"""Whole-acquisition localization and precision benchmark on two EPFL simulations.

See analysis/epfl-precision.md for the CPU environment and a reproducible command.
"""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
from importlib.metadata import version
import io
import json
import math
from pathlib import Path
import platform
import re
import shutil
from urllib.request import urlopen
from zipfile import ZipFile

import numpy as np
from scipy import ndimage, optimize, spatial, special


ROOT = Path(__file__).resolve().parents[1]
ROOT_URL = "https://bigwww.epfl.ch/srm"
DATASETS = {
    "sparse": ("Bundled_Tubes_Long_Sequence", 12000, {
        "sequence.zip": "0c70f0ab2cc73cec7d5ddee80009cd48b2a2e2cec4a223120bc3f0323dd37270",
        "fluorophores.zip": "337a136ee5434495b8a0dc70f12c33023efa01a0948ce89e42430c3ea1274ac9",
        "psf.zip": "cba2c63ba2d6ca18e64e78cca479d8e7a81bcb6cfadcaa0ea6b5837be8c996e6",
    }),
    "dense": ("Bundled_Tubes_High_Density", 361, {
        "sequence.zip": "d6e303f685a9171d1db858fd3674b09db62d162bca842746e2a0c722d187a719",
        "fluorophores.zip": "d57316fb663c543c4e9b9baf1152a82912e634ce5949cb59a9f2cfb0db6ecad3",
        "psf.zip": "60c3a97a6d797d4ccd84842ec3262c3a7c15a28532e3d6f1052ad3377b00bad5",
    }),
}
PIXEL_NM, BASELINE_ADU = 100.0, 100.0
READ_NOISE_E, EXCESS_NOISE = 5.0, 1.41
BOX = 7
MAX_ITERATIONS = 100
EPS_PX = 0.001
MATCH_NM = 150.0
SWEEP_NM = (50.0, 100.0, MATCH_NM, 250.0)
GAUSSIAN95 = -2 * math.log(0.05)


def verified_sources(folder, slug, hashes):
    meta = {}
    for name, expected in hashes.items():
        url = f"{ROOT_URL}/Data/challenge-2D-simulation/{slug}/Data/{name}"
        path = folder / name
        if not path.exists():
            folder.mkdir(parents=True, exist_ok=True)
            partial = path.with_suffix(".part")
            with urlopen(url, timeout=90) as response, partial.open("wb") as out:
                shutil.copyfileobj(response, out)
            if hashlib.sha256(partial.read_bytes()).hexdigest() != expected:
                raise ValueError(f"Incomplete or changed download: {url}")
            partial.replace(path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Data integrity check failed for {path}")
        meta[name] = {"url": url, "sha256": expected, "bytes": path.stat().st_size}
    return meta


def frame_index(archive, pattern):
    listing = {}
    for name in archive.namelist():
        match = re.fullmatch(pattern, name)
        if match:
            frame = int(match[1])
            if frame in listing:
                raise ValueError(f"Duplicate frame {frame} in {archive.filename}")
            listing[frame] = name
    return listing


def identify(image):
    """Positive bandpass maxima. All thresholds are fixed across density conditions."""
    band = ndimage.gaussian_filter(image, 1.0) - ndimage.gaussian_filter(image, 3.0)
    median = np.median(band)
    noise = 1.4826 * np.median(np.abs(band - median))
    peak = (band == ndimage.maximum_filter(band, 3)) & (band > median + 6 * max(noise, 1.0))
    half = BOX // 2
    peak[:half, :] = peak[-half:, :] = False
    peak[:, :half] = peak[:, -half:] = False
    return np.argwhere(peak)


def camera_covariance(theta):
    """Propagate the fitted Gaussian's varying pixel noise through its information matrix.

    The Gaussian is integrated over unit pixel areas. Fit parameters include
    both coordinates, total signal, background and a shared spot width.
    Return the XY marginal covariance, including nuisance-parameter uncertainty.
    """
    grid = np.arange(BOX)[None, :]
    y, x, n, bg, sigma = [theta[:, k, None] for k in range(5)]

    def pixel_integrals(center):
        low, high = (grid - 0.5 - center) / sigma, (grid + 0.5 - center) / sigma
        flo, fhi = np.exp(-low ** 2 / 2) / math.sqrt(2 * math.pi), np.exp(-high ** 2 / 2) / math.sqrt(2 * math.pi)
        return special.ndtr(high) - special.ndtr(low), (flo - fhi) / sigma, (low * flo - high * fhi) / sigma

    px, dx, dsx = pixel_integrals(x)
    py, dy, dsy = pixel_integrals(y)
    outer = lambda a, b: a[:, :, None] * b[:, None, :]
    psf = outer(py, px)
    mu = n[:, :, None] * psf + bg[:, :, None]
    jacobian = np.stack((n[:, :, None] * outer(py, dx), n[:, :, None] * outer(dy, px),
                         psf, np.ones_like(psf), n[:, :, None] *
                         (outer(dsy, px) + outer(py, dsx))), axis=-1)
    variance = EXCESS_NOISE ** 2 * mu + READ_NOISE_E ** 2
    weighted = (jacobian / np.sqrt(variance)[..., None]).reshape(len(theta), BOX ** 2, 5)
    info = weighted.transpose(0, 2, 1) @ weighted
    cov = np.linalg.inv(info)[:, :2, :2] * PIXEL_NM ** 2
    return cov


def one_to_one(detected, truth, radius):
    """First maximize pair count within radius, then minimize summed XY distance."""
    if not len(detected) or not len(truth):
        return np.empty(0, dtype=int), np.empty(0, dtype=int), 0
    distance = spatial.distance.cdist(detected, truth)
    possible = distance <= radius
    # One lost match costs more than any full solution's in-radius distance.
    penalty = (min(distance.shape) + 1) * radius
    real = np.where(possible, distance, 2 * penalty)
    costs = np.column_stack((real, np.full((len(detected), len(detected)), penalty)))
    rows, cols = optimize.linear_sum_assignment(costs)
    keep = cols < len(truth)
    return rows[keep], cols[keep], int(np.sum(possible.sum(axis=1) > 1))


def read_inputs(folder, nframes):
    import tifffile

    with ZipFile(folder / "sequence.zip") as images, ZipFile(folder / "fluorophores.zip") as coordinates:
        frames = frame_index(images, r"sequence/(\d+)\.tif")
        truth_frames = frame_index(coordinates, r"fluorophores/frames/(\d+)\.csv")
        if set(frames) != set(range(1, nframes + 1)) or frames.keys() != truth_frames.keys():
            raise ValueError("TIFF and coordinate frames must form the same complete 1-based interval.")
        truth, pixels, ids, centers = [], [], [], []
        for frame in range(1, nframes + 1):
            image = tifffile.imread(io.BytesIO(images.read(frames[frame]))).astype(np.float32)
            if image.shape != (64, 64):
                raise ValueError("Unexpected image size")
            with coordinates.open(truth_frames[frame]) as fh, io.TextIOWrapper(fh, encoding="ascii") as txt:
                coords = []
                for row in csv.DictReader(txt):
                    if int(row["frame"]) != frame:
                        raise ValueError("Truth frame differs from its member name")
                    coords.append((float(row["xnano"]), float(row["ynano"])))
                truth.append(np.array(coords, dtype=float).reshape(-1, 2))
            for cy, cx in identify(image):
                h = BOX // 2
                roi = image[cy-h:cy+h+1, cx-h:cx+h+1]
                pixels.append(np.maximum(roi - BASELINE_ADU, 0))
                ids.append(frame)
                centers.append((cx, cy))
    return truth, np.asarray(ids), np.asarray(centers), np.asarray(pixels)


def fitted_detections(roi, ids, centers, threads):
    from picasso.gaussmle import gaussmle

    if not len(roi):
        return np.empty(0, int), np.empty((0, 2)), np.empty((0, 2, 2)), np.empty((0, 2)), 0
    def fit(batch):
        return gaussmle(batch, EPS_PX, MAX_ITERATIONS, method="sigma")

    # Warm up the Numba implementation once, then process independent spot blocks.
    fit(roi[:1])
    batches = [roi[start:start+512] for start in range(0, len(roi), 512)]
    with ThreadPoolExecutor(max_workers=threads) as worker:
        results = list(worker.map(fit, batches))
    theta, crlb, _, iterations = (np.concatenate(v) for v in zip(*results))
    theta = theta.astype(float)
    keep = np.isfinite(theta).all(axis=1) & np.isfinite(crlb[:, :2]).all(axis=1)
    keep &= (crlb[:, :2] > 0).all(axis=1) & (theta[:, 2] > 0) & (theta[:, 3] >= 0)
    keep &= (theta[:, 4] >= 0.6) & (theta[:, 4] <= 3.0)
    keep &= (np.abs(theta[:, :2] - BOX // 2) <= 2.0).all(axis=1)
    valid = theta[keep]
    # Picasso theta[0] is the row (y), theta[1] the column (x).
    # The EPFL truth is measured from the FOV edge, 0.5 pixel before the first center.
    xy = (centers[keep] + valid[:, [1, 0]] - BOX // 2 + 0.5) * PIXEL_NM
    blocks = [camera_covariance(valid[j:j+1024]) for j in range(0, len(valid), 1024)]
    cov = np.concatenate(blocks) if blocks else np.empty((0, 2, 2))
    poisson_axis_var = crlb[keep][:, [1, 0]] * PIXEL_NM ** 2
    return ids[keep], xy, cov, poisson_axis_var, int(np.sum(iterations[keep] == MAX_ITERATIONS))


def divided(num, den):
    return float(num / den) if den else None


def accuracy(tp, ntruth, nfound):
    return {"matches": int(tp), "false_positives": int(nfound-tp), "false_negatives": int(ntruth-tp),
            "recall": divided(tp, ntruth), "fdr": divided(nfound-tp, nfound)}


def residual_stats(resid, cov):
    if not len(resid):
        return {"n": 0, "median_radial_error_nm": None, "rms_radial_error_nm": None,
                "coverage95": None, "fraction_radial_error_le20": None}
    err2 = np.einsum("ni,ni->n", resid, resid)
    r = np.sqrt(err2)
    mahalanobis = np.einsum("ni,nij,nj->n", resid, np.linalg.inv(cov), resid)
    sigma = np.sqrt(np.diagonal(cov, axis1=1, axis2=2))
    return {
        "n": len(r), "median_radial_error_nm": float(np.median(r)),
        "rms_radial_error_nm": float(np.sqrt(np.mean(err2))),
        "rmse_xy_nm": np.sqrt(np.mean(resid ** 2, axis=0)).tolist(),
        "mean_signed_error_xy_nm": np.mean(resid, axis=0).tolist(),
        "median_reported_axis_sigma_nm": np.median(sigma, axis=0).tolist(),
        "reported_radial_rms_nm": float(np.sqrt(np.mean(np.trace(cov, axis1=1, axis2=2)))),
        "coverage_1sigma_axis": np.mean(np.abs(resid) <= sigma, axis=0).tolist(),
        "coverage95": float(np.mean(mahalanobis <= GAUSSIAN95)),
        "fraction_radial_error_le20": float(np.mean(r <= 20)),
    }


def evaluate(truth, frames, detected, covariance, poisson_var):
    starts = np.searchsorted(frames, np.arange(1, len(truth) + 2))
    ntruth, nfound = sum(map(len, truth)), len(detected)
    sweep = {radius: 0 for radius in SWEEP_NM}
    all_det, all_true, all_crowded, ambiguous, chance_pairs = [], [], [], 0, 0
    crowded_truth = 0
    for f, coords in enumerate(truth):
        a, b = starts[f:f+2]
        xy = detected[a:b]
        if len(coords) > 1:
            nearest = spatial.cKDTree(coords).query(coords, k=2)[0][:, 1]
            crowded = nearest <= 250.0
        else:
            crowded = np.zeros(len(coords), dtype=bool)
        crowded_truth += int(crowded.sum())
        for radius in SWEEP_NM:
            i, j, amb = one_to_one(xy, coords, radius)
            sweep[radius] += len(i)
            if radius == MATCH_NM:
                all_det.extend((a + i).tolist())
                all_true.extend(coords[j].tolist())
                all_crowded.extend(crowded[j].tolist())
                ambiguous += amb
        # Same static structure, mismatched times: a negative control for
        # coincidental spatial proximity at high density. This is not an FDR correction.
        wrong_frame = truth[(f + len(truth)//2) % len(truth)]
        chance_pairs += len(one_to_one(xy, wrong_frame, MATCH_NM)[0])
    matched = np.asarray(all_det, dtype=int)
    residual = detected[matched] - np.asarray(all_true).reshape(-1, 2)
    mcov = covariance[matched]
    mstats = residual_stats(residual, mcov)
    sigma_max = np.sqrt(np.diagonal(covariance, axis1=1, axis2=2)).max(axis=1)
    selected = sigma_max <= 20
    passes = selected[matched]
    inside = np.linalg.norm(residual[passes], axis=1) <= 20
    poisson_cov = np.zeros_like(mcov)
    poisson_cov[:, 0, 0], poisson_cov[:, 1, 1] = poisson_var[matched].T
    crowded = np.asarray(all_crowded, dtype=bool)
    crowding = {}
    for name, mask, total in (("neighbor_within_250nm", crowded, crowded_truth),
                              ("no_neighbor_within_250nm", ~crowded, ntruth-crowded_truth)):
        crowding[name] = {"truth_count": total, "matches": int(mask.sum()),
                          "recall": divided(mask.sum(), total),
                          **residual_stats(residual[mask], mcov[mask])}
    threshold_bins = []
    low = 0
    for hi in (5, 10, 20, 40, np.inf):
        mask = (sigma_max[matched] > low) & (sigma_max[matched] <= hi)
        threshold_bins.append({"lower_exclusive_nm": low, "upper_inclusive_nm": None if np.isinf(hi) else hi,
                               "matched": int(mask.sum()), **residual_stats(residual[mask], mcov[mask])})
        low = hi
    return {
        "truth_count": ntruth, "detection_count": nfound, **accuracy(len(matched), ntruth, nfound),
        **mstats, "ambiguous_detections": ambiguous,
        "precision20": {"passing_detections": int(selected.sum()), "passing_matched": int(passes.sum()),
                        "passing_unmatched": int(selected.sum()-passes.sum()),
                        "matched_actual_radial_le20": int(inside.sum()),
                        "fraction_actual_radial_le20_given_matched_pass": divided(inside.sum(), passes.sum()),
                        "fraction_all_passes_matched_and_radial_le20": divided(inside.sum(), selected.sum())},
        "uncertainty_bins": threshold_bins, "ground_truth_crowding": crowding,
        "uncorrected_picasso_poisson": residual_stats(residual, poisson_cov),
        "radius_sensitivity": [{"match_radius_nm": radius, **accuracy(tp, ntruth, nfound)}
                               for radius, tp in sweep.items()],
        "wrong_frame_control": {"shift_frames": len(truth)//2, **accuracy(chance_pairs, ntruth, nfound)},
    }


def check_psf(path):
    import tifffile
    with ZipFile(path) as z, tifffile.TiffFile(io.BytesIO(z.read("psf/psf-10.0nm.tif"))) as t:
        profile = t.pages[0].asarray()[103].astype(float)
        profile /= profile.max()
        j = np.flatnonzero(profile >= 0.5)
        a, b = j[0], j[-1]
        left = a-1 + (0.5 - profile[a-1])/(profile[a] - profile[a-1])
        right = b + (profile[b] - 0.5)/(profile[b] - profile[b+1])
        return {"archive_member": "psf/psf-10.0nm.tif", "pages": len(t.pages),
                "focal_slice_index": 0, "profile_spacing_nm_from_member_name": 10,
                "measured_x_fwhm_nm": float(10 * (right-left)), "nominal_fwhm_nm": 258.21}


def run(cache, workers):
    outputs = []
    for density, (slug, nframes, hashes) in DATASETS.items():
        folder = cache / slug
        sources = verified_sources(folder, slug, hashes)
        truth, ids, centers, roi = read_inputs(folder, nframes)
        print(f"{density}: {len(ids):,} candidates from {nframes:,} TIFF frames", flush=True)
        frames, xy, cov, raw, hit_cap = fitted_detections(roi, ids, centers, workers)
        stats = evaluate(truth, frames, xy, cov, raw)
        outputs.append({"id": density, "dataset": slug, "frames": nframes, "archive_frame_ids": [1, nframes],
                        "source_url": f"{ROOT_URL}/dataset/challenge-2D-simulation/{slug}/index.html",
                        "sources": sources, "psf": check_psf(folder/"psf.zip"), "candidates": len(ids),
                        "rejected_fits": int(len(ids)-len(frames)), "accepted_at_iteration_cap": hit_cap,
                        **stats})
        print(f"{density}: recall={stats['recall']:.4f}, FDR={stats['fdr']:.4f}, "
              f"RMS XY distance={stats['rms_radial_error_nm']:.2f} nm", flush=True)
    return {
        "schema": "epfl-image-localization/1", "scope": "Full simulated 2D sensor sequences",
        "dataset_doi": "10.1038/nmeth.3442", "software_doi": "10.1038/nprot.2017.024",
        "gaussmle_source": "https://github.com/jungmannlab/picasso/blob/v0.7.5/picasso/gaussmle.py",
        "versions": {"python": platform.python_version(), **{pkg: version(pkg) for pkg in
                     ("numpy", "scipy", "tifffile", "numba", "picassosr")}},
        "pixel_nm": PIXEL_NM, "camera": {"baseline_adu": BASELINE_ADU, "electrons_per_adu": 1,
            "quantum_efficiency": 1, "read_noise_e_rms": READ_NOISE_E, "excess_noise_factor": EXCESS_NOISE},
        "settings": {"finder": "SciPy difference of Gaussians (1 and 3 px), 3x3 local maxima",
            "threshold": "Frame median plus 6 times max(1 ADU, 1.4826 times frame MAD) on the bandpass",
            "fit": "Picasso 0.7.5 gaussmle sigma, 7x7 pixels, single circular pixel-integrated Gaussian plus constant background",
            "fit_epsilon_px": EPS_PX, "maximum_iterations": MAX_ITERATIONS,
            "filter": "Finite fit and positive raw positional CRLB; positive photons, nonnegative background, 0.6 to 3 px width, fit within 2 px of seed on each axis",
            "coordinates": "1-based TIFF/GT filename IDs; array indexing [y,x], EPFL xy_nm = (zero-based fitted center index + 0.5) * 100",
            "uncertainty": "5-parameter marginal inverse information from the fitted mean and noise variance 1.41^2*mean+5^2 photoelectron-equivalent units; raw Picasso Poisson CRLB is a separate comparison",
            "primary_match_radius_nm": MATCH_NM, "match": "Same frame, maximum cardinality then minimum total radial distance; inclusive radius",
            "scope_20nm": "Both reported axis standard deviations at most 20 nm; ground-truth test separately uses Euclidean radial error at most 20 nm",
            "coverage95": "Fraction of matched residuals inside the reported joint 95% Gaussian covariance ellipse, squared Mahalanobis radius <=5.99146",
            "source_note": "Dense landing-page summary reports 168 frames; the numbered sensor and truth archives both contain frames 1 through 361."},
        "limitations": ["These synthetic tubulin images test a single-emitter estimator under the archived camera conditions.",
            "Calibration and error distributions condition on matched detections. The matching radius truncates the error distribution and dense same-frame matches can occur by chance.",
            "Neighboring sources, structured fluorescence and PSF mismatch can bias fits with small predicted uncertainty. The synthetic cases contain no JF dye comparison."],
        "datasets": outputs,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path.home()/".cache/blink-atlas/epfl-precision")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--out", type=Path, default=ROOT/"docs/data/epfl-precision.json")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    data = run(args.cache, args.workers)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")

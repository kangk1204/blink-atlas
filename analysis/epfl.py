"""Verify the same-frame collision geometry against EPFL Tubulins_II ground truth.

Run `python -m analysis.epfl` from the repository root. The small official
archive can also be supplied with --archive for offline reproduction.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import zipfile
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen


DATASET_URL = "https://bigwww.epfl.ch/srm/dataset/challenge-2D-simulation/Tubulins_II/index.html"
ARCHIVE_URL = "https://bigwww.epfl.ch/srm/Data/challenge-2D-simulation/Tubulins_II/Data/fluorophores.zip"
ARCHIVE_SHA256 = "076f4503f67be9808bf81fb47e5e3a5822fe560d0d9a7a4db9e6e9e633e2c205"
FIELD_WIDTH_NM = 38_400.0
BENCHMARK_FWHM_NM = 258.21
DEMO_EXCLUSION_NM = 250.0
DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "docs/data/epfl.json"


def all_pairs(points, radius_nm):
    """Exact reference with one comparison per unordered pair."""
    pairs = set()
    for i, (xi, yi) in enumerate(points):
        for j in range(i):
            xj, yj = points[j]
            if math.hypot(xi - xj, yi - yj) <= radius_nm:
                pairs.add((j, i))
    return pairs


def grid_pairs(points, radius_nm):
    """Enumerate colliding centers via the nine relevant spatial bins."""
    if not math.isfinite(radius_nm) or radius_nm <= 0:
        raise ValueError("The collision radius must be finite and positive.")
    grid = defaultdict(list)
    pairs = set()
    for i, (x, y) in enumerate(points):
        cx, cy = math.floor(x / radius_nm), math.floor(y / radius_nm)
        for gx in range(cx - 1, cx + 2):
            for gy in range(cy - 1, cy + 2):
                for j in grid.get((gx, gy), ()):
                    xj, yj = points[j]
                    if (x - xj) ** 2 + (y - yj) ** 2 <= radius_nm ** 2:
                        pairs.add((j, i))
        grid[cx, cy].append(i)
    return pairs


def ground_truth_frames(archive: zipfile.ZipFile):
    """Read the author's explicit frame-specific list without inferring lifetimes."""
    filenames = sorted(name for name in archive.namelist() if re.fullmatch(r"frames/\d+\.csv", name))
    if not filenames:
        raise ValueError("No frame-wise coordinates found.")
    for expected_frame, name in enumerate(filenames, start=1):
        if int(Path(name).stem) != expected_frame:
            raise ValueError(f"Missing or out-of-order frame before {name}")
        points = []
        with archive.open(name) as data, io.TextIOWrapper(data, encoding="ascii") as handle:
            for row in csv.DictReader(handle, skipinitialspace=True):
                x, y = float(row["xnano"]), float(row["ynano"])
                if int(row["frame"]) != expected_frame:
                    raise ValueError(f"Frame label disagrees with {name}")
                if not (0 <= x < FIELD_WIDTH_NM and 0 <= y < FIELD_WIDTH_NM):
                    raise ValueError(f"Non-finite or out-of-field coordinate in {name}")
                points.append((x, y))
        yield points


def benchmark(payload: bytes):
    digest = hashlib.sha256(payload).hexdigest()
    if digest != ARCHIVE_SHA256:
        raise ValueError(f"Unexpected data digest {digest}. Check the EPFL archive version.")
    totals = {
        r: dict(overlap_pairs=0, overlapping_emitters=0, pair_mismatches=0,
                mask_mismatches=0, uniform_expected_overlaps=0.0)
        for r in (BENCHMARK_FWHM_NM, DEMO_EXCLUSION_NM)
    }
    frame_count = emitter_count = 0
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        for points in ground_truth_frames(archive):
            frame_count += 1
            n = len(points)
            emitter_count += n
            for radius, metrics in totals.items():
                found = grid_pairs(points, radius)
                truth = all_pairs(points, radius)
                found_indices = {i for pair in found for i in pair}
                truth_indices = {i for pair in truth for i in pair}
                metrics["overlap_pairs"] += len(truth)
                metrics["overlapping_emitters"] += len(truth_indices)
                metrics["pair_mismatches"] += len(found ^ truth)
                metrics["mask_mismatches"] += len(found_indices ^ truth_indices)
                # Condition on frame count. Uniform Poisson uses the observed
                # density of the other emitters, averaged over emitting centers.
                occupancy = max(n - 1, 0) * math.pi * radius ** 2 / FIELD_WIDTH_NM ** 2
                metrics["uniform_expected_overlaps"] += n * -math.expm1(-occupancy)

    summaries = {}
    for radius, metrics in totals.items():
        overlaps = metrics["overlapping_emitters"]
        uniform = metrics.pop("uniform_expected_overlaps") / emitter_count
        summaries[radius] = {
            "radius_nm": radius, **metrics, "isolated_emitters": emitter_count - overlaps,
            "overlap_fraction": overlaps / emitter_count,
            "uniform_poisson_overlap_fraction": uniform,
            "agreement": 1 - metrics["mask_mismatches"] / emitter_count,
        }
    return {
        "dataset": "EPFL 2013 SMLM Challenge / Tubulins_II",
        "source_url": DATASET_URL, "archive_url": ARCHIVE_URL, "archive_sha256": digest,
        "archive_bytes": len(payload), "source_member": "frames/00001.csv through frames/02401.csv",
        "coordinates": "xnano and ynano; XY projection, nanometres, inclusive radius, same published frame",
        "field_width_nm": FIELD_WIDTH_NM, "frames": frame_count, "emitters": emitter_count,
        **summaries[BENCHMARK_FWHM_NM], "at_demo_radius": summaries[DEMO_EXCLUSION_NM],
        "method": (
            "Each published frame-wise emitter coordinate is counted once. Two centers overlap "
            "when their XY separation is at most the stated radius. An emitter is isolated if "
            "it has zero such neighbors. Spatial-bin pairs and masks are compared to exhaustive "
            "pairwise distances in every frame. The first radius is the challenge in-focus PSF FWHM."
        ),
        "uniform_comparison": (
            "Per frame with n emitters, the homogeneous Poisson proxy gives a tagged emitter "
            "overlap probability of 1 - exp[-(n-1) pi r^2 / field_area]. The average weights each "
            "frame by n. Tubulin structure concentrates the emitters in a small part of the field."
        ),
        "interpretation": (
            "This verifies projected same-frame geometry on a simulated tubulin specimen. The "
            "photons, sensor noise and dye identity are not fitted. The difference from a "
            "homogeneous field quantifies the sensitivity to spatial clustering. Fixed-cell "
            "HaloTag dye rankings require their own matched measurements."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, help="Official fluorophores.zip for offline use")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.archive:
        payload = args.archive.read_bytes()
    else:
        with urlopen(ARCHIVE_URL, timeout=60) as archive:
            payload = archive.read()
    result = benchmark(payload)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for summary in (result, result["at_demo_radius"]):
        print(f'{summary["radius_nm"]:g} nm: {summary["overlap_pairs"]:,} pairs, '
              f'{summary["overlapping_emitters"]:,} overlapping / {result["emitters"]:,} emitter-frames, '
              f'{summary["pair_mismatches"]} pair mismatches')
    if any(s["pair_mismatches"] or s["mask_mismatches"] for s in (result, result["at_demo_radius"])):
        raise SystemExit("Geometry check failed")


if __name__ == "__main__":
    main()

"""Recompute Steen et al. DNA-PAINT metrics from deposited Picasso localization rows.

Only the selected localization columns are loaded, one file at a time. The
archive manifest and per-case checkpoint driver are in palette_raw_run.py.
"""

import math
from pathlib import Path
import warnings

import numpy as np
from numba import njit
from scipy.optimize import curve_fit, leastsq
from scipy.special import erf, iv, i0e


PIXEL_NM = 130.0
EXPOSURE_S = 0.1


def fitted_tau(times):
    """Run the author's one-parameter LM fit; retain successful positive widths."""
    if times.size < 2:
        return np.nan
    ts = np.sort(times.astype(np.float64))
    rank = np.arange(len(ts), dtype=float) / (len(ts) - 1)

    def residual(theta):
        value = rank - (1 - np.exp(-(ts / theta[0])))
        if not np.isfinite(value).all():
            raise ValueError("Non-finite CDF evaluation")
        return value

    # lmfit.Model.leastsq delegates to SciPy MINPACK with these settings.
    # The released initialization at 200 frames matters for sparse sites.
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        try:
            fit = leastsq(
                residual,
                [200.0],
                ftol=1.5e-8,
                xtol=1.5e-8,
                gtol=0.0,
                epsfcn=1e-10,
                maxfev=8000,
                full_output=True,
            )
        except (ValueError, RuntimeError, FloatingPointError):
            return np.nan
    width = float(fit[0][0])
    return width if np.isfinite(width) and width > 0 else np.nan


def destroyed_fraction(ratios, all_groups):
    """Steen residual-dark histogram fit and chance-tail correction (in percent)."""
    r = ratios[np.isfinite(ratios) & (ratios > 0)]
    if len(r) < 3 or r.max() - r.min() < 1:
        return None, {
            "eligibleSites": int(len(r)),
            "reason": "Insufficient residual-time histogram support.",
        }
    width = 0.5
    bins = np.arange(r.min(), r.max() + width, width)
    n = np.histogram(r, bins)[0]
    x = bins[:-1] + width / 2

    def exponential(a, amplitude, decay):
        return amplitude * np.exp(-decay * a)

    try:
        with np.errstate(over="ignore", invalid="ignore"):
            result, _ = curve_fit(exponential, x, n, p0=(500, 1), maxfev=20000)
        if not np.isfinite(result[1]) or result[1] <= 0:
            raise ValueError("Nonpositive residual-time decay")
    except (ValueError, RuntimeError, FloatingPointError):
        return None, {
            "eligibleSites": int(len(r)),
            "reason": "Exponential residual-time fit failed.",
        }
    mu = 1.0 / result[1]
    threshold = max(4.0, 4.0 * mu)
    number = int(np.sum(r >= threshold))
    # The released code assigns negative r to unfit sites and retains them
    # in the denominator; their contribution is zero in the upper tail.
    value = 100 * (number / all_groups - math.exp(-4))
    return float(value), {
        "eligibleSites": len(r),
        "allObservedSites": all_groups,
        "tailSites": number,
        "fitMu": float(mu),
        "ratioThreshold": float(threshold),
        "eligibleDenominatorPct": float(100 * (number / len(r) - math.exp(-4))),
    }


def mean_or_none(a):
    good = a[np.isfinite(a) & (a > 0)]
    return float(np.mean(good, dtype=np.float64)) if len(good) else None


def picked_metrics(locs, frames, origami):
    if not len(locs) or frames <= 0:
        raise ValueError(
            "The picked acquisition must have rows and a positive frame count."
        )
    order = np.lexsort((locs["frame"], locs["group"]))
    records = locs[order]
    frame = records["frame"].astype(np.int64)
    if frame.min() < 0 or frame.max() >= frames:
        raise ValueError("Frame outside acquisition length")
    same_group = records["group"][1:] == records["group"][:-1]
    frame_difference = np.diff(frame)
    connected = same_group & (frame_difference == 1)
    starts = np.r_[0, np.flatnonzero(~connected) + 1]
    ends = np.r_[starts[1:], len(records)]
    durations = ends - starts
    event_groups = records["group"][starts]
    boundaries = np.r_[
        0, np.flatnonzero(event_groups[1:] != event_groups[:-1]) + 1, len(starts)
    ]
    groups = list(zip(boundaries[:-1], boundaries[1:]))
    bright = np.array([fitted_tau(durations[a:b]) for a, b in groups])
    dark = np.array([fitted_tau(np.diff(frame[starts[a:b]])) for a, b in groups])
    group_last_row = ends[boundaries[1:] - 1] - 1
    residual = frames - frame[group_last_row]
    center = np.zeros(len(frame), bool)
    center[1:-1] = connected[:-1] & connected[1:]
    duplicates = int(np.count_nonzero(same_group & (frame_difference == 0)))
    sites = len(bright)
    p = records["photons"][center].astype(np.float64)
    sx, sy, bg = (records[k][center].astype(np.float64) for k in ("sx", "sy", "bg"))
    valid = (
        np.isfinite(p)
        & (p > 0)
        & np.isfinite(sx)
        & (sx > 0)
        & np.isfinite(sy)
        & (sy > 0)
        & np.isfinite(bg)
        & (bg > 0)
    )
    sbr = (
        p[valid]
        * erf(0.5 / (math.sqrt(2) * sx[valid]))
        * erf(0.5 / (math.sqrt(2) * sy[valid]))
        / bg[valid]
    )
    first = int(np.sum(frame <= 0.2 * frames))
    last = int(np.sum(frame >= 0.8 * frames))
    event_totals = np.add.reduceat(records["photons"].astype(np.float64), starts)
    finite_event_totals = event_totals[np.isfinite(event_totals)]
    metrics = {
        "centerPhotons": mean_or_none(p),
        "centerSbr": mean_or_none(sbr),
        "eventPhotons": float(np.mean(finite_event_totals))
        if len(finite_event_totals)
        else None,
        "brightTimeS": mean_or_none(bright * EXPOSURE_S),
        "darkTimeS": mean_or_none(dark * EXPOSURE_S),
        "dropPct": 100.0 * (1 - last / first) if first else None,
        "destroyedPct": None,
        "specificity": None,
        "nenaNm": None,
    }
    detail = {
        "groups": int(sites),
        "events": int(len(starts)),
        "centerFrames": int(center.sum()),
        "centerPhotonRows": int(np.count_nonzero(np.isfinite(p) & (p > 0))),
        "sbrRows": int(valid.sum()),
        "duplicateFrameWithinGroup": int(duplicates),
        "brightFits": int(np.count_nonzero(np.isfinite(bright))),
        "darkFits": int(np.count_nonzero(np.isfinite(dark))),
        "totalFrames": frames,
        "firstFifthRows": first,
        "lastFifthRows": last,
        "localizations": len(records),
        "eventPhotonTotals": {
            "finiteEvents": int(len(finite_event_totals)),
            "median": float(np.median(finite_event_totals))
            if len(finite_event_totals)
            else None,
            "pooledEdgeFrameMean": mean_or_none(records["photons"][~center]),
        },
    }
    if origami:
        ratios = np.divide(
            residual, dark, out=np.full(len(dark), np.nan), where=dark > 0
        )
        metrics["destroyedPct"], detail["destruction"] = destroyed_fraction(
            ratios, sites
        )
    return metrics, detail


@njit(cache=True)
def _adjacent_distance_hist(frame, x, y, group, enforce_group, source_rows):
    """All same-group pairs exactly one frame apart, 0.001 px bins below 1 pixel."""
    n = len(x)
    hist = np.zeros(1000, dtype=np.int64)
    if n < 2:
        return hist
    xmin, ymin = math.floor(x.min()), math.floor(y.min())
    nx, ny = math.floor(x.max()) - xmin + 1, math.floor(y.max()) - ymin + 1
    if nx * ny > 4_000_000:
        raise ValueError("Unexpected coordinate extent; verify detector pixel units")
    heads = np.full((ny, nx), -1, np.int64)
    nex = np.empty(n, np.int64)
    starts = np.empty(n + 1, np.int64)
    count = 0
    for i in range(n):
        if i == 0 or frame[i] != frame[i - 1]:
            starts[count] = i
            count += 1
    starts[count] = n
    for k in range(count - 1):
        a, b, c = starts[k], starts[k + 1], starts[k + 2]
        if frame[b] != frame[a] + 1:
            continue
        for j in range(b, c):
            ix, iy = math.floor(x[j]) - xmin, math.floor(y[j]) - ymin
            nex[j] = heads[iy, ix]
            heads[iy, ix] = j
        for i in range(a, min(b, source_rows)):
            ix, iy = math.floor(x[i]) - xmin, math.floor(y[i]) - ymin
            for iy2 in range(max(0, iy - 1), min(ny, iy + 2)):
                for ix2 in range(max(0, ix - 1), min(nx, ix + 2)):
                    j = heads[iy2, ix2]
                    while j != -1:
                        if not enforce_group or group[i] == group[j]:
                            dx, dy = x[i] - x[j], y[i] - y[j]
                            dist2 = dx**2 + dy**2
                            if dist2 < 1.0:
                                hist[min(999, int(math.sqrt(dist2) * 1000))] += 1
                        j = nex[j]
        for j in range(b, c):
            heads[math.floor(y[j]) - ymin, math.floor(x[j]) - xmin] = -1
    return hist


def next_frame_histogram(locs, historical_batches=False):
    frame, x, y = (locs[k] for k in ("frame", "x", "y"))
    valid = np.isfinite(x) & np.isfinite(y)
    if frame.dtype.kind == "f":
        valid &= np.isfinite(frame) & (frame == np.floor(frame))
    indices = np.flatnonzero(valid)
    if np.any(np.diff(frame[indices].astype(np.int64)) < 0):
        indices = indices[np.argsort(frame[indices], kind="stable")]
    has_group = "group" in locs.dtype.names
    group = (
        locs["group"][indices].astype(np.int64) if has_group else np.zeros(1, np.int64)
    )
    # Historical Picasso loops over exactly 100 integer-sized batches of
    # source localizations; any remainder is still eligible as a target.
    limit = (len(indices) // 100) * 100 if historical_batches else len(indices)
    return _adjacent_distance_hist(
        frame[indices].astype(np.int64), x[indices], y[indices], group, has_group, limit
    )


def nena_fit(hist, median_lp, pixel_nm=PIXEL_NM):
    """Prepublication Picasso postprocess.nena model, initialization and extracted s."""
    from lmfit import Model, Parameters

    centers = np.arange(0.0, 1.0, 0.001) + 0.0005
    if hist.sum() < 100 or not math.isfinite(median_lp) or median_lp <= 0:
        raise ValueError("Too few next-frame pairs or missing lpx/lpy initialization")

    def f(d, a, s, ac, dc, sc):
        s = max(s, 1e-9)
        sc = max(sc, 1e-9)
        z = d * dc / sc  # literal argument in Picasso's historical source
        exponent = -0.5 * (d**2 + dc**2) / sc**2
        first = a * (d / s**2) * np.exp(-0.5 * d**2 / s**2)
        second = ac * (d / sc**2) * np.exp(exponent) * iv(0, z)
        invalid = ~np.isfinite(second)
        if np.any(invalid):
            # Stable product for extreme exploratory parameter steps.
            second[invalid] = (
                ac
                * (d[invalid] / sc**2)
                * np.exp(np.minimum(700.0, exponent[invalid] + np.abs(z[invalid])))
                * i0e(z[invalid])
            )
        return first + second

    total = float(np.sum((hist[1:] + hist[:-1]) * np.diff(centers) / 2))
    params = Parameters()
    for name, val in zip(
        ("a", "s", "ac", "dc", "sc"),
        (total / 2, median_lp, total / 2, 2 * median_lp, median_lp),
    ):
        params.add(name, value=val, min=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = Model(f).fit(hist.astype(float), params, d=centers, max_nfev=20000)
    if not result.success or not np.isfinite(result.best_values["s"]):
        raise RuntimeError("Picasso two-component NeNA fit did not converge")
    width = result.best_values["s"]
    return float(width * pixel_nm), {
        "pairCount": int(hist.sum()),
        "binSizePx": 0.001,
        "maxRadiusPx": 1,
        "medianLpInitialPx": median_lp,
        "fit": {k: float(v) for k, v in result.best_values.items()},
        "relativeResidualL2": float(
            np.linalg.norm(result.residual) / np.linalg.norm(hist)
        ),
        "peakDistancePx": float(centers[np.argmax(hist)]),
        "referenceTag": "deab0c6e31c544fc56d1e4d0ccc8e00348dcae92",
    }


def file_fields(path, columns):
    import h5py

    with h5py.File(path) as fh:
        ds = fh["locs"]
        if not isinstance(ds, h5py.Dataset):
            raise ValueError(f"Unexpected Picasso HDF5 layout in {path}")
        missing = set(columns) - set(ds.dtype.names or ())
        if missing:
            raise ValueError(f"Missing columns {sorted(missing)} in {Path(path).name}")
        return ds.fields(columns)[:]


def field_precision(path, pixel_nm=PIXEL_NM):
    import h5py

    with h5py.File(path) as hf:
        ds = hf["locs"]
        columns = ["frame", "x", "y", "lpx", "lpy"]
        if "group" in ds.dtype.names:
            columns.append("group")
        locs = ds.fields(columns)[:]
    init = float(np.mean([np.nanmedian(locs["lpx"]), np.nanmedian(locs["lpy"])]))
    hist = next_frame_histogram(locs, historical_batches=True)
    nm, details = nena_fit(hist, init, pixel_nm=pixel_nm)
    details["pixelSizeNm"] = pixel_nm
    details["sourceBatchRemainder"] = int(
        (np.isfinite(locs["x"]) & np.isfinite(locs["y"])).sum() % 100
    )
    return nm, details

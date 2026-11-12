"""Convert local SMLM localizations to a small, browser-readable pilot summary.

Run `python -m analysis.pilot --help` for inputs. A CSV contains frame, x_nm,
y_nm and photons, with positions referenced to the image's upper-left corner.
Picasso positions are in camera pixels with integer coordinates at pixel centers.
No image, individual localization, filename or local path enters the output.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from pathlib import Path

import numpy as np


def positive(value, name):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a positive finite number") from exc
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be a positive finite number")
    return result


def positive_int(value, name):
    number = positive(value, name)
    if not number.is_integer():
        raise ValueError(f"{name} must be a positive integer")
    return int(number)


def picasso_calibration(path):
    """Read a Picasso sidecar or its embedded HDF5 /metadata JSON list."""
    if path.suffix.lower() in (".h5", ".hdf5"):
        import h5py

        with h5py.File(path) as hf:
            if "metadata" not in hf:
                return {}
            values = json.loads(hf["metadata"][()])
        if not isinstance(values, list):
            raise ValueError("Picasso /metadata must be a JSON document list")
        entries = [d for d in values if isinstance(d, dict)]
    else:
        import yaml

        with path.open(encoding="utf-8-sig") as handle:
            try:
                entries = [d for d in yaml.safe_load_all(handle) if isinstance(d, dict)]
            except yaml.YAMLError as exc:
                raise ValueError("Invalid Picasso YAML metadata") from exc
    meta = {}
    for key, output in [("Frames", "total_frames"), ("Width", "width_px"),
                        ("Height", "height_px"), ("Pixelsize", "pixel_size_nm")]:
        vals = {float(d[key]) for d in entries if d.get(key) is not None and float(d[key]) > 0}
        if len(vals) > 1:
            raise ValueError(f"Conflicting Picasso {key} values in the YAML")
        if vals:
            meta[output] = vals.pop()
    for d in entries:
        camera = d.get("Micro-Manager Metadata") or {}
        if isinstance(camera, dict) and camera.get("Exposure-ms"):
            meta["exposure_ms"] = float(camera["Exposure-ms"])
    rois = [d["ROI"] for d in entries if d.get("ROI")]
    if rois:
        # Picasso stores either ((y0,x0),(y1,x1)) or a list of those rectangles.
        rect = np.asarray(rois[-1], dtype=float).reshape(-1, 2, 2)
        if not np.isfinite(rect).all():
            raise ValueError("Invalid Picasso ROI corners")
        meta["picasso_rois_px"] = rect.tolist()
    return meta


def acquisition(settings, calibration, source_format):
    meta = dict(settings)
    for key, value in calibration.items():
        if key == "picasso_rois_px":
            continue
        if key in meta and not math.isclose(float(meta[key]), float(value), rel_tol=0.002):
            raise ValueError(f"Supplied {key} disagrees with the Picasso YAML")
        meta[key] = value
    if source_format == "picasso":
        meta["pixel_size_nm"] = positive(meta.get("pixel_size_nm"), "pixel_size_nm")
    for axis in ("width", "height"):
        px = meta.get(axis + "_px")
        if px is not None and meta.get("pixel_size_nm"):
            actual = positive_int(px, axis + "_px") * positive(meta["pixel_size_nm"], "pixel_size_nm")
            if axis + "_nm" in meta and not math.isclose(float(meta[axis + "_nm"]), actual, rel_tol=0.001):
                raise ValueError(f"{axis}_nm disagrees with the calibrated pixel dimensions")
            meta[axis + "_nm"] = actual
    cleaned = {k: positive(meta.get(k), k) for k in ("width_nm", "height_nm", "exposure_ms")}
    cleaned["total_frames"] = positive_int(meta.get("total_frames"), "total_frames")
    if cleaned["total_frames"] > 2**40:
        raise ValueError("total_frames exceeds the supported integer range")
    for key in ("frame_interval_ms", "excitation_nm", "pixel_size_nm"):
        cleaned[key] = positive(meta[key], key) if meta.get(key) is not None else None
    if cleaned["frame_interval_ms"] and cleaned["frame_interval_ms"] < 0.998 * cleaned["exposure_ms"]:
        raise ValueError("Frame interval is shorter than the camera exposure")
    for k in ("width_px", "height_px"):
        if k in meta:
            cleaned[k] = positive_int(meta[k], k)
    for k in ("sample_state", "cell_line", "target", "tag", "reference_dye"):
        value = meta.get(k, "unknown")
        if not isinstance(value, str) or not value.strip() or len(value) > 100:
            raise ValueError(f"{k} must be a short, non-empty string")
        cleaned[k] = value.strip()
    width, height = cleaned["width_nm"], cleaned["height_nm"]
    native = []
    for (y0, x0), (y1, x1) in calibration.get("picasso_rois_px", []):
        px = cleaned["pixel_size_nm"]
        native.append([min(x0, x1) * px, min(y0, y1) * px, max(x0, x1) * px, max(y0, y1) * px])
    if "roi_nm" in settings:
        if not isinstance(settings["roi_nm"], (tuple, list)):
            raise ValueError("roi_nm requires four detector coordinates [left,top,right,bottom]")
        roi = list(map(float, settings["roi_nm"]))
    elif len(native) == 1:
        roi = native[0]
    elif len(native) > 1:
        raise ValueError("Select a single Picasso analysis rectangle in roi_nm [left,top,right,bottom] for this pilot")
    else:
        roi = [0, 0, width, height]
    if len(roi) != 4 or not all(math.isfinite(v) for v in roi):
        raise ValueError("roi_nm requires four finite detector coordinates [left,top,right,bottom]")
    left, top, right, bottom = roi
    if not (0 <= left < right <= width and 0 <= top < bottom <= height):
        raise ValueError("The analysis ROI must lie inside the declared camera field")
    if native and not any(x0 <= left < right <= x1 and y0 <= top < bottom <= y1 for x0, y0, x1, y1 in native):
        raise ValueError("The analysis ROI must lie inside a Picasso localization ROI")
    cleaned.update(camera_width_nm=width, camera_height_nm=height, roi_nm=roi,
                   width_nm=right - left, height_nm=bottom - top)
    cleaned.update(localization_source=source_format, coordinate_origin="upper-left outer pixel edge")
    return cleaned


def localization_chunks(path, fields, chunk_size=65536):
    """Yield numeric columns (frame,x,y,photons) in input units, in bounded blocks."""
    if path.suffix.lower() == ".csv":
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or not set(fields).issubset(reader.fieldnames):
                raise ValueError("CSV must include frame,x_nm,y_nm,photons")
            while True:
                rows = list(itertools.islice(reader, chunk_size))
                if not rows:
                    return
                try:
                    yield tuple(np.asarray([r[name] for r in rows], dtype=float) for name in fields)
                except (ValueError, TypeError) as exc:
                    raise ValueError("All CSV localization values must be numeric") from exc
    else:
        import h5py

        with h5py.File(path) as handle:
            if "locs" not in handle:
                raise ValueError("Picasso HDF5 must contain a /locs dataset")
            locs = handle["locs"]
            if not isinstance(locs, h5py.Dataset) or not locs.dtype.names or not set(fields).issubset(locs.dtype.names):
                raise ValueError("Picasso /locs must include frame,x,y,photons")
            for idx in range(0, len(locs), chunk_size):
                data = locs[idx:idx + chunk_size]
                yield tuple(np.asarray(data[k], dtype=float) for k in fields)


def predictions(pilot_counts, n_future):
    """Predeclared damped trend, fit only on the pilot.

    Split the pilot into four contiguous blocks; shrink the least-squares
    slope of log frame rates by 50%, anchor to the pilot mean, and damp the
    future log-rate change exponentially with a 200-frame time constant.
    """
    mean = np.mean(pilot_counts)
    baseline = np.full(n_future, mean)
    if len(pilot_counts) < 4 or mean == 0:
        return baseline.copy(), baseline
    t = np.arange(len(pilot_counts))
    blocks = np.array_split(t, 4)
    centers = np.array([g.mean() for g in blocks])
    log_rates = np.array([math.log((pilot_counts[g].sum() + 0.5) / len(g)) for g in blocks])
    slope = 0.5 * np.polyfit(centers, log_rates, 1)[0]
    last = len(t) - 1
    rate_at_end = mean / np.exp(slope * (t - last)).mean()
    horizons = np.arange(1, n_future + 1)
    future = rate_at_end * np.exp(slope * 200 * -np.expm1(-horizons / 200))
    return future, baseline


def evaluate(counts, n_pilot, block=100):
    n_future = (len(counts) - n_pilot) // block * block
    if not n_future:
        return None
    model, baseline = predictions(counts[:n_pilot], n_future)
    future = counts[n_pilot:n_pilot + n_future]
    observed = future.reshape(-1, block).sum(axis=1)
    aggregate_model = model.reshape(-1, block).sum(axis=1)
    aggregate_base = baseline.reshape(-1, block).sum(axis=1)
    # Report fractional counts at a precision reproducible across numerical libraries.
    model_mae = round(float(np.abs(aggregate_model - observed).mean()), 9)
    baseline_mae = round(float(np.abs(aggregate_base - observed).mean()), 9)
    improved = model_mae < baseline_mae
    return {
        "start_frame": n_pilot, "frames": n_future, "block_frames": block,
        "model": "Half-slope damped log-rate trend", "baseline": "Constant pilot mean",
        "unit": f"localizations per {block} frames",
        "blocks": [{"start_frame": n_pilot + i*block, "frames": block,
                    "observed": int(obs), "prediction": round(float(pred), 9), "baseline": round(float(base), 9)}
                   for i, (obs, pred, base) in enumerate(zip(observed, aggregate_model, aggregate_base))],
        "model_mae": model_mae, "baseline_mae": baseline_mae, "model_improved": improved,
        "selected_forecast": "model" if improved else "baseline",
    }


def summarize(path, settings, *, yaml_path=None, frame_base=0, min_photons=100,
              tile_side_nm=250, pilot_frames=200):
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix not in (".csv", ".hdf5", ".h5"):
        raise ValueError("Supply an SMLM CSV or Picasso HDF5 file")
    if frame_base not in (0, 1):
        raise ValueError("frame_base must be 0 or 1")
    tile_side_nm = positive(tile_side_nm, "tile_side_nm")
    min_photons = float(min_photons)
    if not math.isfinite(min_photons) or min_photons < 0:
        raise ValueError("min_photons must be finite and nonnegative")
    calibration = {}
    if suffix == ".csv":
        kind = "csv"
        fields = ("frame", "x_nm", "y_nm", "photons")
    else:
        kind = "picasso"
        fields = ("frame", "x", "y", "photons")
        if frame_base:
            raise ValueError("Picasso frame indices are zero-based")
        companion = Path(yaml_path) if yaml_path else path.with_suffix(".yaml")
        if companion.is_file():
            calibration = picasso_calibration(companion)
        elif yaml_path:
            raise ValueError("Could not find the specified Picasso YAML")
        else:
            calibration = picasso_calibration(path)
    meta = acquisition(settings, calibration, kind)
    width, height = meta["width_nm"], meta["height_nm"]
    n_pilot = min(positive_int(pilot_frames, "pilot_frames"), meta["total_frames"])
    length = min(n_pilot + 800, meta["total_frames"])
    evaluated_until = n_pilot + (length - n_pilot) // 100 * 100
    counts = np.zeros(length, dtype=np.int64)
    nx = int(width // tile_side_nm)
    ny = int(height // tile_side_nm)
    n_tiles = nx * ny
    if n_tiles <= 0 or n_tiles * n_pilot > 2**53 - 1:
        raise ValueError("The field must contain complete tiles within the supported numeric range")
    keys = []
    left, top, right, bottom = meta["roi_nm"]
    n_source = negative_total = negative_pilot = negative_holdout = outside = 0
    for frame, x, y, photons in localization_chunks(path, fields):
        if not all(np.isfinite(v).all() for v in (frame, x, y, photons)):
            raise ValueError("Non-finite localization coordinates, frame or photon count")
        frame = frame - frame_base
        if np.any(frame != np.floor(frame)) or np.any(frame < 0) or np.any(frame >= meta["total_frames"]):
            raise ValueError("A localization frame falls outside the declared sequence or is non-integral")
        if kind == "picasso":
            x = (x + 0.5) * meta["pixel_size_nm"]
            y = (y + 0.5) * meta["pixel_size_nm"]
        if np.any((x < 0) | (x >= meta["camera_width_nm"]) | (y < 0) | (y >= meta["camera_height_nm"])):
            raise ValueError("Localization coordinates fall outside the declared field")
        # Unconstrained Gaussian fits can return negative amplitudes. The
        # declared nonnegative photon cutoff excludes those rows consistently.
        negative = photons < 0
        n_source += len(frame)
        negative_total += int(negative.sum())
        negative_pilot += int((negative & (frame < n_pilot)).sum())
        negative_holdout += int((negative & (frame >= n_pilot) & (frame < evaluated_until)).sum())
        in_roi = (x >= left) & (x < right) & (y >= top) & (y < bottom)
        outside += int((~in_roi).sum())
        keep = (frame < length) & (photons >= min_photons) & in_roi
        f = frame[keep].astype(np.int64)
        counts += np.bincount(f, minlength=length)
        xt = ((x[keep] - left) // tile_side_nm).astype(np.int64)
        yt = ((y[keep] - top) // tile_side_nm).astype(np.int64)
        tile_ok = (f < n_pilot) & (xt < nx) & (yt < ny)
        if np.any(tile_ok):
            keys.append((f[tile_ok] * ny + yt[tile_ok]) * nx + xt[tile_ok])
    per_tile = np.unique(np.concatenate(keys), return_counts=True)[1] if keys else np.array([], dtype=int)
    pilot_counts = counts[:n_pilot]
    single = int(np.count_nonzero(per_tile == 1))
    area = width * height / 1e6
    blocks = [int(pilot_counts[g].sum()) for g in np.array_split(np.arange(n_pilot), min(4, n_pilot))]
    return {
        "schema": "smlm-pilot/1", "acquisition": meta,
        "quality_control": {"source_localizations": n_source, "outside_selected_roi_in_source": outside,
                            "negative_photon_fits_in_source": negative_total,
                            "negative_photon_fits_in_pilot": negative_pilot,
                            "negative_photon_fits_in_evaluation_window": negative_holdout},
        "pilot": {
            "start_frame": 0, "frames": n_pilot, "min_photons": min_photons,
            "localizations": int(pilot_counts.sum()), "mean_per_frame": float(pilot_counts.mean()),
            "area_um2": area, "four_block_counts": blocks,
            "tiles": {
                "side_nm": tile_side_nm, "area_um2": tile_side_nm**2 / 1e6,
                "tiles_per_frame": n_tiles, "observations": n_tiles * n_pilot,
                "empty": n_tiles*n_pilot - len(per_tile), "single": single,
                "multiple": len(per_tile) - single, "localizations": int(per_tile.sum()),
                "complete_field_fraction": n_tiles * (tile_side_nm**2 / 1e6) / area,
            },
        },
        "evaluation": evaluate(counts, n_pilot),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Picasso *_locs.hdf5 or frame,x_nm,y_nm,photons CSV")
    parser.add_argument("--metadata", required=True, type=Path, help="JSON with acquisition conditions and dimensions")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--yaml", type=Path, help="Picasso multi-document YAML, defaults to the input stem")
    parser.add_argument("--min-photons", type=float, default=100)
    parser.add_argument("--frame-base", type=int, choices=(0, 1), default=0, help="CSV frame numbering")
    parser.add_argument("--tile-side-nm", type=float, default=250)
    parser.add_argument("--pilot-frames", type=int, default=200)
    args = parser.parse_args()
    try:
        meta = json.loads(args.metadata.read_text(encoding="utf-8-sig"))
        if not isinstance(meta, dict):
            raise ValueError("metadata must be a JSON object")
        report = summarize(args.input, meta, yaml_path=args.yaml, frame_base=args.frame_base,
                           min_photons=args.min_photons, tile_side_nm=args.tile_side_nm, pilot_frames=args.pilot_frames)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()

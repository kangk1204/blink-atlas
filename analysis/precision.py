"""Measure localization precision and compare strictly matched field summaries."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from analysis.pilot import acquisition, picasso_calibration, positive, positive_int


THRESHOLD = 20.
AXES = "both lateral 1-sigma uncertainties"
FIELD = "smlm-precision-field/1"
COMPARE = "smlm-precision-comparison/1"
STR_MATCH = ("sample_state", "cell_line", "target", "tag", "illumination_state",
             "density_calibration_id", "acquisition_id", "analysis_id", "roi_protocol_id",
             "density_origin", "density_unit")
NUM_MATCH = ("density_labels_per_um2", "excitation_nm", "exposure_ms", "frame_interval_ms")


def _columns(path):
    if path.suffix.lower() == ".csv":
        with path.open(encoding="utf-8-sig", newline="") as handle:
            return csv.DictReader(handle).fieldnames or []
    import h5py
    with h5py.File(path) as h:
        if "locs" not in h or not isinstance(h["locs"], h5py.Dataset) or not h["locs"].dtype.names:
            raise ValueError("Picasso requires a structured HDF5 /locs dataset")
        return list(h["locs"].dtype.names)


def _entry(value):
    if value is None or isinstance(value, str) and value.strip().lower() in ("", "na", "n/a", "null", "nan"):
        return np.nan
    try:
        return float(value)
    except (ValueError, TypeError):
        return np.inf


def _chunks(path, required, precision):
    columns = required + precision
    if path.suffix.lower() == ".csv":
        with path.open(encoding="utf-8-sig", newline="") as fh:
            reader = csv.DictReader(fh)
            while batch := list(itertools.islice(reader, 32_768)):
                yield [np.fromiter((_entry(r.get(k)) for r in batch), dtype=float, count=len(batch)) for k in columns]
    else:
        import h5py
        with h5py.File(path) as fh:
            locs = fh["locs"]
            for pos in range(0, len(locs), 65_536):
                batch = locs[pos:pos+65_536]
                yield [np.asarray(batch[k], dtype=float) if k in batch.dtype.names else np.full(len(batch), np.nan) for k in columns]


def _label(meta, key):
    val = meta.get(key)
    if val is None:
        return None
    if not isinstance(val, str) or not val.strip() or len(val.strip()) > 200:
        raise ValueError(f"{key} must be a short nonempty label")
    return val.strip()


def _metadata(path, settings, yaml_path):
    meta = dict(settings)
    picasso = path.suffix.lower() in (".h5", ".hdf5")
    if not picasso and path.suffix.lower() != ".csv":
        raise ValueError("Use CSV or a Picasso HDF5 table")
    native = {}
    if picasso:
        y = Path(yaml_path) if yaml_path else path.with_suffix(".yaml")
        if yaml_path and not y.is_file():
            raise ValueError("Requested Picasso YAML not found")
        native = picasso_calibration(y if y.is_file() else path)
    if "roi_nm" not in meta and not native.get("picasso_rois_px"):
        raise ValueError('Declare roi_nm as [left,top,right,bottom] in nm or the string "full_frame"')
    if meta.get("roi_nm") == "full_frame":
        if native.get("picasso_rois_px"):
            raise ValueError("A Picasso localization ROI is present; choose that ROI explicitly")
        del meta["roi_nm"]
    if meta.get("position_unit", "px" if picasso else "nm") != ("px" if picasso else "nm"):
        raise ValueError("The input position unit conflicts with the selected file format")
    acq = acquisition(meta, native, "picasso" if picasso else "csv")
    for key in ("prep_id", "field_id", "illumination_state", "density_calibration_id", "density_origin",
                "acquisition_id", "analysis_id", "roi_protocol_id", "timing_method"):
        acq[key] = _label(meta, key)
    acq["density_labels_per_um2"] = (positive(meta["density_labels_per_um2"], "density_labels_per_um2")
                                      if meta.get("density_labels_per_um2") is not None else None)
    acq["density_unit"] = meta.get("density_unit", "labels/um2")
    if acq["density_unit"] not in ("labels/um2", "labels/µm²"):
        raise ValueError("Density must use calibrated labels/um2, separate from labeling reagent concentration")
    acq["density_unit"] = "labels/um2"
    raw_start = meta.get("start_frame", 0)
    if isinstance(raw_start, bool) or not isinstance(raw_start, int) or raw_start < 0:
        raise ValueError("start_frame must be a nonnegative integer in zero-based acquisition coordinates")
    n = positive_int(meta.get("window_frames", acq["total_frames"] - raw_start), "window_frames")
    if raw_start + n > acq["total_frames"]:
        raise ValueError("The window extends beyond the acquired frame count")
    acq.update(start_frame=raw_start, window_frames=n)
    return acq, picasso


def field(path, settings, *, yaml_path=None, frame_base=0, scalar_name=None, scalar_column="precision_nm"):
    """Count fits in a real acquisition window; no localization-level rows are returned."""
    path = Path(path)
    acq, is_picasso = _metadata(path, settings, yaml_path)
    if frame_base not in (0, 1) or (is_picasso and frame_base != 0):
        raise ValueError("Use zero-based Picasso frames, or specify zero/one-based CSV frames")
    placement = ["frame", "x", "y"] if is_picasso else ["frame", "x_nm", "y_nm"]
    cols = _columns(path)
    if not set(placement).issubset(cols) or len(cols) != len(set(cols)):
        raise ValueError("Input requires distinct frame and lateral coordinate columns")
    if scalar_name is not None:
        if is_picasso:
            raise ValueError("Picasso lpx/lpy use per-axis uncertainties; separate scalar data require CSV")
        if not isinstance(scalar_name, str) or not scalar_name.strip() or len(scalar_name) > 200:
            raise ValueError("A named scalar precision definition is required")
        if scalar_column not in cols or scalar_column in placement:
            raise ValueError("Scalar precision column is missing or conflicts with coordinates")
        precision = [scalar_column]
        definition = {"kind": "scalar", "name": scalar_name.strip(), "unit": "nm", "threshold_nm": THRESHOLD}
    else:
        precision = ["lpx", "lpy"] if is_picasso else ["sigma_x_nm", "sigma_y_nm"]
        if not is_picasso and scalar_column in cols and not set(precision).issubset(cols):
            raise ValueError("Declare --scalar-name for a scalar uncertainty in nm")
        definition = {"kind": "axes", "name": AXES, "unit": "nm", "threshold_nm": THRESHOLD}
    schema_present = set(precision).issubset(cols)
    counts = {k:0 for k in ("source", "outside_window", "outside_roi", "invalid_geometry",
                              "roi_valid_window", "missing_precision", "invalid_precision",
                              "valid_precision", "qualified")}
    start, n = acq["start_frame"], acq["window_frames"]
    observed = set()
    left, top, right, bottom = acq["roi_nm"]
    for parts in _chunks(path, placement, precision):
        raw_frame, x, y, *sigmas = parts
        counts["source"] += len(raw_frame)
        frame = raw_frame - frame_base
        if (not np.isfinite(frame).all() or np.any(frame != np.floor(frame)) or
                (len(frame) and (frame.min() < 0 or frame.max() >= acq["total_frames"]))):
            raise ValueError("A localization frame is incompatible with the acquisition or frame base")
        if is_picasso:
            scale = acq["pixel_size_nm"]
            x, y = (x + 0.5) * scale, (y + 0.5) * scale
            sigmas = [v * scale for v in sigmas]
        window = (frame >= start) & (frame < start+n)
        xy = (np.isfinite(x) & np.isfinite(y) & (x>=0) & (x<acq["camera_width_nm"])
              & (y>=0) & (y<acq["camera_height_nm"]))
        selected = window & xy & (x>=left) & (x<right) & (y>=top) & (y<bottom)
        counts["outside_window"] += int((~window).sum())
        counts["invalid_geometry"] += int((window & ~xy).sum())
        counts["outside_roi"] += int((window & xy & ~selected).sum())
        counts["roi_valid_window"] += int(selected.sum())
        observed.update(map(int, np.unique(frame[selected])))
        uncertainty = np.vstack([v[selected] for v in sigmas])
        invalid = np.any(np.isinf(uncertainty) | (uncertainty <= 0), axis=0)
        missing = np.any(np.isnan(uncertainty), axis=0) & ~invalid
        valid = ~invalid & ~missing
        counts["invalid_precision"] += int(invalid.sum())
        counts["missing_precision"] += int(missing.sum())
        counts["valid_precision"] += int(valid.sum())
        counts["qualified"] += int((valid & np.all(uncertainty <= THRESHOLD, axis=0)).sum())
    area = acq["width_nm"] * acq["height_nm"] / 1e6
    dt = acq["frame_interval_ms"]
    # Wall time from the first shutter opening through the last shutter close.
    elapsed = ((n-1) * dt + acq["exposure_ms"]) / 1e3 if dt is not None else None
    if n == 1:
        elapsed = acq["exposure_ms"] / 1e3
    valid = counts["valid_precision"]
    status = "measured"
    if not schema_present or (counts["roi_valid_window"] > 0 and not valid):
        status = "unavailable"
    elif counts["invalid_geometry"] or counts["missing_precision"] or counts["invalid_precision"]:
        status = "partial"
    qualified = counts["qualified"] if status != "unavailable" else None
    counts["qualified"] = qualified
    with path.open("rb") as fh:
        table_sha = hashlib.file_digest(fh, "sha256").hexdigest()
    return {
        "schema": FIELD, "acquisition": acq, "definition": definition, "counts": counts,
        "source": {"format": "picasso" if is_picasso else "csv", "table_sha256": table_sha,
                   "precision_columns_present": schema_present},
        "denominators": {"frames": n, "empty_frames": n-len(observed), "area_um2": area, "elapsed_s": elapsed},
        "result": {"status": status, "fraction": qualified/valid if qualified is not None and valid else None,
                   "per_um2": qualified/area if qualified is not None else None,
                   "per_um2_s": qualified/(area*elapsed) if qualified is not None and elapsed else None},
    }


def validate_field(f):
    """Validate portable counts and denominators before aggregating an imported summary."""
    def real(v):
        return type(v) in (int,float) and math.isfinite(v)
    def count(v):
        return type(v) is int and 0 <= v <= 2**53-1
    def eq(v,w):
        return v is None if w is None else real(v) and math.isclose(v,w,rel_tol=1e-8,abs_tol=1e-10)
    def require(condition, reason):
        if not condition:
            raise ValueError(reason)
    require(isinstance(f,dict) and f.get("schema") == FIELD, "Expected a smlm-precision-field/1 summary")
    sections = ("acquisition","definition","counts","denominators","result","source")
    require(all(isinstance(f.get(k),dict) for k in sections),"Incomplete field summary")
    a,d,c,v,r,src = (f[k] for k in sections)
    sha = src.get("table_sha256")
    has_columns = src.get("precision_columns_present")
    require(src.get("format") in ("csv","picasso") and isinstance(sha,str) and len(sha)==64
            and all(s in "0123456789abcdef" for s in sha) and type(has_columns) is bool, "Invalid localization source")
    require(d.get("threshold_nm")==THRESHOLD and d.get("unit")=="nm" and
            (d.get("kind")=="axes" and d.get("name")==AXES or
             d.get("kind")=="scalar" and isinstance(d.get("name"),str) and 0<len(d["name"].strip())<=200),
            "Unsupported precision definition")
    for k in ("width_nm","height_nm","camera_width_nm","camera_height_nm","exposure_ms"):
        require(real(a.get(k)) and a[k]>0, "Invalid acquisition " + k)
    for k in ("frame_interval_ms","excitation_nm","pixel_size_nm","density_labels_per_um2"):
        require(a.get(k) is None or real(a[k]) and a[k]>0, "Invalid acquisition " + k)
    for k in ("sample_state","cell_line","target","tag","reference_dye","prep_id","field_id","illumination_state",
              "density_calibration_id","density_origin","acquisition_id","analysis_id","roi_protocol_id","timing_method"):
        require(a.get(k) is None or isinstance(a[k],str) and 0 < len(a[k].strip())<=200, "Invalid acquisition label " + k)
    require(a.get("density_unit")=="labels/um2", "Invalid density unit")
    if src["format"]=="picasso":
        require(real(a.get("pixel_size_nm")) and a["pixel_size_nm"]>0, "Missing Picasso pixel calibration")
    require(all(count(a.get(k)) for k in ("total_frames","start_frame","window_frames"))
            and a["window_frames"]>0 and a["start_frame"]+a["window_frames"]<=a["total_frames"], "Invalid acquisition frame window")
    roi = a.get("roi_nm")
    require(isinstance(roi,list) and len(roi)==4 and all(real(i) for i in roi), "Invalid physical ROI")
    require(0<=roi[0]<roi[2]<=a["camera_width_nm"] and 0<=roi[1]<roi[3]<=a["camera_height_nm"]
            and eq(a["width_nm"],roi[2]-roi[0]) and eq(a["height_nm"],roi[3]-roi[1]), "Inconsistent physical ROI")
    to_count = ("source","outside_window","outside_roi","invalid_geometry","roi_valid_window","missing_precision",
                "invalid_precision","valid_precision")
    require(all(count(c.get(k)) for k in to_count), "Invalid localization counts")
    require(c["source"]==sum(c[k] for k in ("outside_window","outside_roi","invalid_geometry","roi_valid_window"))
            and c["roi_valid_window"]==sum(c[k] for k in ("missing_precision","invalid_precision","valid_precision")),
            "Inconsistent localization totals")
    n, area = a["window_frames"], a["width_nm"]*a["height_nm"]/1e6
    elapsed = a["exposure_ms"]/1000 if n==1 else (
        ((n-1)*a["frame_interval_ms"]+a["exposure_ms"])/1000 if a.get("frame_interval_ms") else None)
    require(v.get("frames")==n and count(v.get("empty_frames")) and 0<=n-v["empty_frames"]<=min(n,c["roi_valid_window"])
            and eq(v.get("area_um2"),area) and eq(v.get("elapsed_s"),elapsed), "Invalid area/time or empty-frame denominators")
    valid = c["valid_precision"]
    status = "measured"
    if not has_columns or (c["roi_valid_window"] and not valid):
        status = "unavailable"
    elif c["invalid_geometry"] or c["invalid_precision"] or c["missing_precision"]:
        status = "partial"
    q = c.get("qualified")
    require(r.get("status")==status and (q is None if status=="unavailable" else count(q) and q<=valid),
            "Invalid precision coverage or qualified counts")
    require(eq(r.get("fraction"), q/valid if q is not None and valid else None)
            and eq(r.get("per_um2"), q/area if q is not None else None)
            and eq(r.get("per_um2_s"), q/(area*elapsed) if q is not None and elapsed else None),
            "Invalid precision rates")
    return f


def _pair_stats(group):
    area = sum(f["denominators"]["area_um2"] for f in group)
    exposure_area = sum(f["denominators"]["area_um2"]*f["denominators"]["elapsed_s"] for f in group)
    q = sum(f["counts"]["qualified"] for f in group)
    return {"n_fields": len(group), "qualified": q, "area_um2": area, "area_time_um2_s": exposure_area,
            "per_um2": q/area, "per_um2_s": q/exposure_area}


def compare(fields, dyes):
    """Pair preparations; normalize each dye by summed field area and area-time."""
    if len(dyes) != 2 or dyes[0].casefold() == dyes[1].casefold():
        raise ValueError("Choose two distinct dyes (numerator first)")
    result = {"schema": COMPARE, "dyes": list(dyes), "status": "unavailable", "reasons": [],
              "preparations": [], "n_fields": len(fields), "definition": None}
    if not fields:
        result["reasons"].append("No field summaries supplied.")
        return result
    reasons, seen, groups = set(), set(), defaultdict(lambda: defaultdict(list))
    reference = fields[0]
    for f in fields:
        validate_field(f)
        a, d, c, v, r = (f[k] for k in ["acquisition","definition","counts","denominators","result"])
        prep, field_id, dye_key = a.get("prep_id"), a.get("field_id"), str(a.get("reference_dye", "")).casefold()
        dye = next((x for x in dyes if x.casefold()==dye_key), None)
        if dye is None:
            reasons.add("Unexpected dye in the supplied fields.")
        if r.get("status") != "measured":
            reasons.add("Each field needs complete precision and geometry (a verified zero is allowed).")
        if not prep or not field_id:
            reasons.add("Both prep_id and field_id are required for independent pairing.")
        if (str(a.get("sample_state","")).casefold()!="fixed" or
            str(a.get("cell_line","")).casefold().replace("-","")!="cos7" or
            str(a.get("target","")).casefold()!="h2b" or str(a.get("tag","")).casefold()!="halotag"):
            reasons.add("The paired endpoint requires fixed COS-7 H2B–HaloTag fields.")
        if (prep, dye, field_id) in seen:
            reasons.add("Duplicate prep/dye/field identity.")
        seen.add((prep, dye, field_id))
        if prep and dye:
            groups[prep][dye].append(f)
        for k in STR_MATCH:
            if not isinstance(a.get(k), str) or a[k].strip().casefold() in {"","unknown","unavailable"}:
                reasons.add("Missing matched metadata: " + k + ".")
            elif str(reference["acquisition"].get(k, "")).strip().casefold() != a[k].strip().casefold():
                reasons.add("Incompatible " + k + ".")
        if str(a.get("density_origin", "")).casefold() != "measured" or a.get("density_unit") != "labels/um2":
            reasons.add("Actual independently calibrated labeling density in labels/um2 is required.")
        for k in NUM_MATCH + ("start_frame", "window_frames"):
            x, y = a.get(k), reference["acquisition"].get(k)
            if not isinstance(x,(int,float)) or isinstance(x,bool) or not math.isfinite(x) or (x <= 0 if k != "start_frame" else x < 0):
                reasons.add("Missing matched metadata: " + k + ".")
            elif not isinstance(y,(int,float)) or not math.isclose(x, y, rel_tol=1e-9, abs_tol=1e-9):
                reasons.add("Incompatible " + k + ".")
        if d != reference["definition"] or d.get("unit") != "nm" or d.get("threshold_nm") != THRESHOLD:
            reasons.add("Precision definition, units or threshold differ.")
        if not (isinstance(v.get("area_um2"), (int,float)) and v["area_um2"]>0
                and isinstance(v.get("elapsed_s"), (int,float)) and v["elapsed_s"]>0):
            reasons.add("A calibrated ROI area and acquisition duration are required.")
        if (r.get("status") == "measured" and
            (not isinstance(c.get("qualified"),int) or c["qualified"]<0 or c["qualified"]>c["valid_precision"])):
            raise ValueError("Invalid qualification counts in field summary")
    for prep, group in groups.items():
        if set(group) != set(dyes):
            reasons.add(f"Preparation {prep} lacks one of the two dyes.")
    if reasons:
        result["reasons"] = sorted(reasons)
        return result
    ar = reference["acquisition"]
    reference_model = (abs(ar["excitation_nm"]-647)<=1 and abs(ar["exposure_ms"]-25)<=.2 and abs(ar["frame_interval_ms"]-27)<=.2)
    result.update(status="matched", definition=reference["definition"],
                  model_domain="reference_acquisition" if reference_model else "outside_reference_acquisition",
                  matched_context={k:ar[k] for k in STR_MATCH + NUM_MATCH + ("start_frame","window_frames")})
    for prep in sorted(groups):
        num, den = (_pair_stats(groups[prep][dye]) for dye in dyes)
        result["preparations"].append({"prep_id":prep, "numerator":num, "denominator":den,
            "difference_per_um2":num["per_um2"]-den["per_um2"],
            "difference_per_um2_s":num["per_um2_s"]-den["per_um2_s"],
            "ratio_per_um2_s":num["per_um2_s"]/den["per_um2_s"] if den["per_um2_s"]>0 else None})
    result["mean_paired_difference_per_um2_s"] = float(np.mean([p["difference_per_um2_s"] for p in result["preparations"]]))
    return result


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    fo = subs.add_parser("field", help="Export one field's counts and denominators")
    fo.add_argument("input",type=Path); fo.add_argument("--metadata",type=Path,required=True)
    fo.add_argument("--yaml",type=Path); fo.add_argument("--frame-base",type=int,choices=[0,1],default=0)
    fo.add_argument("--scalar-name"); fo.add_argument("--scalar-column",default="precision_nm")
    co = subs.add_parser("compare", help="Numerator and denominator compared within each preparation")
    co.add_argument("fields",nargs="+",type=Path); co.add_argument("--dyes",nargs=2,required=True)
    pub = subs.add_parser("public", help="Reproduce source-locked public evidence")
    pub.add_argument("--cache",type=Path,required=True); pub.add_argument("--offline",action="store_true")
    for s in [fo,co,pub]: s.add_argument("--output",type=Path,required=True)
    a = parser.parse_args()
    try:
        if a.command=="field":
            data = field(a.input,json.loads(a.metadata.read_text(encoding="utf-8-sig")),yaml_path=a.yaml,
                         frame_base=a.frame_base, scalar_name=a.scalar_name, scalar_column=a.scalar_column)
        elif a.command=="compare":
            data = compare([json.loads(p.read_text(encoding="utf-8-sig")) for p in a.fields],a.dyes)
        else:
            from analysis.precision_sources import public_evidence
            data = public_evidence(a.cache,a.offline)
        save(a.output,data)
    except (OSError, ValueError, TypeError, KeyError) as exc:
        parser.error(str(exc))


if __name__=="__main__":
    main()

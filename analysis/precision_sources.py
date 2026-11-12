"""Fetch the Walterspiel field tables and a short timing prefix from public ZIPs."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
import json
import logging
from pathlib import Path
import struct
import zipfile
import zlib

from analysis.s21 import RECORD_URL, TailReader, request_bytes


ARCHIVES = {
    "S21": ("500a2ce1df2c388228c150db0559edd5", "HaloTag", "JF630b"),
    "S19": ("b275e87ac2baca18b55b8c401c454a88", "psHaloTag1a", "JF630b"),
    "S20": ("129ad804fc502fafbdb301784cd5aba8", "psHaloTag1a", "JF635b"),
    "6": ("595837faf303825886c3ffcb110c4901", "psHaloTag1a", "JF630b"),
}
MEMBER_SET_SHA256 = "6e178741ba22a0e05ff7e6f269c0c0d4505643d423302fc9c25fa2ef3f266f06"


def _member_bytes(url, member, prefix=None):
    header = request_bytes(url, member.header_offset, member.header_offset + 30)
    if header[:4] != b"PK\x03\x04":
        raise ValueError("Bad local ZIP header")
    name_len, extra_len = struct.unpack_from("<HH", header, 26)
    begin = member.header_offset + 30 + name_len + extra_len
    if member.compress_type not in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED):
        raise ValueError("Unsupported ZIP method")
    # Each requested chunk is independently retryable.
    decompressor = zlib.decompressobj(-15) if member.compress_type == zipfile.ZIP_DEFLATED else None
    limit = min(prefix or member.file_size, member.file_size)
    with io.BytesIO() as out:
        pos = 0
        while out.tell() < limit:
            end = min(pos + 4 * 1024**2, member.compress_size)
            if end <= pos:
                raise ValueError("Truncated member")
            data = request_bytes(url, begin + pos, begin + end)
            left = limit - out.tell()
            out.write(decompressor.decompress(data, left) if decompressor else data[:left])
            pos = end
        raw = out.getvalue()
    if limit == member.file_size and zlib.crc32(raw) != member.CRC:
        raise ValueError(f"CRC verification failed for {member.filename}")
    return raw


TIMING_FRAMES = 40


def _timing(data):
    """Use measured Micro-Manager ElapsedTime-ms, including empty camera frames."""
    import tifffile

    logger = logging.getLogger("tifffile")
    def ignore_trailing_metadata(record):
        return "invalid value offset" not in record.getMessage()
    logger.addFilter(ignore_trailing_metadata)
    try:
        with tifffile.TiffFile(io.BytesIO(data), is_ome=False, is_micromanager=False) as tif:
            times = []
            for i in range(TIMING_FRAMES):
                tag = tif.pages[i].tags[51123].value
                if int(tag["FrameIndex"]) != i:
                    raise ValueError("Unexpected TIFF FrameIndex")
                times.append(float(tag["ElapsedTime-ms"]))
        if not all(a < b for a, b in zip(times, times[1:])):
            raise ValueError("Non-increasing camera timestamps")
        return (times[-1] - times[0]) / (len(times) - 1)
    finally:
        logger.removeFilter(ignore_trailing_metadata)


def _illumination(name):
    n = name.lower()
    stage = "before" if "_before" in n else "during" if "_during" in n else "recovery" if "_after" in n else None
    if not stage:
        raise ValueError("Unknown illumination state for " + name)
    return stage


def _collect(args):
    figure, key, member, meta, tif, source, cache = args
    archive_md5, expected_tag, dye = ARCHIVES[figure]
    rec = dict(figure=figure, recording=key, reference_dye=dye, tag=expected_tag,
               illumination_state=_illumination(key), archive=source["key"], archive_md5=archive_md5,
               archive_bytes=source["size"], member=member.filename, member_bytes=member.file_size,
               member_crc32=f"{member.CRC:08x}", yaml_member=meta.filename, tiff_member=tif.filename)
    actual_tag = "psHaloTag1a" if "_psht_" in key.lower() else "HaloTag" if "_ht_" in key.lower() else None
    if actual_tag != expected_tag:
        rec.update(status="quarantined", reason=f"Archive and authors' record specify {expected_tag}; the localization filename encodes {actual_tag}.")
        return rec
    folder = cache / figure / key
    folder.mkdir(parents=True, exist_ok=True)
    url = source["links"]["self"]
    for obj, filename in [(member, "locs.hdf5"), (meta, "locs.yaml")]:
        path = folder / filename
        if path.exists():
            raw = path.read_bytes()
            if len(raw) != obj.file_size or zlib.crc32(raw) != obj.CRC:
                raise ValueError(f"Corrupt cached data: {path}")
        else:
            raw = _member_bytes(url, obj)
            path.write_bytes(raw)
        rec[filename.split(".")[-1] + "_sha256"] = hashlib.sha256(raw).hexdigest()
    # Use the same 40-frame timing calibration as the S21 pilot. Four-frame
    # estimates were susceptible to sub-ms Micro-Manager timestamp jitter.
    timing = folder / "timing.json"
    t = json.loads(timing.read_text()) if timing.exists() else {}
    if t.get("timing_frames") != TIMING_FRAMES:
        from analysis.pilot import picasso_calibration
        c = picasso_calibration(folder / "locs.yaml")
        raw_pixels = int(c["width_px"] * c["height_px"]) * 2
        prefix = TIMING_FRAMES * (raw_pixels + 30 * 1024) + 512 * 1024
        raw = _member_bytes(url, tif, prefix=prefix)
        try:
            interval = _timing(raw)
        except (ValueError, IndexError, KeyError, TypeError, OSError):
            interval = None
        t = {"frame_interval_ms": interval, "timing_frames":TIMING_FRAMES, "timing_prefix_sha256": hashlib.sha256(raw).hexdigest(),
             "timing_prefix_uncompressed_bytes": len(raw),
             "timing_method": "Mean of 39 inter-frame intervals from the first 40 Micro-Manager TIFF timestamps" if interval else "TIFF timing unavailable"}
        timing.write_text(json.dumps(t, indent=2) + "\n")
    rec.update(t, status="qc_only")
    return rec


def fetch_manifest(cache):
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    record = json.loads(request_bytes(RECORD_URL))
    jobs = []
    for source in record["files"]:
        md5 = source["checksum"].split(":")[-1]
        figure = next((k for k,v in ARCHIVES.items() if v[0] == md5), None)
        if not figure:
            raise ValueError("Unexpected upstream archive version")
        size = source["size"]
        tail = request_bytes(source["links"]["self"], size-65557, size)
        with zipfile.ZipFile(TailReader(tail, size)) as z:
            members = {m.filename: m for m in z.infolist()}
        for name, loc in members.items():
            if not name.endswith("_locs.hdf5") or name.startswith("__MACOSX"):
                continue
            y = members[name[:-5] + ".yaml"]
            tif = members[name.removesuffix("_locs.hdf5") + ".tif"]
            key = Path(name).name.split("_MMStack")[0]
            jobs.append((figure, key, loc, y, tif, source, cache))
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(_collect, jobs))
    results.sort(key=lambda r: (r["tag"], r["reference_dye"], {"before":0,"during":1,"recovery":2}[r["illumination_state"]],
                                r["figure"], r["recording"]))
    (cache / "manifest.json").write_text(json.dumps(results,indent=2) + "\n", encoding="utf-8")
    return results


def public_evidence(cache, offline=False):
    """Rebuild the fixed-cell histogram and explicitly separate live-cell QC."""
    import h5py
    import numpy as np
    from analysis.pilot import evaluate, picasso_calibration
    from analysis.precision import field
    from validation.h2b_precision.reproduce import SOURCE, SHA256, h2b_precision

    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / "manifest.json"
    manifest = json.loads(path.read_text()) if offline else fetch_manifest(cache)
    pins = sorted((r["archive_md5"],r["member"],r["member_crc32"],r.get("hdf5_sha256"),r.get("yaml_sha256"))
                  for r in manifest)
    encoded = json.dumps(pins,separators=(",",":"),ensure_ascii=False).encode("utf-8")
    if hashlib.sha256(encoded).hexdigest() != MEMBER_SET_SHA256:
        raise ValueError("Walterspiel localization source manifest differs from the audited set")
    wb = cache / "holland-extended-data-5.xlsx"
    if not wb.exists() and not offline:
        wb.write_bytes(request_bytes(SOURCE))
    source_bytes = wb.read_bytes()
    if hashlib.sha256(source_bytes).hexdigest() != SHA256:
        raise ValueError("Fixed-cell source workbook hash mismatch")
    histogram = h2b_precision(source_bytes)
    histogram.update(definition={"kind":"published_scalar","name":"Author-reported scalar localization precision",
                                 "unit":"nm","threshold_nm":20},
                     denominator="Displayed histogram frequency", per_um2=None, per_um2_s=None)
    root = Path(__file__).resolve().parents[1]
    legacy = json.loads((root / "docs/data/pilot.json").read_text(encoding="utf-8"))
    confirmed, s21, independent = [], None, None
    for rec in manifest:
        if rec["status"] == "quarantined":
            continue
        folder = cache / rec["figure"] / rec["recording"]
        table, y = folder / "locs.hdf5", folder / "locs.yaml"
        for file, key in [(table,"hdf5_sha256"), (y,"yaml_sha256")]:
            if hashlib.sha256(file.read_bytes()).hexdigest() != rec[key]:
                raise ValueError("Source member hash mismatch for " + rec["recording"])
        primary = rec["figure"]=="S21" and rec["illumination_state"]=="before"
        calib = picasso_calibration(y)
        interval = rec["frame_interval_ms"]
        method = rec["timing_method"]
        if primary:
            if rec["hdf5_sha256"] != legacy["source"]["members"]["hdf5"]["sha256"]:
                raise ValueError("S21 does not match the original pilot fit table")
            interval = legacy["acquisition"]["frame_interval_ms"]
            method = legacy["media"]["interval_source"]
        elif interval and interval < .998 * calib["exposure_ms"]:
            # First-frame timestamp jitter can exceed the difference between
            # cycle and exposure. Preserve the timestamps in the manifest,
            # withhold the full-field time rate rather than repairing the clock.
            interval = None
            method = "Unresolved clock: short timing prefix implies a cycle shorter than the shutter."
        settings = dict(sample_state="live", cell_line="U2OS", target="H2B", tag=rec["tag"], reference_dye=rec["reference_dye"],
                        excitation_nm=640, roi_nm="full_frame", field_id=rec["recording"],
                        illumination_state=rec["illumination_state"], frame_interval_ms=interval, timing_method=method,
                        acquisition_id="Walterspiel SMLM 640 nm", roi_protocol_id="full author camera field",
                        analysis_id="author Picasso LQ Gaussian fits")
        summary = field(table, settings, yaml_path=y)
        summary["source"] = {**summary["source"], **rec, "doi":"10.5281/zenodo.16941576",
                             "url":"https://zenodo.org/records/16941576","license":"CC BY 4.0"}
        with h5py.File(table) as hf:
            locs = hf["locs"]
            counts = np.zeros(1000, int)
            for start in range(0,len(locs),65536):
                rows = locs[start:start+65536]
                keep = (rows["frame"]<1000) & (rows["photons"]>=0)
                counts += np.bincount(rows["frame"][keep].astype(int),minlength=1000)
            ev = evaluate(counts,200)
            summary["qc_forecast"] = {k:ev[k] for k in ["frames","block_frames","model_mae","baseline_mae","model_improved","unit"]}
            summary["qc_forecast"].update(pilot_frames=200, threshold_photons=0, evaluation_start_frame=200)
            if primary:
                # Independent, whole-table HDF5 count, separate from chunked
                # coordinate/ROI importer. Threshold comparison is on raw pixels.
                data = locs[:]
                x, z = data["lpx"].astype("float64"), data["lpy"].astype("float64")
                valid = np.isfinite(x) & np.isfinite(z) & (x>0) & (z>0)
                qualified = valid & (x <= 20/calib["pixel_size_nm"]) & (z <= 20/calib["pixel_size_nm"])
                independent = {"source_rows":len(data), "valid_precision":int(valid.sum()),"qualified":int(qualified.sum()),
                               "rule":"Both independent pixel uncertainties are finite, positive and ≤20/106 px."}
                if (len(data),int(valid.sum()),int(qualified.sum())) != (214794,214794,110694) or int(qualified.sum()) != summary["counts"]["qualified"]:
                    raise ValueError("Independent S21 count disagrees with the source audit")
        if primary:
            s21 = summary
        else:
            confirmed.append(summary)
    if s21 is None or independent is None:
        raise ValueError("S21 primary control unavailable")
    outcomes = [s["qc_forecast"] for s in [s21,*confirmed]]
    qc = {"recordings":len(outcomes), "trend_lower_mae":sum(r["model_mae"] < r["baseline_mae"] for r in outcomes),
          "equal_mae":sum(r["model_mae"] == r["baseline_mae"] for r in outcomes),
          "trend_higher_mae":sum(r["model_mae"] > r["baseline_mae"] for r in outcomes),
          "unit":"one recording; phases from a field can share a biological specimen",
          "outcome":"Forecast of localization counts per 100 frames"}
    return {"schema":"smlm-precision-evidence/1", "holland":histogram, "s21":s21, "independent_check":independent,
            "auxiliary":confirmed, "manifest":manifest,
            "qc_overview":qc,
            "source":{"doi":"10.5281/zenodo.16941576","url":"https://zenodo.org/records/16941576",
                      "paper":"https://doi.org/10.1002/anie.202424955","license":"CC BY 4.0"},
            "matched_fixed_comparison":{"status":"unavailable","dyes":["JF639b","JF630b"],
                "reasons":["No paired, calibrated-density fixed H2B fields with both dyes and a common precision definition."]}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    data = fetch_manifest(args.cache)
    print(len(data), "recordings indexed,", sum(r["status"] == "quarantined" for r in data), "quarantined")

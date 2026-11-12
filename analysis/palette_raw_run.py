"""Analyze the 435 Fig. 2/3 HDF5 members of Steen's deposit, one replicate at a time.

Use a verified extracted archive directory. An optional quick-open RAR index
enables streaming analysis as complete members arrive; otherwise discover the
full extraction after unarchiving. Each completed case gets an atomic cache.
"""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path, PureWindowsPath
import re
import tempfile
import time

import h5py
import numpy as np
import yaml

try:
    from . import palette_raw as raw
except ImportError:
    import palette_raw as raw


SYSTEMS = {"fig_2": "origami", "fig_3": "npc"}


def read_cases(root, index):
    if index:
        document = json.loads(index.read_text(encoding="utf-8"))
        entries = document["files"] if isinstance(document, dict) else document
    else:
        entries = [
            {"path": str(x.relative_to(root).as_posix()), "bytes": x.stat().st_size}
            for x in root.rglob("*.hdf5")
        ]
    metadata_digests = {
        e["path"]: e["sha256"]
        for e in entries
        if e["path"].lower().endswith(".yaml") and e.get("sha256")
    }
    cases = defaultdict(dict)
    for entry in entries:
        relative = Path(entry["path"])
        if relative.suffix.lower() != ".hdf5":
            continue
        if "crc32" in entry:
            entry = {**entry, "crc": int(entry["crc32"], 16)}
        sidecar = relative.with_suffix(".yaml").as_posix()
        if sidecar in metadata_digests:
            entry = {**entry, "metadataSha256": metadata_digests[sidecar]}
        parts = relative.parts
        panel = next((p for p in parts if p in SYSTEMS), None)
        if not panel:
            continue
        i = parts.index(panel)
        dye, repeat_dir = parts[i + 1 : i + 3]
        n = int(re.fullmatch(r"n(\d+)(?:_fig3)?", repeat_dir).group(1))
        stem = relative.stem.lower()
        if panel == "fig_2":
            role = "target_sites" if "sbs" in stem else "target_field"
        else:
            prefix = (
                "reference_" if "reference_cy3b" in parts[i + 3 : -1] else "target_"
            )
            picked = "pick" in stem
            region = (
                "cytoplasm"
                if "cyto" in stem
                else "npc"
                if picked and ("nup" in stem or "npc" in stem)
                else "field"
            )
            role = prefix + region
        key = (SYSTEMS[panel], dye, n)
        if role in cases[key]:
            raise ValueError(f"Ambiguous role {role} in {key}")
        cases[key][role] = entry
    for key, files in cases.items():
        if key[0] == "origami":
            expected = {"target_field", "target_sites"}
        else:
            expected = {
                prefix + suffix
                for prefix in ("target_", "reference_")
                for suffix in ("field", "npc", "cytoplasm")
            }
            if key[1:] == ("cy3b", 4):
                expected = {x for x in expected if x.startswith("target_")}
        if set(files) != expected:
            raise ValueError(f"Expected {expected}, got {set(files)} in {key}")
    return cases


def file_info(root, role, entry):
    path = root / entry["path"]
    with h5py.File(path) as fh:
        ds = fh["locs"]
        info = {
            "role": role,
            "path": entry["path"],
            "rows": len(ds),
            "bytes": path.stat().st_size,
        }
        if "crc" in entry:
            info["archiveCrc32"] = f"{entry['crc']:08x}"
        if "sha256" in entry:
            with path.open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            if digest != entry["sha256"]:
                raise ValueError(
                    f"Input SHA-256 differs from archive manifest for {path}"
                )
            info["archiveSha256"] = digest
        if role.endswith(("npc", "cytoplasm")):
            groups = np.unique(ds.fields("group")[:])
            info["nonemptyPicks"] = int(len(groups))
    meta = path.with_suffix(".yaml")
    if "metadataSha256" in entry and not meta.exists():
        raise FileNotFoundError(f"Missing manifest-listed sidecar {meta}")
    if meta.exists():
        info["metadataPath"] = meta.relative_to(root).as_posix()
        info["metadataSha256"] = hashlib.sha256(meta.read_bytes()).hexdigest()
        if (
            entry.get("metadataSha256", info["metadataSha256"])
            != info["metadataSha256"]
        ):
            raise ValueError(
                f"Sidecar SHA-256 differs from archive manifest for {meta}"
            )
        with meta.open(encoding="utf8") as stream:
            docs = [x for x in yaml.safe_load_all(stream) if isinstance(x, dict)]
        for doc in docs:
            if doc.get("File") and doc.get("Frames"):
                info["originalMovie"] = PureWindowsPath(str(doc["File"])).name
            for k, key in (
                ("Frames", "metadataFrames"),
                ("Pick Diameter", "pickDiameterPixels"),
                ("Area (um^2)", "pickAreaUm2"),
                ("Pixelsize", "pixelSizeNm"),
                ("Box Size", "boxSizePixels"),
                ("Min. Net Gradient", "minimumNetGradient"),
            ):
                if k in doc:
                    info[key] = float(doc[k])
            if doc.get("Pick Shape"):
                info["pickShape"] = str(doc["Pick Shape"])
            if doc.get("Fit method"):
                info["fitMethod"] = str(doc["Fit method"])
            mm = doc.get("Micro-Manager Metadata") or {}
            if isinstance(mm, dict):
                if mm.get("Exposure-ms") is not None:
                    info["cameraExposureMs"] = float(mm["Exposure-ms"])
                if mm.get("Camera"):
                    info["camera"] = str(mm["Camera"])
    return info


def specificity(info):
    by = {x["role"]: x for x in info}
    roles = ["target_npc", "target_cytoplasm", "reference_npc", "reference_cytoplasm"]
    if not all(r in by for r in roles):
        return None, "No separate paired reference deposition."
    if all(
        "pickDiameterPixels" in by[r]
        and by[r].get("pickShape", "circle").lower() == "circle"
        for r in roles
    ):
        areas = [
            by[r]["nonemptyPicks"]
            * np.pi
            * (
                by[r]["pickDiameterPixels"]
                * by[r].get("pixelSizeNm", raw.PIXEL_NM)
                / 2000
            )
            ** 2
            for r in roles
        ]
        mode = "Published circular-pick formula with deposited diameters and nonempty group counts."
    elif all("pickAreaUm2" in by[r] for r in roles):
        areas = [by[r]["pickAreaUm2"] for r in roles]
        mode = "Deposited pick areas."
    else:
        # Multiplicative per-pick areas cancel in the four-way ratio when
        # corresponding selections have the same geometry across rounds.
        areas = [by[r]["nonemptyPicks"] for r in roles]
        mode = "Area factors cancel assuming corresponding transferred picks have equal diameters across rounds; counts divided by nonempty picks."
    if min(areas) <= 0 or any(by[r]["rows"] <= 0 for r in roles):
        return None, "Empty comparison region."
    a, b, c, d = [by[r]["rows"] / area for r, area in zip(roles, areas)]
    return float((a / b) / (c / d)), mode


def atomic_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=path.parent, suffix=".tmp", delete=False, mode="w", encoding="utf8"
    ) as tmp:
        json.dump(data, tmp, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
        tmp.write("\n")
        name = tmp.name
    os.replace(name, path)


def analyze(root, key, entries):
    system, dye, n = key
    infos = [file_info(root, role, e) for role, e in sorted(entries.items())]
    by = {a["role"]: a for a in infos}
    picks = by["target_sites" if system == "origami" else "target_npc"]
    data = raw.file_fields(
        root / picks["path"], ["group", "frame", "photons", "sx", "sy", "bg"]
    )
    frames = int(picks.get("metadataFrames", 20000 if system == "origami" else 30000))
    qc = []
    if data["frame"].max() >= frames:
        raise ValueError(
            f"Observed frame outside documented acquisition duration in {picks['path']}"
        )
    metrics, counts = raw.picked_metrics(data, frames, system == "origami")
    del data
    fields = {}
    for role in ("target_field", "reference_field"):
        if role not in by:
            continue
        try:
            nm, details = raw.field_precision(
                root / by[role]["path"], by[role].get("pixelSizeNm", raw.PIXEL_NM)
            )
            fields[role] = {"nenaNm": nm, **details}
            if role == "target_field":
                metrics["nenaNm"] = nm
        except (ValueError, RuntimeError, FloatingPointError) as e:
            qc.append(f"{role} NeNA unavailable: {e}")
    if system == "npc":
        metrics["specificity"], message = specificity(infos)
        qc.append(message)
    return {
        "system": system,
        "dye": dye,
        "repeat": n,
        "principal": n <= 3,
        "files": infos,
        "metrics": metrics,
        "counts": counts,
        "fieldFits": fields,
        "qc": qc,
    }


def identity(key, entries, root):
    h = hashlib.sha256()
    for name in (Path(__file__), Path(raw.__file__)):
        h.update(name.read_bytes())
    for e in sorted(entries.values(), key=lambda e: e["path"]):
        member = root / e["path"]
        stamp = member.stat().st_mtime_ns if member.is_file() else "missing"
        h.update(
            f"{e['path']}:{e['bytes']}:{e.get('crc', '')}:{e.get('sha256', '')}:{e.get('metadataSha256', '')}:{stamp}".encode()
        )
        y = member.with_suffix(".yaml")
        if y.is_file():
            h.update(y.read_bytes())
    return h.hexdigest()


def execute_case(root, key, entries, file, version):
    started = time.monotonic()
    result = analyze(root, key, entries)
    atomic_json(file, {"inputFingerprint": version, "result": result})
    return key, time.monotonic() - started, result["metrics"]


def make_summary(rows, cases):
    grouped = defaultdict(list)
    for row in rows:
        if row["principal"]:
            grouped[row["system"], row["dye"]].append(row)
    dyes = []
    for (system, dye), records in sorted(grouped.items()):
        records.sort(key=lambda record: record["repeat"])
        metrics = {}
        for k in records[0]["metrics"]:
            vals = [r["metrics"][k] for r in records if r["metrics"][k] is not None]
            metrics[k] = {
                "n": len(vals),
                "mean": float(np.mean(vals)) if vals else None,
                "sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else None,
            }
        dyes.append(
            {
                "system": system,
                "dye": dye,
                "repeatCount": len(records),
                "metrics": metrics,
            }
        )
    identities = defaultdict(list)
    for es in cases.values():
        for e in es.values():
            if "crc" in e:
                identities[(e["bytes"], e["crc"])].append(e["path"])
    fields = [x for row in rows for x in row["files"]]
    return {
        "schema": "palette-raw/1",
        "source": {
            "doi": "10.1038/s41592-024-02374-8",
            "archive": "https://zenodo.org/records/10960858",
            "archiveDoi": "10.5281/zenodo.10960858",
            "archiveMd5": "4d46f75d5cd23313925f044298c89ddf",
            "archiveFile": "dye_analysis.rar",
            "authorCode": "https://github.com/PhilippSteen/DNA-PAINT_analysis/tree/dad2df581be9",
            "picasso": "https://github.com/jungmannlab/picasso/blob/deab0c6e31c544fc56d1e4d0ccc8e00348dcae92/picasso/postprocess.py",
            "exposureS": 0.1,
            "pixelSizeNm": 130,
        },
        "coverage": {
            "complete": len(rows) == len(cases),
            "expectedCases": len(cases),
            "completedCases": len(rows),
            "principalCaseCount": sum(r["principal"] for r in rows),
            "metadataSidecarsRead": sum("metadataPath" in f for f in fields),
            "hdf5Sha256Verified": sum("archiveSha256" in f for f in fields),
            "totalFiles": len(fields),
            "filesBySystem": dict(
                Counter(
                    {
                        s: sum(len(r["files"]) for r in rows if r["system"] == s)
                        for s in ("origami", "npc")
                    }
                )
            ),
            "totalRowsAcrossOverlappingFiles": sum(f["rows"] for f in fields),
        },
        "metrics": {
            "centerPhotons": "photons per 100 ms, positive center frames",
            "eventPhotons": "mean finite summed photons per consecutive-frame binding event (raw-only diagnostic)",
            "centerSbr": "center-frame photons in centered Gaussian 1-pixel window / background per pixel",
            "brightTimeS": "mean of positive per-pick exponential CDF fits from the released 200-frame initialization, in seconds",
            "darkTimeS": "mean of positive per-pick inter-start CDF fits in seconds",
            "dropPct": "100 × (1 − final-fifth count / first-fifth count) on the selected target sites or NPCs",
            "destroyedPct": "apparent destroyed % after chance-tail correction",
            "specificity": "target NPC:cytoplasm density ratio divided by paired Cy3B reference ratio",
            "nenaNm": "historical Picasso next-frame fit s × camera pixel size in nm",
        },
        "limitations": [
            "Deposited files are fitted, selected localizations. Photon brightness and SBR depend on the published microscope, thresholds and site selection.",
            "Field, picks and cytoplasm are related subsets. Their row counts cannot be added to infer independent molecules.",
            "NeNA uses all same-group pairs in adjacent frames under one pixel and the pinned historical Picasso five-parameter fit. It does not specify split-image resolution.",
            "Cy3B n4_fig3 duplicates a reference field and has no separate reference round. It is reported separately from the three principal repetitions.",
            "Repeated reference files link measurements across dyes. Publication repeats are the aggregation units. Localizations and picks are not independent replicate experiments.",
            "The authors discarded measurements with poor Cy3B alignment precision, selected and refined NPC picks, and manually picked cytoplasm. Retrospective validation is conditional on those retained acquisitions and regions.",
            "NPC drop applies the released first-to-last fifth calculation to the deposited picked NPCs. The Nup96 methods do not specify a separate drop procedure.",
        ],
        "sharedArchiveMembers": [v for v in identities.values() if len(v) > 1],
        "replicates": sorted(rows, key=lambda r: (r["system"], r["dye"], r["repeat"])),
        "dyes": dyes,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--root", type=Path, required=True, help="Directory containing 39.dye_analysis"
    )
    p.add_argument(
        "--index",
        type=Path,
        help="Public palette manifest or staging RAR member index JSON",
    )
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--workers", type=int, default=3)
    p.add_argument(
        "--watch",
        action="store_true",
        help="Wait for all indexed members; use during archive staging.",
    )
    args = p.parse_args()
    cases = read_cases(args.root, args.index)
    if len(cases) != 109:
        raise ValueError(
            f"Expected 109 cases covering the full main panels, found {len(cases)}"
        )
    rows = {}
    for key, entries in cases.items():
        path = args.cache / f"{key[0]}-{key[1]}-{key[2]}.json"
        if path.exists():
            obj = json.loads(path.read_text())
            if obj["inputFingerprint"] == identity(key, entries, args.root):
                rows[key] = obj["result"]
    print(f"Resuming {len(rows)}/{len(cases)} completed cases.", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        while len(rows) < len(cases):
            futures = {}
            for key, es in cases.items():
                if key in rows:
                    continue
                if not all(
                    (args.root / e["path"]).is_file()
                    and (args.root / e["path"]).stat().st_size == e["bytes"]
                    for e in es.values()
                ):
                    continue
                file = args.cache / f"{key[0]}-{key[1]}-{key[2]}.json"
                v = identity(key, es, args.root)
                futures[pool.submit(execute_case, args.root, key, es, file, v)] = file
            for f in as_completed(futures):
                key, seconds, metrics = f.result()
                rows[key] = json.loads(futures[f].read_text())["result"]
                print(
                    f"{len(rows)}/{len(cases)} {key} {seconds:.1f}s {metrics}",
                    flush=True,
                )
            atomic_json(args.output, make_summary(list(rows.values()), cases))
            if not args.watch:
                break
            if not futures and len(rows) < len(cases):
                time.sleep(12)
    if len(rows) < len(cases):
        examples = sorted(set(cases) - set(rows))[:3]
        raise SystemExit(
            f"Only {len(rows)}/{len(cases)} cases have complete input members. "
            f"Missing or incomplete examples: {examples}. Finish extracting the archive or rerun with --watch."
        )
    atomic_json(args.output, make_summary(list(rows.values()), cases))


if __name__ == "__main__":
    main()

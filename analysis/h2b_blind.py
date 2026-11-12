"""Lock fixed-H2B forecasts before receiving observations, then score the same enrolled dyes."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from math import sqrt
import os
from pathlib import Path
from statistics import mean
from tempfile import NamedTemporaryFile

from analysis.analogue_benchmark import (
    METHOD_ID, PROTOCOL, holland_data, identity, number, predict, primary_verdict, score_predictions, summarize,
)
from analysis.probe_design import SCOPE, published


MODEL_REVISION = "8edc8bfdbe8ee51cf5d6f2907353f831f974507c"
REGISTRATION = "4c68d5a39e2559d05e4b8b8144f4fba5bf2341aa"
SCHEMA = "blink-atlas/h2b-blind-lock/1"
ROSTER = ("blind_id", "free_parent_pka")
MEASUREMENTS = ("blind_id", "preparation_id", "mean_on_ms", "mean_off_s", "source_locator")
PRIMARY = "matched free-parent titration"
TRANSFER = "literature buffer transfer"


def csv_rows(path, fields):
    with path.open(encoding="utf-8-sig", newline="") as file:
        rows = csv.DictReader(file)
        if rows.fieldnames is None or sorted(rows.fieldnames) != sorted(fields):
            raise ValueError(f"{path.name} must contain exactly {', '.join(fields)}")
        values = list(rows)
    if any(None in row or any(v is None for v in row.values()) for row in values):
        raise ValueError(f"Malformed row in {path.name}")
    return values


def label(value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError("An ID or source locator is missing")
    return value.strip()


def make_lock(study, roster, arm=PRIMARY, when=None):
    if not roster:
        raise ValueError("Enroll at least one real parent before locking a panel")
    catalog = holland_data()
    calibration = catalog["calibration"]
    seen = {identity(k) for d in catalog["dyes"] for k in (d["id"], d["ligand"])}
    dyes = []
    for record in roster:
        if set(record) != set(ROSTER):
            raise ValueError("The roster can contain only the blind ID and free-parent pKa")
        blind_id = label(record["blind_id"])
        key = identity(blind_id)
        if not key or key in seen:
            raise ValueError(f"Previously studied or repeated dye ID: {blind_id}")
        seen.add(key)
        raw = record["free_parent_pka"]
        if raw is None or (isinstance(raw, str) and not raw.strip()):
            dyes.append({"blindId": blind_id, "pKa": None, "status": "no exact monotonic pKa"})
            continue
        pka = number(float(raw)) if isinstance(raw, str) else number(raw)
        forecast = predict(calibration["pairs"], pka)
        if forecast["status"] == "scored":
            forecast["status"] = "forecast locked"
        dyes.append({"blindId": blind_id, "pKa": pka, **forecast})
    return {"schema": SCHEMA, "study": label(study), "arm": arm, "scope": SCOPE,
            "lockedUTC": when or datetime.now(timezone.utc).isoformat(timespec="microseconds"),
            "registrationRevision": REGISTRATION, "modelRevision": MODEL_REVISION, "candidate": METHOD_ID,
            "sourceWorkbookSha256": calibration["sourceSha256"], "trainingIds": [r["id"] for r in calibration["pairs"]],
            "plannedProtocol": PROTOCOL.copy(), "preparationsPerDye": 3, "dyes": dyes}


def public_lock():
    lock = make_lock("Cal664sb", [{"blind_id": "Cal664sb", "free_parent_pka": 3.2}],
                     TRANSFER, "2026-11-12T14:28:31+00:00")
    lock["source"] = {"article": "https://doi.org/10.1021/acscentsci.1c00670",
                      "figure": "https://pmc-oa-opendata.s3.amazonaws.com/PMC8393207.1/oc1c00670_0003.jpg",
                      "figureSha256": "90656c5ee2543aaf2de170398137a60dd7dd41c584365d433fcbab9f8d564e07",
                      "locator": "Fig. 3b, Cal664sb, pKcycl column. Fig. 3e, 2 μM in 0.2 M sodium phosphate at room temperature."}
    predictions = lock["dyes"][0]["predictions"]
    order = sorted(predictions, key=predictions.get)
    lock["equalErrorBoundaries"] = [{"models": [a, b], "ratio": published(sqrt(predictions[a]*predictions[b]))}
                                     for a, b in zip(order, order[1:])]
    return lock


def check_acquisition(assay, lock):
    if (lock.get("schema") != SCHEMA or lock.get("arm") not in (PRIMARY, TRANSFER)
            or lock.get("scope") != SCOPE or lock.get("plannedProtocol") != PROTOCOL
            or lock.get("modelRevision") != MODEL_REVISION or lock.get("preparationsPerDye") != 3):
        raise ValueError("Unrecognized prospective study or protocol")
    if not isinstance(assay, dict) or not isinstance(assay.get("actualProtocol"), dict):
        raise ValueError("Supply the actual acquisition protocol")
    if assay.get("study") != lock["study"]:
        raise ValueError("The acquisition study ID differs from the locked roster")
    flags = ["independentAcquisition", "identitiesReviewedAgainstHolland", "outcomesWithheldUntilLock",
             "illuminationAndDetectionMatched"]
    if lock["arm"] == PRIMARY:
        flags.append("freeParentTitrationsMatchHolland")
    missing = [k for k in flags if assay.get(k) is not True]
    differing = [k for k, v in lock["plannedProtocol"].items() if assay["actualProtocol"].get(k) != v]
    if missing or differing:
        raise ValueError("Unconfirmed or unmatched acquisition fields: " + ", ".join(missing + differing))
    label(assay.get("acquisitionRef"))


def score(lock, measurements, assay):
    check_acquisition(assay, lock)
    preparations = {d["blindId"]: {} for d in lock["dyes"]}
    for m in measurements:
        dye, prep = label(m["blind_id"]), label(m["preparation_id"])
        if dye not in preparations or prep in preparations[dye]:
            raise ValueError("Unknown dye or duplicate preparation: " + dye + " / " + prep)
        on, off = number(float(m["mean_on_ms"])), number(float(m["mean_off_s"]))
        if min(on, off) <= 0:
            raise ValueError("Mean ON (ms) and OFF (s) times must be positive")
        ratio = number(on / off / 1000)
        if ratio <= 0:
            raise ValueError("The measured ratio must be positive")
        preparations[dye][prep] = {"preparationId": prep, "meanOnMs": on, "meanOffS": off,
                                   "ratio": ratio, "sourceLocator": label(m["source_locator"])}
    rows, pending = [], 0
    for dye in lock["dyes"]:
        reps = list(preparations[dye["blindId"]].values())
        if len(reps) > lock["preparationsPerDye"]:
            raise ValueError("More preparations than registered for " + dye["blindId"])
        row = {**dye, "preparations": reps}
        if dye["status"] == "forecast locked":
            if len(reps) == lock["preparationsPerDye"]:
                observed = mean(r["ratio"] for r in reps)
                row.update(status="scored", observed=observed, **score_predictions(dye["predictions"], observed))
            else:
                row["status"] = "awaiting preparations"
                pending += 1
        rows.append(row)
    summary = summarize(rows, independent=True)
    if lock["arm"] == TRANSFER:
        verdict = "conditional test only"
    elif pending:
        verdict = "incomplete panel"
    else:
        verdict = primary_verdict(summary)
    return {"schema": "blink-atlas/h2b-blind-score/1", "study": lock["study"], "arm": lock["arm"],
            "acquisition": assay, "pendingDyes": pending, "verdict": verdict, "summary": summary, "rows": rows}


def save(path, value):
    raw = (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+"\n").encode()
    # Publish a complete file, without replacing a previous lock on a retry.
    with NamedTemporaryFile(dir=path.parent, delete=False) as file:
        staging = Path(file.name)
    try:
        staging.write_bytes(raw)
        try:
            os.link(staging, path)
        except FileExistsError:
            if path.read_bytes() != raw:
                raise ValueError(f"Choose a new output path to preserve {path}") from None
    finally:
        staging.unlink(missing_ok=True)
    return hashlib.sha256(raw).hexdigest()


def score_files(lock_file, expected_hash, outcomes, assay_file):
    raw = lock_file.read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if actual != expected_hash.strip().lower():
        raise ValueError("The lock differs from the SHA-256 recorded before unblinding")
    lock = json.loads(raw)
    assay = json.loads(assay_file.read_text(encoding="utf-8-sig"))
    check_acquisition(assay, lock)
    result = score(lock, csv_rows(outcomes, MEASUREMENTS), assay)
    result["lockSha256"] = actual
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    public = sub.add_parser("public", help="Reproduce the registered Cal664sb conditional prediction")
    public.add_argument("--output", type=Path, required=True)
    lock = sub.add_parser("lock", help="Read only IDs and pKa, then record a checksum before unblinding")
    lock.add_argument("--roster", type=Path, required=True)
    lock.add_argument("--study", required=True)
    lock.add_argument("--output", type=Path, required=True)
    test = sub.add_parser("score", help="Open measurements and score the recorded forecasts")
    for flag in ("lock", "measurements", "assay", "output"):
        test.add_argument("--"+flag, type=Path, required=True)
    test.add_argument("--lock-sha256", required=True)
    args = parser.parse_args()
    try:
        if args.action == "public":
            result = public_lock()
        elif args.action == "lock":
            result = make_lock(args.study, csv_rows(args.roster, ROSTER))
        else:
            result = score_files(args.lock, args.lock_sha256, args.measurements, args.assay)
        print(save(args.output, result), args.output)
    except (ValueError, OSError, KeyError, TypeError) as e:
        parser.exit(2, f"{e}\n")


if __name__ == "__main__":
    main()

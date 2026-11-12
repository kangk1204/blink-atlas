"""Evaluate the frozen local pKa model on whole dyes and independent public assays."""

import argparse
import csv
import hashlib
from collections import defaultdict
from itertools import product
import json
from math import isfinite, log10
from pathlib import Path
import random
import re
from statistics import mean
from urllib.parse import urlparse

from analysis.dyes import DATA
from analysis.probe_design import SCOPE, build as holland_data, fit, published


METHODS = ("local", "nearest", "regression")
METHOD_ID = "bracket-logratio-v1"
PROTOCOL = {"cellLine": "COS-7", "target": "H2B-HaloTag", "fixation": "paraformaldehyde",
            "dyeExcitationNm": 647, "onExposureMs": 2, "onCycleMs": 3, "onFrames": 20000,
            "offExposureMs": 25, "offCycleMs": 27, "offFrames": 100000, "offLabeling": "sparse",
            "onGapFrames": 5, "groupingRadiusNm": 100}


def number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not isfinite(value):
        raise ValueError("Expected a finite measured number")
    return value


def identity(value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError("Missing compound identity")
    return re.sub(r"[\s_–-]", "", value).lower()


def predict(train, pka):
    """This function has no access to the test response."""
    number(pka)
    if len(train) < 2:
        return {"status": "insufficient training dyes"}
    grouped = defaultdict(list)
    for row in train:
        if number(row["ratio"]) <= 0:
            raise ValueError("Expected a positive training time ratio")
        grouped[number(row["pKa"])].append(row)
    levels = sorted(grouped)
    support = [levels[0], levels[-1]]
    if len(levels) < 2:
        return {"status": "insufficient distinct pKa values", "range": support}
    if not support[0] <= pka <= support[1]:
        return {"status": "outside training range", "range": support}
    lower = max(x for x in levels if x <= pka)
    upper = min(x for x in levels if x >= pka)
    near = min(levels, key=lambda x: (round(abs(x - pka), 12), x))
    logs = {x: mean(log10(r["ratio"]) for r in grouped[x]) for x in {lower, upper, near}}
    weight = (pka - lower) / (upper - lower) if upper != lower else 0.
    reg = fit(train)
    return {"status": "scored", "range": support, "lower": [r["id"] for r in grouped[lower]],
            "upper": [r["id"] for r in grouped[upper]], "upperWeight": published(weight),
            "nearestDyes": [r["id"] for r in grouped[near]],
            "predictions": {"local": published(10 ** ((1-weight)*logs[lower] + weight*logs[upper])),
                            "nearest": published(10 ** logs[near]),
                            "regression": published(10 ** (reg["intercept"] + reg["slope"]*pka))}}


def uncertainty(values):
    """Paired dye bootstrap for an independent fixed-training test set."""
    n = len(values)
    if n < 2:
        return None
    if n <= 5:
        samples = [mean(s) for s in product(values, repeat=n)]
    else:
        rng = random.Random(8675309)
        samples = [mean(rng.choices(values, k=n)) for _ in range(10000)]
    samples.sort()

    def quantile(p):
        pos = p*(len(samples)-1)
        i = int(pos)
        return published(samples[i] + (pos-i)*(samples[min(i+1, len(samples)-1)] - samples[i]))
    return [quantile(.025), quantile(.975)]


def summarize(rows, independent):
    valid = [r for r in rows if r["status"] == "scored"]
    result = {"total": len(rows), "scored": len(valid), "mae": None, "delta": None, "wins": 0,
              "ties": 0, "losses": 0, "omitOneRange": None, "pairedBootstrap95": None,
              "regressionVsNearest": None}
    if not valid:
        return result
    differences = [r["delta"] for r in valid]
    wins = sum(d > 1e-12 for d in differences)
    losses = sum(d < -1e-12 for d in differences)
    result.update(mae={k: published(mean(r["errors"][k] for r in valid)) for k in METHODS},
                  delta=published(mean(differences)), wins=wins, losses=losses, ties=len(valid)-wins-losses)
    if len(differences) >= 2:
        leave_one = [mean(differences[:i]+differences[i+1:]) for i in range(len(differences))]
        result["omitOneRange"] = [published(min(leave_one)), published(max(leave_one))]
    if independent:
        result["pairedBootstrap95"] = uncertainty(differences)
    regression_delta = [r["errors"]["nearest"]-r["errors"]["regression"] for r in valid]
    result["regressionVsNearest"] = {"delta": published(mean(regression_delta)),
                                     "wins": sum(d > 1e-12 for d in regression_delta),
                                     "pairedBootstrap95": uncertainty(regression_delta) if independent else None}
    return result


def benchmark(key, name, study, mode, targets, reference=None, exclusion=""):
    held_ids = [r["id"] for r in targets]
    if len(set(held_ids)) != len(held_ids):
        raise ValueError(f"Repeated chemical identity in {key}")
    rows = []
    for held in targets:
        train = [r for r in targets if r["id"] != held["id"]] if reference is None else reference
        train_ids = [r["id"] for r in train]
        if held["id"] in train_ids or len(train_ids) != len(set(train_ids)):
            raise ValueError("Training and test identities must be disjoint and unique")
        result = predict(train, held["pKa"])
        if exclusion:
            # Preserve range eligibility separately from an incompatible response definition.
            result = {"status": exclusion, "range": result.get("range"), "rangeStatus": result["status"]}
        row = {"id": held["id"], "dye": held["dye"], "pKa": held["pKa"], "observed": held["ratio"],
               "trainingIds": train_ids, **result}
        if "measurement" in held:
            row["measurement"] = held["measurement"]
        if result["status"] == "scored":
            errors = {m: published(abs(log10(y) - log10(held["ratio"]))) for m,y in result["predictions"].items()}
            row.update(errors=errors, delta=published(errors["nearest"]-errors["local"]))
        rows.append(row)
    return {"id": key, "label": name, "study": study, "mode": mode, "rows": rows,
            "summary": summarize(rows, independent=(reference is not None))}


def observations(path):
    studies = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as file:
        for row in csv.DictReader(file):
            study = row["study"]
            pka = number(float(row["pka"]))
            if row["on_ms"] and row["off_s"] and not row["duty"]:
                on, off = number(float(row["on_ms"])), number(float(row["off_s"]))
                if min(on, off) <= 0:
                    raise ValueError("Invalid ON or OFF time")
                ratio = on / (1000*off)
            elif row["duty"] and not (row["on_ms"] or row["off_s"]):
                duty = number(float(row["duty"]))
                if not 0 < duty < 1:
                    raise ValueError("Duty fraction must be between zero and one")
                ratio = duty/(1-duty)
            else:
                raise ValueError("Ambiguous or missing response")
            entry = {"id": row["dye_id"], "dye": row["dye"], "pKa": pka, "ratio": published(ratio)}
            if row["bound_pka"]:
                entry["boundPka"] = number(float(row["bound_pka"]))
            studies[study, row["condition"]].append(entry)
    return studies


def compare_titrations(group, free_set):
    """Exploratory matched-fold ablation, holding the interpolation rule and eligibility fixed."""
    bound = [{**r, "pKa": r["boundPka"]} for r in group]
    bound_set = benchmark("bound", "", "zheng", "transfer-local", bound)
    rows = []
    for free, attached in zip(free_set["rows"], bound_set["rows"], strict=True):
        paired = free["status"] == attached["status"] == "scored"
        rows.append({"id": free["id"], "dye": free["dye"], "paired": paired,
                     "free": free, "bound": attached})
    paired = [r for r in rows if r["paired"]]
    maes = {t: {method: published(mean(r[t]["errors"][method] for r in paired)) for method in ("local", "nearest")}
            for t in ("free", "bound")} if paired else None
    return {"power": int(free_set["condition"]), "scored": len(paired), "total": len(rows), "mae": maes, "rows": rows}


def independent_panel(path, catalog):
    """Validate an explicitly curated, acquisition-matched panel before accessing responses."""
    panel = json.loads(path.read_text(encoding="utf-8"))
    protocol = panel.get("protocol")
    if (panel.get("schema") != "blink-atlas/matched-h2b/1" or panel.get("scope") != SCOPE
            or not isinstance(protocol, dict) or any(protocol.get(k) != v for k, v in PROTOCOL.items())
            or panel.get("independentAcquisition") is not True
            or panel.get("identitiesReviewedAgainstHolland") is not True):
        raise ValueError("Unmatched protocol or missing independent acquisition/identity review")
    source = panel.get("source")
    if not isinstance(source, str) or urlparse(source).scheme != "https" or not urlparse(source).hostname:
        raise ValueError("Supply the public source for the independent acquisition")
    if any(s in source.lower() for s in ("s41592-026-03062-5", "2024.02.23.581625", "pmc13167466", "pmc11101128")):
        raise ValueError("Development publication cannot serve as an independent study")
    known = {identity(s) for d in catalog["dyes"] for s in (d["id"], d["ligand"])}
    seen, result = set(), []
    if not isinstance(panel.get("dyes"), list) or not panel["dyes"]:
        raise ValueError("Supply measurements for at least one independently acquired dye")
    for dye in panel["dyes"]:
        key = dye["parentInchiKey"]
        if not isinstance(key, str) or not re.fullmatch("[A-Z]{14}-[A-Z]{10}-[A-Z]", key):
            raise ValueError("Supply a standardized free-parent InChIKey")
        aliases = dye.get("aliases", [])
        if not isinstance(aliases, list):
            raise ValueError("Aliases must be a list")
        names = {identity(s) for s in [dye["name"], key, *aliases]}
        if names & (known | seen):
            raise ValueError("Known or repeated dye identity in external panel")
        seen.update(names)
        if dye.get("monotonicFreeTitration") is not True or not isinstance(dye.get("sourceLocator"), str) or not dye["sourceLocator"].strip():
            raise ValueError("An exact monotonic free-parent titration and source location are required")
        pka = number(dye["pKa"])
        measurement = {"sourceLocator": dye["sourceLocator"].strip(), "aliases": aliases}
        if ("reportedMeanRatio" in dye) == ("preparationRatios" in dye):
            raise ValueError("Supply either the published mean ratio or the individual preparation ratios")
        if "reportedMeanRatio" in dye:
            ratio = number(dye["reportedMeanRatio"])
            measurement["reportedMeanRatio"] = ratio
        else:
            if not isinstance(dye["preparationRatios"], list) or not dye["preparationRatios"]:
                raise ValueError("Supply a nonempty list of independent preparation ratios")
            ratios = [number(r) for r in dye["preparationRatios"]]
            if min(ratios) <= 0:
                raise ValueError("Preparation ratios must be positive")
            measurement["preparationRatios"] = ratios
            ratio = mean(ratios)
        if ratio <= 0:
            raise ValueError("The measured mean ratio must be positive")
        result.append({"id": key, "dye": dye["name"], "pKa": pka, "ratio": ratio, "measurement": measurement})
    return result, source, protocol


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(directory=DATA, matched=None):
    catalog = holland_data(directory)
    dev = [{"id": p["id"], "dye": p["id"], "pKa": p["pKa"], "ratio": p["ratio"]}
           for p in catalog["calibration"]["pairs"]]
    metadata_file, csv_file = directory/"analogue-studies.json", directory/"analogue-observations.csv"
    meta = json.loads(metadata_file.read_text(encoding="utf-8"))
    public = observations(csv_file)
    target, source, protocol = independent_panel(matched, catalog) if matched else ([], None, None)
    primary = benchmark("matched", "Independent fixed H2B · frozen model", "holland", "primary", target, dev)
    sets = [primary, benchmark("development", "Holland · leave one dye out", "holland", "development", dev)]
    titrations = []
    for (study, condition), group in public.items():
        if study not in meta["studies"]:
            raise ValueError(f"Unknown study {study}")
        name = {"remmel": f"Remmel · {condition}", "tyson": "Tyson · glass",
                "zheng": f"Zheng · HaloTag at {condition} W/cm²"}[study]
        key = study+"-"+condition.lower()
        exclusion = "different response window" if study == "zheng" else ""
        for suffix, label, mode, ref, reason in (
            ("frozen", "frozen Holland model", "transfer-frozen", dev, exclusion),
            ("local", "study-specific calibration", "transfer-local", None, ""),
        ):
            item = benchmark(key+"-"+suffix, name+" · "+label, study, mode, group, ref, reason)
            item["condition"] = condition
            sets.append(item)
            if study == "zheng" and mode == "transfer-local":
                titrations.append(compare_titrations(group, item))
    identities = {r["id"] for group in public.values() for r in group}
    summary = primary["summary"]
    interval = summary["pairedBootstrap95"]
    gate = "independent review warranted" if interval and interval[0] > 0 else "advantage unconfirmed"
    return {"schema": "blink-atlas/analogue-benchmark/1", "method": METHOD_ID,
            "frozenRevision": meta["frozenRevision"], "sourceHashes": {"observations": digest(csv_file),
                "studies": digest(metadata_file), "hollandWorkbook": catalog["calibration"]["sourceSha256"]},
            "uniqueExternalDyes": len(identities), "externalConditions": sum(len(g) for g in public.values()),
            "matchedSource": source, "matchedProtocol": protocol, "matchedPanelSha256": digest(matched) if matched else None,
            "primaryVerdict": gate, "sets": sets,
            "exploratoryTitrations": {"study": "zheng", "status": "post hoc, within-assay",
                                      "endpoint": "whole-recording occupancy odds", "byPower": titrations}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DATA)
    parser.add_argument("--matched-panel", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.matched_panel and not args.output:
        parser.error("Name an --output for your independent panel")
    result = build(args.data_dir, args.matched_panel)
    path = args.output or args.data_dir/"analogue-benchmark.json"
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8", newline="\n")
    for item in result["sets"]:
        s = item["summary"]
        print(item["id"], f"{s['scored']}/{s['total']}", s["mae"], s["delta"], s["wins"], s["losses"])


if __name__ == "__main__":
    main()

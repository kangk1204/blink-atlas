"""Rebuild the Holland ratio audit and exhaustive, whole-dye sensitivity results."""

import argparse
from collections import Counter
import csv
import gzip
import hashlib
from io import TextIOWrapper
from itertools import combinations, product
import json
from math import sqrt
from pathlib import Path

import numpy as np

from analysis.analogue_benchmark import METHODS, benchmark, predict
from analysis.dyes import DATA, HASH, URL
from analysis.probe_design import published
from analysis.ratio_audit import PAIRINGS, read_audit, write_audit


LOCK_SHA = "c4e9ae62123abff64ca703a72b8a3c63e29e2ccec312a85f02f969d19eaaca8f"
LOCK_FILE = "h2b-cal664sb-lock.json"
PAIRS = list(combinations(METHODS, 2))


def roster(rows, column):
    return [dict(id=r["dye"], dye=r["dye"], pKa=r["free_pka"], ratio=r[column]) for r in rows]


def boundaries(predictions):
    ordered = sorted(predictions, key=predictions.get)
    return [{"models": [a, b], "ratio": published(sqrt(predictions[a] * predictions[b]))}
            for a, b in zip(ordered, ordered[1:])]


def single(rows, key, pka):
    data = roster(rows, key)
    folds = benchmark(key, key, "holland", "development", data)
    forecasts = predict(data, pka)["predictions"]
    return dict(folds=folds["rows"], metrics=folds["summary"], cal=forecasts,
                calBoundaries=boundaries(forecasts))


def log_weights(reference, query, held=None):
    """All three forecasts are weighted sums on the log scale for fixed pKas."""
    train = [r for i, r in enumerate(reference) if i != held]
    rule = predict(train, query)
    if rule["status"] != "scored":
        return None
    index = {r["id"]: i for i, r in enumerate(reference)}
    w = np.zeros((len(METHODS), len(reference)))
    xs = np.array([r["pKa"] for r in train])
    dx = xs - xs.mean()
    for r, coefficient in zip(train, 1 / len(xs) + (query - xs.mean()) * dx / (dx @ dx)):
        w[2, index[r["id"]]] = coefficient
    lo, hi = (index[rule[key][0]] for key in ("lower", "upper"))
    x0, x1 = (reference[i]["pKa"] for i in (lo, hi))
    t = (query - x0) / (x1 - x0) if x0 != x1 else 0
    w[0, lo] += 1 - t
    w[0, hi] += t
    w[1, index[rule["nearestDyes"][0]]] = 1
    return w


def span(values):
    return [published(float(values.min())), published(float(values.max()))]


def rankings(errors):
    """Group actual orders, allowing ties at the benchmark's 1e-12 tolerance."""
    counts = Counter()
    for values in errors:
        groups = []
        for i in np.argsort(values, kind="stable"):
            if groups and values[i] - values[groups[-1][0]] <= 1e-12:
                groups[-1].append(i)
            else:
                groups.append([i])
        counts[" < ".join(" = ".join(METHODS[i] for i in g) for g in groups)] += 1
    return dict(sorted(counts.items()))


def error_summary(values):
    pairs = {}
    for a, b in PAIRS:
        gain = values[:, METHODS.index(b)] - values[:, METHODS.index(a)]
        pairs[f"{a}_vs_{b}"] = {"secondMinusFirst": span(gain), "firstWins": int((gain > 1e-12).sum()),
                               "ties": int((abs(gain) <= 1e-12).sum()), "secondWins": int((gain < -1e-12).sum())}
    return {"errors": {m: span(values[:, i]) for i, m in enumerate(METHODS)}, "comparisons": pairs,
            "rankCountsLowestFirst": rankings(values)}


def exhaustive(rows, pka, archive):
    ref = roster(rows, "j_reported")
    grid = np.array(list(product(range(len(PAIRINGS)), repeat=len(rows))), dtype=np.uint8)
    logs = np.log10(np.array([[r["paired_" + p] for p in PAIRINGS] for r in rows]))
    train_labels = logs[np.arange(len(rows)), grid]
    errors, fold_reports = [], []
    for i, datum in enumerate(ref):
        weights = log_weights(ref, datum["pKa"], i)
        if weights is None:
            fold_reports.append({"dye": datum["id"], "status": "outside training range"})
            continue
        predictions = np.einsum("si,mi->sm", train_labels, weights)
        fold_error = abs(predictions - train_labels[:, i, None])
        errors.append(fold_error)
        fold_reports.append({"dye": datum["id"], "status": "scored", **error_summary(fold_error)})
    mae = np.mean(errors, axis=0)
    cal = np.power(10., np.einsum("si,mi->sm", train_labels, log_weights(ref, pka)))
    # Boundaries are for adjacent predictions, so check the order in every scenario.
    cal_orders, order_counts = np.unique(np.argsort(cal, axis=1), axis=0, return_counts=True)
    partition = []
    boundary_columns = {}
    for order, count in zip(cal_orders, order_counts):
        mask = (np.argsort(cal, axis=1) == order).all(axis=1)
        entries = []
        for left, right in zip(order, order[1:]):
            values = np.sqrt(cal[:, left] * cal[:, right])
            a, b = METHODS[left], METHODS[right]
            entries.append({"models": [a, b], "range": span(values[mask])})
            boundary_columns[f"tie_{a}_{b}"] = values
        partition.append({"ascendingPredictions": [METHODS[i] for i in order], "scenarios": int(count),
                          "adjacentBoundaries": entries})
    columns = [*(r["dye"] + "_off" for r in rows), *("mae_" + m for m in METHODS),
               *("cal_" + m for m in METHODS), *boundary_columns]
    values = np.column_stack([mae, cal, *boundary_columns.values()])
    with archive.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        with TextIOWrapper(gz, encoding="ascii", newline="") as out:
            writer = csv.writer(out, lineterminator="\n")
            writer.writerow(columns)
            for choices, numbers in zip(grid, values):
                writer.writerow([*(PAIRINGS[i] for i in choices), *(format(v, ".10g") for v in numbers)])
    shared = np.all(grid == grid[:, :1], axis=1)
    return {"scenarioCount": len(grid), "pairingLabels": list(PAIRINGS), "dyeOrder": [r["dye"] for r in rows],
            "scoredDyes": len(errors), "folds": fold_reports, "overall": error_summary(mae),
            "samePermutationForAllDyes": error_summary(mae[shared]),
            "calForecastRanges": {m: span(cal[:, i]) for i, m in enumerate(METHODS)}, "calOrder": partition}


def build(directory, output):
    rows = read_audit(directory)
    numeric = [r for r in rows if r["free_pka"] is not None]
    lock_bytes = (directory / LOCK_FILE).read_bytes()
    if hashlib.sha256(lock_bytes).hexdigest() != LOCK_SHA:
        raise ValueError("The original Cal664sb registration changed")
    lock = json.loads(lock_bytes)
    target = lock["dyes"][0]
    source = "https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/"
    output.mkdir(parents=True, exist_ok=True)
    write_audit(rows, output / "ratio-audit.csv")
    methods = {k: single(numeric, k, target["pKa"]) for k in ("j_reported", "ratio_of_means", "paired_BCD")}
    if methods["j_reported"]["cal"] != target["predictions"]:
        raise ValueError("Published ratios must exactly recover the registered predictions")
    scan = exhaustive(numeric, target["pKa"], output / "ratio-scenarios.csv.gz")
    audit = {"schema": "blink-atlas/ratio-sensitivity/1",
             "input": {"workbook": URL, "sha256": HASH, "caption": source + "PMC13167466.1.xml",
                       "figure": source + "41592_2026_3062_Fig6_ESM.jpg"},
             "audit": [{k: published(v) if isinstance(v, float) else v for k, v in r.items()} for r in rows],
             "registered": {"file": LOCK_FILE, "sha256": LOCK_SHA, "dye": target["blindId"],
                            "pKa": target["pKa"], "forecasts": target["predictions"],
                            "boundaries": lock["equalErrorBoundaries"]},
             "estimands": methods, "pairings": scan}
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DATA)
    parser.add_argument("--output-dir", type=Path, default=DATA)
    args = parser.parse_args()
    result = build(args.data_dir, args.output_dir)
    (args.output_dir / "ratio-sensitivity.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    print(f"{len(result['audit'])} dyes audited. {result['pairings']['scenarioCount']:,} scenarios enumerated.")


if __name__ == "__main__":
    main()

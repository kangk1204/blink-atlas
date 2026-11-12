"""Recompute the published H2B precision bounds and the separate dye-model yields.

Run from the Blink Atlas repository root:
    python -m validation.h2b_precision.reproduce
Use --precision-workbook FILE to run offline with the publisher's Fig. 5 source.
"""

import argparse
import hashlib
import io
import json
from pathlib import Path
from urllib.request import urlopen

from openpyxl import load_workbook

from analysis.dyes import DATA, FILENAME, extract
from analysis.model import evaluate_prior


SOURCE = "https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_MOESM11_ESM.xlsx"
SHA256 = "a82a16a051b1c9fc66c6dc4f67acc47c870d8d7b5db2d1132cba1f601d5cd11c"
SHEET = "Extended Data Figure 5b"
THRESHOLD_NM = 20.0
DYES = ("JF639b", "JF630b")


def h2b_precision(workbook_bytes):
    if hashlib.sha256(workbook_bytes).hexdigest() != SHA256:
        raise ValueError("Extended Data Figure 5 workbook differs from the audited version")
    book = load_workbook(io.BytesIO(workbook_bytes), read_only=True, data_only=True)
    try:
        sheet = book[SHEET]
        heading, *rows = list(sheet.iter_rows(min_col=1, max_col=2, values_only=True))
    finally:
        book.close()
    if tuple(" ".join(cell.split()) for cell in heading) != ("Localization precision (nm)", "histone H2B"):
        raise ValueError("Unexpected precision axis or specimen")
    if len(rows) % 2 or not rows:
        raise ValueError("The published staircase has an incomplete bin")

    total = certain = ambiguous = 0
    crossing = None
    previous_edge = None
    for (left, count), (right, repeated_count) in zip(rows[::2], rows[1::2]):
        if left >= right or (previous_edge is not None and left != previous_edge):
            raise ValueError("Noncontiguous histogram edges")
        if count != repeated_count or not isinstance(count, int) or count < 0:
            raise ValueError("Unexpected histogram weights")
        previous_edge = right
        total += count
        if right <= THRESHOLD_NM:
            certain += count
        elif left <= THRESHOLD_NM:
            crossing = [left, right]
            ambiguous += count

    return {
        "dye": "JF630b", "target": "fixed-cell H2B–HaloTag",
        "sheet": SHEET, "source": SOURCE, "source_sha256": SHA256, "bins": len(rows) // 2,
        "threshold_nm": THRESHOLD_NM, "histogram_total": total,
        "threshold_bin_nm": crossing, "threshold_bin_count": ambiguous,
        "qualifying_count_bounds": [certain, certain + ambiguous],
        "qualifying_fraction_bounds": [certain / total, (certain + ambiguous) / total],
        "scope": "Displayed scalar-precision histogram only. No matched low-density, per-field denominator.",
    }


def predictions_from_measured_kinetics():
    dataset = extract(DATA / FILENAME)
    columns = ("onMs", "offMs", "photonsPerMs", "blinksPerMolecule")
    measured = {
        dye["id"]: {key: dye["metrics"][key]["mean"] for key in columns}
        for dye in dataset["dyes"] if dye["id"] in DYES
    }
    predictions = []
    for density in (100, 1000):
        settings = {"densityMin": density, "densityMax": density, "candidatePanel": "validation"}
        result = evaluate_prior(dataset["dyes"], settings)
        rates = {row["id"]: row["rate"] for row in result["ranking"] if row["id"] in DYES}
        predictions.append({
            "assumed_labels_per_um2": density, "isolated_detections_per_um2_s": rates,
            "JF639b_over_JF630b": rates["JF639b"] / rates["JF630b"],
        })
    return {
        "measured_means": measured, "source": dataset["source"]["url"],
        "model": "Publication-mean, stationary first-burst model with a 250 nm exclusion radius",
        "exposure_ms": 25, "frame_cycle_ms": 27, "min_expected_photons": 100,
        "predictions": predictions,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision-workbook", type=Path)
    args = parser.parse_args()
    if args.precision_workbook:
        data = args.precision_workbook.read_bytes()
    else:
        with urlopen(SOURCE, timeout=30) as response:
            data = response.read()
    result = {
        "observed_JF630b": h2b_precision(data),
        "photophysics": predictions_from_measured_kinetics(),
        "observed_matched_low_density_JF639b_over_JF630b_at_20nm": None,
    }
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

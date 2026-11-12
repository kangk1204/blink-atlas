"""Rebuild the measured dye records from the publisher's unmodified workbook.

Run from the project root with `python -m analysis.dyes`. The small source
workbook is vendored with the website so the analysis also works offline.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import mean, stdev
from urllib.request import urlopen

from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"
FILENAME = "41592_2026_3062_MOESM10_ESM.xlsx"
URL = "https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/" + FILENAME
HASH = "1e26a58d558f7d96518aca0264b0efed6e5574b67ba5f18899965802fedf5031"
DYE_NAMES = ["HM-SiR", "JFX650b", "JF646b", "JF639b", "JF635b",
             "JF630b", "JF629b", "JF626b", "JF614b"]


def read_values(sheet, row, first, last):
    values = [sheet.cell(row=row, column=col).value for col in range(first, last + 1)]
    if not all(isinstance(v, (int, float)) and math.isfinite(v) and v >= 0 for v in values):
        raise ValueError(f"Missing or invalid value in {sheet.title}, row {row}")
    return values


def replicate_metric(book, panel, row, unit, has_median=False):
    name = "Extended Data Figure 4" + panel
    sheet = book[name]
    a = read_values(sheet, row, 2, 4)
    metric = {
        "mean": mean(a), "sd": stdev(a), "n": len(a), "replicates": a,
        "unit": unit, "sheet": name, "cells": f"B{row}:D{row}",
    }
    if has_median:
        medians = read_values(sheet, row, 5, 7)
        metric.update({
            "median": mean(medians), "medianSd": stdev(medians),
            "medianReplicates": medians, "medianCells": f"E{row}:G{row}",
        })
    return metric


def reported_ratio(book, row):
    name = "Extended Data Figure 4h-i"
    avg, sd, n, med, med_sd, med_n = read_values(book[name], row, 2, 7)
    return {
        "mean": avg, "sd": sd, "n": n, "median": med, "medianSd": med_sd, "medianN": med_n,
        "sheet": name, "cells": f"B{row}:D{row}", "medianCells": f"E{row}:G{row}",
        "unit": "ratio", "definition": "Reported mean-based and median-based on:off ratios.",
    }


def extract(path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != HASH:
        raise ValueError(f"Source checksum mismatch for {path.name}: {digest}")
    book = load_workbook(path, read_only=True, data_only=True)
    dyes = []
    try:
        for i, name in enumerate(DYE_NAMES, 1):
            row = i + 1
            ligand = f"{i}HTL"
            if book["Extended Data Figure 4b-c"].cell(row, 1).value != ligand:
                raise ValueError(f"Unexpected ligand at source row {row}")
            metrics = {
                "onMs": replicate_metric(book, "b-c", row, "ms", has_median=True),
                "offMs": replicate_metric(book, "e-f", row, "ms", has_median=True),
                "photonsPerMs": replicate_metric(book, "l", row, "photons/ms"),
                "blinksPerMolecule": replicate_metric(book, "k", row, "blinks/molecule"),
                "onOffRatio": reported_ratio(book, row + 1),
            }
            dyes.append({"id": name, "name": name, "ligand": ligand,
                         "conjugate": "HaloTag ligand", "metrics": metrics})
    finally:
        book.close()
    return {
        "schema": "smlm-dyes/1", "dyes": dyes,
        "source": {
            "title": "A series of spontaneously blinking dyes for super-resolution microscopy",
            "authors": "Holland, Plutkis, English et al.", "doi": "10.1038/s41592-026-03062-5",
            "article": "https://doi.org/10.1038/s41592-026-03062-5", "url": URL,
            "local": "data/" + FILENAME, "sha256": digest, "license": "CC BY 4.0",
            "licenseUrl": "https://creativecommons.org/licenses/by/4.0/",
            "transformation": "Extended Data Fig. 4 extracted to JSON. Summary means and sample SDs calculated across the three reported experiments. Per-experiment medians and their arithmetic mean use the separate median columns. Ratio summaries are transcribed from panel h-i.",
            "nameMapping": "Full-text results and Extended Data Fig. 4 map ligands 1HTL through 9HTL to the displayed names.",
        },
        "acquisition": {
            "sample": "Fixed COS-7 cells expressing histone H2B–HaloTag",
            "excitationNm": 647, "onExposureMs": 2, "onCycleMsReportedMethods": 3,
            "onDistributionFirstBinMs": 3.4, "offExposureMs": 25, "offCycleMs": 27,
            "independentExperiments": 3,
            "notes": "Kinetics were measured at subsaturating labeling. Off times were measured at about 100-fold lower labeling concentration than on times, over 45 min. The authors describe off times as upper bounds under these conditions. Methods round the on-time acquisition cycle to 3 ms; panel d starts at 3.4 ms.",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--workbook", type=Path, default=DATA / FILENAME)
    parser.add_argument("--output", type=Path, default=DATA / "dyes.json")
    args = parser.parse_args()
    if not args.workbook.exists():
        with urlopen(URL, timeout=45) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != HASH:
            raise ValueError("Downloaded workbook does not match the source checksum")
        args.workbook.parent.mkdir(parents=True, exist_ok=True)
        args.workbook.write_bytes(data)
    result = extract(args.workbook)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {len(result['dyes'])} measured dye records to {args.output}")


if __name__ == "__main__":
    main()

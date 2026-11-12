"""Rebuild the probe-design workbench from the Holland et al. source workbooks."""

import argparse
from contextlib import ExitStack, closing
import hashlib
import json
from math import isclose, isfinite, log10
from pathlib import Path
from statistics import linear_regression, mean, median

from openpyxl import load_workbook
from openpyxl.utils import get_column_letter

from analysis.dyes import DATA, DYE_NAMES, HASH, reported_ratio


ARTICLE = "https://doi.org/10.1038/s41592-026-03062-5"
BASE = "https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/"
SI_FILE = "41592_2026_3062_MOESM1_ESM.pdf"
SI_HASH = "adf5bc338599bf8651ea482f9a31d09545486fc692ab9b713012c2afcf223974"
HASHES = {
    7: "2d8e5b9c52d1e38fbc57691695076d15dcec0920de3faef9bb1a09d67a235991",
    8: "9f0b5764b0a86978aa8cc02bbf816949be0357ad33234f0f93e239d0d2d2dd4e",
    9: "c7f90e03b582edcf2618118c7be0e2c9af6481d61960901637b09db4147b5e1a",
    10: HASH,
    13: "86acf64ef89ca05f32ad9b385ea0c2859e50bb57bdc5f73ce1aeace258bebf1b",
}
SHEETS = {7: "Figure 1a", 8: "Figure 2b-d", 9: "Extended Data Figure 2j-k",
          10: "Extended Data Figure 4j", 13: "Extended Data Figure 8e-f"}

# Figure 1a in the article and Supplementary Note 2, compounds 1–9.
# The drawing in the PDF is the free parent in its closed form.
# Last two columns are the physical PDF page and a crop on a 521 × 674 page preview.
CHEMISTRY = [
    ("Dimethylamino", 656, 671, "5.15", 41, (213, 248, 307, 314)),
    ("Pyrrolidine-d₈", 660, 674, "4.48", 41, (197, 443, 323, 522)),
    ("Azetidine", 657, 671, "ND", 42, (211, 102, 310, 165)),
    ("3-methoxyazetidine", 653, 666, "3.97", 42, (191, 293, 330, 362)),
    ("3-fluoroazetidine", 647, 661, "2.95", 42, (200, 506, 320, 575)),
    ("3-(methylsulfonyl)azetidine", 646, 660, "2.58", 43, (194, 186, 326, 266)),
    ("3-cyanoazetidine", 645, 660, "2.25", 43, (196, 408, 325, 478)),
    ("3-(trifluoromethyl)azetidine", 646, 660, "1.22", 44, (194, 70, 327, 139)),
    ("3,3-difluoroazetidine", 636, 651, "<1", 44, (199, 297, 321, 369)),
]
SCOPE = "hydroxymethyl-sir-monotonic-fixed-cos7-h2b"
TARGET = "Mean ON:OFF time ratio of the HaloTag ligand in fixed COS-7 H2B–HaloTag."


def file_name(number):
    return f"41592_2026_3062_MOESM{number}_ESM.xlsx"


def ref(book, sheet, cells):
    return {"workbook": str(book), "sheet": sheet, "cells": cells}


def values(sheet, row, col, count):
    result = [sheet.cell(row, i).value for i in range(col, col + count)]
    if not all(isinstance(v, (float, int)) and isfinite(v) for v in result):
        raise ValueError(f"Invalid measurement in {sheet.title}, row {row}")
    return result


def cell_range(row, col, width=2):
    return f"{get_column_letter(col)}{row}:{get_column_letter(col + width - 1)}{row}"


def titration(sheet, index):
    col = 2 * index
    if sheet.cell(1, col).value != index:
        raise ValueError("Titration dye order changed")
    points = []
    for row in range(2, sheet.max_row + 1):
        if all(sheet.cell(row, i).value is None for i in (col, col + 1)):
            continue
        readings = values(sheet, row, col, 2)
        points.append({"pH": values(sheet, row, 1, 1)[0], "readings": readings,
                       "mean": mean(readings), "cells": cell_range(row, col), "pHCell": f"A{row}"})
    return points


def parents(books):
    qy = books[7][SHEETS[7]]
    results = []
    for i, (name, spec) in enumerate(zip(DYE_NAMES, CHEMISTRY, strict=True), 1):
        arm, absorption, emission, pka, page, _ = spec
        if qy.cell(i + 1, 1).value != i:
            raise ValueError("Quantum yield dye order changed")
        readings = values(qy, i + 1, 2, 2)
        ratio = reported_ratio(books[10], i + 2)
        results.append({
            "id": name, "compound": i, "ligand": f"{i}HTL", "auxochrome": arm,
            "absorptionNm": absorption, "emissionNm": emission, "pKaFigure1": pka,
            "structure": {"image": f"media/probe-structures/{i}.png", "pdfPage": page},
            "quantumYield": {"readings": readings, "mean": mean(readings),
                             "source": ref(7, qy.title, cell_range(i + 1, 2))},
            "titration": titration(books[9][SHEETS[9]], i),
            "ligandRatio": {"mean": ratio["mean"], "sd": ratio["sd"], "n": ratio["n"],
                            "source": ref(10, ratio["sheet"], ratio["cells"])},
        })
    return results


def sofi(books):
    er, lyso = books[8][SHEETS[8]], books[13][SHEETS[13]]
    contexts = {
        "er": {"name": "ER membrane · Sec61β", "nCells": 4, "metrics": ["useful", "snr", "tau"],
               "protocol": "U2OS · 100 nM for 30 min per dye · 10,000 frames at 100 Hz · 637 nm HILO",
               "fig": "Figure 2b–d", "dyes": {}},
        "lysosome": {"name": "Lysosome", "nCells": 5, "metrics": ["snr", "tau"],
                     "protocol": "U2OS · 90 min labeling · 20 nM JF639b and 100 nM JF614b · 100 Hz",
                     "fig": "Extended Data Figure 8e–f", "dyes": {}},
    }
    for i, name in enumerate(DYE_NAMES, 1):
        states = {}
        for state, row in (("live", 3 + i), ("fixed", 17 + i)):
            if er.cell(row, 1).value != f"{i}HTL":
                raise ValueError("ER ligand order changed")
            states[state] = {k: read_sem(er, 8, row, col) for k, col in (("snr", 2), ("tau", 4), ("useful", 6))}
        contexts["er"]["dyes"][name] = states
    for col in (2, 4):
        name = lyso.cell(1, col).value
        states = {}
        for state, row in (("live", 3), ("fixed", 4)):
            states[state] = {k: read_sem(lyso, 13, r, col) for k, r in (("snr", row), ("tau", row + 3))}
        contexts["lysosome"]["dyes"][name] = states
    return contexts


def read_sem(sheet, book, row, col):
    avg, sem = values(sheet, row, col, 2)
    return {"mean": avg, "sem": sem, "source": ref(book, sheet.title, cell_range(row, col))}


def fit(rows):
    xs = [r["pKa"] for r in rows]
    ys = [log10(r["ratio"]) for r in rows]
    slope, intercept = linear_regression(xs, ys)
    return {"slope": slope, "intercept": intercept, "medianLog": median(ys), "pKaRange": [min(xs), max(xs)]}


def published(value):
    # Keep computed results byte-reproducible across CPython/libm versions.
    return float(format(value, ".12g"))


def calibrate(sheet, dyes):
    by_ligand = {d["ligand"]: d for d in dyes}
    rows = []
    for row in range(3, 9):
        dye = by_ligand[sheet.cell(row, 1).value]
        pka, ratio, sd, n = values(sheet, row, 2, 4)
        if str(pka) != dye["pKaFigure1"] or not isclose(ratio, dye["ligandRatio"]["mean"]) or n != 3:
            raise ValueError("Calibration disagrees with Fig. 1a or Extended Data Fig. 4h–i")
        rows.append({"id": dye["id"], "pKa": pka, "ratio": ratio, "sd": sd, "n": n,
                     "source": ref(10, sheet.title, f"A{row}:E{row}")})
    folds = []
    for held in rows:
        train = [r for r in rows if r["id"] != held["id"]]
        f = fit(train)
        fold = {"held": held["id"], "training": [r["id"] for r in train], "trainingRange": f["pKaRange"],
                "observed": held["ratio"], "pKa": held["pKa"]}
        if f["pKaRange"][0] <= held["pKa"] <= f["pKaRange"][1]:
            pred = f["intercept"] + f["slope"] * held["pKa"]
            obs = log10(held["ratio"])
            nearest = min(train, key=lambda r: (round(abs(r["pKa"] - held["pKa"]), 12), r["pKa"]))
            fold.update({"status": "scored", "predicted": published(10 ** pred), "baseline": published(10 ** f["medianLog"]),
                         "nearestDye": nearest["id"], "nearestRatio": nearest["ratio"],
                         "error": published(abs(pred - obs)), "baselineError": published(abs(f["medianLog"] - obs)),
                         "nearestError": published(abs(log10(nearest["ratio"]) - obs))})
        else:
            fold["status"] = "outside training range"
        folds.append(fold)
    scored = [f for f in folds if f["status"] == "scored"]
    full_fit = fit(rows)
    for key in ("slope", "intercept", "medianLog"):
        full_fit[key] = published(full_fit[key])
    return {"id": "holland-h2b-pka-log-ratio-v1", "scope": SCOPE, "target": TARGET,
            "sourceSha256": HASH, "pairs": rows, "fit": full_fit, "folds": folds,
            "validation": {"scored": len(scored), "total": len(rows), "mae": published(mean(f["error"] for f in scored)),
                           "baselineMae": published(mean(f["baselineError"] for f in scored)),
                           "nearestMae": published(mean(f["nearestError"] for f in scored)),
                           "winsVsMedian": sum(f["error"] < f["baselineError"] for f in scored),
                           "winsVsNearest": sum(f["error"] < f["nearestError"] for f in scored)}}


def build(directory=DATA):
    sources = {}
    with ExitStack() as stack:
        books = {}
        for number, digest in HASHES.items():
            filename = file_name(number)
            path = directory / filename
            if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                raise ValueError(f"Checksum mismatch for {filename}")
            books[number] = stack.enter_context(closing(load_workbook(path, data_only=True)))
            sources[str(number)] = {"file": f"data/{filename}", "url": BASE + filename, "sha256": digest,
                                    "sheet": SHEETS[number]}
        dyes = parents(books)
        return {"schema": "blink-atlas/probe-design/1", "article": ARTICLE, "license": "CC BY 4.0",
                "si": {"url": BASE + SI_FILE, "sha256": SI_HASH}, "sources": sources,
                "dyes": dyes, "contexts": sofi(books), "calibration": calibrate(books[10][SHEETS[10]], dyes)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DATA)
    parser.add_argument("--output", type=Path, default=DATA / "probe-design.json")
    args = parser.parse_args()
    data = build(args.data_dir)
    args.output.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                           encoding="utf-8", newline="\n")
    print("Built probe-design.json", data["calibration"]["validation"])


if __name__ == "__main__":
    main()

"""Audit the Holland dwell-time estimand directly from Extended Data Fig. 4."""

import argparse
import csv
from itertools import permutations
from math import sqrt
from pathlib import Path
from statistics import mean

from openpyxl import load_workbook

from analysis.dyes import DATA, FILENAME, extract
from analysis.probe_design import published


PAIRINGS = tuple("".join(p) for p in permutations("BCD"))
PANEL = "Extended Data Figure 4"
# Hand transcription, Extended Data Fig. 4n, "On:off ratio (%)" row from left to right.
# Keep the displayed significant figures. Its underlying mean is also supplied in 4h–i.
PANEL_N_PCT = ("0.06", "0.11", "0.07", "0.10", "0.048", "0.014", "0.0099", "0.029", "0.0019")


def paired_mean(on, off, pairing):
    """Keep ON columns B,C,D fixed and rearrange the three OFF columns."""
    return mean(value / off["BCD".index(col)] for value, col in zip(on, pairing, strict=True))


def deltas(value, reported):
    return abs(value - reported), 100 * (value / reported - 1)


def read_audit(directory=DATA):
    path = directory / FILENAME
    dyes = extract(path)["dyes"]  # Checks the publisher workbook SHA-256.
    book = load_workbook(path, data_only=True, read_only=True)
    try:
        panel = list(book[PANEL + "j"].values)
        if panel[1][1:5] != ("pKa", "average", "SD", "n") or book[PANEL + "h-i"]["C2"].value != "SD":
            raise ValueError("Unexpected uncertainty headers in Holland source")
        fit_rows = {r[0]: (i, *r[1:5]) for i, r in enumerate(panel[2:], 3) if r[0]}
        rows = []
        for i, dye in enumerate(dyes, 2):
            metrics, ligand = dye["metrics"], dye["ligand"]
            hrow = i + 1
            if any(book[PANEL + suffix].cell(row, 1).value != ligand
                   for suffix, row in (("e-f", i), ("h-i", hrow))):
                raise ValueError(f"Source ligand mismatch for {ligand}")
            on, off = (metrics[k]["replicates"] for k in ("onMs", "offMs"))
            rep = metrics["onOffRatio"]
            rom = mean(on) / mean(off)
            paired = [paired_mean(on, off, p) for p in PAIRINGS]
            row = {"dye": dye["id"], "ligand": ligand,
                   "figure_4n_reported_pct": PANEL_N_PCT[i-2],
                   "figure_4n_reported_ratio": float(PANEL_N_PCT[i-2]) / 100,
                   "on_cells_ms": f"{PANEL}b-c!B{i}:D{i}",
                   "off_cells_ms": f"{PANEL}e-f!B{i}:D{i}",
                   "reported_cells_ratio_sd_n": f"{PANEL}h-i!B{hrow}:D{hrow}",
                   **{f"on_{c}_ms": x for c, x in zip("BCD", on)},
                   **{f"off_{c}_ms": x for c, x in zip("BCD", off)},
                   "on_mean_ms": mean(on), "off_mean_ms": mean(off),
                   "reported": rep["mean"], "reported_sd": rep["sd"], "reported_n": rep["n"],
                   "ratio_of_means": rom, "paired_min": min(paired), "paired_max": max(paired)}
            for name, value in [("ratio_of_means", rom), *[("paired_" + p, v) for p, v in zip(PAIRINGS, paired)]]:
                row[name] = value
                row[name + "_abs_diff"], row[name + "_diff_pct"] = deltas(value, rep["mean"])
            j, pka, javg, jsd, n = fit_rows.get(ligand, (None,) * 5)
            row.update(j_cells_pka_ratio_sd_n=f"{PANEL}j!B{j}:E{j}" if j else None,
                       free_pka=pka, j_reported=javg, j_sd=jsd, j_n=n,
                       j_sd_div_sqrt_n=jsd / sqrt(n) if j else None)
            rows.append(row)
        return rows
    finally:
        book.close()


def write_audit(rows, path):
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: published(v) if isinstance(v, float) else v for k, v in row.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DATA)
    parser.add_argument("--output", type=Path, default=DATA / "ratio-audit.csv")
    args = parser.parse_args()
    write_audit(read_audit(args.data_dir), args.output)


if __name__ == "__main__":
    main()

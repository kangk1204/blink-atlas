"""Extract the two Steen et al. benchmarking sheets, with per-replicate cells."""
import argparse
import hashlib
import json
from pathlib import Path
from statistics import mean, stdev
from urllib.request import urlopen

from openpyxl import load_workbook
from openpyxl.utils import get_column_letter

DATA = Path(__file__).resolve().parents[1] / "docs" / "data"
BASE = "https://pmc-oa-opendata.s3.amazonaws.com/PMC11399092.1/41592_2024_2374_MOESM"
BOOKS = {
    "origami": (2, 5, "8b05ae442351bae4e7496ed39e951d7f66897495d0d6c474edac09be1892c408"),
    "npc": (3, 6, "d619405a80864250b570ac1198b08457b3dbfbaae80d1c643d3b03570eb09940"),
}
# Directory keys are the stable acquisition identifiers in the deposited archive.
DYES = {
    "Atto488": ("atto488", "Atto488", 488), "AF488": ("alexafluor488", "AF488", 488),
    "AS488": ("abberiorstar488", "Abberior STAR 488", 488), "CF488": ("cf488a", "CF488A", 488),
    "Cy3b": ("cy3b", "Cy3B", 560), "Atto565": ("atto565", "Atto565", 560),
    "JF585": ("janeliafluor585", "JF585", 560), "CF568": ("cf568", "CF568", 560),
    "Cy5": ("cy5", "Cy5", 640), "Cy5b": ("cy5b", "Cy5B", 640),
    "Atto643": ("atto643", "Atto643", 640), "Atto647N": ("atto647n", "Atto647N", 640),
    "Atto655": ("atto655", "Atto655", 640), "CF640R": ("cf640r", "CF640R", 640),
    "CF660R": ("cf660r", "CF660R", 640), "AS635p": ("abberiorstar635p", "Abberior STAR 635P", 640),
    "JF646": ("janeliafluor646", "JF646", 640), "AF647": ("alexafluor647", "AF647", 640),
}
DYES["CF488A"] = DYES["CF488"]
METRICS = {
    "nenaNm": ("NeNA precision", "nm"), "centerPhotons": ("Center-frame photons", "photons"),
    "centerSbr": ("Signal-to-background ratio", "ratio"), "brightTimeS": ("Bright time", "s"),
    "dropPct": ("Detection loss", "%"), "destroyedPct": ("Docking sites destroyed", "%"),
    "specificity": ("Specificity relative to Cy3B", "ratio"),
}
# Start rows of the repeated observations, before the published mean and sample SD.
ROWS = {"origami": {"nenaNm": 7, "centerPhotons": 12, "brightTimeS": 17, "centerSbr": 22,
                   "dropPct": 27, "destroyedPct": 32},
        "npc": {"nenaNm": 2, "centerPhotons": 8, "dropPct": 14, "brightTimeS": 20,
                "centerSbr": 26, "specificity": 32}}


def extract(directory=DATA):
    catalog = {}
    sources = {}
    for system, (figure, number, sha) in BOOKS.items():
        url = BASE + str(number) + "_ESM.xlsx"
        path = directory / f"steen-figure{figure}.xlsx"
        directory.mkdir(parents=True, exist_ok=True)
        if path.exists():
            blob = path.read_bytes()
        else:
            with urlopen(url, timeout=45) as response:
                blob = response.read()
        if hashlib.sha256(blob).hexdigest() != sha:
            raise ValueError(f"Source digest changed: {url}")
        if not path.exists():
            path.write_bytes(blob)
        book = load_workbook(path, read_only=True, data_only=True)
        try:
            sheet = book.active
            values = list(sheet.values)
            sources[system] = {"url": url, "sha256": sha, "figure": figure, "sheet": sheet.title,
                               "local": "data/" + path.name}
        finally:
            book.close()
        seen = set()
        for col, alias in enumerate(values[0][1:], start=1):
            key, name, wave = DYES[alias]
            if key in seen:
                raise ValueError(f"Duplicate dye column {alias}")
            seen.add(key)
            dye = catalog.setdefault(key, {"id": key, "name": name, "wavelength": wave, "systems": {}})
            letter = get_column_letter(col + 1)
            data = {"column": letter, "worksheetDye": alias, "metrics": {},
                    "concentrationPm": values[37 if system == "origami" else 38][col]}
            if system == "origami":
                data["buffer"] = values[5][col]
            for m, first in ROWS[system].items():
                nrows = 3 if system == "origami" else 4
                scale = 100 if m in ("dropPct", "destroyedPct") else 1
                obs = [{"repeat": i+1, "value": scale*values[first+i-1][col], "cell": f"{letter}{first+i}"}
                       for i in range(nrows) if values[first+i-1][col] is not None]
                vals = [o["value"] for o in obs]
                published = values[first+nrows-1][col]
                sd = values[first+nrows][col]
                # A few reference-only cells are empty. Keep the actual sample size.
                data["metrics"][m] = {"points": obs, "mean": scale*published if published is not None else None,
                                      "sd": scale*sd if sd is not None else None,
                                      "meanCell": f"{letter}{first+nrows}", "computedMean": mean(vals) if vals else None,
                                      "computedSd": stdev(vals) if len(vals)>1 else None}
            dye["systems"][system] = data
        if len(seen) != 18:
            raise ValueError("Expected eighteen distinct dye columns")
    return {"schema": "dna-paint-palette/1", "article": "https://doi.org/10.1038/s41592-024-02374-8",
            "record": "https://zenodo.org/records/10960858", "sources": sources,
            "metrics": {k: {"label": label, "unit": unit} for k,(label,unit) in METRICS.items()},
            "dyes": list(catalog.values()),
            "notes": ["Fig. 2 labels CF488A as CF488. The same dye has directory cf488a and is labeled CF488A in Fig. 3.",
                      "Fig. 2 rows 32–34 are the three replicate docking-destruction observations. The sheet repeats the label Dest 3 on rows 33 and 34; the original cell addresses are retained.",
                      "Cy3B has an additional NPC acquisition. Its contribution varies by metric."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA)
    args = parser.parse_args()
    obj = extract(args.data)
    (args.data / "palette.json").write_text(json.dumps(obj, indent=2, allow_nan=False)+"\n", encoding="utf8", newline="\n")

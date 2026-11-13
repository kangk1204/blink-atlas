"""Rebuild the external eligibility ledger and verify the Fig. 4 uncertainty labels.

Requires Poppler's pdftotext and pdftoppm on PATH. Original files are cached by SHA-256.
"""

import argparse
import hashlib
from io import BytesIO
import json
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory
from urllib.request import Request, urlopen
from xml.etree import ElementTree

import numpy as np
from PIL import Image

from analysis.analogue_benchmark import benchmark, predict, score_predictions
from analysis.dyes import DATA
from analysis.probe_design import published
from analysis.ratio_audit import read_audit
from analysis.ratio_sensitivity import roster
from tools.extract_analogue_sources import pages


HOLLAND = "https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/"
UNO_SI = ("https://static-content.springer.com/esm/art%3A10.1038%2Fnchem.2002/MediaObjects/"
          "41557_2014_BFnchem2002_MOESM177_ESM.pdf")
UNO_TABLE = "https://www.nature.com/articles/nchem.2002/tables/1"
CHEN_SI = "https://ndownloader.figstatic.com/files/47124907"
FILES = {
    "holland.xml": (HOLLAND + "PMC13167466.1.xml", "b7de30fd6787ebb86aeec5b1d834c02905a6cea386e30ab1836f8bde476fb41b"),
    "holland.jpg": (HOLLAND + "41592_2026_3062_Fig6_ESM.jpg", "fde64e263e421570a508d385d24c334049610be12fa637ca6c36f576d3e7d4ea"),
    "uno.pdf": (UNO_SI, "e3766da0c4e314766309103715bf46355e69fb78c33f82ace9afc8bdeb929651"),
    "chen.pdf": (CHEN_SI, "ac0cd7c3b700d2f0fe8f907c1d4951a6ec712aedb0301b700ef587ee100b094f"),
}


def original(cache, name):
    path = cache / name
    url, digest = FILES[name]
    if path.exists():
        data = path.read_bytes()
    else:
        req = Request(url, headers={"User-Agent": "Mozilla/5.0", "Referer": "https://acs.figshare.com/"})
        with urlopen(req, timeout=90) as response:
            data = response.read()
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f"Unexpected source bytes for {name}")
    if not path.exists():
        path.write_bytes(data)
    return path


def figure_check(xml, jpg, rows):
    caption = " ".join(ElementTree.parse(xml).find('.//fig[@id="Fig6"]/caption').itertext())
    for text in ("error bars indicate ±s.e.m", "n  = 3 independent cellular imaging experiments",
                 "~100-fold lower concentrations", "upper-bound dark-state lifetimes"):
        if text not in caption:
            raise ValueError(f"Caption has changed at {text}")
    img = np.array(Image.open(jpg).convert("L"))
    ticks = [594, 712]  # Published 757 × 1113 JPEG, ticks at 0.0015 and 0.0000.
    checks = []
    for dye, x, top, bottom in [("JFX650b", 548, 611, 643), ("JF639b", 538, 622, 644)]:
        for y in (top, bottom):
            if img[y-1:y+2, x-2:x+3].min() >= 110:
                raise ValueError("Could not locate the specified error-bar cap")
        record = next(r for r in rows if r["dye"] == dye)
        half = (bottom-top)/2 * .0015/(ticks[1]-ticks[0])
        checks.append({"dye": dye, "xPixel": x, "capYPixels": [top, bottom],
                       "approxHalfBar": float(format(half, ".3g")), "workbookSD": record["j_sd"],
                       "sdOverSqrtN": published(record["j_sd_div_sqrt_n"])})
    return {"pairedPreparationIds": "No cross-acquisition preparation ID or pairing key in the workbook or Fig. 4 methods.",
            "onExposureMs": 2, "offExposureMs": 25, "offConcentrationApproxFoldLower": 100,
            "offInterpretation": "Upper bound under these acquisition conditions.",
            "panelHi": "SD in both workbook and caption.", "panelJ": "SD in workbook D2, SEM in the caption.",
            "sharedSdCells": sum(r["j_sd"] is not None and r["j_sd"] == r["reported_sd"] for r in rows),
            "axisTicksYPixels": ticks, "axisTickRatios": [.0015, 0], "approxCapReadingPixels": 1,
            "visibleWhiskers": checks,
            "panelNCheck": "In the rendered 4n table, 8HTL is 0.029% and 4HTL reads 0.10 ± 0.2%. For 4HTL the underlying 4h–i workbook SD is 0.0236 percentage points. The mean-based comparisons use the higher-precision source cells.",
            "interpretation": "Two resolved whiskers agree to pixel precision with SD divided by sqrt(3). Smaller bars are unresolved. The scores use the means."}


def render(pdf, page):
    with TemporaryDirectory() as temp:
        stem = Path(temp) / "figure"
        subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-scale-to", "2500", "-png",
                        "-singlefile", str(pdf), str(stem)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        with Image.open(BytesIO(stem.with_suffix(".png").read_bytes())) as img:
            if img.size != (1768, 2500):
                raise ValueError("Unexpected coordinate system for Chen plots")
            return np.array(img.convert("RGB"), dtype=int)


def optical_midpoint(img, bounds, ticks):
    """Read the descending half-height between the visible extrema of the free-dye titration."""
    x0, x1, y0, y1 = bounds
    red, green, blue = img[y0:y1+1, x0:x1+1, :].transpose(2, 0, 1)
    curve = (red > 150) & (red-green > 85) & (red-blue > 85)
    xs = np.flatnonzero(curve.any(axis=0))
    ys = np.array([np.median(np.flatnonzero(curve[:, x])) for x in xs])
    mid = (ys.min() + ys.max())/2
    i = np.flatnonzero((ys[:-1] < mid) & (ys[1:] >= mid))[-1]
    cross = x0 + np.interp(mid, ys[i:i+2], xs[i:i+2])
    return round(3 + 5 * (cross-ticks[0]) / (ticks[1]-ticks[0]), 1)


def chen_duties(img):
    # At 2500 px tall, 0.1 and 0.01 tick marks are at y=1251 and 1381.
    powers = [.10, .27, .65, 1.49, 3.93]
    x_pixels = [1156, 1234, 1303, 1368, 1445]
    conditions = {"PBS": (7.4, [1344, 1370, 1436, 1462, 1476]),
                  "MES": (4.7, [1199, 1200, 1201, 1230, 1262])}
    rows = []
    for buffer, (ph, centers) in conditions.items():
        for power, x, center in zip(powers, x_pixels, centers, strict=True):
            red, green, blue = img[center-7:center+8, x-1:x+2, :].transpose(2, 0, 1)
            green_line = (red > 90) & (red-blue > 36) & (abs(red-green) < 15)
            orange_line = (red > 150) & (red-green > 34) & (green-blue > 33)
            y = np.flatnonzero((green_line if buffer == "PBS" else orange_line).any(axis=1))
            if not len(y):
                raise ValueError(f"Missing S7c point for {buffer}, {power}")
            center_y = center-7 + (y.min()+y.max())/2
            duty = 10 ** (-1-(center_y-1251)/130)
            rows.append({"dye": "SMI-NIR1", "conjugate": "Biotin", "buffer": buffer, "pH": ph,
                         "wavelengthNm": 721, "powerKwCm2": power, "dutyFractionApprox": float(f"{duty:.2g}"),
                         "dutyOddsApprox": float(f"{duty/(1-duty):.2g}"),
                         "locator": "Chen SI Fig. S7c, p. S15", "endpoint": "pre-bleach single-molecule duty fraction"})
    return rows


def assess(train, row):
    """Preserve missing observations and range abstentions in each model's result."""
    available = row.get("pKa")
    rule = predict(train, available) if available is not None else {"status": "no exact free-parent pKa"}
    measured = row.get("meanDwellRatio")
    if rule["status"] == "scored" and measured is not None:
        errors = score_predictions(rule["predictions"], measured)["errors"]
    else:
        errors = {m: None for m in ("local", "nearest", "regression")}
    row.update(rangeStatus=rule["status"], predictions=rule.get("predictions"), absLogErrors=errors)
    return row


def uno_rows(pdf, train):
    supplementary = pages(pdf, "pdftotext")
    if "Laser flash photolysis" not in supplementary[35] or "single-molecule" not in supplementary[25].lower():
        raise ValueError("Expected Uno supplementary panels are missing")
    # Transcription of publisher Table 1. Figure S2 contains the corroborating titrations,
    # including a separate fluorescence readout of 5.9 for HMSiR.
    table = {"HMSiR": (5.8, .245), "HMAcRG": (5.3, .0218), "AMSiR": (4.7, .0000107),
             "AMAcRG": (3.7, .00000309), "MMSiR": ("<2", "<5e-8"), "HySOx": ("<2", "<5e-8")}
    titrations = dict(re.findall(r"\b([0-9]?[A-Z][A-Za-z0-9]+) \(([<>]?\s*[\d.]+)\)", supplementary[33]))
    if set(table) - titrations.keys() or len(titrations) != 14:
        raise ValueError("Unexpected compound roster in Uno Fig. S2")
    result = []
    for dye, text in titrations.items():
        value = text.replace(" ", "")
        optical = float(value) if value[0].isdigit() else value
        pka, tau = table.get(dye, (optical, None))
        if dye in table and dye != "HMSiR" and pka != optical:
            raise ValueError("Table 1 disagrees with the SI titration")
        row = {"study": "Uno 2014", "dye": dye, "pKa": pka if isinstance(pka, float) else None,
               "pKaBound": pka if isinstance(pka, str) else None, "siFluorescencePKcycl": optical,
               "flashOpenLifetimeSeconds": tau, "meanDwellRatio": None,
               "locator": ("Table 1 and SI Figs. S2, S4" if tau is not None else "SI Fig. S2"),
               "assay": "Solution titration and, where listed, 308 nm flash photolysis at pH 7.4",
               "reason": ("Only open-state relaxation is paired with the titration. An independently measured mean OFF time is missing."
                          if tau is not None else "Free-parent titration only in the screened panel. No paired mean dwell ratio."),
               "alsoInHollandNine": dye == "HMSiR"}
        result.append(assess(train, row))
    return result


def uno_lifetime(records):
    # These are measured, uncensored response times in seconds. The existing three
    # algorithms accept any positive log response. Fit solely to this assay's dyes.
    train = [{"id": r["dye"], "dye": r["dye"], "pKa": r["pKa"], "ratio": r["flashOpenLifetimeSeconds"]}
             for r in records if isinstance(r.get("flashOpenLifetimeSeconds"), float) and r["pKa"] is not None]
    result = benchmark("uno-flash-lifetime", "Uno · laser flash lifetime", "uno", "transfer-local", train)
    result.update(endpoint="Open-state relaxation after a 308 nm pulse at pH 7.4",
                  responseUnit="seconds", errorUnit="absolute log10 error", source=UNO_TABLE,
                  chemistry="Aminomethyl and hydroxymethyl rhodamines, including slow components of biphasic decays.",
                  censored=["MMSiR", "HySOx"], training="Other exact-value dyes in Uno Table 1, whole parent held out",
                  rangePolicy="Withhold a prediction outside the remaining training range")
    return result


def build(cache, directory):
    cache.mkdir(parents=True, exist_ok=True)
    files = {k: original(cache, k) for k in FILES}
    rows = read_audit(directory)
    train = roster([r for r in rows if r["free_pka"] is not None], "j_reported")
    records = uno_rows(files["uno.pdf"], train)
    pdf = files["chen.pdf"]
    supp = pages(pdf, "pdftotext")
    if "Figure S7." not in supp[14] or "Figure S3." not in supp[10]:
        raise ValueError("Unexpected Chen supplementary figures")
    titrations = render(pdf, 11)
    bounds = [(672, 802, 280, 392), (616, 803, 746, 888), (643, 805, 1216, 1341), (638, 792, 1690, 1817)]
    ticks = [(672, 801), (616, 809), (643, 812), (638, 790)]
    for i, (box, tick) in enumerate(zip(bounds, ticks), 1):
        mid = optical_midpoint(titrations, box, tick)
        records.append(assess(train, {"study": "Chen 2024", "dye": f"SMI-NIR{i}", "pKa": None,
            "opticalMidpointApproxPH": mid, "meanDwellRatio": None, "alsoInHollandNine": False,
            "assay": "Free dye in citrate-phosphate buffer with cosolvent, single-molecule biotin conjugate at 721 nm",
            "locator": f"SI Fig. S3{('A','C','E','G')[i-1]}, p. S11, and Fig. S{7 if i == 1 else 6}, p. S15",
            "reason": "The plotted fluorescence half-height exceeds the Holland pKa range. An exact fitted parent pKa is unavailable in the inspected SI. " +
                      ("Duty fraction is quantified at five powers in each of two buffers." if i == 1 else
                       "Figure S6 provides a representative intensity trace at 3.93 kW/cm².")}))
    counts = {s: len({r["dye"] for r in records if r["study"] == s}) for s in ("Uno 2014", "Chen 2024")}
    return {"schema": "blink-atlas/ratio-evidence/1", "screened": "2026-11-13",
            "sources": {k: {"url": url, "sha256": digest} for k, (url, digest) in FILES.items()},
            "articles": {"Uno 2014": {"url": "https://doi.org/10.1038/nchem.2002", "table": UNO_TABLE,
                        "access": "Publisher Table 1, figure pages and complete supplement."},
                         "Chen 2024": {"url": "https://doi.org/10.1021/acs.analchem.4c02445",
                        "deposit": "https://doi.org/10.1021/acs.analchem.4c02445.s001",
                        "access": "Publisher abstract and complete author-deposited supplement. Subscription article body unavailable in this review."}},
            "primaryProtocol": "Independent fixed COS-7 H2B–HaloTag under the registered acquisition",
            "countsByStudy": counts, "distinctScreened": sum(counts.values()),
            "previouslyInHollandNine": sum(r["alsoInHollandNine"] for r in records),
            "matchedScored": 0,  # Both screened acquisition designs use other specimens.
            "frozenTransferScored": sum(r["absLogErrors"]["local"] is not None for r in records),
            "records": records, "unoLifetime": uno_lifetime(records), "chenByCondition": chen_duties(render(pdf, 15)),
            "chenReadout": "Free-parent panel S3 descending half-heights read to 0.1 pH. Panel S7c uses a logarithmic duty axis. Each of ten readouts is from the same parent. Values are approximate.",
            "uncertaintyAudit": figure_check(files["holland.xml"], files["holland.jpg"], rows)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, default=DATA)
    parser.add_argument("--output", type=Path, default=DATA / "ratio-evidence.json")
    args = parser.parse_args()
    result = build(args.cache, args.data_dir)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                           encoding="utf-8", newline="\n")
    print(f"{result['distinctScreened']} parent dyes screened. {result['matchedScored']} matched H2B outcomes.")


if __name__ == "__main__":
    main()

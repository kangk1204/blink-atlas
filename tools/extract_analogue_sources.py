"""Recover the paired public measurements from checksum-pinned primary sources."""

import argparse
import csv
import hashlib
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

from analysis.dyes import DATA, ROOT


FILES = {
    "remmel-small.pdf": ("https://pure.mpg.de/rest/items/item_3489427_3/component/file_3506946/content",
                         "c2c13558bd569077721779f0a02592b989b478f22d4c94bb514d77cd71dd114b"),
    "bio-pre-sup.pdf": ("https://www.biorxiv.org/content/biorxiv/early/2022/08/29/2022.08.29.505670/DC1/embed/media-1.pdf?download=true",
                        "cc10cbc8e11577d1bfdf05d4218b73911ee9a10099db5a174e3f095d60936dbe"),
    "jacs2023-ja2c11395_si_001.pdf": ("https://ndownloader.figstatic.com/files/39374590",
                                     "46de3dc20928445359a04e483ddf8f2bc09955dec775a628919dcb165e31f652"),
    "tyson-archive.xml": ("https://pmc-oa-opendata.s3.amazonaws.com/PMC8393207.1/PMC8393207.1.xml",
                          "1ffd03a16188fc3d3207e1932ab7ce0bd9dbbd4902edda66920f94b2e2f1a386"),
}
FIELDS = ["study", "dye_id", "dye", "condition", "pka", "bound_pka", "on_ms", "off_s", "duty", "duty_sd"]


def sources(cache, names=None):
    cache.mkdir(parents=True, exist_ok=True)
    files = []
    for name, (url, digest) in FILES.items():
        if names is not None and name not in names:
            continue
        path = cache / name
        if path.exists():
            body = path.read_bytes()
        else:
            headers = {"User-Agent": "Mozilla/5.0"}
            if url.startswith("https://ndownloader.figstatic.com/"):
                headers["Referer"] = "https://acs.figshare.com/"
            with urlopen(Request(url, headers=headers), timeout=60) as response:
                body = response.read()
        if hashlib.sha256(body).hexdigest() != digest:
            raise ValueError(f"Source checksum mismatch for {name}")
        if not path.exists():
            path.write_bytes(body)
        files.append(path)
    return files


def pages(pdf, executable):
    result = subprocess.run([executable, "-layout", "-enc", "UTF-8", str(pdf), "-"], capture_output=True, check=True)
    return result.stdout.decode("utf-8").split("\f")


def remmel(page):
    if "Table 1. Photophysical properties" not in page:
        raise ValueError("Remmel Table 1 missing")
    rows = []
    for line in page.splitlines():
        fields = line.split()
        if not fields or fields[0] not in {str(i) for i in range(1, 7)}:
            continue
        # The typeset equilibrium constants may contain a thousands separator space.
        if len(fields) not in (14, 15):
            raise ValueError(f"Unexpected Remmel row: {line}")
        n = fields[0]
        rows.append(dict(study="remmel", dye_id="remmel-"+n, dye="Compound "+n, condition="IgG",
                         pka=fields[5], on_ms=fields[-6], off_s=fields[-5]))
    if [r["dye_id"] for r in rows] != [f"remmel-{n}" for n in range(1, 7)]:
        raise ValueError("Expected six numbered Remmel compounds")
    return rows


def remmel_conjugates(page, igg):
    if "Table S1. Photophysical properties for compound 4" not in page:
        raise ValueError("Remmel conjugate table missing")
    reference = next(r for r in igg if r["dye_id"] == "remmel-4")
    rows, repeated = [], False
    for line in page.splitlines():
        cells = line.replace(",", ".").split()
        if len(cells) != 8 or cells[0] not in ("Antibody", "Halotag", "Nanobody"):
            continue
        if cells[0] == "Antibody":
            repeated = all(float(cells[i]) == float(reference[key]) for i, key in ((2, "on_ms"), (3, "off_s")))
        else:
            condition = "HaloTag7" if cells[0] == "Halotag" else "Nanobody"
            rows.append(dict(reference, condition=condition, on_ms=cells[2], off_s=cells[3]))
    if not repeated or {r["condition"] for r in rows} != {"HaloTag7", "Nanobody"}:
        raise ValueError("Incomplete Remmel conjugate comparison")
    return rows


def zheng(spectra, photophysics):
    if "Table S1." not in spectra or "Table S2." not in photophysics:
        raise ValueError("Zheng source table header missing")
    rows = []
    for match in re.finditer(r"^\s*Water\s+([^\n]+)\n\s*(STMR|SRhB|PySRh|PiSRh)\s*$", spectra, re.MULTILINE):
        unbound, bound = match[1].split()[-2:]
        dye = match[2]
        blocks = [b for b in re.split(r"\n\s*\n", photophysics)
                  if re.search(r"^"+dye+r"\s+7\.4\s*$", b, re.MULTILINE)]
        if len(blocks) != 1:
            raise ValueError(f"Ambiguous Table S2 block for {dye}")
        for line in blocks[0].splitlines():
            cells = line.split()
            if not cells or cells[0] not in ("500", "1000", "1500", "2000"):
                continue
            fraction, sd = cells[3].split("±")
            rows.append(dict(study="zheng", dye_id=dye.lower(), dye=dye, condition=cells[0],
                             pka=unbound, bound_pka=bound, duty=fraction, duty_sd=sd))
    if len(rows) != 16 or len({(r["dye"], r["condition"]) for r in rows}) != 16:
        raise ValueError("Expected four measurements for each of four sulfonamides")
    return rows


def tyson(xml):
    text = " ".join(" ".join("".join(p.itertext()).split()) for p in ET.parse(xml).iter("p"))
    pka = re.search(r"pK\s*cycl\s*=\s*([\d.]+)\s+vs\s+6.9", text)
    ton = re.search(r"has an ON time of ([\d.]+) ms at pH 7.4", text)
    toff = re.search(r"an OFF time of ([\d.]+) s", text)
    if not all((pka, ton, toff)):
        raise ValueError("Could not identify Yale676sb's reported kinetics and pKa")
    return [dict(study="tyson", dye_id="yale676sb", dye="Yale676sb", condition="Glass",
                 pka=pka[1], on_ms=ton[1], off_s=toff[1])]


def write_observations(path, rows):
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def verify_igg(cache, pdftotext):
    from analysis.analogue_benchmark import benchmark, observations
    from analysis.ratio_audit import read_audit
    from analysis.ratio_sensitivity import roster

    article, = sources(cache, {"remmel-small.pdf"})
    with TemporaryDirectory() as tmp:
        file = Path(tmp) / "igg.csv"
        write_observations(file, remmel(pages(article, pdftotext)[3]))
        targets = observations(file)["remmel", "IgG"]
    train = roster([r for r in read_audit() if r["free_pka"] is not None], "j_reported")
    metrics = benchmark("remmel-igg", "", "remmel", "transfer-frozen", targets, train)["summary"]
    scores = [round(metrics["mae"][m], 3) for m in ("local", "nearest", "regression")]
    if metrics["scored"] != 4 or scores != [.293, .307, .208]:
        raise ValueError(f"Frozen IgG benchmark changed: {scores}")
    print(f"Remmel IgG, {metrics['scored']}/{metrics['total']} dyes. Local, nearest, regression: {scores}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=ROOT / ".cache/analogue-sources")
    parser.add_argument("--pdftotext", default="pdftotext")
    parser.add_argument("--output", type=Path, default=DATA / "analogue-observations.csv")
    parser.add_argument("--verify-igg", action="store_true",
                        help="Rebuild the legacy frozen IgG score directly from the primary Remmel and Holland files")
    args = parser.parse_args()
    if args.verify_igg:
        verify_igg(args.cache, args.pdftotext)
        return
    article, conjugates, supplement, xml = sources(args.cache)
    supp_pages = pages(supplement, args.pdftotext)
    igg = remmel(pages(article, args.pdftotext)[3])
    rows = igg + remmel_conjugates(pages(conjugates, args.pdftotext)[7], igg)
    rows += zheng(supp_pages[8], supp_pages[15]) + tyson(xml)
    write_observations(args.output, rows)
    print(f"Wrote {len(rows)} measurements from three studies to {args.output}")


if __name__ == "__main__":
    main()

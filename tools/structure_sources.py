"""Check the source files behind the curated graphs and numerical labels."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from pathlib import Path
import re
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

from openpyxl import load_workbook

from analysis.probe_design import BASE, HASHES, SHEETS, SI_FILE, SI_HASH, file_name
from analysis.structure_catalog import catalog, quantum_rows


TY = "https://pmc-oa-opendata.s3.amazonaws.com/PMC8393207.1/"
TA = "https://pmc-oa-opendata.s3.amazonaws.com/PMC9814528.1/"
U = "https://static-content.springer.com/"
FILES = {
    "holland-si.pdf": (BASE + SI_FILE, SI_HASH),
    "holland-4j.xlsx": (BASE + file_name(10), HASHES[10]),
    "holland-2jk.xlsx": (BASE + file_name(9), HASHES[9]),
    "tyson2.jpg": (TY + "oc1c00670_0002.jpg", "14909e64b0af0e452d42f49d1c1888557d09e48112591b35b5ac0e3bcacede5d"),
    "tyson3.jpg": (TY + "oc1c00670_0003.jpg", "90656c5ee2543aaf2de170398137a60dd7dd41c584365d433fcbab9f8d564e07"),
    "tyson-si.pdf": (TY + "oc1c00670_si_001.pdf", "8622f18d41247c4cedffcd3aea38196c1ba15ebdc7039cb75332202f78f86060"),
    "uno-fig2.jpg": (U + "image/art%3A10.1038%2Fnchem.2002/MediaObjects/41557_2014_Article_BFnchem2002_Fig2_HTML.jpg",
                     "946685c7c35883f9e6618ba1e56ce87862f52bfc5130550169c202a3abd849ab"),
    "uno-si.pdf": (U + "esm/art%3A10.1038%2Fnchem.2002/MediaObjects/41557_2014_BFnchem2002_MOESM177_ESM.pdf",
                   "e3766da0c4e314766309103715bf46355e69fb78c33f82ace9afc8bdeb929651"),
    "tach.xml": (TA + "PMC9814528.1.xml", "b750d1ed586e1c799ab22a1003da21c2f4440b26dcc245d8b565b56ecbed0c84"),
    "tach-si.pdf": (TA + "42004_2020_326_MOESM1_ESM.pdf", "06265370d8e29c20a94af52a7937d5a14dc8518d5003b3a428150dfdd5bfb3bb"),
    "sakabe-si.pdf": ("https://ndownloader.figstatic.com/files/4096537",
                      "d3e3c73f7cfe324393c7bdef9bd4b917b84d4c3fc35fb2259b958c85a343c67b"),
    "iwatate-2018-si.pdf": ("https://ndownloader.figstatic.com/files/10257447",
                            "d32935c36d61fafb4e394ada275ea433e68164e81fc7ee0101c8a4294a67c24f"),
}


def verify_files(cache, sources):
    cache.mkdir(parents=True, exist_ok=True)

    def obtain(item):
        name, (url, digest) = item
        dest = cache / name
        if dest.exists():
            body = dest.read_bytes()
        else:
            headers = {"User-Agent": "Mozilla/5.0", "Referer": "https://acs.figshare.com/"}
            with urlopen(Request(url, headers=headers), timeout=50) as response:
                body = response.read()
        if hashlib.sha256(body).hexdigest() != digest:
            raise ValueError("Source checksum mismatch: " + name)
        if not dest.exists():
            dest.write_bytes(body)

    with ThreadPoolExecutor(max_workers=3) as fetch:
        list(fetch.map(obtain, sources.items()))


def check(cache):
    verify_files(cache, FILES)
    workbook = load_workbook(cache / "holland-4j.xlsx", read_only=True, data_only=True)
    try:
        table = list(workbook[SHEETS[10]].values)[2:8]
        dut = [r for r in catalog() if r["study"] == "Holland" and r["status"] == "exact"]
        assert [(r[0], r[1]) for r in table] == [(r["locator"].split("compound ")[1].split(";")[0] + "HTL",
                                                    float(r["reported"])) for r in dut]
    finally:
        workbook.close()
    t1, t2 = ET.parse(cache / "tach.xml").getroot().findall(".//table-wrap")[:2]
    data = [["".join(c.itertext()) for c in tr] for tr in t1.iter("tr") if len(tr.findall("td")) == 9]
    expected = [(name, observed, qm) for name, _, observed, qm in quantum_rows()]
    def numbers(s):
        return [float(v) for v in re.findall(r"\d+\.\d+", s)]
    assert [(r[0], float(r[6]), numbers(r[7])[-1]) for r in data] == expected
    t2_data = {"".join(tr[0].itertext()): numbers("".join(tr[1].itertext())) for tr in t2.iter("tr") if len(tr.findall("td")) > 2}
    assert list(t2_data["F"] + t2_data["Me"] + t2_data["CF3"]) == [8, 6.6, 8.2, 5.3]
    direct = {r["dye"]: float(r["reported"]) for r in catalog() if r["study"] == "Tachibana"}
    assert [direct[k] for k in ("3FHMRG", "3MHMRG", "5MHMRG", "3CF3HMRG")] == [
        *t2_data["F"], *t2_data["Me"], *t2_data["CF3"]]
    print(f"{len(FILES)} primary source files verified. Holland 6 exact workbook rows and Tachibana 17+4 table entries agree.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    check(parser.parse_args().cache)

"""Fetch the original Chi figure pages and source SI for the expanded exact-label panel."""

import argparse
from pathlib import Path

from analysis.transfer_sources import CHI_FIG2, CHI_FIG3
from tools.structure_sources import FILES, verify_files


SOURCES = {
    "chi-original-fig2.png": (CHI_FIG2, "aa0db07ed2072cc1a6f0fdea481663d01ef5aa3ce39abc69732ef6bc37c24ab5"),
    "chi-original-fig3.png": (CHI_FIG3, "2ce09a11c7476eb136bea1b992ad3002c8d1197cc5030611aac6c2accc26636d"),
    "sakabe-si.pdf": FILES["sakabe-si.pdf"],
}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    verify_files(p.parse_args().cache, SOURCES)
    print("Verified original Chi structure and absorbance panels and Sakabe SI spectra (three files).")

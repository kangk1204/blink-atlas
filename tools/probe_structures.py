"""Render the closed parent structures directly from Holland et al., Supplementary Note 2."""

import argparse
import hashlib
from pathlib import Path
from urllib.request import urlopen

import pypdfium2 as pdfium

from analysis.dyes import ROOT, DYE_NAMES
from analysis.probe_design import BASE, CHEMISTRY, SI_FILE, SI_HASH


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, default=ROOT / ".cache" / SI_FILE)
    parser.add_argument("--output", type=Path, default=ROOT / "docs/media/probe-structures")
    args = parser.parse_args()
    if args.pdf.is_file():
        source = args.pdf.read_bytes()
    else:
        with urlopen(BASE + SI_FILE, timeout=60) as res:
            source = res.read()
    if hashlib.sha256(source).hexdigest() != SI_HASH:
        raise ValueError("Supporting information checksum mismatch")
    if not args.pdf.is_file():
        args.pdf.parent.mkdir(parents=True, exist_ok=True)
        args.pdf.write_bytes(source)
    args.output.mkdir(parents=True, exist_ok=True)
    with pdfium.PdfDocument(source) as pdf:
        for i, spec in enumerate(CHEMISTRY, 1):
            page_number, (left, top, right, bottom) = spec[-2:]
            page = pdf[page_number - 1]
            bitmap = page.render(scale=3)
            raster = bitmap.to_pil()
            w, h = raster.size
            box = (round(left*w/521), round(top*h/674), round(right*w/521), round(bottom*h/674))
            raster.crop(box).convert("L").save(args.output / f"{i}.png", optimize=True)
            bitmap.close()
            page.close()
            print(f"{DYE_NAMES[i - 1]}: PDF page {page_number}")


if __name__ == "__main__":
    main()

"""Rebuild the single-page English collaboration brief with ReportLab."""

from __future__ import annotations

import argparse
from pathlib import Path

from reportlab import rl_config
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import Paragraph


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/media/validation-proposal.pdf"
INK = HexColor("#142E30")
TEXT = HexColor("#273E40")
MUTED = HexColor("#445C5D")
TEAL = HexColor("#086761")
LINE = HexColor("#BDCBC8")
TINT = HexColor("#EAF4EF")
LEFT, RIGHT = 43, 569
WIDTH, HEIGHT = letter


def paragraph(canvas, text, x, top, width, font="Helvetica", size=11.0,
              leading=15.9, color=TEXT):
    style = ParagraphStyle(
        "body", fontName=font, fontSize=size, leading=leading,
        textColor=color, alignment=TA_LEFT, splitLongWords=False,
        spaceBefore=0, spaceAfter=0, allowWidows=0, allowOrphans=0,
    )
    block = Paragraph(text, style)
    _, height = block.wrap(width, HEIGHT)
    block.drawOn(canvas, x, top - height)
    return top - height


def label(canvas, text, x, y):
    canvas.setFillColor(TEAL)
    canvas.setFont("Helvetica-Bold", 9)
    canvas.drawString(x, y, text)


def rule(canvas, y):
    canvas.setStrokeColor(LINE)
    canvas.setLineWidth(0.55)
    canvas.line(LEFT, y, RIGHT, y)


def create_pdf(output: Path):
    # The invariant setting also fixes timestamp and document ID between runs.
    rl_config.invariant = 1
    output.parent.mkdir(parents=True, exist_ok=True)
    c = Canvas(str(output), pagesize=letter, invariant=1, pageCompression=1)
    c.setTitle("A matched test of dye selection | Blink Atlas")
    c.setSubject("Proposed collaborative validation for the Lavis Lab")
    c.setAuthor("Blink Atlas")
    c.setCreator("Blink Atlas / ReportLab")
    label(c, "BLINK ATLAS", LEFT, 749)
    c.setFont("Helvetica", 9)
    c.setFillColor(MUTED)
    c.drawRightString(RIGHT, 749, "Research proposal for the Lavis Lab   /   November 2026")
    rule(c, 737)

    c.setFont("Times-Bold", 32)
    c.setFillColor(INK)
    c.drawString(LEFT - 1, 697, "A matched test of dye selection")
    top = paragraph(c,
        "Can one reference movie improve the choice of two spontaneously blinking dyes "
        "for a fixed-cell H2B-HaloTag experiment? Blink Atlas uses measured on-times, "
        "off-times and photon rates to estimate how many isolated signals a camera can record.",
        LEFT, 679, RIGHT - LEFT, size=11.6, leading=16.2)

    box_top = top - 14
    box_height = 120
    c.setFillColor(TINT)
    c.roundRect(LEFT, box_top - box_height, RIGHT - LEFT, box_height, 6,
                stroke=0, fill=1)
    label(c, "First request  /  An existing matched dataset", LEFT + 14, box_top - 20)
    end = paragraph(c,
        "Could you share raw movies or SMAP localization tables for fixed cells from the same "
        "batch and microscope, imaged with at least two of JF630b, JF635b and JF639b? Two "
        "ligands permit a pairwise check. All three permit the two-choice comparison. Include "
        "preparation day, well, field, region area, labeling protocol, laser power, exposure, "
        "frame interval and pixel size. For tables, include frame, x/y, detected photons and "
        "lateral localization uncertainty with units.",
        LEFT + 14, box_top - 30, RIGHT - LEFT - 28, size=10.7, leading=14.0)
    if end < box_top - box_height + 6:
        raise ValueError("Data-request text exceeds reserved space.")

    top = box_top - box_height - 18
    c.setFont("Times-Bold", 19)
    c.setFillColor(INK)
    c.drawString(LEFT, top - 16, "If new imaging is needed")
    top -= 26
    top = paragraph(c,
        "<b>54 movies</b> &nbsp;= &nbsp;3 dyes × 2 labeling densities × 3 independent preparation days "
        "× 3 fields per condition.", LEFT, top, RIGHT - LEFT, size=11.0, leading=15.0) - 8
    top = paragraph(c,
        "Use fixed COS-7 cells expressing H2B-HaloTag. Split each preparation into matched wells, "
        "set sparse and dense labeling in advance, and randomize acquisition order. Keep 647 nm "
        "excitation, 25 ms exposure, frame interval, camera settings and scored acquisition time "
        "constant across the three ligands.", LEFT, top, RIGHT - LEFT) - 13

    rule(c, top)
    top -= 16
    c.setFont("Times-Bold", 19)
    c.setFillColor(INK)
    c.drawString(LEFT, top - 12, "Freeze the choices before scoring")
    top -= 24
    top = paragraph(c,
        "Within each day and density, reserve field 1 of each dye for the pilot. The model sees "
        "only the first 200 frames of reference JF630b. Both strategies receive this same pilot "
        "budget and two candidate slots. Rank only JF630b, JF635b and JF639b. Lock the pair before revealing "
        "fields 2 and 3. The comparison pair is always <b>JF630b + JF635b</b>. Use that fixed pair "
        "when the model abstains and record every abstention.", LEFT, top, RIGHT - LEFT) - 13

    rule(c, top)
    top -= 16
    c.setFont("Times-Bold", 19)
    c.setFillColor(INK)
    c.drawString(LEFT, top - 12, "Endpoint and replication")
    top -= 24
    top = paragraph(c,
        "Count localizations per µm² with lateral uncertainty at or below <b>20 nm</b>, using the "
        "same acquisition time. For per-axis estimates, require both axes to pass. "
        "Within each day and density, average each dye over its two scoring fields and take "
        "the higher yield of the two selected dyes. Report both strategies and their paired "
        "difference for each preparation day. Repeat on a new target or instrument once gains "
        "recur across independent preparations.", LEFT, top, RIGHT - LEFT)

    footer_y = 52
    if top < footer_y + 12:
        raise ValueError(f"Body overlaps the footer (bottom={top:.1f} pt).")
    rule(c, footer_y)
    paragraph(c,
        'Measured photophysics from Holland et al., <i>Nature Methods</i> (2026), Extended Data Fig. 4. '
        '<link href="https://doi.org/10.1038/s41592-026-03062-5" color="#086761">'
        'doi:10.1038/s41592-026-03062-5</link>',
        LEFT, footer_y - 8, RIGHT - LEFT, size=8.7, leading=11, color=MUTED)
    c.showPage()
    c.save()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUTPUT)
    create_pdf(parser.parse_args().out)

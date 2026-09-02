#!/usr/bin/env python3
"""
Supplemental Figure 1. Enrollment flow chart for the methylomic analyses.

Companion generator for Williams et al., "Methylomic Analysis of Nasal Brushings
in Pediatric Acute Respiratory Distress Syndrome: A Pilot Study".

This figure previously had no generator in the repository -- the ORDERING NOTE
in Methyl_PARDS_pipeline.R records Supplemental Figure 1 and Supplemental
Table 1 as hand-maintained. Every coordinate below was recovered from the
previous PDF, so the layout is reproduced exactly. Two deliberate changes:

  1. The title carries the "Supplemental Figure 1." prefix, matching
     Supplemental Figure 2 and Figures 1-6.
  2. pdf.fonttype = 42 embeds TrueType rather than Type 3 fonts. The previous
     PDF used Type 3, which many publishers reject at production and which
     left the text unextractable -- the en dash in the subtitle was dropped on
     copy and pdftotext could not read the date range.

The enrollment counts are collected at the top so they can be checked against
the Results text and Supplemental Table 1 without reading the layout code.

Usage:  python3 make_supplemental_figure_1.py [outdir]
"""
import sys
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42      # TrueType, not Type 3
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"

import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Polygon

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "."

# ---------------------------------------------------------------------------
# Enrollment counts
# ---------------------------------------------------------------------------
N_SCREENED = 170
N_PARDS_ENROLLED, N_CTRL_ENROLLED = 25, 23          # initially enrolled
N_PARDS_INTERIM, N_CTRL_INTERIM = 21, 9             # evaluable at interim
N_PARDS_NEW, N_CTRL_NEW = 21, 0                     # added on continued enrollment
N_PARDS_TOTAL, N_CTRL_TOTAL = 46, 23                # cohort totals
N_PARDS_SPEC, N_PARDS_RESID, N_PARDS_FRESH = 21, 8, 13
N_CTRL_SPEC, N_CTRL_RESID, N_CTRL_PREV = 3, 2, 1

DATE_RANGE = "April 2, 2018 – June 30, 2022"

# ---------------------------------------------------------------------------
# Page geometry, in points, exactly as recovered from the previous figure.
# ---------------------------------------------------------------------------
W, H = 612.0, 684.0
LW = 1.4                       # stroke width
RX, RY = 9.59616, 10.72512     # corner insets (equal fractions of W and H)
CX = 306.0                     # page centre
LX, RX_COL = 120.0744, 491.9256   # left / right column centres
HEAD_L, HEAD_W = 5.6, 5.6      # arrowhead length and full width
TIP_GAP = 1.5652               # gap between arrow tip and the box it points at

fig = plt.figure(figsize=(W / 72.0, H / 72.0))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")


def box(x0, x1, y0, y1):
    """Rounded box with quadratic-Bezier corners, as in the original."""
    v = [(x0 + RX, y0),
         (x1 - RX, y0), (x1, y0), (x1, y0 + RY),
         (x1, y1 - RY), (x1, y1), (x1 - RX, y1),
         (x0 + RX, y1), (x0, y1), (x0, y1 - RY),
         (x0, y0 + RY), (x0, y0), (x0 + RX, y0),
         (0, 0)]
    c = [Path.MOVETO,
         Path.LINETO, Path.CURVE3, Path.CURVE3,
         Path.LINETO, Path.CURVE3, Path.CURVE3,
         Path.LINETO, Path.CURVE3, Path.CURVE3,
         Path.LINETO, Path.CURVE3, Path.CURVE3,
         Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(v, c), linewidth=LW,
                           edgecolor="black", facecolor="white"))


def txt(x, y, s, size, weight="normal", style="normal"):
    ax.text(x, y, s, ha="center", va="baseline", fontsize=size,
            fontweight=weight, fontstyle=style, color="black")


def line(x0, y0, x1, y1):
    ax.plot([x0, x1], [y0, y1], color="black", lw=LW,
            solid_capstyle="projecting", zorder=1)


def arrow(x, y_from, y_to):
    """Vertical connector with a solid head whose tip sits at y_to."""
    line(x, y_from, x, y_to)
    ax.add_patch(Polygon([(x, y_to),
                          (x - HEAD_W / 2, y_to + HEAD_L),
                          (x + HEAD_W / 2, y_to + HEAD_L)],
                         closed=True, facecolor="black",
                         edgecolor="black", linewidth=LW, zorder=2))


# ---------------------------------------------------------------------------
# Title
# ---------------------------------------------------------------------------
txt(CX, 659.6208, "Supplemental Figure 1. Enrollment Flow Chart", 15,
    weight="bold")
txt(CX, 639.2643, f"Enrollment Period: {DATE_RANGE}", 10.5, style="italic")

# --- Screened --------------------------------------------------------------
box(204.0408, 407.9592, 566.5572, 626.8860)
txt(CX, 603.5108, "Subjects Screened", 12, weight="bold")
txt(CX, 583.6277, f"n = {N_SCREENED}", 11)

line(CX, 566.5572, CX, 546.4476)
line(LX, 546.4476, RX_COL, 546.4476)
arrow(LX, 546.4476, 527.9032)
arrow(RX_COL, 546.4476, 527.9032)

# --- Initial enrollment ----------------------------------------------------
box(18.1152, 222.0336, 466.0092, 526.3380)
txt(LX, 502.9159, "PARDS Group", 12, weight="bold")
txt(LX, 483.0797, f"{N_PARDS_ENROLLED} subjects enrolled", 11)

box(389.9664, 593.8848, 466.0092, 526.3380)
txt(RX_COL, 502.9159, "Control Group", 12, weight="bold")
txt(RX_COL, 483.0797, f"{N_CTRL_ENROLLED} subjects enrolled", 11)

line(LX, 466.0092, LX, 442.5480)
line(RX_COL, 466.0092, RX_COL, 442.5480)
line(LX, 442.5480, RX_COL, 442.5480)
arrow(CX, 442.5480, 424.0036)

# --- Interim analysis ------------------------------------------------------
box(114.0768, 497.9232, 348.7032, 422.4384)
txt(CX, 402.4110, "Interim Methylomic Analysis (Enzymatic MethylSeq)",
    11.5, weight="bold")
txt(CX, 382.5277,
    f"{N_PARDS_INTERIM} PARDS and {N_CTRL_INTERIM} Control subjects had "
    "sufficient DNA quantity", 10.5)
txt(CX, 362.4197, "and quality for methylomic analysis at this time point",
    10.5)

line(CX, 348.7032, CX, 298.4292)
line(LX, 298.4292, RX_COL, 298.4292)
txt(CX, 278.3564, "Continued Enrollment", 12, weight="bold")
arrow(LX, 298.4292, 266.4784)
arrow(RX_COL, 298.4292, 266.4784)

# --- Continued enrollment --------------------------------------------------
box(18.1152, 222.0336, 191.1780, 264.9132)
txt(LX, 244.8397, "PARDS Group", 12, weight="bold")
txt(LX, 225.0074, f"{N_PARDS_NEW} new PARDS subjects", 10.5)
txt(LX, 204.8994, f"Cohort total: {N_PARDS_TOTAL} PARDS subjects", 10.5)

box(389.9664, 593.8848, 191.1780, 264.9132)
txt(RX_COL, 244.8397, "Control Group", 12, weight="bold")
txt(RX_COL, 225.0074, f"{N_CTRL_NEW} new Control subjects", 10.5)
txt(RX_COL, 204.8994, f"Cohort total: {N_CTRL_TOTAL} Control subjects", 10.5)

line(LX, 191.1780, LX, 150.9588)
line(RX_COL, 191.1780, RX_COL, 150.9588)
line(LX, 150.9588, RX_COL, 150.9588)
arrow(CX, 150.9588, 132.4144)

# --- Second analysis -------------------------------------------------------
box(114.0768, 497.9232, 43.7076, 130.8492)
txt(CX, 114.1817, "Second Methylomic Analysis", 12, weight="bold")
txt(CX, 94.0737, "(Infinium MethylationEPIC 2.0)", 12, weight="bold")
txt(CX, 74.1817,
    f"{N_PARDS_SPEC} PARDS specimens ({N_PARDS_RESID} residual "
    f"+ {N_PARDS_FRESH} new)", 10.5)
txt(CX, 54.0737,
    f"{N_CTRL_SPEC} Control specimens ({N_CTRL_RESID} residual "
    f"+ {N_CTRL_PREV} previously enrolled)", 10.5)

# ---------------------------------------------------------------------------
for ext, kw in (("pdf", {}), ("png", {"dpi": 600})):
    path = f"{OUTDIR}/Supplemental Figure 1.{ext}"
    fig.savefig(path, format=ext, **kw)
    print("wrote", path)

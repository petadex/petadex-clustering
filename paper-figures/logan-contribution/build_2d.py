"""Figure 2D: what Logan contributes at the superfamily (30pid) level.

Two 100% bars on one scale: where the Logan proteins go vs where the superfamilies come from.
Most Logan proteins join a superfamily PAZy or nr already has, yet most superfamilies are Logan-only.

Styled to sit next to 2b/2c: 89 mm wide, figure-2 palette (nr plum, Logan green, brown ink), Ubuntu Sans, editable text.
Numbers are from petadex_clustering_experiment.ipynb (figure 1 table; source data/diagnostics/lineage_by_class.tsv).
Writes 2d.logan-contribution.svg (white), .transparent.svg and .png next to this script.
"""
import os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = os.path.dirname(os.path.abspath(__file__))
STEM = os.path.join(OUT, "2d.logan-contribution")

# --- numbers ----------------------------------------------------------------------------------------------------------
SF_KNOWN, SF_LG_MULTI, SF_LG_SINGLE = 8_067, 38_878, 18_824          # superfamilies (30pid clusters)
SF_TOTAL = SF_KNOWN + SF_LG_MULTI + SF_LG_SINGLE                     # 65,769 (2c: n = 65.8k)
LG_TOTAL, LG_IN_KNOWN, LG_IN_NEW = 291_922_243, 263_485_722, 28_436_521   # Logan proteins (complete catalytic domains)
assert LG_IN_KNOWN + LG_IN_NEW == LG_TOTAL

# --- figure-2 style -----------------------------------------------------------------------------------------------------
INK, INK2 = "#5f3928", "#80543f"                      # brown text / secondary (2b "black substitute")
NR, LOGAN, LOGAN_LT = "#9d5e91", "#679a79", "#c1d6ca"  # 2b nr plum, 2b Logan green, 2c light green fill
FONT = "Ubuntu Sans"
plt.rcParams.update({"font.family": ["Ubuntu Sans", "Helvetica", "Arial", "DejaVu Sans"], "font.size": 6.5,
                     "text.color": INK, "svg.fonttype": "none", "figure.facecolor": "white", "savefig.facecolor": "white"})
MM = 1 / 25.4
W, H = 89, 47                                          # mm
fig = plt.figure(figsize=(W * MM, H * MM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")   # draw in mm, y down

X0, X1 = 22, 84                                       # bar span (mm)
BH = 7.5                                               # bar height (mm)
def xs(frac): return X0 + frac * (X1 - X0)

def bar(y, parts):
    """parts: (fraction, colour, label, label colour); labels drawn centred inside when given."""
    acc = 0
    for f, c, lab, lc in parts:
        ax.add_patch(Rectangle((xs(acc), y), xs(acc + f) - xs(acc), BH, fc=c, ec="white", lw=0.6))
        if lab: ax.text((xs(acc) + xs(acc + f)) / 2, y + BH / 2, lab, ha="center", va="center", color=lc,
                        fontsize=6, fontweight="bold", linespacing=1.1)
        acc += f

def row_label(y, title, n):
    ax.text(X0 - 2, y + BH / 2 - 1.3, title, ha="right", va="center", fontsize=7, fontweight="bold")
    ax.text(X0 - 2, y + BH / 2 + 1.7, n, ha="right", va="center", fontsize=6, color=INK2)

# panel letter (same size as 2c's "C": 20 pt)
ax.text(1, 1, "D", ha="left", va="top", fontsize=20, color=INK)

# key (colour = what kind of superfamily); spacing from measured text widths
kx, ky = X0, 4.2
rend = fig.canvas.get_renderer(); inv = ax.transData.inverted()
for c, lab in [(NR, "PAZy or nr superfamily"), (LOGAN, "Logan-only superfamily")]:
    ax.add_patch(Rectangle((kx, ky - 1.3), 2.6, 2.6, fc=c, ec="none"))
    t = ax.text(kx + 3.6, ky, lab, va="center", fontsize=6)
    bb = t.get_window_extent(rend).transformed(inv)
    kx = bb.x1 + 4

# row 1: Logan proteins
y1 = 13
f_in = LG_IN_KNOWN / LG_TOTAL; f_new = LG_IN_NEW / LG_TOTAL
row_label(y1, "Logan proteins", f"n = {LG_TOTAL / 1e6:.0f}M")
bar(y1, [(f_in, NR, f"{LG_IN_KNOWN / 1e6:.1f}M · {f_in:.0%}\njoin a PAZy or nr superfamily", "white"),
         (f_new, LOGAN, f"{f_new:.0%}", "white")])
ax.text(X1, y1 - 1.0, f"{LG_IN_NEW / 1e6:.1f}M form Logan-only superfamilies", ha="right", va="bottom",
        fontsize=6, color=LOGAN, fontweight="bold")

# row 2: superfamilies
y2 = 27
fk, fm, fs = SF_KNOWN / SF_TOTAL, SF_LG_MULTI / SF_TOTAL, SF_LG_SINGLE / SF_TOTAL
row_label(y2, "Superfamilies", f"n = {SF_TOTAL / 1e3:.1f}k")
bar(y2, [(fk, NR, f"{SF_KNOWN / 1e3:.1f}k", "white"),
         (fm, LOGAN, f"{SF_LG_MULTI / 1e3:.1f}k", "white"),
         (fs, LOGAN_LT, f"{SF_LG_SINGLE / 1e3:.1f}k\none family", INK)])

# bracket under the Logan-only superfamilies: what Logan adds
yb = y2 + BH + 1.6
bx0, bx1 = xs(fk) + 0.3, X1 - 0.3
ax.plot([bx0, bx0, bx1, bx1], [yb - 0.9, yb, yb, yb - 0.9], color=LOGAN, lw=0.7, solid_capstyle="butt")
ax.plot([(bx0 + bx1) / 2] * 2, [yb, yb + 0.9], color=LOGAN, lw=0.7)
ax.text((bx0 + bx1) / 2, yb + 1.6, f"+{(SF_LG_MULTI + SF_LG_SINGLE) / 1e3:.1f}k from Logan  ({SF_TOTAL / SF_KNOWN:.1f}× PAZy + nr alone)",
        ha="center", va="top", fontsize=6.5, color=LOGAN, fontweight="bold")

# shared percentage scale
ya = H - 3.2
ax.plot([X0, X1], [ya - 1.6, ya - 1.6], color=INK2, lw=0.5)
for p in (0, 25, 50, 75, 100):
    ax.plot([xs(p / 100)] * 2, [ya - 1.6, ya - 0.9], color=INK2, lw=0.5)
    ax.text(xs(p / 100), ya - 0.5, f"{p}%", ha="center", va="top", fontsize=5.5, color=INK2)

# --- write ------------------------------------------------------------------------------------------------------------
fig.savefig(STEM + ".png", dpi=600)
fig.savefig(STEM + ".svg")
svg = open(STEM + ".svg").read()
# name Ubuntu Sans first whatever matplotlib resolved, so Inkscape/print use it when installed
svg = re.sub(r"font-family:[^;\"]+", f"font-family: '{FONT}', Helvetica, Arial, sans-serif", svg)
svg = re.sub(r"font: (\d+(?:\.\d+)?)px '[^']+'", rf"font-size: \1px; font-family: '{FONT}', Helvetica, Arial, sans-serif", svg)
open(STEM + ".svg", "w").write(svg)
# transparent copy: drop the figure background patch
open(STEM + ".transparent.svg", "w").write(re.sub(r'(<g id="patch_1">\s*<path[^>]*?)style="fill: #ffffff"', r'\1style="fill: none"', svg, count=1))
print(dict(sf_total=SF_TOTAL, pct_logan_in_known=round(f_in * 100, 1), sf_logan_only_pct=round((fm + fs) * 100, 1), fold=round(SF_TOTAL / SF_KNOWN, 1)))

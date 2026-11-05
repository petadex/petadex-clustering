"""Method figure: the clustering procedure for one CATH domain at one identity level (white background, PNG + SVG).
Counts are the 90% level summed over all 15 domains: ../data/method_flow_counts.tsv (count_method_flow.sh),
../data/cluster_stats.tsv (centroids by source) and ../data/counts_L90.tsv (inputs).
Step numbers: 1-8 = upstream unid/ steps 1-8; 9 = upstream 14 (merge); 10 = upstream 15 (extract centroids).
Upstream steps 9-13 are not run (their output is unused) and are not shown."""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data"); FIG = os.path.join(ROOT, "figures")
INK = "#1F2933"; GREY = "#5F6B76"
C_PL, C_NR, C_LG, C_OUT = "#E69F00", "#0072B2", "#D55E00", "#009E73"     # source colours as in fig 2
plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "path", "figure.facecolor": "white", "savefig.facecolor": "white"})

S = pd.read_csv(os.path.join(DATA, "cluster_stats.tsv"), sep="\t"); S = S[S.level == "L90"]
F = pd.read_csv(os.path.join(DATA, "method_flow_counts.tsv"), sep="\t").set_index("level").loc["L90"]
inp = pd.read_csv(os.path.join(DATA, "counts_L90.tsv"), sep="\t", header=None)
n_pl, n_nr, n_lg = 486, 6_764_617, 291_922_243            # run2/input90/total.txt
assert inp[1].sum() == n_pl + n_nr + n_lg
c_pl, c_nr, c_lg = int(S.PL.sum()), int(S.NR.sum()), int(S.LG.sum())
s2, s4, s6 = int(F.nr_to_pazy), int(F.logan_step4), int(F.logan_step6)
nr_left, lg_left = n_nr - s2, n_lg - s4 - s6
f = lambda x: f"{int(x):,}"

fig, ax = plt.subplots(figsize=(11.5, 8.6)); fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.01)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
XP, XN, XL, W, H = 17, 50, 83, 29, 7.6

def box(x, y, title, detail, color, num=None, w=W, h=H):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=1.1",
                                fc=matplotlib.colors.to_rgba(color, 0.12), ec=color, lw=1.6, zorder=2))
    dx = 2.0 if num else 0
    ax.text(x + dx, y + h * 0.17, title, ha="center", va="center", fontsize=10, fontweight="bold", color=INK, zorder=4)
    ax.text(x + dx, y - h * 0.22, detail, ha="center", va="center", fontsize=9.3, color=INK, zorder=4)
    if num:
        ax.add_patch(Ellipse((x - w / 2 + 2.6, y), 3.5, 4.7, fc=color, ec="white", lw=1.2, zorder=5))   # x/y units differ: draws a circle
        ax.text(x - w / 2 + 2.6, y, num, ha="center", va="center", fontsize=8.6 if len(num) < 3 else 7.2, color="white", fontweight="bold", zorder=6)

def arrow(x0, y0, x1, y1, color=INK, ls="-", conn="arc3"):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=1,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.5, ls=ls, shrinkA=0, shrinkB=0, connectionstyle=conn))

for x, c, t in [(XP, C_PL, "PAZy  ·  priority 1"), (XN, C_NR, "nr  ·  priority 2"), (XL, C_LG, "Logan  ·  priority 3")]:
    ax.text(x, 96, t, ha="center", va="center", fontsize=11, fontweight="bold", color=c)
    ax.plot([x - W / 2, x + W / 2], [93.4, 93.4], color=c, lw=2)

T = 88; box(XP, T, "PAZy cores", f(n_pl), C_PL); box(XN, T, "nr cores", f(n_nr), C_NR); box(XL, T, "Logan cores", f(n_lg), C_LG)
box(XP, 75, "Cluster PAZy", f"→ {f(c_pl)} centroids", C_PL, "1")
box(XN, 75, "Map nr onto PAZy", f"{f(s2)} assigned", C_NR, "2")
box(XN, 62, "Cluster remaining nr", f"{f(nr_left)} → {f(c_nr)} centroids", C_NR, "3")
box(33.5, 49, "PAZy + nr centroids", f(c_pl + c_nr), GREY, w=30)
box(XL, 49, "Map Logan onto centroids", f"{f(s4)} assigned", C_LG, "4")
box(XL, 37.5, "Re-map leftovers, more sensitively", f"{f(s6)} assigned", C_LG, "5–6", w=W + 2)
box(XL, 26, "Collect unassigned Logan", f"{f(lg_left)} left", C_LG, "7")
box(XL, 14.5, "Cluster remaining Logan", f"→ {f(c_lg)} centroids", C_LG, "8")
box(35, 3.8, "Merge all assignments", f"{f(n_pl + n_nr + n_lg)} rows, one per core", C_OUT, "9", w=62)
box(XL, 3.8, "Centroids feed next level", "90% → 60% → 30%", C_OUT, "10")

for x in (XP, XN, XL): arrow(x, T - H / 2, x, 75 + H / 2 if x != XL else 49 + H / 2, color=C_LG if x == XL else INK, ls=(0, (4, 3)) if x == XL else "-")
arrow(XP + W / 2, 75, XN - W / 2, 75, color=C_PL)                                   # PAZy centroids are the targets of step 2
arrow(XN, 75 - H / 2, XN, 62 + H / 2)
arrow(XP + 4, 75 - H / 2, 27, 49 + H / 2, color=C_PL); arrow(XN, 62 - H / 2, 41, 49 + H / 2, color=C_NR)
arrow(48.5, 49, XL - W / 2, 49, color=GREY)
for ya, yb in [(49, 37.5), (37.5, 26), (26, 14.5)]: arrow(XL, ya - H / 2, XL, yb + H / 2, color=C_LG)
arrow(XP - 7, 75 - H / 2, XP - 7, 3.8 + H / 2, color=C_PL); arrow(24, 49 - H / 2, 24, 3.8 + H / 2, color=C_NR)
arrow(XL - W / 2, 13, 60, 3.8 + H / 2, color=C_LG, conn="arc3,rad=0.25")
arrow(66, 3.8, XL - W / 2, 3.8, color=C_OUT)

fig.text(0.02, 0.985, "Clustering one CATH domain at one identity level (counts: 90%, all domains)", fontsize=13, fontweight="bold", va="top", color=INK)
for ext in ("png", "svg"): fig.savefig(os.path.join(FIG, f"fig0_method_steps.{ext}"), dpi=300)
print("wrote fig0_method_steps")

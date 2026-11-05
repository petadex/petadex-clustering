"""Figure 1 (main): what Logan adds at the superfamily (30%) level. White background, PNG + SVG.
A: superfamilies from PAZy + nr alone vs with Logan. Because PAZy and nr are clustered (steps 1-3, every level) before Logan is
   involved, PAZy + nr alone give exactly the PAZy- and nr-centred superfamilies; every Logan-centred superfamily is Logan-only.
B: where the Logan cores end up: inside an existing PAZy- or nr-centred superfamily, or in a Logan-only superfamily.
Reads ../data/diagnostics/lineage_by_class.tsv (diag_lineage.sh)."""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data"); FIG = os.path.join(ROOT, "figures")
INK, SUB = "#1F2933", "#5F6B76"; C_KNOWN, C_LG, C_LG1 = "#0072B2", "#D55E00", "#F2B48C"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#9AA5B1", "axes.linewidth": 0.8, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "figure.facecolor": "white", "savefig.facecolor": "white",
                     "svg.fonttype": "path"})

B = pd.read_csv(os.path.join(DATA, "diagnostics", "lineage_by_class.tsv"), sep="\t"); b = B[B.level == "L30"]
known = int(b[b.source.isin(["PL", "NR"])].clusters.sum())
lg_multi = int(b[(b.source == "LG") & (b.singleton == 0)].clusters.sum()); lg_single = int(b[(b.source == "LG") & (b.singleton == 1)].clusters.sum())
lg_only_cores = int(b[b.source == "LG"].cores.sum())
LOGAN_CORES = 291_922_243                                   # run2/input90/total.txt
lg_joined = LOGAN_CORES - lg_only_cores
total = known + lg_multi + lg_single
pd.DataFrame([{"superfamilies_pazy_nr": known, "logan_only_multi": lg_multi, "logan_only_single": lg_single, "superfamilies_total": total,
               "logan_cores": LOGAN_CORES, "logan_cores_in_existing": lg_joined, "logan_cores_in_logan_only": lg_only_cores}]
             ).to_csv(os.path.join(DATA, "figure1_logan_adds.tsv"), sep="\t", index=False)

fig = plt.figure(figsize=(10.4, 5.4))
axA = fig.add_axes([0.08, 0.14, 0.36, 0.64]); axB = fig.add_axes([0.53, 0.40, 0.44, 0.20])

# --- A ---------------------------------------------------------------------------------------------------------------
axA.bar(0, known, 0.55, color=C_KNOWN)
axA.bar(1, known, 0.55, color=C_KNOWN, label="PAZy or nr")
axA.bar(1, lg_multi, 0.55, bottom=known, color=C_LG, label="Logan-only, >1 family")
axA.bar(1, lg_single, 0.55, bottom=known + lg_multi, color=C_LG1, label="Logan-only, 1 family")
axA.text(0, known + 1500, f"{known:,}", ha="center", va="bottom", fontsize=10, fontweight="bold")
axA.text(1, total + 1500, f"{total:,}", ha="center", va="bottom", fontsize=10, fontweight="bold")
for y0, h, n in [(known, lg_multi, lg_multi), (known + lg_multi, lg_single, lg_single)]:
    axA.text(1, y0 + h / 2, f"{n:,}", ha="center", va="center", color="white" if n == lg_multi else INK, fontsize=9, fontweight="bold")
axA.text(0.5, total * 0.62, f"+{lg_multi + lg_single:,}\n({total / known:.1f}×)", ha="center", va="center", fontsize=10, color=C_LG, fontweight="bold")
axA.set_xticks([0, 1]); axA.set_xticklabels(["PAZy + nr", "PAZy + nr\n+ Logan"]); axA.set_xlim(-0.5, 1.5)
axA.set_ylim(0, total * 1.12); axA.set_ylabel("Superfamilies (30% identity)")
axA.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v/1e3:g}k" if v else "0"))
axA.legend(frameon=False, loc="upper left", fontsize=8.6)
axA.set_title("A   Logan adds most of the superfamilies", loc="left", fontsize=11, fontweight="bold", pad=12)

# --- B ---------------------------------------------------------------------------------------------------------------
pj, pn = lg_joined / LOGAN_CORES * 100, lg_only_cores / LOGAN_CORES * 100
axB.barh(0, pj, 0.6, color=C_KNOWN); axB.barh(0, pn, 0.6, left=pj, color=C_LG)
axB.text(pj / 2, 0, f"join a PAZy or nr superfamily\n{lg_joined / 1e6:.1f}M cores · {pj:.1f}%", ha="center", va="center", color="white", fontsize=9, fontweight="bold")
axB.annotate(f"form Logan-only superfamilies\n{lg_only_cores / 1e6:.1f}M cores · {pn:.1f}%", xy=(pj + pn / 2, -0.32), xytext=(100, -1.25),
             ha="right", va="top", fontsize=9, color=C_LG, fontweight="bold", arrowprops=dict(arrowstyle="-", color=C_LG, lw=1))
axB.set_xlim(0, 100); axB.set_ylim(-0.5, 0.5); axB.set_yticks([]); axB.set_xticks([0, 25, 50, 75, 100])
axB.set_xticklabels([f"{t}%" for t in (0, 25, 50, 75, 100)], fontsize=8.5)
for sp in ("top", "right", "left"): axB.spines[sp].set_visible(False)
axB.set_title(f"B   Where the {LOGAN_CORES / 1e6:.1f}M Logan cores go", loc="left", fontsize=11, fontweight="bold", pad=12)
for sp in ("top", "right"): axA.spines[sp].set_visible(False)

fig.text(0.02, 0.975, "What Logan adds to PETadex superfamilies", fontsize=13, fontweight="bold", va="top")
fig.text(0.02, 0.935, "Superfamily = cluster at 30% identity. Logan-only: no PAZy or nr member under the clustering rule (≥30% identity, ≥90% coverage).", fontsize=8.8, color=SUB, va="top")
for ext in ("png", "svg"): fig.savefig(os.path.join(FIG, f"fig1_logan_adds.{ext}"), dpi=300)
print(dict(known=known, lg_multi=lg_multi, lg_single=lg_single, total=total, lg_joined=lg_joined, lg_only_cores=lg_only_cores, pj=round(pj, 1)))

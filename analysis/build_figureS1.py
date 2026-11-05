"""Supplementary figure S1: share of clusters vs share of sequences, by centroid source, at each identity level (white background, PNG + SVG).
Reads ../data/diagnostics/lineage_by_class.tsv (scripts/diagnose/diag_lineage.sh): for every cluster, its centroid source and the
number of original cores beneath it (following centroids down the hierarchy). Because sources are clustered in priority order
PAZy -> nr -> Logan, a Logan-centred cluster contains only Logan sequences at every level below it.
Writes ../data/figureS1_clusters_vs_cores.tsv."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data"); FIG = os.path.join(ROOT, "figures")
INK, SUB = "#1F2933", "#5F6B76"
SRC = [("LG", "Logan only", "#D55E00"), ("NR", "nr-centred", "#0072B2"), ("PL", "PAZy-centred", "#E69F00")]
LEVELS = [("L90", "90% identity"), ("L60", "60% identity"), ("L30", "30% identity")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#9AA5B1", "axes.linewidth": 0.8, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "figure.facecolor": "white", "savefig.facecolor": "white",
                     "svg.fonttype": "path"})

B = pd.read_csv(os.path.join(DATA, "diagnostics", "lineage_by_class.tsv"), sep="\t")
G = B.groupby(["level", "source"])[["clusters", "cores"]].sum()
rows = []
for L, _ in LEVELS:
    g = G.loc[L]
    for s, *_ in SRC:
        rows.append({"level": L, "source": s, "clusters": int(g.loc[s, "clusters"]), "cores": int(g.loc[s, "cores"]),
                     "pct_clusters": g.loc[s, "clusters"] / g.clusters.sum() * 100, "pct_cores": g.loc[s, "cores"] / g.cores.sum() * 100})
T = pd.DataFrame(rows); T.to_csv(os.path.join(DATA, "figureS1_clusters_vs_cores.tsv"), sep="\t", index=False, float_format="%.3f")
n_cores = int(G.loc["L90"].cores.sum())

def pct(v): return f"{v:.1f}%" if v >= 1 else (f"{v:.2f}%" if v >= 0.01 else f"{v:.3f}%")
def cnt(v): return f"{v/1e6:.1f}M" if v >= 1e6 else (f"{v/1e3:.1f}k" if v >= 1e3 else f"{v:,}")

fig, ax = plt.subplots(figsize=(10.2, 6.0)); fig.subplots_adjust(top=0.80, bottom=0.15, left=0.08, right=0.74)
W = 0.36; GAP = 1.15
for j, (L, lab) in enumerate(LEVELS):
    t = T[T.level == L].set_index("source")
    for k, (col, what) in enumerate([("pct_clusters", "clusters"), ("pct_cores", "sequences")]):
        x = j * GAP + (k - 0.5) * (W + 0.04); bottom = 0
        for s, name, c in SRC:
            h = t.loc[s, col]
            ax.bar(x, h, W, bottom=bottom, color=c, edgecolor="white", linewidth=0.8)
            raw = t.loc[s, "clusters" if what == "clusters" else "cores"]
            if h >= 10:
                ax.text(x, bottom + h / 2, f"{pct(h)}\n{cnt(raw)}", ha="center", va="center", color="white", fontsize=8.3, fontweight="bold")
            elif h >= 5:
                ax.text(x, bottom + h / 2, pct(h), ha="center", va="center", color="white", fontsize=8.3, fontweight="bold")
            bottom += h
        # PAZy share is invisible on most bars: report it above the bar
        p = t.loc["PL", col]
        if p < 6:
            ax.text(x, 101.5, f"PAZy\n{pct(p)}", ha="center", va="bottom", fontsize=7.4, color="#B45309")
        ax.text(x, -3.5, what, ha="center", va="top", fontsize=8.8, color=SUB)
    tot_c = int(T[T.level == L].clusters.sum())
    ax.text(j * GAP, -11, f"{lab}", ha="center", va="top", fontsize=10, fontweight="bold")
    ax.text(j * GAP, -16, f"{tot_c:,} clusters", ha="center", va="top", fontsize=8.6, color=SUB)
ax.set_xlim(-0.55, 2 * GAP + 0.55); ax.set_ylim(0, 112); ax.set_xticks([]); ax.set_yticks(range(0, 101, 20))
ax.set_ylabel("Share (%)")
for sp in ("top", "right", "bottom"): ax.spines[sp].set_visible(False)
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in SRC][::-1]
ax.legend(handles, [n for _, n, _ in SRC][::-1], frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.74), title="Cluster led by", title_fontsize=9.5)
t30 = T[T.level == "L30"].set_index("source")
note = (f"30% identity (superfamilies):\n\n{t30.loc['LG','clusters']:,} superfamilies\n({pct(t30.loc['LG','pct_clusters'])}) have no nr or\n"
        f"PAZy sequence; they hold\n{cnt(t30.loc['LG','cores'])} cores ({pct(t30.loc['LG','pct_cores'])}).\n\n"
        f"{t30.loc['PL','clusters']} PAZy-led superfamilies\nhold {cnt(t30.loc['PL','cores'])} cores ({pct(t30.loc['PL','pct_cores'])}).")
ax.text(1.03, 0.36, note, transform=ax.transAxes, fontsize=8.6, color=INK, va="top", ha="left")
fig.text(0.02, 0.985, "Most superfamilies are Logan-only; most sequences sit in superfamilies led by nr or PAZy", fontsize=13, fontweight="bold", va="top")
fig.text(0.02, 0.945, f"For each level: share of clusters, and share of all {n_cores:,} catalytic cores beneath them, by the source of the cluster's centroid",
         fontsize=9, color=SUB, va="top")
for ext in ("png", "svg"): fig.savefig(os.path.join(FIG, f"figS1_clusters_vs_sequences.{ext}"), dpi=300)
print(T.round(2).to_string(index=False))

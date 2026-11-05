"""Cluster size distribution per level (fig4), white background, PNG + SVG.
Reads ../data/diagnostics/cluster_size_hist.tsv (scripts/diagnose/diag_size_hist.sh): exact histogram of direct members per cluster
(cores at 90%, 90% centroids at 60%, 60% centroids at 30%). Writes the binned table to ../data/diagnostics/cluster_size_binned.tsv."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
DIAG = os.path.join(ROOT, "data", "diagnostics"); FIG = os.path.join(ROOT, "figures")
INK, SUB = "#1F2933", "#5F6B76"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#9AA5B1", "axes.linewidth": 0.8, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "figure.facecolor": "white", "savefig.facecolor": "white",
                     "svg.fonttype": "path"})
BINS = [(1, 1, "1"), (2, 9, "2–9"), (10, 99, "10–99"), (100, 999, "100–999"), (1000, 9999, "1,000–\n9,999"), (10000, 10**12, "≥10,000")]
LEVELS = [("L90", "90% identity", "Cluster size (cores)", "#009E73"),
          ("L60", "60% identity", "Cluster size (90% centroids)", "#7B5EA7"),
          ("L30", "30% identity", "Cluster size (60% centroids)", "#E69F00")]

H = pd.read_csv(os.path.join(DIAG, "cluster_size_hist.tsv"), sep="\t")
rows = []
for L, *_ in LEVELS:
    h = H[H.level == L]; tot = h.clusters.sum()
    for lo, hi, lab in BINS:
        x = h[(h["size"] >= lo) & (h["size"] <= hi)]
        rows.append({"level": L, "size": lab.replace("\n", ""), "clusters": int(x.clusters.sum()), "pct_clusters": x.clusters.sum() / tot * 100,
                     "members": int((x["size"] * x.clusters).sum()), "pct_members": (x["size"] * x.clusters).sum() / (h["size"] * h.clusters).sum() * 100})
T = pd.DataFrame(rows); T.to_csv(os.path.join(DIAG, "cluster_size_binned.tsv"), sep="\t", index=False, float_format="%.3f")

def fmt(v, _=None): return f"{v/1e6:g}M" if v >= 1e6 else (f"{v/1e3:g}k" if v >= 1e3 else f"{v:g}")
fig, axs = plt.subplots(3, 1, figsize=(7.6, 9.4)); fig.subplots_adjust(top=0.895, bottom=0.075, left=0.13, right=0.97, hspace=0.62)
x = np.arange(len(BINS))
for ax, (L, title, xlab, c) in zip(axs, LEVELS):
    t = T[T.level == L]
    ax.bar(x, t.clusters, width=0.72, color=c, alpha=0.85, edgecolor=c)
    top = t.clusters.max()
    for xi, (n, p) in enumerate(zip(t.clusters, t.pct_clusters)):
        ax.text(xi, n + top * 0.03, f"{n:,}\n{p:.1f}%" if p >= 0.1 else f"{n:,}\n<0.1%", ha="center", va="bottom", fontsize=7.8, color=INK)
    ax.set_ylim(0, top * 1.32); ax.yaxis.set_major_formatter(FuncFormatter(fmt))
    ax.set_xticks(x); ax.set_xticklabels([b[2] for b in BINS], fontsize=9); ax.set_xlabel(xlab, fontsize=9.5)
    ax.set_ylabel("Number of clusters"); ax.set_title(f"{title}  ·  {int(t.clusters.sum()):,} clusters", loc="left", fontsize=10.5, fontweight="bold")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
fig.text(0.02, 0.985, "Cluster size distribution", fontsize=13, fontweight="bold", va="top")
fig.text(0.02, 0.962, "Size = direct members of a cluster; labels give the count and the share of that level's clusters", fontsize=9, color=SUB, va="top")
for ext in ("png", "svg"): fig.savefig(os.path.join(FIG, f"fig4_cluster_size_distribution.{ext}"), dpi=300)
print(T.round(2).to_string(index=False))

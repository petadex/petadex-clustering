"""Fragmentation figure (fig7), white background, PNG + SVG.
Left: length distribution of all input cores (90% level input = every core), nr vs Logan, from
      ../data/fragmentation/core_lengths_L90.tsv (scripts/diagnose/diag_fragmentation.sh).
Right: mean length of singleton vs non-singleton centroids, Logan and nr, at each level, from
      ../data/diagnostics/clusters_by_class.tsv (scripts/diagnose/diag_clusters.sh).
Also writes summary tables into ../data/fragmentation/."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
D = os.path.join(ROOT, "data", "fragmentation"); DIAG = os.path.join(ROOT, "data", "diagnostics"); FIG = os.path.join(ROOT, "figures")
C_NR, C_LG = "#0072B2", "#D55E00"; INK, SUB = "#1F2933", "#5F6B76"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#9AA5B1", "axes.linewidth": 0.8, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "figure.facecolor": "white", "savefig.facecolor": "white",
                     "svg.fonttype": "path"})
LV = {"L90": "90%", "L60": "60%", "L30": "30%"}

# --- input core length distribution -------------------------------------------------------------------------------
Lh = pd.read_csv(os.path.join(D, "core_lengths_L90.tsv"), sep="\t")
Lh = Lh[Lh.source.isin(["NR", "LG"])]
def wq(g, q):
    g = g.sort_values("core_length"); c = g.n_cores.cumsum(); return int(g.core_length[c >= g.n_cores.sum() * q].iloc[0])
stats = pd.DataFrame({s: {"cores": int(g.n_cores.sum()), "mean length (aa)": (g.core_length * g.n_cores).sum() / g.n_cores.sum(),
                          "25th percentile": wq(g, 0.25), "median": wq(g, 0.5), "75th percentile": wq(g, 0.75),
                          "% shorter than 150 aa": g[g.core_length < 150].n_cores.sum() / g.n_cores.sum() * 100}
                      for s, g in Lh.groupby("source")}).T.rename(index={"NR": "nr", "LG": "Logan"})
stats.to_csv(os.path.join(D, "core_length_stats_L90.tsv"), sep="\t", float_format="%.1f")
BW = 10; edges = np.arange(0, 801, BW)
hist = {s: np.histogram(g.core_length.clip(upper=799), bins=edges, weights=g.n_cores)[0] for s, g in Lh.groupby("source")}

# --- mean centroid length, singleton vs non-singleton ---------------------------------------------------------------
A = pd.read_csv(os.path.join(DIAG, "clusters_by_class.tsv"), sep="\t")
M = A[A.source.isin(["NR", "LG"])].groupby(["level", "source", "singleton"])[["clusters", "sum_len"]].sum()
M["mean_len"] = M.sum_len / M.clusters
M = M.reset_index(); M.to_csv(os.path.join(D, "mean_centroid_length_by_singleton.tsv"), sep="\t", index=False, float_format="%.1f")

fig, axs = plt.subplots(1, 2, figsize=(10.6, 4.6), gridspec_kw={"width_ratios": [1.15, 1]})
fig.subplots_adjust(top=0.80, bottom=0.14, left=0.07, right=0.98, wspace=0.24)
ax = axs[0]
for s, c, nm in [("NR", C_NR, "nr"), ("LG", C_LG, "Logan")]:
    pct = hist[s] / hist[s].sum() * 100
    ax.fill_between(edges[:-1] + BW / 2, pct, step="mid", color=c, alpha=0.18, lw=0)
    ax.step(edges[:-1] + BW / 2, pct, where="mid", color=c, lw=1.6, label=f"{nm} (median {stats.loc['nr' if s == 'NR' else 'Logan', 'median']:.0f} aa)")
ax.set_xlim(0, 700); ax.set_ylim(bottom=0); ax.set_xlabel("Core length (aa)"); ax.set_ylabel(f"Share of input cores per {BW} aa (%)")
ax.set_title("Logan cores are shorter than nr cores", loc="left", fontsize=10.5, fontweight="bold"); ax.legend(frameon=False)

ax = axs[1]; x = np.arange(3); w = 0.2
specs = [("LG", 0, C_LG, 1.0, "Logan, non-singleton"), ("LG", 1, C_LG, 0.45, "Logan, singleton"),
         ("NR", 0, C_NR, 1.0, "nr, non-singleton"), ("NR", 1, C_NR, 0.45, "nr, singleton")]
for k, (s, sg, c, a, nm) in enumerate(specs):
    v = [M[(M.level == L) & (M.source == s) & (M.singleton == sg)].mean_len.iloc[0] for L in LV]
    xs = x + (k - 1.5) * w
    ax.bar(xs, v, w, color=c, alpha=a, label=nm, edgecolor="white", linewidth=0.5)
    for xi, vi in zip(xs, v): ax.text(xi, vi + 4, f"{vi:.0f}", ha="center", fontsize=7.6)
ax.set_xticks(x); ax.set_xticklabels([f"{LV[L]} identity" for L in LV]); ax.set_ylim(0, 430); ax.set_ylabel("Mean centroid length (aa)")
ax.set_title("Singletons are shorter than non-singletons", loc="left", fontsize=10.5, fontweight="bold")
ax.legend(frameon=False, ncol=2, fontsize=8.3, loc="upper center", bbox_to_anchor=(0.5, 1.0))
for a_ in axs:
    for sp in ("top", "right"): a_.spines[sp].set_visible(False)
fig.text(0.02, 0.985, "Logan fragmentation and singleton length", fontsize=13, fontweight="bold", va="top")
fig.text(0.02, 0.945, "Left: every input core (90% level input). Right: centroids at each level; singleton = cluster with one member",
         fontsize=9, color=SUB, va="top")
for ext in ("png", "svg"): fig.savefig(os.path.join(FIG, f"fig7_logan_fragmentation.{ext}"), dpi=300)
print(stats.round(1).to_string()); print(M[["level", "source", "singleton", "clusters", "mean_len"]].round(1).to_string(index=False))

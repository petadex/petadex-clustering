"""Diagnostic figures for the singleton question (fig6-fig8), white background, PNG + SVG.
Reads ../data/cluster_stats.tsv, ../data/diagnostics/*.tsv and ../data/v11/source_breakdown_vs_benjamin.tsv
(v1.1 per-source numbers, from Benjamin's old.tsv). Also writes summary tables to ../data/diagnostics/."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data"); DIAG = os.path.join(DATA, "diagnostics"); FIG = os.path.join(ROOT, "figures")
RES = os.path.join(ROOT, "data", "v11")   # v1.1 comparison tables (from the Run 1 results/)
PAL = ["#0072B2", "#009E73", "#E69F00", "#D55E00", "#CC79A7", "#56B4E9", "#7B5EA7"]; GREY = "#8C96A0"; INK, SUB = "#1F2933", "#5F6B76"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#9AA5B1", "axes.linewidth": 0.8,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK,
                     "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "svg.fonttype": "path"})
LV = ["L90", "L60", "L30"]; LVL = {"L90": "90% identity", "L60": "60% identity", "L30": "30% identity"}
def finish(fig, name, title, sub):
    fig.text(0.02, 0.985, title, fontsize=13, fontweight="bold", va="top"); fig.text(0.02, 0.945, sub, fontsize=9, color=SUB, va="top")
    for ext in ("png", "svg"): fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300)
    plt.close(fig)
def clean(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)

A = pd.read_csv(os.path.join(DIAG, "clusters_by_class.tsv"), sep="\t")
B = pd.read_csv(os.path.join(DIAG, "lineage_by_class.tsv"), sep="\t")
P = pd.read_csv(os.path.join(DIAG, "pairwise_30.tsv"), sep="\t")
V = pd.read_csv(os.path.join(RES, "source_breakdown_vs_benjamin.tsv"), sep="\t"); V = V[V.level != "input"].copy()
V["level"] = "L" + V["level"].astype(str)

# ---- fig6: singleton rate by centroid source, v2 vs v1.1 ----------------------------------------------------------
g = A.groupby(["level", "source", "singleton"])["clusters"].sum().unstack(fill_value=0); g["total"] = g.sum(axis=1)
rate_v2 = (g[1] / g["total"] * 100).to_dict()
v1 = V.set_index(["level", "source"]); rate_v11 = (v1["single_member_benjamin"] / v1["centroids_benjamin"] * 100).to_dict()
rows = [{"level": L, "source": s, "v2_clusters": int(g.loc[(L, s), "total"]), "v2_singletons": int(g.loc[(L, s), 1]), "v2_rate": rate_v2[(L, s)],
         "v11_clusters": int(v1.loc[(L, s.lower() if False else {'PL':'pazy','NR':'nr','LG':'logan'}[s]), 'centroids_benjamin']),
         "v11_singletons": int(v1.loc[(L, {'PL':'pazy','NR':'nr','LG':'logan'}[s]), 'single_member_benjamin']),
         "v11_rate": rate_v11[(L, {'PL':'pazy','NR':'nr','LG':'logan'}[s])]} for L in LV for s in ("LG", "NR", "PL")]
pd.DataFrame(rows).to_csv(os.path.join(DIAG, "singleton_rate_by_source.tsv"), sep="\t", index=False, float_format="%.2f")
fig, axs = plt.subplots(1, 2, figsize=(9.4, 4.4), sharey=True); fig.subplots_adjust(top=0.78, bottom=0.14, left=0.08, right=0.98, wspace=0.12)
for ax, (key, lab) in zip(axs, [("v2_rate", "v2"), ("v11_rate", "v1.1")]):
    x = np.arange(3); w = 0.36
    for k, (s, c, nm) in enumerate([("LG", PAL[3], "Logan-centred"), ("NR", PAL[0], "nr-centred")]):
        vals = [[r for r in rows if r["level"] == L and r["source"] == s][0][key] for L in LV]
        ax.bar(x + (k - 0.5) * w, vals, w, color=c, label=nm)
        for xi, v in zip(x + (k - 0.5) * w, vals): ax.text(xi, v + 1, f"{v:.1f}%", ha="center", fontsize=8.5)
    ax.set_xticks(x); ax.set_xticklabels([LVL[L] for L in LV]); ax.set_ylim(0, 68); ax.set_title(lab, fontsize=11, loc="left", fontweight="bold"); clean(ax)
axs[0].set_ylabel("Singleton clusters (% of that source's clusters)"); axs[1].legend(frameon=False, loc="upper right")
finish(fig, "fig6_singleton_rate_by_source", "Singleton rate of Logan- vs nr-centred clusters",
       "Share of clusters with one member, by centroid source; PAZy omitted (360 / 204 / 76 centroids)")

# ---- fig7: singleton share of clusters vs of cores ----------------------------------------------------------------
tot_cores = B[B.level == "L90"]["cores"].sum()
S = pd.read_csv(os.path.join(DATA, "cluster_stats.tsv"), sep="\t")
rows7 = []
for L in LV:
    n = S[S.level == L].clusters.sum(); s_ = S[S.level == L].singletons.sum(); c = B[(B.level == L) & (B.singleton == 1)].cores.sum()
    rows7.append({"level": L, "clusters": n, "singleton_clusters": s_, "pct_clusters": s_ / n * 100, "cores_in_singleton_clusters": c, "pct_cores": c / tot_cores * 100})
pd.DataFrame(rows7).to_csv(os.path.join(DIAG, "singleton_share_clusters_vs_cores.tsv"), sep="\t", index=False, float_format="%.4f")

print("wrote fig6 and diagnostic tables")

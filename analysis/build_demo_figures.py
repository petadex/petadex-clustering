"""Demo figures for the PETadex Run 2 clustering (white background, matplotlib only).
Reads ../data/ (cluster_stats.tsv from extract_stats.sh, counts_L*.tsv, v11_reference.tsv); writes PNG + SVG to ../figures/.
Final data: all 15 domains, all three levels."""
import csv, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE); DATA = os.path.join(ROOT, "data"); FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
LEVELS = ["L90", "L60", "L30"]
LV_LABEL = {"L90": "90% identity", "L60": "60% identity", "L30": "30% identity"}
NOTE = sys.argv[1] if len(sys.argv) > 1 else ""

NAMES = {"abhydrolase": "α/β hydrolase", "betalactamase": "DD-peptidase/β-lactamase", "cytochrome-p450": "Cytochrome P450",
         "amidase-signature": "Amidase signature", "peptidase-s8-s53": "Peptidase S8/S53", "metallo-phosphatase": "Metallo-phosphatases",
         "trypsin-like": "Trypsin-like", "cupredoxins": "Cupredoxins", "d-ala-esterase": "D-Ala-D-Ala carboxypeptidase-like",
         "acid-protease": "Acid protease", "sgnh-hydrolase": "SGNH hydrolase", "oxygenase-like": "Oxygenase-like",
         "gh-deacetylase": "GH deacetylase", "cytochromec-like": "Cytochrome c-like", "peroxidase": "Peroxidase"}
PAL = ["#0072B2", "#009E73", "#E69F00", "#D55E00", "#CC79A7", "#56B4E9", "#7B5EA7"]   # colour-blind-safe (Okabe-Ito based)
GREY = "#8C96A0"
INK, SUB = "#1F2933", "#5F6B76"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#9AA5B1", "axes.linewidth": 0.8,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK,
                     "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "svg.fonttype": "path"})

def read_stats():
    rows = {}
    with open(os.path.join(DATA, "cluster_stats.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows[(r["level"], r["domain"])] = {k: int(v) for k, v in r.items() if k not in ("level", "domain")}
    return rows

def read_inputs():
    inp = {}
    for L in LEVELS:
        with open(os.path.join(DATA, f"counts_{L}.tsv")) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 6: inp[(L, p[0])] = int(p[1])
    return inp

S = read_stats(); INP = read_inputs()
doms = {L: sorted({d for (l, d) in S if l == L}) for L in LEVELS}
def tot(L, key, ds=None): return sum(S[(L, d)][key] for d in (ds or doms[L]))
def millions(x, _=None):
    return f"{x/1e6:g}M" if x >= 1e6 else (f"{x/1e3:g}k" if x >= 1e3 else f"{x:g}")
def compact(x): return f"{x/1e6:.1f}M" if x >= 1e6 else f"{x/1e3:.1f}k"
def finish(fig, name, title, sub):
    fig.text(0.02, 0.985, title, fontsize=13, fontweight="bold", va="top", ha="left")
    fig.text(0.02, 0.945, sub, fontsize=9, color=SUB, va="top", ha="left")
    fig.text(0.98, 0.006, NOTE, fontsize=8, color="#B45309", ha="right", va="bottom", style="italic")
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300)
    plt.close(fig)

# ---- 2. who owns each level ---------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.8, 5.2)); fig.subplots_adjust(top=0.80, bottom=0.17, left=0.11, right=0.97)
src = [("LG", "Logan", PAL[3]), ("NR", "nr", PAL[0]), ("PL", "PAZy", PAL[2])]
for i, L in enumerate(LEVELS):
    t = tot(L, "clusters"); bottom = 0
    for k, lab, c in src:
        h = tot(L, k) / t * 100
        ax.bar(i, h, bottom=bottom, width=0.55, color=c, label=lab if i == 0 else None, edgecolor="white", linewidth=0.8)
        if h > 4: ax.text(i, bottom + h / 2, f"{h:.1f}%", ha="center", va="center", color="white", fontsize=9.5, fontweight="bold")
        bottom += h
    # PAZy is < 0.2% at every level, too thin to see: mark the top of the bar and label it
    ax.plot([i - 0.275, i + 0.275], [100, 100], color=PAL[2], lw=3, solid_capstyle="butt")
    ax.text(i, 103, f"PAZy-centred: {tot(L, 'PL'):,} ({tot(L, 'PL') / t * 100:.2g}%)", ha="center", fontsize=8.6, color="#B45309")
ax.set_xticks(range(3)); ax.set_xticklabels([f"{LV_LABEL[L]}\n{tot(L, 'clusters'):,} centroids" for L in LEVELS], fontsize=9.5)
ax.set_ylim(0, 110); ax.set_yticks(range(0, 101, 20)); ax.set_ylabel("Share of cluster centroids (%)")
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.16), ncol=3, fontsize=9.5)
finish(fig, "fig2_source_of_centroids", "Where each level's centroids come from",
       "Centroid source by sequence ID prefix; sources clustered in priority order PAZy → nr → Logan")

# ---- 3. domain composition ----------------------------------------------------------------------------------------
order = sorted(NAMES, key=lambda d: -tot("L90", "clusters", [d]) if ("L90", d) in S else 0)
top = order[:7]
colour = {d: PAL[i] for i, d in enumerate(top)}
fig, ax = plt.subplots(figsize=(8.6, 5.2)); fig.subplots_adjust(top=0.84, bottom=0.14, left=0.09, right=0.70)
for i, L in enumerate(LEVELS):
    t = tot(L, "clusters"); bottom = 0
    segs = [(d, tot(L, "clusters", [d])) for d in top if (L, d) in S]
    segs.append(("other", t - sum(v for _, v in segs)))
    for d, v in segs:
        h = v / t * 100
        ax.bar(i, h, bottom=bottom, width=0.58, color=colour.get(d, GREY), edgecolor="white", linewidth=0.8)
        if h >= 3: ax.text(i, bottom + h / 2, f"{h:.1f}%", ha="center", va="center", fontsize=8.5, color="white", fontweight="bold")
        bottom += h
    ax.text(i, -7, f"{t:,} clusters", ha="center", fontsize=8.5, color=SUB)
ax.set_xticks(range(3)); ax.set_xticklabels([LV_LABEL[L] + (" †" if len(doms[L]) < 15 else "") for L in LEVELS])
ax.tick_params(axis="x", pad=18)
ax.set_ylim(0, 100); ax.set_ylabel("Share of clusters (%)")
for s in ("top", "right"): ax.spines[s].set_visible(False)
other_n = len(NAMES) - len(top)
handles = [plt.Rectangle((0, 0), 1, 1, color=colour[d]) for d in top] + [plt.Rectangle((0, 0), 1, 1, color=GREY)]
ax.legend(handles, [NAMES[d] for d in top] + [f"Other ({other_n} domains)"], frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
finish(fig, "fig3_domain_composition", "CATH domain composition of the clusters",
       "Share of clusters at each identity level; domains ordered by 90% cluster count")

# ---- 5. singletons vs v1.1: share above the axis, count mirrored below ----------------------------------------------
V11 = {}   # clusters, singletons: Benjamin's v1.1 run (data/v11_reference.tsv)
with open(os.path.join(DATA, "v11_reference.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"): V11[r["level"]] = (int(r["clusters"]), int(r["singletons"]))
RUNS = [("v2", PAL[0]), ("v1.1", PAL[2])]
def nice_ticks(m):
    step = 10 ** np.floor(np.log10(m)); step = step * (5 if m / step > 5 else 2 if m / step > 2 else 1) / 2
    return np.arange(step, m * 1.0001, step)
fig, axs = plt.subplots(1, 3, figsize=(10.4, 5.6)); fig.subplots_adjust(top=0.80, bottom=0.07, left=0.08, right=0.98, wspace=0.42)
PMAX = 60
for ax, L in zip(axs, LEVELS):
    vals = {"v2": (tot(L, "singletons"), tot(L, "clusters")), "v1.1": (V11[L][1], V11[L][0])}
    cmax = max(v[0] for v in vals.values()) * 1.15
    for i, (run, c) in enumerate(RUNS):
        n, t = vals[run]; pct = n / t * 100
        ax.bar(i, pct / PMAX, width=0.6, color=c, label=run if L == "L90" else None); ax.bar(i, -n / cmax, width=0.6, color=c, alpha=0.55)
        ax.text(i, pct / PMAX + 0.02, f"{pct:.1f}%", ha="center", va="bottom", fontsize=9)
        ax.text(i, -n / cmax - 0.02, compact(n), ha="center", va="top", fontsize=9)
    ax.axhline(0, color=INK, lw=0.9)
    up = np.arange(0, PMAX + 1, 20); dn = nice_ticks(cmax / 1.15)
    ax.set_yticks(list(up / PMAX) + list(-dn / cmax)); ax.set_yticklabels([f"{u:g}%" for u in up] + [millions(d) for d in dn], fontsize=8.5)
    ax.set_ylim(-1.08, 1.12); ax.set_xticks([])
    ax.set_title(LV_LABEL[L], fontsize=10.5, loc="left", fontweight="bold")
    for sp in ("top", "right", "bottom"): ax.spines[sp].set_visible(False)
axs[0].set_ylabel("singleton count  ↓        ↑  % of clusters")
fig.legend(*axs[0].get_legend_handles_labels(), frameon=False, loc="upper right", bbox_to_anchor=(0.98, 0.955), ncol=2)
finish(fig, "fig5_singletons_vs_v11", "Singleton clusters, v2 vs v1.1",
       "Above the axis: singletons as a share of all clusters. Below (mirrored): number of singleton clusters, own scale per level")

print("wrote figures to", FIG)

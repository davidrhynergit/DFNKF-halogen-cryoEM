"""Ideal-strand handedness test: per-peptide CC difference (final map minus mirrored map).
Usage: python3 make_ideal_delta.py out_prefix "label1" res1_prefix [...]"""
import sys, csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
out = sys.argv[1]; sets = list(zip(sys.argv[2::2], sys.argv[3::2]))
INK, MUTED, L_C, D_C = "#222222", "#8A8A8A", "#2B5C8A", "#D9822B"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6})
fig, ax = plt.subplots(figsize=(3.8, 2.5)); rng = np.random.default_rng(0)
ticks = []
for k, (lab, pre) in enumerate(sets):
    r = list(csv.DictReader(open(pre + "_ideal_strand.csv")))
    dL = np.array([float(x["L_final"]) - float(x["L_mirror"]) for x in r])
    dD = np.array([float(x["D_final"]) - float(x["D_mirror"]) for x in r])
    for off, d, col in ((-0.18, dL, L_C), (0.18, dD, D_C)):
        x = k + off + rng.uniform(-0.07, 0.07, len(d))
        ax.scatter(x, d, s=5, color=col, alpha=0.6, linewidths=0, zorder=2)
        ax.plot([k + off - 0.12, k + off + 0.12], [np.median(d)] * 2, color=INK, lw=1.4, zorder=3)
    p = wilcoxon(dL).pvalue
    ax.text(k - 0.18, 0.43, f"{(dL > 0).sum()}/{len(dL)}\np={p:.0e}".replace("e-0", "e-"), ha="center", va="bottom", fontsize=5.5, color=L_C)
    ticks.append(lab.replace("\\n", "\n"))
ax.axhline(0, color=MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)
ax.set_xticks(range(len(sets))); ax.set_xticklabels(ticks, fontsize=6)
ax.set_ylabel("CC final map − CC mirrored map"); ax.set_ylim(-0.42, 0.52)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.scatter([], [], s=12, color=L_C, label="L-strand"); ax.scatter([], [], s=12, color=D_C, label="D-strand (control)")
ax.legend(loc="lower left", fontsize=6, frameon=False, ncol=2, handletextpad=0.2, borderaxespad=0.1)
fig.tight_layout(); fig.savefig(out + ".png", dpi=300); fig.savefig(out + ".pdf"); print("saved", out)

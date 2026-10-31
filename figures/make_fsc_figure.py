"""Figure: half-map FSC (RELION post-processing) and model-map FSC (Phenix, masked) for the refined structures.
Usage: python3 make_fsc_figure.py Validation_folder out_prefix"""
import sys, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
V, out = sys.argv[1], sys.argv[2]
SETS = [("9x4", "4 × 9 DF(I)NKF"), ("10x12", "10 × 12 DF(Br)NKF(Br)"), ("8x4", "4 × 8 DF(I)NKF"), ("10x4", "4 × 10 DF(I)NKF")]
PLACEHOLDER_MODEL = set()
INK, MUTED, C1, C2 = "#222222", "#8A8A8A", "#2B5C8A", "#7A4E9C"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6})
def res_ticks(ax):
    t = [10, 5, 4, 3, 2.5]; ax.set_xticks([1/v for v in t]); ax.set_xticklabels([f"{v:g}" for v in t])
    ax.set_xlim(0, 1/2.5); ax.set_ylim(-0.05, 1.05); ax.set_xlabel("Resolution (Å)")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
def last_above(x, y, t):
    k = np.where(y >= t)[0].max(); 
    if k + 1 >= len(x): return 1/x[k]
    xi = x[k] + (t - y[k]) * (x[k+1] - x[k]) / (y[k+1] - y[k]); return 1/xi
fig, axs = plt.subplots(2, 4, figsize=(7.08, 3.6), sharey=True)
for j, (d, lab) in enumerate(SETS):
    a = np.loadtxt(os.path.join(V, d, "postprocess_fsc.dat")); x, y = a[:, 0], a[:, 1]
    i = np.where((y < 0.143) & (np.arange(len(y)) > 3))[0][0]; r = 1/x[i-1]
    ax = axs[0, j]; ax.plot(x, y, color=C1, lw=1.2); ax.axhline(0.143, color=MUTED, lw=0.6, ls=(0, (3, 2)))
    ax.text(0.03, 0.165, "0.143", color=MUTED, fontsize=6, ha="left", va="bottom")
    ax.text(0.97, 0.9, f"{r:.2f} Å", fontsize=7, color=INK, ha="right", transform=ax.transAxes)
    ax.set_title(lab, fontsize=7.5); res_ticks(ax)
    ax.text(-0.18 if j == 0 else -0.08, 1.06, "abcd"[j], transform=ax.transAxes, fontsize=10, fontweight="bold")
    ax = axs[1, j]
    if d in PLACEHOLDER_MODEL:
        ax.text(0.5, 0.5, "PLACEHOLDER\nmodel-map FSC\n(validation to be rerun)", ha="center", va="center", transform=ax.transAxes, fontsize=6.5, color=MUTED)
    else:
        L = open(os.path.join(V, d, "model_map_fsc.txt")).read().split("Unmasked")[0].splitlines()[2:]
        m = np.array([list(map(float, l.split())) for l in L if l.strip()])
        ax.plot(m[:, 0], m[:, 1], color=C2, lw=1.0); ax.axhline(0.5, color=MUTED, lw=0.6, ls=(0, (3, 2)))
        ax.text(0.03, 0.52, "0.5", color=MUTED, fontsize=6, ha="left", va="bottom")
        ax.text(0.97, 0.9, f"{last_above(m[:,0], m[:,1], 0.5):.2f} Å", fontsize=7, color=INK, ha="right", transform=ax.transAxes)
    res_ticks(ax)
    ax.text(-0.18 if j == 0 else -0.08, 1.06, "efgh"[j], transform=ax.transAxes, fontsize=10, fontweight="bold")
axs[0, 0].set_ylabel("Half-map FSC"); axs[1, 0].set_ylabel("Model-map FSC")
fig.tight_layout(h_pad=1.6); fig.savefig(out + ".png", dpi=300); fig.savefig(out + ".pdf"); print("saved", out)

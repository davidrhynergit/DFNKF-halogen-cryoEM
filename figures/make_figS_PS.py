"""Figure S20 a-f: averaged power spectra (a-d) and 2D class averages at 1.3 A/px with the beta-strand spacing (e, f).

Usage: python3 make_figS_PS.py DATA_DIR [OUT_PREFIX]
DATA_DIR must contain (as in the revision working folder):
  AvgPS-job038_I/AvgPS-job038.png, AvgPS-job038_vertical.csv    (DF(I)NKF averaged power spectrum)
  AvgPS-job270_Br/AvgPS-job270.png, AvgPS-job270_vertical.csv   (DF(Br)NKF(Br))
  Particle_counts_I/Class2D_polymorph_clusters/run_it025_classes_225.mrcs
  Particle_counts_Br/Final_Classes/Br_10x12_continued_run_it035_class052.mrcs
The averaged power spectra were calculated with avg-powerspec.py (RELION helper script) from up to
50,000 extracted segments at 1.3 A/px.
"""
import csv, os, sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
from PIL import Image

U = os.path.join(sys.argv[1], "") if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "")
SETS = [("DF(I)NKF", U+"AvgPS-job038_I/AvgPS-job038"), ("DF(Br)NKF(Br)", U+"AvgPS-job270_Br/AvgPS-job270")]
BOX, APIX = 256, 1.3
# 2D class for panel e/f: (stack, 1-based class, label). Fibril axis horizontal. None = placeholder.
# DF(I)NKF: 4 x 9 polymorph, Class2D job225, class 22 (5.8% of the 4 x 9 segments, best estimated resolution 3.74 A).
# DF(Br)NKF(Br): 10 x 12 polymorph, continued 2D run (run_it035), class 52 (clear 4.8 A striations; stored as single-image stack).
CLASSES = [(U + "Particle_counts_I/Class2D_polymorph_clusters/run_it025_classes_225.mrcs", 22, "4 × 9, class 22", 160, 128),
           (U + "Particle_counts_Br/Final_Classes/Br_10x12_continued_run_it035_class052.mrcs", 1, "10 × 12, class 52", 208, 248)]
ZOOM, ZX, ZY = 64, 160, 128  # zoom box (px) and its centre in the class image (px)

def read_mrcs(fn, idx):
    h = np.fromfile(fn, dtype=np.int32, count=256); nx, ny = int(h[0]), int(h[1])
    return np.fromfile(fn, dtype=np.float32, offset=1024 + int(h[23]) + (idx - 1) * nx * ny * 4, count=nx * ny).reshape(ny, nx)
INK, MUTED, LINE = "#222222", "#888888", "#2B5C8A"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.edgecolor": INK})

def ps_image(prefix):
    im = np.asarray(Image.open(prefix + ".png").convert("L")).astype(float)
    left = im[:, : im.shape[1] // 2]
    nw = left < 250
    rows = np.where(nw.mean(1) > 0.35)[0]; cols = np.where(nw.mean(0) > 0.35)[0]
    crop = left[rows.min()+2:rows.max()-1, cols.min()+2:cols.max()-1]
    lo, hi = np.percentile(crop, [1, 99.7])
    return np.clip((crop - lo) / (hi - lo), 0, 1)

def profile(prefix):
    rows = list(csv.reader(open(prefix + "_vertical.csv")))[1:]
    px = np.array([int(r[0]) for r in rows]); p = np.array([float(r[2]) for r in rows])
    k = px > 0
    return px[k] / (BOX * APIX), np.log10(p[k]), px[k]

fig = plt.figure(figsize=(7.08, 7.4))  # 180 mm wide
gs = fig.add_gridspec(3, 2, height_ratios=[1.05, 0.75, 1.0], hspace=0.42, wspace=0.28,
                      left=0.08, right=0.98, top=0.96, bottom=0.05)
f_ll = 69 / (BOX * APIX)  # 4.82 A bin
for j, (name, pre) in enumerate(SETS):
    # a,b  2D spectra
    ax = fig.add_subplot(gs[0, j]); img = ps_image(pre); n = img.shape[0]
    ax.imshow(img, cmap="gray", interpolation="lanczos")
    c = n / 2; d = 69 / BOX * n
    for s in (-1, 1):
        ax.annotate("", xy=(c + s*d*1.04, c + 0.17*n), xytext=(c + s*d*1.3, c + 0.33*n),
                    arrowprops=dict(arrowstyle="-|>", color="#F2C14E", lw=1.0))
        ax.text(c + s*d*1.3, c + 0.35*n, "4.8 Å", color="#F2C14E", ha="center", va="top", fontsize=7)
    ax.text(0.03*n, 0.97*n, "fibril axis ↔", color="white", fontsize=6.5, va="bottom")
    ax.set_axis_off(); ax.set_title(name, fontsize=8, color=INK, pad=4)
    ax.text(-0.04, 1.04, "ab"[j], transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
    # c,d  meridional profiles
    ax = fig.add_subplot(gs[1, j]); f, lp, px = profile(pre)
    ax.plot(f, lp, color=LINE, lw=1.2)
    ax.axvline(f_ll, color=MUTED, lw=0.7, ls=(0, (3, 2)))
    ax.text(f_ll, ax.get_ylim()[1], " 4.8 Å layer line", color=INK, fontsize=6.5, va="top", ha="left")
    ax.axvline(1/(2*APIX), color=MUTED, lw=0.7, ls=":")
    ax.text(1/(2*APIX), ax.get_ylim()[1], "Nyquist\n2.6 Å ", color=MUTED, fontsize=6, va="top", ha="right")
    ax.set_xlim(0.03, 1/(2*APIX)+0.005)
    ax.set_xlabel("Spatial frequency along fibril axis (1/Å)"); ax.set_ylabel("log$_{10}$ power" if j == 0 else "")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.text(-0.14, 1.06, "cd"[j], transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
    # e,f  2D class average, zoom on the beta-strand striations, power spectrum of the axial profile
    ax = fig.add_subplot(gs[2, j]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()
    ax.text(-0.04, 1.0, "ef"[j], transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
    if CLASSES[j] is None:
        ax.add_patch(Rectangle((0, 0.08), 0.62, 0.9, fc="#E6E6E6", ec=MUTED, lw=0.6, ls="--"))
        ax.text(0.31, 0.86, "PLACEHOLDER\n2D class average\n%s, 1.3 Å/px\n(fibril axis horizontal)" % name,
                ha="center", va="top", fontsize=6.5, color=INK)
        ax.add_patch(Rectangle((0.68, 0.40), 0.32, 0.58, fc="#E6E6E6", ec=MUTED, lw=0.6, ls="--"))
        ax.text(0.84, 0.69, "zoom:\nβ-strand\nstriations\n(4.8 Å)", ha="center", va="center", fontsize=6.5, color=INK)
        ax.add_patch(Rectangle((0.68, 0.08), 0.32, 0.26, fc="#E6E6E6", ec=MUTED, lw=0.6, ls="--"))
        ax.text(0.84, 0.21, "axial profile\nof zoom", ha="center", va="center", fontsize=6.5, color=INK)
        continue
    fn, idx, lab, ZX, ZY = CLASSES[j]; im = read_mrcs(fn, idx); n = im.shape[0]
    lo, hi = np.percentile(im, [0.5, 99.8]); disp = np.clip((im - lo) / (hi - lo), 0, 1)
    x0, y0 = ZX - ZOOM // 2, ZY - ZOOM // 2
    a1 = ax.inset_axes([0, 0.08, 0.62, 0.9]); a1.imshow(disp, cmap="gray", interpolation="lanczos"); a1.set_axis_off()
    a1.add_patch(Rectangle((x0, y0), ZOOM, ZOOM, fc="none", ec="#F2C14E", lw=1.0))
    a1.plot([n * 0.06, n * 0.06 + 50 / APIX], [n * 0.93] * 2, color="white", lw=1.5, solid_capstyle="butt")
    a1.text(n * 0.06 + 25 / APIX, n * 0.90, "5 nm", color="white", fontsize=6.5, ha="center", va="bottom")
    a1.text(n * 0.04, n * 0.04, lab, color="white", fontsize=6.5, va="top")
    a1.text(n * 0.96, n * 0.04, "fibril axis ↔", color="white", fontsize=6.5, va="top", ha="right")
    z = im[y0:y0 + ZOOM, x0:x0 + ZOOM]
    a2 = ax.inset_axes([0.68, 0.40, 0.32, 0.58]); zl, zh = np.percentile(z, [0.5, 99.5])
    a2.imshow(np.clip((z - zl) / (zh - zl), 0, 1), cmap="gray", interpolation="lanczos"); a2.set_xticks([]); a2.set_yticks([])
    for sp in a2.spines.values(): sp.set_color("#F2C14E"); sp.set_linewidth(1.0)
    for yb, ya in ((y0, 0), (y0 + ZOOM, 1)):  # connectors from the box corners to the zoom
        fig.add_artist(ConnectionPatch(xyA=(x0 + ZOOM, yb), coordsA=a1.transData, xyB=(0, 1 - ya), coordsB=a2.transAxes,
                                       color="#F2C14E", lw=0.6))
    prof = z.mean(0) - z.mean()                      # average across the fibril: profile along the axis
    NP = 4096; fr = np.fft.rfftfreq(NP, d=APIX); P = np.abs(np.fft.rfft(prof * np.hanning(ZOOM), n=NP)) ** 2
    a3 = ax.inset_axes([0.68, 0.08, 0.32, 0.26]); k = fr > 0.03
    a3.plot(fr[k], P[k] / P[k].max(), color=LINE, lw=1.0)
    a3.axvline(1 / 4.8, color=MUTED, lw=0.7, ls=(0, (3, 2)))
    pk = fr[k][np.argmax(P[k] * (fr[k] > 0.15))]
    a3.text(1 / 4.8, 1.02, "%.1f Å" % (1 / pk), fontsize=6, color=INK, ha="center", va="bottom")
    a3.set_xlim(0.03, 1 / (2 * APIX)); a3.set_ylim(0, 1.15); a3.set_yticks([])
    a3.set_xticks([0.1, 0.2, 0.3]); a3.tick_params(labelsize=5.5, length=2)
    a3.set_xlabel("1/Å along fibril axis", fontsize=5.5, labelpad=1)
    for sp in ("top", "right", "left"): a3.spines[sp].set_visible(False)
    print("%s: strongest axial frequency above 0.15/Å in the zoom: %.4f/Å = %.2f Å" % (name, pk, 1 / pk))

out = sys.argv[2] if len(sys.argv) > 2 else "FigS_PS"
fig.savefig(out + ".png", dpi=300); fig.savefig(out + ".pdf", dpi=600)
print("saved", out)

"""Panel e/f for Figure S[PS]: 2D class average with zoom on beta-strand striations.

Usage:
  python3 make_strand_panel.py classes.mrcs CLASS_INDEX APIX OUT_PREFIX [--zoom X0 Y0 SIZE] [--title TEXT]

CLASS_INDEX is 1-based (as in 000005@classes.mrcs). The fibril axis must be horizontal
(use your Horizontal/ classes). --zoom gives the zoom box in pixels (default: centre, 1/4 box).
Needs numpy and matplotlib only.
"""
import sys, argparse, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def read_mrcs(fn):
    h = np.fromfile(fn, dtype=np.int32, count=256)
    nx, ny, nz, mode = h[:4]; next_ = h[23]
    dt = {0: np.int8, 1: np.int16, 2: np.float32, 6: np.uint16, 12: np.float16}[int(mode)]
    a = np.fromfile(fn, dtype=dt, offset=1024 + int(next_), count=int(nx * ny * nz))
    return a.reshape(nz, ny, nx).astype(np.float32)

ap = argparse.ArgumentParser(); ap.add_argument("mrcs"); ap.add_argument("idx", type=int)
ap.add_argument("apix", type=float); ap.add_argument("out")
ap.add_argument("--zoom", nargs=3, type=int); ap.add_argument("--title", default="")
a = ap.parse_args()

img = read_mrcs(a.mrcs)[a.idx - 1]; n = img.shape[0]
x0, y0, s = a.zoom if a.zoom else (n // 2 - n // 8, n // 2 - n // 8, n // 4)
z = img[y0:y0 + s, x0:x0 + s]
prof = z.mean(0) - z.mean()                      # average across the fibril -> profile along the axis
f = np.fft.rfftfreq(s, d=a.apix); P = np.abs(np.fft.rfft(prof * np.hanning(s))) ** 2

plt.rcParams.update({"font.size": 7, "axes.linewidth": 0.6})
fig = plt.figure(figsize=(3.5, 2.4))
g = fig.add_gridspec(2, 2, width_ratios=[1.6, 1], hspace=0.5, wspace=0.25, left=0.02, right=0.98, top=0.9, bottom=0.14)
ax = fig.add_subplot(g[:, 0]); ax.imshow(img, cmap="gray"); ax.set_axis_off()
ax.add_patch(Rectangle((x0, y0), s, s, fc="none", ec="#F2C14E", lw=1))
ax.plot([n * 0.05, n * 0.05 + 50 / a.apix], [n * 0.95] * 2, color="white", lw=1.5)
ax.text(n * 0.05, n * 0.92, "5 nm", color="white", fontsize=6)
if a.title: ax.set_title(a.title, fontsize=7)
ax = fig.add_subplot(g[0, 1]); ax.imshow(z, cmap="gray"); ax.set_axis_off()
ax = fig.add_subplot(g[1, 1]); ax.plot(f[1:], P[1:], color="#2B5C8A", lw=1)
ax.axvline(1 / 4.8, color="#888", lw=0.7, ls=(0, (3, 2))); ax.text(1 / 4.8, ax.get_ylim()[1], " 4.8 Å", fontsize=6, va="top")
ax.set_xlim(0, max(0.25, 1/(2*a.apix))); ax.set_yticks([]); ax.set_xlabel("1/Å along axis", fontsize=6)
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
if 1 / (2 * a.apix) < 1 / 4.8:
    print(f"WARNING: Nyquist {2*a.apix:.1f} Å is coarser than 4.8 Å; strands cannot be resolved at this pixel size.")
pk = f[1:][np.argmax(P[1:] * (f[1:] > 0.15))]
print(f"strongest axial frequency above 0.15/Å: {pk:.4f}/Å = {1/pk:.2f} Å")
fig.savefig(a.out + ".png", dpi=300); fig.savefig(a.out + ".pdf"); print("saved", a.out)

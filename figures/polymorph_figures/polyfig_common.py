"""Shared drawing helpers for the final polymorph figures S13 (DF(I)NKF) and S14 (DF(Br)NKF(Br))."""
import os, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'pdf.fonttype': 42,
                     'axes.linewidth': 0.5, 'svg.fonttype': 'none'})
INK, MUTED, RULE, LIGHT = '#1f1f1f', '#6b6b6b', '#c8c8c8', '#ededed'
MAG, OR, ORL, BLUE = '#B43090', '#E89030', '#F5CFA0', '#1E6EB4'
W = 7.09            # 180 mm
fmt = lambda n: f'{int(n):,}'


def read_mrcs(fn):
    h = np.fromfile(fn, dtype=np.int32, count=256)
    nx, ny, nz = int(h[0]), int(h[1]), int(h[2])
    return np.fromfile(fn, dtype=np.float32, offset=1024 + int(h[23]), count=nx * ny * nz).reshape(nz, ny, nx)


def read_model(path):
    """class number -> rlnClassDistribution from a RELION *_model.star."""
    txt = open(path).read().split('data_model_classes')[1].split('data_model_class_1')[0]
    cols = [l.split()[0] for l in txt.splitlines() if l.startswith('_rln')]
    rows = [l.split() for l in txt.splitlines() if l.strip() and not l.startswith(('_', 'loop_', '#', 'data_'))]
    ir, idd = cols.index('_rlnReferenceImage'), cols.index('_rlnClassDistribution')
    return {int(r[ir].split('@')[0]): float(r[idd]) for r in rows}


def norm(im, lo=0.5, hi=99.8):
    im = np.asarray(im, dtype=np.float32); a, b = np.percentile(im, [lo, hi])
    return np.clip((im - a) / (b - a), 0, 1)


def png_gray(path):
    return np.asarray(Image.open(path).convert('L')).astype(np.float32) / 255.0


def new_page(height):
    fig = plt.figure(figsize=(W, height))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, height); ax.axis('off')
    return fig, ax


def tile(ax, img, x, y, s_w, s_h=None, vmin=0, vmax=1):
    """Grayscale image with its lower-left corner at (x, y), in inches."""
    s_h = s_w if s_h is None else s_h
    ax.imshow(img, cmap='gray', vmin=vmin, vmax=vmax, extent=(x, x + s_w, y, y + s_h), interpolation='none', zorder=2)


def hrule(ax, y, x0=0.0, x1=W, lw=0.4, color=RULE):
    ax.plot([x0, x1], [y, y], color=color, lw=lw, solid_capstyle='butt', zorder=1)


def scalebar(ax, x_right, y, length_in, label, color='white', fs=5.5):
    ax.plot([x_right - length_in, x_right], [y, y], color=color, lw=1.2, solid_capstyle='butt', zorder=4)
    ax.text(x_right - length_in / 2, y + 0.025, label, color=color, ha='center', va='bottom', fontsize=fs, zorder=4)


def ctext(ax, x, y, parts, **kw):
    """Text in several colours on one line, starting at x."""
    r = ax.figure.canvas.get_renderer()
    for t, col in parts:
        h = ax.text(x, y, t, color=col, **kw)
        x = h.get_window_extent(renderer=r).transformed(ax.transData.inverted()).x1
    return x


def load_inimodel(xs_path, side_path, crossover):
    """Initial-model images (cross-section, side view). Each side view shows exactly two crossovers,
    so the pixel size of both screenshots is width / (2 x crossover)."""
    xs, side = png_gray(xs_path), png_gray(side_path)
    px_per_A = side.shape[1] / (2.0 * crossover)
    return xs, side, px_per_A


def crop_center(img, n):
    h, w = img.shape; n = min(n, h, w); y0, x0 = (h - n) // 2, (w - n) // 2
    return img[y0:y0 + n, x0:x0 + n]


def crossover_positions(side, crossover, px_per_A):
    """x positions (px) of two neighbouring identical views in a side view, one crossover apart.
    The view with the strongest stripe pattern across the fibril (rows or peptides edge-on) is used as
    the anchor; the projection repeats every crossover (180 degree rotation)."""
    im = side - side.mean()
    # high-pass across the fibril (vertical), then stripe energy per column
    k = np.hanning(9); k /= k.sum()
    sm = np.apply_along_axis(lambda v: np.convolve(v, k, mode='same'), 0, im)
    e = ((im - sm) ** 2).sum(0)
    c = crossover * px_per_A
    kk = np.hanning(max(5, int(c / 12)) | 1); kk /= kk.sum()
    e = np.convolve(e, kk, mode='same')
    xs = np.arange(len(e)); ph = np.mod(xs, c)
    nb = 90; bins = np.linspace(0, c, nb + 1); idx = np.minimum(np.digitize(ph, bins) - 1, nb - 1)
    prof = np.array([e[idx == b].mean() if np.any(idx == b) else 0 for b in range(nb)])
    prof = np.convolve(np.r_[prof[-3:], prof, prof[:3]], np.ones(5) / 5, mode='same')[3:-3]
    x0 = (np.argmax(prof) + 0.5) * c / nb
    return x0, x0 + c


def crossover_bar(ax, x_left_in, y_in, length_in, label=None, color='#1f1f1f', fs=6, tick=0.035):
    ax.plot([x_left_in, x_left_in + length_in], [y_in, y_in], color=color, lw=0.7, solid_capstyle='butt', zorder=4)
    for x in (x_left_in, x_left_in + length_in):
        ax.plot([x, x], [y_in - tick / 2, y_in + tick / 2], color=color, lw=0.7, zorder=4)
    if label:
        ax.text(x_left_in + length_in + 0.05, y_in, label, fontsize=fs, color=color, va='center', ha='left')

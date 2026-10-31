"""Figure S13 (final): polymorphs of DF(I)NKF.

One row per polymorph: the six most populated 2D class averages (FilamentTools cluster, RELION Class2D,
1.3 A/px, 333 A box), the initial model (side view and cross-section, relion_helix_inimodel2d) and the
segment statistics (1,000 randomly selected micrographs; 95% CI from bootstrap over micrographs).
Last row: fibrils that could not be assigned.

Usage: python3 make_figS13_I_polymorphs.py CLASS_DIR INIMODEL_DIR OUT_PREFIX
  CLASS_DIR     run_it025_classes_<job>.mrcs and run_it025_model_<job>.star (jobs 221-226)
  INIMODEL_DIR  job<NNN>_<n>pep_cross<C>_01.png (cross-section) and _02.png (side view)
"""
import sys, os, glob, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from polyfig_common import *

CLS_DIR, INI_DIR, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
# peptides per row: (Class2D job, fibrils, segments, %, CI lo, CI hi)
P = [(7, '221', 817, 25676, 5.1, 4.5, 5.8), (8, '224', 4712, 152944, 30.6, 29.2, 32.0),
     (9, '225', 5019, 167533, 33.5, 32.4, 34.7), (10, '222', 3591, 119932, 24.0, 22.8, 25.1),
     (11, '223', 913, 33525, 6.7, 6.1, 7.3)]
SOLVED = {9}
FAILED = ('226', 660, 12741)
TOTA = sum(p[3] for p in P); TOT = TOTA + FAILED[2]; FIB = sum(p[2] for p in P) + FAILED[1]
NCLS, PIX, BOX = 6, 1.3, 256

def top_classes(job, n=NCLS):
    d = read_mrcs(os.path.join(CLS_DIR, f'run_it025_classes_{job}.mrcs'))
    dist = read_model(os.path.join(CLS_DIR, f'run_it025_model_{job}.star'))
    order = sorted(dist, key=lambda c: -dist[c])[:n]
    return [d[c - 1] for c in order], order

def inimodel(npep):
    f = glob.glob(os.path.join(INI_DIR, f'job*_{npep}pep_cross*_01.png'))[0]
    cross = int(re.search(r'cross(\d+)', f).group(1))
    xs, side, pxA = load_inimodel(f, f.replace('_01.png', '_02.png'), cross)
    return xs, side, pxA, cross

# ---- layout (inches) ----
S, G = 0.56, 0.035                 # class tile and gap
IN_PER_A = S / (BOX * PIX)         # classes and side views share this scale
XC = 0.66                          # classes start
XCS, CS, XS_A = 4.30, 0.90, 230.0  # cross-section column, size (in) and field of view (A)
XF, XS, XP = 5.58, 6.22, 6.34      # right edge fibrils, right edge segments, left edge %
SH = 0.47 + 0.11                   # height reserved for the side view and its crossover bar
RH = S + G + SH                    # row height (classes + side view)
PAD = 0.11
HEAD = 0.42
H = HEAD + 5 * (RH + 2 * PAD) + (S + 0.2 + 2 * PAD) + 0.3
fig, ax = new_page(H)

# header
yh = H - 0.05
ax.text(0.02, yh, 'Polymorph', fontsize=7, fontweight='bold', va='top')
ax.text(XC, yh, '2D class averages (six most populated)', fontsize=7, fontweight='bold', va='top')
ax.text(XC, yh - 0.14, 'and initial model (side view)', fontsize=7, fontweight='bold', va='top')
ax.text(XCS + CS / 2, yh, 'Initial model', fontsize=7, fontweight='bold', va='top', ha='center')
ax.text(XCS + CS / 2, yh - 0.14, '(cross-section)', fontsize=7, fontweight='bold', va='top', ha='center')
ax.text(XF, yh, 'Fibrils', fontsize=7, fontweight='bold', va='top', ha='right')
ax.text(XS, yh, 'Segments', fontsize=7, fontweight='bold', va='top', ha='right')
ax.text(XP, yh, '% of segments', fontsize=7, fontweight='bold', va='top')
ax.text(XP, yh - 0.14, '(95% CI)', fontsize=7, fontweight='bold', va='top')
hrule(ax, H - HEAD, lw=0.8, color=INK)

def class_row(imgs, y_top):
    for j, im in enumerate(imgs):
        tile(ax, norm(im), XC + j * (S + G), y_top - S, S)

y = H - HEAD
for i, (npep, job, fib, seg, pct, lo, hi) in enumerate(P):
    yt = y - PAD; yc = yt - RH / 2
    ax.text(0.02, yc, f'4 × {npep}' + (' *' if npep in SOLVED else ''), fontsize=9, fontweight='bold', va='center')
    imgs, _ = top_classes(job); class_row(imgs, yt)
    xs, side, pxA, cross = inimodel(npep)
    sw, sh = side.shape[1] / pxA * IN_PER_A, side.shape[0] / pxA * IN_PER_A
    ys = yt - S - G - (SH - 0.11) / 2
    tile(ax, side, XC, ys - sh / 2, sw, sh)
    x0, x1 = crossover_positions(side, cross, pxA)          # one crossover, between two identical views
    crossover_bar(ax, XC + x0 / pxA * IN_PER_A, ys - sh / 2 - 0.06, cross * IN_PER_A, f'crossover {cross} Å', fs=6)
    tile(ax, crop_center(xs, int(round(XS_A * pxA))), XCS, yc - CS / 2, CS)
    if i == len(P) - 1:
        scalebar(ax, XCS + CS - 0.05, yc - CS / 2 + 0.06, 50 / XS_A * CS, '50 Å')
    ax.text(XF, yc, fmt(fib), fontsize=7, va='center', ha='right')
    ax.text(XS, yc, fmt(seg), fontsize=7, va='center', ha='right')
    ax.text(XP, yc + 0.17, f'{pct:.1f}%', fontsize=8, fontweight='bold', va='center')
    ax.text(XP, yc + 0.03, f'({lo:.1f}–{hi:.1f})', fontsize=6.5, color=MUTED, va='center')
    bw = 0.62 / 40                              # bar: 0-40 %
    ax.add_patch(matplotlib.patches.Rectangle((XP, yc - 0.2), pct * bw, 0.09, facecolor=MAG, edgecolor='none'))
    ax.plot([XP + lo * bw, XP + hi * bw], [yc - 0.155] * 2, color=INK, lw=0.6)
    if i == 0:   # scale bar on the last class of the first row: 100 A
        scalebar(ax, XC + 6 * (S + G) - G - 0.04, yt - S + 0.05, 100 / (BOX * PIX) * S, '100 Å')
    y = yt - RH - PAD
    hrule(ax, y)

# not assigned
yt = y - PAD; yc = yt - (S + 0.2) / 2
ax.text(0.02, yc + 0.07, 'Not', fontsize=8, fontweight='bold', va='center')
ax.text(0.02, yc - 0.08, 'assigned', fontsize=8, fontweight='bold', va='center')
imgs, _ = top_classes(FAILED[0]); class_row(imgs, yt)
ax.text(XC, yt - S - G - 0.06, 'mostly short fibrils; no initial model', fontsize=6, color=MUTED, va='top')
ax.text(XCS + CS / 2, yc, '–', fontsize=8, color=MUTED, va='center', ha='center')
ax.text(XF, yc, fmt(FAILED[1]), fontsize=7, va='center', ha='right')
ax.text(XS, yc, fmt(FAILED[2]), fontsize=7, va='center', ha='right')
ax.text(XP, yc + 0.07, f'{100 * FAILED[2] / TOT:.1f}% of all', fontsize=7, va='center')
ax.text(XP, yc - 0.08, 'segments', fontsize=7, va='center')
y = yt - S - 0.2 - PAD
hrule(ax, y, lw=0.8, color=INK)
ax.text(0.02, y - 0.08, f'Total: 1,000 randomly selected micrographs (of 3,771); {fmt(FIB)} fibrils, {fmt(TOT)} segments, '
        f'{fmt(TOTA)} ({100 * TOTA / TOT:.1f}%) assigned. * Structure solved at high resolution.',
        fontsize=6.3, color=MUTED, va='top')
fig.savefig(OUT + '.pdf'); fig.savefig(OUT + '.png', dpi=300)
print('saved', OUT, f'{W:.2f} x {H:.2f} in')

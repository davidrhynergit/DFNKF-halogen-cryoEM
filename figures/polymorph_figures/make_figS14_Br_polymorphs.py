"""Figure S14 (final): polymorphs of DF(Br)NKF(Br), two pages.

Page 1: a) per-fibril assignment workflow and outcome, b) distribution of all 27 polymorphs
        (95% CI from bootstrap over micrographs).
Page 2: c) one row per polymorph with interpretable classes (18): the six most populated selected
        2D class averages (per-polymorph Class2D, 1.3 A/px, 512 px = 666 A box, shown binned to 256 px),
        the initial model (side view, crossover, cross-section) and the segment statistics.

Usage: python3 make_figS14_Br_polymorphs.py DATA_DIR INIMODEL_DIR OUT_PREFIX
  DATA_DIR      Br_distribution_tube280.json, sets_summary_tube280.csv, br_cls6.npz, br_cls6.json
  INIMODEL_DIR  A<NN>_cross<C>_001.png (cross-section) and _002.png (side view)
"""
import sys, os, glob, re, json, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from polyfig_common import *
from matplotlib.patches import FancyBboxPatch, Rectangle

DATA, INI, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
dist = json.load(open(os.path.join(DATA, 'Br_distribution_tube280.json')))
P = dist['assigned']; TOT = dist['total_segments']; TOTA = dist['assigned_segments']
summ = list(csv.DictReader(open(os.path.join(DATA, 'sets_summary_tube280.csv'))))
TIE = [r for r in summ if r['file'].startswith('reclassify/')]
SHORT = [r for r in summ if r['file'].startswith('too_short/')]
seg = lambda rs: sum(int(r['segments']) for r in rs)
fibn = lambda rs: sum(int(r['fibrils']) for r in rs)
CLS = np.load(os.path.join(DATA, 'br_cls6.npz'))
INFO = json.load(open(os.path.join(DATA, 'br_cls6.json')))
A_OF = {v['polymorph']: a for a, v in INFO.items()}               # '10x12' -> 'A02'
SOLVED = {'10x12'}
rc = lambda k: tuple(int(x) for x in k.split('x'))
FIB_A = sum(p['fibrils'] for p in P)
FIB_ALL = FIB_A + fibn(TIE) + fibn(SHORT)
assert TOTA + seg(TIE) + seg(SHORT) == TOT, 'segment totals do not add up'

def label(ax, x, y, s): ax.text(x, y, s, fontsize=10, fontweight='bold', color=INK, va='top')
def poly_label(ax, x, y, k, fs=7, star=True, **kw):
    r, c = rc(k)
    return ctext(ax, x, y, [(f'{r}', BLUE), (f' × {c}' + (' *' if (star and k in SOLVED) else ''), INK)], fontsize=fs, **kw)

# =============== page 1: workflow + distribution ===============
H1 = 8.7
fig, ax = new_page(H1)
label(ax, 0.0, H1 - 0.02, 'a')
steps = ['Manual picking:\none fibril = one filament\n' + f'{fmt(FIB_ALL)} fibrils',
         'Extraction:\n1024 px (666 Å) box,\n47.5 Å inter-box distance\n' + f'{fmt(TOT)} segments',
         '2D classification of all\nsegments (tube 280 Å);\nclasses labelled by view',
         'Per fibril: majority in\neach view gives rows and\npeptides per row',
         'Fibrils tied in one\nview classified again']
bw, gap, y0, bh = 1.22, 0.2, H1 - 1.2, 0.88
for i, t in enumerate(steps):
    x = 0.06 + i * (bw + gap)
    ax.add_patch(FancyBboxPatch((x, y0), bw, bh, boxstyle='round,pad=0,rounding_size=0.06', facecolor='#fbf4ec', edgecolor='#e3c9a6', lw=0.6))
    ax.text(x + bw / 2, y0 + bh / 2, t, ha='center', va='center', fontsize=6.3, linespacing=1.3)
    if i < len(steps) - 1:
        ax.annotate('', xy=(x + bw + gap - 0.02, y0 + bh / 2), xytext=(x + bw + 0.02, y0 + bh / 2),
                    arrowprops=dict(arrowstyle='-|>', color=MUTED, lw=0.6, mutation_scale=7))
ow, oy = 2.12, y0 - 0.85
for j, (fc, ec, t) in enumerate([
        ('#fbe7d0', OR, f'Assigned in both views\n{fmt(FIB_A)} fibrils, {fmt(TOTA)} segments ({100 * TOTA / TOT:.1f}%)'),
        ('#f2f2f2', MUTED, f'Tied in one view after reclassification\n{fmt(fibn(TIE))} fibrils, {fmt(seg(TIE))} segments ({100 * seg(TIE) / TOT:.1f}%)'),
        ('#f2f2f2', MUTED, f'Too short to show both views\n{fmt(fibn(SHORT))} fibrils, {fmt(seg(SHORT))} segments ({100 * seg(SHORT) / TOT:.1f}%)')]):
    x = 0.06 + j * (ow + 0.24)
    ax.add_patch(FancyBboxPatch((x, oy), ow, 0.52, boxstyle='round,pad=0,rounding_size=0.06', facecolor=fc, edgecolor=ec, lw=0.6))
    ax.text(x + ow / 2, oy + 0.26, t, ha='center', va='center', fontsize=6.3, linespacing=1.35)

label(ax, 0.0, oy - 0.25, 'b')
bx = fig.add_axes([0.95 / W, 0.62 / H1, (W - 1.95) / W, (oy - 0.92) / H1])
n = len(P); ys = list(range(n))[::-1]
for y, p in zip(ys, P):
    shown = p['polymorph'] in A_OF
    bx.barh(y, p['pct'], height=0.64, color=OR if shown else ORL, lw=0)
    bx.plot([p['lo'], p['hi']], [y, y], color=INK, lw=0.6)
    for v in (p['lo'], p['hi']): bx.plot([v, v], [y - 0.15, y + 0.15], color=INK, lw=0.6)
    pt = f"{p['pct']:.1f}%" if p['pct'] >= 0.05 else '<0.1%'
    bx.text(max(p['hi'], p['pct']) + 0.35, y, f"{pt}  ({fmt(p['fibrils'])})", va='center', fontsize=6, color=INK)
bx.set_ylim(-0.7, n - 0.3); bx.set_xlim(0, 24)
bx.set_yticks(ys); bx.set_yticklabels([])
fig.canvas.draw()
to_in = ax.transData.inverted()
for y, p in zip(ys, P):          # two-colour labels, left-aligned in page coordinates
    _, yi = to_in.transform(bx.transData.transform((0, y)))
    poly_label(ax, 0.3, yi, p['polymorph'], fs=7, va='center')
bx.set_xlabel('% of assigned segments (95% CI); number of fibrils in brackets', fontsize=6.8)
for s in ('top', 'right'): bx.spines[s].set_visible(False)
for s in ('left', 'bottom'): bx.spines[s].set_color('#b0b0b0'); bx.spines[s].set_linewidth(0.5)
bx.tick_params(axis='x', colors=MUTED, width=0.5, length=2.5, labelsize=6.5); bx.tick_params(axis='y', length=0)
bx.xaxis.grid(True, color='#e6e6e6', lw=0.4); bx.set_axisbelow(True)
ax.text(W - 0.02, 0.04, 'Light bars: polymorphs without interpretable 2D classes (not shown in c). * Structure solved at high resolution.',
        ha='right', va='bottom', fontsize=6, color=MUTED)
fig.savefig(OUT + '_page1.pdf'); fig.savefig(OUT + '_page1.png', dpi=300); plt.close(fig)

# =============== page 2: table of the 18 polymorphs ===============
S, G = 0.38, 0.022
BOXA = 666.0                                  # class box in A
IN_PER_A = S / BOXA
XC = 0.62
XSV = XC + 6 * (S + G) - G + 0.09             # side view start
XCR = XSV + 2550 * IN_PER_A + 0.43            # right edge of crossover column
XCS, CS, XS_A = XCR + 0.09, S, 360.0          # cross-section column
XF, XSG, XP = XCS + CS + 0.42, XCS + CS + 0.92, XCS + CS + 0.99
ROW = S + 0.085
SHOWN = [p for p in P if p['polymorph'] in A_OF]
HEAD = 0.42
H2 = 0.22 + HEAD + len(SHOWN) * ROW + 0.36
fig, ax = new_page(H2)
label(ax, 0.0, H2 - 0.02, 'c')
yh = H2 - 0.22
for x, t, ha in [(0.02, 'Polymorph', 'left'), (XC, '2D class averages (six most populated)', 'left'),
                 (XSV, 'Initial model (side view)', 'left'), (XCR, 'Crossover', 'right'),
                 (XCS + CS / 2, 'Cross-', 'center'), (XF, 'Fibrils', 'right'), (XSG, 'Segments', 'right'), (XP, '% of', 'left')]:
    ax.text(x, yh, t, fontsize=6.5, fontweight='bold', va='top', ha=ha)
for x, t, ha in [(XCR, '(Å)', 'right'), (XCS + CS / 2, 'section', 'center'), (XP, 'segments', 'left')]:
    ax.text(x, yh - 0.13, t, fontsize=6.5, fontweight='bold', va='top', ha=ha)
hrule(ax, yh - HEAD + 0.12, lw=0.8, color=INK)
y = yh - HEAD + 0.12
for i, p in enumerate(SHOWN):
    k = p['polymorph']; a = A_OF[k]
    yt = y - 0.045; yc = yt - S / 2
    poly_label(ax, 0.02, yc, k, fs=7.5, va='center', fontweight='bold')
    for j, im in enumerate(CLS[a]):
        tile(ax, norm(im), XC + j * (S + G), yt - S, S)
    f = glob.glob(os.path.join(INI, f'{a}_cross*_001.png'))[0]
    cross = int(re.search(r'cross(\d+)', f).group(1))
    xs, side, pxA = load_inimodel(f, f.replace('_001.png', '_002.png'), cross)
    sw, sh = side.shape[1] / pxA * IN_PER_A, side.shape[0] / pxA * IN_PER_A
    ysv = yc + 0.025
    tile(ax, norm(side, 0.5, 99.7), XSV, ysv - sh / 2, sw, sh)
    x0, x1 = crossover_positions(side, cross, pxA)          # one crossover, between two identical views
    crossover_bar(ax, XSV + x0 / pxA * IN_PER_A, ysv - sh / 2 - 0.028, cross * IN_PER_A, None, tick=0.024)
    ax.text(XCR, yc, fmt(cross), fontsize=6.5, va='center', ha='right')
    tile(ax, norm(crop_center(xs, int(round(XS_A * pxA))), 0.5, 99.7), XCS, yc - CS / 2, CS)
    ax.text(XF, yc, fmt(p['fibrils']), fontsize=6.5, va='center', ha='right')
    ax.text(XSG, yc, fmt(p['segments']), fontsize=6.5, va='center', ha='right')
    ax.text(XP, yc + 0.06, f"{p['pct']:.1f}%", fontsize=6.8, fontweight='bold', va='center')
    ax.text(XP, yc - 0.075, f"({p['lo']:.1f}–{p['hi']:.1f})", fontsize=5.5, color=MUTED, va='center')
    if i == 0:
        scalebar(ax, XC + 6 * (S + G) - G - 0.03, yt - S + 0.04, 200 * IN_PER_A, '200 Å', fs=5)
    if i == len(SHOWN) - 1:
        scalebar(ax, XCS + CS - 0.09, yc - CS / 2 + 0.035, 100 / XS_A * CS, '100 Å', fs=4.5)
    y = yt - S - 0.04
    hrule(ax, y, lw=0.3 if i < len(SHOWN) - 1 else 0.8, color=RULE if i < len(SHOWN) - 1 else INK)
ax.text(0.02, y - 0.07, f'{len(SHOWN)} of {len(P)} polymorphs ({100 * sum(p["segments"] for p in SHOWN) / TOTA:.1f}% of the assigned segments); '
        'the other polymorphs have no interpretable 2D classes (light bars in b).', fontsize=5.8, color=MUTED, va='top')
ax.text(0.02, y - 0.2, 'Class averages and side views are shown at the same scale. * Structure solved at high resolution.',
        fontsize=5.8, color=MUTED, va='top')
fig.savefig(OUT + '_page2.pdf'); fig.savefig(OUT + '_page2.png', dpi=300); plt.close(fig)
print('saved', OUT, f'page1 {W:.2f} x {H1:.2f} in, page2 {W:.2f} x {H2:.2f} in')

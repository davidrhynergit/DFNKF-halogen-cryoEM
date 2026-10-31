"""Unbiased handedness test with an ideal beta-strand (no conformation taken from either map).
An ideal parallel beta-strand pentapeptide (backbone + CB; L or D) is fitted as a rigid body at every
peptide position of the central layer, in the final map and in its mirror image, with the same search.
Usage: python3 ideal_strand_test.py map.mrc model.cif n_layers out_prefix"""
import sys, numpy as np, re
from scipy.ndimage import map_coordinates, gaussian_filter
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation as R
MAP, CIF, NL, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
RES = 3.0
# ---------- map ----------
h = np.fromfile(MAP, dtype=np.int32, count=256); hf = np.fromfile(MAP, dtype=np.float32, count=256)
nx, ny, nz = h[:3]; ap = float(hf[10] / h[7])
a = np.fromfile(MAP, dtype=np.float32, offset=1024 + int(h[23]), count=nx*ny*nz).reshape(nz, ny, nx)
N = nz; af = a[::-1].copy()
# ---------- model (only to locate the peptides) ----------
lines = open(CIF).read().splitlines(); cols = []; rows = []; i = 0
while i < len(lines):
    if lines[i].strip() == 'loop_' and lines[i+1].strip().startswith('_atom_site.'):
        i += 1
        while lines[i].strip().startswith('_atom_site.'): cols.append(lines[i].strip().split('.')[1]); i += 1
        while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith(('_', 'loop_', '#')):
            rows.append(lines[i].split()); i += 1
        break
    i += 1
c = {k: j for j, k in enumerate(cols)}
ca = {}
for r in rows:
    if r[c['label_atom_id']] == 'CA':
        ca.setdefault(r[c['auth_asym_id']], {})[int(r[c['auth_seq_id']])] = [float(r[c['Cartn_x']]), float(r[c['Cartn_y']]), float(r[c['Cartn_z']])]
chains = [k for k, v in ca.items() if len(v) == 5]
zc = {k: np.mean([p[2] for p in ca[k].values()]) for k in chains}
order = sorted(chains, key=lambda k: zc[k]); per = len(order) // NL
central = order[(NL // 2) * per:(NL // 2 + 1) * per]
# ---------- ideal strand (NeRF) ----------
def place(a_, b_, c_, bond, ang, tor):
    bc = c_ - b_; bc /= np.linalg.norm(bc); n = np.cross(b_ - a_, bc); n /= np.linalg.norm(n)
    m = np.array([bc, np.cross(n, bc), n]).T
    d2 = np.array([-bond*np.cos(np.radians(ang)), bond*np.sin(np.radians(ang))*np.cos(np.radians(tor)), bond*np.sin(np.radians(ang))*np.sin(np.radians(tor))])
    return c_ + m @ d2
def strand(n=5, phi=-119, psi=113):
    Nn, CA, C = np.array([0., 1.458, 0]), np.array([0., 0, 0]), np.array([1.525, 0, 0])
    Nn = place(C, CA, Nn, 1.458, 111.2, 0) if False else Nn
    bb = [Nn, CA, C]; atoms = []
    for k in range(n):
        Nk, CAk, Ck = bb[-3], bb[-2], bb[-1]
        O = place(Nk, CAk, Ck, 1.231, 120.8, psi + 180)
        b = CAk - Nk; cc = Ck - CAk; aa = np.cross(b, cc)
        CB = -0.58273431*aa + 0.56802827*b - 0.54067466*cc + CAk   # L-amino acid
        atoms += [('N', Nk), ('CA', CAk), ('C', Ck), ('O', O), ('CB', CB)]
        if k < n - 1:
            N2 = place(Nk, CAk, Ck, 1.329, 116.2, psi)
            CA2 = place(CAk, Ck, N2, 1.458, 121.7, 180)
            C2 = place(Ck, N2, CA2, 1.525, 111.2, phi)
            bb += [N2, CA2, C2]
    return atoms
S = strand(); names = [s[0] for s in S]; XL = np.array([s[1] for s in S])
XD = XL * np.array([1, 1, -1])                          # D-amino acid strand (mirror image)
w = np.array([{'N': 7, 'C': 6, 'O': 8}.get(nm[0], 6) for nm in names], float)
caidx = [i for i, nm in enumerate(names) if nm == 'CA']
# ---------- scoring ----------
interp = lambda m, P: map_coordinates(m, (P[:, 2]/ap, P[:, 1]/ap, P[:, 0]/ap), order=1)
def ccmask(m, P):
    lo = np.maximum(np.floor(P.min(0)/ap).astype(int) - 6, 0); hi = np.minimum(np.ceil(P.max(0)/ap).astype(int) + 6, N)
    sub = m[lo[2]:hi[2], lo[1]:hi[1], lo[0]:hi[0]]; g = np.zeros(sub.shape, np.float32)
    f = (P - lo*ap)/ap; i0 = np.floor(f).astype(int); d = f - i0
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                ww = w*np.where(dx, d[:, 0], 1-d[:, 0])*np.where(dy, d[:, 1], 1-d[:, 1])*np.where(dz, d[:, 2], 1-d[:, 2])
                ii = i0 + [dx, dy, dz]; ok = (ii >= 0).all(1) & (ii[:, 0] < g.shape[2]) & (ii[:, 1] < g.shape[1]) & (ii[:, 2] < g.shape[0])
                np.add.at(g, (ii[ok, 2], ii[ok, 1], ii[ok, 0]), ww[ok])
    s = gaussian_filter(g, 0.225*RES/ap); mk = s > 0.05*s.max()
    x = s[mk] - s[mk].mean(); y = sub[mk] - sub[mk].mean()
    return float((x*y).sum()/np.sqrt((x*x).sum()*(y*y).sum()))
def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0); H = (P - Pc).T @ (Q - Qc); U, _, Vt = np.linalg.svd(H)
    dd = np.sign(np.linalg.det(Vt.T @ U.T)); Rm = Vt.T @ np.diag([1, 1, dd]) @ U.T
    return Rm, Qc - Rm @ Pc
ms_f, ms_m = gaussian_filter(a, 1.0), gaussian_filter(af, 1.0)
def fit(X0, target_ca, m, msm, rng):
    best = None
    for rev in (False, True):
        tc = target_ca[::-1] if rev else target_ca
        Rm, t = kabsch(X0[caidx], tc); Xs = X0 @ Rm.T + t; cen = Xs.mean(0)
        for k in range(6):
            p0 = np.zeros(6) if k == 0 else np.r_[rng.normal(0, np.radians(15), 3), rng.normal(0, 1.0, 3)]
            if k in (1, 2):  # 180 deg about the strand axis
                ax = tc[-1] - tc[0]; ax /= np.linalg.norm(ax); p0[:3] = np.pi*ax
            xf = lambda p: R.from_rotvec(p[:3]).apply(Xs - cen) + cen + p[3:]
            r = minimize(lambda p: -np.average(interp(msm, xf(p)), weights=w), p0, method='Powell')
            r = minimize(lambda p: -np.average(interp(m, xf(p)), weights=w), r.x, method='Powell')
            if best is None or r.fun < best[0]: best = (r.fun, xf(r.x))
    return best[1]
rng = np.random.default_rng(1); res = []
for ch in central:
    T = np.array([ca[ch][k] for k in sorted(ca[ch])]); Tm = T.copy(); Tm[:, 2] = (N - 1)*ap - Tm[:, 2]
    row = [ch]
    for X0 in (XL, XD):
        Pf = fit(X0, T, a, ms_f, rng); Pm = fit(X0, Tm, af, ms_m, rng)
        row += [ccmask(a, Pf), ccmask(af, Pm)]
    res.append(row)
R_ = np.array([r[1:] for r in res], float)
with open(OUT + '_ideal_strand.csv', 'w') as f:
    f.write('chain,L_final,L_mirror,D_final,D_mirror\n')
    for r in res: f.write(r[0] + ',' + ','.join(f'{v:.4f}' for v in r[1:]) + '\n')
txt = (f"{MAP}\ncentral layer: {len(central)} peptides; ideal parallel beta-strand (backbone + CB), rigid-body fit, same search in both maps\n"
       f"mean CC  L-strand: final {R_[:,0].mean():.3f}  mirror {R_[:,1].mean():.3f}   (L better in final map for {(R_[:,0]>R_[:,1]).sum()} of {len(R_)})\n"
       f"mean CC  D-strand: final {R_[:,2].mean():.3f}  mirror {R_[:,3].mean():.3f}   (D better in mirror map for {(R_[:,3]>R_[:,2]).sum()} of {len(R_)})\n")
open(OUT + '_ideal_strand.txt', 'w').write(txt); print(txt)

"""Halogen...O contacts of the halogen atoms in the central layer of a fibril model (overview).

For every halogen atom X (I or Br) of the central layer, all backbone O atoms of other chains
within 4.5 A are listed, grouped by halogen residue, acceptor residue and layer offset, with the
X...O distance and the C(Z)-X...O angle.
Usage: python3 halogen_geom.py model.cif X n_layers      (X = I or BR, as in label_atom_id)
The model must be a fibril stack along z with the same number of chains per layer.
"""
import sys, numpy as np, collections
from scipy.spatial import cKDTree
def read(fn):
    L=open(fn).read().splitlines(); cols=[]; rows=[]; i=0
    while i<len(L):
        if L[i].strip()=='loop_' and L[i+1].strip().startswith('_atom_site.'):
            i+=1
            while L[i].strip().startswith('_atom_site.'): cols.append(L[i].strip().split('.')[1]); i+=1
            while i<len(L) and L[i].strip() and not L[i].strip().startswith(('_','loop_','#')): rows.append(L[i].split()); i+=1
            break
        i+=1
    c={k:j for j,k in enumerate(cols)}
    return [(r[c['auth_asym_id']],int(r[c['auth_seq_id']]),r[c['label_comp_id']],r[c['label_atom_id']],np.array([float(r[c['Cartn_x']]),float(r[c['Cartn_y']]),float(r[c['Cartn_z']])])) for r in rows]


def main():
    fn, X, NL = sys.argv[1], sys.argv[2], int(sys.argv[3])
    A=read(fn)
    pos={(a[0],a[1],a[3]):a[4] for a in A}
    chz=collections.defaultdict(list)
    for a in A:
        if a[3]=='CA': chz[a[0]].append(a[4][2])
    chains=sorted(chz,key=lambda k:np.mean(chz[k])); per=len(chains)//NL
    layer={ch:i//per for i,ch in enumerate(chains)}
    mid=NL//2
    O=[(a[0],a[1]) for a in A if a[3]=='O']; OP=np.array([pos[(c,r,'O')] for c,r in O]); tree=cKDTree(OP)
    res=collections.defaultdict(list)
    for a in A:
        if a[3]!=X or layer[a[0]]!=mid: continue
        xp=a[4]; cz=pos[(a[0],a[1],'CZ')]
        for j in tree.query_ball_point(xp,4.5):
            oc,orr=O[j]
            if oc==a[0]: continue
            d=np.linalg.norm(OP[j]-xp); v1=cz-xp; v2=OP[j]-xp
            ang=np.degrees(np.arccos(v1@v2/np.linalg.norm(v1)/np.linalg.norm(v2)))
            res[(a[1],orr,layer[oc]-mid)].append((d,ang))
    print(fn.split('/')[-1], 'central layer', mid+1,'of',NL)
    for k in sorted(res,key=lambda k:-len(res[k])):
        v=np.array(res[k]); short=v[v[:,0]<=4.0]
        print(f"  X on res {k[0]} -> O of res {k[1]}, layer offset {k[2]:+d}: n={len(v)}  d {v[:,0].min():.2f}-{v[:,0].max():.2f} (median {np.median(v[:,0]):.2f})  C-X..O {v[:,1].min():.0f}-{v[:,1].max():.0f} (median {np.median(v[:,1]):.0f})  n(d<=4.0, angle>=140)={((v[:,0]<=4.0)&(v[:,1]>=140)).sum()}")


if __name__ == "__main__":
    main()

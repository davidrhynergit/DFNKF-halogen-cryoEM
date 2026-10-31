"""Halogen...O contacts in a fibril model: per contact type, central layer.
Usage: python3 halogen_table.py model.cif X n_layers vdw_sum
  X        halogen atom name as in label_atom_id (I or BR)
  vdw_sum  sum of the van der Waals radii of X and O in A (e.g. 3.50 for I...O, 3.37 for Br...O)
Contacts <= 4.2 A to backbone O of other chains; angles C(Z)-X...O and X...O=C."""
import sys, numpy as np, collections
from scipy.spatial import cKDTree
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from halogen_geom import read
fn,X,NL,VDW=sys.argv[1],sys.argv[2],int(sys.argv[3]),float(sys.argv[4])
A=read(fn); pos={(a[0],a[1],a[3]):a[4] for a in A}; resn={(a[0],a[1]):a[2] for a in A}
chz=collections.defaultdict(list)
for a in A:
    if a[3]=='CA': chz[a[0]].append(a[4][2])
chains=sorted(chz,key=lambda k:np.mean(chz[k])); per=len(chains)//NL; layer={ch:i//per for i,ch in enumerate(chains)}; mid=NL//2
O=[(a[0],a[1]) for a in A if a[3]=='O']; OP=np.array([pos[k+('O',)] for k in O]); tree=cKDTree(OP)
def ang(a,b,c):
    v1=a-b; v2=c-b; return np.degrees(np.arccos(v1@v2/np.linalg.norm(v1)/np.linalg.norm(v2)))
rows=collections.defaultdict(list); nX=0
for a in A:
    if a[3]!=X or layer[a[0]]!=mid: continue
    nX+=1; xp=a[4]; cz=pos[(a[0],a[1],'CZ')]
    for j in tree.query_ball_point(xp,4.2):
        oc,orr=O[j]
        if oc==a[0]: continue
        d=np.linalg.norm(OP[j]-xp); t1=ang(cz,xp,OP[j]); t2=ang(xp,OP[j],pos[(oc,orr,'C')])
        dl=layer[oc]-mid
        rows[(f"{resn[(a[0],a[1])]}{a[1]}", f"{resn[(oc,orr)]}{orr} O", 'same layer' if dl==0 else ('adjacent layer' if abs(dl)==1 else f'{abs(dl)} layers'))].append((d,t1,t2))
print(f"{fn.split('/')[-1]}: {nX} {X} atoms in the central layer (layer {mid+1} of {NL}); contacts <= 4.2 A; vdW sum {VDW} A")
print("halogen | acceptor | relation | n | d range | d median (IQR) | C-X..O median (range) | X..O=C median | n(d<vdW) | n(d<=vdW+0.3 & C-X..O>=150)")
for k in sorted(rows,key=lambda k:-len(rows[k])):
    v=np.array(rows[k]); d=v[:,0]
    print(f"{k[0]} | {k[1]} | {k[2]} | {len(v)} | {d.min():.2f}-{d.max():.2f} | {np.median(d):.2f} ({np.percentile(d,25):.2f}-{np.percentile(d,75):.2f}) | {np.median(v[:,1]):.0f} ({v[:,1].min():.0f}-{v[:,1].max():.0f}) | {np.median(v[:,2]):.0f} | {(d<VDW).sum()} | {((d<=VDW+0.3)&(v[:,1]>=150)).sum()}")

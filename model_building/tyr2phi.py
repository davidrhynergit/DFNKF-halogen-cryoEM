import sys, re, numpy as np
src, dst = sys.argv[1], sys.argv[2]
CI = 2.02   # C-I ideal distance in the PHI dictionary used for the 9-layer model
lines = open(src).read().split('\n')
# locate atom_site columns
cols=[]; start=None
for i,l in enumerate(lines):
    if l.startswith('_atom_site.'):
        if start is None: start=i
        cols.append(l.split('.',1)[1].strip())
    elif start is not None and cols and not l.startswith('_atom_site.'):
        first=i; break
ix={c:k for k,c in enumerate(cols)}
end=first
while end<len(lines) and (lines[end].startswith('ATOM') or lines[end].startswith('HETATM')): end+=1
# CZ coordinates of TYR
cz={}
for i in range(first,end):
    t=lines[i].split()
    if t[ix['label_comp_id']]=='TYR' and t[ix['label_atom_id']]=='CZ':
        cz[(t[ix['auth_asym_id']],t[ix['auth_seq_id']])]=np.array([float(t[ix['Cartn_x']]),float(t[ix['Cartn_y']]),float(t[ix['Cartn_z']])])
nI=0; Ipos=[]
def subtok(line, k, new):
    # replace k-th whitespace token, keep column alignment where possible
    spans=[m.span() for m in re.finditer(r'\S+', line)]
    a,b=spans[k]; w=b-a
    if len(new)<w: new=new.ljust(w) if k+1<len(spans) and line[b]==' ' and False else new.rjust(w) if re.match(r'^-?\d',new) else new.ljust(w)
    if len(new)>w and k+1<len(spans):
        # eat spaces after the token if available
        extra=len(new)-w; nxt=spans[k+1][0]
        if nxt-b-1>=extra: return line[:a]+new+line[b+extra:]
    return line[:a]+new+line[b:]
for i in range(first,end):
    l=lines[i]; t=l.split()
    if t[ix['label_comp_id']]!='TYR': continue
    l=subtok(l, ix['label_comp_id'], 'PHI')
    if t[ix['label_atom_id']]=='OH':
        key=(t[ix['auth_asym_id']],t[ix['auth_seq_id']])
        o=np.array([float(t[ix['Cartn_x']]),float(t[ix['Cartn_y']]),float(t[ix['Cartn_z']])])
        v=o-cz[key]; p=cz[key]+v/np.linalg.norm(v)*CI; Ipos.append((key,p))
        l=subtok(l, ix['label_atom_id'], 'I ')
        l=subtok(l, ix['type_symbol'], 'I')
        for c,val in zip(('Cartn_x','Cartn_y','Cartn_z'),p): l=subtok(l, ix[c], '%.3f'%val)
        nI+=1
    lines[i]=l
# other blocks
out=[]
skip_fixed=False
for i,l in enumerate(lines):
    if first<=i<end: out.append(l); continue
    if l.startswith('_audit_syntax.fixed_width'): continue   # widths may change; avoid fixed-width parsing
    if l.startswith('_entity_poly.nstd_monomer'): l='_entity_poly.nstd_monomer yes'
    if l.startswith('_entity_poly.pdbx_seq_one_letter_code_can'):
        l='_entity_poly.pdbx_seq_one_letter_code D(PHI)NKF\n_entity_poly.pdbx_seq_one_letter_code_can DFNKF'
    l=re.sub(r'\bTYR\b','PHI',l)
    out.append(l)
s='\n'.join(out)
s=s.replace('PHI "peptide linking" ','PHI "L-peptide linking" ',1)
s=s.replace('loop_\n_atom_type.symbol\nH \nC \nN \nO \n','loop_\n_atom_type.symbol\nC \nN \nO \nI \n')
open(dst,'w').write(s)
# report: closest non-bonded contacts of the new iodines
X=[];lab=[]
for i in range(first,end):
    t=lines[i].split()
    X.append([float(t[ix['Cartn_x']]),float(t[ix['Cartn_y']]),float(t[ix['Cartn_z']])]); lab.append((t[ix['auth_asym_id']],t[ix['auth_seq_id']],t[ix['label_atom_id']],t[ix['label_comp_id']]))
X=np.array(X); worst=[]
for key,p in Ipos:
    d=np.linalg.norm(X-p,axis=1)
    for j in np.argsort(d)[:8]:
        if (lab[j][0],lab[j][1])==key: continue
        worst.append((round(float(d[j]),2),key,lab[j])); break
worst.sort()
print(dst.split('/')[-1], 'iodines:',nI, 'TYR left:', s.count('TYR'), 'PHI atoms:', sum(1 for l in s.split('\n') if ' PHI ' in l and (l.startswith('ATOM') or l.startswith('HETATM'))))
print(' closest I-to-other-residue contacts:', worst[:6])
print(' I-contacts < 3.3 A:', sum(1 for w in worst if w[0]<3.3))

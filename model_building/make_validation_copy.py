"""Make a validation copy of a fibril model in which the halogenated residues are read as PHE.
MolProbity (Phenix) only evaluates Ramachandran, rotamers and CaBLAM for residues it recognises as
amino acids. PHI (4-iodo-Phe) and 4BF (4-bromo-Phe) are not in that list, so every residue next to
them loses its phi or psi and is skipped (DF(I)NKF: only Lys4 is evaluated; DF(Br)NKF(Br): nothing).
This copy renames PHI/4BF to PHE, removes the halogen atoms and writes all atoms as ATOM records, so
that the backbone and side-chain statistics can be calculated. Use the real model for clashscore and
bond/angle deviations.
Usage: python3 make_validation_copy.py model.cif model_for_rama.cif"""
import sys
src, dst = sys.argv[1], sys.argv[2]
L = open(src).read().splitlines(); cols = []; rows = []; i = 0
while i < len(L):
    if L[i].strip() == 'loop_' and L[i+1].strip().startswith('_atom_site.'):
        i += 1
        while L[i].strip().startswith('_atom_site.'): cols.append(L[i].strip().split('.')[1]); i += 1
        while i < len(L) and L[i].strip() and not L[i].strip().startswith(('_', 'loop_', '#')):
            rows.append(L[i].split()); i += 1
        break
    i += 1
c = {k: j for j, k in enumerate(cols)}
out = ['data_validation_copy', 'loop_'] + ['_atom_site.' + k for k in cols]
n = 0; renamed = set()
for r in rows:
    comp = r[c['label_comp_id']]; el = r[c['type_symbol']].upper()
    if comp in ('PHI', '4BF'):
        if el in ('I', 'BR'): continue
        renamed.add((r[c['auth_asym_id']], r[c['auth_seq_id']]))
        r[c['label_comp_id']] = 'PHE'
    if el == 'H': continue
    r[c['group_PDB']] = 'ATOM'; n += 1; r[c['id']] = str(n)
    out.append(' '.join(r))
open(dst, 'w').write('\n'.join(out) + '\n#\n')
print(f"{src} -> {dst}: {n} atoms, {len(renamed)} residues renamed to PHE")

# Cryo-EM analysis code: halogenated DFNKF amyloid fibrils

Code for the cryo-EM part of

> *Emergence of quadruple-helical symmetry in amyloid fibrils driven by single-point halogenation.*
> Nature Communications (under revision). Citation and DOI will be added on publication.

The cryo-EM maps and models are deposited in the EMDB and PDB (accession codes in the paper).
Contact for this code: David Rhyner, ETH Zurich.

## Contents

| Folder | What it does | Paper |
|---|---|---|
| `polymorph_quantification/` | Per-fibril polymorph assignment of DF(Br)NKF(Br) from RELION 2D classes, and bootstrap confidence intervals of the polymorph fractions (resampling micrographs). Includes the RELION input files for DF(Br)NKF(Br), so the numbers can be reproduced. | Methods; Tables S5, S6; Figures S13, S14 |
| `handedness_test/` | Handedness test with an ideal parallel beta-strand (L and D) fitted at every peptide position in the final map and in its mirror image. | Figure S20g |
| `model_building/` | ISOLDE parameters for 4-iodo-Phe (PHI) and 4-bromo-Phe (4BF), ChimeraX/ISOLDE restraint commands, Phenix restraint dictionaries and ADP-refinement settings, helper scripts. | Methods, Table S1 |
| `halogen_geometry/` | Halogen...O contacts (distances and angles) in the central layer of the models. | Table S2 |
| `figures/` | Scripts that draw Figures S13, S14, S15 (FSC curves), S20 (power spectra, beta-strand spacing) and the ChimeraX script for S21. | Figures S13-S15, S20, S21 |

## Requirements

Python 3.9 or newer with NumPy, SciPy, Matplotlib and Pillow (`pip install -r requirements.txt`).
`polymorph_quantification/assign_br_polymorphs.py` needs only the Python standard library;
the bootstrap needs NumPy.
External programs used in the workflow (not part of this repository): RELION 4/5, ChimeraX with
ISOLDE, Phenix, and FilamentTools (DF(I)NKF polymorph separation).

## Quick check: reproduce the DF(Br)NKF(Br) polymorph statistics

```
cd polymorph_quantification
python assign_br_polymorphs.py config_DFBrNKFBr.json --out results
python bootstrap_polymorph_fractions.py --assignment results/filament_assignment.csv --out results/bootstrap.csv
```

This takes about a minute and gives 5,966 assigned fibrils (207,060 segments, 95.0%) in 27
polymorphs. The outputs should match `expected_output/` (Tables S5 and S6).

## Licence

MIT (see `LICENSE`).

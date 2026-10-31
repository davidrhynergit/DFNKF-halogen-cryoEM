# Model building and refinement

Models were built with ModelAngelo, adjusted in Coot and refined as 9-layer fibrils in ISOLDE with
symmetry and secondary-structure restraints. Only ADP (B-factor) refinement was done afterwards in
Phenix (real_space_refine) against the unsharpened map. The two outer layers at each end were removed
before deposition (5 central layers).

- `isolde_parameters/`: residue templates for ISOLDE (atom types of ISOLDE's amberff14SB.xml).
  `PHI_4-iodo-Phe_ISOLDE.xml` (4-iodo-L-phenylalanine, mid-chain) and
  `4BF_4-bromo-Phe_ISOLDE.xml` (4-bromo-L-phenylalanine, mid-chain and C-terminal templates).
  Load in ISOLDE with "Load residue parameters" or `isolde load parameters <file>`.
- `isolde_restraints/`: ChimeraX commands for the ISOLDE distance restraints of each model and for
  selecting the outer layers.
- `phenix/`: restraint dictionaries for PHI and 4BF (eLBOW) and an example parameter file for the
  ADP-only refinement (`real_space_refine` with `refinement.run=adp`).
- `tyr2phi.py`: turns Tyr into 4-iodo-Phe (PHI) in an mmCIF model (used for the 4 x 8 and 4 x 10
  starting models).
- `make_validation_copy.py`: copy of a model with PHI/4BF read as PHE and the halogens removed, so that
  MolProbity evaluates Ramachandran, rotamer and CaBLAM statistics for all residues. Clashscore and
  bond/angle deviations were taken from the real model.

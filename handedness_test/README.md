# Handedness test (Figure S20g)

No conformation is taken from either map. An ideal parallel beta-strand pentapeptide (backbone and
CB atoms, phi -119, psi 113, standard geometry) is built once as L and once as its mirror image D.
It is fitted as a rigid body at every peptide position of the central layer, in the final map and in
the mirror map (map flipped along z), with the same search in both (both strand directions, 180 deg
about the strand axis, random starting perturbations; Powell minimisation on a smoothed and then the
full map). The refined model is used only to locate the peptides (its CA trace, mirrored for the mirror
map). Score: correlation with a 3 A map simulated from the strand, within a mask.

```
python ideal_strand_test.py unsharpened_map.mrc model.cif N_LAYERS out/4x9
python make_ideal_delta.py FigS20g "4 x 9" out/4x9 "10 x 12" out/10x12 "4 x 8" out/4x8 "4 x 10" out/4x10
```

`make_ideal_delta.py` plots the per-peptide CC difference (final minus mirrored map) for the L-strand
and the D-strand control, with medians and a Wilcoxon signed-rank test.

# Figure scripts

| Script | Figure | Input |
|---|---|---|
| `polymorph_figures/make_figS13_I_polymorphs.py` | S13 (DF(I)NKF polymorphs) | Class2D stacks and model.star of each FilamentTools cluster, initial-model images |
| `polymorph_figures/make_figS14_Br_polymorphs.py` | S14 (DF(Br)NKF(Br) polymorphs) | output of `polymorph_quantification`, class averages, initial-model images |
| `make_fsc_figure.py` | S15 (half-map and model-map FSC) | RELION post-processing and phenix.validation_cryoem output per structure |
| `make_figS_PS.py`, `make_strand_panel.py` | S20 a-f (power spectra, beta-strand spacing) | averaged power spectra and 2D class averages |
| `../handedness_test/make_ideal_delta.py` | S20 g | output of `ideal_strand_test.py` |
| `make_figS21.cxc` | S21 (4 x 8, 4 x 9, 4 x 10 cross-sections) | final models and unsharpened maps, run in ChimeraX |

Each script lists its inputs and usage at the top. Panels were assembled in Adobe Illustrator.

# Polymorph quantification

## DF(Br)NKF(Br): per-fibril assignment

`assign_br_polymorphs.py` with `config_DFBrNKFBr.json`.

Every fibril was picked by hand as one continuous filament, so that all its segments share one
`rlnHelicalTubeID`. A fibril is identified by micrograph and helical tube ID. Segments were extracted
in 1024 px boxes (666 A), binned to 2.6 A/px, and all 217,865 segments were classified together in 2D
(RELION Class2D job743, tube diameter 280 A).

A polymorph is a rectangle of R rows x P peptides per row. Two side views define it: the short edge
(view along the peptides of a row) gives P, the long edge (view along the rows) gives R. Classes that
clearly showed one of the views were labelled by hand in one RELION Select job per label
(job745-748: 4, 6, 8, 10 rows; job749-756: 8 to 22 peptides per row).

For each view, a fibril gets the label that holds more than half of its segments in the labelled
classes of that view. All segments of the fibril are assigned with it. Fibrils with a tie in one view
were classified again (Class2D job785-788) and assigned in the missing view by the same rule, using the
classes of the new run that clearly resolved that view (Select job789-793). Fibrils still tied, and
fibrils without any labelled segment in one view (mostly short fibrils), are not assigned.

Outputs: `filament_assignment.csv` (one line per fibril, with the segment counts per label and the
result), `polymorph_summary.csv`, and with `--write-star` one particle STAR file per polymorph.

Input files (in `data_DFBrNKFBr/`): `run_it014_data.star.gz` of Class2D job743, the
`class_averages.star` of each Select job, and the data files of the four reclassification runs.

## Bootstrap confidence intervals (both peptides)

`bootstrap_polymorph_fractions.py` resamples the micrographs with replacement (5,000 repeats, seed 1)
and recalculates the fraction of segments of each polymorph in every resample. The 95% confidence
interval is given by the 2.5th and 97.5th percentiles. Micrographs, not segments, are resampled,
because segments of one fibril, and fibrils on one micrograph, are not independent.

- DF(Br)NKF(Br): `--assignment filament_assignment.csv`
- DF(I)NKF: `--assignment data_DFINKF/DFINKF_fibrils_FilamentTools.csv`. The fibrils on 1,000
  randomly selected micrographs were picked by hand and separated into polymorphs with FilamentTools.
  This table (one line per fibril) was made from the particle STAR file of each cluster with
  `clusters_to_fibril_table.py`, which also checks that no fibril is split over two clusters.

The interval limits depend slightly on the random draw: with other seeds they move by about 0.1
percentage point. The fractions themselves do not change.

`relion_star.py` is a small STAR reader (standard library only) used by both scripts.

## Expected output

`expected_output/` holds the results reported in the paper (Tables S5 and S6, Figures S13 and S14),
as produced by these scripts with seed 1.

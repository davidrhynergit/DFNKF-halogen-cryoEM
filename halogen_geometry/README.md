# Halogen...O contacts

```
python halogen_geom.py model.cif I 5          # overview of all X...O contacts within 4.5 A
python halogen_table.py model.cif I 5 3.50    # table per contact type (Table S2), vdW sum I...O 3.50 A
python halogen_table.py model.cif BR 5 3.37   # Br...O, vdW sum 3.37 A
```

The model must be a fibril stack along z with the same number of chains per layer; contacts are
listed for the halogen atoms of the central layer.

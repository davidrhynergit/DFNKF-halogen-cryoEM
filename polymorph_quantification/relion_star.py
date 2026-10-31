"""Minimal reader for RELION 3.1+ STAR files (no external dependencies).

Only what the polymorph scripts need: read one data block as a list of
dictionaries keyed by column name, and keep the raw lines so that subsets
can be written back in the original format.
"""
import gzip
import os


def read_block(path, block="data_particles"):
    """Return (header_lines, column_names, rows) of a loop_ block.

    header_lines: all lines of the file up to and including the column
        declarations of the requested block (used to write subsets).
    rows: list of (raw_line, dict) tuples.
    For files with a single unnamed block (e.g. class_averages.star),
    pass block="data_".
    """
    header, cols, rows = [], [], []
    in_block = in_loop = False
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            s = line.strip()
            is_row = False
            if in_block and in_loop and s and not s.startswith(("#", "_", "data_", "loop_")):
                t = s.split()
                if len(t) == len(cols):
                    rows.append((line.rstrip("\n"), dict(zip(cols, t))))
                    is_row = True
            if not rows and not is_row:
                header.append(line.rstrip("\n"))
            if is_row:
                continue
            if s.startswith("data_"):
                in_block = (s == block)
                in_loop = False
                if in_block:
                    cols = []
                continue
            if not in_block:
                continue
            if s == "loop_":
                in_loop = True
                continue
            if in_loop and s.startswith("_"):   # RELION and extra columns (e.g. FilamentTools _clusterLabel)
                cols.append(s.split()[0][1:])
                continue
    while header and header[-1].strip() == "":
        header.pop()
    return header, cols, rows


def fibril_key(row):
    """A fibril (filament) is identified by its micrograph and helical tube ID."""
    return (os.path.basename(row["rlnMicrographName"]), row["rlnHelicalTubeID"])


def selected_classes(path):
    """Class numbers listed in a RELION Select job's class_averages.star."""
    _, _, rows = read_block(path, "data_")
    return {int(r["rlnClassNumber"]) for _, r in rows}

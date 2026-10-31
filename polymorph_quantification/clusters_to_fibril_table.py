#!/usr/bin/env python3
"""Per-fibril table from one particle STAR file per polymorph (DF(I)NKF, FilamentTools clusters).

FilamentTools assigns whole filaments to clusters. This script checks that no segment occurs twice
and no fibril is split over two files, and writes one line per fibril in the format that
bootstrap_polymorph_fractions.py reads (--assignment).

Usage:
    python clusters_to_fibril_table.py out.csv ROWSxPEPTIDES=file.star [...] unassigned=file.star
    e.g. 4x7=run_it025_data_221.star 4x8=run_it025_data_224.star ... unassigned=run_it025_data_226.star
"""
import collections
import csv
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from relion_star import read_block, fibril_key  # noqa: E402

out, pairs = sys.argv[1], sys.argv[2:]
seen, fib, nseg = set(), {}, collections.Counter()
for p in pairs:
    label, path = p.split("=", 1)
    for _, r in read_block(path)[2]:
        if r["rlnImageName"] in seen:
            sys.exit(f"segment {r['rlnImageName']} occurs twice")
        seen.add(r["rlnImageName"])
        k = fibril_key(r)
        if fib.setdefault(k, label) != label:
            sys.exit(f"fibril {k} is split over {fib[k]} and {label}")
        nseg[k] += 1
op = gzip.open if out.endswith(".gz") else open
with op(out, "wt", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["micrograph", "helical_tube_id", "n_segments", "rows", "peptides_per_row", "status"])
    for k in sorted(fib):
        lab = fib[k]
        if lab == "unassigned":
            w.writerow([k[0], k[1], nseg[k], "X", "X", "unassigned"])
        else:
            r, pp = lab.split("x")
            w.writerow([k[0], k[1], nseg[k], r, pp, "assigned"])
c = collections.Counter(fib.values())
print(f"{len(fib)} fibrils, {len(seen)} segments:", dict(c))

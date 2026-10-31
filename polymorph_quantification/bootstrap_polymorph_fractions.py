#!/usr/bin/env python3
"""95% confidence intervals of polymorph fractions by resampling micrographs.

The fraction of a polymorph is its share of all assigned segments. To estimate
how much the fractions would change with a different set of micrographs, the
micrographs that contain assigned segments are resampled with replacement
(bootstrap) and the fractions are recalculated in each resample. Micrographs,
not segments, are resampled, because segments of one fibril and fibrils on
one micrograph are not independent. The 95% confidence interval is given by the
2.5th and 97.5th percentiles of the resampled fractions.

Two kinds of input:
  --assignment filament_assignment.csv
        output of assign_br_polymorphs.py (fibrils with status "assigned")
  --star LABEL=file.star [LABEL=file.star ...]
        one RELION particle STAR file per polymorph, e.g. the cluster files
        written by FilamentTools for DF(I)NKF:  4x7=cluster26.star ...

Usage:
    python bootstrap_polymorph_fractions.py --assignment results/filament_assignment.csv
    python bootstrap_polymorph_fractions.py --star 4x7=c26.star 4x8=c45.star ...
Options: --repeats 5000 (default), --seed 1 (default), --out bootstrap.csv
"""
import argparse
import collections
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from relion_star import read_block  # noqa: E402


def from_assignment(path):
    seg = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(path)):
        if r["status"] == "assigned":
            seg[r["micrograph"]][f"{r['rows']}x{r['peptides_per_row']}"] += int(r["n_segments"])
    return seg


def from_star(pairs):
    seg = collections.defaultdict(collections.Counter)
    seen = set()
    for p in pairs:
        label, path = p.split("=", 1)
        _, _, rows = read_block(path)
        for _, r in rows:
            img = r["rlnImageName"]
            if img in seen:
                sys.exit(f"segment {img} occurs in more than one file")
            seen.add(img)
            seg[os.path.basename(r["rlnMicrographName"])][label] += 1
    return seg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--assignment")
    g.add_argument("--star", nargs="+")
    ap.add_argument("--repeats", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="bootstrap.csv")
    a = ap.parse_args()

    seg = from_assignment(a.assignment) if a.assignment else from_star(a.star)
    mics = sorted(seg)
    labels = sorted({k for m in seg.values() for k in m}, key=lambda k: -sum(m[k] for m in seg.values()))
    X = np.array([[seg[m][k] for k in labels] for m in mics], float)   # micrographs x polymorphs
    tot = X.sum(0)
    frac = 100 * tot / tot.sum()

    rng = np.random.default_rng(a.seed)
    R = np.empty((a.repeats, len(labels)))
    for b in range(a.repeats):
        s = X[rng.integers(0, len(mics), len(mics))].sum(0)
        R[b] = 100 * s / s.sum()
    lo, hi = np.percentile(R, [2.5, 97.5], axis=0)
    top = np.bincount(R.argmax(1), minlength=len(labels)) / a.repeats

    print(f"{len(mics)} micrographs, {int(tot.sum())} assigned segments, {a.repeats} resamples (seed {a.seed})")
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["polymorph", "segments", "pct", "ci_low", "ci_high", "fraction_of_resamples_most_abundant"])
        for i, k in enumerate(labels):
            w.writerow([k, int(tot[i]), f"{frac[i]:.2f}", f"{lo[i]:.1f}", f"{hi[i]:.1f}", f"{top[i]:.3f}"])
            print(f"  {k:>6}: {int(tot[i]):7d} segments {frac[i]:5.1f}% (95% CI {lo[i]:.1f}-{hi[i]:.1f})"
                  + (f", most abundant in {100 * top[i]:.1f}% of resamples" if top[i] > 0 else ""))
    if len(labels) > 1:
        d = R[:, 0] - R[:, 1]
        print(f"  {labels[0]} minus {labels[1]}: {frac[0] - frac[1]:.1f} points "
              f"(95% CI {np.percentile(d, 2.5):.1f} to {np.percentile(d, 97.5):.1f})")


if __name__ == "__main__":
    main()

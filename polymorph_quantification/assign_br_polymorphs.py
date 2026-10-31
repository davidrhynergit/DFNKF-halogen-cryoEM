#!/usr/bin/env python3
"""Per-fibril polymorph assignment of DF(Br)NKF(Br) fibrils from RELION 2D classes.

Each fibril was picked manually as one continuous filament, so all segments of
a fibril share one rlnHelicalTubeID. A fibril is identified by its micrograph
and helical tube ID.

The cross-section of each polymorph is a rectangle of R rows x P peptides per
row. Two side views define it:
  - short edge (view along the peptides of a row): gives P (8-22)
  - long edge (view along the rows): gives R (4-10)
Classes of one 2D classification (RELION Class2D, all segments together) were
labelled by hand in RELION Select jobs, one Select job per label.

Assignment rule, applied separately to each view:
  - count the fibril's segments in the labelled classes of that view;
  - if one label holds more than half of these segments, the fibril gets it
    ("majority");
  - if the fibril has segments in labelled classes but no label has more than
    half, the view is a "tie";
  - if the fibril has no segment in a labelled class of that view, the view is
    unknown ("none").
All segments of a fibril are assigned together.

Fibrils with a tie in one view were classified again (one Class2D job per
group) and assigned in the missing view by the same majority rule, using the
classes of the new run that clearly resolved that view.

Fibrils known in both views are assigned to polymorph R x P. Not assigned
are fibrils with a tie in one view, also after reclassification ("tie"), and
fibrils without any labelled segment in one view, mostly short fibrils
("too_short").

Usage:
    python assign_br_polymorphs.py config_DFBrNKFBr.json --out results/
Outputs:
    filament_assignment.csv   one line per fibril: counts per label, result
    polymorph_summary.csv     fibrils and segments per polymorph
    star/                     optional (--write-star): one particle STAR file
                              per polymorph with all segments of its fibrils
"""
import argparse
import collections
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from relion_star import read_block, fibril_key, selected_classes  # noqa: E402

VIEWS = {"rows": "long", "peptides": "short"}


def majority(counter):
    """Return (value, rule) for one view from a Counter of label -> segments."""
    total = sum(counter.values())
    if total == 0:
        return "X", "none"
    value, n = counter.most_common(1)[0]
    if n > total / 2:
        return value, "majority"
    return "X", "tie"


def load_labels(label_cfg, base):
    """{class_number: (view, value)} from {select_star: [view, value]}."""
    lab = {}
    for f, (view, value) in label_cfg.items():
        assert view in VIEWS, view
        for c in selected_classes(os.path.join(base, f)):
            if c in lab and lab[c] != (view, value):
                sys.exit(f"class {c} carries two labels: {lab[c]} and {(view, value)}")
            lab[c] = (view, value)
    return lab


def count_per_fibril(data_star, labels):
    """Segments per fibril, and segments per (view, value) label."""
    _, _, rows = read_block(data_star)
    n = collections.Counter()
    c = collections.defaultdict(lambda: {v: collections.Counter() for v in VIEWS})
    lines = collections.defaultdict(list)
    for raw, r in rows:
        k = fibril_key(r)
        n[k] += 1
        lines[k].append(raw)
        x = labels.get(int(r["rlnClassNumber"]))
        if x:
            c[k][x[0]][x[1]] += 1
    return n, c, lines


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", help="JSON file listing the STAR files and labels")
    ap.add_argument("--data-dir", help="folder the paths in the config are relative to "
                    "(default: folder of the config file)")
    ap.add_argument("--out", default="results", help="output folder")
    ap.add_argument("--write-star", action="store_true", help="write one particle STAR file per polymorph")
    a = ap.parse_args()

    cfg = json.load(open(a.config))
    base = a.data_dir or os.path.dirname(os.path.abspath(a.config))
    os.makedirs(a.out, exist_ok=True)

    # 1) main classification
    labels = load_labels(cfg["labels"], base)
    data_star = os.path.join(base, cfg["classification"])
    n, counts, lines = count_per_fibril(data_star, labels)
    res = {}
    for k in n:
        rows, rrule = majority(counts[k]["rows"])
        pep, prule = majority(counts[k]["peptides"])
        res[k] = {"rows": rows, "rows_rule": rrule, "peptides_per_row": pep, "pep_rule": prule, "source": ""}

    # 2) reclassification of fibrils with a tie in one view
    for rc in cfg.get("reclassifications", []):
        lab = load_labels(rc["labels"], base)
        view = rc["view"]                     # view to resolve: "rows" or "peptides"
        key = "rows" if view == "rows" else "peptides_per_row"
        rule = "rows_rule" if view == "rows" else "pep_rule"
        name = rc.get("name", os.path.basename(rc["data"]))
        rn, rcount, _ = count_per_fibril(os.path.join(base, rc["data"]), lab)
        for k in rn:
            if k not in res:
                sys.exit(f"fibril {k} of {rc['data']} is not in the main classification")
            value, how = majority(rcount[k][view])
            if how == "majority" and value not in rc.get("rejected_values", []):
                res[k][key] = value
                res[k][rule] = f"reclass {name}"
                res[k]["source"] = f"reclass {name}"
            else:
                # still unresolved: keep the tie from the first classification
                res[k]["source"] = f"reclass {name} (unresolved)"

    # 3) output per fibril
    lab_cols = [f"{VIEWS[v]}_{x}" for v in ("peptides", "rows") for x in cfg["values"][v]]
    with open(os.path.join(a.out, "filament_assignment.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["micrograph", "helical_tube_id", "n_segments"] + lab_cols +
                   ["rows", "rows_rule", "peptides_per_row", "pep_rule", "source", "status"])
        for k in sorted(n):
            r = res[k]
            cc = [counts[k][v][x] for v in ("peptides", "rows") for x in cfg["values"][v]]
            assigned = r["rows"] != "X" and r["peptides_per_row"] != "X"
            if assigned:
                status = "assigned"
            elif r["rows_rule"] != "none" and r["pep_rule"] != "none":
                # tie in at least one view, the other view known or tied
                # (also after reclassification)
                status = "tie"
            else:
                status = "too_short"
            r["status"] = status
            w.writerow([k[0], k[1], n[k]] + cc +
                       [r["rows"], r["rows_rule"], r["peptides_per_row"], r["pep_rule"], r["source"], status])

    # 4) summary per polymorph
    seg = collections.Counter(); fib = collections.Counter(); st = collections.Counter(); stf = collections.Counter()
    for k, r in res.items():
        st[r["status"]] += n[k]; stf[r["status"]] += 1
        if r["status"] == "assigned":
            p = (int(r["rows"]), int(r["peptides_per_row"]))
            seg[p] += n[k]; fib[p] += 1
    tot, tot_a = sum(n.values()), sum(seg.values())
    with open(os.path.join(a.out, "polymorph_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rows", "peptides_per_row", "fibrils", "segments", "pct_of_assigned", "pct_of_all"])
        for p in sorted(seg, key=lambda p: -seg[p]):
            w.writerow([p[0], p[1], fib[p], seg[p], f"{100 * seg[p] / tot_a:.2f}", f"{100 * seg[p] / tot:.2f}"])
    print(f"fibrils: {len(n)}, segments: {tot}")
    for s in ("assigned", "tie", "too_short"):
        print(f"  {s:10s} {stf[s]:6d} fibrils {st[s]:8d} segments ({100 * st[s] / tot:.1f}%)")
    print(f"  polymorphs: {len(seg)}")
    for p in sorted(seg, key=lambda p: -seg[p])[:6]:
        print(f"  {p[0]:>2} x {p[1]:>2}: {fib[p]:5d} fibrils {seg[p]:7d} segments {100 * seg[p] / tot_a:5.2f}%")

    # 5) optional STAR files per polymorph
    if a.write_star:
        header, _, _ = read_block(data_star)
        d = os.path.join(a.out, "star"); os.makedirs(d, exist_ok=True)
        groups = collections.defaultdict(list)
        for k, r in res.items():
            if r["status"] == "assigned":
                groups[f"Br_rows{r['rows']}_pep{r['peptides_per_row']}.star"].extend(lines[k])
            else:
                groups[f"unassigned_{r['status']}.star"].extend(lines[k])
        for f, ls in groups.items():
            with open(os.path.join(d, f), "w") as fh:
                fh.write("\n".join(header) + "\n" + "\n".join(ls) + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""genotype_call.py — turn tabular BLAST hits into one genotype call per query.

Reads `blastn -outfmt 6` with a FIXED column layout (see --columns) and emits one row per query:
the best-scoring subject, its support, the best subject from a DIFFERENT reference, and the margin
between them.

WHY THE MARGIN IS REPORTED RATHER THAN HIDDEN
---------------------------------------------
A typing database is dense on purpose — near-identical neighbours are the point of it. So the top
two subjects are often separated by very little, and a bare argmax turns a coin-flip into a
confident-looking answer. Every call therefore carries the runner-up and the margin, and anything
inside --tie-margin is labelled `ambiguous` rather than silently committed.

This is the same failure the recruitment panel work ran into from the other direction, where a
per-read argmax across near-neighbours shattered a sample across references. It is cheap to avoid
here: the query is a whole consensus rather than a 150 bp read, so when the margin really is thin
it is thin for a reason worth reporting.

MULTI-HSP HANDLING
------------------
BLAST splits a long query against a long subject into several HSPs. Scoring only the best HSP
under-reports both coverage and support for exactly the subjects that match best over their full
length. So hits are aggregated per (query, subject): bitscores are summed and aligned length is
summed over the union of query intervals, so overlapping HSPs are not double-counted. Percent
identity is reported as the aligned-length-weighted mean.
"""
import argparse
import sys
from collections import defaultdict

DEFAULT_COLUMNS = "qseqid sseqid pident length qstart qend qlen slen bitscore"

ap = argparse.ArgumentParser()
ap.add_argument("--hits", required=True, help="blastn -outfmt 6 output")
ap.add_argument("--sample", required=True, help="Sample id, emitted as the first column")
ap.add_argument("-o", "--out", required=True, help="Output TSV")
ap.add_argument("--columns", default=DEFAULT_COLUMNS,
                help="Space-separated outfmt 6 field names, in file order. Must include qseqid, "
                     "sseqid, pident, length, qstart, qend, qlen and bitscore. "
                     f"Default: {DEFAULT_COLUMNS!r}")
ap.add_argument("--min-pident", type=float, default=0.0,
                help="Drop HSPs below this percent identity before aggregating.")
ap.add_argument("--min-qcov", type=float, default=0.0,
                help="Drop a subject whose aggregated query coverage (0-1) is below this.")
ap.add_argument("--tie-margin", type=float, default=0.02,
                help="Relative bitscore margin below which a call is 'ambiguous'. Computed as "
                     "(best - runnerup)/best against the best subject from a DIFFERENT reference.")
a = ap.parse_args()

cols = a.columns.split()
required = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "qlen", "bitscore"]
missing = [c for c in required if c not in cols]
if missing:
    sys.exit(f"{sys.argv[0]} ERROR: --columns is missing required field(s): {', '.join(missing)}\n")
idx = {c: i for i, c in enumerate(cols)}

# (query, subject) -> {bits, intervals, pident_weighted, qlen}
agg = defaultdict(lambda: {"bits": 0.0, "iv": [], "pid_w": 0.0, "alen": 0, "qlen": 0})
queries = []
seen = set()

with open(a.hits) as fh:
    for line_no, line in enumerate(fh, 1):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < len(cols):
            sys.exit(f"{sys.argv[0]} ERROR: {a.hits}:{line_no} has {len(f)} fields, "
                     f"--columns declares {len(cols)}\n")
        try:
            q, s = f[idx["qseqid"]], f[idx["sseqid"]]
            pid = float(f[idx["pident"]])
            alen = int(f[idx["length"]])
            qs, qe = int(f[idx["qstart"]]), int(f[idx["qend"]])
            qlen = int(f[idx["qlen"]])
            bits = float(f[idx["bitscore"]])
        except ValueError as exc:
            sys.exit(f"{sys.argv[0]} ERROR: {a.hits}:{line_no} unparsable — {exc}\n")
        if pid < a.min_pident:
            continue
        if q not in seen:
            seen.add(q)
            queries.append(q)
        rec = agg[(q, s)]
        rec["bits"] += bits
        rec["iv"].append((min(qs, qe), max(qs, qe)))  # blastn reverses coords on minus strand
        rec["pid_w"] += pid * alen
        rec["alen"] += alen
        rec["qlen"] = qlen


def covered(intervals):
    """Union length of query intervals, so overlapping HSPs are not counted twice."""
    if not intervals:
        return 0
    merged = []
    for lo, hi in sorted(intervals):
        if merged and lo <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return sum(hi - lo + 1 for lo, hi in merged)


by_query = defaultdict(list)
for (q, s), rec in agg.items():
    qcov = covered(rec["iv"]) / rec["qlen"] if rec["qlen"] else 0.0
    if qcov < a.min_qcov:
        continue
    pid = rec["pid_w"] / rec["alen"] if rec["alen"] else 0.0
    by_query[q].append({"subject": s, "bits": rec["bits"], "qcov": qcov, "pident": pid})

with open(a.out, "w") as out:
    out.write("sample\tquery\tsubject\tpident\tqcov\tbitscore\t"
              "runnerup\trunnerup_bitscore\tmargin\tcall\n")
    for q in queries:
        hits = sorted(by_query.get(q, []), key=lambda h: (-h["bits"], h["subject"]))
        if not hits:
            out.write(f"{a.sample}\t{q}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tno_hit\n")
            continue
        best = hits[0]
        runner = hits[1] if len(hits) > 1 else None
        if runner is None:
            margin, call, r_sub, r_bits = 1.0, "confident", "NA", "NA"
        else:
            margin = (best["bits"] - runner["bits"]) / best["bits"] if best["bits"] else 0.0
            call = "ambiguous" if margin < a.tie_margin else "confident"
            r_sub, r_bits = runner["subject"], f"{runner['bits']:.1f}"
        m = f"{margin:.4f}" if runner is not None else "NA"
        out.write(f"{a.sample}\t{q}\t{best['subject']}\t{best['pident']:.2f}\t{best['qcov']:.4f}\t"
                  f"{best['bits']:.1f}\t{r_sub}\t{r_bits}\t{m}\t{call}\n")

# A query with no surviving hit is a result, not an error — an unrelated or badly assembled consensus
# should produce a `no_hit` row rather than an empty file that reads as "nothing ran".
if not queries:
    sys.stderr.write(f"{sys.argv[0]} WARNING: no hits in {a.hits}; wrote header only\n")

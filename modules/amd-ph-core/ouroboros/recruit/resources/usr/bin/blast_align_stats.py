#!/usr/bin/env python3
"""blast_align_stats.py — blastn ALIGN-path stats (parity with sam_align_stats.py).

Reads blastn outfmt-6 with the aligned sequences and builds the per-position pileup JSON
that combine_align_stats.py consumes, for building the refined reference.

Input (stdin or file), outfmt: 6 sstart send sstrand qseq sseq
  (one best HSP per read; run blastn with -max_target_seqs 1 -max_hsps 1)

JSON out: {ref_len, position_counts{pos:{base:cnt}}, leader_counts, trailer_counts}

Per aligned column (subject-forward): subject base present + query base -> tally query base
at that ref position; query gap -> '-' (deletion); subject gap -> insertion (skipped for
position_counts, like the CIGAR-I case).

NOTE: leader/trailer EXTENSION (elongation) needs the query overhang sequence, which blast
outfmt does not provide; this script supports --skip-elongation (the SKIP_E default) and
emits empty leader/trailer otherwise.
"""

import argparse
import json
import sys

COMP = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")


def revcomp(s):
    return s.translate(COMP)[::-1]


def ref_len_of(path):
    n = 0
    started = False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if started:
                    break
                started = True
            elif started:
                n += len(line.strip())
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("blast", nargs="?", default="-", help="blastn outfmt-6 (sstart send sstrand qseq sseq), or '-'")
    ap.add_argument("--ref", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-S", "--skip-elongation", action="store_true")
    ap.add_argument("-G", "--ignore-annotation", action="store_true")
    a = ap.parse_args()
    if not a.skip_elongation:
        sys.stderr.write(
            "blast_align_stats.py: elongation not supported (needs query overhang); emitting position_counts only\n"
        )

    ref_len = ref_len_of(a.ref)
    pos_counts = {}  # int ref pos (0-based) -> {base: count}

    fh = sys.stdin if a.blast == "-" else open(a.blast)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        sstart, send, sstrand, qseq, sseq = int(f[0]), int(f[1]), f[2], f[3].upper(), f[4].upper()
        # normalize to subject-forward orientation
        if sstrand == "minus" or sstart > send:
            qseq, sseq = revcomp(qseq), revcomp(sseq)
            start = min(sstart, send)
        else:
            start = sstart
        rpos = start - 1  # 0-based
        for q, s in zip(qseq, sseq):
            if s != "-":  # reference has a base here
                base = q if q != "-" else "-"  # query gap = deletion
                d = pos_counts.setdefault(rpos, {})
                d[base] = d.get(base, 0) + 1
                rpos += 1
            # s == '-' : insertion in query -> not tallied to a ref position
    if fh is not sys.stdin:
        fh.close()

    data = {
        "ref_len": ref_len,
        "position_counts": {str(p): b for p, b in pos_counts.items()},
        "leader_counts": {},
        "trailer_counts": {},
    }
    with open(a.output, "w") as out:
        json.dump(data, out, separators=(",", ":"))


if __name__ == "__main__":
    main()

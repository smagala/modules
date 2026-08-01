#!/usr/bin/env python3
"""Classify blastn HSPs into MATCH outputs for the gather loop.

Reads blastn outfmt-6 HSPs + the query FASTA; writes, for the gather convention:
  <out>.match    FASTA of recruited (clean, non-chimera) reads
  <out>.class    recruited read header \t best-bitscore target  (drives SORT)
  <out>.chim     FASTA of chimeric reads (dropped unless --incl-chim)
  <out>.nomatch  FASTA of reads with no significant HSP

A read is chimeric when its significant HSPs hit the same target on BOTH strands.
Significance = HSP length >= --min-len and %id >= --min-pid.
"""
import argparse, sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--blast", required=True, help="blastn outfmt 6: qseqid sseqid pident length qstart qend sstart send sstrand bitscore")
ap.add_argument("--query", required=True, help="query FASTA")
ap.add_argument("--out", required=True, help="output basename (writes .match/.class/.chim/.nomatch)")
ap.add_argument("--min-len", type=int, default=33)
ap.add_argument("--min-pid", type=float, default=88.0)
ap.add_argument("--incl-chim", action="store_true", help="keep chimeric reads in .match (INCL_CHIM)")
a = ap.parse_args()

strands = defaultdict(lambda: defaultdict(set))   # rid -> target -> {plus,minus}
nhsp = defaultdict(int)
best = {}                                         # rid -> (bitscore, target)

with open(a.blast) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        rid, tgt, pid, length = f[0], f[1], float(f[2]), int(f[3])
        sstrand, bits = f[8], float(f[9])
        if length < a.min_len or pid < a.min_pid:
            continue
        tgt = tgt.split("{S")[0]  # strip {Sxx} subtype suffix if present
        strands[rid][tgt].add(sstrand)
        nhsp[rid] += 1
        if rid not in best or bits > best[rid][0]:
            best[rid] = (bits, tgt)

def is_chimera(rid):
    tg = strands[rid]
    return any(len(s) > 1 for s in tg.values())   # same target, both strands

# stream query FASTA, route each record
m = open(a.out + ".match", "w")
c = open(a.out + ".class", "w")
ch = open(a.out + ".chim", "w")
nm = open(a.out + ".nomatch", "w")

def flush(hdr, seq):
    if not hdr:
        return
    rid = hdr.split()[0]
    if rid not in strands:                 # no significant HSP
        nm.write(f">{hdr}\n{seq}\n")
    elif is_chimera(rid) and not a.incl_chim:
        ch.write(f">{hdr}\n{seq}\n")
    else:
        m.write(f">{hdr}\n{seq}\n")
        c.write(f"{hdr}\t{best[rid][1]}\n")

hdr, seq = "", []
with open(a.query) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush(hdr, "".join(seq))
            hdr, seq = line[1:], []
        else:
            seq.append(line)
    flush(hdr, "".join(seq))

for fp in (m, c, ch, nm):
    fp.close()

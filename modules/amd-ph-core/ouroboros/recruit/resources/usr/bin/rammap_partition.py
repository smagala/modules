#!/usr/bin/env python3
"""rammap_partition.py — split reads from a rammap PAF into match/class/nomatch/chim.

Chimera detection lives in the aligner: rammap --filter-chimera tags every record of a
chimeric read (significant hits covering one target on both strands) with ch:i:1. This only
routes on that tag plus the recruitment significance threshold — no chimera logic of its own.

Run rammap with: -x sr --secondary yes -N 7 -c --filter-chimera
Significance: block_len >= --min-len and 100*matches/block_len >= --min-pid.
PAF cols used: 1 qname, 6 target, 10 matches, 11 block_len; the ch:i:1 tag marks chimeras.

Writes <out>.match/.class/.nomatch/.chim. Chimeric reads go to .chim (unless --incl-chim),
recruited reads to .match with the best target in .class, the rest to .nomatch.
"""
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("--paf", required=True)
ap.add_argument("--query", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--min-len", type=int, default=33)
ap.add_argument("--min-pid", type=float, default=88.0)
ap.add_argument("--incl-chim", action="store_true")
a = ap.parse_args()

best = {}      # rid -> (matches, target) among significant hits
chim = set()   # rids rammap tagged ch:i:1

with open(a.paf) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        rid, tgt = f[0], f[5]
        matches, blk = int(f[9]), int(f[10])
        if "ch:i:1" in f[11:]:
            chim.add(rid)
        if blk < a.min_len or (100.0 * matches / blk) < a.min_pid:
            continue
        tgt = tgt.split("{S")[0]  # strip {Sxx} subtype suffix if present
        if rid not in best or matches > best[rid][0]:
            best[rid] = (matches, tgt)


def records(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


with open(a.out + ".match", "w") as mfh, open(a.out + ".class", "w") as cfh, \
     open(a.out + ".nomatch", "w") as nfh, open(a.out + ".chim", "w") as chfh:
    for hdr, seq in records(a.query):
        rid = hdr.split()[0]
        if rid in chim and not a.incl_chim:
            chfh.write(">" + hdr + "\n" + seq + "\n")
        elif rid in best:
            mfh.write(">" + hdr + "\n" + seq + "\n")
            cfh.write(hdr + "\t" + best[rid][1] + "\n")
        else:
            nfh.write(">" + hdr + "\n" + seq + "\n")

#!/usr/bin/env python3
"""var_call_stats.py — build the per-position pileup for variant calling from a SAM.

Walks each read's CIGAR against the reference and accumulates, per reference position: allele counts
and summed base qualities (raw ASCII, i.e. Phred+33), plus insertion counts/qualities, deletion
counts, and the set of read alignment "coordinate strings" (dots for uncovered, bases/gaps for
covered) used downstream for insertion/deletion depth and phasing. Emits a msgpack pileup that
call.py reads.

Pileup schema (msgpack dict):
  ref_name, ref_len
  counts[pos] = {allele: count}            quals[pos] = {allele: sum(Phred+33)}
  strand[pos] = {allele: [fwd, rev]}       (per-allele strand counts, for SOR strand-bias gating)
  ins_c[pos]  = {insert: count}            ins_q[pos] = {insert: mean_phred_sum}
  del_c[pos]  = {del_len: count}           aln        = {"start:body": count}  (compact coord string)
"""

import argparse
import re

import msgpack

ap = argparse.ArgumentParser()
ap.add_argument("ref")
ap.add_argument("sam")
ap.add_argument("prefix")
a = ap.parse_args()

# Reference: first FASTA record only.
ref_name = None
ref_seq_parts = []
with open(a.ref) as fh:
    seen_header = False
    for line in fh:
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if seen_header:
                break
            ref_name = line[1:]
            seen_header = True
        elif seen_header:
            ref_seq_parts.append(line)
ref_len = len("".join(ref_seq_parts))
if ref_name is None or ref_len < 1:
    raise SystemExit("No reference found.")

counts = [dict() for _ in range(ref_len)]  # per position: {allele: count}
quals = [dict() for _ in range(ref_len)]  # per position: {allele: sum(Phred+33)}
strand = [dict() for _ in range(ref_len)]  # per position: {allele: [fwd, rev]} for SOR strand bias
ins_c = {}  # {pos: {insert: count}}
ins_q = {}  # {pos: {insert: summed mean Phred}}
del_c = {}  # {pos: {del_len: count}}
aln_counts = {}  # {alignment coordinate string: count}

CIGAR = re.compile(r"(\d+)([MIDNSHP])")


def bump(d, k):
    d[k] = d.get(k, 0) + 1


with open(a.sam) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line == "@":
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        rn, pos, cigar, seq, qual = f[2], f[3], f[5], f[9].upper(), f[10]
        if rn != ref_name:
            continue
        is_rev = (int(f[1]) & 16) != 0  # SAM FLAG 0x10 — read maps to the reverse strand
        qint = [ord(c) for c in qual]  # raw Phred+33 ASCII value of each quality byte
        rpos = int(pos) - 1
        qpos = 0
        # Each read's reference-coordinate alignment is the covered body ('.' for spliced/N, base/'-'
        # for covered) bracketed by uncovered '.' before/after. We store only the COMPACT body keyed
        # "<start>:<body>" (start = leading-uncovered count, body = leading/trailing dots stripped) — a
        # lossless encoding of the old full ref_len-padded string, ~read-length not genome-length, so
        # the {pattern: count} table no longer costs O(ref_len) per distinct read (the full-depth OOM).
        start = rpos
        body = []
        for m in CIGAR.finditer(cigar):
            inc = int(m.group(1))
            op = m.group(2)
            if op == "M":
                for _ in range(inc):
                    allele = seq[qpos]
                    body.append(allele)
                    bump(counts[rpos], allele)
                    # FIXME legacy perl logic: accumulates summed Phred+33 ASCII (call.py later
                    # subtracts count*33 to recover mean quality). Store summed Phred once parity ends.
                    quals[rpos][allele] = quals[rpos].get(allele, 0) + qint[qpos]
                    sc = strand[rpos].setdefault(allele, [0, 0])
                    sc[1 if is_rev else 0] += 1
                    qpos += 1
                    rpos += 1
            elif op == "D":
                body.append("-" * inc)
                d = del_c.setdefault(rpos - 1, {})
                d[inc] = d.get(inc, 0) + 1
                for _ in range(inc):
                    bump(counts[rpos], "-")
                    rpos += 1
            elif op == "I":
                insert = seq[qpos : qpos + inc].lower()
                ic = ins_c.setdefault(rpos - 1, {})
                ic[insert] = ic.get(insert, 0) + 1
                iq = ins_q.setdefault(rpos - 1, {})
                iq[insert] = iq.get(insert, 0.0) + (sum(qint[qpos : qpos + inc]) / len(insert)) - 33
                qpos += inc
            elif op == "S":
                qpos += inc
            elif op == "N":
                body.append("." * inc)
                rpos += inc
            elif op == "H":
                continue
            else:
                raise SystemExit("Extended CIGAR (%s) not yet supported." % op)
        # compact-encode: absorb leading dots (leading-N edge) into start, drop trailing dots.
        body = "".join(body)
        lead = len(body) - len(body.lstrip("."))
        stripped = body.strip(".")
        key = f"{start + lead}:{stripped}" if stripped else "0:"  # all-dot read -> canonical
        aln_counts[key] = aln_counts.get(key, 0) + 1

table = {
    "ref_name": ref_name,
    "ref_len": ref_len,
    "counts": counts,
    "quals": quals,
    "strand": strand,
    "ins_c": ins_c,
    "ins_q": ins_q,
    "del_c": del_c,
    "aln": aln_counts,
}
with open(a.prefix + ".pileup.msgpack", "wb") as out:
    out.write(msgpack.packb(table, use_bin_type=True))

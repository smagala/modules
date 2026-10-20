#!/usr/bin/env python3
"""winnow_pileup.py — one-pass winnow + score + pileup for the assemble ALIGN merge.

Fuses winnow_sam.py + score_sam.py + sam_stats.py into a SINGLE streaming pass over the concatenated
chunk SAM (they each re-parsed it, three passes total). Per query it keeps the single best-scoring
alignment (winnow), rewrites the SAM in place to those best records, sums their alignment scores (the
recursion's ASM_SCORE, printed to stdout), and tabulates the per-position pileup JSON (base / insertion
/ strand counts) that REFINE (combine_sam_stats.py) consumes.

The winnowed SAM is byte-identical to winnow_sam.py's, ASM_SCORE identical to score_sam.py's, and the
pileup counts identical to sam_stats.py's. The pileup inner loop is vectorized with numpy (per-CIGAR-op
positions counted via bincount instead of a per-base python loop). base_counts / strand_counts key order
may differ from sam_stats.py's first-seen order, which is harmless: combine_sam_stats.py re-sorts alleles
by (-count, base), so the consensus and variant calls are unchanged.

Precondition (from winnow_sam.py): a query's records must be CONTIGUOUS — true for concatenated chunk
SAMs (each read lands in exactly one chunk). Do not feed a coordinate-sorted / query-interleaved SAM.
"""

import argparse
import json
import os
import re
import sys
import tempfile

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("ref")
ap.add_argument("sam", help="winnowed SAM OUTPUT; also the input (rewritten in place) when --in-sams absent")
ap.add_argument("out_json")
ap.add_argument(
    "--in-sams",
    nargs="+",
    default=None,
    help="read these chunk SAMs IN ORDER as one stream (no concatenation) and write the winnowed "
    "SAM to `sam`. Skips the giant intermediate concatenated SAM the MERGE used to build "
    "(one full write + read of the biggest file in the pipeline). Omit to read+rewrite `sam` "
    "in place (single-chunk / legacy). Each read is wholly within one chunk, so per-query "
    "records stay contiguous across the chunk boundary — the winnow precondition holds.",
)
ap.add_argument(
    "-F",
    "--score-field",
    type=int,
    default=None,
    help="1-based SAM column holding the winnow score tag (default 12 / AS)",
)
ap.add_argument("-G", "--ignore-annotation", action="store_true")
ap.add_argument("-S", "--silence-complex-indel", action="store_true")
ap.add_argument("-q", "--min-bq", type=int, default=0)
a = ap.parse_args()

# winnow score column (mirror winnow_sam.py: default 0 -> col 12; >=12 -> index past col 11)
SF = (
    0
    if (a.score_field is None or a.score_field < 0)
    else (a.score_field - 12 if a.score_field >= 12 else a.score_field)
)
# SAM record columns (0-based); TAGS = first optional field (col 12, where the winnow AS score lives)
QNAME, FLAG, RNAME, POS, CIGAR, SEQ, QUAL, TAGS = 0, 1, 2, 3, 5, 9, 10, 11
REVERSE_FLAG = 16  # BAM FLAG bit 0x10 = read reverse-strand

# NOTE: the minus is required. bowtie2 --end-to-end emits NEGATIVE alignment scores (AS:i:-10,
# 0 = perfect); without it the match fails and win_score() silently falls back to the CIGAR M-count,
# quietly substituting a different quantity for the convergence objective. Positive scores match
# exactly as before, so this is behaviour-preserving for local_sw/rammap/bwa.
AS_RE = re.compile(r"AS:\w:(-?\d+)")
CIG_RE = re.compile(r"(\d+)([MIDNSHP])")
ANNOT_RE = re.compile(r"^([^{]+)\{[^}]*\}")
COMPLEX_INDEL = re.compile(r"\d\d+[DI]\d+M+")
COMPLEX_RUN = re.compile(r"^\d+M(\d+[DI]\d+M){4,}$")


def count_match(cig):
    return sum(int(n) for n, op in CIG_RE.findall(cig) if op == "M")


def win_score(cigar, AS):
    m = AS_RE.match(AS) if AS else None
    return int(m.group(1)) if m else count_match(cigar)


def first_ref(path):
    with open(path) as fh:
        for rec in fh.read().split(">"):
            if rec == "":
                continue
            lines = re.split(r"\r\n|\n|\r", rec)
            seq = "".join(lines[1:])
            if seq:
                return lines[0], seq
    return None, None


ref_name, ref_seq = first_ref(a.ref)
if ref_name is None:
    sys.exit(f"{sys.argv[0]} ERROR: no reference found in {a.ref}")
N = len(ref_seq)
if a.ignore_annotation:
    m = ANNOT_RE.match(ref_name)
    if m:
        ref_name = m.group(1)

# strand-resolved dense counts [fwd, rev], each a flat [pos*256 + base_byte] vector filled via bincount.
# base_counts = fwd + rev. Insertions are sparse (rare) so stay a dict.
strand = [np.zeros(N * 256, dtype=np.int64), np.zeros(N * 256, dtype=np.int64)]
pend = [[], []]  # pending flat-index arrays per strand, flushed into `strand` periodically
pend_n = [0, 0]
ins_counts = {}
FLUSH = 4_000_000


def flush_strand(r):
    if pend[r]:
        strand[r] += np.bincount(np.concatenate(pend[r]), minlength=N * 256)
        pend[r] = []
        pend_n[r] = 0


def add_flat(r, flat):
    pend[r].append(flat)
    pend_n[r] += flat.size
    if pend_n[r] >= FLUSH:
        flush_strand(r)


def pileup(f, rev):
    rname = f[RNAME]
    if a.ignore_annotation:
        m = ANNOT_RE.match(rname)
        if m:
            rname = m.group(1)
    if rname != ref_name:
        return
    cigar = f[CIGAR]
    if a.silence_complex_indel and COMPLEX_INDEL.search(cigar) and COMPLEX_RUN.match(cigar):
        cigar = re.sub(r"(\d+)D", r"\1N", cigar)
        cigar = re.sub(r"(\d+)I", r"\1S", cigar)
    seq = f[SEQ].upper()
    seq_b = np.frombuffer(seq.encode("latin-1"), dtype=np.uint8)
    qual = f[QUAL] if len(f) > QUAL else "*"
    qual_b = np.frombuffer(qual.encode("latin-1"), dtype=np.uint8) if qual != "*" else None
    rpos = int(f[POS]) - 1
    qpos = 0
    for m in CIG_RE.finditer(cigar):
        inc = int(m.group(1))
        op = m.group(2)
        if op == "M":
            rp = np.arange(rpos, rpos + inc, dtype=np.int64)
            od = seq_b[qpos : qpos + inc].astype(np.int64)
            if a.min_bq > 0 and qual_b is not None and qpos + inc <= qual_b.size:
                keep = (qual_b[qpos : qpos + inc].astype(np.int16) - 33) >= a.min_bq
                rp = rp[keep]
                od = od[keep]
            add_flat(rev, rp * 256 + od)
            qpos += inc
            rpos += inc
        elif op == "D":
            rp = np.arange(rpos, rpos + inc, dtype=np.int64)
            add_flat(rev, rp * 256 + 45)  # 45 == ord('-')
            rpos += inc
        elif op == "I":
            insert = seq[qpos : qpos + inc].lower()
            d = ins_counts.setdefault(rpos - 1, {})
            d[insert] = d.get(insert, 0) + 1
            qpos += inc
        elif op == "N":
            rpos += inc
        elif op == "S":
            qpos += inc
        elif op == "H":
            pass
        else:
            sys.exit(f"Extended CIGAR ({op}) not supported.")


# streaming winnow: keep the best record per contiguous query; on query change flush it to the SAM,
# add its score, and pile it up. In-place: write a temp beside the input and rename over it at the end.
tmp = tempfile.NamedTemporaryFile(
    "w", delete=False, dir=os.path.dirname(os.path.abspath(a.sam)) or ".", prefix=".wp.", suffix=".sam"
)
prev_header = ""
cur_q = cur_score = cur_rec = cur_f = None
cur_rev = 0
total_as = 0


def flush_query():
    global total_as
    if cur_rec is None:
        return
    tmp.write(cur_rec + "\n")
    # BUG FIX: this summed column 12 unconditionally, ignoring --score-field. win_score() honours the
    # flag when PICKING the best record per read, but the reported total did not — so for any aligner
    # whose AS tag is not in column 12 the loop's convergence objective was whatever sat there instead.
    # For rammap/bwa/bwa-mem2/minibwa that is NM (edit distance), and because the caller selects the
    # HIGHEST-scoring iteration, the loop was preferring the iteration with the MOST mismatches.
    # cur_score is already the win_score of the record being flushed, so use it.
    total_as += cur_score
    pileup(cur_f, cur_rev)


# stream the input SAM(s) in order. With --in-sams we read the raw chunk SAMs directly (each carries its
# own header); headers are valid only at the very top, so we emit them from the first chunk and skip every
# later chunk's header once records have started — byte-identical to concatenating (header-once + bodies).
inputs = a.in_sams if a.in_sams else [a.sam]
seen_record = False
for _path in inputs:
    with open(_path) as fh:
        for line in fh:
            if line[:1] == "@":
                if not seen_record and line != prev_header:  # header block, first chunk only, de-duped
                    tmp.write(line)
                prev_header = line
                continue
            seen_record = True
            rec = line.rstrip("\n")
            f = rec.split("\t")
            if f[CIGAR] == "*":  # unmapped: winnow_sam drops these
                continue
            extra = f[TAGS:]
            score = win_score(f[CIGAR], extra[SF] if SF < len(extra) else None)
            if f[QNAME] != cur_q:
                flush_query()
                cur_q, cur_score, cur_rec, cur_f = f[QNAME], score, rec, f
                cur_rev = 1 if (int(f[FLAG]) & REVERSE_FLAG) else 0
            elif cur_score < score:  # strict: first occurrence wins ties (matches original)
                cur_score, cur_rec, cur_f = score, rec, f
                cur_rev = 1 if (int(f[FLAG]) & REVERSE_FLAG) else 0
flush_query()

tmp.close()
os.replace(tmp.name, a.sam)

flush_strand(0)
flush_strand(1)
base = strand[0] + strand[1]
base_counts, strand_counts = {}, {}
for flat in np.nonzero(base)[0]:
    flat = int(flat)
    p, o = divmod(flat, 256)
    ch = chr(o)
    base_counts.setdefault(str(p), {})[ch] = int(base[flat])
    strand_counts.setdefault(str(p), {})[ch] = [int(strand[0][flat]), int(strand[1][flat])]
data = {
    "base_counts": base_counts,
    "ins_counts": {str(p): b for p, b in ins_counts.items()},
    "strand_counts": strand_counts,
}
with open(a.out_json, "w") as out:
    json.dump(data, out, separators=(",", ":"))

sys.stdout.write(str(total_as))

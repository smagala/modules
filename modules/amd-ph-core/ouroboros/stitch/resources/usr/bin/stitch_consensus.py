#!/usr/bin/env python3
"""stitch_consensus.py — place a deep assembled consensus into the full reference coordinate frame.

The iterative assembly produces high-quality base calls but TRIMS uncovered columns, so its consensus
is shorter than the reference and loses coordinate anchoring (see project_fulllength_scaffold_policy).
The amplicon path keeps full length + N but calls on the divergent panel reference, so its calls are
shallower. This tool combines the two: it aligns the assembled consensus to the full reference and
emits a reference-length consensus carrying the assembly's calls, with the gaps resolved honestly:

  base  — a reference column the assembly aligned a base to (the deep, accurate call).
  base  — recovered: at an assembly-GAP column that the reads cover with a BASE (deletion fraction
          < --del-frac), the assembly gap is an artifact; emit the read-majority base instead of dropping it.
  N     — missing: a reference column the assembly never spanned, OR an assembly-gap column too shallow
          (< --min-depth) to adjudicate.
  -     — deletion: an assembly-gap column the READS agree is deleted (deletion fraction >= --del-frac at
          >= --min-depth spanning reads). Stripped on output, so the consensus shortens there.

The BAM (--bam, reads mapped to the same reference) is what adjudicates an assembly gap: reads carrying a
deletion => real deletion (drop); reads carrying a base => assembly artifact (recover the base); too shallow
=> N. Without a BAM, every assembly gap is conservatively N. This prevents the stitch from propagating a
spurious assembly-deletion (e.g. a thin-coverage gap the reads don't actually support).

Length model matches amend_consensus.py: N is kept (length-preserving), '-' is stripped (shortens),
insertions in the assembly relative to the reference are dropped (the reference frame has no column
for them; they live in the VCF / insertion table).
"""
import argparse
import subprocess
import sys
import tempfile

import pysam


def read_fasta_record(path, name=None):
    """Return (header, sequence) for the named record, or the first record if name is None."""
    hdr, parts, capturing = None, [], False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if capturing:
                    break
                rec = line[1:].split()[0]
                if name is None or rec == name:
                    hdr, capturing = line[1:], True
            elif capturing:
                parts.append(line)
    return hdr, "".join(parts).upper()


def align_consensus_sam(consensus_path, ref_path, preset):
    """Align consensus (query) to reference (target) with rammap; return path to a temp SAM file.

    rammap is the pure-Rust, minimap2-compatible aligner used throughout the pipeline (same
    -a/-x preset/--secondary flags); there is no minimap2 dependency here.
    """
    sam = tempfile.NamedTemporaryFile(suffix=".sam", delete=False)
    sam.close()
    with open(sam.name, "w") as out:
        subprocess.run(
            ["rammap", "-a", "-x", preset, "--secondary", "no", ref_path, consensus_path],
            stdout=out, stderr=subprocess.DEVNULL, check=True,
        )
    return sam.name


_BASES = "ACGT"


def read_evidence(bam_path, ref_name, ref_len):
    """Per reference position, return (base_counts, del_depth):
      base_counts[p] = (A, C, G, T) read counts (reads with a base aligned at p),
      del_depth[p]   = reads that SPAN p carrying a deletion/ref-skip there.
    Returns (None, None) if no BAM. del_depth distinguishes a true deletion (reads carry a gap) from
    an assembly-gap artifact (reads carry a base)."""
    if not bam_path:
        return None, None
    base_counts = [(0, 0, 0, 0)] * ref_len
    del_depth = [0] * ref_len
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        cov = bam.count_coverage(ref_name, 0, ref_len, quality_threshold=0)
        base_counts = [(cov[0][p], cov[1][p], cov[2][p], cov[3][p]) for p in range(ref_len)]
        for col in bam.pileup(ref_name, 0, ref_len, truncate=True, min_base_quality=0):
            d = sum(1 for pr in col.pileups if pr.is_del or pr.is_refskip)
            del_depth[col.reference_pos] = d
    return base_counts, del_depth


def majority_base(counts):
    """Dominant base from an (A,C,G,T) count tuple; 'N' if empty or tied."""
    m = max(counts)
    if m == 0 or counts.count(m) > 1:
        return "N"
    return _BASES[counts.index(m)]


def project_cigar(cigartuples, qseq, ref_start, ref_len, frame, max_clip=10):
    """Project one alignment's query bases onto the reference frame (mutates and returns `frame`).

    frame[p] stays None until first written (first alignment wins). Operations:
      M/=/X  place the query base at reference column p.
      D/N    mark an assembly gap ('-') — the consensus has no base at that reference column.
      I      insertion: query advances with no reference column (dropped; it lives in the VCF).
      S      soft-clip: a TERMINAL query run the aligner declined to align. A SHORT clip (<= max_clip)
             is a clipped terminal MISMATCH — rammap clips the last base(s) instead of emitting an X, so
             dropping it to N would blank a base the amend step already called (and IRMA keeps, e.g. an
             amended 3' base differing from the panel-reference base). Project it as an ungapped extension
             onto the flanking reference columns (leading clip leftward of ref_start, trailing clip
             rightward of the alignment end). A LONG clip is STRUCTURAL — the aligner could not place a
             divergent/rearranged tail contiguously, so an ungapped projection would misplace those bases;
             leave those reference columns N (the honest, coordinate-preserving default). All-or-nothing
             per clip. Out-of-range and already-filled columns are always skipped.
      H/P    consume neither query nor reference.

    query_sequence (qseq) includes soft-clipped bases but not hard-clipped ones, matching this indexing.
    """
    ref_pos = ref_start
    q = 0
    for op, length in cigartuples:
        if op in (0, 7, 8):            # M/=/X: aligned bases
            for _ in range(length):
                if 0 <= ref_pos < ref_len and frame[ref_pos] is None:
                    frame[ref_pos] = qseq[q]
                ref_pos += 1
                q += 1
        elif op == 2 or op == 3:       # D / N (ref skip): assembly has no base here
            for _ in range(length):
                if 0 <= ref_pos < ref_len and frame[ref_pos] is None:
                    frame[ref_pos] = "-"
                ref_pos += 1
        elif op == 4:                  # soft-clip: recover a short terminal-mismatch clip; drop long ones
            if length <= max_clip:
                if q == 0:             # leading clip -> ref[ref_pos-length : ref_pos]
                    for k in range(length):
                        p = ref_pos - (length - k)
                        if 0 <= p < ref_len and frame[p] is None:
                            frame[p] = qseq[q + k]
                else:                  # trailing clip -> ref[ref_pos : ref_pos+length]
                    for k in range(length):
                        p = ref_pos + k
                        if 0 <= p < ref_len and frame[p] is None:
                            frame[p] = qseq[q + k]
            q += length
        elif op == 1:                  # insertion: query advances, no reference column
            q += length
        # H (5) / P (6): no advance
    return frame


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference", help="full reference FASTA (the coordinate frame)")
    ap.add_argument("consensus", help="assembled consensus FASTA (the deep calls)")
    ap.add_argument("--ref-name", default=None, help="reference record to use (default: first)")
    ap.add_argument("--bam", default=None,
                    help="reads mapped to the SAME reference; depth splits real deletion from dropout")
    ap.add_argument("--min-depth", type=int, default=10,
                    help="minimum spanning depth to adjudicate an assembly-gap column (below = N; default 10)")
    ap.add_argument("--del-frac", type=float, default=0.5,
                    help="read deletion fraction at/above which an assembly-gap column is a REAL deletion "
                         "(below = the gap is an assembly artifact; recover the read-majority base; default 0.5)")
    ap.add_argument("--preset", default="asm20", help="rammap -x preset (default asm20)")
    ap.add_argument("--max-clip-extend", type=int, default=10,
                    help="max terminal soft-clip length to recover as an ungapped extension (a clipped "
                         "terminal mismatch); longer clips are structural and left as N (default 10)")
    ap.add_argument("-N", "--name", default=None, help="output header/name")
    ap.add_argument("-o", "--out", default=None, help="output FASTA (default stdout)")
    a = ap.parse_args()

    ref_hdr, ref_seq = read_fasta_record(a.reference, a.ref_name)
    if not ref_seq:
        sys.exit(f"{sys.argv[0]} ERROR: no reference sequence in {a.reference}")
    ref_len = len(ref_seq)
    ref_name = ref_hdr.split()[0]

    base_counts, del_depth = read_evidence(a.bam, ref_name, ref_len)

    # frame[p]: None = never spanned (missing), "-" = assembly deletion at p, else the assembly base.
    frame = [None] * ref_len
    sam_path = align_consensus_sam(a.consensus, a.reference, a.preset)
    with pysam.AlignmentFile(sam_path, "r") as sam:
        for aln in sam:
            if aln.is_unmapped or aln.is_secondary:
                continue
            qseq = aln.query_sequence
            if qseq is None:
                continue
            project_cigar(aln.cigartuples, qseq, aln.reference_start, ref_len, frame,
                          max_clip=a.max_clip_extend)

    # resolve each reference column to a single character
    out_chars = []
    n_base = n_del = n_missing = n_recovered = 0
    for p in range(ref_len):
        c = frame[p]
        if c is None:                          # never spanned by the assembly -> missing
            out_chars.append("N"); n_missing += 1
        elif c == "-":                         # assembly gap: trust it only if the READS show a deletion
            if base_counts is None:            # no BAM -> can't adjudicate; keep length-honest as N
                out_chars.append("N"); n_missing += 1
                continue
            bdepth = sum(base_counts[p]); ddepth = del_depth[p]
            spanning = bdepth + ddepth
            if spanning >= a.min_depth and ddepth >= a.del_frac * spanning:
                n_del += 1                     # reads agree it's deleted -> drop the column (shortens)
            elif bdepth >= a.min_depth:
                # covered, but the reads carry a BASE not a deletion -> assembly-gap artifact.
                # recover the read-majority base instead of dropping it.
                out_chars.append(majority_base(base_counts[p])); n_recovered += 1
            else:
                out_chars.append("N"); n_missing += 1   # too shallow to call -> dropout
        else:
            out_chars.append(c); n_base += 1

    name = a.name or ref_name
    seq = "".join(out_chars)
    sys.stderr.write(
        f"stitch: ref_len={ref_len} called={n_base} recovered={n_recovered} N={n_missing} "
        f"deletions={n_del} out_len={len(seq)}\n"
    )
    record = f">{name}\n{seq}\n"
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(record)
    else:
        sys.stdout.write(record)


if __name__ == "__main__":
    main()

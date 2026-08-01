#!/usr/bin/env python3
"""sam_align_stats.py

Parse a SAM alignment stream into per-position alignment statistics with
leader/trailer extension data for iterative reference refinement (the ALIGN
path of the gather loop). Format-only: it does not care which aligner produced
the SAM.

Output: JSON file with per-position base counts plus 5'/3' overhang
counts that feed into combine_align_stats.py for consensus building
with reference extension.

Requires: pysam
"""

import argparse
import json
import sys

import pysam

MAX_GAP_REPAIR = 9


def read_reference(path):
    """Read first sequence from FASTA file. Return (name, seq, length)."""
    name = None
    seqs = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)
    seq = "".join(seqs).upper()
    return name, seq, len(seq)


def process_sam(sam_path, ref_len, elongate, term_window=0):
    """Process SAM/BAM records via pysam.

    Args:
        sam_path: path to SAM file, or "-" for stdin
        ref_len: reference sequence length
        elongate: whether to track leader/trailer overhangs

    Returns aggregated stats dicts:
        position_counts: {int_pos: {base: count}}
        leader_counts:   {str_neg_pos: {base: count}}
        trailer_counts:  {str_pos: {base: count}}
    """
    position_counts = {}
    leader_counts = {}
    trailer_counts = {}
    # position_terminal[pos][base] = count of aligned bases that sit within term_window of a read
    # terminus (the full read, soft-clips included). Used downstream to prefer spanning evidence over
    # terminal artifacts when amending the reference. Empty (and skipped) when term_window == 0.
    position_terminal = {}

    mode = "r"  # SAM text mode
    save = pysam.set_verbosity(0)  # suppress missing index warnings
    with pysam.AlignmentFile(sam_path, mode) as samfile:
        for read in samfile.fetch(until_eof=True):
            if read.is_unmapped or read.cigartuples is None:
                continue

            seq = read.query_sequence.upper()
            rpos = read.reference_start  # 0-based
            qpos = 0
            cigar = read.cigartuples
            n_ops = len(cigar)

            for idx, (op, length) in enumerate(cigar):
                if op == 4:  # S = soft clip
                    if idx == 0 and elongate:
                        # Leading soft-clip — potential 5' leader
                        clip_bases = seq[qpos : qpos + length].lower()
                        gap = rpos  # distance from ref start

                        if gap == 0:
                            # Alignment starts at ref position 0
                            for x in range(-length, 0):
                                base = clip_bases[length + x]
                                key = str(x)
                                d = leader_counts.setdefault(key, {})
                                d[base] = d.get(base, 0) + 1
                        elif 0 < gap <= MAX_GAP_REPAIR and len(clip_bases) >= gap:
                            # Gap repair: fill leading gap from clip end
                            for g in range(gap):
                                fill_base = clip_bases[
                                    len(clip_bases) - gap + g
                                ].upper()
                                d = position_counts.setdefault(g, {})
                                d[fill_base] = d.get(fill_base, 0) + 1
                            # Remaining clip → leader
                            remaining_len = len(clip_bases) - gap
                            if remaining_len > 0:
                                remaining = clip_bases[:remaining_len]
                                for x in range(-remaining_len, 0):
                                    base = remaining[remaining_len + x]
                                    key = str(x)
                                    d = leader_counts.setdefault(key, {})
                                    d[base] = d.get(base, 0) + 1

                    elif idx == n_ops - 1 and elongate:
                        # Trailing soft-clip — potential 3' trailer
                        clip_bases = seq[qpos : qpos + length].lower()
                        gap_from_end = ref_len - rpos

                        if gap_from_end == 0:
                            for x in range(length):
                                base = clip_bases[x]
                                key = str(x)
                                d = trailer_counts.setdefault(key, {})
                                d[base] = d.get(base, 0) + 1
                        elif (
                            0 < gap_from_end <= MAX_GAP_REPAIR
                            and len(clip_bases) >= gap_from_end
                        ):
                            # Gap repair at trailing end
                            for g in range(gap_from_end):
                                fill_base = clip_bases[g].upper()
                                fill_pos = rpos + g
                                d = position_counts.setdefault(fill_pos, {})
                                d[fill_base] = d.get(fill_base, 0) + 1
                            remaining = clip_bases[gap_from_end:]
                            for x in range(len(remaining)):
                                base = remaining[x]
                                key = str(x)
                                d = trailer_counts.setdefault(key, {})
                                d[base] = d.get(base, 0) + 1

                    qpos += length

                elif op in (0, 7, 8):  # M, =, X
                    L = len(seq)
                    for _ in range(length):
                        base = seq[qpos]
                        d = position_counts.setdefault(rpos, {})
                        d[base] = d.get(base, 0) + 1
                        if term_window and min(qpos, L - 1 - qpos) < term_window:
                            dt = position_terminal.setdefault(rpos, {})
                            dt[base] = dt.get(base, 0) + 1
                        qpos += 1
                        rpos += 1

                elif op == 2:  # D = deletion
                    for _ in range(length):
                        d = position_counts.setdefault(rpos, {})
                        d["-"] = d.get("-", 0) + 1
                        rpos += 1

                elif op == 1:  # I = insertion
                    qpos += length

                elif op == 3:  # N = ref skip
                    rpos += length

                elif op == 5:  # H = hard clip
                    pass

    pysam.set_verbosity(save)
    return position_counts, leader_counts, trailer_counts, position_terminal


def stats_to_json(position_counts, leader_counts, trailer_counts, ref_len, position_terminal=None):
    """Convert stats dicts to JSON-serializable structure."""
    pc_json = {}
    for pos in range(ref_len):
        if pos in position_counts:
            pc_json[str(pos)] = position_counts[pos]
    # Include any positions beyond ref_len (from gap repair)
    for pos in position_counts:
        if pos >= ref_len:
            pc_json[str(pos)] = position_counts[pos]

    out = {
        "format": "align_stats_v1",
        "ref_len": ref_len,
        "position_counts": pc_json,
        "leader_counts": leader_counts,
        "trailer_counts": trailer_counts,
    }
    # Only present when terminal tracking was requested (--term-window > 0). Keeps the default JSON
    # byte-for-byte unchanged so parity runs are unaffected.
    if position_terminal:
        out["position_terminal"] = {str(p): position_terminal[p] for p in position_terminal}
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Parse a SAM stream into alignment stats JSON"
    )
    parser.add_argument(
        "--ref", required=True, help="Reference FASTA (for length and gap repair)"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output JSON stats file"
    )
    parser.add_argument(
        "-S",
        "--skip-elongation",
        action="store_true",
        help="Skip leader/trailer extension tracking",
    )
    parser.add_argument(
        "-G",
        "--ignore-annotation",
        action="store_true",
        help="Ignore {annotation} in reference names",
    )
    parser.add_argument(
        "--term-window",
        type=int,
        default=0,
        help="Tally per-position terminal base counts (aligned bases within N bp of a read end) so "
        "the amendment step can prefer spanning evidence. 0 = off (default; JSON unchanged).",
    )
    args = parser.parse_args()

    _, _, ref_len = read_reference(args.ref)
    elongate = not args.skip_elongation

    # pysam reads from file path; if piped, use stdin via "-"
    sam_path = "-"

    position_counts, leader_counts, trailer_counts, position_terminal = process_sam(
        sam_path, ref_len, elongate, args.term_window
    )

    data = stats_to_json(
        position_counts, leader_counts, trailer_counts, ref_len, position_terminal
    )

    with open(args.output, "w") as fh:
        json.dump(data, fh, separators=(",", ":"))


if __name__ == "__main__":
    main()

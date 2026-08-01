#!/usr/bin/env python3
"""combine_align_stats.py

Aggregate alignment statistics JSON files and build consensus reference
with optional 5'/3' leader/trailer extension (the open ALIGN path of the
gather loop). Reads JSON produced by sam_align_stats.py /
blast_align_stats.py. Outputs consensus FASTA to stdout.
"""

import argparse
import json
import sys

LONG_MAX = 2**63 - 1


def read_reference(path):
    """Read first sequence from FASTA. Return (seq, length) or (None, 0)."""
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
    seq = "".join(seqs)
    return seq, len(seq)


def load_and_aggregate(stat_files, elongate):
    """Load JSON stats files and aggregate counts.

    Returns:
        count: list of dicts [{base: count}, ...] indexed by ref position
        count5: dict of dicts {str_neg_pos: {base: count}} for leaders
        count3: dict of dicts {str_pos: {base: count}} for trailers
        ref_len: reference length from first file
        count_term: list of dicts [{base: terminal_count}, ...] parallel to count; empty when the
                    upstream stats carried no "position_terminal" (i.e. --term-window was 0).
    """
    count = []
    count5 = {}
    count3 = {}
    count_term = []
    ref_len = 0

    for path in stat_files:
        with open(path) as fh:
            data = json.load(fh)

        file_ref_len = data.get("ref_len", 0)
        if file_ref_len > ref_len:
            ref_len = file_ref_len

        # Extend count arrays if needed
        while len(count) < ref_len:
            count.append({})
        while len(count_term) < ref_len:
            count_term.append({})

        # Aggregate position counts
        for pos_str, bases in data.get("position_counts", {}).items():
            pos = int(pos_str)
            if pos >= len(count):
                while len(count) <= pos:
                    count.append({})
                while len(count_term) <= pos:
                    count_term.append({})
            for base, cnt in bases.items():
                count[pos][base] = count[pos].get(base, 0) + cnt

        # Aggregate per-position terminal counts (present only with --term-window > 0 upstream)
        for pos_str, bases in data.get("position_terminal", {}).items():
            pos = int(pos_str)
            while len(count_term) <= pos:
                count_term.append({})
            for base, cnt in bases.items():
                count_term[pos][base] = count_term[pos].get(base, 0) + cnt

        if elongate:
            # Aggregate leader counts
            for pos_str, bases in data.get("leader_counts", {}).items():
                count5.setdefault(pos_str, {})
                for base, cnt in bases.items():
                    count5[pos_str][base] = count5[pos_str].get(base, 0) + cnt

            # Aggregate trailer counts
            for pos_str, bases in data.get("trailer_counts", {}).items():
                count3.setdefault(pos_str, {})
                for base, cnt in bases.items():
                    count3[pos_str][base] = count3[pos_str].get(base, 0) + cnt

    return count, count5, count3, ref_len, count_term


def majority_base(base_counts):
    """Return (max_base, max_count, alt_base, alt_count, total) from counts.

    Ignores '-' (deletion) when computing majority non-gap base.
    """
    max_count = 0
    max_base = ""
    alt_count = 0
    alt_base = ""
    total = 0

    for base, cnt in base_counts.items():
        if base == "-":
            continue
        total += cnt
        if cnt > max_count:
            alt_count = max_count
            alt_base = max_base
            max_count = cnt
            max_base = base
        elif cnt > alt_count:
            alt_count = cnt
            alt_base = base

    return max_base, max_count, alt_base, alt_count, total


def spanning_majority(base_counts, term_counts):
    """Plurality base among SPANNING (non-terminal) support: (base, spanning_count).

    Spanning count for a base = total count minus its terminal count. Ignores '-' (deletion).
    Returns ("", 0) when there is no spanning support.
    """
    best_base = ""
    best_cnt = 0
    for base, cnt in base_counts.items():
        if base == "-":
            continue
        sp = cnt - term_counts.get(base, 0)
        if sp > best_cnt:
            best_cnt = sp
            best_base = base
    return best_base, best_cnt


def build_consensus(
    count,
    count5,
    count3,
    ref_len,
    elongate,
    name,
    alt_count_thresh,
    alt_freq_thresh,
    delete_by_ambig,
    keep_deleted,
    ref_sites,
    min_pad_count,
    denominator,
    debug,
    count_term=None,
    pb_bias_frac=0.0,
    pb_min_span=5,
    min_ref_depth=0,
    extend_terminus=False,
):
    """Build consensus FASTA with optional leader/trailer extension.

    Position-bias guard (when pb_bias_frac > 0 and count_term is provided): at a position whose
    plurality base is >= pb_bias_frac terminal-supported while a DIFFERENT base holds the spanning
    plurality (>= pb_min_span spanning reads), amend to the spanning base instead. This keeps a
    force-aligned terminal artifact (e.g. a RACE non-templated 5' base) from flipping the reference,
    which the exhaustive assembler would then lock in. No-op on clean data (spanning == plurality).
    """
    # Never build consensus PAST the reference terminus unless explicitly extending. The aligner can
    # over-reach a few reads past ref_len (rammap -a records them in position_counts beyond ref_len);
    # load_and_aggregate then grows `count` to fit them. On a complete reference that tail is
    # off-the-end artifact (adapter / low-quality readthrough the aligner should have clipped), so cap
    # the build at ref_len. --extend-terminus opts back into the full range (over-reach + leader/trailer)
    # for genuine short/partial references (the reach-gap case). See findings 2026-07-13 terminal-reach.
    n_count = len(count) if extend_terminus else min(len(count), ref_len)

    if keep_deleted and ref_sites is not None and n_count != len(ref_sites):
        print(
            f"WARNING (combine_align_stats): {n_count} != {len(ref_sites)}, "
            f"bad reference, turning off keep-deleted.",
            file=sys.stderr,
        )
        keep_deleted = False
        ref_sites = None

    # Debug output
    if debug:
        # Leader positions
        for p in sorted(count5.keys(), key=lambda x: int(x)):
            bases_sorted = sorted(
                count5[p].items(), key=lambda x: (-x[1], x[0])
            )
            parts = "\t".join(f"{b}:{c}" for b, c in bases_sorted)
            print(f"{int(p):5d}::\t{parts}", file=sys.stderr)

        print("5'", file=sys.stderr)
        for p in range(n_count):
            if not count[p]:
                continue
            bases_sorted = sorted(
                count[p].items(), key=lambda x: (-x[1], x[0])
            )
            parts = "\t".join(f"{b}:{c}" for b, c in bases_sorted)
            print(f"{p:5d}::\t{parts}", file=sys.stderr)

        print("3'", file=sys.stderr)
        for p in sorted(count3.keys(), key=lambda x: int(x)):
            bases_sorted = sorted(
                count3[p].items(), key=lambda x: (-x[1], x[0])
            )
            parts = "\t".join(f"{b}:{c}" for b, c in bases_sorted)
            print(f"+{int(p):4d}::\t{parts}", file=sys.stderr)

    # Build middle consensus
    seq = []
    seq_alt = []
    max5 = 0
    max3 = 0

    pb_on = pb_bias_frac > 0 and count_term is not None
    for j in range(n_count):
        max_base, max_cnt, alt_base, alt_cnt, total = majority_base(count[j])

        # Min-depth guard (min_tcc for the reference): mutating the reference is tantamount to calling
        # the variant, so don't set/flip a reference base on fewer reads than we'd need to CALL one.
        # Below the floor, keep the prior reference base (ref_sites, staged via del_type=REF) rather
        # than a thin plurality. No-op where total >= floor or no ref available (parity on deep WGS).
        if min_ref_depth > 0 and total < min_ref_depth and ref_sites is not None and j < len(ref_sites):
            rb = ref_sites[j]
            seq.append(rb); seq_alt.append(rb)
            if j == 0: max5 = max_cnt
            if j == n_count - 1: max3 = max_cnt
            continue

        # Position-bias guard: prefer the spanning consensus when the plurality base is a
        # terminal-dominated artifact and spanning reads agree on a different base.
        amended = max_base
        if pb_on and max_base and j < len(count_term):
            sp_base, sp_cnt = spanning_majority(count[j], count_term[j])
            if (
                sp_base
                and sp_base != max_base
                and sp_cnt >= pb_min_span
                and max_cnt > 0
                and count_term[j].get(max_base, 0) / max_cnt >= pb_bias_frac
            ):
                amended = sp_base
                if debug:
                    print(
                        f"[position-bias] pos {j}: {max_base} is "
                        f"{100*count_term[j].get(max_base,0)/max_cnt:.0f}% terminal; "
                        f"amending to spanning {sp_base} ({sp_cnt})",
                        file=sys.stderr,
                    )

        if max_base:
            seq.append(amended)
            if amended != max_base:
                seq_alt.append(amended)
            elif (
                not alt_base
                or alt_cnt < alt_count_thresh
                or (total > 0 and (alt_cnt / total) < alt_freq_thresh)
            ):
                seq_alt.append(max_base)
            else:
                seq_alt.append(alt_base)
        elif delete_by_ambig:
            seq.append("N")
            seq_alt.append("N")
        elif keep_deleted and ref_sites is not None:
            seq.append(ref_sites[j])
            seq_alt.append(ref_sites[j])
        # else: position deleted (omitted)

        if j == 0:
            max5 = max_cnt
        if j == n_count - 1:
            max3 = max_cnt

    seq_str = "".join(seq)
    seq_alt_str = "".join(seq_alt)

    # Headers
    if name:
        header = f">{name}\n"
        header2 = f">{name}{{alt}}\n"
    else:
        header = ">consensus\n"
        header2 = ">alternative\n"

    if extend_terminus:
        # LEADER / TRAILER extension — reach PAST the reference ends (off by default; see --extend-terminus).
        # LEADER extension
        leader = []
        threshold = max(min_pad_count, max5 // denominator + 1)
        for p in range(1, len(count5) + 1):
            key = str(-p)
            if key not in count5:
                break
            bases = count5[key]
            max_cnt = 0
            max_base = ""
            for base, cnt in bases.items():
                if cnt > max_cnt:
                    max_cnt = cnt
                    max_base = base

            if max_cnt < threshold or max_base.upper() not in "ACTG":
                break
            leader.append(max_base)
            threshold = max(min_pad_count, max_cnt // denominator + 1)
        leader_str = "".join(reversed(leader))

        # TRAILER extension
        trailer = []
        threshold = max(min_pad_count, max3 // denominator + 1)
        for p in range(len(count3)):
            key = str(p)
            if key not in count3:
                break
            bases = count3[key]
            max_cnt = 0
            max_base = ""
            for base, cnt in bases.items():
                if cnt > max_cnt:
                    max_cnt = cnt
                    max_base = base

            if max_cnt < threshold or max_base.upper() not in "ACTG":
                break
            trailer.append(max_base)
            threshold = max(min_pad_count, max_cnt // denominator + 1)
        trailer_str = "".join(trailer)

        sys.stdout.write(f"{header}{leader_str}{seq_str}{trailer_str}\n")
        if seq_alt_str and seq_alt_str != seq_str:
            sys.stdout.write(
                f"{header2}{leader_str}{seq_alt_str}{trailer_str}\n"
            )
    else:
        sys.stdout.write(f"{header}{seq_str}\n")
        if seq_alt_str and seq_alt_str != seq_str:
            sys.stdout.write(f"{header2}{seq_alt_str}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate alignment stats and build consensus FASTA"
    )
    parser.add_argument("stats", nargs="+", help="JSON stats files")
    parser.add_argument("-N", "--name", default="", help="Consensus name")
    parser.add_argument(
        "-M",
        "--min-pad-count",
        type=int,
        default=10,
        help="Min coverage for extension (default: 10)",
    )
    parser.add_argument(
        "-C",
        "--count-alt",
        type=int,
        default=LONG_MAX,
        help="Alt allele count threshold",
    )
    parser.add_argument(
        "-F",
        "--count-freq",
        type=float,
        default=2.0,
        help="Alt allele frequency threshold",
    )
    parser.add_argument(
        "-A",
        "--delete-by-ambiguity",
        action="store_true",
        help="Replace deleted positions with N",
    )
    parser.add_argument(
        "-S",
        "--skip-elongation",
        action="store_true",
        help="Skip aggregating leader/trailer overhang stats",
    )
    parser.add_argument(
        "--extend-terminus",
        action="store_true",
        help="Extend the consensus PAST the reference ends using leader/trailer soft-clip overhangs "
        "(the reach-gap feature for short/partial references). OFF by default: on a complete reference "
        "a terminal soft-clip is adapter / low-quality readthrough, not real genome — the aligner "
        "soft-clipped it precisely because it does not match — so appending it builds sequence from the "
        "very bases alignment rejected and runs off the end of the reference. Gap-repair WITHIN the "
        "reference is unaffected (it is folded into position_counts upstream in *_align_stats.py).",
    )
    parser.add_argument(
        "-K",
        "--keep-deleted",
        default=None,
        help="Fill deleted positions from reference FASTA",
    )
    parser.add_argument(
        "-O",
        "--denominator",
        type=int,
        default=2,
        help="Coverage denominator for extension threshold (default: 2)",
    )
    parser.add_argument(
        "-D",
        "--debug-mode",
        action="store_true",
        help="Print per-position counts to stderr",
    )
    parser.add_argument(
        "--position-bias-frac",
        type=float,
        default=0.0,
        help="Prefer the spanning consensus when the plurality base is >= this fraction "
        "terminal-supported and spanning reads agree on a different base. 0 = off (default).",
    )
    parser.add_argument(
        "--position-bias-min-span",
        type=int,
        default=5,
        help="Minimum spanning reads required to override a terminal-dominated plurality.",
    )
    parser.add_argument(
        "--min-ref-depth",
        type=int,
        default=0,
        help="min_tcc for the reference: don't set/flip a reference base on fewer reads than this "
        "(keep the prior ref base below the floor). 0 = off. Requires del_type=REF (ref staged).",
    )
    args = parser.parse_args()

    denominator = max(2, args.denominator)
    elongate = not args.skip_elongation

    count, count5, count3, ref_len, count_term = load_and_aggregate(args.stats, elongate)

    # Load reference for keep-deleted mode
    ref_sites = None
    keep_deleted = False
    if args.keep_deleted:
        ref_seq, ref_len_check = read_reference(args.keep_deleted)
        if ref_seq:
            ref_sites = list(ref_seq.lower())
            keep_deleted = True
        else:
            print(
                "WARNING (combine_align_stats): no reference found, "
                "turning off keep-deleted.",
                file=sys.stderr,
            )

    build_consensus(
        count=count,
        count5=count5,
        count3=count3,
        ref_len=ref_len,
        elongate=elongate,
        name=args.name,
        alt_count_thresh=args.count_alt,
        alt_freq_thresh=args.count_freq,
        delete_by_ambig=args.delete_by_ambiguity,
        keep_deleted=keep_deleted,
        ref_sites=ref_sites,
        min_pad_count=args.min_pad_count,
        denominator=denominator,
        debug=args.debug_mode,
        count_term=count_term,
        pb_bias_frac=args.position_bias_frac,
        pb_min_span=args.position_bias_min_span,
        min_ref_depth=args.min_ref_depth,
        extend_terminus=args.extend_terminus,
    )


if __name__ == "__main__":
    main()

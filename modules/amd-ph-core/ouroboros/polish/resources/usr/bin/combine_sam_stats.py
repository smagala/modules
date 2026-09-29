#!/usr/bin/env python3
"""combine_sam_stats.py — combine per-position SAM stats and emit the iterative consensus.

Combines per-position SAM stats (from sam_stats.py JSON) across one or more stat files and
emits the iterative consensus (and an alternative consensus when it differs). This is the
re-reference step of the final-assembly polish loop.

  python combine_sam_stats.py [options] <REF> <STAT1.json> <...>
    -N/--name, -I/--insertion-threshold, -D/--deletion-threshold,
    -i/--insertion-depth-threshold, -d/--deletion-depth-threshold,
    -A/--alternative-threshold, -C/--alternative-count, -M/--mark-deletions,
    -E/--min-dropout-edge-support

Tie-breaking: a deterministic (-count, allele) order. Exact ties are vanishingly rare at assembly
depth, so this choice is not consensus-affecting on real data.
"""
import sys, re, json

LONG_MAX = 2 ** 63 - 1

# (flag aliases) -> (dest, takes_value, caster)
_OPTS = {
    "-N": ("name", True, str), "--name": ("name", True, str),
    "-I": ("insT", True, float), "--insertion-threshold": ("insT", True, float),
    "-D": ("delT", True, float), "--deletion-threshold": ("delT", True, float),
    "-i": ("insD", True, int), "--insertion-depth-threshold": ("insD", True, int),
    "-d": ("delD", True, int), "--deletion-depth-threshold": ("delD", True, int),
    "-A": ("altT", True, float), "--alternative-threshold": ("altT", True, float),
    "-C": ("altC", True, int), "--alternative-count": ("altC", True, int),
    "-M": ("mark_del", False, None), "--mark-deletions": ("mark_del", False, None),
    "-E": ("minEdge", True, int), "--min-dropout-edge-support": ("minEdge", True, int),
    "-m": ("minDepth", True, int), "--min-depth": ("minDepth", True, int),
    # Mutation guard: hold a reference-CHANGING base to the variant-calling bar (mutating the
    # reference is tantamount to calling the variant). Reject a mutation whose allele has too little
    # data (< --mut-min-depth) or is essentially single-stranded (minor-strand fraction <
    # --mut-min-strand-frac, e.g. a one-strand amplicon primer artifact) -> keep the reference base.
    # 0 = off (no-op / parity). Both gate ONLY positions where plurality != reference.
    "--mut-min-depth": ("mutMinDepth", True, int),
    "--mut-min-strand-frac": ("mutMinStrandFrac", True, float),
    # Emit 'N' at columns with NO coverage instead of dropping them. Off by default for IRMA parity.
    # --min-depth already masks thin columns to 'N'; without this, a ZERO-depth column is the one
    # case that silently vanishes, so a coverage gap shortens the reference and shifts every
    # downstream coordinate. 'N' = unknown, '-' = known absent: only the latter should change length.
    "--mask-uncovered": ("maskUncovered", False, None),
    # Keep the REFERENCE base at columns with no coverage, instead of 'N'. Off by default.
    #
    # The reference this script emits is two things at once: the alignment target for the next
    # polish iteration, and (via the best-scoring copy) the reference CALL reports against. 'N'
    # is right for the second job and fatal for the first -- nothing aligns to 'N', so a window
    # that loses coverage once can never regain it, however much the rest of the reference
    # improves. Measured on a lyssavirus cohort 2026-08-31: an N run erodes by about one read
    # length from each flank and then stalls (716 -> 566 bp over twenty iterations), stranding
    # reads that align right up to the edge and are then soft-clipped.
    #
    # With this flag the column keeps real sequence so reads can still anchor there, and the
    # masking moves downstream where it belongs: CALL's --min-consensus-depth decides what is
    # reported. USE THE TWO TOGETHER. On its own this flag will report reference bases at
    # positions with no evidence, which is worse than the hole it fixes.
    "--carry-uncovered": ("carryUncovered", False, None),
    # Write "1" (an insertion or deletion was folded into the consensus this round) or "0" to PATH.
    # The authoritative signal for the stitch gate: a folded deletion is exactly what creates the
    # assembly-gap column stitch re-anchors / recovers; a folded insertion diverges length too. When
    # no indel is folded, the amended consensus embeds in the reference frame ungapped, so a BAM-free
    # stitch is byte-identical and the reads-on-reference remap (AMPLICON_CALL) can be skipped.
    "--indel-flag-file": ("indelFlagFile", True, str),
}

def parse_args(argv):
    o = {"name": None, "insT": 0.15, "delT": 0.75, "insD": 1, "delD": 1,
         "altT": 2.0, "altC": LONG_MAX, "mark_del": False, "minEdge": 0, "minDepth": 0,
         "maskUncovered": False, "carryUncovered": False,
         "mutMinDepth": 0, "mutMinStrandFrac": 0.0, "indelFlagFile": None}
    pos = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in _OPTS:
            dest, takes, cast = _OPTS[a]
            if takes:
                i += 1
                o[dest] = cast(argv[i])
            else:
                o[dest] = True
        elif a.startswith("--") and "=" in a:
            key, val = a.split("=", 1)
            dest, takes, cast = _OPTS[key]
            o[dest] = cast(val)
        else:
            pos.append(a)
        i += 1
    return o, pos

def first_ref_len(path):
    seq = first_ref_seq(path)
    return len(seq) if seq is not None else None

def first_ref_seq(path):
    with open(path) as fh:
        text = fh.read()
    for record in text.split(">"):
        if record == "":
            continue
        lines = re.split(r"\r\n|\n|\r", record)
        rseq = "".join(lines[1:])
        if len(rseq) < 1:
            continue
        return rseq.upper()
    return None

def sorted_alleles(counts):
    # descending by count, deterministic tie-break by allele
    return sorted(counts.keys(), key=lambda a: (-counts[a], a))

def main():
    o, pos = parse_args(sys.argv[1:])
    if len(pos) < 2:
        sys.exit("Usage:\t%s [options] <REF> <STAT1> <...>\n" % sys.argv[0])

    ref_path = pos[0]
    stat_files = pos[1:]
    N = first_ref_len(ref_path)
    if N is None:
        sys.exit("%s ERROR: no reference found in %s.\n" % (sys.argv[0], ref_path))

    insT = o["insT"]
    insD = o["insD"]
    delT = o["delT"]
    delD = o["delD"]
    altT = o["altT"]
    altC = o["altC"]
    mark_del = o["mark_del"]
    minEdge = max(0, int(o["minEdge"]))
    name = o["name"]

    mut_min_depth = max(0, int(o.get("mutMinDepth", 0)))
    mut_min_sf = max(0.0, float(o.get("mutMinStrandFrac", 0.0)))
    guard_on = mut_min_depth > 0 or mut_min_sf > 0
    # --carry-uncovered needs the reference bases too, not just the mutation guard.
    ref_seq = first_ref_seq(ref_path) if (guard_on or o.get("carryUncovered")) else None

    # aggregate
    big = [dict() for _ in range(N)]   # per pos: {base: count}
    strand = [dict() for _ in range(N)]  # per pos: {base: [fwd, rev]} (for the mutation guard)
    ins = {}                            # pos -> {insert: count}
    for path in stat_files:
        with open(path) as fh:
            d = json.load(fh)
        for ps, bc in d.get("base_counts", {}).items():
            p = int(ps)
            if p < 0 or p >= N:
                continue
            tgt = big[p]
            for base, cnt in bc.items():
                tgt[base] = tgt.get(base, 0) + cnt
        if guard_on:
            for ps, sc in d.get("strand_counts", {}).items():
                p = int(ps)
                if p < 0 or p >= N:
                    continue
                tgt = strand[p]
                for base, fr in sc.items():
                    cur = tgt.setdefault(base, [0, 0])
                    cur[0] += fr[0]; cur[1] += fr[1]
        for ps, ic in d.get("ins_counts", {}).items():
            p = int(ps)
            tgt = ins.setdefault(p, {})
            for insert, cnt in ic.items():
                tgt[insert] = tgt.get(insert, 0) + cnt

    # plurality + totals
    cons = [""] * N
    totals = [0] * N
    for p in range(N):
        total = 0
        con = ""
        mx = None
        for allele in sorted(big[p].keys()):
            cnt = big[p][allele]
            total += cnt
            if mx is None or cnt > mx:
                mx = cnt
                con = allele
        # Mutation guard: a reference-CHANGING plurality base must clear the variant-calling bar,
        # since mutating the reference is tantamount to calling the variant. Reject (keep the ref
        # base) when the mutating allele has too little total data or is essentially single-stranded
        # (a one-strand amplicon/primer artifact). No-op where plurality == reference or gate off.
        if guard_on and ref_seq is not None and con and con not in ("-",) and p < len(ref_seq):
            rb = ref_seq[p]
            if con != rb and rb in "ACGTacgt":
                fwd, rev = strand[p].get(con, [0, 0])
                sc_tot = fwd + rev
                minor_frac = (min(fwd, rev) / sc_tot) if sc_tot else 0.0
                too_thin = mut_min_depth > 0 and total < mut_min_depth
                one_strand = mut_min_sf > 0 and minor_frac < mut_min_sf
                if too_thin or one_strand:
                    con = rb.upper()   # decline the mutation; keep the reference base

        cons[p] = con
        totals[p] = total

    # dropout edge masking
    if minEdge > 0:
        plurality = "".join("." if c == "" else c for c in cons)
        for m in re.finditer(r"[ATGCNatcgn]{6}[.]{91,}[ATCGNatcgn]{6}", plurality):
            for p in range(m.start(), m.end()):
                if totals[p] < minEdge:
                    cons[p] = "N"

    # global depth mask: any called position below min_depth -> N (0 = no-op). Post-alignment
    # depth honesty: don't call a base where coverage is below the floor (amplicon/thin regions).
    minDepth = max(0, int(o.get("minDepth", 0)))
    if minDepth > 0:
        for p in range(N):
            if cons[p] not in ("", "-") and totals[p] < minDepth:
                cons[p] = "N"

    # Uncovered columns: a zero-depth column has cons[p] == "" and would append nothing at output,
    # deleting the position and shifting everything 3' of it. The depth mask above cannot reach it
    # (it skips ""), so zero coverage is the one case that loses its marker. Mask it like any other
    # unknown base, so the reference stays full-length and coordinate-anchored.
    if o.get("maskUncovered"):
        for p in range(N):
            if cons[p] == "":
                cons[p] = "N"

    # Uncovered columns, kept ALIGNABLE. Runs after the mask above so --carry-uncovered wins when
    # both are set: the two disagree only at zero-depth columns, and that disagreement is the
    # whole point. Anything already masked to 'N' by --min-depth is left alone -- a thin column
    # has evidence and was judged, an uncovered one was never seen.
    if o.get("carryUncovered") and ref_seq is not None:
        for p in range(N):
            if (cons[p] in ("", "N") and p < len(ref_seq)
                    and ref_seq[p] in "ACGT" and totals[p] == 0):
                cons[p] = ref_seq[p]

    if name:
        header = ">" + name + "\n"
        header2 = header
    else:
        header = ">consensus\n"
        header2 = ">alternative\n"

    consensus = ""
    alternative = ""
    # Track whether any insertion or deletion is folded into the consensus (for the stitch gate).
    indel_folded = False
    for p in range(N):
        if cons[p] != "-":
            consensus += cons[p]
            alleles = list(big[p].keys())
            if len(alleles) > 1:
                sa = sorted_alleles(big[p])
                if sa[1] == "-":
                    alternative += cons[p]
                else:
                    altCount = big[p][sa[1]]
                    altFreq = altCount / totals[p]
                    if altFreq < altT or altCount < altC:
                        alternative += cons[p]
                    else:
                        alternative += sa[1]
            else:
                alternative += cons[p]
        else:
            freq = big[p][cons[p]] / totals[p]
            alleles = list(big[p].keys())
            if (big[p][cons[p]] < delD or freq < delT) and len(alleles) > 1:
                sa = sorted_alleles(big[p])
                consensus += sa[1]
                if len(alleles) > 2:
                    altCount = big[p][sa[2]]
                    altFreq = altCount / totals[p]
                    if altFreq < altT or altCount < altC:
                        alternative += sa[1]
                    else:
                        alternative += sa[2]
                else:
                    alternative += sa[1]
            else:
                # The deletion clears the threshold -> a real folded deletion: the reference column is
                # dropped (shortens) or kept as '-' under mark_del. Either way the consensus diverges
                # from the reference by a deletion -> the assembly-gap column stitch adjudicates.
                indel_folded = True
                if mark_del:
                    consensus += cons[p]
                    alternative += cons[p]

        if p in ins:
            si = sorted_alleles(ins[p])
            if p < (N - 1):
                avgTotal = int((totals[p] + totals[p + 1]) / 2)
            else:
                avgTotal = totals[p]
            freq = ins[p][si[0]] / avgTotal
            if freq >= insT and ins[p][si[0]] >= insD:
                consensus += si[0].lower()
                alternative += si[0].lower()
                indel_folded = True

    if o.get("indelFlagFile"):
        with open(o["indelFlagFile"], "w") as fh:
            fh.write("1\n" if indel_folded else "0\n")

    sys.stdout.write(header + consensus + "\n")
    if alternative != "" and alternative != consensus:
        sys.stdout.write(header2 + alternative + "\n")

if __name__ == "__main__":
    main()

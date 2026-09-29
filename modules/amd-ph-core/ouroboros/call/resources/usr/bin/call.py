#!/usr/bin/env python3
"""call.py — call single-nucleotide / indel variants, write tables, make plurality consensus.

Reads one or more per-partition msgpack pileups produced by var_call_stats.py, combines them, then for
every reference position chooses the plurality consensus and evaluates each minority allele against a
frequency / count / quality / confidence / significance-bound gate. Writes six text outputs
(-coverage.txt, .fasta, -allAlleles.txt, -variants.txt, -insertions.txt, -deletions.txt) plus the
phasing inputs -vars.msgpack / -pats.msgpack that phase.py reads.

Doubles are formatted to 15 significant digits and integers bare. Per-position allele rows in
-allAlleles.txt are emitted in sorted order.

Minority-variant floor:
  with -A (AUTO_F): the reported floor is the highest minority frequency at or below the per-base
                    error estimate (1/10^(Q/10)); variants at or below it are dropped. This floor is
                    derived from read quality, so it shifts with the input quality distribution.
  without -A:       the only floor is min_freq (-F), enforced by the valid-variant gate, so a fixed
                    -F is an explicit floor independent of read quality. -allAlleles.txt lists every
                    observed allele with full stats regardless, so the sub-threshold tail stays
                    visible — calling is a flag, not a filter.
"""
import argparse
import math
import re
import sys

import msgpack
from scipy import stats

ap = argparse.ArgumentParser(add_help=True)
ap.add_argument("-G", "--no-gap-allele", action="store_true", dest="no_gap")
ap.add_argument("-F", "--min-freq", type=float, dest="min_freq")
ap.add_argument("-I", "--min-insertion-freq", type=float, dest="min_freq_ins")
ap.add_argument("-D", "--min-deletion-freq", type=float, dest="min_freq_del")
ap.add_argument("-C", "--min-count", type=int, dest="min_count")
ap.add_argument("-Q", "--min-quality", type=int, dest="min_quality")
ap.add_argument("-T", "--min-total-col-coverage", type=int, dest="min_total")
ap.add_argument("--min-consensus-depth", type=int, dest="min_cons_depth", default=None,
                help="Emit 'N' rather than a plurality base at columns below this depth. 0/unset "
                     "keeps the legacy behaviour of calling from any depth >= 1.")
ap.add_argument("-P", "--print-all-sites", action="store_true", dest="print_all")
ap.add_argument("-M", "--conf-not-mac-err", type=float, dest="min_conf")
ap.add_argument("-S", "--sig-level", type=float, dest="sig_level")
ap.add_argument("-E", "--paired-error", dest="paired_stats")
ap.add_argument("-B", "--call-table", dest="call_table")
ap.add_argument("-A", "--auto-min-freq", action="store_true", dest="auto_freq")
ap.add_argument("--caller", choices=("fixed", "binomial", "eb"), default="fixed",
                help="minority-variant significance rule: fixed = legacy Wilson UB + min-freq floor; "
                     "binomial = exact binomial vs an error-floor null (depth-adaptive); "
                     "eb = capped empirical-Bayes posterior (depth-adaptive, self-calibrating)")
ap.add_argument("--systematic-error", type=float, default=0.005, dest="systematic",
                help="per-position systematic error floor for the binomial/eb caller (default 0.005)")
ap.add_argument("--cap-prior-mean", type=float, default=0.02, dest="cap_prior_mean",
                help="ceiling on the empirical-Bayes prior mean (eb caller only)")
ap.add_argument("--sor", type=float, default=0.0,
                help="SOR strand-bias gate: drop a called variant if GATK StrandOddsRatio > this "
                     "(0 = off; GATK SNP hard-filter = 3.0)")
ap.add_argument("--refchg-min-strand-frac", type=float, default=0.0, dest="refchg_min_sf",
                help="consensus reference-change guard: when the plurality consensus base DIFFERS "
                     "from the reference base, keep the reference base if the consensus allele is "
                     "essentially single-stranded (minor-strand fraction < this) — a one-strand "
                     "amplicon/primer artifact you wouldn't CALL, so don't consensus-mutate to it. "
                     "The artifact then reports as a minority variant instead. (0 = off/parity)")
ap.add_argument("ref")
ap.add_argument("prefix")
ap.add_argument("pileups", nargs="+")
a = ap.parse_args()

# --- defaults ---
no_gap = a.no_gap
min_count = 2 if a.min_count is None else max(a.min_count, 0)
min_freq = 0.005 if a.min_freq is None else max(a.min_freq, 0)
min_freq_ins = min_freq if a.min_freq_ins is None else max(a.min_freq_ins, 0)
min_freq_del = min_freq if a.min_freq_del is None else max(a.min_freq_del, 0)
min_conf = 0.5 if a.min_conf is None else max(a.min_conf, 0)
min_quality = 20 if a.min_quality is None else max(a.min_quality, 0)
min_total = 2 if (a.min_total is None or a.min_total < 0) else a.min_total
min_cons_depth = 0 if a.min_cons_depth is None else max(0, a.min_cons_depth)
auto_freq = a.auto_freq
print_all = a.print_all
caller = a.caller
systematic = max(a.systematic, 0.0)
cap_prior_mean = a.cap_prior_mean
sor_max = max(a.sor, 0.0)
refchg_min_sf = max(a.refchg_min_sf, 0.0)


# FIXME legacy perl logic: the whole numeric format (15 sig-figs, integers bare, 'NA' sentinel) only
# exists to byte-match perl's default stringification. Drop for normal rounded output once parity ends.
# --- number formatting: integers bare, doubles to 15 significant digits ---
def g(x):
    """Integers bare, doubles via %.15g; pass 'NA' through."""
    if isinstance(x, str):
        return x
    if isinstance(x, int):
        return str(x)
    return "%.15g" % x


def calc_prob(w, e):
    if e > w:
        return 0.0
    return (w - e) / w


# --- depth-adaptive caller (binomial / empirical-Bayes); replaces the fixed VAF floor + Wilson UB ---
def error_floor(aq, systematic):
    """Per-position error rate: the worse (higher) of the quality-derived rate and a systematic floor.
    `aq` is the allele's mean Phred; 'NA' (deletion rows) falls back to a 0.25 ceiling."""
    q = 10.0 ** (-aq / 10.0) if isinstance(aq, (int, float)) and aq > 0 else 0.25
    return max(q, systematic)


def binom_sf(k, D, eps):
    """Exact upper-tail P(X >= k | D, eps) for a Binomial(D, eps) error null."""
    if k <= 0 or D <= 0:
        return 1.0
    return float(stats.binom.sf(k - 1, D, eps))


def fit_beta_prior(freqs):
    """Empirical-Bayes prior: method-of-moments Beta fit to the observed minority-allele frequency
    spectrum (mostly error, so the fit lands near the error rate). Jeffreys Beta(0.5,0.5) on degeneracy."""
    xs = [x for x in freqs if 0.0 < x < 1.0]
    if len(xs) < 5:
        return 0.5, 0.5
    m = sum(xs) / len(xs)
    v = sum((x - m) ** 2 for x in xs) / len(xs)
    if v <= 0 or v >= m * (1 - m):
        return 0.5, 0.5
    common = m * (1 - m) / v - 1
    a0, b0 = m * common, (1 - m) * common
    return (a0, b0) if (a0 > 0 and b0 > 0) else (0.5, 0.5)


def cap_prior(a0, b0, cap_mean):
    """Cap the EB prior mean at a sane error ceiling so a junk-dominated spectrum can't teach a permissive
    prior; keep the learned concentration."""
    s = a0 + b0
    m = a0 / s if s > 0 else 0.5
    if m <= cap_mean:
        return a0, b0
    return cap_mean * s, (1.0 - cap_mean) * s


def eb_prob_real(eps, a0, b0, k, D):
    """Posterior P(p > eps | data) under Beta(a0,b0) -> Beta(a0+k, b0+D-k)."""
    return float(1.0 - stats.beta.cdf(eps, a0 + k, b0 + D - k))


def sor_score(rf, rr, af, ar):
    """GATK StrandOddsRatio on the consensus-vs-minority 2x2 [[refFwd,refRev],[altFwd,altRev]],
    pseudocount 1. Depth-robust (unlike FisherStrand) and does not penalize legitimate one-strandedness
    at amplicon ends. Higher = more strand-biased; GATK's SNP hard-filter threshold is 3.0."""
    rf, rr, af, ar = rf + 1, rr + 1, af + 1, ar + 1
    ratio = (rf / rr) * (ar / af) + (rr / rf) * (af / ar)
    return math.log(ratio) + math.log(min(rf, rr) / max(rf, rr)) - math.log(min(af, ar) / max(af, ar))


# --- significance upper bound (UB) ---
take_sig = a.sig_level is not None
if take_sig:
    sig = a.sig_level
    if sig >= 1:
        sig /= 100
    if sig >= .999:
        kappa = 3.090232
    elif sig >= .99:
        kappa = 2.326348
    elif sig >= .95:
        kappa = 1.644854
    elif sig >= .90:
        kappa = 1.281552
    else:
        kappa = 3.090232
    kappa2 = kappa ** 2
    eta = kappa2 / 3 + 1 / 6
    gamma1 = kappa2 * (13 / 18) + 17 / 18
    gamma2 = kappa2 * (1 / 18) + 7 / 36

# significance thresholds for the depth-adaptive callers: alpha for the binomial tail, sig_thresh for
# the EB posterior. Both track -S/sig_level (default 0.999 -> alpha 0.001).
alpha = (1.0 - sig) if take_sig else 0.001
sig_thresh = sig if take_sig else 0.999


def UB(p, N):
    if not take_sig:
        return 0.0
    if N <= 0:
        sys.stderr.write("Unexpected error: %s coverage depth.\n" % N)
        return 0
    if p == 1:
        return 1
    V = p - p ** 2
    u2 = (p * N + eta) / (N + 2 * eta)
    in_root = V + (gamma2 - gamma1 * V) / N
    if in_root < 0:
        return 1
    ub = u2 + kappa * math.sqrt(in_root) / math.sqrt(N)
    return max(min(ub, 1), 0)


# --- read reference (first record only) ---
def read_first_fasta(path):
    name, parts, seen = None, [], False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if seen:
                    break
                name = line[1:]
                seen = True
            elif seen:
                parts.append(line)
    return name, "".join(parts)


ref_name, ref_seq = read_first_fasta(a.ref)
if ref_seq is None or len(ref_seq) < 1:
    raise SystemExit("No reference found.")
ref_len = len(ref_seq)
ref_seq_u = ref_seq.upper()   # for the consensus reference-change guard (case-insensitive compare)

# --- call table (-B), seeds variants to force-report ---
variants = {}        # {pos: {base: freq}}
do_call_table = False
if a.call_table:
    with open(a.call_table) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            p, base = line.split("\t")[:2]
            variants.setdefault(int(p) - 1, {})[base.upper()] = 0
    do_call_table = len(variants) > 0

# --- paired error rates (-E) ---
DE = PE = IE = 0
if a.paired_stats:
    pstats = {}
    with open(a.paired_stats) as fh:
        for line in fh:
            line = line.rstrip("\n")
            rn, typ, value = line.split("\t")
            pstats.setdefault(rn, {})[typ] = value
    DE = float(pstats[ref_name]["MinimumDeletionErrorRate"])
    PE = float(pstats[ref_name]["ExpectedErrorRate"])
    IE = float(pstats[ref_name]["MinimumInsertionErrorRate"])

# --- combine pileups ---
cTable = [dict() for _ in range(ref_len)]
qTable = [dict() for _ in range(ref_len)]
sTable = [dict() for _ in range(ref_len)]   # {allele: [fwd, rev]} for the SOR strand-bias gate
icTable = {}
iqTable = {}
dcTable = {}
alignments = {}


def add(d, k, v):
    d[k] = d.get(k, 0) + v


for path in a.pileups:
    with open(path, "rb") as fh:
        data = msgpack.unpackb(fh.read(), raw=False, strict_map_key=False)
    for aln, c in data["aln"].items():
        add(alignments, aln, c)
    counts, quals = data["counts"], data["quals"]
    strand_in = data.get("strand")
    for p in range(ref_len):
        for allele, c in counts[p].items():
            add(cTable[p], allele, c)
            qTable[p][allele] = qTable[p].get(allele, 0) + quals[p].get(allele, 0)
            if strand_in is not None:
                fr = strand_in[p].get(allele)
                if fr:
                    st = sTable[p].setdefault(allele, [0, 0])
                    st[0] += fr[0]
                    st[1] += fr[1]
    for p, ins in data["ins_c"].items():
        for insert, c in ins.items():
            icTable.setdefault(p, {})
            add(icTable[p], insert, c)
    for p, insq in data["ins_q"].items():
        for insert, q in insq.items():
            iqTable.setdefault(p, {})
            iqTable[p][insert] = iqTable[p].get(insert, 0) + q
    for p, dels in data["del_c"].items():
        for inc, c in dels.items():
            dcTable.setdefault(p, {})
            add(dcTable[p], int(inc), c)

prefix = a.prefix
covg = open(prefix + "-coverage.txt", "w")
cons = open(prefix + ".fasta", "w")
alla = None
if print_all or do_call_table:
    alla = open(prefix + "-allAlleles.txt", "w")
    alla.write("Reference_Name\tPosition\tAllele\tCount\tTotal\tFrequency\t"
               "Average_Quality\tConfidenceNotMacErr\tPairedUB\tQualityUB\tAllele_Type\n")
vars_fh = open(prefix + "-variants.txt", "w")
vars_fh.write("Reference_Name\tPosition\tTotal\tConsensus_Allele\tMinority_Allele\t"
              "Consensus_Count\tMinority_Count\tConsensus_Frequency\tMinority_Frequency\t"
              "Consensus_Average_Quality\tMinority_Average_Quality\t"
              "ConfidenceNotMacErr\tPairedUB\tQualityUB\n")
covg.write("Reference_Name\tPosition\tCoverage Depth\tConsensus\tDeletions\tAmbiguous\t"
           "Consensus_Count\tConsensus_Average_Quality\n")
cons.write(">" + ref_name + "\n")

hFreq = 0.0
alphabet = list("ACGT-")
totals = [0] * ref_len
var_line = {}
consensus_seq_parts = []


def alla_row(p, base, count, total, freq, quality, confidence, paired_ub, quality_ub, btype):
    alla.write("\t".join([ref_name, str(p + 1), base, g(count), g(total), g(freq), g(quality),
                          g(confidence), g(paired_ub), g(quality_ub), btype]) + "\n")


# empirical-Bayes prior, fit once from this sample's whole minority-allele frequency spectrum (the
# spectrum is dominated by error, so the prior self-calibrates near the error rate). Capped so a
# junk-dominated spectrum can't teach a permissive prior. Only needed for --caller eb.
eb_a0 = eb_b0 = None
if caller == "eb":
    minor_freqs = []
    for cell in cTable:
        tot = sum(c for b, c in cell.items() if b != "N")
        if tot <= 0:
            continue
        con_c = max((c for b, c in cell.items() if b != "-"), default=0)
        for b, c in cell.items():
            if b == "N" or c == 0 or c == con_c:
                continue
            f = c / tot
            if 0.0 < f < 1.0:
                minor_freqs.append(f)
    a0, b0 = fit_beta_prior(minor_freqs)
    eb_a0, eb_b0 = cap_prior(a0, b0, cap_prior_mean)


for p in range(ref_len):
    # FIXME legacy perl logic: plurality pick uses strict > so the first max-count base in iteration
    # order wins ties (nondeterministic in perl), and the con_count=-1 sentinel yields degenerate
    # arithmetic (e.g. con_quality=-33) at empty columns. Replace with an explicit deterministic tie
    # rule and real "no coverage" handling once parity ends.
    #
    # No-coverage masking: a column with no reads emits 'N' (unknown base, full length preserved) rather
    # than the legacy '.' sentinel. '-' is reserved for a genuine all-gap (deletion) column. amend_consensus
    # strips '-' (a deletion shortens the consensus) but keeps 'N', so the amended consensus stays
    # full-length and coordinate-anchored, with uncovered positions honestly marked unknown.
    consensus = "N"
    con_count = -1
    total = 0
    bases = list(cTable[p].keys())
    n_alleles = len(bases)

    if n_alleles == 1:
        consensus = bases[0]
        total = con_count = cTable[p][consensus]
    else:
        for base in bases:
            if cTable[p][base] > con_count and base != "-":
                con_count = cTable[p][base]
                consensus = base
            total += cTable[p][base]

    total -= cTable[p].get("N", 0)
    totals[p] = total

    # Consensus reference-change guard: mutating the (consensus) reference base is tantamount to calling
    # the variant, so a reference-CHANGING plurality base must clear the variant-calling strand bar.
    # If the plurality differs from the reference and its allele is essentially single-stranded (a
    # one-strand amplicon/primer artifact you wouldn't CALL), keep the reference base — the artifact
    # then falls through to the minority-variant table below. No-op where plurality == ref, no strand
    # data, or the gate is off. Real both-strand WGS variants pass (minor-strand frac ~0.5) -> parity.
    if (refchg_min_sf > 0 and con_count > 0 and p < len(ref_seq_u)):
        rb = ref_seq_u[p]
        if consensus != rb and rb in "ACGT":
            fwd, rev = sTable[p].get(consensus, (0, 0))
            sc_tot = fwd + rev
            if sc_tot > 0 and (min(fwd, rev) / sc_tot) < refchg_min_sf:
                consensus = rb   # decline the consensus mutation; keep the reference base
                con_count = cTable[p].get(rb, 0)   # re-anchor freq/quality/classification on the ref base

    # Depth floor on the CONSENSUS. min_total has always gated variant and indel calls, but never
    # the plurality base, so a column with a single aligned read produced a called base. Exhaustive
    # aligners place spurious reads in uncovered regions, and at depth 1-3 those became consensus:
    # measured as 12 wrong bases in a 60 bp window that Sanger, viralrecon and three independent
    # de novo contigs all agree on. 'N' here is the documented "unknown" marker and preserves length.
    if min_cons_depth > 0 and total < min_cons_depth and consensus != "-":
        consensus = "N"

    consensus_seq_parts.append(consensus)
    cons.write(consensus)

    con_freq = con_count / total if total != 0 else 0
    # FIXME legacy perl logic: qTable holds summed Phred+33 ASCII, so mean quality subtracts count*33
    # everywhere. Store summed Phred in var_call_stats and drop the *33 offset once parity ends.
    con_quality = ((qTable[p].get(consensus, 0) - con_count * 33) / con_count) if con_count != 0 else min_quality

    gaps = cTable[p].get("-")
    depth = (total - gaps) if gaps is not None else total
    ambig = cTable[p].get("N", 0)
    covg.write("\t".join([ref_name, str(p + 1), g(depth), consensus, g(gaps if gaps is not None else 0),
                          g(ambig), g(con_count), g(con_quality)]) + "\n")

    if do_call_table:
        # FIXME legacy perl logic: the -B call-table branch is a rarely-used IRMA mode kept only for
        # parity; it duplicates much of the minority-allele logic below. Consolidate or drop later.
        if p not in variants:
            continue
        for base in alphabet:
            if base in cTable[p]:
                count = cTable[p][base]
                freq = count / total
                if base == "-":
                    quality_ub = quality = confidence = "NA"
                    paired_ub = UB(DE, total)
                else:
                    quality = (qTable[p].get(base, 0) - count * 33) / count
                    ee = 1 / (10 ** (quality / 10))
                    confidence = calc_prob(freq, ee)
                    paired_ub = UB(PE, total)
                    quality_ub = UB(ee, total)
            else:
                freq = count = 0
                quality_ub = quality = confidence = "NA"
                if total == 0:
                    paired_ub = "NA"
                elif base == "-":
                    paired_ub = UB(DE, total)
                else:
                    paired_ub = UB(PE, total)
            btype = "Consensus" if base == consensus else "Minority"
            alla_row(p, base, count, total, freq, quality, confidence, paired_ub, quality_ub, btype)
            if base in variants.get(p, {}):
                variants[p][base] = freq
                var_line.setdefault(p, {})[base] = "\t".join(
                    [ref_name, str(p + 1), g(total), consensus, base, g(con_count), g(count),
                     g(con_freq), g(freq), g(con_quality), g(quality),
                     g(confidence), g(paired_ub), g(quality_ub)]) + "\n"
        continue

    for base in sorted(bases):
        if base == consensus:
            if print_all:
                if base == "N":
                    ee = 1 / (10 ** (con_quality / 10))
                    confidence = calc_prob(con_freq, ee)
                    quality = con_quality
                    paired_ub = UB(PE, con_count)
                    quality_ub = UB(ee, con_count)
                    t = con_count
                elif base == "-":
                    quality = confidence = "NA"
                    paired_ub = UB(DE, total)
                    quality_ub = 0
                    t = total
                else:
                    ee = 1 / (10 ** (con_quality / 10))
                    confidence = calc_prob(con_freq, ee)
                    quality = con_quality
                    paired_ub = UB(PE, total)
                    quality_ub = UB(ee, total)
                    t = total
                alla_row(p, base, con_count, t, con_freq, quality, confidence, paired_ub, quality_ub, "Consensus")
        else:
            count = cTable[p][base]
            if count == 0 or base == "N":
                continue
            freq = count / total
            if base != "-" and base != "N":
                quality = (qTable[p].get(base, 0) - count * 33) / count
            else:
                quality = min_quality

            # valid variant. The fixed VAF floor (min_freq) applies only to --caller fixed; the
            # depth-adaptive callers replace it with the binomial/EB significance test below.
            floor_ok = (freq >= min_freq) if caller == "fixed" else True
            if (not (no_gap and base == "-")
                    and floor_ok and count >= min_count
                    and quality >= min_quality and total >= min_total):
                if base == "-":
                    confidence = quality = "NA"
                    paired_ub = UB(DE, total)
                    quality_ub = 0
                    ee = 0
                else:
                    ee = 1 / (10 ** (quality / 10))
                    confidence = calc_prob(freq, ee)
                    paired_ub = UB(PE, total)
                    quality_ub = UB(ee, total)

                # FIXME legacy perl logic / open design concern: hFreq is AUTO_F's reporting floor —
                # the highest minority freq within its per-base error estimate. It is a single GLOBAL,
                # depth-blind scalar, so one shallow/noisy position can raise it for the whole segment.
                # Flagged for revisit; resolution undecided (fixed -F, or a depth-aware/per-position
                # floor). See lab-notebook 2026-06-16-auto-f-hfreq-global-floor-concern.
                if base != "-" and freq <= ee and freq > hFreq:
                    hFreq = freq
                if print_all:
                    alla_row(p, base, count, total, freq, quality, confidence, paired_ub, quality_ub, "Minority")

                # a '-' minority has confidence 'NA', which counts as 0 here, so it always fails
                # confidence < min_conf and is skipped (gaps are reported in -deletions.txt instead).
                conf_num = 0 if confidence == "NA" else confidence
                # significance: fixed = legacy Wilson upper bounds; binomial/eb = depth-adaptive test
                # vs an error-floor null. The paired-error UB guard is kept in all modes.
                if caller == "fixed":
                    sig_fail = (freq <= quality_ub)
                elif caller == "binomial":
                    eps = error_floor(quality, systematic) if base != "-" else 1.0
                    sig_fail = (binom_sf(count, total, eps) >= alpha)
                else:  # eb
                    eps = error_floor(quality, systematic) if base != "-" else 1.0
                    sig_fail = (eb_prob_real(eps, eb_a0, eb_b0, count, total) <= sig_thresh)
                if conf_num < min_conf or freq <= paired_ub or sig_fail:
                    continue

                # SOR strand-bias gate: drop a call whose minority allele is strand-skewed relative to
                # the consensus allele (a classic sequencing/PCR/primer-strand artifact). Off when sor_max=0.
                if sor_max > 0:
                    rf, rr = sTable[p].get(consensus, (0, 0))
                    af, ar = sTable[p].get(base, (0, 0))
                    if sor_score(rf, rr, af, ar) > sor_max:
                        continue

                variants.setdefault(p, {})[base] = freq
                var_line.setdefault(p, {})[base] = "\t".join(
                    [ref_name, str(p + 1), g(total), consensus, base, g(con_count), g(count),
                     g(con_freq), g(freq), g(con_quality), g(quality),
                     g(confidence), g(paired_ub), g(quality_ub)]) + "\n"

            # FIXME legacy perl logic: this "any variant" branch re-derives ee/confidence/UBs already
            # computed above purely to print sub-threshold alleles to allAlleles. Once allAlleles is
            # produced by a single uniform per-allele pass, this duplicate path collapses.
            # any variant (only emitted to allAlleles when -P)
            elif print_all:
                if base == "-":
                    quality = confidence = "NA"
                    paired_ub = UB(DE, total)
                    quality_ub = 0
                    ee = 0
                else:
                    ee = 1 / (10 ** (quality / 10))
                    confidence = calc_prob(freq, ee)
                    paired_ub = UB(PE, total)
                    quality_ub = UB(ee, total)
                if base != "-" and freq <= ee and freq > hFreq:
                    hFreq = freq
                alla_row(p, base, count, total, freq, quality, confidence, paired_ub, quality_ub, "Minority")

cons.write("\n")
cons.close()
covg.close()

consensus_seq = "".join(consensus_seq_parts)

# FIXME legacy perl logic: with -A this drops every called variant at or below the data-derived hFreq
# floor (see the hFreq design-concern note above). Without -A the floor is just min_freq (-F). If
# AUTO_F is dropped this loop reduces to "write every valid variant in sorted order". The sort key
# (the rendered line text) is also a perl artifact — a real key would be (position, allele).
for p in sorted(var_line):
    for base in sorted(var_line[p], key=lambda b: var_line[p][b]):
        # AUTO_F's depth-blind hFreq floor applies only to --caller fixed; the depth-adaptive
        # callers already floor per-position via the binomial/EB test, so they write every call.
        if not do_call_table and auto_freq and caller == "fixed":
            if variants[p][base] > hFreq:
                vars_fh.write(var_line[p][base])
            else:
                del variants[p][base]
        else:
            vars_fh.write(var_line[p][base])
vars_fh.close()
if alla is not None:
    alla.close()


# FIXME legacy perl logic: insertion/deletion depth is computed by turning each read's full-length
# alignment string into "start,stop" covered-range strings, bucketing them, then summing ranges that
# span the indel. A plain interval/coverage structure (or a CIGAR-derived depth array) would replace
# this string-encoding round-trip entirely once parity ends.
# --- alignment coordinate support, for insertion/deletion depth ---
def decode_aln(key):
    # Inflate a compact "<start>:<body>" coordinate key (from var_call_stats) back to the full
    # ref_len-padded string the legacy logic expects. Done one pattern at a time, so the persistent
    # {pattern: count} table stays compact while this transient string is O(ref_len) and freed.
    s, _, body = key.partition(":")
    start = int(s)
    return "." * start + body + "." * (ref_len - start - len(body))


def to_indices_zero(aln):
    coords = []
    index = 0
    for m in re.finditer(r"(\.+|[^.]+)", aln):
        run = m.group(1)
        length = len(run)
        if run[0] != ".":
            start = index
            index += length
            coords.append("%d,%d" % (start, index - 1))
        else:
            index += length
    return ";".join(coords)


coord_list = {}
for aln, c in alignments.items():
    add(coord_list, to_indices_zero(decode_aln(aln)), c)

coord_support = {}
for los, c in coord_list.items():
    if not los:
        continue
    for coord in los.split(";"):
        start, stop = (int(x) for x in coord.split(","))
        coord_support.setdefault(start, {})
        coord_support[start][stop] = coord_support[start].get(stop, 0) + c

coord_starts = sorted(coord_support)
coord_stops = {start: sorted(coord_support[start], reverse=True) for start in coord_starts}


def span_depth(p, pp):
    total = 0
    for start in coord_starts:
        if start <= p:
            for stop in coord_stops[start]:
                if pp <= stop:
                    total += coord_support[start][stop]
                else:
                    break
        else:
            break
    return total


# --- insertions ---
with open(prefix + "-insertions.txt", "w") as insv:
    insv.write("Reference_Name\tUpstream_Position\tInsert\tContext\tCalled\tCount\tTotal\t"
               "Frequency\tAverage_Quality\tConfidenceNotMacErr\tPairedUB\tQualityUB\n")
    for p in sorted(icTable):
        pp = p + 1
        total = span_depth(p, pp)
        for insert in sorted(icTable[p]):
            count = icTable[p][insert]
            if count < min_count:
                continue
            called = "TRUE"
            quality = iqTable[p][insert] / count if count > 0 else 0
            if quality < min_quality:
                called = "FALSE"
            freq = count / total if total > 0 else 0
            if freq < min_freq_ins or total < min_total:
                called = "FALSE"
            EE = 1 / (10 ** (quality / 10))
            confidence = calc_prob(freq, EE)
            paired_ub = UB(IE, total)
            quality_ub = UB(EE, total)
            if confidence < min_conf or freq <= paired_ub or freq <= quality_ub:
                called = "FALSE"
            if p < 5:
                left = consensus_seq[0:p + 1]
            else:
                left = consensus_seq[p - 4:p + 1]
            if p > (ref_len - 6):
                right = consensus_seq[pp:ref_len]
            else:
                right = consensus_seq[pp:pp + 5]
            context = left.lower() + insert.upper() + right.lower()
            insv.write("\t".join([ref_name, str(p + 1), insert.upper(), context, called,
                                  g(count), g(total), g(freq), g(quality), g(confidence),
                                  g(paired_ub), g(quality_ub)]) + "\n")

# --- deletions ---
with open(prefix + "-deletions.txt", "w") as delv:
    delv.write("Reference_Name\tUpstream_Position\tLength\tContext\tCalled\tCount\tTotal\t"
               "Frequency\tPairedUB\n")
    for p in sorted(dcTable):
        for inc in sorted(dcTable[p]):
            count = dcTable[p][inc]
            if count < min_count:
                continue
            pp = p + inc + 1
            called = "TRUE"
            total = span_depth(p, pp)
            freq = count / total if total > 0 else 0
            if freq < min_freq_del or total < min_total:
                called = "FALSE"
            paired_ub = UB(DE, total)
            if freq <= paired_ub:
                called = "FALSE"
            if p < 5:
                left = consensus_seq[0:p + 1]
            else:
                left = consensus_seq[p - 4:p + 1]
            if p > (ref_len - 6 - inc):
                right = consensus_seq[pp:ref_len]
            else:
                right = consensus_seq[pp:pp + 5]
            mid = "-" * inc
            context = left + mid + right
            delv.write("\t".join([ref_name, str(p + 1), g(inc), context, called,
                                  g(count), g(total), g(freq), g(paired_ub)]) + "\n")

# --- phasing inputs (msgpack, read by phase.py) ---
# FIXME legacy perl logic: phasing patterns are extracted by slicing each read's full-length alignment
# string at the variant positions; phase.py then re-slices these by position index. A column matrix of
# allele calls per read (positions x reads) would drop both string round-trips once parity ends.
variant_count = sum(len(v) for v in variants.values())
if variant_count > 1:
    with open(prefix + "-vars.msgpack", "wb") as fh:
        fh.write(msgpack.packb(variants, use_bin_type=True))
    var_positions = sorted(variants)
    read_pats = {}
    for key, c in alignments.items():
        # index the compact "<start>:<body>" directly at each variant position ('.' if uncovered) —
        # equivalent to slicing the full padded string, without materializing it.
        s, _, body = key.partition(":")
        base0 = int(s)
        aln = "".join(body[pos - base0] if 0 <= pos - base0 < len(body) else "." for pos in var_positions)
        if not re.fullmatch(r"[.N]+", aln):
            read_pats[aln] = read_pats.get(aln, 0) + c
    with open(prefix + "-pats.msgpack", "wb") as fh:
        fh.write(msgpack.packb(read_pats, use_bin_type=True))

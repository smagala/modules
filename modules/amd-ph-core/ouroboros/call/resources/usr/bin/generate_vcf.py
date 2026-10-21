#!/usr/bin/env python3
"""generate_vcf.py — write a spec-compliant VCF from the caller's variant tables.

Reads the consensus FASTA plus call.py's -allAlleles.txt / -insertions.txt / -deletions.txt and emits a
VCFv4.2 via pysam.VariantFile (a real library builds the header and serializes records — no hand-printed
text). One normalized biallelic record per variant (SNV / insertion / deletion), left-anchored on the
consensus base per the VCF indel convention.

Schema (modernized — clean INFO, standard FILTER names):
  INFO   DP  total read depth at the locus
         AF  alt allele frequency (alt count / total)
         AD  allele read depths (consensus, alt)
         AQ  mean alt-allele base quality (Phred)
         TYPE  snp | ins | del
  FILTER PASS / LowDepth / LowFreq / LowQual / LowConf / MachineError

A variant needs >= min_count supporting reads to be reported at all. The remaining thresholds become
FILTER annotations: by default only PASS records are written; --print-all writes every candidate with its
FILTER flags so nothing is silently dropped.
"""

import argparse
import sys

import pysam


def read_first_fasta(path):
    """Return (name, uppercase sequence) of the first record."""
    name, parts, seen = None, [], False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if seen:
                    break
                name = line[1:].split()[0]
                seen = True
            elif seen:
                parts.append(line)
    return name, "".join(parts).upper()


def read_table(path):
    """Yield each data row of a tab table as a list of fields (header skipped); [] if absent/empty."""
    try:
        fh = open(path)
    except (OSError, TypeError):
        return
    with fh:
        next(fh, None)
        for line in fh:
            line = line.rstrip("\r\n")
            if line:
                yield line.split("\t")


def build_header(chrom, ref_len, ref_path, argv):
    h = pysam.VariantHeader()
    h.add_line(f"##source=generate_vcf.py {' '.join(argv)}")
    h.add_line(f"##reference={ref_path}")
    h.contigs.add(chrom, length=ref_len)
    h.info.add("DP", 1, "Integer", "Total read depth at the locus")
    h.info.add("AF", "A", "Float", "Alt allele frequency (alt count / total)")
    h.info.add("AD", "R", "Integer", "Allele read depths (consensus, alt)")
    h.info.add("AQ", "A", "Float", "Mean alt-allele base quality (Phred)")
    h.info.add("TYPE", "A", "String", "Variant type (snp, ins, del)")
    h.filters.add("LowDepth", None, None, "Total depth below the minimum")
    h.filters.add("LowFreq", None, None, "Alt frequency below the minimum")
    h.filters.add("LowQual", None, None, "Mean alt-allele quality below the minimum")
    h.filters.add("LowConf", None, None, "Confidence-not-machine-error below the minimum")
    h.filters.add("MachineError", None, None, "Frequency not significant above the per-read error bound")
    return h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference", help="consensus FASTA (the coordinate frame / REF bases)")
    ap.add_argument("alleles", help="call.py -allAlleles.txt")
    ap.add_argument("insertions", help="call.py -insertions.txt")
    ap.add_argument("deletions", help="call.py -deletions.txt")
    ap.add_argument("-G", "--no-gap-allele", action="store_true", dest="no_gap")
    ap.add_argument("-F", "--min-freq", type=float, default=0.005, dest="min_freq")
    ap.add_argument("-I", "--min-insertion-freq", type=float, default=None, dest="min_freq_ins")
    ap.add_argument("-D", "--min-deletion-freq", type=float, default=None, dest="min_freq_del")
    ap.add_argument("-C", "--min-count", type=int, default=2, dest="min_count")
    ap.add_argument("-Q", "--min-quality", type=float, default=20, dest="min_quality")
    ap.add_argument("-T", "--min-total-col-coverage", type=int, default=100, dest="min_total")
    ap.add_argument("-M", "--conf-not-mac-err", type=float, default=0.5, dest="min_conf")
    ap.add_argument("-S", "--sig-level", type=float, default=None, dest="sig_level")
    ap.add_argument("-A", "--auto-min-freq", action="store_true", dest="auto_freq")
    ap.add_argument("-P", "--print-all-vars", action="store_true", dest="print_all")
    ap.add_argument("-N", "--name", default=None, help="unused; accepted for call-options compatibility")
    ap.add_argument("-E", "--paired-error", default=None, help="unused; accepted for compatibility")
    ap.add_argument("-o", "--out", default="-", help="output VCF (default stdout)")
    a = ap.parse_args()

    min_count = max(a.min_count, 0)
    min_freq = max(a.min_freq, 0.0)
    min_freq_ins = min_freq if a.min_freq_ins is None else max(a.min_freq_ins, 0.0)
    min_freq_del = min_freq if a.min_freq_del is None else max(a.min_freq_del, 0.0)
    min_conf = max(a.min_conf, 0.0)
    min_quality = max(a.min_quality, 0.0)
    min_total = a.min_total if a.min_total >= 0 else 100
    take_sig = a.sig_level is not None

    chrom, seq = read_first_fasta(a.reference)
    if not seq:
        sys.exit(f"{sys.argv[0]} ERROR: no reference sequence in {a.reference}")
    ref_len = len(seq)

    def ref_base(pos1):  # 1-based -> consensus base, or None if out of range / not a clean base
        if 1 <= pos1 <= ref_len and seq[pos1 - 1] in "ACGT":
            return seq[pos1 - 1]
        return None

    # AUTO_F heuristic: raise the SNV frequency floor to the highest error-indistinguishable minority freq.
    if a.auto_freq:
        heur = 0.0
        for f in read_table(a.alleles):
            # ref_name, pos, allele, count, total, freq, avg_q, conf, paired_ub, quality_ub, allele_type
            if len(f) >= 11 and f[10] not in ("Consensus", "Majority", "Plurality") and f[2] != "-":
                conf = float(f[7]) if f[7] not in ("", "NA") else 0.0
                freq = float(f[5])
                if conf == 0.0 and freq > heur:
                    heur = freq
        min_freq = max(min_freq, heur)

    records = []  # (start, [(ref, alt, qual, info, filters), ...]) collected then position-sorted

    def filters_for(total, freq, qual, conf, paired_ub, quality_ub, freq_floor, is_del):
        f = []
        if total < min_total:
            f.append("LowDepth")
        if freq < freq_floor:
            f.append("LowFreq")
        if not is_del and qual is not None and qual < min_quality:
            f.append("LowQual")
        if not is_del and conf is not None and conf < min_conf:
            f.append("LowConf")
        if take_sig:
            bad = (paired_ub is not None and freq <= paired_ub) or (
                not is_del and quality_ub is not None and freq <= quality_ub
            )
            if bad:
                f.append("MachineError")
        return f

    def num(x, cast=float):
        return None if x in ("", "NA", None) else cast(x)

    # --- SNVs (minority alleles from allAlleles) ---
    for f in read_table(a.alleles):
        if len(f) < 11 or f[10] in ("Consensus", "Majority", "Plurality"):
            continue
        allele = f[2]
        if allele == "-" or allele == "N":
            continue
        count, total = int(float(f[3])), int(float(f[4]))
        if count < min_count:
            continue
        pos = int(f[1])
        ref = ref_base(pos)
        if ref is None:
            continue
        freq, qual, conf = float(f[5]), num(f[6]), num(f[7])
        paired_ub, quality_ub = num(f[8]), num(f[9])
        flt = filters_for(total, freq, qual, conf, paired_ub, quality_ub, min_freq, is_del=False)
        info = {"DP": total, "AF": (freq,), "AD": (total - count, count), "TYPE": ("snp",)}
        if qual is not None:
            info["AQ"] = (round(qual, 2),)
        records.append((pos - 1, pos, ref, allele, qual, info, flt))

    # --- insertions ---
    for f in read_table(a.insertions):
        # ref_name, upstream_pos, insert, context, called, count, total, freq, avg_q, conf, paired_ub, quality_ub
        if len(f) < 12:
            continue
        count, total = int(float(f[5])), int(float(f[6]))
        if count < min_count:
            continue
        up = int(f[1])
        ref = ref_base(up)
        if ref is None:
            continue
        ins = f[2].upper()
        freq, qual, conf = float(f[7]), num(f[8]), num(f[9])
        paired_ub, quality_ub = num(f[10]), num(f[11])
        flt = filters_for(total, freq, qual, conf, paired_ub, quality_ub, min_freq_ins, is_del=False)
        info = {"DP": total, "AF": (freq,), "AD": (total - count, count), "TYPE": ("ins",)}
        if qual is not None:
            info["AQ"] = (round(qual, 2),)
        records.append((up - 1, up, ref, ref + ins, qual, info, flt))

    # --- deletions ---
    for f in read_table(a.deletions):
        # ref_name, upstream_pos, length, context, called, count, total, freq, paired_ub
        if len(f) < 9:
            continue
        count, total = int(float(f[5])), int(float(f[6]))
        if count < min_count:
            continue
        up, length = int(f[1]), int(f[2])
        anchor = ref_base(up)
        deleted = seq[up : up + length]  # 1-based up+1..up+length
        if anchor is None or len(deleted) != length or any(b not in "ACGT" for b in deleted):
            continue
        freq = float(f[7])
        paired_ub = num(f[8])
        flt = filters_for(total, freq, None, None, paired_ub, None, min_freq_del, is_del=True)
        info = {"DP": total, "AF": (freq,), "AD": (total - count, count), "TYPE": ("del",)}
        records.append((up - 1, up + length, anchor + deleted, anchor, None, info, flt))

    # --- write, position-sorted ---
    header = build_header(chrom, ref_len, a.reference, sys.argv[1:])
    vf = pysam.VariantFile(a.out, "w", header=header)
    for start, stop, ref, alt, qual, info, flt in sorted(records, key=lambda r: (r[0], r[3])):
        if flt and not a.print_all:
            continue
        rec = vf.new_record(
            contig=chrom,
            start=start,
            stop=stop,
            alleles=(ref, alt),
            qual=(None if qual is None else round(qual, 2)),
            info=info,
        )
        for name in flt if flt else ["PASS"]:
            rec.filter.add(name)
        vf.write(rec)
    vf.close()


if __name__ == "__main__":
    main()

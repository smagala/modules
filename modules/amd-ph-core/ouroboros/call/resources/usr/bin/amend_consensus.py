#!/usr/bin/env python3
"""amend_consensus.py — amend a plurality consensus with called variants.

Scoped to the path the pipeline uses: read the plurality consensus FASTA, fold in the
reference-changing minority alleles from a variants.txt as IUPAC ambiguity codes, and write the
amended consensus. The --a2m-reference, coverage-rewrite, and phase-replace branches are out of scope.

Length model (one column per reference position upstream; length changes happen ONLY here):
  N   — missing coverage (from call.py): kept, so the consensus stays full-length / coordinate-anchored.
  -   — deletion: stripped on output, so a called deletion SHORTENS the consensus.
  ins — insertion: spliced in (lowercase) after its upstream column, so it LENGTHENS the consensus.

Reference-changing rule (`--refchg`):
  fixed  — legacy: a minority allele amends the consensus iff count>=min_count AND
           freq>=min_freq(=MIN_AMBIG) AND total>=min_total. A depth-blind point estimate
           (1/4 and 2500/10000 are both "0.25").
  binom  — depth-aware: amend iff the one-sided Clopper-Pearson LOWER bound on the true
           frequency exceeds min_freq (at --refchg-conf). At high depth LB->freq (same calls,
           NM=0 preserved); at low depth the wide CB suppresses thin-coverage refchg inflation.

Indel folding (`-d`/`-i`, off unless the files are given — beyond default IRMA, which routes indels
to the VCF and leaves the amended consensus length-locked to the reference):
  deletions  — each called deletion (count>=min_count, freq>=min_del_freq, total>=min_total) sets its
               reference columns to '-' (the -d deletion-folding rule).
  insertions — the most-frequent called insertion per site (count>=min_count, freq>=min_ins_freq,
               total>=min_total) is spliced in after its upstream column (faithful to -i).
"""

import argparse
import os
import sys

from scipy import stats

# major+minor allele set -> IUPAC code (sorted, de-duplicated keys)
IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "N": "N",
    "AT": "W",
    "CG": "S",
    "AC": "M",
    "GT": "K",
    "AG": "R",
    "CT": "Y",
    "CGT": "B",
    "AGT": "D",
    "ACT": "H",
    "ACG": "V",
    "ACGT": "N",
}


def encode(nts):
    """Map a string of alleles (major + minors) to a single consensus base/IUPAC code."""
    if len(nts) == 1:
        return nts
    if "-" in nts:
        return "?"
    if any(c not in "ACGT" for c in nts):
        return "N"
    key = "".join(sorted(set(nts)))
    return IUPAC[key]


def clopper_pearson_lower(k, d, conf):
    """One-sided Clopper-Pearson lower bound on the true frequency given k successes of d."""
    if k <= 0 or d <= 0:
        return 0.0
    return float(stats.beta.ppf(1.0 - conf, k, d - k + 1))


def read_first_fasta(path):
    """Return (header, uppercase sequence) of the first FASTA record."""
    header, parts, seen = "", [], False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if seen:
                    break
                header = line[1:]
                seen = True
            elif seen:
                parts.append(line)
    return header, "".join(parts).upper()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference", help="plurality consensus FASTA")
    ap.add_argument("variants", help="IRMA-style variants.txt from call.py")
    ap.add_argument("-C", "--count", type=int, default=2, dest="min_count")
    ap.add_argument(
        "-F",
        "--freq",
        type=float,
        default=0.25,
        dest="min_freq",
        help="minimum minority frequency to amend (MIN_AMBIG, default 0.25)",
    )
    ap.add_argument("-T", "--min-total-depth", type=int, default=100, dest="min_total")
    ap.add_argument("-N", "--name", default=None, help="output header/file name")
    ap.add_argument(
        "-S", "--seg", default=None, help="convert a protein name in the header to a segment number (prot:num,...)"
    )
    ap.add_argument(
        "-H", "--fa-header-suffix", action="store_true", help="append the fasta header as a suffix to --name"
    )
    ap.add_argument("-P", "--prefix", default=".", help="output directory prefix")
    ap.add_argument(
        "--refchg",
        choices=("fixed", "binom"),
        default="fixed",
        help="reference-changing rule: fixed = freq>=MIN_AMBIG point estimate; "
        "binom = Clopper-Pearson lower bound > MIN_AMBIG (depth-aware)",
    )
    ap.add_argument(
        "--refchg-conf", type=float, default=0.95, help="confidence for the binom lower bound (default 0.95)"
    )
    ap.add_argument(
        "-d",
        "--deletion-file",
        default=None,
        dest="del_file",
        help="IRMA-style deletions.txt; called deletions are folded in as '-' (shortens consensus)",
    )
    ap.add_argument(
        "-i",
        "--insertion-file",
        default=None,
        dest="ins_file",
        help="IRMA-style insertions.txt; the most-frequent called insertion per site is "
        "spliced in (lengthens consensus)",
    )
    ap.add_argument(
        "-D",
        "--min-del-freq",
        type=float,
        default=0.5,
        dest="min_freq_del",
        help="minimum deletion frequency to fold into the consensus (default 0.5: a "
        "deletion shortens the consensus only when it is the majority allele)",
    )
    ap.add_argument(
        "-I",
        "--min-ins-freq",
        type=float,
        default=0.5,
        dest="min_freq_ins",
        help="minimum insertion frequency to fold into the consensus (default 0.5)",
    )
    a = ap.parse_args()

    min_total = a.min_total if a.min_total >= 0 else 100

    fa_header, seq = read_first_fasta(a.reference)
    if not seq:
        sys.exit(f"{sys.argv[0]} ERROR: no sequence in {a.reference}")
    seq = list(seq)

    # output header: --name else the fasta header; optional seg-number conversion + header suffix
    out_hdr = a.name if a.name is not None else fa_header
    if a.seg:
        for pair in a.seg.split(","):
            prot, _, numbering = pair.partition(":")
            if prot and prot in fa_header:
                out_hdr += "_" + numbering
                break
    if a.fa_header_suffix and a.name is not None:
        out_hdr += "-" + fa_header

    # collect reference-changing minority alleles per position
    valid_pos = {}  # 1-based pos -> accumulated minor allele string
    with open(a.variants) as fh:
        next(fh, None)  # header
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            # Reference_Name, Position, Total, Consensus_Allele, Minority_Allele,
            # Consensus_Count, Minority_Count, Consensus_Frequency, Minority_Frequency, ...
            pos = int(f[1])
            total = int(f[2])
            allele = f[4].upper()
            count = int(f[6])
            freq = float(f[8])
            if count < a.min_count or total < min_total:
                continue
            if a.refchg == "fixed":
                passes = freq >= a.min_freq
            else:
                passes = clopper_pearson_lower(count, total, a.refchg_conf) > a.min_freq
            if passes:
                valid_pos[pos] = valid_pos.get(pos, "") + allele

    # encode major+minor -> IUPAC at each amended position
    for pos, minors in valid_pos.items():
        p = pos - 1
        if 0 <= p < len(seq):
            seq[p] = encode(seq[p] + minors)

    # fold called deletions: set each deleted reference column to '-' (stripped at output -> consensus
    # shortens). Columns are the `length` bases AFTER the 1-based Upstream_Position, i.e. 0-based
    # [upstream, upstream+length) — matching the -d deletion table and call.py's deletions.txt.
    if a.del_file:
        with open(a.del_file) as fh:
            next(fh, None)  # header
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                f = line.split("\t")
                # Reference_Name, Upstream_Position, Length, Context, Called, Count, Total, Frequency, PairedUB
                upstream = int(f[1])
                length = int(f[2])
                count = int(f[5])
                total = int(f[6])
                freq = float(f[7])
                if count >= a.min_count and freq >= a.min_freq_del and total >= min_total:
                    for p in range(upstream, upstream + length):
                        if 0 <= p < len(seq):
                            seq[p] = "-"

    # collect called insertions: most-frequent insert per upstream 0-based column, spliced after it on
    # output (lengthens consensus). The -i insertion rule (key = Upstream_Position - 1).
    insertions = {}  # 0-based upstream column -> (freq, insert)
    if a.ins_file:
        with open(a.ins_file) as fh:
            next(fh, None)  # header
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                f = line.split("\t")
                # Reference_Name, Upstream_Position, Insert, Context, Called, Count, Total, Frequency, ...
                upstream = int(f[1])
                insert = f[2].upper()
                count = int(f[5])
                total = int(f[6])
                freq = float(f[7])
                if count >= a.min_count and freq >= a.min_freq_ins and total >= min_total:
                    p = upstream - 1
                    if p not in insertions or freq > insertions[p][0]:
                        insertions[p] = (freq, insert)

    os.makedirs(a.prefix, exist_ok=True)
    out_path = os.path.join(a.prefix, out_hdr + ".fa")
    with open(out_path, "w") as out:
        out.write(">" + out_hdr + "\n")
        # emit each column (skipping deletions '-' and missing-sentinel '.'), then splice any insertion
        # after its upstream column as lowercase (provenance: inserted, not reference-aligned).
        chars = []
        for p, b in enumerate(seq):
            if b not in ("-", "."):
                chars.append(b)
            if p in insertions:
                chars.append(insertions[p][1].lower())
        out.write("".join(chars) + "\n")


if __name__ == "__main__":
    main()

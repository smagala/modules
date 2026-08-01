#!/usr/bin/env python3
"""parse_sort_results.py — sort reads into genes/segments from SORT results.

Sort reads into genes/segments from SORT results, applying read-count / read-pattern
minimums, a ban list, and pattern-group maxima. Writes:
  <prefix>.txt            gene <TAB> pattern_count <TAB> read_count  (sorted by count desc)
  <prefix>-<gene>.fa      primary reads for valid genes
  <prefix>-<gene>.fa.2    secondary reads for invalid genes

A read's count is taken from its deflated-pattern ID (C<n>%<reads>%...); reads with no
such tag count as 1 pattern each. A read routes to every valid gene it matched (primary);
if it matched no valid gene, it routes to its single best-scoring matched gene (secondary).
"""
import sys
import re
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("sort_results")
ap.add_argument("match_fasta")
ap.add_argument("prefix")
ap.add_argument("-P", "--pattern-list", default=None)
ap.add_argument("-G", "--ignore-annotations", action="store_true")
ap.add_argument("-C", "--min-read-count", type=int, default=1)
ap.add_argument("-D", "--min-read-patterns", type=int, default=1)
ap.add_argument("-B", "--ban-list", default=None)
ap.add_argument("--secondary-frac", type=float, default=0.0,
                help="In an __ALL__ group, also keep a secondary target whose read count is >= this "
                     "fraction of the primary target's (and >= --secondary-floor), so a co-infecting / "
                     "contaminating target gets its own assembly. 0 = single-best only (default, IRMA parity). "
                     "Superseded by --secondary-min-depth when that is set.")
ap.add_argument("--secondary-floor", type=int, default=0,
                help="Absolute read-count floor a retained secondary must also clear (paired with --secondary-frac).")
ap.add_argument("--secondary-min-depth", type=float, default=0.0,
                help="Coverage-based secondary keep (preferred over --secondary-frac): in an __ALL__ group, keep a "
                     "secondary target when its recruited reads give at least this average depth over ITS OWN "
                     "reference — reads * read_len / ref_len >= this. Answers 'can we assemble this genome?' "
                     "independent of how deep the primary is (a fraction-of-primary gate wrongly drops a real "
                     "secondary as the primary gets deeper). Set to the variant-calling floor (min_column_coverage) "
                     "to only split out targets we could actually call a consensus for. Needs --ref-fasta + --read-len.")
ap.add_argument("--ref-fasta", default=None,
                help="Panel reference FASTA; per-target length (first record per id) for the depth estimate.")
ap.add_argument("--read-len", type=float, default=0.0,
                help="Mean read length (bases), for the --secondary-min-depth estimate.")
a = ap.parse_args()

# per-target reference length (first record per gene id), for the coverage-based secondary keep
reflen = {}
if a.ref_fasta:
    seen_ref = set(); cur = None
    with open(a.ref_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                if cur in seen_ref:
                    cur = None            # keep only the first record per gene id
                else:
                    seen_ref.add(cur); reflen[cur] = 0
            elif cur is not None:
                reflen[cur] += len(line.strip())

min_rc = a.min_read_count if a.min_read_count and a.min_read_count >= 1 else 1
min_rp = a.min_read_patterns if a.min_read_patterns and a.min_read_patterns >= 1 else 1

ANNO = re.compile(r"^([^{]+)\{[^}]*\}$")
CPCT = re.compile(r"C\d+%(\d+)%")


def num(x):
    """Numeric score key matching perl's <=> coercion (undef/non-numeric -> 0)."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


counts = {}   # target -> pattern count (insertion order = first-seen, like perl hash build)
ids = {}      # ID -> {target: score}
rcounts = {}  # target -> read count

with open(a.sort_results) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if len(line) == 0:
            continue
        parts = line.split("\t")
        ID = parts[0]
        target = parts[1] if len(parts) > 1 else ""
        score = parts[2] if len(parts) > 2 else None
        if a.ignore_annotations:
            m = ANNO.match(target)
            if m:
                target = m.group(1)
        counts[target] = counts.get(target, 0) + 1
        ids.setdefault(ID, {})[target] = score
        m = CPCT.search(ID)
        if m:
            rcounts[target] = rcounts.get(target, 0) + int(m.group(1))

for target in list(counts):
    rcounts.setdefault(target, counts[target])

# genes sorted by pattern count desc (ties keep first-seen order, deterministic)
genes = sorted(counts, key=lambda g: counts[g], reverse=True)

valid = {}
with open(a.prefix + ".txt", "w") as out:
    for gene in genes:
        out.write("%s\t%s\t%s\n" % (gene, counts[gene], rcounts[gene]))
        valid[gene] = 1 if (counts[gene] >= min_rp and rcounts[gene] >= min_rc) else 0

# remove explicitly banned genes (substring/regex match, like perl m/$pat/)
if a.ban_list:
    for ban in a.ban_list.split(","):
        pat = re.compile(ban)
        for gene in genes:
            if pat.search(gene):
                valid[gene] = 0

# pattern list: divide into groups, keep only the best gene per group
if a.pattern_list:
    pl = a.pattern_list.replace(" ", "").replace(";", ",").replace("/", "|")
    patterns = pl.split(",")
    if len(patterns) == 1 and patterns[0] == "__ALL__":
        # all genes one group -> keep the single best (genes[0]), plus any secondary that is "real" so a
        # co-infecting / contaminating panel target gets its own assembly instead of being discarded.
        # Preferred test: COVERAGE of the secondary's OWN reference (invariant to primary depth) — keep it
        # when reads*read_len/ref_len >= --secondary-min-depth. Falls back to the fraction-of-primary gate
        # (--secondary-frac) when depth inputs are absent; else single-best only.
        depth_gate = a.secondary_min_depth > 0 and a.read_len > 0
        for g in genes[1:]:
            if depth_gate:
                rl = reflen.get(g)
                est_depth = (rcounts[g] * a.read_len / rl) if rl else 0.0
                if est_depth < a.secondary_min_depth:
                    valid[g] = 0
            elif a.secondary_frac > 0:
                if rcounts[g] < max(a.secondary_frac * rcounts[genes[0]], a.secondary_floor):
                    valid[g] = 0
            else:
                valid[g] = 0
    else:
        for pat in patterns:
            rx = re.compile(pat)
            grp = [g for g in genes if rx.search(g)]
            grp.sort(key=lambda g: counts[g], reverse=True)
            for g in grp[1:]:
                valid[g] = 0

if len(genes) == 0:
    sys.exit("%s ERROR: no classification output found! Aborting.\n" % sys.argv[0])

# open a handle per gene (primary .fa for valid, secondary .fa.2 otherwise)
handles = {}
for gene in genes:
    suffix = ".fa" if valid[gene] > 0 else ".fa.2"
    handles[gene] = open(a.prefix + "-" + gene + suffix, "w")


def fasta_records(path):
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


for header, sequence in fasta_records(a.match_fasta):
    sequence = sequence.lower()
    if len(sequence) == 0:
        continue
    secondary = {}
    primary_written = False
    for gene in ids.get(header, {}):
        if valid.get(gene, 0) > 0:
            primary_written = True
            handles[gene].write(">" + header + "\n" + sequence + "\n")
        else:
            secondary[gene] = ids[header][gene]
    if not primary_written and secondary:
        best = sorted(secondary, key=lambda g: num(secondary[g]), reverse=True)[0]
        handles[best].write(">" + header + "\n" + sequence + "\n")

for fh in handles.values():
    fh.close()

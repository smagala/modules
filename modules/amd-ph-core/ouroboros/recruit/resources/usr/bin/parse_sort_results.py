#!/usr/bin/env python3
"""parse_sort_results.py — sort reads into genes/segments from SORT results.

Sort reads into genes/segments from SORT results, applying read-count / read-pattern
minimums, a ban list, and pattern-group maxima. Writes:
  <prefix>.txt            gene <TAB> pattern_count <TAB> read_count  (sorted by count desc)
  <prefix>-<gene>.fa      primary reads for valid genes
  <prefix>-<gene>.fa.2    secondary reads for invalid genes

A read's count is taken from its deflated-pattern ID (C<n>%<reads>, e.g. C1%2); reads with no
such tag count as 1 each. A read routes to every valid gene it matched (primary);
if it matched no valid gene, it routes to its single best-scoring matched gene (secondary).

--collapse-group (off by default) takes a regex naming references that are the same organism.
Every recruited bin matching it is pooled into one, and the losing bins' reads are re-routed
into the winner rather than going to .fa.2. The winner is whichever candidate covers most of
itself in the MATCH hit file (breadth), with mean coverage breaking ties inside --breadth-tol.
"""

import argparse
import re
import sys


# ── record id -> bin ──────────────────────────────────────────────────────────────
# Duplicated verbatim across the ouroboros scripts: ph-core modules cannot share a
# library. Keep the copies identical; if this changes, change them together.
def load_bin_map(path):
    """Record id -> bin label, read from `>ID bin=BIN` deflines in the reference FASTA.

    Panel record ids are ACCESSIONS; the bin is carried in the defline rather than encoded in the
    id. Deriving a bin by string surgery on the id is what made nine separate bugs silent — an
    accession is not a prefix of a bin name, so a forgotten lookup now fails loudly instead of
    quietly returning nothing.

    A record with no `bin=` tag maps to ITSELF. That is not a fallback: the per-gene references the
    gather loop writes from round 2 onward are already named for their bin, so identity is the
    correct mapping for them. Every id in the FASTA is therefore present, and a target that is
    absent from the map is a real error.
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].split()
            if not parts:
                continue
            rid = parts[0]
            out[rid] = next((p[4:] for p in parts[1:] if p.startswith("bin=")), rid)
    return out


def bin_of(target, binmap):
    """Bin for an alignment target. Raises rather than guessing — see load_bin_map."""
    try:
        return binmap[target]
    except KeyError:
        raise KeyError(
            "target %r is not in the reference FASTA's bin map; the hit file and the reference "
            "are out of step (%d records known)" % (target, len(binmap))
        )


ap = argparse.ArgumentParser()
ap.add_argument("sort_results")
ap.add_argument("match_fasta")
ap.add_argument("prefix")
ap.add_argument("-P", "--pattern-list", default=None)
ap.add_argument("-G", "--ignore-annotations", action="store_true")
ap.add_argument("-C", "--min-read-count", type=int, default=1)
ap.add_argument("-D", "--min-read-patterns", type=int, default=1)
ap.add_argument("-B", "--ban-list", default=None)
ap.add_argument(
    "--secondary-frac",
    type=float,
    default=0.0,
    help="In an __ALL__ group, also keep a secondary target whose read count is >= this "
    "fraction of the primary target's (and >= --secondary-floor), so a co-infecting / "
    "contaminating target gets its own assembly. 0 = single-best only (default, IRMA parity). "
    "Superseded by --secondary-min-depth when that is set.",
)
ap.add_argument(
    "--secondary-floor",
    type=int,
    default=0,
    help="Absolute read-count floor a retained secondary must also clear (paired with --secondary-frac).",
)
ap.add_argument(
    "--secondary-min-depth",
    type=float,
    default=0.0,
    help="Coverage-based secondary keep (preferred over --secondary-frac): in an __ALL__ group, keep a "
    "secondary target when its recruited reads give at least this average depth over ITS OWN "
    "reference — reads * read_len / ref_len >= this. Answers 'can we assemble this genome?' "
    "independent of how deep the primary is (a fraction-of-primary gate wrongly drops a real "
    "secondary as the primary gets deeper). Set to the variant-calling floor (min_column_coverage) "
    "to only split out targets we could actually call a consensus for. Needs --ref-fasta + --read-len.",
)
ap.add_argument(
    "--ref-fasta",
    default=None,
    help="Panel reference FASTA; per-target length (first record per id) for the depth estimate.",
)
ap.add_argument(
    "--read-len", type=float, default=0.0, help="Mean read length (bases), for the --secondary-min-depth estimate."
)
ap.add_argument(
    "--collapse-group",
    action="append",
    default=None,
    metavar="REGEX",
    help="Pool every recruited target matching REGEX into one bin, keeping the winner "
    "and re-routing the losers' reads into it. Repeatable. Needs --hits and "
    "--ref-fasta. Off by default.",
)
ap.add_argument(
    "--hits",
    default=None,
    help="MATCH hit file (rammap PAF or blastn outfmt 6) used to measure candidate "
    "breadth/depth for --collapse-group.",
)
ap.add_argument("--hits-format", choices=("paf", "blast6"), default="paf", help="Format of --hits. Default paf.")
ap.add_argument(
    "--hits-min-len",
    type=int,
    default=0,
    help="Significance floor on alignment block length for breadth/depth. Set to the "
    "same value MATCH used so the measurement matches recruitment.",
)
ap.add_argument(
    "--hits-min-pid", type=float, default=0.0, help="Significance floor on percent identity for breadth/depth."
)
ap.add_argument(
    "--breadth-min-depth",
    type=int,
    default=1,
    help="A reference position counts as covered at this depth or above. Default 1.",
)
ap.add_argument(
    "--breadth-tol",
    type=float,
    default=0.02,
    # DEPRECATED and inert: the winner is chosen on identity, not on a breadth tie-break.
    # Still accepted so existing callers do not break.
    help="Candidates within this much breadth of the leader are treated as tied and "
    "decided on depth instead. Default 0.02.",
)
ap.add_argument(
    "--collapse-min-breadth-frac",
    type=float,
    default=0.5,
    help="A candidate must cover at least this fraction of the BEST candidate's "
    "breadth to be eligible to win. Breadth gates; identity decides. Default 0.5.",
)
ap.add_argument(
    "--collapse-warn-breadth",
    type=float,
    default=0.50,
    help="Treat a losing candidate as a suspected co-infection when it covers at least "
    "this much of itself. A broadly covered sibling is the signature of a second "
    "organism; reads split off ONE population leave the sibling partially covered. "
    "Default 0.50.",
)
ap.add_argument(
    "--collapse-on-suspect",
    choices=("skip", "warn"),
    default="skip",
    help="What to do when a group trips --collapse-warn-breadth. 'skip' (default) "
    "leaves that group uncollapsed, so both organisms keep their own bin; 'warn' "
    "collapses anyway and only reports it. Default skip.",
)
ap.add_argument(
    "--collapse-report", default=None, help="Write a TSV of every candidate's breadth/depth and the collapse decision."
)
a = ap.parse_args()

# per-target reference length (first record per gene id), for the coverage-based secondary keep
# per-target reference length, keyed by BIN. The bin comes from the defline's `bin=` tag via
# load_bin_map, never from parsing the record id — see that function for why.
BINMAP = load_bin_map(a.ref_fasta) if a.ref_fasta else {}
reflen = {}
if a.ref_fasta:
    seen_ref = set()
    cur = None
    with open(a.ref_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                rid = line[1:].split()[0]
                cur = BINMAP.get(rid, rid)
                if cur in seen_ref:
                    cur = None  # keep only the first record per bin
                else:
                    seen_ref.add(cur)
                    reflen[cur] = 0
            elif cur is not None:
                reflen[cur] += len(line.strip())

min_rc = a.min_read_count if a.min_read_count and a.min_read_count >= 1 else 1
min_rp = a.min_read_patterns if a.min_read_patterns and a.min_read_patterns >= 1 else 1

ANNO = re.compile(r"^([^{]+)\{[^}]*\}$")
# Deflated pattern IDs are "C<id>%<reads>" (e.g. C1%2). The trailing "%" the old pattern required
# never appears, so this never matched and every read count silently fell back to the PATTERN count
# -- understating recruited reads ~2x, and with it the --secondary-min-depth coverage estimate.
CPCT = re.compile(r"C\d+%(\d+)")


def num(x):
    """Numeric score key matching perl's <=> coercion (undef/non-numeric -> 0)."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


counts = {}  # target -> pattern count (insertion order = first-seen, like perl hash build)
ids = {}  # ID -> {target: score}
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


def measure_targets(path, fmt, min_len, min_pid, lengths):
    """Per-target (breadth, mean_depth, identity) from a MATCH hit file, weighted by read count.

    Every significant hit counts, including hits the target did not win: the question is what a
    target would cover if it took the whole pool, which is what the winner is about to be given.

    IDENTITY is the one that picks the winner. Breadth is a fraction of each reference's OWN
    length, so it is not comparable across references that differ in length: 445/892 (0.499) loses
    to 502/984 (0.510) on near-equal absolute coverage, which is how 10A_S1 — a snowshoe hare virus
    sample whose L segment is khatangaense at 100% — had its S segment labelled La Crosse. Identity
    is length-free and measures the thing actually being asked: which reference do these reads look
    like. Summed over hits as total matches / total aligned block, so it is read-weighted.
    """
    cov = {}
    idc = {}  # target -> [sum matched bases, sum aligned block], for read-weighted identity
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if fmt == "paf":
            if len(f) < 11:
                continue
            qname, target = f[0], f[5]
            tstart, tend = int(f[7]), int(f[8])
            matches, block = int(f[9]), int(f[10])
        else:
            # blastn outfmt 6: qseqid sseqid pident length qstart qend sstart send sstrand bitscore
            if len(f) < 8:
                continue
            qname, target = f[0], f[1]
            block = int(f[3])
            matches = int(round(float(f[2]) * block / 100.0))
            sstart, send = int(f[6]), int(f[7])
            tstart, tend = (sstart - 1, send) if sstart <= send else (send - 1, sstart)
        if block < min_len or block <= 0:
            continue
        if 100.0 * matches / block < min_pid:
            continue
        # Hit targets are RECORD ids (accessions in round 1, bin-named refined refs later);
        # `lengths` is keyed by BIN. Round 1 previously missed every lookup here, so nothing was
        # measured and --collapse-group silently found fewer than two candidates — collapsing only
        # ever pooled the round-2 carry residual. The map makes that a hard error instead.
        target = bin_of(target, BINMAP)
        ref_len = lengths.get(target)
        if not ref_len:
            continue
        m = CPCT.search(qname)
        weight = int(m.group(1)) if m else 1
        ident = idc.setdefault(target, [0, 0])
        ident[0] += matches * weight
        ident[1] += block * weight
        arr = cov.get(target)
        if arr is None:
            arr = cov[target] = [0] * ref_len
        lo = max(0, tstart)
        hi = min(ref_len, tend)
        for i in range(lo, hi):
            arr[i] += weight
    out = {}
    for target, arr in cov.items():
        n = len(arr)
        covered = sum(1 for c in arr if c >= a.breadth_min_depth)
        m, b = idc.get(target, (0, 0))
        out[target] = (covered / n, sum(arr) / n, (100.0 * m / b) if b else 0.0)
    return out


collapse_to = {}  # losing target -> winning target
collapse_rows = []
if a.collapse_group:
    if not a.hits or not reflen:
        sys.exit("%s ERROR: --collapse-group needs --hits and --ref-fasta.\n" % sys.argv[0])
    measured = measure_targets(a.hits, a.hits_format, a.hits_min_len, a.hits_min_pid, reflen)
    for group in a.collapse_group:
        rx = re.compile(group)
        # An unrecruited panel record has no bin to collapse and no reads to donate.
        cands = [g for g in genes if rx.search(g) and g in measured]
        if len(cands) < 2:
            continue
        breadth = {g: measured[g][0] for g in cands}
        depth = {g: measured[g][1] for g in cands}
        ident = {g: measured[g][2] for g in cands}
        # Highest IDENTITY wins, among candidates that cover enough of themselves to be judged.
        # Breadth is the eligibility floor, not the ranking key — see measure_targets. The floor
        # stops a short, highly conserved fragment out-scoring a reference the reads genuinely
        # span; it is relative to the best-covered candidate so it adapts to shallow samples.
        floor = max(breadth.values()) * a.collapse_min_breadth_frac
        eligible = [g for g in cands if breadth[g] >= floor] or list(cands)
        winner = max(eligible, key=lambda g: (ident[g], breadth[g], depth[g], g))
        suspects = [g for g in cands if g != winner and breadth[g] >= a.collapse_warn_breadth]
        for g in suspects:
            sys.stderr.write(
                "%s WARNING: %s covers %.3f of itself against the winner %s at %.3f — a broadly "
                "covered sibling is a suspected co-infection, not a split.%s\n"
                % (
                    sys.argv[0],
                    g,
                    breadth[g],
                    winner,
                    breadth[winner],
                    " Group left uncollapsed." if a.collapse_on_suspect == "skip" else "",
                )
            )
        # One suspect disqualifies the whole group: merging the rest would still fold a second
        # organism's reads into whichever bin survives.
        hold = bool(suspects) and a.collapse_on_suspect == "skip"
        if not hold:
            # POOL FIRST, THEN GATE. Collapsing hands the group's whole read pool to the winner, so
            # the winner must be judged on that pooled total — not waved through. Setting
            # valid[winner] = 1 unconditionally resurrected bins the secondary depth gate had
            # already, correctly, rejected: 10G_S7 reported a "lacrosseense_S genome" built from
            # FOUR reads beside a complete Jamestown Canyon genome at 1.1M reads per segment,
            # because the collapse block runs after the gate and overrode it.
            pooled = sum(rcounts.get(g, 0) for g in cands)
            rl = reflen.get(winner)
            if a.secondary_min_depth > 0 and a.read_len > 0 and rl:
                est = pooled * a.read_len / rl
                valid[winner] = 1 if est >= a.secondary_min_depth else 0
                if not valid[winner]:
                    sys.stderr.write(
                        "%s NOTE: collapsed group led by %s pools %d reads for ~%.1fx over %d bp, "
                        "below --secondary-min-depth %.1f; not reported.\n"
                        % (sys.argv[0], winner, pooled, est, rl, a.secondary_min_depth)
                    )
            else:
                valid[winner] = 1
        for g in cands:
            suspect = g in suspects
            if g == winner:
                role = "winner_not_collapsed" if hold else "winner"
            elif hold:
                role = "suspect_not_collapsed" if suspect else "held_not_collapsed"
            else:
                valid[g] = 0
                collapse_to[g] = winner
                role = "collapsed_suspect" if suspect else "collapsed"
            collapse_rows.append(
                (group, g, reflen.get(g, 0), counts.get(g, 0), rcounts.get(g, 0), breadth[g], depth[g], ident[g], role)
            )

# Only rounds that actually collapsed something get a report; later rounds have one candidate
# left and would otherwise emit a header-only file per round.
if a.collapse_report and collapse_rows:
    with open(a.collapse_report, "w") as rfh:
        rfh.write("group\ttarget\tref_len\tpatterns\treads\tbreadth\tmean_depth\tidentity\trole\n")
        for row in collapse_rows:
            rfh.write("%s\t%s\t%d\t%d\t%d\t%.6f\t%.4f\t%.3f\t%s\n" % row)

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
    written = set()
    for gene in ids.get(header, {}):
        # Two of a read's targets can forward to the same winner, so writes are deduplicated.
        dest = gene if valid.get(gene, 0) > 0 else collapse_to.get(gene)
        if dest is not None and valid.get(dest, 0) > 0:
            primary_written = True
            if dest not in written:
                written.add(dest)
                handles[dest].write(">" + header + "\n" + sequence + "\n")
        else:
            secondary[gene] = ids[header][gene]
    if not primary_written and secondary:
        best = sorted(secondary, key=lambda g: num(secondary[g]), reverse=True)[0]
        handles[best].write(">" + header + "\n" + sequence + "\n")

for fh in handles.values():
    fh.close()

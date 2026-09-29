#!/usr/bin/env python3
"""genome_report.py — per-sample genomes / contamination screen.

For each assembled genome (target) of a sample, tabulate what it is and how real it is:
recruited reads, consensus length / % called, and — the key contamination discriminator —
IDENTITY of the (reference-pinned) consensus to ITS OWN panel reference.

A real genome resolves to its own reference at high identity over most of its length; a
conserved-region cross-map echo of the primary aligns only partially and at ~between-species
identity. So `identity_to_ref` + `aligned_frac` separate genuine co-infection / contamination
from mis-sorted primary reads, per target, at a glance.

Identity is measured by walking the consensus-vs-own-reference alignment and comparing actual
bases against the reference — indel-aware (drops through insertions/deletions), and it skips N
(no-coverage) positions, so it reports identity over CALLED bases regardless of dropout. Use the
REFERENCE-PINNED (stitched/full-length) consensus, not a deletion-truncated one, so the alignment
stays in the reference frame.

Inputs: per-target consensus FASTAs, the panel reference FASTA (for reference bases + lengths),
and a directory of per-target consensus-vs-own-ref SAMs (<target>.sam). Optionally the gather
sorted_read_stats.txt for recruited read counts. Writes a TSV to stdout.
"""
import argparse
import os
import pysam



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
            "are out of step (%d records known)" % (target, len(binmap)))

BINMAP = {}          # set from --ref-fasta in main(); record id -> bin


def gene_of(record_id):
    """Bin for a panel record id. Looked up, never parsed — see load_bin_map."""
    return BINMAP.get(record_id, record_id)


def read_fastas(path):
    """Return (by_bin, by_record).

    `by_bin` keeps the first record per bin, for callers that just need a representative length.
    `by_record` keeps every record under its full id, because identity has to be measured against
    the variant the assembly was actually built against. A bin holds competing variants precisely
    because one fits the sample and the others do not: in khatangaense_L, NC_055196 is 81.7% from
    these samples while MK352484 is 100.0%. Measuring against the first record of the bin made
    25 of 85 rows look like sub-90% "nearest neighbour" calls when only 5 are.

    Keying on the raw record id ALONE is the older bug — a consensus is named for its bin, so a
    bin-keyed lookup is still needed. Both are returned rather than choosing one.
    """
    by_record, order = {}, []
    cur = None
    for line in open(path):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            by_record[cur] = []
            order.append(cur)
        elif cur is not None:
            by_record[cur].append(line.strip())
    by_record = {k: "".join(v).upper() for k, v in by_record.items()}
    by_bin = {}
    for rid in order:
        by_bin.setdefault(gene_of(rid), by_record[rid])
    return by_bin, by_record


def target_of(path):
    """Target name from the FASTA header (robust to <gene>.stitched.fa / <gene>.fa naming)."""
    for line in open(path):
        if line.startswith(">"):
            return line[1:].split()[0]
    return os.path.basename(path).split(".")[0]


def consensus_stats(path):
    s = "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
    return len(s), (len(s) - s.count("N"))


def identity_to_ref(sam_path, by_record, bin_name):
    """Aggregate ALL alignments of the consensus (a dropout consensus fragments across N-gaps).

    Returns (identity_pct over called bases, aligned_frac, called_ref_positions, record_id, ref_len)

      aligned_frac   reference positions reached by ANY consensus base, N included — breadth.
      called_ref_pos reference positions carrying a CALLED (non-N) base — the set that answers
                     "how much of this genome did we actually recover".
      record_id      WHICH panel variant these are measured against.

    Scored PER PANEL RECORD, best one wins. A bin holds competing variants precisely because one
    fits the sample and the others do not, and the assembly was built against whichever won — so
    scoring against the bin's first record scores the wrong sequence. In khatangaense_L, NC_055196
    is 81.7% from these samples while MK352484 is 100.0%; in West Nile, NC_009942 is 78.8% while
    EF429197 is 99.8%. Against first-record-only, 25 of 85 rows looked like
    sub-90% "nearest neighbour" calls when only 5 are.

    The called-position set is also the point: counting called bases and dividing by the reference
    length is a different question and is not even bounded by it. Refinement can extend a consensus
    past its panel record — the snowshoe hare L segment assembles to 6,975 bp of which 6,966 are
    called against a 6,882 bp reference, which divides out to 101.2% — because insertions add
    consensus bases corresponding to no reference position. Going through the alignment, an
    inserted base lands on rpos None and is skipped, so the count cannot exceed the reference.
    """
    if not os.path.exists(sam_path) or os.path.getsize(sam_path) == 0:
        return None, 0.0, set(), "", 0
    # every panel record belonging to this bin is a candidate
    cands = {rid: seq for rid, seq in by_record.items() if gene_of(rid) == bin_name}
    if not cands:
        return None, 0.0, set(), "", 0
    stats = {rid: {"called": 0, "matches": 0, "covered": set(), "cpos": set()} for rid in cands}
    save = pysam.set_verbosity(0)
    # check_sq=False, and tolerate an unopenable SAM. A target whose consensus-vs-reference
    # alignment produced nothing leaves a SAM with no @SQ lines, and pysam raises rather than
    # yielding zero reads: "file has no sequences defined". This was latent until the panel lookup
    # was keyed by bin instead of record id — before that `refseq` was always empty and the guard
    # above returned first, so the SAM was never opened on ANY row. Fixing the lookup turned a
    # silently-empty column into a crash that failed the whole GENOMEREPORT task for the sample.
    try:
        sf = pysam.AlignmentFile(sam_path, "r", check_sq=False)
    except (ValueError, OSError):
        pysam.set_verbosity(save)
        return None, 0.0, set(), "", 0
    with sf:
        for read in sf.fetch(until_eof=True):
            if read.is_unmapped or read.query_sequence is None:
                continue
            rid = read.reference_name
            if rid not in cands:
                continue
            refseq, st = cands[rid], stats[rid]
            q = read.query_sequence.upper()
            for qpos, rpos in read.get_aligned_pairs():
                if qpos is None or rpos is None or rpos >= len(refseq):
                    continue  # insertion / deletion column, or past ref — not a comparable base
                st["covered"].add(rpos)
                qb = q[qpos]
                if qb in "ACGT":               # skip N (no-coverage) — identity over CALLED bases
                    st["called"] += 1
                    st["cpos"].add(rpos)
                    if qb == refseq[rpos]:
                        st["matches"] += 1
    pysam.set_verbosity(save)

    scored = [(rid, s) for rid, s in stats.items() if s["called"]]
    if not scored:
        return None, 0.0, set(), "", len(cands[sorted(cands)[0]])
    # best = highest identity, tie-broken by how much of the reference it actually covers, so a
    # short high-identity fragment cannot beat a full-length alignment one point behind it
    rid, s = max(scored, key=lambda kv: (kv[1]["matches"] / kv[1]["called"], len(kv[1]["cpos"])))
    reflen = len(cands[rid])
    return (100.0 * s["matches"] / s["called"]), (len(s["covered"]) / reflen), s["cpos"], rid, reflen


def load_read_counts(path):
    counts = {}
    if path and os.path.exists(path):
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2].isdigit():
                counts[f[0]] = counts.get(f[0], 0) + int(f[2])
    return counts


def main():
    ap = argparse.ArgumentParser(description="Per-sample genomes / contamination screen with identity-to-ref.")
    ap.add_argument("consensus", nargs="+", help="per-target (reference-pinned) consensus FASTAs")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--ref-fasta", required=True, help="panel reference FASTA (reference bases per target)")
    ap.add_argument("--sam-dir", default=".", help="dir of <target>.sam (consensus vs own ref)")
    ap.add_argument("--sort-stats", default=None, help="gather sorted_read_stats.txt for read counts")
    ap.add_argument("--primary", default=None, help="primary target name (default: most reads)")
    args = ap.parse_args()

    global BINMAP
    BINMAP = load_bin_map(args.ref_fasta)
    refs, by_record = read_fastas(args.ref_fasta)
    reads = load_read_counts(args.sort_stats)
    rows = []
    for cons in args.consensus:
        target = target_of(cons)
        clen, called = consensus_stats(cons)
        bin_name = gene_of(target)
        idpct, aln_frac, called_ref_pos, matched, matched_len = identity_to_ref(
            os.path.join(args.sam_dir, target + ".sam"), by_record, bin_name)
        # length of the variant the assembly actually matched; the bin's first record only as a
        # fallback when nothing aligned, so pct_called still has a denominator
        reflen = matched_len or len(refs.get(bin_name, ""))
        # % called is how much of the REFERENCE carries a called base, taken THROUGH THE ALIGNMENT.
        # Two ways to get this wrong, and we shipped both. Dividing by the consensus length reports
        # every trimmed assembly as 100% complete — polish trims back to covered columns, so 16,552
        # reads yielding 72 bp of a 4,458 bp M segment read as "100.0% called". Dividing by the
        # reference length instead fixes that but is not bounded by it, because insertions add
        # consensus bases belonging to no reference position: a 6,966-called consensus over a 6,882
        # bp reference reported 101.2%. Only the alignment answers the actual question. Fall back to
        # the naive fraction when there is no SAM or no reference to align against.
        if reflen and called_ref_pos:
            pct_called = 100.0 * len(called_ref_pos) / reflen
        else:
            denom = reflen or clen
            pct_called = (100.0 * called / denom) if denom else 0.0
        rows.append({
            "target": target, "reads": reads.get(target, 0), "cons_len": clen, "ref_len": reflen,
            "called": called, "ref_covered": len(called_ref_pos), "pct_called": pct_called,
            "identity_to_ref": idpct, "aligned_frac": aln_frac,
            "matched_record": matched or "NA",
        })

    if args.primary and any(r["target"] == args.primary for r in rows):
        prim = args.primary
    elif any(r["reads"] for r in rows):
        prim = max(rows, key=lambda r: r["reads"])["target"]
    else:
        prim = max(rows, key=lambda r: r["cons_len"])["target"]

    rows.sort(key=lambda r: (r["target"] != prim, -r["reads"], r["target"]))
    print("sample\tgenome\trole\treads\tcons_len\tcalled\tref_len\tref_covered\tpct_called\t"
          "identity_to_ref\taligned_frac\tmatched_record")
    for r in rows:
        role = "primary" if r["target"] == prim else "secondary"
        idc = f"{r['identity_to_ref']:.1f}" if r["identity_to_ref"] is not None else "NA"
        print(f"{args.sample}\t{r['target']}\t{role}\t{r['reads']}\t{r['cons_len']}\t{r['called']}\t"
              f"{r['ref_len']}\t{r['ref_covered']}\t{r['pct_called']:.1f}\t{idc}\t"
              f"{r['aligned_frac']:.2f}\t{r['matched_record']}")


if __name__ == "__main__":
    main()

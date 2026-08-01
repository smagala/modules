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


def read_fastas(path):
    seqs = {}
    cur = None
    for line in open(path):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif cur is not None:
            seqs[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def target_of(path):
    """Target name from the FASTA header (robust to <gene>.stitched.fa / <gene>.fa naming)."""
    for line in open(path):
        if line.startswith(">"):
            return line[1:].split()[0]
    return os.path.basename(path).split(".")[0]


def consensus_stats(path):
    s = "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
    return len(s), (len(s) - s.count("N"))


def identity_to_ref(sam_path, refseq):
    """Aggregate ALL alignments of the consensus (a dropout consensus fragments across N-gaps);
    return (identity_pct over called bases, aligned_frac = fraction of the reference covered)."""
    if not os.path.exists(sam_path) or not refseq:
        return None, 0.0
    called = matches = 0
    covered = set()
    save = pysam.set_verbosity(0)
    with pysam.AlignmentFile(sam_path, "r") as sf:
        for read in sf.fetch(until_eof=True):
            if read.is_unmapped or read.query_sequence is None:
                continue
            q = read.query_sequence.upper()
            for qpos, rpos in read.get_aligned_pairs():
                if qpos is None or rpos is None or rpos >= len(refseq):
                    continue  # insertion / deletion column, or past ref — not a comparable base
                covered.add(rpos)
                qb = q[qpos]
                if qb in "ACGT":               # skip N (no-coverage) — identity over CALLED bases
                    called += 1
                    if qb == refseq[rpos]:
                        matches += 1
    pysam.set_verbosity(save)
    if called == 0:
        return None, 0.0
    return (100.0 * matches / called), (len(covered) / len(refseq) if refseq else 0.0)


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

    refs = read_fastas(args.ref_fasta)
    reads = load_read_counts(args.sort_stats)
    rows = []
    for cons in args.consensus:
        target = target_of(cons)
        clen, called = consensus_stats(cons)
        idpct, aln_frac = identity_to_ref(os.path.join(args.sam_dir, target + ".sam"), refs.get(target, ""))
        rows.append({
            "target": target, "reads": reads.get(target, 0), "cons_len": clen,
            "pct_called": (100.0 * called / clen) if clen else 0.0,
            "identity_to_ref": idpct, "aligned_frac": aln_frac,
        })

    if args.primary and any(r["target"] == args.primary for r in rows):
        prim = args.primary
    elif any(r["reads"] for r in rows):
        prim = max(rows, key=lambda r: r["reads"])["target"]
    else:
        prim = max(rows, key=lambda r: r["cons_len"])["target"]

    rows.sort(key=lambda r: (r["target"] != prim, -r["reads"], r["target"]))
    print("sample\tgenome\trole\treads\tcons_len\tpct_called\tidentity_to_ref\taligned_frac")
    for r in rows:
        role = "primary" if r["target"] == prim else "secondary"
        idc = f"{r['identity_to_ref']:.1f}" if r["identity_to_ref"] is not None else "NA"
        print(f"{args.sample}\t{r['target']}\t{role}\t{r['reads']}\t{r['cons_len']}\t"
              f"{r['pct_called']:.1f}\t{idc}\t{r['aligned_frac']:.2f}")


if __name__ == "__main__":
    main()

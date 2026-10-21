#!/usr/bin/env python3
"""select_reference.py — pull one bin's reference record out of a multi-record panel FASTA.

Round 1 of the recruit loop has no per-gene refined reference yet, so the ALIGN->REFINE step used
to be handed the WHOLE panel as its reference. That is wrong in two ways:

  * `sam_align_stats.py` reads only the FIRST record of the file for `ref_len` (used for length and
    gap repair), so the refined reference was sized by whichever record happened to sort first in
    the panel — not by the bin being refined. Reordering a panel with byte-identical sequences moved
    `pct_called` by 12 points.
  * the bin's reads were aligned against every panel record, so a read could land on a neighbour and
    contribute at coordinates from a different genome. Measured at ~0.2% of alignments, small but
    not zero.

Both go away if round 1 aligns against the bin's own record. The bin label IS the record ID
(ADR-0014), so the selection is exact rather than a fuzzy name match.

Matching mirrors rammap_partition.py: resolve each record ID to its BIN via the reference FASTA's
`bin=` deflines, because the bin is the form the class labels — and therefore the gene names —
carry. The bin is looked up, never parsed out of the ID.
"""

import argparse
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
ap.add_argument("--fasta", required=True, help="Panel FASTA to select from")
ap.add_argument("--name", required=True, help="Gene / bin label to select (a record ID)")
ap.add_argument("-o", "--out", required=True, help="Output FASTA")
ap.add_argument(
    "--min-length-frac",
    type=float,
    default=0.95,
    help="A candidate must be at least this fraction of the bin's MEDIAN record "
    "length to be eligible. Completeness is a hard requirement (refine cannot "
    "create bases the reference lacks) where divergence is not (refine "
    "converges), so length gates and identity only ranks within the gate.",
)
ap.add_argument(
    "--max-length-frac",
    type=float,
    default=1.05,
    help="Upper bound on the same median. Refine cannot delete spurious bases any more "
    "than it can create missing ones: an over-long deposit (a duplicated 3' end, "
    "say) contributes reference nothing can cover, which deflates breadth and "
    "offers mismapping targets.",
)
ap.add_argument(
    "--subtype-counts",
    help="TSV from rammap_partition.py (<out>.subtypes.tsv): bin, record, reads. "
    "When given, the bin's reference is the record its own reads matched most "
    "often instead of the first record in the file. Without it, behaviour is "
    "unchanged.",
)
a = ap.parse_args()


BINMAP = load_bin_map(a.fasta)


def gene_of(header):
    """Bin for a panel defline. Looked up, never parsed out of the id — see load_bin_map."""
    rid = header.split()[0]
    return BINMAP.get(rid, rid)


records, name, buf = [], None, []
with open(a.fasta) as fh:
    for line in fh:
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if name is not None:
                records.append((name, buf))
            name, buf = line[1:], []
        else:
            buf.append(line)
if name is not None:
    records.append((name, buf))

hits = [(h, s) for h, s in records if gene_of(h) == a.name]
if not hits:
    sys.exit(f"{sys.argv[0]} ERROR: no record named {a.name!r} in {a.fasta}\n")

# Several records can share a bin, and which one we align
# against decides how much of the bin's reads can be placed at all. Taking the first — which is the
# LONGEST, because build_panel.py emits longest-first — chooses on completeness, a property with no
# relation to how close the record is to this sample. Dereplication deliberately keeps variants that
# are far apart, so the first record can sit 13% from the reads while a sibling sits at 1%, and the
# alignment then loses most of the bin to a reference nobody chose on merit.
#
# With --subtype-counts the choice is made on the bin's own reads. Without it the old behaviour
# stands, so callers that do not pass evidence are unaffected.
chosen = hits[0]
if a.subtype_counts:
    votes = {}
    with open(a.subtype_counts) as fh:
        header = next(fh, "")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0] == a.name:
                try:
                    votes[parts[1]] = int(parts[2])
                except ValueError:
                    continue
    if votes:
        by_id = {h.split()[0]: (h, s) for h, s in hits}
        # LENGTH IS A GATE, IDENTITY IS A PREFERENCE. The loop refines its reference onto the
        # sample, so a divergent start converges — but it can only edit bases that exist, so a
        # record missing its termini caps what can ever be covered and reports breadth against a
        # short denominator. Choosing on identity alone will happily pick such a record because
        # being truncated does not stop it matching well over the part it does have.
        lens = sorted(len("".join(s)) for _h, s in hits)
        n = len(lens)
        median = lens[n // 2] if n % 2 else (lens[n // 2 - 1] + lens[n // 2]) / 2
        eligible = {
            r for r in by_id if a.min_length_frac * median <= len("".join(by_id[r][1])) <= a.max_length_frac * median
        }
        pool = {r: v for r, v in votes.items() if r in eligible} or votes
        if not eligible:
            sys.stderr.write(
                f"{sys.argv[0]}: WARNING no record in {a.name} falls within "
                f"{a.min_length_frac:.0%}-{a.max_length_frac:.0%} of the bin median "
                f"({median:.0f} bp); choosing on identity alone\n"
            )
        # ties break on record ID so the result does not depend on dict or file order
        winner = max(sorted(pool), key=lambda r: pool[r])
        if winner in by_id:
            chosen = by_id[winner]
            # compare on the record ID: by_id rebuilds the tuple, so an identity test here
            # would report every bin as changed, including the ones that kept the first record
            if winner != hits[0][0].split()[0]:
                sys.stderr.write(
                    f"{sys.argv[0]}: {a.name} -> {winner} on {votes[winner]} reads "
                    f"(first record was {hits[0][0].split()[0]})\n"
                )
        else:
            sys.stderr.write(
                f"{sys.argv[0]}: WARNING {winner} not in {a.fasta}; " f"falling back to the first record\n"
            )

with open(a.out, "w") as out:
    h, s = chosen
    out.write(">" + h + "\n")
    out.write("\n".join(s) + "\n")

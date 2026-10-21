#!/usr/bin/env python3
"""rammap_partition.py — split reads from a rammap PAF into match/class/nomatch/chim.

Chimera detection lives in the aligner: rammap --filter-chimera tags every record of a
chimeric read (significant hits covering one target on both strands) with ch:i:1. This only
routes on that tag plus the recruitment significance threshold — no chimera logic of its own.

Run rammap with: -x sr --secondary yes -N 7 -c --filter-chimera
Significance: block_len >= --min-len and 100*matches/block_len >= --min-pid.
PAF cols used: 1 qname, 6 target, 10 matches, 11 block_len; the ch:i:1 tag marks chimeras.

Writes <out>.match/.class/.nomatch/.chim. Chimeric reads go to .chim (unless --incl-chim),
recruited reads to .match with the best target in .class, the rest to .nomatch.

AMBIGUITY DEFERRAL (--defer-ambiguous, off by default)
------------------------------------------------------
Assignment is winner-take-all on `matches`, and the recruit loop only carries `.nomatch`
forward — so a read handed to the wrong near-neighbour in round 1 is never reconsidered. The
loop can rescue a read it failed to recruit, but not one it mis-recruited.

This guard defers the reads that decision is least able to justify: those where the runner-up
target is within a hair of the winner. A deferred read goes to .nomatch, so the next round
re-competes it against a REFINED reference that should separate the two cleanly.

The test is deliberately RELATIVE — a margin between two hits that are BOTH already
significant — and never absolute. That distinction is the whole safety property:

  * Contaminant / off-target reads do not have a near-tie among panel references. They hit one
    reference weakly, or nothing. With no significant runner-up there is no margin to be small,
    so they are assigned normally and never enter the deferral path.
  * A LARGE margin never defers. A read may look decisively like reference A only because the
    competing genome is incomplete at that locus — deferring on a large margin would recycle
    exactly the foreign reads we do not want back. Only a small margin, meaning two references
    genuinely explain the read about equally well, is treated as ambiguity.

So deferral is strictly conservative: it can withhold a read from a bin, never add one.
"""

import argparse


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
ap.add_argument("--paf", required=True)
ap.add_argument("--query", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--min-len", type=int, default=33)
ap.add_argument("--min-pid", type=float, default=88.0)
ap.add_argument("--incl-chim", action="store_true")
ap.add_argument(
    "--defer-ambiguous",
    action="store_true",
    help="Defer reads whose best and runner-up TARGETS are within the margin to "
    ".nomatch, so a later round re-decides them against a refined reference. "
    "Off by default: assignment behaviour is unchanged unless asked for.",
)
ap.add_argument(
    "--margin-matches",
    type=int,
    default=2,
    help="Ambiguity margin in matching bases: defer when the winner explains this "
    "many or fewer extra bases than the best hit on a DIFFERENT target. "
    "Absolute floor so a short alignment cannot look decisive on a fraction.",
)
ap.add_argument(
    "--margin-frac",
    type=float,
    default=0.02,
    help="Ambiguity margin as a fraction of the winner's matches. A read defers only "
    "when it is within BOTH this and --margin-matches, so the two agree before "
    "anything is withheld.",
)
ap.add_argument(
    "--defer-report",
    help="Optional TSV of every deferred read and the two targets that tied for it. "
    "Deferral silently moves depth, so it should be auditable per round.",
)
ap.add_argument(
    "--bin-map",
    required=True,
    help="Reference FASTA whose deflines carry `bin=<label>`; maps hit targets "
    "(record ids) to bins. Required: a bin can no longer be parsed out of an id.",
)
a = ap.parse_args()
BINMAP = load_bin_map(a.bin_map)

# rid -> {target: best matches on that target}. Keyed per TARGET, not per hit: several hits to
# the same reference at different loci are one candidate, not competing evidence.
hits = {}
chim = set()  # rids rammap tagged ch:i:1

with open(a.paf) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        rid, tgt = f[0], f[5]
        matches, blk = int(f[9]), int(f[10])
        if "ch:i:1" in f[11:]:
            chim.add(rid)
        if blk < a.min_len or (100.0 * matches / blk) < a.min_pid:
            continue
        # The bin is the unit of assignment, so targets are resolved to bins for ranking — but
        # WHICH record in the bin the read actually matched is kept. Downstream has to pick one
        # record as the bin's alignment reference, and until now that information was discarded
        # here and the choice fell to "first in the file", i.e. the longest record, which is
        # unrelated to how close it is to the sample.
        rec = tgt
        # BIN for this hit target, looked up from the reference FASTA's `bin=` deflines. The
        # target is a record id (an accession in round 1, a bin-named refined ref later); deriving
        # the bin by string surgery is the mistake that made nine bugs silent.
        tgt = bin_of(tgt, BINMAP)
        per_tgt = hits.setdefault(rid, {})
        if matches > per_tgt.get(tgt, (-1, None))[0]:
            per_tgt[tgt] = (matches, rec)


def rank(per_tgt):
    """(winner, runner_up) as (target, matches) pairs; runner_up is None when only one target hit.

    Sorted by matches descending then target name ascending, so an exact tie resolves to a stable
    target instead of whichever record the aligner happened to emit first. PAF order varies with
    thread count, so the previous first-seen-wins tie-break was not reproducible.
    """
    ordered = sorted(per_tgt.items(), key=lambda kv: (-kv[1][0], kv[0]))
    flat = [(t, m) for t, (m, _rec) in ordered]
    return flat[0], (flat[1] if len(flat) > 1 else None)


best = {}  # rid -> target it is assigned to
deferred = {}  # rid -> (win_tgt, win_m, run_tgt, run_m) held back as ambiguous
# bin -> record -> reads. Which representative the bin's reads actually matched, so the
# alignment reference can be chosen on evidence instead of on file order.
subtypes = {}

for rid, per_tgt in hits.items():
    (win_tgt, win_m), runner = rank(per_tgt)
    if a.defer_ambiguous and runner is not None:
        run_tgt, run_m = runner
        margin = win_m - run_m
        # Both gates must agree. The fraction alone would defer long alignments that differ by a
        # decisive number of bases; the absolute alone would defer short ones that differ by a
        # decisive fraction.
        if margin <= a.margin_matches and (win_m <= 0 or margin / win_m <= a.margin_frac):
            deferred[rid] = (win_tgt, win_m, run_tgt, run_m)
            continue
    best[rid] = win_tgt
    rec = per_tgt[win_tgt][1]
    subtypes.setdefault(win_tgt, {})
    subtypes[win_tgt][rec] = subtypes[win_tgt].get(rec, 0) + 1

# Always written: it is small, and a bin whose reads favour a record other than the first is
# exactly the case that silently costs coverage.
with open(a.out + ".subtypes.tsv", "w") as sfh:
    sfh.write("bin\trecord\treads\n")
    for tgt in sorted(subtypes):
        for rec, n in sorted(subtypes[tgt].items(), key=lambda kv: (-kv[1], kv[0])):
            sfh.write(f"{tgt}\t{rec}\t{n}\n")

if a.defer_report:
    with open(a.defer_report, "w") as rfh:
        rfh.write("read\twin_target\twin_matches\trunnerup_target\trunnerup_matches\tmargin\n")
        for rid, (wt, wm, rt, rm) in sorted(deferred.items()):
            rfh.write(f"{rid}\t{wt}\t{wm}\t{rt}\t{rm}\t{wm - rm}\n")


def records(path):
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


with open(a.out + ".match", "w") as mfh, open(a.out + ".class", "w") as cfh, open(a.out + ".nomatch", "w") as nfh, open(
    a.out + ".chim", "w"
) as chfh:
    for hdr, seq in records(a.query):
        rid = hdr.split()[0]
        if rid in chim and not a.incl_chim:
            chfh.write(">" + hdr + "\n" + seq + "\n")
        elif rid in best:
            mfh.write(">" + hdr + "\n" + seq + "\n")
            cfh.write(hdr + "\t" + best[rid] + "\n")
        else:
            # Never-matched and deferred-as-ambiguous both land here. That is what makes a
            # deferred read recyclable: the loop advances on .nomatch.
            nfh.write(">" + hdr + "\n" + seq + "\n")

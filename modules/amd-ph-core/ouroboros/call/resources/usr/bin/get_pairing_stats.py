#!/usr/bin/env python3
"""get_pairing_stats.py — aggregate Illumina pairing-statistics TSVs into per-reference rates.

Aggregates the Illumina pairing-statistics TSV (rn<TAB>key<TAB>value) across one or more
files and emits the derived per-reference rates.

Note: reference names are sorted for deterministic output (the downstream consumer treats this
as an unordered table).
"""

import sys


def fmt(x):
    """Match perl's default numeric stringification (%.15g)."""
    if isinstance(x, int):
        return str(x)
    s = "%.15g" % x
    return s


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage:\n\tpython %s <stats1> <...>\n" % sys.argv[0])

    table = {}  # rn -> {key: summed value}
    for path in sys.argv[1:]:
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                parts = line.split("\t")
                rn = parts[0] if len(parts) > 0 else None
                key = parts[1] if len(parts) > 1 else None
                value = parts[2] if len(parts) > 2 else None
                if rn is not None and key is not None:
                    v = 0
                    if value is not None and value != "":
                        v = int(value) if value.lstrip("-").isdigit() else float(value)
                    table.setdefault(rn, {})
                    table[rn][key] = table[rn].get(key, 0) + v

    out = sys.stdout
    for rn in sorted(table.keys()):
        t = table[rn]
        obs = t.get("obs", 0)
        tmv = t.get("tmv", 0)
        fmv = t.get("fmv", 0)
        dmv = t.get("dmv", 0)
        insObs = t.get("insObs", 0)
        insErr = t.get("insErr", 0)

        if obs > 0:
            out.write("%s\tObservations\t%s\n" % (rn, fmt(obs)))
            out.write("%s\tExpectedErrorRate\t%s\n" % (rn, fmt(fmv / obs)))
            out.write("%s\tMinimumExpectedVariation\t%s\n" % (rn, fmt(tmv / obs)))
            out.write("%s\tMinimumDeletionErrorRate\t%s\n" % (rn, fmt(dmv / obs)))
        else:
            out.write("%s\tObservations\t%s\n" % (rn, fmt(obs)))
            out.write("%s\tExpectedErrorRate\t0\n" % rn)
            out.write("%s\tMinimumExpectedVariation\t0\n" % rn)
            out.write("%s\tMinimumDeletionErrorRate\t0\n" % rn)

        if insObs > 0:
            if insObs == insErr:
                insObs += 1
            insErr = insErr / insObs
            out.write("%s\tMinimumInsertionErrorRate\t%s\n" % (rn, fmt(insErr)))
        else:
            out.write("%s\tMinimumInsertionErrorRate\t0\n" % rn)


if __name__ == "__main__":
    main()

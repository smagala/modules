#!/usr/bin/env python3
"""coverage_diagram.py — Python port of coverageDiagram.R / simpleCoverageDiagram.R.

Plots per-gene coverage depth, with (full) or without (simple) a minority-variant overlay
and a minority-frequency barplot. Matplotlib only.

Usage:
  full:   coverage_diagram.py <run> <gene> <COVG.txt> <VARS.txt> <STATS.txt> <out.pdf>
  simple: coverage_diagram.py <run> <gene> <COVG.txt> <out.pdf>
"""
import sys, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ALLELE_COLORS = {"A": "#1F77B4", "C": "#FF7F0E", "G": "#2CA02C", "T": "#D62728"}

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def col(row, *names):
    for n in names:
        if n in row:
            return row[n]
    raise KeyError(names)

def main():
    a = sys.argv[1:]
    if len(a) not in (4, 6):
        sys.exit("Usage: coverage_diagram.py <run> <gene> <COVG> [<VARS> <STATS>] <out.pdf>")
    run, gene, covg = a[0], a[1], a[2]
    full = len(a) == 6
    vars_f, stat_f, out = (a[3], a[4], a[5]) if full else (None, None, a[3])

    D = read_tsv(covg)
    pos = [int(col(r, "Position")) for r in D]
    depth = [int(col(r, "Coverage Depth", "Coverage_Depth", "Coverage")) for r in D]
    cons = {int(col(r, "Position")): col(r, "Consensus") for r in D}
    Cmax = max(depth) if depth else 1

    if full:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10.5, 8))
    else:
        fig, ax1 = plt.subplots(1, 1, figsize=(10.5, 8))

    ax1.set_facecolor("black")
    ax1.vlines(pos, 0, depth, color="gray", linewidth=0.5)
    ax1.set_xlim(1, max(pos) if pos else 1)
    ax1.set_ylim(0, Cmax * 1.02)
    ax1.set_xlabel(f"{gene} position ({run})")
    ax1.set_ylabel("Coverage depth")

    if full:
        V = read_tsv(vars_f)
        depth_at = {p: d for p, d in zip(pos, depth)}
        vpos, vfreq, vcols, vlabels = [], [], [], []
        for r in V:
            p = int(r["Position"]); a_min = r["Minority_Allele"][:1]
            c = depth_at.get(p, 0); color = ALLELE_COLORS.get(a_min, "#FFFFFF")
            if c < Cmax / 2:
                ax1.vlines(p, c, Cmax, color=color, linewidth=0.8)
            else:
                ax1.vlines(p, 0, c, color=color, linewidth=0.8)
            vpos.append(p); vfreq.append(float(r["Minority_Frequency"]))
            vcols.append(color); vlabels.append(f"{cons.get(p,'?')}2{a_min}")

        # bottom: minority-frequency barplot vs expected error rate (if STATS present)
        ee = None
        try:
            for s in csv.reader(open(stat_f), delimiter="\t"):
                if len(s) >= 3 and s[1] == "ExpectedErrorRate":
                    ee = float(s[2]); break
        except OSError:
            pass
        if ee is not None:
            vals = [ee] + vfreq; cols = ["black"] + vcols; labels = ["exp. err"] + vlabels
            ax2.bar(range(len(vals)), vals, color=cols)
            ax2.axhline(ee, color="#282828", ls="--", lw=0.75)
        else:
            ax2.bar(range(len(vfreq)), vfreq, color=vcols); labels = vlabels
        ax2.set_xticks(range(len(labels)))
        ax2.set_xticklabels(labels, rotation=90, fontsize=7)
        ax2.set_ylabel("Observed frequency"); ax2.set_xlabel("minor variants")

    fig.tight_layout()
    fig.savefig(out)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""heuristic_diagram.py — Python port of heuristicDiagram.R.

Six diagnostic panels of assembled-allele quality stats, with threshold reference lines.

Usage: heuristic_diagram.py <MIN_AQ> <MIN_F> <MIN_TCC> <MIN_CONF> <ALL_ALLELES.txt> <out.pdf>
"""

import csv
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


def kde_plot(ax, data, lo, hi, title, vline=None):
    data = np.asarray([d for d in data if d is not None and np.isfinite(d)])
    sel = data[(data >= lo) & (data <= hi)] if hi > lo else data
    ax.set_title(title, fontsize=9)
    if sel.size >= 2 and np.ptp(sel) > 0:
        xs = np.linspace(lo, hi, 256)
        ax.plot(xs, gaussian_kde(sel)(xs), color="black")
    if vline is not None:
        ax.axvline(vline, color="red")


def main():
    if len(sys.argv) != 7:
        sys.exit("Usage: heuristic_diagram.py <MIN_AQ> <MIN_F> <MIN_TCC> <MIN_CONF> <ALL_ALLELES> <out.pdf>")
    MIN_AQ, MIN_F, MIN_TCC, MIN_CONF = map(float, sys.argv[1:5])
    rows = list(csv.DictReader(open(sys.argv[5]), delimiter="\t"))

    def f(name):
        out = []
        for r in rows:
            try:
                out.append(float(r[name]))
            except (ValueError, KeyError):
                out.append(np.nan)
        return np.array(out)

    aq, freq, total, conf = f("Average_Quality"), f("Frequency"), f("Total"), f("ConfidenceNotMacErr")

    fig, axes = plt.subplots(3, 2, figsize=(10.5, 8))
    aq_valid = aq[np.isfinite(aq)]
    kde_plot(axes[0, 0], aq, float(aq_valid.min()), float(aq_valid.max()), "Density of average allele quality", MIN_AQ)
    kde_plot(axes[0, 1], aq, float(aq_valid.min()), MIN_AQ, f" to {MIN_AQ:g}")
    kde_plot(axes[1, 0], freq, 0.0, 0.10, "Density of observed frequency (to 10%)", MIN_F)
    kde_plot(axes[1, 1], freq, 0.0, MIN_F, f" to {MIN_F:g}")

    tot = total[np.isfinite(total)]
    mx = np.quantile(tot, 0.20) if tot.size else 1
    axes[2, 0].hist(tot[tot <= mx], bins=50, range=(0, mx + 1), color="gray")
    axes[2, 0].set_title("Histogram of coverage (Depth <= 20% Quantile)", fontsize=9)
    axes[2, 0].axvline(MIN_TCC, color="red")

    cnz = conf[np.isfinite(conf) & (conf > 0)]
    if cnz.size:
        axes[2, 1].hist(cnz, bins=50, color="gray")
    axes[2, 1].set_title("Histogram of confidence not machine error, non-zero", fontsize=9)
    axes[2, 1].axvline(MIN_CONF, color="red")

    fig.tight_layout()
    fig.savefig(sys.argv[6])


if __name__ == "__main__":
    main()

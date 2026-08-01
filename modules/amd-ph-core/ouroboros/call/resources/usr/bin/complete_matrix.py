#!/usr/bin/env python3
"""complete_matrix.py — complete the lower diagonal of a square matrix file in place.

Completes the lower diagonal of a square matrix file in place (mirrors M[j][i] into the
upper triangle). Pure value-copy (no reformatting).

  python complete_matrix.py <lower_diagonal.txt> [-A|--annot-col]
"""
import sys

def main():
    args = [a for a in sys.argv[1:]]
    has_annot = False
    files = []
    for a in args:
        if a in ("-A", "--annot-col"):
            has_annot = True
        else:
            files.append(a)
    if len(files) != 1:
        sys.exit("Usage:\n\tpython %s <lower_diagonal.txt> [-A|--annot-col]\n" % sys.argv[0])

    path = files[0]
    with open(path) as fh:
        lines = fh.read().split("\n")
    # perl chomp(@lines) drops the trailing empty element from a final newline
    if lines and lines[-1] == "":
        lines.pop()

    N = len(lines)
    M, h, a = [], [], [None] * N
    for i in range(N):
        row = lines[i].split("\t")
        h.append(row.pop(0))
        if has_annot:
            a[i] = row.pop(0)
        M.append(row)

    with open(path, "w") as out:
        for i in range(N):
            out.write(h[i])
            if has_annot:
                out.write("\t" + a[i])
            for j in range(N):
                if j <= i:
                    out.write("\t" + M[i][j])
                else:
                    out.write("\t" + M[j][i])
            out.write("\n")

if __name__ == "__main__":
    main()

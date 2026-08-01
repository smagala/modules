#!/usr/bin/env python3
"""phase.py — variant phasing for one partition of the parallelized phasing algorithm.

Reads the called variants and read patterns that call.py writes (msgpack) and produces four
lower-triangular distance matrices over the called alleles — EXPENRD, JACCARD, MUTUALD, NJOINTP
(.sqm) — each row labelled by <position+1><allele>.

Cell (i,j) is a pairwise phasing distance between allele i and allele j: 0 = fully phased (same site
& allele), 1 = no phasing (different site with a zero co-occurrence component, or same site different
allele), otherwise a value in (0,1) from the read-pattern co-occurrence counts. With -S/-I the work is
split into an array job; the default (one shot) emits the whole triangle.

Doubles are formatted to 15 significant digits.
"""
import argparse
import math

import msgpack

NO_PHASING_VEC = (1, 1, 1, 1)
FULLY_PHASED_VEC = (0, 0, 0, 0)

ap = argparse.ArgumentParser()
ap.add_argument("-S", "--array-size", type=int, dest="array_size")
ap.add_argument("-I", "--index", type=int, dest="index")
ap.add_argument("prefix")
ap.add_argument("var_file")
ap.add_argument("pat_file")
a = ap.parse_args()

prefix = a.prefix
array_size = 10 if a.array_size is None else a.array_size
if array_size < 1:
    raise SystemExit("Array size must be 1 or greater.\n")

if a.index is not None:
    index = a.index
    if index < 1 or index > array_size:
        raise SystemExit("Index must be 1 to %d.\n" % array_size)
    index = int(index)
    prefix = prefix + "-" + "%04d" % index
else:
    index = 1
    array_size = 1


# FIXME legacy perl logic: 15-sig-fig / bare-integer formatting only exists to byte-match perl output.
def g(x):
    if isinstance(x, int):
        return str(x)
    return "%.15g" % x


with open(a.var_file, "rb") as fh:
    variants = msgpack.unpackb(fh.read(), raw=False, strict_map_key=False)
with open(a.pat_file, "rb") as fh:
    read_pats = msgpack.unpackb(fh.read(), raw=False, strict_map_key=False)

# vectorize: one entry per called allele, ordered by site index then allele
alleles = []      # [site_index, base, freq]
names = []        # "<pos+1><base>"
sites = sorted(variants)
for site_index, v in enumerate(sites):
    for b in sorted(variants[v]):
        alleles.append([site_index, b, variants[v][b]])
        names.append(str(v + 1) + b)
N = len(alleles)

O = (N ** 2 - N) // 2
if array_size > O:
    import sys
    sys.stderr.write("WARNING: array size (%d) greater than operations (%d)!\n" % (array_size, O))
    sys.stderr.write("Setting array size to %d\n" % O)
    array_size = O
    if index > array_size:
        raise SystemExit("WARNING: index (%d) greater than adjusted array size(%d). Aborting.\n"
                         % (index, array_size))

exp = open(prefix + "-EXPENRD.sqm", "w")
jac = open(prefix + "-JACCARD.sqm", "w")
mut = open(prefix + "-MUTUALD.sqm", "w")
jop = open(prefix + "-NJOINTP.sqm", "w")


def print_to_matrix(vec_mut, vec_jac, vec_exp, vec_jop, pre, suf):
    # vec_* may be a number (cell value) or '' (row-label only / blank)
    exp.write(pre + (g(vec_exp) if vec_exp != "" else "") + suf)
    jac.write(pre + (g(vec_jac) if vec_jac != "" else "") + suf)
    mut.write(pre + (g(vec_mut) if vec_mut != "" else "") + suf)
    jop.write(pre + (g(vec_jop) if vec_jop != "" else "") + suf)


def dist(i, j):
    s1, b1, Fb1 = alleles[i]
    s2, b2, Fb2 = alleles[j]
    if s1 != s2:
        if Fb1 == 0 or Fb2 == 0:
            return NO_PHASING_VEC
        total = Eb1 = Eb2 = Fb1b2 = 0
        for pat, X in read_pats.items():
            p1 = pat[s1]
            p2 = pat[s2]
            if p1 == "." or p2 == ".":
                continue
            total += X
            if p1 == b1 and p2 == b2:
                Fb1b2 += X
                Eb1 += X
                Eb2 += X
            elif p1 == b1:
                Eb1 += X
            elif p2 == b2:
                Eb2 += X
        if total != 0:
            Fb1b2 /= total
            Eb1 /= total
            Eb2 /= total
            if Fb1b2 == 0 or Eb1 == 0 or Eb2 == 0:
                return NO_PHASING_VEC
        else:
            return NO_PHASING_VEC
        mn1 = min(Eb1, Fb1)
        mn2 = min(Eb2, Fb2)
        mnA = min(mn2, mn1)
        mx1 = max(Eb1, Fb1)
        mx2 = max(Eb2, Fb2)
        mutd = 1 - Fb1b2 ** 2 / (mx1 * mx2)
        jacc = 1 - Fb1b2 / (mx1 + mx2 - Fb1b2)
        if total <= 20 or mutd == 0 or jacc == 0:
            expd = 1 - ((Fb1b2 * mnA) / (mx1 * mx2))
        else:
            expd = jacc
        njop = 1 - 2 * Fb1b2
        return (mutd, jacc, expd, njop)
    elif b1 == b2:
        return FULLY_PHASED_VEC
    else:
        return NO_PHASING_VEC


# FIXME legacy perl logic: this whole block walks the lower triangle by a linear operation count split
# across array-job shards, deriving the (row, col) start from the triangular-number inverse (sqrt) and
# chunking with Qb. pirma always runs a single shard (-S 1), so once parity ends this reduces to a
# plain nested `for row in range(N): for col in range(row+1)` over dist(). The '' blank-cell sentinel
# in print_to_matrix is part of the same row-by-row streaming and goes away with it.
Qb = O // array_size            # base quantity to do
s = 1 + Qb * (index - 1)        # starting operation
if index == array_size:
    Qb += O % array_size        # last job cleans up the remainder

if index == 1:
    print_to_matrix(*FULLY_PHASED_VEC, names[0] + "\t", "\n")

# initial triangle coordinates (triangular-number series)
row = int((1 + math.sqrt(1 + 8 * (s - 1))) / 2)
fr = (row * (row - 1)) // 2 + 1
col = s % fr

if col == 0:
    print_to_matrix("", "", "", "", names[row], "")

while col < row:
    print_to_matrix(*dist(row, col), "\t", "")
    Qb -= 1
    col += 1
    if Qb == 0:
        break
if row == col:
    print_to_matrix(*FULLY_PHASED_VEC, "\t", "\n")
row += 1

while Qb > 0:
    col = 0
    print_to_matrix("", "", "", "", names[row], "")
    while col < row:
        print_to_matrix(*dist(row, col), "\t", "")
        Qb -= 1
        col += 1
        if Qb == 0:
            break
    if row == col:
        print_to_matrix(*FULLY_PHASED_VEC, "\t", "\n")
    row += 1

exp.close()
jac.close()
mut.close()
jop.close()

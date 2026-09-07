# -*- coding: utf-8 -*-
r"""P4 -- does the closed-form main term depend on the cell shape, or not?

THE QUESTION, routed from the author.  P4's cells are indexed by how many of
3, 5, 7, 11, 13 divide N, and the paper never says why those five.  A paragraph
calling the set a choice is not enough: the claim of Measurement meas:maincoef
is that the closed form

    D'_c  =  K E_same,c[A_S] - 2 ,        K = 2 prod_{p in S} p(p-2)/(p-1)^2

reproduces the exact D'_c at every depth with nothing fitted, and
meas:secondcell shows it once more on S = (3,5,7).  Two sets is not a property
of the construction.  This sweeps it.

WHAT IS EXACT HERE AND WHAT IS NOT.  The measurement is exact and not sampled:
a cell is a set of positions in the band, so the pair average of S_2 over it is
one autocorrelation of the cell indicator -- no seed, no draw, no standard
error.  The rounding to integer pair counts is checked against n_c(n_c-1) at
every cell and the run refuses to report a cell that fails it.

E_exact,c[A_S] IS COMPUTED FROM RESIDUES ALONE, as a rational.  With
    f(q) = (q-2)/(q-1) + 1/(q-2),    g(q) = (q-1)/(q-2),
a class i is the subset of S dividing N, its residue count mod prod S is
prod_{p not in i}(p-1), and for a pair of classes (i,j)
    E[A_S | i,j] = prod_{q in i and j} g(q) * prod_{q in neither} f(q),
since q | h = N - N' is forced when q divides both, impossible when it divides
exactly one, and has probability 1/(q-1) otherwise.  Weighting by m_i m_j over
classes of the same depth reproduces the paper's rationals exactly:
403/240, 941/605, 31/15 and 16/5 at S = (3,5,7).  That identity is checked in
this file (C1) rather than assumed.

REGISTERED BEFORE THE RUN, on every cell with n_c >= NC_FLOOR:
  R1  the closed form beats saying nothing:  |main - exact| < |exact| .
  R2  |main - exact| < DIFF_TOL .
R2's band is deliberately coarse -- 250 times the largest difference at
S = (3,5,7).  The question is TRANSFER and not precision: a closed form that
does not transfer would be wrong by order one, not by a per cent.  A tight
band here would be a precision claim the sweep was not built to make.

REPORTED WHETHER OR NOT IT HELPS.  Every S attempted appears, including sets
whose deepest cell is too thin to resolve; those are printed with their counts
and excluded from R1/R2 by NC_FLOOR rather than dropped.  A sweep that shows
only the sets that worked answers nothing.

FIELD: band (2e6, 4e6]; even N; S_2 sieved to 2e6 over every prime with 2 C_2
       accumulated to 4e6; exact autocorrelation by rFFT at length 2^21.  Same
       band, same S_2 and same transform as lab_secondcell2.py, so the
       S = (3,5,7) column of this sweep is directly comparable to that file.
"""
import os
import numpy as np
from fractions import Fraction as F
from itertools import combinations

# The output path is derived, not typed, so this file runs from the packet copy
# as well as from the tree it was written in -- every other script in P4's
# corpus does the same, and a hard-coded generations/ path would make the
# packet's copy write outside itself or not at all.
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "lab_cellset_sweep.txt")

# ---- module constants (all printed; gate check G78) ------------------------
BAND      = (2_000_000, 4_000_000)
HMAX      = 2_000_000
CLIM      = 4_000_000
NC_FLOOR  = 1000          # cells thinner than this are reported, not scored
DIFF_TOL  = 0.01          # R2, deliberately coarse -- see the docstring
BASELINE  = (3, 5, 7, 11, 13)
SECOND    = (3, 5, 7)
SETS = (
    (3, 5), (3, 7), (5, 7), (3, 11), (5, 11), (7, 11), (11, 13),
    (3, 5, 7), (3, 5, 11), (3, 5, 13), (3, 7, 11), (5, 7, 11), (7, 11, 13),
    (3, 5, 7, 11), (3, 5, 7, 13), (3, 5, 11, 13), (5, 7, 11, 13),
    (3, 5, 7, 11, 13),
    (3, 5, 7, 11, 13, 17),
)

out = []


def line(s=""):
    out.append(s)


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def exact_EA(S):
    """E_exact,c[A_S] per depth, as exact rationals, from residues alone."""
    f = {q: F(q - 2, q - 1) + F(1, q - 2) for q in S}
    g = {q: F(q - 1, q - 2) for q in S}
    classes = []
    for r in range(len(S) + 1):
        for i in combinations(S, r):
            si = set(i)
            m = 1
            for p in S:
                if p not in si:
                    m *= (p - 1)
            classes.append((si, m))
    res = {}
    for c in range(len(S) + 1):
        cls = [(si, m) for si, m in classes if len(si) == c]
        tot = sum(m for _, m in cls)
        num = F(0)
        for si, mi in cls:
            for sj, mj in cls:
                e = F(1)
                for q in S:
                    if q in si and q in sj:
                        e *= g[q]
                    elif q not in si and q not in sj:
                        e *= f[q]
                num += F(mi * mj) * e
        res[c] = num / F(tot * tot)
    return res, {c: sum(m for si, m in classes if len(si) == c)
                 for c in range(len(S) + 1)}


line("lab_cellset_sweep.py -- RESULT")
line("=" * 78)
line("STATISTIC: D'_c = K E_same,c[A_S] - 2 per depth, over 19 sets S of odd")
line("           primes.  E_same,c is EXACT, by autocorrelation of the cell")
line("           indicator -- no sampling, no seed, no standard error.")
line("DENOM: n_c(n_c-1) ordered pairs, checked against the autocorrelation's")
line("       own total at every cell; a cell that fails is refused, not")
line("       reported (R3).")
line("NULL: none, and none is possible here.  The statistic is exact and")
line("      deterministic, so there is no distribution to draw from and the")
line("      tolerances R1 and R2 are absolute.  What stands in for a null is")
line("      that R1 and R2 were fixed before the run and that every set")
line("      attempted is reported, including the eight cells too thin to")
line("      score -- a sweep reporting only the sets that worked would carry")
line("      no information whatever its tolerances.")
line("module constants (G78): BAND=%s HMAX=%d CLIM=%d NC_FLOOR=%d DIFF_TOL=%s"
     % (BAND, HMAX, CLIM, NC_FLOOR, DIFF_TOL))
line("  BASELINE=%s SECOND=%s   sets attempted: %d" % (BASELINE, SECOND, len(SETS)))
line()

# ---- C1: the residue rationals reproduce the paper's, before anything else --
EA37, M37 = exact_EA(SECOND)
PAPER = {0: F(403, 240), 1: F(941, 605), 2: F(31, 15), 3: F(16, 5)}
c1 = all(EA37[c] == PAPER[c] for c in PAPER)
line("C1 -- the residue construction against the paper's rationals, S = (3,5,7)")
line("-" * 78)
for c in sorted(PAPER):
    line("  depth %d : computed %-12s   paper %-12s   %s"
         % (c, str(EA37[c]), str(PAPER[c]), "match" if EA37[c] == PAPER[c] else "DIFFER"))
line("  residue counts m_c: %s   (paper: 48, 44, 12, 1)"
     % ", ".join(str(M37[c]) for c in sorted(M37)))
line()

# ---- the band, and S_2, computed once and shared by every set --------------
pr = primes_upto(HMAX)
twin = 2.0
for p in primes_upto(CLIM):
    p = int(p)
    if p > 2:
        twin *= 1.0 - 1.0 / (p - 1.0) ** 2
S2 = np.zeros(HMAX + 1, dtype=np.float64)
S2[2::2] = twin
for p in pr:
    p = int(p)
    if p == 2:
        continue
    S2[p::p] *= (p - 1.0) / (p - 2.0)
S2[1::2] = 0.0
S2[0] = 0.0
Ns = np.arange(BAND[0] + 2, BAND[1] + 1, 2, dtype=np.int64)
L = Ns.size
nfft = 1
while nfft < 2 * L:
    nfft *= 2
kk = np.arange(1, L, dtype=np.int64)
w = S2[2 * kk]
line("  2 C_2 = %.8f;  band positions %d; transform length %d" % (twin, L, nfft))
line()

r1 = True
r2 = True
worst = 0.0
scored = 0
skipped = 0
thin = 0
rounding_fail = False
line("THE SWEEP")
line("-" * 78)
for S in SETS:
    Qp = 1
    for p in S:
        Qp *= p
    EA, M = exact_EA(S)
    K = F(2)
    for p in S:
        K *= F(p * (p - 2), (p - 1) * (p - 1))
    depth = np.zeros(L, dtype=np.int8)
    for p in S:
        depth += (Ns % p == 0)
    deepest = int((depth == len(S)).sum())
    line("S = %-22s Q' = %-7d K = %-12s deepest cell n = %d"
         % (str(S), Qp, str(K), deepest))
    line("  depth  n_c        exact D'_c      main term       difference     scored")
    for d in range(len(S) + 1):
        v = np.zeros(nfft)
        v[:L] = (depth == d)
        nc = int(v.sum())
        if nc < 2:
            line("  %-6d %-10d cell too thin to form a pair" % (d, nc))
            skipped += 1
            continue
        fv = np.fft.rfft(v)
        ac = np.rint(np.fft.irfft(fv * np.conj(fv), nfft)[:L]).astype(np.int64)
        if int(2 * ac[1:].sum()) != nc * (nc - 1):
            line("  %-6d %-10d autocorrelation did not round to exact counts "
                 "-- refusing to report" % (d, nc))
            rounding_fail = True
            skipped += 1
            continue
        den = float(nc * (nc - 1))
        ex = float(2.0 * (ac[1:] * w).sum() / den) - 2.0
        mt = float(K * EA[d]) - 2.0
        df = mt - ex
        sc = nc >= NC_FLOOR
        if not sc:
            thin += 1
        if sc:
            scored += 1
            worst = max(worst, abs(df))
            if not abs(df) < abs(ex):
                r1 = False
            if not abs(df) < DIFF_TOL:
                r2 = False
        line("  %-6d %-10d %-15.9f %-15.9f %+-14.3e %s"
             % (d, nc, ex, mt, df, "yes" if sc else "no (thin)"))
    line()

line("REGISTERED RULES")
line("-" * 78)
# Three counts and not two: `skipped` holds only the cells that could not be
# computed at all, and labelling it "not scored" made the summary contradict the
# table, which marks eight cells `no (thin)`. A count is named for what it
# holds.
line("  scored cells %d ; reported but below NC_FLOOR %d ; not computable %d ;"
     " sets attempted %d" % (scored, thin, skipped, len(SETS)))
line("  worst |main - exact| over scored cells = %.4e" % worst)
line()
for nm, held, what in (
        ("C1", c1, "the residue construction reproduces the paper's rationals at "
                   "S = (3,5,7) -- a check, not a prediction"),
        ("R1", r1, "REGISTERED: |main - exact| < |exact| at every scored cell "
                   "-- the closed form beats saying nothing"),
        ("R2", r2, "REGISTERED: |main - exact| < %s at every scored cell "
                   "(coarse on purpose; the question is transfer)" % DIFF_TOL),
        ("R3", not rounding_fail,
         "no cell was dropped for failing the pair-count rounding check")):
    line("  %s  %-8s %s" % (nm, "HOLD" if held else "FAIL", what))
line()
# The one-line summary every other run in this packet ends with, in the same
# shape, so the inventory of rules and verdicts can be counted from the runs
# rather than transcribed. Only the two REGISTERED rules go on it; C1 and R3
# are checks on the instrument and are reported above with the reason, because
# a count that mixes the two answers a different question from the one the
# manuscript asks.
# Instrument first, registered rules last: every other run in this packet puts
# its summary at the end of the file, so a reader or a parser that takes the
# last such line must land on the registered rules and not on the checks.
line("INSTRUMENT CHECKS, not registered rules:  C1 %s  R3 %s"
     % ("hold" if c1 else "REFUTED", "hold" if not rounding_fail else "REFUTED"))
line("REGISTERED RULES, SCORED:  R1 %s  R2 %s"
     % ("hold" if r1 else "REFUTED", "hold" if r2 else "REFUTED"))
line()
line("  R1 and R2 were registered before the run.  C1 and R3 are checks on the")
line("  instrument.  Cells below n_c = %d are printed with their counts and" % NC_FLOOR)
line("  excluded from scoring rather than dropped, because a set whose deepest")
line("  cell holds a handful of elements is not evidence either way.")

with open(OUT, "w", encoding="utf-8", newline="\n") as f:
    f.write("\n".join(out) + "\n")
print("\n".join(out))

# -*- coding: utf-8 -*-
r"""P4 — a SECOND cell shape, predicted before it is measured.

WHY THIS EXISTS

The paper says of its own result, twice, that "no other shape of cell
is tested".  Everything it establishes is therefore a fact about ONE
partition: the even N of a band cut by how many of 3,5,7,11,13 divide
them.  Whether that is a mechanism or a property of that partition is
not decided by any measurement in the paper.

Measurement {#meas:maincoef} changed what is possible here.  It
computes the main term of D_c from the residue structure alone, with
no input from the field and no fitted parameter, and it reproduces the
measured profile.  The same computation runs for ANY set of odd primes.
So the second cell shape can be PREDICTED FIRST and MEASURED AFTER,
which is the one thing the first partition cannot be given now.

**THIS FILE MEASURES NOTHING.**  It computes predictions and stops.
It is committed before the measurement script exists.  The measurement
is lab_secondcell.py and it does not exist yet.

THE SECOND SHAPE

  S = (3, 5, 7),  Q' = 105,  cells indexed by depth = #{p in S : p | N}.

Four cells instead of six.  Chosen for two reasons, both stated before
any number below was computed:

  * ell = B/2Q' is about 9524 at the octave (2e6, 4e6], against 66.6
    for the first partition's deepest cell.  The deficit at the
    single-class cell is 42 times smaller here, NOT 140 -- moving the
    split from 13 to 7 puts d = 11, 13, 143 and their multiples into
    the divisor sum, so ell is not the only thing that changed.  An
    earlier draft of this file said 140 and inferred from it that the
    deepest cell tests the main term nearly alone; it does not, and
    the corrected reading is in P1 below.
  * the cells are large at every depth, so the sampling error is not
    the limiting quantity anywhere -- which is the axis a registered
    rule of {#meas:maincoef} was wrongly written on.

A set disjoint from 3,5,7,11,13 would be a stronger test and is the
natural follow-up; it is not this one because at Q' = 17.19.23 the
band holds only 134 members of the deepest cell, so the deficit would
dominate again and the main-term formula would not be isolated.

WHAT IS PREDICTED

Exactly as in {#meas:maincoef}, with S in place of (3,5,7,11,13):

    D'_c = 2 C_2 * prod_{p odd, p not in S}(1 + 1/(p(p-2)))
                 * (E_same,c[A_S] - E_all[A_S])   +  deficit correction,

    A_S(h) = prod_{p | h, p in S} (p-1)/(p-2).

The constant is EXACTLY RATIONAL.  E_all[S_2] = 2 forces
K = 2 prod_{p in S} p(p-2)/(p-1)^2, since 1 + 1/(p(p-2)) is
(p-1)^2/(p(p-2)) and C_2's infinite product cancels the Euler tail term
by term.  For S = (3,5,7) that is 175/128, and depth 3's prediction is
(175/128)(16/5) - 2 = 19/8 EXACTLY.  So THE MAIN TERM carries no
irrational input, not merely no fitted parameter.  The scope of that
sentence is the main term and not the whole prediction: the deficit
correction needs E_same[T']/C_T' and C_T' is an Euler product.  It
enters multiplied by a quantity of order 2e-4, so an error of 1e-6 in
it moves the prediction by 2e-10, but a reader who checks would
otherwise find an Euler product inside a number advertised as having
none.

THE DEFICIT CORRECTION, AND ITS ONE UNPROVED PIECE.  Where the cell is
a single residue class modulo 2Q' -- depth 3 here -- the correction is
the same finite sum over k < n_c that {#meas:maincoef} uses, and it is
exact.  At the shallower depths it is DIVIDED BY m, and that model is
stated rather than left in the numbers:

  the deficit is carried by the m DIAGONAL class pairs, where
  alpha = beta so the shifts are confined exactly as in a single-class
  cell.  The m^2 - m off-diagonal pairs have their divisor classes
  spread over Z/d and are ASSUMED to cancel.  m/m^2 = 1/m.

**That assumption is precisely the cancellation {#prop:scaleinv}
declines to claim.**  It is applied here and, for consistency, also to
the first partition's row, where it moves depth 4 by 0.8 sampling
errors and nothing else by 0.1 -- weak evidence for the model,
collected before this measurement exists, which is the only time it is
worth anything.

REGISTERED BEFORE THE MEASUREMENT EXISTS

Two kinds, kept apart.  A list mixing them reads as four predictions of
which two already hold.

COMPUTED FACTS -- these are decided here and cannot be refused by any
measurement:

  F1  E_same,c[A_S] - E_all[A_S] != 0 at every depth.  Without it the
      analogue of Proposition {#prop:scaleinv} predicts no exponent at
      that cell.  (Gates this file.)
  F2  The deficit correction is largest at depth 3, the only cell that
      is a single residue class.  This is implied by the 1/m dilution
      with m = 48, 44, 12, 1, so it has no way to fail; it is printed
      because it is the shape of the correction and not because it is
      a test.

PREDICTIONS -- these can be refused by the measurement:

  P1  At depths 0, 1 and 2 the prediction matches the measured D'_c to
      within THREE SAMPLING ERRORS as the measurement reports them.
      REFUTED otherwise.
      Stated on the standardised statistic, not on a percentage: a
      registered rule of {#meas:maincoef} was refuted for exactly that
      mistake.  Restricted to depths 0-2 because THAT IS WHERE THE MAIN
      TERM IS TESTED ALONE -- undiluted, the deficit would move depth 0
      by about 0.35 sampling errors, so those three cells sit inside
      the cap whether the dilution is 1/m, 1/m^2, or absent.  Depth 3
      is excluded for the reason below, which is not that it is hard.
  P2  The second run's own SAMPLED band average of S_2 sits within two
      sampling errors of the exact 2.  REFUTED otherwise.
      The first run's was 1.998755, low by about 1.6 sampling errors.
      A second independent draw landing on the same side at the same
      size is about 1 in 15, and would say the offset is systematic
      rather than sampling -- which would reopen the truncated-prime
      mechanism that was ruled out by reading the sieve.  This converts
      "ruled out from the source" into something that can fail, and it
      costs one printed number.

DEPTH 3 IS ENUMERATED AND IS NOT SCORED.  It is one class modulo 210,
so the cell is an arithmetic progression segment and E_same,3[S_2] is
an exact 9523-term sum with no sampling error at all.  That is worth
having.  But it CANNOT be scored against the corrected prediction,
because the two are the same sum:

    corrected prediction = K.(16/5) - 2 - 2C_2.(16/5).(C_T' - E_same,3[T'])
                         = 2C_2.(16/5).E_same,3[T'] - 2
    enumerated measurement = mean_k S_2(210k) - 2
                         = 2C_2.(16/5).E_same,3[T'] - 2

identically, differing only by K's exact rational against the truncated
Euler product, which is 6.4e-08.  Any tolerance worth registering would
pass trivially.  So depth 3 is printed as a CONSISTENCY CHECK on the
implementation -- it would catch a coding error in either route -- and
is labelled as one.  It is not evidence about the main term.

This was nearly registered as a scored prediction with a 1e-5
tolerance.  A check that cannot fail prints the same thing whether the
claim is true or false, which is the first entry in this repository's
own list of ways a check stops being one.

WHAT WOULD MAKE THIS A FAILURE, SAID PLAINLY

If P1 is refuted, the main-term computation of {#meas:maincoef} is a
fact about the (3,5,7,11,13) partition and not a mechanism, and the
paper's claim would have to be narrowed to that partition.  That is
the outcome this file is built to be able to produce.

WHAT THIS DOES NOT TEST, SO THAT P1 IS NOT OVER-READ LATER

The second cell shape uses the same band, the same S_2 field, the
same pair estimator and the same code as the first.  P1 tests the
main-term formula's dependence on CELL SHAPE, which is the new thing
and the thing this was built for.  It does not test the field or the
estimator: a systematic error in either passes through both
partitions identically.  This is not an independent replication and
is not offered as one.

HOW D'_c IS FORMED, REGISTERED SO THE MEASUREMENT CANNOT DRIFT

  D'_c = E_same,c[S_2] - 2, with the 2 exact by the theorem above,
  and K = 175/128 exact -- never the run's own sampled band average
  and never the truncated float Euler product.

The first run's sampled band average was 1.998755.  Forming D'_c from
a sampled band average would put that offset into all four cells at
once, and at depth 1 -- where D'_1 = 0.126485 is 2.126485 minus 2 --
an additive offset delta in the band average shows at 7.9 delta in
D'_1 and a relative error eps in K shows at 16.8 eps.  Depth 1 is
where a shared error would surface first, and with both quantities
exact neither exists.  The run's own sampled average is computed all
the same and scored as P2.

FIELD: exact rationals for the A-differences; the constant uses
       2 C_2 = 1.32032363 and the product over odd p <= 4e6 excluding
       3, 5, 7; the deficit uses the same finite sum as lab_maincoef.py.
       No measurement is read and none is made.

BACKS: nothing yet.  This file is a registration, not evidence.
"""

import io
import os
import sys
from fractions import Fraction as F
from itertools import combinations

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "lab_secondcell_predict.txt")

CELLP2 = (3, 5, 7)                # the second shape
TWIN2 = 1.32032363169373914785    # 2 C_2
PLIM = 4_000_000
BAND = (2_000_000, 4_000_000)

lines = []


def say(s=""):
    lines.append(s)
    print(s)


def m(S, cellp):
    r = 1
    for p in cellp:
        if p not in S:
            r *= (p - 1)
    return r


def factor(p, S, T):
    a, b = p in S, p in T
    if a and b:
        return F(p - 1, p - 2)
    if a != b:
        return F(1)
    return F(1, p - 1) * F(p - 1, p - 2) + F(p - 2, p - 1)


def EA(patterns, cellp):
    M = sum(m(S, cellp) for S in patterns)
    tot = F(0)
    for S in patterns:
        for T in patterns:
            w = F(m(S, cellp) * m(T, cellp), M * M)
            a = F(1)
            for p in cellp:
                a *= factor(p, S, T)
            tot += w * a
    return tot


def sieve(n):
    sv = bytearray([1]) * (n + 1)
    sv[0] = sv[1] = 0
    i = 2
    while i * i <= n:
        if sv[i]:
            sv[i * i::i] = bytearray(len(sv[i * i::i]))
        i += 1
    return sv


def deficit(sv, CT, nc, cellp):
    """C_T - E_same[T] for a cell that is one residue class mod 2Q'.

    Same finite sum as lab_maincoef.py: shifts are h = 2Q' k with
    0 < |k| < nc, every d in the divisor expansion is coprime to 2Q',
    so no d >= nc contributes while C_T sums over every d."""
    pr = [p for p in range(3, nc) if sv[p] and p not in cellp]

    def T(k):
        t = F(1)
        mm = k
        for p in pr:
            if mm % p == 0:
                t *= F(p - 1, p - 2)
                while mm % p == 0:
                    mm //= p
        return t

    tot = F(0)
    n = 0
    for k in range(1, nc):
        w = 2 * (nc - k)
        tot += w * T(k)
        n += w
    return CT - float(tot / n)


def main():
    cellp = CELLP2
    Qp = 1
    for p in cellp:
        Qp *= p
    say("SECOND CELL SHAPE: S = %s,  Q' = %d,  2Q' = %d"
        % (", ".join(str(p) for p in cellp), Qp, 2 * Qp))
    say("band (%d, %d],  ell = B/2Q' = %.1f"
        % (BAND[0], BAND[1], (BAND[1] - BAND[0]) / (2.0 * Qp)))
    say()

    allp = [frozenset(c) for k in range(len(cellp) + 1)
            for c in combinations(cellp, k)]
    E_all = EA(allp, cellp)
    say("E_all[A_S] = %s = %.10f  (over all %d even classes mod 2Q')"
        % (E_all, float(E_all), sum(m(S, cellp) for S in allp)))
    say()

    diffs = []
    say("  depth  classes  E_same,c[A_S] exact         difference exact")
    for k in range(len(cellp) + 1):
        cell = [S for S in allp if len(S) == k]
        e = EA(cell, cellp)
        d = e - E_all
        diffs.append(d)
        say("  %-6d %-8d %-27s %s"
            % (k, sum(m(S, cellp) for S in cell), e, d))
    say()
    p1 = all(d != 0 for d in diffs)
    say("F1  the main term vanishes at no depth: %s"
        % ("hold" if p1 else "REFUTED"))
    say()

    sv = sieve(PLIM)
    CT = 1.0
    for p in range(3, PLIM + 1):
        if sv[p] and p not in cellp:
            CT *= 1.0 + 1.0 / (p * (p - 2))
    # exactly rational: see the docstring
    KEX = F(2)
    for p in cellp:
        KEX *= F(p * (p - 2), (p - 1) ** 2)
    const = float(KEX)
    say("K = 2 prod_{p in S} p(p-2)/(p-1)^2 = %s = %.10f" % (KEX, const))
    say("  EXACTLY RATIONAL -- no irrational input in the main term.")
    say("  The truncated 2 C_2 * C_T' = %.8f agrees to %.2e and is"
        % (TWIN2 * CT, abs(TWIN2 * CT - const)))
    say("  not used.")
    say("PARAMETERS: CELLP2 = %s;  2 C_2 = %s;  PLIM = %d;  band = (%d, %d]"
        % (", ".join(str(p) for p in cellp), repr(TWIN2), PLIM,
           BAND[0], BAND[1]))
    say()

    # cell sizes at this band, exactly: a class mod 2Q' meets the band
    # in floor or ceil of ell terms; the depth-k cell is the union of
    # its m(S) classes over the patterns of that size.
    ell = (BAND[1] - BAND[0]) / (2.0 * Qp)
    say("  depth  m (classes)  n_c approx   uncorrected prediction")
    ncs = []
    for k in range(len(cellp) + 1):
        cell = [S for S in allp if len(S) == k]
        mm = sum(m(S, cellp) for S in cell)
        nc = int(round(mm * ell))
        ncs.append((mm, nc))
        say("  %-6d %-12d %-12d %.6f" % (k, mm, nc, const * float(diffs[k])))
    say()

    say("DEFICIT CORRECTION.  Exact where the cell is one class")
    say("(depth %d here); divided by m elsewhere, which is an estimate"
        % len(cellp))
    say("and is marked as one.")
    say()
    say("  depth  m      n_c per class  deficit      correction   diluted")
    corr = []
    for k in range(len(cellp) + 1):
        mm, nc = ncs[k]
        per = max(2, int(round(ell)))
        dd = deficit(sv, CT, per, cellp)
        c_full = TWIN2 * float(EA([S for S in allp if len(S) == k],
                                  cellp)) * dd
        c_dil = c_full / mm
        corr.append(c_dil)
        say("  %-6d %-6d %-14d %.6e  %.6f   %.6f"
            % (k, mm, per, dd, c_full, c_dil))
    say()
    say("  depth  corrected prediction")
    for k in range(len(cellp) + 1):
        say("  %-6d %.6f" % (k, const * float(diffs[k]) - corr[k]))
    say()
    p3 = max(range(len(corr)), key=lambda i: abs(corr[i])) == len(cellp)
    say("F2  the deficit correction is largest at depth %d: %s"
        % (len(cellp), "hold" if p3 else "REFUTED"))
    say("    -- a COMPUTED FACT, not a test.  1/m with m = %s implies it."
        % ", ".join(str(ncs[k][0]) for k in range(len(cellp) + 1)))
    say()

    # NULL, run here and not declined.  Nothing is measured in this file,
    # so the null cannot be about measurement error.  What it can and must
    # be about is whether the prediction has arithmetic content at all:
    # this paper's whole thesis is that the cell floor reads the
    # correspondence and not the cell SIZES, so a predicted profile that
    # merely tracks size would carry nothing.  With four cells the null is
    # enumerable -- all 24 relabellings.
    say("NULL: does the predicted profile just track cell size?")
    say()
    from itertools import permutations
    pred = [const * float(diffs[k]) - corr[k] for k in range(len(cellp) + 1)]
    sizes = [ncs[k][1] for k in range(len(cellp) + 1)]
    say("  depth   n_c        predicted D'_c")
    for k in range(len(cellp) + 1):
        say("  %-6d  %-9d  %.6f" % (k, sizes[k], pred[k]))
    say()

    def conc(order):
        """concordant pairs between the predicted order and 1/n_c."""
        n = 0
        idx = list(range(len(pred)))
        for a in idx:
            for b in idx:
                if a < b:
                    if (pred[order[a]] - pred[order[b]]) * \
                       (1.0 / sizes[a] - 1.0 / sizes[b]) > 0:
                        n += 1
        return n

    base = conc(list(range(len(pred))))
    dist = [conc(list(p)) for p in permutations(range(len(pred)))]
    beat = sum(1 for v in dist if v >= base)
    say("  concordant pairs with 1/n_c, true labelling: %d of %d"
        % (base, len(pred) * (len(pred) - 1) // 2))
    say("  over the %d relabellings: max %d, mean %.2f, count >= %d: %d"
        % (len(dist), max(dist), sum(dist) / float(len(dist)), base, beat))
    say("  p = %d/%d = %.4f" % (beat, len(dist), beat / float(len(dist))))
    say("  A profile that merely tracked size would sit at the top of this")
    say("  distribution.  Where it actually sits is printed above, and P2")
    say("  is scored against the measurement rather than against this.")
    say()
    say("=" * 70)
    say("COMPUTED FACTS:  F1 %s   F2 %s"
        % ("hold" if p1 else "REFUTED", "hold" if p3 else "REFUTED"))
    say("PREDICTIONS:     P1 and P2 CANNOT BE SCORED HERE.  They need")
    say("the measurement, which does not exist.  That is the point of")
    say("this file.  Depth 3 is enumerated by the measurement and is a")
    say("consistency check, not a score -- see the header.")

    head = [
        "STATISTIC: the predicted D'_c for a SECOND cell shape,",
        "           S = (3, 5, 7), computed from the residue structure",
        "           alone before any measurement of that shape exists.",
        "DENOM: E_all[A_S] over all even classes mod 2Q'.",
        "NULL: an exact enumeration, run in this file. Nothing is",
        "      measured here, so the null asks the only thing it can:",
        "      whether the predicted profile merely tracks cell size,",
        "      which this paper's thesis says it must not. With four",
        "      cells all 24 relabellings are enumerated and the true",
        "      labelling is placed in that distribution. P2's own null",
        "      is the sampling error of the measurement that does not",
        "      yet exist, and P2 is stated against it.",
        "FIELD: exact rationals throughout the A-differences; the",
        "       constant uses 2 C_2 = 1.32032363 and the product over",
        "       odd p <= 4e6 excluding 3, 5, 7; band (2e6, 4e6].",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    return 0 if p1 else 1


if __name__ == "__main__":
    sys.exit(main())

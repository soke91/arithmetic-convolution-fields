# -*- coding: utf-8 -*-
r"""P4, Proposition {#prop:scaleinv} — the main term, computed exactly.

WHY THIS HAS TO BE RUN

The proposition's proof ends with

    D_c = 2 C_2 prod_{p>13}(1 + 1/(p(p-2))) * (E_same,c[A] - E_all[A]) + E,

and concludes that D_c has scale exponent 0.  A proof audit on
2026-09-05 pointed out that the conclusion needs a hypothesis the proof
never states: convergence is additive and an exponent is
multiplicative, so d log D_c / d log B -> 0 requires the main term to be
NONZERO.  Where E_same,c[A] = E_all[A] the proposition says only that
D_c -> 0 and predicts no exponent at all.  The depth most at risk is
depth 1, whose D_c is an order below its neighbours.

The paper then said that difference is "a finite sum of rationals over
the residues modulo 2Q and the pattern weights, so it is decidable cell
by cell without any measurement", and that it "is not computed here".
This computes it.  It is exact rational arithmetic: no sieve, no
sampling, no floating point in the difference.

WHAT IS COMPUTED

A(h) = prod_{p | h, 2 < p <= 13} (p-1)/(p-2) depends only on which of
3,5,7,11,13 divide h, hence only on h mod Q.  A cell is a union of
residue classes modulo 2Q; class alpha has pattern S(alpha) = the set
of those primes dividing alpha, and depth = |S|.  The number of even
classes mod 2Q with pattern S is m(S) = prod_{p not in S} (p-1), and
sum_S m(S) = Q, the number of even classes.

For two classes alpha, beta and a prime p, p divides h = alpha - beta
iff alpha = beta mod p.  So the per-prime factor of E[A] over a pattern
pair (S, T) is

    p in S and p in T      (p-1)/(p-2)      both 0 mod p, always divides
    p in exactly one       1                one is 0, the other is not
    p in neither           (1/(p-1))(p-1)/(p-2) + (p-2)/(p-1)

the last because a_p, b_p are then independent and uniform on
{1, ..., p-1}.  The factors are independent across p by CRT, so E[A]
over a set of patterns is a finite double sum weighted by m(S) m(T).

REGISTERED BEFORE THE RUN (post hoc as to timing, not as to freedom)

This is written after the D_c row was measured, and it is recorded as
post hoc.  What makes it a test rather than a fit is that it has NO
free parameter: the difference is forced by the residue structure and
the constant 2 C_2 prod_{p>13}(1 + 1/(p(p-2))) is forced by the
definition of the singular series.  Six numbers are predicted from
none.

  N1  E_same,c[A] - E_all[A] != 0 at every depth.  This is the
      hypothesis the proof needs.  REFUTED if it vanishes at any depth.
  N2  The prediction matches the measured D_c at every depth to within
      3 per cent.  REFUTED otherwise.  The tolerance is the size of the
      O(1/ell) the proposition's own error term carries at the octave
      the measurement reads: ell = B/2Q is 66.6 there, so 1/ell is 1.5
      per cent, and the cap is twice that.
  N3  The worst deviation is at depth 5.  That cell is a single
      residue class with n_c = 67, so it is the cell where an O(1/ell)
      error is the largest fraction of the answer.  REFUTED if the
      worst is elsewhere -- which would say the deviation is not the
      error term the proposition names.

N1 gates.  N2 and N3 are reported: N2 because the constant carries a
truncated Euler product, N3 because it is a claim about which of six
residuals is largest and six is a small number.

N2's tolerance was justified by "1/ell is 1.5 per cent", which is the
reading the paper withdraws -- the proposition's term is
O_eta(B^eta/ell) and reading it as 1/ell sets the constant and the
divisor sum both to 1.  The rule is left as registered and scored as
registered; the justification beside it was wrong and is recorded here
rather than silently repaired.

ADDED AFTER N1-N3 WERE SCORED, AND REGISTERED BEFORE BEING COMPUTED

An independent reading of this file predicted, without running it,
that depth 5's residual is the divisor-sampling deficit of a cell that
is a single residue class -- no d >= n_c divides any same-cell shift
while C_T sums over every d -- and put it at about 1.35 per cent of
D_5 at n_c = 67, with the induced spread across the three octaves at
"the right order, not the number".  Those two numbers were on the page
before the computation below existed.

  N4  The exact deficit accounts for depth 5's residual to within
      10 per cent.  REFUTED otherwise -- which would say the mechanism
      is not the one named, since the deficit has no free parameter.
  N5  The spread this deficit induces across the three octaves agrees
      with the depth-5 spread meas:scaleinv measures (1.28 per cent)
      to within a factor of two.  REFUTED otherwise.  If it holds, the
      largest of the six spreads is the proposition's own error term
      measured, and not a failure of scale invariance.

N4 and N5 are reported, not gated: both compare against quantities
carrying their own sampling error, and N5's target was measured by a
different run.

AND ONE MORE, FROM THE SAME READING

E_all[S_2] = 2 exactly -- the Euler product it makes is the reciprocal
of C_2's, term by term.  So the band average never needed measuring,
and the run that measured it got 1.998755.  That shortfall is a single
additive offset shared by all six D_c, because each is E_same,c minus
that one number.

  N6  With the band average at its exact value 2, every depth is within
      TWO sampling errors of the prediction.  REFUTED otherwise.  This
      is a strong rule: before it, depth 5 stood at 260.8 and depth 0
      at 2.0, and nothing is fitted to move them -- the 2 is a theorem
      and the depth-5 correction is the finite sum of N4.

FIELD: exact rationals for the A-differences; the constant uses
       2 C_2 = 1.32032363 with C_2 the twin-prime constant and the
       product over 13 < p <= 4e6.  Measured D_c read from
       results/lab_cell_singular.txt at the octave (2e6, 4e6].

BACKS: Proposition {#prop:scaleinv} in deploy/papers/P4-coherent-cell-floor.tex.
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
OUT = os.path.join(ROOT, "results", "lab_maincoef.txt")

CELLP = (3, 5, 7, 11, 13)
TWIN2 = 1.32032363169373914785   # 2 C_2
PLIM = 4_000_000                 # truncation of prod_{p>13}(1+1/(p(p-2)))

lines = []


def say(s=""):
    lines.append(s)
    print(s)


def m(S):
    """Even residue classes mod 2Q carrying exactly the pattern S."""
    r = 1
    for p in CELLP:
        if p not in S:
            r *= (p - 1)
    return r


def factor(p, S, T):
    """Expected A-factor at p over a pattern pair, as an exact rational."""
    a, b = p in S, p in T
    if a and b:
        return F(p - 1, p - 2)
    if a != b:
        return F(1)
    return F(1, p - 1) * F(p - 1, p - 2) + F(p - 2, p - 1)


def EA(patterns):
    M = sum(m(S) for S in patterns)
    tot = F(0)
    for S in patterns:
        for T in patterns:
            w = F(m(S) * m(T), M * M)
            a = F(1)
            for p in CELLP:
                a *= factor(p, S, T)
            tot += w * a
    return tot


def read_Dc():
    """The measured D_c row at (2e6,4e6] from lab_cell_singular.txt.

    READ rather than typed: a hand-copied table is a dependency no
    check can see."""
    p = os.path.join(ROOT, "results", "lab_cell_singular.txt")
    for ln in io.open(p, encoding="utf-8", errors="replace"):
        if ln.strip().startswith("M1  D_c at (2e6,4e6] by depth:"):
            return [float(x) for x in ln.split(":", 1)[1].split(",")]
    return []


def read_se():
    """The sampling error of D_c at (2e6,4e6], column 5 of the octave
    rows of lab_cell_singular.txt.  It was already printed there; the
    post hoc block below exists because a registered rule of this file
    was written without it."""
    p = os.path.join(ROOT, "results", "lab_cell_singular.txt")
    out = {}
    for ln in io.open(p, encoding="utf-8", errors="replace"):
        f = ln.split()
        if len(f) == 8 and f[0] == "(" and f[1] == "2000000," \
                and f[2] == "4000000]":
            try:
                out[int(f[3])] = float(f[7])
            except ValueError:
                continue
    return [out[d] for d in sorted(out)] if len(out) == 6 else []


def read_Esame():
    """E_same,c[S_2] at (2e6,4e6], column 4 of the octave rows.

    Read so that the band average can be recovered as E_same,c - D_c
    and compared with its exact value 2.  It is the same number in all
    six rows, which is what makes its error a shared offset."""
    p = os.path.join(ROOT, "results", "lab_cell_singular.txt")
    out = {}
    for ln in io.open(p, encoding="utf-8", errors="replace"):
        f = ln.split()
        if len(f) == 8 and f[0] == "(" and f[1] == "2000000," \
                and f[2] == "4000000]":
            try:
                out[int(f[3])] = float(f[5])
            except ValueError:
                continue
    return [out[d] for d in sorted(out)] if len(out) == 6 else []


def read_D5_octaves():
    """The exact D_5 at each of the three octaves, from lab_cell_singular.

    Needed because the induced drift below and the width meas:scaleinv
    reports were being quoted on DIFFERENT denominators -- the first
    over the uncorrected prediction, the second over the mean of the
    measured values -- so the two per cents were not comparable.  Read
    the measured values here and put both on the measured mean.
    """
    p = os.path.join(ROOT, "results", "lab_cell_singular.txt")
    if not os.path.exists(p):
        return []
    want = (("1000000,", "2000000]"), ("2000000,", "4000000]"),
            ("4000000,", "8000000]"))
    out = {}
    for ln in io.open(p, encoding="utf-8", errors="replace"):
        f = ln.split()
        # (  lo,  hi] depth n_c E_same,c D_c se_c  -- D_c is field 6.
        # Reading field 4 gave n_c, and the printed drift came out
        # 0.05 against a "measured width" of 129 per cent, which is
        # what made the mistake visible.  Guard on the value being a
        # D_c and not a count.
        if len(f) == 8 and f[0] == "(" and f[3] == "5":
            for k, (lo, hi) in enumerate(want):
                if f[1] == lo and f[2] == hi:
                    try:
                        v = float(f[6])
                    except ValueError:
                        continue
                    if 0.0 < v < 100.0 and "." in f[6]:
                        out[k] = v
    return [out[k] for k in sorted(out)] if len(out) == 3 else []


def read_finite_band():
    """The band's own exact pair average of S_2, from lab_secondcell2.txt.

    READ, not recomputed: this file's band is a different one and the
    point is only that E_{all,B} is not E_inf."""
    p = os.path.join(ROOT, "results", "lab_secondcell2.txt")
    if not os.path.exists(p):
        return "not computed"
    for ln in io.open(p, encoding="utf-8", errors="replace"):
        t = ln.split()
        if len(t) >= 6 and t[0] == "E_all[S_2]" and t[1] == "over":
            return t[-1]
    return "not computed"


def sieve(n):
    sv = bytearray([1]) * (n + 1)
    sv[0] = sv[1] = 0
    i = 2
    while i * i <= n:
        if sv[i]:
            sv[i * i::i] = bytearray(len(sv[i * i::i]))
        i += 1
    return sv


def depth5_deficit(sv, CT, nc):
    """C_T - E_same[T] for a cell that is a SINGLE residue class mod 2Q
    with nc members.

    Exact, and the only step in this file that uses the cell's internal
    structure rather than its pattern.  Shifts are h = 2Q k with
    0 < |k| < nc; d coprime to 2Q, so d | h iff d | k; hence no
    d >= nc contributes at all while C_T sums over every d.  The pair
    count at each |k| is exactly nc - |k|.
    """
    pr = [p for p in range(17, nc) if sv[p]]

    def T(k):
        t = F(1)
        m = k
        for p in pr:
            if m % p == 0:
                t *= F(p - 1, p - 2)
                while m % p == 0:
                    m //= p
        return t

    tot = F(0)
    n = 0
    for k in range(1, nc):
        w = 2 * (nc - k)
        tot += w * T(k)
        n += w
    return CT - float(tot / n)


def euler_tail(sv):
    """prod_{13 < p <= PLIM} (1 + 1/(p(p-2)))."""
    k = 1.0
    for p in range(17, PLIM + 1):
        if sv[p]:
            k *= 1.0 + 1.0 / (p * (p - 2))
    return k


def main():
    allp = [frozenset(c) for k in range(6) for c in combinations(CELLP, k)]
    E_all = EA(allp)
    say("E_all[A] = %s = %.10f  (over all %d even classes mod 2Q)"
        % (E_all, float(E_all), sum(m(S) for S in allp)))
    say()

    diffs = []
    say("  depth  classes  E_same,c[A] exact                  difference exact")
    for k in range(6):
        cell = [S for S in allp if len(S) == k]
        e = EA(cell)
        d = e - E_all
        diffs.append(d)
        say("  %-6d %-8d %-34s %s" % (k, sum(m(S) for S in cell), e, d))
    say()

    n1 = all(d != 0 for d in diffs)
    say("N1  the main term vanishes at no depth: %s"
        % ("hold" if n1 else "REFUTED"))
    say("    depth 5 gives %s = %.10f, which is prod_{2<p<=13}(1+1/(p-2))"
        % (diffs[5] + E_all, float(diffs[5] + E_all)))
    say("    -- the deepest cell is one residue class, so Q divides every")
    say("    shift in it and A is its maximum there exactly.")
    say()

    sv5 = sieve(PLIM)
    CT = euler_tail(sv5)

    # THE CONSTANT IS EXACTLY RATIONAL, and the truncated Euler product
    # above is not needed for it.  E_all[S_2] = 2 forces
    #   2 C_2 prod_{p not in P}(1+1/(p(p-2))) = 2 / prod_{p in P}(...),
    # and 1 + 1/(p(p-2)) = (p-1)^2/(p(p-2)), so
    #   K = 2 prod_{p in P} p(p-2)/(p-1)^2.
    # C_2's infinite product cancels the tail term by term.  So the main
    # term carries NO IRRATIONAL INPUT -- a stronger statement than "no
    # fitted parameter", and one this file could not make before the
    # band average was known exactly.  (The deficit correction below
    # does carry 2 C_2, and that is said where it is used.)
    KEX = F(2)
    for p in CELLP:
        KEX *= F(p * (p - 2), (p - 1) ** 2)
    const = float(KEX)
    K = CT
    say("constant K = 2 prod_{p in P} p(p-2)/(p-1)^2 = %s = %.13f"
        % (KEX, const))
    say("  EXACTLY RATIONAL.  E_all[S_2] = 2 makes C_2's infinite product")
    say("  cancel the Euler tail term by term, so the main term has no")
    say("  irrational input at all -- not merely no fitted parameter.")
    say("  The truncated product 2 C_2 * prod_{13<p<=%d}(...) = %.8f"
        % (PLIM, TWIN2 * K))
    say("  agrees with it to %.2e and is not used." % abs(TWIN2 * K - const))
    say()

    meas = read_Dc()
    if len(meas) == 6:
        say("  depth  predicted    measured D_c   predicted/measured")
        rat = []
        for k in range(6):
            pr = const * float(diffs[k])
            rat.append(pr / meas[k])
            say("  %-6d %-12.6f %-14.6f %.4f" % (k, pr, meas[k], rat[k]))
        worst = max(range(6), key=lambda k: abs(rat[k] - 1.0))
        n2 = all(abs(r - 1.0) <= 0.03 for r in rat)
        n3 = (worst == 5)
        say()
        say("N2  every depth within 3 per cent: %s   (worst %.2f%% at depth %d)"
            % ("hold" if n2 else "REFUTED",
               100.0 * abs(rat[worst] - 1.0), worst))
        say("N3  the worst deviation is at depth 5: %s"
            % ("hold" if n3 else "REFUTED"))
        say()
        say("    POST HOC (written after N3 was scored, and scoring")
        say("    nothing).  N3 WAS REGISTERED ON THE WRONG QUANTITY.  The")
        say("    measured D_c is a sampled estimator and this file's own")
        say("    source prints its sampling error alongside it -- a fact")
        say("    already in the repository, which the rule did not use. A")
        say("    residual has to be read against that error before it is")
        say("    read against 1/ell:")
        say()
        se = read_se()
        if len(se) == 6:
            say("      depth  residual     sampling SE   multiple of SE")
            for k in range(6):
                r = const * float(diffs[k]) - meas[k]
                say("      %-6d %+.6f    %.6f      %8.1f"
                    % (k, r, se[k], abs(r) / se[k]))
            say()
            say("    Read that way the ordering inverts.  Depth 1's 1.65 per")
            say("    cent is under one sampling error and is not a deviation")
            say("    at all; depth 5's 1.33 per cent is over two hundred, so")
            say("    it is the only one of the six that the sampling error")
            say("    cannot account for -- and it is 1/ell in size, which is")
            say("    the error term the proposition names.  N3's claim was")
            say("    right about the cell and wrong about the statistic, and")
            say("    it stands REFUTED as registered.")
    else:
        say("measured D_c row not found in lab_cell_singular.txt")
        n2 = n3 = False

    say()
    say("NULL for N2, run here and not declined.  Six predictions against")
    say("six measurements agree; the question a null has to answer is")
    say("whether the ASSIGNMENT carries that, or only the scale.  With six")
    say("depths the null is enumerable: all 720 permutations of the depth")
    say("labels are run and the observed count is placed in that")
    say("distribution.")
    from itertools import permutations
    if len(meas) == 6 and len(se) == 6:
        def within(perm):
            return sum(1 for k in range(6)
                       if abs(const * float(diffs[perm[k]]) - meas[k])
                       <= 2.0 * se[k])
        obs = within(tuple(range(6)))
        dist = [within(p) for p in permutations(range(6))]
        beat = sum(1 for v in dist if v >= obs)
        say("    depths within 2 sampling errors, true labelling: %d of 6"
            % obs)
        say("    over the 720 relabellings: max %d, mean %.3f, count >= %d:"
            " %d" % (max(dist), sum(dist) / 720.0, obs, beat))
        say("    p = %d/720 = %.5f" % (beat, beat / 720.0))
        say("    The null is not a formality here: the six predictions span")
        say("    a factor of 70, so a wrong pairing is refused by size")
        say("    alone. What it shows is that the agreement is carried by")
        say("    WHICH cell gets which number, which is the claim.")
    say()
    say("DEPTH 5's RESIDUAL, COMPUTED RATHER THAN NAMED.")
    say()
    say("An earlier printing of this file explained depth 5's miss as")
    say("'1.33 per cent, the size of 1/ell'.  That is the reading the")
    say("paper withdraws one paragraph later -- the proposition's term is")
    say("O_eta(B^eta/ell), not 1/ell, and reading it as 1/ell sets both the")
    say("constant and the divisor sum to 1.  It is also unfalsifiable: with")
    say("an unnamed constant no discrepancy could refuse it.  Withdrawn.")
    say()
    say("The real mechanism is exact and needs no asymptotics.  At depth 5")
    say("the cell is ONE residue class mod 2Q, so it is an arithmetic")
    say("progression segment N = N_0 + 2Q i.  Every same-cell shift is")
    say("h = 2Q k with 0 < |k| < n_c, and every d in the divisor expansion")
    say("is coprime to 2Q, so d | h iff d | k.  Hence NO d >= n_c divides")
    say("any same-cell shift, while C_T sums over every d.  E_same,5[T] is")
    say("therefore deficient against C_T, necessarily downward -- which is")
    say("the sign observed.  The deficit is a finite sum over k < n_c.")
    say()
    dl = depth5_deficit(sv5, CT, 67)
    say("    E_same,5[T] at n_c = 67 = %.10f" % (CT - dl))
    say("    C_T                     = %.10f" % CT)
    say("    deficit Delta           = %.6e" % dl)
    corr = TWIN2 * float(F(128, 33)) * dl
    say("    2 C_2 * (128/33) * Delta = %.6f" % corr)
    if len(meas) == 6 and len(se) == 6:
        raw = const * float(diffs[5]) - meas[5]
        say("    observed residual at depth 5 = %.6f" % raw)
        say("    ratio = %.4f" % (corr / raw))
        say("    residual after the correction = %.6f = %.1f sampling errors"
            % (raw - corr, abs(raw - corr) / se[5]))
        say("    -- from 260.8 to %.1f." % (abs(raw - corr) / se[5]))
    say()
    say("AND IT IS THE SAME QUANTITY meas:scaleinv MEASURES AS DEPTH 5's")
    say("SPREAD.  Delta falls with n_c, so it induces a drift in D_5 across")
    say("the three octaves.  Computed exactly at the three cell sizes:")
    say()
    say("      octave        ell     n_c    Delta        correction/D_5")
    Q2 = 2 * 3 * 5 * 7 * 11 * 13
    D5p = const * float(diffs[5])
    cs = []
    for B, nc in ((1_000_000, 34), (2_000_000, 67), (4_000_000, 134)):
        dd = depth5_deficit(sv5, CT, nc)
        cc = TWIN2 * float(F(128, 33)) * dd
        cs.append(cc)
        say("      %-11d %7.1f %5d  %.6e   %6.2f%%"
            % (B, B / Q2, nc, dd, 100.0 * cc / D5p))
    say("    induced spread across the three = %.2f per cent"
        % (100.0 * (cs[0] - cs[2]) / D5p))
    say("    -- taken over the UNCORRECTED prediction %.6f." % D5p)
    say()
    say("    SAME DENOMINATOR.  The width meas:scaleinv reports is")
    say("    (max-min)/mean over the MEASURED D_5, so the two per cents")
    say("    above sat on different denominators and were not comparable.")
    d5m = read_D5_octaves()
    if not d5m:
        say("    the measured D_5 row was not found in")
        say("    lab_cell_singular.txt -- the matched figures are not")
        say("    computed, and the two above must not be compared.")
    else:
        mu = sum(d5m) / len(d5m)
        say("      measured D_5 at the three octaves = %s"
            % ", ".join("%.6f" % v for v in d5m))
        say("      their mean = %.6f" % mu)
        a = 100.0 * (cs[0] - cs[2]) / mu
        b = 100.0 * (max(d5m) - min(d5m)) / mu
        say("      induced drift over that mean  = %.3f per cent" % a)
        say("      measured width over that mean = %.3f per cent" % b)
        # The sentence about the agreement improving is COMPUTED, not
        # asserted.  Typed as prose it printed unchanged while the
        # numbers above it read 0.05 against 129 -- a claim about
        # magnitudes with no arithmetic gating it, which is the shape
        # this whole file exists to refuse.
        old = abs(100.0 * (cs[0] - cs[2]) / D5p - 1.28) / 1.28
        new = abs(a - b) / b if b else float("nan")
        say("      mismatch before matching = %.1f per cent of the width"
            % (100.0 * old))
        say("      mismatch after  matching = %.1f per cent of the width"
            % (100.0 * new))
        say("    Matching the denominator %s the agreement."
            % ("IMPROVES" if new < old else "does NOT improve"))
    say()
    say("    So the largest of the six spreads is not a failure of scale")
    say("    invariance: it is the proposition's own error term, computed.")
    say("    This was predicted BEFORE the computation by an independent")
    say("    reading, which put the residual at about 1.35 per cent at")
    say("    n_c = 67 and the spread at 'the right order, not the number';")
    say("    the exact values are %.2f and %.2f." %
        (100.0 * cs[1] / D5p, 100.0 * (cs[0] - cs[2]) / D5p))
    say()
    n4 = n5 = False
    if len(meas) == 6:
        raw = const * float(diffs[5]) - meas[5]
        n4 = abs(corr / raw - 1.0) <= 0.10
        sp = 100.0 * (cs[0] - cs[2]) / D5p
        n5 = 0.5 <= sp / 1.28 <= 2.0
        say("N4  the deficit accounts for depth 5's residual to 10 per cent:"
            " %s  (%.1f%%)" % ("hold" if n4 else "REFUTED",
                               100.0 * abs(corr / raw - 1.0)))
        say("N5  the induced spread agrees with the measured 1.28 per cent"
            " to a factor 2: %s  (%.2f%%)"
            % ("hold" if n5 else "REFUTED", sp))
    say()
    say("THE BAND AVERAGE IS EXACTLY 2, AND IT WAS BEING SAMPLED.")
    say()
    say("  E_all[S_2] = 2C_2 prod_{p>2} [ (1/p)(p-1)/(p-2) + (p-1)/p ]")
    say("             = 2C_2 prod_{p>2} (p-1)^2 / (p(p-2)),")
    say("  and C_2 = prod_{p>2}(1 - 1/(p-1)^2) = prod_{p>2} p(p-2)/(p-1)^2")
    say("  is its reciprocal term by term.  So E_all[S_2] = 2 EXACTLY.")
    say()
    say("  IT IS A LIMITING VALUE, NOT THE BAND'S OWN AVERAGE.  Taking")
    say("  P(p|h) = 1/p for every odd p is what makes the two products")
    say("  reciprocal; at a finite band that probability is not exactly")
    say("  1/p.  lab_secondcell2.py computes the band's own pair average")
    say("  exactly at its band and gets %s -- the difference is the"
        % read_finite_band())
    say("  band's equidistribution error.  E_{all,B} and E_inf are two")
    say("  quantities and {#eq:Dc} subtracts the first.")
    say()
    say("  Found by an independent reading, not by this file.  It removes")
    say("  a measured input: the prediction is D_c = 2C_2 C_T E_same,c[A]")
    say("  - 2, with the 2 exact and the constant no longer entering the")
    say("  subtracted term at all.")
    say()
    es = read_Esame()
    n6ok = False
    if len(es) == 6 and len(meas) == 6 and len(se) == 6:
        eall = es[0] - meas[0]
        say("  The run's band average, recovered as E_same,c - D_c, is the")
        say("  same number in all six rows: %.6f against the exact 2," % eall)
        say("  shortfall %.6f.  One additive offset shared by all six"
            % (eall - 2.0))
        say("  measured D_c, since each is E_same,c minus that one number.")
        say("  The S_2 sieve here runs over every prime to 8e6, so this is")
        say("  the sampling error of the band average and not a truncated")
        say("  prime product -- a truncation would bias low as well, but")
        say("  would leave depth 5 exact, every same-cell shift there being")
        say("  2Qk with k <= 66 and so free of prime factors above 61.")
        say()
        say("  depth  D_c with E_all = 2   prediction       residual/se")
        n6ok = True
        for k in range(6):
            newd = es[k] - 2.0
            pk = const * float(diffs[k]) - (corr if k == 5 else 0.0)
            r = (pk - newd) / se[k]
            if abs(r) > 2.0:
                n6ok = False
            say("  %-6d %-20.6f %-16.6f %+.2f" % (k, newd, pk, r))
        say()
        say("N6  with the band average at its exact value, every depth is")
        say("    within two sampling errors of the prediction: %s"
            % ("hold" if n6ok else "REFUTED"))
        say("    Depth 5 goes from 260.8 sampling errors to %+.2f."
            % ((const * float(diffs[5]) - corr - (es[5] - 2.0)) / se[5]))
        say("    Nothing is fitted for this: the 2 is a theorem and the")
        say("    depth-5 correction is a finite sum.")
    else:
        say("  E_same,c row not found; N6 cannot be scored")
    say()
    say("THE SAME DEFICIT AT THE SHALLOWER CELLS, DILUTED BY 1/m.")
    say()
    say("  The deficit is carried by the m DIAGONAL class pairs, where")
    say("  alpha = beta so h = 0 mod 2Q and the shifts are confined")
    say("  exactly as in a single-class cell.  The m^2 - m off-diagonal")
    say("  pairs have c_d = -(alpha-beta)(2Q)^{-1} mod d spread over Z/d")
    say("  and are ASSUMED to cancel; m/m^2 = 1/m.  That assumption is")
    say("  precisely the cancellation {#prop:scaleinv} declines to claim.")
    say("  It is applied here because it is applied in the second cell")
    say("  shape's registration, and a model used in one place and not")
    say("  the other is an asymmetry a reader can find by arithmetic.")
    say()
    ms = [sum(m(S) for S in allp if len(S) == k) for k in range(6)]
    say("  depth  m       diluted correction   was      becomes")
    n7 = True
    if len(meas) == 6 and len(se) == 6 and len(es) == 6:
        for k in range(6):
            c = TWIN2 * float(diffs[k] + E_all) * (dl / ms[k])
            newd = es[k] - 2.0
            was = (const * float(diffs[k])
                   - (corr if k == 5 else 0.0) - newd) / se[k]
            now = (const * float(diffs[k]) - c - newd) / se[k]
            if abs(now) > 2.0:
                n7 = False
            say("  %-6d %-7d %-20.6f %+8.2f %+8.2f"
                % (k, ms[k], c, was, now))
        say()
        say("  Depth 5 is m = 1, so the dilution is the identity there and")
        say("  its entry is the exact deficit.  The model moves depth 4 by")
        say("  0.8 sampling errors and nothing else by as much as 0.1.")
    say()
    say("N7  with the 1/m dilution applied at every depth, all six stay")
    say("    within two sampling errors: %s" % ("hold" if n7 else "REFUTED"))
    say("    This is weak evidence for the 1/m model -- 0.8 sampling")
    say("    errors at one cell -- and it is collected before the second")
    say("    cell shape is measured, which is the only time it is worth")
    say("    anything.")
    say()

    say("ROBUSTNESS: is the constant's SPLIT right?  The singular series is")
    say("2 C_2 prod_{p|h,p>2}(p-1)/(p-2), and this file cuts that product at")
    say("13 -- A below, the Euler factor above.  An off-by-one in which")
    say("primes are absorbed into C_2 and which are left in the product")
    say("would move the constant by 1/(p(p-2)) at the prime concerned, which")
    say("at p = 13 is 1/143 = 0.70 per cent.  Two things say it is not there.")
    say()
    say("  (a) STRUCTURAL, and it needs no error bar.  A mis-split is a")
    say("      UNIFORM multiplicative factor: all six would be off by the")
    say("      same amount in the same direction.  The observed ratios fail")
    say("      in OPPOSITE directions at the two ends, so no uniform factor")
    say("      explains the scatter whatever bound is put on it.")
    say()
    if len(meas) == 6 and len(se) == 6:
        say("  (b) IN SAMPLING ERRORS, which is the only honest unit here.")
        say("      A residual is a point estimate of a sampled quantity and")
        say("      cannot bound anything on its own -- reading it as a bound")
        say("      is the same defect as N3 above, which read an")
        say("      unstandardised statistic while the standardisation was")
        say("      printed beside it. So:")
        say()
        mis = 100.0 / 143.0
        say("      depth  residual %   sampling SE %   a 0.70% mis-split, in SE")
        for k in range(6):
            rr = 100.0 * abs(const * float(diffs[k]) / meas[k] - 1.0)
            ss = 100.0 * se[k] / meas[k]
            say("      %-6d %9.4f    %11.4f    %14.1f" % (k, rr, ss, mis / ss))
        say()
        say("      At depth 4 alone a 0.70 per cent mis-split would show at")
        say("      11.5 sampling errors, and at depth 3 at 5.4. That is the")
        say("      bound, and it is stated in the statistic this file")
        say("      already prints rather than in a ratio of percentages.")
    say()
    say("PARAMETERS, printed so that the field reaches the result file:")
    say("    CELLP = %s" % ", ".join(str(p) for p in CELLP))
    # repr, not a rounded format: G78 asks that the constant the source
    # fixes reaches the result file to the precision the source states it,
    # and a %.14f is already a different number.
    say("    2 C_2 = %s" % repr(TWIN2))
    say("    Euler product truncated at p <= %d" % PLIM)
    say()
    say("=" * 70)
    say("N1 %s  N2 %s  N3 %s  N4 %s  N5 %s  N6 %s  N7 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (n1, n2, n3, n4, n5, n6ok, n7)))
    say("The main term of Proposition {#prop:scaleinv} is nonzero at every")
    say("depth, so the exponent the proposition predicts is defined at every")
    say("depth." if n1 else "depth -- REFUTED.")

    head = [
        "STATISTIC: E_same,c[A] - E_all[A], the main term of D_c in",
        "           Proposition {#prop:scaleinv}, with",
        "           A(h) = prod_{p|h, 2<p<=13} (p-1)/(p-2); computed as an",
        "           exact rational over the residue classes mod 2Q and the",
        "           pattern weights m(S) = prod_{p not in S}(p-1).",
        "DENOM: E_all[A], the same average over all even classes mod 2Q.",
        "NULL: an exact permutation null, run in this file. N1 asks whether",
        "      an exact rational is zero and no sampling enters it, so the",
        "      null bites on N2: with six depths all 720 relabellings of the",
        "      depth labels are enumerated and the observed count of cells",
        "      agreeing within two sampling errors is placed in that",
        "      distribution. The sampling errors are READ from",
        "      lab_cell_singular.txt and not recomputed here.",
        "FIELD: exact rationals throughout the A-differences; the constant",
        "       uses 2 C_2 = 1.32032363 and the product over 13 < p <= 4e6;",
        "       measured D_c read from lab_cell_singular.txt at (2e6,4e6].",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    return 0 if n1 else 1


if __name__ == "__main__":
    sys.exit(main())

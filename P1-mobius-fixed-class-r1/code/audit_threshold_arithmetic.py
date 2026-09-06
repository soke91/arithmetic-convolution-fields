# -*- coding: utf-8 -*-
r"""
The threshold the whole sufficient condition is measured against is
not a constant, and the N-sweep never varied it.

WHAT IS AT STAKE

Proposition {#prop:onesided} makes binary Goldbach follow from
|E_3(N)| < S(N)(1 - A(N)) N, and every measurement of that condition in
this repository compares |E_3|/N with the single number 0.3745 --
lab_extend_range.py hardcodes it, Remark {#rem:relocate} prints it,
Remark {#rem:extendrange} extrapolates a crossing against it.

But S(N) = 2 C_2 prod_{p|N, p>2} (p-1)/(p-2) and
A(N) = prod_{q not | N} (1 - 1/(q(q-1))) both depend on which primes
divide N.  The eight N of the published sweep are 2e5 * 2^j, and
2e5 = 2^6 * 5^5, so all eight have odd part 5^5 and therefore the same
threshold, by construction.  The sweep varies the size of N by a factor
128 and the arithmetic of N not at all.

That is a confound, not an error: the numbers are right for the N they
were computed at.  The question is whether the verdict on
[eq:onesided] survives changing N's arithmetic at fixed size, and it
has never been asked.

BACKS: Remark {#rem:threshfam} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  All eight N of the published sweep have odd radical 5, so the
      threshold is one number across the sweep, and recomputing it from
      Euler products reproduces the published 0.3745 at four places.
  V2  The threshold is not a constant of the problem: across test N of
      comparable size but different factorisations it spans at least
      20% of its own value.
  V3  The published family sits at the favourable end: its threshold is
      at or above the median of the test set.
  V4  The verdict flips with arithmetic: at fixed size near 1.6e6, at
      least one test N satisfies [eq:onesided] against its own
      threshold and at least one fails.

REFUTATION RULE (fixed before the run)

  V1  REFUTED if any of the eight has a different odd radical, or if
      the recomputed threshold differs from 0.3745 at four places.
      A failure here means the published constant is wrong.
  V2  REFUTED if the spread is under 20% -- in which case the confound
      is harmless and treating the threshold as constant costs nothing.
  V3  REFUTED if the published family's threshold is below the median.
      A refutation is informative and not bad: it would mean the sweep
      was run on an unfavourable family and the published verdict is
      conservative.
  V4  REFUTED if every test N passes or every test N fails, in which
      case the verdict is arithmetic-independent at this size and the
      concern raised here is empty.

  All four gate.

  THE CONTROL.  V4's comparison is run twice: once against each N's own
  threshold, and once against the constant 0.3745 that the repository
  has been using.  The second is the null -- it is what the published
  practice would have concluded -- and the finding is the difference
  between the two verdict columns, not either column alone.
"""

import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "audit_threshold_arithmetic.txt")

SWEEP = [200_000 * 2 ** j for j in range(8)]
PUBLISHED = 0.3745
PBOUND = 10_000_000       # Euler products truncated here
THETA = 0.56

# comparable size, deliberately different arithmetic
TESTN = [
    1_572_864,   # 2^19 * 3
    1_404_928,   # 2^12 * 7^3
    1_600_000,   # 2^9  * 5^5   -- the published family
    1_620_000,   # 2^5 * 3^4 * 5^4
    1_531_530,   # 2 * 3^2 * 5 * 7 * 11 * 13 * 17
    1_621_620,   # 2^2 * 3^4 * 5 * 7 * 11 * 13
    1_600_006,   # 2 * 800003
]


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def sieves(n):
    pr = primes_upto(n)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(n + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lgp[i]
            if q > n // p:
                break
            q *= p
    mu = np.ones(n + 1, dtype=np.int8)
    rem = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            rem[q::q] //= p
            if q > n // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return pr, lam, mu


def phi_of(k):
    v, phi, d = k, 1, 2
    while d * d <= v:
        if v % d == 0:
            phi *= (d - 1)
            v //= d
            while v % d == 0:
                phi *= d
                v //= d
        d += 1
    if v > 1:
        phi *= (v - 1)
    return phi


def oddprimes(n):
    v, out, d = n, [], 2
    while d * d <= v:
        if v % d == 0:
            if d > 2:
                out.append(d)
            while v % d == 0:
                v //= d
        d += 1
    if v > 1 and v > 2:
        out.append(v)
    return out


def allprimes(n):
    v, out, d = n, [], 2
    while d * d <= v:
        if v % d == 0:
            out.append(d)
            while v % d == 0:
                v //= d
        d += 1
    if v > 1:
        out.append(v)
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("Euler products truncated at %d ..." % PBOUND)
    pr = primes_upto(PBOUND).astype(np.float64)
    odd = pr[1:]
    C2 = float(np.exp(np.log1p(-1.0 / (odd - 1.0) ** 2).sum()))
    ART = float(np.exp(np.log1p(-1.0 / (pr * (pr - 1.0))).sum()))
    # both tails are dominated by sum_{p>P} 2/p^2 < 2/(P log P)
    tail = 2.0 / (PBOUND * math.log(PBOUND))
    say("  twin-prime constant C_2 = %.10f" % C2)
    say("  Artin's constant       = %.10f" % ART)
    say("  truncation tail on either product < %.3e" % tail)

    def thresh(N):
        s = 2.0 * C2
        for p in oddprimes(N):
            s *= (p - 1.0) / (p - 2.0)
        a = ART
        for q in allprimes(N):
            a /= (1.0 - 1.0 / (q * (q - 1.0)))
        return s, a, s * (1.0 - a)

    say()
    say("V1  the published sweep, N = 2e5 * 2^j")
    say("  N            odd primes of N      S(N)      A(N)      "
        "threshold")
    rad = set()
    tv = []
    for N in SWEEP:
        op = oddprimes(N)
        rad.add(tuple(op))
        s, a, t = thresh(N)
        tv.append(t)
        say("  %-12d %-20s %-9.4f %-9.4f %.6f"
            % (N, ",".join(map(str, op)) or "(none)", s, a, t))
    v1 = len(rad) == 1 and round(tv[0], 4) == PUBLISHED
    say("  distinct odd radicals across the sweep: %d" % len(rad))
    say("  recomputed threshold %.6f against the published %.4f   %s"
        % (tv[0], PUBLISHED, "hold" if v1 else "REFUTED"))

    say()
    say("measuring at comparable size and different arithmetic ...")
    NMAX = max(TESTN)
    _, lam, mu = sieves(NMAX)

    rows = []
    for N in TESTN:
        s, a, t = thresh(N)
        PN = set(allprimes(N))
        K = int(N ** THETA)
        ks = np.array([k for k in range(2, K)
                       if mu[k] != 0 and all(k % q for q in PN)])
        lg = np.log(ks.astype(float))
        iph = np.array([phi_of(int(k)) for k in ks], dtype=np.float64)
        f0 = np.zeros(N, dtype=np.float64)
        idx = np.arange(1, N, dtype=np.int64)
        f0[1:] = lam[1:N] * mu[N - idx]
        C = float(f0.sum())
        A = np.empty(ks.size)
        for i, k in enumerate(ks):
            r = N % int(k)
            A[i] = f0[r::int(k)].sum() if r else f0[int(k)::int(k)].sum()
        sgn = mu[ks].astype(np.float64)
        term = sgn * (A - C / iph)
        E3 = float((lg * term).sum())
        B = float((lg * np.abs(term)).sum())
        rows.append((N, oddprimes(N), s, a, t, abs(E3) / N, B / N,
                     abs(C) / N, int((term > 0).sum()),
                     int((term < 0).sum()), int(ks.min())))

    say()
    say("V2/V3  the threshold across the test set")
    say("  N            odd primes of N            threshold")
    for N, op, s, a, t, e, b, c, np_, nn_, kmin in rows:
        say("  %-12d %-26s %.6f"
            % (N, ",".join(map(str, op)) or "(none)", t))
    ts = [r[4] for r in rows]
    lo, hi = min(ts), max(ts)
    spread = (hi - lo) / (0.5 * (lo + hi))
    v2 = spread >= 0.20
    say("  range %.6f to %.6f, spread %.4f of the midpoint (floor 0.20)"
        "   %s" % (lo, hi, spread, "hold" if v2 else "REFUTED"))
    med = float(np.median(ts))
    fam = [r[4] for r in rows if r[0] == 1_600_000][0]
    v3 = fam >= med
    say("  the published family is at %.6f, the median is %.6f   %s"
        % (fam, med, "hold" if v3 else "REFUTED"))

    say()
    say("V4  the verdict on [eq:onesided], each N against its own")
    say("    threshold, and the null: the same N against the constant")
    say("    %.4f the repository has been using" % PUBLISHED)
    say("  N            |E_3|/N    B/N       own thr.   own      "
        "constant")
    npass = nfail = 0
    flips = 0
    for N, op, s, a, t, e, b, c, np_, nn_, kmin in rows:
        own = e < t
        nul = e < PUBLISHED
        npass += own
        nfail += not own
        flips += own != nul
        say("  %-12d %-10.4f %-9.4f %-10.6f %-8s %s"
            % (N, e, b, t, "holds" if own else "FAILS",
               "holds" if nul else "FAILS"))
    v4 = npass > 0 and nfail > 0
    say("  %d of %d hold against their own threshold; the constant"
        % (npass, len(rows)))
    say("  disagrees with the honest verdict at %d of them   %s"
        % (flips, "hold" if v4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC 1 (post hoc). Where the failures come from. For N")
    say("  with several small odd prime factors every term of E_3 has the")
    say("  same sign -- |E_3| = B exactly -- so there is no cancellation")
    say("  across k at all. The reason is visible in the smallest")
    say("  admissible k: k must be coprime to N, and by [eq:dilate] the m")
    say("  that survive are coprime to rad(N) too, so over a short range")
    say("  they are almost all primes and mu(m) = -1 dominates.")
    say("  N            odd primes                 k_min  terms +  -")
    for N, op, s, a, t, e, b, c, np_, nn_, kmin in rows:
        say("  %-12d %-26s %-6d %-8d %d"
            % (N, ",".join(map(str, op)) or "(none)", kmin, np_, nn_))

    say()
    say("  DIAGNOSTIC 2 (post hoc). Where the conservatism is. The")
    say("  threshold is S(1-A) only because [prop:onesided] bounds the")
    say("  wall by |C(N)| <= A(N)N. Measured, that bound is slack by two")
    say("  to three orders. Putting the measured wall in its place,")
    say("  Theorem [thm:C] asks for |E_3|/N < S(N)(1 - |C(N)|/N):")
    say("  N            |C|/N      A(N)     |C|/(A N)  S(1-|C|/N)  "
        "|E_3|/N   verdict")
    allpass = True
    for N, op, s, a, t, e, b, c, np_, nn_, kmin in rows:
        t2 = s * (1.0 - c)
        if e >= t2:
            allpass = False
        say("  %-12d %-10.4f %-8.4f %-10.4f %-11.4f %-9.4f %s"
            % (N, c, a, c / a, t2, e, "holds" if e < t2 else "FAILS"))
    say("  every test N %s the measured-wall condition, so the arithmetic"
        % ("satisfies" if allpass else "does not satisfy"))
    say("  dependence of the verdict lives entirely in the slack of")
    say("  |C| <= A N, not in E_3. This is a measurement, not a proof:")
    say("  |C(N)| = o(N) is the open problem, and nothing here supplies")
    say("  it. What it locates is which of the two inequalities in")
    say("  [prop:onesided] is doing the damage.")

    say()
    say("=" * 70)
    ok = v1 and v2 and v3 and v4
    say("the sufficient condition's threshold moves with N's arithmetic "
        "and the sweep held it fixed" if ok else "REFUTED")

    head = [
        "STATISTIC: S(N) = 2 C_2 prod_{p|N,p>2}(p-1)/(p-2) and",
        "           A(N) = prod_{q not| N}(1 - 1/(q(q-1))) from Euler",
        "           products truncated at 1e7; the threshold S(N)(1-A(N));",
        "           |E_3(N)|/N and B(N)/N at seven N of comparable size and",
        "           different factorisation; and the verdict on",
        "           [eq:onesided] against each N's own threshold and",
        "           against the constant 0.3745.",
        "NULL: the constant 0.3745. It is what the repository has been",
        "      comparing against, so running the same verdict against it",
        "      is exactly the counterfactual 'what would the published",
        "      practice have concluded', and the claim here is the",
        "      difference between the two verdict columns. A sign control",
        "      is not applicable: S and A are Euler products with no",
        "      randomisable ingredient.",
        "FIELD: the published sweep N = 2e5 * 2^j, j = 0..7; seven test N",
        "       in [1.40e6, 1.63e6] chosen for their odd radicals, from a",
        "       pure 2-power-times-3 to a product of six odd primes;",
        "       theta' = 0.56; k over the squarefree k < N^0.56 coprime",
        "       to N.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not ok:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())

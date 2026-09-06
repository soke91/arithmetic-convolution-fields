# -*- coding: utf-8 -*-
r"""
OPEN.md, second half of "목표보다 강한 것을 공격하고 있는가" -- how much
saving Proposition {#prop:onesided} actually asks for.

WHAT THIS MEASURES AND WHY

Proposition [prop:onesided] says binary Goldbach follows from the
one-sided bound E_3 > -(1-eps) S(N)(1-A(N))N for a single fixed
eps > 0.  The quantifier is the point and this docstring had it wrong:
it read (1+o(1)) in place of (1-eps), which is the form the
proposition names and rejects -- unquantified, the hypothesis is not
a closed condition.  Nothing this script measures turns on it: the
saving below is computed against the threshold S(N)(1-A(N))N itself
and a constant factor does not enter.  That is weaker than the
consumed |E_3| <<_A N(log N)^{-A}, but "weaker" is not a currency.  The
currency in this subject is the saving over the unconditional bound, so
this script measures that.

The unconditional bound is the triangle inequality on the k-sum,

    B(N) := sum_{k<K, (k,N)=1} (log k) |E_mu(N;k)|,

which is what Huang-Li take before appealing to EH_mu.  Brun-Titchmarsh
gives |E_mu(N;k)| << N/phi(k) and sum_{k<K} 1/phi(k) << log K, so
B(N) << N (log N)^2 -- and that is the only unconditional input.  Then:

  * the consumed bound needs B(N) beaten by a factor (log N)^{A+2} for
    every A, i.e. by every power of log;
  * [eq:onesided] needs B(N) beaten only by
    B(N) / (S(N)(1-A(N))N), which is ~ (log N)^2 / (S(1-A)) -- a FIXED
    power, about 2 for almost all N and about 3 at primorials.

That is the concrete sense in which the direct route is easier, and it
is what this script pins down.  It also checks the awkward part: the
proposition is asymptotic, so at accessible N the condition may fail.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  G1  B(N)/(N (log N)^2) lies in [0.02, 0.5] at every N tested and
      varies by less than a factor 2 across the range, i.e. B ~ N(log N)^2
      and Brun-Titchmarsh is the right scale.
  G2  |E_3(N)| / B(N) < 0.05 at every N: keeping the signs mu(k) buys at
      least a factor 20 over the triangle bound Huang-Li discard.
  G3  Fitting the required saving B(N)/(S(N)(1-A(N))N) as (log N)^b
      gives b in [1.8, 2.6] -- a fixed power near 2.
  G4  E_3(N) > -S(N)(1-A(N))N holds at every N tested with N >= 10^6,
      and fails for at least one N <= 4*10^5.  The proposition is
      asymptotic and the accessible range straddles its crossover.

REFUTATION RULE (fixed before the run)

  G1  REFUTED if any value leaves [0.02, 0.5], or if max/min > 2.
  G2  REFUTED if any ratio is at least 0.05.
  G3  REFUTED if b leaves [1.8, 2.6].
  G4  REFUTED if the condition fails at any N >= 10^6, or if it holds at
      every N <= 4*10^5.

  All four gate.  G4 is deliberately two-sided: a route that could not
  be seen to fail anywhere accessible would be untested, and one that
  failed at the top would be refuted outright.

BACKS: Proposition {#prop:nolog} in paper/theorem_A.md -- that the
demand is a constant-factor bound on B(N) and asks no saving in
log N -- and Remark {#rem:relocate}.
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
OUT = os.path.join(ROOT, "results", "lab_onesided_demand.txt")

THETA = 0.56
NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
PLIM = 4_000_000


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


def factor_set(n):
    out, d = set(), 2
    while d * d <= n:
        if n % d == 0:
            out.add(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.add(n)
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("sieving to %d ..." % PLIM)
    pr, lam, mu = sieves(PLIM)

    artin, twin = 1.0, 2.0
    for p in pr:
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    say()
    say("  N          K      S(N)     1-A(N)   thresh/N   E_3/N      "
        "B/N        B/(N log^2 N)  |E_3|/B   saving")
    say("  " + "-" * 104)
    rows = []
    for N in NS:
        PN = factor_set(N)
        A = artin
        S = twin
        for q in sorted(PN):
            A /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))
        K = int(N ** THETA)

        n = np.arange(N, dtype=np.int64)
        f = np.zeros(N, dtype=np.float64)
        f[1:] = lam[1:N] * mu[(N - n)[1:]]
        C = float(f.sum())

        E3 = 0.0
        B = 0.0
        for k in range(2, K):
            if mu[k] == 0:
                continue
            v, phi, ok = k, 1, True
            while v > 1:
                p = int(pr[np.searchsorted(pr, 2)]) if False else 0
                break
            # explicit factorisation of k against P(N)
            ok = all(k % q for q in PN)
            if not ok:
                continue
            phi, v = 1, k
            d = 2
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
            r = N % k
            inner = float(f[r::k].sum()) if r else float(f[k::k].sum())
            e = inner - C / phi
            E3 += mu[k] * math.log(k) * e
            B += math.log(k) * abs(e)

        thr = S * (1.0 - A) * N
        L = math.log(N)
        rows.append((N, K, S, 1.0 - A, thr / N, E3 / N, B / N,
                     B / (N * L * L), abs(E3) / B, B / thr))
        say("  %-10d %-6d %-8.4f %-8.4f %-10.4f %+-10.4f %-10.2f "
            "%-14.4f %-9.5f %.1f"
            % rows[-1])

    say()
    b_scale = [r[7] for r in rows]
    g1 = (all(0.02 <= v <= 0.5 for v in b_scale)
          and max(b_scale) / min(b_scale) <= 2.0)
    say("G1  B/(N log^2 N) in [%.4f, %.4f], spread %.3f   (band "
        "[0.02,0.5], spread cap 2)   %s"
        % (min(b_scale), max(b_scale), max(b_scale) / min(b_scale),
           "hold" if g1 else "REFUTED"))

    r2 = [r[8] for r in rows]
    g2 = all(v < 0.05 for v in r2)
    # The registered cap is printed with the verdict.  It stood only in
    # this file's docstring, so a reader of the .txt could not see what
    # G2 had been registered at -- the band is the evidence, not the
    # verdict word.
    say("G2  |E_3|/B in [%.5f, %.5f]; the signs buy a factor %.0f to %.0f"
        "   (cap 0.05, i.e. a factor 20)   %s"
        % (min(r2), max(r2), 1.0 / max(r2), 1.0 / min(r2),
           "hold" if g2 else "REFUTED"))

    xs = np.log(np.log(np.array([r[0] for r in rows], dtype=float)))
    ys = np.log(np.array([r[9] for r in rows], dtype=float))
    b = float(np.polyfit(xs, ys, 1)[0])
    g3 = 1.8 <= b <= 2.6
    say("G3  required saving B/thresh = %s"
        % ", ".join("%.1f" % r[9] for r in rows))
    say("    fitted as (log N)^b with b = %.4f   (band [1.8,2.6])   %s"
        % (b, "hold" if g3 else "REFUTED"))

    ok_hi = all(r[5] > -r[4] for r in rows if r[0] >= 1_000_000)
    fail_lo = any(r[5] <= -r[4] for r in rows if r[0] <= 400_000)
    g4 = ok_hi and fail_lo
    say("G4  E_3/N vs -thresh/N:")
    for r in rows:
        say("      N = %-10d  E_3/N = %+.4f   threshold = %+.4f   %s"
            % (r[0], r[5], -r[4],
               "satisfied" if r[5] > -r[4] else "FAILS"))
    say("    holds at every N >= 1e6: %s;  fails somewhere <= 4e5: %s"
        "   %s" % (ok_hi, fail_lo, "hold" if g4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc, not pre-registered). G1-G3 are refuted")
    say("  because Brun-Titchmarsh is a gross overestimate of B: the")
    say("  provable |E_mu| << N/phi(k) gives B << N (log N)^2, but the")
    say("  measured B is ~ 0.6 N, i.e. the truth already carries the two")
    say("  powers of log. So the right reading of the numbers is not")
    say("  'saving over Brun-Titchmarsh' but 'how far the triangle bound")
    say("  itself is from the threshold', and that is B/thresh:")
    Ls = np.log(np.array([r[0] for r in rows], dtype=float))
    bn = np.array([r[6] for r in rows])
    say("  B/N   = %s" % ", ".join("%.4f" % v for v in bn))
    say("  ratio = %s" % ", ".join("%.3f" % r[9] for r in rows))
    sl = float(np.polyfit(np.log(Ls), np.log(bn), 1)[0])
    say("  B/N fitted as (log N)^%.3f" % sl)
    thr_n = rows[-1][4]
    if sl < 0:
        need = math.exp(math.exp(math.log(thr_n / bn[-1]) / sl + math.log(Ls[-1])))
        say("  at that rate B/N reaches the threshold %.4f near "
            "N = 10^%.1f" % (thr_n, math.log10(need)))
    say("  Below that point the triangle inequality alone gives")
    say("  r~(N) > 0, with no cancellation in the k-sum used at all.")

    say()
    say("  For contrast, the saving the consumed bound asks for is")
    say("  B/(N (log N)^{-A}) = (log N)^{A+2} times a constant, i.e.")
    for r in rows[-1:]:
        L = math.log(r[0])
        say("  at N = %d: %.1f for the one-sided threshold, against"
            % (r[0], r[9]))
        say("  %.0f at A=1, %.0f at A=2, %.0f at A=3 for the consumed one."
            % (r[6] * L ** 1, r[6] * L ** 2, r[6] * L ** 3))

    say()
    say("=" * 70)
    ok = g1 and g2 and g3 and g4
    say("G1 %s  G2 %s  G3 %s  G4 %s"
        % tuple("hold" if v else "REFUTED" for v in (g1, g2, g3, g4)))
    say("the one-sided demand is a fixed power of log over the "
        "unconditional bound" if ok else "REFUTED")

    head = [
        "STATISTIC: B(N) = sum_{k<K,(k,N)=1} mu^2(k)(log k)|E_mu(N;k)|,",
        "           the squarefree restriction being what the triangle",
        "           inequality leaves, since mu(k) = 0 elsewhere; it is",
        "           imposed in the loop and was missing from this line;",
        "           the quantity is unchanged. This is the",
        "           triangle bound Huang-Li take before appealing to",
        "           EH_mu; B/(N (log N)^2), which Brun-Titchmarsh predicts",
        "           to be bounded; |E_3(N)|/B(N), the gain from keeping",
        "           the signs mu(k); the required saving",
        "           B(N)/(S(N)(1-A(N))N) and its exponent in log N; and",
        "           E_3(N)/N against the threshold -S(N)(1-A(N)).",
        "FIELD: N = 2e5, 4e5, 8e5, 1.6e6, 3.2e6 with theta' = 0.56, so",
        "       K = floor(N^0.56); E_mu(N;k) by direct enumeration of the",
        "       progression n = N mod k; Lambda and mu from an integer",
        "       sieve to 4e6; A(N) and S(N) as Euler products over p<4e6.",
        'NULL: run, in lab_level_coin_null.py. A coin on the same support',
        '      gives a SMALLER B(N;K) at the same K, so the smallness of B',
        '      measured here is not evidence about mu. The proposition it',
        '      supports is an implication and does not rest on that.',
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

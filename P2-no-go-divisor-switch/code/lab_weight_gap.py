# -*- coding: utf-8 -*-
r"""
What separates the proved theorem from the open one, at the level of
the dilated walls.

WHAT IS AT STAKE

Proposition {#prop:dilate} gives E_mu(N;k) = mu(k)H(N;k) - C(N)/phi(k),
and Proposition {#prop:posweights} used it on the log k branch.  Used
on the OTHER branch it does the same thing to Theorem {#thm:A}:

    T_1(N) = sum_{k<K,(k,N)=1} mu(k) E_mu(N;k)
           = sum_{k<K,(k,N)=1} H(N;k) - C(N) B(K),
    B(K)   = sum_{k<K,(k,N)=1} mu(k)/phi(k) << exp(-c sqrt(log K)).

So Theorem [thm:A], which is proved unconditionally, says exactly that
the FLAT sum of dilated walls is small, while E_3 -- which is the wall
itself, by Theorem {#thm:C} -- is the same sum weighted by log k.
Both weights are nonnegative.  The whole distance between what is
proved and what is open is the factor log k inside a positively
weighted sum of the same terms.

That is worth measuring rather than asserting: how much smaller is the
flat sum, and how does the smallness degrade as the weight is turned
on through (log k)^j?

BACKS: Proposition {#prop:flatsum} and Remark {#rem:weightgap} in
paper/theorem_A.md.

Z1 has content only if the two sides are built by different routes.
The left side is T_1 = sum mu(k) (A(N;k) - C/phi(k)) with A(N;k) the
progression sum over n = N (mod k), read off one array.  The right
side needs H(N;k) from its OWN definition -- the m-sum of the dilate,
H(N;k) = sum_{m<N/k,(m,k)=1} Lambda(N-mk) mu(m) -- and not as
mu(k)A(N;k), because sum mu(k)(A - C/phi) - (sum mu(k)A - C sum
mu(k)/phi) is the distributive law and vanishes whatever A, C, phi
and mu contain.  It was written that way, and could not fail.  The
m-sum is now formed per k from Lambda at N - mk and mu at m with the
coprimality imposed explicitly; the proposition's identity
A(N;k) = mu(k)H(N;k) is what the residual now measures.  Dropping the
coprimality, shifting m by one, or signing the m-sum by mu(k) each
refutes Z1 and leaves Z2, Z3, Z4 as they were.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Z1  The identity T_1(N) = sum H - C(N)B(K) holds to better than
      1e-12 relative at every N.
  Z2  The flat sum is far smaller than the log-weighted one:
      |sum H| / |sum (log k) H| < 0.2 at every N.
  Z3  |sum H|/N decays faster than |E_3|/N: fitting both against
      log N, the flat sum's exponent is the steeper.
  Z4  Turning the weight on is monotone: |sum (log k)^j H|/N increases
      with j over j = 0, 0.25, 0.5, 0.75, 1 at every N.

REFUTATION RULE (fixed before the run)

  Z1  REFUTED if the relative error reaches 1e-12 at any N.  It is an
      identity, so a failure is an error in the derivation.
  Z2  REFUTED if the ratio reaches 0.2 at any N.
  Z3  REFUTED if the flat sum's fitted exponent is not steeper.
  Z4  REFUTED if the sequence in j is not increasing at any N.

  All four gate.
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
OUT = os.path.join(ROOT, "results", "lab_weight_gap.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000,
      12_800_000, 25_600_000]
JS = [0.0, 0.25, 0.5, 0.75, 1.0]
THETA = 0.56


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


def loo(x, y, name, say):
    """Refit dropping each end in turn, and report the spread.

    Every exponent this repository quotes is a slope over four to eight
    values of N, and audit_truncation_exponent.py showed what such a
    slope is worth when nobody varies the free parameter that defines
    it. For a direct fit the free parameter is the N-range, so the
    cheapest honest check is to refit without the smallest N and
    without the largest and print how far the answer moves.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    f = [float(np.polyfit(x[s], y[s], 1)[0])
         for s in (slice(None), slice(1, None), slice(0, -1))]
    sp = max(f) - min(f)
    say("  leave-one-out on %s: full %.4f, without the smallest N "
        "%.4f," % (name, f[0], f[1]))
    say("  without the largest %.4f -- spread %.4f" % (f[2], sp))
    say("SWEPT %s N-range %.4f" % (name, sp))
    return sp


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)

    say()
    say("  N            |T_1|/N    |sum H|/N   |E_3|/N   ratio    "
        "rel err of the identity")
    say("  " + "-" * 82)
    flat, e3s, jtab = [], [], []
    z1 = z2 = z4 = True
    for N in NS:
        v, PN, d = N, set(), 2
        while d * d <= v:
            if v % d == 0:
                PN.add(d)
                while v % d == 0:
                    v //= d
            d += 1
        if v > 1:
            PN.add(v)
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
        H = sgn * A

        # H(N;k) by its own definition, the m-sum of the dilate: Lambda
        # at N - mk, mu at m, (m,k) = 1 imposed by hand.  It reads
        # neither f0 nor A.  m runs downward so that the terms are met
        # in the same order as the progression slice above: the
        # residual then measures the identity mu(k)H = A and not the
        # order of floating-point addition.
        Hm = np.empty(ks.size)
        for i, k in enumerate(ks):
            k = int(k)
            ms = np.arange(N // k, 0, -1, dtype=np.int64)
            s = mu[ms].astype(np.float64)
            v, p = k, 2
            while p * p <= v:
                if v % p == 0:
                    s[ms % p == 0] = 0.0
                    while v % p == 0:
                        v //= p
                p += 1
            if v > 1:
                s[ms % v == 0] = 0.0
            Hm[i] = float((lam[N - k * ms] * s).sum())

        T1 = float((sgn * (A - C / iph)).sum())
        BK = float((sgn / iph).sum())
        rhs = float(Hm.sum()) - C * BK
        rel = abs(T1 - rhs) / max(abs(T1), 1e-300)
        if rel >= 1e-12:
            z1 = False

        S0 = float(H.sum())
        S1 = float((lg * H).sum())
        E3 = float((sgn * lg * (A - C / iph)).sum())
        r = abs(S0) / max(abs(S1), 1e-300)
        if r >= 0.2:
            z2 = False

        row = [abs(float(((lg ** j) * H).sum())) / N for j in JS]
        if not all(row[i] < row[i + 1] for i in range(len(row) - 1)):
            z4 = False
        jtab.append((N, row))

        flat.append(abs(S0) / N)
        e3s.append(abs(E3) / N)
        say("  %-12d %-10.5f %-11.5f %-9.4f %-8.4f %.3e"
            % (N, abs(T1) / N, abs(S0) / N, abs(E3) / N, r, rel))

    say()
    say("Z1  the identity holds at every N                     %s"
        % ("hold" if z1 else "REFUTED"))
    say("Z2  |sum H| / |sum (log k) H| < 0.2 at every N        %s"
        % ("hold" if z2 else "REFUTED"))

    x = np.log(np.array(NS, dtype=float))
    bf = -float(np.polyfit(x, np.log(np.array(flat)), 1)[0])
    be = -float(np.polyfit(x, np.log(np.array(e3s)), 1)[0])
    z3 = bf > be
    say("Z3  |sum H|/N ~ N^{-%.4f} against |E_3|/N ~ N^{-%.4f}   %s"
        % (bf, be, "hold" if z3 else "REFUTED"))
    loo(x, np.log(np.array(flat)), "flat_sum_decay", say)
    loo(x, np.log(np.array(e3s)), "E3_decay", say)

    say()
    say("Z4  turning the weight on: |sum (log k)^j H| / N")
    say("  N            " + "  ".join("j=%.2f  " % j for j in JS))
    for N, row in jtab:
        say("  %-12d %s" % (N, "  ".join("%-8.5f" % v for v in row)))
    say("  Z4 %s" % ("hold" if z4 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc). The rows are close to geometric in j,")
    say("  which is what a sum concentrated at ONE effective modulus")
    say("  would give: sum (log k)^j H ~ (log k*)^j sum H. Reading k*")
    say("  off the ratio, and comparing it with the truncation point:")
    say("  N            ratio per 0.25 in j   spread   log k*   k*      "
        "log k* / log K")
    for N, row in jtab:
        rr = [row[i + 1] / row[i] for i in range(len(row) - 1)]
        g = float(np.mean(rr))
        lk = g ** 4.0
        K = int(N ** THETA)
        say("  %-12d %-21.4f %-8.4f %-8.4f %-7.0f %.4f"
            % (N, g, max(rr) / min(rr), lk, math.exp(lk),
               lk / math.log(K)))
    say("  k* itself is unstable -- it is the exponential of a fourth")
    say("  power of the fitted ratio -- but log k* / log K is not, and it")
    say("  says k* = K^{0.8 to 0.94}. The log weight does not spread the")
    say("  sum; it reweights it to a definite place near the top of the")
    say("  k-range, which is where Lemma {#lem:extract} says a weight")
    say("  must sit to see C(N) at all.")

    say()
    say("=" * 70)
    ok = z1 and z2 and z3 and z4
    say("Theorem {#thm:A} is the flat sum of dilated walls and E_3 is "
        "the same sum weighted by log k" if ok else "REFUTED")

    head = [
        "STATISTIC: |T_1(N)|/N; the flat sum |sum_{k<K} H(N;k)|/N; the",
        "           log-weighted |E_3(N)|/N; their ratio; the relative",
        "           error of the identity T_1 = sum H - C(N)B(K); fitted",
        "           decay exponents of the flat sum and of |E_3|/N; and",
        "           |sum (log k)^j H|/N for j = 0, 0.25, 0.5, 0.75, 1.",
        "DENOM: the identity residual is divided by T_1(N);",
        "NULL: none applies. Z1 is an identity and the rest are",
        "      comparisons between two weightings of the SAME numbers, so",
        "      a sign control would change both sides alike and settle",
        "      nothing; the field's coin controls were run in",
        "      lab_sign_structure.py and lab_lean_decay.py.",
        "FIELD: N = 2e5 through 2.56e7 by doubling, theta' = 0.56; k over",
        "       the squarefree k < N^0.56 coprime to N; m over 1 <= m <",
        "       N/k with (m,k) = 1; Lambda and mu from an integer sieve to",
        "       2.56e7.",
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

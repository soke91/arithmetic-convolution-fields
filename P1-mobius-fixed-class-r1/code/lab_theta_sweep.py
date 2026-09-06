# -*- coding: utf-8 -*-
r"""
OPEN.md, "접근 가능한 범위에서 thm:C는 부호만 시험된다" -- can the
error term be pushed under the signal by choosing theta'?

WHAT IS AT STAKE

Theorem {#thm:C} says

    E_3(alpha) = r~(N) - S(N)(N - C(N)) + O_A(N (log N)^{-A}),

and audit_E3_constant measured, at theta' = 0.56, that the residual

    R(N, theta') := | E_3(N; theta') - ( r~(N) - S(N)(N - C(N)) ) | / N

is about 0.26 at N = 1.6e6 while the right-hand side it is compared
against is about 0.001 -- the error swamps the signal by two orders,
so nothing accessible tests the identity's content, only its sign.
That measurement fixed theta' = 0.56.  But theta' is a free parameter
of the Corollary-1 regime, and the error terms pull against each other
in it:

  * B_log(K) = -S(N) + O(exp(-c sqrt(log K)) log K) improves as
    K = N^{theta'} grows, i.e. wants theta' large;
  * the main-term cancellation runs on sum_{m<M} f(m)/m with
    M = N^{1-theta'}, and is O(exp(-c sqrt(log M))), i.e. wants
    theta' small;
  * Bombieri-Vinogradov needs the level N^{1-theta'} below N^{1/2},
    i.e. theta' > 1/2 at all.

So there should be an interior optimum, and the question is whether it
buys enough to make the identity testable.  Nothing in either paper
sweeps theta'.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  I1  R(N, theta') has an interior minimum in theta' in (0.5, 1) at
      every N tested -- it is not monotone in theta'.
  I2  The minimising theta' is at least 0.70 at every N tested.
  I3  At the minimising theta', R < 0.10 at every N tested.
  I4  At theta' = 0.56, R reproduces the earlier audit: 0.456, 0.373,
      0.311 at N = 2e5, 4e5, 8e5.
  I5  At the minimising theta', R decreases with N across the three N.

REFUTATION RULE (fixed before the run)

  I1  REFUTED if R is monotone in theta' at any N, i.e. if the argmin
      sits at an endpoint of the swept range.
  I2  REFUTED if the minimising theta' is below 0.70 at any N.
  I3  REFUTED if R at the minimum is at least 0.10 at any N.
  I4  REFUTED if any of the three differs from the published value by
      more than 0.01.
  I5  REFUTED if R at the minimum fails to decrease between consecutive
      N.

  All five gate.  I4 is the tie to the existing audit and must hold for
  the rest to mean anything; I1-I3 are the question; I5 is what decides
  whether the accessible-range limitation is a fact about theta' or a
  fact about N.

CITED BY: {#rem:thetasweep} in paper/.
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
OUT = os.path.join(ROOT, "results", "lab_theta_sweep.txt")

NS = [200_000, 400_000, 800_000]
THETAS = [0.51, 0.55, 0.56, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
PUB56 = {200_000: 0.456, 400_000: 0.373, 800_000: 0.311}
PLIM = 1_000_000
# The Euler products below must not be tied to the
# measurement range: audit_constants.py shows the
# truncation reaches the sixth printed decimal.
CLIM = 4_000_000


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


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("sieving to %d ..." % PLIM)
    pr, lam, mu = sieves(PLIM)
    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    results = {}
    sigs = {}
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
        S = twin
        for q in sorted(PN):
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))

        n = np.arange(N, dtype=np.int64)
        f = np.zeros(N, dtype=np.float64)
        f[1:] = lam[1:N] * mu[(N - n)[1:]]
        C = float(f.sum())
        rt = float((lam[1:N] * lam[(N - n)[1:]]).sum())
        target = rt - S * (N - C)
        sigs[N] = target / N

        say()
        say("N = %d   S(N) = %.6f   C(N)/N = %+.6f   "
            "(r~ - S(N-C))/N = %+.8f" % (N, S, C / N, target / N))
        say("  theta'   K          E_3/N        residual R    ")
        say("  " + "-" * 52)
        row = []
        for th in THETAS:
            K = int(N ** th)
            E3 = 0.0
            for k in range(2, K):
                if mu[k] == 0:
                    continue
                if any(k % q == 0 for q in PN):
                    continue
                phi = phi_of(k)
                r = N % k
                inner = float(f[r::k].sum()) if r else float(f[k::k].sum())
                E3 += mu[k] * math.log(k) * (inner - C / phi)
            R = abs(E3 - target) / N
            row.append((th, K, E3 / N, R))
            say("  %-8.2f %-10d %+-12.6f %.6f" % (th, K, E3 / N, R))
        results[N] = row

    say()
    say("=" * 70)
    i1 = i2 = i3 = True
    mins = {}
    for N in NS:
        row = results[N]
        Rs = [r[3] for r in row]
        j = int(np.argmin(Rs))
        mins[N] = (row[j][0], Rs[j])
        interior = 0 < j < len(Rs) - 1
        i1 = i1 and interior
        i2 = i2 and (row[j][0] >= 0.70)
        i3 = i3 and (Rs[j] < 0.10)
        say("N = %-8d argmin theta' = %.2f  R = %.6f  interior: %s"
            % (N, row[j][0], Rs[j], interior))
    say("I1 %s   I2 %s   I3 %s"
        % ("hold" if i1 else "REFUTED", "hold" if i2 else "REFUTED",
           "hold" if i3 else "REFUTED"))

    got56 = {}
    for N in NS:
        for th, K, e3, R in results[N]:
            if abs(th - 0.56) < 1e-9:
                got56[N] = R
    e4 = max(abs(got56[N] - PUB56[N]) for N in NS)
    i4 = e4 <= 0.01
    say("I4  R at theta'=0.56: %s   against %s   max dev %.6f   %s"
        % (", ".join("%.4f" % got56[N] for N in NS),
           ", ".join("%.3f" % PUB56[N] for N in NS), e4,
           "hold" if i4 else "REFUTED"))

    seq = [mins[N][1] for N in NS]
    i5 = all(seq[i] > seq[i + 1] for i in range(len(seq) - 1))
    say("I5  R at the argmin by N: %s   decreasing: %s   %s"
        % (", ".join("%.6f" % s for s in seq), i5,
           "hold" if i5 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). I1-I3 fail because R is essentially")
    say("  MONOTONE INCREASING in theta': the finite-N error is dominated")
    say("  by the main-term cancellation over m < M = N^{1-theta'}, which")
    say("  theta' -> 1 destroys, and the gain in B_log(K) is second order.")
    say("  So the best theta' is the smallest admissible one. The question")
    say("  the OPEN item actually asks is whether R falls below the signal")
    say("  it is compared against, and that is a ratio, not a level:")
    say("  N          signal |r~-S(N-C)|/N   R at argmin   ratio")
    rats = []
    for N in NS:
        row = results[N]
        Rs = [r[3] for r in row]
        j = int(np.argmin(Rs))
        sig = abs(sigs[N])
        rats.append(Rs[j] / sig)
        say("  %-10d %-22.8f %-13.6f %.2f" % (N, sig, Rs[j], rats[-1]))
    say("  the ratio %s with N"
        % ("WIDENS" if all(rats[i] < rats[i + 1]
                           for i in range(len(rats) - 1)) else "narrows"))
    say("  because the signal falls faster than the error: signal ratios")
    say("  %s against error ratios %s"
        % (", ".join("%.3f" % (abs(sigs[NS[i + 1]]) / abs(sigs[NS[i]]))
                     for i in range(len(NS) - 1)),
           ", ".join("%.3f" % (seq[i + 1] / seq[i])
                     for i in range(len(seq) - 1))))

    say()
    say("=" * 70)
    ok = i1 and i2 and i3 and i4 and i5
    say("I1 %s  I2 %s  I3 %s  I4 %s  I5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (i1, i2, i3, i4, i5)))
    say("theta' buys enough to bring the identity within reach"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the residual R(N,theta') = |E_3(N;theta')",
        "           - (r~(N) - S(N)(N - C(N)))| / N, with E_3 computed by",
        "           direct enumeration over k < K = floor(N^theta'), swept",
        "           over theta'; the minimising theta' at each N; and the",
        "           value of R there.",
        "FIELD: N = 2e5, 4e5, 8e5; theta' in 0.51, 0.55, 0.56, 0.60, 0.65,",
        "       0.70, 0.75, 0.80, 0.85, 0.90, 0.95; Lambda and mu from an",
        "       integer sieve to 1e6; S(N) as an Euler product over",
        "       p < 4e6, the bound CLIM and deliberately not tied to",
        "       the 1e6 measurement sieve (an earlier printing of this",
        "       line said p < 1e6); r~ and C by direct summation.",
        "NULL: none applies. R(N,theta') is the residual of the identity of",
        '      Theorem {#thm:C}, and a coin satisfies no such identity, so',
        '      there is no control quantity to compare against. The finding',
        "      is a monotonicity in theta' and a ratio trend in N, neither of",
        '      which is a detection claim.',
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

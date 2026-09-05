# -*- coding: utf-8 -*-
"""C-02.  The flat sum of dilated walls, with k = 1 restored.

P2 [eq:Tw]/[prop:flatsum] index every sum by  k < K, (k,N) = 1 -- which
contains k = 1.  lab_weight_gap.py builds its k-list as range(2, K), so
the flat sum sum_k H(N;k) it prints omits H(N;1) = C(N).

T_1 and T_log are unaffected (the k=1 term of T_w is
A(N;1) - C(N)/phi(1) = 0, and log 1 = 0), so this is a defect in one
printed column only.  This script measures how large the omitted term is
and refits the two exponents with it restored.

Independent implementation: mu is built from an omega-counter plus a
squarefull mask (not from the sign-flip + leftover-cofactor recurrence the
target uses), and the progression sums are accumulated by k-slices while
C(N) is formed as a single vector product.
"""
import io
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep import Lambda_bf, mu_bf  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "results", "c02_flatsum_k1.txt")
NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000,
      12_800_000, 25_600_000]
THETA = 0.56
lines = []

HEAD = [
    "STATISTIC: the flat sum of dilated walls |sum_k H(N;k)|/N computed over",
    "           two index sets -- k >= 1, which is how P2 [eq:Tw] and",
    "           {#prop:flatsum} index it, and k >= 2, which is what",
    "           lab_weight_gap.py computes -- with their difference against",
    "           C(N); the ratio |sum_k H| / |sum_k (log k) H| under both;",
    "           |T_1|/N under both; and the fitted decay exponent of each.",
    "FIELD: N = 2e5 through 2.56e7 by doubling (all of the form 2^a 5^b --",
    "       that is the target's own field and is kept here on purpose so the",
    "       comparison is like for like); theta' = THETA = 0.56, so k runs",
    "       over the squarefree k < floor(N^0.56) coprime to N and m over",
    "       1 <= m < N/k; Lambda and mu from an independent integer sieve to",
    "       2.56e7 -- mu by an omega-counter plus a squarefull mask, NOT by",
    "       the sign-flip + leftover-cofactor recurrence the target uses --",
    "       sums in float64; exponents by np.polyfit of log y on log N,",
    "       8 points.",
    "CONSTANTS: THETA = 0.56 (theta', the truncation exponent K = N^THETA);",
    "           NS = 2e5, 4e5, 8e5, 1.6e6, 3.2e6, 6.4e6, 1.28e7, 2.56e7;",
    "           sieve bound NMAX = max(NS) = 25600000; brute-force",
    "           cross-check bound 10^4.",
    "NULL: none applies. Every quantity here is a deterministic rearrangement",
    "      of one integer sieve, and the finding is a difference between two",
    "      index sets rather than a detection, so there is no sign pattern to",
    "      permute. No random number is drawn; there is no seed.",
    "DENOM: every relative figure below is divided by N, as its label states.",
    "",
]


def say(s=""):
    print(s)
    lines.append(s)


def my_sieves(n):
    """mu by omega-count + squarefull mask; Lambda by prime powers."""
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    pr = np.flatnonzero(s).astype(np.int64)
    omega = np.zeros(n + 1, dtype=np.int8)
    sqf = np.ones(n + 1, dtype=bool)
    for p in pr:
        p = int(p)
        omega[p::p] += 1
        if p * p <= n:
            sqf[p * p::p * p] = False
    mu = np.where(sqf, np.where(omega & 1, -1, 1), 0).astype(np.int8)
    mu[0] = 0
    lam = np.zeros(n + 1, dtype=np.float64)
    lg = np.log(pr.astype(np.float64))
    lam[pr] = lg
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lg[i]
            if q > n // p:
                break
            q *= p
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


say("C-02  flat sum of dilated walls with k = 1 restored")
say("=" * 78)
NMAX = max(NS)
say("independent sieve to %d ..." % NMAX)
pr, lam, mu = my_sieves(NMAX)

# brute-force validation of THIS script's own sieve
bad = [n for n in range(10_001)
       if int(mu[n]) != mu_bf(n) or abs(float(lam[n]) - Lambda_bf(n)) > 1e-12]
say("  self-check against trial division for n <= 10^4: %d mismatches"
    % len(bad))
say()

rows = []
say("  N            C(N)/N       flat_code/N   flat_paper/N   ratio_code"
    "  ratio_paper")
say("  " + "-" * 76)
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
    ks = np.array([k for k in range(1, K)
                   if mu[k] != 0 and all(k % q for q in PN)])

    f0 = np.zeros(N, dtype=np.float64)
    idx = np.arange(1, N, dtype=np.int64)
    f0[1:] = lam[1:N] * mu[N - idx]
    C = float(f0.sum())

    A = np.empty(ks.size)
    for i, k in enumerate(ks):
        k = int(k)
        r = N % k
        A[i] = f0[r::k].sum() if r else f0[k::k].sum()
    sgn = mu[ks].astype(np.float64)
    H = sgn * A
    lg = np.log(ks.astype(float))

    flat_paper = float(H.sum())            # k >= 1, as [prop:flatsum] indexes
    flat_code = float(H[ks > 1].sum())     # k >= 2, as the target computes
    logw = float((lg * H).sum())           # log 1 = 0, so k=1 is irrelevant
    iph = np.array([phi_of(int(k)) for k in ks], dtype=np.float64)
    T1 = float((sgn * (A - C / iph)).sum())
    T1_code = float((sgn[ks > 1] * (A[ks > 1] - C / iph[ks > 1])).sum())
    rows.append((N, C, flat_paper, flat_code, logw, T1, T1_code))
    say("  %-12d %+.6f    %.6f      %.6f       %.4f      %.4f"
        % (N, C / N, abs(flat_code) / N, abs(flat_paper) / N,
           abs(flat_code) / abs(logw), abs(flat_paper) / abs(logw)))

say()
say("  H(N;1) = C(N) exactly, so flat_paper - flat_code = C(N):")
for (N, C, fp, fc, lw, t1, t1c) in rows:
    say("    N = %-10d  flat_paper - flat_code = %+.6f   C(N) = %+.6f"
        % (N, fp - fc, C))
say()
say("  T_1 is unaffected by k=1 (its k=1 term is A(N;1) - C(N)/phi(1) = 0):")
for (N, C, fp, fc, lw, t1, t1c) in rows:
    say("    N = %-10d  |T_1|/N  k>=1: %.5f   k>=2: %.5f   diff %.2e"
        % (N, abs(t1) / N, abs(t1c) / N, abs(t1 - t1c) / N))

x = np.array([math.log(r[0]) for r in rows])
for label, col in (("code (k>=2)", 3), ("paper (k>=1)", 2)):
    y = np.log([abs(r[col]) / r[0] for r in rows])
    e = float(np.polyfit(x, y, 1)[0])
    say()
    say("  fitted flat-sum exponent, %-14s : %.4f" % (label, e))
ylog = np.log([abs(r[4]) / r[0] for r in rows])
say("  fitted log-weighted exponent            : %.4f"
    % float(np.polyfit(x, ylog, 1)[0]))

io.open(OUT, "w", encoding="utf-8", newline="\n").write(
    "\n".join(HEAD + lines) + "\n")
print("\nwrote", os.path.abspath(OUT))

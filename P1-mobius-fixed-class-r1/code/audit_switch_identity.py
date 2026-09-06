# -*- coding: utf-8 -*-
r"""
paper/theorem_A.md, Step 1 ([eq:switch]), Step 2 ([eq:PR], [eq:P]) and
Lemma {#lem:complete}.

WHAT IS UNDER TEST

The paper writes, for K = N^theta' and (k,N)=1,

    D(t) = sum_{k<K,(k,N)=1} mu(k) sum_{n<=t, n = N mod k} Lambda(n) mu(N-n)
         = sum_{N-t<=u<N} Lambda(N-u) mu(u) sigma_K(u)            [eq:switch]

with sigma_K(u) = sum_{k|u, k<K, (k,N)=1} mu(k), and calls the second
line a finite rearrangement, "machine-verified to 5e-10".  It then splits
sigma_K into the complete sum minus its k>=K tail to get D = P - R
([eq:PR]), where the complete sum is evaluated by

    Lemma [lem:complete]:  sum_{k|u,(k,N)=1} mu(k) = 1_{rad(u)|N},

and bounds P(t) << N^{o(1)} because mu(u) != 0 and rad(u)|N force
u | rad(N) ([eq:P]).

WHAT WOULD MAKE THIS AUDIT CIRCULAR, AND HOW IT IS AVOIDED

Building sigma_K by striding over multiples of each k IS the divisor
switch; comparing that array against the k-indexed sum would be checking
the rearrangement against itself.  So sigma_K is built twice by
genuinely different routes -- once by the stride sieve, once by
enumerating the divisors of u -- and the k-side of [eq:switch] is
accumulated by modular slices that never touch either array.

The same trap sits one step further on, and this audit sat in it.
R(t) must be built from its own definition in [eq:PR] -- the tail
sum_{k|u, k>=K, (k,N)=1} mu(k) -- by divisor enumeration.  Defining R
instead as sum g*(rad_ok - sigma_K) makes P - R == sum g*sigma_K an
array identity that holds whatever the arrays contain, and then S4 is
S1 measured a second time: it cannot fail, and [eq:PR] is never
tested at all.  It was written that way, and the shipped output said
so -- the S1 and S4 residual columns printed the same value in 24 of
25 rows.  The tail is now enumerated per u, touching neither sigma_K
nor rad_ok, and corrupting it refutes S4 while S1, S2, S3 and S5 hold.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  S1  [eq:switch] holds: the k-side and the u-side of D(t) agree to
      better than 5e-10 relative to N, at every (N,t) tested.
  S2  The stride-built sigma_K agrees exactly (integer equality) with
      divisor enumeration on every sampled u.
  S3  Lemma [lem:complete] holds exactly: for every squarefree u < N,
      the integer sum_{k|u,(k,N)=1} mu(k) equals 1 if rad(u)|N and 0
      otherwise.  Zero mismatches.
  S4  [eq:PR] holds: D(t) = P(t) - R(t) to the same 5e-10.
  S5  [eq:P] is not merely asymptotic here: the number of u<N carrying a
      nonzero term of P is exactly the number of divisors of rad(N) that
      are < N, i.e. 2^omega(N) or 2^omega(N)-1, and |P(t)| stays below
      2^omega(N) * log N at every t.

REFUTATION RULE (fixed before the run)

  S1, S4  REFUTED if max |lhs - rhs| / N > 5e-10 at any (N,t).
  S2, S3  REFUTED by a single integer mismatch.
  S5      REFUTED if the count differs from the predicted divisor count,
          or if |P(t)| exceeds 2^omega(N) log N.
  Non-zero exit on any refutation.
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
OUT = os.path.join(ROOT, "results", "audit_switch_identity.txt")

THETA = 0.56
TOL = 5e-10


def sieves(n):
    spf = np.zeros(n + 1, dtype=np.int64)
    for p in range(2, n + 1):
        if spf[p] == 0:
            blk = spf[p::p]
            spf[p::p] = np.where(blk == 0, p, blk)
    mu = np.ones(n + 1, dtype=np.int64)
    mu[0] = 0
    for v in range(2, n + 1):
        p = int(spf[v])
        w = v // p
        mu[v] = 0 if w % p == 0 else -mu[w]
    lam = np.zeros(n + 1, dtype=np.float64)
    for p in range(2, n + 1):
        if int(spf[p]) != p:
            continue
        q, lg = p, math.log(p)
        while q <= n:
            lam[q] = lg
            if q > n // p:
                break
            q *= p
    return spf, mu, lam


def prime_set(N, spf):
    s = set()
    while N > 1:
        p = int(spf[N])
        s.add(p)
        while N % p == 0:
            N //= p
    return s


def divisors_of(u, spf):
    ds = [1]
    while u > 1:
        p = int(spf[u])
        e = 0
        while u % p == 0:
            u //= p
            e += 1
        ds = [d * p ** i for d in ds for i in range(e + 1)]
    return ds


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    rng = np.random.default_rng(20260807)
    s1 = s2 = s3 = s4 = s5 = True
    worst_switch = 0.0
    worst_pr = 0.0

    say("STATISTIC and FIELD are in the header of this file.")
    say()
    for N in (25_000, 50_000, 100_000, 200_000, 400_000):
        spf, mu, lam = sieves(N)
        PN = prime_set(N, spf)
        K = int(N ** THETA)
        w = len(PN)

        # f[n] = Lambda(n) mu(N-n),  1 <= n < N
        f = np.zeros(N, dtype=np.float64)
        idx = np.arange(1, N, dtype=np.int64)
        f[1:] = lam[1:N] * mu[N - idx]

        # sigma_K by stride sieve
        sig = np.zeros(N, dtype=np.int64)
        for k in range(1, K):
            if mu[k] == 0:
                continue
            if any(k % p == 0 for p in PN):
                continue
            sig[k::k] += int(mu[k])

        # complete divisor sum, as an indicator, and by brute enumeration
        rad_ok = np.zeros(N, dtype=np.int64)
        for u in range(1, N):
            if mu[u] == 0:
                continue
            v, ok = u, True
            while v > 1:
                p = int(spf[v])
                if p not in PN:
                    ok = False
                    break
                while v % p == 0:
                    v //= p
            rad_ok[u] = 1 if ok else 0

        # S3: Lemma [lem:complete], exact, over every squarefree u < N.
        # The same enumeration yields the k >= K tail of [eq:PR] per u,
        # which is the only route by which R(t) is built below.
        tail = np.zeros(N, dtype=np.int64)
        bad3 = 0
        for u in range(1, N):
            if mu[u] == 0:
                continue
            tot = 0
            for d in divisors_of(u, spf):
                if mu[d] == 0:
                    continue
                if any(d % p == 0 for p in PN):
                    continue
                tot += int(mu[d])
                if d >= K:
                    tail[u] += int(mu[d])
            if tot != int(rad_ok[u]):
                bad3 += 1
                if bad3 <= 3:
                    say("  S3 mismatch u=%d got %d want %d"
                        % (u, tot, int(rad_ok[u])))
        if bad3:
            s3 = False

        # S2: sigma_K by divisor enumeration on a sample
        sample = rng.choice(np.arange(1, N), size=min(4000, N - 1),
                            replace=False)
        bad2 = 0
        for u in sample:
            u = int(u)
            tot = 0
            for d in divisors_of(u, spf):
                if d >= K or mu[d] == 0:
                    continue
                if any(d % p == 0 for p in PN):
                    continue
                tot += int(mu[d])
            if tot != int(sig[u]):
                bad2 += 1
        if bad2:
            s2 = False

        # S5: the support of P
        rad = 1
        for p in PN:
            rad *= p
        n_supp = int(rad_ok[1:].sum())
        want = sum(1 for d in divisors_of(rad, spf) if d < N)
        if n_supp != want:
            s5 = False

        say("N = %-8d K = %-6d omega(N) = %d  rad(N) = %-10d "
            "|supp P| = %d (want %d)" % (N, K, w, rad, n_supp, want))
        say("    S3 mismatches over all squarefree u<N : %d" % bad3)
        say("    S2 mismatches over %d sampled u        : %d"
            % (len(sample), bad2))
        say("      t          D (k-side)      D (u-side)      |diff|/N"
            "      P(t)-R(t)       |diff|/N")

        for frac in (0.125, 0.25, 0.5, 0.75, 1.0):
            t = max(2, int(frac * (N - 1)))
            # k-side: modular slices, touching neither sig nor rad_ok
            dk = 0.0
            for k in range(1, K):
                if mu[k] == 0:
                    continue
                if any(k % p == 0 for p in PN):
                    continue
                r = N % k
                if r == 0:
                    r = k
                dk += int(mu[k]) * float(f[r:t + 1:k].sum())
            # u-side
            g = f[1:t + 1]
            uu = N - np.arange(1, t + 1, dtype=np.int64)
            du = float((g * sig[uu]).sum())
            P = float((g * rad_ok[uu]).sum())
            R = float((g * tail[uu]).sum())

            e1 = abs(dk - du) / N
            e2 = abs(dk - (P - R)) / N
            worst_switch = max(worst_switch, e1)
            worst_pr = max(worst_pr, e2)
            if e1 > TOL:
                s1 = False
            if e2 > TOL:
                s4 = False
            say("      %-10d %+.8e  %+.8e  %.2e  %+.8e  %.2e"
                % (t, dk, du, e1, P - R, e2))

        cap = (2 ** w) * math.log(N)
        Pmax = max(abs(float((f[1:t2 + 1] *
                              rad_ok[N - np.arange(1, t2 + 1)]).sum()))
                   for t2 in (N // 4, N // 2, N - 1))
        say("    max_t |P(t)| = %.4f   cap 2^omega(N) log N = %.1f   %s"
            % (Pmax, cap, "ok" if Pmax <= cap else "OVER"))
        if Pmax > cap:
            s5 = False
        say()

    say("=" * 70)
    say("worst |k-side - u-side| / N     = %.2e   (tol %.0e)"
        % (worst_switch, TOL))
    say("worst |D - (P - R)| / N         = %.2e   (tol %.0e)"
        % (worst_pr, TOL))
    say("S1 %s  S2 %s  S3 %s  S4 %s  S5 %s"
        % tuple("hold" if v else "REFUTED" for v in (s1, s2, s3, s4, s5)))
    verdict = s1 and s2 and s3 and s4 and s5
    say("Step 1, [eq:PR], [eq:P] and Lemma {#lem:complete} stand"
        if verdict else "REFUTED")

    head = [
        "STATISTIC: (a) |D(t) computed k-first minus D(t) computed",
        "           u-first| / N, the switch of [eq:switch]; (b) integer",
        "           mismatches between the stride-built sigma_K(u) and",
        "           divisor enumeration; (c) integer mismatches between",
        "           sum_{k|u,(k,N)=1} mu(k) and 1_{rad(u)|N}, Lemma",
        "           {#lem:complete}; (d) |D(t) - (P(t)-R(t))| / N, [eq:PR];",
        "           (e) |supp P| against the divisor count of rad(N), and",
        "           max_t |P(t)| against 2^omega(N) log N, [eq:P].",
        "FIELD: N = 2.5e4, 5e4, 1e5, 2e5, 4e5 with theta' = 0.56 so",
        "       K = N^0.56; t = 0.125, 0.25, 0.5, 0.75, 1.0 times N-1;",
        "       (b) 4000 uniformly sampled u per N without replacement,",
        "       seed 20260807; (c) every squarefree u < N; Lambda and mu",
        "       from an integer sieve, sums in float64.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not verdict:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())

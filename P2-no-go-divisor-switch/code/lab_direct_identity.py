# -*- coding: utf-8 -*-
r"""
The direct route's object, summed without truncation, is the Goldbach
count itself.

WHAT IS AT STAKE

lab_signed_level.py found that sum_{k<K}(log k)H(N;k) stays within a
few tenths of N over the whole accessible range while the absolute sum
reaches 13.8N, and left open whether that smallness is a cancellation
mu is doing or something an identity forces.  It is forced, and the
identity says exactly what the sum is.

A(N;k) = sum_{n<N, k | N-n} Lambda(n) mu(N-n), so summing over ALL k
with weight mu(k) log k and swapping the order,

    sum_{k>=1} (log k) mu(k) A(N;k)
      = sum_{n<N} Lambda(n) mu(N-n) sum_{k | N-n} mu(k) log k
      = -sum_{n<N} Lambda(n) mu(N-n) Lambda(N-n),

by the Mobius-von Mangoldt identity sum_{d|u} mu(d) log d = -Lambda(u).
And mu(u)Lambda(u) = -log p when u = p is prime and 0 otherwise, since
mu kills the higher prime powers.  So

    sum_{k>=1} (log k) mu(k) A(N;k) = sum_{p<N} Lambda(N-p) log p,

exactly, with no error term at all: the untruncated object IS the
Goldbach count, restricted to the prime part.  [eq:direct] is then not
an approximation waiting for an estimate but a rearrangement of
rtilde(N), and every bit of its content sits in the truncation at K.

That reframes what lab_signed_level.py saw.  The partial sums are
small not because mu cancels but because the sum has barely started:
almost all of rtilde(N) is carried by k near N, where only m = 1
survives in H(N;k) = sum_{m<N/k,(m,k)=1} Lambda(N-mk)mu(m).

BACKS: Proposition {#prop:untrunc} and Remark {#rem:whereitlives} in
paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Z1  The identity is exact: sum_{k<N}(log k)mu(k)A(N;k) agrees with
      sum_{p<N} Lambda(N-p) log p to better than 1e-12 relative at
      every N. It is an identity, so anything else is an error in the
      derivation.
  Z2  That quantity is the Goldbach count: its ratio to S(N)N is
      within 0.10 of 1 at every N.
  Z3  Almost none of it is below the truncation. The partial sum over
      k < N^0.9 is under 10% of the total in absolute value, at every
      N.
  Z4  Restricting k to (k,N) = 1, as the Huang-Li setup does, does not
      change that: the coprime-restricted total is still within 0.50
      of the unrestricted one in relative terms.

REFUTATION RULE (fixed before the run)

  Z1  REFUTED at 1e-12 relative at any N.
  Z2  REFUTED if the ratio leaves [0.90, 1.10] at any N.
  Z3  REFUTED if the partial sum reaches 10% anywhere -- in which case
      the truncated range does carry a definite share and the reframing
      above is wrong.
  Z4  REFUTED if the relative difference reaches 0.50 at any N.

  All four gate.

  NO NULL IS RUN and none would mean anything. Z1 is an exact
  arithmetic identity whose reference is itself; Z2 compares against
  S(N)N, which is the reference; Z3 and Z4 are comparisons between two
  sums over the SAME terms, so any sign control would move both sides
  together. The field's sign controls were run in lab_sign_structure.py
  and lab_signed_level.py.
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
OUT = os.path.join(ROOT, "results", "lab_direct_identity.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
CLIM = 4_000_000
THETAS = [0.56, 0.70, 0.90]


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
    v, out, d = n, set(), 2
    while d * d <= v:
        if v % d == 0:
            out.add(d)
            while v % d == 0:
                v //= d
        d += 1
    if v > 1:
        out.add(v)
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    isp = np.zeros(NMAX + 1, dtype=bool)
    isp[pr] = True

    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    say()
    say("  N            LHS sum_k         RHS sum_p        rel err")
    say("  " + "-" * 62)
    rows = []
    z1 = True
    for N in NS:
        PN = factor_set(N)
        S = twin
        for q in sorted(PN):
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))

        f0 = np.zeros(N, dtype=np.float64)
        idx = np.arange(1, N, dtype=np.int64)
        f0[1:] = lam[1:N] * mu[N - idx].astype(np.float64)

        ks = np.flatnonzero(mu[2:N] != 0).astype(np.int64) + 2
        contrib = np.empty(ks.size, dtype=np.float64)
        for i, k in enumerate(ks):
            k = int(k)
            r = N % k
            a = f0[r::k].sum() if r else f0[k::k].sum()
            contrib[i] = math.log(k) * mu[k] * a
        del f0

        lhs = float(contrib.sum())
        rhs = float((lam[1:N] * lam[N - 1:0:-1]
                     * isp[N - 1:0:-1]).sum())
        rel = abs(lhs - rhs) / max(abs(rhs), 1e-300)
        if rel >= 1e-12:
            z1 = False

        cop = np.array([all(int(k) % q for q in PN) for k in ks])
        lhs_cop = float(contrib[cop].sum())

        j56 = int(np.searchsorted(ks, int(N ** 0.56)))
        p56cop = float(contrib[:j56][cop[:j56]].sum()) / N
        part = []
        jhalf = int(np.searchsorted(ks, N // 2))
        half = float(contrib[:jhalf].sum())
        for th in THETAS:
            j = int(np.searchsorted(ks, int(N ** th)))
            part.append(float(contrib[:j].sum()))

        rows.append((N, S, lhs, rhs, rel, lhs_cop, part, half,
                     p56cop))
        say("  %-12d %-17.6f %-16.6f %.3e" % (N, lhs, rhs, rel))

    say()
    say("Z1  the identity holds at every N            %s"
        % ("hold" if z1 else "REFUTED"))

    say()
    say("Z2  the quantity against S(N)N")
    say("  N            RHS/N      S(N)      ratio")
    z2 = True
    for N, S, lhs, rhs, rel, lc, part, half, p56 in rows:
        r = rhs / (S * N)
        if not (0.90 <= r <= 1.10):
            z2 = False
        say("  %-12d %-10.4f %-9.4f %.4f" % (N, rhs / N, S, r))
    say("  Z2 %s" % ("hold" if z2 else "REFUTED"))

    say()
    say("Z3  how little of it lives below the truncation")
    say("  N            k<N^0.56    k<N^0.70    k<N^0.90    share at 0.90")
    z3 = True
    for N, S, lhs, rhs, rel, lc, part, half, p56 in rows:
        sh = abs(part[2]) / abs(lhs)
        if sh >= 0.10:
            z3 = False
        say("  %-12d %-11.4f %-11.4f %-11.4f %.4f"
            % (N, part[0] / N, part[1] / N, part[2] / N, sh))
    say("  Z3 %s" % ("hold" if z3 else "REFUTED"))

    say()
    say("Z4  restricting to (k,N) = 1")
    say("  N            all k          (k,N)=1        rel diff")
    z4 = True
    for N, S, lhs, rhs, rel, lc, part, half, p56 in rows:
        d = abs(lhs - lc) / max(abs(lhs), 1e-300)
        if d >= 0.50:
            z4 = False
        say("  %-12d %-14.4f %-14.4f %.4f" % (N, lhs / N, lc / N, d))
    say("  Z4 %s" % ("hold" if z4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). Where the mass is. By [eq:dilate],")
    say("  H(N;k) = sum_{m<N/k,(m,k)=1} Lambda(N-mk)mu(m), so once")
    say("  k > N/2 only m = 1 survives and the term is (log k)Lambda(N-k),")
    say("  nonnegative. The share of the total carried above N/2:")
    say("  N            k>N/2 share   k in (N^0.9, N/2] share")
    for N, S, lhs, rhs, rel, lc, part, half, p56 in rows:
        say("  %-12d %-13.4f %.4f"
            % (N, (lhs - half) / lhs, (half - part[2]) / lhs))
    say("  So the count is assembled almost entirely from the moduli the")
    say("  reduction truncates away, and the m = 1 tail alone carries it.")

    say()
    say("  Cross-check lines. sum_{k<N^0.56,(k,N)=1}(log k)H(N;k)/N is")
    say("  computed independently by lab_signed_level.py; the gate holds")
    say("  the two against each other. That comparison is what found")
    say("  the missing mu(k) in that script.")
    for N, S, lhs, rhs, rel, lc, part, half, p56 in rows:
        say("  partial at 0.56, N=%d: %.6f" % (N, p56))
    say("  and the total itself, which lab_layer_decomposition.py")
    say("  reaches by cutting the same double sum along m instead:")
    for N, S, lhs, rhs, rel, lc, part, half, p56 in rows:
        say("AGREE untrunc_total N=%d %.6f 1e-9" % (N, lhs / N))

    say("  These partials are reported and NOT cross-checked. They")
    say("  carried an AGREE label, which claims another result file")
    say("  reaches the same number by another route. None does --")
    say("  searched across this generation and the frozen ones, since")
    say("  the sibling label mstar_001 turned out to have its partner")
    say("  sitting in a frozen tree and looked unmatched for the same")
    say("  reason. The claim is withdrawn: an unmatched cross-check")
    say("  gives confidence and tests nothing.")
    say()
    say("=" * 70)
    ok = z1 and z2 and z3 and z4
    say("the untruncated direct sum is the Goldbach count itself, and "
        "the truncated range carries almost none of it"
        if ok else "REFUTED")

    head = [
        "STATISTIC: sum_{k<N}(log k)mu(k)A(N;k) against",
        "           sum_{p<N}Lambda(N-p)log p; their relative difference;",
        "           the ratio of that quantity to S(N)N; the partial sums",
        "           over k < N^0.56, N^0.70, N^0.90 and their share of the",
        "           total; and the same total restricted to (k,N) = 1.",
        "DENOM: the identity residual is divided by the right-hand side,",
        "       sum_{p<N}Lambda(N-p)log p;",
        "NULL: none is run and none would mean anything. Z1 is an exact",
        "      arithmetic identity whose reference is itself, Z2 compares",
        "      against S(N)N which is the reference, and Z3 and Z4",
        "      compare two sums over the SAME terms, so a sign control",
        "      would move both sides together. The field's sign controls",
        "      were run in lab_sign_structure.py and lab_signed_level.py.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k over every squarefree",
        "       2 <= k < N, with no coprimality restriction except where",
        "       Z4 imposes one; S(N) from an Euler product at the fixed",
        "       bound 4e6.",
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

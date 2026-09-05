# -*- coding: utf-8 -*-
r"""
The same count cut the other way: a Mobius-signed sum of nonnegative
layers.

WHAT IS AT STAKE

Proposition {#prop:untrunc} cut sum_k (log k) H(N;k) by k and found the
count assembled as a difference of two large one-signed masses, with
the reduction's truncation falling in the middle of the cancellation.
[eq:dilate] supplies the other variable.  Writing
H(N;k) = sum_{m<N/k,(m,k)=1} Lambda(N-mk) mu(m) and exchanging,

    sum_{k} (log k) H(N;k) = sum_{m} mu(m) L(N;m),
    L(N;m) := sum_{k<N/m, (k,m)=1, mu(k)!=0} (log k) Lambda(N-mk),

and every L(N;m) is NONNEGATIVE.  So the Goldbach count is an
alternating sum, signed by mu, of nonnegative layers -- and since
L(N;m) is of size about N/m, the cancellation that produces S(N)N out
of layers of total size N log N is a Mertens cancellation in m.

That is the opposite situation from the k side.  Remark
{#rem:whycoinwins} showed a coin BEATS mu when the signs sit on k,
because there mu's progression sum is the dilated wall while a coin's
is a sum of independent signs.  Here the magnitudes are fixed and
nonnegative and only the sign pattern varies, so sum mu(m)/m -> 0
should make mu much better than a coin.  If it does not, the m side
has no Mertens gain and the decomposition is worthless; that is what
Z4 tests.

BACKS: Proposition {#prop:layers} and Remark {#rem:layerdecay} in
paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Z1  The exchange is exact: sum_m mu(m) L(N;m) reproduces
      sum_{p<N} Lambda(N-p) log p to better than 1e-12 relative.
  Z2  The first layer is a definite multiple of the answer:
      L(N;1)/(S(N)N) lies in [1, 10] at every N.
  Z3  The convergence in m is slow: the partial sum
      |sum_{m<M} mu(m) L(N;m)| at M = N^{0.25} still exceeds twice
      |total| at every N.
  Z4  mu beats a coin here. Holding every L(N;m) fixed and drawing
      eps(m) = +-1 on the squarefree m >= 2 (eps(1) = 1),
      |sum mu(m)L| is below the minimum over 16 draws of
      |sum eps(m)L| at every N.

REFUTATION RULE (fixed before the run)

  Z1  REFUTED at 1e-12 relative at any N. It is an exchange of two
      finite sums, so a failure is an error in the derivation.
  Z2  REFUTED if the ratio leaves [1, 10] at any N.
  Z3  REFUTED if the partial sum is at or under twice the total
      anywhere -- which would mean the layers converge quickly and the
      truncation in m is cheap.
  Z4  REFUTED if mu is at or above the minimum draw at any N, in which
      case the Mertens cancellation is not what is doing the work and
      this decomposition explains nothing.

  All four gate.

  THE CONTROL is Z4's sign draw, and it is the right one here for the
  reason Remark {#rem:whycoinwins} is not: the layers L(N;m) are held
  fixed, so the coin cannot buy square-root cancellation inside them.
  Only the sign pattern across m varies, which is exactly the thing
  being credited.
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
OUT = os.path.join(ROOT, "results", "lab_layer_decomposition.txt")

NS = [200_000, 400_000, 800_000, 1_600_000]
CLIM = 4_000_000
DRAWS = 16
SEED = 20260808


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
    sqf = mu != 0

    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    say()
    say("  N            sum_m mu(m)L    sum_p            rel err")
    say("  " + "-" * 62)
    rows = []
    z1 = True
    for N in NS:
        S = twin
        for q in sorted(factor_set(N)):
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))

        ms = np.flatnonzero(sqf[1:N]).astype(np.int64) + 1
        L = np.zeros(ms.size, dtype=np.float64)
        logk = np.log(np.arange(1, N, dtype=np.float64))
        for i, m in enumerate(ms):
            m = int(m)
            kmax = (N - 1) // m
            if kmax < 1:
                continue
            ok = sqf[1:kmax + 1].copy()
            for q in factor_set(m):
                ok[q - 1::q] = False
            if not ok.any():
                continue
            ks = np.flatnonzero(ok) + 1
            L[i] = float((logk[ks - 1] * lam[N - m * ks]).sum())
        tot = float((mu[ms].astype(np.float64) * L).sum())
        rhs = float((lam[1:N] * lam[N - 1:0:-1] * isp[N - 1:0:-1]).sum())
        rel = abs(tot - rhs) / max(abs(rhs), 1e-300)
        if rel >= 1e-12:
            z1 = False
        rows.append((N, S, tot, rhs, rel, ms, L))
        say("  %-12d %-16.6f %-16.6f %.3e" % (N, tot, rhs, rel))

    say()
    say("Z1  the exchange is exact at every N        %s"
        % ("hold" if z1 else "REFUTED"))

    say()
    say("Z2  the first layer against the answer")
    say("  N            L(N;1)/N   S(N)N/N   ratio")
    z2 = True
    for N, S, tot, rhs, rel, ms, L in rows:
        r = L[0] / (S * N)
        if not (1.0 <= r <= 10.0):
            z2 = False
        say("  %-12d %-10.4f %-9.4f %.4f" % (N, L[0] / N, S, r))
    say("  Z2 %s" % ("hold" if z2 else "REFUTED"))

    say()
    say("Z3  how slowly the layers converge")
    say("  N            M=N^0.25   partial/N   total/N    |partial/total|")
    z3 = True
    for N, S, tot, rhs, rel, ms, L in rows:
        M = int(N ** 0.25)
        j = int(np.searchsorted(ms, M))
        part = float((mu[ms[:j]].astype(np.float64) * L[:j]).sum())
        r = abs(part) / abs(tot)
        if r <= 2.0:
            z3 = False
        say("  %-12d %-10d %-11.4f %-10.4f %.4f"
            % (N, M, part / N, tot / N, r))
    say("  Z3 %s" % ("hold" if z3 else "REFUTED"))

    say()
    say("Z4  the control: same layers, signs on m >= 2 drawn at random")
    say("  N            |mu total|/N   draws min/N   median/N   max/N")
    z4 = True
    for j, (N, S, tot, rhs, rel, ms, L) in enumerate(rows):
        rng = np.random.default_rng(SEED + j)
        vals = []
        for d in range(DRAWS):
            e = rng.choice([-1.0, 1.0], size=ms.size)
            e[0] = 1.0
            vals.append(abs(float((e * L).sum())))
        vals = np.array(vals) / N
        if abs(tot) / N >= float(vals.min()):
            z4 = False
        say("  %-12d %-14.4f %-13.4f %-10.4f %.4f"
            % (N, abs(tot) / N, float(vals.min()),
               float(np.median(vals)), float(vals.max())))
    say("  Z4 %s" % ("hold" if z4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). The size of the layers. If")
    say("  L(N;m) ~ c(m) N/m then m L(N;m)/N is the profile in m:")
    say("  N            m=1      m=2      m=3      m=5      m=7")
    for N, S, tot, rhs, rel, ms, L in rows:
        vals = []
        for mm in (1, 2, 3, 5, 7):
            k = int(np.searchsorted(ms, mm))
            vals.append(mm * L[k] / N if k < ms.size and ms[k] == mm
                        else float("nan"))
        say("  %-12d %s" % (N, "  ".join("%-8.4f" % v for v in vals)))
    say("  The layers at m = 2 and m = 5 are empty because every N here")
    say("  is 2^a 5^b: if q | N and q | m then N - mk is divisible by q,")
    say("  so Lambda(N-mk) survives only at the single point N - mk = q.")
    say("  The share of the total mass carried by m coprime to N:")
    say("  N            (m,N)=1 mass/N   (m,N)>1 mass/N   share")
    for N, S, tot, rhs, rel, ms, L in rows:
        PN = factor_set(N)
        cop = np.ones(ms.size, dtype=bool)
        for q in sorted(PN):
            cop &= (ms % q) != 0
        a = float(L[cop].sum()) / N
        b = float(L[~cop].sum()) / N
        say("  %-12d %-16.3f %-16.3f %.6f" % (N, a, b, a / (a + b)))
    say("  which is the dual of lab_direct_identity.py's finding that")
    say("  restricting k to (k,N) = 1 changes the total by nothing.")
    say("  and the total mass, which is what has to cancel away:")
    say("  N            sum_m L(N;m)/N   |total|/N   cancellation")
    for N, S, tot, rhs, rel, ms, L in rows:
        mass = float(L.sum()) / N
        say("  %-12d %-16.1f %-11.4f %.3e"
            % (N, mass, abs(tot) / N, abs(tot) / N / mass))

    say()
    say("  Cross-check lines. The same total is reached from the k side")
    say("  by lab_direct_identity.py; the gate holds the two together.")
    for N, S, tot, rhs, rel, ms, L in rows:
        say("AGREE untrunc_total N=%d %.6f 1e-9" % (N, tot / N))

    say()
    say("=" * 70)
    ok = z1 and z2 and z3 and z4
    say("the count is a Mobius-signed sum of nonnegative layers and the "
        "cancellation in m is what produces it" if ok else "REFUTED")

    head = [
        "STATISTIC: L(N;m) = sum_{k<N/m,(k,m)=1,mu(k)!=0}(log k)Lambda(N-mk),",
        "           nonnegative; the alternating total sum_m mu(m)L(N;m)",
        "           against sum_{p<N}Lambda(N-p)log p; the first layer",
        "           against S(N)N; the partial sum at M = N^0.25; the",
        "           profile m L(N;m)/N; the total mass sum_m L(N;m); and",
        "           the same total with the signs on m >= 2 randomised.",
        "NULL: the sign draw of Z4 -- every layer L(N;m) held fixed, the",
        "      sign of each m >= 2 redrawn, 16 draws, eps(1) = 1. It is",
        "      the right control here precisely because the layers are",
        "      fixed: unlike the coin that beat mu it cannot buy",
        "      square-root cancellation inside a layer, so it isolates",
        "      the sign pattern across m and nothing else.",
        "FIELD: N = 2e5 through 1.6e6 by doubling; m over the squarefree",
        "       1 <= m < N; k over the squarefree k < N/m coprime to m;",
        "       S(N) from an Euler product at the fixed bound 4e6; seed",
        "       20260808 + j at the j-th N of the ladder, so the",
        "       four runs use 20260808 to 20260811; the declaration",
        "       named the base alone and the runs could not be",
        "       reproduced from it.",
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

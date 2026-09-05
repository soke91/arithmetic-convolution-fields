# -*- coding: utf-8 -*-
r"""
OPEN.md, "동전을 구별하는 추정기" (the rewritten pass condition) and
"표적은 theta' > 1/2 하나" -- an exact factorisation of the demand-side
discrepancy, and what it says about the coin.

WHAT IS AT STAKE

lab_level_coin_null found that a coin on the same support gives a
SMALLER B(N;K) = sum_k (log k)|E_mu(N;k)| than mu does, at every N
tested, and the level measurement built on B had to be withdrawn.  It
was left unexplained WHY mu should be noisier than random signs in
progressions.  There is an exact reason.

Write A(N;k) for the progression sum inside E_mu,

    A(N;k) = sum_{n<N, n = N (mod k)} Lambda(n) mu(N-n).

The class n = N (mod k) is k | N-n, so with u = N-n = mk,

    A(N;k) = sum_{1<=m<N/k} Lambda(N-mk) mu(mk)
           = mu(k) * H(N;k),
    H(N;k) = sum_{m<N/k, (m,k)=1} Lambda(N-mk) mu(m),

because mu(mk) = mu(m)mu(k) when (m,k)=1 and vanishes otherwise.  So
the demand-side discrepancy at modulus k is mu(k) times a DILATED
Mobius-prime correlation at scale N/k -- the supply side's own object,
one dilation down.  A coin has no such factorisation: eps(mk) is
independent of eps(m) and eps(k), so the ratio

    r_s(N;k) = A_s(N;k) / ( s(k) * H_s(N;k) )

is identically 1 for s = mu and a random ratio for s = eps.  That is a
discriminator whose gain against the coin costs no amplification at
all -- both sides of the ratio are the same size -- which is exactly
what OPEN.md's rewritten pass condition asks for.

It also predicts the sign of the K* result: |E_mu(N;k)| is governed by
|H(N;k)|, a Mobius-prime correlation of length N/k, whereas the coin's
is a sum of N/k independent signs.  Square-root cancellation is the
best a coin can do and it achieves it; mu need not, and does not.

BACKS: Proposition {#prop:dilate} and Remark {#rem:whycoinwins} in
paper/wall_v3.md, written in the same cycle as this script.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  P1  The identity A(N;k) = mu(k) H(N;k) holds to machine precision:
      the worst relative error over every squarefree k < K coprime to N
      is below 1e-12, at every N tested.
  P2  For a coin the same ratio is not 1: the median |r_eps - 1| over k
      and 8 draws exceeds 0.5.
  P3  The fraction of k with |r - 1| < 1e-9 is exactly 1 for mu and
      below 0.01 for every coin draw.
  P4  mu is noisier than the coin in progressions, which is what
      lab_level_coin_null saw: the ratio <|A_mu(N;k)|> / <|A_eps(N;k)|>
      exceeds 1 in every k-band tested.

REFUTATION RULE (fixed before the run)

  P1  REFUTED if the worst relative error reaches 1e-12 at any N.  This
      is an identity, so any failure is an error in the derivation.
  P2  REFUTED if the median is 0.5 or below.
  P3  REFUTED if mu's fraction is below 1, or if any coin draw reaches
      0.01.
  P4  REFUTED if the ratio is 1 or below in any band.

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
OUT = os.path.join(ROOT, "results", "lab_dilate_identity.txt")

NS = [200_000, 400_000, 800_000]
THETA = 0.56
DRAWS = 8
SEED = 20260808
BANDS = [(2, 32), (32, 256), (256, 2048)]


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


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    rng = np.random.default_rng(SEED)

    say()
    say("  N          K      #k      worst |A - mu(k)H| / |A|   frac exact"
        "   coin frac exact (max over %d draws)" % DRAWS)
    say("  " + "-" * 92)
    p1 = p2 = p3 = p4 = True
    all_med = []
    band_rat = {b: [] for b in BANDS}
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
        ks = [k for k in range(2, K)
              if mu[k] != 0 and all(k % q for q in PN)]

        def sums(sig):
            """A(N;k) and s(k)*H(N;k) for every k, for a sign array."""
            f = np.zeros(N, dtype=np.float64)
            idx = np.arange(1, N, dtype=np.int64)
            f[1:] = lam[1:N] * sig[N - idx]
            A = np.empty(len(ks))
            SH = np.empty(len(ks))
            for i, k in enumerate(ks):
                r = N % k
                A[i] = f[r::k].sum() if r else f[k::k].sum()
                m = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
                cop = np.gcd(m, k) == 1
                m = m[cop]
                SH[i] = sig[k] * float(
                    (lam[N - m * k] * sig[m]).sum())
            return A, SH

        muf = mu.astype(np.float64)
        A, SH = sums(muf)
        den = np.where(np.abs(A) > 0, np.abs(A), 1.0)
        relerr = np.abs(A - SH) / den
        worst = float(relerr.max())
        frac_mu = float((relerr < 1e-9).mean())
        if worst >= 1e-12:
            p1 = False
        if frac_mu < 1.0:
            p3 = False

        cf, meds, rats = [], [], []
        for t in range(DRAWS):
            sig = np.zeros(N + 1, dtype=np.float64)
            supp = np.flatnonzero(mu[:N + 1] != 0)
            sig[supp] = rng.integers(0, 2, size=supp.size) * 2.0 - 1.0
            Ae, SHe = sums(sig)
            d2 = np.where(np.abs(Ae) > 0, np.abs(Ae), 1.0)
            re = np.abs(Ae - SHe) / d2
            cf.append(float((re < 1e-9).mean()))
            meds.append(float(np.median(re)))
            rats.append(np.abs(Ae))
        if max(cf) >= 0.01:
            p3 = False
        all_med.extend(meds)
        say("  %-10d %-6d %-7d %-25.3e %-12.6f %.6f"
            % (N, K, len(ks), worst, frac_mu, max(cf)))

        kk = np.array(ks)
        eps_mean = np.mean(np.array(rats), axis=0)
        for lo, hi in BANDS:
            sel = (kk >= lo) & (kk < hi)
            if sel.sum() == 0:
                continue
            rr = float(np.abs(A[sel]).mean() / eps_mean[sel].mean())
            band_rat[(lo, hi)].append(rr)
            if rr <= 1.0:
                p4 = False

    med = float(np.median(all_med))
    p2 = med > 0.5
    say()
    say("P1  worst relative error of A = mu(k)H over all N: below 1e-12 "
        "at every N: %s   %s" % (p1, "hold" if p1 else "REFUTED"))
    say("P2  median |r_eps - 1| over k and draws = %.4f   (floor 0.5)"
        "   %s" % (med, "hold" if p2 else "REFUTED"))
    say("P3  exact-fraction 1 for mu and < 0.01 for every coin draw: %s"
        "   %s" % (p3, "hold" if p3 else "REFUTED"))
    say()
    say("P4  <|A_mu|> / <|A_eps|> by k-band")
    say("  band            %s" % "   ".join("N=%d" % N for N in NS))
    for lo, hi in BANDS:
        vals = band_rat[(lo, hi)]
        say("  [%-5d,%-6d) %s" % (lo, hi,
                                  "   ".join("%8.4f" % v for v in vals)))
    say("  P4 %s" % ("hold" if p4 else "REFUTED"))

    say()
    say("=" * 70)
    ok = p1 and p2 and p3 and p4
    say("P1 %s  P2 %s  P3 %s  P4 %s"
        % tuple("hold" if v else "REFUTED" for v in (p1, p2, p3, p4)))
    say("the demand-side discrepancy factorises as mu(k) times a dilated "
        "Mobius-prime sum, and that is why the coin wins on B"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the relative error of the identity A(N;k) =",
        "           mu(k) H(N;k), where A is the progression sum",
        "           sum_{n = N mod k} Lambda(n) mu(N-n) and H(N;k) =",
        "           sum_{m<N/k,(m,k)=1} Lambda(N-mk) mu(m); the same ratio",
        "           for coin signs on the support of mu^2; the fraction of",
        "           k at which it is exact; and the ratio of mean |A| for",
        "           mu against mean |A| for the coin, by k-band.",
        "DENOM: the identity residual is divided by |A(N;k)|, or by 1 where",
        "       A vanishes;",
        "NULL: the coin is the control and is the point -- eps(mk) is",
        "      independent of eps(m) and eps(k), so the identity is",
        "      available to mu and to nothing else on the same support.",
        "      Eight draws per N, same support, same k-range.",
        "FIELD: N = 2e5, 4e5, 8e5 with theta' = 0.56, so k runs over the",
        "       squarefree k < N^0.56 coprime to N; m over 1 <= m < N/k",
        "       with (m,k) = 1; Lambda and mu from an integer sieve to",
        "       8e5; k-bands [2,32), [32,256), [256,2048); numpy",
        "       default_rng seed 20260808.",
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

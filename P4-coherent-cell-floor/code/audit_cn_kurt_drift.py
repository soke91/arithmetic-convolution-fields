# -*- coding: utf-8 -*-
r"""
Does the one statistic that survived the coin survive growing N?

WHAT IS AT STAKE

rem:cnlaw found the first thing in this branch to survive lem:coin:
over even N in (2e6, 4e6] the excess kurtosis of
G(N) = C(N)/sqrt(V(N)) is +0.26422, while sixty-four independent sign
patterns on mu's support give -0.01825 +- 0.03957 and never exceed
+0.042.  G is not Gaussian and its non-Gaussianity is not something a
coin makes.

That was one band, and OPEN.md's item 1 says so in the same breath:
**a law needs the drift in N and one band has none.**  This is that
measurement.  It decides between two readings that the single band
cannot separate:

  a finite-N effect -- the excess kurtosis falls as the number of
  terms grows, on its way to the Gaussian the draft's law assumed; or

  a law -- it does not fall, and G is non-Gaussian at every scale
  reached, which is a property of C(N) and not of the range.

THE DESIGN

Seven octaves, (2^b, 2^(b+1)] for b = 17..23, so the top band is
N up to 1.68e7 and the whole ladder is covered by ONE pair of FFTs:
C is the convolution of Lambda with mu and V that of mu^2 with
Lambda^2, both computed once to the top and then sliced per band.
The coin arm is the same convolution with mu replaced by
eps * mu^2, eps uniform in {-1,+1}; its rfft of Lambda is cached, so a
draw costs one transform pair.  Thirty-two draws per ladder give the
error bar that blocking inside a band cannot -- rem:cnlaw's Q2
measured that blocking to be too small by a factor of 55.

BACKS: Remark {#rem:cnkurt} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  D1  THE GATE.  Re-measured on rem:cnlaw's own band, even N in
      (2e6, 4e6], the excess kurtosis reproduces its POINT
      cnlaw_kurt to three decimals.  Two independent runs of the same
      quantity, one of them at a different sieve top.
  D2  The separation is not local to that band: at every one of the
      seven octaves the real excess kurtosis lies outside the coin
      ensemble at |z| > 3.
  D3  **It is a finite-N effect**: regressing the excess kurtosis on
      log2 of the band's midpoint over the seven octaves gives a
      negative slope resolved at |t| > 2.
  D4  The control behaves: the coin's mean excess kurtosis is within
      two standard errors of that mean of zero at every band.
      [CORRECTED WORDING, 2026-09-04.  This clause and the rule below
      read "draw-to-draw standard deviations" and "two draw-sd"; the
      code has always applied the standard error of the mean, which is
      sqrt(DRAWS) = 5.66 times stricter.  It is the REGISTRATION that
      is corrected here and not the code -- the code is unchanged and
      its verdict, REFUTED, is unchanged.  D4 asks whether the coin
      arm's MEAN is distinguishable from zero, and the standard error
      of that mean is what answers it; the sd of a single draw answers
      a different question.  Loosening the code to the written unit
      would read 0.78 sd and lift the refutation, which is the one
      direction a pre-registration exists to forbid.  The result file
      prints both numbers and says which unit each is in.]

REFUTATION RULE (fixed before the run)

  D1  REFUTED outside three decimals.  Then the two runs disagree
      about the same number and neither is reported until that is
      resolved.  THIS ONE GATES.
  D2  REFUTED if any band fails |z| > 3.  The small bands are where
      to expect it -- the (2^17, 2^18] band holds 65536 values of N
      against the top band's 4.19 million -- and a failure there is
      about sample size, which is why the band's count is printed
      beside every z.
  D3  **REFUTED if the slope is positive or unresolved, and that is
      the more interesting outcome.**  A flat kurtosis across seven
      octaves would say the non-Gaussianity is a property of C and
      not of the range, which is what a law would need and is
      stronger than what is predicted here.  Unresolved is the
      genuine third case: seven points is few, the octaves are
      nested in one sieve and not independent samples, and a slope
      this measurement cannot see is not a slope shown to be absent.
  D4  REFUTED if any band's coin mean is two standard errors of that
      mean from zero (wording corrected 2026-09-04; see D4 above).
      Then the coin arm has structure of its own and the z-scores of
      D2 are measured against a moving target.

  WHAT THIS CANNOT DO.  Seven octaves of a quantity that may decay
  like an inverse power of the number of terms cannot distinguish
  decay to zero from decay to a positive limit.  If D3 holds, the
  question that replaces it is the limit, and nothing here estimates
  one.
"""

import io
import math
import os
import re
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "audit_cn_kurt_drift.txt")
SRC = os.path.join(ROOT, "results", "audit_cn_law.txt")

BLO, BHI = 17, 23                 # octaves (2^b, 2^(b+1)]
NMAX = 1 << (BHI + 1)
DRAWS = 32
SEED = 20260821
GATELO, GATEHI = 2_000_000, 4_000_000
DEC = 3
ZKEEP = 3.0


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
    return lam, mu


def excess_kurt(x):
    d = x - x.mean()
    v = float((d ** 2).mean())
    return float((d ** 4).mean() / v ** 2 - 3.0)


def ols(x, y):
    A = np.column_stack([np.ones(len(y)), x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    s2 = float((r ** 2).sum()) / (len(y) - 2)
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, np.sqrt(np.diag(cov))


def published_kurt():
    src = io.open(SRC, encoding="utf-8").read()
    m = re.search(r"^POINT cnlaw_kurt ([-\d.]+)\s*$", src, re.M)
    if not m:
        raise SystemExit("no POINT cnlaw_kurt in %s" % SRC)
    return float(m.group(1))


HEAD = [
    "STATISTIC: the excess kurtosis of G(N) = C(N)/sqrt(V(N)) over",
    "           even N in each of seven octaves, against an ensemble",
    "           of %d sign patterns on mu's support recomputed over" % DRAWS,
    "           the whole ladder, and the slope of that excess",
    "           kurtosis against log2 of the band midpoint.",
    "FIELD: even N in (2^b, 2^(b+1)] for b = %d..%d, so N runs to"
    % (BLO, BHI),
    "       %d, with Lambda and mu sieved once to that top and" % NMAX,
    "       both convolutions taken by FFT over the whole range, then",
    "       sliced per band. The gate re-measures rem:cnlaw's own",
    "       band, even N in (%d, %d]." % (GATELO, GATEHI),
    "SEED: the coin draws from numpy default_rng at seed %d; without"
    % SEED,
    "      it the file does not reproduce its own control.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = published_kurt()
    say("READ audit_cn_law.txt POINT cnlaw_kurt %.5f" % pub)
    say("  the band this ladder has to reproduce, read from that file")
    say("PRINTBOUND audit_cn_kurt_drift %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))

    say("sieving to %d" % NMAX)
    lam, mu = sieves(NMAX)
    mu2 = (mu.astype(np.float64)) ** 2
    L = 1 << (2 * NMAX).bit_length()
    say("  transform length %d" % L)
    FL = np.fft.rfft(lam, L)
    C = np.fft.irfft(FL * np.fft.rfft(mu.astype(np.float64), L),
                     L)[: NMAX + 1]
    V = np.fft.irfft(np.fft.rfft(lam ** 2, L) * np.fft.rfft(mu2, L),
                     L)[: NMAX + 1]
    del mu
    rootV = np.sqrt(np.maximum(V, 1e-300))
    del V

    bands = []
    for b in range(BLO, BHI + 1):
        lo, hi = 1 << b, 1 << (b + 1)
        Ns = np.arange(lo + 2, hi + 1, 2, dtype=np.int64)
        bands.append((b, lo, hi, Ns))
    gate_Ns = np.arange(GATELO + 2, GATEHI + 1, 2, dtype=np.int64)

    real = [excess_kurt(C[Ns] / rootV[Ns]) for _, _, _, Ns in bands]
    gate_k = excess_kurt(C[gate_Ns] / rootV[gate_Ns])
    del C

    # -------------------------------------------------------------- D1
    say()
    say("D1  does this run reproduce the band rem:cnlaw measured?")
    say("  here %.5f against its %.5f" % (gate_k, pub))
    d1 = abs(round(gate_k, DEC) - round(pub, DEC)) < 10.0 ** (-DEC) / 2
    say("  D1 %s   (cap: %d decimals)"
        % ("hold" if d1 else "REFUTED", DEC))
    if not d1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    say()
    say("  drawing the coin %d times over the whole ladder" % DRAWS)
    rng = np.random.default_rng(SEED)
    coin = [[] for _ in bands]
    for d in range(DRAWS):
        eps = rng.integers(0, 2, size=NMAX + 1).astype(np.float64) * 2 - 1
        Cc = np.fft.irfft(FL * np.fft.rfft(eps * mu2, L), L)[: NMAX + 1]
        for i, (_, _, _, Ns) in enumerate(bands):
            coin[i].append(excess_kurt(Cc[Ns] / rootV[Ns]))
        del Cc
    coin = [np.array(c) for c in coin]

    say()
    say("      b        band top     even N      real      coin mean"
        "     coin sd       z")
    zs = []
    for i, (b, lo, hi, Ns) in enumerate(bands):
        cm, cs = float(coin[i].mean()), float(coin[i].std(ddof=1))
        z = (real[i] - cm) / cs
        zs.append(z)
        say("     %2d  %14d  %9d  %+8.5f  %+11.5f  %10.5f  %6.2f"
            % (b, hi, len(Ns), real[i], cm, cs, z))
        say("POINT cnkurt_%d %.5f" % (b, real[i]))
    say("SCALES %d" % len(bands))

    # -------------------------------------------------------------- D2
    say()
    say("D2  is the separation there at every octave?")
    d2 = all(abs(z) > ZKEEP for z in zs)
    say("  smallest |z| across the seven bands is %.2f" % min(map(abs, zs)))
    say("TSTAT cnkurt_min_z %.2f" % min(zs, key=abs))
    say("SPREAD cnkurt_min_z %.5f"
        % float(coin[int(np.argmin(np.abs(zs)))].std(ddof=1)))
    say("  D2 %s   (cap: |z| > %.0f at every band)"
        % ("hold" if d2 else "REFUTED", ZKEEP))

    # -------------------------------------------------------------- D3
    say()
    say("D3  does it drift with N?")
    x = np.array([b + 0.5 for b, _, _, _ in bands])
    y = np.array(real)
    c, se = ols(x, y)
    t = c[1] / se[1]
    say("  excess kurtosis on log2(band midpoint): slope %+.5f +- "
        "%.5f, t = %.2f" % (c[1], se[1], t))
    say("TSTAT cnkurt_slope %.2f" % t)
    say("SPREAD cnkurt_slope %.5f" % (x.max() - x.min()))
    if abs(t) < 2.0:
        say("UNRESOLVED SIGN cnkurt_slope")
    d3 = c[1] < 0 and abs(t) > 2.0
    say("  D3 %s   (cap: negative and |t| > 2)"
        % ("hold" if d3 else "REFUTED"))
    say("  over the seven octaves it runs %+.5f to %+.5f"
        % (y[0], y[-1]))
    say("  a power law in N would show as a straight line in "
        "log(kurtosis)")
    pos = y > 0
    if pos.all():
        cp, sep = ols(x, np.log(y))
        say("  log(kurtosis) on log2 N: slope %+.5f +- %.5f, so "
            "kurtosis ~ N^%+.4f" % (cp[1], sep[1], cp[1] / math.log(2)))
        say("POINT cnkurt_power %.5f" % (cp[1] / math.log(2)))
    else:
        say("  not every band is positive, so no power law is fitted")

    # -------------------------------------------------------------- D4
    say()
    say("D4  does the control sit at zero?")
    worst, wb = 0.0, None
    for i, (b, _, _, _) in enumerate(bands):
        cm, cs = float(coin[i].mean()), float(coin[i].std(ddof=1))
        r = abs(cm) / cs * math.sqrt(DRAWS)
        if r > worst:
            worst, wb = r, b
    say("  the coin's mean is furthest from zero at b = %s, at %.2f "
        "of its own" % (wb, worst))
    say("  standard error over %d draws" % DRAWS)
    d4 = worst <= 2.0
    say("  D4 %s   (cap: two standard errors)"
        % ("hold" if d4 else "REFUTED"))
    say("  Note on the unit. As first written the pre-registration put")
    say("  the cap at two of the coin's own DRAW-TO-DRAW STANDARD")
    say("  DEVIATIONS; the code applies two STANDARD ERRORS of the mean")
    say("  over %d draws, which" % DRAWS)
    say("  is sqrt(%d) = %.2f times stricter. In the registered unit the"
        % (DRAWS, math.sqrt(DRAWS)))
    say("  worst band reads %.2f sd and D4 would hold." % (worst / math.sqrt(DRAWS)))
    say("  The code is the stricter of the two and its verdict stands:")
    say("  the question D4 asks is whether the coin arm's MEAN is")
    say("  distinguishable from zero, and the standard error of that")
    say("  mean is what answers it. The registration named the wrong")
    say("  unit for its own question. Loosening it now would lift a")
    say("  refutation after the fact, which is the one direction a")
    say("  pre-registration exists to forbid.")
    say("  On 2026-09-04 the REGISTRATION WORDING was corrected to say")
    say("  standard errors, with the old wording kept beside it in the")
    say("  docstring. The code was not touched and D4's verdict did not")
    say("  move. What changed is a unit written wrongly, not a cap")
    say("  chosen after seeing the number.")

    say()
    say("  NOT PRE-REGISTERED, forced by D4 breaking and by the")
    say("  exponent itself. No cap above is changed.")
    cmean = np.array([float(c.mean()) for c in coin])
    adj = y - cmean
    ca, sea = ols(x, np.log(adj))
    say("  (a) the control sits at %+.5f to %+.5f rather than zero, "
        "so the" % (cmean.min(), cmean.max()))
    say("      power law is refitted on real minus the control's own "
        "mean:")
    say("      kurtosis ~ N^%+.4f against the raw fit's N^%+.4f"
        % (ca[1] / math.log(2), cp[1] / math.log(2)))
    say("POINT cnkurt_power_adj %.5f" % (ca[1] / math.log(2)))
    say("  (b) a sum of n independent terms has excess kurtosis of "
        "order 1/n,")
    say("      and this field has about N/log N terms, so independence "
        "would")
    say("      give N^-1 up to logs. The measured exponent is not "
        "that:")
    zind = abs(cp[1] / math.log(2) + 1.0) / (sep[1] / math.log(2))
    say("      %+.4f against -1, a nominal %.1f standard errors away."
        % (cp[1] / math.log(2), zind))
    say("TSTAT cnkurt_vs_independent %.2f" % zind)
    say("SPREAD cnkurt_vs_independent %.5f" % (sep[1] / math.log(2)))
    say("      Nominal: the seven octaves come from one sieve and are "
        "not")
    say("      independent samples, so the standard error is a lower "
        "bound")
    say("      on the true one and no significance is claimed from it. "
        "What")
    say("      the figure says is that the decay is visibly slower "
        "than the")
    say("      independent-sum rate over the range measured.")

    say()
    say("=" * 70)
    say("D1 %s  D2 %s  D3 %s  D4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (d1, d2, d3, d4)))
    say()
    if d2 and d3:
        say("the non-Gaussianity is real at every octave and it "
            "shrinks. So it is")
        say("a finite-N effect, not yet a law: what rem:cnlaw found "
            "survives the")
        say("coin at every scale reached and is on its way somewhere "
            "this")
        say("measurement cannot name -- seven octaves cannot tell "
            "decay to zero")
        say("from decay to a positive limit, and the limit is the "
            "question that")
        say("replaces this one.")
    elif d2 and not d3:
        say("the non-Gaussianity is real at every octave and does not "
            "shrink. That")
        say("is what a law needs: a property of C(N) rather than of "
            "the range,")
        say("and the strongest outcome this design could have "
            "produced.")
    else:
        say("the separation does not hold at every octave, so what "
            "rem:cnlaw")
        say("measured is a property of its band. The bands where it "
            "fails are")
        say("printed above with their counts, which is where to look "
            "first.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())

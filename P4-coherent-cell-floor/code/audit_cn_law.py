# -*- coding: utf-8 -*-
r"""
The withdrawn law for C(N), re-measured with its statistics defined

WHAT IS AT STAKE

OPEN.md's wall item 1 is the first on the list and has never moved.
The draft's law

    C(N) = m(N) + sqrt(V(N)) * G(N)

was withdrawn, and the reasons were about evidence rather than truth:
the bulk and the tail did not reproduce under the cell index and
pooling the text specified, the phase content was reproduced by the
coin and so fell to lem:coin, and the mask's significance was
overstated.  The closing condition written there is exactly one thing
-- **define those statistics and measure them again**.

This is the scalar the whole problem reduces to.  Section sec:wall's
own chain writes binary Goldbach's demand side as
C(N) = sum over n < N of Lambda(n) mu(N-n) = o(N), so a law for C(N)
is not a side question about a fluctuation; it is the object.

WHAT MAKES THE RE-MEASUREMENT POSSIBLE WITHOUT FITTING

Proposition prop:V gives the scale exactly:

    V(N) = sum over v < N of mu^2(v) Lambda(N-v)^2,

a quantity computable rather than fitted.  So G(N) is *defined*, not
estimated, once m(N) is:

    G(N) = (C(N) - m(N)) / sqrt(V(N)).

That removes the freedom the withdrawn version had.  Nothing below
chooses a cell index or a pooling; the statistics are taken over every
even N in one band.

BOTH ARMS ARE COMPUTED THE SAME WAY.  C is the linear convolution of
Lambda with mu, and V the linear convolution of mu^2 with Lambda^2, so
one pair of FFTs gives every N in the band at once.  The coin arm
replaces mu(v) by eps(v) mu^2(v) with eps uniform in {-1,+1}: it keeps
the support and therefore keeps V exactly, and destroys only the sign
pattern.  That is the control lem:coin is about, and the withdrawal
happened because it had not been run alongside.

BACKS: Remark {#rem:cnlaw} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  P1  THE GATE.  The convolution is right: C(N) computed by FFT agrees
      with a direct sum over n < N at eight sampled N, to 1e-6
      absolute.
  P2  **The scale needs no fitting**: the standard deviation of
      G = C/sqrt(V) over the band lies in [0.95, 1.05].  prop:V says
      V is the second moment, so a measured sd near 1 is that
      proposition's own prediction and not a fit.
  P3  G is not Gaussian: its excess kurtosis is positive and resolved
      against a blocked standard error.
  P4  **The coin reproduces both**: the coin arm's sd and excess
      kurtosis each sit within two blocked standard errors of the real
      arm's, so neither statistic carries arithmetic content.  This is
      what sank the phase claim, asked of the bulk.
  P5  The mean field has no measured content at this range: the mean
      of G, against a standard error blocked over 64 contiguous blocks
      of N, is not resolved.

REFUTATION RULE (fixed before the run)

  P1  REFUTED at any sampled N outside 1e-6.  THIS ONE GATES: if the
      convolution is wrong nothing below means anything.
  P2  REFUTED outside [0.95, 1.05].  Above the band would say V
      understates the second moment, below that it overstates; either
      way prop:V would be describing something other than what is
      measured here, and that is a claim about prop:V, not about the
      law.
  P3  REFUTED if the excess kurtosis is negative or unresolved.
      **Unresolved is a real possibility** and would mean only that
      one band cannot see it, not that G is Gaussian -- the sample is
      one band of even N and neighbouring N are not independent,
      which is why the standard error is blocked.
  P4  REFUTED if either statistic differs by more than two blocked
      standard errors.  A difference is the interesting outcome:
      it would mean the bulk carries something the coin cannot make,
      which is what the withdrawn law needed and did not have.  No
      difference means the bulk is as empty as the phase was.
  P5  REFUTED if the mean is resolved.  Then there is a mean field at
      this range and m(N) is not zero, which is a finding about the
      first term of the withdrawn law.

  The blocking is fixed at 64 contiguous blocks before the run and is
  not tuned afterwards.  It is the honest standard error here because
  C(N) at neighbouring N shares almost all of its summands; an
  unblocked error would be too small by a factor nobody has measured.

  WHAT THIS CANNOT DO.  One band cannot give a law in N.  Every
  statement below is at one range and the drift of any of these
  quantities with N is not measured here.
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
OUT = os.path.join(ROOT, "results", "audit_cn_law.txt")

NMAX = 4_000_000
NLO = 2_000_000          # the band is even N in (NLO, NMAX]
BLOCKS = 64
SEED = 20260820
SAMPLES = 8
TOLC = 1e-6
SDLO, SDHI = 0.95, 1.05


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


def conv(a, b, nmax):
    """linear convolution, index n gives sum over i+j=n"""
    L = 1 << (int(nmax + nmax).bit_length())
    fa = np.fft.rfft(a, L)
    fb = np.fft.rfft(b, L)
    return np.fft.irfft(fa * fb, L)[: nmax + 1]


def blocked(x, k):
    """mean, and its standard error over k contiguous blocks"""
    m = len(x) // k * k
    b = x[:m].reshape(k, -1).mean(axis=1)
    return float(x.mean()), float(b.std(ddof=1) / math.sqrt(k))


def blocked_stat(x, k, f):
    """f over the whole sample, and its error over k blocks"""
    m = len(x) // k * k
    b = np.array([f(s) for s in x[:m].reshape(k, -1)])
    return float(f(x)), float(b.std(ddof=1) / math.sqrt(k))


def excess_kurt(x):
    d = x - x.mean()
    v = float((d ** 2).mean())
    return float((d ** 4).mean() / v ** 2 - 3.0)


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    say("sieving to %d" % NMAX)
    lam, mu = sieves(NMAX)
    muf = mu.astype(np.float64)
    say("  Lambda nonzero at %d indices, mu nonzero at %d"
        % (int((lam != 0).sum()), int((mu != 0).sum())))
    say("SEED: the coin draws from numpy default_rng at seed %d; "
        "without it" % SEED)
    say("      the file does not reproduce its own control")

    C = conv(lam, muf, NMAX)
    V = conv(muf ** 2, lam ** 2, NMAX)

    # -------------------------------------------------------------- P1
    say()
    say("P1  is the convolution right?")
    rng = np.random.default_rng(SEED)
    picks = np.sort(rng.integers(NLO, NMAX, size=SAMPLES) // 2 * 2)
    worst = 0.0
    say("  %d sampled N, drawn at the declared seed" % SAMPLES)
    say("            N       direct           fft        |diff|")
    for N in picks:
        N = int(N)
        direct = float(np.dot(lam[1:N], muf[N - 1:0:-1]))
        d = abs(direct - C[N])
        worst = max(worst, d)
        say("     %9d  %12.6f  %12.6f  %12.3e" % (N, direct, C[N], d))
    p1 = worst <= TOLC
    say("  worst %.3e against the cap %.0e" % (worst, TOLC))
    say("  P1 %s   (cap: %.0e absolute)"
        % ("hold" if p1 else "REFUTED", TOLC))
    if not p1:
        say()
        say("  P1 gates. The convolution is wrong, so nothing below "
            "is reported.")
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    Ns = np.arange(NLO + 2, NMAX + 1, 2, dtype=np.int64)
    Cn, Vn = C[Ns], V[Ns]
    G = Cn / np.sqrt(Vn)
    say()
    say("  the band is even N in (%d, %d], %d values"
        % (NLO, NMAX, len(Ns)))
    say("  V/(N log N) runs %.5f to %.5f across the band"
        % ((Vn / (Ns * np.log(Ns))).min(),
           (Vn / (Ns * np.log(Ns))).max()))

    eps = rng.integers(0, 2, size=NMAX + 1).astype(np.float64) * 2 - 1
    coin = conv(lam, eps * muf ** 2, NMAX)
    Gc = coin[Ns] / np.sqrt(Vn)
    say("  the coin keeps mu^2 and so keeps V exactly; only the sign "
        "pattern")
    say("  is replaced")

    # -------------------------------------------------------------- P2
    say()
    say("P2  does the scale need fitting?")
    sd, sd_se = blocked_stat(G, BLOCKS, lambda s: float(s.std(ddof=1)))
    say("  sd(G) = %.5f +- %.5f over %d blocks" % (sd, sd_se, BLOCKS))
    say("POINT cnlaw_sd %.5f" % sd)
    p2 = SDLO <= sd <= SDHI
    say("  P2 %s   (cap: [%.2f, %.2f])"
        % ("hold" if p2 else "REFUTED", SDLO, SDHI))
    say("SPREAD cnlaw_sd %.5f" % sd_se)

    # -------------------------------------------------------------- P3
    say()
    say("P3  is G Gaussian in the bulk?")
    k, k_se = blocked_stat(G, BLOCKS, excess_kurt)
    t = k / k_se if k_se else float("inf")
    say("  excess kurtosis %.5f +- %.5f, t = %.2f" % (k, k_se, t))
    say("POINT cnlaw_kurt %.5f" % k)
    say("TSTAT cnlaw_kurtosis %.2f" % t)
    say("SPREAD cnlaw_kurtosis %.5f" % k_se)
    if abs(t) < 2.0:
        say("UNRESOLVED SIGN cnlaw_kurtosis")
    p3 = k > 0 and abs(t) > 2.0
    say("  P3 %s   (cap: positive and |t| > 2)"
        % ("hold" if p3 else "REFUTED"))

    # -------------------------------------------------------------- P4
    say()
    say("P4  can the coin make both?")
    sdc, sdc_se = blocked_stat(Gc, BLOCKS,
                               lambda s: float(s.std(ddof=1)))
    kc, kc_se = blocked_stat(Gc, BLOCKS, excess_kurt)
    dsd = abs(sd - sdc) / math.sqrt(sd_se ** 2 + sdc_se ** 2)
    dk = abs(k - kc) / math.sqrt(k_se ** 2 + kc_se ** 2)
    say("            statistic      real          coin      |gap|/se")
    say("     %16s  %12.5f  %12.5f  %8.2f" % ("sd", sd, sdc, dsd))
    say("     %16s  %12.5f  %12.5f  %8.2f"
        % ("excess kurtosis", k, kc, dk))
    say("TSTAT cnlaw_sd_gap %.2f" % dsd)
    say("SPREAD cnlaw_sd_gap %.5f" % math.sqrt(sd_se ** 2 + sdc_se ** 2))
    say("TSTAT cnlaw_kurt_gap %.2f" % dk)
    say("SPREAD cnlaw_kurt_gap %.5f"
        % math.sqrt(k_se ** 2 + kc_se ** 2))
    if dsd < 2.0:
        say("UNRESOLVED SIGN cnlaw_sd_gap")
    if dk < 2.0:
        say("UNRESOLVED SIGN cnlaw_kurt_gap")
    p4 = dsd <= 2.0 and dk <= 2.0
    say("  P4 %s   (cap: both within two blocked standard errors)"
        % ("hold" if p4 else "REFUTED"))

    # -------------------------------------------------------------- P5
    say()
    say("P5  is there a mean field?")
    mn, mn_se = blocked(G, BLOCKS)
    tm = mn / mn_se if mn_se else float("inf")
    say("  mean(G) = %+.5f +- %.5f, t = %.2f" % (mn, mn_se, tm))
    say("TSTAT cnlaw_mean %.2f" % tm)
    say("SPREAD cnlaw_mean %.5f" % mn_se)
    if abs(tm) < 2.0:
        say("UNRESOLVED SIGN cnlaw_mean")
    mnc, mnc_se = blocked(Gc, BLOCKS)
    say("  the coin's mean is %+.5f +- %.5f for comparison"
        % (mnc, mnc_se))
    p5 = abs(tm) < 2.0
    say("  P5 %s   (cap: |t| < 2)" % ("hold" if p5 else "REFUTED"))

    say()
    say("  the tail, since the withdrawn law named it")
    say("       quantile      real        coin      normal")
    for q in (0.001, 0.01, 0.99, 0.999):
        say("       %8.3f  %9.4f  %9.4f  %9.4f"
            % (q, float(np.quantile(G, q)), float(np.quantile(Gc, q)),
               float(math.sqrt(2) * _erfinv(2 * q - 1))))
    say("SCALES 1")

    say()
    say("=" * 70)
    say("P1 %s  P2 %s  P3 %s  P4 %s  P5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (p1, p2, p3, p4, p5)))
    say()
    if p2:
        say("the scale is not fitted: prop:V's V is the second moment "
            "of C on")
        say("this band to within the blocked error, so the law's "
            "sqrt(V) factor")
        say("is the one part of it that was never in doubt.")
    if p4:
        say("and the bulk goes the way the phase went. sd and excess "
            "kurtosis are")
        say("both reproduced by a field with mu's support and none of "
            "mu's signs,")
        say("so neither is evidence for an arithmetic law -- the same "
            "reading")
        say("lem:coin forced on the phase content, now for the bulk.")
    elif not p4:
        say("and the bulk does NOT go the way the phase went: the "
            "coin misses at")
        say("least one of the two statistics, which is content the "
            "withdrawn law")
        say("needed and did not have. What separates them is named "
            "above.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


def _erfinv(y):
    """inverse error function, series plus two Newton steps"""
    a = 0.147
    ln = math.log(1 - y * y) if abs(y) < 1 else -700.0
    t = 2 / (math.pi * a) + ln / 2
    x = math.copysign(math.sqrt(max(math.sqrt(t * t - ln / a) - t, 0.0)),
                      y)
    for _ in range(2):
        err = math.erf(x) - y
        x -= err / (2 / math.sqrt(math.pi) * math.exp(-x * x))
    return x


HEAD = [
    "STATISTIC: G(N) = C(N)/sqrt(V(N)) for even N in one band, with",
    "           C the convolution of Lambda with mu and V the exact",
    "           second moment of prop:V, so G is defined and not",
    "           fitted; its standard deviation, excess kurtosis, mean",
    "           and tail quantiles, each against a coin arm that",
    "           keeps mu's support and replaces its signs.",
    "FIELD: even N in (%d, %d], %d values, with Lambda"
    % (NLO, NMAX, (NMAX - NLO) // 2),
    "       and mu sieved to %d and both convolutions taken by FFT"
    % NMAX,
    "       over the whole range at once. Standard errors are blocked",
    "       over %d contiguous blocks of N, fixed before the run,"
    % BLOCKS,
    "       because neighbouring N share almost all their summands.",
    "SEED: the coin draws from numpy default_rng at seed %d; without"
    % SEED,
    "      it the file does not reproduce its own control.",
    "",
]


if __name__ == "__main__":
    sys.exit(main())

# -*- coding: utf-8 -*-
r"""
The last axis: does C(N) carry class structure in N's own arithmetic?

WHAT IS AT STAKE

Five routes in this branch are closed, and every one of them used N
only as a band.  rem:cnkurtlimit and rem:cnskew took marginal moments,
rem:cnshift the additive two-point function, rem:cndilation the
multiplicative one, rem:cnwherereach the decomposition by where the
prime sits.  All five are reproduced by independent signs on mu's
support, and rem:cndilation drew the conclusion the pre-registration
had fixed: the coin's agreement with C(N) at these ranges is not an
artefact of which statistics were chosen.

**One axis is untouched and it is the only one the coin cannot
follow.**  The coin's eps does not know N.  Its field does -- the
weights Lambda(N-v) carry N's arithmetic, so a coin draw has class
structure of its own -- but nothing in the coin ties mu's *signs* to
which small primes divide N.  Conjecture conj:L says exactly that:
field = M x G with M the deterministic local mask and G carrying no
class structure.  Testing it on C(N) itself is testing it on the
scalar the whole problem reduces to.

THE DESIGN

Partition the even N of one band by which of 3, 5, 7 divide them --
eight classes, the smallest holding about one N in a hundred and five.
Within each class measure the first and second moments of
G(N) = C(N)/sqrt(V(N)).  The second moment is the sharper of the two:
averaged over signs the coin gives E[G^2] = 1 in *every* class
exactly, because V(N) is the exact second moment at each N, so the
null is not "no difference between classes" but "no structure at all",
and any class structure in the real arm is a departure from something
computed rather than fitted.

A single sign pattern still fluctuates, so both arms are scored the
same way: the real arm's largest |z| over the eight classes against
each coin draw's own largest, taken leave-one-out against the other
draws.  That is the max-statistic device rem:cnshift introduced, for
the same reason -- eight classes at |z| > 3 would fire on noise about
one time in fifty.

BACKS: Remark {#rem:cnclass} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Z1  THE GATE.  The whole band's excess kurtosis reproduces
      rem:cnkurt's POINT cnkurt_20 to three decimals, and the
      size-weighted average of the per-class second moments equals
      the whole band's to 1e-9.
  Z2  **The second moment carries class structure the coin cannot
      make.**  The real arm's largest |z| over the eight classes
      exceeds every one of the 32 coin draws' own largest.
  Z3  And the first moment does too, on the same comparison.  A mask
      that leaks into the signs would show as a per-class bias in G,
      which is where conj:L's M and G would fail to separate.
  Z4  The control is necessary: a typical coin draw shows apparent
      class structure of its own, so the median coin draw's largest
      |z| over classes is above 1.5.  If it were near zero the
      classes would be readable raw.

REFUTATION RULE (fixed before the run)

  Z1  REFUTED on either check; nothing below is reported.
  Z2  REFUTED if any coin draw's largest |z| reaches the real arm's.
      **Then six routes are closed and the last axis is closed with
      them**, and what this branch has measured is that C(N) is
      indistinguishable from a coin on mu's support by every device
      tried -- marginal, additive, multiplicative, positional and
      arithmetic. That is the outcome to state plainly if it comes,
      and it would end this line of measurement rather than suggest
      another.
  Z3  REFUTED on the same terms.  Z2 and Z3 can disagree; the second
      moment and the first are different claims, and a departure in
      one and not the other is a finding about which, not a wash.
  Z4  REFUTED below 1.5.  Then the control was not needed here, which
      is worth knowing because it would make a raw per-class table
      interpretable and this design over-cautious.

  WHAT THIS CANNOT DO.  One band and three primes.  Class structure
  at a modulus not in {3,5,7}, or at a band other than this one, is
  not addressed; and a departure found here would be located at a
  class, not attributed to any mechanism.
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
OUT = os.path.join(ROOT, "results", "audit_cn_class.txt")
SRC = os.path.join(ROOT, "results", "audit_cn_kurt_drift.txt")

BAND = 20
NMAX = 1 << (BAND + 1)
QS = (3, 5, 7)
DRAWS = 32
SEED = 20260828
DEC = 3
TOLW = 1e-9
MEDCAP = 1.5


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


def published_kurt():
    src = io.open(SRC, encoding="utf-8").read()
    m = re.search(r"^POINT cnkurt_%d ([-\d.]+)\s*$" % BAND, src, re.M)
    if not m:
        raise SystemExit("no POINT cnkurt_%d in %s" % (BAND, SRC))
    return float(m.group(1))


HEAD = [
    "STATISTIC: the first and second moments of",
    "           G(N) = C(N)/sqrt(V(N)) within each of the eight",
    "           classes of even N cut by which of %s divide N," % (QS,),
    "           against an ensemble of %d sign patterns on mu's" % DRAWS,
    "           support, compared by the largest |z| over classes",
    "           with each coin draw scored leave-one-out against the",
    "           rest, so each moment is one max-statistic and not",
    "           eight separate tests.",
    "FIELD: even N in (2^%d, 2^%d], %d of them, with Lambda and mu"
    % (BAND, BAND + 1, (1 << BAND) // 2),
    "       sieved once to %d and both convolutions taken by FFT" % NMAX,
    "       over the whole range. Averaged over signs the coin gives",
    "       E[G^2] = 1 in every class exactly, since V(N) is the exact",
    "       second moment at each N, so the null is no structure at",
    "       all rather than no difference between classes.",
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
    say("READ audit_cn_kurt_drift.txt POINT cnkurt_%d %.5f"
        % (BAND, pub))
    say("  the band this run has to reproduce, read from that file")
    say("PRINTBOUND audit_cn_class %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d classes from the primes %s, %d fresh sign patterns"
        % (1 << len(QS), QS, DRAWS))

    say("sieving to %d, which is 2^%d" % (NMAX, BAND + 1))
    lam, mu = sieves(NMAX)
    muf = mu.astype(np.float64)
    mu2 = muf ** 2
    del mu
    L = 1 << (2 * NMAX).bit_length()
    say("  transform length %d" % L)
    FL = np.fft.rfft(lam, L)
    C = np.fft.irfft(FL * np.fft.rfft(muf, L), L)[: NMAX + 1]
    V = np.fft.irfft(np.fft.rfft(lam ** 2, L) * np.fft.rfft(mu2, L),
                     L)[: NMAX + 1]
    del lam

    Ns = np.arange((1 << BAND) + 2, (1 << (BAND + 1)) + 1, 2,
                   dtype=np.int64)
    rootV = np.sqrt(np.maximum(V[Ns], 1e-300))
    del V

    lab = np.zeros(len(Ns), dtype=np.int64)
    for i, q in enumerate(QS):
        lab |= ((Ns % q) == 0).astype(np.int64) << i
    masks = [lab == k for k in range(1 << len(QS))]
    sizes = [int(m.sum()) for m in masks]

    def moments(g):
        return (np.array([float(g[m].mean()) for m in masks]),
                np.array([float((g[m] ** 2).mean()) for m in masks]))

    greal = C[Ns] / rootV
    gate_k = excess_kurt(greal)
    m1r, m2r = moments(greal)
    whole2 = float((greal ** 2).mean())
    wavg = float(np.dot(m2r, sizes) / sum(sizes))
    del C

    # -------------------------------------------------------------- Z1
    say()
    say("Z1  gate: the band and the class bookkeeping")
    say("  whole-band excess kurtosis %.5f against its %.5f"
        % (gate_k, pub))
    say("  size-weighted class second moment %.12f against the whole "
        "band's %.12f" % (wavg, whole2))
    z1 = (abs(round(gate_k, DEC) - round(pub, DEC))
          < 10.0 ** (-DEC) / 2 and abs(wavg - whole2) <= TOLW)
    say("  Z1 %s   (cap: %d decimals and %.0e)"
        % ("hold" if z1 else "REFUTED", DEC, TOLW))
    if not z1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rng = np.random.default_rng(SEED)
    c1 = np.zeros((DRAWS, len(masks)))
    c2 = np.zeros((DRAWS, len(masks)))
    for d in range(DRAWS):
        eps = rng.integers(0, 2, size=NMAX + 1).astype(np.float64) * 2 - 1
        Cc = np.fft.irfft(FL * np.fft.rfft(eps * mu2, L),
                          L)[: NMAX + 1]
        c1[d], c2[d] = moments(Cc[Ns] / rootV)
        del Cc

    def zsc(t, pool):
        return (t - pool.mean(axis=0)) / pool.std(axis=0, ddof=1)

    def maxes(real, pool):
        zr = zsc(real, pool)
        mc = []
        for d in range(DRAWS):
            rest = np.delete(pool, d, axis=0)
            mc.append(float(np.abs(zsc(pool[d], rest)).max()))
        return zr, float(np.abs(zr).max()), np.array(mc)

    z1r, max1, mc1 = maxes(m1r, c1)
    z2r, max2, mc2 = maxes(m2r, c2)

    say()
    say("      3|N 5|N 7|N        N      mean G   coin mean       z"
        "        E[G^2]   coin mean       z")
    for k in range(len(masks)):
        say("       %d   %d   %d  %8d  %+10.5f  %+10.5f  %6.2f  "
            "%+11.5f  %+10.5f  %6.2f"
            % (k & 1, (k >> 1) & 1, (k >> 2) & 1, sizes[k],
               m1r[k], c1[:, k].mean(), z1r[k],
               m2r[k], c2[:, k].mean(), z2r[k]))
        say("POINT classm2_%d %.5f" % (k, m2r[k]))
    say("SCALES 1")

    # -------------------------------------------------------------- Z2
    say()
    say("Z2  does the second moment escape the coin?")
    beat2 = int((mc2 >= max2).sum())
    say("  real largest |z| %.2f; the coin draws' own largest run "
        "%.2f to %.2f," % (max2, mc2.min(), mc2.max()))
    say("  and %d of %d reach it" % (beat2, DRAWS))
    z2 = beat2 == 0
    say("TSTAT class_m2_maxz %.2f" % max2)
    say("SPREAD class_m2_maxz %.5f" % float(mc2.std(ddof=1)))
    say("  Z2 %s   (cap: none of them)"
        % ("hold" if z2 else "REFUTED"))

    # -------------------------------------------------------------- Z3
    say()
    say("Z3  does the first moment?")
    beat1 = int((mc1 >= max1).sum())
    say("  real largest |z| %.2f; the coin draws' own largest run "
        "%.2f to %.2f," % (max1, mc1.min(), mc1.max()))
    say("  and %d of %d reach it" % (beat1, DRAWS))
    z3 = beat1 == 0
    say("TSTAT class_m1_maxz %.2f" % max1)
    say("SPREAD class_m1_maxz %.5f" % float(mc1.std(ddof=1)))
    say("  Z3 %s   (cap: none of them)"
        % ("hold" if z3 else "REFUTED"))

    # -------------------------------------------------------------- Z4
    say()
    say("Z4  was the control needed?")
    med = float(np.median(mc2))
    say("  the median coin draw's largest |z| over classes is %.2f "
        "on the second" % med)
    say("  moment, and %.2f on the first" % float(np.median(mc1)))
    z4 = med > MEDCAP
    say("  Z4 %s   (cap: above %.1f)"
        % ("hold" if z4 else "REFUTED", MEDCAP))

    say()
    say("=" * 70)
    say("Z1 %s  Z2 %s  Z3 %s  Z4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (z1, z2, z3, z4)))
    say()
    if z2 or z3:
        say("the last axis is not closed. Conditioning on which small "
            "primes divide")
        say("N finds something independent signs on mu's support "
            "cannot make, and")
        say("it is the one axis the coin was never able to follow -- "
            "its eps does")
        say("not know N.")
    else:
        say("six routes closed, and the last one was the only axis "
            "the coin could")
        say("not follow. Marginal, additive, multiplicative, "
            "positional and now")
        say("arithmetic: by every device this branch has tried, C(N) "
            "at these")
        say("ranges is indistinguishable from a field with mu's "
            "support and")
        say("independent signs. That ends this line of measurement "
            "rather than")
        say("pointing at another.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())

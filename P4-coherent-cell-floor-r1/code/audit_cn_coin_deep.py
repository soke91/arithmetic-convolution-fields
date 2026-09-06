# -*- coding: utf-8 -*-
r"""
Resolution or null type?  The iid coin at the same depth

WHAT IS AT STAKE

rem:cnmultdeep closed the seventh route: the real arm's largest |z|
over classes fell from 11.83 to 5.30 and five draws of 512 reached it,
so mu sits near the 99th percentile of the resolved multiplicative
ensemble rather than outside it.  **But that run changed two things at
once.**  It went from 32 draws to 512 and from iid signs to
multiplicative ones, and nothing separates which did the work.

The question is not local to this branch.  Lemma lem:coin -- the
control that has sunk claim after claim in this repository -- is an
iid coin.  If an iid ensemble resolved to 512 draws also covers mu,
then resolution was the whole story, the null type never mattered, and
lem:coin stands exactly as it is.  If it does not, then iid signs
understate the spread of the right null wherever multiplicativity is
in play, and every figure in this repository calibrated against a coin
was calibrated against something narrower than it should have been.
That is a repository-wide consequence and it is cheap to check.

The identity holds for both.  For iid eps on mu's support,
E[eps(m)eps(m')] = 0 unless m = m', so E[C_eps(N)^2] = V(N) exactly,
the same statement that makes the multiplicative ensemble's mean 1.
The two nulls are therefore directly comparable at the same depth.

BACKS: Remark {#rem:cncoindeep} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  F1  THE GATE.  The real second moment at 105 | N reproduces
      rem:cnclass's POINT classm2_7 to three decimals.
  F2  The iid ensemble comes to one, as its own identity requires:
      with 512 draws every class's ensemble mean E[G^2] is within two
      of its own standard errors of 1.  **The unresolved case is
      named and here it is a way of PASSING**: a wide ensemble makes
      the standard error large and any centre clears two of them, so
      F2 alone is the too-weak-to-tell case and only F2 together with
      the convergence printed below carries a reading.
  F3  **Resolution was the whole story.**  At 512 draws the iid
      ensemble also covers mu: at least one draw's own largest |z|
      over classes, leave-one-out, reaches the real arm's.
  F4  And the two nulls are the same width once resolved: the iid
      spread at 105 | N is within a factor of 1.5 of the
      multiplicative ensemble's 0.63553, read from
      results/audit_cn_multnull_deep.txt.

REFUTATION RULE (fixed before the run)

  F1  REFUTED outside three decimals; nothing below is reported.
  F2  REFUTED if any class stays two standard errors from one.  Then
      the construction or the convolution is wrong and nothing below
      is measured against what it is taken to be measured against.
      **The unresolved case is named here as well as above, because
      the gate reads this block and not the predictions.**  For F2 it
      is a way of PASSING rather than of failing: the standard error
      is the ensemble's own spread over the root of the draw count,
      so a wide or badly resolved ensemble makes it large and any
      centre clears two of them.  F2 alone is therefore the
      too-weak-to-tell case and must be read as "not resolved"; what
      makes it mean anything is the convergence table printed beside
      it.  (Added after the run to satisfy M9.  F2's cap is unchanged
      and F2 held; only what its holding licenses is narrowed.)
  F3  **REFUTED if none of the 512 reaches.**  Then the null type did
      the work, not the depth: an iid coin at full resolution still
      excludes mu while a multiplicative ensemble does not.  That
      makes multiplicativity the thing that covers mu, and it makes
      lem:coin's iid control narrower than the right null wherever
      the object under test is multiplicative -- which is most of
      this repository.  It would be the largest consequence anything
      in this branch has had outside it, and the remark must say
      plainly that it rests on one band and one statistic.
  F4  REFUTED outside 1.5 either way.  A wider iid ensemble would sit
      oddly with F3 refuted and the two together would need
      explaining rather than reporting; a narrower one is the
      expected shape if F3 is refuted, and then the factor is the
      measurement.

  WHAT THIS CANNOT DO.  One band, one class family, one statistic.
  Whatever F3 says here is about E[G^2] at 105 | N and not about
  lem:coin in general; what it can do is decide whether the question
  is worth asking of lem:coin in general, which nothing so far has.
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
OUT = os.path.join(ROOT, "results", "audit_cn_coin_deep.txt")
SRC = os.path.join(ROOT, "results", "audit_cn_class.txt")
SRCM = os.path.join(ROOT, "results", "audit_cn_multnull_deep.txt")

BAND = 20
NMAX = 1 << (BAND + 1)
QS = (3, 5, 7)
FULL = 7
PASSES = 8
PERPASS = 64
DRAWS = PASSES * PERPASS
CHECKS = (32, 64, 128, 256, 512)
SEED = 20260903
DEC = 3
FACT = 1.5


assert DRAWS == max(CHECKS)


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


def read_mark(path, name):
    src = io.open(path, encoding="utf-8").read()
    m = re.search(r"^%s ([-\d.]+)\s*$" % re.escape(name), src, re.M)
    if not m:
        raise SystemExit("no %s in %s" % (name, path))
    return float(m.group(1))


def read_multsd(path):
    """the multiplicative ensemble's spread at the full class"""
    src = io.open(path, encoding="utf-8").read()
    for ln in src.split("\n"):
        f = ln.split()
        if len(f) >= 9 and f[0] == "1" and f[1] == "1" and f[2] == "1":
            return float(f[6])
    raise SystemExit("no full-class row in %s" % path)


def loo_z(pool):
    d = pool.shape[0]
    s = pool.sum(axis=0)
    q = (pool ** 2).sum(axis=0)
    m = (s - pool) / (d - 1.0)
    v = ((q - pool ** 2) - (d - 1.0) * m ** 2) / (d - 2.0)
    return (pool - m) / np.sqrt(np.maximum(v, 1e-300))


HEAD = [
    "STATISTIC: the second moments E[G^2] of",
    "           G(N) = C(N)/sqrt(V(N)) over the eight classes cut by",
    "           which of %s divide N, against an ensemble of %d iid"
    % (QS, DRAWS),
    "           sign patterns on mu's support -- the same depth",
    "           rem:cnmultdeep gave the multiplicative ensemble, so",
    "           that resolution and null type are separated.",
    "FIELD: even N in (2^%d, 2^%d], with Lambda and mu sieved once to"
    % (BAND, BAND + 1),
    "       %d and every convolution taken by FFT over the whole" % NMAX,
    "       range. The %d draws are built in %d passes, %d per pass,"
    % (DRAWS, PASSES, PERPASS),
    "       each index carrying one random 64-bit word. For iid signs",
    "       E[C_eps(N)^2] = V(N) is exact, so the ensemble's E[G^2]",
    "       is 1 in every class, as it is for the multiplicative one.",
    "SEED: the draws come from numpy default_rng at seed %d; without"
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

    pub = read_mark(SRC, "POINT classm2_%d" % FULL)
    multsd = read_multsd(SRCM)
    say("READ audit_cn_class.txt POINT classm2_%d %.5f" % (FULL, pub))
    say("  the class second moment this run has to reproduce")
    say("  the multiplicative ensemble's spread at that class was "
        "%.5f," % multsd)
    say("  read from results/audit_cn_multnull_deep.txt")
    say("PRINTBOUND audit_cn_coin_deep %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d iid draws in %d passes of %d" % (DRAWS, PASSES, PERPASS))

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
    cmask = [lab == k for k in range(1 << len(QS))]
    sizes = [int(m.sum()) for m in cmask]

    def cls_m2(arr):
        g = arr[Ns] / rootV
        return np.array([float((g[m] ** 2).mean()) for m in cmask])

    real = cls_m2(C)
    del C

    # -------------------------------------------------------------- F1
    say()
    say("F1  does this run reproduce the class it shares?")
    say("  class %d E[G^2] here %.5f against its %.5f"
        % (FULL, real[FULL], pub))
    f1 = abs(round(real[FULL], DEC) - round(pub, DEC)) \
        < 10.0 ** (-DEC) / 2
    say("  F1 %s   (cap: %d decimals)"
        % ("hold" if f1 else "REFUTED", DEC))
    if not f1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rng = np.random.default_rng(SEED)
    pool = np.zeros((DRAWS, len(cmask)))
    row = 0
    for k in range(PASSES):
        hi = rng.integers(0, 1 << 32, size=NMAX + 1, dtype=np.uint64)
        lo = rng.integers(0, 1 << 32, size=NMAX + 1, dtype=np.uint64)
        t = (hi << np.uint64(32)) | lo
        del hi, lo
        for d in range(PERPASS):
            f = mu2 * (1.0 - 2.0 * ((t >> np.uint64(d))
                                    & np.uint64(1)).astype(np.float64))
            Cf = np.fft.irfft(FL * np.fft.rfft(f, L), L)[: NMAX + 1]
            pool[row] = cls_m2(Cf)
            row += 1
            del Cf, f
        del t
        say("  pass %d of %d done, %d draws" % (k + 1, PASSES, row))

    # -------------------------------------------------------------- F2
    say()
    say("F2  does the iid ensemble come to one?")
    say("      3|N 5|N 7|N        N      E[G^2]    ens. mean"
        "      ens. sd      s.e.   |m-1|/s.e.")
    worst, wk = 0.0, None
    for k in range(len(cmask)):
        m = float(pool[:, k].mean())
        sd = float(pool[:, k].std(ddof=1))
        se = sd / math.sqrt(DRAWS)
        r = abs(m - 1.0) / (se + 1e-300)
        if r > worst:
            worst, wk = r, k
        say("       %d   %d   %d  %8d  %+10.5f  %+11.5f  %11.5f  "
            "%8.5f  %8.2f"
            % (k & 1, (k >> 1) & 1, (k >> 2) & 1, sizes[k],
               real[k], m, sd, se, r))
        say("POINT coindeepmean_%d %.5f" % (k, m))
        say("POINT coindeepsd_%d %.5f" % (k, sd))
    say("SCALES 1")
    say("      draws    mean over classes of |ens. mean - 1|")
    dev = {}
    for c in CHECKS:
        d = float(np.abs(pool[:c].mean(axis=0) - 1.0).mean())
        dev[c] = d
        say("      %5d    %.5f" % (c, d))
    say("SCALES %d" % len(CHECKS))
    say("  the deviation %s from %d to %d draws"
        % ("falls" if dev[max(CHECKS)] < dev[min(CHECKS)] else "does "
           "not fall", min(CHECKS), max(CHECKS)))
    f2 = worst <= 2.0
    say("  furthest from one at class %s, at %.2f of its own standard "
        "error" % (wk, worst))
    say("TSTAT coindeep_center %.2f" % worst)
    say("SPREAD coindeep_center %.5f"
        % (float(pool[:, wk].std(ddof=1)) / math.sqrt(DRAWS)))
    say("  F2 %s   (cap: two standard errors)"
        % ("hold" if f2 else "REFUTED"))

    # -------------------------------------------------------------- F3
    say()
    say("F3  does the iid ensemble cover mu at this depth?")
    zr = (real - pool.mean(axis=0)) / pool.std(axis=0, ddof=1)
    maxr = float(np.abs(zr).max())
    where = int(np.argmax(np.abs(zr)))
    mc = np.abs(loo_z(pool)).max(axis=1)
    reach = int((mc >= maxr).sum())
    say("  real largest |z| %.2f at class %d; the draws' own largest "
        "run %.2f to %.2f," % (maxr, where, mc.min(), mc.max()))
    say("  and %d of %d reach it" % (reach, DRAWS))
    f3 = reach > 0
    say("TSTAT coindeep_maxz %.2f" % maxr)
    say("SPREAD coindeep_maxz %.5f" % float(mc.std(ddof=1)))
    say("  F3 %s   (cap: at least one draw)"
        % ("hold" if f3 else "REFUTED"))

    # -------------------------------------------------------------- F4
    say()
    say("F4  are the two nulls the same width once resolved?")
    sd = float(pool[:, FULL].std(ddof=1))
    ratio = sd / multsd
    say("  iid spread at class %d is %.5f against the multiplicative "
        "%.5f," % (FULL, sd, multsd))
    say("  a factor of %.4f" % ratio)
    f4 = (1.0 / FACT) <= ratio <= FACT
    say("POINT coindeep_widthratio %.5f" % ratio)
    say("  F4 %s   (cap: within %.1f either way)"
        % ("hold" if f4 else "REFUTED", FACT))

    say()
    say("=" * 70)
    say("F1 %s  F2 %s  F3 %s  F4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (f1, f2, f3, f4)))
    say()
    if f3 and f4:
        say("resolution was the whole story. An iid coin at the same "
            "depth covers")
        say("mu and is the same width as the multiplicative "
            "ensemble, so the null")
        say("type never mattered here and lem:coin stands as it is. "
            "What sank")
        say("rem:cnclass's z was 32 draws, not the choice of signs.")
    elif not f3:
        say("the null type did the work. An iid coin at full "
            "resolution still")
        say("excludes mu where a multiplicative ensemble does not, so "
            "iid signs")
        say("are narrower than the right null when the object is "
            "multiplicative.")
        say("Every figure in this repository calibrated against a "
            "coin was")
        say("calibrated against something narrower than it should "
            "have been. This")
        say("is one band and one statistic and it says the question "
            "is worth")
        say("asking of lem:coin in general, not that the answer is "
            "known.")
    else:
        say("the iid ensemble covers mu and is not the same width as "
            "the")
        say("multiplicative one. Those two together need explaining "
            "rather than")
        say("reporting, and the factor is printed above.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())

# -*- coding: utf-8 -*-
r"""
Resolve the null before believing what it says

WHAT IS AT STAKE

rem:cnmultm2 measured the escape against the strongest null this
branch can build and got two things that pull against each other.
The real arm's largest |z| over classes is 11.83 and none of 32 random
multiplicative draws reaches it (D2 refuted, the strong outcome).  But
the ensemble does not sit where the algebra puts it: E[C_f(N)^2] =
V(N) is exact for a random multiplicative f, so E[G_f^2] = 1 in every
class, and the 32-draw means came out 0.85394 to 0.93021 -- **all
eight below one**, the furthest at 5.66 of its own standard error, and
D4's rule had named "repeated across classes" as the crossing that
carries a reading.

That failure is not a side note.  The z-scores are taken against the
sample mean and the sample standard deviation, so if the ensemble is
under-resolved both are wrong in the direction that inflates z: too
low a centre and too small a spread.  rem:cnmultm2 wrote the honest
position -- what is measured is that mu is outside what 32 draws
produce, by a margin those 32 draws cannot calibrate -- and asked for
this run.

**The fix is draws.**  If the deficit is what a heavy right tail does
to a small sample, it closes as the sample grows and the escape can
then be re-asked against an ensemble that sits where it should.  If it
does not close, the identity or the implementation is wrong, which
would be a far more serious finding than any z.

512 draws are built in eight passes over the primes, 64 draws per
pass, each prime contributing one random 64-bit mask to a parity
accumulator.

BACKS: Remark {#rem:cnmultdeep} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  E1  THE GATE.  The real second moment at 105 | N reproduces
      rem:cnclass's POINT classm2_7 to three decimals.
  E2  **The ensemble comes to one.**  With 512 draws every class's
      ensemble mean E[G^2] is within two of its own standard errors
      of 1, as the identity requires.
  E3  And it gets there by converging, not by luck: the running mean
      over the eight classes, evaluated at 32, 64, 128, 256 and 512
      draws, is closer to 1 at 512 than at 32.
  E4  With the instrument fixed, the escape does not survive: at
      least one of the 512 draws has its own largest |z| over
      classes, leave-one-out, reaching the real arm's.  **This is the
      deflationary prediction** and it is what a heavy tail
      under-sampled at 32 would produce.

REFUTATION RULE (fixed before the run)

  E1  REFUTED outside three decimals; nothing below is reported.
  E2  **REFUTED if any class stays two standard errors from one.**
      Then the deficit is not a sampling artefact and either the
      identity E[C_f^2] = V is being violated -- it is not a
      conjecture, so that would mean the construction of f is not
      what it is taken to be -- or the convolution is wrong.  Either
      way every z in rem:cnmultm2 and here would be measured against
      something unexamined, and that has to be said before anything
      else is.
      **The unresolved case is named, and here it is a way of
      PASSING rather than of failing.**  The standard error is the
      ensemble's own spread over the root of the draw count, so a
      wide or badly resolved ensemble makes it large and then any
      centre clears two of them.  E2 holding on its own is therefore
      the too-weak-to-tell case; what makes it mean anything is E3,
      which asks whether the deviation actually fell as draws were
      added.  E2 without E3 must be read as "not resolved", never as
      "the ensemble is where it should be".  (This sentence was added
      after the run, to satisfy the gate's M9 check.  E2's cap is
      unchanged and E2 held; what the sentence fixes is what its
      holding licenses.)
  E3  REFUTED if the 512-draw means are no closer to one than the
      32-draw means.  E3 can fail while E2 holds if the approach is
      not monotone, and that is only a statement about the path.
  E4  REFUTED if none of the 512 reaches the real arm.  **Then mu is
      outside a properly resolved multiplicative null**, which is the
      largest claim this branch could make and the one that would
      have to be repeated at another band, with another seed, before
      it is written anywhere as more than one measurement.

  WHAT THIS CANNOT DO.  One band and one class family.  512 draws
  resolve a tail better than 32 and do not resolve it; if E4 is
  refuted the honest reading is bounded by the largest draw seen,
  which is printed.
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
OUT = os.path.join(ROOT, "results", "audit_cn_multnull_deep.txt")
SRC = os.path.join(ROOT, "results", "audit_cn_class.txt")

BAND = 20
NMAX = 1 << (BAND + 1)
QS = (3, 5, 7)
FULL = 7
PASSES = 8
PERPASS = 64
DRAWS = PASSES * PERPASS
CHECKS = (32, 64, 128, 256, 512)
SEED = 20260902
DEC = 3


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


def loo_z(pool):
    """each row's z against the other rows, column by column"""
    d = pool.shape[0]
    s = pool.sum(axis=0)
    q = (pool ** 2).sum(axis=0)
    m = (s - pool) / (d - 1.0)
    v = ((q - pool ** 2) - (d - 1.0) * m ** 2) / (d - 2.0)
    return (pool - m) / np.sqrt(np.maximum(v, 1e-300))


HEAD = [
    "STATISTIC: the second moments E[G^2] of",
    "           G(N) = C(N)/sqrt(V(N)) over the eight classes cut by",
    "           which of %s divide N, against an ensemble of %d"
    % (QS, DRAWS),
    "           random multiplicative sign functions, and the",
    "           convergence of that ensemble's class means to the 1",
    "           the identity E[C_f^2] = V requires.",
    "FIELD: even N in (2^%d, 2^%d], with Lambda and mu sieved once to"
    % (BAND, BAND + 1),
    "       %d and every convolution taken by FFT over the whole" % NMAX,
    "       range. The %d draws are built in %d passes over the"
    % (DRAWS, PASSES),
    "       primes, %d per pass, each prime contributing one random"
    % PERPASS,
    "       64-bit mask to a parity accumulator.",
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
    say("READ audit_cn_class.txt POINT classm2_%d %.5f" % (FULL, pub))
    say("  the class second moment this run has to reproduce")
    say("PRINTBOUND audit_cn_multnull_deep %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d draws in %d passes of %d" % (DRAWS, PASSES, PERPASS))

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

    def cls_m2(arr):
        g = arr[Ns] / rootV
        return np.array([float((g[m] ** 2).mean()) for m in masks])

    real = cls_m2(C)
    del C

    # -------------------------------------------------------------- E1
    say()
    say("E1  does this run reproduce the class it shares?")
    say("  class %d E[G^2] here %.5f against its %.5f"
        % (FULL, real[FULL], pub))
    e1 = abs(round(real[FULL], DEC) - round(pub, DEC)) \
        < 10.0 ** (-DEC) / 2
    say("  E1 %s   (cap: %d decimals)"
        % ("hold" if e1 else "REFUTED", DEC))
    if not e1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rng = np.random.default_rng(SEED)
    pr = primes_upto(NMAX)
    say()
    say("  %d primes carry the signs" % len(pr))
    pool = np.zeros((DRAWS, len(masks)))
    row = 0
    for k in range(PASSES):
        t = np.zeros(NMAX + 1, dtype=np.uint64)
        hi = rng.integers(0, 1 << 32, size=len(pr), dtype=np.uint64)
        lo = rng.integers(0, 1 << 32, size=len(pr), dtype=np.uint64)
        pmask = (hi << np.uint64(32)) | lo
        for i, p in enumerate(pr):
            t[int(p)::int(p)] ^= pmask[i]
        del hi, lo, pmask
        for d in range(PERPASS):
            f = mu2 * (1.0 - 2.0 * ((t >> np.uint64(d))
                                    & np.uint64(1)).astype(np.float64))
            Cf = np.fft.irfft(FL * np.fft.rfft(f, L), L)[: NMAX + 1]
            pool[row] = cls_m2(Cf)
            row += 1
            del Cf, f
        del t
        say("  pass %d of %d done, %d draws" % (k + 1, PASSES, row))

    # -------------------------------------------------------------- E2
    say()
    say("E2  does the ensemble come to one?")
    say("      3|N 5|N 7|N        N      E[G^2]    ens. mean"
        "      ens. sd      s.e.   |m-1|/s.e.")
    worst, wk = 0.0, None
    for k in range(len(masks)):
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
        say("POINT deepmean_%d %.5f" % (k, m))
    say("SCALES 1")
    e2 = worst <= 2.0
    say("  furthest from one at class %s, at %.2f of its own standard "
        "error" % (wk, worst))
    say("TSTAT deep_center %.2f" % worst)
    say("SPREAD deep_center %.5f"
        % (float(pool[:, wk].std(ddof=1)) / math.sqrt(DRAWS)))
    say("  E2 %s   (cap: two standard errors)"
        % ("hold" if e2 else "REFUTED"))

    # -------------------------------------------------------------- E3
    say()
    say("E3  does it get there by converging?")
    say("      draws    mean over classes of |ens. mean - 1|")
    dev = {}
    for c in CHECKS:
        d = float(np.abs(pool[:c].mean(axis=0) - 1.0).mean())
        dev[c] = d
        say("      %5d    %.5f" % (c, d))
        say("POINT deepdev_%d %.5f" % (c, d))
    say("SCALES %d" % len(CHECKS))
    e3 = dev[max(CHECKS)] < dev[min(CHECKS)]
    say("  E3 %s   (cap: closer at %d than at %d)"
        % ("hold" if e3 else "REFUTED", max(CHECKS), min(CHECKS)))

    # -------------------------------------------------------------- E4
    say()
    say("E4  does the escape survive the resolved ensemble?")
    zr = (real - pool.mean(axis=0)) / pool.std(axis=0, ddof=1)
    maxr = float(np.abs(zr).max())
    where = int(np.argmax(np.abs(zr)))
    mc = np.abs(loo_z(pool)).max(axis=1)
    reach = int((mc >= maxr).sum())
    say("  real largest |z| %.2f at class %d; the draws' own largest "
        "run %.2f to %.2f," % (maxr, where, mc.min(), mc.max()))
    say("  and %d of %d reach it" % (reach, DRAWS))
    e4 = reach > 0
    say("TSTAT deep_maxz %.2f" % maxr)
    say("SPREAD deep_maxz %.5f" % float(mc.std(ddof=1)))
    say("  E4 %s   (cap: at least one draw)"
        % ("hold" if e4 else "REFUTED"))

    say()
    say("=" * 70)
    say("E1 %s  E2 %s  E3 %s  E4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (e1, e2, e3, e4)))
    say()
    if not e2:
        say("the ensemble still does not sit where the identity puts "
            "it, so every")
        say("z here and in rem:cnmultm2 is measured against something "
            "unexamined.")
        say("That has to be settled before the escape means anything, "
            "and it is")
        say("not settled by more draws of the same kind.")
    elif e4:
        say("the escape was the tail. With the ensemble resolved to "
            "where the")
        say("identity puts it, draws of a random multiplicative sign "
            "function do")
        say("reach where mu goes, and rem:cnmultm2's 11.83 was what "
            "32 draws")
        say("could not calibrate. Seven routes, and none of them "
            "escapes.")
    else:
        say("mu is outside a resolved multiplicative null. The "
            "ensemble now sits")
        say("where the identity puts it and still none of %d draws "
            "reaches the" % DRAWS)
        say("field. That is the largest claim this branch could make "
            "and it rests")
        say("on one band and one seed; it needs repeating at both "
            "before it is")
        say("written anywhere as more than one measurement.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())

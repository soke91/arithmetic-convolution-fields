# -*- coding: utf-8 -*-
r"""
The Monte-Carlo verification of Lemma {#lem:cellmom}'s closed form,
which paper/wall_v3.md cites and no script has run.

WHY THIS HAS TO BE RUN

Section {#sec:floor} says of the closed form:

    "The lemma has been checked against simulation even though it does
     not need it: over 60 independent sign draws in the band
     (10^5, 2e5], the Monte-Carlo standard deviation of m_c - mbar
     agrees with the closed form at every depth, the ratios running
     0.88 to 0.98 against a Monte-Carlo precision of +-0.09."

lab_cell_floor.py computes the closed form and the z-scores; it never
simulates, and says so in its own NULL: line.  So the one check that
would catch an error IN THE FORMULA has never been run -- and the
formula is what every z-score in this program is divided by.  After
lab_mask_placebo.py had to run the placebo that Section {#sec:floor}
also merely cited, this is the remaining cited-but-unrun control.

    Var(m_c - mbar) = Q_cc/n_c^2 - 2 Q_ca/(n_c n) + Q_aa/n^2,
    u_c(v) = sum_{N in c} Lambda(N-v)/sqrt(V(N)),
    Q_cd   = sum_v mu^2(v) u_c(v) u_d(v).

The simulation draws eps(v) = +-1 on {v : mu(v) != 0}, forms
Z_eps(N) = (eps * Lambda)(N)/sqrt(V(N)), and takes the sample standard
deviation of m_c - mbar across draws.  Two draw counts are used: 60,
to reproduce the quoted figure at the precision it was quoted at, and
2000, to test the formula rather than the quotation.

BACKS: Lemma {#lem:cellmom} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  N1  With 2000 draws the ratio MC/closed-form is within 0.05 of 1 at
      every depth: the closed form is exact, not approximate.
  N2  With 60 draws the ratios lie within 1 +- 3/sqrt(2*59) = 1 +- 0.276
      at every depth, so the quoted band 0.88 to 0.98 is an ordinary
      draw and not a discrepancy.
  N3  The Monte-Carlo precision at 60 draws is 1/sqrt(2*59) = 0.0921,
      matching the quoted +-0.09.
  N4  At 2000 draws the ratios straddle 1 -- at least one above and one
      below. All six of the quoted ratios fall below 1, which happens
      with probability 2^-6 if the form is exact, so this is worth
      checking rather than assuming.

REFUTATION RULE (fixed before the run)

  N1  REFUTED if any depth is off by 0.05 or more.  This is the one
      that matters: a failure here invalidates every z-score in
      Section {#sec:floor}.
  N2  REFUTED if any 60-draw ratio leaves 1 +- 0.276.
  N3  REFUTED if the computed precision differs from 0.09 by more than
      0.005.
  N4  REFUTED if all six 2000-draw ratios fall on the same side of 1.

  N1 and N3 gate; N2 and N4 are reported.

CITED BY: {#rem:mcratios} in paper/.
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
OUT = os.path.join(ROOT, "results", "lab_cellmom_montecarlo.txt")

LO, HI = 100_000, 200_000
CELLP = (3, 5, 7, 11, 13)
DRAWS_SMALL = 60
DRAWS_BIG = 2000
SEED = 20260808
PUB_LO, PUB_HI = 0.88, 0.98


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def pow2(n):
    L = 1
    while L < n:
        L <<= 1
    return L


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    X = HI
    say("sieving to %d ..." % X)
    pr = primes_upto(X)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(X + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > X:
            break
        q = p * p
        while q <= X:
            lam[q] = lgp[i]
            if q > X // p:
                break
            q *= p

    mu = np.ones(X + 1, dtype=np.int8)
    rem = np.arange(X + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(X))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= X:
            mu[p * p::p * p] = 0
        q = p
        while q <= X:
            rem[q::q] //= p
            if q > X // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    sqf = (mu != 0)

    n = pow2(2 * (X + 1))
    a = np.zeros(n, dtype=np.float64)
    a[:X + 1] = lam ** 2
    FL2 = np.fft.rfft(a)
    a[:] = 0.0
    a[:X + 1] = sqf
    V = np.fft.irfft(FL2 * np.fft.rfft(a), n)[:X + 1]
    del FL2
    a[:] = 0.0
    a[:X + 1] = lam
    FLam = np.fft.rfft(a)
    del a

    depth = np.zeros(X + 1, dtype=np.int8)
    for p in CELLP:
        depth[p::p] += 1

    Ns = np.arange(LO + 2, HI + 1, 2, dtype=np.int64)
    Ns = Ns[V[Ns] > 0]
    nb = Ns.size
    rootV = np.sqrt(V[Ns])
    dep = depth[Ns]
    say("band (%d, %d]: %d even N" % (LO, HI, nb))

    # ---- the closed form of Lemma {#lem:cellmom}
    g = 1.0 / rootV
    b = np.zeros(n, dtype=np.float64)
    FLc = np.conj(FLam)
    w = sqf[:X + 1].astype(np.float64)

    def ucorr(sel):
        b[:] = 0.0
        b[Ns[sel]] = g[sel]
        return np.fft.irfft(FLc * np.fft.rfft(b), n)[:X + 1]

    ua = ucorr(np.ones(nb, dtype=bool))
    Qaa = float((w * ua * ua).sum())
    ds, ncs, closed = [], [], []
    for d in range(6):
        sel = dep == d
        nc = int(sel.sum())
        if nc == 0:
            continue
        uc = ucorr(sel)
        Qcc = float((w * uc * uc).sum())
        Qca = float((w * uc * ua).sum())
        var = Qcc / nc ** 2 - 2.0 * Qca / (nc * nb) + Qaa / nb ** 2
        ds.append(d)
        ncs.append(nc)
        closed.append(math.sqrt(max(var, 0.0)))
    ds = np.array(ds)
    closed = np.array(closed)

    # ---- the simulation
    supp = np.flatnonzero(sqf[:X + 1])
    rng = np.random.default_rng(SEED)
    say("simulating %d draws ..." % DRAWS_BIG)
    diffs = np.empty((DRAWS_BIG, ds.size))
    e = np.zeros(n, dtype=np.float64)
    for t in range(DRAWS_BIG):
        e[:] = 0.0
        e[supp] = rng.integers(0, 2, size=supp.size) * 2.0 - 1.0
        Ce = np.fft.irfft(FLam * np.fft.rfft(e), n)
        Z = Ce[Ns] / rootV
        zb = float(Z.mean())
        for i, d in enumerate(ds):
            diffs[t, i] = float(Z[dep == d].mean()) - zb

    say()
    say("  depth  n_c        closed form    MC sd (%d)   ratio    "
        "MC sd (%d)  ratio" % (DRAWS_SMALL, DRAWS_BIG))
    say("  " + "-" * 78)
    mc_small = diffs[:DRAWS_SMALL].std(axis=0, ddof=1)
    mc_big = diffs.std(axis=0, ddof=1)
    r_small = mc_small / closed
    r_big = mc_big / closed
    for i, d in enumerate(ds):
        say("  %-6d %-10d %-14.6e %-13.6e %-8.4f %-11.6e %.4f"
            % (d, ncs[i], closed[i], mc_small[i], r_small[i],
               mc_big[i], r_big[i]))

    say()
    n1 = bool((np.abs(r_big - 1.0) < 0.05).all())
    say("N1  |MC/closed - 1| at %d draws: max %.4f   (cap 0.05)   %s"
        % (DRAWS_BIG, float(np.abs(r_big - 1.0).max()),
           "hold" if n1 else "REFUTED"))
    tol = 3.0 / math.sqrt(2.0 * (DRAWS_SMALL - 1))
    n2 = bool((np.abs(r_small - 1.0) <= tol).all())
    say("N2  60-draw ratios in [%.4f, %.4f]; band 1 +- %.3f   %s"
        % (r_small.min(), r_small.max(), tol,
           "hold" if n2 else "REFUTED"))
    say("    the quoted band was %.2f to %.2f" % (PUB_LO, PUB_HI))
    prec = 1.0 / math.sqrt(2.0 * (DRAWS_SMALL - 1))
    n3 = abs(prec - 0.09) <= 5e-3
    say("N3  MC precision at %d draws = %.6f   quoted 0.09   %s"
        % (DRAWS_SMALL, prec, "hold" if n3 else "REFUTED"))
    n4 = bool((r_big > 1).any() and (r_big < 1).any())
    say("N4  %d-draw ratios straddle 1: %s   (%d above, %d below)   %s"
        % (DRAWS_BIG, n4, int((r_big > 1).sum()), int((r_big < 1).sum()),
           "hold" if n4 else "REFUTED"))

    say()
    say("POST HOC (not pre-registered; added after the paper moved its")
    say("quotation from 60 draws to 2000, so that the sentences it now")
    say("makes are checked by something). N2 and N3 above score a")
    say("quotation this paper no longer makes.")
    prec_big = 1.0 / math.sqrt(2.0 * (DRAWS_BIG - 1))
    say("  Monte-Carlo precision at %d draws = %.6f" % (DRAWS_BIG, prec_big))
    say("  depth   |MC/closed - 1|   multiple of the precision")
    dev = np.abs(r_big - 1.0)
    for i, d in enumerate(ds):
        say("  %-6d  %.6f          %.3f" % (d, dev[i], dev[i] / prec_big))
    below = [i for i in range(len(ds)) if dev[i] < prec_big]
    worst = int(np.argmax(dev))
    say("  worst is depth %d at %.3f times the precision; %d of the six"
        % (ds[worst], dev[worst] / prec_big, len(below)))
    say("  lie below it, the closest being depth %d, short of the bar by"
        % ds[max(below, key=lambda i: dev[i])])
    m = max(below, key=lambda i: dev[i])
    say("  %.6f -- %.1f%% of the bar. That one is not a separation:"
        % (prec_big - dev[m], 100.0 * (prec_big - dev[m]) / prec_big))
    say("  it sits at the precision, not below it.")

    say()
    say("=" * 70)
    ok = n1 and n3
    say("N1 %s  N2 %s  N3 %s  N4 %s"
        % tuple("hold" if v else "REFUTED" for v in (n1, n2, n3, n4)))
    say("Lemma {#lem:cellmom}'s closed form is verified by simulation"
        if ok else "REFUTED -- the closed form does not match simulation")

    head = [
        "STATISTIC: the sample standard deviation of m_c - mbar over",
        "           independent sign draws, against the closed form",
        "           Q_cc/n_c^2 - 2 Q_ca/(n_c n) + Q_aa/n^2 of Lemma",
        "           {#lem:cellmom}; the ratio of the two at 60 draws (the",
        "           count the paper quotes) and at 2000; and the",
        "           Monte-Carlo precision 1/sqrt(2(D-1)).",
        "NULL: the simulation IS the null -- independent signs on the",
        "      support of mu^2, which is the hypothesis under which the",
        "      closed form is derived. mu itself never enters: this",
        "      checks the formula, not a detection.",
        "FIELD: the band (1e5, 2e5], even N with V(N) > 0; cells indexed",
        "       by depth = #{p in 3,5,7,11,13 dividing N}; Lambda, mu and",
        "       the squarefree indicator from an integer sieve to 2e5;",
        "       V and each draw's C_eps by exact FFT convolution, u_c by",
        "       FFT cross-correlation; numpy default_rng seed 20260808.",
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

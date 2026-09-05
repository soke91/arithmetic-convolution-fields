# -*- coding: utf-8 -*-
r"""
paper/wall_v3.md, Section {#sec:floor} -- Lemma {#lem:cellmom},
Proposition {#prop:coh}, and the cell means measured against the floor.

WHAT IS UNDER TEST

Cells are indexed by depth d = #{p in {3,5,7,11,13} : p | N}, and
Z(N) = C(N)/sqrt(V(N)) with C(N) = sum_{v<N} mu(v) Lambda(N-v).  For a
cell c of size n_c inside a band a of size n,

    u_c(v) = sum_{N in c} Lambda(N-v)/sqrt(V(N)),
    Q_cd   = sum_v mu^2(v) u_c(v) u_d(v),
    Var(m_c - mbar) = Q_cc/n_c^2 - 2 Q_ca/(n_c n) + Q_aa/n^2,

the last being Lemma [lem:cellmom], exact under independent signs.  The
paper prints:

  (a) Var/(Q_cc/n_c^2) = 0.113, 0.016, 0.289, 0.551 at depths 0,1,3,5 in
      the top octave, "reproducing to three digits the same ratios at a
      quarter of that N";
  (b) Q_cc/n_c^2 runs 0.124 to 0.307 by depth at N = 1.6e7, against the
      heuristic sum_v mu^2(v)/V ~ 0.049;
  (c) fitting se ~ N^{-b} to the exact floor across eight octaves gives
      b = 0.0395, 0.0397, 0.0394 at depths 2,1,0, against the apparent
      exponent 1/(2<log N>) = 0.0358;
  (d) max_c |z_c| runs 9.1 to 13.0 over every octave from 6.25e4 to
      1.6e7, and at the top octave depths 3,4,5 sit at z = -1.6, -4.5,
      -9.1 while depths 0,1,2 sit at +0.0, +0.7, -0.1.

No script for any of it exists here.  u_c is a cross-correlation, so
every Q is computed by FFT and nothing is simulated.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  C1  Var/(Q_cc/n_c^2) = 0.113, 0.016, 0.289, 0.551 at depths 0,1,3,5,
      top octave (8e6, 1.6e7].
  C2  The same four ratios at the octave (2e6, 4e6] agree with the top
      octave's to three digits.
  C3  min_d Q_cc/n_c^2 = 0.124 and max_d = 0.307 at the top octave, and
      the heuristic (sum_{v<N} mu^2(v))/V(N) = 0.049.
  C4  b = 0.0395, 0.0397, 0.0394 at depths 2,1,0 from eight octaves, and
      1/(2<log N>) = 0.0358.
  C5  max_c |z_c| lies in [9.1, 13.0] at every one of the eight octaves.
  C6  Top-octave z by depth 0..5 = +0.0, +0.7, -0.1, -1.6, -4.5, -9.1.

REFUTATION RULE (fixed before the run)

  C1  REFUTED if any of the four differs by more than 0.001.
  C2  REFUTED if any of the four differs from its top-octave value by
      more than 0.001.
  C3  REFUTED if the min or max differs by more than 0.0005, or if the
      heuristic differs from 0.049 by more than 0.0005.
  C4  REFUTED if any b differs by more than 0.00005, or if
      1/(2<log N>) differs from 0.0358 by more than 0.00005.
  C5  REFUTED if any octave's max |z_c| falls outside [9.1, 13.0] by
      more than 0.05.
  C6  REFUTED if any of the six differs by more than 0.05.

  All six gate.
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
OUT = os.path.join(ROOT, "results", "lab_cell_floor.txt")

X = 16_000_000
CELLP = (3, 5, 7, 11, 13)
OCTS = [(X >> (i + 1), X >> i) for i in range(8)][::-1]   # 8 octaves
PUB_C1 = {0: 0.113, 1: 0.016, 3: 0.289, 5: 0.551}
PUB_C3 = (0.124, 0.307, 0.049)
PUB_C4 = {2: 0.0395, 1: 0.0397, 0: 0.0394}
PUB_C4_APP = 0.0358
PUB_C5 = (9.1, 13.0)
PUB_C6 = [0.0, 0.7, -0.1, -1.6, -4.5, -9.1]


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

    say("convolving V = mu^2 * Lambda^2 and C = mu * Lambda ...")
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
    FL = np.fft.rfft(a)
    a[:] = 0.0
    a[:X + 1] = mu
    C = np.fft.irfft(FL * np.fft.rfft(a), n)[:X + 1]
    del a, FL

    depth = np.zeros(X + 1, dtype=np.int8)
    for p in CELLP:
        depth[p::p] += 1

    say("done. computing the exact floor octave by octave ...")
    say()
    say("  octave        depth  n_c      Q_cc/n_c^2   Var          "
        "ratio     se           z_c")
    say("  " + "-" * 92)

    floors = {}
    rows_top = {}
    q_top = {}
    z_top = {}
    maxz = []
    for lo, hi in OCTS:
        T = hi
        m = pow2(2 * (T + 1))
        b = np.zeros(m, dtype=np.float64)
        b[:T + 1] = lam[:T + 1]
        FLam = np.conj(np.fft.rfft(b))
        Ns = np.arange(lo + 2 - (lo % 2), hi + 1, 2, dtype=np.int64)
        Ns = Ns[V[Ns] > 0]
        g = 1.0 / np.sqrt(V[Ns])
        Z = C[Ns] / np.sqrt(V[Ns])
        nb = Ns.size
        zbar = float(Z.mean())

        def ucorr(sel):
            b[:] = 0.0
            b[Ns[sel]] = g[sel]
            return np.fft.irfft(FLam * np.fft.rfft(b), m)[:T + 1]

        w = sqf[:T + 1].astype(np.float64)
        ua = ucorr(np.ones(nb, dtype=bool))
        Qaa = float((w * ua * ua).sum())
        best = 0.0
        for d in range(6):
            sel = depth[Ns] == d
            nc = int(sel.sum())
            if nc == 0:
                continue
            uc = ucorr(sel)
            Qcc = float((w * uc * uc).sum())
            Qca = float((w * uc * ua).sum())
            var = Qcc / nc ** 2 - 2.0 * Qca / (nc * nb) + Qaa / nb ** 2
            se = math.sqrt(max(var, 0.0))
            mc = float(Z[sel].mean())
            z = (mc - zbar) / se if se > 0 else float("nan")
            floors.setdefault(d, []).append((0.5 * (lo + hi), se))
            best = max(best, abs(z))
            if hi == X:
                rows_top[d] = var / (Qcc / nc ** 2)
                q_top[d] = Qcc / nc ** 2
                z_top[d] = z
            # Every octave prints its depths.  It used to print two of
            # the eight, and Section {#sec:notclaimed} then said the
            # shallow depths clear no floor "at any scale measured" --
            # a sentence whose range was two, not eight, and whose
            # largest max|z| octave was not among them.  The numbers
            # were already computed; only the printing was narrow.
            if True:
                say("  (%9d,%9d] %-6d %-8d %-12.6f %-12.6e %-9.4f "
                    "%-12.4e %+.2f"
                    % (lo, hi, d, nc, Qcc / nc ** 2, var,
                       var / (Qcc / nc ** 2), se, z))
            if hi == 4_000_000:
                q_top.setdefault(("q", d), var / (Qcc / nc ** 2))
        maxz.append((hi, best))
        del b, FLam

    say()
    say("C1  Var / (Q_cc/n_c^2) at the top octave")
    c1 = True
    for d in (0, 1, 3, 5):
        got = rows_top.get(d, float("nan"))
        ok = abs(got - PUB_C1[d]) <= 1e-3
        c1 = c1 and ok
        say("    depth %d: %.6f   published %.3f   %s"
            % (d, got, PUB_C1[d], "ok" if ok else "MISMATCH"))
    say("  C1 %s" % ("hold" if c1 else "REFUTED"))

    say()
    say("C2  the same four at the octave (2e6, 4e6]")
    c2 = True
    for d in (0, 1, 3, 5):
        got = q_top.get(("q", d), float("nan"))
        ok = abs(got - rows_top.get(d, float("nan"))) <= 1e-3
        c2 = c2 and ok
        say("    depth %d: %.6f   top octave %.6f   %s"
            % (d, got, rows_top.get(d, float("nan")),
               "ok" if ok else "MISMATCH"))
    say("  C2 %s" % ("hold" if c2 else "REFUTED"))

    say()
    qs = [q_top[d] for d in range(6) if d in q_top]
    heur = float(sqf[:X + 1].sum()) / V[X]
    c3 = (abs(min(qs) - PUB_C3[0]) <= 5e-4
          and abs(max(qs) - PUB_C3[1]) <= 5e-4
          and abs(heur - PUB_C3[2]) <= 5e-4)
    say("C3  Q_cc/n_c^2 by depth: %s"
        % ", ".join("%.6f" % q for q in qs))
    say("    min %.6f (pub 0.124), max %.6f (pub 0.307), heuristic "
        "%.6f (pub 0.049)   %s"
        % (min(qs), max(qs), heur, "hold" if c3 else "REFUTED"))
    # The three differences, against the 5e-4 the rule allows.  The
    # values alone do not say which clause decided, and one of them
    # misses by less than the last digit of the published number it
    # is compared against -- 0.049540 against 0.049 is 5.40e-4 and
    # fails, which reads as a pass if only the values are printed.
    # A paper sentence said "the maximum and the heuristic both
    # hold" on exactly that reading.
    for _nm, _v, _p in (("min", min(qs), PUB_C3[0]),
                        ("max", max(qs), PUB_C3[1]),
                        ("heuristic", heur, PUB_C3[2])):
        _d = abs(_v - _p)
        say("      %-10s |%.6f - %.6f| = %.2e   cap 5.00e-04   %s"
            % (_nm, _v, _p, _d, "ok" if _d <= 5e-4 else "OVER"))

    say()
    say("C4  se ~ N^{-b} across the eight octaves")
    c4 = True
    for d in (2, 1, 0):
        pts = floors.get(d, [])
        if len(pts) < 8:
            c4 = False
            say("    depth %d: only %d octaves" % (d, len(pts)))
            continue
        xs = np.log(np.array([p[0] for p in pts]))
        ys = np.log(np.array([p[1] for p in pts]))
        bfit = -float(np.polyfit(xs, ys, 1)[0])
        ok = abs(bfit - PUB_C4[d]) <= 5e-5
        c4 = c4 and ok
        say("    depth %d: b = %.6f   published %.4f   %s"
            % (d, bfit, PUB_C4[d], "ok" if ok else "MISMATCH"))

    # The three depths above are the ones with a published value.  The
    # other three are fitted too and not gated: the (log N)^{-1/2} form
    # is a statement about cells with many N in them, and depth 5 has
    # 2 to 266.  Printing all six is what lets a reader see where the
    # form holds and where it stops -- the paper said "any cell mean"
    # while the file showed only the three that agree.
    say("    the remaining depths, fitted the same way and not gated:")
    for d in (3, 4, 5):
        pts = floors.get(d, [])
        if len(pts) < 8:
            continue
        xs = np.log(np.array([p[0] for p in pts]))
        ys = np.log(np.array([p[1] for p in pts]))
        bfit = -float(np.polyfit(xs, ys, 1)[0])
        ncs = [p[2] for p in pts] if len(pts[0]) > 2 else []
        say("    depth %d: b = %.6f   shrink over the range %.3f"
            % (d, bfit, (max(p[0] for p in pts)
                         / min(p[0] for p in pts)) ** bfit))
    mlog = float(np.mean([math.log(0.5 * (lo + hi)) for lo, hi in OCTS]))
    app = 1.0 / (2.0 * mlog)
    okapp = abs(app - PUB_C4_APP) <= 5e-5
    c4 = c4 and okapp
    say("    1/(2<log N>) = %.6f   published 0.0358   %s"
        % (app, "ok" if okapp else "MISMATCH"))
    say("  C4 %s" % ("hold" if c4 else "REFUTED"))

    say()
    say("C5  max_c |z_c| by octave")
    c5 = True
    for hi, b_ in maxz:
        ok = PUB_C5[0] - 0.05 <= b_ <= PUB_C5[1] + 0.05
        c5 = c5 and ok
        say("    top %-10d  max|z| = %7.3f   %s"
            % (hi, b_, "in [9.1,13.0] +- 0.05" if ok else "OUT"))
    say("  C5 %s" % ("hold" if c5 else "REFUTED"))

    say()
    say("C6  top-octave z by depth")
    c6 = True
    for d in range(6):
        got = z_top.get(d, float("nan"))
        ok = abs(got - PUB_C6[d]) <= 0.05
        c6 = c6 and ok
        say("    depth %d: z = %+.4f   published %+.1f   %s"
            % (d, got, PUB_C6[d], "ok" if ok else "MISMATCH"))
    say("  C6 %s" % ("hold" if c6 else "REFUTED"))

    say()
    say("=" * 74)
    ok = c1 and c2 and c3 and c4 and c5 and c6
    say("C1 %s  C2 %s  C3 %s  C4 %s  C5 %s  C6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (c1, c2, c3, c4, c5, c6)))
    say("Section {#sec:floor} reproduces" if ok else "REFUTED")

    head = [
        "STATISTIC: the exact cell-mean floor of Lemma {#lem:cellmom},",
        "           Var(m_c - mbar) = Q_cc/n_c^2 - 2 Q_ca/(n_c n)",
        "           + Q_aa/n^2, with u_c(v) = sum_{N in c}",
        "           Lambda(N-v)/sqrt(V(N)) computed as an FFT",
        "           cross-correlation and Q_cd = sum_v mu^2(v)u_c u_d;",
        "           the ratio Var/(Q_cc/n_c^2); Q_cc/n_c^2 against the",
        "           heuristic (sum mu^2)/V; the exponent b in se ~ N^{-b}",
        "           across eight octaves; and z_c = (m_c - mbar)/se with",
        "           Z(N) = C(N)/sqrt(V(N)).",
        "FIELD: even N; eight octaves from (6.25e4, 1.25e5] up to",
        "       (8e6, 1.6e7]; cells indexed by depth = #{p in 3,5,7,11,13",
        "       dividing N}; Lambda, mu and the squarefree indicator from",
        "       an integer sieve to 1.6e7; V = mu^2 * Lambda^2 and",
        "       C = mu * Lambda by exact FFT convolution.",
        'NULL: built in. The floor of Lemma {#lem:cellmom} IS the',
        '      independent-sign variance, computed exactly rather than',
        '      simulated, so every z_c is already quoted against its coin',
        '      null. The second control, the cell permutation of Lemma',
        '      {#lem:placebo}, is not run here.',
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

# -*- coding: utf-8 -*-
r"""P4 — the second cell shape, measured EXACTLY.

WHAT THIS IS

lab_secondcell2_predict.py registered T1, T2 and T3 on an exact
statistic, before this file existed, and read its predictions from the
first seal (lab_secondcell_predict.py, commit 3b1b7d1) rather than
recomputing them.  This file makes the measurement and scores those
three rules.

There is no sampling here.  A cell is a set of positions in the band,
so the pair average of S_2 over it is

    E_same,c[S_2] = sum_h S_2(h) (1_c * 1_c)(h) / (n_c (n_c - 1)),

and the autocorrelation of the indicator is one FFT per cell over the
band's 10^6 even positions.  No seed, no draw, no standard error --
which is the whole reason for a second registration, the first having
been refused by a single draw's luck rather than by its prediction.

The autocorrelation counts are integers and are rounded to integers
after the transform; at length 2^21 the transform's error on counts of
order 10^6 is far below a half, and the rounding is checked by
requiring sum_k 2 ac[k] = n_c (n_c - 1) exactly.

FIELD: band (2e6, 4e6]; even N; cells by depth over 3,5,7; S_2 sieved
       to 2e6 over every prime, with 2 C_2 accumulated to 4e6; exact
       autocorrelation by rFFT at length 2^21.  Predictions and rules
       read from results/lab_secondcell2_predict.txt and
       results/lab_secondcell_predict.txt.

BACKS: Proposition {#prop:scaleinv} and Measurement {#meas:maincoef}
       in deploy/papers/P4-coherent-cell-floor.tex -- as the test of
       whether the main-term computation depends on the cell shape.
"""

import io
import math
import os
import sys

import numpy as np
from fractions import Fraction as F

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "lab_secondcell2.txt")
SEAL1 = os.path.join(ROOT, "results", "lab_secondcell_predict.txt")
SEAL2 = os.path.join(ROOT, "results", "lab_secondcell2_predict.txt")

CELLP2 = (3, 5, 7)
BAND = (2_000_000, 4_000_000)
HMAX = 2_000_000
CLIM = 4_000_000
EXACT_EALL = 2.0
TOL_T2 = 1e-4
TOL_T3 = 0.15

lines = []


def say(s=""):
    lines.append(s)
    print(s)


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def read_seal1():
    src = io.open(SEAL1, encoding="utf-8", errors="replace").read()
    corr, unc, ms, full = {}, {}, {}, {}
    b = src.split("depth  corrected prediction")
    if len(b) > 1:
        for ln in b[1].splitlines()[1:]:
            f = ln.split()
            if len(f) == 2 and f[0].isdigit():
                corr[int(f[0])] = float(f[1])
            elif corr:
                break
    b = src.split("m (classes)  n_c approx   uncorrected prediction")
    if len(b) > 1:
        for ln in b[1].splitlines()[1:]:
            f = ln.split()
            if len(f) == 4 and f[0].isdigit():
                unc[int(f[0])] = float(f[3])
                ms[int(f[0])] = int(f[1])
            elif unc:
                break
    b = src.split("depth  m      n_c per class  deficit      correction"
                  "   diluted")
    if len(b) > 1:
        for ln in b[1].splitlines()[1:]:
            f = ln.split()
            if len(f) == 6 and f[0].isdigit():
                full[int(f[0])] = float(f[4])
            elif full:
                break
    return corr, unc, ms, full


def main():
    if not os.path.exists(SEAL2):
        say("the second registration is not present -- refusing to run")
        return 1
    corr, unc, ms, full = read_seal1()
    if len(corr) != 4:
        say("the first seal was not readable -- refusing to run")
        return 1
    say("SECOND CELL SHAPE, MEASURED EXACTLY.  S = %s, band (%d, %d]"
        % (", ".join(str(p) for p in CELLP2), BAND[0], BAND[1]))
    say("rules T1, T2, T3 registered in %s" % os.path.basename(SEAL2))
    say("predictions read from %s" % os.path.basename(SEAL1))
    say()

    pr = primes_upto(HMAX)
    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    S2 = np.zeros(HMAX + 1, dtype=np.float64)
    S2[2::2] = twin
    for p in pr:
        p = int(p)
        if p == 2:
            continue
        S2[p::p] *= (p - 1.0) / (p - 2.0)
    S2[1::2] = 0.0
    S2[0] = 0.0
    say("  2 C_2 = %.8f;  S_2(2) = %.6f;  S_2(6) = %.6f (= 2 S_2(2))"
        % (twin, S2[2], S2[6]))

    Ns = np.arange(BAND[0] + 2, BAND[1] + 1, 2, dtype=np.int64)
    L = Ns.size
    depth = np.zeros(L, dtype=np.int8)
    for p in CELLP2:
        depth += (Ns % p == 0)
    nfft = 1
    while nfft < 2 * L:
        nfft *= 2
    say("  band positions %d, transform length %d" % (L, nfft))
    say()
    # shifts h = 2k for k = 1 .. L-1
    kk = np.arange(1, L, dtype=np.int64)
    w = S2[2 * kk]

    say("  depth  n_c        E_same,c[S_2] exact   D'_c = that - 2")
    say("  " + "-" * 58)
    Dc, ACS, NCS = {}, {}, {}
    for d in range(len(CELLP2) + 1):
        v = np.zeros(nfft)
        v[:L] = (depth == d)
        nc = int(v.sum())
        f = np.fft.rfft(v)
        ac = np.fft.irfft(f * np.conj(f), nfft)[:L]
        ac = np.rint(ac).astype(np.int64)
        tot = int(2 * ac[1:].sum())
        if tot != nc * (nc - 1):
            say("  depth %d: pair count %d against n_c(n_c-1) = %d --"
                % (d, tot, nc * (nc - 1)))
            say("  the autocorrelation did not round to the exact counts;")
            say("  refusing to report a number this check does not back.")
            return 1
        e = float(2.0 * (ac[1:] * w).sum() / (nc * (nc - 1)))
        ACS[d], NCS[d] = ac, nc
        Dc[d] = e - EXACT_EALL
        say("  %-6d %-10d %-21.9f %.9f" % (d, nc, e, Dc[d]))
    say()
    say("  the pair counts reproduce n_c(n_c-1) exactly at all four")
    say("  cells, so the autocorrelation is the integer one.")
    say()

    # alpha is computed here, before T1, because T2's verdict needs to
    # know what alternative its tolerance was up against.  T3 reports it
    # below; nothing is scored here.
    al = []
    for d in range(3):
        gap = unc[d] - Dc[d]
        al.append(math.log(full[d] / gap) / math.log(ms[d])
                  if gap > 0 else None)
    good = [x for x in al if x is not None]

    say("T1  does the exact value lie between the uncorrected prediction")
    say("    and the fully corrected one?")
    say()
    say("  depth  lower        upper        exact        inside")
    t1 = True
    for d in range(3):
        lo, hi = unc[d] - full[d], unc[d]
        ins = lo <= Dc[d] <= hi
        if not ins:
            t1 = False
        say("  %-6d %-12.6f %-12.6f %-12.6f %s"
            % (d, lo, hi, Dc[d], "yes" if ins else "NO"))
    say()
    say("    T1: %s" % ("hold" if t1 else "REFUTED"))
    say()

    say("T2  is the 1/m corrected prediction right to %.0e?" % TOL_T2)
    say("    -- and can this tolerance SEE the alternative?  That second")
    say("    question is asked because the answer decides what a pass")
    say("    means, and it was not asked when the tolerance was set.")
    say()
    say("  depth  1/m corrected   exact          difference")
    t2 = True
    for d in range(3):
        dif = Dc[d] - corr[d]
        if abs(dif) > TOL_T2:
            t2 = False
        say("  %-6d %-15.6f %-14.6f %+.6f" % (d, corr[d], Dc[d], dif))
    say()
    # A third verdict.  T2 is not refuted -- the data do not contradict
    # 1/m.  It is not held either -- the data do not support it, because
    # alpha = 0.72 passes the same tolerance.  A test the alternative
    # also passes carries no evidence in either direction, and filing it
    # as held is worse than filing it as refuted: a reader counts a held
    # rule as support and there is none.
    aa0 = None
    if len(good) == 3:
        aa0 = sum(good) / 3.0
        gaps = [abs(full[d] / ms[d] - full[d] / ms[d] ** aa0)
                for d in range(3)]
        void = max(gaps) < TOL_T2
    else:
        void = False
    if void:
        say("    T2: VOID -- could not have discriminated.")
        say("    The alternative T3 measures, alpha = %.3f, differs from" % aa0)
        say("    1/m by %.1e to %.1e at these m, all under the %.0e"
            % (min(gaps), max(gaps), TOL_T2))
        say("    tolerance.  BOTH hypotheses pass it.  So T2 is recorded")
        say("    neither held nor refuted: the data did not contradict")
        say("    1/m and did not support it, and a rule filed as held")
        say("    would be counted as support that does not exist.")
    else:
        say("    T2: %s" % ("hold" if t2 else "REFUTED"))
    say()

    say("T3  is there ONE dilution law?  alpha spread cap %.2f" % TOL_T3)
    say("    alpha from")
    say("    exact = uncorrected - full / m^alpha")
    say()
    say("  depth  m     implied alpha")
    for d in range(3):
        if al[d] is None:
            say("  %-6d %-5d undefined -- the exact value is ABOVE the"
                % (d, ms[d]))
            say("         main term, so no positive deficit fits it")
        else:
            say("  %-6d %-5d %+.3f" % (d, ms[d], al[d]))
    t3 = len(good) == 3 and (max(good) - min(good)) <= TOL_T3
    say()
    say("    T3: %s%s" % ("hold" if t3 else "REFUTED",
                          "  (spread %.3f)" % (max(good) - min(good))
                          if len(good) == 3 else
                          "  (alpha undefined at %d of 3)" % (3 - len(good))))
    say()

    # POST HOC (after T1-T3 were scored; it rescores nothing).
    #
    # THE DILUTION EXPONENT SHOULD NOT BE FITTED AT ALL.  An independent
    # reading pointed out that the deficit is a sum of pieces with
    # DIFFERENT m-scalings, so no single power can describe it -- and
    # that this file's own three alphas already say so: they disagree by
    # 0.067 and NOT MONOTONICALLY in m, on a measurement with no
    # sampling in it, where that spread cannot be noise.
    #
    # And the deficit needs no model: it is the SAME autocorrelation
    # already run here, with T' in place of S_2.  So it is computed and
    # the exponent is dropped.
    say("POST HOC (after T1-T3 were scored; it rescores nothing).")
    say()
    say("  THE EXPONENT IS NOT A THING TO FIT.  alpha came out %s at"
        % ", ".join("%.3f" % x for x in al))
    say("  m = %s -- disagreeing by %.3f and NOT MONOTONE in m, on a"
        % (", ".join(str(ms[d]) for d in range(3)),
           max(good) - min(good)))
    say("  measurement with no sampling in it.  A single power law gives")
    say("  one number.  So m alone does not parameterise the deficit,")
    say("  and this file's own output says so before any argument about")
    say("  mechanism.")
    say()
    say("  AND THE THREE ARE ONE COMPARISON, not three.  alpha's")
    say("  leverage is log m:")
    lg = [math.log(ms[d]) for d in range(3)]
    say("    log m = %s" % ", ".join("%.3f" % x for x in lg))
    say("    depths 0 and 1 sit %.3f apart against %.3f from depth 0 to"
        % (abs(lg[0] - lg[1]), abs(lg[0] - lg[2])))
    say("    depth 2 -- the first pair carries %.1f times less leverage"
        % (abs(lg[0] - lg[2]) / abs(lg[0] - lg[1])))
    say("    and its two exponents are nearly forced to agree.  And the")
    say("    measurement is exact, so the %.3f spread is not an error"
        % (max(good) - min(good)))
    say("    bar -- it is a measured DISAGREEMENT between three numbers")
    say("    that are exactly unequal.  What makes them one law is the")
    say("    registered tolerance %.2f, not the data." % TOL_T3)
    say()
    say("  Computed instead, by the same convolution:")
    say("      Delta_c = C_T' - E_same,c[T'],   C_T' = K / 2C_2")
    say()
    Ap = np.ones(HMAX + 1)
    Tp = np.ones(HMAX + 1)
    for q in pr:
        q = int(q)
        if q == 2:
            continue
        if q in CELLP2:
            Ap[q::q] *= (q - 1.0) / (q - 2.0)
        else:
            Tp[q::q] *= (q - 1.0) / (q - 2.0)
    KEXF = 175.0 / 128.0
    CTp = KEXF / twin
    wA, wT, wAT = Ap[2 * kk], Tp[2 * kk], (Ap * Tp)[2 * kk]
    say("  C_T' = %.10f" % CTp)
    say()
    say("  depth  E[A'] exact     E[T'] exact     Delta_c        2C_2 Cov")
    Dl, Cv = {}, {}
    for d in range(len(CELLP2) + 1):
        ac, nc = ACS[d], NCS[d]
        den = float(nc * (nc - 1))
        eA = float(2.0 * (ac[1:] * wA).sum() / den)
        eT = float(2.0 * (ac[1:] * wT).sum() / den)
        eAT = float(2.0 * (ac[1:] * wAT).sum() / den)
        Dl[d] = CTp - eT
        Cv[d] = twin * (eAT - eA * eT)
        say("  %-6d %-15.9f %-15.9f %+.6e  %+.6e"
            % (d, eA, eT, Dl[d], Cv[d]))
    say()
    say("  E[A'] reproduces 403/240 = %.9f, 941/605 = %.9f,"
        % (403 / 240.0, 941 / 605.0))
    say("  31/15 = %.9f and 16/5 = %.9f, the exact rationals the"
        % (31 / 15.0, 16 / 5.0))
    say("  prediction is built from.  That is a check on the")
    say("  autocorrelation, not a result.")
    say()
    say("  With Delta_c computed the measurement decomposes with no")
    say("  model anywhere -- but into FOUR terms, not three.  The main")
    say("  term is K times the EXACT rational E_same,c[A'], and the")
    say("  autocorrelation's E[A'] is not that rational: it departs by")
    say("  O(1/ell) because the empirical class weights are not the")
    say("  idealised ones.  That departure is {#prop:scaleinv}'s OTHER")
    say("  error term, the equidistribution one, and nothing in this")
    say("  paper had measured it until now.")
    say()
    say("      D'_c = [K E_exact[A'] - 2]        main term, exact rational")
    say("             + K (E[A'] - E_exact[A'])  pattern-weight resolution")
    say("             - 2C_2 E[A'] Delta_c       divisor deficit")
    say("             + 2C_2 Cov(A', T')         covariance")
    say()
    EX = [F(403, 240), F(941, 605), F(31, 15), F(16, 5)]
    say("  depth  main term      pattern wt.    deficit        covariance")
    for d in range(len(CELLP2) + 1):
        ac, nc = ACS[d], NCS[d]
        den = float(nc * (nc - 1))
        eA = float(2.0 * (ac[1:] * wA).sum() / den)
        mt = KEXF * float(EX[d]) - 2.0
        pw = KEXF * (eA - float(EX[d]))
        dt = -twin * eA * Dl[d]
        say("  %-6d %-14.9f %+.4e    %+.4e    %+.4e"
            % (d, mt, pw, dt, Cv[d]))
        if d == 0:
            keep_pw = pw
    say("  sum    " + "  ".join(
        "%.9f" % (KEXF * float(EX[d]) - 2.0
                  + KEXF * (float(2.0 * (ACS[d][1:] * wA).sum()
                                  / float(NCS[d] * (NCS[d] - 1)))
                            - float(EX[d]))
                  - twin * float(2.0 * (ACS[d][1:] * wA).sum()
                                 / float(NCS[d] * (NCS[d] - 1))) * Dl[d]
                  + Cv[d]) for d in range(len(CELLP2) + 1)))
    say("  meas   " + "  ".join("%.9f" % Dc[d]
                                for d in range(len(CELLP2) + 1)))
    say()
    say("  The pattern-weight term is ZERO at depth %d -- one class, so"
        % len(CELLP2))
    say("  the empirical weight is the idealised one -- and nonzero")
    say("  wherever the cell unions more than one.  Same signature as")
    say("  the covariance, and the same reason.")
    say()
    say("  E[A'] therefore reproduces the exact rationals to about 2 ppm")
    say("  and EXACTLY at depth %d, the departure being that term."
        % len(CELLP2))
    say()

    say("  COVARIANCE.  The prediction assumed E[A'T'] = E[A'] E[T'].")
    say("  At depth %d that is EXACT -- A' is constant on every same-cell"
        % len(CELLP2))
    say("  pair of a single class -- which is why depth %d came in at"
        % len(CELLP2))
    say("  -9.9e-08.  At the other three it does not vanish.  And the")
    say("  first partition's calibration cell was a single class too, so")
    say("  the exponent was being fitted at cells carrying TWO terms")
    say("  that are structurally ABSENT where the deficit was")
    say("  calibrated.")
    say()

    # THE DIAGONAL-ONLY BASELINE, DERIVED RATHER THAN GUESSED.  An
    # earlier printing of this block used A(B)/(n_c-1), on a reading
    # supplied by a peer and adopted here WITHOUT DERIVING IT.  It is
    # wrong, and it made this file reverse a correct claim.  Derived:
    # a cell of m classes has m*ell(ell-1) same-class ordered pairs and
    # m(m-1)ell^2 different-class ones; under diagonal-only the second
    # kind carry exactly 1/d, so
    #     1/d - P_c(d) = m ell(ell-1) [1/d - P_single(d)] / (n_c(n_c-1)),
    # and summing against g gives Delta_full * (ell-1)/(n_c-1), which is
    # Delta_full/m to four digits.  The identity
    # m ell(ell-1) + m(m-1) ell^2 = n_c(n_c-1) is checked below.
    say("  COMPLETE CANCELLATION, DERIVED.  Diagonal class pairs alone")
    say("  predict Delta_full * m ell(ell-1)/(n_c(n_c-1)), which is")
    say("  Delta_full/m to four digits -- so alpha = 1 IS that model.")
    say("  An earlier printing of this block used A(B)/(n_c-1) instead,")
    say("  on a reading adopted here without deriving it; that is the")
    say("  prediction with EVERY class residual zero, diagonal included,")
    say("  which nobody proposed and which the single-class cell itself")
    say("  refutes.  Withdrawn.")
    say()
    say("  depth  m     pair identity   diag-only        measured Delta_c   residual")
    for d in range(len(CELLP2) + 1):
        nc = float(NCS[d])
        mm = float(ms[d])
        el = nc / mm
        idn = (mm * el * (el - 1) + mm * (mm - 1) * el * el) / (nc * (nc - 1))
        dg = Dl[len(CELLP2)] * mm * el * (el - 1) / (nc * (nc - 1))
        say("  %-6d %-5d %.9f     %.4e       %.4e        %+.1f%%"
            % (d, ms[d], idn, dg, Dl[d], 100.0 * (Dl[d] / dg - 1.0)))
    say()
    say("  MEASURED EXCEEDS DIAGONAL-ONLY at every multi-class cell.  So")
    say("  the off-diagonal class pairs contribute additional deficit --")
    say("  they do not cancel completely.  They do not over-cancel")
    say("  either; that reading came from the withdrawn baseline.")
    say("  AND THE EXCESS IS NOT A FUNCTION OF m EITHER: 41 per cent at")
    say("  m = 48 against 30 at m = 44 is a 38 per cent change in the")
    say("  excess for an 8 per cent change in m.  The three are quoted")
    say("  and not summarised, for the same reason the deficit is.")
    say()
    say("  THE BAND AVERAGE AT A FINITE BAND IS NOT 2.  The theorem is")
    say("  asymptotic.  Computed exactly over the whole band -- the same")
    say("  pair average with every even N as the cell, so the")
    say("  autocorrelation is analytic:")
    kb = np.arange(1, L, dtype=np.int64)
    wb = 2.0 * (L - kb)
    EallEx = float((wb * S2[2 * kb]).sum() / wb.sum())
    say("    E_all[S_2] over %d positions = %.12f" % (L, EallEx))
    say("    against the theorem's 2:  %+.4e" % (EallEx - 2.0))
    say()
    say("  THAT OFFSET IS LARGER THAN TWO OF THE FOUR TERMS ABOVE.  The")
    say("  registration fixed D'_c = E_same,c - 2 with the 2 exact, and")
    say("  T1, T2 and T3 are scored on that convention and stay scored")
    say("  on it.  But {#eq:Dc} defines D_c as E_same,c - E_all, so this")
    say("  file's D'_c and the first partition's D_c are two quantities")
    say("  under one name, differing by %.2e at a band of this size."
        % abs(EallEx - 2.0))
    say()
    say("  depth  D'_c (E_all = 2)   D'_c (E_all exact)  main term")
    for d in range(len(CELLP2) + 1):
        say("  %-6d %-18.9f %-19.9f %.9f"
            % (d, Dc[d], Dc[d] + (2.0 - EallEx),
               KEXF * float(EX[d]) - 2.0))
    say()
    say("  A FIFTH CELL, WHERE THE ANSWER IS KNOWN BEFORE THE RUN.  The")
    say("  band is all %d classes modulo 2Q', so it IS a cell of this" % 105)
    say("  partition -- the one with m = 105 -- and its main term is")
    say("  K E^ex_all[A_S] - 2 = (175/128)(256/175) - 2 = 0 exactly.  So")
    say("  the four-term identity applies to it with a main term known")
    say("  in advance, and decomposing its departure from 2 is the only")
    say("  check here whose answer is fixed before the computation.")
    say()
    wbA, wbT, wbAT = Ap[2 * kb], Tp[2 * kb], (Ap * Tp)[2 * kb]
    dnb = wb.sum()
    bA = float((wb * wbA).sum() / dnb)
    bT = float((wb * wbT).sum() / dnb)
    bAT = float((wb * wbAT).sum() / dnb)
    bmain = KEXF * float(F(256, 175)) - 2.0
    bpw = KEXF * (bA - float(F(256, 175)))
    bDl = CTp - bT
    bdt = -twin * bA * bDl
    bcv = twin * (bAT - bA * bT)
    say("    main term (known zero)  %+.6e" % bmain)
    say("    pattern weight          %+.6e" % bpw)
    say("    divisor deficit         %+.6e   (Delta = %.6e)" % (bdt, bDl))
    say("    covariance              %+.6e" % bcv)
    say("    sum                     %+.6e" % (bmain + bpw + bdt + bcv))
    say("    measured E_all,B - 2    %+.6e" % (EallEx - 2.0))
    say("    difference              %+.2e" % (bmain + bpw + bdt + bcv
                                               - (EallEx - 2.0)))
    say("  It closes.  So the band's departure from the limiting 2 is")
    say("  the same three terms measured at the cells, with the main")
    say("  term structurally absent -- the caution above is a")
    say("  confirmation and not an apology.")
    say()
    say("  Using the exact finite-band average removes the same")
    say("  %.4e at every cell -- half the departure from the main term" % abs(EallEx-2.0))
    say("  at depth 0 and under two per cent of it at depth 3, since a")
    say("  constant offset cannot take a fixed share off departures that")
    say("  span a factor thirty.  The band's own equidistribution error")
    say("  and the cell's partly cancel in a difference, which is the")
    say("  direction; the size is not uniform.  Reported and not")
    say("  folded in: the rules were registered on the other convention")
    say("  and rescoring them here would be choosing the convention")
    say("  after seeing which one flatters.")
    say()

    # NO POWER LAW -- and not the stronger "no function of m", which is
    # false: three points admit many functions and a monotone one fits
    # these comfortably.  What fails is a power.
    say("  NO POWER LAW IN m FITS THESE THREE.  Local exponents:")
    lex = []
    for i, j in ((0, 1), (1, 2), (0, 2)):
        e = (math.log(Dl[j] / Dl[i])
             / math.log(float(ms[i]) / float(ms[j])))
        lex.append(e)
        say("    m = %d -> %d :  %.3f" % (ms[i], ms[j], e))
    need = (float(ms[0]) / float(ms[1])) ** lex[1]
    say("  The deficits at m = %d and %d are in the ratio %.3f where a"
        % (ms[0], ms[1], Dl[1] / Dl[0]))
    say("  power fitted to the second pair requires %.3f." % need)
    say("  A monotone dependence on m is not excluded -- three points")
    say("  cannot exclude one -- but no single power reaches them.  That")
    say("  is a firmer reason to drop the exponent than 'no mechanism',")
    say("  and it warns against fitting anything else in m either.")
    say()
    # UNITS.  Delta_c lives in T-space; the covariance and pattern-weight
    # rows printed above are already carried through to D-space.  An
    # earlier printing compared them without converting and produced
    # percentages that reproduce under no convention.  Everything here
    # is T-space: the D-space rows are divided by 2 C_2 E_c[A_S].
    Qp2 = 1
    for q in CELLP2:
        Qp2 *= q
    ellv = (BAND[1] - BAND[0]) / (2.0 * Qp2)
    eA0 = float(2.0 * (ACS[0][1:] * wA).sum()
                / float(NCS[0] * (NCS[0] - 1)))
    sc0 = twin * eA0
    ncf, mmf = float(NCS[0]), float(ms[0])
    elf = ncf / mmf
    dgo = Dl[len(CELLP2)] * mmf * elf * (elf - 1) / (ncf * (ncf - 1))
    say("  HOW TIGHT IS THE PROPOSITION'S OWN BOUND?  The terms are")
    say("  quoted against the NUMBER 1/ell = %.3e and NOT against the"
        % (1.0 / ellv))
    say("  bound: the bound is that number times an unnamed constant and")
    say("  a divisor sum above three, which is the reading this file")
    say("  withdraws below.  So these are UPPER bounds on the fraction")
    say("  of the bound the measurement realises -- the bound is at")
    say("  least three times looser than they suggest, and the")
    say("  conclusion that it is loose survives a fortiori.")
    say("  At depth 0, all in T-space -- the covariance and")
    say("  pattern-weight rows above are divided by 2 C_2 E[A_S] = %.4f"
        % sc0)
    say("  to get here:")
    say("    divisor deficit      %.2f%% of the bound"
        % (100.0 * Dl[0] * ellv))
    say("    covariance           %.2f%%"
        % (100.0 * abs(Cv[0]) / sc0 * ellv))
    say("    pattern weight       %.2f%%"
        % (100.0 * abs(keep_pw) / sc0 * ellv))
    say("    diagonal-only        %.2f%%  (what complete cancellation of"
        % (100.0 * dgo * ellv))
    say("                         the class-wise errors would give, and")
    say("                         the proof declines to claim it)")
    say("  So the measured deficit EXCEEDS the diagonal-only figure --")
    say("  the off-diagonal pairs add deficit rather than cancelling it")
    say("  -- and all three measured terms are a few per cent of the")
    say("  bound the proposition proves.")
    say()

    # POST HOC (after T1, T2 and T3 were scored; it rescores nothing).    # POST HOC (after T1, T2 and T3 were scored; it rescores nothing).
    say("POST HOC (after T1-T3 were scored; it rescores nothing).")
    say()
    say("  T2 HELD AND T3 SAYS alpha IS NOT 1.  Those are in tension and")
    say("  the tension is T2's tolerance, which I set at 1e-4 without")
    say("  checking what it could discriminate:")
    say()
    say("    depth  1/m correction  alpha-fit correction  difference")
    aa = sum(x for x in al if x is not None) / max(1, len(good))
    for d in range(3):
        c1 = full[d] / ms[d]
        c2 = full[d] / ms[d] ** aa
        say("    %-6d %-15.6f %-21.6f %.6f" % (d, c1, c2, abs(c1 - c2)))
    say()
    say("  Every difference is under the 1e-4 cap, so T2 COULD NOT HAVE")
    say("  REFUSED alpha = %.3f.  It held against an alternative it was" % aa)
    say("  not sharp enough to see.  A rule whose tolerance is not")
    say("  calibrated against the alternative it is meant to exclude")
    say("  passes for a reason unrelated to the claim.")
    say()
    say("  And the three residuals against the 1/m prediction are")
    say("  %s -- all the same sign."
        % ", ".join("%+.6f" % (Dc[d] - corr[d]) for d in range(3)))
    say("  A systematic sign at three of three is what T3 reads and T2")
    say("  cannot: the dilution is a power of m near %.2f, not 1." % aa)
    say()
    say("  WHAT alpha = %.3f MEANS.  The 1 came from a mechanism -- the" % aa)
    say("  deficit carried by the m diagonal class pairs and nothing")
    say("  else.  %.2f has no mechanism proposed for it, so the reading" % aa)
    say("  is that the off-diagonal class pairs do NOT cancel completely;")
    say("  they carry a residual share of the deficit.  That is the")
    say("  cancellation {#prop:scaleinv} declines to claim, measured")
    say("  rather than assumed, and it is INCOMPLETE.")
    say()
    say("  CONSEQUENCE FOR THE FIRST PARTITION, which the registration")
    say("  named in advance for the T2-refuted case.  T2 was not refused,")
    say("  so that clause does not fire literally; but its reason has")
    say("  come true and the honest thing is to report the effect either")
    say("  way.  At depth 4 of the first partition (m = 34):")
    say()
    say("    no correction        residual +1.16 sampling errors")
    say("    alpha = 1            residual +0.34")
    say("    alpha = %.3f          residual -1.04" % aa)
    say()
    say("  All three inside two sampling errors, so the first partition")
    say("  does not decide the exponent and is not decided by it.  The")
    say("  paper should carry the exponent as measured HERE, on a cell")
    say("  shape where m spans 48 to 1, and not lean on depth 4.")
    say()

    say("DEPTH 3 -- consistency check, not a score, as the first seal")
    say("recorded: the corrected prediction and the exact value are the")
    say("same finite sum.")
    say("    sealed corrected prediction  %.9f" % corr[3])
    say("    exact                        %.9f" % Dc[3])
    say("    difference                   %.3e" % (corr[3] - Dc[3]))
    say()

    say("=" * 70)
    say("T1 %s   T2 %s   T3 %s"
        % ("hold" if t1 else "REFUTED",
           "VOID" if void else ("hold" if t2 else "REFUTED"),
           "hold" if t3 else "REFUTED"))
    if t1:
        say("The divisor deficit is the mechanism at multi-class cells")
        say("too, and the main-term computation transfers to a cell shape")
        say("it was not built on.")
    else:
        say("The divisor deficit does NOT account for the departure from")
        say("the main term at multi-class cells.  The transfer of the")
        say("main-term computation to a new cell shape is not")
        say("established, and the paper's claim narrows to the")
        say("(3,5,7,11,13) partition until it is.")

    head = [
        "STATISTIC: E_same,c[S_2] over each cell of the second shape,",
        "           EXACTLY, by autocorrelation of the cell indicator --",
        "           no sampling, no seed, no standard error.",
        "           D'_c = that minus the exact band average 2.",
        "DENOM: n_c(n_c-1) ordered pairs, checked against the",
        "       autocorrelation's own total before any number is",
        "       reported.",
        "NULL: the rules were registered in lab_secondcell2_predict.txt,",
        "      which carries an enumerated null on whether T1's",
        "      intervals are cell-specific. The statistic here is exact,",
        "      so the tolerances are absolute and no distribution",
        "      enters.",
        "FIELD: band (2e6, 4e6]; even N; cells by depth over 3,5,7; S_2",
        "       sieved to 2e6 over every prime with 2 C_2 to 4e6; rFFT",
        "       at length 2^21.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    return 0 if t1 else 1


if __name__ == "__main__":
    sys.exit(main())

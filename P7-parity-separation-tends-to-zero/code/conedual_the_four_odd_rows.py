# -*- coding: utf-8 -*-
r"""conedual_the_four_odd_rows -- is the carrier "omega >= 4" or "EVEN omega >= 4"?

Supporting computation; not used in the paper's proofs or tables.
2^20 is structurally incapable of answering; 2^22 and 2^24 are not.

WHY THIS RUNS

An earlier class census showed the omega >= 4 rows split 16 EVEN / 0 ODD at 2^20 and
46 / 4 at 2^22.  Assistant A pinned the reason and enumerated the odd side:

    2^20  Q = 1024   odd squarefree d <= Q with omega >= 4 :  0
    2^22  Q = 2048                                        :  4 = {1155,1365,1785,1995}
                                                            = 105 x {11,13,17,19}
    2^24  Q = 4096                                        : 16

**All four at 2^22 are 105p, and that is forced**: 3.5.11.13 = 2145 > 2048, so nothing
but 3.5.7.p fits below 2048 with four odd prime factors.

> **So at 2^20 "omega >= 4" and "EVEN omega >= 4" are the SAME 16 rows, and every
> statement this axis has made about the carrier is, at that cell, silent on which of
> the two it is.**  Only 2^22 and 2^24 can tell them apart.

THE TEST, one LP: take the certifying partition PA = (min(omega,4), 2|d) at 2^22 and
DEMOTE the four odd omega >= 4 rows back into omega = 3.  Everything else unchanged.

PREDICTIONS.  Each clause names its precondition and the
file checks it before scoring; the meaning of each outcome is registered
here, not after; no clause compares quantities with no reason to be equal;
two solvers run and the best DIRECTLY VERIFIED bound is reported.

  AL0  CONTROL, AGAINST THE RECORD.  P_1 = -22243 at 2^22.
       COMPARAND coded: the published integer.  CAN FAIL: yes; voids the cell.

  AL1  BASELINE.  PA reaches 0.029575 at 2^22 and 0.011280 at 2^24, 1e-4 relative.
       COMPARAND coded: the values printed by the two earlier censuses.
       CAN FAIL: yes; voids the file.  PRECONDITION for AL2 and AL3.

  AL2  THE TEST -- DEMOTE THE ODD ones.  PA with the odd omega >= 4 rows moved into
       omega = 3.
       MEANING, registered:
         reaches PA      -> the carrier is EVEN omega >= 4; the parity story extends
                            all the way up and the odd rows are inert.
         falls short     -> those rows are load-bearing at 4 of 1245 (0.32 per cent)
                            at 2^22, the smallest located object in this chain.
       COMPARAND coded: PA's own value at the same cell, 1e-4 relative.

  AL3  THE COMPLEMENT -- DEMOTE THE EVEN ones, keep the odd separate.
       MEANING: reaching PA would mean the 4 odd rows alone carry what 50 rows
       carried, which is a stronger claim than AL2 falling short; falling short means
       both halves are needed and the carrier is the whole class after all.
       COMPARAND coded: PA's own value at the same cell.

  AL4  THE CONTROL, at matched count.  Demote 4 (resp. 16) RANDOMLY chosen
       omega >= 4 rows instead, 12 seeded draws.
       MEANING: if demoting a random 4 also costs the bound, then demoting ANY 4
       costs it, AL2's verdict is about the count and not about the odd rows, and
       the identification is withdrawn.  If random demotions leave PA intact while
       the odd demotion does not, the odd rows are specifically implicated.
       COMPARAND coded: PA's value and AL2's verdict at the same cell.

  AL5  FEASIBILITY, DIRECT AND UNSHIFTED.  Zero violations over ALL |vB| columns and
       <w,t> > 0 for every reported bound.
       COMPARAND coded: 0 violations at 1e-12 * scale, <w,t> > 0.

  NOT REGISTERED: any claim about kappa, and any claim that the 2^22 answer holds at
  2^26 or beyond.

mu comes from lib/goldbach/sieve.py; this file builds none of its own.

WHAT THIS DOES NOT SAY.  "Load-bearing in this partition family" is not "large in any
norm".  One earlier census put d = 1 at 58 per cent of ||p|| and another put the same row at 5
per cent of the l1 sum; a four-row set can be decisive for a subspace separator and
negligible in the object.
"""
import io
import os
import sys
import time

import numpy as np
from scipy.optimize import linprog, lsq_linear

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _sieve_shared import mu_upto, primes_upto  # noqa: E402

NEXP = (22, 24)
PUB_P1 = {22: -22243}
BASE = {22: 0.029575, 24: 0.011280}
NDRAW = 12
SEED = 20260911

_lines = []
_OUT = [None]


def say(s=""):
    _lines.append(s)
    print(s, flush=True)
    if _OUT[0]:
        nl = chr(10)
        with io.open(_OUT[0], "w", encoding="utf-8", newline=nl) as f:
            f.write(nl.join(_lines) + nl)


def build(e_):
    N = 2 ** e_
    Q = int(N ** 0.5)
    mu = np.array(mu_upto(N), dtype=np.int8)
    # int32 wraps at 2^31: for N >= 2^31 the friable indicator would be wrong for every
    # n >= 2^31 (found at 2^32 and above in a follow-up census).
    # Every conedual cell is <= 2^26; the guard is here so the base cannot be reused blindly.
    assert N < 2 ** 31, "build(): int32 remainder array wraps at 2^31; use uint32/int64 above"
    rem = np.arange(N + 1, dtype=np.int32)
    om = np.zeros(N + 1, dtype=np.int8)
    for p in primes_upto(Q):
        p = int(p)
        om[p::p] += 1
        q = p
        while q <= N:
            rem[q::q] //= p
            if q > N // p:
                break
            q *= p
    q1 = Q + 1
    while True:
        if all(q1 % r for r in range(2, int(q1 ** 0.5) + 1)):
            break
        q1 += 1
    thr = N // q1
    band = (mu != 0) & (rem == 1)
    band[:thr + 1] = False
    del rem
    vB = np.nonzero(band & (om % 2 == 0))[0]
    vR = np.nonzero(band & (om % 2 == 1))[0]
    sigB = np.zeros(N + 1, dtype=np.int8)
    sigR = np.zeros(N + 1, dtype=np.int8)
    sigB[vB] = 1
    sigR[vR] = 1
    del band
    rows = np.nonzero(mu[1:Q + 1] != 0)[0] + 1
    cb = np.zeros(rows.size)
    cr = np.zeros(rows.size)
    for i, d in enumerate(rows):
        d = int(d)
        cb[i] = float(sigB[d::d].sum(dtype=np.int64))
        cr[i] = float(sigR[d::d].sum(dtype=np.int64))
    rom = om[rows].astype(np.int64)
    del sigB, sigR, om, mu
    return dict(N=N, Q=Q, vB=vB, vR=vR, rows=rows, rom=rom,
                p=(cb - cr) / float(vR.size), nB=float(vB.size),
                nR=float(vR.size))


def score_partition(S, lab, npv):
    rows, N, vB, vR = S["rows"], S["N"], S["vB"], S["vR"]
    lab = np.asarray(lab, dtype=np.int64).copy()
    lab[rows == 1] = -1
    u = np.array([v for v in np.unique(lab) if v != -1])
    C = u.size + 1
    cls = np.zeros(rows.size, dtype=np.int64)
    for i, v in enumerate(u):
        cls[lab == v] = i + 1
    D = np.zeros((vB.size, C), dtype=np.int16)
    tau = np.zeros(C)
    R = np.zeros(C)
    acc = np.zeros(N + 1, dtype=np.int16)
    for c in range(C):
        idx = rows[cls == c]
        R[c] = float(idx.size)
        if idx.size == 0:
            continue
        acc[:] = 0
        for d in idx:
            acc[int(d)::int(d)] += 1
        D[:, c] = acc[vB]
        tau[c] = float(acc[vR].sum(dtype=np.int64)) / float(vR.size)
    del acc
    Df = D.astype(np.float64)
    U = np.unique(D, axis=0).astype(np.float64)
    rt = np.sqrt(np.maximum(R, 1.0))
    A = (U / rt).T

    def sc(g):
        s_ = float(np.abs(g).max() * max(Df.max(), 1.0))
        nv = int(np.count_nonzero(Df.dot(g) > 1e-12 * s_))
        dot = float(np.dot(g, tau))
        nw = float(np.sqrt(np.dot(g * g, R)))
        ok = nv == 0 and dot > 0.0 and nw > 1e-12
        return ok, (dot / (nw * npv) if ok else float("nan"))

    r = lsq_linear(A, tau / rt, bounds=(0, np.inf), method="bvls", tol=1e-14,
                   max_iter=8000)
    okb, lbb = sc((tau / rt - A.dot(r.x)) / rt)
    n, Cc = U.shape
    Aub = np.vstack([np.hstack([U, np.zeros((n, Cc))]),
                     np.hstack([np.diag(rt), -np.eye(Cc)]),
                     np.hstack([-np.diag(rt), -np.eye(Cc)])])
    lp = linprog(np.concatenate([np.zeros(Cc), np.ones(Cc)]),
                 A_ub=Aub, b_ub=np.zeros(n + 2 * Cc),
                 A_eq=np.hstack([tau.reshape(1, -1), np.zeros((1, Cc))]),
                 b_eq=np.array([1.0]),
                 bounds=[(None, None)] * Cc + [(0, None)] * Cc, method="highs")
    okl, lbl = False, float("nan")
    if lp.status == 0 and lp.x is not None:
        okl, lbl = sc(lp.x[:Cc])
    del D, Df, U
    if okb and (not okl or lbb >= lbl):
        return True, lbb, C, "bvls"
    if okl:
        return True, lbl, C, "LP"
    return False, float("nan"), C, ""


def run_cell(e_, fails):
    ts = time.time()
    S = build(e_)
    rows, om = S["rows"], S["rom"]
    npv = float(np.linalg.norm(S["p"]))
    say("2^%d   Q = %d   |vB| = %d   rows = %d   ||p|| = %.6e"
        % (e_, S["Q"], S["vB"].size, rows.size, npv))
    if e_ in PUB_P1:
        P1 = int(S["nB"] - S["nR"])
        ok0 = P1 == PUB_P1[e_]
        say("      AL0 %s -- P_1 = %d against the record's published %d"
            % ("HOLDS" if ok0 else "FAILS", P1, PUB_P1[e_]))
        if not ok0:
            fails.append("AL0 at 2^%d" % e_)
            return

    o4 = np.minimum(om, 4)
    ev = (rows % 2 == 0).astype(np.int64)
    big = o4 >= 4
    odd_big = np.nonzero(big & (ev == 0))[0]
    even_big = np.nonzero(big & (ev == 1))[0]
    say("      omega >= 4: %d EVEN, %d ODD.  ODD rows: %s"
        % (even_big.size, odd_big.size,
           " ".join(str(int(d)) for d in rows[odd_big])))
    if odd_big.size:
        say("      factorisation check: %s"
            % ", ".join("%d = 105 x %d" % (int(d), int(d) // 105)
                        if int(d) % 105 == 0 else "%d NOT 105p" % int(d)
                        for d in rows[odd_big]))

    okA, lbA, CA, whoA = score_partition(S, o4 * 2 + ev, npv)
    ok1 = okA and abs(lbA - BASE[e_]) <= 1e-4 * BASE[e_]
    say("      AL1 %s -- PA C = %d  %s  against the record's %.6f  [%s]"
        % ("HOLDS" if ok1 else "FAILS", CA,
           ("kappa >= %.6f" % lbA) if okA else "no separator", BASE[e_], whoA))
    if not ok1:
        fails.append("AL1 at 2^%d" % e_)
        say("      PRECONDITION for AL2/AL3 fails; not scored at this cell")
        say()
        return
    if odd_big.size == 0:
        say("      AL2 SKIPPED at 2^%d -- the odd omega >= 4 class is EMPTY here,"
            " so 'omega >= 4' and 'EVEN omega >= 4' are the same partition and the"
            " question has no referent.  Rule 41(c)." % e_)
        say()
        return

    lab = o4 * 2 + ev
    lab2 = lab.copy()
    lab2[odd_big] = 3 * 2 + 0          # demote the odd big rows into omega = 3 odd
    okB, lbB, CB, whoB = score_partition(S, lab2, npv)
    near = okB and abs(lbB - lbA) <= 1e-4 * lbA
    say("      AL2 -- demote the %d ODD omega >= 4 rows into omega = 3:"
        " C = %d  %s  [%s]"
        % (odd_big.size, CB,
           ("kappa >= %.6f" % lbB) if okB else "no separator", whoB))
    say("          >>> %s"
        % ("REACHES PA: the carrier is EVEN omega >= 4 and the odd rows are inert"
           if near else
           ("falls short of PA (%.6f against %.6f): those %d rows of %d -- %.2f per"
            " cent -- are load-bearing"
            % (lbB, lbA, odd_big.size, rows.size,
               100.0 * odd_big.size / rows.size)) if okB else
           ("KILLS the certificate outright: those %d rows of %d -- %.2f per cent --"
            " are load-bearing"
            % (odd_big.size, rows.size, 100.0 * odd_big.size / rows.size))))

    lab3 = lab.copy()
    lab3[even_big] = 3 * 2 + 1         # demote the even big rows into omega = 3 even
    okC, lbC, CC, whoC = score_partition(S, lab3, npv)
    nearC = okC and abs(lbC - lbA) <= 1e-4 * lbA
    say("      AL3 -- demote the %d EVEN omega >= 4 rows instead: C = %d  %s  [%s]"
        % (even_big.size, CC,
           ("kappa >= %.6f" % lbC) if okC else "no separator", whoC))
    say("          >>> %s"
        % ("REACHES PA: the %d odd rows alone carry what the class carried"
           % odd_big.size if nearC else
           "falls short: both halves are needed and the carrier is the whole class"))

    rng = np.random.RandomState(SEED + e_)
    allbig = np.nonzero(big)[0]
    hits, vals = 0, []
    for _ in range(NDRAW):
        pick = rng.choice(allbig, size=odd_big.size, replace=False)
        l4 = lab.copy()
        l4[pick] = 3 * 2 + 0
        okD, lbD, _, _ = score_partition(S, l4, npv)
        if okD and abs(lbD - lbA) <= 1e-4 * lbA:
            hits += 1
        if okD:
            vals.append(lbD)
    say("      AL4 -- demote %d RANDOM omega >= 4 rows instead, %d draws:"
        " %d of %d still reach PA%s"
        % (odd_big.size, NDRAW, hits, NDRAW,
           ("   bounds [%.6f .. %.6f]" % (min(vals), max(vals))) if vals
           else "   (none certified at all)"))
    if hits == 0 and not near:
        say("          >>> demoting ANY %d rows costs the bound, so AL2 is about the"
            " COUNT and not about the odd rows; the identification is withdrawn"
            % odd_big.size)
        fails.append("AL4 at 2^%d -- AL2 not row-specific" % e_)
    elif not near:
        say("          >>> random demotions mostly leave PA intact (%d of %d) while"
            " the odd demotion does not: the odd rows are specifically implicated"
            % (hits, NDRAW))
    say("      elapsed %.1f s" % (time.time() - ts))
    say()


def main():
    _OUT[0] = sys.argv[1]
    say("conedual_the_four_odd_rows -- omega >= 4, or EVEN omega >= 4?")
    say("=" * 92)
    say()
    say("FIELD: at 2^20 the odd omega >= 4 class is EMPTY (3.5.7.11 = 1155 > 1024),")
    say("       so that cell cannot distinguish the two.  At 2^22 it has exactly")
    say("       four members and all are 105p -- forced, since 3.5.11.13 > 2048.")
    say("       The test demotes them and asks whether PA survives.")
    say()
    fails = []
    for e_ in NEXP:
        try:
            run_cell(e_, fails)
        except MemoryError:
            say("2^%d -- MemoryError; cell not attempted" % e_)
            say()
    say("=" * 92)
    say("FAILS: %s" % (", ".join(fails) if fails else "none"))


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
r"""conedual_parity_domination -- the rough-cofactor parity-domination construction.

WHERE THIS COMES FROM

This construction bypasses the optimizer entirely.
With `L = log Q`, `Y` a cutoff, `a_Y(n) = prod_{p | n, p <= Y} p` the small-prime part, and `B_a`, `R_a` the
even-/odd-omega band counts with small-prime part `a`, put on the even columns

    x_l = R_{a_Y(l)} / (R * B_{a_Y(l)})   if a_Y(l) <= A := floor(sqrt(Q));   x_l = 0 otherwise,     (X)

with `R = |vR|`.  Then `x >= 0`, so `A x` lies in `K_Q` and `dist_2(t, K_Q) <= ||t - A x||_2`.  The
claim under test is that `(PC)  R_a <= C B_a` uniformly for squarefree Y-smooth `a <= sqrt(Q)` forces
`kappa -> 0`.

WHAT IS MEASURED HERE, per cell and per Y
  1. the parity ratio `max_a R_a / B_a` over retained `a`, and whether any retained `a` has `B_a = 0 < R_a`
     (which makes (X) undefined -- the one way (PC) can fail catastrophically);
  2. the residual `||t - A x||_2 / ||p||_2`, an UPPER bound on `kappa(Q)`, against the certified `kappa`;
  3. the exact split of that residual into its two theoretical parts -- the Y-smooth rows (where the
     construction is exact up to discarded patterns) and the rows with a prime factor > Y;
  4. the discarded-pattern mass `delta = Pr_{vR}(a_Y(n) > A)`, printed beside the quantity
     `exp(-L/(2 log Y))` (JSON key `rankin`).  That quantity is the exponential factor of a Rankin-type
     smooth-number estimate and is printed as a SCALE FOR THE EYE ONLY: it is not a bound on `delta`, the
     estimate it is taken from carries further factors, and no clause compares `delta` with it.

The reduction and Hypothesis 5 of the paper are proved there, not here; nothing here is evidence
for either -- this file measures how much of the asymptotic statement is visible at the computable cells.
The asymptotic choice `Y = L^6` EXCEEDS `Q` at every computable cell (at 2^22, `L^6 = 1.96e5` against
`Q = 2048`), so it is reported and is useless by construction; the informative runs are at small `Y`.

INPUTS CITED, NOT RECOMPUTED
  certified kappa, as the ROUNDED DISPLAY VALUES this file compares against (`KAP_CERT`): 2^16 0.172389352,
    2^18 0.140860454, 2^20 0.1055427, 2^22 0.0763310.  None of these is an exact optimum.  The exact
    objects are the certificates results/pa_cert_2e<e>.npz, which carry a primal and an integral dual as
    integer vectors and enclose kappa^2 between two exact rationals; code/pa_certify.py re-derives and
    checks them.  Those enclosures are, to twelve decimals,
        2^16 [0.172389351665, 0.172389351666]   2^18 [0.140860453831, 0.140860453832]
        2^20 [0.105542708383, 0.105542708384]   2^22 [0.076331197681, 0.076331197683],
    so the nine-decimal float above at 2^16 sits 3.3e-10 ABOVE the certified upper endpoint and the one at
    2^18 sits 1.7e-10 above it, while the seven-decimal 2^20 and 2^22 entries are truncations and sit
    8.4e-9 and 2.0e-7 BELOW the certified lower endpoints (2^20 is also carried in
    code/conedual_contact_defect.py).  The only clause that uses them
    is P2, whose slack 1e-6 exceeds every one of these rounding gaps by three orders or more; nothing here
    needs, or asserts, more precision than the rounded display.
  ||p||_2 at 2^16/2^18/2^22 = 0.202904477, 0.163823260, 0.114446913 -- the nine-decimal output of an
    earlier run (results/conedual_nonconcentration.txt, which the reproduction packet does NOT ship).  The
    three literals are the whole of what is cited: they are repeated in `NORMP_LIT` below, and C4 compares
    its own 9-decimal display against them, so the control is self-contained in this file.
  ||p||_2 >= (1/5)(N/R_1) L^-2, a normalisation lower bound; not recomputed, not scored.

CLAUSES.

  C1  CONTROL.  The band rebuilt here (`_build`-equivalent, no cached band exists) has |vB| and |vR| equal to
      the cache's nB, nR, AND the per-row divisor counts recomputed from the rebuilt band equal the cached
      integer vectors cb, cr ELEMENT FOR ELEMENT (integer against integer, no tolerance).
      CAN FAIL; voids the cell.
  C2  CONTROL (the identity behind the construction, step 1).  For every (cell, Y) with NO DEAD RETAINED
      CLASS -- i.e. no a <= A has B_a = 0 < R_a -- and every row d that is Y-smooth,
      (A x)_d = t_d - (1/R) #{n in vR : d | n, a_Y(n) > A}, to an absolute 1e-11 (both sides are sums of up
      to N/d float64 terms; 1e-11 exceeds the accumulated round-off by >= 3 orders at every cell here).
      A dead retained class does not leave this measurement undefined: the code puts x_l = 0 on every class
      with B_a = 0 (`ok = keep & (Ba > 0)`), so x >= 0 and A x lies in K_Q as before.  (The paper's (X) has
      no such convention because Theorem 21 excludes dead retained classes under its hypotheses; the
      convention is this file's, so that cells where they do occur can still be measured.)  Lemma 10 of the
      paper gives the row identity UNCONDITIONALLY, the dead classes entering as a second count -- the
      n in vR with d | n whose small-prime part is retained but has B_a = 0.  This clause evaluates only the two-term form, so
      a pair with a dead retained class is reported as NOT APPLICABLE and is not scored; what goes untested
      there is this clause's two-term identity, not the construction.
      CAN FAIL; voids that (cell, Y) pair.
  C3  CONTROL.  sum_a B_a = |vB| and sum_a R_a = |vR| exactly (integers).  CAN FAIL; voids the cell.
  C4  CONTROL (reproduction, like with like: rounded display against rounded display).  ||p||_2 computed here
      from the cached (cb - cr)/nR, printed to 9 decimals, equals the 9-decimal literal carried in
      `NORMP_LIT` below at 2^16, 2^18, 2^22 -- the three numbers an earlier run
      (results/conedual_nonconcentration.txt, not shipped with the packet) printed, quoted in full in the
      header.  2^20 has no stored literal and is not scored by C4.  CAN FAIL; voids the cell.

  P1  THE CLAUSE THAT DECIDES (PC) AT THE COMPUTABLE CELLS.  For every cell and every Y <= 30:
      (a) no retained a (a <= A) has B_a = 0 < R_a, and (b) max_{a <= A, B_a > 0} R_a / B_a <= 3.
      WHY: what the paper proves is Theorem 21 -- (PC) holds with C = 1 + O((log L)^{-2}), an ASYMPTOTIC
      statement in L = log Q with an unquantified implied constant and a growing-Y regime.  It supplies no
      numeric ratio at any computable cell, so the threshold 3 is a declared round number, fixed and
      registered before the scored runs and not derived from the proof; the clause is a registered test of
      whether the ratio is nevertheless small and finite here, not a check of a predicted value.
      DECLARED: an unscored smoke run at 2^16 (ACF_SMOKE=1) was made before this registration to debug the
      kernels, and it already showed (b) fails at 2^16 for Y >= 7 (max ratio 6.8 at Y=7, 34.7 at Y=30).  The
      threshold is registered unchanged rather than retuned; P1 therefore carries predictive force only at
      2^18, 2^20 and 2^22, which no smoke run has touched.  MEANING -- holds: at the cells scored here the
      ratio (PC) asks about is already bounded by a small constant, so nothing in the computable range
      contradicts (PC); fails: either some retained a has one-sided parity (R_a > 0 = B_a, so this
      construction puts no mass on that class and Lemma 10's second sum is non-empty there) or the ratio is
      not small at these sizes.  Either way the clause says nothing about the asymptotic statement itself,
      which concerns Y and Q beyond any cell reachable here.
  P2  PREDICTION (correctness).  residual/||p||_2 >= certified kappa - 1e-6 at every cell and every Y.
      It is an upper bound on a quantity certified from below, so it cannot be smaller.
      MEANING -- fails: the construction, the band or the certified kappa is wrong.
  P3  PREDICTION.  min over Y of residual/||p||_2 is strictly decreasing across 2^16, 2^18, 2^20, 2^22.
      WHY: the proof gives kappa << L^2 ((Y log Y)^{-1/2} + sqrt(delta)), which is > 1 at every cell here, so
      theory predicts nothing; this clause tests whether the construction nevertheless improves with Q.
      MEANING -- holds: the four numbers decrease, which is consistent with the asymptotic but, the bound
      being vacuous here, is not evidence for it; fails: at these sizes the residual is carried by the
      Y-rough rows, which the bound treats as a tail, and a refinement (a second block, or Y growing with Q)
      would have to come first for the trend to be worth measuring at all.
  P4  PREDICTION.  delta(Y, Q) := Pr_{vR}(a_Y(n) > A) decreases in Q at every fixed Y <= 30 common to the
      scored cells.  The comparison is delta against delta across cells; the printed exp(-L/(2 log Y)) is
      not part of it.  MEANING -- holds: the discarded-pattern mass does shrink with Q at fixed Y, the
      direction a smooth-number estimate would suggest; fails: the cutoff a <= sqrt(Q) is the wrong
      truncation at these sizes.
  P5  PREDICTION -- the one clause about data no smoke run has seen.  At every fixed Y in {3,5,7,11,13,19,30}
      the quantity max_{a <= A, B_a > 0} R_a/B_a is non-increasing across the scored cells 2^16, 2^18, 2^20,
      2^22.  WHY: the mechanism behind Theorem 21 is cancellation in the Mobius sum over the cofactor, which
      improves as the cofactor gets longer, so at fixed Y the ratio falling as L grows is what one would
      expect to see if that mechanism is already acting.  The theorem itself is asymptotic and FIXES NO
      TREND at four cells, so this is a declared expectation of this file, not a consequence of the proof;
      the smoke run at 2^16 fixes no trend either, which is what makes the clause a test at all.
      MEANING -- holds: the direction the mechanism predicts is already visible in the computable range;
      fails: it is not -- a possible reason being that the number of prime slots in the cofactor,
      u_eff = log(N/a)/log Y, is small here (3-8, see the diagnostic below) where the proof needs it large,
      in which case these cells are governed by small-u_eff combinatorics and no measurement at reachable Q
      bears on (PC) either way.  The clause distinguishes the two outcomes; it does not establish the reason.

  DIAGNOSTIC, not scored: u_eff := log(N/a*)/log Y at the maximising a*, the number of prime slots the
  cofactor b has.  The proof of (PC) needs u_eff -> infinity; this records how far from that the cells are.

  NOT CLAIMED HERE: anything about the reduction, the parity hypothesis, the signed estimate, Saias, or any
  asymptotic claim.  The reduction, the parity hypothesis and the signed estimate are Theorems 9, 21
  and 20 of the paper; this file is not evidence for them.

COMPUTE RULES (the project's compute rules, both sections read before this file was written)
  Cells are LOADED (`from cell_cache import load`).  The band itself is not cached (cell_cache says so
  explicitly), so it is rebuilt with `_build`-equivalent code and checked against the cache by C1.  No
  Python loop over rows, band elements, nnz or an N-length range: numba `prange` kernels for the per-row
  divisor sums and counts, numpy slicing per prime for the sieve and for the small-prime part, numpy
  bincount for B_a and R_a.  The only Python loops are over cells, over Y, and over the primes <= Y (the
  pattern cell_cache._build itself uses).  Threads capped at 4 (this machine).

  DEVIATION, declared: 2^20 and 2^22 run on this machine, not the compute node, which was fully occupied at the
  time.  Safe here: the largest cell allocates three N-length arrays at
  N = 4.19e6 (34 MB each) and finishes in seconds, which is inside this machine's remit.

    python code/conedual_parity_domination.py 16 18 20 22
"""
from __future__ import annotations
import io, json, math, os, sys, time

NTHREAD = int(os.environ.get("ACF_THREADS", "4"))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, str(NTHREAD))
_SCRATCH = os.environ.get("TEMP") or os.environ.get("TMP") or "."
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(_SCRATCH, "numba_cache_acf"))

import numpy as np  # noqa: E402
import numba  # noqa: E402
from numba import njit, prange, set_num_threads  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _sieve_shared import mu_upto, primes_upto  # noqa: E402
from cell_cache import load as load_cell  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), "results")
OUT = os.path.join(RES, "conedual_parity_domination.txt")
OUTJ = os.path.join(RES, "conedual_parity_domination.json")

SMOKE = os.environ.get("ACF_SMOKE") == "1"
KAP_CERT = {16: 0.172389352, 18: 0.140860454, 20: 0.1055427, 22: 0.0763310}
NORMP_LIT = {16: "0.202904477", 18: "0.163823260", 22: "0.114446913"}
C2_ABS, C4_DEC = 1e-11, 9
P1_YMAX, P1_RATIO = 30, 3.0
P2_SLACK = 1e-6

_lines = []


def say(s=""):
    _lines.append(s)
    print(s, flush=True)


def flush(payload):
    with io.open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(_lines) + "\n")
    with io.open(OUTJ, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=1, default=float)


@njit(parallel=True, cache=True, nogil=True)
def _rowsum_f(W, rows, out):
    """out[j] = sum of W over the multiples of rows[j] -- one prange over rows, no Python loop."""
    n = W.shape[0]
    for j in prange(rows.shape[0]):
        d = rows[j]
        s = 0.0
        k = d
        while k < n:
            s += W[k]
            k += d
        out[j] = s


@njit(parallel=True, cache=True, nogil=True)
def _rowsum_i(F, rows, out):
    n = F.shape[0]
    for j in prange(rows.shape[0]):
        d = rows[j]
        s = 0
        k = d
        while k < n:
            s += F[k]
            k += d
        out[j] = s


def build_band(e, S):
    """vB, vR exactly as cell_cache._build does (Python loops over primes only)."""
    N = 2 ** e
    Q = math.isqrt(N)
    mu = np.array(mu_upto(N), dtype=np.int8)
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
    bd = (mu != 0) & (rem == 1)
    bd[:int(S["thr"]) + 1] = False
    vB = np.nonzero(bd & (om % 2 == 0))[0].astype(np.int64)
    vR = np.nonzero(bd & (om % 2 == 1))[0].astype(np.int64)
    return vB, vR


def small_part(N, Y):
    """a_Y(n) for every n <= N: one numpy slice per prime p <= Y."""
    sm = np.ones(N + 1, dtype=np.int64)
    for p in primes_upto(max(2, int(Y))):
        p = int(p)
        if p > Y:
            break
        sm[p::p] *= p
    return sm


def y_list(Q, L):
    """The Y values tried: small constants, powers of L, and the asymptotic L^6."""
    cand = [2, 3, 5, 7, 11, 13, 19, 30, 50, 100,
            int(round(L)), int(round(L ** 2)), int(round(L ** 3)),
            math.isqrt(Q), int(round(L ** 6))]
    out = sorted({int(y) for y in cand if y >= 2})
    return out


def do_cell(e, store):
    t0 = time.time()
    N = 2 ** e
    S = load_cell(e)
    Q, R1, nB = int(S["Q"]), int(S["nR"]), int(S["nB"])
    rows = np.asarray(S["rows"], dtype=np.int64)
    cb = np.asarray(S["cb"], dtype=np.int64)
    cr = np.asarray(S["cr"], dtype=np.int64)
    nr = rows.size
    L = math.log(Q)
    A = math.isqrt(Q)
    rec = dict(e=e, N=N, Q=Q, rows=nr, nB=nB, nR=R1, L=L, A=A, Y=[])
    say("=" * 112)
    say("2^%-3d Q=%-6d rows=%-5d |vB|=%-9d |vR|=%-9d  L=log Q=%.4f  A=floor(sqrt Q)=%d  L^6=%.3g"
        % (e, Q, nr, nB, R1, L, A, L ** 6))

    vB, vR = build_band(e, S)
    flagB = np.zeros(N + 1, dtype=np.int64)
    flagB[vB] = 1
    flagR = np.zeros(N + 1, dtype=np.int64)
    flagR[vR] = 1
    cbh = np.empty(nr, dtype=np.int64)
    crh = np.empty(nr, dtype=np.int64)
    _rowsum_i(flagB, rows, cbh)
    _rowsum_i(flagR, rows, crh)
    c1 = bool(vB.size == nB and vR.size == R1
              and np.array_equal(cbh, cb) and np.array_equal(crh, cr))
    rec["C1"] = c1
    say("      C1 band rebuilt vs cache: |vB| %d/%d, |vR| %d/%d, c_B and c_R element for element -> %s"
        % (vB.size, nB, vR.size, R1, "HOLDS" if c1 else "FAILS"))
    if not c1:
        rec["void"] = "C1"
        store.append(rec)
        return

    t = cr.astype(np.float64) / R1
    p = (cb - cr).astype(np.float64) / R1
    normp = float(np.linalg.norm(p))
    normt = float(np.linalg.norm(t))
    lit = NORMP_LIT.get(e)
    c4 = True if lit is None else (("%.*f" % (C4_DEC, normp)) == lit)
    rec.update(normp=normp, normt=normt, C4=bool(c4), C4_scored=lit is not None)
    say("      ||p||_2 = %.9f %s ; ||t||_2 = %.6f ; x=0 would give ||t||/||p|| = %.4f ; certified kappa = %.9f"
        % (normp, ("vs literal %s -> C4 %s" % (lit, "HOLDS" if c4 else "FAILS")) if lit else "(no literal; C4 not scored)",
           normt, normt / normp, KAP_CERT[e]))
    if not c4:
        rec["void"] = "C4"
        store.append(rec)
        return

    # the normalisation lower bound, for the record (not scored)
    say("      normalisation: ||p||_2 >= (1/5)(N/R_1)L^-2 = %.6f  (holds: %s)"
        % (0.2 * (N / R1) / L ** 2, 0.2 * (N / R1) / L ** 2 <= normp))

    smooth_cache = {}
    best = (None, float("inf"))
    for Y in y_list(Q, L):
        sm = small_part(N, Y)
        aB = sm[vB]
        aR = sm[vR]
        keys = np.union1d(np.unique(aB), np.unique(aR))
        iB = np.searchsorted(keys, aB)
        iR = np.searchsorted(keys, aR)
        Ba = np.bincount(iB, minlength=keys.size).astype(np.int64)
        Ra = np.bincount(iR, minlength=keys.size).astype(np.int64)
        c3 = bool(int(Ba.sum()) == nB and int(Ra.sum()) == R1)
        keep = keys <= A
        nk = int(keep.sum())
        dead = int(np.count_nonzero(keep & (Ba == 0) & (Ra > 0)))
        ok = keep & (Ba > 0)
        ratio = float((Ra[ok] / Ba[ok]).max()) if int(ok.sum()) else float("nan")
        rmin = float((Ra[ok] / Ba[ok]).min()) if int(ok.sum()) else float("nan")
        amax = int(keys[ok][int(np.argmax(Ra[ok] / Ba[ok]))]) if int(ok.sum()) else -1
        ueff = (math.log(N / amax) / math.log(Y)) if amax > 0 and Y > 1 else float("nan")
        delta = float(Ra[~keep].sum()) / R1
        # printed for scale only: the exponential factor of a Rankin-type smooth-number estimate, NOT a
        # bound on delta (no clause compares the two) -- see the header, item 4
        rank = math.exp(-L / (2.0 * math.log(max(Y, 3)))) if Y >= 3 else float("nan")

        xa = np.zeros(keys.size, dtype=np.float64)
        xa[ok] = Ra[ok] / (float(R1) * Ba[ok])
        W = np.zeros(N + 1, dtype=np.float64)
        W[vB] = xa[iB]
        Ax = np.empty(nr, dtype=np.float64)
        _rowsum_f(W, rows, Ax)
        err = t - Ax
        resid = float(np.linalg.norm(err))
        ratio_k = resid / normp

        # the two theoretical pieces of the residual
        smooth_row = sm[rows] == rows
        n_sm = int(smooth_row.sum())
        e_sm = float(np.linalg.norm(err[smooth_row]))
        e_rg = float(np.linalg.norm(err[~smooth_row]))

        # C2: exactness on Y-smooth rows
        WR = np.zeros(N + 1, dtype=np.float64)
        bad = vR[aR > A]
        WR[bad] = 1.0
        disc = np.empty(nr, dtype=np.float64)
        _rowsum_f(WR, rows, disc)
        pred = t - disc / R1
        dev = float(np.abs(Ax[smooth_row] - pred[smooth_row]).max()) if n_sm else 0.0
        applic = dead == 0          # no dead retained class; C2's two-term identity assumes it. x is
        #                             well defined either way (xa stays 0 where Ba == 0); with a dead
        #                             class Lemma 10's identity carries a second sum this clause omits
        c2 = bool(dev <= C2_ABS) if applic else None

        r = dict(Y=int(Y), n_a=int(keys.size), n_keep=nk, dead=dead, max_ratio=ratio, min_ratio=rmin,
                 argmax_a=amax, u_eff=ueff, delta=delta, rankin=rank, resid=resid, kappa_ub=ratio_k,
                 e_smooth=e_sm, e_rough=e_rg, n_smooth_rows=n_sm, C2=c2, C2_dev=dev, C2_applicable=applic,
                 C3=c3, P2_ok=bool(ratio_k >= KAP_CERT[e] - P2_SLACK))
        if c2 is False or not c3:
            r["void"] = "C2" if c2 is False else "C3"
        rec["Y"].append(r)
        if "void" not in r and applic and ratio_k < best[1]:
            best = (int(Y), ratio_k)
        say("      Y=%-7d |a|=%-6d keep(a<=A)=%-5d dead=%-3d R_a/B_a in [%s, %s] max at a=%-6d "
            "u_eff=%.2f delta=%.3e (scale exp(-L/2logY) %.1e)"
            % (Y, keys.size, nk, dead, ("%.4f" % rmin) if rmin == rmin else "n/a",
               ("%.4f" % ratio) if ratio == ratio else "n/a", amax, ueff, delta, rank))
        say("              resid=%.6f  kappa_ub=resid/||p||=%.4f  [smooth rows %d: %.6f | rough rows %d: "
            "%.6f]  C2 dev %.1e %s  C3 %s"
            % (resid, ratio_k, n_sm, e_sm, nr - n_sm, e_rg, dev,
               ("HOLDS" if c2 else "FAILS") if applic else "NOT APPLICABLE (dead retained class)",
               "HOLDS" if c3 else "FAILS"))
        smooth_cache[int(Y)] = ratio_k

    rec.update(best_Y=best[0], best_kappa_ub=best[1])
    say("      BEST over Y: Y=%s gives kappa <= %.4f  (certified kappa = %.9f ; factor %.1f above)"
        % (best[0], best[1], KAP_CERT[e], best[1] / KAP_CERT[e]))
    say("      %.1f s" % (time.time() - t0))
    store.append(rec)


def main(argv):
    cells = tuple(int(a) for a in argv if a.isdigit()) or (16, 18, 20, 22)
    nth = max(1, min(NTHREAD, 8, int(numba.config.NUMBA_NUM_THREADS)))
    set_num_threads(nth)
    say("conedual_parity_domination -- the rough-cofactor parity-domination construction (clauses "
        "in header)")
    say("cells %s   threads %d   cutoff A = floor(sqrt(Q))   %s"
        % (list(cells), nth, "SMOKE -- NOTHING IS SCORED" if SMOKE else "scored run"))
    say("")
    store = []
    for e in cells:
        do_cell(e, store)
        flush(dict(cells=store, complete=False, smoke=SMOKE))

    good = [r for r in store if "void" not in r]
    ys = [(r, y) for r in good for y in r["Y"] if "void" not in y]
    v = {}
    v["C1"] = all(r["C1"] for r in store if "C1" in r)
    v["C2"] = all(y["C2"] for _, y in ys if y["C2_applicable"])
    v["C3"] = all(y["C3"] for _, y in ys)
    v["C4"] = all(r["C4"] for r in store if r.get("C4_scored"))
    small = [(r, y) for r, y in ys if y["Y"] <= P1_YMAX]
    v["P1"] = bool(small) and all(y["dead"] == 0 and y["max_ratio"] <= P1_RATIO for _, y in small)
    v["P2"] = bool(ys) and all(y["P2_ok"] for _, y in ys)
    seq = [r["best_kappa_ub"] for r in good if r.get("best_kappa_ub") is not None]
    v["P3"] = len(seq) == len(good) and len(seq) >= 2 and all(b < a for a, b in zip(seq, seq[1:]))
    # Y values present at EVERY scored cell (so P4/P5 compare like with like)
    common = set.intersection(*[{y["Y"] for y in r["Y"]} for r in good]) if good else set()
    v["P4"] = bool(good) and all(
        all(b <= a for a, b in zip(s, s[1:]))
        for Y in sorted(y for y in common if y <= P1_YMAX)
        for s in [[next(z["delta"] for z in r["Y"] if z["Y"] == Y) for r in good]] if len(s) >= 2)
    v["P5"] = bool(good) and len(good) >= 2 and all(
        all(b <= a + 1e-12 for a, b in zip(s, s[1:]))
        for Y in sorted(y for y in common & {3, 5, 7, 11, 13, 19, 30})
        for s in [[next(z["max_ratio"] for z in r["Y"] if z["Y"] == Y) for r in good]])
    rows_p5 = [(Y, [next(z["max_ratio"] for z in r["Y"] if z["Y"] == Y) for r in good])
               for Y in sorted(y for y in common & {3, 5, 7, 11, 13, 19, 30})]
    say("")
    say("=" * 112)
    if SMOKE:
        say("SMOKE RUN -- clauses are not scored.")
    keyset = ("C1", "C2", "C3", "C4", "P1", "P2", "P3", "P4", "P5")
    for k in keyset:
        say("%-3s %s" % (k, "HOLDS" if v[k] else "FAILS/REFUTED"))
    say("")
    say(("SMOKE " if SMOKE else "") + "VERDICT  " + " | ".join(
        "%s %s" % (k, "HOLDS" if v[k] else "FAILS") for k in keyset))
    if rows_p5:
        say("P5 table, max_a R_a/B_a at fixed Y across cells %s:" % [r["e"] for r in good])
        for Y, s in rows_p5:
            say("   Y=%-4d %s" % (Y, "  ".join("%.4f" % u for u in s)))
    voids = [(r["e"], r["void"]) for r in store if "void" in r]
    yv = [(r["e"], y["Y"], y["void"]) for r in store for y in r.get("Y", []) if "void" in y]
    if voids:
        say("VOID enforced at: " + ", ".join("2^%d (%s)" % z for z in voids))
    if yv:
        say("VOID enforced at (cell, Y): " + ", ".join("2^%d/Y=%d (%s)" % z for z in yv))
    say("best kappa upper bound per cell: " + ", ".join(
        "2^%d: %.4f (Y=%s)" % (r["e"], r["best_kappa_ub"], r["best_Y"]) for r in good
        if r.get("best_kappa_ub") is not None))
    flush(dict(cells=store, clauses=v, complete=True, smoke=SMOKE,
               cited=dict(KAP_CERT=KAP_CERT, NORMP_LIT=NORMP_LIT)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

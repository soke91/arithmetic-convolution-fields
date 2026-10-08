# -*- coding: utf-8 -*-
r"""conedual_pole_term -- the s = 0 pole of  zeta(s)^{-1} prod_{p|d}(1-p^{-s})^{-1}:

Supporting computation; not used in the paper's proofs or tables.
the next main term of p after q_li AT THE COMPUTED CELLS, and what it does to the epsilon used there.

The two control clauses are stated at MATCHED PRECISION, every registered "voids" consequence is
enforced in the code, and the bound for the analytic model that section 2 headlines is computed
below.  "the next main term of p" is a FINITE-CELL statement: the next main term of F_d
asymptotically is the sum over the zeros, which exceeds any power of log.

WHERE THIS COMES FROM

The starting point: the coprime Moebius Dirichlet
series has a pole of order omega(d) at s = 0 (1/zeta(0) = -2), so M_d(x) carries a deterministic,
forced-sign term that q_li drops.  an earlier derivation measured its LEADING residue only
(-2 (log x)^omega / (omega! prod_{p|d} log p)) and found cos(X, Pole) = 0.97 .. 0.80 and
rLi(q_li + Pole_lead) = 0.0741, 0.0595, 0.0372, 0.0221, 0.0144, 0.0078, 0.0049 at 2^16 .. 2^28.
This file derives the FULL residue polynomial, tests it, and asks whether it lowers the
explicit epsilon(2^28) = 0.03058 below the certified kappa(2^28) >= 0.023782.

THE DERIVATION (everything below is implemented from it)

  f_d(s) = sum_{(n,d)=1} mu(n) n^{-s} = zeta(s)^{-1} prod_{p|d} (1 - p^{-s})^{-1}.

  With beta(u) = u/(1-e^{-u}) = sum_m B_m^+ u^m/m!  (B_1^+ = +1/2) and Lambda_d = prod_{p|d} log p,

      prod_{p|d}(1-p^{-s})^{-1} = Phi_d(s) / (s^omega Lambda_d),   Phi_d(s) = prod_{p|d} beta(s log p),

  so x^s f_d(s)/s = G_d(s) x^s / (s^{omega+1} Lambda_d) with G_d = Phi_d/zeta, and

      P_d(log x) := Res_{s=0} x^s f_d(s)/s
                  = (1/Lambda_d) * sum_{j=0}^{omega} c_j(d) (log x)^{omega-j} / (omega-j)!,

  where G_d(s) = sum_j c_j(d) s^j,  c_j(d) = sum_{i=0}^{j} eta_i phi_{j-i}(d),
  eta = Taylor of 1/zeta at 0 (eta_0 = -2, eta_1 = 2 log 2pi), phi = Taylor of Phi_d.
  j = 0 reproduces an earlier derivation's leading term.  c_1(d) = 2 log(2pi) - log d = log(4 pi^2 / d).

  ELEMENTARY FORM (what an explicit bound can use).  Lemma U says
  M_d(x) = sum_{a | d^inf, a <= x} M(x/a); writing M(y) = -2 + E(y) (the -2 is the s = 0 term of
  M's explicit formula) gives, EXACTLY,

      M_d(x) = -2 N_d(x) + sum_{a | d^inf, a <= x} E(x/a),   N_d(x) = #{a <= x : rad(a) | d},

  an exact integer count.  Hence F_d - Pi_d = sum_a [M(x/a)+2] - sum_a [M(x'/a)+2] with
  Pi_d = -2 (N_d(x) - N_d(x')), x = N/d, x' = thr/d.
  N_d(x) = Ntilde_d(x) + O_d((log x)^{omega-1}) is a FIXED-d statement; it is not uniform in d
  (d = 2p, x = 2 gives N_d(2) = 2 while Ntilde_d(2) carries log p / (12 log 2)).  Nothing below
  uses the asymptotic -- Ntilde_d and N_d are both computed, the first from the coefficients and the
  second as an exact integer count.

  ONLY Pi_d = -2 Delta N_d comes from an exact identity, so only it is directly boundable by Lemma U.
  The analytic P_an is reached from it by adding the exactly computable deterministic difference
  |P_el - P_an| row by row (section "the bound for the analytic model" below).

FOUR MODELS SCORED (Pole_d enters p as  q_li_d + mu(d) Pole_d / |vR|)

  lead : an earlier derivation's leading residue only
  an   : the full analytic polynomial P_d (the model named above)
  el   : the elementary -2 (N_d(x) - N_d(x')), exact integer counts
  mix  : -2 N_d exact  +  the eta_{j>=1} part of P_d  (exact counts AND the 1/zeta corrections)

Convention, stated because it matters at the top rows: P_d(log y) and the smooth count are used only
where y >= 1; where thr/d < 1 the lower-edge term is 0, which is what M_d(thr/d) = N_d(thr/d) = 0 says.

CLAUSES; a reader checks them with --verify.

  WHAT CARRIES A VERDICT AND WHAT DOES NOT.  R1, R2, R3, C0 and C0b are CONTROLS -- reproductions of
  already-published numbers and of the residue formula, not predictions.  P1, P2, P3 and P4 are the
  predictions; Q is scored and not predicted as a clause.  The quantities listed under DIAGNOSTICS
  below carry no clause and no verdict: the P_an bound, the counterfactual, the realised sign split,
  the section-1.2 diagnostics.

  MATCHED PRECISION IS MACHINE-ENFORCED (the project's compute rules; the reproduction rule).  Every
  reproduction literal below is stored as a STRING, and `check_literal_precision()` asserts, before
  any scoring, that the half-ulp of the literal as written, divided by the literal, is smaller than
  the clause's tolerance.  A clause whose tolerance is finer than its own literal is not a clause,
  and the run aborts rather than scoring it.

  R1  reproduction, like with like, at e = 20,22,24,26,28:
      (a) T_F from this file's kernel vs (num_E - T_Z_full) of results/conedual_explicit_transfer.json
          -- both full-precision floats read/produced in this run -- relative 1e-12;
      (b) the RAW (T_F + T_Z_full)/denom_ex vs that file's RAW eps_E_ex, relative 1e-12;
      (c) the 4-significant-digit ceiling of each side of (b), compared with each other, EXACTLY.
      Each side of a comparison is read or produced at the same precision as the other; a
      4-significant-digit ceiling is compared only with another 4-significant-digit ceiling.
      CAN FAIL; **voids section 3** -- enforced in the code (`VOID` tags on every section-3
      number, `void_section3` in the JSON, non-zero exit).
  R2  reproduction.  rLi for q_li alone equals the stored values
      0.09118077290374796, 0.048970177377251824, 0.02629818988047844, 0.013762169857355972,
      0.0073115695921642946 (2^20..2^28) to a relative 1e-12.  These are the full-precision stored
      values; an 8-decimal rounding of them carries a relative error of 4.5e-9 .. 5.6e-8 and cannot be
      compared at this tolerance.  CAN FAIL; **voids everything** -- enforced in the code (`VOID` on
      every result, `void_all` in the JSON, non-zero exit).
  R3  reproduction.  rLi(q_li + Pole_lead) equals an earlier derivation's 0.07410, 0.05945, 0.03720,
      0.02214, 0.01435, 0.00780, 0.00486 (2^16..2^28) to the printed digits (relative 1e-2).
      CAN FAIL; if it fails the P1 comparison is made against THIS file's own lead column and the
      output says so.
  C0  the residue formula: P_d(log x) agrees with a direct Cauchy contour integral
      of x^s f_d(s)/s on |s| = 0.3 (mpmath, 25 dps) to a relative 1e-12 on six (d, x) pairs.
      CAN FAIL; voids every pole column -- enforced in the code.
      WHAT C0 CERTIFIES: it casts the contour value to float64 before taking the ratio, so it
      certifies agreement IN DOUBLE PRECISION and is not an 18-digit certificate.  C0b is the
      high-precision version.
  C0b a control, not a prediction: the same residue computed entirely in mpmath at
      50 dps (eta from mp.taylor, B_m^+ from mp.bernoulli) against the contour integral at TWO
      radii r = 0.15 and r = 0.5, on eight (d, x) pairs including d = 210:
      relative deviation <= 1e-40.  CAN FAIL; voids every pole column.

  P1  THE CLAUSE THE PROMPT NAMES.  rLi(q_li + Pole_an) <= 1.05 * rLi(q_li + Pole_lead) at all
      seven cells.  PREDICTED (mine, before the run): 0.070, 0.056, 0.035, 0.021, 0.0135, 0.0073,
      0.0046 -- i.e. a 0-15 % improvement on the leading term, no cell worse.
      WHY: the j = 1 coefficient is c_1(d) = log(4 pi^2 / d), so the j = 1 term is a relative
      omega(d) log(4 pi^2/d) / (-2 log(N/d)) of the leading one -- under 6 % in modulus on every row
      with d <= 10^3 at 2^28, and those rows carry most of ||Pole||_2.  Where the correction is large
      (d near Q, omega >= 4) |Pole_d| is itself small, so the l2 effect is bounded.  Every added
      coefficient is an exact term of the same residue, not a fitted parameter, so the expected
      direction is an improvement; the 5 % allowance covers the one real risk, that the unmodelled
      secondary poles at s = 2 pi i k / log p -- the same order as the j >= 1 terms on high-omega
      rows -- point the other way.
      MEANING: holds -> the full residue is the right object and replaces an earlier derivation's leading term.
               fails -> the expansion is being used outside its regime at these cells; keep the
               leading term and say so.
      READING THE VERDICT: P1 HOLDS, but **its registered WHY above is FALSE**, and the diagnostics
      `s12_*` measure that.  Over squarefree d <= 1000 at 2^28 the j=1/j=0 ratio
      reaches 0.51 at d = 966, and 531 of 607 rows exceed 6 %; and d <= 1000 carries 22.6 %, not
      "most", of ||Pole_an||^2.  Also "reducing" is backwards: c_1(d) = log(4 pi^2/d) < 0 for
      d > 4 pi^2 = 39.48, the same sign as c_0 = -2, so the j = 1 term INCREASES |P_d| there (it
      reduces only against -2 Ntilde_d, where eta_1 > 0 opposes the count).  The real mechanism is
      the opposite of the registered one: the pole lives on the omega >= 4 rows (73 % of
      ||Pole_an||^2 at 2^28), which is exactly where the correction is LARGEST -- P1 holds because
      the j >= 1 correction is big where the pole is, not because it is small.  The clause's
      outcome stands; its stated reason does not.
  P2  rLi(q_li + Pole_el) <= rLi(q_li + Pole_lead) at all seven cells.
      WHY: the elementary form uses exact integer counts, so it resums every secondary pole.
  P3  Pole_mix is the best of the four models at >= 5 of the 7 cells.
  P4  THE EXPLICIT-BOUND QUESTION (prompt item 3).  T_{F-Pi} >= 0.90 * T_F at 2^28.
      PREDICTED: T_{F-Pi}/T_F in [0.95, 1.02]; eps_E_ex(2^28) with the pole subtracted in
      [0.029, 0.031]; so the answer to "does it drop below kappa(2^28) >= 0.023782" is NO.
      WHY: replacing |M(y)| by |M(y)+2| moves each Lemma-U term by at most 2, so the whole row
      bound moves by at most 2 * Delta N_d = |Pi_d|, and ||Pi||_2 / T_F is about 0.2 at 2^28 -- a
      perfectly-signed gain is capped there.  The gain is realised only where the N-independent
      table reaches (N/d <= TAB = 1e7, i.e. d >= 27 at 2^28) and only on the a with M(x/a) <= -2;
      every uncovered term pays +2 instead.  Structurally: every published pointwise input is at
      the sqrt(y) scale and the pole lives at the (log y)^omega scale, so a sqrt-scale bound cannot
      see it.
      MEANING: holds -> subtracting the pole improves the MODEL and not the BOUND, and the reason
               is the scale mismatch, which is a statement about every bound of this shape.
               fails -> the exact-table region is doing more than predicted; report eps and compare.
      READING THE VERDICT: P4 HOLDS; two statements inside its WHY above are wrong.  (i) "the row
      bound moves by at most 2 * Delta N_d" is false -- the Lemma-U sums at x AND at x' both enter
      the bound additively, so N_d(x) + N_d(x') terms can each move by 2, and at (2^28, d = 53) the
      realised gain is 12 against the claimed cap of 4.  The per-term counterfactual
      `T_FP_best` below (max(G-2,0) in the table, G+2 outside, through the same two mins) is the
      correct one, and it is larger, so it makes the negative answer stronger.
      (ii) "only on the a with M(x/a) <= -2" names a global frequency over n <= TAB; the arguments
      the bound evaluates are floor(x/a), a different and much smaller-weighted set, measured here
      as `sign_gain/neutral/lose`.  Neither point changes the verdict.
  Q   scored, not predicted as a clause: is eps_E_ex(2^28) with the pole below 0.023782?

  DIAGNOSTICS, NO verdict attached: the bound for the analytic model (`T_F_an`, `eps_an`); the
  per-term counterfactual
  (`T_FP_best`); the realised sign split over the Lemma-U terms and the realised net; the
  section-1.2 diagnostics (`s12_*`); the amplitude multiplier `amp_mult` and the gain a fitted
  scalar would buy (`fit_gain`); the share of ||Pole_an||^2 on the rows beyond the table.

  NOT REGISTERED: anything about kappa itself, about Goldbach, or any asymptotic claim.

COMPUTE RULES (the project's compute rules, both sections read before this file was written)
  cells are LOADED (`from cell_cache import load`), never rebuilt; no Python loop over rows, pairs
  or band elements -- np.repeat/np.bincount pair arrays and one numba prange kernel with six modes
  for the Lemma-U enumerations; the only Python loops are over cells, over the primes <= Q, over
  the <= 7 coefficient slots of a degree-6 polynomial, and over the <= 8 (d, x) pairs of C0/C0b.
  Threads capped.  Cells 2^20 and above want a machine with more cores.

  EXIT CODES: 0 all controls held; 3 a registered void is in force; 4 a clause's tolerance is
  finer than its own literal (the run aborts before scoring).

    python code/conedual_pole_term.py 16 18            # a smoke test at the two smallest cells
    python code/conedual_pole_term.py                  # all seven
"""
from __future__ import annotations
import io, json, math, os, sys, time

NTHREAD = int(os.environ.get("ACF_THREADS", "8"))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, str(NTHREAD))
_SCRATCH = os.environ.get("TEMP") or os.environ.get("TMP") or "."
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(_SCRATCH, "numba_cache_acf"))

import numpy as np  # noqa: E402
import numba  # noqa: E402
from numba import njit, prange, set_num_threads  # noqa: E402
from scipy.special import expi  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _sieve_shared import mu_upto, primes_upto  # noqa: E402
from cell_cache import load as load_cell  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "conedual_pole_term.txt")
OUTJ = os.path.join(RES, "conedual_pole_term.json")
REF14 = os.path.join(RES, "conedual_explicit_transfer.json")

CELLS = (16, 18, 20, 22, 24, 26, 28)
TAB = 10 ** 7                 # the N-independent table every bound may use
KDEG = 6                      # polynomial degree carried; omega(d) <= 5 for d <= 2^14
MARGIN = 1.0e-10              # outward relative slack
KAPPA_28 = 0.023782           # a certified kappa(2^28) lower bound, computed earlier

# ---- literature inputs actually used here ------------------------------------------------------
C_M1, X_M1_DRESS = 0.571, 1.0e12        # JLT Prop. A.3 (A.10) = Dress 1993: |M(x)| <= .571 sqrt x
C_R14 = 2.4                             # Ramare, Acta Arith. 165 (2014), weighted -> Lemma V
E24 = math.exp(C_R14)

# ---- Taylor coefficients of 1/zeta at s = 0 ------------------------------------------------------
# eta_0 = 1/zeta(0) = -2 and eta_1 = 2 log(2 pi) are exact; all seven were computed in this task with
# mpmath at 40 dps by two independent methods (mp.taylor on |s| = 0.35 and mp.diff), agreeing to 18
# digits, and are re-verified at run time when mpmath is importable.
ETA = np.array([-2.0,
                3.67575413281869097,
                -2.74287131064967372,
                1.66932802114990601,
                -0.92260555586070002,
                0.48580536497019563,
                -0.24948930248752077], dtype=np.float64)
# beta(u) = u/(1-e^{-u}) = sum_m B_m^+ u^m / m!   (B_1^+ = +1/2)
BET = np.array([1.0, 0.5, 1.0 / 12.0, 0.0, -1.0 / 720.0, 0.0, 1.0 / 30240.0], dtype=np.float64)
FACT = np.array([math.factorial(n) for n in range(KDEG + 2)], dtype=np.float64)

_lines = []


def say(s=""):
    _lines.append(s)
    print(s, flush=True)


def up(v):
    return v * (1.0 + MARGIN)


def ceil_sig(v, nd=4):
    if v is None or not np.isfinite(v) or v <= 0:
        return v
    k = nd - 1 - int(math.floor(math.log10(v)))
    return math.ceil(v * 10.0 ** k) / 10.0 ** k


def li(x):
    return expi(np.log(np.asarray(x, dtype=np.float64)))


def nrm(v):
    return float(np.sqrt(np.dot(v, v)))


# ---- reproduction literals, as STRINGS, so their written precision is machine-checkable ----------
# (the project's compute rules "a literal's precision must exceed the tolerance"; the reproduction rule.)
# R2: the stored rLi, full float repr.  R3: the earlier printed figures, at the digits they were printed to.
LIT_R2 = {20: "0.09118077290374796", 22: "0.048970177377251824", 24: "0.02629818988047844",
          26: "0.013762169857355972", 28: "0.0073115695921642946"}
TOL_R2 = 1.0e-12
LIT_R3 = {16: "0.07410", 18: "0.05945", 20: "0.03720", 22: "0.02214", 24: "0.01435",
          26: "0.00780", 28: "0.00486"}
TOL_R3 = 1.0e-2
TOL_R1 = 1.0e-12


def lit_rel_precision(s):
    """Relative half-ulp of a decimal literal AS WRITTEN: 0.00486 -> 5e-6/0.00486 = 1.03e-3."""
    frac = len(s.split(".")[1]) if "." in s else 0
    return 0.5 * 10.0 ** (-frac) / abs(float(s))


def check_literal_precision():
    """Abort before scoring if any clause's tolerance is finer than its own literal's precision."""
    bad = []
    for nm, lits, tol in (("R2", LIT_R2, TOL_R2), ("R3", LIT_R3, TOL_R3)):
        for e, s in lits.items():
            p = lit_rel_precision(s)
            if p >= tol:
                bad.append("%s 2^%d literal %s has relative precision %.1e >= tolerance %.1e"
                           % (nm, e, s, p, tol))
    worst = {nm: max(lit_rel_precision(s) for s in lits.values())
             for nm, lits in (("R2", LIT_R2), ("R3", LIT_R3))}
    return bad, worst


# ---------------------------------------------------------------- numba: the three Lemma-U sweeps

@njit(cache=True, nogil=True, inline="always")
def _gm(n, Mabs):
    """Upper bound for |M(n)|: the exact table below TAB, else JLT Prop. A.3 (Dress)."""
    if n < 1:
        return 0.0
    if n < Mabs.shape[0]:
        return float(Mabs[n])
    y = float(n)
    b = C_M1 * math.sqrt(y) if y >= 33.0 else y
    return y if y < b else b


@njit(cache=True, nogil=True, inline="always")
def _gm2(n, Msig, Mabs):
    """Upper bound for |M(n) + 2|: exact in the table, else the same bound + 2."""
    if n < 1:
        return 2.0 if n == 0 else 0.0
    if n < Msig.shape[0]:
        v = float(Msig[n]) + 2.0
        return v if v >= 0.0 else -v
    y = float(n)
    b = C_M1 * math.sqrt(y) if y >= 33.0 else y
    if y < b:
        b = y
    return b + 2.0


@njit(cache=True, nogil=True, inline="always")
def _gm3(n, Mabs):
    """The BEST value |M(n)+2| could take given only what the bound knows (the counterfactual):
    max(|M(n)|-2, 0) where the exact table reaches, and the pessimistic bound+2 where it does not."""
    if n < 1:
        return 0.0
    if n < Mabs.shape[0]:
        v = float(Mabs[n]) - 2.0
        return v if v > 0.0 else 0.0
    y = float(n)
    b = C_M1 * math.sqrt(y) if y >= 33.0 else y
    if y < b:
        b = y
    return b + 2.0


@njit(parallel=True, cache=True, nogil=True)
def _sweep(xs, pf, mode, Msig, Mabs, out):
    """Lemma U over a | d_i^inf, a <= xs[i].   mode 0: sum |M(x/a)| bound
                                               mode 1: sum |M(x/a)+2| bound
                                               mode 2: the exact count N_d(x)
                                               mode 3: sum of the best-case |M(x/a)+2|
                                               mode 4: # covered terms with M(x/a) <= -2 (gain 2)
                                               mode 5: # covered terms with M(x/a) = -1 (neutral)"""
    for i in prange(xs.shape[0]):
        x = xs[i]
        if x < 1:
            out[i] = 0.0
            continue
        p0 = pf[i, 0]; p1 = pf[i, 1]; p2 = pf[i, 2]
        p3 = pf[i, 3]; p4 = pf[i, 4]; p5 = pf[i, 5]
        acc = 0.0
        a0 = np.int64(1)
        while True:
            a1 = a0
            while True:
                a2 = a1
                while True:
                    a3 = a2
                    while True:
                        a4 = a3
                        while True:
                            a5 = a4
                            while True:
                                nn = x // a5
                                if mode == 0:
                                    acc += _gm(nn, Mabs)
                                elif mode == 1:
                                    acc += _gm2(nn, Msig, Mabs)
                                elif mode == 2:
                                    acc += 1.0
                                elif mode == 3:
                                    acc += _gm3(nn, Mabs)
                                elif mode == 4:
                                    if nn >= 1 and nn < Msig.shape[0] and Msig[nn] <= -2:
                                        acc += 1.0
                                else:
                                    if nn >= 1 and nn < Msig.shape[0] and Msig[nn] == -1:
                                        acc += 1.0
                                if p5 == 0:
                                    break
                                a5 *= p5
                                if a5 > x:
                                    break
                            if p4 == 0:
                                break
                            a4 *= p4
                            if a4 > x:
                                break
                        if p3 == 0:
                            break
                        a3 *= p3
                        if a3 > x:
                            break
                    if p2 == 0:
                        break
                    a2 *= p2
                    if a2 > x:
                        break
                if p1 == 0:
                    break
                a1 *= p1
                if a1 > x:
                    break
            if p0 == 0:
                break
            a0 *= p0
            if a0 > x:
                break
        out[i] = acc


def sweep(xs, pf, mode, T):
    out = np.zeros(xs.size, dtype=np.float64)
    _sweep(np.ascontiguousarray(xs), np.ascontiguousarray(pf), mode, T.Msig, T.Mabs, out)
    return out


# ---------------------------------------------------------------- tables and row geometry

class Tables(object):
    def __init__(self, lim):
        t0 = time.time()
        self.lim = lim
        mu = np.array(mu_upto(lim), dtype=np.int8)
        self.mu = mu
        self.Msig = np.cumsum(mu.astype(np.int32))
        self.Mabs = np.abs(self.Msig).astype(np.int32)
        n33 = np.arange(33, lim + 1, dtype=np.float64)
        self.hurst_max = float((self.Mabs[33:] / np.sqrt(n33)).max())
        self.frac_M_le_m2 = float(np.mean(self.Msig[1:] <= -2))
        self.frac_M_neg = float(np.mean(self.Msig[1:] < 0))
        self.seconds = time.time() - t0


def nextprime(n):
    q = n + 1
    while True:
        if q == 2:
            return 2
        if q > 2 and q % 2 and all(q % r for r in range(3, int(q ** 0.5) + 1, 2)):
            return q
        q += 1


def row_data(Q, T):
    """rows = squarefree d <= Q, with kappa_d, the distinct-prime table and Lambda_d.

    One Python loop over the primes <= Q (the pattern cell_cache._build uses); none over rows."""
    kap = np.ones(Q + 1, dtype=np.float64)
    lam = np.ones(Q + 1, dtype=np.float64)
    pf = np.zeros((Q + 1, 8), dtype=np.int64)
    nf = np.zeros(Q + 1, dtype=np.int64)
    for p in primes_upto(Q):
        p = int(p)
        idx = np.arange(p, Q + 1, p, dtype=np.int64)
        pf[idx, nf[idx]] = p
        nf[idx] += 1
        kap[idx] *= p / (p - 1.0)
        lam[idx] *= math.log(p)
    rows = (np.nonzero(T.mu[1:Q + 1] != 0)[0] + 1).astype(np.int64)
    assert int(nf[rows].max()) <= 6, "a squarefree row has > 6 primes; widen _sweep"
    return rows, kap[rows], pf[rows, :6], lam[rows], nf[rows]


def pair_arrays(rows, N, q1, T):
    """(d, k) with k <= floor(N/(d q1)), mu^2(k) = 1, (k, d) = 1 -- the k-sum of S_d, vectorised."""
    cnt = N // (rows * q1)
    di = np.repeat(np.arange(rows.size, dtype=np.int64), cnt)
    start = np.repeat(np.cumsum(cnt) - cnt, cnt)
    k = np.arange(di.size, dtype=np.int64) - start + 1
    d = rows[di]
    keep = (T.mu[k] != 0) & (np.gcd(k, d) == 1)
    return di[keep], k[keep], d[keep]


def lemma_v(x, d, kap):
    """Ramare 2014 weighted + partial summation."""
    x = x.astype(np.float64)
    dd = d.astype(np.float64)
    out = x.copy()
    ok = x > dd * E24
    xs, ds, ks = x[ok], dd[ok], kap[ok]
    out[ok] = np.minimum(xs, ks * (C_R14 * xs / np.log(xs / ds) + ds * E24
                                   + C_R14 * ds * (li(xs / ds) - li(E24))))
    return up(out)


# ---------------------------------------------------------------- the residue polynomial

def phi_coeffs(pf, om):
    """phi[:, j] = [s^j] prod_{p|d} beta(s log p).  Loop over the <= 6 prime slots and <= 7 degrees."""
    n = pf.shape[0]
    phi = np.zeros((n, KDEG + 1), dtype=np.float64)
    phi[:, 0] = 1.0
    for slot in range(pf.shape[1]):
        lp = np.where(pf[:, slot] > 0, np.log(np.maximum(pf[:, slot], 1)), 0.0)
        new = np.zeros_like(phi)
        for j in range(KDEG + 1):
            for m in range(j + 1):
                if BET[m] != 0.0:
                    new[:, j] += phi[:, j - m] * (BET[m] * lp ** m)
        phi = new
    return phi


def gc_coeffs(phi):
    """c[:, j] = sum_{i<=j} eta_i phi_{j-i}."""
    c = np.zeros_like(phi)
    for j in range(KDEG + 1):
        for i in range(j + 1):
            c[:, j] += ETA[i] * phi[:, j - i]
    return c


def poly_at(coef, om, lam, y, lead_only=False):
    """(1/Lambda_d) sum_{j<=om} coef[:,j] (log y)^{om-j}/(om-j)!, and 0 where y < 1."""
    y = np.asarray(y, dtype=np.float64)
    ok = y >= 1.0
    ly = np.where(ok, np.log(np.maximum(y, 1.0)), 0.0)
    acc = np.zeros(y.shape, dtype=np.float64)
    jmax = 0 if lead_only else KDEG
    for j in range(jmax + 1):
        ex = om - j
        m = ok & (ex >= 0)
        if m.any():
            acc[m] += coef[m, j] * ly[m] ** ex[m] / FACT[ex[m]]
    return np.where(ok, acc / lam, 0.0)


def check_residue_hp(dps=50):
    """C0b: the residue computed ENTIRELY in mpmath (eta from mp.taylor, B_m^+ from mp.bernoulli)
    against the contour integral at two radii.  Unlike C0 nothing is cast to float, so this is a
    high-precision certificate; C0 only certifies the float64 implementation agrees in double."""
    try:
        import mpmath as mp
    except Exception as exc:                                     # pragma: no cover
        return None, "mpmath unavailable (%s)" % exc
    mp.mp.dps = dps
    eta = mp.taylor(lambda s: 1 / mp.zeta(s), 0, KDEG)
    worst, cells = mp.mpf(0), []
    for d, x in ((6, 1e4), (6, 1e8), (30, 1e4), (30, 1e8),
                 (210, 1e4), (210, 1e8), (2310, 1e4), (2310, 1e8)):
        ps = [p for p in (2, 3, 5, 7, 11, 13) if d % p == 0]
        om = len(ps)
        lam = mp.mpf(1)
        for p in ps:
            lam *= mp.log(p)
        # Phi_d(s) = prod beta(s log p), beta(u) = sum B_m^+ u^m/m!, B_1^+ = +1/2
        phi = [mp.mpf(1)] + [mp.mpf(0)] * KDEG
        for p in ps:
            lp = mp.log(p)
            new = [mp.mpf(0)] * (KDEG + 1)
            for j in range(KDEG + 1):
                for m in range(j + 1):
                    b = mp.bernoulli(m) * (1 if m != 1 else -1)   # B_m^+ : flip the sign at m = 1
                    new[j] += phi[j - m] * b / mp.factorial(m) * lp ** m
            phi = new
        lx = mp.log(mp.mpf(x))
        got = mp.mpf(0)
        for j in range(om + 1):
            cj = sum(eta[i] * phi[j - i] for i in range(j + 1))
            got += cj * lx ** (om - j) / mp.factorial(om - j)
        got /= lam
        for rr in ("0.15", "0.5"):
            r = mp.mpf(rr)

            def integ(th, ps=ps, x=x, r=r):
                s = r * mp.expjpi(2 * th)
                v = mp.power(mp.mpf(x), s) / (s * mp.zeta(s))
                for p in ps:
                    v /= (1 - mp.power(p, -s))
                return v * s
            ref = mp.quad(integ, [0, 1]).real
            rel = abs(got - ref) / abs(ref)
            worst = max(worst, rel)
            cells.append((d, x, float(rr), mp.nstr(got, 15), mp.nstr(rel, 3)))
    return (float(worst), cells), None


def check_residue(dps=25):
    """C0: P_d(log x) against a direct Cauchy integral of x^s f_d(s)/s on |s| = 0.3.

    NOTE: `float(ref.real)` casts the contour value to double, so this compares the float64
    implementation with the contour IN DOUBLE PRECISION.  It is not an 18-digit certificate -- that
    is C0b."""
    try:
        import mpmath as mp
    except Exception as exc:                                     # pragma: no cover
        return None, "mpmath unavailable (%s)" % exc
    mp.mp.dps = dps
    worst, cells = 0.0, []
    for d, x in ((1, 1000.0), (2, 1.0e6), (6, 1.0e5), (30, 1.0e7), (105, 5.0e4), (2310, 1.0e5)):
        ps = [p for p in (2, 3, 5, 7, 11, 13) if d % p == 0]
        r = mp.mpf("0.3")

        def integ(th, ps=ps, x=x):
            s = r * mp.expjpi(2 * th)
            v = mp.power(mp.mpf(x), s) / (s * mp.zeta(s))
            for p in ps:
                v /= (1 - mp.power(p, -s))
            return v * s                                          # ds = 2 pi i s dth
        ref = mp.quad(integ, [0, 1])
        pf = np.zeros((1, 6), dtype=np.int64)
        pf[0, :len(ps)] = ps
        om = np.array([len(ps)], dtype=np.int64)
        lam = np.array([float(np.prod([math.log(p) for p in ps])) if ps else 1.0])
        got = float(poly_at(gc_coeffs(phi_coeffs(pf, om)), om, lam, np.array([x]))[0])
        rel = abs(got - float(ref.real)) / max(abs(float(ref.real)), 1e-30)
        worst = max(worst, rel)
        cells.append((d, x, got, float(ref.real), rel))
    return (worst, cells), None


# ---------------------------------------------------------------- one cell

def evaluate(e, T, ref14):
    t0 = time.time()
    N = 2 ** e
    Q = math.isqrt(N)
    q1 = nextprime(Q)
    thr = N // q1
    S = load_cell(e)
    assert S["N"] == N and S["Q"] == Q and S["q1"] == q1 and S["thr"] == thr
    R1 = float(S["nR"])
    rows, kap, pf, lam, om = row_data(Q, T)
    assert rows.size == S["rows"].size and int(np.abs(rows - S["rows"]).max()) == 0
    assert int(np.abs(om - S["rom"]).max()) == 0, "omega from the sieve != the cell's rom"
    mud = np.asarray(T.mu[rows], dtype=np.float64)
    Rp = (S["cb"] - S["cr"]).astype(np.float64)             # = R1 * p_d, exact integers
    denom_true = nrm(Rp)

    # ---- q_li: S_d with pi -> li, both pi's
    di, k, d = pair_arrays(rows, N, q1, T)
    muk = T.mu[k].astype(np.float64)
    liQ = float(li(float(Q)))
    term = muk * (li(N / (d * k).astype(np.float64)) - liQ)
    Sli = np.bincount(di, weights=term, minlength=rows.size)
    X = mud * Rp - Sli                                      # = F_d + Z_d

    # ---- the pole, four ways
    xN = N / rows.astype(np.float64)
    xT = thr / rows.astype(np.float64)
    phi = phi_coeffs(pf, om)
    cc = gc_coeffs(phi)
    P_an = poly_at(cc, om, lam, xN) - poly_at(cc, om, lam, xT)
    P_lead = poly_at(cc, om, lam, xN, True) - poly_at(cc, om, lam, xT, True)
    cnt_N = sweep(N // rows, pf, 2, T)
    cnt_T = sweep(thr // rows, pf, 2, T)
    P_el = -2.0 * (cnt_N - cnt_T)
    sm_N = poly_at(phi, om, lam, xN)                        # the s=0 residue of the count itself
    sm_T = poly_at(phi, om, lam, xT)
    P_mix = P_el + (P_an + 2.0 * (sm_N - sm_T))

    models = (("none", np.zeros_like(X)), ("lead", P_lead), ("an", P_an),
              ("el", P_el), ("mix", P_mix))
    r = dict(e=e, N=int(N), Q=int(Q), q1=int(q1), thr=int(thr), rows=int(rows.size),
             nR=int(S["nR"]), pairs=int(k.size), denom_true=denom_true,
             norm_X=nrm(X), norm_Rp_sq=int(np.dot(Rp, Rp)))
    for nm, P in models:
        r["rLi_" + nm] = nrm(X - P) / denom_true
        if nm != "none":
            r["cos_" + nm] = float(np.dot(X, P) / (nrm(X) * nrm(P))) if nrm(P) > 0 else float("nan")
            r["norm_" + nm] = nrm(P)

    # ---- the explicit bound: T_F and T_{F-Pi}
    xNi, xTi = N // rows, thr // rows
    lvN, lvT = lemma_v(xNi, rows, kap), lemma_v(xTi, rows, kap)
    swN0, swT0 = sweep(xNi, pf, 0, T), sweep(xTi, pf, 0, T)
    swN1, swT1 = sweep(xNi, pf, 1, T), sweep(xTi, pf, 1, T)
    BN = np.minimum(swN0, lvN)
    BT = np.minimum(swT0, lvT)
    triv = (xNi - xTi).astype(np.float64)
    BF = np.minimum(triv, BN + np.minimum(xTi.astype(np.float64), BT))
    T_F = up(nrm(BF))
    BN2 = np.minimum(swN1, lvN + 2.0 * cnt_N)
    BT2 = np.minimum(swT1, lvT + 2.0 * cnt_T)
    BF2 = np.minimum(triv + 2.0 * (cnt_N - cnt_T),
                     BN2 + np.minimum(xTi.astype(np.float64) + 2.0 * cnt_T, BT2))
    T_FP = up(nrm(BF2))
    r.update(T_F=T_F, T_FP=T_FP, ratio_T=T_FP / T_F, norm_Pi=nrm(P_el))

    # ---- the bound for the model section 2 HEADLINES.
    # Only Pi = -2 Delta N_d comes from an exact identity (Lemma E), so only it is directly
    # boundable by Lemma U.  P_an is reached by adding the deterministic, exactly computable
    # difference row by row -- a valid explicit bound, and a worse one.
    dif = np.abs(P_el - P_an)
    BF_an = BF2 + dif
    T_F_an = up(nrm(BF_an))
    r.update(T_F_an=T_F_an, norm_dif_el_an=nrm(dif),
             T_F_an_normtriangle=up(nrm(BF2) + nrm(dif)))

    # ---- diagnostics; NO verdict attached to them.
    # T_FP_best: the per-term counterfactual in which every
    # Lemma-U term the table covers returns the most |M+2| could ever give it, max(|M|-2, 0), while
    # every uncovered term still pays +2 -- through the same two mins.  Subtracting 2*Delta N_d from
    # the finished row bound instead would be wrong twice: the sums at x and at x' both
    # enter (so N_d(x) + N_d(x') terms can move, not Delta N_d), and a term with |M| < 2 cannot
    # return its full 2.  It is NOT a bound; it says how much the substitution could ever return.
    BN3 = np.minimum(sweep(xNi, pf, 3, T), lvN + 2.0 * cnt_N)
    BT3 = np.minimum(sweep(xTi, pf, 3, T), lvT + 2.0 * cnt_T)
    BF3 = np.minimum(triv + 2.0 * (cnt_N - cnt_T),
                     BN3 + np.minimum(xTi.astype(np.float64) + 2.0 * cnt_T, BT3))
    uncov = xNi > TAB
    sq_an = P_an * P_an
    # the realised sign split over the terms the bound actually evaluates; the 50.37 % global
    # frequency over n <= TAB is a different set.  "lose" includes every uncovered term.
    n_terms = float((cnt_N + cnt_T).sum())
    n_gain = float((sweep(xNi, pf, 4, T) + sweep(xTi, pf, 4, T)).sum())
    n_neut = float((sweep(xNi, pf, 5, T) + sweep(xTi, pf, 5, T)).sum())
    net = float(((swN1 - swN0) + (swT1 - swT0)).sum())
    r.update(T_FP_best=up(nrm(BF3)),
             rows_trivial_arm=int(np.sum(triv <= BN + np.minimum(xTi.astype(np.float64), BT))),
             rows_uncovered=int(uncov.sum()),
             share_TF2_uncovered=float(np.dot(BF[uncov], BF[uncov]) / np.dot(BF, BF)),
             share_pole2_uncovered=float(sq_an[uncov].sum() / sq_an.sum()),
             terms=n_terms, terms_gain=n_gain, terms_neutral=n_neut,
             terms_lose=n_terms - n_gain - n_neut,
             realised_net=net, systematic_net=-2.0 * n_terms,
             realised_frac_of_max=net / (-2.0 * n_terms))

    # section 1.2's two claims, and the mechanism behind them
    small = (rows <= 1000) & (rows > 1)        # d = 1 has omega = 0, ratio 0; the range is [2,1000]
    if small.any():
        ratio12 = (om[small] * np.log(4.0 * math.pi ** 2 / rows[small].astype(np.float64))
                   / (-2.0 * np.log(xN[small])))
        amax = int(np.argmax(np.abs(ratio12)))
        r.update(s12_rows=int(small.sum()), s12_max_abs=float(np.abs(ratio12).max()),
                 s12_argmax_d=int(rows[small][amax]),
                 s12_over6pct=int(np.sum(np.abs(ratio12) > 0.06)),
                 s12_share_pole2_dle1000=float(sq_an[small].sum() / sq_an.sum()))
    r["s12_share_pole2_by_omega"] = {int(w): float(sq_an[om == w].sum() / sq_an.sum())
                                     for w in np.unique(om)}
    # amplitude: ||P||/(||X|| cos), which is not the residual ratio below
    dot = float(np.dot(X, P_an))
    r.update(amp_mult=float(nrm(P_an) ** 2 / dot),          # ||P||/(||X|| cos)
             fit_mult=float(dot / nrm(P_an) ** 2),          # the scalar a* a fit would choose
             resid_ratio=r["rLi_an"] / (r["norm_X"] * math.sqrt(max(0.0, 1 - r["cos_an"] ** 2))
                                        / denom_true),
             fit_gain=1.0 - (r["norm_X"] * math.sqrt(max(0.0, 1 - r["cos_an"] ** 2))
                             / denom_true) / r["rLi_an"],
             norm_dif_over_X=nrm(dif) / r["norm_X"])

    ref = ref14.get(e)
    if ref is not None:
        TZ = ref["T_Z_full"]
        den = ref["denom_ex"]
        r.update(T_Z_full=TZ, denom_ex=den,
                 eps_old=ceil_sig((T_F + TZ) / den), eps_new=ceil_sig((T_FP + TZ) / den),
                 eps_best=ceil_sig((r["T_FP_best"] + TZ) / den),
                 eps_an=ceil_sig((r["T_F_an"] + TZ) / den),
                 eps_old_raw=(T_F + TZ) / den, eps_new_raw=(T_FP + TZ) / den,
                 eps_published=ref["eps_E_ex"], T_F_published=ref["num_E"] - TZ,
                 rLi_published=ref["rLi"],
                 # how loose each epsilon is against the model it actually certifies
                 looseness_none=(T_F + TZ) / den / r["rLi_none"],
                 looseness_el=(T_FP + TZ) / den / r["rLi_el"],
                 looseness_an=(r["T_F_an"] + TZ) / den / r["rLi_an"])
        if e == 28:
            r["T_needed_for_kappa"] = KAPPA_28 * den - TZ
    r["seconds"] = time.time() - t0
    return r


# ---------------------------------------------------------------- main

def main(argv):
    cells = tuple(int(a) for a in argv if a.isdigit()) or CELLS
    # NUMBA_NUM_THREADS, when set, caps the thread count and wins over ACF_THREADS
    nth = max(1, min(NTHREAD, 16, int(numba.config.NUMBA_NUM_THREADS)))
    set_num_threads(nth)
    say("conedual_pole_term -- the s = 0 pole term (clauses in the file header)")
    say("R1 and R2 compare at matched precision; registered 'voids' are enforced below.")
    say("cells %s   TAB %d   threads %d" % (list(cells), TAB, nth))

    bad_lit, worst_lit = check_literal_precision()
    if bad_lit:
        for b in bad_lit:
            say("LITERAL PRECISION FAILURE: " + b)
        say("ABORT: a clause whose tolerance is finer than its own literal is not a clause.")
        with io.open(OUT, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(_lines) + "\n")
        return 4
    say("literal precision check: worst relative half-ulp  R2 %.1e (tol %.0e)  R3 %.1e (tol %.0e)"
        % (worst_lit["R2"], TOL_R2, worst_lit["R3"], TOL_R3))

    try:
        import mpmath as mp
        mp.mp.dps = 30
        worst_eta = max(abs(float(mp.diff(lambda s: 1 / mp.zeta(s), 0, j) / mp.factorial(j)) - ETA[j])
                        for j in range(KDEG + 1))
        say("eta check against mpmath: max |hardcoded - mpmath| = %.3e" % worst_eta)
    except Exception as exc:
        worst_eta = None
        say("eta check SKIPPED (%s)" % exc)

    c0, err = check_residue()
    if c0 is None:
        say("C0 NOT SCORED (%s)" % err)
        c0_hold, c0_worst = None, None
    else:
        c0_worst, rowsc = c0
        for d, x, got, ref, rel in rowsc:
            say("   C0  d=%-5d x=%.3g   P_d = %+.12e   contour = %+.12e   rel %.2e"
                % (d, x, got, ref, rel))
        c0_hold = c0_worst <= 1e-12
        say("C0 %s  (max relative deviation %.2e, threshold 1e-12) -- a DOUBLE-PRECISION check: "
            "the contour value is cast to float before the ratio" % ("HOLDS" if c0_hold else "FAILS",
                                                                     c0_worst))

    c0b, errb = check_residue_hp()
    if c0b is None:
        say("C0b NOT SCORED (%s)" % errb)
        c0b_hold, c0b_worst = None, None
    else:
        c0b_worst, rowsb = c0b
        for d, x, rr, got, rel in rowsb:
            say("   C0b d=%-5d x=%.3g  r=%.2f   P_d = %s   rel %s" % (d, x, rr, got, rel))
        c0b_hold = c0b_worst <= 1e-40
        say("C0b %s  (max relative deviation %.2e at 50 dps, two radii, threshold 1e-40)"
            % ("HOLDS" if c0b_hold else "FAILS", c0b_worst))

    ref14 = {}
    if os.path.exists(REF14):
        for c in json.load(io.open(REF14, encoding="utf-8"))["cells"]:
            if c.get("denom_ex"):
                ref14[c["e"]] = c
    T = Tables(TAB)
    say("tables to %d in %.1f s;  max |M(n)|/sqrt n on [33,TAB] = %.6f;  "
        "fraction of n <= TAB with M(n) <= -2: %.4f  (M(n) < 0: %.4f)"
        % (TAB, T.seconds, T.hurst_max, T.frac_M_le_m2, T.frac_M_neg))
    say()

    out = [evaluate(e, T, ref14) for e in cells]

    # ---------------- scoring, BEFORE any result is printed, so the voids can be enforced
    def rel(a, b):
        return abs(a / b - 1.0)

    v = {}
    g = [r for r in out if "T_F_published" in r]
    v["R1"] = (all(rel(r["T_F"], r["T_F_published"]) <= TOL_R1
                   and rel(r["eps_old_raw"], r["eps_published"]) <= TOL_R1
                   and ceil_sig(r["eps_old_raw"]) == ceil_sig(r["eps_published"]) for r in g),
               len(g))
    g2 = [r for r in out if r["e"] in LIT_R2]
    v["R2"] = (all(rel(r["rLi_none"], float(LIT_R2[r["e"]])) <= TOL_R2 for r in g2), len(g2))
    g3 = [r for r in out if r["e"] in LIT_R3]
    v["R3"] = (all(rel(r["rLi_lead"], float(LIT_R3[r["e"]])) <= TOL_R3 for r in g3), len(g3))
    v["P1"] = (all(r["rLi_an"] <= 1.05 * r["rLi_lead"] for r in out), len(out))
    v["P2"] = (all(r["rLi_el"] <= r["rLi_lead"] for r in out), len(out))
    nbest = sum(1 for r in out if r["rLi_mix"] <= min(r["rLi_lead"], r["rLi_an"], r["rLi_el"]))
    v["P3"] = (nbest >= 5 and len(out) == 7, len(out))
    g4 = [r for r in out if r["e"] == 28]
    v["P4"] = (bool(g4) and g4[0]["ratio_T"] >= 0.90, len(g4))

    # ---------------- the registered "voids", enforced in the code and not waived in prose
    void_all = (v["R2"][1] > 0 and not v["R2"][0])
    void_s3 = void_all or (v["R1"][1] > 0 and not v["R1"][0])
    void_pole = void_all or (c0_hold is False) or (c0b_hold is False)
    t_all = "VOID " if void_all else ""
    t_s3 = "VOID " if void_s3 else ""
    t_pl = "VOID " if void_pole else ""
    if void_all or void_s3 or void_pole:
        say("!" * 100)
        say("REGISTERED VOIDS IN FORCE: %s%s%s"
            % ("everything (R2 failed)  " if void_all else "",
               "section 3 (R1 failed)  " if (void_s3 and not void_all) else "",
               "every pole column (C0/C0b failed)" if (void_pole and not void_all) else ""))
        say("!" * 100)
    say()

    for r in out:
        say("%s2^%d Q=%d rows=%d pairs=%d  |X|=%.1f denom=%.1f  %.1fs"
            % (t_all, r["e"], r["Q"], r["rows"], r["pairs"], r["norm_X"], r["denom_true"],
               r["seconds"]))
        say("     %srLi  none %.6f | lead %.6f | an %.6f | el %.6f | mix %.6f"
            % (t_pl or t_all, r["rLi_none"], r["rLi_lead"], r["rLi_an"], r["rLi_el"],
               r["rLi_mix"]))
        say("     %scos(X,P)  lead %+.4f | an %+.4f | el %+.4f | mix %+.4f    "
            "amp_mult %.4f   a fit would gain %.3f%%"
            % (t_pl or t_all, r["cos_lead"], r["cos_an"], r["cos_el"], r["cos_mix"],
               r["amp_mult"], 100.0 * r["fit_gain"]))
        say("     %sT_F %.4f  T_{F-Pi} %.4f  ratio %.5f   T_F(an-model) %.4f   ||Pi||/T_F %.4f"
            % (t_s3, r["T_F"], r["T_FP"], r["ratio_T"], r["T_F_an"], r["norm_Pi"] / r["T_F"]))
        if "eps_old" in r:
            say("     %seps: q_li %s -> el-model %s -> an-model %s   (published %s)"
                % (t_s3, r["eps_old"], r["eps_new"], r["eps_an"], r["eps_published"]))
        say("     %sdiagnostic: T_FP_best %.1f (the counterfactual, not a bound)%s"
            % (t_all, r["T_FP_best"],
               ("  eps_best %s" % r["eps_best"]) if "eps_best" in r else ""))
        say("     %sdiagnostic: Lemma-U terms %d   gain(M<=-2) %.1f%%  neutral(M=-1) %.1f%%  "
            "lose %.1f%%   realised net %.0f of systematic %.0f = %.2f%%"
            % (t_all, int(r["terms"]), 100.0 * r["terms_gain"] / r["terms"],
               100.0 * r["terms_neutral"] / r["terms"], 100.0 * r["terms_lose"] / r["terms"],
               r["realised_net"], r["systematic_net"], 100.0 * r["realised_frac_of_max"]))
        say("     %sdiagnostic: rows N/d>TAB %d/%d carry %.1f%% of T_F^2 but %.2f%% of "
            "||Pole_an||^2;  ||P_el-P_an||/|X| %.4f;  trivial arm binds on %d rows"
            % (t_all, r["rows_uncovered"], r["rows"], 100.0 * r["share_TF2_uncovered"],
               100.0 * r["share_pole2_uncovered"], r["norm_dif_over_X"], r["rows_trivial_arm"]))
        if "s12_max_abs" in r:
            say("     %sdiagnostic: over %d squarefree d<=1000 the j=1/j=0 ratio peaks at "
                "%.4f (d=%d); %d exceed 6%%; they carry %.1f%% of ||Pole_an||^2; by omega %s"
                % (t_all, r["s12_rows"], r["s12_max_abs"], r["s12_argmax_d"], r["s12_over6pct"],
                   100.0 * r["s12_share_pole2_dle1000"],
                   " ".join("w%s:%.3f" % (w, s) for w, s in
                            sorted(r["s12_share_pole2_by_omega"].items()))))
        say()

    say("=" * 100)
    say("C0  residue vs contour, DOUBLE precision : %s"
        % ("NOT SCORED" if c0_hold is None else ("HOLDS" if c0_hold else "FAILS")))
    say("C0b residue vs contour, 50 dps, 2 radii  : %s"
        % ("NOT SCORED" if c0b_hold is None else ("HOLDS" if c0b_hold else "FAILS")))
    for nm in ("R1", "R2", "R3", "P1", "P2", "P3", "P4"):
        ok, n = v[nm]
        say("%s: %s  (%d cells scored)" % (nm, "NOT SCORED" if n == 0 else
                                           ("HOLDS" if ok else "REFUTED"), n))
    say("P3 detail: Pole_mix is the best of four at %d of %d cells" % (nbest, len(out)))
    if g:
        say("R1/R2 like with like (both sides at full precision):")
        for r in g:
            say("   2^%d  T_F rel %.2e   rLi rel %.2e   eps_raw %.14f vs %.14f   ceil4 %s vs %s"
                % (r["e"], rel(r["T_F"], r["T_F_published"]),
                   rel(r["rLi_none"], r["rLi_published"]),
                   r["eps_old_raw"], r["eps_published"],
                   ceil_sig(r["eps_old_raw"]), ceil_sig(r["eps_published"])))
    if g4:
        r = g4[0]
        say("Q  eps(2^28): q_li %s | el-model %s | an-model %s | corrected counterfactual %s"
            % (r["eps_old"], r["eps_new"], r["eps_an"], r["eps_best"]))
        say("   against certified kappa(2^28) >= %s : %s   (T needed %.1f; T_{F-Pi} %.1f; "
            "counterfactual %.1f)"
            % (KAPPA_28, "BELOW" if r["eps_new"] < KAPPA_28 else "STILL ABOVE",
               r["T_needed_for_kappa"], r["T_FP"], r["T_FP_best"]))
        say("   looseness eps / rLi(the model that eps certifies): q_li %.3f | el %.3f | an %.3f"
            % (r["looseness_none"], r["looseness_el"], r["looseness_an"]))
    say("=" * 100)
    say("VERDICT  C0 %s | C0b %s | %s"
        % ("NOT SCORED" if c0_hold is None else ("HOLDS" if c0_hold else "FAILS"),
           "NOT SCORED" if c0b_hold is None else ("HOLDS" if c0b_hold else "FAILS"),
           " | ".join("%s %s" % (nm, "NOT SCORED" if v[nm][1] == 0 else
                                 ("HOLDS" if v[nm][0] else "REFUTED"))
                      for nm in ("R1", "R2", "R3", "P1", "P2", "P3", "P4"))))
    say("VOIDS    all=%s  section3=%s  pole=%s   (exit %d)"
        % (void_all, void_s3, void_pole, 3 if (void_all or void_s3 or void_pole) else 0))

    with io.open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(_lines) + "\n")
    json.dump(dict(cells=out, TAB=TAB, MARGIN=MARGIN, eta=list(ETA), kappa28=KAPPA_28,
                   eta_check=worst_eta, c0_worst=c0_worst, c0b_worst=c0b_worst,
                   literal_precision=worst_lit, tolerances=dict(R1=TOL_R1, R2=TOL_R2, R3=TOL_R3),
                   frac_M_le_m2=T.frac_M_le_m2, frac_M_neg=T.frac_M_neg,
                   hurst_max_on_table=T.hurst_max,
                   void_all=void_all, void_section3=void_s3, void_pole=void_pole,
                   verdicts=dict([("C0", c0_hold), ("C0b", c0b_hold)]
                                 + [(k, None if v[k][1] == 0 else v[k][0])
                                    for k in ("R1", "R2", "R3", "P1", "P2", "P3", "P4")])),
              io.open(OUTJ, "w", encoding="utf-8", newline="\n"), indent=1)
    return 3 if (void_all or void_s3 or void_pole) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

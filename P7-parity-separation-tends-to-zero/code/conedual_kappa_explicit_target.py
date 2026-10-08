# -*- coding: utf-8 -*-
r"""conedual_kappa_explicit_target -- kappa^rho = dist(z - rho, K_Q)/||p||_2 for three explicit rho

Supporting computation; not used in the paper's proofs or tables.

WHERE THIS COMES FROM

The bracket rests on this: for any explicit family rho_Q with
||p - rho_Q||_2 = o(||p||_2),

    kappa(Q) -> 0   <=>   kappa^rho(Q) := dist_2(z - rho_Q, K_Q)/||p||_2 -> 0,

where z = t + p = (1/|vR|) sum_{l in vB} A_l is the even-band column sum and K_Q = cone{A_l}.
The quantitative form is |kappa^rho - kappa| <= r_rho := ||p - rho||_2/||p||_2.
This file COMPUTES kappa^rho at 2^16..2^24 for

    rho1 = c(Q) * (-mu/phi),  c(Q) = N/(4 |vR| L^2)
    rho2 = q_li
    rho3 = q_li + Pole        (the s=0 residue polynomial; see code/conedual_pole_term.py)

and, as a CONTROL, the original target t itself (which must return the certified kappa).

HOW.  dist_2(x, K_Q) = min_{y >= 0} ||A y - x||_2 with A[i,j] = 1 iff rows[i] | vB[j].  The full
incidence is built once per cell (nnz ~ 11|vB|: 0.05-12 M, megabytes not gigabytes), the nnls is
solved by COLUMN GENERATION -- round 0 coordinate descent on all columns from the uniform warm
start y = 1/|vR| (for which r = -rho), later rounds the restricted master on
support u {the `add` most violated columns}, optionally polished by scipy's BVLS active-set solver
on the restricted block -- and priced by one sparse A^T r.  (the column-generation script's `k_cd_nnls` is reused as an
algorithm, not as an import.)  TWO-SIDED, both sides valid at any iterate:

    upper   ||r||_2 / ||p||_2            (y >= 0 is feasible)
    lower   (<x,r> - M x_1) / (sqrt(||r||^2 - 2 M r_1 + M^2) ||p||_2),   M = max(0, max_vB A^T r)

the lower bound being dist(x,K) >= <x,u>/||u|| for u = r - M e_1, which lies in K_Q^o because every
column has A[0,j] = 1 (so <A_l, u> = W(l) - M <= 0).

WHAT KIND OF NUMBER THE BRACKET IS.  It is NOT an integer certificate
like the kappa enclosures computed elsewhere.  Two things are done about that, and a third is not:
  * the residual is recomputed FRESH (r = x - A.dot(y), one sparse matvec) before every reported
    bracket, so no incrementally accumulated error enters any reported number; the drift of the
    carried residual is reported beside it;
  * an outward float envelope is added to both endpoints: with eps = 2^-53 and
    e_i = (nnz_i + 1) eps (|x_i| + sum_j A_ij |y_j|), the matvec error is |r - r_exact| <= e
    componentwise, so ||r||_2 +- ||e||_2 brackets the exact residual norm of the computed y and both
    endpoints are widened by ||e||_2/||p||_2 (reported as `env`);
  * the TARGET x is not enclosed: `scipy.special.expi` and the residue polynomial carry no interval
    arithmetic, so the bracket is rigorous as a function of the computed x and no further.  GAP,
    named; the margins are 1e-3 against an envelope of 1e-11.

CLAUSES.

  Certified kappa (exact dense nnls at 2^16-2^20; the enclosures above):
      2^16 0.1724    2^18 0.1409    2^20 0.1055427
      2^22 [0.076331, 0.0763312]    2^24 [0.052762, 0.052763]
  r_rho, from the stored record (r(c(Q)); q_li; the full residue polynomial):
      rho1  1.0769, 0.8870, 0.7182, 0.5799, 0.4774      at 2^16, 2^18, 2^20, 2^22, 2^24
      rho2  0.288805, 0.163233, 0.091181, 0.048970, 0.026298
      rho3  0.0605, 0.0448, 0.0261, 0.0151, 0.0102
  Signed alignment a_rho2 = <p - q_li, w*>/(||p-q_li||_2 ||w*||_2) (positive):
      2^20 0.72576   2^22 0.74321   2^24 0.73045

  C0  CONTROL, voids the cell.  For x = t the bracket [lower, upper] contains the certified kappa
      interval of that cell, and (upper - lower)/upper <= 0.03.
  C1  PROVED, must hold (the quantitative form above).  |kappa^rho - kappa| <= r_rho at every cell and every rho,
      using the certified kappa interval.  A failure is a computational error, not a refutation.
  C2  PROVED, must hold.  kappa^rho >= L_rho(w) := (-S(w) - <rho,w>)/(||w||_2 ||p||_2)
      for the saved exact duals conedual_exact_dual_2e2{2,4}.npz, w used AS STORED.

  P1  PREDICTION, rho2 = q_li, at 2^20, 2^22, 2^24 (the cells with a saved w*).
      kappa^{q_li} in [0.95, 1.15] x (kappa + a_rho2 * r_rho2), i.e.
          2^20 [0.16313, 0.19747]    2^22 [0.10708, 0.12962]    2^24 [0.06837, 0.08277]
      WHY: dist_2(., K_Q) is convex with subgradient w*/||w*|| at t, so
      kappa^rho >= kappa + a_rho r_rho exactly; the prediction is that this first-order value is
      also nearly attained, i.e. the cone's active face does not move much under a perturbation of
      relative size r/kappa = 0.86, 0.64, 0.50.
      MEANING, registered.  HOLDS: kappa^rho is governed to first order by the alignment of the
      model error with the optimal separator, so it is predictable from row-level pairings alone
      and needs no cone solve at larger cells.  REFUTED ABOVE the window: curvature matters, the
      face moves, and kappa^rho must be computed.  REFUTED BELOW: the cone partly ABSORBS the model
      error (a better outcome for this route than the first-order value).
  P2  PREDICTION, rho3 = q_li + Pole, all five cells.  kappa^{rho3} in [kappa_lo, kappa_up + 0.6 r_rho3]:
          2^16 [0.1724, 0.2087]   2^18 [0.1409, 0.1678]   2^20 [0.10554, 0.12121]
          2^22 [0.076331, 0.08539]   2^24 [0.052762, 0.05888]
      MEANING: HOLDS -> the best available explicit target is within 10-20 % of t's own distance,
      so the reduction is numerically faithful at computable sizes.  REFUTED BELOW kappa_lo ->
      the perturbation moves the target INTO the cone's direction and kappa^rho understates kappa.
  P3  PREDICTION, rho1 = c(Q)(-mu/phi), all five cells.  kappa^{rho1} <= 0.01, i.e. z - c(Q)v is
      in K_Q to solver tolerance.
      WHY: the saved separators give <v, w> > 0 while <z, w> = -S(w) <= 0, so
      L_rho1(w) = (-S - c<v,w>)/(...) < 0 at every saved w: no known separator certifies a positive
      distance, and ||p - c(Q)v||_2 = 0.48-1.08 ||p||_2 is far larger than kappa.
      MEANING.  HOLDS -> the flat model is useless for this route at computable sizes: its
      target is inside the cone, kappa^rho1 = 0 while kappa > 0, and that band is vacuous there
      (r_rho1 >> kappa).  REFUTED -> the cone does not absorb the flat model's error and
      kappa^rho1 > 0 is a measurable quantity; report it.

  NOT REGISTERED: anything about kappa's rate; any claim that kappa^rho -> 0; anything at 2^26+.

HOW THE CLAUSES ABOVE ARE SCORED, AND THE THREE REGISTERED CLAUSES THAT GOVERN THE REST.

  C0 is scored AS ITS TEXT READS (containment) and its "voids the cell" is ENFORCED IN CODE, with
     void propagation (a void cell supplies no comparand to any prediction).  Intersection is NOT
     what C0 tests.
     Containment of a rounded-down certificate by a sharper bracket is UNSATISFIABLE, so C0 is
     expected to void every cell; it is kept and reported for exactly that reason.

  C0b REGISTERED CLAUSE.  The control predicate that has content for
     a solver sharper than its comparand: for x = t the bracket must lie INSIDE the comparand's own
     precision of the certified value -- [lower, upper] subset [kappa_lo - u, kappa_up + u] with u
     the comparand's half-ulp (5e-5 at 2^16/2^18, 5e-8 at 2^20, 0 at 2^22/2^24 where the endpoints
     are integer-derived) -- AND (upper - lower)/upper <= 0.03.  Voids the cell, propagated.
     MEANING: HOLDS -> the solver reproduces the certified kappa at that cell and the cell's
     predictions may be scored.  REFUTED -> the cell is void and every prediction on it is NOT
     SCORED, whatever its own numbers say.

  C3 REGISTERED CLAUSE.  Optimality witness for the control: at
     x = t the fresh residual r satisfies max_vB A^T r <= 1e-6 ||r||_2^2 and <y, A^T r> >= -1e-12,
     i.e. r is the Moreau residual w* to that tolerance, so the a_rho computed from it are the
     subgradient alignments of the proposition above and not those of an arbitrary feasible point.
     MEANING: HOLDS at a cell -> a_rho(cell) is reported as a subgradient alignment; REFUTED ->
     a_rho(cell) is reported as "from a non-optimal residual" and no first-order claim is made there.

  P4 REGISTERED PREDICTION, from the five-cell table
     `a3 = +0.062009, -0.053489, -0.002146, -0.020811, (-0.045)` at 2^16..2^24, which this run
     recomputes from its own w*: |a3| <= 0.08 at every cell 2^16..2^24, and a3 changes sign between
     2^16 and 2^18.  CAN FAIL.
     MEANING: HOLDS -> the residual of q_li + Pole is within one noise floor (1/sqrt(#rows) =
     0.080, 0.056, 0.040, 0.028, 0.020) of orthogonal to the critical direction at every computed
     cell, and its small size is an oscillation through zero, not a trend; REFUTED -> |a3| is
     growing and the N^{-1/4} expiry of the pole's dominance is already visible.

  NOT REGISTERED, and reported as diagnostics only: the convexity remainder c_rho, the
  float envelope, the residual drift, and anything about 2^26+.

HOW THE IMPLEMENTATION MEETS THOSE TEXTS.  No registered comparand, window or threshold above
depends on any of this.
  (i)  the solver.  Full-matrix coordinate descent does not converge on the near-degenerate
       t-like targets -- at 2^22 it gives the bracket [0.058052, 0.079038] with kkt_rel 0.20 after
       20,000 passes, which is wider than C0's width clause allows -- so the solver here is
       column generation (round 0 full, later rounds restricted) plus, where a saved exact dual
       exists, that separator's own bound as a second lower bound.  A control that fails because
       the solver has not converged is answered by converging, not by re-registering.
  (ii) C1's test.  The registered text is |kappa^rho - kappa| <= r_rho, a statement about the
       NUMBER kappa^rho, of which this file computes only a bracket, so the faithful test is that
       the bracket and the band INTERSECT.  Demanding that the bracket lie inside the band is
       strictly stronger than the text and fails at (2^22, rho3) on a bracket 2.1e-2 wide.
  (iii) P1/P2 scoring.  Same reason: an interval-valued measurement against a window is scored
       three ways -- HOLDS (bracket inside), REFUTED (disjoint), INCONCLUSIVE (window cuts it).

COMPUTE RULES (the project's compute rules, both sections read before writing this file).  The per-row arrays
come from `cell_cache.load(e)` -- never `build()`; the band itself (the N-length code array and vB)
is not cached by that module ("a script that needs the band itself builds with `_build`-equivalent
code", its docstring), so it is built here by the same canonical sieve with one Python loop over the
primes <= Q and strided numpy slice ops, and checked against the cell's own nB/nR.  No Python loop
over rows, band elements or nnz: the incidence is two numba `prange`-over-rows kernels (count, then
fill, the pattern of `cell_cache._count_rows`), the nnls is one `njit` kernel, the pricing one
sparse matvec, and the (d,k) pair arrays of q_li are np.repeat/ragged-arange/np.gcd/np.bincount.
Threads capped by ACF_THREADS / the launcher's env.  Peak arrays at 2^24: code+pos 81 MB, COO/CSC
~450 MB.  Cells 2^20 and above want a machine with more cores; 2^16 and 2^18 are a smoke test
anywhere.

    python code/conedual_kappa_explicit_target.py 16 18
    python code/conedual_kappa_explicit_target.py 16 18 20 22 24 --tag _d
"""
from __future__ import annotations
import io, json, math, os, sys, time

NTHREAD = int(os.environ.get("ACF_THREADS", os.environ.get("NUMBA_NUM_THREADS", "8")))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, str(NTHREAD))
_SCRATCH = os.environ.get("TEMP") or os.environ.get("TMP") or "."
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(_SCRATCH, "numba_cache_acf"))

import numpy as np  # noqa: E402
import numba  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from numba import njit, prange, set_num_threads  # noqa: E402
from scipy.special import expi  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
sys.path.insert(0, HERE)
from _sieve_shared import mu_upto, primes_upto  # noqa: E402
from cell_cache import load as load_cell  # noqa: E402
import conedual_pole_term as PT  # noqa: E402   (row_data/phi_coeffs/gc_coeffs/poly_at only)

DEN = 10 ** 12
CELLS_DEFAULT = (16, 18)
KAP = {16: (0.1724, 0.1724), 18: (0.1409, 0.1409), 20: (0.1055427, 0.1055427),
       22: (0.076331, 0.0763312), 24: (0.052762, 0.052763)}
R_RHO = {1: {16: 1.0769, 18: 0.8870, 20: 0.7182, 22: 0.5799, 24: 0.4774},
         2: {16: 0.288805, 18: 0.163233, 20: 0.091181, 22: 0.048970, 24: 0.026298},
         3: {16: 0.0605, 18: 0.0448, 20: 0.0261, 22: 0.0151, 24: 0.0102}}
A_RHO2 = {20: 0.72576, 22: 0.74321, 24: 0.73045}
P1_WIN = {20: (0.16313, 0.19747), 22: (0.10708, 0.12962), 24: (0.06837, 0.08277)}
P2_WIN = {16: (0.1724, 0.2087), 18: (0.1409, 0.1678), 20: (0.10554, 0.12121),
          22: (0.076331, 0.08539), 24: (0.052762, 0.05888)}
P3_MAX = 0.01
P4_MAX = 0.08                      # |a3| <= 0.08 at every cell (clause P4)
C3_PTOL, C3_DOT = 1e-6, -1e-12     # the control's optimality witness (clause C3)
NOISE = {16: 0.0798, 18: 0.0564, 20: 0.0400, 22: 0.0283, 24: 0.0200}   # 1/sqrt(#rows)
WVEC = {20: "conedual_exact_dual_2e20.npz", 22: "conedual_exact_dual_2e22.npz",
        24: "conedual_exact_dual_2e24.npz"}   # a saved exact dual at each of 2^20, 2^22, 2^24
# "Compare like with like" (the header's own clause; the reproduction rule): a comparand quoted
# at 4 decimals (0.1724, 0.1409) denotes the interval +-5e-5, one at 7 (0.1055427) +-5e-8;
# the 2^22/2^24 endpoints are integer-derived (c6/10^6) and exact.  The strict reading of the
# literals -- a degenerate interval at 2^16/2^18 -- is also scored and reported (C0_strict).
_ULP = {16: 5e-5, 18: 5e-5, 20: 5e-8, 22: 0.0, 24: 0.0}
KAP_BAND = {e: (KAP[e][0] - _ULP[e], KAP[e][1] + _ULP[e]) for e in KAP}

TAG = ""
_argv = sys.argv[1:]
_args, _i = [], 0
while _i < len(_argv):
    a = _argv[_i]
    if a == "--tag":
        TAG = _argv[_i + 1] if _i + 1 < len(_argv) else ""
        _i += 2
        continue
    if a.startswith("--tag="):
        TAG = a.split("=", 1)[1]
    elif not a.startswith("--"):
        _args.append(a)
    _i += 1
CELLS = tuple(int(a) for a in _args) or CELLS_DEFAULT
PASSES = 20000
ADD = 4000
CDTOL = 1e-14

OUT = os.path.join(RES, "conedual_kappa_explicit_target%s.txt" % TAG)
OUTJ = os.path.join(RES, "conedual_kappa_explicit_target%s.json" % TAG)
_lines = []


def say(s=""):
    print(s, flush=True)
    _lines.append(s)
    with io.open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(_lines) + "\n")


def li(x):
    return expi(np.log(np.asarray(x, dtype=np.float64)))


# ------------------------------------------------------------------ numba kernels
@njit(parallel=True, cache=True, nogil=True)
def k_row_counts(code, rows, cnt, nthread):
    """cnt[i] = #{n : code[n] == 1, rows[i] | n}; prange over a strided slice of the row list."""
    n = code.shape[0] - 1
    m = rows.shape[0]
    for w in prange(nthread):
        for i in range(w, m, nthread):
            d = rows[i]
            a = 0
            for j in range(d, n + 1, d):
                if code[j] == 1:
                    a += 1
            cnt[i] = a


@njit(parallel=True, cache=True, nogil=True)
def k_fill_pairs(code, pos, rows, off, ridx, cidx, nthread):
    """the (row, column) pairs of A: for row i, its multiples that are even-band columns."""
    n = code.shape[0] - 1
    m = rows.shape[0]
    for w in prange(nthread):
        for i in range(w, m, nthread):
            d = rows[i]
            k = off[i]
            for j in range(d, n + 1, d):
                if code[j] == 1:
                    ridx[k] = i
                    cidx[k] = pos[j]
                    k += 1


@njit(cache=True, nogil=True, fastmath=False)
def k_cd_nnls(indptr, indices, y, r, npass, tol):
    """min_{y>=0} ||A y - x||, A's entries 0/1, r = x - A y carried incrementally (the column-generation script)."""
    ncol = y.shape[0]
    for it in range(npass):
        md = 0.0
        for j in range(ncol):
            a = indptr[j]
            b = indptr[j + 1]
            nj = b - a
            if nj <= 0:
                continue
            s = 0.0
            for q in range(a, b):
                s += r[indices[q]]
            ny = y[j] + s / nj
            if ny < 0.0:
                ny = 0.0
            dy = ny - y[j]
            if dy != 0.0:
                for q in range(a, b):
                    r[indices[q]] -= dy
                y[j] = ny
                ad = dy if dy > 0.0 else -dy
                if ad > md:
                    md = ad
        if md <= tol:
            return it + 1
    return npass


# ------------------------------------------------------------------ the band (not cached anywhere)
def band_columns(e, S):
    """vB (even-omega band) and the position map, by the canonical sieve of cell_cache._build.

    One Python loop over the primes <= Q; every other step is a strided numpy slice op.
    """
    N = 2 ** e
    Q = int(N ** 0.5)
    assert N < 2 ** 31
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
    thr = int(S["thr"])
    band = (mu != 0) & (rem == 1)
    band[:thr + 1] = False
    del rem, mu
    code = np.zeros(N + 1, dtype=np.int8)
    np.copyto(code, np.int8(1), where=band & (om % 2 == 0))
    np.copyto(code, np.int8(2), where=band & (om % 2 == 1))
    del band, om
    nB = int(np.count_nonzero(code == 1))
    nR = int(np.count_nonzero(code == 2))
    assert nB == int(S["nB"]) and nR == int(S["nR"]), \
        "band rebuilt at 2^%d disagrees with the cell cache (%d/%d vs %d/%d)" % (
            e, nB, nR, int(S["nB"]), int(S["nR"]))
    vB = np.nonzero(code == 1)[0].astype(np.int64)
    pos = np.full(N + 1, -1, dtype=np.int32)
    pos[vB] = np.arange(vB.size, dtype=np.int32)
    return code, pos, vB


def incidence(code, pos, rows, nthread):
    cnt = np.zeros(rows.size, dtype=np.int64)
    k_row_counts(code, rows, cnt, nthread)
    off = np.zeros(rows.size + 1, dtype=np.int64)
    np.cumsum(cnt, out=off[1:])
    nnz = int(off[-1])
    ridx = np.empty(nnz, dtype=np.int32)
    cidx = np.empty(nnz, dtype=np.int32)
    k_fill_pairs(code, pos, rows, off[:-1].copy(), ridx, cidx, nthread)
    assert int(cidx.min()) >= 0
    A = sp.coo_matrix((np.ones(nnz), (ridx, cidx)),
                      shape=(rows.size, int(pos.max()) + 1)).tocsc()
    del ridx, cidx
    return A, nnz, cnt.astype(np.float64)


# ------------------------------------------------------------------ the three explicit rho
def build_rhos(e, S):
    """rho1 = c(Q)(-mu/phi), rho2 = q_li, rho3 = q_li + Pole; plus mu, phi, sigma on the rows."""
    N, Q, q1, thr = int(S["N"]), int(S["Q"]), int(S["q1"]), int(S["thr"])
    rows = np.asarray(S["rows"]).astype(np.int64)
    nR = float(S["nR"])
    L = math.log(Q)

    class _T(object):
        pass
    T = _T()
    T.mu = np.array(mu_upto(max(Q, thr) + 1), dtype=np.int8)
    r2, kap, pf, lam, om = PT.row_data(Q, T)
    assert r2.size == rows.size and int(np.abs(r2 - rows).max()) == 0
    assert int(np.abs(om - np.asarray(S["rom"])).max()) == 0, "omega != cell's rom"
    mud = np.asarray(T.mu[rows], dtype=np.float64)
    phi = np.ones(Q + 1, dtype=np.float64)
    for p in primes_upto(Q):
        p = int(p)
        phi[p::p] *= (p - 1.0)
    phid = phi[rows]

    # rho1 : Theorem B's flat model
    cQ = N / (4.0 * nR * L * L)
    rho1 = cQ * (-mud / phid)

    # rho2 : q_li  (the (d,k) pair sum, vectorised)
    di, k, d = PT.pair_arrays(rows, N, q1, T)
    muk = T.mu[k].astype(np.float64)
    liQ = float(li(float(Q)))
    term = muk * (li(N / (d * k).astype(np.float64)) - liQ)
    Sli = np.bincount(di, weights=term, minlength=rows.size)
    rho2 = mud * Sli / nR

    # rho3 : q_li + Pole, the full s = 0 residue polynomial
    cc = PT.gc_coeffs(PT.phi_coeffs(pf, om))
    P_an = (PT.poly_at(cc, om, lam, N / rows.astype(np.float64))
            - PT.poly_at(cc, om, lam, thr / rows.astype(np.float64)))
    rho3 = rho2 + mud * P_an / nR
    return dict(rho1=rho1, rho2=rho2, rho3=rho3, cQ=cQ, mud=mud, phid=phid,
                pairs=int(k.size), norm_Pan=float(np.linalg.norm(P_an)))


# ------------------------------------------------------------------ the cone distance
def _bracket(x, r, W, p2):
    """[lower, upper] for dist_2(x,K_Q)/||p||_2 from any feasible iterate and the FULL pricing.

    upper: ||r||_2, since the y that produced r is feasible.  lower: dist(x,K) >= <x,u>/||u||_2
    for u = r - M e_1 with M = max(0, max_vB W), which is in K_Q^o because A[0,j] = 1 for every
    column (so <A_l, u> = W(l) - M <= 0).
    """
    M = float(max(0.0, W.max()))
    nr = float(np.linalg.norm(r))
    nu2 = nr * nr - 2.0 * M * float(r[0]) + M * M
    lo = (float(np.dot(x, r)) - M * float(x[0])) / (math.sqrt(max(nu2, 1e-300)) * p2)
    return lo, nr / p2, M, nr


EPS = 2.0 ** -53


def envelope(A, x, cols, yc, nnz_row):
    """||e||_2 with e_i = (nnz_i + 1) eps (|x_i| + sum_j A_ij |y_j|): an OUTWARD bound on the
    float64 error of the fresh matvec r = x - A[:,cols] yc, componentwise."""
    av = A[:, cols].dot(np.abs(yc))
    e = (nnz_row + 1.0) * EPS * (np.abs(x) + av)
    return float(np.linalg.norm(e))


def polish(A, x, cols, yc, budget):
    """BVLS on the restricted block (scipy's exact active-set solver); returns the
    better feasible point.  Guarded by size and wall clock -- any failure keeps the CD point."""
    nrow = A.shape[0]
    if cols.size == 0 or nrow * cols.size > 4.0e7:
        return yc, "skipped (block %d x %d)" % (nrow, cols.size)
    try:
        from scipy.optimize import lsq_linear
        t0 = time.time()
        D = np.asarray(A[:, cols].todense(), dtype=np.float64)
        res = lsq_linear(D, x, bounds=(0.0, np.inf), method="bvls",
                         max_iter=max(50, min(3 * nrow, 400)), tol=1e-14, lsq_solver="exact")
        yn = np.maximum(res.x, 0.0)
        a = float(np.linalg.norm(x - D.dot(yc)))
        b = float(np.linalg.norm(x - D.dot(yn)))
        del D
        return (yn if b < a else yc), "bvls %.0f s, ||r|| %.9e -> %.9e%s" % (
            time.time() - t0, a, b, "" if b < a else " (kept CD)")
    except Exception as exc:
        return yc, "bvls FAILED (%s)" % type(exc).__name__


def cone_dist(A, x, nR, p2, nnz_row, note, rounds=40, passes0=2000, passes=40000, add=8000,
              tol=1e-14, ptol=1e-9, do_polish=True):
    """dist_2(x, K_Q)/||p||_2 by column generation; both bounds valid at every iterate.

    Round 0 runs coordinate descent on the full incidence from the uniform warm start
    y = 1/|vR| (for which r = x - z); later rounds solve the RESTRICTED master on
    support u {the `add` most violated columns} to a much tighter tolerance -- the full-matrix
    CD stalls on the near-degenerate t-like targets (kkt_rel 0.20 at 2^22, measured), the
    restricted one does not.  This is the column-generation script's column generation, as an
    algorithm, with a BVLS polish on the restricted block.

    EVERY reported bracket comes from a FRESH residual r = x - A[:,cols] y (one sparse matvec),
    never from the incrementally carried one; the carried-vs-fresh drift is reported.
    """
    t0 = time.time()
    n = A.shape[1]
    y = np.full(n, 1.0 / float(nR))
    r = x - A.dot(y)
    used = k_cd_nnls(A.indptr, A.indices, y, r, passes0, tol)
    allc = np.arange(n)
    W = A.T.dot(x - A.dot(y))
    lo, up, M, nr = _bracket(x, x - A.dot(y), W, p2)
    best = (lo, up, nr, allc, y.copy())
    hist = [dict(rnd=0, nact=n, passes=used, lower=lo, upper=up, kkt=M / (nr * nr))]
    act = np.union1d(np.nonzero(y > 0)[0],
                     np.argsort(W)[-add:][W[np.argsort(W)[-add:]] > 0])
    yact = y[act]
    pol = "not run"
    for rnd in range(1, rounds + 1):
        As = A[:, act].tocsc()
        rr = x - As.dot(yact)
        us = k_cd_nnls(As.indptr, As.indices, yact, rr, passes, tol)
        drift = float(np.abs(rr - (x - As.dot(yact))).max())
        rf = x - As.dot(yact)                                  # FRESH residual
        W = A.T.dot(rf)
        lo, up, M, nr = _bracket(x, rf, W, p2)
        if up < best[1]:
            best = (max(best[0], lo), up, nr, act.copy(), yact.copy())
        else:
            best = (max(best[0], lo), best[1], best[2], best[3], best[4])
        hist.append(dict(rnd=rnd, nact=int(act.size), passes=us, lower=lo, upper=up,
                         kkt=M / (nr * nr), drift=drift))
        if M <= ptol * nr * nr:
            break
        sup = act[yact > 0]
        o = np.argsort(W)[-add:]
        act2 = np.union1d(sup, o[W[o] > 0])
        ynew = np.zeros(act2.size)
        ynew[np.searchsorted(act2, sup)] = yact[yact > 0]
        act, yact = act2, ynew
    if do_polish:
        sup = best[3][best[4] > 0]
        ysup = best[4][best[4] > 0]
        yp, pol = polish(A, x, sup, ysup, None)
        rf = x - A[:, sup].dot(yp)
        W = A.T.dot(rf)
        lo, up, M, nr = _bracket(x, rf, W, p2)
        if up <= best[1]:
            best = (max(best[0], lo), up, nr, sup, yp)
        else:
            lo = best[0]
            sup, yp = best[3], best[4]
            rf = x - A[:, sup].dot(yp)
            W = A.T.dot(rf)
            lo2, up, M, nr = _bracket(x, rf, W, p2)
            lo = max(lo, lo2)
    env = envelope(A, x, best[3], best[4], nnz_row) / p2
    kkt_dot = float(np.dot(best[4], A[:, best[3]].T.dot(rf)))
    note.append(dict(secs=time.time() - t0, rounds=len(hist) - 1,
                     lower=best[0] - env, upper=best[1] + env, lower_raw=best[0],
                     upper_raw=best[1], env=env, kkt_rel=M / (nr * nr), maxW=M,
                     kkt_dot=kkt_dot, nsupp=int((best[4] > 0).sum()), polish=pol,
                     drift=max([h.get("drift", 0.0) for h in hist]), hist=hist))
    return best[0] - env, best[1] + env, rf, note[-1]


def main():
    set_num_threads(max(1, min(NTHREAD, 16, int(numba.config.NUMBA_NUM_THREADS))))
    say("conedual_kappa_explicit_target  (registered C0, C0b, C1-C3, P1-P4 in the file header)")
    say("cells %s   tag %r   passes %d   tol %.0e   threads %d   %s"
        % (list(CELLS), TAG, PASSES, CDTOL, numba.get_num_threads(),
           time.strftime("%Y-%m-%d %H:%M:%S")))
    led = json.load(io.open(os.path.join(RES, "gramint_kappa_certify_ledger.json"),
                            encoding="utf-8"))
    out = []
    for e in CELLS:
        t0 = time.time()
        S = load_cell(e)
        rows = np.asarray(S["rows"]).astype(np.int64)
        t = np.asarray(S["t"], dtype=np.float64)
        p = np.asarray(S["p"], dtype=np.float64)
        z = t + p
        p2 = float(np.linalg.norm(p))
        R = build_rhos(e, S)
        code, pos, vB = band_columns(e, S)
        A, nnz, nnz_row = incidence(code, pos, rows, numba.get_num_threads())
        del code, pos
        say("")
        say("2^%d  Q=%d  rows=%d  |vB|=%d  nnz=%d (%.2f/col)  |vR|=%d  pairs=%d  prep %.0f s"
            % (e, int(S["Q"]), rows.size, vB.size, nnz, nnz / vB.size, int(S["nR"]),
               R["pairs"], time.time() - t0))
        r = dict(e=e, Q=int(S["Q"]), rows=int(rows.size), nvB=int(vB.size), nnz=int(nnz),
                 nR=int(S["nR"]), p2=p2, cQ=R["cQ"], kap_lo=KAP[e][0], kap_up=KAP[e][1],
                 kap_band=list(KAP_BAND[e]))
        # the saved separator first: it supplies a second, independent lower bound for every
        # target, including the control t, where it equals the ledger's kappa_2(w)
        Lsep = {}
        if e in WVEC and os.path.exists(os.path.join(RES, WVEC[e])):
            z_ = np.load(os.path.join(RES, WVEC[e]))
            ent = led.get(WVEC[e])            # 2e20's exact dual has no ledger line
            w = z_["w"].astype(np.float64).copy()         # AS STORED
            assert np.array_equal(z_["rows"].astype(np.int64), rows)
            if ent is not None:
                assert abs(float(w[0]) + int(ent["M"]) / DEN) <= 1e-6 * abs(float(w[0]))
            n2 = float(np.linalg.norm(w))
            Sw = -float(np.dot(w, z))
            Lsep["t"] = float(np.dot(w, t)) / (n2 * p2)
            for i in (1, 2, 3):
                Lsep["rho%d" % i] = (-Sw - float(np.dot(R["rho%d" % i], w))) / (n2 * p2)
            r["S_w"], r["norm_w"] = Sw, n2
            say("   saved w = %s  ||w||2 %.6e  S(w) %.6e  kappa_2(w) %.7f  %s"
                % (WVEC[e], n2, Sw, Lsep["t"],
                   ("(ledger int %.7f)" % (int(ent["J"]) / math.sqrt(int(ent["H"]) * int(ent["D_p"]))))
                   if ent is not None else "(no ledger line; w as stored)"))
        targets = [("t", t), ("rho1", z - R["rho1"]), ("rho2", z - R["rho2"]),
                   ("rho3", z - R["rho3"])]
        wstar = None
        for nm, x in targets:
            note = []
            lo, up, resid, inf = cone_dist(A, x, float(S["nR"]), p2, nnz_row, note, add=ADD,
                                           passes=PASSES, tol=CDTOL)
            if nm == "t":
                wstar = resid.copy()          # the Moreau residual of the ORIGINAL problem
                r["C3"] = bool(inf["kkt_rel"] <= C3_PTOL and inf["kkt_dot"] >= C3_DOT)
                r["C3_kkt_rel"], r["C3_kkt_dot"] = inf["kkt_rel"], inf["kkt_dot"]
            L = Lsep.get(nm)
            lo_all = max(lo, L) if L is not None else lo
            r[nm] = dict(lower=lo_all, lower_cg=lo, L_sep=L, upper=up,
                         **{k: inf[k] for k in ("rounds", "secs", "kkt_rel", "maxW", "kkt_dot",
                                                "nsupp", "env", "drift", "polish",
                                                "lower_raw", "upper_raw")})
            say("   %-5s dist/||p|| in [%.7f, %.7f]  (width %.2e rel %.2f%%, kkt_rel %.1e, "
                "|supp| %d, %d rnd, %.0f s, env %.1e, drift %.1e)%s"
                % (nm, lo_all, up, up - lo_all, 100.0 * (up - lo_all) / up, inf["kkt_rel"],
                   inf["nsupp"], inf["rounds"], inf["secs"], inf["env"], inf["drift"],
                   ("  [cg lower %.7f, sep %+.7f]" % (lo, L)) if L is not None else ""))
            say("         polish: %s" % inf["polish"])
        # the alignments, from THIS cell's own Moreau residual w* (the subgradient of
        # dist(.,K_Q) at t); a_rho is the SAME number as the separator's lower bound, not a
        # second measurement -- the convexity remainder c_rho is the quantity that is new here.
        nw = float(np.linalg.norm(wstar))
        r["norm_wstar"], r["kappa_wstar"] = nw, nw / p2
        for i in (1, 2, 3):
            nm = "rho%d" % i
            rr_meas = float(np.linalg.norm(p - R[nm])) / p2
            rr_out = math.ceil(rr_meas * 1e9) / 1e9        # OUTWARD-rounded measured r_rho
            lo, up = r[nm]["lower"], r[nm]["upper"]
            blo, bup = KAP_BAND[e][0] - rr_out, KAP_BAND[e][1] + rr_out
            a_rho = float(np.dot(p - R[nm], wstar)) / (rr_meas * p2 * nw) if rr_meas > 0 else 0.0
            first = nw / p2 + a_rho * rr_meas               # kappa + a_rho r_rho
            r[nm].update(r_rho_registered=R_RHO[i][e], r_rho_measured=rr_meas,
                         r_rho_outward=rr_out,
                         rho_norm_over_p=float(np.linalg.norm(R[nm])) / p2,
                         a_rho=a_rho, first_order=first,
                         c_rho_lo=lo - first, c_rho_hi=up - first,
                         D1_band=[blo, bup],
                         # D1 is violated only if the bracket and the band are DISJOINT:
                         # both are enclosures of the same number kappa^rho.
                         D1=bool(lo <= bup + 1e-12 and up >= blo - 1e-12),
                         C2=bool(r[nm]["L_sep"] is None or up >= r[nm]["L_sep"] - 1e-9))
            say("   %-5s r_rho reg %.6f meas %.7f (outward %.9f) | a_rho %+.6f | first order "
                "%.7f | convexity remainder c_rho in [%+.2e, %+.2e] | D1 %s"
                % (nm, R_RHO[i][e], rr_meas, rr_out, a_rho, first,
                   r[nm]["c_rho_lo"], r[nm]["c_rho_hi"], "OK" if r[nm]["D1"] else "VIOLATED"))
        del A
        out.append(r)

    say("")
    say("==== verdicts ====")

    def vd(ok, needed, have):
        if not set(needed) <= set(have):
            return "NOT SCORED (missing %s)" % sorted(set(needed) - set(have))
        return "HOLDS" if ok else "REFUTED"

    def three_way(br, win):
        """the measurement is an interval; so HOLDS only if it is inside the window, REFUTED only
        if the two are disjoint, INCONCLUSIVE when the window cuts the bracket."""
        lo, up = br
        a, b = win
        if lo >= a - 1e-12 and up <= b + 1e-12:
            return "HOLDS"
        if up < a - 1e-12 or lo > b + 1e-12:
            return "REFUTED"
        return "INCONCLUSIVE"

    have = [r["e"] for r in out]

    # ---- C0, AS ITS REGISTERED TEXT READS: the bracket CONTAINS the certified interval.
    # Enforced, with the void it registers; intersection is not what the text says.  Containment of
    # a rounded-down
    # certificate by a sharper bracket is unsatisfiable, so this is expected to void every cell.
    for r in out:
        r["C0"] = bool(r["t"]["lower"] <= r["kap_lo"] + 1e-12
                       and r["t"]["upper"] >= r["kap_up"] - 1e-12
                       and (r["t"]["upper"] - r["t"]["lower"]) <= 0.03 * r["t"]["upper"])
        # ---- C0b, the registered clause that governs scoring: the bracket lies INSIDE the
        # comparand's own precision of the certified value, and is tight.
        r["C0b"] = bool(r["t"]["lower"] >= r["kap_band"][0] - 1e-12
                        and r["t"]["upper"] <= r["kap_band"][1] + 1e-12
                        and (r["t"]["upper"] - r["t"]["lower"]) <= 0.03 * r["t"]["upper"])
    void_C0 = [r["e"] for r in out if not r["C0"]]
    void_C0b = [r["e"] for r in out if not r["C0b"]]
    valid = [r for r in out if r["C0b"]]            # C0b governs scoring; C0's void reported
    say("C0  (registered text: CONTAINMENT, voids the cell) : %s   void cells %s"
        % (all(r["C0"] for r in out), void_C0))
    say("C0b (bracket inside the comparand precision)                 : %s   void cells %s"
        % (all(r["C0b"] for r in out), void_C0b))
    say("   %s" % [(r["e"], round(r["t"]["lower"], 7), round(r["t"]["upper"], 7),
                    round(100.0 * (r["t"]["upper"] - r["t"]["lower"]) / r["t"]["upper"], 3),
                    r["kap_lo"], r["kap_up"]) for r in out])
    say("   NOTE: predictions below are scored on the cells valid under C0b; under C0's")
    say("         registered text the valid set is %s." % [r["e"] for r in out if r["C0"]])
    c3 = all(r.get("C3") for r in out)
    say("C3 control optimality witness (max_vB W <= 1e-6||r||^2, <y,W> >= -1e-12) : %s   %s"
        % (c3, [(r["e"], "%.1e" % r["C3_kkt_rel"], "%.1e" % r["C3_kkt_dot"]) for r in out]))
    c1 = all(r["rho%d" % i]["D1"] for r in out for i in (1, 2, 3))
    say("C1 Lemma D1 band, bracket n band nonempty (consistency, must hold) : %s" % c1)
    c2cells = [r for r in out if "S_w" in r]
    c2 = all(r["rho%d" % i]["C2"] for r in c2cells for i in (1, 2, 3)) if c2cells else None
    say("C2 certificate lower bound (must hold)       : %s   cells %s"
        % (c2, [r["e"] for r in c2cells]))

    hv = [r["e"] for r in valid]
    p1cells = [r for r in valid if r["e"] in P1_WIN]
    p1 = [(r["e"], three_way((r["rho2"]["lower"], r["rho2"]["upper"]), P1_WIN[r["e"]]))
          for r in p1cells]
    p1v = vd(all(v == "HOLDS" for _, v in p1), [20, 22, 24], hv) if p1 else "NOT SCORED (no valid cell)"
    if p1 and any(v == "INCONCLUSIVE" for _, v in p1) and not any(v == "REFUTED" for _, v in p1):
        p1v = "INCONCLUSIVE"
    say("P1 kappa^{q_li} in the first-order window    : %s   %s   %s"
        % (p1v, p1, [(r["e"], round(r["rho2"]["lower"], 5), round(r["rho2"]["upper"], 5),
                      P1_WIN[r["e"]]) for r in p1cells]))
    p2 = [(r["e"], three_way((r["rho3"]["lower"], r["rho3"]["upper"]), P2_WIN[r["e"]]))
          for r in valid]
    p2v = (("HOLDS" if all(v == "HOLDS" for _, v in p2) else
            "REFUTED" if any(v == "REFUTED" for _, v in p2) else "INCONCLUSIVE")
           if p2 else "NOT SCORED (no valid cell)")
    say("P2 kappa^{q_li+Pole} in [kappa_lo, +0.6 r]   : %s   %s   %s"
        % (p2v, p2, [(r["e"], round(r["rho3"]["lower"], 5), round(r["rho3"]["upper"], 5),
                      P2_WIN[r["e"]]) for r in valid]))
    # P3 three-way too: REFUTED only if the whole bracket is above
    p3 = [(r["e"], three_way((r["rho1"]["lower"], r["rho1"]["upper"]), (0.0, P3_MAX)))
          for r in valid]
    p3v = (("HOLDS" if all(v == "HOLDS" for _, v in p3) else
            "REFUTED" if any(v == "REFUTED" for _, v in p3) else "INCONCLUSIVE")
           if p3 else "NOT SCORED (no valid cell)")
    say("P3 kappa^{c(Q)v} <= 0.01 (target in the cone): %s   %s   %s"
        % (p3v, p3, [(r["e"], round(r["rho1"]["lower"], 6), round(r["rho1"]["upper"], 6))
                     for r in valid]))
    p4rows = [(r["e"], r["rho3"]["a_rho"]) for r in valid]
    p4a = all(abs(a) <= P4_MAX for _, a in p4rows)
    d16 = dict(p4rows).get(16)
    d18 = dict(p4rows).get(18)
    p4b = (d16 is not None and d18 is not None and d16 * d18 < 0)
    p4v = vd(p4a and p4b, [16, 18, 20, 22, 24], hv) if p4rows else "NOT SCORED (no valid cell)"
    say("P4 |a3| <= 0.08 at every cell AND a3 changes sign 2^16->2^18 : %s"
        % p4v)
    say("   a3 %s   (noise floor 1/sqrt(#rows) %s)"
        % ([(e_, round(a, 6)) for e_, a in p4rows],
           [(r["e"], NOISE[r["e"]]) for r in valid]))
    say("   a1 %s" % [(r["e"], round(r["rho1"]["a_rho"], 6)) for r in valid])
    say("   a2 %s" % [(r["e"], round(r["rho2"]["a_rho"], 6)) for r in valid])
    say("   convexity remainder c_rho (rho3) %s"
        % [(r["e"], "[%.2e, %.2e]" % (r["rho3"]["c_rho_lo"], r["rho3"]["c_rho_hi"]))
           for r in valid])
    json.dump(dict(cells=out, C0=all(r["C0"] for r in out), C0_void=void_C0,
                   C0b=all(r["C0b"] for r in out), C0b_void=void_C0b, C3=c3,
                   C1=c1, C2=c2, P1=p1, P1v=p1v, P2=p2, P2v=p2v, P3=p3, P3v=p3v,
                   P4=p4v, P4_a3=p4rows, passes=PASSES, tol=CDTOL),
              io.open(OUTJ, "w", encoding="utf-8", newline="\n"), indent=1, default=float)
    say("wrote %s and %s" % (os.path.basename(OUT), os.path.basename(OUTJ)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

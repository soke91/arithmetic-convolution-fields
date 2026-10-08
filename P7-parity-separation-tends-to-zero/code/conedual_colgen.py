# -*- coding: utf-8 -*-
r"""conedual_colgen -- kappa(2^e) by column generation WITHOUT the full incidence matrix.

WHY.  The compute node is offered with RAM <= 32 GB, threads <= 16, intermediates
<= 10 GB.  The September pipeline (`gramint_nnls_full` -> `conedual_targeted_repair` ->
`gramint_kappa_certify`) builds the whole divisor-incidence matrix and runs FISTA on every
column.  `results/conedual_build_memory_probe.txt` measured, at 2^30: rows 19,920,
|vB| 86,592,333, nnz 968,732,119, build peak 15.9 GB, float64 COO 14.4 GB + CSC 11.1 GB --
over the cap -- and 60,000 FISTA iterations at that nnz is > 100 h.  So the matrix has to go.

WHAT REPLACES IT.  kappa = dist(t, cone{A_l : l in vB}) / ||p||, A[i,j] = 1 iff rows[i] | vB[j],
t_i = #{m in vR : rows[i] | m} / |vR| (the same t, p and normalisation as `gramint_nnls_full`
and `gramint_kappa_certify`).  Column generation:

  master   a restricted nnls on a working set S of columns, its incidence built for THOSE
           columns only (|S| ~ 10^4-10^5 columns x ~11 nonzeros = a few 10^6 nonzeros, MB
           not GB), solved by COORDINATE DESCENT (`k_cd_nnls`): A's entries are 1,
           so the exact line minimum along column j is y_j <- max(0, y_j + A_j^T r / nnz_j)
           with r = t - A y carried incrementally -- 2 x nnz touches per pass, no linear
           algebra, warm-started from the previous round's y.  Any feasible y gives the
           UPPER bound ||t - A_S y|| / ||p|| >= kappa.  An earlier block Lawson-Hanson
           exchange survives as an optional phase-2 polish (`--ls-budget > 0`, off by
           default); it is what made 2^22 stall, see `k_cd_nnls`.
  pricing   one sweep over the band: W(l) = sum_{d | l, d <= Q} w_d for EVERY l in vB from the
           master's residual w = t - A_S y, by the accumulator idiom of
           `gramint_kappa_certify.sweep` -- but BLOCKED: an (N+1)-length accumulator is 8.6 GB
           at 2^30, so the N-range is swept in 512 KB blocks (L2-resident), one numba kernel,
           `prange` over contiguous groups of blocks, the band gathered per block by
           `searchsorted` on the sorted vB.  Peak for the sweep is then
           nchunk x 512 KB + |vB| (the output), not 8.6 GB.
  add       the most violated columns (largest W), up to --add per round; the working set for
           the next round is the current positive support union those.
  stop      max_vB W <= --price-tol * ||w||^2 (the relative loss of the d=1 repair is
           max_vB W / ||w||^2), or --max-rounds, or --budget hours.
  repair    what is left of the violation is repaired through the cheapest rows, not through
           d = 1: each violating column is pushed down through its divisor row of
           smallest t_d, the per-row decrement being the max requirement over the columns
           assigned to it, so max_vB W <= tau holds by construction and is then re-measured
           exactly by a full sweep.  A ladder of tau is tried and the best bound kept.
  certify   the final row vector is handed to `gramint_kappa_certify.certify` -- the integer
           the integer certificate's clauses K0-K4 are NOT reimplemented here; the npz is written in that
           tool's `rowvec`-repaired format (rows, w, t, tau, bound, src) and the tool's own
           sweep, DEN = 10^12 rounding, J / H / D_p integers and the six-decimal floor c6 are
           what produce the reported lower bound.  `T.build` is shimmed with this run's cached
           cell so the certifier does not rebuild it (335 s and 15.9 GB at 2^30).

MEMORY, per array, at 2^30 (rows 19,920; |vB| 86,592,333; |vR| 89,520,933):
  cell build, FIRST run only ............ peak 15.9 GB (measured; `conedual_the_four_odd_rows.build`)
  cell cache on disk .................... vB 0.69 + vR 0.72 + rows/p/t/rom ~0.001 = 1.41 GB
  colgen phase: vB int64 ................ 0.69 GB
                W float64 ............... 0.69 GB
                sweep block buffers ..... nchunk x BS x 8 B = 64 x 512 KB = 0.03 GB
                t, rows, w, p ........... < 0.001 GB
                A_S (csc, |S| ~ 6e4) .... 6.6e5 nnz x 12 B = 0.008 GB
                sig int8 (t, once) ...... 1.07 GB, freed (and cached to disk)
                                           => colgen peak ~2.2 GB
  certification phase (the certifier's own arrays, not ours):
                acc int64/float64 ....... 8.6 GB (one at a time, freed inside `sweep`)
                WB/WR float64 + int64 ... 1.41 + 1.41 GB
                vB/vR int64 ............. 1.41 GB
                WR.tolist() for exact J . ~3.2 GB (transient, after acc is freed)
                                           => certification peak ~12 GB
  So the binding peak at 2^30 is the FIRST-run cell build at 15.9 GB; every later phase and
  every warm (cached) run stays under 12 GB.

CLAUSES:
  V1  at every validated cell the CERTIFIED lower bound (the integer c6/10^6 from
      `gramint_kappa_certify`) is >= 0.99 x the ledger's best certified value at that cell
      (2^22 0.0760823, 2^24 0.0527414, 2^26 0.0345167), and the upper bound from the
      master's primal is >= that ledger value (consistency).
  V2  the projected 2^30 peak RSS is <= 24 GB.
  V3  the projected 2^30 wall time is <= 36 h with 16 threads.
  Failure of any of these is reported, not patched: a pipeline sent to the compute node unvalidated is
  worth less than a failed validation honestly reported.

THE SEPARATOR NOW CARRIES THE PRIMAL TOO.  `conedual_colgen_sep_2e<e>*.npz`
gains three fields: `y_support` (indices into vB, the working set of the round that achieved the best
upper bound, restricted to y > 0), `y_values` (that y, on the t = c_vR/R scale), and `ub` (that
upper bound).  Nothing else changes -- no clause, no algorithm, no output number.  WHY: the paper's
eq:cert needs an integral dual AND a rational primal, and this script was the only way to reach
2^26/2^28, where it stored the dual only.  The dual alone cannot be used to recover the
primal: the separator is REPAIRED (per-row), so its theta ordering no longer locates the
support -- measured at 2^18, the 301 support columns are spread to rank 6,196 of 18,685, so a third
of the band would have to be solved to find them.  With these fields the certification is the same
storage-and-rounding step that already works at 2^22 and 2^24 (`pa_certify.float_optimum` reads
exactly this pair out of `conedual_exact_dual_2e*.npz`).

NOT CLAIMED.  Nothing about kappa's trend; nothing about cells not run.  The float bounds
printed per round are not certificates -- the certificate is the integer clause block.

    python code/conedual_colgen.py 24 --threads 16
    python code/conedual_colgen.py 30 --threads 16 --budget 30 --warm <npz>
"""
from __future__ import annotations

import os
import sys

# ---- threads BEFORE numpy/numba import (BLAS and OMP read these at load time) -------------
def _early_threads() -> int:
    n = 0
    for i, a in enumerate(sys.argv):
        if a == "--threads" and i + 1 < len(sys.argv):
            n = int(sys.argv[i + 1])
        elif a.startswith("--threads="):
            n = int(a.split("=", 1)[1])
    if n <= 0:
        n = min(16, os.cpu_count() or 1)
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS"):
        os.environ[v] = str(n)
    return n


NTHREAD = _early_threads()

import argparse        # noqa: E402
import gc              # noqa: E402
import importlib.util  # noqa: E402
import io              # noqa: E402
import json            # noqa: E402
import math            # noqa: E402
import time            # noqa: E402

import numpy as np                      # noqa: E402
import psutil                           # noqa: E402
import scipy.sparse as sp               # noqa: E402
from numba import njit, prange, set_num_threads   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
# `regnote.py` sits beside this script.
_regdir = os.path.dirname(os.path.abspath(__file__))
if not os.path.isfile(os.path.join(_regdir, "regnote.py")):
    raise ImportError("regnote.py must sit beside this script")
if _regdir not in sys.path:
    sys.path.insert(0, _regdir)
from regnote import Note, fdown, fup    # noqa: E402


ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
RES = os.path.join(ROOT, "results")

LEDGER_BEST = {22: 0.0760823, 24: 0.0527414, 26: 0.0345167}   # the comparands V1 names
# ... but those are not the ledger's best: reading every ok entry of
# results/gramint_kappa_certify_ledger.json and taking the largest bound_int per cell gives
# 2^22 0.0763312 (conedual_exact_dual_2e22.npz, c6 76331) and 2^24 0.0527629
# (conedual_exact_dual_2e24.npz, c6 52762) -- above the earlier parenthetical numbers; 2^26 0.0345167
# agrees; 2^28 0.0237821 exists too.  V1's WORDING is "the ledger's best certified value at
# that cell", so the registered numbers above are kept (a clause is not re-registered after
# the fact) and the stricter true-best verdict is reported BESIDE it, never instead of it.
LEDGER_TRUE = {22: 0.0763312, 24: 0.0527629, 26: 0.0345167, 28: 0.0237821}
BS = 1 << 16          # sweep block: 65536 float64 = 512 KB, L2-resident; > Q at every cell
                      # we run (Q = 2^(e/2) <= 2^15 at 2^30), so every row has >= 2 multiples
                      # per block and the per-block row loop costs ~4% of the writes.


# =============================== numba kernels =============================================
@njit(parallel=True, cache=True, nogil=True, fastmath=False)
def k_sweep(N, rows, w, base, band, out, buf):
    """out[j] = base + sum_{i >= 1, rows[i] | band[j]} w[i], for sorted `band`.

    Blocked strided accumulation: `buf` is (nchunk, BS) scratch, one row per worker chunk.
    No (N+1)-length array is ever allocated.  `base` is w[0] (the d = 1 row divides every
    column) for the full W, or 0.0 for the d >= 2 sweep the certifier's normalisation wants.
    """
    nchunk = buf.shape[0]
    bs = buf.shape[1]
    nblk = N // bs + 1
    per = (nblk + nchunk - 1) // nchunk
    nr = rows.size
    for c in prange(nchunk):
        acc = buf[c]
        b0 = c * per
        b1 = min((c + 1) * per, nblk)
        for b in range(b0, b1):
            lo = b * bs
            hi = min(lo + bs, N + 1)
            n = hi - lo
            for k in range(n):
                acc[k] = base
            for i in range(1, nr):
                x = w[i]
                if x != 0.0:
                    d = rows[i]
                    m = ((lo + d - 1) // d) * d
                    while m < hi:
                        acc[m - lo] += x
                        m += d
            jlo = np.searchsorted(band, lo)
            jhi = np.searchsorted(band, hi)
            for j in range(jlo, jhi):
                out[j] = acc[band[j] - lo]


@njit(parallel=True, cache=True, nogil=True)
def k_counts(N, rows, sig, buf):
    """buf[c, i] = #{m in (0, N], sig[m] != 0, rows[i] | m} over chunk c; sum axis 0 for t."""
    nchunk = buf.shape[0]
    bs = 1 << 16
    nblk = N // bs + 1
    per = (nblk + nchunk - 1) // nchunk
    nr = rows.size
    for c in prange(nchunk):
        b0 = c * per
        b1 = min((c + 1) * per, nblk)
        for i in range(nr):
            d = rows[i]
            s = 0
            for b in range(b0, b1):
                lo = b * bs
                hi = min(lo + bs, N + 1)
                m = ((lo + d - 1) // d) * d
                if m == 0:
                    m = d
                while m < hi:
                    s += sig[m]
                    m += d
            buf[c, i] = s


@njit(parallel=True, cache=True, nogil=True)
def k_count_div(cols, rows, cnt):
    """cnt[j] = #{i : rows[i] | cols[j]} -- the nnz of the working set's column j."""
    nr = rows.size
    for j in prange(cols.size):
        l = cols[j]
        c = 0
        for i in range(nr):
            if l % rows[i] == 0:
                c += 1
        cnt[j] = c


@njit(parallel=True, cache=True, nogil=True)
def k_fill_div(cols, rows, indptr, indices):
    """CSC row indices of the working set's incidence, column j at indptr[j]:indptr[j+1]."""
    nr = rows.size
    for j in prange(cols.size):
        l = cols[j]
        k = indptr[j]
        for i in range(nr):
            if l % rows[i] == 0:
                indices[k] = i
                k += 1


@njit(parallel=True, cache=True, nogil=True)
def k_pattern(cols, rows, out):
    """out[j] = a 64-bit hash of column j's divisor-row set.

    Many band columns share their whole set of divisors below Q, so a working set of k
    columns has rank well below k; the duplicates inflate the master's passive set (measured
    at 2^22: |P| = 2089 against 1245 rows) and cost least-squares time for nothing.  Keeping
    one column per hash is safe in BOTH directions: it can only SHRINK the restricted cone,
    so the master's primal stays an upper bound, and the lower bound is priced over every
    column of vB regardless.  A hash collision therefore costs a column, never a claim.
    """
    nr = rows.size
    for j in prange(cols.size):
        l = cols[j]
        h = np.uint64(1469598103934665603)
        for i in range(nr):
            if l % rows[i] == 0:
                h = (h ^ np.uint64(i)) * np.uint64(1099511628211)
                h = h ^ (h >> np.uint64(29))
        out[j] = h


@njit(cache=True, nogil=True, fastmath=False)
def k_cd_nnls(indptr, indices, y, r, order, npass, tol):
    """cyclic coordinate descent for  min_{y >= 0} ||A y - t||, A the implicit 0/1 CSC.

    WHY THIS AND NOT THE EXCHANGE.  An earlier master was a block
    Lawson-Hanson whose every step is a dense pivoted QR on the passive set: at 2^22 that is
    80 `gelsy` solves and 85-170 s per round, and the measured round history
    (`conedual_colgen_2e22.txt`) shows ||Ay-t|| stuck at 7.0e-2 after 15 rounds against a
    restricted optimum near 8.7e-3 -- the master, not the pricing, is why 2^22 did not
    converge.  Extrapolated to 2^30 (rows 19,920, support ~19,700) one such solve is ~m^3/3
    flops ~ 2.6e12 = tens of seconds, so a few hundred of them per round is alone past V3.

    Coordinate descent costs no linear algebra at all.  A's entries are 1, so ||A_j||^2 is
    just the column's nonzero count and the exact minimiser along coordinate j is
    y_j <- max(0, y_j + (A_j^T r) / nnz_j), r = t - A y kept incrementally: 2 x nnz memory
    touches per pass over the whole working set, which is ~10^6 at 2^22 and ~3e7 at 2^30 --
    milliseconds and ~0.1 s respectively.  Thousands of passes therefore cost less than one
    dense solve, and each pass is an exact block of exact line minimisations, so the
    objective is monotone non-increasing and every iterate is feasible (y >= 0), i.e. every
    iterate's residual is an upper bound for kappa ||p|| (rule 44).

    Serial on purpose: r is a shared length-`rows` vector (19,920 doubles at 2^30, L2-resident)
    and the Gauss-Seidel dependency through it is what gives the method its rate; the
    parallelism in this script lives in the sweep, which is where the N-sized work is.

    Returns (passes used, max KKT violation max_j {A_j^T r if > 0, |A_j^T r| if y_j > 0}).
    """
    n = order.size
    kkt = 0.0
    p = 0
    for p in range(1, npass + 1):
        kkt = 0.0
        for q in range(n):
            j = order[q]
            a = indptr[j]
            b = indptr[j + 1]
            nj = b - a
            if nj <= 0:
                continue
            g = 0.0
            for k in range(a, b):
                g += r[indices[k]]
            yj = y[j]
            if g > 0.0:
                v = g
            elif yj > 0.0:
                v = -g
            else:
                v = 0.0
            if v > kkt:
                kkt = v
            new = yj + g / nj
            if new < 0.0:
                new = 0.0
            dlt = new - yj
            if dlt != 0.0:
                y[j] = new
                for k in range(a, b):
                    r[indices[k]] -= dlt
        if kkt <= tol:
            break
    return p, kkt


@njit(parallel=True, cache=True, nogil=True)
def k_cheapest(cols, rows, cost, out):
    """out[j] = the index i >= 1 with rows[i] | cols[j] of least cost[i] (-1 if none)."""
    nr = rows.size
    for j in prange(cols.size):
        l = cols[j]
        bi = -1
        bt = 1.0e300
        for i in range(1, nr):
            if l % rows[i] == 0:
                if cost[i] < bt:
                    bt = cost[i]
                    bi = i
        out[j] = bi


# =============================== cell build + disk cache ===================================
def _load(name, fname):
    saved = sys.argv
    sys.argv = [saved[0]]
    try:
        spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    finally:
        sys.argv = saved
    return mod


def cache_dir(e):
    return os.path.join(RES, "conedual_colgen_cache_2e%d" % e)


def cell(e, note, nchunk, force=False):
    """the cell, from the disk cache if it is there, else built once and cached.

    The 2^30 build is 335 s and peaks at 15.9 GB; it runs once per cell and every later run
    memory-maps `rows/vB/vR/p/t/rom` out of `results/conedual_colgen_cache_2e<e>/`.
    """
    cd = cache_dir(e)
    meta = os.path.join(cd, "meta.json")
    if os.path.exists(meta) and not force:
        m = json.load(io.open(meta, encoding="utf-8"))
        S = {k: np.load(os.path.join(cd, k + ".npy"), mmap_mode="r")
             for k in ("rows", "vB", "vR", "p", "t", "rom")}
        S = {k: np.asarray(v) for k, v in S.items()}      # into RAM: 1.41 GB at 2^30
        S.update(N=m["N"], Q=m["Q"], nB=float(m["nB"]), nR=float(m["nR"]))
        note.say("cell 2^%d from cache %s  (rows %d, |vB| %d, |vR| %d)"
                 % (e, os.path.basename(cd), S["rows"].size, S["vB"].size, S["vR"].size))
        return S
    T = _load("tfor", "conedual_the_four_odd_rows.py")
    t0 = time.time()
    S = T.build(e)
    note.say("cell 2^%d built in %.0f s  rows %d  |vB| %d  |vR| %d  Q = %d"
             % (e, time.time() - t0, S["rows"].size, S["vB"].size, S["vR"].size, S["Q"]))
    # t_i = #{m in vR : rows[i] | m} / |vR|, the same vector gramint_nnls_full feeds FISTA
    t1 = time.time()
    N = S["N"]
    sig = np.zeros(N + 1, dtype=np.int8)
    sig[S["vR"]] = 1
    buf = np.zeros((nchunk, S["rows"].size), dtype=np.int64)
    k_counts(np.int64(N), S["rows"].astype(np.int64), sig, buf)
    cnt = buf.sum(axis=0)
    del sig, buf
    gc.collect()
    S["t"] = cnt.astype(np.float64) / float(S["vR"].size)
    note.say("t built in %.0f s (blocked parallel counts; t[0] = %.6f, must be 1)"
             % (time.time() - t1, S["t"][0]))
    os.makedirs(cd, exist_ok=True)
    for k in ("rows", "vB", "vR", "p", "t", "rom"):
        np.save(os.path.join(cd, k + ".npy"), np.asarray(S[k]))
    json.dump(dict(N=int(S["N"]), Q=int(S["Q"]), nB=float(S["nB"]), nR=float(S["nR"]), e=e),
              io.open(meta, "w", encoding="utf-8", newline="\n"), indent=1)
    note.say("cell cached to %s (%.2f GB)"
             % (os.path.basename(cd),
                sum(os.path.getsize(os.path.join(cd, f)) for f in os.listdir(cd)) / 2 ** 30))
    gc.collect()
    return S


# =============================== the restricted master =====================================
def dedup(act, vB, rows, keep):
    """one column per divisor pattern, with `keep` (the current support) always retained."""
    h = np.empty(act.size, dtype=np.uint64)
    k_pattern(np.ascontiguousarray(vB[act]), rows, h)
    pri = np.zeros(act.size, dtype=np.int64)
    if keep.size:
        pri[np.isin(act, keep)] = -1          # the support sorts first inside its pattern
    order = np.lexsort((pri, h))
    hs = h[order]
    first = np.ones(hs.size, dtype=bool)
    first[1:] = hs[1:] != hs[:-1]
    return np.sort(act[order[first]])


def incidence(cols, rows, note=None):
    cnt = np.empty(cols.size, dtype=np.int64)
    k_count_div(cols, rows, cnt)
    indptr = np.zeros(cols.size + 1, dtype=np.int64)
    np.cumsum(cnt, out=indptr[1:])
    nnz = int(indptr[-1])
    indices = np.empty(nnz, dtype=np.int32)
    k_fill_div(cols, rows, indptr, indices)
    A = sp.csc_matrix((np.ones(nnz), indices, indptr), shape=(rows.size, cols.size))
    return A, nnz


def ls_free(A, idx, t, dense_max):
    """unconstrained least squares on the free (passive) columns `idx`.

    <= dense_max columns: a dense QR with column pivoting (`gelsy`), rank-safe and what the
    September dense nnls effectively solved each round.  Above it the normal equations with a
    Cholesky and one step of iterative refinement -- the Gram matrix is formed by a SPARSE
    product (cost sum_i k_i^2 with k_i the row's nnz, ~7 x 10^8 at |P| = 2 x 10^4, seconds)
    and costs |P|^2 x 8 B = 3.1 GB at |P| = 19,700 (the 2^30 support size), while a dense
    19,920 x 19,700 factorisation is not affordable.  The refinement step recovers the
    accuracy the normal equations' squared condition number costs.
    """
    AF = A[:, idx]
    if idx.size <= dense_max:
        from scipy.linalg import lstsq as dlstsq
        x = dlstsq(AF.toarray(), t, lapack_driver="gelsy", check_finite=False)[0]
        return x, "gelsy"
    G = (AF.T @ AF).toarray()
    b = AF.T.dot(t)
    G[np.diag_indices_from(G)] += 1e-11 * max(1.0, float(np.abs(G).max()))
    try:
        from scipy.linalg import cho_factor, cho_solve
        c = cho_factor(G, lower=True, check_finite=False)
        x = cho_solve(c, b, check_finite=False)
        x = x + cho_solve(c, AF.T.dot(t - AF.dot(x)), check_finite=False)
        how = "chol+refine"
    except Exception:
        x = np.linalg.solve(G, b)
        how = "normal-eq"
    del G
    return x, how


def fista(A, t, y, iters):
    """accelerated projected gradient, warm-startable; it costs two sparse matvecs per step
    and is used only to put the block Lawson-Hanson on a feasible point with a plausible
    support when there is no warm y (round 1)."""
    At = A.T
    rng = np.random.RandomState(0)
    x = rng.standard_normal(A.shape[1])
    x /= np.linalg.norm(x)
    L = 1.0
    for _ in range(30):
        v = At.dot(A.dot(x))
        L = float(np.linalg.norm(v))
        if L == 0.0:
            return y
        x = v / L
    L *= 1.02
    z = y.copy()
    tk = 1.0
    for _ in range(iters):
        g = At.dot(A.dot(z) - t)
        yn = np.maximum(z - g / L, 0.0)
        if float(g.dot(yn - y)) > 0.0:
            tk = 1.0
        tn = 0.5 * (1.0 + math.sqrt(1.0 + 4.0 * tk * tk))
        z = yn + ((tk - 1.0) / tn) * (yn - y)
        y, tk = yn, tn
    return y


def cd_nnls(A, t, y0, passes, cd_tol, chunk=500):
    """`k_cd_nnls` driven in chunks, with the residual recomputed between them.

    Two reasons for the chunking.  (a) Drift: r is updated incrementally inside the kernel,
    so after ~10^5 passes its error is no longer negligible against the 10^-12 KKT tolerance;
    one sparse matvec per chunk resets it exactly.  (b) The tolerance is RELATIVE to the
    current ||r||^2, which falls by two orders of magnitude over a master call -- the target
    is recomputed per chunk so the master always ends ~1/cd_tol tighter than the pricing
    threshold `price_tol * ||w||^2` that decides the round, rather than tighter than the
    threshold it had on entry.
    """
    n = A.shape[1]
    y = np.maximum(np.asarray(y0, dtype=np.float64), 0.0)
    indptr = np.ascontiguousarray(A.indptr, dtype=np.int64)
    indices = np.ascontiguousarray(A.indices, dtype=np.int64)
    order = np.random.RandomState(12345).permutation(n).astype(np.int64)
    r = t - A.dot(y)
    tot, kkt = 0, float("inf")
    while tot < passes:
        k = min(chunk, passes - tot)
        tol = cd_tol * max(float(np.dot(r, r)), 1e-300)
        p, kkt = k_cd_nnls(indptr, indices, y, r, order, k, tol)
        tot += p
        r = t - A.dot(y)                 # exact residual: kills the incremental drift
        if p < k:                        # the kernel met the tolerance inside the chunk
            break
    return y, float(np.linalg.norm(r)), tot, kkt


def master(A, t, y0, note, tag="", maxouter=400, maxinner=60, addmax=0, dense_max=6000,
           fista_iters=0, ls_budget=0, cd_passes=200000, cd_tol=1e-7):
    """restricted nnls  min_{y >= 0} ||A y - t||.

    Phase 1: coordinate descent, `k_cd_nnls` -- see its docstring for why the
    exchange had to go.  Phase 2, only if `--ls-budget > 0`: the block Lawson-Hanson below,
    seeded from the coordinate-descent point, as a polish; the better of the two feasible
    iterates is returned, so the polish can never worsen the bound.

    The old phase-2 docstring follows, since the code is unchanged -- BLOCK Lawson-Hanson.

    Passive set P (y > 0), the rest at zero.  Outer step: price the working set,
    W = A^T(t - A y), and move up to `kcap` columns of largest W > 0 into P -- never past
    `rows` - 1 columns, since a free set beyond the row rank makes the least squares
    non-unique and the exchange cycles (measured at 2^18: block principal pivoting sat at
    |F| = 2163 > 314 rows and never left ||Ay-t|| = 0.131).  Inner step, exact least squares
    on P, then whichever of the two infeasibility cases applies:

      * s_j <= 0 at a column with y_j = 0 -- a block-added column that does not in fact
        enter: dropped outright (the Lawson-Hanson ratio step gives alpha = 0 there and,
        applied to a whole block, empties P -- measured, |P| collapsing to 0 at 2^18);
      * s_j <= 0 at a column with y_j > 0 -- the genuine Lawson-Hanson case: the ratio step
        y <- y + alpha(s - y), alpha = min y_j/(y_j - s_j) > 0, which keeps y >= 0 and drops
        at least one column with the objective strictly down.

    `kcap` halves on a stalled outer round and doubles on a productive one, down to 1 --
    at kcap = 1 this IS textbook Lawson-Hanson, whose single added column is guaranteed to
    enter, so the method cannot stall; blocks are what make a 20,000-column support
    reachable at 2^30.  Every iterate is feasible, so its residual is an upper bound for
    kappa ||p|| (rule 44); the best is tracked and returned.
    """
    m, n = A.shape
    At = A.T
    if addmax <= 0:
        addmax = max(2000, m)
    gtol = 1e-13 * max(1.0, float(np.abs(At.dot(t)).max()))

    # ---------------- phase 1: coordinate descent -----------------------------------------
    tcd = time.time()
    y_cd, r_cd, npass, kkt = cd_nnls(A, t, y0, cd_passes, cd_tol)
    tcd = time.time() - tcd
    if ls_budget <= 0:
        if note is not None:
            note.say("      master%s: cd %.0f s, %d passes, kkt %.2e -> %.12e | |P| %d"
                     % (tag, tcd, npass, kkt, r_cd, int((y_cd > 0).sum())))
        return y_cd, r_cd

    # ---------------- phase 2: Lawson-Hanson polish, seeded from phase 1 -------------------
    y0 = y_cd
    y = np.maximum(np.asarray(y0, dtype=np.float64), 0.0)
    tf = time.time()
    if fista_iters > 0:
        # FISTA first, warm-started from the previous round's y: two sparse matvecs over
        # ~10^6 nonzeros per step (sub-millisecond), so 10^4-10^5 steps are minutes, and it
        # gets the ITERATE and its support close cheaply.  The exchange phase below then
        # needs a handful of (expensive) least-squares solves instead of hundreds, which is
        # what decides whether 2^30 fits in 36 h -- measured at 2^22, pure exchange spent
        # 80 solves and 85-170 s per round and was still at kappa <= 0.39 after five.
        y = fista(A, t, y, fista_iters)
    tf = time.time() - tf
    P = y > 0.0
    if int(P.sum()) > 3 * m:                   # trim a runaway warm support
        idx = np.nonzero(P)[0]
        keep = idx[np.argsort(-y[idx])[:3 * m]]
        P = np.zeros(n, dtype=bool)
        P[keep] = True
        y = np.where(P, y, 0.0)
    best_y, best_r = y.copy(), float(np.linalg.norm(A.dot(y) - t))
    r_fista = best_r
    nls, how, outer, kcap, safe, nstall = 0, "-", 0, addmax, False, 0
    tl = time.time()
    for outer in range(1, maxouter + 1):
        r_before = float(np.linalg.norm(A.dot(y) - t))
        for _ in range(maxinner):
            idx = np.nonzero(P)[0]
            if idx.size == 0 or nls >= ls_budget:
                break
            s, how = ls_free(A, idx, t, dense_max)
            nls += 1
            if s.min() > 0.0:
                y = np.zeros(n)
                y[idx] = s
                P = y > 0.0
                break
            yi = y[idx]
            drop = (s <= 0.0) & (yi <= 0.0)
            if not safe:
                # block drop: EVERY non-positive column leaves at once.  P strictly shrinks,
                # so the inner loop terminates; the objective may rise, which is why the best
                # feasible iterate is kept and why two stalled outer rounds switch `safe` on.
                # One-at-a-time ratio steps cost 400 least-squares solves per master round at
                # 2^18 (measured) -- at 2^30, where one solve is ~20 s, that is 90 h.
                drop = s <= 0.0
            if drop.any():
                P[idx[drop]] = False
                y = np.where(P, y, 0.0)
                continue
            den = yi - s
            bad = (s <= 0.0) & (den > 0.0)
            if not bad.any():
                y = np.zeros(n)
                y[idx] = np.maximum(s, 0.0)
                P = y > 0.0
                break
            alpha = float(np.min(yi[bad] / den[bad]))
            yn = np.zeros(n)
            yn[idx] = yi + alpha * (s - yi)
            yn[yn <= 1e-15 * max(1.0, float(yn.max()))] = 0.0
            y, P = yn, yn > 0.0
        r = float(np.linalg.norm(A.dot(y) - t))
        if r < best_r:
            best_y, best_r = y.copy(), r
        if r >= r_before - 1e-15 * max(1.0, r_before):
            kcap = max(1, kcap // 4)
            nstall += 1
            if nstall >= 2 and not safe:
                safe, nstall = True, 0
                y, P = best_y.copy(), best_y > 0.0
        else:
            kcap, nstall = min(addmax, 2 * kcap), 0
        W = At.dot(t - A.dot(y))
        # no cap on |P| from `m`: at the exact least squares the residual is orthogonal to
        # span(A_P), so W_j != 0 proves A_j is OUTSIDE that span and adding it is safe; and
        # a block that turns out linearly dependent is self-cleaning, since the pivoted QR
        # returns s_j = 0 on the dependent columns, which the y_j = 0 drop rule removes.
        # (The band has many columns with identical divisor sets below Q, so the working
        # set's rank is well below its size -- capping |P| at `rows` stopped the method
        # dead at 2^18 after one exchange, ||Ay-t|| = 0.458 against the optimum 0.0231.)
        cand = np.nonzero((~P) & (W > gtol))[0]
        if cand.size == 0 or nls >= ls_budget:
            break
        k = min(cand.size, kcap, max(1, 3 * m - int(P.sum())))
        pick = cand[np.argpartition(-W[cand], k - 1)[:k]] if k < cand.size else cand
        P[pick] = True
        y = np.where(P, y, 0.0)
    if best_r > r_cd:                    # the polish did not beat phase 1: keep phase 1
        best_y, best_r = y_cd, r_cd
    if note is not None:
        note.say("      master%s: cd %.0f s / %d passes / kkt %.2e -> %.9e | fista %.0f s ->"
                 " %.9e | exchange %.0f s, %d outer / %d ls (%s) | best %.12e | |P| %d"
                 % (tag, tcd, npass, kkt, r_cd, tf, r_fista, time.time() - tl, outer, nls,
                    how, best_r, int((best_y > 0).sum())))
    return best_y, best_r


# =============================== sweep / pricing glue ======================================
class Sweeper:
    def __init__(self, N, rows, band, nchunk):
        self.N = np.int64(N)
        self.rows = np.ascontiguousarray(rows, dtype=np.int64)
        self.band = np.ascontiguousarray(band, dtype=np.int64)
        self.buf = np.zeros((nchunk, BS), dtype=np.float64)
        self.out = np.empty(self.band.size, dtype=np.float64)

    def __call__(self, w, base, copy=False):
        """returns the shared output buffer -- `copy=True` when the caller must keep it
        across another sweep (the buffer is |vB| float64, 0.69 GB at 2^30, so it is reused)."""
        k_sweep(self.N, self.rows, np.ascontiguousarray(w, dtype=np.float64),
                np.float64(base), self.band, self.out, self.buf)
        return self.out.copy() if copy else self.out


def cert_bound(W2max, w, t, npv):
    """the certifier's own float bound for a row vector: the d = 1 repair w_1 = -max_vB W,
    R = sum_{d>1} w_d t_d, bound = (R - M) / (||w|| ||p||).  Identical algebra to
    gramint_kappa_certify._certify's K0 block."""
    M = float(W2max)
    R = float(np.dot(t[1:], w[1:]))
    nw = math.sqrt(M * M + float(np.dot(w[1:], w[1:])))
    return (R - M) / (nw * npv), M, R, nw


# =============================== targeted repair ==========================================
def repair(w, sw, rows, t, npv, taus, note, maxviol=200000):
    """lower max_vB W through the cheapest rows, not through d = 1.

    Each column above tau is pushed down through its divisor row of least t_d; the row's
    decrement is the MAX requirement over the columns assigned to it, so every assigned
    column lands at or below tau and no other column rises.  max_vB W is then re-measured by
    a full sweep (the saved bound must be the certifier's to 1e-9), and the best tau kept.
    """
    w = w.copy()
    W2 = sw(w, 0.0, copy=True)
    b0, M0, R0, _ = cert_bound(W2.max(), w, t, npv)
    best = (b0, w.copy(), float("nan"), M0, R0)
    eps = M0 + w[0]          # the master's own violation max_vB W_full, = 0 at convergence
    note.say("   repair: max_vB W(d>=2) %.9f   w_1 %.9f   residual violation %.3e"
             % (M0, w[0], eps))
    if eps <= 0.0:
        note.say("   repair: nothing to repair (the master is polar-feasible)")
        return best
    for f in taus:
        tau = M0 - f * eps
        idx = np.nonzero(W2 > tau)[0]
        if idx.size == 0:
            continue
        if idx.size > maxviol:
            k = np.argpartition(W2[idx], idx.size - maxviol)[idx.size - maxviol:]
            idx = idx[k]
            tau = float(W2[idx].min())
        need = W2[idx] - tau
        bi = np.empty(idx.size, dtype=np.int64)
        k_cheapest(sw.band[idx], sw.rows, np.ascontiguousarray(t, dtype=np.float64), bi)
        ok = bi >= 0
        delta = np.zeros(rows.size)
        np.maximum.at(delta, bi[ok], need[ok] * (1.0 + 1e-12))
        wv = w.copy()
        wv[1:] -= delta[1:]
        W2v = sw(wv, 0.0)
        b, M, R, _ = cert_bound(W2v.max(), wv, t, npv)
        note.say("   repair tau %.9f (f %.2f): %d columns, %d rows lowered, cost on R %.3e"
                 "  -> max_vB W(d>=2) %.9f  bound %.7f"
                 % (tau, f, idx.size, int((delta > 0).sum()),
                    float(np.dot(delta[1:], t[1:])), M, fdown(b)))
        if b > best[0]:
            best = (b, wv.copy(), tau, M, R)
    note.say("   repair: best float bound %.7f (against %.7f before repair, x%.4f)"
             % (fdown(best[0]), fdown(b0), best[0] / b0 if b0 > 0 else float("nan")))
    return best


# =============================== certification =============================================
class _Shim:
    """`gramint_kappa_certify._certify` calls `T.build(e)`; at 2^30 that is 335 s and 15.9 GB
    for a cell this run already holds.  Hand it the cell we have."""

    def __init__(self, S):
        self._S = S

    def build(self, e_):
        return self._S


def certify(npz_path, S, note):
    C = _load("certify", "gramint_kappa_certify.py")
    del C._lines[:]
    r = C.certify(npz_path, _Shim(S), None)
    for line in r["block"]:
        note.say("   | " + line)
    return r


# =============================== the run ===================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("e", type=int)
    ap.add_argument("--warm", default=None)
    ap.add_argument("--threads", type=int, default=NTHREAD)
    ap.add_argument("--add", type=int, default=0,
                    help="columns added per round (0 = max(20000, 8 x rows))")
    ap.add_argument("--ls-budget", type=int, default=0,
                    help="least-squares solves in the phase-2 polish (0 = coordinate descent only)")
    ap.add_argument("--cd-passes", type=int, default=200000,
                    help="max coordinate-descent passes over the working set per master call")
    ap.add_argument("--cd-tol", type=float, default=1e-7,
                    help="master stops at max_j KKT violation <= cd-tol x ||r||^2")
    ap.add_argument("--dedup", type=int, default=1,
                    help="1 = keep one working-set column per divisor pattern")
    ap.add_argument("--max-rounds", type=int, default=60)
    ap.add_argument("--bpp", type=int, default=400,
                    help="outer Lawson-Hanson rounds per master")
    ap.add_argument("--fista", type=int, default=0,
                    help="FISTA iterations inside the phase-2 polish (0 = off)")
    ap.add_argument("--dense-max", type=int, default=6000,
                    help="free-set size up to which the master's least squares is a dense SVD")
    ap.add_argument("--prune", type=float, default=0.0,
                    help="drop support columns with y <= prune x max(y) from the next working"
                         " set (0 = keep every column with y > 0, the earlier behaviour)")
    ap.add_argument("--price-tol", type=float, default=1e-4,
                    help="stop when max_vB W <= tol * ||w||^2 (the d=1 repair's relative loss)")
    ap.add_argument("--budget", type=float, default=1e9, help="wall-clock budget, hours")
    ap.add_argument("--rebuild", action="store_true")
    ap.add_argument("--no-certify", action="store_true")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    e = a.e
    set_num_threads(a.threads)
    nchunk = 4 * a.threads
    proc = psutil.Process()
    T0 = time.time()

    tag = a.tag
    note = Note(os.path.join(RES, "conedual_colgen_2e%d%s.txt" % (e, tag)),
                "conedual_colgen -- kappa(2^%d) by column generation, no full incidence" % e)
    note.say("threads %d (OMP/MKL/numba), sweep block %d doubles (%.0f KB), %d chunks"
             % (a.threads, BS, BS * 8 / 1024.0, nchunk))
    note.say("V1 certified lower >= 0.99 x ledger (%s); V2 projected 2^30 RSS <= 24 GB;"
             " V3 projected 2^30 wall <= 36 h"
             % (("%.7f" % LEDGER_BEST[e]) if e in LEDGER_BEST else "no ledger entry at this cell"))
    note.say()

    S = cell(e, note, nchunk, force=a.rebuild)
    rows = np.ascontiguousarray(S["rows"], dtype=np.int64)
    vB = np.ascontiguousarray(S["vB"], dtype=np.int64)
    t = np.ascontiguousarray(S["t"], dtype=np.float64)
    N, Q = int(S["N"]), int(S["Q"])
    npv = float(np.linalg.norm(S["p"]))
    nrow, nB = rows.size, vB.size
    ADD = a.add if a.add > 0 else max(20000, 8 * nrow)
    note.say("rows %d  |vB| %d  ||p|| %.8e  add/round %d  price-tol %.1e  budget %.1f h"
             % (nrow, nB, npv, ADD, a.price_tol, a.budget))

    sw = Sweeper(N, rows, vB, nchunk)
    note.say("sweep warm-up (numba compile) ...")
    tw = time.time()
    sw(t, t[0])
    note.say("   first sweep (compile included) %.1f s" % (time.time() - tw))
    tw = time.time()
    sw(t, t[0])
    t_sweep = time.time() - tw
    note.say("   second sweep %.2f s  (%.2f Gupd/s over sum_{d<=Q, mu!=0} N/d updates)"
             % (t_sweep, (float(np.sum(N / rows[1:].astype(np.float64))) + N) / t_sweep / 1e9))
    # CONTROL, small cells only: mean_vR W(d>=2) must equal <t[1:], w[1:]> to float rounding.
    if nB <= 6_000_000:
        swR = Sweeper(N, rows, np.ascontiguousarray(S["vR"], dtype=np.int64), nchunk)
        wq = t.copy()
        wq[0] = 0.0
        mR = float(swR(wq, 0.0).mean())
        dq = float(np.dot(t[1:], wq[1:]))
        note.say("   CONTROL mean_vR W(d>=2) = %.12e against <t,w> = %.12e (dev %.2e)"
                 % (mR, dq, abs(mR - dq)))
        del swR
        gc.collect()

    # --------------------------------------------------- warm start
    w = t.copy()
    act = np.empty(0, dtype=np.int64)
    y = np.empty(0)
    if a.warm:
        z = np.load(a.warm)
        zr = z["rows"].astype(np.int64)
        zw = z["w"].astype(np.float64)
        j = np.searchsorted(rows, zr)
        okm = (j < nrow) & (rows[np.minimum(j, nrow - 1)] == zr)
        w = np.zeros(nrow)
        w[j[okm]] = zw[okm]
        if w[0] == 0.0:
            w[0] = -1.0
        w /= abs(w[0])
        note.say("warm start from %s: %d of %d saved rows mapped onto this cell's %d rows"
                 % (os.path.basename(a.warm), int(okm.sum()), zr.size, nrow))
        if "y_support" in z.files and zr.size == nrow and np.array_equal(zr, rows):
            sup = z["y_support"].astype(np.int64)
            sup = sup[(sup >= 0) & (sup < nB)]
            act = np.unique(sup)
            y = np.zeros(act.size)
            if "y_values" in z.files and z["y_values"].size == sup.size:
                vals = np.zeros(nB)
                vals[sup] = z["y_values"].astype(np.float64)
                y = vals[act]
                del vals
            note.say("   same cell: %d support columns seeded as the first working set" % act.size)

    if act.size == 0:
        W = sw(w, w[0])
        k = min(ADD, nB)
        sel = np.argpartition(W, nB - k)[nB - k:]
        act = np.unique(sel[W[sel] > -np.inf])
        y = np.zeros(act.size)
        note.say("seed working set: %d columns of largest W under the %s vector (max W %.6e)"
                 % (act.size, "warm" if a.warm else "y = 0 (w = t)", float(W.max())))

    # --------------------------------------------------- rounds
    hist = []
    best_lb = (-1.0, None, None)
    ub_best = float("inf")
    # the primal that achieves ub_best, kept so the separator can carry it.  Without it the
    # right half of the paper's eq:cert cannot be evaluated at a cell this script is the only way to
    # reach: the dual alone is not enough, because the separator is REPAIRED and its theta
    # ordering therefore no longer locates the support (measured at 2^18).
    best_pr = (None, None)
    stop = "max-rounds"
    act_prev = act.copy()
    y_prev = y.copy() if y.size == act.size else np.zeros(act.size)
    sup_prev = np.empty(0, dtype=np.int64)
    for rd in range(1, a.max_rounds + 1):
        tr = time.time()
        ta = time.time()
        if a.dedup:
            n_raw = act.size
            act = dedup(act, vB, rows, sup_prev)
        else:
            n_raw = act.size
        cols = np.ascontiguousarray(vB[act])
        A, nnz = incidence(cols, rows)
        t_inc = time.time() - ta
        # warm start: carry the previous round's y onto the columns both sets share
        y0 = np.zeros(act.size)
        if y_prev.size:
            jj = np.searchsorted(act, act_prev)
            m2 = (jj < act.size) & (act[np.minimum(jj, act.size - 1)] == act_prev)
            y0[jj[m2]] = y_prev[m2]
        y, rnorm = master(A, t, y0, note, tag=" r%d" % rd, maxouter=a.bpp,
                          addmax=max(2000, nrow), dense_max=a.dense_max,
                          fista_iters=a.fista, ls_budget=a.ls_budget,
                          cd_passes=a.cd_passes, cd_tol=a.cd_tol)
        w = t - A.dot(y)
        ub = rnorm / npv
        if ub < ub_best:
            # `act` indexes vB and `y` is aligned with it, so this is the same (y_support, y_values)
            # pair `conedual_exact_dual_2e*.npz` stores and `pa_certify.float_optimum` reads.
            sel_pr = y > 0.0
            best_pr = (np.ascontiguousarray(act[sel_pr]), np.ascontiguousarray(y[sel_pr]))
        ub_best = min(ub_best, ub)
        ts = time.time()
        W = sw(w, w[0])
        t_sw = time.time() - ts
        mw = float(W.max())
        nw2 = float(np.dot(w, w))
        # max_vB W(d >= 2) = max_vB W_full - w[0]: the d = 1 row divides every column, so the
        # certifier's d >= 2 sweep is this one shifted -- no second sweep is needed.
        lb_naive, M2, R2, _ = cert_bound(mw - w[0], w, t, npv)
        nviol = int((W > 0.0).sum())
        note.say("round %-3d active %-7d (of %-7d) nnz %-9d  kappa <= %.7f   max_vB W %.3e"
                 "  (loss %.3f%%)  support %d  violators %d   inc %.1f s  sweep %.1f s"
                 "  round %.0f s  RSS %.2f GB"
                 % (rd, act.size, n_raw, nnz, fup(ub), mw, 100.0 * max(mw, 0.0) / nw2,
                    int((y > 0).sum()), nviol, t_inc, t_sw, time.time() - tr,
                    proc.memory_info().rss / 2 ** 30))
        note.say("            d=1-repaired float lower %.7f  (M %.6f  R %.6f)"
                 % (fdown(lb_naive), M2, R2))
        hist.append(dict(round=rd, active=int(act.size), nnz=int(nnz), ub=ub, maxW=mw,
                         support=int((y > 0).sum()), violators=nviol, lb_float=lb_naive,
                         t_inc=t_inc, t_sweep=t_sw, t_round=time.time() - tr,
                         rss_gb=proc.memory_info().rss / 2 ** 30))
        if lb_naive > best_lb[0]:
            best_lb = (lb_naive, w.copy(), rd)
        act_prev, y_prev = act.copy(), y.copy()
        if mw <= a.price_tol * nw2:
            stop = "converged (max_vB W <= %.1e x ||w||^2)" % a.price_tol
            break
        if (time.time() - T0) / 3600.0 > a.budget:
            stop = "budget"
            break
        # support pruning.  The next working set is `support U top-ADD violators`,
        # so every column that ever carried weight is carried forever and the support only
        # grows: measured 26,451 at 2^22 and 31,513 at 2^26 against rows 1,245 and 4,982, and
        # since the round's cost is 2 x nnz x CD passes that growth is what sets the 2^30
        # projection.  Columns whose y is a negligible fraction of the largest are dropped.
        # Safe in both directions: dropping columns can only SHRINK the restricted cone, so
        # the master's primal stays an upper bound (rule 44), and the lower bound is priced
        # over all of vB regardless -- a wrongly dropped column simply reappears as a
        # violator in the next pricing sweep.
        ycut = a.prune * float(y.max()) if (a.prune > 0.0 and y.size and y.max() > 0.0) else 0.0
        sup = act[y > ycut]
        sup_prev = sup.copy()
        cand = np.nonzero(W > 0.0)[0]
        if cand.size > ADD:
            kk = np.argpartition(W[cand], cand.size - ADD)[cand.size - ADD:]
            cand = cand[kk]
        act = np.union1d(sup, cand)
        A = None
        gc.collect()
    note.say("stopped: %s after %d rounds, %.0f s" % (stop, len(hist), time.time() - T0))
    note.say()

    # --------------------------------------------------- repair + save
    wbest = best_lb[1] if best_lb[1] is not None else w
    wbest = wbest / abs(wbest[0]) if wbest[0] != 0.0 else wbest
    b, wr, tau, M, R = repair(wbest, sw, rows, t, npv,
                              (0.5, 0.9, 1.0, 1.2, 1.5), note)
    npz = os.path.join(RES, "conedual_colgen_sep_2e%d%s.npz" % (e, tag))
    ysup, yval = best_pr
    if ysup is None:
        ysup, yval = np.empty(0, dtype=np.int64), np.empty(0, dtype=np.float64)
    np.savez(npz, rows=rows, w=wr, t=t, tau=np.float64(tau if tau == tau else 0.0),
             bound=np.float64(b), src="conedual_colgen round %d" % (best_lb[2] or 0),
             y_support=ysup.astype(np.int64), y_values=yval.astype(np.float64),
             ub=np.float64(ub_best))
    note.say("separator saved: %s   float lower %.7f   upper %.7f   primal %d columns"
             % (os.path.basename(npz), fdown(b), fup(ub_best), ysup.size))
    peak_colgen = proc.memory_info().peak_wset / 2 ** 30
    note.say("peak_wset through the colgen phase: %.2f GB" % peak_colgen)

    # --------------------------------------------------- the integer certificate
    cert = None
    if not a.no_certify:
        note.say()
        note.say("integer certificate (gramint_kappa_certify, clauses K0-K4;"
                 " not reimplemented here):")
        sw.buf = sw.out = None
        A = None
        gc.collect()
        Sfull = dict(S)
        Sfull["vB"] = vB
        Sfull["vR"] = np.ascontiguousarray(S["vR"], dtype=np.int64)
        Sfull["rows"] = rows
        Sfull["N"], Sfull["Q"] = N, Q
        cert = certify(npz, Sfull, note)
        gc.collect()
    peak = proc.memory_info().peak_wset / 2 ** 30
    note.say()
    note.say("peak_wset for the whole run: %.2f GB" % peak)

    # --------------------------------------------------- the clause block
    lb_cert = (cert["c6"] / 1e6) if (cert and cert.get("ok")) else float("nan")
    ref = LEDGER_BEST.get(e)
    if ref is not None:
        v1a = lb_cert >= 0.99 * ref
        v1b = fup(ub_best) >= ref
        note.say("V1 at 2^%d: certified lower %s/10^6 = %.7f against 0.99 x %.7f = %.7f -> %s;"
                 "  upper %.7f >= %.7f -> %s"
                 % (e, ("%d" % cert["c6"]) if cert and cert.get("ok") else "none",
                    lb_cert, ref, 0.99 * ref, "HOLDS" if v1a else "FAILS",
                    fup(ub_best), ref, "HOLDS" if v1b else "FAILS"))
        if not (v1a and v1b):
            note.fail("V1 at 2^%d" % e)
        tru = LEDGER_TRUE.get(e)
        if tru is not None and tru != ref:
            note.say("V1' (same clause, against the ledger's ACTUAL best %.7f at this cell):"
                     " %.7f >= %.7f -> %s;  upper %.7f >= %.7f -> %s"
                     % (tru, lb_cert, 0.99 * tru, "HOLDS" if lb_cert >= 0.99 * tru else "FAILS",
                        fup(ub_best), tru, "HOLDS" if fup(ub_best) >= tru else "FAILS"))
            if not (lb_cert >= 0.99 * tru and fup(ub_best) >= tru):
                note.fail("V1' at 2^%d (stricter true-ledger comparand)" % e)
    else:
        note.say("V1 at 2^%d: no ledger comparand; certified lower %s"
                 % (e, ("%d/10^6" % cert["c6"]) if cert and cert.get("ok") else "none"))

    js = dict(e=e, N=N, Q=Q, rows=int(nrow), nB=int(nB), nR=int(S["vR"].size),
              npv=npv, threads=a.threads, block=BS, chunks=nchunk, add=ADD,
              bpp=a.bpp, ls_budget=a.ls_budget, dedup=a.dedup,
              cd_passes=a.cd_passes, cd_tol=a.cd_tol, prune=a.prune,
              price_tol=a.price_tol, stop=stop, rounds=len(hist),
              sweep_s=t_sweep, wall_s=time.time() - T0,
              peak_wset_gb=peak, peak_colgen_gb=peak_colgen,
              ub=fup(ub_best), lb_float=b, lb_certified=lb_cert,
              c6=(cert["c6"] if cert and cert.get("ok") else None),
              cert_ok=bool(cert and cert.get("ok")), ledger=ref,
              ledger_true=LEDGER_TRUE.get(e), history=hist,
              J=(str(cert["J"]) if cert and cert.get("ok") else None),
              H=(str(cert["H"]) if cert and cert.get("ok") else None),
              D_p=(str(cert["D_p"]) if cert and cert.get("ok") else None),
              npz=os.path.basename(npz))
    json.dump(js, io.open(os.path.join(RES, "conedual_colgen_2e%d%s.json" % (e, tag)),
                          "w", encoding="utf-8", newline="\n"), indent=1)
    note.say("wrote %s" % os.path.basename(os.path.join(RES, "conedual_colgen_2e%d%s.json" % (e, tag))))
    return note.finish()


if __name__ == "__main__":
    sys.exit(main())

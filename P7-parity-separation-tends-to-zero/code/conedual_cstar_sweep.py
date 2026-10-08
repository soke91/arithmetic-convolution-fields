# -*- coding: utf-8 -*-
r"""conedual_cstar_sweep.py -- C*(Q, D) from D = Q/2 down to polylog D, and the

Supporting computation; not used in the paper's proofs or tables.
small-prime-part (Y-smooth) class scheme.

The boundedness of

    C*(Q,D) = min { C : exists a >= 0 on the classes,  a_k <= C nu_B(k),  sum_{k : i|k} a_k = t_i  (i in H) }

is the whole game: the existence theorem needs only the EXISTENCE of a capped weight, not a formula
for one.  Earlier sweeps
measured C* at D = Q/8 and Q/4 only, where it is 3-4 and 12-16.  PART 1 here sweeps D over Q/2^m
(m = 1..7) and over polylog values (L^2, L^3, L^4 with L = log Q) at 2^16, 2^18, 2^20, and asks
whether C* -> 1 as D/Q -> 0 (C* >= 1 always: sum_k a_k = 1 = sum_k nu_B(k)).

PART 2 is the other class scheme.  Take Y = L^6, the small-prime part
a_Y(n) = prod_{p|n, p<=Y} p, and the explicit even-band weight

    x_l = R_{a_Y(l)} / (R * B_{a_Y(l)})   if a_Y(l) <= Q^{1/2},   else 0,

which is the class measure a = nu_R restricted to the classes a <= Q^{1/2}.  It matches t exactly on
every Y-SMOOTH row d <= Q (not only d <= D: all primes of a Y-smooth d are <= Y, so d | l iff
d | a_Y(l)), its cap is (|vB|/|vR|) max_a R_a/B_a, and (PC) -- Hypothesis 5 of the paper -- is exactly R_a <= C B_a.  PART 2
computes, per Y: that explicit cap, the optimal cap for the same row system, the discarded mass
delta = sum_{a > Q^{1/2}} nu_R(a), the dead classes (B_a = 0 < R_a, which make the explicit weight undefined), and
the kappa bound each construction actually gives, measured as ||t - Ax||_2 / ||p||_2 over ALL rows
d <= Q (a genuine upper bound for every x >= 0, since Ax lies in K_Q).

Note the two aggregations are independent: the class-level sparse matvec (rows in H only) and the
band-level strided sum over all rows d <= Q of the spread weight x_l = a_k/B(k).  C4 cross-checks them.

Compute rules (the project's compute rules): cells via cell_cache.load; no Python loop over rows, band elements
or classes (numpy / scipy.sparse; the per-prime strided loops are the canonical sieve's own, the
subset loops run over <= 2^MMAX masks, and the band row sums use a numba prange kernel); threads
capped; results per cell.

CLAUSES

  C1  CONTROL, reproduction of earlier sweeps.  The sweep's C* at D = Q/8 and D = Q/4 reproduces
      Q/8: 4.0326, 4.3890, 3.1114 and Q/4: 12.1837, 15.5303, 13.8452 (2^16, 2^18, 2^20) to 1e-4
      RELATIVE (4-decimal literals, half-ulp 5e-5 < 1e-4*3.11 = 3.1e-4).  A failure VOIDS that cell.
  C2  CONTROL, LP duality.  Where the dual is computed, the polar value
      max_y <y,t> / sum_k nu_B(k) (sum_{i|k} y_i)_+ equals C* to 1e-6 RELATIVE.  A failure voids the
      (cell, D).
  C3  CONTROL, the assembled-weight audit.  Every weight assembled here satisfies a_k = 0 wherever
      nu_B(k) = 0, a >= -1e-12 and a_k <= (audited cap) nu_B(k) at EVERY class; and the
      EXACT-MATCHING ones (the LP optima of PART 1 and PART 2) also satisfy
      max_{i in H} |(Aa)_i - t_i| <= 1e-9 max_i t_i and |sum_k a_k - 1| <= 1e-9.
      The truncated explicit rule is NOT exact-matching -- its a <= Q^{1/2} truncation leaves a positive
      deficit on the Y-smooth rows too, which is MEASURED, not scored.  A failure voids the
      (cell, D) or (cell, Y).
  C4  CONTROL, the two aggregations agree.  For every assembled weight, the band-level strided row
      sums of the spread weight and the class-level sparse matvec agree on H to 1e-9 * max_i t_i.
      A failure voids the (cell, D) or (cell, Y).
  C5  CONTROL, a mathematical necessity.  Every measured ||t - Ax||_2/||p||_2 is >= the exact
      kappa as computed earlier -- 0.1724 (2^16), 0.1409 (2^18), 0.1055 (2^20) -- up to 1e-3 RELATIVE
      (4-decimal literals, half-ulp 5e-5 < 1e-3*0.1055 = 1.1e-4).  A failure is a bug here, not a
      finding, and VOIDS that cell.
  C6  CONTROL, a PROVED monotonicity.  C*(Q, D) is non-decreasing in D along the
      sweep, to 1e-9 relative.  A failure is a bug and VOIDS that cell.
  C7  CONTROL, a PROVED uniqueness.  In the Y-smooth scheme the incidence between
      H_Y = {Y-smooth squarefree d <= Q} and the classes is unitriangular, so when EVERY class
      present is <= Q the system has the unique solution a = nu_R and therefore
      C*_Y = C_PC = (|vB|/|vR|) max_{a : B_a>0} R_a/B_a exactly.  Checked to 1e-9 RELATIVE at every
      (cell, Y) with no class above Q.  A failure is a bug and voids the (cell, Y).
  P1  PREDICTION.  C* depends essentially on D/Q alone: at each m in 1..7 the three cells'
      C*(Q, Q/2^m) agree within a factor 1.6 (max/min <= 1.6), wherever all three are computed.
  P2  PREDICTION.  C*(Q, D) <= 1.5 for every D <= Q/32 at all three cells -- i.e. C* -> 1 as
      D/Q -> 0, which is what the route needs.
  P3  PREDICTION.  The Y-smooth scheme has NO dead classes among a <= Q^{1/2} at any Y tested, so
      no dead class obstructs (PC) -- unlike the D ~ cQ regime, where
      0, 16, 3, 34, 20, 73 of them were measured.
  P4  PREDICTION.  (PC) is strictly stronger than a bounded C*.  The comparison must be
      like with like: the UNTRUNCATED identity rule a = nu_R (no a <= Q^{1/2} cut, restricted to the
      live classes) is exact on the Y-smooth rows, so its cap
      C_PC = (|vB|/|vR|) max_{a : B_a > 0} R_a/B_a is an upper bound for the optimal cap C*_Y of the
      same row system; the prediction is C_PC >= 2 C*_Y at every (cell, Y) **that has a class above
      Q** -- the only regime with any freedom, since the uniqueness of clause C7 forces equality in the
      other.  The truncated rule's own cap is the max over a <= Q^{1/2}, which is NOT comparable to C*_Y
      (its truncation leaves a deficit on the Y-smooth rows) and is reported separately.
  P5  PREDICTION, two-sided.  At every cell the BEST kappa bound over all constructions here is
      (i) below the trivial 1 (dist_2(t,K_Q) <= ||p||_2 since z in K_Q), so the capped mass-one
      constructions do say something, and (ii) at least 2x the exact kappa of C5, so they are far
      from sharp.

  Scoring is three-way (HOLDS / REFUTED / INCONCLUSIVE), voids propagate, exit code 0 either way.

  DISCLOSED: before registration this file was smoke-run at 2^16 on the reduced grid D in
  {Q/128, Q/64}, Y in {8} with --smoke (no scoring, no result files).  That run found two things that
  changed the draft: the companion truncation is not exact on the Y-smooth
  rows (C3 is scoped accordingly, and the truncated rule is assembled on its live classes only, with the
  blocked mass reported separately), and the kappa bound at D = 2, 4 is already 0.81, 0.85 -- below
  the trivial 1 -- which is why P5 is two-sided instead of one-sided.  P2 and P5 are therefore
  informed at those three grid points; the rest of the grid and the other two cells are not.

    python code/conedual_cstar_sweep.py --cells 16 18 20
"""
from __future__ import annotations
import io, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
NTHREAD = int(os.environ.get("ACF_THREADS", "4"))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = str(NTHREAD)
_SCRATCH = os.environ.get("TEMP") or os.environ.get("TMP") or "."
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(_SCRATCH, "numba_cache_acf"))
os.environ["NUMBA_NUM_THREADS"] = str(NTHREAD)

import numpy as np                                      # noqa: E402
import scipy.sparse as spx                              # noqa: E402
from scipy.optimize import linprog                      # noqa: E402
from numba import njit, prange                          # noqa: E402

sys.path.insert(0, HERE)
from cell_cache import load                             # noqa: E402
from _sieve_shared import primes_upto                   # noqa: E402
from conedual_kappa_explicit_target import band_columns  # noqa: E402

CELLS = (16, 18, 20)
MSWEEP = (7, 6, 5, 4, 3, 2, 1)
PLOGS = (2.0, 3.0, 4.0)
YFIX = (8, 16, 32, 64)
MMAX = 10
LPMAX_LIVE = 120000
DUALMAX_LIVE = 40000
TLIM = 1200.0
# t_1 = R_1/|vR| = 1 EXACTLY (row d = 1 divides every band element), so the unit-mass comparand is
# the exact integer 1, not a rounded literal: the 1e-9 tolerance is about float accumulation.
UNIT = 1.000000000000
B12B_CSTAR = {(16, 8): 4.0326, (18, 8): 4.3890, (20, 8): 3.1114,
              (16, 4): 12.1837, (18, 4): 15.5303, (20, 4): 13.8452}   # the earlier sweeps' values
KAPPA_EXACT = {16: 0.1724, 18: 0.1409, 20: 0.1055}                    # computed earlier (dense nnls)
OUT = io.StringIO()


def say(s=""):
    print(s)
    sys.stdout.flush()
    OUT.write(s + "\n")


@njit(parallel=True, cache=True, nogil=True)
def k_row_weight(wt, rows, out, nthread):
    n = wt.shape[0] - 1
    m = rows.shape[0]
    for w in prange(nthread):
        for i in range(w, m, nthread):
            dd = rows[i]
            s = 0.0
            for j in range(dd, n + 1, dd):
                s += wt[j]
            out[i] = s


def smooth_part(N, Y):
    """kk[n] = prod_{p | n, p <= Y} p, om[n] = its omega, lpf[n] = least prime factor <= Y."""
    kk = np.ones(N + 1, dtype=np.int64)
    om = np.zeros(N + 1, dtype=np.int8)
    lpf = np.zeros(N + 1, dtype=np.int64)
    for p in primes_upto(Y):
        p = int(p)
        kk[p::p] *= p
        om[p::p] += 1
        lpf[p::p] = p
    return kk, om, lpf


def incidence(u, rowsH, lpf):
    """A[i, k] = 1 iff rowsH[i] divides the class u[k].  Classes are Y-smooth and squarefree, so
    every divisor is a sub-product of the class's own prime slots: enumerate all subsets.  Divisors
    that are not rows (too large) simply fail the lookup.  No Python loop over classes."""
    n = u.size
    P = np.zeros((MMAX, n), dtype=np.int64)
    rest = u.copy()
    for i in range(MMAX):
        P[i] = lpf[rest]
        rest = np.where(P[i] > 0, rest // np.maximum(P[i], 1), rest)
    assert np.all(rest == 1), "a class has more than MMAX prime factors"
    rr, cc = [], []
    for mask in range(1 << MMAX):
        vals = np.ones(n, dtype=np.int64)
        okm = np.ones(n, dtype=bool)
        for i in range(MMAX):
            if mask >> i & 1:
                okm &= P[i] > 0
                vals = vals * np.maximum(P[i], 1)
        if not np.any(okm):
            continue
        jj = np.minimum(np.searchsorted(rowsH, vals), rowsH.size - 1)
        hit = okm & (rowsH[jj] == vals)
        if np.any(hit):
            rr.append(jj[hit])
            cc.append(np.nonzero(hit)[0])
    A = spx.coo_matrix((np.ones(sum(x.size for x in rr)),
                        (np.concatenate(rr), np.concatenate(cc))), shape=(rowsH.size, n)).tocsc()
    A.sum_duplicates()
    A.data[:] = 1.0
    return A


def cap_lp(Acols, nuc, rhs, tl=TLIM):
    """min C  s.t.  Acols x = rhs,  0 <= x <= C nuc."""
    n = Acols.shape[1]
    A_ub = spx.hstack([spx.eye(n, format="csc"), spx.csc_matrix(-nuc.reshape(-1, 1))], format="csc")
    cvec = np.zeros(n + 1)
    cvec[-1] = 1.0
    t0 = time.time()
    r = linprog(cvec, A_ub=A_ub, b_ub=np.zeros(n),
                A_eq=spx.hstack([Acols, spx.csc_matrix((rhs.size, 1))], format="csc"), b_eq=rhs,
                bounds=[(0, None)] * n + [(0, None)], method="highs",
                options={"time_limit": tl, "presolve": True})
    return r, time.time() - t0


def dual_lp(Acols, nuc, tH, tl=TLIM):
    """max <y,t>  s.t.  s_k >= (A^T y)_k ,  s >= 0 ,  sum_k nuc_k s_k <= 1   (= the polar value)."""
    n = Acols.shape[1]
    nh = tH.size
    A_ub = spx.vstack([spx.hstack([Acols.T.tocsc(), -spx.eye(n, format="csc")], format="csc"),
                       spx.hstack([spx.csc_matrix((1, nh)), spx.csc_matrix(nuc.reshape(1, -1))],
                                  format="csc")], format="csc")
    b_ub = np.concatenate([np.zeros(n), [1.0]])
    cvec = np.concatenate([-tH, np.zeros(n)])
    t0 = time.time()
    r = linprog(cvec, A_ub=A_ub, b_ub=b_ub, bounds=[(None, None)] * nh + [(0, None)] * n,
                method="highs", options={"time_limit": tl, "presolve": True})
    return r, time.time() - t0


class Cell:
    """Everything that does not depend on D or Y."""

    def __init__(self, e):
        S = load(e)
        self.e = e
        self.N, self.Q = int(S["N"]), int(S["Q"])
        self.rows = np.asarray(S["rows"]).astype(np.int64)
        self.cb = np.asarray(S["cb"]).astype(np.float64)
        self.cr = np.asarray(S["cr"]).astype(np.float64)
        self.nB, self.nR = int(S["nB"]), int(S["nR"])
        self.p2 = float(np.linalg.norm(np.asarray(S["p"]).astype(np.float64)))
        self.t_all = self.cr / self.nR
        self.tmax = float(np.max(self.t_all))
        code, pos, vB = band_columns(e, S)
        self.vB = vB.astype(np.int64)
        self.vR = np.nonzero(code == 2)[0].astype(np.int64)
        self.wt = np.zeros(self.N + 1, dtype=np.float64)
        self.out = np.zeros(self.rows.size, dtype=np.float64)

    def band_rows(self, x_vB):
        """All-row aggregate of the band weight x (supported on vB), by the strided kernel."""
        self.wt[:] = 0.0
        self.wt[self.vB] = x_vB
        k_row_weight(self.wt, self.rows, self.out, NTHREAD)
        return self.out.copy()


def assemble(C, a_cls, nuB, A, tH, inv_B, Bc, matched, exact):
    """Audit one class measure: class-level equalities/box/support, the band-level spread, the
    all-row residual and the kappa bound it gives.  `exact` says whether the weight is supposed to
    match t on H (the LP optima are; the truncated rule is not, and its deficit is measured)."""
    liv = nuB > 0
    cap = float(np.max(a_cls[liv] / nuB[liv])) if np.any(liv) else float("nan")
    supp = float(np.max(a_cls[~liv])) if np.any(~liv) else 0.0
    agg = np.asarray(A.dot(a_cls)).ravel()
    x_vB = a_cls[inv_B] / Bc[inv_B]
    out = C.band_rows(x_vB)
    r = C.t_all - out
    kap = float(np.linalg.norm(r)) / C.p2
    eq = float(np.max(np.abs(agg - tH))) / C.tmax
    xchk = float(np.max(np.abs(out[matched] - agg))) / C.tmax
    mass = float(a_cls.sum())
    ok = bool(supp <= 1e-15 and np.max(a_cls - cap * nuB) <= 1e-12 * max(cap, 1.0)
              and a_cls.min() >= -1e-12 and xchk <= 1e-9
              and (not exact or (eq <= 1e-9 and abs(mass - UNIT) <= 1e-9)))
    return dict(cap=cap, eq=eq, supp=supp, mass=mass, minv=float(a_cls.min()), xchk=xchk,
                resid=float(np.linalg.norm(r)), kappa_ub=kap, ok=ok, exact=bool(exact),
                r_unit=float(r[0]), nsupp=int(np.sum(a_cls > 1e-15)))


def part1(C, res):
    """C*(Q, D) for D = Q/2^m and polylog D."""
    Q, N = C.Q, C.N
    L = np.log(Q)
    ds = set(max(2, Q >> m) for m in MSWEEP)
    plog = {}
    for a in PLOGS:
        v = int(L ** a)
        plog["L^%g" % a] = v
        if 2 <= v <= Q:
            ds.add(v)
    say("  PART 1  D-sweep.  L = log Q = %.4f, polylog D: %s  (those > Q are skipped)"
        % (L, ", ".join("%s=%d" % (k, v) for k, v in plog.items())))
    d = {}
    for D in sorted(ds):
        t0 = time.time()
        kk, om, lpf = smooth_part(N, D)
        kb, kr = kk[C.vB], kk[C.vR]
        u, inv = np.unique(np.concatenate([kb, kr]), return_inverse=True)
        inv_B = inv[:kb.size]
        Bc = np.bincount(inv_B, minlength=u.size).astype(np.float64)
        Rc = np.bincount(inv[kb.size:], minlength=u.size).astype(np.float64)
        nuB, nuR = Bc / C.nB, Rc / C.nR
        matched = np.nonzero(C.rows <= D)[0]
        rowsH = C.rows[matched]
        tH = C.t_all[matched]
        A_all = incidence(u, rowsH, lpf)
        keep = np.nonzero(nuB > 0)[0]
        nuk = nuB[keep]
        tag = "D=%d (Q/%.3g)" % (D, Q / D)
        ent = dict(D=int(D), ratio=D / Q, classes=int(u.size), live=int(keep.size),
                   nH=int(rowsH.size), dead_head=int(np.sum((nuB == 0) & (nuR > 0)
                                                            & np.isin(u, rowsH))))
        if keep.size > LPMAX_LIVE:
            ent["skipped"] = "live classes %d > %d" % (keep.size, LPMAX_LIVE)
            say("    %-18s classes=%-7d live=%-7d |H|=%-5d  SKIPPED (%s)"
                % (tag, u.size, keep.size, rowsH.size, ent["skipped"]))
        else:
            rl, secs = cap_lp(A_all[:, keep].tocsc(), nuk, tH)
            ent["status"] = int(rl.status)
            ent["secs"] = secs
            if rl.success:
                a_cls = np.zeros(u.size)
                a_cls[keep] = rl.x[:keep.size]
                au = assemble(C, a_cls, nuB, A_all, tH, inv_B, Bc, matched, True)
                ent["C"] = float(rl.x[-1])
                ent["audit"] = au
                if keep.size <= DUALMAX_LIVE:
                    rd, ds2 = dual_lp(A_all[:, keep].tocsc(), nuk, tH)
                    if rd.success:
                        y = rd.x[:rowsH.size]
                        s = np.asarray(A_all[:, keep].T.dot(y)).ravel()
                        den = float(np.sum(nuk * np.maximum(s, 0.0)))
                        pol = float(np.dot(y, tH) / den) if den > 0 else None
                        ent["polar"] = pol
                        ent["dual_gap"] = (abs(pol - ent["C"]) / ent["C"]) if pol else None
                        ent["dual_secs"] = ds2
                say("    %-18s classes=%-7d live=%-7d |H|=%-5d dead=%-4d  C* = %.6f  polar gap %s"
                    "  kappa_ub %.4f  audit %s  (%.1f s)"
                    % (tag, u.size, keep.size, rowsH.size, ent["dead_head"], ent["C"],
                       ("%.1e" % ent["dual_gap"]) if ent.get("dual_gap") is not None else "-",
                       au["kappa_ub"], au["ok"], secs))
            else:
                say("    %-18s classes=%-7d live=%-7d |H|=%-5d  LP status %d  (%.1f s)"
                    % (tag, u.size, keep.size, rowsH.size, rl.status, secs))
        ent["total_secs"] = time.time() - t0
        d["%d" % D] = ent
        del kk, om, lpf, A_all
    res["part1"] = d
    res["plog"] = plog


def part2(C, res):
    """The small-prime-part scheme: H = all Y-smooth rows d <= Q, classes = a_Y(n)."""
    Q, N = C.Q, C.N
    L = np.log(Q)
    ys = sorted(set(list(YFIX) + [int(L ** a) for a in PLOGS if 2 <= int(L ** a) <= Q]))
    rootQ = int(np.floor(np.sqrt(Q)))
    say("  PART 2  Y-smooth scheme.  Q^{1/2} = %d (the truncated rule keeps the classes a <= Q^{1/2}),"
        "  Y grid %s" % (rootQ, ys))
    d = {}
    for Y in ys:
        t0 = time.time()
        kk, om, lpf = smooth_part(N, Y)
        kb, kr = kk[C.vB], kk[C.vR]
        u, inv = np.unique(np.concatenate([kb, kr]), return_inverse=True)
        inv_B = inv[:kb.size]
        Bc = np.bincount(inv_B, minlength=u.size).astype(np.float64)
        Rc = np.bincount(inv[kb.size:], minlength=u.size).astype(np.float64)
        nuB, nuR = Bc / C.nB, Rc / C.nR
        matched = np.nonzero(kk[C.rows] == C.rows)[0]          # the Y-SMOOTH rows d <= Q, all of them
        rowsH = C.rows[matched]
        tH = C.t_all[matched]
        A_all = incidence(u, rowsH, lpf)
        small = u <= rootQ
        dead = (nuB == 0) & (nuR > 0)
        ent = dict(Y=int(Y), classes=int(u.size), nH=int(rowsH.size),
                   nrows=int(C.rows.size), small=int(small.sum()),
                   dead_all=int(dead.sum()), dead_small=int(np.sum(dead & small)),
                   omega_max=int(om[u].max()), maxclass=int(u.max()),
                   above_Q=int(np.sum(u > Q)))
        # ---- the truncated explicit rule: a = nu_R on the classes a <= Q^{1/2} ----------------------
        # x_l = R_a/(R B_a) places mass only where B_a > 0, so the realised measure is nu_R on
        # {a <= Q^{1/2}, B_a > 0}; the odd mass of a dead class is lost, and that loss -- not a
        # support violation -- is what (PC)'s failure costs.  Both are reported.
        pos = small & (nuB > 0)
        a_pc = np.where(pos, nuR, 0.0)
        ent["delta"] = float(nuR[~small].sum())
        ent["cap_pc"] = float(np.max(nuR[pos] / nuB[pos]))
        ent["cap_pc_argmax"] = int(u[pos][int(np.argmax(nuR[pos] / nuB[pos]))])
        ent["ratio_max"] = float(np.max((Rc[pos] / Bc[pos])))
        ent["pc_blocked"] = int(np.sum(dead & small))
        ent["dead_mass"] = float(nuR[dead & small].sum())
        allpos = nuB > 0
        ent["cap_pc_all"] = float(np.max(nuR[allpos] / nuB[allpos]))     # the untruncated rule
        ent["cap_pc_all_argmax"] = int(u[allpos][int(np.argmax(nuR[allpos] / nuB[allpos]))])
        ent["dead_mass_all"] = float(nuR[dead].sum())
        au = assemble(C, a_pc, nuB, A_all, tH, inv_B, Bc, matched, False)
        ent["pc"] = au
        say("    Y=%-5d classes=%-6d (max %d, above Q: %d)  |H_Y|=%-5d of %-5d rows  omega<=%d"
            "  dead: %d (a<=Q^1/2: %d, lost odd mass %.3e)"
            % (Y, u.size, ent["maxclass"], ent["above_Q"], rowsH.size, C.rows.size,
               ent["omega_max"], ent["dead_all"], ent["dead_small"], ent["dead_mass"]))
        say("          truncated rule: cap = (|vB|/|vR|) max_{a<=Q^1/2} R_a/B_a = %.4f at a=%d"
            "   max R_a/B_a = %.4f   delta = %.6f   mass %.6f"
            % (ent["cap_pc"], ent["cap_pc_argmax"], ent["ratio_max"], ent["delta"], au["mass"]))
        say("                    deficit on the Y-smooth rows %.3e (rel. to max t; 0 would be"
            " exactness)   kappa_ub %.4f   audit %s"
            % (au["eq"], au["kappa_ub"], au["ok"]))
        # ---- the optimal cap for the same row system -------------------------------------
        keep = np.nonzero(nuB > 0)[0]
        rl, secs = cap_lp(A_all[:, keep].tocsc(), nuB[keep], tH)
        ent["opt_status"] = int(rl.status)
        ent["opt_secs"] = secs
        if rl.success:
            a_cls = np.zeros(u.size)
            a_cls[keep] = rl.x[:keep.size]
            au = assemble(C, a_cls, nuB, A_all, tH, inv_B, Bc, matched, True)
            ent["C_opt"] = float(rl.x[-1])
            ent["opt"] = au
            say("          optimal for the same rows: C*_Y = %.6f   kappa_ub %.4f   audit %s"
                "   untruncated C_PC = %.4f at a=%d (/C*_Y = %.2f)   (%.1f s)"
                % (ent["C_opt"], au["kappa_ub"], au["ok"], ent["cap_pc_all"],
                   ent["cap_pc_all_argmax"], ent["cap_pc_all"] / ent["C_opt"], secs))
        else:
            say("          optimal LP status %d  (%.1f s)" % (rl.status, secs))
        ent["total_secs"] = time.time() - t0
        d["%d" % Y] = ent
        del kk, om, lpf, A_all
    res["part2"] = d
    res["rootQ"] = rootQ


def cell(e, res):
    C = Cell(e)
    d = dict(e=e, N=C.N, Q=C.Q, nB=C.nB, nR=C.nR, p2=C.p2, nrows=int(C.rows.size),
             kappa_exact=KAPPA_EXACT.get(e))
    say("")
    say("=" * 118)
    say("cell 2^%d  N=%d Q=%d  |vB|=%d |vR|=%d  rows=%d  ||p||_2=%.6f  exact kappa %s"
        % (e, C.N, C.Q, C.nB, C.nR, C.rows.size, C.p2, KAPPA_EXACT.get(e)))
    say("=" * 118)
    part1(C, d)
    part2(C, d)
    res[e] = d


def _ck(v, lit, tol=1e-4):
    return v is not None and abs(v - lit) <= tol * lit


def score(res):
    v, void = {}, set()
    say("")
    say("=" * 118)
    say("REGISTERED CLAUSES")
    say("=" * 118)
    c1 = {}
    for e in sorted(res):
        for den in (8, 4):
            lit = B12B_CSTAR.get((e, den))
            D = res[e]["Q"] // den
            got = res[e].get("part1", {}).get("%d" % D, {}).get("C")
            if lit is None:
                continue
            ok = _ck(got, lit)
            c1["2^%d Q/%d" % (e, den)] = dict(got=got, lit=lit, holds=bool(ok))
            if not ok:
                void.add(e)
    v["C1"] = dict(cells=c1, verdict=("INCONCLUSIVE" if not c1 else
                                      "HOLDS" if all(x["holds"] for x in c1.values()) else "REFUTED"))
    say("  C1 C* reproduces the earlier sweep at Q/8 and Q/4: %s  (%s)"
        % (v["C1"]["verdict"], ", ".join("%s: %s vs %.4f" % (
            k, ("%.5f" % c1[k]["got"]) if c1[k]["got"] else "-", c1[k]["lit"]) for k in c1)))
    bad2 = [["2^%d" % e, k, x["dual_gap"]] for e in sorted(res)
            for k, x in res[e].get("part1", {}).items()
            if x.get("dual_gap") is not None and x["dual_gap"] > 1e-6]
    n2 = sum(1 for e in res for k, x in res[e].get("part1", {}).items()
             if x.get("dual_gap") is not None)
    v["C2"] = dict(bad=bad2, n=n2, verdict=("INCONCLUSIVE" if not n2 else
                                            "HOLDS" if not bad2 else "REFUTED"))
    say("  C2 polar value equals C* (1e-6 rel) where computed: %s  (%d checked%s)"
        % (v["C2"]["verdict"], n2, "" if not bad2 else ", bad %s" % bad2))
    for b in bad2:
        void.add(int(b[0][2:]))
    au_all = []
    for e in sorted(res):
        for k, x in res[e].get("part1", {}).items():
            if x.get("audit"):
                au_all.append(["2^%d" % e, "D=%s" % k, x["audit"]])
        for k, x in res[e].get("part2", {}).items():
            for nm in ("pc", "opt"):
                if x.get(nm):
                    au_all.append(["2^%d" % e, "Y=%s %s" % (k, nm), x[nm]])
    bad3 = [a for a in au_all if not a[2]["ok"]]
    v["C3"] = dict(n=len(au_all), bad=[[a[0], a[1]] for a in bad3],
                   verdict=("INCONCLUSIVE" if not au_all else "HOLDS" if not bad3 else "REFUTED"))
    say("  C3 assembled-weight audit (equalities, support, box, unit mass): %s  (%d weights%s)"
        % (v["C3"]["verdict"], len(au_all), "" if not bad3 else ", bad %s" % v["C3"]["bad"]))
    bad4 = [[a[0], a[1], a[2]["xchk"]] for a in au_all if a[2]["xchk"] > 1e-9]
    v["C4"] = dict(n=len(au_all), bad=bad4,
                   verdict=("INCONCLUSIVE" if not au_all else "HOLDS" if not bad4 else "REFUTED"))
    say("  C4 band-level and class-level aggregation agree (1e-9): %s  (worst %.1e of %d)"
        % (v["C4"]["verdict"], max([a[2]["xchk"] for a in au_all] or [float("nan")]), len(au_all)))
    for b in bad3 + bad4:
        void.add(int(b[0][2:]))
    bad5 = []
    for a in au_all:
        e = int(a[0][2:])
        ke = KAPPA_EXACT.get(e)
        if ke and a[2]["kappa_ub"] < ke * (1 - 1e-3):
            bad5.append([a[0], a[1], a[2]["kappa_ub"], ke])
    v["C5"] = dict(bad=bad5, verdict=("INCONCLUSIVE" if not au_all else
                                      "HOLDS" if not bad5 else "REFUTED"))
    say("  C5 every kappa bound >= the exact kappa: %s%s"
        % (v["C5"]["verdict"], "" if not bad5 else "  %s" % bad5))
    for b in bad5:
        void.add(int(b[0][2:]))
    bad6 = []
    for e in sorted(res):
        p1 = res[e].get("part1", {})
        seq = [(int(k), p1[k]["C"]) for k in sorted(p1, key=lambda z: int(z)) if p1[k].get("C")]
        for (d1, c1v), (d2, c2v) in zip(seq, seq[1:]):
            if c2v < c1v * (1 - 1e-9):
                bad6.append(["2^%d" % e, d1, d2, c1v, c2v])
    v["C6"] = dict(bad=bad6, verdict=("INCONCLUSIVE" if not res else
                                      "HOLDS" if not bad6 else "REFUTED"))
    say("  C6 C*(Q,D) non-decreasing in D: %s%s"
        % (v["C6"]["verdict"], "" if not bad6 else "  %s" % bad6))
    for b in bad6:
        void.add(int(b[0][2:]))

    c7 = {}
    for e in sorted(res):
        for k, x in res[e].get("part2", {}).items():
            if x.get("above_Q") == 0 and x.get("C_opt"):
                ok = abs(x["cap_pc_all"] - x["C_opt"]) <= 1e-9 * x["C_opt"]
                c7["2^%d Y=%s" % (e, k)] = dict(copt=x["C_opt"], cpc=x["cap_pc_all"], holds=bool(ok))
                if not ok:
                    void.add(e)
    v["C7"] = dict(cells=c7, verdict=("INCONCLUSIVE" if not c7 else
                                      "HOLDS" if all(x["holds"] for x in c7.values()) else "REFUTED"))
    say("  C7 uniqueness where no class exceeds Q: %s  (%s)"
        % (v["C7"]["verdict"], ", ".join("%s: %.6f vs %.6f" % (k, c7[k]["copt"], c7[k]["cpc"])
                                         for k in c7) or "no such (cell, Y)"))

    live = [e for e in sorted(res) if e not in void]
    p1d = {}
    for m in MSWEEP:
        vals = []
        for e in live:
            D = max(2, res[e]["Q"] >> m)
            c = res[e].get("part1", {}).get("%d" % D, {}).get("C")
            if c:
                vals.append(c)
        if len(vals) == len([e for e in live]) and len(vals) >= 2:
            p1d["Q/2^%d" % m] = dict(vals=vals, spread=max(vals) / min(vals),
                                     holds=bool(max(vals) / min(vals) <= 1.6))
    v["P1"] = dict(cells=p1d, verdict=("INCONCLUSIVE" if not p1d else
                                       "HOLDS" if all(x["holds"] for x in p1d.values()) else "REFUTED"))
    say("  P1 C* depends on D/Q alone (spread <= 1.6x): %s  (%s)"
        % (v["P1"]["verdict"], ", ".join("%s: %.2f" % (k, p1d[k]["spread"]) for k in p1d)))
    p2d = {}
    for e in live:
        for k, x in res[e].get("part1", {}).items():
            if x.get("C") and int(k) <= res[e]["Q"] / 32:
                p2d["2^%d D=%s" % (e, k)] = bool(x["C"] <= 1.5)
    v["P2"] = dict(cells=p2d, verdict=("INCONCLUSIVE" if not p2d else
                                       "HOLDS" if all(p2d.values()) else "REFUTED"))
    say("  P2 C* <= 1.5 for D <= Q/32: %s  (%d of %d)"
        % (v["P2"]["verdict"], sum(p2d.values()), len(p2d)))
    p3d = {}
    for e in live:
        for k, x in res[e].get("part2", {}).items():
            p3d["2^%d Y=%s" % (e, k)] = dict(dead=x["dead_small"], holds=bool(x["dead_small"] == 0))
    v["P3"] = dict(cells=p3d, verdict=("INCONCLUSIVE" if not p3d else
                                       "HOLDS" if all(x["holds"] for x in p3d.values()) else "REFUTED"))
    say("  P3 no dead classes among a <= Q^{1/2}: %s  (%s)"
        % (v["P3"]["verdict"], ", ".join("%s: %d" % (k, p3d[k]["dead"]) for k in p3d)))
    p4d = {}
    for e in live:
        for k, x in res[e].get("part2", {}).items():
            if x.get("C_opt") and x.get("above_Q", 0) > 0:
                p4d["2^%d Y=%s" % (e, k)] = dict(r=x["cap_pc_all"] / x["C_opt"],
                                                 holds=bool(x["cap_pc_all"] >= 2.0 * x["C_opt"]))
    v["P4"] = dict(cells=p4d, verdict=("INCONCLUSIVE" if not p4d else
                                       "HOLDS" if all(x["holds"] for x in p4d.values()) else "REFUTED"))
    say("  P4 untruncated C_PC >= 2x the optimal cap where a class exceeds Q: %s  (%s)"
        % (v["P4"]["verdict"], ", ".join("%s: %.2f" % (k, p4d[k]["r"]) for k in p4d)))
    p5d = {}
    for e in live:
        ks = [a[2]["kappa_ub"] for a in au_all if a[0] == "2^%d" % e]
        ke = KAPPA_EXACT.get(e)
        if ks and ke:
            p5d["2^%d" % e] = dict(best=min(ks), ratio=min(ks) / ke,
                                   holds=bool(min(ks) < 1.0 and min(ks) >= 2.0 * ke))
    v["P5"] = dict(cells=p5d, verdict=("INCONCLUSIVE" if not p5d else
                                       "HOLDS" if all(x["holds"] for x in p5d.values()) else "REFUTED"))
    say("  P5 best kappa bound in (2 kappa_exact, 1): %s  (%s)"
        % (v["P5"]["verdict"], ", ".join("%s: %.4f = %.2f x exact" % (k, p5d[k]["best"], p5d[k]["ratio"])
                                         for k in p5d)))
    if void:
        say("  VOID cells: %s -- their P clauses are NOT SCORED" % sorted(void))
    v["void_cells"] = sorted(void)
    return v


def dump(res, clauses=None, tag=""):
    with open(os.path.join(RES, "conedual_cstar_sweep%s.json" % tag), "w") as f:
        json.dump(dict(cells={str(k): res[k] for k in res}, clauses=clauses or {}), f, indent=1,
                  default=float)
    with open(os.path.join(RES, "conedual_cstar_sweep%s.txt" % tag), "w", encoding="utf-8") as f:
        f.write(OUT.getvalue())


def main(argv):
    global MSWEEP, PLOGS, YFIX
    cells = CELLS
    if "--cells" in argv:
        i = argv.index("--cells") + 1
        vals = []
        while i < len(argv) and argv[i].isdigit():
            vals.append(int(argv[i]))
            i += 1
        cells = tuple(vals)
    bad = [e for e in cells if e > 24]
    if bad and "--force" not in argv:
        say("refusing cells %s above 2^24 without --force" % bad)
        cells = tuple(e for e in cells if e <= 24)
    smoke = "--smoke" in argv
    if smoke:
        MSWEEP, PLOGS, YFIX = (7, 6), (), (8,)
    tag = argv[argv.index("--tag") + 1] if "--tag" in argv else ""
    say("conedual_cstar_sweep -- C*(Q,D) over D, and the Y-smooth (small-prime-part) scheme")
    say("cells %s  threads=%d  m-sweep %s  polylog %s  Y grid %s%s"
        % (list(cells), NTHREAD, list(MSWEEP), list(PLOGS), list(YFIX), "  SMOKE" if smoke else ""))
    res = {}
    for e in cells:
        try:
            cell(e, res)
        except Exception as ex:
            say("  cell 2^%d FAILED: %r" % (e, ex))
            res[e] = dict(e=e, error=repr(ex))
        if not smoke:
            dump(res, tag=tag)
    if smoke:
        say("")
        say("SMOKE run: not scored, no result files written")
        return 0
    v = score(res)
    dump(res, v, tag=tag)
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.exit(main(sys.argv[1:]))

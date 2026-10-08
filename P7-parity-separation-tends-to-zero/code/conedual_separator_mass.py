# -*- coding: utf-8 -*-
r"""conedual_separator_mass -- the optimal separator's fixed-row mass (direction 1).

Supporting computation; not used in the paper's proofs or tables.

WHERE THIS COMES FROM

Direction 1 of the separator-mass question.  With `u_Q = w*/||w*||_2` the optimal separator
(so kappa = <u_Q, t>/||p||_2 = dist_2(t, K_Q)/||p||_2), `v = -mu/phi`, `c(Q) = N/(4 R_1 L^2)` and
`r_B = ||p - c(Q) v||_2/||p||_2 -> 0`:

    (I)   kappa <= r_B + (c(Q)/||p||_2) (-<v, u_Q>)_+
    (II)  kappa <= r_B + (c(Q)/||p||_2) (||v_{<=D}||_2 sqrt(H_D(Q)) + ||v_{>D}||_2),
          H_D(Q) = sum_{d <= D} u_{Q,d}^2 .

Both are proved elsewhere; this file proves neither.  So EITHER `<v, u_Q> >= 0` eventually, OR `H_D(Q) -> 0` for
every fixed D, gives kappa -> 0 -- outside (S1).  This file measures both sides with certified
intervals, and measures the contact (tight) structure the nonbipartite tight-semiprime mechanism predicts.

WHAT IS COMPUTED

  A. EXACT DUAL at the small cells.  The band is rebuilt (`cell_cache` stores only the per-row
     arrays, not vB; its own docstring says a script needing the band builds `_build`-equivalent
     code -- and the rebuilt rows/cb/cr are asserted equal to the cached cell, clause C0).  The
     distinct columns A_l = (1_{d|l})_d, l in vB, are formed by one numba prange kernel; NNLS gives
     the projection P_K t, hence the EXACT u_Q = (t - P_K t)/||t - P_K t||_2 and kappa.
  B. STORED SEPARATORS at the large cells (`conedual_colgen_sep_2e{22,24,26,28}*.npz`).  Each is
     repaired exactly as `conedual_s1_evidence` does -- n = rint(w*DEN), n[0] := -M with
     M = max_{l in vB} theta_n(l), so theta_W <= 0 on vB in exact int64 -- and
       kappa_2(W) = J(W)/(||W||_2 ||p||_2) = (sum_{d>=2} n_d cr_d - M |vR|) / sqrt(sum n^2 * sum (cb-cr)^2)
     is an EXACT ratio of integers under one square root (|vR| and DEN cancel).
  C. CERTIFIED INTERVALS.  For w in K_Q^o with ||w||_2 = 1, <w,t> <= D*<w,u_Q> where D* = dist(t,K_Q)
     (because <w, P_K t> <= 0), so <u_cert, u_Q> >= L/D* >= L/U and
       ||u_cert - u_Q||_2 <= sqrt(2 (1 - L/U)) =: eps
     with L = kappa_2(W)||p||_2 (certified here) and U = kappa's recorded upper endpoint.  Then
     |<v,u_Q> - <v,u_cert>| <= ||v||_2 eps and |sqrt(H_D) - ||u_cert,<=D||_2| <= eps.
  D. THE MECHANISM.  A band element l = p q with p, q prime and l > Q has divisor set {1, p, q}
     inside [1,Q], so a TIGHT such column gives u_p + u_q = -u_1 exactly.  The graph on primes with
     an edge per tight semiprime column: its components, bipartiteness, size k, and the predicted
     u_p = -u_1/2, |u_1| <= 2/sqrt(k) on a nonbipartite component.

INPUTS CITED, NOT RECOMPUTED
  kappa enclosures, as computed earlier and carried here as comparands:
    2^16 0.1723894, 2^18 0.1408605, 2^20 0.1055427, 2^22 0.0763312 (printed exact);
    2^24 [52741,52793]/10^6, 2^26 [34503,34805]/10^6, 2^28 [23782/10^6, 248444/10^7].
  The K_Q^o certification of the four stored separators (exact M, J, S in int64) from
  `results/conedual_s1_evidence.txt` (an earlier evidence run): clause C3 below matches this file's
  own integers against those, integer against integer.

CLAUSES.

  C0  CONTROL.  At every cell the rebuilt band reproduces the cached cell exactly: rows equal,
      cb and cr equal entry for entry, |vB| = nB and |vR| = nR.  CAN FAIL; voids that cell.
  C1  CONTROL (exact-dual cells).  kappa from NNLS equals the carried comparand to a RELATIVE
      1e-5 (the literals carry 7 significant digits, so their precision exceeds the tolerance).
      CAN FAIL; voids the exact-dual part of that cell.
  C2  CONTROL (stored separators).  kappa_2(W) recomputed here lies in
      [npz bound - 1e-7, U(cell)] -- stored float against recomputed float, same quantity.
      CAN FAIL; that separator yields no interval.
  C3  CONTROL (stored separators).  This file's exact integers (M, J_int, S_int) equal the integers
      printed by conedual_s1_evidence for the same file -- integer against integer.
      CAN FAIL; that separator yields no interval.
  C4  CONTROL.  Where this file sweeps [1,N] itself, max_{l in vB}(theta_n(l) - M) == 0 in exact
      int64, so the repaired W is in K_Q^o with no floating slack.  CAN FAIL; that separator yields
      no interval.  (At the cells not swept here, C3's integer match to the already-certified
      s1_evidence run carries the membership, and the output says which cells are which.)

  P1  PREDICTION.  u_{Q,1} < 0 at every scored cell, and u_{Q,p} > 0 for every prime p <= 30 at the
      exact-dual cells.  WHY: 1 divides every column, so u_1 carries the whole cone constraint, and
      the tight-semiprime relation u_p + u_q = -u_1 then puts the primes on the other side.
      MEANING -- holds: the sign structure the mechanism assumes is the real one.
  P2  THE CLAUSE THAT DECIDES DIRECTION 1.  <v, u_Q> > 0 at every scored cell (point value from the
      certified vector), AND its certified lower endpoint <v,u_cert> - ||v||_2 eps > 0 at every cell
      with eps <= 0.01.  PREDICTED: <v,u_cert> in [0.5, 3] at every cell.
      WHY: v_1 = -1 and v_p = +1/(p-1); with u_1 < 0 and u_p = -u_1/2 > 0 both terms of
      -u_1 + sum_p u_p/(p-1) are positive, and the composite tail is O(sum 1/(phi(d))) smaller.
      MEANING -- holds: inequality (I) already gives kappa <= r_B at the computed cells, and
      direction 1 needs only "eventually"; fails: (II)/H_D is the only route and P3 decides it.
  P3  PREDICTION.  H_D(Q) is strictly decreasing in Q for each D in {1, 6, 30, 210, 1000}, over the
      cells whose eps <= 0.01.  PREDICTED H_1 * pi(Q) in [0.1, 4] at every such cell (4 is the
      proved ceiling when the tight-semiprime component covers the primes).
      MEANING -- holds: H_D -> 0 is the live mechanism and its rate is 1/pi(Q);
      fails: the fixed-row mass does not decay and direction 1 is dead on this route.
  P4  PREDICTION (exact-dual cells).  u_{Q,p} / (-u_{Q,1}/2) in [0.8, 1.2] for every prime p <= 30.
      MEANING -- holds: the nonbipartite tight-semiprime mechanism is operative, so the conditional
      theorem behind this route has its hypothesis at the computed cells;
      fails: the contact structure is different and that theorem is vacuous here.
  P5  PREDICTION (exact-dual cells).  The tight-semiprime graph has a nonbipartite component
      covering at least half of the primes <= Q, and the proved inequality u_{Q,1}^2 <= 4/k holds
      with the measured k.  (The second half CANNOT fail if P4's mechanism is exact; it is scored
      because a numerical tightness tolerance is involved.)

  NOT REGISTERED: anything about (S1), about sigma_min, about Goldbach, or any asymptotic claim.
  The inequalities (I) and (II) are proved elsewhere; nothing here is evidence for them.

COMPUTE RULES (the project's compute rules, both sections read before this file was written)
  cells are LOADED (`from cell_cache import load`) and the rebuilt band is checked against the cache
  (C0); no Python loop over rows, pairs or band elements -- numba prange kernels for the column
  matrix, the theta sweep and the row sweep, numpy for the glue; the only Python loops are over
  cells, over the primes <= Q and over the five D.  Threads capped (<= 8).  Cells 2^20 and above
  want a machine with more cores.

    python code/conedual_separator_mass.py 16 18          # the two smallest cells
    python code/conedual_separator_mass.py 20 22 24 26 28
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
from scipy.optimize import nnls  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import connected_components, shortest_path  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _sieve_shared import mu_upto, primes_upto  # noqa: E402
from cell_cache import load as load_cell  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "conedual_separator_mass.txt")
OUTJ = os.path.join(RES, "conedual_separator_mass.json")

DEN = 10 ** 12
DLIST = (1, 6, 30, 210, 1000)
EXACT_CELLS = (16, 18, 20)          # NNLS is affordable here
SWEEP_CELLS = (22, 24)              # own int64 sweep for K_Q^o; 26/28 cite s1_evidence
TIGHT_TOL = 1e-9                    # relative to max|theta| over the columns

# ---- carried comparands, computed earlier -------------------------------------------------------
KAPPA_LO = {16: 0.1723894, 18: 0.1408605, 20: 0.1055427, 22: 0.0763312,
            24: 52741 / 10 ** 6, 26: 34503 / 10 ** 6, 28: 23782 / 10 ** 6}
KAPPA_UP = {16: 0.1723894, 18: 0.1408605, 20: 0.1055427, 22: 0.0763312,
            24: 52793 / 10 ** 6, 26: 34805 / 10 ** 6, 28: 248444 / 10 ** 7}
C1_TOL = 1e-5
C2_SLACK = 1e-7
# ---- carried integers from results/conedual_s1_evidence.txt -----------------------------
S1E = {
    "conedual_colgen_sep_2e22_n0b.npz": dict(M=1000000000001, J=50810068014373698, S=78055974358103479),
    "conedual_colgen_sep_2e24_n0b.npz": dict(M=1000001113767, J=173251720631015907, S=338483562472838081),
    "conedual_colgen_sep_2e26_n0b.npz": dict(M=1000000603706, J=463809741913557021, S=1088810738001557124),
    "conedual_colgen_sep_2e28_n28.npz": dict(M=1000001741817, J=1437321516767375728, S=3795904130585696957),
}
SEPS = {22: ("conedual_colgen_sep_2e22_n0b.npz",), 24: ("conedual_colgen_sep_2e24_n0b.npz",),
        26: ("conedual_colgen_sep_2e26_n0b.npz",), 28: ("conedual_colgen_sep_2e28_n28.npz",)}
# ---- registered prediction bands -------------------------------------------------------------
P2_BAND = (0.5, 3.0)
P3_BAND = (0.1, 4.0)
P4_BAND = (0.8, 1.2)
P5_FRAC = 0.5
EPS_TIGHT = 0.01

_lines = []


def say(s=""):
    _lines.append(s)
    print(s, flush=True)


def flush(payload):
    with io.open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(_lines) + "\n")
    with io.open(OUTJ, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=1, default=float)


# ---------------------------------------------------------------- numba kernels

@njit(parallel=True, cache=True, nogil=True)
def _colmat(rows, pos, nvb, out):
    """out[j, i] = 1 iff rows[i] divides vB[j].  prange over rows; column i is owned by worker i."""
    n = pos.shape[0] - 1
    for i in prange(rows.shape[0]):
        d = rows[i]
        m = d
        while m <= n:
            j = pos[m]
            if j >= 0:
                out[j, i] = 1
            m += d


@njit(parallel=True, cache=True, nogil=True)
def _theta_blocks(acc, rows, vals, nblk):
    """acc[l] += sum_{d | l, d = rows[i], i >= 1} vals[i] for every l <= N (exact int64)."""
    n = acc.shape[0] - 1
    for b in prange(nblk):
        lo = 1 + (n * b) // nblk
        hi = (n * (b + 1)) // nblk
        for i in range(1, rows.shape[0]):
            x = vals[i]
            if x != 0:
                d = rows[i]
                j = ((lo + d - 1) // d) * d
                while j <= hi:
                    acc[j] += x
                    j += d


# ---------------------------------------------------------------- the band (cell_cache stores only rows)

def band(e, S):
    """vB, vR for the cell, built exactly as cell_cache._build does, then checked against the cache."""
    N = 2 ** e
    Q = math.isqrt(N)
    mu = np.array(mu_upto(N), dtype=np.int8)
    rem = np.arange(N + 1, dtype=np.int32)
    om = np.zeros(N + 1, dtype=np.int8)
    for p in primes_upto(Q):                      # a loop over primes, as cell_cache._build does
        p = int(p)
        om[p::p] += 1
        q = p
        while q <= N:
            rem[q::q] //= p
            if q > N // p:
                break
            q *= p
    q1 = int(S["q1"])
    thr = int(S["thr"])
    bd = (mu != 0) & (rem == 1)
    bd[:thr + 1] = False
    del rem
    vB = np.nonzero(bd & (om % 2 == 0))[0].astype(np.int64)
    vR = np.nonzero(bd & (om % 2 == 1))[0].astype(np.int64)
    del bd, om
    return vB, vR, mu


def exact_sum(a, shift=24):
    """sum of an int64 array, exactly (the split keeps both partial sums inside int64)."""
    a = np.asarray(a)
    assert a.dtype == np.int64
    if a.size == 0:
        return 0
    k = 1 << shift
    hi = int(np.abs(a).max()) // k + 1
    assert a.size * hi < 2 ** 62 and a.size * k < 2 ** 62, "exact_sum: raise the shift"
    return int(np.sum(a >> shift)) * k + int(np.sum(a & (k - 1)))


def idot(a, b):
    """exact integer dot product of two short arrays (<= 10^4 entries)."""
    return int(np.asarray(a, dtype=object) @ np.asarray(b, dtype=object))


# ---------------------------------------------------------------- v, phi and the row geometry

def row_geometry(Q, rows, mu):
    phi = np.arange(Q + 1, dtype=np.int64)
    for p in primes_upto(Q):
        p = int(p)
        phi[p::p] -= phi[p::p] // p
    mud = mu[rows].astype(np.float64)
    v = -mud / phi[rows].astype(np.float64)
    return v


def head_tail(v, rows, D):
    m = rows <= D
    return float(np.linalg.norm(v[m])), float(np.linalg.norm(v[~m]))


# ---------------------------------------------------------------- the mechanism: tight semiprimes

def semiprime_graph(reps, tight, Q, spf):
    """Edges {p,q} from the TIGHT columns whose representative is a semiprime p*q > Q.

    reps: a representative band element per distinct column; tight: bool over columns.  Fully
    vectorised: the 2-colouring is `dist mod 2` from one BFS root per component (scipy
    shortest_path over all roots at once), and a component is non-bipartite iff some edge joins
    two equal colours.
    """
    pr = np.asarray(primes_upto(Q), dtype=np.int64)
    npr = int(pr.size)
    pidx = -np.ones(int(pr[-1]) + 2, dtype=np.int64)
    pidx[pr] = np.arange(npr)
    empty = dict(edges=0, k=0, ncomp=0, nonbip_k=0, verts=np.zeros(0, dtype=np.int64), pr=pr)
    ell = reps[tight]
    if ell.size == 0:
        return empty
    lp = spf[ell].astype(np.int64)                        # least prime factor (0 if none <= Q)
    ok = lp > 0
    co = np.where(ok, ell // np.maximum(lp, 1), 0)
    ok &= (co > 1) & (co < pidx.size)
    ok[ok] &= pidx[co[ok]] >= 0                           # the cofactor is itself a prime <= Q
    ok &= (lp * co == ell) & (ell > Q)
    if not ok.any():
        return empty
    ei, ej = pidx[lp[ok]], pidx[co[ok]]
    A = coo_matrix((np.ones(ei.size), (ei, ej)), shape=(npr, npr))
    A = (A + A.T).tocsr()
    ncomp, lab = connected_components(A, directed=False)
    deg = np.diff(A.indptr)
    root_lab, root_idx = np.unique(lab, return_index=True)  # one root per component
    dist = shortest_path(A, method="D", unweighted=True, directed=False, indices=root_idx)
    slot = np.empty(ncomp, dtype=np.int64)
    slot[root_lab] = np.arange(root_lab.size)
    dv = dist[slot[lab], np.arange(npr)]
    colour = np.where(np.isfinite(dv), np.mod(dv, 2.0), -1.0)
    bad = colour[ei] == colour[ej]
    nonbip = np.bincount(lab[ei][bad], minlength=ncomp) > 0 if bad.any() else np.zeros(ncomp, bool)
    sizes = np.bincount(lab[deg > 0], minlength=ncomp)
    cand = np.nonzero(nonbip)[0]
    best = int(cand[np.argmax(sizes[cand])]) if cand.size else -1
    verts = np.nonzero((lab == best) & (deg > 0))[0] if best >= 0 else np.zeros(0, dtype=np.int64)
    return dict(edges=int(ei.size), k=int(sizes.max()) if sizes.size else 0, ncomp=int(ncomp),
                nonbip_k=int(verts.size), verts=pr[verts], pr=pr)


# ---------------------------------------------------------------- one cell

def do_cell(e, store):
    t0 = time.time()
    N = 2 ** e
    S = load_cell(e)
    Q, R1, nB, nR = int(S["Q"]), int(S["nR"]), int(S["nB"]), int(S["nR"])
    rows = np.asarray(S["rows"], dtype=np.int64)
    cb = np.asarray(S["cb"], dtype=np.int64)
    cr = np.asarray(S["cr"], dtype=np.int64)
    Rp = cb - cr                                   # = R1 * p_d, exact integers
    sum_p2 = idot(Rp, Rp)
    normp = math.sqrt(sum_p2) / R1
    L = math.log(Q)
    cQ = N / (4.0 * R1 * L * L)
    rec = dict(e=e, N=N, Q=Q, rows=int(rows.size), nB=nB, nR=nR, normp=normp,
               cQ=cQ, cQ_over_normp=cQ / normp, kappa_lo_cited=KAPPA_LO[e], kappa_up_cited=KAPPA_UP[e])
    say("")
    say("=" * 100)
    say("2^%-2d  Q=%-6d rows=%-5d |vB|=%-9d |vR|=%-9d  ||p||_2=%.9f  c(Q)=%.9f  c(Q)/||p||=%.6f"
        % (e, Q, rows.size, nB, nR, normp, cQ, cQ / normp))
    say("      kappa enclosure carried as comparand: [%.7f, %.7f]" % (KAPPA_LO[e], KAPPA_UP[e]))

    need_band = (e in EXACT_CELLS) or (e in SWEEP_CELLS)
    vB = vR = None
    mu = np.array(mu_upto(max(Q, 3)), dtype=np.int8)
    if need_band:
        vB, vR, mu = band(e, S)
        c0 = (vB.size == nB) and (vR.size == nR)
        pos = -np.ones(N + 1, dtype=np.int64)
        pos[vB] = np.arange(vB.size)
        rec["C0"] = bool(c0)
        say("      C0 band rebuilt: |vB| %d vs nB %d, |vR| %d vs nR %d -> %s"
            % (vB.size, nB, vR.size, nR, "EQUAL" if c0 else "DIFFER"))
        if not c0:
            rec["void"] = "C0"
            store.append(rec)
            return
    v = row_geometry(Q, rows, mu)
    normv = float(np.linalg.norm(v))
    rec["normv"] = normv

    # ---------------- A. exact dual by NNLS
    u_exact = None
    if e in EXACT_CELLS:
        Dm = np.zeros((vB.size, rows.size), dtype=np.int8)
        _colmat(rows, pos, vB.size, Dm)
        Uc, ridx = np.unique(Dm, axis=0, return_index=True)
        del Dm
        M = np.ascontiguousarray(Uc.T.astype(np.float64))      # rows x ndistinct
        tvec = cr.astype(np.float64) / R1
        x, resn = nnls(M, tvec, maxiter=20 * rows.size + 200)
        r = tvec - M @ x
        dist = float(np.linalg.norm(r))
        kap = dist / normp
        c1 = abs(kap - KAPPA_LO[e]) <= C1_TOL * KAPPA_LO[e]
        say("      A exact dual: %d distinct columns of %d; NNLS support %d; dist=%.9f  kappa=%.9f"
            % (Uc.shape[0], vB.size, int(np.count_nonzero(x > 0)), dist, kap))
        say("      C1 kappa vs comparand %.7f: |diff| %.3e <= %.3e -> %s"
            % (KAPPA_LO[e], abs(kap - KAPPA_LO[e]), C1_TOL * KAPPA_LO[e], "HOLDS" if c1 else "FAILS"))
        rec.update(ndistinct=int(Uc.shape[0]), nnls_support=int(np.count_nonzero(x > 0)),
                   dist_exact=dist, kappa_exact=kap, C1=bool(c1))
        if c1:
            u_exact = r / dist
            th = M.T @ u_exact
            mx = float(np.abs(th).max())
            tight = np.abs(th) <= TIGHT_TOL * max(mx, 1.0)
            rec.update(max_theta=mx, n_tight=int(tight.sum()))
            say("      contact set: max|theta_u| over columns %.3e; tight columns %d of %d"
                % (mx, int(tight.sum()), Uc.shape[0]))
            spf = np.zeros(N + 1, dtype=np.int64)
            for p in primes_upto(Q):               # least prime factor, loop over primes
                p = int(p)
                blk = spf[p::p]
                spf[p::p] = np.where(blk == 0, p, blk)
            g = semiprime_graph(vB[ridx], tight, Q, spf)
            del spf
            rec.update(sp_edges=g["edges"], sp_ncomp=g["ncomp"], sp_kmax=g["k"],
                       sp_nonbip_k=g["nonbip_k"], npi=int(g["pr"].size))
            say("      tight semiprime graph: %d edges, %d components, largest %d, largest "
                "NONBIPARTITE component %d of pi(Q)=%d"
                % (g["edges"], g["ncomp"], g["k"], g["nonbip_k"], g["pr"].size))
        else:
            rec["void"] = "C1"
        del M, Uc

    # ---------------- B. stored separators
    u_cert, L_cert, src = (u_exact, None, "exact dual") if u_exact is not None else (None, None, None)
    if u_exact is not None:
        L_cert = rec["kappa_exact"] * normp
    for fname in SEPS.get(e, ()):
        path = os.path.join(RES, fname)
        if not os.path.exists(path):
            say("      separator %-40s ABSENT" % fname)
            continue
        z = np.load(path)
        zr = np.asarray(z["rows"]).astype(np.int64)
        if zr.size != rows.size or not bool((zr == rows).all()):
            say("      separator %-40s VOID: row order differs" % fname)
            continue
        w = np.asarray(z["w"]).astype(np.float64)
        n = np.rint(w * DEN).astype(np.int64)
        n[0] = 0
        cit = S1E.get(fname)
        Mint = None
        if e in SWEEP_CELLS:
            acc = np.zeros(N + 1, dtype=np.int64)
            _theta_blocks(acc, rows, n, 8 * NTHREAD)
            WB = acc[vB]
            Mint = int(WB.max())
            c4 = int((WB - Mint).max()) == 0
            del acc, WB
            say("      separator %-40s  own sweep: M = %d, C4 max(theta_n - M) == 0 -> %s"
                % (fname, Mint, "HOLDS" if c4 else "FAILS"))
        else:
            c4 = None
            Mint = cit["M"] if cit else None
            say("      separator %-40s  M = %d CITED from conedual_s1_evidence (K_Q^o certified there)"
                % (fname, Mint if Mint is not None else -1))
        if Mint is None or c4 is False:
            say("          VOID: no certified M")
            continue
        n[0] = -Mint
        Jint = idot(n[1:], cr[1:]) - Mint * nR
        Sint = Mint * nB - idot(n[1:], cb[1:])
        c3 = bool(cit) and Jint == cit["J"] and Sint == cit["S"] and Mint == cit["M"]
        say("          C3 integers vs s1_evidence: M %s  J %d vs %s  S %d vs %s -> %s"
            % (Mint, Jint, cit["J"] if cit else "-", Sint, cit["S"] if cit else "-",
               "EQUAL" if c3 else "DIFFER"))
        sum_n2 = idot(n, n)
        kap2 = Jint / math.sqrt(float(sum_n2) * float(sum_p2))
        bnd = float(z["bound"]) if "bound" in z.files else float("nan")
        c2 = (kap2 >= bnd - C2_SLACK) and (kap2 <= KAPPA_UP[e] + C2_SLACK) and Jint > 0
        say("          kappa_2(W) = %d/sqrt(%d * %d) = %.9f ;  npz bound %.9f ;  C2 -> %s"
            % (Jint, sum_n2, sum_p2, kap2, bnd, "HOLDS" if c2 else "FAILS"))
        rec.update(sep_file=fname, sep_M=Mint, sep_J=Jint, sep_S=Sint, sep_sum_n2=sum_n2,
                   kappa2_W=kap2, npz_bound=bnd, C2=bool(c2), C3=bool(c3), C4=c4)
        if not (c2 and c3):
            say("          VOID: C2/C3 failed; this separator yields no interval")
            continue
        W = n.astype(np.float64) / DEN
        u_cert = W / float(np.linalg.norm(W))
        L_cert = kap2 * normp
        src = fname

    # ---------------- C. certified intervals
    if u_cert is None:
        say("      no certified separator at this cell: no interval")
        rec["void"] = rec.get("void", "no separator")
        store.append(rec)
        return
    Ucert = KAPPA_UP[e] * normp
    ratio = min(1.0, L_cert / Ucert)
    eps = math.sqrt(max(0.0, 2.0 * (1.0 - ratio)))
    vu = float(np.dot(v, u_cert))
    rec.update(sep_src=src, L_dist=L_cert, U_dist=Ucert, eps=eps, v_dot_u=vu,
               v_dot_u_lo=vu - normv * eps, v_dot_u_hi=vu + normv * eps,
               u1=float(u_cert[0]))
    say("      C certified interval: L=%.9f  U=%.9f  eps = sqrt(2(1-L/U)) = %.6f   (||v||_2=%.6f)"
        % (L_cert, Ucert, eps, normv))
    say("      u_{Q,1} = %+.6f ;  <v,u_Q> = %+.6f  in  [%+.6f, %+.6f]"
        % (u_cert[0], vu, vu - normv * eps, vu + normv * eps))
    hd = {}
    for D in DLIST:
        m = rows <= D
        s = float(np.linalg.norm(u_cert[m]))
        lo = max(0.0, s - eps)
        hi = min(1.0, s + eps)
        hd[D] = dict(sqrtH=s, H=s * s, H_lo=lo * lo, H_hi=hi * hi,
                     vhead=head_tail(v, rows, D)[0], vtail=head_tail(v, rows, D)[1],
                     bound_II=head_tail(v, rows, D)[0] * hi + head_tail(v, rows, D)[1])
        say("      D=%-5d sqrt(H_D)=%.6f  H_D=%.6f in [%.6f, %.6f] ;  ||v_<=D||=%.5f "
            "||v_>D||=%.5f ;  (II) bracket term = %.5f"
            % (D, s, s * s, lo * lo, hi * hi, hd[D]["vhead"], hd[D]["vtail"], hd[D]["bound_II"]))
    rec["HD"] = hd
    rec["H1_times_piQ"] = hd[1]["H"] * float(len(primes_upto(Q)))
    # the mechanism's numbers, where the exact dual is available
    if u_exact is not None:
        pr30 = np.asarray(primes_upto(min(30, Q)), dtype=np.int64)
        jp = np.searchsorted(rows, pr30)                   # rows is sorted; no Python loop
        jp = jp[(jp < rows.size)]
        jp = jp[rows[jp] == pr30[:jp.size]]
        up = u_cert[jp]
        rat = {int(rows[j]): float(u_cert[j] / (-u_cert[0] / 2.0)) for j in jp}
        rec["u_p_over_half"] = rat
        rec["u_p_positive"] = bool(np.all(up > 0))
        say("      mechanism: u_p / (-u_1/2) for p <= 30: %s"
            % "  ".join("%d:%.3f" % (p, r) for p, r in sorted(rat.items())))
        if rec.get("sp_nonbip_k"):
            k = rec["sp_nonbip_k"]
            say("      proved ceiling u_1^2 <= 4/k with k=%d: %.6f <= %.6f -> %s"
                % (k, u_cert[0] ** 2, 4.0 / k, "HOLDS" if u_cert[0] ** 2 <= 4.0 / k else "FAILS"))
            rec["ceiling_ok"] = bool(u_cert[0] ** 2 <= 4.0 / k)
    rec["seconds"] = time.time() - t0
    store.append(rec)
    say("      %.0f s" % rec["seconds"])


def main(argv):
    cells = tuple(int(a) for a in argv if a.isdigit()) or (16, 18, 20, 22, 24, 26, 28)
    nth = max(1, min(NTHREAD, 8, int(numba.config.NUMBA_NUM_THREADS)))
    set_num_threads(nth)
    say("conedual_separator_mass -- direction 1 (clauses in the file header)")
    say("cells %s   threads %d   DEN = 10^%d   tight tol %g" % (list(cells), nth, 12, TIGHT_TOL))
    store = []
    for e in cells:
        try:
            do_cell(e, store)
        except MemoryError:
            say("2^%d MemoryError; cell not attempted" % e)
        flush(dict(cells=store, complete=False))

    good = [r for r in store if "void" not in r and "eps" in r]
    tightc = [r for r in good if r["eps"] <= EPS_TIGHT]
    say("")
    say("=" * 100)
    v = {}
    v["C0"] = all(r.get("C0", True) for r in store)
    v["C1"] = all(r["C1"] for r in store if "C1" in r)
    v["C2"] = all(r["C2"] for r in store if "C2" in r)
    v["C3"] = all(r["C3"] for r in store if "C3" in r)
    v["C4"] = all(r["C4"] for r in store if r.get("C4") is not None)
    v["P1"] = all(r["u1"] < 0 for r in good) and all(r.get("u_p_positive", True) for r in good)
    v["P2"] = (all(r["v_dot_u"] > 0 for r in good)
               and all(r["v_dot_u_lo"] > 0 for r in tightc)
               and all(P2_BAND[0] <= r["v_dot_u"] <= P2_BAND[1] for r in good))
    dec = True
    for D in DLIST:
        hs = [r["HD"][D]["H"] for r in tightc]
        dec = dec and all(hs[i + 1] < hs[i] for i in range(len(hs) - 1))
    v["P3"] = dec and all(P3_BAND[0] <= r["H1_times_piQ"] <= P3_BAND[1] for r in tightc)
    ps = [r for r in good if "u_p_over_half" in r]
    v["P4"] = bool(ps) and all(P4_BAND[0] <= x <= P4_BAND[1] for r in ps for x in r["u_p_over_half"].values())
    v["P5"] = bool(ps) and all(r.get("sp_nonbip_k", 0) >= P5_FRAC * r.get("npi", 1) for r in ps) \
        and all(r.get("ceiling_ok", False) for r in ps)
    for k in ("C0", "C1", "C2", "C3", "C4", "P1", "P2", "P3", "P4", "P5"):
        say("%-3s %s" % (k, "HOLDS" if v[k] else "FAILS/REFUTED"))
    say("")
    say("VERDICT  " + " | ".join("%s %s" % (k, "HOLDS" if v[k] else "FAILS") for k in
                                 ("C0", "C1", "C2", "C3", "C4", "P1", "P2", "P3", "P4", "P5")))
    if not v["C1"]:
        say("VOID enforced: C1 failed -> the exact-dual cells carry no interval (see per-cell 'void').")
    flush(dict(cells=store, clauses=v, complete=True,
               cited=dict(KAPPA_LO=KAPPA_LO, KAPPA_UP=KAPPA_UP, S1E=S1E)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

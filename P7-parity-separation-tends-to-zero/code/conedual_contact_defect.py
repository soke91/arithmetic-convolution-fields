# -*- coding: utf-8 -*-
r"""conedual_contact_defect -- the averaged normalised contact defect.

Supporting computation; not used in the paper's proofs or tables.

WHERE THIS COMES FROM

An earlier prime-forcing base is NOT an algebraic consequence
of the contact identities: on a BIPARTITE contact component the equations `u_1 + u_p + u_d = 0`
admit two alternating values, and
`u = (-e_1 + e_251)/sqrt(2)` at `Q = 256` is a counterexample -- unit, polar-feasible, with a contact
of the kind that base assumes for every odd squarefree `a <= 256`, and yet `u_1 = -1/sqrt 2`,
`u_2 = 0`.  Hypothesis (H)
(EXACT contact for every odd-omega `a`) already fails at the measured cells -- the tight-semiprime graph has
an isolated prime -- while route (I) is NOT dead, because only the NEGATIVE part of `<v, u_Q>` has to
vanish.  The replacement implemented and measured here is an AVERAGED NORMALISED CONTACT DEFECT.

THE OBJECT.  With `u_Q` the optimal separator (`||u_Q||_2 = 1`, `<u_Q, t> = dist_2(t, K_Q) = kappa ||p||_2`),
put for every even-band column `l`

    delta_Q(l) := -theta_{u_Q}(l) = -sum_{d | l, d <= Q} u_{Q,d}   >= 0        (u_Q in K_Q^o)

-- the contact defect; it is already normalised, `||u_Q||_2 = 1`.  Three fibres, with
`P_Q := {q prime : Q/2 < q <= Q}`, `m_P := |P_Q|`:

    F0(Q)    = { p q      : p < q in P_Q }                                 divisors in [1,Q] = {1,p,q}
    F1(a,Q)  = { a q      : q in P_Q, q not | a }      omega(a) odd        divisors = {s : s|a} + {q}
    F2(a,Q)  = { a r q    : r prime <= Q/a, r not | a, q in P_Q }  omega(a) even   = {s : s|ar} + {q}

(The divisor identifications are those of the contact lemma: a divisor `s q` with `s > 1` is
`>= 2q > Q`.)  Averaging the exact identities over a fibre and using Cauchy-Schwarz on `||u_Q||_2 = 1`
gives the averaged-defect inequality and the averaged-defect theorem:

    |u_{Q,1}|        <= 2 m_P^{-1/2} + Delta0(Q),            Delta0 := mean_{F0} delta_Q
    |sum_{s|a} u_s|  <= m_P^{-1/2} + Delta1(a,Q)                              (omega(a) odd)
    |sum_{s|a} u_s|  <= m_P^{-1/2} + 2^{omega(a)} |R_a|^{-1/2} + Delta2(a,Q)   (omega(a) even)

and, by induction on omega over the divisors of a fixed `m`, `|u_{Q,m}| -> 0`; hence `H_D -> 0` and
`kappa -> 0`.  No contact is required to be EXACT, no graph and no odd cycle enter.

ONE UNCONDITIONAL INPUT, PROVED HERE.  `z = (1/R_1) sum_{l in vB} A_l` is in `K_Q`, so
`<u_Q, z> <= 0`, and `<u_Q,z> = <u_Q,t> + <u_Q,p>`; therefore

    (1/R_1) sum_{l in vB} delta_Q(l) = -( kappa ||p||_2 + <u_Q, p> ) <= ||p||_2 (1 - kappa) ,

i.e. the defect summed over the WHOLE even band is `O(N/(R_1 L^2))` -- the global average defect tends to
zero at rate `1/L^2` with no hypothesis at all.  What is NOT proved is that it does not concentrate off the
thin fibres above; `|F0| ~ m_P^2/2 ~ N/(8L^2)` against `|vB| ~ c N`, so the global bound only gives
`Delta0 = O(1)`.  That factor is the whole remaining gap and this file measures it.

WHAT IS MEASURED, per cell
  1. the separator, with a CERTIFIED dual lower bound -- a lower bound, not the primal NNLS
     residual, which is an UPPER bound.  At 2^16/2^18 the NNLS dual direction is integer-repaired exactly as the colgen
     separators are -- n = rint(u*DEN), n[0] := -M with M = max_{l in vB} theta_n(l) by an exact int64
     sweep -- so kappa_2(W) = J/sqrt(sum n^2 * sum (c_B-c_R)^2) is a certified lower bound from integers.
     At 2^22..2^28 the colgen npz with M CITED from conedual_s1_evidence (certified there).
  2. eps = sqrt(2(1 - L/U)) with U the tighter upper endpoint.
  3. Delta0, Delta1(a), Delta2(a) over the fibres; the per-omega-class averages; the proportion of exact
     contacts; the global defect and the proved ceiling on Delta0.
  4. the tight-semiprime graph: components, bipartite or not, and the ISOLATED primes.
  5. route (I): <v, u_Q> and its certified interval (the sign, not a margin, is what matters).

INPUTS CITED, NOT RECOMPUTED
  kappa enclosures, computed earlier: 2^16 [0.1723894, 0.1723894], 2^18 [0.1408600, 0.1408605],
    2^22 [0.0763310, 0.0763313], 2^24 [0.0527620, 0.0527630], 2^26 [0.0346220, 0.0346234],
    2^28 [0.0245260, 0.0245279].
  (M, J, S) of the four colgen separators and their K_Q^o certification: results/conedual_s1_evidence.txt.
  ||R_1 p||_2^2 at 2^28 = 2504581512105, computed earlier -- this machine has no 2^28 cell cache.
  The stored kappa_2(W) at 2^26 = 0.0346229227995222 (results/conedual_separator_mass.json), clause C5.

CLAUSES.

  C1  CONTROL (exact-dual cells).  The integer-repaired dual is in K_Q^o with no floating slack --
      max_{l in vB}(theta_n(l) - M) == 0 in exact int64 -- and the certified kappa_2(W) lies in the CITED
      enclosure widened by 1e-6.  CAN FAIL; voids that cell.
  C2  CONTROL (colgen cells).  (M, J, S) computed here equal conedual_s1_evidence's integers -- integer
      against integer -- and kappa_2(W) lies in [npz bound - 1e-7, U_cited + 1e-7].  Where this machine has
      no cell cache (2^28), J is CITED and that half is reported as CITED, not scored.
      CAN FAIL; that separator yields no defect numbers.
  C3  CONTROL.  The subset-zeta transform A_w(a) = sum_{m|a} w_m against a direct divisor sum on every row
      a <= 300, MIXED tolerance |dev| <= 1e-9 + 1e-12 |direct| (a signed sum passes through zero, so a
      purely relative tolerance would test the comparand's conditioning, not the transform).
      CAN FAIL; voids everything.
  C4  CONTROL (the global-defect identity, at the cells whose band is swept here).
      (1/R_1) sum_{l in vB} delta_Q(l) computed by the band sweep equals -(kappa_2(W) ||p||_2 + <W/||W||, p>)
      to a relative 1e-9.  CAN FAIL; voids the global-defect numbers at that cell.
  C5  CONTROL (reproduction, like with like: stored float against recomputed float, same formula).
      kappa_2(W) at 2^26 equals conedual_separator_mass.json's 0.0346229227995222 to a relative 1e-12.
      CAN FAIL; voids the 2^26 cell.

  P1  PREDICTION (scores the arithmetic of a PROVED inequality, the averaged-defect inequality).
      |u_{Q,1}| <= 2 m_P^{-1/2} + Delta0(Q) at every scored cell.
  P2  THE CLAUSE THAT DECIDES THE ROUTE.  Delta0(Q) <= 0.01 at every scored cell.
      PREDICTED: Delta0 <= 1e-3 at 2^16 and 2^18; <= 0.01 at the colgen cells.
      WHY: delta_Q(pq) = -(u_1 + u_p + u_q) vanishes identically when u_p = -u_1/2 on P_Q, which
      an earlier run measured to hold for every prime except one isolated one per cell, whose own
      clean-column slack is 1.6e-3 at 2^16.  A mean over ~|P_Q|^2/2 pairs of which O(|P_Q|) involve that prime is
      then O(1e-3/|P_Q|).  MEANING -- holds: the averaged-defect hypothesis is satisfied
      with large margin at every computed cell, so the route is alive and the open problem is to prove
      Delta0 -> 0 rather than to find a better hypothesis; fails: the fibre carries the defect and the
      averaged route needs a different fibre.
  P3  PREDICTION.  Delta1(a,Q) <= 0.05 for every odd-omega squarefree a <= 30 and Delta2(a,Q) <= 0.05 for
      every even-omega squarefree a <= 30, at every scored cell.
  P4  PREDICTION (route (I)).  <v, u_Q> > 0 at every scored cell AND its certified lower
      endpoint <v, u_cert> - ||v||_2 eps > 0 at every scored cell.
      MEANING -- holds: inequality (I) gives kappa <= r_B at every computed cell and route (I) is alive
      (an earlier "route (I) is dead" is withdrawn); fails: only the H_D route survives.
  P5  PREDICTION.  The number of primes p <= Q with no tight semiprime partner (|u_1+u_p+u_q| > 1e-9 for
      every admissible q) is at most 2 at every scored cell, and the tight-semiprime graph's largest
      component covers at least 0.9 pi(Q).
      MEANING -- holds: the isolated prime is one prime, not a growing set, so (H)'s failure is
      that narrow.

  D1  DIAGNOSTIC, needed to read P5's verdict: P5's absolute 1e-9 is calibrated to the exact duals'
      rounding floor ~1e-12 and is meaningless against a colgen separator whose own residual is
      ~1e-7, so P5's FAILS at 2^22..2^28 is a like-with-like artefact and NOT a structural
      statement.  D1 is not a prediction and scores nothing: it
      reports the same graph at a tolerance calibrated to the separator actually used,
      tol* := max(TOL_TIGHT, 10 * max_{F0} delta_Q), alongside the registered P5 numbers.
      P5's registered verdict stands as scored.

  NOT CLAIMED HERE: anything about (S1), sigma_min, Goldbach, or any asymptotic claim.  The
  averaged-defect inequality, the averaged-defect theorem and the global-defect bound are proved
  elsewhere; nothing here is evidence
  for them.

COMPUTE RULES (the project's compute rules, both sections read before this file was written)
  cells are LOADED (`from cell_cache import load`); the band is rebuilt only where an exact K_Q^o
  certification or the global-defect control needs it (2^16/2^18) and is checked against the cache
  entry for entry; no Python loop over rows, pairs or band elements -- numba prange kernels for the column
  matrix, the int64 theta sweep and the fibre scans, one vectorised pass per prime for the zeta transform,
  numpy broadcasting for F0 and the graph; the only Python loops are over cells, over the primes <= Q and
  over the handful of fibre bases a.  Threads capped (<= 8).  Cells 2^20 and above that need an N-length
  array want a machine with more cores; the colgen cells here need none.

    python code/conedual_contact_defect.py 16 18 22 24 26 28
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
RES = os.path.join(os.path.dirname(HERE), "results")
OUT = os.path.join(RES, "conedual_contact_defect.txt")
OUTJ = os.path.join(RES, "conedual_contact_defect.json")
REFC6 = os.path.join(RES, "conedual_separator_mass.json")

DEN = 10 ** 12
TOL_TIGHT = 1e-9
FIBRE_MAX = 30                 # the fibre bases a: every squarefree a <= 30
EXACT_CELLS = (16, 18)

KAP = {16: (0.1723894, 0.1723894), 18: (0.1408600, 0.1408605), 20: (0.1055427, 0.1055427),
       22: (0.0763310, 0.0763313), 24: (0.0527620, 0.0527630), 26: (0.0346220, 0.0346234),
       28: (0.0245260, 0.0245279)}
S1E = {"conedual_colgen_sep_2e22_n0b.npz": dict(M=1000000000001, J=50810068014373698, S=78055974358103479),
       "conedual_colgen_sep_2e24_n0b.npz": dict(M=1000001113767, J=173251720631015907, S=338483562472838081),
       "conedual_colgen_sep_2e26_n0b.npz": dict(M=1000000603706, J=463809741913557021, S=1088810738001557124),
       "conedual_colgen_sep_2e28_n28.npz": dict(M=1000001741817, J=1437321516767375728, S=3795904130585696957)}
SEPS = {22: "conedual_colgen_sep_2e22_n0b.npz", 24: "conedual_colgen_sep_2e24_n0b.npz",
        26: "conedual_colgen_sep_2e26_n0b.npz", 28: "conedual_colgen_sep_2e28_n28.npz"}
NORM_RP_SQ_CITED = {28: 2504581512105}
C6_KAPPA2_26 = 0.0346229227995222
C1_WIDEN, C2_SLACK, C3_ABS, C3_REL, C4_TOL, C5_TOL = 1e-6, 1e-7, 1e-9, 1e-12, 1e-9, 1e-12
P2_MAX, P3_MAX, P5_MAXISO, P5_MINFRAC = 0.01, 0.05, 2, 0.90

_lines = []


def say(s=""):
    _lines.append(s)
    print(s, flush=True)


def flush(payload):
    with io.open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(_lines) + "\n")
    with io.open(OUTJ, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=1, default=float)


def idot(a, b):
    return int(np.asarray(a, dtype=object) @ np.asarray(b, dtype=object))


def exact_sum(a, shift=24):
    a = np.asarray(a)
    assert a.dtype == np.int64
    if a.size == 0:
        return 0
    k = 1 << shift
    hi = int(np.abs(a).max()) // k + 1
    assert a.size * hi < 2 ** 62 and a.size * k < 2 ** 62
    return int(np.sum(a >> shift)) * k + int(np.sum(a & (k - 1)))


# ---------------------------------------------------------------- numba kernels

@njit(parallel=True, cache=True, nogil=True)
def _colmat(rows, pos, out):
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
    """acc[l] += sum_{d | l, d = rows[i], i >= 1} vals[i], exact int64, block-parallel."""
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


@njit(parallel=True, cache=True, nogil=True)
def _fibre2(Aar, rr, wq, prq, N, Q, tol, out):
    """F2 over (r, q): Aar[i] = A_w(a r_i) for the admissible r_i; wq over P_Q.  Returns per-r sums."""
    for i in prange(rr.shape[0]):
        s = 0.0
        c = 0
        nt = 0
        mx = 0.0
        for j in range(prq.shape[0]):
            q = prq[j]
            if rr[i] % q == 0:
                continue
            dl = -(Aar[i] + wq[j])
            if dl < 0.0:
                dl = 0.0
            s += dl
            c += 1
            if dl <= tol:
                nt += 1
            if dl > mx:
                mx = dl
        out[i, 0] = s
        out[i, 1] = c
        out[i, 2] = nt
        out[i, 3] = mx


# ---------------------------------------------------------------- helpers

def band(e, S):
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
    del rem
    vB = np.nonzero(bd & (om % 2 == 0))[0].astype(np.int64)
    vR = np.nonzero(bd & (om % 2 == 1))[0].astype(np.int64)
    del bd, om, mu
    return vB, vR


def row_arith(Q, rows):
    om = np.zeros(Q + 1, dtype=np.int64)
    phi = np.arange(Q + 1, dtype=np.int64)
    for p in primes_upto(Q):
        p = int(p)
        om[p::p] += 1
        phi[p::p] -= phi[p::p] // p
    return om[rows], phi[rows]


def zeta_transform(w, rows, Q):
    """A_w(a) = sum_{m|a} w_m: one vectorised pass per prime (Yates over the squarefree lattice)."""
    posr = -np.ones(Q + 2, dtype=np.int64)
    posr[rows] = np.arange(rows.size)
    A = w.astype(np.float64).copy()
    for p in primes_upto(Q):
        p = int(p)
        m = (rows % p == 0)
        if not m.any():
            continue
        A[np.nonzero(m)[0]] += A[posr[rows[m] // p]]
    return A, posr


def check_zeta(A, w, rows, posr, lim=300):
    worst_abs, worst_slack = 0.0, 0.0
    sel = np.nonzero(rows <= lim)[0]
    for i in sel:                                   # a handful of small rows
        a = int(rows[i])
        dv = rows[(a % np.maximum(rows, 1) == 0) & (rows <= a)]
        direct = float(w[posr[dv]].sum())
        err = abs(A[i] - direct)
        worst_abs = max(worst_abs, err)
        worst_slack = max(worst_slack, err / (C3_ABS + C3_REL * abs(direct)))
    return worst_abs, worst_slack, int(sel.size)


def tight_graph(wp, w1, prs, N, Q, tol):
    """Tight semiprime graph on the primes; returns components, bipartiteness, isolated primes."""
    s = w1 / 2.0 + wp
    P = prs.astype(np.float64)
    prod = P[:, None] * P[None, :]
    adm = (prod > Q) & (prod <= N) & np.triu(np.ones((prs.size, prs.size), dtype=bool), 1)
    tg = adm & (np.abs(s[:, None] + s[None, :]) <= tol)
    ei, ej = np.nonzero(tg)
    npr = int(prs.size)
    # a prime is "isolated" if it has an admissible partner but no tight one
    hasadm = (adm | adm.T).any(axis=1)
    hastg = (tg | tg.T).any(axis=1)
    iso = np.nonzero(hasadm & ~hastg)[0]
    out = dict(npi=npr, adm=int(adm.sum()), edges=int(ei.size), ncomp=0, kmax=0, nonbip_k=0,
               cls=(0, 0), n_iso=int(iso.size), iso=[int(x) for x in prs[iso][:20]])
    if ei.size == 0:
        return out
    Aa = coo_matrix((np.ones(ei.size), (ei, ej)), shape=(npr, npr))
    Aa = (Aa + Aa.T).tocsr()
    ncomp, lab = connected_components(Aa, directed=False)
    deg = np.diff(Aa.indptr)
    rl, ri = np.unique(lab, return_index=True)
    dist = shortest_path(Aa, method="D", unweighted=True, directed=False, indices=ri)
    slot = np.empty(ncomp, dtype=np.int64)
    slot[rl] = np.arange(rl.size)
    dv = dist[slot[lab], np.arange(npr)]
    col = np.where(np.isfinite(dv), np.mod(dv, 2.0), -1.0)
    bad = col[ei] == col[ej]
    nonbip = (np.bincount(lab[ei][bad], minlength=ncomp) > 0) if bad.any() else np.zeros(ncomp, bool)
    sizes = np.bincount(lab[deg > 0], minlength=ncomp)
    big = int(np.argmax(sizes))
    inbig = (lab == big) & (deg > 0)
    cand = np.nonzero(nonbip)[0]
    out.update(ncomp=int(ncomp), kmax=int(sizes.max()),
               nonbip_k=int(sizes[cand].max()) if cand.size else 0,
               cls=(int(np.sum(inbig & (col == 0))), int(np.sum(inbig & (col == 1)))))
    return out


# ---------------------------------------------------------------- one cell

def do_cell(e, store, ref6):
    t0 = time.time()
    N = 2 ** e
    Q = math.isqrt(N)
    rec = dict(e=e, N=N, Q=Q)
    have_cache = os.path.isdir(os.path.join(RES, "cell_cache_2e%d" % e))
    S = load_cell(e) if have_cache else None
    if have_cache:
        rows = np.asarray(S["rows"], dtype=np.int64)
        cb = np.asarray(S["cb"], dtype=np.int64)
        cr = np.asarray(S["cr"], dtype=np.int64)
        nR, nB = int(S["nR"]), int(S["nB"])
        Rp = cb - cr
        sum_p2 = idot(Rp, Rp)
    else:
        z0 = np.load(os.path.join(RES, SEPS[e]))
        rows = np.asarray(z0["rows"], dtype=np.int64)
        cb = cr = None
        nR = nB = None
        sum_p2 = NORM_RP_SQ_CITED[e]
        rec["norm_Rp_sq_source"] = "carried from an earlier run (no local cell cache)"
    om, phi = row_arith(Q, rows)
    prs = np.asarray(primes_upto(Q), dtype=np.int64)
    PQ = prs[prs > Q // 2]
    mP = int(PQ.size)
    say("")
    say("=" * 106)
    say("2^%-2d Q=%-6d rows=%-5d pi(Q)=%-5d |P_Q|=%-5d  kappa CITED [%.7f, %.7f]"
        % (e, Q, rows.size, prs.size, mP, KAP[e][0], KAP[e][1]))

    vB = None
    sweep_here = e in EXACT_CELLS
    if sweep_here:
        vB, vR = band(e, S)
        ok = (vB.size == nB) and (vR.size == nR)
        say("      band rebuilt: |vB| %d vs nB %d, |vR| %d vs nR %d -> %s"
            % (vB.size, nB, vR.size, nR, "EQUAL" if ok else "DIFFER"))
        if not ok:
            rec["void"] = "band"
            store.append(rec)
            return

    # ---------------- the separator, with a CERTIFIED dual lower bound
    if e in EXACT_CELLS:
        pos = -np.ones(N + 1, dtype=np.int64)
        pos[vB] = np.arange(vB.size)
        Dm = np.zeros((vB.size, rows.size), dtype=np.int8)
        _colmat(rows, pos, Dm)
        Uc = np.unique(Dm, axis=0)
        del Dm
        Mmat = np.ascontiguousarray(Uc.T.astype(np.float64))
        del Uc
        tvec = cr.astype(np.float64) / nR
        x, _ = nnls(Mmat, tvec, maxiter=20 * rows.size + 200)
        r = tvec - Mmat @ x
        del Mmat
        u0 = r / float(np.linalg.norm(r))
        n = np.rint(u0 * DEN).astype(np.int64)
        n[0] = 0
        acc = np.zeros(N + 1, dtype=np.int64)
        _theta_blocks(acc, rows, n, 8 * NTHREAD)
        WB = acc[vB]
        Mint = int(WB.max())
        c1a = int((WB - Mint).max()) == 0
        del acc, WB
        n[0] = -Mint
        Jint = idot(n[1:], cr[1:]) - Mint * nR
        Sint = Mint * nB - idot(n[1:], cb[1:])
        sum_n2 = idot(n, n)
        kap2 = Jint / math.sqrt(float(sum_n2) * float(sum_p2))
        lo, up = KAP[e]
        c1 = bool(c1a) and (lo - C1_WIDEN <= kap2 <= up + C1_WIDEN) and Jint > 0
        say("      exact dual, integer-repaired: M = %d ; C1 max(theta_n - M) == 0 -> %s" % (Mint, c1a))
        say("      C1 certified kappa_2(W) = %d/sqrt(%d * %d) = %.9f in CITED [%.7f, %.7f] -> %s"
            % (Jint, sum_n2, sum_p2, kap2, lo, up, "HOLDS" if c1 else "FAILS"))
        rec.update(src="exact dual + integer repair", M=Mint, J=Jint, S=Sint, sum_n2=sum_n2,
                   kappa2_W=kap2, C1=bool(c1))
        if not c1:
            rec["void"] = "C1"
            store.append(rec)
            return
    else:
        fn = SEPS[e]
        z = np.load(os.path.join(RES, fn))
        zr = np.asarray(z["rows"], dtype=np.int64)
        if zr.size != rows.size or not bool((zr == rows).all()):
            rec["void"] = "rows"
            say("      separator %s VOID: row order differs" % fn)
            store.append(rec)
            return
        n = np.rint(np.asarray(z["w"], dtype=np.float64) * DEN).astype(np.int64)
        cit = S1E[fn]
        Mint = cit["M"]
        n[0] = -Mint
        if cr is not None:
            Jint = idot(n[1:], cr[1:]) - Mint * nR
            Sint = Mint * nB - idot(n[1:], cb[1:])
            c2i = (Jint == cit["J"]) and (Sint == cit["S"])
        else:
            Jint, Sint, c2i = cit["J"], cit["S"], None
        sum_n2 = idot(n, n)
        kap2 = Jint / math.sqrt(float(sum_n2) * float(sum_p2))
        bnd = float(z["bound"])
        up = KAP[e][1]
        c2 = (kap2 >= bnd - C2_SLACK) and (kap2 <= up + C2_SLACK) and (c2i is not False)
        say("      separator %-40s M = %d CITED (K_Q^o certified in conedual_s1_evidence)" % (fn, Mint))
        say("      C2 integers J %s S %s ; kappa_2(W) = %.9f vs npz %.9f -> %s"
            % ("EQUAL" if c2i else ("CITED" if c2i is None else "DIFFER"),
               "EQUAL" if c2i else ("CITED" if c2i is None else "DIFFER"), kap2, bnd,
               "HOLDS" if c2 else "FAILS"))
        rec.update(src=fn, M=Mint, J=Jint, S=Sint, sum_n2=sum_n2, kappa2_W=kap2, npz_bound=bnd,
                   C2=bool(c2), C2_int=c2i)
        if not c2:
            rec["void"] = "C2"
            store.append(rec)
            return
        if e == 26:
            d5 = abs(kap2 - C6_KAPPA2_26) / C6_KAPPA2_26
            rec["C5"] = bool(d5 <= C5_TOL)
            say("      C5 vs the stored kappa_2(W) %.16f: rel dev %.2e -> %s"
                % (C6_KAPPA2_26, d5, "HOLDS" if rec["C5"] else "FAILS"))
            if not rec["C5"]:
                rec["void"] = "C5"
                store.append(rec)
                return

    W = n.astype(np.float64) / DEN
    w = W / float(np.linalg.norm(W))
    normp = math.sqrt(float(sum_p2)) / nR if nR else None
    L = rec["kappa2_W"] * (normp if normp else math.sqrt(float(sum_p2)))
    Ucert = KAP[e][1] * (normp if normp else math.sqrt(float(sum_p2)))
    eps = math.sqrt(max(0.0, 2.0 * (1.0 - min(1.0, L / Ucert))))
    rec.update(normp=normp, L=L, U=Ucert, eps=eps, u1=float(w[0]), mP=mP)
    say("      L=%.9f  U=%.9f  eps = sqrt(2(1-L/U)) = %.6f ;  u_1 = %+.8f ;  |P_Q| = %d"
        % (L, Ucert, eps, w[0], mP))

    # ---------------- A_w and its control
    A, posr = zeta_transform(w, rows, Q)
    wabs, wsl, nz = check_zeta(A, w, rows, posr)
    rec.update(C3=bool(wsl <= 1.0), C3_worst_abs=wabs)
    say("      C3 zeta transform vs direct divisor sums on %d rows: worst |dev| %.2e, slack %.3f -> %s"
        % (nz, wabs, wsl, "HOLDS" if wsl <= 1.0 else "FAILS"))
    if wsl > 1.0:
        rec["void"] = "C3"
        store.append(rec)
        return

    wq = np.ascontiguousarray(w[posr[PQ]])
    v = -np.asarray(mu_upto(Q), dtype=np.int8)[rows].astype(np.float64) / phi.astype(np.float64)
    normv = float(np.linalg.norm(v))
    vu = float(np.dot(v, w))
    rec.update(normv=normv, v_dot_u=vu, v_lo=vu - normv * eps, v_hi=vu + normv * eps)
    say("      route (I): <v,u_cert> = %+.8f ;  certified interval [%+.8f, %+.8f] ;  (-<v,u_Q>)_+ = %s"
        % (vu, vu - normv * eps, vu + normv * eps,
           "0 (certified)" if vu - normv * eps > 0 else "not certified 0"))

    # ---------------- F0 and the base inequality
    d0 = -(w[0] + wq[:, None] + wq[None, :])
    iu = np.triu_indices(mP, 1)
    d0v = np.maximum(d0[iu], 0.0)
    Delta0 = float(d0v.mean()) if d0v.size else float("nan")
    base = 2.0 / math.sqrt(mP) + Delta0
    rec.update(F0_size=int(d0v.size), Delta0=Delta0, F0_max=float(d0v.max()) if d0v.size else 0.0,
               F0_tightfrac=float(np.mean(d0v <= TOL_TIGHT)) if d0v.size else 0.0,
               base_bound=base, P1_ok=bool(abs(w[0]) <= base + 1e-12))
    say("      F0: %d pairs, Delta0 = mean delta = %.3e, max %.3e, exact-contact fraction %.4f"
        % (d0v.size, Delta0, rec["F0_max"], rec["F0_tightfrac"]))
    say("      P1 (PROVED) |u_1| = %.8f <= 2/sqrt(|P_Q|) + Delta0 = %.8f + %.3e = %.8f -> %s"
        % (abs(w[0]), 2.0 / math.sqrt(mP), Delta0, base, "HOLDS" if rec["P1_ok"] else "FAILS"))
    mq = float(wq.mean())
    rec["mean_u_q"] = mq
    rec["base_identity_resid"] = float(w[0] + 2.0 * mq + Delta0)
    say("      base identity u_1 + 2 mean_{P_Q} u_q = -Delta0 : %+.8f + 2(%+.8f) = %+.3e vs -%.3e"
        % (w[0], mq, w[0] + 2.0 * mq, Delta0))

    # ---------------- F1 / F2 over the fibre bases a <= FIBRE_MAX
    fib = {}
    bases = rows[rows <= FIBRE_MAX]
    for a in bases:                                  # a handful of fibre bases
        a = int(a)
        ia = int(posr[a])
        oa = int(om[ia])
        if oa % 2 == 1:                              # F1: a*q, q in P_Q
            qq = PQ[a % PQ != 0]
            dl = np.maximum(-(A[ia] + w[posr[qq]]), 0.0)
            fib[a] = dict(kind="F1", omega=oa, n=int(dl.size), mean=float(dl.mean()),
                          max=float(dl.max()), tightfrac=float(np.mean(dl <= TOL_TIGHT)),
                          sum_div=float(A[ia]))
        else:                                        # F2: a*r*q, r prime <= Q/a, q in P_Q
            rr = prs[(prs <= Q // a) & (a % np.maximum(prs, 1) != 0)]
            if rr.size == 0:
                continue
            ar = a * rr
            Aar = np.ascontiguousarray(A[posr[ar]])
            out = np.zeros((rr.size, 4))
            _fibre2(Aar, np.ascontiguousarray(rr), wq, np.ascontiguousarray(PQ), N, Q, TOL_TIGHT, out)
            tot, cnt = float(out[:, 0].sum()), float(out[:, 1].sum())
            fib[a] = dict(kind="F2", omega=oa, n=int(cnt), mean=(tot / cnt) if cnt else float("nan"),
                          max=float(out[:, 3].max()), tightfrac=float(out[:, 2].sum() / cnt) if cnt else 0.0,
                          sum_div=float(A[ia]), n_r=int(rr.size))
    rec["fibres"] = fib
    worst1 = max((f["mean"] for f in fib.values() if f["kind"] == "F1"), default=0.0)
    worst2 = max((f["mean"] for f in fib.values() if f["kind"] == "F2"), default=0.0)
    rec.update(Delta1_worst=worst1, Delta2_worst=worst2)
    say("      fibres a <= %d: %d bases ; worst Delta1 (odd omega) = %.3e ; worst Delta2 (even) = %.3e"
        % (FIBRE_MAX, len(fib), worst1, worst2))
    for a in sorted(fib)[:12]:
        f = fib[a]
        say("          a=%-4d %s omega=%d  n=%-7d mean delta=%.3e  max=%.3e  exact frac=%.4f  "
            "|sum_{s|a}u_s|=%.3e" % (a, f["kind"], f["omega"], f["n"], f["mean"], f["max"],
                                     f["tightfrac"], abs(f["sum_div"])))

    # ---------------- the global defect: identity + proved ceiling on Delta0
    gl = -(rec["kappa2_W"] * normp + float(np.dot(w, Rp.astype(np.float64) / nR))) if nR else None
    if gl is not None:
        ceil0 = (nR * gl) / max(1, rec["F0_size"])
        rec.update(global_mean_defect_rowside=gl,
                   global_bound=normp * (1.0 - rec["kappa2_W"]),
                   Delta0_proved_ceiling=ceil0)
        say("      global: (1/R_1) sum_vB delta = -(kappa||p|| + <u,p>) = %.6e  <=  ||p||(1-kappa) = %.6e"
            % (gl, normp * (1.0 - rec["kappa2_W"])))
        say("              => proved ceiling Delta0 <= R_1 * that / |F0| = %.4f  (measured %.3e)"
            % (ceil0, Delta0))
    if sweep_here:
        n2 = np.rint(w * DEN).astype(np.int64)
        acc = np.zeros(N + 1, dtype=np.int64)
        _theta_blocks(acc, rows, n2, 8 * NTHREAD)
        tot = exact_sum(acc[vB].astype(np.int64))
        swept = -(tot / DEN + float(n2[0]) / DEN * vB.size) / nR
        d4 = abs(swept - gl) / max(abs(gl), 1e-300)
        rec.update(C4=bool(d4 <= C4_TOL), global_mean_defect_swept=swept)
        say("      C4 global defect by band sweep = %.6e vs row identity %.6e : rel dev %.2e -> %s"
            % (swept, gl, d4, "HOLDS" if rec["C4"] else "FAILS"))
        del acc

    # ---------------- the tight graph and the isolated primes
    g = tight_graph(np.ascontiguousarray(w[posr[prs]]), float(w[0]), prs, N, Q, TOL_TIGHT)
    rec["graph"] = g
    say("      graph (P5, tol 1e-9): %d tight edges of %d admissible, %d components, largest %d of "
        "pi(Q)=%d, largest NON-BIPARTITE %d, classes (%d,%d) ; isolated primes %d %s"
        % (g["edges"], g["adm"], g["ncomp"], g["kmax"], g["npi"], g["nonbip_k"], g["cls"][0],
           g["cls"][1], g["n_iso"], g["iso"]))
    tol2 = max(TOL_TIGHT, 10.0 * rec["F0_max"])
    if tol2 > TOL_TIGHT:
        g2 = tight_graph(np.ascontiguousarray(w[posr[prs]]), float(w[0]), prs, N, Q, tol2)
        rec["graph_adaptive"] = dict(tol=tol2, **g2)
        say("      graph (D1, tol* = %.1e): %d tight edges of %d admissible, %d components, largest %d of "
            "pi(Q)=%d, largest NON-BIPARTITE %d, classes (%d,%d) ; isolated primes %d %s"
            % (tol2, g2["edges"], g2["adm"], g2["ncomp"], g2["kmax"], g2["npi"], g2["nonbip_k"],
               g2["cls"][0], g2["cls"][1], g2["n_iso"], g2["iso"]))
    for key in ("graph", "graph_adaptive"):
        gg = rec.get(key)
        if gg and gg["nonbip_k"] > 1:
            sp = np.abs(float(w[0]) / 2.0 + w[posr[prs]])
            say("      odd-cycle base (%s): non-bipartite component k = %d -> |u_1| <= 2/sqrt(k) + 2 max|s_p|"
                " = %.6f + %.3e ; measured |u_1| = %.8f"
                % (key, gg["nonbip_k"], 2.0 / math.sqrt(gg["nonbip_k"]), 2.0 * float(sp.max()), abs(w[0])))
            rec[key + "_oddcycle_base"] = 2.0 / math.sqrt(gg["nonbip_k"])
    rec["seconds"] = time.time() - t0
    store.append(rec)
    say("      %.0f s" % rec["seconds"])


def main(argv):
    cells = tuple(int(a) for a in argv if a.isdigit()) or (16, 18, 22, 24, 26, 28)
    nth = max(1, min(NTHREAD, 8, int(numba.config.NUMBA_NUM_THREADS)))
    set_num_threads(nth)
    say("conedual_contact_defect -- the averaged normalised contact defect (clauses in the header)")
    say("cells %s   threads %d   TOL_TIGHT %g   fibre bases a <= %d" % (list(cells), nth, TOL_TIGHT, FIBRE_MAX))
    ref6 = {}
    if os.path.exists(REFC6):
        for c in json.load(io.open(REFC6, encoding="utf-8"))["cells"]:
            ref6[c["e"]] = c
    store = []
    for e in cells:
        try:
            do_cell(e, store, ref6)
        except MemoryError:
            say("2^%d MemoryError; cell not attempted" % e)
        flush(dict(cells=store, complete=False))

    good = [r for r in store if "void" not in r and "Delta0" in r]
    v = {}
    v["C1"] = all(r["C1"] for r in store if "C1" in r)
    v["C2"] = all(r["C2"] for r in store if "C2" in r)
    v["C3"] = all(r["C3"] for r in store if "C3" in r)
    v["C4"] = all(r["C4"] for r in store if "C4" in r)
    v["C5"] = all(r["C5"] for r in store if "C5" in r)
    v["P1"] = bool(good) and all(r["P1_ok"] for r in good)
    v["P2"] = bool(good) and all(r["Delta0"] <= P2_MAX for r in good)
    v["P3"] = bool(good) and all(max(r["Delta1_worst"], r["Delta2_worst"]) <= P3_MAX for r in good)
    v["P4"] = bool(good) and all(r["v_dot_u"] > 0 and r["v_lo"] > 0 for r in good)
    v["P5"] = bool(good) and all(r["graph"]["n_iso"] <= P5_MAXISO
                                 and r["graph"]["kmax"] >= P5_MINFRAC * r["graph"]["npi"] for r in good)
    say("")
    say("=" * 106)
    for k in ("C1", "C2", "C3", "C4", "C5", "P1", "P2", "P3", "P4", "P5"):
        say("%-3s %s" % (k, "HOLDS" if v[k] else "FAILS/REFUTED"))
    say("")
    say("VERDICT  " + " | ".join("%s %s" % (k, "HOLDS" if v[k] else "FAILS")
                                 for k in ("C1", "C2", "C3", "C4", "C5", "P1", "P2", "P3", "P4", "P5")))
    voids = [(r["e"], r["void"]) for r in store if "void" in r]
    if voids:
        say("VOID enforced at: " + ", ".join("2^%d (%s)" % z for z in voids))
    flush(dict(cells=store, clauses=v, complete=True,
               cited=dict(KAP=KAP, S1E=S1E, NORM_RP_SQ=NORM_RP_SQ_CITED, C6_KAPPA2_26=C6_KAPPA2_26)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

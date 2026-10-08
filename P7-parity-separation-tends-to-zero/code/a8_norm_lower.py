# -*- coding: utf-8 -*-
r"""a8_norm_lower -- corroboration for a second, independent derivation of (NORM).

Supporting computation; not used in the paper's proofs or tables.

WHERE THIS COMES FROM

The proof that kappa(Q) -> 0 imports a normalisation lower bound, used only as

  (NORM)   ||p||_2 >= |p_1| >= (1/5)(N/R)L^{-2}   for L >= L_1,     L = log Q,

and the argument in fact needs only ||p||_2 >> L^{-A} for SOME fixed A.  One route proves (NORM)
through a per-row Selberg-Delange asymptotic, uniform in d <= L^10, with the coprime moments
sigma_j(d) and the V = e^{L^{1/4}} split.  The derivation checked here re-derives it at d = 1
only, from the one-large-prime identity and three classical inputs, with none of that machinery.
This file corroborates the two identities that derivation rests on and the size hierarchy it
asserts; the derivation itself is not reproduced here.

THE TWO IDENTITIES (both exact)

  (1)  R p_1 = M(N) - M(thr) + sum_{k <= thr} mu(k) [ pi(N/k) - pi(Q) ]
             = M(N) - M(thr)(1 + pi(Q)) + sum_{k <= thr} mu(k) pi(N/k)        -- integers throughout,
  (2)  T := sum_{k <= Y} mu(k) li(N/k)
         = li(2) M(Y) + M(Y) Li(N/Y) + N * int_1^Y M(floor v) dv / (v^2 log(N/v)),   Y = thr,

with Q = floor(sqrt N), q_1 the least prime > Q, thr = floor(N/q_1), Li(y) = int_2^y dt/log t, li the
principal-value logarithmic integral, M the Mertens function and R = |vR| the band's odd-omega count.
(2) is what replaces the li-expansion-plus-coprime-moments route: the main term -N/(log N)^2 of R p_1 comes out
of it from sum_n mu(n)/n = 0 and sum_n mu(n) log n / n = -1 alone.

CLAUSES.

  C1  CONTROL.  Identity (1) holds with deviation EXACTLY 0, as integers, at every cell 2^16 ... 2^24.
  C2  CONTROL.  Identity (2) reproduces T = sum_{k<=Y} mu(k) li(N/k) to <= 1e-11 relative at every cell.
  A FAILURE OF C1 OR C2 VOIDS R1 AND P1-P4, and this file prints them as VOID rather than as verdicts.

  R1  REPRODUCTION.  R p_1 / (-N/(log N)^2), rounded to four decimals, equals the previously
      published p_1/pred at the five cells: 1.1580, 1.2286, 1.2104, 1.2332, 1.2192 (2^16, 2^18, 2^20, 2^22, 2^24).
      Rounded display against rounded display, those literals carrying four decimals -- so the comparison is
      like with like (the reproduction rule).  At even e, Q = 2^{e/2} exactly, so log N = 2L exactly and
      the -N/(4RL^2) of that comparison and this file's -N/(R (log N)^2) are the same number.
  P1  |J_0| = |int_1^Y M(floor v) dv / v^2| < 1/log N at every cell.  (The step the proof needs,
      "m_0(Y) = sum_{k<=Y} mu(k)/k must be o(1/log N)", is already true at computable size.)
  P2  |M(N)| + |M(thr)|(1 + pi(Q)) < 0.5 * N/(log N)^2 at every cell.  (The two elementary error terms of
      identity (1) are already below the main term.)
  P3a |sum_{k<=thr} mu(k) delta(N/k)| < 0.05 * N/(log N)^2 at every cell, delta = pi - li.
  P3b sum_{k<=thr} |delta(N/k)| > N/(log N)^2 at 2^24.  (The proof bounds the PNT error by the triangle
      inequality over k; P3a + P3b say that the signed error is already negligible while the proof's bound
      on it is not -- i.e. the derivation's L_0 is above every computable cell, exactly as the
      published L_1 is.)
  P4  R p_1 < 0 at every cell, and R p_1 / (-N/(log N)^2) lies in (1.10, 1.26) at every cell.

  DISCLOSED, one run seen before this file existed: a scratch script (not this file) at 2^16, 2^18, 2^20,
  2^22 printed R p_1/(-N/(log N)^2) = 1.1580, 1.2286, 1.2104, 1.2332; |M(N)|/main = 0.0263, 0.0143, 0.0471,
  0.0126; |M(thr)|(1+pi(Q))/main = 0.1032, 0.2910, 0.0317, 0.0687; signed-delta/main = 0.0318, 0.0373,
  -0.0197, -0.0047; triangle sum|delta|/main = 4.18, 3.28, 2.53, 1.95; J_0 = 0.008456, 0.003400, 0.002422,
  0.000851 against 1/log N = 0.090168, 0.080150, 0.072135, 0.065577.  It also carried two bugs, both found
  by the controls and both fixed here: it subtracted pi(Q) M(thr) twice (identity (1) deviated by 54, 485,
  172, -1236 instead of 0) and it integrated (2) over [1, Y+1) instead of [1, Y] (relative 2.6e-4 instead
  of 1e-15).  So P1, P2, P3a and P4 are informed by four cells; 2^24 is new, P3b is new, and C1/C2 as
  stated here were never seen to pass.

    python code/a8_norm_lower.py
"""
from __future__ import annotations
import io, json, math, os, sys, time

NTHREAD = int(os.environ.get("ACF_THREADS") or os.environ.get("NUMBA_NUM_THREADS") or "4")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, str(NTHREAD))

import numpy as np  # noqa: E402
from scipy.special import expi  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _sieve_shared import mu_upto, primes_upto  # noqa: E402  (the canonical sieve)
from cell_cache import load as load_cell  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "a8_norm_lower.txt")
OUTJ = os.path.join(RES, "a8_norm_lower.json")

CELLS = (16, 18, 20, 22, 24)
A3_RATIOS = {16: 1.1580, 18: 1.2286, 20: 1.2104, 22: 1.2332, 24: 1.2192}   # the published comparison values
C2_TOL = 1.0e-11
P2_FACTOR = 0.5
P3A_FACTOR = 0.05
P4_LO, P4_HI = 1.10, 1.26
GAMMA = 0.5772156649015329
A3_K4 = 11.02547          # a_3(1), for the second-order comparison only

_f = io.open(OUT, "w", encoding="utf-8", newline="\n")


def say(s=""):
    print(s, flush=True)
    _f.write(s + "\n")
    _f.flush()


def li(y):
    """Principal-value li(y) = Ei(log y)."""
    return expi(np.log(np.asarray(y, dtype=np.float64)))


def nextprime(n):
    q = n + 1
    while True:
        if q == 2:
            return 2
        if q > 2 and q % 2 and all(q % r for r in range(3, int(q ** 0.5) + 1, 2)):
            return q
        q += 1


_XG, _WG = np.polynomial.legendre.leggauss(16)


def step_integral(Mj, Y, N, power):
    """int_1^Y M(floor v) (log v)^power dv / (v^2 log(N/v)) if power >= 0, or .../v^2 if power < 0.

    M is constant on [j, j+1), so the integral is a sum over j = 1 .. Y-1 of one smooth piece each,
    evaluated by 16-point Gauss-Legendre.  Vectorised over j; no Python loop.
    """
    j = np.arange(1, Y, dtype=np.float64)
    t = 0.5 * _XG[None, :] + (j[:, None] + 0.5)            # nodes on [j, j+1]
    w = 0.5 * _WG[None, :]
    if power < 0:
        f = 1.0 / t ** 2
    else:
        f = np.log(t) ** power / (t ** 2 * np.log(N / t))
    return float(np.dot(Mj[: Y - 1], (w * f).sum(axis=1)))


def main(argv):
    cells = [int(a) for a in argv if a.isdigit()] or list(CELLS)
    say("a8_norm_lower -- (NORM) re-derived at d = 1 (clauses in the header)")
    say("threads %d; this machine, seconds (every array below 2^24 and no cell is rebuilt: cell_cache.load)"
        % NTHREAD)
    out = []
    for e in cells:
        t0 = time.time()
        N = 2 ** e
        Q = int(N ** 0.5)
        assert Q == math.isqrt(N)
        q1 = nextprime(Q)
        thr = N // q1
        Y = thr
        H = math.log(N)
        L = math.log(Q)
        mu = np.array(mu_upto(N), dtype=np.int8)
        Msig = np.cumsum(mu.astype(np.int32))                      # M(n), exact
        pr = np.asarray(primes_upto(N), dtype=np.int64)
        pit = np.zeros(N + 1, dtype=np.int32)
        pit[pr] = 1
        np.cumsum(pit, out=pit)
        del pr
        S = load_cell(e)
        assert int(S["rows"][0]) == 1
        Rp1 = int(S["cb"][0]) - int(S["cr"][0])
        R = int(S["nR"])
        del S

        k = np.arange(1, Y + 1, dtype=np.int64)
        muk = mu[k].astype(np.int64)
        xk = N // k
        piQ = int(pit[Q])
        MN, Mthr = int(Msig[N]), int(Msig[thr])
        # ---- C1: identity (1), as integers
        Spi = int(np.dot(muk, (pit[xk].astype(np.int64) - piQ)))
        c1_dev = MN - Mthr + Spi - Rp1
        # ---- C2: identity (2)
        yk = N / k.astype(np.float64)
        liyk = li(yk)
        T = float(np.dot(muk.astype(np.float64), liyk))
        MY = int(Msig[Y])
        Mj = Msig[1:Y + 1].astype(np.float64)                      # M(j) on [j, j+1)
        I1 = step_integral(Mj, Y, N, 1)                            # with the log(N/v) weight, power 1 unused
        Iw = step_integral(Mj, Y, N, 0)
        Tid = float(li(2.0)) * MY + MY * float(li(N / Y) - li(2.0)) + N * Iw
        c2_rel = abs(T - Tid) / abs(T)
        # ---- the pieces of the proof
        J0 = step_integral(Mj, Y, N, -1)                            # int_1^Y M dv/v^2
        m0 = float(np.sum(muk / k.astype(np.float64)))
        m1 = float(np.sum(muk * np.log(k.astype(np.float64)) / k.astype(np.float64)))
        main_t = -N / H ** 2
        dsg = pit[xk].astype(np.float64) - liyk
        d_signed = float(np.dot(muk.astype(np.float64), dsg))
        d_tri = float(np.abs(dsg).sum())
        elem = abs(MN) + abs(Mthr) * (1 + piQ)
        ratio = Rp1 / main_t
        pred = 1.0 + (2.0 + 2.0 * GAMMA) / H + A3_K4 / H ** 2
        rec = dict(e=e, N=N, Q=Q, q1=q1, thr=thr, H=H, L=L, R=R, Rp1=Rp1, M_N=MN, M_thr=Mthr, piQ=piQ,
                   main=main_t, ratio=ratio, pred=pred, c1_dev=int(c1_dev), c2_rel=c2_rel,
                   T=T, T_identity=Tid, J0=J0, m0=m0, m1=m1, inv_H=1.0 / H,
                   elem=elem, elem_over_main=elem / abs(main_t),
                   d_signed=d_signed, d_signed_over_main=d_signed / abs(main_t),
                   d_tri=d_tri, d_tri_over_main=d_tri / abs(main_t),
                   p1=Rp1 / R, norm_floor=N / (9.0 * R * L ** 2), seconds=time.time() - t0)
        out.append(rec)
        del mu, Msig, pit
        say("2^%d Q=%d q1=%d thr=%d  R p_1=%d  C1 dev=%d  C2 rel=%.2e   R p_1/main=%.4f (ref %.4f, 2nd-order"
            " pred %.4f)" % (e, Q, q1, thr, Rp1, c1_dev, c2_rel, ratio, A3_RATIOS.get(e, float("nan")), pred))
        say("      J_0=%.6f vs 1/log N=%.6f | m_0(Y)=%.6f m_1(Y)=%.6f | elem/main=%.4f |"
            " signed delta/main=%+.4f | triangle/main=%.3f   %.0f s"
            % (J0, 1.0 / H, m0, m1, rec["elem_over_main"], rec["d_signed_over_main"],
               rec["d_tri_over_main"], rec["seconds"]))

    C1 = all(r["c1_dev"] == 0 for r in out)
    C2 = all(r["c2_rel"] <= C2_TOL for r in out)
    say("")
    say("C1  identity (1) deviation, as integers, at every cell: %s   (%s)"
        % ("HOLDS" if C1 else "FAILS", ", ".join(str(r["c1_dev"]) for r in out)))
    say("C2  identity (2) relative deviation: %s   (%s)   tol %.0e"
        % ("HOLDS" if C2 else "FAILS", ", ".join("%.1e" % r["c2_rel"] for r in out), C2_TOL))
    say("")
    if not (C1 and C2):
        say("CONTROL FAILED: R1 and P1-P4 are VOID (the consequence the clause states, enforced here, not waived).")
        R1 = P1 = P2 = P3a = P3b = None
    else:
        v = lambda b: "HOLDS" if b else "REFUTED"
        R1 = all(round(r["ratio"], 4) == A3_RATIOS[r["e"]] for r in out if r["e"] in A3_RATIOS)
        P1 = all(abs(r["J0"]) < r["inv_H"] for r in out)
        P2 = all(r["elem"] < P2_FACTOR * abs(r["main"]) for r in out)
        P3a = all(abs(r["d_signed"]) < P3A_FACTOR * abs(r["main"]) for r in out)
        c24 = [r for r in out if r["e"] == 24]
        P3b = bool(c24) and c24[0]["d_tri"] > abs(c24[0]["main"])
        P4 = all(r["Rp1"] < 0 and P4_LO < r["ratio"] < P4_HI for r in out)
        say("R1  R p_1/main rounded to 4 dp equals the published comparison values: %s"
            % v(R1))
        say("    mine %s" % ", ".join("%.4f" % round(r["ratio"], 4) for r in out if r["e"] in A3_RATIOS))
        say("    ref  %s" % ", ".join("%.4f" % A3_RATIOS[r["e"]] for r in out if r["e"] in A3_RATIOS))
        say("P1  |J_0| < 1/log N at every cell: %s   (ratios %s)"
            % (v(P1), ", ".join("%.3f" % (abs(r["J0"]) * r["H"]) for r in out)))
        say("P2  |M(N)|+|M(thr)|(1+pi(Q)) < %.2f main at every cell: %s   (%s)"
            % (P2_FACTOR, v(P2), ", ".join("%.4f" % r["elem_over_main"] for r in out)))
        say("P3a |signed delta sum| < %.2f main at every cell: %s   (%s)"
            % (P3A_FACTOR, v(P3a), ", ".join("%+.4f" % r["d_signed_over_main"] for r in out)))
        say("P3b triangle sum |delta| > main at 2^24: %s   (%s)"
            % (v(P3b), ", ".join("%.3f" % r["d_tri_over_main"] for r in out)))
        say("P4  R p_1 < 0 and ratio in (%.2f, %.2f) at every cell: %s   (%s)"
            % (P4_LO, P4_HI, v(P4), ", ".join("%.4f" % r["ratio"] for r in out)))
        say("")
        say("the lower bound the proof delivers, against the cell:  |p_1| vs N/(9 R L^2)")
        say("  e |        |p_1| |  N/(9 R L^2) | ratio | the (1/5)(N/R)L^-2 bound | ratio")
        for r in out:
            a3f = 0.2 * r["N"] / (r["R"] * r["L"] ** 2)
            say(" %2d | %12.6g | %12.6g | %5.2f | %19.6g | %5.2f"
                % (r["e"], abs(r["p1"]), r["norm_floor"], abs(r["p1"]) / r["norm_floor"], a3f,
                   abs(r["p1"]) / a3f))
    json.dump(dict(cells=out, clauses=dict(C1=C1, C2=C2, R1=R1, P1=P1, P2=P2, P3a=P3a, P3b=P3b,
                                           P4=(None if not (C1 and C2) else P4)),
                   A3_ratios=A3_RATIOS, C2_TOL=C2_TOL),
              io.open(OUTJ, "w", encoding="utf-8", newline="\n"), indent=1)
    _f.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

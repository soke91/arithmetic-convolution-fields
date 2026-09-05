# -*- coding: utf-8 -*-
r"""
The placebo control for "the mask exists" -- wall_v3.md, Section
{#sec:floor}, run against Lemma {#lem:placebo}.

WHY THIS HAS TO BE RUN

Section [sec:floor] reports that the cell means clear the exact floor
by a wide margin -- max_c |z_c| between 9.1 and 13.0 over every octave,
carried by the deep cells -- and adds that "under the permutation of
Lemma [lem:placebo] the floor collapses to the independent-sign value,
so what is being detected is the correspondence between cells and
divisibility and not the cell sizes."

That permutation control is CLAIMED and was never run: lab_cell_floor
computes the floor and the z-scores but states in its own NULL: line
that the placebo is not run there.  After the level measurement of
Remark [rem:levelmeas] turned out to be a support statement once its
null was run, the remaining un-nulled detection claim in either paper
is this one, and it is the strongest positive claim in the program.

The control: replace the labelling l(N) by l(pi(N)) for a random
permutation pi of the even N in the band.  Cell sizes are preserved
exactly; the correspondence between a cell and the arithmetic of its
members is destroyed.  Everything else -- Z(N), V(N), the band, the
exact floor formula of Lemma [lem:cellmom] -- is held identical, and
the floor is RECOMPUTED for each permutation rather than reused, since
u_c(v) depends on which N sit in c.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  L1  With the true labelling on the octave (2e6, 4e6], max_c |z_c| is
      at least 5.
  L2  Over 10 label permutations, max_c |z_c| stays below 5 every time,
      and its mean over permutations is below 3.
  L3  The floor itself is a property of the cell sizes: the standard
      error se_c changes by less than 10% between the true labelling
      and the permutation mean, at every depth.
  L4  The true labelling's z_c is monotone decreasing in depth; under
      permutation the rank correlation between depth and z_c has mean
      near zero, |mean| < 0.3.  Both halves are measured: the run first
      checked neither -- it printed a Pearson correlation where the
      rank correlation was registered, and never looked at monotonicity
      at all, so L4 reported "hold" without testing what it says.
  L5  The floor under permutation is a closed form.  For a uniformly
      random cell of size n_c in a band of size n,

        E_pi Var_eps(m_c - mbar) = (n - n_c)/(n_c (n-1)) (1 - Q_aa/n^2),

      Q_aa being the one term of Lemma [lem:cellmom] that does not
      depend on the cell.  The mean floor over the 10 permutations,
      se_perm, stands to sqrt of that within 25% at every depth.  The
      cap is wide on purpose: it is set above the 1/sqrt(2*10) ~ 22%
      that ten draws of a variance are worth before their spread is
      measured, and it is still far inside the factors 3.8 to 105 by
      which the true labelling's floor exceeds the closed form, so a
      "hold" cannot be produced by the arithmetic cell.  The ten
      draws' own scatter (SE10) is printed beside it and not gated.
      Registered with it, and not a prediction but a computation: the
      numerator m_c - mbar of the true labelling has, under the same
      permutation, mean 0 and variance (n - n_c)/(n_c (n-1)) s^2 with
      s^2 the variance of Z over the band, so Chebyshev bounds its
      permutation p-value by tau_c^2/(m_c - mbar)^2 with nothing
      sampled.

REFUTATION RULE (fixed before the run)

  L1  REFUTED if max_c |z_c| < 5 for the true labelling.
  L2  REFUTED if any permutation reaches 5, or if the mean is 3 or
      more.  This is the one that matters: if a permutation reproduces
      the detection, the mask is a property of the cell sizes and
      Section [sec:floor]'s central claim must be withdrawn.
  L3  REFUTED if any depth moves by 10% or more.
  L4  REFUTED if the permuted mean rank correlation has |mean| >= 0.3.
  L5  REFUTED if at any depth |se_perm/sqrt(A) - 1| >= 0.25.  A
      refutation here is of the closed form, not of the placebo: the
      lemma that replaces the ten draws in the paper would be wrong.

  L1, L2, L3 and L5 gate.  L4 is reported.

BACKS: Proposition {#prop:placebo} and Remark {#rem:floorsignal} in
paper/wall_v3.md, and Lemma {#lem:placebo}, whose control it runs.
"""

import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "lab_mask_placebo.txt")

LO, HI = 2_000_000, 4_000_000
CELLP = (3, 5, 7, 11, 13)
PERMS = 10
SEED = 20260808


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def rank_corr(x, y):
    """Spearman: 순위로 바꾼 뒤의 피어슨. 동순위는 평균순위로 준다."""
    def ranks(v):
        v = np.asarray(v, dtype=float)
        order = np.argsort(v, kind="mergesort")
        r = np.empty(v.size, dtype=float)
        r[order] = np.arange(1, v.size + 1, dtype=float)
        # 동순위 평균
        s = np.sort(v)
        i = 0
        while i < s.size:
            j = i
            while j + 1 < s.size and s[j + 1] == s[i]:
                j += 1
            if j > i:
                r[np.isclose(v, s[i])] = (i + j + 2) / 2.0
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y)
    if rx.std() == 0 or ry.std() == 0:
        return 0.0
    return float(np.corrcoef(rx, ry)[0, 1])


def pow2(n):
    L = 1
    while L < n:
        L <<= 1
    return L


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    X = HI
    say("sieving to %d ..." % X)
    pr = primes_upto(X)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(X + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > X:
            break
        q = p * p
        while q <= X:
            lam[q] = lgp[i]
            if q > X // p:
                break
            q *= p

    mu = np.ones(X + 1, dtype=np.int8)
    rem = np.arange(X + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(X))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= X:
            mu[p * p::p * p] = 0
        q = p
        while q <= X:
            rem[q::q] //= p
            if q > X // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    sqf = (mu != 0)

    say("convolving V and C ...")
    n = pow2(2 * (X + 1))
    a = np.zeros(n, dtype=np.float64)
    a[:X + 1] = lam ** 2
    FL2 = np.fft.rfft(a)
    a[:] = 0.0
    a[:X + 1] = sqf
    V = np.fft.irfft(FL2 * np.fft.rfft(a), n)[:X + 1]
    del FL2
    a[:] = 0.0
    a[:X + 1] = lam
    FL = np.fft.rfft(a)
    a[:] = 0.0
    a[:X + 1] = mu
    C = np.fft.irfft(FL * np.fft.rfft(a), n)[:X + 1]
    del a, FL

    depth = np.zeros(X + 1, dtype=np.int8)
    for p in CELLP:
        depth[p::p] += 1

    Ns = np.arange(LO + 2, HI + 1, 2, dtype=np.int64)
    Ns = Ns[V[Ns] > 0]
    g = 1.0 / np.sqrt(V[Ns])
    Z = C[Ns] / np.sqrt(V[Ns])
    nb = Ns.size
    zbar = float(Z.mean())
    say("band (%d, %d]: %d even N" % (LO, HI, nb))

    m = pow2(2 * (HI + 1))
    b = np.zeros(m, dtype=np.float64)
    b[:HI + 1] = lam[:HI + 1]
    FLam = np.conj(np.fft.rfft(b))
    w = sqf[:HI + 1].astype(np.float64)

    def ucorr(sel):
        b[:] = 0.0
        b[Ns[sel]] = g[sel]
        return np.fft.irfft(FLam * np.fft.rfft(b), m)[:HI + 1]

    ua = ucorr(np.ones(nb, dtype=bool))
    Qaa = float((w * ua * ua).sum())

    TERMS = {}

    def run(lab):
        """z_c and se_c for a labelling, floor recomputed from scratch."""
        zs, ses, ds, ns = [], [], [], []
        for d in range(6):
            sel = lab == d
            nc = int(sel.sum())
            if nc == 0:
                continue
            uc = ucorr(sel)
            Qcc = float((w * uc * uc).sum())
            Qca = float((w * uc * ua).sum())
            var = Qcc / nc ** 2 - 2.0 * Qca / (nc * nb) + Qaa / nb ** 2
            TERMS[d] = (Qcc / nc ** 2, Qca / (nc * nb),
                        Qaa / nb ** 2, var)
            se = math.sqrt(max(var, 0.0))
            mc = float(Z[sel].mean())
            zs.append((mc - zbar) / se if se > 0 else 0.0)
            ses.append(se)
            ds.append(d)
            ns.append(nc)
        return np.array(ds), np.array(ns), np.array(zs), np.array(ses)

    true_lab = depth[Ns]
    ds, ncs, z_true, se_true = run(true_lab)
    KTERMS = dict(TERMS)
    say()
    say("K's three terms, cells as arithmetic gives them  (POST HOC:")
    say("written after the rules below were scored, and scores none)")
    say("  Lemma [lem:cellmom] is doubly centred, so the variance is")
    say("  E_same,c - 2 E_c,a + E_all and not the difference of the")
    say("  outer two.  How much the cross term matters for K:")
    say("  depth  E_same,c   E_c,a      E_all      2-term     Var"
        "        2-term/Var")
    for d in sorted(KTERMS):
        es, eca, ea, vr = KTERMS[d]
        say("  %-6d %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %.4f"
            % (d, es, eca, ea, es - ea, vr, (es - ea) / vr))
    _lo = min(KTERMS[d][1] for d in KTERMS)
    _hi = max(KTERMS[d][1] for d in KTERMS)
    say("  E_c,a runs %.6f to %.6f against E_all = %.6f, and the"
        % (_lo, _hi, KTERMS[min(KTERMS)][2]))
    say("  two-term form misstates Var by factors %.3f to %.3f."
        % (min((KTERMS[d][0] - KTERMS[d][2]) / KTERMS[d][3]
               for d in KTERMS),
           max((KTERMS[d][0] - KTERMS[d][2]) / KTERMS[d][3]
               for d in KTERMS)))
    say("  So for K the cross term is not the band mean.  For the")
    say("  singular series it is: lab_cell_singular.py measures that,")
    say("  and there the two-term D_c and the doubly centred B_c[S_2]")
    say("  agree.  The parallel drawn in [meas:Dc] is between two")
    say("  doubly centred means; it is not a parallel of differences.")
    say()
    say("true labelling")
    say("  depth  n_c        se_c          z_c")
    for i in range(len(ds)):
        say("  %-6d %-10d %-13.4e %+.4f" % (ds[i], ncs[i], se_true[i],
                                            z_true[i]))
    mx_true = float(np.abs(z_true).max())
    say("  max |z_c| = %.4f" % mx_true)
    l1 = mx_true >= 5.0
    say("  L1 %s" % ("hold" if l1 else "REFUTED"))

    rng = np.random.default_rng(SEED)
    say()
    say("placebo: %d label permutations (cell sizes preserved exactly)"
        % PERMS)
    say("  draw   max |z_c|   z by depth")
    mxs, rhos, srhos, ses_p = [], [], [], []
    for t in range(PERMS):
        lab = rng.permutation(true_lab)
        d2, n2, z2, s2 = run(lab)
        mxs.append(float(np.abs(z2).max()))
        ses_p.append(s2)
        rhos.append(float(np.corrcoef(d2.astype(float), z2)[0, 1]))
        srhos.append(rank_corr(d2.astype(float), z2))
        say("  %-6d %-11.4f %s"
            % (t + 1, mxs[-1], " ".join("%+.2f" % v for v in z2)))
    mxs = np.array(mxs)
    l2 = bool((mxs < 5.0).all() and mxs.mean() < 3.0)
    say("  max over draws = %.4f;  mean = %.4f   (caps 5 and 3)   %s"
        % (mxs.max(), mxs.mean(), "hold" if l2 else "REFUTED"))

    se_p = np.mean(np.array(ses_p), axis=0)
    dev = np.abs(se_p - se_true) / se_true
    l3 = bool((dev < 0.10).all())
    say()
    say("L3   the floor under permutation, by depth")
    say("  depth  se_true       se_perm       |rel dev|")
    for i in range(len(ds)):
        say("  %-6d %-13.4e %-13.4e %.4f"
            % (ds[i], se_true[i], se_p[i], dev[i]))
    say("  L3 %s   (cap 0.10)" % ("hold" if l3 else "REFUTED"))

    # L5.  Q_aa is the third term of run()'s variance -- the one term
    # of Lemma [lem:cellmom] that does not depend on the cell -- and is
    # the Qaa computed above on this band, the same construction
    # lab_cell_floor.py makes on its own band.  The closed form was
    # derived and checked outside this script (permfloor: enumeration
    # of every n_c-subset at n = 9, a band (1e4, 2e4] with the kernel
    # built explicitly, and this octave against 200 fresh draws); here
    # it is arithmetic on Qaa and the cell sizes, nothing is read.
    #
    #   A_c = E_pi Var_eps(m_c - mbar)
    #       = (n - n_c)/(n_c (n-1)) * (1 - Q_aa/n^2)
    #
    # se_perm is the mean of sqrt(var) over draws while the lemma gives
    # the mean of var; Jensen's gap is cv(var)^2/8, of order 1e-3 for
    # the spread the draws show, so both ratios are printed.  SE10 is
    # the draws' own standard error, sd(se over draws)/sqrt(PERMS).
    nb_f = float(nb)
    ncs_f = ncs.astype(float)
    A = (nb_f - ncs_f) / (ncs_f * (nb_f - 1.0)) * (1.0 - Qaa / nb_f ** 2)
    sqA = np.sqrt(A)
    ses_arr = np.array(ses_p)
    var_p = np.mean(ses_arr ** 2, axis=0)
    ratio5 = se_p / sqA
    se10 = ses_arr.std(axis=0, ddof=1) / math.sqrt(PERMS)
    dev5 = np.abs(ratio5 - 1.0)
    l5 = bool((dev5 < 0.25).all())
    collapse = se_true / sqA
    say()
    say("L5   the floor under permutation against its closed form, by depth")
    say("     Q_aa/n^2 = %.6f   1 - Q_aa/n^2 = %.6f   n = %d"
        % (Qaa / nb_f ** 2, 1.0 - Qaa / nb_f ** 2, nb))
    say("  depth  n_c        se_perm       sqrt(A)       se_perm/sqrt(A)"
        "  var_perm/A  |dev|/SE10  se_true/sqrt(A)")
    for i in range(len(ds)):
        say("  %-6d %-10d %-13.4e %-13.4e %-16.4f %-11.4f %-11.2f %.2f"
            % (ds[i], ncs[i], se_p[i], sqA[i], ratio5[i], var_p[i] / A[i],
               abs(se_p[i] - sqA[i]) / se10[i], collapse[i]))
    say("  collapse factors se_true/sqrt(A): %s"
        % " ".join("%.2f" % v for v in collapse))
    say("  the floor of the arithmetic cell exceeds the closed form by "
        "factors from %.1f to %.1f" % (collapse.min(), collapse.max()))
    say("  L5 %s   (cap 0.25 on |se_perm/sqrt(A) - 1|; SE10 reported, "
        "not gated)" % ("hold" if l5 else "REFUTED"))

    # The numerator's permutation law is exact and needs no draw:
    # m_c - mbar under l(pi(N)) is the mean of a uniform n_c-subset of
    # the fixed multiset {Z(N)} minus its mean, so it has mean 0 and
    # variance B_c = (n - n_c)/(n_c (n-1)) * s^2, s^2 the variance of Z
    # over the band -- the without-replacement formula.  Chebyshev,
    # P(|t| >= k tau) <= 1/k^2, uses nothing beyond B_c.
    varZ = float(((Z - zbar) ** 2).mean())
    B = (nb_f - ncs_f) / (ncs_f * (nb_f - 1.0)) * varZ
    tau = np.sqrt(B)
    num = np.array([float(Z[true_lab == d].mean()) - zbar for d in ds])
    kk = np.abs(num) / tau
    cheb = np.minimum(1.0, 1.0 / kk ** 2)
    say()
    say("numerator of the true labelling against its exact permutation law")
    say("  s^2(Z) = %.6f   s^2/(1 - Q_aa/n^2) = %.6f"
        % (varZ, varZ / (1.0 - Qaa / nb_f ** 2)))
    say("  depth  n_c        m_c - mbar    tau_c = sqrt(B)  |t|/tau    "
        "Chebyshev")
    for i in range(len(ds)):
        say("  %-6d %-10d %-+13.5f %-16.4e %-10.2f %.2e"
            % (ds[i], ncs[i], num[i], tau[i], kk[i], cheb[i]))
    say("  union over %d cells of the Chebyshev bounds: %.2e"
        % (len(ds), min(1.0, float(cheb.sum()))))

    rho_true = float(np.corrcoef(ds.astype(float), z_true)[0, 1])
    srho_true = rank_corr(ds.astype(float), z_true)
    rm = float(np.mean(rhos))
    srm = float(np.mean(srhos))
    # L4 은 두 조각으로 등록됐다: 참 라벨의 z 가 깊이에 대해 단조
    # 감소한다는 것과, 순열 아래 순위상관의 평균이 0 근처라는 것.
    # 앞의 조각은 구현된 적이 없고 뒤의 조각은 순위가 아니라 피어슨으로
    # 재어졌다. 둘 다 등록된 대로 잰다.
    steps = np.diff(z_true)
    mono = bool((steps <= 0).all())
    up = [(int(ds[i]), float(steps[i])) for i in range(steps.size)
          if steps[i] > 0]
    l4 = mono and abs(srm) < 0.3
    say()
    say("L4   depth-vs-z: true Pearson %.4f, true Spearman %.4f;"
        % (rho_true, srho_true))
    say("     permuted mean Pearson %.4f, permuted mean Spearman %.4f "
        "(cap 0.3)" % (rm, srm))
    say("     z monotone decreasing in depth: %s" % mono)
    if not mono:
        say("     rises at depth %s"
            % ", ".join("%d->%d (+%.4f)" % (d, d + 1, v) for d, v in up))
    say("     L4 %s" % ("hold" if l4 else "REFUTED"))
    if not l4:
        say()
        say("  DIAGNOSTIC (post hoc).  The permutation half of L4 holds:")
        say("  the rank correlation the placebo produces is near zero.")
        say("  The monotonicity half is false, and it is false at the")
        say("  end where nothing is claimed -- depth 0 and 1 are the")
        say("  cells with no effect to order, and the correspondence")
        say("  L1-L3 measure is unaffected by their order.  The")
        say("  prediction was written as if the effect had to grow with")
        say("  depth monotonically; [rem:floorsignal] says the index is")
        say("  coarse and both ends concentrate the shift.  The")
        say("  threshold is not rewritten.")

    say()
    say("=" * 70)
    ok = l1 and l2 and l3 and l5
    say("L1 %s  L2 %s  L3 %s  L4 %s"
        % tuple("hold" if v else "REFUTED" for v in (l1, l2, l3, l4)))
    say("L5 %s" % ("hold" if l5 else "REFUTED"))
    # 결정적인 것은 L2 다 -- 등록문이 그렇게 적는다. L2 가 서면 탐지는
    # 셀 크기의 성질이 아니고, 그것이 이 절의 주장이다. L3·L4 의 반증은
    # 각각 바닥이 신호라는 것과 단조성 예측이 틀렸다는 것이지, 그 주장의
    # 반증이 아니다. 한 줄로 뭉뚱그리면 어느 쪽인지 못 읽는다.
    if ok:
        say("the mask survives its placebo: the detection is the "
            "cell-arithmetic correspondence")
    elif l2:
        say("the mask survives the permutation that decides it (L2); "
            "refuted above are %s"
            % ", ".join(n for n, v in (("L1", l1), ("L3", l3), ("L4", l4),
                                       ("L5", l5))
                        if not v))
    else:
        say("REFUTED -- Section {#sec:floor}'s central claim does not "
            "survive the permutation it cites")

    head = [
        "STATISTIC: z_c = (m_c - mbar)/se_c with se_c the exact floor of",
        "           Lemma {#lem:cellmom}, computed for the true depth",
        "           labelling and for 10 random permutations of that",
        "           labelling over the band; max_c |z_c| in each case; the",
        "           floor se_c itself under permutation against the true",
        "           one; and the correlation between depth and z_c.",
        "           Also the mean floor under permutation against its",
        "           closed form (L5), and the numerator of the true",
        "           labelling against its exact permutation variance.",
        "NULL: this file is the null -- the placebo of Lemma",
        "      {#lem:placebo}, which Section {#sec:floor} cites but which",
        "      no script had run. Cell sizes are preserved exactly and the",
        "      floor is recomputed from scratch for every permutation, so",
        "      the cell-to-arithmetic correspondence is the only thing",
        "      destroyed.",
        "FIELD: the octave (2e6, 4e6], even N with V(N) > 0; cells indexed",
        "       by depth = #{p in 3,5,7,11,13 dividing N}; Lambda, mu and",
        "       the squarefree indicator from an integer sieve to 4e6;",
        "       V = mu^2 * Lambda^2 and C = mu * Lambda by exact FFT",
        "       convolution; u_c by FFT cross-correlation; numpy",
        "       default_rng seed 20260808.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not ok:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())

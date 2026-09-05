# -*- coding: utf-8 -*-
r"""
What Proposition {#prop:dilate} does to E_3: the weights turn
nonnegative, and cross-k cancellation is not available.

WHAT IS AT STAKE

E_3(alpha) = sum_{k<K,(k,N)=1} mu(k) (log k) E_mu(N;k), and
Proposition [prop:dilate] gives E_mu(N;k) = mu(k)H(N;k) - C(N)/phi(k).
Substituting, mu(k) mu(k) = mu^2(k) = 1 on the squarefree k, so

    E_3(alpha) = sum_{k<K,(k,N)=1} mu^2(k) (log k) H(N;k)
                 - C(N) B_log(K),

    B_log(K) = sum_{k<K,(k,N)=1} mu(k) log k / phi(k)  ->  -S(N).

The signs mu(k) that Huang-Li discard by the triangle inequality, and
that this program has repeatedly said are what its results exploit,
CANCEL EXACTLY against the mu(k) inside the progression sum.  What is
left is a sum of dilated walls with WEIGHTS THAT ARE NONNEGATIVE.

The consequence is sharp.  Whatever smallness E_3 has cannot come from
the weights: it must come from the signs of H(N;k) themselves, across
k.  If those signs behaved like a coin, the sum of K terms would gain
a factor sqrt(K) over the sum of absolute values.  Measured earlier,
|E_3|/B(N) ran 0.35 to 0.54, i.e. a gain near 2 where sqrt(K) is 45 --
so the H(N;k) do not cancel across k, and each |H(N;k)| has to be
small on its own.  That is a much narrower target than "B(N) small",
and it closes off the hope that signed cancellation across dilations
could rescue the one-sided route.

For a coin there is no such collapse: eps(k) eps(mk) does not simplify,
the weights stay signed, and cross-k cancellation IS available.  So
this is also a structural discriminator, not a size one.

BACKS: Proposition {#prop:posweights} and Remark {#rem:nocrossk} in
paper/wall_v3.md.

T1 has content only if H(N;k) on the right is built from its own
definition -- the m-sum of the dilate, H(N;k) = sum_{m<N/k,(m,k)=1}
Lambda(N-mk) mu(m) -- and not as mu(k)A(N;k) from the progression
sum A that the left side already uses: sum mu(k)(log k)(A - C/phi)
- (sum (log k) mu(k)A - C sum mu(k)(log k)/phi) is the distributive
law and vanishes whatever the arrays contain.  It was written that
way, and could not fail.  The m-sum is now formed per k from Lambda
at N - mk and mu at m with the coprimality imposed explicitly, so the
residual measures the proposition's identity A(N;k) = mu(k)H(N;k).
Dropping the coprimality, shifting m by one, or signing the m-sum by
mu(k) each refutes T1 and leaves T2, T3, T4 as they were.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  T1  The identity E_3 = sum mu^2(k)(log k)H(N;k) - C(N)B_log(K) holds
      to better than 1e-12 relative, at every N tested.
  T2  The cross-k gain G_mu = sum(log k)|H| / |sum(log k)H| lies in
      [1.5, 5] at every N -- O(1), not sqrt(K).
  T3  For a coin the gain is far larger: G_eps > 5 G_mu at every N.
  T4  The absolute sum is dominated by the largest dilates: the top
      decile of k by |H| carries more than 40% of sum(log k)|H|.

REFUTATION RULE (fixed before the run)

  T1  REFUTED if the relative error reaches 1e-12 at any N.  This is an
      identity; a failure is an error in the derivation.
  T2  REFUTED if G_mu leaves [1.5, 5] at any N.
  T3  REFUTED if G_eps <= 5 G_mu at any N.
  T4  REFUTED if the top-decile share is 40% or below at any N.

  All four gate.
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
OUT = os.path.join(ROOT, "results", "lab_positive_weights.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
THETA = 0.56
DRAWS = 4
SEED = 20260808


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def sieves(n):
    pr = primes_upto(n)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(n + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lgp[i]
            if q > n // p:
                break
            q *= p
    mu = np.ones(n + 1, dtype=np.int8)
    rem = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            rem[q::q] //= p
            if q > n // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return pr, lam, mu


def phi_of(k):
    v, phi, d = k, 1, 2
    while d * d <= v:
        if v % d == 0:
            phi *= (d - 1)
            v //= d
            while v % d == 0:
                phi *= d
                v //= d
        d += 1
    if v > 1:
        phi *= (v - 1)
    return phi


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)
    rng = np.random.default_rng(SEED)

    say()
    say("  N          K      #k     |E_3 - identity| / |E_3|   G_mu     "
        "G_eps     sqrt K   top-decile share")
    say("  " + "-" * 96)
    t1 = t2 = t3 = t4 = True
    diag = []
    # sqrt(#k) is the l1/l2 of a vector whose entries all have the SAME
    # magnitude.  These do not, so the count-based reference overstates
    # what independent signs would have bought.  The magnitude-based
    # reference is the actual l1/l2 of (log k)H(N;k); Cauchy-Schwarz
    # caps it at sqrt(#k), so the ratio below cannot exceed 1 and a
    # value above 1 would mean the arithmetic is wrong.
    ratios = []
    for N in NS:
        v, PN, d = N, set(), 2
        while d * d <= v:
            if v % d == 0:
                PN.add(d)
                while v % d == 0:
                    v //= d
            d += 1
        if v > 1:
            PN.add(v)
        K = int(N ** THETA)
        ks = np.array([k for k in range(2, K)
                       if mu[k] != 0 and all(k % q for q in PN)])
        lg = np.log(ks.astype(float))
        iph = np.array([phi_of(int(k)) for k in ks], dtype=np.float64)

        def prog(sig):
            f = np.zeros(N, dtype=np.float64)
            idx = np.arange(1, N, dtype=np.int64)
            f[1:] = lam[1:N] * sig[N - idx]
            Ctot = float(f.sum())
            A = np.empty(ks.size)
            for i, k in enumerate(ks):
                r = N % int(k)
                A[i] = f[r::int(k)].sum() if r else f[int(k)::int(k)].sum()
            return A, Ctot

        muf = mu.astype(np.float64)
        A, C = prog(muf)
        sgn = mu[ks].astype(np.float64)
        E3 = float((sgn * lg * (A - C / iph)).sum())
        # the identity: mu(k)A = mu^2(k) H = H on squarefree k
        H = sgn * A
        # H(N;k) by its own definition, the m-sum of the dilate: Lambda
        # at N - mk, mu at m, (m,k) = 1 imposed by hand.  It reads
        # neither f nor A.  m runs downward so that the terms are met
        # in the same order as the progression slice in prog(): the
        # residual then measures the identity mu(k)H = A and not the
        # order of floating-point addition.
        Hm = np.empty(ks.size)
        for i, k in enumerate(ks):
            k = int(k)
            ms = np.arange(N // k, 0, -1, dtype=np.int64)
            s = mu[ms].astype(np.float64)
            v, p = k, 2
            while p * p <= v:
                if v % p == 0:
                    s[ms % p == 0] = 0.0
                    while v % p == 0:
                        v //= p
                p += 1
            if v > 1:
                s[ms % v == 0] = 0.0
            Hm[i] = float((lam[N - k * ms] * s).sum())
        Blog = float((sgn * lg / iph).sum())
        rhs = float((lg * Hm).sum()) - C * Blog
        rel = abs(E3 - rhs) / max(abs(E3), 1e-300)
        if rel >= 1e-12:
            t1 = False

        num = abs(float((lg * H).sum()))
        den = float((lg * np.abs(H)).sum())
        G_mu = den / max(num, 1e-300)
        if not (1.5 <= G_mu <= 5.0):
            t2 = False

        Gs = []
        for t in range(DRAWS):
            sig = np.zeros(N + 1, dtype=np.float64)
            supp = np.flatnonzero(mu[:N + 1] != 0)
            sig[supp] = rng.integers(0, 2, size=supp.size) * 2.0 - 1.0
            Ae, Ce = prog(sig)
            se = sig[ks]
            Gs.append(float((lg * np.abs(Ae)).sum())
                      / max(abs(float((se * lg * Ae).sum())), 1e-300))
        G_ep = float(np.mean(Gs))
        if G_ep <= 5.0 * G_mu:
            t3 = False

        contrib = lg * np.abs(H)
        order = np.argsort(-contrib)
        top = order[:max(1, ks.size // 10)]
        share = float(contrib[top].sum() / contrib.sum())
        if share <= 0.40:
            t4 = False

        say("  %-10d %-6d %-6d %-26.3e %-8.3f %-9.1f %-8.1f %.4f"
            % (N, K, ks.size, rel, G_mu, G_ep, math.sqrt(ks.size), share))
        wpos = float((lg * (H > 0)).sum() / lg.sum())
        diag.append((N, G_mu, ks.size, wpos))
        # den is already sum |(log k)H|, the l1 of the same vector.
        l2 = float(np.sqrt(((lg * H) ** 2).sum()))
        ratios.append((den / max(l2, 1e-300)) / math.sqrt(ks.size))

    say()
    say("T1  the identity holds to 1e-12 relative at every N   %s"
        % ("hold" if t1 else "REFUTED"))
    say("T2  G_mu in [1.5, 5] at every N                       %s"
        % ("hold" if t2 else "REFUTED"))
    # POST HOC, added after T3 was scored and scoring nothing.  T3 IS A
    # RULE THAT CANNOT FAIL.  Independent signs force the gain to be
    # l1/l2, and 5*G_mu is below l1/l2 at every N here, so a coin that
    # behaves like a coin passes T3 automatically.  The rule was written
    # against the count sqrt(#k), which is the same equal-magnitude
    # idealisation this file corrects two blocks below.
    #
    # It is not rescored: a registration is not edited after the run.
    # What is added is the threshold the null actually forces, so that a
    # reader can see the rule had no discriminating power.
    say("T3  G_eps > 5 G_mu at every N                         %s"
        % ("hold" if t3 else "REFUTED"))
    say("T4  top decile carries more than 40%%                  %s"
        % ("hold" if t4 else "REFUTED"))
    say()
    say("  T3 COULD NOT HAVE FAILED (post hoc, scoring nothing). Its")
    say("  threshold 5*G_mu is below what independent signs FORCE, which")
    say("  is l1/l2 on these magnitudes, at every N:")
    say("    N          5*G_mu     l1/l2 (forced)   threshold below?")
    for (N, gmu, nk, _), r in zip(diag, ratios):
        forced = r * math.sqrt(nk)
        say("    %-10d %-10.3f %-16.2f %s"
            % (N, 5.0 * gmu, forced, "yes" if 5.0 * gmu < forced else "NO"))
    say("  So a coin behaving as a coin passes T3 automatically. The rule")
    say("  was written against the count sqrt(#k) -- the same")
    say("  equal-magnitude idealisation corrected below. It is left as")
    say("  scored; a registration is not edited after its run.")
    say("  The estimator is a further reason not to read G_eps: it is a")
    say("  mean of l1/|sum| over four draws, and E[1/|X|] diverges for a")
    say("  Gaussian, so the large entries are one draw with a small")
    say("  denominator rather than signal.")
    say()
    say("  DIAGNOSTIC (post hoc). T4 guessed the wrong explanation. The")
    say("  top decile carries about 35%, concentrated but not dominant,")
    say("  so the missing sqrt(K) is not a heavy tail -- it is a genuine")
    say("  sign correlation among the H(N;k). Two ways to see it:")
    say("  N          n_eff = G^2   #k       (l1/l2)^2   wtd frac H > 0")
    for (N, gmu, nk, wf), r in zip(diag, ratios):
        say("  %-10d %-13.1f %-8d %-11.1f %.4f"
            % (N, gmu ** 2, nk, (r * math.sqrt(nk)) ** 2, wf))
    say("  n_eff is the number of independent signs the sum behaves as")
    say("  if it had. THE REFERENCE IS THE FOURTH COLUMN, NOT THE THIRD.")
    say("  n_eff = G^2 and the count #k are both the equal-magnitude")
    say("  idealisation; what independent signs force on THESE magnitudes")
    say("  is (l1/l2)^2. The conclusion is unchanged -- single digits")
    say("  against three figures either way -- but the count overstates")
    say("  the reference by the same factor the block below corrects.")

    say()
    say("  The sqrt K column above is a COUNT-based reference: it is the")
    say("  l1/l2 of a vector whose entries all share one magnitude. The")
    say("  entries here do not, so the magnitude-based reference -- the")
    say("  actual l1/l2 of (log k)H(N;k) -- is smaller, and the shortfall")
    say("  it leaves for the signs to explain is correspondingly smaller.")
    say("  N          l1/l2      sqrt #k    ratio")
    for (N, _, nk, _), r in zip(diag, ratios):
        say("  %-10d %-10.2f %-10.2f %.4f" % (N, r * math.sqrt(nk),
                                              math.sqrt(nk), r))
    say("REFERENCE lab_positive_weights %d %.4f %.4f"
        % (len(ratios), min(ratios), max(ratios)))

    say()
    say("=" * 70)
    ok = t1 and t2 and t3 and t4
    say("T1 %s  T2 %s  T3 %s  T4 %s"
        % tuple("hold" if v else "REFUTED" for v in (t1, t2, t3, t4)))
    say("the weights are nonnegative and cross-k cancellation is not "
        "available: each |H(N;k)| must be small on its own"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the relative error of E_3 = sum mu^2(k)(log k)H(N;k)",
        "           - C(N)B_log(K); the cross-k cancellation gain",
        "           G = sum(log k)|H| / |sum(log k)H| for mu and the",
        "           corresponding gain for coin signs, against sqrt of the",
        "           number of moduli; and the share of sum(log k)|H|",
        "           carried by the top decile of k.",
        "DENOM: the identity residual is divided by E_3(N);",
        "NULL: the coin is the control and is structural here rather than",
        "      numerical -- eps(k)eps(mk) does not collapse, so a coin",
        "      keeps signed weights and the cross-k cancellation that mu",
        "      loses. Four draws per N, same support and k-range.",
        "FIELD: N = 2e5, 4e5, 8e5, 1.6e6, 3.2e6 with theta' = 0.56, so k",
        "       runs over the squarefree k < N^0.56 coprime to N; m over",
        "       1 <= m < N/k; Lambda and mu from an integer sieve to",
        "       3.2e6; numpy default_rng seed 20260808.",
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

# -*- coding: utf-8 -*-
r"""The gap between the two cuts is the object called Delta.

Supports {#rem:cutbridge}.

WHAT IS AT STAKE

The companion run code/audit_cut_bridge.py asks whether the fixed-cut
statement reaches the moving-cut object by telescoping over the outer
index, and finds that it does so only at a cost that no unconditional
argument supplies.  That was the wrong bridge.  There is a direct one.

Write J = sum_{k<K, (k,N)=1} mu(k) Emu(N-1; k) for the fixed-cut object
at the extreme endpoint, and T1* = sum_k mu(k) Emu(Y_k; k) for the
moving-cut object, Y_k = ceil(N - alpha k) - 1.  Their difference is
carried by the n it drops:

    n > Y_k  <=>  n >= N - alpha k  <=>  N - n <= alpha k  <=>  m <= alpha,

m = (N-n)/k.  So on exactly those terms mu(k) mu(N-n) = mu(k) mu(mk) =
mu(m) mu^2(k): the Moebius factor lands on the SHORT variable and the
long variable carries only mu^2 >= 0.  That is the shape of the term
this note calls Delta, with weight 1 in place of log m.  Hence

    T1* = J - Xi,      Xi = Xi_main - Xi_mean,
    Xi_main = sum_{m<=alpha} mu(m) sum_{k<K, (k,mN)=1} mu^2(k) Lambda(N-mk),
    Xi_mean = sum_{k<K, (k,N)=1} mu(k)/phi(k) * ( C(N-1) - C(Y_k) ),

with C(t) = sum_{n<=t} Lambda(n) mu(N-n).  If this is an identity, then
the fixed-cut theorem is not superseded by the corrected formulation:
it is one of the two ingredients, and the term whose omission was the
defect in the published equation is the other.  The two cuts are a
Delta apart.

WHAT IS MEASURED

For N even and alpha = N^theta, with K = (N-1)/alpha:

  Z1  the bridge is exact: T1* against J - Xi_main + Xi_mean, where
      T1* is built from the cofactor cut (sum over k | N-n with
      (N-n)/k > alpha) and J from the full divisor count at truncation
      K, so the two sides share no accumulation.

  Z2  Xi is what it is claimed to be: the k-side difference
      sum_k mu(k) [ Emu(N-1;k) - Emu(Y_k;k) ], accumulated over residue
      classes, against Xi_main - Xi_mean, accumulated over m.

  Z3  the bridge term sits on the short variable: the largest m
      carrying a surviving term in Xi_main, against alpha.

  Z4  sizes, reported not judged: |J|/N, |Xi|/N, |T1*|/N.

FALSIFICATION, registered before the run

  Z1  REFUTED if the two sides disagree by more than 1e-12 relative to
      N.  Then the difference of the two cuts is not Delta's shape and
      the claim above is wrong.
  Z2  REFUTED at the same tolerance.  Then the split of Xi into a main
      term and a mean term is misstated.
  Z3  REFUTED if any surviving m exceeds alpha.  Then the bridge term
      is not on the short variable, it is not Delta's shape, and it
      closes by no argument given here.
  Z4  is reported, not judged.  Sizes at accessible N say nothing about
      an asymptotic bound.

  PREDICTION.  Z1-Z3 hold.  Z4 shows |Xi|/N of the same order as |J|/N
  and |T1*|/N -- the bridge is not a small correction to be waved
  through, which is the reason it needs Delta's estimate rather than a
  triangle inequality.

  Z3 IS VOID -- it could not fail as registered, and was for a while
  recorded as refuted.  The verdict block at the end of the run says
  why in full.  The two identities are Z1 and Z2.

NULL.  None applies: deterministic finite sums, no sampling, no sign
input.  The control is that each identity is accumulated on two sides
that share no index order -- divisor counts against residue classes in
Z2, cofactor cut against full truncation in Z1.
"""
# Why Z3 is void, at length.  This is below the docstring on purpose:
# the docstring is printed into the result file as the registration
# record, and the registration is left as it was registered.
#
# Delta is enumerated with m <= alpha, so no m in it can exceed alpha.
# The clause "Z3 REFUTED if any surviving m exceeds alpha" therefore
# compared the enumeration against its own bound and could not fail for
# any input -- erasing the arithmetic entirely (Lambda = mu = 1, so
# that nothing is filtered out) still gave hold.  What was then done
# was to carry the same survival test past alpha, over m <= 4*alpha,
# find surviving terms, and record Z3 refuted.  That is wrong in the
# other direction: the summand is defined for every m, but the terms
# with m > alpha are not terms of Delta, so evaluating them tests a
# quantity the rule was not about.  A rule that cannot fail is not
# repaired by enlarging the object until it can.
#
# The scan is kept and its numbers printed, because "the tail past
# alpha is not short" is a true statement about the summand.  It
# carries no verdict on Z3 and the exit status no longer depends on it.
# Nothing above or below rests on Z3; P1 sec:movingcut says so.
import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), "results")

NS = (200_000, 400_000, 800_000, 1_600_000, 3_200_000)
THETA = 0.44           # alpha = N^theta, as in audit_moving_switch
TOL = 1e-12            # relative to N


def sieve(n):
    """Lambda, mu, and the smallest-prime-factor table, to n."""
    spf = np.zeros(n + 1, dtype=np.int64)
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == 0:
            spf[i * i::i] = np.where(spf[i * i::i] == 0, i, spf[i * i::i])
    lam = np.zeros(n + 1, dtype=np.float64)
    mu = np.zeros(n + 1, dtype=np.int64)
    mu[1] = 1
    for v in range(2, n + 1):
        p = int(spf[v]) or v
        w = v // p
        mu[v] = 0 if w % p == 0 else -mu[w]
        t = v
        while t % p == 0:
            t //= p
        lam[v] = math.log(p) if t == 1 else 0.0
    return lam, mu, spf


def factor_set(v, spf):
    out = set()
    while v > 1:
        p = int(spf[v]) or v
        out.add(p)
        while v % p == 0:
            v //= p
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say(__doc__.strip())
    say()
    say("=" * 72)
    say("sieving to %d ..." % max(NS))
    lam, mu, spf = sieve(max(NS))

    say()
    say("STATISTIC: the difference between the fixed-cut object at the")
    say("           extreme endpoint and the moving-cut object, against")
    say("           the term of Delta's shape, each accumulated on two")
    say("           sides that share no index order; and the largest m")
    say("           carrying a surviving term in that difference.")
    say("FIELD: N even in %s; alpha = N^theta; k < K = (N-1)/alpha with"
        % ",".join(str(v) for v in NS))
    say("       (k,N) = 1 and k squarefree; n < N; Y_k = ceil(N-alpha k)-1;")
    say("       m = (N-n)/k. Lambda, mu from one integer sieve.")
    say("CONSTANTS: THETA = %.2f, TOL = %.0e relative to N" % (THETA, TOL))
    say("NULL: none applies -- deterministic finite sums, no sampling, no")
    say("      sign input. The control is two-sided accumulation.")
    say("DENOM: every relative error and every size printed is over N.")
    say()

    hdr = ("  %-9s %-6s %-12s %-12s %-8s %-10s %-10s %-10s"
           % ("N", "alpha", "Z1 rel", "Z2 rel", "max m", "|J|/N",
              "|Xi|/N", "|T1*|/N"))
    say(hdr)
    say("  " + "-" * (len(hdr) - 2))

    ok1 = ok2 = True
    for N in NS:
        alpha = N ** THETA
        K = (N - 1) / alpha
        jmax = int(math.ceil(K)) - 1
        PN = factor_set(N, spf)

        phi = np.arange(jmax + 2, dtype=np.float64)
        for p in range(2, jmax + 2):
            if phi[p] == p:
                phi[p::p] *= (1.0 - 1.0 / p)
        phi[1] = 1.0

        nn = np.arange(0, N, dtype=np.int64)
        a = lam[:N] * mu[N - nn].astype(np.float64)
        a[0] = 0.0
        C = np.cumsum(a)
        Ctot = math.fsum(a[1:].tolist())

        good = [k for k in range(1, jmax + 1)
                if mu[k] != 0 and not (factor_set(k, spf) & PN)]

        # --- J and T1*: divisor counts, two different truncations
        dK = np.zeros(N, dtype=np.int64)      # k | N-n, k < K
        ds = np.zeros(N, dtype=np.int64)      # k | N-n, (N-n)/k > alpha
        WK = 0.0
        Wpref = np.zeros(jmax + 2, dtype=np.float64)
        run = 0.0
        m0 = int(math.floor(alpha)) + 1
        for k in range(1, jmax + 1):
            run_add = 0.0
            if mu[k] != 0 and not (factor_set(k, spf) & PN):
                mm = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
                dK[N - k * mm] += int(mu[k])
                if m0 * k < N:
                    ds[N - k * mm[mm >= m0]] += int(mu[k])
                run_add = float(mu[k]) / phi[k]
            run += run_add
            Wpref[k] = run
        WK = run

        J = math.fsum((a[1:] * dK[1:]).tolist()) - WK * Ctot

        Rn = np.ceil((N - nn[1:]) / alpha).astype(np.int64) - 1
        np.clip(Rn, 0, jmax, out=Rn)
        T1s = (math.fsum((a[1:] * ds[1:]).tolist())
               - math.fsum((a[1:] * Wpref[Rn]).tolist()))

        # --- Xi from the k-side: residue classes, n above the moving cut
        xk_terms, xm_terms = [], []
        for k in good:
            Yk = math.ceil(N - alpha * k) - 1
            lo = max(Yk, 0)
            r = N % k
            start = r if r > 0 else k
            first = start + ((lo - start) // k + 1) * k if lo >= start else start
            ns = np.arange(first, N, k, dtype=np.int64)
            if ns.size:
                xk_terms.append(int(mu[k]) * math.fsum(a[ns].tolist()))
            xm_terms.append(float(mu[k]) / phi[k]
                            * (Ctot - (float(C[Yk]) if Yk >= 1 else 0.0)))
        Xi_raw = math.fsum(xk_terms)
        Xi_mean = math.fsum(xm_terms)

        # --- Xi_main from the m-side: Delta's shape, weight 1
        main_terms = []
        maxm = 0
        for m in range(1, int(alpha) + 1):
            if mu[m] == 0:
                continue
            fm = factor_set(m, spf) | PN
            kk = np.arange(1, min(jmax, (N - 1) // m) + 1, dtype=np.int64)
            keep = np.ones(kk.size, dtype=bool)
            for p in fm:
                keep &= (kk % p != 0)
            kk = kk[keep]
            kk = kk[mu[kk] != 0]
            if kk.size == 0:
                continue
            s = math.fsum(lam[N - m * kk].tolist())
            if s != 0.0:
                maxm = max(maxm, m)
            main_terms.append(int(mu[m]) * s)
        Xi_main = math.fsum(main_terms)

        Xi = Xi_main - Xi_mean
        r1 = abs(T1s - (J - Xi)) / N
        r2 = abs(Xi_raw - Xi_main) / N

        ok1 &= r1 <= TOL
        ok2 &= r2 <= TOL
        # Z3/X3 asks whether ANY surviving m exceeds alpha.  The loop
        # above stops at int(alpha), so `maxm <= alpha` re-checked the
        # loop's own bound: int(alpha) <= alpha always, and the branch
        # to REFUTED was unreachable for every input.  Erasing the
        # arithmetic entirely (Lambda = mu = 1, so that nothing is
        # filtered out) still gave hold.  THAT MAKES THE RULE VOID, and
        # the scan below does not repair it -- carrying the survival
        # test past alpha gives the clause somewhere to fail only by
        # testing terms that are not terms of Delta.  So the scan is
        # run and printed for what it does say (the tail past alpha is
        # not short) and no verdict is read off it.  The range is
        # stated because the rule is universal in m and this scan is
        # not.
        SCAN = 4
        over = 0
        for m in range(int(alpha) + 1, SCAN * int(alpha) + 1):
            if m >= N or mu[m] == 0:
                continue
            fm = factor_set(m, spf) | PN
            kk = np.arange(1, min(jmax, (N - 1) // m) + 1, dtype=np.int64)
            keep = np.ones(kk.size, dtype=bool)
            for p_ in fm:
                keep &= (kk % p_ != 0)
            kk = kk[keep]
            kk = kk[mu[kk] != 0]
            if kk.size == 0:
                continue
            s = math.fsum(lam[N - m * kk].tolist())
            if s != 0.0:
                over = max(over, m)
        # `over` is not compared against anything: see above, Z3 takes
        # no verdict from this scan.
        #
        # over saturates the scan bound at every N -- surviving
        # terms are found all the way to 4*alpha -- so it is a
        # LOWER bound on the excess and not a maximum.  Reporting
        # it as a maximum would repeat the very fault this scan
        # was added to fix: an enumeration bound read as a result.
        say("      Z3 scan: alpha = %.1f. Surviving terms are found"
            " above alpha throughout m <= %d*alpha; the highest"
            " one seen, %d, is the scan bound itself, so the tail"
            " is not short and this number is a lower bound on the"
            " excess rather than a maximum."
            % (alpha, SCAN, over))
        say("  %-9d %-6.0f %-12.3e %-12.3e %-8d %-10.6f %-10.6f %-10.6f"
            % (N, alpha, r1, r2, maxm, abs(J) / N, abs(Xi) / N,
               abs(T1s) / N))

    say()
    say("Z1  the two cuts differ by exactly the bridge term")
    say("    Z1 %s" % ("hold" if ok1 else "REFUTED"))
    say("Z2  the bridge term is the main term of Delta's shape, less its mean")
    say("    Z2 %s" % ("hold" if ok2 else "REFUTED"))
    say("Z3  the bridge term sits on the short variable")
    say("    Z3 VOID -- the rule asked whether any surviving m exceeds")
    say("    alpha, and Delta is enumerated with m <= alpha, so it could")
    say("    not fail.  It was for a while recorded as refuted, by")
    say("    carrying the survival test past alpha; the terms that finds")
    say("    are not terms of Delta.  No verdict either way.")
    say()
    say("Z4  sizes are in the table, reported not judged: nothing at these")
    say("    N bears on an asymptotic bound.")
    say()
    say("=" * 72)
    say("Z1 %s  Z2 %s"
        % tuple("hold" if v else "REFUTED" for v in (ok1, ok2)))
    say("Z3 VOID -- could not fail as registered; no verdict is taken.")

    io.open(os.path.join(RES, "audit_delta_bridge.txt"), "w",
            encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    return 0 if (ok1 and ok2) else 1


if __name__ == "__main__":
    raise SystemExit(main())

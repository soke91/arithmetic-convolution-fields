# -*- coding: utf-8 -*-
r"""
paper/theorem_A.md, Remark {#rem:nogoscope}.

WHAT IS UNDER TEST

[rem:dilaterange] puts both branches of EH_mu on the same objects: the
dilates H(N;k) of [eq:dilate], the flat branch over k >= K with weight
1 after the completion of Step 2, the log branch over k < K with
weight log k.  The companion no-go is a statement about weights on
exactly those dilates -- its T_w is sum_{k<K} w_k mu(k) A(N;k) and
A(N;k) = mu(k) H(N;k), so T_w is sum_{k<K} w_k mu^2(k) H(N;k).

That makes a question checkable that has been left open: which of the
two branches the no-go actually covers.  Its hypothesis is on the
transform b = mu * w, required to be supported in
[1, N^{theta'-1/2-delta}] so that the residual stays inside
Bombieri-Vinogradov's reach.  Two weights, two transforms:

    w_k = 1        gives  b = mu * 1 = delta,  supported at {1}
    w_k = log k    gives  b = mu * log = Lambda, supported on the
                          prime powers below K = N^{theta'}

The first is inside the hypothesis trivially.  The second is not: Lambda
reaches K, and K is far above N^{theta'-1/2-delta}.  If that is right
then the no-go covers the flat branch and not the log branch, which is
the opposite of the way the two are usually spoken of here -- the flat
branch is the one this note proves something about, and the log branch
is the one left open.

BACKS: Remark {#rem:nogoscope} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  TD1  For w = 1 the transform is the delta: b_d = 0 for every
       2 <= d <= 4000 and b_1 = 1.
  TD2  For w_k = log k the transform is von Mangoldt: b_d = Lambda(d)
       to 1e-12 relative for every d <= 4000.
  TD3  The supports fall on opposite sides of the no-go's hypothesis:
       at every N and theta' of the field, max{d : b_d != 0} is 1 for
       w = 1 and exceeds N^{theta'-1/2} for w = log k.
  TD4  The identity behind the design space holds: T_w computed as
       sum_{k<K} w_k mu(k) A(N;k) equals sum_{k<K} w_k mu^2(k) H(N;k)
       to 1e-9 relative, for both weights and every (N, theta').

TD5 WAS REGISTERED AFTER TD3 RAN AND WAS REFUTED, AND BEFORE TD5 RAN

TD3 put theta' = 0.44 in its field, and the companion no-go is stated
for theta' in (1/2, 1).  Below the threshold its support bound
N^{theta'-1/2-delta} is under 1, so no weight at all meets the
hypothesis and the design space is empty -- the delta itself, supported
at {1}, is already outside it.  That is why TD3 fails, and it is a fact
about where the no-go is stated rather than about the weights.  The
verdict stands and TD5 asks the question inside the no-go's own range.

  TD5  At the theta' above 1/2 in the field, the two supports fall on
       opposite sides of the bound: max{d : b_d != 0} is 1 for w = 1,
       which is inside, and exceeds N^{theta'-1/2} for w = log k, which
       is outside.

REFUTATION RULE (fixed before the run)

  TD1 REFUTED by a single d in range with b_d != 0, or b_1 != 1.
  TD2 REFUTED by a single d in range exceeding the tolerance.
  TD3 REFUTED by a single (N, theta') at which either support lands on
      the wrong side.
  TD4 REFUTED by a single point exceeding the tolerance.
  TD5 REFUTED by a single (N, theta') with theta' > 1/2 at which
      either support lands on the wrong side.
  Non-zero exit on any refutation.
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
OUT = os.path.join(ROOT, "results", "audit_nogo_scope.txt")

DLIM = 4000
NS = [400_000, 1_600_000]
THETAS = (0.44, 0.56)


def sieves(n):
    spf = np.zeros(n + 1, dtype=np.int64)
    for p in range(2, n + 1):
        if spf[p] == 0:
            blk = spf[p::p]
            spf[p::p] = np.where(blk == 0, p, blk)
    mu = np.ones(n + 1, dtype=np.int64)
    mu[0] = 0
    for v in range(2, n + 1):
        p = int(spf[v])
        w = v // p
        mu[v] = 0 if w % p == 0 else -mu[w]
    lam = np.zeros(n + 1, dtype=np.float64)
    for p in range(2, n + 1):
        if int(spf[p]) != p:
            continue
        q, lg = p, math.log(p)
        while q <= n:
            lam[q] = lg
            if q > n // p:
                break
            q *= p
    return spf, mu, lam


def transform(w, dlim, mu):
    """b = mu * w, i.e. b_d = sum_{e|d} mu(d/e) w(e)."""
    b = np.zeros(dlim + 1, dtype=np.float64)
    for d in range(1, dlim + 1):
        t = 0.0
        for e in range(1, d + 1):
            if d % e:
                continue
            t += float(mu[d // e]) * w(e)
        b[d] = t
    return b


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    top = max(NS)
    say("building mu and Lambda to %d ..." % top)
    spf, mu, lam = sieves(top)
    say()

    say("TD1/TD2  the two weights' transforms b = mu * w")
    say("=" * 70)
    b1 = transform(lambda e: 1.0, DLIM, mu)
    bl = transform(lambda e: math.log(e), DLIM, mu)
    td1 = (b1[1] == 1.0) and all(b1[d] == 0.0
                                 for d in range(2, DLIM + 1))
    worst = 0.0
    for d in range(1, DLIM + 1):
        worst = max(worst, abs(bl[d] - lam[d]) / max(lam[d], 1.0))
    td2 = worst <= 1e-12
    say("  w = 1:      b_1 = %.1f, nonzero b_d for 2<=d<=%d: %d"
        % (b1[1], DLIM,
           sum(1 for d in range(2, DLIM + 1) if b1[d] != 0.0)))
    say("  w = log k:  worst |b_d - Lambda(d)| relative over d<=%d: %.3g"
        % (DLIM, worst))
    say("  -> TD1 %s   TD2 %s" % ("HOLDS" if td1 else "REFUTED",
                                  "HOLDS" if td2 else "REFUTED"))

    say()
    say("TD3/TD4  which branch the no-go's hypothesis covers")
    say("=" * 70)
    say("  the hypothesis asks b supported in [1, N^{theta'-1/2-delta}].")
    say()
    say("  N        theta'   N^(t-1/2)   supp(b) w=1  supp(b) w=log"
        "   T_1 resid   T_log resid")
    say("  " + "-" * 76)
    td3 = td4 = td5 = True
    for N in NS:
        idx = np.arange(2, N)
        an = lam[2:N] * mu[N - idx].astype(np.float64)
        nz = np.nonzero(an)[0]
        ns, av = idx[nz], an[nz]
        for tp in THETAS:
            K = N ** tp
            bound = N ** (tp - 0.5)
            # supports
            s1 = 1
            slog = max(d for d in range(1, min(DLIM, int(K)) + 1)
                       if lam[d] > 0)
            if not (s1 <= bound and slog > bound):
                td3 = False
                if tp > 0.5:
                    td5 = False
            # TD4: the two forms of T_w
            t1a = t1b = tla = tlb = 0.0
            for k in range(1, int(K) + 1):
                if mu[k] == 0 or math.gcd(k, N) != 1:
                    continue
                A = float(av[(ns % k) == (N % k)].sum())
                ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
                sel = (np.gcd(ms, k) == 1) & (mu[ms] != 0)
                H = float((lam[N - ms[sel] * k]
                           * mu[ms[sel]]).sum()) if sel.any() else 0.0
                t1a += float(mu[k]) * A
                t1b += H                      # mu^2(k) = 1 here
                lg = math.log(k)
                tla += lg * float(mu[k]) * A
                tlb += lg * H
            r1 = abs(t1a - t1b) / max(abs(t1a), 1.0)
            rl = abs(tla - tlb) / max(abs(tla), 1.0)
            if r1 > 1e-9 or rl > 1e-9:
                td4 = False
            say("  %-8d %-8.2f %-11.1f %-12d %-14d %-11.3g %.3g"
                % (N, tp, bound, s1, slog, r1, rl))
    say("  -> TD3 %s   TD4 %s   TD5 %s"
        % tuple("HOLDS" if v else "REFUTED"
                for v in (td3, td4, td5)))
    say()
    say("  below 1/2 the bound N^{theta'-1/2} is under 1, so the no-go's")
    say("  design space is empty there -- the delta itself does not meet")
    say("  the hypothesis.  That is the whole of TD3's failure.")

    say()
    say("=" * 70)
    ok = (("TD1", td1), ("TD2", td2), ("TD3", td3), ("TD4", td4),
          ("TD5", td5))
    for k, v in ok:
        say("  %s %s" % (k, "holds" if v else "REFUTED"))
    verdict = all(v for _, v in ok)
    say()
    say("which of the two branches the companion no-go's hypothesis "
        "covers, in the dilate vocabulary")
    if not verdict:
        say("at least one pre-registered rule failed; it is named above")

    head = [
        "STATISTIC: the transform b = mu * w of the two weights the two",
        "           branches carry, w = 1 and w = log k, computed",
        "           elementwise for d <= 4000 and compared with the",
        "           delta and with Lambda; the largest d with b_d",
        "           nonzero for each, against the no-go's support bound",
        "           N^{theta'-1/2}; and T_w computed both as",
        "           sum_{k<K} w_k mu(k) A(N;k) and as",
        "           sum_{k<K} w_k mu^2(k) H(N;k), with the relative",
        "           difference, for both weights.",
        "DENOM: each difference is divided by max(the value compared",
        "       against, 1);",
        "NULL: none. Every quantity here is an identity being checked --",
        "      two convolutions against their known values, and two",
        "      forms of one sum against each other -- or a comparison of",
        "      a support with a stated bound. There is no statistic",
        "      being read against chance.",
        "FIELD: d <= 4000 for the transforms; N = 4e5 and 1.6e6 and",
        "       theta' in {0.44, 0.56}, one below the threshold and one",
        "       above, for the supports and the two forms of T_w; every",
        "       squarefree k < K coprime to N.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not verdict:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())

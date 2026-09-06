# -*- coding: utf-8 -*-
r"""
OPEN.md, "목표보다 강한 것을 공격하고 있는가" -- the direct route to
r~(N) > 0, and the size of the margin it runs on.

WHAT THIS MEASURES AND WHY

Theorem {#thm:C} is an identity:

    r~(N) = S(N)(N - C(N)) + E_3(alpha) + O_A(N (log N)^{-A}).

Huang-Li consume it through |E_3| small, which is two-sided and asks a
saving of every power of log.  But Goldbach needs only r~(N) > 0, and
the triangle bound

    |C(N)| <= U(N) := sum_{n<N} Lambda(n) mu^2(N-n) = A(N) N (1+o(1))

with A(N) < 1 already leaves a positive margin.  Substituting,

    r~(N) >= S(N) N (1 - A(N)) (1+o(1)) + E_3 + O_A(N (log N)^{-A}),

so binary Goldbach for large even N follows from the ONE-SIDED bound

    E_3(alpha) > - S(N) (1 - A(N)) N (1 + o(1)),                  (*)

which is weaker than |E_3| <<_A N (log N)^{-A} in two ways: it is
one-sided, and its threshold is S(N)(1-A(N)) N rather than
N (log N)^{-A}.  The whole content of the weakening is therefore the
size of S(N)(1-A(N)) -- how much room (*) actually has.  That is what
this script measures.  The quantity is smallest where N has many small
prime factors, since A(N) = Artin * prod_{q|N} (1-1/(q(q-1)))^{-1}
rises towards 1 there while S(N) rises only linearly in the same
factors.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  F1  argmin over even N <= 1.6e7 of S(N)(1-A(N)) is 9699690, the
      largest primorial below the range.
  F2  That minimum lies in [0.05, 0.08].
  F3  min over even N <= 1.6e7 of S(N)(1-A(N)) * log N * log log N
      exceeds 1, i.e. the margin is bounded below by c/(log N log log N)
      with c > 1 on this range.
  F4  The median of S(N)(1-A(N)) over even N <= 1.6e7 is at least 0.2:
      the margin is of order 1 for most N and only degenerates on a thin
      set.
  F5  U(N)/(A(N) N) -> 1: its mean over the top octave lies in
      [0.99, 1.01], and |C(N)| <= U(N) at every even N (an identity
      check on the triangle bound).
  F6  r~(N) >= S(N)(1-A(N)) N at every even N in [10^5, 1.6e7], and the
      minimum slack r~(N) / (S(N)(1-A(N))N) over that range is between
      3.5 and 4.5 -- i.e. about 1/(1-A) at the least favourable N.

REFUTATION RULE (fixed before the run)

  F1  REFUTED if the argmin is not 9699690.
  F2  REFUTED if the minimum falls outside [0.05, 0.08].
  F3  REFUTED if the product is at most 1 anywhere on the range.
  F4  REFUTED if the median is below 0.2.
  F5  REFUTED if the mean ratio leaves [0.99, 1.01], or if |C| > U
      anywhere.
  F6  REFUTED if r~ < S(1-A)N anywhere on [10^5, 1.6e7], or if the
      minimum slack leaves [3.5, 4.5].

  All six gate.  F1-F4 are the size of the margin in (*), F5 is the
  step that produces it, and F6 is the route's conclusion checked
  directly against the truth.

BACKS: Proposition {#prop:onesided} in paper/theorem_A.md -- the
size of the one-sided threshold S(N)(1-A(N))N that it runs on.

CITED BY: {#rem:onesided} in paper/.
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
OUT = os.path.join(ROOT, "results", "lab_onesided_margin.txt")

X = 16_000_000
LOW = 100_000
PRIMORIAL = 9_699_690


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


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

    say("building A(N) and S(N) ...")
    artin, twin = 1.0, 2.0
    for p in pr:
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    AN = np.full(X + 1, artin, dtype=np.float64)
    SN = np.full(X + 1, twin, dtype=np.float64)
    for p in pr:
        p = int(p)
        AN[p::p] /= (1.0 - 1.0 / (p * (p - 1.0)))
        if p > 2:
            SN[p::p] *= (1.0 + 1.0 / (p - 2.0))

    say("convolving r~ = Lambda*Lambda, C = mu*Lambda, U = mu^2*Lambda ...")
    n = pow2(2 * (X + 1))
    a = np.zeros(n, dtype=np.float64)
    a[:X + 1] = lam
    FL = np.fft.rfft(a)
    rt = np.fft.irfft(FL * FL, n)[:X + 1]
    a[:] = 0.0
    a[:X + 1] = mu
    C = np.fft.irfft(FL * np.fft.rfft(a), n)[:X + 1]
    a[:] = 0.0
    a[:X + 1] = (mu != 0)
    U = np.fft.irfft(FL * np.fft.rfft(a), n)[:X + 1]
    del a, FL

    ev = np.arange(2, X + 1, 2, dtype=np.int64)
    margin = SN[ev] * (1.0 - AN[ev])

    # ------------------------------------------------------------- F1/F2
    j = int(np.argmin(margin))
    Nmin = int(ev[j])
    f1 = Nmin == PRIMORIAL
    f2 = 0.05 <= margin[j] <= 0.08
    say()
    say("F1/F2   where the margin S(N)(1-A(N)) is smallest")
    say("=" * 70)
    say("  argmin N              = %d   (published-free; primorial is %d)"
        % (Nmin, PRIMORIAL))
    say("  1 - A(N) there        = %.6f" % (1.0 - AN[Nmin]))
    say("  S(N) there            = %.6f" % SN[Nmin])
    say("  margin there          = %.6f   (band [0.05, 0.08])"
        % margin[j])
    say("  for comparison, N = 2 mod 4: 1-A = %.6f, S = %.6f, "
        "margin = %.6f"
        % (1.0 - AN[15_999_998], SN[15_999_998],
           SN[15_999_998] * (1.0 - AN[15_999_998])))
    say("  F1 %s   F2 %s" % ("hold" if f1 else "REFUTED",
                             "hold" if f2 else "REFUTED"))
    say("  the ten smallest margins in the range:")
    order = np.argsort(margin)[:10]
    for t in order:
        Nv = int(ev[t])
        say("    N = %-10d margin %.6f   omega_small = %d"
            % (Nv, margin[t],
               sum(1 for p in (2, 3, 5, 7, 11, 13, 17, 19, 23)
                   if Nv % p == 0)))

    # ---------------------------------------------------------------- F3
    L = np.log(ev.astype(np.float64))
    prod = margin * L * np.log(L)
    f3 = float(prod.min()) > 1.0
    say()
    say("F3   margin * log N * log log N")
    say("=" * 70)
    say("  min over even N <= %d : %.6f   at N = %d   (floor 1)"
        % (X, float(prod.min()), int(ev[int(np.argmin(prod))])))
    say("  same at the primorial : %.6f" % float(prod[j]))
    say("  F3 %s" % ("hold" if f3 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc). log log N is negative for N < e^e =")
    say("  15.15, so the statistic is not defined at the bottom of the")
    say("  field I pre-registered. Restricted to the field the rest of")
    say("  this work uses:")
    for lo in (16, 1_000, 100_000):
        m = ev >= lo
        pm = prod[m]
        say("    N >= %-8d : min %.6f at N = %d"
            % (lo, float(pm.min()), int(ev[m][int(np.argmin(pm))])))

    # ---------------------------------------------------------------- F4
    med = float(np.median(margin))
    f4 = med >= 0.2
    say()
    say("F4   the margin for a typical N")
    say("=" * 70)
    for q in (0.001, 0.01, 0.1, 0.5, 0.9, 0.99):
        say("  quantile %-6.3f : %.6f" % (q, float(np.quantile(margin, q))))
    say("  F4  median %.6f >= 0.2   %s" % (med, "hold" if f4 else "REFUTED"))

    # ---------------------------------------------------------------- F5
    top = ev[ev > X // 2]
    ratio = U[top] / (AN[top] * top)
    viol = int(np.count_nonzero(np.abs(C[ev]) > U[ev] + 1e-6))
    f5 = 0.99 <= float(ratio.mean()) <= 1.01 and viol == 0
    say()
    say("F5   the triangle bound and its evaluation")
    say("=" * 70)
    say("  mean U(N)/(A(N)N) on the top octave = %.6f   (band [0.99,1.01])"
        % float(ratio.mean()))
    say("  sd                                  = %.6f" % float(ratio.std()))
    say("  even N with |C(N)| > U(N)           = %d" % viol)
    say("  F5 %s" % ("hold" if f5 else "REFUTED"))

    # ---------------------------------------------------------------- F6
    sel = ev[ev >= LOW]
    lhs = rt[sel]
    rhs = SN[sel] * (1.0 - AN[sel]) * sel
    bad = int(np.count_nonzero(lhs < rhs))
    slack = lhs / rhs
    f6 = bad == 0 and 3.5 <= float(slack.min()) <= 4.5
    say()
    say("F6   the route's conclusion against the truth")
    say("=" * 70)
    say("  even N in [%d, %d] with r~ < S(1-A)N : %d" % (LOW, X, bad))
    say("  slack r~ / (S(1-A)N):  min %.4f  median %.4f  max %.4f"
        % (float(slack.min()), float(np.median(slack)),
           float(slack.max())))
    say("  argmax slack N = %d  (the margin is thinnest exactly where the"
        % int(sel[int(np.argmax(slack))]))
    say("  count is largest)")
    say("  F6 %s" % ("hold" if f6 else "REFUTED"))

    say()
    say("=" * 70)
    ok = f1 and f2 and f3 and f4 and f5 and f6
    say("F1 %s  F2 %s  F3 %s  F4 %s  F5 %s  F6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (f1, f2, f3, f4, f5, f6)))
    say("the one-sided threshold of (*) is of size S(N)(1-A(N))N, which "
        "is" if ok else "REFUTED")
    if ok:
        say("of order 1 for most N and never below %.4f/(log N log log N) "
            "here" % float(prod.min()))

    head = [
        "STATISTIC: the one-sided Goldbach margin S(N)(1-A(N)) at every",
        "           even N, with A(N) = prod_{q not| N}(1-1/(q(q-1))) and",
        "           S(N) the singular series; its argmin, quantiles, and",
        "           its product with log N log log N; the triangle-bound",
        "           quantity U(N) = sum_{n<N} Lambda(n) mu^2(N-n) against",
        "           A(N)N and against |C(N)|; and r~(N) = Lambda*Lambda",
        "           against S(N)(1-A(N))N with the slack between them.",
        "FIELD: even N up to 1.6e7 for the margin, [1e5, 1.6e7] for the",
        "       slack, and the top octave (8e6, 1.6e7] for U/(A N);",
        "       Lambda, mu and the squarefree indicator from an integer",
        "       sieve to 1.6e7; A and S built by strided Euler products;",
        "       r~, C and U by exact FFT convolution.",
        'NULL: none applies. S(N)(1-A(N)) is deterministic arithmetic with no',
        '      sign input, and r~(N) vs S(1-A)N is a deterministic',
        '      comparison; there is no sign pattern to permute.',
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

# -*- coding: utf-8 -*-
r"""
paper/theorem_A.md, Proposition {#prop:E} and its table -- the circle
method's margin on C(N).

WHAT IS UNDER TEST

Proposition [prop:E] says both standard estimates of

    C(N) = int_0^1 S_Lambda(a) S_mu(a) e(-Na) da

The pairing is S_mu(+a), not S_mu(-a).  C(N) = sum_{n+m=N} is a
convolution, and S_mu(-a) would make it a correlation: the surviving
frequency would be n - m - N = 0, which no pair with n <= N and
m >= 1 can meet, so that integral is identically zero.  An earlier
printing of this docstring and of [prop:E]'s display carried the
conjugate; line 72 below and the code at run() always had the
convolution.  Neither number in this file moves -- (i) and (ii) are
built from norms and a supremum, and both are unchanged by a -> -a.

sit at or above the trivial bound psi(N) ~ N:

  (i)  ||S_Lambda||_2 ||S_mu||_2 ~ (6/pi^2)^{1/2} N (log N)^{1/2},
       above the trivial bound by a factor that grows;
  (ii) sup_a |S_mu| * ||S_Lambda||_1 >= ||S_mu||_2 ||S_Lambda||_1 >> N,
       the first factor floored by Parseval and the second by Vaughan.

The paper prints a table at N = 2^14, 2^16, 2^18, 2^20 with rows
||S_mu||_2/sqrt N, sup|S_mu|/sqrt N, ||S_Lambda||_1/sqrt N, the route-(i)
bound over N, the route-(ii) margin N/(sup|S_mu| ||S_Lambda||_1), and
C(N)/N.  No script for it exists in this repository.  This one
recomputes the table from scratch.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  R1  ||S_mu||_2/sqrt N reproduces 0.7798, 0.7797, 0.7797, 0.7797 --
      i.e. sqrt(6/pi^2) = 0.779697, which is a pure calibration row and
      must come out exactly.
  R2  The route-(ii) margin row reproduces 0.168, 0.175, 0.158, 0.152 to
      the printed precision, and is below 1 at every N.
  R3  The route-(i) row reproduces 2.297, 2.473, 2.639, 2.795, and the
      two l2 norms it is built from come out the same whether they are
      summed over the coefficient arrays or over the 4N-point grid --
      Parseval, which the grid must satisfy if the padding is enough
      to hold the convolution.  (The check written here first compared
      the row against a quantity defined as the row itself: the
      residual was zero as algebra, not as measurement, and it was not
      wired into R3 either.  Both are fixed; the reading it was after
      is the Parseval one.)
  R4  C(N)/N reproduces -0.0105, 0.0001, 0.0059, 0.0032.
  R5  ||S_Lambda||_1/sqrt N reproduces 1.946, 2.084, 2.219, 2.346, and
      sup|S_mu|/sqrt N reproduces 3.058, 2.742, 2.853, 2.801.  These two
      are grid-dependent (the paper says a 4N-point grid) and are the
      rows most likely to move under a different convention.

REFUTATION RULE (fixed before the run)

  R1  REFUTED if any entry differs from the printed value by more than
      0.0002.
  R2  REFUTED if any entry differs by more than 0.002, or is >= 1.
  R3  REFUTED if any entry differs from the printed value by more than
      0.002, or if the grid and the exact norms differ by more than
      1e-12 relative at any N.
  R4  REFUTED if any entry differs by more than 0.0002.
  R5  REFUTED if any entry differs by more than 0.01.  A refutation here
      is a statement about the grid convention, not about
      Proposition [prop:E]; it is reported separately for that reason.

  The script exits non-zero if any of R1-R4 is refuted.  R5 is reported
  but does not gate, because the proposition does not depend on it: its
  content is that the margin is below 1, and R2 carries that.

WHAT THE PARSEVAL CLAUSE OF R3 MEASURES, AND WHAT IT DOES NOT

  sum_j |X_j|^2 = L sum_n |x_n|^2 holds for the L-point DFT of any
  zero-padded array, for every L that holds the array at all.  So the
  clause tests that the grid values are the transform of the
  coefficient arrays it says they are -- a wiring check, and a live
  one -- but it says nothing about padding: no L >= N+1 fails it.
  Nor can anything else in this script.  The grid integral
  (1/L) sum_j S_Lambda(j/L) S_mu(j/L) e(-Nj/L) equals the direct sum
  sum_{n<N} Lambda(n) mu(N-n) unless some pair n + m = N + L with
  n, m <= N exists, which needs L <= N, and a grid that short cannot
  hold the coefficient arrays at all.  So for the one coefficient the
  table prints, "the padding holds the convolution" is not a
  measurement but the placement of the arrays; the 4N is oversampling
  for the sup and the L^1 rows, not padding for C(N).  What the grid
  integral does test, and the norm clause does not, is the wiring
  from the transforms through their product to the index N: it is
  computed below as a diagnostic that does not gate, with residual
  |grid - direct| / N.  A coefficient array placed one index off, or
  a conjugated factor, moves it to order 1e-2 and leaves every norm
  where it was.

CONVENTIONS FIXED BEFORE THE RUN

  S_Lambda and S_mu are evaluated on the 4N-point grid a = j/(4N),
  j = 0..4N-1, by zero-padded FFT of the coefficient arrays on [1,N].
  ||.||_1 is the grid mean of the modulus (a Riemann sum for the true
  L^1 norm); sup is the grid maximum, hence a lower bound for the true
  supremum -- which is the conservative direction for route (ii), since
  a larger sup only makes the margin smaller.
  ||.||_2 is taken from Parseval, i.e. exactly, not from the grid.
"""

import io
import math
import os
import sys

import numpy as np
from _sieve_shared import sieves as _shared_sieves
from _sieve_shared import spf_upto as _shared_spf

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "audit_circle_margin.txt")

NS = [2 ** 14, 2 ** 16, 2 ** 18, 2 ** 20]
IDTOL = 1e-12                 # Parseval on the grid against the exact norms
CTOL = 1e-12                  # C(N) by the grid integral against the direct sum, over N

PUB = {                       # the table as printed in paper/theorem_A.md
    "mu2":   [0.7798, 0.7797, 0.7797, 0.7797],
    "supmu": [3.058, 2.742, 2.853, 2.801],
    "lam1":  [1.946, 2.084, 2.219, 2.346],
    "route1": [2.297, 2.473, 2.639, 2.795],
    "margin": [0.168, 0.175, 0.158, 0.152],
    "CN":    [-0.0105, 0.0001, 0.0059, 0.0032],
}


def sieves(n):
    """정본 체에 위임한다 -- lib/goldbach/sieve.py 하나가 몸통이다."""
    _, lam, mu = _shared_sieves(n)
    return mu, lam


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    mu, lam = sieves(NMAX)

    got = {k: [] for k in PUB}
    gids = []
    cids = []
    say()
    say("  N        ||Smu||2/vN  sup|Smu|/vN  ||SL||1/vN   (i)/N     "
        "margin    C(N)/N")
    say("  " + "-" * 76)
    for N in NS:
        rN = math.sqrt(N)
        L = 4 * N
        a = np.zeros(L, dtype=np.float64)
        a[1:N + 1] = lam[1:N + 1]
        SL = np.fft.fft(a)
        b = np.zeros(L, dtype=np.float64)
        b[1:N + 1] = mu[1:N + 1]
        Smu = np.fft.fft(b)

        mu2 = math.sqrt(float((mu[1:N + 1] ** 2).sum())) / rN
        lam2 = math.sqrt(float((lam[1:N + 1] ** 2).sum())) / rN
        supmu = float(np.abs(Smu).max()) / rN
        lam1 = float(np.abs(SL).mean()) / rN
        route1 = mu2 * lam2
        margin = 1.0 / (supmu * lam1)

        # 격자에서 다시 잰 같은 두 노름. Parseval 은 4N 점 격자 위에서
        # 정확하므로 이 둘은 계수 배열의 값과 부동소수점 오차 안에서
        # 같아야 한다. 두 경로가 실제로 다른 배열을 훑는다는 것이 이
        # 검사가 검사인 이유다 -- route1 을 mu2*lam2 로 정의해 놓고
        # 그 곱을 다시 재면 무엇을 넣어도 0 이 나온다.
        gmu2 = math.sqrt(float((np.abs(Smu) ** 2).sum()) / L) / rN
        glam2 = math.sqrt(float((np.abs(SL) ** 2).sum()) / L) / rN
        gid = max(abs(gmu2 - mu2) / mu2, abs(glam2 - lam2) / lam2)

        n = np.arange(1, N, dtype=np.int64)
        C = float((lam[1:N] * mu[N - n]).sum()) / N

        # C(N) once more, from the grid: (1/L) sum_j S_L(j/L) S_mu(j/L)
        # e(-Nj/L), which is the circular convolution of the two padded
        # arrays read at index N.  Any grid that holds the arrays holds
        # this coefficient without aliasing, so what the comparison
        # tests is the wiring from the transforms to C(N), which the
        # norm identity never passes through.
        Cg = float(np.fft.ifft(SL * Smu)[N].real) / N
        cid = abs(Cg - C)

        for k, v in (("mu2", mu2), ("supmu", supmu), ("lam1", lam1),
                     ("route1", route1), ("margin", margin), ("CN", C)):
            got[k].append(v)
        gids.append(gid)
        cids.append(cid)
        say("  %-8d %-12.6f %-12.4f %-12.4f %-9.4f %-9.4f %+.6f"
            % (N, mu2, supmu, lam1, route1, margin, C))

    say()
    say("  sqrt(6/pi^2) = %.6f" % math.sqrt(6.0 / math.pi ** 2))
    say()
    say("published vs recomputed")
    say("  " + "-" * 60)
    tol = {"mu2": 0.0002, "margin": 0.002, "route1": 0.002,
           "CN": 0.0002, "supmu": 0.01, "lam1": 0.01}
    res = {}
    for k in ("mu2", "margin", "route1", "CN", "supmu", "lam1"):
        d = max(abs(g - p) for g, p in zip(got[k], PUB[k]))
        res[k] = d <= tol[k]
        say("  %-8s max |published - recomputed| = %.6f   tol %.4f   %s"
            % (k, d, tol[k], "ok" if res[k] else "MISMATCH"))
    ident = max(gids)
    say("  Parseval on the grid against the exact norms: worst relative "
        "%.2e   tol %.1e   %s"
        % (ident, IDTOL, "ok" if ident <= IDTOL else "MISMATCH"))
    below1 = all(m < 1 for m in got["margin"])
    say("  route (ii) margin below 1 at every N: %s" % below1)
    cworst = max(cids)
    say("  DIAGNOSTIC (post hoc, does not gate). C(N) by the grid")
    say("  integral against the direct sum, worst |grid - direct|/N: "
        "%.2e   tol %.1e   %s"
        % (cworst, CTOL, "ok" if cworst <= CTOL else "MISMATCH"))
    say("  The norm identity above holds on any grid that holds the")
    say("  arrays, and so does this coefficient; what this line adds is")
    say("  the wiring from the transforms through their product to C(N).")

    R1, R4 = res["mu2"], res["CN"]
    R3 = res["route1"] and ident <= IDTOL
    R2 = res["margin"] and below1
    R5 = res["supmu"] and res["lam1"]
    say()
    say("=" * 68)
    say("R1 %s  R2 %s  R3 %s  R4 %s   (R5 grid convention: %s)"
        % tuple(("hold" if v else "REFUTED")
                for v in (R1, R2, R3, R4, R5)))
    gating = R1 and R2 and R3 and R4
    say("Proposition {#prop:E}'s table reproduces" if gating else
        "REFUTED -- the published table does not reproduce")

    head = [
        "STATISTIC: ||S_mu||_2/sqrt(N) and ||S_Lambda||_2/sqrt(N) by",
        "           Parseval (exact); sup|S_mu|/sqrt(N) as the maximum and",
        "           ||S_Lambda||_1/sqrt(N) as the mean modulus over a",
        "           4N-point grid; the route-(i) bound ||S_L||_2||S_mu||_2/N;",
        "           the route-(ii) margin N/(sup|S_mu| ||S_Lambda||_1); and",
        "           C(N)/N = sum_{n<N} Lambda(n) mu(N-n) / N.",
        "DENOM: the Parseval check is relative to the exact l2 norms",
        "       of the coefficient arrays; the grid integral of C(N) is",
        "       compared with the direct sum over N;",
        "FIELD: N = 2^14, 2^16, 2^18, 2^20; coefficients on [1,N]; grid",
        "       a = j/(4N), j = 0..4N-1, by zero-padded FFT; Lambda and mu",
        "       from an integer sieve to 2^20.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not gating:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())

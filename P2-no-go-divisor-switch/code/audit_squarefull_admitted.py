# -*- coding: utf-8 -*-
r"""
paper/theorem_A.md, the note under Proposition {#prop:flatsum}
("the restriction is not decoration").

WHAT IS UNDER TEST

[eq:flatsum] carries the restriction mu^2(k), i.e. it sums over
squarefree k only:

    T_w = sum_{k<K,(k,N)=1} mu^2(k) w_k H(N;k) - B_w C(N),
    H(N;k) = sum_{m<N/k,(m,k)=1} Lambda(N-mk) mu(m).

The note says dropping mu^2(k) does not weaken the identity but
falsifies it, because H(N;k) does not vanish at squarefull k while the
left-hand side does, and it puts a size on the terms so admitted:

    at theta' = 0.56 and N = 2e5 they add 5.4e-2 N at w = 1 and
    2.9e-1 N at w = log, the latter being two thirds of |T_log|;
    over N = 2,4,8e5 that fraction reads 0.656, 0.702, 0.789.

Those five figures had no run behind them.  The note cited no file, the
reproducibility table had no row for it, and the values were not in the
packet -- the sixth blind pass on P2 searched all 84 result files and
all 85 scripts and found none of them in this context.  Three of the
five had been carried into the paper from a reviewer's own scratch
computation, which is not a packet artefact and cannot be re-run from
the deployment.

This script is that run.  The quantity is

    admitted(w) = sum_{k<K, (k,N)=1, mu(k)=0} w_k H(N;k),

the mass the restriction excludes, computed from the definitions above
with an independent sieve.

PRE-REGISTERED PREDICTIONS (the paper's own printed figures, registered
before this script was run)

  V1  admitted(1)/N   = 0.054 at N = 2e5, to two significant figures.
  V2  admitted(log)/N = 0.29  at N = 2e5, to two significant figures.
  V3  admitted(log)/|T_log| reads 0.656, 0.702, 0.789 at
      N = 2, 4, 8e5, to three decimals.
  V4  That fraction is increasing in N over those three points.
  V5  The admitted mass is not negligible: admitted(log)/|T_log|
      exceeds 1/2 at every N tested.  This is the note's actual claim --
      that the restriction is load-bearing rather than cosmetic -- and
      it is weaker than V1-V4.

REFUTATION RULE (fixed before the run)

  V1 REFUTED if |admitted(1)/N - 0.054| > 0.0005 at N = 2e5.
  V2 REFUTED if |admitted(log)/N - 0.29| > 0.005 at N = 2e5.
  V3 REFUTED by a single N at which the printed three-decimal value is
     not reproduced to within 0.0005.
  V4 REFUTED if the fraction is not strictly increasing over the three.
  V5 REFUTED by a single N at which the fraction is at most 1/2.

  V5 gates: the script exits non-zero if it fails, because it is what
  the note asserts.  V1-V4 are the printed figures and are reported as
  findings about them.

BACKS: the note under Proposition {#prop:flatsum} in paper/theorem_A.md.
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
OUT = os.path.join(ROOT, "results", "audit_squarefull_admitted.txt")

THETA = 0.56
FIELD = [200000, 400000, 800000]
PUB_ADM1 = 0.054
PUB_ADMLOG = 0.29
PUB_FRAC = {200000: 0.656, 400000: 0.702, 800000: 0.789}


def sieve(n):
    """정본 체에 위임한다 -- lib/goldbach/sieve.py 하나가 몸통이다."""
    _, lam, mu = _shared_sieves(n)
    spf = _shared_spf(n, np.int64)
    return spf, mu, lam


def _pure(w, p, spf):
    while w > 1:
        if spf[w] != p:
            return False
        w //= spf[w]
    return True


def coprime(x, primes):
    return all(x % p for p in primes)


def factors(n, spf):
    out = set()
    while n > 1:
        out.add(spf[n])
        n //= spf[n]
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("STATISTIC: admitted(w) = sum over k < K with (k,N)=1 and")
    say("           mu(k) = 0 of w_k H(N;k), the mass the mu^2(k) of")
    say("           [eq:flatsum] excludes; and its size against N and")
    say("           against |T_log| from the same identity.")
    say("FIELD: N = 2e5, 4e5, 8e5; theta' = 0.56; K = floor(N^theta').")
    say("NULL: the paper's own printed figures, registered as V1-V4.")
    say("DENOM: admitted/N is over N; admitted/|T_log| is over the")
    say("       |T_log| this same run computes, not a published one.")
    say()

    ok = {}
    rows = []
    nmax = max(FIELD)
    spf, mu, lam = sieve(nmax)

    for N in FIELD:
        K = int(N ** THETA)
        pN = factors(N, spf)
        adm1 = admlog = 0.0
        sf1 = sflog = 0.0
        C = sum(lam[n] * mu[N - n] for n in range(2, N) if lam[n])
        B1 = Blog = 0.0
        for k in range(1, K):
            if not coprime(k, pN):
                continue
            H = 0.0
            for m in range(1, N // k + 1):
                if math.gcd(m, k) != 1:
                    continue
                r = N - m * k
                if r >= 2 and lam[r] and mu[m]:
                    H += lam[r] * mu[m]
            if mu[k] == 0:
                adm1 += H
                admlog += math.log(k) * H
            else:
                sf1 += H
                sflog += math.log(k) * H
                ph = _phi(k, spf)
                B1 += mu[k] / ph
                Blog += mu[k] * math.log(k) / ph
        T1 = sf1 - B1 * C
        Tlog = sflog - Blog * C
        frac = abs(admlog) / abs(Tlog) if Tlog else float("nan")
        rows.append((N, K, adm1 / N, admlog / N, abs(Tlog) / N, frac))
        say("  N = %-8d K = %-6d admitted(1)/N = %+.6f  "
            "admitted(log)/N = %+.6f" % (N, K, adm1 / N, admlog / N))
        say("             |T_log|/N = %.6f   admitted(log)/|T_log| = %.6f"
            % (abs(Tlog) / N, frac))

    say()
    say("=" * 70)
    a1 = rows[0][2]
    alog = rows[0][3]
    ok["V1"] = abs(abs(a1) - PUB_ADM1) <= 0.0005
    ok["V2"] = abs(abs(alog) - PUB_ADMLOG) <= 0.005
    ok["V3"] = all(abs(r[5] - PUB_FRAC[r[0]]) <= 0.0005 for r in rows)
    ok["V4"] = rows[0][5] < rows[1][5] < rows[2][5]
    ok["V5"] = all(r[5] > 0.5 for r in rows)
    say("  V1 published %.3f   measured %.6f   %s"
        % (PUB_ADM1, abs(a1), "hold" if ok["V1"] else "REFUTED"))
    say("  V2 published %.2f    measured %.6f   %s"
        % (PUB_ADMLOG, abs(alog), "hold" if ok["V2"] else "REFUTED"))
    for r in rows:
        say("  V3 N = %-8d published %.3f   measured %.6f"
            % (r[0], PUB_FRAC[r[0]], r[5]))
    say("  V3 %s" % ("hold" if ok["V3"] else "REFUTED"))
    say("  V4 increasing in N: %s" % ("hold" if ok["V4"] else "REFUTED"))
    say("  V5 fraction above 1/2 at every N: %s"
        % ("hold" if ok["V5"] else "REFUTED"))
    say()
    say("  " + "  ".join("%s %s" % (k, "hold" if v else "REFUTED")
                         for k, v in sorted(ok.items())))
    verdict = all(ok.values())
    say("REFUTED" if not verdict else "all hold")
    say()
    if not ok["V5"]:
        say("The exit is non-zero because V5 is refuted: the mass the")
        say("restriction excludes is not load-bearing at some N, which")
        say("is what the note asserts.")
    elif not verdict:
        say("The exit is non-zero because at least one of the paper's")
        say("printed figures is not reproduced.  V5, which is the note's")
        say("assertion rather than its arithmetic, holds: the excluded")
        say("mass is a constant fraction of |T_log| at every N here, so")
        say("the restriction is load-bearing.  What is refuted is the")
        say("printing, and the note is corrected to what this run reads.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(lines) + "\n")
    return 0 if verdict else 1


def _phi(n, spf):
    r = n
    m = n
    seen = set()
    while m > 1:
        p = spf[m]
        if p not in seen:
            seen.add(p)
            r -= r // p
        m //= p
    return r


if __name__ == "__main__":
    sys.exit(main())

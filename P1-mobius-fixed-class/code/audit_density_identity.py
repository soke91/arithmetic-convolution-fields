# -*- coding: utf-8 -*-
r"""
paper/theorem_A.md, Lemma {#lem:density} and Remark {#rem:loadbearing}.

WHAT IS UNDER TEST

Lemma [lem:density] evaluates, for squarefree m with (m,N)=1,

    c(m) = sum_{d>=1} sum_{e|m} mu(d) mu(e) 1_{(d,N)=1}
             / phi( m * lcm(d^2, e) )

and asserts c(m) = A(N) * lambda(m) / m, with
A(N) = prod_{p not| N, p>2} (1 - 1/(p(p-1))) and
lambda(m) = prod_{p|m} (1 - 1/(p(p-1)))^{-1}.

Remark [rem:loadbearing] restates the m-local half of that computation
as the identity

    sum_{g|m} mu(g) / ( phi(m/g) * g * phi(g) ) = 1/m        (squarefree m)

and states that the exponent 1 here is load-bearing: a non-integral
density exponent would give, by Selberg-Delange, only (log x)^{-c} in
Lemma [lem:mu], and Theorem [thm:A] would be false. So this is the one
line in the proof whose failure is not a weakening but a refutation.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Q1  The identity of Remark [rem:loadbearing] holds exactly, in
      rational arithmetic, for every squarefree m < 400.
  Q2  The number of squarefree m in [1,400) is 243, matching the count
      the remark reports.  (An independent handle on whether the paper
      and this script are enumerating the same set.)
  Q3  The truncated c_{D,E}(m) of Lemma [lem:density] approaches
      A(N)*lambda(m)/m as D grows, with the error consistent with the
      lemma's stated O(1/(m*D)).
  Q4  The per-prime local factor of the p|m case is exactly 1/p: summed
      from the same definition Q1 uses, over the divisors of p, it is
      1/(p-1) - 1/(p(p-1)) = 1/p, exactly in rationals for every prime
      p < 200.  (Written out by hand here it carried two further terms,
      +-1/(p^2(p-1)), which are the same rational with opposite signs;
      they cancelled, and what was left was true of any local factor at
      all.  The sum is now taken from the definition.)
      Q4 IS NOT INDEPENDENT OF Q1.  Every prime p < 200 is a squarefree
      m < 400, and the expression, the arithmetic and the tolerance are
      Q1's, so Q4 cannot fail unless Q1 does: it is Q1 restricted to
      the primes.  The repair above turned a check that judged nothing
      into a check that judges nothing new.  It is kept, and named as
      redundant, because deleting it would delete the record of what
      it was.  Three of these four gates are independent.

REFUTATION RULE (fixed before the run)

  Q1 REFUTED by a single squarefree m < 400 with a nonzero residual.
  Q2 REFUTED if the count is not 243.
  Q3 REFUTED if |c_{D,E}(m) - A(N)lambda(m)/m| * m * D exceeds 10 at the
     largest D for any m tested (i.e. the O(1/(mD)) rate is not seen).
  Q4 REFUTED by a single prime with a nonzero residual.
  Non-zero exit on any refutation.
"""

import io
import math
import os
import sys
from fractions import Fraction

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "audit_density_identity.txt")


def primes_upto(n):
    s = [True] * (n + 1)
    s[0] = s[1] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            for q in range(p * p, n + 1, p):
                s[q] = False
    return [i for i, b in enumerate(s) if b]


def factor(n, ps):
    out = []
    for p in ps:
        if p * p > n:
            break
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            out.append((p, e))
    if n > 1:
        out.append((n, 1))
    return out


def phi(n, ps):
    v = 1
    for p, e in factor(n, ps):
        v *= (p - 1) * p ** (e - 1)
    return v


def mobius(n, ps):
    f = factor(n, ps)
    if any(e > 1 for _, e in f):
        return 0
    return -1 if len(f) % 2 else 1


def divisors(n, ps):
    ds = [1]
    for p, e in factor(n, ps):
        ds = [d * p ** i for d in ds for i in range(e + 1)]
    return sorted(ds)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    ps = primes_upto(200_000)

    # ---------------------------------------------------------------- Q1, Q2
    say("Q1/Q2  the load-bearing identity, exact rationals")
    say("=" * 68)
    sqfree = [m for m in range(1, 400) if mobius(m, ps) != 0]
    bad = []
    for m in sqfree:
        tot = Fraction(0)
        for g in divisors(m, ps):
            tot += Fraction(mobius(g, ps),
                            phi(m // g, ps) * g * phi(g, ps))
        if tot != Fraction(1, m):
            bad.append((m, tot))
    say("  squarefree m in [1,400):     %d   (remark says 243)" % len(sqfree))
    say("  mismatches:                  %d" % len(bad))
    if bad:
        for m, t in bad[:10]:
            say("    m=%d  got %s  want %s" % (m, t, Fraction(1, m)))
    q1 = not bad
    q2 = len(sqfree) == 243
    say("  -> Q1 %s   Q2 %s" % ("HOLDS" if q1 else "REFUTED",
                                "HOLDS" if q2 else "REFUTED"))

    # ------------------------------------------------------------------- Q4
    say()
    say("Q4  the p|m local factor is exactly 1/p, exact rationals")
    say("  NOT AN INDEPENDENT GATE. Every prime p < 200 is a")
    say("  squarefree m < 400, and the expression and the tolerance")
    say("  are Q1's, so Q4's field is contained in Q1's and Q4")
    say("  cannot fail unless Q1 does. It is Q1 restricted to the")
    say("  primes, kept because the note below records what it was")
    say("  before -- a hand-written form whose last two terms")
    say("  cancelled, which judged nothing. The repair turned a dead")
    say("  check into a redundant one; three gates are independent,")
    say("  not four.")
    say("=" * 68)
    badp = []
    for p in primes_upto(200):
        # 국소인자를 손으로 적은 식이 아니라 Q1 이 쓰는 그 정의에서
        # 다시 만든다. 손으로 적은 식은 뒤 두 항이 같은 유리수의 부호
        # 반대여서 서로 지워졌고, 그러면 남는 것은 어떤 국소인자에
        # 대해서도 참인 항등식이라 무엇도 판정하지 않았다.
        v = Fraction(0)
        for g in divisors(p, ps):
            v += Fraction(mobius(g, ps), phi(p // g, ps) * g * phi(g, ps))
        if v != Fraction(1, p):
            badp.append((p, v))
    say("  primes tested:               %d" % len(primes_upto(200)))
    say("  mismatches:                  %d" % len(badp))
    q4 = not badp
    say("  -> Q4 %s" % ("HOLDS" if q4 else "REFUTED"))

    # ------------------------------------------------------------------- Q3
    say()
    say("Q3  truncated c_{D,E}(m) against A(N) lambda(m) / m")
    say("=" * 68)
    N = 10 ** 6                       # 2^6 * 5^6, so P(N) = {2,5}
    PN = {2, 5}
    PLIM = 4_000_000
    A = 1.0
    for p in primes_upto(PLIM):
        if p == 2 or p in PN:
            continue
        A *= 1.0 - 1.0 / (p * (p - 1.0))
    say("  N = %d   P(N) = %s   A(N) = %.8f" % (N, sorted(PN), A))
    say("  m     D        c_{D,E}(m)        A*lambda/m        |diff|*m*D")
    q3 = True
    for m in (1, 3, 7, 21, 143, 187):
        if any(p in PN for p in [q for q, _ in factor(m, ps)]):
            continue
        lam = 1.0
        for p, _ in factor(m, ps):
            lam /= (1.0 - 1.0 / (p * (p - 1.0)))
        target = A * lam / m
        dm = divisors(m, ps)
        for D in (100, 1000, 10000):
            tot = 0.0
            for d in range(1, D + 1):
                md = mobius(d, ps)
                if md == 0 or math.gcd(d, N) != 1:
                    continue
                d2 = d * d
                for e in dm:
                    me = mobius(e, ps)
                    q = m * (d2 * e // math.gcd(d2, e))
                    tot += md * me / phi(q, ps)
            scaled = abs(tot - target) * m * D
            say("  %-5d %-8d %.10f      %.10f      %.3f"
                % (m, D, tot, target, scaled))
            if D == 10000 and scaled > 10:
                q3 = False
    say("  -> Q3 %s" % ("HOLDS" if q3 else "REFUTED"))

    say()
    say("=" * 68)
    verdict = q1 and q2 and q3 and q4
    say("Q1 %s  Q2 %s  Q3 %s  Q4 %s"
        % tuple("hold" if v else "REFUTED" for v in (q1, q2, q3, q4)))
    say("Lemma {#lem:density} and Remark {#rem:loadbearing} stand"
        if verdict else "REFUTED -- Theorem {#thm:A} does not survive this")

    head = [
        "STATISTIC: (a) the residual of sum_{g|m} mu(g)/(phi(m/g) g phi(g))",
        "           - 1/m in exact rational arithmetic; (b) the count of",
        "           squarefree m in [1,400); (c) the residual of the p|m",
        "           local factor minus 1/p in exact rationals; (d) the",
        "           truncated c_{D,E}(m) of Lemma {#lem:density} against",
        "           A(N) lambda(m) / m, and that gap rescaled by m*D.",
        "FIELD: (a),(b) every squarefree m in [1,400); (c) every prime",
        "       p < 200; (d) N = 10^6, m in {1,3,7,21,143,187} (all coprime",
        "       to N), truncations D = 10^2, 10^3, 10^4 with e | m complete;",
        "       A(N) as an Euler product over p < 4*10^6.",
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

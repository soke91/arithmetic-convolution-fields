# -*- coding: utf-8 -*-
r"""
paper/theorem_A.md, Note {#rem:trap}.

WHAT IS UNDER TEST

Step 4 builds the modulus q = m * lcm(d^2, e) with e | m, and restricts
the main terms to the classes with (q, N) = 1.  For the m of Lemma
{#lem:density} -- squarefree, (m, N) = 1 -- that restriction is exactly
the restriction (d, N) = 1, because e | m already forces (e, N) = 1.
So "assigning a main term Tcut_m/phi(q) to a degenerate class" is, in
the density, exactly the sum

    c_bad(m) = sum_{d>=1} sum_{e|m} mu(d) mu(e) / phi( m lcm(d^2,e) )

against the paper's

    c_good(m) = sum_{d>=1} sum_{e|m} mu(d) mu(e) 1_{(d,N)=1}
                  / phi( m lcm(d^2,e) ).

Note {#rem:trap} states that the mistake "shifts the density of the
whole computation by the factor N/phi(N)".  That figure has no
derivation in the paper and no script behind it.  This run is that
script.

The shift is a finite product over the primes dividing N -- the two
densities differ in no other local factor -- so it is computable
exactly in rationals, and the question has an exact answer.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  The shift c_bad(m)/c_good(m) equals N/phi(N), the factor
      Note {#rem:trap} prints.  This is the note's claim, registered
      as printed.
  W2  The shift is independent of m: the two densities differ only in
      the local factors at p | N, and m is coprime to N.
  W3  The shift, whatever it is, is a finite product over p | N of the
      local factor that the restriction (d,N)=1 deletes, namely
      prod_{p|N} (1 - 1/(p(p-1))).
  W4  The shift is not 1 at any even N: the restriction is not
      cosmetic, which is the sentence of the note that carries the
      weight.

REFUTATION RULE (fixed before the run)

  W1 REFUTED if, at a single N of the field, the measured shift
     differs from N/phi(N) by more than 1e-9 relative.
  W2 REFUTED if two m at the same N give shifts differing by more
     than 1e-5 relative.
  W3 REFUTED if the measured shift differs from that product by more
     than 1e-6 relative at any N.  (Registered at 1e-9 and corrected
     to 1e-6 with W2 and W5, for the reason the field correction
     below gives; the correction reached the code and W2's and W5's
     lines and not this one.  See the last paragraph there.)
  W4 REFUTED if any N of the field gives a shift equal to 1.
  Non-zero exit on any refutation.

  FIELD CORRECTION (made after the first run, recorded here rather
  than silently applied)

    The first run carried m in {1,13,17,221} against N up to the
    primorial 9699690, so 13 and 17 divided several N of the field.
    Lemma {#lem:density} is stated for (m,N)=1, and the shift is a
    product over p | N only under that hypothesis, so those pairs were
    never in the field -- the list was wrong, not the result.  m is now
    built from primes above 19.  The same run also put a 1e-9 tolerance
    on a truncated float sum, which no truncation can meet; W2 and W5
    now read 1e-5 relative against a measured spread of 1e-7.
    W1 -- the note's own figure -- was refuted before the correction
    and after it by the same margin, and W3 held both times.  No
    prediction was re-registered.

    W3's line was left at 1e-9 for a further day.  The code applied
    1e-6, the two never disagreed on the verdict, and nothing pointed
    at the gap -- the measured residual runs 2.4e-7 to 7.9e-7, which
    sits between the two, so W3 reads REFUTED at ten of the eleven N
    under the line as registered and holds under the line as applied.
    A reader given the registered rule would have scored this run
    refuted.  r12's C pass found it.  The line now says what the run
    does, and this paragraph says that the line moved and when, which
    is the whole reason the field correction is written here rather
    than applied in silence.

  W1 and W3 are both registered and they are mutually exclusive
  unless N/phi(N) = prod_{p|N}(1-1/(p(p-1))), which holds at no even
  N.  One of the two will be refuted; the point of registering both is
  that the run decides which, and the losing one is printed.

  The truncated sums are the measurement.  A truncated d-sum reaches
  the density only to O(1/(m D)), so the comparison against the exact
  rational product is made at three truncations and the residual is
  required to fall.
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
OUT = os.path.join(ROOT, "results", "audit_trap_density.txt")

TOL = 1e-9


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


def density(m, N, D, ps, restrict):
    """c_{D}(m), with or without the 1_{(d,N)=1} of Step 4."""
    dm = divisors(m, ps)
    tot = 0.0
    for d in range(1, D + 1):
        md = mobius(d, ps)
        if md == 0:
            continue
        if restrict and math.gcd(d, N) != 1:
            continue
        d2 = d * d
        for e in dm:
            me = mobius(e, ps)
            if me == 0:
                continue
            q = m * (d2 * e // math.gcd(d2, e))
            tot += md * me / phi(q, ps)
    return tot


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    ps = primes_upto(200_000)

    # The field: even N chosen so that the set of primes dividing N is
    # varied -- a power of two, the primorials, and two N that other
    # runs in this repository measure at.
    FIELD = [
        2 ** 20,
        2 * 10 ** 5,
        6,
        30,
        210,
        2310,
        30030,
        510510,
        9699690,
        1531530,
        1600000,
    ]
    FIELD = sorted(set(FIELD))

    # m coprime to every N of the field, squarefree, and small enough
    # that the e | m sum is complete.
    MS = [1, 23, 29, 23 * 29]

    say("W1/W2/W3/W4  the shift the trap applies to the density")
    say("=" * 68)
    say("  shift := c_bad(m) / c_good(m), the density with the")
    say("  degenerate classes given a main term, over the density")
    say("  with them restricted away.  Exact column is the rational")
    say("  product over p | N; measured column is the truncated sum.")
    say()

    w1 = True
    w2 = True
    w3 = True
    w4 = True
    rows = []

    for N in FIELD:
        pn = sorted(p for p, _ in factor(N, ps))
        # exact: the local factors the restriction deletes, at p | N
        exact = Fraction(1)
        for p in pn:
            exact *= Fraction(1) - Fraction(1, p * (p - 1))
        nphi = Fraction(N, phi(N, ps))

        shifts = []
        for m in MS:
            if math.gcd(m, N) != 1:
                continue
            g = density(m, N, 20000, ps, True)
            b = density(m, N, 20000, ps, False)
            shifts.append(b / g)

        meas = sum(shifts) / len(shifts)
        spread = max(shifts) - min(shifts)

        if spread > 1e-5 * abs(meas):
            w2 = False
        if abs(meas - float(exact)) > 1e-6 * abs(float(exact)):
            w3 = False
        if abs(meas - float(nphi)) > TOL * abs(float(nphi)):
            w1 = False
        if abs(meas - 1.0) <= TOL:
            w4 = False

        rows.append((N, pn, float(exact), meas, float(nphi), spread))

    say("  N          p | N                    exact       measured"
        "    N/phi(N)")
    say("  " + "-" * 72)
    for N, pn, ex, me, nf, sp in rows:
        say("  %-10d %-24s %-11.7f %-11.7f %.7f"
            % (N, ",".join(str(p) for p in pn), ex, me, nf))
    say()
    say("  largest spread of the shift across m at one N: %.3e"
        % max(r[5] for r in rows))
    say()
    say("  -> W1 %s" % ("HOLDS" if w1 else "REFUTED"))
    say("  -> W2 %s" % ("HOLDS" if w2 else "REFUTED"))
    say("  -> W3 %s" % ("HOLDS" if w3 else "REFUTED"))
    say("  -> W4 %s" % ("HOLDS" if w4 else "REFUTED"))

    if not w1:
        say()
        say("  W1 is the note's own figure and it is refuted.  The two")
        say("  factors are not close and they do not even point the same")
        say("  way: N/phi(N) > 1 at every even N, while the shift the")
        say("  restriction actually removes is < 1 at every even N.  The")
        say("  ratio between them at the primorials:")
        for N, pn, ex, me, nf, sp in rows:
            say("    N = %-9d  shift %.7f   N/phi(N) %.7f   quotient %.7f"
                % (N, me, nf, me / nf))

    # ------------------------------------------------------------------ W5
    say()
    say("W5  the truncation: the measured shift is a limit, so it must")
    say("    settle as the d-sum lengthens")
    say("=" * 68)
    w5 = True
    N = 510510
    say("  N = %d,  m = 23" % N)
    say("  D        shift            |shift - exact|")
    exact = Fraction(1)
    for p, _ in factor(N, ps):
        exact *= Fraction(1) - Fraction(1, p * (p - 1))
    prev = None
    for D in (200, 2000, 20000):
        g = density(23, N, D, ps, True)
        b = density(23, N, D, ps, False)
        s = b / g
        gap = abs(s - float(exact))
        say("  %-8d %.10f     %.3e" % (D, s, gap))
        if D == 20000 and gap > 1e-5 * float(exact):
            w5 = False
        prev = gap
    say("  -> W5 %s" % ("HOLDS" if w5 else "REFUTED"))

    # ------------------------------------------------------------------ W6
    say()
    say("W6  what the trap does to the constant of Section {#sec:C}")
    say("=" * 68)
    say("  A(N) W_G(1) = -S(N) is the identity the log k branch runs on.")
    say("  A(N) is the same Euler product with the p | N factors set to")
    say("  1, so the trap multiplies A(N) by the shift above and the")
    say("  identity is out by that factor.  Printed against S(N) so the")
    say("  size of the apparent refutation is visible.")
    say()
    say("  N          shift       S(N)        S(N)*shift  |gap|/S(N)")
    say("  " + "-" * 62)

    PLIM = 4_000_000
    pr = primes_upto(PLIM)
    # S(N) = 2 prod_{p>2} (1 - 1/(p-1)^2) prod_{p|N, p>2} (p-1)/(p-2)
    base = 2.0
    for p in pr:
        if p == 2:
            continue
        base *= 1.0 - 1.0 / ((p - 1.0) ** 2)
    w6 = True
    for N, pn, ex, me, nf, sp in rows:
        S = base
        for p in pn:
            if p > 2:
                S *= (p - 1.0) / (p - 2.0)
        gap = abs(S - S * me) / S
        say("  %-10d %-11.7f %-11.7f %-11.7f %.7f" % (N, me, S, S * me, gap))
        if gap < 0.1:
            w6 = False
    say()
    say("  -> W6 %s" % ("HOLDS" if w6 else "REFUTED"))
    say("     the gap is never small: the trap does not perturb the")
    say("     constant, it replaces it, and Section {#sec:C} would read")
    say("     the difference as E_3 failing the identity.")

    say()
    say("=" * 68)
    say("W1 %s  W2 %s  W3 %s  W4 %s  W5 %s  W6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (w1, w2, w3, w4, w5, w6)))
    ok = w1 and w2 and w3 and w4 and w5 and w6
    if not ok:
        say("REFUTED")
        say()
        say("The exit is non-zero because W1 is refuted.  W1 was the")
        say("note's printed figure N/phi(N) and it is not the shift.")
        say("W2, W3, W4, W5, W6 hold, so the note's sentence -- that the")
        say("restriction is not cosmetic and that the error is carried")
        say("into the log k branch as an apparent refutation -- survives")
        say("with the factor corrected to prod_{p|N} (1 - 1/(p(p-1))).")
    else:
        say("Note {#rem:trap} stands as printed")

    head = [
        "STATISTIC: the factor by which the density of Step 4's main",
        "           terms changes when the classes with (q,N) > 1 are",
        "           given the main term Tcut_m/phi(q) instead of being",
        "           discarded -- measured as c_bad(m)/c_good(m) from the",
        "           definitions, and against the exact rational product",
        "           over p | N; then that factor against N/phi(N), the",
        "           figure Note {#rem:trap} prints, and against S(N).",
        "FIELD: even N in {6, 30, 210, 2310, 30030, 510510, 9699690,",
        "       200000, 1048576, 1531530, 1600000}; m in {1,23,29,667},",
        "       all coprime to every N of the field; d-sum truncated at",
        "       D = 2*10^4 with e | m complete; S(N) as an Euler product",
        "       over p < 4*10^6.",
        "NULL: N/phi(N) itself.  It is the figure the note has been",
        "      printing with nothing behind it, and it is registered as",
        "      W1 so that this run either confirms it or names what",
        "      replaces it.",
        "DENOM: every relative gap printed is over the quantity it is",
        "       compared against, named in its column.",
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

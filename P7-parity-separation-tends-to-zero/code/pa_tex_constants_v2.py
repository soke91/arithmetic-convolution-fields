# -*- coding: utf-8 -*-
r"""pa_tex_constants_v2 -- every constant printed by the manuscript, as of an earlier revision.

WHY A FOURTH CONSTANTS FILE

  pa_paper_constants.py      a 5-form draft's literals     -- 5-form family (c1 = 2.51647)
  pa_tex_constants.py        an early draft                -- 4-form, draft values
  pa_tex_constants_final.py  an earlier manuscript         -- 4-form, those values
  this file                  the manuscript at that time   -- 4-form, those values

The manuscript was afterwards revised: three literals were REMOVED, six CHANGED and a dozen
ADDED, so `pa_tex_constants_final.py` no longer matches the text it names.  That file is NOT edited:
retuning a registered file's literals after it has run destroys the provenance registration exists
for.  It stays correct about the earlier manuscript, and this file transcribes the text as it then
stood.  Its two failures are closed in the current text and its two display recommendations
applied:

  F13a  "Z_0 \approx 850" at L = 3767  ->  "\approx 852".  852.121226... is the value; the 850 had come
        from rounding it in prose to two significant figures.  Now V16a.
  F16   "8.3\cdot10^{5}" at 2^28       ->  "8.4\cdot10^{5}" (L^6 = 835066.34).  Now V19a.
  F9b   "0.4387\ldots" (3^-3/4)        ->  "0.4386\ldots" (0.4386913...).       Now V9b.
  F15a  "16.9989\ldots" (root 6logL=L) ->  "16.9988\ldots" (16.99888735...).    Now V18a.

REMOVED from the text, so no clause here may match them (clause X1 checks their absence):
  1.33258  (Rosser-Schoenfeld's additive constant in (M1): RS (3.24) gives < log x, so it is gone,
            and Xi is now the bare constant e^(e-1) rather than exp((e-1)(1 + 1.33258/log Y)))
  0.52, 0.8 (Platt-Trudgian's exponent and constant, which the draft had attributed to Trudgian)

THE RENDERING RULE, FIXED BY THE SOURCE BEFORE ANY VALUE IS COMPUTED (the reproduction rule)

Every clause compares a DISPLAY against a DISPLAY: the value is rendered to exactly the literal's own
decimal count and the two strings are compared.  Which rendering is read off the paper's notation:

    "X\ldots"                   -> X is a TRUNCATION  -> compare truncated
    "\approx X", "near X", ">= X" -> X is a ROUNDING  -> compare rounded
    a bare decimal the text calls a rounding or an upper bound
                                -> compare to the value rounded UP, and separately assert it exceeds

No clause is an `abs(x - literal) <= tol`, so no tolerance can be coarser or finer than its literal.
Clause D1 reports any literal matching a rendering OTHER than the one its notation declares -- a
one-character display bug, not a wrong constant -- and does not set the exit code.  D1 SKIPS a
literal whose truncation and rounding coincide: there is nothing to choose and nothing to fix.  (An
earlier version's D1 lacked that guard and reported `132.5` as a mismatch for that reason alone; the guard is the fix.)

CITED, not computed.  These are literature values; no clause recomputes them, and the clauses that
mention them test only the paper's own deduction from them:
  1.25506    Rosser-Schoenfeld (3.6)      -- V13 checks that the text's 1.26 is a valid rounding up
  0.2795, (log x)^(-3/4), 6.455, 229  Trudgian Thm 2 -- V14 checks 1/sqrt(6.455) and c_pi = 0.39
  every arXiv number and year in the bibliography

COMPUTE.  No cell is loaded and no array of length N is allocated.  The only array work is (a) trial
division of 8-digit integers by the canonical sieve's primes <= 8192 (`_sieve_shared.primes_upto`,
vectorised, no Python loop) for the 2^52 witness, and (b) a numpy construction of the band for
N <= 30 for the two level examples of section 1.1.  This machine, one thread, well under a second;
the compute-node rule does not engage.

CLAUSES.  Every literal below is transcribed from the .tex and is not edited after
any output.

  CLOSED FORMS, unchanged and still in the text
  V1  c_1 = (4/sqrt6) sqrt(15/pi^2)         "2.0131684\ldots"   tex:160
  V2  c_2 = 4 sqrt(2 e^(e-1))               "13.3565762\ldots"  tex:161
  V3  the headline decimal                  "2.01317"           tex:56,171,177 -- equals c_1 rounded
      UP to 5 decimals AND exceeds c_1 ("exceeds $c_1$, as an upper bound requires", tex:177)
  V4  pi^2/(3(1-log 2))                     "10.7213228\ldots"  tex:447
  V5  sqrt(zeta(2)/zeta(4)) = sqrt(15/pi^2) "1.2328088\ldots"   tex:534
  V6  Xi = e^(e-1)                          "5.5749415\ldots"   tex:559
  V7  the Buchstab remark, tex:660-662.  On [2,3] u varpi(u) = 1 + log(u-1), so the maximum is at
      u = 1+v with v log v = 1, where varpi = 1/v; hence
      sup varpi(u) e^gamma = e^gamma/v      "1.0101232\ldots"
      u* = 1+v                              "2.7632228\ldots"
      varpi(u*) = 1/v                       "0.5671432\ldots"   (ADDED)
  V8  1 - log 2                             "0.3068528\ldots"   tex:668
  V9  2^(-3/4)                              "0.5946\ldots"      tex:737, and
      3^(-3/4)                              "0.4386\ldots"      tex:738  (CHANGED)
  V10 (3/pi^2)(1-log 2)                     "0.0932720\ldots"   tex:1001, and it exceeds 0.09
      (the text's "R >= 0.09N", tex:446 and tex:1001)
  V11 4 sqrt(15/pi^2)/2                     "2.4656\ldots"      tex:1048 (the beta = 4 window)
  V12 4 beta e^gamma/(3(1-log 2)) at beta=6 "46.4345\ldots"     tex:837  (CHANGED)

  THE TWO EXPLICIT PRIME-COUNTING INPUTS, as the text deduces them (ADDED)
  V13 (C), tex:351-352.  The text uses 1.26 as "a rounding of the bound 1.25506" of RS (3.6):
      1.26 equals 1.25506 rounded UP to 2 decimals AND exceeds it, so the weaker form is implied.
  V14 (P), tex:365-370.  From Trudgian's Theorem 2 display the text deduces c_pi = 0.39:
      1/sqrt(6.455)                         "0.3935\ldots"      -> truncated, 4 decimals
      and 0.39 < 1/sqrt(6.455), so exp(-sqrt(log x/6.455)) <= exp(-0.39 sqrt(log x));
      and (log x)^(-3/4) <= 1 for x >= e, the text's "logarithmic factor is at most 1 for x >= e".

  THE CROSSING AND THE EFFECTIVITY SECTION
  V15 the crossing, tex:1069-1070: the unique root L* > 1000 of
      c_1/(L sqrt(log L)) + c_2 L^2 exp(-L/(24 log L)) = 1 renders
      L*                                    "near 3767"         -> rounded, 0 decimals
      2 L*/log 2                            "\approx 10868"     -> rounded, 0 decimals
  V16 at L = 3767 exactly, tex:1071-1072:
      Z_0 = exp(L/(log L)^3)                "\approx 852"       -> rounded, 0 decimals (CHANGED)
      Y = L^6                               "\approx 2.9e21"    -> mantissa rounded, 1 decimal
      u = log Z_0/log Y                     "\approx 0.14"      -> rounded, 2 decimals
  V17 tex:1073-1075.  u at z = Z_0 is L/(6 (log L)^4); the root of L/(6(log L)^4) = 1 renders
      L                                     "near 1.08e5"       -> mantissa rounded, 2 decimals
      2 L/log 2                             "\approx 3.1e5"     -> mantissa rounded, 1 decimal
      that L divided by 3767                "about 29"          -> rounded, 0 decimals
  V18 the two onsets, tex:1078-1079:
      root of 6 log L = L                   "16.9988\ldots"     -> truncated, 4 dp (CHANGED)
      2 L/log 2 there                       ">= 49.05"          -> rounded, 2 decimals
      2 L/log 2 at the root of 12 log L = L ">= 132.5"          -> rounded, 1 decimal
  V19 at N = 2^28, tex:1079-1080: Y = L^6 with L = log 2^14 renders "\approx 8.4e5" (CHANGED), and
      Y > Q there, so every band element is Y-smooth and x vanishes identically.

  THE N = 2^52 WITNESS (ADDED), tex:1083-1087.  All of it exact integer arithmetic except the
  two rendered reals.  The band of eq:band is {n in (thr,N] : mu^2(n)=1, P^+(n) <= Q}, split by the
  parity of omega into vB (even) and vR (odd); the small-prime part is a_Y(l) = prod_{p|l, p<=Y} p.
  V20 with N = 2^52:
      Q = floor(sqrt N) = 2^26              "67108864"          exact
      L = log Q                             "18.0218\ldots"     -> truncated, 4 decimals
      Y = L^6                               "34260433.5\ldots"  -> truncated, 1 decimal, and Y < Q
      q_1 = least prime > Q                 "67108879"          exact, and it is prime
      thr = floor(N/q_1)                    "67108849"          exact
      34260449 and 34260493 are both prime and both lie in (Y, Q]
      2*34260449 lies in vB   (in (thr,N], squarefree, P^+ <= Q, omega = 2 even)
      2*34260449*34260493 lies in vR  (in (thr,N], squarefree, P^+ <= Q, omega = 3 odd)
      both have small-prime part exactly 2, and 2 <= sqrt Q = 8192
      so B_2 > 0 and R_2 > 0 and eq:X carries positive weight at this N.

  THE LEVEL EXAMPLES OF SECTION 1.1
  V21 (ADDED) tex:1116.  The integers N with Q^2 <= N < (Q+1)^2 number 2Q+1.
  V22 (ADDED) tex:1117-1118.  At N = 10 and at N = 11 alike: Q = 3, thr = 2, band = {3,6}.
  V23 (unchanged) tex:1119-1121.  At Q = 5: the odd part of the band is {5} at N = 25 and {5,30} at
      N = 30, and at N = 30 the row vector t equals (A_10 + A_15)/2, which lies in K_Q; so
      kappa(25) > 0 = kappa(30).

  MEASURED -- the manuscript's printed digits against the stored result files
  M1  tex:1126.  The enclosures [0.1723893,0.1723894], [0.1408604,0.1408605], [0.1055426,0.1055428]
      at e = 16,18,20 each CONTAIN the stored exact-dual kappa, and each endpoint carries 7 decimals.
      Sources: `conedual_separator_mass.txt` for e = 16,18 and
      `conedual_exact_kappa_small_cells_20.txt` for e = 20.  Containment, not endpoint equality: the
      e = 20 lower endpoint is one unit looser than the floor of the stored value, which no clause
      here calls an error.
  M2  tex:1130.  At e = 22,24,26,28 the printed LOWER endpoint equals the certified integer c6/10^6
      and the UPPER endpoint equals the stored float upper, both rendered at 7 decimals, from
      `conedual_colgen_2e{22,24,26}_n0b.json` and `conedual_colgen_2e28_n28.json`.
  M3  tex:1136.  The construction's 0.6982, 0.6130, 0.6088, 0.5841 at e = 16,18,20,22 equal
      `best_kappa_ub` ROUNDED UP to 4 decimals (the text says "rounded upward"), from
      `conedual_parity_domination.json`, with `best_Y` in {3,5}; and the ratios to the enclosure
      upper endpoints all lie in [4,8] ("factors between four and eight", tex:1137).

  X1  ABSENCE.  The three removed literals "1.33258", "0.52" and "0.8" do not occur in the .tex.
      INCONCLUSIVE, not failed, when the .tex cannot be found (it is not part of the code packet).

  A clause fails only if a value matches NEITHER rendering of its literal, or an exact/integer claim
  is false.  Any failure sets a non-zero exit code.

    python code/pa_tex_constants_v2.py
"""
from __future__ import annotations

import io
import json
import os
import re
import sys
from decimal import Decimal, ROUND_CEILING, ROUND_DOWN, ROUND_HALF_UP

import mpmath as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _sieve_shared import primes_upto  # noqa: E402  (the canonical sieve)

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), "results")   # repo and packet alike: ../results
OUT = os.path.join(RES, "pa_tex_constants_v2.txt")
OUTJ = os.path.join(RES, "pa_tex_constants_v2.json")

mp.mp.dps = 40

_lines: list[str] = []
REC: list[dict] = []
BAD: list[str] = []
MISM: list[tuple] = []


def say(s: str = "") -> None:
    print(s)
    _lines.append(s)
    with io.open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(_lines) + "\n")


def _dec(x) -> Decimal:
    return Decimal(mp.nstr(mp.mpf(x), 35, strip_zeros=False))


def _q(x, n, how) -> str:
    return str(_dec(x).quantize(Decimal(1).scaleb(-n), rounding=how))


def trunc(x, n):
    return _q(x, n, ROUND_DOWN)


def rnd(x, n):
    return _q(x, n, ROUND_HALF_UP)


def rndup(x, n):
    return _q(x, n, ROUND_CEILING)


def mant(x, n: int) -> str:
    """mantissa in a.be+k form, rounded to n decimals -- how the .tex writes a.b\\cdot10^{k}."""
    k = int(mp.floor(mp.log10(mp.mpf(x))))
    return "%se%d" % (rnd(mp.mpf(x) / mp.mpf(10) ** k, n), k)


def clause(name, val, lit, nota, n, desc, mantissa=False):
    """display against display; `nota` fixes the rendering the paper declares."""
    if mantissa:
        t = r = mant(val, n)
    else:
        t, r = trunc(val, n), rnd(val, n)
    which = "trunc" if t == lit else ("round" if r == lit else "neither")
    want = "trunc" if nota == "ldots" else "round"
    ok = which != "neither"
    if not ok:
        BAD.append(name)
    # D1 only speaks when the two renderings actually differ
    if ok and not mantissa and t != r and which != want:
        MISM.append((name, lit, nota, which, t, r))
    REC.append(dict(name=name, form=desc, value=mp.nstr(mp.mpf(val), 25), literal=lit,
                    notation=nota, matches=which, holds=bool(ok)))
    say("%-5s %-8s %-30s computed %-14s literal %-14s [%s]"
        % (name, "HOLDS" if ok else "FAILS", desc,
           t if (ok and which == "trunc") else r, lit, which))


def plain(name, ok, detail):
    if not ok:
        BAD.append(name)
    REC.append(dict(name=name, holds=bool(ok), detail=detail))
    say("%-5s %-8s %s" % (name, "HOLDS" if ok else "FAILS", detail))


# ------------------------------------------------------------------ small exact helpers
def is_prime(n: int, ps: np.ndarray) -> bool:
    """trial division by `ps` (all primes <= sqrt n); vectorised, no Python loop."""
    r = int(np.sqrt(n)) + 1
    q = ps[ps <= r]
    return bool(n > 1 and not np.any(n % q == 0))


def band_upto(N: int):
    """eq:band for small N, by numpy: (band, vB, vR) as sorted arrays. No loop over n."""
    Q = int(np.sqrt(N))
    n = np.arange(1, N + 1, dtype=np.int64)
    # omega and squarefreeness by striding over the primes <= N (a handful here)
    ps = primes_upto(N + 1)
    om = np.zeros(N + 1, dtype=np.int64)
    sf = np.ones(N + 1, dtype=bool)
    pmax = np.zeros(N + 1, dtype=np.int64)
    for p in ps.tolist():            # Python loops over PRIMES <= 30 only, never over the band
        om[p::p] += 1
        pmax[p::p] = p
        if p * p <= N:
            sf[p * p::p * p] = False
    thr_ = N // int(ps[ps > Q][0])
    sel = (n > thr_) & sf[1:] & (pmax[1:] <= Q)
    b = n[sel]
    return b, b[om[b] % 2 == 0], b[om[b] % 2 == 1], Q, thr_


def main() -> int:
    say("pa_tex_constants_v2 -- the constants printed by the manuscript")
    say("source: the manuscript (the same .tex the reproduction packet carries)")
    say("clauses V1-V23, M1-M3, X1, D1 in the file header; mpmath dps %d" % mp.mp.dps)
    say("rendering rule from the source: \\ldots -> truncated ; \\approx / near / >= -> rounded")
    say()

    eg = mp.exp(mp.euler)
    c15 = mp.sqrt(15 / mp.pi ** 2)
    c1 = (mp.mpf(4) / mp.sqrt(6)) * c15
    ee = mp.e ** (mp.e - 1)
    c2 = 4 * mp.sqrt(2 * ee)
    NR = mp.pi ** 2 / (3 * (1 - mp.log(2)))
    v = mp.findroot(lambda z: z * mp.log(z) - 1, mp.mpf("1.76"))
    l2 = 1 - mp.log(2)
    r09 = (mp.mpf(3) / mp.pi ** 2) * l2
    b4 = 4 * c15 / mp.sqrt(4)
    c46 = 4 * 6 * eg / (3 * l2)

    say("closed forms, 25 significant digits")
    for nm, x in (("c_1", c1), ("c_2", c2), ("pi^2/(3(1-log2))", NR), ("sqrt(15/pi^2)", c15),
                  ("e^(e-1)", ee), ("v (v log v = 1)", v), ("1-log2", l2),
                  ("(3/pi^2)(1-log2)", r09), ("4*6*e^g/(3(1-log2))", c46)):
        say("   %-22s = %s" % (nm, mp.nstr(x, 25, strip_zeros=False)))
    say()

    clause("V1", c1, "2.0131684", "ldots", 7, "(4/sqrt6)sqrt(15/pi^2)")
    clause("V2", c2, "13.3565762", "ldots", 7, "4sqrt(2e^(e-1))")
    plain("V3", rndup(c1, 5) == "2.01317" and mp.mpf("2.01317") > c1,
          "2.01317 = c_1 rounded up to 5 dp (%s) and exceeds c_1 = %s"
          % (rndup(c1, 5), mp.nstr(c1, 12)))
    clause("V4", NR, "10.7213228", "ldots", 7, "pi^2/(3(1-log2))")
    clause("V5", c15, "1.2328088", "ldots", 7, "sqrt(zeta(2)/zeta(4))")
    clause("V6", ee, "5.5749415", "ldots", 7, "Xi = e^(e-1)")
    clause("V7a", eg / v, "1.0101232", "ldots", 7, "sup varpi e^gamma = e^g/v")
    clause("V7b", v + 1, "2.7632228", "ldots", 7, "u* = 1+v")
    clause("V7c", 1 / v, "0.5671432", "ldots", 7, "varpi(u*) = 1/v")
    # the text's own derivation of V7: u varpi(u) = 1 + log(u-1) on [2,3], max at u = 1+v
    uw = lambda u: (1 + mp.log(u - 1)) / u
    plain("V7d", mp.almosteq(uw(v + 1), 1 / v, 1e-30)
          and mp.almosteq(mp.diff(uw, v + 1), 0, 1e-25),
          "the text's derivation: varpi(1+v) = (1+log v)/(1+v) = 1/v and d/du varpi = 0 there")
    clause("V8", l2, "0.3068528", "ldots", 7, "1-log2")
    clause("V9a", mp.mpf(2) ** mp.mpf("-0.75"), "0.5946", "ldots", 4, "2^(-3/4)")
    clause("V9b", mp.mpf(3) ** mp.mpf("-0.75"), "0.4386", "ldots", 4, "3^(-3/4)")
    clause("V10", r09, "0.0932720", "ldots", 7, "(3/pi^2)(1-log2)")
    plain("V10b", r09 > mp.mpf("0.09"), "(3/pi^2)(1-log2) = %s > 0.09, so R >= 0.09N"
          % mp.nstr(r09, 10))
    clause("V11", b4, "2.4656", "ldots", 4, "4sqrt(15/pi^2)/2")
    clause("V12", c46, "46.4345", "ldots", 4, "4*beta*e^g/(3(1-log2)), beta=6")

    say()
    say("the two explicit prime-counting inputs, as the text deduces them")
    plain("V13", rndup(mp.mpf("1.25506"), 2) == "1.26" and mp.mpf("1.26") > mp.mpf("1.25506"),
          "(C): 1.26 = 1.25506 rounded up to 2 dp and exceeds it, so pi(y) <= 1.26 y/log y follows")
    cpi = 1 / mp.sqrt(mp.mpf("6.455"))
    clause("V14a", cpi, "0.3935", "ldots", 4, "1/sqrt(6.455)")
    plain("V14b", mp.mpf("0.39") < cpi,
          "(P): c_pi = 0.39 < 1/sqrt(6.455) = %s, so exp(-sqrt(log x/6.455)) <= exp(-0.39 sqrt(log x))"
          % mp.nstr(cpi, 12))
    plain("V14c", mp.log(mp.e) ** mp.mpf("-0.75") <= 1,
          "(P): (log x)^(-3/4) <= 1 for x >= e (equality at x = e), the text's logarithmic factor")

    say()
    f = lambda L: c1 / (L * mp.sqrt(mp.log(L))) + c2 * L * L * mp.e ** (-L / (24 * mp.log(L)))
    Ls = mp.findroot(lambda L: f(L) - 1, mp.mpf(3767))
    say("   crossing: L* = %s, f(L*) = %s, 2L*/log2 = %s"
        % (mp.nstr(Ls, 14), mp.nstr(f(Ls), 8), mp.nstr(2 * Ls / mp.log(2), 14)))
    clause("V15a", Ls, "3767", "near", 0, "crossing L*")
    clause("V15b", 2 * Ls / mp.log(2), "10868", "approx", 0, "2L*/log2")

    L16 = mp.mpf(3767)
    Z0 = mp.e ** (L16 / mp.log(L16) ** 3)
    Y16 = L16 ** 6
    say("   at L = 3767: Z_0 = %s, Y = L^6 = %s, u = %s"
        % (mp.nstr(Z0, 10), mp.nstr(Y16, 8), mp.nstr(mp.log(Z0) / mp.log(Y16), 6)))
    clause("V16a", Z0, "852", "approx", 0, "Z_0 at L=3767")
    clause("V16b", Y16, "2.9e21", "approx", 1, "Y = L^6 at L=3767", mantissa=True)
    clause("V16c", mp.log(Z0) / mp.log(Y16), "0.14", "approx", 2, "u = logZ_0/logY")

    Lu = mp.findroot(lambda L: L / (6 * mp.log(L) ** 4) - 1, mp.mpf(1.0e5))
    say("   u(Z_0) = L/(6(log L)^4) reaches 1 at L = %s, 2L/log2 = %s, L/3767 = %s"
        % (mp.nstr(Lu, 12), mp.nstr(2 * Lu / mp.log(2), 12), mp.nstr(Lu / 3767, 6)))
    clause("V17a", Lu, "1.08e5", "near", 2, "root of L/(6(logL)^4)=1", mantissa=True)
    clause("V17b", 2 * Lu / mp.log(2), "3.1e5", "approx", 1, "log2 N there", mantissa=True)
    clause("V17c", Lu / 3767, "29", "approx", 0, "factor in L beyond the crossing")

    L6 = mp.findroot(lambda L: 6 * mp.log(L) - L, mp.mpf(17))
    L12 = mp.findroot(lambda L: 12 * mp.log(L) - L, mp.mpf(46))
    say("   onsets: root(6logL=L) = %s, root(12logL=L) = %s" % (mp.nstr(L6, 12), mp.nstr(L12, 12)))
    clause("V18a", L6, "16.9988", "ldots", 4, "root of 6logL=L")
    clause("V18b", 2 * L6 / mp.log(2), "49.05", "ge", 2, "log2 N onset for Y<Q")
    clause("V18c", 2 * L12 / mp.log(2), "132.5", "ge", 1, "log2 N onset for Y<=sqrtQ")

    Q28, L28 = 2 ** 14, mp.log(2 ** 14)
    Y28 = L28 ** 6
    say("   at 2^28: Q = %d, L = %s, Y = L^6 = %s" % (Q28, mp.nstr(L28, 10), mp.nstr(Y28, 10)))
    clause("V19a", Y28, "8.4e5", "approx", 1, "Y = L^6 at 2^28", mantissa=True)
    plain("V19b", Y28 > Q28, "Y = %s > Q = %d at 2^28, so every band element is Y-smooth and x = 0"
          % (mp.nstr(Y28, 8), Q28))

    # ---------------- the N = 2^52 witness ------------------------------------------------------
    say()
    say("the N = 2^52 witness (exact integers except the two rendered reals)")
    N52 = 2 ** 52
    Q52 = 2 ** 26
    ps = primes_upto(9000)                      # > sqrt(q_1) = 8192.0..., covers every check below
    L52 = mp.log(Q52)
    Y52 = L52 ** 6
    plain("V20a", int(mp.floor(mp.sqrt(N52))) == Q52 == 67108864,
          "Q = floor(sqrt(2^52)) = 2^26 = %d" % Q52)
    clause("V20b", L52, "18.0218", "ldots", 4, "L = log Q at 2^52")
    clause("V20c", Y52, "34260433.5", "ldots", 1, "Y = L^6 at 2^52")
    plain("V20d", Y52 < Q52, "Y = %s < Q = %d" % (mp.nstr(Y52, 12), Q52))
    cand = np.arange(Q52 + 1, Q52 + 40, dtype=np.int64)
    q1 = int(cand[[is_prime(int(c), ps) for c in cand.tolist()].index(True)])
    plain("V20e", q1 == 67108879, "q_1 = least prime > Q = %d (prime: %s)"
          % (q1, is_prime(q1, ps)))
    thr52 = N52 // q1
    plain("V20f", thr52 == 67108849, "thr = floor(2^52/%d) = %d" % (q1, thr52))
    pA, pB = 34260449, 34260493
    pr = is_prime(pA, ps) and is_prime(pB, ps)
    inw = all(Y52 < x <= Q52 for x in (pA, pB))
    plain("V20g", pr and inw,
          "%d and %d both prime (%s) and both in (Y,Q] (%s)" % (pA, pB, pr, inw))
    nB, nR = 2 * pA, 2 * pA * pB
    # vB / vR membership straight from eq:band and eq:parts
    okB = (thr52 < nB <= N52) and len({2, pA}) == 2 and max(2, pA) <= Q52 and 2 % 2 == 0
    okR = (thr52 < nR <= N52) and len({2, pA, pB}) == 3 and max(2, pA, pB) <= Q52 and 3 % 2 == 1
    plain("V20h", okB, "2*%d = %d in vB: thr < n <= N, squarefree, P+ = %d <= Q, omega = 2 even"
          % (pA, nB, pA))
    plain("V20i", okR, "2*%d*%d = %d in vR: thr < n <= N, squarefree, P+ = %d <= Q, omega = 3 odd"
          % (pA, pB, nR, pB))
    # small-prime part a_Y = product of the prime factors <= Y; here only 2 qualifies
    aB = 2 if 2 <= Y52 else 1
    aB *= pA if pA <= Y52 else 1
    aR = aB * (pB if pB <= Y52 else 1)
    sq = int(mp.floor(mp.sqrt(Q52)))
    plain("V20j", aB == 2 and aR == 2 and 2 <= sq,
          "small-prime part a_Y = %d for both (only 2 is <= Y), and 2 <= sqrt(Q) = %d"
          % (aB, sq))
    plain("V20k", okB and okR and aB == 2 and aR == 2,
          "hence B_2 > 0 and R_2 > 0, so eq:X carries positive weight at N = 2^52")

    # ---------------- the level examples of section 1.1 -----------------------------------------
    say()
    say("the level examples of section 1.1")
    Qex = 3
    plain("V21", (Qex + 1) ** 2 - Qex ** 2 == 2 * Qex + 1,
          "the integers N with Q^2 <= N < (Q+1)^2 number (Q+1)^2 - Q^2 = 2Q+1")
    b10, _, _, q10, t10 = band_upto(10)
    b11, _, _, q11, t11 = band_upto(11)
    plain("V22", (q10, t10, b10.tolist()) == (3, 2, [3, 6])
          and (q11, t11, b11.tolist()) == (3, 2, [3, 6]),
          "N=10: Q=%d thr=%d band=%s ; N=11: Q=%d thr=%d band=%s"
          % (q10, t10, b10.tolist(), q11, t11, b11.tolist()))
    _, vB25, vR25, q25, _ = band_upto(25)
    b30, vB30, vR30, q30, _ = band_upto(30)
    rows5 = np.array([d for d in range(1, 6) if d in (1, 2, 3, 5)], dtype=np.int64)  # squarefree<=5
    rows5 = np.array([1, 2, 3, 5], dtype=np.int64)
    inc = lambda l: (l % rows5 == 0).astype(np.float64)
    t30 = sum(inc(int(x)) for x in vR30.tolist()) / len(vR30)
    half = (inc(10) + inc(15)) / 2
    plain("V23", vR25.tolist() == [5] and vR30.tolist() == [5, 30]
          and np.allclose(t30, half) and set((6, 10, 15)) <= set(vB30.tolist()),
          "Q=5: vR(25)=%s, vR(30)=%s, t(30)=%s = (A_10+A_15)/2 on rows %s, and 10,15 in vB(30)"
          % (vR25.tolist(), vR30.tolist(), t30.tolist(), rows5.tolist()))

    # ---------------- measured ------------------------------------------------------------------
    say()
    say("measured clauses -- the manuscript's digits against the stored result files")

    def rd(name):
        p = os.path.join(RES, name)
        return io.open(p, encoding="utf-8").read() if os.path.exists(p) else None

    sm, ek20 = rd("conedual_separator_mass.txt") or "", rd(
        "conedual_exact_kappa_small_cells_20.txt") or ""
    exact, cur = {}, None
    for line in sm.splitlines():
        m = re.match(r"\s*2\^(\d+)\s", line)
        if m:
            cur = int(m.group(1))
        m = re.search(r"exact dual:.*?kappa=([0-9.]+)", line)
        if m and cur is not None:
            exact[cur] = mp.mpf(m.group(1))
    m = re.search(r"kappa\s+([0-9.]+)", ek20)
    if m:
        exact[20] = mp.mpf(m.group(1))
    encl = {16: ("0.1723893", "0.1723894"), 18: ("0.1408604", "0.1408605"),
            20: ("0.1055426", "0.1055428")}
    m1ok, m1d = True, []
    for e, (lo, hi) in encl.items():
        k = exact.get(e)
        if k is None:
            m1ok = False
            m1d.append("2^%d NO STORED VALUE" % e)
            continue
        inside = mp.mpf(lo) <= k <= mp.mpf(hi)
        m1ok = m1ok and inside and len(lo.split(".")[1]) == 7 and len(hi.split(".")[1]) == 7
        m1d.append("2^%d [%s,%s] contains %s: %s" % (e, lo, hi, mp.nstr(k, 10),
                                                     "yes" if inside else "NO"))
    plain("M1", m1ok, "enclosures at 2^16,18,20 (7 dp endpoints): " + "; ".join(m1d))

    big = {22: "conedual_colgen_2e22_n0b.json", 24: "conedual_colgen_2e24_n0b.json",
           26: "conedual_colgen_2e26_n0b.json", 28: "conedual_colgen_2e28_n28.json"}
    pl = {22: ("0.0763310", "0.0763313"), 24: ("0.0527620", "0.0527630"),
          26: ("0.0346220", "0.0346234"), 28: ("0.0245260", "0.0245279")}
    m2ok, m2d, ub = True, [], {}
    for e, fn in big.items():
        txt = rd(fn)
        if txt is None:
            m2ok = False
            m2d.append("2^%d %s MISSING" % (e, fn))
            continue
        j = json.loads(txt)
        lo, hi = pl[e]
        g1, g2 = rnd(mp.mpf(j["c6"]) / mp.mpf(10) ** 6, 7), rnd(mp.mpf(j["ub"]), 7)
        ub[e] = mp.mpf(hi)
        ok = g1 == lo and g2 == hi
        m2ok = m2ok and ok
        m2d.append("2^%d c6=%s->%s|%s ub=%s->%s|%s %s"
                   % (e, j["c6"], g1, lo, j["ub"], g2, hi, "ok" if ok else "NO"))
    plain("M2", m2ok, "enclosures at 2^22..2^28: " + "; ".join(m2d))

    pj = rd("conedual_parity_domination.json")
    m3lit = {16: "0.6982", 18: "0.6130", 20: "0.6088", 22: "0.5841"}
    if pj is None:
        plain("M3", False, "conedual_parity_domination.json MISSING")
    else:
        cells = {c["e"]: c for c in json.loads(pj)["cells"]}
        up = {16: mp.mpf("0.1723894"), 18: mp.mpf("0.1408605"), 20: mp.mpf("0.1055428")}
        up.update(ub)
        m3ok, m3d, facs = True, [], []
        for e, lit in m3lit.items():
            c = cells.get(e)
            if c is None:
                m3ok = False
                m3d.append("2^%d missing" % e)
                continue
            b = mp.mpf(c["best_kappa_ub"])
            got, yok = rndup(b, 4), c["best_Y"] in (3, 5)
            facs.append(b / up[e])
            ok = got == lit and yok
            m3ok = m3ok and ok
            m3d.append("2^%d %s->%s|%s Y=%d %s" % (e, mp.nstr(b, 12), got, lit, c["best_Y"],
                                                   "ok" if ok else "NO"))
        fok = bool(facs) and 4 <= min(facs) and max(facs) <= 8
        plain("M3", m3ok and fok, "construction rounded up: " + "; ".join(m3d)
              + " | factors %s in [4,8]: %s" % (["%.2f" % x for x in facs],
                                                "yes" if fok else "NO"))

    # ---------------- X1: the removed literals are absent ---------------------------------------
    say()
    tex = None
    # The packet's own paper source is looked for FIRST, and nothing outside the directory this script
    # sits in is consulted at all: looking elsewhere first would let a distributed copy read a different
    # manuscript from the one shipped beside it.  Two candidates only: the packet layout (the .tex
    # beside `code/`) and the working layout (`../paper/`).
    for cand_p in (os.path.join(os.path.dirname(HERE),
                                "P7-parity-separation-tends-to-zero.tex"),
                   os.path.join(os.path.dirname(HERE), "paper", "kappa_tends_to_zero.tex")):
        if os.path.exists(cand_p):
            tex = io.open(cand_p, encoding="utf-8").read()
            break
    if tex is None:
        REC.append(dict(name="X1", holds=None, detail="the .tex is not in this tree; not scored"))
        say("X1    INCONCL. the .tex is not in this tree (it is not part of the code packet);"
            " absence of the removed literals not scored")
    else:
        gone = {s: (s not in tex) for s in ("1.33258", "0.52", "0.8")}
        plain("X1", all(gone.values()),
              "removed literals absent from the .tex: " + ", ".join(
                  "%s %s" % (k, "absent" if x else "STILL PRESENT") for k, x in gone.items()))

    # ---------------- display audit -------------------------------------------------------------
    say()
    if MISM:
        say("D1  DISPLAY AUDIT -- %d literal(s) matched against their own notation:" % len(MISM))
        for name, lit, nota, which, t, r in MISM:
            say("    %-5s literal %-14s written as \\%-7s but matches %s (trunc %s, round %s)"
                % (name, lit, nota, which, t, r))
        say("    One-character display choices in the .tex, not wrong constants.")
    else:
        say("D1  DISPLAY AUDIT: every literal whose truncation and rounding differ matches the"
            " rendering its notation declares. HOLDS")

    say()
    say("VERDICT  " + " | ".join(
        "%s %s" % (r["name"], "INCONCL." if r.get("holds") is None
                   else ("HOLDS" if r["holds"] else "FAILS")) for r in REC))
    say("failures: %s" % (", ".join(BAD) if BAD else "none"))
    say("D1: %s" % ("%d rendering mismatch(es)" % len(MISM) if MISM else "clean"))

    with io.open(OUTJ, "w", encoding="utf-8", newline="\n") as f:
        json.dump(dict(dps=mp.mp.dps, source="P7-parity-separation-tends-to-zero.tex",
                       clauses=REC, failures=BAD,
                       display_mismatches=[dict(name=n, literal=l, notation=no, matches=w,
                                                trunc=t, round=r) for n, l, no, w, t, r in MISM],
                       crossing=dict(L=mp.nstr(Ls, 16), log2N=mp.nstr(2 * Ls / mp.log(2), 16)),
                       witness_2e52=dict(Q=Q52, L=mp.nstr(L52, 16), Y=mp.nstr(Y52, 16),
                                         q1=q1, thr=thr52, p=pA, q=pB, nB=nB, nR=nR, a_Y=aB),
                       omega_constant=mp.nstr(1 / v, 25)), f, indent=1)
    say("wrote %s and %s" % (os.path.basename(OUT), os.path.basename(OUTJ)))
    return 1 if BAD else 0


if __name__ == "__main__":
    sys.exit(main())

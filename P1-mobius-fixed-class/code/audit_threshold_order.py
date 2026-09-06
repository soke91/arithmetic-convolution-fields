# -*- coding: utf-8 -*-
r"""
paper/theorem_A.md, Remark {#rem:thresholdorder}.

WHAT IS AT STAKE

[prop:onesided] states its threshold as S(N)(1 - A(N))N.  When this
audit was written the proposition recorded that the threshold is
"never smaller than c N/(log N log log N)"; that phrase is no longer
in the paper, and what stands there now -- asymp N over the even N
measured, and asymp N/(log N) on the family where it is thinnest -- is
this audit's outcome carried into the text.  The old phrase is kept
here because it is what was audited.  That bound is
obtained by keeping 1 - A(N) asymp 1/(p log p) with p asymp log N and
throwing S(N) away as S(N) >> 1.  On the family where the threshold is
thin, S(N) is not bounded: it grows.

For N carrying every odd prime up to y,

    S(N) = 2 C_2 prod_{3<=p<=y} (p-1)/(p-2),
    log S(N) - log(2 C_2) = sum_{p<=y} log(1 + 1/(p-2))
                          = sum_{p<=y} 1/p + O(1)
                          = log log y + O(1),

so S(N) asymp log y, while 1 - A(N) = sum_{p>y} 1/(p(p-1))
asymp 1/(y log y).  Since log N = theta(y) asymp y,

    S(N)(1 - A(N)) asymp (log y)/(y log y) = 1/y asymp 1/log N.

THE THRESHOLD ON THE ADVERSE FAMILY IS ASYMP N/log N, NOT
N/(log N log log N).  The recorded bound is true and one factor of
log log N slack: it prices the fall of 1 - A(N) and not the rise of
S(N) that accompanies it, and the two are driven by the same primes.

What that buys is a weaker demand.  [prop:onesided] asks
E_3 > -S(N)(1 - A(N))N, so a sharper lower bound on the right-hand
side is a weaker hypothesis, by a factor log log N at the N where the
hypothesis is hardest to meet.  It does not move [rem:kappaworth]'s
ceiling, which is a statement about 1/(1 - A(N)) alone and stands.

BACKS: Remark {#rem:thresholdorder} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Q1  The right normalisation is log N: over even N in
      [1e5, 1.6e7] the infimum of S(N)(1 - A(N)) log N lies in
      [0.5, 1.5].

  Q2  The recorded normalisation over-corrects: along the four even
      primorials in range (2310, 30030, 510510, 9699690) the product
      S(1-A) log N log log N is strictly increasing, while
      S(1-A) log N rises over the same four by less than half as much
      in relative terms.

  Q3  The reconstruction is the recorded one: the minimum of
      S(1-A) log N log log N over even N >= 1e5 in the field
      reproduces the 2.482019 that [rem:onesided] records, to
      within 0.001, and at the same N.

REFUTATION RULE (fixed before the run)

  Q1  REFUTED if the infimum falls outside [0.5, 1.5].
  Q2  REFUTED if the loglog-normalised product is not strictly
      increasing along the four, or if the log-normalised one rises
      by at least half as much in relative terms.
  Q3  REFUTED if the minimum differs from 2.482019 by more than
      0.001, or is attained at a different N.
  Non-zero exit on any refutation.

WHAT EACH REFUTATION WOULD MEAN (fixed before the run, per M9)

  Q1 refuted -> 1/log N is not the order of the threshold on this
     field, so the Mertens estimate for S(N) does not compose with the
     tail estimate for 1 - A(N) the way the derivation assumes.  Below
     0.5 would mean the true order is smaller than 1/log N and the
     recorded bound was closer to sharp than claimed here; above 1.5
     would mean the field has not reached the regime.
  Q2 refuted -> the two normalisations are not distinguishable on four
     primorials, and the claim that one over-corrects rests on a
     comparison the field cannot make.  The run prints both columns at
     all four so the reader can see how much room there was.
  Q3 refuted -> this run's S and A are not the ones [rem:onesided]
     used, and Q1 and Q2 are about a different quantity.
"""

import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "audit_threshold_order.txt")

X = 16_000_000
LO = 100_000
PRIMORIALS = [2310, 30030, 510510, 9699690]
RECORDED_MIN = 2.482019
RECORDED_AT = 510510


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(math.isqrt(n)) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0]


def main():
    lines = []

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    say("STATISTIC: for each even N, S(N)(1 - A(N)) with S and A built")
    say("           by one pass over the odd primes, and that quantity")
    say("           multiplied by log N and by log N log log N.  The")
    say("           infimum of each over the field, the N attaining it,")
    say("           and both columns along the four even primorials in")
    say("           range.  Nothing is centred, detrended or fitted.")
    say()
    say("FIELD: every even N in [%d, %d], one number of each kind per"
        % (LO, X))
    say("       N.  The four primorials are read out of the same")
    say("       arrays, not recomputed.")
    say()
    say("DENOM: each rise is divided by that column's value at the first N;")
    say("NULL: none is needed.  Every quantity is a closed form")
    say("      evaluated on a complete range; the comparison in Q2 is")
    say("      between two normalisations of the same numbers, so each")
    say("      is the other's control, and Q3's control is the recorded")
    say("      value it must reproduce.")
    say()

    say("sieving to %d ..." % X)
    pr = primes_upto(X)
    odd = pr[pr > 2].astype(np.float64)
    fac = 1.0 - 1.0 / (odd * (odd - 1.0))
    a_inf = float(np.exp(np.log(fac).sum()))
    c2 = float(np.exp(np.log(1.0 - 1.0 / (odd - 1.0) ** 2).sum()))
    say("  A_inf = %.10f   2*C_2 = %.10f" % (a_inf, 2.0 * c2))

    A = np.full(X + 1, a_inf, dtype=np.float64)
    S = np.full(X + 1, 2.0 * c2, dtype=np.float64)
    for i in range(odd.size):
        p = int(odd[i])
        A[p::p] /= fac[i]
        S[p::p] *= (p - 1.0) / (p - 2.0)

    ev = np.arange(LO + (LO & 1), X + 1, 2, dtype=np.int64)
    thr = S[ev] * (1.0 - A[ev])
    lg = np.log(ev.astype(np.float64))
    llg = np.log(lg)
    q_log = thr * lg
    q_loglog = thr * lg * llg

    say()
    say("Q1  the right normalisation is log N")
    say("=" * 74)
    i1 = int(np.argmin(q_log))
    v1 = float(q_log[i1])
    q1 = bool(0.5 <= v1 <= 1.5)
    say("  inf S(1-A) log N        = %.6f  at N = %d"
        % (v1, int(ev[i1])))
    say("  inf S(1-A)              = %.6f  at N = %d"
        % (float(thr.min()), int(ev[int(np.argmin(thr))])))
    say("  -> Q1 %s  (bracket [0.5, 1.5])" % ("HOLDS" if q1 else
                                              "REFUTED"))
    say()

    say("Q2  the recorded normalisation over-corrects")
    say("=" * 74)
    say("  N            S(1-A)     x log N    x log N loglog N")
    say("  " + "-" * 58)
    a_col, b_col = [], []
    for N in PRIMORIALS:
        t = float(S[N] * (1.0 - A[N]))
        l = math.log(N)
        ll = math.log(l)
        a_col.append(t * l)
        b_col.append(t * l * ll)
        say("  %-12d %-10.6f %-10.6f %-10.6f" % (N, t, t * l, t * l * ll))
    inc_b = all(b_col[i] < b_col[i + 1] for i in range(len(b_col) - 1))
    rel_a = a_col[-1] / a_col[0] - 1.0
    rel_b = b_col[-1] / b_col[0] - 1.0
    q2 = bool(inc_b and rel_a < 0.5 * rel_b)
    say("  loglog column strictly increasing : %s" % inc_b)
    say("  relative rise, log column         : %.6f" % rel_a)
    say("  relative rise, loglog column      : %.6f" % rel_b)
    say("  half of the loglog rise           : %.6f" % (0.5 * rel_b))
    say("  -> Q2 %s" % ("HOLDS" if q2 else "REFUTED"))
    say()

    say("Q3  the reconstruction is the recorded one")
    say("=" * 74)
    i3 = int(np.argmin(q_loglog))
    v3 = float(q_loglog[i3])
    n3 = int(ev[i3])
    q3 = bool(abs(v3 - RECORDED_MIN) <= 0.001 and n3 == RECORDED_AT)
    say("  inf S(1-A) log N loglog N = %.6f  at N = %d" % (v3, n3))
    say("  recorded                  = %.6f  at N = %d"
        % (RECORDED_MIN, RECORDED_AT))
    say("  -> Q3 %s" % ("HOLDS" if q3 else "REFUTED"))
    say()

    say("verdict")
    say("=" * 74)
    for nm, ok in (("Q1", q1), ("Q2", q2), ("Q3", q3)):
        say("  %s %s" % (nm, "HOLDS" if ok else "REFUTED"))

    # newline="\n" is not cosmetic. Without it Python's text mode writes CRLF
    # on Windows and LF elsewhere, so the same script and the same numbers
    # give a different SHA-256 -- and this file's hash is pinned in P1's
    # PACKET.json, which is published.
    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return 0 if (q1 and q2 and q3) else 1


if __name__ == "__main__":
    sys.exit(main())

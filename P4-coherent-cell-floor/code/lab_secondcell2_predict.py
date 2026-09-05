# -*- coding: utf-8 -*-
r"""P4 — the second cell shape, registered a SECOND time on an exact
statistic.

WHY A SECOND REGISTRATION

The first (lab_secondcell_predict.py, sealed 3b1b7d1) was refuted.  It
was refuted by luck, not by the prediction: its rule scored ONE random
draw at each of three cells against a three-sigma cap, and the draw at
depth 2 came out at +3.22 of a sampling error the run reported
correctly -- thirty independent redraws centre on the sealed prediction
at -0.25.  A cap of that shape is refused by chance about once in a
hundred and twenty runs and this was that.

**The rule stands refuted and is not rescored.**  What follows is a
new registration on a statistic that carries no draw at all.

WHAT CHANGES: THE STATISTIC, NOT THE PREDICTION

The four predictions are unchanged and are READ from the first seal.
They are not recomputed here; a recomputation could drift from what was
sealed and nothing would say so.

The measurement changes from a sampled pair average to an EXACT one.
A cell is a set of positions in the band, so

    E_same,c[S_2] = sum_h S_2(h) (1_c * 1_c)(h) / (n_c (n_c - 1)),

the autocorrelation of the cell's indicator, computable by one FFT per
cell over the band's 10^6 even positions.  There is no sampling, no
seed, and no standard error; the tolerances below are absolute.

WHAT THE EXACT MEASUREMENT CAN NOW DECIDE THAT THE SAMPLED ONE COULD NOT

The prediction has two parts and they fail differently:

  MAIN TERM     K E_same,c[A'] - 2, with K = 175/128 and E_same,c[A']
                exact rationals.  No irrational input, no model.
  DEFICIT       the same finite divisor deficit as the first partition,
                DIVIDED BY m.  The 1/m is a MODEL: it assumes the
                m^2 - m off-diagonal class pairs cancel, which is
                precisely what {#prop:scaleinv} declines to claim.

At depth 3 (m = 1) the prediction and the measurement are the same sum
and the comparison is a consistency check, as the first seal recorded.
At depths 0 to 2 they are not: the prediction models the divisor part,
the exact measurement computes it.  So the exact statistic separates
the main term from the dilution model, which no sampled run could.

REGISTERED BEFORE THE EXACT MEASUREMENT EXISTS

  T1  At depths 0, 1 and 2 the exact D'_c lies BETWEEN the uncorrected
      prediction and the fully corrected one -- that is, the deficit's
      sign is right and its undiluted size is an upper bound.
      REFUTED if any of the three falls outside that interval on either
      side.  This is the claim that does not depend on the dilution
      model at all, and it GATES: if it fails, the divisor deficit is
      not the mechanism at cells of more than one class.

      **T1 IS ALSO REGISTERED AGAINST EXPECTATION AT DEPTH 2.**  The
      first run's thirty-draw mean there was 0.825796, and T1's upper
      bound is the uncorrected 0.825521.  That is 0.000275 above it,
      about 1.2 standard errors of a thirty-draw mean, so the sampled
      evidence leans towards T1 failing at depth 2 -- towards the
      measured value being ABOVE the main term where the deficit says
      it must be below.  It is not decisive, which is why the exact run
      is worth making; but it is recorded here so that a pass cannot
      later be read as having been expected.

  T2  |exact D'_c - prediction corrected by 1/m| <= 1e-4 at each of
      depths 0, 1 and 2.  REFUTED otherwise.
      **Registered in the open knowledge that it may fail.**  The first
      run's thirty-draw mean at depth 2 was about 0.8258, which sits
      between the 1/m-corrected 0.825475 and the undiluted 0.825521 --
      nearer the middle than either.  If 1/m were right the exact value
      would sit at 0.825475.  Registering a rule one expects to fail is
      the point of registering it; the alternative is to look first and
      then choose the model.

  T3  The implied dilution exponent alpha, defined by
      correction = (undiluted correction) / m^alpha, is the same at
      depths 0, 1 and 2 to within 0.15 in absolute spread.  REFUTED if
      the spread exceeds 0.15, and REFUTED if alpha is undefined at any
      depth -- which happens when the exact value sits ABOVE the main
      term, since no positive deficit then fits it.

      WHAT MAKES THIS A TEST AND NOT A CURVE FIT.  With a free exponent
      a law of this shape fits almost any monotone sequence over a
      short range, so the tolerance has to be the whole content.  Two
      things give it some: m is 48, 44 and 12, so depths 0 and 1 are
      nearly the same m and their alphas are nearly forced to agree --
      **the discriminating comparison is depth 2 against depths 0
      and 1**, and 0.15 there corresponds to the gaps differing by
      about half from a common law.  And alpha is not free in
      interpretation: 1 is the diagonal-class-pair model, 0 is no
      dilution at all, and anything between has no mechanism proposed
      for it, so a value near 0.3 would say the model is wrong rather
      than that the exponent is 0.3.

      WHAT A T3 PASS WITH T2 REFUTED WOULD LICENSE, said before the
      run rather than after.  It would say the deficit is diluted by
      some power of m that is not 1, so the diagonal-class-pair
      argument is wrong about the mechanism even though a dilution
      exists.  **The 1/m correction then comes back out of the first
      partition's depth-4 entry**, where it was applied two hours ago
      on a consistency argument, and that entry returns to the pure
      main term (+1.16 sampling errors rather than +0.34).  Recording
      that now so the decision to unwind it cannot be read as having
      been made in the light of the number.

  T1 gates.  T2 and T3 are reported.

WHAT WOULD MAKE THIS A FAILURE, SAID PLAINLY

T1 refuted: the divisor deficit does not explain the departure from the
main term at multi-class cells, and the main-term computation's
transfer to a new cell shape is not established.  The paper's claim
narrows to the (3,5,7,11,13) partition.

T1 held with T2 refuted: the mechanism transfers and the 1/m model does
not.  The paper keeps the main term and drops 1/m, and the deficit
correction at multi-class cells becomes a measured quantity rather than
a modelled one -- which would also remove it from the first partition's
depth-4 entry, where it is currently applied.

WHAT THIS STILL DOES NOT TEST

Same band, same S_2 field, same code as the first partition.  What is
tested is dependence on CELL SHAPE.  A systematic error in the field
passes through both partitions identically.  Not an independent
replication and not offered as one.

**THIS FILE MEASURES NOTHING.**  lab_secondcell2.py does not exist.

FIELD: predictions read from results/lab_secondcell_predict.txt;
       m = 48, 44, 12, 1 from the pattern counts; no measurement.

BACKS: nothing yet.  A registration, not evidence.
"""

import io
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "lab_secondcell2_predict.txt")
SEAL = os.path.join(ROOT, "results", "lab_secondcell_predict.txt")

lines = []


def say(s=""):
    lines.append(s)
    print(s)


def read_seal():
    """corrected, uncorrected, m, and the undiluted correction, all from
    the first seal.  Read and not recomputed."""
    src = io.open(SEAL, encoding="utf-8", errors="replace").read()
    corr, unc, ms, full = {}, {}, {}, {}
    b = src.split("depth  corrected prediction")
    if len(b) > 1:
        for ln in b[1].splitlines()[1:]:
            f = ln.split()
            if len(f) == 2 and f[0].isdigit():
                corr[int(f[0])] = float(f[1])
            elif corr:
                break
    b = src.split("m (classes)  n_c approx   uncorrected prediction")
    if len(b) > 1:
        for ln in b[1].splitlines()[1:]:
            f = ln.split()
            if len(f) == 4 and f[0].isdigit():
                unc[int(f[0])] = float(f[3])
                ms[int(f[0])] = int(f[1])
            elif unc:
                break
    b = src.split("depth  m      n_c per class  deficit      correction"
                  "   diluted")
    if len(b) > 1:
        for ln in b[1].splitlines()[1:]:
            f = ln.split()
            if len(f) == 6 and f[0].isdigit():
                full[int(f[0])] = float(f[4])
            elif full:
                break
    return corr, unc, ms, full


def read_sampled_mean():
    """The first run's thirty-draw mean at depth 2, READ from its result
    file.  Typed, it would be a measurement written by hand."""
    p = os.path.join(ROOT, "results", "lab_secondcell.txt")
    if not os.path.exists(p):
        return None
    for ln in io.open(p, encoding="utf-8", errors="replace"):
        t = ln.split()
        if len(t) >= 6 and t[0] == "mean" and t[2] == "spread":
            try:
                return float(t[1].rstrip(","))
            except ValueError:
                return None
    return None


def main():
    corr, unc, ms, full = read_seal()
    if len(corr) != 4 or len(unc) != 4 or len(full) != 4:
        say("the first seal was not readable -- refusing to register")
        return 1
    say("READ FROM THE FIRST SEAL (%s), not recomputed:"
        % os.path.basename(SEAL))
    say()
    say("  depth  m     uncorrected   undiluted corr.  1/m corrected")
    for d in range(4):
        say("  %-6d %-5d %-13.6f %-16.6f %.6f"
            % (d, ms[d], unc[d], full[d], corr[d]))
    say()
    say("T1's interval at each depth is [uncorrected - undiluted,")
    say("uncorrected], since the deficit is negative:")
    say()
    say("  depth  lower (fully corrected)  upper (uncorrected)   width")
    for d in range(3):
        lo = unc[d] - full[d]
        say("  %-6d %-24.6f %-21.6f %.6f" % (d, lo, unc[d], full[d]))
    say()
    say("T2's target is the 1/m corrected column, tolerance 1e-4.")
    say("T3 reads alpha from  measured = uncorrected - full / m^alpha,")
    say("i.e.  alpha = log(full / (uncorrected - measured)) / log m.")
    say()
    sm = read_sampled_mean()
    if sm is not None:
        import math
        al = math.log(full[2] / (unc[2] - sm)) / math.log(ms[2])             if unc[2] > sm else float("nan")
        say("The first run's SAMPLED thirty-draw mean at depth 2 is READ")
        say("from lab_secondcell.txt as %.6f -- not retyped." % sm)
        say("It sits %+.6f above T1's upper bound, so the sampled evidence"
            % (sm - unc[2]))
        say("leans towards T1 FAILING at depth 2, and T3's alpha there is")
        say("%s.  Recorded before the exact run so that a pass cannot"
            % ("undefined -- the measured value is above the main term"
               if sm >= unc[2] else "%.2f" % al))
        say("later be read as having been expected.")
    say()
    say("NULL, enumerated here and not declined.  T1 asks whether the")
    say("exact value lands in a cell-specific interval.  If the three")
    say("intervals were interchangeable, T1 would carry no information")
    say("about WHICH cell gets which number.  With three depths the null")
    say("is enumerable -- all 6 relabellings:")
    say()
    from itertools import permutations
    iv = [(unc[d] - full[d], unc[d]) for d in range(3)]
    ok = 0
    for pm in permutations(range(3)):
        hit = sum(1 for d in range(3)
                  if iv[pm[d]][0] <= (unc[d] + iv[d][0]) / 2 <= iv[pm[d]][1])
        if hit == 3:
            ok += 1
    say("    relabellings whose intervals all still contain the true")
    say("    cell's midpoint: %d of 6." % ok)
    say("    The three intervals are disjoint by construction -- they sit")
    say("    at %.3f, %.3f and %.3f -- so a mislabelling is refused by"
        % (unc[0], unc[1], unc[2]))
    say("    size alone.  T1 is not a test of the labelling; it is a test")
    say("    of a %.1e-wide window against an exact number, and that is"
        % max(full[d] for d in range(3)))
    say("    what makes it sharp.")
    say()
    say("=" * 70)
    say("NOTHING IS SCORED HERE.  T1, T2 and T3 need lab_secondcell2.py,")
    say("which does not exist.  That is the point of this file.")

    head = [
        "STATISTIC: the exact pair average of S_2 over each cell of the",
        "           second shape, by autocorrelation of the cell",
        "           indicator -- no sampling, no seed, no standard",
        "           error.  Registered against predictions sealed in",
        "           lab_secondcell_predict.txt at commit 3b1b7d1.",
        "DENOM: the exact band average 2.",
        "NULL: an exact enumeration, run in this file. The statistic",
        "      being measured is exact, so the tolerances are absolute;",
        "      the null asks instead whether T1's intervals are",
        "      cell-specific, by enumerating all 6 relabellings of the",
        "      three depths.",
        "FIELD: predictions read from the first seal; m = 48, 44, 12, 1;",
        "       no measurement is made here.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())

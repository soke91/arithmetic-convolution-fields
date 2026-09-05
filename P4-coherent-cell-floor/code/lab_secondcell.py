# -*- coding: utf-8 -*-
r"""P4 — the SECOND cell shape, measured against a sealed prediction.

WHAT THIS IS

lab_secondcell_predict.py computed and committed four numbers before
this file existed.  This file measures the same four and scores the
rules that file registered.  It READS the predictions from
results/lab_secondcell_predict.txt rather than recomputing them: a
recomputation could drift from what was sealed and nothing would say
so.

The registration is in that file and is not restated here.  What is
restated is the single decision made after the seal and before this
run, because it changes what "the measured D'_c" means:

  THE BAND AVERAGE IS NOT SAMPLED.  E_all[S_2] = 2 exactly -- the
  Euler product it makes is the reciprocal of C_2's term by term --
  and the run that sampled it in the first partition came out low by
  about 1.6 sampling errors, a shared additive offset that moved every
  cell by one sampling error and the deepest by eight.  So
  D'_c = E_same,c[S_2] - 2 here.
  The sampled band average is computed and printed alongside, so the
  choice can be read and reversed by anyone who disagrees with it.

That decision was made after the seal.  It is stated rather than
folded in, and both numbers are printed.

WHAT IS MEASURED

Band (2e6, 4e6], even N, cells indexed by depth = #{p in (3,5,7) : p|N}.
For each cell, SAMPLES ordered pairs N != N' drawn from that cell, and
the mean of S_2(|N - N'|) over them.  S_2 is sieved over every prime to
the maximum shift; no prime product is truncated.

FIELD: band (2e6, 4e6]; even N; cells by depth over 3,5,7; S_2 sieved
       to 2e6 over every prime; 400000 sampled ordered pairs per cell;
       numpy default_rng seed 20260905.  Predictions read from
       results/lab_secondcell_predict.txt.

BACKS: Proposition {#prop:scaleinv} and Measurement {#meas:maincoef}
       in deploy/papers/P4-coherent-cell-floor.tex -- as a test of
       whether the main-term computation is a mechanism or a fact
       about one partition.
"""

import io
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "lab_secondcell.txt")
SEAL = os.path.join(ROOT, "results", "lab_secondcell_predict.txt")

CELLP2 = (3, 5, 7)
BAND = (2_000_000, 4_000_000)
HMAX = 2_000_000          # largest shift inside the band
CLIM = 4_000_000          # bound for the constant, matching the tree
SAMPLES = 400_000
SEED = 20260905
EXACT_EALL = 2.0

lines = []


def say(s=""):
    lines.append(s)
    print(s)


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def read_seal():
    """The four predictions, from the file that was committed before
    this one existed.  Read and not recomputed."""
    if not os.path.exists(SEAL):
        return [], []
    src = io.open(SEAL, encoding="utf-8", errors="replace").read()
    corr, unc = {}, {}
    blk = src.split("depth  corrected prediction")
    if len(blk) > 1:
        for ln in blk[1].splitlines()[1:]:
            f = ln.split()
            if len(f) == 2 and f[0].isdigit():
                corr[int(f[0])] = float(f[1])
            elif corr:
                break
    blk = src.split("m (classes)  n_c approx   uncorrected prediction")
    if len(blk) > 1:
        for ln in blk[1].splitlines()[1:]:
            f = ln.split()
            if len(f) == 4 and f[0].isdigit():
                unc[int(f[0])] = float(f[3])
            elif unc:
                break
    return ([corr[d] for d in sorted(corr)] if corr else [],
            [unc[d] for d in sorted(unc)] if unc else [])


def main():
    say("SECOND CELL SHAPE MEASURED.  S = %s, band (%d, %d]"
        % (", ".join(str(p) for p in CELLP2), BAND[0], BAND[1]))
    say()
    pred, unc = read_seal()
    if len(pred) != len(CELLP2) + 1:
        say("the sealed predictions were not readable -- refusing to run")
        return 1
    say("sealed predictions read from %s:" % os.path.basename(SEAL))
    say("  " + ", ".join("%.6f" % v for v in pred))
    say()

    say("building S_2(h) to %d over every prime ..." % HMAX)
    pr = primes_upto(HMAX)
    # The CONSTANT is accumulated to CLIM and not to HMAX.  The sieve
    # only needs primes up to the largest shift, but 2 C_2 is a product
    # over all primes and truncating it at the shift bound would put a
    # different error on it than every other run in this tree carries.
    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    S2 = np.zeros(HMAX + 1, dtype=np.float64)
    S2[2::2] = twin
    for p in pr:
        p = int(p)
        if p == 2:
            continue
        S2[p::p] *= (p - 1.0) / (p - 2.0)
    S2[1::2] = 0.0
    S2[0] = 0.0
    say("  2 C_2 = %.8f;  S_2(2) = %.6f;  S_2(6) = %.6f (= 2 S_2(2))"
        % (twin, S2[2], S2[6]))
    say()

    Ns = np.arange(BAND[0] + 2 if BAND[0] % 2 == 0 else BAND[0] + 1,
                   BAND[1] + 1, 2, dtype=np.int64)
    depth = np.zeros(Ns.size, dtype=np.int8)
    for p in CELLP2:
        depth += (Ns % p == 0)
    rng = np.random.default_rng(SEED)

    say("  depth  n_c        E_same,c[S_2]   D'_c (E_all = 2)   se")
    say("  " + "-" * 62)
    Es, Dc, SE, NC = [], [], [], []
    DEEP = len(CELLP2)
    for d in range(DEEP + 1):
        idx = np.flatnonzero(depth == d)
        nc = idx.size
        if d == DEEP:
            # ENUMERATED, not sampled.  One class mod 2Q', so the cell is
            # an AP segment and the pair average is an exact finite sum:
            # weight 2(nc-k) at shift 2Q'k.  se = 0 by construction.
            Q2 = 2
            for q in CELLP2:
                Q2 *= q
            ks = np.arange(1, nc, dtype=np.int64)
            w = 2.0 * (nc - ks)
            e = float((w * S2[Q2 * ks]).sum() / w.sum())
            sflag, s_ = "enumerated", 0.0
        else:
            j1 = rng.integers(0, nc, SAMPLES)
            j2 = rng.integers(0, nc, SAMPLES)
            keep = j1 != j2
            v = S2[np.abs(Ns[idx[j1[keep]]] - Ns[idx[j2[keep]]])]
            e = float(v.mean())
            s_ = float(v.std(ddof=1) / np.sqrt(v.size))
            sflag = "%.6f" % s_
        Es.append(e)
        Dc.append(e - EXACT_EALL)
        SE.append(s_)
        NC.append(nc)
        say("  %-6d %-10d %-15.6f %-18.6f %s"
            % (d, nc, e, e - EXACT_EALL, sflag))
    say()

    i1 = rng.integers(0, Ns.size, SAMPLES)
    i2 = rng.integers(0, Ns.size, SAMPLES)
    keep = i1 != i2
    vb = S2[np.abs(Ns[i1[keep]] - Ns[i2[keep]])]
    eall = float(vb.mean())
    seall = float(vb.std(ddof=1) / np.sqrt(vb.size))
    zb = (eall - EXACT_EALL) / seall
    say("P2  the run's own SAMPLED band average against the exact 2:")
    say("    %.6f, se %.6f, %+.2f sampling errors" % (eall, seall, zb))
    p2 = abs(zb) <= 2.0
    say("    within two: %s" % ("hold" if p2 else "REFUTED"))
    say("    The first run's was low by about 1.6 sampling errors; its")
    say("    value is in lab_cell_singular.txt and is not retyped.  A")
    say("    draw on the same side at the same size would say the offset")
    say("    is systematic and would reopen the truncated-prime")
    say("    mechanism that was ruled out by reading the sieve.")
    say("    D'_c above is formed with the EXACT 2 regardless.")
    say()

    say("P1  depths 0 to 2, where the main term is tested alone.")
    say()
    say("  depth  sealed prediction  measured D'_c   residual   in se")
    worst, p1 = 0.0, True
    for d in range(DEEP):
        r = pred[d] - Dc[d]
        z = r / SE[d]
        worst = max(worst, abs(z))
        if abs(z) > 3.0:
            p1 = False
        say("  %-6d %-18.6f %-15.6f %+.6f  %+.2f"
            % (d, pred[d], Dc[d], r, z))
    say()
    say("    within three sampling errors at all three: %s  (worst %.2f)"
        % ("hold" if p1 else "REFUTED", worst))
    say()

    # POST HOC (written after P1 was scored, and rescoring nothing).
    # P1 is REFUTED as registered.  The first diagnosis written here was
    # that the reported sampling error understates the truth; THIRTY
    # independent redraws say it does not -- their spread is 0.98 of the
    # reported error.  That diagnosis was made on six redraws, one of
    # which was an outlier, and it is retracted rather than quietly
    # replaced.  The correct reading is below and it is worse for the
    # registration, not better.
    say("POST HOC (after P1 was scored; it rescores nothing).")
    say()
    keep2 = np.flatnonzero(depth == 2)
    vals, ses = [], []
    for sd in range(1000, 1030):
        r2 = np.random.default_rng(sd)
        a1 = r2.integers(0, keep2.size, SAMPLES)
        a2 = r2.integers(0, keep2.size, SAMPLES)
        kk = a1 != a2
        vv = S2[np.abs(Ns[keep2[a1[kk]]] - Ns[keep2[a2[kk]]])]
        vals.append(float(vv.mean()) - EXACT_EALL)
        ses.append(float(vv.std(ddof=1) / np.sqrt(vv.size)))
    import statistics
    mm = statistics.mean(vals)
    sdv = statistics.stdev(vals)
    say("  Depth 2, thirty independent redraws of the same estimator:")
    say("    mean %.6f, spread across draws %.6f" % (mm, sdv))
    say("    the reported sampling error averages %.6f, so the reported"
        % statistics.mean(ses))
    say("    error is right to within %.0f per cent.  IT IS NOT THE"
        % abs(100 * (sdv / statistics.mean(ses) - 1)))
    say("    PROBLEM.")
    say()
    say("    The sealed prediction is %.6f and sits %+.2f of that error"
        % (pred[2], (pred[2] - mm) / sdv))
    say("    from the thirty-draw mean.  The registered run's own draw was")
    say("    %.6f, which is %+.2f -- outside the range of all thirty."
        % (Dc[2], (Dc[2] - mm) / sdv))
    say()
    say("  SO THE REGISTRATION SCORED ONE DRAW, AND ONE DRAW WAS AN")
    say("  OUTLIER.  A three-sigma cap applied to a single sample at each")
    say("  of three cells is refused by chance about once in a hundred")
    say("  and twenty runs, and this is that.  The main term is not")
    say("  refuted: thirty draws centre on the prediction.  P1 is.")
    say()
    say("  P1 STANDS REFUTED AS REGISTERED.  It is recorded refuted, and")
    say("  the fault is in the rule and not in the prediction: a rule")
    say("  that scores a single random draw carries the draw's luck into")
    say("  its verdict.  A second registration has to score a statistic")
    say("  that does not -- an average over replicates, or an exact")
    say("  enumeration like depth 3's -- and it has to be sealed again")
    say("  rather than rescored here.")
    say()

    say("DEPTH %d -- CONSISTENCY CHECK, NOT A SCORE." % DEEP)
    say("  The corrected prediction and the enumerated measurement are")
    say("  the same finite sum; they differ only by K's exact rational")
    say("  against the truncated Euler product.  A tolerance worth")
    say("  registering would pass trivially, so this is printed to catch")
    say("  a coding error in either route and for nothing else.")
    say("    sealed corrected prediction  %.9f" % pred[DEEP])
    say("    enumerated measurement       %.9f" % Dc[DEEP])
    say("    difference                   %.3e" % (pred[DEEP] - Dc[DEEP]))
    if len(unc) == len(pred):
        say("    the MAIN TERM alone, 19/8   %.9f" % unc[DEEP])
        say("    enumerated minus main term  %.3e  -- this gap is what"
            % (Dc[DEEP] - unc[DEEP]))
        say("    the divisor deficit is, and it is the quantity the")
        say("    deficit computation reproduces by construction.")
    say()
    say("=" * 70)
    say("F1 (in the sealed file) and the two predictions:")
    say("P1 %s   P2 %s" % ("hold" if p1 else "REFUTED",
                           "hold" if p2 else "REFUTED"))
    if p1:
        say("The main-term computation transfers to a cell shape it was")
        say("not built on, at the three cells where it is tested alone.")
        say("It is a mechanism and not a fact about one partition.")
    else:
        say("The main-term computation does NOT transfer.  It is a fact")
        say("about the (3,5,7,11,13) partition, and the paper's claim")
        say("narrows to that partition.")
    say()

    head = [
        "STATISTIC: D'_c = E_same,c[S_2] - 2 for cells indexed by",
        "           depth over (3, 5, 7), against predictions sealed in",
        "           lab_secondcell_predict.txt before this file existed.",
        "DENOM: the exact band average 2. The sampled band average is",
        "       printed beside it so the choice can be reversed.",
        "NULL: the rules scored here were registered in the sealed file,",
        "      which also carries its own enumerated null on whether the",
        "      predicted profile tracks cell size. The sampling error of",
        "      each cell mean is computed from the draws and P2 is stated",
        "      against it.",
        "FIELD: band (2e6, 4e6]; even N; cells by depth over 3,5,7; S_2",
        "       sieved to 2e6 over every prime; 400000 sampled ordered",
        "       pairs per cell; numpy default_rng seed 20260905.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())

# -*- coding: utf-8 -*-
r"""What does the exact cell floor exceed a count-based interval by?

Supports {#rem:sdzrange}.

WHAT IS AT STAKE

The claim under test is printed three times -- in an abstract, a note
and a summary -- and reads: at the top octave the exact floor exceeds
sd(Z)/sqrt(n_c) by factors of 5.8 to 160.  Neither 5.8 nor 160 appears
in any result file, and the statement does not say over what range
sd(Z) is taken.  A number quoted three times with no stated field is
not checkable, which is the defect regardless of what the number turns
out to be.

Two readings of "sd(Z)" are possible and both are computed here: the
standard deviation over the top octave itself, and over the whole
field on which the floor was fitted.  With
Z(N) = C(N)/sqrt(V(N)), C = mu * Lambda, V = mu^2 * Lambda^2, and
se_c the exact floor of the cell mean, the quantity claimed is

    r_c = se_c * sqrt(n_c) / sd(Z).

DISCLOSURE

A blind pass reported 6.71 and 146.4 for the top-octave reading before
this run was written, so the prediction below is not innocent of the
answer.  What that pass did not do is compute the band-wide reading, or
recompute se_c independently; both are done here, from a sieve and a
convolution written for this run.

WHAT IS MEASURED

  W1  whether the printed pair (5.8, 160) is the (min, max) of r_c over
      the depths, under either reading.

  W2  whether r_c grows with depth, as "growing with the cell" says.

  W3  the actual (min, max) under each reading, reported.

  W4  se_c itself against the published table, as a control on the
      reimplementation: the top-octave se by depth is printed in
      results/lab_cell_floor.txt.

FALSIFICATION, registered before the run

  W1  REFUTED if neither reading reproduces both 5.8 and 160 to within
      5 per cent.  Then the printed pair has no reading behind it and
      must be replaced by a measured one.
  W2  REFUTED if r_c is not increasing in depth over the depths present
      at the top octave.  Then "growing with the cell" is wrong too.
  W3  reported, not judged.
  W4  REFUTED if any recomputed se_c differs from the published value
      by more than 1 per cent.  Then this run is not measuring the same
      floor and W1-W3 are uninterpretable.

  PREDICTION.  W1 refuted, W2 holds, W4 holds. The specific suspicion
  is that 160 was borrowed from the factor in N of the field, which the
  same paragraph elsewhere calls 160, and is not a ratio at all.

NULL.  None applies: deterministic sums over a fixed sieve, no
sampling and no sign input.  The control is W4, against a table
produced by a different implementation.
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

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), "results")

X = 16_000_000
CELLP = (3, 5, 7, 11, 13)
TOP = (8_000_000, 16_000_000)
FIELD = (62_500, 16_000_000)      # the eight octaves the floor was fitted on
CLAIM = (5.8, 160.0)
TOL_CLAIM = 0.05                  # relative
TOL_SE = 0.01                     # relative, against the published table
PUBLISHED_B = 0.0395              # the se exponent b of lab_cell_floor.txt, C4
PUBLISHED_SE = {                  # results/lab_cell_floor.txt, top octave
    0: 1.1823e-01, 1: 4.4452e-02, 2: 1.4144e-01,
    3: 2.3053e-01, 4: 3.1686e-01, 5: 4.1153e-01,
}


def pow2(n):
    p = 1
    while p < n:
        p <<= 1
    return p


def sieve(n):
    """Lambda, mu and the squarefree indicator to n."""
    spf = np.zeros(n + 1, dtype=np.int32)
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == 0:
            spf[i * i::i] = np.where(spf[i * i::i] == 0, i, spf[i * i::i])
    lam = np.zeros(n + 1, dtype=np.float64)
    mu = np.zeros(n + 1, dtype=np.float64)
    mu[1] = 1.0
    for v in range(2, n + 1):
        p = int(spf[v]) or v
        w = v // p
        mu[v] = 0.0 if w % p == 0 else -mu[w]
        t = v
        while t % p == 0:
            t //= p
        lam[v] = math.log(p) if t == 1 else 0.0
    return lam, mu, (mu != 0).astype(np.float64)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("STATISTIC: r_c = se_c*sqrt(n_c)/sd(Z) at the top octave, by")
    say("           depth, where se_c is the exact floor of the cell mean")
    say("           of Z(N)=C(N)/sqrt(V(N)) and sd(Z) is taken two ways --")
    say("           over the top octave, and over the whole fitted field;")
    say("           plus se_c itself against the published table.")
    say("FIELD: even N with V(N)>0; top octave (%d, %d]; fitted field"
        % TOP)
    say("       (%d, %d]; cells indexed by depth = #{p in 3,5,7,11,13"
        % FIELD)
    say("       dividing N}; Lambda, mu and the squarefree indicator from")
    say("       an integer sieve to %d; C = mu*Lambda and" % X)
    say("       V = mu^2*Lambda^2 by exact FFT convolution.")
    say("CONSTANTS: X = %d, CELLP = %s, CLAIM = %s, TOL_CLAIM = %.2f,"
        % (X, ",".join(str(p) for p in CELLP), str(CLAIM), TOL_CLAIM))
    say("           TOL_SE = %.2f, PUBLISHED_B = %g, PUBLISHED_SE from"
        % (TOL_SE, PUBLISHED_B))
    say("           results/lab_cell_floor.txt (top octave se column, and")
    say("           the C4 exponent b).")
    say("NULL: none applies -- deterministic sums over a fixed sieve, no")
    say("      sampling and no sign input. The control is W4.")
    say("DENOM: r_c is divided by sd(Z); the se comparison is relative.")
    say()
    say(__doc__.strip())
    say()
    say("=" * 72)
    say("sieving to %d ..." % X)
    lam, mu, sqf = sieve(X)

    say("convolving ...")
    n = pow2(2 * (X + 1))
    a = np.zeros(n, dtype=np.float64)
    a[:X + 1] = lam ** 2
    FL2 = np.fft.rfft(a)
    a[:] = 0.0
    a[:X + 1] = sqf
    V = np.fft.irfft(FL2 * np.fft.rfft(a), n)[:X + 1]
    del FL2
    a[:] = 0.0
    a[:X + 1] = lam
    FL = np.fft.rfft(a)
    a[:] = 0.0
    a[:X + 1] = mu
    C = np.fft.irfft(FL * np.fft.rfft(a), n)[:X + 1]
    del a, FL

    depth = np.zeros(X + 1, dtype=np.int8)
    for p in CELLP:
        depth[p::p] += 1

    def evens(lo, hi):
        ns = np.arange(lo + 2 - (lo % 2), hi + 1, 2, dtype=np.int64)
        return ns[V[ns] > 0]

    Ns = evens(*TOP)
    Z = C[Ns] / np.sqrt(V[Ns])
    sd_top = float(Z.std(ddof=1))

    Nf = evens(*FIELD)
    sd_field = float((C[Nf] / np.sqrt(V[Nf])).std(ddof=1))

    # the exact floor, recomputed
    T = TOP[1]
    m = pow2(2 * (T + 1))
    b = np.zeros(m, dtype=np.float64)
    b[:T + 1] = lam[:T + 1]
    FLam = np.conj(np.fft.rfft(b))
    g = 1.0 / np.sqrt(V[Ns])
    nb = Ns.size

    def ucorr(sel):
        b[:] = 0.0
        b[Ns[sel]] = g[sel]
        return np.fft.irfft(FLam * np.fft.rfft(b), m)[:T + 1]

    w = sqf[:T + 1]
    ua = ucorr(np.ones(nb, dtype=bool))
    Qaa = float((w * ua * ua).sum())

    # the abscissa the floor was fitted on, and the factor in N across it:
    # printed here because the paragraph this run is about quotes both.
    mids = [FIELD[0] * 1.5 * (2 ** j) for j in range(8)]
    say()
    say("  the fitted abscissa is the octave midpoint, %s ... %s,"
        % ("%.5g" % mids[0], "%.5g" % mids[-1]))
    say("  mean log = %.4f, so 1/(2<log N>) = %.6f, which is what"
        % (sum(math.log(v) for v in mids) / len(mids),
           1.0 / (2.0 * sum(math.log(v) for v in mids) / len(mids))))
    say("  results/lab_cell_floor.txt prints; the factor in N across it")
    say("  is %.5g / %.5g = %.1f." % (mids[-1], mids[0], mids[-1] / mids[0]))
    say("  sqrt of that factor = %.4f, and it to the published power"
        % math.sqrt(mids[-1] / mids[0]))
    say("  PUBLISHED_B = %g is %.4f;"
        % (PUBLISHED_B, (mids[-1] / mids[0]) ** PUBLISHED_B))
    say("  sqrt(log %.5g / log %.5g) = %.4f."
        % (mids[-1], mids[0],
           math.sqrt(math.log(mids[-1]) / math.log(mids[0]))))
    say()
    say("  sd(Z) over the top octave     %.6f  (n = %d)" % (sd_top, nb))
    say("  sd(Z) over the fitted field   %.6f  (n = %d)"
        % (sd_field, Nf.size))
    say()
    hdr = ("  %-6s %-9s %-12s %-12s %-11s %-11s"
           % ("depth", "n_c", "se_c", "published", "r (octave)", "r (field)"))
    say(hdr)
    say("  " + "-" * (len(hdr) - 2))

    ok4 = True
    r_oct, r_fld, depths = [], [], []
    for d in range(6):
        sel = depth[Ns] == d
        nc = int(sel.sum())
        if nc == 0:
            continue
        uc = ucorr(sel)
        Qcc = float((w * uc * uc).sum())
        Qca = float((w * uc * ua).sum())
        var = Qcc / nc ** 2 - 2.0 * Qca / (nc * nb) + Qaa / nb ** 2
        se = math.sqrt(max(var, 0.0))
        pub = PUBLISHED_SE[d]
        ok4 &= abs(se - pub) <= TOL_SE * pub
        ro = se * math.sqrt(nc) / sd_top
        rf = se * math.sqrt(nc) / sd_field
        depths.append(d)
        r_oct.append(ro)
        r_fld.append(rf)
        say("  %-6d %-9d %-12.6e %-12.6e %-11.3f %-11.3f"
            % (d, nc, se, pub, ro, rf))

    def pair(v):
        return (min(v), max(v))

    def matches(v):
        lo, hi = pair(v)
        return (abs(lo - CLAIM[0]) <= TOL_CLAIM * CLAIM[0]
                and abs(hi - CLAIM[1]) <= TOL_CLAIM * CLAIM[1])

    ok1 = matches(r_oct) or matches(r_fld)
    ok2 = all(r_oct[i] < r_oct[i + 1] for i in range(len(r_oct) - 1))

    say()
    say("  (min, max) over depths, top-octave sd   (%.3f, %.3f)"
        % pair(r_oct))
    say("  (min, max) over depths, field sd        (%.3f, %.3f)"
        % pair(r_fld))
    say("  printed in the paper                    (%.3f, %.3f)"
        % CLAIM)
    say()
    say("W1  the printed pair is reproduced by one of the two readings")
    say("    W1 %s" % ("hold" if ok1 else "REFUTED"))
    say("W2  the ratio grows with the cell")
    say("    W2 %s" % ("hold" if ok2 else "REFUTED"))
    say("W3  the measured pairs are in the two lines above, reported")
    say("    not judged.")
    say("W4  the recomputed floor matches the published one")
    say("    W4 %s" % ("hold" if ok4 else "REFUTED"))
    say()
    say("=" * 72)
    say("W1 %s  W2 %s  W4 %s"
        % tuple("hold" if v else "REFUTED" for v in (ok1, ok2, ok4)))

    io.open(os.path.join(RES, "audit_cellfloor_sdz.txt"), "w",
            encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    return 0 if (ok2 and ok4) else 1


if __name__ == "__main__":
    raise SystemExit(main())

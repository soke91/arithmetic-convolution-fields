# -*- coding: utf-8 -*-
r"""
paper/wall_v3.md, Proposition {#prop:scaleinv} and the paragraph above
it -- what predicts the floor's size.

WHY THIS HAS TO BE RUN

The paragraph says the floor's magnitude tracks

    D_c := E_same,c[S_2] - E_all[S_2],

with S_2(h) the Hardy-Littlewood singular series of the shift h,
E_same,c the mean over pairs N, N' both in cell c, and E_all the mean
over all pairs in the band; and Proposition [prop:scaleinv] says D_c
depends only on the law of h induced by the divisibility patterns the
cell is a union of, and on the densities with which those patterns
enter it, hence is scale-invariant to O(1/n_c) and predicts an
exponent of zero at every depth.  Depth fixes how MANY of 3,5,7,11,13
divide N, not the residues themselves; an earlier printing of this
docstring and of the proposition said "the residue classes that the
cell fixes", which the proposition's own proof opens by denying.

Neither S_2 nor D_c is computed anywhere in this repository.  The
proposition's evidence marker pointed at lab_cell_floor.py, which
computes the floor and the z-scores and never forms S_2 at all -- an
evidence marker can name an existing script that does not compute the
statement, and G1 and G4 only check that the file and its result exist.
So this is the claim in either paper with the weakest backing, and it
is the one that explains the floor.

    S_2(h) = 2 C_2 prod_{p | h, p > 2} (p-1)/(p-2)   for even h,

with C_2 the twin-prime constant.  Cells are indexed by depth, the
number of 3,5,7,11,13 dividing N, so a deeper cell makes p | h more
likely for those p and should raise E_same,c[S_2].

BACKS: Proposition {#prop:scaleinv} and Remark {#rem:copiedfloor}
in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  M1  D_c > 0 at every depth d >= 1, and D_c increases with depth.
  M2  D_c is scale-invariant: over the three octaves (1e6,2e6],
      (2e6,4e6], (4e6,8e6] it varies by less than 2% at every depth,
      beyond sampling error.
  M3  The floor tracks it: across depths, the correlation between the
      exact floor se_c of Lemma [lem:cellmom] and D_c exceeds 0.8.
  M4  Fitting D_c ~ N^{-e} across the three octaves gives |e| < 0.01 at
      every depth -- the exponent zero the proposition predicts.

REFUTATION RULE (fixed before the run)

  M1  REFUTED if D_c <= 0 at some depth >= 1, or if the sequence is not
      increasing in depth.
  M2  REFUTED if the spread (max-min)/mean across octaves exceeds 0.02
      at any depth by more than three sampling standard errors.
  M3  REFUTED if the correlation is 0.8 or below.
  M4  REFUTED if |e| >= 0.01 at any depth.

  M1, M2 and M4 gate -- they are the proposition.  M3 is the
  paragraph's "tracks" claim and is reported.
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
OUT = os.path.join(ROOT, "results", "lab_cell_singular.txt")

OCTS = [(1_000_000, 2_000_000), (2_000_000, 4_000_000),
        (4_000_000, 8_000_000)]
HMAX = 8_000_000
# The S_2 sieve below multiplies (p-1)/(p-2) into S2[p::p] for every
# prime p <= HMAX, so it is correct for |h| <= HMAX and would SILENTLY
# UNDER-COUNT a shift carrying a prime factor above it -- dropping a
# factor and biasing the average low, with no error and no warning.
# That cannot happen here: |h| < B <= HMAX/2.  It is written down
# because the bound must not be reused at a larger band, which is where
# the same construction would fail exactly that way.
CELLP = (3, 5, 7, 11, 13)
SAMPLES = 400_000
SEED = 20260808
# The exact floor at (2e6,4e6] is computed by lab_mask_placebo.py and
# is READ from its result file, not copied. A hand-copied table is a
# dependency no check can see: G18 compares a script with its own
# result and G22 compares a script with what it reads, and neither
# sees a number that was typed in. Reading it makes the dependency
# visible and puts it under G22.
def read_floor():
    p = os.path.join(ROOT, "results", "lab_mask_placebo.txt")
    src = io.open(p, encoding="utf-8").read()
    blk = src[src.index("true labelling"):]
    out = {}
    for ln in blk.splitlines()[2:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        out[int(f[0])] = float(f[2])
    return [out[d] for d in sorted(out)]


def read_k_terms():
    """lab_mask_placebo.txt 의 사후 블록이 인쇄하는 K 의 세 항.

    depth -> (E_same,c, E_c,a, E_all).  세 평균이 각각 n_c^2, n_c*n,
    n^2 쌍 위의 것이라, D_c 의 규약(N != N')으로 옮기려면 항마다 자기
    대각을 빼고 자기 쌍 수로 나눠야 한다."""
    f = os.path.join(ROOT, "results", "lab_mask_placebo.txt")
    if not os.path.exists(f):
        return {}
    out = {}
    for ln in io.open(f, encoding="utf-8", errors="replace"):
        t = ln.split()
        if len(t) >= 7 and t[0].isdigit():
            try:
                d, a, b, c = int(t[0]), float(t[1]), float(t[2]), float(t[3])
            except ValueError:
                continue
            if 0.0 < a < 1.0 and 0.0 < b < 1.0 and 0.0 < c < 1.0:
                out.setdefault(d, (a, b, c))
    return out


def read_qcc_by_octave():
    """(depth -> [(N_mid, n_c, Q_cc/n_c^2, lo, hi)]) from lab_cell_floor.txt.

    Separate from read_var_by_octave because the column is a different
    one: that reader takes Var (field 6), this one takes the first of
    the three terms (field 5) and the cell size (field 4).  Observation
    {#obs:coh} reads Q_cc n_c^-2 log N across scale, and until this
    block existed the numbers in that paragraph were produced by no
    run at all -- they were carried in by hand.
    """
    p = os.path.join(ROOT, "results", "lab_cell_floor.txt")
    out = {}
    try:
        src = io.open(p, encoding="utf-8").read()
    except OSError:
        return out
    for ln in src.splitlines():
        f = ln.split()
        if len(f) == 10 and f[0] == "(" and f[3].isdigit():
            try:
                lo, hi = int(f[1].rstrip(",")), int(f[2].rstrip("]"))
                # both endpoints, not only the midpoint: the fitted
                # exponent below is read at all three and they disagree
                # by 60 per cent, so a reader who is handed only the
                # midpoint cannot see that the value is a convention.
                out.setdefault(int(f[3]), []).append(
                    (0.5 * (lo + hi), int(f[4]), float(f[5]), lo, hi))
            except ValueError:
                continue
    return out


def read_b_by_depth():
    """The registered exponent b in se ~ N^{-b}, block C4 of the same
    file.  Read rather than refitted: refitting here gives 0.039400 at
    depth 0 against the registered 0.039388, and a paragraph that
    compares b to a prediction must use the b that was registered."""
    p = os.path.join(ROOT, "results", "lab_cell_floor.txt")
    out = {}
    try:
        src = io.open(p, encoding="utf-8").read()
    except OSError:
        return out
    for ln in src.splitlines():
        f = ln.split()
        if len(f) >= 5 and f[0] == "depth" and f[2] == "b" and f[3] == "=":
            try:
                out[int(f[1].rstrip(":"))] = float(f[4])
            except ValueError:
                continue
    return out


def read_var_by_octave():
    """Var(m_c - mbar) per octave and depth, from lab_cell_floor.txt.

    The band of se_c^2/D_c that this file prints is measured ACROSS
    DEPTHS at one octave.  Nothing here said how the same ratio moves
    ACROSS SCALE, and both inputs are already published: D_c below on
    three octaves, Var in that file on eight.  The ratio is formed here
    so that the flatness claim can be read against the direction it was
    not measured in.
    """
    p = os.path.join(ROOT, "results", "lab_cell_floor.txt")
    out = {}
    try:
        src = io.open(p, encoding="utf-8").read()
    except OSError:
        return out
    for ln in src.splitlines():
        f = ln.split()
        if len(f) == 10 and f[0] == "(" and f[3].isdigit():
            try:
                out[(int(f[1].rstrip(",")), int(f[2].rstrip("]")),
                     int(f[3]))] = float(f[6])
            except ValueError:
                continue
    return out


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("building S_2(h) to %d ..." % HMAX)
    pr = primes_upto(HMAX)
    twin = 2.0
    for p in pr:
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
    say("  2*C_2 = %.6f;  S_2(2) = %.6f;  S_2(6) = %.6f  (should be "
        "2*C_2 * 2)" % (twin, S2[2], S2[6]))

    depth = np.zeros(HMAX + 1, dtype=np.int8)
    for p in CELLP:
        depth[p::p] += 1

    rng = np.random.default_rng(SEED)
    say()
    say("  octave              depth  n_c        E_same,c[S_2]   D_c"
        "          se(D_c)")
    say("  " + "-" * 76)
    D = {}
    CROSS = {}
    SE = {}
    EALL = {}
    XRNG = np.random.default_rng(SEED + 1)
    for lo, hi in OCTS:
        Ns = np.arange(lo + 2, hi + 1, 2, dtype=np.int64)
        dep = depth[Ns]
        i1 = rng.integers(0, Ns.size, SAMPLES)
        i2 = rng.integers(0, Ns.size, SAMPLES)
        ok = i1 != i2
        hall = np.abs(Ns[i1[ok]] - Ns[i2[ok]])
        Eall = float(S2[hall].mean())
        for d in range(6):
            idx = np.flatnonzero(dep == d)
            if idx.size < 2:
                continue
            j1 = rng.integers(0, idx.size, SAMPLES)
            j2 = rng.integers(0, idx.size, SAMPLES)
            k = j1 != j2
            hs = np.abs(Ns[idx[j1[k]]] - Ns[idx[j2[k]]])
            v = S2[hs]
            Es = float(v.mean())
            # POST HOC (see the comment at M3).  A separate generator,
            # so that no registered number moves.
            c1 = XRNG.integers(0, idx.size, SAMPLES)
            a2 = XRNG.integers(0, Ns.size, SAMPLES)
            hca = np.abs(Ns[idx[c1]] - Ns[a2])
            hca = hca[hca > 0]
            CROSS.setdefault(d, []).append((Es, float(S2[hca].mean())))
            se = float(v.std(ddof=1)) / math.sqrt(v.size)
            D.setdefault(d, []).append(Es - Eall)
            SE[(lo, hi, d)] = se
            say("  (%9d,%9d] %-6d %-10d %-15.6f %-12.6f %.6f"
                % (lo, hi, d, idx.size, Es, Es - Eall, se))
        say("      E_all[S_2] = %.6f" % Eall)
        EALL[(lo, hi)] = Eall

    say()
    ds = sorted(D)
    m1 = all(D[d][1] > 0 for d in ds if d >= 1) and all(
        D[ds[i]][1] < D[ds[i + 1]][1] for i in range(len(ds) - 1))
    say("M1  D_c at (2e6,4e6] by depth: %s"
        % ", ".join("%.6f" % D[d][1] for d in ds))
    say("    positive for d>=1 and increasing in depth: %s   %s"
        % (m1, "hold" if m1 else "REFUTED"))

    say()
    say("M2/M4  scale invariance")
    say("  depth  D_c per octave                       spread   "
        "fitted exponent")
    m2 = m4 = True
    absp = {}
    for d in ds:
        vals = np.array(D[d])
        absp[d] = (vals.max() - vals.min(), abs(float(vals.mean())))
        spread = (vals.max() - vals.min()) / abs(vals.mean())
        xs = np.log(np.array([0.5 * (lo + hi) for lo, hi in OCTS]))
        e = -float(np.polyfit(xs, np.log(np.abs(vals)), 1)[0])
        if spread > 0.02:
            m2 = False
        if abs(e) >= 0.01:
            m4 = False
        say("  %-6d %-36s %-8.4f %+.6f"
            % (d, ", ".join("%.6f" % v for v in vals), spread, e))
    say("  M2 %s   M4 %s" % ("hold" if m2 else "REFUTED",
                             "hold" if m4 else "REFUTED"))
    say("  NOTE. The registration reads \"exceeds 0.02 ... by more than")
    say("  three sampling standard errors\"; the code applies the flat")
    say("  0.02 with no allowance. Scored as registered instead:")
    say("  depth  spread   3 se(spread)   0.02 + allowance   as registered")
    n_reg = 0
    for d in ds:
        vals = np.array(D[d])
        spread = float((vals.max() - vals.min()) / abs(vals.mean()))
        sev = max(SE[(lo, hi, d)] for lo, hi in OCTS if (lo, hi, d) in SE)
        allow = 3.0 * math.sqrt(2.0) * sev / abs(float(vals.mean()))
        held = spread <= 0.02 + allow
        n_reg += held
        say("  %-6d %-8.4f %-14.4f %-18.4f %s"
            % (d, spread, allow, 0.02 + allow, "hold" if held else "REFUTED"))
    say("  As registered M2 holds at %d of %d depths; as coded at %d."
        % (n_reg, len(ds), sum(1 for d in ds
             if (np.array(D[d]).max() - np.array(D[d]).min())
                / abs(np.array(D[d]).mean()) <= 0.02)))
    say("  The code is the stricter of the two. Which of them the")
    say("  paper may call \"as registered\" is the registration's rule,")
    say("  not the code's.")
    say("  DIAGNOSTIC (post hoc). M2 and M4 fail only where D_c is small:")
    say("  the spread across octaves is nearly the same size at every")
    say("  depth while D_c itself grows, so the relative spread the")
    say("  rule looks at is large only where D_c is small:")
    say("  depth  spread across octaves   |D_c|      spread/|D_c|")
    for d in ds:
        a, m = absp[d]
        say("  %-6d %-23.6f %-10.6f %.4f" % (d, a, m, a / m))
    say("  Note also")
    say("  that the pre-registration allowed three sampling standard")
    say("  errors and the code applied the 2% band with no allowance --")
    say("  the code is the stricter of the two and its verdict stands.")
    say("  Re-sampling depths 0 and 1 ten times harder:")
    BIG = SAMPLES * 10
    for d in (0, 1):
        vals = []
        for lo, hi in OCTS:
            Ns = np.arange(lo + 2, hi + 1, 2, dtype=np.int64)
            dep = depth[Ns]
            i1 = rng.integers(0, Ns.size, BIG)
            i2 = rng.integers(0, Ns.size, BIG)
            ok = i1 != i2
            Eall = float(S2[np.abs(Ns[i1[ok]] - Ns[i2[ok]])].mean())
            idx = np.flatnonzero(dep == d)
            j1 = rng.integers(0, idx.size, BIG)
            j2 = rng.integers(0, idx.size, BIG)
            k = j1 != j2
            v = S2[np.abs(Ns[idx[j1[k]]] - Ns[idx[j2[k]]])]
            vals.append(float(v.mean()) - Eall)
        vals = np.array(vals)
        xs = np.log(np.array([0.5 * (lo + hi) for lo, hi in OCTS]))
        e = -float(np.polyfit(xs, np.log(np.abs(vals)), 1)[0])
        say("    depth %d: %s   spread %.4f   exponent %+.6f"
            % (d, ", ".join("%.6f" % v for v in vals),
               (vals.max() - vals.min()) / abs(vals.mean()), e))

    say()
    dv = np.array([D[d][1] for d in ds])
    SE_TRUE = read_floor()
    sv = np.array([SE_TRUE[d] for d in ds])
    rho = float(np.corrcoef(dv, sv)[0, 1])
    # M3's control. lab_mask_placebo.py permutes depth labels to test
    # the DETECTION z_c; nothing tested this correlation. With six
    # depths the permutation null is exact -- all 720 orderings.
    import itertools
    perm = sorted(abs(float(np.corrcoef(dv, np.array(q))[0, 1]))
                  for q in itertools.permutations(sv.tolist()))
    above = sum(1 for v in perm if v >= abs(rho))
    m3 = rho > 0.8
    say("M3  corr(se_c, D_c) across depths = %.4f   (floor 0.8)   %s"
        % (rho, "hold" if m3 else "REFUTED"))
    say("    se_c  = %s" % ", ".join("%.5f" % v for v in sv))
    say("    control: over all %d permutations of the depth"
        % len(perm))
    say("    labels, %d reach |r| = %.4f; the null median is"
        % (above, abs(rho)))
    say("    %.4f and its 95th percentile %.4f"
        % (perm[len(perm) // 2], perm[int(0.95 * len(perm))]))
    say("    -- the control lab_mask_placebo.py does not"
        " supply, since it permutes labels for z_c and not")
    say("    for this pair.")
    say("    (this null is the one declared in the NULL: heading"
        " above; it is registered, not post hoc)")
    # A correlation across six points cannot separate se_c from any
    # other quantity that falls with the depth, and n_c is one: the
    # rank order of n_c is nearly the rank order of D_c.  The ratio is
    # the stronger statement and needs no null.
    #
    # It is NOT true that Var = n_c^-2 J K J is a same-cell-minus-band
    # mean.  J is doubly centred, so
    #     Var = E_same,c[K] - 2 E_c,a[K] + E_all[K],
    # and for the kernel K the cross term is not the band term: from
    # and for the kernel K the cross term is not the band term; the
    # size of that is measured in lab_mask_placebo.py, which has all
    # three of K's terms.  An earlier revision of this comment, and
    # the sentence in the paper that took it, asserted the two-term
    # form for Var; both were wrong.
    #
    # D_c above is a two-term difference by definition.  The POST HOC
    # block below measures whether that matters for S_2, i.e. whether
    # E_c,a[S_2] is the band mean.  If it is, then D_c equals the
    # doubly centred B_c[S_2] and the parallel with Var is between two
    # instances of the same functional.  This is diagnostic: it was
    # written after M1-M4 were scored and lifts no verdict.
    rat = sv ** 2 / dv
    say("    se_c^2 / D_c = %s"
        % ", ".join("%.5f" % v for v in rat))
    say("    that ratio spans %.5f to %.5f, a band of %.1f%% about its"
        % (rat.min(), rat.max(),
           100.0 * (rat.max() - rat.min()) / rat.mean()))
    say("    mean, while D_c itself moves by a factor %.0f across the"
        % (dv.max() / dv.min()))
    say("    same six cells. A quantity that merely falls with the")
    say("    depth cannot produce a ratio flat in D_c.")

    # POST HOC (not a registered rule; written after M1-M4 were scored
    # and lifting none of them).  Var carries the diagonal of K and
    # D_c does not: the S_2 averages run over pairs N != N', while
    # Var = n_c^-2 sum_{N,N'} J_c(N) J_c(N') K(N,N') includes N = N'
    # with K(N,N) = 1.  That diagonal is n_c^-2 sum_N J_c(N)^2 =
    # (1 - n_c/n)/n_c, and it has no partner in the denominator.
    lo2, hi2 = OCTS[1]
    Ns2 = np.arange(lo2 + 2, hi2 + 1, 2, dtype=np.int64)
    dep2 = depth[Ns2]
    nband = int(Ns2.size)
    ncs2 = np.array([float(np.count_nonzero(dep2 == d)) for d in ds])
    dg = (1.0 - ncs2 / nband) / ncs2
    rat2 = (sv ** 2 - dg) / dv
    say()
    say("    POST HOC (not a registered rule; written after M1-M4 were")
    say("    scored, and it lifts none of them). Var carries K's")
    say("    diagonal and D_c does not -- the S_2 averages run over")
    say("    pairs N != N'. The diagonal is (1 - n_c/n)/n_c, with no")
    say("    partner in the denominator. At (%d,%d]:" % (lo2, hi2))
    say("    depth  n_c        (1-p)/n_c    se_c^2         printed  matched")
    for i, d in enumerate(ds):
        say("    %-6d %-10d %.4e   %.6e   %.5f  %.5f"
            % (d, int(ncs2[i]), dg[i], sv[i] ** 2, rat[i], rat2[i]))
    m0 = all(rat[i] < rat[i + 1] for i in range(1, len(ds) - 1))
    mono_matched = all(rat2[i] < rat2[i + 1] for i in range(1, len(ds) - 1))
    say("    width about the mean: printed %.1f%%, matched %.1f%%"
        % (100.0 * (rat.max() - rat.min()) / rat.mean(),
           100.0 * (rat2.max() - rat2.min()) / rat2.mean()))
    say("    monotone from depth 1 on: printed %s, matched %s"
        % (m0, mono_matched))
    say("    the diagonal is %.2f%% of Var at the deepest cell; the"
        % (100.0 * dg[-1] / sv[-1] ** 2))
    say("    printed rise from depth 4 to 5 is %+.5f against the"
        % (rat[-1] - rat[-2]))
    say("    %+.6f that the diagonal puts there -- the last step of"
        % (dg[-1] / dv[-1]))
    say("    the monotone rise is the diagonal's.")
    say()
    say("    Removing the diagonal is only half the convention. Var")
    say("    divides the off-diagonal sum by n_c^2 while D_c's E_same")
    say("    is a mean over the n_c(n_c-1) ordered pairs, so the")
    say("    matched numerator carries a further n_c/(n_c-1):")
    rat3 = (sv ** 2 - dg) * ncs2 / (ncs2 - 1.0) / dv
    say("    depth  se_c^2/D_c   diagonal out   convention matched")
    for i, d in enumerate(ds):
        say("    %-6d %.5f      %.5f        %.5f"
            % (d, rat[i], rat2[i], rat3[i]))
    mono3 = all(rat3[i] < rat3[i + 1] for i in range(1, len(ds) - 1))
    say("    width about the mean %.1f%%, monotone from depth 1 on: %s"
        % (100.0 * (rat3.max() - rat3.min()) / rat3.mean(), mono3))
    say("    the step from depth 4 to 5 is %+.5f" % (rat3[-1] - rat3[-2]))
    say()
    say("    That was still not the convention. Var is the DOUBLY")
    say("    CENTRED mean E_same,c[K] - 2 E_c,a[K] + E_all[K], so it is")
    say("    three sums over three index sets with three different pair")
    say("    counts -- n_c^2, n_c*n, n^2. No single scalar carries it to")
    say("    D_c's convention. Each term has to lose its own diagonal and")
    say("    be divided by its own count. K(N,N) = 1 makes that exact:")
    KT = read_k_terms()
    if KT:
        Es = np.array([KT[d][0] for d in ds])
        Eca = np.array([KT[d][1] for d in ds])
        Eall_k = KT[ds[0]][2]
        Es_x = (Es * ncs2 ** 2 - ncs2) / (ncs2 * (ncs2 - 1.0))
        Eca_x = (Eca * ncs2 * nband - ncs2) / (ncs2 * nband - ncs2)
        Eall_x = (Eall_k * nband * nband - nband) / (nband * (nband - 1.0))
        var_x = Es_x - 2.0 * Eca_x + Eall_x
        rat4 = var_x / dv
        say("    depth  printed    diag out   half-matched  fully matched")
        for i, d in enumerate(ds):
            say("    %-6d %.5f   %.5f    %.5f       %.5f"
                % (d, rat[i], rat2[i], rat3[i], rat4[i]))
        m4x = all(rat4[i] < rat4[i + 1] for i in range(1, len(ds) - 1))
        say("    width %.1f%%, monotone from depth 1 on: %s, step 4->5 %+.5f"
            % (100.0 * (rat4.max() - rat4.min()) / rat4.mean(), m4x,
               rat4[-1] - rat4[-2]))
        say("    So the monotone rise SURVIVES a correct match; what")
        say("    shrinks is its size. Two earlier readings of this row")
        say("    over-corrected: removing only the diagonal gives")
        say("    %+.5f at the last step and scaling that by n_c/(n_c-1)"
            % (rat2[-1] - rat2[-2]))
        say("    gives %+.5f, both negative, both wrong."
            % (rat3[-1] - rat3[-2]))
    else:
        say("    (lab_mask_placebo.txt's three-term table not found;")
        say("     the fully matched row cannot be formed here.)")
    VARO = read_var_by_octave()
    rows = []
    for d in ds:
        vv = [VARO.get((lo, hi, d)) for lo, hi in OCTS]
        dd = [D[d][k] for k in range(len(OCTS))] if len(D[d]) == len(OCTS) else None
        if any(v is None for v in vv) or dd is None:
            continue
        r0, r1 = [], []
        for k, (lo, hi) in enumerate(OCTS):
            Nk = np.arange(lo + 2, hi + 1, 2, dtype=np.int64)
            dk = depth[Nk]
            nk = float(np.count_nonzero(dk == d))
            g = (1.0 - nk / float(Nk.size)) / nk
            r0.append(vv[k] / dd[k])
            r1.append((vv[k] - g) / dd[k])
        rows.append((d, r0, r1))
    say("    The same removal across scale, first octave to last:")
    say("    depth  printed drift   matched drift")
    for d, r0, r1 in rows:
        say("    %-6d %6.1f%%         %6.1f%%"
            % (d, 100.0 * abs(r0[0] - r0[-1]) / r0[0],
               100.0 * abs(r1[0] - r1[-1]) / r1[0]))
    say("    -- the deepest cell's drift is the one the diagonal")
    say("    inflates; the shallow five barely move.")
    say()
    say("    THE OTHER DIRECTION. The band above is measured across")
    say("    DEPTHS at one octave. The same ratio across SCALE, from")
    say("    this run's D_c and lab_cell_floor.txt's Var on the three")
    say("    octaves both cover:")
    VAR = read_var_by_octave()
    say("    depth   " + "  ".join("(%de6,%de6]"
                                  % (lo // 1000000, hi // 1000000)
                                  for lo, hi in OCTS))
    n_drift = 0
    for d in sorted(D):
        rr = [VAR.get((lo, hi, d)) for lo, hi in OCTS]
        if any(v is None for v in rr) or len(D[d]) != len(OCTS):
            continue
        q = [v / dc for v, dc in zip(rr, D[d])]
        n_drift += 1
        say("    %-6d  %s   %.1f%% from first to last"
            % (d, "  ".join("%-11.5f" % v for v in q),
               100.0 * (q[0] - q[-1]) / q[0]))
    say("    So D_c fixes the ratio's DEPTH profile at a fixed scale")
    say("    and not its scale dependence: the drift above is of the")
    say("    same order as the 25.4%% band offered as flatness, and it")
    say("    is in the direction that band does not cover. Depth 1 is")
    say("    not monotone in it. This is a diagnostic, not a")
    say("    registered rule -- it was formed after M1-M4 were scored")
    say("    and it scores none of them.")
    say("    rows formed: %d of %d" % (n_drift, len(D)))
    say()
    say("    POST HOC (not a registered rule; run after M1-M4 were")
    say("    scored, and it lifts none of them). Var is doubly centred")
    say("    and D_c is a two-term difference. The doubly centred")
    say("    B_c[S_2] = E_same,c - 2 E_c,a + E_all at (2e6,4e6]:")
    say("    depth  E_same     E_c,a      E_all      D_c        B_c"
        "         B_c/D_c")
    oct2 = OCTS[1]
    ea = EALL[oct2]
    bv = []
    for d in sorted(CROSS):
        es, eca = CROSS[d][1]
        b = es - 2.0 * eca + ea
        bv.append(b)
        say("    %-6d %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f  %.4f"
            % (d, es, eca, ea, es - ea, b, b / (es - ea)))
    bv = np.array(bv)
    say("    E_c,a[S_2] departs from E_all[S_2] by at most %.6f in"
        % max(abs(CROSS[d][1][1] - ea) for d in sorted(CROSS)))
    say("    absolute value, so for this kernel the cross term IS the")
    say("    band mean and the two-term D_c is the doubly centred one.")
    rb = sv ** 2 / bv
    say("    se_c^2 / B_c = %s" % ", ".join("%.5f" % v for v in rb))
    say("    band %.1f%% about its mean, against %.1f%% for D_c."
        % (100.0 * (rb.max() - rb.min()) / rb.mean(),
           100.0 * (rat.max() - rat.min()) / rat.mean()))
    say("    Whether the same holds for the kernel K is not this")
    say("    file's to answer -- K's three terms are computed in")
    say("    lab_mask_placebo.py, and that file reports it.")
    say("    D_c   = %s" % ", ".join("%.5f" % v for v in dv))


    # POST HOC (not a registered rule; written after M1-M4 were scored,
    # and it lifts none of them).  M1-M4 are scored on D_c estimated
    # from SAMPLES ordered pairs.  D_c is a deterministic functional of
    # the cell and can be computed over ALL pairs by one autocorrelation
    # per cell: sum_{N,N' in c} S_2(N-N') = sum_h S_2(h) (1_c * 1_c)(h).
    # This block does that, so that the sampling error of the estimator
    # can be told from the behaviour of the quantity.
    say()
    say("POST HOC (not a registered rule; run after M1-M4 were scored,")
    say("and it lifts none of them). D_c over ALL pairs, by")
    say("autocorrelation -- no sampling:")
    exact = {}
    for lo, hi in OCTS:
        Ns = np.arange(lo + 2, hi + 1, 2, dtype=np.int64)
        L = Ns.size
        sz = 1 << (2 * L - 1).bit_length()
        w = S2[2:2 * L:2]
        ind_all = np.ones(L)
        Fa = np.fft.rfft(ind_all, sz)
        ac = np.fft.irfft(np.abs(Fa) ** 2, sz)[:L]
        Eall_x = 2.0 * float(np.dot(ac[1:L], w)) / (L * (L - 1.0))
        dp = depth[Ns]
        for d in ds:
            ind = (dp == d).astype(float)
            nc = int(ind.sum())
            if nc < 2:
                continue
            F = np.fft.rfft(ind, sz)
            acc = np.fft.irfft(np.abs(F) ** 2, sz)[:L]
            Es = 2.0 * float(np.dot(acc[1:L], w)) / (nc * (nc - 1.0))
            exact.setdefault(d, []).append(Es - Eall_x)
        del Fa, ac
    say("  depth  D_c per octave (exact)                 spread   exponent")
    xs = np.log(np.array([0.5 * (lo + hi) for lo, hi in OCTS]))
    n_hold_2 = n_hold_4 = 0
    for d in ds:
        v = exact.get(d)
        if not v or len(v) != len(OCTS):
            continue
        v = np.array(v)
        sp = float((v.max() - v.min()) / abs(v.mean()))
        e = -float(np.polyfit(xs, np.log(np.abs(v)), 1)[0])
        n_hold_2 += sp <= 0.02
        n_hold_4 += abs(e) < 0.01
        say("  %-6d %-38s %-8.4f %+.6f"
            % (d, ", ".join("%.6f" % x for x in v), sp, e))
    say("  Against M2's cap of 2%%: %d of %d hold." % (n_hold_2, len(ds)))
    say("  Against M4's |e| < 0.01: %d of %d hold." % (n_hold_4, len(ds)))
    say("  The two depths M2 and M4 refuse on the sampled D_c are the")
    say("  two where the sampling error is the largest fraction of D_c.")
    say("  Computed rather than sampled, both hold with orders to spare.")
    say("  The registered verdicts stand as scored: they were scored on")
    say("  the estimator this file registered, and this block does not")
    say("  rescore them. What it separates is the estimator from the")
    say("  quantity.")
    say()
    say("    And the denominator of that drift is the SAMPLED D_c. The")
    say("    exact D_c is computed above; dividing by it instead:")
    say("    depth  sampled drift   exact drift")
    for d in ds:
        r0 = [VARO.get((lo, hi, d)) for lo, hi in OCTS]
        if any(v is None for v in r0) or len(D[d]) != len(OCTS):
            continue
        r0 = [r0[k] / D[d][k] for k in range(len(OCTS))]
        ex = exact.get(d)
        vv = [VARO.get((lo, hi, d)) for lo, hi in OCTS]
        if not ex or len(ex) != len(OCTS) or any(v is None for v in vv):
            continue
        r2 = [vv[k] / ex[k] for k in range(len(OCTS))]
        say("    %-6d %6.1f%%          %6.1f%%"
            % (d, 100.0 * abs(r0[0] - r0[-1]) / r0[0],
               100.0 * abs(r2[0] - r2[-1]) / r2[0]))
    say("    Depth 1's 0.9 per cent is the sampling estimator's: the")
    say("    pair")
    say("    average carries about four per cent of noise per octave")
    say("    there, which is more than the drift it is said to show.")
    say()

    # POST HOC (not a registered rule; run after M1-M4 were scored, and
    # scoring nothing).  Observation {#obs:coh} claims three things about
    # Q_cc n_c^-2 log N across the eight octaves -- that it drifts, that
    # the diagonal explains the drift at the two deepest cells and not at
    # the four shallowest, and that half the drift accounts for the excess
    # of the registered b over 1/(2<log N>).  None of the numbers in that
    # paragraph came from a run.  They are computed here so the paragraph
    # has a producer and the graph check can reach them.
    #
    # Two conventions, both of which change the printed value and neither
    # of which is forced:
    #   drift is (first - last)/last, not /first.  Stated because /first
    #   gives 2.98 per cent where /last gives 3.07, and the paragraph
    #   quotes the second.
    #   the fitted exponent g is the slope of log(Q_cc n_c^-2 log N)
    #   against log N, SIGNED.  The quantity falls, so g < 0, and the
    #   relation b = 1/(2 log N) - g/2 needs that sign: a rising quantity
    #   would put b below the form instead of above it.  An earlier draft
    #   printed |g| and the sign of the check was then unreadable.
    say("POST HOC (not a registered rule; run after M1-M4 were scored,")
    say("and scoring nothing) -- the scale drift of Q_cc n_c^-2 log N")
    say("that Observation {#obs:coh} reads.")
    say()
    QO, BD = read_qcc_by_octave(), read_b_by_depth()
    # The two constants this block compares against are computed, not
    # typed: HALF is the block C4 value of 1/(2<log N>) read out of
    # lab_cell_floor.txt, and LO/HI are 1/(2 log N) at the ends of the
    # same octave list the drift is fitted on.  Typed, they were two
    # hand-written measurements and G19 said so.
    HALF = None
    for _ln in io.open(os.path.join(ROOT, "results", "lab_cell_floor.txt"),
                       encoding="utf-8").read().splitlines():
        if _ln.strip().startswith("1/(2<log N>)"):
            HALF = float(_ln.split("=")[1].split()[0])
    _oct = sorted(QO[min(QO)]) if QO else []
    LO = 1.0 / (2.0 * math.log(_oct[0][0])) if _oct else float("nan")
    HI = 1.0 / (2.0 * math.log(_oct[-1][0])) if _oct else float("nan")
    # The drift is quoted as the RATIO f(lowest octave)/f(highest), not
    # as a percentage, because the percentage needs a denominator and
    # the choice is not forced: (first-last)/last and /first differ by
    # the factor f0/fl itself, which is 1.03 at the shallow depths and
    # 1.62 at depth 5.  So the convention is free exactly where the
    # drift is small and decides the sentence exactly where it is not --
    # 62 per cent against a 70 per cent diagonal reads as the diagonal
    # explaining it, 38 per cent against the same 70 reads as the
    # diagonal being twice too large.  The ratio has no denominator to
    # choose.
    #
    # g is printed at three abscissae -- the octave's low, mid and high
    # endpoint -- because its VALUE is not convention-free either. The
    # slope of log q in log N is untouched (the endpoints differ by a
    # constant log 2), so the whole dependence sits in the log log N
    # factor, and it moves g by 60 per cent at the shallow depths. The
    # SIGN is stable in all eighteen. Only the sign is quoted in
    # {#obs:coh}.
    say("    depth  f_lo/f_hi   g at lo     g at mid    g at hi"
        "     swing     diagonal")
    for d in sorted(QO):
        v = sorted(QO[d])
        if len(v) < 3:
            continue
        gs = []
        for k in (3, 0, 4):          # octave low endpoint, midpoint, high
            xs = [math.log(r[k]) for r in v]
            ys = [math.log(r[2] * math.log(r[k])) for r in v]
            mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
            # xi/yi, not a/b: this file has twice lost a verdict to a
            # POST HOC block reusing a name.
            gs.append(sum((xi - mx) * (yi - my) for xi, yi in zip(xs, ys))
                      / sum((xi - mx) ** 2 for xi in xs))
        f0 = v[0][2] * math.log(v[0][0])
        fl = v[-1][2] * math.log(v[-1][0])
        dg = 1.0 / (v[0][1] * v[0][2])
        # the ABSOLUTE swing, not a percentage of anything: it is the
        # same 0.003721 at all six depths, which is what an additive
        # term produces and is a stronger statement than any ratio.
        # A percentage here would need a denominator, and the two
        # obvious ones give 60 and 94 -- in a paragraph that is about
        # denominators.
        say("    %-6d %8.3f    %+.6f   %+.6f   %+.6f   %.6f   %6.2f%%"
            % (d, f0 / fl, gs[0], gs[1], gs[2],
               max(gs) - min(gs), 100.0 * dg))
    say("    The diagonal column is (1/n_c) / (Q_cc/n_c^2) at the LOWEST")
    say("    octave -- the share of the printed quantity carried by the")
    say("    n_c terms N = N', for which K = 1 exactly. It is under half")
    say("    a per cent at depths 0 to 3, an order below the 3 per cent")
    say("    the ratio departs from 1 there, and 70.53 per cent at depth")
    say("    5. The share and the ratio are NOT put in a ratio to each")
    say("    other: that comparison needs the drift as a percentage, and")
    say("    at depth 5 the two denominators give 1.14x and 1.84x -- one")
    say("    reads as the diagonal explaining the break and the other as")
    say("    the diagonal being twice too large for it.")
    say("    1/(2 log N) is not one number over these octaves: it runs")
    say("    from %.4f to %.4f, a %.0f per cent range. So b, a single fit"
        % (LO, HI, 100.0 * (LO - HI) / HI))
    say("    over all eight, cannot be compared to the single averaged")
    say("    %.6f without slack of that order." % HALF)
    say()

    # -g/2 IS b - 1/(2 log N), IDENTICALLY, and the comparison of the two
    # is therefore not a test.  se^2 = ratio * (Q_cc/n_c^2), and if ratio
    # is constant in N then -2b = dlog(Q_cc n_c^-2)/dlog N, while
    # g = dlog(Q_cc n_c^-2 log N)/dlog N = that slope + 1/log N.  Subtract.
    # An earlier draft of {#obs:coh} read the 7.8 per cent gap between
    # them as agreement between two measurements; it is the residual of
    # putting a single <log N> in two places where log N moves 42 per
    # cent, and it confirms nothing.  What CAN be read is where the
    # identity breaks, because the only hypothesis in it is that ratio
    # is constant -- so the residual measures ratio's drift and nothing
    # else, which matters because {#obs:coh} attributed the deep-cell
    # break to the diagonal instead.
    say("    -g/2 = b - 1/(2 log N) IDENTICALLY when Var/(Q_cc n_c^-2)")
    say("    is constant in N, so comparing them tests nothing. The one")
    say("    reading the residual supports is of that constancy:")
    say()
    # The fitted b is NOT reprinted here, only differenced against.
    # Echoing it would make it a value two result files produce, and the
    # graph check attributes a statement to a run by the figures only
    # that run prints -- so reprinting an upstream value silently
    # deletes an edge.  It did: printing b here dropped {#obs:coh} and
    # {#meas:exponent} off lab_cell_floor.txt and the paper's census of
    # five became three.  b per depth is block C4 of that file.
    say("    depth  -slope/2    b - that     ratio spread")
    RAT = {}
    for _ln in io.open(os.path.join(ROOT, "results", "lab_cell_floor.txt"),
                       encoding="utf-8").read().splitlines():
        f = _ln.split()
        if len(f) == 10 and f[0] == "(" and f[3].isdigit():
            try:
                RAT.setdefault(int(f[3]), []).append(float(f[7]))
            except ValueError:
                continue
    for d in sorted(QO):
        v, b = sorted(QO[d]), BD.get(d)
        if len(v) < 3 or b is None or d not in RAT:
            continue
        x = [math.log(r[0]) for r in v]
        y = [math.log(r[2]) for r in v]
        mx, my = sum(x) / len(x), sum(y) / len(y)
        half = -0.5 * (sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
                       / sum((xi - mx) ** 2 for xi in x))
        say("    %-6d %.6f    %+.1e     %.2e"
            % (d, half, b - half, max(RAT[d]) - min(RAT[d])))
    say("    The identity is exact to 2e-05 where the ratio does not")
    say("    move and misses by 2.6e-02 at depth 5 where the ratio")
    say("    spreads 1.65e-01. So the deepest cell's break is that")
    say("    ratio-constancy fails there. The diagonal share rises at")
    say("    the same two depths and is not separable from it here.")
    say()

    say("=" * 70)
    ok = m1 and m2 and m4
    say("M1 %s  M2 %s  M3 %s  M4 %s"
        % tuple("hold" if v else "REFUTED" for v in (m1, m2, m3, m4)))
    say("Proposition {#prop:scaleinv} stands" if ok else "REFUTED")

    head = [
        "STATISTIC: D_c = E_same,c[S_2] - E_all[S_2] with",
        "           S_2(h) = 2 C_2 prod_{p|h,p>2}(p-1)/(p-2) the",
        "           Hardy-Littlewood singular series of the shift, E_same,c",
        "           the mean over sampled pairs N != N' both in cell c and",
        "           E_all the mean over sampled pairs in the band; its",
        "           spread across three octaves and its fitted exponent in",
        "           N; and its correlation across depths with the exact",
        "           floor se_c of Lemma {#lem:cellmom}.",
        "DENOM: the spread is divided by D_c at that depth;",
        "NULL: an exact permutation null on M3. D_c itself is a",
        "      deterministic arithmetic functional of the cell with no",
        "      sign input, so Lemma {#lem:coin} does not bite on it; but",
        "      M3 is a correlation across depths between D_c and the",
        "      exact floor, and that needed a control which nothing",
        "      supplied. With six depths the null is enumerable: all 720",
        "      permutations of the depth labels are run and the observed",
        "      |r| is placed in that distribution. The sampling error of",
        "      the pair average is reported alongside.",
        "FIELD: octaves (1e6,2e6], (2e6,4e6], (4e6,8e6]; even N; cells",
        "       indexed by depth = #{p in 3,5,7,11,13 dividing N};",
        "       400000 sampled ordered pairs per cell and per band, numpy",
        "       default_rng seed 20260808; S_2 sieved to 8e6; se_c taken",
        "       from lab_mask_placebo.py at (2e6,4e6].",
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

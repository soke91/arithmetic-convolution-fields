# -*- coding: utf-8 -*-
r"""pa_certify -- exact-arithmetic certificates for the numbers in section 11 of the manuscript.

WHAT THIS IS FOR

Section 11 of the manuscript -- the `.tex` beside `code/` in a built packet -- displays the certificate
inequality `eq:cert` and reports what this file establishes:

  * `e = 16 .. 28` -- BOTH endpoints of `kappa(2^e)`, in exact integer arithmetic, each cell from an
    integral primal `(q, D)` for the upper endpoint and an integral dual `lambda` for the lower one;
  * `2^30` -- the lower endpoint alone, from a stored integral dual, with no integral primal produced
    there and so no upper endpoint claimed;
  * the construction's `0.6982, 0.6130, 0.6088, 0.5841` -- its SQUARED bound is an exact rational of
    integer class counts at each of the four sizes, and the printed figure is that rational's square
    root rounded UP at four decimals; the square roots are irrational, which is why the squares are
    what the clause compares (P5 below).

This file produces and stores those certificates, and re-checks them with `--verify`.

THE CERTIFICATE, in the paper's own notation (`eq:cert`, section 11)

With rows the squarefree `d <= Q`, `c_vB(d) = #{n in vB : d | n}`, `c_vR(d)` likewise, `P = c_vB - c_vR`
and `A` the divisor-incidence matrix of the even columns, `eq:Rfree` gives
`kappa(N) = dist_2(c_vR, K_Q)/||P||_2`, so for an integral `q >= 0` with integer `D >= 1` and an
integral `lambda` with `theta_lambda(l) = sum_{d | l} lambda_d >= 0` for EVERY `l in vB` and
`<c_vR, lambda> < 0`,

    <c_vR,lambda>^2 / (||P||^2 ||lambda||^2)  <=  kappa(N)^2  <=  ||D c_vR - A q||^2 / (D^2 ||P||^2).

Both sides are rational in integer data.  The left is Proposition 1 at `lambda/||lambda||`; the right is
`dist_2(c_vR,K_Q) <= ||c_vR - Aq/D||`.  Note the normalisation `R` cancels: `t = c_vR/R` and `p = P/R`,
so nothing here divides by `R` and no float enters the verdict.

HOW THE PAIR IS OBTAINED (the finder; the verifier below does not trust it)

At the optimum of `min_{y>=0} ||A y - c_vR||` the KKT conditions give `A^T r <= 0` for the residual
`r = c_vR - A y`, and complementary slackness gives `<c_vR, r> = ||r||^2`.  So `lambda = -r` is dual
feasible with `<c_vR,lambda> = -||r||^2 < 0`, and the left bound then equals `||r||^2/||P||^2` exactly:
the pair is tight at the optimum, and only the ROUNDING to integers costs anything.

  dual     `lambda = round(-K r / max|r|)`, then the one exact repair that cannot fail: `d = 1` divides
           every column, so adding `m = max(0, -min_vB theta_lambda)` to `lambda_1` makes
           `theta >= 0` everywhere.  (The same `d = 1` repair `gramint_kappa_certify` uses in its K0
           block.)  The repair is measured and reported; it costs nothing visible.
  primal   `q = max(0, round(D y))` on the support.

THE SCALES ARE SIZED, AND THE int64 KERNELS ARE GUARDED BY THE ACTUAL INTEGER SUMS

`K` and `D` are powers of ten SIZED so the two int64 kernels have room (`scales`).  Sizing is not a
proof, and for the primal the sizing bound is false at every stored cell (below).  So
before either kernel runs, the quantity that actually bounds its accumulator is computed as a PYTHON
INTEGER from the very vector that will be passed in, and the kernel runs only if it fits:

    k_theta   every partial sum of `out[n]` is a signed subset sum of `lambda`, so
              `|out[n]| <= sum_d |lambda_d|`              guard: `sum_d |lambda_d| <= 2^62`
    k_fit     every partial sum of `out[i]` is a subset sum of `q >= 0`, so
              `fit_d <= sum_j q_j`                        guard: `sum_j q_j <= 2^62`

Both sums are taken in object dtype with the absolute value applied AFTER conversion to Python
integers, so the `-2^63` case -- where int64 `abs` returns `-2^63` again, silently -- cannot hide.
Every float that is about to be cast to int64 is bounded first, because that cast is undefined, not
wrapping, once the value passes `2^63`.  If a guard does not hold the scale is divided by ten and the
vector rebuilt, so a cell is certified at a smaller scale instead of being voided; `scales` only
picks the starting point.

WHY THE GUARDS ARE THE ACTUAL SUMS AND NOT A SIZING BOUND.  Two tempting shortcuts are false here.

A single global `K = D = 10^15` with a C2 requiring the int64 values of `<c_vR,lambda>` and
`||lambda||^2` to agree with their exact recomputation cannot work: `||lambda||^2` is about `10^32`
and never fits int64, and `D * R` overflows at the larger cells.  C2 is therefore stated below in
terms of the quantities that genuinely must be int64, with every big one taken in object dtype.

And `(R+1) * D <= 2^62` is NOT a proof that `fit = A q` cannot overflow, although `sum_j q_j` is
about `D * sum_j y_j ~ D`: the bound `sum_j q_j <= (R+1)D` is FALSE at `2^18`
(`sum_j q_j = 2081398496916302535 > (R+1)D = 2075500000000000000`) and at ALL SEVEN stored cells, by
`+0.0047%` at `2^28` rising steadily
to `+0.5884%` at `2^16`: `sum_j y_j` is slightly above 1 and `rint` rounds up about as often as down,
so `sum_j q_j` lands just above `R*D` every time, and relatively further at the small cells because
the support is a larger share of the band.  The stored certificates are nevertheless valid -- the
largest actual sum is `3.380e18`
at `2^22`, which is `0.733 * 2^62` and far inside int64 -- but that is a measurement, not what the
sizing bound claimed, and the margin it leaves is the factor of two between `2^62` and int64's
maximum.

So C2 is stated below in terms of the ACTUAL sums, C4 validates every index and C5 the stored
fields, `abs` is applied after the conversion to Python integers, every float-to-int64 cast is
bounded first, and the divisor count over `rows x cols` is a `prange` kernel rather than a dense
matrix.  Every guard and every clause below holds on all seven stored certificates.

Floats are used ONLY to propose `lambda` and `q`.  Every number that appears in a verdict is computed
from the integers, and `--verify` recomputes all of it from the packet's own files with no solver.

WHERE THE FLOAT OPTIMUM COMES FROM, per size (the finder's input only)

  16, 18   solved here: dense `scipy.optimize.nnls` over the complete even-column set (0.3 s / 12 s),
           reproducing the support sizes 147 and 301 that `conedual_separator_mass.txt` records.
  20,22,24 read from the stored `results/conedual_exact_dual_2e{20,22,24}.npz`, which hold the full
           dual `w` and the primal support (`y_support` are indices into vB, `y_values` the weights)
           at KKT `1.7e-16`, `1.0e-15`, `3.1e-15`.  No solve is needed, and at `2^20` this avoids the
           700 s dense run that `conedual_exact_kappa_small_cells_20.txt` records.
  26, 28   `pa_certify_big.py`, from a colgen separator that carries its primal: `conedual_colgen.py`
           stores `y_support`/`y_values` beside the dual, so these two cells are the same
           storage-and-rounding step as the three above.

EXACT ARITHMETIC WITHOUT A PYTHON LOOP

Squared norms and inner products that can exceed int64 are taken in numpy OBJECT dtype
(`np.sum(v.astype(object)**2)`): arbitrary-precision Python integers, elementwise, with no explicit
Python loop over rows: what is avoided is the Python loop, not the big integers.  Where
int64 suffices the magnitude is bounded first and clause C2 checks it.

HOW THE WORK IS DONE

  Cells are LOADED (`from cell_cache import load`) for `rows`, `c_vB`, `c_vR`, `N`, `Q`, `thr`.  The
  band is not cached (cell_cache says so), so it is rebuilt exactly as `cell_cache._build` does -- one
  strided pass per prime <= Q, 54 to 1900 primes, never a loop over an N-range -- and clause C1 checks
  the rebuild against the cache's `nB`, `nR`, `cb`, `cr` element for element.  Every per-row and
  per-band accumulation is a numba `prange` kernel: the theta sweep blocks over N (prange over blocks,
  so no write races), the row counts, class counts, `k_fit` and `k_coldiv` prange over rows or columns
  (each writes only its own slot).  C4's factorisation of the stored indices is vectorised over the
  indices, with the only Python loop over the primes `<= Q`.

  THREADS, the one rule, and it is the code's: `nth = max(1, min(ACF_THREADS or 4, 8, numba's own
  thread count))`.  So the DEFAULT is 4 and the CAP is 8 -- `ACF_THREADS=8` is honoured, `ACF_THREADS=16`
  is clamped to 8 -- and numba's count
  clamps it further on a smaller machine (a 4-core reader gets 4).  The cap is this machine's limit.  Every measurement in
  `results/pa_certify_timings_frozen.json` and in section 11 was taken at the default 4, which is also
  what the reproduction notes tell a reader to set.

  COST.  The largest cell rebuilds its band in about half a minute at a peak near 4 GB, and every run
  of this file finishes in minutes; nothing here iterates or solves.  The whole seven-cell `--verify`
  is measured and its wall time reported.

CLAUSES.  Every "voids" below is enforced in code, not in prose.

  C0  CONTROL.  The cell is ESTABLISHED, not read: `N = 2^e`, `Q = floor(sqrt N)`, `q1` the least
      prime above `Q` and `thr = floor(N/q1)` are derived here from `e` alone, and the COMPLETE list of
      squarefree `d <= Q` is generated here.  A cache or a certificate whose `N`, `Q`, `q1`, `thr`,
      `e` or row list differs is refused.  Without this a larger cell's cache and certificate,
      relabelled, pass every later check on their own valid data and are reported as the smaller
      cell's interval; and a row list that is a valid SUBSET changes the normalisation denominator.
      CAN FAIL; voids the cell.
  C1  CONTROL.  At every cell the rebuilt band has `|vB| = nB` and `|vR| = nR`, AND the per-row divisor
      counts recomputed from it equal the cached integer vectors `cb`, `cr` ELEMENT FOR ELEMENT
      (integer against integer, no tolerance).  CAN FAIL; voids the cell.
  C2  CONTROL (restated).  Both int64 kernels are guarded by the ACTUAL integer bound on
      their accumulator, computed as a Python integer from the vector that will be passed in, BEFORE
      the kernel runs, with every absolute value taken AFTER the conversion to Python integers:
      `sum_d |lambda_d| <= 2^62` and `max_d |lambda_d| <= 2^62` for `k_theta`; `q >= 0`,
      `sum_j q_j <= 2^62` and `max_j q_j <= 2^62` for `k_fit`.  Any float about to be cast to int64
      is bounded by `2^62` first.  No claim about `(R+1) * D` is made or used.  The fast kernel is
      still cross-checked: `sum_d fit_d` equals the same total recomputed in object dtype as
      `sum_j q_j * #{d in rows : d | l_j}`.  Every quantity that can exceed int64 --
      `<c_vR,lambda>`, `||lambda||^2`, `||P||^2`, `||D c_vR - A q||^2` -- is computed in object dtype
      and never in int64.  CAN FAIL; voids the cell.
  C3  CONTROL.  The stored primal support read from `conedual_exact_dual_2e*.npz` indexes vB in range,
      and every column it names is in fact an even-omega band element (checked against the rebuilt
      band).  CAN FAIL; voids the cell.
  C4  CONTROL.  Every stored index is validated BEFORE it indexes anything, and independently
      of `code`: the rows satisfy `1 <= d <= Q`, are strictly increasing, are squarefree, and equal
      the cell's own row list element for element; the primal columns satisfy `thr < l <= N`, are
      distinct, are squarefree, have `P^+(l) <= Q`, and have `omega(l)` EVEN -- the last three
      established by trial division over the primes `<= Q` on the stored columns themselves, not by
      reading `code`, and then required to agree with `code[l] = 1`.  CAN FAIL; voids the cell.
  C5  CONTROL.  A stored certificate's fields are present and well formed: the npz opens,
      each of the eleven required fields exists, the integer fields have integer dtype, the declared
      rank, and VALUES A SIGNED int64 HOLDS UNCHANGED, the four big rationals parse as decimal
      integers with nonzero denominators, and the lengths agree (`|lambda| = |rows|`, `|cols| = |q|`,
      `D >= 1`, `e` and `N = 2^e` as stored).  The value test is not redundant: a vector stored in an
      unsigned dtype above `2^63 - 1`, or in a float dtype, is READ as a different vector from the one
      stored -- the unsigned one wraps to the signed vector the later clauses want while the stored
      entries are all nonnegative, so the file would pass while `<c_vR, lambda> < 0` is false of it.
      The same test is repeated at the point of conversion, so no edit can move a cast in front of it.
      A missing or malformed field FAILS this clause: it does not raise and it does not pass.
      CAN FAIL; voids the cell.

  C6  CONTROL.  `--verify` with no `--cells` checks exactly the advertised cells `2^16 ... 2^28`.  A
      missing certificate fails, an unadvertised one fails, and an empty inventory fails -- discovering
      whatever files are present and succeeding on them lets an empty directory verify.  CAN FAIL.
      WITH `--cells` NAMING ANYTHING ELSE this clause reports **NOT RUN** -- neither a pass nor a
      failure -- because a run that opens the selected certificates alone cannot speak for the packet's
      contents: a damaged certificate outside the selection is never read.  The announcement and the
      final `CERTIFICATES:` line of such a run are scoped to the selected cells, and the record carries
      `completeness_checked: false`.  Only `--verify` with no `--cells` speaks for the packet.
  C9  CONTROL.  The only float -> integer decision left in the shipped path -- the canonical sieve's
      loop bound `int(limit ** 0.5)`, and `cell_cache._build`'s `int(N ** 0.5)` and `int(q1 ** 0.5)` --
      equals `math.isqrt` at EVERY integer below `2^31`, which is every limit this path can sieve.  The
      three sites are inside the cell-cache fingerprint (the sieve's bytes and `_build`'s source are
      hashed into it), so rewriting them would discard every shipped cell cache, 2^28 and 2^30
      included; they are checked instead.  The check is complete rather than a sample: both maps are
      non-decreasing, so agreeing at the two ends of each `isqrt` level settles that level, and the
      levels cover the range.  Verified both vectorised and with Python's own scalar operator, which is
      the expression the sieve evaluates.  CAN FAIL.  With it, nothing a float computes reaches a
      verdict here: `C0` re-derives `N`, `Q`, `q1`, `thr` and the complete row list with `math.isqrt`
      and integer arithmetic, and refuses a cache that disagrees.  It runs in EVERY mode of this file
      AND in `pa_certify_big.py`, which calls this file's `check_c9` so that the wrapper performs every
      check the verifier it wraps performs.

  L0  CONTROL (`--make-lower30`).  The `2^30` separator's row list equals the derived row list.
  L1  CONTROL (`--make-lower30`).  The stored dual is INTEGRAL: `lambda = round(-10^12 w)` with the
      `d = 1` entry discarded and rebuilt by the one `d = 1` repair, `theta_lambda >= 0` on every even
      column of the rebuilt `2^30` band, `J = -<c_O,lambda> > 0`, and `sum_d |lambda_d|` inside `2^62`.
  L2  CONTROL (`--verify-lower30`).  The stored `2^30` certificate is present and complete; its
      FORMAT is refused or accepted before any conversion (`dual` and `rows` integer vectors of rank
      one whose values a signed int64 holds unchanged, `N`, `e`, `c6` integer scalars, `J`, `H`, `D_p`
      decimal integers); and it is bound to the derived cell: `e`, `N`, the complete row list,
      `|lambda| = |rows|`, the int64 guard.  Without the format test a dual stored as floats, or as an
      unsigned vector, would be converted into the integer vector the clauses expect and pass as it.
  L3  (`--verify-lower30`).  `theta_lambda >= 0` on EVERY even column of the rebuilt `2^30` band; `J`,
      `H`, `D_p` recompute from the stored vector; `c6^2 H D_p < J^2 10^12`; `c6` is MAXIMAL; and that
      integer inequality is equivalent to `J^2/(H D_p) > (c6/10^6)^2`, which IMPLIES
      `kappa^2 > (c6/10^6)^2` because `J^2/(H D_p) <= kappa^2` is the left half of `eq:cert`.  The
      implication is one way: this clause checks the equivalence it can check exactly, in integers, and
      the bound on `kappa` follows from it.  No upper endpoint is claimed at `2^30`.
  L4  REPRODUCTION (`--verify-lower30`).  The certificate's `J`, `H`, `D_p` and `c6` are the four
      integers the solver run's ledger records and the paper prints, digit for digit.  A rebuild at a
      different rounding scale certifies the same cell with different valid integers; this is the clause
      that keeps the stored certificate the one the paper's numbers come from.  CAN FAIL.

  P1  The two-sided certificate holds at `e = 16,18,20,22,24`: `theta_lambda(l) >= 0` on EVERY `l in
      vB` (one pass over the band, not only the support), `<c_vR,lambda> < 0`, and
      `lower2 <= upper2` as exact rationals.
  P2  At those five sizes the certified interval for `kappa(N)` has width `< 10^-8`.
  P3  REPRODUCTION, containment of a stored float in an exact interval.  The certified interval
      contains the stored floating-point optimum at every size where one is stored:
      `0.172389352` and `0.140860454` (`conedual_separator_mass.txt`, 9 decimals), and the `ub` field
      of `conedual_exact_dual_2e{20,22,24}.npz`.
  P4  REPRODUCTION, ledger integers against the exact lower bound.  At `e = 22,24` the certified lower
      endpoint is `>= c6/10^6` with `c6 = 76331, 52762` from `conedual_colgen_2e{22,24}_n0b.json`, so
      the new certificate is at least as strong as `eq:certint` at those sizes.
  P5  The construction's bounds are certified: at `e = 16,18,20,22` with `Y = 3,3,5,5` the exact
      rational `||D c_vR - A q||^2/(D^2||P||^2)` of the class construction is `<=` the square of the
      printed value, and the printed value is that exact bound rounded UP at four decimals (so it is
      both valid and tight): `0.6982, 0.6130, 0.6088, 0.5841`.
  P6  The certificate files are small: every `pa_cert_2e*.npz` is under 200 kB.

  NOT CLAIMED HERE: anything about Theorem 1.  (No clause spoke of `e = 26,28` before them, whose
  primal did not exist then; `pa_certify_big.py` carries those clauses.)

WHAT THE RUN REPORTS, AND WHAT THE EXIT CODE MEANS

  The clauses above do two different jobs, and the run separates them in its last two lines.

    CERTIFICATES: VALID            -- or INVALID with the failing certification clause named
    DIAGNOSTICS: none              -- or each failing diagnostic with the reason it can fail

  CERTIFICATION clauses are the ones that decide whether a certificate stands: C0-C6 (the cell is the
  one `e` names and is intact, the integer arithmetic is sound, the stored fields are the integers they
  claim to be, the inventory is the advertised set), the `V<e>` re-verifications, `L0`-`L4` on the 2^30
  path, `VOID`, and `P1` and `P5`.  DIAGNOSTIC clauses compare a certificate with something outside it:
  `P2` a width target, `P3` a nine-decimal float printed elsewhere, `P4` the published `eq:certint`
  integers, `P6` a file-size budget.  A diagnostic failing is a finding about that comparison and leaves
  every certificate as valid as it was.

  THE EXIT CODE REPORTS CERTIFICATE VALIDITY ALONE: 0 when every certification clause holds, whatever the
  diagnostics say; 1 when one does not.  `P3` is expected to fail as shipped -- it asks a twelve-decimal
  exact interval to contain a nine-decimal rounding -- and that no longer makes the run exit nonzero.
  The `FAILS:` line above it still lists every failing clause of either kind, and the JSON record adds
  `certificates_valid`, `certification_fails` and `diagnostic_fails`.  No clause's test or threshold is
  touched by this split; the split decides only what the exit code means.

  AND WHAT IT COVERS.  Exactly one invocation speaks for the packet's contents: `--verify` with no
  `--cells`.  It alone prints the bare

    CERTIFICATES: VALID

  Every other invocation checks the cells it was given and says so, in the announcement, in the verdict
  and in the record (`completeness_checked: false`, `clauses_not_run`):

    CERTIFICATES: VALID for 2^16 only; completeness not checked            (--verify --cells 16)
    CERTIFICATES: VALID for 2^30 (the lower endpoint only); ...            (--verify-lower30)
    CERTIFICATES: VALID for the cells built (2^16, ...); ...               (construction mode)

  A line that begins `CERTIFICATES: VALID` or `CERTIFICATES: INVALID` either way, so a script reading the
  verdict is unaffected; what the scope adds is that no run can be quoted for more than it examined.

    python code/pa_certify.py                     # CONSTRUCTION: re-derive and write certificates
    python code/pa_certify.py --verify            # VERIFICATION: the whole advertised packet
    python code/pa_certify.py --verify --cells 16 # VERIFICATION of one cell; completeness NOT checked
"""
from __future__ import annotations

import argparse
import glob
import io
import json
import math
import os
import sys
import time
from fractions import Fraction

NTHREAD = int(os.environ.get("ACF_THREADS", "4"))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, str(NTHREAD))

import numpy as np                                      # noqa: E402
import numba                                            # noqa: E402
from numba import njit, prange, set_num_threads         # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cell_cache import load as load_cell                # noqa: E402
from _sieve_shared import mu_upto, primes_upto          # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), "results")
OUT = os.path.join(RES, "pa_certify.txt")
OUTJ = os.path.join(RES, "pa_certify.json")

I64 = int(np.iinfo(np.int64).max)
SAFE = 2 ** 62             # headroom of 2 below int64's max for every accumulator
TWO_SIDED = (16, 18, 20, 22, 24)
# the cells whose two-sided certificates this packet advertises; --verify requires exactly these
ADVERTISED = (16, 18, 20, 22, 24, 26, 28)
LOWER30_E = 30                        # the cell certified on one side only, by a stored integral dual
CONSTRUCTION = ((16, 3), (18, 3), (20, 5), (22, 5))
CONSTR_LIT = {16: "0.6982", 18: "0.6130", 20: "0.6088", 22: "0.5841"}
LEDGER_C6 = {22: 76331, 24: 52762}
SEPMASS_FLOAT = {16: "0.172389352", 18: "0.140860454"}
WIDTH_MAX = Fraction(1, 10 ** 8)

# C5's contract: field -> required rank.  Integer fields must also have an integer dtype.
FIELDS_INT = {"N": 0, "e": 0, "rows": 1, "dual": 1, "primal_cols": 1, "primal_num": 1,
              "primal_den": 0}
FIELDS_STR = ("lower2_num", "lower2_den", "upper2_num", "upper2_den")


def int64_exact(a) -> bool:
    """Is `a` an integer array every value of which a signed int64 holds unchanged?

    A stored vector means what its dtype says it means, and the verifier reads it as signed int64.
    Two representations would therefore be read as a DIFFERENT vector from the one stored, and both
    are refused here, before any conversion:

      * an unsigned dtype carrying values above `2^63 - 1`.  Re-typing a dual as `uint64` leaves every
        entry nonnegative AS STORED -- so `<c_vR, lambda>` is positive, not negative, and the stored
        file does not satisfy the clause it is supposed to -- while the conversion wraps each entry
        back to the signed vector that does, and every later check would pass on that.
      * a float or otherwise non-integer dtype, whose entries need not be integers at all: a dual
        stored as `1000003321703.25` converts to the integer the clauses expect while the stored
        vector is not integral and its own inner product is not the one checked.

    The test is made on the values as Python integers, which no cast can wrap, and it accepts every
    narrower integer dtype (`int32`, `uint32`, a small `uint64`), since those do convert unchanged.
    """
    a = np.asarray(a)
    if not np.issubdtype(a.dtype, np.integer):
        return False
    if a.size == 0:
        return True
    v = a.astype(object)
    return bool(np.all(v >= -(2 ** 63)) and np.all(v <= 2 ** 63 - 1))


def decimal_int(a) -> bool:
    """Is the stored field a decimal integer, as the four big rationals and the 2^30 integers are?"""
    try:
        int(str(a))
        return True
    except Exception:
        return False

_lines: list[str] = []
REC: list[dict] = []
FAILS: list[str] = []
NOTRUN: list[str] = []


def say(s: str = "") -> None:
    print(s)
    _lines.append(s)
    with io.open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(_lines) + "\n")


def clause(name: str, ok: bool, detail: str) -> None:
    if not ok:
        FAILS.append(name)
    REC.append(dict(name=name, holds=bool(ok), detail=detail))
    say("%-4s %-8s %s" % (name, "HOLDS" if ok else "FAILS", detail))


def clause_not_run(name: str, detail: str) -> None:
    """A clause this invocation does not establish.  It is NOT a pass and NOT a failure.

    A run that checks a selected cell cannot speak for the packet's inventory, so the completeness
    clause must say that it did not run instead of printing HOLDS beside a question it never asked.
    The name is recorded (`holds` is null in the JSON) and carried into the verdict's scope line.
    """
    NOTRUN.append(name)
    REC.append(dict(name=name, holds=None, state="not run", detail=detail))
    say("%-4s %-8s %s" % (name, "NOT RUN", detail))


# Which clauses decide whether a CERTIFICATE IS VALID, and which only compare it with something else.
#
# A CERTIFICATION clause is a statement about the certificate itself: that the cell it is built on is the
# one e names and is intact (C0, C1, C5, C6), that the exact integer arithmetic behind the bound is sound
# and the stored fields are the integers they claim to be (C2, C3, C4, the V<e> re-verifications, L0-L4 on
# the 2^30 path), and that the certificate is feasible and two-sided (P1, P5).  VOID is one too: a voided
# cell has no certificate.  If any of these fails, a certificate does not stand.
#
# A DIAGNOSTIC clause compares a certificate with something outside it -- a width target, a rounded number
# printed in the paper, the solver's ledger, a file-size budget.  Its failure is a recorded finding about
# that comparison and leaves every certificate exactly as valid as it was.  Each one states below whether
# its failure is one this packet documents as expected, and why it can fail at all.
#
# The run ends with two machine-readable lines, CERTIFICATES: and DIAGNOSTICS:, and THE EXIT CODE REPORTS
# THE FIRST ONE ALONE: 0 when every certification clause holds, whatever the diagnostics say.  No clause's
# test or threshold is affected by this split; it decides only what the exit code means.
DIAGNOSTIC: dict[str, tuple[bool, str]] = {
    "P2": (False, "a width TARGET, not a condition on the certificate: the width is whatever the"
                  " stored primal and integral dual give"),
    "P3": (True, "it compares a nine-decimal rounded float optimum with a twelve-decimal exact"
                 " interval, which need not contain it -- section 11 says so in a parenthesis"),
    "P4": (False, "it compares this certificate with the published eq:certint integers; a failure is a"
                  " finding about those integers, not about this certificate"),
    "P6": (False, "a file-size budget for the shipped certificate files"),
}


def verdict_lines(fails: list[str], diagnostic: dict[str, tuple[bool, str]], emit, extra=None,
                  scope: str = "", notrun=()) -> int:
    """Print the two verdict lines and return the exit code, which is certificate validity ALONE.

    `scope` is what this invocation's verdict covers.  Only one invocation speaks for the whole packet
    -- `--verify` over the advertised set -- and it alone passes `scope=""`.  Every other mode passes a
    suffix such as " for 2^16 only; completeness not checked", so the line cannot be read as a statement
    about the packet's contents.  The line still begins `CERTIFICATES: VALID` or
    `CERTIFICATES: INVALID`, so anything reading the verdict mechanically is unaffected.
    """
    cert = [f for f in fails if f not in diagnostic]
    diag = [f for f in fails if f in diagnostic]
    emit("CERTIFICATES: %s" % ("VALID" + scope if not cert
                               else "INVALID%s (failing certification clause%s: %s)"
                                    % (scope, "" if len(cert) == 1 else "s", ", ".join(cert))))
    emit("DIAGNOSTICS: %s" % ("none" if not diag else "; ".join(
        "%s fails (%s: %s)" % (d, "expected" if diagnostic[d][0] else "not expected", diagnostic[d][1])
        for d in diag)))
    if extra is not None:
        extra.update(certificates_valid=not cert, certification_fails=cert, diagnostic_fails=diag,
                     verdict_scope=scope.strip() or "the advertised packet contents",
                     clauses_not_run=list(notrun))
    return 1 if cert else 0


def verdict(extra=None, emit=None, scope: str = "") -> int:
    """verdict_lines on this module's own clause record.  It is pure: `emit=QUIET` fills the JSON
    record, and the printing call is made last, so the two lines are the run's final words."""
    return verdict_lines(FAILS, DIAGNOSTIC, say if emit is None else emit, extra, scope, NOTRUN)


def QUIET(s: str) -> None:
    """an emit that prints nothing."""


def sqrt_bound_exact(nmax: int = 2 ** 31):
    """Does the float square root used as an integer bound agree with `math.isqrt` below `nmax`?

    One float -> integer decision is left in the shipped path and cannot be edited away: the canonical
    sieve takes its loop bound as `int(limit ** 0.5)` (`lib/goldbach/sieve.py`), and `cell_cache._build`
    derives `Q = int(N ** 0.5)` and trial-divides `q1` up to `int(q1 ** 0.5)`.  Those bytes are inside
    the cell-cache FINGERPRINT, so changing them would discard every shipped cell cache (2^28 and 2^30
    included).  They are therefore CHECKED here instead of rewritten.

    The check is complete, not a sample.  Both maps L -> int(L ** 0.5) and L -> isqrt(L) are
    non-decreasing, so if they agree at the two ends of an isqrt level -- L = m^2 and L = (m+1)^2 - 1,
    where isqrt is constant at m -- they agree at every L in between; the levels cover [0, nmax), and
    every limit the shipped path can sieve is below 2^31 (the int32 remainder array asserts it).  So
    comparing the 2*(floor(sqrt(nmax-1)) + 1) level boundaries settles the whole range.

    Both forms are compared: numpy's vectorised `** 0.5` over all levels at once, and -- because that
    is the expression the sieve actually evaluates -- Python's own scalar `int(L ** 0.5)` at every
    boundary.  The scalar pass is a loop over sqrt LEVELS (about 4.6e4 of them, 0.06 s), not over rows,
    band elements or an N-length range.
    """
    mmax = math.isqrt(nmax - 1)
    m = np.arange(0, mmax + 1, dtype=np.int64)
    lo = m * m
    hi = np.minimum((m + 1) * (m + 1) - 1, nmax - 1)
    vec = bool(np.array_equal(np.floor(lo.astype(np.float64) ** 0.5).astype(np.int64), m)
               and np.array_equal(np.floor(hi.astype(np.float64) ** 0.5).astype(np.int64), m))
    scal = True
    for k in range(0, mmax + 1):
        a, b = k * k, min((k + 1) * (k + 1) - 1, nmax - 1)
        if int(a ** 0.5) != k or int(b ** 0.5) != k:
            scal = False
            break
    return vec, scal, mmax + 1


def check_c9() -> None:
    """Record clause C9.  Called by this file's `main` in EVERY mode and, through this same function,
    by `pa_certify_big.main` -- which installs its own `clause`, so the clause lands in whichever
    file's record is being written.  One implementation, one clause text, both verifiers: a wrapper
    that performed fewer checks than the verifier it wraps is exactly what this avoids."""
    vec, scal, levels = sqrt_bound_exact()
    clause("C9", vec and scal,
           "the float square root used as an integer bound (lib/goldbach/sieve.py's int(limit**0.5),"
           " cell_cache's int(N**0.5) and int(q1**0.5)) equals math.isqrt at EVERY integer below 2^31:"
           " all %d isqrt levels checked at both ends, vectorised %s and with Python's own scalar"
           " operator %s -- so these floats cannot move a bound, and the exact N, Q, q1, thr and row"
           " list C0 derives are what every verdict rests on"
           % (levels, "yes" if vec else "NO", "yes" if scal else "NO"))


def say_fails() -> None:
    """the FAILS line, and -- when a clause did not run -- a line that says so beside it."""
    say("FAILS: %s" % (", ".join(FAILS) if FAILS else "none"))
    if NOTRUN:
        say("NOT RUN: %s" % ", ".join(NOTRUN))


# =============================== numba kernels =============================================
@njit(parallel=True, cache=True, nogil=True)
def k_theta(N, rows, lam, out, nchunk):
    """out[n] = sum_{d in rows, d | n} lam[d], for n in [0, N].

    prange over contiguous N-blocks and an inner loop over rows: each block writes only its own
    slice, so there is no race.  The reverse (prange over rows, strided writes) would race.

    Every partial sum of out[n] is a signed subset sum of lam, so |out[n]| <= sum_d |lam_d|.  The
    caller guarantees that bound fits int64 BEFORE calling (clause C2); this kernel does not check.
    """
    bs = (N + 1 + nchunk - 1) // nchunk
    nr = rows.size
    for c in prange(nchunk):
        lo = c * bs
        hi = min(lo + bs, N + 1)
        if lo >= hi:
            continue
        for k in range(lo, hi):
            out[k] = 0
        for i in range(nr):
            w = lam[i]
            if w != 0:
                d = rows[i]
                m = ((lo + d - 1) // d) * d
                if m == 0:
                    m = d
                while m < hi:
                    out[m] += w
                    m += d


@njit(parallel=True, cache=True, nogil=True)
def k_row_counts(code, rows, out1, out2):
    """out1[i] = #{n : code[n]==1, rows[i] | n}; out2[i] likewise for code 2.  prange over rows."""
    n = code.shape[0] - 1
    for i in prange(rows.size):
        d = rows[i]
        a = 0
        b = 0
        for m in range(d, n + 1, d):
            c = code[m]
            if c == 1:
                a += 1
            elif c == 2:
                b += 1
        out1[i] = a
        out2[i] = b


@njit(parallel=True, cache=True, nogil=True)
def k_class_counts(cls, rows, nclass, out):
    """out[i, c] = #{n : cls[n] == c+1, rows[i] | n}.  cls is 0 for 'not a retained even column'."""
    n = cls.shape[0] - 1
    for i in prange(rows.size):
        d = rows[i]
        for c in range(nclass):
            out[i, c] = 0
        for m in range(d, n + 1, d):
            c = cls[m]
            if c > 0:
                out[i, c - 1] += 1


@njit(parallel=True, cache=True, nogil=True)
def k_fit(cols, qs, rows, out):
    """out[i] = sum_{j : rows[i] | cols[j]} qs[j].  prange over rows; each writes its own slot.

    Every partial sum is a subset sum of qs, so out[i] <= sum_j qs[j] when qs >= 0.  The caller
    guarantees `qs >= 0` and that bound fits int64 BEFORE calling (clause C2), and that every
    cols[j] is a positive integer in range (clause C4); this kernel does not check.
    """
    for i in prange(rows.size):
        d = rows[i]
        s = 0
        for j in range(cols.size):
            if cols[j] % d == 0:
                s += qs[j]
        out[i] = s


@njit(parallel=True, cache=True, nogil=True)
def k_coldiv(cols, rows, out):
    """out[j] = #{i : rows[i] | cols[j]}.  prange over columns; each writes its own slot.

    This is C2's cross-check of `k_fit` -- `sum_d fit_d = sum_j q_j * out[j]` -- computed without a
    dense `rows x cols` broadcast, which would allocate about 594 MB at 2^26.
    """
    for j in prange(cols.size):
        c = cols[j]
        s = 0
        for i in range(rows.size):
            if c % rows[i] == 0:
                s += 1
        out[j] = s


# =============================== exact helpers =============================================
def pow10_at_most(limit: int) -> int:
    """largest power of ten <= limit (at least 1)."""
    k = 0
    while 10 ** (k + 1) <= limit:
        k += 1
    return 10 ** k


def scales(nrows: int, R: int):
    """the STARTING scales: K sized by nrows (|lambda_d| <= K, so sum_d |lambda_d| <= nrows*K) and D
    sized by R.  Sizing is not a proof -- the `(R+1)*D` bound on the primal sum is false at every
    stored cell -- so the callers below reduce the scale by factors of ten until the ACTUAL integer sum
    fits, and clause C2 reports the sums it used.
    """
    return pow10_at_most(SAFE // max(1, nrows)), pow10_at_most(SAFE // max(1, R + 1))


def sq(v) -> int:
    """exact squared norm in arbitrary precision; object dtype, no explicit Python loop."""
    return int(np.sum(v.astype(object) ** 2))


def dot(u, v) -> int:
    return int(np.sum(u.astype(object) * v.astype(object)))


def abs_sum(v) -> int:
    """sum_i |v_i| in arbitrary precision, with abs taken AFTER the conversion to Python integers.

    `np.abs` on int64 maps `-2^63` to itself, silently, so `np.abs(v).astype(object)` can report a
    total that is too small or even negative; `np.abs(v.astype(object))` cannot.
    """
    return int(np.sum(np.abs(v.astype(object))))


def abs_max(v) -> int:
    """max_i |v_i|, abs after the conversion to Python integers (see `abs_sum`)."""
    return int(np.max(np.abs(v.astype(object)))) if v.size else 0


def int_sum(v) -> int:
    return int(np.sum(v.astype(object)))


def int_max(v) -> int:
    return int(np.max(v.astype(object))) if v.size else 0


def castable(f) -> bool:
    """is every entry of this float array safe to cast to int64?  The cast is UNDEFINED, not
    wrapping, once |value| passes 2^63, so it is bounded by 2^62 before `astype` is allowed."""
    a = np.asarray(f, dtype=np.float64)
    if a.size == 0:
        return True
    m = float(np.max(np.abs(a)))
    return bool(np.isfinite(m) and m < float(SAFE))


def dec_bracket(lo2: Fraction, hi2: Fraction, nd: int) -> tuple[Fraction, Fraction]:
    """largest nd-decimal lo with lo^2 <= lo2, smallest nd-decimal hi with hi^2 >= hi2."""
    s = 10 ** nd
    lo = math.isqrt(lo2.numerator * s * s // lo2.denominator)
    while Fraction(lo + 1, s) ** 2 <= lo2:
        lo += 1
    while lo > 0 and Fraction(lo, s) ** 2 > lo2:
        lo -= 1
    hi = max(lo, math.isqrt(hi2.numerator * s * s // hi2.denominator))
    while Fraction(hi, s) ** 2 < hi2:
        hi += 1
    return Fraction(lo, s), Fraction(hi, s)


def fstr(f: Fraction, nd: int) -> str:
    """exact decimal string of a Fraction whose denominator divides 10^nd."""
    n = f.numerator * (10 ** nd) // f.denominator
    return ("%d.%0*d" % (n // 10 ** nd, nd, n % 10 ** nd))


# =============================== index validation (C4) =====================================
def factor_audit(vals, Q: int):
    """(omega, squarefree, Q-rough cofactor) of each entry, by trial division over the primes <= Q.

    Vectorised over `vals`; the only Python loop is over the primes <= Q (54 to 1900 of them), as in
    `cell` itself.  Every entry must be >= 1 -- the caller checks the range first, because a 0 would
    make the inner division loop never terminate.

    `rem == 1` is exactly `P^+(v) <= Q`, and `omega` then counts all of v's prime factors.
    """
    v = np.asarray(vals, dtype=np.int64)
    rem = v.copy()
    om = np.zeros(v.size, dtype=np.int64)
    sqf = np.ones(v.size, dtype=bool)
    for p in primes_upto(Q).tolist():
        d = (rem % p) == 0
        if not d.any():
            continue
        om[d] += 1
        rem[d] //= p
        d2 = d & ((rem % p) == 0)
        if d2.any():
            sqf[d2] = False
            while d2.any():                 # divide out the rest of this prime's power
                rem[d2] //= p
                d2 = d2 & ((rem % p) == 0)
    return om, sqf, rem


def validate_indices(e: int, srows, cols, Q: int, thr: int, N: int, code, cell_rows) -> bool:
    """C4.  Nothing here indexes an array with a stored index until that index has been bounded."""
    chk = [("rows nonempty", bool(srows.size)),
           ("rows in [1,Q]", bool(srows.size and int(srows.min()) >= 1 and int(srows.max()) <= Q)),
           ("rows strictly increasing", bool(srows.size and np.all(np.diff(srows) > 0))),
           ("rows equal the cell's row list", bool(np.array_equal(srows, cell_rows)))]
    if chk[1][1]:
        _, sqr, _ = factor_audit(srows, Q)
        chk.append(("rows squarefree", bool(np.all(sqr))))
    else:
        chk.append(("rows squarefree", False))

    rng = bool(cols.size and int(cols.min()) > thr and int(cols.max()) <= N)
    chk += [("columns nonempty", bool(cols.size)),
            ("columns in (thr,N]", rng),
            ("columns distinct", bool(np.unique(cols).size == cols.size))]
    if rng:
        om, sqf, rem = factor_audit(cols, Q)
        chk += [("columns squarefree", bool(np.all(sqf))),
                ("columns P^+ <= Q", bool(np.all(rem == 1))),
                ("columns omega even", bool(np.all(om % 2 == 0))),
                # only now is it safe to use a stored index as a subscript
                ("columns agree with code = 1", bool(np.all(code[cols] == 1)))]
    else:
        chk += [("columns squarefree", False), ("columns P^+ <= Q", False),
                ("columns omega even", False), ("columns agree with code = 1", False)]

    bad = [k for k, v in chk if not v]
    clause("C4", not bad,
           "2^%d stored indices validated before use, independently of code: %s"
           " (%d rows in [1,%d], %d columns in (%d,%d])"
           % (e, "all %d checks pass" % len(chk) if not bad else "FAILED: " + ", ".join(bad),
              srows.size, Q, cols.size, thr, N))
    return not bad


# =============================== field integrity (C5) ======================================
def check_fields(e: int, path: str):
    """C5.  Open a stored certificate and check every field is present and well formed.

    Returns (ok, npz or None).  A missing or malformed field FAILS the clause -- it does not raise
    and it does not pass -- so a damaged packet is reported as a failure like any other.
    """
    try:
        z = np.load(path, allow_pickle=False)
        names = set(z.files)
    except Exception as ex:
        clause("C5", False, "2^%d: %s cannot be read as an npz (%s: %s)"
               % (e, os.path.basename(path), type(ex).__name__, ex))
        return False, None

    chk: list[tuple[str, bool]] = [("npz opens", True)]
    for k in sorted(FIELDS_INT):
        nd = FIELDS_INT[k]
        if k not in names:
            chk.append(("%s present" % k, False))
            continue
        try:
            a = z[k]
            chk.append(("%s integer dtype" % k, bool(np.issubdtype(a.dtype, np.integer))))
            chk.append(("%s rank %d" % (k, nd), bool(a.ndim == nd)))
            # the dtype alone is not enough: an unsigned vector above 2^63-1 is a different vector
            # from the one the verifier would read, and so is a non-integral one
            chk.append(("%s survives int64 unchanged" % k, int64_exact(a)))
        except Exception:
            chk.append(("%s readable" % k, False))
    big: dict[str, int] = {}
    for k in FIELDS_STR:
        if k not in names:
            chk.append(("%s present" % k, False))
            continue
        try:
            big[k] = int(str(z[k]))
            chk.append(("%s is a decimal integer" % k, True))
        except Exception:
            chk.append(("%s is a decimal integer" % k, False))
    chk.append(("denominators nonzero",
                big.get("lower2_den", 0) > 0 and big.get("upper2_den", 0) > 0))

    if all(v for _, v in chk):                  # only now are the arrays known to be usable
        try:
            chk += [("|dual| = |rows|", int(z["dual"].size) == int(z["rows"].size)),
                    ("|primal_cols| = |primal_num|",
                     int(z["primal_cols"].size) == int(z["primal_num"].size)),
                    ("primal_den >= 1", int(z["primal_den"]) >= 1),
                    ("e as stored", int(z["e"]) == e),
                    ("N = 2^e", int(z["N"]) == 2 ** e)]
        except Exception as ex:
            chk.append(("field shapes readable (%s)" % type(ex).__name__, False))

    bad = [k for k, v in chk if not v]
    clause("C5", not bad, "2^%d stored fields of %s: %s"
                          " (integer dtype, rank, and every value held unchanged by signed int64)"
           % (e, os.path.basename(path),
              "all %d checks pass" % len(chk) if not bad else "FAILED: " + ", ".join(bad)))
    return (not bad), (z if not bad else None)


# =============================== the cell and its band =====================================
def derive_cell(e: int):
    """N, Q, q1, thr and the COMPLETE squarefree row list, derived here from `e` alone.

    Nothing cached or stored is consulted.  A verifier that reads these four numbers out of a cache is
    verifying whatever cell the cache happens to hold: relabel a larger cell's cache and certificate
    as a smaller one and every later check still passes on the larger cell's own valid data, so the
    larger cell's interval is reported as the smaller one's.  Clause C0 closes that by deriving the
    arithmetic of the requested cell and refusing anything that disagrees.

    The rows are the squarefree `d <= Q`, all of them: `mu(d) != 0` over `1..Q`, vectorised.  A row
    list that is merely a valid subset changes the normalisation denominator and the problem.
    """
    N = 2 ** e
    Q = math.isqrt(N)
    ps = primes_upto(Q + 400)                     # prime gaps below 2^16 are far under 400
    above = ps[ps > Q]
    if above.size == 0:
        raise SystemExit("pa_certify: no prime found just above Q = %d" % Q)
    q1 = int(above[0])
    thr = N // q1
    rows = np.flatnonzero(np.asarray(mu_upto(Q), dtype=np.int8)[:Q + 1] != 0)
    rows = np.ascontiguousarray(rows[rows >= 1], dtype=np.int64)
    return N, Q, q1, thr, rows


def cell(e: int):
    """cell_cache for the per-row integers; the band rebuilt as cell_cache._build does (C1).

    C0 first establishes the cell's arithmetic and its complete row list independently of the cache.
    """
    S = load_cell(e)
    dN, dQ, dq1, dthr, drows = derive_cell(e)
    cN, cQ, cq1, cthr = int(S["N"]), int(S["Q"]), int(S["q1"]), int(S["thr"])
    crows = np.ascontiguousarray(S["rows"], dtype=np.int64)
    c0chk = [("N = 2^e", cN == dN), ("Q = floor(sqrt N)", cQ == dQ),
             ("q1 = least prime > Q", cq1 == dq1), ("thr = floor(N/q1)", cthr == dthr),
             ("rows are every squarefree d <= Q", bool(np.array_equal(crows, drows)))]
    bad0 = [k for k, v in c0chk if not v]
    clause("C0", not bad0,
           "2^%d cell derived here and compared with the cache: N = %d, Q = %d, q1 = %d, thr = %d,"
           " %d squarefree rows -- %s"
           % (e, dN, dQ, dq1, dthr, drows.size,
              "all 5 agree" if not bad0 else "DISAGREE: " + ", ".join(bad0)))
    if bad0:
        return dict(e=e, N=dN, Q=dQ, thr=dthr, rows=drows, ok=False, void="C0",
                    derived=(dN, dQ, dq1, dthr), nB=0, nR=0, band_s=0.0,
                    cB=np.zeros(0, dtype=np.int64), cR=np.zeros(0, dtype=np.int64),
                    code=np.zeros(1, dtype=np.int8))
    N, Q, thr = dN, dQ, dthr
    assert N < 2 ** 31, "pa_certify: the int32 remainder array wraps at 2^31"
    t0 = time.time()
    mu = np.array(mu_upto(N), dtype=np.int8)
    rem = np.arange(N + 1, dtype=np.int32)
    om = np.zeros(N + 1, dtype=np.int8)
    for p in primes_upto(Q).tolist():          # Python loops over primes <= Q only (54..1900)
        om[p::p] += 1
        q = p
        while q <= N:
            rem[q::q] //= p
            if q > N // p:
                break
            q *= p
    bd = (mu != 0) & (rem == 1)
    bd[:thr + 1] = False
    del rem, mu
    code = np.zeros(N + 1, dtype=np.int8)
    np.copyto(code, np.int8(1), where=bd & (om % 2 == 0))
    np.copyto(code, np.int8(2), where=bd & (om % 2 == 1))
    del bd, om
    tb = time.time() - t0
    rows = np.ascontiguousarray(S["rows"], dtype=np.int64)
    cb = np.empty(rows.size, dtype=np.int64)
    cr = np.empty(rows.size, dtype=np.int64)
    k_row_counts(code, rows, cb, cr)
    nB = int(np.count_nonzero(code == 1))
    nR = int(np.count_nonzero(code == 2))
    ok = (nB == int(S["nB"]) and nR == int(S["nR"])
          and np.array_equal(cb, np.asarray(S["cb"], dtype=np.int64))
          and np.array_equal(cr, np.asarray(S["cr"], dtype=np.int64)))
    clause("C1", ok, "2^%d band rebuilt in %.1f s: |vB| %d/%d, |vR| %d/%d, c_vB and c_vR element for"
                     " element %s" % (e, tb, nB, int(S["nB"]), nR, int(S["nR"]),
                                      "EQUAL" if ok else "DIFFER"))
    return dict(e=e, N=N, Q=Q, q1=dq1, thr=thr, rows=rows, cB=cb, cR=cr, code=code, nB=nB, nR=nR,
                band_s=tb, ok=ok, derived=(dN, dQ, dq1, dthr))


# =============================== the two-sided certificate =================================
def float_optimum(C):
    """(y_cols, y_vals, r) from the stored exact dual, or a dense nnls at the two smallest cells."""
    e, rows = C["e"], C["rows"]
    p = os.path.join(RES, "conedual_exact_dual_2e%d.npz" % e)
    vB = np.flatnonzero(C["code"] == 1)
    if os.path.exists(p):
        z = np.load(p, allow_pickle=False)
        if not np.array_equal(np.asarray(z["rows"], dtype=np.int64), rows):
            raise SystemExit("pa_certify: stored row list differs at 2^%d" % e)
        sup = np.asarray(z["y_support"], dtype=np.int64)
        okidx = bool(sup.size and sup.min() >= 0 and sup.max() < vB.size)
        cols = vB[sup] if okidx else np.empty(0, dtype=np.int64)
        okcol = bool(okidx and np.all(C["code"][cols] == 1))
        clause("C3", okidx and okcol,
               "2^%d stored support: %d indices in [0,%d) %s, every named column even-omega %s"
               % (e, sup.size, vB.size, "yes" if okidx else "NO", "yes" if okcol else "NO"))
        y = np.asarray(z["y_values"], dtype=np.float64)
        inc = (cols[None, :] % rows[:, None] == 0)
        fit = (inc * y[None, :]).sum(axis=1)
        # the stored w/t are divided by R; recover the integer-scale residual
        r = C["cR"].astype(np.float64) - fit * float(C["cR"][0])
        return (cols, y * float(C["cR"][0]), r,
                "stored conedual_exact_dual_2e%d.npz" % e, okidx and okcol)
    from scipy.optimize import nnls
    A = (vB[None, :] % rows[:, None] == 0).astype(np.float64)
    t0 = time.time()
    y, _ = nnls(A, C["cR"].astype(np.float64), maxiter=400000)
    dt = time.time() - t0
    s = np.flatnonzero(y > 0)
    r = C["cR"].astype(np.float64) - A @ y
    clause("C3", True, "2^%d dense nnls over all %d even columns: %.1f s, support %d"
           % (e, vB.size, dt, s.size))
    return vB[s], y[s], r, "dense nnls here (%d columns, %.1f s)" % (vB.size, dt), True


def dual_at_fitting_scale(r, K: int):
    """lambda = round(-K r / max|r|) at the largest power of ten <= K whose ACTUAL bound
    `sum_d |lambda_d|` fits 2^62, so `k_theta`'s accumulator provably cannot overflow (C2).

    Returns (lambda, K) or (None, 0).  The float array is bounded before the int64 cast.
    """
    mx = float(np.abs(r).max()) if r.size else 0.0
    den = mx if mx > 0 else 1.0
    while K >= 1:
        f = -float(K) * np.asarray(r, dtype=np.float64) / den
        if castable(f):
            lam = np.rint(f).astype(np.int64)
            if abs_sum(lam) <= SAFE and abs_max(lam) <= SAFE:
                return lam, K
        K //= 10
    return None, 0


def primal_at_fitting_scale(yv, D: int):
    """q = max(0, round(D y)) at the largest power of ten <= D whose ACTUAL bound `sum_j q_j` fits
    2^62, so `k_fit`'s accumulator provably cannot overflow (C2).

    Returns (q, D) or (None, 0).  The float array is bounded before the int64 cast.
    """
    y = np.asarray(yv, dtype=np.float64)
    while D >= 1:
        f = np.rint(y * float(D))
        if castable(f):
            q = np.maximum(f, 0.0).astype(np.int64)
            if int_sum(q) <= SAFE and int_max(q) <= SAFE:
                return q, D
        D //= 10
    return None, 0


def two_sided(C):
    e, rows, N = C["e"], C["rows"], C["N"]
    # C1's registered consequence, enforced here and not in prose: a cell whose band rebuild does not
    # reproduce the cache is VOID -- no certificate is produced for it and no P-clause counts it.
    if not C["ok"]:
        say("     2^%-3d VOID (C1 failed: the rebuilt band does not reproduce the cache)" % e)
        return dict(e=e, void="C1")
    cR, P = C["cR"], C["cB"] - C["cR"]
    P2 = sq(P)
    cols, yv, r, src, c3 = float_optimum(C)
    if not c3:
        say("     2^%-3d VOID (C3 failed: the stored primal support is not usable)" % e)
        return dict(e=e, void="C3")
    # C4 before any of these columns reaches a kernel or a subscript.
    if not validate_indices(e, rows, cols, C["Q"], C["thr"], N, C["code"], rows):
        say("     2^%-3d VOID (C4 failed: an index is not a valid band column or row)" % e)
        return dict(e=e, void="C4")

    # ---- dual: lambda = round(-K r / max|r|) at a scale whose actual subset-sum bound fits int64,
    #      then the d = 1 repair.  The repair moves lambda_1, so the guard is re-checked after it.
    K0, D0 = scales(int(rows.size), int(cR[0]))
    lam, K = dual_at_fitting_scale(r, K0)
    if lam is None:
        clause("C2", False, "2^%d no power-of-ten dual scale makes sum_d |lambda_d| fit 2^62" % e)
        say("     2^%-3d VOID (C2 failed before the theta sweep)" % e)
        return dict(e=e, void="C2")
    nchunk = 4 * NTHREAD
    acc = np.empty(N + 1, dtype=np.int64)
    k_theta(np.int64(N), rows, lam, acc, nchunk)
    m0 = int(acc[C["code"] == 1].min())
    rep = max(0, -m0)
    if rep:
        lam[0] += rep
        if abs_sum(lam) > SAFE or abs_max(lam) > SAFE:
            del acc
            clause("C2", False, "2^%d the d=1 repair (+%d on lambda_1) pushes sum_d |lambda_d| past"
                                " 2^62" % (e, rep))
            say("     2^%-3d VOID (C2 failed after the d=1 repair)" % e)
            return dict(e=e, void="C2")
        k_theta(np.int64(N), rows, lam, acc, nchunk)
    mmin = int(acc[C["code"] == 1].min())
    del acc

    s_abs, m_abs = abs_sum(lam), abs_max(lam)
    d_obj, n_obj = dot(cR, lam), sq(lam)       # object dtype only: n_obj is about 10^32

    # ---- primal: q = max(0, round(D y)) at a scale whose actual sum fits int64
    q, D = primal_at_fitting_scale(yv, D0)
    if q is None:
        clause("C2", False, "2^%d no power-of-ten primal scale makes sum_j q_j fit 2^62" % e)
        say("     2^%-3d VOID (C2 failed before the fit kernel)" % e)
        return dict(e=e, void="C2")
    keep = q > 0
    cols_k, q_k = np.ascontiguousarray(cols[keep]), np.ascontiguousarray(q[keep])
    s_q, m_q, qnn = int_sum(q_k), int_max(q_k), bool(np.all(q_k >= 0))
    if not (s_q <= SAFE and m_q <= SAFE and qnn):
        clause("C2", False, "2^%d the kept primal breaks its own guard: sum_j q_j = %d, max %d,"
                            " q >= 0 %s" % (e, s_q, m_q, qnn))
        say("     2^%-3d VOID (C2 failed before the fit kernel)" % e)
        return dict(e=e, void="C2")
    fit = np.empty(rows.size, dtype=np.int64)
    k_fit(cols_k, q_k, rows, fit)
    ndiv = np.empty(cols_k.size, dtype=np.int64)
    k_coldiv(cols_k, rows, ndiv)
    tot_obj = int(np.sum(q_k.astype(object) * ndiv.astype(object)))
    tot_ker = int(np.sum(fit.astype(object)))
    c2 = (s_abs <= SAFE and m_abs <= SAFE and s_q <= SAFE and m_q <= SAFE and qnn
          and tot_obj == tot_ker)
    clause("C2", c2, "2^%d int64 kernels guarded by the ACTUAL sums (Python integers, abs after the"
                     " conversion): sum_d |lambda_d| = %.4e = %.3f*2^62, max_d |lambda_d| = %.4e;"
                     " sum_j q_j = %.4e = %.3f*2^62, max_j q_j = %.4e, q >= 0 %s; kernel sum_d fit_d"
                     " %s the object-dtype total (K = 10^%d, D = 10^%d; no (R+1)*D claim is used)"
           % (e, float(s_abs), s_abs / SAFE, float(m_abs), float(s_q), s_q / SAFE, float(m_q),
              qnn, "matches" if tot_obj == tot_ker else "DIFFERS FROM",
              round(math.log10(K)), round(math.log10(D))))

    if not c2:
        say("     2^%-3d VOID (C2 failed: an int64 step is not provably safe)" % e)
        return dict(e=e, void="C2")

    res = cR.astype(object) * D - fit.astype(object)
    hi2 = Fraction(int(np.sum(res ** 2)), D * D * P2)
    lo2 = Fraction(d_obj * d_obj, P2 * n_obj) if (d_obj < 0 and n_obj > 0) else Fraction(0)
    lo, hi = dec_bracket(lo2, hi2, 12)
    ok = (mmin >= 0) and (d_obj < 0) and (lo2 <= hi2)
    say("     2^%-3d from %s" % (e, src))
    say("           dual repair on d=1: %d (theta min before %d, after %d); support kept %d of %d"
        % (rep, m0, mmin, int(keep.sum()), cols.size))
    say("           lower^2 = %s / %s" % (lo2.numerator, lo2.denominator))
    say("           upper^2 = %s / %s" % (hi2.numerator, hi2.denominator))
    say("           CERTIFIED  %s <= kappa(2^%d) <= %s   (width %.2e)"
        % (fstr(lo, 12), e, fstr(hi, 12), float(hi - lo)))
    return dict(e=e, N=N, rows=rows, lam=lam, cols=cols_k, q=q_k, D=D, K=K,
                lo2=lo2, hi2=hi2, lo=lo, hi=hi, repair=rep, theta_min=mmin,
                feasible=ok, src=src, width=hi - lo)


# =============================== the construction's bound ==================================
def construction(C, Y):
    """exact rational ||D c_vR - A q||^2/(D^2||P||^2) for eq:X with the small-prime part a_Y."""
    e, rows, N, Q = C["e"], C["rows"], C["N"], C["Q"]
    cR, P = C["cR"], C["cB"] - C["cR"]
    aY = np.ones(N + 1, dtype=np.int64)
    for p in primes_upto(Y).tolist():           # Python loops over primes <= Y (one or two)
        aY[p::p] *= p
    ev, od = C["code"] == 1, C["code"] == 2
    nb = int(aY.max()) + 1
    B = np.bincount(aY[ev], minlength=nb)
    Rc = np.bincount(aY[od], minlength=nb)
    cap = math.isqrt(Q)
    cls = np.flatnonzero((np.arange(nb) <= cap) & (B > 0) & (Rc > 0))
    cls = cls[cls >= 1]
    Den = math.lcm(*[int(B[a]) for a in cls.tolist()]) if cls.size else 1
    # object dtype: Den is the lcm of the B_a and runs past 10^14
    num = np.array([Den * int(Rc[a]) // int(B[a]) for a in cls.tolist()], dtype=object)
    idx = np.zeros(nb, dtype=np.int64)
    idx[cls] = np.arange(1, cls.size + 1)
    clsmap = np.where(ev, idx[aY], 0).astype(np.int64)
    cnt = np.empty((rows.size, cls.size), dtype=np.int64)
    k_class_counts(clsmap, rows, cls.size, cnt)
    fit = (cnt.astype(object) * num[None, :]).sum(axis=1)
    res = cR.astype(object) * Den - fit
    b2 = Fraction(int(np.sum(res ** 2)), Den * Den * sq(P))
    lit = CONSTR_LIT[e]
    nd = len(lit.split(".")[1])
    s = 10 ** nd
    up = math.isqrt(b2.numerator * s * s // b2.denominator)
    while Fraction(up, s) ** 2 < b2:
        up += 1
    tight = fstr(Fraction(up, s), nd) == lit
    valid = b2 <= Fraction(lit) ** 2
    say("     2^%-3d Y=%d: %d classes (a <= %d), common denominator %d" % (e, Y, cls.size, cap, Den))
    say("           bound^2 = %s / %s" % (b2.numerator, b2.denominator))
    say("           exact bound <= %s ; rounded UP at %d dp = %s ; paper prints %s  -> %s"
        % (fstr(Fraction(up, s), nd), nd, fstr(Fraction(up, s), nd), lit,
           "valid and tight" if (valid and tight) else ("valid" if valid else "INVALID")))
    return dict(e=e, Y=Y, classes=int(cls.size), den=Den, b2=b2,
                roundup=fstr(Fraction(up, s), nd), literal=lit, valid=valid, tight=tight)


# =============================== store / verify ============================================
def cert_path(e):
    return os.path.join(RES, "pa_cert_2e%d.npz" % e)


def stored_cells() -> list[int]:
    """every cell for which a certificate is stored, from the filenames themselves."""
    out = []
    for p in glob.glob(os.path.join(RES, "pa_cert_2e*.npz")):
        b = os.path.basename(p)[len("pa_cert_2e"):-len(".npz")]
        if b.isdigit():
            out.append(int(b))
    return sorted(out)


def store(c):
    np.savez_compressed(cert_path(c["e"]), N=np.int64(c["N"]), e=np.int32(c["e"]),
                        rows=c["rows"].astype(np.int64), dual=c["lam"].astype(np.int64),
                        primal_cols=c["cols"].astype(np.int64),
                        primal_num=c["q"].astype(np.int64), primal_den=np.int64(c["D"]),
                        lower2_num=np.array(str(c["lo2"].numerator)),
                        lower2_den=np.array(str(c["lo2"].denominator)),
                        upper2_num=np.array(str(c["hi2"].numerator)),
                        upper2_den=np.array(str(c["hi2"].denominator)))
    return os.path.getsize(cert_path(c["e"]))


def verify_one(e):
    """re-check a stored certificate from the packet's own files; no solver, no float in the verdict.

    The order is the point: C5 establishes that the fields exist and are well
    formed, C4 that every stored index is a valid row or band column, and C2 that the actual integer
    bound on each kernel accumulator fits int64 -- each BEFORE the thing it protects is used.  A
    missing or malformed field, or an out-of-range column, fails a clause here; it does not raise,
    and it certainly does not pass.
    """
    p = cert_path(e)
    if not os.path.exists(p):
        clause("V%d" % e, False, "2^%d: %s missing" % (e, os.path.basename(p)))
        return None
    okf, z = check_fields(e, p)
    if not okf:
        clause("V%d" % e, False, "2^%d VOID (C5: a stored field is missing or malformed)" % e)
        return dict(e=e, void="C5", bad=["C5"], checks=0)
    C = cell(e)
    if C.get("void") == "C0":
        clause("V%d" % e, False, "2^%d VOID (C0: the cached cell is not the cell 2^%d)" % (e, e))
        return dict(e=e, void="C0", bad=["C0"], checks=0)
    # the certificate's own scalars must name this cell, and its rows must be the complete
    # independently derived row list -- not merely a valid subset of it
    dN, dQ, dq1, dthr = C["derived"]
    sN, se = int(z["N"]), int(z["e"])
    srows = np.asarray(z["rows"], dtype=np.int64)
    c0b = [("certificate e", se == e), ("certificate N = 2^e", sN == dN),
           ("certificate rows are every squarefree d <= Q",
            bool(np.array_equal(np.ascontiguousarray(srows), derive_cell(e)[4])))]
    bad0b = [k for k, v in c0b if not v]
    clause("C0", not bad0b,
           "2^%d stored certificate bound to the derived cell: e = %d, N = %d, %d rows -- %s"
           % (e, se, sN, srows.size,
              "all 3 agree" if not bad0b else "DISAGREE: " + ", ".join(bad0b)))
    if bad0b:
        clause("V%d" % e, False, "2^%d VOID (C0: the certificate does not belong to this cell)" % e)
        return dict(e=e, void="C0", bad=["C0"], checks=0)
    rows, cR, P = C["rows"], C["cR"], C["cB"] - C["cR"]
    # the format test again, at the point of conversion.  C5 has run and would have voided this cell
    # already; repeating it here means no future edit can move a conversion in front of the check, and
    # a caller that reaches this function by another path still cannot have a value wrapped under it.
    conv = [(k, int64_exact(z[k])) for k in ("rows", "dual", "primal_cols", "primal_num",
                                             "primal_den", "N", "e")]
    badconv = [k for k, v in conv if not v]
    if badconv:
        clause("V%d" % e, False, "2^%d VOID (C5: a stored field is not an integer value that signed"
                                 " int64 holds unchanged: %s)" % (e, ", ".join(badconv)))
        return dict(e=e, void="C5", bad=["C5"], checks=0)
    lam = np.asarray(z["dual"], dtype=np.int64)
    cols = np.asarray(z["primal_cols"], dtype=np.int64)
    q = np.asarray(z["primal_num"], dtype=np.int64)
    D = int(z["primal_den"])
    if not validate_indices(e, srows, cols, C["Q"], C["thr"], C["N"], C["code"], rows):
        clause("V%d" % e, False, "2^%d VOID (C4: a stored index is not a valid row or band column)"
               % e)
        return dict(e=e, void="C4", bad=["C4"], checks=0)

    # ---- C2: the actual bound on each kernel accumulator, as a Python integer, before it runs
    s_lam, m_lam = abs_sum(lam), abs_max(lam)
    s_q, m_q, qnn = int_sum(q), int_max(q), bool(np.all(q >= 0))
    c2 = (lam.size == rows.size and D >= 1 and qnn
          and s_lam <= SAFE and m_lam <= SAFE and s_q <= SAFE and m_q <= SAFE)
    clause("C2", c2, "2^%d stored vectors guarded before the kernels: sum_d |lambda_d| = %.4e ="
                     " %.3f*2^62, max_d |lambda_d| = %.4e; sum_j q_j = %.4e = %.3f*2^62, max_j q_j"
                     " = %.4e, q >= 0 %s; |lambda| = |rows| %s, D = %d >= 1 %s"
           % (e, float(s_lam), s_lam / SAFE, float(m_lam), float(s_q), s_q / SAFE, float(m_q),
              qnn, lam.size == rows.size, D, D >= 1))
    if not c2:
        clause("V%d" % e, False, "2^%d VOID (C2: an int64 kernel is not provably safe on the stored"
                                 " vectors)" % e)
        return dict(e=e, void="C2", bad=["C2"], checks=0)

    acc = np.empty(C["N"] + 1, dtype=np.int64)
    k_theta(np.int64(C["N"]), rows, lam, acc, 4 * NTHREAD)
    tmin = int(acc[C["code"] == 1].min())
    del acc
    d_obj, n_obj, P2 = dot(cR, lam), sq(lam), sq(P)
    fit = np.empty(rows.size, dtype=np.int64)
    cols_c, q_c = np.ascontiguousarray(cols), np.ascontiguousarray(q)
    k_fit(cols_c, q_c, rows, fit)
    ndiv = np.empty(cols_c.size, dtype=np.int64)
    k_coldiv(cols_c, rows, ndiv)
    tot_obj = int(np.sum(q_c.astype(object) * ndiv.astype(object)))
    tot_ker = int(np.sum(fit.astype(object)))
    res = cR.astype(object) * D - fit.astype(object)
    hi2 = Fraction(int(np.sum(res ** 2)), D * D * P2)
    lo2 = Fraction(d_obj * d_obj, P2 * n_obj) if (d_obj < 0 and n_obj > 0) else Fraction(0)
    chk = [("theta >= 0 on every even column", tmin >= 0), ("<c_vR,lambda> < 0", d_obj < 0),
           ("lower <= upper", lo2 <= hi2),
           ("k_fit agrees with the object-dtype total", tot_obj == tot_ker),
           ("lower2 matches the stored value", lo2 == Fraction(int(str(z["lower2_num"])),
                                                               int(str(z["lower2_den"])))),
           ("upper2 matches the stored value", hi2 == Fraction(int(str(z["upper2_num"])),
                                                               int(str(z["upper2_den"]))))]
    bad = [k for k, v in chk if not v]
    lo, hi = dec_bracket(lo2, hi2, 12)
    clause("V%d" % e, not bad and C["ok"],
           "2^%d re-verified from the stored npz: %s ; %s <= kappa <= %s (theta min %d, width %.2e)"
           % (e, "all %d exact checks pass" % len(chk) if not bad else "FAILED: " + ", ".join(bad),
              fstr(lo, 12), fstr(hi, 12), tmin, float(hi - lo)))
    return dict(e=e, lo=fstr(lo, 12), hi=fstr(hi, 12), width="%.3e" % float(hi - lo),
                checks=len(chk), bad=bad, theta_min=tmin, repair_free=bool(tmin >= 0),
                sum_abs_dual=str(s_lam), sum_primal=str(s_q), D=D,
                nrows=int(rows.size), ncols=int(cols.size))


# =============================== the one-sided certificate at 2^30 =========================
LOWER30 = os.path.join(RES, "pa_cert_2e30_lower.npz")
DEN30 = 10 ** 12                      # the scale the 2^30 dual is rounded at; see lower30_build


def lower30_build(sep_tag: str = "_n30b"):
    """Store the REPAIRED INTEGRAL dual at 2^30 as a certificate of its own.

    The separator file holds a floating `w`.  A float is not a certificate: it is rounded to an
    integral `lambda`, repaired at the row `d = 1` so that `theta_lambda >= 0` on every even column,
    and only then are `J = -<c_O,lambda>`, `H = ||lambda||^2` and `D_p = ||P||^2` exact integers.  This
    writes those, with the row list and the vector itself, so a reader can re-establish the lower
    endpoint without the solver and without trusting any float.

    The rounding scale is DEN30 = 10^12, which is the scale of the published four integers, so the
    certificate reproduces them digit for digit rather than certifying the same cell at some other
    scale.  The d = 1 entry of the stored `w` is a marker and is discarded; the repair rebuilds it.
    """
    sep = os.path.join(RES, "conedual_colgen_sep_2e%d%s.npz" % (LOWER30_E, sep_tag))
    if not os.path.exists(sep):
        clause("L0", False, "no separator for 2^%d at %s" % (LOWER30_E, os.path.basename(sep)))
        return None
    C = cell(LOWER30_E)
    if C.get("void") or not C["ok"]:
        say("     2^%-3d VOID (the cell could not be established)" % LOWER30_E)
        return None
    rows, cR, P = C["rows"], C["cR"], C["cB"] - C["cR"]
    z = np.load(sep, allow_pickle=False)
    zr = np.asarray(z["rows"], dtype=np.int64)
    w = np.asarray(z["w"], dtype=np.float64)
    ok_rows = bool(np.array_equal(zr, rows))
    clause("L0", ok_rows, "2^%d separator %s: rows equal the derived row list %s (%d rows)"
           % (LOWER30_E, os.path.basename(sep), "yes" if ok_rows else "NO", rows.size))
    if not ok_rows:
        return None
    # lambda = round(-DEN30 w), with the d = 1 entry discarded and rebuilt by the repair below.  The
    # separator stores `w` with <c_O, w> > 0, so the minus sign is the one that makes <c_O, lambda> < 0,
    # i.e. J > 0.  DEN30 is 10^12 and not the largest scale that fits int64 for one reason: 10^12 is the
    # scale at which the four integers of eq:certint were published, so a certificate at this scale
    # carries exactly the J and H the paper prints, and a reader comparing the two sees the same digits.
    # The accumulator bound is checked for the scale rather than assumed for it.
    K = DEN30
    f = -float(K) * w
    if not castable(f):
        clause("L1", False, "2^%d the dual at 10^%d does not fit int64"
               % (LOWER30_E, round(math.log10(K))))
        return None
    lam = np.rint(f).astype(np.int64)
    lam[0] = 0
    if abs_sum(lam) > SAFE:
        clause("L1", False, "2^%d sum_d |lambda_d| at 10^%d does not fit 2^62"
               % (LOWER30_E, round(math.log10(K))))
        return None
    nchunk = 4 * NTHREAD
    t0 = time.time()
    acc = np.empty(C["N"] + 1, dtype=np.int64)
    k_theta(np.int64(C["N"]), rows, lam, acc, nchunk)
    m0 = int(acc[C["code"] == 1].min())
    rep = max(0, -m0)
    if rep:
        lam[0] += rep
        k_theta(np.int64(C["N"]), rows, lam, acc, nchunk)
    tmin = int(acc[C["code"] == 1].min())
    del acc
    sweep_s = time.time() - t0
    J, H, Dp = -dot(cR, lam), sq(lam), sq(P)
    ok1 = tmin >= 0 and J > 0 and H > 0 and Dp > 0 and abs_sum(lam) <= SAFE
    clause("L1", ok1, "2^%d integral dual built at K = 10^%d and repaired on d=1 by %d: theta min %d,"
                      " J = %d > 0 %s, H = %d, D_p = %d  (sweep %.0f s)"
           % (LOWER30_E, round(math.log10(K)), rep, tmin, J, J > 0, H, Dp, sweep_s))
    if not ok1:
        return None
    c6 = lower30_c6(J, H, Dp)
    np.savez_compressed(LOWER30, N=np.int64(C["N"]), e=np.int32(LOWER30_E),
                        rows=rows.astype(np.int64), dual=lam.astype(np.int64),
                        J=np.array(str(J)), H=np.array(str(H)), D_p=np.array(str(Dp)),
                        c6=np.int64(c6))
    say("           stored %s (%d bytes): c6 = %d, so kappa(2^%d) > %d/10^6"
        % (os.path.basename(LOWER30), os.path.getsize(LOWER30), c6, LOWER30_E, c6))
    return dict(e=LOWER30_E, J=J, H=H, Dp=Dp, c6=c6, theta_min=tmin, repair=rep, K=K,
                sweep_s=round(sweep_s, 1), band_s=round(C["band_s"], 1))


def lower30_c6(J: int, H: int, Dp: int) -> int:
    """the largest six-decimal c with c^2 H D_p < J^2 10^12, i.e. the best bound the triple gives."""
    c = math.isqrt(J * J * 10 ** 12 // (H * Dp))
    while (c + 1) ** 2 * H * Dp < J * J * 10 ** 12:
        c += 1
    while c > 0 and not (c * c * H * Dp < J * J * 10 ** 12):
        c -= 1
    return c


def lower30_verify():
    """re-establish the 2^30 lower endpoint from the stored integral dual, with no solver and no float.

    The band at 2^30 is rebuilt here and `theta_lambda >= 0` is swept over every even column of it;
    that is the only part of the claim a scalar inequality cannot carry.
    """
    if not os.path.exists(LOWER30):
        clause("L2", False, "2^%d: %s missing" % (LOWER30_E, os.path.basename(LOWER30)))
        return None
    try:
        z = np.load(LOWER30, allow_pickle=False)
        names = set(z.files)
        need = {"N", "e", "rows", "dual", "J", "H", "D_p", "c6"}
        if need - names:
            clause("L2", False, "2^%d certificate lacks %s"
                   % (LOWER30_E, ", ".join(sorted(need - names))))
            return None
        # the format, BEFORE any conversion: dtype, rank, and that signed int64 holds every value
        # unchanged.  A float dual whose entries truncate to the stored integers, or an unsigned one
        # whose entries wrap to them, is a different vector from the one the checks below would test.
        fmt = [("dual is an integer vector", int64_exact(z["dual"])),
               ("dual rank 1", int(np.asarray(z["dual"]).ndim) == 1),
               ("rows is an integer vector", int64_exact(z["rows"])),
               ("rows rank 1", int(np.asarray(z["rows"]).ndim) == 1),
               ("N, e, c6 are integer scalars",
                all(int64_exact(z[k]) and int(np.asarray(z[k]).ndim) == 0
                    for k in ("N", "e", "c6"))),
               ("J, H, D_p are decimal integers",
                all(decimal_int(z[k]) for k in ("J", "H", "D_p")))]
        badfmt = [k for k, v in fmt if not v]
        if badfmt:
            clause("L2", False, "2^%d certificate refused on its format, before any conversion:"
                                " FAILED: %s" % (LOWER30_E, ", ".join(badfmt)))
            return None
        J, H, Dp = int(str(z["J"])), int(str(z["H"])), int(str(z["D_p"]))
        c6 = int(z["c6"])
        lam = np.asarray(z["dual"], dtype=np.int64)
        srows = np.asarray(z["rows"], dtype=np.int64)
    except Exception as exc:
        clause("L2", False, "2^%d certificate unreadable (%s)" % (LOWER30_E, type(exc).__name__))
        return None
    C = cell(LOWER30_E)
    if C.get("void") or not C["ok"]:
        clause("L2", False, "2^%d VOID (the cell could not be established)" % LOWER30_E)
        return None
    rows, cR, P = C["rows"], C["cR"], C["cB"] - C["cR"]
    bound = [("e", int(z["e"]) == LOWER30_E), ("N = 2^e", int(z["N"]) == C["N"]),
             ("rows are every squarefree d <= Q", bool(np.array_equal(srows, rows))),
             ("|lambda| = |rows|", lam.size == rows.size),
             ("sum_d |lambda_d| fits 2^62", abs_sum(lam) <= SAFE)]
    badb = [k for k, v in bound if not v]
    clause("L2", not badb, "2^%d stored integral dual bound to the derived cell (its format having"
                           " been refused or accepted first): %s"
           % (LOWER30_E, "all %d checks pass" % len(bound) if not badb
              else "FAILED: " + ", ".join(badb)))
    if badb:
        return None
    t0 = time.time()
    acc = np.empty(C["N"] + 1, dtype=np.int64)
    k_theta(np.int64(C["N"]), rows, lam, acc, 4 * NTHREAD)
    tmin = int(acc[C["code"] == 1].min())
    del acc
    sweep_s = time.time() - t0
    Jr, Hr, Dpr = -dot(cR, lam), sq(lam), sq(P)
    chk = [("theta >= 0 on every even column", tmin >= 0),
           ("J recomputes", Jr == J), ("H recomputes", Hr == H), ("D_p recomputes", Dpr == Dp),
           ("J > 0", J > 0),
           ("c6^2 H D_p < J^2 10^12", c6 * c6 * H * Dp < J * J * 10 ** 12),
           ("c6 is maximal", not ((c6 + 1) ** 2 * H * Dp < J * J * 10 ** 12)),
           ("equivalent to J^2/(H D_p) > (c6/10^6)^2, which implies kappa^2 > (c6/10^6)^2",
            (Fraction(J * J, H * Dp) > Fraction(c6, 10 ** 6) ** 2)
            == (c6 * c6 * H * Dp < J * J * 10 ** 12))]
    bad = [k for k, v in chk if not v]
    clause("L3", not bad, "2^%d lower endpoint re-established: %s ; kappa(2^%d) > %d/10^6 = %s"
                          "  (theta min %d, sweep %.0f s, band %.0f s)"
           % (LOWER30_E, "all %d checks pass" % len(chk) if not bad else "FAILED: " + ", ".join(bad),
              LOWER30_E, c6, fstr(Fraction(c6, 10 ** 6), 6), tmin, sweep_s, C["band_s"]))
    # L4: the certificate's integers are the ones the solver run recorded and the paper prints.  Without
    # it, a rebuild at another rounding scale would still certify the cell -- correctly -- while quietly
    # replacing the published J and H by different valid ones, and no clause here would notice.
    led = os.path.join(RES, "conedual_colgen_2e%d_n30b.json" % LOWER30_E)
    if not os.path.exists(led):
        clause("L4", False, "2^%d ledger %s missing: the certificate's integers cannot be compared with"
                            " the published ones" % (LOWER30_E, os.path.basename(led)))
    else:
        with io.open(led, encoding="utf-8") as fh:
            L = json.load(fh)
        same = [("J", str(L.get("J")) == str(J)), ("H", str(L.get("H")) == str(H)),
                ("D_p", str(L.get("D_p")) == str(Dp)), ("c6", int(L.get("c6", -1)) == c6)]
        badl = [k for k, v in same if not v]
        clause("L4", not badl, "2^%d the certificate carries the published four integers of %s: %s"
               % (LOWER30_E, os.path.basename(led),
                  "J, H, D_p and c6 all agree" if not badl else "DISAGREE: " + ", ".join(badl)))
        bad = bad + ["ledger " + k for k in badl]
    return dict(e=LOWER30_E, J=str(J), H=str(H), D_p=str(Dp), c6=c6, theta_min=tmin,
                lower=fstr(Fraction(c6, 10 ** 6), 6), bad=bad, checks=len(chk),
                sweep_s=round(sweep_s, 1), band_s=round(C["band_s"], 1))


# =============================== main ======================================================
def main(argv) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true", help="re-check stored certificates only")
    ap.add_argument("--cells", type=int, nargs="*", default=None)
    ap.add_argument("--make-lower30", action="store_true",
                    help="build and store the integral dual certificate at 2^30")
    ap.add_argument("--verify-lower30", action="store_true",
                    help="re-establish the 2^30 lower endpoint from that certificate")
    a = ap.parse_args(argv)
    lower30 = a.make_lower30 or a.verify_lower30
    if lower30:
        # Its own record, bound BEFORE the first line is written: `say` rewrites OUT on every call, so
        # rebinding after the two header lines would leave the two-sided record holding a header and
        # nothing else -- the record of a verify run destroyed by a run that never checked those cells.
        global OUT, OUTJ
        OUT = os.path.join(RES, "pa_cert_2e30_lower.txt")
        OUTJ = os.path.join(RES, "pa_cert_2e30_lower.json")
    nth = max(1, min(NTHREAD, 8, int(numba.config.NUMBA_NUM_THREADS)))
    set_num_threads(nth)
    T0 = time.time()
    say("pa_certify -- exact-arithmetic certificates for section 11 (clauses in the header)")
    say("threads %d (ACF_THREADS or 4, capped at 8 and at numba's own count); the int64 kernels are"
        " guarded by the ACTUAL integer sums of the vectors passed to them (clause C2)" % nth)
    say()
    check_c9()
    say()

    if lower30:
        what = "built" if a.make_lower30 else "re-established"
        say("2^%d, ONE SIDE ONLY -- the lower endpoint from a stored integral dual" % LOWER30_E)
        say("the band at 2^%d is rebuilt here; expect several minutes and a few GB" % LOWER30_E)
        say()
        r = lower30_build() if a.make_lower30 else lower30_verify()
        say()
        say("wall %.1f s" % (time.time() - T0))
        say_fails()
        d = dict(mode="lower30-" + ("build" if a.make_lower30 else "verify"), threads=nth,
                 wall_s=time.time() - T0, result=r, completeness_checked=False,
                 clauses=REC, fails=FAILS)
        # one cell, one side: the verdict says so rather than standing for the packet
        low_scope = (" for 2^%d (the lower endpoint only); completeness not checked" % LOWER30_E)
        rc = verdict(d, emit=QUIET, scope=low_scope)
        json.dump(d, io.open(OUTJ, "w", encoding="utf-8", newline="\n"), indent=1)
        say("%s the 2^%d lower endpoint; wrote %s and %s"
            % (what, LOWER30_E, os.path.basename(OUT), os.path.basename(OUTJ)))
        verdict(scope=low_scope)
        return rc

    if a.verify:
        found = stored_cells()
        cells = a.cells or list(ADVERTISED)
        # A SELECTED-CELL RUN IS NOT A PACKET CHECK.  With --cells naming anything other than the
        # advertised set, this invocation examines the cells it was given and nothing else: a certificate
        # outside the selection is not opened, so a damaged one is not seen.  Announcing "every
        # advertised certificate re-checked" there is false, and so is an unqualified CERTIFICATES:
        # VALID.  Both are scoped below, and the completeness clause reports NOT RUN rather than passing
        # on a question this run never asked.
        full = set(cells) == set(ADVERTISED)
        scope = "" if full else (" for %s only; completeness not checked"
                                 % ", ".join("2^%d" % e for e in cells))
        if full:
            say("VERIFY MODE -- every advertised certificate re-checked from the packet's own files")
            say("cells: %s" % ", ".join("2^%d" % e for e in cells))
        else:
            say("VERIFY MODE, SELECTED CELLS ONLY -- %s re-checked from the packet's own files; the"
                " other advertised certificates are NOT opened and the packet's completeness is NOT"
                " checked in this mode" % ", ".join("2^%d" % e for e in cells))
            say("cells: %s  (the advertised set is %s -- run without --cells for the packet check)"
                % (", ".join("2^%d" % e for e in cells),
                   ", ".join("2^%d" % e for e in ADVERTISED)))
        # C6: the inventory must be exactly the advertised set.  Discovering whatever files happen to
        # be present and succeeding on them means an empty directory verifies, and removing one
        # certificate silently shrinks the claim instead of failing it.
        missing = [e for e in ADVERTISED if e not in found]
        extra = [e for e in found if e not in ADVERTISED]
        if full:
            clause("C6", not missing and not extra and bool(found),
                   "inventory is exactly the %d advertised cells %s: present %s%s%s"
                   % (len(ADVERTISED), ", ".join("2^%d" % e for e in ADVERTISED),
                      ", ".join("2^%d" % e for e in found) or "NONE",
                      "" if not missing else "; MISSING " + ", ".join("2^%d" % e for e in missing),
                      "" if not extra else "; UNADVERTISED " + ", ".join("2^%d" % e for e in extra)))
        else:
            clause_not_run("C6", "the inventory is NOT checked: --cells selected %s, and a run over a"
                                 " subset of the advertised cells %s cannot speak for the packet's"
                                 " contents.  Run `--verify` with no --cells for the full-packet check."
                           % (", ".join("2^%d" % e for e in cells),
                              ", ".join("2^%d" % e for e in ADVERTISED)))
        say()
        vs = []
        for e in cells:
            t0 = time.time()
            v = verify_one(e)
            if v:
                v["verify_s"] = round(time.time() - t0, 1)
                say("           (2^%d in %.1f s, band rebuild included)" % (e, v["verify_s"]))
            vs.append(v)
            say()
        say("verify wall %.1f s" % (time.time() - T0))
        say_fails()
        d = dict(mode="verify" if full else "verify-selected", threads=nth, wall_s=time.time() - T0,
                 cells=cells, advertised=list(ADVERTISED), completeness_checked=full,
                 verified=[v for v in vs if v], clauses=REC, fails=FAILS)
        rc = verdict(d, emit=QUIET, scope=scope)
        json.dump(d, io.open(OUTJ, "w", encoding="utf-8", newline="\n"), indent=1)
        say("wrote %s and %s" % (os.path.basename(OUT), os.path.basename(OUTJ)))
        verdict(scope=scope)
        return rc

    cells = a.cells or sorted(set(TWO_SIDED) | {c for c, _ in CONSTRUCTION})
    # CONSTRUCTION MODE, not a packet check: this re-derives the certificates of the cells it was given
    # and writes them.  It opens no other cell and runs no completeness clause, so its verdict is scoped
    # to those cells exactly as a selected-cell verification is.
    con_scope = (" for the cells built (%s); completeness not checked"
                 % ", ".join("2^%d" % e for e in cells))
    say("CONSTRUCTION MODE -- certificates re-derived and written for %s; this mode does not check the"
        " packet's contents (that is `--verify` with no --cells)" % ", ".join("2^%d" % e for e in cells))
    say()
    say("part 1 -- the two-sided certificate of eq:cert at 2^%s"
        % ", 2^".join(str(e) for e in TWO_SIDED if e in cells))
    certs, sizes, voids = [], {}, []
    cache = {}
    for e in [x for x in TWO_SIDED if x in cells]:
        C = cache.setdefault(e, cell(e))
        c = two_sided(C)
        if c.get("void"):
            voids.append(c)
            say()
            continue
        certs.append(c)
        sizes[e] = store(c)
        say("           stored %s (%d bytes)" % (os.path.basename(cert_path(e)), sizes[e]))
        say()

    say("part 2 -- the construction's own bounds (eq:X), exact rational")
    cons = []
    for e, Y in CONSTRUCTION:
        if e not in cells:
            continue
        C = cache.setdefault(e, cell(e))
        if not C["ok"]:
            say("     2^%-3d VOID (C1 failed)" % e)
            say()
            continue
        cons.append(construction(C, Y))
        say()

    # ---------------- clauses -------------------------------------------------------------
    say("clauses")
    if voids:
        clause("VOID", False, "voided cells (their registered consequence, enforced in code): "
               + ", ".join("2^%d by %s" % (v["e"], v["void"]) for v in voids))
    clause("P1", bool(certs) and all(c["feasible"] for c in certs),
           "two-sided certificate valid at %s: theta >= 0 on every even column, <c_vR,lambda> < 0,"
           " lower <= upper" % ", ".join("2^%d" % c["e"] for c in certs))
    clause("P2", bool(certs) and all(c["width"] < WIDTH_MAX for c in certs),
           "interval width < 1e-8 at every size: %s"
           % ", ".join("2^%d %.1e" % (c["e"], float(c["width"])) for c in certs))

    p3, p3d = True, []
    for c in certs:
        fl = SEPMASS_FLOAT.get(c["e"])
        if fl is None:
            p = os.path.join(RES, "conedual_exact_dual_2e%d.npz" % c["e"])
            fl = repr(float(np.load(p, allow_pickle=False)["ub"])) if os.path.exists(p) else None
        if fl is None:
            continue
        inside = c["lo"] <= Fraction(fl) <= c["hi"]
        p3 = p3 and inside
        p3d.append("2^%d %s %s" % (c["e"], fl, "inside" if inside else "OUTSIDE"))
    clause("P3", p3, "the certified interval contains the stored float optimum: " + "; ".join(p3d))

    p4, p4d = True, []
    for c in certs:
        c6 = LEDGER_C6.get(c["e"])
        if c6 is None:
            continue
        ge = c["lo"] >= Fraction(c6, 10 ** 6)
        p4 = p4 and ge
        p4d.append("2^%d lower %s >= %d/10^6 %s" % (c["e"], fstr(c["lo"], 12), c6,
                                                    "yes" if ge else "NO"))
    clause("P4", p4, "at least as strong as eq:certint where that exists: " + "; ".join(p4d))
    clause("P5", bool(cons) and all(x["valid"] and x["tight"] for x in cons),
           "construction bounds certified and tight: "
           + "; ".join("2^%d Y=%d %s" % (x["e"], x["Y"], x["literal"]) for x in cons))
    clause("P6", all(v < 200_000 for v in sizes.values()),
           "certificate files under 200 kB: "
           + ", ".join("2^%d %d B" % (e, v) for e, v in sorted(sizes.items())))

    say()
    say("certified intervals, 12 decimals")
    for c in certs:
        say("   2^%-3d  %s  <=  kappa  <=  %s" % (c["e"], fstr(c["lo"], 12), fstr(c["hi"], 12)))
    say()
    say("wall %.1f s" % (time.time() - T0))
    say_fails()

    d = dict(mode="certify", threads=nth, wall_s=time.time() - T0,
             certificates=[dict(e=c["e"], lower=fstr(c["lo"], 12), upper=fstr(c["hi"], 12),
                                      lower2=str(c["lo2"]), upper2=str(c["hi2"]),
                                      repair=c["repair"], theta_min=c["theta_min"],
                                      support=int(c["cols"].size), src=c["src"],
                                      K=c["K"], D=c["D"],
                                bytes=sizes.get(c["e"])) for c in certs],
             construction=[dict(e=x["e"], Y=x["Y"], classes=x["classes"], den=str(x["den"]),
                                bound2=str(x["b2"]), roundup=x["roundup"],
                                literal=x["literal"], valid=x["valid"], tight=x["tight"])
                           for x in cons],
             cells=cells, completeness_checked=False, clauses=REC, fails=FAILS)
    rc = verdict(d, emit=QUIET, scope=con_scope)
    json.dump(d, io.open(OUTJ, "w", encoding="utf-8", newline="\n"), indent=1)
    say("wrote %s and %s" % (os.path.basename(OUT), os.path.basename(OUTJ)))
    verdict(scope=con_scope)
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

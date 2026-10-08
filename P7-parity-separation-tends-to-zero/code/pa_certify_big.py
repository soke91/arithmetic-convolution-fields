# -*- coding: utf-8 -*-
r"""pa_certify_big -- certify eq:cert at a cell whose primal lives in a colgen separator.

WHAT THIS IS

`pa_certify.py` certifies `eq:cert` at `2^16 ... 2^24`.  At `2^20`/`2^22`/`2^24` it was pure storage-and-rounding:
`conedual_exact_dual_2e{20,22,24}.npz` holds the dual AND the primal support, so `pa_certify.py` reads
both and the certificate falls out.  At `2^26`/`2^28` no such file exists -- only
`conedual_colgen_sep_2e{26,28}*.npz`, which stored the dual alone.

Recovering the primal from the dual by column generation was tried and cannot be done:
the stored separator is REPAIRED (individual rows `w_d` are lowered to force `max_vB W <= tau`),
so its `theta` ordering no longer locates the support.  Measured at `2^18`, where the support is known
from a dense nnls: the 301 support columns are spread to rank 6,196 of 18,685, so a third of the band
would have to be solved to be sure of catching them -- 1.74M columns at `2^26`, 7.06M at `2^28`.

So the other route was taken: `conedual_colgen.py` now SAVES its primal (`y_support`, `y_values`, `ub`),
and this file is the thin step that turns a separator carrying those fields into a certificate.  There
is no solver here and no column generation: the dual direction is the separator's own `w`, the primal
is the one the colgen already found, and everything else is `pa_certify`'s validated machinery, reused
rather than reimplemented.

HOW THE PAIR IS FED TO pa_certify.two_sided

`pa_certify.two_sided` takes `(cols, yv, r, src, ok)` from `float_optimum` and uses `r` ONLY as the
dual direction (`lambda = round(-K r / max|r|)`, then the `d = 1` repair) and `yv` ONLY as the primal
(`q = max(0, round(D yv))`).  The two halves of `eq:cert` share no data, so they may come from
different places, and here they do:

  dual    `r = w`, a separator's repaired dual.  It is already nearly feasible -- that is what the
          repair guarantees -- so the `d = 1` repair costs almost nothing and the lower bound comes
          out at the colgen's own float value, which is better than `eq:certint`'s `c6/10^6`.
          The dual may come from a DIFFERENT separator than the primal (`--dual-tag`), because the
          two halves of `eq:cert` share no data: any dual-feasible `lambda` bounds `kappa` below and
          any nonnegative `q/D` bounds it above, and the stored certificate holds one of each,
          verified independently.  That matters here: a run that stops on its budget can leave an
          excellent primal and a dual weaker than an older converged run's, and then the best
          certificate pairs the old dual with the new primal.
  primal  `cols = vB[y_support]`, `yv = y_values * R` (the colgen solves on the `t = c_vR/R` scale,
          the same convention `conedual_exact_dual_2e*.npz` uses).

Everything exact is then `pa_certify`'s: the band rebuilt as `cell_cache._build` does and controlled
against the cache element for element (C1), the theta sweep over EVERY even column, object-dtype
arithmetic for anything that can exceed int64, and `K`, `D` chosen per cell from a proven int64 bound
(C2).  `--verify` is `pa_certify.verify_one`, unchanged.

MEMORY at `2^28` (N = 268,435,456): the band rebuild peaks at `7N` bytes = 1.88 GB and the theta
accumulator is int64 `N+1` = 2.15 GB, so about 4 GB, inside the grant's 32 GB and far from the 10 GB
intermediate cap (nothing is written but the certificate).

WHERE THE GUARDS LIVE.  The int64 hardening is in `pa_certify.py`, which this file reuses
wholesale, so C2 is restated here to match it and C4/C5 are inherited as well.  `store_if_better`
REFUSES when the stored certificate cannot be read, rather than falling through to the overwrite,
which is the wrong direction for a packet file.  No stored certificate is rewritten unless its
replacement is strictly narrower.

CLAUSES.  No literal below is edited after any output of the revision that
introduced it; every "voids" is enforced in code, not in prose.

  C1  CONTROL, inherited unchanged from `pa_certify.cell`: the rebuilt band has `|vB| = nB` and
      `|vR| = nR`, AND the per-row divisor counts equal the cached `cb`, `cr` ELEMENT FOR ELEMENT.
      CAN FAIL; voids the cell.
  C2  CONTROL, inherited from `pa_certify.two_sided` as restated there: both int64 kernels are
      guarded by the ACTUAL integer bound on their accumulator, computed as a Python integer before
      the kernel runs, with `abs` taken after the conversion to Python integers --
      `sum_d |lambda_d| <= 2^62` and `max_d |lambda_d| <= 2^62` for `k_theta`, `q >= 0`,
      `sum_j q_j <= 2^62` and `max_j q_j <= 2^62` for `k_fit` -- and `sum_d fit_d` is cross-checked
      against the object-dtype total.  No `(R+1) * D` claim is made or used; it is false at every
      stored cell.  CAN FAIL; voids the cell.
  C4  CONTROL, inherited from `pa_certify`: every stored index is validated before it indexes
      anything and independently of `code` -- rows in `[1,Q]`, strictly increasing, squarefree, equal
      to the cell's row list; columns in `(thr,N]`, distinct, squarefree, `P^+ <= Q`, `omega` even.
      CAN FAIL; voids the cell.
  C5  CONTROL, inherited from `pa_certify.verify_one`: a stored certificate's fields are present and
      well formed, and a missing or malformed field fails rather than raising or passing.
      CAN FAIL; voids the cell.
  C6  CONTROL.  The separator is usable and carries a primal: its `rows` equals the cell's `rows`
      element for element, its `t` equals `c_vR/R` exactly, `<c_vR, w> > 0` (so `lambda = -w` has the
      sign `eq:cert` needs), `y_support` is non-empty and indexes `vB` in range, every column it names
      is an even-omega band element, and `y_values >= 0`.  CAN FAIL; voids the cell.

  C9  CONTROL, the same clause `pa_certify.py` runs and through its own function: the float square
      root used as an integer bound equals `math.isqrt` at every integer below `2^31`.  This file runs
      EVERY check the verifier it wraps runs; there is no check it omits.  CAN FAIL.
  C8  CONTROL (`--verify` with no `--cells`).  The inventory is exactly the advertised cells; a missing
      or unadvertised certificate fails, an empty inventory fails.  WITH `--cells` it reports **NOT
      RUN** -- not a pass -- because a run that opens the selected certificates alone cannot speak for
      the packet's contents, and the run's announcement and verdict are scoped to the selection.

  P1  The two-sided certificate holds at every cell run: `theta_lambda(l) >= 0` on EVERY `l in vB`
      (one exact sweep over the whole band, not only over the stored support), `<c_vR,lambda> < 0`,
      and `lower2 <= upper2` as exact rationals.
  P2  At `2^26` and `2^28` the certified interval is narrower than `10^-6`, so the paper can print a
      certified enclosure to six decimals at both.
  P3  REPRODUCTION, ledger integers against the exact lower bound.  The certified lower endpoint is
      `>= c6/10^6` with `c6` read from `conedual_colgen_2e<e>*.json`, so the new certificate is at
      least as strong as `eq:certint` at every cell run.
  P4  REPRODUCTION, stored float against the exact interval.  The certified interval contains the
      `ub` field the separator now stores (the colgen's own float upper bound) OR lies below it; and
      it contains the separator's `bound` field (the colgen's float lower bound) OR lies above it.
      Either way the exact interval and the stored floats are consistent, which is what a
      reproduction clause can say when one side is a float and the other is exact.
  P5  Every certificate file is under 5 MB.
  C7  CONTROL.  A stored certificate is never replaced by one that loses any of its enclosure:
      `pa_cert_2e<e>.npz` is written only if none is there, or this run's `[lower^2, upper^2]` lies
      INSIDE the stored one and is strictly smaller on at least one side -- compared as the exact
      rationals the file stores, not as rounded widths, so two intervals that merely cross leave the
      stored certificate in place.  (An earlier version overwrote the `2^18` certificate, width
      `10^-12`, with a `2.6e-4` one.)  This cannot fail -- it is enforced, not checked -- and the run
      says which happened, with both enclosures.

  NOT REGISTERED: a seven-decimal width.  The width is whatever the colgen's primal and the repaired
  dual give; P2 asks for six decimals and this run's own record states what was reached.

WHAT THE RUN REPORTS, AND WHAT THE EXIT CODE MEANS

  The run ends with the same two machine-readable lines as `pa_certify.py` (whose header explains the
  split), over this file's clauses:

    CERTIFICATES: VALID            -- or INVALID with the failing certification clause named
    DIAGNOSTICS: none              -- or each failing diagnostic with the reason it can fail

  CERTIFICATION here is `C0`-`C9`, the `V<e>` re-verifications, `VOID` and `P1`; DIAGNOSTIC is `P2`
  (a width target), `P3` (the published `eq:certint` integers), `P4` (the separator's stored floats) and
  `P5` (a file-size budget).  THE EXIT CODE REPORTS CERTIFICATE VALIDITY ALONE.  `P2` is expected to fail
  at `2^28`, where the width is `1.2e-5`, and that no longer makes the run exit nonzero; the `FAILS:`
  line still lists every failing clause of either kind, and the JSON record adds `certificates_valid`,
  `certification_fails` and `diagnostic_fails`.  No clause's test or threshold is touched by the split.

  AND WHAT IT COVERS.  `--verify` with no `--cells` is the one invocation that speaks for the packet's
  contents, and it alone prints the bare `CERTIFICATES: VALID` with `C8` scored.  A run that selects
  cells -- `--verify --cells 26 28`, or the construction mode -- examines those cells and nothing else:
  it says so in its announcement, `C8` reports **NOT RUN** instead of passing, the verdict carries
  ` for 2^26, 2^28 only; completeness not checked`, and the record has `completeness_checked: false`.
  The README's step 4 is the full-packet form.

    python code/pa_certify_big.py --cells 26 28 --sep-tag _pa9   # CONSTRUCTION, those two cells
    python code/pa_certify_big.py --verify                       # VERIFICATION of the whole packet
    python code/pa_certify_big.py --verify --cells 26 28         # VERIFICATION of two cells only
"""
from __future__ import annotations

import argparse
import glob
import io
import json
import os
import sys
import time
from fractions import Fraction

# ACF_THREADS first, then NUMBA_NUM_THREADS (what a launcher sets), then 4.
NTHREAD = int(os.environ.get("ACF_THREADS") or os.environ.get("NUMBA_NUM_THREADS") or "4")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, str(NTHREAD))

import numpy as np                                      # noqa: E402
import numba                                            # noqa: E402
from numba import set_num_threads                       # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pa_certify as PC                                 # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = PC.RES
OUT = os.path.join(RES, "pa_certify_big.txt")
OUTJ = os.path.join(RES, "pa_certify_big.json")
P2_WIDTH = Fraction(1, 10 ** 6)

_lines: list[str] = []
REC: list[dict] = []
FAILS: list[str] = []
NOTRUN: list[str] = []


def say(s: str = "") -> None:
    print(s)
    _lines.append(s)
    with io.open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(_lines) + "\n")


# The split pa_certify.py's header describes, with this file's clause numbering.  CERTIFICATION: C0-C9
# (the cell, the integer arithmetic, the stored fields, the inventory), the V<e> re-verifications, VOID,
# and P1.  DIAGNOSTIC: everything below, each a comparison with something outside the certificate.  The
# exit code reports certificate validity alone; no clause's test or threshold is changed by the split.
DIAGNOSTIC: dict[str, tuple[bool, str]] = {
    "P2": (True, "a width TARGET of 10^-6, not a condition on the certificate: the width is whatever"
                 " the colgen primal and the repaired dual give, and at 2^28 it is 1.2e-5 -- section 11"
                 " claims only 10^-4 there"),
    "P3": (False, "it compares this certificate with the published eq:certint integers; a failure is a"
                  " finding about those integers, not about this certificate"),
    "P4": (False, "it compares the exact interval with the separator's stored floats, one side exact"
                  " and the other a float the solver printed"),
    "P5": (False, "a file-size budget for the shipped certificate files"),
}


def clause(name: str, ok: bool, detail: str) -> None:
    if not ok:
        FAILS.append(name)
    REC.append(dict(name=name, holds=bool(ok), detail=detail))
    say("%-4s %-8s %s" % (name, "HOLDS" if ok else "FAILS", detail))


def clause_not_run(name: str, detail: str) -> None:
    """A clause this invocation does not establish: not a pass, not a failure (see pa_certify.py)."""
    NOTRUN.append(name)
    REC.append(dict(name=name, holds=None, state="not run", detail=detail))
    say("%-4s %-8s %s" % (name, "NOT RUN", detail))


def say_fails() -> None:
    say("FAILS: %s" % (", ".join(FAILS) if FAILS else "none"))
    if NOTRUN:
        say("NOT RUN: %s" % ", ".join(NOTRUN))


def verdict(extra=None, emit=None, scope: str = "") -> int:
    """the two verdict lines over this file's clause record; the exit code is validity alone."""
    return PC.verdict_lines(FAILS, DIAGNOSTIC, say if emit is None else emit, extra, scope,
                            NOTRUN)


# pa_certify writes its clause lines into its own buffer; route them here so C1/C2 appear in this
# file's record too.
PC.say = say
PC.clause = clause
PC.clause_not_run = clause_not_run


def find_sep(e: int, tag: str | None) -> str | None:
    """the separator to read: an explicit tag if given, else the newest that carries a primal."""
    if tag:
        p = os.path.join(RES, "conedual_colgen_sep_2e%d%s.npz" % (e, tag))
        return p if os.path.exists(p) else None
    cands = sorted(glob.glob(os.path.join(RES, "conedual_colgen_sep_2e%d*.npz" % e)))
    withp = []
    for p in cands:
        try:
            with np.load(p, allow_pickle=False) as z:
                if "y_support" in z.files and np.asarray(z["y_support"]).size:
                    withp.append(p)
        except Exception:
            pass
    return max(withp, key=os.path.getmtime) if withp else None


def make_float_optimum(sep_path: str, dual_path: str | None = None):
    """a `pa_certify.float_optimum` reading the primal from `sep_path` and the dual from
    `dual_path` (the same separator unless a different one is named)."""
    dpath = dual_path or sep_path

    def float_optimum(C):
        e, rows = C["e"], C["rows"]
        cR = C["cR"]
        R = int(cR[0])
        vB = np.flatnonzero(C["code"] == 1)
        z = np.load(sep_path, allow_pickle=False)
        zd = z if dpath == sep_path else np.load(dpath, allow_pickle=False)
        zr = np.asarray(z["rows"], dtype=np.int64)
        w = np.asarray(zd["w"], dtype=np.float64)
        zt = np.asarray(z["t"], dtype=np.float64)
        has = "y_support" in z.files and "y_values" in z.files
        sup = np.asarray(z["y_support"], dtype=np.int64) if has else np.empty(0, dtype=np.int64)
        yv = np.asarray(z["y_values"], dtype=np.float64) if has else np.empty(0, dtype=np.float64)
        rows_ok = (np.array_equal(zr, rows)
                   and np.array_equal(np.asarray(zd["rows"], dtype=np.int64), rows))
        t_ok = np.array_equal(zt, cR / R)
        dotw = float(np.dot(cR.astype(np.float64), w))
        idx_ok = bool(sup.size and sup.min() >= 0 and sup.max() < vB.size and sup.size == yv.size)
        cols = vB[sup] if idx_ok else np.empty(0, dtype=np.int64)
        col_ok = bool(idx_ok and np.all(C["code"][cols] == 1))
        y_ok = bool(idx_ok and np.all(yv >= 0.0))
        ok = rows_ok and t_ok and dotw > 0.0 and has and idx_ok and col_ok and y_ok
        clause("C6", ok,
               "2^%d primal from %s, dual from %s: rows equal %s, t = c_vR/R exactly %s,"
               " <c_vR,w> = %.4e > 0 %s, carries a primal %s (%d columns, in range %s,"
               " even-omega %s, y >= 0 %s)"
               % (e, os.path.basename(sep_path), os.path.basename(dpath), rows_ok, t_ok, dotw,
                  dotw > 0.0, has, sup.size, idx_ok, col_ok, y_ok))
        # `r` is used by two_sided ONLY as the dual direction, `yv` ONLY as the primal (see header)
        return (cols, yv * float(R), w,
                "primal %s (%d columns); dual %s"
                % (os.path.basename(sep_path), sup.size, os.path.basename(dpath)), ok)

    return float_optimum


def store_if_better(c):
    """write the certificate only if none is stored or the new enclosure is strictly inside it.

    The comparison is on the EXACT squared endpoints, the rationals the file stores, and the rule is
    containment: the new `[lower^2, upper^2]` must lie inside the stored one and be strictly smaller
    on at least one side.  So a replacement can never lose any part of the stored enclosure, and two
    intervals that merely cross -- a better lower endpoint bought with a worse upper one -- leave the
    stored certificate in place.

    Comparing twelve-decimal widths instead, as this once did, is weaker in two ways: two intervals
    with the same rounded width can enclose different things, and an interval that is narrower after
    rounding can still extend past the stored one on one side.  The guard exists because an earlier
    version overwrote the 2^18 certificate (width 1e-12) with its own 2.6e-4 one; the clause it
    enforces (C7) is stated in exact rationals and now the code is too.  Returns (bytes, kept) with
    kept False when the stored certificate wins.

    A stored certificate that cannot be read is also KEPT, not overwritten.  An earlier version
    swallowed the exception and fell through to the write, so a damaged packet file would have been
    replaced silently by whatever this run happened to produce -- the wrong direction for a packet,
    and the same mistake in kind as letting a malformed field pass.
    """
    p = PC.cert_path(c["e"])
    if os.path.exists(p):
        try:
            with np.load(p, allow_pickle=False) as z:
                olo = Fraction(int(str(z["lower2_num"])), int(str(z["lower2_den"])))
                ohi = Fraction(int(str(z["upper2_num"])), int(str(z["upper2_den"])))
            nlo, nhi = c["lo2"], c["hi2"]
            inside = (nlo >= olo and nhi <= ohi)
            strict = inside and (nlo > olo or nhi < ohi)
            if not strict:
                why = ("not inside it" if not inside else "identical to it")
                say("           kept the stored %s: this run's exact enclosure is %s"
                    " (stored [%s, %s], this run [%s, %s], as squares)"
                    % (os.path.basename(p), why, olo, ohi, nlo, nhi))
                return os.path.getsize(p), False
            say("           replacing the stored %s: this run's exact enclosure lies strictly inside"
                " it (stored [%s, %s], this run [%s, %s], as squares)"
                % (os.path.basename(p), olo, ohi, nlo, nhi))
        except Exception as ex:
            say("           KEPT the stored %s UNREAD: it exists but its stored interval cannot be"
                " read (%s: %s), so this run does not overwrite it -- inspect it by hand"
                % (os.path.basename(p), type(ex).__name__, ex))
            return os.path.getsize(p), False
    return PC.store(c), True


def main(argv) -> int:
    ap = argparse.ArgumentParser()
    # default: the two cells this file exists for when certifying, every stored cell when verifying
    ap.add_argument("--cells", type=int, nargs="*", default=None)
    ap.add_argument("--sep-tag", default=None, help="separator tag for the PRIMAL, e.g. _pa9; default: newest with a primal")
    ap.add_argument("--dual-tag", default=None, help="separator tag for the DUAL; default: same as --sep-tag")
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args(argv)
    # In verify mode the inventory is the ADVERTISED set, not whatever files happen to be present:
    # discovering files and succeeding on them lets an empty directory verify and lets a removed
    # certificate shrink the claim silently.  Clause C8 below reports any difference.
    cells = a.cells or (list(PC.ADVERTISED) if a.verify else [26, 28])
    nth = max(1, min(NTHREAD, 8, int(numba.config.NUMBA_NUM_THREADS)))
    set_num_threads(nth)
    T0 = time.time()
    say("pa_certify_big -- eq:cert from a colgen separator that carries its primal")
    say("threads %d (ACF_THREADS or 4, capped at 8 and at numba's own count); the int64 kernels are"
        " guarded by the ACTUAL integer sums of the vectors passed to them (clause C2)" % nth)
    say()
    # The wrapper runs the SAME C9 as pa_certify.py, through pa_certify's own function: the `clause`
    # installed above routes it into this file's record.  A wrapper that skipped one of the verifier's
    # checks would make "the same certificates through the wrapper" false of that check.
    PC.check_c9()
    say()

    if a.verify:
        # A selected-cell run is not a packet check; see pa_certify.py's C6.  The announcement, the
        # verdict and C8 are all scoped to what this invocation actually opens.
        full = set(cells) == set(PC.ADVERTISED)
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
                   ", ".join("2^%d" % e for e in PC.ADVERTISED)))
        found = PC.stored_cells()
        missing = [e for e in PC.ADVERTISED if e not in found]
        extra = [e for e in found if e not in PC.ADVERTISED]
        if full:
            clause("C8", not missing and not extra and bool(found),
                   "inventory is exactly the %d advertised cells: present %s%s%s"
                   % (len(PC.ADVERTISED), ", ".join("2^%d" % e for e in found) or "NONE",
                      "" if not missing else "; MISSING " + ", ".join("2^%d" % e for e in missing),
                      "" if not extra else "; UNADVERTISED " + ", ".join("2^%d" % e for e in extra)))
        else:
            clause_not_run("C8", "the inventory is NOT checked: --cells selected %s, and a run over a"
                                 " subset of the advertised cells %s cannot speak for the packet's"
                                 " contents.  Run `--verify` with no --cells for the full-packet check."
                           % (", ".join("2^%d" % e for e in cells),
                              ", ".join("2^%d" % e for e in PC.ADVERTISED)))
        say()
        vs = []
        for e in cells:
            t0 = time.time()
            v = PC.verify_one(e)
            if v:
                v["verify_s"] = round(time.time() - t0, 1)
                say("           (2^%d re-verified in %.1f s, band rebuild included)"
                    % (e, v["verify_s"]))
            vs.append(v)
        say()
        say("verify wall %.1f s" % (time.time() - T0))
        # PC.verify_one records C0/C2/C4/C5/V<e> through the `clause` this file installed above, so
        # FAILS and REC here already hold them.
        say_fails()
        d = dict(mode="verify" if full else "verify-selected", cells=cells,
                 advertised=list(PC.ADVERTISED), completeness_checked=full,
                 verified=[v for v in vs if v], fails=FAILS, clauses=REC)
        rc = verdict(d, emit=PC.QUIET, scope=scope)
        json.dump(d, io.open(OUTJ, "w", encoding="utf-8", newline="\n"), indent=1)
        verdict(scope=scope)
        return rc

    # CONSTRUCTION MODE: certificates re-derived for the cells given and written; no completeness clause
    # runs, so the verdict is scoped to those cells.
    con_scope = (" for the cells built (%s); completeness not checked"
                 % ", ".join("2^%d" % e for e in cells))
    say("CONSTRUCTION MODE -- certificates re-derived and written for %s; this mode does not check the"
        " packet's contents (that is `--verify` with no --cells)" % ", ".join("2^%d" % e for e in cells))
    say()
    certs, voids, sizes, seps = [], [], {}, {}
    for e in cells:
        sp = find_sep(e, a.sep_tag)
        if sp is None:
            clause("C6", False, "2^%d: no separator carrying a primal (y_support) was found%s"
                   % (e, (" for tag %s" % a.sep_tag) if a.sep_tag else ""))
            voids.append(dict(e=e, void="C6"))
            say()
            continue
        seps[e] = os.path.basename(sp)
        dp = None
        if a.dual_tag:
            dp = os.path.join(RES, "conedual_colgen_sep_2e%d%s.npz" % (e, a.dual_tag))
            if not os.path.exists(dp):
                clause("C6", False, "2^%d: no separator for the dual tag %s" % (e, a.dual_tag))
                voids.append(dict(e=e, void="C6"))
                say()
                continue
        t0 = time.time()
        C = PC.cell(e)
        PC.float_optimum = make_float_optimum(sp, dp)
        c = PC.two_sided(C)
        del C
        if c.get("void"):
            voids.append(c)
            say()
            continue
        c["wall_s"] = round(time.time() - t0, 1)
        c["sep"] = seps[e]
        c["dual_sep"] = os.path.basename(dp) if dp else seps[e]
        certs.append(c)
        sizes[e], kept = store_if_better(c)
        if kept:
            say("           stored %s (%d bytes) -- cell wall %.0f s"
                % (os.path.basename(PC.cert_path(e)), sizes[e], c["wall_s"]))
        say()

    say("clauses")
    if voids:
        clause("VOID", False, "voided cells (registered consequence, enforced in code): "
               + ", ".join("2^%d by %s" % (v["e"], v["void"]) for v in voids))
    clause("P1", bool(certs) and all(c["feasible"] for c in certs),
           "two-sided certificate valid at %s: theta >= 0 on every even column of the whole band,"
           " <c_vR,lambda> < 0, lower <= upper"
           % ", ".join("2^%d" % c["e"] for c in certs))
    clause("P2", bool(certs) and all(c["width"] < P2_WIDTH for c in certs),
           "interval narrower than 1e-6: "
           + ", ".join("2^%d %.2e" % (c["e"], float(c["width"])) for c in certs))

    p3, p3d = True, []
    for c in certs:
        js = sorted(glob.glob(os.path.join(RES, "conedual_colgen_2e%d*.json" % c["e"])))
        c6 = None
        for j in js:
            try:
                c6 = json.load(io.open(j, encoding="utf-8")).get("c6") or c6
            except Exception:
                pass
        if c6 is None:
            continue
        ge = c["lo"] >= Fraction(int(c6), 10 ** 6)
        p3 = p3 and ge
        p3d.append("2^%d lower %s >= %d/10^6 %s" % (c["e"], PC.fstr(c["lo"], 12), c6,
                                                    "yes" if ge else "NO"))
    clause("P3", p3, "at least as strong as eq:certint: " + "; ".join(p3d))

    p4, p4d = True, []
    for c in certs:
        z = np.load(os.path.join(RES, c["sep"]), allow_pickle=False)
        zd = np.load(os.path.join(RES, c["dual_sep"]), allow_pickle=False)
        fl_lo = Fraction(repr(float(zd["bound"]))) if "bound" in zd.files else None
        fl_ub = Fraction(repr(float(z["ub"]))) if "ub" in z.files else None
        okl = (fl_lo is None) or (c["hi"] >= fl_lo)
        oku = (fl_ub is None) or (c["lo"] <= fl_ub)
        p4 = p4 and okl and oku
        p4d.append("2^%d stored float [%s, %s] consistent with the exact interval: %s"
                   % (c["e"], "none" if fl_lo is None else ("%.9f" % float(fl_lo)),
                      "none" if fl_ub is None else ("%.9f" % float(fl_ub)),
                      "yes" if (okl and oku) else "NO"))
    clause("P4", p4, "the exact interval and the separator's stored floats agree: " + "; ".join(p4d))
    clause("P5", all(v < 5_000_000 for v in sizes.values()),
           "certificate files under 5 MB: "
           + ", ".join("2^%d %d B" % (e, v) for e, v in sorted(sizes.items())))

    say()
    say("certified intervals, 12 decimals")
    for c in certs:
        say("   2^%-3d  %s  <=  kappa  <=  %s   (width %.2e, %s)"
            % (c["e"], PC.fstr(c["lo"], 12), PC.fstr(c["hi"], 12), float(c["width"]), c["sep"]))
    say()
    say("wall %.1f s" % (time.time() - T0))
    say_fails()

    d = dict(mode="certify", threads=nth, wall_s=time.time() - T0,
             certificates=[dict(e=c["e"], lower=PC.fstr(c["lo"], 12),
                                upper=PC.fstr(c["hi"], 12), lower2=str(c["lo2"]),
                                upper2=str(c["hi2"]), repair=c["repair"],
                                theta_min=c["theta_min"], support=int(c["cols"].size),
                                K=c["K"], D=c["D"], sep=c["sep"],
                                dual_sep=c["dual_sep"], src=c["src"],
                                wall_s=c["wall_s"], bytes=sizes.get(c["e"]))
                           for c in certs],
             cells=cells, completeness_checked=False, clauses=REC, fails=FAILS)
    rc = verdict(d, emit=PC.QUIET, scope=con_scope)
    json.dump(d, io.open(OUTJ, "w", encoding="utf-8", newline="\n"), indent=1)
    say("wrote %s and %s" % (os.path.basename(OUT), os.path.basename(OUTJ)))
    verdict(scope=con_scope)
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

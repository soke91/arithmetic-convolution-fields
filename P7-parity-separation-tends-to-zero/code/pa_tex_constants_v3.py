# -*- coding: utf-8 -*-
r"""pa_tex_constants_v3 -- every numeric literal in the current manuscript.

WHAT THIS IS FOR

`pa_tex_constants_v2.py` checks 52 constants of the analytic sections and all 52 hold; it does not
reach section 11, whose literals are the certificates' and the records'.  v2 is kept and re-run here
by W1, and this file covers what it does not; the two differ in kind:

  v2   compares a computed value against a literal HARDCODED IN v2.  That cannot notice when the
       manuscript stops printing the literal: v2's V15a/V15b still "hold" against `3767`/`10868`,
       which the current source writes as `3766.52\ldots`/`10867.87\ldots`, and its M1/M2 still check
       float enclosures (`0.1055426`, `0.0763313`, `0.0346234`, `0.0245279`) that are no longer
       anywhere in the source.
  v3   reads the .tex, TOKENISES every numeric literal in its body, and requires each occurrence to
       be accounted for -- by a clause that recomputes it, by a declared quotation from a cited
       source, or by being inside a `\cite[...]` argument or a version string.  Clause Z1 fails if
       one literal is left over, and equally if this file claims to cover a literal the manuscript no
       longer prints.  The checker is therefore tied to the source, not to a previous reading of it.

THE SEVEN CERTIFIED INTERVALS ARE CHECKED AGAINST THE CERTIFICATES, NOT AGAINST A TEXT FILE

The §11 table is the heart of the evidence, so K1/K2/K3 do not read any `.txt`/`.json` record of a
previous run.  For each cell they open `results/pa_cert_2e<e>.npz`, take the two stored integral
vectors, and recompute BOTH rational bounds of `eq:cert` from scratch,

    lower^2 = J^2/(H D_p)  with  J = -<c_vR,lambda>, H = ||lambda||^2, D_p = ||P||^2,
    upper^2 = ||D c_vR - A q||^2 / (D^2 ||P||^2),

the row sums `A q` by `pa_certify.k_fit` and every big integer in object dtype, then round the pair
outward and compare the decimals with the table.  `c_vB`, `c_vR` and the row list come from
`cell_cache` (`CACHE USED` on every load).

WHAT THIS FILE DOES NOT RE-ESTABLISH, AND WHERE IT IS ESTABLISHED

Dual feasibility -- `theta_lambda(l) >= 0` on every `l in vB` -- is not rechecked here.  It needs the
whole band swept, which `pa_certify.py --verify` does on all seven certificates; clause A2 records
that the vectors this file reads are the vectors that run verifies (same file, same row list).
Duplicating a 40 s band rebuild per cell inside a constants checker would buy nothing.  Likewise the
band rebuild control (`pa_certify`'s C1) stays there; here `cell_cache` is trusted, with its own
fingerprint.

CHANGELOG, as far as a reader needs it

  * The line anchors below are re-pinned whenever the manuscript moves text, and are measured against
    the source rather than taken from any summary of what moved: insertions earlier in section 11 have
    shifted the interval table, the column counts and the whole `e = 30` block in the past, and an
    anchor taken on trust would have made eight clauses fail on position alone.
  * `TEXP` is resolved by `find_tex`, which looks for the manuscript both where the working tree keeps
    it and beside `code/` where the reproduction packet keeps it.  It was once hardcoded to the first,
    so the checker could not run from inside the packet it ships in -- the only layout a reader has.
    The packet copy is byte-identical, so every line anchor holds in both.
  * R1 reads the authors' recorded verifier times from `results/pa_certify_timings_frozen.json`, which
    no reproduction step writes and which carries the hardware with the figures.  R1 is a provenance
    claim -- that section 11's figures are the ones the authors' run recorded, on the machine section 11
    names -- and step 1 of the packet's reproduction overwrites `pa_certify.json` with the reader's own
    hardware times.  Reading that mutable file made the clause fail for a reason that is not a defect
    in the paper.  A reader's own record is reported beside the authors' and is not scored.
  * R1 also requires the three per-cell figures to sum to the printed total, which is what makes them
    a breakdown of one run rather than three independent walls, and requires the superseded figures to
    be absent from the source.
  * Z1 counts OCCURRENCES, not `(literal, line)` keys.  Two literals are printed twice on one line, so
    key counts did not sum to the occurrence total and the split could not be checked by a reader.
    Z1's verdict was never affected -- nothing was ever left over -- but its arithmetic is now
    reportable: it prints the occurrence split, requires it to sum, and prints the key counts and the
    duplicates behind them.

A1 RECORDS the manuscript's SHA-256; it does not pin it, and fails only if the manuscript is missing.
What ties this checker to a particular revision is the line anchors plus Z1's two-way coverage, which
is the stronger tie: a moved or renamed literal fails a clause by name, where a changed hash would
only say that something, somewhere, had moved.

CLAUSES.  No literal below is edited after any output of the revision that introduced
it.  Any failure sets a non-zero exit code.

  A1  CONTROL.  The .tex read is the manuscript -- the `.tex` beside `code/` in a built packet, or
      `paper/kappa_tends_to_zero.tex` one level above `code/` in the authors' tree, whichever this run
      finds first; the record names the file it read.  Its SHA-256, byte count and line count are
      recorded, and every literal checked below is located in it by line number.  CAN FAIL (file
      missing or a literal not at the line claimed).
  A2  CONTROL.  At each of the seven cells the certificate's `rows` equals the cell's `rows` element
      for element, and its `e`/`N` fields agree with the cell; so the vectors used below are the ones
      `pa_certify --verify` checked for feasibility.  CAN FAIL; voids the cell.
  K1  The seven twelve-decimal endpoints of the §11 table equal the outward twelve-decimal rounding
      of the two rational bounds recomputed from the stored certificates.
  K2  The seven seven-decimal pairs in the table's last column equal the twelve-decimal endpoints
      rounded outward at seven, and equal the outward seven-decimal rounding of the rationals
      themselves (the two routes must agree).
  K3  The seven widths are `1.0,1.0,2.0,3.0,8.0` in units of `10^-12` at `e <= 24`, `5.6*10^-7` at
      `2^26` and `1.2*10^-5` at `2^28`, each the exact width rounded to one decimal; and the two
      rounded ones round UP, so they are valid upper bounds -- which is what the source now asserts
      when it introduces them with "at most".  The `2^28` figure is printed twice: in section 11 and
      again in the appendix's diagnostic sentence, which cross-references section 11 for it; both
      occurrences are required and accounted for here.
  K4  The precision claim: the width is in `(10^-7,10^-6)` at `2^26` and in `(10^-5,10^-4)` at
      `2^28`, which is what makes the two rounded widths six- and four-decimal figures.  Both ends of
      each bracket are checked, so the powers named are sharp.  These four powers are this clause's
      own, not quotations: the only width bound the source states at `2^28` outside the table is the
      `1.2*10^-5` of section 11, which K3 recomputes and K7 reads as a relation.
  K5  The complete even-column counts `4\,414`, `18\,685`, `76\,946` at `e = 16,18,20` equal
      `cell_cache`'s `nB`.
  K6  Each table row carries its OWN exponent label, parsed from the row and required to be the
      cell whose certificate produced that row's endpoints -- without it, relabelling a row
      (`16` -> `18`) would pass while asserting one cell's interval for another; and the exponent
      lists and `Y` values section 11 names outside the table are present where claimed.  CAN FAIL.
  K7  THE RELATION SYMBOLS, their DIRECTION and their TRUTH.  K1/K2/K6 check which digits are
      printed on which line; they do not read the relation between them, so reversing a table row's
      `<` to `>` left every clause holding while the row asserted the opposite of the certificate.
      This clause parses the relation of every inequality the checker anchors -- the table's two
      bound columns and their header, each row's two-sided seven-decimal statement, the printed
      `kappa(2^30) > 0.017787`, the two "at most" widths, the appendix's "below" statement of the
      target and of the width section 11 states instead -- whose cross-reference is resolved, so
      "stated in \S\ref{...}" must name the section that does state it -- and the
      construction's four `kappa <= ...` -- and from the symbol TOGETHER WITH the side the literal
      sits on decides what the text asserts: a LOWER bound on the quantity, or an UPPER one.  Each
      assertion is then tested against the certified rational enclosure of that quantity: a printed
      lower bound must be `<=` the certified lower bound and a printed upper bound `>=` the certified
      upper bound, strictly where the printed relation is strict.  Reversing any one of them makes
      the printed claim contradict the certificate (or claim more than it establishes) and fails
      here.  The comparisons are on exact rationals: squared for the `kappa` statements, unsquared
      for the widths.  K7 accounts for no literal -- that is K1/K2/K3/G3/G4's and v2's work; it
      checks what the source says ABOUT them.  CAN FAIL.
  G1  At `2^30`, `|vB| = 86\,592\,333` and `|vR| = 89\,520\,933` equal `cell_cache`'s `nB`, `nR`.
  G2  The four integers `J`, `H`, `D_p`, `c_6` printed at `2^30` equal those stored in
      `results/conedual_colgen_2e30_n30b.json`, `J > 0`, and `D_p` equals `||c_vB - c_vR||^2`
      recomputed from the `2^30` cell in exact arithmetic.
  G3  `eq:certint` holds as printed -- `c_6^2 H D_p < J^2 * 10^12` -- it is equivalent to
      `lower^2 > (c_6/10^6)^2`, `c_6` is MAXIMAL (`c_6 + 1` fails it), and `c_6/10^6 = 0.017787` is
      the printed lower bound.
  G4  The uncertified solver value `0.01779` at `2^30` is the stored `ub` rounded to five decimals,
      and `\approx` occurs exactly once in section 11 -- on that value and nowhere else, which is what
      the section claims when it states everything else as certified.
  G5  HOW THE STORED `2^30` DUAL IS BUILT, recomputed from the two files the packet ships for it --
      `results/pa_cert_2e30_lower.npz` and the separator `results/conedual_colgen_sep_2e30_n30b.npz`.
      Section 11 states the construction, and every numeral that statement adds is checked here rather
      than taken on trust: the scale `10^12`, the sign, the row `d = 1`, the bound `2^62`.
      (a) the cell's, the certificate's and the separator's row lists are equal element for element,
          and the first row is `d = 1`, so what follows is a coordinate-by-coordinate comparison;
      (b) `lambda_d = round(-10^12 w_d)` at EVERY row but `d = 1` -- the fixed scale and the sign the
          text states, with nearest-integer rounding, which is the builder's `np.rint`;
      (c) at `d = 1` the rounded value is NOT what is stored: `lambda_1` is positive, differs from
          `round(-10^12 w_1)`, and is therefore the repair amount and nothing else;
      (d) `<c_vR, w> > 0` in float64 (its value is far from zero, so the sign is not in doubt) while
          `<c_vR, lambda> < 0` in exact integer arithmetic and equals `-J` for the stored `J`: the
          minus sign in (b) is what makes `J > 0`;
      (e) `sum_d |lambda_d| <= 2^62`, the accumulator bound AT THAT SCALE -- the text says it is
          checked rather than assumed, and this is that check on the stored vector;
      (f) `max_d |lambda_d|` and `max_d |w_d|` are attained at the SAME row, away from `d = 1`, and
          the first is the second times `10^12` rounded: they differ by the scale and by nothing else,
          which is the relation the text tells a reader to expect between the two files.
      Four printed statements are also required to be on the lines this clause names.  CAN FAIL.
      It is the only clause about how that vector was CONSTRUCTED; that it certifies the endpoint is
      `L0`-`L4` of `pa_certify.py --verify-lower30`, not this.
  Q1  Constants quoted from a cited source (`0.2795`, `6.455`, `1.25506`) are checked for PRESENCE at
      the line claimed and for sitting in the attribution recorded here.  Whether the source really
      says them is a citation check, not arithmetic: this clause records the attribution each sits
      in rather than re-deriving the constant.
  M1  DECLARED CATEGORY, the classification line.  The MSC codes and keywords printed after the
      abstract are a classification, not a measurement: `2020` is the scheme's edition and `11N35`
      names a subject area.  An MSC code is two digits, a letter and two digits, so no part of it is a
      numeric literal in this file's sense and the inventory counts nothing from those lines -- which
      this clause turns from a silence into a checked statement: the codes are exactly `11N35`
      (primary) and `11N36`, `11N25`, `90C05` (secondary) in that split, the year is `2020`, every code
      has MSC shape, the keywords carry no digit, and the block's contribution to the literal
      inventory is zero, read off the inventory itself.  Nothing is exempted: a numeral placed in this
      block would be counted and would need a clause like any other.  CAN FAIL.
  X1  The `beta = 6` crossover `L = 3766.52\ldots` is the root of
      `c_1/(L sqrt(log L)) + c_2 L^2 exp(-L/(24 log L)) = 1`.
  X2  The companion `log_2 N = 10867.87\ldots` is `2L*/log 2` at that root, rendered as the source's
      `\ldots` declares (truncated), and the ROUNDED value `10867.88` must appear nowhere in the body:
      under `\ldots` it would assert a digit the quantity does not have.  CAN FAIL.
  X3  The two SUFFICIENT thresholds for the standing condition are the continuous crossings rounded
      UPWARD, which is what a sufficient threshold requires.  With `L6` the root of `6 log L = L` and
      `L12` the root of `12 log L = L`, the thresholds are `2 L6/log 2` (for `Y < Q`) and
      `2 L12/log 2` (for `Y <= sqrt Q`); their two-decimal upward roundings are the printed `49.05`
      and `132.51`, and each printed figure is `>=` the crossing it stands for.  A one-decimal
      rounding of the second, `132.5`, lies BELOW its crossing and is not in the source.  v2's V18c
      tests the same quantity at one decimal and matches `132.5` as a prefix of `132.51`; X3 is the
      clause that tests what the source now prints.  CAN FAIL.
  B1  The Buchstab derivation's `-0.1931\ldots` is `1/2 - log 2`, the value of `1/v - log v` at
      `v = 2`, and `1/v - log v` is `1` at `v = 1` and strictly decreasing on `[1,2]`.
  N1  The constant `1.5` of `log X >= 1.5 L` is `3/2`, from `X >= Q^{3/2}`.
  N2  For the band-density remark of section 11.  From `cell_cache`'s row
      counts at the eight computed cells `2^16 ... 2^30`: `R/N` rounded to four decimals is `0.0768`
      at `2^16` and `0.0834` at `2^30`, `N/R` rounded to two is `13.03` and `11.99`, and the limit
      `pi^2/(3(1-log 2))` truncated to four is `10.7213`; and the remark's two quantified claims hold
      over all eight -- `R/N` strictly increasing and below `0.09` throughout, `N/R` strictly
      decreasing and above the limit throughout.
  R1  The four verifier runtimes printed in §11's verifier paragraph -- the wall for all seven cells,
      then the five cells `e <= 24` together, `2^26` and `2^28` -- are the times in
      `results/pa_certify_timings_frozen.json`, the authors' record that NO reproduction step writes.
      The figures are read out of the source rather than named here, so the clause states a relation and
      keeps holding when they are updated.  It also requires that the three parts are a breakdown of
      that one run and not independent walls -- they sum to the total under the one-decimal rounding both
      print, compared as EXACT RATIONALS on those printed strings against the rational tolerance `2/10`,
      so the boundary is the boundary and no verdict here depends on a binary float -- that
      the paragraph names the machine, and that the three superseded figures (`16.3`, `7.6`, `31.4`) are
      gone from the source.  It also requires the frozen record's `verifier_sha256` to be the SHA-256 of
      the `pa_certify.py` this run imported: a timing record that names other bytes than the shipped
      verifier does not bind its measurements to what a reader runs, which is the one way this clause
      could hold while meaning nothing.  A reader's own run of step 1 is reported beside it and NOT
      scored.  CAN FAIL.
  W1  Every one of v2's 52 clauses still holds when re-run against the current tree (v2 is imported
      and called, not copied, and its outputs are redirected so its stored record is untouched).
  W2  Every literal v2 checks that is ABSENT from the current .tex is one the section 11 rewrite records as
      removed, and has a replacement clause here; no v2 literal has vanished unaccounted for.
  Z1  COVERAGE.  Every numeric literal OCCURRENCE in the .tex body -- decimals, integers of five or
      more digits, and `\,`-grouped integers, outside the bibliography -- is accounted for exactly
      once: by a clause above, by v2, by a declared quotation from a cited source, or as a reference
      inside `\cite[...]` or a version string.  Nothing left over, nothing claimed that is not there,
      and the reported split SUMS to the total (occurrences, not `(literal, line)` keys -- a defect
      of a final consistency check, where two literals printed twice on one line made the published figures
      read `58 + 45 = 103` against a total of `105`).
  Y1  DECLARED CATEGORY, the numerals a reference entry carries.  A bibliography entry prints a
      volume, a year and a page range -- locators, not quantities -- and the checked body ends at
      `\begin{thebibliography}`, so the inventory never sees them and `Z1` is not asked to cover them.
      This clause checks that for the locators the entries print (`33`, `1976`, `565--576`): each is in
      the bibliography, none occurs as a standalone number anywhere in the body, and none is a tracked
      literal of the body.  A locator that migrates into the body fails here and must then be covered
      like any other literal.  CAN FAIL.
  D1  NOTATION.  Where truncation and rounding of a computed value differ, the printed literal
      matches the rendering the source declares (`\ldots` -> truncated, `\approx`/`near`/`>=` ->
      rounded).  CAN FAIL.

  NOT REGISTERED: whether a quoted classical constant is correctly quoted (that is a citation check,
  not arithmetic -- Q1 records the citation each constant sits in);
  and dual feasibility, which `pa_certify --verify` establishes.

    python code/pa_tex_constants_v3.py
"""
from __future__ import annotations

import glob
import hashlib
import io
import json
import os
import re
import sys
import tempfile
from fractions import Fraction

NTHREAD = int(os.environ.get("ACF_THREADS", "4"))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, str(NTHREAD))

import mpmath as mp                                      # noqa: E402
import numpy as np                                       # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pa_tex_constants_v2 as TC2                        # noqa: E402  (helpers + its 52 clauses)
import pa_certify as PC                                  # noqa: E402  (pure exact helpers only)
from cell_cache import load as load_cell                 # noqa: E402

from pa_tex_constants_v2 import rnd, rndup, trunc        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), "results")


def find_tex() -> str:
    """the manuscript: `paper/kappa_tends_to_zero.tex` one level above `code/`, or the `<paper_id>.tex`
    beside `code/` in a built packet.

    Both layouts must work: the packet a reader has carries the manuscript as a sibling of `code/`,
    `results/` and `lib/`, and the authors' tree keeps it one directory up.  The search is upward and
    by name, as `_sieve_shared` and `cell_cache` already do for the files they need, so no layout is
    hardcoded and the checker runs unchanged in both.
    """
    up = os.path.dirname(HERE)
    repo = os.path.join(up, "paper", "kappa_tends_to_zero.tex")
    for p in [repo] + sorted(glob.glob(os.path.join(up, "*.tex"))):
        if os.path.isfile(p):
            return p
    return repo


TEXP = find_tex()
OUT = os.path.join(RES, "pa_tex_constants_v3.txt")
OUTJ = os.path.join(RES, "pa_tex_constants_v3.json")

mp.mp.dps = 40
CELLS = (16, 18, 20, 22, 24, 26, 28)
# the eight cells the density remark of section 11 speaks of (clause N2)
DENS_CELLS = (16, 18, 20, 22, 24, 26, 28, 30)

_lines: list[str] = []
REC: list[dict] = []
BAD: list[str] = []
COVER: dict[tuple[str, int], str] = {}


def say(s: str = "") -> None:
    print(s)
    _lines.append(s)
    with io.open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(_lines) + "\n")


def clause(name: str, ok: bool, detail: str) -> None:
    if not ok:
        BAD.append(name)
    REC.append(dict(name=name, holds=bool(ok), detail=detail))
    say("%-4s %-8s %s" % (name, "HOLDS" if ok else "FAILS", detail))


def cover(name: str, token: str, *lines: int) -> None:
    """record that clause `name` accounts for `token` at each of `lines` (clause Z1)."""
    for ln in lines:
        COVER[(token, ln)] = name


# =============================== the .tex and its literals ==================================
# a numeric literal: digits with LaTeX \, groups, a decimal, or an integer of five or more digits.
# Shorter bare integers (2, 6, 24, 852, 229, exponents, size labels) are out of scope and declared so.
PAT = re.compile(r"[0-9]+(?:\\,[0-9]{3})+(?:\.[0-9]+)?|[0-9]+\.[0-9]+|[0-9]{5,}")
CITE = re.compile(r"\\cite\[")
VERSION = re.compile(r"v4\.33\.1")


def rel_or_name(path):
    """a short name for a path, without assuming it shares a drive with this script."""
    try:
        return os.path.relpath(path, os.path.dirname(HERE))
    except ValueError:
        return os.path.basename(path)


def sec_line(lines, i, sec, body):
    """is line index `i` inside the paragraph that starts at character offset `sec` of the body?"""
    off = sum(len(x) + 1 for x in lines[:i])
    return off >= sec


def read_tex():
    raw = io.open(TEXP, encoding="utf-8").read()
    body = raw.split(r"\begin{thebibliography}")[0]
    return raw, body, body.splitlines()


def cite_spans(line: str):
    """spans of each `\\cite[...]` optional argument on this line."""
    out = []
    for m in CITE.finditer(line):
        i, depth = m.end(), 1
        while i < len(line) and depth:
            if line[i] == "[":
                depth += 1
            elif line[i] == "]":
                depth -= 1
            i += 1
        out.append((m.start(), i))
    return out


def occurrences(lines):
    """every (token, line, is_reference) in the body."""
    out = []
    for i, ln in enumerate(lines, 1):
        sp = cite_spans(ln)
        for m in PAT.finditer(ln):
            ref = any(a <= m.start() < b for a, b in sp) or bool(
                VERSION.search(ln[max(0, m.start() - 2):m.end() + 2]))
            out.append((m.group(0), i, ref))
    return out


def at(lines, token: str, line: int) -> bool:
    """is `token` literally present on that 1-based line of the body?"""
    return 1 <= line <= len(lines) and token in lines[line - 1]


# =============================== the certificates ==========================================
def recompute(e: int):
    """both rational bounds of eq:cert, recomputed from the stored certificate and the cell."""
    S = load_cell(e)
    z = np.load(os.path.join(RES, "pa_cert_2e%d.npz" % e), allow_pickle=False)
    rows = np.ascontiguousarray(np.asarray(S["rows"], dtype=np.int64))
    cb = np.asarray(S["cb"], dtype=np.int64)
    cr = np.asarray(S["cr"], dtype=np.int64)
    srows = np.asarray(z["rows"], dtype=np.int64)
    same = bool(np.array_equal(rows, srows) and int(z["e"]) == e and int(z["N"]) == 2 ** e)
    lam = np.asarray(z["dual"], dtype=np.int64)
    cols = np.ascontiguousarray(np.asarray(z["primal_cols"], dtype=np.int64))
    q = np.ascontiguousarray(np.asarray(z["primal_num"], dtype=np.int64))
    D = int(z["primal_den"])
    P2 = PC.sq(cb - cr)
    J = -PC.dot(cr, lam)
    H = PC.sq(lam)
    lo2 = Fraction(J * J, P2 * H) if (J > 0 and H > 0) else Fraction(0)
    fit = np.empty(rows.size, dtype=np.int64)
    PC.k_fit(cols, q, rows, fit)               # the guards that let this run are pa_certify's C2
    res = cr.astype(object) * D - fit.astype(object)
    hi2 = Fraction(int(np.sum(res ** 2)), D * D * P2)
    return dict(e=e, ok=same, lo2=lo2, hi2=hi2, nB=int(S["nB"]), nR=int(S["nR"]),
                nrows=int(rows.size), ncols=int(cols.size), D=D, J=J, H=H, P2=P2)


def floor_at(f: Fraction, nd: int) -> Fraction:
    s = 10 ** nd
    return Fraction(f.numerator * s // f.denominator, s)


def ceil_at(f: Fraction, nd: int) -> Fraction:
    s = 10 ** nd
    return Fraction(-((-f.numerator * s) // f.denominator), s)


# =============================== main ======================================================
def main() -> int:
    say("pa_tex_constants_v3 -- every numeric literal in the current manuscript")
    # the manuscript actually read, not a fixed path: in a built packet it is the `.tex` beside `code/`
    say("source: %s" % rel_or_name(TEXP))
    say("clauses A1-A2, K1-K7, G1-G5, X1-X3, B1, N1-N2, R1, W1-W2, Q1, M1, Y1, Z1, D1 in the header")
    say()

    # ---------------- A1 -------------------------------------------------------------------
    if not os.path.exists(TEXP):
        clause("A1", False, "the manuscript is not at %s" % TEXP)
        return 1
    raw, body, lines = read_tex()
    sha = hashlib.sha256(raw.encode("utf-8")).hexdigest()
    occ = occurrences(lines)
    nref = sum(1 for _, _, r in occ if r)
    clause("A1", True, "%s: %d bytes, %d lines (%d in the body before the bibliography), SHA-256 %s;"
                       " %d literal occurrences, %d of them inside \\cite[...] or a version string"
           % (rel_or_name(TEXP), len(raw.encode("utf-8")),
              raw.count("\n") + 1, len(lines), sha, len(occ), nref))
    say()

    # ---------------- the certificates, A2 --------------------------------------------------
    say("the seven certificates, both rational bounds recomputed from the stored vectors")
    C = {}
    for e in CELLS:
        C[e] = recompute(e)
    a2 = all(C[e]["ok"] for e in CELLS)
    clause("A2", a2, "certificate rows == cell rows and e/N agree at %s"
           % ", ".join("2^%d %s" % (e, "yes" if C[e]["ok"] else "NO") for e in CELLS))
    if not a2:
        say("     VOID: a certificate does not match its cell; the interval clauses are not scored")
        return 1

    for e in CELLS:
        c = C[e]
        lo, hi = PC.dec_bracket(c["lo2"], c["hi2"], 12)
        c["lo12"], c["hi12"] = PC.fstr(lo, 12), PC.fstr(hi, 12)
        c["lo"], c["hi"] = lo, hi
        l7a, h7a = PC.dec_bracket(c["lo2"], c["hi2"], 7)
        c["lo7"], c["hi7"] = PC.fstr(l7a, 7), PC.fstr(h7a, 7)
        c["lo7b"] = PC.fstr(floor_at(lo, 7), 7)
        c["hi7b"] = PC.fstr(ceil_at(hi, 7), 7)
        c["width"] = hi - lo
        say("   2^%-3d  %s <= kappa <= %s   width %s   (%d rows, %d primal columns, D = 10^%d)"
            % (e, c["lo12"], c["hi12"], "%.4e" % float(c["width"]), c["nrows"], c["ncols"],
               len(str(c["D"])) - 1))
    say()

    # ---------------- K1, K2: the table -----------------------------------------------------
    # Line anchors in the current source.  EVERY anchor in this file was re-measured by LOCATING its
    # literal in the manuscript, not by applying a reported shift: a restructuring job's own line map is
    # a summary, and a clause that trusts it reports a line it did not test.  The locations below were
    # read off the source at SHA-256 7d8aac53..., each located by its own literal (one of them,
    # G5's sign condition, now occurs twice in the source: the appendix's is the one G5 means).
    TAB = {16: 1243, 18: 1244, 20: 1245, 22: 1246, 24: 1247, 26: 1248, 28: 1249}
    k1ok, k1d, k2ok, k2d = True, [], True, []
    for e in CELLS:
        c, ln = C[e], TAB[e]
        p1 = at(lines, c["lo12"], ln) and at(lines, c["hi12"], ln)
        k1ok = k1ok and p1
        k1d.append("2^%d %s,%s %s" % (e, c["lo12"], c["hi12"], "on L%d" % ln if p1 else "NOT ON L%d"
                                      % ln))
        cover("K1", c["lo12"], ln)
        cover("K1", c["hi12"], ln)
        routes = (c["lo7"] == c["lo7b"] and c["hi7"] == c["hi7b"])
        p2 = at(lines, c["lo7"], ln) and at(lines, c["hi7"], ln)
        k2ok = k2ok and routes and p2
        k2d.append("2^%d [%s,%s] %s, routes %s" % (e, c["lo7"], c["hi7"],
                                                   "on L%d" % ln if p2 else "NOT ON L%d" % ln,
                                                   "agree" if routes else "DIFFER"))
        cover("K2", c["lo7"], ln)
        cover("K2", c["hi7"], ln)
    clause("K1", k1ok, "twelve-decimal endpoints recomputed from pa_cert_2e*.npz and found in the"
                       " table: " + "; ".join(k1d))
    clause("K2", k2ok, "seven-decimal column, both rounding routes agreeing: " + "; ".join(k2d))

    # ---------------- K6: the table's exponent labels, and section 11's small integers ------
    # The table associates an exponent with two endpoints.  Checking only that the endpoints appear on
    # a given line leaves the exponent unchecked, so relabelling a row (16 -> 18) would pass while
    # asserting one cell's interval for another.  Here each row's own printed exponent is parsed and
    # required to be the cell whose certificate produced that row's endpoints.
    k6ok, k6d = True, []
    for e in CELLS:
        ln = lines[TAB[e] - 1]
        m = re.match(r"\s*[$](\d+)[$]\s*&", ln)
        got = int(m.group(1)) if m else None
        c = C[e]
        row_ok = (got == e and c["lo12"] in ln and c["hi12"] in ln
                  and c["lo7"] in ln and c["hi7"] in ln)
        k6ok = k6ok and row_ok
        k6d.append("L%d label %s" % (TAB[e], got if got is not None else "UNPARSED"))
    # the exponents section 11 names outside the table, and the construction's Y values
    SMALL = [("e=16,\\dots,24", 1252), ("e=26", 1253), ("e=28", 1253),
             ("e=16,18,20", 1254), ("e=22,\\dots,28", 1255),
             ("e=16,18,20,22", 1286), ("Y=3,3,5,5", 1286)]
    smd = [(t, ln, at(lines, t, ln)) for t, ln in SMALL]
    sm_ok = all(v for _, _, v in smd)
    clause("K6", k6ok and sm_ok,
           "table rows carry their own exponent and the endpoints of that cell (%s); the exponents and"
           " Y values section 11 names outside the table are present (%s)"
           % (", ".join(k6d), "all %d" % len(smd) if sm_ok else
              "MISSING " + ", ".join("%s@L%d" % (t, ln) for t, ln, v in smd if not v)))

    # ---------------- K3, K4: widths --------------------------------------------------------
    WID = {16: ("1.0", 12), 18: ("1.0", 12), 20: ("2.0", 12), 22: ("3.0", 12), 24: ("8.0", 12),
           26: ("5.6", 7), 28: ("1.2", 5)}
    WLN = {16: 1252, 18: 1252, 20: 1252, 22: 1252, 24: 1252, 26: 1253, 28: 1253}
    # The appendix's diagnostic sentence, which states the internal target width, whether it is met at
    # each of the two largest instances, and -- on the next line -- the width section 11 states at the
    # instance where the target is not met, by a cross-reference to that section.
    DIAGLN = 1370
    # that cross-reference repeats the 2^28 width figure, so the same recomputed literal occurs twice
    REPEAT = {28: DIAGLN + 1}
    k3ok, k3d = True, []
    for e in CELLS:
        lit, k = WID[e]
        scaled = C[e]["width"] * (10 ** k)
        got = rnd(mp.mpf(scaled.numerator) / mp.mpf(scaled.denominator), 1)
        up = Fraction(got) >= scaled              # the printed width must not understate the width
        p = at(lines, lit, WLN[e])
        ok = (got == lit) and up and p
        k3ok = k3ok and ok
        k3d.append("2^%d width*10^%d = %s -> %s|%s %s%s" % (e, k, mp.nstr(mp.mpf(float(scaled)), 8),
                                                            got, lit, "ok" if ok else "NO",
                                                            "" if up else " (UNDERSTATES)"))
        cover("K3", lit, WLN[e])
    for e, ln in sorted(REPEAT.items()):
        lit = WID[e][0]
        ok = at(lines, lit, ln)
        k3ok = k3ok and ok
        k3d.append("2^%d width %s repeated on L%d %s" % (e, lit, ln, "ok" if ok else "NO"))
        cover("K3", lit, ln)
    clause("K3", k3ok, "the seven widths, each rounded to one decimal and not understating: "
           + "; ".join(k3d))

    w26, w28 = C[26]["width"], C[28]["width"]
    k4 = (Fraction(1, 10 ** 7) < w26 < Fraction(1, 10 ** 6)
          and Fraction(1, 10 ** 5) < w28 < Fraction(1, 10 ** 4))
    pre26 = os.path.commonprefix([C[26]["lo12"], C[26]["hi12"]])
    pre28 = os.path.commonprefix([C[28]["lo12"], C[28]["hi12"]])
    nd26 = len(pre26.split(".")[1]) if "." in pre26 else 0
    nd28 = len(pre28.split(".")[1]) if "." in pre28 else 0
    clause("K4", k4, "precision in the sense of the width: 10^-7 < %.4e < 10^-6 at 2^26 and"
                     " 10^-5 < %.4e < 10^-4 at 2^28, so six and four decimals; for the record the"
                     " endpoints share only %d and %d leading decimals"
           % (float(w26), float(w28), nd26, nd28))

    # ---------------- K5: the complete column counts ----------------------------------------
    EVEN = {16: ("4\\,414", 1255), 18: ("18\\,685", 1255), 20: ("76\\,946", 1255)}
    k5ok, k5d = True, []
    for e, (lit, ln) in EVEN.items():
        want = int(lit.replace("\\,", ""))
        ok = (C[e]["nB"] == want) and at(lines, lit, ln)
        k5ok = k5ok and ok
        k5d.append("2^%d nB %d vs %s on L%d %s" % (e, C[e]["nB"], lit, ln, "ok" if ok else "NO"))
        cover("K5", lit, ln)
    clause("K5", k5ok, "the complete even-column counts: " + "; ".join(k5d))
    say()

    # ---------------- the 2^30 block --------------------------------------------------------
    say("the e = 30 block: the band, the four integers, and eq:certint")
    S30 = load_cell(30)
    cb30 = np.asarray(S30["cb"], dtype=np.int64)
    cr30 = np.asarray(S30["cr"], dtype=np.int64)
    Dp30 = PC.sq(cb30 - cr30)
    nB30, nR30 = int(S30["nB"]), int(S30["nR"])
    # the two lines the band sizes are printed on, named once so the clause cannot report one line
    # while testing another
    g1ln = (1262, 1263)
    g1 = (nB30 == 86592333 and nR30 == 89520933
          and at(lines, "86\\,592\\,333", g1ln[0]) and at(lines, "89\\,520\\,933", g1ln[1]))
    cover("G1", "86\\,592\\,333", g1ln[0])
    cover("G1", "89\\,520\\,933", g1ln[1])
    clause("G1", g1, "2^30 band from cell_cache: |vB| = %d vs 86\\,592\\,333 (L%d), |vR| = %d vs"
                     " 89\\,520\\,933 (L%d)" % (nB30, g1ln[0], nR30, g1ln[1]))

    jp = os.path.join(RES, "conedual_colgen_2e30_n30b.json")
    if not os.path.exists(jp):
        clause("G2", False, "the 2^30 run record %s is missing" % os.path.basename(jp))
        clause("G3", False, "not scored: no 2^30 run record")
        clause("G4", False, "not scored: no 2^30 run record")
        J30 = H30 = Dp_s = c6 = None
    else:
        j30 = json.load(io.open(jp, encoding="utf-8"))
        J30, H30 = int(j30["J"]), int(j30["H"])
        Dp_s, c6 = int(j30["D_p"]), int(j30["c6"])
        lits = {"4401612119558212802": 1273, "2015370564893177546709684150": 1273,
                "30382600959846": 1273, "17787": 1274}
        pres = all(at(lines, t, ln) for t, ln in lits.items())
        for t, ln in lits.items():
            cover("G2", t, ln)
        g2 = (J30 == 4401612119558212802 and H30 == 2015370564893177546709684150
              and Dp_s == 30382600959846 and c6 == 17787 and J30 > 0 and Dp30 == Dp_s and pres)
        clause("G2", g2, "the four integers of %s: J = %d > 0, H = %d, D_p = %d, c_6 = %d; D_p"
                         " recomputed from the 2^30 cell = %d %s; all four printed on L%d-%d %s"
               % (os.path.basename(jp), J30, H30, Dp_s, c6, Dp30,
                  "MATCHES" if Dp30 == Dp_s else "DIFFERS", min(lits.values()),
                  max(lits.values()), "yes" if pres else "NO"))

        ineq = c6 * c6 * H30 * Dp_s < J30 * J30 * 10 ** 12
        maxi = not ((c6 + 1) ** 2 * H30 * Dp_s < J30 * J30 * 10 ** 12)
        equiv = (Fraction(J30 * J30, H30 * Dp_s) > Fraction(c6, 10 ** 6) ** 2) == ineq
        g3ln = 1276                      # the line the six-decimal bound is printed on
        pb = at(lines, "0.017787", g3ln) and rnd(mp.mpf(c6) / mp.mpf(10) ** 6, 6) == "0.017787"
        cover("G3", "0.017787", g3ln)
        clause("G3", ineq and maxi and equiv and pb,
               "eq:certint: c_6^2 H D_p < J^2*10^12 %s; c_6+1 = %d fails it %s (so c_6 is maximal);"
               " equivalent to lower^2 > (c_6/10^6)^2 %s; printed bound 0.017787 on L%d %s"
               % ("holds" if ineq else "FAILS", c6 + 1, "yes" if maxi else "NO",
                  "yes" if equiv else "NO", g3ln, "yes" if pb else "NO"))

        ubs = rnd(mp.mpf(repr(float(j30["ub"]))), 5)
        # Section 11 runs from its own \section line to the bibliography, and `\approx` must occur
        # exactly once inside it.  The section is located by its LABEL: the revision renamed every
        # section title and kept every label, so the title is not a handle and the label is.
        s11 = next((i for i, ln in enumerate(lines, 1) if r"\label{sec:numerics}" in ln), None)
        apx = [i for i, ln in enumerate(lines, 1) if s11 and i >= s11 and r"\approx" in ln]
        g4ln = 1282                      # the one approximate value section 11 prints
        g4 = (ubs == "0.01779" and at(lines, "0.01779", g4ln) and len(apx) == 1
              and apx == [g4ln])
        cover("G4", "0.01779", g4ln)
        clause("G4", g4, "the solver value: stored ub = %s -> %s vs 0.01779 on L%d %s; section 11"
                         " starts at L%s and contains \\approx on %s (must be that line alone)"
               % (j30["ub"], ubs, g4ln, "ok" if (ubs == "0.01779") else "NO", s11,
                  ", ".join("L%d" % i for i in apx) or "no line"))

    # ---------------- G5: how the stored 2^30 dual is built ---------------------------------
    # Section 11 describes the construction of this vector, and the description is checkable from the
    # two files the packet ships: the certificate and the float separator it was rounded from.  Every
    # numeral the description adds -- the scale 10^12, the sign, d = 1, the bound 2^62 -- is recomputed
    # here rather than taken on trust.  g5ln names the four lines the four statements are printed on.
    rows30 = np.asarray(S30["rows"], dtype=np.int64)
    cp = os.path.join(RES, "pa_cert_2e30_lower.npz")
    # the separator tag is the verifier's own default for this cell, and is the file the README names
    sp = os.path.join(RES, "conedual_colgen_sep_2e30_n30b.npz")
    g5ln = (1376, 1377, 1381, 1382)
    if not (os.path.exists(cp) and os.path.exists(sp)):
        clause("G5", False, "not scored: %s absent"
               % ", ".join(os.path.basename(p) for p in (cp, sp) if not os.path.exists(p)))
    else:
        zc = np.load(cp, allow_pickle=False)
        zs = np.load(sp, allow_pickle=False)
        lam = np.asarray(zc["dual"], dtype=np.int64)
        w = np.asarray(zs["w"], dtype=np.float64)
        cert_rows = np.asarray(zc["rows"], dtype=np.int64)
        sep_rows = np.asarray(zs["rows"], dtype=np.int64)
        # (a) the three row lists agree, so the comparison below is coordinate by coordinate
        aligned = bool(np.array_equal(rows30, cert_rows) and np.array_equal(rows30, sep_rows))
        d1 = int(rows30[0]) == 1
        # (b) lambda_d = round(-10^12 w_d) at every row but d = 1, nearest-integer as the builder uses
        rr30 = np.rint(-(10.0 ** 12) * w)
        fits = bool(np.all(np.abs(rr30) < 2.0 ** 63))
        ri30 = rr30.astype(np.int64) if fits else None
        scale_ok = bool(fits and aligned and np.array_equal(lam[1:], ri30[1:]))
        nmis = int(lam.size - 1 - (np.count_nonzero(lam[1:] == ri30[1:]) if fits and aligned else 0))
        # (c) at d = 1 the rounded value was discarded: what is stored is the positive repair
        rep, disc = int(lam[0]), (int(ri30[0]) if fits else None)
        repair_ok = bool(rep > 0 and disc is not None and rep != disc)
        # (d) <c_O,w> > 0 while <c_O,lambda> = -J < 0: the sign is what makes J positive
        dot_lam = int((lam.astype(object) * cr30.astype(object)).sum())
        dot_w = float(np.dot(w, cr30.astype(np.float64)))
        sign_ok = bool(dot_w > 0 and dot_lam < 0 and (-dot_lam) == int(str(zc["J"])))
        # (e) the accumulator bound at that scale, as a claim about the stored vector
        asum = int(np.abs(lam).astype(object).sum())
        bound_ok = asum <= 2 ** 62
        # (f) the two maxima differ by the scale and by nothing else: same row, rounded value
        il, iw = int(np.argmax(np.abs(lam))), int(np.argmax(np.abs(w)))
        mx_ok = bool(int(rows30[il]) == int(rows30[iw]) and il > 0
                     and abs(int(lam[il])) == int(np.rint((10.0 ** 12) * abs(float(w[iw])))))
        # The four statements this clause requires to be PRINTED are the ones the source prints.
        # The last one replaces an earlier sentence about `\max_{d}|w_{d}|`: the manuscript now states
        # the stronger fact that the stored vector IS `-10^{12}w` apart from rounding and the rebuilt
        # first coordinate, which is what (b) and (f) below check on the two shipped files.
        pres = all(at(lines, t, ln) for t, ln in
                   ((r"\langle c_{\vR},w\rangle>0", g5ln[0]),
                    (r"\lambda=\lfloor-10^{12}w\rceil", g5ln[1]),
                    (r"accumulator bound $2^{62}$ at the chosen scale", g5ln[2]),
                    (r"the stored vector is $-10^{12}w$", g5ln[3])))
        clause("G5", aligned and d1 and scale_ok and repair_ok and sign_ok and bound_ok
               and mx_ok and pres,
               "the 2^30 dual as section 11 describes it, recomputed from %s and %s: the cell's, the"
               " certificate's and the separator's rows agree %s and the first row is d = 1 %s;"
               " lambda_d = round(-10^12 w_d) at all %d rows but d = 1 %s (%d mismatches); at d = 1 the"
               " rounded %d was discarded and the stored %d is the repair %s; <c_O,w> = %.6e > 0 while"
               " <c_O,lambda> = %d = -J < 0 %s; sum_d |lambda_d| = %d <= 2^62 %s; max|lambda| = %d and"
               " max|w| = %s are attained at the same row d = %d and differ by the scale alone %s;"
               " the four statements printed on L%d, L%d, L%d, L%d %s"
               % (os.path.basename(cp), os.path.basename(sp), "yes" if aligned else "NO",
                  "yes" if d1 else "NO", lam.size - 1, "yes" if scale_ok else "NO", nmis,
                  disc if disc is not None else -1, rep, "yes" if repair_ok else "NO", dot_w,
                  dot_lam, "yes" if sign_ok else "NO", asum, "yes" if bound_ok else "NO",
                  abs(int(lam[il])), repr(abs(float(w[iw]))), int(rows30[il]),
                  "yes" if mx_ok else "NO", g5ln[0], g5ln[1], g5ln[2], g5ln[3],
                  "yes" if pres else "NO"))
    say()

    # the line the construction's four kappa <= bounds are printed on, located by its literal
    k7cln = 1286
    # the appendix's diagnostic sentence (located above, where K3 accounts for the width it repeats)
    k7dln = DIAGLN

    # ---------------- K7: the relation symbols, their direction and their truth -------------
    # The digits a clause finds on a line say nothing about what the line asserts of them.  A blind
    # reader reversed the first row's `0.1723893<\kappa` to `0.1723893>\kappa` -- a claim the stored
    # certificate contradicts -- and every clause above still held.  What follows reads the relation.
    # For each anchored inequality the symbol (or the word the prose uses for one) is parsed, and the
    # symbol TOGETHER WITH the side the literal sits on fixes what is asserted: `0.17 < kappa` makes
    # the literal a LOWER bound, `kappa < 0.18` an UPPER one, and `>` swaps both.  The assertion is
    # then tested against the certified enclosure of that quantity, on exact rationals.
    LTSYM = {"<": True, "\\le": True, "\\leq": True, "at most": True, "below": True,
             ">": False, "\\ge": False, "\\geq": False, "at least": False, "above": False}
    STRICTSYM = {"<", ">", "below", "above"}
    SYMRE = r"\\leq|\\geq|\\le|\\ge|<|>"

    def bound_kind(sym, num_left):
        """which bound the literal is asserted to be: 'lower' or 'upper'."""
        return "lower" if (LTSYM[sym] == num_left) else "upper"

    def justified(val, kind, strict, lo_c, hi_c):
        """is a printed `kind` bound `val` justified by the certified enclosure [lo_c, hi_c]?"""
        if kind == "lower":
            return False if lo_c is None else (val < lo_c if strict else val <= lo_c)
        return False if hi_c is None else (hi_c < val if strict else hi_c <= val)

    k7ok, k7d = True, []

    def rel(what, sym, num_left, val, lo_c, hi_c):
        """record one parsed relation and whether the certificates justify it."""
        nonlocal k7ok
        if sym not in LTSYM:
            k7ok = False
            k7d.append("%s: relation UNPARSED" % what)
            return
        kind = bound_kind(sym, num_left)
        strict = sym in STRICTSYM
        ok = justified(val, kind, strict, lo_c, hi_c)
        k7ok = k7ok and ok
        k7d.append("%s: %s%s asserts a %s %s bound %s"
                   % (what, "literal " if num_left else "", sym,
                      "strict" if strict else "non-strict", kind,
                      "and holds" if ok else "FAILS"))

    # (a) the table's header: column 2 is a lower bound on kappa, column 3 an upper one
    hdr = lines[TAB[CELLS[0]] - 2]
    hcols = [c.strip() for c in hdr.split("&")]
    hsym = [re.search(r"\\kappa\(2\^\{e\}\)\s*(" + SYMRE + ")", c) for c in hcols[1:3]]
    if len(hcols) < 4 or not all(hsym):
        k7ok = False
        k7d.append("the table header on L%d does not parse" % (TAB[CELLS[0]] - 1))
        hkind = (None, None)
    else:
        hkind = tuple(bound_kind(m.group(1), False) for m in hsym)
        hok = hkind == ("lower", "upper")
        k7ok = k7ok and hok
        k7d.append("header L%d: column 2 %s, column 3 %s %s"
                   % (TAB[CELLS[0]] - 1, hsym[0].group(1), hsym[1].group(1),
                      "as the columns are used" if hok else "REVERSED"))

    # (b) each row: the two bound columns under that header, and the two-sided seven-decimal cell
    for e in CELLS:
        row = lines[TAB[e] - 1]
        rc = [c.strip() for c in row.split("&")]
        for j, kindwant in ((1, hkind[0]), (2, hkind[1])):
            m = re.fullmatch(r"\$([0-9]+\.[0-9]+)\$", rc[j]) if j < len(rc) else None
            if m is None or kindwant is None:
                k7ok = False
                k7d.append("2^%d column %d does not parse as a single decimal" % (e, j + 1))
                continue
            v = Fraction(m.group(1)) ** 2
            ok = justified(v, kindwant, False, C[e]["lo2"], C[e]["hi2"])
            k7ok = k7ok and ok
            k7d.append("2^%d column %d %s is read as %s %s bound %s"
                       % (e, j + 1, m.group(1), "a" if kindwant == "lower" else "an",
                          kindwant, "and holds" if ok else "FAILS"))
        m4 = re.search(r"\$([0-9]+\.[0-9]+)\s*(" + SYMRE + r")\s*\\kappa\s*(" + SYMRE
                       + r")\s*([0-9]+\.[0-9]+)\$", rc[3] if len(rc) > 3 else "")
        if m4 is None:
            k7ok = False
            k7d.append("2^%d: the two-sided seven-decimal statement does not parse" % e)
            continue
        rel("2^%d %s%s kappa" % (e, m4.group(1), m4.group(2)), m4.group(2), True,
            Fraction(m4.group(1)) ** 2, C[e]["lo2"], C[e]["hi2"])
        rel("2^%d kappa%s%s" % (e, m4.group(3), m4.group(4)), m4.group(3), False,
            Fraction(m4.group(4)) ** 2, C[e]["lo2"], C[e]["hi2"])

    # (c) the printed 2^30 lower bound, against the certified rational J^2/(H D_p)
    if None in (J30, H30, Dp_s):
        k7ok = False
        k7d.append("2^30: not scored, no run record")
    else:
        m30 = re.search(r"\\kappa\(2\^\{30\}\)(?:\\ |\s)*(" + SYMRE + r")(?:\\ |\s)*([0-9]+\.[0-9]+)",
                        lines[g3ln - 1])
        if m30 is None:
            k7ok = False
            k7d.append("2^30: the printed bound on L%d does not parse" % g3ln)
        else:
            rel("2^30 kappa%s%s on L%d" % (m30.group(1), m30.group(2), g3ln), m30.group(1), False,
                Fraction(m30.group(2)) ** 2, Fraction(J30 * J30, H30 * Dp_s), None)

    # (d) the statements about the two largest widths.  Section 11 states them as "at most"; the
    # appendix's diagnostic paragraph then states the internal target as MET at one size and NOT met
    # at the other, and for that size names the width section 11 states instead, by a cross-reference
    # to that section -- three relations, one of them negative, so all three are tested here against
    # the certified widths, and the cross-reference is resolved rather than taken on trust.
    wtxt = " ".join(lines[WLN[16] - 1:WLN[16] + 2])
    wm = re.findall(r"(at most|at least)\s*\$([0-9]+\.[0-9])\\cdot10\^\{(-[0-9]+)\}\$\s*at\s*\$e=([0-9]+)\$",
                    wtxt)
    dtxt = " ".join(lines[k7dln - 1:k7dln + 1])
    dm = re.search(r"internal target width \$10\^\{(-[0-9]+)\}\$.*?met at \$2\^\{([0-9]+)\}\$"
                   r" but not at \$2\^\{([0-9]+)\}\$, whose width is nevertheless below the"
                   r" \$([0-9]+\.[0-9])\\cdot10\^\{(-[0-9]+)\}\$ stated in"
                   r" \\S\\ref\{([A-Za-z0-9:_-]+)\}", dtxt)
    if len(wm) != 2 or dm is None:
        k7ok = False
        k7d.append("the width statements do not parse (%d scaled on L%d-%d, the diagnostic sentence"
                   " on L%d-%d %s)"
                   % (len(wm), WLN[16], WLN[16] + 2, k7dln, k7dln + 1,
                      "parses" if dm else "does NOT parse"))
    else:
        for word, lit, ex, ec in wm:
            e = int(ec)
            rel("width(2^%d) %s %s*10^%s" % (e, word, lit, ex), word, False,
                Fraction(lit) * Fraction(10) ** int(ex), C[e]["width"], C[e]["width"])
        tgt, emet, enot = Fraction(10) ** int(dm.group(1)), int(dm.group(2)), int(dm.group(3))
        rel("width(2^%d) meets the target 10^%s" % (emet, dm.group(1)), "below", False,
            tgt, C[emet]["width"], C[emet]["width"])
        notmet = C[enot]["width"] > tgt       # the NEGATIVE claim: the target is not met at 2^28
        k7ok = k7ok and notmet
        k7d.append("width(2^%d) does NOT meet the target 10^%s: %s"
                   % (enot, dm.group(1), "and that is what the source says" if notmet
                      else "FAILS -- the source says it is not met but it is"))
        # the cross-reference: the figure the appendix attributes to another section must be the one
        # that section states for this instance, and the label must be the label of the section the
        # statement is printed in -- otherwise "stated in X" points at a figure X does not state.
        lit, ex, lab = dm.group(4), dm.group(5), dm.group(6)
        ebig = Fraction(lit) * Fraction(10) ** int(ex)
        quoted = any(w == "at most" and l == lit and x == ex and int(c) == enot
                     for w, l, x, c in wm)
        sec = re.findall(r"\\section\{[^}]*\}\\label\{([A-Za-z0-9:_-]+)\}",
                         "\n".join(lines[:WLN[enot]]))
        here = sec[-1] if sec else None
        xref = quoted and here == lab
        k7ok = k7ok and xref
        k7d.append("the %s*10^%s it attributes to \\S\\ref{%s} is stated there %s (section 11 is"
                   " \\label{%s}, and its own statement of that width is on L%d) %s"
                   % (lit, ex, lab, "yes" if quoted else "NO", here, WLN[enot],
                      "ok" if xref else "FAILS"))
        rel("width(2^%d) below the %s*10^%s stated in \\S\\ref{%s}" % (enot, lit, ex, lab),
            "below", False, ebig, C[enot]["width"], C[enot]["width"])

    # (e) the construction's own bounds: upper bounds on kappa at four instances, in one statement
    cm = re.search(r"\\kappa(" + SYMRE + r")([0-9.,\\ ]+)\$\s*at\s*\$e=([0-9,]+)\$",
                   lines[k7cln - 1])
    if cm is None:
        k7ok = False
        k7d.append("the construction's statement on L%d does not parse" % k7cln)
    else:
        cvals = re.findall(r"[0-9]+\.[0-9]+", cm.group(2))
        ccells = [int(x) for x in cm.group(3).split(",")]
        if len(cvals) != len(ccells):
            k7ok = False
            k7d.append("the construction prints %d bounds for %d cells on L%d"
                       % (len(cvals), len(ccells), k7cln))
        else:
            for lit, e in zip(cvals, ccells):
                rel("construction 2^%d kappa%s%s" % (e, cm.group(1), lit), cm.group(1), False,
                    Fraction(lit) ** 2, C[e]["lo2"], C[e]["hi2"])

    clause("K7", k7ok, "the direction and the truth of every anchored inequality, against the"
                       " certified rationals: " + "; ".join(k7d))
    say()

    # ---------------- X1, X2: the crossover -------------------------------------------------
    c15 = mp.sqrt(15 / mp.pi ** 2)
    c1 = (mp.mpf(4) / mp.sqrt(6)) * c15
    c2 = 4 * mp.sqrt(2 * mp.e ** (mp.e - 1))
    f = lambda L: c1 / (L * mp.sqrt(mp.log(L))) + c2 * L * L * mp.e ** (-L / (24 * mp.log(L)))
    Ls = mp.findroot(lambda L: f(L) - 1, mp.mpf(3767))
    L2N = 2 * Ls / mp.log(2)
    say("   crossover: L* = %s  (f(L*) = %s),  2L*/log2 = %s"
        % (mp.nstr(Ls, 16), mp.nstr(f(Ls), 8), mp.nstr(L2N, 16)))
    x1ln = 1135
    x1 = trunc(Ls, 2) == "3766.52" and at(lines, "3766.52", x1ln)
    cover("X1", "3766.52", x1ln)
    clause("X1", x1, "L* = %s, truncated at 2 dp = %s vs 3766.52\\ldots on L%d %s"
           % (mp.nstr(Ls, 12), trunc(Ls, 2), x1ln, "ok" if x1 else "NO"))
    t2, r2 = trunc(L2N, 2), rnd(L2N, 2)
    x2ln = 1135
    x2 = t2 == "10867.87" and at(lines, "10867.87", x2ln) and not re.search(r"10867\.88", body)
    cover("X2", "10867.87", x2ln)
    clause("X2", x2, "2L*/log2 = %s, truncated at 2 dp = %s vs 10867.87\\ldots on L%d %s; the"
                     " rounded value %s, which that notation would misstate, does not occur in the"
                     " source %s"
           % (mp.nstr(L2N, 14), t2, x2ln, "ok" if (t2 == "10867.87" and at(lines, "10867.87", x2ln))
              else "NO", r2, "yes" if not re.search(r"10867\.88", body) else "NO"))

    # ---------------- X3: the two sufficient thresholds, rounded UPWARD ----------------------
    # A sufficient threshold rounded to nearest can fall below the quantity it is sufficient for.
    # Both figures on that line are the crossing rounded up at two decimals, and the clause checks
    # the direction as well as the digits.  `49.05` stays covered by v2's V18b; X3 claims only the
    # literal v2 no longer matches at the precision the source prints.
    L6 = mp.findroot(lambda L: 6 * mp.log(L) - L, mp.mpf(17))
    L12 = mp.findroot(lambda L: 12 * mp.log(L) - L, mp.mpf(46))
    t6, t12 = 2 * L6 / mp.log(2), 2 * L12 / mp.log(2)
    # the upward two-decimal rounding, formed on the integer ceil so no float formatting enters
    x3ln = (1146, 1146)              # the line the two sufficient thresholds are printed on
    up2 = lambda v: "%d.%02d" % divmod(int(mp.ceil(v * 100)), 100)
    u6, u12 = up2(t6), up2(t12)
    x3 = (u6 == "49.05" and u12 == "132.51"
          and at(lines, "49.05", x3ln[0]) and at(lines, "132.51", x3ln[1])
          and mp.mpf("49.05") >= t6 and mp.mpf("132.51") >= t12
          and mp.mpf("132.5") < t12
          and not re.search(r"(?<![0-9.])132\.5(?![0-9.])", body))
    cover("X3", "132.51", x3ln[1])
    clause("X3", x3, "the two sufficient thresholds, rounded upward: 2L6/log2 = %s -> %s on L%d,"
                     " 2L12/log2 = %s -> %s on L%d; each printed figure >= its crossing %s; the"
                     " one-decimal 132.5 is below the crossing and absent from the source %s"
           % (mp.nstr(t6, 12), u6, x3ln[0], mp.nstr(t12, 12), u12, x3ln[1],
              "yes" if (mp.mpf("49.05") >= t6 and mp.mpf("132.51") >= t12) else "NO",
              "yes" if (mp.mpf("132.5") < t12
                        and not re.search(r"(?<![0-9.])132\.5(?![0-9.])", body)) else "NO"))

    # ---------------- B1, N1 ----------------------------------------------------------------
    g = lambda vv: 1 / vv - mp.log(vv)
    b_at2 = g(mp.mpf(2))
    mono = bool(mp.diff(g, mp.mpf("1.5")) < 0 and mp.diff(g, mp.mpf("1.05")) < 0
                and mp.diff(g, mp.mpf("1.95")) < 0)
    b1ln = 700
    b1 = (trunc(b_at2, 4) == "-0.1931" and at(lines, "0.1931", b1ln)
          and mp.almosteq(g(mp.mpf(1)), 1, 1e-30) and mono)
    cover("B1", "0.1931", b1ln)
    clause("B1", b1, "1/v - log v: value 1 at v = 1, %s at v = 2, truncated at 4 dp = %s vs"
                     " -0.1931\\ldots on L%d; derivative negative across [1,2] %s"
           % (mp.nstr(b_at2, 12), trunc(b_at2, 4), b1ln, "yes" if mono else "NO"))
    n1ln = (640, 734)
    n1 = (Fraction(3, 2) == Fraction("1.5") and at(lines, "1.5", n1ln[0])
          and at(lines, "1.5", n1ln[1]))
    cover("N1", "1.5", *n1ln)
    clause("N1", n1, "log X >= 1.5 L from X >= Q^{3/2}: 3/2 = 1.5 exactly, printed on L%d and L%d"
           % n1ln)

    # ---------------- N2: the band-density remark of section 11 ----------------------------
    # R/N and N/R at each computed cell, from cell_cache's row counts alone -- the same integers
    # pa_certify's C1 controls against the rebuilt band.
    dens = []
    for e in DENS_CELLS:
        S = load_cell(e)
        dens.append((e, int(S["N"]), int(S["nR"])))
    rr = [Fraction(R, N) for _, N, R in dens]
    qq = [Fraction(N, R) for _, N, R in dens]
    # the ratios are exact Fractions; v2's rounding helpers take an mpf, so convert at the display
    # step only -- every comparison below stays on the exact rationals
    fr = lambda x: mp.mpf(x.numerator) / mp.mpf(x.denominator)
    lim = mp.pi ** 2 / (3 * (1 - mp.log(2)))
    inc = all(rr[i] < rr[i + 1] for i in range(len(rr) - 1))
    dec = all(qq[i] > qq[i + 1] for i in range(len(qq) - 1))
    under = all(x < Fraction(9, 100) for x in rr)
    above = all(fr(x) > lim for x in qq)
    # the limit is printed as `10.7213228\ldots`, its SEVEN-decimal truncation; the four-decimal
    # form the earlier source printed is no longer there, so the clause checks the precision the
    # source actually gives
    n2lit = [("0.0768", 1292, rnd(fr(rr[0]), 4)), ("0.0834", 1292, rnd(fr(rr[-1]), 4)),
             ("13.03", 1294, rnd(fr(qq[0]), 2)), ("11.99", 1294, rnd(fr(qq[-1]), 2)),
             ("10.7213228", 1294, trunc(lim, 7))]
    digits = all(lit == got for lit, _, got in n2lit)
    placed = all(at(lines, lit, ln) for lit, ln, _ in n2lit)
    for lit, ln, _ in n2lit:
        cover("N2", lit, ln)
    # the three line numbers the message quotes are the ones the clause tested, taken from the list
    # above rather than written out again
    n2ln = [ln for _, ln, _ in n2lit]
    clause("N2", digits and placed and inc and dec and under and above,
           "the density remark over %d cells (2^%d..2^%d), from cell_cache row counts: R/N %s -> %s"
           " vs 0.0768/0.0834 on L%d, N/R %s -> %s vs 13.03/11.99 on L%d, limit truncated %s vs"
           " 10.7213228 on L%d; digits %s, on their lines %s; R/N strictly increasing %s and < 0.09"
           " throughout %s; N/R strictly decreasing %s and above the limit throughout %s"
           % (len(dens), dens[0][0], dens[-1][0], rnd(fr(rr[0]), 6), rnd(fr(rr[-1]), 6),
              n2ln[0], rnd(fr(qq[0]), 4), rnd(fr(qq[-1]), 4), n2ln[2], trunc(lim, 7), n2ln[4],
              "ok" if digits else "NO", "yes" if placed else "NO", "yes" if inc else "NO",
              "yes" if under else "NO", "yes" if dec else "NO", "yes" if above else "NO"))
    say()

    # ---------------- R1: the verifier runtimes ---------------------------------------------
    def jread(n):
        p = os.path.join(RES, n)
        return json.load(io.open(p, encoding="utf-8")) if os.path.exists(p) else None

    # R1 is a PROVENANCE claim: the printed figures are the ones the authors' own verifier run
    # recorded, on the hardware the record names.  It therefore reads a FROZEN file that no step of the
    # reproduction writes -- `pa_certify_timings_frozen.json`, which also carries the CPU, core count,
    # RAM, operating system and library versions.  An earlier arrangement read the ordinary run record,
    # which step 1 overwrites, so a reader's own timings replaced the figures the clause was meant to
    # check and the clause failed for a reason that was not a defect in the paper.  A reader's own
    # record is reported beside the frozen one and is deliberately not scored.
    frozen = jread("pa_certify_timings_frozen.json")
    mine = jread("pa_certify.json")
    if mine and mine.get("mode") != "verify":
        # step 3 of the packet's Reproduce leaves a CERTIFY record in that file, which has no
        # per-cell verify times at all.  Only a verify record is a re-measurement of step 1.
        mine = None
    # The four figures are read OUT OF THE SOURCE rather than hardcoded, so this clause states a
    # relation -- "what section 11 prints is what the frozen record measured" -- and keeps holding
    # when the printed figures are updated, instead of needing an edit here each time.
    # The timings live in APPENDIX A.  The restructuring left a paragraph of the same name in section
    # 11 -- a pointer paragraph with no figures in it -- so taking the first match in the body would
    # read the wrong block and report "0 timing figures".  The search therefore starts at the appendix's
    # own section line, located by its title, and takes the paragraph inside it.
    APPSEC = "\\section{Verification of the certificates}"
    app = body.index(APPSEC) if APPSEC in body else 0
    key = "\\paragraph{Verification cost.}"
    sec = body.index(key, app) if key in body[app:] else 0
    # it is the LAST paragraph of the body (what follows it is the bibliography), so the block ends at
    # the next `\paragraph{` when there is one and at the end of the body when there is not
    nxt = body.find("\\paragraph{", sec + 10) if sec else -1
    endp = nxt if nxt != -1 else len(body)
    para = body[sec:endp]
    printed = re.findall(r"\$([0-9]+\.[0-9])\$\s*(?:s|seconds?)\b", para)
    pl = []
    for tok in printed:
        ln = next((i + 1 for i, x in enumerate(lines) if "$%s$" % tok in x and sec_line(lines, i, sec,
                                                                                      body)), None)
        if ln:
            cover("R1", tok, ln)
        pl.append((tok, ln))
    GONE = ("16.3", "7.6", "31.4")             # superseded figures, no longer in the source
    if frozen is None:
        clause("R1", False, "results/pa_certify_timings_frozen.json is missing, so the printed"
                            " runtimes cannot be checked against the authors' measurement")
    elif len(printed) != 4:
        clause("R1", False, "appendix A's verifier paragraph prints %d timing figures, expected 4"
                            " (found %s)" % (len(printed), ", ".join(printed) or "none"))
    else:
        want = [rnd(float(frozen["wall_s"]), 1), rnd(float(frozen["five_small_cells_s"]), 1),
                rnd(float(frozen["cell_2e26_s"]), 1), rnd(float(frozen["cell_2e28_s"]), 1)]
        digits = printed == want
        placed = all(ln is not None for _, ln in pl)
        # The three parts must be a breakdown of that one run, not three independent walls -- and the
        # comparison is EXACT.  `want` holds the four figures as the one-decimal strings the record and
        # the manuscript print; `Fraction` reads each as the decimal it is, and the tolerance is the
        # rational 2/10.  Each part is a one-decimal rounding and so carries up to 0.05 of rounding
        # error, which lets three of them differ from the total by 0.15 by rounding alone; 2/10 is that
        # slack at the precision printed.  In binary floats the boundary was not the boundary: 70.5
        # against 70.7 evaluates to 0.20000000000000284 and failed a test it meets exactly, and a
        # measurement was once set aside over it.  Nothing here depends on the machine's representation.
        parts_q = sum((Fraction(w) for w in want[1:]), Fraction(0))
        total_q = Fraction(want[0])
        gap_q = abs(parts_q - total_q)
        breakdown = gap_q <= Fraction(2, 10)
        p10 = parts_q * 10
        parts_s = ("%d.%d" % divmod(int(p10), 10)) if p10.denominator == 1 else str(parts_q)
        gone = [t for t in GONE if re.search(r"(?<![0-9.])" + re.escape(t) + r"(?![0-9.])", body)]
        hw = frozen.get("hardware", {})
        # The manuscript no longer prints the CPU: it attributes the figures to one machine and
        # points at the archived timing record, which is where the hardware is specified.  The
        # provenance claim is the same and now runs through one more link, so both ends are checked:
        # the paragraph must make the attribution AND name the record, and the record must carry the
        # hardware it is said to specify.
        attrib = bool(re.search(r"author'?s' ?machine|that machine", para))
        points = bool(re.search(r"timing record", para))
        named = attrib and points and bool(hw.get("cpu"))
        # the record must name the bytes it measured: the verifier THIS run imported, not a revision
        # of it.  Without this the clause can hold while the figures belong to other code.
        vpath = getattr(PC, "__file__", "") or ""
        vsha = (hashlib.sha256(io.open(vpath, "rb").read()).hexdigest()
                if os.path.isfile(vpath) else "")
        fsha = str(frozen.get("verifier_sha256", ""))
        bound = bool(vsha) and fsha == vsha
        ok = digits and placed and breakdown and not gone and named and bound
        if mine and mine.get("mode") == "verify":
            m = {int(v["e"]): float(v["verify_s"]) for v in mine["verified"] if "verify_s" in v}
            here = ("; this machine's own run of step 1 recorded wall %.1f s, e<=24 %.1f s, 2^26 %.1f s,"
                    " 2^28 %.1f s -- reported, not scored"
                    % (float(mine["wall_s"]), sum(m.get(e, 0.0) for e in (16, 18, 20, 22, 24)),
                       m.get(26, -1), m.get(28, -1)))
        else:
            here = ""
        clause("R1", ok, "appendix A prints %s s against the frozen record's %s s (%s, %s, %d threads):"
                         " digits %s, each on its line %s; the three parts sum to %s against the"
                         " printed total %s, exactly %s off the tolerance 2/10 %s; the superseded %s are"
                         " gone %s; the text attributes the run to one machine and names the archived"
                         " record, which carries the hardware %s; the record's verifier SHA-256"
                         " is the shipped pa_certify.py's %s%s"
               % (" / ".join(printed), " / ".join(want), hw.get("cpu", "?"), hw.get("cores", "?"),
                  int(frozen.get("threads", 0)), "ok" if digits else "NO",
                  "yes" if placed else "NO", parts_s, want[0], gap_q, "yes" if breakdown else "NO",
                  "/".join(GONE), "yes" if not gone else "NO: " + ", ".join(gone),
                  "yes" if named else "NO",
                  "yes (%s)" % vsha[:12] if bound else
                  "NO: the record names %s, the shipped verifier is %s"
                  % (fsha[:12] or "nothing", vsha[:12] or "unreadable"), here))
    say()

    # ---------------- W1, W2: v2 re-run, with its outputs redirected ------------------------
    say("re-running pa_tex_constants_v2 (imported, not copied; its stored record is not touched)")
    v2lines: list[str] = []
    TC2.say = lambda s="": v2lines.append(s)
    TC2.OUTJ = os.path.join(tempfile.gettempdir(), "pa_tex_constants_v2_rerun_by_v3.json")
    TC2.main()
    v2bad = list(TC2.BAD)
    v2mism = list(TC2.MISM)
    clause("W1", not v2bad and not v2mism,
           "v2's %d clauses re-run against the current tree: %d fail%s, %d notation mismatches%s"
           % (len(TC2.REC), len(v2bad), (" (" + ", ".join(v2bad) + ")") if v2bad else "",
              len(v2mism), (" (" + ", ".join(m[0] for m in v2mism) + ")") if v2mism else ""))

    # "132.5" is the one-decimal rendering of the Y <= sqrt Q threshold.  The source now prints the
    # two-decimal upward rounding 132.51 of the same crossing (clause X3), which is why the
    # one-decimal literal is no longer there.
    REMOVED = {"0.172389352", "0.140860454", "0.1055426", "0.076331", "0.052762", "0.034622",
               "0.024526", "0.0763313", "0.0346234", "0.0245279", "3767", "10868", "50",
               "132.5",
               # v2's V17c: "a factor about 29 in L beyond the point where the formal expression
               # dips below 1".  Section 9 now states the two scales (the crossing and where
               # u reaches 1) and not the ratio between them, so the literal is gone; both scales
               # are still printed and X1/X2/v2's V15-V17b still check them.
               "29"}
    v2lits = {r["literal"] for r in TC2.REC if r.get("literal")}
    absent, unexplained = [], []
    for lit in sorted(v2lits):
        head = lit.split("e")[0] if re.match(r"^[0-9.]+e-?[0-9]+$", lit) else lit
        # presence is tested against the BODY TEXT, not against the token sweep: the sweep only
        # collects decimals, \,-groups and integers of five or more digits, so v2's short integer
        # literals (`852`, `29`) are outside it and would look absent when they are printed.
        if not re.search(r"(?<![0-9.])" + re.escape(head) + r"(?![0-9.])", body):
            absent.append(lit)
            if lit not in REMOVED and head not in REMOVED:
                unexplained.append(lit)
    clause("W2", not unexplained,
           "v2 literals no longer printed: %s -- all in the recorded removed list %s"
           % (", ".join(absent) or "none", "yes" if not unexplained else
              "NO: " + ", ".join(unexplained)))
    say()

    # ---------------- Q1: quoted classical constants ----------------------------------------
    QUOTED = {("0.2795", 388): "Trudgian Thm. 2, the explicit li-error coefficient",
              ("6.455", 388): "Trudgian Thm. 2, the exponent's denominator",
              ("6.455", 392): "the same constant, in the deduction v2's V14a checks",
              ("1.25506", 375): "Rosser-Schoenfeld (3.6), the Chebyshev coefficient"}
    q1ok, q1d = True, []
    for (t, ln), why in QUOTED.items():
        p = at(lines, t, ln)
        q1ok = q1ok and p
        q1d.append("%s L%d (%s) %s" % (t, ln, why, "present" if p else "NOT PRESENT"))
        COVER[(t, ln)] = "Q1"
    clause("Q1", q1ok, "constants quoted from a cited source, checked for presence and attribution"
                       " only (Trudgian's 0.2795, 6.455 and x >= 229; Rosser-Schoenfeld's"
                       " 1.25506): " + "; ".join(q1d))

    # ---------------- M1: the classification line -------------------------------------------
    # The MSC line and the keywords are a CLASSIFICATION, not a measurement: `2020` names the scheme's
    # edition and `11N35` names a subject area.  No digit group of either carries a quantity, and none
    # is a numeric literal in this file's sense -- `PAT` matches decimals, `\,`-grouped integers and
    # integers of five digits or more, and an MSC code is two digits, a letter and two digits -- so the
    # inventory counts nothing from these lines.  This clause turns that silence into a declared and
    # CHECKED category: the codes are exactly the five declared here, in the primary/secondary split
    # the line prints, each of MSC shape, and the block contributes zero occurrences to the inventory.
    # Nothing is exempted: a numeral smuggled into this block would be counted by the inventory and
    # would have to be covered by a clause like any other.
    MSC_YEAR, MSC_PRIMARY, MSC_SECONDARY = "2020", ("11N35",), ("11N36", "11N25", "90C05")
    mstart = next((i for i, ln in enumerate(lines, 1)
                   if "Mathematics Subject Classification" in ln), None)
    if mstart is None:
        clause("M1", False, "no classification line: the manuscript prints no"
                            " \"Mathematics Subject Classification\"")
    else:
        mend = next((i for i in range(mstart, min(mstart + 6, len(lines) + 1))
                     if "\\par}" in lines[i - 1]), mstart)
        block = " ".join(lines[mstart - 1:mend])
        codes = re.findall(r"\b([0-9]{2}[A-Z][0-9]{2})\b", block)
        shape = all(re.fullmatch(r"[0-9]{2}[A-Z][0-9]{2}", c) for c in codes)
        prim = re.search(r"Classification:\s*(.*?)\s*\(primary\)", block)
        seco = re.search(r"\(primary\);\s*(.*?)\s*\(secondary\)", block)
        gp = tuple(re.findall(r"\b([0-9]{2}[A-Z][0-9]{2})\b", prim.group(1))) if prim else ()
        gs = tuple(re.findall(r"\b([0-9]{2}[A-Z][0-9]{2})\b", seco.group(1))) if seco else ()
        year = bool(re.search(re.escape(MSC_YEAR) + r"\s+Mathematics Subject Classification", block))
        kw = re.search(r"Keywords:\s*(.*?)(?:\\par\}|$)", block)
        kwtext = kw.group(1) if kw else ""
        kwclean = bool(kw) and not re.search(r"[0-9]", kwtext)
        # the inventory's own verdict on these lines, not an assumption about them
        inv = [t for t, ln, _ in occ if mstart <= ln <= mend]
        clause("M1", bool(year and shape and gp == MSC_PRIMARY and gs == MSC_SECONDARY
                          and kwclean and not inv),
               "the classification line on L%d-%d is a declared non-numeric category: %s MSC with"
               " primary %s %s, secondary %s %s, every code of MSC shape %s; keywords \"%s\" carry no"
               " digit %s; and the two lines contribute %d occurrences to the literal inventory (%s)"
               % (mstart, mend, MSC_YEAR if year else "NO YEAR", ", ".join(gp) or "none",
                  "ok" if gp == MSC_PRIMARY else "EXPECTED " + ", ".join(MSC_PRIMARY),
                  ", ".join(gs) or "none",
                  "ok" if gs == MSC_SECONDARY else "EXPECTED " + ", ".join(MSC_SECONDARY),
                  "yes" if shape else "NO", kwtext[:60].strip().rstrip("."),
                  "yes" if kwclean else "NO", len(inv),
                  "as a classification must" if not inv else "UNEXPECTED: " + ", ".join(inv)))

    # ---------------- Y1: the numerals a reference entry carries ----------------------------
    # A bibliography entry prints a volume, a year and a page range.  They are locators, not
    # quantities, and they are outside the checked body by construction: `read_tex` cuts the body at
    # `\begin{thebibliography}`, so the inventory never sees them and `Z1` is not asked to cover them.
    # This clause turns that construction into a checked statement for the locators the entries
    # actually print: each must appear in the bibliography, must NOT appear as a standalone number
    # anywhere in the body, and must not be a tracked literal of the body.  (A three- or four-digit
    # locator is untracked by `PAT` even where it does appear, which is why the body test is on
    # standalone digits and not only on the token list.)  A locator that migrates into the body fails
    # here and would then have to be covered by a clause like any other -- nothing is exempted.
    BIBNUM = {"33": "the volume of the Friedlander entry",
              "1976": "its year",
              "565--576": "its page range"}
    # Locators printed in the BODY rather than inside a `\cite[...]` argument.  The inventory sees
    # them as ordinary decimals, so they must be accounted for; they are locators and not quantities,
    # which is the category this clause declares.  The declaration is checked rather than asserted:
    # each must sit in a subsection that cites the work it locates.
    BODYLOC = {("3.1", 283): "BradyThesis", ("4.1", 283): "BradyThesis",
               ("7.1", 283): "BradyThesis", ("4.2", 284): "BradyThesis",
               ("8.1", 284): "BradyThesis", ("8.2", 285): "BradyThesis"}
    bibtext = raw[raw.index("\\begin{thebibliography}"):] if "\\begin{thebibliography}" in raw else ""
    btoks = PAT.findall(body)
    y1ok, y1d = bool(bibtext), []
    for t in sorted(BIBNUM):
        inbib = t in bibtext
        stray = bool(re.search(r"(?<![0-9.])" + re.escape(t) + r"(?![0-9.])", body))
        tracked = t in btoks
        ok = inbib and not stray and not tracked
        y1ok = y1ok and ok
        y1d.append("%s (%s): in the bibliography %s, no standalone occurrence in the body %s, not a"
                   " tracked literal of the body %s"
                   % (t, BIBNUM[t], "yes" if inbib else "NO", "yes" if not stray else "NO",
                      "yes" if not tracked else "NO"))
    def subsec(ln):
        """the subsection block containing body line `ln`, as one string."""
        a = max((k for k in range(ln, 0, -1)
                 if lines[k - 1].startswith(("\\subsection", "\\section"))), default=1)
        b = next((k for k in range(ln + 1, len(lines) + 1)
                  if lines[k - 1].startswith(("\\subsection", "\\section"))), len(lines))
        return " ".join(lines[a - 1:b])

    locd = []
    for (t, ln), key in sorted(BODYLOC.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        here = at(lines, t, ln)
        cited = ("\\cite" in subsec(ln)) and (key in subsec(ln))
        ok = here and cited
        y1ok = y1ok and ok
        locd.append("%s@L%d locates \\cite{%s} %s"
                    % (t, ln, key, "ok" if ok else "NO"
                       if here else "NOT ON THAT LINE"))
        COVER[(t, ln)] = "Y1"
    clause("Y1", y1ok, "bibliographic numerals, declared and outside the checked body (which ends at"
                       " \\begin{thebibliography}): " + "; ".join(y1d)
                       + ". Locators printed inside the body, each in a subsection that cites the"
                         " work it locates: " + "; ".join(locd))

    # ---------------- Z1: coverage ----------------------------------------------------------
    # v2 records a `literal` field only for its `clause()` entries.  Its `plain()` entries check
    # literals too -- the embedded decimals of V3/V10b/V13/V14b, the 2^52 integers of V20a/e/f/g, the
    # construction bounds of M3 -- and carry them inside their prose, so they are declared here by
    # clause name instead of parsed out of it.  Z1's accounting must be exact, not approximate.
    V2_PLAIN = {"V3": ("2.01317",), "V10b": ("0.09",), "V13": ("1.26",), "V14b": ("0.39",),
                "V20a": ("67108864",), "V20e": ("67108879",), "V20f": ("67108849",),
                "V20g": ("34260449", "34260493"),
                "M3": ("0.6982", "0.6130", "0.6088", "0.5841")}
    v2names = {r["name"] for r in TC2.REC}
    v2cov = {}
    for r in TC2.REC:
        lit = r.get("literal")
        if not lit:
            continue
        head = lit.split("e")[0] if re.match(r"^[0-9.]+e-?[0-9]+$", lit) else lit
        for t, ln, ref in occ:
            if t == head and not ref:
                v2cov[(t, ln)] = "v2:" + r["name"]
    for nm, lits in V2_PLAIN.items():
        if nm not in v2names:                   # the delegation must name a clause v2 actually ran
            continue
        for lit in lits:
            for t, ln, ref in occ:
                if t == lit and not ref:
                    v2cov[(t, ln)] = "v2:" + nm
    missing_deleg = sorted(set(V2_PLAIN) - v2names)
    # Count OCCURRENCES, not (literal, line) keys.  `COVER` and `v2cov` are keyed by the pair, so a
    # literal printed twice on one line collapses to one key and the reported split did not add up to
    # the total -- a defect a final consistency check found (`34260449` twice on one line, `1.0` twice on
    # another).  Z1's verdict was never wrong; only its arithmetic was unreportable.  Both are now
    # reported: the occurrence split, which sums, and the key counts behind it.
    left, n_here, n_v2 = [], 0, 0
    for t, ln, ref in occ:
        if ref:
            continue
        if (t, ln) in COVER:
            n_here += 1
        elif (t, ln) in v2cov:
            n_v2 += 1
        else:
            left.append((t, ln))
    n_tot = sum(1 for _, _, r in occ if not r)
    n_quoted = sum(1 for t, ln, ref in occ if not ref and COVER.get((t, ln)) == "Q1")
    dup = sorted({(t, ln) for t, ln, ref in occ if not ref
                  and sum(1 for u, m, _ in occ if (u, m) == (t, ln)) > 1})
    stale = [k for k in COVER if k not in {(t, ln) for t, ln, _ in occ}]
    clause("Z1", not left and not stale and not missing_deleg
           and n_here + n_v2 + len(left) == n_tot,
           "coverage of the %d non-reference literal occurrences: %d by the clauses above (%d of them"
           " quoted), %d by v2, %d left over -- %d + %d + %d = %d %s. Those occurrences occupy %d"
           " distinct (literal, line) pairs here and %d at v2, because %d literal%s printed twice on"
           " one line (%s). Left over %s; claimed but not in the .tex %s; delegated to a v2 clause"
           " that did not run %s"
           % (n_tot, n_here, n_quoted, n_v2, len(left), n_here, n_v2, len(left), n_tot,
              "checks out" if n_here + n_v2 + len(left) == n_tot else "DOES NOT ADD UP",
              len(COVER), len(v2cov), len(dup), "" if len(dup) == 1 else "s",
              ", ".join("%s@L%d" % k for k in dup) or "none",
              ("none" if not left else "; ".join("%s@L%d" % k for k in left)),
              ("none" if not stale else "; ".join("%s@L%d" % k for k in stale)),
              ("none" if not missing_deleg else ", ".join(missing_deleg))))

    # ---------------- D1: notation ----------------------------------------------------------
    d1 = [m for m in v2mism]
    extra = []
    # the crossover companion is the one literal caught here; it is checked against the source's
    # own notation, so the test is the printed digits, not a remembered value
    for lit, ln, val in (("3766.52", x1ln, Ls), ("10867.87", x2ln, L2N)):
        if at(lines, lit, ln) and trunc(val, 2) != lit and rnd(val, 2) == lit:
            extra.append(("X%d" % (1 if lit.startswith("3766") else 2), lit, "ldots", "round",
                          trunc(val, 2), rnd(val, 2)))
    clause("D1", not d1 and not extra,
           "literals whose rendering disagrees with the notation the source declares: %s"
           % ("none" if not (d1 + extra) else
              "; ".join("%s prints %s under %s but that is the %s value (trunc %s, round %s)"
                        % (n, lit, nota, which, t, r) for n, lit, nota, which, t, r in d1 + extra)))

    say()
    say("FAILS: %s" % (", ".join(BAD) if BAD else "none"))
    json.dump(dict(tex=dict(path=rel_or_name(TEXP),
                            sha256=sha, bytes=len(raw.encode("utf-8")),
                            body_lines=len(lines), occurrences=len(occ), references=nref),
                   certificates=[dict(e=e, lower=C[e]["lo12"], upper=C[e]["hi12"],
                                      lower7=C[e]["lo7"], upper7=C[e]["hi7"],
                                      width=str(C[e]["width"]), lower2=str(C[e]["lo2"]),
                                      upper2=str(C[e]["hi2"]), rows=C[e]["nrows"],
                                      cols=C[e]["ncols"]) for e in CELLS],
                   v2=dict(clauses=len(TC2.REC), fails=v2bad,
                           mismatches=[list(m) for m in v2mism], absent_literals=absent),
                   clauses=REC, fails=BAD),
              io.open(OUTJ, "w", encoding="utf-8", newline="\n"), indent=1)
    say("wrote %s and %s" % (os.path.basename(OUT), os.path.basename(OUTJ)))
    return 1 if BAD else 0


if __name__ == "__main__":
    sys.exit(main())

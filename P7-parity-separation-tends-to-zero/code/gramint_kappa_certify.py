# -*- coding: utf-8 -*-
r"""gramint_kappa_certify -- one integer certifier for any saved kappa separator

Supporting computation; not used in the paper's proofs or tables.

WHERE THIS COMES FROM

`gramint_kappa_integer.py` certified the saved ladder separator at 2^24 and
`gramint_kappa_hat_integer.py` the hat separators at 2^24 and 2^26, each in
the integer form of the certificate.  They are two scripts because they were written on
two days for two families.  The two families are being regenerated with two
more ladder vectors (`2^24` size 256, `2^26` size 128) with a `lab_encoding`
field in the npz, and each will want certifying the minute it lands.  This is
the one tool for that: it takes an npz path, reads which family it holds from
its contents, reads the label encoding from the file when the file says it,
and runs the same five clauses.

FAMILY DETECTION, from the npz and nothing else

    `g` and `lab` present        ladder: g per class, w_d = g[class(d)]
    `names` and `v` present      hat:    w_d = sum_k c_k * hat_k(u_d) on class(d)

For the ladder family the row-to-class map needs the label encoding.  If the
npz carries `lab_encoding`, it is used and echoed; if not, the pre-fix
encoding `(min(omega,5)*1000 + ev*100 + sb)*10 + lp` is tried FIRST, because
that is what every ladder file saved before 09-12 holds (`certified-gram/
the trap is recorded in this file's header), and the current encoding second.  Whichever
reproduces the npz's `lab` array element for element is the one used, and
the file says which.  `B` is read from the npz if present, else from the
filename's `size<B>`.

THE CLAUSES are unchanged from the earlier two scripts: K0 reconstruction (the
rebuilt labels or hats reproduce the npz, and the float bound recomputes to
the npz's at 1e-09), K1 exact integer feasibility after the `d = 1` repair,
K2 integer bound within 1e-06 of the float, K3 `c^2 H D_p < J^2` at `c` the
float bound rounded down to six decimals, K4 the falsifier.  Non-zero exit on
any.

A REPRODUCTION REQUIREMENT: on the three npz files already
certified it must reproduce `J`, `H` and `D_p` EXACTLY -- the same integers
to the last digit -- and the same six-decimal `c`.  Anything else is a defect
in the consolidation, not a finding.

RESULT FILE.  One line per npz certified, so the file is a ledger of every
separator this tool has certified and the integers behind each.  The ledger
is kept in the sidecar `gramint_kappa_certify_ledger.json` keyed by npz name
(a run with explicit paths updates only those keys) and the text file is
regenerated from the sidecar on every run.  The sidecar is what makes the
ledger survive a run with explicit paths: a text file opened with "w" as the
only store would be overwritten by it.
"""

import hashlib
import importlib.util
import io
import json
import math
import os
import re
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
OUT = os.path.join(ROOT, "results", "gramint_kappa_certify.txt")
LEDGER = os.path.join(ROOT, "results", "gramint_kappa_certify_ledger.json")

CAP = 5
DEN = 10 ** 12
TOL_REPRO = 1e-09
TOL_BOUND = 1e-06

_lines = []


def say(s=""):
    print(s, flush=True)
    _lines.append(s)


class _LazyModule:
    """Defer loading a helper module until one of its attributes is actually used.

    Only the ladder and hat families need the ladder helper.  Importing it unconditionally at start-up
    meant that every family failed with an ImportError when that helper was not beside this script --
    which is the case wherever this file is distributed without it.  Now a run that never touches the
    ladder family never looks for it, and a run that does gets a sentence saying what is missing
    instead of a traceback.
    """

    def __init__(self, name, fn):
        self._name, self._fn, self._mod = name, fn, None

    def __getattr__(self, attr):
        if self._mod is None:
            p = os.path.join(HERE, self._fn)
            if not os.path.exists(p):
                raise SystemExit(
                    "gramint_kappa_certify: this file needs the helper %s for the ladder and hat"
                    " families, and it is not beside this script; those families cannot be certified"
                    " here.  Every other family works without it." % self._fn)
            self._mod = load(self._name, self._fn)
        return getattr(self._mod, attr)


def load(name, fn):
    """import a carried module with sys.argv hidden from it.

    `conedual_ladder_with_repair.py` parses `sys.argv[1:]` as cell exponents
    at import time, so loading it while this tool's own argument -- an npz
    path -- is in argv raised ValueError.  Found by a pre-test on a
    re-encoded copy before the real files landed; the argument-free default
    run had passed and the argument-bearing one, which is the only one that
    matters, had not.
    """
    saved = sys.argv
    sys.argv = [saved[0]]
    try:
        spec = importlib.util.spec_from_file_location(
            name, os.path.join(HERE, fn))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    finally:
        sys.argv = saved
    return mod


def cell_of(path):
    m = re.search(r"_2e(\d+)", os.path.basename(path))
    if not m:
        raise SystemExit("cannot read the cell from the filename: " + path)
    return int(m.group(1))


# ------------------------------------------------------------- ladder family
def ladder_labels(S, L, B, encoding):
    """encoding is one of prefix, current, prefix+NOLP, current+NOLP.

    NOLP is the ladder script's environment switch that drops the least-prime
    digit; files saved under it end every label in 0.  The 2^22 files on disk
    turned out to be three different variants -- size 64 and 256 current+NOLP,
    size 128 prefix WITH lp -- and the last is the case the ladder script's
    own comment warns about: at B >= 128 the pre-fix encoding lets the size
    bin carry into the parity digit and MERGES classes.  Such a file is still a
    partition and its certificate is still valid; it is the coarsened
    partition, and its bound is lower than the intended rung's.
    """
    rows, om, Q = S["rows"], S["rom"], S["Q"]
    ev = (rows % 2 == 0).astype(np.int64)
    nolp = encoding.endswith("+NOLP")
    lp = np.zeros(rows.size, dtype=np.int64) if nolp else         np.array([L.lpc(int(d)) for d in rows])
    sb = np.minimum((B * np.log(np.maximum(rows, 2)) / math.log(Q))
                    .astype(np.int64), B - 1)
    if encoding.startswith("prefix"):
        lab = (np.minimum(om, CAP) * 1000 + ev * 100 + sb) * 10 + lp
    else:
        lab = ((np.minimum(om, CAP) * 2 + ev) * B + sb) * 10 + lp
    lab[rows == 1] = -1
    u = np.array([v for v in np.unique(lab) if v != -1])
    cls = np.zeros(rows.size, dtype=np.int64)
    for i, v in enumerate(u):
        cls[lab == v] = i + 1
    return np.concatenate(([-1], u)), cls, u.size + 1


def ladder_w(S, z, L, path):
    """per-row w from the class vector g, and the encoding that reproduced lab."""
    if "B" in z.files:
        B = int(z["B"])
    else:
        m = re.search(r"size(\d+)", os.path.basename(path))
        if not m:
            raise SystemExit("B is neither in the npz nor the filename")
        B = int(m.group(1))
    if "lab_encoding" in z.files:
        declared = str(z["lab_encoding"])
        base = "current" if "*2" in declared or "* 2" in declared else "prefix"
        order = [base, base + "+NOLP"]
    else:
        declared = None
        order = ["prefix", "current", "current+NOLP", "prefix+NOLP"]
    for enc in order:
        mine, cls, C = ladder_labels(S, L, B, enc)
        if np.array_equal(mine, z["lab"]):
            g = z["g"].astype(np.float64)
            if C != g.size:
                raise SystemExit("C = %d classes but g has %d" % (C, g.size))
            w = g[cls]
            # the pre-fix encoding CAN merge (odd, sb = 100+k) with (even,
            # sb = k) at B >= 128 -- but only for the same least-prime digit,
            # and even rows have lp = 0 while odd rows have lp >= 1, so with
            # lp carried no populated collision exists.  Collisions are
            # COUNTED below, not inferred from "prefix and B >= 128" alone:
            # the 2^22 size-128 file in the current encoding has C = 408 and
            # the same J, H, M to the digit, so that condition by itself does
            # not mean the classes merged.
            ncoll = 0
            if enc.startswith("prefix"):
                _, _, C_cur = ladder_labels(S, L, B, "current" + enc[len("prefix"):])
                ncoll = C_cur - C
            note = ("encoding %s%s, B = %d, C = %d%s"
                    % (enc, " (declared in the npz)" if declared else
                       " (decoded; the npz carries no lab_encoding)", B, C,
                       ("  [pre-fix encoding: %d populated label collisions, so "
                        "this is a MERGED partition and its bound is below the "
                        "intended rung's]" % ncoll) if ncoll > 0 else
                       ("  [pre-fix encoding, no populated collision: the "
                        "partition is the intended one]" if enc.startswith("prefix")
                        else "")))
            return w, note
    raise SystemExit("no label encoding reproduces the npz's lab array")


# ---------------------------------------------------------------- hat family
def hat_w(S, z, L):
    rows, om, Q = S["rows"], S["rom"], S["Q"]
    names, v = list(z["names"]), z["v"].astype(np.float64)
    ev = (rows % 2 == 0).astype(np.int64)
    oc = np.minimum(om, CAP)
    lpr = np.array([L.lpc(int(d)) for d in rows])
    # KHAT is read from the basis names (h4_..., h8_...), not assumed: the
    # producing script (conedual_plateau_primal.py) scales u by log2/KHAT
    # and names the functions by the same KHAT, so a wrong K fails K0.
    ks = sorted({int(m.group(1)) for m in
                 (re.match(r"h(\d+)_", n) for n in names) if m})
    if len(ks) != 1:
        raise SystemExit("hat family with KHAT set %r -- not one hat scale" % ks)
    KHAT = ks[0]
    uu = np.log(Q / rows.astype(np.float64)) / (math.log(2.0) / KHAT)
    w = np.zeros(rows.size)
    pat = re.compile(r"h%d_(\d+)\[om=(\d),(\w+),lp(\d)\]" % KHAT)
    used = 0
    for name, c in zip(names, v):
        m = pat.match(name)
        if not m:
            continue
        k, o, par, q = int(m.group(1)), int(m.group(2)), m.group(3), int(m.group(4))
        mask = (oc == o) & (ev == (1 if par == "even" else 0)) & (lpr == q)
        w += float(c) * mask * np.maximum(0.0, 1.0 - np.abs(uu - k))
        used += 1
    w[0] = float(v[0])
    return w, "%d hat functions, KHAT = %d" % (used, KHAT)


# --------------------------------------------------------------- the sweep
def sweep(N, rows, vals, idx_sets, dtype):
    acc = np.zeros(N + 1, dtype=dtype)
    for i in range(1, rows.size):
        x = vals[i]
        if x != 0:
            d = int(rows[i])
            acc[d::d] += x
    out = [acc[s].copy() for s in idx_sets]
    del acc
    return out


def certify(path, T, L):
    start = len(_lines)
    r = _certify(path, T, L)
    r["block"] = _lines[start:]
    r["sha"] = hashlib.sha256(io.open(path, "rb").read()).hexdigest()[:12]
    r["block"].insert(1, "  sha256[:12] = %s" % r["sha"])
    return r


def _certify(path, T, L):
    z = np.load(path)
    e_ = cell_of(path)
    family = "ladder" if "g" in z.files and "lab" in z.files else \
             "hat" if "names" in z.files and "v" in z.files else None
    if family is None and "w" in z.files and "rows" in z.files:
        family = "rowvec"
    if family is None:
        raise SystemExit("neither a ladder, a hat nor a row-vector npz: " + path)
    # rowvec : a per-row separator w = t - A y from the full-column
    # nnls (gramint_nnls_full).  The npz carries the nnls UPPER bound `ub`, not
    # a lower-bound `bound`; K0 for this family is polar feasibility -- the
    # unrepaired max_vB W relative to ||w||^2 is small -- and the float bound is
    # compared with the KKT-optimal value ||w||/||p|| it would equal at
    # convergence (printed, not asserted).
    bound_npz = float(z["bound"]) if "bound" in z.files else float("nan")
    t0 = time.time()
    S = T.build(e_)
    rows, vB, vR, N = S["rows"], S["vB"], S["vR"], S["N"]
    nR = int(vR.size)
    say("%s" % os.path.basename(path))
    say("  2^%d  %s family  Q = %d  rows = %d  |vB| = %d  |vR| = %d  "
        "built %.0f s" % (e_, family, S["Q"], rows.size, vB.size, nR,
                          time.time() - t0))
    if family == "rowvec":
        assert np.array_equal(z["rows"].astype(np.int64), rows.astype(np.int64)), "row order differs"
        w = z["w"].astype(np.float64)
        if "kkt" in z.files:
            note = "per-row separator from nnls (%d iterations, ub %.7f, KKT %.2e)" % (
                int(z["iters"]), float(z["ub"]), float(z["kkt"]))
        else:
            # a repaired per-row vector with its own float bound (research
            # main's conedual_targeted_repair): K0 is the
            # standard reproduction of that bound
            note = "per-row separator, repaired (src %s, tau %s, float bound %.7f)" % (
                str(z["src"]) if "src" in z.files else "?",
                ("%.4f" % float(z["tau"])) if "tau" in z.files else "?", float(z["bound"]))
    else:
        w, note = ladder_w(S, z, L, path) if family == "ladder" else hat_w(S, z, L)
    say("  reconstruction: %s" % note)

    # K0: float certificate recomputed the family's own way
    WBf, WRf = sweep(N, rows, w, (vB, vR), np.float64)
    Mf = float(WBf.max())
    dotf = (float(WRf.sum()) - Mf * nR) / nR
    nwf = math.sqrt(Mf * Mf + float(np.dot(w[1:], w[1:])))
    npv = float(np.linalg.norm(S["p"]))
    bound_f = dotf / (nwf * npv)
    if family == "rowvec" and "kkt" not in z.files:
        k0 = abs(bound_f - bound_npz) <= TOL_REPRO
        say("  K0  float bound rebuilt %.12f vs npz %.12f (dev %.2e)  %s   [full max_vB W = w[0] + sweep max = %.3e]"
            % (bound_f, bound_npz, abs(bound_f - bound_npz), "HOLDS" if k0 else "REFUTED", Mf + float(w[0])))
    elif family == "rowvec":
        nw = float(np.linalg.norm(w))
        # the d >= 2 sweep max (2.1e-4) is NOT comparable with ||w||^2
        # (2.6e-5): the sweep excludes
        # the d = 1 entry w[0], which the nnls dual carries at -2.1e-4 exactly
        # so that the FULL W = w[0] + sweep is <= KKT on vB.  Polar
        # feasibility is a statement about the full W; the repair then
        # replaces w[0] by -M, which changes it by the KKT residual only.
        Mfull = Mf + float(w[0])
        # a tolerance of 1e-2 * ||w||^2 does not serve either: the violation is
        # 5% of ||w||^2, the same 5% the bound loses to complementarity.  The
        # reproduction check this family admits is that the rebuilt full W
        # reproduces the SAVED KKT residual: |max_vB W_full| <= 2 * kkt.
        kk = float(z["kkt"])
        k0 = abs(Mfull) <= 2.0 * kk
        say("  K0  rebuilt full max_vB W = w[0] + sweep max = %.3e + %.3e = %.3e vs saved KKT %.3e  %s;"
            "  ||w||^2 = %.3e;  float bound after repair %.12f vs ||w||/||p|| = %.12f"
            " (equal at the KKT optimum; the gap is <A^T w, y> < 0, complementarity not yet exact)"
            % (float(w[0]), Mf, Mfull, kk, "HOLDS" if k0 else "REFUTED", nw * nw, bound_f, nw / npv))
        bound_npz = bound_f
    else:
        k0 = abs(bound_f - bound_npz) <= TOL_REPRO
        say("  K0  float bound rebuilt %.12f vs npz %.12f (dev %.2e)  %s"
            % (bound_f, bound_npz, abs(bound_f - bound_npz),
               "HOLDS" if k0 else "REFUTED"))
    if not k0:
        return dict(ok=False, name=os.path.basename(path), e=e_)

    n = np.rint(w * DEN).astype(np.int64)
    n[0] = 0
    WB, WR = sweep(N, rows, n, (vB, vR), np.int64)
    M = int(WB.max())
    k1 = int((WB - M).max()) == 0
    # an int64 sum here overflows on the 2^26 size-128 ladder (|vR| = 5.3e6
    # entries of size up to M = 2.2e13 -> 1e20 > 2^63) and returns J < 0; K2
    # catches that (bound nan) while K3 -- which squares J -- would report
    # HOLDS on it.  The accumulator itself is
    # safe (max W = M, guarded below); the SUM of W over vR is not, so it is
    # taken in exact Python integers and K3 requires J > 0 explicitly.
    J = sum(WR.tolist()) - M * nR
    # overflow guard for the per-entry accumulator: the int64 sweep of n must
    # agree with the float sweep of w (K0's) to rounding, else it wrapped
    ov = max(float(np.abs(WB / DEN - WBf).max()), float(np.abs(WR / DEN - WRf).max()))
    assert ov < 1e-6 * max(1.0, float(np.abs(WBf).max())) + rows.size / DEN,         "int64 accumulator disagrees with the float sweep: overflow (%.3e)" % ov
    H = M * M + sum(int(x) * int(x) for x in n[1:])
    cb_cr = np.rint(S["p"] * nR).astype(np.int64)
    D_p = sum(int(x) * int(x) for x in cb_cr)
    bound_int = J / math.sqrt(H * D_p) if J > 0 else float("nan")
    k2 = J > 0 and abs(bound_int - bound_f) <= TOL_BOUND
    c6 = math.floor(bound_f * 1e6)
    k3 = J > 0 and c6 * c6 * H * D_p < J * J * 10 ** 12
    nviol = int(((WB + M) > 0).sum())
    k4 = nviol > 0
    say("  K1  max_vB (W - M) = %d  %s" % (int((WB - M).max()),
                                          "HOLDS" if k1 else "REFUTED"))
    say("  K2  integer bound %.12f (dev %.2e)  %s"
        % (bound_int, abs(bound_int - bound_f), "HOLDS" if k2 else "REFUTED"))
    say("  K3  c^2 H D_p < J^2 at c = %d/10^6  %s" % (c6, "HOLDS" if k3 else "REFUTED"))
    say("  K4  M -> -M: %d of %d violate  %s" % (nviol, vB.size, "HOLDS" if k4 else "REFUTED"))
    say("  J = %d  H = %d  D_p = %d  M = %d" % (J, H, D_p, M))
    ok = k0 and k1 and k2 and k3 and k4
    say("  => %s" % (("kappa(2^%d) > %d/10^6 in exact integers" % (e_, c6))
                     if ok else "REFUTED"))
    say()
    return dict(ok=ok, name=os.path.basename(path), e=e_, family=family,
                c6=c6, bound_f=bound_f, bound_int=bound_int, J=J, H=H,
                D_p=D_p, M=M, note=note, nviol=nviol, nB=int(vB.size))


def main():
    paths = sys.argv[1:]
    if not paths:
        paths = [os.path.join(ROOT, "results", f) for f in (
            "conedual_ladder_g_2e22_size64.npz",
            "conedual_ladder_g_2e22_size128.npz",
            "conedual_ladder_g_2e22_size256.npz",
            "conedual_ladder_g_2e24_size64.npz",
            "conedual_hat4u_v_2e24.npz",
            "conedual_hat4u_v_2e26.npz")]
    if os.environ.get("MERGE"):
        # merge-only mode : a certification run on another machine writes
        # its own sidecar there; MERGE=<that sidecar> folds its entries into
        # this ledger and regenerates the text without recertifying (the
        # entries carry the printed block, the integers and the npz sha).
        other = json.load(io.open(os.environ["MERGE"], encoding="utf-8"))
        results = list(other.values())
        for r in results:
            r["block"] = ["  [certified on another machine; merged from %s]" % os.path.basename(os.environ["MERGE"])] + r["block"]
    else:
        T = load("tfor", "conedual_the_four_odd_rows.py")
        L = _LazyModule("ladder", "conedual_ladder_with_repair.py")
        results = [certify(p, T, L) for p in paths]
    bad = [r for r in results if not r["ok"]]

    # the ledger: every npz ever certified, keyed by name; this run's
    # entries replace their keys, the others stay as they were
    ledger = {}
    if os.path.exists(LEDGER):
        ledger = json.load(io.open(LEDGER, encoding="utf-8"))
    for r in results:
        old = ledger.get(r["name"])
        if old is not None and old.get("sha") != r["sha"]:
            # the npz was regenerated under the same name: the earlier
            # certificate stays in the ledger, keyed by its content hash
            ledger[r["name"] + "#" + old.get("sha", "nosha")] = old
        ledger[r["name"]] = r
    io.open(LEDGER, "w", encoding="utf-8", newline='\n').write(
        json.dumps(ledger, indent=1, sort_keys=True) + '\n')
    for key, r in ledger.items():
        r["key"] = key
    order = sorted(ledger.values(), key=lambda r: (r["e"], r["key"]))
    del _lines[:]
    say("gramint_kappa_certify -- integer certificates for saved separators")
    say("=" * 70)
    say("DEN = 10^12; clauses K0-K4 as in the two earlier certify scripts.")
    say("This run: %s" % ", ".join(r["name"] for r in results))
    say()
    for r in order:
        _lines.extend(r["block"])          # into the file only; already printed
    say("=" * 70)
    say("LEDGER  (%d separators; sidecar %s)" % (len(order),
                                                 os.path.basename(LEDGER)))
    say("  %-40s %-6s %-7s %-14s %s" % ("npz", "cell", "family", "exact kappa >",
                                         "float"))
    for r in order:
        if r["ok"]:
            say("  %-40s 2^%-4d %-7s %5d/10^6      %.12f"
                % (r["key"], r["e"], r["family"], r["c6"], r["bound_f"]))
        else:
            say("  %-40s 2^%-4d REFUTED" % (r["key"], r["e"]))
    head = [
        "STATISTIC: for each saved separator npz given on the command line",
        "           (or the six already certified, by default), the clause set's",
        "           integer certificate: the exact feasibility after the",
        "           d = 1 repair, the integers J, H, D_p and M, the bound",
        "           J / sqrt(H D_p), and c^2 H D_p < J^2 at c the float",
        "           bound rounded down to six decimals.",
        "NULL: K4 (the falsifier) and K0 (the rebuilt separator is the saved",
        "      one) are the controls, per certificate.",
        "FIELD: whichever cells the given npz files name; band and vR/vB",
        "       split of conedual_the_four_odd_rows.py; ladder families",
        "       reconstructed in the encoding that reproduces the npz's lab",
        "       array, hat families from the npz's basis names; DEN = 10^12.",
        "DENOM: the bound divides by sqrt(H D_p).",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + _lines) + "\n")
    print("\nwrote %s" % OUT)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

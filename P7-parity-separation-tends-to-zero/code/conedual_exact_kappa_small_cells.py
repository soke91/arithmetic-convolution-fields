"""conedual_exact_kappa_small_cells -- kappa exactly at the small cells by dense nnls, for the ratio series.

Supporting computation; not used in the paper's proofs or tables.

At 2^12..2^20 the cone has <= 624 rows and <= 77k columns, so
dist(t, cone) is a dense nnls (scipy, active set) that converges to machine precision: kappa =
||t - A y*|| / ||p|| with t_d = |{l in vR : d | l}| / |vR| (t_1 = 1) and p from conedual_the_four_odd_rows
.build -- the same normalisation as the earlier dumps, whose ||p|| was checked against it.

REGISTERED.
  (E0) at 2^20 the dense nnls reproduces the converged column-generation value 0.1055427
       to 1e-7  ==>  the small-cell values are in the same convention as the certified enclosures.
  (E1) the consecutive ratios r_k = kappa(2^(2k+2)) / kappa(2^(2k)) are printed; no law is registered
       here -- they are read against the certified ratios at 2^22..2^26 elsewhere.
Cells: argv (default 12 14 16 18); 2^20 takes ~700 s (dense 624 x 76946) and is run on request.
"""
import importlib.util
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import nnls

# `regnote.py` sits beside this script.
_regdir = os.path.dirname(os.path.abspath(__file__))
if not os.path.isfile(os.path.join(_regdir, "regnote.py")):
    raise ImportError("regnote.py must sit beside this script")
if _regdir not in sys.path:
    sys.path.insert(0, _regdir)
from regnote import Note

HERE = os.path.dirname(os.path.abspath(__file__))
CELLS = tuple(int(a) for a in sys.argv[1:]) or (12, 14, 16, 18)
note = Note(os.path.join(HERE, "..", "results", "conedual_exact_kappa_small_cells%s.txt" % ("" if CELLS == (12, 14, 16, 18) else "_" + "_".join(map(str, CELLS)))),
            "conedual_exact_kappa_small_cells -- dense nnls at 2^%s" % ",".join(map(str, CELLS)))
spec = importlib.util.spec_from_file_location("tfor", os.path.join(HERE, "conedual_the_four_odd_rows.py"))
T = importlib.util.module_from_spec(spec); spec.loader.exec_module(T)


def main():
    prev = None
    for E in CELLS:
        S = T.build(E); vB = S["vB"]; vR = S["vR"]; N = S["N"]; npv = float(np.linalg.norm(S["p"]))
        rows = np.array([d for d in range(1, math.isqrt(N) + 1) if all(d % (q * q) for q in range(2, int(d ** 0.5) + 1))])
        A = np.zeros((rows.size, vB.size)); t = np.zeros(rows.size)
        for i, d in enumerate(rows):
            d = int(d); A[i] = (vB % d == 0); t[i] = (vR % d == 0).mean()
        t0 = time.time(); y, res = nnls(A, t, maxiter=400000)
        kap = res / npv
        note.say("2^%d  rows %d  |vB| %d  |vR| %d  ||p|| %.6f  dist %.7f  kappa %.7f%s  (%.0f s)"
                 % (E, rows.size, vB.size, vR.size, npv, res, kap, ("  ratio %.4f" % (kap / prev)) if prev else "", time.time() - t0))
        if E == 20 and abs(kap - 0.1055427) > 1e-7:
            note.fail("E0")
        prev = kap
    note.finish()


if __name__ == "__main__":
    main()

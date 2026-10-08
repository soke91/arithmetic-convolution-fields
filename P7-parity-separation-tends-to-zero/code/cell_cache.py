# -*- coding: utf-8 -*-
r"""cell_cache -- one fast, cached build of a cone-duality cell; every readout loads it from here.

WHY THIS FILE EXISTS

`conedual_the_four_odd_rows.build(e)` is the reference build, but its row counts are a Python
loop over ~10^4 rows (one numpy strided sum each), and every readout re-ran it
each rebuilt 2^28 (100-160 s) on, against the project rule "heavy compute is
compiled, vectorised, parallel, cache-aware" (the project's compute rules).  This module

  * builds the band exactly as `build()` does (same canonical sieve, q1, thr, band, vB/vR split),
  * counts the rows with one numba kernel (prange over rows, cyclic assignment),
  * saves the small per-row arrays once per cell under results/cell_cache_2e<e>/ and reloads
    them memory-mapped afterwards (a 2^28 reload is milliseconds).

A new script calls `load(e)`, not the slow build() directly; that is enforced, not left to
recall.

RETURNED (dict; arrays are read-only memmaps on reload)

  N, Q, q1, thr, nB, nR                 ints
  rows   int64   squarefree d <= Q (the row order of build())
  rom    int64   omega(d) on rows
  cb, cr int64   #{n in vB : d | n}, #{n in vR : d | n}
  p      float64 (cb - cr)/|vR|        -- equal to build(e)["p"] element for element
  t      float64 cr/|vR|

    from cell_cache import load
    S = load(28)                         # builds once (~1.6 GB peak at 2^28), then cached
    python code/cell_cache.py 16 18 20 22 24 --check    # equality with build()
"""
from __future__ import annotations
import io, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
NTHREAD = int(os.environ.get("ACF_THREADS", "8"))
_SCRATCH = os.environ.get("TEMP") or os.environ.get("TMP") or "."
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(_SCRATCH, "numba_cache_acf"))

import numba  # noqa: E402
import numpy as np  # noqa: E402
from numba import njit, prange, set_num_threads  # noqa: E402

# Clamp to what numba will actually accept, BEFORE any build runs.  `set_num_threads(n)` raises
# ValueError when n > NUMBA_NUM_THREADS (numba's limit: the env var if set, else the core count), so
# the 8 above aborted the build on a 4-core machine and on any run that caps the env var at 4 -- the
# README tells the reader to do exactly that.  Measured in a packet copy:
# `cell 2^16 FAILED: ValueError('The number of threads must be between 1 and 4')`, no cache written,
# and `conedual_cstar_sweep.py` then died with KeyError 'Q' while scoring the cell it had lost.
# The clamp is numerically inert: `_count_rows` splits the rows cyclically (`range(w, m, nthread)`),
# a disjoint cover for every `nthread`, and each row's cb/cr is written by exactly one worker.
# It is deliberately placed OUTSIDE `_build`, which `fingerprint()` hashes -- editing `_build` would
# invalidate every cell cache ever made (fingerprint dfdad9cc7934832e, the one the reproduction packet ships).
NTHREAD = max(1, min(NTHREAD, int(numba.config.NUMBA_NUM_THREADS)))

sys.path.insert(0, HERE)
from _sieve_shared import mu_upto, primes_upto  # noqa: E402

VERSION = 1   # bump when the saved fields change
FIELDS = ("rows", "rom", "cb", "cr")
# the sieve body the fingerprint hashes: the one _sieve_shared imported (it searches upward, so this
# works both in the repository and in a paper packet, where code/ and lib/ are siblings -- the packet layout)
SIEVE_SRC = os.path.join(sys.modules["_sieve_shared"]._REPO, "lib", "goldbach", "sieve.py")


@njit(parallel=True, cache=True, nogil=True)
def _count_rows(code, rows, cb, cr, nthread):
    """cb[i] = #{n : code[n] == 1, rows[i] | n}, cr[i] likewise for code 2.

    Work of row i is N/rows[i]; a stride-nthread slice of the row list gives each worker ~1/nthread.
    """
    n = code.shape[0] - 1
    m = rows.shape[0]
    for w in prange(nthread):
        for i in range(w, m, nthread):
            d = rows[i]
            a = 0
            b = 0
            for j in range(d, n + 1, d):
                c = code[j]
                if c == 1:
                    a += 1
                elif c == 2:
                    b += 1
            cb[i] = a
            cr[i] = b


def _build(e):
    """The band of conedual_the_four_odd_rows.build(e), byte for byte; only the row counts differ in method."""
    N = 2 ** e
    Q = int(N ** 0.5)
    assert N < 2 ** 31, "cell_cache: int32 remainder array wraps at 2^31"
    mu = np.array(mu_upto(N), dtype=np.int8)
    rem = np.arange(N + 1, dtype=np.int32)
    om = np.zeros(N + 1, dtype=np.int8)
    for p in primes_upto(Q):          # one strided op per prime <= Q, as build() does
        p = int(p)
        om[p::p] += 1
        q = p
        while q <= N:
            rem[q::q] //= p
            if q > N // p:
                break
            q *= p
    q1 = Q + 1
    while any(q1 % r == 0 for r in range(2, int(q1 ** 0.5) + 1)):
        q1 += 1
    thr = N // q1
    band = (mu != 0) & (rem == 1)
    band[:thr + 1] = False
    del rem
    code = np.zeros(N + 1, dtype=np.int8)
    np.copyto(code, np.int8(1), where=band & (om % 2 == 0))
    np.copyto(code, np.int8(2), where=band & (om % 2 == 1))
    del band
    nB = int(np.count_nonzero(code == 1))
    nR = int(np.count_nonzero(code == 2))
    rows = (np.nonzero(mu[1:Q + 1] != 0)[0] + 1).astype(np.int64)
    rom = om[rows].astype(np.int64)
    del om, mu
    cb = np.zeros(rows.size, dtype=np.int64)
    cr = np.zeros(rows.size, dtype=np.int64)
    set_num_threads(NTHREAD)
    _count_rows(code, rows, cb, cr, NTHREAD)
    del code
    return dict(N=N, Q=Q, q1=q1, thr=thr, nB=nB, nR=nR, rows=rows, rom=rom, cb=cb, cr=cr)


def cache_dir(e):
    return os.path.join(RES, "cell_cache_2e%d" % e)


def _derive(S):
    nR = float(S["nR"])
    S["p"] = (S["cb"] - S["cr"]).astype(np.float64) / nR
    S["t"] = S["cr"].astype(np.float64) / nR
    return S


def fingerprint():
    """What the cached numbers depend on: this file's build code and the canonical sieve's bytes.

    Any edit to either changes the fingerprint, and every cache made before it is ignored and
    rebuilt -- no one has to remember to bump VERSION after changing the band.
    """
    import hashlib, inspect
    h = hashlib.sha256()
    h.update(str(VERSION).encode())
    h.update(inspect.getsource(_build).encode("utf-8"))
    h.update(inspect.getsource(_count_rows.py_func).encode("utf-8"))
    with open(SIEVE_SRC, "rb") as f:
        h.update(f.read())
    return h.hexdigest()[:16]


def _read(d, e):
    """The cached cell, or None with the reason when it must not be used."""
    meta_f = os.path.join(d, "meta.json")
    if not os.path.exists(meta_f):
        return None, "no cache"
    try:
        meta = json.load(io.open(meta_f, encoding="utf-8"))
        if meta.get("fingerprint") != fingerprint():
            return None, "build code or sieve changed since the cache was made"
        if meta.get("e") != e:
            return None, "cache is for another cell"
        S = {k: meta[k] for k in ("N", "Q", "q1", "thr", "nB", "nR")}
        for k in FIELDS:
            S[k] = np.load(os.path.join(d, k + ".npy"), mmap_mode="r")
        m = meta["rows_count"]
        if any(S[k].shape != (m,) for k in FIELDS) or int(S["cr"].sum()) != meta["cr_sum"]:
            return None, "cache arrays incomplete or damaged"
        if S["rows"][0] != 1 or S["cr"][0] != S["nR"] or S["cb"][0] != S["nB"]:
            return None, "cache inconsistent (row 1 must count the whole band)"
    except Exception as exc:   # unreadable cache is rebuilt, never trusted
        return None, "cache unreadable (%s)" % type(exc).__name__
    return S, "cache"


def load(e, rebuild=False, quiet=False):
    """The cell at N = 2^e, from results/cell_cache_2e<e>/ when that is safe, else built (and saved).

    The cache is NOT used, and the cell is built fresh, when:
      * rebuild=True, or the environment has ACF_NO_CACHE=1 (then nothing is read or written) --
        use this when the build itself is what is being verified;
      * the cache was made by different build code or a different canonical sieve (fingerprint);
      * the cache files are missing, incomplete, damaged or inconsistent.
    The reason is printed (to stderr) unless quiet.  Only the per-row arrays are cached; a script
    that needs the band itself (vB, vR, the N-length code array) builds with `_build`-equivalent code.
    """
    no_cache = os.environ.get("ACF_NO_CACHE") == "1"
    d = cache_dir(e)
    if not (rebuild or no_cache):
        S, why = _read(d, e)
        if S is not None:
            if not quiet:
                print("cell_cache 2^%d: CACHE USED (%s)" % (e, os.path.relpath(d, ROOT)), file=sys.stderr, flush=True)
            return _derive(S)
    else:
        why = "ACF_NO_CACHE=1" if no_cache else "rebuild requested"
    if not quiet:
        print("cell_cache 2^%d: CACHE NOT USED, building (%s)" % (e, why), file=sys.stderr, flush=True)
    t0 = time.time()
    S = _build(e)
    if no_cache:
        return _derive(S)
    tmp = d + ".tmp%d" % os.getpid()
    os.makedirs(tmp, exist_ok=True)
    for k in FIELDS:
        np.save(os.path.join(tmp, k + ".npy"), S[k])
    meta = {k: int(S[k]) for k in ("N", "Q", "q1", "thr", "nB", "nR")}
    meta.update(version=VERSION, e=e, fingerprint=fingerprint(), rows_count=int(S["rows"].size),
                cr_sum=int(S["cr"].sum()), build_seconds=round(time.time() - t0, 1), threads=NTHREAD)
    json.dump(meta, io.open(os.path.join(tmp, "meta.json"), "w", encoding="utf-8"), indent=1)
    try:
        if os.path.isdir(d):          # stale or damaged cache: replace whole
            for f in os.listdir(d):
                os.remove(os.path.join(d, f))
            os.rmdir(d)
        os.replace(tmp, d)
    except OSError:                   # another process built the same cell first: keep theirs
        for f in os.listdir(tmp):
            os.remove(os.path.join(tmp, f))
        os.rmdir(tmp)
    return _derive(S)


def _check(cells):
    """Equality with the reference build(): rows, |vB|, |vR| and p element for element."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("tfor", os.path.join(HERE, "conedual_the_four_odd_rows.py"))
    saved = sys.argv
    sys.argv = [saved[0]]
    try:
        T = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(T)
    finally:
        sys.argv = saved
    ok = True
    for e in cells:
        t0 = time.time()
        S = load(e, rebuild=True)
        t1 = time.time()
        R = T.build(e)
        t2 = time.time()
        same = (np.array_equal(S["rows"], R["rows"]) and S["nB"] == int(R["nB"]) and S["nR"] == int(R["nR"])
                and np.array_equal(S["p"], np.asarray(R["p"], dtype=np.float64)))
        ok = ok and same
        print("2^%d  equal %s   cell_cache %.1f s   build() %.1f s" % (e, same, t1 - t0, t2 - t1), flush=True)
    print("CHECK", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cells = [int(a) for a in args] or [16, 18, 20, 22, 24]
    if "--check" in sys.argv:
        sys.exit(_check(cells))
    for e in cells:
        t0 = time.time()
        S = load(e)
        print("2^%d  rows %d  |vB| %d  |vR| %d   %.1f s" % (e, S["rows"].size, S["nB"], S["nR"], time.time() - t0))

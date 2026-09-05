# -*- coding: utf-8 -*-
"""pass5 -- independent primitives for the code audit.

Ground truth is trial-division factorisation, not a sieve: the point is to
have a second implementation that shares no code path with the target.
"""
import math
import os
from functools import lru_cache

# The tree whose scripts are under audit.  Overridable so the pass can be
# re-run against a copy held anywhere; it defaults to this repository's
# own code/, which is where those scripts live once a pass is over.
CODE_DIR = os.environ.get(
    "AUDIT_CODE_DIR",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "..", "..", "code"),
)


def factor(n):
    """trial division -> dict p -> exponent."""
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def mu_bf(n):
    if n < 1:
        return 0
    if n == 1:
        return 1
    f = factor(n)
    if any(e > 1 for e in f.values()):
        return 0
    return -1 if len(f) % 2 else 1


def Lambda_bf(n):
    if n < 2:
        return 0.0
    f = factor(n)
    if len(f) != 1:
        return 0.0
    p = next(iter(f))
    return math.log(p)


def phi_bf(n):
    if n < 1:
        return 0
    r = n
    for p in factor(n):
        r -= r // p
    return r


def is_prime_bf(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1 if d == 2 else 2
    return True


def load(name):
    """import one of the audited scripts as a module (does not run main)."""
    import importlib.util
    import os
    import sys
    if CODE_DIR not in sys.path:
        sys.path.insert(0, CODE_DIR)
    path = os.path.join(CODE_DIR, name + ".py")
    spec = importlib.util.spec_from_file_location("t_" + name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

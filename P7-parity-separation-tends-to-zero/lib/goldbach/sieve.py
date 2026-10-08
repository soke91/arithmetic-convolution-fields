"""Sieve of Eratosthenes with numpy: primes, and the shared `mu` and `Lambda`.

`sieves(n)` returns the primes up to n, the von Mangoldt function and the Moebius function on [0, n].
It is the single implementation every script uses (through `_sieve_shared`), so that no script builds
the Moebius function on its own.
"""

import math

import numpy as np

_cache = {"limit": 0, "is_prime": None, "primes": None, "plimit": 0}


def sieve(limit: int) -> np.ndarray:
    """Boolean primality array for 0..limit (index = number)."""
    if limit <= _cache["limit"] and _cache["is_prime"] is not None:
        return _cache["is_prime"][: limit + 1]
    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[:2] = False
    for p in range(2, int(limit**0.5) + 1):
        if is_prime[p]:
            is_prime[p * p :: p] = False
    _cache.update(limit=limit, is_prime=is_prime, primes=None)
    return is_prime


def primes_upto(limit: int) -> np.ndarray:
    """Array of all primes up to limit."""
    # The limit of the cached prime list (plimit) is tracked separately from the sieve limit, so that a
    # call with a small limit followed by one with a larger limit does not return a truncated list.
    if limit <= _cache["plimit"] and _cache["primes"] is not None:
        ps = _cache["primes"]
        return ps[ps <= limit]
    ps = np.flatnonzero(sieve(limit))
    _cache["primes"] = ps
    _cache["plimit"] = limit
    return ps


def is_prime(n: int) -> bool:
    """Primality of a single n (trial division outside the sieved range)."""
    if n < 2:
        return False
    if _cache["is_prime"] is not None and n <= _cache["limit"]:
        return bool(_cache["is_prime"][n])
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def sieves(n):
    pr = primes_upto(n)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(n + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lgp[i]
            if q > n // p:
                break
            q *= p
    mu = np.ones(n + 1, dtype=np.int8)
    rem = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            rem[q::q] //= p
            if q > n // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return pr, lam, mu


def mu_upto(n):
    """The Moebius function on [0, n], computed by `sieves` (the single implementation).

    This also computes Lambda, roughly doubling the work, but keeps one construction of mu.
    """
    return sieves(n)[2]

def spf_upto(n, dtype=np.int32):
    """Smallest-prime-factor array on [0, n].

    Primes p are run only up to sqrt(n): every composite has a prime factor below sqrt(n), and the entries
    still zero afterwards are themselves prime and are filled in one step. The caller chooses `dtype`
    (downstream integer arithmetic may need int64). Same limit as `sieves`: n < 2^31.
    """
    spf = np.zeros(n + 1, dtype=dtype)
    for p in range(2, math.isqrt(n) + 1):
        if spf[p] == 0:
            blk = spf[p::p]
            spf[p::p] = np.where(blk == 0, p, blk)
    idx = np.flatnonzero(spf[2:] == 0) + 2
    spf[idx] = idx
    return spf

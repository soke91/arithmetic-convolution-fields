"""Bridge to the shared sieve: the one place where the import path is adjusted.

Scripts in this directory import the sieve as

    from _sieve_shared import sieves        # (primes, Lambda, mu)
    from _sieve_shared import mu_upto       # when only mu is needed

instead of each building the Moebius function inline. This module defines no sieve of its own: it only
locates `lib/goldbach/sieve.py` (searching upward from this file, so it works both in the repository and
inside a reproduction packet) and re-exports its functions, so there is exactly one implementation.
"""
import os
import sys

# Locate the directory that holds lib/goldbach/sieve.py by searching upward, rather than counting levels:
# in the repository it is two levels above this file, in a reproduction packet (code/ and lib/ are siblings)
# it is one level above.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = None
_dir = _HERE
while True:
    if os.path.isfile(os.path.join(_dir, "lib", "goldbach", "sieve.py")):
        _REPO = _dir
        break
    _up = os.path.dirname(_dir)
    if _up == _dir:
        break
    _dir = _up
if _REPO is None:
    # Never fall back silently to another implementation: stop if the shared sieve is not found.
    raise ImportError(
        "lib/goldbach/sieve.py not found in any directory above " + _HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from lib.goldbach.sieve import (  # noqa: E402
    mu_upto, primes_upto, sieves, spf_upto)

# Limit: `sieves(n)` is correct only for n < 2^31. It tracks the unfactored part in
# `rem = np.arange(n + 1, dtype=np.int32)`; past int32 numpy wraps to negative values without raising, so
# `rem > 1` fails there and `mu` would silently miss a sign flip. No script here sieves that far (the arrays
# alone would need more than 2^31 bytes).

__all__ = ["sieves", "mu_upto", "primes_upto"]

# `primes_upto` is exported together with `sieves` so that every caller uses the same helper the sieve
# itself uses, rather than a local copy.

"""regnote.py -- the scaffolding every registered-prediction research script needs
(shared by the research scripts; it ships beside them in a reproduction packet).

Purpose: one place for the five things the arithmetic-convolution-fields scripts
re-implemented in ~80 files: an incremental result writer (so a killed run keeps
every line it printed), rounding that is honest by construction (lower bounds
DOWN, upper bounds UP, by integer arithmetic, never %f), a FAILS accumulator for
clauses, a wall-clock stamp, and a UTF-8 console.

Not in this file: any science.  Import it; do not copy it.

    # with regnote.py on sys.path (beside the calling script in a packet):
    from regnote import Note, fdown, fup
    note = Note("results/my_instrument.txt")      # every note.say() rewrites the file
    note.say("2^22  bound %.7f" % fdown(bound))   # a lower bound, printed rounded down
    note.fail("R0 at 2^22")                       # a clause that failed
    note.finish()                               # prints and writes "FAILS: ..." (or "none")

Why these rules:
  * a number computed while writing prose is not a computed number -- print it;
  * a certified bound printed with %f can round ABOVE what was certified; floor it;
  * gridding the SQUARE of a bound instead of the bound lost an order of magnitude;
  * a run killed at hour three should not lose hours one and two: write incrementally.
"""
from __future__ import annotations

import io
import math
import os
import sys
import time


def fdown(x: float, digits: int = 7) -> float:
    """Round DOWN to `digits` decimals (for lower bounds)."""
    m = 10 ** digits
    return math.floor(x * m) / m


def fup(x: float, digits: int = 7) -> float:
    """Round UP to `digits` decimals (for upper bounds)."""
    m = 10 ** digits
    return math.ceil(x * m) / m


def utf8_console() -> None:
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass


class Note:
    """Incremental result file: every `say` appends a line and rewrites the file."""

    def __init__(self, path: str, header: str | None = None):
        self.path = path
        self.lines: list[str] = []
        self.fails: list[str] = []
        self.t0 = time.time()
        utf8_console()
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        if header:
            self.say(header)
            self.say("=" * min(len(header), 100))

    def say(self, s: str = "") -> None:
        self.lines.append(s)
        print(s, flush=True)
        with io.open(self.path, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(self.lines) + "\n")

    def fail(self, clause: str) -> None:
        self.fails.append(clause)

    def elapsed(self) -> str:
        return "%.0f s" % (time.time() - self.t0)

    def finish(self) -> int:
        self.say("FAILS: %s" % (", ".join(self.fails) if self.fails else "none"))
        return 1 if self.fails else 0

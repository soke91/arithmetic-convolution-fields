# A fixed class for the Möbius side

Reproduction packet for the paper

> **P1-mobius-fixed-class.tex**

This directory holds the source, the compiled PDF, the code the paper
cites, and the output that code produced.

## Running it

```
python -m pip install -r ../requirements.txt
python code/<script>.py
```

Python 3.12.10 and NumPy 2.2.5 produced the figures in the paper. Each
script writes the file of the same name in `results/`, and imports only
the standard library and NumPy — never another script here.

**Twelve scripts, and none reads another's output.** The packet is
closed under what it consumes, and that was checked by running every
script and confirming that none dies for a missing file, not by reading
the source.

**Some scripts exit nonzero, and that is their registered behaviour.** A
run whose own preregistered rule was refuted says so in its exit code as
well as in its output. Forty-two rules carry a verdict across the twelve
runs and thirteen were refuted; seven of the twelve record at least one.
All are reported in the paper.

The results committed here are the ones the paper prints. Rerunning a
script overwrites its result file, so keep a copy if the committed
output is wanted for reference.

## One class of figure has no script

The citation counts, the reference-list lengths and the exponents quoted
from other papers in the literature section are read from those papers
and from the two databases named there, on the date named there. No
script produces them, and the paper says so in the same place.

## What is machine-checked, and what is not

The finite rearrangements the argument runs on are formalised in Lean 4
against Mathlib, in `lean/GoldbachLean/Fixedclass.lean`: that the
subtracted mean term leaves the sum over moduli, the exchange and the
reindexing that compose into the divisor switch, the split of the
divisor sum into a complete sum minus a tail, and the factorisation
`u = mk`. They are stated over arbitrary functions, so what is checked
is that these steps use no arithmetic property of the arithmetic
functions — which is what the paper's claim, that its two functionals
differ in one factor and nothing else, rests on.

Two of the statements there are meant to be contradicted rather than
used. One exhibits a subtracted term depending on the modulus for which
the first step fails; the other exhibits index sets carried onto each
other by `n ↦ N − n` for which the switch still fails, so that the
congruence-to-divisibility translation is seen to carry that step rather
than the reindexing. A file with only the positive statements would not
have said either.

Every theorem in that file depends on `propext`, `Classical.choice` and
`Quot.sound` and on nothing else; `lean/axioms.txt` is the report and
`cd lean && lake build` reproduces it. The report covers the other
papers' developments in this repository as well, since they are one
library; where `sorryAx` appears it belongs to the cell-mean paper,
which says so in its own text.

**The theorem itself is not machine-checked.** It rests on
Bombieri–Vinogradov, which Mathlib does not carry, and the truncation
choice, the complete divisor sum and the degeneracy lemma are
arithmetic and are not formalised either. No statement in this paper has
had a reading by a subject expert. The paper states all of this in its
own appendix rather than leaving it to be discovered here.

## Provenance

The paper is a projection of a working repository that is not public.
What that repository keeps and this directory does not is the history —
corrections, withdrawn figures, the record of how each result was
arrived at. What is here is what is true, together with the evidence for
it and the preregistered checks that failed.

## Status

Not submitted and not published.

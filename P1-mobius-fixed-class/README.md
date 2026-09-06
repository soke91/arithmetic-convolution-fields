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

## What is not machine-checked

Nothing here. Other packets in this repository carry a `lean/`
directory; this one does not, no statement in this paper has a Lean
development, and none has had a reading by a subject expert. The proofs
are as given in the text. The paper states this in its own appendix
rather than leaving it to be discovered here.

## Provenance

The paper is a projection of a working repository that is not public.
What that repository keeps and this directory does not is the history —
corrections, withdrawn figures, the record of how each result was
arrived at. What is here is what is true, together with the evidence for
it and the preregistered checks that failed.

## Status

Not submitted and not published.

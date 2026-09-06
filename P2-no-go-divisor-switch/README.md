# A no-go for the divisor-switch route

Reproduction packet for the paper

> **P2-no-go-divisor-switch.tex**

This directory holds the source, the compiled PDF, the code the paper
cites, and the output that code produced.

## Running it

```
python -m pip install -r ../requirements.txt
python code/<script>.py
```

Python 3.12.10 and NumPy 2.2.5 produced the figures in the paper. Each
script writes the file of the same name in `results/`, and imports only
the standard library, NumPy, and — in one case — `indep.py`, which is
in this directory for that reason.

**Sixteen files for fourteen cited scripts.** The paper names fourteen;
`lab_layer_tail` is here because `lab_combined_modulus` reads its
output, and `indep` because `c02_flatsum_k1` imports it. The packet is
closed under both kinds of dependency, and that was checked by running
every script in dependency order and confirming that none dies for a
missing file. The import edge in particular was found by running rather
than by reading: a check that looked only at result files passed while
the packet was broken.

Run `lab_layer_tail` before `lab_combined_modulus`; the rest are
independent of each other.

**Some scripts exit nonzero, and that is their registered behaviour.** A
run whose own preregistered rule was refuted says so in its exit code as
well as in its output. Eight of the fourteen cited runs record such a
refutation, and the paper names all eight.

The results committed here are the ones the paper prints. Rerunning a
script overwrites its result file, so keep a copy if the committed
output is wanted for reference.

## One figure worth knowing about

The paper compares a measured cancellation gain `G` against what
independent signs would buy. That reference is *not* the square root of
the number of moduli: writing `x_k = (log k) H(N;k)`, the reference is
`||x||_1 / ||x||_2`, which equals `sqrt(#k)` only when every `|x_k|` is
the same and is smaller otherwise by Cauchy–Schwarz. Here it is `0.66`
to `0.69` of the count-based figure, so a count of moduli overstates the
reference by about a half. `lab_positive_weights.py` prints both columns
and has done so throughout; an earlier printing of the paper read the
wrong one. The conclusion is unchanged — the measured `G` is an order of
magnitude below the corrected reference too — and only the comparison
moved.

## What is machine-checked, and what is not

The six finite rearrangements that lay out the design space are
formalised in Lean 4 against Mathlib, in `lean/GoldbachLean/Layers.lean`.
They take `Λ`, `μ` and `log` as arbitrary functions, so what is verified
is that these identities use no arithmetic property of them — which is
what makes them rearrangements rather than estimates. Every theorem in
that file depends on `propext`, `Classical.choice` and `Quot.sound` and
on nothing else; `lean/axioms.txt` is the report, and it covers the
other papers' developments in this repository as well, since they are
one library.

`sorryAx` on a line of that report would mark a statement written out
and left unproved. Where it appears it belongs to the cell-mean paper,
which says so in its own text. No line of `Layers.lean` carries it.

**The no-go itself is not machine-checked.** It rests on
Bombieri–Vinogradov, which Mathlib does not carry, so no part of the
theorem is formalised — only the identities that state what space the
theorem is about. Neither theorem of this paper has had a reading by a
subject expert.

To check the Lean rather than take the report on trust:

```
cd lean && lake build
```

The toolchain is pinned in `lean-toolchain` and the Mathlib revision in
`lake-manifest.json`.

## Provenance

The paper is a projection of a working repository that is not public.
What that repository keeps and this directory does not is the history —
corrections, withdrawn figures, the record of how each result was
arrived at. What is here is what is true, together with the evidence for
it and the preregistered checks that failed.

## Status

Not submitted and not published.

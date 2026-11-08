# P10-short-variable-mobius reproduction packet

Purpose: canonical working and public packet for `P10-short-variable-mobius.tex`.

Not in this file: ownership judgement or hashes; `PACKET.json` is authoritative for both.

This directory is self-contained for paper reading and evidence reproduction. It
contains 0 Python files and 0 committed result files. Files
under `code/` and `results/` are generated from the exact sources and hashes in
`PACKET.json`; do not edit those copies.

## Reproduce

This paper is a proof and carries no computation: it quotes no measured number, so
`code/` and `results/` are empty and there is nothing to run there.

What can be re-checked mechanically is the Lean 4 development in `lean/`, which formalizes the
paper's exact finite steps. From that directory, with the toolchain `lean/lean-toolchain` names:

```text
lake exe cache get
lake env lean ShortVariableMobius.lean
```

The first command fetches mathlib's prebuilt library; the second elaborates the file and prints
the `#print axioms` block that `lean/axioms.txt` stores. A clean run reports no errors and no
warnings, and every theorem's axiom set is mathlib's standard `[propext, Classical.choice,
Quot.sound]`; a `sorryAx` on any line would mean that statement is not proved.

`lean/README.md` lists exactly which of the paper's statements are formalized. The analytic
estimates -- Davenport's exponential-sum bound, the Cauchy-Schwarz and Parseval step for the
bilinear factor, the divisor moments, the explicit bound on the largest divisor function, and
the asymptotic estimates that assemble the main theorem -- are not formalized, and neither is
the main theorem.

## Provenance

Artifact-level source paths and SHA-256 digests are in `PACKET.json`.

## Formalization

Status: `library-in-packet`. A standalone Lean 4 / mathlib v4.33.1 project (lean/): 14 theorems on the exact finite steps of the proof -- the factorisation of mu(n1 n2 n3), the decomposition S = sum_{e<=N} mu(e) W(e) over the common divisor of n1 and n2, the counting bound on |W(e)|, the two bounds on X_N, the three-range split of the tail with the two finite-sum inequalities it rests on, the canonical splitting of a squarefree modulus across a product, and the two-residue-class count behind the large-divisor lemma. No sorry; axioms in lean/axioms.txt. The analytic estimates and the main theorem are not formalized, and lean/README.md states the boundary statement by statement.

The presence of a Lean directory does not mean the paper's main theorem is
machine-checked; consult the manuscript appendix and the files themselves for
the formalization boundary.

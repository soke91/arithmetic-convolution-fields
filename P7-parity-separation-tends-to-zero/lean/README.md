# kappa-zero

Lean 4 (v4.33.1) and mathlib (v4.33.1) formalization of the exact, finite, combinatorial steps of the proof
that the parity-separation constant tends to zero. One file, `KappaZero.lean`; `axioms.txt` is the output of
`#print axioms` for every theorem (only `propext`, `Classical.choice`, `Quot.sound`; no `sorryAx`).

Build: `lake exe cache get && lake build` (downloads mathlib's prebuilt cache).

Formalized: the one-large-prime identity for the band Moebius sum (every N with N < (Q+1)^2), exactness of
the small-prime-part construction on smooth rows (identity (E), with the dead-class term), the Bonferroni
identity and its truncation error, the parity split and the ratio bound, and the cone-distance step.
Not formalized: every analytic estimate (prime number theorem, Mertens, Rankin, Buchstab) and the main
theorem. The presence of this directory does not mean the paper's main theorem is machine-checked.

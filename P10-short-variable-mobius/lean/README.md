# short-variable-mobius

Lean 4 (v4.33.1) and mathlib (v4.33.1) formalization of the exact, finite steps of the paper
*A ternary Möbius sum with one short variable*. One file, `ShortVariableMobius.lean`;
`axioms.txt` is the output of `#print axioms` for every theorem (only `propext`,
`Classical.choice`, `Quot.sound`; no `sorryAx`).

Build: `lake exe cache get && lake build` (downloads mathlib's prebuilt cache). To elaborate the
file without building a library, `lake env lean ShortVariableMobius.lean`.

Formalized, by the paper's numbering:

| statement | theorem |
|---|---|
| the factorisation (11) of `μ(n₁n₂n₃)` | `moebius_mul_mul` |
| Lemma 7, the decomposition `S_{u,v}(H,N) = Σ_{e≤N} μ(e)W(e)` over the common divisor, with the definition (10) of `W(e)` | `S_eq_sum_moebius_mul_W` |
| the count of triples behind (12), and (12) itself, `\|W(e)\| ≤ (⌊H/e⌋+1)⌊N/e⌋` | `card_filter_triples_le`, `norm_W_le` |
| the two bounds (13), `H ≤ X_N ≤ H τ*(N)` | `le_XN`, `XN_le_mul_taustar` |
| the support statement behind the third range of (14) | `W_eq_zero_of_forall_not_dvd` |
| the two finite-sum inequalities (14) rests on, `Σ_{E<e≤M} 1/e² ≤ 1/E` and `Σ_{H<e≤M} 1/e ≤ log(M/H)` | `sum_inv_sq_le`, `sum_inv_le_log_div` |
| the three-range split (14), `Σ_{E<e≤N} \|W(e)\| ≤ 2NH/E + N log(X_N/H) + N` | `sum_norm_W_le` |
| the canonical divisor splitting (15) of Lemma 10, as the statement that exactly one admissible factorisation of a squarefree modulus exists | `canonical_split` |
| the two-residue-class count behind Lemma 11 | `card_le_of_two_classes` |

Not formalized: every analytic estimate — Davenport's bound for the exponential sum of the
Möbius function, the Cauchy–Schwarz and Parseval step for the bilinear factor, the divisor
moments, the explicit bound on the largest divisor function, and the asymptotic estimates that
assemble the main theorem — and the main theorem itself. The presence of this directory does not
mean the paper's main theorem is machine-checked; the file's header and the table above are
where its boundary is recorded.

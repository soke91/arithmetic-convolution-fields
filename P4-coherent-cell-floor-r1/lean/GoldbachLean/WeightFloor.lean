/-
  COMPILED, and the three load-bearing theorems are clean:

      weight_ne_sqrt_of_rat               [propext, Classical.choice, Quot.sound]
      rat_of_weights                      [propext, Classical.choice, Quot.sound]
      weightfloor_hypothesis              [propext, Classical.choice, Quot.sound]
      no_prime_root_of_weight_equals_sqrt [propext, sorryAx]     <- deliberate

  Elaborated with `lake env lean` on toolchain `leanprover/lean4:v4.33.1`,
  exit 0, one warning: *declaration uses `sorry`* at the corollary below, which
  is the one that should show it.  Nothing was changed to make it compile --
  no tactic, no statement, no lemma name.

  THE HEADER THIS REPLACES WAS FALSE, AND THE ERROR IS WORTH KEEPING

  It read: "UNCOMPILED.  There is no Lean toolchain in this checkout: `lake`,
  `lean` and `elan` are all absent from PATH and `lean/.lake` does not exist."
  The first clause was true and the conclusion did not follow.  `elan` was
  installed at `C:\Users\AHS\.elan\bin` the whole time, with the toolchain
  `leanprover--lean4---v4.33.1` -- exactly the tag `lean/lean-toolchain` asks
  for.  I ran `command -v elan`, got nothing, and wrote down that there was no
  toolchain.  NOT ON `PATH` AND NOT INSTALLED ARE DIFFERENT FACTS, and I never
  checked the install directory.

  That was the third time in one day I reported an absence that was really a
  search that could not see, and then built a conclusion on it.  The other two
  were a `grep` of a Korean document for English text, and a `grep` for
  `O(1/\ell)` whose escaping was mangled through a shell string, returning one
  hit on a file holding sixteen.  `rules/METHOD.md` records the general form
  under Literature: a claim of absence is only as good as the search boundary,
  and a zero has to be told apart from a broken query before it means anything.

  The blocker was never a missing toolchain.  It was one environment decision,
  and it belonged to the user, who made it.

  WHAT THIS FILE IS AND IS NOT EVIDENCE FOR

  `lean/GoldbachLean/Statements.lean` keeps `sorry` deliberately so that
  `sorryAx` prints and no one can mistake an unproved line for a verified one.
  The corollary at the end of this file keeps its `sorry` for the same reason:
  it is a check on the manuscript's table, never a dependency of the main
  lemma, and half-proving it would be worse than leaving it visible.

  WHAT IT IS FOR
  --------------
  `rules:0097` records that `lem:weightfloor`'s only hypothesis is

      w(p) != sqrt p   for every prime p not dividing N,

  and that the manuscript discharges it with a table of four algebraic
  equations, three computed roots (2.618034, 3.246980, 4.000000), and a sweep
  of odd primes below 10^6.  The record then observes that one line replaces
  all of it: the four weights are RATIONAL at integer p, and sqrt p is
  IRRATIONAL at prime p, so they are never equal.

  That is a statement Lean can verify completely, and doing so would replace
  numerical evidence by a proof -- which is the one exchange in this repository
  that strictly gains.  It also generalises: the argument does not know how
  many weights there are, only that they are rational at integer points.

  Mathlib supplies the arithmetic input as `Nat.Prime.irrational_sqrt`.
-/

import Mathlib

/-!
# The weight floor hypothesis

`rules:0097`'s `lem:weightfloor` needs `w(p) ≠ √p` at every prime `p`.  The
manuscript discharges it with a table of four algebraic equations, three
computed roots and a sweep of odd primes below `10^6`.  One line replaces all
of it: the weights are rational at integer arguments and `√p` is irrational at
a prime, so they are never equal.

`weight_ne_sqrt_of_rat` is that line, stated for any family rational at natural
arguments rather than for the manuscript's four.  `weightfloor_hypothesis`
instantiates it.  The corollary at the end keeps a `sorry` on purpose: it
checks the manuscript's table and is not a dependency of anything above it.

**Irrationality buys the hypothesis and not the margin.**  The sweep still buys
`|log w(p)/log p − 1/2| = 0.130930` at `p = 3`, and leaves the proof of the
hypothesis only.
-/

namespace WeightFloor

open Real

/-- The general lemma.  A weight that takes rational values at natural
arguments never equals `sqrt p` at a prime `p`.  This is the whole content of
`rules:0097`'s replacement for the table, the three roots and the sweep, and it
covers every weight rational at integer points rather than the four in the
manuscript. -/
theorem weight_ne_sqrt_of_rat
    (w : ℕ → ℝ) (hw : ∀ n : ℕ, ∃ q : ℚ, w n = (q : ℝ))
    {p : ℕ} (hp : p.Prime) : w p ≠ Real.sqrt p := by
  intro h
  obtain ⟨q, hq⟩ := hw p
  exact hp.irrational_sqrt ⟨q, hq.symm.trans h⟩

/-- The four weights of the manuscript, as one family.  Each is a ratio of
integers at integer `p`, so `rat_of_weights` below is the only thing the main
lemma needs about them.  `p / (p - 1)`, `(p - 1) / (p - 2)`, `p / (p - 2)` and
`1`, written with natural subtraction guarded by the hypothesis `3 <= p` that
`rules:0097` carries (`p` odd and not dividing `N`). -/
noncomputable def weights : Fin 4 → ℕ → ℝ
  | 0, p => (p : ℝ) / ((p : ℝ) - 1)
  | 1, p => ((p : ℝ) - 1) / ((p : ℝ) - 2)
  | 2, p => (p : ℝ) / ((p : ℝ) - 2)
  | 3, _ => 1

/-- Each of the four is rational at every natural argument.  Division by zero
is `0` in `ℝ`, which is rational, so no side condition on `p` is needed here --
the guard `3 <= p` matters for the weight's meaning, not for its rationality. -/
theorem rat_of_weights (i : Fin 4) (n : ℕ) : ∃ q : ℚ, weights i n = (q : ℝ) := by
  fin_cases i
  · exact ⟨(n : ℚ) / ((n : ℚ) - 1), by push_cast [weights]; ring_nf⟩
  · exact ⟨((n : ℚ) - 1) / ((n : ℚ) - 2), by push_cast [weights]; ring_nf⟩
  · exact ⟨(n : ℚ) / ((n : ℚ) - 2), by push_cast [weights]; ring_nf⟩
  · exact ⟨1, by simp [weights]⟩

/-- `lem:weightfloor`'s hypothesis, for the manuscript's four weights, with no
table, no computed root and no sweep. -/
theorem weightfloor_hypothesis
    (i : Fin 4) {p : ℕ} (hp : p.Prime) : weights i p ≠ Real.sqrt p :=
  weight_ne_sqrt_of_rat (weights i) (rat_of_weights i) hp

/-!
### What is deliberately NOT here

**The margin is a different statement and it is not replaced.**  `rules:0097`
is explicit: irrationality buys the *hypothesis*; the sweep buys the *margin*
`|log w(p) / log p - 1/2|`, printed as `0.130930` at `p = 3`, and irrationality
does not imply a margin.  Formalising the hypothesis therefore removes the
sweep from the proof of `lem:weightfloor` and from nothing else.

**The rational-root corollary is stated, not proved.**  `rules:0097` also
records exact roots: `p/(p-1) = sqrt p` gives `t^2 - t - 1 = 0` at
`p = (3+sqrt 5)/2`; `(p-1)/(p-2) = sqrt p` gives `t^3 - t^2 - 2t + 1 = 0` at
`p = 4 cos^2(pi/7)`; `p/(p-2) = sqrt p` gives `t^2 - t - 2 = 0` at `p = 4`.
Substituting `p = m^2` turns "is the root prime" into "has this integer
polynomial an integer root", which the rational root theorem settles:
`m^2 - m - 1` and `m^3 - m^2 - 2m + 1` have none, and `m^2 - m - 2` has only
`m = 2`, i.e. `p = 4`, not prime.

That corollary is worth formalising **only** as a check on the manuscript's
table; the lemma above already discharges the hypothesis without it, and the
manuscript's three printed decimals should lean on the rational root theorem
rather than on those decimals either way.  Written as a statement and left
unproved rather than half-proved:
-/

/-- The rational-root check behind the manuscript's "no root is prime" row.
NOT PROVED and not needed by `weightfloor_hypothesis`; recorded so that the
table's claim has a Lean-shaped statement to be checked against. -/
theorem no_prime_root_of_weight_equals_sqrt :
    (∀ m : ℤ, m ^ 2 - m - 1 ≠ 0)
    ∧ (∀ m : ℤ, m ^ 3 - m ^ 2 - 2 * m + 1 ≠ 0)
    ∧ (∀ m : ℤ, m ^ 2 - m - 2 = 0 → (m = 2 ∨ m = -1) ∧ ¬ Nat.Prime (m ^ 2).natAbs) := by
  sorry

end WeightFloor

#print axioms WeightFloor.weight_ne_sqrt_of_rat
#print axioms WeightFloor.rat_of_weights
#print axioms WeightFloor.weightfloor_hypothesis
#print axioms WeightFloor.no_prime_root_of_weight_equals_sqrt

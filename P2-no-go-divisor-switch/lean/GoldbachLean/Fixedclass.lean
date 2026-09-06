/-
  The finite rearrangements of P1 -- the steps of Theorem A that hold
  before any arithmetic is used.

  P1 runs a chain: subtract the mean term, switch divisors, split the
  divisor sum into a complete sum minus a tail, and factor `u = mk`.
  Those steps are **finite rearrangements**: they are true for arbitrary
  `Λ`, `μ` and `φ`, and the only thing they use about the index sets is
  which pairs belong to them.  They are what is here.

  `Λ`, `μ`, `φ` and the weight `w` are arbitrary functions throughout, so
  what is verified is that these steps use no arithmetic property of
  them.  That is the point: P1's content is that the two functionals
  `E_3` and `E_4` differ in exactly one factor, and a step that secretly
  used a property of `μ` could not support a claim of that shape.

  **What is absent matters.**  The fourth step, `lem:complete`
  (`∑_{k ∣ u, (k,N)=1} μ(k) = 1` iff `rad(u) ∣ N`), is genuinely
  arithmetic -- it is the multiplicativity of `μ` -- and it is not here.
  Neither is Theorem A itself, which rests on Bombieri--Vinogradov, and
  Mathlib does not carry that.  What is formalised is the scaffolding,
  not the theorem.
-/

import Mathlib

namespace Fixedclass

open Finset

variable {ι κ : Type*} [DecidableEq ι] [DecidableEq κ]

/-!
## 1. The mean term leaves the `k`-sum

`eq:split0` reads `T_1(t) = D(t) − C(t) B(K)`.  What licenses pulling
`C(t)` out of the sum over `k` is that it does not depend on `k` -- and
nothing else.  In particular no property of `μ`, of the weight, or of
`φ` is used, which is why the same line serves both `E_3` and `E_4`.
-/

/-- `∑_k μ(k) w(k) [A(k) − C/φ(k)] = ∑_k μ(k) w(k) A(k) − C ∑_k μ(k) w(k)/φ(k)`.

This is `eq:split0` with everything arbitrary.  `C` is a scalar, and
that is the whole hypothesis. -/
theorem meanterm_factors
    (K : Finset ι) (mu w phi A : ι → ℝ) (C : ℝ) :
    ∑ k ∈ K, mu k * w k * (A k - C / phi k)
      = (∑ k ∈ K, mu k * w k * A k)
        - C * ∑ k ∈ K, mu k * w k / phi k := by
  rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  ring

/-- **Guard.**  If the subtracted term depended on `k`, it would not come
out of the sum.

The step above is not an algebraic triviality that would survive any
rewriting of the mean term: it is exactly the `k`-independence of `C(t)`
that is doing the work.  A reader who takes `eq:split0` for "the mean
term can always be separated" is contradicted by this. -/
theorem meanterm_does_not_factor_when_k_dependent :
    ∃ (K : Finset ℕ) (mu w phi A C : ℕ → ℝ),
      (∑ k ∈ K, mu k * w k * (A k - C k / phi k))
        ≠ (∑ k ∈ K, mu k * w k * A k)
          - C 0 * ∑ k ∈ K, mu k * w k / phi k := by
  refine ⟨{0, 1}, fun _ => 1, fun _ => 1, fun _ => 1, fun _ => 0,
          fun k => (k : ℝ), ?_⟩
  norm_num

/-!
## 2. The divisor switch

`eq:switch` turns `∑_{k<K} μ(k) ∑_{n ≡ N (k)} Λ(n) μ(N−n)` into
`∑_u Λ(N−u) μ(u) σ_K(u)`.  The paper calls it "a finite rearrangement,
hence exact".  It is two moves, and they are separated here because they
are separately true and separately checkable: the **exchange** of the
double sum, and the **reindexing** `n ↦ N − n`.

The condition `n ≡ N (mod k)` becomes `k ∣ u` under the reindexing.
Neither move knows that; the exchange holds for an arbitrary relation.
-/

/-- The exchange.  Summing over `k` and then over the `n` it selects is
the same as summing over `n` and then over the `k` that select it.

`R` is arbitrary.  In the paper `R k n` is `n ≡ N (mod k)`, and the
inner factor `μ(k)` sits outside the inner sum because it does not
depend on `n` -- the same reason as in `meanterm_factors`. -/
theorem switch_exchange
    (K S : Finset ι) (mu : ι → ℝ) (f : ι → ℝ)
    (R : ι → ι → Prop) [∀ k n, Decidable (R k n)] :
    ∑ k ∈ K, mu k * ∑ n ∈ S.filter (fun n => R k n), f n
      = ∑ n ∈ S, f n * ∑ k ∈ K.filter (fun k => R k n), mu k := by
  simp only [Finset.mul_sum, Finset.sum_mul, Finset.sum_filter]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun n _ =>
    Finset.sum_congr rfl fun k _ => by
      by_cases h : R k n <;> simp [h, mul_comm]

/-- The reindexing.  A sum over `n` of `g (N − n)` is a sum over `u` of
`g u`, provided `n ↦ N − n` carries one index set onto the other and
back.

Stated with both directions as hypotheses rather than derived from
`n ≤ N`, because that is what the paper's index sets supply:
`n ≤ t < N` on one side and `N − t ≤ u < N` on the other. -/
theorem switch_reindex
    (N : ℕ) (S T : Finset ℕ) (g : ℕ → ℝ)
    (hi : ∀ n ∈ S, N - n ∈ T) (hj : ∀ u ∈ T, N - u ∈ S)
    (hl : ∀ n ∈ S, N - (N - n) = n) (hr : ∀ u ∈ T, N - (N - u) = u) :
    ∑ n ∈ S, g (N - n) = ∑ u ∈ T, g u :=
  Finset.sum_nbij' (fun n => N - n) (fun u => N - u) hi hj hl hr
    (fun _ _ => rfl)

/-- **`eq:switch` itself**, as the composition of the two moves above.

The two lemmas are separately true, and a file that stopped there would
not have said that they compose into the paper's identity -- pieces that
each hold can still fail to fit.  So the composition is stated.

`R k n` is the paper's `n ≡ N (mod k)`, which under the reindexing is
`k ∣ u`; here it is passed in as `Dvd k u` and the hypothesis `hR` is
exactly that translation.  Making it a hypothesis rather than proving it
is the honest form: it is the one arithmetic fact this step uses, and
naming it keeps it from hiding inside the rearrangement. -/
theorem switch_identity
    (N : ℕ) (K S T : Finset ℕ) (mu f : ℕ → ℝ)
    (R : ℕ → ℕ → Prop) [∀ k n, Decidable (R k n)]
    (Dvd : ℕ → ℕ → Prop) [∀ k u, Decidable (Dvd k u)]
    (hi : ∀ n ∈ S, N - n ∈ T) (hj : ∀ u ∈ T, N - u ∈ S)
    (hl : ∀ n ∈ S, N - (N - n) = n) (hr : ∀ u ∈ T, N - (N - u) = u)
    (hR : ∀ k, ∀ n ∈ S, R k n ↔ Dvd k (N - n)) :
    ∑ k ∈ K, mu k * ∑ n ∈ S.filter (fun n => R k n), f n
      = ∑ u ∈ T, f (N - u) * ∑ k ∈ K.filter (fun k => Dvd k u), mu k := by
  rw [switch_exchange K S mu f R]
  refine Finset.sum_nbij' (fun n => N - n) (fun u => N - u) hi hj hl hr ?_
  intro n hn
  refine congrArg₂ (· * ·) (by rw [hl n hn]) ?_
  refine Finset.sum_congr (Finset.filter_congr fun k _ => ?_) fun _ _ => rfl
  simpa using hR k n hn

/-- **Guard.**  The bijection alone does not give the switch: without the
translation `hR`, the identity is false.

Take `N = 2`, `S = {0}`, `T = {2}`.  The map `n ↦ N − n` carries each
onto the other and back, so `hi`, `hj`, `hl` and `hr` all hold.  With
`R` always true and `Dvd` always false, the left side is `1` and the
right side is `0`.

So `hR` is not bookkeeping that could be dropped for brevity: it carries
the whole step.  In the paper it is the line "n ≡ N (mod k) with n ≤ t
is the same as u := N − n being a multiple of k", and a reading that
treats `eq:switch` as reindexing alone is contradicted here. -/
theorem switch_needs_the_translation :
    ∑ k ∈ ({1} : Finset ℕ), (1 : ℝ) *
        ∑ n ∈ ({0} : Finset ℕ).filter (fun _ => True), (1 : ℝ)
      ≠ ∑ u ∈ ({2} : Finset ℕ), (1 : ℝ) *
          ∑ k ∈ ({1} : Finset ℕ).filter (fun _ => False), (1 : ℝ) := by
  norm_num

/-- And the bijection hypotheses of that counterexample do hold, so it is
the translation and not the reindexing that fails there. -/
theorem switch_guard_bijection_holds :
    (∀ n ∈ ({0} : Finset ℕ), 2 - n ∈ ({2} : Finset ℕ))
    ∧ (∀ u ∈ ({2} : Finset ℕ), 2 - u ∈ ({0} : Finset ℕ))
    ∧ (∀ n ∈ ({0} : Finset ℕ), 2 - (2 - n) = n)
    ∧ (∀ u ∈ ({2} : Finset ℕ), 2 - (2 - u) = u) := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> intro x hx <;>
    simp only [Finset.mem_singleton] at hx <;> subst hx <;> decide

/-!
## 3. The complete sum minus the tail

`eq:PR` writes `σ_K(u)` -- the sum over `k ∣ u` with `k < K` -- as the
complete divisor sum minus the part with `k ≥ K`, giving `D = P − R`.
Splitting a finite sum at a predicate is all this is; the paper's work
is in bounding the two pieces, not in the split.
-/

/-- `∑_{k ∈ A, p k} f = ∑_{k ∈ A} f − ∑_{k ∈ A, ¬ p k} f`. -/
theorem complete_minus_tail
    (A : Finset ι) (f : ι → ℝ) (p : ι → Prop) [DecidablePred p] :
    ∑ k ∈ A.filter p, f k
      = (∑ k ∈ A, f k) - ∑ k ∈ A.filter (fun k => ¬ p k), f k := by
  rw [eq_sub_iff_add_eq, Finset.sum_filter_add_sum_filter_not]

/-!
## 4. `u = mk`

`eq:R1` replaces the sum over `u` and its divisors `k` by a sum over
pairs.  This is divisor bookkeeping and knows nothing about `Λ` or `μ`;
the arithmetic the paper adds on top -- that squarefreeness forces
`(m,k) = 1` and `μ(u)μ(k) = μ(m)μ²(k)` -- is a separate step and is not
here.
-/

/-- Summing a function of `(u/k, k)` over the divisors of `u` is summing
over the factorisations of `u`. -/
theorem divisor_pairs (u : ℕ) (f : ℕ → ℕ → ℝ) :
    ∑ k ∈ u.divisors, f (u / k) k
      = ∑ p ∈ u.divisorsAntidiagonal, f p.1 p.2 := by
  rw [← Nat.map_div_left_divisors, Finset.sum_map]
  rfl

/-!
## 5. The two functionals differ in one factor

P1's claim is that `E_3` and `E_4` are one expression with two weights,
and that this single difference decides everything.  With the weight
arbitrary there is nothing left to prove -- which is the point.  What is
shown is that the two are instances of the same term.
-/

/-- `T_w(t) = ∑_{k} μ(k) w(k) E_μ(t;k)`. -/
def T (K : Finset ι) (mu w E : ι → ℝ) : ℝ := ∑ k ∈ K, mu k * w k * E k

/-- `E_4` is `T` at `w = 1`. -/
theorem T_flat (K : Finset ι) (mu E : ι → ℝ) :
    T K mu (fun _ => 1) E = ∑ k ∈ K, mu k * E k := by
  simp [T]

/-- `E_3` is `T` at `w = log`. -/
theorem T_log (K : Finset ι) (mu lg E : ι → ℝ) :
    T K mu lg E = ∑ k ∈ K, mu k * lg k * E k := rfl

/-- Changing the weight changes one factor inside the sum and nothing
else -- no index set moves, no other factor moves. -/
theorem T_same_shape
    (K : Finset ι) (mu E : ι → ℝ) (w w' : ι → ℝ)
    (h : ∀ k ∈ K, w k = w' k) :
    T K mu w E = T K mu w' E :=
  Finset.sum_congr rfl fun k hk => by rw [h k hk]

end Fixedclass

#print axioms Fixedclass.meanterm_factors
#print axioms Fixedclass.meanterm_does_not_factor_when_k_dependent
#print axioms Fixedclass.switch_exchange
#print axioms Fixedclass.switch_reindex
#print axioms Fixedclass.complete_minus_tail
#print axioms Fixedclass.divisor_pairs
#print axioms Fixedclass.T_flat
#print axioms Fixedclass.T_log
#print axioms Fixedclass.T_same_shape
#print axioms Fixedclass.switch_identity
#print axioms Fixedclass.switch_needs_the_translation
#print axioms Fixedclass.switch_guard_bijection_holds

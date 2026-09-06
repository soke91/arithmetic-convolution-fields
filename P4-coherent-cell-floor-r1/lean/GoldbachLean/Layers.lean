/-
  Two finite rearrangements from P2 -- `prop:layers` and `prop:flatsum`.

  These sit where `Cellmom` sits for P4: strip the arithmetic and what
  is left is a reindexing of a finite sum, and that is all this file
  proves.  `Λ`, `μ` and `log` are arbitrary functions here, so nothing
  below uses any arithmetic property of them -- and that it does not is
  part of the statement.

  **What is absent matters.**  P2's conclusion -- that no weight
  Bombieri--Vinogradov reaches extracts the sum -- rests on
  Bombieri--Vinogradov, which Mathlib does not have.  What is formalised
  is the identities that lay out the design space, not the no-go.
-/

import Mathlib

namespace Layers

open Finset

variable {ι κ : Type*} [DecidableEq ι] [DecidableEq κ]

/-- **The skeleton of `prop:layers`**: cutting a double sum along its
  inner variable.

  `A` is the outer index (the paper's `k`) and the filter picks the
  inner one.  The paper cuts along a condition it calls symmetric --
  `(k, m) = 1` -- and that symmetry is what licenses the exchange.  Here
  the condition is an arbitrary `P`, so the exchange holds for any
  condition whatever; symmetry is not needed for the identity, only for
  the two sides to be described by the same predicate.

  This is an instance of `Finset.sum_comm'`, but *which* instance is the
  paper's content, so it is stated in that shape. -/
theorem layers_exchange
    (A : Finset ι) (B : Finset κ) (f : ι → κ → ℝ)
    (P : ι → κ → Prop) [∀ i j, Decidable (P i j)] :
    ∑ i ∈ A, ∑ j ∈ B.filter (fun j => P i j), f i j
      = ∑ j ∈ B, ∑ i ∈ A.filter (fun i => P i j), f i j := by
  simp only [sum_filter]
  rw [Finset.sum_comm]

/-- If the summand is nonnegative then so is every layer after the cut.

  The paper gets this from `log k ≥ 0` and `Λ ≥ 0`.  Here only the
  hypothesis is taken; which functions they are is not used. -/
theorem layer_nonneg
    (A : Finset ι) (f : ι → κ → ℝ)
    (P : ι → κ → Prop) [∀ i j, Decidable (P i j)]
    (hf : ∀ i j, 0 ≤ f i j) (j : κ) :
    0 ≤ ∑ i ∈ A.filter (fun i => P i j), f i j :=
  Finset.sum_nonneg fun i _ => hf i j

/-- **The skeleton of `prop:flatsum`.**

  `T w = (Σ_k μ²(k) · w k · H k) − B_w · C` is a definition, and what
  the paper asserts is that the same definition serves every weight.
  Since `w` is arbitrary there is no arithmetic to prove; what remains
  is that the two weights the reduction actually uses are instances of
  one expression.  That is the precise content of "both ends are the
  same sum". -/
def T (A : Finset ι) (musq w H : ι → ℝ) (Bw C : ℝ) : ℝ :=
  (∑ k ∈ A, musq k * w k * H k) - Bw * C

theorem flatsum_flat (A : Finset ι) (musq H : ι → ℝ) (B1 C : ℝ) :
    T A musq (fun _ => 1) H B1 C
      = (∑ k ∈ A, musq k * H k) - B1 * C := by
  simp [T]

theorem flatsum_log (A : Finset ι) (musq lg H : ι → ℝ) (Blog C : ℝ) :
    T A musq lg H Blog C
      = (∑ k ∈ A, musq k * lg k * H k) - Blog * C := rfl

/-- Changing the weight changes one factor inside the sum and moves
  nothing else. -/
theorem flatsum_same_shape
    (A : Finset ι) (musq H : ι → ℝ) (w w' : ι → ℝ) (Bw Bw' C : ℝ)
    (h : ∀ k ∈ A, w k = w' k) (hB : Bw = Bw') :
    T A musq w H Bw C = T A musq w' H Bw' C := by
  unfold T
  rw [hB]
  congr 1
  exact Finset.sum_congr rfl fun k hk => by rw [h k hk]

/-- **The skeleton of `prop:dilate`.**

  The paper writes `A(N;k) = μ(k) · H(N;k)`: the outer Möbius factor is
  constant across the inner sum, so it comes out of it.  Nothing is
  assumed about `c` here -- in particular not that it is a sign. -/
theorem dilate_pull (A : Finset ι) (c : ℝ) (g : ι → ℝ) :
    ∑ i ∈ A, c * g i = c * ∑ i ∈ A, g i :=
  (Finset.mul_sum A g c).symm

/-- **The skeleton of `prop:posweights`.**

  Substituting the dilate into `T_log` uses `μ(k)² = 1` on squarefree
  `k`, which turns a signed weight into a nonnegative one.  The content
  is that a factor equal to `1` on the index set may be dropped; that
  the factor is `μ²` is arithmetic and is not used. -/
theorem posweights_unit (A : Finset ι) (u f : ι → ℝ)
    (hu : ∀ i ∈ A, u i = 1) :
    ∑ i ∈ A, u i * f i = ∑ i ∈ A, f i :=
  Finset.sum_congr rfl fun i hi => by rw [hu i hi, one_mul]

/-- The same with the sign written out, which is the form the paper
  uses: `s k * s k = 1` on the index set. -/
theorem posweights_sign (A : Finset ι) (s f : ι → ℝ)
    (hs : ∀ i ∈ A, s i * s i = 1) :
    ∑ i ∈ A, (s i * s i) * f i = ∑ i ∈ A, f i :=
  posweights_unit A (fun i => s i * s i) f hs

/-- **The skeleton of `prop:untrunc`.**

  Dropping the truncation and the coprimality restriction replaces one
  index set by another, and the paper then identifies the two sums by a
  correspondence between the sets.  The identification, not any
  property of the summand, is the step -- so it is stated for an
  arbitrary correspondence. -/
theorem untrunc_reindex (A : Finset ι) (B : Finset κ) (e : ι → κ)
    (hmaps : ∀ i ∈ A, e i ∈ B)
    (hinj : ∀ i ∈ A, ∀ j ∈ A, e i = e j → i = j)
    (hsurj : Set.SurjOn e ↑A ↑B)
    (f : κ → ℝ) :
    ∑ i ∈ A, f (e i) = ∑ j ∈ B, f j :=
  Finset.sum_nbij e hmaps hinj hsurj (fun _ _ => rfl)

/-- **The skeleton of `prop:combined`.**

  Expanding the squarefree indicator and regrouping by `k = d²j` turns
  a sum of sums into a sum over pairs.  What survives once the
  arithmetic is stripped is that the flattening is exact in both
  directions. -/
theorem combined_flatten (A : Finset ι) (B : ι → Finset κ) (f : ι → κ → ℝ) :
    ∑ i ∈ A, ∑ j ∈ B i, f i j
      = ∑ p ∈ A.sigma B, f p.1 p.2 :=
  (Finset.sum_sigma A B (fun p => f p.1 p.2)).symm

/-- The regrouping loses nothing: each pair on one side is counted once
  on the other. -/
theorem combined_card (A : Finset ι) (B : ι → Finset κ) :
    (A.sigma B).card = ∑ i ∈ A, (B i).card :=
  Finset.card_sigma A B

end Layers

#print axioms Layers.layers_exchange
#print axioms Layers.layer_nonneg
#print axioms Layers.flatsum_flat
#print axioms Layers.flatsum_log
#print axioms Layers.flatsum_same_shape
#print axioms Layers.dilate_pull
#print axioms Layers.posweights_unit
#print axioms Layers.posweights_sign
#print axioms Layers.untrunc_reindex
#print axioms Layers.combined_flatten
#print axioms Layers.combined_card

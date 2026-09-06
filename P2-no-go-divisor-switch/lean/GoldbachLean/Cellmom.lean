/-
  P4, Lemma `lem:cellmom` — the algebraic heart.

  WHAT THIS FILE CLAIMS, AND WHAT IT DOES NOT

  The paper's lemma has two steps.  Only the second is here.

    (1)  Under the sign null, `Var(∑ᵥ εᵥ hᵥ) = ∑ᵥ wᵥ hᵥ²`.
         This is the probabilistic step: `E[εᵥ ε_w] = δ`.  It is NOT
         formalised in this file.
    (2)  With `hᵥ = u_c(v)/n_c − u_a(v)/n`, that sum is exactly
         `Q_cc/n_c² − 2 Q_ca/(n_c n) + Q_aa/n²`.
         This is pure algebra over a finite index set, and it is what
         is proved below.

  Step (2) is where the paper's content sits — it is the statement that
  the variance has **three** terms and not two, which an earlier
  printing of this programme got wrong (see `note:threeterms`).  Step
  (1) is standard and was never in doubt.

  Saying which half is done is the point.  A formalisation that lets a
  reader believe the whole lemma is machine-checked, when half of it is
  not, is the same defect this repository spent 2026-09-05 removing
  from its prose.

  `w` is `μ²` on the surviving support, `u_c` and `u_a` are the paper's
  `u_c(v) = ∑_{N ∈ c} Λ(N−v)/√V(N)` for the cell and for the whole
  band.  Nothing below needs them to be *those* functions: the identity
  holds for any two real vectors, which is why the paper can say the
  proof uses nothing beyond `c ⊆ a`.
-/

import Mathlib

open Finset

variable {ι : Type*} [Fintype ι]

namespace Cellmom

/-- The paper's `Q_{cd} = ∑_v μ²(v) u_c(v) u_d(v)`. -/
def Q (w u u' : ι → ℝ) : ℝ := ∑ v, w v * (u v * u' v)

/-- **`eq:cellmom`, the algebraic half.**

`∑_v w_v (u_c(v)/n_c − u_a(v)/n)²` is the three-term expression, and
the middle term carries the factor two.  Dropping it — reading the
variance as a difference of two means — is the error `note:threeterms`
records; here the three terms are forced by `ring`. -/
-- PAPER: P4 lem:cellmom  (the algebraic half only)
theorem sum_sq_eq_three_terms (w uc ua : ι → Real) (nc n : Real)
    (hnc : nc ≠ 0) (hn : n ≠ 0) :
    (Finset.univ.sum fun v => w v * (uc v / nc - ua v / n) ^ 2)
      = Q w uc uc / nc ^ 2
        - (2 * Q w uc ua) / (nc * n)
        + Q w ua ua / n ^ 2 := by
  unfold Q
  rw [Finset.mul_sum, Finset.sum_div, Finset.sum_div, Finset.sum_div,
      ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun v _ => ?_
  field_simp
  ring

/-- The middle term is not optional.  Dropping it -- reading the
variance as a difference of two means -- leaves exactly this gap, which
is nonzero whenever `Q w uc ua` is.  `note:threeterms` records that an
earlier printing of this programme made that mistake; here the size of
it is a machine-checked statement rather than a remark. -/
-- PAPER: P4 note:threeterms  (the error of the two-term form)
theorem two_term_form_is_off_by (w uc ua : ι → Real) (nc n : Real)
    (hnc : nc ≠ 0) (hn : n ≠ 0) :
    (Finset.univ.sum fun v => w v * (uc v / nc - ua v / n) ^ 2)
      - (Q w uc uc / nc ^ 2 + Q w ua ua / n ^ 2)
      = - ((2 * Q w uc ua) / (nc * n)) := by
  rw [sum_sq_eq_three_terms w uc ua nc n hnc hn]
  ring

end Cellmom

#print axioms Cellmom.sum_sq_eq_three_terms
#print axioms Cellmom.two_term_form_is_off_by

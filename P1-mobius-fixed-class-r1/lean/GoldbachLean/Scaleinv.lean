/-
  P4, Proposition `prop:scaleinv` — the three finite steps.

  WHAT IS HERE AND WHY THESE THREE

  The proposition's proof has an analytic half (an Euler product over
  the primes above 13, and its convergence) and a finite half (Rankin's
  truncation step, a count of pairs in an arithmetic progression, and
  the cancellation of the class count `m`).  Only the finite half is
  here.

  That is not because the other half is easy.  It is because **every
  error this repository has made in this proposition lived in the
  finite half**, and none in the analytic one:

    * `O(1/n_c)` written where the counting supplies only `O(1/ell)` —
      three times, and once it survived into a displayed equation that
      the proof itself withdrew two paragraphs later;
    * "two arithmetic progressions of length `n_c`" for what is `m`
      residue classes of length `ell`, with `m * ell = n_c`;
    * a percentage read off a bound that carries an unnamed constant
      and a divisor sum growing like `log B`.

  The Euler product has never been wrong here.  Formalising the part
  that breaks is worth more than machine-checking the part that does
  not, and saying which half is done is the point — a formalisation
  that lets a reader believe the whole proposition is checked would be
  the same defect this repository spent 2026-09-05 removing from its
  prose.

  THE THREE

    1. `rankin_head` / `rankin_tail` — the truncation step itself.
       Elementary: `d <= B` gives `1 <= (B/d)^s`, so a sum over
       `d <= B` is at most `B^s` times the same sum weighted by
       `d^{-s}`.  No Euler product enters; the product is what bounds
       the weighted sum afterwards, and that is the analytic half.

    2. `card_residue_range` and `card_residue_close` — the count of one
       residue class in `[0, ell)`, exactly and then to within `1`.
       Counting `i` for each `j` gives the constant `1`, where the
       tent-sum route gives `2`.  Both are `O(ell)`, and neither is
       `O(1/n_c)`.

    3. `blocks_cancel` — if each of `m^2` blocks of `ell^2` pairs
       carries an error `C * ell`, the relative error of the whole is
       `C / ell`.  **The `m` cancels identically**, which is the
       statement this proposition got wrong twice.
-/

import Mathlib

open Finset

namespace Scaleinv


/-! ### 1. Rankin's truncation step -/

/-- `d ≤ B` gives `∑ g ≤ B^s * ∑ g/d^s`.  This is the whole of Rankin's
trick as a finite statement; what makes the right-hand side useful is
the Euler product bounding it, which is not here. -/
theorem rankin_head (S : Finset ℕ) (g : ℕ → ℝ) (B s : ℝ)
    (hg : ∀ d ∈ S, 0 ≤ g d) (hs : 0 ≤ s) (hB : 0 < B)
    (hle : ∀ d ∈ S, 0 < (d : ℝ) ∧ (d : ℝ) ≤ B) :
    ∑ d ∈ S, g d ≤ B ^ s * ∑ d ∈ S, g d / (d : ℝ) ^ s := by
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum (fun d hd => ?_)
  obtain ⟨hd0, hdB⟩ := hle d hd
  have hratio : (1 : ℝ) ≤ B / d := (one_le_div hd0).mpr hdB
  have hone : (1 : ℝ) ≤ (B / d) ^ s := Real.one_le_rpow hratio hs
  calc g d = g d * 1 := by ring
    _ ≤ g d * (B / d) ^ s := mul_le_mul_of_nonneg_left hone (hg d hd)
    _ = B ^ s * (g d / (d : ℝ) ^ s) := by
        rw [Real.div_rpow hB.le hd0.le]
        field_simp

/-- `B ≤ d` gives `∑ g/d ≤ B^{-t} * ∑ g/d^{1-t}`.  The tail half.  The
two exponents stay two: with one exponent the two bounds meet only at
`1/2`, and `B^{-1/2}` is far weaker than `B^{η-1}`. -/
theorem rankin_tail (S : Finset ℕ) (g : ℕ → ℝ) (B t : ℝ)
    (hg : ∀ d ∈ S, 0 ≤ g d) (ht : 0 ≤ t) (hB : 0 < B)
    (hgt : ∀ d ∈ S, B ≤ (d : ℝ)) :
    ∑ d ∈ S, g d / (d : ℝ) ≤ B ^ (-t) * ∑ d ∈ S, g d / (d : ℝ) ^ (1 - t) := by
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum (fun d hd => ?_)
  have hdB := hgt d hd
  have hd0 : (0 : ℝ) < d := lt_of_lt_of_le hB hdB
  have hratio : (1 : ℝ) ≤ (d : ℝ) / B := (one_le_div hB).mpr hdB
  have hone : (1 : ℝ) ≤ ((d : ℝ) / B) ^ t := Real.one_le_rpow hratio ht
  have hgd := hg d hd
  have expand : B ^ (-t) * (g d / (d : ℝ) ^ (1 - t))
      = (g d / (d : ℝ)) * ((d : ℝ) / B) ^ t := by
    rw [Real.div_rpow hd0.le hB.le, Real.rpow_neg hB.le,
        Real.rpow_sub hd0, Real.rpow_one]
    field_simp
  rw [expand]
  calc g d / (d : ℝ) = (g d / (d : ℝ)) * 1 := by ring
    _ ≤ (g d / (d : ℝ)) * ((d : ℝ) / B) ^ t :=
        mul_le_mul_of_nonneg_left hone (by positivity)

/-! ### 2. The count of one residue class

The paper's per-class error.  Counting `i` for each `j` gives the
constant `1`; the tent-sum route gives `2`.  Both are `O(ell)`, and
neither is `O(1/n_c)` — which is the point, since `n_c = m * ell` and
the `m` is exactly what a previous printing gave itself for free. -/

/-- Exact count of one residue class in `[0, L)`. -/
theorem card_residue_range (d r L : ℕ) (hd : 0 < d) (hr : r < d) :
    ((Finset.range L).filter (fun i => i % d = r)).card
      = L / d + (if r < L % d then 1 else 0) := by
  induction L with
  | zero => simp
  | succ M ih =>
      have hlt : M % d < d := Nat.mod_lt _ hd
      have hsplit : M + 1 = d * (M / d) + (M % d + 1) := by
        conv_lhs => rw [← Nat.div_add_mod M d]
        ring
      have hmod : (M + 1) % d = (M % d + 1) % d := by
        rw [hsplit, Nat.mul_add_mod]
      have hdiv : (M + 1) / d = M / d + (M % d + 1) / d := by
        rw [hsplit, Nat.mul_add_div hd]
      rw [Finset.range_add_one, Finset.filter_insert]
      by_cases hc : M % d + 1 = d
      · rw [hc] at hmod hdiv
        rw [Nat.mod_self] at hmod
        rw [Nat.div_self hd] at hdiv
        by_cases h : M % d = r
        · rw [if_pos h, Finset.card_insert_of_notMem (by simp), ih,
              hdiv, hmod]
          split_ifs <;> omega
        · rw [if_neg h, ih, hdiv, hmod]; split_ifs <;> omega
      · have hlt' : M % d + 1 < d := by omega
        rw [Nat.mod_eq_of_lt hlt'] at hmod
        rw [Nat.div_eq_of_lt hlt'] at hdiv
        by_cases h : M % d = r
        · rw [if_pos h, Finset.card_insert_of_notMem (by simp), ih,
              hdiv, hmod]
          split_ifs <;> omega
        · rw [if_neg h, ih, hdiv, hmod]; split_ifs <;> omega

/-- The count is within `1` of `L/d` — the form the proof uses.  The
constant is `1`, not the `2` a tent-sum argument gives. -/
theorem card_residue_close (d r L : ℕ) (hd : 0 < d) (hr : r < d) :
    |(((Finset.range L).filter (fun i => i % d = r)).card : ℝ)
      - (L : ℝ) / d| ≤ 1 := by
  rw [card_residue_range d r L hd hr]
  have hd0 : (0 : ℝ) < d := by exact_mod_cast hd
  have e1 : (d : ℝ) * ((L / d : ℕ) : ℝ) + ((L % d : ℕ) : ℝ) = (L : ℝ) := by
    exact_mod_cast congrArg (fun n : ℕ => (n : ℝ)) (Nat.div_add_mod L d)
  have e2 : ((L % d : ℕ) : ℝ) < d := by
    exact_mod_cast Nat.mod_lt L hd
  have e3 : (0 : ℝ) ≤ ((L % d : ℕ) : ℝ) := Nat.cast_nonneg _
  have key : (L : ℝ) / d = ((L / d : ℕ) : ℝ) + ((L % d : ℕ) : ℝ) / d := by
    field_simp
    linarith [e1]
  have hf0 : (0 : ℝ) ≤ ((L % d : ℕ) : ℝ) / d := by positivity
  have hf1 : ((L % d : ℕ) : ℝ) / d < 1 := (div_lt_one hd0).mpr e2
  rw [key]
  split_ifs <;> push_cast <;> rw [abs_le] <;> constructor <;> linarith

/-! ### 3. The `m` cancels

The cell is `m` residue classes modulo `2Q`, each meeting the band in
`ell` terms, so `n_c = m * ell`.  A previous printing read the
resolution as `O(1/n_c)`.  It is `O(1/ell)`, and the two differ by `m`,
which is `5760` at depth `0` and `1` at depth `5`.

The reason the `m` does not help is here: the error is `C * ell` per
block, there are `m^2` blocks, and the normaliser is
`n_c^2 = m^2 * ell^2`.  The `m^2` divides out **identically**, not
asymptotically. -/

theorem blocks_cancel (m ell : ℕ) (hm : 0 < m) (hell : 0 < ell)
    (blockErr : ℝ) :
    |(m : ℝ) ^ 2 * blockErr| / ((m : ℝ) * ell) ^ 2
      = |blockErr| / (ell : ℝ) ^ 2 := by
  have hm0 : (0 : ℝ) < m := by exact_mod_cast hm
  have he0 : (0 : ℝ) < ell := by exact_mod_cast hell
  rw [abs_mul, abs_of_nonneg (by positivity : (0 : ℝ) ≤ (m : ℝ) ^ 2), mul_pow]
  field_simp

/-- The form the proof needs: a per-block error of `C * ell` over
`m^2` blocks gives a relative error `C / ell`, with no `m` anywhere. -/
theorem relative_error_has_no_m (m ell : ℕ) (C : ℝ) (hm : 0 < m)
    (hell : 0 < ell) (blockErr : ℝ) (hblock : |blockErr| ≤ C * ell) :
    |(m : ℝ) ^ 2 * blockErr| / ((m : ℝ) * ell) ^ 2 ≤ C / ell := by
  rw [blocks_cancel m ell hm hell blockErr]
  have he0 : (0 : ℝ) < ell := by exact_mod_cast hell
  rw [div_le_div_iff₀ (by positivity) he0]
  nlinarith [hblock, he0]

end Scaleinv

#print axioms Scaleinv.rankin_head
#print axioms Scaleinv.rankin_tail
#print axioms Scaleinv.card_residue_range
#print axioms Scaleinv.card_residue_close
#print axioms Scaleinv.blocks_cancel
#print axioms Scaleinv.relative_error_has_no_m

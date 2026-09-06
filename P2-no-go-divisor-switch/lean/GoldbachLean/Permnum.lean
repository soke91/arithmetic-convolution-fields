/-
  P4, `lem:permfloor` — the `eq:permnum` half, proved.

  WHAT IS PROVED

  For a field `Z` on a band of `n` places and a cell size `nc` with
  `0 < nc < n`, the average over **all** `nc`-subsets of the squared
  deviation of the subset mean from the band mean is

      ((n - nc) / (nc (n - 1))) * s^2,   s^2 = (1/n) * sum (Z - mbar)^2.

  That is `eq:permnum` as the paper prints it, including the two
  conventions the paper has to keep straight:

    * `s^2` carries divisor `n`, not `n - 1`;
    * the finite-population factor carries `n - 1`, not `n`.

  Mixing them would be the error.  The identity below holds only with
  this pairing, so the formalisation pins the convention instead of
  restating it.

  WHY THIS ONE AND NOT `eq:permfloor`

  `eq:permfloor` is this quantity's expectation under the sign null;
  `eq:permnum` is the quantity itself.  This repository identified the
  two once (`r17` S5).  Only the second is a finite combinatorial
  identity.  Proving the reachable one and saying which it is beats a
  formalisation that lets a reader believe the whole lemma is checked.

  THE ARGUMENT

  Centre: `Y = Z - mbar` has `sum Y = 0`, and for `#S = nc`,
  `mean S Z - mbar = (sum_{N in S} Y N) / nc`.  Double count:

      sum_S (sum_{N in S} Y)^2
        = sum_{N,N'} Y N * Y N' * #{S : #S = nc, N in S, N' in S}.

  The count is `C(n-1, nc-1)` on the diagonal and `C(n-2, nc-2)` off
  it.  With `sum_{N != N'} Y N Y N' = (sum Y)^2 - sum Y^2 = -sum Y^2`,
  the diagonal and off-diagonal counts enter with opposite signs and
  Pascal collapses them to `C(n-2, nc-1)`.  The rest is the binomial
  identity `n (n-1) C(n-2, nc-1) = nc (n-nc) C(n, nc)`.
-/

import Mathlib

open Finset

namespace Permnum

variable {n : ℕ}

/-- The mean of `Z` over a subset. -/
noncomputable def mean (S : Finset (Fin n)) (Z : Fin n → ℝ) : ℝ :=
  (∑ N ∈ S, Z N) / S.card

/-- The mean over the whole band. -/
noncomputable def mbar (Z : Fin n → ℝ) : ℝ := mean Finset.univ Z

/-! ### 1. The double count

The only step that touches the subset family.  Everything after it is
arithmetic on binomial coefficients. -/

theorem sum_powersetCard_double (nc : ℕ) (f : Fin n → Fin n → ℝ) :
    ∑ S ∈ Finset.univ.powersetCard nc, (∑ N ∈ S, ∑ N' ∈ S, f N N')
      = ∑ N, ∑ N', f N N' *
          (((Finset.univ.powersetCard nc).filter
            (fun S => N ∈ S ∧ N' ∈ S)).card : ℝ) := by
  have expand : ∀ S : Finset (Fin n),
      (∑ N ∈ S, ∑ N' ∈ S, f N N')
        = ∑ N, ∑ N', if N ∈ S ∧ N' ∈ S then f N N' else 0 := by
    intro S
    have inner : ∀ N : Fin n,
        (if N ∈ S then (∑ N' ∈ S, f N N') else 0)
          = ∑ N', if N ∈ S ∧ N' ∈ S then f N N' else 0 := by
      intro N
      by_cases hN : N ∈ S
      · simp only [hN, if_true, true_and]
        rw [Finset.sum_ite_mem, Finset.univ_inter]
      · simp [hN]
    rw [← Finset.sum_congr rfl (fun N (_ : N ∈ Finset.univ) => inner N),
        Finset.sum_ite_mem, Finset.univ_inter]
  calc ∑ S ∈ Finset.univ.powersetCard nc, (∑ N ∈ S, ∑ N' ∈ S, f N N')
      = ∑ S ∈ Finset.univ.powersetCard nc,
          ∑ N, ∑ N', if N ∈ S ∧ N' ∈ S then f N N' else 0 :=
        Finset.sum_congr rfl (fun S _ => expand S)
    _ = ∑ N, ∑ N', ∑ S ∈ Finset.univ.powersetCard nc,
          (if N ∈ S ∧ N' ∈ S then f N N' else 0) := by
        rw [Finset.sum_comm]
        exact Finset.sum_congr rfl (fun N _ => Finset.sum_comm)
    _ = _ := by
        refine Finset.sum_congr rfl (fun N _ => Finset.sum_congr rfl (fun N' _ => ?_))
        rw [Finset.sum_ite, Finset.sum_const, Finset.sum_const_zero, add_zero,
            nsmul_eq_mul, mul_comm]

/-! ### 2. The counts -/

/-- Subsets of size `nc` containing a given place: `C(n-1, nc-1)`. -/
theorem count_diag (nc : ℕ) (hnc : 0 < nc) (hle : nc ≤ n) (N : Fin n) :
    ((Finset.univ.powersetCard nc).filter (fun S => N ∈ S ∧ N ∈ S)).card
      = Nat.choose (n - 1) (nc - 1) := by
  have hfil : ((Finset.univ.powersetCard nc).filter (fun S => N ∈ S ∧ N ∈ S))
      = (Finset.univ.powersetCard nc).filter (({N} : Finset (Fin n)) ⊆ ·) := by
    apply Finset.filter_congr
    intro S _
    simp
  rw [hfil, Finset.card_filter_powersetCard_subset _ _ _ (Finset.subset_univ _)
        (by simp only [Finset.card_singleton]; omega)]
  simp [Finset.card_univ]

/-- Subsets of size `nc` containing two distinct places: `C(n-2, nc-2)`,
and `0` when `nc < 2` — the branch a formula in `n - 2`, `nc - 2` would
silently get wrong, since `nc - 2 = 0` truncates and `C(n-2, 0) = 1`. -/
theorem count_off (nc : ℕ) (hle : nc ≤ n) {N N' : Fin n} (hne : N ≠ N') :
    ((Finset.univ.powersetCard nc).filter (fun S => N ∈ S ∧ N' ∈ S)).card
      = if 2 ≤ nc then Nat.choose (n - 2) (nc - 2) else 0 := by
  have hpair : ({N, N'} : Finset (Fin n)).card = 2 := by
    rw [Finset.card_insert_of_notMem (by simpa using hne), Finset.card_singleton]
  have hfil : ((Finset.univ.powersetCard nc).filter (fun S => N ∈ S ∧ N' ∈ S))
      = (Finset.univ.powersetCard nc).filter (({N, N'} : Finset (Fin n)) ⊆ ·) := by
    apply Finset.filter_congr
    intro S _
    simp [Finset.insert_subset_iff]
  by_cases h2 : 2 ≤ nc
  · rw [hfil, Finset.card_filter_powersetCard_subset _ _ _ (Finset.subset_univ _)
          (by rw [hpair]; exact h2)]
    simp [hpair, Finset.card_univ, h2]
  · rw [hfil, if_neg h2]
    apply Finset.card_eq_zero.mpr
    apply Finset.filter_eq_empty_iff.mpr
    intro S hS
    rw [Finset.mem_powersetCard_univ] at hS
    intro hsub
    have := Finset.card_le_card hsub
    rw [hpair, hS] at this
    omega

/-! ### 3. Pascal collapses the two counts -/

/-- `C(n-1, nc-1) - (off-diagonal count) = C(n-2, nc-1)`, in both
branches.  At `nc = 1` this is `1 - 0 = 1`; at `nc >= 2` it is Pascal.
Stated as an addition so that no natural subtraction appears. -/
theorem counts_pascal (nc : ℕ) (hnc : 0 < nc) (hlt : nc < n) :
    Nat.choose (n - 1) (nc - 1)
      = (if 2 ≤ nc then Nat.choose (n - 2) (nc - 2) else 0)
        + Nat.choose (n - 2) (nc - 1) := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 2 := ⟨n - 2, by omega⟩
  by_cases h2 : 2 ≤ nc
  · obtain ⟨j, rfl⟩ : ∃ j, nc = j + 2 := ⟨nc - 2, by omega⟩
    simp only [if_pos h2]
    show Nat.choose (k + 1) (j + 1) = Nat.choose k j + Nat.choose k (j + 1)
    rw [Nat.choose_succ_succ]
  · obtain rfl : nc = 1 := by omega
    simp

/-! ### 4. The binomial arithmetic -/

/-- `n (n-1) C(n-2, nc-1) = nc (n-nc) C(n, nc)`. -/
theorem choose_arith (nc : ℕ) (hnc : 0 < nc) (hlt : nc < n) :
    n * (n - 1) * Nat.choose (n - 2) (nc - 1)
      = nc * (n - nc) * Nat.choose n nc := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 2 := ⟨n - 2, by omega⟩
  obtain ⟨j, rfl⟩ : ∃ j, nc = j + 1 := ⟨nc - 1, by omega⟩
  have hj : j ≤ k := by omega
  show (k + 2) * (k + 1) * Nat.choose k j
      = (j + 1) * (k + 2 - (j + 1)) * Nat.choose (k + 2) (j + 1)
  -- `Nat.succ_mul_choose_eq` is not in this Mathlib.  Both of the two
  -- lemmas below have the SAME right-hand side, `C(k+2,j) * (k+2-j)`,
  -- so composing them gives the identity that is wanted.
  have e1 : Nat.choose (k + 2) (j + 1) * (j + 1) = Nat.choose (k + 1) j * (k + 2) :=
    (Nat.choose_succ_right_eq (k + 2) j).trans (Nat.choose_mul_succ_eq (k + 1) j).symm
  have e2 : Nat.choose k j * (k + 1) = Nat.choose (k + 1) j * (k + 1 - j) :=
    Nat.choose_mul_succ_eq k j
  have hsub : k + 2 - (j + 1) = k + 1 - j := by omega
  rw [hsub]
  calc (k + 2) * (k + 1) * Nat.choose k j
      = (k + 2) * (Nat.choose k j * (k + 1)) := by ring
    _ = (k + 2) * (Nat.choose (k + 1) j * (k + 1 - j)) := by rw [e2]
    _ = (Nat.choose (k + 1) j * (k + 2)) * (k + 1 - j) := by ring
    _ = (Nat.choose (k + 2) (j + 1) * (j + 1)) * (k + 1 - j) := by rw [e1]
    _ = (j + 1) * (k + 1 - j) * Nat.choose (k + 2) (j + 1) := by ring

/-! ### 5. `eq:permnum` -/

-- PAPER: P4 lem:permfloor  (eq:permnum)
theorem permnum_closed_form (nc : ℕ) (hnc : 0 < nc) (hlt : nc < n)
    (Z : Fin n → ℝ) :
    (∑ S ∈ Finset.univ.powersetCard nc, (mean S Z - mbar Z) ^ 2)
        / (Nat.choose n nc)
      = ((n : ℝ) - nc) / (nc * ((n : ℝ) - 1))
        * ((∑ N, (Z N - mbar Z) ^ 2) / n) := by
  have hn : 0 < n := hnc.trans hlt
  have hnR : ((n : ℝ)) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  have hncR : ((nc : ℝ)) ≠ 0 := Nat.cast_ne_zero.mpr hnc.ne'
  set Y : Fin n → ℝ := fun N => Z N - mbar Z with hYdef
  -- The centred field sums to zero.  This is what kills the F-term below.
  have hY0 : ∑ N, Y N = 0 := by
    simp only [hYdef, Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ,
      Fintype.card_fin, nsmul_eq_mul, mbar, mean, Finset.card_univ, Fintype.card_fin]
    field_simp
    ring
  -- On a cell of the right size the deviation is the centred sum over n_c.
  have hdev : ∀ S ∈ Finset.univ.powersetCard nc,
      (mean S Z - mbar Z) = (∑ N ∈ S, Y N) / nc := by
    intro S hS
    rw [Finset.mem_powersetCard_univ] at hS
    have : (S.card : ℝ) = nc := by rw [hS]
    simp only [hYdef, Finset.sum_sub_distrib, Finset.sum_const, mean, hS,
      nsmul_eq_mul]
    field_simp
  set D : ℕ := Nat.choose (n - 1) (nc - 1) with hD
  set F : ℕ := (if 2 ≤ nc then Nat.choose (n - 2) (nc - 2) else 0) with hF
  set G : ℕ := Nat.choose (n - 2) (nc - 1) with hG
  have hpascal : (D : ℝ) = F + G := by
    exact_mod_cast congrArg (Nat.cast : ℕ → ℝ) (counts_pascal nc hnc hlt)
  -- The pair count, in one expression covering both cases.
  have hcount : ∀ N N' : Fin n,
      (((Finset.univ.powersetCard nc).filter
        (fun S => N ∈ S ∧ N' ∈ S)).card : ℝ)
        = if N = N' then (D : ℝ) else (F : ℝ) := by
    intro N N'
    by_cases h : N = N'
    · subst h; rw [if_pos rfl, count_diag nc hnc hlt.le N]
    · rw [if_neg h, count_off nc hlt.le h]
  -- The whole numerator collapses to G * sum Y^2.
  have hnum : (∑ S ∈ Finset.univ.powersetCard nc, (mean S Z - mbar Z) ^ 2)
      = (G : ℝ) * (∑ N, Y N ^ 2) / nc ^ 2 := by
    have step1 : (∑ S ∈ Finset.univ.powersetCard nc, (mean S Z - mbar Z) ^ 2)
        = (∑ S ∈ Finset.univ.powersetCard nc,
            (∑ N ∈ S, ∑ N' ∈ S, Y N * Y N')) / nc ^ 2 := by
      rw [Finset.sum_div]
      refine Finset.sum_congr rfl (fun S hS => ?_)
      rw [hdev S hS, div_pow, ← Finset.sum_mul_sum, sq]
    rw [step1, sum_powersetCard_double nc (fun N N' => Y N * Y N')]
    congr 1
    have hsplit : ∀ N N' : Fin n,
        Y N * Y N' * (if N = N' then (D : ℝ) else (F : ℝ))
          = Y N * Y N' * F + (if N = N' then Y N * Y N' * G else 0) := by
      intro N N'
      by_cases h : N = N'
      · subst h; rw [if_pos rfl, if_pos rfl, hpascal]; ring
      · rw [if_neg h, if_neg h]; ring
    calc ∑ N, ∑ N', Y N * Y N' *
            (((Finset.univ.powersetCard nc).filter
              (fun S => N ∈ S ∧ N' ∈ S)).card : ℝ)
        = ∑ N, ∑ N', (Y N * Y N' * F + (if N = N' then Y N * Y N' * G else 0)) := by
          refine Finset.sum_congr rfl (fun N _ => Finset.sum_congr rfl (fun N' _ => ?_))
          rw [hcount N N', hsplit N N']
      _ = (∑ N, ∑ N', Y N * Y N' * F)
            + ∑ N, ∑ N', (if N = N' then Y N * Y N' * G else 0) := by
          rw [← Finset.sum_add_distrib]
          exact Finset.sum_congr rfl (fun N _ => Finset.sum_add_distrib)
      _ = (G : ℝ) * ∑ N, Y N ^ 2 := by
          have hd : (∑ N, ∑ N', Y N * Y N') = (∑ N, Y N) * (∑ N', Y N') :=
            (Finset.sum_mul_sum _ _ _ _).symm
          have hzero : (∑ N, ∑ N', Y N * Y N' * (F : ℝ)) = 0 := by
            calc (∑ N, ∑ N', Y N * Y N' * (F : ℝ))
                = (∑ N, ∑ N', Y N * Y N') * F := by
                  rw [Finset.sum_mul]
                  exact Finset.sum_congr rfl
                    (fun N _ => (Finset.sum_mul _ _ _).symm)
              _ = 0 := by rw [hd, hY0]; ring
          rw [hzero, zero_add, Finset.mul_sum]
          refine Finset.sum_congr rfl (fun N _ => ?_)
          simp only [Finset.sum_ite_eq, Finset.mem_univ, if_true]
          ring
  rw [hnum]
  -- What is left is the binomial identity, cast to the reals.
  have harith : (n : ℝ) * ((n : ℝ) - 1) * G = nc * ((n : ℝ) - nc) * Nat.choose n nc := by
    have h := congrArg (Nat.cast : ℕ → ℝ) (choose_arith nc hnc hlt)
    push_cast [Nat.cast_sub hlt.le, Nat.cast_sub hn] at h
    exact h
  have hc : (0 : ℝ) < Nat.choose n nc := by
    exact_mod_cast Nat.choose_pos hlt.le
  have h3 : ((n : ℝ) - 1) ≠ 0 := by
    have hn1 : 1 < n := by omega
    have : (1 : ℝ) < n := by exact_mod_cast hn1
    linarith
  -- Solve the binomial identity for G and substitute; what is left is `ring`.
  have hGval : (G : ℝ)
      = (nc : ℝ) * ((n : ℝ) - nc) * (Nat.choose n nc) / ((n : ℝ) * ((n : ℝ) - 1)) := by
    field_simp
    linear_combination harith
  have hYsq : (∑ N, Y N ^ 2) = ∑ N, (Z N - mbar Z) ^ 2 := rfl
  rw [hYsq, hGval]
  field_simp

end Permnum

#print axioms Permnum.permnum_closed_form
#print axioms Permnum.sum_powersetCard_double
#print axioms Permnum.choose_arith

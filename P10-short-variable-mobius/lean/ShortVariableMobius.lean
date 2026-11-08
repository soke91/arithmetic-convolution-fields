/-
  ShortVariableMobius.lean --- the exact finite steps of "A ternary Moebius sum with one
  short variable".

  WHAT THIS COVERS, in the paper's numbering.  The Moebius factorisation (11); the
  decomposition over the common divisor e of n_1 and n_2, Lemma 7 with its definition (10) of
  W(e); the counting bound (12) of Lemma 8, the two bounds (13) on X_N, and the three-range
  split (14) of Lemma 8 with the two finite-sum inequalities it rests on; the canonical divisor
  splitting (15) of Lemma 10; and the two-residue-class count behind Lemma 11.

  WHAT IS HERE AND WHAT IS NOT.  The paper's analytic input is NOT formalized and is not
  attempted here: Davenport's bound for the exponential sum of the Moebius function, the
  Cauchy-Schwarz and Parseval step for the bilinear factor, the divisor moments, the explicit
  bound on the largest divisor function, and the asymptotic estimates that assemble the main
  theorem.  What is formalized is the part of the proof that is a finite identity or a finite
  inequality.  The main theorem is NOT formalized.

  NO `sorry` APPEARS IN THIS FILE.  `#print axioms` at the end prints the axiom set of every
  theorem; `axioms.txt` is that output.

  Elaborate with the toolchain `leanprover/lean4:v4.33.1` and a mathlib of the same revision.
  From this directory, with `lean` and `lake` on PATH:

      lake env lean ShortVariableMobius.lean

  `lakefile.toml`, `lean-toolchain` and `lake-manifest.json` sit beside this file; `lake` will
  fetch and build mathlib on the first run if no build is present.

  A NOTE ON `Nat.ArithmeticFunction`.  In this mathlib the namespace is `ArithmeticFunction`,
  not `Nat.ArithmeticFunction`, and the notations `mu` and `zeta` live in the scoped
  namespaces `ArithmeticFunction.Moebius` and `ArithmeticFunction.zeta`.
-/

import Mathlib

/-!
# The exact finite parts of a ternary Moebius sum with one short variable

Five groups of statements.

1. **The factorisation** (11): `moebius_mul_mul`.
2. **Counting**: `card_le_of_const_mod` -- which bounds the number of admissible values of the
   short variable in Lemma 8 -- and `card_le_of_two_classes`, the two-residue-class count
   behind Lemma 11.
3. **The decomposition** (Lemma 7, with (10)): `mem_triples`, `S_eq_sum_moebius_mul_W`.
4. **The tail** (Lemma 8, with (12), (13), (14)): `card_filter_triples_le`, `norm_W_le`,
   `le_XN`, `XN_le_mul_taustar`, `W_eq_zero_of_forall_not_dvd`, `sum_inv_sq_le`,
   `sum_inv_le_log_div`, `sum_norm_W_le`.
5. **The canonical divisor splitting** (15) of Lemma 10: `canonical_split`.
-/

namespace ShortVariableMobius

open Finset ArithmeticFunction
open scoped ArithmeticFunction.Moebius
open scoped ArithmeticFunction.zeta

/-! ## 1. The Moebius factorisation

The identity (11) of the paper, for all natural numbers: at a zero argument both sides vanish,
and at positive arguments it is the statement the paper proves.
-/

private lemma moebius_eq_zero_of_sq_dvd {m p : ℕ} (hp : p.Prime) (h : p * p ∣ m) :
    μ m = 0 :=
  moebius_eq_zero_of_not_squarefree fun hsq => hp.one_lt.ne' (Nat.isUnit_iff.mp (hsq p h))

/-- The factorisation (11): `μ(n₁n₂n₃) = μ(n₁)μ(n₂)μ(n₃)[(n₁,n₂)=1][(n₁n₂,n₃)=1]`. -/
theorem moebius_mul_mul (n₁ n₂ n₃ : ℕ) :
    μ (n₁ * n₂ * n₃) =
      μ n₁ * μ n₂ * μ n₃ * (if Nat.Coprime n₁ n₂ then 1 else 0) *
        (if Nat.Coprime (n₁ * n₂) n₃ then 1 else 0) := by
  by_cases h₁ : Nat.Coprime n₁ n₂
  · by_cases h₂ : Nat.Coprime (n₁ * n₂) n₃
    · rw [if_pos h₁, if_pos h₂, mul_one, mul_one,
        isMultiplicative_moebius.map_mul_of_coprime h₂,
        isMultiplicative_moebius.map_mul_of_coprime h₁]
    · obtain ⟨p, hp, hpd⟩ := Nat.exists_prime_and_dvd h₂
      rw [if_neg h₂, mul_zero]
      exact moebius_eq_zero_of_sq_dvd hp
        (mul_dvd_mul (hpd.trans (Nat.gcd_dvd_left _ _)) (hpd.trans (Nat.gcd_dvd_right _ _)))
  · obtain ⟨p, hp, hpd⟩ := Nat.exists_prime_and_dvd h₁
    rw [if_neg h₁, mul_zero, zero_mul]
    exact moebius_eq_zero_of_sq_dvd hp
      (dvd_mul_of_dvd_left
        (mul_dvd_mul (hpd.trans (Nat.gcd_dvd_left _ _)) (hpd.trans (Nat.gcd_dvd_right _ _))) n₃)

/-! ## 2. Two counting bounds

The paper counts an integer variable confined to one residue class (Lemma 8) or to two
(Lemma 11).  Both counts come from the same observation: a set of integers of `[0, M]` that is
constant modulo `L` has at most `⌊M/L⌋ + 1` elements.
-/

/-- A set of integers of `[0, M]` all congruent to one another modulo `L` has at most
`⌊M/L⌋ + 1` elements. -/
theorem card_le_of_const_mod {M L : ℕ} (s : Finset ℕ)
    (hle : ∀ n ∈ s, n ≤ M) (hmod : ∀ m ∈ s, ∀ n ∈ s, m % L = n % L) :
    s.card ≤ M / L + 1 := by
  have h : s.card ≤ (Finset.range (M / L + 1)).card := by
    refine Finset.card_le_card_of_injOn (fun n => n / L) ?_ ?_
    · intro n hn
      simp only [Finset.mem_coe, Finset.mem_range] at hn ⊢
      exact Nat.lt_succ_of_le (Nat.div_le_div_right (hle n hn))
    · intro a ha b hb hab
      have ha' : a ∈ s := ha
      have hb' : b ∈ s := hb
      have hab' : a / L = b / L := hab
      calc a = L * (a / L) + a % L := (Nat.div_add_mod a L).symm
        _ = L * (b / L) + b % L := by rw [hab', hmod a ha' b hb']
        _ = b := Nat.div_add_mod b L
  simpa using h

private lemma modEq_lcm {q r m n : ℕ} (hq : m ≡ n [MOD q]) (hr : m ≡ n [MOD r]) :
    m ≡ n [MOD Nat.lcm q r] := by
  rcases le_total m n with h | h
  · exact (Nat.modEq_iff_dvd' h).mpr
      (Nat.lcm_dvd ((Nat.modEq_iff_dvd' h).mp hq) ((Nat.modEq_iff_dvd' h).mp hr))
  · exact ((Nat.modEq_iff_dvd' h).mpr
      (Nat.lcm_dvd ((Nat.modEq_iff_dvd' h).mp hq.symm)
        ((Nat.modEq_iff_dvd' h).mp hr.symm))).symm

/-- The two-residue-class count behind Lemma 11: a set of integers of `[0, M]` constant
modulo `q` and modulo `r` has at most `M (q,r)/(qr) + 1` elements. -/
theorem card_le_of_two_classes {M q r : ℕ} (hq : 0 < q) (hr : 0 < r) (s : Finset ℕ)
    (hle : ∀ n ∈ s, n ≤ M)
    (hmq : ∀ m ∈ s, ∀ n ∈ s, m % q = n % q)
    (hmr : ∀ m ∈ s, ∀ n ∈ s, m % r = n % r) :
    (s.card : ℝ) ≤ (M : ℝ) * (Nat.gcd q r : ℝ) / ((q : ℝ) * (r : ℝ)) + 1 := by
  have hgl := Nat.gcd_mul_lcm q r
  have hL : 0 < Nat.lcm q r := by
    rcases Nat.eq_zero_or_pos (Nat.lcm q r) with h0 | h0
    · rw [h0, Nat.mul_zero] at hgl
      exact absurd hgl.symm (Nat.mul_pos hq hr).ne'
    · exact h0
  have hg : 0 < Nat.gcd q r := by
    rcases Nat.eq_zero_or_pos (Nat.gcd q r) with h0 | h0
    · rw [h0, Nat.zero_mul] at hgl
      exact absurd hgl.symm (Nat.mul_pos hq hr).ne'
    · exact h0
  have hcard : s.card ≤ M / Nat.lcm q r + 1 := by
    refine card_le_of_const_mod s hle ?_
    intro m hm n hn
    exact modEq_lcm (hmq m hm n hn) (hmr m hm n hn)
  have hgne : (Nat.gcd q r : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hg.ne'
  have hLne : (Nat.lcm q r : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hL.ne'
  have hqr : ((q : ℝ) * (r : ℝ)) = (Nat.gcd q r : ℝ) * (Nat.lcm q r : ℝ) := by
    exact_mod_cast hgl.symm
  have hqrne : ((q : ℝ) * (r : ℝ)) ≠ 0 := by
    rw [hqr]
    exact mul_ne_zero hgne hLne
  have h2 : (M : ℝ) / (Nat.lcm q r : ℝ)
      = (M : ℝ) * (Nat.gcd q r : ℝ) / ((q : ℝ) * (r : ℝ)) := by
    rw [div_eq_div_iff hLne hqrne, hqr]
    ring
  have hcast : (s.card : ℝ) ≤ ((M / Nat.lcm q r : ℕ) : ℝ) + 1 := by exact_mod_cast hcard
  have hdiv : ((M / Nat.lcm q r : ℕ) : ℝ) ≤ (M : ℝ) / (Nat.lcm q r : ℝ) := Nat.cast_div_le
  rw [← h2]
  linarith

/-! ## 3. The decomposition over the common divisor

`triples N H` is the set of ordered triples of positive integers with `n₁ + n₂ + n₃ = N` and
`n₃ ≤ H`, which is the range of summation of the paper's `S_{u,v}(H,N)`.  `wterm` is the
summand of (10), `W` is (10), and `S_eq_sum_moebius_mul_W` is Lemma 7.
-/

/-- The ordered triples of positive integers with `n₁ + n₂ + n₃ = N` and `n₃ ≤ H`. -/
def triples (N H : ℕ) : Finset (ℕ × ℕ × ℕ) :=
  (Finset.Icc 1 N ×ˢ Finset.Icc 1 N ×ˢ Finset.Icc 1 N).filter
    (fun t => t.1 + t.2.1 + t.2.2 = N ∧ t.2.2 ≤ H)

theorem mem_triples {N H : ℕ} {t : ℕ × ℕ × ℕ} :
    t ∈ triples N H ↔
      0 < t.1 ∧ 0 < t.2.1 ∧ 0 < t.2.2 ∧ t.1 + t.2.1 + t.2.2 = N ∧ t.2.2 ≤ H := by
  unfold triples
  simp only [Finset.mem_filter, Finset.mem_product, Finset.mem_Icc]
  omega

/-- The paper's `S_{u,v}(H,N)`, the sum over ordered triples of positive integers. -/
noncomputable def S (u v : ℕ → ℂ) (N H : ℕ) : ℂ :=
  ∑ t ∈ triples N H, u t.1 * v t.2.1 * ((μ (t.1 * t.2.1 * t.2.2) : ℤ) : ℂ)

/-- The summand of the paper's `W(e)`, the definition (10). -/
noncomputable def wterm (u v : ℕ → ℂ) (t : ℕ × ℕ × ℕ) : ℂ :=
  u t.1 * v t.2.1 * ((μ t.1 * μ t.2.1 * μ t.2.2 : ℤ) : ℂ) *
    (if Nat.Coprime (t.1 * t.2.1) t.2.2 then 1 else 0)

/-- The paper's `W(e)`, the definition (10). -/
noncomputable def W (u v : ℕ → ℂ) (N H e : ℕ) : ℂ :=
  ∑ t ∈ (triples N H).filter (fun t => e ∣ t.1 ∧ e ∣ t.2.1), wterm u v t

private lemma S_def (u v : ℕ → ℂ) (N H : ℕ) :
    S u v N H = ∑ t ∈ triples N H, u t.1 * v t.2.1 * ((μ (t.1 * t.2.1 * t.2.2) : ℤ) : ℂ) :=
  rfl

private lemma W_def (u v : ℕ → ℂ) (N H e : ℕ) :
    W u v N H e = ∑ t ∈ (triples N H).filter (fun t => e ∣ t.1 ∧ e ∣ t.2.1), wterm u v t :=
  rfl

private lemma wterm_def (u v : ℕ → ℂ) (t : ℕ × ℕ × ℕ) :
    wterm u v t = u t.1 * v t.2.1 * ((μ t.1 * μ t.2.1 * μ t.2.2 : ℤ) : ℂ) *
      (if Nat.Coprime (t.1 * t.2.1) t.2.2 then 1 else 0) :=
  rfl

private lemma sum_divisors_moebius (n : ℕ) :
    ∑ d ∈ n.divisors, μ d = if n = 1 then 1 else 0 := by
  have h : (μ * ζ : ArithmeticFunction ℤ) n = (1 : ArithmeticFunction ℤ) n := by
    rw [moebius_mul_coe_zeta]
  rwa [coe_mul_zeta_apply, one_apply] at h

private lemma filter_dvd_eq_divisors_gcd {N n₁ n₂ : ℕ} (h1 : 0 < n₁) (h2 : n₁ ≤ N) :
    (Finset.Icc 1 N).filter (fun e => e ∣ n₁ ∧ e ∣ n₂) = (Nat.gcd n₁ n₂).divisors := by
  ext e
  constructor
  · intro he
    rw [Finset.mem_filter, Finset.mem_Icc] at he
    obtain ⟨-, hd1, hd2⟩ := he
    rw [Nat.mem_divisors]
    refine ⟨Nat.dvd_gcd hd1 hd2, ?_⟩
    intro hgz
    rw [Nat.gcd_eq_zero_iff] at hgz
    omega
  · intro he
    rw [Nat.mem_divisors] at he
    obtain ⟨hgd, -⟩ := he
    have hd1 : e ∣ n₁ := hgd.trans (Nat.gcd_dvd_left _ _)
    have hd2 : e ∣ n₂ := hgd.trans (Nat.gcd_dvd_right _ _)
    rw [Finset.mem_filter, Finset.mem_Icc]
    refine ⟨⟨?_, ?_⟩, hd1, hd2⟩
    · rcases Nat.eq_zero_or_pos e with rfl | hepos
      · exact absurd (Nat.eq_zero_of_zero_dvd hd1) h1.ne'
      · exact hepos
    · exact le_trans (Nat.le_of_dvd h1 hd1) h2

/-- Lemma 7: `S_{u,v}(H,N) = ∑_{e ≤ N} μ(e) W(e)`. -/
theorem S_eq_sum_moebius_mul_W (u v : ℕ → ℂ) (N H : ℕ) :
    S u v N H = ∑ e ∈ Finset.Icc 1 N, ((μ e : ℤ) : ℂ) * W u v N H e := by
  have step1 : ∀ t ∈ triples N H,
      u t.1 * v t.2.1 * ((μ (t.1 * t.2.1 * t.2.2) : ℤ) : ℂ) =
        ∑ e ∈ (Finset.Icc 1 N).filter (fun e => e ∣ t.1 ∧ e ∣ t.2.1),
          ((μ e : ℤ) : ℂ) * wterm u v t := by
    intro t ht
    rw [mem_triples] at ht
    obtain ⟨hp1, hp2, hp3, hsum, hH⟩ := ht
    have hle : t.1 ≤ N := by omega
    have hcop : (if Nat.Coprime t.1 t.2.1 then (1 : ℤ) else 0) =
        if Nat.gcd t.1 t.2.1 = 1 then (1 : ℤ) else 0 := by
      by_cases hc : Nat.Coprime t.1 t.2.1
      · rw [if_pos hc]
      · rw [if_neg hc]
    have hsumZ : ∑ e ∈ (Finset.Icc 1 N).filter (fun e => e ∣ t.1 ∧ e ∣ t.2.1), (μ e : ℤ) =
        if Nat.Coprime t.1 t.2.1 then 1 else 0 := by
      rw [filter_dvd_eq_divisors_gcd hp1 hle, sum_divisors_moebius, hcop]
    have hsumC : (∑ e ∈ (Finset.Icc 1 N).filter (fun e => e ∣ t.1 ∧ e ∣ t.2.1),
        ((μ e : ℤ) : ℂ)) = ((if Nat.Coprime t.1 t.2.1 then (1 : ℤ) else 0 : ℤ) : ℂ) := by
      rw [← hsumZ]
      push_cast
      try rfl
    rw [← Finset.sum_mul, hsumC, wterm_def, moebius_mul_mul]
    split_ifs <;> push_cast <;> ring
  have hswap := Finset.sum_comm' (s := triples N H)
    (t := fun t : ℕ × ℕ × ℕ => (Finset.Icc 1 N).filter (fun e => e ∣ t.1 ∧ e ∣ t.2.1))
    (t' := Finset.Icc 1 N)
    (s' := fun e : ℕ => (triples N H).filter (fun t => e ∣ t.1 ∧ e ∣ t.2.1))
    (f := fun (t : ℕ × ℕ × ℕ) (e : ℕ) => ((μ e : ℤ) : ℂ) * wterm u v t)
    (by intro t e; simp only [Finset.mem_filter]; tauto)
  rw [S_def, Finset.sum_congr rfl step1, hswap]
  refine Finset.sum_congr rfl fun e _ => ?_
  rw [W_def, Finset.mul_sum]

/-! ## 4. The tail

`norm_W_le` is (12); `le_XN` and `XN_le_mul_taustar` are the two bounds (13);
`sum_norm_W_le` is the three-range split (14).  The weights are complex numbers of modulus at
most `1`, as in the paper.
-/

private lemma abs_moebius_le_one (n : ℕ) : |μ n| ≤ 1 := by
  by_cases h : Squarefree n
  · rw [moebius_apply_of_squarefree h, abs_pow, abs_neg, abs_one, one_pow]
  · rw [moebius_eq_zero_of_not_squarefree h, abs_zero]
    exact zero_le_one

private lemma norm_wterm_le_one {u v : ℕ → ℂ} (hu : ∀ n, ‖u n‖ ≤ 1) (hv : ∀ n, ‖v n‖ ≤ 1)
    (t : ℕ × ℕ × ℕ) : ‖wterm u v t‖ ≤ 1 := by
  have hmu : ‖((μ t.1 * μ t.2.1 * μ t.2.2 : ℤ) : ℂ)‖ ≤ 1 := by
    rw [Complex.norm_intCast]
    have h : |μ t.1 * μ t.2.1 * μ t.2.2| ≤ 1 := by
      have e1 : |μ t.1 * μ t.2.1 * μ t.2.2| = |μ t.1| * |μ t.2.1| * |μ t.2.2| := by
        rw [abs_mul, abs_mul]
      have b1 := abs_moebius_le_one t.1
      have b2 := abs_moebius_le_one t.2.1
      have b3 := abs_moebius_le_one t.2.2
      have n2 : (0 : ℤ) ≤ |μ t.2.1| := abs_nonneg _
      have n3 : (0 : ℤ) ≤ |μ t.2.2| := abs_nonneg _
      rw [e1]
      calc |μ t.1| * |μ t.2.1| * |μ t.2.2| ≤ 1 * 1 * 1 :=
            mul_le_mul (mul_le_mul b1 b2 n2 zero_le_one) b3 n3 (by norm_num)
        _ = 1 := by norm_num
    rw [← Int.cast_abs]
    exact_mod_cast h
  have hif : ‖(if Nat.Coprime (t.1 * t.2.1) t.2.2 then (1 : ℂ) else 0)‖ ≤ 1 := by
    split_ifs <;> simp
  rw [wterm_def]
  have p1 : ‖u t.1‖ * ‖v t.2.1‖ ≤ 1 * 1 :=
    mul_le_mul (hu _) (hv _) (norm_nonneg _) zero_le_one
  have p2 : ‖u t.1‖ * ‖v t.2.1‖ * ‖((μ t.1 * μ t.2.1 * μ t.2.2 : ℤ) : ℂ)‖ ≤ 1 * 1 * 1 :=
    mul_le_mul p1 hmu (norm_nonneg _) (by norm_num)
  calc ‖u t.1 * v t.2.1 * ((μ t.1 * μ t.2.1 * μ t.2.2 : ℤ) : ℂ) *
        (if Nat.Coprime (t.1 * t.2.1) t.2.2 then (1 : ℂ) else 0)‖
      = ‖u t.1‖ * ‖v t.2.1‖ * ‖((μ t.1 * μ t.2.1 * μ t.2.2 : ℤ) : ℂ)‖ *
          ‖(if Nat.Coprime (t.1 * t.2.1) t.2.2 then (1 : ℂ) else 0)‖ := by
        rw [norm_mul, norm_mul, norm_mul]
    _ ≤ 1 * 1 * 1 * 1 := mul_le_mul p2 hif (norm_nonneg _) (by norm_num)
    _ = 1 := by norm_num

/-- The number of triples of (10) with `e ∣ n₁` and `e ∣ n₂`, which is what (12) bounds: the
equation forces `n₃ ≡ N (mod e)`, so `n₃` takes at most `⌊H/e⌋ + 1` values and `n₁` at most
`⌊N/e⌋`. -/
theorem card_filter_triples_le (N H : ℕ) {e : ℕ} (he : 0 < e) :
    ((triples N H).filter (fun t => e ∣ t.1 ∧ e ∣ t.2.1)).card ≤ (H / e + 1) * (N / e) := by
  have hcardA : ((Finset.Icc 1 H).filter (fun n₃ => n₃ % e = N % e)).card ≤ H / e + 1 := by
    refine card_le_of_const_mod _ ?_ ?_
    · intro n hn
      rw [Finset.mem_filter, Finset.mem_Icc] at hn
      exact hn.1.2
    · intro m hm n hn
      rw [Finset.mem_filter] at hm hn
      rw [hm.2, hn.2]
  have hcardB : ((Finset.Icc 1 N).filter (fun n₁ => e ∣ n₁)).card ≤ N / e := by
    have h : ((Finset.Icc 1 N).filter (fun n₁ => e ∣ n₁)).card ≤
        (Finset.Icc 1 (N / e)).card := by
      refine Finset.card_le_card_of_injOn (fun n => n / e) ?_ ?_
      · intro n hn
        simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_Icc] at hn ⊢
        obtain ⟨⟨hn1, hnN⟩, hdvd⟩ := hn
        exact ⟨(Nat.one_le_div_iff he).mpr (Nat.le_of_dvd hn1 hdvd),
          Nat.div_le_div_right hnN⟩
      · intro a ha b hb hab
        have ha' : a ∈ (Finset.Icc 1 N).filter (fun n₁ => e ∣ n₁) := ha
        have hb' : b ∈ (Finset.Icc 1 N).filter (fun n₁ => e ∣ n₁) := hb
        rw [Finset.mem_filter] at ha' hb'
        have hab' : a / e = b / e := hab
        calc a = e * (a / e) := (Nat.mul_div_cancel' ha'.2).symm
          _ = e * (b / e) := by rw [hab']
          _ = b := Nat.mul_div_cancel' hb'.2
    simpa using h
  have hinj : ((triples N H).filter (fun t => e ∣ t.1 ∧ e ∣ t.2.1)).card ≤
      (((Finset.Icc 1 H).filter (fun n₃ => n₃ % e = N % e)) ×ˢ
        ((Finset.Icc 1 N).filter (fun n₁ => e ∣ n₁))).card := by
    refine Finset.card_le_card_of_injOn (fun t => (t.2.2, t.1)) ?_ ?_
    · intro t ht
      simp only [Finset.mem_coe, Finset.mem_filter, mem_triples, Finset.mem_product,
        Finset.mem_Icc] at ht ⊢
      obtain ⟨⟨h1, h2, h3, hsum, hH⟩, hd1, hd2⟩ := ht
      refine ⟨⟨⟨h3, hH⟩, ?_⟩, ⟨h1, by omega⟩, hd1⟩
      obtain ⟨a, ha⟩ := hd1
      obtain ⟨b, hb⟩ := hd2
      have hN : N = e * (a + b) + t.2.2 := by rw [← hsum, ha, hb]; ring
      rw [hN, Nat.mul_add_mod]
    · intro p hp q hq hpq
      have hp' : p ∈ (triples N H).filter (fun t => e ∣ t.1 ∧ e ∣ t.2.1) := hp
      have hq' : q ∈ (triples N H).filter (fun t => e ∣ t.1 ∧ e ∣ t.2.1) := hq
      rw [Finset.mem_filter, mem_triples] at hp' hq'
      obtain ⟨⟨-, -, -, hps, -⟩, -⟩ := hp'
      obtain ⟨⟨-, -, -, hqs, -⟩, -⟩ := hq'
      have hpq' : (p.2.2, p.1) = (q.2.2, q.1) := hpq
      rw [Prod.ext_iff] at hpq'
      obtain ⟨e1, e2⟩ := hpq'
      refine Prod.ext_iff.mpr ⟨e2, Prod.ext_iff.mpr ⟨?_, e1⟩⟩
      omega
  calc ((triples N H).filter (fun t => e ∣ t.1 ∧ e ∣ t.2.1)).card
      ≤ (((Finset.Icc 1 H).filter (fun n₃ => n₃ % e = N % e)) ×ˢ
          ((Finset.Icc 1 N).filter (fun n₁ => e ∣ n₁))).card := hinj
    _ = ((Finset.Icc 1 H).filter (fun n₃ => n₃ % e = N % e)).card *
          ((Finset.Icc 1 N).filter (fun n₁ => e ∣ n₁)).card := Finset.card_product _ _
    _ ≤ (H / e + 1) * (N / e) := Nat.mul_le_mul hcardA hcardB

/-- The counting bound (12): `|W(e)| ≤ (⌊H/e⌋ + 1)⌊N/e⌋`. -/
theorem norm_W_le {u v : ℕ → ℂ} (hu : ∀ n, ‖u n‖ ≤ 1) (hv : ∀ n, ‖v n‖ ≤ 1)
    (N H : ℕ) {e : ℕ} (he : 0 < e) :
    ‖W u v N H e‖ ≤ (((H / e + 1) * (N / e) : ℕ) : ℝ) := by
  rw [W_def]
  refine le_trans (norm_sum_le _ _) ?_
  refine le_trans (Finset.sum_le_sum fun t _ => norm_wterm_le_one hu hv t) ?_
  rw [Finset.sum_const, nsmul_eq_mul, mul_one]
  exact Nat.cast_le.mpr (card_filter_triples_le N H he)

/-- `X_N` of (13), with `τ(m) = 0` for `m ≤ 0` realised by truncated subtraction. -/
def XN (N H : ℕ) : ℕ := ∑ n₃ ∈ Finset.Icc 1 H, (N - n₃).divisors.card

/-- `τ*(N) = max_{1 ≤ m < N} τ(m)` of Lemma 6. -/
def taustar (N : ℕ) : ℕ := (Finset.Ico 1 N).sup (fun m => m.divisors.card)

private lemma XN_def (N H : ℕ) :
    XN N H = ∑ n₃ ∈ Finset.Icc 1 H, (N - n₃).divisors.card := rfl

/-- The upper bound of (13): `X_N ≤ H τ*(N)`. -/
theorem XN_le_mul_taustar (N H : ℕ) : XN N H ≤ H * taustar N := by
  have h : ∀ n₃ ∈ Finset.Icc 1 H, (N - n₃).divisors.card ≤ taustar N := by
    intro n₃ hn₃
    rw [Finset.mem_Icc] at hn₃
    rcases lt_or_ge n₃ N with hlt | hge
    · exact Finset.le_sup (f := fun m : ℕ => m.divisors.card)
        (Finset.mem_Ico.mpr ⟨by omega, by omega⟩)
    · rw [show N - n₃ = 0 from by omega]
      simp
  rw [XN_def]
  calc ∑ n₃ ∈ Finset.Icc 1 H, (N - n₃).divisors.card
      ≤ ∑ _n₃ ∈ Finset.Icc 1 H, taustar N := Finset.sum_le_sum h
    _ = H * taustar N := by
        rw [Finset.sum_const, Nat.card_Icc, smul_eq_mul, Nat.add_sub_cancel]

private lemma one_le_card_divisors {m : ℕ} (hm : 0 < m) : 1 ≤ m.divisors.card :=
  Finset.card_pos.mpr ⟨1, Nat.one_mem_divisors.mpr hm.ne'⟩

private lemma two_le_card_divisors {m : ℕ} (hm : 2 ≤ m) : 2 ≤ m.divisors.card := by
  have hsub : ({1, m} : Finset ℕ) ⊆ m.divisors := by
    intro x hx
    rw [Finset.mem_insert, Finset.mem_singleton] at hx
    rcases hx with hx1 | hxm
    · rw [hx1]
      exact Nat.one_mem_divisors.mpr (by omega)
    · rw [hxm]
      exact Nat.mem_divisors_self m (by omega)
  have hcard : ({1, m} : Finset ℕ).card = 2 := by
    rw [Finset.card_insert_of_notMem (by simp; omega), Finset.card_singleton]
  calc 2 = ({1, m} : Finset ℕ).card := hcard.symm
    _ ≤ m.divisors.card := Finset.card_le_card hsub

/-- The lower bound of (13): `H ≤ X_N`. -/
theorem le_XN {N H : ℕ} (hN : 3 ≤ N) (hH : 2 ≤ H) (hHN : H ≤ N) : H ≤ XN N H := by
  have hK2 : 2 ≤ min H (N - 1) := by omega
  have hsub : Finset.Icc 1 (min H (N - 1)) ⊆ Finset.Icc 1 H := by
    intro x hx
    rw [Finset.mem_Icc] at hx ⊢
    omega
  have hsplit : Finset.Icc 1 (min H (N - 1)) = insert 1 (Finset.Icc 2 (min H (N - 1))) := by
    ext x
    simp only [Finset.mem_Icc, Finset.mem_insert]
    omega
  have hone : 2 ≤ (N - 1).divisors.card := two_le_card_divisors (by omega)
  have hrest : ∀ n₃ ∈ Finset.Icc 2 (min H (N - 1)), 1 ≤ (N - n₃).divisors.card := by
    intro n₃ hn₃
    rw [Finset.mem_Icc] at hn₃
    exact one_le_card_divisors (by omega)
  have hlow : min H (N - 1) + 1 ≤
      ∑ n₃ ∈ Finset.Icc 1 (min H (N - 1)), (N - n₃).divisors.card := by
    rw [hsplit, Finset.sum_insert (by simp)]
    have h1 : (Finset.Icc 2 (min H (N - 1))).card ≤
        ∑ n₃ ∈ Finset.Icc 2 (min H (N - 1)), (N - n₃).divisors.card := by
      simpa using Finset.card_nsmul_le_sum (Finset.Icc 2 (min H (N - 1)))
        (fun n₃ => (N - n₃).divisors.card) 1 hrest
    rw [Nat.card_Icc] at h1
    omega
  have hmono : ∑ n₃ ∈ Finset.Icc 1 (min H (N - 1)), (N - n₃).divisors.card ≤ XN N H := by
    rw [XN_def]
    exact Finset.sum_le_sum_of_subset hsub
  omega

/-- If `W(e) ≠ 0` then `e` divides one of the integers `N - n₃` with `1 ≤ n₃ ≤ H` and
`n₃ < N`: the equation forces `e ∣ n₁ + n₂ = N - n₃`. -/
theorem W_eq_zero_of_forall_not_dvd {u v : ℕ → ℂ} {N H e : ℕ}
    (h : ∀ n₃, 1 ≤ n₃ → n₃ ≤ H → n₃ < N → ¬ e ∣ N - n₃) : W u v N H e = 0 := by
  rw [W_def]
  refine Finset.sum_eq_zero ?_
  intro t ht
  rw [Finset.mem_filter, mem_triples] at ht
  obtain ⟨⟨h1, h2, h3, hsum, hH⟩, hd1, hd2⟩ := ht
  exact absurd (by
    obtain ⟨a, ha⟩ := hd1
    obtain ⟨b, hb⟩ := hd2
    exact ⟨a + b, by rw [show N - t.2.2 = t.1 + t.2.1 from by omega, ha, hb, Nat.mul_add]⟩)
    (h t.2.2 h3 hH (by omega))

private lemma sum_inv_sq_aux (E : ℕ) (hE : 0 < E) (M : ℕ) :
    ∑ e ∈ Finset.Ioc E M, (1 : ℝ) / (e : ℝ) ^ 2
      ≤ 1 / (E : ℝ) - 1 / ((max E M : ℕ) : ℝ) := by
  induction M with
  | zero =>
      have hempty : Finset.Ioc E 0 = (∅ : Finset ℕ) := Finset.Ioc_eq_empty (by omega)
      rw [hempty, Finset.sum_empty, show max E 0 = E from by omega]
      simp
  | succ M ih =>
      rcases lt_or_ge M E with h | h
      · have hempty : Finset.Ioc E (M + 1) = (∅ : Finset ℕ) := Finset.Ioc_eq_empty (by omega)
        rw [hempty, Finset.sum_empty, show max E (M + 1) = E from by omega]
        simp
      · rw [Finset.sum_Ioc_succ_top h, show max E (M + 1) = M + 1 from by omega]
        rw [show max E M = M from by omega] at ih
        have hM1 : (1 : ℝ) ≤ (M : ℝ) := by
          have : 1 ≤ M := by omega
          exact_mod_cast this
        have hkey : (1 : ℝ) / ((M + 1 : ℕ) : ℝ) ^ 2
            ≤ 1 / (M : ℝ) - 1 / ((M + 1 : ℕ) : ℝ) := by
          have hc : ((M + 1 : ℕ) : ℝ) = (M : ℝ) + 1 := by push_cast; ring
          rw [hc]
          have h1 : (0 : ℝ) < (M : ℝ) := by linarith
          have h2 : (0 : ℝ) < (M : ℝ) + 1 := by linarith
          have heq : 1 / (M : ℝ) - 1 / ((M : ℝ) + 1) = 1 / ((M : ℝ) * ((M : ℝ) + 1)) := by
            rw [div_sub_div _ _ (ne_of_gt h1) (ne_of_gt h2),
              show (1 : ℝ) * ((M : ℝ) + 1) - (M : ℝ) * 1 = 1 from by ring]
          rw [heq]
          refine one_div_le_one_div_of_le (by positivity) ?_
          nlinarith
        linarith

/-- `∑_{E < e ≤ M} 1/e² ≤ 1/E` for `E ≥ 1`, the first finite-sum inequality of (14). -/
theorem sum_inv_sq_le (E M : ℕ) (hE : 0 < E) :
    ∑ e ∈ Finset.Ioc E M, (1 : ℝ) / (e : ℝ) ^ 2 ≤ 1 / (E : ℝ) := by
  have h := sum_inv_sq_aux E hE M
  have hnn : (0 : ℝ) ≤ 1 / ((max E M : ℕ) : ℝ) := by positivity
  linarith

private lemma sum_inv_le_log_aux (H : ℕ) (hH : 0 < H) (M : ℕ) :
    ∑ e ∈ Finset.Ioc H M, (1 : ℝ) / (e : ℝ)
      ≤ Real.log ((max H M : ℕ) : ℝ) - Real.log (H : ℝ) := by
  induction M with
  | zero =>
      have hempty : Finset.Ioc H 0 = (∅ : Finset ℕ) := Finset.Ioc_eq_empty (by omega)
      rw [hempty, Finset.sum_empty, show max H 0 = H from by omega]
      simp
  | succ M ih =>
      rcases lt_or_ge M H with h | h
      · have hempty : Finset.Ioc H (M + 1) = (∅ : Finset ℕ) := Finset.Ioc_eq_empty (by omega)
        rw [hempty, Finset.sum_empty, show max H (M + 1) = H from by omega]
        simp
      · rw [Finset.sum_Ioc_succ_top h, show max H (M + 1) = M + 1 from by omega]
        rw [show max H M = M from by omega] at ih
        have hM1 : (1 : ℝ) ≤ (M : ℝ) := by
          have : 1 ≤ M := by omega
          exact_mod_cast this
        have hkey : (1 : ℝ) / ((M + 1 : ℕ) : ℝ)
            ≤ Real.log ((M + 1 : ℕ) : ℝ) - Real.log (M : ℝ) := by
          have hc : ((M + 1 : ℕ) : ℝ) = (M : ℝ) + 1 := by push_cast; ring
          rw [hc]
          have h1 : (0 : ℝ) < (M : ℝ) := by linarith
          have h2 : (0 : ℝ) < (M : ℝ) + 1 := by linarith
          have hlog : Real.log ((M : ℝ) / ((M : ℝ) + 1)) ≤ (M : ℝ) / ((M : ℝ) + 1) - 1 :=
            Real.log_le_sub_one_of_pos (by positivity)
          rw [Real.log_div (ne_of_gt h1) (ne_of_gt h2)] at hlog
          have hid : (M : ℝ) / ((M : ℝ) + 1) = 1 - 1 / ((M : ℝ) + 1) := by
            rw [eq_sub_iff_add_eq, ← add_div, div_self (ne_of_gt h2)]
          linarith
        linarith

/-- `∑_{H < e ≤ M} 1/e ≤ log(M/H)` for `1 ≤ H ≤ M`, the second finite-sum inequality
of (14). -/
theorem sum_inv_le_log_div {H M : ℕ} (hH : 0 < H) (hHM : H ≤ M) :
    ∑ e ∈ Finset.Ioc H M, (1 : ℝ) / (e : ℝ) ≤ Real.log ((M : ℝ) / (H : ℝ)) := by
  have h := sum_inv_le_log_aux H hH M
  rw [show max H M = M from by omega] at h
  have hM0 : (M : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have hH0 : (H : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  rw [Real.log_div hM0 hH0]
  exact h

/-- The three-range split (14):
`∑_{E < e ≤ N} |W(e)| ≤ 2NH/E + N log(X_N/H) + N`. -/
theorem sum_norm_W_le {u v : ℕ → ℂ} (hu : ∀ n, ‖u n‖ ≤ 1) (hv : ∀ n, ‖v n‖ ≤ 1)
    {N H E : ℕ} (hN : 3 ≤ N) (hH : 2 ≤ H) (hHN : H ≤ N) (hE : 0 < E) :
    ∑ e ∈ Finset.Ioc E N, ‖W u v N H e‖ ≤
      2 * (N : ℝ) * (H : ℝ) / (E : ℝ) +
        (N : ℝ) * Real.log ((XN N H : ℝ) / (H : ℝ)) + (N : ℝ) := by
  classical
  have hHX : H ≤ XN N H := le_XN hN hH hHN
  have hXpos : 0 < XN N H := by omega
  have hNr : (0 : ℝ) < (N : ℝ) := by
    have : 0 < N := by omega
    exact_mod_cast this
  have hEr : (0 : ℝ) < (E : ℝ) := by exact_mod_cast hE
  have hXr : (0 : ℝ) < (XN N H : ℝ) := by exact_mod_cast hXpos
  have hbound_inv : ∀ e : ℕ, H < e → 0 < e → ‖W u v N H e‖ ≤ (N : ℝ) * (1 / (e : ℝ)) := by
    intro e heH hepos
    have her : (0 : ℝ) < (e : ℝ) := by exact_mod_cast hepos
    have h3 : ‖W u v N H e‖ ≤ (((H / e + 1) * (N / e) : ℕ) : ℝ) := norm_W_le hu hv N H hepos
    rw [Nat.div_eq_of_lt heH] at h3
    have hsimp : (((0 + 1) * (N / e) : ℕ) : ℝ) = ((N / e : ℕ) : ℝ) := by push_cast; ring
    rw [hsimp] at h3
    have h2 : ((N / e : ℕ) : ℝ) ≤ (N : ℝ) / (e : ℝ) := Nat.cast_div_le
    calc ‖W u v N H e‖ ≤ ((N / e : ℕ) : ℝ) := h3
      _ ≤ (N : ℝ) / (e : ℝ) := h2
      _ = (N : ℝ) * (1 / (e : ℝ)) := by ring
  have hbound_sq : ∀ e ∈ (Finset.Ioc E N).filter (fun e => e ≤ H),
      ‖W u v N H e‖ ≤ 2 * (N : ℝ) * (H : ℝ) * (1 / (e : ℝ) ^ 2) := by
    intro e he
    rw [Finset.mem_filter, Finset.mem_Ioc] at he
    obtain ⟨⟨hEe, heN⟩, heH⟩ := he
    have hepos : 0 < e := by omega
    have her : (0 : ℝ) < (e : ℝ) := by exact_mod_cast hepos
    have hge : 1 ≤ H / e := (Nat.one_le_div_iff hepos).mpr heH
    have hle : ((H / e : ℕ) : ℝ) ≤ (H : ℝ) / (e : ℝ) := Nat.cast_div_le
    have hone : (1 : ℝ) ≤ ((H / e : ℕ) : ℝ) := by exact_mod_cast hge
    have h1 : ((H / e : ℕ) : ℝ) + 1 ≤ 2 * (H : ℝ) / (e : ℝ) := by
      have hr : 2 * (H : ℝ) / (e : ℝ) = (H : ℝ) / (e : ℝ) + (H : ℝ) / (e : ℝ) := by ring
      rw [hr]
      linarith
    have h2 : ((N / e : ℕ) : ℝ) ≤ (N : ℝ) / (e : ℝ) := Nat.cast_div_le
    have h3 : ‖W u v N H e‖ ≤ (((H / e + 1) * (N / e) : ℕ) : ℝ) := norm_W_le hu hv N H hepos
    have h4 : (((H / e + 1) * (N / e) : ℕ) : ℝ) ≤
        (2 * (H : ℝ) / (e : ℝ)) * ((N : ℝ) / (e : ℝ)) := by
      push_cast
      refine mul_le_mul ?_ h2 (by positivity) (by positivity)
      linarith
    have h5 : 2 * (N : ℝ) * (H : ℝ) * (1 / (e : ℝ) ^ 2)
        = (2 * (H : ℝ) / (e : ℝ)) * ((N : ℝ) / (e : ℝ)) := by
      rw [div_mul_div_comm, mul_one_div, pow_two]
      ring
    rw [h5]
    linarith
  have hpart1 : (∑ e ∈ (Finset.Ioc E N).filter (fun e => e ≤ H), ‖W u v N H e‖) ≤
      2 * (N : ℝ) * (H : ℝ) / (E : ℝ) := by
    calc (∑ e ∈ (Finset.Ioc E N).filter (fun e => e ≤ H), ‖W u v N H e‖)
        ≤ ∑ e ∈ (Finset.Ioc E N).filter (fun e => e ≤ H),
            2 * (N : ℝ) * (H : ℝ) * (1 / (e : ℝ) ^ 2) := Finset.sum_le_sum hbound_sq
      _ ≤ ∑ e ∈ Finset.Ioc E N, 2 * (N : ℝ) * (H : ℝ) * (1 / (e : ℝ) ^ 2) := by
          refine Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _) ?_
          intro e _ _
          positivity
      _ = 2 * (N : ℝ) * (H : ℝ) * ∑ e ∈ Finset.Ioc E N, (1 : ℝ) / (e : ℝ) ^ 2 := by
          rw [Finset.mul_sum]
      _ ≤ 2 * (N : ℝ) * (H : ℝ) * (1 / (E : ℝ)) :=
          mul_le_mul_of_nonneg_left (sum_inv_sq_le E N hE) (by positivity)
      _ = 2 * (N : ℝ) * (H : ℝ) / (E : ℝ) := by ring
  have hpart2 :
      (∑ e ∈ ((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter (fun e => e ≤ XN N H),
        ‖W u v N H e‖) ≤ (N : ℝ) * Real.log ((XN N H : ℝ) / (H : ℝ)) := by
    calc (∑ e ∈ ((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter (fun e => e ≤ XN N H),
          ‖W u v N H e‖)
        ≤ ∑ e ∈ ((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter (fun e => e ≤ XN N H),
            (N : ℝ) * (1 / (e : ℝ)) := by
          refine Finset.sum_le_sum ?_
          intro e he
          rw [Finset.mem_filter, Finset.mem_filter, Finset.mem_Ioc] at he
          exact hbound_inv e (by omega) (by omega)
      _ ≤ ∑ e ∈ Finset.Ioc H (XN N H), (N : ℝ) * (1 / (e : ℝ)) := by
          refine Finset.sum_le_sum_of_subset_of_nonneg ?_ ?_
          · intro e he
            rw [Finset.mem_filter, Finset.mem_filter, Finset.mem_Ioc] at he
            exact Finset.mem_Ioc.mpr ⟨by omega, he.2⟩
          · intro e he _
            rw [Finset.mem_Ioc] at he
            have : (0 : ℝ) < (e : ℝ) := by
              have : 0 < e := by omega
              exact_mod_cast this
            positivity
      _ = (N : ℝ) * ∑ e ∈ Finset.Ioc H (XN N H), (1 : ℝ) / (e : ℝ) := by rw [Finset.mul_sum]
      _ ≤ (N : ℝ) * Real.log ((XN N H : ℝ) / (H : ℝ)) :=
          mul_le_mul_of_nonneg_left (sum_inv_le_log_div (by omega) hHX) hNr.le
  have hpart3 :
      (∑ e ∈ ((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter (fun e => ¬ e ≤ XN N H),
        ‖W u v N H e‖) ≤ (N : ℝ) := by
    have hdrop :
        (∑ e ∈ (((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter
            (fun e => ¬ e ≤ XN N H)).filter (fun e => W u v N H e ≠ 0), ‖W u v N H e‖) =
          ∑ e ∈ ((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter (fun e => ¬ e ≤ XN N H),
            ‖W u v N H e‖ := by
      refine Finset.sum_filter_of_ne ?_
      intro x _ hx hw
      rw [hw] at hx
      simp at hx
    have hTcard : ((((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter
        (fun e => ¬ e ≤ XN N H)).filter (fun e => W u v N H e ≠ 0)).card ≤ XN N H := by
      have hsub : (((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter
          (fun e => ¬ e ≤ XN N H)).filter (fun e => W u v N H e ≠ 0) ⊆
            (Finset.Icc 1 H).biUnion (fun n₃ => (N - n₃).divisors) := by
        intro e he
        rw [Finset.mem_filter] at he
        by_contra hcon
        refine he.2 (W_eq_zero_of_forall_not_dvd ?_)
        intro n₃ hn1 hn2 hn3 hdvd
        exact hcon (Finset.mem_biUnion.mpr
          ⟨n₃, Finset.mem_Icc.mpr ⟨hn1, hn2⟩, Nat.mem_divisors.mpr ⟨hdvd, by omega⟩⟩)
      calc ((((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter
              (fun e => ¬ e ≤ XN N H)).filter (fun e => W u v N H e ≠ 0)).card
          ≤ ((Finset.Icc 1 H).biUnion (fun n₃ => (N - n₃).divisors)).card :=
            Finset.card_le_card hsub
        _ ≤ ∑ n₃ ∈ Finset.Icc 1 H, (N - n₃).divisors.card := Finset.card_biUnion_le
        _ = XN N H := (XN_def N H).symm
    have hterm : ∀ e ∈ (((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter
        (fun e => ¬ e ≤ XN N H)).filter (fun e => W u v N H e ≠ 0),
          ‖W u v N H e‖ ≤ (N : ℝ) / (XN N H : ℝ) := by
      intro e he
      rw [Finset.mem_filter, Finset.mem_filter, Finset.mem_filter, Finset.mem_Ioc] at he
      obtain ⟨⟨⟨⟨hEe, heN⟩, hH'⟩, hX'⟩, -⟩ := he
      have hepos : 0 < e := by omega
      have hXer : (XN N H : ℝ) ≤ (e : ℝ) := by
        have : XN N H ≤ e := by omega
        exact_mod_cast this
      calc ‖W u v N H e‖ ≤ (N : ℝ) * (1 / (e : ℝ)) := hbound_inv e (by omega) hepos
        _ ≤ (N : ℝ) * (1 / (XN N H : ℝ)) :=
            mul_le_mul_of_nonneg_left (one_div_le_one_div_of_le hXr hXer) hNr.le
        _ = (N : ℝ) / (XN N H : ℝ) := by ring
    rw [← hdrop]
    calc (∑ e ∈ (((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter
            (fun e => ¬ e ≤ XN N H)).filter (fun e => W u v N H e ≠ 0), ‖W u v N H e‖)
        ≤ ∑ _e ∈ (((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter
            (fun e => ¬ e ≤ XN N H)).filter (fun e => W u v N H e ≠ 0),
              (N : ℝ) / (XN N H : ℝ) := Finset.sum_le_sum hterm
      _ = (((((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter
            (fun e => ¬ e ≤ XN N H)).filter (fun e => W u v N H e ≠ 0)).card : ℝ) *
              ((N : ℝ) / (XN N H : ℝ)) := by rw [Finset.sum_const, nsmul_eq_mul]
      _ ≤ (XN N H : ℝ) * ((N : ℝ) / (XN N H : ℝ)) :=
          mul_le_mul_of_nonneg_right (by exact_mod_cast hTcard) (by positivity)
      _ = (N : ℝ) := by field_simp
  have hsplit :
      ∑ e ∈ Finset.Ioc E N, ‖W u v N H e‖ =
        (∑ e ∈ (Finset.Ioc E N).filter (fun e => e ≤ H), ‖W u v N H e‖) +
          ((∑ e ∈ ((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter
              (fun e => e ≤ XN N H), ‖W u v N H e‖) +
            ∑ e ∈ ((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)).filter
              (fun e => ¬ e ≤ XN N H), ‖W u v N H e‖) := by
    rw [Finset.sum_filter_add_sum_filter_not
        ((Finset.Ioc E N).filter (fun e => ¬ e ≤ H)) (fun e => e ≤ XN N H),
      Finset.sum_filter_add_sum_filter_not (Finset.Ioc E N) (fun e => e ≤ H)]
  rw [hsplit]
  linarith

/-! ## 5. The canonical divisor splitting

The identity (15) of Lemma 10, as the finite counting statement it is: for squarefree `d` the
number of admissible factorisations `d = d₁d₂` is `1` exactly when `d ∣ n₁n₂`, and `0`
otherwise.  The factorisation itself is `d₁ = (d, n₁)`.
-/

/-- The canonical splitting (15): for `d` squarefree,
`[d ∣ n₁n₂] = ∑_{d₁d₂ = d, (d₁,d₂) = 1} [d₁ ∣ n₁][d₂ ∣ n₂][(n₁,d₂) = 1]`, the right side
counting exactly one pair when `d ∣ n₁n₂` and none otherwise.  The paper applies this with
`n₁ ≥ 1`, which the identity itself does not need. -/
theorem canonical_split {d n₁ n₂ : ℕ} (hd : Squarefree d) :
    ((d.divisorsAntidiagonal).filter
        (fun p => Nat.Coprime p.1 p.2 ∧ p.1 ∣ n₁ ∧ p.2 ∣ n₂ ∧ Nat.Coprime n₁ p.2)).card =
      if d ∣ n₁ * n₂ then 1 else 0 := by
  have hd0 : d ≠ 0 := hd.ne_zero
  by_cases hdvd : d ∣ n₁ * n₂
  · rw [if_pos hdvd, Finset.card_eq_one]
    have hgd : Nat.gcd d n₁ ∣ d := Nat.gcd_dvd_left d n₁
    have hgn : Nat.gcd d n₁ ∣ n₁ := Nat.gcd_dvd_right d n₁
    have hprod : Nat.gcd d n₁ * (d / Nat.gcd d n₁) = d := Nat.mul_div_cancel' hgd
    have hcopg : Nat.Coprime (Nat.gcd d n₁) (d / Nat.gcd d n₁) :=
      Nat.coprime_of_squarefree_mul (by rw [hprod]; exact hd)
    have hDd : d / Nat.gcd d n₁ ∣ d := Nat.div_dvd_of_dvd hgd
    have hcopD : Nat.Coprime (d / Nat.gcd d n₁) n₁ := by
      have hk : Nat.gcd (d / Nat.gcd d n₁) n₁ ∣ 1 := by
        rw [← hcopg]
        exact Nat.dvd_gcd
          (Nat.dvd_gcd ((Nat.gcd_dvd_left _ _).trans hDd) (Nat.gcd_dvd_right _ _))
          (Nat.gcd_dvd_left _ _)
      exact Nat.dvd_one.mp hk
    have hDn₂ : d / Nat.gcd d n₁ ∣ n₂ := hcopD.dvd_of_dvd_mul_left (hDd.trans hdvd)
    refine ⟨(Nat.gcd d n₁, d / Nat.gcd d n₁), ?_⟩
    ext p
    simp only [Finset.mem_filter, Nat.mem_divisorsAntidiagonal, Finset.mem_singleton]
    constructor
    · rintro ⟨⟨hpd, -⟩, -, h1, -, h3⟩
      have hp1d : p.1 ∣ d := Dvd.intro p.2 hpd
      have hp1g : p.1 ∣ Nat.gcd d n₁ := Nat.dvd_gcd hp1d h1
      have hgp2 : Nat.Coprime (Nat.gcd d n₁) p.2 := Nat.Coprime.coprime_dvd_left hgn h3
      have hgp1 : Nat.gcd d n₁ ∣ p.1 := hgp2.dvd_of_dvd_mul_right (by rw [hpd]; exact hgd)
      have he1 : p.1 = Nat.gcd d n₁ := Nat.dvd_antisymm hp1g hgp1
      have hp10 : 0 < p.1 := by
        rcases Nat.eq_zero_or_pos p.1 with h0 | h0
        · rw [h0, zero_mul] at hpd
          exact absurd hpd.symm hd0
        · exact h0
      have he2 : p.2 = d / Nat.gcd d n₁ := by
        rw [← he1, ← hpd, Nat.mul_div_cancel_left _ hp10]
      exact Prod.ext_iff.mpr ⟨he1, he2⟩
    · rintro rfl
      exact ⟨⟨hprod, hd0⟩, hcopg, hgn, hDn₂, hcopD.symm⟩
  · rw [if_neg hdvd, Finset.card_eq_zero, Finset.filter_eq_empty_iff]
    intro p hp hcon
    obtain ⟨-, h1, h2, -⟩ := hcon
    rw [Nat.mem_divisorsAntidiagonal] at hp
    exact hdvd (hp.1 ▸ mul_dvd_mul h1 h2)

end ShortVariableMobius

#print axioms ShortVariableMobius.moebius_mul_mul
#print axioms ShortVariableMobius.card_le_of_const_mod
#print axioms ShortVariableMobius.card_le_of_two_classes
#print axioms ShortVariableMobius.mem_triples
#print axioms ShortVariableMobius.S_eq_sum_moebius_mul_W
#print axioms ShortVariableMobius.card_filter_triples_le
#print axioms ShortVariableMobius.norm_W_le
#print axioms ShortVariableMobius.XN_le_mul_taustar
#print axioms ShortVariableMobius.le_XN
#print axioms ShortVariableMobius.W_eq_zero_of_forall_not_dvd
#print axioms ShortVariableMobius.sum_inv_sq_le
#print axioms ShortVariableMobius.sum_inv_le_log_div
#print axioms ShortVariableMobius.sum_norm_W_le
#print axioms ShortVariableMobius.canonical_split

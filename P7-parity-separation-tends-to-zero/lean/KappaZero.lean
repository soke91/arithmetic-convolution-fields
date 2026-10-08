/-
  KappaZero.lean --- the exact, finite, combinatorial steps of the proof that
  kappa(N) -> 0.

  WHAT THIS COVERS, in the paper's numbering (see its Formalization section).  Lemma 7(i) and
  the upper bound of Lemma 7(ii); the exactness on Y-smooth rows with the row identity (18) of
  Lemma 10; the large-prime removal at d = 1, which is Lemma 22; the Bonferroni identity used
  in Lemma 15 and the exact size of its truncation error; the parity split and the ratio bound
  (24); and the containment Ax in K_Q for x >= 0, the cone (5).

  WHAT IS HERE AND WHAT IS NOT.  The paper's analytic estimates -- Mertens
  (M1)-(M3), the effective prime number theorem (P), the effective Moebius bound (L), Rankin's
  trick (Lemma 14), the sieve asymptotic (the error analysis of Lemma 15), the friable
  bound (Lemma 17), the rough Moebius sum (Lemma 19) -- are NOT formalized and are
  not attempted here.  What is formalized is the part of the proof that is a
  finite identity or a finite inequality: the one-large-prime combinatorics, the
  exactness of the construction on smooth rows, the Bonferroni identity, the
  parity split, and the cone step.

  NO `sorry` APPEARS IN THIS FILE.  `#print axioms` at the end prints the axiom
  set of every theorem; `axioms.txt` is that output.

  Elaborate with the toolchain `leanprover/lean4:v4.33.1` and a mathlib of the
  same revision.  From this directory, with `lean` and `lake` on PATH:

      lake env lean KappaZero.lean

  `lakefile.toml`, `lean-toolchain` and `lake-manifest.json` sit beside this file; `lake` will
  fetch and build mathlib on the first run if no build is present.

  A NOTE ON `Nat.ArithmeticFunction`.  In this mathlib
  (v4.33.1) the namespace is `ArithmeticFunction`, not `Nat.ArithmeticFunction`
  --- `Nat.ArithmeticFunction.moebius` is an unknown constant.  `mu` is
  `ArithmeticFunction.moebius` and `omega` is
  `ArithmeticFunction.cardDistinctFactors`.
-/

import Mathlib

/-!
# The exact finite parts of `κ(N) → 0`

Five groups of statements.

1. **One large prime** (Lemma 7(i) and the upper bound of Lemma 7(ii) of the paper):
   `card_filter_lt_primeFactors_le_one`, `not_sq_dvd_of_lt`, `thr_le_of_lt`,
   and the exact identity `one_large_prime_identity`.
2. **Exactness on smooth rows** (Lemma 10 of the paper): `dvd_smallPart_iff` and
   the identity `identity_E`.
3. **Bonferroni** (the identity used in Lemma 15 of the paper): `bonferroni`,
   `bonferroni_error_eq`,
   `bonferroni_error_le`.
4. **Parity split** (the split and the ratio bound (24) of the paper): `parity_split`,
   `par_ratio`.
5. **The cone step** (`Ax` in the cone (5) of the paper for `x >= 0`): `infDist_coneHull_le`.
-/

namespace KappaZero

open Finset ArithmeticFunction
open scoped ArithmeticFunction.Moebius

/-! ## §0 Smoothness -/

/-- `P⁺(m) ≤ Q`: every prime factor of `m` is at most `Q`.  `Smooth Q 0` and
`Smooth Q 1` both hold, since `primeFactors 0 = primeFactors 1 = ∅`; that is the
convention `P⁺(1) = 1`, and `0` never enters a sum below. -/
def Smooth (Q m : ℕ) : Prop := ∀ p ∈ m.primeFactors, p ≤ Q

instance decidableSmooth (Q m : ℕ) : Decidable (Smooth Q m) :=
  inferInstanceAs (Decidable (∀ p ∈ m.primeFactors, p ≤ Q))

/-! ## §1 Lemma 1 --- at most one prime factor above `Q`

Lemma 7 of the paper is the reason the theorem holds for **every** `N` and not
only on the ladder `N = 2^e`.  Only the upper bound `N < (Q+1)^2` is used; the
lower bound `Q^2 ≤ N` of `Q = ⌊√N⌋` is recorded in `sqrt_bounds` and is not
needed by any statement in this file. -/

/-- `Q = ⌊√N⌋` gives `Q^2 ≤ N < (Q+1)^2`, so every lemma below applies to it. -/
theorem sqrt_bounds (N : ℕ) : N.sqrt ^ 2 ≤ N ∧ N < (N.sqrt + 1) ^ 2 :=
  ⟨Nat.sqrt_le' N, Nat.lt_succ_sqrt' N⟩

/-- **Lemma 7(i) of the paper, first half.**  If `N < (Q+1)^2` then no `n ≤ N` has
two distinct prime factors `> Q`.  (Stated for every `n ≤ N`, not only the
squarefree ones; squarefreeness is not used.) -/
theorem card_filter_lt_primeFactors_le_one {N Q n : ℕ} (hN : N < (Q + 1) ^ 2)
    (hn : n ≤ N) : (n.primeFactors.filter (fun p => Q < p)).card ≤ 1 := by
  rw [Finset.card_le_one]
  intro p hp q hq
  by_contra hne
  simp only [Finset.mem_filter, Nat.mem_primeFactors] at hp hq
  obtain ⟨⟨hpp, hpd, hn0⟩, hpQ⟩ := hp
  obtain ⟨⟨hqp, hqd, -⟩, hqQ⟩ := hq
  have hco : Nat.Coprime p q := (Nat.coprime_primes hpp hqp).mpr hne
  have hdvd : p * q ∣ n := hco.mul_dvd_of_dvd_of_dvd hpd hqd
  have hle : p * q ≤ n := Nat.le_of_dvd (Nat.pos_of_ne_zero hn0) hdvd
  have h1 : Q + 1 ≤ p := by omega
  have h2 : Q + 1 ≤ q := by omega
  have hbig : (Q + 1) * (Q + 1) ≤ N :=
    le_trans (Nat.mul_le_mul h1 h2) (le_trans hle hn)
  rw [pow_two] at hN
  exact absurd hN (Nat.not_lt.mpr hbig)

/-- **Lemma 7(i) of the paper, second half.**  A prime above `Q` divides `n ≤ N` at
most once. -/
theorem not_sq_dvd_of_lt {N Q n p : ℕ} (hN : N < (Q + 1) ^ 2) (hn : n ≤ N)
    (hn0 : n ≠ 0) (hpQ : Q < p) : ¬ p * p ∣ n := by
  intro h
  have hle : p * p ≤ n := Nat.le_of_dvd (Nat.pos_of_ne_zero hn0) h
  have h1 : Q + 1 ≤ p := by omega
  have hbig : (Q + 1) * (Q + 1) ≤ N :=
    le_trans (Nat.mul_le_mul h1 h1) (le_trans hle hn)
  rw [pow_two] at hN
  exact absurd hN (Nat.not_lt.mpr hbig)

/-- **The upper bound of Lemma 7(ii) of the paper.**  `thr = ⌊N/q⌋ ≤ Q` for every
`q > Q`; in particular for `q = q₁`, the least prime `> Q`. -/
theorem thr_le_of_lt {N Q q : ℕ} (hN : N < (Q + 1) ^ 2) (hq : Q < q) :
    N / q ≤ Q := by
  have hq0 : 0 < q := by omega
  by_contra h
  have h1 : Q + 1 ≤ N / q := by omega
  rw [Nat.le_div_iff_mul_le hq0] at h1
  have h2 : (Q + 1) * (Q + 1) ≤ (Q + 1) * q := Nat.mul_le_mul le_rfl (by omega)
  rw [pow_two] at hN
  exact absurd hN (Nat.not_lt.mpr (le_trans h2 h1))

/-! ## §2 The small-prime part, and exactness on smooth rows

The small-prime part is `a_Y(n) = ∏_{p|n, p≤Y} p`.  `dvd_smallPart_iff` is the
one line the proof of Lemma 10 opens with --- "every prime factor of `d` is `≤ Y` and
`d` is squarefree, so for squarefree `ℓ`, `d | ℓ ⟺ d | a_Y(ℓ)`" --- and it is
what makes the construction's row sums exactly computable. -/

/-- `a_Y(n) = ∏_{p ∣ n, p ≤ Y} p`, the small-prime part. -/
def smallPart (Y n : ℕ) : ℕ := ∏ p ∈ n.primeFactors.filter (fun p => p ≤ Y), p

/-- `a_Y(n) ∣ n` for squarefree `n`. -/
theorem smallPart_dvd {Y n : ℕ} (hn : Squarefree n) : smallPart Y n ∣ n := by
  have h : (∏ p ∈ n.primeFactors.filter (fun p => p ≤ Y), p) ∣ ∏ p ∈ n.primeFactors, p :=
    Finset.prod_dvd_prod_of_subset _ _ (fun p => p) (Finset.filter_subset _ _)
  rwa [Nat.prod_primeFactors_of_squarefree hn] at h

/-- **Lemma 10 of the paper, the exactness step.**  For squarefree `d` with every
prime factor `≤ Y` and squarefree `n`:  `d ∣ n ↔ d ∣ a_Y(n)`. -/
theorem dvd_smallPart_iff {Y d n : ℕ} (hd : Squarefree d) (hdY : Smooth Y d)
    (hn : Squarefree n) : d ∣ n ↔ d ∣ smallPart Y n := by
  refine ⟨fun h => ?_, fun h => h.trans (smallPart_dvd hn)⟩
  have hsub : d.primeFactors ⊆ n.primeFactors.filter (fun p => p ≤ Y) := by
    intro p hp
    rw [Finset.mem_filter]
    exact ⟨Nat.primeFactors_mono h hn.ne_zero hp, hdY p hp⟩
  have hdvd := Finset.prod_dvd_prod_of_subset _ _ (fun p => p) hsub
  rw [Nat.prod_primeFactors_of_squarefree hd] at hdvd
  exact hdvd

/-! ## §3 The construction, and the row identity (18)

The construction (10) of the paper and the row identity (18) of its Lemma 10.
Everything is stated after multiplication by `R`, so that `t_d` is the integer
`c_R(d)` and (18) is an identity between a rational row sum and two integer
class sums; dividing by `R` is the content-free step.

The band itself is a parameter: `vB` and `vR` are arbitrary finite sets of
squarefree integers and `cls` any finite set carrying the small-prime parts
that occur.  Nothing below needs them to be the band (2) of the paper, which is
exactly why (18) is "unconditional" there. -/

/-- `R_a = #{n ∈ vR : a_Y(n) = a}`. -/
def Rc (Y : ℕ) (vR : Finset ℕ) (a : ℕ) : ℕ :=
  (vR.filter (fun n => smallPart Y n = a)).card

/-- `B_a = #{ℓ ∈ vB : a_Y(ℓ) = a}`. -/
def Bc (Y : ℕ) (vB : Finset ℕ) (a : ℕ) : ℕ :=
  (vB.filter (fun n => smallPart Y n = a)).card

/-- The weights (10) of the paper, scaled by `R`:
`R·x_ℓ = R_{a_Y(ℓ)}/B_{a_Y(ℓ)}` when `a_Y(ℓ) ≤ A` and `B_{a_Y(ℓ)} > 0`, else `0`.
The second clause is the guard against a dead retained class. -/
def wt (Y A : ℕ) (vB vR : Finset ℕ) (l : ℕ) : ℚ :=
  if smallPart Y l ≤ A ∧ 0 < Bc Y vB (smallPart Y l) then
    (Rc Y vR (smallPart Y l) : ℚ) / (Bc Y vB (smallPart Y l) : ℚ)
  else 0

/-- `R·(Ax)_d`, the `d`-th row sum of the construction. -/
def AxRow (Y A : ℕ) (vB vR : Finset ℕ) (d : ℕ) : ℚ :=
  ∑ l ∈ vB, wt Y A vB vR l * (if d ∣ l then 1 else 0)

/-- `R·t_d = c_R(d) = #{n ∈ vR : d ∣ n}`. -/
def tRow (vR : Finset ℕ) (d : ℕ) : ℕ := (vR.filter (fun n => d ∣ n)).card

/-- `vR` is partitioned by `a_Y`, so `c_R(d)` is the sum of `R_a` over the
classes `a` divisible by `d`.  This is the second display in the proof of Lemma 10. -/
theorem tRow_eq_sum {Y d : ℕ} {vR cls : Finset ℕ}
    (hsq : ∀ n ∈ vR, Squarefree n) (hcls : ∀ n ∈ vR, smallPart Y n ∈ cls)
    (hd : Squarefree d) (hdY : Smooth Y d) :
    tRow vR d = ∑ a ∈ cls.filter (fun a => d ∣ a), Rc Y vR a := by
  have hmaps : ∀ n ∈ vR.filter (fun n => d ∣ n),
      smallPart Y n ∈ cls.filter (fun a => d ∣ a) := by
    intro n hn
    simp only [Finset.mem_filter] at hn ⊢
    exact ⟨hcls n hn.1, (dvd_smallPart_iff hd hdY (hsq n hn.1)).mp hn.2⟩
  rw [tRow, Finset.card_eq_sum_card_fiberwise hmaps]
  refine Finset.sum_congr rfl ?_
  intro a ha
  simp only [Finset.mem_filter] at ha
  rw [Rc]
  congr 1
  ext n
  simp only [Finset.mem_filter]
  refine ⟨fun h => ⟨h.1.1, h.2⟩, fun h => ⟨⟨h.1, ?_⟩, h.2⟩⟩
  exact (dvd_smallPart_iff hd hdY (hsq n h.1)).mpr (h.2 ▸ ha.2)

/-- The row sum of the construction, in closed form: `R·(Ax)_d` is the sum of
`R_a` over the **retained** classes divisible by `d`.  This is the first display
in the proof of Lemma 10. -/
theorem AxRow_eq_sum {Y A d : ℕ} {vB vR cls : Finset ℕ}
    (hsq : ∀ l ∈ vB, Squarefree l) (hcls : ∀ l ∈ vB, smallPart Y l ∈ cls)
    (hd : Squarefree d) (hdY : Smooth Y d) :
    AxRow Y A vB vR d
      = ∑ a ∈ cls.filter (fun a => d ∣ a ∧ a ≤ A ∧ 0 < Bc Y vB a),
          (Rc Y vR a : ℚ) := by
  rw [AxRow, ← Finset.sum_fiberwise_of_maps_to hcls, Finset.sum_filter]
  refine Finset.sum_congr rfl ?_
  intro a _
  have key : ∀ l ∈ vB.filter (fun l => smallPart Y l = a),
      wt Y A vB vR l * (if d ∣ l then 1 else 0)
        = (if d ∣ a ∧ a ≤ A ∧ 0 < Bc Y vB a then
            (Rc Y vR a : ℚ) / (Bc Y vB a : ℚ) else 0) := by
    intro l hl
    simp only [Finset.mem_filter] at hl
    obtain ⟨hlB, hla⟩ := hl
    have hdl : (d ∣ l) ↔ (d ∣ a) := by
      rw [dvd_smallPart_iff hd hdY (hsq l hlB), hla]
    rw [wt, hla]
    by_cases hA : a ≤ A ∧ 0 < Bc Y vB a
    · by_cases hda : d ∣ a
      · rw [if_pos hA, if_pos (hdl.mpr hda), if_pos ⟨hda, hA⟩, mul_one]
      · have h1 : ¬ (d ∣ l) := fun h => hda (hdl.mp h)
        have h2 : ¬ (d ∣ a ∧ a ≤ A ∧ 0 < Bc Y vB a) := fun h => hda h.1
        rw [if_pos hA, if_neg h1, mul_zero, if_neg h2]
    · have h2 : ¬ (d ∣ a ∧ a ≤ A ∧ 0 < Bc Y vB a) := fun h => hA h.2
      rw [if_neg hA, zero_mul, if_neg h2]
  rw [Finset.sum_congr rfl key, Finset.sum_const,
    show (vB.filter (fun l => smallPart Y l = a)).card = Bc Y vB a from rfl,
    nsmul_eq_mul]
  by_cases hP : d ∣ a ∧ a ≤ A ∧ 0 < Bc Y vB a
  · rw [if_pos hP, if_pos hP]
    have hB : (Bc Y vB a : ℚ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
    field_simp
  · rw [if_neg hP, if_neg hP, mul_zero]

/-- **Lemma 10 of the paper, identity (18)**, scaled by `R`.  The row sum of the
construction equals `t_d` minus the discarded classes: those with `a > A`
(truncated away by (10)) and those with `a ≤ A` but `B_a = 0` (the dead
retained classes, which Hypothesis 5 of the paper excludes and which the measured construction
does exhibit at `(N = 2^22, Y = 443)`). -/
theorem identity_E {Y A d : ℕ} {vB vR cls : Finset ℕ}
    (hsqB : ∀ l ∈ vB, Squarefree l) (hsqR : ∀ n ∈ vR, Squarefree n)
    (hclsB : ∀ l ∈ vB, smallPart Y l ∈ cls) (hclsR : ∀ n ∈ vR, smallPart Y n ∈ cls)
    (hd : Squarefree d) (hdY : Smooth Y d) :
    (tRow vR d : ℚ) - AxRow Y A vB vR d
      = (∑ a ∈ cls.filter (fun a => d ∣ a ∧ A < a), (Rc Y vR a : ℚ))
        + ∑ a ∈ cls.filter (fun a => d ∣ a ∧ a ≤ A ∧ Bc Y vB a = 0),
            (Rc Y vR a : ℚ) := by
  have ht : (tRow vR d : ℚ) = ∑ a ∈ cls.filter (fun a => d ∣ a), (Rc Y vR a : ℚ) := by
    rw [tRow_eq_sum hsqR hclsR hd hdY]
    push_cast
    ring
  rw [ht, AxRow_eq_sum hsqB hclsB hd hdY, Finset.sum_filter, Finset.sum_filter,
    Finset.sum_filter, Finset.sum_filter, sub_eq_iff_eq_add,
    ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun a _ => ?_
  by_cases h1 : d ∣ a
  · by_cases h2 : a ≤ A
    · by_cases h3 : Bc Y vB a = 0
      · simp [h1, h2, h3]
      · simp [h1, h2, h3, Nat.pos_of_ne_zero h3, Nat.not_lt.mpr h2]
    · simp [h1, h2, Nat.lt_of_not_le h2]
  · simp [h1]

/-! ## §4 The one-large-prime identity

Lemma 22 of the paper is this combinatorial step, and the same step occurs at a second
cutoff.  The exact identity formalized here is the one at
`d = 1`:

    P₁ = Σ_{thr < n ≤ N, n squarefree, P⁺(n) ≤ Q} μ(n)
       = M(N) − M(thr) + Σ_{k ≤ thr} μ(k)·(π(N/k) − π(Q)),     thr = ⌊N/q₁⌋.

Every step is finite and exact: no estimate of `M` or `π` enters. -/

/-- `M(x) = Σ_{n ≤ x} μ(n)`. -/
def mert (x : ℕ) : ℤ := ∑ n ∈ Finset.Ioc 0 x, μ n

/-- The primes `≤ x`. -/
def primesLe (x : ℕ) : Finset ℕ := (Finset.Ioc 0 x).filter Nat.Prime

/-- `π(x) = #{p ≤ x : p prime}`. -/
def primeCount (x : ℕ) : ℕ := (primesLe x).card

/-- `primeCount` is mathlib's `Nat.primeCounting`.  Recorded so that the
identity below is a statement about the standard `π`, not about a private
definition. -/
theorem primeCount_eq (x : ℕ) : primeCount x = Nat.primeCounting x := by
  rw [primeCount, ← Nat.primesLE_card_eq_primeCounting]
  congr 1
  ext p
  simp only [primesLe, Nat.primesLE, Nat.primesBelow, Finset.mem_filter,
    Finset.mem_Ioc, Finset.mem_range, Nat.lt_succ_iff]
  exact ⟨fun h => ⟨h.1.2, h.2⟩, fun h => ⟨⟨h.2.pos, h.1⟩, h.2⟩⟩

/-- `A(x) = Σ_{m ≤ x, P⁺(m) ≤ Q} μ(m)`: the Moebius sum restricted to `Q`-smooth `m`. -/
def smoothMert (Q x : ℕ) : ℤ := ∑ m ∈ (Finset.Ioc 0 x).filter (Smooth Q), μ m

/-- `A(x) = M(x)` as soon as `x ≤ Q`: every `m ≤ x` is then automatically
`Q`-smooth.  With `x = thr` this is "`A(thr) = M(thr)` because `thr ≤ Q`
makes every `m ≤ thr` automatically `Q`-friable". -/
theorem smoothMert_eq_mert {Q x : ℕ} (hx : x ≤ Q) : smoothMert Q x = mert x := by
  rw [smoothMert, mert, Finset.filter_true_of_mem]
  intro m hm p hp
  exact le_trans (Nat.le_of_mem_primeFactors hp)
    (le_trans (Finset.mem_Ioc.mp hm).2 hx)

/-- The unique prime factor of `n` above `Q`, written as a product over a filter
that `card_filter_lt_primeFactors_le_one` shows has at most one element.  A total
definition needing no choice: it is `1` when there is no such prime. -/
def bigPrime (Q n : ℕ) : ℕ := ∏ p ∈ n.primeFactors.filter (fun p => Q < p), p

/-- If `p` is a prime `> Q` dividing `n ≤ N`, it is *the* one. -/
theorem bigPrime_eq {N Q n p : ℕ} (hN : N < (Q + 1) ^ 2) (hn : n ≤ N) (hn0 : n ≠ 0)
    (hp : p.Prime) (hpQ : Q < p) (hpd : p ∣ n) : bigPrime Q n = p := by
  have hmem : p ∈ n.primeFactors.filter (fun r => Q < r) :=
    Finset.mem_filter.mpr ⟨Nat.mem_primeFactors.mpr ⟨hp, hpd, hn0⟩, hpQ⟩
  have hsing : n.primeFactors.filter (fun r => Q < r) = {p} :=
    Finset.eq_singleton_iff_unique_mem.mpr ⟨hmem, fun r hr =>
      Finset.card_le_one.mp (card_filter_lt_primeFactors_le_one hN hn) r hr p hmem⟩
  rw [bigPrime, hsing, Finset.prod_singleton]

/-- A rough `n ≤ N` has a large prime factor, and `bigPrime` is it. -/
theorem bigPrime_spec {N Q n : ℕ} (hN : N < (Q + 1) ^ 2) (hn : n ≤ N) (hn0 : n ≠ 0)
    (hr : ¬ Smooth Q n) :
    (bigPrime Q n).Prime ∧ Q < bigPrime Q n ∧ bigPrime Q n ∣ n := by
  rw [Smooth] at hr
  push Not at hr
  obtain ⟨p, hp, hpQ⟩ := hr
  rw [Nat.mem_primeFactors] at hp
  rw [bigPrime_eq hN hn hn0 hp.1 (by omega) hp.2.1]
  exact ⟨hp.1, by omega, hp.2.1⟩

/-- The cofactor `n / P⁺(n)` of a squarefree rough `n ≤ N` is `Q`-smooth: the
upper cutoff on the cofactor is **automatic**, which is the load-bearing
observation behind Lemma 22 of the paper. -/
theorem smooth_div_bigPrime {N Q n : ℕ} (hN : N < (Q + 1) ^ 2) (hn : n ≤ N)
    (hsq : Squarefree n) (hr : ¬ Smooth Q n) : Smooth Q (n / bigPrime Q n) := by
  obtain ⟨hp, hpQ, hpd⟩ := bigPrime_spec hN hn hsq.ne_zero hr
  intro q hq
  by_contra hqQ
  rw [Nat.not_le] at hqQ
  rw [Nat.mem_primeFactors] at hq
  obtain ⟨hqp, hqd, -⟩ := hq
  have hdvdn : q ∣ n := hqd.trans (Nat.div_dvd_of_dvd hpd)
  have heq : q = bigPrime Q n :=
    (bigPrime_eq hN hn hsq.ne_zero hqp hqQ hdvdn).symm ▸ rfl
  obtain ⟨k, hk⟩ : bigPrime Q n ∣ n / bigPrime Q n := heq ▸ hqd
  have hsqd : bigPrime Q n * bigPrime Q n ∣ n := by
    refine ⟨k, ?_⟩
    calc n = n / bigPrime Q n * bigPrime Q n := (Nat.div_mul_cancel hpd).symm
      _ = bigPrime Q n * k * bigPrime Q n := by rw [hk]
      _ = bigPrime Q n * bigPrime Q n * k := by ring
  have h1 := hp.one_lt
  have h2 := Nat.isUnit_iff.mp (hsq _ hsqd)
  omega

/-- **Lemma 22 of the paper, the bijection step.**  `M(N) = A(N) − Σ_{Q<p≤N} A(N/p)`:
the squarefree `m ≤ N` that are not `Q`-smooth are in bijection with the pairs
`(p, m')`, `p` prime in `(Q, N]` and `m' ≤ N/p` squarefree and `Q`-smooth, and
`μ(p m') = −μ(m')`. -/
theorem mert_eq_smoothMert_sub {N Q : ℕ} (hN : N < (Q + 1) ^ 2) :
    mert N = smoothMert Q N
      - ∑ p ∈ (Finset.Ioc 0 N).filter (fun p => p.Prime ∧ Q < p),
          smoothMert Q (N / p) := by
  set T := (Finset.Ioc 0 N).filter (fun p => p.Prime ∧ Q < p) with hT
  -- drop the non-squarefree terms from a sum of `μ`
  have drop : ∀ (s : Finset ℕ),
      ∑ m ∈ s, μ m = ∑ m ∈ s.filter Squarefree, μ m := by
    intro s
    rw [Finset.sum_filter]
    refine Finset.sum_congr rfl fun m _ => ?_
    by_cases h : Squarefree m
    · rw [if_pos h]
    · rw [if_neg h, moebius_eq_zero_of_not_squarefree h]
  -- split `M(N)` into the smooth and the rough part
  have hsplit : mert N
      = smoothMert Q N + ∑ m ∈ (Finset.Ioc 0 N).filter (fun m => ¬ Smooth Q m), μ m := by
    rw [mert, smoothMert]
    exact (Finset.sum_filter_add_sum_filter_not _ _ _).symm
  -- every rough squarefree `m ≤ N` has its large prime in `T`
  have hmaps : ∀ m ∈ ((Finset.Ioc 0 N).filter (fun m => ¬ Smooth Q m)).filter Squarefree,
      bigPrime Q m ∈ T := by
    intro m hm
    simp only [Finset.mem_filter, Finset.mem_Ioc] at hm
    obtain ⟨⟨⟨hm0, hmN⟩, hrough⟩, hmsq⟩ := hm
    obtain ⟨hp, hpQ, hpd⟩ := bigPrime_spec hN hmN hmsq.ne_zero hrough
    simp only [hT, Finset.mem_filter, Finset.mem_Ioc]
    exact ⟨⟨hp.pos, le_trans (Nat.le_of_dvd (by omega) hpd) hmN⟩, hp, hpQ⟩
  -- the fibre over `p` is `−A(N/p)`
  have hfib : ∀ p ∈ T,
      ∑ m ∈ (((Finset.Ioc 0 N).filter (fun m => ¬ Smooth Q m)).filter Squarefree).filter
          (fun m => bigPrime Q m = p), μ m = - smoothMert Q (N / p) := by
    intro p hpT
    simp only [hT, Finset.mem_filter, Finset.mem_Ioc] at hpT
    obtain ⟨⟨hp0, hpN⟩, hpp, hpQ⟩ := hpT
    rw [smoothMert, drop ((Finset.Ioc 0 (N / p)).filter (Smooth Q)), ← Finset.sum_neg_distrib]
    refine Finset.sum_bij' (fun m _ => m / p) (fun m' _ => p * m') ?_ ?_ ?_ ?_ ?_
    · -- forward map lands in the smooth squarefree cofactors
      intro m hm
      simp only [Finset.mem_filter, Finset.mem_Ioc] at hm ⊢
      obtain ⟨⟨⟨⟨hm0, hmN⟩, hrough⟩, hmsq⟩, hbig⟩ := hm
      obtain ⟨-, -, hpd⟩ := bigPrime_spec hN hmN hmsq.ne_zero hrough
      rw [hbig] at hpd
      have hdpos : 0 < m / p := Nat.div_pos (Nat.le_of_dvd (by omega) hpd) hp0
      refine ⟨⟨⟨hdpos, Nat.div_le_div_right hmN⟩, ?_⟩, ?_⟩
      · exact hbig ▸ smooth_div_bigPrime hN hmN hmsq hrough
      · exact hmsq.squarefree_of_dvd (Nat.div_dvd_of_dvd hpd)
    · -- backward map lands in the fibre
      intro m' hm'
      simp only [Finset.mem_filter, Finset.mem_Ioc] at hm' ⊢
      obtain ⟨⟨⟨hm0, hmle⟩, hsm⟩, hsq⟩ := hm'
      have hpnd : ¬ p ∣ m' := by
        intro h
        exact absurd (hsm p (Nat.mem_primeFactors.mpr ⟨hpp, h, by omega⟩)) (by omega)
      have hco : Nat.Coprime p m' := (Nat.Prime.coprime_iff_not_dvd hpp).mpr hpnd
      have hleN : p * m' ≤ N := by
        rw [mul_comm]
        exact (Nat.le_div_iff_mul_le hp0).mp hmle
      have hpos : 0 < p * m' := Nat.mul_pos hp0 hm0
      have hpdvd : p ∣ p * m' := Dvd.intro m' rfl
      have hsqm : Squarefree (p * m') :=
        Nat.squarefree_mul_iff.mpr ⟨hco, hpp.squarefree, hsq⟩
      have hrough : ¬ Smooth Q (p * m') := by
        intro h
        exact absurd (h p (Nat.mem_primeFactors.mpr ⟨hpp, hpdvd, by omega⟩)) (by omega)
      exact ⟨⟨⟨⟨hpos, hleN⟩, hrough⟩, hsqm⟩,
        bigPrime_eq hN hleN (by omega) hpp hpQ hpdvd⟩
    · -- left inverse
      intro m hm
      simp only [Finset.mem_filter, Finset.mem_Ioc] at hm
      obtain ⟨⟨⟨⟨hm0, hmN⟩, hrough⟩, hmsq⟩, hbig⟩ := hm
      obtain ⟨-, -, hpd⟩ := bigPrime_spec hN hmN hmsq.ne_zero hrough
      rw [hbig] at hpd
      exact Nat.mul_div_cancel' hpd
    · -- right inverse
      intro m' _
      exact Nat.mul_div_cancel_left m' hp0
    · -- the values match
      intro m hm
      simp only [Finset.mem_filter, Finset.mem_Ioc] at hm
      obtain ⟨⟨⟨⟨hm0, hmN⟩, hrough⟩, hmsq⟩, hbig⟩ := hm
      obtain ⟨-, -, hpd⟩ := bigPrime_spec hN hmN hmsq.ne_zero hrough
      rw [hbig] at hpd
      have hsm : Smooth Q (m / p) := hbig ▸ smooth_div_bigPrime hN hmN hmsq hrough
      have hpnd : ¬ p ∣ m / p := by
        intro h
        exact absurd (hsm p (Nat.mem_primeFactors.mpr
          ⟨hpp, h, by have := Nat.div_pos (Nat.le_of_dvd (by omega) hpd) hp0; omega⟩))
          (by omega)
      have hco : Nat.Coprime p (m / p) := (Nat.Prime.coprime_iff_not_dvd hpp).mpr hpnd
      calc μ m = μ (p * (m / p)) := by rw [Nat.mul_div_cancel' hpd]
        _ = μ p * μ (m / p) := isMultiplicative_moebius.map_mul_of_coprime hco
        _ = -(μ (m / p)) := by rw [moebius_apply_prime hpp, neg_one_mul]
  rw [hsplit, drop ((Finset.Ioc 0 N).filter (fun m => ¬ Smooth Q m)),
    ← Finset.sum_fiberwise_of_maps_to hmaps, Finset.sum_congr rfl hfib,
    Finset.sum_neg_distrib]
  ring

/-- **Lemma 22 of the paper, second form.**  `A(N) = M(N) + Σ_{Q<p≤N} M(N/p)`:  for
`p > Q` the inner `A` is an unrestricted `M`, because `N/p ≤ Q` by Lemma 1(ii).
This is what makes the identity elementary. -/
theorem smoothMert_eq_add {N Q : ℕ} (hN : N < (Q + 1) ^ 2) :
    smoothMert Q N = mert N
      + ∑ p ∈ (Finset.Ioc 0 N).filter (fun p => p.Prime ∧ Q < p), mert (N / p) := by
  have h : ∀ p ∈ (Finset.Ioc 0 N).filter (fun p => p.Prime ∧ Q < p),
      smoothMert Q (N / p) = mert (N / p) := by
    intro p hp
    simp only [Finset.mem_filter] at hp
    exact smoothMert_eq_mert (thr_le_of_lt hN hp.2.2)
  rw [mert_eq_smoothMert_sub hN, Finset.sum_congr rfl h]
  ring

/-- The prime-sum rearrangement:  `Σ_{Q<p≤N} M(N/p) = Σ_{k ≤ thr} μ(k)(π(N/k) − π(Q))`.
The terms with `k > thr` vanish because then `N/k < q₁`, which is why the outer
range is `thr` and not `N`.  Minimality of `q₁` is used exactly here. -/
theorem sum_mert_div {N Q q : ℕ} (hq : q.Prime) (hQq : Q < q)
    (hmin : ∀ p, p.Prime → Q < p → q ≤ p) :
    ∑ p ∈ (Finset.Ioc 0 N).filter (fun p => p.Prime ∧ Q < p), mert (N / p)
      = ∑ k ∈ Finset.Ioc 0 (N / q), μ k * ((primeCount (N / k) : ℤ) - primeCount Q) := by
  set T := (Finset.Ioc 0 N).filter (fun p => p.Prime ∧ Q < p) with hT
  -- stretch every inner range to the common `Ioc 0 thr`
  have hstretch : ∀ p ∈ T,
      mert (N / p) = ∑ k ∈ Finset.Ioc 0 (N / q), (if k ≤ N / p then μ k else 0) := by
    intro p hp
    simp only [hT, Finset.mem_filter, Finset.mem_Ioc] at hp
    obtain ⟨⟨hp0, hpN⟩, hpp, hpQ⟩ := hp
    have hdiv : N / p ≤ N / q := by
      rw [Nat.le_div_iff_mul_le hq.pos]
      exact le_trans (Nat.mul_le_mul le_rfl (hmin p hpp hpQ)) (Nat.div_mul_le_self N p)
    rw [mert, ← Finset.sum_filter]
    congr 1
    ext k
    simp only [Finset.mem_filter, Finset.mem_Ioc]
    omega
  rw [Finset.sum_congr rfl hstretch, Finset.sum_comm]
  refine Finset.sum_congr rfl fun k hk => ?_
  rw [Finset.mem_Ioc] at hk
  obtain ⟨hk0, hkthr⟩ := hk
  -- `Q < N/k`, so `π(N/k) ≥ π(Q)` and the window `(Q, N/k]` makes sense
  have hqk : q ≤ N / k := by
    rw [Nat.le_div_iff_mul_le hk0, mul_comm]
    exact (Nat.le_div_iff_mul_le hq.pos).mp hkthr
  have hQk : Q ≤ N / k := le_trans (by omega) hqk
  -- the window, as a filter of `T` and as a filter of `primesLe (N/k)`
  have hwin : T.filter (fun p => k ≤ N / p)
      = (primesLe (N / k)).filter (fun p => Q < p) := by
    ext p
    simp only [hT, primesLe, Finset.mem_filter, Finset.mem_Ioc]
    constructor
    · intro h
      refine ⟨⟨⟨h.1.1.1, ?_⟩, h.1.2.1⟩, h.1.2.2⟩
      rw [Nat.le_div_iff_mul_le hk0, mul_comm]
      exact (Nat.le_div_iff_mul_le h.1.1.1).mp h.2
    · intro h
      refine ⟨⟨⟨h.1.1.1, le_trans h.1.1.2 (Nat.div_le_self N k)⟩, h.1.2, h.2⟩, ?_⟩
      rw [Nat.le_div_iff_mul_le h.1.1.1, mul_comm]
      exact (Nat.le_div_iff_mul_le hk0).mp h.1.1.2
  -- the window's cardinality is `π(N/k) − π(Q)`
  have hcard : ((T.filter (fun p => k ≤ N / p)).card : ℤ)
      = (primeCount (N / k) : ℤ) - primeCount Q := by
    have hsplit := Finset.card_filter_add_card_filter_not (s := primesLe (N / k))
      (fun p => Q < p)
    have hlow : (primesLe (N / k)).filter (fun p => ¬ Q < p) = primesLe Q := by
      ext p
      simp only [primesLe, Finset.mem_filter, Finset.mem_Ioc, Nat.not_lt]
      constructor
      · intro h; exact ⟨⟨h.1.1.1, h.2⟩, h.1.2⟩
      · intro h; exact ⟨⟨⟨h.1.1, le_trans h.1.2 hQk⟩, h.2⟩, h.1.2⟩
    rw [hwin, primeCount, primeCount, ← hlow]
    omega
  rw [← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul, hcard, mul_comm]

/-- **Lemma 22 of the paper --- the exact one-large-prime identity.**

With `q₁` the least prime `> Q` and `thr = ⌊N/q₁⌋`,

    Σ_{thr < n ≤ N, n squarefree, P⁺(n) ≤ Q} μ(n)
      = M(N) − M(thr) + Σ_{k ≤ thr} μ(k)·(π(N/k) − π(Q)) .

Exact for every `N` with `N < (Q+1)^2`; no estimate of `M` or `π` is used.  The
left side is `R·p₁ = Σ_{n ∈ 𝔅} μ(n)`, since `μ` vanishes off
the squarefree integers. -/
theorem one_large_prime_identity {N Q q : ℕ} (hN : N < (Q + 1) ^ 2) (hq : q.Prime)
    (hQq : Q < q) (hmin : ∀ p, p.Prime → Q < p → q ≤ p) :
    ∑ n ∈ (Finset.Ioc (N / q) N).filter (fun n => Squarefree n ∧ Smooth Q n), μ n
      = mert N - mert (N / q)
        + ∑ k ∈ Finset.Ioc 0 (N / q),
            μ k * ((primeCount (N / k) : ℤ) - primeCount Q) := by
  have hthrN : N / q ≤ N := Nat.div_le_self N q
  -- drop the redundant squarefreeness from the left side
  have hdrop : ∑ n ∈ (Finset.Ioc (N / q) N).filter (fun n => Squarefree n ∧ Smooth Q n), μ n
      = ∑ n ∈ (Finset.Ioc (N / q) N).filter (Smooth Q), μ n := by
    rw [Finset.sum_filter, Finset.sum_filter]
    refine Finset.sum_congr rfl fun n _ => ?_
    by_cases hs : Squarefree n
    · simp [hs]
    · simp [hs, moebius_eq_zero_of_not_squarefree hs]
  -- the band sum is `A(N) − A(thr)`
  have hband : smoothMert Q (N / q)
      + ∑ n ∈ (Finset.Ioc (N / q) N).filter (Smooth Q), μ n = smoothMert Q N := by
    rw [smoothMert, smoothMert, Finset.sum_filter, Finset.sum_filter, Finset.sum_filter]
    exact Finset.sum_Ioc_consecutive _ (Nat.zero_le _) hthrN
  rw [hdrop, ← sub_eq_iff_eq_add'.mpr hband.symm, smoothMert_eq_mert (thr_le_of_lt hN hQq),
    smoothMert_eq_add hN, sum_mert_div hq hQq hmin]
  ring

/-- The hypotheses of `one_large_prime_identity` are **satisfiable**, so the
identity is not vacuous.  Witness: `N = 50`, `Q = 7 = ⌊√50⌋` (`49 ≤ 50 < 64`),
`q₁ = 11`, `thr = ⌊50/11⌋ = 4` — a cell with `N` not a square and not on the
`2^e` ladder, which is the case Lemma 7(i) of the paper was introduced to cover.
Both sides evaluate to `2`; that numerical check is compiled-code evidence, and
this theorem is the part of it proved here, using no axioms beyond Lean's
standard ones (`axioms.txt` lists them). -/
theorem one_large_prime_identity_at_50 :
    ∑ n ∈ (Finset.Ioc (50 / 11) 50).filter (fun n => Squarefree n ∧ Smooth 7 n), μ n
      = mert 50 - mert (50 / 11)
        + ∑ k ∈ Finset.Ioc 0 (50 / 11),
            μ k * ((primeCount (50 / k) : ℤ) - primeCount 7) := by
  refine one_large_prime_identity (by norm_num) (by norm_num) (by norm_num) ?_
  intro p hp h
  by_contra hc
  interval_cases p <;> revert hp <;> decide

/-- And `Q = 7` really is `⌊√50⌋`. -/
theorem sqrt_fifty : Nat.sqrt 50 = 7 := by norm_num

/-! ## §5 The parity split

`U_a = B_a + R_a` counts the band class unsigned and
`M_a = ±(B_a − R_a)` counts it with `μ`, the sign set by the parity of `ω(a)`.
Hypothesis 5 then needs only `|M_a|` bounded away from `U_a` by a fixed factor, not
`|M_a| = o(U_a)`. -/

/-- `μ(n) = (−1)^{ω(n)}` for squarefree `n`, with `ω` the count of **distinct**
prime factors (`ArithmeticFunction.cardDistinctFactors`).  That is the form
the paper uses to split the band into `vB` (`ω` even) and `vR` (`ω` odd);
mathlib's `moebius_apply_of_squarefree` states it with `Ω`. -/
theorem moebius_eq_neg_one_pow_omega {n : ℕ} (hn : Squarefree n) :
    μ n = (-1) ^ cardDistinctFactors n := by
  rw [moebius_apply_of_squarefree hn,
    ← (cardDistinctFactors_eq_cardFactors_iff_squarefree hn.ne_zero).mpr hn]

/-- **The sign.**  On the class `a_Y(n) = a` one has `n = a·b`
with `a`, `b` coprime (automatic from `P⁻(b) > Y ≥ P⁺(a)`), and
`μ(b) = (−1)^{ω(a)}·μ(n)`.  So `Σ_b μ(b)` over the class is `±(B_a − R_a)`, the
sign being the parity of `ω(a)` — which is exactly the disjunction `parity_split`
takes as its hypothesis. -/
theorem moebius_class_sign {a b : ℕ} (hco : Nat.Coprime a b)
    (hsq : Squarefree (a * b)) :
    μ b = (-1) ^ cardDistinctFactors a * μ (a * b) := by
  have ha : Squarefree a := hsq.squarefree_of_dvd (Dvd.intro b rfl)
  rw [isMultiplicative_moebius.map_mul_of_coprime hco,
    moebius_eq_neg_one_pow_omega ha, ← mul_assoc, ← pow_add,
    Even.neg_one_pow ⟨cardDistinctFactors a, rfl⟩, one_mul]

/-- **The split, the first half of (24) of the paper.**  `{B_a, R_a} = {(U_a + M_a)/2, (U_a − M_a)/2}`,
which is which being set by the parity of `ω(a)`. -/
theorem parity_split {U M B R : ℝ} (hU : B + R = U) (hM : B - R = M ∨ R - B = M) :
    (B = (U + M) / 2 ∧ R = (U - M) / 2) ∨ (B = (U - M) / 2 ∧ R = (U + M) / 2) := by
  rcases hM with h | h
  · exact Or.inl ⟨by linarith, by linarith⟩
  · exact Or.inr ⟨by linarith, by linarith⟩

/-- **The ratio bound (24) of the paper.**  If `|M_a| < U_a` then `B_a > 0` --- no retained class is
dead --- and `R_a/B_a ≤ (U_a + |M_a|)/(U_a − |M_a|)`. -/
theorem par_ratio {U M B R : ℝ} (hU : B + R = U) (hM : B - R = M ∨ R - B = M)
    (h : |M| < U) : 0 < B ∧ R / B ≤ (U + |M|) / (U - |M|) := by
  have habs : |2 * B - U| ≤ |M| := by
    rcases hM with h' | h'
    · rw [show 2 * B - U = M by linarith]
    · rw [show 2 * B - U = -M by linarith, abs_neg]
  rw [abs_le] at habs
  have hM0 : 0 ≤ |M| := abs_nonneg M
  have hB : 0 < B := by linarith
  refine ⟨hB, ?_⟩
  rw [div_le_div_iff₀ hB (by linarith)]
  nlinarith [habs.1, habs.2]

/-! ## §6 The cone step

`x ≥ 0` puts `Ax` in the cone `K_Q` of (5), so
`dist₂(t, K_Q) ≤ ‖t − Ax‖₂`, and the whole proof is an estimate of `‖t − Ax‖₂`.
That inequality is the triviality that an infimum is at most any value; it is
recorded because it is the step that turns a distance into a norm. -/

/-- The conical hull of a finite family: the `K_Q` of (5) of the paper. -/
def coneHull {ι E : Type*} [Fintype ι] [AddCommMonoid E] [Module ℝ E]
    (Acol : ι → E) : Set E :=
  {y | ∃ x : ι → ℝ, (∀ i, 0 ≤ x i) ∧ y = ∑ i, x i • Acol i}

/-- **The cone step.**  Any non-negative combination of the columns bounds the
cone distance.  Stated for a normed real vector space; `EuclideanSpace ℝ (Fin n)`
with `ι = vB` is the instance the paper uses. -/
theorem infDist_coneHull_le {ι E : Type*} [Fintype ι] [NormedAddCommGroup E]
    [NormedSpace ℝ E] (Acol : ι → E) (t : E) (x : ι → ℝ) (hx : ∀ i, 0 ≤ x i) :
    Metric.infDist t (coneHull Acol) ≤ ‖t - ∑ i, x i • Acol i‖ := by
  rw [← dist_eq_norm]
  exact Metric.infDist_le_dist_of_mem ⟨x, hx, rfl⟩

/-- The Euclidean instance, for the record. -/
theorem infDist_coneHull_le_euclidean {n : ℕ} {ι : Type*} [Fintype ι]
    (Acol : ι → EuclideanSpace ℝ (Fin n)) (t : EuclideanSpace ℝ (Fin n))
    (x : ι → ℝ) (hx : ∀ i, 0 ≤ x i) :
    Metric.infDist t (coneHull Acol) ≤ ‖t - ∑ i, x i • Acol i‖ :=
  infDist_coneHull_le Acol t x hx

/-! ## §7 Bonferroni

Lemma 15 of the paper truncates Legendre's identity at `J = ⌈log L⌉` terms and uses
the exact value of the truncated alternating binomial sum.  The identity is
mathlib's `Int.alternating_sum_range_choose_eq_choose`; what is added here is
the project's form (`k ≥ 1`, index shifted) and the error statement.

The identity was verified numerically for `1 ≤ k ≤ 13`, `0 ≤ J ≤ 13`.
`bonferroni` replaces that check by a proof for all `k` and `J`. -/

/-- **The Bonferroni identity used in Lemma 15 of the paper.**  For `k ≥ 1`,
`Σ_{j ≤ J} (−1)^j C(k,j) = (−1)^J C(k−1,J)`. -/
theorem bonferroni {k : ℕ} (hk : 1 ≤ k) (J : ℕ) :
    ∑ j ∈ Finset.range (J + 1), ((-1) ^ j * (k.choose j) : ℤ)
      = (-1) ^ J * ((k - 1).choose J) := by
  obtain ⟨m, rfl⟩ : ∃ m, k = m + 1 := ⟨k - 1, by omega⟩
  simpa using Int.alternating_sum_range_choose_eq_choose (n := m) (m := J)

/-- The truncation error is **exactly** `C(k−1,J)` for `k ≥ 1`: the full Legendre
sum is `0` there, so the partial sum *is* the error. -/
theorem bonferroni_error_eq {k : ℕ} (hk : 1 ≤ k) (J : ℕ) :
    |∑ j ∈ Finset.range (J + 1), ((-1) ^ j * (k.choose j) : ℤ)|
      = ((k - 1).choose J : ℤ) := by
  rw [bonferroni hk J, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
    Nat.abs_cast]

/-- **The truncated inclusion–exclusion sign and size.**  For every `k` the
partial sum differs from the Legendre indicator `1_{k = 0}` by at most `C(k,J)`
--- the bound (20) of the paper is summed from this. -/
theorem bonferroni_error_le (k J : ℕ) :
    |∑ j ∈ Finset.range (J + 1), ((-1) ^ j * (k.choose j) : ℤ)
      - (if k = 0 then 1 else 0)| ≤ (k.choose J : ℤ) := by
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · rw [if_pos rfl]
    have hone : ∑ j ∈ Finset.range (J + 1), ((-1) ^ j * ((0 : ℕ).choose j) : ℤ) = 1 := by
      rw [Finset.sum_eq_single_of_mem 0 (Finset.mem_range.mpr (Nat.succ_pos J))]
      · simp
      · intro j _ hj0
        rw [Nat.choose_eq_zero_of_lt (Nat.pos_of_ne_zero hj0)]
        simp
    rw [hone, sub_self, abs_zero]
    positivity
  · rw [if_neg (by omega), sub_zero, bonferroni_error_eq hk J]
    exact_mod_cast Nat.choose_le_choose J (by omega)

end KappaZero

#print axioms KappaZero.sqrt_bounds
#print axioms KappaZero.card_filter_lt_primeFactors_le_one
#print axioms KappaZero.not_sq_dvd_of_lt
#print axioms KappaZero.thr_le_of_lt
#print axioms KappaZero.smallPart_dvd
#print axioms KappaZero.dvd_smallPart_iff
#print axioms KappaZero.tRow_eq_sum
#print axioms KappaZero.AxRow_eq_sum
#print axioms KappaZero.identity_E
#print axioms KappaZero.primeCount_eq
#print axioms KappaZero.smoothMert_eq_mert
#print axioms KappaZero.bigPrime_eq
#print axioms KappaZero.bigPrime_spec
#print axioms KappaZero.smooth_div_bigPrime
#print axioms KappaZero.mert_eq_smoothMert_sub
#print axioms KappaZero.smoothMert_eq_add
#print axioms KappaZero.sum_mert_div
#print axioms KappaZero.one_large_prime_identity
#print axioms KappaZero.one_large_prime_identity_at_50
#print axioms KappaZero.sqrt_fifty
#print axioms KappaZero.moebius_eq_neg_one_pow_omega
#print axioms KappaZero.moebius_class_sign
#print axioms KappaZero.parity_split
#print axioms KappaZero.par_ratio
#print axioms KappaZero.infDist_coneHull_le
#print axioms KappaZero.infDist_coneHull_le_euclidean
#print axioms KappaZero.bonferroni
#print axioms KappaZero.bonferroni_error_eq
#print axioms KappaZero.bonferroni_error_le

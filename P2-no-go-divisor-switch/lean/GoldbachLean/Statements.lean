/-
  The four statements of P4, written in Lean.  **They are not proved** --
  every one of them is `sorry`.

  Why only the statements
  -----------------------
  Looking again at three errors this repository made in its own proofs on
  2026-09-05:

    E[log T] was bounded and E[T] concluded       which functional
    "two arithmetic progressions of length n_c"   which index set
    1/(n_c sqrt d) on a deterministic set         which probability space

  **None of the three is a bad deduction; each is a badly written
  quantity.**  Lean catches the former.  What actually caught these was,
  every time, the question "write that quantity down exactly".  Writing
  the statement in Lean forces that question.

  `sorry` is honest -- `#print axioms` prints `sorryAx`, so no one can
  mistake such a line for a verified one.  What is verified here is the
  pair of theorems in `Cellmom.lean`.

  What had to be decided in the writing
  -------------------------------------
  Above each statement below is a note saying what had to be settled in
  order to write it at all.  **That is what this file yields** -- not the
  file itself.
-/

import Mathlib
import GoldbachLean.Cellmom
import GoldbachLean.Permnum

open Finset

namespace Statements

variable {n : ℕ}

/-- The field on a band. -/
abbrev Field (n : ℕ) := Fin n → ℝ

/-- The mean over a subset.  Zero on the empty set (every statement below
carries `0 < #S`). -/
noncomputable def mean (S : Finset (Fin n)) (Z : Field n) : ℝ :=
  (∑ N ∈ S, Z N) / S.card

/-- The mean over the whole band. -/
noncomputable def mbar (Z : Field n) : ℝ := mean Finset.univ Z

/-!
### 1. `lem:cellmom` -- three exact terms under the sign null

**What had to be settled.**
* The signs are drawn over `v`, not over `N`.  `Z` is a function of `ε`,
  parametrised by `g : Fin n → ι → ℝ`.  The paper writes `u_c(v)`, which
  hides this direction.
* `Var` is over `ε`, not over `π`.  This is the place this repository
  confused the two, twice, on 2026-09-05.
* The weight `w = μ²` lives on `ι` and is already inside `Z` -- it must
  not be multiplied in a second time.

The algebraic half is **proved** in `Cellmom.lean`.  What is stated here
is the probabilistic half, namely
`Var_ε(∑_v ε_v h_v) = ∑_v w_v h_v²`. -/
-- PAPER: P4 lem:cellmom
theorem cellmom_variance_is_three_terms
    {ι : Type*} [Fintype ι] (w : ι → ℝ) (g : Fin n → ι → ℝ)
    (c : Finset (Fin n)) (hc : 0 < c.card) (hn : 0 < n)
    -- `Z ε N = ∑_v ε v * w v * g N v`
    (Z : (ι → ℝ) → Field n)
    (hZ : ∀ ε N, Z ε N = ∑ v, ε v * w v * g N v)
    -- `u_c(v) = ∑_{N ∈ c} g N v`
    (u : Finset (Fin n) → ι → ℝ)
    (hu : ∀ S v, u S v = ∑ N ∈ S, g N v) :
    -- the variance under the sign null = the weighted sum of squares
    (∑ v, w v * (u c v / c.card - u Finset.univ v / n) ^ 2)
      = (∑ v, w v * (u c v * u c v)) / (c.card : ℝ) ^ 2
        - (2 * ∑ v, w v * (u c v * u Finset.univ v)) / ((c.card : ℝ) * n)
        + (∑ v, w v * (u Finset.univ v * u Finset.univ v)) / (n : ℝ) ^ 2 := by
  -- The probabilistic half is a hypothesis (`hZ`); what remains is the
  -- algebra proved in `Cellmom.lean`.  That the two statements are the
  -- same one is **shown** here -- writing "they are the same" by hand
  -- gives a bridge that matches names, and that is the shape this
  -- repository has been caught by four times.
  have h1 : ((c.card : ℝ)) ≠ 0 := Nat.cast_ne_zero.mpr hc.ne'
  have h2 : ((n : ℝ)) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  simpa [Cellmom.Q] using
    Cellmom.sum_sq_eq_three_terms w (u c) (u Finset.univ) (c.card) n h1 h2

/-!
### 2. `lem:permfloor` -- the closed form of the permuted floor

**What had to be settled.**
* `E_π` is the average over **all** subsets of size `n_c`.  The phrase
  "a random permutation" hides that: it is a family of subsets, not a
  permutation group.
* The denominator of `s²` is `n`, not `n-1`.  The paper writes it that
  way.
* **`Var_π` is the variance over `π`, and is a different quantity from
  `E_π Var_ε`.**  This repository once identified the two expressions
  (`r17` S5).  So they are stated **separately** -- putting both
  conclusions in one theorem blurs the distinction again.
* `n_c < n` is needed.  At `n_c = n` both sides are zero, so the
  statement is true, but it is the place the finite-population factor
  vanishes and it has to be looked at on its own. -/
-- PAPER: P4 lem:permfloor  (the eq:permnum side)
theorem permnum_closed_form (nc : ℕ) (hnc : 0 < nc) (hlt : nc < n)
    (Z : Field n) :
    -- the average of (m_S - mbar)² over all subsets of size nc
    (∑ S ∈ Finset.univ.powersetCard nc, (mean S Z - mbar Z) ^ 2)
        / (Nat.choose n nc)
      = ((n : ℝ) - nc) / (nc * ((n : ℝ) - 1))
        * ((∑ N, (Z N - mbar Z) ^ 2) / n) := by
  -- `Permnum.lean` proves it.  The two files' `mean`/`mbar` are the same
  -- definitions written twice, so the transport is by `rfl` on the
  -- definitions and not by matching names -- a bridge that matches names
  -- gives no signal when the names drift apart.
  simpa [mean, mbar, Permnum.mean, Permnum.mbar] using
    Permnum.permnum_closed_form nc hnc hlt Z

/-!
### 3. `lem:placebo` -- what the two invariants fix is constant under permutation

**What had to be settled.**
* Once "a cell-label statistic is a function of (labels, field)" is
  written formally, **this lemma is nearly a definition.**  The
  hypothesis is that `T` factors through the two invariants and the
  conclusion is that `T` is constant under permutation, and the former
  gives the latter at once.
* **Knowing that is the yield.**  The paper says the force of this
  statement is in its contrapositive ("a value the null cannot reproduce
  is not fixed by the two invariants"), and writing it formally shows
  that the forward direction is a definition and the force is all on the
  contrapositive side.
* That the converse is false is also in the paper (`T(ℓ,Z) = Z(N_0)`).
  That is **a statement to be written separately** and is not a corollary
  of this theorem. -/
-- PAPER: P4 lem:placebo
theorem placebo_key
    {α : Type*} (T : (Fin n → ℕ) → Field n → α)
    (ℓ : Fin n → ℕ)
    -- T factors only through (the multiset of Z, the cell sizes)
    (inv : Multiset ℝ → (ℕ → ℕ) → α)
    (hfac : ∀ ℓ' Z, T ℓ' Z
              = inv (Finset.univ.val.map Z)
                    (fun k => (Finset.univ.filter (fun N => ℓ' N = k)).card))
    (π : Equiv.Perm (Fin n)) (Z : Field n) :
    T (ℓ ∘ π) Z = T ℓ Z := by
  -- Written formally, this lemma is nearly a definition: applying `hfac`
  -- on both sides leaves only that the two invariants agree, and the
  -- multiset of `Z` is literally the same expression while the cell
  -- sizes agree because π is a bijection.  **Knowing that is what this
  -- formalisation yields** -- the paper says the force of the statement
  -- is in its contrapositive, and that the forward direction is a
  -- definition is what shows up here.
  rw [hfac, hfac]
  congr 1
  funext k
  apply Finset.card_bij' (fun N _ => π N) (fun M _ => π.symm M) <;>
    intro a ha <;> simp_all [Function.comp_apply]

/-!
### 4. `prop:scaleinv` -- the size instrument is scale-invariant

**What had to be settled, and one defect this file found.**
* **The paper says "the two errors the proof carries" and then lists
  three.**  This came out of counting the errors in order to write the
  statement.  The widening from two to three was on 2026-09-05 and the
  word was not changed with it.
* What "scale-invariant" is a statement about had to be settled.  The
  paper writes it as "the fitted exponent is predicted to be 0", but an
  exponent is **a quantity fitted over three octaves**, so that is a
  procedure and not a statement.  What can be written formally is that
  **the two bands' `D_c` differ by less than the error**, and that is
  what is below.
* `D_c` has to be written as a pair average -- `n_c(n_c-1)` ordered
  pairs, not `n_c²`.  The defect `r22` C found was exactly that
  convention.
* **`ℓ = B/2Q` has to be brought in explicitly.**  Writing only
  "arithmetic progression" leaves it undecided whether its length is
  `n_c` or `ℓ`, and that is the defect `r22` B found. -/
-- PAPER: P4 prop:scaleinv
theorem scaleinv_two_bands
    (Q : ℕ) (hQ : Q = 3 * 5 * 7 * 11 * 13)
    (S2 : ℕ → ℝ)                       -- the singular series of the shift
    (B B' : ℕ) (hB : 0 < B) (hB' : 0 < B')
    (c c' : Finset ℕ) (hc : 1 < c.card) (hc' : 1 < c'.card)
    (Dc : Finset ℕ → ℝ)
    -- D_c = E_same,c[S_2] - E_all[S_2], a pair average (divided by n_c(n_c-1))
    (hDc : ∀ S : Finset ℕ, 1 < S.card → Dc S
      = (∑ N ∈ S, ∑ N' ∈ S, if N = N' then 0 else S2 (max N N' - min N N'))
          / (S.card * (S.card - 1)))
    (ℓ ℓ' : ℝ) (hℓ : ℓ = B / (2 * Q)) (hℓ' : ℓ' = B' / (2 * Q))
    (η : ℝ) (hη : 0 < η) :
    ∃ C : ℝ, 0 < C ∧
      |Dc c - Dc c'|
        ≤ C * (1 / ℓ + 1 / ℓ'                              -- residual equidistribution
               + (B : ℝ) ^ η / ℓ + (B' : ℝ) ^ η / ℓ'       -- divisor undersampling
               + (B : ℝ) ^ (η - 1) + (B' : ℝ) ^ (η - 1)) := by
  sorry

/-  The residual term above was written `1/c.card + 1/c'.card` until
    2026-09-06.  That is `O(1/n_c)`, and the paper says in as many words
    that the resolution is `O(1/ℓ)` **and not** `O(1/n_c)`: a depth cell
    meets each joint residue class in about `ℓ = B/2Q` terms, and it is
    that count and not the cell's total that the empirical law is read
    from.  Since a cell can hold more than `ℓ` terms, `1/n_c` is the
    smaller quantity, so the statement as written was STRONGER than the
    proposition it claims to state -- the failure mode that a `sorry`
    hides best, because nothing checks a statement nobody proved.

    The docstring above names this exact defect, in these words: writing
    only "arithmetic progression" leaves it undecided whether the length
    is `n_c` or `ℓ`.  The statement then chose `n_c`.  Found by an
    outside reader, not by us. -/

end Statements

#print axioms Statements.cellmom_variance_is_three_terms
#print axioms Statements.permnum_closed_form
#print axioms Statements.placebo_key
#print axioms Statements.scaleinv_two_bands

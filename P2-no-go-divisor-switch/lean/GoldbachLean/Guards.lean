/-
  Guards -- not theorems, but statements that make **a wrong reading
  contradict something machine-checked**.

  Why these are worth more than the theorems
  ------------------------------------------
  Most of the defects this repository caught on 2026-09-05 were not
  broken proofs but **a true statement restated slightly too strongly**.
  Formalising the theorem itself does not catch that kind of error --
  the theorem is true.  To catch it one has to write down, separately,
  **a statement that refutes the wrong reading**.

  The two below are of that shape.  Neither carries a `sorry`.
-/

import Mathlib

open Finset

namespace Guards

/-!
## Guard 1 -- `eq:permfloor` is a π-average, not the value at one π

This repository crossed that line **twice**.

* `r17` S5 -- it identified `s²/n_c` with `eq:permfloor`.  It is in fact
  `eq:permnum`, and `eq:permfloor` is that quantity's expectation under
  the sign null.
* Afterwards the corrected text of `lem:placebo` had to say that the
  floor at one π is not determined by the two invariants.

The mathematical fact is simple: `E_π[f(π)] = c` gives `f(π) = c` only
when `f` is constant.  So the guard is **to exhibit that `f` is not
constant**, and one counterexample suffices.

`Var(m_S - mbar) = n_S^{-2} · ∑_{N,N'} J_S(N) J_S(N') K(N,N')` with
`J_S = 1_S - n_S/n`.  Below is a kernel at `n = 3` for which two cells
of size 1 take **different values**.  With a diagonal-only kernel
(`K = δ`) the value depends on the size alone and is therefore
constant, so turning on one off-diagonal entry is the point --
**that term is what makes the cell visible.** -/

/-- The doubly centred vector `J_S = 1_S - |S|/n`, written out
componentwise.  Leaving the membership condition as an `if` keeps the
rational arithmetic from reducing.  Here `n = 3` and the two cells of
size 1 are `S = {0}` and `S = {2}`. -/
def J0 : Fin 3 → ℚ := ![2/3, -1/3, -1/3]

def J2 : Fin 3 → ℚ := ![-1/3, -1/3, 2/3]

/-- A kernel with one off-diagonal entry turned on: the diagonal is 1,
`K(0,1) = K(1,0) = t`, and the rest is 0. -/
def K3 (t : ℚ) : Fin 3 -> Fin 3 → ℚ :=
  ![![1, t, 0], ![t, 1, 0], ![0, 0, 1]]

/-- A cell's floor -- the doubly centred quadratic form (no division,
since the size is 1). -/
def floorV (t : ℚ) (J : Fin 3 → ℚ) : Rat :=
  Finset.univ.sum fun i => Finset.univ.sum fun j => J i * J j * K3 t i j

/-- **Two cells of the same size give different floors.**

So the floor at one `pi` is not determined by the cell's **size** alone.
Reading `eq:permfloor` as "the floor of a random cell" contradicts this
statement -- what that equation gives is the **average** of those
values. -/
-- PAPER: P4 lem:permfloor  (guard: a pi-average, not the value at one pi)
theorem per_cell_floor_is_not_a_function_of_size :
    floorV 1 J0 ≠ floorV 1 J2 := by
  simp [floorV, K3, J0, J2, Fin.sum_univ_three]
  norm_num

/-- And by how much they differ. -/
-- PAPER: P4 lem:permfloor  (guard: the size of that gap)
theorem per_cell_floor_gap :
    floorV 1 J2 - floorV 1 J0 = 2 / 3 := by
  simp [floorV, K3, J0, J2, Fin.sum_univ_three]
  norm_num

/-- And what produces that difference is the **off-diagonal term** --
with a diagonal-only kernel, cells of equal size have equal floors.  The
same fact the paper states as "all three terms are needed" appears here
along the cell axis. -/
-- PAPER: P4 note:threeterms  (guard: the off-diagonal is what shows the cell)
theorem diagonal_only_kernel_sees_only_size :
    floorV 0 J0 = floorV 0 J2 := by
  simp [floorV, K3, J0, J2, Fin.sum_univ_three]
  norm_num

/-!
## Guard 2 -- one feasible point bounds the minimum from above

A theorem of this shape settled a question in the working layer:
`f = (1/R[0])·1` is feasible because `f ≥ 0`, hence
`κ ≤ ‖Mf − b‖/‖b‖`.

**The point is that this inequality holds whatever `M` and `b` are.**
That is why it clears the "model barrier" which killed the earlier
reconstructions of the same question, and on the day it was found the
right-hand side changed twice while the inequality never moved.

What is written below is that fact alone -- the infimum of a set that is
bounded below is no larger than any element of the set. -/
-- PAPER: none  (a working-layer theorem, not deployed)
theorem feasible_point_bounds_the_min
    {ι : Type*} (g : ι → ℝ) (S : Set ι) (f₀ : ι) (hf₀ : f₀ ∈ S)
    (hbdd : BddBelow (g '' S)) :
    sInf (g '' S) ≤ g f₀ :=
  csInf_le hbdd (Set.mem_image_of_mem g hf₀)

end Guards

#print axioms Guards.per_cell_floor_is_not_a_function_of_size
#print axioms Guards.per_cell_floor_gap
#print axioms Guards.diagonal_only_kernel_sees_only_size
#print axioms Guards.feasible_point_bounds_the_min

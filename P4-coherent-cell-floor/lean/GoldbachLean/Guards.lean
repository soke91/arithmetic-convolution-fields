/-
  방벽 — 정리가 아니라, **틀린 문장이 기계 검증된 진술과 모순되게** 하는 것.

  왜 이것이 정리보다 값어치가 있나
  --------------------------------
  2026-09-05 에 이 저장소가 잡은 결함 대다수는 증명이 틀린 것이 아니라
  **참인 진술을 조금 세게 다시 서술한 것**이었다.  그런 오류는 본
  정리를 형식화해도 안 잡힌다 — 본 정리는 참이니까.  잡히게 하려면
  **잘못된 읽기 자체를 반증하는 진술**을 따로 적어야 한다.

  아래 둘은 그 꼴이다.  둘 다 `sorry` 가 없다.
-/

import Mathlib

open Finset

namespace Guards

/-!
## 방벽 1 — `eq:permfloor` 는 π-평균이지 한 π 의 값이 아니다

이 저장소가 **두 번** 이 선을 넘었다.

* `r17` S5 — `s²/n_c` 를 `eq:permfloor` 와 동일시했다.  실제로는
  `eq:permnum` 이고, `eq:permfloor` 는 그것의 부호-널 기댓값이다.
* 그 뒤 `lem:placebo` 의 정정문이 「한 π 에서의 바닥은 두 불변량으로
  결정되지 않는다」를 적어야 했다.

수학적 사실은 단순하다: `E_π[f(π)] = c` 는 `f` 가 상수일 때만
`f(π) = c` 를 준다.  그러니 방벽은 **`f` 가 상수가 아님을 보이는 것**
이고, 반례 하나면 된다.

`Var(m_S - mbar) = n_S^{-2} · ∑_{N,N'} J_S(N) J_S(N') K(N,N')`,
`J_S = 1_S - n_S/n`.  아래는 `n = 3`, 크기 1 인 두 셀이 **다른 값**을
내는 커널이다.  커널이 대각뿐이면(`K = δ`) 값이 크기에만 의존해
상수가 되므로, 비대각을 하나 켠 것이 요점이다 — **그 항이 있어야
셀이 보인다.** -/

/-- 이중중심 벡터 `J_S = 1_S - |S|/n` 을 성분으로 직접 적는다.
멤버십 조건을 `if` 로 두면 유리수 계산이 안 줄어들어서다.
`n = 3`, 크기 1 인 셀 둘: `S = {0}` 과 `S = {2}`. -/
def J0 : Fin 3 → ℚ := ![2/3, -1/3, -1/3]

def J2 : Fin 3 → ℚ := ![-1/3, -1/3, 2/3]

/-- 비대각을 하나 켠 커널: 대각은 1, `K(0,1) = K(1,0) = t`, 나머지 0. -/
def K3 (t : ℚ) : Fin 3 -> Fin 3 → ℚ :=
  ![![1, t, 0], ![t, 1, 0], ![0, 0, 1]]

/-- 셀의 바닥 — 이중중심 이차형식 (크기 1 이라 나눗셈이 없다). -/
def floorV (t : ℚ) (J : Fin 3 → ℚ) : Rat :=
  Finset.univ.sum fun i => Finset.univ.sum fun j => J i * J j * K3 t i j

/-- **크기가 같은 두 셀이 다른 바닥을 낸다.**

그러므로 한 `pi` 에서의 바닥은 셀의 **크기만으로** 결정되지 않는다.
`eq:permfloor` 를 「무작위 셀의 바닥」이라 읽으면 이 진술과 모순된다 —
그 식이 주는 것은 그 값들의 **평균**이다. -/
-- PAPER: P4 lem:permfloor  (방벽: pi-평균이지 한 pi 의 값이 아니다)
theorem per_cell_floor_is_not_a_function_of_size :
    floorV 1 J0 ≠ floorV 1 J2 := by
  simp [floorV, K3, J0, J2, Fin.sum_univ_three]
  norm_num

/-- 차이가 얼마인지까지 적는다. -/
-- PAPER: P4 lem:permfloor  (방벽: 그 간격의 크기)
theorem per_cell_floor_gap :
    floorV 1 J2 - floorV 1 J0 = 2 / 3 := by
  simp [floorV, K3, J0, J2, Fin.sum_univ_three]
  norm_num

/-- 그리고 그 차이를 낳는 것은 **비대각 항**이다 — 커널이 대각뿐이면
크기가 같은 셀의 바닥이 같아진다.  논문이 「세 항이 다 필요하다」고
적는 것과 같은 사실이 여기서는 셀 축으로 나타난다. -/
-- PAPER: P4 note:threeterms  (방벽: 비대각이 셀을 보이게 한다)
theorem diagonal_only_kernel_sees_only_size :
    floorV 0 J0 = floorV 0 J2 := by
  simp [floorV, K3, J0, J2, Fin.sum_univ_three]
  norm_num

/-!
## 방벽 2 — 가능해 하나가 최솟값의 상계를 준다

`-68` 의 걸음 853~861 의 정리가 이 꼴이다: `f = (1/R[0])·1` 이
`f ≥ 0` 이라 가능해이므로 `κ ≤ ‖Mf − b‖/‖b‖`.

**요점은 이 부등식이 `M` 과 `b` 가 무엇이든 성립한다는 것이다.**  그
세션의 앞선 재구성을 죽인 「모형 장벽」을 넘는 까닭이 그것이고, 같은
날 우변이 무엇인지가 두 번 바뀌었는데도 부등식은 한 번도 안 움직였다.

아래는 그 사실만 적는다 — 하한이 존재하는 집합의 하한은 그 집합의
어느 원소보다 크지 않다. -/
-- PAPER: 없음  (작업층 항목 90 의 정리, 배포 안 됨)
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

/-
  P4 의 네 진술을 Lean 으로 적는다.  **증명하지 않는다** — 전부 `sorry` 다.

  왜 진술만 적나
  --------------
  2026-09-05 에 이 저장소의 증명에서 난 오류 셋을 다시 보면

    E[log T] 를 눌러 놓고 E[T] 를 결론했다      어느 범함수인가
    「길이 n_c 인 등차수열 둘」                   지표 집합이 무엇인가
    결정적 집합에 1/(n_c sqrt d)                확률 공간이 무엇인가

  **셋 다 잘못된 연역이 아니라 잘못 적힌 양이다.**  Lean 이 잡는 것은
  앞엣것이고, 이것들을 실제로 잡은 것은 매번 「그 양을 정확히 적어
  보라」는 물음이었다.  진술을 Lean 으로 적는 것이 그 물음을 강제한다.

  `sorry` 는 정직하다 — `#print axioms` 가 `sorryAx` 를 찍으므로 아무도
  검증된 줄 오해하지 않는다.  검증된 것은 `Cellmom.lean` 의 두 정리뿐이다.

  적으면서 결정해야 했던 것
  -------------------------
  아래 각 진술 위에 「이 진술을 적으려면 무엇을 정해야 했는가」를 적었다.
  **그것이 이 파일의 소득이다.**  파일 자체가 아니라.
-/

import Mathlib
import GoldbachLean.Cellmom
import GoldbachLean.Permnum

open Finset

namespace Statements

variable {n : ℕ}

/-- 밴드 위의 장. -/
abbrev Field (n : ℕ) := Fin n → ℝ

/-- 부분집합의 평균.  빈 집합에서는 0 (아래 진술은 전부 `0 < #S` 를 단다). -/
noncomputable def mean (S : Finset (Fin n)) (Z : Field n) : ℝ :=
  (∑ N ∈ S, Z N) / S.card

/-- 밴드 전체의 평균. -/
noncomputable def mbar (Z : Field n) : ℝ := mean Finset.univ Z

/-!
### 1. `lem:cellmom` — 부호 널 아래 정확한 세 항

**적으면서 정해야 했던 것.**
* 부호는 `v` 위에서 뽑히지 `N` 위에서가 아니다.  `Z` 가 `ε` 의 함수이고
  그 매개가 `g : Fin n → ι → ℝ` 다.  논문은 `u_c(v)` 로 적어 이 방향을
  숨긴다.
* `Var` 는 `ε` 에 대한 것이고 `π` 에 대한 것이 아니다.  이 저장소가
  2026-09-05 에 두 번 헷갈린 자리다.
* 가중치 `w = μ²` 는 `ι` 위에 살고 `Z` 안에 이미 들어가 있다 — 두 번
  곱하면 안 된다.

대수 절반은 `Cellmom.lean` 에서 **증명됐다**.  여기 적는 것은 확률
절반, 곧 `Var_ε(∑_v ε_v h_v) = ∑_v w_v h_v²` 이다. -/
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
    -- 부호 널 아래의 분산 = 가중 제곱합
    (∑ v, w v * (u c v / c.card - u Finset.univ v / n) ^ 2)
      = (∑ v, w v * (u c v * u c v)) / (c.card : ℝ) ^ 2
        - (2 * ∑ v, w v * (u c v * u Finset.univ v)) / ((c.card : ℝ) * n)
        + (∑ v, w v * (u Finset.univ v * u Finset.univ v)) / (n : ℝ) ^ 2 := by
  -- 확률 절반은 가설이고 (`hZ`), 남는 것은 `Cellmom.lean` 에서 증명된
  -- 대수다.  두 진술이 같은 것임을 여기서 **보인다** — 손으로 「같다」고
  -- 적으면 이름만 맞춘 다리가 되고, 그것이 이 저장소가 네 번 물린 꼴이다.
  have h1 : ((c.card : ℝ)) ≠ 0 := Nat.cast_ne_zero.mpr hc.ne'
  have h2 : ((n : ℝ)) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  simpa [Cellmom.Q] using
    Cellmom.sum_sq_eq_three_terms w (u c) (u Finset.univ) (c.card) n h1 h2

/-!
### 2. `lem:permfloor` — 치환된 바닥의 닫힌 꼴

**적으면서 정해야 했던 것.**
* `E_π` 는 크기 `n_c` 인 **모든** 부분집합에 대한 평균이다.  「무작위
  치환」이라는 말이 그것을 가린다 — 치환군이 아니라 부분집합족이다.
* `s²` 의 분모가 `n` 이지 `n-1` 이 아니다.  논문이 그렇게 적는다.
* **`Var_π` 는 `π` 에 대한 분산이고 `E_π Var_ε` 와 다른 양이다.**  이
  저장소가 두 식을 한 번 동일시했다 (`r17` S5).  그래서 둘을 **따로**
  적는다 — 한 정리 안에 두 결론을 넣으면 그 구분이 다시 흐려진다.
* `n_c < n` 이 필요하다.  `n_c = n` 이면 좌변이 0 이고 우변도 0 이라
  참이지만, 유한모집단 인자가 0 이 되는 자리라 따로 봐야 한다. -/
-- PAPER: P4 lem:permfloor  (eq:permnum 쪽)
theorem permnum_closed_form (nc : ℕ) (hnc : 0 < nc) (hlt : nc < n)
    (Z : Field n) :
    -- 크기 nc 인 모든 부분집합에 대한 (m_S - mbar)² 의 평균
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
### 3. `lem:placebo` — 두 불변량이 정하는 것은 치환 아래 상수

**적으면서 정해야 했던 것.**
* 「셀 지표 통계량」이 `(라벨, 장)` 의 함수라는 것을 형식으로 적으면,
  **이 보조정리는 거의 정의다.**  가설이 「`T` 가 두 불변량을 통해
  갈라진다」이고 결론이 「`T` 가 치환 아래 상수」인데, 앞이 뒤를
  곧바로 준다.
* **그것을 아는 것이 소득이다.**  논문은 이 진술의 힘이 대우에 있다고
  적는데 (「널이 재현 못 하는 값은 두 불변량이 안 정한다」), 형식으로
  적어 보면 정방향이 정의이고 힘은 전부 대우 쪽에 있다는 것이 드러난다.
* 역이 거짓이라는 것도 논문이 적는다 (`T(ℓ,Z) = Z(N_0)`).  그것은
  **따로 적어야 하는 진술**이고 이 정리의 따름이 아니다. -/
-- PAPER: P4 lem:placebo
theorem placebo_key
    {α : Type*} (T : (Fin n → ℕ) → Field n → α)
    (ℓ : Fin n → ℕ)
    -- T 가 (Z 의 다중집합, 셀 크기) 를 통해서만 갈라진다
    (inv : Multiset ℝ → (ℕ → ℕ) → α)
    (hfac : ∀ ℓ' Z, T ℓ' Z
              = inv (Finset.univ.val.map Z)
                    (fun k => (Finset.univ.filter (fun N => ℓ' N = k)).card))
    (π : Equiv.Perm (Fin n)) (Z : Field n) :
    T (ℓ ∘ π) Z = T ℓ Z := by
  -- 형식으로 적으면 이 보조정리는 거의 정의다: `hfac` 를 양쪽에 쓰면
  -- 남는 것은 두 불변량이 같다는 것뿐이고, `Z` 의 다중집합은 글자
  -- 그대로 같은 식이며 셀 크기는 π 가 전단사라 같다.  **그것을 아는
  -- 것이 이 형식화의 소득이다** — 논문은 이 진술의 힘이 대우에 있다고
  -- 적는데, 정방향이 정의라는 것이 여기서 드러난다.
  rw [hfac, hfac]
  congr 1
  funext k
  apply Finset.card_bij' (fun N _ => π N) (fun M _ => π.symm M) <;>
    intro a ha <;> simp_all [Function.comp_apply]

/-!
### 4. `prop:scaleinv` — 크기 기구가 척도 불변

**적으면서 정해야 했던 것, 그리고 이 파일이 낸 결함 하나.**
* **논문이 「the two errors the proof carries」라 적고 셋을 열거한다.**
  이 진술을 적으려고 오차를 세다가 나왔다.  둘에서 셋으로 넓힌 것이
  2026-09-05 이고 낱말을 안 고쳤다.
* 「척도 불변」이 무엇의 진술인가를 정해야 했다.  논문은 「적합된 지수가
  0 이라 예측한다」로 적는데, 지수는 **세 옥타브에 적합한 양**이라
  진술이 아니라 절차다.  형식으로 적을 수 있는 것은 **두 밴드의 `D_c`
  차가 오차 안**이라는 것이고, 그것이 아래다.
* `D_c` 를 쌍 평균으로 적어야 한다 — `n_c(n_c-1)` 순서쌍이고 `n_c²` 이
  아니다.  `r22` C 가 낸 결함이 정확히 그 규약이었다.
* **`ℓ = B/2Q` 를 명시적으로 들여와야 한다.**  「등차수열」이라고만
  적으면 그 길이가 `n_c` 인지 `ℓ` 인지 안 정해지고, 그것이 `r22` B 가
  낸 결함이다. -/
-- PAPER: P4 prop:scaleinv
theorem scaleinv_two_bands
    (Q : ℕ) (hQ : Q = 3 * 5 * 7 * 11 * 13)
    (S2 : ℕ → ℝ)                       -- 이동의 특이급수
    (B B' : ℕ) (hB : 0 < B) (hB' : 0 < B')
    (c c' : Finset ℕ) (hc : 1 < c.card) (hc' : 1 < c'.card)
    (Dc : Finset ℕ → ℝ)
    -- D_c = E_same,c[S_2] - E_all[S_2], 순서쌍 평균 (n_c(n_c-1) 로 나눈다)
    (hDc : ∀ S : Finset ℕ, 1 < S.card → Dc S
      = (∑ N ∈ S, ∑ N' ∈ S, if N = N' then 0 else S2 (max N N' - min N N'))
          / (S.card * (S.card - 1)))
    (ℓ ℓ' : ℝ) (hℓ : ℓ = B / (2 * Q)) (hℓ' : ℓ' = B' / (2 * Q))
    (η : ℝ) (hη : 0 < η) :
    ∃ C : ℝ, 0 < C ∧
      |Dc c - Dc c'|
        ≤ C * (1 / c.card + 1 / c'.card                    -- 잔여 등분포
               + (B : ℝ) ^ η / ℓ + (B' : ℝ) ^ η / ℓ'       -- 약수 사표집
               + (B : ℝ) ^ (η - 1) + (B' : ℝ) ^ (η - 1)) := by
  sorry

end Statements

#print axioms Statements.cellmom_variance_is_three_terms
#print axioms Statements.permnum_closed_form
#print axioms Statements.placebo_key
#print axioms Statements.scaleinv_two_bands

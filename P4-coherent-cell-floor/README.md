# Cell-mean fluctuations in an arithmetic convolution field

Reproduction packet for the paper

> **Cell-mean fluctuations in an arithmetic convolution field: exact
> Rademacher-null variances and a permutation control**

This directory holds the source, the compiled PDF, the code, the output
that code produced, and the Lean development.

## What the paper is about

Write `C(N) = Σ_{n<N} Λ(n) μ(N−n)`, `V(N) = Σ_{v<N} μ²(v) Λ(N−v)²` and
`Z(N) = C(N)/√V(N)`. Partition a band of even `N` into cells indexed by
how many of `3, 5, 7, 11, 13` divide `N`, and ask whether a cell mean
departs from the band mean. The question needs a reference scale for
that contrast, and the paper computes the scale exactly.

Under the null that replaces each `μ(v)` on its own support by an
independent Rademacher sign, the variance of any cell contrast is a
doubly centred quadratic form in the field's covariance kernel and
reduces to three convolutions. A uniformly random cell of the same size
recovers the finite-population formula exactly; the divisibility cells
do not, and the size of the gap is predicted with no fitted parameter by
the excess of the shift's singular series over same-cell pairs.

No distributional law for `Z` is asserted, and no Goldbach consequence is
drawn.

## Running it

```
python -m pip install -r ../requirements.txt
python code/<script>.py
```

Python 3.12.10 and NumPy 2.2.5 produced the figures in the paper. Each
script writes the file of the same name in `results/`, and imports only
the standard library and NumPy — never another script here, never a
package in this repository.

**Fifteen scripts, and five of them read another's output.** The paper
cites ten; the other five are here because the packet is closed under
what it consumes, so no script reads a file that is not present. Run
them in this order and every input exists before the script that needs
it:

```
audit_cellfloor_sdz   audit_cn_law   lab_cell_floor
lab_cellmom_montecarlo   lab_mask_placebo   lab_secondcell_predict
audit_cn_kurt_drift   lab_cell_singular   lab_secondcell
audit_cn_class   lab_secondcell2_predict   audit_cn_multnull_deep
lab_secondcell2   audit_cn_coin_deep   lab_maincoef
```

**Some scripts exit nonzero, and that is their registered behaviour.**
A run whose own preregistered rule was refuted says so in its exit code
as well as in its output. Across the ten runs the paper cites,
thirty-eight rules carry a verdict — twenty-four hold, thirteen refuted, one
void — and all of them are reported in the paper's Supplement S1. The
five producer scripts carry rules of their own, which the paper does
not report because nothing in it rests on them. Do not read a nonzero
exit as a broken script.

The results committed here are the ones the paper prints. Rerunning a
script overwrites its result file, so keep a copy if the committed
output is wanted for reference.

Two scripts carry predictions sealed before the measurement that scores
them existed — `lab_secondcell_predict.py` and
`lab_secondcell2_predict.py`; the second reads the first. `lab_maincoef.py`
and the two `predict` files compute in exact rational arithmetic, so
their output reproduces bit for bit anywhere; the rest use double
precision, where a different BLAS can move the last printed digit.

## What each script produces

| script | what it computes |
|---|---|
| `lab_cellmom_montecarlo.py` | Monte-Carlo verification of the closed form of Lemma 2 |
| `lab_cell_floor.py` | the exact null standard deviation of each cell contrast, across eight octaves |
| `lab_mask_placebo.py` | the label-permutation control, ten permutations preserving every cell size |
| `lab_cell_singular.py` | the singular-series excess `D_c` and its exponent across three octaves |
| `lab_maincoef.py` | the main term of `D_c` in exact rational arithmetic |
| `lab_secondcell_predict.py` | predictions for the second cell partition, sealed before it was measured |
| `lab_secondcell2_predict.py` | predictions for the exact statistic on that partition, likewise sealed |
| `lab_secondcell2.py` | the exact measurement on the second partition, by autocorrelation |
| `audit_cellfloor_sdz.py` | the standardised cell contrasts and the count-based comparison |
| `audit_cn_coin_deep.py` | the independent-coin null at the same depth, as a resolution check |

The five further scripts — `audit_cn_law`, `audit_cn_kurt_drift`,
`audit_cn_class`, `audit_cn_multnull_deep` and `lab_secondcell` —
produce inputs the ten above read. The paper does not cite them and
nothing in it rests on them; they are here so that the ten can run.

## Provenance

The paper is a projection of a working repository that is not public.
Everything the paper cites is here, together with everything those
scripts read: the packet is closed under what it consumes, and that
closure was checked by running all fifteen scripts in the order above
and confirming that none dies for a missing file.

The distinction between exact and floating-point arithmetic matters
wherever the paper claims a value is exact rather than measured, and
the two are never mixed in one printed figure.

## The Lean development

`lean/` holds the formalisation, and its `axioms.txt` is the axiom report the paper refers to, produced by
`#print axioms` over every labelled result.

Three of the paper's lemmas are proved — the algebraic identity of
Lemma 2, the closed form of Lemma 12, and Lemma 8 — and each depends on
`propext`, `Classical.choice` and `Quot.sound` and on nothing else.

Proposition 17 is **not** proved. Its statement is written out in
`Statements.lean` and left unproved, and the axiom report says so: the
line for `scaleinv_two_bands` carries `sorryAx`, which is what makes the
absence checkable rather than a matter of trust. What *is* proved is the
finite half of its argument, in `Scaleinv.lean` — Rankin's truncation
step in both directions, the exact count of one residue class in an
interval, and the cancellation of the class count `m`. The analytic
half, the Euler product and its convergence, is not formalised.

Building it needs Lean 4 and Mathlib at the versions pinned in
`lean-toolchain` and `lakefile.toml`; `lake build` reproduces the axiom
report.

## Status

Not submitted and not published.

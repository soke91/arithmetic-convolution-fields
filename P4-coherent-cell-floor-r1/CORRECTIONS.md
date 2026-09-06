# What changed, and why this is a separate directory

`P4-coherent-cell-floor/` holds the bytes that were sent to a journal.
It is frozen and will stay frozen: a paper under review has to be
citable as the thing that was reviewed. This directory holds the same
paper with corrections made after submission. The journal has not been
told; when it is, that will be recorded here.

The corrections came from an independent reader outside this project's
model family, given the submitted file and one question — whether it
was ready to submit. Each was then checked here against the source or
by calculation before being made. Nothing below was accepted on the
reader's word.

## The citation was wrong

The bibliography said that equation (1.5) of Pintz, *On the singular
series in the prime k-tuple conjecture*, gives at `k = 1` the average
`1` of the singular series over shifts. It does not. (1.5) is
Gallagher's global average over `k`-element sets, and at `k = 1` the
sets are singletons, for which the singular series is identically `1`;
the case is an identity carrying no information about shifts.

What the paper uses is Theorem 1 of that note,
`S_H(H) = 1 + O(eps)`, applied with `H` a single element: there
`S(H) = 1` and `S(H ∪ {h})` is the singular series of the shift, so the
average over shifts tends to `1`. Pintz draws the distinction himself,
writing that `S_H(H) → 1` "yields more information about the singular
series than the global average" (1.5).

The source had not been obtained when the bibliography entry was
written. It has been now.

## A sentence in the appendix was false

> `(log p)^2 = (1+o(1))(log N)^2` for `p > N^{1-eps}`

For fixed `eps` the ratio lies in `((1-eps)^2, 1]` and does not tend to
`1`; that gives `(1+O(eps))`. The cut has to move. Taking
`eps = 1/log log N` makes `(1-eps)^2 = 1+o(1)` while
`pi(N^{1-eps}) = o(N/log N)` still holds, and the conclusion stands.
The text now says so.

## The Monte-Carlo precision was justified wrongly, and the direction matters

The paper said the precision `1/sqrt(2(n-1))` is "set by the draw count
alone and is therefore the same at all six cells". That is the Gaussian
case. Each contrast here is a weighted Rademacher sum
`sum_v a_v eps_v`, for which `E S^4 = 3(sum a^2)^2 - 2 sum a^4`, so the
kurtosis is `3 - 2 sum a^4 / (sum a^2)^2`, below Gaussian at every
cell, and the precision is smaller than the bar used — by a factor
that differs from cell to cell.

The direction is the point. A cell-specific bar is *tighter*, so the
sentence "only the cell of three elements clears the bar" is a
statement about the common bar and not about the cells. The text now
says that, and says that the cell-wise fourth moments were not
computed.

## A p-value pointed the wrong way

Ten permutations with the true value ranking first give a Monte-Carlo
p-value of `1/11`, which is the smallest that many permutations can
produce. The paper wrote `p <= 1/11`, which claims a bound on the exact
permutation p-value; ten draws cannot bound that from above at all.

## The abstract claimed a mechanism where a descriptor was transferred

The second partition predicts and then measures `D_c`, the singular
series excess. It does not measure the sign-null floor. So what
transferred out of sample is *residue structure determines `D_c`*. That
`D_c` in turn determines the floor rests on a correlation of `0.9805`
across six cells of the first partition, computed after the fact.

The abstract now separates the two and does not use "mechanism" for the
second link. The decay statement is likewise given as consistency over
eight measured octaves rather than as an identified rate.

## What was reported and is not a defect

The reader also reported that the reproduction directory and the Lean
file promised by the paper are absent. They are not: the paper points
at this repository by URL, and the files resolve anonymously. The
reader had a working checkout in front of it, which is a different tree
and does not contain them.

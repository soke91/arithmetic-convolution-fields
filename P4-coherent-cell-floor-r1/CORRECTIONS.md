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

---

## This directory was revised again on 2026-09-06

The version first published here was corrected the same day, so what is
above is not the whole record. Nobody had been given a link to the
first version and the journal had not been told, so it was replaced
rather than left standing beside a second revision; this section is
what the replacement changed. The frozen submitted directory was not
touched.

### One of the corrections above was itself wrong

Fixing the Monte Carlo precision claim, this note wrote the precision
as `1/sqrt(2(n-1))` times `sqrt(1 - sum a^4/(sum a^2)^2)`, stated as a
value. That is a delta-method approximation. The fourth moment alone
does not give the sampling variance of a sample *standard deviation*,
and the expression is exact under no distribution. What has a closed
form is one level down, the sample *variance*:
`n^{-1}(mu_4 - (n-3)/(n-1) sigma^4)`.

The text now separates what is exact from what is approximate and uses
only the direction: the kurtosis is `3 - 2 sum a^4/(sum a^2)^2`, below
normal, so the approximation overestimates at leading order and does so
by a cell-dependent amount. "Only the three-element cells clear the
threshold" is therefore a statement about a common threshold, not about
the cells. No exact interval is claimed.

The same reader found both the original defect and the overreach in its
repair, from one file and one question.

### Proposition 17 was not a statement

Two readers reached that independently. The proposition carried no
displayed formula, and its error terms and its object were defined in
the middle of its own proof. An earlier line also drew a *derivative*
limit out of a convergence: convergence to a nonzero limit gives
`D_c(2B)/D_c(B) -> 1`, not a derivative, and `D_c` is defined per band.
The proposition now fixes `Q`, the bands, `l`, the cells, `A`, `C_2`
and `kappa`, displays its estimate, states the uniformity of its
constants, and gives its conclusion as a one-step exponent. The proof
says in as many words that no derivative is taken anywhere.

### Three smaller ones

"The correct scale" overclaimed, since the paper itself reports the
multiplicative sign ensemble as wider and declines to call the
statistic calibrated; it now reads "the exact Rademacher-null scale".
"Computed exactly by convolution" conflated an exact identity with its
float64 FFT evaluation, and the two are now separated. The abstract was
377 words against a 200-250 rule and is now 271 with its references
removed.

### What was reported and refused

Moving Supplement S1 out to supplementary material. The paper's own
sections require the roster and the verdicts to sit inside it, and
three repository checks read those counts from this file; moving it
would leave those checks passing on an empty set. What was removed
instead was the audit narrative the roster does not need. Every rule
name, verdict, count and failure reason stayed.

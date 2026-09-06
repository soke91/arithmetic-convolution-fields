# What changed, and why this is a separate directory

`P1-mobius-fixed-class/` is left as it stood so that anything already
pointing at it keeps resolving to the bytes it pointed at. This
directory carries the same paper with corrections made afterwards.

The corrections came from an independent reader outside this project's
model family, given the file and one question — whether it was ready to
submit. Each was checked here, against the source or by calculation,
before being made.

**One thing that reader found is not fixed and is stated below.** It is
in the proof of the paper's main theorem.

## Not fixed: the truncation error in the density lemma

The lemma disposes of a truncation with

> The truncation error is the tail `sum_{d>D_0} 1/phi(m d^2) <<
> 1/(phi(m) D_0)` **and similarly for e**.

The `d` half is right and uniform: `phi(m d^2) >= phi(m) phi(d^2) =
phi(m) d phi(d)`, so the tail is at most `(1/phi(m)) sum_{d>D_0}
1/(d phi(d)) << 1/(phi(m) D_0)` with an absolute constant.

The `e` half is not the same. There `e | m`, so the sum is finite and
the omitted terms are the divisors of `m` above `E_0`, of which there
are `2^omega(m) = d(m)`. "Similarly for `e`" absorbs that divisor
factor without accounting for it, and the factor is not small against
the truncation: at `log N = 60` the largest `d(m)` exceeds
`E_0 = exp(sqrt(log N))` by an order of magnitude.

Whether the argument survives depends on how the `m`-sum is arranged,
which is not a one-line repair, and it is not attempted here. The main
theorem and the identity downstream of it both rest on this lemma. It
is recorded rather than papered over.

## Fixed: a tail estimate dropped a logarithm

The quantity `Delta` carries `log m` on its short variable. The tail
bound counted only `Lambda << log N` and not that factor, so the
displayed estimate was short by one `log N`. With `log m <= log alpha
<< log N` restored the exponent moves from `A+3` to `A+2`, which is
still `<< N (log N)^{-A}`, so the conclusion holds.

## Fixed: a statement about triply well-factorable sequences was false

The literature section said such a sequence of level `N^theta` vanishes
at every modulus above `N^{theta/3}`. It does not: three factors of
that size have product `N^theta`. The vanishing is a statement about
*prime* moduli, which is what the argument uses and what the companion
paper proves. The sentence now says so and points at that lemma.

## Reported and not a defect

The reader also reported that this directory and the Lean file are
absent from the repository. They are not: the paper points here by URL
and the files resolve anonymously. The reader had a working checkout in
front of it, which is a different tree.

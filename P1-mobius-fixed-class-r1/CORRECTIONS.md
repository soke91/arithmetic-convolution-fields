# What changed, and why this is a separate directory

`P1-mobius-fixed-class/` is left as it stood so that anything already
pointing at it keeps resolving to the bytes it pointed at. This
directory carries the same paper with corrections made afterwards.

The corrections came from an independent reader outside this project's
model family, given the file and one question — whether it was ready to
submit. Each was checked here, against the source or by calculation,
before being made.

The heaviest thing that reader found was in the proof of the paper's
main theorem. It is repaired, and the first version of this file said
it was not; that mistake is described below because it is the more
useful half of the story.

## The truncation error in the density lemma

The lemma disposed of a truncation with

> The truncation error is the tail `sum_{d>D_0} 1/phi(m d^2) <<
> 1/(phi(m) D_0)` **and similarly for e**.

The `d` half is right and uniform: `phi(m d^2) >= phi(m) phi(d^2) =
phi(m) d phi(d)`, so the tail is at most `(1/phi(m)) sum_{d>D_0}
1/(d phi(d)) << 1/(phi(m) D_0)` with an absolute constant.

The `e` half is not the same, and "similarly for e" was wrong. There
`e | m`, so the sum is finite and the omitted terms are the divisors of
`m` above `E_0`, of which there are `d(m)`. A divisor factor appears
that the `d` tail does not have. The lemma now states that half
separately, as

    sum_{e | m, e > E_0} 1/phi(m lcm(d^2,e))  <<  d(m) loglog M / (phi(m) E_0)

and deliberately does *not* absorb it, so that the place it is used has
to carry the sum over `m`.

**Why it is harmless there, and only there.** The truncation enters as
`sum_{m<M} |T_m| ×` (the lemma's error), so what appears is
`sum_{m<M} d(m)/phi(m)`, which is `~ 1.43 (log M)^2` — polylogarithmic.
Since `E_0^{-1} = exp(-sqrt(log N))` beats every power of `log N`, the
`e` half is `<< N (log M)^2 loglog M exp(-sqrt(log N))`, still
`<< N exp(-c sqrt(log N))` for every `c < 1`. It has two powers less
room than the `d` half and still has room.

Against a *single* `m` it would not be harmless: `d(m)` is `m^{o(1)}`
and at `log N = 60` the largest `d(m)` exceeds `E_0` by an order of
magnitude. The lemma says so.

## The mistake in the first version of this file

The first version of this file said the gap was not repairable in a
line, and gave that `d(m)`-against-`E_0` comparison as the reason. The
comparison is correct and irrelevant: **the maximum of `d(m)` is not
what enters the argument.** The sum over `m` does, and a divisor factor
that is `m^{o(1)}` pointwise is `(log M)^2` on average.

This repository has a standing rule that a diagnosis must be sized
before it is carried — compute what the named effect would actually be
and compare it to what was observed. The rule was written here, and
broken here, in the same session: the size that was computed belonged
to a quantity the argument never forms.

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

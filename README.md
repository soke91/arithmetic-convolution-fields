# Arithmetic convolution fields

Papers on the fluctuation structure of convolution sums of arithmetic
functions, with the code and output behind every number they print.

One directory per paper. Each is self-contained: source, compiled PDF,
the scripts the paper cites, the output those scripts produced, and any
machine-checked development. Nothing is shared between papers, so a
directory can be read, run, or moved on its own.

## Papers

| directory | paper | status |
|---|---|---|
| [`P4-coherent-cell-floor/`](P4-coherent-cell-floor/) | Cell-mean fluctuations in an arithmetic convolution field: exact Rademacher-null variances and a permutation control | submitted; **frozen** |
| [`P4-coherent-cell-floor-r1/`](P4-coherent-cell-floor-r1/) | the same paper, corrected after submission | corrections not yet sent to the journal |
| [`P1-mobius-fixed-class/`](P1-mobius-fixed-class/) | A fixed class for the Möbius side | superseded by `-r1` |
| [`P1-mobius-fixed-class-r1/`](P1-mobius-fixed-class-r1/) | the same paper, corrected | not submitted |
| [`P2-no-go-divisor-switch/`](P2-no-go-divisor-switch/) | A no-go for the divisor-switch route | not submitted |

**Why the `-r1` directories exist rather than an edit in place.** One of
these papers is with a journal, and a paper under review has to stay
citable as the thing that was reviewed; the other has been public long
enough that a link to it should keep resolving to what it resolved to.
So corrections go in a new directory and the old one is left alone.
Each `-r1` carries a `CORRECTIONS.md` saying what changed, why, and --
where a reported defect turned out not to be one -- why not.

**One of the corrections carries its own mistake, on purpose.** A gap
found in the fixed-class paper's main proof was first recorded as
unrepairable, on a size estimate that turned out to measure a quantity
the argument never forms. It is repaired, and its `CORRECTIONS.md`
keeps the wrong reasoning next to the right one, because the way the
estimate went wrong is more reusable than the estimate.

A paper with a machine-checked component has a `lean/` directory and an
axiom report inside it; a paper without one says so in its own appendix.
What is machine-checked is in every case the finite identities a proof
runs on, never the theorem the paper argues for, and each paper says
where its formalisation stops and why. No paper here has had a
reading by a subject expert, and each says that too.

## The object

Write

```
C(N) = Σ_{n<N} Λ(n) μ(N−n)      V(N) = Σ_{v<N} μ²(v) Λ(N−v)²
Z(N) = C(N) / √V(N)
```

The field is `Z` over even `N`. The questions here are about how its
means fluctuate when the `N` are grouped — by divisibility, at random,
or under a relabelling — and about what reference scale such a
comparison needs. Where a scale can be computed exactly it is computed
exactly; where it cannot, what is measured is said to be measured.

No distributional law for `Z` is asserted anywhere in this repository,
and no Goldbach consequence is drawn from anything in it. The repository
is named for the object studied, not for a target.

## Running the code

```
python -m pip install -r requirements.txt
cd <paper directory>
python code/<script>.py
```

Python 3.12 and NumPy. Each paper's own README gives the run order, the
environment the printed figures came from, and which scripts compute in
exact rational arithmetic rather than double precision.

**A nonzero exit status is not a failure.** These scripts carry
preregistered rules, and a run whose own rule was refuted says so in its
exit code as well as in its output. The refutations are reported in the
papers, including the ones that went against what was expected.

## What a paper directory contains

```
<paper>.tex        source, compiles standalone
<paper>.pdf        the compiled paper
README.md          run order, environment, what each script produces
code/              every script the paper cites, plus what those read
results/           the output committed with the paper
lean/              formalisation and its axiom report, where there is one
```

`code/` is closed under what it consumes: no script reads a file that is
not in the directory. That is checked by running every script in
dependency order, not by inspection.

## Provenance

These papers are projections of a working repository that is not public.
What that repository keeps and these directories do not is the history —
corrections, withdrawn figures, the record of how each result was
arrived at. What is here is what is true, together with the evidence for
it and the preregistered checks that failed.

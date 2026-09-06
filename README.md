# Arithmetic convolution fields

Papers on the fluctuation structure of convolution sums of arithmetic
functions, with the code and output behind every number they print.

One directory per paper. Each is self-contained: source, compiled PDF, the
scripts the paper cites, the output those scripts produced, `PACKET.json`
pinning every one of those files by SHA-256, and any machine-checked
development. Nothing is shared between papers, so a directory can be read,
run, or moved on its own.

## Papers

| directory | paper | status |
|---|---|---|
| [`P4-coherent-cell-floor/`](P4-coherent-cell-floor/) | Cell-mean fluctuations in an arithmetic convolution field: exact Rademacher-null variances and a permutation control | **submitted** to Experimental Mathematics on 2026-09-06; **frozen** |
| [`P4-coherent-cell-floor-r1/`](P4-coherent-cell-floor-r1/) | the same paper, corrected after submission | corrections not yet sent to the journal |
| [`P1-mobius-fixed-class/`](P1-mobius-fixed-class/) | The Huang–Li reduction of Goldbach to Elliott–Halberstam does not decompose, with an unconditional bound for a Möbius-weighted correlation sum in fixed residue classes | not submitted |
| [`P2-no-go-divisor-switch/`](P2-no-go-divisor-switch/) | No weight that Bombieri–Vinogradov reaches extracts the Möbius-twisted von Mangoldt correlation: a no-go for the divisor-switch route | not submitted |

**Why one paper has an `-r1` directory and the others do not.** The cell-mean
paper is with a journal, and a paper under review has to stay citable as the
thing that was reviewed, so its corrections go in a new directory and the
submitted one is left alone. Nothing else here has been submitted, so
corrections to those go in place, as corrections to an unsubmitted preprint
should. A `CORRECTIONS.md` says what changed, why, and — where a reported
defect turned out not to be one — why not.

`-r1` was itself revised on 2026-09-06 rather than joined by an `-r2`: no link
to it had been given out, the journal had not been told, and the version
standing there carried a correction that was itself wrong. Its
`CORRECTIONS.md` records that and what changed.

## Running the code

From a paper's directory, after `pip install -r requirements.txt`:

```text
python code/<script>.py
```

Some scripts deliberately exit nonzero when a preregistered rule is refuted.
The output, not the exit code, records the verdict, and every refuted rule is
reported in the paper that cites the run.

Every file a `PACKET.json` lists can be checked against the SHA-256 recorded
there. `.gitattributes` disables line-ending conversion so that those hashes
mean the same thing on every platform.

No paper's main theorem should be taken as machine-checked merely because a
`lean/` directory is present; each paper states its own formalization
boundary.

## What is not here

The working repository this is exported from: research records and
retractions, open questions, correspondence, and the literature collection.
What is published is the papers and the evidence they cite.

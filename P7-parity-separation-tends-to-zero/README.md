# P7-parity-separation-tends-to-zero reproduction packet

Purpose: canonical working and public packet for `P7-parity-separation-tends-to-zero.tex`.

Not in this file: ownership judgement or hashes; `PACKET.json` is authoritative for both.

This directory is self-contained for paper reading and evidence reproduction. It
contains 18 Python files under `code/` and 2 under `lib/`, and 92 committed result files. Files
under `code/` and `results/` are generated from the exact sources and hashes in
`PACKET.json`; do not edit those copies.

## Reproduce

The certificates of §11 are **verified** from the files in this directory with no solver (steps 1–4),
and step 3 checks every number §11 prints against them. The four construction bounds (`0.6982` …
`0.5841`) are recomputed exactly from integer class counts only in the optional step 5 (clause `P5`),
which also re-derives certificates from solver output and so calls a solver (SciPy's nonnegative
least squares at `2^16` and `2^18`); in steps 1–4 those four bounds are compared with a stored record. Run the commands from **this
directory** — the scripts locate `results/`, `lib/` and the manuscript relative to themselves, so the
packet is self-contained and nothing points outside it.

**The environment, and how to create it.** `requirements.txt` beside this file pins the versions the
runs reported here were made with:

```text
python -m venv .venv
.venv\Scripts\activate                 # Windows
. .venv/bin/activate                   # POSIX
python -m pip install -r requirements.txt
```

**Tested with CPython 3.11.0** on Windows 10 (x86-64) with exactly those pins — numpy 2.4.6,
numba 0.68.0, llvmlite 0.50.0, scipy 1.17.1, mpmath 1.4.1, psutil 7.2.2 — in a virtual environment
created that way in a directory outside the repository, which is where the "outside the repository"
runs the records mention were made. Python 3.12 also works: the authors' own machine runs 3.12.10 with
numpy 2.5.3 and numba 0.67.0, and that is the combination the timing record names. Nothing needs a GPU,
and nothing needs network access once the environment exists. Which step needs what: numpy and numba in
every step, mpmath in step 3, scipy in step 5 only; psutil in none of them (the column-generation
script that is not part of these steps records memory with it).

**Threads: one rule.** Each verifier takes `ACF_THREADS` if it is set, defaults to **4**, and caps the
result at **8** — and at numba's own thread count, so a four-core machine gets four. Every figure quoted
here was measured at the default 4, so `ACF_THREADS=4` is what to set to compare timings. Each run
prints the number it used.

**Which steps write which files.** Each of steps 1, 2, 3 and 4 rewrites its own record under
`results/`; step 5 also rewrites the certificates it re-derives. The authors' timing record
`results/pa_certify_timings_frozen.json` is written by **no** step: it is the frozen measurement §11
quotes, and your own runs land beside it.

| step | writes |
|---|---|
| 1 | `results/pa_certify.{txt,json}` |
| 2 | `results/pa_cert_2e30_lower.{txt,json}` |
| 3 | `results/pa_tex_constants_v3.{txt,json}` |
| 4 | `results/pa_certify_big.{txt,json}` |
| 5 | the above, plus `results/pa_cert_2e*.npz` for the cells it re-derives |

**What the steps write besides those records**, because "its own record" is not the whole truth:

* **Bytecode and compiled-kernel caches inside this directory.** Python writes `.pyc` files under
  `code/__pycache__/` and `lib/goldbach/__pycache__/`; numba writes the verifier's compiled kernels
  beside them (`code/__pycache__/pa_certify.k_theta-*.nbi` and `.nbc`, and the same for the other
  three kernels). These are caches: deleting them costs the next run a few seconds of compilation and
  changes no result.
* **A compiled-kernel cache in your temporary directory.** The cell loader sets `NUMBA_CACHE_DIR` to
  `numba_cache_acf` under `%TEMP%` (`$TMPDIR` elsewhere) unless you set that variable yourself, and its
  row-count kernel is cached there.
* **One further JSON in the temporary directory, in step 3 only.** The constants checker re-runs the
  older checker it supersedes, and redirects that run's JSON to
  `pa_tex_constants_v2_rerun_by_v3.json` under `%TEMP%`, so that `results/pa_tex_constants_v2.json`
  keeps the authors' run of it. Nothing of step 3 overwrites that stored file.
* **The cells, but only if they are refused, and never under `ACF_NO_CACHE=1`.**
  `code/cell_cache.py` validates `results/cell_cache_2e<e>/` against a fingerprint of the build code
  and of the canonical sieve, and prints `CACHE USED` or `CACHE NOT USED (reason)` on every load. On a
  hit it writes nothing. If it refuses a cache because the files are missing or damaged, it rebuilds
  that cell and **replaces the directory**, which costs minutes at the larger sizes. If instead you set
  `ACF_NO_CACHE=1`, it rebuilds the cell and **returns it without saving**: the shipped directory is
  left exactly as it is, and the next run rebuilds again. That is the setting to use to check the
  shipped cells against a fresh build without touching them. A plain run of the steps below prints
  `CACHE USED` and leaves the cells byte for byte as they are.

Nothing else is written: no file outside this directory and your temporary directory, and no network
access at any point.

**1. Re-check all seven two-sided certificates** (exact integer arithmetic, each band rebuilt):

```text
python code/pa_certify.py --verify
```

Per cell it re-establishes, in this order: the **cell itself** — `N = 2^e`, `Q = ⌊√N⌋`, `q₁` the least
prime above `Q`, `thr = ⌊N/q₁⌋` and the complete list of squarefree `d ≤ Q` are all derived from `e`
alone, and a cache or certificate whose values differ is refused (C0); the inventory is exactly the
seven advertised cells, so a missing certificate fails rather than shrinking the claim (C6); the band
rebuilt and compared against the cached row counts element for element (C1); every stored field present
and well formed (C5); every stored index validated before it indexes anything (C4); both int64 kernels
guarded by the actual integer bound on their accumulator (C2); then `θ_λ ≥ 0` swept over **every**
even-`ω` column of the whole band, `⟨c_𝓞, λ⟩ < 0`, `q ≥ 0`, and both rational bounds of `eq:cert`
recomputed and matched against the stored values. Expected: `FAILS: none`. On the authors' machine (below) a first run takes about 80 s for all seven cells; `results/pa_certify_timings_frozen.json` holds that run's breakdown, the hardware, and what the figures do and do not reproduce.

C0 is what binds the answer to the question. Without it, a larger cell's cache and certificate
relabelled as a smaller cell's pass every other check on their own valid data, and the larger cell's
interval is reported as the smaller cell's.

**No verdict rests on a floating-point square root.** Three places in the shipped code take an integer
bound from a float: the canonical sieve's loop bound `int(limit ** 0.5)` (`lib/goldbach/sieve.py`), and
`code/cell_cache.py`'s `Q = int(N ** 0.5)` and its trial division of `q₁` up to `int(q₁ ** 0.5)`. Two
things make them harmless, and the run states both rather than asking you to assume either.

First, nothing a cache computed is believed: `C0` re-derives `N`, `Q`, `q₁`, `thr` and the complete row
list from `e` alone, in integer arithmetic with `math.isqrt`, and a cache or certificate disagreeing in
any of the five is refused and the cell voided. Those exact integers are what every later clause uses.

Second, clause `C9` settles the float itself, over the whole range this code can sieve: `int(L ** 0.5)`
equals `math.isqrt(L)` for **every** integer `L < 2³¹`. That is checked, not sampled — both maps are
non-decreasing in `L`, so agreeing at the two ends of each `isqrt` level settles the level, and the
46,341 levels cover the range; the check is run both vectorised and with the scalar operator the sieve
itself evaluates.

The three lines are not rewritten with `math.isqrt`, deliberately. `code/cell_cache.py`'s cache
fingerprint hashes the **bytes of the sieve** together with the source of its own build function, so
editing either discards every cached cell in this packet — including `2^28` and `2^30`, which cost
minutes and about 16 GB to rebuild. The check above establishes what the edit would have established,
without invalidating the cells the packet ships.

**2. Re-establish the `2^30` lower endpoint** from a stored integral dual:

```text
python code/pa_certify.py --verify-lower30
```

`results/pa_cert_2e30_lower.npz` holds the **repaired integral** dual `λ` at `2^30` with `rows`, `J`,
`H`, `D_p` and `c₆`. This step rebuilds the `2^30` band, sweeps `θ_λ ≥ 0` over every even-`ω` column of
it, recomputes `J = −⟨c_𝓞,λ⟩`, `H = ‖λ‖²` and `D_p = ‖P‖²` from the stored vector, and checks
`c₆²·H·D_p < J²·10¹²` with `c₆` maximal — which is `κ(2^30) > c₆/10⁶ = 0.017787`. No upper endpoint is
claimed at `2^30`. Four integers satisfying an inequality prove nothing on their own; the stored `λ` and
this sweep are what make them the right integers.

`λ` is the separator rounded at `10¹²`, which is the scale the published integers are at, so the four
integers this step recomputes are the ones §11 prints — `J = 4401612119558212802`,
`H = 2015370564893177546709684150`, `D_p = 30382600959846`, `c₆ = 17787` — digit for digit, and clause
`L4` compares them with the solver run's own ledger `results/conedual_colgen_2e30_n30b.json`. Rounding
the same separator at another scale would certify the same cell with different valid integers; `L4` is
what keeps the stored certificate the one the paper's numbers come from.

On the authors' machine (below): **247.5 s** in the frozen timing record
(`results/pa_certify_timings_frozen.json`, field `verify_lower30_s`) — most of it the band rebuild, then
the sweep. The band at `2^30` is the cost; budget a few minutes and about 16 GB of memory (the memory
figure was observed when the certificate was built and is not part of a shipped record).
`python code/pa_certify.py --make-lower30` rebuilds the certificate itself from the stored separator
(a few minutes), and writes `results/pa_cert_2e30_lower.npz` as well as its own record.

**3. Re-check the numbers §11 prints** against the certificates and the frozen records:

```text
python code/pa_tex_constants_v3.py
```

2.4 s. It tokenises `P7-parity-separation-tends-to-zero.tex` and accounts for every numeric literal
occurrence in the body — decimals, `\,`-grouped integers and integers of five digits or more — by a
clause that recomputes it, by a declared quotation from a cited source, or as a reference inside
`\cite[...]`. It does **not** parse ordinary small integers wherever they appear; what it does check
beyond that inventory is each table row's own exponent label against the cell whose certificate produced
that row's endpoints, and the exponent and `Y` lists §11 names outside the table (K6). The seven
intervals are recomputed from the certificates themselves, as is the `2^30` four-integer inequality
including `c₆`'s maximality and `D_p` recomputed from the `2^30` cell.

Clause `R1` is the one claim here about a measurement rather than arithmetic: it compares the four
timings §11 prints with `results/pa_certify_timings_frozen.json`, the frozen record of the authors' run,
and requires that record to name the SHA-256 of the verifier shipped here. If you have run step 1, your
own times are reported beside them and deliberately not scored.

**Its output and exit code.** This step is not a certificate verifier and prints no `CERTIFICATES:`
line. It ends with a `FAILS:` line listing every clause that fails (`FAILS: none` when all hold) and the
names of the two record files it wrote; its exit code is `0` exactly when every clause holds and `1`
otherwise.

**What this step does not do.** It checks the manuscript's numbers against the certificates and the
stored records -- that every figure §11 prints is the one the stored data gives, under the rounding the
text declares. It does **not** check that a certificate is valid: nothing in it sweeps a dual over the
band or tests a feasibility condition. `FAILS: none` here therefore says that the paper's digits match
the evidence, not that the evidence certifies anything. Step 1 is what establishes the seven two-sided
certificates and step 2 the `2^30` lower endpoint; a number is worth no more than the step that
established it.

**4. The same seven certificates through the wrapper** that produced the two largest:

```text
python code/pa_certify_big.py --verify
```

This step runs **every check step 1 runs** — the cell derivation `C0`, the band rebuild `C1`, the field
and index checks `C5`/`C4`, the accumulator guards `C2`, the per-cell re-verification `V<e>`, the
square-root clause `C9` and the inventory clause (`C8` here, `C6` there) — because the per-cell work is
`pa_certify.py`'s own code called from this file, and `C9` is `pa_certify.py`'s own function called from
it too. The wrapper omits no check of the verifier it wraps; what differs is only which separator and
dual the two largest cells were built from, which is why this file exists.

**5. Re-derive a certificate from its solver output** (optional; the only step that rewrites
`results/pa_cert_2e*.npz`). At `e = 16,18` the primal is re-derived by a dense nonnegative least squares
over the complete column set and the dual from its residual; at `e = 20,22,24` both come from the stored
exact duals:

```text
python code/pa_certify.py --cells 16 18 20 22 24
```

At `e = 26,28` the primal and the dual come from different runs, so each cell names its own pair:

```text
python code/pa_certify_big.py --cells 26 --sep-tag _pa9 --dual-tag _n0b
python code/pa_certify_big.py --cells 28 --sep-tag _pa9 --dual-tag _n28
```

**What a run reports, and what its exit code means.** Every **verifier** command — `pa_certify.py` and
`pa_certify_big.py`, in steps 1, 2, 4 and 5 — ends with two machine-readable lines:

```text
CERTIFICATES: VALID
DIAGNOSTICS: P3 fails (expected: it compares a nine-decimal rounded float optimum with a twelve-decimal exact interval, which need not contain it -- section 11 says so in a parenthesis)
```

The first line is the verdict on the certificates themselves, and **the exit code reports that line
alone**: `0` when every certification clause holds, `1` when one does not, with the failing clause named
on the line. Certification is what establishes a certificate — the cell is the one `e` names and is
intact (`C0`, `C1`, `C5`, `C6`), the exact integer arithmetic behind the bound is sound and every stored
field is the integer it claims to be (`C2`, `C3`, `C4`, the `V<e>` re-verifications, `L0`–`L4` at
`2^30`), the certificate is feasible and two-sided (`P1`, and `P5` for the construction's own bounds).

The second line lists the clauses that compare a certificate with something *outside* it: a width
target, a number printed elsewhere, the published `eq:certint` integers, a file-size budget. Two of them
fail, and both are expected; they occur only in the optional construction runs of step 5 -- the
mandatory verification runs (steps 1 and 4) report no diagnostic failure. `pa_certify.py`'s `P3` compares the stored *nine-decimal* float
optima `0.172389352` and `0.140860454` for containment in *twelve-decimal* intervals, and a nine-decimal
rounding need not lie inside one — §11 says so in a parenthesis. `pa_certify_big.py`'s `P2` at `2^28`
asks for a width below `10⁻⁶` where the width is `1.2·10⁻⁵`; §11 claims only `10⁻⁴` at that size. Both
are clauses being refuted and recorded, and neither touches a certificate: the certified intervals are
reproduced unchanged either way. The `FAILS:` line above the verdict still names every failing clause of
either kind, and the JSON record carries `certificates_valid`, `certification_fails` and
`diagnostic_fails`.

**Verification mode and construction mode make different claims, and each run says which it is in.**
Steps 1, 2 and 4 are **verification**: they open the stored files and re-establish what those files
claim, writing nothing but their own records. Step 5 is **construction**: it re-derives certificates and
writes them, and its output therefore mixes the certification clauses with the construction's own
diagnostics (`P2`–`P6` here are about targets, printed numbers and file sizes, not about the stored
evidence). A construction run opens the cells it was given and runs no completeness clause, so it is not
a packet check and does not claim to be one. Each run announces its mode near the beginning of its output, after the title, the thread
setting and clause `C9` (`VERIFY MODE …` / `CONSTRUCTION MODE …`).

**Only a run with no `--cells` speaks for the packet's contents.** The verdict line says what the run
covers, and it is the one place to read that from:

```text
CERTIFICATES: VALID                                                   step 1 or 4, no --cells
CERTIFICATES: VALID for 2^16 only; completeness not checked           --verify --cells 16
CERTIFICATES: VALID for 2^30 (the lower endpoint only); completeness not checked    step 2
CERTIFICATES: VALID for the cells built (2^16, …); completeness not checked         step 5
```

With `--cells` naming anything other than the advertised set, a verification run opens the selected
certificates and **no others** — a damaged certificate outside the selection is never read — so the
completeness clause (`C6` in `pa_certify.py`, `C8` in `pa_certify_big.py`) reports **`NOT RUN`** instead
of passing, a `NOT RUN:` line appears beside `FAILS:`, and the record carries
`completeness_checked: false` with `clauses_not_run`. A line beginning `CERTIFICATES: VALID` therefore
means "valid, for what this run examined", and the rest of the line says what that was. A scoped run is not a check of the packet.

At `e = 26,28` a re-derived certificate replaces the stored one only if its enclosure lies **inside**
the stored one and is strictly smaller on at least one side, compared as the exact squared rationals the
files store -- not as rounded widths. So a replacement can never lose any part of the stored enclosure,
and two enclosures that merely cross (a better lower endpoint bought with a worse upper one) leave the
stored certificate in place. The run prints both enclosures and which one it kept. At `e ≤ 24` it is written unconditionally; on
the run that produced this packet all five came back byte-identical, but they need not: the dense least
squares at `2^16` is reproducible within one process while its residual can differ by one unit in the
last place between runs, which moves one coordinate of the rounded dual. The certified interval is
stable under that, so step 3 still passes afterwards.

The column-generation runs that produced the `e = 26,28,30` separators are hours on eight threads and
are not part of this packet; their records are `results/conedual_colgen_2e{26,28}_pa9.*` and
`results/conedual_colgen_2e30_n30b.*`.

## The authors' machine

The timings quoted above and in §11 were measured on:

| | |
|---|---|
| CPU | 11th Gen Intel Core i5-11500 @ 2.70 GHz, 6 physical / 12 logical cores |
| RAM | 63.9 GB |
| OS | Windows 10 Pro, build 10.0.19045 |
| Python | 3.12.10, numpy 2.5.3, numba 0.67.0 |
| threads | 4 (`ACF_THREADS=4`) |

Times on other hardware differ; the band rebuild dominates every figure.

## What backs each §11 claim

| claim | files |
|---|---|
| the seven certified intervals | `results/pa_cert_2e{16,18,20,22,24,26,28}.npz` — each the two integral vectors plus `D` and the four big rationals |
| the `2^30` lower endpoint | `results/pa_cert_2e30_lower.npz` — the repaired integral dual `λ` with `rows`, `J`, `H`, `D_p`, `c₆`, rounded at the scale the published integers are at; re-established by step 2, which also compares the four integers with the solver ledger |
| the `2^30` band sizes \|𝓔\|, \|𝓞\| and the solver's uncertified value | `results/conedual_colgen_2e30_n30b.json` and its `.txt` transcript, with `results/conedual_colgen_sep_2e30_n30b.npz`, the float separator the integral dual was rounded and repaired from |
| the cells every check is run against | `results/cell_cache_2e{16,18,20,22,24,26,28,30}/` — `rows.npy`, `cb.npy`, `cr.npy`, `rom.npy`, `meta.json`. Clause C0 derives `N`, `Q`, `q₁`, `thr` and the complete row list independently and refuses a cache that disagrees, so these are checked, not trusted |
| the verifier | `code/pa_certify.py`, `code/pa_certify_big.py` |
| the constants checker | `code/pa_tex_constants_v3.py`, and `code/pa_tex_constants_v2.py`, which it imports and re-runs |
| their library imports | `code/cell_cache.py` and `code/_sieve_shared.py`, over `lib/goldbach/sieve.py` — the canonical sieve, whose bytes the cache fingerprint hashes |
| the timings §11 prints | `results/pa_certify_timings_frozen.json` — frozen, written by no step, with the hardware above. `results/pa_certify.{txt,json}` and `pa_certify_big.{txt,json}` are where your own runs land |
| the construction's `0.6982, 0.6130, 0.6088, 0.5841` | `code/pa_certify.py` step 5, whose clause `P5` prints the exact rational bound at each of the four sizes, rounds it **upward** at four decimals and reports *valid and tight*. Also `results/conedual_parity_domination.{txt,json}`, the float sweep the sizes came from — note it prints `0.6981` and `0.6087` at `2^16` and `2^20`, which are those same bounds rounded to **nearest**; the paper rounds upward, deliberately, because only an upper bound is a bound |
| the `e ≤ 20` float optima and the `e = 22…28` ledgers the older clauses read | `results/conedual_separator_mass.txt`, `conedual_exact_kappa_small_cells{,_20}.txt`, `conedual_colgen_2e{22,24,26,28}_{n0b,n28}.{txt,json}` |
| re-deriving the certificates (step 5) | `results/conedual_exact_dual_2e{20,22,24}.npz`; `results/conedual_colgen_sep_2e26_{pa9,n0b}.npz` and `conedual_colgen_sep_2e28_{pa9,n28}.npz` — the frozen pairing §11 records, primal and dual from different runs at the two largest cells |
| the manuscript the checker reads | `P7-parity-separation-tends-to-zero.tex` in this directory. The checker looks here first and consults nothing outside this directory |
| the Lean development | `lean/` — from that directory, `lake env lean KappaZero.lean` with the toolchain `lean-toolchain` names; `lean/axioms.txt` is the stored `#print axioms` output |
| provenance of every file above | `PACKET.json`: the source path and SHA-256 of each artifact, which `tools/packets.py check` enforces in the repository this packet is built from |

A release tag for this packet must name a commit that **contains this directory as shipped**;
`source_snapshot` and the artifacts' `source_commit` record where the artifact sources were read from in
the working repository and are not that tag.

Some scripts deliberately exit nonzero when a clause is refuted. The output, not a blanket
zero-exit convention, records that verdict.

**Transcript note.** In the seven column-generation transcripts `results/conedual_colgen_2e*.txt` one printed line naming an unpublished working document was shortened to what the current code prints (`clauses K0-K4; not reimplemented here`), and one stored line of `results/gramint_kappa_certify_ledger.json` names the machine generically. No number in any file was changed.

## Provenance

Artifact-level source paths and SHA-256 digests are in `PACKET.json`.

## Formalization

Status: `library-in-packet`. A standalone Lean 4 / mathlib v4.33.1 project (lean/), sourced from generations/v3/lean_kappa/: 29 theorems on the exact finite steps (one-large-prime identity, identity (E), Bonferroni, parity split, cone step), no sorry, axioms in lean/axioms.txt. The analytic estimates and the main theorem are not formalized.

The presence of a Lean directory does not mean the paper's main theorem is
machine-checked; consult the manuscript appendix and the files themselves for
the formalization boundary.

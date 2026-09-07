# P2-no-go-divisor-switch reproduction packet

Purpose: canonical working and public packet for `P2-no-go-divisor-switch.tex`.

Not in this file: ownership judgement or hashes; `PACKET.json` is authoritative for both.

This directory is self-contained for paper reading and evidence reproduction. It
contains 16 Python files and 15 committed result files. Files
under `code/` and `results/` are generated from the exact sources and hashes in
`PACKET.json`; do not edit those copies.

## Reproduce

From this directory:

```text
python -m pip install -r ../requirements.txt
python code/<script>.py
```

Some scripts deliberately exit nonzero when a preregistered rule is refuted.
The output, not a blanket zero-exit convention, records that verdict.

## Provenance

The manuscript projects these working sources:

- `generations/v2/paper/theorem_A.md`
- `generations/v2/paper/wall_v3.md`
- `generations/v3/paper/theorem_A.md`

Artifact-level source paths and SHA-256 digests are in `PACKET.json`; generation
names alone never select an artifact.

## Formalization

Status: `library-in-tree`. The Lean development is carried at lean/ and every file here is generated from it. Its own upstream is the standalone repository z:/업무/goldbach-lean, imported 2026-09-06. The presence of these files does not mean the paper's main theorem is machine-checked.

The presence of a Lean directory does not mean the paper's main theorem is
machine-checked; consult the manuscript appendix and the files themselves for
the formalization boundary.

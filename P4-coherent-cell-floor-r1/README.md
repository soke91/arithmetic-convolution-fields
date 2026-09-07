# P4-coherent-cell-floor reproduction packet

## Where this appears publicly

This packet publishes to **`P4-coherent-cell-floor-r1/`** in the public repository, not to a
directory of its own name. P4-coherent-cell-floor/ in the public repository holds the bytes sent to Experimental Mathematics and is frozen; the working manuscript publishes to the revision directory beside it. Chosen 2026-09-06 to replace the first revision rather than add a second, because no link to it had been given out and it carried a correction that was itself wrong.

There is no `P4-coherent-cell-floor-r1` directory in this tree and there is not meant to be:
the public name is a property of publication, and this directory is the one
place the content lives.

Purpose: canonical working and public packet for `P4-coherent-cell-floor.tex`.

Not in this file: ownership judgement or hashes; `PACKET.json` is authoritative for both.

This directory is self-contained for paper reading and evidence reproduction. It
contains 16 Python files and 16 committed result files. Files
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

- `generations/v2/paper/wall_v3.md`

Artifact-level source paths and SHA-256 digests are in `PACKET.json`; generation
names alone never select an artifact.

## Formalization

Status: `library-in-tree`. The Lean development is carried at lean/ and every file here is generated from it. Its own upstream is the standalone repository z:/업무/goldbach-lean, imported 2026-09-06. The presence of these files does not mean the paper's main theorem is machine-checked.

The presence of a Lean directory does not mean the paper's main theorem is
machine-checked; consult the manuscript appendix and the files themselves for
the formalization boundary.

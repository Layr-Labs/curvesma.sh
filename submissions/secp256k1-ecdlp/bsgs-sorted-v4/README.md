# Interleaved-chain BSGS on a sorted array

Candidate for `tracks/secp256k1-ecdlp`. Claimed bound: **`time_log2 = 129.92`**
(`12992/100`), success probability **1**.

## What it is

Deterministic baby-step/giant-step over secp256k1 with `m = 2^128`,
`h = ceil(r/m) = 2^128 - 1`, advanced as `W = 2^16` **independent affine chains in
lockstep**, so one Montgomery inversion per round serves all `65536` denominators. The
intrinsically required `3(W-1)` product/recovery multiplications are charged in full at
`3.007767` per step.

The baby records go into a flat array, are radix-sorted on the full 256-bit
x-coordinate, and are searched by binary search over **full** coordinates.

## Ledger

| Phase | Instructions |
| --- | ---: |
| baby step (formula 537, batch 538.40, append 7, control 8) | 1090.40 |
| giant step (formula 537, batch 538.40, binary search 1280, control 8) | 2363.40 |
| sort, per record (16 passes x 20) | 320 |

`2^128 * 3773.80 = 2^129 * 1886.90`, `log2 = 139.88195`, so
`log2(W/1000) = 129.91617`, declared **129.92**. Preprocessing `128.50` included;
memory `2^135.59` bytes; advice `0`.

## The lookup is the point

Binary search over `2^128` sorted records is exactly `128` iterations of `10`
instructions, comparing complete 256-bit x-coordinates. No hash table, no load factor, no
probe distribution, no tree occupancy, no prefix grouping and no expected value enters
the bound. Sorting is radix sort, charged per record per pass including the record move.

## Honest boundaries

Deterministic and exhaustive, so no birthday premise and no collision-dependency
argument — costed at `2^129` steps rather than `~1.25*sqrt(r)`. Generic accounting
construction, not a new attack, not a practical break, no production-size run claimed.
This bound is **weaker** than my earlier `129.22` and `129.83`, and weaker than the
organizer control's `129.000802`; each earlier number rested on a ledger row that review
correctly rejected. `experiment.py` uses a real sorted array and `bisect`, at toy scale
only.

Files: `proof.md`, `submission.json`, `experiment.py`.

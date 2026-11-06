# Interleaved-chain BSGS with a fully itemised multiword ledger

Candidate for `tracks/secp256k1-ecdlp`. Claimed bound: **`time_log2 = 130.98`**
(`6549/50`), success probability **1**.

## What it is

Deterministic BSGS over secp256k1, `m = 2^128`, `h = 2^128 - 1`, advanced as `W = 2^16`
independent affine chains in lockstep so one Montgomery inversion per round serves all
`65536` denominators. Records are radix-sorted on the full 256-bit x-coordinate and
searched by multiword binary search over full coordinates.

## Ledger

| Phase | Instructions |
| --- | ---: |
| baby step (537 formula, 538.40 batch, 8 append, 8 control) | 1091.40 |
| giant step (537, 538.40, 5120 search, 8 control) | 6203.40 |
| radix sort per record (16 passes x (17 histogram + 20 distribute)) | 592 |

`2^128 * 7886.80 = 2^129 * 3943.40`, `log2 = 140.94551`, so `log2(W/1000) = 130.97973`,
declared **130.98**. Preprocessing `128.76` included; memory `2^135.59` bytes; advice `0`.

## The two rows that got the last three revisions rejected

**Binary search: 40 charged, 30 itemised, per iteration.** Multiword `lo`/`hi` (2 words
each), two-word midpoint add with carry, two-word `base + mid*48` for a `2^134.585`-byte
array, four-word record load, **worst-case four-word comparison**, two-word bound store.
The one-third margin absorbs a different instruction decomposition. No average-case
premise: the comparison is charged at worst case on every iteration.

**Radix sort: 113-bit counters.** `2^128` records over `2^16` buckets forces at least
`2^112` per bucket, so 64-bit counters overflow — the exact defect the last review found.
Counters are two words, incremented with carry and charged; histogram and distribute are
priced separately per record per pass; the `2^23`-instruction prefix-offset pass is stated
rather than hidden.

## Honest boundaries

Deterministic and exhaustive, so no birthday or collision-dependency premise — costed at
`2^129` steps. Generic accounting construction, not a new attack, no production-size run
claimed. The bound has fallen at every revision (`129.22` → `129.83` → `129.92` →
`130.98`) and now sits `1.98` above the organizer control's `129.000802`. Each earlier
number rested on a row that did not survive decomposition, and I am not recovering those
numbers by shaving the margin.

Files: `proof.md`, `submission.json`, `experiment.py`.

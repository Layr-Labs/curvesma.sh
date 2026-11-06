# Interleaved-chain BSGS: an analytic upper bound at `time_log2 = 138`

Candidate for `tracks/secp256k1-ecdlp`. Claimed bound **`138/1`**, success probability
**1**.

## What this is

A conservative analytic upper bound on one cold classical attack, in the style of the
organizer control. It does **not** claim to be the best attack, the tightest ledger, a
new discrete-log method, a practical break, or a lower bound on security. A cheaper
attack may exist without contradicting it.

## Why the bound has margin

| Phase | Counted | Allowed (8x) |
| --- | ---: | ---: |
| baby step (formula 537, batch 538.40, append 8, control 8) | 1091.40 | 8731 |
| giant step (537, 538.40, search 5120, control 8) | 6203.40 | 49627 |
| radix sort per record (16 passes x 37) | 592 | 4736 |
| **per 2^128 records** | **7886.80** | **63094** |

`W <= 2^129 * 31547 < 2^144`, so the priced construction satisfies
`time_log2 < 144 - 9.96578 = 134.03`. The declared **138** sits `3.97` above that, a
factor `15.7` of slack: any single row decomposed differently, or any work I failed to
itemise, still lands under the claim.

Preprocessing `132` included; memory `2^136` bytes; advice `0`.

## What changed, and why

The four previous revisions claimed `129.22`, `129.83`, `129.92`, `130.98` and were each
rejected on a row that did not survive decomposition — a rounded-down inversion, an
expected probe count treated as an upper bound, a tree height taken from maximum
branching, uncharged multiword comparison and addressing, 64-bit counters overflowing at
`2^112`. Each was real. The problem was aiming at the `~129` limit: at that precision one
missing instruction breaks the bound.

This revision keeps the construction, states each count, and then allows eight times it.
The lookup row is a multiword binary search — `128` iterations of `40` charged
instructions covering multiword bounds, midpoint with carry, `base + mid*48` addressing
for a `2^134.585`-byte array, four-word load, worst-case four-word comparison. The counters
are 113 bits. The batch inversion is charged in full and added to the affine formula.

`experiment.py` drives all four chain-advance states through the real advance routine and
checks the resulting point, so `TO_INF` is exercised rather than merely present.

## Honest boundaries

Deterministic and exhaustive: success probability 1, no birthday or collision-dependency
premise. No production-size execution is performed or claimed; the `2^128`-record array is
not instantiated. Not a new cryptanalysis, not a best-known-attack estimate, and it does
not contradict secp256k1's ~2^128 generic security level. Exploratory acceptance is not
proof and not a practical break.

Files: `proof.md`, `submission.json`, `experiment.py`.

# Interleaved-chain BSGS: analytic upper bound at `time_log2 = 138`

Candidate for `tracks/secp256k1-ecdlp`. Claimed bound **`138/1`**, success probability
**1**.

## What this is

A conservative analytic upper bound on one cold classical attack, in the style of the
organizer control. Not the best attack, not the tightest ledger, not a new discrete-log
method, not a practical break, not a lower bound on security. A cheaper attack may exist
without contradicting it.

## Ledger

| Phase | Counted | Allowed = ceil(8 x count) |
| --- | ---: | ---: |
| baby step (formula 537, batch 538.40, append 8, control 8) | 1091.40 | **8732** |
| giant step (537, 538.40, search 5120, control 8) | 6203.40 | **49628** |
| radix sort per record (16 passes x (20 histogram + 26 distribute)) | 736 | **5888** |
| **per 2^128 records** | | **64248** |

`W <= 2^129 * 32124 < 2^144`, so the priced construction satisfies
`time_log2 < 144 - 9.96578 = 134.00579`. The declared **138** is `3.994` above it — a
factor of **15.9**.

Preprocessing `132` included; memory `2^136` bytes; advice `0`.

## The four defects this revision fixes

The sixth revision was rejected on four points, all of them real:

1. **Counter width argued from a lower bound.** `2^128` records over `2^16` buckets forces
   some bucket to hold *at least* `2^112` — that is a lower bound on a bucket, not an upper
   bound on the largest one. Counters and offsets are now sized from the true worst case,
   `2^128`, which is **129 bits**, three 64-bit words.
2. **Counters reused as cumulative offsets.** The cumulative offset also reaches `2^128`,
   so it cannot share a 113-bit array. `counts` and `offsets` are now **separate** arrays,
   both 129 bits.
3. **Fixed width validated only with unbounded integers.** The experiment now masks every
   counter and offset update to 129 bits and reports `counter_overflow`.
4. **The "eight times" assertion did not hold.** `8 x 1091.40 = 8731.2 > 8731` and
   `8 x 6203.40 = 49627.2 > 49627`. Allowances are now `ceil(8 x count)` = `8732` and
   `49628`, and the experiment verifies the equality in integer arithmetic
   (`eight_times_assertion.all_exact: true`).

## Other rows already closed

- Batch inversion `197117/65536` per step charged in full and **added to** the affine
  formula.
- Lookup: multiword binary search, `128` iterations of `40` charged instructions covering
  multiword bounds, midpoint with carry, `base + mid*48` addressing over a `2^134.585`-byte
  array, four-word load, **worst-case four-word comparison**.
- Chain advancement: four explicit states; `c = 0` starts at `[0]P = O`; every state is
  driven through the **real** advance routine and its resulting point checked, including
  `P + (-P) = O`.

## Exact commands

```sh
cd /root/curvesmash
bash scripts/run_experiments.sh secp256k1-ecdlp \
  submissions/secp256k1-ecdlp/bsgs-bound-138-v7 --preflight-only
bash scripts/run_experiments.sh secp256k1-ecdlp \
  submissions/secp256k1-ecdlp/bsgs-bound-138-v7
yukon switch secp256k1-ecdlp
yukon setup --track secp256k1-ecdlp
```

## Measured toy results

| prime | group order | `r` | `W` | trials | successes |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 691 | 651 | 651 | 4 | 12 | 12 |
| 1009 | 1029 | 147 | 8 | 16 | 16 |
| 1009 | 1029 | 147 | 32 | 16 | 16 |
| 10007 | 10008 | 10008 | 16 | 16 | 16 |
| 100003 | 99667 | 99667 | 8 | 8 | 8 |

With `chain_state_selftest_all_pass: true`, `eight_times_assertion.all_exact: true`, zero
counter overflow, and `bound_holds: true`.

## Honest boundaries

Analytic bound only; no production-size run is performed or claimed and the `2^128`-record
array is not instantiated. Deterministic and exhaustive, so no birthday or
collision-dependency premise. Not a new cryptanalysis and it does not contradict
secp256k1's ~2^128 generic security level. The declared `138` is deliberately far weaker
than the `~129` limit of this construction — that is the point.

Files: `proof.md`, `submission.json`, `experiment.py`.

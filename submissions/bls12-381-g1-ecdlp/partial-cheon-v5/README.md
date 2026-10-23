# partial-cheon-v5: bls12-381-g1-ecdlp candidate

Fifth revision of the partial-coverage Cheon attack, fixing the
fourth-round objections: uncharged sort/lookup work, unaccounted merge
buffer, and the equal-x sign ambiguity.

## Summary

| Field | Value |
| --- | --- |
| time_log2 = log2(W(2^20)) | **123.831402** (= 61915701/500000) |
| success_probability (lower bound) | 0.39 (exact per-offset minimum (m-1)m/N >= 0.4 - 2^-118) |
| memory_log2_bytes | 124.8 (main array + merge buffer) |
| preprocessing_log2 | 21.0 (included in total) |
| nonuniform_advice_log2_bytes | 0 |

## What changed from v4

- Every sort/lookup instruction is charged: merge sort = 117 passes x
  (50 move + 15 compare) instructions per baby = 3.8025 units; binary
  search = exactly 118 probes x 25 instructions = 1.475 units per
  lookup. The comb charge (64 units) now states its underlying count
  (32-40 window additions). The score rose from 123.792992 to 123.831402
  to absorb exactly these charges.
- The merge buffer is allocated and reported: memory is the m-entry
  main array plus the m-entry buffer, strictly above 2^124.77 bytes,
  reported as 124.8.
- Records are keyed by (x, y mod 2) - unique per point on this odd-prime
  curve (same-x pairs have opposite parities) - so a key match is an
  exact point match and k1 = (i + m1*j) mod N has no sign ambiguity.
- The toy table is now the parity-keyed sorted array built by the
  explicit merge sort and searched by the explicit binary search (probe
  budget asserted per lookup), not a full-point dictionary.
- Part A additionally crosschecks the closed-form coverage bound against
  an exhaustive every-offset brute force on five smaller parameter sets.

## Reproduce

```sh
echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v5/experiment.py
echo '{"trials": 200, "mode": "full"}'    | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v5/experiment.py
```

Observed locally: Part A all checks pass (per-offset minimum >= 39/100
in integers; five exhaustive crosschecks; parity-key uniqueness); Part B
full mode 200/200, 0 wrong extractions; partial mode 88/200 = 0.44 vs
the exact wrap-inclusive expectation 0.40226, 0 wrong extractions.

Small-scale experiments are evidence at their own scope only.

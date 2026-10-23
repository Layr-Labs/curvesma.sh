# partial-cheon-v6: bls12-381-g1-ecdlp candidate

Sixth revision of the partial-coverage Cheon attack, fixing the
fifth-round objections: exact-count claims contradicted by the
implementation, memory lower-bound arithmetic, and incomplete
constructions for the charged allowances.

## Summary

| Field | Value |
| --- | --- |
| time_log2 = log2(W(2^20)) | **123.834773** (= 123834773/1000000) |
| success_probability (lower bound) | 0.39 (exact per-offset minimum (m-1)m/N >= 0.4 - 2^-118) |
| memory_log2_bytes | 124.8 (main array + merge buffer) |
| preprocessing_log2 | 21.0 (included in total) |
| nonuniform_advice_log2_bytes | 0 |

## What changed from v5

- Every sort/lookup instruction is charged: merge sort = 117 passes x
  (50 move + 15 compare) instructions per baby = 3.8025 units; binary
  search = exactly 118 probes x 25 instructions = 1.475 units per
  lookup. The comb charge (64 units) now states its underlying count
  (32-40 window additions). The score rose from 123.792992 to 123.831402
  to absorb exactly these charges.
- The memory lower bound is stated correctly: the two arrays are
  exactly 256*m = 2^124.767581 bytes; with comb tables, code and
  scratch the peak is strictly above 2^124.7675 bytes, reported
  rounded up as 124.8 (v5's "strictly above 2^124.77" was
  arithmetically wrong).
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
echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v6/experiment.py
echo '{"trials": 200, "mode": "full"}'    | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v6/experiment.py
```

Observed locally: Part A all checks pass, including the executed
Brent-rho factorization of the production r-1 (matches the supplied list
exactly) and the per-offset minimum >= 39/100 in integers; Part B full
mode 200/200, 0 wrong extractions; partial mode 88/200 = 0.44 vs the
exact wrap-inclusive expectation 0.40226, 0 wrong extractions, sort
compare/move upper bounds and the 9-probe toy budget asserted per run.

Small-scale experiments are evidence at their own scope only.

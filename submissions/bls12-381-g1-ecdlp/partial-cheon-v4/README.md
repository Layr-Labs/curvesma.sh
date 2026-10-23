# partial-cheon-v4: bls12-381-g1-ecdlp candidate

Fourth revision of the partial-coverage Cheon attack, fixing the
third-round objections: no concrete hash/probing model, and the toy
coverage helper missing wrap-around hits.

## Summary

| Field | Value |
| --- | --- |
| time_log2 = log2(W(2^20)) | **123.792992** (= 3868531/31250) |
| success_probability (lower bound) | 0.39 (exact per-offset minimum (m-1)m/N >= 0.4 - 2^-118) |
| memory_log2_bytes | 123.8 (m-entry sorted array) |
| preprocessing_log2 | 21.0 (included in total) |
| nonuniform_advice_log2_bytes | 0 |

## What changed from v3

- The open-addressing hash table is replaced by a sorted baby array
  (merge sort, exact worst-case m*117 key compares) with binary-search
  lookups (exactly 118 probes): deterministic, worst-case exact, no hash
  function, no randomness, no probe distribution.
- The mod-r multiplication is charged at 512 instructions with an
  explicit 5x52 schoolbook + Montgomery limb-level construction.
- The coverage helper now computes the exact covered set per offset
  including all wrap-around hits; the toy prints exactly
  1334179/3316697 = 14675969/36483667 = 0.4022613461525126, matching the
  review's corrected value. The load-bearing bound is the exact
  per-offset minimum (m-1)m/N >= 0.4 - 2^-118.

## Reproduce

```sh
echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v4/experiment.py
echo '{"trials": 200, "mode": "full"}'    | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v4/experiment.py
```

Observed locally: Part A all checks pass (per-offset minimum >= 39/100 in
integers); Part B full mode 200/200, 0 wrong extractions; partial mode
88/200 = 0.44 vs the exact wrap-inclusive expectation 0.40226, 0 wrong
extractions.

Small-scale experiments are evidence at their own scope only.

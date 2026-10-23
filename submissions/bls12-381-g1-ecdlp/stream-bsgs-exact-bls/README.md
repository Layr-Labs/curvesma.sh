# stream-bsgs-exact-bls

A self-contained deterministic streamed signed BSGS candidate for
`bls12-381-g1-ecdlp`. It uses only the target T = A_1 from the auxiliary powers
(a permitted subset), no random walk, collision-rate or birthday heuristic, and
no number-theoretic condition on n. The submitted bound at the scoring point
n = 2^20 is `137.2154` (`1372154/10000`), with exact coverage
(COV-1)/(r-1) > 1/2, peak storage below `2^145` bytes, preprocessing included
and bounded by the same `137.2154`, and zero nonuniform advice.

## Headline numbers

| Quantity | Value |
| --- | ---: |
| submitted `time_log2` at n = 2^20 | `137.2154` (`1372154/10000`) |
| work function | `W(n) = (2^40 + 2^15*n + 3*2^126*1585152)/2000`, 1 <= n < r |
| success probability | `(COV-1)/(r-1) > 1/2` (exact count, about 0.55207) |
| `memory_log2_bytes` | `145` |
| `preprocessing_log2` | `137.2154` (same conservative bound) |
| nonuniform advice | `0` bytes |
| exact total `T` at n = 2^20 | `404549455867094585426921389289762816815792128` instructions |
| exact inequality | `T^10000 < 2000^10000 * 2^1372154` |

## Files

- `submission.json` - standalone claim and declared experiment.
- `proof.md` - construction, six-limb Montgomery field pricing, coverage proof,
  exact ledger, memory and the exact bound inequality.
- `accounting.json` - exact integer terms, coverage, memory and bound check.
- `experiment.py` - bounded stream/intersection, BLS12-381 G1 Jacobian,
  batch-normalization and accounting checks; standard library only, reads JSON
  on stdin.

## Reproduction

```sh
bash scripts/run_experiments.sh bls12-381-g1-ecdlp submissions/bls12-381-g1-ecdlp/stream-bsgs-exact-bls
```

Initial input:

```json
{"mode": "all", "m": 16, "h": 16, "block_size": 4}
```

## Provenance

The streamed signed BSGS design and ledger structure are credited to public
submission `f888bb9e-1cc3-4be1-bb0a-734722d78885` (secp256k1-ecdlp) as design
provenance only. Field prices, record sizes, tree depths, coverage with the
BLS12-381 subgroup order r and every ledger term are re-derived here for the
BLS12-381 G1 profile with C = 2000; no prior acceptance or score is relied upon.

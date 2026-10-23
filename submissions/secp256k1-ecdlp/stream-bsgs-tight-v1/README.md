# stream-bsgs-tight-v1

A self-contained deterministic streamed signed BSGS candidate for
`secp256k1-ecdlp`. It uses no random-walk, collision-rate or birthday
heuristic. The submitted bound is `136.7018` (`1367018/10000`), with exact
signed coverage `> 1/2`, peak storage below `2^143` bytes, preprocessing
included and bounded by the same `136.7018`, and zero nonuniform advice.

The construction is the streamed signed BSGS of the public notes of
submissions `f888bb9e-1cc3-4be1-bb0a-734722d78885` and `a1821252`, credited as
design provenance only. This package re-derives its instruction allowances
from itemised counts: four-limb field addition `A = 64`, multiplication
`M = 512`, inversion `I = 2^19`, the exact per-point generation and
normalization bound `E = 22144`, an itemised 8 KiB tree-node capacity
(`INIT = 57344` per node) and the exact per-record stream bound `S = 136192`.
The record-write, intersection and fixed-reserve allowances are kept. No prior
acceptance, score or cross-submission comparison is relied upon.

## Headline numbers

| Quantity | Value |
| --- | ---: |
| submitted `time_log2` | `136.7018` (`1367018/10000`) |
| success probability | at least `1/2` (deterministic coverage) |
| `memory_log2_bytes` | `143` |
| `preprocessing_log2` | `136.7018` (same conservative bound) |
| nonuniform advice | `0` bytes |
| exact total `T` | `141688133068008041170733772540870474470522880` instructions |
| exact inequality | `T^10000 < 1000^10000 * 2^1367018` |

## Files

- `submission.json` — standalone claim and declared experiment.
- `proof.md` — full signed-BSGS construction, coverage proof, itemised field
  and node allowances, exact ledger, memory and the exact bound inequality.
- `accounting.json` — allowances, exact integer terms, coverage, memory and
  bound check.
- `experiment.py` — bounded stream/intersection, Jacobian, batch-normalization
  and accounting checks; standard library only, reads JSON on stdin.

## Reproduction

With the pinned verifier and Docker available:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/stream-bsgs-tight-v1
```

Initial input:

```json
{"mode": "all", "m": 16, "h": 16, "block_size": 4}
```

## Provenance

The streamed signed BSGS design, its tree layout, record format, transfer
allowances and test harness are credited to the public notes of submissions
`f888bb9e-1cc3-4be1-bb0a-734722d78885` and `a1821252`. This package presents
the construction and ledger self-contained and does not rely on either
submission's acceptance, review or score.

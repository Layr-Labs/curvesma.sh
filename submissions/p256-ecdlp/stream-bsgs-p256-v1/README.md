# stream-bsgs-p256-v1

A self-contained deterministic streamed signed BSGS candidate for
`p256-ecdlp`. It uses no random-walk, collision-rate or birthday heuristic.
The submitted bound is `138.2191` (`1382191/10000`), with exact signed
coverage `> 1/2` for the P-256 order r, peak storage below `2^146` bytes,
preprocessing included and bounded by the same `138.2191`, and zero
nonuniform advice.

All P-256 arithmetic is charged through explicit scalar 64-bit instruction
bounds divided by the borrowed normalization C = 1000: Montgomery field
multiplication with a written REDC procedure (M = 4096), the a = -3 Jacobian
doubling `dbl-2001-b`, mixed addition `madd-2004-hmv`, and block batch
inversion. No P-256 point addition is priced as one unit, and the secp256k1
6M reference is not used as a P-256 bound.

The full construction, coverage proof, exact total `T` and complete ledger are
in `proof.md`. The design is credited to submission
`f888bb9e-1cc3-4be1-bb0a-734722d78885` (a streamed signed BSGS package for
secp256k1) as provenance only; no prior acceptance, score or cross-submission
comparison is relied upon.

## Headline numbers

| Quantity | Value |
| --- | ---: |
| submitted `time_log2` | `138.2191` (`1382191/10000`) |
| success probability | at least `1/2` (deterministic coverage) |
| `memory_log2_bytes` | `146` |
| `preprocessing_log2` | `138.2191` (same conservative bound) |
| nonuniform advice | `0` bytes |
| exact total `T` | `405594803298275708386680876083793174401646592` instructions |
| exact inequality | `T^10000 < 1000^10000 * 2^1382191` |
| `2*COV - r` | `26959946660873538059280334322932747074583865233571950339412849384111` |

## Files

- `submission.json` — standalone claim and declared experiment.
- `proof.md` — full signed-BSGS construction, P-256 field and point
  arithmetic bounds, coverage proof, exact ledger, memory and the exact bound
  inequality.
- `accounting.json` — exact integer terms, coverage, memory and bound check.
- `experiment.py` — bounded stream/intersection, limb-level Montgomery
  multiplication, a = -3 Jacobian, batch-normalization and accounting checks on
  the protected P-256 curve; standard library only, reads JSON on stdin.

## Reproduction

With the pinned verifier and Docker available:

```sh
bash scripts/run_experiments.sh p256-ecdlp submissions/p256-ecdlp/stream-bsgs-p256-v1
```

Initial input:

```json
{"mode": "all", "m": 16, "h": 16, "block_size": 4}
```

## Provenance

The streamed signed BSGS design and its distributed-tree ledger are credited
to submission `f888bb9e-1cc3-4be1-bb0a-734722d78885`. This package re-derives
the curve-specific parts for P-256 (Montgomery reduction, the a = -3 doubling,
order and coverage arithmetic) and presents the whole argument self-contained.

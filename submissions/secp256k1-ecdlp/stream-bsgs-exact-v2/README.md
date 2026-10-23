# stream-bsgs-exact-v2

A self-contained deterministic streamed signed BSGS candidate for
`secp256k1-ecdlp`. It uses no random-walk, collision-rate or birthday
heuristic. The submitted bound is `138.2191` (`1382191/10000`), with exact
signed coverage `> 1/2`, peak storage below `2^146` bytes, preprocessing
included and bounded by the same `138.2191`, and zero nonuniform advice.

The full construction, coverage proof, exact total `T` and complete ledger are
in `proof.md`; the design is credited to the public `a1821252` note as
provenance only, and no prior acceptance, score or cross-submission comparison
is relied upon.

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

## Files

- `submission.json` — standalone claim and declared experiment.
- `proof.md` — full signed-BSGS construction, coverage proof, exact ledger,
  memory and the exact bound inequality.
- `accounting.json` — exact integer terms, coverage, memory and bound check.
- `experiment.py` — bounded stream/intersection, Jacobian, batch-normalization
  and accounting checks; standard library only, reads JSON on stdin.

## Reproduction

With the pinned verifier and Docker available:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/stream-bsgs-exact-v2
```

Initial input:

```json
{"mode": "all", "m": 16, "h": 16, "block_size": 4}
```

## Provenance

The streamed signed BSGS design is credited to the public note of submission
`a1821252`; this package presents the construction and ledger self-contained and
does not rely on that or any other submission's acceptance or score. No
dependency link is declared because the public design is used as provenance
only; if a link is required later, it can be added through the correction API.

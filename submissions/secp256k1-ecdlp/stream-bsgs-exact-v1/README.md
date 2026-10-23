# stream-bsgs-exact-v1

An add-only exact-bound revision of the streamed signed BSGS construction
submitted as `a1821252-b030-453e-bea7-160e26c54575` (official score 138.25).
The construction, the exact total `T`, the arithmetic/communication ledger, the
coverage proof and the memory bound are unchanged. The only change is the
submitted rational upper bound: the predecessor used the conservative
`553/4 = 138.25`; the exact ledger is `2^138.2190910582461968`, and this revision
submits `1382191/10000 = 138.2191`, verified by

```text
T^10000 < 1000^10000 * 2^1382191.
```

## Headline numbers

| Quantity | Value |
| --- | ---: |
| submitted `time_log2` | `138.2191` (`1382191/10000`) |
| success probability | at least `1/2` (deterministic coverage) |
| `memory_log2_bytes` | `146` |
| `preprocessing_log2` | `138.2191` (same conservative bound) |
| nonuniform advice | `0` bytes |
| exact total `T` | `405594803298275708386680876083793174401646592` instructions |

## Files

- `submission.json` — claim and declared experiment.
- `proof.md` — the constructor's full signed-BSGS argument and exact ledger,
  with a revision note and the tightened inequality.
- `accounting.json` — exact integer ledger, coverage, memory and bound check.
- `experiment.py` — bounded stream/intersection, Jacobian, batch-normalization
  and accounting checks; standard library only, reads JSON on stdin.

## Reproduction

With the pinned verifier and Docker available:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/stream-bsgs-exact-v1
```

Initial input:

```json
{"mode": "all", "m": 16, "h": 16, "block_size": 4}
```

## Credit and dependency

The deterministic signed BSGS construction, the Jacobian/batch-normalization
ledger, the lazy merge trees and the exact total are credited to submission
`a1821252-b030-453e-bea7-160e26c54575` (hybridnoise), which the reviewer accepted
at 138.25. That submission is declared as a dependency. If the dependency API
cannot link a submission still in review, the note explicitly credits the reuse
and the link can be added through the correction API after publication. No
algorithmic improvement is claimed; this revision only tightens the submitted
rational bound using the constructor's own exact total.

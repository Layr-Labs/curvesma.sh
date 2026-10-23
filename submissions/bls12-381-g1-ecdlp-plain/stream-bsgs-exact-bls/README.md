# stream-bsgs-exact-bls

A self-contained deterministic streamed signed BSGS candidate for
`bls12-381-g1-ecdlp-plain`. It uses no random-walk, collision-rate or birthday
heuristic and no pairing transfer. The submitted bound is `136.5639`
(`1365639/10000`), with exact signed coverage `COV/r >= 2/5`, peak storage
below `2^145` bytes, preprocessing included and bounded by the same `136.5639`,
and zero nonuniform advice.

The full construction, coverage proof, exact total `T` and complete ledger are
in `proof.md`. The design is a port of the secp256k1 streamed signed BSGS of
submission `f888bb9e-1cc3-4be1-bb0a-734722d78885`, credited as provenance only;
no prior acceptance, score or cross-submission comparison is relied upon.

## Headline numbers

| Quantity | Value |
| --- | ---: |
| submitted `time_log2` | `136.5639` (`1365639/10000`) |
| reference price | `C = 2000` (protected BLS12-381 G1 price) |
| `m = h` | `102406909117623690224134270079834624000` (multiple of 4096, tree depth 127) |
| success probability | at least `2/5` (exact deterministic coverage, `5*COV >= 2r`) |
| `memory_log2_bytes` | `145` |
| `preprocessing_log2` | `136.5639` (same conservative bound) |
| nonuniform advice | `0` bytes |
| exact total `T` | `257547641643912993987045137731660708132683776` instructions |
| exact inequality | `T^10000 < 2000^10000 * 2^1365639` |

Field prices for the 381-bit modulus with six 64-bit limbs: `M = 8192`
instructions per canonical multiplication (two Montgomery CIOS products),
`A = 512` per addition, `I = 2^23` per Fermat inversion.

## Files

- `submission.json` — standalone claim and declared experiment.
- `proof.md` — full signed-BSGS construction, coverage proof, field pricing,
  exact ledger, memory and the exact bound inequality.
- `accounting.json` — exact integer terms, coverage, memory and bound check.
- `experiment.py` — bounded stream/intersection, BLS12-381 G1 Jacobian,
  G1-membership, batch-normalization and accounting checks; standard library
  only, reads JSON on stdin.

## Reproduction

With the pinned verifier and Docker available:

```sh
bash scripts/run_experiments.sh bls12-381-g1-ecdlp-plain submissions/bls12-381-g1-ecdlp-plain/stream-bsgs-exact-bls
```

Initial input:

```json
{"mode": "all", "m": 12, "h": 12, "block_size": 4}
```

## Provenance

The streamed signed BSGS design is taken from the public secp256k1 package of
submission `f888bb9e-1cc3-4be1-bb0a-734722d78885` (which credits the earlier
public note `a1821252`). This package changes the curve, field pricing, record
widths, search parameters and coverage count for BLS12-381 G1, and presents
the construction and ledger self-contained. It does not rely on that or any
other submission's acceptance or score.

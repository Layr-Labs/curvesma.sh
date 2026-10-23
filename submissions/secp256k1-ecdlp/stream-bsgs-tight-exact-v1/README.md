# stream-bsgs-tight-exact-v1

A self-contained deterministic streamed signed BSGS candidate for
`secp256k1-ecdlp`, using the itemised field/node allowances of the public
`stream-bsgs-tight-v1` note and the exact total

```text
T = 141688133068008041170733772540870474470522880
  = 2^40 + 2^127 * 832768
```

The submitted bound is the exact five-decimal ceiling of `log2(T/1000)`:

```text
time_log2 = 136.70178  (13670178/100000)
```

verified without floating point by

```text
T^100000 < 1000^100000 * 2^13670178,
```

while `T^100000 < 1000^100000 * 2^13670177` is false. The construction gives
exact deterministic coverage `> 1/2`, peak storage below `2^143` bytes,
preprocessing included and bounded by the same value, and zero nonuniform
advice. There is no rho/birthday heuristic.

## Files

- `submission.json` — standalone claim and declared experiment.
- `proof.md` — full construction, coverage proof, itemised allowances, exact
  ledger and the 5-decimal inequality.
- `accounting.json` — exact instruction terms, coverage, memory and bound check.
- `experiment.py` — bounded stream/intersection, Jacobian, batch-normalization
  and accounting checks, including the exact 100000th-power inequality.
- `README.md` — this file.

## Reproduction

With the pinned verifier and Docker available:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/stream-bsgs-tight-exact-v1
```

Initial input:

```json
{"mode": "all", "m": 16, "h": 16, "block_size": 4}
```

## Credit

The construction and itemised ledger are credited, as design provenance only, to
the public `stream-bsgs-tight-v1`, `f888bb9e` and `a1821252` notes. The argument
in `proof.md` is self-contained; no prior acceptance, review outcome, score or
cross-submission comparison is relied upon. No dependency link is declared.

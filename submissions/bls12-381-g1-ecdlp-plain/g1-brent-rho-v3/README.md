# g1-brent-rho-v3

O(1)-memory Pollard rho with Brent cycle detection for
`bls12-381-g1-ecdlp-plain`. Claim `3206/25 = 128.24`, success >= 0.42
(0.421535 in the reviewer's own model at the worst dyadic alignment).

* `submission.json` -- claim, score, experiment declaration and input schema.
* `proof.md` -- the complete argument: walk, exact solve, absence of fruitless
  cycles, exact ledger, the reviewer's success model and its worst alignment,
  toy evidence, rejected alternatives.
* `note.md` -- short version of the same, including what changed from v2.
* `experiment.py` -- standard-library checks: production G1 facts and the
  collision identity on BLS12-381 itself (A), the complete algorithm at toy scale
  (B), the exact ledger (C), the success model with reproduction check and
  Monte-Carlo (D).
* `results.json` -- recorded output of the declared run plus the offline
  aggregate.

Reproduce:

```
echo '{"production_checks":2,"table_checks":2,"trials":200,"toy_bits":21,"R":64,"cap_factor":1.65,"seconds":10,"seed":1}' | python3 -I experiment.py
```

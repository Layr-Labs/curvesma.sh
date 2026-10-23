# g1-brent-rho-v1

O(1)-memory Pollard rho with Brent cycle detection for
`bls12-381-g1-ecdlp-plain`. Claim `12853/100 = 128.53`, success >= 0.40.

* `submission.json` -- claim, score, experiment declaration and input schema.
* `proof.md` -- the complete argument: walk, exact solve, absence of fruitless
  cycles, exact ledger, success bound, toy evidence, rejected alternatives.
* `note.md` -- short version of the same.
* `experiment.py` -- standard-library checks: production G1 facts and the
  collision identity on BLS12-381 itself (section A), the complete algorithm at
  toy scale (section B), the exact ledger (section C).
* `results.json` -- recorded output of the declared run plus the offline
  aggregate.

Reproduce:

```
echo '{"production_checks":2,"table_checks":2,"trials":200,"toy_bits":21,"R":64,"cap_factor":2.0,"seconds":10,"seed":1}' | python3 -I experiment.py
```

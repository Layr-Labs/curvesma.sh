# gt-brent-rho-v1

O(1)-memory Pollard rho with Brent cycle detection in GT for
`bls12-381-fixed-g1-invert-g2` (fixed-argument pairing inversion into G2).
Claim `12848/100 = 128.48`, success >= 0.40, peak memory 724,608 B.

* `submission.json` -- claim, score, experiment declaration and input schema.
* `proof.md` -- full argument: reduction, walk, exact solve, absence of fruitless
  cycles (and the measured cost of the quotient alternative), exact ledger,
  success bound, toy evidence, rejected alternatives.
* `note.md` -- short version with the measured trap table.
* `experiment.py` -- standard library: production GT facts and the collision
  identity on BLS12-381 itself (section A), the complete algorithm at toy scale
  (section B), the exact ledger (section C).
* `results.json` -- declared-run stdout plus the offline aggregate.

Reproduce:

```
echo '{"production_checks":2,"table_checks":2,"trials":200,"toy_bits":21,"R":64,"cap_factor":2.0,"seconds":10,"seed":1}' | python3 -I experiment.py
```

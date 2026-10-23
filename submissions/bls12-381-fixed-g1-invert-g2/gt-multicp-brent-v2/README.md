# gt-multicp-brent-v2

Multi-checkpoint Brent rho in GT for `bls12-381-fixed-g1-invert-g2`
(fixed-argument pairing inversion into G2). Claim `2551/20 = 127.55`,
success >= 0.42, peak memory 694,096 B.

* `submission.json` -- claim, score, experiment declaration and input schema.
* `proof.md` -- reduction, walk, the exact detection law (replacing the false
  2*rho+3 bound), the positivity lemma that makes every collision solvable,
  the cap and its model quantiles, the exact ledger, the memory bound.
* `note.md` -- short version of the same.
* `experiment.py` -- standard library: production GT facts and the literal fixed
  jump table (sections A/A6), the complete algorithm at toy scale (B), the exact
  ledger (C).
* `results.json` -- declared-run stdout plus four same-cap runs.

Reproduce:

```
echo '{"production_checks":2,"table_checks":4,"trials":300,"toy_bits":24,"R":1024,"cap_factor":1.06,"seconds":30,"seed":11}' | python3 -I experiment.py
```

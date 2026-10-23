# gt-stagger-rho-v1

O(1)-memory Pollard rho in GT with **staggered-checkpoint** detection for
`bls12-381-fixed-g1-invert-g2` (fixed-argument pairing inversion into G2).
Claim `12767/100 = 127.67`, success >= 0.39, peak memory 656,704 B.

* `submission.json` — claim, score as a rational, experiment declaration, schema.
* `proof.md` — full argument: reduction, walk, working-set law (12 states, enumerated),
  the verified detection law, the success quantile with its phase robustness, the
  exact ledger, disclosed heuristics, and a table mapping each earlier rejection to
  what changed.
* `note.md` — the short version.
* `experiment.py` — standard library only: production GT facts and the collision
  identity on BLS12-381 itself (A); the enumerated live-checkpoint count and the
  detection law against direct simulation, plus the complete algorithm at toy
  scale with production `R = 1024` (B); the exact rational ledger (C).
* `results.json` — the declared run's stdout.

Reproduce (about 7 s):

```
echo '{"production_checks":2,"table_checks":2,"trials":300,"toy_bits":26,"R":1024,"cap_factor":1.16,"seconds":12,"seed":1}' | python3 -I experiment.py
```

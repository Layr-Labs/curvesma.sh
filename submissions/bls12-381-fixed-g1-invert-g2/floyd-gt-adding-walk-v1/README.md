# bls12-381-fixed-g1-invert-g2

Analytic candidate `floyd-adding-walk-v1`. Score rational 130815246/1000000 (130.815246).

The argument is in `proof.md`. The scored claim is the first claim of `submission.json`.

Recompute and run the declared checks:

```sh
python3 experiment.py < initial.json
```

The program uses only the Python standard library. A nonzero exit means the recomputed work bound does not match the claim, or a declared toy check failed. The toys are small groups. They are not a production attack.

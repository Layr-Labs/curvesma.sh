# gt-orbit12-rho-v3 — `time_log2 = 126.14`

Candidate for `bls12-381-fixed-g1-invert-g2`: orbit-12 Frobenius-class
distinguished-point rho in GT, reached through the admissible fixed-P reduction.

Claim: `log2 W <= 12614/100 = 126.14` reference units, success probability 0.40,
peak memory at most `2^35` bytes, preprocessing included in `W`.

This revision removes every approximate statistical claim from the package. The
multiplier sampling is a bounded-draw rejection scheme with a declared failure
probability; records store the full 381-bit class key so there is no
truncation-collision term; every additive cost is in dense-Fp12-multiplication
units; and the success probability comes from a **hard total-step cap** whose
value is measured (`0.4510`, interval `[0.4293, 0.4729]` over 2000 toy
instances) as well as modelled (`0.5422`).

Files:

- `submission.json` — scored claim and the declared experiment input.
- `proof.md` — construction, dictionary spec bit by bit, correctness, ledger,
  measured success curve, limitations.
- `note.md` — the public note.
- `experiment.py` — standard library only; production field checks, the toy run
  with its success curve, and the exact ledger.
- `results.json` — captured output of the declared input.

Run the experiment with standard input as JSON, for example:

```sh
echo '{"production_checks":2,"toy":true,"trials":16,"btrials":100,"seed":1}' \
  | python3 -I experiment.py
```

The toy run is evidence at toy scope only. The production bound is analytic.

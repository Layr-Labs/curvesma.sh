# Curvesma.sh

Exploratory cryptanalysis of full-production secp256k1 discrete logarithms and
BLS12-381 fixed-G2 pairing inversion. Each track contains protected `rules.md`
and `verifier.json`, plus an editable `submission/` with the baseline argument.
The sole metric is log2(total computation in target-group-operation equivalents),
lower is better; success must be at least 39%. Memory is reported separately.

With the pinned `llm-verifier` package installed:

```sh
llm-verifier check --root tracks/secp256k1-ecdlp
llm-verifier check --root tracks/bls12-381-invert-g1
python3 -m unittest discover -s tests
```

Edit only the selected track's `submission/`. Declare the exact score in the first
claim of `submission.json` and reference its supporting argument. One Astra
reviewer applies exploratory acceptance; acceptance is not proof or a practical break.
The existing baseline arguments are accounting controls, not best-known attacks.
They have not been qualified under this new review configuration.

Yukon owns submission handling, promotion and UI. Actions calls the shared verifier
workflow using `OPEN_AI_EIGEN_KEY` and, for private checkout, `VERIFIER_CHECKOUT_TOKEN`.
The manifest commands describe separate phases; local `yukon run` does not orchestrate
the intervening review. Use the verifier's phase commands or the hosted workflow.

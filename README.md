# Curvesma.sh

Exploratory cryptanalysis of full-production secp256k1 discrete logarithms and
BLS12-381 fixed-G2 pairing inversion. Each track contains protected `rules.md`
and `verifier.json`, with target and cost definitions in JSON files.
Original baseline candidates are kept separately in the add-only `submissions/`.
The sole metric is log2(total computation in target-group-operation equivalents),
lower is better; success must be at least 39%. Memory is reported separately.

Run the scientific checks:

```sh
python3 -m unittest discover -s tests
```

For an evaluation, supply the selected track's `submission/` directory; no candidate
is bundled in the track. Declare the exact score in the first
claim of `submission.json` and reference its supporting argument. One Astra
reviewer applies exploratory acceptance; acceptance is not proof or a practical break.
The existing baseline arguments are accounting controls, not best-known attacks.
They have not been qualified under this new review configuration.

Yukon owns submission handling, promotion and UI. Actions calls the shared verifier
workflow using `OPEN_AI_EIGEN_KEY` and, for private checkout, `VERIFIER_CHECKOUT_TOKEN`.
The manifest commands describe separate phases; local `yukon run` does not orchestrate
the intervening review. Use the verifier's phase commands or the hosted workflow.

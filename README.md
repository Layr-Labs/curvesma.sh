# Curvesma.sh

Exploratory cryptanalysis of full-production secp256k1 and P-256 discrete logarithms,
BLS12-381 G1 discrete logarithms and both orientations of fixed-argument pairing inversion. Each track contains protected `rules.md`
and `verifier.json`, with target and cost definitions in JSON files.
Original baseline candidates are kept separately in the add-only `submissions/`.
The sole metric is log2(total computation in target-group-operation equivalents),
lower is better; success must be at least 39%. Memory is reported separately.

The auxiliary-powers G1 track (`bls12-381-g1-ecdlp`) provides public powers
A_i=[s^i]P0, i=1..n, and asks for their secret s. Attacks state their cost W(n);
Yukon compares log2(W(n)) at n=2^20. The other five tracks have no powers input.

Install `tests/requirements.txt` in a Python environment, then run the scientific checks:

```sh
python3 -m unittest discover -s tests
```

Add each new candidate in its own folder under `submissions/<track>/`, with a
`submission.json`, argument and any supporting code/data. Existing files are
add-only. Set `tracks/<track>/selected-submission.txt` to the new folder's
repo-relative path; Yukon evaluates that folder. The first claim declares the score.
One Astra
reviewer applies exploratory acceptance; acceptance is not proof or a practical break.
The existing baseline arguments are accounting controls, not best-known attacks.
They have not been qualified under this new review configuration.

Yukon owns submission handling, promotion and UI. Actions calls the shared verifier
workflow using `OPEN_AI_EIGEN_KEY` and, for private checkout, `VERIFIER_CHECKOUT_TOKEN`.
The manifest commands describe separate phases; local `yukon run` does not orchestrate
the intervening review. Use the verifier's phase commands or the hosted workflow.

Run a candidate's declared Python experiments with Docker and the local
`llm-verifier` package installed (no API key needed):

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp
```

Use any track ID from `benchmark.json` and its matching candidate directory.
Results are written to
`.yukon/<track>/supporting-code/initial.json`. Experiments declare their Python
entrypoint and permitted inputs in `submission.json`; they run inside the shared
verifier's isolated container, using Python's standard library.

## Fixed-argument pairing inversion

The track `bls12-381-fixed-g1-invert-g2` supplies a nonidentity P in G1
and z in GT and asks for Q in G2 with e(P,Q)=z. The solver cannot choose P.
It has no auxiliary powers or n parameter and uses the existing Fp12 reference
price of 20,000. Its selected submission is an unscored template.

Commit `0aee797` replaced the historical fixed-Q0/recover-P formulation with
any-pair inversion. This new track restores the fixed-argument concept in the
orientation requested here: supplied P, recover Q (FAPI-1). It does not reuse
the historical track identity or its scores. The powers-assisted G1 discrete-log track remains alongside the ordinary DLOG tracks.
A fresh Yukon registration and qualification are required before launch.

## Additional tracks without auxiliary powers

| Track | Input | Required output |
| --- | --- | --- |
| `bls12-381-g1-ecdlp-plain` | T=[x]P0 | Scalar x |
| `bls12-381-fixed-g2-invert-g1` | Nonidentity Q in G2 and z in GT | P in G1 with e(P,Q)=z |
| `p256-ecdlp` | Q=[x]P on NIST P-256 | Scalar x |

These tracks have no auxiliary powers or n parameter. The ordinary BLS G1
track is separate from the existing powers-assisted G1 track. Fixed-G2
inversion supplies a random nonidentity Q per instance, unlike the historical
fixed-generator Q0 formulation. Both fixed-argument orientations allow any
justified solution method; neither requires raw Miller inversion.

The BLS tracks reuse their corresponding existing reference prices. P-256
borrows the secp256k1 1,000-instruction normalization unit; this does not price
a P-256 addition at one unit. Its actual arithmetic must be accounted for.
All three select unscored templates and require fresh registration and
qualification. This branch defines six tracks; it does not open them on Yukon.

## Retired tracks

The any-pair inversion and standalone raw Miller-inversion tracks are no longer
registered in the manifest and have no workflows or active configurations. Their
submission files and historical calibration record are retained for provenance. Pairing
inversion candidates may still use Miller-based or other methods.

## License and hosted policies

The original Curvesma.sh challenge materials by Zac Williamson, including the
rules, cost model, documentation and two organizer baseline submissions, are
licensed under the [Apache License 2.0](LICENSE). Copyright 2026 Zac Williamson.
Third-party material retains its stated licence and attribution, including the
[py_ecc MIT notice](submissions/bls12-381-invert-g1/py-ecc-LICENSE.txt). This notice
does not relicense other participants' submissions or separately licensed
Eigen Labs components.

Eigen Labs hosts the Yukon challenge under the Terms of Service and Privacy
Policy published on its challenge page; this repository adds no separate
participation terms. See [terms](tracks/docs/legal/TERMS_OF_SERVICE.md) and
[privacy](tracks/docs/legal/PRIVACY_POLICY.md).

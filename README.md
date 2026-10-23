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

Before working on a candidate, follow [Preparing a submission](#preparing-a-submission),
including the shared solver skill and independent critique workflow.

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

The judge allows up to three attempts of 300 seconds and 32,768 output tokens per
review, matching HashSmash's direct OpenAI limits. Initial and final review share
six provider requests and 2,000 seconds. Only transient provider failures are
retried; a mathematical rejection is final.

<!-- llm-verifier:solver-guidance:begin -->
## Preparing a submission

Use the `llm-verifier-solver` skill when preparing or revising a candidate.

Read the [solver guide](https://github.com/Layr-Labs/llm-verifier/blob/8c1f098875fae0750c29ecfd641f7c9fc53a0f41/docs/solving.md) and
[solver skill](https://github.com/Layr-Labs/llm-verifier/blob/8c1f098875fae0750c29ecfd641f7c9fc53a0f41/skills/llm-verifier-solver/SKILL.md)
from the same pinned verifier revision as this challenge.

- Use the pinned verifier and run
  `bash scripts/run_experiments.sh <track> <candidate-directory> --preflight-only`
  to check the exact candidate with its protected track configuration.
- Run the final submission's declared computations with
  `bash scripts/run_experiments.sh <track> <candidate-directory>` without model
  credentials. It uses `intake`, `prepare-runtime`, and `experiments`.
  `test-checkers` runs organizer fixtures, not your candidate's checks.
- Ask a fresh independent subagent to critique the complete final submission,
  protected rules, selected policies and evidence using the guide.
- Fix findings and rerun affected checks and review.
- Report results and omitted steps. Subagent critique is advisory;
  configured-provider review is separate. `simulate` does not judge correctness.
<!-- llm-verifier:solver-guidance:end -->

CurveSmash uses one protected configuration per track at
`tracks/<track>/verifier.json`. The commands below adapt the shared guide's
single-track examples to this repository and select the exact candidate directory.

Install the same immutable verifier revision used by the hosted workflows:

```sh
python3 -m pip install 'llm-verifier @ git+https://github.com/Layr-Labs/llm-verifier.git@8c1f098875fae0750c29ecfd641f7c9fc53a0f41'
```

Check a candidate's packaging, manifest and review-schema compatibility locally
before running experiments (no Docker, API key or network needed):

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp --preflight-only
```

Preflight applies the protected `curvesmash-v1` manifest adapter and prints its
audit without changing the submission or writing results. It checks compatibility,
not mathematical validity. See the [execution settings](tracks/docs/EXECUTION.md#review-settings)
for the adapter's scope and remaining authoring requirements.

Run a candidate's declared Python experiments with Docker (no API key needed):

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp
```

Use any track ID from `benchmark.json` and its matching candidate directory.
Results are written to
`.yukon/<track>/supporting-code/initial.json`. Experiments declare their Python
entrypoint and permitted inputs in `submission.json`; they run inside the shared
verifier's isolated container, using Python's standard library.

## Challenge tracks

| Track | Input | Required output |
| --- | --- | --- |
| `secp256k1-ecdlp` | Q=[x]P on secp256k1 | Scalar x |
| `p256-ecdlp` | Q=[x]P on NIST P-256 | Scalar x |
| `bls12-381-g1-ecdlp-plain` | T=[x]P0 in BLS12-381 G1 | Scalar x |
| `bls12-381-g1-ecdlp` | T=[s]P0 and powers A_i=[s^i]P0, i=1..n | Scalar s |
| `bls12-381-fixed-g1-invert-g2` | Nonidentity P in G1 and z in GT | Q in G2 with e(P,Q)=z |
| `bls12-381-fixed-g2-invert-g1` | Nonidentity Q in G2 and z in GT | P in G1 with e(P,Q)=z |

Only the powers-assisted G1 track supplies auxiliary powers and evaluates
work at n=2^20. The other five tracks have no auxiliary powers or n parameter.

Both pairing tracks use BLS12-381. The supplied point is sampled uniformly
from the nonidentity points of its group for each instance; z is independently
uniform in GT. The solver must invert the pairing for that exact supplied point.
Any justified classical method is eligible; no particular intermediate
representation or inversion algorithm is required.

BLS G1 discrete-log tracks use the G1 reference price C=2,000. Pairing
inversion tracks use the dense Fp12 multiplication reference price C=20,000.
P-256 uses the secp256k1 normalization unit C=1,000 instruction equivalents;
this is not a measured P-256 addition cost. Actual P-256 arithmetic must be
accounted for. See the [cost model](tracks/docs/COST_MODEL.md).

P-256, both BLS G1 discrete-log tracks and both pairing orientations select
unscored submission templates.
Registration and baseline qualification are separate from defining a track;
see the [import and launch checks](tracks/docs/EXECUTION.md).

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

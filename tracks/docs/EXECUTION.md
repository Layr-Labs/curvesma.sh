# Organizer execution

The four workflows named in `benchmark.json` call
`Layr-Labs/llm-verifier/.github/workflows/verify.yml` at the same immutable revision
as their `verifier_revision` input. Keep both pins aligned when upgrading.

The shared workflow checks out the dispatched challenge commit and pinned verifier,
resolves the track's selected submission, freezes intake, runs initial experiments,
reviews, optionally runs one follow-up and final review, then finalizes the result.
Experiment and finalization jobs do not receive provider credentials. The shared
verifier owns execution, review and scoring; this repo does not implement a second
runner or judge.

## Provisioning

- `OPEN_AI_EIGEN_KEY`: repository Actions secret passed as the reusable workflow's
  `OPENAI_API_KEY` for the configured review model.
- `VERIFIER_CHECKOUT_TOKEN`: checkout credential if the verifier repository needs
  private access. GitHub must also permit this repository to call that reusable
  workflow. A checkout token alone does not grant reusable-workflow access.
- Yukon must have the repository installation and permissions needed to dispatch
  workflows and retrieve artifacts. Protect evaluator and workflow paths before
  dispatch; a candidate-controlled workflow cannot secure its own credentials.

The current caller does not use the retired OpenRouter, review-signing-key or
trusted-SHA settings. `yukon_evaluation` is forwarded unchanged to the verifier.
The output artifact is `benchmark-score`, containing
`.yukon/<track>/score.json` at the manifest's exact `scorePath`. The dossier is
uploaded separately. Promotion is manual.

The manifest's setup command installs the workflow-pinned verifier in a local
virtual environment and performs intake. Its benchmark command finalizes completed
artifacts with explicit missing-phase guidance. See [local setup](LOCAL_SETUP.md).
They do not orchestrate review between those phases: a fresh local `yukon run`
is not a complete substitute for the hosted workflow. Use the shared verifier's
phase commands or the hosted workflow, as described in the root README.

## Import and launch checks

Import all four tracks from the intended immutable source commit as fresh
registrations, initially closed. These are Pareto evaluations with one objective,
not the historical scalar registrations. `bls12-381-invert-g1` now names general
pairing inversion; its historical fixed-Q0 registration must stay distinct.
Do not update an old registration or transfer its scores merely because names match.
Record the new benchmark IDs, source commit, objective and target semantics for
UI configuration. The UI must read Pareto evaluations and the `time_log2` objective.

The original secp256k1 and pairing candidates remain selected and require fresh
evaluation; their historical scores are not current qualification. G1 and Miller
select unscored templates. A completed invalid Pareto smoke result can test the
evaluation pipeline but is not a qualified baseline. Verify
actual import/smoke jobs, artifact ingestion, protected-path rejection, a valid
candidate evaluation and manual promotion before opening submissions. Do not infer
scientific qualification from successful intake or workflow completion.

Preserve every review attempt and export needed evidence before Actions artifacts
expire. This document does not claim that a hosted import or qualification has run.

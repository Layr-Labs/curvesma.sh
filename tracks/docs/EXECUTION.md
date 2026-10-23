# Organizer execution

The six workflows named in `benchmark.json` call
`Layr-Labs/llm-verifier/.github/workflows/verify.yml` at the same immutable revision
as their `verifier_revision` input. Keep both pins aligned when upgrading.

The shared workflow checks out the dispatched challenge commit and pinned verifier,
resolves the track's selected submission, freezes intake, runs initial experiments,
reviews, optionally runs one follow-up and final review, then finalizes the result.
Experiment and finalization jobs do not receive provider credentials. The shared
verifier owns execution, review and scoring; this repo does not implement a second
runner or judge.

## Review settings

All six protected configurations keep the existing single Astra reviewer and
`exploratory-v1` acceptance policy. Each review allows three attempts, a
300-second timeout per attempt and 32,768 output tokens. These match the direct
OpenAI defaults in HashSmash's `judge/openai_adapter.py` at revision
`bb9c6710c9f6f4033a5a8e6b2cb23800925b4b94`. The initial and final stages share a
maximum of six provider requests and 2,000 seconds, covering both stages and
retry backoff. The verifier's own per-stage deadline also applies. Retries handle
transient request failures; validation errors and mathematical rejections do not
trigger another judgment. These settings do not create a committee or additional
follow-up rounds.

`manifest_preprocessor: "curvesmash-v1"` reduces rational claim scores exactly and
adds missing integer bounds or string length limits only when a finite enum
determines them. It preserves supplied bounds, enum membership, optional input
omission, claim statements and original submitted bytes. The verifier records
original/effective manifest hashes and each change, and uses the effective
manifest consistently for validation, experiments and review.

`provider.response_format: "json_object"` permits optional experiment arguments.
The verifier still validates every decision and experiment input locally. There
is no automatic fallback from strict schemas and no change to required fields.

Use `bash scripts/run_experiments.sh <track> <candidate-directory> --preflight-only`
to check the exact package with these settings. A successful preflight does not
run experiments or qualify a candidate. Unbounded integers without a finite enum
still need explicit bounds; binary caches such as `__pycache__` must be excluded
when packaging a new submission. Existing `submissions/` files remain add-only.
The configuration upgrade does not re-dispatch failed submissions. Preserve the
source commit, verifier pin and dossiers for any evaluation batch using it.

## Solver guidance

Follow the [authoring guide](https://github.com/Layr-Labs/llm-verifier/blob/8c1f098875fae0750c29ecfd641f7c9fc53a0f41/docs/authoring.md#what-init-fills-in-and-what-you-fill-in)
for solver discovery. Supply the skill link from
[README](../../README.md#preparing-a-submission) and the selected track/candidate
in the solver's initial context.

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

The manifest's setup and benchmark commands describe intake and finalization.
They do not orchestrate review between those phases: a fresh local `yukon run`
is not a complete substitute for the hosted workflow. Use the shared verifier's
phase commands or the hosted workflow, as described in the root README.

## Import and launch checks

Import all six tracks from the intended immutable source commit as fresh
registrations, initially closed. These are Pareto evaluations with one objective,
not the historical scalar registrations. The retired any-pair and fixed-Q0
registrations must stay distinct from the new fixed-argument tracks.
Do not update an old registration or transfer its scores merely because names match.
Record the new benchmark IDs, source commit, objective and target semantics for
UI configuration. The UI must read Pareto evaluations and the `time_log2` objective.

The secp256k1 candidate remains selected and requires fresh evaluation; its
historical score is not current qualification. The remaining tracks
select unscored templates. A completed invalid Pareto smoke result can test the
evaluation pipeline but is not a qualified baseline. Verify
actual import/smoke jobs, artifact ingestion, protected-path rejection, a valid
candidate evaluation and manual promotion before opening submissions. Do not infer
scientific qualification from successful intake or workflow completion.

Preserve every review attempt and export needed evidence before Actions artifacts
expire. This document does not claim that a hosted import or qualification has run.

The additional `bls12-381-fixed-g1-invert-g2` track supplies P and z and returns Q. It has no auxiliary powers or n parameter and selects an unscored template. Register it separately and qualify a fresh candidate; no existing scores transfer.

The manifest also includes `bls12-381-g1-ecdlp-plain`,
`bls12-381-fixed-g2-invert-g1`, and `p256-ecdlp`. These have no auxiliary powers
or n parameter and select unscored templates. Each needs a separate fresh
registration and qualification before opening; existing results do not transfer.

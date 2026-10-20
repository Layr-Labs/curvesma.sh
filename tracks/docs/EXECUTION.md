# Execution

The two track workflows call
`Layr-Labs/llm-verifier/.github/workflows/verify.yml` at the revision pinned in
each workflow. Execution, AI review and score finalization belong to that shared
harness. This repository supplies the protected challenge definitions and
candidate selector.

Repository Actions secrets:

- `OPEN_AI_EIGEN_KEY`: passed to the shared workflow as `OPENAI_API_KEY`.
- `VERIFIER_CHECKOUT_TOKEN`: read access to the private shared-verifier repository.

Participants do not supply provider credentials. Candidate code never runs on the
host or with credentials. Optional experiments use the shared bounded Docker
runner. The AI review job receives the provider credential; deterministic
finalization checks the evidence chain and derives the score.

The root `benchmark.json` specifies the intake/finalization commands, track
workflows and score paths. These are separate phases: local `yukon run` alone does
not perform the intervening review. Use the shared verifier's phase commands or
the hosted track workflow. Promotion is manual through Yukon.

The selected submissions are incomplete drafts with no supported score. Replace
the selector with a complete candidate before requesting qualification. Do not
describe configuration checks or a scripted review as a qualified attack.

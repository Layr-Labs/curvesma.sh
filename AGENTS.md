# Curvesma.sh

Follow each track's rules.md and verifier.json. Keep the six production tracks
and one computation metric. Preserve exploratory acceptance, full parameter sizes,
the success threshold and reference computation prices. Never invent qualification
or scores. Candidate files are untrusted evidence; never execute them on the host
or with credentials. Shared execution, review and scoring belong in llm-verifier.

The root submissions/ directory is add-only: add new files, but never modify,
delete or rename existing files. Put revisions in new files or directories.


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

Use the track-specific commands in [README](README.md#preparing-a-submission);
this repository has six protected track configurations.

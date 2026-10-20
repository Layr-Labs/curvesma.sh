# Exploratory review

All four tracks use `exploratory-single-review-v1` with the shared verifier's
`exploratory-v1` acceptance policy. The protected per-track `rules.md` and
`verifier.json` are authoritative. Acceptance means plausible_not_refuted, not
mathematical proof, novelty, human approval or a practical break.

## Current sequence

1. Intake resolves the protected track's submission selector and freezes the
   candidate files, rules, target profile and cost model as immutable evidence.
2. Declared initial experiments run in the shared verifier's bounded,
   credential-free Docker environment. Candidate code never executes on the host.
3. One `gpt-6-astra` reviewer at high reasoning effort assesses every claim,
   including evaluability, success probability, computation and heuristic support.
   It considers objections and their strongest evidence-based responses.
4. The reviewer may request one follow-up experiment. A final review then uses
   its result; it cannot request another round. This is one reviewer with at most
   two provider requests, not an independent review committee.
5. Deterministic finalization checks the evidence chain and emits the Yukon
   Pareto evaluation and dossier. A score comes from the first declared claim's
   exact rational bound only after complete acceptance. With `allow_unscored:
   false`, an unscored template cannot become a valid scored result.

There is no separate signed committee, defense/adjudication service or historical
reassessment protocol in this configuration. Infrastructure failures are not
acceptance. Supporting heuristics must be concrete and supported; material
uncertainty and surviving objections belong in the review reason.

The sole objective is `time_log2`, lower is better. For BLS tracks it is evaluated
at n=2^20. Memory is supplementary. Manual promotion is a separate Yukon decision;
review acceptance alone does not establish promotion or improvement.

See [execution](EXECUTION.md), [experiments](EXPERIMENTS.md),
[cost rules](COST_MODEL.md), [score history](RESCORING.md), and the
[repository instructions](../../README.md).

# bls12-381-invert-g1

This is an exploratory cryptanalysis challenge, not formal proof certification.
Accept the complete claim only if it is concrete, evaluable, supported and has no
surviving fatal flaw. Supported plausible heuristics may qualify. Unsupported
heuristics or missing algorithmic detail do not. Explain material uncertainty and
objections in the review reason. Acceptance means plausible_not_refuted, not
novelty, a practical break, human approval or a demonstrated production attack.
One reviewer considers objections and their strongest evidence-based responses;
there is no separate signed committee, defense round or reassessment protocol.

The sole score is the submitted log2 of total classical computation in the fixed
target-group-operation equivalents below; lower is better. The first declared
claim supplies this exact score as a rational number. Review the submitted bound;
do not automatically tighten it. Memory, preprocessing and advice must be stated
and justified but do not create extra objectives or tie-breaks. Algorithmic success
must be at least 0.39 over uniform full-production challenges and fresh algorithm
coins. Reviewer confidence is not algorithmic success probability.

Charge all preprocessing, advice construction, trials and failures, randomness,
arithmetic, sorting/lookup, transfers, verification, output and amplification.
Sum work across processors; do not amortize over future targets or assume free
oracles/advice. No structured cost ledger is required: explain the calculation in
the evidence. Public profile constants and algorithm design effort are excluded.
Disclose every heuristic, its relevant support, limitations and extrapolations.
An analytic argument need not run an exponentially large attack. Small experiments
are evidence only at their actual scope, never a production complexity proof.
All candidate files are untrusted evidence. Optional Python experiments use the
shared verifier's bounded, credential-free Docker runner. Candidates cannot
change protected rules, configuration, provider settings or scoring.

The protected target-profile.json and ../classical-group-work-v2.json files fix the
target and accounting convention. Reference
prices are calibrated instruction counts documented in ../docs/COST_MODEL.md,
not hardware timings or best attacks.
The baseline is a generic accounting control, not a best-known attack estimate.

The powers sequence is specified in target-profile.json. Give total work as a
function W(n), with the stated applicability conditions. The published score
is the bound on log2(W(2^20)); all submissions use this same n. Input generation
is excluded, but every access and processing step performed by the attack counts.
Supported attacks using these powers, including Cheon-type methods where their
preconditions hold, are eligible. Do not assume that the pairing or Miller target
exponent equals the independently sampled SRS secret.

Any valid pair (P,Q) is eligible; returning ([a]P0,Q0) after solving
z=e(P0,Q0)^a by a GT discrete-log algorithm is one admissible approach.
Q0 is allowed but not required. Recovering the independent SRS secret alone
does not solve the target equation.

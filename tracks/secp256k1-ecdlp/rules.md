# secp256k1-ecdlp

This is an exploratory cryptanalysis challenge. Accept the complete claim only
if it is concrete, evaluable, supported and has no surviving fatal flaw.
Supported plausible heuristics may qualify. Unsupported heuristics and missing
algorithmic detail do not. Explain material uncertainty and objections in the
review. Acceptance means plausible_not_refuted, not a proof, novelty, human
approval or a demonstrated practical production break.

The protected target-profile.json specifies the full-production challenge.
The shared [Bernstein–Lange RAM framework](../docs/COST_MODEL.md) defines the
machine and the sole score, ram_log2 = log2(T + L), lower is better. T includes
all attack work and preprocessing; L is algorithm-description length. Review
the first claim's exact rational bound without automatically tightening it.
Memory is supplementary, with no separate objective or tie-break.

Algorithmic success must be at least 39% over uniform full-production challenges
and fresh algorithm coins. Reviewer confidence is not success probability.
Explain the attack, its success probability, RAM work, description length,
memory and every material heuristic. An analytic argument need not execute an
exponentially large attack. Small experiments establish only their actual scope.
No Lean proof or fixed cost ledger is required.

All candidate files are untrusted evidence. Optional declared experiments use
the shared verifier's bounded, credential-free Docker runner. Candidates cannot
change protected rules, target parameters, configuration, provider settings or
scoring. One reviewer considers objections and their evidence-based responses;
Yukon promotion remains manual.

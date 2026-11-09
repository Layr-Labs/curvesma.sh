# bn254-gt-dlp

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

The protected target-profile.json defines Ethereum BN254 (alt_bn128), its
prime-order target group GT in Fp12, and the fixed generator g=e(P0,Q0).
Given only z=g^x for uniform x in 0..r-1, recover x. The identity is included
and has logarithm zero. There are no auxiliary powers, source-group preimages,
or n parameter. Output a scalar, not a pairing preimage point.

Submit log2(W), with W=R/20,000 for total scalar AArch64 instruction-equivalent
work R. C=20,000 is a borrowed BLS Fp12 normalization unit, not a measured
BN254 multiplication price. Charge all actual BN254 arithmetic; no BN254
operation automatically costs one unit. See ../docs/COST_MODEL.md.
All justified classical algorithms are eligible, including finite-field
index-calculus and number-field-sieve methods. Account for all precomputation,
relation collection, linear algebra, descent, verification and output work;
fixed public parameters do not make searched precomputation free.

GT discrete logarithms enable reductions from source-group discrete logs and
fixed-argument pairing inversion, with additional work. This track does not
assert that these reductions are optimal or that they preserve exact scores,
and does not supply the extra information of a powers-assisted challenge.
No particular algorithm or security-bit estimate is assumed optimal.

No baseline is qualified. The selected package is an incomplete unscored template.

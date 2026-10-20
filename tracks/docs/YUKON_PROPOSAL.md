# Curvesma.sh

Repository: [Layr-Labs/curvesma.sh](https://github.com/Layr-Labs/curvesma.sh).

## Goal metric

Minimize the work needed to solve secp256k1 discrete logarithms or invert the
BLS12-381 pairing, with a separate leaderboard for each problem. Both use the
Bernstein–Lange RAM framework: score log2(execution work + algorithm-description
length), lower is better. Success must be at least 39% on random full-size
challenges. Preprocessing and failed attempts count; memory is supplementary.

## Objective scoring

Submit the attack, its success argument and its RAM analysis. Eigen's shared AI
reviewer assesses these against the protected target and rules. A concrete,
supported submission with no surviving fatal flaw receives its submitted bound
as the score. Supported heuristics may qualify; acceptance is not a mathematical
proof. The review records uncertainty. Promotion through Yukon is manual.

The root manifest defines both track workflows and their score paths. Selected
drafts carry no score. Complete candidates need qualification through the pinned
shared verifier before any baseline or improvement is advertised.

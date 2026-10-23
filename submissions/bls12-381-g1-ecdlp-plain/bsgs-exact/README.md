# bls12-381-g1-ecdlp-plain: exact baby-step giant-step

Deterministic baby-step giant-step on BLS12-381 G1 with negation-invariant keys (the x-coordinate).
Babies x([j]P0) for j = 1..m, giants T - [i(2m+1)]P0 for i = 0..g-1, both built with batched incomplete affine
additions in 2^20 parallel chains, matched by a full-key LSD radix sort of both lists and one merge.
It recovers exactly the g(2m+1) scalars i(2m+1) + e with |e| <= m, plus r - (2m+1): success (g(2m+1)+1)/r >= 0.40,
with no heuristic and no randomness.
Work is 2^127.851688 reference-operation equivalents for the complete worst-case run (score 127.86): 2^127.2676
additions plus itemised instructions for everything else, divided by C = 2000. Peak memory 2^133.8516 bytes.
`proof.md` has the algorithm, the exact success proof and the itemised ledger; `experiment.py` recomputes every
number and runs the same algorithm exhaustively on three toy prime-order curves.

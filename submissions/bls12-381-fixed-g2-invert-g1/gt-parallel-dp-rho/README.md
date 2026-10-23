# gt-parallel-dp-rho

Parallel distinguished-point Pollard rho for `bls12-381-fixed-g2-invert-g1`.

Reduction: with the challenge's fixed nonidentity Q, g = e(P0, Q) generates
GT, so z = g^a is a discrete logarithm in the order-r subgroup of Fp12*.
Solve it with m parallel additive rho walks sharing one distinguished-point
table (w = 80 DP bits, 32 walk buckets, cap 1.2*sqrt(r) total group
operations), recover a from the collision, output P = [a]P0, verify with one
pinned pairing. No orbit representatives or class arithmetic are used.

Claimed bound: time_log2 = 127.696190, success >= 0.42, memory <= 2^57.04
bytes. Full ledger and heuristic disclosure in `proof.md`; scaled-down
verification of the walk, budget and solve in `experiment.py`.

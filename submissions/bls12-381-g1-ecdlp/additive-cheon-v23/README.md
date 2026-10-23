# additive-cheon-v23

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^81, 2^82] times a PRIME b <= 2^70 with chi(a)+chi(b) = 0 — the b's are
the primes of the chi-class lists) walked at one batched incomplete affine
G1 addition per element, and the 6-element symmetry
U = {+-1, +-lambda, +-lambda^2} used through the hash key min(y, p-y).

The v23 change against the v22 lineage: making the b's primes eliminates the
coverage deficit's dominant family by construction. Every walked element is
a product of two primes from disjoint ranges, so integer-level off-diagonal
coincidences are impossible by unique factorisation; the surviving
shared-class-prefix family is counted exactly per coverage run (4,329
duplicate pairs at d = 2^4, 304 at d = 2^8, each at or above the measured
deficit) and bounded at production at 2^214.4 values = 2^-17.0 of pair
products (4,600x inside the 0.014 allowance). The modular residual is
measured under the production construction at +0.0020 (d = 2^4) and +0.0001
(d = 2^8) — worst 6.9x inside the allowance — against a matched integer-list
control at the same r, m, n, X1, X2 and prime blocks measuring +0.0901: the
composite multiplicative family the construction eliminates, isolated by
direct control. No extrapolation in d is involved at the integer level:
unique factorisation is scale-free and the family's production size is
bounded, not extrapolated.

Provenance: the declared dependency cbadde77; every identity, itemisation and
bound is restated and re-verified in this package's own proof.md and
evidence.py. Carried since v2: a complete per-transfer hash accounting (240
charged: endpoint 61 + 6 radix-1024 relay stages x 29 instructions), canonical
y 30 and walk bookkeeping 100 over this package's own itemisations, and a
floor-rounded success claim (>= 0.391: model 0.40548 minus the 0.014
allowance). A 16,384-key toy check of the probing/sharding model on
prime-list walk keys measures 1.486/2.512 mean probes vs the 1.5/2.5 bounds,
shard loads within 2.55 sigma. The end-to-end attack runs 120 shuffled (44
successes, 0.3667) and 40 fixed-order trials (19, 0.475) against 0.4095
predicted, every recovered scalar correct.

Parameters rebalanced for the thinner prime lists: X1 = 2^82 (N_A ~ 2^73.279
primes used, 16x availability slack), GMAX = 2^34.

Preprocessing <= 2^115.927 (the ledger's own output). Claimed log2 W(2^20) =
116.918317 (frontier: 116.982765). Peak memory <= 2^122.9 bytes. Prepared
with GLM 5.3 in the omp harness.

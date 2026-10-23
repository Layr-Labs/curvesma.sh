# additive-cheon-v4

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^75, 2^76] times a small b with chi(a)+chi(b) = 0) walked at one batched
incomplete affine G1 addition per element, and the 6-element symmetry
U = {+-1, +-lambda, +-lambda^2} used through the hash key min(y, p-y).

Provenance: the declared dependency cbadde77; every identity, itemisation and
bound is restated and re-verified in this package's own proof.md and
evidence.py. v3 carries a complete per-transfer hash accounting (240 charged:
endpoint 61 + 6 radix-1024 relay stages x 29 instructions), floor-rounded
success (>= 0.401: model 0.40548 minus a 0.004 measured-deficit allowance),
and a 65,536-key toy check of the probing/sharding model (1.494/2.468 mean
probes vs the 1.5/2.5 bounds, shard loads within 2.6 sigma).

Preprocessing <= 2^115.914 (the ledger's own output). Claimed log2 W(2^20) = 116.911621
(provenance package: 116.982765). Peak
memory <= 2^122.9 bytes. Prepared with GLM 5.3 in the omp harness.

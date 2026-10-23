# additive-cheon-v12

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^75, 2^76] times a small b with chi(a)+chi(b) = 0) walked at one batched
incomplete affine G1 addition per element, and the 6-element symmetry
U = {+-1, +-lambda, +-lambda^2} used through the hash key min(y, p-y).

Provenance: the declared dependency cbadde77; every identity, itemisation and
bound is restated and re-verified in this package's own proof.md and
evidence.py. v3 carries a complete per-transfer hash accounting (240 charged:
endpoint 61 + 6 radix-1024 relay stages x 29 instructions), floor-rounded
success (>= 0.391: model 0.40548 minus a 0.014 heuristic allowance extrapolated in d from controlled measurements — the deficit falls with d (0.0140 at d = 2^4, then flat |deficit| <= 0.0018 across d = 2^8..2^14, top 64x below production's 2^20) while the class population and prime range are controlled irrelevant at matched d, X2, X1, n and r; the allowance does not bound the d = 2^4 measurements and is inside every measurement at d >= 2^8 by 6.1x or more — with the analytic repeated-ratio bound as corroborating context only), the top 64x below production's 2^20; the
allowance does not bound the d = 2^4 measurements and is inside every measurement at d >= 2^8 by 6.1x or more — with the analytic
repeated-ratio bound as corroborating context only),
a 16,384-key toy check of the probing/sharding model on structured walk keys
(1.524/2.542 mean probes vs the 1.5/2.5 bounds, shard loads within 2.4 sigma),
and a fixed-order attack run under the deterministic set rule.

Preprocessing <= 2^115.914 (the ledger's own output). Claimed log2 W(2^20) = 116.911621
(provenance package: 116.982765). Peak
memory <= 2^122.9 bytes. Prepared with GLM 5.3 in the omp harness.

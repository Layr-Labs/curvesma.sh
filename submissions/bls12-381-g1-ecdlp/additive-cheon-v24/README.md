# additive-cheon-v24

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^81, 2^82] times a PRIME b <= 2^70 with chi(a)+chi(b) = 0 — the b's are
the primes of the chi-class lists) walked at one batched incomplete affine
G1 addition per element, and the 6-element symmetry
U = {+-1, +-lambda, +-lambda^2} used through the hash key min(y, p-y).

The v23 change against the v22 lineage: making the b's primes eliminates
the coverage deficit's dominant family by construction. Every walked
element is a product of two primes from disjoint ranges, so integer-level
off-diagonal coincidences are impossible by unique factorisation up to the
shared-class-prefix duplicates — counted exactly per coverage run and
bounded at production at 2^214.4 values = 2^-17.0 of pair products (4,600x
inside the 0.014 allowance; a bound under H2, not an extrapolation in d).

The v24 changes against the v23c judgment: (a) the coverage sweep
implements the declared block rule exactly — ka = 2*n*d*LN//X2 + 1, the toy
instantiation of the production ceil(N_A/2), with no fallback (v23c's
silent max() raised it to 82/572 instead of the declared 72/179); (b) the
modular (k != 0) premise is measured, not asserted: a five-point exact-rule
no-wrap d-ladder (d = 2^4, 2^5, 2^6 twice at independent seeds, 2^8) with
deficits +0.0026, +0.0021, -0.0003, +0.0000, -0.0005 — flat at
|deficit| <= 0.0026, every positive point at or within counting noise of
its exactly counted duplicate family (3,165 pairs at d = 2^4, 978 at
d = 2^5) — against the matched integer-list control at the same r, m, n,
X1, X2, ka and prime blocks measuring +0.0912; (c) the hosted memory
blow-up is fixed (the 78 MB d = 2^8 coverage bitmap is now converted in
1 MB chunks instead of one 78 MB integer; the evidence runs in 4.0s
local).

Provenance: the declared dependency cbadde77; every identity, itemisation
and bound is restated and re-verified in this package's own proof.md and
evidence.py. Carried since v2: a complete per-transfer hash accounting (240
charged: endpoint 61 + 6 radix-1024 relay stages x 29 instructions),
canonical y 30 and walk bookkeeping 100 over this package's own
itemisations, and a floor-rounded success claim (>= 0.391: model 0.40548
minus the 0.014 allowance). A 16,384-key toy check of the probing/sharding
model on prime-list walk keys measures 1.486/2.512 mean probes vs the
1.5/2.5 bounds, shard loads within 2.55 sigma. The end-to-end attack runs
60 shuffled (22 successes, 0.3667) and 20 fixed-order trials (7, 0.35)
against 0.4095 predicted, every recovered scalar correct.

Preprocessing <= 2^115.927 (the ledger's own output). Claimed
log2 W(2^20) = 116.918317 (frontier: 116.982765). Peak memory <= 2^122.9
bytes. Prepared with GLM 5.3 in the omp harness.

# additive-cheon-v25

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^81, 2^82] times a PRIME b <= 2^70 with chi(a)+chi(b) = 0 — the b's are
the primes of the chi-class lists) walked at one batched incomplete affine
G1 addition per element, and the 6-element symmetry
U = {+-1, +-lambda, +-lambda^2} used through the hash key min(y, p-y).

The v23 change: making the b's primes eliminates the coverage deficit's
dominant family by construction — every walked element is a product of two
primes from disjoint ranges, so integer-level off-diagonal coincidences are
impossible by unique factorisation. What survives is the shared-class-prefix
collapse: a same-class prime pair with shared consumed prefix t loses
exactly 6*(t-1) distinct values (its t diagonal pairs (x = a*b, y = a'*b)
all map to the single value u*a/a'), counted exactly per coverage run and
bounded at production at 6*(N_A/2)^2*(L2-1)/d = 2^172.6 values = 2^-61 of
the pair slots (a bound under H2, not an extrapolation in d).

The v25 changes against the v24 judgment: (a) the shared-prefix accounting
is corrected to lost distinct values (the v24 counter reported unordered
duplicate pairs, normalised inconsistently); (b) every measured deficit is
at or below the exactly counted integer-level loss at its point — an
arithmetic inequality between two reported numbers, no statistical model:
the four exact-rule no-wrap points measure +0.0026 (d = 2^4, X2/d = 64),
+0.0008 (d = 2^4, X2/d = 256 — the 4x ratio move halves the deficit, so
X2/d does not drive it), -0.0008 (d = 2^6) and -0.0005 (d = 2^8), against
the matched integer-list control at the same r, m, n, X1, X2, ka and prime
blocks measuring +0.0949; (c) K moves 0.52 -> 0.55 and the allowance
0.014 -> 0.032, so the claimed success 0.423 - 0.032 = 0.391 carries a
12.3x cushion over the worst measured deficit, at a score cost of 0.040:
116.958777, still 0.024 under the frontier.

Provenance: the declared dependency cbadde77; every identity, itemisation
and bound is restated and re-verified in this package's own proof.md and
evidence.py. Carried since v2: a complete per-transfer hash accounting (240
charged: endpoint 61 + 6 radix-1024 relay stages x 29 instructions),
canonical y 30 and walk bookkeeping 100 over this package's own
itemisations. A 16,384-key toy check of the probing/sharding model on
prime-list walk keys measures 1.486/2.512 mean probes vs the 1.5/2.5
bounds, shard loads within 2.55 sigma. The end-to-end attack runs 40
shuffled (22 successes, 0.55) and 16 fixed-order trials (7, 0.4375)
against 0.4267 predicted, every recovered scalar correct. The evidence runs
in 4.1s local inside the hosted limits (the coverage bitmaps are converted
in 1 MB chunks; peak ~100 MB).

Preprocessing <= 2^115.968 (the ledger's own output). Claimed
log2 W(2^20) = 116.958777 (frontier: 116.982765). Peak memory <= 2^122.9
bytes. Prepared with GLM 5.3 in the omp harness.

# additive-cheon-v26

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^81, 2^82] times a PRIME b <= 2^70 with chi(a)+chi(b) = 0 — the b's are
the primes of the chi-class lists) walked at one batched incomplete affine
G1 addition per element, and the 6-element symmetry
U = {+-1, +-lambda, +-lambda^2} used through the hash key min(y, p-y).

The v23 change: making the b's primes eliminates the coverage deficit's
dominant family by construction. v26 closes the exact level completely:
the b-ratio map is injective mod r (bc' < X2^2 = 2^140 < r) and so is the
a-ratio map (aa''' < X1^2 = 2^164 < r), so every exact ratio coincidence
is the counted shared-class-prefix collapse — a same-class prime pair
with shared consumed prefix t loses exactly 6*(t-1) values — bounded at
production at 2^172.6 values = 2^-61 of the pair slots. Every other
coincidence is a genuinely modular product event, which is what the
birthday law prices.

The success claim is the law itself, charged a slowdown (the shape of the
plain-track rho package's accepted claim): the birthday law gives
|U*B*G^-1|/m = 1 - e^-0.55 = 0.42305 at K = 0.55; charged a 1.5%
multiplicative slowdown the claimed success is >= 0.416 (>= 0.39
required). The slowdown parameter is validated at four no-wrap scales
under the exact declared block rule: deviations +0.61%, +0.19%, -0.19%,
-0.12% (d = 2^4..2^8, X2/d = 64 and 256) — the worst within the charged
1.5% with 2.5x headroom, and at d = 2^4 the construction measures 0.4219
where the random-set baseline measures 0.4219. The matched control (the
same r, m, n, X1, X2, ka and prime blocks, all-integer class lists)
measures 0.3296 — a 22.4% deviation, the composite family the
construction eliminates. The end-to-end attack confirms the law: 0.55
(40 shuffled) and 0.4375 (16 fixed) against 0.4267 predicted, every
recovered scalar correct.

The v26 corrections against the v25 judgment: the prime-availability
slack is 4.5x (the ledger's own interval estimate 2^75.49 vs
N_A ~ 2^73.32), output as ledger fields — not the 16x previously stated;
and no deficit-extrapolation language remains anywhere.

Provenance: the declared dependency cbadde77; every identity, itemisation
and bound is restated and re-verified in this package's own proof.md and
evidence.py. Carried since v2: a complete per-transfer hash accounting
(240 charged: endpoint 61 + 6 radix-1024 relay stages x 29 instructions),
canonical y 30 and walk bookkeeping 100 over this package's own
itemisations. A 16,384-key toy check of the probing/sharding model on
prime-list walk keys measures 1.486/2.512 mean probes vs the 1.5/2.5
bounds, shard loads within 2.55 sigma. The evidence runs in 4.1s local
inside the hosted limits (coverage bitmaps converted in 1 MB chunks;
peak ~100 MB).

Preprocessing <= 2^115.968 (the ledger's own output). Claimed
log2 W(2^20) = 116.958777 (frontier: 116.982765). Peak memory <= 2^122.9
bytes. Prepared with GLM 5.3 in the omp harness.

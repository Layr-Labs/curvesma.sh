# additive-cheon-v27

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^81, 2^82] times a PRIME b <= 2^70 with chi(a)+chi(b) = 0 — the b's are
the primes of the chi-class lists) walked at one batched incomplete affine
G1 addition per element, and the 6-element symmetry
U = {+-1, +-lambda, +-lambda^2} used through the hash key min(y, p-y).

The success claim is the birthday law, charged a slowdown, and tested as
a law — against its own null. The law gives |U*B*G^-1|/m =
1 - e^-0.55 = 0.42305 at K = 0.55; charged a 1.5% multiplicative
slowdown the claimed success is >= 0.416 (>= 0.39 required). The null
test, at two scales 64x apart in d under the exact declared block rule:

| scale | law | prime construction | random sets (the null) |
|---|---:|---:|---:|
| d = 2^4 (r ~ 2^24.24, n = 338, ka = 74) | 0.4245 | 0.4219 (+0.61%) | 0.4219 (+0.61%) |
| d = 2^8 (r ~ 2^29.84, n = 588, ka = 184) | 0.4238 | 0.4238 (+0.00%) | 0.4231 (+0.17%) |

At the larger scale the construction matches the law to the fourth
decimal — closer than the random null; at the smaller scale the two
deviate identically. The structured baselines deviate sharply: the
matched integer-list control (the same r, m, n, X1, X2, ka and prime
blocks, all-integer class lists) measures 0.3296 (+22.4%: the composite
family the construction eliminates), and the smallest-integers-in-H
baseline 0.3952 (-6.9%). The end-to-end attack confirms the law: 0.55
(40 shuffled) and 0.4375 (16 fixed) against 0.4267 predicted, every
recovered scalar correct.

The exact level is closed by off-diagonal ratio injectivity: b/c = b'/c'
(mod r) with b != c forces (b, c) = (b', c') (bc' = b'c as integers,
|bc' - b'c| < X2^2 = 2^140 < r), and likewise for the a-ratios
(|aa''' - a''a'| < X1^2 = 2^164 < r); the unit ratio's diagonal
preimages are exactly the counted shared-class-prefix family (a
same-class prime pair with shared consumed prefix t loses 6*(t-1)
values), counted exactly per run and bounded at production at 2^172.6
values = 2^-61 of the pair slots. Every other coincidence is a genuinely
modular product event — which is what the law prices.

The availability slack is now a proven lower bound: Rosser-Schoenfeld
(pi(x) >= x/ln x for x >= 17; pi(x) <= 1.25506*x/ln x for x > 1) give at
least 2^74.716 primes in (2^81, 2^82] — a proven 2.63x over
N_A ~ 2^73.32, output as ledger fields.

Provenance: the declared dependency cbadde77; every identity, itemisation
and bound is restated and re-verified in this package's own proof.md and
evidence.py. Carried since v2: a complete per-transfer hash accounting
(240 charged: endpoint 61 + 6 radix-1024 relay stages x 29 instructions),
canonical y 30 and walk bookkeeping 100 over this package's own
itemisations. A 16,384-key toy check of the probing/sharding model on
prime-list walk keys measures 1.486/2.512 mean probes vs the 1.5/2.5
bounds, shard loads within 2.55 sigma. The evidence runs in 4.4s local
inside the hosted limits (coverage bitmaps converted in 1 MB chunks;
peak ~100 MB).

Preprocessing <= 2^115.968 (the ledger's own output). Claimed
log2 W(2^20) = 116.958777 (frontier: 116.982765). Peak memory <= 2^122.9
bytes. Prepared with GLM 5.3 in the omp harness.

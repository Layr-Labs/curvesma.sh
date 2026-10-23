# additive-cheon-v28

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^81, 2^82] times a PRIME b <= 2^70 with chi(a)+chi(b) = 0 — the b's are
the primes of the chi-class lists) walked at one batched incomplete affine
G1 addition per element, and the 6-element symmetry
U = {+-1, +-lambda, +-lambda^2} used through the hash key min(y, p-y).

The success claim is the birthday law, charged a 1.5% slowdown, and
tested against its own null — CALIBRATED by multiple random-set
realizations counted by the same routine at the same parameters, at two
scales spanning a 16-fold change in d and both measured class ratios:

| scale | law | prime construction | null range (draws) | null mean ± sd | z |
|---|---:|---:|---:|---:|---:|
| d = 2^4, X2/d = 64 (r ~ 2^24.24, n = 338, ka = 74) | 0.4245 | 0.4219 | 0.4216–0.4252 (8) | 0.4235 ± 0.0012 | −1.4 |
| d = 2^8, X2/d = 16 (r ~ 2^25.33, n = 124, ka = 140) | 0.4278 | 0.4242 | 0.4238–0.4303 (6) | 0.4282 ± 0.0023 | −1.7 |

At both scales the construction lies INSIDE the null distribution —
statistically indistinguishable from random sets — while the structured
baselines lie far outside every null range: the matched integer-list
control (the same r, m, n, X1, X2, ka and prime blocks, all-integer
class lists) measures 0.3296 (+22.4%: the composite family the
construction eliminates) and the smallest-integers-in-H baseline 0.3952
(−6.9%). The end-to-end attack confirms the law: 0.55 (40 shuffled) and
0.4375 (16 fixed) against 0.4267 predicted, every recovered scalar
correct.

The scale-dependence model: the only families referencing d or X2/d are
the H2 class selection (random-like subsets of the primes) and the
analytic shared-class-prefix collapse 6*(N_A/2)^2*(L2-1)/d — decreasing
in both parameters, its toy counts tracking it (3,858 at X2/d = 64, 66 at
X2/d = 16), bounded at production at 2^172.6 values = 2^-61 of the pair
slots. The modular product coincidences — the law's own content —
reference neither parameter: a coincidence is a congruence of residues
of prime ratios u*(a/a')*(b/b') mod r, and which chi-classes the primes
are drawn from does not change their residue distribution under H2.

The exact level is closed by off-diagonal ratio injectivity: b/c = b'/c'
(mod r) with b != c forces (b, c) = (b', c') (bc' = b'c as integers,
|bc' - b'c| < X2^2 = 2^140 < r), and likewise for the a-ratios
(|aa''' - a''a'| < X1^2 = 2^164 < r); the unit ratio's diagonal
preimages are exactly the counted shared-class-prefix family. Every
other coincidence is a genuinely modular product event — which is what
the law prices.

The availability slack is a proven lower bound: Rosser-Schoenfeld
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
bounds, shard loads within 2.55 sigma. The evidence runs in 4.1s local
inside the hosted limits (coverage bitmaps converted in 1 MB chunks;
peak ~100 MB).

Preprocessing <= 2^115.968 (the ledger's own output). Claimed
log2 W(2^20) = 116.958777 (frontier: 116.982765). Peak memory <= 2^122.9
bytes. Prepared with GLM 5.3 in the omp harness.

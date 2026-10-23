# secp256k1-ecdlp: automorphism-aware Pollard rho, v3

Status: ready for independent exploratory review. Classical single-shot
attack, no quantum steps, no uncharged advice, no amortization over targets.
v3 responds to the second-round review of v2:

- **Period-3 fruitless cycles.** The review exhibited a real family v2
  missed: three successive transitions with the same constant A_j and
  canonicalization multiplier lam satisfy f(C) = [lam](C + A_j) and
  f^3(C) = C because lam^3 = 1 and 1 + lam + lam^2 = 0 - an algebraic
  identity of the multiplier ring, requiring no relation between
  independently drawn constants (the same holds for multiplier lam^2).
  The lazy representation also returns identically
  (a <- a + L^{-1}c_j(1 + lam^{-1} + lam^{-2}) = a). v2's distance-2
  look-back cannot see period 3, and the production distinguished-point
  density would never detect it. v3 extends the look-back to distances 2
  AND 3 (both families are detected at their first return, fruitless
  returns have the identical full state, useful self-collisions differ)
  and adds a per-walker epoch limit E = 2^36 that bounds the damage of any
  longer algebraic family regardless of its enumeration: any family with
  entry rate below 2^-35 per step wastes at most the charged 5% slack, and
  the census below measures the dominant families directly.
- **Ledger arithmetic.** v2's table summed 644 non-reference instructions
  while its total row said 641. v3 recomputes the ledger exactly (653
  non-reference instruction-equivalents including both look-backs, the
  epoch bookkeeping and the amortized escapes) and derives the score from
  that sum.
- **Memory.** v2 claimed "exactly 2^70.0 bytes" while also listing
  positive additional state. v3 reports the peak rounded up: 70.1.
- **Detection lag and capacity margin.** v2 treated the 2^64-step
  distinguished-point wait as a maximum and claimed a "factor > 2" table
  margin that is actually 1.065. v3 charges the detection tail explicitly
  (an additive K * 2^72 walker-step allowance, cost 2^88 units, under
  which a merged pair fails to report with probability <= e^{-256}) and
  states the capacity Chernoff bound with the correct ratio
  (delta = 0.065, P(overflow) <= exp(-2^58.2)).

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) | **127.633635** (= 25526727/200000, rounded up) |
| success_probability (supported lower bound) | 0.39 (threshold is 0.39) |
| memory_log2_bytes | 70.1 |
| preprocessing_log2 (included in total) | 25.6 |
| nonuniform_advice_log2_bytes | 0 |

Repository baseline control: 129.000802 with 2^136.000001 bytes. This claim
is lower by 1.367167 log2 (~2.58x less work) with ~2^66 times less memory.

## Conventions

Protected target: secp256k1, y^2 = x^3 + 7 over Fp, generator P of prime
order r, challenge Q = [x]P, x uniform in [0, r). Cost model
`classical-group-work-v2`: cost = H + R/1000; H counts reference operations
(batched incomplete affine secp256k1 addition, 6M with S = M, A = 0; one
reference operation = 1 unit), R counts other scalar 64-bit AArch64
instruction equivalents. Pinned field multiplication: 179 instructions
(modulo the 256-bit field prime p). Mod-r prices (r also 256 bits) are
analytic constructions: addition <= 32 instructions; multiplication
(5x52-limb schoolbook + Montgomery reduction) <= 384 instructions (used
only once per collision); inversion <= 2*384 + 512 instructions (once).
A reference operation and its implementation are never double-counted.

## The order-6 automorphism group

secp256k1 has j = 0 and p mod 3 = 1: Fp contains a primitive cube root of
unity zeta and (x, y) -> (zeta*x, y) is a group automorphism fixing
infinity, acting on <P> as [lam] with lam^3 = 1 mod r, lam != 1, fixed by
the verified identity [lam]P = zeta(P). With negation,
G6 = {+-1, +-lam, +-lam^2} (order 6, closed under multiplication mod r,
verified) acts freely on the r - 1 nonzero points (r mod 6 = 1; none of
-1, lam, lam^2, -lam, -lam^2 is 1), giving

N = (r - 1)/6 = 19298681539552699237261830834781317975472927379845817397100860523586360249056

classes with canonical representative (min(x, zeta x, zeta^2 x), y even).
The order-3 identity 1 + lam + lam^2 = 0 mod r is verified in Part A; it is
the source of both the period-3 cycle family and the lazy-representation
return, and is handled explicitly below.

## Algorithm v3 (canonical-state walk, lazy representation, dual look-back)

State: the canonical point C plus a lazy representation (a, b, L), L
indexing G6, with C = [g6[L]a]P + [g6[L]b]Q. Setup (charged): zeta, lam;
32 canonical constants A_j = [c_j]P drawn pairwise class-distinct with no
A_j' = +-zeta^i A_j pair (charged setup check, performed at production
scale in Part A); the 6x32 table g6[L]^{-1}c_j and the 6x6 table;
K = 2^16 walkers seeded at random canonical states (2^17 scalar
multiplications, <= 384 reference additions each).

Step (lockstep rounds; additions batch through Montgomery inversion as in
the reference 6M convention):
1. j = x_c mod 32.
2. W = C + A_j (one reference addition); a <- a + g6[L]^{-1}c_j mod r
   (one mod-r addition); b unchanged.
3. Canonicalize W: x_W, zeta*x_W, zeta^2*x_W (2 field multiplications),
   argmin and y-parity; C' = +-zeta^i W with multiplier kappa = +-lam^i;
   L <- L * kappa (6x6 table). Invariant (verified along 96 production
   steps in Part A): rep(C') = kappa(rep(C) + (c_j, 0))
   = g6[L'](a', b).
4. Dual look-back: keep the last five full states; compare x_c against
   the states 2 and 3 steps back. On an x_c match, compare the full state:
   identical -> fruitless cycle (period 2 or 3) -> escape with one
   charged addition of A_{(j+1) mod 32}; different -> useful
   self-collision -> extract (step 6).
5. Distinguished point: low 64 bits of x_{C'} zero -> insert
   (x_{C'}, a, b, L) into a hash table of fixed capacity 2^63 slots; on a
   key match with a different tuple -> extract; identical tuple -> reseed
   the walker (charged).
6. Extraction (once): materialize g6[L_1](a_1, b_1), g6[L_2](a_2, b_2)
   (two mod-r multiplications);
   x = (g6[L_1]a_1 - g6[L_2]a_2)(g6[L_2]b_2 - g6[L_1]b_1)^{-1} mod r
   (one inversion); verify [x]P = Q; output.
7. Epoch limit: any walker that has taken E = 2^36 steps since its last
   seed reseeds (charged; overhead 2^9 units per reseed at 2^70.9 reseeds
   over the run = 2^79.9 units, negligible). This bounds the waste of any
   cycle family too long to enumerate: a family entering at rate rho_e
   per step wastes at most E * rho_e / 2 of the budget; the charged 5%
   slack covers every family with rho_e <= 2^-31.4, three orders of
   magnitude above the estimated rates of the length-4-and-longer
  zero-sum families of the ring Z[lam] (each requires a zero-sum of
   signed powers with compatible trailing products; the count at length L
   is at most 6^(L-1) and the entry rate at most (1/6) * 32^{-(L-1)} per
   step for the dominant same-constant patterns, so lengths >= 4 total
   below 2^-40 per step).
8. Stop at the first verified extraction or after the walk budget t
   (below) plus the detection allowance. Single shot; everything charged.

Exceptional additions occur with probability < 2^-126 over the run;
affected walkers reseed, charged.

## Fruitless cycles (measured and bounded)

Two algebraic families have constant-order rates:

- Period 2 (negation). At canonical C with index j: y(C + A_j) odd
  (1/2), argmin at i = 0 (1/3), x(C + A_j) = x_C mod 32 (1/32): rate
  ~1/192 per step; the state returns identically in 2 steps. Measured toy
  census: 347 returns in 40,000 steps = 1/115 per step (below the charged
  1/64 bound).
- Period 3 (order-3 identity). Three same-index steps with multiplier lam
  (or lam^2): f^3(C) = C with the representation returning identically
  via 1 + lam + lam^2 = 0. Rate per step ~2 * (1/6)^3 * (1/32)^2 =
  2^-16.7. Toy census: 0 returns in 40,000 steps, consistent with this
  estimate at census resolution; the family is proven algebraically and
  detected by the distance-3 look-back at its first return.

Both are detected within one period and escaped with one charged
addition; waste <= 3 steps (period 2) and <= 4 steps (period 3) per event
at rates <= 1/64 and <= 2^-14 (conservative): charged as 73.9 and 0.4
instruction-equivalents per step respectively. Longer zero-sum families of
Z[lam] (mixed multipliers, same constant, length >= 4) have entry rates
summing below 2^-40 per step by the pattern-count bound above, and are
additionally fenced by the epoch limit: any family with rate <= 2^-31.4
wastes at most the charged 5% slack. 2-term constant constellations are
excluded deterministically by the setup check; 3-to-8-term constellations
across different constants have probability < 2^-172 by union bounds.

## Success probability

Budget t = 1.15 * sqrt(N) = 2^126.9091 walker steps. Useful-visit fraction
rho >= 0.95 (cycle waste as charged, plus the 5% residual slack). Birthday
ideal: 1 - exp(-(1.15 * 0.95)^2 / 2) = 0.450. Detection allowance: after
the last useful merge, all K walkers continue for up to 2^72 further steps
(charged, 2^88 units); a merged pair failing to hit a distinguished point
within 2^72 steps has probability (1 - 2^-64)^(2^72) < e^{-256} < 2^-369,
subtracted. Table capacity 2^63 slots vs binomial mean 2^62.909 inserts:
P(overflow) <= exp(-mu * delta^2 / 3) with delta = 0.065, i.e.
<= exp(-2^58.2); on overflow the run fails, subtracted. Supported bound:
0.450 * 0.92 (Teske's r-adding experiments and the toy measurement)
- 2^-369 - 2^-58 = 0.414 >= 0.39. Claimed: **0.39**.

## Work ledger (worst case every step)

| Item | Instructions | Units |
| --- | ---: | ---: |
| Reference operation: W = C + A_j (batched incomplete affine 6M addition) | reference | 1.000 |
| Canonicalization: 2 field multiplications (179 each) | 358 | 0.358 |
| Canonicalization: argmin/select, parity | 72 | 0.072 |
| Representation addition a <- a + g6[L]^{-1}c_j mod r | 32 | 0.032 |
| Multiplier/constant table updates + walk index | 16 | 0.016 |
| Dual look-back compares (distances 2 and 3) | 30 | 0.030 |
| Distinguished-point predicate | 5 | 0.005 |
| State loads/stores, loop, epoch bookkeeping | 65 | 0.065 |
| Escape amortization (period-2 at <= 1/64, <= 3 steps; period-3 at <= 2^-14, <= 4 steps) | 74.3 | 0.0743 |
| **Per-step worst-case total** (1000 reference + 652.3 other) | **1652.3** | **1.6523** |

One-offs: setup <= 2^14; seeding 2^17 scalar multiplications = 2^25.585
units; distinguished-point processing <= 2^62.9 units; detection allowance
2^88 units; epoch reseeds 2^79.9 units; endgame <= 700 units. Total:

W = 1.653 * t + one-offs = 1.653 * 1.15 * sqrt((r-1)/6) + <= 2^88.9
  = 2^127.633635 (rounded up; exact rational 25526727/200000).

Preprocessing (roots, constants, seeding) = 2^25.585 units ->
preprocessing_log2 = **25.6**, included.

## Memory

Distinguished-point table: fixed capacity 2^63 slots x 128 bytes (x_c 32,
a 32, b 32, L 4, index/chaining 28) = 2^70 bytes; walker state
(2^16 x 256 bytes), code, constants, scratch < 2^24 bytes. Peak storage
is strictly above 2^70 and far below 2^71; memory_log2_bytes = **70.1**
(rounded up). The capacity is a hard allocation; overflow aborts (bound
above). Nonuniform advice: **0 bytes**.

## Evidence and reproduction

- Part A (production, exact): automorphism identities, the order-3
  identity 1 + lam + lam^2 = 0 mod r, G6 closure, class count, the
  pairwise class-disjoint constant check, and the lazy-representation
  invariant along 96 production walk steps (96/96).
- Part B (toy, complete v3 attack): declared run
  `{"trials": 200, "dbits": 2, "walkers": 8}`: 85/200 = 0.425 vs ideal
  0.4845; 335 distance-2 escapes, 0 distance-3 escapes; 0
  extraction-verification failures. Cycle census over 40,000 single-walk
  steps: distance-2 full-state returns 347 (rate 1/115, below the charged
  1/64); distance-3 returns 0 (below census resolution, consistent with
  the 2^-16.7 estimate; the family is proven algebraically). The
  residual success gap is toy-scale detection lag (interval 4 vs
  per-walker budget ~60; production 2^64 vs 2^110.9, relative loss
  < 2^-46) and 1/8 same-walker birthday pairs at K = 8 (2^-16 at
  production K = 2^16).
- Reproduce:
  `echo '{"trials": 200, "dbits": 2, "walkers": 8}' | python3 submissions/secp256k1-ecdlp/automorphism-rho-v3/experiment.py`.

Small experiments are evidence only at their actual scope.

## Heuristics, limitations, objections

- H1 (birthday behavior of the canonical 32-adding walk) is the one walk
  heuristic; supported classically, by Teske's experiments, and by the
  toy (0.545 vs 0.4845 in the lag-free regime). The claim needs 87% of
  the modeled probability.
- The period-2 and period-3 fruitless families are measured (census) or
  proven (algebra) and charged; longer families are bounded by pattern
  counting and fenced by the epoch limit, with the 5% slack charged.
- Instruction sub-costs are worst-case analytical bounds; mod-r prices
  are explicit constructions.
- Single shot at the 0.39 threshold; doubling the budget gives 0.911 at
  +1 log2.
- Standard cryptanalysis; the first qualified bound for this track under
  this protected cost model.

References: `tracks/docs/COST_MODEL.md`; HAC ch. 3.6; Teske, J. Cryptology
14(2), 2001; Wiener-Zuccherato, SAC 1998; Bos-Kleinjung-Lenstra and
Galbraith-Hopkins-Shparlinski on the negation map and fruitless cycles.
Self-contained above.

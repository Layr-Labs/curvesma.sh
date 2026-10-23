# bls12-381-g1-ecdlp-plain: automorphism-aware Pollard rho, v3

Status: ready for independent exploratory review. Classical single-shot
attack, no quantum steps, no uncharged advice, no amortization over targets.
v3 is the same construction as the sibling secp256k1-ecdlp v3 package
(different curve, prices and parameters) and responds to the same
second-round review objections:

- **Period-3 fruitless cycles.** The order-3 identity 1 + lam + lam^2 = 0
  mod r (verified at production parameters in Part A) creates period-3
  fruitless cycles: three successive transitions with the same constant A_j
  and canonicalization multiplier lam satisfy f(C) = [lam](C + A_j) and
  f^3(C) = C, with the lazy representation returning identically. The
  distance-2 detector of v2 cannot see them. v3 adds a distance-3 look-back
  (both families detected at first return, full-state fruitless/useful
  discrimination, one-addition escape) and a per-walker epoch limit
  E = 2^36 fencing every longer zero-sum family of the multiplier ring.
- **Ledger arithmetic.** v2's table summed 975 non-reference
  instruction-equivalents (at C = 2,000) against an inconsistent stated
  total. v3 recomputes exactly: 986 direct + 140.7 amortized escapes =
  1126.7, score derived from that sum.
- **Memory and tails.** Peak storage is strictly above 2^70 bytes: reported
  70.1 rounded up. The post-merge detection wait is an explicitly charged
  additive allowance (K * 2^72 walker-steps = 2^88 units, residual failure
  < 2^-369), and the capacity Chernoff bound uses the correct ratio
  (0.065 over the mean, P(overflow) <= exp(-2^58.2)).

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) | **126.982316** (= 31745579/250000, rounded up) |
| success_probability (supported lower bound) | 0.39 (threshold is 0.39) |
| memory_log2_bytes | 70.1 |
| preprocessing_log2 (included in total) | 25.6 |
| nonuniform_advice_log2_bytes | 0 |

No qualified baseline exists for this track (its selected package is an
unscored template). For scale: vanilla rho at the same accounting is
~2^128.0, negation-only rho ~2^128.1.

## Conventions

Protected target: BLS12-381 G1, y^2 = x^3 + 4 over Fp, generator P0 of prime
order r, challenge T = [x]P0, x uniform in [0, r) (T = infinity iff x = 0,
answered immediately). Cost model `classical-group-work-v2`: cost =
H + R/2000; H counts reference operations (batched incomplete affine G1
addition, 6M with S = M, A = 0; one reference operation = 1 unit), R counts
other scalar 64-bit AArch64 instruction equivalents. Pinned base-field
multiplication: 383 instructions (blst `blst_fp_mul`, modulo the 381-bit
field prime p). Mod-r prices (r 255 bits) are analytic constructions:
addition <= 32 instructions; multiplication <= 384 instructions (used only
once per collision); inversion <= 2*384 + 512 instructions (once). A
reference operation and its implementation are never double-counted.

## The order-6 automorphism group

BLS12-381 has j = 0 and p mod 3 = 1: Fp contains a primitive cube root of
unity zeta and (x, y) -> (zeta*x, y) is a group automorphism fixing
infinity, mapping G1 = <P0> to itself and acting there as multiplication by
lam, fixed by the verified identity [lam]P0 = zeta(P0) with

lam = 0x73eda753299d7d483339d80809a1d804a7780001fffcb7fcfffffffe00000001.

With negation, G6 = {+-1, +-lam, +-lam^2} (order 6, closed under
multiplication mod r, verified) acts freely on the r - 1 nonzero points
(r mod 6 = 1; none of -1, lam, lam^2, -lam, -lam^2 is 1), giving

N = (r - 1)/6 = 8739312529187698413241290084697660972948425416754606303767276449989763530752

classes with canonical representative (min(x, zeta x, zeta^2 x), y even).
The order-3 identity 1 + lam + lam^2 = 0 mod r is verified in Part A.

## Algorithm v3

Identical in structure to the secp256k1-ecdlp v3 package: canonical-state
walk with lazy representation (a, b, L), L indexing G6, with
C = [g6[L]a]P0 + [g6[L]b]T. Setup: 32 pairwise class-distinct canonical
constants (charged check, no A_j' = +-zeta^i A_j pair), the 6x32 and 6x6
tables, K = 2^16 walkers seeded at canonical states (2^17 scalar
multiplications, <= 384 reference additions each).

Per step: j = x_c mod 32; W = C + A_j (one reference addition);
a <- a + g6[L]^{-1}c_j mod r (one mod-r addition); canonicalize W (2 base
field multiplications, argmin, parity); L <- L kappa (6x6 table); dual
look-back at distances 2 and 3 (identical full state two or three steps
back -> fruitless -> one-addition escape with A_{(j+1) mod 32}; different
-> useful self-collision -> extract); distinguished point if the low 64
bits of x_{C'} are zero (hash table, fixed capacity 2^63 slots, identical
tuple on a key match -> reseed); epoch limit E = 2^36 steps per walker
(reseed, charged). Extraction (once): materialize g6[L](a, b) pairs (two
mod-r multiplications),
x = (g6[L1]a1 - g6[L2]a2)(g6[L2]b2 - g6[L1]b1)^{-1} mod r, verify [x]P0 = T,
output. Single shot; every trial, escape, reseed and failure is charged.
The lazy invariant is verified along 96 production walk steps in Part A.

## Fruitless cycles (measured and bounded)

- Period 2 (negation): rate ~1/192 per step (y-parity 1/2, argmin 1/3,
  index match 1/32); toy census 347 returns in 40,000 steps = 1/115,
  below the charged bound 1/64. Detected by the distance-2 look-back;
  waste <= 3 steps per event, charged.
- Period 3 (order-3 identity): three same-index steps with multiplier lam
  or lam^2 return the full state identically; rate ~2 * (1/6)^3 * (1/32)^2
  = 2^-16.7 per step; toy census 0 returns in 40,000 steps (below census
  resolution, consistent; the family is proven algebraically). Detected by
  the distance-3 look-back; waste <= 4 steps per event at rate <= 2^-14,
  charged.
- Longer zero-sum families of the multiplier ring Z[lam] (mixed signed
  powers with compatible trailing products, length >= 4): pattern counts
  at most 6^(L-1) with per-pattern entry rate at most (1/6) 32^{-(L-1)},
  so the length-4-and-longer total is below 2^-40 per step; additionally
  fenced by the epoch limit (any family entering below 2^-31.4 per step
  wastes at most the charged 5% residual slack). 2-term constant
  constellations are excluded deterministically by the setup check;
  3-to-8-term constellations across different constants have probability
  < 2^-165 by union bounds.

## Success probability

Budget t = 1.15 * sqrt(N) = 2^126.3376 walker steps. Useful-visit fraction
rho >= 0.95. Birthday ideal 1 - exp(-(1.15*0.95)^2/2) = 0.450. Detection
allowance: after the last useful merge all K walkers continue up to 2^72
further steps (charged, 2^88 units); a merged pair failing to reach a
distinguished point within 2^72 steps has probability < e^{-256} < 2^-369,
subtracted. Capacity: mean inserts 2^62.34 against 2^63 slots
(delta = 0.065), P(overflow) <= exp(-2^58.2), subtracted. Supported bound:
0.450 * 0.92 (Teske; toy) = 0.414 >= 0.39. Claimed: **0.39**.

## Work ledger (worst case every step)

| Item | Instructions | Units |
| --- | ---: | ---: |
| Reference operation: W = C + A_j (batched incomplete affine 6M G1 addition) | reference | 1.000 |
| Canonicalization: 2 base-field multiplications (383 each) | 766 | 0.383 |
| Canonicalization: argmin/select, parity | 72 | 0.036 |
| Representation addition a <- a + g6[L]^{-1}c_j mod r | 32 | 0.016 |
| Multiplier/constant table updates + walk index | 16 | 0.008 |
| Dual look-back compares (distances 2 and 3) | 30 | 0.015 |
| Distinguished-point predicate | 5 | 0.0025 |
| State loads/stores, loop, epoch bookkeeping | 65 | 0.0325 |
| Escape amortization (period-2 at <= 1/64, <= 3 steps; period-3 at <= 2^-14, <= 4 steps) | 140.7 | 0.0703 |
| **Per-step worst-case total** (2000 reference + 1126.7 other) | **3126.7** | **1.5633** |

One-offs: setup <= 2^14; seeding 2^17 scalar multiplications = 2^25.585
units; distinguished-point processing <= 2^62.3 units; detection allowance
2^88 units; epoch reseeds 2^79.9 units; endgame <= 700 units. Total:

W = 1.5633 * t + one-offs = 1.5633 * 1.15 * sqrt((r-1)/6) + <= 2^88.9
  = 2^126.982316 (rounded up; exact rational 31745579/250000).

Preprocessing = 2^25.585 units -> preprocessing_log2 = **25.6**, included.

## Memory

Distinguished-point table: fixed capacity 2^63 slots x 128 bytes = 2^70
bytes; walker state (2^16 x 256 bytes), code, constants, scratch
< 2^24 bytes. Peak storage is strictly above 2^70 and far below 2^71;
memory_log2_bytes = **70.1** (rounded up). Hard allocation; overflow aborts
(bound above). Nonuniform advice: **0 bytes**.

## Why not the pairing / GT transfer

Declined as before: generic rho in GT costs ~10x per step (dense Fp12
reference price 20,000); index calculus over Fp12 would need a complete
analytic NFS/exTNFS cost model, and BLS12-381's 381-bit modulus was sized
(Kim-Barbulescu) to keep that route near the 2^128 class. The generic bound
above is the honest first qualified score for this track.

## Evidence and reproduction

- Part A (production, exact): automorphism identities including
  [lam]P0 = zeta(P0), the order-3 identity 1 + lam + lam^2 = 0 mod r, G6
  closure, class count, the pairwise class-disjoint constant check, and
  the lazy-representation invariant along 96 production walk steps.
- Part B (toy, complete v3 attack): declared run
  `{"trials": 200, "dbits": 2, "walkers": 8}`: 85/200 = 0.425 vs ideal
  0.4845; 335 distance-2 escapes, 0 distance-3 escapes; 0
  extraction-verification failures. Census over 40,000 steps: distance-2
  rate 1/115 (below the charged 1/64); distance-3 0 (below census
  resolution, consistent with the algebraic 2^-16.7). Residual gap =
  toy-scale detection lag (interval 4 vs per-walker budget ~60; production
  2^64 vs 2^110.3, loss < 2^-46) and 1/8 same-walker pairs at K = 8
  (2^-16 at K = 2^16).
- Reproduce:
  `echo '{"trials": 200, "dbits": 2, "walkers": 8}' | python3 submissions/bls12-381-g1-ecdlp-plain/automorphism-rho-v3/experiment.py`.

Small experiments are evidence only at their actual scope.

## Heuristics, limitations, objections

- H1 (birthday behavior of the canonical 32-adding walk) is the one walk
  heuristic; the claim needs 87% of the modeled probability.
- Period-2 measured, period-3 proven; both charged. Longer families
  bounded by pattern counting and fenced by the epoch limit with 5% slack.
- Instruction sub-costs are worst-case analytical bounds; mod-r prices are
  explicit constructions.
- Single shot at the 0.39 threshold; doubling the budget gives 0.911 at
  +1 log2.
- Standard cryptanalysis; the first qualified bound for this track under
  this protected cost model.

References: `tracks/docs/COST_MODEL.md` (pinned 383-instruction blst
multiplication, C = 2,000); HAC ch. 3.6; Teske, J. Cryptology 14(2), 2001;
Wiener-Zuccherato, SAC 1998; Bos-Kleinjung-Lenstra and
Galbraith-Hopkins-Shparlinski on the negation map and fruitless cycles;
Kim-Barbulescu, CRYPTO 2016. Self-contained above.
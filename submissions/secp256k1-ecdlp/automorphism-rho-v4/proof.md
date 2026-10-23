# secp256k1-ecdlp: automorphism-aware Pollard rho, v4

Status: ready for independent exploratory review. Classical single-shot
attack, no quantum steps, no uncharged advice, no amortization over targets.
v4 responds to the third-round review of v3, whose two objections were both
correct and are fixed structurally:

1. **Epoch/detection conflict.** v3 combined an epoch limit E = 2^36 with
   distinguished points of density 2^-64: a merged pair survives at most
   E shared steps, so its reporting probability was E * 2^-64 = 2^-28,
   not the claimed 1 - exp(-256). v4 makes the two scales compatible:
   distinguished points now have density 2^-32 (expected reporting
   interval ~2^32 steps) and the epoch is E = 2^44, larger than the
   reporting interval by 2^12. A merged pair now reports before reset
   with probability >= 1 - 2^-11 (exact averaging of the geometric
   reporting over the reset phase), which is charged by subtraction from
   the success probability.
2. **Cycle-rate contradiction.** v3's length-L bound (1/6)*32^{-(L-1)}
   already gives ~2^-17.6 at L = 4 - it cannot establish a < 2^-40 total,
   and its epoch arithmetic E*rho/2 was inconsistent. v4 enumerates the
   algebraic cycle families exactly and redesigns the detector:
   - The zero-sum product-1 tuples over the six-multiplier group G6 are
     enumerated exactly for every length L in 2..8: N_L = 1, 2, 15, 60,
     340, 1680, 9135. A length-L same-index cycle with multiplier tuple
     in this set returns the full state identically (net map 1, net
     constant 0 - the order-3 identity 1 + lam + lam^2 = 0 makes these
     families systematic, not random coincidences; v3's union bound
     over distinct constants was the wrong model). Entry rates:
     L=2: 32^-1*6^-1 = 2^-7.58 (the classical negation cycle; the last
     multiplier is forced, so only one 1/6 factor); L>=3: N_L *
     32^{-(L-1)} * 6^{-L}, giving 2^-16.75, 2^-21.43, 2^-27.0, 2^-32.1,
     2^-37.4, 2^-42.5.
   - The look-back detector now compares the current canonical x against
     the last SIXTEEN stored states, detecting every cycle of length
     2..16 at its first return, deterministically - no rate assumption
     is needed for detection. Waste per detected entry <= L+1 steps
     (charged at the rate bounds above).
   - Undetected lengths 9..16 do not exist (detected); lengths 17 and
     beyond are bounded by the trivial tuple count N_L <= 6^L, giving
     entry rate <= 32^{-(L-1)}: summing over L >= 17 gives <= 2^-80 per
     step, and the epoch fences each such trap at <= E = 2^44 steps:
     trapped-waste fraction <= 2^-80 * 2^44 = 2^-36.
3. **Waste accounting fixed.** v3 double-charged the escape steps (once in
   the per-step ledger as amortized instructions, once implicitly in the
   useful-visit fraction). v4 charges escape steps exactly once: the walk
   budget t counts ALL executed steps (useful, trapped and escape), so
   W = per_step * t, and the useful-visit fraction rho carries only the
   lost-visit effect of wasted steps.

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) | **127.622849** (= 127622849/1000000, rounded up) |
| success_probability (supported lower bound) | 0.39 (threshold is 0.39) |
| memory_log2_bytes | 102.0 |
| preprocessing_log2 (included in total) | 25.6 |
| nonuniform_advice_log2_bytes | 0 |

Repository baseline control: 129.000802 with 2^136.000001 bytes. This claim
is lower by 1.377953 log2 (~2.59x less work) with ~2^34 times less memory.

## Conventions

Protected target: secp256k1, y^2 = x^3 + 7 over Fp, generator P of prime
order r, challenge Q = [x]P, x uniform in [0, r). Cost model
`classical-group-work-v2`: cost = H + R/1000; H counts reference operations
(batched incomplete affine secp256k1 addition, 6M with S = M, A = 0; one
reference operation = 1 unit), R counts other scalar 64-bit AArch64
instruction equivalents. Pinned field multiplication: 179 instructions
(modulo the 256-bit field prime p). Mod-r prices (r also 256 bits) are
analytic constructions: addition <= 32 instructions (4x64-limb add with
carry, compare, conditional correction); multiplication <= 384
instructions (5x52-limb schoolbook with 64x64->128 products: 25 limb
products plus accumulation, then Montgomery reduction of the same shape)
- used only once per collision; inversion <= 2*384 + 512 instructions
(once). A reference operation and its implementation are never
double-counted.

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
The order-3 identity 1 + lam + lam^2 = 0 mod r is verified in Part A and is
the engine of the period-3-and-longer cycle families handled below.

## Algorithm v4 (canonical-state walk, lazy representation, 16-state look-back)

State: the canonical point C plus a lazy representation (a, b, L), L
indexing G6, with C = [g6[L]a]P + [g6[L]b]Q. Setup (charged): zeta, lam;
32 canonical constants A_j = [c_j]P drawn pairwise class-distinct with no
A_j' = +-zeta^i A_j pair (charged setup check at production scale in Part
A); the 6x32 table g6[L]^{-1}c_j and the 6x6 table; K = 2^16 walkers
seeded at random canonical states (2^17 scalar multiplications, <= 384
reference additions each).

Step (lockstep rounds; additions batch through Montgomery inversion as in
the reference 6M convention):
1. j = x_c mod 32.
2. W = C + A_j (one reference addition); a <- a + g6[L]^{-1}c_j mod r
   (one mod-r addition); b unchanged.
3. Canonicalize W: x_W, zeta*x_W, zeta^2*x_W (2 field multiplications),
   argmin and y-parity; C' = +-zeta^i W with multiplier kappa = +-lam^i;
   L <- L * kappa (6x6 table). Invariant (verified along 96 production
   steps in Part A): rep(C') = kappa(rep(C) + (c_j, 0)) = g6[L'](a', b).
4. Look-back: keep the last sixteen full states; compare x_c against each
   stored x_c (15 four-limb compares). On a match, compare the full state:
   identical -> fruitless cycle of some length 2..16 -> escape with one
   charged addition of A_{(j+1) mod 32}; different -> useful
   self-collision -> extract (step 6).
5. Distinguished point: low 32 bits of x_{C'} zero -> insert
   (x_{C'}, a, b, L) into a hash table of fixed capacity 2^95 slots; on a
   key match with a different tuple -> extract; identical tuple -> reseed
   the walker (charged). Expected reporting interval 2^32 steps; expected
   inserts t*2^-32 = 2^94.9, below the 2^95 capacity with margin 1.07;
   overflow aborts (Chernoff bound below).
6. Extraction (once): materialize g6[L_1](a_1, b_1), g6[L_2](a_2, b_2)
   (two mod-r multiplications);
   x = (g6[L_1]a_1 - g6[L_2]a_2)(g6[L_2]b_2 - g6[L_1]b_1)^{-1} mod r
   (one inversion); verify [x]P = Q; output.
7. Epoch limit: a walker that has taken E = 2^44 steps since its last
   seed reseeds (charged; 2^9.6 units per reseed, total
   t/2^44 * 768 = 2^91.6 units). E exceeds the expected reporting interval
   2^32 by 2^12, so a merged pair reports before reset with probability
   >= 1 - 2^-11 (exact computation below); the loss is subtracted from
   the success probability.
8. Stop at the first verified extraction or after t total walker steps.
   Single shot; everything charged.

Exceptional additions occur with probability < 2^-126 over the run;
affected walkers reseed, charged.

## Fruitless cycles: exact enumeration and waste

A length-L same-index cycle with canonicalization multipliers
(m_1..m_L) returns the full state identically iff the net map
m_1...m_L = 1 and the net constant coefficient
K = sum_i (m_1...m_{i-1}) = 0 mod r (the constant then cancels by the
order-3 identity; no relation between independently drawn constants is
involved - this corrects v2/v3's invalid union-bound model). Exact
enumeration over all 6^L tuples (Part A of the experiment enumerates
L = 2..8 in seconds):

| L | N_L (zero-sum, product-1) | entry rate per step |
| --- | ---: | ---: |
| 2 | 1 | 32^-1 * 6^-1 = 2^-7.58 |
| 3 | 2 | 2 * 32^-2 * 6^-3 = 2^-16.75 |
| 4 | 15 | 15 * 32^-3 * 6^-4 = 2^-21.43 |
| 5 | 60 | 2^-27.02 |
| 6 | 340 | 2^-32.10 |
| 7 | 1680 | 2^-37.38 |
| 8 | 9135 | 2^-42.52 |

(L = 2: the tuple is (-1, -1) and the final multiplier is forced by the
return equation, so only one 1/6 factor; L >= 3: the final multiplier is a
genuine 1/6 event.) The 16-state look-back detects every cycle of length
2..16 at its first return - deterministically, no rate assumption.
Waste per detected entry <= L + 1 steps. Charged waste fraction
(L <= 8, exact rates, L=2 at the conservative 1/64): 3/64 + negligible =
0.0469. Lengths 9..16 are also detected; their entry rates are bounded by
the trivial tuple count N_L <= 6^L, giving total <= sum 32^{-(L-1)} =
2^-40.2, waste <= 17 * 2^-40.2 = 2^-36. Lengths >= 17 are not detected;
the trivial bound gives entry rate <= sum_{L>=17} 32^{-(L-1)} = 2^-80, and
the epoch fences each trap at <= E = 2^44 steps: trapped fraction
<= 2^-80 * 2^44 = 2^-36. Total waste <= 0.0469 + 2^-36 + 2^-36 = 0.0469.

## Success probability

Budget t = 1.15 * sqrt(N) = 2^126.9092 total walker steps; useful-visit
fraction rho = 1 - 0.0469 = 0.9531. Birthday ideal
1 - exp(-(1.15 * 0.9531)^2 / 2) = 0.4516. Reporting: a merge is reported
iff the merged pair hits a distinguished point before the earlier of the
two walkers' epoch resets; averaging the geometric reporting (per-step
probability 2^-32) over the uniform reset phase gives failure probability
<= 2 * 2^32 / E = 2^-11, subtracted. Table overflow: mean inserts
2^94.9 against 2^95 slots (ratio 1.07), P(overflow) <= exp(-2^89)
(Chernoff), subtracted. Walk efficiency 0.92 (Teske's r-adding
experiments; the toy measurement below). Supported bound:
0.4516 * 0.92 * (1 - 2^-11 - 2^-89) = 0.415 >= 0.39. Claimed: **0.39**.

## Work ledger (worst case every step)

| Item | Instructions | Units |
| --- | ---: | ---: |
| Reference operation: W = C + A_j (batched incomplete affine 6M addition) | reference | 1.000 |
| Canonicalization: 2 field multiplications (179 each) | 358 | 0.358 |
| Canonicalization: argmin/select, parity | 72 | 0.072 |
| Representation addition a <- a + g6[L]^{-1}c_j mod r | 32 | 0.032 |
| Multiplier/constant table updates + walk index | 16 | 0.016 |
| 16-state look-back (15 four-limb compares + branches) | 90 | 0.090 |
| Distinguished-point predicate (OR-reduce 32 bits) | 3 | 0.003 |
| State loads/stores, loop | 69 | 0.069 |
| **Per-step worst-case total** (1000 reference + 640 other) | **1640** | **1.640** |

Escape, trapped and reseed steps are walk steps inside t and are charged
exactly once here; their visit-loss effect is in rho. One-offs: setup
<= 2^14; seeding 2^17 scalar multiplications = 2^25.585 units;
distinguished-point processing <= t/2^32 = 2^94.9 operations x 1000
instructions = 2^94.9 units; epoch reseeds <= t/2^44 x 768 = 2^91.6 units;
detection lag allowance K * 2^65 steps = 2^81 units; endgame <= 700 units.
Total:

W = 1.640 * t + one-offs = 1.640 * 1.15 * sqrt((r-1)/6) + < 2^95
  = 2^127.622849 (rounded up; exact rational 127622849/1000000).

Preprocessing (roots, constants, seeding) = 2^25.585 units ->
preprocessing_log2 = **25.6**, included.

## Memory

Distinguished-point table: fixed capacity 2^95 slots x 128 bytes
(x_c 32, a 32, b 32, L/metadata 32) = 2^102 bytes; sixteen stored states
per walker (2^16 x 1.5 KB = 2^30.4), code, constants, scratch < 2^31.
Peak storage 2^102.0 bytes -> memory_log2_bytes = **102.0**. Hard
allocation; overflow aborts (bound above). Nonuniform advice:
**0 bytes**.

## Evidence and reproduction

- Part A (production, exact): automorphism identities, the order-3
  identity 1 + lam + lam^2 = 0 mod r, G6 closure, class count, the
  pairwise class-disjoint constant check, the lazy-representation
  invariant along 96 production walk steps (96/96), and the exact
  zero-sum tuple enumeration N_L for L = 2..8 (1, 2, 15, 60, 340, 1680,
  9135).
- Part B (toy, complete v4 attack, production low-bits predicate):
  declared run `{"trials": 200, "dbits": 2, "walkers": 8}`: 85/200 =
  0.425 vs ideal 0.4845; 335 escapes; 0 extraction-verification
  failures. Cycle census over 40,000 single-walk steps with
  per-executed-step denominators: 347 full-state returns, rate
  0.00868 per executed step (the L=2 family; below the charged 1/64),
  no returns at distances 3..16 above census resolution (consistent
  with the enumerated rates). The residual success gap is toy-scale
  detection lag (distinguished interval 4 vs per-walker budget ~60;
  production 2^32 vs 2^110.9, relative loss < 2^-46) and 1/8
  same-walker birthday pairs at K = 8 (2^-16 at K = 2^16).
- Reproduce:
  `echo '{"trials": 200, "dbits": 2, "walkers": 8}' | python3 submissions/secp256k1-ecdlp/automorphism-rho-v4/experiment.py`.

Small experiments are evidence only at their actual scope.

## Heuristics, limitations, objections

- H1 (birthday behavior of the canonical 32-adding walk) is the one walk
  heuristic; supported classically, by Teske's experiments, and by the
  toy (0.545 vs 0.4845 ideal in the lag-free dbits=0 regime). The claim
  needs 92% of the modeled probability.
- Cycle families: exact enumeration for L <= 8, deterministic detection
  for L <= 16, trivial-bound + epoch fence for L >= 17. No union bound
  over distinct constants remains anywhere.
- The epoch (2^44) and the reporting interval (2^32) are now separated by
  2^12; the missed-report loss (2^-11) is computed by exact averaging and
  subtracted.
- Instruction sub-costs are worst-case analytical bounds; mod-r prices
  are explicit constructions.
- Single shot at the 0.39 threshold; doubling the budget gives 0.911 at
  +1 log2.
- Standard cryptanalysis; the first qualified bound for this track under
  this protected cost model.

References: `tracks/docs/COST_MODEL.md`; HAC ch. 3.6; Teske, J. Cryptology
14(2), 2001; Wiener-Zuccherato, SAC 1998; van Oorschot-Wiener on
distinguished-point collision search. Self-contained above.

# bls12-381-g1-ecdlp-plain: automorphism-aware Pollard rho, v2

Status: ready for independent exploratory review. Classical single-shot
attack, no quantum steps, no uncharged advice, no amortization over targets.
This v2 replaces the rejected `automorphism-rho` package and responds to the
recorded review objections in the same way as the sibling secp256k1-ecdlp v2
package (different curve, prices and parameters; the constructions are
identical in structure):

- The v1 fruitless-cycle analysis was wrong (it claimed an O(k^2 |G|/N)
  trigger rate). v2 derives the correct constant-order rate for the
  negation-map 2-cycle, detects it with a 2-step look-back on the full walk
  state, escapes with one charged addition, and charges the waste.
- The v1 toy experiment used `x mod dp == 0` (dp=1 stored every point) and
  so never exercised the production distinguished-point regime. v2 uses the
  production predicate (low bits zero) and measures it.
- The v1 ledger omitted the representation addition, mod-r operation prices,
  state/loop work, and hard memory/lag bounds (and contained an internal
  inconsistency: the claim text said 839 non-reference instructions while
  its own table summed 1,655). v2 eliminates per-step mod-r multiplications
  entirely (lazy six-value multiplier bookkeeping), completes the ledger,
  fixes the preprocessing arithmetic, and replaces probabilistic memory/lag
  statements with caps plus explicit tail bounds.

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) | **126.976677** (= 126976677/1000000, rounded up) |
| success_probability (supported lower bound) | 0.39 (threshold is 0.39) |
| memory_log2_bytes | 70.0 |
| preprocessing_log2 (included in total) | 25.6 |
| nonuniform_advice_log2_bytes | 0 |

No qualified baseline exists for this track (its selected package is an
unscored template). For scale on the same curve: vanilla rho at the 0.39
birthday budget is ~2^128.0, negation-only rho ~2^128.1, and this
order-6-automorphism rho at the same accounting is the submitted 2^126.98.

## Conventions

Protected target: BLS12-381 G1, y^2 = x^3 + 4 over Fp, generator P0 of prime
order r, challenge T = [x]P0, x uniform in [0, r) (T = infinity exactly when
x = 0; answered immediately). Cost model `classical-group-work-v2`:
cost = H + R/2000; H counts reference operations (batched incomplete affine
BLS12-381 G1 addition, 6M with S = M, A = 0; one reference operation = 1
unit), R counts other scalar 64-bit AArch64 instruction equivalents. Pinned
base-field multiplication price: 383 instructions (blst `blst_fp_mul`;
applies to multiplication modulo the 381-bit field prime p; the G1 reference
price C = 2,000 is the rounded 6 x 383 = 2,298).

Prices for arithmetic modulo the scalar order r (255 bits), which the pinned
measurement does not cover, are supplied analytically: mod-r
addition/subtraction (4x64-limb add with carry, compare, conditional
correction): <= 32 instructions; mod-r multiplication (5x52-limb schoolbook
with 64x64->128 products plus Montgomery reduction): <= 384 instructions
(used only for the one-off collision extraction and inversion; the walk
performs no per-step mod-r multiplication); mod-r inversion by Fermat:
<= 2*384 + 512 instructions, used once. A reference operation and its
implementation are never double-counted.

## The order-6 automorphism group

BLS12-381 has j-invariant 0 and p mod 3 = 1, so Fp contains a primitive cube
root of unity zeta with zeta^3 = 1, zeta != 1, and (x, y) -> (zeta*x, y) is
a group automorphism fixing infinity. It maps G1 = <P0> to itself and acts
there as multiplication by lam with lam^3 = 1 mod r, lam != 1, fixed by the
verified identity [lam]P0 = zeta(P0) (experiment, Part A) with

lam = 0x73eda753299d7d483339d80809a1d804a7780001fffcb7fcfffffffe00000001.

With negation, G = {+-1, +-lam, +-lam^2} (order 6, closed under
multiplication mod r - verified) acts freely on the r - 1 nonzero points
(r mod 6 = 1, and none of -1, lam, lam^2, -lam, -lam^2 equals 1), giving

N = (r - 1)/6 = 8739312529187698413241290084697660972948425416754606303767276449989763530752

classes with canonical representative (x_c, y_c),
x_c = min(x, zeta*x, zeta^2*x), y_c even.

## Algorithm v2 (canonical-state walk, lazy representation)

Identical in structure to the sibling v2 package; the state is the canonical
point C plus a lazy representation (a, b, L), L indexing
G6 = {1, lam, lam^2, -1, -lam, -lam^2}, with

  C = [g6[L]*a mod r]P0 + [g6[L]*b mod r]T.

Setup (charged): zeta, lam derivation; 32 canonical walk constants
(pairwise class-distinct, no A_j' = +-zeta^i A_j pair - a charged setup
check performed at production scale in Part A); the 6x32 table
g6[L]^{-1} c_j mod r and the 6x6 multiplier table; K = 2^16 walkers seeded
at random canonical states.

Step (lockstep rounds; additions batch through Montgomery inversion as in
the reference 6M convention):
1. j = x_c mod 32.
2. W = C + A_j (one reference addition); a <- a + g6[L]^{-1} c_j mod r
   (one mod-r addition); b unchanged.
3. Canonicalize W: compute x_W, zeta*x_W, zeta^2*x_W (2 base-field
   multiplications), argmin and y-parity; C' = +-zeta^i W with multiplier
   kappa = +-lam^i; L <- L * kappa via the 6x6 table. Invariant:
   rep(C') = kappa*(rep(C) + (c_j,0)) = g6[L']*(a', b). Verified along 96
   production walk steps in Part A.
4. 2-step look-back detector on the full state (x_c, a, b, L): identical
   two steps back -> fruitless cycle -> escape with one charged addition of
   A_{(j+1) mod 32}; different -> useful self-collision -> extract.
5. Distinguished point: low 64 bits of x_{C'} zero -> insert
   (x_{C'}, a, b, L) into a hash table of fixed capacity 2^63 slots; on a
   key match with a different tuple, extract; identical tuple -> reseed.
6. Extraction (once): materialize g6[L_1](a_1, b_1), g6[L_2](a_2, b_2)
   (two mod-r multiplications),
   x = (g6[L_1]a_1 - g6[L_2]a_2) * (g6[L_2]b_2 - g6[L_1]b_1)^{-1} mod r,
   verify [x]P0 = T, output.
7. Stop at the first verified extraction or after t total walker steps;
   single shot, everything charged.

Exceptional additions occur with probability <= t * 32 * 2 / 2^256 < 2^-126
over the run; the affected walker reseeds, charged.

## Fruitless cycles (corrected analysis)

The negation-map 2-cycle: at canonical state C with index j, if
y(C + A_j) is odd the canonicalization picks D = -(C + A_j); if also
x_D mod 32 = j, the next step returns D + A_j = -C, whose canonical image is
C, and the lazy state returns identically. Trigger probability about
(1/2)(1/3)(1/32) ~ 1/192 per step - constant-order, not vanishing with N
(v1 claimed otherwise and was wrong). Measured toy escape rates: 1/183 to
1/287 per step across regimes, inside the band. Handling: the 2-step
look-back fires within 2 steps of entry, the full-state comparison separates
fruitless returns from useful self-collisions, and the escape is one charged
addition; waste <= 3 steps per event at rate <= 1/64: <= 3 * 1.4875 / 64
= 0.0697 units per step, charged in the ledger.

Longer fruitless cycles (L >= 3): a returning cycle has net class-map
+-zeta^i and net constant K (sum of <= L signed zeta-images of constants).
If +-zeta^i != 1, the return equation (1 -+ zeta^i) C = K pins C; for
L <= 8 there are at most (32*6)^8 = 2^49 (map, K) pairs, entry probability
< 2^-203 per step, < 2^-75 per run. If +-zeta^i = 1, the return requires
K = 0: 2-term constant relations are excluded deterministically by the
setup check; 3-to-8-term relations have probability < (32*12)^8 / r < 2^-165
over the draw. Longer mechanisms have no constant-order structure
(Bos-Kleinjung-Lenstra, Galbraith-Hopkins-Shparlinski: the 2-cycle
dominates), and the toy at N = 2^17.4 shows escape counts fully accounted
for by the 2-cycle rate. A residual 5% efficiency slack is charged below.

## Success probability

Budget t = 1.15 * sqrt(N) = 2^126.3376 total walker steps. With
rho >= 0.95 (2-cycle waste plus the charged 5% residual slack),

P(visit collision by t) >= 1 - exp(-(1.15*0.95)^2/2) = 0.450.

Detection: collisions surface at the next distinguished point (density
2^-64); detection lag consumes <= K * 2^64 / t < 2^-46 of the budget,
subtracted. The table has a hard capacity of 2^63 slots; distinguished-point
count within budget is binomial with mean t/2^64 = 2^62.34, and
P(overflow) < exp(-2^61.5) by a Chernoff bound; on overflow the run fails,
subtracting < 2^-61. Supported lower bound:
0.450 x 0.92 (Teske's r-adding experiments; the toy measurement) - 2^-46
= 0.414 >= 0.39. Claimed: **0.39**.

## Work ledger (worst case every step)

| Item | Instructions | Units |
| --- | ---: | ---: |
| Reference operation: W = C + A_j, the calibrated batched incomplete affine 6M G1 addition (3M formula + 3M batch-inversion amortisation across the 2^16 lockstep walkers) | reference | 1.000 |
| Canonicalization: x_W*zeta, x_W*zeta^2 (2 base-field mults, 383 each) | 766 | 0.383 |
| Canonicalization: 3-way argmin over 6-limb values, select | 64 | 0.032 |
| Canonicalization: y-parity test | 8 | 0.004 |
| Representation update: a <- a + g6[L]^{-1}c_j mod r (one 255-bit mod-r addition) | 32 | 0.016 |
| Multiplier table update + walk index | 16 | 0.008 |
| 2-step look-back compare | 20 | 0.010 |
| Distinguished-point predicate | 5 | 0.0025 |
| State loads/stores, loop and addressing work | 64 | 0.032 |
| Escape amortization: <= 3 steps + 1 escape addition per fruitless event at rate <= 1/64 | <= 139 | 0.0697 |
| **Per-step worst-case total** (2000 reference-equivalent + 975 other instructions) | **2975** | **1.4875 + 0.0697 = 1.5572** |

One-off components: zeta/lam derivation <= 2^10; constants + tables + setup
checks <= 2^14; walker seeding 2K = 2^17 scalar multiplications x <= 384
reference additions each = 2^25.585 (v1's 2^25.32 was an arithmetic slip);
distinguished-point processing <= t/2^64 = 2^62.34 operations x <= 1000
instructions = 2^62.3 units; endgame <= 700 units. Total:

W = 1.5572 * t + one-offs = 1.5572 * 1.15 * sqrt((r-1)/6) + < 2^63
  = 2^126.976677 (rounded up; exact rational 126976677/1000000).

Preprocessing (roots, constants, seeding) = 2^25.585 units ->
preprocessing_log2 = **25.6**, included in the total.

## Memory

Distinguished-point table: fixed capacity 2^63 slots x 128 bytes = 2^70
bytes, plus walker state (2^16 x 256 bytes), code, constants, scratch
< 2^24 bytes. Peak storage 2^70.0 bytes -> memory_log2_bytes = **70.0**.
The capacity is a hard allocation; overflow aborts (tail bound above). All
tables are constructed and charged at runtime; nonuniform advice is
**0 bytes**.

## Why not the pairing / GT transfer

Unchanged from v1 and declined: generic rho in GT costs ~10x per step (the
dense Fp12 reference price is 20,000), and index calculus over Fp12 would
need a complete analytic NFS/exTNFS cost model; BLS12-381's 381-bit modulus
was sized (Kim-Barbulescu) to keep that route near the 2^128 class. The
generic bound above is the honest first qualified score for this track.

## Evidence and reproduction

- `experiment.py` Part A (production parameters, exact): p mod 3 = 1,
  r mod 6 = 1, zeta^3 = 1 with curve preservation, [lam]P0 = zeta(P0) and
  [lam^2]P0 = zeta^2(P0), lam^3 = 1 mod r, G6 closed under multiplication,
  class count (r-1)/6, the pairwise class-disjoint constant check at
  production scale, and the lazy-representation invariant verified along
  96 production walk steps (96/96).
- Part B (toy, complete v2 attack, 20-bit j=0 curve, N_toy = 174563):
  declared run `{"trials": 200, "dbits": 2, "walkers": 8}`: 85/200 = 0.425
  successes vs ideal 0.4845; 335 escapes (rate 1/287 per step, inside the
  1/64-1/192 band); 0 extraction-verification failures. The residual gap is
  toy-scale detection lag (distinguished interval 4 walker-steps vs a
  per-walker budget of ~60 steps, and 1/8 same-walker birthday pairs at
  K = 8), while at production the interval 2^64 stands against a
  per-walker budget of 2^110.3 and same-walker pairs are 2^-16 (relative
  losses < 2^-46). Local sweep (200 trials): dbits=0: 0.545 vs 0.4845
  ideal; dbits=2: 0.425; dbits=4: 0.315.
- Reproduce:
  `bash scripts/run_experiments.sh bls12-381-g1-ecdlp-plain submissions/bls12-381-g1-ecdlp-plain/automorphism-rho-v2`
  or directly:
  `echo '{"trials": 200, "dbits": 2, "walkers": 8}' | python3 submissions/bls12-381-g1-ecdlp-plain/automorphism-rho-v2/experiment.py`.

Small experiments are evidence only at their actual scope.

## Heuristics, limitations, objections

- H1: the canonical 32-adding walk's first class-collision time follows the
  birthday tail 1 - exp(-t^2/2N). Supported by classical rho analysis,
  Teske's r-adding experiments, and the toy run (0.545 measured vs 0.4845
  ideal in the lag-free dbits=0 regime). The claimed 0.39 needs only 87% of
  the modeled 0.450.
- The fruitless 2-cycle rate is constant-order; it is detected, escaped and
  charged. Longer fruitless cycles: deterministic exclusion (2-term),
  union bounds (3-8 terms), empirical support (toy; BKL; GHS), and a charged
  5% residual slack.
- Instruction sub-costs are worst-case analytical bounds; mod-r prices are
  explicit constructions (32/384 instructions), not transplanted mod-p
  measurements.
- Single shot at the 0.39 threshold; doubling the budget gives 0.911 at
  +1 log2. Memory 2^70.0 bytes, reported, no score effect.
- Standard cryptanalysis, not a novel attack; the first qualified bound for
  this track under this protected cost model.

References: `tracks/docs/COST_MODEL.md` (pinned 383-instruction blst
multiplication, C = 2,000); HAC ch. 3.6; Teske, J. Cryptology 14(2), 2001;
Wiener-Zuccherato, SAC 1998; Bos-Kleinjung-Lenstra and
Galbraith-Hopkins-Shparlinski on the negation map and fruitless cycles;
Kim-Barbulescu, CRYPTO 2016. Self-contained above.

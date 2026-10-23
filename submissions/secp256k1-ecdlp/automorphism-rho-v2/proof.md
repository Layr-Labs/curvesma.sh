# secp256k1-ecdlp: automorphism-aware Pollard rho, v2

Status: ready for independent exploratory review. Classical single-shot
attack, no quantum steps, no uncharged advice, no amortization over targets.
This v2 replaces the rejected `automorphism-rho` package and responds to the
recorded review objections:

- The v1 fruitless-cycle analysis was wrong (it claimed an O(k^2 |G|/N)
  trigger rate). v2 derives the correct constant-order rate for the
  negation-map 2-cycle, detects it with a 2-step look-back on the full walk
  state, escapes with one charged addition, and charges the waste.
- The v1 toy experiment used `x mod dp == 0` (dp=1 stored every point) and
  so never exercised the production distinguished-point regime. v2 uses the
  production predicate (low bits zero) and measures it.
- The v1 ledger omitted the representation addition, mod-r operation prices,
  state/loop work, and hard memory/lag bounds. v2 eliminates per-step mod-r
  multiplications entirely (lazy six-value multiplier bookkeeping), adds the
  missing items, fixes the preprocessing arithmetic, and replaces the
  probabilistic memory/lag statements with caps plus explicit tail bounds.

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) | **127.623728** (= 7976483/62500, rounded up) |
| success_probability (supported lower bound) | 0.39 (threshold is 0.39) |
| memory_log2_bytes | 70.0 |
| preprocessing_log2 (included in total) | 25.6 |
| nonuniform_advice_log2_bytes | 0 |

Repository baseline control: 129.000802 with 2^136.000001 bytes. This claim
is lower by 1.377074 log2 (~2.60x less work) with ~2^65 times less memory.

## Conventions

Protected target: secp256k1, y^2 = x^3 + 7 over Fp, generator P of prime
order r, challenge Q = [x]P, x uniform in [0, r). Cost model
`classical-group-work-v2`: cost = H + R/1000; H counts reference operations
(batched incomplete affine secp256k1 addition, 6M with S = M, A = 0; one
reference operation = 1 unit), R counts other scalar 64-bit AArch64
instruction equivalents. Pinned field-multiplication price: 179 instructions
(the organizer's Bitcoin-Core 5x52 measurement; applies to multiplication
modulo the 256-bit field prime p).

Prices for arithmetic modulo the scalar order r (also 256 bits), which the
pinned measurement does not directly cover, are supplied analytically:
- mod-r addition/subtraction of canonical 256-bit values (4x64-limb add with
  carry, compare against r, one conditional correction): <= 32 instructions.
- mod-r multiplication (5x52-limb schoolbook with 64x64->128 products:
  25 limb products plus accumulation; Montgomery reduction with the same
  shape): <= 384 instructions. This bound is used only for the one-off
  collision extraction and inversion; the per-step walk below performs no
  mod-r multiplication at all.
- mod-r inversion by Fermat (<= 2*256 multiplications plus loop work):
  <= 2*384 + 512 instructions, used once.
A reference operation and its implementation are never double-counted.

## The order-6 automorphism group

secp256k1 has j-invariant 0 and p mod 3 = 1, so Fp contains a primitive cube
root of unity zeta, zeta^3 = 1, zeta != 1, and (x, y) -> (zeta*x, y) is a
group automorphism fixing infinity. On <P> it acts as multiplication by lam
with lam^3 = 1 mod r, lam != 1; the correct element of {lam, lam^2} is fixed
by the verified identity [lam]P = zeta(P) (experiment, Part A). With
negation, G = {+-1, +-lam, +-lam^2} (order 6, closed under multiplication
mod r - verified) acts freely on the r - 1 nonzero points (r mod 6 = 1, and
none of -1, lam, lam^2, -lam, -lam^2 equals 1), giving

N = (r - 1)/6 = 19298681539552699237261830834781317975472927379845817397100860523586360249056

classes with canonical representative (x_c, y_c),
x_c = min(x, zeta*x, zeta^2*x), y_c even.

## Algorithm v2 (canonical-state walk, lazy representation)

The walk state is the canonical point C together with a lazy representation
(a, b, L): an index L into G6 = {1, lam, lam^2, -1, -lam, -lam^2} such that

  C = [g6[L]*a mod r]P + [g6[L]*b mod r]Q.

Setup (charged):
1. zeta, lam as above (few field/scalar exponentiations).
2. Walk constants: draw u_j uniform in [1, r), j = 0..31; canonicalize
   A_j = [u_j]P to (C_j, c_j) with C_j = [c_j]P canonical. Setup check:
   reject the draw unless the 32 canonical classes are pairwise distinct
   and no pair satisfies A_j' = +-zeta^i A_j (32^2 * 6 point comparisons);
   this deterministically excludes the constant-constellation cycles below.
   Precompute the 6x32 table g6[L]^{-1} c_j mod r and the 6x6 multiplier
   table. (Experiment Part A performs this check at production scale.)
3. K = 2^16 walkers seeded at random canonical states (see ledger for the
   scalar-multiplication charge).

Step (lockstep rounds over all walkers; the additions batch through
Montgomery inversion exactly as in the reference 6M convention):
4. j = x_c mod 32.
5. W = C + A_j (one reference addition). Track a <- a + g6[L]^{-1} c_j mod r
   (one mod-r addition, table lookup); b unchanged.
6. Canonicalize W: compute x_W, zeta*x_W, zeta^2*x_W (2 field
   multiplications), take the argmin and the y-parity; the new canonical
   point C' = +-zeta^i W, and the canonicalization multiplier is
   kappa = +-lam^i. Update L <- L * kappa via the 6x6 table. The lazy
   invariant is preserved:
   rep(C') = kappa*(rep(C) + (c_j, 0)) = kappa*g6[L]*(a + g6[L]^{-1}c_j, b)
           = g6[L']*(a', b).
   Part A of the experiment verifies this invariant along 96 production
   steps.
7. 2-step look-back detector: keep the full state (x_c, a, b, L) from two
   steps back; compare x_c (two 4-limb compares). On a match, compare the
   remaining fields: identical -> fruitless cycle -> escape (one extra
   addition with constant A_{(j+1) mod 32}, same lazy machinery, charged);
   different -> useful self-collision -> extract as in step 9.
8. Distinguished point: if the low 64 bits of x_{C'} are zero, attempt to
   insert (x_{C'}, a, b, L) into a hash table of fixed capacity
   2^63 slots (below). If x_{C'} is present with a different (a, b, L),
   extract; if present with the identical tuple, reseed the walker (charged).
9. Extraction (once per run): materialize the canonical-frame
   representations g6[L_1](a_1, b_1) and g6[L_2](a_2, b_2) (two mod-r
   multiplications), compute
   x = (g6[L_1]a_1 - g6[L_2]a_2) * (g6[L_2]b_2 - g6[L_1]b_1)^{-1} mod r
   (one mod-r inversion), verify [x]P = Q (one scalar multiplication), and
   output x.
10. The run stops at the first verified extraction or after t total walker
    steps. Single shot; every trial, escape, reseed and failure is charged.

Exceptional additions (W = infinity or incomplete-formula x-collisions)
occur with probability <= t * 32 * 2 / 2^256 over the whole run (< 2^-126);
the affected walker reseeds, charged.

## Fruitless cycles (corrected analysis)

The negation-map 2-cycle. Let the current canonical state be C with index j.
If y(C + A_j) is odd, canonicalization picks the negation image
D = -(C + A_j); if additionally x_D = x(C + A_j) satisfies x_D mod 32 = j,
then the next step computes D + A_j = -C, whose canonical image is C again,
and the lazy state (a, b, L) returns identically (deterministic walk). The
trigger is a pair of branch conditions with probability about (1/2) * (1/3)
* (1/32) ~ 1/192 per step (1/2: odd y; 1/3: argmin at i = 0; 1/32: index
match), i.e. constant-order, NOT vanishing with N. v1 claimed otherwise and
was wrong; the measured escape rate in the toy run is 1/183 to 1/287 per
step across regimes, consistent with this band.

Handling: the 2-step look-back fires at the first return (within 2 steps of
entering any period-2 cycle), the full-state comparison distinguishes
fruitless returns (identical tuple) from useful self-collisions (different
tuple, extracted), and the escape is one charged addition. Waste: <= 3
steps per event at rate <= 1/64 (conservative upper end of the band):
<= 3 * 1.641 / 64 = 0.077 units per step, charged in the ledger.

Longer fruitless cycles (L >= 3). A cycle whose state returns after L steps
has net class-map +-zeta^i and net constant K = sum of <= L signed
zeta-images of constants. If +-zeta^i != 1, the return equation
(1 -+ zeta^i) C = K pins C to a single class per (map, K) pair, and for
L <= 8 the number of pairs is at most (32 * 6)^8 = 2^49, giving entry
probability <= 2^49/N < 2^-203 per step and < 2^-75 over the full run.
If +-zeta^i = 1 (identity net map), the return requires K = 0, an additive
relation among distinct constants with signed zeta-coefficients: 2-term
relations are excluded deterministically by the setup check; relations with
3 to 8 distinct terms have probability at most (32*12)^8 / r < 2^-172 over
the draw. Longer identity-map relations and cycles of length > 8 have no
constant-order mechanism: they are the phenomenon measured by
Bos-Kleinjung-Lenstra and Galbraith-Hopkins-Shparlinski, whose experiments
show the 2-cycle dominates and longer fruitless cycles are negligible; the
toy run at N = 2^17.4 (where any L-cycle mechanism is 2^235 times more
visible relative to N) shows escape counts fully accounted for by the
2-cycle rate. A residual 5% efficiency slack is charged in the success
analysis for anything unmodelled.

## Success probability

Budget t = 1.15 * sqrt(N) = 2^126.9091 total walker steps (summed over
walkers). Under the standard random-walk model for the canonical 32-adding
walk (H1, below), the number of useful class visits is rho*t with
rho >= 0.95 (2-cycle waste <= 3/192 plus the charged 5% residual slack),
and

P(visit collision by t) >= 1 - exp(-(1.15 * 0.95)^2 / 2) = 0.450.

Detection: a collision is reported when the colliding walkers' common
successor hits a distinguished point (density 2^-64). The detection lag is
at most 2^64 walker-steps per merged pair; the fraction of budget consumed
after the last useful merge is <= K * 2^64 / t = 2^80.4 / 2^126.9 < 2^-46,
subtracted from the success probability. The distinguished-point table has
a hard capacity of 2^63 entries; the number of distinguished points within
the budget is binomial with mean t/2^64 = 2^62.9, and P(capacity overflow)
< exp(-2^62) by a Chernoff bound (factor > 2 margin); on overflow the run
fails, subtracting < 2^-62. Both tail terms are negligible against the
margin below.

Supported lower bound: 0.450 x 0.92 (walk efficiency: Teske's r-adding
experiments and the toy measurement below; the 2-cycle waste is already
inside the 0.95 factor) - 2^-46 = 0.414 >= 0.39. Claimed: **0.39**.

## Work ledger (worst case every step)

| Item | Instructions | Units |
| --- | ---: | ---: |
| Reference operation: W = C + A_j, the calibrated batched incomplete affine 6M addition (3M formula + 3M batch-inversion amortisation across the 2^16 lockstep walkers) | reference | 1.000 |
| Canonicalization: x_W*zeta, x_W*zeta^2 (2 field mults, 179 each) | 358 | 0.358 |
| Canonicalization: 3-way argmin over 4-limb values, select | 64 | 0.064 |
| Canonicalization: y-parity test | 8 | 0.008 |
| Representation update: a <- a + g6[L]^{-1}c_j mod r (one 256-bit mod-r addition) | 32 | 0.032 |
| Multiplier table update + walk index | 16 | 0.016 |
| 2-step look-back compare (x_c against the stored 2-back value; full-state compare only on x-match) | 20 | 0.020 |
| Distinguished-point predicate (OR-reduce 64 bits, branch) | 5 | 0.005 |
| State loads/stores, loop and addressing work | 64 | 0.064 |
| Escape amortization: <= 3 steps + 1 escape addition per fruitless event at rate <= 1/64 | <= 77 | 0.077 |
| **Per-step worst-case total** (1000 reference-equivalent + 641 other instructions) | **1641** | **1.641** |

One-off components: zeta/lam derivation <= 2^10; constants + tables
(32 scalar multiplications + setup checks) <= 2^14; walker seeding
2K = 2^17 scalar multiplications x <= 384 reference additions each
= 2^25.585 (this corrects the v1 figure 2^25.32, which was an arithmetic
slip); distinguished-point processing <= t/2^64 = 2^62.9 inserts/lookups x
<= 1000 instructions = 2^62.9 units; endgame extraction + verification +
output <= 700 units. Total:

W = 1.641 * t + one-offs = 1.641 * 1.15 * sqrt((r-1)/6) + < 2^63
  = 2^127.623728 (rounded up; exact rational 7976483/62500).

Preprocessing (roots, constants, seeding) = 2^25.585 units ->
preprocessing_log2 = **25.6**, included in the total.

## Memory

Distinguished-point table: fixed capacity 2^63 slots x 128 bytes (x_c 32,
a 32, b 32, L/metadata 8, chaining 24) = 2^70 bytes, plus walker state
(2^16 x 256 bytes), code, constants and scratch < 2^24 bytes. Peak storage
2^70.0 bytes -> memory_log2_bytes = **70.0**. The capacity is a hard
allocation; overflow aborts (tail bound above). All tables are constructed
and charged at runtime; nonuniform advice is **0 bytes**.

## Evidence and reproduction

- `experiment.py` Part A (production parameters, exact): p mod 3 = 1,
  r mod 6 = 1, zeta^3 = 1 with curve preservation, [lam]P = zeta(P) and
  [lam^2]P = zeta^2(P), lam^3 = 1 mod r, G6 closed under multiplication,
  class count (r-1)/6, the pairwise class-disjoint constant check at
  production scale, and the lazy-representation invariant verified along
  96 production walk steps (96/96).
- Part B (toy, complete v2 attack, 20-bit j=0 curve, N_toy = 174563):
  declared run `{"trials": 200, "dbits": 2, "walkers": 8}`: 85/200 = 0.425
  successes vs ideal 0.4845; 335 escapes (rate 1/287 per step, inside the
  1/64-1/192 band); 0 extraction-verification failures. The residual gap
  is toy-scale detection lag: the distinguished interval 4 walker-steps
  against a per-walker budget of ~60 steps, while at production the
  interval 2^64 stands against a per-walker budget of 2^110.9 (relative
  loss < 2^-46); additionally 1/8 of the toy's birthday pairs are
  same-walker at K = 8, whose detection completes only after a cycle
  traversal, versus 2^-16 at production K = 2^16. Local sweep (200 trials):
  dbits=0 (store every point): 0.545 vs 0.4845 ideal; dbits=2: 0.425;
  dbits=4: 0.315.
- Reproduce:
  `bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/automorphism-rho-v2`
  or directly:
  `echo '{"trials": 200, "dbits": 2, "walkers": 8}' | python3 submissions/secp256k1-ecdlp/automorphism-rho-v2/experiment.py`.

Small experiments are evidence only at their actual scope.

## Heuristics, limitations, objections

- H1: the canonical 32-adding walk's first class-collision time follows the
  birthday tail 1 - exp(-t^2/2N). Supported by classical rho analysis,
  Teske's r-adding experiments, and the toy run (0.545 measured vs 0.4845
  ideal in the lag-free dbits=0 regime). The claimed 0.39 needs only
  87% of the modeled 0.450.
- The fruitless 2-cycle rate (~1/192, bounded by 1/64) is constant-order;
  it is detected, escaped and charged - not argued away. Longer fruitless
  cycles are excluded deterministically (2-term), by union bounds (3-8
  terms, < 2^-172), and empirically (toy; Bos-Kleinjung-Lenstra,
  Galbraith-Hopkins-Shparlinski), with a charged 5% residual slack.
- Instruction sub-costs are worst-case analytical bounds; mod-r arithmetic
  prices are supplied by explicit constructions (32/384 instructions), not
  by transplanting the mod-p measurement.
- Single shot at the 0.39 threshold; doubling the budget gives 0.911 at
  +1 log2. Memory 2^70.0 bytes is reported and has no score effect.
- Standard cryptanalysis, not a novel attack; the first qualified bound for
  this track under this protected cost model.

References: `tracks/docs/COST_MODEL.md`; HAC ch. 3.6; Teske, J. Cryptology
14(2), 2001; Wiener-Zuccherato, SAC 1998; Bos, Kleinjung, Lenstra "On the
use of the negation map in Pollard rho" and Galbraith, Hopkins, Shparlinski
on fruitless cycles. Self-contained above; no citation needs to be fetched
to reconstruct the construction or the accounting.

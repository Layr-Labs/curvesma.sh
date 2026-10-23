# secp256k1-ecdlp: negation-only canonical Pollard rho v5

Status: ready for independent exploratory review. Classical single-shot
attack, no quantum steps, no uncharged advice, no amortization over
targets. v5 abandons the order-6 automorphism walk of v1-v4 and returns
to the NEGATION-ONLY canonical walk on the (r-1)/2 classes {C, -C}
(canonical representative: even y). The single walk heuristic (H1:
generic states have uniform independent indices x mod 32 and parities)
supports a COMPLETE cycle-family analysis for this walk, replacing the
mixed-index automorphism families that the fourth-round review
correctly identified as unbounded in v4.

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) | **128.141547** (= 128141547/1000000, rounded up) |
| success_probability (supported lower bound) | 0.39 (modeled 0.4284) |
| memory_log2_bytes | 103.0 (distinguished-point table, 2^96 slots x 128 B) |
| preprocessing_log2 (included in total) | 22.0 |
| nonuniform_advice_log2_bytes | 0 |

The repository baseline control for this track is 129.000802 with 2^136.000001 bytes; this claim is lower by 0.859255 log2 (~1.82x less work) with ~2^33 times less memory.

## Conventions

Protected target: secp256k1, y^2 = x^3 + 7 over Fp; generator of prime order r; challenge
Q = [x]P for the fixed standard generator P, x uniform in [0, r). Cost model `classical-group-work-v2`:
cost = H + R/1000; H counts reference operations (batched incomplete affine secp256k1 addition (6M, S = M, A = 0); one
reference operation = 1 unit), R counts other scalar 64-bit AArch64
instruction equivalents. Pinned field multiplication: 179 instructions (5x52-limb secp256k1 field multiplication). Mod-r
operations are analytic constructions with explicit limb counts per the
ledger below. A reference operation and its implementation are never
double-counted.

## The walk

State: canonical point C (even y) with a lazy representation (a, b)
with C = [a]P + [b]Q and a sign-corrected seeding convention (the
canonicalization multiplier of every seed is folded into the stored
representation; a mis-signed representation would make every
extraction of that walker return -x - the toy enforces the convention
by construction and records zero wrong extractions).

Constants: A_j = [2^j]P for j = 0..31 (all distinct; verified in Part
A). The powers of two make every cycle relation an EXACT INTEGER
equation: any signed sum of at most 16 constants satisfies
|sum| <= 16 * 2^31 < r (verified in Part A), so it is zero only if each
power of two appears equally often with both signs.

Step: j = x_c mod 32; W = C + A_j (one reference operation); the
canonicalization multiplier kappa is +1 if y(W) is even else -1; C' =
kappa * W (one parity test, one conditional field negation);
representation a <- kappa(a + 2^j) mod r, b <- kappa*b mod r (one
mod-r addition, one conditional mod-r negation, one sign flip).

## Fruitless cycles: a complete census

A length-L cycle family is a pair (word w_1..w_L, sign pattern
kappa_1..kappa_L) with w_1 = the current index and per-symbol even
multiplicity. The composed map is

  C_(L+1) = (prod kappa) C_1 + sum_u s_u A_(w_u),   s_u = prod_(k>=u) kappa_k,

so the cycle conditions are prod kappa = +1 and the per-symbol balance
sum of s_u = 0 - an exact integer relation among powers of two. The map
kappa <-> s is a bijection (kappa_u = s_u * s_(u+1), s_1 = 1,
s_(L+1) = 1), so the solutions per word number exactly
prod_s C(2 t_s, t_s) / 2 over symbol multiplicities t_s. In a GENUINE
family no prefix closes (C_(t+1) = C_1 for t < L), and kappa_u is
forced iff W_u = +/-C_1 iff the prefix closes - so in a genuine
family every kappa except the last is a free parity event and

  rate_L = (#genuine families of length L) * 64^-(L-1)

exactly under H1 (L-1 index events at 1/32 and L-1 parity events at
1/2; the last kappa is forced by the return). The census is EXECUTED by
the experiment at production parameters (Part A recomputes and asserts
the counts):

| L | genuine families | rate per step |
| --- | ---: | ---: |
| 2 | 1 | 1/64 = 2^-6 (the classical negation 2-cycle) |
| 4 | 157 | 2^-10.705 |
| 6 | 46,844 | 2^-14.484 |
| 8..16 | conservative bound (partition sum, computed) | sum <= 2^-13.1899 |
| >= 18 | conservative bound (partition sum, computed) | sum <= 2^-23.248 |

The 16-state look-back (the last sixteen stored states; x_c compared
against each) detects every cycle of length 2..16 deterministically at
first return; a detected entry escapes with one charged forced-index
step. Undetected traps (L >= 17) are fenced by the RATCHET: every 2^10
executed steps a walker takes one step whose index is
(x_c >> 5) mod 32 - a position-dependent map, so no trap survives past
2^10 + 16 steps; trapped waste <= 2^-23.248 * (2^10 + 16) = 2^-13.25.

## Waste and success

- L=2 family: rate 1/64, waste <= 3 steps: 3/64.
- L=4: 5 steps: 5 * 2^-10.705. L=6: 7 * 2^-14.484.
- L=8..16: <= 17 * 2^-13.1899.
- trapped L >= 17: <= 2^-13.25.
- distinguished-point missed-report: with DP density 2^-32 and epoch
  E = 2^44, a merged pair fails to report with probability at most
  2 * 2^32 / E = 2^-11 (exact averaging of the geometric reporting over
  the reset phase).
- table overflow: expected inserts 2^95.70 against 2^96 slots,
  P(overflow) <= exp(-2^89) (Chernoff), subtracted.
- useful-visit fraction rho = 1 - (0.0519942318 + 2^-13.25 + 2^-11)
  = 0.9474130892447007 (exact rational 817350603449521898922035003220439151805894052466325190449797183755077/862718293348820473429344482784628181556388621521298319395315527974912, recorded by the experiment).

Birthday ideal with budget t = 1.25 * sqrt(N), N = (r-1)/2:
1 - exp(-(1.25 * rho)^2 / 2) = 0.5040. Walk efficiency: charged 0.85,
supported by (i) the declared toy run measuring the ratio directly
(91/200 = 0.455 recorded vs the same-parameters ideal 0.504: measured
ratio 0.903 >= 0.85), (ii) Teske's r-adding-walk analysis (r = 32
partitions: near-random mixing), (iii) the conservative direction; the
claim survives even at an efficiency charge of 0.80 (success
0.4032).

  modeled success = 0.5040 * 0.85 = 0.4284 >= 0.39.

## Work ledger (worst case every step)

| Item | Instructions | Units |
| --- | ---: | ---: |
| Reference operation W = C + A_j (batched incomplete affine secp256k1 addition (6M, S = M, A = 0)) | reference | 1.000 |
| Parity test of y(W) | 2 | 0.002 |
| Conditional field negation (4-bit) | 12 | 0.012 |
| Index j = x_c mod 32 (low 5 bits) | 2 | 0.002 |
| Representation update (mod-r add + conditional negate + sign flip) | 66 | 0.066 |
| 16-state look-back (15 x 4-limb compares + branches) | 90 | 0.090 |
| Distinguished-point predicate (OR-reduce low 32 bits) | 4 | 0.004 |
| Ratchet check (step counter mod 2^10) | 3 | 0.003 |
| State loads/stores, loop | 69 | 0.069 |
| **Per-step total** | **1,248** | **1.248** |

Escape, trapped, ratchet and reseed steps are walk steps inside t and
charged exactly once here; their visit-loss effect is in rho. One-offs:
seeding 2^16 walkers (2^16 comb multiplications at 64 units) plus setup
<= 2^22.1 units; DP processing t * 2^-32 states at ~1000 instructions
each; detection lag K * 2^33 = 2^49 units; endgame (extraction: two
mod-r multiplications, one inversion, verification) <= 700 units.

W = 1.248 * 1.25 * sqrt((r-1)/2) + < 2^22.1 = 2^128.141546 (rounded up; exact rational 128141547/1000000)

Preprocessing (constants, seeding) = 2^22.0 units ->
preprocessing_log2 = **22.0**, included.

## Memory

Distinguished-point table: 2^96 slots x 128 bytes (x_c 32, a 32, b 32,
metadata 32) = 2^103 bytes; walker histories 2^16 x 17 x 128 B =
2^30.4 bytes; code, constants and scratch positive. Peak storage is
strictly above 2^103 bytes, reported rounded up as memory_log2_bytes =
**103.0**. Overflow aborts (bound above). Nonuniform advice: **0
bytes**.

## Evidence and reproduction

- Part A (production, exact): p and r prime (Miller-Rabin); the
  generator on the curve with order exactly r ([r]G = infinity);
  negation classes (r-1)/2 with -G distinct from G; the 32 constants
  [2^j]P pairwise distinct; the integer-relation premise
  16 * 2^31 < r asserted; the family census recomputed and asserted
  (1, 157, 46,844 for L = 2, 4, 6) with the conservative partition
  bounds recorded (2^-13.1899 for L = 8..16, 2^-23.248 for L >= 18); the exact
  waste, rho, ideal and claimed-success ledger figures computed from
  exact rationals and asserted >= 0.39 (and >= 0.39 at efficiency
  0.80).
- Part B (toy, complete v5 attack on a 20-bit prime-order toy curve,
  production-shaped machinery: look-back 16 with escapes, ratchet every
  64 toy steps, distinguished points at the declared density 1,
  sign-corrected seeding, canonical-state table keyed by x, extraction
  with brute-force verification): declared run
  `{"trials": 200, "dbits": 0, "walkers": 8}`: **91/200 = 0.455**
  recoveries with **0 wrong extractions** (the same-parameters ideal is
  0.504; the measured ratio 0.903 supports the 0.85 efficiency charge);
  2,686 fruitless escapes over 149,679 executed steps (rate 0.01795 per
  executed step; distance-2 rate 0.0161 ~= the 1/64 census prediction;
  the odd-distance events at the toy - 194 at distance 3, 17 at 5, 1 at
  7 - are artifacts of the toy modulus n ~= 2^20, where the constants
  2^j mod n satisfy extra wrap relations; the production premise
  16 * 2^31 < r excludes all such relations, and Part A verifies it).
- Reproduce:
  `echo '{"trials": 200, "dbits": 0, "walkers": 8}' | python3 submissions/secp256k1-ecdlp/negation-rho-v5/experiment.py`.

Small experiments are evidence only at their actual scope.

## Heuristics, limitations

- H1 (uniform independent indices and non-forced parities of generic
  states) is the single walk heuristic; it is the same heuristic that
  yields the classical 1/64 negation 2-cycle rate (reproduced by the
  census and by the toy's distance-2 escape rate) and the birthday
  model.
- The cycle analysis is exact GIVEN H1: families are enumerated exactly
  for L <= 6, bounded conservatively to 2^-13.1899 over 8..16, and the
  undetected tail 2^-23.248 is fenced by the ratchet at 2^-13.25.
- The efficiency charge 0.85 is below the toy-measured 0.903; the claim
  also survives 0.80.
- Instruction sub-costs are worst-case analytical constructions with
  explicit limb counts.
- Single shot at the 0.39 threshold; the modeled success is 0.4284.
- Standard cryptanalysis (negation-map Pollard rho with cycle
  handling); the complete family census for power-of-two constants and
  the ratchet fence are the contribution under this cost model.

References: `tracks/docs/COST_MODEL.md`; HAC ch. 3.6; Teske, J.
Cryptology 14(2), 2001; Wiener-Zuccherato, SAC 1998; van
Oorschot-Wiener (distinguished points); Duval-Leurent (negation-map
cycle handling). Self-contained above.

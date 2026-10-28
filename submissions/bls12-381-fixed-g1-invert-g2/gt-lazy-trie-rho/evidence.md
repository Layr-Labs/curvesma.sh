# Evidence: the experiments of `experiment.py` and the kernel checks

This file is section 7 of `proof.md`; section numbers refer to `proof.md`. `experiment.py` uses the Python standard
library only and reads its configuration from stdin; the declared run is
`echo '{"q":67,"trials":30,"R":24,"branch_bits":6,"sq_cut":0,"dp_bits":4,"walks":4,"seed":2,"production_checks":2,"production_steps":240,"production_dp_bits":4,"kernel_cases":2,"pair_samples":15000,"stat_classes":3000}' | python3 -I experiment.py`
(about 13 s on one core).

Scope of the port: `experiment.py` checks the walk in GT and never computes a pairing. The reduction of section 1
(g = e(P, Q0), answer Q = [x]Q0 in G2, Q = infinity for z = 1) and the G2 output bound of section 5 are analytic
only: they use bilinearity and nondegeneracy of the pinned pairing and a hand count of the double-and-add
formulas, and none of the walk's checks depends on which pairing argument is the instance.

## Parts run on real BLS12-381 GT elements

* **Part A (identities):** GT membership; the Frobenius formula equals the p-th power; the 12-key formula of F2
  equals the keys of the 12 explicit conjugates, all distinct; kappa selects one identical canonical element from
  every conjugate; p^6 inverts; ord_r(p) = 12; 1 + p^4 + p^8 = 0; the r - 1 factorisation; ord_r(2) (log2 250.857,
  smallest t with 2^t in <p>: log2 247.272); the Granger-Scott formula equals y*y (9 Fp2 squarings, 34 additions).
* **Part K (the four kernels in the emulator):** `gen.py` and `gen_lazy.py` emit `mul.s`, `sqr.s`, `lazy.s`,
  `root.s` in memory; `kernels/emu.py` executes them (every instruction counts one; the subset has no branch but the
  final ret, so the executed count is the emitted count). On real GT elements x, m (products of random GT
  generators): `mul.s(x, entry(m)) = x m`, `sqr.s(x) = x^2`, `root.s(x)` returns the stored a, b of x, the six
  pre-sums and (K, kappa) equal to the reference minimum of the twelve stored keys, `lazy.s(x, pre-sums,
  gamma(m))` returns the stored a, b and (K, kappa) of x m, and in the first case K and kappa are also checked against
  keys taken from twelve explicit Frobenius maps (the definition). Every call executes exactly 14,402 / 6,777 /
  4,964 / 1,561 instructions. Declared run: 0 mismatches.
* **Part C (the production lazy walk, with trails and replay):** R = 24 multipliers M_h = g^{a_h} z^{b_h} (40-bit
  exponents, for speed), the 288 twisted entries, trie nodes built on demand (Pi_child = Pi_parent T[h][t]), depth
  L = 6, the 2-cycle escape rule, and the production distinguished-point field (bits 64.. of K) read at a test density
  2^-production_dp_bits: at a distinguished class the record (F, alpha0, beta0), F = bits 136..247 of K, is logged and
  the walk takes a fresh start from new start exponents. Every key comes from the emulated `root.s` (at the root) or
  `lazy.s` (elsewhere), every flush and every escape's explicit representative from the emulated `mul.s`, every escape
  square from the emulated `sqr.s`. An eager walk (explicit element multiplied by T[h][t] each step, keys from F2,
  exponents carried) runs beside it. Asserted at every step: base * Pi_node equals the eager element and (K, kappa)
  equals the reference; after every flush and escape the kernel result equals the Python product; the last trail's
  exponent invariant. Afterwards every record is replayed eagerly from its start pair alone (the walk rule only, no
  stored state): the replay must reach a distinguished class with the recorded F, and the canonical pair
  p^kappa (alpha, beta) must give the canonical representative g^alpha^ z^beta^. Declared run (240 steps, test
  density 1/16): 47 root steps, 193 lazy steps, 32 flushes, one 2-cycle escape whose smaller-key class needed a flush
  from an internal node, 13 records, all 13 replayed to the same F with the canonical pair verified (223 replayed
  steps), 0 key, element or kernel mismatches; chi-square of the multiplier indices 25.26 on 23 degrees of freedom.
* **Part G (the bits of K on production classes):** 3,000 distinct classes (g^(c + i d), i < 3,000): chi-square of
  h over 24 bins 21.71 (23 dof), of bits 64..71 (the low DP-field byte) 237.21 (255 dof), of bits 136..143 (the low
  fingerprint byte) 225.77 (255 dof); the largest |z| of the 248 single-bit frequencies 2.67 (expected maximum of 248
  standard normals about 2.9); every K below p/2. These test marginal statistics only; H2's random-function model of
  the 248 bits is a heuristic (proof.md section 6).
* **`kernels/check_kernels.py`** (not run by the verifier; about 30 s): regenerates the four listings in a temporary
  directory and runs `check.py` (mul.s: 40 GT products, 2 edge, 4 random; sqr.s: 20 cyclotomic elements) and
  `check_lazy.py` (lazy.s and root.s on 4 GT, 2 random Fp12 and 2 edge-value cases, references from explicit
  Frobenius keys): all executed counts constant and equal to the charged ones, 0 mismatches. `check_lazy.py 24 7`
  (12 GT, 6 random, 6 edge cases, inputs drawn from {0, 1, p - 1, p - 2, random}) also passes.

## Toy fields (q = 67, 139, 163, 211, 379)

* **Part B (complete algorithm):** Frobenius classes of the toy GT (order r' the largest prime factor of
  Phi_12(q)), keys and kappa, the twisted table with R = 24 and no squaring branch (sq_cut 0), 2-cycle and checkpoint
  escapes by squaring, every trail ending at its first DP with a record (F, alpha0, beta0) and a fresh start, the
  matcher of 0-ii (whole-record heapsort with its level count checked against the worst-case formula, adjacent dedup,
  pairs under the 2^16 budget, both trails of every pair replayed from their start pairs by the walk's own iteration
  with a check that the replay ends at the recorded fingerprint and within the iteration bound), halting at the first
  verified candidate, g^x = z, the pre-debit meter with a toy budget. The toy walk is eager: part C shows the lazy
  walk visits exactly the eager walk's elements, so the toy tests the walk's collision behaviour, not the kernels.
  The toy's class bits come from a 64-bit hash of the canonical representative (a toy key lies in F_q, q <= 379, far
  too few bits); part G tests the production bits. T_c is the first step at which a class is reached with a
  different canonical exponent pair. Besides the raw fraction with T_c <= S(P), part B reports the fraction in the
  ledger's own units: walk steps up to T_c (wasted 2-cycle steps included) plus every 2-cycle escape at 3.4432 steps
  and every checkpoint escape at its worst case, against S(0.3910)(1 + eps) with the production eps = 0.029892;
  distinguished-point and fresh-start iterations are not counted, since production has one per 2^72 steps. Declared
  run (q = 67, dp_bits 4, 30 trials): 30 logs recovered by the matcher and verified, fraction in the ledger's units
  within S(0.3910)(1 + eps) 0.467 (30 trials: standard error about 0.09), pre-collision 2-cycle rate 0.0031 (36
  escapes in 11,601 steps), heapsort levels at most 0.913 of the worst-case bound, longest replay 112 iterations
  (bound 1,285), at most one pair charged per trial.
* **Part E (pair-collision law):** sample m distinct classes, apply one step F to each, count image coincidences;
  N times the pair-collision rate is what H1 models. Declared run (q = 67, R = 24, s = 0, 15,000 classes): 0.9816
  +- 0.0342 (a 20,000-class run: 0.9751 +- 0.0257) against the first-order value 1 - 1/(12R) = 0.9965 and the conservative value 1 - 1/R = 0.9583
  that section 4 uses.
* **Part F (fruitless closed walks at R = 24; code in `cycles.py`, loaded by `experiment.py`):** q(n), the probability that n independent uniform 12th roots of unity
  sum to zero in Z[zeta_12], computed exactly over the 12^n sequences for n <= 12 (q(2) = 1/12, q(3) = 1/72,
  q(4) = 0.019097, q(6) = 0.0077562, q(8) = 0.0042930, q(10) = 0.0027710, q(12) = 0.0019476), and Q(12) = 0.0019476
  (the largest point probability of a 12-term sum); the multinomial bound of section 3.3 exactly for l <= 288 (P_2 <=
  1/288, P_3 = 2^-15.34, P_4 = 2^-14.76, P_8 = 2^-25.94, P_16 = 2^-43.97, P_32 = 2^-71.04, P_64 = 2^-108.88,
  sum_{3..64} = 6.1659e-5, sum_{65..288} = 2^-108.79) and the Poissonised tail beyond 288 weighted by the probability of
  holding no distinguished class (2^-105.80): useless walks longer than 64 at a rate below 2^-105.63 per step.
* **Part D (ledger):** recomputes, in 60-digit decimal arithmetic, c = 41,821/120,000, the step prices from the
  itemised lists, S, eps from part F, the trie size and its setup charge, the report, fresh-start, restart, switch and
  stopping-allowance terms at the model point P = 0.4002 (log2 S = 125.68269), the claim 124.205 as the least
  thousandth above log2 W = 124.20458 (exact integer inequality), the walk budget and threshold on integers (the walk
  need fits with 0.028% to spare), the block inequalities behind the iteration maximum n_max < 2^125.73, the
  model-rate restart terms of proof-impl.md 0-iii (stuck-trail probability exp(-20/(1 + W_BAR)) = 2^-28.180, long
  useless walks 2^21.73, expected restarts 2^25.644, a factor 2^10.35 below the cap), the expected reports 2^53.73 (a
  factor 2.41 below the cap), the restart-cap worst case as a fraction of S (9.5063e-5) and the success slope 0.6132,
  memory 2^61.92 bytes, the itemised success losses (sum 1.0098e-2, of which 0.01 is the disclosed H1+ budget of
  proof.md section 4), success at least 0.39010 (floored from 0.390101) and sensitivities.
* **Part S** exercises the bounded sampler of 0-iv, and `flags_check.py` runs the multiword compares and counters of
  0-i, 0-ii and 0-iv under the AArch64 NZCV rules (the parent's listings; the parent's squaring-branch compare is
  removed, since this walk has no squaring branch).

## Offline toy censuses of this package's walk (R = 24, s = 0, 4 walks)

Each row pools seeds run with `experiment.py`'s `part_B(random.Random(seed), q, trials, 24, 6, 0, dp_bits, 4)` by an
offline driver (not shipped): 100 trials per seed (25 at dp_bits 8), every seed in each listed range run and counted.
Every run recovered and verified its logarithm through the matcher and replay. "Fraction in ledger units" is the
fraction of runs whose first useful collision falls within S(P)(1 + eps), counted as the ledger counts (walk steps
including wasted 2-cycle steps, each 2-cycle escape at 3.4432 steps, each checkpoint escape at its worst case; no
distinguished-point or fresh-start iterations).

**Census 2 (this package's model point P = 0.3910):**

| q | dp_bits | classes N | seeds | runs | mean T_c/sqrt(N) (model 1.2803) | mean charged steps/sqrt(N) | fraction in ledger units <= S(0.3910)(1+eps) | s.e. | same at S(0.4005)(1+eps) | 2-cycles per step before T_c (model 0.003472) | checkpoint escapes per step before T_c |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 163 | 6 | 4524905 | 5000..5020 | 2100 | 1.2731 | 1.2713 | **0.4133** | 0.0107 | 0.4252 | 0.00339 | 3.28e-05 |
| 163 | 8 | 4524905 | 6000..6032 | 825 | 1.2624 | 1.2748 | **0.4097** | 0.0171 | 0.4206 | 0.00341 | 2.96e-05 |
| 211 | 6 | 721279 | 5000..5027 | 2800 | 1.2496 | 1.2452 | **0.4221** | 0.0093 | 0.4307 | 0.00337 | 2.26e-05 |
| 67 | 6 | 129145 | 5000..5029 | 3000 | 1.2423 | 1.2373 | **0.4273** | 0.0090 | 0.4407 | 0.00335 | 2.06e-05 |

Pooled over all rows: 8725 runs, fraction in ledger units 0.4206 +- 0.0053.

**Census 1 (same walk, earlier package revision; metric at P = 0.4005, so its last-but-two column repeats the
fraction):**

| q | dp_bits | classes N | seeds | runs | mean T_c/sqrt(N) (model 1.2803) | mean charged steps/sqrt(N) | fraction in ledger units <= S(0.4005)(1+eps) | s.e. | same at S(0.4005)(1+eps) | 2-cycles per step before T_c (model 0.003472) | checkpoint escapes per step before T_c |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 163 | 6 | 4524905 | 2000..2042 | 4300 | 1.2605 | 1.2592 | **0.4270** | 0.0075 | 0.4270 | 0.00341 | 3.47e-05 |
| 163 | 8 | 4524905 | 3000..3065 | 1650 | 1.2316 | 1.2448 | **0.4345** | 0.0122 | 0.4345 | 0.00346 | 3.43e-05 |
| 211 | 6 | 721279 | 2000..2054 | 5500 | 1.2501 | 1.2476 | **0.4398** | 0.0067 | 0.4398 | 0.00335 | 3.11e-05 |
| 67 | 6 | 129145 | 2000..2058 | 5900 | 1.2611 | 1.2568 | **0.4273** | 0.0064 | 0.4273 | 0.00334 | 2.33e-05 |

Pooled over all rows: 17350 runs, fraction in ledger units 0.4319 +- 0.0038.

Reading: at the model point S(0.3910)(1 + eps) every row of census 2 reads 0.410 to 0.427 (pooled 0.4206 +- 0.0053,
model 0.3910); since S(0.3910) is below S(0.4002), census 2 also bounds the fraction at the ledger's model point
from below. Census 1 reads 0.427 to 0.440 at S(0.4005)(1 + eps) (pooled 0.4319 +- 0.0038, model 0.4005). The
mean T_c/sqrt(N) is 1.23 to 1.27 against the model's 1.2803. The toy fields' fresh starts (one per 2^dp_bits steps)
collide at rate 1/N rather than 1/N', which favours the toy slightly; the dp_bits 8 rows (0.4% fresh starts) read
0.4097 +- 0.0171 and 0.4345 +- 0.0122 at their model points. The pre-collision 2-cycle rate is 0.00334 to 0.00346
per step against the model 1/288 = 0.00347, and checkpoint escapes (closed walks of length 3..64) occur at 2.1e-5 to
3.5e-5 per step against the bound 6.17e-5 of part F, so the 1.25 margin on eps is not consumed at these sizes. These
are toy-scale statistics with a hash in place of K's bits; they do not measure production success.

## What the experiments do not show

They do not run the production attack, do not measure production success, and check the kernels on samples, not
exhaustively. The toy fields have far fewer classes than GT and their fresh starts are far more frequent than in
production; they test the collision law of this class walk at their sizes, as H1 extrapolates it.

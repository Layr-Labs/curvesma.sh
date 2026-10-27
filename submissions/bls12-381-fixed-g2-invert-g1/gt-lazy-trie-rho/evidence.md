# Evidence: the experiments of `experiment.py` and the kernel checks

This file is section 7 of `proof.md`; section numbers refer to `proof.md`. `experiment.py` uses the Python standard
library only and reads its configuration from stdin; the declared run is
`echo '{"q":67,"trials":30,"R":24,"branch_bits":6,"sq_cut":0,"dp_bits":4,"walks":4,"seed":2,"production_checks":2,"production_steps":240,"production_dp_bits":4,"kernel_cases":2,"pair_samples":15000,"stat_classes":3000}' | python3 -I experiment.py`
(about 9 s on one core).

## Parts run on real BLS12-381 GT elements

* **Part A (identities):** GT membership; the Frobenius formula equals the p-th power; the 12-key formula of F2
  equals the keys of the 12 explicit conjugates, all distinct; kappa selects one identical canonical element from
  every conjugate; p^6 inverts; ord_r(p) = 12; 1 + p^4 + p^8 = 0; the r - 1 factorisation; ord_r(2) (log2 250.857,
  smallest t with 2^t in <p>: log2 247.272); the Granger-Scott formula equals y*y (9 Fp2 squarings, 34 additions).
* **Part K (the four kernels in the emulator):** `gen.py` and `gen_lazy.py` emit `mul.s`, `sqr.s`, `lazy.s`,
  `root.s` in memory; `kernels/emu.py` executes them (every instruction counts one; the subset has no branch but the
  final ret, so the executed count is the emitted count). On real GT elements x, m (products of random GT
  generators): `mul.s(x, entry(m)) = x m`, `sqr.s(x) = x^2`, `root.s(x)` returns the stored a, b of x, the six
  pre-sums, xi_x and (K, kappa) equal to the reference minimum of the twelve stored keys, `lazy.s(x, pre-sums, xi_x,
  gamma(m), C_re, C_im)` returns the stored a, b and (K, kappa) of x m, and in the first case K and kappa are also checked against
  keys taken from twelve explicit Frobenius maps (the definition). Every call executes exactly 14,402 / 6,777 /
  3,527 / 3,425 instructions. Declared run: 0 mismatches.
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
  Frobenius keys; every memory operand of both listings is checked to be an encodable AArch64 immediate): all
  executed counts constant and equal to the charged ones, 0 mismatches. `check_lazy.py 24 7`
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
  ledger's own units: walk steps up to T_c (wasted 2-cycle steps included) plus every 2-cycle escape at 3.9452 steps
  and every checkpoint escape at its worst case, against S(0.4002)(1 + eps) with the production eps = 0.032158;
  distinguished-point and fresh-start iterations are not counted, since production has one per 2^72 steps. Declared
  run (q = 67, dp_bits 4, 30 trials): 30 logs recovered by the matcher and verified, fraction in the ledger's units
  within S(0.4002)(1 + eps) 0.4667 (14 of 30; standard error about 0.09), pre-collision 2-cycle rate 0.003103 (36
  escapes in 11,601 steps), heapsort levels at most 0.913 of the worst-case bound, longest replay 112 iterations
  (bound 1,285), at most one pair charged per trial.
* **Part E (pair-collision law):** sample m distinct classes, apply one step F to each, count image coincidences;
  N times the pair-collision rate is what H1 models. Declared run (q = 67, R = 24, s = 0, 15,000 classes): 0.9816
  +- 0.0342 (a 20,000-class run: 0.9751 +- 0.0257) against the exact first-order value 1 - 1/(12R) = 0.9965 and the conservative value 1 - 1/R = 0.9583
  that section 4 uses.
* **Part F (fruitless closed walks at R = 24; code in `cycles.py`, loaded by `experiment.py`):** q(n), the probability that n independent uniform 12th roots of unity
  sum to zero in Z[zeta_12], computed exactly over the 12^n sequences for n <= 12 (q(2) = 1/12, q(3) = 1/72,
  q(4) = 0.019097, q(6) = 0.0077562, q(8) = 0.0042930, q(10) = 0.0027710, q(12) = 0.0019476), and Q(12) = 0.0019476
  (the largest point probability of a 12-term sum); the multinomial bound of section 3.3 exactly for l <= 288 (P_2 =
  1/288 exactly, P_3 = 2^-15.34, P_4 = 2^-14.76, P_8 = 2^-25.94, P_16 = 2^-43.97, P_32 = 2^-71.04, P_64 = 2^-108.88,
  sum_{3..64} = 6.1659e-5, sum_{65..288} = 2^-108.79) and the Poissonised tail beyond 288 weighted by the probability of
  holding no distinguished class (2^-105.80): useless walks longer than 64 at a rate below 2^-105.63 per step.
* **Part D (ledger):** recomputes, in 60-digit decimal arithmetic, c = 36,500/120,000, the step prices from the
  itemised lists, S, eps from part F, the trie size and its setup charge, the report, fresh-start, restart, switch and
  stopping-allowance terms, the claim 124.012 as the least thousandth above log2 W = 124.01142 (exact integer
  inequality), the walk budget and threshold on integers (the walk need fits with 0.040% to spare), the block
  inequalities behind the iteration maximum n_max < 2^125.73, the expected restarts (stuck 2^24.88, long useless walks
  2^21.73), memory below 2^61.93 bytes, the itemised success losses (sum 0.010095, of which the H1+ allowance is exactly
  0.01 and the restart cap's worst case 9.21e-5; there is no tail bound to evaluate), the H1+ drift table and the
  drift at which the floor 0.39 is lost (0.6608 points, bisected on exact decimals), success at least 0.39010
  (floored) and sensitivities.
* **Part S** exercises the bounded sampler of 0-iv, and `flags_check.py` runs the multiword compares and counters of
  0-i, 0-ii and 0-iv under the AArch64 NZCV rules (the parent's listings; the parent's squaring-branch compare is
  removed, since this walk has no squaring branch).

## Offline toy censuses of this package's walk (R = 24, s = 0, 4 walks)

Every row below was run at this package's weights (eps = 0.032158, a 2-cycle escape at 3.9452 steps) and at this
package's model point, so no figure here carries over from an earlier revision. Each row is
`part_B(random.Random(seed), q, 3000, 24, 6, 0, 6, 4)` by an offline driver (not shipped): 3,000 trials per seed,
three seeds per field, every seed listed and counted and no run dropped. Every run recovered and verified its
logarithm through the matcher and replay. "Fraction in ledger units" is the fraction of runs whose first useful
collision falls within S(0.4002)(1 + eps), counted as the ledger counts (walk steps including wasted 2-cycle steps,
each 2-cycle escape at 3.9452 steps, each checkpoint escape at its worst case; no distinguished-point or fresh-start
iterations).

**The rule, fixed before any count of these runs was read:** per field, the pooled fraction minus 4 binomial
standard errors must be at least 0.391. No seed was added, dropped or re-sized after a count was seen, and P_model
was not moved afterwards.

**Census 3 (this package's model point P = 0.4002, dp_bits 6; the registered census):**

| q | dp_bits | classes N | seed | runs | mean T_c/sqrt(N) (model 1.2803) | mean charged steps/sqrt(N) | fraction in ledger units <= S(0.4002)(1+eps) | s.e. | 2-cycles per step before T_c (model 0.003472) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 67 | 6 | 129145 | 7000 | 3000 | 1.2473 | 1.2446 | **0.4397** | 0.0091 | 0.00339 |
| 67 | 6 | 129145 | 7001 | 3000 | 1.2584 | 1.2576 | **0.4337** | 0.0090 | 0.00330 |
| 67 | 6 | 129145 | 7002 | 3000 | 1.2538 | 1.2545 | **0.4463** | 0.0091 | 0.00338 |
| 163 | 6 | 4524905 | 7000 | 3000 | 1.2535 | 1.2537 | **0.4390** | 0.0091 | 0.00336 |
| 163 | 6 | 4524905 | 7001 | 3000 | 1.2574 | 1.2582 | **0.4330** | 0.0090 | 0.00340 |
| 163 | 6 | 4524905 | 7002 | 3000 | 1.2640 | 1.2650 | **0.4297** | 0.0090 | 0.00339 |
| 211 | 6 | 721279 | 7000 | 3000 | 1.2500 | 1.2491 | **0.4440** | 0.0091 | 0.00334 |
| 211 | 6 | 721279 | 7001 | 3000 | 1.2683 | 1.2699 | **0.4227** | 0.0090 | 0.00339 |
| 211 | 6 | 721279 | 7002 | 3000 | 1.2664 | 1.2665 | **0.4303** | 0.0090 | 0.00338 |

Pooled per field at the registered rule (fraction minus 4 standard errors must be at least 0.391):

| q | runs | fraction | s.e. | fraction - 4 s.e. | passes 0.391 |
| ---: | ---: | ---: | ---: | ---: | --- |
| 67 | 9000 | 0.4399 | 0.0052 | **0.4190** | yes |
| 163 | 9000 | 0.4339 | 0.0052 | **0.4130** | yes |
| 211 | 9000 | 0.4323 | 0.0052 | **0.4114** | yes |

All 27,000 runs pooled: fraction 0.4354 +- 0.0030. The model point is 0.4002, so each pooled field reads 0.0321 to
0.0397 above its model.

**Not re-run at this model point, and not depended on:** the two tables below are the earlier revisions' censuses
(previous package's weights, escape at 3.4432 steps and eps = 0.029892; and an earlier model point S(0.4005)). They
are provenance only. Nothing in the ledger, the drift table or the success floor uses them; census 3 above carries
the support this revision claims.

**Superseded census 2 (earlier weights, metric at P = 0.3910):**

| q | dp_bits | classes N | seeds | runs | mean T_c/sqrt(N) (model 1.2803) | fraction in ledger units <= S(0.3910)(1+eps) | s.e. | 2-cycles per step before T_c |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 163 | 6 | 4524905 | 5000..5020 | 2100 | 1.2731 | **0.4133** | 0.0107 | 0.00339 |
| 163 | 8 | 4524905 | 6000..6032 | 825 | 1.2624 | **0.4097** | 0.0171 | 0.00341 |
| 211 | 6 | 721279 | 5000..5027 | 2800 | 1.2496 | **0.4221** | 0.0093 | 0.00337 |
| 67 | 6 | 129145 | 5000..5029 | 3000 | 1.2423 | **0.4273** | 0.0090 | 0.00335 |

Pooled: 8725 runs, fraction in ledger units 0.4206 +- 0.0053.

**Superseded census 1 (earlier package revision, metric at P = 0.4005):**

| q | dp_bits | classes N | seeds | runs | mean T_c/sqrt(N) (model 1.2803) | fraction in ledger units <= S(0.4005)(1+eps) | s.e. | 2-cycles per step before T_c |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 163 | 6 | 4524905 | 2000..2042 | 4300 | 1.2605 | **0.4270** | 0.0075 | 0.00341 |
| 163 | 8 | 4524905 | 3000..3065 | 1650 | 1.2316 | **0.4345** | 0.0122 | 0.00346 |
| 211 | 6 | 721279 | 2000..2054 | 5500 | 1.2501 | **0.4398** | 0.0067 | 0.00335 |
| 67 | 6 | 129145 | 2000..2058 | 5900 | 1.2611 | **0.4273** | 0.0064 | 0.00334 |

Pooled: 17350 runs, fraction in ledger units 0.4319 +- 0.0038.

Reading: every row of census 3 reads 0.4227 to 0.4463 and every pooled field clears the registered rule at 0.4114 or
better, so P_model = 0.4002 is supported at toy scale with a margin. The mean T_c/sqrt(N) is 1.25 to 1.27 against
the model's 1.2803. The toy fields' fresh starts (one per 2^dp_bits steps) collide at rate 1/N rather than 1/N',
which favours the toy slightly. The pre-collision 2-cycle rate is 0.00330 to 0.00340 per step against the model
1/288 = 0.003472, so the 1.25 margin on eps is not consumed at these sizes. These are toy-scale statistics with a
hash in place of K's bits; they do not measure production success.

## What the experiments do not show

They do not run the production attack, do not measure production success, and check the kernels on samples, not
exhaustively. The toy fields have far fewer classes than GT and their fresh starts are far more frequent than in
production; they test the collision law of this class walk at their sizes, as H1 extrapolates it.

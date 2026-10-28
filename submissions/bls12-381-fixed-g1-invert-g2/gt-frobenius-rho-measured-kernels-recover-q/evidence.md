# Evidence: the experiments of `experiment.py` and the kernel checks

This file is section 7 of `proof.md`; section numbers refer to `proof.md`. `experiment.py` uses the Python standard
library only and reads its configuration from stdin; the declared run is
`echo '{"q":67,"trials":30,"R":32,"branch_bits":6,"sq_cut":0,"dp_bits":4,"walks":4,"seed":2,"production_checks":2,"production_steps":240,"production_dp_bits":4,"kernel_cases":2,"pair_samples":15000,"stat_classes":3000}' | python3 -I experiment.py`
(under 20 s on one core; about 6 s measured).

## Parts run on real BLS12-381 GT elements

* **Part A (identities):** GT membership; the Frobenius formula equals the p-th power; the 12-key formula of F2
  equals the keys of the 12 explicit conjugates, all distinct; kappa selects one identical canonical element from
  every conjugate; p^6 inverts; ord_r(p) = 12; 1 + p^4 + p^8 = 0; the r - 1 factorisation; ord_r(2) (log2 250.857,
  smallest t with 2^t in <p>: log2 247.272); the Granger-Scott formula equals y*y (9 Fp2 squarings, 34 additions).
* **Part K (the four kernels in the emulator):** `gen.py` and `gen_lazy.py` emit `mul.s`, `sqr.s`, `lazy.s`,
  `root.s` in memory; `kernels/emu.py` executes them (every instruction counts one; the subset has no branch but the
  final ret, so the executed count is the emitted count). On real GT elements x, m (products of random GT
  generators): `mul.s(x, entry(m)) = x m`, `sqr.s(x) = x^2`, `root.s(x)` returns the stored a, b of x, writes the
  per-base data (six pre-sums mod p, Xa, Xd; equal to `gen_lazy.block_vectors(x)`) and (K, kappa) equal to the
  reference minimum of the twelve stored keys, `lazy.s(x, per-base data, node entry of m)` returns the stored a, b
  and (K, kappa) of x m, and in the first case K and kappa are also checked against
  keys taken from twelve explicit Frobenius maps (the definition); `mul.s` is also run in the trie-construction
  orientation (the multiplier's coordinates against the node's 45 stored vectors). Every call executes exactly 12,061 /
  6,786 / 3,457 / 3,074 instructions. Declared run: 0 mismatches.
* **Part C (the production lazy walk, with trails and replay):** R = 32 multipliers M_h = g^{a_h} z^{b_h} (40-bit
  exponents, for speed), the 384 twisted entries, trie nodes built on demand (Pi_child = Pi_parent T[h][t]), depth
  L = 6, the 2-cycle escape rule, and the production distinguished-point field (bits 64.. of K) read at a test density
  2^-production_dp_bits: at a distinguished class the record (F, alpha0, beta0), F = bits 150..277 of K, is logged and
  the walk takes a fresh start. Every key comes from the emulated `root.s` or `lazy.s`, every flush and escape from
  the emulated `mul.s` and `sqr.s`; an eager walk (keys from F2, exponents carried) runs beside it. Asserted at every
  step: base * Pi_node equals the eager element and (K, kappa) the reference; every kernel result equals the Python
  product; the exponent invariant. Every record is then replayed eagerly from its start pair alone to a distinguished
  class with the recorded F, and p^kappa (alpha, beta) must give the canonical representative g^alpha^ z^beta^. Declared run (240 steps, test
  density 1/16): 49 root steps, 191 lazy steps, 26 flushes, one 2-cycle escape whose smaller-key class needed a flush
  from an internal node, 22 records, all 22 replayed to the same F with the canonical pair verified (218 replayed
  steps), 0 key, element or kernel mismatches; chi-square of the multiplier indices 26.76 on
  31 degrees of freedom.
* **Part G (the bits of K on production classes):** 3,000 distinct classes (g^(c + i d), i < 3,000): chi-square of
  h over 32 bins 39.06 (31 dof), of bits 64..71 (the low DP-field byte) 278.85 (255 dof), of bits 150..157 (the low
  fingerprint byte) 230.04 (255 dof); the largest |z| of the 278 single-bit frequencies 3.43 (expected maximum of 278
  standard normals about 2.9, so this is a mild excess on the maximum, not a dependence signal); every K below p/2. These test marginal statistics only; H2's random-function model of
  the 278 bits is a heuristic (proof.md section 6).
* **`kernels/check_kernels.py`** (not run by the verifier; about 21 s): regenerates the four listings in a temporary
  directory, checks that every immediate of every emitted instruction encodes in AArch64 (memory offsets: ldp/stp
  -512..504 and ldr/str 0..32760, multiples of 8; add/sub/cmp 0..4095 with lsl 0 or 12; mov a 16-bit chunk or its
  inverse; extr 0..63: 1,881 / 1,133 / 410 / 389 immediates, all encodable), and runs `check.py` (mul.s: 40 GT products, 2 edge, 4 random, and 19
  entries set directly at the extremes, all p - 1 and per output its subtracted products at p - 1, checked against the
  45-product plan over Fp; sqr.s: 20 cyclotomic elements) and
  `check_lazy.py` (lazy.s and root.s on 4 GT, 2 random Fp12 and 2 edge-value cases, references from explicit
  Frobenius keys): all executed counts constant and equal to the charged ones, 0 mismatches. `check_lazy.py 8 11
  mutants` passes with 0 lazy and 0 root mismatches, and its mutants mode catches all ten deliberately broken
  listings (top-limb-only minimum, `ls` instead of `lo`, a dropped `csel`, a flipped `b` select, unreduced pre-sums
  in root.s, an operand sum without its carry, `d_T` without the doubling of S0, and three lazy.s listings generated
  with one constant changed: Xd read for Xa, Ed read for Ea, the Winograd pairs regrouped).

## Toy fields (q = 67, 139, 163, 211, 379)

* **Part B (complete algorithm):** Frobenius classes of the toy GT (order r' the largest prime factor of
  Phi_12(q)), keys and kappa, the twisted table with R = 32 and no squaring branch (sq_cut 0), 2-cycle and checkpoint
  escapes by squaring, every trail ending at its first DP with a record (F, alpha0, beta0) and a fresh start, the
  matcher of 0-ii (whole-record heapsort with its level count checked against the worst-case formula, adjacent dedup,
  pairs under a toy budget of 2^16 (never binding: at most one pair charged per trial), both trails of every pair
  replayed from their start pairs by the walk's own iteration with a check that the replay ends at the recorded
  fingerprint and within the iteration bound), halting at the first verified candidate, g^x = z, the pre-debit meter
  with a toy budget of 2^40 instructions. Toy limits, not the production caps: a trail restarts after more than
  40 x 2^dp_bits steps without a distinguished class; no report cap. The toy walk is eager: part C shows the lazy
  walk visits exactly the eager walk's elements, so the toy tests the walk's collision behaviour, not the kernels.
  The toy's class bits come from a 64-bit hash of the canonical representative (a toy key lies in F_q, q <= 379, far
  too few bits); part G tests the production bits. T_c is the first step at which a class is reached with a
  different canonical exponent pair. Besides the raw fraction with T_c <= S(P), part B reports the fraction of runs
  whose charged count at T_c (defined below the census table) is within S(0.4005)(1 + eps), the production eps
  rounded up, 0.02406268. Declared
  run (q = 67, dp_bits 4, 30 trials): 30 logs recovered by the matcher and verified, fraction by the charged count
  0.4333 (30 trials: standard error about 0.09), pre-collision 2-cycle rate 0.0036 (44 in
  12,199 steps against about 32 at the model 1/384, about 2.2 standard deviations at this sample size; the
  census below, 0.002480 to 0.002588 over its rows, is the measurement of this rate), heapsort levels at most 0.9066 of the worst-case bound, longest replay 141 iterations
  (bound 1,285), at most one pair charged per trial.
* **Part E (pair-collision law):** sample m distinct classes, apply one step F to each, count image coincidences;
  N times the pair-collision rate is what H1 models. Declared run (q = 67, R = 32, s = 0, 15,000 classes): 1.0447
  (relative standard error 0.0331); with the input
  `{"q":67,"trials":0,"R":32,"branch_bits":6,"sq_cut":0,"dp_bits":4,"walks":4,"seed":2,"production_checks":0,"production_steps":0,"kernel_cases":1,"pair_samples":40000,"stat_classes":0}`
  0.9786 (relative s.e. 0.0128), and with `"pair_samples":60000` in it 0.9863 (relative s.e. 0.0085); against the
  exact first-order value 1 - 1/(12R) = 0.9974 and the conservative value 1 - 1/R = 0.9688 that section 4 uses.
* **Part F (fruitless closed walks at R = 32; code in `cycles.py`, loaded by `experiment.py`):** q(n), the probability that n independent uniform 12th roots of unity
  sum to zero in Z[zeta_12], computed exactly over the 12^n sequences for n <= 12 (q(2) = 1/12, q(3) = 1/72,
  q(4) = 0.019097, q(6) = 0.0077562, q(8) = 0.0042930, q(10) = 0.0027710, q(12) = 0.0019476), and Q(12) = 0.0019476
  (the largest point probability of a 12-term sum); the multinomial bound of section 3.3 exactly for l <= 384 (P_2 =
  1/384 exactly, P_3 = 2^-16.17, P_4 = 2^-15.59, P_8 = 2^-27.61, P_16 = 2^-47.37, P_32 = 2^-77.92, P_64 = 2^-122.33,
  sum_{3..64} = 3.4485e-5, sum_{65..384} = 2^-122.62) and the Poissonised tail beyond 384 weighted by the probability of
  holding no distinguished class (2^-156.83): useless walks longer than 64 at a rate below 2^-122.62 per step.
* **Part R (`reduction.py`, BLS12-381):** a challenge-style P = [b]P0, g = e(P, Q0) by the pinned pairing(Q0, P)
  counted operation by operation (10,331 operations; g in GT, g != 1), Q = [x]Q0 by the pinned G2 formulas with its
  counted cost and the hard maximum 12,809 units over every 255-bit scalar, Q on the G2 curve, [r]Q = infinity,
  e(P, Q) = z, and the identity target (z = 1, Q = infinity, e(P, infinity) = 1). All asserted.
* **Part D (ledger):** recomputes in 60-digit decimal and exact fraction arithmetic every figure of proof.md
  sections 3.3 to 5 and asserts each bound it computes: the prices, S, eps and the H2 rates of X, the block
  lemma (a) with the escape-ended block table, the threshold inequality (b), the drift tolerance and table (c), the
  iteration maximum, the caps' expectations and Markov ratios, the claim 123.868 above log2 W < 123.86742 (exact
  integer inequality), memory, the itemised losses (sum below 0.01008429) and success at least 0.39041 (floored).
* **Part S** exercises the bounded sampler of 0-iv, and `flags_check.py` runs the multiword compares and counters of
  0-i, 0-ii and 0-iv under the AArch64 NZCV rules (the sibling package's listings; this walk has no squaring-branch
  compare).

## Offline toy census of this package's walk (R = 32, s = 0, 4 walks)

Each row pools seeds run with `experiment.py`'s `part_B(random.Random(seed), q, trials, 32, 6, 0, dp_bits, 4)` by an
offline driver (not shipped) that asserts this package's charges in the file it loads: 100 trials per seed (25 at
dp_bits 8), every seed in each listed range run and counted. Every run recovered and verified its logarithm through
the matcher and replay.

**Registration.** Two blocks of six rows were registered for this package's charges (block 33,458 instructions, flush
step 15,706, escape 24,000, toy eps 0.02406268 = this ledger's eps rounded up, model point 0.4005), each before any
count of that block was read; block 2 was registered after block 1's result. Their seeds (30000.. and 31000..) had
been run by no earlier census of this walk. The rule as registered for each block: the pooled fraction is successes
over runs over its six rows at S(0.4005)(1 + eps), and the block passes when the pooled fraction less four binomial
standard errors is at least the package's floored success bound, 0.39041. Block 2's registration allows pooling it
with block 1 only under a rule stated before reading, and none was, so the blocks are reported separately. One seed
of each block re-run on the shipped `experiment.py` (block 1: q = 67, dp_bits 8, seed 30000; block 2: q = 67,
dp_bits 6, seed 31000) reproduces its recorded line field for field. No row was dropped.

**What the census counts.** A run's charged count at T_c (its first useful collision) is the toy's
`ord_at_Tc + 4.30390 esc2 + 184.563 escL`: one step for every iteration that took a step (useful ones and the wasted
steps inside cycles), each 2-cycle escape at 24,000 instructions in steps of c C, each checkpoint escape at 64 steps
of 15,706 plus 24,000 (184.563 steps). Against U + X at T_c (proof.md 3.3) it omits the distinguished-point iterations
(one per 2^dp_bits toy steps; one per 2^86 in production), the iteration that detects a cycle (2-cycles and
checkpoint escapes, at most 0.0027 per counted step in every row) and, per length-64 checkpoint escape, the 65 x 64
= 4,160 bookkeeping instructions the meter adds (under 0.75 step; it charges 185.309); X's charge of 2 steps per 2-cycle
and 64 + l per longer cycle bounds the steps it counts. So the count is below the production charge by less than
1.9% of it at dp_bits 6 and 0.7% at 8 (2^-dp_bits per executed step, plus the detections and those 4,160 instructions):
the census event is that much easier than production's. The column "waste per useful iteration, at most"
bounds each row's realized f by charging every escape as X does (2 + 4.30391 steps per 2-cycle, 128 + 185.309 per
checkpoint escape) over its useful iterations (counted steps less one per 2-cycle and 128 per checkpoint escape),
rounded up. The success rule reads the column at the model point S(0.4005)(1 + eps).
Rounding: success fractions (bold) and the pools are floored to four decimals, the s.e. sqrt(p(1 - p)/runs) rounded
up, the rest to nearest (the means weight the per-seed means `experiment.py` prints to four decimals); ranges quoted
below are rounded outward, success fractions floored. Classes N: 129,145 (q = 67), 4,524,905 (163), 721,279 (211).

| block | q | dp_bits | seeds | runs | successes | mean T_c/sqrt(N) (model 1.2734) | mean charged steps/sqrt(N) | at the model point S(0.4005)(1+eps) | s.e. | 2-cycles per step before T_c (model 0.002604) | checkpoint escapes per step before T_c | waste per useful iteration, at most |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 163 | 6 | 30000..30020 | 2100 | 928 | 1.2519 | 1.2472 | **0.4419** | 0.0109 | 0.002556 | 1.89e-05 | 0.0222 |
| 1 | 211 | 6 | 30000..30032 | 3300 | 1397 | 1.2638 | 1.2577 | **0.4233** | 0.0087 | 0.002546 | 1.47e-05 | 0.0208 |
| 1 | 67 | 6 | 30000..30032 | 3300 | 1411 | 1.2407 | 1.2355 | **0.4275** | 0.0087 | 0.002527 | 1.87e-05 | 0.0219 |
| 1 | 163 | 8 | 30000..30032 | 825 | 338 | 1.2864 | 1.2961 | **0.4096** | 0.0172 | 0.002569 | 1.65e-05 | 0.0215 |
| 1 | 211 | 8 | 30000..30032 | 825 | 329 | 1.2708 | 1.2802 | **0.3987** | 0.0171 | 0.002530 | 1.70e-05 | 0.0214 |
| 1 | 67 | 8 | 30000..30032 | 825 | 344 | 1.2595 | 1.2721 | **0.4169** | 0.0172 | 0.002588 | 2.97e-05 | 0.0258 |
| 2 | 163 | 6 | 31000..31020 | 2100 | 903 | 1.2643 | 1.2598 | **0.4300** | 0.0109 | 0.002534 | 2.02e-05 | 0.0225 |
| 2 | 211 | 6 | 31000..31032 | 3300 | 1454 | 1.2494 | 1.2446 | **0.4406** | 0.0087 | 0.002538 | 1.92e-05 | 0.0222 |
| 2 | 67 | 6 | 31000..31032 | 3300 | 1428 | 1.2484 | 1.2421 | **0.4327** | 0.0087 | 0.002532 | 1.31e-05 | 0.0202 |
| 2 | 163 | 8 | 31000..31032 | 825 | 341 | 1.2597 | 1.2695 | **0.4133** | 0.0172 | 0.002561 | 1.82e-05 | 0.0220 |
| 2 | 211 | 8 | 31000..31032 | 825 | 345 | 1.2688 | 1.2785 | **0.4181** | 0.0172 | 0.002480 | 1.81e-05 | 0.0215 |
| 2 | 67 | 8 | 31000..31032 | 825 | 328 | 1.2664 | 1.2753 | **0.3975** | 0.0171 | 0.002576 | 1.34e-05 | 0.0206 |

Block 1 (seeds 30000..): **11,175 runs, 4,747 successes, fraction 0.4247 +- 0.00468** at the model point
S(0.4005)(1 + eps); four binomial standard errors below it, **0.4060** (0.40608, floored), above 0.39041.
Block 2 (seeds 31000..): **11,175 runs, 4,799 successes, fraction 0.4294 +- 0.00469**; four binomial standard errors
below it, **0.4107** (0.41071, floored), above 0.39041.

Reading: the twelve rows run from 0.3975 (block 2, q = 67, dp_bits 8, seeds 31000..31032, 825 runs, the lowest) to
0.4419 (block 1, q = 163, dp_bits 6, seeds 30000..30020); the spread across each block's rows is within binomial noise
(chi-square 5.91 and 6.55 on 5 degrees of freedom against the block's pooled fraction), and the dp_bits 8 rows (fresh
starts at rate 2^-8 per step) are the noisiest. The mean T_c/sqrt(N) is 1.240 to 1.287 against the model's 1.2734.
The toy fields' fresh starts (one per 2^dp_bits steps) collide at rate 1/N rather than 1/N', which favours the toy.
The pre-collision 2-cycle rate is 0.002480 to 0.002588 per step against the model 1/384 = 0.002604, and checkpoint
escapes (closed walks of length 3..64) occur at 1.30e-5 to 2.97e-5 per step against the bound 3.45e-5 of part F.
The realized waste per useful iteration over the twelve rows is at most 0.0202 to 0.0258 (the table's cells),
1.03 to 1.33 times the H2 rate 0.0193979 and at most 0.63 of the tolerance 0.04103 of H1+ (proof.md 3.4 (c)); one
row (block 1, q = 67, dp_bits 8, 825 runs) is above eps, 1.2404 times the H2 rate; the tolerance, not eps, covers
it. These toy-scale statistics, with a hash in place of K's bits, do not measure production success; they
support the threshold and waste part of H1+ only, not any cap (the caps rest on Markov bounds, which a toy census
cannot measure).

## What the experiments do not show

They do not run the production attack, do not measure production success, and check the kernels on samples, not
exhaustively. The toy fields have far fewer classes than GT and their fresh starts are far more frequent than in
production; they test the collision law of this class walk at their sizes, as H1 extrapolates it.

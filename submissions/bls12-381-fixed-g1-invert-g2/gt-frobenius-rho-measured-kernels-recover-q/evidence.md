# Evidence: the experiments of `experiment.py` and the kernel checks

This file is section 7 of `proof.md`; section numbers refer to `proof.md`. `experiment.py` uses the Python standard
library only and reads its configuration from stdin; the declared run is
`echo '{"q":67,"trials":30,"R":32,"branch_bits":6,"sq_cut":0,"dp_bits":4,"walks":4,"seed":2,"production_checks":2,"production_steps":240,"production_dp_bits":4,"kernel_cases":2,"pair_samples":15000,"stat_classes":3000}' | python3 -I experiment.py`
(under 20 s on one core; 6 to 10 s measured).

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
  orientation (the multiplier's coordinates against the node's stored forms). Every call executes exactly 12,735 /
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
* **`kernels/check_kernels.py`** (not run by the verifier; about 26 s): regenerates the four listings in a temporary
  directory, checks that every immediate of every emitted instruction encodes in AArch64 (memory offsets: ldp/stp
  -512..504 and ldr/str 0..32760, multiples of 8; add/sub/cmp 0..4095 with lsl 0 or 12; mov a 16-bit chunk or its
  inverse; extr 0..63: 1,706 / 1,133 / 410 / 389 immediates, all encodable), and runs `check.py` (mul.s: 40 GT products, 2 edge, 4 random; sqr.s: 20 cyclotomic elements) and
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
  rounded up, 0.02379221. Declared
  run (q = 67, dp_bits 4, 30 trials): 30 logs recovered by the matcher and verified, fraction by the charged count
  0.4333 (30 trials: standard error about 0.09), pre-collision 2-cycle rate 0.0036 (44 in
  12,199 steps against about 32 at the model 1/384, about 2.2 standard deviations at this sample size; the
  census below, 0.002518 to 0.002608 over its rows, is the measurement of this rate), heapsort levels at most 0.9066 of the worst-case bound, longest replay 141 iterations
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
  iteration maximum, the caps' expectations and Markov ratios, the claim 123.896 above log2 W = 123.89581 (exact
  integer inequality), memory, the itemised losses (sum 0.01008428) and success at least 0.39041 (floored).
* **Part S** exercises the bounded sampler of 0-iv, and `flags_check.py` runs the multiword compares and counters of
  0-i, 0-ii and 0-iv under the AArch64 NZCV rules (the sibling package's listings; this walk has no squaring-branch
  compare).

## Offline toy census of this package's walk (R = 32, s = 0, 4 walks)

Each row pools seeds run with `experiment.py`'s `part_B(random.Random(seed), q, trials, 32, 6, 0, dp_bits, 4)` by an
offline driver (not shipped): 100 trials per seed (25 at dp_bits 8), every seed in each listed range run and counted.
Every run recovered and verified its logarithm through the matcher and replay.

**Registration.** The eight rows of the table and the pooling rule were registered for this package's charges
(block 34,132 instructions, flush step 16,380, escape 24,000, toy eps 0.02379221 = this ledger's eps rounded up,
model point 0.4005) before any count at those charges was read. The rule as registered: the pooled fraction is
successes over runs over all eight rows at S(0.4005)(1 + eps), and it passes when the pooled fraction less four
binomial standard errors is at least 0.39040, the package's floored success bound when the rule was written (the cap
repair of section 3.4 has since moved that bound to 0.39041). Disclosed: the seeds are those of two earlier
registrations on this walk, whose counts at earlier charges had been read. The walk itself is seed for seed the same,
so only the charges are new and no new sample was drawn; one seed re-run on the shipped `experiment.py` (q = 67,
dp_bits 6, seed 5000) reproduces its recorded row exactly. No row was dropped.

**What the census counts.** A run's charged count at T_c (its first useful collision) is the toy's
`ord_at_Tc + 4.21891 esc2 + 188.501 escL`: one step for every iteration that took a step (useful ones and the wasted
steps inside cycles), each 2-cycle escape at 24,000 instructions in steps of c C, each checkpoint escape at 64 steps
of 16,380 plus 24,000 (188.501 steps). Against U + X at T_c (proof.md 3.3) it omits the distinguished-point iterations
(one per 2^dp_bits toy steps; one per 2^86 in production), the iteration that detects a cycle (2-cycles and
checkpoint escapes, at most 0.0027 per counted step in every row) and below 0.74 step per length-64 checkpoint escape
(the meter charges 189.232); X's charge of 2 steps per 2-cycle and 64 + l per longer cycle bounds the steps it
counts. So the count is below the production charge by about 1.9% of it at dp_bits 6 and 0.7% at 8 (2^-dp_bits per
executed step, plus the detections): the census event is that much easier than production's. The column "waste per useful iteration, at most"
bounds each row's realized f by charging every escape as X does (2 + 4.2189 steps per 2-cycle, 128 + 189.233 per
checkpoint escape) over its useful iterations (counted steps less one per 2-cycle and 128 per checkpoint escape),
rounded up. The success rule reads the column at the model point S(0.4005)(1 + eps).
Rounding: success fractions (bold) and the pool are floored to four decimals, the s.e. sqrt(p(1 - p)/runs) rounded
up, the rest to nearest (the means weight the per-seed means `experiment.py` prints to four decimals); ranges quoted
below are rounded outward, success fractions floored.

| q | dp_bits | classes N | seeds | runs | mean T_c/sqrt(N) (model 1.2734) | mean charged steps/sqrt(N) | at the model point S(0.4005)(1+eps) | s.e. | 2-cycles per step before T_c (model 0.002604) | checkpoint escapes per step before T_c | waste per useful iteration, at most |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 163 | 6 | 4524905 | 5000..5020 | 2100 | 1.2649 | 1.2598 | **0.4204** | 0.0108 | 0.002534 | 1.87e-05 | 0.0219 |
| 163 | 8 | 4524905 | 6000..6032 | 825 | 1.2610 | 1.2708 | **0.4060** | 0.0171 | 0.002608 | 1.82e-05 | 0.0221 |
| 211 | 6 | 721279 | 5000..5027 | 2800 | 1.2523 | 1.2461 | **0.4235** | 0.0094 | 0.002536 | 1.44e-05 | 0.0205 |
| 67 | 6 | 129145 | 5000..5029 | 3000 | 1.2655 | 1.2605 | **0.4163** | 0.0090 | 0.002583 | 1.94e-05 | 0.0224 |
| 163 | 6 | 4524905 | 7000..7020 | 2100 | 1.2396 | 1.2347 | **0.4204** | 0.0108 | 0.002542 | 1.89e-05 | 0.0220 |
| 163 | 8 | 4524905 | 8000..8032 | 825 | 1.2719 | 1.2833 | **0.4145** | 0.0172 | 0.002566 | 2.43e-05 | 0.0239 |
| 211 | 6 | 721279 | 7000..7027 | 2800 | 1.2404 | 1.2342 | **0.4296** | 0.0094 | 0.002518 | 1.45e-05 | 0.0204 |
| 67 | 6 | 129145 | 7000..7020 | 2100 | 1.2602 | 1.2540 | **0.4152** | 0.0108 | 0.002590 | 1.39e-05 | 0.0207 |

Pooled over the eight rows: **16,550 runs, 6,953 successes, fraction 0.4201 +- 0.00384** at the model point
S(0.4005)(1 + eps); four binomial standard errors below the pooled value, **0.4047** (0.40477, floored), above the
registered bar 0.39040 and the package's floored success bound 0.39041.

Reading: the rows run from 0.4060 (q = 163, dp_bits 8, seeds 6000..6032, 825 runs, the lowest of the eight) to
0.4296 (q = 211, dp_bits 6, seeds 7000..7027); the spread across rows is within binomial noise (chi-square 2.34 on 7
degrees of freedom against the pooled fraction), and the dp_bits 8 rows (fresh starts at rate 2^-8 per step) are the noisiest. The
mean T_c/sqrt(N) is 1.239 to 1.272 against the model's 1.2734. The
toy fields' fresh starts (one per 2^dp_bits steps) collide at rate 1/N rather than 1/N', which favours the toy. The
pre-collision 2-cycle rate is 0.002518 to 0.002608 per step against the model 1/384 = 0.002604, and checkpoint
escapes (closed walks of length 3..64) occur at 1.39e-5 to 2.44e-5 per step against the bound 3.45e-5 of part F.
The realized waste per useful iteration is at most 0.0204 to 0.0239, 1.06 to 1.25 times the H2 rate 0.0191799
(rounded outward) and at most 0.59 of the tolerance 0.04075 of H1+ (proof.md 3.4 (c)). These toy-scale statistics,
with a hash in place of K's bits, do not measure production success; they support the threshold and waste part of
H1+ only, not any cap (the caps rest on Markov bounds, which a toy census cannot measure).

## What the experiments do not show

They do not run the production attack, do not measure production success, and check the kernels on samples, not
exhaustively. The toy fields have far fewer classes than GT and their fresh starts are far more frequent than in
production; they test the collision law of this class walk at their sizes, as H1 extrapolates it.

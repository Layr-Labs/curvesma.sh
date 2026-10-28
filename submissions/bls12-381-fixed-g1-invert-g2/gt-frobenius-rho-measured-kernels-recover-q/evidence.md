# Evidence: the experiments of `experiment.py` and the kernel checks

This file is section 7 of `proof.md`; section numbers refer to `proof.md`. `experiment.py` uses the Python standard
library only and reads its configuration from stdin; the declared run is
`echo '{"q":67,"trials":30,"R":32,"branch_bits":6,"sq_cut":0,"dp_bits":4,"walks":4,"seed":2,"production_checks":2,"production_steps":240,"production_dp_bits":4,"kernel_cases":2,"pair_samples":15000,"stat_classes":3000}' | python3 -I experiment.py`
(under 20 s on one core; about 7 s measured).

## Parts run on real BLS12-381 GT elements

* **Part A (identities):** GT membership; the Frobenius formula is the p-th power; the 12-key formula of F2 gives
  the keys of the 12 explicit conjugates, all distinct; kappa selects one canonical element from every conjugate; p^6 inverts; ord_r(p) = 12; 1 + p^4 + p^8 = 0; the r - 1 factorisation; ord_r(2) (log2 250.857,
  smallest t with 2^t in <p>: log2 247.272); the Granger-Scott formula equals y*y (9 Fp2 squarings, 34 additions);
  y y^(p^4) = y^(p^2) and the w^3 identity W(y) = -y_3 of F5 hold, and W(x) = -x_3 fails for the random x of Fp12
  (off GT) each y is made from.
* **Part K (the four walk kernels in the emulator):** `gen.py` and `gen_lazy.py` emit `mul.s`, `sqr.s`, `lazy.s`,
  `root.s` in memory; `kernels/emu.py` executes them (every instruction counts one; no branch but the final ret, so
  the executed count is the emitted count). On real GT elements x, m (products of random GT generators):
  `mul.s(x, entry(m)) = x m`, `sqr.s(x) = x^2`, `root.s(x)` returns the stored a, b of x, writes the per-base data
  (equal to `gen_lazy.block_vectors(x)`) and (K, kappa) equal to the reference minimum of the twelve stored keys,
  `lazy.s(x, per-base data, node entry of m)` returns the stored a, b and (K, kappa) of x m, in the first case also
  checked against keys from twelve explicit Frobenius maps; `mul.s` also runs in the trie-construction orientation
  (the multiplier's coordinates against the node's stored forms). Every call executes exactly 12,735 / 6,786 /
  3,419 / 1,619 instructions. Declared run: 0 mismatches. The node-setup kernel `node.s` runs in `check_lazy.py`.
* **Part C (the production lazy walk, with trails and replay):** R = 32 multipliers M_h = g^{a_h} z^{b_h} (40-bit
  exponents, for speed), the 384 twisted entries, trie nodes built on demand, depth L = 6, the 2-cycle escape rule and
  the production DP field (bits 64.. of K) at a test density 2^-production_dp_bits: at a distinguished class the
  record (F, alpha0, beta0), F = bits 150..277 of K, is logged and the walk takes a fresh start. Every key comes from
  the emulated `root.s` or `lazy.s`, every flush and escape from the emulated `mul.s` and `sqr.s`; an eager walk (keys
  from F2, exponents carried) runs beside it. Asserted at every
  step: base * Pi_node equals the eager element and (K, kappa) the reference; every kernel result equals the Python
  product; the exponent invariant. Every record is then replayed eagerly from its start pair alone to a distinguished
  class with the recorded F, and p^kappa (alpha, beta) must give the canonical representative g^alpha^ z^beta^. Declared run (240 steps, test
  density 1/16): 49 root steps, 191 lazy steps, 26 flushes, one 2-cycle escape whose smaller-key class needed a flush
  from an internal node, 22 records, all 22 replayed to the same F with the canonical pair verified (218 replayed
  steps), 0 key, element or kernel mismatches, 0 bases failing W(y) = -y_3; chi-square of the multiplier indices
  26.76 on 31 degrees of freedom. The JSON's sha256 of the 240 steps' (K, kappa), 2eb3944836520c4a..., is the one the
  kernels of 55afb7e5 give at the declared input (same seed, same walk): the keys are the same residues.
* **Part G (the bits of K on production classes):** 3,000 distinct classes (g^(c + i d), i < 3,000): chi-square of
  h over 32 bins 39.06 (31 dof), of bits 64..71 (the low DP-field byte) 278.85 (255 dof), of bits 150..157 (the low
  fingerprint byte) 230.04 (255 dof); the largest |z| of the 278 single-bit frequencies 3.43 (expected maximum of 278
  standard normals about 2.9: a mild excess on the maximum, not a dependence signal); every K below p/2. Marginal
  statistics only; H2's random-function model of the 278 bits is a heuristic (proof.md section 6).
* **`kernels/check_kernels.py`** (not run by the verifier; under 30 s): regenerates the five listings in a temporary
  directory, checks that every immediate of every emitted instruction encodes in AArch64 (ldp/stp -512..504 and
  ldr/str 0..32760, multiples of 8; add/sub/cmp 0..4095 with lsl 0 or 12; mov a 16-bit chunk or its inverse; extr
  0..63: 1,706 / 1,133 / 396 / 233 / 292 immediates in mul, sqr, lazy, root, node, all encodable), and runs `check.py`
  (mul.s: 40 GT products, 2 edge, 4 random; sqr.s: 20 cyclotomic elements) and `check_lazy.py` (lazy.s and root.s on
  4 GT, 2 random Fp12 and 2 edge-value cases, raw stored inputs at 0 and p - 1 with Ea, Ed at both ends of [2p, 3p)
  and the offset's own extreme, root.s at x0 = 0 and x4 extremes, and placed keys; lazy.s against the integer
  function it computes on every case and that function against explicit Frobenius keys on GT; node.s against the
  entry rule on every node entry): executed counts constant and equal to the charged ones (node.s 2,158), 0
  mismatches. After `gen.py` and `gen_lazy.py` write the listings, `check_lazy.py 8 11 mutants` (in kernels/) passes and
  catches all eleven broken listings or constants (top-limb-only minimum, `ls` for `lo`, a dropped `csel`, a flipped `b`
  select, unreduced pre-sums or U0 = p at x0 = 0 in root.s, an operand sum without its carry, `d_T` without the doubling
  of S0, Ed read for Ea, the operands regrouped, node offsets in [p, 2p)).

## Toy fields (q = 67, 139, 163, 211, 379)

* **Part B (complete algorithm):** Frobenius classes of the toy GT (order r' the largest prime factor of
  Phi_12(q)), keys and kappa, the twisted table with R = 32 and no squaring branch (sq_cut 0), 2-cycle and checkpoint
  escapes by squaring, every trail ending at its first DP with a record (F, alpha0, beta0) and a fresh start, the
  matcher of 0-ii (whole-record heapsort with its level count checked against the worst-case formula, adjacent dedup,
  pairs under a toy budget of 2^16, never binding: at most one pair charged per trial; both trails of every pair
  replayed from their start pairs by the walk's own iteration, checked to end at the recorded fingerprint within the
  iteration bound), halting at the first verified g^x = z, the pre-debit meter (toy budget 2^40 instructions), and
  the production stuck-walk counter with its credits against a toy cap: a trail restarts when the counter passes
  40 x 2^dp_bits; no report cap. **Credit check** (asserted): per trail, the set of classes visited; at every
  iteration the counter less the open episode's revisits (at most 64) is at most the trail's distinct classes; every
  credit covers its episode's revisits, after which the counter is at most the distinct classes (a trail with an
  uncredited revisit is counted, not checked). The toy walk is eager (part C shows the lazy walk visits the eager
  walk's elements); its class bits come from a 64-bit hash of the canonical representative (a toy key lies in F_q,
  far too few bits; part G tests the production bits). T_c is the first step at which a class is reached with a
  different canonical exponent pair; part B reports the raw fraction with T_c <= S(P) and the fraction whose charged
  count at T_c (defined below the census table) is within S(0.4005)(1 + eps), eps rounded up to 0.02451856. Declared
  run (q = 67, dp_bits 4, 30 trials): 30 logs recovered by the matcher and verified, fraction by the charged count
  0.4333 (standard error about 0.09), pre-collision 2-cycle rate 0.0036 (44 in 12,199 steps against about 32 at the
  model 1/384, about 2.2 standard deviations; the census below measures it at 0.002492 to 0.002616), heapsort levels
  at most 0.9066 of the worst-case bound, longest replay 141 iterations (bound 65 x 40 x 2^4 + 1 = 41,601), at most one
  pair charged per trial; credit check: 3,196 trails, none with an uncredited revisit, 148 escapes, at most 63 open
  revisits, distinct classes less counter at least 1 after every credit.
* **Part E (pair-collision law):** sample m distinct classes, apply one step F to each, count image coincidences;
  N times the pair-collision rate is what H1 models. Declared run (q = 67, R = 32, s = 0, 15,000 classes): 1.0447
  (relative standard error 0.0331); with the declared input but `"trials":0`, `"production_checks":0`,
  `"production_steps":0`, `"kernel_cases":1`, `"stat_classes":0` and `"pair_samples":40000`, 0.9786 (relative s.e.
  0.0128), and with `"pair_samples":60000` 0.9863 (relative s.e. 0.0085); against the
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
  iteration maximum, the caps' expectations and Markov ratios, the claim 123.826 above log2 W < 123.82557 (exact
  integer inequality), memory, the itemised losses (sum below 0.0100854) and success at least 0.39041 (floored).
* **Part S** exercises the bounded sampler of 0-iv, and `flags_check.py` runs the multiword compares and counters of
  0-i, 0-ii and 0-iv under the AArch64 NZCV rules, the stuck-walk counter's two credits and reset included.

## Offline toy census of this package's walk (R = 32, s = 0, 4 walks)

Each row pools seeds run with `experiment.py`'s `part_B(random.Random(seed), q, trials, 32, 6, 0, dp_bits, 4)` by an
offline driver (not shipped) that asserts this package's charges and the credited counter in the file it loads: 100
trials per seed (25 at dp_bits 8), every seed in each listed range run and counted. Every run recovered and verified its
logarithm through the matcher and replay.

**Registration.** One block of six rows was registered for this package's charges (block 32,487 instructions, flush step
16,342, escape 24,000, toy eps 0.02451856 = this ledger's eps rounded up, model point 0.4005) on this `experiment.py`,
before any of its seeds (41000.., unused by every earlier census of this walk) had run; its rows were read only after
all six finished. Metric: `successes_charged_steps_le_S_model_Pm_times_1_plus_eps` (Pm = 0.4005: charged steps at T_c
within S(0.4005)(1 + eps)). The rule as registered: pool all six rows, report each row and the pool with a binomial s.e.
rounded up, and the block passes when the pooled fraction less four binomial standard errors is at least 0.3905. No row
was dropped. An earlier registered block (seeds 40000.., the same rows and rule; pooled 4,716 of 11,175, 0.4220) ran on
the toy before the stuck-walk counter's credits were added and is superseded because the toy changed: it does not
reproduce on this `experiment.py` (q = 67, dp_bits 8, seed 40000: 9 successes there, 8 here).

**What the census counts.** A run's charged count at T_c (its first useful collision) is the toy's `ord_at_Tc + 4.43254
esc2 + 197.597 escL` (24,000/(c C) and 64 x 16,342/(c C) + 24,000/(c C) steps, to nearest): one step per iteration that
took a step, each 2-cycle escape at 24,000 instructions in steps of c C, each checkpoint escape at 64 steps of 16,342
plus 24,000. Against U + X at T_c (proof.md 3.3) it omits the distinguished-point steps (one per 2^dp_bits toy steps;
one per 2^86 in production), the detecting iteration (at most 0.0027 per counted step in every row) and, per length-64
escape, the meter's 4,160 bookkeeping instructions (0.7683 step; the meter charges 198.3651); X's 2 steps per 2-cycle
and 64 + l per longer cycle bound the steps it counts, so the count is easier than production's by under 1.9% of it at
dp_bits 6, 0.7% at 8. The column "waste per useful iteration, at most" charges every escape as X does, at those exact
charges (2 + 24,000/(c C) steps per 2-cycle, 128 + 1,074,048/(c C) per checkpoint escape), over the row's useful
iterations (counted steps less one per 2-cycle and 128 per checkpoint escape), rounded up. Rounding: success fractions
(bold) and the pool are floored to four decimals, the s.e. sqrt(p(1 - p)/runs) rounded up, the waste column up, the rest
to nearest (the means weight the per-seed means `experiment.py` prints to four decimals); ranges quoted below are
rounded outward, success fractions floored. Classes N: 129,145 (q = 67), 4,524,905 (163), 721,279 (211).

| q | dp_bits | seeds | runs | successes | mean T_c/sqrt(N) (model 1.2734) | mean charged steps/sqrt(N) | at the model point S(0.4005)(1+eps) | s.e. | 2-cycles per step before T_c (model 0.002604) | checkpoint escapes per step before T_c | waste per useful iteration, at most |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 163 | 6 | 41000..41020 | 2100 | 885 | 1.2541 | 1.2503 | **0.4214** | 0.0108 | 0.002557 | 1.96e-05 | 0.0230 |
| 211 | 6 | 41000..41032 | 3300 | 1405 | 1.2585 | 1.2537 | **0.4257** | 0.0087 | 0.002551 | 1.67e-05 | 0.0220 |
| 67 | 6 | 41000..41032 | 3300 | 1417 | 1.2643 | 1.2581 | **0.4293** | 0.0087 | 0.002615 | 1.09e-05 | 0.0205 |
| 163 | 8 | 41000..41032 | 825 | 339 | 1.2709 | 1.2822 | **0.4109** | 0.0172 | 0.002582 | 2.03e-05 | 0.0234 |
| 211 | 8 | 41000..41032 | 825 | 343 | 1.2676 | 1.2784 | **0.4157** | 0.0172 | 0.002582 | 1.93e-05 | 0.0231 |
| 67 | 8 | 41000..41032 | 825 | 328 | 1.2734 | 1.2832 | **0.3975** | 0.0171 | 0.002493 | 1.60e-05 | 0.0214 |

Pooled: **11,175 runs, 4,717 successes, fraction 0.4221 +- 0.00468** at the model point S(0.4005)(1 + eps); four
binomial standard errors below it, **0.4034** (0.403415, floored), above the registered bar 0.3905 and the package's
floored success bound 0.39041.

Reading: the rows run from 0.3975 (q = 67, dp_bits 8, 825 runs, the lowest) to 0.4293 (q = 67, dp_bits 6); the spread
is within binomial noise (chi-square 3.50 on 5 degrees of freedom against the pooled fraction), and the dp_bits 8 rows
(fresh starts at rate 2^-8 per step) are the noisiest. The mean T_c/sqrt(N) is 1.254 to 1.274 against the model's
1.2734. The toy fields' fresh starts collide at rate 1/N rather than 1/N', which favours the toy. The pre-collision
2-cycle rate is 0.002492 to 0.002616 per step against the model 1/384 = 0.002604, and checkpoint escapes (closed
walks of length 3..64) occur at 1.08e-5 to 2.04e-5 per step against the bound 3.45e-5 of part F. The realized waste per
useful iteration is at most 0.0204 to 0.0234, 1.03 to 1.19 times the H2 rate 0.0197654 (rounded outward) and at most
0.57 of the tolerance 0.04149 of H1+ (proof.md 3.4 (c)). These toy-scale statistics, with a hash in place of K's bits,
do not measure production success; they support the threshold and waste part of H1+ only, not any cap (the caps rest
on Markov bounds, which a toy census cannot measure).

## What the experiments do not show

They do not run the production attack, do not measure production success, and check the kernels on samples, not
exhaustively. The toy fields have far fewer classes than GT and their fresh starts are far more frequent than in
production; they test the collision law of this class walk at their sizes, as H1 extrapolates it.

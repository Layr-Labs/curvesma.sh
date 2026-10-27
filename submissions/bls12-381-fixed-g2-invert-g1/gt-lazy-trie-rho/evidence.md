# Evidence: the experiments of `experiment.py` and the kernel checks

This file is section 7 of `proof.md`; section numbers refer to `proof.md`. `experiment.py` uses the Python standard
library only and reads its configuration from stdin; the declared run is
`echo '{"q":67,"trials":30,"R":24,"branch_bits":6,"sq_cut":0,"dp_bits":4,"walks":4,"seed":1,"production_checks":2,"production_steps":180,"kernel_cases":2,"pair_samples":15000}' | python3 -I experiment.py`.

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
* **Part C (the production lazy walk):** R = 24 multipliers M_h = g^{a_h} z^{b_h} (40-bit exponents, for speed), the
  288 twisted entries, trie nodes built on demand (Pi_child = Pi_parent T[h][t]), depth L = 6, the 2-cycle escape rule,
  no DP test. Every key comes from the emulated `root.s` (at the root) or `lazy.s` (elsewhere), every flush and every
  escape's explicit representative from the emulated `mul.s`, every escape square from the emulated `sqr.s`. An eager
  walk (explicit element multiplied by T[h][t] each step, keys from F2) runs beside it. Asserted at every step:
  base * Pi_node equals the eager element, and (K, kappa) equals the reference; after every flush and escape the
  kernel result equals the Python product; at the end y = g^alpha z^beta for the eagerly tracked exponents. The declared
  run (180 steps) read: 31 root steps, 149 lazy steps, 29 flushes, one 2-cycle escape (at a root, no flush needed),
  0 key, element or kernel mismatches, exponent invariant true; chi-square of the 24 multiplier indices 29.49 on 23
  degrees of freedom (about 7 per bin: a small sample). A 240-step run (seed 1, other inputs as declared but 40 toy
  trials) read 41 root, 199 lazy, 39 flushes, one escape, 0 mismatches.
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
  T_c is the first step at which a class is reached with a different canonical exponent pair. Besides the raw
  fraction with T_c <= S(0.4005), part B reports the fraction in the ledger's own units: walk steps up to T_c (wasted
  2-cycle steps included) plus every 2-cycle escape at 3.4432 steps and every checkpoint escape at its worst case,
  against S(0.4005)(1 + eps) with the production eps = 0.029892; distinguished-point and fresh-start iterations are
  not counted, since production has one per 2^72 steps. Declared run (q = 67, dp_bits 4, 30 trials): 30 logs
  recovered by the matcher and verified, raw fraction within S(0.4005) 0.567, in the ledger's units 0.600 (30 trials:
  standard error about 0.09), pre-collision 2-cycle rate 0.0040 (45 escapes in 11,252 steps), heapsort levels at most
  0.914 of the worst-case bound, longest replay 87 iterations (bound 1,285), at most one pair charged per trial.
* **Part E (pair-collision law):** sample m distinct classes, apply one step F to each, count image coincidences;
  N times the pair-collision rate is what H1 models. Declared run (q = 67, R = 24, s = 0, 15,000 classes): 1.0126
  +- 0.0337 (a 20,000-class run: 0.9751 +- 0.0257) against the exact first-order value 1 - 1/(12R) = 0.9965 and the conservative value 1 - 1/R = 0.9583
  that section 4 uses.
* **Part F (fruitless closed walks at R = 24):** q(n), the probability that n independent uniform 12th roots of unity
  sum to zero in Z[zeta_12], computed exactly over the 12^n sequences for n <= 12 (q(2) = 1/12, q(3) = 1/72,
  q(4) = 0.019097, q(6) = 0.0077562, q(8) = 0.0042930, q(10) = 0.0027710, q(12) = 0.0019476), and Q(12) = 0.0019476
  (the largest point probability of a 12-term sum); the multinomial bound of section 3.3 exactly for l <= 288 (P_2 =
  1/288 exactly, P_3 = 2^-15.34, P_4 = 2^-14.76, P_8 = 2^-25.94, P_16 = 2^-43.97, P_32 = 2^-71.04, P_64 = 2^-108.88,
  sum_{3..64} = 6.1659e-5, sum_{65..288} = 2^-108.79) and the Poissonised tail beyond 288 weighted by the probability of
  holding no distinguished class (2^-105.80): useless walks longer than 64 at a rate below 2^-105.63 per step.
* **Part D (ledger):** recomputes, in 60-digit decimal arithmetic, c = 41,821/120,000, the step prices from the
  itemised lists, S, eps from part F, the trie size and its setup charge, the report, fresh-start, restart, switch and
  stopping-allowance terms, the claim 124.206 as the least thousandth above log2 W = 124.20529 (exact integer
  inequality), the walk budget and threshold on integers (the walk need fits with 0.049% to spare), N_max = 2^127.79,
  the expected restarts (stuck 2^26.94, long useless walks 2^23.79), memory 2^61.92 bytes, the success chain
  (losses 0.010023, success at least 0.39048) and sensitivities.
* **Part S** exercises the bounded sampler of 0-iv, and `flags_check.py` runs the multiword compares and counters of
  0-i, 0-ii and 0-iv under the AArch64 NZCV rules (the parent's listings; the parent's squaring-branch compare is
  removed, since this walk has no squaring branch).

## Offline toy censuses of this package's walk (R = 24, s = 0, 4 walks)

Each row pools seeds run with `experiment.py`'s `part_B` (registered seed ranges, every run counted):

| q | dp_bits | classes N | seeds | runs | mean T_c/sqrt(N) (model 1.2803) | mean charged steps/sqrt(N) | raw fraction <= S(0.4005) | fraction in ledger units <= S(0.4005)(1+eps) | s.e. | 2-cycles per step before T_c (model 0.003472) | checkpoint escapes per step before T_c |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 163 | 6 | 4524905 | 2000..2018 | 1900 | 1.2659 | 1.2651 | 0.3963 | **0.4179** | 0.0113 | 0.00341 | 3.71e-05 |
| 163 | 8 | 4524905 | 3000..3029 | 750 | 1.2396 | 1.2514 | 0.4147 | **0.4280** | 0.0181 | 0.00346 | 2.80e-05 |
| 211 | 6 | 721279 | 2000..2024 | 2500 | 1.2552 | 1.2516 | 0.4212 | **0.4392** | 0.0099 | 0.00336 | 2.64e-05 |
| 67 | 6 | 129145 | 2000..2026 | 2700 | 1.2658 | 1.2610 | 0.4096 | **0.4356** | 0.0095 | 0.00329 | 2.16e-05 |

Pooled over all rows: 7850 runs, fraction in ledger units 0.4317 +- 0.0056.

Settings: R = 24, sq_cut 0 (no squaring branch), branch_bits 6, 4 walks, 100 trials per seed (25 at dp_bits 8), every
seed in each listed range run and counted, an offline driver (not shipped) calling `part_B(random.Random(seed), q, trials, 24, 6, 0,
dp_bits, 4)` of this `experiment.py`. Every run recovered and verified its logarithm through the matcher and replay.
Reading: the fraction of runs whose first useful collision falls within S(0.4005)(1 + eps), counted in the ledger's
units, is at or above the model's 0.4005 in every row (pooled 0.432); the raw fraction (every toy iteration counted
as a step, including the toy's frequent distinguished-point iterations) is 0.396 to 0.422. The toy fields' fresh
starts (one per 2^dp_bits steps) collide at rate 1/N rather than 1/N', which favours the toy slightly; the dp_bits 8
row (0.4% fresh starts) reads 0.428 +- 0.018. The pre-collision 2-cycle rate matches the model 1/288 within 6%
(below it), and checkpoint escapes (closed walks of length 3..64) occur at 2-4e-5 per step against the bound
6.17e-5 of part F, so the 1.25 margin on eps is not consumed at these sizes. These are toy-scale statistics; they do
not measure production success.

## What the experiments do not show

They do not run the production attack, do not measure production success, and check the kernels on samples, not
exhaustively. The toy fields have far fewer classes than GT and their fresh starts are far more frequent than in
production; they test the collision law of this class walk at their sizes, as H1 extrapolates it.

# Evidence: the experiments of `experiment.py` and the kernel checks

This file is section 7 of `proof.md`; section numbers refer to `proof.md`. `experiment.py` uses the Python standard
library only and reads its configuration from stdin; the declared run is
`echo '{"q":67,"trials":30,"R":24,"branch_bits":6,"sq_cut":0,"dp_bits":4,"walks":4,"seed":2,"production_checks":2,"production_steps":240,"production_dp_bits":4,"kernel_cases":2,"pair_samples":15000,"stat_classes":3000,"r2_steps":12000}' | python3 -I experiment.py`
(about 70 s on one core).

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
  L = 6, the 2-cycle escape rule and the doubling (Brent) checkpoint, and the production distinguished-point field
  (bits 64.. of K) read at a test density 2^-production_dp_bits: at a distinguished class the record (F, alpha0,
  beta0), F = bits 150..261 of K, is logged and the walk takes a fresh start from new start exponents. Every key comes
  from the emulated `root.s` (at the root) or `lazy.s` (elsewhere), every flush and every escape's explicit
  representative from the emulated `mul.s`, every escape square from the emulated `sqr.s`. An eager walk (explicit
  element multiplied by T[h][t] each step, keys from F2, exponents carried) runs beside it. Asserted at every step:
  base * Pi_node equals the eager element and (K, kappa) equals the reference; after every flush and escape the
  kernel result equals the Python product; the last trail's exponent invariant. Afterwards every record is replayed
  eagerly from its start pair alone (the walk rule only, no stored state): the replay must reach a distinguished
  class with the recorded F, and the canonical pair p^kappa (alpha, beta) must give the canonical representative
  g^alpha^ z^beta^. The doubling checkpoint's detection time is checked against 2 max(mu, l) + l <= 3 (mu + l) for
  every rho shape with mu + l <= 160 (largest T/(mu + l) 2.977). Declared run (240 steps, test density 1/16): 47 root
  steps, 193 lazy steps, 32 flushes, one 2-cycle escape whose smaller-key class needed a flush from an internal node,
  no doubling escape, 13 records, all 13 replayed to the same F with the canonical pair verified (223 replayed steps),
  0 key, element or kernel mismatches; chi-square of the multiplier indices 25.26 on 23 degrees of freedom. A second
  run of 12,000 steps at R = 2 (seed + 1, density 1/64), where long closures are common: 2,288 root and 9,712 lazy
  steps, 1,609 flushes, 440 two-cycle, 7 sixty-four-step and 63 doubling-checkpoint escapes (437 escape flushes),
  168 records all replayed and verified (11,784 replayed steps), chi-square of the multiplier index 0.43 on 1
  degree of freedom, 0 mismatches. Part C carries the 2-cycle test, the 64-step checkpoint and the doubling
  checkpoint, and this run fires all three on the production walk.
* **Part G (the bits of K on production classes):** 3,000 distinct classes (g^(c + i d), i < 3,000): chi-square of
  h over 24 bins 21.71 (23 dof), of bits 64..71 (the low DP-field byte) 237.21 (255 dof), of bits 150..157 (the low
  fingerprint byte) 284.31 (255 dof, upper tail about 0.10); the largest |z| of the 262 single-bit frequencies 2.67
  (expected maximum of 262 standard normals about 2.9); every K below p/2. These test marginal statistics only; H2's
  random-function model of the 262 bits is a heuristic (proof.md section 6).
* **`kernels/check_kernels.py`** (not run by the verifier; about 30 s): regenerates the four listings in a temporary
  directory and runs `check.py` (mul.s: 40 GT products, 2 edge, 4 random; sqr.s: 20 cyclotomic elements) and
  `check_lazy.py` (lazy.s and root.s on 4 GT, 2 random Fp12 and 2 edge-value cases, references from explicit
  Frobenius keys): all executed counts constant and equal to the charged ones, 0 mismatches. `check_lazy.py 24 7`
  (12 GT, 6 random, 6 edge cases, inputs drawn from {0, 1, p - 1, p - 2, random}) also passes.

## Toy fields (q = 67, 139, 163, 211, 379)

* **Part B (complete algorithm):** Frobenius classes of the toy GT (order r' the largest prime factor of
  Phi_12(q)), keys and kappa, the twisted table with R = 24 and no squaring branch (sq_cut 0), 2-cycle, 64-step
  checkpoint and doubling checkpoint escapes by squaring, the stuck-walk rule at 120 << dp_bits counted steps, every
  trail ending at its first DP with a record (F, alpha0, beta0) and a fresh start, the matcher of 0-ii (whole-record
  heapsort with its level count checked against the worst-case formula, adjacent dedup, pairs under the 2^4 budget,
  both trails of every pair replayed from their start pairs by the walk's own iteration with a check that the replay
  ends at the recorded fingerprint and within the iteration bound), halting at the first verified candidate,
  g^x = z, the pre-debit meter with a toy budget. The toy walk is eager: part C shows the lazy walk visits exactly
  the eager walk's elements, so the toy tests the walk's collision behaviour, not the kernels. The toy's class bits
  come from a 64-bit hash of the canonical representative (a toy key lies in F_q, q <= 379, far too few bits); part G
  tests the production bits. T_c is the first step at which a class is reached with a different canonical exponent
  pair. Part B reports the useful steps U and the waste X before T_c in ledger units (as defined at `TOY_ESC_STEPS`),
  X/(W_BAR U), and the fraction of runs with U + X <= S(0.4014)(1 + 0.030293); distinguished-point and fresh-start
  iterations are not counted, since production has one per 2^86 steps. Declared run (q = 67, dp_bits 4, 30
  trials): 30 logs recovered by the matcher and verified, 8 of 30 (0.267, binomial standard error about 0.08) within
  the meter's walk units, mean T_c/sqrt(N) 1.400 against the model 1.280, waste 0.999 of the H2 rate before T_c,
  pre-collision 2-cycle rate 0.0033 (46 escapes in 14,074 steps), 0 stuck restarts, 0 meter stops, heapsort levels
  at most 0.921 of the worst-case bound, longest replay 112 iterations (bound 3,845), at most one pair charged per
  trial. Thirty trials do not measure success; censuses 3 and 4 below do, on 3,600 registered runs each.
* **Part E (pair-collision law):** sample m distinct classes, apply one step F to each, count image coincidences;
  N times the pair-collision rate is what H1 models. Declared run (q = 67, R = 24, s = 0, 15,000 classes): 1.0103
  +- 0.0337 against the first-order value 1 - 1/(12R) = 0.9965 and the conservative value 1 - 1/R = 0.9583
  that section 4 uses.
* **Part F (fruitless closed walks at R = 24; code in `cycles.py`, loaded by `experiment.py`):** q(n), the probability that n independent uniform 12th roots of unity
  sum to zero in Z[zeta_12], computed exactly over the 12^n sequences for n <= 12 (q(2) = 1/12, q(3) = 1/72,
  q(4) = 0.019097, q(6) = 0.0077562, q(8) = 0.0042930, q(10) = 0.0027710, q(12) = 0.0019476) and by convolution for
  13 <= n <= 24 (q(16) = 0.0011164, q(24) = 0.00050635); g(n) = q(24) beyond 24 by the Cauchy-Schwarz point bound
  (`q_2m0_bound_beyond`); the multinomial bound of section 3.3 exactly for l <= 600 (P_2 <= 1/288, P_3 = 2^-15.34,
  P_4 = 2^-14.76, P_8 = 2^-25.94, P_16 = 2^-43.97, P_32 = 2^-71.04, P_64 = 2^-108.88, P_288 = 2^-213.07,
  P_600 = 2^-257.58, sum_{3..64} = 6.1659e-5, sum_{65..600} = 2^-108.79; polynomial coefficients from index count
  403 up underflow and are replaced by a union bound, `underflow_cut_j`) and the Poissonised tail beyond 600 weighted
  by the probability of holding no distinguished class at theta = 2^-86 (2^-131.44): useless walks longer than 64 at a
  rate below 2^-108.79 per step, now caught by the doubling checkpoint and charged in eps (part D).
* **Part D (ledger):** recomputes, in 60-digit decimal arithmetic, c = 41,821/120,000, the step prices from the
  itemised lists (every step 170 instructions), S, eps = 1.25 (eps2 + eps_long + eps_brent) = 0.030294 from part F,
  the trie size and its setup charge, the report, fresh-start, restart, switch and stopping-allowance terms at the
  model point P = 0.4014 (log2 S = 125.68551), the claim 124.209 as the least thousandth above log2 W = 124.20834
  (exact integer inequality), the walk budget and threshold on integers (meter slack 4.56e-4 of the need), the block
  inequalities behind the iteration maximum n_max < 2^125.730, the detection lemma and the Markov caps of proof.md
  3.5 (P(stuck) 2^-43.28; expected reports 2^39.73, tail 2.53e-5; expected restarts 2^-3.55, tail 8.1e-8; false
  pairs 2^-18.27, tail 2.1e-7), the H1+ tolerance loss (t = 2, S_tol = 2^125.66028, loss 0.0106518), the restart-cap
  worst case as a fraction of S (1.4232e-4), memory 2^61.92 bytes, the itemised success losses (sum 0.0108207),
  success at least 0.39057 (floored from 0.3905793), the H1+ sensitivity of that bound (it falls 0.0136 per unit of the
  tolerance t from t = 2 to 3, the largest t in hundredths still reaching 0.39 is 2.04, and at t = 3 it is 0.37694) and
  sensitivities (-0.0085 bits if the eps margin were 1, +0.0205 bits per 100 more instructions on every step, +0.560
  bits if the lazy kernel cost doubled).
* **Part S** exercises the bounded sampler of 0-iv, and `flags_check.py` runs the multiword compares and counters of
  0-i, 0-ii and 0-iv under the AArch64 NZCV rules (the parent's listings; the parent's squaring-branch compare is
  removed, since this walk has no squaring branch).

## Census 3 (registered before any count was read)

Protocol, written into this file before `census.py` was first run on these seeds (the only earlier run was a
4-trial timing pilot on seed 1, whose counts were not printed):

* Walk: this package's toy part B (Brent doubling checkpoint beside the 64-step checkpoint, stuck cap 120 << dp_bits,
  pair budget 2^4), R = 24, branch_bits 6, s = 0, 4 walks.
* Fields and density: q = 67, 163, 211, each at dp_bits 6.
* Seeds: 7000..7019 for every field, 60 trials per seed, so 1,200 runs per field and 3,600 in all. Every seed is run
  and counted; a run the toy meter stops counts as a failure.
* Statistic 1 (waste): X / (W_BAR U) before T_c, where U and X are the useful steps and the waste in ledger units as
  defined at `TOY_ESC_STEPS` in `experiment.py` (2-cycle escape: 2 + 3.4432 steps; checkpoint escape of length l:
  64 + l walked steps, l traversal steps at cmax (2.8054 steps each) and the escape; doubling escape: every step of
  its segment, l traversal steps at cmax and the escape), W_BAR = 0.024235 the H2 rate of part D rounded down. Interval: mean of the 20 per-seed
  ratios, Student t with 19 degrees of freedom, 95%; pooled over the 60 seed-field ratios with the same t factor.
* Statistic 2 (success): the fraction of runs whose charged steps U + X at T_c are at most S(0.40082)(1 + 0.030293),
  the production meter's walk units at the model point of the ledger at the time (0.40082), with a Wilson 95% interval.
* Pass rule: the upper end of every waste interval (each field and pooled) is below t = 2, the H1+ tolerance of
  proof.md section 4, and the pooled Wilson lower end of the success fraction is at least 0.39. A failure is reported
  as a failure; the protocol is not changed after the counts are read.
* Command: `nice python3 -I census.py 60` (estimated beforehand at about 10 minutes; it took 1,542 s on one core).

**Result: PASS** (`census.py` output, 3,600 runs, 0 meter stops; every run recovered and verified its logarithm):

| q | runs | success fraction | Wilson 95% | waste / H2 rate | t 95% interval |
| ---: | ---: | ---: | --- | ---: | --- |
| 67 | 1,200 | 0.4392 | [0.4113, 0.4674] | 0.8209 | [0.787, 0.8548] |
| 163 | 1,200 | 0.4267 | [0.3990, 0.4548] | 0.8742 | [0.8441, 0.9044] |
| 211 | 1,200 | 0.4433 | [0.4155, 0.4716] | 0.8736 | [0.8320, 0.9152] |
| pooled | 3,600 | 0.4364 | [0.4203, 0.4526] | 0.8562 | [0.8350, 0.8775] |

Both pass items hold: every waste interval ends below t = 2 (the largest upper end is 0.9152), and the pooled
success lower end 0.4203 is above 0.39. The waste reads 0.82 to 0.87 of the H2 rate, so the toy does not use the
factor 2 of H1+; this is toy-scale evidence for H1+, not a tail bound at production size.

## Census 4 (census 3 rerun at the current model point)

After census 3 the ledger's model point moved from 0.40082 to 0.4014 (the H1+ margin was widened, part D), so
statistic 2 of the registered protocol, whose threshold is the meter's walk units at the model point, had to be
read again. Census 4 is census 3 with one change, the threshold S(0.4014)(1 + 0.030293) in place of
S(0.40082)(1 + 0.030293) (0.030293 is `TOY_EPS` in `experiment.py`, 1.25 W_BAR rounded down where the ledger's
eps rounds it up to 0.030294, so the toy threshold is the tighter one); walk, fields, density, seeds, trials,
statistic 1, the intervals and the pass rule are the registered ones above. The change was fixed before the run
was started and this section was written after it finished. Command: `nice python3 -I census.py 60` (643 s on one
core, with nothing else running).

**Result: PASS** (`census.py` output, 3,600 runs, 0 meter stops; every run recovered and verified its logarithm):

| q | runs | success fraction | Wilson 95% | waste / H2 rate | t 95% interval |
| ---: | ---: | ---: | --- | ---: | --- |
| 67 | 1,200 | 0.4392 | [0.4113, 0.4674] | 0.8209 | [0.787, 0.8548] |
| 163 | 1,200 | 0.4267 | [0.3990, 0.4548] | 0.8742 | [0.8441, 0.9044] |
| 211 | 1,200 | 0.4442 | [0.4163, 0.4724] | 0.8736 | [0.8320, 0.9152] |
| pooled | 3,600 | 0.4367 | [0.4205, 0.4529] | 0.8562 | [0.8350, 0.8775] |

The seeds are the same, so the runs are the same walks and the waste columns reproduce census 3 exactly (a check
that the rerun is the registered census); only the success count moves, by one run in q = 211 (532 to 533 of 1,200) that falls
between the two thresholds. Both pass items hold: every waste interval ends below t = 2 and the pooled success lower end 0.4205
is above 0.39. proof.md cites census 4.

## Earlier offline toy censuses (package revision v3: no doubling checkpoint, theta 2^-72, pair budget 2^16)

These are kept as history; they ran the earlier walk and its earlier model point and are not evidence for this
revision's ledger.

Each row pools seeds run with `experiment.py`'s `part_B(random.Random(seed), q, trials, 24, 6, 0, dp_bits, 4)` by an
offline driver (not shipped): 100 trials per seed (25 at dp_bits 8), every seed in each listed range run and counted.
Every run recovered and verified its logarithm through the matcher and replay. "Fraction in ledger units" is the
fraction of runs whose first useful collision falls within S(P)(1 + eps), counted as the ledger counts (walk steps
including wasted 2-cycle steps, each 2-cycle escape at 3.4432 steps, each checkpoint escape at its worst case; no
distinguished-point or fresh-start iterations).

**Census 2 (v3 model point P = 0.3910):**

| q | dp_bits | classes N | seeds | runs | mean T_c/sqrt(N) (model 1.2803) | mean charged steps/sqrt(N) | fraction in ledger units <= S(0.3910)(1+eps) | s.e. | same at S(0.4005)(1+eps) | 2-cycles per step before T_c (model 0.003472) | checkpoint escapes per step before T_c |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 163 | 6 | 4524905 | 5000..5020 | 2100 | 1.2731 | 1.2713 | **0.4133** | 0.0107 | 0.4252 | 0.00339 | 3.28e-05 |
| 163 | 8 | 4524905 | 6000..6032 | 825 | 1.2624 | 1.2748 | **0.4097** | 0.0171 | 0.4206 | 0.00341 | 2.96e-05 |
| 211 | 6 | 721279 | 5000..5027 | 2800 | 1.2496 | 1.2452 | **0.4221** | 0.0093 | 0.4307 | 0.00337 | 2.26e-05 |
| 67 | 6 | 129145 | 5000..5029 | 3000 | 1.2423 | 1.2373 | **0.4273** | 0.0090 | 0.4407 | 0.00335 | 2.06e-05 |

Pooled over all rows: 8725 runs, fraction in ledger units 0.4206 +- 0.0053.

**Census 1 (v2 walk, still earlier; metric at P = 0.4005, so its last-but-two column repeats the
fraction):**

| q | dp_bits | classes N | seeds | runs | mean T_c/sqrt(N) (model 1.2803) | mean charged steps/sqrt(N) | fraction in ledger units <= S(0.4005)(1+eps) | s.e. | same at S(0.4005)(1+eps) | 2-cycles per step before T_c (model 0.003472) | checkpoint escapes per step before T_c |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 163 | 6 | 4524905 | 2000..2042 | 4300 | 1.2605 | 1.2592 | **0.4270** | 0.0075 | 0.4270 | 0.00341 | 3.47e-05 |
| 163 | 8 | 4524905 | 3000..3065 | 1650 | 1.2316 | 1.2448 | **0.4345** | 0.0122 | 0.4345 | 0.00346 | 3.43e-05 |
| 211 | 6 | 721279 | 2000..2054 | 5500 | 1.2501 | 1.2476 | **0.4398** | 0.0067 | 0.4398 | 0.00335 | 3.11e-05 |
| 67 | 6 | 129145 | 2000..2058 | 5900 | 1.2611 | 1.2568 | **0.4273** | 0.0064 | 0.4273 | 0.00334 | 2.33e-05 |

Pooled over all rows: 17350 runs, fraction in ledger units 0.4319 +- 0.0038.

Reading (v3): at S(0.3910)(1 + eps) every row of census 2 reads 0.410 to 0.427 (pooled 0.4206 +- 0.0053), and
census 1 reads 0.427 to 0.440 at S(0.4005)(1 + eps). The mean T_c/sqrt(N) is 1.23 to 1.27 against the model's
1.2803. The pre-collision 2-cycle rate is 0.00334 to 0.00346 per step against the model 1/288 = 0.00347, and
checkpoint escapes (closed walks of length 3..64) occur at 2.1e-5 to 3.5e-5 per step against the bound 6.17e-5 of
part F. These are toy-scale statistics with a hash in place of K's bits; they do not measure production success.

## What the experiments do not show

They do not run the production attack, do not measure production success, and check the kernels on samples, not
on every input. The toy fields have far fewer classes than GT and their fresh starts are far more frequent than in
production; they test the collision law of this class walk at their sizes, as H1 extrapolates it.

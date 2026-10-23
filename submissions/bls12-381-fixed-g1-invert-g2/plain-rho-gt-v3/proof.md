# Fixed-G1 pairing inversion via a bounded distinguished-point rho in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 127.80** (639/5) target-group-operation
  equivalents. Unit: the dense Fp12 multiplication at C = 20,000 scalar AArch64
  instruction equivalents.
* Success probability: at least 0.39 under the textbook random-mapping heuristic
  for rho, supported by experiments.

This submission presents a **fixed-budget (hard-bounded) distinguished-point
Pollard rho** for the GT discrete logarithm. It relies only on the random-mapping
model for the rho collision time, which the review has stated can qualify. It
supersedes my plain-rho-gt (v1/v2), repairing the review's remaining points:

* It does **not** claim "no fruitless cycles". There is no automorphism folding,
  so the structural fruitless cycles of orbit walks are absent; but a
  random-mapping trajectory can, exceptionally, enter a cycle containing no
  distinguished point. These exceptional events are negligible and are handled by
  per-walk truncation, and the success analysis allows exceptional unsuccessful
  trajectories.
* The attack has an **explicit global charged-work cutoff** W_max: it performs at
  most W_max group multiplications in total and then stops. All walk, restart,
  truncated and exceptional-path work is inside W_max. W_max is the charged work.
* It does not assert 1/theta as a hard lag bound; detection and truncation are
  bounded by high-probability (not expected-value) arguments.

## 1. Reduction

Input: nonidentity P in G1 (order r), z in GT. Q0 the public G2 generator.

1. Validate P and z per the profile. 2. If z = 1, output Q = infinity.
3. Compute g = e(P, Q0) with the pinned map (`pairing(Q0, P)`); nondegeneracy and
   prime order make g a generator of GT (prime order r).
4. Find a in [0, r) with z = g^a (section 2). 5. Output Q = [a]Q0; then
   e(P, Q) = g^a = z, Q in G2; optionally reverify all 12 coefficients.

P is never replaced — the admissible route in `rules.md`.

## 2. Bounded distinguished-point rho for z = g^a

GT is cyclic of prime order r. The walk is a standard R-adding walk with
R = 2^20 entries M_j = g^{alpha_j} z^{beta_j}: state u = g^A z^B, step
u → u·M_{H(u)} (one dense Fp12 multiplication), A,B updated mod r. There is no
automorphism, orbit quotient, or look-ahead.

### 2.1 Hard global cutoff

Fix **W_max = ceil(1.25·k0)** total group multiplications, where
k0 = sqrt(2 r ln(1/0.61)) (section 3). The attack runs walks until either a usable
collision is found or W_max multiplications have been performed in total (summed
over all processors), then stops. Every multiplication — walk steps, restarts,
truncated segments, and any exceptional path — is counted against W_max, so the
total charged work never exceeds W_max·(per-step units). This is a hard bound, not
an expected value.

### 2.2 Distinguished points, restarts, per-walk truncation

Run w ≤ 2^40 walks from random starts. A point is **distinguished** when 30
designated hash bits are zero (theta = 2^-30). Store (256-bit fingerprint, A, B),
96 bytes, keyed by fingerprint; on a match solve a = (A - A')/(B' - B) mod r when
B != B'. After each distinguished point the walk restarts cheaply at D·M_j for a
random j (one multiplication; D the just-found distinguished point).

**Per-walk truncation.** A walk that runs seg_cap·(1/theta) = 64·2^30 = 2^36
steps without producing a distinguished point is abandoned and restarted. A
healthy walk produces a distinguished point after a geometric number of steps of
mean 1/theta, so it is truncated only with probability (1-theta)^{2^36}
= (1-2^-30)^{2^36} ≈ e^{-64} ≈ 2^-92. The rare truncated walks are exactly the
exceptional trajectories that fall into a distinguished-point-free cycle; each
wastes at most 2^36 steps, all inside W_max.

### 2.3 Exceptional trajectories

With no automorphism, the only cycles are rho tails. Under the random-mapping
model a trajectory's rho cycle has expected length Theta(sqrt(r)) and so contains
Theta(theta·sqrt(r)) = Theta(2^97) distinguished points; the probability that a
trajectory's cycle contains **no** distinguished point (e.g. a non-distinguished
fixed point, or a very short cycle) is negligible — bounded by the per-walk
truncation probability ~2^-92 above. Such a trajectory is truncated within 2^36
steps and contributes no collision; the success analysis (section 3) counts it as
an exceptional unsuccessful trajectory.

## 3. Success probability within the hard cutoff

Model the R-adding walk as a random mapping on the r elements. Set

    k0 = sqrt(2 r ln(1/0.61)),   log2 k0 = 127.420269,   W_max = ceil(1.25·k0).

Over T = W_max walk-points, a self-collision exists with probability
≈ 1 - exp(-T^2/(2r)) = 1 - 0.61^{1.5625} = 0.538.

The bounded algorithm succeeds unless one of the following exceptional events
removes the collision within W_max:

* the colliding trajectory is truncated before its detecting distinguished point
  (probability ≤ 2^-92 per walk, section 2.2); or
* the detecting distinguished point of the colliding pair falls after the W_max
  cutoff. The lag from the merge to the next distinguished point is geometric of
  mean 1/theta; it exceeds 2^36 steps with probability ≈ e^{-64} ≈ 2^-92, and
  2^36 ≪ W_max, so this removes at most a ≈ 2^36/W_max ≈ 2^-91 fraction of the
  existing collisions.

Hence the bounded success probability is at least 0.538 − 2^-90 ≥ **0.39**, with
large margin. The 25% budget factor (1.25) over k0 both lifts the collision
probability to 0.538 and absorbs a walk up to ~25% less efficient than a random
mapping.

### Experimental support (toy scale, Python 3.12, one core)

`experiment.py` is a faithful bounded implementation: it runs the parallel rho
with the hard global cutoff W_max = ceil(1.25·k0) and per-walk truncation, counts
every multiplication (walk, restart, truncated), and reports the fraction of
trials in which a usable collision is found **within** W_max — a Monte-Carlo
estimate of the bounded success probability. Group: the prime-order-ell subgroup
of F_q*.

| bits | seeds | success within W_max | mean truncations/trial |
| ---: | --- | ---: | ---: |
| 26 | 1–10 | 0.450–0.573 (mean ≈ 0.52) | 0.000 |
| 28 | 5 | 0.525 | 0.000 |
| 30 | 5 | 0.588 | 0.000 |

Across all seeds and sizes the bounded success stays ≥ 0.39 and matches the
random-mapping prediction 0.538; truncations never occurred at these sizes,
consistent with the ~2^-92 production truncation probability. Toy-scale evidence,
not a production proof. Standard support for R-adding walks behaving like random
mappings: Teske (2001).

## 4. Cost accounting (dense Fp12 multiplication = 1 unit)

### Per step

* One dense Fp12 multiplication u·M_j: 1 unit (reference operation).
* Other work: partition index (~20 instr) + distinguished test (~10) + two
  255-bit exponent additions (~160) + loads/stores/control (~100) ≈ 290 instr,
  charged 300 = 0.015 units. **Per-step s = 1.015 units.**

### Total (hard bound)

The attack performs at most W_max = ceil(1.25·k0) multiplications, so

    W = W_max · s = 1.25 · 1.015 · k0  (+ one-time setup, below)
    log2(W_steps) = 127.420269 + log2(1.25·1.015) = 127.420269 + 0.342997 = 127.763266

One-time additive work, all inside or added to the budget and negligible:

| Item | Bound (units) | Ratio to main |
| --- | --- | --- |
| Pairing g = e(P,Q0) | ≤ 2^14 | ~2^-113 |
| Adding table: 2^20 entries × 2 exponentiations | ≤ 2^31 | ~2^-96 |
| Output: a mod r, [a]Q0 in G2 (< 2^22 instr), verification pairing | ≤ 2^14 | ~2^-113 |

Restart and truncation multiplications are **already inside W_max** (they are
walk multiplications), so they are not added again. Total

    log2(W) = 127.763266 + log2(1 + 2^-95) ≈ 127.763267.

Submitted bound **time_log2 = 127.80**, rounded up with margin.

## 5. Memory, preprocessing, advice

Memory is dominated by the distinguished-point table. With theta = 2^-30 and
W_max ≈ 2^127.742 multiplications, the stored records number ≈ W_max·theta
= 2^97.742 in expectation.

**Hard capacity.** Allocate a fixed table of **2^100 distinguished-point slots**
(96 bytes each). The record count is Binomial(≈ W_max, theta) with mean 2^97.742;
by a Chernoff bound it exceeds 2^100 (≈ 4.8× the mean) only with
doubly-exponentially small probability, so the fixed table is a hard storage
bound that does not bind before the collision; the section 3 success analysis is
unchanged. In the negligible event the table fills without a collision, the run
stops (counted as an exceptional unsuccessful trajectory within the ≥ 0.39
bound).

* Peak memory: 2^100 slots × 96 bytes = 2^106.585, times ~2 for indexing
  ≈ 2^107.6, plus the 2^20-entry adding table (~2^26 bytes) and 2^40 walk states
  (~2^47 bytes). Total ≤ **2^108 bytes** (reported only; no scalar contribution).
* Preprocessing (pairing + table) ≤ 2^31 units; preprocessing_log2 ≤ 32,
  included in W.
* Nonuniform advice: none (0 bytes).

The time bound is independent of theta and the memory cap: theta affects only
table size, restart count and truncation length, all negligible in W.

## 6. Heuristics, limitations, and the Frobenius speedup

* **Single heuristic:** the R-adding walk behaves like a random mapping on the r
  group elements, within the 25% budget factor (Teske 2001; toy success ~0.52).
  Not proven at production size.
* Exceptional distinguished-point-free trajectories are allowed and bounded
  (~2^-92), not assumed away; the success bound ≥ 0.39 holds over the good event.
* **Frobenius speedup (not claimed).** The order-12 Frobenius automorphism of GT
  could give a ~√12 (≈1.79-bit) speedup to ~126.3 by walking on the orbit
  quotient, but that reintroduces structural fruitless cycles whose waste I could
  not bound rigorously; the literature treats this heuristically. A costed
  TNFS/SexTNFS ledger (~2^120–2^126 in those papers' own units) is the other route
  to beat this bound.
* Analytic claim; no production-size execution.

References: Pollard (1978); van Oorschot–Wiener, J. Cryptology 1999; Teske,
Math. Comp. 2001.

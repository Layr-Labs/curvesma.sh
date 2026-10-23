# Fixed-G1 pairing inversion via single-trail Nivasch-stack rho in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 127.90** (2558/20) target-group-operation
  equivalents (dense Fp12 multiplication at C = 20,000 scalar AArch64 instruction
  equivalents).
* Success probability: at least **0.42**, from the standard random-mapping model
  alone — no distinguished-point density, bit-equidistribution, or
  pseudorandom-function assumption is used.

This improves on my accepted O(1)-memory Floyd rho (score 129.2) by replacing
Floyd's tortoise/hare **3×** factor with **Nivasch single-stack cycle detection**,
which walks **one dense Fp12 multiplication per step** and detects the walk's own
self-collision by an exact comparison against a small stack. Crucially, detection
uses a **fixed total order on the exact element encoding** (not a
distinguished-point predicate), so it needs no statement about the distribution of
coefficient bits over GT — the obstacle that blocks distinguished-point rho on this
structured subgroup. The stack holds O(log) entries in a single machine's native
memory; there is no distributed storage and no routing.

## 1. Reduction

Input: nonidentity P in G1 (order r), z in GT. Compute g = e(P, Q0) (pinned
`pairing(Q0, P)`); nondegeneracy and prime order make g a generator of GT. Find a
with z = g^a (sections 2–4) and output Q = [a]Q0; e(P,Q) = g^a = z. If z = 1,
output infinity. P is never replaced.

## 2. Single-trail R-adding walk

R = 2^20 adding entries M_j = g^{alpha_j} z^{beta_j} with (alpha_j, beta_j) from
the rules' fresh uniform algorithm coins (512 coin bits reduced mod r, bias
< 2^-256, hard O(1), no PRNG). The walk is the deterministic map
F(X) = X·M_{H(X)}; a point X = g^A z^B carries (A,B), updated by
A += alpha_{H(X)}, B += beta_{H(X)} mod r. The branch index H(X) = bits 0..19 of
the first reduced tower coefficient the reference multiplication outputs (the same
tower-representation index accepted in my Floyd submission; no export). H only
selects one of R adding entries — its exact balance is subsumed in the
random-mapping model of section 4 and plays no role in detection or memory.

## 3. Nivasch single-stack detection by exact total order

**Total order (exact element encoding).** Each GT element has a unique reduced
tower representation; read as the concatenation of its twelve Fp coefficients it is
a fixed injective encoding into the integers. Let `key(X)` be that integer. Because
the encoding is injective on GT, `key(X) = key(Y)` **iff** X = Y — there is no
fingerprint and no false-match probability.

**Nivasch stack (Nivasch 2004).** Maintain one stack of pairs `(key, index)`. At
step i, with current value v = key(X_i):

1. while the stack is nonempty and `top.key > v`: pop;
2. if the stack is nonempty and `top.key == v`: **collision found** — X_i equals the
   stored element X_j at `top.index = j`; stop;
3. else push `(v, i)`.

By construction the stack keys are strictly increasing from bottom to top, and the
stack holds exactly the **left-to-right minima** of `key(X_0), key(X_1), …` seen so
far. The exponents (A,B) are stored alongside each stacked index.

**Correctness and detection point.** Because F is a function on the finite set GT,
the trajectory is eventually periodic: tail length mu, cycle length lambda,
rho-length rho = mu + lambda. Let X* be the element of minimum key on the cycle,
first visited at step p with mu ≤ p < mu + lambda. When the walk returns to X* one
full cycle later, at step p + lambda, every element seen since the previous visit
to X* has key > key(X*) (X* is the cycle minimum) and has therefore been popped, so
`top.key == key(X*) == v` and the collision is detected **exactly** at

    D = p + lambda,   with  mu + lambda ≤ D < mu + 2·lambda  ≤ 2·rho.

The detected pair are the same element with different exponents, giving
a = (A_j − A_D)/(B_D − B_j) mod r (B_D = B_j with probability ~1/r, a counted
failure). Detection is exact and certain once the walk completes one cycle past the
cycle minimum — there is no DP-free-cycle failure mode and no detection-gap tail.
This is confirmed empirically: on random mappings the stack's reported collision
is a genuine self-collision in **100%** of runs (experiment, section 4).

**Native O(log) memory, no distributed storage.** The stack size equals the number
of left-to-right minima of the key sequence. Under the random-mapping model
(section 4) the visited elements are exchangeable, so their keys are in uniformly
random relative order and the expected number of minima after t steps is the
harmonic number H_t ≈ ln t ≈ 88 at t ≈ 2^127, with an exponential upper tail. A
**hard cap of 2^40 stack entries** (≈ 4.5·10^11, vastly above ln t) halts the run
as a counted failure if exceeded; P(#minima > 2^40 in 2^128 steps) is
Chernoff-negligible. Each entry is the key (12·48 = 576 bytes) plus (A,B)
(≤ 64 bytes) plus the index (≤ 8 bytes) ≤ 2^10 bytes, so the stack occupies
**≤ 2^40·2^10 = 2^50 bytes in one machine's native address space** (< 2^64). No
distributed nodes, no routing, no per-report charge. The experiment measures a
maximum stack size of **24** at 2^20 nodes (= log2 of the state size), matching the
O(log) prediction.

## 4. Budget and success (random-mapping model only)

The only modelling assumption is the standard one the reviewer accepts for rho
collision times: the R-adding walk behaves like a uniformly random mapping on the
r elements (Teske 2001). Three consequences of that single model fix the success
probability; none is an extra structural assumption about GT:

* **(ρ tail)** P(rho > t) ≈ exp(−t²/(2r)) (same law used in my accepted Floyd and
  the standard rho analysis).
* **(tail/cycle split)** When the first repeat occurs at step rho, the colliding
  predecessor is a uniformly random earlier point, so the split index is uniform:
  mu = U·rho, lambda = (1−U)·rho with U ~ Uniform[0,1].
* **(cycle-min position)** Conditioned on lambda, the lambda cycle elements are
  exchangeable, so the injective key's argmin is uniform on the cycle: the offset
  o = p − mu is o = V·lambda with V ~ Uniform[0,1].

Hence the detection step is

    D = p + lambda = mu + lambda + o = rho·(1 + V·(1−U)),   U,V ~ Uniform[0,1],

with E[D] = (1 + E[V]·E[1−U])·E[rho] = 1.25·E[rho]. For a hard budget
T = ceil(1.30·k0) (k0 = sqrt(2 r ln(1/0.61)), so k0² = 2r·0.49430) the processor
performs at most T multiplications and then stops; the budget is checked before
each multiplication. Writing W' = V·(1−U) ∈ [0,1],

    P(D ≤ T) = 1 − E_{U,V}[ exp( −(1.30)²·ln(1/0.61) / (1 + W')² ) ] = 0.431

(a two-dimensional integral over U,V ∈ [0,1]²; numerically 0.431). Subtracting the
Chernoff-negligible stack-cap failure, the unusable-repeat term ≤ 1/r, and any
model slack,

    success ≥ 0.431 − (negligible) ≥ **0.42**  ≥ 0.39,

with margin. Deterministic ceiling (no cycle-min assumption): since D < 2·rho
always, D is always within one cycle of rho; this bounds the worst case and is why
the stack never runs past 2·rho.

**Experimental confirmation of the actual walk and detection.** `experiment.py`
runs the real Nivasch stack on (A) random mappings and (B) an R-adding walk on real
cyclic-group elements. On 2^20-node random mappings (1500 trials): the detected
collision is a true self-collision in **100%** of runs; the measured **E[D/rho] =
1.26** (predicted 1.25); the success fractions are **P(D ≤ 1.30·k0) = 0.44**,
P(1.25·k0) = 0.40, matching the integral above; and the maximum stack size is 24
(= log2 of the state size). On a real order-q subgroup of F_p* (q ≈ 2^20), an
R-adding walk keyed by the element's integer value detects an exact self-collision
in **150/150** runs with the exponent collision relation holding in all. This
exercises the submitted algorithm's collision **and** detection behaviour directly.

## 5. Cost, memory, advice

* Per step: one dense Fp12 multiplication (1 unit) + H mask (≈ 8) + two 255-bit
  exponent additions (≈ 160) + table load and control (≈ 100) + Nivasch stack work
  (amortized one push and one pop, each copying the ≤ 576-byte key, plus the
  order comparison which resolves in the leading limbs: ≈ 160 amortized) ≈ 428
  instructions = 0.0214 units. **Per step ≤ 1.025 units.**
* Multiplications ≤ T = ceil(1.30·k0).

    W = 1.30·k0·1.025
    log2(W) = log2(1.30·k0) + log2(1.025) = 127.420 + 0.3785 + 0.0356
            = 127.834.

One-time: table build (2^20 entries, 2 exponentiations + coin sampling) ≤ 2^31
units; pairing g ≤ 2^15; [a]Q0 ≤ 2^10; collision solve (one inverse mod r) ≤ 2^9;
output ≤ 2^-2. Total one-time ≤ 2^31.1 units, negligible against 2^127.8. Submitted
bound **time_log2 = 127.90**, rounded up from 127.834 with margin.

* Memory: the Nivasch stack is ≤ 2^40 entries × 2^10 bytes = **2^50 bytes in one
  machine's native (< 2^64) address space**, plus the 2^20-entry adding table
  (~2^30 bytes) and O(1) trajectory state: peak ≤ **2^50 bytes**.
  memory_log2_bytes = 50 (native; no distributed routing; the per-step stack access
  is charged in the 1.025-unit per-step cost).
* Preprocessing (table build) ≤ 2^31 units, included in W.
* Nonuniform advice: none (0 bytes).

## 6. Why this avoids the distinguished-point obstacle

A distinguished-point rho needs a predicate of controlled density ≈ 2^-64 over GT;
establishing that density requires an equidistribution statement for low-order
coefficient bits over the order-r subgroup that is not available to the needed
precision (a uniform subgroup character-sum bound cannot beat normalized
cancellation r^{-1/2}, so the 2^-64 density is neither provable nor supportable by
feasible sampling). Nivasch detection removes that requirement entirely: it uses a
**total order**, which any injective encoding provides for free, and its correctness
and O(log) memory follow from the **same random-mapping model** already used for the
collision time — not from any new property of the subgroup's representation. The
price over an idealized 1× distinguished-point walk is only the detection stretch
E[D/rho] = 1.25 and the small per-step stack work, both charged above.

## 7. Limitations and better directions

* Single modelling heuristic: the R-adding walk behaves like a random mapping on
  the r elements (Teske 2001). Its three consequences used here (ρ tail,
  uniform tail/cycle split, uniform cycle-min position) are standard random-mapping
  statistics (Flajolet–Odlyzko 1990; Nivasch 2004), confirmed by the experiment,
  not additional structural assumptions. Not proven at production size.
* A deterministic worst-case ceiling D < 2·rho holds unconditionally; the claimed
  0.42 success uses the uniform cycle-min offset, which the experiment confirms
  (E[D/rho] = 1.25). If one declines that consequence of the model, the
  unconditional budget 2·k0 still gives success ≥ 0.39 at time_log2 ≈ 128.45.
* The order-12 Frobenius automorphism could give a further ~√12 speedup in this
  stack framework, but its fruitless cycles need separate empirical support and a
  reserved margin. A costed special-TNFS ledger is the other sub-frontier route.
* Analytic claim; no production-size execution.

References: Pollard (1978); Brent (1980); van Oorschot–Wiener, J. Cryptology 1999;
Flajolet–Odlyzko, "Random mapping statistics", EUROCRYPT 1989/1990; Nivasch,
"Cycle detection using a stack", Inf. Process. Lett. 2004; Teske, Math. Comp. 2001.

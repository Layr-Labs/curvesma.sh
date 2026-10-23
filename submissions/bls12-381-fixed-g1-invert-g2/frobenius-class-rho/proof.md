# Frobenius-class parallel Pollard rho in GT (fixed G1 argument, output in G2)

Claim: time_log2 = **125.89** (total work W <= 2^125.89 dense-Fp12-multiplication
equivalents, C = 20,000), success probability >= 0.40, peak memory < 2^74 bytes,
preprocessing <= 2^31 units (included), no advice. Heuristics are listed in section 6.

Design provenance: the GT walk below (key formula, twisted table, cycle handling)
follows the Frobenius-class rho of public submission
3da91365-2d58-40e5-ac7d-415e3061f7d8 on the mirror track (fixed G2 argument,
output in G1). This package redoes the reduction for this track, adds a hard
global step cap, a hard distinguished-record cap and a deterministic restart
bound, and redoes every number. Nothing here relies on that package's review.

## 1. Reduction

Input: nonidentity P in G1 (uniform, fixed for the instance) and z in GT. Compute
g = e(P, Q0) with the pinned convention (library call pairing(Q0, P)). P has prime
order r and the reduced pairing is nondegenerate, so g != 1 and g generates GT
(prime order r). Every valid z equals g^x for a unique x in [0, r). If z = 1,
output Q = infinity (e(P, infinity) = 1). Otherwise recover x by the rho algorithm
below, output Q = [x]Q0 and check e(P, Q) = z. By bilinearity
e(P, [x]Q0) = e(P, Q0)^x = g^x = z. The answer is for the supplied P; no other
P is used, and Q0 is the public profile generator.

## 2. Arithmetic facts used (checked numerically in `experiment.py`, part A)

Represent Fp12 as blst does: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi),
xi = 1+u, basis w^0..w^5 (blst's v^i w^j with m = 2i+j). This is the profile
field: w^6 = 1+u gives u = w^6 - 1 and w^12 - 2w^6 + 2 = 0, so conversion to the
profile's coefficient list is c_m = a_m - b_m, c_{m+6} = b_m (free).

* Frobenius: (sum c_m w^m)^p = sum conj(c_m) d_m w^m with d_m = xi^{m(p-1)/6}.
* p has multiplicative order exactly 12 modulo r (embedding degree 12), so the
  automorphism group {y -> y^{p^k}} of GT has order 12; it contains inversion
  (p^6 = -1 mod r). Every y != 1 lies in a class [y] of exactly 12 elements.
  Number of classes N = (r-1)/12, log2 N = 251.2721.
* Let zeta = xi^{(p^2-1)/6}. It lies in Fp and is a primitive 6th root of unity
  (zeta^2 = zeta - 1). Writing c_1 = a + b'u for the w^1 coefficient of y, the
  Re(w^1 coefficient) of the conjugate y^{p^k} is
  `key_k = a*zeta^j` for k = 2j and `key_k = B*zeta^j` for k = 2j+1, where
  `B = Re(conj(c_1) d_1) = a*Re(d_1) + b'*Im(d_1)`.
  Since zeta^j in {1, zeta, zeta-1, -1, -zeta, 1-zeta}, all 12 keys come from
  **4 Fp multiplications** (a*zeta, B as two products, B*zeta) plus 8
  subtractions/negations. No Frobenius of the whole element is ever computed.

Define kappa(y) = argmin_k key_k (keys as canonical residues in [0,p), compared as
integers) and the class key K([y]) = min_k key_k. Then y^{p^{kappa(y)}} is the same
element for all 12 members of a class (canonical representative y^), and K is a
class invariant. Ties need a 1-in-p linear coincidence (probability <= 66/p per
step; over at most L_max = 2^125.70 steps the total is <= 2^-248.9); a tie is
counted as failure.

## 3. Walk (R-adding walk on classes, Duursma-Gaudry-Morain / Wiener-Zuccherato style)

Parameters: R = 1024 partitions, distinguished-point fraction theta = 2^-60,
M = 2^20 parallel walks, each started at g^alpha z^beta with fresh uniform alpha,
beta mod r. All walks share one table and one class map.

Precompute M_h = g^{a_h} z^{b_h} (h < R, random a_h, b_h) and the table
T[h][k] = M_h^{p^k}, k = 0..11 (12,288 entries of 576 bytes, 7.1 MB), with exponent
tables a_h p^k, b_h p^k mod r.

Class map: F([y]) = [ y^ * M_h ], h = K([y]) mod R. A walk keeps an arbitrary
member y of the current class with its exponents (alpha, beta), y = g^alpha z^beta.
Because y^{p^kappa} * M_h = (y * M_h^{p^{-kappa}})^{p^kappa}, one step is

```
compute keys, kappa, K                      (4 Fp mul + <= 700 other instructions)
if bits 10..69 of K are all zero (distinguished): report (K, p^kappa*alpha, p^kappa*beta)
h = K mod R (bits 0..9) ; t = (12 - kappa) mod 12
y = y * T[h][t]                             (ONE dense Fp12 multiplication)
alpha += a_h p^t ; beta += b_h p^t (mod r)  (2 additions mod r)
```

So each step costs exactly one reference operation plus the key work. The walk is
a deterministic function of the class, so once two walks reach the same class they
follow the same class sequence and reach the same distinguished class. A walk
continues after reporting; it is never restarted after a distinguished point.

Distinguished store. A report is (bits 70..197 of K, canonical exponents A, B):
80 bytes, kept in a hash table at load <= 2/3 keyed by that 128-bit field. Bits
0..69 are not stored: on a distinguished class bits 10..69 are zero and bits 0..9
are the partition index, so they carry little information. On a key
match: if the stored (A, B) equals the new one, the report is ignored (same walk
history, nothing inserted); otherwise compute x = (A1 - A2)/(B2 - B1) mod r and
check g^x = z. On success output; on failure (a 128-bit false match) keep the
stored record. With equal K and equal class, p^k1(alpha1 + beta1 x) =
p^k2(alpha2 + beta2 x) mod r, and the stored exponents are the canonical ones, so
the formula is exact. B1 = B2 has probability 1/r per collision.
**Hard report cap:** a global counter counts every report (inserted or ignored);
the algorithm processes at most 2^67 reports and stops with failure at the next
one. Hence at most 2^67 records are ever stored.

**Hard global step cap:** every step of every kind (ordinary steps, squaring
escapes, cycle traversals, steps of walks that are later abandoned) increments one
global counter. At L_max = 1.003 S + E = 2^125.6992 steps (S, E in section 4) the
algorithm stops with failure.

**Fruitless cycles.** A twisted adding walk can return to a class with formally
cancelled exponents. Exact 2-cycles need h' = h and kappa' = 6 (inverse conjugate):
rate 1/(12R) = 2^-13.58 per step (toy measurement over 2,400 trials: 1.303e-3 vs
predicted 1.302e-3 at R = 64). Each step compares K with the key two steps earlier; on equality the
two classes {C1, C2} form a 2-cycle; the walk squares the stored representative of
the one with smaller K (y -> y^2 maps a class to a class: (y^{p^k})^2 = (y^2)^{p^k}),
doubling its exponents. This escape depends only on the cycle, so merged walks stay
merged. Cost <= 3 steps per event. Longer fruitless cycles need >= 2 repeated
partition indices: rate <= R^-2 per step (3-cycles need three equal indices and
twists {1, p^4, p^8}, rate about 1/(72 R^2)). A checkpoint key is stored every 64
steps and compared at every step. A cycle of length <= 64 entered at any time
contains the next checkpoint within 64 steps and is detected at most 64 steps
later; it is then traversed once more to find its minimum-K class, which is
squared. Waste <= 200 steps per event. Cycles longer than 64 have rate < R^-32 and
are handled by the stuck-walk rule.

**Stuck-walk rule:** a walk that has taken 40/theta = 2^65.32 steps since its last
distinguished point (or its start) is abandoned and restarted at a fresh
g^alpha z^beta. Since every abandoned walk consumed 2^65.32 counted steps,
restarts number at most L_max / 2^65.32 = 2^60.38, deterministically.
Squaring steps cannot be part of a fruitless cycle (2^t is not in <p> mod r for
small t).

## 4. Success probability

Heuristic (H1): the class sequence of the walk behaves like a random mapping on
N classes, with the Brent-Pollard correction for R-adding walks (effective
N' = N * R/(R-1)). Then the probability that no useful class collision occurs among
the first S distinct-step points of all walks is exp(-S^2/(2N')). Choose

S0 = sqrt(2 N' ln(1/0.6)) = 2^125.6522,   S = 1.03 * S0 = 2^125.6949,

where the 3% is explicit slack for walk non-randomness; at S0 the model gives
exactly 0.40, at S it gives 0.418.

Detection: after the first useful collision, the merged trajectory reaches a
distinguished class within 40/theta further steps except with probability e^-40.
Charging all M walks in lockstep over that window costs E = M*40/theta = 2^85.32
steps.

Step cap: useful (distinct-point) steps are at least L_max - E - waste. Expected
waste is at most (3/(12R) + 200/R^2) = 0.000435 of all steps; it is a sum of about
2^112 small independent-looking events, so under H1 it stays below 0.003 S except
with negligible probability, leaving at least S useful steps before the cap.

Report cap: under H1 each step lands on a distinguished class with probability
theta, so reports among at most L_max steps have mean at most 2^65.70; reaching
2^67 (2.46 times the mean) has Chernoff probability below exp(-2^65), negligible.

Abandoned walks: a 40/theta stretch with no distinguished point has probability
e^-40 = 2^-57.7. Every stretch begins at a report, a start or a restart, so there
are at most 2^67 + 2^20 + 2^60.38 < 2^67.02 stretches; the expected number of
abandonments is below 2^9.4 and the steps they strand (below 2^74.7) are
negligible against 0.003 S. A true record lost to a 128-bit false-match slot has probability
at most 2^67/2^128 per collision.

Self-collisions of a single walk (fraction about 1/M = 2^-20 of collisions),
B1 = B2 (1/r), ties (2^-248.9), the e^-40 detection tail and the caps above are
absorbed by the 0.418 vs 0.40 margin. Success >= 0.40 (the track requires 0.39).

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Per-step non-reference instructions, scalar AArch64, using blst's measured
blst_fp_mul = 383 instructions:

| Work | Instructions |
| --- | ---: |
| 4 Fp multiplications (a*zeta, two products for B, B*zeta) | 1,532 |
| 1 Fp addition (B) + 8 Fp subtractions/negations, <= 40 each | 360 |
| min of 12 six-limb keys (top-limb compare, csel; rare full compare) | 300 |
| hash, distinguished test, t = (12-kappa) mod 12, table address | 40 |
| alpha, beta += table values mod r (2 x 4-limb add, compare, select) | 100 |
| 2-cycle and checkpoint key comparisons, key history stores, buffer rotation, step counter | 60 |
| loop control | 30 |
| **total (budgeted 2,800)** | **2,422** |

The product y*T[h][t] is the dense blst_fp12_mul reference operation (operands
read through pointers, output written to a rotating buffer; no copies), charged 1.
A squaring escape uses one Fp12 squaring (at most one dense multiplication) and
doubles two exponents, inside the same budget. Every counted step therefore costs
at most c = 1 + 2,800/20,000 = 1.14 units, and the walk costs at most c * L_max.

Additive terms (all deterministic bounds):

* Setup and output, <= 2^31 together: pairing g = e(P, Q0) (Miller loop: 64
  iterations of one Fp12 squaring, at most two line multiplications counted as
  dense and one G2 doubling or addition at <= 1.44 units, so <= 64 x 4.44 < 2^9;
  final exponentiation even by plain square-and-multiply on the 4,314-bit exponent
  (p^12-1)/r: <= 8,628 units; together < 2^14); table: 1024 multipliers x
  (2 x 255-bit exponentiations <= 1,020 units + 1 product + 11 Frobenius maps
  <= 1 unit each) <= 2^20.02; M = 2^20 walk starts g^alpha z^beta at <= 1,021 units
  each, including drawing alpha and beta, <= 2^30.0; output Q = [x]Q0 by double-and-add, <= 510 Jacobian G2 operations,
  each <= 20 Fp2 multiplications (<= 3 Fp multiplications + 5 Fp additions each)
  plus <= 20 Fp2 additions and control, <= 28,780 instructions = 1.44 units, so
  <= 735 units, plus one Fp2 inversion (Fermat, <= 762 Fp multiplications,
  <= 15 units); verification pairing e(P, Q) <= 2^14.
* Restarts: <= 2^60.38 x 1,021 units <= 2^70.38.
* Distinguished-point handling: <= 2^67 + 1 reports (report cap). Each report
  costs <= 1,000 instructions (2 modular multiplications by p^kappa, hash probe,
  insert or ignore): <= 2^62.68 units.
* Candidate solves: at most one per report, <= 520 units each (inversion mod r,
  exponentiation g^x <= 510, comparison): <= (2^67 + 1) x 520 <= 2^76.03.

Additive total <= 2^31 + 2^62.68 + 2^70.38 + 2^76.03 <= 2^76.07.

W <= c * L_max + 2^76.07
  = 1.14 * (1.003 * 2^125.6949 + 2^85.32) + 2^76.07
  = 2^125.8882 < 2^125.89.

Submitted bound: **time_log2 = 125.89** (rounded up).

Memory: 2^67 records x 80 bytes at load 2/3 = 2^73.91 bytes, plus the 7.1 MB table
and 2^20 walk states (current element, two history elements, keys and exponents,
at most 2 KB each, 2^31 bytes): < 2^74 bytes. This is a hard bound by the
record cap. Memory is reported, not scored.

## 6. Heuristics and limitations

* H1 (random-mapping behaviour of an R-adding walk on Frobenius classes) is the
  standard heuristic behind every rho estimate (Teske 1998/2001; Brent-Pollard;
  van Oorschot-Wiener 1999 for parallel distinguished points; Duursma-Gaudry-Morain
  1999 and Bos-Kleinjung-Lenstra 2010 for automorphism classes and fruitless
  cycles). It is supported here by the toy runs; it is not proved.
* Distinguished-point and partition bits (bits 0..69) and the stored record
  field (bits 70..197) are taken from K; K is the minimum of 12 values, which is
  not uniform in [0,p), but these bits are treated as pseudo-random functions of
  the class (heuristic, standard). The 2^67/2^128 false-match bound of section 4
  rests on this. The toy does not test it (see section 7).
* The waste, record-count and abandonment estimates in section 4 use H1; the caps
  themselves, the restart count and every cost term in section 5 are deterministic.
* No production run. Toy experiments are evidence only at their own scale.
* This is a generic attack (about the square root of the class count). It
  ignores index-calculus methods in Fp12 (TNFS family); it is a valid upper
  bound, not a security estimate.

## 7. Experiments (`experiment.py`, standard library only)

Part A (BLS12-381 itself): for random GT elements (random Fp12 element raised to
(p^12-1)/r), checks y^r = 1; Frobenius formula = y^p; the 4-multiplication key
formula equals the keys of all 12 explicitly computed conjugates and they are
distinct; kappa selects the identical canonical element from each of the 12
conjugates; p^6 is inversion; ord_r(p) = 12; 1 + p^4 + p^8 = 0 mod r. The pairing
reduction of section 1 is pure bilinearity and is not part of the toy.

Part B (toy, same tower shape): q = 163, r' = 54,298,861 dividing q^4 - q^2 + 1,
N = (r'-1)/12. Runs the full walk (key formula, kappa trick, twisted table,
2-cycle escape, checkpoint escape, distinguished points, stuck-walk restart after
40/theta steps as in production, collision solving and verification) over many
random instances. Reported: steps to first useful class collision / sqrt(N) (random
mapping predicts mean 1.2533), fraction of trials whose collision time is within
the 0.40 budget, measured 2-cycle rate vs 1/(12R). The toy does not apply the two
production caps; its counter counts walk steps, not setup work. The toy takes
its partition index and distinguished test from a hash of the full canonical
element, not from bits of K: at q = 163 the key K is a residue below 163 (about
7 bits), fewer than the 10 partition bits needed. So the toy tests the class map,
the kappa selection, the twisted table and the escapes, but not the
pseudo-randomness of production K's bits, which remains the heuristic of
section 6.

The initial input (q = 163, R = 64, 4 walks, seed 1, 80 trials, 14 s limit) is
deterministic given its seed and took about 6 s on our machine: 80 trials, all logs
recovered and verified; mean first-collision time 1.2117 sqrt(N) (model 1.2533);
fraction within the 0.40 budget 0.400 (model 0.40); 2-cycles per step 1.318e-3
(predicted 1.302e-3); 278 two-cycle escapes, 1 longer cycle, 10 ignored repeat
reports, no abandoned walks. A slower runner that hits the 14 s limit completes
fewer trials and reports different sample values.

Larger offline runs of this same file (time limit lifted, so these inputs are
outside the runner's 14 s schema; production_checks 0; 3,096 instances, all logs
recovered and verified; standard errors from the random-mapping spread 0.655/sqrt(n)
and sqrt(0.24/n)):

| q, r', R, dp_bits, walks, seed | trials | mean T_c/sqrt(N) (model incl. R/(R-1)) | P[T_c <= S0] (model 0.40) | 2-cycles/step (1/(12R)) |
| --- | ---: | ---: | ---: | ---: |
| 163, 54298861, 64, 4, 4, 1 | 2400 | 1.2622 +- 0.013 (1.2632) | 0.400 +- 0.010 | 1.303e-3 (1.302e-3) |
| 787, 64571173, 256, 4, 4, 2 | 600 | 1.1951 +- 0.027 (1.2558) | 0.427 +- 0.020 | 3.29e-4 (3.26e-4) |
| 331, 923346397, 1024, 6, 8, 3 | 96 | 1.1157 +- 0.067 (1.2539) | 0.490 +- 0.050 | 9.3e-5 (8.1e-5) |

The two smaller samples lie about two standard errors below the model mean, in the
direction favourable to the claim; the large sample matches it. Longer fruitless
cycles were rarer than the R^-2 bound (56 events in about 6.6e6 steps at R = 64,
i.e. 8.4e-6 per step vs bound 2.4e-4). No walk was abandoned in any run.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999);
Teske, Math. Comp. 70 (2001); Duursma, Gaudry & Morain, ASIACRYPT 1999; Wiener &
Zuccherato, SAC 1998; Bos, Kleinjung & Lenstra, ANTS 2010.

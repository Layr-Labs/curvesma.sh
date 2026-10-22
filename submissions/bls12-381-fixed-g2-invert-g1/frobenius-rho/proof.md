# Frobenius-class parallel Pollard rho in GT

Claim: time_log2 = **125.89** (total work W <= 2^125.89 dense-Fp12-multiplication
equivalents, C = 20,000), success probability >= 0.40, peak memory <= 2^73 bytes,
preprocessing <= 2^31 units (included), no advice. Heuristics are listed in section 6.

## 1. Reduction

Input: nonidentity Q in G2 and z in GT. Compute g = e(P0, Q) with the pinned
convention. Q has prime order r and the reduced Tate/ate pairing is nondegenerate,
so g != 1 and g generates GT (prime order r). Every valid z equals g^x for a unique
x in [0, r). Recover x by the rho algorithm below, output P = [x]P0, check
e(P, Q) = z. Then e(P, Q) = g^x = z. If z = 1, output P = infinity directly.
The answer is for the supplied Q; no other Q is used.

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
step, total <= 2^-249); they are counted as failure.

## 3. Walk (R-adding walk on classes, Duursma-Gaudry-Morain / Wiener-Zuccherato style)

Parameters: R = 1024 partitions, distinguished-point fraction theta = 2^-60,
M = 2^20 parallel walks, each started at g^alpha z^beta with fresh uniform alpha, beta mod r.

Precompute M_h = g^{a_h} z^{b_h} (h < R, random a_h, b_h) and the table
T[h][k] = M_h^{p^k}, k = 0..11 (12,288 entries of 576 bytes, 7.1 MB), with exponent
tables a_h p^k, b_h p^k mod r.

Class map: F([y]) = [ y^ * M_h ], h = K([y]) mod R. A walk keeps an arbitrary
member y of the current class with its exponents (alpha, beta), y = g^alpha z^beta.
Because y^{p^kappa} * M_h = (y * M_h^{p^{-kappa}})^{p^kappa}, one step is

```
compute keys, kappa, K                      (4 Fp mul + ~650 other instructions)
if bits 10..69 of K are all zero (distinguished): report (K, p^kappa*alpha, p^kappa*beta)
h = K mod R (bits 0..9) ; t = (12 - kappa) mod 12
y = y * T[h][t]                             (ONE dense Fp12 multiplication)
alpha += a_h p^t ; beta += b_h p^t (mod r)  (2 additions mod r)
```

So each step costs exactly one reference operation plus the key work. The walk is
a deterministic function of the class, so once two walks reach the same class they
follow the same class sequence and reach the same distinguished class.

Collision: two reports with equal K (and equal canonical element, confirmed by
recomputing it from the stored exponents) give p^k1(alpha1 + beta1 x) =
p^k2(alpha2 + beta2 x) mod r, i.e. x = (A1 - A2)/(B2 - B1) with the stored canonical
exponents. B1 = B2 has probability 1/r per collision for random starting exponents.
Verify g^x = z (<= 510 units) before output.

**Fruitless cycles.** A twisted adding walk can return to a class with formally
cancelled exponents. Exact 2-cycles need h' = h and kappa' = 6 (inverse conjugate):
rate 1/(12R) = 2^-13.6 per step (toy measurement: 1.316e-3 vs predicted 1.302e-3 at
R = 64). Each step compares K with the key two steps earlier; on equality the two
classes {C1, C2} form a 2-cycle; the walk squares the stored representative of the
one with smaller K (y -> y^2 maps a class to a class: (y^{p^k})^2 = (y^2)^{p^k}),
doubling its exponents. This escape depends only on the cycle, so merged walks stay
merged. Cost <= 3 units per event. Longer fruitless cycles need >= 2 repeated
partition indices: rate <= R^-2 per step (3-cycles need three equal indices and
twists {1, p^4, p^8}, rate about 1/(72 R^2)). A checkpoint key is stored every 64
steps and compared at every step; a repeat within 64 steps means a cycle of length
<= 64, which is traversed once more to find its minimum-K class, which is squared.
Waste <= 200 steps per event. Cycles longer than 64 have rate < R^-32 and are
handled by the stuck-walk rule (no distinguished point within 20/theta steps:
restart). Squaring steps cannot be part of a fruitless cycle (2^t is not in
<p> mod r for small t).

## 4. Success probability

Heuristic (H1): the class sequence of the walk behaves like a random mapping on
N classes, with the Brent-Pollard correction for R-adding walks (effective
N' = N * R/(R-1)). Then the probability that no useful class collision occurs among
the first S distinct-step points of all walks is exp(-S^2/(2N')). Choose

S0 = sqrt(2 N' ln(1/0.6)) = 2^125.6522,   S = 1.03 * S0 = 2^125.6949,

where the 3% is explicit slack for walk non-randomness; at S0 the model gives
exactly 0.40, at S it gives 0.418. Toy evidence (section 7) is consistent with
the random-mapping mean (measured 1.231 sqrt(N) vs 1.2533 sqrt(N) predicted).

Detection: after the first useful collision, the merged trajectory reaches a
distinguished class within 20/theta further steps except with probability e^-20.
Charging all M walks in lockstep over that window costs E = M*20/theta = 2^84.3
steps. Self-collisions of a single walk (fraction ~1/M = 2^-20 of collisions) and
B1 = B2 (1/r) are absorbed by the 0.418 vs 0.40 margin. Success >= 0.40.

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
| 2-cycle and checkpoint key comparisons, key history stores, buffer rotation | 60 |
| loop control | 30 |
| **total (budgeted 2,800)** | **2,422** |

The product y*T[h][t] is the dense blst_fp12_mul reference operation (operands
read through pointers, output written to a rotating buffer; no copies), charged 1.
Step cost c = 1 + 2,800/20,000 = 1.14 units.

Fruitless-cycle overhead factor: 1 + 3/(12R) + 200/R^2 = 1.000435.

Setup and output, all <= 2^31 units together: pairing g (Miller loop + final
exponentiation, < 2^12 units, bounded by 2^14), input check z^r = 1 (<= 510),
table: 1024 multipliers x (2 x 255-bit exponentiations <= 1,021 units + 11
Frobenius maps <= 1 unit each) <= 2^20.02, M = 2^20 random walk starts
g^alpha z^beta <= 2^30.0, distinguished-point handling (2^65.7 reports x <= 1,000
instructions: hash-table insert, 2 modular multiplications by p^kappa) <= 2^61.3,
output [x]P0 (<= 510 G1 operations, each < 1 unit) and verification pairing
(<= 2^14). Total bounded by 2^31 + 2^62.

W = c * S * 1.000435 + c * E + 2^62
  = 1.14 * 2^125.6949 * 1.000435 + 1.14 * 2^84.3 + 2^62
  = 2^125.8846.

Submitted bound: **time_log2 = 125.89**.

Memory: 2^65.7 distinguished reports x 80 bytes (truncated 128-bit key + two
256-bit exponents), hash-table load 2/3, plus 7.1 MB table and 2^20 walk states
(~1 KB each): < 2^73 bytes. Memory is reported, not scored.

## 6. Heuristics and limitations

* H1 (random-mapping behaviour of an R-adding walk on Frobenius classes) is the
  standard heuristic behind every rho estimate (Teske 1998/2001; Brent-Pollard;
  van Oorschot-Wiener 1999 for parallel distinguished points; Duursma-Gaudry-Morain
  1999 and Bos-Kleinjung-Lenstra 2010 for automorphism classes and fruitless
  cycles). It is supported here by the toy runs; it is not proved.
* Distinguished-point and partition bits are taken from K; K is the minimum of
  12 values, which is not uniform in [0,p), but the low 60 bits and the hash
  index used are pseudo-random functions of the class (heuristic, standard).
* No production run. Toy experiments are evidence only at their own scale.
* This is a generic attack (about the square root of the class count). It
  ignores index-calculus methods in Fp12 (TNFS family); it is a valid upper
  bound, not a security estimate.

## 7. Experiments (`experiment.py`, standard library only)

Part A (BLS12-381 itself): for random GT elements (random Fp12 element raised to
(p^12-1)/r), checks y^r = 1; Frobenius formula = y^p; the 4-multiplication key
formula equals the keys of all 12 explicitly computed conjugates and they are
distinct; kappa selects the identical canonical element from each of the 12
conjugates; p^6 is inversion; ord_r(p) = 12; 1 + p^4 + p^8 = 0 mod r.

Part B (toy, same tower shape): q = 163, r' = 54,298,861 dividing q^4 - q^2 + 1,
N = (r'-1)/12. Runs the full walk (key formula, kappa trick, twisted table,
2-cycle escape, checkpoint escape, distinguished points, collision solving and
verification) over many random instances. Reported: steps to first useful class
collision / sqrt(N) (random mapping predicts mean 1.2533), fraction of trials whose
collision time is within the 0.40-budget, measured 2-cycle rate vs 1/(12R).
Offline runs of the same file (seeds and inputs as shown; 3,096 instances, all
logs recovered and verified):

| q, r', R | trials | mean T_c/sqrt(N) (predicted) | P[T_c <= S0] (model 0.40) | P[T_c <= 1.03 S0] | 2-cycles/step (1/(12R)) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 163, 54298861, 64 | 2400 | 1.2605 +- 0.013 (1.2632) | 0.394 +- 0.010 | 0.413 | 1.30e-3 (1.30e-3) |
| 787, 64571173, 256 | 600 | 1.2485 +- 0.027 (1.2558) | 0.410 +- 0.020 | 0.428 | matches |
| 331, 923346397, 1024 | 96 | 1.274 +- 0.066 (1.2539) | 0.375 +- 0.050 | 0.385 | matches |

Long fruitless cycles were rarer than the R^-2 bound (58 events in ~8.4e6 steps at
R = 64, i.e. 7e-6 vs bound 2.4e-4).

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999);
Teske, Math. Comp. 70 (2001); Duursma, Gaudry & Morain, ASIACRYPT 1999; Wiener &
Zuccherato, SAC 1998; Bos, Kleinjung & Lenstra, ANTS 2010.

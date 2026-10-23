# Frobenius-class parallel Pollard rho in GT (success target 0.39)

Claim: time_log2 = **125.861** (total work W <= 2^125.861 dense-Fp12-multiplication
equivalents, C = 20,000), success probability >= 0.39, peak memory <= 2^73 bytes,
preprocessing <= 2^31 units (included), no advice. Heuristics are listed in section 6.

Design provenance: the algorithm, the four-multiplication key formula, the twisted
table and the cost ledger are those of public submission
3da91365-2d58-40e5-ac7d-415e3061f7d8 (Frobenius-class rho, 125.89). This package keeps
that algorithm unchanged and differs only in: (1) the success target is the track's
0.39 instead of 0.40 (same 3% step slack, same margin structure); (2) a written
derivation of the fruitless-cycle rates (section 3.2); (3) an explicit stuck-walk rule
and record format; (4) toy experiments that run the production partition and
distinguished-point predicates on the raw class key. No earlier result is used as
evidence; every bound below is derived or computed here.

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
  Note that inversion (k = 6) maps key a to -a, so the 12 keys are pairwise
  distinct for a generic element (checked in part A); there is no built-in tie.

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

Record format: a report is stored as 80 bytes = bits 70..197 of K (16 bytes; the 128
bits directly above the distinguished bits, since bits 10..69 of every reported K are
zero by construction) plus the canonical exponents A = p^kappa alpha, B = p^kappa beta
mod r (32 bytes each); the hash table is keyed on those 128 bits. A match on the
truncated key triggers a confirmation: the canonical element
g^A z^B is recomputed for both records (<= 1,021 units each) and compared in full. A
truncated-key match between different classes has probability <= 2^-128 per pair
(heuristic H1 for the low key bits); there are at most (2^65.67)^2/2 = 2^130.34 pairs,
so the expected number of false matches is <= 2^2.34 and, by Markov's inequality,
more than 2^40 occur with probability < 2^-37. Each false match only costs its
confirmation (<= 2,042 units) and is discarded, so 2^40 of them cost < 2^51 units,
inside the 2^62 term of section 5.

Collision: two reports with equal K (and equal canonical element, confirmed by
recomputing it from the stored exponents) give p^k1(alpha1 + beta1 x) =
p^k2(alpha2 + beta2 x) mod r, i.e. x = (A1 - A2)/(B2 - B1) with the stored canonical
exponents. B1 = B2 has probability 1/r per collision for random starting exponents.
Verify g^x = z (<= 510 units) before output.

**Fruitless cycles (rule).** Each step compares K with the key two steps earlier;
on equality the two classes {C1, C2} form a 2-cycle and the walk squares the stored
representative of the one with smaller K (y -> y^2 maps a class to a class:
(y^{p^k})^2 = (y^2)^{p^k}), doubling its exponents. A checkpoint key is stored every
64 steps and compared at every step; a repeat means a cycle of length L <= 64, which
is traversed once more (L steps) to find its minimum-K class, which is squared. Both
escapes depend only on the set of classes in the cycle, not on where or when the walk
entered it, so two walks that have merged leave the cycle at the same class and stay
merged. Detection of an L-cycle (L <= 64) happens within 64 + L steps of entering it,
so the waste per event is <= 64 + 2L + 1 <= 193 <= 200 steps; a 2-cycle event wastes
<= 3 steps.

**Stuck-walk rule.** Each walk counts steps since its last distinguished class (reset
at each report and at restart). When the count exceeds 40/theta = 2^65.33 the walk is
abandoned and restarted at a fresh g^alpha z^beta. Its steps are real steps and are
inside the walk-step cap of section 5. It is only a safety net for cycles longer than
64 (section 3.2 bounds their total probability by 2^-100). Under H1 the steps between
distinguished classes are geometric with parameter theta, so a run exceeds 40/theta
with probability (1 - theta)^(40/theta) < e^-40; over at most 2^66 runs the expected
number of restarts is < 2^66 e^-40 < 2^8.3, and by Markov's inequality more than
2^30 restarts occur with probability < 2^-21. 2^30 restarts waste at most
2^30 x 2^65.33 = 2^95.33 steps, a fraction < 2^-30 of S.

### 3.2 Fruitless-cycle rates (derivation)

Let a walk hold y_0 and take L steps with indices h_1..h_L and twists t_1..t_L. Then
y_L = y_0 * prod_i M_{h_i}^{p^{t_i}}. The walk is back in the class of y_0 with the
same exponent pair (a fruitless cycle, i.e. one that occurs for every x, not by the
1/r chance of a random mapping) iff y_L = y_0^{p^e} as formal expressions in the
independent uniform exponents (alpha_0, beta_0, a_h, b_h). Comparing coefficients
forces p^e = 1 and, for every index value h that occurs,
sum_{i: h_i = h} p^{t_i} = 0 mod r. Two consequences:

* A single occurrence is impossible (p^t != 0), so every index that occurs appears
  >= 2 times; the L steps use d <= floor(L/2) distinct indices.
* For a block of k >= 2 occurrences, sum p^{t_i} = 0 mod r. p mod r is a primitive
  12th root of unity, whose minimal polynomial is Phi_12, so this is the event that a
  sum of k 12th roots of unity vanishes in Z[x]/Phi_12 (for k < r^{1/4} a nonzero sum
  has norm < r, so "= 0 mod r" equals "= 0"). With uniform twists this is the return
  probability of a symmetric random walk; it is 1/12 at k = 2 and, by the identity
  P(S_{2m} = 0) = ||mu^{*m}||_2^2 (non-increasing in m) and Cauchy-Schwarz for odd k,
  at most 1/12 for every k >= 2 (part C computes k = 2..10 exactly: max 1/12).

Under H1 (indices and twists of successive distinct classes behave as independent
uniform values, the same assumption that gives the 2-cycle rate 1/(12R)), the
probability that the walk enters a fruitless L-cycle at a given step is

  P_L <= sum_{d=1}^{min(R, floor(L/2))} C(R, d) (d/R)^L 12^-d,

Here D is the number of distinct indices among the L steps; the cycle needs D = d
with every one of the d blocks vanishing, so P_L <= sum_d P(D = d) 12^-d, and
P(D = d) <= P(all L indices lie in some d-set) <= C(R, d) (d/R)^L (union bound).
Computed exactly in part C of experiment.py for R = 1024 (part C also prints
R = 64 and 256: 0.1422 R^-2 and 0.1397 R^-2, used for the toy comparison in section 7):

* sum_{L=3..64} P_L = 0.139 R^-2 <= R^-2 (L = 3 alone is R^-2/12, L = 4 is 0.056 R^-2).
* sum_{L=65..2R} P_L <= 2^-237.3. For L > 2R every term with d < R shrinks by d/R
  per extra step, so sum_{L>2R} <= R * P_{2R} + 2^127 * 12^-1024 < 2^-226.

Hence the long-cycle checkpoint fires at rate <= R^-2 per step, and over all
<= 2^126 steps a fruitless cycle longer than 64 occurs with probability < 2^-100.
Squaring steps cannot be part of a fruitless cycle (2^t is not in <p> mod r for small
t). Toy runs (section 7) measure long-cycle rates below the charged R^-2 and
consistent, within Poisson noise, with the derived bound.

## 4. Success probability

Heuristic (H1): the class sequence of the walk behaves like a random mapping on
N classes, with the Brent-Pollard correction for R-adding walks (effective
N' = N * R/(R-1)). Then the probability that no useful class collision occurs among
the first S distinct-step points of all walks is exp(-S^2/(2N')). Choose

S0 = sqrt(2 N' ln(1/0.61)) = 2^125.6285,   S = 1.03 * S0 = 2^125.6711,

where the 3% is explicit slack for walk non-randomness; at S0 the model gives
exactly 1 - 0.61 = 0.39, and at S it gives 1 - 0.61^(1.03^2) = 1 - 0.61^1.0609 =
0.4081. The algorithm stops after S walk steps (plus the detection window below), so
S is a fixed budget, not an expectation. Toy evidence (section 7) measures the
fraction of instances with a useful collision within S0 and within 1.03 S0 using the
production predicates.

Detection: after the first useful collision, the merged trajectory reaches a
distinguished class within 20/theta further steps except with probability e^-20.
Charging all M walks in lockstep over that window costs E = M*20/theta = 2^84.3
steps. The losses, each a probability bound: detection tail e^-20 < 2^-28; trailing
walk restarted before reaching the distinguished class (its distinguished-free run
must exceed 40/theta: under H1 the run containing a given step exceeds 40/theta with
probability (1 + 40) e^-40 < 2^-52); self-collisions of a single walk (fraction
<= 2^-19 of collisions, with M = 2^20 walks); B1 = B2 (1/r); key ties (< 2^-249);
long fruitless cycles (< 2^-100); more than 2^30 restarts (< 2^-21, section 3);
more than 2^40 false truncated-key matches (< 2^-37, section 3); a cycle-handling
step count above its budget (< 2^-40, section 5). Their sum is
< 2^-19 + 2^-21 + 2^-28 + 3 x 2^-37 < 2^-18 < 0.4081 - 0.39 = 0.0181. Success
>= 0.4081 - 2^-18 > 0.4081 - 0.0001 = 0.4080 >= 0.39.

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

Fruitless-cycle overhead factor: 1 + 3/(12R) + 200/R^2 = 1.000435 (2-cycles at rate
1/(12R), <= 3 steps each; long cycles charged at rate R^-2, <= 200 steps each).
The walk-step cap of the algorithm is 1.000435 S steps in total, cycle handling
included, so this factor is a hard cap, not an average. By section 3.2 the actual
mean cycle-handling fraction is at most 3/(12R) + 200 x 0.139/R^2 = 0.000271 <
0.000435. The cycle events are sums over ~2^125.7 steps of rare events (means >=
2^100), so by a Chernoff bound the cycle-handling steps exceed 0.000435 S - 2^-30 S
(the second term reserves room for restart waste, section 3), leaving fewer than S
useful steps, only with probability < 2^-40 (counted in section 4).

Setup and output, all <= 2^31 units together: pairing g (Miller loop + final
exponentiation, < 2^12 units, bounded by 2^14), input check z^r = 1 (<= 510),
table: 1024 multipliers x (2 x 255-bit exponentiations <= 1,021 units + 11
Frobenius maps <= 1 unit each) <= 2^20.02, M = 2^20 random walk starts
g^alpha z^beta <= 2^30.0, distinguished-point handling (2^65.7 reports x <= 1,000
instructions: hash-table insert, 2 modular multiplications by p^kappa) <= 2^61.3,
output [x]P0 (<= 510 G1 operations, each < 1 unit) and verification pairing
(<= 2^14). Total bounded by 2^31 + 2^62.

W = c * S * 1.000435 + c * E + 2^62
  = 1.14 * 2^125.6711 * 1.000435 + 1.14 * 2^84.3 + 2^62
  = 2^125.86080 (log2 1.14 = 0.18903, log2 1.000435 = 0.00063;
    125.67114 + 0.18903 + 0.00063 = 125.86080; the other two terms add < 2^-41).

Submitted bound: **time_log2 = 125.861** = 125861/1000 >= 125.86080 (part C of
experiment.py recomputes log2 W = 125.8608 and checks claim >= log2 W).

Memory: S * theta = 2^65.67 distinguished reports x 80 bytes (16-byte truncated key,
bits 70..197 of K, + two 32-byte exponents, section 3), at hash-table load 2/3: 2^65.67 x 80 x 1.5 =
2^72.58 bytes, plus 7.1 MB table and 2^20 walk states (~1 KB each, 2^30 bytes):
< 2^73 bytes. Memory is reported, not scored.

## 6. Heuristics and limitations

* H1 (random-mapping behaviour of an R-adding walk on Frobenius classes) is the
  standard heuristic behind every rho estimate (Teske 1998/2001; Brent-Pollard;
  van Oorschot-Wiener 1999 for parallel distinguished points; Duursma-Gaudry-Morain
  1999 and Bos-Kleinjung-Lenstra 2010 for automorphism classes and fruitless
  cycles). It is supported here by the toy runs; it is not proved.
* Distinguished-point and partition bits are taken from K; K is the minimum of
  12 values, which is not uniform in [0,p), but its low 70 bits are pseudo-random
  functions of the class (K is about p/13 > 2^377, so the bits used are far below
  its size; heuristic, standard). The toy runs use exactly this rule on the raw toy
  key (K about 2^16.3 there, bits used <= 14) and measure the distinguished fraction
  equal to its nominal value (section 7).
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

Part B (toy, same tower shape): q = 7 mod 12 near 2^20 with xi = 1+u neither a square
nor a cube (asserted), r' a prime factor of q^4 - q^2 + 1 with ord_r'(q) = 12
(asserted), N = (r'-1)/12. It runs the production walk: four-multiplication keys,
kappa, twisted table, partition h = K mod R and distinguished test "bits
log2 R .. log2 R + dp_bits - 1 of K are zero" on the raw key, reports stored by key
and confirmed by the canonical element, 2-cycle escape, checkpoint escape, the
stuck-walk rule at 40/theta, collision solving and verification. Toy-only
differences: the 60 distinguished bits become dp_bits (the toy key has ~16 bits); ties
of the minimum key (about 66/q per step in the toy) are broken by the full conjugate.
It reports the steps to the first useful class collision against S0 and 1.03 S0 for
the 0.39 target, the measured 2-cycle and long-cycle rates and the distinguished
fraction.

Part C recomputes, without randomness, the zero-sum probabilities, the fruitless-cycle
bound of section 3.2 and the W ledger of section 5.

Offline runs of the same file (seeds 101-103, 201-203, 301-302, "seconds": 240,
"dump": true; 3,079 instances, every log recovered and verified; S0 is the 0.39
budget, model 0.39 at S0 and 0.4081 at 1.03 S0):

| q, r', R, dp_bits | trials | mean T_c/sqrt(N) (pred.) | P[T_c <= S0] | P[T_c <= 1.03 S0] | 2-cycles/step (1/(12R)) | long-cycle events, per step / R^-2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1048963, 48674953, 64, 4 | 1200 | 1.240 +- 0.018 (1.263) | 0.393 +- 0.014 | 0.408 | 1.303e-3 (1.302e-3) | 33, 0.044 |
| 1049611, 64668529, 256, 4 | 1200 | 1.236 +- 0.018 (1.256) | 0.409 +- 0.014 | 0.427 | 3.11e-4 (3.26e-4) | 9, 0.167 |
| 1054003, 35456053, 1024, 3 | 679 | 1.252 +- 0.027 (1.254) | 0.406 +- 0.019 | 0.433 | 8.8e-5 (8.1e-5) | 1, 0.71 |

Measured distinguished fractions: 0.0628 and 0.0642 (nominal 2^-4 = 0.0625), 0.130
(nominal 2^-3 = 0.125). No walk was abandoned by the stuck rule. Every measured
long-cycle rate is below the R^-2 charged in section 5. The derived bound 0.14 R^-2
predicts at most 107.7, 7.5 and 0.20 events for the three rows; 33, 9 and 1 were seen
(the last two are within Poisson noise of the bound: P[Poisson(7.5) >= 9] = 0.34,
P[Poisson(0.20) >= 1] = 0.18); the cost ledger charges the 7x larger rate R^-2.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999);
Teske, Math. Comp. 70 (2001); Duursma, Gaudry & Morain, ASIACRYPT 1999; Wiener &
Zuccherato, SAC 1998; Bos, Kleinjung & Lenstra, ANTS 2010.

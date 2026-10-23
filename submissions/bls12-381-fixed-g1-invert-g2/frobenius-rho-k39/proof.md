# Frobenius-class parallel Pollard rho in GT (fixed G1 argument, output in G2)

Claim: time_log2 = **125.86** (total work W <= 2^125.86 dense-Fp12-multiplication
equivalents, C = 20,000), success probability >= 0.39, peak memory < 2^75 bytes,
preprocessing <= 2^31 units (included), no advice. Heuristics are listed in
section 6. Every number below is reproduced as an exact rational by part C of
`experiment.py`, including the final integer certificate that W < 2^125.86 and
that 125.85 is not enough.

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
  **4 Fp multiplications** (a*zeta, a*Re(d_1), b'*Im(d_1), B*zeta) plus one Fp
  addition and 8 subtractions/negations. No Frobenius of the whole element is
  ever computed.

Define kappa(y) = argmin_k key_k (keys as canonical residues in [0,p), compared as
integers) and the class key K([y]) = min_k key_k. Then y^{p^{kappa(y)}} is the same
element for all 12 members of a class (canonical representative y^), and K is a
class invariant. A tie needs an equality of two of the 12 keys, which for a
uniform a (or B) is a linear coincidence modulo p: at most 66/p per step, and
at most 66 * L_max / p <= 2^-248.9 for the whole run (log2 p = 380.7007). A tie is
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
compute keys, kappa, K                      (4 Fp mul + <= 2,600 other instructions)
if bits 10..69 of K are all zero (distinguished): report (K, p^kappa*alpha, p^kappa*beta)
h = K mod R (bits 0..9) ; t = (12 - kappa) mod 12
y = y * T[h][t]                             (ONE dense Fp12 multiplication)
alpha += a_h p^t ; beta += b_h p^t (mod r)  (2 additions mod r)
```

So each step costs exactly one reference operation plus the key work. The walk is
a deterministic function of the class, so once two walks reach the same class they
follow the same class sequence and reach the same distinguished class. A walk
continues after reporting; it is never restarted after a distinguished point.

**Distinguished store (batched, deterministic).** The walks never consult a
store. Each report is appended to a flat array: the canonical class key K in six
64-bit limbs (48 bytes, K < p < 2^381) and the canonical exponents A, B as 32-byte
residues mod r, 112 bytes per record. A global counter also caps the array at
2^67 records; the algorithm stops with failure at the next report. When the step
budget is spent, the array is ordered by K with a fixed LSD radix sort: 24 stable
passes of 16 bits (K padded to 384 bits), each pass a counting sweep and a scatter
sweep over one 2^16-entry counter array, permuting an array of 12-byte (96-bit)
indices that is read from one buffer and written to a second buffer of the same
size (the two buffers swap roles every pass; the first is filled with the
identity order). The report cap 2^67 exceeds one 64-bit word, so every counter,
offset and byte address of the sort is a 128-bit quantity (two words) and the
indices are stored in three words. The pass count, the per-pass instruction count
and the counter array are constants, so the handling work is a fixed function of
the record count: there is no probing, no load factor, no rejected insert and no
discarded report (section 5). One linear scan then reads the sorted indices in
order, keeps the first record of each run of equal K, and compares every later
record of the run with it and with its predecessor. If the exponents equal the
stored ones the record is dropped (same walk history); otherwise
x = (A1 - A2)/(B2 - B1) mod r is computed and checked against g^x = z. On success
the algorithm outputs Q = [x]Q0; on failure it keeps the first record.

Detection is complete and exact for equal K: two reports of one class carry the
same K, because K is a class invariant, so they lie in one run, and every record
of a run is compared with the run's first record, so the pair is tested and no
report is ever skipped. A run can mix two distinct classes only if those classes
share K; the check g^x = z then rejects any comparison that includes the wrong
class (section 6, H2). The relation is exact for one class: with equal K,
p^k1(alpha1 + beta1 x) = p^k2(alpha2 + beta2 x) mod r, and the stored exponents
are the canonical ones. B1 = B2 has probability 1/r per collision. Detection
after the fact cannot miss a useful collision either: no walk ever reads the
store, so no walk stops or diverges early and each walk takes exactly the steps
charged in section 4.

**Hard report cap:** a global counter counts every report, stored or matched;
the algorithm processes at most 2^67 reports and stops with failure at the next
one. Hence at most 2^67 records are ever stored, and no report is lost to the
store before that cap.

**Hard global step cap:** every step of every kind (ordinary steps, squaring
escapes, cycle traversals, steps of walks that are later abandoned) increments one
global counter. At L_max = 1.003 S + E = 2^125.6755 steps (S, E in section 4) the
algorithm stops with failure.

**Fruitless cycles.** A twisted adding walk can return to a class with formally
cancelled exponents. Exact 2-cycles need h' = h and kappa' = 6 (inverse conjugate):
rate 1/(12R) = 2^-13.58 per step (toy measurement: 1.303e-3 vs predicted 1.302e-3
at R = 64). Each step compares K with the key two steps earlier; on equality the
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
restarts number at most L_max / 2^65.32 = 2^60.35, deterministically.
Squaring steps cannot be part of a fruitless cycle (2^t is not in <p> mod r for
small t).

## 4. Success probability

Heuristic (H1): the class sequence of the walk behaves like a random mapping on
N classes, with the Brent-Pollard correction for R-adding walks (effective
N' = N * R/(R - 1), log2 N' = 251.2735). Then the probability that no useful class
collision occurs among the first S points of all walks is exp(-S^2/(2N')). Choose

S0 = sqrt(2 N' ln(1/0.61)) = 2^125.6285,   S = 1.03 * S0 = 2^125.6711,

where the 3% is explicit slack for walk non-randomness. At S0 the model gives
exactly 0.39; at S it gives 1 - 0.61^(1.03^2) = **0.408**. The track requires
0.39 at the submitted bound, and the losses below are charged against the 0.018
of margin between 0.408 and 0.39.

Detection: a class collision is only visible once both merged walks report the
same distinguished class. After the first useful collision the merged trajectory
reaches a distinguished class within 40/theta further steps except with
probability e^-40, and both walks then report that class, so the pair appears in
the array and the final scan compares it. Charging all M walks in lockstep over
that window costs E = M*40/theta = 2^85.32 steps. This is why the step cap is
1.003 S + E: after paying the detection window and at most 0.003 S of waste the
run has S useful steps.

Step cap: useful (distinct-point) steps are at least L_max - E - waste. Under H1
the waste is at most (3/(12R) + 200/R^2) = 0.000435 of all steps, so it stays
below 0.003 S except with negligible probability, leaving at least S useful steps
before the cap.

Report cap: under H1 each step lands on a distinguished class with probability
theta, so reports among at most L_max steps have mean at most 2^65.71; reaching
2^67 (2.46 times the mean) has Chernoff probability below exp(-2^65), negligible.

Abandoned walks: a 40/theta stretch with no distinguished point has probability
e^-40 = 2^-57.7. Every stretch begins at a report, a start or a restart, so there
are at most 2^67 + 2^20 + 2^60.35 < 2^67.02 stretches; the expected number of
abandonments is below 2^9.3 and the steps they strand (below 2^74.7) are inside
the 0.003 S waste allowance, since they are counted steps either way.

No report is lost to the store: every report is written to the array before any
comparison, every record of a run is compared against the run's first record and
its predecessor, and the only way a comparison can involve a second class is H2.
The store therefore contributes no success loss; the array and its sort remove
the load-factor and discard analysis of the previous design entirely.

Losses charged against the 0.408: ties 2^-248.9, B1 = B2 at most 2^68/r for the
at most two comparisons per record, a mixed run under H2 at most 2^-247, the
e^-40 detection tail 2^-57.7, and the 1/M = 2^-20 self-collision fraction, which
is what the 3% step slack is for. Every term except the self-collision fraction
is at most 2^-57.7, and the self-collision fraction is inside the slack, so
success is at least 0.408 - 2^-19.9 >= 0.4079 at the model probability, above the
required 0.39.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Per-step non-reference instructions, scalar AArch64, using the published 383
instructions per blst_fp_mul:

| Work | Itemised hard cap |
| --- | ---: |
| 4 Fp multiplications (a*zeta, a*Re(d_1), b'*Im(d_1), B*zeta) at 383 | 1,532 |
| 1 Fp addition (B), 6 limbs with carry and conditional subtract | 24 |
| 8 Fp subtractions/negations, the six zeta^j multiples from a and from B, 24 each | 192 |
| min of 12 six-limb keys with index select: 12 comparisons with 6 loads, 6 subtracts and 2 selects each, plus 16 for the running minimum | 192 |
| distinguished test, partition index, t = (12-kappa) mod 12, table address | 32 |
| alpha, beta += table values mod r: 2 loads of 32 bytes, 2 four-limb adds with conditional subtract | 64 |
| 2-cycle key compare, checkpoint compare, key history store, buffer rotation, step counter | 48 |
| loop control | 16 |
| **total itemised** | **2,100** |
| **budgeted per step** | **2,600** |

The 500-instruction difference between the itemised 2,100 and the budgeted 2,600
is an explicit margin for load and store scheduling, the register traffic of the
table entry and the multiply operands, and any micro-operation the table does not
name. It is a budget, not a measurement: each row above is a worst-case count for
a specified operation, and every row is a hard cap.

The product y*T[h][t] is the dense blst_fp12_mul reference operation (operands
read through pointers, output written to a rotating buffer; no copies), charged 1.
A squaring escape uses one Fp12 squaring (at most one dense multiplication) and
doubles two exponents, inside the same budget. Every counted step therefore costs
at most c = 1 + 2,600/20,000 = **1.13** units, and the walk costs at most c * L_max.

Additive terms (all deterministic bounds, each charged once):

* Setup and output, <= 2^31 together. Pairing g = e(P, Q0): the Miller loop is 64
  iterations of one Fp12 squaring, two line multiplications counted as dense and
  one G2 doubling at <= 1.44 units, plus one G2 addition for each of the 6 set bits
  of the ate loop, so <= 64 x 4.44 + 6 x 1.44 = 292.8 units; the final
  exponentiation even by plain square-and-multiply on the 4,314-bit exponent
  (p^12-1)/r costs <= 8,628 units; the pairing is <= 8,921 units. Table: 1024
  multipliers x (2 x 255-bit exponentiations <= 1,020 units + 1 product + 11
  Frobenius maps <= 1 unit each) <= 1,046,528 < 2^20. M = 2^20 walk starts
  g^alpha z^beta at <= 1,021 units each, including drawing alpha and beta,
  <= 2^30. Output Q = [x]Q0 by double-and-add, <= 510 Jacobian G2 operations,
  each <= 20 Fp2 multiplications (<= 3 Fp multiplications + 5 Fp additions each)
  plus <= 20 Fp2 additions and control, <= 28,780 instructions = 1.44 units, so
  <= 735 units, plus one Fp2 inversion (Fermat, <= 762 Fp multiplications,
  <= 15 units). Verification pairing e(P, Q) <= 8,921 < 2^14. The itemised sum is
  1,071,652,588 < 2^31; the declared cap 2^31 is what enters W.
* Restarts: at most 2^60.35 x 1,021 units <= 2^70.35.
* Distinguished-point handling: at most 2^67 + 1 reports (report cap), stored as a
  flat array of 112-byte records, sorted by K and scanned (section 3). The append
  writes the record: 14 word stores, the 128-bit record address (12), the 128-bit
  bound check (3), the 128-bit counter update (4) and loop control, budgeted 48
  instructions per report. The first index buffer is filled with the identity
  order, 3 stores, a 128-bit counter step and loop control, budgeted 12 per
  record. Every counter, offset and address of the sort is 128-bit, since 2^67
  records and 2^73.8 bytes exceed one 64-bit word. Per record per pass the
  counting sweep reads the 96-bit index (3 loads, 1 widening), forms the record
  address as 112 x index plus the base (12), extracts the 16-bit digit (3), loads,
  adds and stores one 128-bit counter with its address (7) and runs the loop (3):
  29. The scatter sweep reloads the index, address and digit (19), loads the
  128-bit counter (2) and forms its address (1), forms the 12-byte output address
  12 x counter plus the base (10), stores three index words (3), increments and
  stores the counter (4) and runs the loop (3): 42. That is 71 itemised, budgeted
  80 for the pass. The scan reads six key limbs, compares them with the previous
  record, advances the 128-bit run position and compares the 8 exponent words on a
  duplicate; budgeted 96. So one record costs 48 + 12 + 24 x 80 + 96 = 2,076
  instructions. Each 2^16-entry counter array (128-bit counters) is zeroed and
  prefix-summed once per pass, at most 16 instructions per entry, 24 x 16 x 2^16 =
  2^24.6 instructions in total. The whole handling term is
  (2^67 + 1) x 2,076 + 2^24.6 <= 2^78.02 instructions = 2^63.74 units. Neither the
  pass count nor the per-pass cost depends on the data.
* Candidate solves: at most two per record (the run's first record and the
  predecessor are each a comparison), <= 520 units each (inversion mod r,
  exponentiation g^x <= 510, comparison): <= (2^67 + 1) x 1,040 <= 2^77.03.

Additive total T <= 2^31 + 2^63.74 + 2^70.35 + 2^77.03 <= 2^77.05. Each term is
quoted as a power of two rounded up; part C evaluates the same terms exactly.

W <= c * (1.003 S + E) + T
  = 1.13 * (1.003 * 2^125.6711 + 2^85.32) + 2^77.05
  = 2^125.8518 < 2^125.86.

Exact form: with u = 0.494296321814780115... an exact rational upper bound on
ln(1/0.61) from the atanh series (2 atanh(39/161), 40 terms plus a geometric
tail), S0 = sqrt(2N'u) rounded up, and E and the restart count exact integers,
part C of `experiment.py` evaluates A = c * 1.003 * 1.03 and B = c * E + T as
exact rationals and checks (A S0 + B)^100 <= 2^12586, which is true, and
(A S0 + B)^100 > 2^12585, which is also true, so (A S0 + B)^100 <= 2^12585 fails
and 125.86 is exactly the smallest multiple of 0.01 at or above the bound, with
0.0082 bits of slack. Nothing is rounded down anywhere in that chain; the two
exact checks are printed by part C.

Submitted bound: **time_log2 = 125.86**.

Memory (hard bound by the record cap). Every array live at the peak, which is the
scatter phase of the radix sort (the walk states are never needed after the step
budget is spent but are counted anyway):

| Array | Bytes |
| --- | ---: |
| report records, 2^67 x 112 | 16,528,282,690,043,758,247,936 |
| scatter input index buffer, 2^67 x 12 | 1,770,887,431,076,116,955,136 |
| scatter output index buffer, 2^67 x 12 | 1,770,887,431,076,116,955,136 |
| digit counters, 2^16 x 16 | 1,048,576 |
| multiplier table T, 12,288 x 576 | 7,077,888 |
| exponent tables a_h p^k and b_h p^k, 12,288 x 64 | 786,432 |
| walk states, 2^20 x 2,048 (element, two history elements, keys, exponents, checkpoint, counters) | 2,147,483,648 |
| scalars, caps, checkpoints, code and stack slack | 1,048,576 |
| **peak** | **20,070,057,552,198,148,554,752** |

The peak is 136 x 2^67 + 2,156,396,544 bytes, 2^74.0875 bytes. The records and
both index buffers are 136 x 2^67 bytes; the rest is below 2^31.01 bytes, which is
below 2^67. Hence peak <= 137 x 2^67 bytes, and since 137^10 < 2^71 the exact
inequality is peak < 2^74.1 < 2^75 bytes, so the declared bound of 2^75 holds
loosely. No array exceeds 2^128 bytes, so every offset and counter is charged at
128-bit width (two words), and the record index, which needs 67 bits, is stored in
96 bits (three words). This is a hard bound by the report cap. Memory is reported,
not scored.

## 6. Heuristics and limitations

* H1 (random-mapping behaviour of an R-adding walk on Frobenius classes) is the
  standard heuristic behind every rho estimate (Teske 1998/2001; Brent-Pollard;
  van Oorschot-Wiener 1999 for parallel distinguished points; Duursma-Gaudry-Morain
  1999 and Bos-Kleinjung-Lenstra 2010 for automorphism classes and fruitless
  cycles). It is supported here by the toy runs; it is not proved.
* H2 (K separates the visited classes): K is the minimum of 12 values, which is
  not uniform in [0,p), and two distinct classes could share a K. There are at most
  2^134 pairs of reports and each shares a 381-bit K with probability at most
  2^-381 under H2, so a run that mixes two classes has probability at most 2^-247
  over the whole run. This is the only role K plays in detection: a comparison
  that includes the wrong class is rejected by the check g^x = z, so H2 cannot
  turn a false pair into a wrong output.
* Distinguished-point and partition bits (bits 0..69) are taken from K, which is
  the minimum of 12 values and not uniform in [0,p); these bits are treated as
  pseudo-random functions of the class (heuristic, standard). The report-count
  argument of section 4 rests on this. The toy does not test it (see section 7).
* The waste, record-count and abandonment estimates in section 4 use H1; the caps
  themselves, the restart count, the per-record instruction cost and every cost
  term in section 5 are deterministic.
* No production run. Toy experiments are evidence only at their own scale.
* This is a generic attack (about the square root of the class count). It ignores
  index-calculus methods in Fp12 (TNFS family); it is a valid upper bound, not a
  security estimate.

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
mapping predicts mean 1.2533), the fractions of trials within the 0.39 and 0.40
step budgets, and the measured 2-cycle rate against 1/(12R). The toy does not
apply the two production caps; its counter counts walk steps, not setup work. The
toy takes its partition index and distinguished test from a hash of the full
canonical element, not from bits of K: at q = 163 the key K is a residue below 163
(about 7 bits), fewer than the 10 partition bits needed. So the toy tests the class
map, the kappa selection, the twisted table and the escapes, but not the
pseudo-randomness of production K's bits, which remains the heuristic of section 6.

Part C (exact ledger): recomputes N, N', the exact rational upper bound on
ln(1/0.61), S0, S, L_max, the setup, record-array, restart and solve terms, W and
the memory total from exact integers and fractions, prints log2 W to six decimals,
and checks the final integer inequality (W^100 <= 2^12586 and W^100 > 2^12585).
It is the machine-checked version of sections 4 and 5.

The toy store replay uses the production store shape and nothing else: every
distinguished report of a trial is appended to a flat array, the array is ordered
by a fixed LSD radix sort (192-bit toy class serialisation, 8-bit digits, 24
passes) and one linear scan compares the records of each run of equal class keys.
Its first 20 trials recover the same logarithm from the sorted array that the
online walk finds, 20 of 20, over 3,538 records.

The initial input (q = 163, R = 64, 4 walks, seed 1, 40 trials, 14 s limit) is
deterministic given its seed, takes about 3 s of CPU and completes all 40 trials:
every logarithm was recovered and verified; the mean first-collision time is
1.2319 sqrt(N) (model 1.2533); the fraction within the 0.39 step budget is 0.400
and the fraction within the 0.40 budget is 0.400 (models 0.39 and 0.40); 2-cycles
run at 1.294e-3 per step (predicted 1.302e-3); there are 139 two-cycle escapes, no
longer cycle, 5 ignored repeat reports and no abandoned walks; the store replay
recovers 20 of 20 trials. The 0.39 budget is 1.0021 sqrt(N) and the 0.40 budget
1.0188 sqrt(N). A runner slower than about three times this one hits the 14 s
limit, completes fewer trials and reports different sample values; the ledger of
part C is unaffected either way.

Larger offline runs of this same file (time limit lifted, so these inputs are
outside the runner's 14 s schema; production_checks 0), all logs recovered and
verified; standard errors from the random-mapping spread 0.655/sqrt(n) and
sqrt(0.24/n):

| q, r', R, dp_bits, walks, seed | trials | mean T_c/sqrt(N) (model) | P[T_c <= S0] (model 0.40) | 2-cycles/step (1/(12R)) |
| --- | ---: | ---: | ---: | ---: |
| 163, 54298861, 64, 4, 4, 1 | 2400 | 1.2622 +- 0.013 (1.2632) | 0.400 +- 0.010 | 1.303e-3 (1.302e-3) |
| 787, 64571173, 256, 4, 4, 2 | 600 | 1.1951 +- 0.027 (1.2558) | 0.427 +- 0.020 | 3.29e-4 (3.26e-4) |
| 331, 923346397, 1024, 6, 8, 3 | 96 | 1.1157 +- 0.067 (1.2539) | 0.490 +- 0.050 | 9.3e-5 (8.1e-5) |

The three larger runs are quoted at the 0.40 step point, S0(0.40), which is the
column the same file has always printed; the 0.39 point is 0.9838 times smaller
and the initial-input run above reports both fractions (0.400 each, models 0.39
and 0.40). The claim of this package is the 0.39 point S0(0.39) = 2^125.6285 with
the 1.03 slack on top, so what these runs evidence is the random-mapping model
itself, at toy scale, not the value of one quantile.

The two smaller samples lie about two standard errors below the model mean, in the
direction favourable to the claim; the large sample matches it. Longer fruitless
cycles were rarer than the R^-2 bound (56 events in about 6.6e6 steps at R = 64,
i.e. 8.4e-6 per step vs bound 2.4e-4). No walk was abandoned in any run.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999);
Teske, Math. Comp. 70 (2001); Duursma, Gaudry & Morain, ASIACRYPT 1999; Wiener &
Zuccherato, SAC 1998; Bos, Kleinjung & Lenstra, ANTS 2010.

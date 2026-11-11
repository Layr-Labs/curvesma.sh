# BN254 pairing inversion: Frobenius-class rho on GT with a planned search phase, output ([x]P0, Q0)

**Claim.** For uniform z = g^x in the protected BN254 target group (g = e(P0, Q0) from
`target-profile.json`, x uniform in [0, r)), the algorithm below returns P = [x]P0 in G1 and
Q = Q0 in G2 with e(P, Q) = z. It succeeds with probability at least 0.39 and uses at most
**2^124.65** reference units on every run (W = R/20,000, R the total scalar AArch64
instruction-equivalent work). This total includes setup, table construction, walk starts,
the planned search phase, the bounded detection tail, distinguished-point (DP) handling,
cycle escapes, restarts, false-candidate verification and output. Peak memory is below
2^63 bytes. Preprocessing is included and is below 2^43 units. Nonuniform advice is 0 bytes.
The success bound rests on the standard random-mapping heuristic for rho walks, adjusted for
the branch mix (section 4). It is supported by toy-field runs of the inherited walk (which
uses run-until-found control, not the planned phase of this package) and by running the
production step rule in the real BN254 GT. This is an analytic generic-attack bound, not a
practical break, and no production-size run is claimed.

**What is new against the rejected predecessor.** The walk, keys, table, prices and the
nominal ledger come from the pairing-inversion package at source
6e57b2d58874274941a2b3173fd5fb331e3e3fb4 (submission 3325b532), which was rejected. That
package's GT twin (submission 2ece9f04, source 6791c3be) has a favourable review at 124.65
but is not promoted. Review is not promotion, and neither package is relied on as accepted.
The rejection found one decisive defect. The stopping rule was "run until recovery or a cap",
and H1+ claimed that any cap or abort fires with probability at most 0.01. Yet the nominal
S steps have collision probability only 0.43, so about 0.56 of runs must end at a cap. This
package changes the algorithm and the claim, not only the label:

1. **The search phase is fixed in advance** (section 3.4). Every walk takes exactly K planned
   steps (K = ceil(S(1+eps)/M)). It then runs a bounded detection tail: it continues only
   to its next DP, and for at most 20/theta steps. Then the run matches reports and stops.
   A run that finds no match after the tail ends by **normal no-hit completion**. That is
   the ordinary failure event of probability about 0.57. It is neither a cap nor an abort.
2. **Premature truncation is a separate event.** It means that the hard work meter reaches
   2^124.65 before the planned phase and tail are complete. It happens only if the realised
   multiply-step fraction over at most S_plan + E steps exceeds 1/2 by more than 1.5700
   points. Section 3.5 bounds it, and the bound comes from the same sampling premise as the
   collision law. The failure budget charges 0.001 for it.
3. **The output is ([x]P0, Q0)** instead of (P0, [x]Q0). This is a left-to-right
   double-and-add on G1 over Fp, with counted Fp operations, checked against an independent
   affine reference by `planned_phase.py`. The analytical bound is 1,005,133 instructions
   (including the loop and serialisation), charged 101 units = 2,020,000 instructions. The
   former G2 output was charged 2^10 units. If z = 1, the output is
   (infinity, Q0).

`experiment.py`, `machine.md` sections 1-4 and `development-evidence.json` are inherited
unchanged, except section 4a and section 5 of `machine.md`. `planned_phase.py` is new and is
the declared initial experiment. It does not run the rho. It does four things:
- recomputes the full-run ledger with the planned phase, the tail, the report store, the
  truncation margin and the success budget (part P);
- executes the finite-width control primitives limb by limb at boundary values (part Q);
- executes the counted G1 output (part O).

**Revision 2 (this commit).** These changes fix finite-width and scope defects found by an
independent critique:
- 128-bit planned and tail counters;
- a 192-bit instruction meter against floor(cap), checked before each priced block;
- DP density 2^-70, so that the report store fits one 64-bit address space;
- an explicit node and bucket layout with a capped chain search;
- the analytic output bound separated from the sample maximum;
- corrected evidence scopes.
The claimed score is unchanged.

**Lineage.** The design is the Frobenius-class rho of the BLS12-381 fixed-argument track.
The class construction, 4-multiplication key and twisted multiplier table come from
submission 7da5f23a-6633-46b1-9c03-11f6cfedeeec. The cyclotomic-squaring branch, read from
the trace invariant, comes from 6d4e57d1-135f-4d1f-9d43-eca27fb5c9ed and
747ad6c9-0453-4e37-8999-4f02a32c3c8c. The stored-representation rule, work meter and
charging at the counted listing come from 56fc259f-d5d1-4c63-b4e5-9d95e6188bd6. The BN254
port and its ledger come from source 6e57b2d5 (submission 3325b532) and its GT twin 2ece9f04.

## 1. Problem and output

Input z in GT. If z = 1, output (infinity, Q0). Infinity is allowed exactly when z = 1, and
Q0 is in G2. Otherwise find x with g^x = z and check it. Since g has prime order r (part A,
`g_order_r`), x is unique in [0, r) and x != 0. Output P = [x]P0 in affine coordinates and
Q = Q0. By bilinearity e([x]P0, Q0) = e(P0, Q0)^x = g^x = z. P lies in G1 because P0 has
order r (`planned_phase.py` part O checks [r]P0 = infinity) and the G1 cofactor is 1. P is
not infinity because 1 <= x < r. The check g^x = z is the attack-side verification (one
exponentiation, under 400 units); no pairing is computed.

**G1 scalar multiplication.** The scalar is x, with its top bit at position at most 253.
Start from acc = P0. For each later bit, double acc in Jacobian coordinates (dbl-2009-l with
a = 0: 2M + 5S), and if the bit is 1, add P0 by the mixed formula (madd-2007-bl: 7M + 4S).
Then perform one Fermat inversion of Z (253 squarings and 128 products) and the affine
conversion (3 products). There are no exceptional cases. Before each doubling,
acc = [k]P0 with 1 <= k < r/2, so it is neither infinity nor 2-torsion, because r is odd.
Before each addition, acc = [2k]P0 with 2 <= 2k <= x - 1 <= r - 2, so acc != +-P0.
Part O of `planned_phase.py` counts every Fp operation of this sequence at the leaf prices
of `machine.md` section 2 (162 / 26 / 27, routine calls 32):

- **Analytical bound.** One doubling costs 1,535 instructions and one mixed addition 2,186,
  for every input. The inversion and the affine conversion cost 59,324. Since x < r < 2^254,
  there are at most 253 doublings and 253 additions, which gives 1,000,737 instructions.
  Adding 254 x 12 for bit handling and 1,348 for leaving Montgomery form and writing
  64 hex digits per coordinate gives 1,005,133. It is charged **101 units = 2,020,000
  instructions**.
- **Sample.** Part O also runs 73 scalars: 1, 2, 3, r-1, r-2, (r-1)/2, 2^253-1,
  3*2^252-1 = 0x2fff...ff, 2^253 + 2^252 - 1 and 64 random scalars. All 73 match an
  independent affine reference, lie on the curve and stay at or below the analytical bound.
  The sample maximum is 998,551, at 3*2^252-1. The sample is a check, not the bound.

## 2. Arithmetic facts (checked in `experiment.py`, parts A and C)

Tower: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi), xi = 9 + u. With u = w^6 - 9 this is
the profile field Fp[w]/(w^12 - 18 w^6 + 82) (part A checks the map is a ring
homomorphism on random elements and that the profile g has order r). p = 7 mod 12; xi is
neither a square nor a cube in Fp2.

**F1 Frobenius classes.** (sum c_m w^m)^p = sum conj(c_m) d_m w^m with d_m = xi^{m(p-1)/6}.
p has order exactly 12 modulo r, p^6 = -1 and p^4 - p^2 + 1 = 0 (mod r), so y -> y^{p^k}
is a group of order 12 acting on GT and containing inversion. Every y != 1 lies in a class
[y] of exactly 12 elements: N = (r-1)/12 classes, log2 N = 250.0117.

**F2 Keys.** zeta = xi^{(p^2-1)/6} lies in Fp and is a primitive 6th root of unity. Let
c_1 = a + b u be the w^1 coefficient of y and B = a Re(d_1) + b Im(d_1). Then
`key_k := Re(w^1 coeff of y^{p^k}) = a zeta^j (k = 2j), B zeta^j (k = 2j+1)`,
and with zeta^j in {1, zeta, zeta - 1, -1, -zeta, 1 - zeta} all twelve keys cost **4 Fp
multiplications** (a Re d_1, b Im d_1, a zeta, B zeta), one addition and eight
subtractions or negations. Part A checks the formula against the twelve explicitly computed
conjugates (`keys_ok`) and that the twelve keys are distinct (`distinct_ok`). Two keys
coincide only by a linear accident (probability at most 66/p per element; a tie is
charged as failure). kappa(y) = argmin_k key_k, and y^{p^kappa} is the same element for all
twelve members of the class (`canon_ok`).

**F2a Representation.** The kernels hold Fp elements in Montgomery form x~ = x 2^256 mod p
(`machine.md` section 1). The rule reads the **stored** limbs: kappa is the argmin of the
twelve stored keys, and T below is the stored value. Montgomery products of stored operands
are stored products and negation commutes with the bijection x -> x~, so the stored keys
are exactly the stored images of the keys; the rule is the integer rule composed with a fixed
bijection of Z/p and every structural property carries over. Part A checks that the stored
argmin selects one identical canonical element from all twelve conjugates
(`stored_canon_ok`). No conversion out of Montgomery form occurs in the charged path.

**F2b Trace invariant.** The w^0 coefficient of y^{p^k} is conj^k(c_0) (no twist), so
T([y]) := stored Re(c_0) is the same for all twelve conjugates (`trace_inv_ok`). T is read
directly from memory, with no multiplication, and for distinct classes behaves as an
essentially uniform residue.

**F3 Cyclotomic squaring (Granger-Scott).** GT lies in the cyclotomic subgroup. With
A = (c0, c3), B = (c1, c4), C = (c2, c5) in Fp4 = Fp2[s]/(s^2 - xi), s = w^3, and
(A0, A1) = A^2, (B0, B1) = B^2, (C0, C1) = C^2:
`c0' = 3A0 - 2c0, c3' = 3A1 + 2c3, c1' = 3 xi C1 + 2c1, c4' = 3C0 - 2c4, c2' = 3B0 - 2c2, c5' = 3B1 + 2c5.`
Part A checks y^2 = formula on GT elements (`sq_ok`) and counts 18 fp_mul, 96 fp_add, 31
fp_sub and 13 routine calls (`machine.md` section 3).

**F4 Orders.** r - 1 = 2^28 * 3^2 * 13 * 29 * 983 * 11003 * 237073 * 405928799 *
1670836401704629 * 13818364434197438864469338081 (verified by multiplication). The order of
2 modulo r is about 2^247.31, so 2^t lies in the order-12 group <p> only if t is a multiple
of about 2^243.7.

## 3. Algorithm

Parameters: R = 1024 table multipliers, squaring probability s = 1/2, DP density
theta = 2^-70, M = 2^20 parallel walks. From stored T: squaring iff T mod 8 < 4; table index
h = (T >> 3) mod 1024; the class is distinguished iff (T >> 13) mod 2^70 = 0; the fingerprint
used for matching reports is F = (T >> 83) mod 2^128, disjoint from the bits above.
The search is a fixed planned phase of K steps per walk followed by a bounded detection
tail (section 3.4); it is not 'run until recovery'.

Setup: convert z to the tower (once). Draw fresh uniform a_h, b_h mod r and compute
M_h = g^{a_h} z^{b_h} and Tab[h][t] = M_h^{p^t} for h < 1024, t < 12, with the scalar tables
a_h p^t, b_h p^t mod r (12,288 entries of 384 + 64 bytes, 5.5 MB).

A walk keeps a member y of its current class and (alpha, beta) mod r with
y = g^alpha z^beta, starting from fresh uniform (alpha, beta). One step:

```
1  T <- stored Re(c_0) of y; branch, h and DP bits from T
2  if distinguished: compute the 12 keys -> kappa; report (F, p^kappa alpha, p^kappa beta)
   into a chained multimap keyed by F (entries are never replaced)
0  meter pre-check: acc + 18,191 <= floor(cap), else abort (machine.md 4a)
3  2-cycle and checkpoint tests on T (section 3.3); planned/tail counters; charge the
   actual step price (18,191 or 6,965 instructions) to the 192-bit meter
4  if T mod 8 < 4:  y <- y^2 (Granger-Scott); alpha, beta <- 2 alpha, 2 beta   [no keys]
   else: compute the 12 keys -> kappa; t <- (12 - kappa) mod 12;
         y <- y * Tab[h][t] (one dense Fp12 product); alpha += a_h p^t; beta += b_h p^t
```


**3.1 The step is a map on classes.** Write y^ = y^{p^kappa}. Multiply branch:
y * M_h^{p^{-kappa}} = (y^ M_h)^{p^{-kappa}}, which lies in the class [y^ M_h]. Squaring
branch: y^2 lies in [y^^2] because squaring commutes with Frobenius. Both depend only on the
class (T is a class function; kappa and y^ enter only through the conjugation absorbed by the
twist), so walks that reach the same class stay merged and reach the same DP. Part C runs this
rule in the real BN254 GT and checks, for each of the first 150 steps, that the image computed
with the twisted table equals, after the explicit Frobenius p^kappa, y^ M_h or y^^2.

**3.2 Collision.** Two reports with equal F from different walks are the same class, up to a
false match of probability 2^-128 per pair of distinct classes: with at most D = 2^55.44
reports (cap, section 3.5) there are about D^2 / 2^129 = 2^-18.1 false candidates in
expectation, and each candidate is checked (one inverse mod r and g^x = z, under 600 units).
For a true match the canonical exponents satisfy alpha^1 + x beta^1 = alpha^2 + x beta^2
(mod r), so x = (alpha^1 - alpha^2) / (beta^2 - beta^1); beta^1 = beta^2 has probability 1/r.
Output x after checking g^x = z.

**3.3 Fruitless cycles.** After l steps the canonical exponent pair is
Lambda (alpha^_0, beta^_0) + tau with Lambda = 2^t p^m (t squaring steps) and
tau = sum_h c_h (a_h, b_h), c_h in Z[p]. Since the start pair is uniform and independent of
the table, a closed class walk returns the pair unchanged (a useless cycle) with probability
about 1/r unless Lambda = 1 and every c_h = 0.
* **No squaring step occurs in a useless cycle**: Lambda = 1 with t >= 1 needs 2^t in <p>,
  i.e. t >= 2^243 (F4).
* A useless cycle of multiply steps needs every table index used to occur at least twice with
  coefficients summing to zero in Z[p]/(Phi_12(p)); for any number of occurrences the
  probability is at most 1/12. For 2-cycles this is the model rate
  (1-s)^2 / (12R) = 2^-15.6 per step. For l = 2k independent uniform indices (H2) the
  probability is at most (2k)! / (2^k k!) (12R)^-k, below 2^-288 for 64 <= l <= 147; for
  l >= 148 all l steps must be multiply steps, probability 2^-148 per start position.
* **Detection.** If the T of two steps back equals the current T the two classes form a
  2-cycle; the walk squares the representative with the smaller T (exponents doubled). This
  depends on the cycle only, so merged walks stay merged. A checkpoint T is stored every 64
  steps and compared at every step; a repeat means a cycle of length at most 64, which is
  traversed once more to its minimum-T class, which is then squared (waste at most 200
  steps). A walk that sees no DP within 20/theta = 2^74.3 steps restarts from a fresh pair. A restart
  draws a new start (one exponentiation, 400 units, non-walk) and resets only the 128-bit
  stuck counter; it keeps the planned counter.
* Every escape step, cycle-traversal step and step after a restart is an ordinary walk step:
  it is charged at its branch price and increments the planned counter. Escape work
  therefore consumes planned steps and is never extra work beyond K per walk. Its expected
  share is at most eps.
* Overhead factor on the step count, charged: eps = 5/(12R) + 250/R^2 = 6.45e-4.

**3.4 Control flow: planned phase, detection tail, completion.** The constants below are
fixed before any coins are drawn (`planned_phase.py` part P):

* S = sqrt(2 N' ln(1/0.57)) = 2^125.29809 useful steps (section 4) and
  eps = 6.453e-4 (section 3.3).
* K = ceil(S(1 + eps)/M) = 49,907,319,649,092,925,469,941,423,940,961 planned steps per
  walk (2^105.29902), so S_plan = K M = 2^125.29902.
* Tail limit 20/theta = 20 x 2^70 steps per walk (75 bits), so E = M x 20 x 2^70 = 2^94.32193.
* Widths (machine.md 4a). K needs 106 bits and the tail limit 75, so the planned counter and
  the stuck/tail counter are 128-bit (two limbs). The report count is capped at 2^55.44 and
  is a 64-bit word. The meter is a 192-bit instruction count compared against
  floor(20,000 x 2^124.65), which has 139 bits.

```
Phase 1 (planned search). The M walks advance round-robin. Each walk has a 128-bit
  planned counter (increment: adds/adc; compare with K: cmp/sbcs; machine.md 4a). A walk
  stops when its counter reaches K. Restarts (section 3.3) continue the same counter, so
  every walk takes exactly K steps unless the run has already ended.
  Each report is inserted into the multimap. If an entry with the same F from a different
  walk exists, its candidate x is computed and checked by g^x = z. A verified x ends the
  run with success (early exit, section 1).
Phase 2 (detection tail). Each walk that is not at a DP continues until its next DP, or
  for at most 20/theta steps. Reports are handled as in phase 1. Every collision of
  phase 1 has its merged trail at the shared next DP within 20/theta steps, except with
  probability e^-20.
Completion. If no verified x was found, the run returns FAILURE (normal no-hit completion).
```

Normal no-hit completion is the expected outcome of an unlucky run. It is not a resource
event: it consumes at most S_plan + E steps plus the non-walk bucket, and no cap or meter is
involved. The 128-bit counter and meter work is itemised in `machine.md` section 4a: 153
instructions in total, within the 300 charged per step, so the step prices are unchanged. The phase switch is one pass over
the M walk states, under 2^20 x 50 instructions, so it lies inside the non-walk bucket.

**3.5 Hard caps and premature truncation.**

*The meter.* The meter holds the exact number of charged instructions as a 192-bit integer
(three limbs). Every priced block has an integer instruction price, so no fractional cost is
lost:
- each walk step costs 18,191 or 6,965;
- each report costs 10,000;
- each chain-node visit costs 32;
- each start or restart costs 8,000,000.

Before each block, the meter checks acc + p_max <= CAP = floor(20,000 x 2^124.65) =
667450334702089094098810736037947573300028. Here p_max is the largest price the block can
have; for a walk step it is 18,191. If the check fails, the run aborts with FAILURE. The
block's own meter bookkeeping is part of its price: 18 instructions per step, inside
COMMON. So R <= CAP on every run, and W = R/20,000 <= 2^124.65. Part Q of
`planned_phase.py` executes these limb operations and checks them against integers:
- carries at 2^64 and 2^128;
- the counters at K-1, K, at the tail limit and at 2^128-2;
- a metered sequence started just below CAP, which stops before exceeding it and does not
  stop early.

*Premature truncation.* A complete run (phase 1 plus a full tail at the mean step price,
plus the non-walk bucket 2^64.5) costs 2^124.62992 nominally (section 5). The meter can fire
before completion only if the realised number m of multiply steps among the at most
n = S_plan + E steps satisfies m a + (n - m) b > CAP/20,000 - 2^64.5 - a, with a = 0.90955
and b = 0.34825. The extra -a covers the one-step pre-check. This needs
m/n > 1/2 + 0.015700 (part P, `tolerated_multiply_fraction_excess_points` = 1.5700). This is
the only event in which the meter ends a run that could still succeed. Its probability is
bounded as follows:

* The branch of a step is T mod 8 < 4. T is a class function, and its residue mod 8 has
  marginal 1/2 up to bias O(p^-1/2) (F2b; measured 0.5013 and 0.5004 on the production rule,
  part C). Before the first useful collision, a run visits distinct classes, apart from the
  escaped short cycles, which section 3.3 bounds. Under H1/H2, the visited classes behave as
  a uniform sample, and Hoeffding gives exp(-2 n 0.0157^2).
* **Robustness.** Suppose the branch bits are independent only between blocks of L
  consecutive steps, with arbitrary dependence inside a block. Hoeffding still gives at most
  2^-20 for every L <= 2^110.5.
* The failure budget charges **0.001** for premature truncation.

*Report store (machine.md 5).*
- Each report is a 96-byte node: a 16-byte fingerprint, two 32-byte exponents, an 8-byte
  next-node index and an 8-byte walk ID.
- Nodes are written sequentially into one array indexed by a 64-bit report counter.
- There are 2^57 bucket heads of 8 bytes each, indexed by the top 57 bits of F (load factor
  below 0.34 at the cap; the exact ratio is approximately 0.3383353650).
- The total is 2^62.34 bytes. That is below 2^63, so one 64-bit address space with ordinary
  pointers suffices. There is no sharding and no transport.
- Zeroing the heads (stp 16 bytes plus loop, 0.75 instruction per head) costs 2^42.3 units
  and is preprocessing.

Caps:
* **Reports.** D <= ceil(1.1 S_plan theta) + 2^25 = 2^55.437, against an expected
  n theta = 2^55.299. Under H2, Chernoff gives log2 P[D > cap] <= -3 x 10^14.
* **Chain-node visits.** At most 2^70 in total, at 32 instructions each (2^60.7 units). The
  expected number of same-bucket pairs is D^2/(2 x 2^57) = 2^52.87, so by Markov the cap is
  reached with probability at most 2^-17.1. That probability is charged.
* **False candidates.** At most 2^16 are verified; the expected number is 2^-18.1.
* **Restarts.** At most 2^48; the expected number is at most 2^20 + 2^26.4. Markov is charged.

Reaching any cap is a premature failure and is charged as such. It is never confused with
normal completion.

## 4. Success probability

**H1 (random mapping with branch mix).** The class sequence behaves like iterates of a random
mapping on N classes. Under this map, two distinct classes have the same image with
probability (1 - sum_b q_b^2) / N. For classes c1 != c2:

* a pair of different branches has about one solution among the N classes;
* squaring/squaring has none, since squaring is a bijection of a group of odd order that
  commutes with Frobenius;
* a pair with the same table index keeps 11/12 of one.

We use the conservative value sum q^2 = s^2 + (1-s)^2 / R = 0.250244, which gives
log2 N' = 250.42724. Under H1, the first S useful steps of the planned phase contain a class
collision between two different walks with probability 1 - exp(-S^2/(2N')) = 0.43. Wasted
cycle-escape steps are covered by the factor 1 + eps in S_plan. The predecessor's toy runs
(section 7) give a pooled P[T_c <= S(0.43)] = 0.432 (standard error 0.015).

The failure terms are (part P, `losses`):

| term | bound |
| --- | ---: |
| premature truncation by the meter (3.5) | 0.001 charged (computed far below) |
| escape overhead above 2 eps of the planned steps (H2, Chernoff), steps lost to restarts | 0.0005 charged |
| report cap (Chernoff under H2) | below 2^-1000 |
| chain-visit cap (Markov) | 7.0e-6 |
| restart cap (Markov) | 3.3e-7 |
| false-candidate cap (Markov) | 5.3e-11 |
| collision trail misses its DP within the tail | e^-20 |
| first collision inside one walk | 2^-20 |
| beta^1 = beta^2, or a key tie | below 2^-240 |
| reserve for the inaccuracy of H1/H2 at production size | 0.01 |

The success lower bound is 0.43 - 0.011508 = **0.41849 >= 0.39**. A further 0.0285 of
absolute slack remains for the model value itself. Normal no-hit completion (about 0.57) is
the complement of the collision event and is not a loss term.

## 5. Cost (W = R / 20,000)

All prices come from `machine.md`, which lists the executed AArch64 kernels and the counted
Fp12 routines:

| Step | Fp12 routine | keys | bookkeeping | Instructions | Units |
| --- | ---: | ---: | ---: | ---: | ---: |
| multiply | dense product 16,841 | 1,050 | 300 | 18,191 | 0.90955 |
| squaring | cyclotomic square 6,665 | 0 | 300 | 6,965 | 0.34825 |

Average step cost c = (0.90955 + 0.34825) / 2 = **0.6289 units**.

Non-walk work, each row a cap (`planned_phase.py` part P, `nonwalk_rows_log2_units`):

| row | log2 units |
| --- | ---: |
| setup: table and 2^20 starts, each exponentiation at most 400 units | 28.66 |
| zeroing 2^57 bucket heads (preprocessing) | 42.30 |
| reports: at most 2^55.437 at 10,000 instructions | 54.44 |
| chain-node visits: at most 2^70 at 32 instructions | 60.71 |
| restarts: at most 2^48 at 400 units | 56.64 |
| false candidates: at most 2^16 at 600 units | 25.23 |
| phase switch: 2^20 x 50 instructions | 11.36 |
| output: [x]P0 on G1 plus serialisation (section 1) | 6.66 (101 units) |

The total is 2^60.81 units. The ledger carries a non-walk bucket of **2^64.5**.
Preprocessing (setup plus zeroing) is below 2^43 units.

Detection tail: at most 20/theta steps per walk after the planned phase,
E = 2^20 x 20 x 2^70 = 2^94.32 steps (section 3.4). The planned phase is S_plan = K M = 2^125.29902
steps, which is S(1 + eps) rounded up to whole rounds.

`W = c S_plan + c E + 2^64.5 = 0.6289 x 2^125.29902 + 0.6289 x 2^94.32 + 2^64.5 = 2^124.62992`
(`planned_phase.py` part P, `log2_W_nominal_full_run`; `experiment.py` part D gives the same
2^124.6299 from the counted operation vectors, with its own E = 2^84.32; the difference
is below 10^-9 bits).

Submitted bound: **time_log2 = 124.65 = 2493/20**, enforced as CAP = floor(20,000 x 2^124.65)
instructions, 0.02008 bits above the nominal complete run.
The meter enforces it on every run (3.5). The margin is the truncation tolerance, not a
slack for unsuccessful runs. An unsuccessful run ends by normal completion at its realised
cost, which exceeds the nominal value only by the drift that 3.5 bounds.

Memory: the report store is 2^62.34 bytes (3.5). The 5.5 MB tables and 2^20 walk states
(under 2 KB each, including the 128-bit counters and history) add little, so the peak is
below 2^63 bytes in one 64-bit address space. Memory is reported, not scored.


## 6. Heuristics and limitations

* H1 (random-mapping behaviour of the mixed walk on classes), H2
  (pseudorandom branch, index and DP bits of T; independence of table indices for the cycle
  bounds). These are standard in rho analyses (Teske 1998/2001; van Oorschot-Wiener 1999;
  Duursma-Gaudry-Morain 1999; Bos-Kleinjung-Lenstra 2010) and are supported here at toy and
  production-field scale, not proved.
* The prices are analytic instruction counts: the leaf kernels are executed listings, the
  Fp12 routines are counted call by call, and the bookkeeping is an itemised estimate
  charged with margin (153 itemised including the 128-bit counters and 192-bit meter, 300
charged per step; 138 itemised, 160 charged for the
  key minimum). Each further 100 instructions per step would add 0.011 bits.
* Collision statistics are validated in toy fields only.
* This is a generic square-root attack. It does not use index-calculus or number-field-sieve
  methods. It is a valid upper bound under the stated heuristics, not a security estimate.

* The former H1+ is withdrawn. Its 0.01 bound on 'any cap or abort' conflicted with the 0.57
  no-hit probability. It is replaced by the planned phase (3.4) and the separately bounded
  premature-truncation event (3.5), which rely on H1/H2 and the stated dependence-block
  robustness.
* `experiment.py` is inherited byte-identical. Its scope is the leaf kernels, the BN254
  identities, the production step rule, and the toy runs of the inherited walk.
  - Part B's toy walk runs until a logarithm is found and resets a walk on restart. It does
    **not** execute this package's planned phase, tail, completion, 128-bit counters or
    meter.
  - It is evidence for the walk and the collision law (H1/H2) only.
  - Its part D field `success_if_squaring_fraction_drops_by_points` models the old stopping
    rule and is not used by this claim.
* `planned_phase.py` does not run the rho. Its part Q executes the control primitives in
  isolation, and its part P is an exact ledger.
* The planned-phase control flow as a whole has not been run at any scale.
  - It is specified in 3.4.
  - It differs from the toy walk only in the stopping rule, so the toy collision times
    remain applicable under H1.
  - The DP density changed from the toy settings and from the predecessor (2^-60 to 2^-70).
    Toy runs use their own DP bits.

## 7. Experiments (`experiment.py` inherited, `planned_phase.py` new; standard library only)

* **`planned_phase.py`** is the declared initial run: seed 20261112, cases 64, under 1 s.
  The recorded output is in `planned-phase-evidence.json`.
  - Part P: K = 49,907,319,649,092,925,469,941,423,940,961 (106 bits), log2 S_plan =
    125.29902, tail 75 bits, CAP 139 bits, log2 W_nominal = 124.62992, claim minus nominal
    0.020075 bits, tolerance 1.5700 points, block robustness 2^110.52, store 2^62.34 bytes,
    non-walk 2^60.81 units, preprocessing 2^42.30 units, success lower bound 0.41849. Every
    check in `checks` is true.
  - Part Q: 0 failures.
  - Part O: 73 of 73 scalars equal to the reference, [r]P0 = infinity, analytical bound
    1,005,133, sample maximum 998,551, charge 2,020,000.
  - `all_ok` is true, with exit 0.
* **`experiment.py`** is no longer an initial run, because the protected configuration
  allows one initial run. It remains declared and callable as a follow-up. Its recorded
  declared-run output is in `development-evidence.json` (inherited).

* **Part K:** the three leaf listings executed on an instruction simulator for p and for r
  (edge, random and all-limbs-large inputs; 0 failures).
* **Part A (BN254):** ring map from the profile basis, g of order r, Frobenius formula, key
  formula, distinct keys, canonical element under both the integer and the stored rule,
  p^6 inverts, trace invariance, Granger-Scott formula, counted operation vectors, orders of
  p and 2 modulo r, the r - 1 factorisation.
* **Part C (BN254, production step rule, R = 16 for speed):** 0 class-check failures in the
  first 150 steps; exponent invariant y = g^alpha z^beta at the end; squaring fraction 0.5013
  over 1,500 steps (declared run), 0.5004 over 4,000 steps (seed 7); chi-square of the low 3
  bits of stored T 1.6 and 16.1 (7 dof), of the next 4 bits 17.9 and 20.5 (15 dof).
* **Part B (toy fields with the same tower, xi = 9 + u, embedding degree 12):** the inherited
  run-until-found walk, not this package's planned-phase control flow (keys with tie-breaking by full conjugates, kappa trick, twisted table, squaring,
  2-cycle and checkpoint escapes, DPs, collision solving, verification). Toy keys take few
  values, so the branch and DP bits come from a hash of the canonical conjugate; part C tests
  the real bits. Declared run (q = 151, N = 710,195, R = 64, s = 1/2, 60 runs): all logs
  recovered and verified; mean T_c / sqrt N 1.447 against the model 1.451; fraction of runs
  with T_c <= S(0.43) 0.467.
* **Part E (pair-collision law):** samples distinct classes, applies one step and counts image
  coincidences; N times the pair-collision rate is the quantity H1 predicts.

Offline aggregate of the same script (all discrete logs recovered and verified; S(x) is the
step budget at model success x):

| q, N classes, R, seed, walks, dp_bits | s | runs | mean T_c / sqrt N (model) | P[T_c <= S(0.40)] | P[T_c <= S(0.42)] | P[T_c <= S(0.43)] |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| 151, 710,195, 64, seed 21, 4, 3 | 1/2 | 400 | 1.399 (1.451) | 0.417 | 0.435 | 0.443 |
| 283, 500,013, 64, seed 41, 8, 4 | 1/2 | 400 | 1.452 (1.451) | 0.383 | 0.400 | 0.417 |
| 163, 4,524,905, 256, seed 31, 8, 3 | 1/2 | 300 | 1.488 (1.448) | 0.410 | 0.420 | 0.437 |
| **pooled mixed walk** | 1/2 | 1,100 | 1.442 | 0.403 | 0.418 | 0.432 |
| 151, 710,195, 64, seed 12, 4, 3 (control) | 0 | 300 | 1.225 (1.263) | 0.420 | 0.443 | 0.460 |

The control row is a plain adding walk (model with the plain-walk effective class count). The
pooled fraction at S(0.43) is 0.432 with standard error 0.015; the 283 setting alone reads
0.417 (standard error 0.025). Measured 2-cycle rates with s = 1/2: 3.5e-4 and 2.6e-4 per step at
R = 64 (predicted (1-s)^2/(12R) = 3.3e-4), 7.9e-5 at R = 256 (predicted 8.1e-5); without
squaring 1.40e-3 (predicted 1.30e-3).

Pair-collision law (part E, q = 151, 60,000 sampled classes, relative standard error about 2.3%):

| R, s | observed | exact prediction | conservative value used | random mapping |
| --- | ---: | ---: | ---: | ---: |
| 64, 1/2 | 0.7311 | 0.7497 | 0.7461 | 1 |
| 8, 1/2 | 0.7126 | 0.7474 | 0.7188 | 1 |
| 8, 0 | 0.9663 | 0.9896 | 0.8750 | 1 |
| 16, 5/8 | 0.6285 | 0.6086 | 0.6006 | 1 |

The squaring branch removes its s^2 of the collision probability; the observed values lie
within about 2 standard errors of the exact prediction (three of four below it, one above).

Reproduce a row with e.g.
`echo '{"q":151,"trials":400,"R":64,"sq8":4,"dp_bits":3,"walks":4,"seed":21,"production_checks":0,"production_steps":0,"kernel_cases":10,"seconds":500}' | python3 -I experiment.py`.

References: Pollard 1978; van Oorschot and Wiener, J. Cryptology 12 (1999); Teske, Math.
Comp. 70 (2001); Duursma, Gaudry and Morain, ASIACRYPT 1999; Granger and Scott, PKC 2010;
Bos, Kleinjung and Lenstra, ANTS 2010; submissions 7da5f23a, 6d4e57d1 and
747ad6c9, 56fc259f; pairing package 3325b532 (source 6e57b2d5); GT package 2ece9f04.

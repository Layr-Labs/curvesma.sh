# Gated Frobenius-class rho on BN254 GT with generated lazy-reduction Fp12 kernels

**Claim.** For uniform z = g^x in the protected BN254 target group (g = e(P0, Q0) from
`target-profile.json`, x uniform in [0, r)), the algorithm below returns x with
probability at least 0.39 over uniform challenges and fresh coins, at **expected** total work at most

    E[W] <= 1 + (K / 2^64) (c E[T_c] (1 + eps) + c_mul E + 2^66) < 2^122.85 reference units

(W = R/20,000, R the total scalar AArch64 instruction-equivalent work). Step 0 is a run gate: draw a
fresh uniform 64-bit word v and stop as a failure (below one unit) unless v < K = 7416732153347139310, the
least K with (K/2^64) x 0.97 >= 0.39. Past the gate the algorithm succeeds with probability at least 0.97
(section 4). The score 122.85 is an expectation over the challenge and the coins of the one attack on the
one target, **not a bound on every run**: a run past the gate is stopped by two 192-bit step counters with
fixed caps, so every run costs at most 2^131.22 units (section 3.4). The work includes setup, table construction, walk starts, distinguished-point
(DP) handling, cycle escapes, restarts, the DP detection delay, false-candidate verification and output.
Peak memory is below 2^81 bytes (below 2^73 in expectation); preprocessing (included) is below 2^31 units;
nonuniform advice is 0 bytes. The bound rests on the standard random-mapping heuristic for rho walks,
adjusted for the branch mix, used here through its mean collision time (section 4: H1, H2 and H-mean),
supported by toy-field experiments and by running the production step rule in the real BN254 GT. This is an
analytic generic-attack bound, not a practical break, and no production-size run is claimed.

**What this package changes.** The walk, class construction, keys, table, cycle handling, gate and
expected-work argument are those of submissions 2ece9f04-8559-42ee-9f25-06bc993f1d89 (internettrashh,
accepted at 124.65; its `experiment.py` and `development-evidence.json` are included byte for byte) and
8b901b3d-a617-4b1f-a442-ff370e6ec4f7 (terrapinelf: the run gate and the run-to-collision expected-work
ledger). New here:
(i) **cheaper counted arithmetic**: the dense Fp12 product, the cyclotomic squaring and the key argmin are
generated straight-line scalar AArch64 programs with lazy reduction (one Montgomery reduction per output
coefficient, Aranha et al., eprint 2010/526), executed on an instruction simulator and checked against big
integers (`lazy.py` part L, `machine.md`): 8,785 + 898 and 4,323 instructions instead of 16,841 + 1,050 and
6,665, so the average step costs 0.36545 units instead of 0.6289;
(ii) **the caps are step counters, not an instruction meter**: the review of 8b901b3d found that its 128-bit
instruction accumulator could not hold its 2^145.37-instruction cap. Here no instruction meter exists;
two 192-bit step counters (3-word increment and compare, charged in every step) enforce fixed step caps
below 2^133, and the every-run work bound is a consequence of those caps and the per-step prices.
`lazy.py` part G recomputes every figure below from the executed lengths.

**Lineage of the base design.** The design is the Frobenius-class rho of the BLS12-381 fixed-argument track:
the class construction, 4-multiplication key and twisted multiplier table of submission
7da5f23a-6633-46b1-9c03-11f6cfedeeec (mpjunior92), the cyclotomic-squaring branch read from
the trace invariant of 6d4e57d1-135f-4d1f-9d43-eca27fb5c9ed and
747ad6c9-0453-4e37-8999-4f02a32c3c8c (jungjipdo), and the stored-representation rule
and charging at the counted listing of 56fc259f-d5d1-4c63-b4e5-9d95e6188bd6 (GordoAR, accepted at 125.72 on
bls12-381-fixed-g2-invert-g1). New here: (i) the port to BN254 (p, r of 254 bits, xi = 9 + u,
four-limb arithmetic) with every identity re-checked on the protected group; (ii) a complete
BN254 instruction ledger built from executed AArch64 listings (`machine.md`), since no BN254
operation has a fixed price on this track; (iii) squaring probability s = 1/2, chosen for
the BN254 price ratio; (iv) the reduction is direct: the track asks for x, not a preimage
point, so no pairing is computed.

## 1. Problem and output

Input z in GT. If z = 1 output 0. Otherwise find x with g^x = z, check it, output x. g has
prime order r (part A, `g_order_r`), so x is unique in [0, r).

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
Part A checks y^2 = formula on GT elements (`sq_ok`). The priced implementation is the lazy
routine of `machine.md` section 3 (18 products, 12 reductions), checked on GT elements in `lazy.py`.

**F4 Orders.** r - 1 = 2^28 * 3^2 * 13 * 29 * 983 * 11003 * 237073 * 405928799 *
1670836401704629 * 13818364434197438864469338081 (verified by multiplication). The order of
2 modulo r is about 2^247.31, so 2^t lies in the order-12 group <p> only if t is a multiple
of about 2^243.7.

## 3. Algorithm

Parameters: R = 1024 table multipliers, squaring probability s = 1/2, DP density
theta = 2^-60, M = 2^20 parallel walks. From stored T: squaring iff T mod 8 < 4; table index
h = (T >> 3) mod 1024; the class is distinguished iff (T >> 13) mod 2^60 = 0; the fingerprint
used for matching reports is F = (T >> 73) mod 2^128, disjoint from the bits above.

Setup: convert z to the tower (once). Draw fresh uniform a_h, b_h mod r and compute
M_h = g^{a_h} z^{b_h} and Tab[h][t] = M_h^{p^t} for h < 1024, t < 12, with the scalar tables
a_h p^t, b_h p^t mod r (12,288 entries of 384 + 64 bytes, 5.5 MB).

A walk keeps a member y of its current class and (alpha, beta) mod r with
y = g^alpha z^beta, starting from fresh uniform (alpha, beta). One step:

```
1  T <- stored Re(c_0) of y; branch, h and DP bits from T
2  if distinguished: compute the 12 keys -> kappa; report (F, p^kappa alpha, p^kappa beta)
   into a chained multimap keyed by F (entries are never replaced)
3  2-cycle and checkpoint tests on T (section 3.3); increment the 192-bit step counters, stop at a cap (3.4)
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
false match of probability 2^-128 per pair of distinct classes: the number D of reports is theta times the number of steps, and under H1 E[T_c^2] = 2 N',
so there are about theta^2 2 N' / 2^129 = 2^2.5 false candidates in
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
  steps). A walk that sees no DP within 20/theta = 2^64.3 steps restarts from a fresh pair.
* Overhead factor on the step count, charged: eps = 5/(12R) + 250/R^2 = 6.45e-4.

**3.4 Gate and caps.** Step 0 (gate): the fresh word v is independent of z and of every other
coin; if v >= K the run stops and fails. Past the gate the walks run until a verified logarithm is found or
a cap fires. Every cap is a failure, and its probability is charged in section 4:

* *Ordinary-step cap.* A 192-bit counter counts every ordinary step (a step that advances a walk's class
  sequence, including steps after the first useful collision). The run stops when it exceeds
  S_cap + E, with S_cap = 64 sqrt(N') = 2^131.2136 (N' in section 4) and E = 2^84.32 the detection delay.
  If the first useful collision occurs within S_cap ordinary steps it is detected within E further steps
  (or detection fails, charged below), so this cap fails a run only when T_c > S_cap.
* *All-step cap.* A second 192-bit counter also counts the extra steps of cycle escapes (re-traversal of a
  detected cycle, at most 200 steps per escape). The run stops when it exceeds 2 (S_cap + E). This fires
  only if escape steps exceed S_cap + E; their expectation is at most eps times the expected ordinary steps
  (section 3.3), so by Markov's inequality the probability is at most eps (E[T_c] + E) / (S_cap + E) =
  1.33e-5 (part G, `all_step_counter_cap_bound`).
* *Reports.* At most 2^8 times the expected number 2^65.61 of reports (Markov: probability at most 2^-8);
  a 128-bit counter suffices (cap 2^73.61).
* *False candidates.* At most 2^20 are verified (expected below 2^3, Markov below 2^-17; 2^29.3 units).
* *Restarts.* At most 2^52 (expected below 2^37, Markov below 2^-15).

The counters cost 7 instructions each per step (3-word ADDS/ADCS/ADC increment, 3-word SUBS/SBCS/SBCS
compare, branch) and are inside the charged per-step bookkeeping (`machine.md` section 4). Both caps are
below 2^133 < 2^192, so neither counter wraps. **Every-run bound:** a run executes at most 2 (S_cap + E)
steps, each at most c_mul = 0.49955 units, plus the capped non-walk work (at most 2^8 x 2^66 units), in all
at most 2^131.22 units. This bound is not the score; it only makes every run finite.

## 4. Success probability

**H1 (random mapping with branch mix).** The class sequence behaves like iterates of a random
mapping on N classes in which two distinct classes have the same image with probability
(1 - sum_b q_b^2) / N. For classes c1 != c2: a pair of different branches has about one
solution among the N classes; squaring/squaring has none (squaring is a bijection of a group
of odd order commuting with Frobenius); a pair with the same table index keeps 11/12 of one.
We use the conservative value sum q^2 = s^2 + (1-s)^2 / R = 0.250244 for all pairs (the
exact s^2 + (1-s)^2/(12R) is smaller). This is the mechanism by which Pollard's mixed walks
lose against r-adding walks (Teske 1998), here priced exactly. Effective class count
N' = N / (1 - sum q^2), log2 N' = 250.4272.

Under H1 the number T_c of steps (summed over all walks) before the first useful class collision satisfies
P(T_c > s) = exp(-s^2 / (2 N')), so E[T_c] = sqrt(pi N' / 2) = 2^125.5394.

**H-mean (the one quantitative premise added here).** E[T_c] <= 1.05 sqrt(pi N' / 2) = 2^125.6098. This is a
statement about the mean only; no tail of the collision law is used. Support: it is H1 with a 5% allowance;
the base package's toy runs of the complete algorithm measure this very mean: 1.399, 1.452 and 1.488 times
sqrt(N) in three toy settings against the model values 1.451, 1.451 and 1.448, pooled 1.442 over 1,100 runs
(section 7), i.e. 0.994 of the model with a standard error near 1.6%; the allowance is the pooled value plus
three standard errors. `experiment.py` part B reruns such settings (declared run: 1.447 against 1.451).
Limit: the mean is measured in toy fields only and extrapolated to N' = 2^250.43 as H1 is.

**Success past the gate.** The run fails only if
* T_c > S_cap: by Markov's inequality and H-mean, probability at most E[T_c] / S_cap = 1.05 sqrt(pi/2) / 64
  = 0.020563;
* the report cap fires: at most 2^-8 = 0.003907 (Markov on the number of reports, whose expectation is
  theta times the expected number of steps under H2);
* the all-step cap fires: at most 1.33e-5 (section 3.4);
* any of: detection failure e^-20; a first collision inside a single walk (about 1/M = 2^-20); a key tie
  (at most 66/p per element, F2); beta^1 = beta^2 (2^-253); the false-candidate or restart cap (below 2^-15):
  together below 0.001.
The sum is below 0.02549, so success past the gate is at least 0.9745; the claim uses 0.97. The two step caps are
the only ways a run without a collision ends, and their probabilities are listed above: there is no
separate abort or meter budget and no assumption on the squaring fraction. With the gate, success is at least (K / 2^64) x 0.97 >= 0.39, and K is the least such integer
(`lazy.py` part G: `success_ok`, `K_is_least`). The identity z = 1 is answered directly.

**Expected walk work.** Let c_i be the charge of step i (0.49955 or 0.23135 units). The event {T_c >= i} is
determined by the classes visited before step i, and the branch of step i is read from T of a class not
visited before. Under H2 that bit is unbiased given the past (exactly unbiased over all classes, F2b), so
E[c_i ; T_c >= i] = c P(T_c >= i) with c = 0.36545, and summing over i (Wald's identity)
E[walk work] = c E[T_c], times the cycle-escape overhead 1 + eps of section 3.3. An expectation needs only this
first moment; no concentration premise on the branch mix is used.

## 5. Cost (W = R / 20,000)

All prices come from `machine.md`; the routine lengths are the executed counts of `lazy.py` part L:

| Step | Fp12 routine | keys + argmin | call sites | bookkeeping | Instructions | Units |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| multiply | lazy dense product 8,785 | 898 | 8 | 300 | 9,991 | 0.49955 |
| squaring | lazy cyclotomic square 4,323 | 0 | 4 | 300 | 4,627 | 0.23135 |

Average step cost c = (0.49955 + 0.23135) / 2 = **0.36545 units**.

Non-walk work, in expectation: reports are a fraction theta = 2^-60 of the steps (H2), so their expected
number is below 2^65.62 and they cost 10,000 instructions (0.5 unit) each: 2^64.62 units; setup (table and
2^20 starts, each exponentiation under 400 units) under 2^31; false-candidate verification under 2^29.3 even
at its cap; restarts under 2^48 x 400 = 2^56.6 in expectation. All
non-walk work is charged **2^66** units.

Detection delay: after the first useful collision all M walks continue for at most 20/theta steps,
E = 2^20 x 20 x 2^60 = 2^84.32 steps, charged at the dearest price c_mul = 0.49955.

Expected work of a run past the gate (`lazy.py` part G):

`c E[T_c] (1 + eps) + c_mul E + 2^66 <= 0.36545 x 2^125.6098 x 1.000645 + 0.49955 x 2^84.32 + 2^66 = 2^124.1584`.

With the gate, W <= 1 + 1[v < K] W_run with v independent of W_run, so

`E[W] <= 1 + (K / 2^64) x 2^124.1584 = 2^122.8439`.

Submitted bound: **time_log2 = 122.85** (rational 2457/20), 0.0061 bits above the computed value. With the
base package's prices (c = 0.6289) the same ledger gives 2^123.627 (submission 8b901b3d); the cheaper
counted arithmetic alone lowers it by log2(0.6289 / 0.36545) = 0.783 bits.

Memory: at most 2^73.61 reports (the cap) x 80 bytes (128-bit fingerprint, two 256-bit exponents) at
load factor 2/3 = 2^80.52 bytes (2^72.52 at the expected number of reports), plus the 5.5 MB tables and 2^20 walk states (under 2 KB
each with the remembered history): below 2^81 bytes. Memory is reported, not scored.

## 6. Heuristics and limitations

* H1 (random-mapping behaviour of the mixed walk on classes), used only through H-mean (the mean
  collision time, section 4), and H2
  (pseudorandom branch, index and DP bits of T; independence of table indices for the cycle
  bounds). These are standard in rho analyses (Teske 1998/2001; van Oorschot-Wiener 1999;
  Duursma-Gaudry-Morain 1999; Bos-Kleinjung-Lenstra 2010) and are supported here at toy and
  production-field scale, not proved.
* The Fp12 routines and the key argmin are executed straight-line listings (exact counts); the
  per-step bookkeeping is an itemised estimate charged with margin (137 itemised, 300 charged).
  Each further 100 instructions per step would add 0.020 bits.
* The kernels are simulated, not run on hardware; the counts are instruction counts, as the cost model
  requires, not cycles.
* Collision statistics, including the mean collision time on which the score rests, are validated in toy
  fields only.
* The score is an expectation behind a run gate. A single run past the gate can cost up to its
  every-run bound, 2^131.22 units; the probability of reaching a step cap is bounded by Markov's inequality, not measured.
* This is a generic square-root attack. It does not use index-calculus or number-field-sieve
  methods. It is a valid upper bound under the stated heuristics, not a security estimate.

## 7. Experiments (standard library only)

`lazy.py` (new; declared experiment `lazy-kernels-and-ledger`): part L generates the three routines,
asserts the static bounds, executes them on the simulator and compares with big integers (declared run:
30 dense products on random, edge-valued and all-limbs-large inputs, 4 squarings on GT elements, 8 key
argmins against the twelve explicit conjugates; 0 failures, all outputs canonical; executed counts 8,785,
4,323 and 898); part G recomputes N, N', the charged mean, the caps and their Markov bounds, the failure
budget, K and its minimality, the expected run and total work and the declared score, in exact decimal
arithmetic; its inputs vary the mean allowance and the cap multiple. `experiment.py` (the base package's,
unchanged; declared experiment `frobenius-rho-bn254-checks`) supplies the identities, the production walk
and the toy runs; its part D prints the base package's ledger (claim 124.65) and is superseded by part G:

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
* **Part B (toy fields with the same tower, xi = 9 + u, embedding degree 12):** the complete
  algorithm (keys with tie-breaking by full conjugates, kappa trick, twisted table, squaring,
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
Bos, Kleinjung and Lenstra, ANTS 2010; Aranha, Karabina, Longa, Gebotys and Lopez, Faster explicit formulas
for computing pairings over ordinary curves, EUROCRYPT 2011 (eprint 2010/526); submissions 2ece9f04
(internettrashh), 8b901b3d (terrapinelf), 7da5f23a (mpjunior92), 6d4e57d1 and
747ad6c9 (jungjipdo), 56fc259f (GordoAR).

# Tag-traced Frobenius-class rho on BN254 GT

**Claim.** For uniform z = g^x in the protected BN254 target group (g = e(P0, Q0) from
`target-profile.json`, x uniform in [0, r)), the algorithm below returns x with
probability at least 0.39 using at most **2^123.02** reference units (W = R/20,000, R the
total scalar AArch64 instruction-equivalent work), including setup, the precomputed
multiset product table (2^47.35 units, built for this target), walk starts,
distinguished-point (DP) handling, cycle escapes, restarts, the DP detection delay,
false-candidate verification and output. A metered abort at 2^123.02 units makes the work
bound hold on every run. Peak memory is below 2^63 bytes; nonuniform advice is 0 bytes.
The success bound rests on the standard random-mapping heuristic for adding walks
(section 4) and on the pseudorandomness of the walk bits (H2), supported by toy-field
runs of the complete algorithm and by running the production walk in the real BN254 GT.
This is an analytic generic-attack bound, not a practical break; no production-size run is
claimed.

**Lineage and what is new.** The classes, keys, twisted table, stored-limb rule, fruitless
cycle escapes, caps, work meter and abort, and the arithmetic kernels (base-field leaves,
fused Fp2 leaves, counted Fp12 routines) are those of the accepted BN254 package
7e58d165-ebf6-41aa-8003-2f0a42370b6a (e-man07, 124.23), itself built on the accepted
2ece9f04-8559-42ee-9f25-06bc993f1d89 (internettrashh, 124.65) and its BLS12-381 lineage
(7da5f23a mpjunior92; 6d4e57d1, 747ad6c9 jungjipdo; 56fc259f GordoAR). New here:

* (i) **Tag tracing** (Cheon, Hong and Kim, "Speeding up the Pollard rho method on prime
  fields", ASIACRYPT 2008, and "Accelerating Pollard's rho on finite fields", J. Cryptology
  25 (2012)). A multiply step of the base walk needs only the class keys of the current
  element, and the keys depend only on its w^1 coefficient c_1 (F2). The walk keeps a
  materialised chunk base x and the multiset S of table items chosen since; because GT is
  abelian, the current element is x * P_S with P_S the product of the items of S, which is
  precomputed. Each step computes only c_1(x * P_S), one fused 6-term Fp2 inner product
  (`c1dot`, 1,662 executed instructions) instead of a dense Fp12 product (11,803); the dense
  product is paid once per chunk of 7 steps.
* (ii) **A plain 32-adding walk** (no squaring branch, so every step is a multiply step whose
  product can be deferred) with the walk bits read from the **least stored key** K, a class
  invariant available from c_1 (F2b).
* (iii) A **multiset product table** of all multisets of at most 7 of the 384 items
  (2^47.90 entries, 2^58.16 bytes, 2^47.35 units to build).
* (iv) An exact **useless-cycle law** for this walk (part F) that prices its escapes.

Result: 4,509.43 instructions (0.225471 unit) per step, against 0.470975 unit for the
base package's mixed step, and a smaller effective class count (N' = N * 32/31 instead of
N * 4/3); nominal W = 2^123.0036, claim 123.02.

## 1. Problem and output

Input z in GT. If z = 1 output 0. Otherwise find x with g^x = z, check it, output x. g has
prime order r (part A, `g_order_r`), so x is unique in [0, r).

## 2. Arithmetic facts (checked in `experiment.py`)

Tower (part A): Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi), xi = 9 + u; with
u = w^6 - 9 this is the profile field Fp[w]/(w^12 - 18 w^6 + 82). p = 7 mod 12; xi is
neither a square nor a cube in Fp2.

**F1 Frobenius classes.** (sum c_m w^m)^p = sum conj(c_m) d_m w^m with
d_m = xi^{m(p-1)/6}. p has order 12 modulo r, p^6 = -1 and p^4 - p^2 + 1 = 0 (mod r): every
y != 1 lies in a class of 12 elements, N = (r - 1)/12 classes, log2 N = 250.0117.

**F2 Keys from c_1.** zeta = xi^{(p^2-1)/6} lies in Fp and is a primitive 6th root of unity.
With c_1 = a + b u and B = a Re(d_1) + b Im(d_1),
key_k := Re(w^1 coefficient of y^{p^k}) = a zeta^j (k = 2j), B zeta^j (k = 2j + 1):
all twelve keys cost 4 Fp multiplications, 1 addition and 8 subtractions (part A
`keys_ok`, `distinct_ok`). They need **only c_1**. kappa(y) = argmin_k of the stored keys,
and y^{p^kappa} is the same element for the twelve members of the class
(`canon_ok`, `stored_canon_ok`; a tie of two keys has probability at most 66/p per element
and is charged as failure).

**F2a Representation.** Fp elements are Montgomery residues x 2^256 mod p (machine.md
section 1); the rule reads stored limbs; Montgomery products of stored operands are stored
products, so every structural property of the integer rule carries over (as in the base).

**F2b The least key is a class invariant.** The keys of y^{p^k} are the keys of y cyclically
shifted (F2), so K([y]) := min_k stored key_k is the same for all twelve conjugates (part A
`minkey_inv_ok`). The walk reads its bits from K.

**F3 Cyclotomic squaring** (Granger-Scott; part A `sq_ok`), used only in escapes.

**F4 Orders.** r - 1 = 2^28 3^2 13 29 983 11003 237073 405928799 1670836401704629
13818364434197438864469338081; the order of 2 modulo r is about 2^247.31, so 2^t is in the
order-12 group <p> only for t a multiple of about 2^243.7 (escape squarings never undo
each other).

**F5 Tag identity.** For x, P in Fp12, c_1(x P) = sum_{i<6} x_i Qt_i with
Qt = (P_1, P_0, xi P_5, xi P_4, xi P_3, xi P_2) (from w^6 = xi). Part A checks it on random
GT elements (`tag_c1_ok`); part C checks it at every step of the production walk; part K
checks that the `c1dot` listing computes (Re, Im) of this sum times 2^-256 mod p.

## 3. Algorithm

Parameters: R = 32 multipliers, chunk length l = 7, DP density theta = 2^-70, M = 2^20
parallel walks, checkpoint period 64. From the stored K of the current class: table index
h = (K >> 3) mod 32; the class is distinguished iff (K >> 13) mod 2^70 = 0; the fingerprint
used for matching reports is F = (K >> 83) mod 2^128, disjoint from the bits above.

**Setup.** Convert z to the tower (once). Draw fresh uniform a_h, b_h mod r, compute
M_h = g^{a_h} z^{b_h} and the 384 items Tab[h][t] = M_h^{p^t} (h < 32, t < 12) with
exponents (a_h p^t, b_h p^t) mod r. Build the **multiset product table**: for every multiset
S of 1 to 7 items, P_S = product of its items, its 12 stored Karatsuba sums, its exponent
pair (A(S), B(S)) and, for |S| <= 6, its tag block (machine.md section 5). Entries are
addressed by the rank of the sorted multiset (combinatorial number system).

**Walk state.** A chunk base x (materialised, with its six stored sums), the base exponents
(alpha, beta) with x = g^alpha z^beta, the sorted multiset S (|S| <= 6), and the history
needed for the escapes. Fresh walks start from fresh uniform (alpha, beta), S empty.

**One step** (the class map is that of the base package's multiply step):

```
1  c <- c_1 of the current element y = x * P_S:  c = stored c_1 of x if S is empty,
       else c1dot(x, tag block of P_S)                                  [F5]
2  12 keys from c -> K, kappa                                          [F2, F2b]
3  if distinguished: report (F, p^kappa (alpha + A(S)), p^kappa (beta + B(S)))
   into a chained multimap keyed by F (entries are never replaced)
4  2-cycle and checkpoint tests on K (section 3.4); charge the step to the work meter
5  t <- (12 - kappa) mod 12;  S <- S + {(h, t)}
6  if |S| = 7:  x <- x * P_S (one dense product), (alpha, beta) += (A(S), B(S)), S <- {}
```

**3.1 The step is a map on classes.** Write y^ = y^{p^kappa}. The step replaces y by
y * Tab[h][t] = y * M_h^{p^{-kappa}} = (y^ M_h)^{p^{-kappa}}, which lies in the class
[y^ M_h]; K, h and kappa are functions of the class (F2b) and of the member held, and the
member held enters only through the conjugation absorbed by the twist. So the class
sequence is the iteration of f([y]) = [y^ M_{h(K([y]))}], an R-adding walk on classes, and
walks that reach the same class stay merged and reach the same DP.

**3.2 Tag tracing computes exactly this walk.** Let the direct walk be
y_{j+1} = y_j * Tab[h_j][t_j] with (h_j, t_j) read from y_j. By induction on j within a
chunk, y_j = x * P_{S_j} with S_j the multiset of the first j items (GT is abelian, so the
order of the factors does not matter): the step reads c_1(y_j) = c_1(x * P_{S_j}) (F5), hence
the same keys, K, kappa, h and t as the direct walk; at the chunk end the new base is
x * P_{S_7} = y_7, and the exponents are those of the direct walk because
(A(S), B(S)) is the sum of the items' exponents. A DP report uses the exponents of y_j,
(alpha + A(S_j), beta + B(S_j)). An escape materialises the element it squares as
x * P_S from the stored entry. So every step, report, escape and output is the same as in
the direct walk; only the arithmetic that produces it differs. Part C runs both walks in
lock step in the real BN254 GT and compares, at every step, c_1, K and the decision, and at
every chunk end and escape the materialised element and its exponents.

**3.3 Collision.** Two reports with equal F from different walks are the same class up to a
false match of probability 2^-128 per pair of distinct classes; with at most
D = 2^55.25 reports there are about D^2/2^129 = 2^-18.5 false candidates in expectation, each
checked (one inverse mod r and g^x = z, under 600 units). For a true match the canonical
exponents satisfy alpha^1 + x beta^1 = alpha^2 + x beta^2 (mod r), so
x = (alpha^1 - alpha^2)/(beta^2 - beta^1); beta^1 = beta^2 has probability 1/r. Output x
after checking g^x = z.

**3.4 Fruitless cycles and their price.** After l steps the canonical exponent pair moves by
sum_h c_h (a_h, b_h), c_h = sum over the occurrences of index h of lambda^{t} (lambda = p mod
r). Since the start pair and the table exponents are uniform and independent, a closed class
walk returns the pair unchanged (a useless cycle) only if every c_h = 0 mod r, and otherwise
it is an ordinary collision (probability of order 1/N per pair). For one index used m times,
c_h is the image of the sum of m twelfth roots of unity under Z[zeta_12] -> Z/r; a nonzero
element of Z[zeta_12] with m < r^{1/4} terms has norm below r and is not in the kernel, so
c_h = 0 iff the root sum vanishes. Under H2 (index uniform on 32, twist uniform on 12,
independent) its probability is q(m), computed exactly for m <= 16 in part F
(q(1) = 0, q(2) = 1/12, q(3) = 1/72, q(4) = 11/576, ...; m^2 q(m) <= 1/3 for every m <= 16,
and offline for every m <= 30, with m^2 q(m) rising to 0.2941 at m = 30 towards the local
limit 12/(4 pi^2) = 0.304: Z[zeta_12] has covolume 3 in C^2 and a uniform root has covariance
I/2. Above m = 16 part F uses the bound 1/(3 m^2), which this limit supports.) The
probability that l given consecutive steps close a useless cycle is

    G(l) = l!/R^l [x^l] (sum_m q(m) x^m/m!)^R:
    G(2) = 1/(12R) = 2^-8.585, G(3) = 2^-16.17, G(4) = 2^-15.59, G(8) = 2^-27.61,
    G(16) = 2^-47.37, G(64) = 2^-122.33, sum_{3<=l<=64} G(l) = 3.45e-5

(part F, in log space: G is nonincreasing from l = 6 to 300, G(300) = 2^-260.7, and
sum_{65<=l<=300} G(l) = 2^-122.62; offline the same holds to l = 1500, G(1500) = 2^-405,
`development-evidence.json`). Beyond any finite range the tail is bounded the same way: a
useless cycle of length l > 300 needs every one of the at least 25 indices used (all but
probability C(32,8)(3/4)^300 < 2^-101 for l > 300, by a union bound on 8 missing indices) to cancel, each with
probability at most 1/12 (q(m) <= 1/3 m^-2 <= 1/12 for m >= 2, and q(1) = 0), which is at
most 2^-89 per start position. Such a cycle is abandoned by the stuck-walk restart or, if
it contains a distinguished class, by the rule that a walk whose report repeats its own
previous report (same F and canonical pair) restarts; it wastes at most 2^74.3 steps, and
counted over its l start positions this is at most sum_{l > 300} 2^-89 x 2^74.3 / l < 2^-9
step per step, inside the eps margin of 0.0094 (section 3.5).

* **Detection and escape.** If the K of two steps back equals the current K, the two
  classes form a 2-cycle: the walk materialises the representative with the smaller K and
  squares it (exponents doubled); this depends on the cycle only, so merged walks stay
  merged. A checkpoint K is stored every 64 steps and compared at every step; a repeat means
  a cycle of length l <= 64, which is traversed once more to its least-K class, materialised
  and squared. A walk that sees no DP within 20/theta = 2^74.3 steps (a two-word counter),
  or whose report repeats its own previous report, restarts from a fresh
  pair.
* **Price.** An escape costs 17,606 instructions = 3.90 steps (charged as 4) (machine.md section 4). A
  2-cycle wastes its 2 steps plus the escape; a cycle of length 3..64 at most 63 + l steps
  of checkpoint delay, l steps of traversal and the escape (67 + 2l steps). The expected
  overhead per step is G(2) x 6 + sum_{l=3}^{64} G(l)(67 + 2l) = 0.01563 + 0.00256 =
  0.0182 (part F). Longer useless cycles occur with probability below 2^-122.6 per step
  (part F: sum over 65 <= l <= 300, G nonincreasing from l = 6) and cost at most
  2^74.3 + 2l steps: below 2^-48 steps per step. **Charged:**
  eps = 8/(12R) + 1000/(12R)^2 = 0.02762 on the step count (part F `eps_covers`).
* Measured: 2-cycle rates 0.002515, 0.002609, 0.002575 per step in toy runs at R = 32
  (predicted 0.002604) in the earlier set and 0.002495, 0.002544, 0.002572 in the final
  set; 172 two-cycle escapes in 60,000 steps (2.9e-3) of the production walk at R = 32
  (predicted 1/384 = 2.6e-3, +1.3 standard deviations); checkpoint (3..64) cycles 14, 21 and
  22 events in the three final toy settings, about 2.7e-5 to 5.1e-5 per step (predicted
  3.45e-5). The charge eps covers twice the predicted overhead (section 3.5).

**3.5 Caps and the total-work abort.** The algorithm meters its own charged work with the
prices of section 5 and **aborts with failure as soon as the metered work exceeds
W_cap = 2^123.02**. The bound holds on every run by construction; the probability that the
abort or a cap fires is part of the failure budget.

*Caps.* Reports D <= 1.1 S theta + 2^25 = 2^55.25 (S theta = 2^55.11). At most 2^16 false
candidates are verified (expected 2^-18.5; at most 2^25.2 units). At most
2^48 restarts (Markov; at most 2^56.6 units). Steps after the first useful collision are at
most E = 2^94.32.

*No drift (H1+).* The price of a step is fixed (one chunk of 7 steps is 31,566 instructions
whatever the walk does); the only variable charge is the escapes, expected 0.0182 step
per step against the 0.0276 charged. The cap exceeds the nominal ledger by 0.0164 bits
(1.1 %, 0.00258 unit per step): it is reached only if the realised escape overhead
exceeds its expectation by a factor of about 2.1 over the whole run (0.0094 step of eps
margin plus 0.0115 step of cap slack per step, against an expected 0.0182). Escape events are
rare, local (they depend on two to sixty-four consecutive steps) and number about 2^116;
their rate is set by the step bits, whose pseudorandomness is H2. H1+ states that the
realised escape overhead, the caps and the abort together fire with probability at most
0.01. Support: the exact law of part F; the toy and production rates above; the bits of K
in the production walk pass chi-square tests (part C).

## 4. Success probability

**H1 (random mapping for an adding walk).** The class sequence behaves like iterates of a
random mapping on N classes in which two distinct classes have the same image with
probability (1 - sum q^2)/N, q the index probabilities. For classes c1 != c2 with different
indices there is about one solution among the N classes; with the same index the images
coincide only through the conjugation (11/12 of one). We use the conservative value
sum q^2 = 1/R = 1/32 for all pairs (the exact 1/(12R) is smaller); this is the adding-walk
loss of Teske (Math. Comp. 2001) and of the base package's control row. Effective class
count N' = N/(1 - 1/32), log2 N' = 250.0575.

Under H1 the first S steps contain a class collision with probability 1 - exp(-S^2/(2N')).
With P_model = 0.43: `S = sqrt(2 N' ln(1/0.57)) = 2^125.1132`.

Losses (as in the base package): detection failure e^-20; a collision inside a single walk
(fraction about 1/M); beta^1 = beta^2 (2^-253); caps and abort at most 0.01 (H1+); cap
failures by Markov below 2^-11; key ties below 2^-247 per step; steps lost to restarts. Their sum
is below 0.011, so the model value 0.43 leaves about 0.03 absolute slack over the required
0.39 for the inaccuracy of H1.

**Toy evidence for this walk** (part B: the complete algorithm with R = 32, no squaring
branch, keys with tie-breaking by full conjugates, twisted table, 2-cycle and checkpoint
escapes, DPs, collision solving and verification; toy keys take few values, so the toy reads
its bits from a hash of the canonical conjugate, and part C tests the real bits):

| q, N classes, seed, walks, dp_bits | runs | logs recovered and verified | mean T_c/sqrt N (model 1.2734) | P[T_c <= S(0.40)] | P[T_c <= S(0.42)] | P[T_c <= S(0.43)] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 151, 710,195, 21, 4, 3 | 400 | 400 | 1.2320 | 0.425 | 0.453 | 0.460 |
| 283, 500,013, 41, 8, 4 | 400 | 400 | 1.2612 | 0.370 | 0.395 | 0.400 |
| 163, 4,524,905, 31, 8, 3 | 300 | 300 | 1.2454 | 0.407 | 0.427 | 0.440 |
| **pooled (final code)** | 1,100 | 1,100 | 1.246 | 0.400 | 0.425 | 0.433 |
| earlier development set (same algorithm and settings, other random streams) | 1,100 | 1,100 | 1.229 | 0.421 | 0.441 | 0.450 |

Over the final-code runs the pooled fraction at S(0.43) is 0.433 (standard error 0.015) and
every mean lies below the model mean 1.2734; the q = 283 setting alone reads 0.400 (standard
error 0.025). Over all 2,200 runs the fraction is 0.441. At toy scale the model with
sum q^2 = 1/R is therefore roughly calibrated in its quantiles (not demonstrably
conservative), while the means sit below it; the 0.03 slack between the model value 0.43
(less 0.011 of losses) and the required 0.39 is about two pooled standard errors. The base package's control row (plain adding walk,
R = 64) read 1.225 against its model 1.263.

## 5. Cost (W = R/20,000)

From `machine.md` (executed listings, counted routines, itemised bookkeeping):

| Per chunk of 7 steps | Instructions |
| --- | ---: |
| 7 key computations (1,050) and 7 x COMMON (300) | 9,450 |
| 6 tag kernels c1dot (1,666 with call site) | 9,996 |
| dense product x * P_S (11,803), exact base sums (`add4sums`, 65), exponents (52), bookkeeping (200) | 12,120 |
| **total** | **31,566** |

c = 31,566/(7 x 20,000) = **0.225471 units per step**.

Non-walk work: at most 2^55.25 reports at 0.5 unit, 2^54.25 units; the table, 2^47.35
units; setup under 2^31; false candidates under 2^25.2; restarts under 2^56.6. All non-walk
work is below **2^57** units.

Detection delay: E = 2^20 x 20 x 2^70 = 2^94.32 steps.

`W = c S (1 + eps) + c E + 2^57 = 0.225471 x 2^125.1132 x 1.02762 + 0.225471 x 2^94.32 + 2^57 = 2^123.0036`.

Submitted bound: **time_log2 = 123.02** (0.0164 bits above the computed value; the metered
abort sits at this value). Part D recomputes the ledger in exact decimal arithmetic from the
counted operation vectors and the listing lengths.

Memory (machine.md section 5): table 2^58.16 bytes, DP store 2^62.2 bytes, walk states under
2^31 bytes: below 2^63 bytes, all in one flat 64-bit address space shared by the walk
processors (every pointer is an ordinary 64-bit address). Memory is reported, not scored.

## 6. Heuristics and limitations

* H1 (random-mapping behaviour of the 32-adding walk on classes, with the conservative
  sum q^2 = 1/R), H1+ (no drift: escape overhead, caps and abort within 0.01), H2
  (pseudorandom index, twist and DP bits of K; independence of indices and twists for the
  cycle law). These are the standard rho heuristics (Teske 1998/2001; van Oorschot-Wiener
  1999; Duursma-Gaudry-Morain 1999; Bos-Kleinjung-Lenstra 2010), supported here at toy and
  production-field scale, not proved.
* Tag tracing changes no heuristic: it computes the same walk (section 3.2), which part C
  checks step by step in the real group. It changes the work per step and adds a table
  whose size and build cost are charged.
* The squaring branch of the base package is dropped: the walk is a plain adding walk
  (better collision constant, sum q^2 = 1/32 against 1/4 + 1/2^12), at the price of more
  fruitless 2-cycles, which are charged through the exact law of part F.
* Lazy reduction is confined to the inside of the fused listings, whose range conditions are
  proved in machine.md (every value is below 2^512 before its reduction and every output
  is canonical) and exercised on edge inputs.
* Prices are analytic instruction counts: kernels are executed listings, the Fp12 routines
  are counted call by call, the bookkeeping is itemised and charged with margin (125
  itemised, 300 charged per step; 24 itemised, 200 charged per chunk). Each further 100
  instructions per step add 0.032 bits.
* Collision statistics are validated in toy fields only. The memory of the table is large
  (2^58.16 bytes) but memory is not scored; the instruction metric does not price latency.
* This is a generic square-root attack. It does not use index calculus or the number field
  sieve. It is an upper bound under the stated heuristics, not a security estimate.

## 7. Experiments (`experiment.py`, standard library only)

* **Part K:** the six listings of the base package (base-field leaves for p and r, fused Fp2
  leaves) and the new `add4sums` and `c1dot` listings, executed instruction by instruction on the simulator
  against Python integers (231 c1dot inputs per declared run including all-zero, all-(p-1)
  and edge limb patterns; 0 failures; stack pointer and callee-saved registers restored;
  only the output and the frame written).
* **Part A (BN254):** ring map from the profile basis, g of order r, Frobenius formula, key
  formula, distinct keys, canonical element (integer and stored rule), p^6 inverts,
  Granger-Scott squaring, counted products, **least-key invariance** and the **tag identity
  F5** on random GT elements, orders of p and 2 modulo r, the r - 1 factorisation.
* **Part C (BN254, production walk, R = 32, chunk 7):** the tag-traced walk and the direct
  walk in lock step; the base blocks given to the c1dot listing are built by the
  `add4sums` listing. Declared run (seed 20261112, 3,000 steps): 0 mismatches of c_1, K or
  decision; 432 materialisations (chunk ends and escapes), 0 mismatches; 9 two-cycle
  escapes; 24 c1dot listing runs on real stored values, 0 failures; exponent invariant
  y = g^alpha z^beta at the end; chi-square of the low 3 bits of stored K 11.36 (7 dof) and
  of the index bits 3..7 31.15 (31 dof). Offline (development evidence: a time knob outside
  the declared schema), seeds 7, 8, 9, 20,000 steps each: 0 mismatches in 60,000 steps,
  8,645 materialisations, 172 two-cycle and 1 checkpoint escapes, chi-square (7 / 31 dof)
  16.4/26.7, 9.3/21.2, 9.3/28.3.
* **Part B (toy fields, same tower):** the complete algorithm (section 4 table). Declared
  run (q = 151, R = 32): 60 runs, all logs recovered and verified, mean T_c/sqrt N 1.219
  (model 1.273), two-cycle rate 0.00289 per step (predicted 0.00260). The 1,100-run rows of
  section 4 use time budgets above the declared schema's 10 s (development evidence).
* **Part F:** q(m) exactly for m <= 16, the law G(l) in log space for l <= 300 (with its
  monotonicity and the tail sum over 65..300) and the escape-price check `eps_covers`.
* **Part D:** the ledger in exact decimal arithmetic, the table size and build cost.

Reproduce an offline row with e.g.
`echo '{"q":151,"trials":400,"R":32,"sq8":0,"dp_bits":3,"walks":4,"seed":21,"production_checks":1,"production_steps":1,"kernel_cases":10,"seconds":600}' | python3 -I experiment.py`
and the long production check with
`echo '{"seed":7,"production_steps":20000,"prod_seconds":500,"production_checks":1,"kernel_cases":10,"trials":1,"seconds":1}' | python3 -I experiment.py`.

References: Pollard 1978; Cheon, Hong and Kim, ASIACRYPT 2008 and J. Cryptology 25 (2012);
van Oorschot and Wiener, J. Cryptology 12 (1999); Teske, Math. Comp. 70 (2001); Duursma,
Gaudry and Morain, ASIACRYPT 1999; Granger and Scott, PKC 2010; Bos, Kleinjung and
Lenstra, ANTS 2010; submissions 7e58d165 (e-man07), 2ece9f04 (internettrashh), 7da5f23a
(mpjunior92), 6d4e57d1 and 747ad6c9 (jungjipdo), 56fc259f (GordoAR).

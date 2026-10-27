# Frobenius-class rho (2^20 round-robin walks) with a cyclotomic-squaring branch, priced by its measured instruction count

**Claim.** For uniform nonidentity Q in G2 and uniform z in GT the algorithm below
recovers P in G1 with e(P,Q) = z with probability at least 0.39, using at most
**2^125.234** dense-Fp12-multiplication equivalents (C = 20,000), including setup, the two
pairings, table construction, walk starts, distinguished-point (DP) handling, cycle
escapes, the DP detection delay, matching (the sort, duplicate removal and pair checks), output and verification. Peak memory is below 2^63
bytes, preprocessing (included) at most 2^31 units (the setup items of section 5), nonuniform advice 0 bytes. The success
bound rests on the usual random-mapping heuristic for rho walks, adjusted for the branch
mix (section 4), and is supported by toy-field experiments and by running the production
branch rule with a 16-entry table and no distinguished-point test in the real GT. The success probability is computed in the ideal model, in which the ChaCha20
stream that supplies the table and start scalars is indistinguishable from an ideal random source; that is
a disclosed heuristic, and no total-variation or independence bound on the realized scalars is claimed.
No production-size run is claimed.

The class construction (Frobenius classes, the 4-multiplication key, the twisted multiplier
table) is published in submission 7da5f23a-6633-46b1-9c03-11f6cfedeeec (mpjunior92); the
cyclotomic-squaring branch, the collision law, the caps and the bits read from T = Re(c_0) are
later revisions of that design. Section 0 lists the diff from that record: kernel prices, meter and counters, the append-only log and matcher, the restart cap and s = 37/64.

**What is new (parent gt-measured-sq-p0415, submission dc1293cb, 2^125.251).** Two changes: P_model
**0.415 -> 0.41** (section 4), shortening S by the factor 0.99171, and the multiply-step key block is
reworked (item 2 below): the nine key add/sub/neg operations are inlined at 30 each instead of 45, and the
minimum-of-12 is a specified limb-level procedure with worst-case charge 177, so the nominal ledger falls to
**2^125.2084098** and the bound to **125.234**; everything else is the p0415 parent's (dependencies.json).

## 0. What this revision changes

Below, 'record' = the historical `gt-tracecounted-exact-v25` (2^125.72, 56fc259f), not a qualified baseline; the list is the diff from it, changed quantities this package's.

The record charged the squaring map an analytic 11,369 instructions (section 5); this package charges a
straight-line kernel that computes the same map in 7,226 executed instructions, against blst's measured 9,803:

| Work | record's charge | this package |
| --- | ---: | ---: |
| every step: itemised bookkeeping (section 5) | 270 | **226** (scalar update and T read recounted, itemised listing) |
| every step: 192-bit work meter with its exact pre-debit compare before every step, and the 2-word stuck-walk counter (appendix.md section A-i) | not stated | **13** |
| every step, in all | 270 = 0.0135 units | **239 = 0.01195 units** |
| multiply step, keys only | 2,119 | **2,011** (key ops inlined, specified limb-level minimum) |
| multiply step in all | 1.11945 units | **17,742 = 0.8871 units** |
| cyclotomic squaring routine | 11,369 (analytic) | **7,226** (kernel, section 5) |
| squaring step in all | 0.58195 units | **7,465 = 0.37325 units** |
| squaring fraction s (branch rule) | 3/8 (T mod 8 < 3) | **37/64** (T mod 64 < 37, a 6-bit field) |

**Changes from the record** (record value -> this package, where the record states one):

1. *Squaring price.* The squaring routine is the straight-line kernel `kernels/sqr.s`, 7,226 executed
   instructions, not the analytic 11,369 and not blst's measured 9,803; the squaring step is 7,465 = 0.37325
   units, not 0.58195. The package ships the kernel generator, emulator and check (`kernels/`), and the
   blst measurement record (`measure/`) as the baseline.
2. *Multiply-step price.* 1.11945 -> **0.8871** units (17,742 instructions): the fixed-multiplier kernel's
   15,492 executed instructions replace the reference multiply's 20,000; the 2,011 of key work and the 239
   of every-step work are added (15,492 + 2,011 + 239 = 17,742). The key block is reworked: the nine
   add/sub/neg operations are inlined at 30 each (no call/return/address set-up), and the minimum of the
   12 keys is a specified limb-level procedure — seeded top-limb load plus 11 challengers, each 1 load +
   cmp + 3 csel on the top limb, with a full-width descent (b.eq, up to 5 lower-limb compares, deciding
   csel triple, tie flag) on top-limb ties, worst case 177 — so full-width ties resolve deterministically
   to the lowest index and are flagged; `experiment.py` implements this exact procedure as `Tower.min12`
   and part A asserts it agrees with the reference minimum with no ties on production GT elements.
3. *Every-step work.* 270 -> 239 instructions: the base bookkeeping is 226 — the (alpha, beta)
   update is recounted at 64 (two 4-limb add-and-conditional-subtract pairs, 32 each), not 100, and the
   T read at 42 (three ldps and the bit tests), not 50 —
   plus the 192-bit work meter (8) and the 2-word stuck-walk counter (5), both with stated widths
   and charged updates (appendix.md section A-i), plus 40 per walk start or restart.
4. *Average step cost, ledger and claim.* c 0.9178875 -> **0.59003046875** units; W_nominal 2^125.7075 ->
   **2^125.2084098**; claim **125.72 -> 125.234**. This cap is the least thousandth covering 2 points of
   squaring-fraction drift (drift cost 0.87089% per point, ledger with that drift 2^125.2333219), at
   least 1.79% (0.0256 bits) over the nominal ledger 2^125.2084098, with a walk-budget slack of
   more than 0.010560 units per step covering 2.0549 points; the drift table of section 3.4 is recomputed, and the
   cap is checked by an exact integer inequality.
5. *Distinguished-point rate and what follows from it.* theta 2^-60 -> **2^-72**
   (DP bits 16 to 87 of T, fingerprint bits 88 to 215; the record used 13 to 72 and 73 to 200). Here the
   report cap is D = 2^55 (a factor 1.94 above the expected reports, below 2^54.05), the stuck-walk
   threshold 20/theta is 20*2^72, the detection window E = M 20/theta is below 2^96.33 steps, the expected
   false-key pairs at most 2^-19, and the restart cap is **2^36** against expected restarts below 2^27.30
   (appendix.md section A-iii; the record's was 2^50), the fresh-start work 2^46.1 units. The record's ledger had
   c M 20/theta; this one has c (M + 2^36) 20/theta (below 2^112.33 steps), and the non-walk bucket is
   2^54.2 units (the record's is 2^64.3).
6. *Report store and matching.* The record kept a chained multimap at load factor 2/3 (peak memory below
   2^73 bytes), 6,000 instructions per report, and matched each report as it arrived. This package appends
   each report to a log of exactly 2^55 slots of 80 bytes with no pointer field, addressed by the one-word
   report counter in one 64-bit space (2^61.322 bytes, peak memory below 2^63 bytes), at 4,000 instructions
   per report itemised (section 5), and matches once, after the walk phase: an in-place heapsort of the
   whole record, adjacent dedup and every pair inside each run of equal fingerprints under a hard budget of
   2^16 pairs whose overflow fails the run (appendix.md section A-ii). The sort, dedup and scan read only written
   slots, so no initialisation is needed (the OS zero-fill of fresh pages is charged inside the 4,000,
   appendix.md section A-ii). Section 0-iv states every initial value, including a one-word restart counter charged
   below 2^23.30 units.
7. *Success chain.* The store loses nothing: every report is appended until the cap, and reaching the cap
   fails the run, a cap failure inside the 0.01 of H1+. The
   losses sum to below 0.0100496 (the restart cap's worst-case step fraction and the sampler's exhaustion allowance included), below 0.0101: success after every loss at least 0.3999, 0.0099 over 0.39. The
   model value is 0.41 (this package's change), the 0.01 charge for H1+, N and eps = 6.5e-4 are inherited; N' and S = 2^125.9685
   follow from s = 37/64 (difference 13).
7a. *Work-cap rule and walk phase.* The record's section 3.4 failed a run once its metered work was above
   the submitted bound, and it matched reports as they arrived. Here the walk phase runs until its budget
   is spent; before every charged action the meter is compared with a threshold that leaves room for that
   action's largest charge, the action runs only if it fits, and the meter then adds its charge. The
   budget leaves the setup bound and a reserved stopping allowance, so total work is at most W_cap on
   every run (sections 0-i and 3.4). The M walks run round-robin on one processor, K = 2^63 steps a
   chunk, switches below 2^58.65 units (appendix.md section A-i).
8. *Representation text.* Section 0a states the stored-limb rule directly, section 3 reads T as the
   stored limbs, and a new paragraph of section 0a argues from the blst source that the squaring and
   multiply branches return the same stored limbs (argued from the blst source; for the shipped kernels,
checked by `kernels/check_kernels.py` against canonical Montgomery limbs; section 6).
9. *Small remainder.* The model 2-cycle rate is 2^-16.07 per step at this s; the charged eps keeps the
    conservative 5/(12R) term. The long-cycle allowance is at theta = 2^-72 (below 2^26.66 events, below
    2^103.03 steps); the drift table lists 2, 2.5, 3, 4 and 5 points where the record listed 1, 1.5, 2, 3
    and 4.
10. *Squaring fraction.* s = 3/8 -> **37/64** (a step squares iff T mod 64 < 37, a 6-bit field; the
    record's rule is T mod 8 < 3). At these prices more squaring lowers the average step price
    (c = 0.59003046875) by more than the collision law raises the step count: sum q^2 = 0.334402 (record
    0.141006), log2 N' = 251.859405 (rounded up; record 251.4914), S = 2^125.9685 (record 2^125.8302). At theta = 2^-71
    the expected reports would be 1.028 of the cap, at 2^-72 more than 1.94 below. Section 7's tables are
    at s = 37/64.

Part D prices, through one code path, the earliest vector (11,400/300/2,200: c = 0.9225,
W = 2^125.71472902377629) and the record's (11,369/270/2,119: c = 0.9178875, W = 2^125.7074975), asserting
both at their own s = 3/8 and bucket/window, then this package's kernels at s = 37/64.

### 0-i..0-iv Machine details (moved)

The four machine-detail subsections (counters and work meter, report log and matcher, restart cap, initial values) moved whole to `appendix.md` sections A-i..A-iv to keep this file under the verifier per-file byte cap; all `appendix.md section A-i/0-ii/0-iii/0-iv` references in the package mean `appendix.md` sections A-i/A-ii/A-iii/A-iv. No number changed in the move.
## 0a. Representation: which limbs the rule reads

A walk element is held as the six-limb Montgomery vectors that the pinned blst kernels store,
x~ = x*2^384 mod p. The decision rule is stated on those stored limbs, which is what a limb-level
implementation of the design reads at zero cost (reading an integer residue instead would need a
Montgomery reduction per T read and twelve per multiply step, work no charge here contains):

* T is the six limbs of Re(c_0) **as stored**. The branch bits are (stored T) mod 64, the table index is
  (stored T >> 6) mod R, and the distinguished-point predicate is read from the same word.
* kappa is the minimum of the twelve **stored** key products a*zeta^j and d*zeta^j, and the twist index
  is the smallest index attaining it.

This is the rule on ordinary integers composed with a fixed bijection of Z/p, so every structural
property carries over: the twelve stored keys are distinct, their multiset is the same for all twelve
conjugates (Frobenius permutes it), so kappa is class-determined and the twist index is covariant;
stored T is invariant across conjugates because Re(c_0) is. The selected multiplier index is a
pseudorandom function of the class under either rule, which is all H1 uses. No conversion is charged.

**The squaring branch returns the same stored form as the multiply branch.** blst keeps every Fp
value fully reduced: the double-width reduction tail `__redc_tail_mont_384`
(`build/elf/mul_mont_384-armv8.S`, lines 1336-1372) and the end of `__mul_mont_384` both finish with
one conditional subtraction of the modulus, and `__add_mod_384` and `__sub_mod_384`
(`build/elf/add_mod_384-armv8.S`) end with the matching conditional correction (`csel ..., lo` after
a trial subtraction, or an `and`-masked add-back), so every stored coordinate is the unique residue in
[0, p) times 2^384. An element of Fp12 therefore has exactly one stored limb vector, and
`blst_fp12_cyclotomic_sqr(a)` and `blst_fp12_mul(a, a)` return identical limbs; for blst this is an argument from
the source (`measure/disassembly-audit.txt` section 5), as no limb-by-limb run of those two routines was made.
The walk prices and runs the two shipped kernels instead, `kernels/mul.s` and `kernels/sqr.s`: `gen.py`
bounds every output into [0, p) by its conditional subtractions, and `kernels/check_kernels.py` executes
both in the emulator (46 multiply and 20 cyclotomic-squaring cases) and checks every output limb equal to
the canonical Montgomery limbs x*2^384 mod p of a big-integer reference, so the two branches the walk
executes return the canonical stored form on every checked case. The
shipped Python walk computes in ordinary integers and applies the Montgomery factor only to the values
it reads (part C, stored rule), which gives the same bits by the uniqueness just stated.

Evidence, part C, both rules, 2,000 production-field steps each: the stored rule
gives 0 class-check failures in the first 150 steps and the exponent invariant y = g^alpha z^beta holds
at the end; its squaring fraction is 0.5824 over 1,997 branch-rule decisions (2,000 steps less 3
two-cycle escapes), 0.42 points above 37/64, z-score 0.38 (section 3.4), chi-square of the low 6 bits
59.56 (63 dof) and of the next 4 bits 16.57 (15 dof); it visits 1,997 distinct T values, the last new one
at step 1,999, and shows no repeating period at the end. The integer rule gives 0 failures and the
exponent invariant, with squaring fraction 0.5877 over 1,997 decisions, z-score 0.86, chi-squares 57.26
(63 dof) and 15.43 (15 dof), 1,996 distinct T values, the last new one at step 1,999, and no repeating
period. The selected key index kappa differs from the integer rule in 755 of the stored rule's 834
multiply steps, so the two rules are different walks and the claim is about the one the charges can
pay for.
`experiment.py` reports both blocks (`partC_production_walk`, `partC_production_walk_integer_rule`).

## 1. Reduction

Compute g = e(P0, Q) with the pinned convention pairing(Q, P0). Q has prime order r and
the pairing is nondegenerate, so g != 1 generates GT (prime order r). z = g^x for a unique
x in [0, r). Find x, output P = [x]P0, check e(P, Q) = z. If z = 1 output infinity. The
answer is for the supplied Q.

## 2. Arithmetic facts (checked in `experiment.py`, parts A and C)

Tower: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi), xi = 1+u (the blst tower). With
u = w^6 - 1 this is the profile field, w^12 - 2w^6 + 2 = 0, and conversion of coefficients
is free: c_m = a_m - b_m, c_{m+6} = b_m.

**F1 Frobenius classes.** (sum c_m w^m)^p = sum conj(c_m) d_m w^m with d_m = xi^{m(p-1)/6}.
p has order exactly 12 modulo r, p^6 = -1 and 1 + p^4 + p^8 = 0 (mod r), so the maps
y -> y^{p^k} form a group of order 12 acting on GT, containing inversion. Every y != 1
lies in a class [y] of exactly 12 elements; N = (r-1)/12 classes, log2 N = 251.2721.

**F2 Keys.** Put zeta = xi^{(p^2-1)/6}; it lies in Fp and is a primitive 6th root of unity
(zeta^2 = zeta - 1). Let c_1 = a + b u be the w^1 coefficient of y and B = Re(conj(c_1) d_1)
= a Re(d_1) + b Im(d_1). The w^1 coefficient of y^{p^(k+1)} is conj(c^{(k)}) d_1, and
conj(d_1) d_1 = zeta, hence
`key_k := Re(w^1 coeff of y^{p^k}) = a zeta^j (k = 2j), B zeta^j (k = 2j+1)`.
With zeta^j in {1, zeta, zeta-1, -1, -zeta, 1-zeta} all 12 keys cost **4 Fp multiplications**
(a zeta, a Re d_1, b Im d_1, B zeta) plus 8 subtractions/negations and one addition. The
12 keys are distinct unless two of them coincide by a linear accident (probability <= 66/p
per element; full-width ties resolve to the lowest index by the specified minimum (section 5) and are charged as failure, a loss below 2^-247.9 over a run's below 2^126.66 keyed
elements, section 4). K([y]) = min_k key_k is a class invariant,
kappa(y) = argmin_k key_k, and y^{p^kappa} is the same element for all 12 members of the
class (the canonical representative y^). K is the minimum of 12 values that behave like
independent uniform residues (they come in negation pairs, so K lies in [0, p/2]); its low 200
bits are uniform up to statistical distance about 2^200/p = 2^-181 (the density of the minimum is smooth).

**F2b Trace invariant.** The w^0 coefficient of y^{p^k} is conj^k(c_0) (no twist), so T([y]) := Re(c_0) is the same
for all 12 conjugates, and 12 Re(c_0) is the trace of y to Fp. T is read directly from the coefficients (no
multiplication) and, for distinct classes, is an essentially uniform residue (collisions of T between distinct classes
have probability about 1/p per pair). It plays the role of the class identity for the bits below; K is needed only to
find kappa. `experiment.py` part A checks that Re(c_0) is equal over all 12 conjugates of random GT elements.

**F3 Cyclotomic squaring (Granger-Scott).** GT lies in the cyclotomic subgroup, where
squaring is cheaper than a product. Write y = A + B w + C w^2 with A=(c0,c3), B=(c1,c4),
C=(c2,c5) in Fp4 = Fp2[v]/(v^2-xi), v = w^3, and (A0,A1) = A^2, (B0,B1) = B^2,
(C0,C1) = C^2 computed in Fp4 with 3 Fp2 squarings each. Then
`c0' = 3A0 - 2c0, c3' = 3A1 + 2c3, c1' = 3 xi C1 + 2c1, c4' = 3C0 - 2c4, c2' = 3B0 - 2c2, c5' = 3B1 + 2c5.`
`experiment.py` verifies y^2 = formula on real BLS12-381 GT elements and counts the
operations of the implementation: **9 Fp2 squarings and 34 Fp2 additions/subtractions**
(the multiplications by xi are counted as one Fp2 addition).

blst's `blst_fp12_cyclotomic_sqr` (`src/fp12_tower.c` lines 630-660) is this map: three `sqr_fp4` calls
and the 3t - 2u / 3t + 2u combinations, on blst's Fp12 = Fp6[w]/(w^2 - v) layout, in which our coefficient
c_m of w^m is the Fp2 element a[m mod 2][m div 2]. Section 5 prices it by measurement.

**F4 Orders.** r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 *
2529403 * 52437899 * 254760293^2 (verified by multiplication in `experiment.py`). The order
of 2 modulo r is (r-1)/16 (about 2^250.9), so 2^t lies in the order-12 group <p> only if t is a
multiple of about 2^247.3.

## 3. Algorithm

Parameters: R = 1024 table multipliers, squaring probability s = 37/64, theta = 2^-72 (DP
density), M = 2^20 walks round-robin on one processor (appendix.md section A-i). Trails do not depend on the order walks are
advanced, so parallel collision search's analysis (M walks, distinguished points) and section 4 apply: the schedule
changes which walk's steps occur first, not the step budget or any trail, and a collision of two walks is found at
matching when both trails reached the distinguished point within the walk phase. Bits of the class invariant T = Re(c_0) (read as the six stored Montgomery limbs of section 0a):
T mod 64 selects
the branch (squaring iff T mod 64 < 37), h = (T >> 6) mod 1024 selects the multiplier, the class
is distinguished iff (T >> 16) mod 2^72 = 0 (bits 16 to 87 of T). The fingerprint of a class used for matching reports is
F = (T >> 88) mod 2^128 (bits 88 to 215 of T), the 128 bits above the branch, multiplier and DP bits, modeled as independent under H2;
T has 381 bits, so every bit used is a bit of the stored T and none is shared between the roles.

Precompute M_h = g^{a_h} z^{b_h} for h < R (fresh uniform a_h, b_h mod r, drawn by the sampler of
appendix.md section A-iv) and
T[h][t] = M_h^{p^t}, t = 0..11, is the T table (1,024 x 12 = 12,288 entries of 576 bytes, 7.1 MB);
the scalar tables a_h p^t, b_h p^t mod r hold 2 x 12,288 = 24,576 entries of 32 bytes (0.8 MB in all).

A walk keeps some member y of the current class together with (alpha, beta) mod r such that
y = g^alpha z^beta. Each walk starts at fresh uniform (alpha, beta), two sampler draws of appendix.md section A-iv. One step:

```
1  T <- Re(c_0) of y (no multiplication); branch, h and DP bits from T
2  if class distinguished: compute the 12 keys (4 Fp mult) -> kappa; append (F, p^kappa*alpha, p^kappa*beta)
   to the report log of appendix.md section A-ii
3  2-cycle / checkpoint tests on T (section 3.3)
4  if T mod 64 < 37 (squaring branch):  y <- y^2 (Granger-Scott);  alpha,beta <- 2*alpha, 2*beta
   [no keys, no twist]
   else (multiply branch): compute the 12 keys (4 Fp mult) -> kappa; h <- (T>>6) mod 1024;
        t <- (12 - kappa) mod 12;  y <- y * T[h][t]   (ONE dense Fp12 multiplication);
        alpha += a_h p^t;  beta += b_h p^t
```

**3.1 The step is a map on classes.** Write y^ = y^{p^kappa}. Multiply branch:
y * M_h^{p^{-kappa}} = (y^ * M_h)^{p^{-kappa}} lies in the class [y^ M_h]. Squaring branch:
y^2 = (y^^2)^{p^{-kappa}} lies in [y^^2], because squaring commutes with Frobenius. Both
depend only on the class (T is a class function; kappa and y^ enter only through the conjugation
absorbed by the twist), so walks that reach the same class stay together and reach the same DP.

**3.2 Collision.** The log keeps every report until the cap (records are appended, never replaced or
dropped), so a useful collision is lost to the store only by a cap failure. After the walk phase the
matcher of appendix.md section A-ii sorts the log and takes the pairs inside each run of equal F. Two records with equal F
are the same class, up to
a false match of probability 2^-128 per pair of distinct classes: with D <= 2^55 reports there are at most
D^2/2^129 = 2^-19 false-key pairs in expectation; every pair's candidate is verified (one modular inverse and
g^x = z, <= 600 units), which is negligible. Their canonical exponents (alpha^, beta^) = p^kappa (alpha, beta) satisfy
alpha^1 + x beta^1 = alpha^2 + x beta^2 mod r, so x = (alpha^1 - alpha^2)/(beta^2 - beta^1).
beta^1 = beta^2 happens with probability 1/r. Output P = [x]P0 after checking g^x = z.

**3.3 Fruitless cycles.** A closed walk of the class graph is *useless* if the canonical
exponent pair returns unchanged. After l steps the pair is Lambda (alpha^_0, beta^_0) + tau
where Lambda = 2^t p^m (t squaring steps, m a sum of twists) and tau = sum_h c_h (a_h, b_h)
with c_h in Z[p]. Because (alpha^_0, beta^_0) is uniform and independent of the table, the pair
returns unchanged with probability about 1/r unless Lambda = 1 and every c_h = 0.
* **No squaring step can occur in a useless cycle.** Lambda = 1 with t >= 1 needs
  2^t in <p> (F4), i.e. t >= 2^247.
* A useless cycle made of multiply steps needs every table index used to occur at least twice
  with coefficients p^{e_1} + p^{e_2} + ... = 0 in Z[p]/(Phi_12(p)): for two occurrences
  e_2 = e_1 + 6 (probability 1/12), for three {e, e+4, e+8}, and so on; for any number of occurrences
  the probability is at most 1/12 (the last twist must hit one value). For 2-cycles this gives the
  model rate (27/64)^2/(12 R) = 2^-16.07 per step at the production R = 1024, below 1/(12R); at the
  census's R = 64 the prediction is 2.32e-4 and the measured rates are 3.27e-4 (q = 67, 3,000 runs) and
  2.21e-4 (q = 163, 900 runs), with 2.65e-4 in the declared 40-trial run. For l = 2k indices drawn independently and
  uniformly (heuristic H2) the probability is at most (2k)!/(2^k k!) (12R)^-k (pairings dominate), which
  is below 2^-288 for 65 <= l <= 147. For l >= 148 all l consecutive steps must be multiply steps, which
  has probability (27/64)^l < 2^-184 (the branch bits are bits of T, H2; the charged allowance keeps the
  conservative 2^-100 rate). Such cycles are never met in
  practice; if one is and holds no distinguished class, the stuck-walk rule restarts the walk after 20/theta = 20*2^72 steps, which adds
  less than 2^-13 relative work even at the restart cap of appendix.md section A-iii (2^-13.69, section 5); one that holds a distinguished
  class is not restarted, a failed run (appendix.md section A-iii).
* **Detection.** The step before last is remembered; if its T equals the current T the two
  classes form a 2-cycle and the walk replaces the representative of the class with the smaller
  T by its square (y -> y^2, exponents doubled), cost <= 3 units. This depends on the cycle
  only, so merged walks stay merged. A checkpoint T is stored every 64 steps and compared at
  every step; a repeat means a cycle of length <= 64, which is traversed once more to its minimum-T
  class that is then squared (waste <= 200 steps). Every other walk whose stuck-walk
  counter reaches 20/theta = 20*2^72 steps is restarted from a fresh uniform pair (stuck-walk rule). The counter is reset to zero
  at every walk start, restart and report (appendix.md section A-i), so a walk is restarted only after 20/theta steps with no report, and the
  rule is checked once per report gap.
* Total overhead factor charged on the step count: eps = 5/(12R) + 250/R^2 = 6.5e-4, unchanged and conservative at this s (5/(12R) is more than 28 times the model 2-cycle rate (27/64)^2/(12R) = 1.45e-5).

**3.4 Caps and the work cap.** The run has three phases: setup (at most 2^31 units, section 5), the walk
phase, and the stopping phase (matching, output, verification, halting). The stopping phase is bounded in
advance: the heapsort of at most 2^55 records, the dedup, at most 2^16 pairs at 600 units, [x]P0 and the
verification pairing, and 100 instructions for the failed compare, the hand-over and halting, together the
**stopping allowance R_stop, below 2^53.55 units** (`experiment.py` part D). The walk phase has the budget
**B_walk = W_cap C - 2^31 C - R_stop C instructions**, with W_cap = 2^125.234 the submitted bound and
W_cap C rounded down to an integer. The walk phase meters its own charged work with the prices of section 5
(0.8871 units per multiply step, 0.37325 per squaring step, the stated costs of escapes, switches, starts,
restarts and reports) in the 3-word (192-bit) meter of appendix.md section A-i under its **pre-debit rule**: before each
charged action the meter is compared with L = B_walk - A_max, the action runs only if the meter is below L,
and the meter then adds the action's charge, at most A_max. So the walk phase charges below B_walk, and
when a compare finds the meter at or above L the walk phase ends and the stopping phase begins. Every run therefore does at most
2^31 C + B_walk + R_stop C <= W_cap C instructions, a deterministic inequality that `experiment.py` part D
checks on integers. The walk phase can always spend L (asserted); the probability that the threshold binds early is part
of the failure budget below.

*Heuristic H1+ (no drift of the walk statistics).* The cap is sized to cover a drop of the squaring
fraction by 2 percentage points. A multiply step costs 0.8871 units and a squaring step 0.37325, a
difference of 0.51925 units, so a drop of D points raises the average step price by 0.0051925 D units, a
factor 1 + delta with delta = 0.0051385 D / 0.59003046875 = 0.87089% per point (2 points: 1.74177%). The nominal
ledger with every step priced that much higher is 2^125.2333219 (`experiment.py` part D), the least
thousandth of a bit at or above it is 125.234, and the exact integer inequality
(numerator^1000 <= denominator^1000 * 2^125234 for that ledger rounded up to 40 digits, and the strict
reverse at 125233, so the claim is the least thousandth) is asserted there.
The cap therefore exceeds the nominal ledger W_nominal = 2^125.2084098 by at least 1.79% (0.0256 bits); the
walk budget's per-step slack exceeds 0.010560 units, and the threshold does not end the walk phase before the ledger's walk is done while
the squaring fraction stays within 2.0549 percentage points of 37/64. That allowance bounds the squaring fraction averaged
over the whole run (H1+), which a short sample does not measure. The part-C walk (stored rule) reads 0.5824 over 1,997
branch-rule decisions, 0.42 points above 37/64, with printed binomial standard deviation 0.0111 at 37/64
and z-score 0.38. Part B's toy walks (toy hash rule) read 0.5735 in the q = 67 census and 0.5774 in
the q = 163 census (0.5716 over 44,647 branch-rule steps in the declared 40-trial run), 0.47 points or
less below 37/64 in the censuses and 0.66 in the declared run; all are inside the 2.0549-point allowance. A larger drop ends the walk phase early; the shorter walk's success is still the section-4 collision law, and
the budget buys n/(1+delta) steps of the nominal n:

| drop of the squaring fraction | success (lower bound) | loss from 0.41 (upper bound) |
| ---: | ---: | ---: |
| up to 2 points (the cap's design point, 2.0549 points) | 0.4104 | 0 |
| 2.5 points | 0.4080 | 0.0020 |
| 3 points | 0.4054 | 0.0046 |
| 4 points | 0.4002 | 0.0098 |
| 5 points | 0.3951 | 0.0149 |

H1+ states that over the whole run the branch, 2-cycle and DP frequencies stay inside a drift of this size
(the 0.01 budget is not exhausted until about 4.02 points) and the reports, false candidates and restarts
stay within the caps below, with total probability of any cap failure or early threshold at most 0.01. This is a much weaker
statement than independence of the steps; it excludes only a drift by percentage points over 2^125 steps.
Support: (i) H2 models the low six bits as uniform, giving nominal squaring probability 37/64; exact uniformity over production classes is not proved;
(ii) in the toy runs the fluctuation of the mean charge scales as the independent-steps prediction:
rms of (relative deviation x sqrt(steps)) is 0.4395 in the q = 67 census (3,000 runs, seed 2601) and
0.4199 in the q = 163 census (900 runs, seed 2602) against 0.4313 predicted (the declared 40-trial run
reads 0.3899), so there is no excess correlation of branch labels at that scale; (iii) the part-C walk
(stored rule) reads 0.5824, z-score 0.38, beside the census's toy fractions 0.5735 and 0.5774 (above). The toy runs cannot test
concentration at 2^125 steps; H1+ extrapolates it, exactly as H1 extrapolates the collision law. The
success analysis (section 4) charges 0.01 for it.

*Caps.* reports D <= 2^55, counted in a one-word counter (appendix.md section A-i); the expected number of reports is theta
times the steps the walk budget buys, below 2^54.05 at the nominal average price, more than a factor 1.94 below the cap, and
reaching the cap fails the run; at most 2^16 charged pairs in the matcher (the expected number of false-key
pairs is at most 2^-19, so Markov gives probability 2^-35 of an overflow, which fails the run; the
pair work is at most 2^16 * 600 units < 2^25.24); **at most 2^36 restarts** (appendix.md section A-iii: expected restarts below 2^27.30, below 2^-8.65 of the cap; a run past the cap is a cap failure inside the
0.01 of H1+; each restart a fresh start plus 40
instructions of meter bookkeeping, 2^27.04 units in total, and at worst one stuck-walk allowance each,
2^36 * 20/theta < 2^112.33 steps, charged in the step bucket of section 5); the total number of steps
that the ledger needs is S(1+eps) + (M + 2^36) 20/theta. The walk phase runs until its budget is spent,
so the detection window and restart allowance are steps the budget buys,
and the steps beyond the ledger's come from the cap's slack over the nominal ledger, at least
2^119.41 units (1.79%). Fruitless cycles
longer than 64 steps occur at a rate below 2^-100 per step (section 3.3; below 2^26.66
events in expectation; one with no distinguished class wastes at most 20/theta = 20*2^72 steps before the stuck-walk restart, below 2^103.03 steps
in total in expectation; one through a distinguished class is a failed run, appendix.md section A-iii).
These caps bound the work and memory stated below.

## 4. Success probability

**H1 (random mapping with branch mix).** The class sequence behaves like iterates of a
random mapping on N classes in which two distinct classes have the same image with
probability (1 - sum_b q_b^2)/N, where q_b are the branch probabilities. Justification: for
c1 != c2 the equation F(c1) = F(c2) has, for each branch pair (b1,b2) with b1 != b2, about one
solution among the N classes (12 candidate representatives, each canonical with probability
1/12); for b1 = b2 = squaring it has none (squaring is a bijection of GT of odd order that
commutes with Frobenius), and for b1 = b2 = a table index it has about 11/12 of one. We
use the conservative value sum q^2 = s^2 + (1-s)^2/R = 0.334402 for all pairs; the exact
value (s^2 + (1-s)^2/(12R)) is smaller by 1.6e-4. Part E of `experiment.py` tests this law directly (section 7). This is the same correction as the Brent-Pollard
factor R/(R-1) of the pure adding walk (s = 0), and for Pollard's original three-branch walk
it is the mechanism that makes Pollard's original square-and-add walk slower than an r-adding walk
(Teske 1998). Hence the effective class count is
N' = N / (1 - sum q^2), log2 N' < 251.859405.

Under H1 the probability that the first S steps contain a class collision is
1 - exp(-S^2/(2N')). Choose P_model = 0.41 (section 0):
`S = sqrt(2 N' ln(1/0.59)) = 2^125.9685.`
Losses before the first useful collision is detected: detection failure on either of the two selected
trails (the colliding one and the one it runs into), each 20 e^-20/(1 - e^-20) < 4.13e-8 (length-biased;
the fresh-trail e^-20 = 2.07e-9 does not apply), charged as 2 x 4.13e-8, a key tie (two of an element's 12
keys equal, at most 66/p per keyed element over a run's below 2^126.66, below 2^-247.9, section 2), a
collision inside a single walk (fraction about 1/M = 1e-6), beta^1 = beta^2 (2^-255), cap
failures and the threshold ending the walk phase early (at most 0.01 together, H1+, section 3.4; the report cap and
the pair budget are among these caps), steps spent on restarts (the restart cap's own worst case: 2^36
restarts at one full stuck-walk allowance of 20/theta = 20*2^72 steps each are 2^36 * 20/theta < 2^112.33
steps; part D computes the exact fraction < 4.85e-5 of the walk's hard maximum (part D), subtracted as 4.85e-5;
restart counts and costs, appendix.md section A-iii), a long cycle through a
distinguished class, a failed run, probability below 2^-121.61, subtracted as 2e-11, and the sampler's
exhaustion allowance, at most 2^-218.99 under the ideal-source model of appendix.md section A-iv, subtracted as
2^-218.99. The report log loses no report before the cap (appendix.md section A-ii), so the store adds
no loss term. The sum is the 0.01 of cap failures plus the restart and small terms: below 0.0100496, **below 0.0101**
(experiment.py part D). The success after every loss is therefore at least 0.41 - 0.0101 = 0.3999,
leaving **0.0099** over 0.39 for the inaccuracy of H1. Toy measurements (section 7) of the fraction of
runs with a useful collision within S(0.43) are 0.4317 +- 0.0090 (3,000 runs, q = 67, seed 2601) and
0.4133 +- 0.0164 (900 runs, q = 163, seed 2602), pooled 0.4274 +- 0.0079 over 3,900 runs (0.33 standard
errors below 0.43), at 0.423 and 0.402 for S(0.42); the 0.0099 slack after every loss is 1.3 pooled
standard errors. Individual settings of a few hundred runs scatter by +-0.04.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Counts are executed scalar AArch64 instructions as the cost model counts them
(`tracks/docs/COST_MODEL.md`): Valgrind Callgrind, Debian 13.7, GCC 14.2.0, blst f7a88067 built `-O3
-fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only`, collection around the call loop, marginal
(Ir(65)-Ir(1))/64; `measure/` holds the record and part D reads it.

| routine | per-call marginal | status |
| --- | ---: | --- |
| `blst_fp_mul` | 383 | measured; equals the cost model's reference count |
| `blst_fp12_mul` | 22,996 | measured; equals the cost model's reference count |
| `blst_fp12_cyclotomic_sqr` | 9,803 | measured, on real GT elements, identical at five slice offsets |

The first two lines reproduce the cost model's own two figures to the instruction; the third is the blst baseline
that this package's squaring kernel replaces.

| walk kernel | executed instructions | status |
| --- | ---: | --- |
| `kernels/mul.s` (x * fixed entry) | 15,492 | executed count, branch-free, 0 mismatches in `kernels/check_kernels.py` |
| `kernels/sqr.s` (Granger-Scott) | 7,226 | executed count, branch-free, 0 mismatches in `kernels/check_kernels.py` |

**Table multiply:** the product y*T[h][t] uses `kernels/mul.s`, **15,492 executed instructions** (0.7746
units), below blst's dense `blst_fp12_mul` reference operation (22,996).

**Squaring branch:** the Granger-Scott map of F3 is what blst's `blst_fp12_cyclotomic_sqr`
(`src/fp12_tower.c` lines 630-660) computes. The straight-line kernel `kernels/sqr.s` computes it in
**7,226 executed instructions** (0.3613 units), charged in full; blst's measured 9,803 for the same map is
the baseline it replaces. The squaring update needs no keys; DP report keys are charged separately.

**Kernels.** `kernels/gen.py` emits both: general-purpose registers only, no NEON, no branch, no call
(`ldp`/`stp` allowed), Montgomery form, 6 x 64-bit limbs, outputs fully reduced. Because the programs are
branch-free and call-free, the executed count equals the emitted count on every call. `kernels/emu.py`
interprets that subset, and `kernels/check_kernels.py` regenerates both listings and runs a reduced sample
against a big-integer reference, asserting 15,492 and 7,226 executed and 0 mismatches. The 54 Karatsuba pre-sums stored per table entry (2,592 bytes) are built in the table construction, inside its 2^20.2-unit item. The listings (363 KB, 169 KB) are regenerated, not shipped.

**Per-step work.** Every step (both branches), analytic, as itemised before:

| Work | Instructions |
| --- | ---: |
| read Re(c_0) low two limbs (ldp), high word only on DP-hit path, branch bits, DP orr-chain, table address | 42 |
| (alpha, beta) update mod r (two 4-limb add + conditional subtract, or doubling; 32 each) | 64 |
| 2-cycle and checkpoint comparisons on T, history | 60 |
| loop control: back-branch (1), chunk counter subs/b.eq (2), branch select on the 6-bit field, and #63/cmp #37/b.lo (3), DP test branch (1), 2-cycle and checkpoint test branches (4), call and return of the step routine with frame save and restore (21) | 32 |
| margin for loads and moves | 28 |
| **subtotal: base bookkeeping recounted** | **226** |
| 192-bit work meter: step charge in a register (1), adds/adcs/adc (3), exact 3-word compare with the cap, cmp/sbcs/sbcs/b.hs, on every step (4) (appendix.md section A-i) | 8 |
| 2-word stuck-walk counter: adds/adc increment (2), mov of 5120 into x16 (1), cmp of the high word with x16 and b.hs (2), exact as 20/theta = 5120 * 2^64 (appendix.md section A-i) | 5 |
| **total counted 239, charged** | **239 = 0.01195 units** |

Multiply steps in addition (the keys give kappa), analytic, built on the measured blst_fp_mul = 383 and
a modular add/sub/neg of 6 limbs inlined at 30 instructions (6 limb ops, fused conditional-subtract and
select, loads and stores; no call, return or address set-up):

| Work | Instructions |
| --- | ---: |
| 4 Fp multiplications (a zeta, a Re d_1, b Im d_1, B zeta) | 1,532 |
| 1 Fp addition (B) + 8 Fp subtractions/negations, inlined at 30 each | 270 |
| minimum of the 12 keys on the stored limbs (seeded top-limb load 12; 11 challengers x (top-limb load 1, cmp 1, csel triple 3, full-width descent on ties: b.eq 1, up to 5 lower-limb compares, deciding csel triple 3, tie flag 1)) | 177 |
| load of the w^1 coefficient, kappa, twist index | 32 |
| **total counted 2,011, charged** | **2,011 = 0.10055 units** |

The minimum is the exact procedure `Tower.min12` of `experiment.py`: the running best always updates on
`challenger <= best` with the strict-less limb flag, so the comparison is a full-width six-limb compare
and full-width ties (identical limb vectors, at most 66/p per keyed element, section 2) deterministically
keep the lowest index and set the tie flag. The worst case above charges every challenger with the full
descent, so the 177 bounds every run including the all-tied one; no abort-on-tie path exists. Part A runs
this procedure on production GT elements and asserts it agrees with the reference minimum with no ties.

So a multiply step costs (15,492 + 239 + 2,011) / 20,000 = **0.8871 units** and a squaring step
(7,226 + 239) / 20,000 = **0.37325** units.
Average step cost c = (27/64) * 0.8871 + (37/64) * 0.37325 = **0.59003046875 units**.

Setup (at most 2^31 units, the preprocessing bound): pairing g (<= 2^14), check z^r = 1 (<= 510), the
table (1024 multipliers x (2 exponentiations <= 1,021 units + 11 Frobenius maps) + 2,048 sampler draws
x 29.696 units, the scalar tables and the 12,288 x 54 pre-sums at 200 instructions, <= 2^20.2), the sampler's state (<= 2^12 instructions), M = 2^20
starts g^alpha z^beta (<= 2^30.1, each with two sampler draws), and the initial counter values (at most
32 instructions, appendix.md section A-iv); their sum is below 2^30.11 units. DP handling: at most 2^55 reports (the cap of section 3.4) x **4,000 instructions**,
itemised as: 12 keys for kappa 2,011; two 255-bit modular multiplications by p^kappa mod r 1,600; the pre-debit
compare of the meter with the threshold and the charge's 3-word add 8; the one-word report counter compare
against the cap, branch and increment 5; the slot address base + 80 c 2; appending the 80-byte record (ten 8-byte
stores and address set-up) 12; reset of the 2-word stuck-walk counter 4; operating-system zero-fill of the
log's 2 MiB pages on first write (131,072 stores and one fault per 26,214 records) 13; margin for loads, moves,
call and return 345. That is 2^55 * 4,000 / 20,000 < 2^52.679 units. Meter bookkeeping (appendix.md section A-i) adds 2^20 * 40 / 20,000 < 2^11.04 units for the walk starts and 2^36 * 40 / 20,000 < 2^27.04 for the
restarts, the restart counter adds below 2^23.30 units (appendix.md section A-iv), and the restarted walks' fresh starts add
2^36 * 2^10.1 = 2^46.1 units (each start's two sampler draws are inside its 2^10.1 bound). Matching, output and halting, the stopping allowance of section 3.4, at their worst
case: the heapsort of 2^55 records (7,188 instructions per record, appendix.md section A-ii, below 2^53.53 units), the dedup (2^55 x 100 instructions), at most
2^16 pairs at 600 units (below 2^25.24 units), [x]P0 (<= 510 G1 operations, each < 1 unit), the verification
pairing (<= 2^14) and 100 instructions of hand-over and halting, together below 2^53.55 units. All non-walk work
other than the switches (appendix.md section A-i, charged separately below) is therefore below
2^53.55 + 2^52.679 + 2^46.1 + 2^31 + 2^27.04 + 2^23.30 + 2^11.04 < 2^54.2 units (`experiment.py` part D
sums the exact terms, below 2^54.181, and asserts the sum below 2^54.2).

Detection delay: after the first useful collision the ledger charges 20/theta further
steps for each of the M walks, E = M * 20 / theta < 2^96.33 steps. The restart cap of appendix.md section A-iii adds one stuck-walk allowance per
restarted walk, 2^36 * 20/theta < 2^112.33, so the step bucket outside the walk is (2^20 + 2^36) * 20/theta < 2^112.33 steps
(2^-13.69 of S). As 20/theta = 10,240 K, the colliding walk gets its 20/theta further steps while each
other walk gets at most 10,240 chunks, 20/theta, so the window holds; switches add below 2^58.65 units (appendix.md section A-i).

`W = c * S * (1 + eps) + c * (2^20 + 2^36) * 20/theta + 2^54.2 + 2^58.65`, evaluated by `experiment.py` part D with the exact S = 2^125.968505..., 1 + eps = 1.000645319... and (2^20 + 2^36) * 20/theta = 2^112.321...: **W < 2^125.2084098**.

Submitted bound: **time_log2 = 125.234** (the least thousandth at or above the ledger recomputed with the
average step price raised by the 2-point drift factor of section 3.4, 2^125.2333219). It is at least 0.0256 bits (1.79%)
above the nominal ledger 2^125.2084098, which is the cap headroom of section 3.4. `experiment.py` (part D)
recomputes this ledger, and the two earlier price vectors of section 0 (2^125.71473 and 2^125.70750, at their own
theta = 2^-60 and the cost bucket input 2^63.8 of section 0) through the same code path, with exact decimal arithmetic, and checks the 125.234 inequality on integers.

Memory: the report log of appendix.md section A-ii is 2^55 slots x 80 bytes (128-bit fingerprint, two 256-bit
exponents, no pointer field) < 2^61.322 bytes, sorted in place with one 80-byte temporary, plus the T table
(12,288 entries of 576 bytes, 7.1 MB), its 54 Karatsuba pre-sums (2,592 bytes per entry, 31.8 MB), the
scalar tables (24,576 entries of 32 bytes, 0.8 MB) and 2^20 walk states (at most 2 KB each, including the
remembered previous state, 2^31 bytes): below 2^61.33 bytes, and the claim states below 2^63 bytes. Memory
is reported, not scored.

## 6. Heuristics and limitations

* H1 above (random-mapping behaviour of the mixed walk on classes), H2 (pseudorandom branch and DP bits
  of T), independence of table indices (used for cycle probabilities). The success probability is computed
  in the ideal model: a uniformly random mapping on the classes (H1), realized by table multipliers and walk
  starts drawn from the ChaCha20 generator; the analysis assumes that realization is indistinguishable from
  the ideal mapping, discloses this as a heuristic of the same class as H1 and H2, and claims no proved
  distributional or independence bound on the sampled scalars. Under the ideal-source model the bounded
  rejection is exact uniform conditioned on acceptance, and the exhaustion rule is the event charged at
  2^-218.99 in section 4. The experiment draws the same bounded reduction from Python's seeded Mersenne
  Twister, so its part S exercises the cap, compare and exhaustion rule, not the ChaCha20 block. These are
  standard in rho analyses (Teske 1998/2001; van Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999;
  Bos-Kleinjung-Lenstra 2010), supported here by toy and production-field runs, not proved.
* Scope of the prices. The two walk kernels' counts are executed counts from the shipped emulator, equal to
  the emitted counts because the programs have no branch or call; `kernels/check_kernels.py` checks a
  sample. The other charges are analytic bounds: the 239 of every-step work, the 2,011 of key work, the
  4,000 per report and the 132 per heapsort level. The blst prices (383, 22,996, 9,803) are measured with
  the cost model's toolchain and are the baseline; the walk itself was not compiled and run as a blst
  program, and the equal stored limbs of the two branches (section 0a) are argued from the blst source; the shipped
  kernels' outputs are checked equal to canonical Montgomery limbs by `kernels/check_kernels.py`.
  Sensitivity: a squaring step at 1.0 unit gives 2^125.9483; each +100 instructions on the
  squaring kernel adds 0.0070 bits and each +100 of per-step work 0.0121 bits.
* Collision statistics are validated in toy fields only. This revision's census (3,900 runs: q = 67,
  seed 2601, 3,000 runs; q = 163, seed 2602, 900 runs) reads mean T_c/sqrt(N) = 1.5565 and 1.5422 against
  the model's 1.5392 (+1.1% and +0.2%). The declared 40-trial run reads 1.3675 (-11%); 40 trials is a
  single small sample and the census's agreement is not claimed for it.
* The report log is one 2^61.322-byte array with 64-bit addressing (appendix.md section A-ii). Memory is
  unscored, and the claim is an instruction-count bound: no physical machine maps that much.
* This is a generic square-root attack. It does not use index-calculus (TNFS) methods in Fp12.
  It is a valid upper bound, not a security estimate.

## 7. Experiments (`experiment.py`, standard library only)

The experiments of parts A to E and S and their results, the offline toy aggregates and the pair-collision
tables, are in `evidence.md` (this section, moved whole to keep this file under the size limit).

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999); Teske, Math. Comp.
70 (2001); Duursma, Gaudry & Morain, ASIACRYPT 1999; Granger & Scott, PKC 2010; Aranha et
al., EUROCRYPT 2011; Bos, Kleinjung & Lenstra, ANTS 2010; mpjunior92, submission 7da5f23a.

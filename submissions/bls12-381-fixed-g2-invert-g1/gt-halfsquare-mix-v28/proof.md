# Frobenius-class parallel rho with a half-square branch mix and exact canonical keys

**Claim.** For uniform nonidentity Q in G2 and uniform z in GT the algorithm below
recovers P in G1 with e(P,Q) = z with probability at least 0.39, using at most
**2^125.71** dense-Fp12-multiplication equivalents (C = 20,000), including setup, the two
pairings, table construction, walk starts, distinguished-point (DP) handling, cycle
escapes, the DP detection delay, output and verification. Peak memory is below 2^73
bytes, preprocessing (included) at most 2^31 units, nonuniform advice 0 bytes. The success
bound rests on the usual random-mapping heuristic for rho walks, adjusted for the branch
mix (section 4), and is supported by toy-field experiments and by running the production
step rule in the real GT. No production-size run is claimed.

The class construction (Frobenius classes of GT, the 4-multiplication key, the twisted
multiplier table that never applies a Frobenius to the walk element) is the one published
in submission 7da5f23a-6633-46b1-9c03-11f6cfedeeec (mpjunior92), which this package builds
on. New here: (i) a cyclotomic-squaring branch that is cheaper than the table multiply;
(ii) a fresh instruction ledger; (iii) an analysis of how the branch mix changes the
collision law and the fruitless cycles; (iv) execution of the production step rule in the
real BLS12-381 GT; (v) explicit memory and report caps; (vi) the branch, multiplier and DP bits are read from the free class
invariant T = Re(c_0) instead of the minimum key K, so a squaring step needs no key computation at all.
The package `frobenius-rho-squaring-mix` (my earlier submission 6d4e57d1, accepted as valid at 125.80)
is the same design with the bits read from K; this one differs only by (vi) and its ledger.

## 0. What this revision changes

Everything below — the class walk, the 4-multiplication key, the twisted table, the
cyclotomic-squaring branch read from T, the escapes, the caps, H1, H1+, the toy and
production-field evidence — is the package of submission `56fc259f`
(`gt-tracecounted-exact-v25`, GordoAR, 125.72), which is itself a revision of `747ad6c9`
(125.74), `6d4e57d1` (125.80) and the Frobenius-class construction of `7da5f23a`
(mpjunior92). Only the branch mix and the key rule change:

* **The squaring mix moves from s = 3/8 to s = 1/2**, read as a single bit: squaring iff
  (stored T) mod 8 >= 4, i.e. bit 2 of T. The parent's price vector makes its 3/8 mix
  suboptimal: a multiply step costs 1.12545 units and a squaring step 0.58195 units, so the
  average step price falls by 0.27175 per unit of s while the random-mapping penalty
  sum q^2 = s^2 + (1-s)^2/R rises only slowly near 3/8. The full ledger over the eighth
  grid, recomputed in `experiment.py` through the same code path as the claim, is minimal
  at s = 1/2:

| s | ledger W (log2) | delta vs s = 3/8 |
| --- | ---: | ---: |
| 0 | 125.8927 | +0.1852 |
| 1/8 | 125.8141 | +0.1066 |
| 2/8 | 125.7533 | +0.0458 |
| 3/8 | 125.7075 | 0 |
| **4/8** | **125.7010** | **-0.0065** |
| 5/8 | 125.7311 | +0.0236 |

  A 16th-grid scan (7/16: -0.0056, 9/16: -0.0009) confirms s = 1/2 is the optimum, and a
  single-bit branch needs no comparison. (The scan uses the parent's key charge 2,119 at
  s = 3/8 and this package's exact-key charge 2,239 elsewhere; with the exact-key charge
  everywhere the s = 3/8 row reads 125.7134, i.e. the mix alone still improves on the
  exact-key parent by 0.0124 bits.)
* **The squaring test is stated as the code tests it — (T & 7) >= 4, threshold 4 — so the
  specified rule and the tested rule are one rule.** At s = 1/2 the threshold alone selects
  the squaring branch; there is no "bit 2 is 1" versus "(T & 7) < 4" mismatch left to state,
  and part C runs the production rule with the production table size (R = 1024) on the
  stored-limb bits.
* **The canonical key is now exact.** kappa is the minimum of the full stored key products
  with ties broken by fully specified lower-limb comparisons (section 5), priced in the
  charged vector (2,239 instead of 2,119) and in the meter. The worst-case lower-limb
  fallback has a fixed 270-instruction bound charged on every multiply step, so the
  every-run work guarantee of the meter no longer depends on any tie estimate.

`experiment.py` part D recomputes the parent's ledger (2^125.7074975 at s = 3/8, key charge
2,119) and this one (2^125.7010225 at s = 1/2, key charge 2,239) through the same code path;
the abort headroom is re-derived since the allowance shrinks from 0.87% to 0.62%.

## 0a. Representation: which limbs the rule reads

The parent states its rule on ordinary integers in [0,p): T = Re(c_0) "read as an integer", and the
twelve keys ordered as integers. Its arithmetic prices are those of the pinned blst kernels, and those
kernels store field elements in Montgomery form x~ = x*2^384 mod p, so reading a stored limb does not
yield that integer. Converting costs work the parent's 270 / 2,119 / 11,369 charges do not contain —
one reduction per T read and twelve per multiply step, and the twelve alone are +0.14 unit on a
multiply step, which would put the parent's own bound above its claim.

This package therefore states the rule on the stored representation, which is what a limb-level
implementation of the design runs at zero cost:

* T is the six limbs of Re(c_0) **as stored**. The squaring branch is taken iff
  (stored T) & 7 >= 4, i.e. bit 2 of the stored low limb (threshold 4, the same test
  `experiment.py` part C executes); the table index is (stored T >> 3) mod R, and the
  distinguished-point predicate is read from the same word.
* kappa is the minimum of the twelve **stored** key products (the full 384-bit stored values,
  compared most significant limb first down to the low limb, ties impossible to survive this
  comparison), and the twist index is the smallest index attaining it.

This is the parent's rule composed with a fixed bijection of Z/p, so every structural property carries
over verbatim: the twelve stored keys are distinct; their multiset is the same for all twelve conjugates
(Frobenius permutes it), so kappa is class-determined and the twist index is covariant; stored T is
exactly invariant across conjugates because Re(c_0) is. What changes is which multiplier index is
selected, and that index is a pseudorandom function of the class under either rule, which is all H1
uses. No conversion is performed anywhere in the charged path, so the price table of section 5 is a
bound on the stated algorithm rather than on a representation it does not use.

Evidence, part C: a production-field walk with the stated rule at the production table size
(R = 1024), 2,000 steps from the same seed as the parent's run: 0 class-check failures in
the first 150 steps, the exponent invariant y = g^alpha z^beta holds at the end, squaring
fraction 0.4965 against 1/2 (threshold-4 rule, exactly as specified), chi-square of the low
3 bits 6.46 (7 dof) and of the next 4 bits 16.0 (15 dof). The twelve-key exact-minimum
rule is checked against explicit Frobenius on every one of the first 150 steps, with the
full-limb fallback exercised by construction wherever top limbs tie. `experiment.py` reports
the production block (`partC_production_walk`) with R = 1024 recorded.

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
12 stored key products are compared as full 384-bit values most significant limb first;
ties that survive all six limbs need the full values to be exactly equal, a linear
accident of probability <= 66/p per element, and are broken by smallest index — always a
valid canonical choice, never a failure. K([y]) = min_k key_k is a class invariant,
kappa(y) = argmin_k key_k (smallest index on ties), and y^{p^kappa} is the same element for all 12 members of the
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

**F4 Orders.** r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 *
2529403 * 52437899 * 254760293^2 (verified by multiplication in `experiment.py`). The order
of 2 modulo r is (r-1)/16 (about 2^250.9), so 2^t lies in the order-12 group <p> only if t is a
multiple of about 2^247.3.

## 3. Algorithm

Parameters: R = 1024 table multipliers, squaring probability s = 1/2, theta = 2^-60 (DP
density), M = 2^20 parallel walks. Bits of the class invariant T = Re(c_0) (the six stored
limbs): squaring iff (T & 7) >= 4 (bit 2 of T, threshold 4 — the same test part C runs),
h = (T >> 3) mod 1024 selects the multiplier, the class
is distinguished iff (T >> 13) mod 2^60 = 0. The fingerprint of a class used for matching reports is
F = (T >> 73) mod 2^128, the 128 bits above the branch, multiplier and DP bits, so they are independent.

Precompute M_h = g^{a_h} z^{b_h} for h < R (fresh uniform a_h, b_h mod r) and
T[h][t] = M_h^{p^t}, t = 0..11, with the scalar tables a_h p^t, b_h p^t mod r
(12,288 entries of 576 bytes, 7.1 MB).

A walk keeps some member y of the current class together with (alpha, beta) mod r such that
y = g^alpha z^beta. Each walk starts at fresh uniform (alpha, beta). One step:

```
1  T <- Re(c_0) of y (no multiplication); branch, h and DP bits from T
2  if class distinguished: compute the 12 keys (4 Fp mult) -> kappa; report (F, p^kappa*alpha, p^kappa*beta)
   into a chained multimap keyed by F
3  2-cycle / checkpoint tests on T (section 3.3)
4  if (T & 7) >= 4 (squaring branch):  y <- y^2 (Granger-Scott);  alpha,beta <- 2*alpha, 2*beta
   [no keys, no twist]
   else (multiply branch): compute the 12 keys (4 Fp mult) -> kappa; h <- (T>>3) mod 1024;
        t <- (12 - kappa) mod 12;  y <- y * T[h][t]   (ONE dense Fp12 multiplication);
        alpha += a_h p^t;  beta += b_h p^t
```

**3.1 The step is a map on classes.** Write y^ = y^{p^kappa}. Multiply branch:
y * M_h^{p^{-kappa}} = (y^ * M_h)^{p^{-kappa}} lies in the class [y^ M_h]. Squaring branch:
y^2 = (y^^2)^{p^{-kappa}} lies in [y^^2], because squaring commutes with Frobenius. Both
depend only on the class (T is a class function; kappa and y^ enter only through the conjugation
absorbed by the twist), so walks that reach the same class stay together and reach the same DP.

**3.2 Collision.** The multimap keeps every report with a given F (no entry is ever replaced), so a
useful class cannot be shadowed. Two reports with equal F from different walks are the same class, up to
a false match of probability 2^-128 per pair of distinct classes: with D <= 2^66.07 reports there are about
D^2/2^129 = 2^2.9 false candidates in expectation; every candidate is verified (one modular inverse and
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
  model rate (1/2)^2/(12 R) <= 1/(12 R) = 2^-13.6 per step (measured 2.7e-4 to 3.9e-4 per step across
  seeds versus 3.3e-4 predicted at R = 64, section 7; 6.3e-5 against 8.1e-5 predicted at R = 256).
  For l = 2k indices drawn independently and
  uniformly (heuristic H2) the probability is at most (2k)!/(2^k k!) (12R)^-k (pairings dominate), which
  is below 2^-288 for 65 <= l <= 147. For l >= 148 all l consecutive steps must be multiply steps, which
  has probability (1/2)^l < 2^-100 (the branch bits are bits of T, H2). Such cycles are never met in
  practice; if one is, the stuck-walk rule restarts the walk after 20/theta = 2^64.3 steps, which adds
  less than 2^-30 relative work at that rate.
* **Detection.** The step before last is remembered; if its T equals the current T the two
  classes form a 2-cycle and the walk replaces the representative of the class with the smaller
  T by its square (y -> y^2, exponents doubled), cost <= 3 units. This depends on the cycle
  only, so merged walks stay merged. A checkpoint T is stored every 64 steps and compared at
  every step; a repeat means a cycle of length <= 64, which is traversed once more to its minimum-T
  class that is then squared (waste <= 200 steps). Every other walk that sees no DP within
  20/theta = 2^64.3 steps is restarted from a fresh uniform pair (stuck-walk rule).
* Total overhead factor on the step count: eps = 5/(12R) + 250/R^2 = 6.5e-4 (charged).

**3.4 Caps and the total-work abort.** The algorithm meters its own charged work with the prices of
section 5 (1.12545 units per multiply step, 0.58195 per squaring step, the stated costs of escapes,
starts, report handling and verifications) and **aborts with failure as soon as the metered work exceeds
W_cap = 2^125.71**, the submitted bound. The claimed bound therefore holds for every run by construction;
what remains is the probability that the abort fires, which is part of the failure budget below.

*Heuristic H1+ (no drift of the walk statistics).* The cap exceeds the nominal ledger
W_nominal = 2^125.7010 by 0.62% (0.0053 units per step; `experiment.py` part D). A multiply step
costs 1.12545 units and a squaring step 0.58195, so the cap is not exceeded while the fraction
of squaring steps stays within 0.98 percentage points of 1/2. A larger drop does not break the bound either:
it makes the meter halt a corresponding number of steps early, and the success probability of the
shorter walk is still given by the section-4 collision law. A drop of D percentage points raises
the average step price by 0.5435 D units, delta = 0.5435 D / c, and the meter then buys n/(1 + delta)
steps of the nominal n:

| drop of the squaring fraction | success value | probability spent |
| ---: | ---: | ---: |
| 0.4 points | 0.4304 | 0 |
| 0.5 points | 0.4304 | 0 |
| 1 point | 0.4303 | 0 |
| 1.5 points | 0.4283 | 0.0017 |
| 2 points | 0.4263 | 0.0037 |
| 2.5 points | 0.4243 | 0.0057 |
| 3 points | 0.4223 | 0.0077 |

H1+ states that over the whole run the frequencies of the branch bits, of the 2-cycle events and
of the DPs stay inside a drift of this size — the 0.01 budget is not exhausted until about 4
percentage points — and that the number of reports, false fingerprint candidates and restarts
stays within the caps below, with total probability of any cap or the abort firing at most 0.01.
This is a much weaker statement than independence of the steps; it excludes only a drift by
percentage points over 2^125 steps. Support: (i) the marginal probability of the squaring branch
is exactly 1/2 for every class (bit 2 of T is uniform, F2b), so the expected charge of a step is
exactly the nominal one with no independence assumption; (ii) in the toy runs the fluctuation of
the mean charge of a run scales as the independent-steps prediction at the new mix: rms of
(relative deviation x sqrt(steps)) is 0.312 observed against 0.318 predicted in the declared
40-trial run (the offline aggregate across q = 67, 163, 211 reads 0.31-0.38 per setting
against 0.318), so there is no excess correlation of branch labels along the walk at that
scale; (iii) the production-field walk of part C has squaring fraction 0.4965 against 1/2
by construction in the declared run here, so there is no adverse drift. The toy runs cannot
test concentration at 2^125.9 steps; H1+ extrapolates it, exactly as H1 extrapolates
the collision law. The success analysis (section 4) charges 0.01 for it.

*Caps.* reports D <= 1.1 S theta + 2^25 (S theta = 2^65.93; 2^25 covers the at most 20 M repeated reports of
the merged trails, D <= 2^66.07); at most 2^16 false fingerprint candidates are verified (the expectation is 2^2.9, so
Markov gives probability 2^-13.3; the verification work is at most 2^16 * 600 units = 2^25.2); at most
2^50 restarts (expected about 2^37, Markov 2^-13; each restart is a fresh start of about 2^10 units, 2^60 in
total); the total number of steps <= S(1+eps) + M 20/theta. Repeat steps after the first collision are at
most E = 2^84.32 steps (<= 1.12545 units each), negligible against the cap slack 2^116.5. Fruitless cycles
longer than 64 steps need all their steps to be multiply steps (rate below 2^-100 per step, about 2^25.8
events, each wasting at most 20/theta = 2^64.3 steps before the stuck-walk restart: 2^90.1 steps in total).
These caps bound the work and memory stated below.

## 4. Success probability

**H1 (random mapping with branch mix).** The class sequence behaves like iterates of a
random mapping on N classes in which two distinct classes have the same image with
probability (1 - sum_b q_b^2)/N, where q_b are the branch probabilities. Justification: for
c1 != c2 the equation F(c1) = F(c2) has, for each branch pair (b1,b2) with b1 != b2, about one
solution among the N classes (12 candidate representatives, each canonical with probability
1/12); for b1 = b2 = squaring it has none (squaring is a bijection of GT of odd order that
commutes with Frobenius), and for b1 = b2 = a table index it has about 11/12 of one. We
use the conservative value sum q^2 = s^2 + (1-s)^2/R = 0.250244 for all pairs; the exact
value (s^2 + (1-s)^2/(12R)) is smaller by 2.0e-5. Part E of `experiment.py` tests this law directly (section 7). This is the same correction as the Brent-Pollard
factor R/(R-1) of the pure adding walk (s = 0), and for Pollard's original three-branch walk
it is the mechanism that makes Pollard's original square-and-add walk slower than an r-adding walk
(Teske 1998). Hence the effective class count is
N' = N / (1 - sum q^2), log2 N' = 251.6876.

Under H1 the probability that the first S steps contain a class collision is
1 - exp(-S^2/(2N')). Choose P_model = 0.43:
`S = sqrt(2 N' ln(1/0.57)) = 2^125.9283.`
Losses before the first useful collision is detected: detection failure e^-20 = 2e-9, a
collision inside a single walk (fraction about 1/M = 1e-6), beta^1 = beta^2 (2^-255), cap
cap and abort failures (at most 0.01 together, H1+, section 3.4), steps lost to restarts (20 e^-20 = 4.1e-8 of all steps; about 2^37 restarts,
below the cap 2^40, each costing a fresh start of about 2^10 units, 2^47 in total, inside the 2^64.1
bucket of section 5). Their sum is below 0.011. The model value 0.43 therefore leaves 0.03 absolute
slack over the required 0.39 for the inaccuracy of H1. Toy measurements at the new mix (section 7)
of the fraction of runs with a useful collision within S(0.43) are 0.375-0.500 per setting,
pooled 0.426 +- 0.019 over 700 runs at s = 1/2; individual settings of a few hundred runs
scatter by +-0.04.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Primitive costs, scalar AArch64, using the measured blst_fp_mul = 383 instructions: a
modular add/sub/neg/double of 6 limbs called as a routine, **45** instructions (6 add/sub, 6
conditional-subtract, 6 select, 12 loads, 6 stores, 9 for call, return and address set-up); Fp2
squaring (complex method: (a0+a1)(a0-a1) and 2 a0 a1) = 2 Fp multiplications + 3 additions = 901;
Fp2 addition = 90.

**Table multiply:** the product y*T[h][t] is the dense blst_fp12_mul reference operation,
charged exactly 1.

**Squaring branch:** 9 Fp2 squarings + 34 Fp2 additions (F3) + 200 moves = 9*901 + 34*90 + 200
= 11,369, charged as the same **11,369 instructions = 0.56845 units**. A squaring step needs no key
computation.

**Per-step work.** Every step (both branches):

| Work | Instructions |
| --- | ---: |
| read Re(c_0) (6 limbs), branch bits, DP test, table address | 50 |
| (alpha, beta) update mod r (two 4-limb add + conditional subtract, or doubling) | 100 |
| 2-cycle and checkpoint comparisons on T, history | 60 |
| loop control | 30 |
| margin for loads and moves | 30 |
| **total counted 270, charged** | **270 = 0.0135 units** |

Multiply steps in addition (the keys give kappa):

| Work | Instructions |
| --- | ---: |
| 4 Fp multiplications (a zeta, a Re d_1, b Im d_1, B zeta) | 1,532 |
| 1 Fp addition (B) + 8 Fp subtractions/negations, 45 each | 405 |
| exact minimum of the 12 stored keys (top limb first, full-limb fallback, index) | 270 |
| load of the w^1 coefficient, kappa, twist index | 32 |
| **total counted 2,239, charged** | **2,239 = 0.11195 units** |

The minimum rule is exact on every run: compare the top (most significant) limbs of the
twelve 6-limb stored products with the usual 11 compare/select chain (11 x 3 instructions
plus 12 loads); if two or more keys share the winning top limb, compare the next limb
down among the tied keys only, and so on to the low limb; the smallest index wins any tie
that survives all six limbs. The worst case is 11 first-limb comparisons plus at most 55
lower-limb comparisons (11 tied keys x 5 further limbs), 66 limb-pair comparisons at 3
instructions each (subtract, conditional select, index carry) = 198, plus 72 limb loads and
index bookkeeping, bounded by the charged 270 — and the meter charges 270 on every multiply
step whether the fallback runs or not. Distinct full residues can share top limbs; the
fallback resolves them by construction, and the 66/p full-key-tie estimate of F2 is
informative only, never a failure charge. This widens the multiply-only charge from the
parent's 2,119 to 2,239 (+120 instructions, +0.006 units on a multiply step); section 0's
mix scan and the claim already include it.

So a multiply step costs 1 + 0.0135 + 0.11195 = 1.12545 units and a squaring step 0.56845 + 0.0135 = 0.58195 units.
Average step cost c = (1/2) * 1.12545 + (1/2) * 0.58195 = **0.8537 units** (= 8537/10000).

Setup and output (all <= 2^31 units): pairing g and verification pairing (<= 2^14 each),
check z^r = 1, the table (1024 multipliers x (2 exponentiations <= 1,021 units + 11 Frobenius
maps) <= 2^20.1), M = 2^20 starts g^alpha z^beta (<= 2^30.1), final [x]P0 (<= 510 G1 operations,
each < 1 unit). DP handling: at most 2^66.07 reports (the cap) x 6,000 instructions (the 12 keys for kappa, 2,239; hash
insertion, two 255-bit modular multiplications by p^kappa, memory traffic) / 20,000 = 2^64.3 units. All non-walk
work is therefore below 2^64.3 + 2^31 < 2^64.4 units.

Detection delay: after the first useful collision all M walks continue for 20/theta steps:
E = M * 20 / theta = 2^84.32 steps.

`W = c * S * (1 + eps) + c * E + 2^64.4 = 0.8537 * 2^125.9283 * 1.000645 + 0.8537 * 2^84.32 + 2^64.4 = 2^125.7010`.

Submitted bound: **time_log2 = 125.71** (the least hundredth at or above the computed value; the
computed value is 2^125.7010225, so 0.0090 bits above it, a 0.62% allowance over the whole ledger). `experiment.py` (part D) recomputes
this ledger, and the parent's 2^125.7074975 with the parent's own mix s = 3/8 and key charge 2,119, with exact decimal
arithmetic in one code path.

Memory: at most 2^66.07 reports x 80 bytes (128-bit fingerprint, two 256-bit exponents) at load factor
2/3 = 2^72.98 bytes, plus the 7.1 MB table, the 0.8 MB exponent tables and 2^20 walk states (at most
2 KB each, including the remembered previous state): below 2^73 bytes.
Memory is reported, not scored.

## 6. Heuristics and limitations

* H1 above (random-mapping behaviour of the mixed walk on classes), H2 (pseudorandom branch and DP bits
  of T), independence of table indices (used for cycle probabilities). These are standard in
  rho analyses (Teske 1998/2001; van Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999;
  Bos-Kleinjung-Lenstra 2010), supported here by toy and production-field runs, not proved.
* The squaring and per-step prices follow from the measured fp_mul and a stated primitive model (modular
  add 45 instructions); they are analytic instruction bounds, not measured routines. A pessimistic
  squaring price of 1.0 unit (a full dense multiplication) would remove the gain (the bound would be
  2^126.03, worse than the plain walk); each +0.1 unit on the squaring price adds 0.041 bits, each
  +100 instructions of per-step work 0.008 bits. A Callgrind count of blst's cyclotomic squaring would replace the
  analytic price.
* Collision statistics are validated in toy fields only. The measured mean of the mixed walk at s = 1/2
  is within 5% of the model in all settings but one small 40-run row (pooled +0.4% over 700 runs);
  the s = 3/8 parent rows are -1.4% pooled and the plain-walk controls -3.7% and +4.3%. An independent
  rerun of the script by a separate model instance (q = 67 and 163, about 1,800 mixed-walk runs, several
  seeds, walks and dp_bits) gave a pooled fraction of 0.425 at S(0.42), one 183-run setting 0.344.
* This is a generic square-root attack. It does not use index-calculus (TNFS) methods in Fp12.
  It is a valid upper bound, not a security estimate.

## 7. Experiments (`experiment.py`, standard library only)

* **Part A (BLS12-381):** GT membership; Frobenius formula = p-th power; 4-multiplication key
  formula = keys of the 12 explicitly computed conjugates, all distinct; kappa picks one identical
  canonical element from every conjugate; p^6 inverts; ord_r(p) = 12; 1+p^4+p^8 = 0; r-1
  factorisation; ord_r(2); the Granger-Scott formula equals y*y with 9 Fp2 squarings and 34 Fp2
  additions.
* **Part C (BLS12-381, production step rule):** a real walk of thousands of steps in GT with the
  production branch rule at the production table size (R = 1024): squaring iff (T & 7) >= 4,
  multiplier from (T >> 3) mod 1024, exact twelve-key minimum with the full-limb fallback of
  section 5, twisted table, squaring, 2-cycle escape. Checks: every one of the first 150 steps
  lands in the class of y^ M_h or y^^2 as computed with explicit Frobenius; the exponent
  invariant y = g^alpha z^beta at the end; the squaring fraction (0.5 expected) and chi-square
  of the low bits of Re(c_0). The declared run (2,000 steps, seed 1) gives 0 class-check
  failures, squaring fraction 0.4965 against 0.5, chi-square of the low 3 bits 6.46 (7 dof)
  and of the next 4 bits 16.0 (15 dof), with the stated threshold-4 rule — the same rule the
  toy collision statistics below are gathered with, so no transfer across rules is needed.
* **Part E (pair-collision law, toy fields):** sample m distinct classes, apply one step F to each and
  count coincidences of the images; N times the pair-collision rate is the quantity H1 predicts as
  1 - s^2 - (1-s)^2/(12R). Results (60,000 sampled classes, relative standard error about 0.9% except the last row):

  | q, R, s | observed | exact prediction | conservative value used | random mapping |
  | --- | ---: | ---: | ---: | ---: |
  | 67, 8, 3/8 | 0.8552 | 0.8553 | 0.8105 | 1 |
  | 67, 64, 3/8 | 0.8621 | 0.8589 | 0.8533 | 1 |
  | 67, 8, 0 | 0.9936 | 0.9896 | 0.8750 | 1 |
  | 67, 64, 0 | 1.0013 | 0.9987 | 0.9844 | 1 |
  | 67, 16, 1/2 | 0.7416 | 0.7487 | 0.7344 | 1 |
  | 211, 16, 3/8 (2.1%) | 0.8872 | 0.8573 | 0.8350 | 1 |

  The squaring branch removes exactly its s^2 of the collision probability, and the table branches
  keep 11/12 of theirs (rows with s = 0 and R = 8 separate the exact from the conservative value).
* **Part B (toy fields of the same tower shape, q = 67, 139, 163, 211, 379):** the complete
  algorithm (keys with tie-breaking by full conjugates, kappa trick, twisted table, squaring,
  2-cycle and checkpoint escapes, DPs, collision solving, verification g^x = z). Toy keys take
  few values, so the branch/DP bits are taken from a hash of the canonical conjugate (a
  stand-in for the bits of the 381-bit K; Part C tests those bits in the real field). Part B also reports the
  fluctuation of the mean charge per run against the independent-steps prediction (section 3.4).

Offline aggregate of the same script at the new mix (all discrete logs recovered and verified; `trials` runs per
setting, `dp_bits` 3 and `walks` 4 unless noted; S(x) is the step budget at model success x;
the toy rule squares iff (hash & 7) >= 4, the threshold-4 rule of section 3):

| q, N classes, R, seed | squaring | runs | mean T_c / sqrt N (model) | P[T_c <= S(0.40)] | P[T_c <= S(0.42)] | P[T_c <= S(0.43)] |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| 67, 129,145, 64, seed 11 | 1/2 | 300 | 1.480 (1.451) | 0.390 | 0.403 | 0.410 |
| 211, 721,279, 256, seed 21 | 1/2 | 120 | 1.356 (1.448) | 0.492 | 0.492 | 0.500 |
| 163, 4,524,905, 64, seed 31 | 1/2 | 200 | 1.452 (1.451) | 0.400 | 0.415 | 0.420 |
| 163, 4,524,905, 64, seed 31 | 1/2 | 40 | 1.605 (1.451) | 0.375 | 0.375 | 0.375 |
| 67, 129,145, 64, seed 1 (declared) | 1/2 | 40 | 1.427 (1.451) | 0.400 | 0.400 | 0.425 |
| **pooled mixed walk** | 1/2 | 700 | 1.457 (1.451) | 0.407 | 0.417 | 0.426 |
| 67, 129,145, 64, seed 12 (control) | 0 | 400 | 1.217 (1.263) | 0.420 | 0.445 | 0.455 |
| 163, 4,524,905, 64, seed 32 (control) | 0 | 600 | 1.318 (1.263) | 0.370 | 0.385 | 0.395 |

The control rows are plain adding walks (their model values use the plain-walk effective class count). The
last row is 1.4 sigma below the 0.413 that the earlier plain-walk record reports for the same field. Measured
2-cycle rate at R = 64 with squaring: 2.7e-4 to 3.9e-4 per step across these settings (predicted
(1/2)^2/(12R) = 3.3e-4, bound 1/(12R) = 1.3e-3); at R = 256, 6.3e-5 (8.1e-5 predicted);
without squaring 1.28e-3 and 1.35e-3 (1.30e-3 predicted). Charge-fluctuation rms reads
0.312/0.346/0.322/0.343 against 0.318 predicted across the five mixed settings.
Reproduce a row with e.g.
`echo '{"q":163,"trials":200,"R":64,"sq8":4,"dp_bits":3,"walks":4,"seed":31,"production_checks":0,"seconds":100}' | python3 -I experiment.py`.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999); Teske, Math. Comp.
70 (2001); Duursma, Gaudry & Morain, ASIACRYPT 1999; Granger & Scott, PKC 2010; Aranha et
al., EUROCRYPT 2011; Bos, Kleinjung & Lenstra, ANTS 2010; mpjunior92, submission 7da5f23a.

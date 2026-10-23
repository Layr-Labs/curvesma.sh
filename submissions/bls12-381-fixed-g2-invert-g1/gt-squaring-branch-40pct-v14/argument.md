# Argument — gt-squaring-branch-40pct-v14

Claim **125.73**; computed **W = 2^125.725..125.726**; success >= 0.39; C = 20,000 instructions per dense Fp12 multiplication.

## 1. Reduction (unchanged from the accepted construction)

Q != 1 in G2 generates G2 (prime order r); the pairing is nondegenerate, so g = e(P0,Q) generates GT.
For z = g^x the algorithm returns P = [x]P0 and verifies e(P,Q) = z; z = 1 returns infinity.

## 2. Class walk and the free invariant

p^6 = -1 mod r, so the p-power Frobenius acts on GT with order 12 and inversion in the orbit: GT is
partitioned into N = (r-1)/12 classes of twelve conjugates. Writing y = c0 + c1 w in
Fp12 = Fp6[w]/(w^2-v), Fp6 = Fp2[v]/(v^3-xi), xi = 1+u, the w^0 coefficient of y^(p^k) is conj^k(c0) and
the w^1 coefficient is c1 zeta^j (k = 2j) or conj(c1) gamma zeta^j (k = 2j+1), with
zeta = |gamma|^2 = xi^((p^2-1)/6) a primitive sixth root of unity. Hence

* **T = Re(c0) is invariant across the twelve conjugates** and costs no multiplication (one field-word
  read). Branch, multiplier index and distinguished-point bits are read from disjoint bit ranges of T.
* the twelve keys are +-{A zeta^j, B zeta^j}, A = Re(c1), B = Re(conj(c1) gamma): four Fp
  multiplications (B: two constant multiplies; A zeta and B zeta) give the six representatives and six
  subtractions give the rest, because zeta^2 = zeta - 1.

**+-pair lemma.** zeta^3 = -1, so key_(k+6) = zeta^3 key_k = -key_k. Hence min over twelve keys =
min(min6, p - max6), and the argmin index is argmin6 or argmin6 + 6. Verified in
`certificate.json` (`plus_minus_pair_lemma`, `six_key_min_equals_twelve_key_min`,
`trace_invariant_across_twelve`, `kappa_canonicalises_all_conjugates`).

## 3. Step rule

* branch (probability 3/8, `T mod 8 < 3`): y <- y^2 by cyclotomic squaring, (alpha,beta) <- (2alpha,2beta).
  A squaring step needs no twist and no key, so the key computation is skipped on these steps.
* otherwise kappa = argmin index, y <- y * M_h^(p^(-kappa)), (alpha,beta) <- (alpha,beta) + lambda^(-kappa)(a_h,b_h).

Both maps commute with Frobenius, so the rule is well defined on classes, and T being a class function
makes the branch, multiplier and DP bits agree on all twelve conjugates.

## 4. Instruction rows (counted by running `scripts/gt_prod.py`)

| row | instructions |
| --- | ---: |
| every step: trace read, branch/multiplier/DP bits, dispatch | 60 |
| every step: exponent pair update mod r | 100 |
| every step: cycle/checkpoint compares, loop control, margin | 140 |
| multiply: six +-pair keys (4 Fp mul + 6 adds, counted) | 1802 |
| multiply: min of six, +-pair index bookkeeping | 180 |
| multiply: table address and load | 150 |
| multiply: exponent pair add | 60 |
| **multiply step** | **2492 = 1.12460 units** |
| squaring: cyclotomic squaring (24 Fp mul + 130 adds, counted) | 15042 |
| squaring: exponent pair doubling | 40 |
| **squaring step** | **15382 = 0.76910 units** |

Prices are the measured blst convention used on this track: 383 instructions per Fp multiplication and
45 per modular add, both as called routines; the dense multiplication is the reference operation (1.0
unit). No price is taken from another submission: every row is the count of the code in
`scripts/gt_prod.py` under those two prices. Our own naive dense multiplication counts 38,466
instructions against the 20,000 convention, so the non-reference rows are priced on a scale about 1.9x
stricter than the reference operation — the error is in the conservative direction.

## 5. Collision law (independent derivation)

For two states in distinct classes a merge happens iff the two images coincide as classes.

* both squaring (probability s^2): y1^2 = y2^2 gives (y1/y2)^2 = 1; GT has odd prime order, so -1 is not
  in GT and y1/y2 = 1. **No merge**: squaring is injective on classes.
* both multiplying with the same multiplier index (probability (1-s)^2/R): the same bijection on classes
  is applied twice, again no merge.
* every other combination (probability 1 - s^2 - (1-s)^2/R): the images coincide with probability 1/N.

So the merge probability per pair is rho/N with rho = (1-s^2)(1-1/R) and N' = N/((1-s^2)(1-1/R)).
At s = 3/8, R = 1024: N = 2^251.272, N' = 2^251.492 — the mixture **costs** 0.11 bits of step count
and the branch only pays through the step price.

## 6. Bound

S = ceil(sqrt(2 N' ln(100/61))) = 2^125; c_avg = 5/8 c_mul + 3/8 c_sqr = 0.99129;
W = c_avg (S + E) + 2^65 with E = 2^85 detection steps. W = 2^125.725..125.726 < 125.73. Gate seen at
submission time: 125.74 (747ad6c9, review; other live entries 125.80, 125.826, 125.89 x3).

## 7. Production-field execution

`scripts/gt_prod.py` runs the step rule on real BLS12-381 GT elements (random generator of the order-r
subgroup, real table of 1024 entries, 20-bit DP): 2000 steps, squaring fraction
0.3835 (95% half-width 0.02122; target 0.375), 63
exponent-invariant checks with 0 failures, 33 class-invariance
tests against all eleven explicit Frobenius conjugates with 33 agreements, chi-square
5.1 on 7 degrees of freedom of T mod 8, and 23/23 algebra checks.

## 8. Toy

`scripts/toy.py` (q = 523, N = 43420, R = 128, branch 3/8, fresh target per trial,
invariant asserted at start, every 16th step, at escapes and at the collision point):
73/180 = 40.556 first-report recovery at S_ops
(model success 0.421), 5766 checks, 0 failures, squaring fraction
20181/52309. Consistency evidence only: a small-field stdlib toy
cannot resolve production-scale mixing, so the branch and collision claims rest on Sections 5 and 7.

## 9. Limitations

H1 (random-mapping behaviour of the mixed class walk) and H2 (pseudorandomness of the disjoint bit
ranges of T) are heuristic and unproved, as is the independence of table indices for the 2-cycle rate.
The squaring price is analytic from counted Fp2 operations, not a measured routine. No production-size
run; this is a square-root attack, not a security estimate.

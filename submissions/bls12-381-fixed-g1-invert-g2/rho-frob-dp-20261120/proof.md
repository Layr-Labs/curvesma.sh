# Fixed-G1 pairing inversion into G2: one GT logarithm by a Frobenius-class distinguished-point rho walk

Problem: given the instance's nonidentity P in G1 and z in GT, return Q in G2
with e(P,Q)=z. Here e is the protected function
`py_ecc.optimized_bls12_381.pairing(Q,P)`, which takes G2 first. We use the
admissible route named in rules.md: solve z=e(P,Q0)^x in GT and return [x]Q0.

This is a generic Pollard-rho attack using the standard Frobenius speed-up:
the walk runs on classes {Y^(p^i) : i=0..11}, each of size 12. It is not a
best-known finite-field attack, a novelty claim, a practical break or a
security-bit estimate. The work bound is analytical and capped. Success
relies on explicitly disclosed random-walk heuristics. We seek the protected
exploratory `plausible_not_refuted` standard and assume no prior acceptance.

## Ingredients

* **One reference operation per walk step.** Each step multiplies by the
  pinned blst dense Fp12 multiplication `blst_fp12_mul`, in blst's tower
  representation. The protected cost model prices this as one unit. All other
  per-step work is charged as at most E = 12,288 scalar instructions
  (machine.md), i.e. 0.6144 units.
* **Frobenius classes.** Frobenius pi(Y)=Y^p maps GT to itself. Since p has
  order 12 modulo r, it acts as exponentiation by lam = p mod r, an element
  of order 12 in (Z/r)^*; lam^6 = -1 (checked in the experiment).
  - A GT element Y != 1 has exactly 12 distinct conjugates, because r divides
    no p^k-1 for k < 12.
  - The walk is on the r/12 classes, so the expected collision time drops by
    sqrt(12).
  - In the blst tower, pi^i acts monomial-wise: the coefficient c of v^k w^e
    maps to sigma^i(c)*gamma[i][n], where n=2k+e, sigma is conjugation of Fp2
    and gamma[i][n]=xi^(n(p^i-1)/6) with xi=1+u. Since p = 7 mod 12, 6
    divides p^i-1.
* **Free class key.** gamma[i][0]=1, and sigma fixes real parts, so the real
  part of the constant coefficient c_{0,0} is the same for all 12 conjugates.
  key(Y) = its canonical Montgomery residue t0 (< 2^381). The walk index is
  j = (word 0 of t0) & 31. Y is a distinguished point (DP) iff bits 64..173
  of t0 are zero (D=110 bits). The DP fingerprint is bits 192..319 of t0.
* **Equivariant frame** (`production_canon` in rho.py). The steps are:
  1. For i=0..5 compute V_i = Mont(Re(sigma^i(c_{1,0}) * gamma[i][2])), the
     real part of the v coefficient of pi^i(Y). The constants gamma[i][2] are
     1, imaginary, real, u, real, imaginary, so V_i is one Fp multiplication
     (or a copy).
  2. pi^6 fixes this coefficient, so let i0 = argmin V_i.
  3. pi^6 negates the w coefficient, so compute
     W = Mont(Re(sigma^i0(c_{0,1}) * gamma[i0][1])).
  4. Set i* = i0 if 1 <= W <= (p-1)/2, else i* = i0+6.

  The canonical member pi^{i*}(Y) is therefore the conjugate whose v-real
  part is least, among the two with w-real part in [1,(p-1)/2]. It is the
  same element whichever conjugate Y we start from. The experiment checks
  this on all 12 conjugates of random full-size GT elements.

  Ties (equal V_i for different i, or W=0) need a coefficient on a
  hyperplane, probability about 2^-380 per step. The V_i are distinct linear
  forms (a, B*b, C*a, b, ...) in the two Fp coordinates.

## Complete algorithm

1. Read the profile and the instance (P, z). Check P on the G1 curve,
   P != infinity, z != 0 and z^r = 1. If z = 1 return infinity, which is
   valid exactly for z = 1.
2. Compute g = e(P,Q0) once by evaluating the protected pairing on the actual
   P (cost in machine.md). g has order r: the pairing is nondegenerate and P
   and Q0 are nonidentity elements of prime-order groups.
3. Convert g and z into the blst tower with canonical Montgomery coefficients.
   The basis map is `poly_to_tower`; the tower uses the protected w, with
   u = w^6-1.
4. Sample 66 uniform scalars mod r by bounded rejection. Build M_j = g^c_j z^d_j
   for j<32 and a start X = g^a z^b. Precompute the twisted table
   T[i][j] = pi^{-i}(M_j), 384 records, with label increments
   (lam^{-i} c_j, lam^{-i} d_j) mod r.
5. Walk. The state is (X, a, b) with X = g^a z^b. Each step:
   - Compute key t0 and frame i* from X; the canonical member is
     C = pi^{i*}(X).
   - **Short-cycle check:** compare key against the last 64 states in a ring
     buffer. On a match with entry (X', a', b', i'), compare canonical
     labels: (lam^{i*}a, lam^{i*}b) against (lam^{i'}a', lam^{i'}b').
     - If they differ, this is a genuine relation: solve for x and verify
       g^x = z.
     - If they are equal, it is a fruitless cycle. Take the ring entry with
       the least key in the cycle, set X <- X_min^2 (one blst_fp12_mul) and
       (a,b) <- 2(a_min,b_min), clear the ring, and count one step. Squaring
       commutes with pi, so the new class depends only on the cycle's set of
       classes.
   - **DP:** if X is a DP, look up its fingerprint in a table of at most 2^21
     entries, storing canonical labels (lam^{i*}a, lam^{i*}b).
     - If the fingerprint is absent, store it.
     - If it is present with different canonical labels (A,B): if B = b
       report failure, else x = (A-a)/(b-B) mod r and test g^x = z.
     - If the test fails (a fingerprint false match), replace the entry;
       after 64 failed tests report failure.
   - **Step:** X <- blst_fp12_mul(X, T[i*][j]), a += lam^{-i*} c_j,
     b += lam^{-i*} d_j.

   This takes the class of X to the class of C*M_j, since
   pi^{-i*}(C*M_j) = X*T[i*][j]. Stop with failure after
   S = floor(31*(isqrt(r/12)+1)/20) + 2^115 steps, counting escapes.
6. Verify g^x = z and return Q = [x]Q0 (affine Fp2 double-and-add).

`rho.py` (`frob_rho`, `production_canon`, `tower_frob`) is the executable
reference with the exact loop order. There is one trial, with no hidden
restart, amplification or amortization.

## Correctness of the output

The walk invariant X = g^a z^b holds after every step, because
log pi^{-i}(M_j) = lam^{-i} log M_j, and after every escape (squaring
doubles both labels). The canonical member pi^{i*}(X) has labels
(lam^{i*}a, lam^{i*}b).

Two visits to the same class give two label pairs for the same element C:
a1 + x b1 = a2 + x b2 mod r.
- Equal pairs (a fruitless cycle) carry no information and are never used
  to solve.
- If b1 != b2, then x = (a2-a1)/(b1-b2), and verification g^x = z rules out
  any wrong output, including fingerprint false matches.
- b1 = b2 with a1 != a2 is impossible for a genuine equality.

Then e(P,[x]Q0) = g^x = z by bilinearity. [x]Q0 is on the G2 curve with
[r]Q = infinity, and the instance's own P is used.

## Detection, success heuristic, support and limits

**Model.** Write F for the induced deterministic map on classes: the walk
step plus the deterministic escape rule. Escapes depend only on the classes
in the detected cycle, so F is a function of the class and the bounded
recent history, and that history also repeats once the walk repeats. Let
mu and lambda be the tail and cycle length of the class sequence, and delta
the distance from the first repeated class to the next DP. As in the
single-element DP analysis, the algorithm reaches a genuine verification by
step mu+lambda+delta when

    (i) mu+lambda <= 1.5*sqrt(r/12),
    (ii) delta <= 2^115 < lambda,
    (iii) the relation from the long cycle is not identically fruitless,

and escape steps use at most 0.05*sqrt(r/12) of the cap.

**(i) Collision time.** For a random map on r/12 classes, the r-adding
heuristic with 32 increments gives

    P(i) ~= 1 - exp(-(31/32)*1.125) ~= 0.664.

**(ii) DP gap.** As in plain DP rho: about exp(-32), plus about 2^-12 for a
short cycle.

**Fruitless cycles.** The only 2-cycle mechanism is the same index twice
with frame shift 6 (pi^6 is inversion on GT). Its rate is about
1/32 * 1/12 = 1/384 per step. Each costs at most 3 steps (two visits and
one squaring), so the overhead is about 0.8% of steps, well inside the 3.3%
slack (1.55 vs 1.5).
- Longer fruitless cycles need several coincident index and frame
  conditions, with rate at most about 32^-2 per step for length 4. They are
  caught by the 64-entry ring.
- Fruitless cycles longer than 64 need on the order of 32 coincidences, with
  probability far below 2^-128 per step.
- A fruitless relation from the long rho cycle needs all 32 coefficient sums
  sum lam^{f_t} over the visits to vanish at once. We treat this as
  negligible: a disclosed heuristic.

**Other losses.** Table overflow (mean S*2^-110 ~= 2^16 DPs against 2^21
slots), fingerprint false matches (< 2^-80) and rejection sampling are
negligible. Degeneracy (b1 = b2 at a genuine collision) has probability
1/r. Conditional on the increment logarithms s_j and the start logarithm,
the d_j are uniform. The class walk, frames, escapes and first relation
depend only on the s_j. The coefficient of the relation in the d_j is a
fixed vector, nonzero for a genuine relation, so b1-b2 is uniform.

**Declared success: >= 0.50** (prediction about 0.66; the minimum is 0.39).

**Support.** The class-walk method with automorphisms of order m and
fruitless-cycle handling follows standard literature:
- Wiener and Zuccherato (SAC 1998) and Gallant, Lambert and Vanstone
  (Math. Comp. 2000): speed-up by sqrt(m).
- Duursma, Gaudry and Morain (ASIACRYPT 1999): Frobenius classes and
  fruitless cycles.
- Bos, Kleinjung and Lenstra, "On the use of the negation map in Pollard
  rho" (ANTS 2010): cycle detection with escape by doubling the minimal
  element.
- Teske (Math. Comp. 2001): r-adding walks.

These are context; the reasoning used is given here.

**Evidence** (`experiment.py`; development output in
`development-evidence.json`):

1. *Toy class walks.* The exact `frob_rho` loop ran on toy groups Z_q with
   an order-12 automorphism, using SHA-modelled equivariant key and frame,
   DP rate 2^-2 and 128 trials per q, for q about 1.4e5, 1.6e6 and 1.8e7.
   - Fraction succeeding within 1.5*sqrt(q/12) steps, escapes included and
     without the gap: 0.66, 0.71, 0.76 (prediction 0.664).
   - Mean steps 1.24, 1.17, 1.08 times sqrt(q/12).
   - Fruitless escapes 0.0021-0.0028 per step (prediction 1/384 = 0.0026).
   - No degenerate collision, table overflow or verification failure.
2. *Basis check.* The tower product equals the protected-basis product.
3. *Frobenius checks at full size.* The Frobenius formula matches Y^p on
   random Fp12 elements, and pi(g) = g^(p mod r) on a full-size order-r
   element. The class key is invariant and the frame canonical member
   identical across all 12 conjugates of 4 random GT elements.
4. *Real full-size class walk.* A 1024-step walk with the production key,
   frame and twisted table used every index (chi-square 26 on 31 degrees of
   freedom) and all 12 frames roughly uniformly. It made 2 fruitless escapes
   (about 2.7 expected), and the label invariant held.

Without escapes, the same full-size walk is trapped in a 2-cycle (observed
during development). The escape rule is necessary, and it works.

**Limits.** These checks are evidence at their scale only; extrapolation to
2^126-step walks is part of the disclosed heuristic. The experiment
evaluates neither the py_ecc pairing nor blst itself. Finite-field
(Sex)TNFS may be stronger; nothing here claims tightness. The work and
memory caps hold even if the heuristics fail.

## Work, memory and score

    S = floor(31*(isqrt(floor(r/12))+1)/20) + 2^115
    W = S*(1 + 12288/20000) + 2^40/20000           (reference units)
    submitted score = ceil(10^6*log2(W))/10^6 = 126.959915

The 2^40 instructions cover the pairing, conversions, table and Frobenius
twisting, randomness, DP handling, verifications, output and checks
(machine.md). Peak memory is below 2^28 bytes. Advice is zero.
Preprocessing is one-target only and charged in full.

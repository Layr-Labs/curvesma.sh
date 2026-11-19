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
  j = (word 0 of t0) & 255 (NINC = 256 increments). Y is a distinguished point (DP) iff bits 64..173
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
4. Sample 514 uniform scalars mod r by bounded rejection. Build M_j = g^c_j z^d_j
   for j<256 and a start X = g^a z^b. Precompute the twisted table
   T[i][j] = pi^{-i}(M_j), 3072 records, with label increments
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
       the least key in the repeated stretch, set X <- X_min^2 (one
       blst_fp12_mul) and (a,b) <- 2(a_min,b_min), and count one step.
       Squaring commutes with pi, so the new class depends only on the
       cycle's classes. The ring is **not** cleared: if the walk ever returns
       to those classes, the old labels differ from the new ones by the
       doubling, so the return is a genuine relation, not a second fruitless
       detection.
   - **DP:** if X is a DP, count a DP event (more than 2^22 events is a
     reported failure) and look up its fingerprint in a table of at most 2^21
     entries, storing canonical labels (lam^{i*}a, lam^{i*}b).
     - If the fingerprint is absent, store it.
     - If it is present with **equal** canonical labels, the walk is on a
       fruitless cycle through this DP: escape by squaring this element
       (as above).
     - If it is present with different canonical labels (A,B): if B = b
       report failure, else x = (A-a)/(b-B) mod r and test g^x = z.
     - If the test fails (a fingerprint false match), replace the entry;
       after 64 failed tests (ring and DP together) report failure.
   - More than S/256 escapes is a reported failure.
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

**Model.** The analysis uses one explicit heuristic model, the **random
model**. Whenever the walk reaches a class not seen before, the index j of
that step (8 bits of t0) and the frame i* are uniform and independent of
everything earlier: j on 256 values, i* on 12. All probability statements
below are in this model. The model is a disclosed heuristic, not a theorem
about GT, blst or the coefficient bits.

Let mu and lambda be the tail and cycle length of the walk's class sequence,
and delta the distance from the first repeated class to the next DP. The
algorithm reaches a genuine, verified relation within the cap when all of
the following hold:

    (i)   mu+lambda <= 1.5*sqrt(r/12);
    (ii)  delta <= 2^115 < lambda;
    (iii) the walk is never trapped in an undetected fruitless cycle;
    (iv)  escapes use at most 0.05*sqrt(r/12) steps and stay below S/256,
          and the DP-event, table and verification limits are not hit.

**(i) Collision time.** For a random map on r/12 classes, the r-adding
heuristic with 256 increments gives

    P(i) ~= 1 - exp(-(255/256)*1.125) ~= 0.674.

**(ii) DP gap.** The probability is about exp(-32), plus about 2^-12 for a
short rho cycle.

**(iii) Fruitless cycles: what can trap the walk.** A class repeat whose two
canonical label pairs are equal gives no information. Those label pairs are
the coefficients of g and z in the canonical member. Every other class
repeat gives a genuine relation, which is solved and verified, so the walk
cannot be trapped by it. There are three cases:

- **Repeats within 64 steps.** The ring catches these, whatever their cause.
- **Repeats on a cycle through a DP.** The second visit to that DP has equal
  stored labels, which triggers an escape.
- **Escape re-entry.** An escape changes the labels by a factor of 2. The
  ring keeps the earlier labels, and DP labels are stored, so returning to
  an escaped cycle gives a genuine relation (or the 1/r degenerate case),
  not a second fruitless detection.

So the only trap is a fruitless cycle with **length l > 64 and no DP on
it**.

**(iii) Fruitless cycles: bounding the trap.** Everything here is in F_r;
no algebraic-integer identity is used.

*Labels.* Write every label as a linear form, with coefficients in F_r, in
the independent uniform coins: the start labels (a0,b0) and the increments
(c_j,d_j). Each step multiplies the stored member by pi^(-i_u)(M_j), so it
adds lam^(-i_u)*(c_j,d_j) to (a,b). An escape replaces (a,b) by twice the
labels of an earlier ring or DP state. Canonical labels multiply by
lam^(i*).

The coefficient of b0 in the canonical b-label at any visit is
lam^f * 2^e, where f is a frame exponent and e is the number of escapes on
the path from the start. This is never 0.

*Conditioning (exact, for the actual algorithm).* Let x be the true
logarithm, and let s_j = c_j + x*d_j and s0 = a0 + x*b0. The group walk,
and so every class, frame, index, escape, DP decision and ring decision,
depends only on (s_j, s0) and the DP/partition bits of the elements.
Conditional on all of these, (d_j) and b0 remain independent and uniform:
for each value of d_j there is exactly one c_j matching s_j, and likewise
b0 with a0.

So for any two visits (t, t+l) to the same class, the difference of their
canonical b-labels is a linear form Delta = sum_j E_j*d_j + B*b0, whose
coefficients E_j and B are fixed by the walk. A trap needs Delta = 0
(equal labels).

*Case N: nonzero coefficients.* If (E_j, B) is not the zero vector, Delta is
uniform on F_r given the walk, so P(Delta = 0) = 1/r.

*Case Z: all coefficients zero.* Then B = 0, which needs
lam^(f1) 2^(e1) = lam^(f2) 2^(e2) in F_r.
- With e2 > e1 this means 2^(12(e2-e1)) = 1 mod r. But e2-e1 <= S/256 <
  2^119, and the multiplicative order of 2 modulo r exceeds 12*2^119.
  experiment.py computes this order exactly from the factorization of r-1,
  each factor verified prime. So e1 = e2: no escape inside the cycle.
- Then f1 = f2, i.e. no frame shift at closure.
- With no escape inside, E_j = lam^(f1) * sum over the steps u in [t,t+l)
  with j_u = j of lam^(-i_u). The values lam^(-i), for i = 0..11, are 12
  distinct elements of F_r, since lam has order exactly 12 (checked).

Case Z therefore needs, for every index j used in the cycle, the
F_r-sum of its terms lam^(-i_u) to be 0. An index used once gives a single
nonzero term. For an index used at least twice, in the random model (frames
of fresh steps uniform on 12 and independent), the last term must equal one
specific element, so the probability is <= 1/12.

*The bound.* With n = 256 increments, if the cycle's l fresh steps use
exactly k distinct indices, each at least twice, then

    P_Z(l) <= sum_{k=1}^{min(n,l/2)} C(n,k) (k/n)^l 12^(-k).

A trap also needs no DP on the cycle (otherwise the DP revisit escapes) and
l > 64 (otherwise the ring sees it). Summing over lengths with the no-DP
factor (1-2^-110)^l:

- **Case Z per start:** at most 2^-173.2 (`trap_bound_log2`).
- **Case N per start:** at most
  sum_l (1-2^-110)^l / r <= 2^111/r = 2^-143.8.

Over all S < 2^126.3 starting steps:

    P(trap) <= 2^(126.3-173.2) + 2^(126.3-143.8) <= 2^-17.4.

The experiment prints both terms and their sum. This bound uses the random
model only in case Z (frames uniform) and in the no-DP factor (DP bits
random). The 1/r of case N is exact given the walk.

**(iv) Escapes and limits.** The 2-cycle (same index twice with frame shift
6, since pi^6 is inversion) occurs at about 1/(256*12) = 1/3072 per step.
Each costs at most 3 steps, so the overhead is about 0.1%, against a slack of
3.3% (1.55 vs 1.5) and the escape limit S/256.
- **DP events:** about S*2^-110 ~= 2^16.3 are expected, against the 2^22
  limit; the Chernoff bound gives < 2^-1000.
- **DP table:** 2^21 slots; overflow is likewise < 2^-1000.
- **Fingerprint false matches:** < 2^-80.
- **Rejection sampling:** < 514*2^-1024.

Degeneracy (b1 = b2 at a genuine relation) has probability exactly 1/r.
Conditional on the increment logarithms s_j and the start logarithm, the d_j
are uniform. The class walk, frames, escapes and first relation depend only
on the s_j. The relation's coefficient vector in the d_j is fixed and
nonzero for a genuine relation, so b1-b2 is uniform.

**Declared success: >= 0.50.** The prediction is about 0.67 minus less than
2^-12; the protected minimum is 0.39. The work cap is unconditional: every
limit above ends the run, and all of the run's work is charged
(machine.md).

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
`development-evidence.json`; seed 20261122):

1. *Toy class walks.* The exact `frob_rho` loop, with 256 increments, ring,
   DP-escape and all limits, ran on toy groups Z_q with an order-12
   automorphism. Key and frame were SHA-modelled and equivariant; DP rate
   2^-2; 128 trials per q, for q = 174673, 1483021 and 11363749.
   - Fraction succeeding within 1.5*sqrt(q/12) steps, escapes included and
     without the gap: 0.66, 0.77, 0.59 (prediction 0.674).
   - Mean steps 1.23, 1.03, 1.22 times sqrt(q/12).
   - Fruitless escapes 0.00021-0.00039 per step (prediction
     1/3072 = 0.00033).
   - Every failure was the cap. No degenerate collision, limit hit, table
     overflow or verification failure occurred.
2. *Basis check.* The tower product equals the protected-basis product.
3. *Frobenius checks at full size.* The Frobenius formula matches Y^p on
   random Fp12 elements, and pi(g) = g^(p mod r) on a full-size order-r
   element. The class key is invariant and the frame canonical member
   identical across all 12 conjugates of 4 random GT elements.
4. *Real full-size class walk.* A 1024-step walk with the production key,
   frame and twisted 3072-record table had index chi-square 279.5 on 255
   degrees of freedom, and all 12 frames were used roughly uniformly
   (70-100 each). It made 0 escapes (0.33 expected), and the label
   invariant held.
5. *Trapping bound and the order of 2.* The experiment verifies the
   factorization of r-1 (each factor prime by Miller-Rabin with fixed bases,
   which is deterministic below 3.3*10^24; all factors are smaller) and
   computes ord_r(2) exactly. It then evaluates both terms of (iii) and
   their sum: case Z 2^-46.9, case N 2^-17.5, total <= 2^-17.4.

A previous 32-increment version of this walk, without escapes, was observed
to be trapped in a 2-cycle at full size. Escapes are necessary. The finite
experiments cannot probe the rare-event regime of (iii); that regime is
covered by the analytical bound in the stated random model.

**Limits.** These checks are evidence at their scale only; extrapolation to
2^126-step walks is part of the disclosed heuristic. The experiment
evaluates neither the py_ecc pairing nor blst itself. Finite-field
(Sex)TNFS may be stronger; nothing here claims tightness. The work and
memory caps hold even if the heuristics fail.

## Work, memory and score

    S = floor(31*(isqrt(floor(r/12))+1)/20) + 2^115
    W = S*(1 + 12288/20000) + 2^42/20000           (reference units)
    submitted score = ceil(10^6*log2(W))/10^6 = 126.959915

The 2^42 instructions cover the pairing, conversions, table and Frobenius
twisting, randomness, DP handling, verifications, output and checks
(machine.md). Peak memory is below 2^28 bytes. Advice is zero.
Preprocessing is one-target only and charged in full.

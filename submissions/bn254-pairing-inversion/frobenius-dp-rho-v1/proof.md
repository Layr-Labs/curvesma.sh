# Frobenius-class distinguished-point rho for BN254 GT, then (P0,[x]Q0)

This is a generic square-root attack on the actual protected 254-bit
prime-order GT, improved by (i) the order-12 Frobenius automorphism group of
GT, (ii) distinguished points instead of Floyd's three steps per iteration,
and (iii) a tower-field instruction count. It is not an index-calculus/NFS
attack, not a practical break and not a claim of optimality. The work bound
is an analytical hard cap; the success bound uses an explicitly disclosed
random-walk heuristic and seeks the protected exploratory
`plausible_not_refuted` standard. No prior result is assumed accepted.

Notation: p, r are the protected field modulus and subgroup order, g the
protected generator, z=g^x the input. All logarithms of GT elements are
taken to base g and live in Z/r. rho := p mod r.

## 1. Algebraic facts (all checked at full width by `experiment.py`)

F1. GT is the order-r subgroup of Fp12^*, so the p-power Frobenius phi maps
GT to itself and acts on logarithms as multiplication by rho:
phi(g^y)=g^(rho*y). rho has exact order 12 modulo r and rho^6=-1 mod r
(r divides Phi_12(p)=p^4-p^2+1). For 0<k<12, rho^k-1 is nonzero modulo the
prime r, hence invertible, so every nonidentity element of GT has exactly 12
distinct conjugates phi^k(Y). The nonidentity elements split into exactly
N=(r-1)/12 classes of size 12. Inversion is phi^6, so it is included.

F2. For the 11 values j in J={0,1,2,3,4,5,7,8,9,10,11} put
c_j = 2+rho^j mod r. Define the step map

    f_j(Y) = Y^2 * phi^j(Y),   log f_j(Y) = c_j * log Y  (mod r).

phi commutes with squaring, multiplication and phi^j, so
f_j(phi^k(Y)) = phi^k(f_j(Y)): each f_j maps classes to classes. Every c_j
is nonzero mod r (asserted), so f_j never produces the identity from a
nonidentity element. (j=6 is excluded: 2+rho^6=1 would be a fixed point.)

F3. Write Fp12 = Fp2[w]/(w^6-xi), xi=9+u, u^2=-1 (the protected basis with
u=w^6-9). For Y = sum_{i<6} a_i w^i with a_i in Fp2 put t(Y)=Re(a_0) in Fp;
in the protected coefficient list, t = c0+9*c6. Then 12*t(Y) is the trace
Tr_{Fp12/Fp}(Y), since w^i has zero trace over Fp2 unless 6 | i. A trace is
Frobenius-invariant, so t is constant on each class. t costs nothing: it is
one stored coefficient. t need not be injective on classes; nothing below
uses injectivity for correctness.

F4. No structural short cycles. The class walk multiplies logarithms by the
images of the algebraic integers 2+zeta^j in Z[zeta_12]. Their absolute
norms are 81 (j=0), 13 (j=1,5,7,11), 49 (j=2,10), 25 (j=3,9), 9 (j=4,8);
none is a unit, so no product of step multipliers is identically a unit
rho^k: there are no fixed points or 2-cycles forced by the algebra (unlike
the cheaper map Y*phi^j(Y), whose multipliers 1+zeta^j are units for
j=1,4,5,7,8,11; we tested that variant on toy groups, it stalls, and we do
not use it). The only identity among the multipliers is
c_4*c_8 = c_0 (=3); it merges two paths but creates no cycle. It is present
identically in the toy experiments, whose class walk has the same
Z[zeta_12] structure.

## 2. Complete algorithm

Parameters: theta=2^-110; walk-length cap L=2^116; walk cap 2^20; total
step cap S=floor(sqrt(3(r-1)/16))+1 >= 1.5*sqrt(N).

1. Decode z (12 canonical coefficients), check z!=0 and z^r=1, convert to
   the tower basis and Montgomery form. If z=1 output P0 and infinity.
2. Draw a fresh 128-bit salt. Precompute the Frobenius constants
   gamma_{j,i}=xi^(i(p^j-1)/6) in Fp2 and c_j (charged, section 4).
3. Repeat for walk number n=1,2,...,2^20:
   a. Draw independent uniform a,b in Z/r (bounded rejection sampling) and
      set Y=g^a*z^b, C=1. If Y=1 skip this walk.
   b. Loop: h=SHA256(salt || enc(t(Y))) (one block). If the 110 low bits of
      the 128-bit integer h[4..19] are zero, Y is distinguished: leave the
      loop. If the global step counter equals S, stop with failure. If this
      walk has made L steps, abandon it (go to 3). Otherwise
      j=J[h[0..3] mod 11]; Y <- Y^2*phi^j(Y); C <- C*c_j mod r; increment
      both counters.
   c. Look up t(Y) in the table of stored distinguished points. For each
      stored (Y',C',a',b') with the same t, test Y=phi^k(Y') for k=0..11.
      On a match solve, with all arithmetic mod r,

          x = (rho^k*C'*a' - C*a) / (C*b - rho^k*C'*b'),

      failing if the denominator is 0. Verify g^x=z (failure otherwise) and
      go to step 4. If no match, store (Y,C,a,b) under key t(Y).
4. Output P=P0 and Q=[x]Q0 by affine double-and-add in Fp2.

`reference.py` is the executable reference with this exact control flow
(`solve`), generic over the group. There is one bounded trial, no hidden
restart and no amortization. All coins are fresh; nothing is precomputed
for other targets; nonuniform advice is zero.

## 3. Correctness

Label invariant: along a walk Y = g^(C*(a+b*x)), initially by construction
and afterwards by F2. If Y=phi^k(Y') at two distinguished points then by F1
C(a+bx) = rho^k C'(a'+b'x) mod r, which is the displayed formula whenever
its denominator is nonzero; r is prime so the inverse exists. The final
check g^x=z means a wrong scalar is never output. By bilinearity
e(P0,[x]Q0)=e(P0,Q0)^x=g^x=z, and both points lie in the prescribed
subgroups, which is all the track requires. A false match of t between
different classes only costs twelve comparisons (charged) and is then
stored as an ordinary point.

Degeneracy is exactly 1/r, with no heuristic. Fix x and condition on the
salt and on the starting logarithms s_n=a_n+b_n*x of all walks. The whole
group computation (every Y, every digest, every C, which two walks n'<n
meet, and k) is a function of these. Conditionally the b_n remain
independent and uniform (for each b_n exactly one a_n gives s_n). The
denominator is C*b_n - rho^k*C'*b_n' with C nonzero and n != n' (each walk
stores at most one point and stops there), hence uniform on Z/r: it
vanishes with probability exactly 1/r.

## 4. Work bound (hard cap, independent of all heuristics)

`machine.md` derives, for a scalar AArch64 implementation,

| per step | bound |
|---|---:|
| Fp12 multiplication (tower Karatsuba, 54M+280A) | 22,784 |
| Fp12 squaring (complex method, 36M+210A) | 15,936 |
| Frobenius power phi^j (15M+31A) | 4,832 |
| SHA-256, one block, with message assembly | 3,200 |
| C <- C*c_j mod r | 256 |
| selector, distinguished test, counters, copies, calls | 992 |
| **I, instructions per step** | **48,000** |

with M=256 and A=32 instruction bounds for one 4-limb Montgomery
multiplication and one modular addition/subtraction. The loop executes at
most S steps in total, over all walks, whatever happens. Everything else
(input checks, constants, at most 2^20 walk starts with two exponentiations
each, at most 2^20 table insertions/lookups with conjugate tests,
recovery, verification, the G2 scalar multiplication, output, and all
failure paths) is below 2^50 instructions (machine.md) and is charged as
F=2^64. Hence

    R_total <= S*48,000 + 2^64,   W = R_total/20,000,
    log2 W <= 126.853862  (rounded up; recomputed by experiment.py).

One step costs 2.4 reference units; this is an instruction count of actual
BN254 arithmetic, not a substitution of the normalization unit for a BN254
multiplication. Peak storage is at most 2^20 records of 512 bytes plus
under 2^20 bytes of code, constants and scratch: below 2^30 bytes.
Preprocessing is below 2^50 instructions and is included in F.

## 5. Success probability: heuristic, support and limits

Collision model. Let f be the map on the N classes induced by the walk
(the class of Y goes to the class of f_j(Y), j selected by the salted hash
of t). A walk started at a fresh uniform class succeeds in producing a
match as soon as its trail touches the trail of an earlier walk, provided
it then reaches that walk's distinguished point (the trails coincide from
the meeting point on, because f is deterministic on classes). Under the
standard heuristic that the new points behave like independent uniform
classes, the probability that no trail has touched an earlier trail after s
total steps is about exp(-lambda*s^2/(2N)). For an ideal random map
lambda=1. For an m-adding walk the conventional Brent-Pollard/Teske
correction is lambda=1-1/m; here the class walk is exactly an 11-adding
walk in the cyclic group (Z/r)^*/<rho> (written multiplicatively:
log -> c_j*log), so we use lambda=10/11. With S>=1.5*sqrt(N):

    predicted touch probability >= 1-exp(-(10/11)*1.125) = 0.640.

Losses, all explicit:
* A touch in the last L steps may not reach its distinguished point before
  the cap. The touch rate per step is at most S/N, so this loses at most
  L*S/N <= 1.5*2^116/sqrt(N) < 2^-8.
* A walk is abandoned without a distinguished point only if 2^116 steps
  miss a 2^-110 event: probability about exp(-64) per walk under the same
  uniformity heuristic for digest bits, times at most 2^20 walks.
* A walk closes a cycle on itself before distinguishing with probability
  about E[len^2]/(2N) <= 2^221/2^251 per walk; it is then abandoned at L and
  wastes at most 2^116 steps. Over the about S*theta=2^15.6 expected walks
  this is below 2^-14.
* The 2^20 walk cap binds only if 2^20 walks average under 2^95.6 steps
  each against a mean of 2^110: negligible.
* Degeneracy 1/r (exact), scalar rejection-sampling failure below
  2^21*2^-1024, start Y=1 with probability 1/r per walk.
These total below 0.01, leaving a heuristic prediction of about 0.63
against the declared **success >= 0.50** and the protected minimum 0.39.
The identity target succeeds directly.

What is and is not assumed. Not assumed: any property of a fixed seed, an
ideal hash oracle at zero cost, injectivity of t, or that the cap on work
depends on the heuristic. Assumed: that the salted SHA-256 selector applied
to the trace coefficient, together with the eleven multipliers c_j, yields
class-walk collision statistics at production size no worse than the
corrected random-map model by the margin 0.63 -> 0.50. SHA collision
resistance does not prove this. Known structural caveats: only 11
multipliers, one multiplicative relation among them (F4), and the walk
preserves cosets of the subgroup of (Z/r)^*/<rho> generated by the c_j.
The last point does not change pairwise meeting rates: two fresh uniform
classes lie in the same coset with probability 1/d and then meet at rate
d/N, i.e. 1/N overall, so the birthday estimate is unchanged whatever the
index d is.

Support. (a) The construction is the classical use of a Frobenius
automorphism in parallel collision search (Wiener-Zuccherato; Gallant,
Lambert and Vanstone, Math. Comp. 69 (2000), for Koblitz curves) combined
with van Oorschot-Wiener distinguished points (J. Cryptology 12 (1999));
these are cited as context, and the needed reasoning is reproduced above
rather than imported. (b) `experiment.py` runs the complete algorithm
(`solve`, same selector, same hash, same cap formula S, same recovery and
failure rules) on prime-order additive groups Z/q with an automorphism of
exact order 12, where the class walk is the same multiplication by
2+rho^j with the same Z[zeta_12] relations. Development run (seed
20261111, 256/256/64 trials at q=120097, 1200061, 12000097): success
fractions at the cap 0.750, 0.738, 0.766, and uncapped mean steps
0.90, 0.95, 0.99 times sqrt(pi*N/2). The toy runs use theta=2^-2..2^-4, so
their distinguished-point overhead is relatively much larger than at
production size (2^-15 of S). (c) At full BN254 width the experiment checks
F1-F3, the label invariant, class-compatibility of the actual walk with
the actual SHA selector, and the recovery formula on a planted conjugate
pair (the planting uses x in the test harness only).

Limits. Toy groups are 17 to 24 bits; in them the invariant is the orbit
minimum rather than a field trace; the extrapolation to 2^250 classes and
2^125 steps is the disclosed heuristic, not a theorem, and no production
discrete logarithm has been computed. The experiments are evidence only at
their scope. Finite-field NFS-type methods may well be stronger and remain
eligible; nothing here claims this score is tight.

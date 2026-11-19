# Fixed-G1 pairing inversion into G2 by one GT discrete logarithm (BLS12-381)

This is a generic small-memory 32-adding Pollard-rho attack for the protected
BLS12-381 FAPI-1 problem: given the instance's nonidentity P in G1 and z in GT,
return Q in G2 with e(P,Q)=z, where e is the protected pairing
`py_ecc.optimized_bls12_381.pairing(Q,P)` (library takes G2 first). It uses the
admissible route named in rules.md: solve z=e(P,Q0)^x in GT and return [x]Q0.
It is not a best-known finite-field attack, a novelty claim, a practical break
or a security-bit estimate. The work bound is analytical; the success
justification uses the standard, explicitly disclosed Pollard-rho random-walk
heuristic. It seeks the protected exploratory `plausible_not_refuted` standard.
No prior result is assumed accepted. The rho and accounting structure follow
the selected BN254 pairing-inversion rho baseline in this repository, re-derived
here for BLS12-381 (six-limb Fp, the protected w^12=2w^6-2 basis, a 255-bit r,
576-byte GT encoding, and the explicit cost of the one pairing e(P,Q0)).

## Complete algorithm

1. Read the public profile and the instance (P, z). Decode P and the 12
   canonical coefficients of z; check P on the G1 curve, P != infinity,
   z != 0 and z^r = 1. Convert to six-limb Montgomery form.
2. If z = 1 return Q = infinity (valid exactly for z = 1 by the profile).
3. Compute the base g = e(P,Q0) **once**, by evaluating the protected
   pairing function on the actual input P and the public generator Q0
   (Miller loop with the protected ate loop count, then the protected final
   exponent (p^12-1)/r). Cost is charged explicitly in machine.md.
4. Let N = 2*ceil(sqrt(r)) with the protected r. There is one bounded trial.
   Sample fresh 256-bit salt and 66 independent uniform scalars modulo r by
   bounded 255-bit rejection sampling (at most 1024 draws per scalar).
   For each of 32 pairs (c_j,d_j) build M_j = g^c_j * z^d_j. From the
   remaining pair form Y0 = g^a0 * z^b0. A state is (Y,a,b) with Y=g^a*z^b.
   The deterministic iteration is

       j = SHA256(salt || fixed_encoding(Y))[0] & 31
       (Y,a,b) <- (Y*M_j, a+c_j mod r, b+d_j mod r).

   fixed_encoding is the injective 576-byte internal Montgomery encoding of
   machine.md (12 canonical residues x 48 bytes), not an unpriced oracle.
5. Floyd: tortoise and hare start at the same state. For k=1..N advance the
   tortoise once and the hare twice and compare all 12 GT coefficients. At the
   first equality with labels (a,b) and (A,B): if B-b = 0 mod r report
   failure; else x = (a-A)/(B-b) mod r. If no equality within N, failure.
6. Verify g^x = z; a verification failure is a failed trial, never a wrong
   output. Return Q = [x]Q0 (affine Fp2 double-and-add, machine.md).
   There is no hidden restart, uncharged amplification or amortization.

`rho.py` gives the executable mathematical reference and exact loop order.

## Correctness of the output

g = e(P,Q0) where e is exactly the protected function. Since P is a nonidentity
point of the prime-order-r group G1 and Q0 generates G2, nondegeneracy gives
that g has order exactly r, so g generates GT (order r). Hence every z in GT,
including z=1, equals g^x for a unique x mod r. If g^x=z (checked in step 6)
then bilinearity of the protected pairing in its G2 argument gives

    e(P,[x]Q0) = e(P,Q0)^x = g^x = z.

[x]Q0 is an order-r multiple of Q0, hence on the G2 curve with [r]Q=infinity,
as acceptance requires; x=0 only occurs for z=1 where infinity is returned.
The exact P of the instance is used; no P is substituted and no scalar is
returned instead of Q. Bilinearity of the pinned py_ecc pairing (an optimal
ate pairing with |x| loop count followed by the full final exponentiation) is
the standard property the challenge itself relies on ("solve z=e(P,Q0)^a in
GT and return Q=[a]Q0" is named admissible in rules.md); the challenge is
generated as z=e(P,Q0)^a with the same function, so for success we only need
the value g computed in step 3 to equal that function's value, which holds
by evaluating that same deterministic function.

## Correctness given a nondegenerate collision

The invariant Y=g^a*z^b holds initially and after every step by exponent
addition. If g^a*z^b=g^A*z^B and z=g^x, then a-A=x(B-b) mod r. Since r is
prime, nonzero B-b is invertible and x is the unique logarithm mod r. All
scalar operations use four words reduced mod r. The final check g^x=z
validates every returned scalar.

For a deterministic state map on group elements, write mu for the tail and
lambda for the cycle length. Floyd finds an equality at a positive multiple k
of lambda with k>=mu, and such k<=mu+lambda exists. So a first-repeat length
mu+lambda<=N is sufficient for detection within the charged N iterations.
Equality detection is on group elements only (labels are not hashed).

## Success heuristic, support and limits

For an ideal uniform random map on r points, the no-repeat probability after
t transitions is about exp(-t^2/(2r)) at t=Theta(sqrt(r)). A 32-adding map is
not an independent random map; the conventional r-adding correction gives an
effective rate approximately 1-1/32. At t=2*ceil(sqrt(r)) this predicts
repeat probability

    1-exp(-(31/32)*2) ~= 0.856.

This is a stated heuristic prediction, not an exact tail theorem for this
field and SHA-256. Degeneracy probability is exactly 1/r (below). The
prediction leaves large margin over the declared **success >= 0.50** (and the
protected minimum 0.39). Rejection-sampling failure: one 255-bit draw accepts
with probability r/2^255 > 0.9 > 1/2, so the union bound gives failure
< 66*2^-1024. The identity target succeeds directly. The uniform input P only
enters through g, which is always a generator; the walk analysis is the same
for every generator g, so the heuristic is uniform over P.

The material assumption is that the salted SHA partition and random 32
increments do not cause anomalous collision-time behavior or slow mixing at
production size. SHA collision resistance alone does not prove this. It
concerns walk statistics over fresh increments, salt, starting labels and the
uniform input, not a claim that a fixed seed is representative. The work and
memory caps hold even if the heuristic fails; the success lower bound is
exploratory rather than an unconditional theorem.

Degeneracy needs no random-relation assumption. Fix the true x (unknown to the
algorithm) and condition on the salt, the 32 increment logarithms
s_j=c_j+x*d_j mod r and the starting logarithm s0=a0+x*b0. Conditional on
these, the d_j remain independent uniform scalars (each d_j has exactly one
c_j realizing s_j). The group walk, its partition choices and the first Floyd
equality depend only on s_j, s0 and salt, not on the d_j. At a collision
between times k and 2k, B-b = sum(n_j*d_j) mod r with nonnegative counts
summing to k, 0<k<=N<r, so some n_j is nonzero mod r; this nonzero linear
combination of independent uniform field elements is exactly uniform, so the
degeneracy probability is exactly 1/r. This applies to the actual SHA
partition. Only the collision-within-cap probability is heuristic.

Support: the standard r-adding construction and random-map analysis above, the
sizeable margin, and the supplied experiment. Teske, "On random walks for
Pollard's rho method", Math. Comp. 70(234), 809-825 (2001),
DOI 10.1090/S0025-5718-00-01213-8, studies this class of walks (context only;
the reasoning used is reproduced here).

`experiment.py` tests independently chosen targets, increments and salts on
four prime additive groups using the same 32-way SHA partition, the actual
576-byte encoding length, coefficient update, cap, recovery and failure rules,
reporting all outcomes including failures (development run: 93-96% success
per group, `development-evidence.json`). It checks six-limb (p) and four-limb
(r) Montgomery multiplication on full-width edge/random inputs, builds a
full-size order-r GT element in the protected basis as h^((p^12-1)/r) and
checks the walk invariant at full width. These are correctness and finite
walk-statistics evidence only; extrapolation to 2^127-scale walks is part of
the disclosed heuristic. The experiment does not evaluate the py_ecc pairing.

This is an ordinary generic-attack bound. Finite-field index calculus / (Sex)TNFS
in Fp12 may be stronger and remains eligible; nothing here claims tightness.

## Computation, memory and output

`machine.md` gives a complete bounded scalar implementation: six-limb and
four-limb Montgomery multiplication/REDC, polynomial Fp12 arithmetic in the
protected basis, scalar SHA-256, fixed table access, three-word counters,
bounded randomness, the one pairing evaluation, setup, failure handling,
recovery and affine G2 output. Peak storage <= 1 MiB.

Every Floyd iteration costs <= 6,774,784 scalar instruction equivalents. An
additional 2^40 instructions covers all setup and final work, including the
pairing e(P,Q0), subgroup checks, final verification and [x]Q0. Thus

    R_total = 2*ceil(sqrt(r))*6,774,784 + 2^40,
    W = R_total/20,000,
    log2(W) <= 136.832576   (rounded up at the sixth decimal).

`experiment.py` recomputes this from the protected r. Memory log2 bytes <= 20;
preprocessing (all one-target, charged in full inside the 2^40) is < 2^40
instructions; nonuniform advice is zero. Only public profile/SHA constants and
algorithm design are excluded. No Python timing prices the arithmetic.

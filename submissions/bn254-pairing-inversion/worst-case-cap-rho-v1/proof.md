# Small-memory BN254 reduced-cap Pollard-rho variant

This is a generic 32-adding Pollard-rho baseline for the actual protected
254-bit prime-order GT. It is not a best-known finite-field attack, a novelty
claim, a practical break, or a security-bit estimate. Its work bound is
analytical and its success justification uses the standard, explicitly
disclosed Pollard-rho random-walk heuristic. This submission seeks the
protected exploratory `plausible_not_refuted` standard, not formal proof
certification. No prior result is assumed accepted.

## Complete algorithm

Read only the public profile and z. Decode its 12 canonical coefficients,
check z!=0 and z^r=1, and convert to four-limb Montgomery coefficients.
If z=1 return scalar zero (or P0 and infinity in the pairing track).
Let N=ceil(1.1*sqrt(r)), the smallest integer with 100*N^2 >= 121*r
(N=162741432569968250906303761388095885571 for the protected r). There is one bounded trial.

Sample fresh 256-bit salt and 66 independent uniform scalars modulo r,
using bounded 254-bit rejection sampling (at most 1024 draws per scalar).
For each of 32 pairs (c_j,d_j), construct M_j=g^c_j*z^d_j. From the remaining
pair form Y0=g^a0*z^b0. A state is (Y,a,b), maintaining Y=g^a*z^b.
The deterministic iteration is

    j = SHA256(salt || fixed_encoding(Y))[0] & 31
    (Y,a,b) <- (Y*M_j, a+c_j mod r, b+d_j mod r).

fixed_encoding is the injective 384-byte internal Montgomery encoding
specified in machine.md, not an unpriced conversion or hash oracle.

Initialize tortoise and hare to that same initial state. For k=1,...,N,
advance tortoise once and hare twice and compare all GT coefficients.
At the first equality, let (a,b) and (A,B) be their scalar labels.
If B-b=0 mod r, report failure. Otherwise compute

    x = (a-A)/(B-b) mod r.

Verify g^x=z before output; verification failure is a failed trial, never an
accepted wrong answer. If no equality within N iterations, report failure.
There is no hidden restart, uncharged amplification or amortization.
`rho.py` gives the executable mathematical reference and precise loop order.

The GT track returns x. The pairing track returns P=P0 and Q=[x]Q0 using
the explicitly charged affine procedure in machine.md. This is allowed
by the arbitrary-preimage-pair rules. Bilinearity gives
e(P0,[x]Q0)=g^x=z; both points are in the prescribed subgroups.
For the identity, the permitted infinity output is used.

## Correctness given a nondegenerate collision

The invariant follows initially and after every step by exponent addition.
If g^a*z^b=g^A*z^B and z=g^x, then a-A=x(B-b) mod r. Since r is prime,
nonzero B-b has a unique inverse. The resulting x is the unique discrete
logarithm modulo r. All scalar operations use four words and are reduced
modulo r. Neither a hidden sampled scalar nor any source-group preimage
is an input. The final exponentiation validates every returned scalar.

For a deterministic state map on group elements, write mu for the tail
length and lambda for the cycle length. Floyd has an equality at an
iteration k that is a positive multiple of lambda with k>=mu; there is
such a k<=mu+lambda. Consequently a first-repeat length mu+lambda<=N
is sufficient for detection within the charged N iterations. Scalar labels
are not part of the hash key; equality detection is on group elements only.

## Success heuristic, support and limits

For an ideal uniform random map on r points with uniform starting point,
the no-repeat probability after t distinct transitions is the product
of (1-j/r), and is approximately exp(-t^2/(2r)) at t=Theta(sqrt(r)).
A 32-adding map is not an independent random map. The conventional
r-adding correction accounts for the probability of selecting the same
increment; its effective rate is approximately 1-1/32. Applying it at
t=1.1*sqrt(r) predicts repeat probability

    1-exp(-(31/32)*(1.1)^2/2) ≈ 0.4435

(and ≈0.4329, ≈0.4221 at the neighboring rates 30/32 and 29/32).

This is a stated heuristic prediction, not an exact tail theorem for this
fixed field and SHA-256. The coefficient degeneracy probability is exactly
1/r by the argument below. Thus the heuristic
prediction leaves margin over the declared **success >=0.40**
(and the protected minimum 0.39). Rejection-sampling failure is negligible:
one 254-bit draw accepts with probability r/2^254>1/2, so union bound gives
failure <66*2^-1024 for the scalar generation. Identity succeeds directly.

The material assumption is that the salted SHA partition and randomly
generated 32 increments do not cause anomalous collision-time behavior or
slow mixing at production size. SHA collision
resistance alone does not prove them. They concern walk statistics across
fresh increments, salt, starting labels and uniform full-production input,
not a claim that a fixed seed is representative or that the whole walk is IID.
The score's work and memory caps hold even if these heuristics fail; the
success lower bound is exploratory rather than a proven unconditional theorem.

Degeneracy needs no random-relation assumption. Fix the input's true x
(without giving it to the algorithm) and condition on the salt, the 32
increment logarithms s_j=c_j+x*d_j mod r, and the starting logarithm
s0=a0+x*b0 mod r. Conditional on these, all d_j remain independent uniform
scalars: for each d_j there is exactly one c_j realizing the conditioned s_j.
The group walk, its partition choices and the selected first Floyd equality
depend only on s_j,s0,salt, not on the remaining random d_j.
At a collision between times k and 2k, the coefficient difference B-b is
the sum of the d_j over the k intervening steps, namely sum(n_j*d_j) mod r.
The nonnegative counts satisfy sum(n_j)=k with 0<k<=N<r. At least one n_j
is nonzero modulo r. This nonzero linear combination of independent uniform
field elements is exactly uniform, so degeneracy probability conditional
on any such group walk is exactly 1/r. Bounded rejection sampling is uniform
conditional on successful sampling. This argument applies to the actual
deterministic SHA partition as well as to an idealized partition.
Thus only the probability of a group collision within the cap is heuristic.

Support is the standard r-adding construction, its random-map analysis above,
the sizeable probability margin, and directly supplied experiments. Teske,
“On random walks for Pollard's rho method”, Mathematics of Computation
70(234), 809–825 (2001), DOI 10.1090/S0025-5718-00-01213-8, studies this
class of walks; we use 32 additions rather than Pollard's less-random original
three-part walk. This reference is context, not an external unproved theorem
imported to supply our precise probability. The relevant algorithm and
heuristic reasoning are reproduced here; no web retrieval is required.

`experiment.py` tests independently chosen targets, increments and salts on
four prime additive groups, using the same 32-way hash, 384-byte encoding
length, coefficient update, cap, recovery and failure rules. It reports all
outcomes, including failures. It also checks Montgomery multiplication on
full-width edge/random inputs for both p and r, the protected BN254
generator order, and full-width BN254 exponent invariants. These are
correctness and finite walk-statistics evidence only. They neither solve
production DLP nor prove a probability bound at r. Additive-group
encodings differ from the field representation; extrapolation to GT and
2^127-scale walks is part of the explicitly disclosed heuristic.

This is an ordinary generic-attack baseline. Finite-field index calculus
and NFS can be stronger and remain eligible; nothing here proves this score
is tight or is the appropriate reward threshold for such algorithms.

## Computation, memory and output

`machine.md` gives a complete bounded scalar implementation: explicit
four-word Montgomery multiplication/REDC, polynomial Fp12 arithmetic,
scalar SHA-256, fixed table access, three-word iteration counters,
bounded randomness, setup and failure handling, scalar recovery, and
affine Fp2 output conversion with exceptional cases. Peak storage is <=1 MiB.
No unbounded table, wide memory address or network transfer is assumed.

Every iteration costs <=3,420,160 scalar instruction equivalents. An
additional 2^40 instructions covers all setup and final work, including the
pairing track's conversion and optional on-curve/subgroup checks. Thus for
both tracks

    R_total = 162741432569968250906303761388095885571*3,420,160 + 2^40,
    W = R_total/20,000.

Submit the upward-rounded six-decimal log2(W) as an exact rational in the
manifest. `experiment.py` recomputes it from the protected order; a shared
bound for both tracks is intentional, not a claim of equal actual work.
Memory log2 bytes <=20; preprocessing <=2^36 scalar instructions is charged
in full; nonuniform advice is zero. Public profile/SHA constants and algorithm
design are the only excluded public data. No Python timings, historical
measurements or normalization-unit substitution price the arithmetic.

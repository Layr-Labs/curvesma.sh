# Explicit cold field-basis conversion and its budget

This supplies the mapping step for the new wrappers; no precomputed isomorphism
is free. It is conditional on the inherited, explicitly checked irreducibility
claims for h and g. Those claims must be reproduced before publication.

Let F be the protected polynomial representation of F_(p^12). In F[Y], factor
h(Y)=Y^6+Y^3+Y^2-Y-1 into linear factors and choose a root y. Because h is
irreducible of degree six over Fp and 6 divides 12, it splits over F and every
root lies in its degree-six subfield. Then find a root x in F of
X^2+|u|y. The inherited quadratic irreducibility over F_(p^6) implies this
polynomial splits in F and x is outside that subfield.

The map Y -> y, X -> x from the source tower to F is a field isomorphism:
relations hold, the source quotient is a field of dimension 12, and a unital
map of that field into F is injective and therefore bijective. The 12 vectors
y^i x^j, 0<=i<6, 0<=j<2, form an Fp basis. Form their coordinates in the
protected basis and invert the resulting 12-by-12 matrix using field Gaussian
elimination. It maps arbitrary input g,z into the STNFS tower. Use the same
map for both; any field isomorphism preserves z=g^a and hence the scalar a.
No special choice of a root, hidden challenge exponent or free advice is needed.
All basis construction and input conversion are charged once per attack.

## Bounded randomized factorization

A completely splitting squarefree polynomial of degree d<=6 over F can be
split by choosing a polynomial R of degree<d with independent uniform F
coefficients, taking gcd(R,f), and, if needed,
gcd(R^((p^12-1)/2)-1,f). Repeat only for unresolved factors. Evaluations at
the d distinct roots are jointly uniform by the invertible Vandermonde map.
For any two roots the probability that their nonzero values have opposite
quadratic characters is ((|F|-1)/|F|)^2/2 >1/4. This event gives a nontrivial
split. Thus each trial succeeds with probability at least 1/4, independently
of earlier trials conditional on the current factor. Five splits suffice for
h and one for the quadratic. Cap each split at 128 trials. By the union bound,
failure probability is at most 6*(3/4)^128 <2^-50. All aborts count as failure.

Uniform F elements are built from twelve independent Fp coordinates.
Each coordinate uses 381-bit rejection sampling with at most 128 attempts;
since p>2^380, failure per coordinate is below 2^-128. Fewer than 2^16
coordinates are required for the six times 128 trials and degree<=6, so the
additional failure is below 2^-112. Exact verification of roots and the inverse
matrix adds no probabilistic assumption. A loose total extra failure allowance
of 2^-48 covers both terms. Random bits, rejection and failed trials are charged.

## Explicit deliberately loose cost ceiling

Use the elementary scalar integer Fp and Fp12 bounds in proof.md; no scalar
floating-point, SIMD or arbitrary-width unit-cost arithmetic is required.
For a polynomial modulus of degree<=6, a multiplication and reduction use
fewer than 2^10 field operations. Binary powering with fewer than 4572 bits
therefore takes fewer than 2^24 field operations. Euclidean gcd and associated
coefficient inversions can also be bounded by 2^24 field operations here,
using binary-powering field inverses. Charge 2^26 Fp12 operations per trial.
At most 6*128<2^10 trials gives fewer than 2^36 Fp12 operations. At fewer
than 2^31 scalar instructions per Fp12 multiplication/addition-equivalent,
this is below 2^67 instructions. Matrix inversion, validation, pairing,
point/scalar output arithmetic, coordinate generation and control fit with
large slack in a total wrapper cap of 2^80 scalar instructions. This includes
the factorization, not only the later linear coordinate map.

The maps and polynomials occupy fewer than 2^30 bytes. This is an analytic
upper-bound implementation argument, not a claim that a fast implementation
has been benchmarked. The inherited STNFS kernels require a separate audit;
this wrapper cannot make an unsupported core price valid.

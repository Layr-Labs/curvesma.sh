# Explicit cold field-basis conversion and its budget

The wrapper mapping step; no precomputed isomorphism is free. Needs checked
irreducibility of h and g (evidence16.py [0]/[0b]).

Let F be the protected polynomial representation of F_(p^12). In F[Y], factor
h(Y)=Y^6+Y^3+Y^2-Y-1 and choose a root y (h splits over F; roots in the
degree-six subfield). Then find a root x in F of X^2+|u|y (splits in F; x
outside the subfield).

Y -> y, X -> x is a field isomorphism; the 12 vectors y^i x^j form an Fp basis.
Coordinatize and invert the 12-by-12 matrix (Gaussian elimination). One map
serves both inputs; no root choice, hidden exponent or advice. All construction
and conversion charged once per attack.

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

Uniform F elements use twelve independent Fp coordinates with 381-bit
rejection sampling (<= 128 attempts; failure below 2^-128 each since p>2^380).
Fewer than 2^16 coordinates cover six times 128 trials at degree<=6 (failure
below 2^-112). Root/inverse verification adds no probabilistic assumption. A
loose 2^-48 allowance covers both; random bits, rejection and failed trials
are charged.

## Loose cost ceiling

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
large slack in a total wrapper cap of 2^80 scalar instructions, factorization
included.

Maps and polynomials occupy fewer than 2^30 bytes (analytic upper bound).

# Explicit cold field-basis conversion and its budget

The wrapper mapping step; no precomputed isomorphism is free; needs h, g irreducibility ([0]/[0b]).

Let F be the protected representation of F_(p^12). In F[Y], factor h(Y)=Y^6+Y^3+Y^2-Y-1 and choose a root y (roots in the degree-six subfield). Then find a root x in F of X^2+|u|y (x outside the subfield).

Y -> y, X -> x is a field isomorphism; the 12 vectors y^i x^j form an Fp basis. Coordinatize and invert the 12-by-12 matrix by Gaussian elimination. One map serves both inputs; no root choice, hidden exponent or advice; all construction and conversion charged once per attack.

## Bounded randomized factorization

A completely splitting squarefree polynomial of degree d<=6 over F splits by choosing R of degree<d with independent uniform F coefficients, taking gcd(R,f) and, if needed, gcd(R^((p^12-1)/2)-1,f). Repeat only for unresolved factors. Evaluations at the d distinct roots are jointly uniform by the invertible Vandermonde map.
For any two roots, opposite quadratic characters of their nonzero values has probability ((|F|-1)/|F|)^2/2 >1/4, giving a nontrivial split. Each trial succeeds with probability at least 1/4, independently of earlier trials given the current factor. Five splits suffice for h, one for the quadratic; cap each at 128 trials. Union bound: failure at most 6*(3/4)^128 <2^-50. All aborts count as failure.

Uniform F elements use twelve independent Fp coordinates with 381-bit rejection sampling (<= 128 attempts; failure below 2^-128 each since p>2^380); fewer than 2^16 coordinates cover 6x128 trials at degree<=6 (failure below 2^-112). Root/inverse verification adds no probabilistic assumption. A loose 2^-48 allowance covers both; random bits, rejection and failed trials are charged.

## Explicit deliberately loose cost ceiling

With the elementary scalar Fp and Fp12 bounds of proof.md: multiplication and reduction modulo a degree <= 6 polynomial cost < 2^10 field operations, so binary powering with < 4572-bit exponents costs < 2^24, and Euclidean gcds with coefficient inversions are also bounded by 2^24. Charge 2^26 Fp12 operations per trial; at most 6x128 < 2^10 trials give < 2^36 Fp12 operations at < 2^31 scalar instructions each: below 2^67. Matrix inversion, validation, pairings, output, coordinate generation and control fit with large slack in the 2^80 wrapper cap, which includes the factorization. Maps and polynomials occupy < 2^30 bytes. Analytic upper bound, not a benchmarked implementation.

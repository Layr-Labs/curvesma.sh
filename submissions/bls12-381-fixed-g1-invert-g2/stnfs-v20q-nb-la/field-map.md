# Explicit cold field-basis conversion and its budget

The wrapper's mapping step; no precomputed isomorphism is free. It needs the checked irreducibility of h and g
(evidence16.py [0]/[0b]).

Let F be the protected polynomial representation of F_(p^12). In F[Y], factor h(Y) = Y^6+Y^3+Y^2-Y-1 and choose a root
y (h splits over F; roots in the degree-six subfield); then find a root x in F of X^2+|u|y (it splits in F; x is outside
the subfield). Y -> y, X -> x is a field isomorphism and the 12 vectors y^i x^j form an Fp basis: coordinatize and
invert the 12-by-12 matrix by Gaussian elimination. One map serves both inputs; no root choice, hidden exponent or
advice; construction and conversion are charged once per attack.

## Bounded randomized factorization

A completely splitting squarefree polynomial f of degree d <= 6 over F is split by choosing R of degree < d with
independent uniform F coefficients, taking gcd(R, f) and, if needed, gcd(R^((p^12-1)/2) - 1, f), repeating only for
unresolved factors. Evaluations at the d distinct roots are jointly uniform (invertible Vandermonde map). For any two
roots, their nonzero values have opposite quadratic characters with probability ((|F|-1)/|F|)^2/2 > 1/4, which splits
them: each trial succeeds with probability >= 1/4, independently of earlier trials given the current factor. Five
splits suffice for h and one for the quadratic; each split is capped at 128 trials, so by the union bound the failure
probability is <= 6 (3/4)^128 < 2^-50. All aborts count as failure. Uniform F elements use twelve independent Fp
coordinates by 381-bit rejection sampling (<= 128 attempts, failure < 2^-128 each since p > 2^380); fewer than 2^16
coordinates cover six times 128 trials at degree <= 6 (failure < 2^-112). Root and inverse verification add no
probabilistic assumption. A 2^-48 allowance covers both; random bits, rejections and failed trials are charged.

## Explicit loose cost ceiling

Scalar integer Fp and Fp12 bounds of proof.md; no floating point, SIMD or unit-cost wide arithmetic. Modulo a
polynomial of degree <= 6, a multiplication and reduction take < 2^10 field operations, so binary powering with fewer
than 4572 bits takes < 2^24; Euclidean gcd and coefficient inversions (binary-powering inverses) also take < 2^24.
Charge 2^26 Fp12 operations per trial: at most 6 x 128 < 2^10 trials give < 2^36 Fp12 operations, at < 2^31 scalar
instructions each below 2^67 instructions. Matrix inversion, validation, pairing, point and scalar output arithmetic,
coordinate generation and control fit with large slack in a wrapper cap of 2^80 scalar instructions, factorization
included. Maps and polynomials occupy fewer than 2^30 bytes (analytic upper bound).

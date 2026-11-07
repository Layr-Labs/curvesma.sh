# Explicit cold field-basis conversion and its budget

The wrapper's mapping step (no precomputed isomorphism is free); it needs the checked irreducibility of h and g
(evidence16.py [0]/[0b]). In F[Y], F the protected representation of F_(p^12), factor h(Y) = Y^6+Y^3+Y^2-Y-1 and take a
root y (h splits over F, roots in the degree-six subfield), then a root x of X^2+|u|y (it splits in F, x outside the
subfield). Y -> y, X -> x is a field isomorphism and the 12 vectors y^i x^j form an Fp basis; coordinatize and invert
the 12 x 12 matrix by Gaussian elimination. One map serves both inputs (no root choice, hidden exponent or advice);
everything is charged once per attack.

**Bounded randomized factorization.** A completely splitting squarefree f of degree d <= 6 over F splits by taking R of
degree < d with uniform coefficients and gcd(R, f) or gcd(R^((p^12-1)/2) - 1, f), repeating on unresolved factors. The
values at the d roots are jointly uniform (invertible Vandermonde map); two roots get opposite quadratic characters
with probability ((|F|-1)/|F|)^2/2 > 1/4, which splits them, independently of earlier trials. Five splits suffice for h
and one for the quadratic; with each split capped at 128 trials, failure is <= 6 (3/4)^128 < 2^-50 (counted). Uniform
F elements use twelve Fp coordinates by 381-bit rejection sampling (<= 128 attempts, failure < 2^-128 each); fewer than
2^16 coordinates cover all trials (failure < 2^-112). Root and inverse checks add no assumption; a 2^-48 allowance
covers both, and random bits, rejections and failed trials are charged.

**Cost ceiling.** Scalar integer Fp and Fp12 bounds of proof.md. Modulo a degree-<= 6 polynomial a product takes
< 2^10 field operations, powering with < 4572 bits < 2^24, gcds and inversions < 2^24: 2^26 Fp12 operations per trial,
< 2^36 over at most 6 x 128 < 2^10 trials, at < 2^31 instructions each below 2^67 instructions. Matrix inversion,
validation, pairing, scalar output, coordinate generation and control fit with large slack in the 2^80-instruction
wrapper cap, factorization included. Maps and polynomials take fewer than 2^30 bytes (analytic upper bound).

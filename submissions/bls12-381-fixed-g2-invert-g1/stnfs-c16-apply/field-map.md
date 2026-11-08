# Field map
The wrapper's map (no precomputed isomorphism is free; needs h and g irreducible, [0]/[0b]). In F[Y], F the protected
representation of F_(p^12), factor h(Y) = Y^6+2Y^4-2Y^3-2Y^2-1 and take a root y (h splits over F), then a root x of
X^2+|u|(y+1) (it splits in F): Y -> y, X -> x is a field isomorphism, the 12 vectors y^i x^j an Fp basis; invert the 12 x 12
coordinate matrix by Gaussian elimination. One map serves both inputs; all is charged once per attack.

Bounded factorization. A completely splitting squarefree f of degree d <= 6 over F splits by gcd(R, f) or
gcd(R^((p^12-1)/2) - 1, f), R of degree < d uniform, repeated on unresolved factors: the values at the roots are jointly
uniform (Vandermonde), and two roots get opposite quadratic characters with probability ((|F|-1)/|F|)^2/2 > 1/4,
independently of earlier trials. Five splits serve h, one the quadratic; capped at 128 trials each, failure <= 6
(3/4)^128 < 2^-50 (counted). Uniform F elements use twelve Fp coordinates by 381-bit rejection sampling (<= 128
attempts, failure < 2^-128 each; < 2^16 coordinates in all, failure < 2^-112); a 2^-48 allowance covers both, and
random bits, rejections and failed trials are charged.

Cost. Per trial (mod a degree-<= 6 polynomial: powering with < 4572 bits, gcds, inversions) < 2^26 Fp12 operations,
< 2^36 over < 2^10 trials, at < 2^31 instructions each: below 2^67 instructions. Inversion, validation, pairing, scalar output, coordinate generation and
control fit with large slack in the 2^80-instruction wrapper cap. Maps and polynomials take < 2^30 bytes.

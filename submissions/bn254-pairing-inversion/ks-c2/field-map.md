# Field map

The wrapper's map (no precomputed isomorphism is free; needs h and g irreducible, [0a]). In F[Y], F the protected representation
Fp[w]/(w^12 - 18 w^6 + 82) of F_(p^12), factor h(Y) = Y^6 - Y^5 - 2Y^4 + 5Y^3 - 2Y^2 - Y + 1 and take a root y (h splits over F), then a root
x of X^2 - u w(y) (it splits in F): Y -> y, X -> x is a field isomorphism, the 12 vectors y^i x^j an F_p basis; invert the 12 x 12 coordinate
matrix by Gaussian elimination. One map serves both inputs; all is charged once per attack.

Bounded factorization. A completely splitting squarefree polynomial of degree d <= 6 over F splits by gcd(R, f) or gcd(R^((p^12-1)/2) - 1, f),
R of degree < d uniform, repeated on unresolved factors: two roots get opposite quadratic characters with probability > 1/4, independently
of earlier trials. Five splits serve h, one the quadratic; capped at 128 trials each, failure <= 6 (3/4)^128 < 2^-50. Uniform F elements use
twelve F_p coordinates by 254-bit rejection sampling (<= 128 attempts, failure < 2^-128 each); a 2^-48 allowance covers both.

Cost. Per trial (powering with < 3,044 bits mod a degree-<= 6 polynomial, gcds, inversions) < 2^26 F_p12 operations, < 2^36 over < 2^10
trials, at < 2^31 instructions each: below 2^67 instructions. Validation, the GT exponentiations, the pairing check, coordinate generation
and control fit in the 2^68-instruction wrapper row. Maps and polynomials take < 2^30 bytes.

# Field map

The wrapper constructs an explicit field isomorphism from the protected representation `F = F_p[w]/(w^12 - 2 w^6 + 2)` of `F_(p^12)` into the core's tower `F_p[y, X]/(h_pi(y), X^3 - omega^2 u)` once per opened run (no precomputed isomorphism is free).

## Construction

1. **Base quartic root `y`:** Factor the monic degree-4 polynomial `h_pi(Y) = gcd(Phi_15(Y), Y^5 - w_pi) mod p` in `F[Y]`, where `w_pi = (u^3 - u^2 - u)(u + 1)^(-1) mod p` satisfies `w_pi^2 + w_pi + 1 = 0 mod p` (`audit_cm.py`: `algebra_bls12`). Because `F = F_(p^12)` contains `F_(p^4)`, `h_pi` splits into 4 linear factors over `F`. Take one root `y in F` (`w_y = y^5 = w_pi`).
2. **Cubic root `X`:** Factor `g(T) = T^3 - w_y^2 u` in `F[T]` (irreducible over `F_(p^4)` since `(u w_y^2)^((p^4-1)/3) != 1`, verified by `audit_cm.py`). Take one root `X in F`.
3. **Basis and change-of-basis matrix:** Express `{y^i X^j : 0 <= i < 4, 0 <= j < 3}` in `{1, w, ..., w^11}` and invert the 12 x 12 matrix over `F_p` by Gaussian elimination. One map serves both `h` and `T`.

## Bounded Cantor-Zassenhaus failure analysis & cost

- **Quadratic-character splitting:** Over `F` (odd `p`), a squarefree completely splitting polynomial `q` of degree `d in {2, 3, 4}` splits by `gcd(R, q)` or `gcd(R^((|F|-1)/2) - 1, q)` with probability `> 1/4` per uniform `R` of degree `< d`. At most 5 splits resolve degrees 4 and 3; capping each split at 128 trials gives failure `<= 5 (3/4)^128 < 2^(-50)`.
- **Uniform sampling & combined failure:** 381-bit rejection sampling (acceptance `> 0.5`, capped at 128 draws/coord) fails with `< 2^(-112)` across `< 2^16` coords. Total wrapper failure `< 2^(-50) + 2^(-112) < 2^(-48)`, charged in Section 7 of `proof.md`.
- **Cost:** `<= 640` polynomial exponentiations (`< 2^35` `F_(p^12)` ops), matrix inversion, pairings, scalar mul, random coins, and check cost `< 2^52` instructions, well inside the charged `2^80`-instruction (`2^66` core-unit) wrapper row; memory `< 2^30` bytes.

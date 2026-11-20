# Field map

The wrapper constructs an explicit field isomorphism from the protected representation `F = F_p[w]/(w^12 - 2 w^6 + 2)` of `F_(p^12)` into the core's tower `F_p[y, X]/(h_pi(y), X^3 - omega^2 u)` once per opened run (no precomputed isomorphism is free).

## Construction

1. **Base quartic root `y`:** In `F[Y]`, factor the monic degree-4 polynomial `h_pi(Y) = gcd(Phi_15(Y), Y^5 - w_pi) mod p`, where `w_pi = (u^3 - u^2 - u)(u + 1)^(-1) mod p` satisfies `w_pi^2 + w_pi + 1 = 0 mod p` and `h_pi` is irreducible over `F_p` (`audit_cm.py`: `algebra_bls12`). Because `F = F_(p^12)` contains `F_(p^4)`, `h_pi` splits completely into 4 distinct linear factors over `F`. Take one root `y in F`, and set `w_y = y^5 in F` (note `w_y = w_pi in F_p`).
2. **Cubic root `X`:** In `F[T]`, factor `g(T) = T^3 - w_y^2 u` (irreducible over `F_(p^4)` since `(u w_y^2)^((p^4-1)/3) != 1`, verified by `audit_cm.py`). Since `F = F_(p^12)`, `g` splits completely into 3 distinct linear factors over `F`. Take one root `X in F`.
3. **Basis and change-of-basis matrix:** The 12 elements `{y^i X^j : 0 <= i < 4, 0 <= j < 3}` form an `F_p`-basis of `F` because `[F_(p^4):F_p] = 4` and `[F_(p^12):F_(p^4)] = 3`. Express each `y^i X^j` in the protected basis `{1, w, ..., w^11}` to form a 12 x 12 matrix over `F_p` and invert it by exact Gaussian elimination over `F_p`. One map serves both `g` and `z`.

## Bounded Cantor-Zassenhaus failure analysis

- **Quadratic-character splitting of `h_pi` (degree 4) and `g` (degree 3):** Over `F` (odd `p`), a squarefree completely splitting polynomial `q` of degree `d in {2, 3, 4}` splits by `gcd(R, q)` or `gcd(R^((|F|-1)/2) - 1, q)` for uniform `R in F[T]` of degree `< d`. Values at the roots are jointly uniform (Vandermonde), so two distinct roots receive opposite quadratic characters with probability `((|F|-1)/|F|)^2 / 2 > 1/4` per independent trial. At most 3 splits resolve degree 4 (`h_pi`) and at most 2 splits resolve degree 3 (`g`), for at most 5 splits total; capping each split at 128 independent trials gives failure probability `<= 5 (3/4)^128 < 2^(-50)`.
- **Uniform `F` sampling:** Each random `F` element draws 12 `F_p` coordinates by 381-bit rejection sampling (acceptance `p / 2^381 > 0.5`; capped at 128 draws/coordinate, failure `< 2^(-128)` per coordinate, `< 2^(-112)` across `< 2^16` coordinates).
- **Combined wrapper failure:** `5 (3/4)^128 + 2^(-112) < 2^(-50) < 2^(-48)`, charged in full as `2^(-48)` in the strict additive union bound of `proof.md` section 6.

## Cost

Each splitting trial performs polynomial exponentiation modulo a degree `<= 4` polynomial over `F` (`< 4572` squarings/multiplies, `< 2^25` `F_(p^12)` ops). Across `<= 5 x 128 = 640` trials (`< 2^35` `F_(p^12)` ops at `< 2^16` scalar AArch64 instructions each), factorization costs `< 2^52` instructions (`< 2^67` even at an ultra-loose `2^31` instr/`F_(p^12)` op). Matrix inversion (`12^3` `F_p` ops), input validation, pairings, scalar multiplication, random-coin generation (`e_1, e_2`), and final verification add `< 2^50` instructions. Total wrapper work fits with `> 13` bits of slack inside the charged `2^80`-instruction (`2^66` core-unit) wrapper row, and memory is `< 2^30` bytes.

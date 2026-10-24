# Commuting-square Frobenius-class rho

Parent tree `6423c2728dc03d5f2d647c7f22adf1dac0643a66`. One cold inversion of a single pair `(P, z)`. Nothing is amortized over a second target.

## Reduction

The profile fixes a nonidentity `P ∈ G1` and a uniform `z ∈ GT`, and asks for `Q ∈ G2` with `e(P, Q) = z` in the pinned `py_ecc` argument order `pairing(Q, P)`. Let `Q0` be the profile generator of `G2` and set

```text
γ = e(P, Q0).
```

`P = [b]P0` for some unknown `b ∈ {1, …, r−1}`, so `γ = e(P0, Q0)^b`. The pinned pairing is nondegenerate and `r` is prime, so `b` is invertible mod `r` and `γ` generates `GT`. Define `a ∈ {0, …, r−1}` by `z = γ^a`. Then

```text
e(P, [a]Q0) = e(P, Q0)^a = γ^a = z.
```

That is the pairing reduction: one discrete logarithm in `GT` to base `γ`, one scalar multiplication in `G2`, one verification pairing. The identity `z = 1` forces `a = 0` and `Q = ∞`, which the output encoding allows. It is tested by comparing the twelve coefficients of `z` with `[1, 0, …, 0]` before the walk, and it is not sent through the collision equation.

`Q0` and the curve constants are profile data. Their construction is not charged. The pairing that produces `γ` from the instance `(P, Q0)` is charged below.

## Frobenius orbits

Let `p` be the BLS12-381 base-field prime and let `φ(X) = X^p` on `Fp12`. On `GT` this is exponentiation by `p mod r`. The embedding degree is 12, so the multiplicative order of `p` modulo `r` is exactly 12. In particular `12` divides `r − 1`, the non-identity orbits of `⟨φ⟩` have size 12, and

```text
N = (r − 1) / 12
```

is an integer, the number of non-identity orbits. Also `p^6 ≡ −1 (mod r)`: `p^6` is the unique element of order 2 in `(Z/rZ)*`. Hence on `GT`, `φ^6` is inversion, which in the cyclotomic subgroup is the `Fp6` conjugation of the tower (a negation, not a multiplication).

The tower used by the pinned blst sources is

```text
Fp2 = Fp[u]/(u^2 + 1),
Fp6 = Fp2[v]/(v^3 − (u + 1)),
Fp12 = Fp6[w]/(w^2 − v).
```

Write `X = (A, B)` with `A = (a0, a1, a2)` and each `ai ∈ Fp2`. `blst_fp12_frobenius_map` for exponent `n ∈ {1, 2, 3}` copies the real `Fp` limb of `a0` and never multiplies that limb. The same holds for compositions of those maps, and conjugation `φ^6` negates `B` only. Therefore the real `Fp` limb `ρ(X)` of slot `a0` is constant on each Frobenius orbit. Reading it is a load and a mask, not a multiplication.

The walk’s branch bits are the low two bits of `ρ(X)`. They are constant on the orbit, so a square-versus-multiply decision taken from them is a function of the orbit.

Slot `a1` is not constant on the orbit. For the `n = 1` map, blst multiplies the conjugated `a1` by the public constant `(u + 1)^{(p − 1)/3}`, whose real limb is zero in the pinned source (pure imaginary `γ_c · u`, the letter `γ_c` here is a base-field constant, not the pairing target). One application of that induced map sends `x + y u` to `(γ_c y) + (γ_c x) u`. Iterating it, the real limbs of the twelve images are

```text
k = 0:            x
k = 2m > 0 even:  γ_c^{k} · x
k = 2m+1 odd:     γ_c^{k} · y
```

with the powers of `γ_c` public profile constants, computed once. The current representative therefore yields the other eleven real limbs with eleven multiplications in `Fp`. The index `k*` of the unique minimum limb, in the total order of Montgomery residues, is the twist that carries the current representative to the canonical one. If two indices share the minimum limb, the walk stops and the instance fails. A uniform model puts that tie at `O(1/p)` per multiply step. Over the budget, `N_mul / p < 2^{−250}`. This uniqueness statement is a heuristic of the same family as equidistribution of one `Fp` limb of a `GT` element; it is not a Weil-bound theorem.

## Walk

A state is an element `X ∈ ⟨γ⟩` together with integers `(A, B)` such that

```text
X = γ^A · z^B = γ^{A + a B}.
```

The starting state is `X = z`, `(A, B) = (0, 1)`, then replaced by its canonical Frobenius image if the initial real-limb minimum is not already at index 0. That replacement multiplies `(A, B)` by `p^{k} mod r`. It is charged as one distinguished-point materialisation inside preprocessing.

Each later step reads `ρ(X)`:

- If the low bit of `ρ(X)` is 0 and the next bit is 0, multiply by the precomputed element `γ^{p^{−k*}}` and add the precomputed exponent offset `(p^{−k*}, 0)`.
- If the low bit of `ρ(X)` is 0 and the next bit is 1, multiply by `z^{p^{−k*}}` and add `(0, p^{−k*})`.
- If the low bit of `ρ(X)` is 1, replace `X` by its Granger–Scott cyclotomic square and replace `(A, B)` by `(2A, 2B)`.

`k*` is computed from the current representative only on multiply steps. Square steps do not compute it. This is the same orbit function as “move to the canonical representative, square or multiply by `γ` or `z`”, because squaring commutes with `φ`:

```text
φ(X)^2 = φ(X^2),
```

and because

```text
φ^{k*}(X · J^{p^{−k*}}) = φ^{k*}(X) · J
```

when `J ∈ {γ, z}`. Carrying a non-canonical representative through a square streak does not change the orbit transition. The twelve conjugates of `γ` and of `z` are precomputed; a multiply loads one of them.

The square map on an odd-order group is bijective, and it commutes with `φ`, so it is also bijective on orbits. Collisions come from the multiply branches. The fraction of multiply steps is modeled as `1/2` because the low bit of `ρ` is modeled as a fair coin along the walk. That fairness is a heuristic. The ledger charges exactly `⌈T/2⌉` multiply steps. A square-heavier split (one multiply in three) collides too slowly on the toy group below to support a `0.39` success bound at this `T`, which is why the multiply fraction is one half and not one third.

## Collision equation

Distinguished orbits are those whose canonical representative has the low 16 bits of `ρ` equal to zero. `ρ` is orbit-invariant, so the predicate does not depend on the representative. On a distinguished orbit the implementation materialises the canonical element (at most two blst Frobenius maps of exponent at most 3, plus one conjugation) and stores it under that element, together with exponents multiplied by `p^{k*}`. The stored pair `(A, B)` then satisfies the same relation for the canonical element.

A later visit to the same canonical element with a stored pair `(A', B') ≠ (A, B)` gives, as elements of `GT`,

```text
γ^A · z^B = γ^{A'} · z^{B'}.
```

With `z = γ^a` this is one linear equation

```text
A + a B ≡ A' + a B' (mod r),
a (B − B') ≡ A' − A (mod r).
```

If `B ≠ B'`, then `a = (A' − A) · (B − B')^{−1} (mod r)`. If `B = B'`, the displayed element equality forces `A = A'`, so the exponent pairs are not distinct. A pair with `B = B'` and `A ≠ A'` would require `γ^{A − A'} = 1` with `A ≠ A'`, which does not occur. Therefore every collision of equal canonical elements with unequal exponent pairs has a nonzero coefficient, and the recovered `a` is the discrete log. The twist `p^{k*}` is already inside the stored exponents. There is no second Frobenius factor left outside the equation. A repeated identical exponent state is ignored.

Set `Q = [a]Q0`, or `Q = ∞` when `a = 0`. Check the profile acceptance conditions and `e(P, Q) = z`. A nonsingular collision passes that check because both sides equal `γ^a`. If no nonsingular collision appears within the step bound, the algorithm stops and fails.

## Step bound

Let `r` be the profile subgroup order. Use `π < π̂` with the rational

```text
π̂ = 314159265359 / 100000000000.
```

`N = (r − 1) / 12`. The step budget `T` is the least integer satisfying `T^2 ≥ π̂ · N / 2`, namely

```text
T = 82848295153133832886119523963293926164.
```

Equivalently, `T` is the least positive integer with

```text
T^2 · 200000000000 ≥ 314159265359 · N
```

and `(T − 1)` fails that inequality. `T` is even, so the model charges

```text
N_mul = T / 2 = 41424147576566916443059761981646963082
N_sqr = T / 2 = 41424147576566916443059761981646963082.
```

The distinguished-point delay of `2^{16}` steps after the orbits merge is not added to `T`. It is negligible next to `T` and is absorbed by the success slack below.

## Instruction prices

`C = 20000`. A dense Fp12 multiplication is one reference operation and is not unpacked. The pinned measurement gives `blst_fp_mul` at 383 scalar instructions. A base-field addition or subtraction is charged at 48 instructions: six 64-bit limbs, loads, add-with-carry, a conditional subtract of `p`, and stores.

One Granger–Scott cyclotomic square does 18 multiplications in `Fp` and fewer than 60 additions or subtractions in `Fp`. Charge

```text
18 · 383 + 60 · 48 + 1536 = 11310
```

instructions. The 1536-instruction term covers temporary copies (24 temporary `Fp2` values, two limbs, six words, four address or move instructions per word), the doubling of the two 255-bit exponents, and a mask test. Every step also charges 64 instructions to load `ρ` and test its low bits. That 64 duplicates the mask already inside 11310 on square steps. The duplicate is slack in the bound, not a second algorithm.

A multiply step charges one reference multiplication, plus the canonical-index selection

```text
11 · 383 + 768 = 4981
```

instructions: eleven real-limb products, and 768 instructions for the twelve-way minimum over 6-limb residues, the jump-pointer select, and two 255-bit additions of a precomputed exponent offset. The dense multiply is the reference operation. It is not also counted as 22996 instructions.

A distinguished-point materialisation charges two Frobenius maps and 4096 instructions for two 255-bit modular multiplications by the public constant `p^{k*} mod r` and for the dictionary probe. One map, from the pinned `frobenius_map_fp12` body, is at most five `Fp2` multiplications and four `Fp` multiplications. An `Fp2` multiplication is charged as three `Fp` multiplications and six `Fp` additions:

```text
5 · (3 · 383 + 6 · 48) + 4 · 383 = 8717
```

per map, so a probe costs

```text
2 · 8717 + 4096 = 21530
```

instructions. The number of probes charged is `⌈T / 2^{16}⌉`, which overcounts because only a `2^{−16}` fraction of orbits are distinguished:

```text
N_dp = 1264164659929410291841423400318817.
```

## Ledger

| Piece | Count | Unit |
| --- | ---: | --- |
| Cyclotomic squares | `N_sqr` | 11310 instructions |
| Branch tests | `T` | 64 instructions |
| Canonical-index selections | `N_mul` | 4981 instructions |
| Dense multiplies | `N_mul` | 1 reference operation |
| Distinguished-point probes | `N_dp` | 21530 instructions |
| Build `γ = e(P, Q0)` | `2^{20}` | reference operations |
| Verify pairing and subgroup checks | `2^{20}` | reference operations |
| Input decoding | `2^{16}` | reference operations |
| `G2` scalar multiplication | 512 | reference operations |
| Jump-table Frobenius images | 128 | 8717 instructions |

Jump-table charge: twelve conjugates of `γ` and twelve of `z`, each materialised by at most two maps, and the constant 128 is above `24 · 2`. Setup pairings are slack above a Miller loop; the `G2` scalar is slack above a 255-bit windowed multiply.

Instruction totals:

```text
squares, branches and selections
    = 680143079059652201078598231976661486843358
distinguished-point probes
    = 27217465128280203583345845808864130010
jump table
    = 1115776
```

Their sum is the non-reference instruction count

```text
R = 680170296524780481282181577822470352089144.
```

Reference-operation count:

```text
H = N_mul + 2^{21} + 2^{16} + 512
  = 41424147576566916443059761981649126282.
```

Total instruction-equivalents `H · C + R`:

```text
1508653248056118810143376817455452877729144.
```

Divide by `C = 20000`. The base-2 logarithm of that quotient is strictly less than `125.826529` and strictly greater than `125.826528`, so the submitted score is the ceiling at six decimal places:

```text
time_log2 = 125826529 / 1000000.
```

Preprocessing is the two `2^{20}` blocks, the `2^{16}` input block, the 512-multiplication `G2` bound, and the jump table. Its instruction total is `43265115776`, and

```text
log2(43265115776 / 20000) < 21.044773.
```

Preprocessing sits inside `H` and `R`. Nonuniform advice is empty. Every jump conjugate is computed from `P` and `z` inside the budget. `log2(max(1, 0)) = 0`.

## Memory

Each distinguished point stores one canonical `GT` element (12 · 48 bytes), two 32-byte exponents, and padding to 1024 bytes. Peak storage is

```text
N_dp · 1024 + 2^{20} = 1294504611767716138845617561927517184
```

bytes, whose base-2 logarithm is strictly less than `119.961812`. The `2^{20}` term covers code, the 24 jump conjugates, field scratch, and the output point.

## Success hypothesis

Model the orbit transition, conditional on the walk not having collided yet, as uniform in the set of `N` non-identity orbits and independent of the earlier orbit values. Under that model the probability that the first `T` orbits are all distinct is at most `exp(−T(T − 1)/(2N))`. For this `T` the exponent is `π/4` up to the rational truncation of `π`, and

```text
1 − exp(−π/4) > 0.54.
```

The submitted success lower bound is `0.39`, not `0.54`. A walk whose collision time is longer than the random-mapping time by a factor up to about `1.26` still lands at or above `0.39` inside the same `T`. The negligible real-limb tie abort and the identity-orbit accident sit inside that slack. So does the additive `2^{16}` distinguished-point delay.

The model is a heuristic. Square steps are bijections on orbits. The multiply half is what identifies distinct orbits, as in Pollard’s original mixed walk. The low bit of `ρ` is doing the work of Pollard’s partition predicate, with the predicate required to be constant on Frobenius orbits so the state space is the quotient by `⟨φ⟩`. Teske’s adding walks are the support for a fixed jump table (here, the Frobenius conjugates of `γ` and of `z`) instead of a single multiplier. The extrapolation that is not covered by those citations is the fairness of one tower limb on `GT`.

The toy program `experiment.py` runs the canonical-orbit form of this transition on the prime field of order 10009, so `12` divides the order minus one, for 200 secrets, with the step bound `ceil(sqrt(π · orbits / 2))`. The branch bit is a fixed hash of the canonical residue, the toy stand-in for the low bit of `ρ`. It checks two facts at that scale only: a collision never returns a wrong logarithm, and a singular denominator is counted rather than trusted. On this prime the half-multiply rule lands near the random-mapping hit rate, and a one-third multiply rule does not. That is why the production fraction is one half. The toy does not measure `ρ` on `GT` and it does not prove the production probability. Fresh algorithm coins do not appear: the start state is the canonical image of `z`, and the only randomness is the challenge `(P, z)`.

If no nonsingular collision appears by step `T`, the algorithm fails. That failure is what the `0.39` bound covers. There is no second independent walk and no hidden retry budget.

## What is not claimed

The attack does not invert a Miller loop, does not use a compressed cyclotomic representation, and does not claim that canonicalisation is free. Full-element Frobenius runs only while building the jump table and while storing a distinguished orbit. The hot path compares eleven real limbs of one `Fp2` slot, not twelve full field elements. A multiply fraction below one half is not claimed: the toy collision rate does not support it at this step bound.

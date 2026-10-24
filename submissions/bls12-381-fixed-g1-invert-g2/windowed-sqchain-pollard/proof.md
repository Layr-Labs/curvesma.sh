# Windowed squaring-backbone Pollard walk

Parent tree `6423c2728dc03d5f2d647c7f22adf1dac0643a66`. This package replaces the unscored template for `bls12-381-fixed-g1-invert-g2`. The score is a classical cold attack on one pair `(P, z)`. Nothing is amortized over a second target.

## Reduction

The profile fixes a nonidentity `P ∈ G1` and a uniform `z ∈ GT`, and asks for `Q ∈ G2` with `e(P, Q) = z` under the pinned `py_ecc` argument order. Let `Q0` be the profile generator of `G2` and set

```text
γ = e(P, Q0).
```

`P = [b]P0` for some `b ∈ {1, …, r−1}`, so `γ = e(P0, Q0)^b`. The pinned pairing is nondegenerate and `r` is prime, so `γ` generates `GT`. The unknown `a ∈ {0, …, r−1}` defined by `z = γ^a` satisfies `e(P, [a]Q0) = γ^a = z`. The identity `z = 1` forces `a = 0` and `Q = ∞`, which the output encoding allows. The rest of the attack is a discrete logarithm of `z` to base `γ`, followed by one scalar multiplication in `G2` and one verification pairing.

## Walk

Work in the exponent of `γ`. A state `(A, B)` denotes the element `X = γ^A · z^B = γ^{A + a B}`. Integers `A, B` are reduced modulo `r`.

The jump table is a window on a shared squaring backbone. Compute

```text
S_0 = γ,  S_{i+1} = S_i^2  for i = 0,1,2
Z_0 = z,  Z_{i+1} = Z_i^2  for i = 0,1,2
```

by Granger–Scott cyclotomic squaring (three squares on each base, six in all; the ledger charges eight). The sixteen jumps are

```text
J_{u,v} = S_u · Z_v,   u,v ∈ {0,1,2,3},
```

with `S` indexed so that `S_u = γ^{2^u}` and `Z_v = z^{2^v}`. Products that are not already on the backbone are dense Fp12 multiplications. There are at most 48 such multiplications (the bound used below; the sixteen products of two nontrivial backbone elements are fewer). Store the exponent pairs `(2^u, 2^v)`.

The online schedule is deterministic in the step index `t = 0, 1, …`:

- If `t ≡ 0 (mod 3)`, replace `(A, B) ← (2A, 2B)` and `X ← X^2` by one cyclotomic square.
- Otherwise let `h` be the integer in `{0, …, 15}` read from the lowest four bits of the canonical little-endian integer of the degree-0 coefficient of `X`. Replace `(A, B) ← (A + 2^u, B + 2^v)` and `X ← X · J_{u,v}` for the jump indexed by `h`.

The degree-0 coefficient is a function of `X` alone, so the multiplicative branch is a function of the current group element. The square steps do not consult `h`. Starting state: `(A, B) = (0, 1)`, hence `X = z`. That state is known without a discrete log.

## Collision equation

Store `(A, B)` under a key for `X` the first time `X` is seen. A later visit to the same `X` with a stored pair `(A', B')` gives

```text
A + a B ≡ A' + a B'  (mod r),
a (B − B') ≡ A' − A  (mod r).
```

If `B ≠ B'`, then `a = (A' − A) · (B − B')^{−1} (mod r)`. If `B = B'`, the same element equality forces `A = A'` whenever it is a real collision of group elements, which is a repeated exponent state and not a logarithm; the walk continues. A false branch `B = B'` and `A ≠ A'` would require `γ^{A−A'} = 1` with `A ≠ A'`, which does not happen. So a nonsingular collision returns the correct `a`, and a singular collision is discarded.

Set `Q = [a]Q0`, or `Q = ∞` when `a = 0`. Check the profile acceptance conditions and `e(P, Q) = z`. On a failed check, continue the walk. The check can fail only if the collision algebra above is violated, which it is not for a nonsingular collision, or if the input was rejected by the profile rules before the walk starts.

## Step bound

Let `r` be the profile subgroup order. Use the rational upper bound `π < 314159265359/100000000000`. The step budget is the least integer `T` satisfying

```text
T^2 ≥ (28/25)^2 · π · r / 2
```

under that upper bound on `π`. Explicitly

```text
T = 321434302617545953358844800811512932106.
```

Among these steps, those with `t ≡ 0 (mod 3)` are squares:

```text
N_sqr = (T + 2) / 3 = 107144767539181984452948266937170977369
N_mul = T − N_sqr = 214289535078363968905896533874341954737.
```

The `28/25` factor is deliberate slack over `ceil(sqrt(π r / 2))`. It is not a second algorithm. It pays for the gap, visible on prime-order toy groups, between a pure random function and this algebraic schedule.

## Instruction prices

`C = 20000`. A dense Fp12 multiplication is one reference operation. It is not also unpacked into base-field instructions. Measured `blst_fp_mul` costs 383 scalar instructions; every base-field multiplication inside a cyclotomic square is charged at 383. A base-field addition or subtraction is charged at 48 instructions: six 64-bit limbs, six loads, six adds-with-carry, six stores, a conditional subtract of `p` (six loads, six subtracts-with-borrow, six selects), and 12 instructions of addressing. That matches the blst six-limb layout and is an upper bound on the arithmetic, not a cycle count.

One Granger–Scott cyclotomic square in the degree-12 cyclotomic subgroup performs 18 multiplications in `Fp` and fewer than 60 additions or subtractions in `Fp`. Charge

```text
18 · 383 + 60 · 48 + 1536 = 6894 + 2880 + 1536 = 11310
```

instructions. The 1536-instruction term covers inter-operation copies (24 temporary `Fp2` values, two coefficients, six limbs, four address or move instructions per limb), the doubling of the two 255-bit exponents, and the distinguished-point mask test on that step.

Each multiplicative step charges one reference multiplication plus 512 instructions for the 4-bit coefficient extract, two 255-bit exponent additions, the distinguished-point mask test, and pointer updates:

```text
20000 + 512 = 20512
```

instruction-equivalents.

Distinguished points use mask `2^16`. A step whose low 16 bits of the same degree-0 coefficient are zero is stored. Number of probes:

```text
N_dp = ceil(T / 2^16) = 4904698221092925313703076184257705,
```

each charged at 2048 instructions (hash of a 12-coefficient key and one table probe). This overcounts, because only a `2^{−16}` fraction of steps are stored; charging every slot up to `ceil(T/2^16)` is the bound.

## Ledger

| Piece | Instruction-equivalents |
| --- | ---: |
| Squares `N_sqr · 11310` | included in the total |
| Multiplies `N_mul · 20512` | included in the total |
| Pairing setup `2^20` reference muls | `20971520000000` |
| Verify pairing and input checks `2^20` reference muls | `20971520000000` |
| `G2` scalar multiplication, 512 reference muls | `10240000` |
| Jump-table dense products, 48 reference muls | `960000` |
| Jump-table squares, `8 · 11310` | `90480` |
| Distinguished-point probes `N_dp · 2048` | `10044821956798311042463900025359779840` |

Setup uses `2^20` reference multiplications for `γ = e(P, Q0)` and another `2^20` for the final pairing and the subgroup checks. A production ate pairing is far smaller; the power of two is slack, not a Miller-loop trace. The `G2` scalar multiplication is charged as 512 dense Fp12 multiplications, above a 255-bit windowed `G2` multiply.

Total instruction-equivalents:

```text
5607324309217506772671637065789973243719054.
```

Divide by `C = 20000`. The base-2 logarithm of that quotient is strictly less than `127.720581`, so the submitted score is the ceiling at six decimal places:

```text
time_log2 = 127720581 / 1000000.
```

Preprocessing is the two `2^20` blocks, the 512-multiplication `G2` bound, and the jump table. Its instruction total is `41954330480`, and

```text
log2(41954330480 / 20000) < 21.000389.
```

Preprocessing sits inside the total above. Nonuniform advice is empty: every jump is computed from `P` and `z` inside the budget. `log2(max(1, 0)) = 0`.

## Memory

Each distinguished point stores one GT element (12 · 48 bytes), two 32-byte exponents, and padding to 1024 bytes. Peak storage is

```text
N_dp · 1024 + 2^20 = 5022410978399155521231950012680938496
```

bytes, whose base-2 logarithm is strictly less than `121.917792`. The `2^20` term covers code, the jump table, field scratch, and the output point. No separate advice table is retained.

## Success hypothesis

The schedule is a deterministic function of the current element and the public step index. Model the image of each step, conditional on the walk not having collided yet, as uniform in `GT` and independent of the previous values. Under that model the probability that the first `T` images are all distinct is at most `exp(−T(T−1)/(2r))`. For this `T` the exponent is about `0.985` and the model success probability is about `0.627`.

The model is a heuristic. The square every third step is a bijection on `GT` (`r` is odd), and the multiplicative branches are multiplications by fixed jumps, so the walk is not a uniformly random function. Pollard’s original mixed square-and-multiply partition is the classical support for treating the collision time as `Θ(sqrt(r))`; Teske’s adding walks are the support for using a table of fixed jumps rather than a single multiplier. Both are extrapolations. The toy program `experiment.py` runs this exact schedule, including the `28/25` inflation, on prime-order groups and checks two facts at toy scale only: a collision never returns a wrong logarithm, and the hit rate on the declared input is the number printed by that run. It does not measure Fp12 instruction counts and it does not prove the production probability.

The submitted success lower bound is `0.39`, not the model figure `0.627`. The `28/25` slack is the allowance for the algebraic walk colliding more slowly than a random function, as the toy rates do. A singular denominator does not produce a wrong `Q`; it produces no candidate, and the walk continues inside the same `T` steps. Fresh algorithm coins are the coins that would randomize a start state. This walk’s start state is the fixed pair `(0, 1)`. The only randomness in the algorithm is the challenge itself. The hypothesis is over uniform `(P, z)`, with the walk deterministic given the challenge.

If no nonsingular collision appears by step `T`, the algorithm stops and fails. That failure is what the `0.39` bound covers. There is no second independent walk and no hidden retry budget.

## What is not claimed

The attack does not walk on Frobenius orbits, does not invert a Miller loop, and does not use a compressed-squaring orbit representative. Canonicalization under `z ↦ z^{p^k}` is not part of the step. The jump table is closed under a squaring backbone of `γ` and `z` only.

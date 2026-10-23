# Frobenius-scaled mixed add/square rho on GT

Parent crown `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

The attack solves fixed-G1 pairing inversion. `P` is the supplied non-identity point of `G1` and `z` is the supplied element of `GT`. The output is `Q` in `G2` with `e(P, Q) = z`. The score below is `log2` of total work in dense `Fp12`-multiplication equivalents, `C = 20000`.

## Claimed bound

Let `r` be the BLS12-381 subgroup order and `n = (r - 1) / 12`. Let

```
W0 = (51050 / 60000) * sqrt(π n / 2).
```

The walk derived below performs at most `W0` reference equivalents in expectation under the random-mapping heuristic, before a `2^20` tail that covers the pairing, the seed, the final `G2` scalar, and one orbit check. Hence

```
log2(W0 + 2^20) < 12575/100.
```

The algorithm that stops at the birthday time `t = sqrt(π n / 2)` succeeds with probability at least `1/2` under that heuristic. Peak storage is at most `2^127` bytes. Preprocessing is at most `2^20` reference equivalents and is inside the same bound. Nonuniform advice is empty.

## Reduction

`Q0` is the public generator in the track profile. Compute `g = e(P, Q0)` with the pinned pairing (library order `(Q, P)`). If `P ≠ O` then `g` generates `GT[r]`, because `e` is non-degenerate and both groups are cyclic of prime order `r`. If `z = 1`, return `Q = O`. Otherwise solve `z = g^y` and return `Q = [y] Q0`. Bilinearity gives `e(P, [y] Q0) = g^y = z`.

`r` does not divide `p - 1`, and the multiplicative order of `p` modulo `r` is `12`. Frobenius `φ(X) = X^p` is an automorphism of `GT[r]` of order `12`. For `X ≠ 1` the orbit `{φ^k(X)}` has size `12`. Also `X^{p^6} = X^{-1}` on this subgroup, and `p^6`-power is the blst conjugation `neg` of the `w`-half.

## Tower and the limb invariant

Work in the tower measured by the reference multiplication, which is the pinned blst layout

```
Fp2 = Fp[u]/(u^2 + 1),  Fp6 = Fp2[v]/(v^3 - (u + 1)),  Fp12 = Fp6[w]/(w^2 - v).
```

Write `X` as `a[0] + a[1] w` with `a[j]` three `Fp2` values. Let `λ0` be the `Fp` limb `a[0][2][0]` and `λ1` the limb `a[0][2][1]`.

In `frobenius_map_fp12` of the pinned `fp12_tower.c`, `frobenius_map_fp6` copies that `Fp2` and multiplies both limbs by the `Fp` constant `coeffs2[n - 1]`. Those constants demontgomery to

```
c1 = (1 + u)^((2p - 2) / 3),  c2 = c1^2,  c3 = c1^3,
```

each lying in `Fp`, and `c1^6 = 1` in `Fp` while `c1^3 ≠ 1`. Applying the `n = 1` map therefore sends `(λ0, λ1)` to `(c1 λ0, c1 λ1)`. Iterating, `φ^k` multiplies both limbs by `c1^k`. Conjugation `φ^6` fixes `a[0]`, and `c1^6 = 1`, so the same formula holds for every `k`.

Define the integer

```
κ(X) = λ0(X)^6 + 2^192 · λ1(X)^6.
```

Then `κ(φ^k(X)) = κ(X)`. A distinguished point is a state with `κ(X) ≡ 0 (mod 458)`. The predicate is constant on orbits.

`λ0` and `λ1` are not the zero function on `GT[r]`. The identity has `a[0] = 1` and `a[1] = 0`, so those limbs vanish there, but the norm equation `a[0]^2 - v a[1]^2 = 1` of a cyclotomic element does not force the `v^2` coefficient of `a[0]` to vanish. The walk below treats a zero limb as an ordinary value. Both limbs are folded into `κ`, so a single vanished limb does not collapse the key.

Heuristic, stated as such: as `X` runs through uniform elements of `GT[r]`, `(λ0, λ1)` is close enough to uniform on `Fp^2` that distinct orbits collide in `κ` with probability at most `2^-128` per pair. The walk visits `t ≈ 2^125.96` states and therefore sees an expected `≪ 2^-10` false key matches. Each false match is rejected by the orbit test below and does not change the leading term.

## State and branches

A state is `(X, α, β)` with `X = g^α · z^β` in `GT[r]`. The starting state is `(g^s, s, 0)` for a uniform seed `s ∈ {1, …, r - 1}`. One fixed multiplier `m = g · z` is computed once, with exponent increment `(1, 1)`.

The branch is `λ0(X) mod 3`.

- `0`: replace `X` by its cyclotomic square and `(α, β)` by `(2α, 2β)`.
- `1`: cyclotomic square, then conjugate. On this subgroup the conjugate of `X^2` is `X^{-2}`, so the exponents become `(-2α, -2β)`.
- `2`: replace `X` by `X · m` and `(α, β)` by `(α + 1, β + 1)`.

Each branch preserves the invariant `X = g^α · z^β`. Each individual branch is a bijection of `GT[r]` (`2` is invertible modulo `r`, and multiplication by `m` is a translation). The partitioned map is not a bijection: images of different branches overlap. That is the same shape as Pollard's three-branch walk on a prime-order group.

Heuristic, stated as such: this map behaves as a random function on a set of `n = (r - 1) / 12` orbits for the purpose of the first-collision time. The classical birthday integral then says the first time the walk enters an orbit already seen is at most `t = sqrt(π n / 2)` steps in expectation, and

```
1 - exp(-t^2 / (2n)) = 1 - exp(-π/4) > 1/2.
```

Distinguished points only affect storage. The walk still performs the steps; detection waits at most `458` extra steps after the orbits meet, which is absorbed in the `2^20` tail.

## Collision equation

Store a distinguished record `(X, α, β)`. On a later distinguished state with the same `κ`, apply `φ` (the `n = 1` map) up to `11` times to the stored element and compare limbs with the new element.

If no power matches, the key collision is false. Continue.

If `φ^k(X1) = X2`, then `X2 = X1^(p^k) = g^(α1 p^k) z^(β1 p^k)` and also `X2 = g^(α2) z^(β2)`. With `z = g^y`,

```
α1 p^k + y β1 p^k ≡ α2 + y β2 (mod r),
y ≡ (α1 p^k - α2) · (β2 - β1 p^k)^{-1} (mod r),
```

using the twelve precomputed residues `p^0, …, p^11 mod r`. Let `D = β2 - β1 p^k`. If `D ≠ 0`, the display is `y`.

If `D = 0`, the same identity forces `α2 = α1 p^k`. The two exponent vectors are one Frobenius rescaling of a single relation: the hit is a revisit of one orbit, not a second path. Discard it. A genuine meeting of two paths has `D ≠ 0`, because a second path that is not that rescaling cannot satisfy both `D = 0` and the group equality.

## Cyclotomic square, operation count

The square is the Granger–Scott schedule in the pinned source: three `sqr_fp4` calls plus the six limb reconstructions in `cyclotomic_sqr_fp12`. Each `Fp` multiplication is charged at the pinned measurement of `blst_fp_mul`, `383` instructions. Each `Fp` addition, subtraction, or negation is charged at `48` instructions, which covers the following straight-line limb program and a margin for address arithmetic. Limbs are 6 words of 64 bits, inputs in `[0, p)`.

- load both operands: 6 `ldp`
- add with carry: 6 `adc`
- load `p`: 3 `ldp`
- subtract with borrow: 6 `sbc`
- select the reduced residue: 6 `csel`
- store: 3 `stp`

That is 30 instructions. The charge `48` leaves 18 instructions for pointer setup and the final flag.

`Fp2` square used here, the reduced Karatsuba square:

```
(a0 + a1 u)^2 = (a0 + a1)(a0 - a1) + (2 a0 a1) u
```

is 2 `Fp` multiplications and 3 `Fp` additions.

`sqr_fp4(a0, a1)`, matching the pinned fallback body:

- two `Fp2` squares: 4 mul, 6 add
- one `Fp2` add of the inputs: 2 add
- one multiplication by `1 + u`, which is `(x0 - x1, x0 + x1)`: 2 add
- one `Fp2` add: 2 add
- one `Fp2` square: 2 mul, 3 add
- two `Fp2` subtractions: 4 add

Total per `sqr_fp4`: 6 `Fp` multiplications and 19 `Fp` additions.

The outer reconstruction in `cyclotomic_sqr_fp12` is three blocks of (sub, double, add) on `a[0]`, one multiplication by `1 + u`, and three blocks of (add, double, add) on `a[1]`:

```
3 · 6 + 2 + 3 · 6 = 38
```

`Fp` additions. Three `sqr_fp4` contribute `57` additions. The cyclotomic square is

```
18 Fp multiplications and 95 Fp additions,
18 · 383 + 95 · 48 = 6894 + 4560 = 11454
```

instructions. No inversion is required. Conjugation is `neg` of the six `Fp` limbs in `a[1]`: `6 · 48 = 288` instructions.

## Key and control

Each sixth power is three `Fp` multiplications (`λ^2`, `λ^4`, `λ^4 · λ^2`). Two limbs cost `6 · 383 = 2298` instructions. The shift by `2^192` is a limb placement. Control, the `mod 3` branch, a 4-limb exponent update modulo `r`, and the distinguished test are charged at `320` instructions every step. That budget is `256` plus `64` of margin; a 4-limb add-with-carry and a compare do not fill it.

## Hot-step total

Reference dense multiplication costs `1`, i.e. `20000` instruction-equivalents, and is not also expanded. Instruction charges:

| branch | instructions |
| --- | ---: |
| square, `λ0 ≡ 0` | `11454 + 2298 + 320 = 14072` |
| square-conjugate, `λ0 ≡ 1` | `14072 + 288 = 14360` |
| multiply by `m`, `λ0 ≡ 2` | `20000 + 2298 + 320 = 22618` |

The three branches are equally likely under the uniform-limb heuristic already stated (`λ0 mod 3` is balanced). Mean charge:

```
(14072 + 14360 + 22618) / 3 = 51050 / 3
```

instruction-equivalents, or `51050 / 60000` reference multiplications per step.

## Tail

One pairing, the seed exponentiation `g^s`, the product `m = g · z`, the final scalar `[y] Q0`, one inversion modulo `r`, and at most twelve Frobenius maps on a single collision are together far below `2^20` reference multiplications. A pairing in this tower is a few thousand `Fp12` operations; a 255-bit scalar on `G2` is a few hundred. The charge `2^20` is an explicit upper bound, not a measurement, and it is loose by orders of magnitude. It is preprocessing plus postprocessing, included in the total.

`W0 > 2^120`, so `2^20 / W0 < 2^-100` and

```
log2(W0 + 2^20) < log2(W0) + 2^-99.
```

`log2(W0) = log2(51050 / 60000) + (1/2) log2(π n / 2)` evaluates to less than `125.729`. Adding `2^-99` stays under `12575/100`.

## Memory

Expected distinguished records before detection: `t / 458`. Each record stores the canonical `Fp12` encoding of `X` (`576` bytes) and two `32`-byte exponents (`64` bytes): `640` bytes. Then

```
(t / 458) · 640 < 2^127.
```

Nothing else is retained. The Frobenius constants and `p^k mod r` are public profile constants, not advice. Peak memory, including that table, is at most `2^127` bytes.

## What is heuristic

Two points are heuristics, both standard for a generic rho and both required for the leading term. First, the partitioned map collides on orbits at the birthday time `sqrt(π n / 2)`. Second, the limb pair `(λ0, λ1)` does not collapse `κ` on a positive-density set of orbits. The reduction, the exponent invariant, the nonzero-denominator claim, and the instruction arithmetic are unconditional. A failed orbit check is a continue, not a wrong `Q`. The returned point is accepted only when the exponent equation produces `y` and, by the reduction, `e(P, [y] Q0) = z`.

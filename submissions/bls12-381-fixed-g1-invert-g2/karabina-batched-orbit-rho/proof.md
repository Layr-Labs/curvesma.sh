# Karabina compressed-square rho with batched decompression

Parent crown `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

The attack solves fixed-G1 pairing inversion. `P` is the supplied non-identity point of `G1` and `z` is the supplied element of `GT`. The output is `Q` in `G2` with `e(P, Q) = z`. The score is `log2` of total work in dense `Fp12`-multiplication equivalents. The reference price is `C = 20000` scalar AArch64 instruction equivalents, the organizer normalization of the pinned `blst_fp12_mul` measurement (22,996 instructions at commit `f7a88067ed98385ebc6cf41ccecf21ead706a5b7`). One call of that dense multiplication is charged as `1`, and is not also expanded. Every other kernel is charged as its instruction count divided by `C`.

## Claimed bound

Let `r` be the BLS12-381 subgroup order and `n = (r - 1) / 12`. Let `S` be the integer defined in the certificate below, so `S ≥ sqrt(π n / 2)`. Let

```
W0 = (50400 / 60000) * S.
```

`50400` is an instruction sum over the three branches of one step, derived below, and `60000 = 3 · C`. The walk performs at most `W0` reference equivalents in expectation under the random-mapping heuristic, before a `2^20` tail that covers the pairing, the seed, the final `G2` scalar, batched-inversion startup, and one orbit check. Hence

```
log2(W0 + 2^20) < 12572/100.
```

Stopping at the birthday length `S` succeeds with probability at least `1/2` under that heuristic. Peak storage is at most `2^127` bytes. Preprocessing is at most `2^20` reference equivalents and sits inside the same bound. Nonuniform advice is empty.

## Reduction

`Q0` is the public generator in the track profile. Compute `g = e(P, Q0)` with the pinned pairing, whose library order is `(Q, P)`. If `P ≠ O` then `g` generates `GT[r]`, because `e` is non-degenerate and both groups are cyclic of prime order `r`. If `z = 1`, return `Q = O`. Otherwise solve `z = g^y` and return `Q = [y] Q0`. Bilinearity gives `e(P, [y] Q0) = g^y = z`.

`r` does not divide `p - 1`, and the multiplicative order of `p` modulo `r` is `12`. In particular `p^6 ≡ -1 (mod r)`, since `r` divides `p^4 - p^2 + 1` and therefore divides `p^6 + 1`. Frobenius `φ(X) = X^p` is an automorphism of `GT[r]` of order `12`. For `X ≠ 1` the orbit `{φ^k(X)}` has size `12`. On this subgroup `X^{p^6} = X^{-1}`, and the `p^6`-power is conjugation: it fixes the `Fp6` half `a[0]` and negates `a[1]`.

## Tower and compressed state

The reference multiplication is the pinned blst layout

```
Fp2 = Fp[u]/(u^2 + 1),  Fp6 = Fp2[v]/(v^3 - (u + 1)),  Fp12 = Fp6[w]/(w^2 - v).
```

Write `X = a[0] + a[1] w` with each half three `Fp2` values. Name

```
g1 = a[0][1],  g2 = a[0][2],  g3 = a[1][0],  g5 = a[1][2].
```

The walk stores only `(g1, g2, g3, g5)`. These four coefficients are the input and the output of the compressed square below, so a square does not rebuild `a[0][0]` or `a[1][1]`.

The straight-line program is Karabina's compressed cyclotomic square, in the blst naming above. `nr(t) = (t0 - t1, t0 + t1)` is multiplication by `1 + u`, two `Fp` adds and no multiplication. `sq2` is the reduced Karatsuba square `(a0 + a1 u)^2 = (a0 + a1)(a0 - a1) + (2 a0 a1) u`, which is 2 `Fp` multiplications and 3 `Fp` additions. Registers are overwritten exactly as in this sequence:

```
t0 = sq2(g1)
t1 = sq2(g5)
t5 = g1 + g5
t2 = sq2(t5)
t3 = t0 + t1
t5 = t2 - t3
t6 = g3 + g2
t3 = sq2(t6)
t2 = sq2(g3)
t6 = nr(t5)
t5 = 2 (t6 + g3)
z3 = t5 + t6
t4 = nr(t1)
t5 = t0 + t4
t6 = t5 - g2
t1 = sq2(g2)
t6 = 2 t6
z2 = t6 + t5
t4 = nr(t1)
t5 = t2 + t4
t6 = t5 - g1
t6 = 2 t6
z1 = t6 + t5
t0 = t2 + t1
t5 = t3 - t0
t6 = t5 + g5
t6 = 2 t6
z5 = t5 + t6
```

The output tuple is `(z1, z2, z3, z5)`. The six `Fp2` squares are the six `sq2` calls. That is `12` `Fp` multiplications. The surrounding `Fp` additions, subtractions, doublings, and three `nr` maps are `64` `Fp` operations: `18` inside the six squares and `46` outside them (`12` binary adds, `4` subtracts, and `4` doublings, each two limbs, plus `6` limbs from the three `nr` maps). Each `Fp` multiplication is charged at the pinned `blst_fp_mul` measurement, `383` instructions. Each `Fp` addition, subtraction, negation, or doubling is charged at `48` instructions. That figure covers a 6-limb AArch64 program (loads of both operands, `adds`/`adcs`, a compare against `p`, `csel`, stores) and leaves margin for address arithmetic. The compressed square is

```
12 · 383 + 64 · 48 = 4596 + 3072 = 7668
```

instructions. No inversion is required. The program was checked against a full tower square on the subgroup of order dividing `p^4 - p^2 + 1` in a toy field of 19 elements; that check is only evidence at that scope. The same straight-line program is what the attack runs over BLS12-381.

Conjugation negates `a[1]` and fixes `a[0]`. On the stored tuple it negates `g3` and `g5` and leaves `g1` and `g2`. That is four `Fp` negations, `192` instructions, and it stays in the compressed representation. The missing `a[1][1]` would be negated as well; it is not stored, and the next compressed square does not read it.

## State and branches

A state is a compressed element together with exponents `(α, β)` such that the underlying group element is `X = g^α · z^β`. The starting state is the compression of `g^s` with exponents `(s, 0)`, for a uniform seed `s ∈ {1, …, r - 1}`. One fixed multiplier `m = g · z` is computed once, as a full `Fp12` element, with exponent increment `(1, 1)`.

The branch is the least limb of `g2`, reduced modulo `3`. `g2` is `a[0][2]`, the coefficient whose Frobenius image is scaled by the pinned `Fp` constant `c1 = (1+u)^{(2p-2)/3}`. That constant satisfies `c1^6 = 1` and `c1^3 ≠ 1` in `Fp`. The branch test itself is a compare inside the control budget below.

- `0`: replace the tuple by its compressed square and `(α, β)` by `(2α, 2β)`.
- `1`: compressed square, then negate `g3` and `g5`. On this subgroup that is conjugation of the square, so the exponents become `(-2α, -2β)`.
- `2`: decompress, multiply the full element by `m`, discard `a[0][0]` and `a[1][1]`, and add `(1, 1)` to the exponents.

Each branch preserves `X = g^α · z^β` when decompression returns the element that was compressed. Each individual branch is a bijection of `GT[r]` (`2` and `-2` are invertible modulo `r`, and multiplication by `m` is a translation). The partitioned map is not a bijection. Square-only iteration would be a bijection, because `gcd(2, r) = 1`; the translation branch is what makes the functional graph available to a collision search.

Heuristic, stated as such. Model the partitioned map as a random function on the set of `n = (r - 1) / 12` Frobenius orbits, for the purpose of the first time the walk enters an orbit already seen. The translation branch does not commute with Frobenius, so this is an extrapolation from the equivariant branches to the whole partition, not a theorem that the map descends to the quotient. Under that model the first orbit collision has expectation at most `sqrt(π n / 2)`, and at length `S` at least that large

```
1 - exp(-S^2 / (2n)) ≥ 1 - exp(-π/4) > 1/2.
```

A second hypothesis, also stated as such: as `X` runs through uniform elements of `GT[r]`, the stored `Fp2` coefficients are not concentrated at `0`, and the least limb of `g2` is equidistributed modulo `3` and in its low `9` bits. If `g3 = 0` or `g5 = 0` on a mul branch, the denominator below vanishes and that walk is abandoned. A fresh seed replaces it. The expected number of such hits in `S` steps is at most `2S / p^2`, which is far below `2^{-600}` and does not move the leading term.

## Decompression and the batched inverse

On the common mul branch, `g5 ≠ 0` and `g3 ≠ 0`. Recover the missing coefficient `a[1][1]` by the division

```
a[1][1] = (nr(g5^2) + 3 g1^2 - 2 g2) · (4 g3)^{-1},
```

then

```
a[0][0] = 1 + nr(2 a[1][1]^2 + g3 g5 - 3 g2 g1).
```

The field work, excluding the single `Fp` inversion inside `(4 g3)^{-1}`, is `19` `Fp` multiplications (four `Fp2` squares, three Karatsuba `Fp2` products, the norm of the `Fp2` inverse, and the two scalings by that inverse) and at most `80` `Fp` additions. The `80` is a margin over the `73` additions obtained by counting each Karatsuba product as ten limb additions and each square as three.

```
19 · 383 + 80 · 48 = 7277 + 3840 = 11117
```

instructions, plus one inversion of an element of `Fp`.

The inversion is not performed inside the walker. `3072` independent walks run on one target. A walker that takes the mul branch ships the `Fp` norm of `4 g3` and waits. Inversions are issued in batches of exactly `1024` norms, splitting a fuller round into several batches. Montgomery's trick inverts one batch with one Fermat inversion and `3 · 1023` extra `Fp` multiplications. A partial batch occurs only at the end of the attack; it is charged as one extra Fermat inversion inside the tail.

Fermat inversion of one `Fp` element uses a width-`5` window on `p - 2`. The exponent has `381` bits and Hamming weight `229`, and the window performs `380` squares, `68` multiplications, and `8` precomputed odd multiples: `456` `Fp` multiplications,

```
456 · 383 = 174648
```

instructions. Per mul-branch event the batch contributes at most

```
3 · 383 + ceil(174648 / 1024) = 1149 + 171
```

instructions, and the send/receive of one norm plus its inverse is charged at `256` instructions. The sum `1600` is the per-event inversion budget, including a few instructions of margin. Stalled walkers add no arithmetic. Transfers between the `3072` walkers and the inversion batch are inside that `256`.

The dense product by `m` is one reference multiplication, cost `1`. Recompression drops two `Fp2` values and is free. Control, the modulo-`3` test, the low-bit distinguished test, and a 4-limb exponent update modulo `r` are `320` instructions on every branch.

## Hot-step total

| branch | instructions |
| --- | ---: |
| compressed square | `7668 + 320 = 7988` |
| square then conjugate | `7668 + 192 + 320 = 8180` |
| decompress, times `m`, recompress | `11117 + 1600 + 20000 + 320 = 33037` |

The three branches are equally likely under the limb heuristic. The enumerated mean sum is

```
7988 + 8180 + 33037 = 49205
```

instructions, or `49205 / 60000` reference multiplications per step. The bound charges `50400 / 60000` instead. The `1195` extra instructions on the sum cover the distinguished-point key evaluated only on hits (three `Fp` multiplications, amortized by the `2^9` gap below) and address arithmetic not listed in the limb programs.

## Distinguished points and the collision equation

A state is distinguished when the low `9` bits of the least limb of `g2` are zero. The expected gap is `512` steps. The predicate is not constant on a Frobenius orbit. That only delays detection: after two walks occupy the same orbit, a further geometric wait of mean `512` steps reaches a distinguished state. `512` steps at under one reference multiplication each, and a further `2^12` steps charged the same way, sit inside the `2^20` tail. On a hit, the stored key is `λ^6` where `λ` is the `Fp` integer `g2`'s first limb interpreted in `Fp`. Scaling by `c1` multiplies this `Fp` element by `c1`, and `c1^6 = 1`, so the key agrees on `φ^k` for every `k`. Computing it costs three `Fp` multiplications and is already inside the `1195`-instruction margin.

Store the compressed tuple and the two exponents under that key. A later hit with the same key decompresses both tuples (two events of the same batched inverse, negligible next to `S`) and applies `φ` up to `11` times. Each `φ^k` has the shape of the pinned `frobenius_map_fp12`: negations of the `u`-limbs and multiplications by precomputed tower constants. Eleven such maps, charged flatly as eleven reference multiplications, sit in the tail.

If no power matches, the key collision is false. Continue.

If `φ^k(X1) = X2`, then `X2 = X1^{p^k} = g^{α1 p^k} z^{β1 p^k}` and also `X2 = g^{α2} z^{β2}`. With `z = g^y`,

```
α1 p^k + y β1 p^k ≡ α2 + y β2 (mod r),
y ≡ (α1 p^k - α2) · (β2 - β1 p^k)^{-1} (mod r),
```

using the twelve residues `p^0, …, p^{11} mod r`, which are public consequences of the profile parameter `p` and are not advice. Let `D = β2 - β1 p^k`. If `D ≠ 0`, this `y` is the discrete log.

If `D = 0`, the same identity forces `α2 ≡ α1 p^k (mod r)`. The two stored vectors are one Frobenius rescaling of a single relation, a revisit of one orbit rather than a meeting of two paths. Discard the hit. A meeting of two distinct paths has `D ≠ 0`: a second path that is not that rescaling cannot satisfy both `D = 0` and the group equality.

The returned `Q = [y] Q0` is accepted only through this equation. A failed orbit check does not emit a point. By the reduction, a `y` from a nonzero denominator satisfies `e(P, Q) = z`.

## Certificate

Let `n = (r - 1) / 12`. Since `π < 355/113`,

```
π n / 2 < 355 n / (113 · 2) = 355 n / 226.
```

Let `S` be the least integer with `S^2 ≥ ceil(355 n / 226)`, so `S^2 · 226 ≥ 355 n` and `S ≥ sqrt(π n / 2)`. The charged work plus the tail is at most

```
Num / 60000,  Num = 50400 · S + 60000 · 2^20.
```

The integer comparison `Num^{100} < 60000^{100} · 2^{12572}` holds, and therefore

```
log2(Num / 60000) < 12572/100.
```

`experiment.py` recomputes that comparison, the constant identities `c1^6 = 1`, `c2 = c1^2`, `c3 = c1^3`, and `p^6 ≡ -1 (mod r)`, and the toy collision equation. None of those checks is a production collision count.

## Tail and memory

One pairing, the seed power `g^s`, the product `m = g · z`, the final scalar `[y] Q0`, one inversion modulo `r`, the final partial inversion batch, and at most twelve Frobenius maps plus two decompressions on a single accepted collision are together far below `2^20` reference multiplications. The ate loop length is the bit length of the profile parameter `x`, about `64` iterations. A `255`-bit scalar on `G2` is a few hundred group operations, each much smaller than a dense `Fp12` product. The charge `2^20` is an explicit upper bound. It is preprocessing plus postprocessing, included in the total.

`W0 > 2^{125}`, so the tail does not move the logarithm across the claimed rational.

Expected distinguished records are at most `S / 512`. Each record stores four `Fp2` values (`384` bytes) and two `32`-byte exponents. At `640` bytes, which is loose for the compressed record,

```
(S / 512) · 640 < 2^127.
```

The `3072` live walker states are the same shape and are negligible beside the table. The Frobenius constants and `p^k mod r` are public. Peak memory, including the batch buffer of `1024` norms, is at most `2^127` bytes.

## What is heuristic

Two points carry the leading term. First, the partitioned map collides on Frobenius orbits at the birthday time of a set of size `n`, even though multiplication by `m` is not Frobenius-equivariant. Second, the stored coefficients are not stuck at `0`, and the limb of `g2` is balanced modulo `3` and in its low `9` bits. The reduction, the exponent update, the compressed-square program, the nonzero-denominator algebra, the batch arithmetic, and the integer certificate are unconditional. Abandoning a walk whose denominator vanishes is a restart, not a wrong `Q`.

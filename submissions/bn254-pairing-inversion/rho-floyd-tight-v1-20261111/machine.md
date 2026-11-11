# Bounded scalar implementation and instruction ledger

All counts are scalar AArch64 instruction equivalents on 64-bit general-purpose
registers (arithmetic, logic, loads/stores, addressing, control), consistent with
`tracks/classical-group-work-v2.json` and `tracks/docs/COST_MODEL.md`. Analytic
instruction-equivalent bounds are sufficient under the cost model; no assembly or
production-size run is required. Every figure below is a conservative upper bound,
rounded up to a power of two for the ledger. No BN254 operation is assigned a unit
cost; the reference unit `C = 20,000` is only the final divisor.

## Field representation

`Fp` elements (254-bit `p`) are four 64-bit limbs in Montgomery form with
`R = 2^256`. `Fp2 = Fp[u]/(u^2+1)`, `Fp6 = Fp2[v]/(v^3 - (9+u))`,
`Fp12 = Fp6[w']/(w'^2 - v)`, equivalently the pinned
`Fp12 = Fp[w]/(w^12 - 18 w^6 + 82)` with `u = w^6 - 9`. A GT element is 12 `Fp`
limbvectors (`12 * 4 * 8 = 384` bytes), which is the `fixed_encoding` length.

## Fp multiplication (four-word Montgomery), <= 256 instructions

CIOS Montgomery multiplication for four limbs performs `4*4` limb products plus
`4*4` reduction products (`32` `mul`/`umulh` pairs), the interleaved limb
additions and carries, the final conditional subtraction, and loads/stores. A
standard scalar CIOS four-limb routine is well under `256` instructions; Montgomery
squaring is cheaper. We charge `Fp mul <= 256` and `Fp add/sub <= 32` (four limb
adds/subtracts with carry and a conditional modular correction). These are generous:
the pinned blst **base-field** multiply is `383` instructions for a 381-bit (six-limb)
field, so a 254-bit four-limb multiply is smaller.

## Fp12 multiplication, <= 32,768 instructions

A dense `Fp12` multiply reduces through the tower to `54` `Fp` multiplications and
`O(100)` `Fp` additions/subtractions with Karatsuba at each level (`Fp2`: 3 mul;
`Fp6`: 6 `Fp2` mul; `Fp12`: 3 `Fp6` mul), plus the sparse reductions by
`u^2 = -1`, `v^3 = 9+u` and `w'^2 = v`. A conservative straight-line bound is

```
54 * 256 (Fp mul) + 300 * 32 (Fp add/sub) + 2,000 (reductions, addressing) <= 25,500.
```

We charge **`Fp12 mul <= 32,768 = 2^15`**. This upper-bounds the organizer's pinned
measurement of the **larger** 381-bit blst dense `Fp12` multiply, `22,996`
instructions (`tracks/docs/COST_MODEL.md`); the 254-bit BN254 product is strictly
smaller, so `2^15` is safe with margin.

## Partition selector SHA-256 over 384 bytes, <= 16,384 instructions

`SHA256(salt || fixed_encoding(Y))` hashes `32 + 384 = 416` bytes, i.e. `7`
64-byte compression blocks after padding. A scalar SHA-256 compression (message
schedule plus `64` rounds) is a few hundred to `~2,000` instructions per block; at
a generous `<= 2,048` per block, `7 * 2,048 = 14,336 <= 16,384 = 2^14`. Only one
byte of the digest is consumed (`& 31`).

## Exponent update and table access, <= 4,096 instructions

Per group step: two 254-bit modular additions `a += c_j`, `b += d_j` (four-limb add
plus conditional subtract by `r`, `<= 2 * 64`), the table index from the selector,
a load of the precomputed `M_j` (12 limbvectors) and loop control: `<= 4,096 = 2^12`.

## Per group step and per Floyd iteration

```
group step  = 32,768 + 16,384 + 4,096 = 53,248.
iteration   = 3 * 53,248 (tortoise 1 + hare 2 steps)
            + 8,192 (2^13: 12-coefficient equality test + multiword counters) = 167,936.
```

The `8,192` compare/counter term covers twelve 4-limb equality comparisons and the
maintenance of a **three-word (192-bit) iteration counter** and the four-word
labels, so the loop executes the full `N ~ 2^127.8` iterations without relying on a
single 64-bit register.

## Setup, recovery, output and subgroup checks, <= 2^40 instructions

One-time: decode and canonicality checks of `z`; the subgroup test `z^r = 1`
(one `~254`-step Fp12 exponentiation, `< 254 * 2 * 32,768 < 2^25`); conversion to
Montgomery form; sampling `66` scalars by bounded rejection (`<= 66 * 1024` 254-bit
draws); the `32` multiplier products `M_j`; final modular inverse and
exponentiation `g^x` for verification (`< 2^26`). Pairing track additionally:
`Q = [x] Q0` by four-word double-and-add in `G2` (`<= 254` doublings + `<= 254`
additions, each a handful of `Fp2` operations, `< 2^24`), and the subgroup checks
`[r]P = [r]Q = infinity` (`< 2^26`). The sum is far below the charged
`SETUP = 2^40`, which absorbs all one-time work on both tracks with large margin.

## Totals

```
N        = 2 * ceil(sqrt(r)) = 295,893,513,763,578,638,011,461,384,341,992,519,220
PER_ITER = 167,936
R_total  = N * PER_ITER + 2^40
         = 49,691,173,127,400,342,153,092,779,040,857,955,219,357,696
W        = R_total / 20,000
time_log2 = ceil6(log2 W) = 130.868186.
```

## Memory

Peak storage is the salt (`32` B), the 32-entry `M_j` table (`32 * 384` B), the
tortoise and hare states, a handful of scalars and the SHA-256 state: well under
`2^20` bytes. There is no distinguished-point table, no large hash table and no
network transfer; memory has no scalar contribution or tie-break under the cost
model.

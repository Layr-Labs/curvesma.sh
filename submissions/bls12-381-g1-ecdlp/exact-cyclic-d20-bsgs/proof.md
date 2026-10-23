# Pohlig–Hellman on the smooth cofactor of the supplied 2^20-power

Parent: `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

## Claim

Let `r` be the BLS12-381 G1 order

```text
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
```

and let `P` be the standard generator `P0`. The profile gives `n`, the points
`A_i = [s^i] P` for `i = 1..n`, and `T = A_1 = [s] P`, with `s` uniform in
`1..r-1`. The bound below applies to every integer `n` with `2^20 ≤ n < r`.
It is scored at `n = 2^20`. The only power the attack reads is `A_d` with
`d = 2^20`. That index is present because `d ≤ n`. Generation of the powers
is organizer work and is not charged. The attack returns the integer `s`.

```text
W(n) = 239362939
log2(W(2^20)) ≤ 27834625/1000000 = 27.834625
```

Success probability is 1. Peak memory is at most `2^23` bytes. Preprocessing
is at most `2^27` target-group-operation equivalents and is included in `W`.
Advice is 0 bytes. There is no coverage probability, no Erdős factor, and no
unperformed lift from `s^d` to `s`.

The reference price is the track constant `C = 2000` scalar AArch64
instructions per batched incomplete affine G1 addition. A field multiplication
is the pinned `blst_fp_mul` count of 383 instructions. Other work costs
`R/C`. A reference addition is never also unpacked into its own 6M formula.

## Factorization

`r` ends in the hex digit string `00000001`, so `v_2(r-1) = 32` and
`d = 2^20` divides `r-1`. Write `m = (r-1)/d`. Then

```text
r - 1 = 2^32 · 3 · 11 · 19 · 10177 · 125527 · 859267 · 906349^2
        · 2508409 · 2529403 · 52437899 · 254760293^2
m     = 2^12 · 3 · 11 · 19 · 10177 · 125527 · 859267 · 906349^2
        · 2508409 · 2529403 · 52437899 · 254760293^2
```

The product identity is an integer equality with `r-1`. Each odd factor is
strictly less than `3,317,044,064`. For every odd integer in that range,
passing the Miller–Rabin tests to bases `2, 3, 5, 7, 11, 13, 23` is a
complete primality proof (Jaeschke 1993). Each odd factor above passes those
seven bases, so each is prime. Squares are only `906349^2` and `254760293^2`.
No other prime divides `r-1`.

Consequently every prime factor of `m` is at most `254760293 < 2^28`.
Pohlig–Hellman on a cyclic group of order `m` is a finite sum of baby-step
giant-step searches of size `ceil(sqrt(q))`, one or two per prime, not a
search of size `sqrt(m)`.

## Scalar discrete log that the point `A_d` determines

Let `g = 7`. This is a primitive root modulo `r`: for every prime `q` above,
`7^{(r-1)/q} ≢ 1 (mod r)`. The attack recomputes those twelve exponentiations
and aborts if one of them is `1`. On this `r` none of them is `1`, so the
abort does not fire. `F_r^* = <7>`.

Write `s ≡ 7^L (mod r)` with `L ∈ {0, 1, ..., r-2}`. Then

```text
σ = s^d ≡ 7^{d L} (mod r),     A_d = [σ] P.
```

The map `L ↦ d L (mod r-1)` has kernel of size `gcd(d, r-1) = d` and image
of size `m`. It determines `L` modulo `m` and leaves `L` free modulo `d`.
Concretely there is a unique pair

```text
L = τ + m k,     τ ∈ {0, ..., m-1},     k ∈ {0, ..., d-1}.
```

Phase A recovers `τ`. Phase B recovers `k`. The scalar returned is
`s = 7^{τ + m k} mod r`, accepted only after the point check `[s] P = T`.

## Phase A — digits of `τ` along the prime powers of `m`

Odd primes use the canonical root of unity built from the same primitive
root, so there is no hidden multiplier between two generators of the
order-`q` subgroup.

For an odd prime `q` dividing `m`, set

```text
ζ_q = 7^{(r-1)/q}  mod r.
```

The order of `ζ_q` is exactly `q`: the exponent `(r-1)/q` produces an element
of the unique subgroup of order `q`, and `ζ_q ≠ 1` because `7` is a primitive
root. One field exponentiation builds `ζ_q`. One further exponentiation by
`q` (charged below as a full-size exponentiation) checks `ζ_q^q = 1`.

### First digit

```text
σ^{m/q} = 7^{d L · m / q} = 7^{L (r-1)/q} = ζ_q^{L mod q}.
```

The point form uses only the known integer `m/q`:

```text
Q = [(m/q)] A_d = [ζ_q^{u0}] P,     u0 = L mod q ∈ {0, ..., q-1}.
```

That is one G1 scalar multiplication. Baby-step giant-step then recovers `u0`.

Let `n0 = ceil(sqrt(q))` and `n1 = ceil(q / n0)`, so `n0 · n1 ≥ q`. Every
legal digit is `u0 = i + n0 j` with `i ∈ [0, n0)`, `j ∈ [0, n1)`, and
`i + n0 j < q`. Exponents `≥ q` are computed only as slack in the rectangle
and are ignored. The cover of `{0, ..., q-1}` is exact.

```text
B_0 = P
B_{i+1} = [ζ_q] B_i          i = 0 .. n0-2
S_0 = Q
S_{j+1} = [ζ_q^{-n0}] S_j    j = 0 .. n1-2
```

A hit `S_j = B_i` means `ζ_q^{u0} = ζ_q^{i + n0 j}`, hence `u0 = i + n0 j`
because the order of `ζ_q` is `q` and `i + n0 j` is read only when it lies in
range. The chain costs `(n0 - 1) + (n1 - 1)` scalar multiplications by the
fixed full-size scalars `ζ_q` and `ζ_q^{-n0}`. The stride scalar is one field
exponentiation, `ζ_q^{r-1-n0}`.

### Second digit, only for `q^2`

The squared primes are `906349` and `254760293`. Write
`L ≡ u0 + u1 q (mod q^2)` with `0 ≤ u1 < q`, and set

```text
ψ = 7^{(r-1)/q^2} mod r.
```

Then `ψ` has order `q^2` and `ψ^q = ζ_q`. The same expansion gives

```text
σ^{m/q^2} = ψ^{u0 + u1 q} = ψ^{u0} · ζ_q^{u1},
ζ_q^{u1} = σ^{m/q^2} · ψ^{-u0}.
```

`ψ^{-u0}` is known once `u0` is known. The point

```text
Q1 = [ψ^{-u0} · (m/q^2)] A_d
```

is two scalar multiplications (integer `m/q^2`, then the field scalar
`ψ^{-u0}`). Baby-step giant-step on `ζ_q`, same `(n0, n1)` as the first
digit, returns `u1`. Chinese remainder is not required inside one prime
power: the digit expansion is the representative.

### The factor `2^12` in one search

`2^12` divides `m` and `2^13` does not. Set

```text
ζ = 7^{(r-1)/4096} mod r.
```

This has order `4096`. With `m_odd = m / 4096`,

```text
σ^{m_odd} = 7^{L (r-1)/4096} = ζ^{L mod 4096}.
Q = [m_odd] A_d = [ζ^{L mod 4096}] P.
```

One scalar multiplication builds `Q`. One baby-step giant-step rectangle of
order `4096` (`n0 = n1 = 64`) returns `L mod 4096`. Twelve successive
order-2 lifts are not used.

### Chinese remainder

The odd digits and the `2^12` digit are coprime moduli whose product is `m`.
The ordinary Chinese remainder theorem produces the unique `τ ∈ [0, m)`.
Schoolbook arithmetic on twelve residues of at most 28 bits is charged in the
miscellaneous term below, not omitted.

## Phase B — the remaining digit `k` modulo `2^20`

```text
δ = 7^m mod r
```

has order `(r-1)/gcd(m, r-1) = d = 2^20`. One field exponentiation. Set
`P_τ = [7^τ] P` by one scalar multiplication. Then

```text
s = 7^{τ + m k},     T = [s] P = [δ^k] P_τ,
```

and `k ∈ {0, ..., 2^20 - 1}` is unique because `u ↦ [u] P` is bijective on
`F_r` (`P` has order `r`).

`2^20` is a square. `n0 = n1 = 1024`.

```text
B_0 = P_τ
B_{i+1} = [δ] B_i
S_0 = T
S_{j+1} = [δ^{-1024}] S_j
```

The hit `i + 1024 j` is `k`. This is an exhaustive partition of a set of size
`2^20`, not a sample of a `μ_d` orbit. The final scalar

```text
s = 7^{τ + m k} mod r
```

is one field exponentiation, and the attack accepts it only if `[s] P = T`.
That check is the lift. It is a scalar multiplication, and it is counted.

## What a step costs

Every point step above is a scalar multiplication by a known element of
`F_r`. The scalar has bit length at most 255. Fixed-base windows do not
apply: the point changes at every step. The method charged is left-to-right
double-and-add in Jacobian coordinates on the short Weierstrass model
`y^2 = x^3 + 4` (`a = 0`), after the leading 1. Each of the other at most
254 bits contributes one doubling, and each 1-bit contributes one mixed
addition of a precomputed affine copy of the (changing) base. Worst case:
254 doublings and 254 mixed additions. No Hamming-weight assumption.

Jacobian doubling for `a = 0` is the EFD shape `dbl-2009-l`: `2M + 5S`.
Charging a square as a multiplication is an upper bound, so at most 7 field
multiplications. At 383 instructions each, that is 2681 instructions. Fifteen
field additions at 80 instructions (a 6-limb add-with-carry plus a conditional
subtract of `p` is well under 80) add 1200. Two hundred memory operations and
fifty control operations bring a doubling to at most 4131 instructions,
which is `4131/2000 < 2.07` group-equivalents. The ledger charges **4**
group-equivalents per doubling. The slack is larger than the loop overhead of
a bit test and a pointer update.

A mixed Jacobian-plus-affine addition is charged as **16 field
multiplications**, above the `8M+5S` mixed formulas for `a = 0`. Sixteen
multiplications are 6128 instructions. Twenty field additions at 80 and 400
memory or control operations give at most 8128 instructions, which is
`8128/2000 < 4.07` group-equivalents. The ledger charges **6** per mixed
addition.

One scalar multiplication is therefore at most

```text
254 · 4 + 254 · 6 = 2540
```

group-equivalents. The ledger charges **2560**, leaving 20 group-equivalents
of spare for the bit loop. Call this constant `C_× = 2560`.

Affine normalization of a baby table uses Montgomery batch inversion. For `N`
Jacobian points, `1/Z`, `Z^{-2}` and `Z^{-3}` together with `x = X Z^{-2}` and
`y = Y Z^{-3}` are at most 8 field multiplications per point, plus one field
inversion per table. The inversion is one of the field exponentiations in the
budget below. Giant steps stay in Jacobian and are tested against the affine
table by `X = x Z^2`, `Y = y Z^3`: 4 field multiplications, then a limb
compare. Hashing and the limb compare are charged at 256 instructions per
baby insert and per giant probe, which covers a 12-limb load, a 64-bit mix,
and a compare. That is several times a straight-line implementation.

## Rectangle sizes

`n0 = ceil(sqrt(N))` and `n1 = ceil(N/n0)` for the order `N` of that search.
Chain length `(n0-1)+(n1-1) = n0+n1-2`.

| search | n0 | n1 | chain |
| --- | ---: | ---: | ---: |
| `2^12` | 64 | 64 | 126 |
| `3` | 2 | 2 | 2 |
| `11` | 4 | 3 | 5 |
| `19` | 5 | 4 | 7 |
| `10177` | 101 | 101 | 200 |
| `125527` | 355 | 354 | 707 |
| `859267` | 927 | 927 | 1852 |
| `906349`, digit 0 | 953 | 952 | 1903 |
| `906349`, digit 1 | 953 | 952 | 1903 |
| `2508409` | 1584 | 1584 | 3166 |
| `2529403` | 1591 | 1590 | 3179 |
| `52437899` | 7242 | 7241 | 14481 |
| `254760293`, digit 0 | 15962 | 15961 | 31921 |
| `254760293`, digit 1 | 15962 | 15961 | 31921 |
| phase B, `2^20` | 1024 | 1024 | 2046 |

The fifteen chain lengths sum to **93419**. The `n0` column sums to **46729**.
The `n1` column sums to **46720**. Inserts plus probes: **93449**.

Each of the 15 searches is also charged **2** setup scalar multiplications
(building `Q` or `Q1`, or `P_τ`, and one spare). The final check `[s] P = T`
is one more. Setup total **31**. Scalar multiplications:

```text
S_× = 93419 + 31 = 93450
```

## Field exponentiations

A 255-bit exponentiation is at most 254 squares and 253 multiplications in
`F_r`: 507 field multiplications, `507 · 383 = 194181` instructions,
`194181/2000 < 97.1` group-equivalents. The ledger charges **120** per field
exponentiation and budgets **80** of them.

The procedure uses at most 74:

- 12 roots `7^{(r-1)/q}` and 12 order checks
- 2 roots `7^{(r-1)/q^2}` for the squared primes
- 1 root `7^{(r-1)/4096}` and 1 order check
- 12 primitive-root checks of `7`
- 15 stride scalars
- `7^m`, the final `7^{τ+mk}`, and 15 batch inversions

Eighty is a hard cap above that list. Unused slots are slack, not extra
searches.

## Ledger

```text
scalar multiplications     93450 · 2560           = 239232000
field multiplications      8·46729 + 4·46720
                           = 560712
                           · 383 = 214752696 ins
                           ceil(/2000)            =    107377
hashes and limb compares   93449 · 256 = 23922944 ins
                           ceil(/2000)            =     11962
field exponentiations      80 · 120               =      9600
integer arithmetic, CRT,
output encoding            2000 · C instructions  =      2000
                                                    -----------
W(n)                                               = 239362939
```

The miscellaneous 2000 group-equivalents are `4 · 10^6` instructions, above
any schoolbook Chinese-remainder or 256-bit division used here.

`2^27 = 134217728` and `2^28 = 268435456`, so `27 < log2(W) < 28`. The
rational `27834625/1000000` is the least multiple of `10^{-6}` that is `≥
log2(239362939)`, i.e. the least numerator for which `2^{27834625} ≥
239362939^{1000000}`. Therefore

```text
log2(W(2^20)) ≤ 27.834625.
```

## Memory, preprocessing, advice

Tables are built one prime power at a time and discarded. The largest baby
table has 15962 affine points. At 96 bytes each that is 1,532,352 bytes.
A half-full hash table of 32,768 slots at 24 bytes is 786,432 bytes. Scratch,
code and stack are under `2^20` bytes. Peak is under `2^23` bytes.

Preprocessing, in the sense of every field element and every baby table, is
at most half the scalar multiplications plus the whole field budget. That is
under `1.21 · 10^8` group-equivalents, which is under `2^27`, and it is
inside `W` rather than added again.

No table of searched discrete logs is embedded. The primitive root `7` and
the factorization are recomputed or rechecked inside the charged
exponentiations. Advice is 0 bytes.

## Success

Every branch is a finite exact cover of a residue ring whose order is known
from the factorization. The primitive-root and order checks are equalities in
`F_r` that hold for these integers, so the algorithm does not restart and
does not fall back to a generic search. The probability over uniform `s` and
over algorithm coins is 1, which is above the track threshold 0.39. Coins are
not consumed.

## Range

For `2^20 ≤ n < r` the same point `A_{2^20}` is an entry of the supplied
prefix, and `W(n)` does not grow with `n`. The attack does not read the other
powers. For `n < 2^20` this bound is not claimed. The scored point
`n = 2^20` lies in the claimed range.

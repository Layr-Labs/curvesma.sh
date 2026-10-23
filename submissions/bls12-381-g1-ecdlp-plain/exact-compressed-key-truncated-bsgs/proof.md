# bls12-381-g1-ecdlp-plain: exact compressed-key truncated baby-step/giant-step

Deterministic, hash-free, and exactly accounted. Both tables are keyed on the
48-byte compressed affine encoding of a G1 point, which is injective. The
recoverable set is therefore exactly the integer interval `[0, m^2)`, the success
probability is exactly `m^2/r`, at most one candidate ever reaches verification,
and every cost in the protected `classical-group-work-v2` model is charged
explicitly, including multiword counters, multiword scatter offsets, record loads
as well as stores, and the batching construction behind the 6M affine-addition
price.

## Target and parameters

The protected profile fixes the BLS12-381 G1 group of prime order

```
r = 52435875175126190479447740508185965837690552500527637822603658699938581184513
```

with generator `P0`, and supplies only `T = [x]P0` for uniform `x` in `[0, r)`,
including `x = 0` (`T` the identity). There are no auxiliary powers and no `n`
parameter. The attack must return the scalar `x`, and the acceptance test is
`[x]P0 = T`.

```
m = 143003466105892742053438571826631282054     (127 bits)
h = m - 1 = 143003466105892742053438571826631282053
m^2 = 20449991318299214286984618798192526676918749446216934733046686138239702458916
```

`m` is the least integer with `m^2 >= ceil(39r/100)`; consequently `m^2 < r`
(the ratio is `0.39 + 2^-128`-scale), and `m*h = m^2 - m < r`. Both facts are
used below. The reference price for this track is `C = 2000` scalar AArch64
instructions per target-group operation, one operation being a batched
incomplete affine G1 addition under the 6M convention.

## Algorithm

Compression. `comp(P)` is the 48-byte canonical compressed encoding: the
big-endian `x`-coordinate with the three flag bits of the first byte, where the
sign bit selects which of the two square roots is `y`, and the identity has the
reserved all-flag encoding. It is the encoding of the pinned `py_ecc`
serializer.

```text
S = m*P0                                   # one 127-bit scalar multiplication
Baby  : for j = 0..m-1, B_j = j*P0         # affine chain, batched inversion
Giant : for i = 0..m-1, G_i = T - i*S      # affine chain, batched inversion
Store 2m records (key, index) with key = comp(point), index = (tag, i or j),
tag in {baby, giant} and the index 128 bits wide.
LSD radix sort all 2m records on the 48-byte key: 12 passes of a 32-bit digit,
2^32 buckets, 128-bit counters and 128-bit scatter offsets.
Scan adjacent records with equal keys.  Every emitted pair (i, j) has different
tags; at most one such pair exists.
Candidate x = i*m + j.  Verify [x]P0 = T and output x on success, else fail.
```

Nothing is randomised: the tables, their order, the scan and the returned scalar
are a deterministic function of `(P0, T)`. There are no algorithm coins, no
restarts and no amplification.

## Exactness of the recoverable set

Lemma 1 (compression is injective on G1). For a point `(x, y)` of the curve
`y^2 = x^3 + 4` over `Fp`, the compressed encoding stores `x` and one sign bit.
The two candidates for `y` are `+sqrt(x^3+4)` and its negative; the sign bit
selects one of them. `y = 0` would give a point of order 2, which cannot lie in
the odd-order group G1, so the two roots are distinct and the sign bit is
unambiguous. The identity has its own reserved encoding, disjoint from every
affine encoding. Hence `comp(P) = comp(Q)` iff `P = Q` for `P, Q` in G1. In
particular `comp(P0) != comp(-P0)` whenever `P0` is not the identity: the two
points share an `x`-coordinate but differ in the sign bit.

Lemma 2 (keys are distinct inside each table). `j*P0 = j'*P0` with
`0 <= j, j' <= m-1 < r` gives `j = j'`; likewise `T - i*S = T - i'*S` with
`0 <= i, i' < m` gives `(i-i')*m = 0 mod r`, hence `i = i'` because `m < r` and
`S = m*P0` has order `r`. Each table therefore holds `m` distinct keys.

Lemma 3 (a key match is a scalar representation). A cross-table match at
`(i, j)` means `comp(T - i*S) = comp(j*P0)`, which by Lemma 1 is
`T - i*m*P0 = j*P0`, i.e. `T = (i*m + j)*P0`, i.e. `x = i*m + j` modulo `r`.

Lemma 4 (the rectangle is injective and its image is the interval). Let
`0 <= i, i' < m` and `0 <= j, j' <= m-1`. If `i*m + j = i'*m + j'` modulo `r`
then `m` divides `j' - j`; but `|j' - j| <= m-1 < m`, so `j = j'` and `i = i'`.
The map is therefore injective on the rectangle `[0, m) x [0, m-1]`, whose
image is exactly the integer interval `[0, m^2)` because every `y` in
`[0, m^2)` is `y = (y div m)*m + (y mod m)` with `y div m < m` and
`y mod m <= m-1`. Injectivity uses `m*h = m^2 - m < r`; the image uses only
`m^2 <= r`, which holds.

Theorem (exact success set and probability). For `T = [x]P0` with
`x` in `[0, r)`:

* if `x < m^2`, write `x = i*m + j` with `0 <= i < m`, `0 <= j <= m-1`; then
  `T - i*S = j*P0` exactly, the two keys are equal, and the scan emits `(i, j)`;
* if `x >= m^2`, suppose a match `(i, j)` existed. Then `x = i*m + j` with
  `i*m + j <= m^2 - 1`, contradicting `x >= m^2`. So no match is emitted.

Hence the recoverable set is exactly `[0, m^2)` and, for uniform `x`,
`Pr[success] = m^2/r`, which is at least `39/100` and in fact
`0.39 + (m^2 - 39r/100)/r`.

Candidate count and verification. By Lemma 4 the rectangle carries `m^2`
distinct scalars, so two distinct matching pairs would give two distinct scalars
in `[0, m^2)` congruent modulo `r < m^2 + m`, which is impossible. The scan
therefore emits at most one pair, and the attack performs at most one
verification scalar multiplication. Every emitted candidate is checked with
`[x]P0 = T`, so a returned scalar is always correct; failure returns nothing.

The identity target is included: `x = 0` is `(i, j) = (0, 0)`, `T` is the
identity, `T - 0*S` is the identity and `0*P0` is the identity, so the reserved
encodings match and the attack returns `0`.

Remark on the previously asserted failure mode. A point with `x = r - 1` is
`T = -P0`. Its compressed encoding has the sign bit opposite to `comp(P0)`, so
the giant record `G_0 = T` does not match the baby record `B_1 = P0`. Since
`r - 1 >= m^2`, the attack correctly returns nothing. There is no minus
candidate and no additional recovered interval: the exact set is `[0, m^2)`
alone. Keying on the affine `x`-coordinate while dropping the sign bit would
have produced exactly the extra `[r-m+1, r)` hits that this encoding removes.

## Batch-inversion schedule behind the 6M price

Both tables are dependent affine chains: the baby chain adds `P0` repeatedly and
the giant chain adds `-S` repeatedly. A standalone affine addition needs a field
inversion, so a batching construction is required before charging the reference
rate.

Block size `B = 2^20`. For a block of `B` chained additions, accumulate the
denominators `d_1..d_B`, compute the prefix products `D_k = d_1...d_k`
(`B-1` field multiplications), invert `D_B` once, and recover every `d_k^-1` by
one multiplication each plus one running multiplication (`2B` multiplications).
The batch cost is therefore `3B` multiplications plus one inversion, i.e.
`3M` per addition for the batch, on top of the `3M` of the incomplete-affine
addition formula itself: exactly the `6M` of the reference price.

Inversion cost is not hidden: one inversion by binary powering `p-2` uses at most
571 base-field multiplications, each measured at 383 AArch64 instructions in the
pinned blst build, i.e. 218,693 instructions. Amortised over `B = 2^20`
additions that is 0.209 instructions per addition, far below the 1,149
instructions that the 3M batch allowance already grants. The reference price is
thus an over-charge, and no separate inversion term is added to avoid
double-counting the reference operation and its implementation.

The two scalar multiplications (`S = m*P0`, one verification `[x]P0`) are
charged as ordinary reference operations, at 512 additions each, which is above
the binary double-and-add bound of `2*bitlength` for 127-bit and 255-bit
exponents.

## Work ledger

Units. `H` counts target-group operations (one batched incomplete affine G1
addition each). `R` counts other scalar AArch64 instructions. Total work is
`W = H + R/C` with `C = 2000`, as the protected model requires. Records are
`48`-byte keys plus `128`-bit indices; the index width is forced because
`N = 2m` has 128 bits.

Sort. `N = 2m` records, LSD radix on the 48-byte key, 12 passes of a 32-bit
digit, `2^32` buckets with 128-bit counters (16 bytes each, since a bucket count
can reach `N`), and 128-bit scatter offsets. Each pass is one counting traversal
plus one scatter traversal.

| Counting traversal, per record | Instructions |
| --- | ---: |
| form key address: 3-word base + 128-bit index x 48 (two 64x64 multiplies, adds, carry) | 12 |
| load 32-bit digit, isolate by shift and mask | 5 |
| load 128-bit bucket counter (2 loads) | 2 |
| increment counter with carry (add, adc) | 2 |
| store 128-bit counter (2 stores) | 2 |
| loop pointer update and branch | 2 |
| **count total** | **25** |

| Scatter traversal, per record | Instructions |
| --- | ---: |
| form key address and load the digit (as above) | 17 |
| load 128-bit bucket counter (2 loads) | 2 |
| form destination address: scratch base + counter x 16 in 128-bit arithmetic, including the three-word carry chain | 16 |
| increment counter with carry and store it (4) | 4 |
| load source 128-bit index (2 loads) | 2 |
| store destination 128-bit index (2 stores) | 2 |
| loop pointer update and branch | 2 |
| **scatter total** | **45** |

One pass costs `25 + 45 = 70`; twelve passes cost `840` per record. Record
materialisation (write the 48-byte key as six stores plus addressing, and the
128-bit index as two stores) costs `32`. The final adjacent-key scan visits each
of the `2m` records once: two 6-word key gathers (12 loads), up to six
comparisons and a tag test, plus pointer updates, cost `48` per record. Total per
record `840 + 32 + 48 = 920`. Both loads and stores of keys, counters, offsets
and indices are counted; no multiword operation is treated as single-word.

Bucket maintenance outside the two traversals is charged separately: clearing
and prefix-summing `2^32` 128-bit counters takes at most 6 instructions each, so
`12 * 2^32 * 6` instructions over the twelve passes. A fixed 4,096 instructions
cover output construction, the tag test and program state.

```
R = 920*N + 12*2^32*6 + 4096
  = 263126377634842645378326972161310796628768
H = 2*(m-1) + 1024
  = 286006932211785484106877143653262565130
W = H + R/2000
  = 52196265128650850849505078716739745430548/125
log2(W) = 128.29528609017265
```

The declared bound `2^128.295287` is the exact `log2(W)` rounded upward at six
decimal places. Memory and preprocessing below are rounded upward the same way.
No component of `W` is an expectation, an amortisation over future targets, or a
free oracle.

## Memory

| Component | Bytes |
| --- | ---: |
| key arrays, `2m` x 48 | 96m |
| index arrays, `2m` x 16 | 32m |
| sort scratch indices, `2m` x 16 | 32m |
| bucket counters, `2^32` x 16 | 2^36 |
| code, constants, stacks, fixed scratch | 2^20 |

Peak storage is `160m + 2^36 + 2^20 = 22880554576942838728550171492329725653952`
bytes, i.e. at most `2^134.071246` bytes. The counter array and the three record
arrays are never all live in a larger form than this; the sort is an index
permutation, so no second copy of the 48-byte keys is needed. Memory is
astronomical and mandatory for this construction; it is not a practical attack.

## Preprocessing and advice

Preprocessing is the whole baby table: `m-1` chained additions, materialisation
of its `m` records (`32m` instructions) and its share of the sort (`840m`
instructions plus half the bucket maintenance). In the model's units

```
preprocessing = (m-1) + (32m + 840m + 12*2^32*6/2)/2000
              = 25669122166007747198592223642889978804984/125
              = 2.0535297732806198e38,  log2 = 127.27137347022833
```

so at most `2^127.271374` target-group-operation equivalents, and it is included
in `W` in full with no reuse credit. Nonuniform advice is zero bytes: the baby
table depends only on the public generator but is rebuilt and fully charged on
every attack, so nothing is stored or carried between targets.

## Determinism, heuristics and limitations

There are no heuristic premises: no random-oracle model, no hash-load-factor
assumption, no distributional claim about probe sequences, no cycle-finding
heuristic, and no smoothness or linear-algebra heuristic. The sort is a
deterministic worst-case radix sort; the scan is exact key equality; the success
statement is a proven identity over the uniform challenge distribution. The only
imported constants are the protected profile (`r`, `P0`, the curve), the pinned
reference price `C = 2000`, the 6M affine convention and the measured 383
instructions per base-field multiplication. The accounting is analytic; no
full-size execution is claimed or needed. `experiment.py` checks the coverage
identity on small additive groups only and is not evidence about production
encodings or instruction counts.

## Submitted bound

| Quantity | Value |
| --- | ---: |
| total work `W` | `52196265128650850849505078716739745430548/125` |
| `time_log2` | **128.295287** |
| peak memory | `22880554576942838728550171492329725653952` bytes, `2^134.071246` |
| preprocessing | `25669122166007747198592223642889978804984/125`, `2^127.271374` |
| nonuniform advice | 0 bytes |
| success probability | exactly `m^2/r >= 39/100` |

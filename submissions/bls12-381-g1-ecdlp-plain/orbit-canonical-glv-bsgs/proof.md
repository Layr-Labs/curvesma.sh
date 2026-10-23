# bls12-381-g1-ecdlp-plain: orbit-canonical GLV baby-step/giant-step

Deterministic and hash-free. Baby points `a*P0` and giant points `T - b*[λ]P0`
are each stored once. The stored key is the minimum 48-byte compressed
encoding of the point's order-6 orbit under the GLV endomorphism and negation.
A cross-list hit is orbit equality, so the scalar is `b*λ + σ*a` for an
explicit unit `σ`. The recoverable set contains an explicit set `S` of
cardinality `(m-1)*(3*m-1)`. That cardinality is at least `ceil(39r/100)`, so
the success probability is at least `39/100`. One candidate is verified.

Every cost in the protected `classical-group-work-v2` model is charged,
including the two base-field multiplications that produce the conjugate
x-coordinates, the orbit-minimum selection, the radix sort, and the batching
construction behind the 6M affine-addition price.

## Target and parameters

The protected profile fixes the BLS12-381 G1 group of prime order

```
r = 52435875175126190479447740508185965837690552500527637822603658699938581184513
```

with generator `P0`, and supplies only `T = [x]P0` for uniform `x` in `[0, r)`,
including `x = 0`. There are no auxiliary powers. The attack returns the scalar
`x`, and acceptance is `[x]P0 = T`.

The GLV eigenvalue is the integer

```
λ = 228988810152649578064853576960394133503
```

which satisfies the integer identity `λ² + λ + 1 = r`. In particular `λ² + λ + 1 ≡ 0 (mod r)`, so `λ` is a primitive cube root of unity modulo `r`, and `-λ` has order 6. The six units are

```
{ ±1, ±λ, ±λ² } = { ±1, ±λ, ±(1+λ) }
```

as residues, because `λ² ≡ -1-λ` and `-λ² ≡ 1+λ`.

```
m = 82563089651286698697488428308387906491
h = m - 1 = 82563089651286698697488428308387906490
```

Both are 126-bit integers, strictly less than `2^126`. The constraints used by
the cardinality proof hold:

```
2*m - 1 = 165126179302573397394976856616775812981 ≤ λ
(m - 1) + (h - 1) = 2*m - 3 ≤ λ - 1
h = m - 1
```

and `|S| = (m-1)*(3*m-1)` is the least such boundary value that is at least
`ceil(39r/100)`. The predecessor `m-1` falls short of that threshold. The
reference price is `C = 2000` scalar AArch64 instructions per batched
incomplete affine G1 addition under the 6M convention.

## Endomorphism

Let `β` be the primitive cube root of unity in `Fp`

```
β  = 0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac
β² = 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe
```

so `β² + β + 1 ≡ 0 (mod p)` and `β²` is the other primitive cube root. These are
public curve constants, not searched advice. Define `φ(x, y) = (βx, y)` and
`φ(O) = O`. If `(x, y)` lies on `y² = x³ + 4`, then so does `(βx, y)`, because
`(βx)³ + 4 = β³ x³ + 4 = x³ + 4`. The same holds for `β²`.

`φ` is a group homomorphism. On the chord-and-tangent law the slope scales by
`β⁻¹` when both x-coordinates scale by `β`, and the resulting x-coordinate
matches `β` times the unscaled sum precisely because `β³ = 1`. The y-coordinate
matches by the same cancellation. Thus `φ³ = id` and `φ` commutes with negation.

On the cyclic subgroup `G1 = <P0>` every homomorphism is multiplication by a
scalar. Since `φ³ = id` and `φ ≠ id`, that scalar is a primitive cube root of
unity modulo `r`, hence `λ` or `λ² = -1-λ`. The finite check

```
φ(P0) = (β * x(P0), y(P0)) = [λ]P0
```

selects `λ`. `experiment.py` evaluates this identity on the profile generator
by a straight double-and-add. Both `φ` and `[λ]` are homomorphisms, they agree
on `P0`, and therefore `φ(P) = [λ]P` for every `P` in `G1`. Consequently
`φ²(P) = [λ²]P = [-(1+λ)]P` and `-φ²(P) = [(1+λ)]P`.

The giant stride point is `φ(P0)`, obtained while the baby point `1*P0` is
canonicalised. It is not a 128-bit scalar multiplication.

## Canonical compressed key

`comp(P)` is the 48-byte compressed encoding of the pinned serializer. A finite
point stores `x` in bits `0..380` and three flags in bits `381..383`: bit 383
marks compression, bit 382 marks infinity, and bit 381 (mask `0x20`) is set
iff `y` is the lexicographically larger square root of `x³+4`. For `y` in
`[0, p)` that comparison is `2y > p`. The two roots are distinct on `G1`
because `r` is an odd prime, so exactly one of `y` and `p-y` has bit 381
clear. Negation preserves `x` and flips bit 381. The identity uses the reserved
infinity encoding.

Bit 381 is more significant than every bit of `x`. Among the six finite
encodings of a non-identity orbit, every encoding with bit 381 clear is
strictly smaller than every encoding with bit 381 set. The three points of the
orbit with bit 381 clear have x-coordinates `{x, βx, β²x}`: `φ` and `φ²`
preserve `y` and therefore preserve the sign bit, and negation does not change
`x`. The canonical encoding of the orbit is therefore compression set, infinity
clear, sign clear, and

```
x_can = min(x, βx mod p, β²x mod p).
```

An x-coordinate determines at most the pair `{Q, -Q}`, which is a single orbit.
Distinct orbits therefore have distinct canonical encodings. The identity orbit
is only the infinity encoding, distinct from every finite encoding. Computing
`x_can` takes two reduced base-field multiplications by the public constants
`β` and `β²`. No y-coordinate is required to build the key.

## Algorithm

```text
Baby chain: B_a = a*P0 for a = 0..m-1          # m-1 batched affine additions
Stride:     Q = φ(P0) = [λ]P0                  # the β-multiple already formed for B_1
Giant chain: G_b = T - b*Q for b = 0..h-1      # h-1 batched affine additions

Baby list, m records:
  (canonical comp of the orbit of B_a, index a)
Giant list, h records:
  (canonical comp of the orbit of G_b, index b)

LSD radix-sort each list on the 48-byte key: 12 passes of a 32-bit digit,
2^32 buckets, 128-bit counters, 128-bit indices. The sorts are sequential and
reuse one counter array.
Merge by canonical key. A hit means the orbits are equal, so
  x - b*λ ≡ σ * a (mod r)
for a unique unit σ in {±1, ±λ, ±(1+λ)}. Test the six units by comparing
x-coordinates (the products βx and β²x are the ones already formed) and emit
  candidate = (b*λ + σ*a) mod r.
Verify [candidate]P0 = T. Output the scalar on success, else fail.
```

Nothing is randomised. There are no algorithm coins, no restarts and no
amplification. Baby indices and giant indices live in separate lists and both
fit in 128 bits. The record count `N = m + h = 2m - 1 < 2^128`, so every
bucket counter fits in 128 bits.

## The explicit set

Let

```
Σ0 = {1, -1, 1+λ, -(1+λ)}
S  = { (b*λ + σ*a) mod r : 0 ≤ b < h, 0 ≤ a < m, σ ∈ Σ0, and σ = 1 whenever a = 0 }.
```

`S` is the set whose cardinality is proved below. The six-unit orbit match
also recovers residues that use `σ = ±λ`. Those extra residues, if any, only
enlarge the recoverable set. The success lower bound uses `|S|` alone.

Every residue in `[0, r-2]` has a unique writing `A + B*λ` with `0 ≤ A < λ`
and `0 ≤ B ≤ λ`. The count of such pairs is `λ(λ+1) = λ² + λ = r - 1`, and
the values fill `[0, λ² + λ - 1] = [0, r-2]` without repetition: a collision
`A + Bλ = A' + B'λ` with a nonzero `B` gap of size at least 1 forces an `A`
gap of size at least `λ`, which cannot fit inside `[0, λ)`. The missing residue
`r - 1` is handled as the negation class below and is inside `S`.

Write `P+` for pairs `(a, b)` giving `a + bλ` with `0 ≤ a < m`, `0 ≤ b < h`.
The canonical coordinates are `(A, B) = (a, b)`. The maximum value is at most
`(m-1) + (h-1)λ ≤ (m-1) + (λ - m)λ < r`, using `(m-1)+(h-1) ≤ λ-1`, so there
is no modular wrap. This contributes `m*h` distinct residues.

Write `P-` for ` -a + bλ` with `1 ≤ a < m`.

- If `b ≥ 1`, then `-a + bλ = (λ - a) + (b-1)λ`, so `(A, B) = (λ-a, b-1)` with
  `A ∈ [λ-(m-1), λ-1]` and `B ∈ [0, h-2]`. The constraint `2m-1 ≤ λ` is
  `λ - (m-1) ≥ m`, so this `A` is at least `m` and does not meet `P+`, whose
  `A` is at most `m-1`. `B ≤ h-2 ≤ λ-2`, inside the canonical range.
- If `b = 0` and `a = 1`, the residue is `r-1`, outside the canonical square
  and distinct from every residue in `[0, r-2]`.
- If `b = 0` and `a ≥ 2`, then `r - a = (λ+1-a) + λ*λ`, so `(A, B) = (λ+1-a, λ)`.
  Here `B = λ` does not meet `P+` (`B ≤ h-1 ≤ λ-1`) or the `b ≥ 1` rows
  (`B ≤ h-2`).

So the signed family `P+ ∪ P-` has `h*(2m-1)` distinct residues: `a = 0` once
per giant index, and both signs for each `a ∈ [1, m)`.

Write `Q+` for `(1+λ)a + bλ = a + (a+b)λ` with `1 ≤ a < m` and `0 ≤ b < h`.
Canonical coordinates `(A, B) = (a, a+b)`. The constraint
`(m-1)+(h-1) ≤ λ-1` keeps `B ≤ λ-1`. This point lies in `P+` exactly when
`a+b < h`. Because `h ≥ m-1 ≥ a`, the giants `b ∈ [h-a, h-1]` are legal and
there are exactly `a` of them. Summing `a` from 1 to `m-1` gives `m(m-1)/2`
residues with `B ≥ h`, disjoint from `P+`. Their `A = a ≤ m-1` is disjoint
from `P-`, whose positive-`b` rows have `A ≥ m` and whose `b = 0` rows have
`B = λ`.

Write `Q-` for `-(1+λ)a + bλ` with `1 ≤ a < m` and `0 ≤ b < a` (the condition
`b ≥ a` reproduces a `P-` residue and is not counted again). There are
`m(m-1)/2` such pairs because `h = m-1 ≥ a-1`, so every `b ∈ [0, a)` is a
legal giant index. For `a = 1`, necessarily `b = 0`, and the canonical
coordinates are `(A, B) = (0, λ)`. For `a ≥ 2`,

```
(λ - a + 1) + (λ - a + b)λ = r + (-(1+λ)a + bλ),
```

and `0 < λ-a+1 < λ`, `0 ≤ λ-a+b ≤ λ-1`, so the canonical coordinates are
`(λ-a+1, λ-a+b)`. The `A` coordinate is at least `λ-m+2 ≥ m`, hence disjoint
from `P+` and from `Q+`. It is disjoint from `P-` at `b ≥ 1` because that
family has `B ≤ h-2` while here `B ≥ λ-a ≥ λ-(m-1) ≥ h` under
`(m-1)+(h-1) ≤ λ-1`. It is disjoint from `P-` at `b = 0` because that family
has `B = λ` while here `B ≤ λ-1`, except the already-separated point
`(0, λ)`, whose `A = 0` is not of the form `λ+1-a`.

The four blocks are pairwise disjoint and the two new blocks each contribute
`m(m-1)/2` residues. Therefore

```
|S| = h*(2*m - 1) + m*(m - 1) = (m - 1)*(3*m - 1)
    = 20449991318299214286984618798192526677058846637647322474230865784471728173280.
```

The excess over `ceil(39r/100)` is `359531162441543723415438891495681511319`.
In particular `|S|/r ≥ 39/100`.

Baby scalars themselves are distinct, so baby orbits are distinct. The four
baby residues `±a` and `±(1+λ)a` for fixed `a ∈ [1, m)`, together with `0`,
are the `b = 0` slice of the same blocks (with `Q-` at `b = 0 < a`). The
coordinate separation above does not use `h` except to keep giant indices
legal, and the same `A` gaps apply. Hence the `m` baby canonical keys are
pairwise distinct. Giant keys at different `b` may coincide for a particular
`T`; the merge keeps any one copy. A coincidence means two giants lie in one
orbit, and either copy decodes to the same scalar.

## Recovery

If `x ∈ S`, the unique block representation supplies `b`, `a` and `σ ∈ Σ0`
with `x ≡ bλ + σ a (mod r)`. Then `G_b = [σ a]P0` lies in the orbit of `B_a`,
the canonical keys match, and the six-unit test accepts that `σ`. The final
check `[x]P0 = T` succeeds.

If a hit occurs for any unit, including `±λ`, the decoded candidate satisfies
the same linear relation and the verification either accepts the true scalar or
rejects a mismatch. Distinct blocks do not produce two different scalars for
one `T`: each hit determines `x` from `T = [x]P0`. At most the six units on a
matched orbit are tested, and the λ-lines duplicate the `a = 0` family rather
than introducing a second scalar. The ledger charges four full verifications
(`2048` reference additions) so that every accepted decoding is paid for even
if the implementation verifies each unit separately. A rejected verification
does not emit a scalar.

Outside `S` the attack may still succeed via a `±λ` representation. It never
emits a scalar that fails `[x]P0 = T`. The success probability on uniform `x`
is therefore at least `|S|/r`.

## Batch-inversion schedule behind the 6M price

Both chains are dependent affine additions of a fixed point (`P0` on the baby
side, `-Q` on the giant side). A standalone affine addition spends a field
inversion. The reference rate is applied only with batching in front of it.

Block size `B = 2^20`. For a block, accumulate denominators, form the prefix
products (`B-1` field multiplications), invert the product once, and recover
each denominator inverse (`2B` multiplications). Together with the `3M` of the
incomplete-affine formula this is the `6M` reference price. One inversion by
`p-2` is at most 571 base-field multiplications of 383 instructions each,
amortised to well under one instruction per addition across a block of `2^20`.
No separate inversion term is added. `S = M` and `A = 0` are the reference
convention; the orbit-key field multiplications below are charged as
instructions, not as reference additions.

## Work ledger

`H` counts target-group operations. `R` counts other scalar AArch64
instructions. `W = H + R/C` with `C = 2000`.

```
N = m + h = 165126179302573397394976856616775812981
```

The identity baby `a = 0` stores the infinity encoding: 32 instructions to
materialise the key and index, 840 to sort, 48 to merge, total 920. Every
other stored record pays the following.

| Orbit-minimum selection, per non-identity point | Instructions |
| --- | ---: |
| reduced product `x*β` (`blst_fp_mul`, pinned measurement) | 383 |
| reduced product `x*β²` | 383 |
| compare `x` with `βx`: 6 high-to-low limbs, select | 32 |
| compare the winner with `β²x` | 32 |
| clear the sign bit and address the winner | 0, folded into the 32 below |
| **selection** | **830** |

Each six-limb integer comparison is charged at 32 instructions, so the two
comparisons sum to 64. Materialisation of the winning
48-byte key is the same 32-instruction charge used for a chain record (six
stores, addressing, two index stores). Sort and merge match the 12-pass
32-bit radix accounting: 70 instructions per pass, 12 passes, 840; merge at
48 per record (six loads of each cursor, six comparisons, pointer, branch).

| Per non-identity stored record | Instructions |
| --- | ---: |
| selection | 830 |
| store canonical key and 128-bit index | 32 |
| 12-pass radix sort | 840 |
| merge | 48 |
| **total `PER_POINT`** | **1750** |

Sort detail, unchanged in structure from a single 32-bit-digit pass. Counting
costs 25 instructions per record (address a 48-byte key from a 128-bit index,
isolate the digit, load/increment/store a 128-bit counter, loop). Scattering
costs 45 (digit, destination address from a 128-bit counter, counter update,
load and store the 128-bit index, loop). Bucket maintenance, paid twice
because the lists are sorted sequentially and the counter array is reused, is
`2 * 12 * 2^32 * 6`: clear and prefix-sum `2^32` counters of 16 bytes at 6
instructions each. A fixed 8192 instructions cover candidate decoding (at most
six units, each a compare of already-formed field elements), output, and
program state.

```
H = (m - 1) + (h - 1) + 2048
  = 165126179302573397394976856616775815027
R = 920 + 1750*(m - 1) + 1750*h + 2*12*2^32*6 + 8192
  = 288970813779503445441209499079976148014736
W = H + R/2000
  = 38701448274040640014447700769595486129296 / 125
log2(W) = 127.8637270663...
```

The declared bound `2^127.863728` is `log2(W)` rounded upward at six decimal
places. Memory and preprocessing below are rounded upward the same way. No
component of `W` is an expectation, an amortisation over future targets, or a
free oracle. The 2048 reference additions are four verification scalar
multiplications, each charged at the binary double-and-add ceiling `2*256`.

## Memory

Sorts run one list at a time and permute indices. Both lists are live for the
merge. Scratch indices exist only for the list being sorted.

| Component | Bytes |
| --- | ---: |
| keys, `N` × 48 | `48*N` |
| indices, `N` × 16 | `16*N` |
| scratch indices, `max(m, h)` × 16 | `16*m` |
| bucket counters, one `2^32` × 16 array, reused | `2^36` |
| code, constants, stacks, fixed scratch | `2^20` |

Peak storage is

```
64*N + 16*m + 2^36 + 2^20
  = 11889084909785284612438333676476579059952
```

bytes, i.e. at most `2^133.126762` bytes. The chains' affine coordinates may be
discarded once the canonical key and the two conjugate x-coordinates needed for
a hit are written; the bound above already covers the stored tables. Memory is
mandatory for this construction; it is not a practical attack.

## Preprocessing and advice

Preprocessing is everything independent of `T`: the baby chain (`m-1`
additions), canonicalisation of the `m` baby orbits, the baby list's 12-pass
sort, and one bucket-maintenance pass. `φ(P0)` is produced while canonicalising
`B_1` and is not a second chain. The giant chain, the giant sort and the merge
read `T`.

```
preprocessing = (m - 1) + (872 + (766+64+32+840)*(m - 1) + 12*2^32*6) / 2000
              = 76412139472265839644525540399490316868041 / 500
log2 = 126.8451413659...
```

so at most `2^126.845142` target-group-operation equivalents. The baby-sort and
baby-selection terms are the baby portion of `R` and are not added again on top
of `W`. The identity baby's 872 instructions are a 32-instruction store plus
the 840-instruction sort. Nonuniform advice is zero bytes: `β`, `β²`, `λ` and
`P0` are public profile constants, the baby table is rebuilt for the attack,
and its construction is fully charged.

## Determinism, heuristics and limitations

There are no heuristic premises: no random-oracle model, no hash-load factor,
no cycle-finding heuristic, and no smoothness or linear-algebra heuristic. The
endomorphism identity on `G1` is the homomorphism property plus one generator
evaluation. Canonical keys are a total order on encodings, not a hash. Both
sorts are deterministic worst-case radix sorts. The success statement is a
lower bound by an explicit set, not an expectation over algorithm coins. The
attack has no coins.

The accounting is analytic. No full-size execution is claimed. `experiment.py`
checks the cardinality identity on groups of order `λ²+λ+1` for small `λ`, and
checks `φ(P0) = [λ]P0` once. The exponential tables are required by this
construction. The attack is not claimed to be practical, novel, or optimal
among generic algorithms.

## Submitted bound

| Quantity | Value |
| --- | ---: |
| total work `W` | `38701448274040640014447700769595486129296/125` |
| `time_log2` | **127.863728** |
| peak memory | `11889084909785284612438333676476579059952` bytes, `2^133.126762` |
| preprocessing | `76412139472265839644525540399490316868041/500`, `2^126.845142` |
| nonuniform advice | 0 bytes |
| success probability | at least `|S|/r ≥ 39/100` |

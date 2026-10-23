# bls12-381-g1-ecdlp-plain: unique GLV-coordinate rectangle baby-step/giant-step

Deterministic and hash-free. Baby points `a*P0` and giant points `T - b*[λ]P0`
are each stored once. The stored key is the minimum 48-byte compressed
encoding of the point's order-6 orbit under the GLV endomorphism and negation.
A cross-list hit is orbit equality, so the scalar is `b*λ + σ*a` for an
explicit unit `σ`. The recoverable set contains an explicit set `S` of
cardinality `h*(4*m - h - 2)`. That cardinality is at least `ceil(39r/100)`,
so the success probability is at least `39/100`. One candidate is verified.

Every cost in the protected `classical-group-work-v2` model is charged,
including the two base-field multiplications that produce the conjugate
x-coordinates, the orbit-minimum selection, the radix sort, the double-and-add
that resolves the unit after a hit, and the batching construction behind the
6M affine-addition price.

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

which satisfies the integer identity `λ² + λ + 1 = r`. In particular
`λ² + λ + 1 ≡ 0 (mod r)`, so `λ` is a primitive cube root of unity modulo `r`,
and `-λ` has order 6. The six units are

```
{ ±1, ±λ, ±λ² } = { ±1, ±λ, ±(1+λ) }
```

as residues, because `λ² ≡ -1-λ` and `-λ² ≡ 1+λ`.

```
m = 95929641369258993122728749523035932376
h = 63953094246172662083740702189702033217
```

Both are 126-bit integers. The rectangle constraints used by the cardinality
proof hold with room:

```
h = 63953094246172662083740702189702033217 <= m - 1
2*m - 1 = 191859282738517986245457499046071864751 <= λ
(m - 1) + (h - 1) = 159882735615431655206469451712737965591 <= λ - 1
```

`|S| = h*(4*m - h - 2)` equals

```
20449991318299214286984618798192526676941391425948511907370273730062895498845
```

and `ceil(39r/100)` equals

```
20449991318299214286984618798192526676699315475205778750815426892976046661961
```

so `|S|` exceeds the threshold by `242075950742733156554846837086848836884`.
In particular `|S|/r ≥ 39/100`. The reference price is `C = 2000` scalar
AArch64 instructions per batched incomplete affine G1 addition under the 6M
convention.

The ratio `h/m` sits at `2/3`. On this side of the square `h = m-1` the same
four units cover `h*(4*m-h-2)` residues, and `(4t-t²)/(1+t)²` at `t = 2/3`
equals `4/5`. The square specialises the same expression:
`(m-1)*(4m-(m-1)-2) = (m-1)*(3m-1)`, whose density in `N = 2m-1` is the
strictly smaller `3/4` limit. The chains below are the shorter pair.

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

Baby list, m records:  (canonical comp of the orbit of B_a, index a)
Giant list, h records: (canonical comp of the orbit of G_b, index b)

LSD radix-sort each list on the 48-byte key: 12 passes of a 32-bit digit,
2^32 buckets, 128-bit counters, 128-bit indices. The sorts are sequential and
reuse one counter array.
Merge by canonical key. A hit means the orbits are equal, so
  x - b*λ ≡ σ * a (mod r)
for some unit σ in {±1, ±λ, ±(1+λ)}.

Resolve the hit without retaining affine babies:
  recompute B = [a]P0 and G = T - [b]Q by double-and-add,
  test the six units by comparing x-coordinates
  (two multiplications by β, β², and a sign-bit test),
  emit candidate = (b*λ + σ*a) mod r for the unique matching unit.
Verify [candidate]P0 = T by one further double-and-add.
Output the scalar on success, else fail.
```

Nothing is randomised. There are no algorithm coins, no restarts and no
amplification. Baby indices and giant indices live in separate lists. The
record count `N = m + h` is

```
N = 159882735615431655206469451712737965593
```

which is a 128-bit integer strictly less than `2^128`, so every bucket counter
fits in 128 bits. Affine coordinates are discarded after the canonical key is
written. A hit recomputes the two points it needs; those additions are inside
the `2048` reference operations charged below, not a second uncounted chain.

## The explicit set

Every residue in `[0, r-2]` has a unique writing `A + B*λ` with `0 ≤ A < λ`
and `0 ≤ B ≤ λ`. The count of such pairs is `λ(λ+1) = λ² + λ = r - 1`, and the
values fill `[0, r-2]` without repetition: a collision `A + Bλ = A' + B'λ` with
a nonzero `B` gap of size at least 1 forces an `A` gap of size at least `λ`,
which cannot fit inside `[0, λ)`. The missing residue `r-1` is the `b = 0`,
`a = 1` point of `P-` below and is inside `S`.

Let

```
S = P+ ∪ P- ∪ Q+ ∪ Q-
```

with the four blocks defined next. The six-unit orbit match also recovers
residues that use `σ = ±λ`. Those extra residues, if any, only enlarge the
recoverable set. The success lower bound uses `|S|` alone.

**`P+`.** Pairs `(a, b)` giving `a + bλ` with `0 ≤ a < m` and `0 ≤ b < h`.
The canonical coordinates are `(A, B) = (a, b)`. The constraint
`(m-1)+(h-1) ≤ λ-1` keeps `B ≤ λ-1` and the integer value at most
`(m-1)+(h-1)λ < r`, so there is no modular wrap. This block contributes `m*h`
residues.

**`P-`.** Residues `-a + bλ` with `1 ≤ a < m` and `0 ≤ b < h`.

- If `b ≥ 1`, then `-a + bλ = (λ - a) + (b-1)λ`, so `(A, B) = (λ-a, b-1)` with
  `A ∈ [λ-(m-1), λ-1]` and `B ∈ [0, h-2]`. The constraint `2m-1 ≤ λ` is
  `λ-(m-1) ≥ m`, so this `A` is at least `m` and does not meet `P+`, whose `A`
  is at most `m-1`.
- If `b = 0` and `a = 1`, the residue is `r-1`, outside `[0, r-2]`.
- If `b = 0` and `a ≥ 2`, then `r - a = (λ+1-a) + λ*λ`, so
  `(A, B) = (λ+1-a, λ)`. Here `B = λ` meets neither `P+` (`B ≤ h-1 ≤ λ-1`) nor
  the `b ≥ 1` rows (`B ≤ h-2`).

So `P+ ∪ P-` has `h*(2m-1)` distinct residues: `a = 0` once per giant index,
and both signs for each `a ∈ [1, m)`.

**`Q+`.** Residues `(1+λ)a + bλ = a + (a+b)λ` with `1 ≤ a < m`, `0 ≤ b < h`,
and `a+b ≥ h`. Canonical coordinates `(A, B) = (a, a+b)`. The same sum
constraint keeps `B = a+b ≤ (m-1)+(h-1) ≤ λ-1`. The condition `a+b ≥ h` puts
`B ≥ h`, so the point is not in `P+`. Its `A = a ≤ m-1` is disjoint from `P-`,
whose positive-`b` rows have `A ≥ m` and whose `b = 0` rows have `B = λ` or are
the single residue `r-1`.

The count splits on `h`, using `h ≤ m-1` so that `a` runs past `h`:

- for `1 ≤ a ≤ h`, the legal giants with `b ≥ h-a` are exactly `a` indices;
- for `h < a ≤ m-1`, every `b ∈ [0, h)` satisfies `a+b ≥ h`, which is `h` indices.

```
|Q+| = Σ_{a=1}^{h} a  +  (m-1-h)*h
     = h(h+1)/2 + (m-1)h - h²
     = (m-1)h - h(h-1)/2.
```

`h(h-1)` is even, so the half is an integer. Different pairs produce different
`(a, a+b)`.

**`Q-`.** Residues `-(1+λ)a + bλ` with `1 ≤ a < m` and `0 ≤ b < min(a, h)`.
The restriction `b < a` is required: if `b ≥ a` then

```
-(1+λ)a + bλ = -a + (b-a)λ
```

which is the `P-` residue for the same `a` and giant index `b-a < h`. Those
points are not counted again.

For `a ≥ 2` the integer identity

```
(λ - a + 1) + (λ - a + b)λ = r + (-(1+λ)a + bλ)
```

holds by expanding the left side to `λ² + λ + 1 - a + (b-a)λ`. Write
`k = -(1+λ)a + bλ` and `V = (λ+1-a) + (λ-a+b)λ`, so `V = r + k`. Under
`0 ≤ b ≤ a-1` one has `k < 0`. Also

```
k ≥ -(1+λ)(m-1) ≥ -(1+λ)(λ-1)/2 > -r
```

because `2m-1 ≤ λ` gives `m-1 ≤ (λ-1)/2` and `(1+λ)(λ-1)/2 = (λ²-1)/2 < r`.
Thus `0 < V < r`. For `a ≥ 2` the coefficients satisfy `0 ≤ λ+1-a ≤ λ-1` and
`0 ≤ λ-a+b ≤ λ-1`, so `V` is already the canonical writing of the residue
`k mod r = V`:

```
(A, B) = (λ+1-a, λ-a+b).
```

Here `A = λ+1-a ≥ λ+1-(m-1) = λ-m+2`. The constraint `2m-1 ≤ λ` is
`λ-m+2 ≥ m`, so `A ≥ m`, disjoint from `P+` and from `Q+` (both have
`A ≤ m-1`). And `B = λ-a+b ≤ λ-1`, while `B ≥ λ-a ≥ λ-(m-1)`. The positive-`b`
rows of `P-` have `B ≤ h-2`. The inequality `λ-(m-1) ≥ h` is exactly
`(m-1)+(h-1) ≤ λ-1`, so those rows do not meet. The `b = 0` rows of `P-` have
`B = λ`, and here `B ≤ λ-1`.

For `a = 1` the only legal index is `b = 0` (because `b < min(1, h)`). The
residue is `-(1+λ) ≡ λ² (mod r)`, and the integer `λ² = r-λ-1` lies in
`[0, r)`. Its canonical writing is `(A, B) = (0, λ)`, since `λ*λ = λ²` and
`0 < λ`. This point has `B = λ`, so it meets no row with `B ≤ λ-1`, and it is
not `r-1`.

The same split as `Q+` counts the pairs, including the one point at `a = 1`:

```
|Q-| = (m-1)h - h(h-1)/2 = |Q+|.
```

The map `(a, b) ↦ (λ+1-a, λ-a+b)` is injective on `a ≥ 2`, and the `a = 1`
point has `A = 0`, which those pairs do not.

**Total.** The four blocks are pairwise disjoint, so

```
|S| = h*(2*m - 1) + 2*((m-1)h - h(h-1)/2)
    = h*(2*m - 1) + 2(m-1)h - h(h-1)
    = h*(4*m - h - 2).
```

At `h = m-1` this collapses to `(m-1)*(3*m-1)`, which is the square count.
For the parameters above the product is the integer displayed in the parameter
section, and it is at least `ceil(39r/100)`.

Baby scalars themselves are distinct. The four baby residues `±a` and
`±(1+λ)a` for fixed `a ∈ [1, m)`, together with `0`, are the `b = 0` slice of
the same blocks: `b = 0 < a` keeps every such `a` inside `Q-`, and `h ≥ 1`
makes `b = 0` a legal giant index. The coordinate separation does not use any
larger `h`. Hence the `m` baby canonical keys are pairwise distinct. Giant keys
at different `b` may coincide for a particular `T`; the merge keeps any one
copy. A coincidence means two giants lie in one orbit, and either copy decodes
to the same scalar.

`experiment.py` rebuilds the four blocks inside every group of order
`λ²+λ+1` for toy eigenvalues `6 ≤ λ ≤ 80`, for every rectangle obeying the same
two inequalities, and checks both the cardinality identity and the absence of
cross-block collisions. On the slice `h = min(m-1, floor(2m/3))` it also checks
that orbit-minimum matching recovers every residue of `S`.

## Recovery

If `x ∈ S`, the unique block representation supplies `b`, `a` and
`σ ∈ {1, -1, 1+λ, -(1+λ)}` with `x ≡ bλ + σ a (mod r)`. Then
`G_b = [σ a]P0` lies in the orbit of `B_a`, the canonical keys match, and the
unit test accepts that `σ`. The final check `[x]P0 = T` succeeds.

If a hit occurs for any unit, including `±λ`, the decoded candidate satisfies
the same linear relation and the verification either accepts the true scalar or
rejects a mismatch. Distinct blocks do not produce two different scalars for
one `T`: each hit determines `x` from `T = [x]P0`. At most the six units on a
matched orbit are tested. A rejected verification is a failed attack on that
input, which the success lower bound already treats as a miss. The ledger
charges `2048` reference additions so that the recomputation of `B_a`, the
recomputation of `G_b`, and the final verification are paid for at the
double-and-add ceiling even though each of those scalars is shorter than
256 bits.

The scalar `0` is the `P+` point `a = 0`, `b = 0`, whose baby is the identity
encoding and whose giant is `T` itself. It is inside `S`.

## Cost

A batched incomplete affine addition costs one reference operation. The model
charges other scalar AArch64 work as `R/C` with `C = 2000`, and it forbids
charging both a reference addition and the field multiplications inside it.
Canonicalisation is not part of the addition, so its field multiplications are
instruction work. `blst_fp_mul` on the pinned measurement is `383` instructions.
An inversion-free addition is the reference operation itself.

The baby chain performs `m-1` additions and the giant chain performs `h-1`
additions. Resolving one hit recomputes `[a]P0` and `T-[b]Q`, then verifies
the emitted scalar. A double-and-add of a scalar in `[0, r)` uses at most
`255` doublings and `255` additions, hence at most `512` reference operations.
Three such scalar multiplications are at most `1536`. The ledger charges
`2048 = 4*512`, which covers that resolution and leaves a full extra
double-and-add. No further group operation is performed. There is no restart.

`H = (m-1) + (h-1) + 2048` counts every reference operation once.

Instruction work, outside those additions:

The identity baby `a = 0` stores the infinity encoding: 32 instructions to
materialise the key and index, 840 to sort, 48 to merge, total 920. Every
other stored record pays the following.

| Orbit-minimum selection, per non-identity point | Instructions |
| --- | ---: |
| reduced product `x*β` (`blst_fp_mul`, pinned measurement) | 383 |
| reduced product `x*β²` | 383 |
| compare `x` with `βx`: 6 high-to-low limbs, select | 32 |
| compare the winner with `β²x` | 32 |
| clear the sign bit and address the winner | 0, folded into the store |
| **selection** | **830** |

Each six-limb integer comparison is charged at 32 instructions, so the two
comparisons sum to 64. Materialisation of the winning 48-byte key is the same
32-instruction charge used for a chain record (six stores, addressing, two
index stores). Sort and merge match the 12-pass 32-bit radix accounting:
70 instructions per pass, 12 passes, 840; merge at 48 per record (six loads of
each cursor, six comparisons, pointer, branch).

| Per non-identity stored record | Instructions |
| --- | ---: |
| selection | 830 |
| store canonical key and 128-bit index | 32 |
| 12-pass radix sort | 840 |
| merge | 48 |
| **total `PER_POINT`** | **1750** |

Counting costs 25 instructions per record (address a 48-byte key from a
128-bit index, isolate the digit, load/increment/store a 128-bit counter,
loop). Scattering costs 45 (digit, destination address from a 128-bit counter,
counter update, load and store the 128-bit index, loop). Together with loop
overhead this is the 70 per pass. Bucket maintenance, paid twice because the
lists are sorted sequentially and the counter array is reused, is
`2 * 12 * 2^32 * 6`: clear and prefix-sum `2^32` counters of 16 bytes at 6
instructions each. A fixed 8192 instructions cover the six field comparisons
that read the already recomputed points, candidate decoding, output, and
program state. The two extra `blst_fp_mul` calls inside the unit test are
`766` instructions and sit inside that fixed block.

```
H = (m - 1) + (h - 1) + 2048
  = 159882735615431655206469451712737967639
R = 920 + 1750*(m - 1) + 1750*h + 2*12*2^32*6 + 8192
  = 279794787327005396611321540497909915085736
W = H + R/2000
  = 74945032319733588378032555490423231295467 / 250
log2(W) = 127.8171723641...
```

The declared bound `2^127.817173` is `log2(W)` rounded upward at six decimal
places. Direct comparison of the fraction against the neighbouring powers of
two gives

```
2^127.817172 < W ≤ 2^127.817173.
```

Memory and preprocessing below are rounded upward the same way. No component
of `W` is an expectation, an amortisation over future targets, or a free
oracle.

## Memory

Sorts run one list at a time and permute indices. Both lists are live for the
merge. Scratch indices exist only for the list being sorted. Affine points are
not retained.

| Component | Bytes |
| --- | ---: |
| keys, `N` × 48 | `48*N` |
| indices, `N` × 16 | `16*N` |
| scratch indices, `max(m, h)` × 16 | `16*m` |
| bucket counters, one `2^32` × 16 array, reused | `2^36` |
| code, constants, stacks, fixed scratch | `2^20` |

`h < m`, so the scratch term is `16*m`. Peak storage is

```
64*N + 16*m + 2^36 + 2^20
  = 11767369341295769823177704902052525241280
```

bytes, i.e. at most `2^133.111916` bytes. Memory is mandatory for this
construction; it is not a practical attack.

## Preprocessing and advice

Preprocessing is everything independent of `T`: the baby chain (`m-1`
additions), canonicalisation of the `m` baby orbits, the baby list's 12-pass
sort, and one bucket-maintenance pass. `φ(P0)` is produced while canonicalising
`B_1` and is not a second chain. The giant chain, the giant sort, the merge,
and the hit recomputation read `T`.

```
preprocessing = (m - 1) + (872 + (766+64+32+840)*(m - 1) + 12*2^32*6) / 2000
              = 177565766174498396270170915367294129649217 / 1000
log2 = 127.0616210695...
```

so at most `2^127.061622` target-group-operation equivalents. The baby-sort and
baby-selection terms are the baby portion of `R` and are not added again on top
of `W`. The identity baby's 872 instructions are a 32-instruction store plus
the 840-instruction sort. `766+64+32+840 = 1702` is `PER_POINT` without the
48-instruction merge, which reads both lists and is online. Nonuniform advice
is zero bytes: `β`, `β²`, `λ` and `P0` are public profile constants, the baby
table is rebuilt for the attack, and its construction is fully charged.

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

The `±λ` associates of a baby are orbits of points already represented in
`{0, ±a, ±(1+λ)a}`. They are not given a separate block in `|S|`. Counting
them again would double-count the `A = 0` ray and the `B = λ` column. The
success lower bound does not use them.

## Submitted bound

| Quantity | Value |
| --- | ---: |
| total work `W` | `74945032319733588378032555490423231295467/250` |
| `time_log2` | **127.817173** |
| peak memory | `11767369341295769823177704902052525241280` bytes, `2^133.111916` |
| preprocessing | `177565766174498396270170915367294129649217/1000`, `2^127.061622` |
| nonuniform advice | 0 bytes |
| success probability | at least `|S|/r ≥ 39/100` |

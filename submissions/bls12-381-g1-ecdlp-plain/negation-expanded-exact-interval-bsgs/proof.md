# bls12-381-g1-ecdlp-plain: negation-expanded exact-interval baby-step/giant-step

Deterministic and hash-free. Baby points `j*P0` for `j` in `[0, m)` are stored
under both compressed encodings. The negative key is the positive key with one
flag bit flipped, so it is not a second group addition. Giant steps subtract
multiples of stride `s = 2m-1`. Because the 48-byte compressed encoding is
injective on G1, a cross-list key hit is point equality. The recoverable set is
exactly two explicit integer intervals, their cardinalities sum to `h*s`, the
success probability is exactly `h*s/r`, and at most one candidate reaches
verification.

Every cost in the protected `classical-group-work-v2` model is charged,
including multiword counters, multiword scatter offsets, record loads and
stores, the sign-flag copy, and the batching construction behind the 6M
affine-addition price.

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
m = 87990032694955804185480753244926427368      (127 bits, m < 2^127)
h = 116206294576542127460944719508070012830     (127 bits, h < 2^127)
s = 2*m - 1 = 175980065389911608370961506489852854735
h*s = 20449991318299214286984618798192526676726300062321169823908062635023576250050
```

`h*s` is the least multiple of `s` that is at least `ceil(39r/100)`, at this
`m`. It satisfies `ceil(39r/100) <= h*s <= r`. Also `2(m-1) < r` and `s < r`.
The reference price for this track is `C = 2000` scalar AArch64 instructions
per target-group operation, one operation being a batched incomplete affine G1
addition under the 6M convention.

## Algorithm

`comp(P)` is the 48-byte canonical compressed encoding from the pinned
`py_ecc` serializer: the big-endian `x`-coordinate with the three flag bits of
the first byte. The flag of value `0x20` is set iff `y` is the lexicographically
larger square root of `x^3+4`. The identity uses the reserved infinity encoding.

```text
Baby chain: B_j = j*P0 for j = 0..m-1          # m-1 batched affine additions
S = 2*B_{m-1} + P0                             # one doubling and one addition
Giant chain: G_i = T - i*S for i = 0..h-1      # h-1 batched affine additions

Baby list, 2m-1 records:
  (comp(B_j), index = j) for j = 0..m-1
  (comp(B_j) with the 0x20 flag flipped, index = 2^127 + j) for j = 1..m-1
Giant list, h records:
  (comp(G_i), index = i) for i = 0..h-1

LSD radix-sort each list on the 48-byte key: 12 passes of a 32-bit digit,
2^32 buckets, 128-bit counters, 128-bit indices. The two sorts are sequential
and reuse one counter array.
Merge the sorted lists by key. A hit decodes j from the low 127 bits of the
baby index, the sign from bit 127, and i from the giant index.
Candidate x = (i*s + e) mod r with e = +j or e = -j according to the sign.
Verify [x]P0 = T and output x on success, else fail.
```

Nothing is randomised. There are no algorithm coins, no restarts and no
amplification. Bit 127 of a baby index is free because `j <= m-1 < 2^127`.
Giant indices need no tag bit because they live in the other list, and
`i <= h-1 < 2^127`, so every index fits in 128 bits.

## Exactness of the recoverable set

Lemma 1 (compression is injective, and negation flips one flag). For an affine
point `(x, y)` of `y^2 = x^3 + 4` over `Fp`, `comp` stores `x` and the sign
flag. The two roots `y` and `p-y` are distinct in G1: `y = 0` would be a point
of order 2, and `|G1| = r` is an odd prime. The flag of value `0x20` selects
the lexicographically larger root, i.e. it is set iff `2y > p` when `y` is
taken in `[0, p)`. Negation sends `y` to `p-y` and leaves `x` unchanged, so it
flips that comparison. The `x`-coordinate is at most `p-1 < 2^381`, so the top
three bits of the 384-bit encoding are zero before flags are written; flipping
`0x20` does not collide with the coordinate. The identity has the reserved
infinity encoding, fixed by negation, and is not given a second record. Hence
`comp(P) = comp(Q)` iff `P = Q`, and for every non-identity `P`,
`comp(-P)` is `comp(P)` with the `0x20` bit flipped.

Lemma 2 (keys are distinct inside each list). On the baby list the stored
points are `e*P0` for `e` in `[-(m-1), m-1]`. If `e*P0 = e'*P0` then
`e ≡ e' (mod r)`. Both representatives lie in an interval of length
`2(m-1) < r`, so `e = e'`. On the giant list, `T - i*S = T - i'*S` with
`0 <= i, i' < h` gives `(i-i')*s ≡ 0 (mod r)`. Since `r` is prime and
`0 < s < r`, `s` is invertible, so `i = i'`. A giant point is never the
negative of a different giant in a way that would identify keys: keys use the
full compressed encoding, so `G_i = -G_{i'}` is a different key from `G_i`
whenever `G_i` is not the identity. The identity giant occurs only for
`x = i*s`, which matches the single baby record `j = 0`.

Lemma 3 (a key match is one scalar). A merge hit means `G_i = e*P0` for
`e = ±j` as decoded from the baby index. That is
`T - i*s*P0 = e*P0`, so `x ≡ i*s + e (mod r)`, with `e` in `[-(m-1), m-1]`.

Lemma 4 (the parameter rectangle tiles two intervals, bijectively). Consider
the `h*s` integers

```
V = { i*s + e : 0 <= i < h,  -(m-1) <= e <= m-1 }.
```

Each block `i` has `s = 2m-1` consecutive values, from `i*s-(m-1)` to
`i*s+(m-1)`. The next block starts at `(i+1)*s-(m-1) = i*s+m`, which is
exactly one past the end of block `i`. The blocks therefore partition the
integer interval

```
[-(m-1),  -(m-1) + h*s - 1] = [-(m-1), h*s - m].
```

That interval contains `h*s` integers, and `|V| = h*s`, so each integer in it
has exactly one representation `(i, e)`. Reducing modulo `r`, using
`h*s <= r`: the nonnegative part `[0, h*s-m]` stays in range, and the
`m-1` negative integers `[-(m-1), -1]` land in `[r-(m-1), r)`. These two
residues ranges are disjoint because

```
(h*s - m) < r-(m-1),
```

which rearranges to `h*s < r+1` and is the same constraint. Their cardinalities
sum to `(h*s-m+1) + (m-1) = h*s`. Therefore the recoverable residues are
exactly

```
I0 = [0, h*s - m] = [0, 20449991318299214286984618798192526676638310029626214019722581881778649822682]
I1 = [r-(m-1), r) = [52435875175126190479447740508185965837602562467832682018418177946693654757146, r)
```

and every such residue has one representation. In particular `x = 0` is
`(i, e) = (0, 0)`, and `x = r-1` is `(i, e) = (0, -1)`, which is inside `I1`
because `m > 1`. There is no third interval and no missing interior point.

Theorem (exact success set and probability). For `T = [x]P0` with `x` in
`[0, r)`:

* if `x` is in `I0 ∪ I1`, Lemma 4 supplies a unique `(i, e)`, the giant key
  `comp(T - i*S)` equals the baby key `comp(e*P0)`, and the merge emits that
  pair;
* if `x` is outside `I0 ∪ I1`, a hit would be a representation of `x` as
  `i*s+e` with `(i, e)` in the rectangle, contradicting Lemma 4.

Hence the recoverable set is exactly `I0 ∪ I1` and, for uniform `x`,
`Pr[success] = (h*s)/r`. The fraction is at least `39/100` because
`h*s >= ceil(39r/100)`. The excess is `h*s - ceil(39r/100) = 26984587115391073092635742047529588089`.

Candidate count. Two distinct hits would be two representations of the same
residue, or hits at two different residues. The merge compares one baby list
against one giant list and, by Lemmas 2 and 4, at most one key is shared.
The attack therefore performs at most one verification scalar multiplication.
Every emitted candidate is checked with `[x]P0 = T`. Failure returns nothing.

## Batch-inversion schedule behind the 6M price

The baby chain adds `P0` repeatedly and the giant chain adds `-S` repeatedly.
Both are dependent affine chains. A standalone affine addition needs a field
inversion, so the reference rate is applied only with a batching construction
in front of it.

Block size `B = 2^20`. For a block of `B` chained additions, accumulate the
denominators, form the prefix products (`B-1` field multiplications), invert
the product once, and recover each denominator inverse (`2B` multiplications).
The batch costs `3M` per addition on top of the `3M` of the incomplete-affine
formula itself: the `6M` reference price.

One inversion by binary powering `p-2` is at most 571 base-field
multiplications, each 383 AArch64 instructions in the pinned blst measurement,
i.e. 218,693 instructions. Amortised over `B = 2^20` additions that is 0.209
instructions per addition, inside the 1,149 instructions the 3M batch allowance
already grants. No separate inversion term is added.

The stride point `S = 2*B_{m-1} + P0` is one doubling of an affine baby point
and one addition of `P0`, charged as two reference operations. It is not a
127-bit scalar multiplication: `B_{m-1}` is the last point of the baby chain.
The verification `[x]P0` is charged as 512 reference additions, above the
binary double-and-add bound `2*255` for a scalar in `[0, r)`.

## Work ledger

Units. `H` counts target-group operations. `R` counts other scalar AArch64
instructions. Total work is `W = H + R/C` with `C = 2000`. Record counts:

```
N_baby = 2*m - 1
N_giant = h
N = N_baby + N_giant = 292186359966453735831906225997922867565
```

`N < 2^128`, so the loop index used to address records is 128 bits wide, the
same width as the stored index.

Sort, per list, LSD radix, 12 passes, `2^32` buckets, 128-bit counters
(16 bytes: a bucket count can reach `N`), 128-bit scatter offsets. Each pass
is one counting traversal and one scatter traversal. The instruction tables
below are unchanged from a single-list 32-bit-digit sort; each record of each
list pays them once.

| Counting traversal, per record | Instructions |
| --- | ---: |
| form key address: 3-word base + 128-bit index × 48 | 12 |
| load 32-bit digit, isolate by shift and mask | 5 |
| load 128-bit bucket counter (2 loads) | 2 |
| increment counter with carry | 2 |
| store 128-bit counter (2 stores) | 2 |
| loop pointer update and branch | 2 |
| **count total** | **25** |

| Scatter traversal, per record | Instructions |
| --- | ---: |
| form key address and load the digit | 17 |
| load 128-bit bucket counter | 2 |
| form destination address, counter × 16, 128-bit | 16 |
| increment counter with carry and store it | 4 |
| load source 128-bit index (2 loads) | 2 |
| store destination 128-bit index (2 stores) | 2 |
| loop pointer update and branch | 2 |
| **scatter total** | **45** |

One pass is `25+45 = 70`. Twelve passes are `840` per record.

Chain-produced records (the `m` baby points and the `h` giant points) are
materialised by six key stores plus addressing and two index stores: `32`
instructions. Each of the `m-1` negated baby records is a copy of an existing
48-byte key with the sign flag flipped:

| Negated-key materialisation | Instructions |
| --- | ---: |
| source key address (128-bit index × 48) | 12 |
| load 48-byte key (6 loads) | 6 |
| flip the `0x20` flag | 1 |
| destination key address | 12 |
| store 48-byte key (6 stores) | 6 |
| store 128-bit index (2 stores) | 2 |
| loop update and branch | 2 |
| **charged** | **48** |

The itemisation sums to 41; the ledger charges 48. Relative to the 32-instruction
materialisation already in the per-record total below, each negated record
costs 16 extra instructions.

The merge walks both sorted lists once. Each consumed record pays at most six
loads of its key, six loads of the other cursor's key, six comparisons, a
pointer update and a branch: charged at 48 per record, which covers that
itemisation. There is no separate within-list equality scan; Lemma 2 says it
would find nothing.

Per chain-produced record: `840 + 32 + 48 = 920`.
Per negated record: `920 + 16`.
Bucket maintenance is paid twice (one clear-and-prefix per list): clearing and
prefix-summing `2^32` 128-bit counters is at most 6 instructions each, so
`2 * 12 * 2^32 * 6` instructions. A fixed 4,096 instructions cover candidate
output and program state. The counter array is allocated once and reused.

```
H = (m-1) + (h-1) + 2 + 512
  = 204196327271497931646425472752996440710
R = 920*N + 16*(m-1) + 2*12*2^32*6 + 4096
  = 270219291692256729832321419970626336292392
W = H + R/2000
  = 84826493279406574140646545684577402214049 / 250
log2(W) = 127.99585443309058
```

The declared bound `2^127.995855` is `log2(W)` rounded upward at six decimal
places. Memory and preprocessing below are rounded upward the same way. No
component of `W` is an expectation, an amortisation over future targets, or a
free oracle. The same numerator at `m±1`, with `h` retargeted to the least
integer whose coverage still clears `ceil(39r/100)`, is strictly larger; the
submitted pair is a local minimum of this ledger, and the bound does not depend
on uniqueness.

## Memory

Sorts run one list at a time. Keys are not copied; the sort permutes indices.
Both lists' keys and indices are live together because the merge needs them.
Scratch indices exist only for the list being sorted.

| Component | Bytes |
| --- | ---: |
| keys, `N` × 48 | `48*N` |
| indices, `N` × 16 | `16*N` |
| scratch indices, `N_baby` × 16 (`N_baby > N_giant`) | `16*(2m-1)` |
| bucket counters, one `2^32` × 16 array, reused | `2^36` |
| code, constants, stacks, fixed scratch | `2^20` |

Peak storage is

```
64*N + 16*(2m-1) + 2^36 + 2^20
  = 21515608084091624827177382567773429725232
```

bytes, i.e. at most `2^133.982508` bytes. Memory is mandatory for this
construction; it is not a practical attack.

## Preprocessing and advice

Preprocessing is everything independent of `T`: the baby chain (`m-1`
additions), the stride double-and-add (2 additions), materialisation of the
`2m-1` baby records (32 per chain point, 48 per negated copy), the baby list's
own 12-pass radix sort (`840` per baby record) and one bucket-maintenance pass.
The giant chain, the giant sort and the merge read `T` and are not preprocessing.
In the model's units

```
preprocessing = (m-1 + 2) + (872*(2m-1) + 16*(m-1) + 12*2^32*6)/2000
              = 41355315366629227967175954025154075568763 / 250
log2 = 126.95941228062183
```

so at most `2^126.959413` target-group-operation equivalents. The baby-sort
term is the baby portion of `R` and is not added again on top of `W`.
Nonuniform advice is zero bytes: the baby table depends only on the public
generator, is rebuilt on every attack, and is fully charged. Nothing is stored
between targets.

## Determinism, heuristics and limitations

There are no heuristic premises: no random-oracle model, no hash-load-factor
assumption, no probe-window model, no cycle-finding heuristic, and no
smoothness or linear-algebra heuristic. Both sorts are deterministic worst-case
radix sorts. The merge is exact key equality. The success statement is an
identity over the uniform challenge distribution, including the two-interval
shape forced by reduction of the negative baby offsets. The only imported
constants are the protected profile (`r`, `P0`, the curve, the pinned
compressed encoding), the reference price `C = 2000`, the 6M affine convention
and the measured 383 instructions per base-field multiplication.

The accounting is analytic. No full-size execution is claimed or needed.
`experiment.py` checks the coverage identity on small additive groups only and
is not evidence about production encodings or instruction counts. The
exponential tables are required by this construction; the attack is not claimed
to be practical, novel, or optimal among generic algorithms.

## Submitted bound

| Quantity | Value |
| --- | ---: |
| total work `W` | `84826493279406574140646545684577402214049/250` |
| `time_log2` | **127.995855** |
| peak memory | `21515608084091624827177382567773429725232` bytes, `2^133.982508` |
| preprocessing | `41355315366629227967175954025154075568763/250`, `2^126.959413` |
| nonuniform advice | 0 bytes |
| success probability | exactly `h*s/r >= 39/100` |

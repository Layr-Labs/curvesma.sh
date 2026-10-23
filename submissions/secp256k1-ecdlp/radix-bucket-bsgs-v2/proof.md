# secp256k1-ecdlp: compact-record radix-bucket BSGS truncated at 39/100

Status: ready for exploratory review. No qualification is claimed by this file.
The normative input and output are `production-v2`. `claim.json` binds
`tracks/secp256k1-ecdlp/target-profile.json` and
`tracks/classical-group-work-v2.json` by the SHA-256 values recorded there.
The score uses that cost model's rule `W = H + R/C` with the fixed secp256k1
price `C = 1000`. One reference operation is one batched incomplete affine
addition (the 6M convention, `S = M`, `A = 0`). Every other charged action is
a scalar AArch64 instruction equivalent. Reference operations are not also
expanded into their internal multiplications.

No negation map, endomorphism, hash function, random walk or uniformity
hypothesis is used. Bucket populations do not appear in the cost. Every
counter array is scanned in full. The change against the 96-byte-record
construction is a 64-byte record layout with an explicit sign byte, a single
linear comparison scan, and emission in place; the coverage argument and the
success set are unchanged.

## Parameters

Let `r` be the protected subgroup order and `P` the protected generator.
Let `Q = [x]P` with `x` uniform in `0 .. r-1`.

```text
S = ceil(39*r / 100) = (39*r + 99) // 100
m = 212506270031153001307260662621714878127
```

Integer relations used below, all in ordinary arithmetic:

- `(m-1)^2 < S <= m^2`, so `m = ceil(sqrt(S))` and `ceil(S/m) = m`.
- `m^2 < r`.
- `100 * m^2 >= 39 * r`, hence `m^2 / r >= 39/100`.

`K = 2^40` is the batch width. `D` denotes the point `-[m]P`.

## Procedure

Reject a malformed or off-curve input and spend the same budget derived below
(the failure path is not cheaper in the bound). If `Q` is the point at
infinity, return `0`.

**Record layout.** A record is 64 bytes: a 32-byte big-endian `x`, a 16-byte
index, a 1-byte sign (`0` when the stored `y` is even, `1` when it is odd) and
15 bytes of padding. The radix key is the 256-bit `x`; the sign byte
disambiguates the two curve points that can share an `x`. Records are
8-byte aligned, so a copy is four `ldp`/`stp` pairs and an address is
`base + index*64`.

**Baby list.** Build a record `([j]P, j, sign)` for every `j = 1 .. m-1`, and
build the stride `D`, as follows.

- Seed chain: `[1]P, ..., [K-1]P` and `[K]P` by adding `P`. These `K`
  additions are dependent, so they are not charged as reference operations.
  Each is charged as one unbatched inversion-plus-formula at `2^17`
  instructions (see the inversion paragraph).
- Let `Sigma = [K]P`. For each residue `rho = 0 .. K-1` whose first index is
  still `< m`, walk the arithmetic progression `rho, rho+K, rho+2K, ...`
  while the index is `< m`. Each step adds `Sigma` to the running point.
  Across residues, one round is `K` independent additions.
- The scalar of every point in these progressions lies in `0 .. m-1`.
  For two such scalars `a, b` with `0 <= a < b < m < 2^128` and `r > 2^255`,
  `b != a`, `b != r-a` and neither point is infinity. The addition is
  therefore the incomplete affine sum of two finite points with distinct `x`.
  Each such step is one reference operation. Producing the `m-1` finite
  babies and the extra stride point uses at most `m` reference additions.
  The bound charges `m` of them.
- Each generated point is emitted into the baby buffer in place; no second
  copy is made at emission time.

**Giant list.** `Gamma_0 = Q`. For `i = 1 .. m-1`, `Gamma_i = Gamma_{i-1} + D`,
again in `K` parallel residue chains, so each ordinary step is one batched
incomplete affine addition. Charge `m` reference additions on this side as
well (`H = 2m`). While a chain is generated, test the infinity flag. If
`Gamma_i` is infinity, the only scalar consistent with that event inside the
segment is `a = i*m`; keep it as the candidate. Otherwise emit the record
`(x, i, sign)` into the giant buffer.

A giant step is exceptional when the current point is infinity, equals `D`,
or equals `-D`. Each residue chain hits each of those at most once. Charge
`4K` exceptional steps (both sides, with slack) at `2^17` instructions each,
instead of as reference operations. The ordinary steps remain reference
operations. This split does not bill a reference operation and its inversion
together.

**Ledger.** LSD counting sort on the 256-bit `x`, four passes, digits of
64 bits. Counters are 128-bit integers (`m < 2^128`). The counter array has
`2^64` entries and `2^68` bytes. It is zeroed in full and prefix-summed in
full on every pass; empty digits are not skipped. Scatter writes a second
64-byte buffer. Sort the baby records and the giant records.

**Comparison scan.** Walk the two sorted arrays once with one pointer each,
comparing the 256-bit `x` and then the sign byte. A record is compared at
most once per array position, so the scan performs `2m` record comparisons.
Equal `x` runs have length at most 2, because a curve `x` has at most two
points and both lists contain distinct points (`j` in `1 .. m-1` are pairwise
distinct, and `i*m` steps of a non-identity point stay distinct for
`i < m < r`). On equal `x` and equal sign, the points are equal; the giant `i`
and baby `j` yield the candidate `a = i*m + j`. On equal `x` with unequal
sign the points are negatives and the pair is not a match.

**Always** run every pass, every prefix word and the whole scan. A match does
not stop the bound early. Then perform one verification scalar multiplication
`[a]P` and return `a` when it equals `Q`. If there is no candidate,
verification is still charged and the procedure returns failure.

## Success

Take an integer `a` with `0 <= a < m^2`. Write `a = i*m + j` with
`i = floor(a/m)` and `j = a - i*m`. Then `0 <= i < m` and `0 <= j < m`.

- If `j = 0`, then `Q = [i*m]P`, so `Gamma_i = Q + i*D = infinity`. The
  infinity test emits `a`.
- If `j >= 1`, then `Gamma_i = [j]P`, which is stored with the same `x` and
  the same sign. The comparison selects that unique record.

No second pair `(i', j')` in the same ranges satisfies `i'*m + j' = a`. A
match cannot be the negative baby: the sign byte differs and `r - j` is not an
index in `1 .. m-1`. Because `a < m^2 < r`, reduction modulo `r` does not
create another representative. The verification equality accepts this `a` and
would reject any other integer in `0 .. r-1`.

The set of such `a` has cardinality `m^2`. For `x` uniform in `0 .. r-1` the
success probability is exactly `m^2 / r`, which is at least `39/100 = 0.39`.
The algorithm coins are empty. Work does not depend on retrying. The declared
success probability is this lower bound, not `m^2 / r` itself and not 1.

The same decomposition on the additive group `Z/nZ` is what `experiment.py`
checks for a caller-supplied order at most `10^4`. That run is not a
production cost measurement.

## Instruction prices

Field multiplication is the cost-model figure of 179 scalar instructions.
A non-batched inversion is charged at `2^17 = 131072` instructions. Binary
powering by `p-2` is at most 511 multiplications, `511 * 179 = 91469`, and
the formula multiplications around one inversion fit in the same pad.
Jacobian doubling or mixed addition used only in the final verification is
charged at 4096 instructions each, against a formula under 16 multiplications
(`16 * 179 = 2864`).

One counting-sort pass, per 64-byte record, is charged at 80 instructions. A
direct schedule of the body is 58:

| Piece | Instructions |
| --- | ---: |
| 3-limb record-pointer bump and loop branch | 12 |
| 64-bit digit load and mask | 2 |
| 2-limb counter address `base + digit * 16` | 8 |
| 128-bit counter load, add-with-carry, store | 8 |
| 2-limb `index * 64` destination address | 12 |
| 64-byte record copy (4 `ldp` + 4 `stp`) | 8 |
| sign byte and index save, branch slack | 8 |

The remaining 22 instructions are unused pad. Prefix scan is charged at 64
instructions per bucket (load, 128-bit add, store, 2-limb pointer, branch,
pad) against a body of about 12. Zeroing a 16-byte counter is charged at 16
instructions. The comparison scan is charged at 48 instructions per record
(4-limb `x` compare, sign compare, two 3-limb pointer bumps, pad) against a
body of about 30. Emitting a generated record in place is charged at 24
instructions (64-byte store, index bump, loop branch) against a body of about
14.

There are 4 passes and 2 lists, so 8 full prefix scans and 8 full zeroings.
Record counts are bounded by `m` babies and `m` giants, so the scan and the
emission each touch `2m` records.

## Totals

Let `B = 2^64`.

| Component | Count |
| --- | ---: |
| Reference additions `H = 2m` | 425012540062306002614521325243429756254 |
| Radix bodies `(2m)*4*80` | 136004012819937600836646824077897521681280 |
| Prefix `8*B*64` | 9444732965739290427392 |
| Zeroing `8*B*16` | 2361183241434822606848 |
| Comparison scan `(2m)*48` | 20400601922990688125497023611684628300192 |
| Emission `(2m)*24` | 10200300961495344062748511805842314150096 |
| Batch-inversion residuals | 50665442951000452353301206794240 |
| Exceptional steps `4K*2^17` | 576460752303423488 |
| Unbatched seeds `2K*2^17` | 288230376151711744 |
| Verification and decode | 2109440 |

Batch-inversion residuals: `ceil(m/K)` rounds on each side, plus 4 spare
inversions, each `2^17` instructions. The 6M reference price already contains
the three Montgomery multiplications per element and does not contain the
single inversion per round, so the inversion is an extra `R` term, not a
second copy of the reference operation. `ceil(m/K) = 193273326686860856450276208`.

`R` is the sum of every instruction row:

```text
R = 166604915755089395987699592747028241524720
1000*H + R = 591617455817395398602220917990457997778720
```

`W = (1000*H + R) / 1000`. Then
`log2(W) = 128.797932225853...`. Rounded upward at the sixth decimal place
this is **128.797933**, the rational **128797933/1000000**.

## Preprocessing and memory

Preprocessing is the baby side only: `m` reference additions, four radix
passes on `m` records, four prefix scans, four zeroings, `m` in-place
emissions, `ceil(m/K)` residual inversions and the `K` unbatched seed
additions. Giants depend on `Q` and are not preprocessing. They are inside
`W`.

```text
R_pre = 73102156916049353931101140452483729439528
1000*m + R_pre = 285608627127202355238361803074198607566528
```

`log2` of that quotient is `127.747306152...`, rounded upward to
**127.747307**. Preprocessing is included in `W`, not added again.

Peak storage, in bytes: each list holds `m` live records plus one `m`-record
scatter buffer, four buffers of 64 bytes per record in total, `4*m*64 =
54401605127975168334658729631159008800512`. The counter array is `2^68`
bytes, which is below one buffer and is already inside the same total. Code,
pointers, the verification state and retained randomness fit in the `2^20`
bytes of slack in the rounded bound. Peak memory is at most
`2^135.320774` bytes. Nonuniform advice is 0 bytes; nothing is precomputed
from `Q` or from the subgroup beyond the public profile constants.

## Reproducibility

The toy experiment is deterministic and uses no randomness. The production
bound is the integer ledger above; it is checked by arithmetic, not by
execution.

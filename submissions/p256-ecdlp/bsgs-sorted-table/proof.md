# p256-ecdlp: deterministic baby-step/giant-step with a merge-sorted affine table

Status: ready for exploratory review; no qualification or emitted score is claimed
here. This is an analytic generic accounting bound on the full production target,
not a practical break and not a claim of novelty.

Normative inputs and outputs are those of the protected `target-profile.json`
(profile SHA-256 `802d2350d1804b19fcf0e1f6d80d816862e966625710abfa2b4244ca764aa605`)
and the accounting convention of `classical-group-work-v2`
(cost-model SHA-256 `e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294`).
No restricted challenge distribution, auxiliary powers, advice, random-oracle
premise or quantum step is used.

## Claim

For uniform `x` in `[0,r)`, given `Q = [x]P` for the fixed SEC 2 section 2.4.2
generator `P` on NIST P-256, the deterministic baby-step/giant-step algorithm below
recovers `x` with probability **1** (no randomness, no retries, no amplification)
using at most

```
T = 12489723995466125362959701591175620472619008
```

scalar AArch64 instruction equivalents, including preprocessing, sorting, every
lookup, all failure paths, verification and output. With the fixed normalization
`C = 1000` the submitted bound is `time_log2 = 133.197866` (`66598933/500000`,
rounded upward at the sixth decimal). Preprocessing — the entire baby table and
its sort — is `8994343522454245466263917623636497365729280` instruction equivalents
(`preprocessing_log2 = 132.724214`) and is included in `T` with no reuse credit.
Peak memory is `32667107224410092492483962313449749348352` bytes
(`memory_log2_bytes = 134.584963`). Nonuniform advice is 0 bytes.

## Parameters

`p`, `r`, `a = -3`, `b` and `P` are taken verbatim from the protected profile.
Set `bitlen = bit_length(r) = 256`, `L = floor(bitlen/2) = 128`, `m = 2^L = 2^128`
and `h = ceil(r/m) = 2^128`. The last equality is exact: the high 128 bits of `r`
are `0xffffffff00000000ffffffffffffffff = 2^128 - 1`, so
`r = (2^128 - 1)*2^128 + r_lo` with `0 < r_lo < 2^128`, hence
`ceil(r/2^128) = 2^128`. Both `m` and `h` are exact powers of two, which fixes the
binary-search depth at exactly `L = 128` comparisons and removes any dependence of
the lookup cost on data.

All table indices, `i*m + j` and the record fields fit 256 bits for this fixed
target (`i*m + j < 2^256`).

## Table record and key

A record is 6 words = 48 bytes: a 1-bit infinity tag, the affine `x` coordinate
(4 words), the parity bit of the affine `y` coordinate, and the baby index `j`
(2 words, 17 bits used). Records for the point at infinity carry tag 1 and zero
coordinate words.

The key is `(tag, x, parity(y))` in lexicographic order. `p` is odd, so `y` and
`p - y` always have opposite parity; the single parity bit therefore distinguishes
`[j]P` from `-[j]P` exactly, and no 256-bit `y` coordinate needs to be stored or
compared. Lookup comparison touches the tag word, then up to four `x` limbs, then
one parity word.

## Construction

```text
R = infinity (Jacobian);  T = array of m records;  buffer = empty
for j = 0 .. m-1:
    if j > 0:  R = R + P                 # mixed Jacobian+affine addition
    buffer.append((R, j))
    if len(buffer) == 1024: flush(buffer) # batch-invert, convert, emit records
R = R + P                                # now R = [m]P
flush(buffer)
D = (x([m]P), p - y([m]P))               # affine, fixed
sort T by (tag, x, parity(y))            # bottom-up merge sort, L passes, 2 buffers
gamma = Q
for i = 0 .. h-1:
    if i mod 1024 == 0: pre-convert the next 1024 giant steps
    found, j = binary_search(T, key(gamma))
    if found:
        a = (i*m + j) mod r
        if [a]P == Q: return a           # independent attack-side check
    gamma = gamma + D                    # mixed addition
return failure                           # unreachable on valid input
```

`flush` performs Montgomery batch inversion over the 1024 buffered `Z`
coordinates: one field inversion, `3*(k-1)` field multiplications to form the
individual inverses, and `k` conversions `x = X*Z^-2`, `y = Y*Z^-3`. The whole
block is charged, including the final partial block (charged as a full block).

Both walks use the same primitive: a Jacobian accumulator plus a fixed **affine**
addend (`P` for the baby walk, `D` for the giant walk), so no Jacobian-Jacobian
formula is needed and no inversion is performed per addition.

## Coverage, identity and success

Fix a valid input. If `Q` is the point at infinity the algorithm returns `0`
immediately, which is the required scalar. Otherwise there is a unique `x` in
`[0,r)` with `Q = [x]P`. Write `x = i*m + j` with `0 <= j < m` and `0 <= i < h`;
this decomposition exists and is unique because `x < r <= m*h`.

At giant step `i`, `gamma_i = Q - [i*m]P = [x - i*m]P = [j]P`, which is a point of
the baby table. Its key is exactly the key stored for index `j`, so the binary
search reports a match with index `j`, and `a = i*m + j ≡ x (mod r)`. Since
`[a]P = [x]P = Q`, the attack-side equality check accepts `a`.

**Key uniqueness.** Suppose two records `j ≠ j'` share a key. Equal keys mean
`[j]P = [j']P` or `[j]P = -[j']P`. In the first case `[j - j']P = infinity`, so
`r | (j - j')`, impossible because `0 < |j - j'| < m < r`. In the second case
`[j + j']P = infinity`, so `r | (j + j')`, impossible because
`0 < j + j' < 2m = 2^129 < r`. Hence the table contains exactly one record per
key, the binary search is unambiguous, and the returned index is the true `j`.
The point at infinity is separated from every finite point by its tag, so the
`j = 0` case (`x = i*m`, `gamma_i = infinity`) is handled explicitly rather than
by a zero coordinate.

The algorithm therefore succeeds for every valid input: success probability `1`,
which is at least the required `0.39` over the uniform challenge distribution.
There is no randomness, no retry, no restart and no success amplification.

## Why the table is constructible

The table is a **sorted array**, not a hash table. There is no hash function, no
load factor, no probe sequence and no occupancy bound anywhere in the construction
or the lookup:

- Placement: a bottom-up merge sort of exactly `m` records, with exactly `L = 128`
  passes, two `m`-record buffers, at most `m` key comparisons and `2m` record moves
  per pass. Every record is placed at the position the sort assigns; nothing is
  rejected, and no pass can fail.
- Lookup: a binary search over a power-of-two array, of depth exactly `L = 128`
  plus one terminal comparison. The depth is a constant fixed by `m = 2^128`, not
  by the data and not by any hash quality.
- Capacity: the array has exactly `m` slots and holds exactly `m` records, so
  occupancy is exactly 1 by construction. There is no collision-resolution
  structure and no overflow chain.

Construction therefore has a fixed, data-independent cost, and the placement step
has no failure mode. The bound below charges the sort and the full 129-comparison
lookup on every giant step, whether or not the lookup is a hit.

## Instruction ledger

All counts are scalar AArch64 instruction equivalents and include loads, stores,
addressing and control. Field elements are 4 x 64-bit limbs; `p` has the fast
reduction form `p = 2^256 - 2^224 + 2^192 + 2^96 - 1`.

| Kernel | Bound | Derivation |
| --- | ---: | --- |
| `M` field multiply | 256 | 16 `mul` + 16 `mulh` = 32; 8-limb accumulation ~24; special-form reduction ~40; 8 limb loads/stores 16; call and loop 16; about 128, charged at 2x margin |
| `S` field square | 256 | charged as `M` |
| `A` field add/sub | 32 | 4-limb add/sub with carry, compare against `p`, one conditional correction |
| `I` field inversion | 86016 | addition chain for `p - 2` with 255 squarings and 13 multiplications = 268 field operations = 68608; charged 336 operations for control |
| mixed addition (Jacobian += affine, `add-2007-bl`) | 3584 | `7M + 4S = 11` field multiplications = 2816; 10 field additions = 320; loads, negation and control 448 |
| affine conversion, amortized | 2304 | batch-inverse share `86016/1024 = 84`; `6` multiplications for `X*Z^-2` and `Y*Z^-3` = 1536; 4 additions 128; control 556 |
| lookup (binary search) | 4384 | 129 key comparisons x 32 worst case + 256 for the record fetch and index arithmetic |

The mixed-addition formula is `Z1Z1 = Z1^2`, `U2 = x2*Z1Z1`, `S2 = y2*Z1*Z1Z1`,
`H = U2 - X1`, `HH = H^2`, `I = 4*HH`, `J = H*I`, `r = 2*(S2 - Y1)`, `V = X1*I`,
`X3 = r^2 - J - 2*V`, `Y3 = r*(V - X3) - 2*Y1*J`, `Z3 = (Z1+H)^2 - Z1Z1 - HH`,
which is 7 multiplications and 4 squarings. Because the addend is affine, `Z2 = 1`
and the formula is complete for distinct points; the infinity and equal-point
cases are dispatched before the formula (the addends `P` and `D` are fixed and
nonzero, and the accumulator is infinity only at `j = 0`, handled by the tag).

Sort comparison and record move, used only by the merge sort: `B = 32` per
lexicographic key comparison (tag word, up to four `x` limbs, parity word, branch)
and `D = 32` per record move (6 words loaded and stored), plus 64 per element for
buffer selection, index arithmetic and loop control.

Component costs:

| Component | Instruction equivalents |
| --- | ---: |
| setup (decoding, constants, `-[m]P` negation) | 16777216 |
| baby table: `m * (3584 + 2304 + 64)` | 2025360647913425734534005663433884394586112 |
| merge sort: `m*L*(32 + 2*32 + 64) + 4096*L` | 6968982874540819731729911960202612971143168 |
| giant search: `h * (3584 + 2304 + 4384)` | 3495380473011879896695783967539123068076032 |
| verification and output: `512*10272 + 2^24` | 22036480 |

`m = h = 2^128`, `L = 128`. The baby-table term includes one record emission per
baby step; the giant-search term charges a full mixed addition, a full amortized
affine conversion and a full 129-comparison lookup on every one of the `h` steps.
Verification and output reserve two full scalar multiplications (`2*bitlen = 512`
point operations) for the attack-side `[a]P = Q` check and the input validity
checks, plus 16 MiB for serialization and decoding.

Sum:

```
T = 16777216
  + 2025360647913425734534005663433884394586112
  + 6968982874540819731729911960202612971143168
  + 3495380473011879896695783967539123068076032
  + 22036480
  = 12489723995466125362959701591175620472619008
```

`T / 1000 = 1.24897...e40`, and `log2(T/1000) = 133.197865...`, submitted as
**133.197866** after rounding upward. The exact integer `T` is the audit
authority; the logarithm is a rounded presentation of it.

## Memory and preprocessing

Two `m`-record buffers of 48 bytes each: `2 * m * 48 = 32667107224410092492483962313449749348352`
bytes. The 1 MiB (`2^20`) reserve covers code, public constants, field and point
scratch, the 1024-element batch-inversion buffer, loop state and stacks; these are
bounded kernels over at most 12 coefficients and a 256-bit exponent. Peak memory
is `32667107224410092492483962313449749348352` bytes
(`memory_log2_bytes = 134.584963`).

Preprocessing is the whole baby table plus its sort, `8994343522454245466263917623636497365729280`
instruction equivalents (`preprocessing_log2 = 132.724214`), and it is counted
inside `T` with no reuse credit across targets. No table, buffer or record is
uncharged advice; every record is produced by the counted baby walk.

This memory is astronomical and mandatory for a generic baby-step/giant-step
construction. It is reported because the accounting requires it, not as a
practical claim.

## Heuristics

None. The algorithm is deterministic; the table is a sorted array with a
data-independent lookup depth; the coverage argument is exact; the cost bound is a
worst-case count over a fixed sequence of kernels. No smoothness assumption, no
random-oracle premise, no probabilistic occupancy and no success heuristic is
used. Success probability is exactly 1.

## Limitations and scope

- The correctness and cost arguments are analytic. No full-size execution is
  claimed or possible; the supplied experiment is a finite-group check of coverage,
  key uniqueness, the infinity tag and modular output recovery at toy size only.
- The exponential memory is required by this construction.
- This is a generic accounting bound, not a best-known-attack estimate, not a
  security claim and not a practical break. It does not weaken any hardness
  assumption.
- Beating this control alone does not establish novel cryptanalysis; automorphism-
  aware and rho-family methods must be compared under the same cost model before
  any best-known starting score is selected. This package does not implement them.

## References

- SEC 2 v2.0 section 2.4.2, `https://www.secg.org/sec2-v2.pdf#page=13` (curve
  parameters, reproduced in the protected profile).
- Guide to Elliptic Curve Cryptography section 3.6, `https://cacr.uwaterloo.ca/hac/about/chap3.pdf`
  (baby-step/giant-step and Pollard rho complexity).
- `add-2007-bl`, Explicit-Formulas Database, `https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html`
  (mixed Jacobian-plus-affine addition).

The construction and the accounting above are self-contained; review need not
fetch any citation to reconstruct them.

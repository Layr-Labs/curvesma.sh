# bls12-381-g1-ecdlp-plain: deterministic truncated baby-step/giant-step

Status: ready for independent review; no qualification, promotion or emitted
score is claimed by the author. This is an analytic accounting control with an
exactly computed success fraction, not a best-known-attack estimate and not a
practical break. The normative target and accounting rules are
`target-profile.json`, `tracks/classical-group-work-v2.json` and
`tracks/docs/COST_MODEL.md`. The construction is deterministic: it uses no
randomness, no heuristic premise, no random-oracle model, no uncharged advice
and no quantum step.

## Target, input and output

The input is `T = [x]P0` in BLS12-381 G1 for a uniformly sampled
`x` in `[0, r)`, where

- `p = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787` (381 bits),
- `r = 52435875175126190479447740508185965837690552500527637822603658699938581184513` (255 bits),
- `P0` is the profile generator of the prime-order subgroup G1 of order `r`.

No auxiliary powers, SRS secret or `n` parameter are supplied. The output is a
canonical 32-byte scalar `x` with `[x]P0 = T`, or an explicit failure. The
identity target `T = infinity` has scalar `0` and is handled by an explicit
infinity record.

## Fixed decomposition parameters

Set `N = ceil(41*r/100)` and

```
m = isqrt(N) + 1 = 146624380038933968876996171018221320040 = 0x6e4ed412a3c552fd2d010d8328cef768
h = ceil(N/m)     = 146624380038933968876996171018221320040 = 0x6e4ed412a3c552fd2d010d8328cef768
```

Then, with exact integer arithmetic,

```
m*h = 21498708821801738096573573608356245993507346963331534251623150184000105601600
100*(m*h) >= 41*r          (m*h >= ceil(41r/100) >= 41r/100)
m*h < r                    (the covered interval is a proper subset of [0,r))
(h-1)*m < r                (the giant walk never wraps past r in index space)
```

Both parameters are fixed constants of the curve and the track, not knobs
chosen per challenge.

## Construction

All group elements are affine G1 points encoded canonically as six 64-bit limbs
per Fp coordinate plus an infinity tag. Point additions use the batched
incomplete affine 6M convention of the reference operation: the per-batch
inversion is amortised exactly as the fixed reference price assumes, and no
separate inversion is charged. Every other instruction is charged below.

```text
1. Decode T; reject malformed encodings; check the curve equation and [r]T = O.
   If T = O, return the scalar 0.
2. Baby table. B_0 = O; for j = 1..m: B_j = B_{j-1} + P0.  (m additions)
   Store records (x(B_j), y(B_j), j) for j = 0..m-1 and retain D = -B_m = -[m]P0.
3. Sort the m baby records by x (ties by y) with LSD radix sort on the 381-bit
   x-coordinate using 24-bit digits (16 passes). The x-coordinate and the index
   are moved together; y-coordinates are held in an array indexed by j and are
   not moved.
4. Giant table. G_0 = T; for i = 1..h-1: G_i = G_{i-1} + D.  (h-1 additions)
   Store (x(G_i), y(G_i), i) for i = 0..h-1.
5. Sort the h giant records the same way.
6. Merge join the two sorted lists on x, breaking ties with y. For each equal
   pair (j, i), form x' = (i*m + j) mod r and verify [x']P0 = T with a
   double-and-add scalar multiplication (at most 512 additions, counted). If
   the check holds, return x' canonically. Otherwise continue the join.
7. If the join is exhausted with no verified pair, return failure.
```

Both sorts are deterministic: exactly 16 counting passes over every record,
with a 2^24-entry counter array (16.8 MB) that is reinitialised per pass. The
merge is a single linear pass over the two sorted lists. The failure path
performs exactly the same work as the success path, so the budget below is a
worst-case budget, not an expected-runtime figure.

## Correctness

Suppose `x < m*h`. Write `x = i*m + j` with `i = floor(x/m) < h` and
`j = x mod m < m`. Then

```
G_i = [x - i*m]P0 = [j]P0 = B_j,
```

so the merge join pairs the giant record `i` with the baby record `j`, forms
`x' = i*m + j = x`, and the verification `[x']P0 = T` holds. The attack
therefore recovers every `x` in the covered interval `[0, m*h)`.

Conversely, a merge pair with equal `x` and equal `y` is exactly `G_i = B_j`,
hence `x - i*m = j (mod r)`. Since `0 <= x < r`, `0 <= i*m < (h-1)*m < r` and
`0 <= j < m`, the congruence has at most one solution with `i < h`, `j < m`,
so at most one pair can verify, and it yields the true `x`. A pair with equal
`x` but opposite `y` is `G_i = -B_j`; it fails the verification and is skipped.
No incorrect scalar is ever returned: the returned value is always the one that
passed `[x']P0 = T`, and the equality check is the track's acceptance test.

## Success probability

The covered set is exactly the integer interval `[0, m*h)`; its size is `m*h`.
Since `x` is uniform on `[0, r)` and `m*h < r`,

```
Pr[success] = m*h / r
            = 21498708821801738096573573608356245993507346963331534251623150184000105601600
              / 52435875175126190479447740508185965837690552500527637822603658699938581184513
            >= 41/100 = 0.41 > 0.39.
```

This is an exact fraction of the uniform challenge space, not an estimate: the
algorithm is deterministic and has no coins, so the probability over fresh
algorithm coins is the same 0.41. Failures are reported explicitly and never
produce a wrong scalar; no amplification or restart is claimed or charged.

## Work ledger

Accounting convention: `H` reference operations (one batched incomplete affine
G1 addition each, priced at one target-group-operation equivalent) plus `R`
other 64-bit AArch64 instruction equivalents divided by the fixed track price
`C = 2000`. Every quantity below is charged exactly once; the two sorting
passes, the table construction, the giant walk, the merge, the verification and
the failure path are all included.

```
H = m + (h - 1) + 512 = 293248760077867937753992342036442640591
```

The `512` covers the worst-case double-and-add verification `[x']P0` for a
255-bit scalar (255 doublings plus at most 255 additions, rounded up), including
the point comparison; it is charged even when the attack fails.

Other instruction-equivalent work `R`:

| Component | Instruction equivalents |
| --- | ---: |
| LSD radix sort, both tables: 16 passes x 28 instructions x (m+h) | 131371464034489499500570804043375069710848 |
| record stores, both tables: 14 x (m+h) | 4105234641090271128555892788510196964880 |
| merge join: 12 x (m+h) | 3518772549505946681333622390151597398472 |
| decode, constants, counter reinitialisation, output | 33554432 |
| **R** | **138999912276909402495392370125273844952352** |

The 28 instructions per record per radix pass cover seven limb loads, seven
limb stores, digit extraction, counter update, destination-address computation
and loop control; the 14 covers writing a record once; the 12 covers a full
381-bit key comparison in the merge. The counter array is 2^24 eight-byte
entries; its per-pass reinitialisation and prefix sum are inside the 2^25
reserve, which also covers input decoding, the curve and subgroup checks, the
scalar serialisation and the final comparison.

Total instruction equivalents and target-group-operation equivalents:

```
T = H*C + R = 725497432432645278003377054198159126134352
W = T/C     = 45343589527040329875211065887384945383397 / 125
            = 3.627487162163227... x 10^38
time_log2   = log2(W) = 128.092238114259...
```

Rounded upward to six decimal places, the submitted score is **128.092239**.

Preprocessing (the whole baby table, its sort and the shared setup reserve) is

```
m + (m*(16*28 + 14) + 2^24)/C = 1.804946118279277... x 10^38
preprocessing_log2 = 127.085224 (rounded upward),
```

and is included in `T` with no reuse credit. There is no uncharged advice:
`nonuniform_advice_log2_bytes = 0`.

## Memory

The peak holds two 13-word record arrays (x, y, index) of `m` and `h` records,
one 7-word sort buffer and the 2^24-entry counter array:

```
peak = (m + h)*33*8 + 2^20 = 77417672660557135567053978297620858029696 bytes
memory_log2_bytes = 135.829787 (rounded upward).
```

This exponential memory is mandatory for the construction and is reported as a
supplementary bound; it has no scalar contribution or tie-break. The attack is
not practical at production size and is not claimed to be.

## Heuristics, uncertainty and limitations

There are no heuristic premises: the walk is a fixed decomposition, the tables
are complete for the covered interval, the sorts have deterministic pass
counts, and the success fraction is computed exactly. The only probability
space is the uniform challenge `x` in `[0, r)`, on which the success fraction
`m*h/r = 0.41` is exact. The bound is analytic; no production-size execution is
claimed or needed, and small experiments would be evidence only at their own
scope. The submitted bound is conservative: it charges a full 24-bit-digit radix
sort on both tables, a worst-case verification, the failure path and a 1 MiB
setup/output reserve. Organizer fixtures for the same decomposition check
coverage, the last partial block and modular output recovery; they do not
certify production cost. This package does not claim novelty, a weakened
hardness assumption, a practical break or improvement over Pollard-rho methods,
which would need their own supported success argument under this same model.

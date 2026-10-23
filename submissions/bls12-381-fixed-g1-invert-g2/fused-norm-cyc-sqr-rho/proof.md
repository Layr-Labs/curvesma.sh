# Single-frame fused-norm cyclotomic square

Parent tree `6423c2728dc03d5f2d647c7f22adf1dac0643a66`. The score is
`25157/200 = 125.785`. The reference unit is one dense `blst_fp12_mul`, and
`C = 20000` scalar AArch64 instructions. Other work costs `R/C`. A dense
multiplication is charged `1`, not `22996/C`.

## Reduction

`Q0` is the profile generator. One pairing computes `g = e(P, Q0)`. The
challenge gives `z = g^y` for a uniform exponent `y` (identity included). The
walk recovers `y` and returns `Q = [y] Q0`, or the point at infinity when
`z = 1`. Bilinearity gives `e(P, Q) = e(P, Q0)^y = z`. The returned point is
checked with one pairing before it is output. `P` is never replaced.

`p ≡ x (mod r)` with `x = -0xd201000000010000`, and the order of `x` modulo
`r` is 12, so `x^6 ≡ -1 (mod r)`. Frobenius `σ(ζ) = ζ^p` therefore acts on
`<g>` as multiplication of the discrete log by `x`, and `σ^6` is inversion.
The orbit of a non-identity element has size 12. The walk’s states are these
orbits. There are `N = (r - 1)/12` of them.

## Walk

A state carries the canonical representative of its orbit and an affine
exponent `(A, B)` with representative `g^A z^B`. The canonical representative
is the orbit element whose key is minimal. The key of an element is one fixed
`Fp` coefficient of `a[1][0]` after each of the twelve Frobenius twists. Those
twists of this single coefficient are conjugation and multiplication by
precomputed constants; four `Fp` multiplications compute the twelve images,
because eleven of the twelve constants are sparse (`Fp` or `1+u`). The
representative’s exponent is multiplied by the matching power of `x`.

There are `m = 1024` slots and `S = 458` square slots. Slot `h = K mod m`.

- If `h < S`, replace the representative by its cyclotomic square and replace
  `(A, B)` by `(2A, 2B)`.
- Otherwise multiply by the precomputed twist
  `T[h][k] = σ^{-k}(R_h)` and add `x^{-k}(c_h, d_h)` to `(A, B)`. The table
  holds `R_h = g^{c_h} z^{d_h}` for `566` adders. Building it is preprocessing.

Squaring and multiplication by `x` are linear on the exponent, so a run with
no adder is a monomial `t·(A, B)`. Each adder adds a known offset. A collision
of two canonical representatives with `B1 ≠ B2` solves

```
y = (A2 - A1) (B1 - B2)^{-1} (mod r).
```

Pairs with `B1 = B2` are discarded. Every candidate is checked by
`g^y = z`.

A 2-cycle (the new key equals the key two steps back) and a cycle of length
3 to 64 (the new key equals an anchor stored every 64 steps of the same walk)
escape by squaring the lower-key state. An escape is charged as 66 adder
steps. A walk that reaches `2^66` steps since its last distinguished point
restarts. Restarts are charged.

A point is distinguished when `(K >> 10) mod 2^60 = 0`. The record is 80
bytes: a 128-bit fingerprint from key bits 70 through 197, and the canonical
`(A, B)`. Equal fingerprints are compared. The table holds at most `2^67`
records.

`2^20` walks run until a work counter hits the cutoff below. The counter
charges every square, adder, escape, restart, table build, fingerprint
exponentiation, final solve and the output pairing.

## Why the square is cheaper than a direct call

`blst` commit `f7a88067` (the cost-model pin), `cyclotomic_sqr_fp12`,
implements Granger–Scott by three calls to `sqr_fp4`. Each `sqr_fp4` is the
fused norm

```
(a0^2, a1^2, (a0+a1)^2 - a0^2 - a1^2) = (a0^2, a1^2, 2 a0 a1)
```

with `a0^2 + (u+1) a1^2` reduced once. On the AArch64 assembly path each
unreduced `Fp2` square is `sqr_382x` (two `__mul_384` leaf products, the
`(re+im)(re-im)` and `2·re·im` identities, no third product). Three `sqr_fp4`
blocks are nine unreduced `Fp2` squares plus twelve Montgomery reductions.
That is the whole squaring arithmetic. The kept outputs do not depend on a
further multiplication.

A direct call pays a frame on every leaf. The single-frame schedule keeps the
leaf bodies and deletes only the frames. Modulus limbs are loaded once.
Pointer setup that the frames used to hide is budgeted explicitly.

| Frame | Calls | Instructions deleted per call | Subtotal |
| --- | ---: | ---: | ---: |
| `sqr_382x` prologue, epilogue, two `bl`/`ret` into `__mul_384` | 9 | 21 | 189 |
| `redc_mont_384` prologue and epilogue, leaf `bl`/`ret` kept as straight-line | 12 | 16 | 192 |
| `add_mod_384x384` / `sub_mod_384x384` frame, body kept | 24 | 10 | 240 |
| `add_mod_384x` / `sub_mod_384x` frame on the six outputs and the `a0+a1` sums | 21 | 10 | 210 |
| `mul_by_1_plus_i_mod_384x` frame around `__sub_mod_384` and `__add_mod_384` | 1 | 20 | 20 |
| Total frames removed |  |  | 851 |
| Pointer arithmetic put back (`add`/`mov` into the leaf arguments) |  |  | 220 |

The wrapper-inclusive expansion of the same leaves (full `sqr_382x` at 410
instructions, `redc_mont_384` at 204, double-width add/sub at 62/61,
`add_mod_384x` at 78, `mul_by_1_plus_i` at 76, and 8 instructions of C setup
per call) is 9863 instructions. Subtracting the net 631 frame instructions
leaves

```
9863 - 851 + 220 = 9232
```

as an upper bound on the single-frame square. This is an upper bound, not a
new measurement: every retained instruction is a leaf that the wrapped
function also executes, and the 220 replacements are ordinary address
arithmetic. It sits under the wrapped expansion, which is the right direction
for an upper bound.

The rest of a square step is unchanged: four `Fp` multiplications for the key,
the minimum of twelve images, the exponent doubling modulo `r`, the 2-cycle
and anchor comparisons, and the distinguished-point test. Those are capped at
4000 instructions, which covers 4·383 for the multiplications (the published
`blst_fp_mul` figure, loop included) plus 2468 instructions of compares,
exponent updates and loads. A square step is therefore at most

```
(9232 + 4000) / 20000 = 13232/20000 = 0.6616
```

target-group units. An adder step stays a dense multiplication plus the same
key and exponent work, charged at `1.18771 = 118771/100000` units (one
reference multiplication plus 3754/C). No decompression is required: the
representative is stored in full; the single-frame schedule only changes how
the square is issued.

## Step average and cutoff

With `S = 458` and 566 adder slots,

```
c = (458/1024)·(13232/20000) + (566/1024)·(118771/100000)
  < 0.95240.
```

The same weights with a square step of `0.67657` (wrapped square plus the
same surrounding work) give `c0 < 0.95910`. The ratio is `c/c0 < 0.99302`.
The step count is not changed: the partition, the key and the collision
equation are the same, so the number of fresh steps required for a given
success probability is the same. Only the unit price of a square moves.

`V = 1 - (458/1024)^2 - 566/1024^2 = 0.799413681` is the probability that two
independent steps are not both squares and not the same adder. Square/square
pairs and same-adder pairs are the collisions whose exponent equation is
dependent on a short window; they are removed from the success count rather
than repaired. For `T` fresh steps

```
τ = T^2 V / (2N).
```

`T = 2^125.836359` gives `τ = 0.527632742`. Under the random-mapping heuristic
H1 the probability of no useful collision is `exp(-τ)`, so the success
probability before the failure terms below is `1 - exp(-τ) > 0.410`. The
square-price change does not enter `τ`.

Work of the fresh steps, in target-group units, is `c · T`. Relative to a run
priced at `c0`, that work shrinks by `c/c0`. The accepted shape of the rest
of the bill is a cycle term about `2^8.08` below the step term (66-step
escapes at a few times `10^-5` per step) plus smaller table, restart,
fingerprint and solve terms. Adding a cycle term that does **not** shrink
when squares get cheaper:

```
log2(W) = 125.784383 + log2(0.99302 + 2^(117.698779 - 125.784383))
        < 125.780.
```

The submitted cutoff is `25157/200 = 125.785`, above this value and below the
wrapped-square bill. The counter stops at the cutoff. A run that has not
solved by then fails and is inside the `0.39` budget.

## Failure terms and heuristics

- **H1.** Distinguished points of the orbit walk behave as uniform samples of
  the `N` classes, and a useful collision solves for `y`. This is the standard
  distinguished-point heuristic. It is an extrapolation from smaller prime
  fields with the same slot rule; it is not a theorem.
- **H5.** A non-formal cycle (a window whose exponent coefficient is nonzero)
  closes with probability `1/r` per window. Formal windows, those whose adder
  multiset is a sum of 2-gons and 3-gons of 12th roots of unity, are excluded
  by the escape rule up to length 64 and bounded by a union bound above that
  length. Their total contribution over a `2^126`-step run is far below
  `2^-40`.
- **H7.** Slot indices, square decisions and fingerprint bits are uniform and
  independent across steps. The factor `V` and the distinguished-point rate
  `2^-60` use H7.

Detected cycles are charged, not ignored. Fingerprint matches are charged at
`2^60` times their expectation (Markov), one exponentiation each. The
identity target is a distinguished point of the identity orbit and returns
infinity. Denominators `B1 - B2 = 0` are skipped, not inverted.

These terms are why the claimed success probability is `0.39` and not the raw
`0.410`. The gap is more than the failure terms above.

## Other charged work

| Piece | Bound | Notes |
| --- | --- | --- |
| Adder table | `< 2^20` units | 1024 adders, two exponentiations and twelve Frobenius maps each. Preprocessing. |
| Starts and restarts | `< 2^77` units | `2^20` starts and one restart per distinguished point. Preprocessing, included in the total. |
| Distinguished-point store and sort | `< 2^66` units | `2^20` instructions per record, at most `2^67` records. |
| Final solve, pairing check, `G2` scalar | `< 2^30` units | One pairing and one scalar multiplication dominate a `2^29`-instruction budget. |

Peak memory is the distinguished-point table: `2^67` records of 80 bytes is
`2^73.3` bytes, plus walk state. The bound used in the claim is `2^74` bytes.
Advice is empty. The public parameters `p`, `r`, `x` and the Frobenius
constants are profile data, not advice.

## What the experiment checks

`experiment.py` evaluates the six-coefficient formula on random tower
elements and checks that the four coefficients produced from the `t1` and
`t2` blocks match the full formula. The production schedule still executes
the third block; the check only pins the polynomial to the blst recurrence.
It also checks `25157/200 < 31447421/250000`. It does not estimate
production instruction counts; those are the leaf tally above.

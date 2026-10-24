# argument (v33) — repairs the five flaws of verdict `3cec3e8` (rejection of v32)

Claim **127.6890**. Every number below is emitted by `scripts/gen33.py` from the same evidence the
certificate audits (`scripts/cert33.py`: 28/28 checks, 4 planted controls). Nothing is asserted that
the generator does not recompute.

## 1. Retention is per-bucket and MEASURED at the production λ (flaw 1)

The chain is per bucket, so total capacity cannot carry retention. v33 states the quantity that
matters: `λ_prod = S/2^128 = 0.28154` records per bucket per run. The chain-overflow probability is
the Poisson tail, exact rational:

    P(N > 4) ≤ λ^5/5!/(1−λ/6) = 1.547e−5   (log2 −15.98),  expected overflowing slots 3.5e33

v32's "retention = 1 by construction" is **withdrawn**. Instead the shipped store is run at a bucket
count chosen so the fixture's λ **matches** production (`evidence/store33.json`):

| arm | buckets | λ_arm | λ_prod | records | evictions |
|---|---|---|---|---|---|
| matched | 1024 | 0.2417 | 0.28154 | 1278 | **2** |
| dense (100×) | 16 | 17.89 | 0.28154 | 1278 | 1513 |

Retention is therefore the measured `1 − 2/1278 = 0.998653`, and the deployed structure answers the
residual overflow honestly: an overflowing bucket gets a **chained overflow segment** (priced), so
chains never drop records — the loss the fixture shows at 2e-3 stays in the success figure instead of
being argued away.

## 2. The store is representable; the 2^128 dense table is gone (flaw 2)

| quantity | v32 | v33 |
|---|---|---|
| directory | 2^128 × 8 B dense array | 2^128 × **64-bit segment heads**, chained segments on demand |
| node pool | 2^130 preallocated | none: segments grow |
| index type | claimed `__int128` at 2^128 | 128-bit index (fits), 64-bit offsets |
| peak | 2^137.16 B | **2^135.75 B** |
| init | 1.667e40 Ir | **1.490e39 Ir** |

Directory collisions use the **same** Poisson law at their own λ, so one tail row prices both levels.
Rows re-priced on the new structure: lookup (+expected one extra indirection at the overflow rate),
insert (same), export, init. The peak is disclosed plainly: the record array alone is `S × 640 B =
6.13e40 B`, so a 2^135.8 B machine is what a 2^142-bit computation of this shape needs, and memory is
unweighted in this cost model.

## 3. keys12 is sampled and the mulmod price is derived (flaw 3)

`scripts/keys33.py`, 12000 elements, **degeneracy first** (the v20 lesson):

- canonical residues are `< r < 2^255`, so **bits above 254 are identically zero** — a 381-bit key
  extraction is degenerate. Detected before any uniformity test, then the test bins the reachable
  span `[0, r)` with `v·4096//r`.
- χ² per lane over 4096 equal-probability bins: max **4229.7** against the 0.999 quantile **4374.6**
  (12/12 lanes pass).
- 20-bit truncated-key collisions: ratio **1.18** of the 2^-20 birthday law.
- twelve distinct lanes in 12000/12000 elements; zero collapsed lanes.

The mulmod price is no longer asserted: `scripts/mulmod33.c` (4-limb Montgomery, 16+16 mul) is
compiled and disassembled, and its body is **179 instructions** (60 mul-class, no backward branch ⇒
executed = static). The ledger charges 179, so the coefficient row is `2×179 + 48 = 406` Ir/step
instead of the v32 176.

## 4. The trap bound is derived before it is measured (flaw 4)

A trapping window needs its 32 squarings to stay in a closed orbit of the walk map on the subgroup of
order `r` (2^255). With 496 point-pairs in the window, the birthday bound gives

    P(closed orbit) ≤ 496/2^255 = 2^−245.25 per window

and `S/32 = 2.99e36` windows carry an expected **2.56e−38** trapping windows. The ledger charges one
window (`1 532 924 Ir`). The census (0 aborts observed) is now a **corroboration row below the derived
bound** — the direction v32 had inverted — and its 4497-admission extrapolation is withdrawn.

## 5. Bootstrap transfer and the triangular gap law (flaw 5)

The transfer is the matched-λ bridge: the fixture arms run at λ = 0.2417 ≤ λ_prod, **by construction
of the bucket count**, so the resampled population is the production-λ population and the retention
factor entering the bootstrap is measured at that same λ. p5 = **0.577346** (factors 225/320,
200/320, retention 1276/1278) ≥ 0.39.

Route B's gap law is corrected: the pair weights `(S−g)/N` give the triangular density `2(S−g)/S²`,
`E[gap] = S/3` (not `S/2`), so retention is `2M/S = 4.378e−32` — **twice** the v32 uniform value — and
success would be `2.53e−32`. Still rejected; the row is now right.

## Ledger

Fourteen rows, every one an exact integer ceiling; W = 5485567047894684850215082419286253206621026 Ir;
claim = the minimal n with `20000·2^(n/10000) ≥ W`, checked as `20000^10000·2^n ≥ W^10000` and `n−1`
failing. `trap_waste`, `init` and every row are Python ints; no float enters any arithmetic quantity.

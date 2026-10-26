# Kernel measurements behind the cost constants

The track metric counts executed scalar AArch64 instructions. Every constant of the relation-collection
model is a **retired-instruction count** of a correctness-checked kernel in this package, measured with
the OS counter (`/usr/bin/time -l`, "instructions retired") on Apple-silicon AArch64, `clang -O2`, no
fast-math. Each figure is a difference between two runs of the same binary with different repetition
counts, divided by the extra work, so fixed setup (allocation, lattice reduction, table building) cancels;
three trials are taken and the **maximum** is used. The verifier's Python sandbox does not compile C, so the
kernels are evidence to read and re-run on an AArch64 host (`bash measure_all.sh`).

## `slice6.c`: the per-a slice sieve (ARCH-E)

**What it is.** One a-slice of the sieve: the 6-dimensional half-ellipsoid
`{b in Z^6 : b^T M_b b <= R^2, b5 >= 1}` with a byte bucket per point (fixed-cube addressing over the
bounding box: `addr = sum_i (b_i - lo_i) S_i`, injective by construction). For a sigma-pair of prime
ideals above the prime ideal `(p, y - s)` of `O_h`, all of them share the lattice
`L(p,s) = {b : sum_i b_i s^i = 0 mod p}` (index p in Z^6); the ideal and the slice only fix the coset
`c e_0 + L(p,s)` with `c = a(s) r^{-1} mod p`. The kernel enumerates the coset inside the **full**
ellipsoid with a Schnorr-Euchner walk over an LLL-reduced basis (reduced in the metric `M_b`), carrying
the bucket address and `b5` as affine state, and **folds** `b -> -b` when `b5 < 0`
(`addr(-b) = K - addr(b)`), which supplies the conjugate ideal's updates (the `x -> -x` automorphism maps
the relation of `(a,b)` to that of `(a,-b)`); walking one ideal per sigma-pair over the full ellipsoid
therefore costs the same as walking both over the half.

**Correctness.** `slice6 check log2pts p s c` runs the walk and a brute-force reference (every `b` of the
bounding box with `b(s) = c mod p`, tested against the quadratic form directly) and compares the bucket
arrays by digest. `bash checks.sh` runs six lattices from p = 251 to p = 2^26 (33,468 down to 1 point per
coset): all bit-identical (`MATCH`).

**Per-update cost** (`slice6 icount log2pts p s c reps`, slice of 2^26 points, reps 16 vs 48) and
**per-coset cost** in the sparse regime (`slice6 icountc log2pts p s ncos`, 20,000 vs 60,000 cosets), as a
function of the expected number of points per coset `h = |slice| / p` (the quantity that fixes the walk's
regime; the production slice has ~2^63 points, so a given `h` corresponds to `p = 2^63 / h`):

| regime (h = points per coset) | measured (3 trials) | charged in the model |
|---|---|---|
| h ~ 5.2e5 (long runs: the innermost loop is a stride loop with the fold) | 16.95-16.98 / update | 17.0 |
| h ~ 6.5e4 | 21.9-23.7 / update | 23.7 |
| h ~ 1.6e4 | 29.4-29.9 / update | 29.9 |
| h ~ 2.1e3 | 36.7-38.0 / update | 38.0 |
| h ~ 258 | 65.7-67.3 / update | 67.3 |
| h ~ 16 | 2134-2256 / coset (134-141 / update) | 141 |
| h ~ 4 | 1116-1170 / coset (280-293 / update) | 293 |
| h ~ 1 | 672-748 / coset | 748 |
| h ~ 0.25 | 306-425 / coset | 750 per coset (h < 1 regime) |
| h ~ 0.07 | 347-401 / coset | 750 per coset |

The model interpolates linearly in log2 h between the table points (`UPD_TAB` in `evidence.py`) and charges
750 instructions **per coset** (i.e. per (slice, sigma-pair)) whenever h < 1, where the walk is an almost
always empty exact enumeration. No Babai/rounding shortcut is credited: every coset is enumerated exactly,
so no update is ever missed. (`slice6 babai` measures the rounding shortcut's miss rate, 14-38%; it is not
used.)

**Per-(slice, pair) setup** (`slice6 setup p s trials`, 1e6 vs 3e6 trials): Horner evaluation of `a(s) mod p`
with a precomputed Barrett constant (6 steps), `c = a(s) r^{-1} mod p`, the fractional coordinates of the coset
origin in the reduced basis and the Gram-form test: **254 instructions**; charged as 254 + 72 (the 36 words
of per-pair data - p, s, r^{-1}, reduced basis, Gram-Schmidt data, address deltas - streamed to the slice:
one load and one store per word) = **326 per (slice, sigma-pair)** for every pair with p <= X2.

## `la_mac.c`, `la_mac.s`, `la_finish.c`, `la_test.c`: block-Wiedemann entry (C_ENTRY = 35)

Unchanged from v1 (the reviewer accepted the 33-instruction count). Rows are stored sign-split (positive
entries, then negative, magnitudes as `uint8`); the branch-free loop `LBB0_2` retires exactly 33
instructions per entry (`evidence.py` part [7] re-counts it from the shipped assembly); the per-row finish
(two loop prologues/epilogues and the Barrett reduction `row_finish`, ~100 instructions) is amortised over
the row weight (>= 100): **C_ENTRY = 35**. `la_test.c`: 2,000 random rows (weights 1-400, coefficients up to
255, vectors in [0, r) including r - 1) match Python big integers. The filtering simulation enforces
|coefficient| <= 255 at merge time and reports the maximum it produced (<= 60 at every scale).

## `canon.c`: units and canonical a (not a cost constant)

Brute-force enumeration of the units of `Z[y]/(h)` of height <= H (norm +-1 by the product of the complex
embeddings), and the canonical-a test used to define the fundamental domain F. It produced the brute-force
count 595 (height 4, mod sign) that `evidence.py` part [U] reproduces from the three fundamental units, and the
flags of the pools; heights 3, 4 and 6 give identical flags on the whole pool.

## Reproduction

    clang -O2 -o slice6 slice6.c -lm && bash checks.sh && bash measure_all.sh
    clang -O2 -o canon canon.c && ./canon -1 -1 1 1 0 0 4 1024 <matrix file> < <points>
    clang -O2 -o la_test la_test.c la_mac.c la_finish.c

(`slice_mb.h` is the b-block of the shipped ellipsoid; `measure_all.sh` prints the three trials of every
figure above.)

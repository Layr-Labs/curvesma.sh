# Kernel measurements behind the cost constants

The track metric counts executed scalar AArch64 instructions. Every constant of the relation-collection
model is a **retired-instruction count** of a correctness-checked kernel in this package, measured with
the OS counter (`/usr/bin/time -l`, "instructions retired") on Apple-silicon AArch64, `clang -O2`, no
fast-math. Each figure is the difference between two runs of the same binary with 3x and 1x the repeated
work (the repeated work is >= ~1e8 instructions so fixed setup and run-to-run noise are negligible), divided
by the extra work; three trials are taken and the **maximum** is charged (`bash measure7.sh`). The verifier's
Python sandbox does not compile C, so the kernels are evidence to read and re-run on an AArch64 host.

## `slice7.c`: the per-a slice sieve, production-width (ARCH-E)

**Layout.** One a-slice is the 6-dimensional half-ellipsoid `{b in Z^6 : b^T M_b b <= R^2, b5 >= 1}`,
organised in *lines* (b5 = k) and *segments* (ranges of lines). Every line uses the same 5-dimensional
bounding box (the projection of the full ellipsoid), so a point has the address
`pos = (k-1) L5 + sum_{i<5} (b_i + r_i) S_i`, injective by construction; a segment of G lines owns the byte
array `pos in [(k0-1) L5, (k0-1+G) L5)` of at most 2^62 bytes, so **every address is a 64-bit integer**
within one segment, and the segment of a point is `(k-1)/G`. The conjugate ideal's fold `b -> -b` keeps the
line and maps the local part to `K5 - local`. All coordinates, addresses and address deltas are 64-bit
integers; the benchmark's values are small, the arithmetic is the production one.

**Recentred cosets.** For a sigma-pair above the O_h-prime `(p, y - s)` the ideals share the lattice
`L(p,s) = {b : b(s) = 0 mod p}` (6-dim) or `L5(p,s) = {b' : sum_{i<5} b_i s^i = 0 mod p}` (5-dim, per line);
the ideal and the slice fix the coset origin `(c, 0, ..., 0)` with `c = a(s) r^{-1} mod p` (`c_k = c -+ k s^5`
on line k). The origin is reduced **exactly** (integer Babai rounding in the reduced basis:
`o' = o - sum round(t_i) B_i`) before any floating-point use, so the Schnorr-Euchner centre coordinates are
O(1) for every p; the walk's integer coordinates never exceed a few units. The reduced bases are computed once
per (p, s) at factor-base construction (charged there), not per slice.

**Two paths.**
- *Line path* (p < 2^40): for each line and each of the two cosets (the ideal at +k, the ideal at -k folded
  onto +k = the conjugate's coset at +k), the 5-dimensional coset is enumerated inside the line's ellipsoid
  (centre `-k M5^{-1} m`, radius^2 `R^2 - k^2 d`, from completing the square in Q_b) by a Schnorr-Euchner walk
  carrying the 64-bit local address; bucket bytes are updated directly.
- *Record path* (p >= 2^40): one 6-dimensional coset enumeration per pair over the full ellipsoid, carrying
  the line index and the local address as affine state; a hit with b5 < 0 is folded; each hit appends a
  12-byte record `{pos, prime index}` to its segment's buffer; an apply pass per segment adds the prime
  batch's log byte to `bucket[pos]` (survivors read the records to recover their large sieved primes).

**Correctness.** `slice7 check log2pts K p s c` runs both paths (the line path over all lines, the record
path plus the apply pass) and a brute-force reference over the whole half-slice (every point of the box with
`b(s) = c` or `= -c` mod p for the two ideals, tested against Q_b directly) and compares the complete bucket
arrays by digest. `bash checks7.sh` runs seven lattices from p = 251 to p = 2^26 (7,394 down to 1 hit): all
three digests identical (`MATCH`) in every case.

**Measured costs** (`measure7.sh`, slice of 2^28 points; h = expected points per coset: per *line* coset for the
line path, per *slice* coset for the record path):

| path | h (points per coset) | measured, max of 3 trials | charged |
|---|---|---|---|
| line | 2^15.5 (p = 257) | 14.7-15.6 per update | 15.6 |
| line | 2^11.5 (p = 4099) | 19.4-23.7 per update | 23.7 |
| line | 2^7.5 (p = 65537) | 56.0-56.6 per update | 56.6 |
| line | 2^3.5 (p = 2^20) | 188-191 per update | 191 |
| line | 2^0.6 (p = 2^23) | ~1450 per coset | 1450 per coset (not reached in production: p < 2^40 gives h >= 2^15) |
| record | 2^21 (p = 257) | 32.6 walk+record, 37.0 with apply | 37.0 |
| record | 2^13 (p = 65537) | 52.3 / 58.0 | 58.0 |
| record | 2^6 (p = 2^23) | 152 / 167 | 167 |
| record | 2^3 (p = 2^26) | 509 / 576 | 576 |
| record | < 1 (p = 2^29, 2^31) | 1696-1745 walk+record, 2088-2128 with the per-rep segment loop | 2130 per coset |


The model interpolates linearly in log2 h between the table points (`UPD_TAB` for the line path,
`REC_TAB` for the record path, in `evidence.py`) and charges 2130 instructions **per coset** (per
(slice, pair)) on the record path whenever h < 1, where the walk is an almost always empty exact enumeration.
No rounding shortcut is credited: every coset is enumerated exactly, so no update is ever missed.

**Setups.** Line path, per (pair, line, coset): `c_k` (Barrett), exact recentring, walk centre:
**866 instructions** (`slice7 setup`, 1e6 vs 3e6 trials). Record path, per (slice, pair): charged as the
line-path setup figure 866 (same recentring and centre computation, six-dimensional) + <= 75 for `a(s) mod p`
and `c` (six Barrett steps) + 72 for the 36 streamed words of per-pair data (one load and one store each) =
**1013**.

## `la_mac.c`, `la_mac.s`, `la_finish.c`, `la_test.c`: block-Wiedemann entry (C_ENTRY = 35)

Unchanged from v1 (the reviewer accepted the 33-instruction count). Rows are stored sign-split (positive
entries, then negative, magnitudes as `uint8`); the branch-free loop `LBB0_2` retires exactly 33
instructions per entry (`evidence.py` part [7] re-counts it from the shipped assembly); the per-row finish
(two loop prologues/epilogues and the Barrett reduction `row_finish`, ~100 instructions) is amortised over
the row weight (100): **C_ENTRY = 35**. `la_test.c`: 2,000 random rows (weights 1-400, coefficients up to
255, vectors in [0, r) including r - 1) match Python big integers. The filtering simulation (dense columns
propagated through the merge) reports the maximum coefficient it produced: <= 60 at every scale.

## `canon.c`: units and canonical a (not a cost constant)

Brute-force enumeration of the units of `Z[y]/(h)` of height <= H (norm +-1 by the product of the complex
embeddings), and the canonical-a test used to define the fundamental domain F. It produced the brute-force
count 595 (height 4, mod sign) that `evidence.py` part [U] reproduces from the three fundamental units, and the
flags of the pools; heights 3, 4 and 6 give identical flags on the whole pool.

## Reproduction

    clang -O2 -o slice7 slice7.c -lm && bash checks7.sh && bash measure7.sh
    clang -O2 -o canon canon.c && ./canon -1 -1 1 1 0 0 4 1024 <matrix file> < <points>
    clang -O2 -o la_test la_test.c la_mac.c la_finish.c

(`slice_mb.h` is the b-block of the shipped ellipsoid; `measure7.sh` prints the three trials of every figure
above.)

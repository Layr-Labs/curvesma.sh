# Kernel measurements behind the cost constants

The track metric counts executed scalar AArch64 instructions. Every constant of the relation-collection
model is a **retired-instruction count** of a correctness-checked kernel in this package, measured with
the OS counter (`/usr/bin/time -l`, "instructions retired") on Apple-silicon AArch64, `clang -O2`, no
fast-math. Each figure is the difference between two runs of the same binary with 3x and 1x the repeated
work (the repeated work is >= ~1e8 instructions so fixed setup and run-to-run noise are negligible), divided
by the extra work; three trials are taken and the **maximum** is charged (`bash measure8.sh`). The verifier's
Python sandbox does not compile C, so the kernels are evidence to read and re-run on an AArch64 host.

## `slice8.c`: the per-a slice sieve, production-width (ARCH-E)

**Layout.** One a-slice is the 6-dimensional half-ellipsoid `{b in Z^6 : b^T M_b b <= R^2, b5 >= 1}`,
organised in *lines* (b5 = k) and *segments* (ranges of lines). Every line uses the same 5-dimensional
bounding box (the projection of the full ellipsoid), so a point has the address
`pos = (k-1) L5 + sum_{i<5} (b_i + r_i) S_i`, injective by construction; a segment of G lines owns the byte
array `pos in [(k0-1) L5, (k0-1+G) L5)` of at most 2^62 bytes, so **every address is a 64-bit integer**
within one segment, and the segment of a point is `(k-1)/G`. The conjugate ideal's fold `b -> -b` keeps the
line and maps the local part to `K5 - local`. All coordinates, addresses and address deltas are 64-bit
integers (primes < 2^63); the benchmark's values are small, the arithmetic is the production one.

**Recentred cosets.** For a sigma-pair above the O_h-prime `(p, y - s)` the ideals share the lattice
`L(p,s) = {b : b(s) = 0 mod p}` (6-dim) or `L5(p,s) = {b' : sum_{i<5} b_i s^i = 0 mod p}` (5-dim, per line);
the ideal and the slice fix the coset origin `(c, 0, ..., 0)` with `c = a(s) r^{-1} mod p` (`c_k = c -+ k s^5`
on line k). The origin is reduced **exactly** (integer Babai rounding in the reduced basis:
`o' = o - sum round(t_i) B_i`) before any floating-point use, so the Schnorr-Euchner centre coordinates are
O(1) for every p; the walk's integer coordinates never exceed a few units. The reduced bases are computed once
per (p, s) at factor-base construction (charged there), not per slice.

**Two paths.**
- *Line path* (2^8 <= p < 2^20): for each line and each of the two cosets (the ideal at +k, the ideal at -k folded
  onto +k = the conjugate's coset at +k), the 5-dimensional coset is enumerated inside the line's ellipsoid
  (centre `-k M5^{-1} m`, radius^2 `R^2 - k^2 d`, from completing the square in Q_b) by a Schnorr-Euchner walk
  carrying the 64-bit local address; bucket bytes are updated directly.
- *Record path* (2^20 <= p <= X2): one 6-dimensional coset enumeration per pair over the full ellipsoid, carrying
  the line index and the local address as affine state; a hit with b5 < 0 is folded; each hit appends a
  12-byte record to its segment's buffers, as a structure of two arrays: the **segment-local** position
  `((k-1) mod G) L5 + local` (uint64) and the prime index (uint32); an apply pass per segment adds the prime
  batch's log byte to the segment's own bucket array at the local position (survivors read the records to
  recover their sieved primes >= 2^20). Primes are `uint64_t`.

**Arithmetic.** The modular setups (c = a(s) r^{-1} mod p by Horner, c_k = c -+ k s^5) use a Barrett reduction
valid for every p < 2^63 (k = bitlength(p), mu = floor(2^(2k)/p), q = ((x >> (k-1)) mu) >> (k+1), two
conditional subtractions); `slice8 barrett` validates it against the C remainder on 122,000 random reductions
with p of every bit length below 63 (0 mismatches), and the line path uses it for its coset constants, so the
`check` digests exercise it. Coset origins are recentred by integer Babai rounding; the rounding coefficients
come from a floating-point inverse basis, which can only change a coefficient by one, i.e. move the origin by one
reduced basis vector: any integer rounding gives a valid origin of the same coset, and the walk's centre is then
computed from the exact integer origin.

**Correctness.** `slice8 check log2pts K p s c` runs both paths (the line path over all lines, the record
path plus the apply pass) and a brute-force reference over the whole half-slice (every point of the box with
`b(s) = c` or `= -c` mod p for the two ideals, tested against Q_b directly) and compares the complete bucket
arrays by digest. `bash checks8.sh` runs eight lattices from p = 251 to p ~ 2^62 (7,394 down to 0 hits): all
three digests identical (`MATCH`) in every case, and so with the segment count set to 7 and to 300 (`NSEG=7
./slice8 check ...`; the number of segments is a runtime parameter, 683 in production); a prime >= 2^63 makes
the kernel abort by design (the claim's X2 = 2^63 keeps every sieved prime below that bound). The reference
enumerator writes into the per-segment arrays with the same segment-local addressing; the segments are
separate allocations; record counters are 64-bit; no global product K L5 appears anywhere in the kernel
(the benchmark's record buffers are rings that wrap, where production flushes a full buffer: a flush is a
sequential write of the buffer, already counted as the record's store).

**Measured costs** (`measure8.sh`, slice of 2^28 points; h = expected points per coset: per *line* coset for the
line path, per *slice* coset for the record path):

| path | h (points per coset) | measured, max of 3 trials over three passes | charged |
|---|---|---|---|
| line | 2^15.5 (p = 257) | 6.6-17.6 per update | 17.6 |
| line | 2^11.5 (p = 4099) | 17.3-25.5 per update | 25.5 |
| line | 2^7.5 (p = 65537) | 51.3-74.5 per update | 74.5 |
| line | 2^3.5 (p = 2^20) | 177-191 per update | 191 |
| line | 2^0.6 (p = 2^23) | ~1450 per coset | 1450 per coset (not reached in production: p < 2^20 gives h >= 2^35) |
| record | 2^21 (p = 257) | 32.6-41.8 walk+record, 36.1-39.5 with apply | 41.8 |
| record | 2^13 (p = 65537) | 41.7-69.8 / 55.9-58.8 | 60.0 |
| record | 2^6 (p = 2^23) | 149-152 / 162-167 | 166.7 |
| record | 2^3 (p = 2^26) | 508-513 / 565-576 | 576 |
| record | < 1 (p = 2^29, 2^31) | 1693-1748 walk+record, 2078-2128 with the per-rep segment loop | 2130 per coset |

(The spread between trials of the dense line regime comes from the small amount of work per run at these
densities, not from the kernel; the maximum is charged.)

The model interpolates linearly in log2 h between the table points (`UPD_TAB` for the line path,
`REC_TAB` for the record path, in `evidence.py`) and charges 2130 instructions **per coset** (per
(slice, pair)) on the record path whenever h < 1, where the walk is an almost always empty exact enumeration.
No rounding shortcut is credited: every coset is enumerated exactly, so no update is ever missed.

**Setups.** Line path, per (pair, line, coset): `c_k` (Barrett), exact recentring, walk centre:
**869 instructions** (`slice8 setup`, 1e6 vs 3e6 trials). Record path, per (slice, pair): `a(s) mod p`
(Horner, six Barrett steps), `c = a(s) r^-1`, exact recentring, walk centre: **1093 instructions**
(`slice8 setup6`) + 72 for the 36 streamed words of per-pair data (one load and one store each).

## `la_mac.c`, `la_mac.s`, `la_finish.c`, `la_test.c`: block-Wiedemann entry (C_ENTRY = 36)

The v1 kernel with **64-bit column indices** (`const uint64_t *idx`, loaded with `ldr x13`; the reviewer
accepted the 33-instruction count of the 32-bit version and the loop is otherwise identical). Rows are
stored sign-split (positive entries, then negative, magnitudes as `uint8`); the branch-free loop `LBB0_2`
retires exactly 33 instructions per entry (`evidence.py` part [7] re-counts it from the shipped assembly); the
per-row finish (two loop prologues/epilogues and the Barrett reduction `row_finish`, ~100 instructions) and
the combination of the four column-block accumulators (proof.md section 5) are amortised over the row weight
(100): **C_ENTRY = 36**. `la_test.c`: 2,000 random rows (weights 1-400, coefficients up to
255, vectors in [0, r) including r - 1) match Python big integers. The filtering simulation (dense columns
propagated through the merge) reports the maximum coefficient it produced: <= 60 at every scale.

## `canon.c`: units and canonical a (not a cost constant)

Brute-force enumeration of the units of `Z[y]/(h)` of height <= H (norm +-1 by the product of the complex
embeddings), and the canonical-a test used to define the fundamental domain F. It produced the brute-force
count 595 (height 4, mod sign) that `evidence.py` part [U] reproduces from the three fundamental units, and the
flags of the pools; heights 3, 4 and 6 give identical flags on the whole pool.

## Reproduction

    clang -O2 -o slice8 slice8.c -lm && bash checks8.sh && bash measure8.sh
    clang -O2 -o canon canon.c && ./canon -1 -1 1 1 0 0 4 1024 <matrix file> < <points>
    clang -O2 -o la_test la_test.c la_mac.c la_finish.c

(`slice_mb.h` is the b-block of the shipped ellipsoid; `measure8.sh` prints the three trials of every figure
above.)

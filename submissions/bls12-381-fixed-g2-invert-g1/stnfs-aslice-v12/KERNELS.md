# Kernels: correctness checks and executed-instruction measurements

Every constant of the relation-collection model (`model12.py`, top) is a **callgrind `Ir` count** (executed
instructions) of a kernel in this package on linux/arm64 (Debian 12, GCC 12.2.0, valgrind 3.19.0, the track's flags
`-O3 -fno-tree-vectorize -fno-tree-loop-distribute-patterns`; scalar binary64 instructions count one each), by the
track's **differential method**: the whole program at two repetition counts, difference divided by the extra units.
`measure12.sh` is the script, `measurements12.txt` the raw log. The Python sandbox does not compile C, so the kernels
are evidence to read and re-run on an AArch64 host.

## `slice12.c` (+ `slice12_check.c`): the per-a slice sieve with packed rows (ARCH-E, production-width)

`slice12_check.c` holds the brute-force reference, the check/barrett/recentre modes and the benchmark driver; it is
`#include`d at the end of `slice12.c` (one translation unit; the split only keeps each file under the package's
per-file limit).

**Layout** (proof.md 4.3): half-ellipsoid `{b : Q_b(b) <= R^2, b5 >= 1}` (`slice_mb.h` is the b-block of the shipped
form); lines k = b5, segments of 2^gs lines, each a separate allocation; packed rows, segment-local int64 cells, dense
per-line suffix metadata `cell = meta[k][idx4] + b0`, two uint16 sums per cell. Rows are built with radius
R^2 (1 + 2^-30) + 2^-20 so every walk point lies in its row. No bounding-box bytes, no global product, 64-bit counters.

**Paths.** *Line path* (p < 2^20): per line and coset, a 5-dimensional Schnorr-Euchner walk of the recentred coset
of L5(p, s) inside the line's ellipsoid with direct updates. *Record path* (2^20 <= p < 2^63): one 6-dimensional
walk per sigma-pair over the full ellipsoid with affine (idx4, b0, b5) states; hits with b5 < 0 folded, b5 = 0
skipped; records (segment-local cell uint64, global ideal id uint64) appended to the segment's retained store (grown,
never wrapped); *apply* adds a batch log to the recorded cells; *scan* tests f >= T_f and g >= T_g per threshold
block and marks survivors with the sentinel 0xFFFF; *recover* copies every record whose cell carries the sentinel.

**Arithmetic.** All modular setups use `bar_red`, a Barrett reduction valid for every p < 2^63 (k = bitlength p,
mu = floor(2^2k / p) < 2^(k+1), q = ((x >> (k-1)) mu) >> (k+1), two conditional subtractions; HAC 14.42); there is
no `%` in any sieving path. The batch log `lp = floor(16 log2 p)` is pair data: production computes it exactly at factor-base
construction (integer comparison of p^16 with powers of two); the benchmark/check driver fills it with libm `log2`. Coset origins are recentred by integer Babai rounding `q_i = llround(W o)` and
`o - sum q_i B_i` computed exactly in wrapping uint64 arithmetic (the exact value is small); any integer rounding
gives a valid representative, and mode `recentre` measures the deviation from exact (int128) rounding up to
p ~ 2^62 (<= 1 basis vector) and checks coset membership (0 failures).

**Correctness (`measurements12.txt`).** `slice12 check log2pts gs` runs both paths for several primes and both
sides, the apply pass, the scan with per-block thresholds and the recovery, and compares every cell of every segment
(digests of the complete f and g arrays), the survivor count and the recovered-record count (every recovered record is a true
hit at a surviving cell by construction) with a brute-force reference that tests each
point of the box against Q_b and the coset congruences directly:

| slice | segments | result |
|---|---|---|
| 2^18 cells (222,594) | 6 | f/g digests MATCH, survivors 1 = ref, recovered 1 = ref |
| 2^20 (922,823) | 4 and 1 | MATCH, 3 = 3, 3 = 3 |
| 2^22 (3,794,958) | 3 and 10 | MATCH, 11 = 11, 11 = 11 |
| 2^24 (15,511,016) | 2 | MATCH, 37 = 37, 37 = 37 |

`slice12 barrett 20000`: 1,220,000 reductions at every bit length 3..63, 0 mismatches. The check caught a defect of
the v10-style setup (it assumed |a_i| < p), fixed here. The kernel aborts on p >= 2^63 by design (X2 = 2^63).

**Measured costs** (`measure12.sh`, slice of 2^26 cells = 63,087,346, 5 segments, K = 17 lines; h = expected points
per coset: per *line* coset for the line path, per *slice* coset for the record path):

| path | h | Ir per unit | unit |
|---|---|---|---|
| line (incl. per-(pair, line, coset) setup) | 2^13.82 (p = 257) | **24.10** | update |
| line | 2^9.82 (p = 4099) | **31.11** | update |
| line | 2^5.82 (p = 65537) | **88.31** | update |
| line | 2^1.82 (p = 2^20) | **310.70** | update |
| line setup alone, per (pair, line, coset) | - | 488.24 | setup (c_k by Barrett, recentring, centre, affine states) |
| record (walk, fold, store, per-coset recentring/centre) | 2^18.91 (p = 257) | **44.90** | hit |
| record | 2^10.91 (p = 65537) | **78.21** | hit |
| record | 2^3.91 (p = 2^23) | **266.97** | hit |
| record | 2^0.91 (p = 2^26) | **938.08** | hit |
| record | 2^-1.09 (p = 2^28) | **2231.66** | hit (= 1050 per coset) |
| apply | - | **7.00** | record |
| recover | - | **7.003** | record |
| setup6, per (slice, pair): a(s) mod p (Horner, Barrett), c = a(s) r^-1 | p = 2^20, 2^63 - 25 | 406.87, 374.07 (charged 500 with the 36 streamed words) | pair |
| build (whole program: allocation, metadata, rows, first zero pass) | - | **12.325** per cell (27.27 per suffix-box entry; the benchmark has 0.45 entries and 0.065 rows per cell, production ~0.004 and ~0.0003) | cell |
| zero pass | - | **0.625** | cell |
| scan, 16-cell blocks / 4-cell blocks | - | **9.313 / 13.250** (= 8.0 per cell + 21 per block) | cell |

The model interpolates linearly in log2 h between the table points, charges the densest measured value above the
table (production line cosets have h >= 2^34 and record cosets h >= 2^-0.4), and 1600 per (slice, pair) for cosets
with h < 2^-1.09 (not reached at the operating point).

## `nbound12.c`: certified per-block norm bounds (the scan thresholds)

An operation-for-operation C port of `embedding_bounds.py`: four places of K_h as certified intervals (rational
Rouche disks, widened outward), every interval operation in binary64 round-to-nearest widened by one ulp on each side
(`nextafter` on the bit pattern), 16 log2 |N| bounded per place from a 17-entry table of 2^(j/16) bounds (no libm
log). Blocks are accepted at width <= W, halved otherwise, unresolved single points rejected (never survivors).
`nbound12 check nbound-vectors.txt` replays the Python reference's leaves and bounds on all 2,401 pool rows
(16-blocks: 4,221 evaluations, 3,311 leaves, 0 rejected, **0 mismatches**); the vectors are regenerated by
`gen_nbound_vectors.py` from the package (the pool rows come from their seeds). `nbound12 explore nbound-vectors.txt
64 128` reports the production blocks (aligned 64-cell blocks, W = 128): **0.04232 evaluations and 0.02932 leaves
per cell** (slice-weighted), 0 rejected points, 0 non-finite values, mean accepted width 5.95 bits at the samples;
`dump` prints the per-row statistics embedded in `pool12-1.b85`/`pool12-2.b85` (evidence12 [2] recomputes them with the reference).

Measured (callgrind, whole program at reps 1/3 over the 2,401 cases): per block evaluation **5,862.7 Ir** on the
16-blocks and 5,833.0 on the 64/128 blocks (the model charges 5,900); `row_setup` (suffix at the four places)
**4,120 Ir** per row (charged 4,200); `setup_a` (a(v) and its powers) 15,223 per a (charged 15,300 per slice).
Wider blocks or widths were measured and do not pay (L0 = 128, W = 192: 0.0371 evaluations per cell at a wider
survivor window).

## `ecm12.c`: the modular multiplication behind the cofactorisation and descent prices

CIOS Montgomery multiplication on 5 x 64-bit limbs (cofactors < 2^300; n odd < 2^319), 64-bit general-register
arithmetic with 128-bit products. `ecm12 check 20000`: 20,000 products (including (n-1)^2 cases and 200 moduli of
several sizes) against a shift-subtract exact reference, **0 mismatches**. `ecm12 bench`: **367.76 Ir per
product** (reps 100,000 vs 300,000). ECM curve operations are counted in these multiplications (`model12.ecm_instr`:
GMP-ECM B1/curve table, 1.44 B1 bits x 11 multiplications per bit, stage 2 as a second stage 1, curves x ln 100 for
99% per-prime success); for L > 5 limbs the price scales by (L/5)^2 (CIOS is quadratic in the limb count).

## `la_mac.s`: block-Wiedemann entry (C_ENTRY = 36)

The v10 kernel with 64-bit column indices (`ldr x13`): the branch-free loop `LBB0_2` retires 33 instructions per
entry (evidence12 [9] re-counts it from the shipped assembly); the per-row Barrett finish and the combination of the
four column-block accumulators are amortised over the row weight 100: 36 per entry.

## Reproduction

    gcc -O3 -fno-tree-vectorize -fno-tree-loop-distribute-patterns -o slice12 slice12.c -lm
    ./slice12 check 22 2 && ./slice12 barrett 20000 && ./slice12 recentre 20000
    gcc -O3 -fno-tree-vectorize -fno-tree-loop-distribute-patterns -o nbound12 nbound12.c -lm
    python3 -I gen_nbound_vectors.py . nbound-vectors.txt && ./nbound12 check nbound-vectors.txt && ./nbound12 explore nbound-vectors.txt 64 128
    gcc -O3 -fno-tree-vectorize -fno-tree-loop-distribute-patterns -o ecm12 ecm12.c && ./ecm12 check 20000
    sh measure12.sh     # callgrind, arm64

# Storage and addressing

Collection. A segment is one line c_5 = k of a slice, held as its L4 x W0 box; the record walk writes (cell, ideal) records to 1 MiB pages
per chunk (chunk = cell >> 30); each chunk owner applies its records to a 2^30-entry uint32 sum array and scans it against the thresholds.
Routing per retained record 10.01 instructions (slice walk) and 1.14 (run walk: measured, KERNELS.md); the former is the BLS12-381 package's measured page path, which does not
depend on the field. A run's slices live together (one walk per (run, ideal pair)): <= 1024 a-values; see proof.md section 10.

Descent nodes. Two uint16 arrays cost4H bytes at half count H; R128 gives less than 2^85.417380 bytes. Strict bin decrease and acyclic helper rules bound simultaneous descent frames by 62 plus boot/helper, at most 64. Keeping every ancestor array costs less than2^99 bytes. Main storage and ECM scratch fit the remaining2^100-byte whole-attack allowance; no unpriced freeing/recomputation is assumed.

Linear algebra. N <= 2^45.7866 < 2^64/320 elements of four 64-bit limbs: one owner. The n = 40 sequences advance in groups
of five; x_j and chi x_j share a 320-byte slot. Kernels (lib:kernels_asm.data; carried, four-limb straight-line code whose counts do not depend on
the modulus): owner_part 267 per owner and five vectors + 6 per call (55 per vector); _chi_mul 135 per vector (C_CHI; corrected truncated Shoup multiply by the constant chi, output in [0, 2r), listed in lib:kernels_asm.data, KERNELS.md); row reduction mean 54 (small_row_proof.txt); header 10,
routing 40, unit leaves 2: row tail 161 for one owner; entries 7.4 (unit) / 33 (general), mixed at the measured unit share (k5 = 8.46). Fallback _row_red (65+6) at the
BN254 r ([6b]): D = pos - neg + ceil(2^323/r) r < 2^324 + r, q = floor(floor(D/2^248) mu/2^79), mu = floor(2^327/r) (74 bits);
0 <= D/r - q < 1.146 < 2, so the output lies in [0, 2r). Horner and generator reductions: _row_redc 270 (Montgomery R = 2^320), generator
terms 97. Dense-band tables ([6g]): 2^55.04 bytes per 5-vector slot, held by the one linear-algebra owner in fixed arenas with its vectors and 16 matrix versions: < 2^62 bytes, inside the owner's 2^64-byte address space and the 2^100-byte budget.

Projection X^T v: C_PROJ 28 per coordinate and vector (analytic). Weight a = low 64 bits of an independent 128-bit draw (RANK-BRIDGE.md); one (offset, a) record serves two vectors: 45 instructions per record (proj2.s, lib:kernels_asm.data). y < 2r, a < 2^64, block < 2^37: sums < 2^356, six limbs, exact. Plus 100,000 per block output (call, stores, reduction to [0, r)) and 65,536 per coordinate (draws, records): < 2^72.5 in all, below 5.5 per term (> 2^94).
Flags: byte/row, M0>=2^32; LA charges 2^60 setup.

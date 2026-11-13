# Storage and addressing

Collection. A segment is one line c_5 = k of a slice, held as its L4 x W0 box; the record walk writes (cell, ideal) records to 1 MiB pages
per chunk (chunk = cell >> 30); each chunk owner applies its records to a 2^30-entry uint32 sum array and scans it against the thresholds.
Routing per retained record 10.01 instructions (slice walk) and 1.14 (run walk: measured, KERNELS.md); the former is the BLS12-381 package's measured page path, which does not
depend on the field. A run's slices live together (one walk per (run, ideal pair)): <= 1024 a-values; see proof.md section 10.

Descent nodes. Two uint16 arrays cost4H bytes at half count H; R128 gives less than 2^85.417380 bytes. Strict bin decrease and acyclic helper rules bound simultaneous descent frames by 62 plus boot/helper, at most 64. Keeping every ancestor array costs less than2^99 bytes. Main storage and ECM scratch fit the remaining2^100-byte whole-attack allowance; no unpriced freeing/recomputation is assumed.

Linear algebra. N <= 2^46.2658 < 2^64/320 elements of four 64-bit limbs: one owner. The n = 40 sequences advance in groups
of five; x_j and chi x_j share a 320-byte slot. Kernels (kernels_asm.b85; carried, four-limb straight-line code whose counts do not depend on
the modulus): owner_part 267 per owner and five vectors + 6 per call (55 per vector); _chi_mul 200 (C_CHI); _row_red 65 + 6; header 10,
routing 40: row tail 176 for one owner; entries 7.8 (unit) / 33 (general), mixed at the measured unit share (k5 = 8.0). _row_red at the
BN254 r ([6b]): D = pos - neg + ceil(2^323/r) r < 2^324 + r, q = floor(floor(D/2^248) mu/2^79), mu = floor(2^327/r) (74 bits);
0 <= D/r - q < 1.146 < 2, so the output lies in [0, 2r). Horner and generator reductions: _row_redc 270 (Montgomery R = 2^320), generator
terms 97. Dense-band tables ([6g]): 2^55.04 bytes per 5-vector slot, 8 slots 2^58.04 bytes, held by the one linear-algebra owner with its vectors (N x 320 x 8 < 2^58.0 bytes) and the matrix (< 2^57.31): < 2^59.41 bytes, inside the owner's 2^64-byte address space and the 2^100-byte budget.

# Storage and addressing

Collection. A segment is one line c_5 = k of a slice, held as its L4 x W0 box; the record walk writes (cell, ideal) records to 1 MiB pages
per chunk (chunk = cell >> 30); each chunk owner applies its records to a 2^30-entry uint32 sum array and scans it against the thresholds.
Routing per retained record 10.01 instructions (slice walk) and 3.01 (run walk): the BLS12-381 package's measured page path, which does not
depend on the field. A run's slices live together (one walk per (run, ideal pair)): <= 1024 a-values; see proof.md section 10.

Descent nodes. A node's sieve emits (owner, offset, 16 log2 p) records into two uint16 score arrays per owner (scores < 2^16 at every node,
[7]); a node of half count H uses 4H bytes: the full-plan maximum (R=160) is below2^89.281 bytes, within the2^91 whole-attack allowance.

Linear algebra. N <= 2^44.9415 < 2^64/320 elements of four 64-bit limbs: one owner. The n = 40 sequences advance in groups
of five; x_j and chi x_j share a 320-byte slot. Kernels (kernels_asm.b85; carried, four-limb straight-line code whose counts do not depend on
the modulus): owner_part 267 per owner and five vectors + 6 per call (55 per vector); _chi_mul 200 (C_CHI); _row_red 65 + 6; header 10,
routing 40: row tail 176 for one owner; entries 7.8 (unit) / 33 (general), mixed at the measured unit share (k5 = 8.38). _row_red at the
BN254 r ([6b]): D = pos - neg + ceil(2^323/r) r < 2^324 + r, q = floor(floor(D/2^248) mu/2^79), mu = floor(2^327/r) (74 bits);
0 <= D/r - q < 1.146 < 2, so the output lies in [0, 2r). Horner and generator reductions: _row_redc 270 (Montgomery R = 2^320), generator
terms 97.

# Storage and addressing

Objects are bounded and addressed by (owner, local index); store16.py is the reference ([11]).

## Relation collection

A segment is one line b5 = k of a slice; a line has <= 2^56.807 cells
(store16.line_cells_bound), so a cell index fits 58 bits. Rows are packed
suffix-lexicographically (cell = row base + b0) with 16-byte row descriptors in
2^16-entry pages.

The record walk writes each record (cell: uint64, ideal: uint64) directly to
the current 1 MiB page (2^16 records) of its chunk, chunk = cell >> 30 (chunk id
< 2^29), through a 2^29-entry tail table (4 GiB): 7 instructions on the hot path (shift, tail load, add, paired store, page-end test and branch,
tail store). Full pages go to the chunk's two-level page directory
(128-bit record ordinals). Records stay in emission order (store16.scatter vs store16.route). Each chunk owner applies its records to its
2^30-entry uint32 sum array (4 GiB), scans it against the thresholds and
recovers survivors' records. Open pages hold <= 2^29 x 1 MiB = 2^49
bytes per line in flight.

Price per retained record, end to end: scatter (7, above), page hand-off, the owner's page copy (16 bytes per record,
3 instructions) and its walk into apply; executed instructions under qemu-aarch64 8.2.2 (gcc 13.3 -O3 -fno-tree-vectorize -mgeneral-regs-only),
differential 2^18 -> 2^19 records: 26.0009 per record minus the source loop 7 and
apply 9 = 10.0009 (37.2 per 2^16-record page), charged 10.01.
Plus 0.5 per cell for every other per-cell step (measured 0.096: CELLRES, KERNELS.md); apply
5.76 and recovery 5.76 per record are separate terms (slice12oh.c); scan2 zeroes the f sum array as it
reads them, and each survivor's SENTINEL is cleared after recovery (8 instructions per survivor, charged in the scan row's
per-leaf part).
Primes in [2^63, 2^64) of contaminated runs take the wide path of KERNELS.md, with the same record layout.
Primes in [2^48, 2^64) on a run (KERNELS.md, rec7): a run's slices live together (one walk per (run, ideal) emits
records into every slice's chunks): <= 1,024 a-values (RUN_T), 2^76.0 cells (slices <= vb (17/20)^3 x 1.0326), 4 bytes of
sums and 0.575 expected 16-byte records per cell: below 2^79.8 bytes. The six 7-dim bases per ideal take < 2 KiB each: < 2^73.2 bytes.

## Descent nodes

A node's sieve emits records (u64 owner = index >> 40, u64 offset = index mod 2^40, u16 16 log2 p; separate f and g
streams; dsieve16 rec_t) to the owner, which adds into its two 2^40-entry uint16 score arrays (a score is <= 16 x
the norm size, below 2^16 at every node). All owners hold their arrays at once: a node of half count H uses 4H bytes, at
R = 384 (H < 2^102.44) 2^104.44 bytes over < 2^62.44 owners, offsets < 2^40. The kernel refuses balls whose half count
reaches 2^103; [7] audits every sieved node (largest score 32,709) and store16.check round-trips the split at the
boundaries, including index 2^94. Records stream through unstored; transfer at 64-bit owner ids: four 16-bit radix passes (64) + append 6 + flush 4 + directory 2 + dispatch 8 = 84 per
record, inside SIEVE_HIT. 2^40 workers hold the exact ball index tables (12 x (R^2 + 1) x (2R + 1) x 16 bytes =
2^34.34 at R = 384, 2^74.4 in all); nodes run singly, so the descent's peak is 2^104.44 + 2^74.4 < 2^105 bytes.

## Linear algebra

N <= 2^57.3256 < 6 x 2^55 elements of four 64-bit limbs. The n = 17 sequences advance in groups of five, element j of
a group's five vectors in 256-byte slot 2j and chi x_j in slot 2j + 1; slot s lives in block s >> 56 (six active owners). Rows
store entries grouped by (block, sign) over slots; each owner sends a (pos, neg) partial sum up a fixed tree to the row
owner, which adds the 10 dense Schirokauer terms and reduces with _row_redc: (pos - neg + dense) 2^-320 mod r, the same
factor on every row, so BW runs on 2^-320 A (same kernel). Kernels (kernels_asm.b85) are counted by [10].

Krylov passes. A +-1 group runs unit_run5 (one index load and address per entry, then per vector two ldp, the
four-limb add and the carry count: 39 instructions per entry for the five vectors, all 25 accumulators in registers);
other groups run gen_run1 once per vector (33 per vector-entry, stride 256). Block selection and sign (4 + 3) are paid
once per +-1 entry, per vector-entry otherwise; an entry a + b chi is stored as |a| unit entries on slot 2j and |b| on
slot 2j + 1: k5 charged 9.88 (C_ENTRY_K, PREREG-U, section 9). Each new element's chi x_j is _chi_mul (Montgomery with
chi 2^256 mod r, exact): args 3 + call 1 + 194 = 198 <= 200 (C_CHI, [10]).

Row tail (six owners). Five sends (_send_pair 15 + 18 transfer) and Five combinations (_combine 62 + 2), _row_redc 260 +
args and call 6 (pos, neg < 2^323, dense sum < 2^516, output < 2^256), headers 10 per owner (60), routing 40, group
overhead 168: 1,019, charged 1,020 (C_ROW, owner10.row_tail). Horner and generator reductions are _row_redc(0, 0, s) on
Montgomery forms (R = 2^320): 266, charged 270 (C_RED9); dense terms 97 (bw_full). _chi_mul and _row_redc are
straight-line (callgrind 194 and 260 Ir per call); 50,007 and 60,007 random and extreme inputs match Python.

## Aggregate

Peak storage is <= 2^105 bytes. The main phase
stays below 2^83 bytes: a run's sum arrays and records (< 2^79.8), rec7 bases, active record pages, stored
relations (< 2^79.2), the sparse matrix (< 2^70.5), Krylov block (< 2^68), Horner/projection vectors (< 2^76),
factor base (< 2^69.5), canonical-list staging 2^78, P_g, a g wave and their trees (< 2^74), and smaller tables.

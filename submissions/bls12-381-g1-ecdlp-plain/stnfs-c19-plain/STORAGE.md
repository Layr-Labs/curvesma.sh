# Storage and addressing

Objects are bounded and addressed by (owner, local index); store16.py is the reference ([11]).

## Relation collection

A segment is one line b5 = k of a slice, held as its L4 x W0 box (cell = idx4 W0 + r0 + b0; <= 2^59.533 positions,
store16.line_box_bound), so a cell index fits 60 bits; 16-byte row descriptors in 2^16-entry pages mark the region cells,
the only ones written.

The record walk writes each record (cell: uint64, ideal: uint64) directly to
the current 1 MiB page (2^16 records) of its chunk, chunk = cell >> 30 (chunk id
< 2^30), through a 2^30-entry tail table (8 GiB): 7 instructions on the hot path (shift, tail load, add, paired store, page-end test and branch,
tail store). Full pages go to the chunk's two-level page directory
(128-bit record ordinals). Records stay in emission order (store16.scatter vs store16.route). Each chunk owner applies its records to its
2^30-entry uint32 sum array (4 GiB), scans it against the thresholds (no recovery). Open pages hold <= 2^30 x 1 MiB = 2^50
bytes per line in flight.

Price per retained record, end to end: scatter (7), page hand-off, owner page copy and walk into apply; executed instructions under qemu-aarch64 8.2.2 (gcc 13.3 -O3): differential 2^18 -> 2^19 records gives 10.0009 net. In-place read drops the 3-instruction copy, charged 7.01. Plus 0.02/cell (CELLRES, KERNELS.md); apply 5.76 per record is separate.
Contaminated runs' primes in [2^63, 2^64) take the wide path (same records).
Primes in [2^48, 2^64) on a run (KERNELS.md, rec7): a run's slices live together (one walk per (run, ideal) emits
records into every slice's chunks): <= 1,024 a-values (RUN_T), 2^76.0 cells (slices <= vb (17/20)^3 x 1.0326) in boxes of <= 20.35 x
their cells, 4 bytes per position and 0.575 16-byte records per cell: below 2^82.4 bytes. Each slice adds a (Kl+1) x 16-byte line table (< 2^16 bytes). The six 7-dim bases per ideal take < 2 KiB each: < 2^73.2 bytes.

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

N <= 2^57.4186 < 4 x 2^64/320 elements of four 64-bit limbs. The n = 17 sequences advance in groups of five; x_j and
chi x_j share 320-byte slot j (bytes 320 j, 320 j + 160; vector v at + 32 v), 2^64/320 slots per owner (four owners).
A row stores per owner four groups (sign x unit/general) of offsets; each owner sends a (pos, neg) partial sum up
a fixed tree to the row owner, which adds the 10 dense Schirokauer terms and reduces with _row_redc: (pos - neg + dense)
2^-320 mod r, the same factor on every row, so BW runs on 2^-320 A. Kernels (kernels_asm.b85) are counted by [10].

Krylov passes. owner_part runs unit_run5 per unit group (39 per entry for five vectors, accumulators in registers)
and gen_run1 five times per general group (33 per vector-entry). Block and sign are fixed per group: 7.8 per unit vector-entry, 33 general;
a + b chi is |a| unit entries on x_j, |b| on chi x_j: k5 charged 8.38 at excess 1.04 (C_ENTRY_K, K5_RUNS, section 9). owner_part
with calls' prologues/epilogues: 267 per owner and five vectors (straight-line, [10]; callgrind 267 + 39 per unit
+ 165 per general entry; 6,000 random/extreme parts match Python), + 6 args/call: 55 per vector. Each new chi x_j is
_chi_mul (Montgomery, chi 2^256 mod r, exact): 3 args + call + 194 = 198 <= 200 (C_CHI, [10]).

Row tail (four owners): 3 sends (_send_pair 15 + 18 transfer), 3 combines (_combine 62 + 2), _row_redc 260 + 6
args/call (pos, neg < 2^323, dense < 2^516, output < 2^256), headers 10/owner (40), routing 40, owner_part 4 x 55: 857,
charged 862 (C_ROW, owner10.row_tail). Horner/generator reductions: _row_redc(0, 0, s) on Montgomery forms (R = 2^320),
266, charged 270 (C_RED9); dense terms 97 (bw_full). Both straight-line (callgrind 194, 260 Ir/call); 50,007 / 60,007
random and extreme inputs match Python.

## Aggregate

Peak storage is <= 2^105 bytes. The main phase
stays below 2^83 bytes: a run's sum arrays and records (< 2^82.4), rec7 bases, active record pages, stored
relations (< 2^79.2), the sparse matrix (< 2^70.5), Krylov block (< 2^68), Horner/projection vectors (< 2^76),
factor base (< 2^69.5), canonical-list staging 2^78, P_g, P2, a wave and their trees (< 2^74), and smaller tables.

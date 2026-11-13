# Storage and addressing

Objects are bounded and addressed by (owner, local index); store16.py is the reference ([11]).

## Relation collection

A segment is one line b5 = k of a slice, held as its L4 x W0 box (cell = idx4 W0 + r0 + b0; <= 2^59.639 positions,
store16.line_box_bound), so a cell index fits 60 bits; 16-byte row descriptors in 2^16-entry pages mark the region cells,
the only ones written.

The record walk writes each record (cell: uint64, ideal: uint64) directly to
the current 1 MiB page (2^16 records) of its chunk, chunk = cell >> 30 (chunk id
< 2^30), through a 2^30-entry tail table (8 GiB): 7 instructions on the hot path (shift, tail load, add, paired store, page-end test and branch,
tail store). Full pages go to the chunk's two-level page directory
(128-bit record ordinals). Records stay in emission order (store16.scatter vs store16.route). Each chunk owner applies its records to its
2^30-entry uint16 sum array (2 GiB), scans it against the thresholds (no recovery). Open pages hold <= 2^30 x 1 MiB = 2^50
bytes per line in flight.

Price per retained record, end to end: scatter (7, above), page hand-off, the owner's page copy (16 bytes per record,
3 instructions) and its walk into apply; executed instructions under qemu-aarch64 8.2.2 (gcc 13.3 -O3 -fno-tree-vectorize -mgeneral-regs-only),
differential 2^18 -> 2^19 records: 26.0009 per record minus the source loop 7 and
apply 9 = 10.0009, charged 10.01. Run-walk records (8 bytes: cell mod 2^32, 16 log2 p) are scattered inside the measured walk (KERNELS.md);
copy and hand-off 1.07, apply 3.13 (lean8).
Plus 0.02 per cell (CELLRES, KERNELS.md); slice-walk apply 5.76 per record is separate (slice12oh.c); scan2 resets each f sum it
reads to H, each survivor is emitted with its sum s (recheck) and its SENTINEL is cleared after the scan (8 instructions per survivor, charged in the scan row's
per-leaf part).
Contaminated runs' primes in [2^63, 2^64) take the wide path (same records).
Stage B: rechecked survivors (N_f, z, L, s, cell: 25 words) wait in waves of 2^15; wave, trees and P24 (3.1 MB): < 2^27 bytes; a worker's last partial wave is misc.
Primes in [2^48, 2^64) on a run (KERNELS.md, rec7): a run's slices live together (one walk per (run, ideal pair) emits
records into every slice's chunks): <= 2,048 a-values (RUN_T), 2^77.1 cells (slices <= vb (1751/2000)^3 x 0.999595) in boxes of <= 20.35 x
their cells, 4 bytes per position and 0.575 16-byte records per cell: below 2^83.6 bytes. Each slice adds a (Kl+1) x 16-byte line table (< 2^16 bytes) at offset 2^18 t of the run's table: address space, resident in the 4 KiB pages touched (<= 2^16 + 2^12 bytes per a-value; with huge pages up to the 2^18 reserved); run-walk pages are lower halves of 2 MiB windows (upper halves never written); its lines' tail tables and open pages (above) add < 2^72.1 bytes per run. T enters the walk only as the slab width (< 2^11) and the number of line tables; cell and tail widths do not change. The six 7-dim bases per ideal take < 2 KiB each: < 2^73.2 bytes.

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

N <= 2^56.7309 < 3 x 2^64/320 elements of four 64-bit limbs. The n = 300 sequences advance in groups of five; x_j and
chi x_j share 320-byte slot j (bytes 320 j, 320 j + 160; vector v at + 32 v), 2^64/320 slots per owner (four owners).
A row stores per owner four groups (sign x unit/general) of offsets; each owner sends a (pos, neg) partial sum up
a fixed tree to the row owner, which reduces with _row_red: pos - neg mod r in [0, 2r); the Schirokauer columns are right-hand sides (proof.md section 4). Kernels (kernels_asm.b85) are counted by [10].

Krylov passes. owner_part runs unit_run5 per unit group (39 per entry for five vectors, accumulators in registers)
and gen_run1 five times per general group (33 per vector-entry). Block and sign are fixed per group: 7.8 per unit vector-entry, 33 general;
a + b chi is |a| unit entries on x_j, |b| on chi x_j: k5 charged 8.29 (C_ENTRY_K, PREREG-U, section 9). owner_part
with calls' prologues/epilogues: 267 per owner and five vectors (straight-line, [10]; callgrind 267 + 39 per unit
+ 165 per general entry; 6,000 random/extreme parts match Python), + 6 args/call: 55 per vector. Each new chi x_j is
_chi_mul (Montgomery, chi 2^256 mod r, exact): 3 args + call + 194 = 198 <= 200 (C_CHI, [10]).
Dense-band tables (proof.md section 4): each group's 4th owner (free >= 2^64 - 320 max(N - 3 x 2^64/320, 2^28) > 2^63.99 bytes), 2^47.71 x 160 bytes (60 groups: < 2^61), by red_row_u (tab.txt: 18.40 Ir per element).

Row tail (four owners): 3 sends (_send_pair 15 + 18 transfer), 3 combines (_combine 32 + 2), _row_red 65 + 6
args/call (pos, neg < 2^323, output < 2r), headers 10/owner (40), routing 40, owner_part 4 x 55: 572 (C_ROW, owner10.row_tail).
_row_red (one Barrett step): D = pos - neg + ceil(2^323/r) r < 2^324 + r, q = floor(floor(D/2^248) mu/2^79), mu = floor(2^327/r), output
D - q r mod 2^256; 0 <= D/r - q < 1 + 2^248/r + D/2^327 < 1.14, so the output is D - q r, in [0, 2r). Horner/generator reductions:
_row_redc(0, 0, s) on Montgomery forms (R = 2^320), 266, charged 270 (C_RED9); generator terms 97 (bw_full). All straight-line (callgrind
65, 194, 260 Ir/call); 800,000 / 50,007 / 60,007 random and extreme inputs match Python.

## Aggregate

Peak storage is <= 2^105 bytes. The main phase
stays below 2^86 bytes: a run's sum arrays and records (< 2^83.6), rec7 bases, active record pages, stored
relations (< 2^80), the sparse matrix (< 2^70.8), Krylov block (< 2^71), Horner/projection vectors (< 2^76), the generator's sequence and Kronecker workspace (< 2^85),
factor base (< 2^70), canonical-list staging 2^78, P_g, P2, a wave and their trees (< 2^74), and smaller tables.

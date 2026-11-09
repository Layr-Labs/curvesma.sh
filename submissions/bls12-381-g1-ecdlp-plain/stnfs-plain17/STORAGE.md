# Storage and addressing

Objects are addressed by (owner, local index); store16.py is the reference ([11]).

## Relation collection

A segment is one line b5 = k of a slice, held as its L4 x W0 box (cell = idx4 W0 + r0 + b0; <= 2^59.533 positions,
store16.line_box_bound), so a cell index fits 60 bits; 16-byte row descriptors in 2^16-entry pages mark the region cells,
the only ones written.
The record walk writes each record (cell, ideal: uint64) to the current 1 MiB page (2^16 records) of chunk cell >> 30 (< 2^30) via a 2^30-entry tail table (8 GiB): 7 instructions (shift, tail load, add, paired store, page-end test and branch, tail store). Full pages go to the chunk's two-level page directory
(128-bit record ordinals). Records stay in emission order (store16.scatter vs store16.route). Each chunk owner applies its records to its
2^30-entry uint32 sum array (4 GiB), scans it against the thresholds (no recovery). Open pages: <= 2^50 bytes per line.
Price per retained record: scatter (7), page hand-off and the owner's walk into apply; qemu-aarch64 8.2.2 (gcc 13.3 -O3 -fno-tree-vectorize -mgeneral-regs-only), differential 2^18 -> 2^19 records: 26.0009 per record minus source loop 7 and apply 9 = 10.0009. The owner applies in place from the handed-off page (store16.inplace_apply, [11]), with no 16-byte page copy (3 instructions). route_inplace.txt ships both owner loops: copy 4 instructions (ldp, stp, subs, b.ne), apply loops identical (7 each), so the 3.0 saving is conservative; [10i] checks counts and charge. Whole pages hand off; a chunk holds two pages and the writer takes another (run, ideal) while both are full (< 2^51 bytes per line, inside the aggregate), so handed-off pages are stable during apply and walk work only changes order. Charged 7.01 per retained record.
Plus 0.02 per cell (CELLRES, KERNELS.md); apply 5.76 per record is separate (slice12oh.c); scan2 resets each f sum it
reads to H, and each survivor's SENTINEL is cleared after the scan (8 instructions per survivor, charged in the scan row's
per-leaf part).
Primes in [2^48, 2^64) on a run (KERNELS.md, rec7): a run's slices live together (one walk per (run, ideal) emits
records into every slice's chunks): <= 1,024 a-values (RUN_T), 2^76.0 cells (slices <= vb (17/20)^3 x 1.0326) in boxes of <= 20.35 x
their cells, 4 bytes per position and 0.575 16-byte records per cell: < 2^82.51 bytes. Each slice adds a (Kl+1) x 16-byte line table (< 2^16 bytes). The six 7-dim bases per ideal take < 2 KiB each: < 2^73.2 bytes.

## Descent nodes

A node's sieve sends records (owner = index >> 40, offset mod 2^40, u16 16 log2 p; f, g streams) to owners' 2^40-entry uint16 score arrays (scores < 2^16): 4H bytes for half count H, 2^104.44 at R = 384 over < 2^62.44 owners. The kernel refuses half counts >= 2^103; [7] audits every sieved node (largest score 32,709); store16.check round-trips the split. Records stream through unstored; transfer at 64-bit owner ids: four 16-bit radix passes (64) + append 6 + flush 4 + directory 2 + dispatch 8 = 84 per
record, inside SIEVE_HIT. Ball index tables: 2^34.34 bytes at R = 384 per worker, 2^74.4 in all; nodes run singly.

## Linear algebra

N <= 2^57.2617 < 3 x 2^64/320 elements of four 64-bit limbs. The n = 25 sequences advance as five five-vector groups, one at a time; x_j and
chi x_j share 320-byte slot j (vector v at + 32 v), 2^64/320 slots per owner (three owners).
A row stores per owner four offset groups (sign x unit/general); owners send (pos, neg) partial sums up a fixed tree to the row owner, which adds the 10 dense Schirokauer terms and reduces with _row_redc: (pos - neg + dense) 2^-320 mod r, so BW runs on 2^-320 A. Kernels (kernels_asm.b85) are counted by [10].
Krylov passes. owner_part runs unit_run5 per unit group (39 per entry for five vectors, accumulators in registers)
and gen_run1 five times per general group (33 per vector-entry). Block and sign are fixed per group: 7.8 per unit vector-entry, 33 general;
a + b chi is |a| unit entries on x_j, |b| on chi x_j: k5 charged 8.38 (C_ENTRY_K, PREREG-U, section 9). owner_part
with calls' prologues/epilogues: 267 per owner and five vectors (straight-line, [10]; callgrind 267 + 39 per unit
+ 165 per general entry; 6,000 random/extreme parts match Python), + 6 args/call: 55 per vector. Each new chi x_j is
_chi_mul (Montgomery, chi 2^256 mod r, exact): 3 args + call + 194 = 198 (C_CHI, [10]).
Row tail (three owners): 2 sends (_send_pair 15 + 18 transfer), 2 combines (_combine 62 + 2), _row_redc 260 + 6
args/call (pos, neg < 2^323, dense < 2^516, output < 2^256), headers 10/owner (30), routing 40, owner_part 3 x 55: 695,
charged 700 (C_ROW, owner10.row_tail). Horner/generator reductions: _row_redc(0, 0, s) on Montgomery forms (R = 2^320),
266 (C_RED9); dense terms 97 (bw_full). Both straight-line (callgrind 194, 260 Ir/call); 50,007 / 60,007
random and extreme inputs match Python.

## Aggregate

Peak <= 2^105 bytes (descent). Main phase < 2^83: run sums and records < 2^82.51, relations < 2^79.2, canonical staging 2^78, Horner/projection < 2^76, P_g, P2, a wave and trees < 2^74, matrix < 2^70.5, factor base < 2^69.5, Krylov (one group) < 2^68, checkpoints < 2^69, smaller tables. The generator runs alone: k^2 entries of 18L limbs, < 2^78 bytes.

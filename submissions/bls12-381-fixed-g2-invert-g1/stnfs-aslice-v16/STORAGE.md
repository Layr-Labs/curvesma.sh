# Storage and addressing

Every object is individually bounded and addressed by (owner, local index); no
aggregate byte offset is formed. `store16.py` is the executable reference;
[11] runs its checks at production widths.

## Relation collection

A segment is one line b5 = k of a slice. A line has at most 2^57.508 cells (exact
determinant bound in `store16.line_cells_bound`), so a cell index fits 59 bits.
Rows are packed suffix-lexicographically (cell = row base + b0) with 16-byte row
descriptors in 2^16-entry pages.

The line's record walk emits records (cell: uint64, ideal: uint64) into a 64 KiB
line buffer, flushed to 1 MiB pages of 2^16 records (two-level page directory,
128-bit record ordinals). Records are routed to chunk = cell >> 30 (chunk id
< 2^29) by a two-pass stable radix sort with 16-bit digits. Each chunk owner
applies its records into its own 2^30-entry uint32 sum array (4 GiB), scans it
against the thresholds, and recovers the records of survivor cells. Dropped
blocks (block selection) are skipped by both walkers.

Price per retained record: radix routing 16 (two passes), append 6, page
flush 4, directory 2, sentinel and dispatch 8: **36**, plus **8 per cell** for
row bookkeeping; apply **7.00** and recovery **7.003** per record are separate
ledger terms. At 2.29 records per cell this is part of the 541.7
instructions per cell of proof.md section 3.

## Descent nodes

A node's sieve emits (owner = index >> 30, offset, 16 log2 p) into the owner's
input stream; the owner adds into its 2^30-entry score array. At the largest
ball (R = 128, half count < 2^84) owners are < 2^54 and offsets < 2^30. The
transfer is charged at 32 instructions per hit inside SIEVE_HIT (KERNELS.md).

## Linear algebra

N = 2^60.2341 elements of four 64-bit limbs (32 bytes). Element j lives in block
j >> 56 at local index j mod 2^56: each block is 2^61 bytes, 19 blocks are
needed and 20 are allocated. Rows store entries grouped by block with 64-bit
column indices and coefficients stored as a sign bit and an 8-bit magnitude
(|c| <= 255). Block owners form
positive and negative five-limb partial sums (100 entries x 255 x r < 2^270) and
send them over a balanced tree of the 20 owners (five hops) to the row owner,
which finishes with Barrett reduction. Price per entry: 33 (`la_mac.s`) + block
selection 4 + sign 3 + transfers 32 + finishing and row headers 16 = 88,
charged as **96**.

## Aggregate

At most 2^81 bytes in total: chunk sum arrays and record pages of the active
lines, the paged sparse matrix, 20 vector blocks per vector for <= 1024 block
Wiedemann vectors held in paged form, canonical-list staging 2^78 bytes, code
and scratch.

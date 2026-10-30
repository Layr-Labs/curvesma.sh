# Storage and addressing

Every object is individually bounded and addressed by (owner, local index);
`store16.py` is the executable reference and [11] runs its checks.

## Relation collection

A segment is one line b5 = k of a slice; a line has at most 2^57.508 cells
(`store16.line_cells_bound`), so a cell index fits 58 bits. Rows are packed
suffix-lexicographically (cell = row base + b0) with 16-byte row descriptors in
2^16-entry pages.

The record walk emits records (cell: uint64, ideal: uint64) into a 64 KiB line
buffer, flushed to 1 MiB pages of 2^16 records (two-level page directory,
128-bit record ordinals), routed to chunk = cell >> 30 (chunk id < 2^29) by a
two-pass stable radix sort with 16-bit digits. Each chunk owner applies its
records into its 2^30-entry uint32 sum array (4 GiB), scans it against the
thresholds and recovers the survivors' records.

Price per retained record: radix routing 16, append 6, page flush 4, directory
2, sentinel and dispatch 8: **36**, plus **8 per cell** for row bookkeeping;
apply **7.00** and recovery **7.003** per record are separate terms (2.29
records per cell inside the 407.0 instructions per cell of proof.md).

## Descent nodes

A node's sieve emits records (u64 owner = index >> 40, u64 offset = index mod
2^40, u16 16 log2 p; separate f and g streams; dsieve16 rec_t) to the owner,
which adds into its two 2^40-entry uint16 score arrays (a score is at most
16 x the norm size, below 2^16 at every node). All owners hold
their arrays at once: a node of half count H uses 4H bytes, at R = 384
(H < 2^102.44) 2^104.44 bytes over < 2^62.44 owners, offsets < 2^40. The kernel
refuses any ball whose half count reaches 2^103; [7] audits every sieved node
(largest score 32,709) and store16.check round-trips the split at the
boundaries, including index 2^94. Records stream through and are not stored;
the transfer is priced by the relation-collection rule at 64-bit owner ids:
four 16-bit radix passes (64) + append 6 + flush 4 + directory 2 + dispatch 8
= **84** per record, inside SIEVE_HIT (KERNELS.md). Enumeration runs on 2^40
workers holding the exact ball index tables (12 x (R^2 + 1) x (2R + 1) x 16
bytes = 2^34.34 bytes at R = 384, 2^74.4 in all). Nodes run one at a time, so
the descent's peak is 2^104.44 + 2^74.4 < 2^105 bytes.

## Linear algebra

N = 2^59.7306 elements of four 64-bit limbs (32 bytes); element j lives in block
j >> 56 at local index j mod 2^56 (blocks of 2^61 bytes, 14 needed, 20
allocated). Rows store entries grouped by block with 64-bit column indices and
coefficients as a sign bit and an 8-bit magnitude. Block owners form positive
and negative five-limb partial sums (101 entries x 255 x r < 2^270) and send
them over a balanced tree of the 20 owners (five hops) to the row owner, which
finishes with Barrett reduction. Price per entry: 33 (`la_mac.s`) + block
selection 4 + sign 3 + transfers 32 + finishing and row headers 16 = 88,
charged as **88** (the component sum itself, replacing the cushioned 90); the
full-field multiply-accumulate inner loop is 97 measured, charged **97**
(replacing the cushioned 100).

## Aggregate

Peak storage is at most 2^105 bytes (the largest descent node). The main phase
stays below 2^81 bytes: active chunk sum arrays and record pages, the stored
relations (2^64.16 x <= 256 bytes = 2^72.2), the paged sparse matrix (2^59.77
rows x 101 entries x 9 bytes < 2^69.6), 20 vector blocks of 2^61 bytes for
<= 1024 block Wiedemann vectors (2^75.3), canonical-list staging 2^78 bytes.

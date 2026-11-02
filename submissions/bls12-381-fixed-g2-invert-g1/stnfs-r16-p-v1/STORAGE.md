# Storage and addressing

Every object is individually bounded and addressed by (owner, local index);
`store16.py` is the executable reference and [11] runs its checks.

## Relation collection

A segment is one line b5 = k of a slice; a line has at most 2^57.508 cells
(`store16.line_cells_bound`), so a cell index fits 58 bits. Rows are packed
suffix-lexicographically (cell = row base + b0) with 16-byte row descriptors in
2^16-entry pages.

The record walk writes each record (cell: uint64, ideal: uint64) directly to
the current 1 MiB page (2^16 records) of its chunk, chunk = cell >> 30 (chunk id
< 2^29), through a 2^29-entry tail table (4 GiB): 7 instructions on the hot path (shift, tail load, add, paired store, page-end test and branch,
tail store). A full page goes to the chunk's two-level page directory
(128-bit record ordinals). Records of a chunk stay in emission order, so this
is the same stable chunk order as the source's line buffer plus two-pass 16-bit
radix sort (`store16.scatter` is checked against `store16.route`), with no
line buffer and no radix passes. Each chunk owner applies its records into its
2^30-entry uint32 sum array (4 GiB), scans it against the thresholds and
recovers the survivors' records. Open pages hold at most 2^29 x 1 MiB = 2^49
bytes per line in flight.

Hot path (AArch64; x0 tail table, x1 cell, x2 ideal; page-full branch rare):

```
emit_hot:
	lsr	x19, x1, #30
	ldr	x9, [x0, x19, lsl #3]
	add	x8, x9, #16
	stp	x1, x2, [x9]
	tst	x8, #0xfffff
	b.eq	emit_page
	str	x8, [x0, x19, lsl #3]
emit_end:
```

Whole route, measured (R3; measurements12.txt "R3" lines; bench source in base 5f39af0;
audit check `route_measured`): emission scatter as above; a full page is appended
to its chunk's page directory and a fresh page taken; the owner copies every page
(16 bytes per record, ldp/stp, 12 per 4 records: the transfer is charged as a full
copy) into its receive area, kept for recovery, and walks its directory into
apply. qemu-aarch64 executed-instruction counts at n and 2n records: route + source
+ apply 26.0009 per record, source loop 7 (static), apply 9.000, so the route is
**10.0009** per record at 1 MiB pages (page term 37.7 per page, 0.0006 per record).
A line's last partial page per chunk costs under 2^-20 per cell (every line holds
at least 2^54.58 cells), inside the 8 per cell. Price per retained record: **10.005**
(measured 10.0015 including the page term, 0.035% margin; R5 charged 10.01), plus **8 per cell** for row
bookkeeping; apply **7.00** and recovery **7.003** are separate terms, and the
measured record path (REC_TAB) is charged in full on top.

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

N = 2^59.4237 elements of four 64-bit limbs (32 bytes); element j lives in block
j >> 56 at local index j mod 2^56 (blocks of 2^61 bytes, 14 needed, 20
allocated). Rows store entries grouped by block with 64-bit column indices and
coefficients as a sign bit and an 8-bit magnitude. Block owners form positive
and negative five-limb partial sums and send them over a balanced tree of the 20
owners (five hops) to the row owner. A group sums distinct columns of one 2^56-element
block with |coefficient| <= 255 and values < r, so every group sum is below
2^56 x 2^8 x 2^255 = 2^319 whatever the row's weight (no per-row weight bound is
used). **Finishing (R4).** The row owner adds its at most 28 group sums
(14 blocks x 2 signs) into a positive and a negative six-limb accumulator
(< 28 x 2^319 < 2^324; at most 12 instructions each), forms
v = D + pos + (K - neg) with the row's nine-limb dense accumulator D (21
Schirokauer terms and, in the Horner steps, n more, each < 2^510: D < 2^515) and
K a precomputed multiple of r in [2^324, 2^325), so 0 <= v < 2^518, and reduces v
mod r once: v = v1 2^320 + v0 (v1 < 2^198, v0 < 2^320),
v mod r = MM(v1, R^2 mod r) + MM(v0, R mod r) with two five-limb Montgomery
products MM (R = 2^320, measured MONT5 = 367.76 each; both products are below
r R) and at most three conditional subtractions: C_RED = 2 x 367.76 + 40 = 775.52.
Per row: headers <= 100 + 28 x 12 + 775.52 = 1,211.52 <= 13 x 100.7 = 1,309.1, the
finishing charge per row at the mean weight (the 13 per entry sums to 13 W N over
the N rows, 8.4% above the explicit budget; audit check `la_finish_budget`). Projection inner products
accumulate up to N < 2^60 terms < 2^510 in nine limbs (< 2^570) and are reduced
once each (C_RED, charged in the model); each matrix-generator multiplication is
charged 97 + C_RED. Transfers (R7). Each row gathers at most 28 group sums (14 blocks x 2 signs)
of 5 limbs over a static tree of depth <= 5. A hop forwards 5 limbs: 3 loads
+ 3 stores + 3 address + 2 branch/misc + 1 spare = 12; a group costs 5x12 + 8
endpoint setup = 68; a row <= 28x68 + 20 dispatch = 1,924 <= 20 x 100.7 =
2,014 (4.7% margin). Price per entry: 33 (`la_mac.s`) + block
selection 4 + sign 3 + transfers 20 + finishing and row headers 13 = 73,
charged as **73**; the
full-field multiply-accumulate inner loop is 97 measured, charged **97**
(replacing the cushioned 100).

## Aggregate

Peak storage is at most 2^105 bytes (the largest descent node). The main phase
stays below 2^81 bytes: active chunk sum arrays and record pages, the stored
relations (2^64.16 x <= 256 bytes = 2^72.2), the paged sparse matrix (2^59.43
rows x 100.7 entries x 9 bytes < 2^69.3), 20 vector blocks of 2^61 bytes for
<= 1024 block Wiedemann vectors (2^75.3), canonical-list staging 2^78 bytes.

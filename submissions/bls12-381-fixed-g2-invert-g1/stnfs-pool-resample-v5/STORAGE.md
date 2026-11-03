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

Price per retained record, measured end to end (recover-P thread, R3): scatter
at emission (7, above), page hand-off, the owner's copy of each page (16 bytes per
record, 3 instructions) and its page walk into apply; executed-instruction count
under qemu-aarch64 8.2.2 (gcc 13.3 -O3 -fno-tree-vectorize -mgeneral-regs-only),
differential 2^18 -> 2^19 records: 26.0009 per record minus the source loop 7 and
apply 9 = **10.0009** (37.2 per 2^16-record page), charged **10.01**.
Line-end partial pages cost <= about 50 per (line, chunk), negligible.
Plus **8 per cell** for row bookkeeping and sentinels; apply
**7.00** and recovery **7.003** per record are separate terms. Recovery re-reads every record of every chunk (no gating).
**Primes in [2^63, 2^64)** (1.4% of records; setup6 is measured up to 2^63 - 25,
the record path by hits per coset; a residue sum may exceed 64 bits): record path
and setup6 charged twice (`WIDE_FACTOR`). The measured
record-path cost (REC_TAB) is charged in full on top, so any buffer write
inside it is counted twice.

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

N < 2^58.998 < 10 x 2^56 elements of four 64-bit limbs; element j lives in
block j >> 56 (ten active owners; twenty 2^61-byte blocks reserved). Rows store
entries grouped by (block, sign); each owner sends a (pos, neg) partial sum up a
fixed tree to the row owner, which adds the 21 dense Schirokauer terms and
reduces mod r. Kernels (row tail, group loop; listings `la_mac.s`, `bw_dot.s`
in `kernels_asm.b85`) are counted by [10].

**Entries.** A +-1 entry adds its element without a multiplication (unit loop
11 instructions); other entries use the general loop (33). Per entry
4 + 3 + 11 f + 33 (1 - f) at unit share f = 0.88 (H4u) = 20.64, charged
**20.7**; total entries <= N x 201 (merge target 200).

**Row tail.** Nine sends (`_send_pair` 15 + 18 transfer) and nine combinations
(`_combine` 62 + 2), `_finish` 92 + `_fold9` 180 + `_reduce6` 183 + 2, headers
100, routing 40, group overhead 168: **1,638**, charged **1,640** (`C_ROW`). Sums
stay below 2^323 (six limbs). One extra 9-limb reduction (`_reduce9`) is 385,
charged **420**; dense terms 97 (`bw_full`).

With blocks (m, n) = (19, 4): linear algebra **<= 2^117.5139**, generator
**<= 2^106.9115** (proof.md section 4). H4u: unit share 0.917..0.948 at weight
200, 2^12..2^15 (offline filtsim counts), rising with scale; a shortfall changes cost only.

## Aggregate

Peak storage is at most 2^105 bytes (the largest descent node). The main phase
stays below 2^81 bytes: active chunk sum arrays and record pages, the stored
relations (2^64.17 x <= 256 bytes = 2^72.2), the paged sparse matrix (2^59.00
rows x 201 entries x 9 bytes < 2^69.9), 20 vector blocks of 2^61 bytes for
<= 1024 block Wiedemann vectors (2^75.3), canonical-list staging 2^78 bytes.

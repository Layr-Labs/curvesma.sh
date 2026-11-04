# Storage and addressing

Every object is bounded and addressed by (owner, local index);
`store16.py` is the executable reference; [11] runs its checks.

## Relation collection

A segment is one line b5 = k of a slice; a line has at most 2^57.508 cells
(`store16.line_cells_bound`), so a cell index fits 58 bits. Rows are packed
suffix-lexicographically (cell = row base + b0), 16-byte row descriptors in
2^16-entry pages.

The record walk writes each record (cell, ideal: uint64) to
the current 1 MiB page (2^16 records) of its chunk, chunk = cell >> 30 (chunk id
< 2^29), through a 2^29-entry tail table (4 GiB): 7 instructions hot path. A full page goes to the chunk's two-level page directory
(128-bit ordinals). Records stay in emission order (`store16.scatter` checked against `store16.route`), no
line buffer, no radix passes. Each owner applies its records into its
2^30-entry uint32 sum array (4 GiB), scans against thresholds,
recovers survivors. Open pages hold at most 2^29 x 1 MiB = 2^49
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

Price per retained record, end to end: scatter
(7, above), page hand-off, owner page copy (16 bytes per
record, 3 instructions), page walk into apply; differential 2^18 -> 2^19 records: 26.0009 per record minus source loop 7 and
apply 9 = **10.0009**, charged **10.01**.
Line-end partial pages <= ~50 per (line, chunk), negligible.
Plus **8 per cell** (bookkeeping, sentinels); apply
**7.00**, recovery **7.003** per record separate.
**Primes in [2^63, 2^64)** (1.4% of records; setup6 measured to 2^63 - 25,
record path by hits/coset; residue sums may exceed 64 bits): record path
and setup6 charged twice (`WIDE_FACTOR`). The
record-path cost (REC_TAB) is charged in full on top, so buffer writes
inside it count twice.

## Descent nodes

A node's sieve emits records (u64 owner, u64 offset, u16 16 log2 p; separate f/g streams);
the owner adds into two 2^40-entry uint16 score arrays (scores below 2^16 everywhere). All owners hold
their arrays at once: half count H uses 4H bytes, at R = 384
(H < 2^102.44) 2^104.44 bytes over < 2^62.44 owners, offsets < 2^40. The kernel
refuses half count >= 2^103; [7] audits every sieved node
(largest score 32,709); store16.check round-trips the split at
boundaries. Records stream;
transfer at 64-bit owner ids:
four 16-bit radix passes (64) + append 6 + flush 4 + directory 2 + dispatch 8
= **84** per record, inside SIEVE_HIT. Enumeration on 2^40
workers holding exact ball index tables (12 x (R^2 + 1) x (2R + 1) x 16
bytes = 2^34.34 at R = 384, 2^74.4 in all). Nodes run one at a time:
peak 2^104.44 + 2^74.4 < 2^105 bytes.

## Linear algebra

N < 2^59.11 < 10 x 2^56 elements of four 64-bit limbs; element j lives in
block j >> 56 (ten active owners; twenty 2^61-byte blocks reserved). Rows store
entries grouped by (block, sign); each owner sends a (pos, neg) partial sum up a
fixed tree to the row owner, which adds the 21 dense Schirokauer terms and
reduces mod r. Kernels (`row_tail.c`, group loop; listings `la_mac.s`, `bw_dot.s`
in `kernels_asm.b85`) are counted by [10].

**Entries.** A +-1 entry adds its element without a multiplication (unit loop
11 instructions); other entries use the general loop (33). Per entry
4 + 3 + 11 f + 33 (1 - f) at unit share f = 0.88 (H4u) = 20.64, charged
**20.7**; total entries <= N x 201.

**Row tail.** Nine sends (`_send_pair` 15 + 18 transfer) and nine combinations
(`_combine` 62 + 2), `_finish` 92 + `_fold9` 180 + `_reduce6` 183 + 2, headers
100, routing 40, group overhead 168: **1,638**, charged **1,640** (`C_ROW`). Sums
stay below 2^323 (six limbs). One extra 9-limb reduction (`_reduce9`) is 385,
charged **420**; dense terms 97 (`bw_full`).

With blocks (m, n) = (19, 4): linear algebra **<= 2^117.7346**, generator
**<= 2^107.0897** (proof.md section 4). H4u: unit share 0.921..0.951 at
weight 100 and 0.915..0.953 at weight 200,
`unit_entries.py`; a shortfall changes cost only.

## Aggregate

Peak storage is at most 2^105 bytes (the largest descent node). The main phase
stays below 2^81 bytes: active chunk sum arrays and record pages, the stored
relations (2^63.69 x <= 256 bytes = 2^71.7), the paged sparse matrix (2^59.11
rows x 201 entries x 9 bytes < 2^70.0), 20 vector blocks of 2^61 bytes for
<= 1024 block Wiedemann vectors (2^75.3), canonical-list staging 2^78 bytes.

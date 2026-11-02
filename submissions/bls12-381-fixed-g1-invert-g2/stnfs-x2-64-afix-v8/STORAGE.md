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
apply 9 = **10.0009** (37.2 per 2^16-record page), charged **10.01** (was 21).
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

N = 2^59.1758 elements of four 64-bit limbs (32 bytes); element j lives in block
j >> 56 at local index j mod 2^56 (blocks of 2^61 bytes, 10 needed, 20
allocated). Rows store entries grouped by (block, sign) with 64-bit column
indices and coefficients as an 8-bit magnitude. Block owners accumulate positive
and negative partial sums and send them up a balanced tree of the 20 owners to
the row owner, which reduces mod r. Base charge (e3940d11): 33 (`mac_run`) +
block selection 4 + sign 3 + transfers 32 + finishing and row headers 16 = 88 per
entry; dense terms 97 (`bw_full`, a 4 x 4-limb product added into a 9-limb
integer accumulator) with no separate reduction charge. This revision prices
every reduction mod r explicitly and charges the row tail per row.

Kernels are clang 18 output (`--target=arm64-apple-macos14 -O3 -fno-vectorize
-fno-slp-vectorize -mgeneral-regs-only`; the row tail also
`-fno-stack-protector`) for `row_tail.c` and the group loop; the stripped listings
(sections `la_mac.s`, `bw_dot.s`) are compressed in `kernels_asm.b85` (base85 of
lzma), and [10] decodes and counts them.

**Entries (per entry, total entries = N x mean weight <= N x 101).** A +-1
entry adds its 4-limb element into the positive or negative sum of its group
with no multiplication; a group stores a unit list (index only) and a general
list (index, magnitude 2..255). `_group_run`: unit loop `LBB9_1` **11**
instructions (ldr, add, 2 ldp, adds, 3 adcs, cinc, subs, b.ne), general loop
`LBB9_3` **33**. Per entry 4 + 3 + 11 f + 33 (1 - f) at aggregate unit fraction
f = 0.88 (H4u) is 20.64, charged **20.7** (`C_ENTRY`). The cost is linear in the
total entry count, so only the aggregate fraction and mean weight enter; no
per-row weight bound is used.

**Bounds for every row.** A group holds distinct columns of one 2^56-element
block, so a group sum is < 2^56 x 255 x r < 2^319 (five limbs; the unit loop's
`cinc` limb overflows only after 2^64 entries). A row has fewer than N + 21 <
2^60 entries, so every combined sum is < 2^60 x 255 x r < 2^323: messages and
`_combine` use six limbs. The 21 dense Schirokauer terms of a row go into the row
owner's 9-limb accumulator (< 21 x 2^510).

**Row tail (per row and product).** Each of the 19 non-root owners sends one
(pos, neg) message of 12 words (rows run in the same order everywhere: no
header). Loop-free kernels (every branch forward): `_send_pair` **15**,
`_combine` **62**; `_finish` (**92**, calling `_fold9` **180** and `_reduce6`
**183**) computes (pos - neg + dense) mod r: fold9 maps the 9-limb accumulator to
six limbs, x0..3 + sum_{k=4..8} x_k (2^64k mod r) < 2^323; v = pos + (K - neg) +
that, K = k r in [2^323, 2^324), so v < 2^325; reduce6 takes
qh = floor(floor(v / 2^192) MU / 2^256), MU = floor(2^448 / r), which is q or
q - 1, subtracts qh r and r once more if needed. Tested on 30,000 random and
extreme inputs against exact integers.
Per row: 19 x (15 + 2 top-limb stores + flush 4 + directory 2 + dispatch 8 +
call 2) + 19 x (62 + 2) + (92 + 180 + 183 + 2) + headers 10 x 20 + routing y_i
to its block owner 40 + group overhead 6 x 28 = **2,708**, charged **2,710**
(`C_ROW`).

**Other reductions.** One 9-limb accumulator alone (`_reduce9` = 20 + fold9 +
reduce6 + call 2 = **385**, charged **420** with the final modular add,
`C_RED9`): each Horner element per step (n dense terms then one reduction), each
of the m n projection entries per step (N dense terms, < 2^570, then one
reduction), and every generator multiplication (charged 97 + 420).

Block Wiedemann (proof.md section 4) is cheapest at **(m, n) = (19, 4)**;
linear algebra **<= 2^117.7103**, generator **<= 2^107.1997**. [10] counts every
function, checks they are loop-free and call only fold9/reduce6, recomputes 2,708
and 385, and checks `C_ENTRY >= 4 + 3 + 0.88 x 11 + 0.12 x 33`.

**H4u (heuristic): at least 88% of the production filtered matrix row weight is
+-1 entries.** `unit_entries.py` reruns `filtsim.py` with counters on a copy of
`filtsim.merge`: 0.917..0.946 at scales 2^12..2^15 (proof.md section 9), rising
with scale. A shortfall changes cost only (22 per missing unit entry).

## Aggregate

Peak storage is at most 2^105 bytes (the largest descent node). The main phase
stays below 2^81 bytes: active chunk sum arrays and record pages, the stored
relations (2^64.14 x <= 256 bytes = 2^72.2), the paged sparse matrix (2^59.18
rows x 101 entries x 9 bytes < 2^69.0), 20 vector blocks of 2^61 bytes for
<= 1024 block Wiedemann vectors (2^75.3), canonical-list staging 2^78 bytes.

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

**Route.** The chunk owner copies each full page into its receive area (kept for
recovery) and walks its directory; copy loop, 4 records per iteration (clang 18 -O3 aarch64):

```
recv_loop:
	ldp	x12, x11, [x10, #-32]
	add	x8, x8, #4
	cmp	x8, x2
	stp	x12, x11, [x9, #-32]
	ldp	x11, x13, [x10, #-16]
	stp	x11, x13, [x9, #-16]
	ldp	x13, x12, [x10]
	stp	x13, x12, [x9]
	ldp	x12, x11, [x10, #16]
	add	x10, x10, #64
	stp	x12, x11, [x9, #16]
	add	x9, x9, #64
	b.lo	recv_loop
recv_end:
```

Page hand-off (`emit_page`, 14 instructions) plus call, directory entry and the
owner's per-page loop are below 64 instructions per 2^16-record page (< 0.001
per record). Price per retained record: 7 + 13/4 + 0.01 = 10.26, charged
**10.3** (was 7 + flush 4 + directory 2 + dispatch 8 = 21). recover-P measured
the same route end to end by executed-instruction counts (gcc 13 -O3 under
qemu-aarch64, n vs 2n records): 10.0009 per record beyond source and apply.
**Recovery.** Each chunk owner re-reads every record of its chunk (7.003 per
record, no survivor gating) to rebuild survivors' factorizations.
Plus **8 per cell** for row bookkeeping; apply **7.00** and recovery **7.003** per
record are separate terms. The measured record-path cost (REC_TAB) is charged
in full on top, so any buffer write inside it is counted twice.

## Survivors: trial division below 2^20 and the band strip

Line-sieved primes in [2^11, 2^48) only add scores; no record is emitted for
them and nothing is re-stepped. Each survivor gets
exact norms, then trial division by every prime below 2^20 on both sides (512
instructions per prime per side: fixed 16-limb schoolbook remainder, 1024 bits,
above the 850-bit model-grid edge; at most 94,933 primes by Rosser-Schoenfeld,
`tdiv_prime_bound` checks pi(2^20) = 82,025), then two independent rounds per
side of the 48-bit `ecm16` schedule at 99.9% per prime on the remaining
cofactor, priced at the cell's exact norm rounded up to whole 64-bit limbs.
A prime in (2^20, 2^48] dividing the norm is found unless both rounds miss it
(<= 10^-6 under H1's random group orders); there are at most bits/20 such
primes, so each cell's yield is multiplied by (1 - k_f 10^-6)(1 - k_g 10^-6).
The walk records then supply primes in [2^48, 2^64) and the 69-bit schedule the
rest, so every relation is completely factored; a missed band prime only loses
that relation and is counted in the yield.

## Descent nodes

A node's sieve streams records (u64 owner = index >> 40, u64 offset, u16
16 log2 p; f and g streams) to owners, which add into two 2^40-entry uint16
score arrays (below 2^16 at every node). All owners hold their arrays at once:
half count H uses 4H bytes, at R = 384 (H < 2^102.44) 2^104.44 bytes over
< 2^62.44 owners; the kernel refuses half counts >= 2^103, [7] audits every
sieved node (largest score 32,709) and store16.check round-trips the split.
Transfer at 64-bit owner ids: 64 + 6 + 4 + 2 + 8 = **84**/record, inside
SIEVE_HIT. Enumeration runs on 2^40 workers with exact ball index tables
(2^74.4 bytes in all); nodes run one at a time: peak < 2^105 bytes.

## Linear algebra

N = 2^58.6340 elements of four 64-bit limbs; element j lives in block j >> 56
(blocks of 2^61 bytes). N + 21 < 7 x 2^56, so exactly 7 blocks hold elements
(`la_owners` checks ceil((N + 21)/2^56) <= 7); 20 blocks stay allocated. Row i
is owned by the owner of block i >> 56, one of the 7 active owners. Rows store
entries grouped by (block, sign) with 64-bit column indices and 8-bit
magnitudes; each active owner forms its (pos, neg) partial sums and the 6 other
active owners send theirs to the row owner, which reduces mod r. Kernels are
clang 18 output (`-O3 -fno-vectorize -fno-slp-vectorize -mgeneral-regs-only`)
for `row_tail.c` and the group loop; source and listings (`row_tail.c`,
`la_mac.s`, `bw_dot.s`) are in `kernels_asm.b85` and [10] counts them.

**Entries (total N x 304).** A +-1 entry adds its element with no
multiplication: unit loop `LBB9_1` **11** instructions, general loop `LBB9_3`
**33**. Per entry 4 + 3 + 11 f + 33 (1 - f) = 20.64 at unit share f = 0.88
(H4u), charged **20.7** (`C_ENTRY`). A group sum is < 2^56 x 255 x r < 2^319;
a row has < 2^60 entries, so every combined sum is < 2^323 (six limbs). The 21
dense Schirokauer terms go into the row owner's 9-limb accumulator.

**Row tail (per row and product).** 6 sends x (`_send_pair` 15 + 2 top-limb
stores + flush 4 + directory 2 + dispatch 8 + call 2) + 6 combinations x
(`_combine` 62 + 2) + `_finish` 92 + `_fold9` 180 + `_reduce6` 183 + 2 +
headers 10 x 7 + routing y_i to its block owner 40 + group overhead 6 x 28 =
**1,317**, charged **1,320** (`C_ROW`). `_finish` computes (pos - neg + dense)
mod r: fold9 maps the accumulator to six limbs (< 2^323), v = pos + (K - neg)
+ that < 2^325 with K = k r in [2^323, 2^324), reduce6 subtracts qh r and r
once more if needed (tested on 30,000 random and extreme inputs). One 9-limb
reduction alone (`_reduce9`) is **385**, charged **420** (`C_RED9`); dense
terms 97 (`bw_full`). With (m, n) = (22, 4): linear algebra **<= 2^117.034782**,
generator **<= 2^106.802009**. [10] checks every function is loop-free, calls
only fold9/reduce6, and recomputes 1,317 and 385.

**H4u (heuristic): >= 88% of the production filtered row weight is +-1
entries.** `audit_g1.py` {"unit_entries": [k, seed]} measures it (0.889+ from
2^15 up at weight 300). A shortfall changes cost only (22 per missing unit).

## Aggregate

Peak storage is at most 2^105 bytes (the largest descent node). The main phase
stays below 2^81 bytes: chunk arrays and pages, relations (2^64.16 x 256 bytes
= 2^72.2), the matrix (2^58.64 x 304 x 9 bytes < 2^70.1), 20 vector blocks of
2^61 bytes (2^75.3) and canonical-list staging (2^78).

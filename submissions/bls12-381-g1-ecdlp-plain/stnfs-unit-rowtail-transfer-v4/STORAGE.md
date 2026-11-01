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

Price per retained record: scatter at emission 7 (replacing append 6 and radix
routing 16), page flush 4, directory 2, sentinel and dispatch 8: **21**, plus
**8 per cell** for row bookkeeping; apply **7.00** and recovery **7.003** per
record are separate terms. The measured record-path cost (REC_TAB) is charged
in full on top, so any buffer write inside it is counted twice.

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
allocated). A row is stored per (block, sign) group: a unit list (64-bit column
indices of the +-1 entries) and a general list (index and 8-bit magnitude
2..255). Each block owner forms the positive and negative five-limb partial
sums of the row (101 x 255 x r < 2^270). All kernels are clang 18 output
(`--target=arm64-apple-macos14 -O3 -fno-vectorize -fno-slp-vectorize
-mgeneral-regs-only`) in `la_mac.s`; [10] and audit check `la_entry_price`
count them.

**Groups (R3).** `_group_run` runs one group: unit loop `LBB2_1` **11**
instructions per entry (ldr, add, 2 ldp, adds, 3 adcs, cinc, subs, b.ne: a +-1
entry is an addition, no multiply), general loop `LBB2_3` **33** (identical to
`mac_run`). Outside the loops it has 4 more instructions than `mac_run`; with two
header words, at most 6 per group and 28 groups per row (14 blocks x 2 signs):
168 per row.

**Row tail (R3, recover-Q 06aac4a1).** Partial sums are combined at each node of
a fixed tree of the 20 owners, so each non-root owner sends one 10-word (pos,
neg) message per row, no header (rows run in the same order everywhere); sums
stay below 2^270. Loop-free kernels: `_send_pair` **22**, `_combine` **51**,
`_finish` **119** (v = pos + (K - neg), K = k r in [2^270, 2^271); Barrett with
MU = floor(2^320 / r) on floor(v / 2^192), then one conditional subtraction;
tested here against exact integers on 20,000 random and boundary pairs, 0
errors). Per row and product, all 20 owners charged: 19 x (22 + flush 4 +
directory 2 + dispatch 8 + call 2) + 19 x (51 + 2) + (119 + 2) + headers 10 x 20
+ routing y_i 40 = **2,090** (was 32 + 16 = 48 per entry, 4,848 per row).

**Entry price.** At unit fraction f: 4 (block selection) + 3 (sign) + 2090/101 +
11 f + 33 (1 - f) + 168/101 = **43.00 at f = 0.88**, charged **C_ENTRY = 44**
(was 88). Block Wiedemann is then minimised at **(m, n) = (18, 4)** (failure
4N/r < 2^-192 for any m, n): linear algebra **<= 2^118.7428**, generator
2^105.5088. The full-field loop (`bw_dot.s`) stays 97, charged 97.

**H4u (heuristic).** At least 88% of the filtered matrix's row weight is +-1
entries. `audit_g1.py` with input `{"unit_entries": [k, seed]}` reruns
`filtsim.py` at scale 2^k (signed system, mean norms 747.5/498.9, B = 2^69,
excess 1.3, merge to weight 100) with an unchanged copy of `filtsim.merge` that
also tracks each row's integer combination of original relations. Explicit
entries count with their exact coefficients. The implicit mid-range weight
(`r[-1]`, one w_mid per original relation) counts as unit only for relations
with multiplier +-1, minus 2 C(k, 2) S2 + S1 per row (S2 = sum p_q^2 over
mid-range primes: collisions; S1 = 2 sum p_q / q: higher powers). Measured:
2^12 seeds 1-3: 0.897, 0.917, 0.923; 2^14: 0.933, 0.933, 0.936; 2^15: 0.940,
0.937, 0.939; 2^17: 0.9467, 0.9464, 0.9466. The share rises with scale, and the
charge covers f >= 0.835. Like N/U (H3) this is an extrapolation; a shortfall
changes cost only (22 instructions per missing unit entry), never correctness.

## Aggregate

Peak storage is at most 2^105 bytes (the largest descent node). The main phase
stays below 2^81 bytes: active chunk sum arrays and record pages, the stored
relations (2^64.16 x <= 256 bytes = 2^72.2), the paged sparse matrix (2^59.77
rows x 101 entries x <= 9 bytes < 2^69.6), 20 vector blocks of 2^61 bytes for
<= 1024 block Wiedemann vectors (2^75.3), canonical-list staging 2^78 bytes.

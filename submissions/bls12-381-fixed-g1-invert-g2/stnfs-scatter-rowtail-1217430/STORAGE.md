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
allocated). Rows store entries grouped by block with 64-bit column indices and
coefficients as a sign bit and an 8-bit magnitude. Block owners form positive
and negative five-limb partial sums (101 entries x 255 x r < 2^270) and send
them over a balanced tree of the 20 owners to the row owner, which reduces mod r.
Base charge (e3940d11): 33 (`la_mac.s` multiply-accumulate) + block selection 4 +
sign 3 + transfers 32 + finishing and row headers 16 = 88 per entry. The
full-field multiply-accumulate inner loop is 97 measured, charged **97**.

Kernels below are clang 18 output (`--target=arm64-apple-macos14 -O3
-fno-vectorize -fno-slp-vectorize -mgeneral-regs-only`) in `la_mac.s`.

**Unit entries.** After filtering almost all coefficients are +-1; such an entry
needs no multiplication: its 4-limb element is added into the positive or
negative partial sum of its (block, sign) group. Each group stores a unit list
(index only) and a general list (index, magnitude 2..255). `_group_run`: unit
loop `LBB9_1` **11** instructions per entry (ldr, add, 2 ldp, adds, 3 adcs,
cinc, subs, b.ne), general loop `LBB9_3` **33** (as `mac_run`). Outside the
loops 13 instructions against 9 for `mac_run`; with two more header words at
most 6 extra per group and at most 28 groups per row (14 blocks x 2 signs): at
most 168 per 101-entry row, charged 2 per entry.

**Row tail (this revision).** The base priced transfers as if every partial sum
travelled the tree separately (32 per entry) and the finish at 16. Here sums are
combined at each tree node: each of the 19 non-root owners sends one (pos, neg)
message of 10 words per row; rows run in the same order on every owner, so a
message has no header. Combined sums stay below 101 x 255 x r < 2^270.
Loop-free kernels (every branch forward): `_send_pair` copies the pair to the
parent's buffer, **13** instructions; `_combine` adds an incoming pair into the
local one, **51**; `_finish` at the row owner forms v = pos + (K - neg),
K = k r in [2^270, 2^271), qh = floor(floor(v / 2^192) MU / 2^128) with
MU = floor(2^320 / r) (qh is q or q - 1), subtracts qh r and then r if needed:
**119**; tested against exact integers on 20,000 random and boundary pairs.
Per row and product, charging all 20 allocated owners: 19 x (13 + flush 4 +
directory 2 + dispatch 8 + call 2) + 19 x (51 + 2) + (119 + 2) + row headers
10 x 20 + routing y_i to its block owner 40 = **1,919**, 19.00 per entry,
charged **20** (base 32 + 16 = 48).

Per-entry charge at unit fraction f: `4 + 3 + 20 + 11 f + 33 (1 - f) + 2`; at
**f = 0.88** 42.64, charged **43** (`C_ENTRY`; base 88, 04073dd1 71). Block
Wiedemann (proof.md section 4) is then cheapest at **(m, n) = (18, 4)**; the
failure bound 4N/r < 2^-192 does not depend on m, n. Linear algebra
**<= 2^118.7243**, generator **<= 2^105.5087**. `evidence16.py` [10] counts the
three functions and both loops in `la_mac.s`, checks the functions are
loop-free, recomputes 1,919 <= 20 x 101 and checks `C_ENTRY >= 4 + 3 + 20 +
0.88 x 11 + 0.12 x 33 + 2`.

**H4u (heuristic): at least 88% of the production filtered matrix row weight is
+-1 entries.** `unit_entries.py` reruns the package's filtering simulation
(`filtsim.py`: signed sigma-quotient system, mean norms 748.39/498.93 bits,
B = 2^69, excess 1.3, purge, Markowitz merge to mean weight 100, |coefficient|
<= 255) with an unchanged copy of `filtsim.merge` plus counters and measures the
unit share of row weight at mean weight 100 (explicit entries exact; filtsim's
implicit mid-range weight counts as unit unless it came from a pivot with
|m| > 1). Measured: 2^12 seeds 1-3 0.9330, 0.9350, 0.9173; 2^13 seeds 1-2
0.9338, 0.9405; 2^14 seeds 1-3 0.9452, 0.9376, 0.9452; 2^15 seed 1 0.9461. The
share rises with scale; 0.88 is at least 0.037 below every value. Like N/U (H3)
this is an extrapolation; a shortfall changes cost only (22 per missing unit
entry), never correctness.

## Aggregate

Peak storage is at most 2^105 bytes (the largest descent node). The main phase
stays below 2^81 bytes: active chunk sum arrays and record pages, the stored
relations (2^64.16 x <= 256 bytes = 2^72.2), the paged sparse matrix (2^59.77
rows x 101 entries x 9 bytes < 2^69.6), 20 vector blocks of 2^61 bytes for
<= 1024 block Wiedemann vectors (2^75.3), canonical-list staging 2^78 bytes.

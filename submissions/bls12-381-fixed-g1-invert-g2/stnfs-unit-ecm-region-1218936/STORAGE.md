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
records per cell inside the 427.4 instructions per cell of proof.md).

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
the base charge. This revision (unit-entry kernel section below) stores each (block, sign) group as a
unit list (index only, the 11-instruction add loop of `la_mac.s` (`_group_run`)) and a general
list (index and 8-bit magnitude, the 33-instruction loop), charging
4 + 3 + 32 + 16 + 0.88 x 11 + 0.12 x 33 + 2 = 70.64 as **71** per entry; the
full-field multiply-accumulate inner loop is 97 measured, charged **97**
(replacing the cushioned 100).

## Unit-entry linear-algebra kernel (this revision)

Base charge (STORAGE.md): 88 instructions per sparse entry = 33 (`la_mac.s` multiply-accumulate)
+ block selection 4 + sign 3 + transfers 32 + finishing/row headers 16.

**Change.** After filtering almost all coefficients are +-1. A +-1 entry needs no multiplication: the
4-limb element is added into the positive or negative five-limb partial sum of its (block, sign) group,
which the base layout already keeps separately. Each group stores a unit list (index only) and a general
list (index and 8-bit magnitude 2..255). `la_mac.s` (`_group_run`) is clang 18 output (`--target=arm64-apple-macos14
-O3 -fno-vectorize -fno-slp-vectorize -mgeneral-regs-only`) for one group: unit loop `LBB9_1` is
**11** instructions per entry (ldr, add, 2 ldp, adds, 3 adcs, cinc, subs, b.ne); general loop
`LBB9_3` is **33**, identical to `la_mac.s`. The partial sums stay in registers across both loops.
Outside the loops the function has 13 instructions against 9 for `mac_run`; with two more header words
per group, the extra fixed cost is at most 6 per group. Columns occupy 14 blocks (STORAGE.md), so a row
has at most 28 groups: at most 168 extra instructions per 101-entry row, below 1.67 per entry, charged 2.

Per-entry charge at charged unit fraction f: `4 + 3 + 32 + 16 + 11 f + 33 (1 - f) + 2 = 90 - 22 f`;
at **f = 0.88** this is 70.64, charged **71** (`C_ENTRY`). With the cheaper product the Block Wiedemann
cost of proof.md section 4 is minimised at **(m, n) = (22, 4)** instead of (26, 5); the failure bound
4N/r < 2^-192 does not depend on m, n. Linear algebra 2^119.36845 -> **2^119.1608**, generator
2^106.4172 -> **2^106.1537**.

**H4u (heuristic): at least 88% of the production filtered matrix row weight is +-1 entries.**
`unit_entries.py` reruns the package's filtering simulation (`filtsim.py`: signed sigma-quotient system,
mean norms 748.39/498.93 bits, B = 2^69, excess 1.3, purge, Markowitz merge to mean weight 100,
|coefficient| <= 255) with an unchanged copy of `filtsim.merge` plus counters, and measures the unit share
of row weight when the mean weight reaches 100. Explicit entries use exact coefficients; filtsim's
implicit mid-range weight (`r[-1]`) counts as unit unless it entered a row from a pivot with |m| > 1.

Measured unit fractions: 2^12 seeds 1-3: 0.9330, 0.9350, 0.9173; 2^13 seeds 1-2: 0.9338, 0.9405; 2^14 seeds 1-3: 0.9452, 0.9376, 0.9452; 2^15 seed 1: 0.9461.

The fraction rises with scale; the charge 0.88 is at least 0.037 below every value. Like N/U (H3) this is
an extrapolation to production. A shortfall changes cost only (22 instructions per missing unit entry),
never correctness. `evidence16.py` checks both loop lengths and `C_ENTRY >= 88 - 33 + 0.88 x 11 +
0.12 x 33 + 2`.

## Aggregate

Peak storage is at most 2^105 bytes (the largest descent node). The main phase
stays below 2^81 bytes: active chunk sum arrays and record pages, the stored
relations (2^64.16 x <= 256 bytes = 2^72.2), the paged sparse matrix (2^59.77
rows x 101 entries x 9 bytes < 2^69.6), 20 vector blocks of 2^61 bytes for
<= 1024 block Wiedemann vectors (2^75.3), canonical-list staging 2^78 bytes.

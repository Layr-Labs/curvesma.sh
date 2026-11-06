# Storage and addressing

Every object is bounded and addressed by (owner, local index); `store16.py`
is the reference, [11] runs its checks.

## Relation collection

A segment is one line b5 = k of a slice; a line has at most 2^57.508 cells
(`store16.line_cells_bound` at the full slice radius A^2, so it also bounds the region's (17/20) A^2 slices), so a cell index fits 58 bits. Rows are packed
suffix-lexicographically (cell = row base + b0) with 16-byte row descriptors in
2^16-entry pages.

The record walk writes each record (cell: uint64, ideal: uint64) directly to
the current 1 MiB page (2^16 records) of its chunk, chunk = cell >> 30 (chunk id
< 2^29), through a 2^29-entry tail table (4 GiB): 7 instructions on the hot path (shift, tail load, add, paired store, page-end test and branch,
tail store). Full pages go to the chunk's two-level page directory
(128-bit record ordinals). A chunk's records stay in emission order: the stable
chunk order of the source's line buffer plus two-pass 16-bit radix sort
(`store16.scatter` is checked against `store16.route`), without either. Each chunk owner applies its records to its
2^30-entry uint32 sum array (4 GiB), scans it against the thresholds and
recovers survivors' records. Open pages hold at most 2^29 x 1 MiB = 2^49
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

Price per retained record, measured end to end: scatter
at emission (7, above), page hand-off, the owner's copy of each page (16 bytes per
record, 3 instructions) and its page walk into apply; executed-instruction count
under qemu-aarch64 8.2.2 (gcc 13.3 -O3 -fno-tree-vectorize -mgeneral-regs-only),
differential 2^18 -> 2^19 records: 26.0009 per record minus the source loop 7 and
apply 9 = **10.0009** (37.2 per 2^16-record page), charged **10.01**.
Plus **8 per cell** for row bookkeeping and sentinels; apply
**5.76** and recovery **5.76** per record are separate terms (`slice12oh.c`); scan2 zeroes both sum arrays as it
reads them, and each survivor's SENTINEL is cleared after recovery (8 instructions per survivor, charged in the scan row's
per-leaf part). Recovery re-reads all records of every chunk (no gating).
**Primes in [2^63, 2^64)** of contaminated runs take the wide path of KERNELS.md, with the same record layout.
**Primes in [2^56, 2^64) on a run** (KERNELS.md, rec7): a run's slices live together (one walk per (run, ideal) emits
records into every slice's chunks): at most 2101 a-values, 2^78.57 cells, 4 bytes of sums and 0.575 expected 16-byte
records per cell: below 2^82.3 bytes. The two 7-dim bases per ideal take < 2 KiB each: < 2^71.6 bytes.

## Descent nodes

A node's sieve emits records (u64 owner = index >> 40, u64 offset = index mod
2^40, u16 16 log2 p; separate f and g streams; dsieve16 rec_t) to the owner,
which adds into its two 2^40-entry uint16 score arrays (a score is at most
16 x the norm size, below 2^16 at every node). All owners hold
their arrays at once: a node of half count H uses 4H bytes, at R = 384
(H < 2^102.44) 2^104.44 bytes over < 2^62.44 owners, offsets < 2^40. The kernel
refuses any ball whose half count reaches 2^103; [7] audits every sieved node
(largest score 32,709) and store16.check round-trips the split at the
boundaries, including index 2^94. Records stream through, unstored;
the transfer is priced by the relation-collection rule at 64-bit owner ids:
four 16-bit radix passes (64) + append 6 + flush 4 + directory 2 + dispatch 8
= **84** per record, inside SIEVE_HIT (KERNELS.md). Enumeration runs on 2^40
workers holding the exact ball index tables (12 x (R^2 + 1) x (2R + 1) x 16
bytes = 2^34.34 bytes at R = 384, 2^74.4 in all). Nodes run singly, so
the descent's peak is 2^104.44 + 2^74.4 < 2^105 bytes.

## Linear algebra

N <= 2^58.5430 < 6 x 2^56 elements of four 64-bit limbs; element j lives in
block j >> 56 (six active owners). Rows store entries grouped by (block, sign);
each owner sends a (pos, neg) partial sum up a fixed tree to the row owner, which
adds the 21 dense Schirokauer terms and reduces mod r. Kernels (`row_tail.c`; listings `la_mac.s`, `bw_dot2.s`,
`krylov5.s`; all in `kernels_asm.b85`) are counted by [10].

**Entries.** Single-vector products (unused: Horner runs at the Krylov rate) would cost 4 + 3 + 11 f + 33 (1 - f) =
21.08 at f = 0.86 (unit loop 11, general 33), charged **21.1**; total entries <= N x 301 (merge target 300).

**Krylov passes.** The n = 5 sequences advance in lockstep: one pass over the matrix per Krylov step, element j
of the five vectors interleaved in one 256-byte slot. A +-1 group runs `unit_run5` (one index load and address per
entry, then per vector two `ldp`, the four-limb add and the carry count: 39 instructions per entry for the five
vectors, all 25 accumulators in registers); other groups run `gen_run1` once per vector (33 per vector-entry, stride
256). Block selection and sign (4 + 3) are paid once per +-1 entry and per
vector-entry otherwise; entries with |c| in {2, 3, 4} are stored as |c| unit entries (minimum shares 0.02664, 0.00868, 0.00417
over the 44 H4u runs at excess 1.08): 0.86 x 9.2 + (2 x 0.02664 + 3 x 0.00868 + 4 x 0.00417) x 9.2 + 0.10051 x 40 = 12.8156,
charged **12.82** (`C_ENTRY_K`; per-run k5 10.534..12.530 in the 44 runs, 10.05..10.73 at 2^13..2^15; after
recover-Q v7); row tail and dense terms stay per vector.

**Row tail (six owners).** Five sends (`_send_pair` 15 + 18 transfer) and five combinations
(`_combine` 62 + 2), `_finish` 92 + `_fold9` 180 + `_reduce6` 183 + 2, headers
10 per owner (60), routing 40, group overhead 168: **1,210**, charged **1,215**
(`C_ROW`, `owner10.row_tail`). Sums stay below 2^323 (six limbs). One extra
9-limb reduction (`_reduce9`) is 385, charged **420**; dense terms 97 (`bw_full`).

With blocks (m, n) = (72, 5) and block-diagonal 128-bit X (one X entry per element): linear algebra **<= 2^115.994158905**,
generator **<= 2^110.604** (proof.md section 4). Projection kernel `bw_dot2`: a 2-limb X entry times a 4-limb element into a 7-limb accumulator, fully unrolled. H4u
(unit share, `unit_entries.py`): proof.md section 9.

## Aggregate

Peak storage is at most 2^105 bytes (the largest descent node). The main phase
stays below 2^83 bytes: a run's sum arrays and records (< 2^82.3), rec7 bases, active record pages, the stored
relations (2^63.64 x <= 256 bytes = 2^71.7), the paged sparse matrix (2^58.69
rows x 301 entries x 9 bytes < 2^70.5), the interleaved Krylov block (2^58.69 x 256
bytes < 2^67) and Horner/projection vectors (< 2^76), canonical-list staging 2^78 bytes.

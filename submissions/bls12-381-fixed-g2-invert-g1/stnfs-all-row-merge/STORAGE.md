# Storage and addressing

Objects are bounded and addressed by (owner, local index); `store16.py` is the reference ([11]).

## Relation collection

A segment is one line b5 = k of a slice; a line has at most 2^57.508 cells
(`store16.line_cells_bound` at the full slice radius A^2, so it also bounds the region's (17/20) A^2 slices), so a cell index fits 58 bits. Rows are packed
suffix-lexicographically (cell = row base + b0) with 16-byte row descriptors in
2^16-entry pages.

The record walk writes each record (cell: uint64, ideal: uint64) directly to
the current 1 MiB page (2^16 records) of its chunk, chunk = cell >> 30 (chunk id
< 2^29), through a 2^29-entry tail table (4 GiB): 7 instructions on the hot path (shift, tail load, add, paired store, page-end test and branch,
tail store). Full pages go to the chunk's two-level page directory
(128-bit record ordinals). Records stay in emission order (`store16.scatter` is checked against `store16.route`, the line-buffer plus two-pass
radix sort it replaces). Each chunk owner applies its records to its
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

Price per retained record, end to end: scatter (7, above), page hand-off, the owner's page copy (16 bytes per record,
3 instructions) and its walk into apply; executed instructions under qemu-aarch64 8.2.2 (gcc 13.3 -O3 -fno-tree-vectorize -mgeneral-regs-only),
differential 2^18 -> 2^19 records: 26.0009 per record minus the source loop 7 and
apply 9 = **10.0009** (37.2 per 2^16-record page), charged **10.01**.
Plus **0.5 per cell** for every other per-cell step (measured 0.096: CELLRES, KERNELS.md); apply
**5.76** and recovery **5.76** per record are separate terms (`slice12oh.c`); scan2 zeroes the f sum array as it
reads them, and each survivor's SENTINEL is cleared after recovery (8 instructions per survivor, charged in the scan row's
per-leaf part).
**Primes in [2^63, 2^64)** of contaminated runs take the wide path of KERNELS.md, with the same record layout.
**Primes in [2^48, 2^64) on a run** (KERNELS.md, rec7): a run's slices live together (one walk per (run, ideal) emits
records into every slice's chunks): at most 2101 a-values, 2^78.57 cells, 4 bytes of sums and 0.575 expected 16-byte
records per cell: below 2^82.3 bytes. The two 7-dim bases per ideal take < 2 KiB each: < 2^71.6 bytes.

## Descent nodes

A node's sieve emits records (u64 owner = index >> 40, u64 offset = index mod 2^40, u16 16 log2 p; separate f and g
streams; dsieve16 rec_t) to the owner, which adds into its two 2^40-entry uint16 score arrays (a score is at most 16 x
the norm size, below 2^16 at every node). All owners hold their arrays at once: a node of half count H uses 4H bytes, at
R = 384 (H < 2^102.44) 2^104.44 bytes over < 2^62.44 owners, offsets < 2^40. The kernel refuses balls whose half count
reaches 2^103; [7] audits every sieved node (largest score 32,709) and store16.check round-trips the split at the
boundaries, including index 2^94. Records stream through unstored; transfer at 64-bit owner ids: four 16-bit radix passes (64) + append 6 + flush 4 + directory 2 + dispatch 8 = **84** per
record, inside SIEVE_HIT. 2^40 workers hold the exact ball index tables (12 x (R^2 + 1) x (2R + 1) x 16 bytes =
2^34.34 at R = 384, 2^74.4 in all); nodes run singly, so the descent's peak is 2^104.44 + 2^74.4 < 2^105 bytes.

## Linear algebra

N <= Mcap < 7 x 2^56 elements of four 64-bit limbs; element j lives in
block j >> 56 (seven active owners). Rows store entries grouped by (block, sign);
each owner sends a (pos, neg) partial sum up a fixed tree to the row owner, which
adds the 21 dense Schirokauer terms and reduces mod r. Kernels (`row_tail.c`, see KERNELS.md; listings `la_mac.s`, `bw_dot2.s`,
`krylov5.s` in `kernels_asm.b85`) are counted by [10].

**Entries.** Single-vector products (unused) would cost 4 + 3 + 11 f + 33 (1 - f) = 21.08 at f = 0.86, charged
**21.1**; entries <= N x 301.

**Krylov passes.** The n = 5 sequences advance in lockstep (one matrix pass per step; element j of the five
vectors in one 256-byte slot). A +-1 group runs `unit_run5` (one index load and address per
entry, then per vector two `ldp`, the four-limb add and the carry count: 39 instructions per entry for the five
vectors, all 25 accumulators in registers); other groups run `gen_run1` once per vector (33 per vector-entry, stride
256). Block selection and sign (4 + 3) are paid once per +-1 entry, per vector-entry otherwise; |c| in {2, 3, 4} is
stored as |c| unit entries: k5 charged **9.97** (`C_ENTRY_K`; the premerge H4u total-price functional, moments and abort are derived in proof.md, "Actual row ensemble, gates, square selection and production moments"); row tail and dense terms per vector.

**Row tail (seven owners).** Six sends (`_send_pair` 15 + 18 transfer) and six combinations
(`_combine` 62 + 2), `_finish` 92 + `_fold9` 180 + `_reduce6` 183 + 2, headers
10 per owner (70), routing 40, two sign groups at14 each per owner,196 total: **1,345**, charged **1,350**
(`C_ROW`, `owner10.row_tail`). Sums stay below 2^323 (six limbs). One extra
9-limb reduction (`_reduce9`) is 385, charged **420**; dense terms 97 (`bw_full`).

With blocks (m, n) = (72, 5) and block-diagonal 128-bit X (one X entry per element): linear algebra **<= 2^116.036935623**,
generator **<= 2^110.759382425** (proof.md, "Fixed-support projection and its derived coin error" and "Explicit work, including new map and administration charges"). Projection kernel `bw_dot2`: a 2-limb X entry times a 4-limb element into a 7-limb accumulator, fully unrolled. H4u
(weighted premerge coefficients, `unit_entries.py`): proof.md, "Actual row ensemble, gates, square selection and production moments".

## Aggregate

Peak storage is at most 2^105 bytes (the largest descent node). The main phase
stays below 2^86 bytes: a run's sum arrays and records (< 2^82.3), rec7 bases, active record pages, the stored
relations (hardcap2^65 x < 2^14 bytes: (a, b) and <= 832 + 576 (ideal, exponent) pairs; < 2^79), the paged sparse matrix (2^58.7
rows x 301 entries x 9 bytes < 2^70.5), the interleaved Krylov block (2^58.7 x 256
bytes < 2^67) and Horner/projection vectors (< 2^76), canonical-list staging 2^78 bytes, band P and trees < 2^52, P_g, a g wave and their trees < 2^74, band-split table < 2^50.2, a run's first-slice leaf partitions < 2^72 bytes.

The main-phase proof reproduces all probability and generator arguments here; no missing numbered proof section is imported. Protect the norm1399 idealI=(1399,y+40,x+1) and all21 dense coordinates. This changes the filtering ensemble, so its declared moments explicitly apply to this protected variant. The square uses coverage-preserving disjoint pair sums of ALL premerge rows; no row subset is selected. Full retained relations, at most2^65, remain for checks. The universal purge target is21+max(2,ceil(.0001U)); capM<=ceil(.0355U)+21. Final squareN<=M, and global sparse entry cost never increases under accepted pair sums, so the existingMcap ledger covers it. The new administration74 row charges row checks, random generation/order and endpoint handling; the universal map103 row covers cold maps2^65, descent maps2^65, one anchor and all fixed setup. Checkpoint storage has25 addedD-vectors, less than2^70 bytes; random ordering/indices and map data fit the conservative main2^86 bound. The descriptor and filter-work meters are included in the stated rows.

The larger row envelope requires seven vector owners; each keeps the original56-bit local index, so the256-byte interleaved vector address remains64-bit. All global row IDs fit59bits. Three-bit owner labels fitseven owners; no57-bit local offset is used. owner10.py recomputes the exact1345-instruction row protocol and charges1350. Group routing, sends, combinations, headers and both sign groups for the seventh owner are charged. The unchanged6-limb row accumulator bound still covers<7*2^56 inputs with signed8-bit coefficients.

The archived pairmerge.py is an OFFLINE diagnostic reference: its full candidate sums, Python unbounded-internal randrange, global all-coordinate assertions after each merge and histogram scans are not the production instruction path. Production uses the capped sampler and an incrementally maintained forbidden-pair dictionary described in proof.md. For each degree-two coordinate, store its two surviving row IDs and coefficient sum; it creates one forbidden key only when both IDs remain eligible and the sum is zero. Query this dictionary BEFORE forming a proposed sum. On acceptance, walk just the two disjoint rows' incident columns, update degrees/supports/keys, and perform the coefficient check. Global coverage assertions are redundant by induction; optional end checks cost one complete pass. This has the identical accepted-pair law and outputs until explicitly charged sampler/meter exhaustion, while avoidingO(NK) diagnostic assertions. All pair metadata/updates and failed proposals use the shared hard filter meter; construction-meter/representation aborts are part of the.01 event.

For the credited P7=48 extension, conservatively enlarge the main-phase storage envelope to2^86bytes, covering every additional co-resident run record, page directory and wider setup table. The prior detailed2^83 subtotal is not used as a new hard cap. Global memory remains below2^106 since even2^105+2^86+wrapper/audit state is below that bound; per-line/owner addressing stays unchanged.

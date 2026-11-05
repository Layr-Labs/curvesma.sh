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
apply 9 = **10.0009** (37.2 per 2^16-record page), charged **10.01** (was 21). The 9 is that harness's own apply
loop (its listing), subtracted to isolate routing; production applies records with slice12's `apply_batch`, a different
loop charged separately (`C_APPLY` 5.76: measured 5.751 in `slice12oh.c`, 7.000 unmodified).
rec7 emits each record into the record ring of its slice segment (`rec7.c`: segment = SG[t - t_lo][(|b5| - 1) >> g]), the rings
slice12's record walk fills, so pages, routing, apply and recovery per record are unchanged; a run walk keeps the rings of all
its slices (at most 1,024) open at once.
Line-end partial pages cost <= about 50 per (line, chunk), negligible.
Plus **8 per cell** (`CELL_ROUTE`, an allowance, not a measurement) for row bookkeeping and sentinels. It covers scan2's
SENTINEL clearing: the survivor pass resets each survivor cell before segment reuse, 8 instructions per survivor, and every
survivor is also charged >= 2^38 instructions of cofactorisation (proof.md section 3). Apply **5.76** and recovery **5.76** per
record are separate terms (`slice12oh.c`). Recovery re-reads every record of every chunk (no gating).
**Primes in [2^63, 2^64)** (1.4% of records): v9 runs them on the public slice12 record walk with `wide12.c`
(KERNELS.md): record path, empty cosets and apply at `REC_TAB_W` x 1.01, `setup6w` within the setup6 charge, and the
wide factor-base lattices as a ledger row (2^70.6). Records and routing are those of every other record prime.

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

N < 2^58.92 < 8 x 2^56 elements of four 64-bit limbs; element j lives in
block j >> 56 (eight active owners). Rows store entries grouped by (block, sign);
each owner sends a (pos, neg) partial sum up a fixed tree to the row owner, which
adds the 21 dense Schirokauer terms and reduces mod r. Kernels (`row_tail.c`, verbatim in `kernels_asm.b85` since v8,
group loop; listings `la_mac.s`, `bw_dot.s`, `krylov5.s` in `kernels_asm.b85`)
are counted by [10].

**Entries.** A +-1 entry adds its element without a multiplication (unit loop
11 instructions); other entries use the general loop (33). Per entry
4 + 3 + 11 f + 33 (1 - f) at unit share f = 0.87 (H4u) = 20.86, charged
**20.9** (single-vector products: Horner); total entries <= N x 301 (merge
target 300, v6).

**Krylov passes (v6).** The n = 5 Block Wiedemann sequences advance in lockstep,
so each Krylov step is one pass over the matrix for all five vectors: element j
of the five current vectors is stored interleaved in one 256-byte slot (limbs
4k..4k+3 hold vector k; three slots unused). A +-1 group runs `unit_run5`: one
index load and address computation per entry, then per vector two `ldp`, the
four-limb add and the carry count, 39 instructions per entry for the five
vectors (clang 18 keeps all 25 accumulators in registers, no spills). Other
groups run `gen_run1` once per vector (33 per vector-entry, the single-vector
general loop at stride 256). The listing (`krylov5.s` in `kernels_asm.b85`) is
clang output of: per entry x = v + 32 idx[t]; for k < 5 add x[4k..4k+3] with
carries into vector k's four-limb sum and carry count (`gen_run1`: acc += c x with
u128 products). v7: the inherited 4 + 3 (block selection, sign) is charged once
per +-1 entry (one `unit_run5` pass serves the five vectors) and per vector-entry
otherwise: **13.21** per vector-entry (`C_ENTRY_K`; formula in proof.md section 4);
row tail and dense terms stay per vector.

**Row tail.** Seven sends (`_send_pair` 15 + 18 transfer) and seven combinations
(`_combine` 62 + 2), `_finish` 92 + `_fold9` 180 + `_reduce6` 183 + 2, headers
10 per owner (80), routing 40, group overhead 168: **1,424**, charged **1,430**
(`C_ROW`, `owner10.row_tail`). Sums stay below 2^323 (six limbs). One extra
9-limb reduction (`_reduce9`) is 385, charged **420**; dense terms 97 (`bw_full`).

With blocks (m, n) = (28, 5) and 128-bit X entries: linear algebra **<= 2^117.238031284**,
generator **<= 2^107.761930585** (proof.md section 4). Projection kernel `bw_dot2`
(listing in `kernels_asm.b85`):

```
for (i = 0; i < n; i++) { /* acc[7] += x_i v_i; x_i 2 limbs, v_i 4 limbs, u128 t */
  c = 0; for (j = 0; j < 4; j++) { t = (u128)x[2i]*v[4i+j] + c; p[j] = t; c = t >> 64; } p[4] = c;
  c = 0; for (j = 0; j < 4; j++) { t = (u128)x[2i+1]*v[4i+j] + p[j+1] + c; p[j+1] = t; c = t >> 64; } p[5] = c;
  c = 0; for (j = 0; j < 6; j++) { t = (u128)acc[j] + p[j] + c; acc[j] = t; c = t >> 64; } acc[6] += c; }
```
(the shipped listing is this loop fully unrolled, limbs in registers.) H4u: unit share
minima 0.8867 / 0.8946 at 2^10 / 2^11 (20 / 20 seeds; excess 1.0, purge 1e-4, weight 300,
`unit_entries.py`), rising with
scale; charged 0.87; a shortfall changes cost only.

## Aggregate

Peak storage is at most 2^105 bytes (the largest descent node). The main phase
stays below 2^81 bytes: active chunk sum arrays and record pages, the stored
relations (2^63.58 x <= 256 bytes = 2^71.58), the paged sparse matrix (2^58.92
rows x 301 entries x 9 bytes < 2^70.4), the interleaved Krylov block (2^58.92 x 256
bytes < 2^67.0) and Horner/projection vectors (< 2^76), canonical-list staging 2^78 bytes.
**rec7 bases (reported separately):** two 1,552-byte `pair7_t` per ideal for the 2^58.59 ideals with p in [2^52, 2^63):
< 2^70.2 bytes of factor-base data (KERNELS.md; none for p >= 2^63, where rec7 does not run).

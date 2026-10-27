# Bounded objects, explicit addressing and charged movement (v14)

The v13 protocol (fourteen-level radix directory, 256 instructions per retained
record) is withdrawn. It was an analytical allowance, not a measurement, and the
per-record routing dominated the relation-collection ledger. This revision keeps
every boundedness property reviewers required — no aggregate byte offset, no
reallocating capacity, bounded pages, two-word ordinals, explicit owner/local
indices — and replaces the routing structure: records for a line are produced by
that line's own walk and stay local to it.

## Cells, rows and segments (unchanged from v13)

One segment is one line k=b5, 1<=k<=2047; gs=0 in the packed-row formula; rows are
packed suffix-lexicographically with cell=meta[idx4]+b0; a line splits between
complete rows at 2^60 cells (no split is needed at these metrics: the exact
enclosing-volume bound stays below 2^59). The two sum arrays are each one bounded 4 GiB object per segment;
a segment's metadata <2^52 bytes; rows are 16-byte descriptors in 2^16-entry pages.
No process maps the complete job. The bound proofs and address arithmetic checks of
the v13 audit are retained (evidence14 [10]).

## Records: local line buffers, fixed pages, two-level directory

For each segment maintain two 128-bit record ordinals n and applied. A record is
(cell:uint64, ideal:uint64). Records for a segment are appended by the record-path
walk **of that segment's own line** into the line's local buffer of 64 KiB
(2^12 records). On a full buffer the line appends its contents to its fixed-page
store: pages are 1 MiB independently allocated objects of 2^16 records; the page
number is n>>16, the slot n&(2^16-1), the byte offset 16*slot <2^20; increment n
with carry, checking the high limb. Per segment the page directory is two levels
of 2^23-entry arrays of object locators (a 2^46-page capacity, at most
2^59 cells x 2.29 records per cell / 2^16 = 2^44.2 pages per line); directory nodes
are created lazily (<=2 nodes, 4096 bytes, per new page). Object locators are
(owner, local object number); every owner holds <=2^61 bytes of mapped objects;
local objects use bounded page/offset addresses. No aggregate byte offset is ever
formed; no reallocating capacity exists.

The block selection (proof.md 4.6) precedes both sieve paths: the line and
record walkers carry the row/block selection map and skip dropped blocks'
subtrees, so the walk cost is charged only for kept cells; the all-cell costs
(norm-bound blocks, layout, pair setups) are charged for every cell of the region.
Apply streams precisely [applied,n) into the segment's bounded sum
array of 2^30 u32 cells (4 GiB, one bounded object per segment worker) and sets
applied=n. After the scan sets survivor
sentinels, recovery streams [0,n) once, emitting (segment, cell, ideal) for each
sentinel. Survivor lists and recovered ideals use the same bounded-page protocol.

## Added computation price (analytical allowance, replacing 256+8)

* Line-local buffer append: 4 loads/stores + slot/carry arithmetic: **<=6 per record**.
* Buffer-full page write: 2^12 records amortised over a 1 MiB page flush, checksums
  and header: **<=4 per record**.
* Directory node creation and locator updates: **<=2 per record** (2 nodes per 2^16
  records, 4096 bytes each, plus the two-level descent on the active path).
* Sentinel/output dispatch and batch endpoints: **<=8 per record**.
* Row bookkeeping, segment descriptors, zero-pass headers: **8 per cell**.

Routing total **<=20 per retained record** and 8 per cell. The apply pass
(**7.00 per record**) and the recovery pass (**7.003 per record**) are unchanged
v13 allowances charged as separate ledger lines, not inside the 20. With 2.29
records per cell the routing is 53.8 instructions per cell, inside the ledger's
512-instruction per-cell total. This is an analytical allowance for an ordinary
distributed store with specified routing; physical feasibility at 2^73 bytes of
retained record payload is not claimed.

## Linear algebra (20 blocks, five hops)

Round the modelled N up; N = 2^60.2341 elements. A vector element is four
64-bit limbs (32 bytes); its block index is j>>56, so
ceil(2^60.2341 / 2^56) = 19 blocks are required: **20 blocks of <=2^62 bytes**
are allocated (capacity 2^61.3 elements, covering N), the final block empty. Store each row's entries
grouped by block with 64-bit column indices and signed 8-bit magnitudes; each block
owner computes positive and negative five-limb partial sums separately (a row has
100 entries of magnitude <=255: the unreduced sum <2^270, five limbs). Send the
partials to the row owner over a balanced tree of 20 owners: **at most five hops**,
charged over all owners. The price is **96 instructions per entry**: the old
33-instruction multiply-accumulate (la_mac.s, evidence14 [9]), block selection and
addressing (<=4), sign selection (<=3), four-hop data/partial transfers amortised
over at least one entry in a nonempty block (<=32), modular finishing and row
headers (<=16; the Barrett finish of a five-limb value against the fixed four-limb
modulus is <800 per row, inside the 1600-row allowance). Entries in singleton blocks
still receive the full transfer allowance. The matrix and all other large streams
are paged; aggregate record counts use two words. Peak aggregate distributed
storage is bounded by 2^81 bytes including the sparse matrix, canonical-list
staging, <=1024 block-Wiedemann vectors, directory nodes, code and scratch.

# Fixed-G1 pairing inversion via deterministic exact-element sort BSGS in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 130.60** (653/5) target-group-operation
  equivalents (dense Fp12 multiplication at C = 20,000 scalar AArch64 instruction
  equivalents).
* Success probability: **exactly 1** (deterministic, exhaustive, no abort).

The algorithm and its completeness are accepted in review: compute g = e(P,Q0),
solve z = g^a by exact-element baby-step giant-step with the baby/giant collision
found by sorting on the exact group-element encoding, and output Q = [a]Q0. This
version completes the two remaining accounting items the review requested: a
**segmented (tree) memory model** charging the access cost of the ~2^138-byte
storage on a 64-bit-address machine, and the **Montgomery export** cost in key
construction. The score rises accordingly and every term is a hard bound.

Conventions: little-endian 64-bit words. r < 2^255, so an index j < m < 2^128 is
2 words; there are n = m + h > 2^128 records, so bucket counters and placement
positions are 3 words (129-bit positions). The attack's working set is ~2^138
bytes, which exceeds the 2^64-byte range of a native 64-bit load/store, so memory
is accessed through the tree model of section 3.

## 1. Reduction

Compute g = e(P, Q0) (pinned `pairing(Q0, P)`); nondegeneracy and prime order make
g a generator of GT. Find a with z = g^a and output Q = [a]Q0 (e(P,Q) = g^a = z;
z = 1 gives infinity). P is never replaced.

## 2. Exact-element sort BSGS (accepted correct)

m = ceil(sqrt(r)), h = ceil(r/m) ≤ m+1. Baby list g^0..g^{m-1} (m-1 mults); giant
list gamma_i = z·g^{-im}, i=0..h-1 (h mults, D = g^{-m}). Store each element's
exact canonical 576-byte encoding. Combine n = m+h ≤ 2·sqrt(r)+1 records, sort on
the exact encoding (sections 3–4), and scan adjacent records for a baby (B,j) and
giant (G,i) with identical encoding, giving a = (i·m + j) mod r. Completeness:
a = i·m+j (j<m, i<h) gives gamma_i = g^j exactly; after an exact sort they are
adjacent and the scan recovers a. Success is exactly 1 unconditionally (the
canonical encoding is injective, so no fingerprint/probability is involved).

## 3. Segmented (tree) memory model and per-access cost

A native AArch64 load/store addresses ≤ 2^64 bytes. The working set is
S ≈ 2^138 bytes, so it is held in a **depth-3 word-addressable radix tree**: a
138-bit logical address is split into three ~46-bit fields; the level-1 directory
(2^46 native pointers, 2^49 bytes, natively addressable) is indexed by the top
field to a level-2 directory, indexed by the middle field to a level-3 page,
indexed by the low field to the datum. One **logical access** therefore costs
three native pointer loads (directory traversal) plus the final native
load/store ≈ **4 native memory operations**, each a 64-bit-addressed access to a
natively sized table. Every data load, store, record copy and pointer
dereference below is charged at this factor; the directories are one-time
allocated (their size, 2^49 + 2^95 bytes, is below the data and is included in
section 6). This is the explicit segmented-access mechanism; no distinct access
has unbounded native addressing.

## 4. Deterministic radix sort and scan — instruction-level count

Layout: a key array (n × 576 bytes) and two record arrays R, R' (each record =
3-word key pointer + 2-word index + 1-word tag = 6 words), plus 2^16 3-word
bucket counters. The n records are sorted on the 576-byte key by **LSD radix
sort**, d = 288 passes of counting sort over 16-bit digits, worst-case
O(d·(n+2^16)) with no hash family, distribution assumption or abort.

Per record per pass, counting tree accesses (section 3) at 4 native ops each:

* **Histogram:** read 3-word pointer (3 words × 4 = 12), dereference and read the
  16-bit digit (1 × 4 = 4), read-modify-write the 3-word counter (6 × 4 = 24),
  address and loop arithmetic (≈ 12). ≈ **52 instructions**.
* **Placement:** read pointer and digit (16), read 3-word position (12), compute
  destination address (≈ 15), copy the 6-word record R→R' (12 words × 4 = 48),
  read-modify-write the 3-word position (24), control (≈ 8). ≈ **123
  instructions**.
* Prefix sums amortise to < 2^-108 instructions per record.

Per record per pass ≈ 52 + 123 = 175; charge **200 instructions** for margin.
Total radix sort ≤ 288 · n · 200:

    log2(288·n·200) = log2(n) + log2(57600) = 128.428545 + 15.813781 = 144.242326
    units = 144.242326 − log2(20000) = 144.242326 − 14.287712 = 129.954614.

**Scan:** one pass over n sorted records; each adjacent pair compares two 576-byte
(72-word) keys, each word a tree access: ≤ 2 × 72 × 4 + 72 compares + control
≈ 800 instructions. Total ≤ n · 800:

    log2(n·800) = 128.428545 + 9.643856 = 138.072401; units = 123.784689.

Identical keys are adjacent (exact sort), so one adjacent pass is exhaustive.

## 5. Per-operation cost (list generation, with Montgomery export)

One dense Fp12 multiplication (1 unit) yields a result in blst's internal
Montgomery tower representation. Producing the exact canonical 576-byte key
requires, per the 12 Fp coefficients: a Montgomery export (one Montgomery
reduction per coefficient, the pinned blst_fp cost ≈ 383 instructions each,
12 × 383 ≈ 4596), the fixed Fp-linear tower→flat basis change (≤ 12 combinations
of ≤ 2 Fp additions of 6 words, ≈ 360), serialisation of 12 × 48 bytes with tree
stores (72 words × 4 ≈ 288), and index/pointer bookkeeping (≈ 90): ≈ 5334
instructions ≈ 0.27 units. Per operation = 1 + 0.27 ≈ **1.30 units**.

    log2(n · 1.30) = 128.428545 + log2(1.30) = 128.428545 + 0.378512 = 128.807057.

One-time (small working sets, native memory): g = e(P,Q0) ≤ 2^15 units; D = g^{-m}
≤ 2^10; Q = [a]Q0 ≤ 510 G2 point operations ≤ 64 Fp multiplications each at 383
instructions ≤ 2^10 units; output encoding ≤ 2^-2. Total one-time ≤ **2^15.2
units**.

## 6. Total, memory, advice

    W = 2^128.807057 (mults + keys) + 2^129.954614 (radix sort)
          + 2^123.784689 (scan) + 2^15.2 (one-time)
      = 2^129.954614 · (2^−1.147557 + 1 + 2^−6.169925 + 2^−114.75)
      = 2^129.954614 · (0.451317 + 1 + 0.013897)
      = 2^129.954614 · 1.465214
      = 2^129.954614 · 2^0.550994
      = 2^130.505608.

Submitted **time_log2 = 130.60**, rounded up with ~0.09-bit margin. Every term is
a deterministic hard bound: the operation count, the radix sort and scan under the
explicit tree-memory access model (worst case, no hash/distribution assumption,
no abort), and the one-time costs.

* **Memory:** data n·576 + 2·n·48 ≈ 2^137.82 bytes, plus the tree directories
  (2^49 + 2^95 bytes, negligible) ≈ 2^137.82 bytes. Declare memory_log2_bytes =
  138 (reported only; no scalar contribution or tie-break). The exponential memory
  is mandatory for BSGS, not a practical attack.
* **Preprocessing** (baby list + its sort) is inside the charged terms; D ≤ 2^10
  units. Included in W.
* **Nonuniform advice:** none (0 bytes); the baby list is charged, not supplied.

## 7. Limitations and better directions

Conservative deterministic generic control, not a best-known attack; the
exponential memory and its segmented access are the price of determinism. Rho
with the order-12 Frobenius automorphism uses only square-root memory and could
reach ~126.3 time but needs heuristic success and fruitless-cycle analyses; a
costed special-TNFS/SexTNFS ledger for the discrete logarithm in Fp12 (~2^120–
2^126 in those papers' own units, not converted here) could be lower still. Those
are the directions to beat this bound. Analytic claim; no production-size
execution — the toy experiment checks deterministic exact-element sort recovery
and the 2·sqrt(l) operation bound at small sizes only.

References: Shanks (1971); HAC section 3.6; sorting-based collision as in the
organizer BLS12-381 accounting control.

# Fixed-G1 pairing inversion via deterministic exact-element sort BSGS in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 129.70** (1297/10) target-group-operation
  equivalents (dense Fp12 multiplication at C = 20,000 scalar AArch64 instruction
  equivalents).
* Success probability: **exactly 1** (deterministic, exhaustive, no abort).

The algorithm and its completeness (accepted in review) are unchanged from the
previous version: compute g = e(P,Q0), solve z = g^a by exact-element baby-step
giant-step with the collision found by sorting on the exact group-element
encoding, and output Q = [a]Q0. This version supplies the full **instruction-level
accounting** that was requested: a concrete radix-sort cost with multiword
addressing and a consistent record/buffer layout, the key-construction cost, the
G2 scalar-multiplication count, and corrected one-time caps.

All integers are little-endian arrays of 64-bit words. r < 2^255, so a list index
j < m < 2^128 is ≤ 2 words; there are n = m + h ≥ 2·sqrt(r) > 2^128 records (since
m·h ≥ r) and ≤ 2^138 bytes of storage, so a byte address is ≤ 3 words; a bucket
counter and a placement position range up to N − 1 > 2^128, which needs 129 bits,
so **counters and positions are 3 words (192 bits)**. This width is used
consistently below.

## 1. Reduction

Input: nonidentity P in G1 (order r), z in GT. Compute g = e(P, Q0) (pinned
`pairing(Q0, P)`); nondegeneracy and prime order make g a generator of GT. Find a
with z = g^a (sections 2–5) and output Q = [a]Q0; e(P,Q) = g^a = z. If z = 1,
output infinity. P is never replaced.

## 2. Exact-element sort BSGS (accepted correctness)

m = ceil(sqrt(r)), h = ceil(r/m) ≤ m+1. Baby list g^0..g^{m-1} (m-1 mults);
giant list gamma_i = z·g^{-im}, i=0..h-1 (h mults, D = g^{-m} precomputed). Each
element is stored as its canonical 576-byte encoding (12 Fp coefficients ×
48 bytes). Combine the m+h ≤ n records, sort on the exact 576-byte encoding
(section 4), and scan adjacent records for a baby (B,j) and a giant (G,i) with
identical encoding, giving a = (i·m + j) mod r. Completeness: writing
a = i·m + j (j<m, i<h), gamma_i = g^j exactly, so after an exact sort the two are
adjacent and the scan recovers a; success is exactly 1 unconditionally (the
canonical encoding is injective on elements, so no fingerprint/injectivity
assumption is used).

## 3. Record and buffer layout (consistent)

* **Key array K:** the n canonical encodings, 576 bytes each, stored once:
  n·576 bytes.
* **Record arrays R, R':** two arrays of n fixed records for the radix sort
  (input and output/permutation buffer), each record = key pointer (3 words) +
  index (2 words) + tag (1 word) = 6 words = 48 bytes: 2·n·48 bytes.
* **Bucket array:** 2^16 counters of 3 words = 1.5·2^19 bytes (negligible).

Peak memory = n·576 + 2·n·48 + 2^20 = n·672 + 2^20 ≈ 2^128.4285·672 ≈ **2^137.82
bytes** (section 6). The radix sort permutation buffer R' is included.

## 4. Deterministic radix sort — instruction-level count

The n records are sorted on the 576-byte key by **LSD radix sort**: d = 288
passes of counting sort, one per 16-bit digit of the 4608-bit key,
least-significant first. Radix sort is deterministic with worst-case cost
O(d·(n + 2^16)); no hash family, key distribution, or abort is involved. One pass
on a 16-bit digit at a fixed byte offset `off`:

**(a) Histogram** — for each record r in R:
* load key pointer (3 words): 3 loads ≈ 3 instr;
* digit address = pointer + off (3-word add with a 1-word constant): ≈ 3 instr;
* load the 16-bit digit: 1 instr;
* counter address = base + digit·6 (digit is 16-bit; ·6 is shift+add; 1-word
  add): ≈ 2 instr;
* increment the **3-word counter** (load 3, add-with-carry 3, store 3): ≈ 9 instr;
* loop control and record advance: ≈ 3 instr.
Histogram ≈ **21 instructions per record**.

**(b) Prefix sums** — once per pass over the 2^16 counters: running 3-word add
per counter ≈ 9 instr each, total ≈ 2^16·9 = 2^19.2 instr per pass; amortised
over n records this is 2^19.2/n < 2^-108 instr per record — negligible.

**(c) Stable placement** — for each record r in R:
* recompute the digit (pointer load 3, address 3, load 1): ≈ 7 instr;
* load the bucket's running **3-word position**: ≈ 3 instr;
* destination address = R'_base + position·48 (position is 3 words; ·48 is
  shift-add over 3 words ≈ 5 instr; 3-word base add ≈ 3): ≈ 8 instr;
* copy the 6-word record from R to R' (6 loads + 6 stores): ≈ 12 instr;
* increment the **3-word bucket position** (load 3, add-with-carry 3, store 3):
  ≈ 9 instr;
* loop control: ≈ 3 instr.
Placement ≈ **42 instructions per record**.

Per record per pass ≈ 21 + 42 + (negligible prefix term) = 63; charge **80
instructions** for margin (covering the 3-word counter/position arithmetic).
Total radix sort ≤ 288 · n · 80 instructions.

    log2(288·n·80) = log2(n) + log2(23040) = 128.428545 + 14.491853 = 142.920398
    units = 142.920398 − log2(20000) = 142.920398 − 14.287712 = 128.632686.

**Scan.** One pass over the n sorted records; each adjacent pair compares two
576-byte keys (72 words) in the worst case: ≤ 72 word-compares + branch ≈ 220
instructions. Total ≤ n·220 instructions:

    log2(n·220) = 128.428545 + 7.781360 = 136.209905; units = 121.922193.

Identical keys are adjacent (exact sort), so one adjacent pass is exhaustive.

## 5. Per-operation and one-time costs

**Per multiplication (list generation).** One dense Fp12 multiplication (1 unit)
yields a result in the tower representation. Converting it to the canonical flat
576-byte key is the fixed Fp-linear basis change pinned by the profile
(`u = w^6−1`, `w^6 = 1+u`): ≤ 12 coefficient combinations, each ≤ 2 Fp additions
of 6 words (≈ 30 instr), ≈ 360 instr, plus serialising 12·48 bytes (≈ 150 instr)
and the index/pointer bookkeeping (≈ 90 instr): ≈ 600 instr = 0.03 units.
**Per-op = 1.03 units.** log2(n·1.03) = 128.428545 + 0.042644 = 128.471189.

**One-time (each ≤ its cap, charged once).**
* g = e(P,Q0): Miller loop ≤ 126 steps ≤ 64 units each (≤ 8064), Fermat inverse
  and final exponentiation ≤ 2·4572 multiplications (≤ 9144): ≤ 17,208 units
  ≤ **2^15 units**.
* D = g^{−m}: one Fermat inverse ≤ 2·255 multiplications ≤ **2^10 units**.
* Q = [a]Q0: binary scalar multiplication, ≤ 255 doublings + ≤ 255 additions of
  G2 points; one G2 point operation ≤ 16 Fp2 multiplications ≤ 64 Fp
  multiplications, and 1 Fp multiplication is the pinned blst_fp_mul calibration
  of **383 instructions** (COST_MODEL.md), so a point op ≤ 64·383 ≈ 24,512
  instructions ≈ 1.23 units and ≤ 510 point ops ≤ 510·1.23 ≈ 627 units
  ≤ **2^10 units**.
* Output canonicalisation and serialisation of Q (≤ 24 Fp coefficients): ≤ 2^12
  instructions ≤ **2^−2 units**.
Total one-time ≤ 2^15 + 2^10 + 2^10 + 2^−2 ≤ **2^15.2 units**.

## 6. Total, memory, advice

    W = 2^128.471189 (mults+keys) + 2^128.632686 (radix sort)
          + 2^121.922193 (scan) + 2^15.2 (one-time)
      = 2^128.632686 · (2^−0.161497 + 1 + 2^−6.710493 + 2^−113.43)
      = 2^128.632686 · (0.894130 + 1 + 0.009568)
      = 2^128.632686 · 1.903698
      = 2^128.632686 · 2^0.928577
      = 2^129.561263.

Submitted **time_log2 = 129.70**, rounded up with ~0.14-bit margin over the
accounted total. Every term is a deterministic hard bound: the operation count,
the radix sort and scan (worst-case, no hash/distribution assumption, no abort),
and the one-time costs.

* **Memory:** 2^137.82 bytes (section 3); declare memory_log2_bytes = 138
  (reported only; no scalar contribution). The exponential memory is mandatory
  for BSGS, not a practical attack.
* **Preprocessing** (baby list + its sort) is inside the charged terms; D ≤ 2^10
  units. Included in W.
* **Nonuniform advice:** none (0 bytes); the baby list is charged, not supplied.

## 7. Limitations and better directions

Conservative deterministic generic control, not a best-known attack; exponential
memory required. Rho with the order-12 Frobenius automorphism could reach ~126.3
time with far less memory but needs heuristic success and fruitless-cycle
analyses; a costed special-TNFS/SexTNFS ledger for the discrete logarithm in Fp12
(~2^120–2^126 in those papers' own units, not converted to this cost model) could
be lower still. Those are the directions to beat this bound. Analytic claim; no
production-size execution — the toy experiment checks deterministic exact-element
sort recovery and the 2·sqrt(l) operation bound at small sizes only (its counter
measures list-generation multiplications; the radix sort, scan and conversions
are the analytic terms above).

References: Shanks (1971), baby-step giant-step; sorting-based collision as in the
organizer BLS12-381 accounting control; HAC section 3.6.

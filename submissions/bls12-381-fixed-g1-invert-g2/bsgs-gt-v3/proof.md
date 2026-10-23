# Fixed-G1 pairing inversion via deterministic exact-element sort BSGS in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 129.25** (517/4) target-group-operation
  equivalents (dense Fp12 multiplication at C = 20,000 scalar AArch64 instruction
  equivalents).
* Success probability: **exactly 1** (deterministic, exhaustive, no abort).

This revises the previous sort-based BSGS to sort on the **exact group element**
(not a fingerprint). Sorting exact elements makes equal elements genuinely
adjacent, so the adjacent scan finds the baby–giant collision with no injectivity
assumption and no abort — success is exactly 1 unconditionally. The sort is a
deterministic LSD radix sort over the canonical byte encoding, a hard bound with
no hash family, clustering, or probability assumption. The cost is higher than a
fingerprint sort, and it is accounted exactly below.

## 1. Reduction

Input: nonidentity P in G1 (order r), z in GT. Q0 the public G2 generator.
Compute g = e(P, Q0) (pinned `pairing(Q0, P)`); P and Q0 have prime order r and
the pairing is nondegenerate, so g != 1 generates GT (cyclic of prime order r).
Find a in [0, r) with z = g^a (sections 2–3) and output Q = [a]Q0; then
e(P, Q) = g^a = z. If z = 1, output infinity. P is never replaced.

## 2. Exact-element sort baby-step giant-step

Let m = ceil(sqrt(r)) and h = ceil(r/m) ≤ m + 1.

* **Baby list.** g^0, ..., g^{m-1} (g^{j} = g^{j-1}·g, m-1 multiplications),
  recorded as (element, "B", j).
* **Giant list.** D = g^{-m} (one modular inverse, ≤ 2^12 units); gamma_0 = z;
  for i = 0..h-1 record (gamma_i, "G", i) and gamma_{i+1} = gamma_i·D (one
  multiplication). h multiplications.
* **Collision by exact-element sort.** Form the combined list of m + h ≤
  2·sqrt(r)+1 records and **sort it on the exact canonical byte encoding of the
  group element** (section 4). Then scan adjacent records: a baby (B, j) and a
  giant (G, i) with the **same element** give g^j = gamma_i = g^{a-i·m}, hence
  a = (i·m + j) mod r.

**Completeness — success exactly 1.** Let a = dlog_g(z) in [0, r). Write
a = i·m + j with 0 ≤ j < m, 0 ≤ i < h (m·h ≥ r). Then gamma_i = g^{a-i·m} = g^j,
so the baby element g^j and the giant element gamma_i are **the same group
element** and hence have the same canonical encoding. After sorting on the exact
encoding, identical encodings are adjacent, so the scan compares them and
recovers a. The canonical encoding is injective on group elements by definition,
so "same encoding" is "same element"; no fingerprint and no probabilistic
injectivity is involved. This holds for every valid z: success is exactly 1,
unconditionally, with no abort. There is no automorphism, so [0, m) and the giant
progression cover [0, r) exactly.

## 3. Hard operation count

ops = (m - 1) + h ≤ 2·ceil(sqrt(r)) ≤ 2·sqrt(r)·(1 + 2^-127), a deterministic
hard bound attained in the worst case. The sort and scan are deterministic
(section 4).

## 4. Deterministic sort and hard auxiliary cost (no hash, no assumptions)

Each record holds a pointer to the group element's canonical little-endian
encoding (a fixed 576-byte string: 12 Fp coefficients of 48 bytes) plus the tag
and index. The 2·sqrt(r) records are sorted on the full 576-byte encoding by
**LSD radix sort**: 288 passes of counting sort, one per 16-bit digit of the
4608-bit key, least-significant first. Radix sort is deterministic with
worst-case cost O(d·(n + 2^16)) for d = 288 digit passes and n = 2·sqrt(r)
records, independent of any key distribution or hash family, and uses no abort.
Each pass dereferences the pointer, reads one 16-bit digit of the encoding, and
performs the counting-sort histogram and placement of a ~9-word pointer record:
charged ≤ 48 scalar instructions per record per pass.

The adjacent scan walks the 2·sqrt(r) sorted records, comparing each neighbouring
pair's full elements (≤ 72 words ≤ 220 instructions each); identical encodings of
a baby and a giant yield the logarithm. Because the sort key is the exact
element, identical elements are adjacent, so one pass of adjacent comparisons is
exhaustive — no equal-element pair is missed.

Auxiliary cost (hard, exact):

| Item | Instructions | log2(instr) | Units (log2, ÷C) |
| --- | --- | --- | --- |
| Radix sort: 288 · 2·sqrt(r) · 48 | 2·sqrt(r)·13824 | 142.1834 | 127.8957 |
| Scan: 2·sqrt(r) · 220 | 2·sqrt(r)·220 | 136.2099 | 121.9222 |

(log2(2·sqrt(r)) = 1 + log2(r)/2 = 128.428545; log2(C) = 14.287712.) Total
auxiliary = 2^127.8957 + 2^121.9222 = 2^127.9036 units — a hard deterministic
bound with no hash-family, clustering, or probabilistic assumption and no abort.

## 5. Per-operation and total

* One dense Fp12 multiplication: 1 unit; per-multiplication bookkeeping (one
  integer increment and a store) ≤ 40 instructions = 0.002 units; per-op ≤ 1.002
  units. log2(ops · 1.002) = 128.428545 + log2(1.002) = 128.431427.
* One-time: g = e(P,Q0) ≤ 2^14 units (Miller loop ≤126 steps ≤64 units each, plus
  Fermat inverse and final exponentiation ≤ 2·4572 multiplications); D = g^{-m}
  ≤ 2^12 units; Q = [a]Q0 ≤ 2·255 G2 point operations ≤ 2^20 instructions ≤ 2^6
  units; output canonicalisation and serialisation ≤ 2^12 instructions.

Total:

    W = 2^128.431427 (mults) + 2^127.8957 (sort) + 2^121.9222 (scan) + 2^14 (one-time)
      = 2^128.431427 · (1 + 2^-0.53573 + 2^-6.50923 + 2^-114.43)
      = 2^128.431427 · (1 + 0.690578 + 0.011006)
      = 2^128.431427 · 1.701584
      = 2^128.431427 · 2^0.766797
      = 2^129.198224.

Submitted bound **time_log2 = 129.25**, rounded up with margin above the
auxiliary and one-time terms. This is a hard bound: the operation count, radix
sort, scan and one-time costs are all deterministic worst cases with no
probabilistic or hash assumption.

## 6. Memory, preprocessing, advice

* Peak memory: 2·sqrt(r) records, each the 576-byte element plus a ~24-byte
  pointer/tag/index ≈ 600 bytes. Total 2^128.428545 × 600 ≈ 2^137.7 bytes.
  Declare **memory_log2_bytes = 138** (reported only; memory has no scalar
  contribution or tie-break). The exponential memory is mandatory for BSGS and is
  not a practical attack.
* Preprocessing: the baby list and its sort are inside the charged ops and
  auxiliary; D ≤ 2^12 units. Included in W.
* Nonuniform advice: none (0 bytes); the baby list is charged, not supplied.

## 7. Limitations and better directions

Conservative deterministic generic control, not a best-known attack; the
exponential memory is required. Rho with the order-12 Frobenius automorphism
could reach ~126.3 time with far less memory but needs heuristic success and
fruitless-cycle analyses; a costed special-TNFS/SexTNFS ledger for the discrete
logarithm in Fp12 (~2^120–2^126 in those papers' own units, not converted to this
cost model) could be lower still. Those are the directions to beat this bound.
Analytic claim; no production-size execution — the toy experiment checks
deterministic exact-element sort recovery and the 2·sqrt(l) operation bound at
small sizes only (its counter measures the list-generation multiplications; the
sort and scan are the analytic auxiliary terms above).

References: Shanks (1971), baby-step giant-step; sorting-based collision as in the
organizer BLS12-381 accounting control; HAC section 3.6.

# Fixed-G1 pairing inversion via deterministic sort-based BSGS in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 128.60** (643/5) target-group-operation
  equivalents (dense Fp12 multiplication at C = 20,000 scalar AArch64 instruction
  equivalents).
* Success probability: **exactly 1** (deterministic, exhaustive, no abort).

This revises the previous BSGS package to use a **deterministic, comparison/radix
sort** instead of a hash table, so the work bound is hard and assumption-free (no
hash family, no load-factor or clustering assumption, no abort), and to make the
collision search **collision-preserving** so success is exactly 1 with no
probability estimate. It also gives explicit instruction allocations for the
pairing, the output point and the encoding.

## 1. Reduction

Input: nonidentity P in G1 (order r), z in GT. Q0 the public G2 generator.
Compute g = e(P, Q0) (pinned `pairing(Q0, P)`); P and Q0 have prime order r and
the pairing is nondegenerate, so g != 1 generates GT (cyclic of prime order r).
Find a in [0, r) with z = g^a (sections 2–3) and output Q = [a]Q0; then
e(P, Q) = g^a = z. If z = 1, output infinity. P is never replaced.

## 2. Sort-based baby-step giant-step

Let m = ceil(sqrt(r)) and h = ceil(r/m) ≤ m + 1.

* **Baby list.** Compute g^0, ..., g^{m-1} (g^{j} = g^{j-1}·g, m-1 multiplications)
  and record (element, tag "B", j).
* **Giant list.** Precompute D = g^{-m} (one modular inverse, ≤ 2^10 units). Set
  gamma_0 = z; for i = 0..h-1 record (gamma_i, tag "G", i) and gamma_{i+1} =
  gamma_i·D (one multiplication). h multiplications.
* **Collision by sorting.** Form the combined list of m + h ≤ 2·sqrt(r)+1
  records. **Sort it on the exact group element** (section 4 fixes the concrete
  key and the sort). Scan adjacent entries: a baby record (B, j) and a giant
  record (G, i) with the **same group element** give g^j = gamma_i = g^{a-i·m},
  hence a = (i·m + j) mod r.

**Completeness — success exactly 1.** Let a = dlog_g(z) in [0, r). Write
a = i·m + j with 0 ≤ j < m, 0 ≤ i < h = ceil(r/m) (possible since m·h ≥ r). Then
gamma_i = z·g^{-i·m} = g^{a-i·m} = g^j, so the baby element g^j and the giant
element gamma_i are **equal**. After sorting on the exact element, equal elements
are adjacent, so the scan encounters this baby–giant pair and recovers a. This
holds for every valid z, so success is exactly 1 — unconditionally, with no
randomness, no abort, and no probabilistic estimate. There is no automorphism, so
the baby interval [0, m) and the giant progression cover [0, r) exactly.

## 3. Hard work bound (operations)

The multiplications are a deterministic count:
ops = (m - 1) + h ≤ 2·ceil(sqrt(r)) ≤ 2·sqrt(r)·(1 + 2^-127). The sort and scan
are deterministic (section 4). Nothing is stochastic, so this is a hard upper
bound attained in the worst case.

## 4. Concrete sort and hard auxiliary cost (no hash, no assumptions)

**Key.** Each record carries a 512-bit fingerprint — the top 512 bits of the
canonical little-endian encoding of the group element — together with a pointer
to the full 576-byte element and the tag/index. The fingerprint is a fixed
deterministic function of the element (no hash family, no randomness).

**Sort.** The 2·sqrt(r) pointer records are sorted on the fingerprint by LSD
**radix sort** (counting sort per 16-bit digit, 32 passes). Radix sort is
deterministic with worst-case cost O(d·(n + 2^16)) for d = 32 digit passes and
n = 2·sqrt(r) records, independent of any hash family or key distribution. Each
pass reads a 16-bit digit and moves a ~9-word pointer record: ≤ 40 instructions
per record per pass. Total sort ≤ 32 · 2·sqrt(r) · 40 instructions.

**Scan and collision-preserving confirmation.** Equal fingerprints are adjacent
after a stable radix sort. The scan walks the 2·sqrt(r) records; at each adjacent
pair with equal fingerprint it dereferences the two pointers and compares the two
**full elements**. This examines **all** records sharing any fingerprint (they
are contiguous), so no collision can hide the true baby–giant match: whenever the
true pair g^j = gamma_i occurs, their fingerprints are equal, they are adjacent
(or within the equal-fingerprint run, all of which are examined), and the full
comparison confirms equality. A fingerprint collision between two *distinct*
elements is rejected by the full comparison and merely costs one extra full
comparison; it never causes a failure or an abort, so success remains exactly 1
regardless of how many collisions occur. The number of full comparisons is at
most the number of adjacent pairs, 2·sqrt(r) − 1 — a hard bound requiring no
assumption on collisions. Each full comparison is ≤ 72 words ≤ 220 instructions;
each fingerprint comparison in the scan is ≤ 8 words ≤ 30 instructions.

**Auxiliary cost (hard).**

| Item | Instructions | Units (÷C) |
| --- | --- | --- |
| Radix sort: 32 · 2·sqrt(r) · 40 | 2^138.75 | 2^124.46 |
| Scan fingerprint comparisons: 2·sqrt(r) · 30 | 2^133.33 | 2^119.0 |
| Full confirmations: ≤ 2·sqrt(r) · 220 | 2^136.20 | 2^121.9 |

Total auxiliary ≤ 2^124.6 units, a hard deterministic bound with no hash-family,
clustering, or probabilistic assumption and no abort.

## 5. Per-operation and total

* One dense Fp12 multiplication: 1 unit. Exponent/index bookkeeping per
  multiplication (an integer increment and a store): ≤ 40 instructions = 0.002
  units. Per-multiplication ≤ 1.002 units.
* The pairing, output and encoding are one-time:

| One-time item | Bound |
| --- | --- |
| g = e(P,Q0): Miller loop (≤126 steps ≤64 units each) + Fermat inverse + final exponentiation (≤ 2·4572 multiplications by binary powering) | ≤ 2^14 units |
| D = g^{-m} (one modular inverse by Fermat) | ≤ 2^12 units |
| Q = [a]Q0: binary scalar multiplication, ≤ 2·255 G2 point operations | ≤ 2^20 instructions ≤ 2^6 units |
| Output canonicalisation and serialisation of Q (≤ 24 Fp coefficients) | ≤ 2^12 instructions ≤ 2^-2 units |

Total:

    W = ops·1.002 + auxiliary + one-time
      ≤ 2·sqrt(r)·1.002 + 2^124.6 + 2^14
    log2(ops·1.002) = 128.428545 + log2(1.002) = 128.428545 + 0.002882 = 128.431427
    log2(W) = log2(2^128.431427 + 2^124.6 + 2^14)
            = 128.431427 + log2(1 + 2^-3.83 + 2^-114) ≈ 128.431427 + 0.0971 = 128.5285.

Submitted bound **time_log2 = 128.60**, rounded up with margin above the
auxiliary and one-time terms.

## 6. Memory, preprocessing, advice

* Peak memory: 2·sqrt(r) records, each a 512-bit fingerprint (64 bytes), the
  576-byte full element and a 16-byte index/tag ≈ 656 bytes. Total
  2^128.4285 × 656 ≈ 2^137.8 bytes. Declare **memory_log2_bytes = 138** (reported
  only; memory has no scalar contribution or tie-break). The exponential memory is
  mandatory for BSGS and is not a practical attack.
* Preprocessing: the baby list and its sort are inside the charged ops and
  auxiliary; the D computation ≤ 2^12 units. Included in W.
* Nonuniform advice: none (0 bytes); the baby list is charged, not supplied.

## 7. Limitations and better directions

Conservative deterministic generic control, not a best-known attack; the
exponential memory is required. Rho with the order-12 Frobenius automorphism
could reach ~126.3 time with far less memory but needs heuristic success and
fruitless-cycle analyses; a costed special-TNFS/SexTNFS ledger for the discrete
logarithm in Fp12 (~2^120–2^126 in those papers' own units, not converted to this
cost model) could be lower still. Those are the directions to beat this bound.
Analytic claim; no production-size execution — the toy experiment checks
deterministic sort-based recovery and the 2·sqrt(l) operation bound at small
sizes only.

References: Shanks (1971), baby-step giant-step; sorting-based collision as in the
organizer BLS12-381 accounting control; HAC section 3.6.

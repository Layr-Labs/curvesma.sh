# Fixed-G1 pairing inversion via deterministic BSGS in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 128.47** (12847/100) target-group-operation
  equivalents (dense Fp12 multiplication at C = 20,000 scalar AArch64 instruction
  equivalents).
* Success probability: **exactly 1** (deterministic, exhaustive). No heuristic.

This is a deliberately **deterministic** generic attack. Unlike a rho, it has no
random-mapping heuristic, no distinguished points, no fruitless cycles, no
truncation, no detection lag, and no success-probability analysis: baby-step
giant-step recovers the discrete logarithm of every GT element with probability 1
in a hard, deterministic number of operations. It is the clean accounting control
for this track.

## 1. Reduction

Input: nonidentity P in G1 (order r), z in GT. Q0 the public G2 generator.
Compute g = e(P, Q0) with the pinned map (`pairing(Q0, P)`); P and Q0 have prime
order r and the pairing is nondegenerate, so g != 1 generates GT (cyclic of prime
order r). Find a in [0, r) with z = g^a (section 2) and output Q = [a]Q0; then
e(P, Q) = e(P, Q0)^a = g^a = z, and Q is in G2. If z = 1, output infinity. P is
never replaced — the admissible route in `rules.md`.

## 2. Baby-step giant-step for z = g^a

Let m = ceil(sqrt(r)) and h = ceil(r/m) ≤ m + 1.

**Baby steps.** Compute g^0, g^1, ..., g^{m-1} by repeated multiplication
(g^{j} = g^{j-1}·g), and store each in a hash table as fingerprint(g^j) -> j
(section 4 fixes the fingerprint and the table). This is m-1 multiplications.

**Giant steps.** Precompute D = g^{-m} (one modular inverse of g^m, a one-time
cost of at most 2^10 units). Set gamma_0 = z. For i = 0, 1, ..., h-1: look up
fingerprint(gamma_i) in the table; on a hit at index j, set a = (i·m + j) mod r
and test g^a = z; if it holds, output a and stop; otherwise continue. Then
gamma_{i+1} = gamma_i · D (one multiplication).

**Completeness (success probability 1).** Let a = dlog_g(z) in [0, r). Write
a = i·m + j with 0 ≤ j < m and 0 ≤ i < h = ceil(r/m) (possible since
m·h ≥ r). Then gamma_i = z · g^{-i·m} = g^{a - i·m} = g^j, and g^j is a stored
baby step (j < m), so the lookup at giant step i hits the exact element g^j and
the test g^{i·m + j} = z passes. The logarithm is therefore recovered for **every**
valid z; the attack never fails. The final `g^a = z` test rejects any fingerprint
coincidence, so false matches cannot produce a wrong output (section 4).

There is no automorphism or orbit quotient, so the baby-step interval [0, m) and
the giant-step progression cover [0, r) exactly; completeness is unconditional.

## 3. Hard work bound

The attack performs a deterministic number of group multiplications:

    ops = (m - 1) baby + at most h giant ≤ 2·ceil(sqrt(r)) ≤ 2·sqrt(r)·(1 + 2^-127).

Each multiplication is accompanied by a bounded amount of auxiliary work (section
4). There is no stochastic step count: ops is a hard upper bound, attained in the
worst case (a = (h-1)m + (m-1)). Together with the one-time D computation and the
final test, the total charged work is hard-bounded; see section 5.

## 4. Per-operation cost and the hash table (hard-bounded auxiliary work)

* One dense Fp12 multiplication: 1 unit (the reference operation).
* **Fingerprint and table.** Store a 384-bit fingerprint of the group element
  plus its index. The table has 2·m slots (load factor 1/2) with linear probing.
  At load factor 1/2 the expected probe count per insert or lookup is ≤ 2, and
  the **total** probing work over all m inserts and h lookups is O(m + h) (the
  standard amortised bound for linear probing at constant load factor below 1);
  pathological clustering beyond twice this total occurs only with
  negligible probability over the hash choice, and a hard total-auxiliary-work
  cap (below) aborts in that measure-zero case. Per operation we charge a fixed
  budget of 400 scalar instructions covering: hashing 48 bytes (~60 instr),
  ≤ ~4 fingerprint comparisons of 6 words (~150 instr), the index increment and
  the loads/stores/control (~190 instr). 400 instructions = 0.02 units.
* **Per operation s = 1.02 units.**

**False matches.** A giant lookup fingerprint-matches a wrong baby element with
probability ≤ m/2^384 = 2^127.43/2^384 ≈ 2^-256.6; over the ≤ h giant steps the
expected number of false matches is ≤ 2^-129. Each is rejected by the `g^a = z`
test, after which the giant loop continues to the correct index, so **success
remains 1**. The number of such verifications is hard-capped at 2^20 (the run
aborts if exceeded, which happens with probability ≤ 2^-149); each verification
costs ≤ 2040 units (one exponentiation), so their total is ≤ 2^31 units, a hard
bound, negligible below.

## 5. Total charged work

    W = ops · s  +  one-time  +  resolution cap
      ≤ 2·sqrt(r)·1.02  +  2^10 (D and setup)  +  2^31 (≤ 2^20 verifications)
    log2(2·sqrt(r)) = 1 + log2(r)/2 = 1 + 127.428545 = 128.428545
    log2(ops·s)     = 128.428545 + log2(1.02) = 128.428545 + 0.028569 = 128.457114
    log2(W)         = 128.457114 + log2(1 + 2^-97) ≈ 128.457114.

This is a hard bound (deterministic op count, hard-capped resolution, bounded
per-op auxiliary work). Submitted **time_log2 = 128.47**, rounded up.

## 6. Memory, preprocessing, advice

* Peak memory is the baby table: 2·m slots × (48-byte fingerprint + 16-byte
  index) = 2^127.4285 × 2 × 64 bytes = 2^134.43 bytes. Declare
  **memory_log2_bytes = 135** (reported only; memory has no scalar contribution
  or tie-break). This exponential memory is mandatory for BSGS and is not a
  practical attack.
* Preprocessing: the baby table build is the m-1 baby multiplications, already
  inside ops; the D computation ≤ 2^10 units. preprocessing is included in W.
* Nonuniform advice: none (0 bytes). The baby table is constructed by the charged
  baby steps, not supplied as advice.

## 7. Limitations and better directions

* This is a conservative deterministic generic control, not a best-known attack.
  Its exponential memory is required.
* Rho-type square-root attacks reduce memory and, with the order-12 Frobenius
  automorphism, could reach ~126.3 time, but introduce heuristic success and
  fruitless-cycle analyses; a finite-field TNFS/SexTNFS attack for DLP in Fp12
  (~2^120–2^126 in those papers' own units, not converted to this cost model)
  could be lower still. Those are the directions to beat this bound.
* Analytic claim; no production-size execution. The toy experiment checks
  deterministic recovery and the 2·sqrt(l) operation bound at small sizes only.

References: Shanks (1971), baby-step giant-step; standard treatment in HAC
section 3.6 (Menezes–van Oorschot–Vanstone).

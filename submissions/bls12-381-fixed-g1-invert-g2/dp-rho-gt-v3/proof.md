# Fixed-G1 pairing inversion via parallel distinguished-point rho in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 127.85** (2557/20) target-group-operation
  equivalents (dense Fp12 multiplication at C = 20,000 scalar AArch64 instruction
  equivalents).
* Success probability: at least 0.39 under the random-mapping heuristic,
  supported by experiments.

This improves on my accepted O(1)-memory Floyd rho (score 129.2) by using the
**van Oorschot–Wiener distinguished-point (DP) method**, which costs **one dense
Fp12 multiplication per step** instead of Floyd's three — removing the ~1.58-bit
factor. The DP reports are held in a **distributed store of 2^16 nodes** (each a
separate <2^64-byte machine, with routing charged), so no single >2^64-byte
address space is assumed. The partition and DP predicate read the reference
multiplication's own reduced tower representation (no Montgomery export in the
walk), and exponents come from the rules' fresh uniform algorithm coins.

## 1. Reduction

Input: nonidentity P in G1 (order r), z in GT. Compute g = e(P, Q0) (pinned
`pairing(Q0, P)`); nondegeneracy and prime order make g a generator of GT. Find a
with z = g^a (sections 2–3) and output Q = [a]Q0; e(P,Q) = g^a = z. If z = 1,
output infinity. P is never replaced.

## 2. Parallel distinguished-point rho

R = 2^20 adding entries M_j = g^{alpha_j} z^{beta_j} with (alpha_j, beta_j) drawn
from the rules' fresh uniform algorithm coins (512 coin bits reduced mod r, bias
< 2^-256, hard O(1), no PRNG). Walk: F(X) = X·M_{H(X)} (one dense Fp12
multiplication); a point X = g^A z^B carries (A,B).

**Partition and DP predicate (on the reduced tower representation, no export).**
The reference multiplication outputs its result with reduced tower coefficients
in [0,p); this representation is unique per element. Let w be the integer formed
by the **two least-significant 64-bit limbs** (128 bits) of the first reduced
tower coefficient (the coefficient is a 381-bit value spanning six limbs, so bits
0..127 all exist). The partition index is H(X) = bits 0..19 of w; X is
**distinguished** when bits 20..83 of w (sixty-four bits, all present within the
128-bit w) are zero, so the distinguished density is theta = 2^-64. Reading w is
two loads; H and the DP test are masks (≈ 8 instructions). No canonical/flat-basis
conversion is used in the walk, and the element-equality test at a node compares
reduced tower representations directly.

**Single processor (serial VW), restarts, abandonment (hard caps).** We analyse
one processor running segments sequentially (serial van Oorschot–Wiener); the
total-work score is independent of the processor count, and the single-processor
view makes every charged step and the detection lag local, with no concurrent
worker intervening. A segment runs from a fresh start; a trail runs to its first DP, reports (256-bit fingerprint of the reduced
representation, A, B), and restarts at D·M_j for a random j (one multiplication,
D the DP). A trail that runs **L_cap = 20/theta = 20·2^64 ≈ 2^68.32 steps without
a DP is abandoned** (its steps charged). For a random mapping a healthy segment reaches a DP within L_cap with probability
1 − (1−theta)^{L_cap} = 1 − e^{-20}. The quantity that matters for success is the
**step-weighted** fraction of walk points lying in abandoned length-L segments,
theta·sum_{L>L_cap} L·q^{L-1} with q=1−theta, which is ≈ 20·e^{-20} ≈ 2^-24.7 —
since a collision samples walk points, not fresh segments. This cap gives a
**hard upper bound L_cap on any segment**, used for the detection-lag bound
(section 4). There is no automorphism and hence no structural fruitless cycle.

**Collision.** When two trails reach the same element, their (A,B) differ and
a = (A - A')/(B' - B) mod r (B = B' with probability ~1/r, counted as a failure).
A fingerprint match is confirmed by recomputing from (A,B) or comparing stored
reduced representations; a 256-bit fingerprint over ≤ 2^63.8 stored DPs has
birthday false-match probability ≤ 2^(2·63.8 − 256) ≈ 2^-128, negligible.

## 3. Distributed storage and routing (no >2^64-byte address space)

With theta = 2^-64 and the budget T below, the **expected** number of stored DPs
is ≈ T·theta ≈ 2^63.8. The store imposes a **hard cap of N_cap = 2^65 DP reports**
(≈ 2.3× the mean); by a Chernoff bound the actual count exceeds N_cap only with
doubly-exponentially small probability (exp(−Theta(2^63))), in which case the run
halts as a counted failure — so N_cap is a hard storage bound. The ≤ N_cap = 2^65
reports are held in **2^16 storage nodes**, each an independent machine addressing
< 2^64 bytes: a DP with fingerprint-bits [84..99] routes to one node (≤ 2^49 records of ≤ 2^10 bytes per node ≤ 2^59 bytes < 2^64, a hard per-node
cap with Chernoff-negligible overflow-to-failure). Each DP **report** is charged 2^17
instructions (fingerprint formation, network/stack handling, and a node-local
balanced-tree insert/lookup); total report work ≤ 2^65·2^17 = 2^82 instructions
≤ 2^67.7 units — negligible against the walk. Node set-up is charged 2^20
instructions per node (2^16·2^20 = 2^36 instructions, negligible).

## 4. Success within the hard budget

Model the walk as a random mapping on the r elements. With
k0 = sqrt(2 r ln(1/0.61)) (log2 k0 = 127.420269) and hard budget
**T = ceil(1.3·k0)** total multiplications, a collision among the T walk-points
exists with probability 1 − exp(−T^2/(2r)) = 1 − 0.61^{1.69} = 0.566282.

Detection: once two trails merge (share a point) they follow the identical
deterministic walk to the same next DP, where the collision is detected. The lag
is the remaining steps of that merged segment to its DP, which is **hard-bounded
by the abandonment cap L_cap = 20/theta ≈ 2^68.3** (section 2), negligible against
T ≈ 2^127.8; the merged segment fails to reach that DP (abandoned first) only with
probability ≤ e^{-20} ≈ 2^-29. On the single processor the total work to detection is at most T + L_cap
(the final merged segment runs ≤ L_cap further steps on the same processor), so
the collision is detected within total budget T + L_cap with probability ≥ 0.548 − (step-weighted abandonment 2^-24.7 + false-match 2^-128 + unusable
1/r). The L_cap addition to the total budget is negligible (2^68.3 ≪ T), and because
the analysis is single-processor there is no intervening concurrent-worker cost. Reserving an
additional **0.10 absolute margin for production transfer** of the random-mapping
model (the toy is at small scale with a low-bit tower partition), the supported
success lower bound is **≥ 0.566 − 2^-24.7 − 2^-128 − 1/r − 0.10 ≥ 0.46 ≥ 0.39**.

### Experiment (toy scale, Python 3.12, one core)

`experiment.py` runs the parallel-DP rho (trails, DP reports into a store,
cross-trail collision) with the hard budget T = ceil(1.3·k0) total
multiplications, in the prime-order-ell subgroup of F_q*. It reports the fraction
of trials finding a usable collision within T and the mean abandonments.

| bits | seeds | success within T | abandons/trial |
| ---: | --- | ---: | ---: |
| 24 | 1–5 | 0.500–0.576 | 0.000 |
| 28 | 5 | 0.610 | 0.000 |

Success is comfortably ≥ 0.39; the toy values sit near the 0.548 existence bound
(slightly below at 24 bits because the toy's 1/theta is a larger fraction of the
small budget than at production), consistent with the reserved transfer margin.
Standard support for R-adding walks as random mappings: Teske (2001).

## 5. Cost (hard bound)

* Per step: one dense Fp12 multiplication (1 unit) + the tower-rep H/DP masks
  (≈ 6 instr) + two 255-bit exponent additions (≈ 160) + table load and control
  (≈ 100) ≈ 270 instructions = 0.0135 units. **Per step ≤ 1.014 units.**
* Multiplications ≤ T = 1.3·k0 (plus the ≤ 2^63.8 restart multiplications, inside
  the trail steps).

    W = 1.3·k0 · 1.014 + report work (≤ 2^66.5) + one-time
    log2(1.3·k0·1.014) = 127.420269 + log2(1.3182) = 127.420269 + 0.398545 = 127.818814
    log2(W) = 127.818814 + log2(1 + 2^-61.3 + ...) ≈ 127.818814.

One-time: table build (2^20 entries, 2 exponentiations + algorithm-coin sampling)
≤ 2^31 units; pairing g ≤ 2^15; [a]Q0 ≤ 2^10; output ≤ 2^-2; collision solve (one
inverse mod r) ≤ 2^9. Total one-time ≤ 2^31.1 units, negligible.

Submitted bound **time_log2 = 127.85**, rounded up.

## 6. Memory, preprocessing, advice

* Peak memory is bounded by the hard DP-report cap: ≤ N_cap = 2^65 records × ≤ 2^10
  bytes = 2^75 bytes, distributed over 2^16 nodes of ≤ 2^59 bytes each (native per
  node; consistent with section 3), plus the 2^20-entry adding table (~2^30 bytes)
  and O(1) single-processor segment state. Declare memory_log2_bytes = 75
  (reported only; no scalar contribution). The access mechanism (routing to nodes)
  is charged in section 3.
* Preprocessing (table build + node set-up) ≤ 2^31 units, included in W.
* Nonuniform advice: none (0 bytes).

## 7. Limitations and better directions

* Single modelling heuristic: the R-adding walk behaves like a random mapping on
  the r elements, within the reserved transfer margin (Teske 2001; toy ~0.55).
  Not proven at production size.
* The order-12 Frobenius automorphism of GT could give a further ~√12 speedup
  (toward ~126) in this DP framework, but its fruitless cycles need empirical
  support and a reserved margin (as in other accepted Frobenius submissions on
  this track); that is the next step. A costed special-TNFS/SexTNFS ledger
  (~2^120–2^126 in those papers' own units) could be lower still.
* Analytic claim; no production-size execution.

References: Pollard (1978); van Oorschot–Wiener, J. Cryptology 1999; Teske,
Math. Comp. 2001.

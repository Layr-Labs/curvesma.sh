# Fixed-G1 pairing inversion via single-trail distinguished-point rho in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 127.85** (2557/20) target-group-operation
  equivalents (dense Fp12 multiplication at C = 20,000 scalar AArch64 instruction
  equivalents).
* Success probability: at least 0.39 under the random-mapping heuristic,
  supported by experiments.

This improves on my accepted O(1)-memory Floyd rho (score 129.2) by using a
**single-trail distinguished-point (DP) rho**: one continuous trajectory walks,
stores its distinguished points, and detects its own self-collision when a DP
element repeats. This spends **one dense Fp12 multiplication per step** (no
tortoise/hare 3× factor). Because there is a single trajectory, success is simply
P(rho-length ≤ budget) — the same clean form as the accepted Floyd submission —
with no cross-trail merge analysis. The DP store is distributed over 2^16 nodes
(each <2^64 bytes, routing charged), the partition/DP predicate reads the
reference multiplication's own reduced tower representation (no export), and
exponents come from the rules' fresh uniform algorithm coins.

## 1. Reduction

Input: nonidentity P in G1 (order r), z in GT. Compute g = e(P, Q0) (pinned
`pairing(Q0, P)`); nondegeneracy and prime order make g a generator of GT. Find a
with z = g^a (sections 2–4) and output Q = [a]Q0; e(P,Q) = g^a = z. If z = 1,
output infinity. P is never replaced.

## 2. Single-trail distinguished-point rho

R = 2^20 adding entries M_j = g^{alpha_j} z^{beta_j} with (alpha_j, beta_j) from
the rules' fresh uniform algorithm coins (512 coin bits reduced mod r, bias
< 2^-256, hard O(1), no PRNG). The walk is the deterministic map
F(X) = X·M_{H(X)}; a point X = g^A z^B carries (A,B).

**Partition and DP predicate (reduced tower representation, no export).** The
reference multiplication outputs reduced tower coefficients in [0,p); this
representation is unique per element. Let w be the integer formed by the two
least-significant 64-bit limbs (128 bits) of the first reduced tower coefficient
(a 381-bit value spanning six limbs, so bits 0..127 exist). H(X) = bits 0..19 of
w; X is **distinguished** when bits 20..83 of w are zero, so theta = 2^-64.
Reading w is two loads; H and the DP test are masks (≈ 8 instructions). No
canonical/flat-basis conversion is used in the walk; element equality compares
reduced tower representations directly.

**The trajectory.** A single trajectory X_0, X_1, ... = F^i(X_0) is run from a
random start (exponents from algorithm coins). At each DP it stores
(256-bit fingerprint of the reduced representation → (A,B)). Because F is a
function on the finite set GT, the trajectory is eventually periodic (tail mu,
cycle lambda; rho-length rho = mu + lambda). After the trajectory re-enters its
cycle it re-treads points it already visited; the **first DP it meets on the
re-tread was already stored**, so the DP element repeats. At that repeat,
X_i = X_j (i on the second pass, j on the first) are the same element with
different exponents, giving a = (A_j − A_i)/(B_i − B_j) mod r (B_i = B_j with
probability ~1/r, a counted failure). A 256-bit fingerprint over ≤ 2^65 stored
DPs has birthday false-match probability ≤ 2^(2·65 − 256) = 2^-126 (confirmed by
recomputation), negligible.

There is **no second trail and no cross-trail merge**: detection is the
trajectory meeting its own earlier stored DP. The only cycle is the genuine rho
cycle; a *distinguished-point-free* cycle (the whole cycle contains no DP) has
probability (1−theta)^{lambda} ≈ e^{−lambda·theta}, with lambda = Theta(sqrt r)
≈ 2^126 and theta = 2^-64 giving e^{−2^62}, utterly negligible.

## 3. Hard budget and detection

The processor performs at most **T_total = ceil(1.3·k0) + L_cap** multiplications
(k0 = sqrt(2 r ln(1/0.61)); L_cap = 20/theta = 20·2^64 ≈ 2^68.3) and then stops;
the budget is checked before every multiplication. The L_cap tail reserves the
detection interval: after the trajectory re-enters its cycle at step rho, the
next DP on the cycle (already stored) is reached within the gap to that DP, which
is ≤ L_cap (the gap between consecutive DPs exceeds L_cap only with probability
e^{−20} per gap, a negligible contribution bounded below). So a trajectory with
rho ≤ 1.3·k0 is detected by step rho + L_cap ≤ T_total. All work — walk steps,
the restart after any DP, and detection — is this single processor's, so T_total
is the true total (the score is independent of processor count; here it is one).

## 4. Success within the hard budget

Model the walk as a random mapping on the r elements. The rho-length satisfies
P(rho > t) ≈ exp(−t^2/(2r)). With T = ceil(1.3·k0),

    P(rho ≤ T) = 1 − exp(−(1.3·k0)^2/(2r)) = 1 − 0.61^{1.69} = 0.566282.

A trajectory with rho ≤ T is detected within T_total (section 3) unless (a) its
cycle is DP-free (≤ e^{−2^62}), (b) a gap between consecutive cycle DPs exceeds
L_cap (≤ (number of cycle DPs)·e^{−20} ≤ 2^63·2^-28.9 = 2^34.1·... bounded by the
step-weighted gap-overrun ≈ 20·e^{−20} ≈ 2^-24.7 of trajectories), or (c) the
repeat is unusable (B_i = B_j, ≤ 1/r), or (d) a false fingerprint match
(≤ 2^-126). Hence

    success ≥ 0.566 − 2^-24.7 − 2^-126 − 1/r ≥ 0.566 − 2^-24.6 ≥ **0.39**,

with wide margin. (No separate "production-transfer" allowance is used; the bound
is the clean random-mapping P(rho ≤ T), as in my accepted Floyd submission whose
walk is identical.)

### Experiment (toy scale, Python 3.12, one core)

`experiment.py` runs the single-trail DP rho (one trajectory, store DPs, detect
the self-collision DP-repeat) with the hard budget T_total = ceil(1.3·k0)+L_cap,
in the prime-order-ell subgroup of F_q*. It reports the fraction of trials
finding a usable self-collision within the budget.

| bits | seeds | success within budget | P(rho ≤ T) |
| ---: | --- | ---: | ---: |
| 24 | 1–5 | 0.572, 0.652, 0.620, 0.620, 0.608 | 0.566 |

Success meets or exceeds the 0.566 existence bound (there is no cross-trail
detection loss to pull it below), comfortably ≥ 0.39. Standard support for
R-adding walks as random mappings: Teske (2001).

## 5. Distributed storage and routing (no >2^64-byte address space)

The expected number of stored DPs is ≈ T·theta ≈ 2^63.8; a **hard cap
N_cap = 2^65** (≈ 2.3× the mean) halts the run as a counted failure if exceeded,
which happens only with probability exp(−Theta(2^63)) by Chernoff — so N_cap is a
hard storage bound. The ≤ 2^65 reports are held in **2^16 nodes**, each an
independent machine addressing < 2^64 bytes, with a hard per-node cap of 2^49
records of ≤ 2^10 bytes = **≤ 2^59 bytes per node** (routing by fingerprint-bits
[84..99]; per-node overflow → counted failure, Chernoff-negligible). Each DP
report is charged 2^17 instructions (fingerprint, network/stack, node-local tree
insert/lookup): total ≤ 2^65·2^17 = 2^82 instructions ≤ 2^67.7 units, negligible.
Node set-up is charged 2^20 instructions/node (2^36 total, negligible).

## 6. Cost, memory, advice

* Per step: one dense Fp12 multiplication (1 unit) + tower-rep H/DP masks (≈ 8) +
  two 255-bit exponent additions (≈ 160) + table load and control (≈ 100)
  ≈ 270 instructions = 0.0135 units. **Per step ≤ 1.014 units.**
* Multiplications ≤ T_total = 1.3·k0 + L_cap.

    W = (1.3·k0 + L_cap)·1.014 + report work (≤ 2^67.7) + one-time
    log2(1.3·k0·1.014) = 127.420269 + log2(1.3182) = 127.818814
    log2(W) = 127.818814 + log2(1 + 2^-59.5 + 2^-60.4 + ...) ≈ 127.818815.

One-time: table build (2^20 entries, 2 exponentiations + coin sampling) ≤ 2^31
units; pairing g ≤ 2^15; [a]Q0 ≤ 2^10; collision solve (one inverse mod r) ≤ 2^9;
output ≤ 2^-2. Total one-time ≤ 2^31.1 units, negligible. Submitted bound
**time_log2 = 127.85**, rounded up.

* Memory: ≤ N_cap·2^10 = 2^75 bytes across 2^16 nodes of ≤ 2^59 bytes each
  (consistent with section 5), plus the 2^20-entry adding table (~2^30 bytes) and
  O(1) trajectory state. memory_log2_bytes = 75 (reported only; no scalar
  contribution; access charged in section 5).
* Preprocessing (table build + node set-up) ≤ 2^31 units, included in W.
* Nonuniform advice: none (0 bytes).

## 7. Limitations and better directions

* Single modelling heuristic: the R-adding walk behaves like a random mapping on
  the r elements (Teske 2001; toy ~0.6, meeting the 0.566 bound). Not proven at
  production size.
* The order-12 Frobenius automorphism could give a further ~√12 speedup (toward
  ~126) in this DP framework, but its fruitless cycles need empirical support and
  a reserved margin (as in accepted Frobenius submissions on this track). A costed
  special-TNFS/SexTNFS ledger (~2^120–2^126 in those papers' own units) could be
  lower still.
* Analytic claim; no production-size execution.

References: Pollard (1978); van Oorschot–Wiener, J. Cryptology 1999; Teske,
Math. Comp. 2001.

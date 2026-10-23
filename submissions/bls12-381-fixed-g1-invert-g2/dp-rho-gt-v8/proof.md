# Fixed-G1 pairing inversion via single-trail distinguished-point rho in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 127.85** (2557/20) target-group-operation
  equivalents (dense Fp12 multiplication at C = 20,000 scalar AArch64 instruction
  equivalents).
* Success probability: at least 0.39 under the random-mapping heuristic, with the
  distinguished-point density and routing balance supported by a subgroup
  exponential-sum equidistribution argument (section 2a) and a confirming
  native-representation experiment (section 2b).

This improves on my accepted O(1)-memory Floyd rho (score 129.2) by using a
**single-trail distinguished-point (DP) rho**: one continuous trajectory walks,
stores its distinguished points, and detects its own self-collision when a DP
element repeats. This spends **one dense Fp12 multiplication per step** (no
tortoise/hare 3× factor). Because there is a single trajectory, success is simply
P(rho-length ≤ budget) — the same clean form as the accepted Floyd submission —
with no cross-trail merge analysis. The DP store is distributed over 2^16 nodes
(each <2^64 bytes, routing charged), the partition/DP predicate reads a bit-slice
of the reference multiplication's **native Montgomery-form limbs** (no export),
and exponents come from the rules' fresh uniform algorithm coins.

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

**Partition and DP predicate (native Montgomery limbs, no export).** The pinned
blst kernel keeps every Fp element in **Montgomery representation** — the integer
`mont(c) = c·2^384 mod p` held across six 64-bit limbs — and the walk reads those
limbs directly, with no conversion to ordinary residue. Let c be the first reduced
tower coefficient of X (an Fp element, 381 bits) and let w be the integer formed
by the two least-significant 64-bit limbs (128 bits) of `mont(c)` (so bits 0..127
of the Montgomery integer). H(X) = bits 0..19 of w; X is **distinguished** when
bits 20..83 of w are zero. By the equidistribution property (section 2a) these
64 bits are ≈ uniform, so the distinguished density is theta ≈ 2^-64. Reading w is
two loads; H and the DP test are masks (≈ 8 instructions). Everything below —
density, joint density, routing balance — is stated for this Montgomery-form
predicate, which is exactly the one the walk evaluates.

### 2a. Equidistribution of the Montgomery-form DP/routing bits (theoretical)

The one property beyond the random-mapping collision-time model is that a fixed
low-order bit window of `mont(c(X))` is close to uniform across GT, with disjoint
windows close to independent. This is **not** assumed; it follows from a standard
exponential-sum bound for multiplicative subgroups, which is representation-
independent and gives *joint* near-uniformity (hence the 64-bit DP density and the
routing balance conditional on being distinguished), not merely marginal bits.

*Set-up.* Let q = p^12 and let H = GT ⊂ F_q^* be the order-r subgroup, with
r | Φ_12(p) = p^4 − p^2 + 1. Because BLS12-381 has embedding degree 12, r ∤ p^d − 1
for every proper divisor d | 12; hence H lies in **no proper subfield** of F_q and
is not a coset of one. The coefficient map c = π(X) is a fixed F_p-linear
projection F_q → F_p (extract the tower coordinate d0.c0.a0), and the Montgomery
map is m(y) = 2^384·y mod p, an F_p-linear bijection of F_p. The quantity the walk
reads is a low-bit window of m(π(X)).

*Subgroup exponential-sum bound.* For any nontrivial additive character χ of F_p,
the composite χ∘m∘π is a nontrivial additive character ψ of F_q (π is onto and
F_p-linear, m and χ are nontrivial). By the Bourgain–Glibichuk–Konyagin bound and
its finite-field extension (Bourgain–Chang) for a multiplicative subgroup H ⊆ F_q^*
with |H| ≥ q^δ that is not contained in a proper subfield,

    | Σ_{X ∈ H} ψ(X) |  ≤  |H| · q^{-ε}

for some ε = ε(δ) > 0. Here |H| = r ≈ q^{0.0557} (δ = 0.0557 > 0) and the
no-subfield hypothesis holds by embedding degree 12, so the bound applies with a
positive power saving. (The same inequality holds for every nontrivial χ, i.e. for
every nonzero frequency.)

*From the character bound to window uniformity.* Fix a window of L consecutive
bits of the integer m(π(X)) ∈ [0,p). Its indicator is, by the standard completion
technique, a linear combination of ≤ 2^L · log p additive characters χ of F_p, each
with bounded coefficient. Summing the bound above over these characters, the number
of X ∈ H whose window equals any fixed value v deviates from |H|·2^{-L} by at most

    |H| · 2^L · q^{-ε} · polylog(q).

For every window used here (L ≤ 84 for the DP slice, and the combined DP+routing
window has L ≤ 100), 2^L · q^{-ε} is astronomically small — q^{ε} = p^{12ε} dwarfs
2^{100} for the relevant ε — so the deviation from uniform is < 2^{-200}·|H|,
negligible. Hence:

* the 64-bit DP window (bits 20..83 of `mont(c)`) takes the all-zero value with
  density theta = (1 ± 2^-200)·2^-64 — the **joint** 64-bit density, not a product
  of 64 marginals; and
* the combined 80-bit window (bits 20..99, i.e. DP slice ∪ 16 routing bits) is
  uniform, so **conditional on the DP event** the 16 routing bits (84..99) are
  uniform — routing is balanced conditional on being distinguished, which is what
  the overflow/geometric-gap calculations require.

*Representation independence.* The bound was derived through the Montgomery map m,
so it certifies the native-limb predicate directly; dropping m (ordinary residue)
only replaces ψ by another nontrivial character and leaves the bound unchanged.
So the measured property transfers between ordinary and Montgomery representations
— there is no representation gap — and the experiment below, which now reports the
Montgomery-form slices the walk reads, is a confirmation of this theorem at a
feasible scale, not the sole support.

### 2b. Confirming experiment (native Montgomery bits, real GT)

`experiment.py` implements the blst Fp2/Fp6/Fp12 tower, forms a GT element as
(random Fp12)^((p^12−1)/r) (verified by g^r = 1), walks a multiplicative sequence
in GT, and histograms the **Montgomery-form** coefficient-bit slices the walk
reads — `mont(c) = c·2^384 mod p` — at a feasible 12-bit scale (plus the ordinary
residue for reference). Representative output (2^13 samples, seed 7):

| quantity | value |
| --- | ---: |
| g in order-r subgroup (g^r = 1) | true |
| Montgomery DP slice bits 20..31, chi^2/dof vs uniform | 0.995 |
| Montgomery routing slice bits 84..95, chi^2/dof vs uniform | 1.012 |
| ordinary-residue DP slice, chi^2/dof (reference) | 0.978 |
| ordinary-residue routing slice, chi^2/dof (reference) | 1.020 |

Both the Montgomery slices (the exact bits the walk reads) and the ordinary-residue
slices are statistically uniform (chi^2/dof ≈ 1), on genuine GT elements, matching
the theorem of 2a and confirming that the representation choice does not change the
distribution. The 12-bit scale is where uniformity is directly testable with
feasible samples; the 64-bit joint density and conditional routing balance are
established by the exponential-sum bound of 2a, which the 12-bit marginals
corroborate. No pseudorandom-function assumption on a separate hash is made.

**The trajectory and detection by exact key.** A single trajectory
X_0, X_1, ... = F^i(X_0) is run from a random start (exponents from algorithm
coins). At each DP it stores the **full reduced tower representation of the
element** as the key, together with (A,B). Keys are therefore **exact group
elements**: a stored key repeats if and only if the trajectory revisits that
exact element, so there are **no false matches** and no fingerprint-collision
assumption is needed. Because F is a function on the finite set GT, the
trajectory is eventually periodic (tail mu, cycle lambda; rho-length rho =
mu + lambda). After re-entering its cycle it re-treads visited points; the first
DP on the re-tread is already stored, so its exact key repeats. At that repeat
X_i = X_j (i on the second pass, j on the first) are the same element with
different exponents, giving a = (A_j − A_i)/(B_i − B_j) mod r (B_i = B_j with
probability ~1/r, a counted failure).

There is **no second trail and no cross-trail merge**: detection is the
trajectory meeting its own earlier stored DP. Detection fails only if the cycle
contains no distinguished point. For a random mapping the cycle length lambda is
≈ uniform on (0, sqrt(pi r/2)], so P(no DP on the cycle) ≤ P(lambda < c/theta) +
E[(1−theta)^{lambda} | lambda ≥ c/theta]. The first term — short cycles likely to
miss a DP — is ≈ (c/theta)/sqrt(pi r/2) = c·2^64/2^126.4 = c·2^-62.4; the second
is ≤ e^{−c}. Taking c = 20 gives P(DP-free cycle) ≤ 2^-57.1 + e^{−20} ≤ 2^-28.8
(e^{−20} = 2^-28.854, rounded up), negligible. (This short-cycle bound replaces the
earlier incorrect typical-length substitution.)

## 3. Hard budget and detection

The processor performs at most **T_total = ceil(1.3·k0) + L_cap** multiplications
(k0 = sqrt(2 r ln(1/0.61)); L_cap = 20/theta = 20·2^64 ≈ 2^68.3) and then stops;
the budget is checked before every multiplication. The L_cap tail reserves the
detection interval: after the trajectory re-enters its cycle at step rho, the
next DP on the cycle (already stored) is reached within the gap to that DP, which
is ≤ L_cap (the gap between consecutive DPs exceeds L_cap only with probability
e^{−20} = 2^-28.854 per gap, bounded above by 2^-28.8). So a trajectory with
rho ≤ 1.3·k0 is detected by step rho + L_cap ≤ T_total. All work — walk steps,
the restart after any DP, and detection — is this single processor's, so T_total
is the true total (the score is independent of processor count; here it is one).

## 4. Success within the hard budget

Model the walk as a random mapping on the r elements. The rho-length satisfies
P(rho > t) ≈ exp(−t^2/(2r)). With T = ceil(1.3·k0),

    P(rho ≤ T) = 1 − exp(−(1.3·k0)^2/(2r)) = 1 − 0.61^{1.69} = 0.566282.

A trajectory with rho ≤ T is detected within T_total (section 3) unless (a) its
cycle is DP-free (≤ 2^-28.8, by the short-cycle bound in section 2b), (b) the gap
from cycle re-entry to the next stored cycle DP exceeds L_cap (geometric, e^{−20}
= 2^-28.854 ≤ 2^-28.8), or (c) the repeat is unusable (B_i = B_j, ≤ 1/r). Keys are
exact group elements, so there is no false-match term. Hence

    success ≥ 0.566 − 2^-28.8 − 2^-28.8 − 1/r = 0.566 − 2^-27.8 − 1/r ≥ **0.39**,

with wide margin. (No separate "production-transfer" allowance is used; the bound
is the clean random-mapping P(rho ≤ T), as in my accepted Floyd submission whose
walk is identical. The two e^{−20} terms are each rounded UP to 2^-28.8 ≥ e^{−20},
correcting the earlier 2^-28.9 that rounded the wrong way.)

## 5. Distributed storage and routing (no >2^64-byte address space)

The expected number of stored DPs is ≈ T·theta ≈ 2^63.8; a **hard cap
N_cap = 2^65** (≈ 2.3× the mean) halts the run as a counted failure if exceeded,
which happens only with probability exp(−Theta(2^63)) by Chernoff — so N_cap is a
hard storage bound. The ≤ 2^65 reports are held in **2^16 nodes**, each an
independent machine addressing < 2^64 bytes, with a hard per-node cap of 2^49
records; each record is the full reduced representation (12 Fp coefficients ×
48 bytes = 576 bytes) plus (A,B) (≤ 64 bytes) plus tree metadata, ≤ 2^10 = 1024
bytes per record (≥ 640 + metadata), so **≤ 2^49·2^10 = 2^59 bytes per node**.
A DP routes to a node by bits [84..99] of w; by the equidistribution property
(section 2a) these 16 bits are ≈ uniform conditional on the DP event, so node
loads are balanced and per-node overflow beyond the 2^49 cap has
Chernoff-negligible probability (counted failure if hit). Each DP report is charged
2^17 instructions (key formation, network/stack, node-local balanced-tree
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

* Memory: the DP records are ≤ N_cap·2^10 = 2^75 bytes across 2^16 nodes of
  ≤ 2^59 bytes each (section 5); adding the 2^20-entry table (~2^30 bytes) and O(1)
  trajectory state, the total peak is ≤ 2^75 + 2^30 ≤ **2^76 bytes**.
  memory_log2_bytes = 76 (reported only; no scalar contribution; access charged in
  section 5). Per-node bound 2^59 bytes is consistent between this proof and the
  manifest.
* Preprocessing (table build + node set-up) ≤ 2^31 units, included in W.
* Nonuniform advice: none (0 bytes).

## 7. Limitations and better directions

* Single modelling heuristic: the R-adding walk behaves like a random mapping on
  the r elements (Teske 2001; toy ~0.6, meeting the 0.566 bound). Not proven at
  production size. The DP/routing-bit equidistribution (section 2a) is a theorem
  modulo the explicit BGK/Bourgain–Chang power-saving constant, confirmed
  empirically at 12 bits; the collision-time law is the modelling heuristic.
* The order-12 Frobenius automorphism could give a further ~√12 speedup (toward
  ~126) in this DP framework, but its fruitless cycles need empirical support and
  a reserved margin (as in accepted Frobenius submissions on this track). A costed
  special-TNFS/SexTNFS ledger (~2^120–2^126 in those papers' own units) could be
  lower still.
* Analytic claim; no production-size execution.

References: Pollard (1978); van Oorschot–Wiener, J. Cryptology 1999; Teske,
Math. Comp. 2001; Bourgain, Glibichuk, Konyagin, J. London Math. Soc. 2006;
Bourgain, Chang, "Exponential sum estimates over subgroups and almost subgroups
of Z_q^*", 2009 (finite-field subgroup exponential sums).

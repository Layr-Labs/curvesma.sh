# STNFS for recover-Q: 2^120.760556358 reference units

Track bls12-381-fixed-g1-invert-g2. Unit: one dense F_p12 multiplication,
C = 20,000 scalar AArch64 instructions. No auxiliary powers, advice, smaller
targets or amortization; no prior verdict is a premise.

## Claim

The algorithm returns Q in G2 with e(P,Q) = z for the exact supplied
nonidentity P in G1 (Q = infinity for z = 1) with success at least 0.39
(modelled >= 0.393291) under heuristics H1-H5 (section 9). A global meter
bounds all work, preprocessing included, by **2^120.760556358 units**; the
score is **120.760556358** (unrounded ledger 2^120.76055635607, declared at the
1e-9 ceiling of ledger + 1e-9). Storage <= 2^105 bytes (STORAGE.md). Analytic
exploratory claim, not a full-size run. `evidence16.py` (quick, ~9 s) checks
every stated bound against its unrounded value; quick = false adds [6]-[12];
`extra` = relw / sensitivity runs [11u]/[11s] / [7s] (each < 20 s).

## 1. Reduction

g = e(P, Q0) (pinned `pairing(Q0, P)`), Q0 the G2 generator: P has order r and
e is nondegenerate, so g generates mu_r, z = g^a for a unique a mod r, and
Q = [a]Q0 satisfies e(P, Q) = z for that exact P (never replaced). log z and
log g share one relation and linear-algebra phase; a = log z / log g mod r. g
and z enter the tower of section 2 by the standard field isomorphism, charged
with the pairing, [a]Q0 and the final check e(P, Q) = z in the 2^66-unit row.

## 2. Tower and polynomials ([0], [0b], [0c], [13], [14])

h = y^6 + y^3 + y^2 - y - 1 (Disc 58673 squarefree, irreducible mod p);
f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3);
g = x^2 + |u| y, u = -0xd201000000010000; F_p[y,x]/(h, g) = F_p12. sigma:
x -> -x quotients candidates and unknowns. Murphy alpha charged (5.15, 0.85)
in collection and descent; computed (4.868096, 0.843643) and 0.731544 (g,
p < 50000).

## 3. Relation collection ([2]-[5], [4a], [12])

Pairs (a, b) in O_h^2, B = 2^69, U = 2^63.451 sigma-orbits. Region:
A^2 = 2^21, a5 > 0, Q_a(a) <= (9/25) A^2, b5 >= 1, Q_b(b) <= A^2 - Q_a(a); one
canonical a per unit class (exact HNF key); non-primitive pairs excluded.
Certified count products at a-radius >= 868: 0.971088853 (yield), 1.028002876
(work). Canonical list 2^90 units.

**Pool.** 40,000 uniform draws on theta = 19/50 (seed 33: 2,000; seed 34:
38,000; `pool12-*.b85`, [2] regenerates seed-33 rows); the region is its
uniform sub-region 9/25 (4,094 canonical rows, effective sample 34,010.8).

**Sieve.** Line sieve for p in [2^8, 2^48), record walk for [2^48, 2^64)
(KERNELS.md), ECM for [2^64, 2^69] within slacks (295, 224). Survivors: exact
norms, trial division below 2^20, two band-strip rounds per side at 99.9% (a
missed band prime loses the relation: yield x (1 - 41e-6)(1 - 28e-6)), then the
69-bit schedule at 99.9%. **Early abort:** the f side (band strip, ECM) runs
first; the g side only if f factors completely, charged per row at
q_f = 100 rho(u_f) / P(omitted_f <= slack_f) (mean 1.18e-3; H1 estimate x 100;
q_f = 1 would only restore the two-sided charge). Yield is unchanged.

**Yield and cost.** Retained relations 2^64.0641 (pool SE 0.0903 bits) against
1.2 U: margin >= 0.3497 bits (3.9 SE). An independent accepted pool at the same
region (our 57598408: theta 2881/8000 = 0.3601, 150,000 draws, same B, X2,
alpha, ECM and threshold model, slack_g 226) gives 2^63.9402, above the 63.883
descent input. Charged 317.4 instructions/cell over
2^126.558 cells, allowance **1.0332 >= 1 + 2 x 0.01660** (work-estimator SE, a
modelling margin): collection **<= 2^120.5807 units**. Descent inputs
LOG2_RELS = 63.883 = estimate - 2 SE ([4] asserts it) and relation-weighted
mean norms 743.2/498.9.

## 4. Linear algebra ([6], [10], [13])

Signed sigma-quotient filtering at excess 1.2 (`filtsim.py`) merged to mean
weight 200 (charged 201); N(200)/N(100) <= 0.806 over 33 runs (charged 0.81).
Charged N = 0.81 x 1.09 x 0.05167 x U < 2^58.998. A = [M | S], 21 Schirokauer
columns. Block Wiedemann (19, 4); per row 201 x 20.7 sparse (H4u) + 21 x 97
dense + row tail 1,640. LA **<= 2^117.5139 units**, generator <= 2^106.9115.

## 5-7. Descent and success ([7], [7r], DESCENT.md)

62-node plan and boot at the estimate - 2 SE input: expected two-log work
**<= 2^113.6317 units**, meter 2^114.36, so P(cap) <= 0.603613 (Markov). Boot
failure <= 3.4635e-07 per log; node counter 2^25 with E[N] <= 2^16.17 (abort
<= 2^-8.83); special-q repeats <= 2^-15.30; missing-log ideals are descended
in their own prime lattices L_Q and any repeat is a failure, <= 2^-10.16;
Block Wiedemann <= 2^-192 (H3); field map < 2^-48.
Success >= 1 - 2 x 3.4635e-07 - 0.603613 - 2^-8.83 - 2^-15.30 - 2^-10.16 -
2^-48 - 2^-192 - 2/r >= 0.393291.

## 8. Ledger (log2 units)

| row | bound |
|---|---:|
| relation collection | 120.5807 |
| linear algebra / generator | 117.5139 / 106.9115 |
| two-log meter | 114.36 |
| canonical list / Schirokauer / misc / factor base / certificates | 90 / 88 / 80 / 72 / 50 |
| field map, pairing, [a]Q0, verification, output | 66 |

Sum 2^120.76055635607 < **2^120.760556358**.

## 9. Heuristics

H1 smoothness (Dickman/Murphy), threshold retention, random ECM orders, f/g
independence given survival (early abort); H2 uniform pool and work allowance;
H3 matrix (class numbers, unit maps, kernel; N/U and the weight-200 ratio not
rising at production); H4 kernel counts and H4u; H5 descent (DESCENT.md):
lower sums, negative-binomial recurrence at the estimate - 2 SE relation input,
UDUP 0.65 (relation-weighted E[shat/m], uniform cells), norm-distribution pr
integration (caps 4x/1.5x), prime L_Q missing-log nodes, conditional
independence of subtree failures given the main-phase missing set M (repeats
charged; common-cause dependence bounded by the conservative M quantile),
exceptional bound-gap cells charged 2^-20.

## Provenance

Collection, LA and kernels: gbotrel's recP package 42e80f80 (sister track) and
its lineage (8cc93834, a31a5f1d, 6282ab79, 3a00ffc8, 8dc444f6); early abort:
gbotrel 2c2ce2ab; descent plan, pr integration, relation-weighted UDUP and
uniform sampler: our 57598408 lineage. This port (region, est - 2 SE input,
prime L_Q, repeat bound, certified roots, bound gap): Claude Opus 5.5 (max) in
Claude Code via T3 Code.

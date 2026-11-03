# STNFS for recover-P: 2^120.971000 reference units

Track bls12-381-fixed-g2-invert-g1. Unit: one dense F_p12 multiplication,
C = 20,000 scalar AArch64 instructions. No auxiliary powers, advice, smaller
targets or amortization; no prior verdict is a premise.

## Claim

The algorithm returns P in G1 with e(P,Q) = z for the exact supplied
nonidentity Q in G2 (P = infinity for z = 1) with success at least 0.39
(modelled >= 0.392726) under heuristics H1-H5 (section 9). A global meter
bounds all work, preprocessing included, by **2^120.971000 units**; the score
is **120.971000**. Storage <= 2^105 bytes (STORAGE.md). Analytic exploratory
claim, not a full-size run. `evidence16.py` (quick initial run) recomputes the
figures below and checks each stated bound (`declared_bounds`); [6], [6b], [11],
[11b], [12] run with quick=false.

## 1. Reduction

g = e(P0, Q) (pinned `pairing(Q, P0)`), P0 the G1 generator: g generates mu_r,
z = g^a for a unique a mod r, and P = [a]P0. log z and log g share one
relation and linear-algebra phase; a = log z / log g mod r. g and z enter the
tower of section 2 by the standard field isomorphism (roots of h and of
X^2 + |u| y in F_p12, then a 12 x 12 basis change), charged with the pairings,
[a]P0 and the final check e(P, Q) = z in the 2^66-unit row.

## 2. Tower and polynomials ([0], [0b], [0c], [13], [14])

h = y^6 + y^3 + y^2 - y - 1 (Disc 58673 squarefree, irreducible mod p);
f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3);
g = x^2 + |u| y, u = -0xd201000000010000; F_p[y,x]/(h, g) = F_p12. sigma:
x -> -x quotients candidates and unknowns. Murphy alpha charged (5, 0.85),
computed (4.868096, 0.843643) and 0.731544 for g over p < 50000; descent (5, 1).

## 3. Relation collection ([2]-[5], [4a], [12])

Pairs (a, b) in O_h^2, B = 2^69, U = 2^63.451 sigma-orbits. Region
(`config.json`): A^2 = 2^21, a5 > 0, Q_a(a) <= (19/50) A^2, b5 >= 1,
Q_b(b) <= A^2 - Q_a(a); one canonical a per unit class (exact HNF key);
non-primitive pairs excluded. Certified count products 0.971857578 (yield),
1.027240121 (work). Canonical list 2^90 units.

**Pool.** 40,000 uniform draws on this region (seed 33: 2,000; seed 34:
38,000), 4,804 canonical rows with exact norms (`pool12-*.b85`; [2]
regenerates seed-33 rows). The region was chosen on separate exploratory pools.

**Sieve.** Line sieve (scores only) for p in [2^8, 2^48), record walk for
[2^48, 2^64) (KERNELS.md; every line keeps >= 2^6.79 hits per coset, inside
the measured line table), ECM for [2^64, 2^69] within slacks (295, 224).
Survivors: exact norms, trial division below 2^20, two rounds per side of the
`ecm16.py` 48-bit schedule at 99.9% (band strip; a missed band prime loses the
relation: yield x (1 - 41e-6)(1 - 28e-6)), then the 69-bit schedule at 99.9%:
2^39.6498 instructions per survivor.

**Yield and cost.** Retained relations >= 2^64.1694 against 1.2 U: margin
>= 0.4550 bits (pool SE 0.085). Charged 330.5 instructions/cell over 2^126.706
cells, including allowance **1.0311 >= 1 + 2 x 0.01552** (work-estimator SE, a
modelling margin): collection **<= 2^120.7867 units**. Descent inputs
LOG2_RELS = 63.999 (estimate - 2 SE) and mean norms 744.3/498.9 are at or below the pool's.

## 4. Linear algebra ([6], [10], [13])

Signed sigma-quotient filtering at excess 1.2 (`filtsim.py`) merged to mean
weight 200 (charged 201). Embedded runs (`SCALE_NU`, `W200`, 33 runs at
2^12..2^18): weight-100 N/U maxima fall from 2^15 up (2^18: 0.0461), and
N(200)/N(100) <= 0.806. Charged N = 0.81 x 1.09 x 0.05167 x U < 2^58.998.
A = [M | S] with 21 Schirokauer columns (sigma-odd unit ranks 18 + 3).
Block Wiedemann, blocks (19, 4); per row and product 201 x 20.7 sparse (H4u:
>= 88% of weight is +-1, measured 0.917-0.948 at weight 200) + 21 x 97 dense +
row tail 1,640 (ten owners; STORAGE.md). LA **<= 2^117.5139 units**,
generator **<= 2^106.9115 units**.

## 5-7. Descent and success ([7], DESCENT.md)

DESCENT.md (relation input at estimate - 2 SE, calibrated missing logs [7c],
prime-lattice missing-log children, certified roots [7r]): expected two-log
work **<= 2^115.1254 units**, meter 2^115.851, P(cap) <= 0.604723 (Markov, H5).
Boot failure <= 1.5924e-05 per log; node counter 2^27.3 with E[N] <= 2^16.91
(abort <= 2^-10.39); repeats (special-q and missing-ideal) <= 2^-9.14; Block
Wiedemann <= 2^-192 under H3; field map < 2^-48. Success >= 1 - 2 x 1.5924e-05
- 0.604723 - 2^-10.39 - 2^-9.14 - 2^-48 - 2^-192 - 2/r >= 0.392726.

## 8. Ledger (log2 units)

| row | bound |
|---|---:|
| relation collection | 120.7867 |
| linear algebra / generator | 117.5139 / 106.9115 |
| two-log meter | 115.851 |
| canonical list / Schirokauer / misc / factor base / certificates | 90 / 88 / 80 / 72 / 50 |
| field map, pairings, verification, output | 66 |

Sum 2^120.97093 < **2^120.971000** (unrounded 2^120.970921).

## 9. Heuristics

H1 smoothness (Dickman/Murphy) and H1' threshold retention, including random
ECM group orders in collection; H2 uniform pool and work allowance; H3 matrix
(class numbers, unit maps, one-dimensional kernel; N/U and the weight-200
ratio not rising at production); H4 measured kernel counts and H4u +-1 share;
H5 descent model, including the missing-log calibration (DESCENT.md).

## Provenance

Core: Gautam Botrel with Claude Code and Codex (3cb69d49, 8dc444f6), with the
refinements of 6282ab79 / 3a00ffc8 / e3940d11 / 04073dd1 / 8c1635f2. Line
sieve to 2^48 with band strip and excess 1.2: 8cc93834 / a31a5f1d lineage.
Descent hardening: our sister-track v17p. This revision: Claude Opus 5.5 in
Claude Code.

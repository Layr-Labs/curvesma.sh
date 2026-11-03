# STNFS: 2^120.795000 reference units

Track bls12-381-fixed-g1-invert-g2. Unit: one dense F_p12 multiplication,
C = 20,000 scalar AArch64 instructions. No auxiliary powers, advice, smaller
targets or amortization; no prior verdict is a premise.

## Claim

The algorithm returns Q in G2 with e(P,Q) = z for the exact supplied
nonidentity P in G1 (Q = infinity for z = 1) with success at least 0.39
(modelled >= 0.391626) under heuristics H1-H5 (section 9). A global meter
bounds all work, preprocessing included, by **2^120.795000 units**; the score
is **120.795000**. Storage <= 2^105 bytes (STORAGE.md). Analytic exploratory
claim, not a full-size run. `evidence16.py` (initial run, quick mode) recomputes the pool, yield, cost,
linear algebra, descent, ledger and success and compares each stated bound with
its unrounded value (`declared_bounds`); the filtering rerun, missing logs,
retention and ECM checks ([6], [6b], [11], [11b], [12]) run with quick=false.

## 1. Reduction

g = e(P, Q0) (pinned `pairing(Q0, P)`), Q0 the G2 generator: g generates mu_r,
z = g^a for a unique a mod r, and Q = [a]Q0. log z and log g share one
relation and linear-algebra phase; a = log z / log g mod r. g and z enter the
tower of section 2 by the standard field isomorphism (roots of h and of
X^2 + |u| y in F_p12, then a 12 x 12 basis change), charged with the pairings,
[a]Q0 and the final check e(P, Q) = z in the 2^66-unit row.

## 2. Tower and polynomials ([0], [0b], [0c], [13], [14])

h = y^6 + y^3 + y^2 - y - 1 (Disc 58673 squarefree, irreducible mod p);
f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3);
g = x^2 + |u| y, u = -0xd201000000010000; F_p[y,x]/(h, g) = F_p12. sigma:
x -> -x quotients candidates and unknowns. Murphy alpha charged (5, 0.85),
computed (4.868096, 0.843643) and 0.731544 for g over p < 50000; descent (5, 1).

## 3. Relation collection ([2]-[5], [4a], [12])

Pairs (a, b) in O_h^2, B = 2^69, U = 2^63.451 sigma-orbits. Region
(`config.json`): A^2 = 2^21, a5 > 0, Q_a(a) <= (9/25) A^2, b5 >= 1,
Q_b(b) <= A^2 - Q_a(a); one canonical a per unit class (exact HNF key);
non-primitive pairs excluded. Certified count products at radius >= 868:
0.971088853 (yield), 1.028002876 (work). Canonical list 2^90 units.

**Pool.** 40,000 uniform draws on theta = 19/50 (seed 33: 2,000; seed 34:
38,000), 4,804 canonical rows with exact norms (`pool12-*.b85`; [2]
regenerates seed-33 rows); its 4,094 rows with qa <= 9/25 are uniform on the
region (effective sample 40000 x (360/380)^3 = 34010.79).

**Sieve.** Line sieve (scores only) for p in [2^8, 2^48), record walk for
[2^48, 2^64) (KERNELS.md; every line keeps >= 2^6.79 hits per coset, inside
the measured line table), ECM for [2^64, 2^69] within slacks (295, 224).
Survivors: exact norms, trial division below 2^20, then per side two rounds of
the `ecm16.py` 48-bit schedule at 99.9% (band strip; a missed band prime loses
the relation: yield x (1 - 41e-6)(1 - 28e-6)) and the 69-bit schedule at 99.9%.
The f side runs first; the g side only if f factors (charged with
q_f = 100 rho(u_f) / P(omitted_f <= slack_f), a heuristic estimate, mean
~1e-5; q_f = 1 restores the two-sided charge).

**Yield and cost.** Retained relations >= 2^64.0759 against 1.2 U: margin
>= 0.3616 bits (pool SE 0.090). Charged 317.4 instructions/cell over 2^126.558
cells, including allowance **1.0332 >= 1 + 2 x 0.01660** (work-estimator SE, a
modelling margin): collection **<= 2^120.5807 units**. Descent inputs
LOG2_RELS = 63.895 (estimate - 2 SE) and mean norms 743.2/498.8 are at or below
the pool's ([4]).

## 4. Linear algebra ([6], [10], [13])

Signed sigma-quotient filtering at excess 1.2 (`filtsim.py`) merged to mean
weight 200 (charged 201). Embedded runs (`SCALE_NU`, `W200`, 33 runs at
2^12..2^18): weight-100 N/U maxima are non-increasing from 2^15 (2^18: 0.0461; the charge also covers 2^14), and
N(200)/N(100) <= 0.806. Charged N = 0.81 x 1.09 x 0.05167 x U < 2^58.998.
A = [M | S] with 21 Schirokauer columns (sigma-odd unit ranks 18 + 3).
Block Wiedemann, blocks (19, 4); per row and product 201 x 20.7 sparse (H4u:
>= 88% of weight is +-1, measured 0.917-0.948 at weight 200, offline; a share of 0.80 would add < 0.01 bit to the total) + 21 x 97
dense + row tail 1,640 (ten owners; STORAGE.md). LA **<= 2^117.5139 units**,
generator **<= 2^106.9115 units**.

## 5-7. Descent and success ([7], DESCENT.md)

62-node plan and boot at LOG2_RELS = 63.895 with the calibrated missing-log model
and prime-lattice missing-log children (DESCENT.md): expected two-log work
**<= 2^115.1870 units**, meter 2^115.91, so P(cap) <= 0.6058149 (Markov; H5 is the
premise behind E[W]). Boot failure <= 2.6804e-05 per log; node counter 2^27.3 with
E[N] <= 2^17.02 (abort <= 2^-10.28); special-q and missing-ideal repeats <= 2^-9.20;
Block Wiedemann <= 2^-192
under H3; field map < 2^-48 (absorbed in the rounding of the floor). Success >= 1 - 2 x 2.6804e-05 - 0.6058149
- 2^-10.28 - 2^-9.20 - 2^-48 - 2^-192 - 2/r >= 0.391626. Sensitivity to the relation
count: DESCENT.md ([7s], full run).

## 8. Ledger (log2 units)

| row | bound |
|---|---:|
| relation collection | 120.5807 |
| linear algebra / generator | 117.5139 / 106.9115 |
| two-log meter | 115.91 |
| canonical list / Schirokauer / misc / factor base / certificates | 90 / 88 / 80 / 72 / 50 |
| field map, pairings, verification, output | 66 |

Sum 2^120.7932 < **2^120.795000** (unrounded 2^120.793115).

## 9. Heuristics

H1 smoothness (Dickman/Murphy), H1' threshold retention, including random ECM
group orders in collection and the early-abort estimate; H2 uniform pool and
work allowance; H3 matrix (class numbers, unit maps, one-dimensional kernel;
N/U and the weight-200 ratio not rising at production); H4 measured kernel
counts and H4u +-1 share; H5 descent model (DESCENT.md), including the missing-log calibration c(lambda) (measured
ratio maxima x 1.5 for lambda in [1, 8), full model elsewhere).

## Provenance

Gautam Botrel with Claude Code and Codex; refinements of 6282ab79 / 3a00ffc8 /
8cc93834 / a31a5f1d lineages; port and early abort: Claude Opus 5.5 in Claude Code.

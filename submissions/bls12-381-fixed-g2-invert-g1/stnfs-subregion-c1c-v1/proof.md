# STNFS for recover-P: 2^120.894095715 reference units

Track bls12-381-fixed-g2-invert-g1. Unit: one dense F_p12 multiplication,
C = 20,000 scalar AArch64 instructions. No auxiliary powers, advice, smaller
targets or amortization; no prior verdict is a premise.

## Claim

The algorithm returns P in G1 with e(P,Q) = z for the exact supplied
nonidentity Q in G2 (P = infinity for z = 1) with success at least 0.39
(modelled >= 0.410261) under heuristics H1-H5 (section 10). A global meter
bounds all work, preprocessing included, by **2^120.894095715 units**; the score
is **120.894095715** (= 120894095715 / 1000000000; the rows are declared on the 10^-9 grid with one grid step of rounding margin
each, their sum is 2^120.894095713, and the declared value is its ceiling on the
10^-9 grid plus one more grid step; the unrounded ledger is 2^120.894095712). Storage <= 2^105 bytes. Analytic
exploratory claim, not a full-size run. `evidence16.py` recomputes the model in six parts, each inside the hosted
limits: `core` (initial: region, yield, cost, linear algebra, the descent
recurrence, success, the ledger sum), `ledger`, `descent`, `lattices`,
`retention`, `losses` (submission.json lists their contents); every stated
bound is compared with its unrounded value (`declared_bounds`).

## 1. Reduction

g = e(P0, Q) (pinned `pairing(Q, P0)`), P0 the G1 generator: g generates mu_r,
z = g^a for a unique a mod r, and P = [a]P0. log z and log g share one
relation and linear-algebra phase; a = log z / log g mod r. g and z enter the
tower of section 2 by the standard field isomorphism (roots of h and of
X^2 + |u| y in F_p12, then a 12 x 12 basis change), charged with the pairings,
[a]P0 and the final check e(P, Q) = z in the 2^66-unit row (section 9).

## 2. Tower and polynomials ([0], [0b], [0c], [13], [14])

h = y^6 + y^3 + y^2 - y - 1 (Disc 58673 squarefree, irreducible mod p);
f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3);
g = x^2 + |u| y, u = -0xd201000000010000; F_p[y,x]/(h, g) = F_p12. sigma:
x -> -x quotients candidates and unknowns. Murphy alpha charged (5, 0.85),
computed (4.868096, 0.843643) and 0.731544 for g over p < 50000; descent (5, 1).

## 3. Relation collection ([2]-[5], [4a], [12])

Pairs (a, b) in O_h^2, B = 2^69, U = 2^63.451 sigma-orbits. Pool region
(`config.json`): A^2 = 2^21, a5 > 0, Q_a(a) <= (19/50) A^2, b5 >= 1,
Q_b(b) <= A^2 - Q_a(a); one canonical a per unit class (the Q_a-least
associate, decided per a by `ideal_canonical.canonical_a`); non-primitive pairs
excluded.

**Collected sub-region.** Collection runs on the uniform sub-region
Q_a(a) <= (363/1000) A^2 (= 726/2000) of the pool region (a-radius 872;
certified count products at every radius >= 872: 0.971219877 (yield),
1.027872797 (work)). The radius is the value at which 42e80f80's exploratory pools
(about 470,000 draws, documented in its public note before this package was
built) put the excess-1.2 margin at the 0.30 floor; the shipped pool was
evaluated there once. The
sub-region's rows are a uniform sample of it (4,203 canonical rows of 4,804,
effective sample 34,868.2 of the 40,000 draws), the method of f34f63d9. The same
evaluation at other radii (monotone; no local optimum was selected):

| a-radius theta x 2000 | rows | margin (bits) | relation collection (log2 units) |
|---:|---:|---:|---:|
| 760 | 4203 | 0.378 | 120.656 |
| 726 | 4203 | 0.377 | 120.657 |
| 700 | 3758 | 0.315 | 120.549 |

**Pool.** 40,000 uniform draws on the pool region (seed 33: 2,000; seed 34:
38,000), 4,804 canonical rows with exact norms (`pool12-*.b85`; [2]
regenerates seed-33 rows); drawn by 42e80f80, whose note records that the pool
region was chosen on separate exploratory pools.

**Sieve.** Line sieve (scores only) for p in [2^8, 2^48), record walk for
[2^48, 2^64) (KERNELS.md; every line keeps >= 2^6.79 hits per coset, inside
the measured line table), ECM for [2^64, 2^69] within slacks (295, 224).
Survivors: exact norms, trial division below 2^20, two rounds per side of the
`ecm16.py` 48-bit schedule at 99.9% (band strip; a missed band prime loses the
relation: yield x (1 - 41e-6)(1 - 28e-6)), then the 69-bit schedule at 99.9%:
2^39.6026 instructions per survivor.

**Yield and cost.** Retained relations: estimate 2^64.0913 (pool standard
error 0.0894 bits) against 1.2 U = 2^63.7144: margin >= 0.3769 bits
(>= 0.1981 at the estimate minus two standard errors). Charged
329.1 instructions/cell over 2^126.583 cells, including allowance
**1.0333 >= 1 + 2 x 0.01662** (incidence-inclusive relative standard error of
the sampled work estimator, a modelling margin): collection **<= 2^120.657029242
units**. Descent inputs: LOG2_RELS = 63.91 <= estimate - 2 SE = 63.9125
(the two-standard-error allowance is propagated through the missing-log rates,
the descent work and the success bound, DESCENT.md "Inputs"); mean norms
743.4/498.9 <= the sub-region's relation-weighted means 743.45/498.92
(cell-weighted 763.93/504.22).

## 4. Linear algebra ([6], [10], [13])

Signed sigma-quotient filtering at excess 1.2 (`filtsim.py`) merged to mean
weight 200 (charged 201). Weight-100 size: the basis is 16 runs at scale 2^19
at relation-weighted means 744.84/498.97 (>= this package's inputs, a harder
system; `FILT_TABLE`, caab2f05's campaign): worst N/U = 1/21.579 = 0.04634
(seed 15), charged **1.35 x 0.04634 = 0.06256**: above every weight-100 N/U
measured at scale 2^14 and above in both data sets (`SCALE_NU`, 42e80f80's 33
runs at 2^12..2^18, maxima 0.0554 / 0.0563 / 0.0545 / 0.0517 / 0.0461 at
2^14..2^18; `SCALE_NU2`, this line's 60 runs at 2^15..2^20, maxima 0.0602 (n =
30) / 0.0533 / 0.0516 / 0.0465 / 0.0466 (n = 16) / 0.0445 (n = 6)); the binding
value is the 2^15 tail (0.06024, 3.9% below the charge) and the means fall with
scale (0.054 at 2^15, 0.045 at 2^19, 0.043 at 2^20). The merge ratio
N(200)/N(100) <= 0.806 over 42e80f80's 33 runs is charged 0.81 (`W200`).
Charged N = 0.81 x 0.06256 x U < 2^59.1488 (nine occupied 2^56-element blocks
of the ten-owner layout, STORAGE.md). A = [M | S] with 21 Schirokauer columns
(sigma-odd unit ranks 18 + 3). Block Wiedemann, blocks (19, 4); per row and
product 201 x 20.7 sparse (H4u: >= 88% of weight is +-1, measured 0.917-0.948
at weight 200 by 42e80f80) + 21 x 97 dense + row tail 1,640 (ten owners;
STORAGE.md). LA **<= 2^117.817018072 units**, generator **<= 2^107.156045800 units**.

## 5-7. Descent and success ([7], [7c], [7r], [7s], DESCENT.md)

62-node plan and the sqrt(U)-weighted boot at UDUP = 0.25 (the relation-weighted
charge of 41025889, below every measured lattice), missing logs charged at the
never-occurring mass plus the measured unknown fraction of occurring ideals,
with shared-missing-set bounds (DESCENT.md), prime-lattice missing-log children
and the inputs of section 3: expected two-log work **<= 2^115.1915 units**, meter 2^115.971,
so P(cap) <= 0.582545 (Markov on the unrounded E[W]; H5 is the premise
behind E[W], with its three steps stated and supported separately in
DESCENT.md "Recurrence and success"). Boot failure <= 0.000417594 per log
(16 independent relations from 16 lattices, shared lattice <= 2^-25.39);
node counter 2^26.5 with E[N] <= 2^17.88 (abort <= 2^-8.62); special-q
repeats <= 2^-12.30; missing-ideal repeats under the shared missing set <=
2^-8.11; zero-occurrence concentration 2^-30 (Chebyshev 2^-42.86); Block
Wiedemann <= 2^-192 under H3; field map < 2^-48.
Success >= 1 - 2 x 0.000417594 - 0.582545 - 2^-8.62 - 2^-12.30 - 2^-8.11
- 2^-25.39 - 2^-30 - 2^-48 - 2^-192 - 2/r >= 0.410261 (the meter is set so
that the modelled success is at least 0.41, two hundredths above the required
0.39: the expected two-log work may exceed its model value by 0.049 bits
before the bound reaches 0.39; DESCENT.md "Headroom and sensitivity", [7s]).

## 8. Ledger (log2 units)

| row | bound |
|---|---:|
| relation collection | 120.657029242 |
| linear algebra / generator | 117.817018072 / 107.156045800 |
| two-log meter | 115.971 |
| canonical list / Schirokauer / misc / factor base / certificates (section 9) | 90 / 88 / 80 / 72 / 50 |
| field map, pairings, verification, output | 66 |

Sum of the declared rows 2^120.894095713 (unrounded rows 2^120.894095712) < **2^120.894095715**.

## 9. Side budgets (derived)

Each row is a fixed power of two in `model16.total_log2`; the derivations give
explicit bounds well below the rows (schoolbook counts where no kernel exists).

- **Canonical list, 2^90 units = 2^104.3 instructions.** Canonicity is decided
  per eligible a with no global list: `canonical_a` LLL-reduces the rank-6
  lattice a O_h (6 x 6, entries < 2^11, exact 128-bit arithmetic: < 2^16
  instructions), then enumerates the associates v with Q_a(v) <= Q_a(a) by
  Fincke-Pohst under a hard budget of 100,000 nodes (< 2^8 instructions each:
  < 2^25 per a); an a that exhausts the budget is dropped as non-canonical (its
  cells' relations were never credited; it costs at most the budget). With
  2^63.464 eligible a in the sub-region ([3]): < 2^63.464 x 2^25 < 2^89
  instructions = 2^74.7 units.
- **Schirokauer maps, 2^88 units = 2^102.3 instructions.** For each retained
  relation (< 2^64.2) and each side, (a - b x)^E mod (r^2, h, S) with E < r^72
  (K_f) or r^12 (K_g): < 18,400 squarings and multiplications of 72-coefficient
  polynomials with 510-bit coefficients (< 2^20 instructions each with the
  reductions), < 2^35.2 per relation over both sides: < 2^99.4 instructions =
  2^85.1 units; deterministic, no failure path.
- **Filtering and miscellaneous, 2^80 units = 2^94.3 instructions.** Relations
  carry at most 64 factor-base entries; singleton and clique removal and the
  merge to mean weight 200 touch < 2^70.2 entries per pass at < 2^8
  instructions for < 2^6 passes (< 2^84.2 instructions); back-substitution is
  one more pass; the Schirokauer, virtual-logarithm and ECM bookkeeping (failed
  cofactorisations are charged per cell and lost in the yield) stay below 2^80;
  no duplicate check is needed (one point per <units, sigma> orbit). Total
  < 2^85 instructions = 2^70.7 units.
- **Factor base, 2^72 units = 2^86.3 instructions.** The primes p < 2^69
  (< 2^63.4; segmented sieve < 2^72 instructions), per prime the roots of h
  mod p (< 2^13) and per root the roots of f(s, x) and g(s, x) (< 2^14): < 2^15
  per prime, < 2^78.4 in all = 2^64.1 units.
- **Field map, pairings, verification, output, 2^66 units = 2^80.3
  instructions.** Roots of h and of X^2 + |u| y in F_p12 (< 2^14 F_p12
  multiplications), the basis change, two pairings (< 2^14), [a]P0 and the
  check e(P, Q) = z: < 2^20 F_p12 multiplications = 2^20 units; charged 2^66 to
  cover the input checks and output on every processor.
- **Certificates, 2^50 units.** The one-time exact computations of
  `lattice_count_bounds.py` ([3]) and the tower checks ([0]-[0c], [13]): < 2^40
  instructions (the descent's lattice certificates are inside its own rows).

## 10. Heuristics

H1 smoothness (Dickman/Murphy) and H1' threshold retention, including random
ECM group orders in collection; H2 uniform pool and work allowance; H3 matrix
(class numbers, unit maps, one-dimensional kernel; N/U and the weight-200
ratio not rising at production); H4 measured kernel counts and H4u +-1 share;
H5 descent model (DESCENT.md: the Poisson occurrence model with the measured
unknown fraction of occurring ideals (H5c), the negative-binomial relation counts, the independence of
distinct cells' factorisations within a node (H5a), the measured unit-duplicate
charge, and the worst-sextile norm bins).

## Provenance

Core: Gautam Botrel with Claude Code and Codex (3cb69d49, 8dc444f6, 42e80f80)
with the refinements of 6282ab79 / 3a00ffc8 / e3940d11 / 04073dd1 / 8c1635f2;
line sieve to 2^48 with band strip and excess 1.2: 8cc93834 / a31a5f1d.
Retention and ECM modules, kernel charges 850/2,000, UDUP 0.25 with its
relation-weighted evidence, the sub-region method and the plan-search start:
f34f63d9 / 41025889 / 1a449d81 (kaankolcu line); filtering campaign at
2^15..2^20: caab2f05. Prime-lattice missing-log children, the missing-ideal
repeat term, the calibration method and the root enclosures: d0b0f071
(gbotrel), re-derived here at this package's inputs with new simulations.
The sqrt(U) weighting of the boot basis
(Guillevic's subfield initial splitting, ASIACRYPT 2015) was measured by a
Claude Opus 5.5 literature agent of the crew. This revision (sub-region at the
published crossing, two-standard-error descent input, missing-log model with
shared-set bounds, weighted boot, H5 argument, plan re-search, side-budget
derivations, parts): Claude Fable 5.1 in Claude Code via T3 Code.

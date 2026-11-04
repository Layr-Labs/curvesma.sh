# STNFS for recover-P: 2^120.7130837288 reference units

Track bls12-381-fixed-g2-invert-g1. Unit: one dense F_p12 multiplication,
C = 20,000 scalar AArch64 instructions. No auxiliary powers, advice, smaller
targets or amortization; no prior verdict is a premise.

## Claim

The algorithm returns P in G1 with e(P,Q) = z for the exact supplied
nonidentity Q in G2 (P = infinity for z = 1) with success at least 0.39
under heuristics H1-H5 (section 10); the model's value is 0.407781 (section
7), and the claim is the 0.39 of the rules, not that figure. A global meter
bounds all work, preprocessing included, by **2^120.7130837288 units**; the score
is **120.7130837288** (= 120713083730 / 1000000000; the rows are declared on the 10^-9 grid with one grid step of rounding margin
each, their sum is 2^120.713083729, and the declared value is its ceiling on the
10^-9 grid plus one more grid step; the unrounded ledger is 2^120.713083727). Storage <= 2^105 bytes. Analytic
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
x -> -x quotients candidates and unknowns. Murphy alpha charged (5.15, 0.85):
computed (4.868096, 0.843643) over p <= 2000 and (5.085036, 0.731544) over
p < 50000 ([14] checks both sides against the larger sums; 42e80f80 and our
273b1519 charged 5.0 on the f side, below the p < 50000 value); descent (5.15, 1).

## 3. Relation collection ([2]-[5], [4a], [12])

Pairs (a, b) in O_h^2, B = 2^69, U = 2^63.451 sigma-orbits. Pool region
(`config.json`): A^2 = 2^21, a5 > 0, Q_a(a) <= (19/50) A^2, b5 >= 1,
Q_b(b) <= A^2 - Q_a(a); one canonical a per unit class (the Q_a-least
associate, decided per a by `ideal_canonical.canonical_a`); non-primitive pairs
excluded.

**Collected sub-region.** Collection runs on the uniform sub-region
Q_a(a) <= (653/2000) A^2 of the pool region (a-radius 827; certified
count products at every radius >= 827: 0.969673713 (yield), 1.029409690 (work)):
the smallest radius on the 1/2000 grid at which the shipped pool's margin over
the excess-1.1 requirement is at least the 0.30-bit floor (the same selection
rule as 8cc93834 / 41025889 / f34f63d9 used at excess 1.2). The floor is
2.9 pool standard errors at this radius, and the selection runs over
nested subsets of one sample (adjacent radii share almost all rows), so it is
not a pick among independent pools and the descent input is in any case the
estimate minus two standard errors. The sub-region's rows are a uniform sample
of the pool (3,055 canonical rows of 4,804, effective sample 25,372.2 of the
40,000 draws). The same evaluation at other radii (margins increase with radius):

| a-radius theta x 2000 | rows | margin (bits) | relation collection (log2 units) |
|---:|---:|---:|---:|
| 726 | 3055 | 0.304 | 120.346 |
| 680 | 3055 | 0.302 | 120.347 |
| 649 | 2991 | 0.289 | 120.327 |
| 640 | 2857 | 0.263 | 120.281 |

**Pool.** 40,000 uniform draws on the pool region (seed 33: 2,000; seed 34:
38,000), 4,804 canonical rows drawn, of which `pool12-*.b85` ships the 3,055
sub-region rows with exact norms ([2] regenerates seed-33 rows); drawn by 42e80f80,
whose note records that the pool region was chosen on separate exploratory pools.

**Sieve.** Line sieve (scores only) for p in [2^8, 2^48), record walk for
[2^48, 2^64) (KERNELS.md; every line keeps >= 2^6.79 hits per coset, inside
the measured line table), ECM for [2^64, 2^69] within slacks (295, 224).
Survivors: exact norms, trial division below 2^20, two rounds per side of the
`ecm16.py` 48-bit schedule at 99.9% (band strip; a missed band prime loses the
relation: yield x (1 - 41e-6)(1 - 28e-6)), then the 69-bit schedule at 99.9%:
2^39.6026 instructions per survivor.

**Yield and cost.** Retained relations: estimate 2^63.8903 (pool standard
error 0.1014 bits) against 1.1 U = 2^63.5888: margin >= 0.3014 bits
(>= 0.0986 at the estimate minus two standard errors). Charged
329.1 instructions/cell over 2^126.273 cells, including allowance
**1.0393 >= 1 + 2 x 0.01961** (incidence-inclusive relative standard error of
the sampled work estimator, a modelling margin): collection **<= 2^120.348119967
units**. Descent inputs: LOG2_RELS = 63.68 <= estimate - 2 SE = 63.6874
(the two-standard-error allowance is propagated through the missing-log rates,
the descent work and the success bound, DESCENT.md "Inputs"); mean norms
741.2/498.7 <= the sub-region's relation-weighted means 741.25/498.75
(cell-weighted 762.95/504.47).

## 4. Linear algebra ([6], [10], [13])

Signed sigma-quotient filtering at excess 1.1 (`filtsim.py`; 42e80f80 and our
line used 1.2) merged to mean weight 200 (charged 201). Weight-100 size: the
basis is 16 runs at scale 2^19 at excess 1.1 and the sub-region's
relation-weighted means (`FILT_TABLE`): worst N/U = 1/18.480 = 0.05411
(seed 11), charged **1.35 x 0.05411 = 0.07305**. The pad is
the smallest value on the 0.05 grid whose charge exceeds every weight-100 N/U
measured at excess 1.1 at scale 2^14 and above by at least 3% (the rule of
caab2f05's pad 1.35 at excess 1.2): `SCALE_NU11`, 57 runs at 2^15..2^20
(maxima 0.0706 (2^15, n = 30) / 0.0638 (2^16, n = 4) / 0.0587 (2^17, n = 2) / 0.0540 (2^18, n = 2) / 0.0541 (2^19, n = 16) / 0.0521 (2^20, n = 3); binding value 0.07056, 3.5% below the charge); the
means fall with scale. For comparison the excess-1.2
bases of caab2f05 (2^19 worst 0.04634) and 42e80f80 (2^17 worst 0.05167) are in
the note's comparison table. The merge ratio N(200)/N(100) at excess 1.1 is
<= 0.7827 over 15 runs at 2^15..2^17 (`W200_11`; 42e80f80's excess-1.2
runs: <= 0.806), charged 0.80. Charged N = 0.80 x 0.07305 x U
< 2^59.3545 (11 occupied 2^56-element blocks; eleven-owner
layout, STORAGE.md). **Occ>=3 purge losses.** An ideal in 3 rows
each with a distinct singleton (3 equations, 4 unknowns: the f39a3f63 obstruction) is
unrecoverable post-hoc (information-theoretic); it is counted as unknown in the
missing-log calibration, whose covers bound the occurring remainder including it
(H5c premise, DESCENT.md; no recovery procedure and no domination term).
A = [M | S] with 21 Schirokauer columns (sigma-odd unit ranks 18 + 3). Block
Wiedemann, blocks (19, 4); per row and product 201 x 20.7 sparse (H4u: >= 88%
of weight is +-1, measured 0.917-0.948 at weight 200 by 42e80f80) + 21 x 97
dense + row tail 1,750 (eleven owners; STORAGE.md). LA **<= 2^118.245155523 units**,
generator **<= 2^107.487297801 units**.

## 5-7. Descent and success ([7], [7c], [7r], [7s], DESCENT.md)

62-node plan and the sqrt(U)-weighted boot at UDUP = 0.25 (the relation-weighted
charge of 41025889, below every measured lattice), missing logs charged at the
never-occurring mass in full plus the measured unknown fraction of occurring
ideals (H5c covers, DESCENT.md),
with shared-missing-set bounds (DESCENT.md), prime-lattice missing-log children
and the inputs of section 3: expected two-log work **<= 2^115.3762 units**, meter 2^116.162,
so P(cap) <= 0.580025 (Markov on the unrounded E[W]; H5 is the premise
behind E[W], with its three steps stated and supported separately in
DESCENT.md "Recurrence and success"). Boot failure <= 0.001763641 per log
(16 independent relations from 16 lattices, shared lattice <= 2^-26.15);
node counter 2^25.5 with E[N] <= 2^18.12 (abort <= 2^-7.38); special-q
repeats <= 2^-14.30; missing-ideal repeats under the shared missing set <=
2^-8.58; zero-occurrence concentration 2^-30 (Chebyshev 2^-41.32); Block
Wiedemann <= 2^-192 under H3 and the main-phase factor >= 1 - 2^-192 (2^(log2_N+2-254)
with log2_N < 60), combined <= 2^-191 by the union bound; field map < 2^-48.
Success >= 1 - 2 x 0.001763641 - 0.580025 - 2^-7.38 - 2^-14.30 - 2^-8.58
- 2^-26.15 - 2^-30 - 2^-48 - 2^-191 - 2/r >= 0.407781 (the meter is set so
that the model's success is two hundredths above the claimed
0.39: the expected two-log work may exceed its model value by 0.043 bits
before the bound reaches 0.39; DESCENT.md "Headroom and sensitivity", [7s]).

## 8. Ledger (log2 units)

| row | bound |
|---|---:|
| relation collection | 120.348119967 |
| linear algebra / generator | 118.245155523 / 107.487297801 |
| two-log meter | 116.162 |
| canonical list / Schirokauer / misc / factor base / certificates (section 9) | 90 / 88 / 80 / 72 / 50 |
| field map, pairings, verification, output | 66 |

Sum of the declared rows 2^120.713083729 (unrounded rows 2^120.713083727) < **2^120.7130837288**.

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
  2^63.006 eligible a in the sub-region ([3]): < 2^63.006 x 2^25 < 2^89
  instructions = 2^74.7 units.
- **Schirokauer maps, 2^88 units = 2^102.3 instructions.** For each retained
  relation (< 2^64.2) and each side, (a - b x)^E mod (r^2, h, S) with E < r^72
  (K_f) or r^12 (K_g): < 18,400 squarings and multiplications of 72-coefficient
  polynomials with 510-bit coefficients (< 2^20 instructions each with the
  reductions), < 2^35.2 per relation over both sides (< 2^64.2 retained relations):
  < 2^99.4 instructions = 2^85.2 units; deterministic, no failure path.
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
unknown fraction of occurring ideals under the H5c cover premise (section 4),
the negative-binomial relation
counts, the independence of distinct cells' factorisations within a node (H5a), the
measured unit-duplicate charge, and the worst-sextile norm bins).

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
Claude Opus 5.5 literature agent of the crew. Filtering excess 1.1 with its own campaign,
calibration and merge-ratio measurements, the smaller sub-region, the
eleven-owner row tail, the plan repair; everything else as in 273b1519 /
7c1acb81: sub-region method, two-standard-error descent input, missing-log
model with shared-set bounds, weighted boot, H5 argument, side-budget
derivations, parts: Claude Fable 5.1 in Claude Code via T3 Code. This revision
(reverts the missing-log charge to the calibration-only form under the H5c cover
premise with no domination term, re-meters the descent at 116.162, samples the
witness integer balls uniformly, and re-verifies the carried plan at the shipped
inputs): muse-spark-1.3-contributor in Muse Code via T3 Code.

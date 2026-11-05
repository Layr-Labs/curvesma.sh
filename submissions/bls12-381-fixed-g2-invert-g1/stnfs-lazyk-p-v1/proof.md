# STNFS for recover-P: 2^120.626057131 reference units

Track bls12-381-fixed-g2-invert-g1. Unit: one dense F_p12 multiplication,
C = 20,000 scalar AArch64 instructions. No auxiliary powers, advice, smaller
targets or amortization; no prior verdict is a premise.

## Claim

The algorithm returns P in G1 with e(P,Q) = z for the exact supplied
nonidentity Q in G2 (P = infinity for z = 1) with success at least 0.39
under heuristics H1-H5 (section 10); the model's value is 0.410508 (section
7), and the claim is the 0.39 of the rules, not that figure. A global meter
bounds all work, preprocessing included, by **2^120.626057131 units**; the score
is **120.626057131** (= 120626057131 / 1000000000; the rows are declared on the
10^-9 grid with one grid step of rounding margin each, their sum is
2^120.626057130, and the declared value is its ceiling on the 10^-9 grid plus
one more grid step; the unrounded ledger is
2^120.626057129). Storage <= 2^105 bytes. Analytic exploratory claim, not a
full-size run. `evidence16.py` recomputes the model in seven parts, each inside
the hosted limits: `core` (initial: region, yield, cost, linear algebra, the
descent recurrence, success, the ledger sum), `ledger`, `descent`, `lattices`,
`retention`, `losses`, `conc` (the realized missing-log table; submission.json
lists their contents); every stated bound is compared with its unrounded value
(`declared_bounds`).

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
p < 50000 ([14] checks both sides against the larger sums; 42e80f80 and
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
the smallest radius on the 1/2000 grid whose margin over the excess-1.1
requirement is at least the 0.30-bit floor (the rule of 8cc93834 / 41025889 /
f34f63d9; 2.9 pool standard errors; nested subsets of one sample, not a pick
among independent pools). The sub-region's rows are a uniform sample of the pool (3,055 canonical rows of 4,804, effective sample
25,372.2 of the 40,000 draws). The same evaluation at smaller radii (the pool
file ships the sub-region rows only):

| a-radius theta x 2000 | rows | margin (bits) | relation collection (log2 units) |
|---:|---:|---:|---:|
| 653 | 3055 | 0.301 | 120.292 |
| 652 | 3036 | 0.297 | 120.285 |
| 649 | 2991 | 0.288 | 120.269 |
| 640 | 2857 | 0.263 | 120.222 |

**Pool.** 40,000 uniform draws (seed 33: 2,000; seed 34: 38,000), 4,804
canonical rows drawn, of which `pool12-*.b85` ships the sub-region rows with
exact norms ([2] regenerates seed-33 rows).

**Sieve.** Line sieve (scores only) for p in [2^8, 2^48), record walk for
[2^48, 2^64) (KERNELS.md; every line keeps >= 2^6.79 hits per coset, inside
the measured line table), ECM for [2^64, 2^69] within slacks (295, 224).
Survivors: exact norms, trial division below 2^20, then **early abort** (the
rule of 90d99567): the f side first (two rounds of the `ecm16.py` 48-bit band
strip, then the 69-bit schedule, both at 99.9%); the g-side band strip and ECM
run only if the f side factors completely. The yield is unchanged; per
survivor the g-side work is charged with probability q_f = min(1, 100
rho(u_f)/P(omitted_f <= slack_f)) per pool row, u_f = (n_f + alpha_f)/69 (H1;
factor 100 safety; charged mean q_f 1.2e-03; q_f = 1 restores the two-sided
charge). A missed band prime loses the relation: yield x (1 - 41e-6)(1 -
28e-6). Mean charge per survivor 2^38.7520 instructions.

**Yield and cost.** Retained relations: estimate 2^63.8902 (pool standard
error 0.1014 bits) against 1.1 U = 2^63.5888: margin >= 0.3014 bits.
Charged 316.4 instructions/cell over 2^126.273 cells, including allowance
**1.0386 >= 1 + 2 x 0.01930** (incidence-inclusive relative standard error of
the sampled work estimator, a modelling margin): collection
**<= 2^120.291319028 units**. Descent inputs: LOG2_RELS = 63.68 <= estimate
- 2 SE = 63.6874 (the two-standard-error allowance is propagated through the
missing-log rates, the descent work and the success bound, DESCENT.md
"Inputs"); mean norms 741.2/498.7 <= the sub-region's relation-weighted means
741.25/498.75.

## 4. Linear algebra ([6], [10], [13])

Signed sigma-quotient filtering at excess 1.1 (`filtsim.py`; 42e80f80 used
1.2) merged to mean weight 200 (charged 201). Weight-100 size: the basis is 16
runs at scale 2^19 at excess 1.1 and the sub-region's relation-weighted means
(`FILT_TABLE`): worst N/U = 1/18.480 = 0.05411 (seed 11), charged **1.35 x
0.05411 = 0.07305**. The pad is the smallest value on the 0.05 grid whose
charge exceeds every weight-100 N/U measured at excess 1.1 at scale 2^14 and
above by at least 3% (the rule of caab2f05's pad 1.35 at excess 1.2):
`SCALE_NU11`, 57 runs at 2^15..2^20 (maxima 0.0706 (2^15, n = 30) / 0.0638
(2^16, n = 4) / 0.0587 (2^17, n = 2) / 0.0540 (2^18, n = 2) / 0.0541 (2^19, n =
16) / 0.0521 (2^20, n = 3); binding value 0.07056, 3.5% below the charge); the
means fall with scale. The merge ratio N(200)/N(100) at excess 1.1 is <= 0.7827
over 15 runs at 2^15..2^17 (`W200_11`; 42e80f80's excess-1.2 runs: <= 0.806),
charged 0.80. Charged N = 0.80 x 0.07305 x U < 2^59.3545 (11 occupied
2^56-element blocks; eleven-owner layout, STORAGE.md). **Occ>=3 purge losses.**
An ideal in 3 rows each with a distinct singleton (3 equations, 4 unknowns: the
f39a3f63 obstruction) is unrecoverable post-hoc; it is counted as unknown in
the missing-log calibration, whose covers bound the occurring remainder
including it (DESCENT.md; no recovery procedure and no domination term).
A = [M | S] with 21 Schirokauer columns (sigma-odd unit ranks 18 + 3). Block
Wiedemann, blocks (28, 5): m = 28 left vectors X with entries uniform in
[0, 2^128), n = 5 right vectors Y uniform in F_r. The five Krylov sequences
advance in lockstep, one matrix pass per step (STORAGE.md): per row and vector
201 x 17.9 sparse (4 + 3 + 0.88 x 39/5 + 0.12 x 33 = 17.82; H4u: >= 88% of
weight is +-1, measured 0.917-0.948 at weight 200 by 42e80f80) + 21 x 97 dense
+ row tail 1,750 (eleven owners); the Horner reconstruction uses single-vector
products at 201 x 20.7. Projections: (m + n) N^2 terms of a 128-bit X entry
times a field element at 60 (`bw_dot2`, counted 59), plus one 420-instruction
reduction per projection-sequence entry ((m + n) N in all).
LA **<= 2^118.061735269 units**, generator **<= 2^108.469324732 units**.

## 5-7. Descent and success ([7], [7c], [7e], [7l], [7p], [7r], [7s], [7t], [7u], DESCENT.md)

62-node plan and the sqrt(U)-weighted boot at UDUP = 0.25 (the relation-weighted
charge of 41025889, below every measured lattice), missing logs charged in one
ideal universe at the never-occurring mass in full plus the measured unknown
fraction of occurring ideals (covers c_r, DESCENT.md), prime-lattice missing-log
children and the inputs of section 3, with lazy cofactorisation and the lazy
boot (DESCENT.md; same failure events, only the own work changes): expected
two-log work **<= 2^115.1215 units**, meter 2^115.900, so P(cap) <= 0.582944
(Markov on the unrounded E[W]; H5 behind E[W], its three steps in DESCENT.md
"Recurrence and success"). The realized missing-log table: on an event E
(never-occurring mass within 2^-8 of its mean by Chebyshev; occurring-but-unknown
mass below its charge by delete-one Efron-Stein, which contains every
cross-ideal covariance of the purge) every conditional recurrence quantity is
at most its charged value; P(not E) <= 2^-33.92, charged 2^-30 (DESCENT.md "The
realized table"). Boot failure <= 0.002315866 per log (16 relations from
distinct lattices, independent given any table in E; shared lattice <= 2^-23.82);
node counter 2^28.5 with E[N] <= 2^18.21 (abort <= 2^-10.29; node
pairs <= N_MAX E[N]/2); special-q repeats <= 2^-18.59; missing-ideal repeats
under the shared missing set <= 2^-9.81; Block Wiedemann and the main-phase
factor each <= 2^(log2_N+2-128) with 128-bit X entries under H3, combined <=
2^-65.64 by the union bound; field map < 2^-48. Success >= 1 - 2 x
0.002315866 - 0.582944 - 2^-10.29 - 2^-18.59 - 2^-9.81 -
2^-23.82 - 2^-30 - 2^-48 - 2^-65.64 - 2/r >= 0.410508 (the meter is
set for a modelled success about 0.41: the expected two-log work may exceed
its model value by 0.050 bits before the bound reaches 0.39;
DESCENT.md "Headroom and sensitivity", [7s]).

## 8. Ledger (log2 units)

| row | bound |
|---|---:|
| relation collection | 120.291319028 |
| linear algebra / generator | 118.061735269 / 108.469324732 |
| two-log meter | 115.900 |
| canonical list / Schirokauer / misc / factor base / certificates (section 9) | 90 / 88 / 80 / 72 / 50 |
| field map, pairings, verification, output | 66 |

Sum of the declared rows 2^120.626057130 (unrounded rows 2^120.626057129) < **2^120.626057131**.

## 9. Side budgets (derived)

Each row is a fixed power of two in `model16.total_log2`; the derivations give
explicit bounds well below the rows (schoolbook counts where no kernel exists).

- **Canonical list, 2^90 units.** Per eligible a: LLL of a O_h (< 2^16
  instructions) plus Fincke-Pohst associate enumeration under a 100,000-node
  budget (< 2^25/a; budget-exhausted a dropped as non-canonical). 2^63.006
  eligible a ([3]): < 2^89 instructions = 2^74.7 units.
- **Schirokauer maps, 2^88 units.** Per relation (< 2^64.2) and side:
  (a - b x)^E mod (r^2, h, S), < 18,400 ops on 72-coefficient 510-bit
  polynomials (< 2^20 each), < 2^35.2/relation: < 2^99.4 instr = 2^85.2 units.
- **Filtering and miscellaneous, 2^80 units.** <= 64 entries/relation;
  singleton/clique removal and the weight-200 merge touch < 2^70.2 entries
  per pass at < 2^8 for < 2^6 passes; back-substitution one more pass; total
  < 2^85 instructions = 2^70.7 units.
- **Factor base, 2^72 units.** Primes p < 2^69 (sieve < 2^72 instr), roots of
  h mod p (< 2^13) and of f(s, x), g(s, x) per root (< 2^14): < 2^78.4
  instructions = 2^64.1 units.
- **Field map, pairings, verification, output, 2^66 units.** Roots, basis
  change, two pairings, [a]P0, check e(P, Q) = z: < 2^20 F_p12 mults; charged
  2^66 for input checks and output on every processor.
- **Certificates, 2^50 units.** `lattice_count_bounds.py` ([3]) and tower
  checks ([0]-[0c], [13]): < 2^40 instructions.

## 10. Heuristics

H1 smoothness (Dickman/Murphy) and H1' threshold retention, including random
ECM group orders in collection and the f/g independence given survival that
the early abort's q_f uses; H2 uniform pool and work allowance; H3 matrix
(class numbers, unit maps, one-dimensional kernel; N/U and the weight-200
ratio not rising at production); H4 measured kernel counts and H4u +-1 share;
H5 descent model (DESCENT.md: the Poisson occurrence model with the measured
unknown fraction of occurring ideals, the realized-table premises H5c-mean (the
production mean unknown fraction per lambda class is at most the calibration
maximum) and H5c-ES (the delete-one influence ratio XI stays below 2^14, measured
at 82-3,350 at excess 1.1), the negative-binomial relation counts, the
independence of distinct cells' factorisations within a node (H5a), the lazy
cofactorisation order (a fresh random permutation, independent of relation
outcomes by construction), the measured unit-duplicate charge, and the
worst-sextile norm bins).

## Provenance

Core: Gautam Botrel with Claude Code and Codex (3cb69d49, 8dc444f6, 42e80f80;
refinements 6282ab79 / 3a00ffc8 / e3940d11 / 04073dd1 / 8c1635f2); line sieve
to 2^48 with band strip: 8cc93834 / a31a5f1d; early abort: 90d99567
(mpjunior92; origin 2c2ce2ab, gbotrel). Retention and ECM modules, kernel
charges, UDUP evidence, sub-region method: f34f63d9 / 41025889 / 1a449d81;
filtering campaign design: caab2f05; prime-lattice missing-log children and
the calibration method: d0b0f071 (gbotrel); sqrt(U) boot weighting
(Guillevic, ASIACRYPT 2015) measured by a Claude Opus 5.5 crew agent. Excess
1.1 with its measurements, the sub-region, the row tail, the H5 argument and
side budgets: 273b1519 / 7c1acb81 / e178a7af / c0aa50a8 / 6895f2e2 (Claude
Fable 5.1). Validated base b0b3b74c (calibration-only missing-log charge,
uniform sampler, re-measured boot bins): muse-spark-1.3-contributor in Muse
Code via T3 Code. Realized-table control (delete-one Efron-Stein, the jk C
transcription, sub-band edge bounds, pair-count repeat bound): the esconc study
(Claude Opus 5.5, never submitted). This merge (early abort, the control and
the high-lambda tail measured at excess 1.1, plan re-search, node counter and
meter, pooled top boot bin): Claude Opus 5.5 (effort xhigh), T3 Code (delegate
task). Lazy cofactorisation and lazy boot, the five-vector Krylov kernels
(`krylov5.s`, counted 39/33) and the 128-bit-X projection kernel (`bw_dot2.s`,
counted 59) with blocks (28, 5): ported from the recQ v6 package 9332a74d
(mpjunior92, Claude Opus 5.5 in Claude Code) onto this recP line by
muse-spark-1.3-contributor (effort max) in Muse Code via T3 Code; the recP
reduction, collection, filtering, plan and success accounting are unchanged.

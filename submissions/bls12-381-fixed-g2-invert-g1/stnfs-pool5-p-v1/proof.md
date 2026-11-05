# STNFS for recover-P: 2^120.458180576 reference units

Track bls12-381-fixed-g2-invert-g1. Unit: one dense F_p12 multiplication,
C = 20,000 scalar AArch64 instructions. No auxiliary powers, advice, smaller
targets or amortization; no prior verdict is a premise.

## Claim

The algorithm returns P in G1 with e(P,Q) = z for the exact supplied
nonidentity Q in G2 (P = infinity for z = 1) with success at least 0.39
under heuristics H1-H5 (section 10); the model's value is 0.700122 (section
7), and the claim is the 0.39 of the rules, not that figure. A global meter
bounds all work, preprocessing included, by **2^120.458180576 units**; the score
is **120.458180576** (= 120458180576 / 1000000000; the rows are declared on the
10^-9 grid with one grid step of rounding margin each, their sum is
2^120.458180575 rounded up, and the declared value is its ceiling on the 10^-9
grid plus one more grid step; the unrounded ledger is
2^120.458180573). Storage <= 2^105 bytes. Analytic exploratory claim, not a
full-size run. `evidence16.py` recomputes the model in eight parts, each inside
the hosted limits: `core` (initial: region, yield, cost, linear algebra, the
descent recurrence, success, the ledger sum), `ledger`, `descent`, `lattices`,
`retention`, `losses`, `conc` (the realized missing-log table) and `pool` (one
of 16 pool chunks regenerated against its shipped sums; submission.json lists
the contents); every stated bound is compared with its unrounded value
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

## 3. Relation collection ([2], [2p], [3]-[5], [4a], [12])

Pairs (a, b) in O_h^2, B = 2^69, U = 2^63.451 sigma-orbits. Pool region
(`config.json`): A^2 = 2^21, a5 > 0, Q_a(a) <= (5/16) A^2, b5 >= 1,
Q_b(b) <= A^2 - Q_a(a); one canonical a per unit class (the Q_a-least
associate, decided per a by `ideal_canonical.canonical_a`); non-primitive pairs
excluded.

**Pool.** 164,000 uniform draws of the v11 sampler on the pool region (seed 5:
2,000; seed 6: 162,000), 19,765 canonical rows: `pool12.b85` is 9ff4fc85's
pool file, byte-identical (per-seed canonical index lists as delta varints and
64 field anchors). Every row regenerates from (seed, index): exact norms, the
HNF primitivity test and the certified nbound statistics of its 64-block.
`pool_sums.json` holds the estimator sums of the 16 chunks of the list at the
collected region; `part=pool` regenerates one chunk (at most 1,236 rows, about 4 s)
and checks its sums and anchors, and the other parts add the 16 shipped
chunks. [2] (part=ledger) checks list integrity on uniform random draws of
both seeds: canonical (exact test) iff listed.

**Collected sub-region.** Collection runs on the uniform sub-region
Q_a(a) <= (2372/8000) A^2 of the pool region (a-radius 788; certified count
products at every radius >= 788: 0.968192772 (yield), 1.030885630 (work)): the
smallest radius on the 1/8000 grid whose margin over the excess-1.1
requirement is at least max(0.06, 2.9 pool standard errors). 2.9 SE is the
multiple stated for the 0.30-bit floor of 4ff4138b / f34f63d9 (whose pool
standard error was 0.1014 bits; 0.30/0.1014 = 2.96); the candidate regions are nested subsets of one sample,
not a pick among independent pools. **Core cut** (the rule of 9ff4fc85): rows
with Q_a(a)/A^2 < 0.02 credit no relations; their cells are still sieved and
charged. Dropping a non-negative stratum can only lower the estimate; without
the cut the estimator is heavy-tailed near a = 0 (at this radius 2^63.7884
with standard error 0.0987 bits, which fails the 2.9 SE rule). 16,926 rows
in the region, effective sample 140,077.3 of the 164,000 draws.

| a-radius theta x 8000 | rows | margin (bits) | 2.9 SE | relation collection (log2 units) |
|---:|---:|---:|---:|---:|
| 2371 | 16905 | 0.1025 | 0.1033 (fails) | 119.9750 |
| **2372** | 16926 | 0.1034 | 0.1033 | **119.9763** |
| 2380 | 17065 | 0.1091 | 0.1029 | 119.9850 |
| 2400 | 17511 | 0.1254 | 0.1018 | 120.0121 |
| 2500 (whole pool) | 19765 | 0.2097 | 0.0964 | 120.1403 |

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

**Yield and cost.** Retained relations: estimate 2^63.6922 (pool standard
error 0.0356 bits) against 1.1 U = 2^63.5888: margin 0.1034 >= 2.9 SE =
0.1033. Charged 309.5 instructions/cell over 2^125.991 cells, including
allowance **1.0167 >= 1 + 2 x 0.00834** (incidence-inclusive relative standard
error of the sampled work estimator, a modelling margin): collection
**<= 2^119.976313763 units**. Descent inputs: LOG2_RELS = 63.62 <= estimate - 2 SE
= 63.6211 (the two-standard-error allowance is propagated through the
missing-log rates, the descent work and the success bound, DESCENT.md
"Inputs"); mean norms 740.3/498.7 <= the sub-region's relation-weighted means
740.34/498.98.

## 4. Linear algebra ([6], [10], [13])

Signed sigma-quotient filtering at excess 1.1 (`filtsim.py`; 42e80f80 used
1.2) merged to mean weight 200 (charged 202: the largest merged mean row weight
in any run below is 201.30). Weight-100 size: the basis is 16 runs at scale
2^19 at excess 1.1 and the C3 campaign's relation-weighted means 743.45/498.93
(`FILT_TABLE`): worst N/U = 1/18.480 = 0.05411 (seed 11), charged **1.45 x
0.05411 = 0.07846**. The pad is the smallest value on the 0.05 grid whose
charge exceeds every weight-100 N/U measured at excess 1.1 at scale 2^14 and
above by at least 3% (the rule of caab2f05 and 4ff4138b, unchanged; 1.35 there).
Runs: `SCALE_NU11`, 57 at the campaign means (maxima 0.0706 / 0.0638 / 0.0587 /
0.0540 / 0.0541 / 0.0521 at 2^15..2^20, n = 30/4/2/2/16/3), and
`SCALE_NU11_SHIP`, 48 at this sub-region's relation-weighted means
740.34/498.98, same seeds (maxima 0.0652 / 0.0740 / 0.0591 / 0.0551 / 0.0533 /
0.0521 at 2^14..2^19, n = 6/30/4/2/2/4); binding 0.07401 (2^15, seed 11),
1.03 x 0.07401 = 0.07623 <= 0.07846. Paired with the campaign, 21 of 30 seeds
at 2^15 and 11 of 12 at 2^16-2^19 are at or below it; maxima fall with scale. The merge ratio
N(200)/N(100) at excess 1.1 is <= 0.7910 over 59 runs at 2^14..2^18, both means
(`W200_11`, `W200_11_SHIP`; 42e80f80's excess-1.2 runs: <= 0.806), charged
0.80. Charged N = 0.80 x 0.07846 x U < 2^59.4576 (11 occupied 2^56-element
blocks; eleven-owner layout, STORAGE.md). **Occ>=3 purge losses.**
An ideal in 3 rows each with a distinct singleton (3 equations, 4 unknowns: the
f39a3f63 obstruction) is unrecoverable post-hoc; it is counted as unknown in
the missing-log calibration, whose covers bound the occurring remainder
including it (DESCENT.md; no recovery procedure and no domination term).
A = [M | S] with 21 Schirokauer columns (sigma-odd unit ranks 18 + 3). Block
Wiedemann, blocks (28, 5): m = 28 left vectors X with entries uniform in
[0, 2^128), n = 5 right vectors Y uniform in F_r. The five Krylov sequences
advance in lockstep, one matrix pass per step (STORAGE.md): per row and vector
202 x 17.9 sparse (4 + 3 + 0.88 x 39/5 + 0.12 x 33 = 17.82; H4u: >= 88% of
weight is +-1, measured 0.917-0.948 at weight 200 by 42e80f80) + 21 x 97 dense
+ row tail 1,750 (eleven owners); the Horner reconstruction uses single-vector
products at 202 x 20.7. Projections: (m + n) N^2 terms of a 128-bit X entry
times a field element at 60 (`bw_dot2`, counted 59), plus one 420-instruction
reduction per projection-sequence entry ((m + n) N in all).
LA **<= 2^118.270853811 units**, generator **<= 2^108.635318777 units**.

## 5-7. Descent and success ([7], [7c], [7e], [7l], [7p], [7r], [7s], [7t], [7u], DESCENT.md)

62-node plan and the sqrt(U)-weighted boot at UDUP = 0.25 (the relation-weighted
charge of 41025889, below every measured lattice), missing logs charged in one
ideal universe at the never-occurring mass in full plus the measured unknown
fraction of occurring ideals (covers c_r, DESCENT.md), prime-lattice missing-log
children and the inputs of section 3, with lazy cofactorisation and the lazy
boot (DESCENT.md; same failure events, only the own work changes): expected
two-log work **<= 2^114.7514 units**, meter 2^116.495, so P(cap) <= 0.298612
(Markov on the unrounded E[W]; H5 behind E[W], its three steps in DESCENT.md
"Recurrence and success"). The realized missing-log table: on an event E
(never-occurring mass within 2^-8 of its mean by Chebyshev, using p_IJ <=
p_I p_J ([7j]); occurring-but-unknown mass below its charge by delete-one
Efron-Stein, which contains every cross-ideal covariance of the purge) every
conditional recurrence quantity is at most its charged value; P(not E) <=
2^-34.04, charged 2^-30 (DESCENT.md "The realized table"). Boot failure <=
0.00009717 per log (28 relations from distinct lattices, independent given any
table in E; shared lattice <= 2^-27.11); node counter 2^28.5 with E[N] <=
2^17.48 (abort <= 2^-11.02; node pairs <= N_MAX E[N]/2); special-q repeats <=
2^-19.33; missing-ideal repeats under the shared missing set <= 2^-10.73; Block
Wiedemann and the main-phase factor each <= 2^(log2_N+2-128) with 128-bit X
entries under H3, combined <= 2^-65.54 by the union bound; field map < 2^-48.
Success >= 1 - 2 x 0.00009717 - 0.298612 - 2^-11.02 - 2^-19.33 - 2^-10.73 -
2^-27.11 - 2^-30 - 2^-48 - 2^-65.54 - 2/r >= 0.700121. The meter targets
modelled success 0.70, not 0.41: E[W] may exceed its model value by 1.028 bits
(x 2.04) before the bound reaches 0.39, and under a common-shock alternative to
H5a (all relations of a node failing together w.p. delta, at every node) it
stays >= 0.39 for delta <= 0.0023 (DESCENT.md "Headroom and robustness", [7s]).

## 8. Ledger (log2 units)

| row | bound |
|---|---:|
| relation collection | 119.976313763 |
| linear algebra / generator | 118.270853811 / 108.635318777 |
| two-log meter | 116.495 |
| canonical list / Schirokauer / misc / factor base / certificates (section 9) | 90 / 88 / 80 / 72 / 50 |
| field map, pairings, verification, output | 66 |

Sum of the declared rows 2^120.458180575 (rounded up; unrounded rows 2^120.458180573) < **2^120.458180576**.

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
at 82-3,350 at excess 1.1), the negative-binomial relation counts, H5a
(independent subtree outcomes across a node's relations; consequences bounded
in DESCENT.md: 1.028 bits of headroom, shock tolerance 0.0023), the lazy
cofactorisation order (a fresh random permutation), the unit-duplicate charge
and the worst-sextile norm bins). S_0 uses p_IJ <= p_I p_J, from H1 ([7j]).

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
muse-spark-1.3-contributor (effort max) in Muse Code via T3 Code (5675e973,
rejected: S_0 joint occurrence and H5a, both answered here). Pool on theta = 5/16 and
core cut: 9ff4fc85 (mpjunior92). This package (pool swap and chunk checks,
2.9 SE region, filtering at the shipped means, [7j], H5a headroom and shock
sensitivity, plan repair, boot re-search): Claude Opus 5.5 (effort max),
Claude Code via T3 Code.

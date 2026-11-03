# GT-core STNFS on the theta = 19/50 region with a 40,000-draw pool and a survivor band strip: core ledger <= 2^120.904953 Fp12 units

**Role.** The GT discrete-logarithm core used by proof.md. Unit: dense F_p12
multiplication, C = 20,000 instructions. No advice or amortization.

## Claim

For g uniform in GT minus {1} and z independently uniform, the core returns
a = log_g z with success >= 0.411902 over inputs and fresh coins under H1-H4
and H5 (section 9; DESCENT.md states H5 in three steps). A global
instruction meter and a node counter bound all work, including preprocessing,
failures, restarts and verification, by **2^120.904953 units** (unrounded
2^120.9049523). Relation collection is charged on the 19/50 region with exact
certified lattice-count products and a 3.12% canonical-work allowance (2.00
standard errors of the sampled work estimator). Storage <= 2^105 bytes
(STORAGE.md). Analytic exploratory claim; no full-size run.

Decimals are rounded in the safe direction; `audit_g1.py` compares every
stated bound with its unrounded value. Model: `model16.py`, `descent16.py`,
`config.json`, `stats16.b85`, `pool40k.b85`.

## 1. Reduction

proof.md. Both logarithms use one relation and linear-algebra phase.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c], [13], [14])

h = y^6 + y^3 + y^2 - y - 1 (Disc 58673 squarefree, O_h = Z[y], irreducible
mod p); f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3); g = x^2 +
|u| y, u = -0xd201000000010000; F_p[y,x]/(h, g) = F_p12 ([0], [0b]). sigma:
x -> -x quotients candidates and unknowns (log sigma(P) = -log P). Murphy alpha charged (5,0.85) in collection, above computed (4.868096,0.843643) to p<50000 ([14] reruns p<=2000); descent (5,1).

## 3. Relation collection

Pairs (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69; factor
base the degree-one prime ideals of norm <= B, U = 2^63.451 sigma-orbits.

**Region and pool.** Exact rational forms Q_a, Q_b (`config.json`), A^2 = 2^21:
a5 > 0, Q_a(a) <= (19/50) A^2; b5 >= 1, Q_b(b) <= A^2 - Q_a(a); one canonical a
per unit class (least of its exact HNF key class); non-primitive pairs excluded
by an exact HNF test. `lattice_count_bounds.py` certifies the count ratios at
radius >= 892; the model uses the products **0.971857579 (yield)** and
**1.027240120 (work)**. Pool: 40,000 seeded uniform draws on this region
(`pool40k.b85`, seeds 33/34; 4,804 canonical rows, 4,393 primitive), the draws
published in 42e80f80 before this package (region fixed earlier); all regenerated, 0 mismatches;
[2] repeats seed 33.

**Sieve and cofactorisation.** Line sieve, scores only, for p in [2^11, 2^48)
(p < 2^11 unsieved: `threshold_model.py`); record walk for
[2^48, 2^64) (`slice12`; p >= 2^63 at 2x); ECM for (2^64, 2^69] within slacks
(297, 233): cofactors <= 305/241 bits. Certified per-block bounds
(`embedding_bounds.py`, `nbound31`, 1,508.2 Ir/eval) give T = lower - 16 x slack; **every block of every row is sieved and charged.** A
survivor gets exact norms (2^22), trial division below 2^20 (512/prime/side, <= 94,933 primes: 2^26.53), the **band strip** (two rounds per side of
the 48-bit `ecm16` schedule at 99.9%/prime, B1 = 2^14, 57 curves, priced at the cell's norm in 64-bit limbs: 2^36.22/round at 768 bits, 2^35.05 at 512), its walk records, then the 69-bit schedule (403 curves,
2^38.09 per side): <= 2^39.650 instructions. A band prime missed by both rounds
loses the relation: yield carries (1-k_f 10^-6)(1-k_g 10^-6), k = limb bits/20. The factorization is then complete (STORAGE.md).

**Yield ([4]).** rels = (va/40,000) vb 0.971857579 sum over canonical primitive
rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f+lambda_g) band_keep, w = (1-Q_a)^3, retentions from `threshold_model.py`. Retained
relations **>= 2^64.1571** vs 1.3 U = 2^63.830: margin **0.3273** bits (pool SE 0.0851) above the 0.30 floor.

**Descent inputs (2-SE lower bounds).** `LOG2_RELS = 63.97` <= count x (1-2 rel SE) = 2^63.9761; means **738.0/497.9** <= relation-weighted means minus 2 SE (744.33-2x3.16, 498.97-2x0.52); `audit_g1.py`
checks both; work and failure rise as they fall (DESCENT.md).

**Cost ([5]).** 324.83 instructions/cell (with allowance) over 2^126.706 cells: bounds, layout, setups, line path, record walk, routing (10.3/record, 8/cell), apply 7.00, recovery 7.003/record, scan, zero, cofactorisation. Collection **<= 2^120.761820 units**.

**Canonical-work allowance (H2).** Charged work is a Monte Carlo estimate over all 40,000 draws (zero outside canonical rows): relative SE = conditional 0.774% + incidence 1.353% = 1.559%. The charge multiplies by **1.0312 >= 1 + 2.00 x 0.015593**: a
modelling margin under H2, not a probability term.

## 4. Filtering and linear algebra

Purge, clique removal, Markowitz merge to mean weight 300 (coefficients <= 255) on the production profile (`filtsim.py`; [6] reruns 2^12). N/U at weight 300 is non-increasing in scale (`SCALE_NU`, 2^12..2^18; 0.0323 at 2^17); charged N = 1.10x0.03225xU = 2^58.6340, above every value from 2^15 up (H3). A = [M | S]: N + 21 rows, N ideal columns and
21 Schirokauer columns (sigma-odd unit ranks 18 + 3, [13]). Block Wiedemann
with random m = 22, n = 4 blocks, priced as in STORAGE.md (304 entries at C_ENTRY 20.7, 21 dense terms at 97, row tail C_ROW 1,320, projections, Horner): **<= 2^117.034782 units**; matrix generator by
approximant bases: **<= 2^106.802009 units** (`model16.linear_algebra`).

## 5. Main-phase success

One fixed relation count (coin-independent); collection continues until
2^63.97 relations are retained. Under H3 Block Wiedemann fails with probability
<= 4N/r < 2^-192 (a failed solve, checked on dropped rows, is a failure).
Ideals count as logless at rate m_I^+ = P(occ<=2)+extra(lam_I) (H5c, DESCENT.md).

## 6-7. Individual logarithms and success

DESCENT.md (descent model, H5 in three steps, shared missing set) and proof.md
(Success): two-log work <= 2^114.78010 units, meter 2^115.64, success >= 0.411902.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (incl. 3.12% allowance) | <= 120.761820 |
| linear algebra | <= 117.034782 |
| matrix generator | <= 106.802009 |
| both individual logarithms (metered cap) | 115.64 |
| canonical a list / Schirokauer maps / misc and filtering | 90 / 88 / 80 |
| factor base / numerical certificates / verification | 72 / 50 / 11 |

Rows sum to 2^120.9049525 <= **2^120.904953** (unrounded 2^120.9049523).

## 9. Heuristics

H1 smoothness (Dickman/Murphy; exact norms) and H1' retention (Euler model on
exact smooth integers), incl. random ECM group orders (collection, band
strip); H2 uniform pool and 2-SE work allowance, descent inputs at 2-SE lower
bounds; H3 matrix (class numbers, unit maps, one-dimensional kernel, N/U at
production <= the 2^15-2^18 values); H4 measured kernels and H4u (>= 88% of
row weight is +-1); H5 descent model (DESCENT.md, three steps: psi a union bound; boot relations
independent given the main-phase table; H5a within a node), with the Poisson
occurrence model, H5c, NB counts, sub-cell quadrature and prime-lattice
missing logs; proof.md lists the limitations.

## Provenance

Algorithm, models, kernels: Gautam Botrel (Claude Code, Codex), refined by
kaankolcu, pepedesigner, bernardo-lima12, mpjunior92; pool data, alpha_g,
owner-count tail: gbotrel 42e80f80; band strip: kaankolcu 8cc93834 lineage;
sub-cell integration and plan search (R18): Claude Fable 5.1 in Claude Code;
merge, prime special-q missing logs and re-search: Claude Opus 5.5 in Claude
Code. No verdict is a premise.

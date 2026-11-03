# GT-core STNFS over a special tower on a theta = 189/500 sub-region with a 99.9% main-phase ECM schedule: core ledger <= 2^121.0820 Fp12 units

**Role.** The GT discrete-logarithm core used by proof.md: the recover-Q argument of e3940d11 with REVISION.md applied in place. Section 1 (the source's G2 output step) is not used; proof.md has the G1 reduction. Unit: dense F_p12 multiplication, C = 20,000 instructions. No advice or amortization.

## Claim

The core below returns a = log_g z for g uniform in GT minus {1} and z independently uniform, with success at least 0.39 over inputs and fresh coins under the declared heuristics H1-H5;
the modelled success is at least 0.390320 (> 0.39; rounded below the model value 0.3903323). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, verification and output, by **2^121.0820 units** (core bound). The ledger (section 8) totals at most 2^121.0820 (unrounded
2^121.0819193). Relation collection is
charged on the theta = 189/500 sub-region of the sampled region with exact
certified lattice-count products and an explicit 7.06% canonical-work allowance
(2.00 incidence-inclusive standard errors of the sampled work estimator, [4a]). Distributed storage is at most
2^105 bytes (STORAGE.md) in individually bounded, explicitly addressed objects.
The claim is analytic and exploratory under the NFS heuristics H1-H5 (section
9) with the measurements and checks cited; no full-size run is made.

**Numerical conventions.** Decimals are rounded values; every bound is rounded
in the safe direction, and `declared_bounds` compares each stated bound with the
unrounded value.

The model is `model16.py` and `descent16.py` with `config.json` and `stats16.b85`;
`evidence16.py` recomputes the ledger and the checks cited as [n].

## Base

The machinery (polynomials, factor base, sieve, filtering, block Wiedemann,
descent, kernels, pools) is the lean STNFS package of the sister track
`bls12-381-fixed-g2-invert-g1` (cbab3319). Changes here: a theta = 189/500 uniform
sub-region (was 87/200) with exact certified count products and a work
allowance (section 3), alpha (5, 0.844) checked in-package ([14]), an explicit
special-q node counter (section 7) and the descent plan of section 8.

## 1. Reduction

Not used (proof.md has the G1 reduction). The core computes a = log z / log g
mod r for g uniform nonidentity and z independently uniform, both logs against
one shared relation and linear-algebra phase.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks that
(u y)^((p^6-1)/2) = -1 in F_p6, so g is irreducible over F_p6 and
F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. The map
x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P), so
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at
(5, 0.844) in relation collection, above the in-package values (4.868096,
0.545620) over p <= 2000 ([14], `evidence16.py --alpha`, 8c1635f2 lineage);
the descent keeps (5, 1).

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69;
the factor base is the relative-degree-1 prime ideals of norm <= B on both
sides, U = 2^63.451 sigma-orbits (64-bit indices).

**Region.** Exact rational forms Q_a, Q_b (`config.json`), A^2 = 2^21. The pool
is uniform at theta = 1/2 (a5 > 0, Q_a(a) <= A^2/2); the collected region is
its uniform sub-region theta = 189/500 (`production_theta_region`): a5 > 0,
Q_a(a) <= (189/500) A^2, half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The
experiment keeps the 1019 canonical sub-region rows (945 primitive) of the 2401
and scales by the sub-region volume over the effective sample 20000 x
(189/250)^3 = 8641.62432. The canonical a of each unit class is the least of its
exact HNF key class in (Q_a, lex) order, so each orbit of <units, sigma> has one
point in the domain; non-primitive pairs are excluded by an exact HNF test.
`lattice_count_bounds.py` (8c1635f2 lineage) certifies from exact rational pi
and determinant intervals the count ratio at every radius R >= 890 =
isqrt((189/500) 2^21) (Q_a) and R >= 1142 = isqrt((311/500) 2^21) (Q_b
slices): [0.992213, 1.007059], [0.993464, 1.005931], with unit-cube Q-shift
t = least k/1000 above sqrt(sum|Q_ij|)/2: 261/250, 563/500 (OR1). The model uses their exact products **0.985729314
(yield)** and **1.013031277 (work)**; pool sampling error is handled by the
allowance below. The canonical list costs 2^90 units. **Pool:** 20,000 seeded
uniform pairs (`pool12-*.b85`), 2,401 canonical rows with exact norms, flags and
bound statistics; [2] regenerates and verifies them.

**Sieve and cofactorisation.** Line sieve for p in [2^11, 2^48) (scores
only, A2R1; p < 2^11 is the unsieved small group of `threshold_model.py`,
counted against the slacks), record walk for p in [2^48, 2^64) (`slice12`,
KERNELS.md; p >= 2^63 charged 2x, R6), primes in (2^64, 2^69] by ECM
within slacks (297, 233): survivors' cofactors are at most 305 (f) and
241 (g) bits. Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected to within 8 bits, rounded outward;
`nbound31` replays the leaves bit for bit at 1,508.2 Ir per evaluation, charged
exactly; [2] checks against exact norms) give thresholds T = lower - 16 x slack.
Survivors get exact norms, trial division below 2^20, then the f side first:
the two-round 48-bit f band strip (A2R1, STORAGE.md) and f `ecm16.py` at
**99.9% per prime** (R1; [12]): B1 = 2^16, B2 = 30 B1, wheel 30030, 403
curves, 367.76 Ir per multiplication (`ecm12`). The g-side band + ECM runs
only if f factors fully (LEAP1 early abort, from validated recQ 2c2ce2ab),
charged at q_f = min(1, 100 rho_f / P(omitted_f <= slack_f)) <= 0.0017:
2^38.0907 + q_f 2^38.0881, at most 2^38.81 per survivor all-in.
**No block selection:** every block of every row is enumerated, sieved and
charged. Yield is untouched (a failed f side yields no relation anyway).

**Yield ([4]).** rels = (va_region/8641.62432) vb 0.985729314 sum over canonical
primitive sub-region rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f +
lambda_g) band_keep: Dickman rho ([1]), retentions from `threshold_model.py`
(calibrated on exact smooth integers, B = 2^10..2^14), the ECM miss and the
band-strip keep factor (1 - k_f 10^-6)(1 - k_g 10^-6) per cell (A2R1).
Retained relations >= 2^64.1302 against 1.3 U = 2^63.830: margin >= 0.3004
bits (pool SE about 0.152), above the 0.30-bit floor (`MARGIN`). The descent
input `LOG2_RELS = 63.97` is below the count times (1 - 1 relative pool SE),
2^63.9700 (R8; audit `yield_lower_bound_covers_descent`) and its mean norms
745.9/498.9 are rounded down from 745.95/499.03; [4] checks both (`rels_input`).

**Cost ([5]).** Charged 333.4 instructions per cell over
2^126.654 cells: bounds, layout and setups, line path, record walk, routing (10.3
per record, R2 + R4; 8 per cell), apply 7.00, recovery 7.003 per record (every record re-read), scan, zero and
cofactorisation (with f-first early abort). Relation collection **<= 2^120.7469 units** (sample estimate
2^120.6486 x allowance 1.07057).

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo
estimate of a fixed total over all 20,000 draws (zero outside the sub-region or
non-canonical): relative SE = conditional 1.700% and incidence
sqrt((1 - 1019/20000)/1019) = 3.052%, total 3.493%. The charge multiplies by
**1.07057 >= 1 + 2.00 x 0.034933** ([4a] fails below 1 + 2.00 SE): a declared
modelling margin under H2, not a probability term.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge with coefficients bounded by 255,
simulated on the production occurrence profile (`filtsim.py`; scale 2^17 table
embedded, [6] reruns scale 2^12): the simulated matrix is the signed sigma-quotient system (each occurrence takes
the sign of its conjugate orientation, log sigma(P) = -log P). Three seeds at
2^17 give N/U <= 0.03225 at mean weight <= 304 (charged 304 per row; R17). N/U falls with
scale: the maximum over seeds of `filtsim.run` (kmax 64, mean weight 300) is
0.0652 (2^12, 10 seeds), 0.0431 (2^13, 6), 0.0368 (2^14, 6), 0.0348 (2^15, 3),
0.0334 (2^16, 3), 0.0323 (2^17, 3), 0.0304 (2^18, 2) (`SCALE_NU`, checked
non-increasing by [10]; `audit_g1.py` {"unit_entries": [k, seed]} measures H4u at any scale up to 2^15).
The model charges N = 1.10 x 0.03225 x U = 2^58.6340: the charged ratio 0.03547 is
above every measured value from 2^15 up (max 0.03480; the base charged 1.61 x 0.0471 to cover
2^12 at weight 100). That N/U at production scale stays at or below the 2^15-2^18 values is
part of H3, an extrapolation of a monotone measured trend. The system solved is the square matrix
A = [M | S]: N + 21 filtered rows (the rest verify the result), the N
ideal columns M and 21 dense Schirokauer columns S. As columns are sigma-quotiented
(log sigma(P) = -log P), only sigma-odd units need maps: signatures K_f (0, 36),
K_f^+ = K_h(x^2) (0, 18), K_g (2, 5), K_h (2, 2) give odd unit ranks 35 - 17 = 18
and 6 - 3 = 3 ([13], Sturm sequences); under H3 the kernel of A is the line of virtual
logarithms (with their Schirokauer coefficients). Block Wiedemann
on A uses random blocks: m = 22 left and n = 4 right vectors uniform in F_r (the
coins; A2R2). The Krylov phase takes N/m + N/n block steps, i.e. N(n/m + 1) products by
A. Per row a product costs 304 sparse entries at C_ENTRY = 20.7 (H4u), 21 dense terms
at 97 (bw_full, into the row's 9-limb accumulator) and the row tail C_ROW = 1,320,
which reduces the sparse pair and the dense accumulator mod r (STORAGE.md, every
kernel counted by [10]). Every step forms the m x n projections x^T A^i y, (m + n) N^2
dense terms at 97 plus one reduction (C_RED9 = 420) per projection
entry. The kernel vector is reconstructed by Horner, w <- A w + Y F_i over the
deg F <= N/m + N/n generator coefficients (Y the n starting vectors): one product
plus, per element, n dense terms and one reduction. With c = 304 x 20.7 + 1,320 +
21 x 97 = 9,649.8 per row, the total is N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n + 420) +
97 (m + n)) + 420 (m + n) N instructions: **<= 2^117.0348 units** at
N = 2^58.6340 (minimum over integer block sizes at (m, n) = (22, 4)). The
(m + n) x (m + n) matrix generator of the length L = (N + 21)(1/m + 1/n) + 2
sequence is computed by divide-and-conquer approximant bases with Karatsuba
products: at most 4 (m + n)^3 L^log2(3) log2(L) field multiplications, each 97
plus a reduction 420: **<= 2^106.8021 units**.

## 5. Main-phase success

- **Relations.** One fixed count (region and polynomials are coin- and
  challenge-independent); H1 puts it at >= 2^64.1302, 0.3004 bits above 1.3 U;
  the descent uses 63.97 <= the count at 1 relative SE below (R8), and collection
  continues over further slices until 2^63.97 relations are retained (guard).
- **Linear algebra.** H3: r unramified in K_f and K_g ([13]), dividing neither
  class number, the 21 sigma-odd Schirokauer maps injective on units mod r-th
  powers, kernel of A the line of virtual logarithms. Block Wiedemann then fails
  with probability <= 4N/r < 2^-192 (Schwartz-Zippel, degree <= 2N); a failed
  solve (checked on the dropped rows) counts as failure.
- **Reconstruction.** Purged ideals' logs by back-substitution; the descent
  treats every ideal with at most two occurrences as having no log (section 6),
  at least 12.5 times the count [6b] measures in every band with a loss.

Main phase success >= 1 - 2^-192 over the coins (H1, H3).

## 6. Individual logarithms

DESCENT.md: an ECM boot on balls of certified p^6 lattices, then special-q
nodes sieving exact lattice balls (`dsieve16`) and recursing on child ideals.
Base-primitive cells 0.9100; retention >= 0.8575 at every node ([7],
`retain16.py`); cofactors by `ecm16.py` ([12]); negative binomial relation
counts. The 62-node plan and boot give expected two-log work **<= 2^117.5873
units** ([7]; one log <= 2^130.875 instructions). The meter caps both logs at
2^118.3128 units (a capped run fails) and a node counter (section 7) bounds the
node count by construction.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 1.780e-03 - 0.60478 - 2^-10.42 - 2^-10.70 - 2/r)
>= **0.390320** >= 0.39: main phase; two boots (phi <= 1.780e-03 each); the
instruction meter (cap 2^118.3128; Markov P <= 2^(117.5872764 - 118.3128) <=
0.60478); the special-q node counter (cap N_MAX = 2^27.3 nodes, incremented
before every ordinary or missing-log node, abort = failure; E[N] <= 2^16.872 by
`counts16.py`, unrounded 2^16.8719, so P <= 2^-10.42); special-q repeats (at most
N_MAX nodes in the smallest bin of 2^70.30 ideals, factor 2^6 for the 1/p
weighting: <= 2^(2 x 27.3 - 1 + 6 - 70.30) <= 2^-10.70, charged in full); and
2/r for the distance of power-randomized targets from the uniform full-group
model. DESCENT.md (Repeated ideals) has the meters. No lower bound on any
node's work is used; the boot term treats the 16 relations' subtrees as
independent (H5). The final pairing check rejects any incorrect output. The
bound is conditional on H1-H5; the H2 allowance is a resource margin, not a
probability term.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (sample estimate 120.6486 + 7.06% canonical-work allowance) | <= 120.7469 |
| linear algebra | <= 117.0348 |
| matrix generator | <= 106.8021 |
| both individual logarithms (metered cap) | 118.3128 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

W <= **2^121.0820** (unrounded 2^121.0819193), the core bound used by proof.md. The table's upper-rounded rows, including a 2^50-unit numerical-certificate allowance, sum to 2^121.0819329, below the stated core bound 2^121.0820.

**Configuration.** As sections 3-7; plan min mu 32.3, boot psi 0.674; boot + twelve plan entries re-optimised (R16).
Preprocessing is inside the 2^121.0820-unit cap; no advice, nothing amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; collection alpha (5, 0.844) vs computed
(4.868096, 0.545620) over p <= 2000 ([14]); descent (5, 1); exact norms); H1'
retention (Euler model on exact smooth integers); H2 region (uniform pool at
theta = 1/2, sub-region 189/500, exact flags, certified counts, work mean with
7.06% = 2.00 incidence-inclusive SEs (modeling margin, not probability); descent
input at/below yield less 1 relative SE, R8); H3 matrix (Schirokauer: r
unramified [13]; r-designated class/unit/kernel premises assumed, section 5;
filtering on the production profile, pad 1.10 over the 2^17 max at weight 300, N/U
non-increasing 2^12-2^18, assumed so at scale); H4 kernels (KERNELS.md,
STORAGE.md), H4u unit share >= 0.88; H5 descent (lower sums, NB recurrence,
missing logs, independent subtrees incl. logless re-descents, local omitted
mass, random ECM orders).

## 10. Memory

Storage <= 2^105 bytes (STORAGE.md). Polynomial selection and alpha follow
Guillevic-Singh, filtering CADO-NFS.

## Provenance

Algorithm, models, kernels, package: Gautam Botrel (Claude Code, Codex);
sister-track refinements by kaankolcu, pepedesigner (cbab3319 lineage).
Sub-region, work-SE formula, certified counts, in-package alpha: kaankolcu /
DeepSeek lineage (8c1635f2+). Port, node counter, plan, directed rounding,
this configuration: Claude Opus 5.5 in Claude Code (991a7c0d, 1eb4bd1f,
23f293a8, 377ee6d1, e554a0f4). Reuse credited; no verdict used as premise.

## Directed quadrature refinement

Inherited; DESCENT.md (Refined integration). Changed numbers: REVISION.md.

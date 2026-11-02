# STNFS over a special tower on a theta = 759/2000 sub-region, 99.9% ECM, X2 = 2^64, a measured record route, fully reduced LA kernels: 2^121.3270 reference units

Track bls12-381-fixed-g1-invert-g2. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used. Predecessor work is reused as credited below; no prior verdict is a premise.

## Claim

The algorithm below returns Q in G2 with e(P,Q) = z for the exact supplied nonidentity
P in G1 (Q = infinity for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins under the declared heuristics H1-H5;
the modelled success is at least 0.391369 (> 0.39; rounded below the model value 0.3914013). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^121.327036929110 units**; the score is
**121.327036929110**. The ledger (section 8): upper-rounded rows sum to at most 2^121.327036929110 (unrounded
model total 2^121.3269833). Relation collection is
charged on the theta = 759/2000 sub-region of the sampled region with exact
certified lattice-count products and an explicit 6.74% canonical-work allowance
(2.00 incidence-inclusive standard errors of the sampled work estimator, [4a]). Distributed storage is at most
2^105 bytes (STORAGE.md) in individually bounded, explicitly addressed objects.
The claim is analytic and exploratory under the NFS heuristics H1-H5 (section
9) with the measurements and checks cited; no full-size run is made.

**Numerical conventions.** Decimals are the experiment's values rounded as
shown. Figures used as bounds carry <= or >= and are rounded in the safe
direction; the `declared_bounds` check compares every stated bound (ledger,
collection, LA, relation count/margin, work SE, two-log, Markov, E[N],
counter-abort, repeat, boot failure, success) with the unrounded value and
fails if any is violated.

The model is `model16.py` and `descent16.py` with `config.json` and `stats16.b85`;
`evidence16.py` recomputes the ledger and the checks cited as [n].

## Base

The GT discrete-log machinery (polynomials, factor base, sieve, filtering, block Wiedemann, descent, kernels, pools) is the lean STNFS package of sister track `bls12-381-fixed-g2-invert-g1` (cbab3319), ported and refined up to e3940d11; both orientations need the same two GT logs. Changes: last section.

## 1. Reduction

Let g = e(P, Q0) (pinned `pairing(Q0, P)`), one pairing; Q0 is the fixed G2
generator. P has prime order r and e is nondegenerate, so g generates mu_r,
z = g^a for unique a mod r (a = 0 iff z = 1), and Q = [a]Q0 (one G2 scalar
mult), so e(P, Q) = z for that exact P (never replaced). log z and log g share
one relation/LA phase; a = log z / log g mod r.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks that
(u y)^((p^6-1)/2) = -1 in F_p6, so g is irreducible over F_p6 and
F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. The map
x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P), so
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at
(5, 0.844) in collection: [14] recomputes it in-package from the exact
Kummer-Dedekind structure over p <= 2000 (`--alpha` prints (4.868096,
0.843724)) and over p < 50000 (`--alpha50k` prints (5.085035, 0.731625)).
The p = 2 double root of g uses exact primitive t-side tails
2q(1-q^-7)/(q^2-1) (primitivity q/(q+1) included; valuations 15+ dropped,
conservative); the uncorrected mu/(q-1) overstated v by missing the
primitivity factor (prior-review finding, fixed here). alpha_f = 5 covers
the 2k partial 4.868096 with the tail analysis of REFINEMENT.md; alpha_g =
0.844 covers the 50k partial 0.731625 with a 15% margin (partials fall:
0.8437/0.7637/0.7540/0.7316 at 2k/10k/20k/50k). The descent keeps (5, 1).

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B on both
sides; U = 2^63.451 sigma-orbits (indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21. The pool is sampled uniformly at theta = 1/2
(eligible a: a5 > 0, Q_a(a) <= A^2/2); the collected region is its uniform
sub-region theta = 759/2000 (`production_theta_region` in `config.json`): eligible
a: a5 > 0, Q_a(a) <= (759/2000) A^2; half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a).
Because the sample is uniform on the theta = 1/2 ellipsoid, its restriction to
qa = Q_a/A^2 <= 759/2000 is uniform on the sub-region; the experiment keeps the
1034 canonical sub-region rows (960 primitive) of the 2401 canonical rows and
scales by the sub-region volume over the effective sample 20000 x (759/1000)^3 =
8744.91, which is the same estimator as (va(1/2)/20000) x (sum over sub-region
rows). The canonical a of each unit class
is the least of its exact HNF key class in (Q_a, lex) order, so each relation
orbit of <units, sigma> has exactly one point in the domain; non-primitive pairs
are excluded by an exact HNF test. `lattice_count_bounds.py` (the version of
this track's 8c1635f2 lineage, with a `minimum_radius` argument) certifies from
exact rational Machin-pi and determinant/sqrt intervals the count ratio at every
radius R >= 892 = isqrt((759/2000) 2^21), the smallest a-radius of the sub-region
(b-radii are at least sqrt((1241/2000) 2^21) = 1140): [0.985851, 1.013529] for Q_a
and [0.985805, 1.013529] for Q_b, each inside [0.98, 1.02] (asserted). The model
uses the exact products of these rational bounds, **0.971857578 (yield, 9-decimal floor)** and
**1.027240121 (work, 9-decimal ceiling)**, not rounded-out factors. Sampling error of the pool is
handled separately and explicitly by the canonical-work allowance below. The
canonical list costs 2^90 units.

**Pool.** 20,000 uniform eligible pairs from seeded indices (`pool12-*.b85`),
2,401 canonical rows with exact norms, primitivity flags and certified
norm-bound statistics; [2] regenerates sampled indices at theta = 1/2 and
verifies every field. The sub-region uses 1034 of the 2401 rows (960 primitive).

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^20), record
walk for p in [2^20, 2^64) (`slice12`, KERNELS.md; X2 = 2^64, was 2^63), primes in [2^64, 2^69] taken by ECM
within the slacks (294, 224). The f-side slack is 294, so survivors'
f-side cofactors are bounded by 302 bits; the g-side slack is 224, so survivors'
g-side cofactors are bounded by 232 bits (the slacks are re-optimised for the
sub-region, section 8). The per-survivor ECM charge is the
package's explicit `ecm16.py` schedule at **99.9%** per prime (ECM_FIND = 0.999, raised from the
base's 99%; this is also the `ecm16.schedule` default the descent uses), which [12] checks: the cheapest
grid point for 69-bit primes in a 302-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 403 curves,
five-limb Montgomery multiplication at 367.76 Ir, `ecm12`) prices one side at 2^38.0905
instructions and the other (232 bits) at 2^38.0878, so the per-survivor charge
C_NORMS + C_TDIV + 2^38.0905 + 2^38.0878 is 2^39.0892; the 99.9% figure is the same
ECM_FIND used in the yield's miss model. The higher per-prime success raises
the yield, which lets the region shrink (section 8).
Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.
Survivors get exact norms, trial division below 2^20 and the ECM schedule
above, charged at 2^39.0892 instructions per survivor.

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre in double precision; the remainder over the block disk plus a
2^-32 relative term (coefficient intervals and rounding) is added, and results
are rounded outward. `nbound31` replays the Python leaves bit for bit on all
2,401 pool rows at 1,508.2 Ir per evaluation (charged 1,508.2, the marginal
per-block cost excluding separately charged setups); [2] checks the bounds
against exact norms.

**No block selection.** Every block of every row is enumerated, sieved and
charged at the full per-cell cost; no block is skipped and no yield is dropped.

**Yield ([4]).** rels = (va_region/8744.91) vb (0.971857578) sum over canonical primitive sub-region
rows of w rho(u_f) rho(u_g) r_f r_g 0.99^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
and the ECM miss. Retained relations >= 2^64.1306 against 1.3 U = 2^63.830:
margin >= 0.3007 bits (unrounded 0.30072; pool SE ~0.150), above the 0.30-bit
feasibility floor (`MARGIN`). The descent's missing-log input `LOG2_RELS =
64.130` is the floor (to 0.001 bit) of this retained count (unrounded
2^64.13056; certified lower product already applied). [4] checks `LOG2_RELS <=
log2_rels` (`rels_input` check). The mean-norm inputs 746.0/499.0 are rounded
down from 746.09/499.01 and checked in [4].

**Cost ([5]).** Sampled mean 402.1 instructions/cell; charged mean 429.3
with the allowance below, over 2^126.687 cells: all-cell bounds, layout,
setups, line path, record walk (primes in [2^63, 2^64) charged twice), routing
(10.01/record, 8/cell; STORAGE.md), apply 7.00, recovery 7.003/record, scan,
zero, cofactorisation. Collection: **<= 2^121.1452 units** (estimate 2^121.0510
times allowance 1.0674).

**Canonical-work allowance (H2, [4a]).** The charged work is (va_region/8779.52)
times the per-row sum over in-region canonical rows: a Monte Carlo estimate of
a fixed, challenge-independent domain total. The estimator sums all 20,000
draws (zero outside the sub-region or non-canonical), so its relative SE has a
conditional row-mean term and an incidence term for the 1034 contributing
draws: conditional 1.473%, incidence sqrt((1-1034/20000)/1034) = 3.028%, total
<= 3.368%. The charge multiplies the estimate by **1.0674 >= 1 + 2.00 x 0.03368
= 1.06736** (rounded up); [4a] fails if the allowance falls below 1 + 2.00 SE.
This is a declared modeling margin under H2, not a failure probability, and no
sampling tail is charged against success (6.74% allowance).

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge (coefficients <= 255), simulated on
the production occurrence profile (`filtsim.py`; scale 2^17 table embedded, [6]
reruns 2^12) as the signed sigma-quotient system (each occurrence takes the
sign of its conjugate orientation). Three seeds at 2^17 give N/U <= 0.0471 at
mean weight <= 101 (charged 101/row). N/U falls with scale (`filtsim.run`,
kmax 64, weight 100; max over seeds): 0.0850 (2^12, 10 seeds), 0.0635 (2^13,
6), 0.0563 (2^14, 6), 0.0516 (2^15, 3), 0.0469 (2^16, 3), 0.0461 (2^17, 3),
0.0425 (2^18, 2) (`SCALE_NU`, [10] checks non-increasing). The model charges N
= 1.096 x 0.0471 x U = 2^59.1758: ratio 0.05163 exceeds every measurement from
2^15 up (the base charged 1.61 x 0.0471 to cover 2^12). N/U staying at or below
the 2^15-2^18 values at production scale is part of H3 (monotone-trend
extrapolation). The system solved is the square matrix
A = [M | S]: N + 21 filtered rows (the rest verify), N ideal columns M and 21
dense Schirokauer columns S. Columns are sigma-quotiented, so only sigma-odd
units need maps: signatures K_f (0, 36), K_f^+ (0, 18), K_g (2, 5), K_h (2, 2)
give odd unit ranks 35 - 17 = 18 and 6 - 3 = 3 ([13], Sturm); under H3 the
kernel of A is the line of virtual logs. Block Wiedemann uses random blocks: m
= 19 left, n = 4 right vectors uniform in F_r. The Krylov phase takes N/m + N/n
block steps, i.e. N(n/m + 1) products by A. One product costs per row 101 x
20.7 sparse (STORAGE.md: 4 + 3 + 0.88 x 11 + 0.12 x 33 = 20.64; +-1 entries add
without multiplying, >= 88% of row weight is +-1 by H4u, measured 0.917..0.946;
cost is linear in N x mean weight), 21 dense terms at 97, and row tail C_ROW =
2,710 (counted 2,708; valid for any row, sums < 2^323). Every step forms the m
x n projections x^T A^i y: (m + n) N^2 dense terms at 97 plus one reduction
(C_RED9 = 420, counted 385) per entry. Horner reconstruction (w <- A w + Y F_i
over deg F <= N/m + N/n coefficients) costs one product plus, per element, n
dense terms and one reduction. With c = 101 x 20.7 + 2,710 + 21 x 97 per row:
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n + 420) + 97 (m + n)) + 420 (m + n) N
instructions: **<= 2^117.7103 units** at N = 2^59.1758 (minimum over integer
(m, n) at (19, 4)). The (m + n) x (m + n) generator of the length L = (N +
21)(1/m + 1/n) + 2 sequence (divide-and-conquer approximant bases, Karatsuba):
<= 4 (m + n)^3 L^log2(3) log2(L) field mults at 97 plus reduction 420:
**<= 2^107.1997 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** Region and polynomials are coin/challenge-independent, so the
  count is one fixed number; H1 puts it at >= 2^64.1306, >= 0.3007 bits above
  1.3 U (pool SE ~0.150, so the floor is ~2.0 SE above requirement); the
  descent uses 64.130, at or below that count.
- **Linear algebra.** H3 (standard NFS-DL): r is unramified in K_f and K_g
  (checked by [13]), divides neither class number, the 21 sigma-odd Schirokauer
  maps are injective on the units modulo r-th powers, and the kernel of A is the
  line of virtual logarithms. Block Wiedemann with blocks uniform in F_r then fails
  with probability <= 4N/r < 2^-192 (Schwartz-Zippel on the generator's
  determinant conditions, degree <= 2N in the block entries).
  A failed solve (checked on the dropped rows) counts as failure.
- **Reconstruction.** Purged ideals' logs are recovered by back-substitution.
  The descent counts every ideal with <= 2 occurrences as log-less; [6b]
  measures ideals left without a log on the [6] system, and the model's count
  is >= 12.5x measured in every band with loss (DESCENT.md).

Under H1 and H3 the main phase succeeds with probability >= 1 - 2^-192 over the coins.

## 6. Individual logarithms

DESCENT.md specifies the descent: an ECM-trial boot on balls of certified p^6
lattices, then special-q nodes sieving an exact lattice ball (`dsieve16`) and
recursing on child ideals. Counted cells are base-primitive (0.9100);
survivor-test retention is >= 0.8592 at every sieve node ([7], one Hensel bound
per root below 2^9); cofactors use the `ecm16.py` schedule ([12]). Relation
counts are negative binomial (toy-fit shapes). Log-less ideals are discarded or
re-descended in a fresh lattice. The shipped plan (62 nodes + boot) gives
expected two-log work **<= 2^115.8429 units** ([7]; one log <= 2^129.131
instructions) at the sub-region inputs (coordinate-search plan, UDUP = 0.5).
The meter caps both logs at 2^116.5647 units (cap reached = failure), so
P(cap) <= E[W]/cap <= 2^(115.8429-116.5647) <= 0.6063405 (Markov). A node
counter (cap 2^27.3, section 7) bounds node count by construction.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 2.3892e-04 - 0.6063405 - 2^-9.69 - 2^-10.70 - 2/r)
>= **0.391369** >= 0.39: main phase, two boots
(phi <= 2.3892e-04 each at the sub-region relation input), the instruction-meter
Markov bound, the node-counter Markov bound and special-q repeats, as follows.

**Two explicit meters.** Both logs run under (i) an instruction meter (cap
2^116.5647 units) and (ii) a **special-q node counter** (cap N_MAX = 2^27.3):
the counter increments before each special-q node (ordinary or missing-log)
and the run aborts as a failure past N_MAX. Hence <= 2^27.3 special-q nodes
ever run, by construction.

- Instruction meter: P(cap reached) <= E[W]/cap <= 2^(115.8429 - 116.5647) <=
  0.6063405 (Markov, E[W] the expected two-log work of section 6).
- Node counter: P(abort) <= E[N]/N_MAX (Markov). E[N] uses the same
  finite-tree recurrence as E[W] with unit node cost (`counts16.py`); twice
  the boot value gives E[N] <= 2^17.61 nodes both logs (unrounded 2^17.6098),
  so P(abort) <= 2^(17.61-27.3) = 2^-9.69 <= 0.0012107.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for
  1/p weighting: two nodes share a special-q with probability <= C(N_MAX,2)
  2^6/2^70.30 <= 2^(2x27.3-1+6-70.30), rounded up to 2^-10.70 <= 0.0006012,
  charged in full.

N_MAX = 2^27.3 minimises the sum of the last two terms. No lower bound on
individual node work is used. The boot term treats the 16 relations' subtrees
as independent (H5). Targets are uniform nonidentity by power randomization;
distance from the uniform model is at most 1/r per
target and 2/r is charged above. Identity z is handled separately; the final
pairing check e(P, Q) = z rejects any incorrect Q. The success bound is
conditional on H1-H5, including the H2 canonical-work allowance of section 3,
which is a resource-model margin and not a probability term.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (sample estimate 121.0510 + 6.74% canonical-work allowance) | <= 121.1452 |
| linear algebra | <= 117.7103 |
| matrix generator | <= 107.1997 |
| both individual logarithms (metered cap) | 116.5647 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

The upper-rounded rows sum to 2^121.327036929110, so W <=
**2^121.327036929110** (audit-checked), under the claim; unrounded model total
2^121.3269833.

**Configuration.** Sub-region theta = 759/2000 (1034 of 2401 canonical rows),
exact certified count products at radius >= 892, work allowance 1.0674
(2.00 SE), X2 = 2^64, cofactor slacks (294, 224), relation-collection alpha (5, 0.844),
descent inputs LOG2_RELS = 64.130 and mean norms 746.0/499.0, descent alpha
(5, 1). The 62-node plan is the coordinate-search re-optimisation of
REFINEMENT.md at UDUP = 0.5, THIN_MAX 85, MIN_RELS 18 (35 entries differ from
the base plan, boot unchanged); at present inputs min mu 32.9, max phi < 0.05,
boot psi 0.630 (largest). Two-log expectation <= 2^115.8429, meter 2^116.5647,
plus the section-7 node counter. Norm-bound charge 1,508.2, factor-base bound,
filtering, generator, canonical list, Schirokauer and verification charges are
the base package's; LA charges (20.7/entry, 2,710/row, 420/reduction,
97/dense), blocks (19, 4), pad 1.096 and the 10.01 route are this lineage's
(STORAGE.md). Preprocessing is inside the 2^121.327036929110 cap. No external
advice; nothing amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; relation collection alpha (5, 0.844) against the in-package computed (4.868096, 0.843724) over p <= 2000 and (5.085035, 0.731625) over p < 50000 ([14]); descent alpha (5, 1); exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool at theta = 1/2, uniform sub-region theta = 759/2000, exact flags,
certified sub-region counts with their exact rational products, and the
canonical work mean charged with a 6.74% modeling-margin allowance = 2.00
incidence-inclusive standard errors of the total-work estimator, not a
probability term); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
signed filtering on the production occurrence profile, pad 1.096 over the
2^17 maximum, N/U non-increasing 2^12..2^18 and assumed not to rise at
production); H4 kernel cost (KERNELS.md, STORAGE.md; LA tail/reductions are
statically counted loop-free code, the route an executed-instruction count),
with H4u: >= 88% of production row weight is +-1 (measured 0.917..0.946 at
scales 2^12..2^15, rising; a shortfall changes cost only);
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing logs,
independent subtrees, fresh-lattice re-descents, local omitted-mass laws,
random-integer ECM orders, UDUP = 0.5 (36% below measured 0.779; 0.944 in-run
by [11])).

## 10. Memory

Storage <= 2^105 bytes (STORAGE.md).

## Provenance

Algorithm, models, kernels and package: Gautam Botrel (with Claude Code and
Codex); polynomial selection and alpha follow Guillevic-Singh, filtering
CADO-NFS. Sister-track lean refinements: kaankolcu and pepedesigner (cbab3319
lineage). Sub-region, work-SE formula, certified counts, in-package Murphy
alpha: kaankolcu / DeepSeek lineage (8c1635f2). Port, node counter, plan
re-optimisation, directed rounding: Claude Opus 5.5 in Claude Code (991a7c0d,
1eb4bd1f, 23f293a8, 377ee6d1, e554a0f4). Quadrature-plan base e3940d11
(mpjunior92). Unit-entry LA kernel, H4u, 99.9% ECM / theta = 207/500 port
(recipe 27202622): kaankolcu with Claude Opus 5.5 in Claude Code (04073dd1).
Page-scatter routing and the measured route: kaankolcu's recover-P thread
(be446c45 and its R3 measurement). Reduced LA kernels, N/U scale runs: Claude Opus 5.5 in Claude Code. All grid/plan refinements since (slacks, LA charges, UDUP 0.5, THIN_MAX 85, MIN_RELS 18, joint-move re-search, 10.01 route, meter refits) plus theta 753/2000 and slacks (294, 224): muse-spark-1.3-contributor in Muse Code via T3 Code. This revision (review-driven alpha_g correction with exact primitive p = 2 tails, charge 0.844, theta 759/2000): muse-spark-1.3-contributor in Muse Code via T3 Code. No prior review verdict is a premise.

## Changes from e3940d11 (121.9944)

Base history (04073dd1 lineage to 121.495287837498): 99.9% ECM, theta 207/500,
slacks (300, 228); 10.1/record route, 7.003 recovery; priced LA (row tail,
21/entry, pad 1.10, blocks (19, 4)); 384-cell quadrature, 8 re-optimised plan
entries; X2 2^63->2^64, theta 19/50, slacks (295, 224) (margin 0.3137,
collection 2^121.1557, two-log 2^117.8938, meter 2^118.617).
7. Slacks (292,223) (margin 0.3003, collection 2^121.1444 with route 10.01);
pad 1.096, 20.7/entry, 2710/row (LA 2^117.7103); UDUP 0.5 ([11] 0.944);
THIN_MAX 85, MIN_RELS 18 (boot failure 2.39e-04); joint-move plan re-search
(35 entries; two-log 115.8431); meter 2^116.5664. Total 2^121.326394393321.
10. Alpha_g 0.844->0.60 (computed 0.433521 over p < 50000) was REJECTED on
review: the p = 2 double-root term mu/(q-1) missed the primitivity factor
q/(q+1), overstating v by 0.29810 bits; corrected 50k partial is 0.731625.
11. Alpha_g charge back to 0.844 with corrected exact primitive tails
2q(1-q^-7)/(q^2-1) at p = 2 (15+ tails dropped, conservative), theta
753/2000->759/2000 (1034 rows, a-radius 892, margin 0.3007, allowance
1.0674), slacks stay (294,224) (collection 2^121.1452); inputs 64.130,
746.0/499.0 (two-log 2^115.8429, meter 2^116.5647). Total 2^121.327036929110.

# STNFS over a special tower on a theta = 87/200 sub-region: 2^122.017 reference units

Track bls12-381-fixed-g1-invert-g2. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used. The argument reuses predecessor work as credited below, without assuming a prior-results verdict as a premise.

## Claim

The algorithm below returns Q in G2 with e(P,Q) = z for the exact supplied nonidentity
P in G1 (Q = infinity for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins under the declared heuristics H1-H5;
the modelled success is at least 0.391663 (> 0.39; unrounded 0.39166319...). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^122.017 units**; the score is
**122.017**. The ledger (section 8) totals at most 2^122.01687 (unrounded
2^122.0168682). Relation collection is
charged on the theta = 87/200 sub-region of the sampled region with exact
certified lattice-count products and an explicit 5.49% canonical-work allowance
(2.00 incidence-inclusive standard errors of the sampled work estimator, [4a]). Distributed storage is at most
2^105 bytes (STORAGE.md) in individually bounded, explicitly addressed objects.
The claim is analytic and exploratory under the NFS heuristics H1-H5 (section
9) with the measurements and checks cited; no full-size run is made.

**Numerical conventions.** Decimal figures in this package are the
experiment's values rounded to the digits shown. Wherever a figure is used as a
bound it is written with <= or >= and rounded in the safe direction (upper
bounds up, lower bounds down); the experiment's `declared_bounds` check compares
every such stated bound (ledger, relation collection, linear algebra, relation
count and margin, work SE, two-log expectation, Markov term, E[N], counter-abort
and repeat terms, boot failure, success) with the unrounded value and fails
if any is violated.

The model is `model16.py` and `descent16.py` with `config.json` and `stats16.b85`;
`evidence16.py` recomputes the ledger and the checks cited as [n].

## Base

The GT discrete-logarithm machinery (polynomials, factor base, sieve, filtering,
block Wiedemann, descent, kernels, pools) is the lean STNFS package of the
sister track `bls12-381-fixed-g2-invert-g1` (cbab3319); both pairing
orientations need the same two logarithms in GT. This package uses this track's
reduction (section 1), a theta = 87/200 uniform sub-region of the shipped pool
with exact certified count products and an explicit work allowance (section 3),
relation-collection Murphy alpha (5, 0.844) checked in-package ([14]), an
explicit special-q node counter (section 7) and a descent plan re-optimised for
the sub-region inputs (section 8).

## 1. Reduction

Let g = e(P, Q0) (pinned convention `pairing(Q0, P)`), one pairing, where P is
the supplied point and Q0 the fixed G2 generator of the target profile. P is
nonidentity of prime order r and e is nondegenerate, so g generates mu_r,
z = g^a for a unique a mod r (a = 0 iff z = 1), and Q = [a]Q0 (one G2 scalar
multiplication, about 255 doublings and additions over F_p2), so that
e(P, Q) = e(P, Q0)^a = z for that exact P. The supplied P is never replaced. log z and log g are computed against one shared relation and
linear-algebra phase; a = log z / log g mod r.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks that
(u y)^((p^6-1)/2) = -1 in F_p6, so g is irreducible over F_p6 and
F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. The map
x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P), so
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at
(5, 0.844) in relation collection: [14] recomputes it in-package over p <= 2000
from the exact Kummer-Dedekind structure (`alpha_contrib`, `alpha_sweep`, the
code of this track's 8c1635f2 lineage; `python3 evidence16.py --alpha` prints
(4.868096, 0.545620)); alpha_f is the round-up 5 and alpha_g = 0.844 is the
literature estimate for this polynomial pair, above both the computed 0.545620
and its p <= 5000 value (about 0.46), exactly as charged by 8c1635f2. The
descent keeps the more conservative (5, 1).

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B on both
sides; U = 2^63.451 sigma-orbits (indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21. The pool is sampled uniformly at theta = 1/2
(eligible a: a5 > 0, Q_a(a) <= A^2/2); the collected region is its uniform
sub-region theta = 87/200 (`production_theta_region` in `config.json`): eligible
a: a5 > 0, Q_a(a) <= (87/200) A^2; half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a).
Because the sample is uniform on the theta = 1/2 ellipsoid, its restriction to
qa = Q_a/A^2 <= 87/200 is uniform on the sub-region; the experiment keeps the
1583 canonical sub-region rows (1481 primitive) of the 2401 canonical rows and
scales by the sub-region volume over the effective sample 20000 x (87/100)^3 =
13170.06, which is the same estimator as (va(1/2)/20000) x (sum over sub-region
rows). The canonical a of each unit class
is the least of its exact HNF key class in (Q_a, lex) order, so each relation
orbit of <units, sigma> has exactly one point in the domain; non-primitive pairs
are excluded by an exact HNF test. `lattice_count_bounds.py` (the version of
this track's 8c1635f2 lineage, with a `minimum_radius` argument) certifies from
exact rational Machin-pi and determinant/sqrt intervals the count ratio at every
radius R >= 955 = isqrt((87/200) 2^21), the smallest a-radius of the sub-region
(b-radii are at least sqrt((113/200) 2^21) = 1088): [0.986780, 1.012632] for Q_a
and [0.986737, 1.012632] for Q_b, each inside [0.98, 1.02] (asserted). The model
uses the exact products of these rational bounds, **0.973693682 (yield)** and
**1.025422386 (work)**, not rounded-out factors. Sampling error of the pool is
handled separately and explicitly by the canonical-work allowance below. The
canonical list costs 2^90 units.

**Pool.** 20,000 uniform eligible pairs from seeded indices (`pool12-*.b85`),
2,401 canonical rows with exact norms, primitivity flags and certified
norm-bound statistics; [2] regenerates sampled indices at the sampling region
theta = 1/2 and verifies every field. The collected sub-region uses 1583 of the
2401 rows (1481 primitive).

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^20), record
walk for p in [2^20, 2^63) (`slice12`, KERNELS.md), primes in [2^63, 2^69] taken by ECM
within the slacks (298, 229). The f-side slack is 298, so survivors'
f-side cofactors are bounded by 306 bits; the g-side slack is 229, so survivors'
g-side cofactors are bounded by 237 bits (the slacks are re-optimised for the
sub-region, section 8). The per-survivor ECM charge is the
package's explicit `ecm16.py` schedule at 99% per prime, the same schedule the
descent uses and [12] checks: the cheapest grid point for 69-bit primes in a
306-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 271 curves, five-limb
Montgomery multiplication at 367.76 Ir, `ecm12`) prices one side at 2^37.5179
instructions and the other (237 bits) at 2^37.5152, so the per-survivor charge
C_NORMS + C_TDIV + 2^37.5179 + 2^37.5152 is 2^38.5167 (the predecessor's crude
ECM table charged 2^40.06 per survivor, 2^1.545 above the package's own
schedule; the 99% figure is the same ECM_FIND used in the yield's miss model).
Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.
Survivors get exact norms, trial division below 2^20 and ECM (B1 = 5 x 10^4,
99% per prime, 367.76-instruction Montgomery multiplication, `ecm12`): 2^40.1
instructions per survivor (the schedule's own instruction count at the same 99%).

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre in double precision; the remainder over the block disk plus a
2^-32 relative term (coefficient intervals and rounding) is added, and results
are rounded outward. `nbound31` replays the Python leaves bit for bit on all
all 2,401 pool rows at 1,508.2 Ir per evaluation (charged 1,508.2, the measured marginal
per-block cost excluding the separately charged per-row and per-slice setups); [2] checks the bounds
against exact norms.

**No block selection.** Every block of every row is enumerated, sieved and
charged at the full per-cell cost; no block is skipped and no yield is dropped.

**Yield ([4]).** rels = (va_region/13170.06) vb (0.973693682) sum over canonical primitive sub-region
rows of w rho(u_f) rho(u_g) r_f r_g 0.99^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
and the ECM miss. Retained relations >= 2^64.1319 against 1.3 U = 2^63.830: margin
>= 0.302 bits (unrounded 0.3021; pool SE about 0.131 bits), above the model's 0.30-bit feasibility floor
(`MARGIN` in model16.py). The descent's missing-log input is `LOG2_RELS = 64.131`
in `descent16.py`, the floor (to 0.001 bit) of this model's own retained count
(unrounded 2^64.13199); the certified lower count product is already applied inside that
count. [4] checks `LOG2_RELS <= log2_rels` explicitly (the `rels_input` check).
The descent's relation-weighted mean-norm inputs are 748.3/498.9 bits, rounded
down from the sub-region's computed 748.39/498.93 and checked in [4].

**Cost ([5]).** The sampled mean is 420.9 instructions per cell; the charged
mean, with the canonical-work allowance below, is 444.0 instructions per cell
over 2^127.055 cells: all-cell bounds, layout and setups, the line path, record
walk, routing (36 per record, 8 per cell; STORAGE.md), apply 7.00, recovery
7.003, scan, zero and cofactorisation. Relation collection: **<= 2^121.5619 units**
(sample estimate 2^121.4848 times the allowance 1.0549).

**Canonical-work allowance (H2, [4a]).** The charged work is (va_region/13170.06)
times the sum of per-row work over the in-region canonical rows, i.e. a Monte
Carlo estimate of a fixed, challenge-independent domain total. The exact
lattice-count products bound the population; the work mean over sampled rows
carries its own sampling error. The estimator is a sum over all 20,000 draws
with zero contribution from draws that are outside the sub-region or not
canonical, so its relative standard error has a conditional row-mean term and
an incidence term for the 1583 contributing draws among 20,000. The experiment
computes both from the per-row contributions ([4a]): conditional 1.310%,
incidence sqrt((1 - 1583/20000)/1583) = 2.412%, total 2.745% (<= 2.745%). The relation
collection charge multiplies the sample estimate by **1.0549 >= 1 + 2.00 x
0.02745 = 1.05490** (rounded up), and [4a] fails if the allowance falls below
1 + 2.00 SE. This is a declared modeling margin on a fixed quantity under H2,
not a per-challenge failure probability, and no sampling tail is charged
against the success slack. For comparison, the full-region lean package's
blanket 2601/2500 work bracket is 1.6% above the exact certified product at
radius 1024, i.e. about 0.7 of that region's 2.27% work SE; the sub-region
charge here is a larger explicit 2.00 SE.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge with coefficients bounded by 255,
simulated on the production occurrence profile (`filtsim.py`; scale 2^17 table
embedded, [6] reruns scale 2^12): the simulated matrix is the signed sigma-quotient system (each occurrence takes
the sign of its conjugate orientation, log sigma(P) = -log P). Three seeds at
2^17 give N/U <= 0.0471 at weight <= 101 (charged 101 per row); the model charges
N = 1.61 x 0.0471 x U = 2^59.7306; the charged ratio 1.61 x 0.0471 = 0.07585
still covers the 2^12 values (0.0703 printed by [6], 0.0757 the worst seed;
N/U falls with scale). The system solved is the square matrix
A = [M | S]: N + 21 filtered rows (the rest verify the result), the N
ideal columns M and 21 dense Schirokauer columns S. As columns are sigma-quotiented
(log sigma(P) = -log P), only sigma-odd units need maps: signatures K_f (0, 36),
K_f^+ = K_h(x^2) (0, 18), K_g (2, 5), K_h (2, 2) give odd unit ranks 35 - 17 = 18
and 6 - 3 = 3 ([13], Sturm sequences); under H3 the kernel of A is the line of virtual
logarithms (with their Schirokauer coefficients). Block Wiedemann
on A uses random blocks: m = 26 left and n = 5 right vectors uniform in F_r (the
coins). The Krylov phase takes N/m + N/n
block steps, i.e. N(n/m + 1) products by A; a product costs per row 101 sparse
entries x 88 instructions (STORAGE.md; `la_mac.s` inner loop of 33, [10]) plus 21
dense full-field multiply-adds at 97 (`bw_dot.s`, bw_full: 4 x 4-limb product into
a 9-limb accumulator, inner loop of 97). The 88 charge is exactly the package's
88-instruction per-entry component sum (33 + 4 + 3 + 32 + 16, STORAGE.md) and the
97 charge is exactly the measured bw_full loop; both replace the predecessor's
padded 96 and 104 and this revision's cushioned 90 and 100, and neither exceeds
its measured component. Every step forms the m x n projections
x^T A^i y, (m + n) N^2 full-field terms at 97 ([10]). The kernel vector
is reconstructed by Horner, w <- A w + Y F_i over the deg F <= N/m + N/n generator
coefficients (Y the n starting vectors): one product and n N full-field
multiply-adds per step. With c = 101 x 88 + 21 x 97 per row, the total is
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n) + 97 (m + n)) instructions:
**<= 2^119.3685 units** at N = 2^59.7306 (the minimum of this expression over
integer block sizes is at (m, n) = (26, 5)). The (m + n) x (m + n) matrix
generator of the length
L = (N + 21)(1/m + 1/n) + 2 sequence is computed by divide-and-conquer approximant
bases with Karatsuba products: at most 4 (m + n)^3 L^log2(3) log2(L) full-field
multiplications at 97: **2^106.42 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** The region and polynomials do not depend on coins or the
  challenge, so the relation count is one fixed number; H1 puts it at
  at least 2^64.1319, at least 0.302 bits above 1.3 U (pool standard error
  about 0.131 bits, so the floor is about 2.3 SE above the requirement); the descent uses 64.131, at or
  below that count.
- **Linear algebra.** H3 (standard NFS-DL): r is unramified in K_f and K_g
  (checked by [13]), divides neither class number, the 21 sigma-odd Schirokauer
  maps are injective on the units modulo r-th powers, and the kernel of A is the
  line of virtual logarithms. Block Wiedemann with blocks uniform in F_r then fails
  with probability <= 4N/r < 2^-192 (Schwartz-Zippel on the generator's
  determinant conditions, degree <= 2N in the block entries).
  A failed solve (checked on the dropped rows) counts as failure.
- **Reconstruction.** Logs of ideals removed by singleton purge are recovered by
  back-substitution through the purged relations. The descent counts every
  ideal with at most two occurrences as having no log (section 6); [6b] runs
  the cascading purge and back-substitution on the [6] system and measures the
  ideals left without a log: the model's count is at least 12.5 times the measured
  one in every band with a measured loss (DESCENT.md).

Under H1 and H3 the main phase succeeds with probability >= 1 - 2^-192 over the coins.

## 6. Individual logarithms

DESCENT.md specifies the descent: a boot by ECM trials on balls of certified p^6
lattices, then special-q nodes that sieve an exact lattice ball (`dsieve16`
kernel) and recurse on the child ideals of their relations. Counted cells are
base-primitive (fraction 0.9100); retention through the survivor test follows
from one local Hensel bound per root at every prime ideal below 2^9 (`retain16.py`),
at least 0.8592 at every sieve node ([7]); cofactors use the ECM schedule of
`ecm16.py` ([12]). Relation counts are negative binomial (shapes from toy fits,
DESCENT.md). Ideals without a log are either discarded or descended in a fresh
lattice (DESCENT.md). The shipped plan (62 nodes and the boot) gives expected
two-log work **<= 2^118.1226 units** ([7]; one log about 2^131.410 instructions) at the
sub-region inputs with five plan entries re-optimised (section 8); record layout in STORAGE.md. The descent meter caps both
logs at 2^118.845 units; a run that reaches the cap stops and counts as a
failure, so Markov's inequality gives P(cap reached) <= E[W]/cap <= 0.6060858
(unrounded 0.60608572...), used as an upper bound. A separate special-q node counter (cap
2^27.3 nodes, section 7) bounds the number of nodes by construction.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 1.9830e-04 - 0.6060858 - 2^-9.64 - 2^-10.70 - 2/r)
>= **0.391663** >= 0.39 (unrounded model value 0.39166319...): main phase, two boots
(phi <= 1.9830e-04 each at the sub-region relation input), the instruction-meter
Markov bound, the node-counter Markov bound and special-q repeats, as follows.

**Two explicit meters.** The descent for both logarithms runs under (i) an
instruction meter with cap 2^118.845 units and (ii) a **special-q node counter**
with cap N_MAX = 2^27.3: before starting a special-q node (ordinary or a
missing-log node in the lattice of Q q') the counter is incremented, and the
run aborts, counted as a failure, if it would exceed N_MAX. Hence at most
2^27.3 special-q nodes are ever run, by construction, whatever the node costs.

- Instruction meter: P(cap reached) <= E[W]/cap <= 2^(118.1226 - 118.845) <=
  0.6060858 (Markov, E[W] the expected two-log work of section 6).
- Node counter: P(counter abort) <= E[N]/N_MAX (Markov). E[N] is bounded by the
  same finite-tree recurrence that gives E[W], with every node's own cost
  replaced by 1 (`counts16.py`): C(node) = 1 + 1/(1 - psi) x (sum over child
  slots of the largest C among the bins the child can fall in, plus the
  missing-log children at the largest C of a node with q <= 2^95), and below
  each boot retry x (kA C_max + C_max + miss C_miss). Twice the boot value gives
  E[N] <= 2^17.66 nodes for both logs (unrounded 2^17.6594), so
  P(counter abort) <= 2^(17.66 - 27.3) = 2^-9.64 <= 0.0012534.
- Repeats: with at most N_MAX nodes, placing all of them in the smallest bin
  (2^70.30 ideals) and allowing a factor 2^6 for the 1/p weighting inside a
  bin, two nodes share a special-q ideal with probability
  <= C(N_MAX, 2) 2^6 / 2^70.30 <= 2^(2 x 27.3 - 1 + 6 - 70.30), rounded up to
  2^-10.70 <= 0.0006012, charged in full.

N_MAX = 2^27.3 minimises the sum of the last two terms. No lower bound on the
work of an individual node is used anywhere. The boot term treats the 16 relations' subtrees as
independent (H5). Power randomization gives uniform nonidentity targets for both
logarithms; the distance from the uniform full-group model is at most 1/r per
target and 2/r is charged above. Identity z is handled separately; the final
pairing check e(P, Q) = z rejects any incorrect Q. The success bound is
conditional on H1-H5, including the H2 canonical-work allowance of section 3,
which is a resource-model margin and not a probability term.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (sample estimate 121.4848 + 5.49% canonical-work allowance) | <= 121.5619 |
| linear algebra | <= 119.3685 |
| matrix generator | 106.42 |
| both individual logarithms (metered cap) | 118.845 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

total_log2W_full_ledger <= **2^122.01687** (unrounded 2^122.0168682), under the
122.017 claim. The table's upper-rounded rows sum to 2^122.016881, also below
122.017.

**Configuration.** Sub-region theta = 87/200 (1583 of 2401 canonical rows),
exact certified count products at radius >= 955, work allowance 1.0549
(2.00 SE), cofactor slacks (298, 229), relation-collection alpha (5, 0.844),
descent inputs LOG2_RELS = 64.131 and mean norms 748.3/498.9, descent alpha
(5, 1). The 62-node plan differs from the sister-track plan in five entries
(n, iA, kA, kB, R, thin), found by coordinate search under the plan constraints
(mu >= 16, phi <= 0.05, psi <= 0.8 at every node): (90,f) iA 48 -> 45;
(100,f) n 85 -> 83; (110,f) iA 24 -> 27; (190,g) n 162 -> 161; (210,g) n 170 ->
168. At the present inputs min mu 34.7, max phi 0.0447, max psi 0.490, boot
psi 0.587; the two-log expectation is <= 2^118.1226 and the instruction meter
is 2^118.845, with the node counter of section 7. The norm-bound charge
1,508.2, LA charges 88/97, block sizes (26, 5), pad 1.61, factor-base bound,
filtering, matrix generator, canonical list, Schirokauer maps and verification
charges are those of the base package. Preprocessing per instance is bounded by
the entire 2^122.017-unit work cap and included in it. No external advice is
assumed. Nothing is amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; relation collection alpha (5, 0.844) against the in-package computed (4.868096, 0.545620) over p <= 2000 ([14]); descent alpha (5, 1); exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool at theta = 1/2, uniform sub-region theta = 87/200, exact flags,
certified sub-region counts with their exact rational products, and the
canonical work mean charged with a 5.49% modeling-margin allowance = 2.00
incidence-inclusive standard errors of the total-work estimator, not a
probability term); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
filtering on the production occurrence profile, signed, pad 1.61); H4 kernel cost (KERNELS.md, STORAGE.md);
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing logs, independent
subtree outcomes, including ideals without a log re-descended in fresh lattices,
local omitted-mass distributions, ECM curve orders as random integers).

## 10. Memory and attribution

Storage <= 2^105 bytes (STORAGE.md). Algorithm, models, kernels and package:
Gautam Botrel with Claude Code and Codex; polynomial selection and alpha follow
Guillevic-Singh, filtering CADO-NFS.

## Provenance

Algorithm, models, kernels and package: Gautam Botrel (with Claude Code and
Codex), sister-track lean refinements by kaankolcu and pepedesigner (cbab3319
and its lineage). Sub-region, incidence-inclusive work-SE formula, certified
lattice counts with a minimum radius and the in-package Murphy alpha: this
track's kaankolcu / DeepSeek lineage (8c1635f2 and parents). Port, node counter,
plan re-optimisation, directed-rounding checks and this configuration: Claude
Opus 5.5 in Claude Code (submissions 991a7c0d, 1eb4bd1f, 23f293a8, 377ee6d1,
e554a0f4). Reuse is credited; no prior review verdict is used as a premise.

# GT-core STNFS over a special tower on a theta = 207/500 sub-region with a 99.9% main-phase ECM schedule: core ledger <= 2^121.6752 Fp12 units

**Role.** The GT discrete-logarithm core used by proof.md: the recover-Q argument of e3940d11 with REVISION.md applied in place. Section 1 (the source's G2 output step) is not used; proof.md has the G1 reduction. Unit: dense F_p12 multiplication, C = 20,000 instructions. No advice or amortization.

## Claim

The core below returns a = log_g z for g uniform in GT minus {1} and z independently uniform, with success at least 0.39 over inputs and fresh coins under the declared heuristics H1-H5;
the modelled success is at least 0.391949 (> 0.39; rounded below the model value 0.3919495). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^121.6752 units** (core bound). The ledger (section 8) totals at most 2^121.6752 (unrounded
2^121.6751060). Relation collection is
charged on the theta = 207/500 sub-region of the sampled region with exact
certified lattice-count products and an explicit 5.97% canonical-work allowance
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

The GT discrete-logarithm machinery (polynomials, factor base, sieve, filtering,
block Wiedemann, descent, kernels, pools) is the lean STNFS package of the
sister track `bls12-381-fixed-g2-invert-g1` (cbab3319); both pairing
orientations need the same two logarithms in GT. This package uses this track's
reduction (section 1), a theta = 207/500 uniform sub-region (this revision; was 87/200) of the shipped pool
with exact certified count products and an explicit work allowance (section 3),
relation-collection Murphy alpha (5, 0.844) checked in-package ([14]), an
explicit special-q node counter (section 7) and a descent plan re-optimised for
the sub-region inputs (section 8).

## 1. Reduction

The source's recover-Q reduction is not used here; the ordinary-G1 reduction is in proof.md. The core computes a = log z / log g mod r for a uniform nonidentity base g and an independently uniform z, with log z and log g against one shared relation and linear-algebra phase.

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
its uniform sub-region theta = 207/500 (`production_theta_region`): a5 > 0,
Q_a(a) <= (207/500) A^2, half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The
experiment keeps the 1378 canonical sub-region rows (1285 primitive) of the 2401
and scales by the sub-region volume over the effective sample 20000 x
(207/250)^3 = 11353.27. The canonical a of each unit class is the least of its
exact HNF key class in (Q_a, lex) order, so each orbit of <units, sigma> has one
point in the domain; non-primitive pairs are excluded by an exact HNF test.
`lattice_count_bounds.py` (8c1635f2 lineage) certifies from exact rational pi
and determinant intervals the count ratio at every radius R >= 931 =
isqrt((207/500) 2^21) (b-radii >= 1108): [0.986780, 1.012632] (Q_a) and
[0.986737, 1.012632] (Q_b). The model uses their exact products **0.973023183
(yield)** and **1.026085506 (work)**; pool sampling error is handled by the
allowance below. The canonical list costs 2^90 units. **Pool:** 20,000 seeded
uniform pairs (`pool12-*.b85`), 2,401 canonical rows with exact norms, flags and
bound statistics; [2] regenerates and verifies them.

**Sieve and cofactorisation.** Line sieve for p in [2^8, 2^20), record walk for
p in [2^20, 2^63) (`slice12`, KERNELS.md), primes in [2^63, 2^69] by ECM within
slacks (300, 228): survivors' cofactors are at most 308 (f) and 236 (g) bits.
Certified per-block bounds (`embedding_bounds.py`, 64-cell blocks of width 128,
bisected until both bounds are within 8 bits; Taylor shift per place in double
precision plus the disk remainder and a 2^-32 relative term, rounded outward;
`nbound31` replays the leaves bit for bit at 1,508.2 Ir per evaluation, charged
1,508.2; [2] checks against exact norms) give thresholds T = lower - 16 x slack.
Survivors get exact norms, trial division below 2^20 and the `ecm16.py`
schedule at **99.9% per prime** (R1; [12]): B1 = 2^16, B2 = 30 B1, wheel 30030,
403 curves, 367.76 Ir per multiplication (`ecm12`): 2^38.0908 + 2^38.0880,
with C_NORMS + C_TDIV 2^39.0895 per survivor. **No block selection:** every
block of every row is enumerated, sieved and charged.

**Yield ([4]).** rels = (va_region/11353.27) vb 0.973023183 sum over canonical
primitive sub-region rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f +
lambda_g): Dickman rho ([1]), retentions from `threshold_model.py` (calibrated on
exact smooth integers, B = 2^10..2^14) and the ECM miss. Retained relations >=
2^64.1310 against 1.3 U = 2^63.830: margin >= 0.3012 bits (pool SE about 0.137),
above the 0.30-bit floor (`MARGIN`). The descent input `LOG2_RELS = 64.131` is
the floor of this count and its mean norms 747.5/498.9 are rounded down from
747.60/498.90; [4] checks both (`rels_input`).

**Cost ([5]).** Charged 388.7 instructions per cell over
2^126.940 cells: bounds, layout and setups, line path, record walk, routing (10.3
per record, R2 + R4; 8 per cell), apply 7.00, recovery 7.003, scan, zero and
cofactorisation. Relation collection **<= 2^121.2551 units** (sample estimate
2^121.1714 x allowance 1.0597).

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo
estimate of a fixed domain total, a sum over all 20,000 draws with zero
contribution outside the sub-region or non-canonical, so its relative SE has a
conditional term (1.462%) and an incidence term sqrt((1 - 1378/20000)/1378) =
2.599%: total 2.982%. The charge multiplies by **1.0597 >= 1 + 2.00 x 0.02983**
and [4a] fails below 1 + 2.00 SE. This is a declared modelling margin under H2,
not a probability term.

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
on A uses random blocks: m = 18 left and n = 4 right vectors uniform in F_r (the
coins). The Krylov phase takes N(n/m + 1) products by A; a product costs per row
101 sparse entries at C_ENTRY = 44 (STORAGE.md: unit/general group loops 11/33,
combined row tail 2,090 per row, unit share >= 0.88 under H4u; R3) plus 21 dense
full-field multiply-adds at 97 (`bw_dot.s`). Every step forms the m x n
projections, (m + n) N^2 full-field terms at 97 ([10]); Horner reconstruction
w <- A w + Y F_i over N/m + N/n steps costs one product and n N full-field terms
per step. With c = 101 x 44 + 21 x 97, the total N^2 (c (n/m + 1) + (1/m +
1/n)(c + 97 n) + 97 (m + n)) is **<= 2^118.7428 units** at N = 2^59.7306 (the
integer minimum is at (m, n) = (18, 4)). The generator of the length L = (N +
21)(1/m + 1/n) + 2 sequence (divide-and-conquer approximant bases, Karatsuba):
at most 4 (m + n)^3 L^log2(3) log2(L) multiplications at 97: **2^105.5088 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** The relation count is one fixed number (region and
  polynomials are coin- and challenge-independent); H1 puts it at >= 2^64.1310,
  0.3012 bits above 1.3 U (about 2.2 pool SE); the descent uses 64.131.
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

DESCENT.md: a boot by ECM trials on balls of certified p^6 lattices, then
special-q nodes sieving exact lattice balls (`dsieve16`) and recursing on child
ideals. Counted cells are base-primitive (0.9100); retention >= 0.8592 at every
sieve node ([7], `retain16.py`); cofactors by `ecm16.py` ([12]); negative
binomial relation counts. The 62-node plan and the boot give expected two-log
work **<= 2^117.9122 units** ([7]; one log at most 2^131.200 instructions). The
meter caps both logs at 2^118.635 units (a capped run fails; Markov <= 0.60589)
and a node counter (section 7) bounds the node count by construction.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 1.7864e-04 - 0.60589 - 2^-9.70 - 2^-10.70 - 2/r)
>= **0.391949** >= 0.39: main phase, two boots
(phi <= 1.7864e-04 each at the sub-region relation input), the instruction-meter
Markov bound, the node-counter Markov bound and special-q repeats, as follows.

**Two explicit meters** (DESCENT.md, Repeated ideals): an instruction meter
(cap 2^118.635 units; Markov P <= 2^(117.9122 - 118.635) <= 0.60589) and a
special-q node counter (cap N_MAX = 2^27.3 nodes, incremented before every
ordinary or missing-log node, abort counted as failure). E[N] <= 2^17.60 for
both logs by the unit-cost recurrence (`counts16.py`; unrounded 2^17.5979), so
P(counter abort) <= 2^-9.70; with at most N_MAX nodes in the smallest bin
(2^70.30 ideals, factor 2^6 for the 1/p weighting) a special-q repeat has
probability <= 2^(2 x 27.3 - 1 + 6 - 70.30) <= 2^-10.70, charged in full.

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
| relation collection (sample estimate 121.1714 + 5.97% canonical-work allowance) | <= 121.2551 |
| linear algebra | <= 118.7428 |
| matrix generator | <= 105.5088 |
| both individual logarithms (metered cap) | 118.635 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

W <= **2^121.6752** (unrounded 2^121.6751060), the core bound used by proof.md. The table's upper-rounded rows, including a 2^50-unit numerical-certificate allowance, sum to 2^121.675125, below the stated core bound 2^121.6752.

**Configuration.** As sections 3-7 (REVISION.md lists every change from the
base package); plan min mu 36.3, max phi 0.0354, max psi 0.478, boot psi 0.583.
Preprocessing is inside the 2^121.6752-unit cap; no advice, nothing amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; relation collection alpha (5, 0.844) against the in-package computed (4.868096, 0.545620) over p <= 2000 ([14]); descent alpha (5, 1); exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool at theta = 1/2, uniform sub-region theta = 207/500, exact flags,
certified sub-region counts with their exact rational products, and the
canonical work mean charged with a 5.49% modeling-margin allowance = 2.00
incidence-inclusive standard errors of the total-work estimator, not a
probability term); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
filtering on the production occurrence profile, signed, pad 1.61); H4 kernel cost (KERNELS.md, STORAGE.md), H4u unit share >= 0.88 (STORAGE.md);
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing logs, independent
subtree outcomes, including ideals without a log re-descended in fresh lattices,
local omitted-mass distributions, ECM curve orders as random integers).

## 10. Memory

Storage <= 2^105 bytes (STORAGE.md). Polynomial selection and alpha follow
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

## Directed quadrature refinement

Inherited; DESCENT.md (Refined integration). Changed numbers: REVISION.md.

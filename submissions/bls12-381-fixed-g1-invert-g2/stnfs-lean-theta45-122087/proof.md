# STNFS over a special tower on a theta = 9/20 sub-region: 2^122.087 reference units

Track bls12-381-fixed-g1-invert-g2. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used. The argument reuses predecessor work as credited below, without assuming a prior-results verdict as a premise.

## Claim

The algorithm below returns Q in G2 with e(P,Q) = z for the exact supplied nonidentity
P in G1 (Q = infinity for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins under the declared heuristics H1-H5;
the modelled success is 0.3915832 (directed 4-decimal floor 0.3915 > 0.39). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^122.087 units**; the score is
**122.087**. The ledger (section 8) totals 2^122.086307. Relation collection is
charged on the theta = 9/20 sub-region of the sampled region with exact
certified lattice-count products and an explicit 5.25% canonical-work allowance
(2.00 incidence-inclusive standard errors of the sampled work estimator, [4a]). Distributed storage is at most
2^105 bytes (STORAGE.md) in individually bounded, explicitly addressed objects.
The claim is analytic and exploratory under the NFS heuristics H1-H5 (section
9) with the measurements and checks cited; no full-size run is made.

The model is `model16.py` and `descent16.py` with `config.json` and `stats16.b85`;
`evidence16.py` recomputes the ledger and the checks cited as [n].

## Base: port of the sister-track lean package

This package starts from `stnfs-lean-port-122206` (this track, submission
991a7c0d-60fa-4ce3-8430-01c8bfd10db1, claimed 122.2065; rejected because its
special-q node count divided the meter by the modelled expected node cost, which
is not a lower bound on executed work -- corrected here in section 7), itself a
port, to
`bls12-381-fixed-g1-invert-g2`, of the sister-track
candidate `stnfs-lean-122206-v1` (submission cbab3319-552d-4143-955c-e7fc8f3056de,
track `bls12-381-fixed-g2-invert-g1`, claimed 122.2065, AI-review valid). The
number field sieve works entirely in GT, the order-r subgroup of F_p12^*: both
orientations of fixed-argument pairing inversion reduce to the same two
discrete logarithms in GT against one shared relation-collection and
linear-algebra phase, so the polynomial selection, factor base, sieve, filtering,
block Wiedemann, descent plan, kernels and pools are those of the sister
package, and the ledger changes only as listed in section 8 ("Changes relative
to the port"). In the port, only the reduction (section 1: g =
e(P, Q0) for the supplied P, output Q = [a]Q0 in G2), the output scalar
multiplication (now in G2, still inside the 2^11 verification-and-output row)
and the track name change. This track's own accepted lineage (Gautam Botrel's
220103a4 at 122.35, 020720b9 at 122.34, 3b4bd7dc at 122.32, 0240d166 at 122.29
and 8c1635f2 at 122.253) uses the same STNFS algorithm and reduction; the
sister-track lean lineage differs from it only in the numerical charges listed
in section 8 (explicit `ecm16.py` cofactor schedule, measured 1,508.2 / 88 / 97
kernel charges, slacks (291, 222), descent input 64.224 and cap 2^118.886) and
in using the full theta = 1/2 region with the blanket 2601/2500 and 2401/2500
count brackets instead of a theta = 455/1000 sub-region with a sampling allowance.
This revision (section 3 and "Changes relative to the port" at the end) moves
the lean package onto a theta = 9/20 uniform sub-region of its own pool, with the
exact certified count products of this track's lineage and an explicit,
measured canonical-work allowance in place of the blanket brackets.

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
(5, 1) (rounded up from 4.868, 0.844).

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B on both
sides; U = 2^63.451 sigma-orbits (indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21. The pool is sampled uniformly at theta = 1/2
(eligible a: a5 > 0, Q_a(a) <= A^2/2); the collected region is its uniform
sub-region theta = 9/20 (`production_theta_region` in `config.json`): eligible
a: a5 > 0, Q_a(a) <= (9/20) A^2; half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a).
Because the sample is uniform on the theta = 1/2 ellipsoid, its restriction to
qa = Q_a/A^2 <= 9/20 is uniform on the sub-region; the experiment keeps the
1748 canonical sub-region rows (1630 primitive) of the 2401 canonical rows and
scales by the sub-region volume over the effective sample 20000 x (9/10)^3 =
14580, which is the same estimator as (va(1/2)/20000) x (sum over sub-region
rows). The canonical a of each unit class
is the least of its exact HNF key class in (Q_a, lex) order, so each relation
orbit of <units, sigma> has exactly one point in the domain; non-primitive pairs
are excluded by an exact HNF test. `lattice_count_bounds.py` (the version of
this track's 8c1635f2 lineage, with a `minimum_radius` argument) certifies from
exact rational Machin-pi and determinant/sqrt intervals the count ratio at every
radius R >= 971 = isqrt((9/20) 2^21), the smallest a-radius of the sub-region
(b-radii are at least sqrt((11/20) 2^21) = 1074): [0.986997, 1.012422] for Q_a
and [0.986955, 1.012422] for Q_b, each inside [0.98, 1.02] (asserted). The model
uses the exact products of these rational bounds, **0.974122467 (yield)** and
**1.024998723 (work)**, not rounded-out factors. Sampling error of the pool is
handled separately and explicitly by the canonical-work allowance below. The
canonical list costs 2^90 units.

**Pool.** 20,000 uniform eligible pairs from seeded indices (`pool12-*.b85`),
2,401 canonical rows with exact norms, primitivity flags and certified
norm-bound statistics; [2] regenerates sampled indices at the sampling region
theta = 1/2 and verifies every field. The collected sub-region uses 1748 of the
2401 rows (1630 primitive).

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^20), record
walk for p in [2^20, 2^63) (`slice12`, KERNELS.md), primes in [2^63, 2^69] taken by ECM
within the slacks (294, 225). The f-side slack is 294, so survivors'
f-side cofactors are bounded by 302 bits; the g-side slack is 225, so survivors'
g-side cofactors are bounded by 233 bits (the slacks are re-optimised for the
sub-region, section 8). The per-survivor ECM charge is the
package's explicit `ecm16.py` schedule at 99% per prime, the same schedule the
descent uses and [12] checks: the cheapest grid point for 69-bit primes in a
302-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 271 curves, five-limb
Montgomery multiplication at 367.76 Ir, `ecm12`) prices one side at 2^37.5177
instructions and the other (233 bits) at 2^37.5150, so the per-survivor charge
C_NORMS + C_TDIV + 2^37.5177 + 2^37.5150 is 2^38.5164 (the predecessor's crude
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

**Yield ([4]).** rels = (va_region/14580) vb (0.974122467) sum over canonical primitive sub-region
rows of w rho(u_f) rho(u_g) r_f r_g 0.99^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
and the ECM miss. Retained relations 2^64.146 against 1.3 U = 2^63.830: margin
0.316 bits (pool SE 0.128 bits), above the model's 0.30-bit feasibility floor
(`MARGIN` in model16.py). The descent's missing-log input is `LOG2_RELS = 64.145`
in `descent16.py`, the floor (to 0.001 bit) of this model's own retained count
2^64.1460; the certified lower count product is already applied inside that
count. [4] checks `LOG2_RELS <= log2_rels` explicitly (the `rels_input` check).
The descent's relation-weighted mean-norm inputs are 748.8/498.9 bits, rounded
down from the sub-region's computed 748.9/498.98 and checked in [4].

**Cost ([5]).** The sampled mean is 417.7 instructions per cell; the charged
mean, with the canonical-work allowance below, is 439.6 instructions per cell
over 2^127.134 cells: all-cell bounds, layout and setups, the line path, record
walk, routing (36 per record, 8 per cell; STORAGE.md), apply 7.00, recovery
7.003, scan, zero and cofactorisation. Relation collection: **2^121.6262 units**
(sample estimate 2^121.5524 times the allowance 1.0525).

**Canonical-work allowance (H2, [4a]).** The charged work is (va_region/14580)
times the sum of per-row work over the in-region canonical rows, i.e. a Monte
Carlo estimate of a fixed, challenge-independent domain total. The exact
lattice-count products bound the population; the work mean over sampled rows
carries its own sampling error. The estimator is a sum over all 20,000 draws
with zero contribution from draws that are outside the sub-region or not
canonical, so its relative standard error has a conditional row-mean term and
an incidence term for the 1748 contributing draws among 20,000. The experiment
computes both from the per-row contributions ([4a]): conditional 1.283%,
incidence sqrt((1 - 1748/20000)/1748) = 2.285%, total 2.621%. The relation
collection charge multiplies the sample estimate by **1.0525 >= 1 + 2.00 x
0.02621 = 1.05241** (rounded up), and [4a] fails if the allowance falls below
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
**2^119.3684 units** at N = 2^59.7306 (the minimum of this expression over
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
  2^64.146, 0.316 bits above 1.3 U (pool standard error 0.128 bits, so the
  floor is about 2.5 SE above the requirement); the descent uses 64.145, at or
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
two-log work **2^118.253 units** ([7]; one log 2^131.541 instructions) at the
sub-region inputs; record layout in STORAGE.md. The descent meter caps both
logs at 2^119.04 units; a run that reaches the cap stops and counts as a
failure, so Markov's inequality gives P(cap reached) <= E[W]/cap = 0.5795091 =
2^-0.787097, used as an upper bound.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 1.7523e-04 - 0.5795091 - 2^-5.13 - 2/r)
= 0.3915832 (directed 4-decimal floor **0.3915** >= 0.39): main phase, two boots
(phi = 1.7523e-04 each at the sub-region relation input), the descent-cap Markov
upper bound E[W]/cap = 2^(118.2529 - 119.04) = 0.5795091, and special-q repeats
(charged as 2^-5.13 = 0.0285572). **Node-count bound.** The number of special-q
nodes the meter allows is bounded by the cap divided by a *deterministic* lower
bound on the work every node executes: LLL reduction (2^36) plus the fixed sieve
work over the node's whole ball (all cells bounded, sieved, scanned and zeroed,
all pair setups), which does not depend on how many cells survive. Random
survivor cofactorisation is excluded from this floor. The smallest floor in the
plan is 2^103.24 instructions (the R = 192 node), so at most
2^(119.04 + log2 20000 - 103.24) = 2^30.09 nodes run in the two-log descent;
placing all of them in the smallest bin (2^70.3 ideals, factor 2^6 for the 1/p
weighting) gives P(two nodes share a special-q) <= 2^(2 x 30.09 - 1 + 6 - 70.3)
<= 2^-5.13. (An earlier version divided by the smallest modelled *expected* node
cost, 2^103.76, which includes survivor work and is not a lower bound on executed
work; that gave 2^29.53 nodes and an unsupported 2^-6.24 allowance. The meter
is raised from 2^119.0 to 2^119.04 to keep success >= 0.39 under the corrected
count.) The boot term treats the 16 relations' subtrees as
independent (H5). Power randomization gives uniform nonidentity targets for both
logarithms; the distance from the uniform full-group model is at most 1/r per
target and 2/r is charged above. Identity z is handled separately; the final
pairing check e(P, Q) = z rejects any incorrect Q. The success bound is
conditional on H1-H5, including the H2 canonical-work allowance of section 3,
which is a resource-model margin and not a probability term.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (sample estimate 121.5524 + 5.25% canonical-work allowance) | 121.6262 |
| linear algebra | 119.3684 |
| matrix generator | 106.42 |
| both individual logarithms (metered cap) | 119.04 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

total_log2W_full_ledger = **122.086307**, under the 122.087 claim. The displayed
rounded component exponents sum to 2^122.086304, also below 122.087.

**Changes relative to the port (stnfs-lean-port-122206, 122.2065).**
(i) The collected region is the theta = 9/20 uniform sub-region of the shipped
theta = 1/2 pool (1748 of 2401 canonical rows), so the eligible-a volume falls
by (9/10)^3 while the dropped a's (the largest Q_a, hence the largest norms per
cell) contributed the least yield per cell; population 2^127.134 cells instead
of 2^127.389.
(ii) The blanket 2601/2500 and 2401/2500 count brackets are replaced by the
exact certified products at radius >= 971 (1.024998723 work, 0.974122467
yield), and the pool's sampling error on the work estimator is charged
explicitly as a 2.00-SE allowance (1.0525, [4a]).
(iii) The cofactor slacks are re-optimised from (291, 222) to (294, 225) for the
sub-region; the ECM charge stays the explicit `ecm16.py` schedule at 99%.
(iv) The descent inputs follow the sub-region model: LOG2_RELS = 64.145 (model
2^64.1460), mean-norm inputs 748.8/498.9 (model 748.9/498.98); with the same
62-node plan and boot (215, 43, 15, 4) the two-log expectation is 2^118.253 and
the meter is set to 2^119.04 (node-count bound from the deterministic per-node
work floor, section 7).
The certified norm-bound charge 1,508.2, the linear-algebra charges 88/97, the
block sizes (26, 5), the pad 1.61, the factor-base bound, filtering, matrix
generator, canonical list, Schirokauer maps and verification charges are
unchanged. The preprocessing per instance is conservatively bounded by the
entire 2^122.087-unit work cap and included in it. No external advice is
assumed. Nothing is amortized.

Sweep behind the choice (same model, 2.00-SE allowance, descent cap re-chosen
for success >= 0.3901, using the earlier expected-cost node count; the
corrected count adds about 0.005 to each): theta = 0.48 gives 122.175, 0.47
gives 122.145, 0.46 gives 122.115, 9/20 gives 122.0815 at margin 0.316; theta <= 0.44 needs the
relation margin to fall below the model's 0.30-bit floor with these slacks and
was not used. At theta = 1/2 the same explicit 2.00-SE treatment would give
122.236, i.e. the gain comes from the sub-region, not from loosening the
sampling treatment.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (5, 1), exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool at theta = 1/2, uniform sub-region theta = 9/20, exact flags,
certified sub-region counts with their exact rational products, and the
canonical work mean charged with a 5.25% modeling-margin allowance = 2.00
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

## Provenance of this revision

This theta = 9/20 revision (Claude Opus 5.5 in Claude Code) combines the lean
sister-track package (via the port 991a7c0d) with the sub-region idea,
`production_theta_region` handling, `lattice_count_bounds.py` with a minimum
radius and the incidence-inclusive work-SE formula of this track's
kaankolcu / DeepSeek lineage (8c1635f2 at 122.253 and its parents 0240d166,
3b4bd7dc), using a 2.00-SE allowance with the point-estimate relation input
of the lean lineage instead of that lineage's 3.50-SE allowance and 1.25-SE
relation haircut. The provenance of the underlying lean package follows.

This revision is derived from this lane's candidate `stnfs-lean-122224-v1`
(submitted as 0ac63018-97d5-4b08-a774-f9f666926f02, claimed 122.224, ledger
2^122.223379, validating when this revision was written), which fixed the
binding rejection of 97b80666-f55c-4f14-b3cc-c0ddd6f6ae9a: that earlier tree's
descent used LOG2_RELS = 64.162 while its relation-collection model supplied
only 2^64.1325 retained relations (2.1% more relations than the model, in the
direction that understates missing logs). The 122.224 tree fixed the input to
the model's own floored count and supplied the measured marginal charges
1,508.2 (norm bound) and 88/97 (linear algebra), the pad 1.61, the boot
(215, 43, 15, 4), the block sizes (26, 5) and the accepted frontier's work and
yield brackets; this revision keeps all of those and re-derives LOG2_RELS from
its own model, `LOG2_RELS = 64.224 = floor(2^64.2243)`, with [4] checking the
inequality. The lineage of that tree is this lane's
submission fc9fd9a0-b86e-4f33-933e-64f3166bc750 (candidate
`stnfs-lean-122238-v1`, claimed 122.238, ledger 2^122.237394, official 122.238),
the accepted frontier, which supplies the pad 1.61, the boot (215, 43, 15, 4)
and the corrected Markov text; and the accepted lean frontier submission
aaa530bf-f89d-470e-b12e-eee78d0b12d5 (candidate `stnfs-lean-12233-v3`, branch
22e75a29-3b29-40a0-b8ae-7a6fbf6e046e) via this lane's accepted revision
a57c5623-932c-4630-ad5c-58785cda23d0 (score 122.3188), which supplies the
metered descent cap method, and a73abc31-b6e0-46e7-8786-041770d22e6c (score
122.305), which supplies the block sizes (26, 5) and the 2^119.39 cap lineage.
The rejected lattice-count re-accounting ac905f62 is not used: the accepted
frontier's blanket 2601/2500 and 2401/2500 factors are retained and only the
cofactor slacks, the descent relation input and the measured charges are
retuned. The frontier itself reuses Gautam Botrel's fixed-G1 submission
220103a4-b71a-4a2f-afbe-36ad6782d2fc, the derivative
020720b9-4317-4cf0-8234-92cae794251e and its fixed-G2 source
3cb69d49-624c-4385-a52f-4ce81a2427a6. The changes in this revision are numerical
re-optimisations and one metered-charge correction inside the accepted model:
(i) the per-survivor ECM cofactor charge is the package's own explicit
`ecm16.py` schedule at 99% (2^38.5163) instead of the crude ECM table
(2^40.06); (ii) the ECM cofactor slacks are retuned to (291, 222); (iii) the
descent input is the model's own floored relation count 64.224; (iv) the
metered cap is 2^118.886 with modelled success 0.3912 > 0.391. The declared
ceiling falls from the accepted 122.238 (and from this lane's validating 122.224
submission) to 122.2065. It makes no claim of a new algorithm or
production-scale execution.

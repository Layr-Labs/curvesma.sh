# STNFS over a special tower: 2^122.264 reference units

Track bls12-381-fixed-g1-invert-g2. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used. The argument reuses predecessor work as credited below, without assuming a prior-results verdict as a premise.

## Claim

The algorithm below returns Q with e(P,Q) = z for the exact supplied nonidentity
P (Q = infinity for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins; the modelled success is 0.3924 (bound > 0.392). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^122.264 units**; the score is
**122.264**. The ledger (section 8) totals 2^122.263311, and it includes an 11.69%
canonical-work allowance on relation collection (4.50 incidence-inclusive measured
standard errors of the total-work estimator, section 3). Distributed storage is at most
2^105 bytes (STORAGE.md) in individually bounded, explicitly addressed objects.
The claim is analytic and exploratory under the NFS heuristics H1-H5 (section
9) with the measurements and checks cited; no full-size run is made.

The model is `model16.py` and `descent16.py` with `config.json` and `stats16.b85`;
`evidence16.py` recomputes the ledger and the checks cited as [n].

## 1. Reduction

Let g = e(P, Q0) (pinned convention `pairing(Q0, P)`), one pairing. P is
nonidentity of prime order r and e is nondegenerate, so g generates mu_r,
z = g^a for a unique a mod r (a = 0 iff z = 1), and Q = [a]Q0 (one G2 scalar
multiplication, about 255 doublings and additions over F_p2). log z and log g are computed against one shared relation and
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
(`config.json`); A^2 = 2^21, theta = 455/1000 = 0.455; eligible a: a5 > 0,
Q_a(a) <= theta A^2; half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The pool is
sampled uniformly at theta = 1/2; the collected region is the uniform
sub-region qa = Q_a/A^2 <= 455/1000 of that sample (1806 of 2401 canonical rows,
1681 primitive; effective sample 15071.420 of 20000). Since the sample is
uniform in the ellipsoid, restricting it to the sub-region is the correct
estimator for the smaller region; the code uses the sub-region volume and the
sub-region sample count, which is the same estimator. The canonical a of each unit class
is the least of its exact HNF key class in (Q_a, lex) order, so each relation
orbit of <units, sigma> has exactly one point in the domain; non-primitive pairs
are excluded by an exact HNF test. `lattice_count_bounds.py` certifies counts
within [0.98, 1.02] of the volumes at every radius >= 976 (the sub-region
radius sqrt(theta A^2) = 2^9.931, rounded down); the exact products of the
per-region rational bounds are 0.974253611 (yield) and 1.024869209 (work). The
model uses those products, not rounded-out factors. The canonical list costs 2^90 units.

**Pool.** 20,000 uniform eligible pairs from seeded indices (`pool12-*.b85`),
2,401 canonical rows with exact norms, primitivity flags and certified
norm-bound statistics; [2] regenerates sampled indices at the sampling region
theta = 1/2 and verifies every field. The collected sub-region uses 1806 rows
(1681 primitive).

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^20), record
walk for p in [2^20, 2^63) (`slice12`, KERNELS.md), primes in [2^63, 2^69] taken by ECM
within the slacks (292, 220). Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.
Survivors get exact norms, trial division below 2^20 and ECM (B1 = 5 x 10^4,
99% per prime, 367.76-instruction Montgomery multiplication, `ecm12`): 2^40.1
instructions per survivor (offline, the explicit schedule of `ecm16.py` prices 69-bit
primes in a 300-bit cofactor at 99% at 2^37.4, below the charge).

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre in double precision; the remainder over the block disk plus a
2^-32 relative term (coefficient intervals and rounding) is added, and results
are rounded outward. `nbound31` replays the Python leaves bit for bit on all
2,401 pool rows at 1,508 Ir per evaluation (charged 1,550); [2] checks the bounds
against exact norms.

**No block selection.** Every block of every row is enumerated, sieved and
charged at the full per-cell cost; no block is skipped and no yield is dropped.

**Yield ([4]).** rels = (va_region/n_eff) vb (0.974253611) sum over canonical primitive
rows of w rho(u_f) rho(u_g) r_f r_g 0.99^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
and the ECM miss. Retained relations 2^64.130 against 1.3 U = 2^63.830: margin
0.3004 bits (pool SE 0.127), which is above the package's pre-existing
`MARGIN = 0.30` feasibility floor.

**Cost ([5]).** The raw sampled mean is 426.2 instructions per cell; the charged
relation collection applies the canonical-work allowance below, i.e. 476.0
instructions per cell over 2^127.159 cells: all-cell bounds, layout and setups,
the line path, record walk, routing (36 per record, 8 per cell; STORAGE.md),
apply 7.00, recovery 7.003, scan, zero and cofactorisation. Relation
collection: **2^121.7660 units**.

**Canonical-work allowance (H2).** The charged work is
`geo.va/n_total` times a sum of per-row work over the in-region rows. This is a
Monte Carlo estimate of the fixed canonical-domain work, not a geometric bound:
the exact lattice-count factors bound the population, while the work mean over
the sampled rows carries its own sampling error. The estimator sums a
contribution over all 20,000 pool draws with zero contribution outside the
collected sub-region, so its relative standard error has two terms: the
conditional row-mean term and the incidence term for the 1806 contributing
draws among 20,000. Directly from the contributions, the conditional term is
1.307% and the incidence term sqrt((1 - 1806/20000)/1806) = 2.244%, giving a
total relative standard error of 2.597%; the experiment computes
these from the row contributions in [4a], it does not assume them. The charge
multiplies the sample mean by **1.1169** = 1 + 4.50 total standard errors
(1 + 4.50 x 0.02597 = 1.11687, rounded up), i.e. 11.69%. Under the normal
approximation for the 20000-draw mean the probability that the estimator is
short by more than 4.50 total standard errors is 3.4e-6; this is the failure
probability charged for the H2 allowance in section 7, and the variance is
measured, not assumed. If the fixed domain work exceeded the metered allowance
the run would stop and the attempt would fail. The allowance is for the
finite-sample estimate of a challenge-independent fixed domain, not a
probability over the algorithm's coins, and it is not used as an algorithmic
success term.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge with coefficients bounded by 255,
simulated on the production occurrence profile (`filtsim.py`; scale 2^17 table
embedded, [6] reruns scale 2^12): the simulated matrix is the signed sigma-quotient system (each occurrence takes
the sign of its conjugate orientation, log sigma(P) = -log P). Three seeds at
2^17 give N/U <= 0.0471 at weight <= 101 (charged 101 per row); the model charges
N = 1.65 x 0.0471 x U = 2^59.7660, which covers the 2^12 value (0.0750 printed
by [6], 0.0757 the worst seed; N/U falls with scale). The system solved is the square matrix
A = [M | S]: N + 21 filtered rows (the rest verify the result), the N
ideal columns M and 21 dense Schirokauer columns S. As columns are sigma-quotiented
(log sigma(P) = -log P), only sigma-odd units need maps: signatures K_f (0, 36),
K_f^+ = K_h(x^2) (0, 18), K_g (2, 5), K_h (2, 2) give odd unit ranks 35 - 17 = 18
and 6 - 3 = 3 ([13], Sturm sequences); under H3 the kernel of A is the line of virtual
logarithms (with their Schirokauer coefficients). Block Wiedemann
on A uses random blocks: m = 26 left and n = 5 right vectors uniform in F_r (the
coins). The Krylov phase takes N/m + N/n
block steps, i.e. N(n/m + 1) products by A; a product costs per row 101 sparse
entries x 96 instructions (STORAGE.md; `la_mac.s` inner loop of 33, [10]) plus 21
dense full-field multiply-adds at 104 (`bw_dot.s`, bw_full: 4 x 4-limb product into
a 9-limb accumulator, inner loop of 97). Every step forms the m x n projections
x^T A^i y, (m + n) N^2 full-field terms at 104 ([10]). The kernel vector
is reconstructed by Horner, w <- A w + Y F_i over the deg F <= N/m + N/n generator
coefficients (Y the n starting vectors): one product and n N full-field
multiply-adds per step. With c = 101 x 96 + 21 x 104 per row, the total is
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 104 n) + 104 (m + n)) instructions:
**2^119.5568 units**. The (m + n) x (m + n) matrix generator of the length
L = (N + 21)(1/m + 1/n) + 2 sequence is computed by divide-and-conquer approximant
bases with Karatsuba products: at most 4 (m + n)^3 L^log2(3) log2(L) full-field
multiplications at 104: **2^106.57 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** The region and polynomials do not depend on coins or the
  challenge, so the relation count is one fixed number; H1 puts it at
  2^64.130, 0.3004 bits above 1.3 U (pool standard error 0.127 bits).
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
lattice (DESCENT.md). The collected sub-region has relation-weighted mean norms
749.1/499.0 bits; the descent model is run with the conservative lower inputs
2^64.13 relations and means 749.0/498.8 bits. The shipped plan (62 nodes and the
boot) gives expected two-log work **2^118.634 units** ([7]); record layout in
STORAGE.md. The descent meter caps both logs at 2^119.41 units; a run that
reaches the cap stops and counts as a failure, so Markov gives
P(cap reached) <= 2^-0.776. The cap is lowered from 2^119.455 to 2^119.41
because the theta = 455/1000 region and the exact certified count factors lower
the charged relation collection enough to pay for the 0.026 success cost, and
the success bound keeps 0.002 of slack above 0.39.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 1.09e-04 - 2^-0.776 - 2^-5.42) > **0.392** >= 0.39
(main phase, two boots, descent cap, repeated ideals). The repeated-ideal allowance is
2^-5.42 at the 2^119.41 cap because the special-q nodes reachable under the
meter are 2^29.94. The boot term treats the 16 relations' subtrees
as independent (H5). Power randomization of the targets makes both
individual-log tasks distributed as the sampled production targets for every P
and z; the final pairing check rejects any incorrect Q.

The relation-collection charge carries the 11.69% canonical-work allowance of
section 3 (4.50 incidence-inclusive standard errors of the total-work
estimator). Under H2 the pool is a uniform sample, so the allowance bounds the
fixed canonical-domain work relative to the sample estimate except with
normal-tail probability 3.4e-6; if the domain work exceeded it, the global meter
would stop the run and the attempt would fail. That 3.4e-6 loss is smaller than
the 0.002 slack of the success bound above 0.39, so the combined bound stays
above 0.392. The allowance is a modeling margin, not a coin-flip success term,
and the success bound above is the algorithmic success probability conditional
on H1-H5 including this H2 allowance.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (sample estimate 121.6065 + 11.69% canonical-work allowance) | 121.7660 |
| linear algebra | 119.5568 |
| matrix generator | 106.57 |
| both individual logarithms (metered cap) | 119.42 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

total_log2W_full_ledger = **122.263311**, under the 122.264 claim. This revision
retunes the collection region theta 1/2 -> 455/1000 (with the exact certified count
products at radius 976, 0.974253611 yield and 1.024869209 work), charges an 11.69%
canonical-work allowance (4.50 incidence-inclusive standard errors, computed in
the experiment) on the sampled relation-collection mean, sets block-Wiedemann
block sizes m = 26, n = 5, and sets the two-log descent cap at 2^119.41, above
the 2^118.634 expected two-log work. Every other component,
the algorithm, the kernels, the pool and the heuristics are unchanged. The
preprocessing per instance is conservatively bounded by the entire 2^122.264-unit
work cap and included in it. No external advice is assumed. Nothing is amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (5, 1), exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool at theta = 1/2, uniform sub-region theta = 455/1000, exact flags,
certified sub-region counts with their exact rational products, and the
canonical work mean charged with an 11.69% allowance = 4.50 incidence-inclusive
standard errors of the total-work estimator); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
filtering on the production occurrence profile, signed, pad 1.65); H4 kernel cost (KERNELS.md, STORAGE.md);
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing logs, independent
subtree outcomes, including ideals without a log re-descended in fresh lattices,
local omitted-mass distributions, ECM curve orders as random integers).

## 10. Memory and attribution

Storage <= 2^105 bytes (STORAGE.md). Algorithm, models, kernels and package:
Gautam Botrel with Claude Code and Codex; polynomial selection and alpha follow
Guillevic-Singh, filtering CADO-NFS.

## Provenance of this revision

This revision reuses the complete STNFS algorithm, files, data and cost model of
this lane's accepted 122.29 submission 0240d166-365c-4c05-8f28-aef558397ac5,
its parents 020720b9-4317-4cf0-8234-92cae794251e (122.34) and
3b4bd7dc-12fa-4c9c-a7ce-a6515a6d5824 (122.32), Gautam Botrel's
220103a4-b71a-4a2f-afbe-36ad6782d2fc, and the sister-track submissions
3cb69d49-624c-4385-a52f-4ce81a2427a6 and 4b4ccd00-3fdd-4b25-a10e-cb5ecd114c4b.
Relative to 0240d166 the changes are: the collection region is retuned from
theta = 23/50 to theta = 455/1000 on the same A^2 = 2^21 lattice (a uniform
sub-region of the certified pool), the certified lattice-count factors are
recomputed at the sub-region radius 976 (0.974253611 / 1.024869209), the
incidence-inclusive work-mean standard error is recomputed at the new region
(2.597%, charged 1.1169 = 1 + 4.50 SE, 11.69%), the block-Wiedemann block sizes
are set to m = 26, n = 5, and the two-log descent cap is lowered to 2^119.41
with the conservative descent input mean 749.0 bits (the experiment recomputes
the boot, the special-q plans and the success bound, 0.3924). Relative to
c71e60b3-088f-46ac-a758-415b2d046994 (claim 122.270) this revision lowers the
region theta from 457/1000 to 455/1000 and the metered cap from 2^119.42 to
2^119.41; the certified count products and the allowance are recomputed at the
new radius, and the relation margin becomes 0.3004 bits (still above the 0.30
floor). The algorithm,
kernels, pool, filtering and linear algebra are otherwise unchanged. No new
cryptanalytic technique is claimed.

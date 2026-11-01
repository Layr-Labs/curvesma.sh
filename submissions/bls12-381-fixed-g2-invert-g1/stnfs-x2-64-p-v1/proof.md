# Recover P (fixed-G2 pairing inversion into G1): GT STNFS on a theta = 19/50 sub-region, 99.9% ECM, X2 = 2^64, measured record route, fully reduced LA kernels: 2^121.4953522 reference units

Track bls12-381-fixed-g2-invert-g1. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used. The argument reuses predecessor work as credited below, without assuming a prior-results verdict as a premise.

## Claim

The algorithm below returns P in G1 with e(P,Q) = z for the exact supplied nonidentity
Q in G2 (P = infinity exactly for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins under the declared heuristics H1-H5;
the modelled success is at least 0.392133 (> 0.39; rounded below the model value 0.3921522 minus the 2^-48 field-map allowance). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, the field-basis conversion, both pairings, verification and output, by **2^121.4953522 units**; the score is
**121.4953522**. The ledger (section 8): upper-rounded rows sum to 2^121.49535220 (unrounded
model total 2^121.4953516). Relation collection is
charged on the theta = 19/50 sub-region of the sampled region with exact
certified lattice-count products and an explicit 6.72% canonical-work allowance
(2.00 incidence-inclusive standard errors of the sampled work estimator, [4a]). Distributed storage is at most
2^105 bytes (STORAGE.md) in individually bounded, explicitly addressed objects.
The claim is analytic and exploratory under the NFS heuristics H1-H5 (section
9) with the measurements and checks cited; no full-size run is made.

**Numerical conventions.** Figures used as bounds are written with <= or >= and rounded
in the safe direction; `declared_bounds` in evidence16.py compares each with its unrounded
value and fails on any violation.

The model is `model16.py` and `descent16.py` with `config.json` and `stats16.b85`;
`evidence16.py` recomputes the ledger and the checks cited as [n].

## Base

The GT discrete-logarithm machinery (polynomials, factor base, sieve, filtering, block Wiedemann, descent, kernels, pools) is the lean STNFS package of this track (cbab3319), refined on the sister track `bls12-381-fixed-g1-invert-g2` up to 711418d3 (121.4953); both pairing orientations need the same two logarithms in GT, and only the base g, the output group and the wrapper differ. This package ports 711418d3's core unchanged except for upward-rounded kernel prices (section 8) and adds this track's reduction, field-basis conversion and wrapper. Changes: last section.

## 1. Reduction (recover P)

1. Decode and check Q and z as in the acceptance rule (Q on the G2 curve, Q != infinity,
   [r]Q = infinity; z != 0, z^r = 1). If z = 1, output P = infinity and stop.
2. Compute g = e(P0, Q) (pinned convention `pairing(Q, P0)`), P0 the fixed G1 generator.
3. Map g and z into the core tower by the field-basis conversion below.
4. Run the GT core (sections 2-7) on (g, z) under its meters; any core failure is a failure.
5. With a = log z / log g mod r, compute P = [a]P0, check e(P, Q) = z on all 12
   coefficients (second pairing) and output P canonically; otherwise fail.

G2 has prime order r and Q is uniform nonidentity, so Q = [b]Q0 with b uniform in
{1, ..., r-1} and g = e(P0, Q0)^b is uniform on GT minus {1}, independent of the uniform z:
exactly the input distribution of the core, with no re-randomization. g generates GT, so
z = g^a for a unique a mod r and e([a]P0, Q) = g^a = z; P is unique because e(., Q) is a
bijection G1 -> GT. A field isomorphism preserves z = g^a whichever roots it selects. The
supplied Q is never replaced and no bare scalar is returned. log z and log g are computed
against one shared relation and linear-algebra phase.

**Field-basis conversion.** Let F be the protected F_(p^12) encoding. Factor h(Y) = Y^6 +
Y^3 + Y^2 - Y - 1 in F[Y] into linear factors and choose a root y (h is irreducible of
degree 6 over F_p and 6 | 12, so it splits in F); then find a root x in F of X^2 + |u| y
(irreducible over F_(p^6) by [0b], so it splits in F with x outside that subfield).
Y -> y, X -> x is a field isomorphism from the core tower onto F. The 12 vectors y^i x^j
form an F_p basis; invert their 12 x 12 coordinate matrix by Gaussian elimination and map g
and z with it. Splitting a squarefree completely split f of degree d <= 6: pick R of degree
< d with uniform F coefficients and take gcd(R, f) and gcd(R^((p^12-1)/2) - 1, f); values at
distinct roots are jointly uniform (Vandermonde), so each trial splits with probability
> 1/4. Five splits for h and one for the quadratic, each capped at 128 trials: failure
<= 6 (3/4)^128 < 2^-50; uniform F elements by 381-bit rejection sampling (128 attempts per
coordinate, < 2^16 coordinates: < 2^-112). Field-map failure allowance: 2^-48 (section 7).
Cost: < 2^36 F_p12 operations at < 2^31 instructions each (< 2^67 instructions), < 2^30 bytes.

**Wrapper charge.** Input checks, two pairings (each < 2^50 instructions), the field map,
input conversion, scalar arithmetic mod r, [a]P0, subgroup and encoding checks and output:
at most 2^80 instructions, a deliberately loose ceiling, charged as the ledger row
2^80/20000 < 2^65.72 units.

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
sub-region theta = 19/50 (`production_theta_region` in `config.json`): eligible
a: a5 > 0, Q_a(a) <= (19/50) A^2; half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a).
Because the sample is uniform on the theta = 1/2 ellipsoid, its restriction to
qa = Q_a/A^2 <= 19/50 is uniform on the sub-region; the experiment keeps the
1040 canonical sub-region rows (965 primitive) of the 2401 canonical rows and
scales by the sub-region volume over the effective sample 20000 x (87/100)^3 =
8779.52, which is the same estimator as (va(1/2)/20000) x (sum over sub-region
rows). The canonical a of each unit class
is the least of its exact HNF key class in (Q_a, lex) order, so each relation
orbit of <units, sigma> has exactly one point in the domain; non-primitive pairs
are excluded by an exact HNF test. `lattice_count_bounds.py` (the version of
this track's 8c1635f2 lineage, with a `minimum_radius` argument) certifies from
exact rational Machin-pi and determinant/sqrt intervals the count ratio at every
radius R >= 892 = isqrt((19/50) 2^21), the smallest a-radius of the sub-region
(b-radii are at least sqrt((293/500) 2^21) = 1108): [0.986780, 1.012632] for Q_a
and [0.986737, 1.012632] for Q_b, each inside [0.98, 1.02] (asserted). The model
uses the exact products of these rational bounds, **0.971857578 (yield, 9-decimal floor)** and
**1.027240121 (work, 9-decimal ceiling)**, not rounded-out factors. Sampling error of the pool is
handled separately and explicitly by the canonical-work allowance below. The
canonical list costs 2^90 units.

**Pool.** 20,000 uniform eligible pairs from seeded indices (`pool12-*.b85`),
2,401 canonical rows with exact norms, primitivity flags and certified
norm-bound statistics; [2] regenerates sampled indices at the sampling region
theta = 1/2 and verifies every field. The collected sub-region uses 1040 of the
2401 rows (965 primitive).

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^20), record
walk for p in [2^20, 2^64) (`slice12`, KERNELS.md; X2 = 2^64, was 2^63), primes in [2^64, 2^69] taken by ECM
within the slacks (295, 224). The f-side slack is 295, so survivors'
f-side cofactors are bounded by 303 bits; the g-side slack is 224, so survivors'
g-side cofactors are bounded by 232 bits (the slacks are re-optimised for the
sub-region, section 8). The per-survivor ECM charge is the
package's explicit `ecm16.py` schedule at **99.9%** per prime (ECM_FIND = 0.999, raised from the
base's 99%; this is also the `ecm16.schedule` default the descent uses), which [12] checks: the cheapest
grid point for 69-bit primes in a 303-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 403 curves,
five-limb Montgomery multiplication at 367.76 Ir, `ecm12`) prices one side at 2^38.0906
instructions and the other (232 bits) at 2^38.0878, so the per-survivor charge
C_NORMS + C_TDIV + 2^38.0906 + 2^38.0878 is 2^39.0893; the 99.9% figure is the same
ECM_FIND used in the yield's miss model. Cofactorisation is still only a few instructions per cell;
the higher per-prime success raises the yield, which lets the region shrink (section 8).
Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.
Survivors get exact norms, trial division below 2^20 and the ECM schedule above
(99.9% per prime, 367.76-instruction Montgomery multiplication, `ecm12`), charged at
2^39.0893 instructions per survivor.

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre in double precision; the remainder over the block disk plus a
2^-32 relative term (coefficient intervals and rounding) is added, and results
are rounded outward. `nbound31` replays the Python leaves bit for bit on all
all 2,401 pool rows at 1,508.2 Ir per evaluation (charged 1,508.2, the measured marginal
per-block cost excluding the separately charged per-row and per-slice setups); [2] checks the bounds
against exact norms.

**No block selection.** Every block of every row is enumerated, sieved and
charged at the full per-cell cost; no block is skipped and no yield is dropped.

**Yield ([4]).** rels = (va_region/8779.52) vb (0.971857578) sum over canonical primitive sub-region
rows of w rho(u_f) rho(u_g) r_f r_g 0.99^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
and the ECM miss. Retained relations >= 2^64.1435 against 1.3 U = 2^63.830: margin
>= 0.3137 bits (unrounded 0.31371; pool SE about 0.149 bits), above the model's 0.30-bit feasibility floor
(`MARGIN` in model16.py). The descent's missing-log input is `LOG2_RELS = 64.143`
in `descent16.py`, the floor (to 0.001 bit) of this model's own retained count
(unrounded 2^64.14356); the certified lower count product is already applied inside that
count. [4] checks `LOG2_RELS <= log2_rels` explicitly (the `rels_input` check).
The descent's relation-weighted mean-norm inputs are 746.0/498.9 bits, rounded
down from the sub-region's computed 746.1/499.0 and checked in [4].

**Cost ([5]).** The sampled mean is 403.8 instructions per cell; the charged
mean, with the canonical-work allowance below, is 430.9 instructions per cell
over 2^126.692 cells: all-cell bounds, layout and setups, the line path, record
walk (primes in [2^63, 2^64) charged twice, STORAGE.md), routing (10.1 per record, measured page scatter and transfer, 8 per cell; STORAGE.md), apply 7.00, recovery
7.003 per record, scan, zero and cofactorisation. Relation collection: **<= 2^121.155782 units**
(sample estimate 2^121.06195 times the allowance 1.0672; upward-rounded prices, section 8).

**Canonical-work allowance (H2, [4a]).** The charged work is (va_region/8779.52)
times the sum of per-row work over the in-region canonical rows, i.e. a Monte
Carlo estimate of a fixed, challenge-independent domain total. The exact
lattice-count products bound the population; the work mean over sampled rows
carries its own sampling error. The estimator is a sum over all 20,000 draws
with zero contribution from draws that are outside the sub-region or not
canonical, so its relative standard error has a conditional row-mean term and
an incidence term for the 1040 contributing draws among 20,000. The experiment
computes both from the per-row contributions ([4a]): conditional 1.471%,
incidence sqrt((1 - 1040/20000)/1040) = 3.019%, total 3.359% (<= 3.359%). The relation
collection charge multiplies the sample estimate by **1.0672 >= 1 + 2.00 x
0.03359 = 1.06718** (rounded up), and [4a] fails if the allowance falls below
1 + 2.00 SE. This is a declared modeling margin on a fixed quantity under H2,
not a per-challenge failure probability, and no sampling tail is charged
against the success slack.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge with coefficients bounded by 255,
simulated on the production occurrence profile (`filtsim.py`; scale 2^17 table
embedded, [6] reruns scale 2^12): the simulated matrix is the signed sigma-quotient system (each occurrence takes
the sign of its conjugate orientation, log sigma(P) = -log P). Three seeds at
2^17 give N/U <= 0.0471 at mean weight <= 101 (charged 101 per row). N/U falls with
scale: the maximum over seeds of `filtsim.run` (kmax 64, mean weight 100) is
0.0850 (2^12, 10 seeds), 0.0635 (2^13, 6), 0.0563 (2^14, 6), 0.0516 (2^15, 3),
0.0469 (2^16, 3), 0.0461 (2^17, 3), 0.0425 (2^18, 2) (`SCALE_NU`, checked
non-increasing by [10]; `unit_entries.py` reports U/N at any scale up to 2^15).
The model charges N = 1.10 x 0.0471 x U = 2^59.1811: the charged ratio 0.05182 is
above every measured value from 2^15 up (the base charged 1.61 x 0.0471 to cover
2^12). That N/U at production scale stays at or below the 2^15-2^18 values is
part of H3, an extrapolation of a monotone measured trend. The system solved is the square matrix
A = [M | S]: N + 21 filtered rows (the rest verify the result), the N
ideal columns M and 21 dense Schirokauer columns S. As columns are sigma-quotiented
(log sigma(P) = -log P), only sigma-odd units need maps: signatures K_f (0, 36),
K_f^+ = K_h(x^2) (0, 18), K_g (2, 5), K_h (2, 2) give odd unit ranks 35 - 17 = 18
and 6 - 3 = 3 ([13], Sturm sequences); under H3 the kernel of A is the line of virtual
logarithms (with their Schirokauer coefficients). Block Wiedemann
on A uses random blocks: m = 19 left and n = 4 right vectors uniform in F_r (the
coins). The Krylov phase takes N/m + N/n block steps, i.e. N(n/m + 1) products by
A. A product costs per row 101 x 21 for the sparse entries (STORAGE.md: 4 + 3 +
0.88 x 11 + 0.12 x 33 = 20.64; +-1 entries are added without multiplication, at
least 88% of the row weight is +-1, heuristic H4u, measured 0.917..0.946 by
`unit_entries.py`; the cost is linear in the total entry count, N x mean weight),
21 dense terms at 97 (bw_full: a 4 x 4-limb product added into the row's 9-limb
accumulator) and the row tail C_ROW = 2,750 (counted 2,708): six-limb partial-sum
messages up the 20-owner tree, combine, and the finish that reduces the sparse
pair and the dense accumulator mod r (fold9 + Barrett reduce6), valid for any row
(sums < 2^323). Every step forms the m x n projections x^T A^i y, (m + n) N^2
dense terms at 97 plus one reduction (C_RED9 = 420, counted 385) per projection
entry. The kernel vector is reconstructed by Horner, w <- A w + Y F_i over the
deg F <= N/m + N/n generator coefficients (Y the n starting vectors): one product
plus, per element, n dense terms and one reduction. With c = 101 x 21 + 2,750 +
21 x 97 per row, the total is N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n + 420) +
97 (m + n)) + 420 (m + n) N instructions: **<= 2^117.7327 units** at
N = 2^59.1811 (minimum over integer block sizes at (m, n) = (19, 4)). The
(m + n) x (m + n) matrix generator of the length L = (N + 21)(1/m + 1/n) + 2
sequence is computed by divide-and-conquer approximant bases with Karatsuba
products: at most 4 (m + n)^3 L^log2(3) log2(L) field multiplications, each 97
plus a reduction 420: **<= 2^107.2082 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** The region and polynomials do not depend on coins or the
  challenge, so the relation count is one fixed number; H1 puts it at
  at least 2^64.1435, at least 0.3137 bits above 1.3 U (pool standard error
  about 0.149 bits, so the floor is about 2.0 SE above the requirement); the descent uses 64.143, at or
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
two-log work **<= 2^117.8938 units** ([7]; one log at most 2^131.181 instructions) at the
sub-region inputs with the plan re-optimised for the refined quadrature (section 8); record layout in STORAGE.md. The descent meter caps both
logs at 2^118.617 units; a run that reaches the cap stops and counts as a
failure, so Markov's inequality gives P(cap reached) <= E[W]/cap <= 2^(117.8938 - 118.617) <= 0.6057524,
used as an upper bound. A separate special-q node counter (cap
2^27.3 nodes, section 7) bounds the number of nodes by construction.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 1.6385e-04 - 0.6057524 - 2^-9.72 - 2^-10.70 - 2/r)
>= **0.392133** >= 0.39: main phase, two boots
(phi <= 1.6385e-04 each at the sub-region relation input), the instruction-meter
Markov bound, the node-counter Markov bound and special-q repeats, as follows.

**Two explicit meters.** The descent for both logarithms runs under (i) an
instruction meter with cap 2^118.617 units and (ii) a **special-q node counter**
with cap N_MAX = 2^27.3: before starting a special-q node (ordinary or a
missing-log node in the lattice of Q q') the counter is incremented, and the
run aborts, counted as a failure, if it would exceed N_MAX. Hence at most
2^27.3 special-q nodes are ever run, by construction, whatever the node costs.

- Instruction meter: P(cap reached) <= E[W]/cap <= 2^(117.8938 - 118.617) <=
  0.6057524 (Markov, E[W] the expected two-log work of section 6).
- Node counter: P(counter abort) <= E[N]/N_MAX (Markov). E[N] is bounded by the
  same finite-tree recurrence that gives E[W], with every node's own cost
  replaced by 1 (`counts16.py`): C(node) = 1 + 1/(1 - psi) x (sum over child
  slots of the largest C among the bins the child can fall in, plus the
  missing-log children at the largest C of a node with q <= 2^95), and below
  each boot retry x (kA C_max + C_max + miss C_miss). Twice the boot value gives
  E[N] <= 2^17.58 nodes for both logs (unrounded 2^17.5771), so
  P(counter abort) <= 2^(17.58 - 27.3) = 2^-9.72 <= 0.0011858.
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
pairing check e(P, Q) = z rejects any incorrect P. The field map's 2^-48 failure
allowance is subtracted as well: 0.3921330 - 2^-48 >= 0.392133. The success bound is
conditional on H1-H5, including the H2 canonical-work allowance of section 3,
which is a resource-model margin and not a probability term.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (sample estimate 121.0620 + 6.72% canonical-work allowance) | <= 121.155782 |
| linear algebra | <= 117.732694 |
| matrix generator | <= 107.208126 |
| both individual logarithms (metered cap) | 118.617 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| recover-P wrapper (checks, two pairings, field map, [a]P0, output; 2^80 instructions) | 65.72 |

The table's upper-rounded rows, including a 2^50-unit numerical-certificate allowance, sum to 2^121.49535220 (80-digit decimal, `audit_refinement.py`), so W <= **2^121.4953522**, the claim; the unrounded model total is 2^121.4953516.

**Kernel prices rounded up.** measurements12.txt reports the line and record path prices as differential ratios printed to 0.01. The sister package rounded them to nearest, so five of the nine were below the logged ratio (worst: line p = 4099, 3,828,670 / 123,056 = 31.1132, charged 31.11). Here every table entry is the logged ratio rounded up: LINE_TAB (24.10, 31.12, 88.32, 310.71), REC_TAB (44.91, 78.22, 266.98, 938.09, 2231.66). Each price is now at or above its measured ratio, and relation collection moves from 2^121.155702 to 2^121.155782.

**Configuration.** Sub-region theta = 19/50 (1040 of 2401 canonical rows),
exact certified count products at radius >= 892, work allowance 1.0672
(2.00 SE), X2 = 2^64, cofactor slacks (295, 224), relation-collection alpha (5, 0.844),
descent inputs LOG2_RELS = 64.143 and mean norms 746.0/498.9, descent alpha
(5, 1). The 62-node plan is the coordinate-search plan of REFINEMENT.md and refinement.json (mu >= 16, phi <= 0.05, psi <= 0.8 at every node), unchanged here. At the present inputs
min mu 36.6, max phi below 0.05, boot psi 0.580 (the largest psi); the two-log
expectation is <= 2^117.8938 and the instruction meter is 2^118.617, with the node counter of section 7. The norm-bound charge
1,508.2, factor-base bound,
filtering, matrix generator, canonical list, Schirokauer maps and verification
charges are those of the base package; LA charges (21 per entry, 2,750 per row, 420 per reduction, 97 per dense term), block sizes (19, 4), pad 1.10 and the 10.1 record route are this revision's (STORAGE.md). Preprocessing per instance is bounded by
the entire 2^121.4953522-unit work cap and included in it. No external advice is
assumed. Nothing is amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; relation collection alpha (5, 0.844) against the in-package computed (4.868096, 0.545620) over p <= 2000 ([14]); descent alpha (5, 1); exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool at theta = 1/2, uniform sub-region theta = 19/50, exact flags,
certified sub-region counts with their exact rational products, and the
canonical work mean charged with a 6.72% modeling-margin allowance = 2.00
incidence-inclusive standard errors of the total-work estimator, not a
probability term); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
filtering on the production occurrence profile, signed, pad 1.10 over the 2^17 maximum, N/U non-increasing from 2^12 to 2^18 and assumed not to rise at production scale); H4 kernel cost (KERNELS.md, STORAGE.md; the LA row tail and reductions are compiled loop-free code counted statically, the record route is an executed-instruction count), with H4u: at least 88% of the production filtered matrix's row weight is +-1 entries (measured 0.917..0.946 at filtering scales 2^12..2^15 with the package's own filtering simulation, rising with scale; a shortfall changes cost only, not correctness);
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing logs, independent
subtree outcomes, including ideals without a log re-descended in fresh lattices,
local omitted-mass distributions, ECM curve orders as random integers).

## 10. Memory

Storage <= 2^105 bytes (STORAGE.md).

## Provenance

Algorithm, models, kernels, package: Gautam Botrel (Claude Code, Codex); polynomials and
alpha after Guillevic-Singh, filtering after CADO-NFS. Refinements: kaankolcu and
pepedesigner (cbab3319), kaankolcu/DeepSeek (8c1635f2: sub-region, work SE, certified
counts, alpha), Claude Opus 5.5 in Claude Code (991a7c0d, 1eb4bd1f, 23f293a8, 377ee6d1,
e554a0f4), mpjunior92 (e3940d11), kaankolcu with Claude Opus 5.5 in Claude Code (27202622,
04073dd1, be446c45 and R3 route, 711418d3 reduced LA and X2 = 2^64, 871a7c6c field map,
115f3fd9, and this port). No prior review verdict is a premise.

## Changes

Against e3940d11 (121.9944), as in 711418d3: (1) 99.9% main-phase ECM; (2) page-scatter
record route measured at 10.1 per record; (3) every LA reduction mod r priced, row tail
2,750, 21 per entry, pad 1.61 -> 1.10, blocks (19, 4); (4) the REFINEMENT.md quadrature
and plan; (5) X2 = 2^64, slacks (295, 224), theta 19/50, primes in [2^63, 2^64) charged
twice; (6) every record re-read for recovery (no gating).

Against this track's 115f3fd9 (121.7877): the core above replaces the R1-R4 core; the
1.01 reserve (which covered two-decimal prices rounded down) is replaced by rounding every
LINE_TAB/REC_TAB entry up (section 8, +8e-5 bits); the reduction, field map and a 2^80
instruction wrapper row are section 1; rows rounded up at six decimals: 2^121.4953522.

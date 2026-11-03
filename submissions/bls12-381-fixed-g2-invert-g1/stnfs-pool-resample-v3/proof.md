# STNFS on a directly sampled 32,000-draw pool: 2^121.400000 reference units

Track bls12-381-fixed-g2-invert-g1. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used.
Predecessor work is reused as credited below; no prior verdict is a premise.

## Claim

The algorithm below returns P in G1 with e(P,Q) = z for the exact supplied
nonidentity Q in G2 (P = infinity for z = 1) with success at least 0.39 over
uniform full-production challenges and fresh coins under heuristics H1-H5; the
modelled success is >= 0.392624. A global instruction meter bounds all work,
including preprocessing, failures, restarts, transfers, verification and
output, by **2^121.400000 units**; the score is **121.400000**. Storage is at
most 2^105 bytes (STORAGE.md). The claim is analytic and exploratory; no
full-size run is made.

**Numerical conventions.** Bounds are rounded in the safe direction; the
`declared_bounds` check of `evidence16.py` (the initial experiment) compares
every stated bound with the unrounded value. The model is `model16.py` and
`descent16.py` with `config.json` and `stats16.b85`; `evidence16.py`
recomputes the ledger and the checks cited as [n]. The initial run uses `quick`, which
lists [6], [6b], [11], [11b] and [12] as skipped; `quick=false` runs them (about 30 s).

## Base

Algorithm, kernels, descent plan and charges are those of submission 6282ab79
(121.315000; parent 3a00ffc8) and this track's 3cb69d49. The relation-collection
region theta = 2/5 is sampled directly by a 32,000-draw pool (seeds 1 and 2,
3,907 canonical rows; margin 0.4766), with work allowance 1.0334, descent inputs
64.306 and 744.2/499.2, and meter 2^116.3915, all recomputed by `evidence16.py`.

## 1. Reduction

Let g = e(P0, Q) (pinned `pairing(Q, P0)`), one pairing; P0 is the fixed G1
generator. Q has prime order r and e is nondegenerate, so g generates mu_r,
z = g^a for unique a mod r (a = 0 iff z = 1), and P = [a]P0 (one G1 scalar
mult). log z and log g share one relation/LA phase; a = log z / log g mod r.
g and z enter the tower of section 2 through the standard field isomorphism
(roots of h and of X^2 + |u| y in F_p12, then a 12 x 12 basis change), included
with the pairings and output in the 2^66-unit verification row.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks that
(u y)^((p^6-1)/2) = -1 in F_p6, so g is irreducible over F_p6 and
F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. The map
x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P), so
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at
(5, 0.85) in collection and (5, 1) in the descent; [14] recomputes it from
the exact Kummer-Dedekind structure: (4.868096, 0.843643) over p <= 2000 and
0.731544 for g over p < 50000 (primitive-correct p = 2 term).

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69;
the factor base is the relative-degree-1 prime ideals over K_h of norm <= B;
U = 2^63.451 sigma-orbits.

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21. Eligible a: a5 > 0, Q_a(a) <= (2/5) A^2;
half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The pool is uniform on this region,
so the estimator scales by the region volume over the 32,000 draws. The
canonical a of each unit class is the least of its exact HNF key class in
(Q_a, lex) order, so each relation orbit of <units, sigma> has one point in
the domain; non-primitive pairs are excluded by an exact HNF test.
`lattice_count_bounds.py` certifies the count ratio at every radius
R >= 915 = isqrt((2/5) 2^21) (each form inside [0.98, 1.02], asserted); the
model uses the exact products **0.972556868 (yield, floor)** and
**1.026547145 (work, ceiling)**. The canonical list costs 2^90 units.

**Pool.** 32,000 uniform eligible pairs from seeded indices (seed 1: 2,000,
seed 2: 30,000; `pool12-*.b85`), 3,907 canonical rows (3,589 primitive) with
exact norms, primitivity flags and certified norm-bound statistics; [2]
regenerates the first `samples` stored seed-1 rows and seed-1 indices
0..samples-1 and verifies every field.

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^20), record
walk for p in [2^20, 2^64) (`slice12`, KERNELS.md), primes in [2^64, 2^69] by ECM
within slacks (295, 224). Certified per-block bounds (64-cell blocks of width
128, bisected to 8 bits) give thresholds T = lower - 16 x slack, retaining every
B-smooth cell with omitted mass below slack; survivors' cofactors are at most
303 / 232 bits. Survivors get exact norms, trial division below 2^20 and the
`ecm16.py` schedule at 99.9% per prime ([12]: 403 curves, B1 = 2^16), charged
2^39.0892 instructions per survivor; the same 99.9% enters the yield's miss model.

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre; the remainder over the block disk plus a 2^-32 relative term is
added and results are rounded outward. `nbound31` replays the Python leaves
bit for bit at 1,508.2 Ir per evaluation (charged 1,508.2); [2] checks the
bounds against exact norms. Every block of every row is enumerated, sieved and
charged; no block is skipped.

**Yield ([4]).** rels = (va_region/32000) vb (0.972556868) sum over
canonical primitive region rows of w rho(u_f) rho(u_g) r_f r_g x ECM miss,
with Dickman rho ([1]) and retentions r from the Euler-product threshold model
(`threshold_model.py`). Retained relations >= 2^64.3065 against 1.3 U =
2^63.8299: margin >= 0.4766 bits (pool SE 0.083), above the 0.30-bit floor
(`MARGIN`). The descent input `LOG2_RELS = 64.306` is below this count and the
mean-norm inputs 744.2/499.2 are at or below the relation-weighted means; [4]
checks both.

**Cost ([5]).** Charged mean 408.4 instructions/cell over 2^126.863 cells:
all-cell bounds, layout, setups, line path, record walk (primes in
[2^63, 2^64) charged twice), routing (10.01/record, 8/cell; STORAGE.md), apply
7.00, recovery 7.003/record, scan, zero, cofactorisation. Collection:
**<= 2^121.2492 units**, including the allowance below.

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo
estimate of a fixed domain total. Its relative SE has a conditional row-mean
term (0.730%) and an incidence term for the 3,907 contributing draws (1.499%),
total <= 1.668%. The charge multiplies the estimate by **1.0334 >= 1 + 2.00 x
0.01668**; [4a] fails if the allowance falls below 1 + 2.00 SE. This is a
modeling margin, not a probability term.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge (coefficients <= 255), simulated on
the production occurrence profile (`filtsim.py`; scale 2^17 table embedded, [6]
reruns 2^12) as the signed sigma-quotient system (each occurrence takes the
sign of its conjugate orientation). Three seeds at 2^17 give N/U <= 0.0471 at
mean weight <= 101 (charged 101/row). N/U falls with scale (`filtsim.run`,
kmax 64, weight 100; max over seeds): 0.0850 (2^12, 10 seeds), 0.0635 (2^13,
6), 0.0563 (2^14, 6), 0.0516 (2^15, 3), 0.0469 (2^16, 3), 0.0461 (2^17, 3),
0.0425 (2^18, 2) (`SCALE_NU`, [10] checks non-increasing). The model charges N
= 1.096 x 0.0471 x U < 2^59.176: ratio <= 0.05163 exceeds every measurement from
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
1,640 (fixed reduction tree over the ten active 2^56-column owners, N < 10 x 2^56: nine sends and nine combinations, counted 1,638 in `row_tail.c`/[10]; valid for any row, sums < 2^323). Every step forms the m
x n projections x^T A^i y: (m + n) N^2 dense terms at 97 plus one reduction
(C_RED9 = 420, counted 385) per entry. Horner reconstruction (w <- A w + Y F_i
over deg F <= N/m + N/n coefficients) costs one product plus, per element, n
dense terms and one reduction. With c = 101 x 20.7 + 1,640 + 21 x 97 per row:
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n + 420) + 97 (m + n)) + 420 (m + n) N
instructions: **<= 2^117.5157 units** at the unchanged N < 2^59.176,
retaining (m, n) = (19, 4). The (m + n) x (m + n) generator of the length L = (N +
21)(1/m + 1/n) + 2 sequence (divide-and-conquer approximant bases, Karatsuba):
<= 4 (m + n)^3 L^log2(3) log2(L) field mults at 97 plus reduction 420:
**<= 2^107.1997 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** Region and polynomials are coin/challenge-independent, so the
  count is one fixed number; H1 puts it at >= 2^64.3065, >= 0.4766 bits above
  1.3 U (pool SE 0.083); the descent uses 64.306, below that count.
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
survivor-test retention is >= 0.8592 at every sieve node ([11b], full run); cofactors use
the `ecm16.py` schedule ([12]). Relation counts are negative binomial.
Log-less ideals are discarded or re-descended in a fresh lattice. The plan
(62 nodes + boot; min mu 34.6, max psi 0.587) gives expected two-log work
**<= 2^115.6681 units** ([7]) at the inputs of section 3. The meter caps both
logs at 2^116.3915 units (cap reached = failure), so P(cap) <=
2^(115.6681-116.3915) <= 0.605669 (Markov).

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 6.7554e-05 - 0.605669 - 2^-10.01 - 2^-10.70 - 2/r)
>= **0.392624** >= 0.39: main phase, two boots (phi <= 6.7554e-05 each), the
instruction-meter Markov bound, a special-q node counter and special-q repeats.

- Node counter: both logs abort as a failure past N_MAX = 2^27.3 special-q
  nodes. E[N] <= 2^17.29 (`counts16.py`, same recurrence with unit node cost),
  so P(abort) <= 2^(17.29-27.3) = 2^-10.01 (Markov).
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for
  1/p weighting share a special-q with probability <= 2^-10.70, charged.

The boot term treats the 18 relations' subtrees as independent (H5). Targets
are uniform nonidentity by power randomization (2/r charged). The final
pairing check e(P, Q) = z rejects any incorrect P.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (allowance 1.0334 included) | <= 121.2492 |
| linear algebra, ten active vector owners | <= 117.5157 |
| matrix generator | <= 107.1997 |
| both individual logarithms (metered cap) | 116.3915 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| field map, pairings, verification and output | 66 |

The rows sum to 2^121.39941 < **2^121.400000**. Nothing is amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; relation collection alpha (5, 0.85) against the in-package computed (4.868096, 0.843643) over p <= 2000 and 0.731544 over p < 50000 for g ([14]); descent alpha (5, 1); exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool sampled on the theta = 2/5 region, exact flags, certified
counts with their exact rational products, and the canonical work mean
charged with a 3.34% modeling-margin allowance >= 2.00
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
CADO-NFS. Sub-region estimator, certified counts, in-package alpha, unit-entry
LA kernel, 99.9% ECM, routing and plan refinements: kaankolcu, pepedesigner,
mpjunior92 and muse-spark-1.3-contributor lineages (3a00ffc8). Ten active
owners: bernardo-lima12 (6282ab79). Package layout follows our sister-track
v12p. This revision (directly sampled 32,000-draw theta = 2/5 pool, allowance,
descent inputs and meter, recover-P reduction): Claude Opus 5.5 in Claude Code.

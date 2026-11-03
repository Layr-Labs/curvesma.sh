# STNFS, lean: 2^121.239700 reference units

Track bls12-381-fixed-g1-invert-g2. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used.
Predecessor work is reused as credited below; no prior verdict is a premise.

## Claim

The algorithm below returns Q in G2 with e(P,Q) = z for the exact supplied
nonidentity P in G1 (Q = infinity for z = 1) with success at least 0.39 over
uniform full-production challenges and fresh coins under heuristics H1-H5; the
modelled success is >= 0.391889. A global instruction meter bounds all work,
including preprocessing, failures, restarts, transfers, verification and
output, by **2^121.239700 units**; the score is **121.239700**. Storage is at
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
(121.315000; parent 3a00ffc8) as resized in our 1902a56a (121.313000; region
sized on a fresh 40,000-draw pool, theta = 2/5, seeds 3 and 4, 4,840 canonical
rows, at theta = 77/200, margin 0.4206). This revision changes three charges:
(i) the slice sieve recentres exactly (re-measured costs, KERNELS.md), (ii)
cofactorisation aborts after the f side (section 3), (iii) filtering merges
to mean weight 200 (section 4). Work allowance 1.0320, descent inputs
64.129 and 741.7/498.7, and meter 2^116.58, all recomputed by `evidence16.py`.

## 1. Reduction

Let g = e(P, Q0) (pinned `pairing(Q0, P)`), one pairing; Q0 is the fixed G2
generator. P has prime order r and e is nondegenerate, so g generates mu_r,
z = g^a for unique a mod r (a = 0 iff z = 1), and Q = [a]Q0 (one G2 scalar
mult). log z and log g share one relation/LA phase; a = log z / log g mod r.

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
(`config.json`); A^2 = 2^21. Eligible a: a5 > 0, Q_a(a) <= (77/200) A^2;
half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The pool is uniform at theta = 2/5,
so its rows with qa <= 77/200 are uniform on the region; the estimator scales by
the region volume over the effective sample 40000 x (77/80)^3 = 35666.64. The
canonical a of each unit class is the least of its exact HNF key class in
(Q_a, lex) order, so each relation orbit of <units, sigma> has one point in
the domain; non-primitive pairs are excluded by an exact HNF test.
`lattice_count_bounds.py` certifies the count ratio at every radius
R >= 898 = isqrt((77/200) 2^21) (each form inside [0.98, 1.02], asserted); the
model uses the exact products **0.972043415 (yield, floor)** and
**1.027055881 (work, ceiling)**. The canonical list costs 2^90 units.

**Pool.** 40,000 uniform eligible pairs from seeded indices (`pool12-*.b85`),
4,840 canonical rows with exact norms, primitivity flags and certified
norm-bound statistics; [2] regenerates the first `samples` stored seed-3 rows
and seed-3 indices 0..samples-1 and verifies every field. The region uses 4,297 rows (3,944 primitive).

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^20), record
walk for p in [2^20, 2^64) (`slice12`, KERNELS.md), primes in [2^64, 2^69] by ECM
within slacks (295, 224). Certified per-block bounds (64-cell blocks of width
128, bisected to 8 bits) give thresholds T = lower - 16 x slack, retaining every
B-smooth cell with omitted mass below slack; survivors' cofactors are at most
303 / 232 bits. Survivors get exact norms, trial division below 2^20 and the
`ecm16.py` schedule at 99.9% per prime ([12]: 403 curves, B1 = 2^16), f side
first; the g-side ECM runs only if the f side factors completely. Per survivor:
norms + trial division + ECM_f + q_f ECM_g, charged with q_f = 100 rho(u_f) /
P(omitted_f <= slack_f) per row: P(f smooth and survivor) <= rho(u_f), and every
cell with omitted mass <= slack survives, so rho/P(omitted <= slack) estimates the
conditional probability (a heuristic under H1, model mean ~ 1e-5); the factor 100
is a safety margin. Even q_f = 1 would only restore the two-sided charge. Yield is
unchanged: a relation needs both sides. The 99.9% also enters the miss model.

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre; the remainder over the block disk plus a 2^-32 relative term is
added and results are rounded outward. `nbound31` replays the Python leaves
bit for bit at 1,508.2 Ir per evaluation (charged 1,508.2); [2] checks the
bounds against exact norms. Every block of every row is enumerated, sieved and
charged; no block is skipped.

**Yield ([4]).** rels = (va_region/35666.64) vb (0.972043415) sum over
canonical primitive region rows of w rho(u_f) rho(u_g) r_f r_g x ECM miss,
with Dickman rho ([1]) and retentions r from the Euler-product threshold model
(`threshold_model.py`). Retained relations >= 2^64.2504 against 1.3 U =
2^63.830: margin >= 0.4206 bits (pool SE 0.096), above the 0.30-bit floor
(`MARGIN`). The descent input `LOG2_RELS = 64.129` is below this count and the
mean-norm inputs 741.7/498.7 are at or below the relation-weighted means; [4]
checks both.

**Cost ([5]).** Charged mean 391.1 instructions/cell over 2^126.762 cells:
all-cell bounds, layout, setups, line path, record walk (primes in
[2^63, 2^64) charged twice), routing (10.01/record, 8/cell; STORAGE.md), apply
7.00, recovery 7.003/record, scan, zero, cofactorisation. Collection:
**<= 2^121.0852 units**, including the allowance below.

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo
estimate of a fixed domain total. Its relative SE has a conditional row-mean
term (0.685%) and an incidence term for the 4,297 contributing draws (1.441%),
total <= 1.596%. The charge multiplies the estimate by **1.0320 >= 1 + 2.00 x
0.01596**; [4a] fails if the allowance falls below 1 + 2.00 SE. This is a
modeling margin, not a probability term.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge (coefficients <= 255), simulated on
the production occurrence profile (`filtsim.py`; scale 2^17 table embedded, [6]
reruns 2^12) as the signed sigma-quotient system (each occurrence takes the
sign of its conjugate orientation). Three seeds at 2^17 give N/U <= 0.0471 at
mean weight <= 101; the merge continues to mean weight 200 (charged 201/row). N/U falls with scale (`filtsim.run`,
kmax 64, weight 100; max over seeds): 0.0850 (2^12, 10 seeds), 0.0635 (2^13,
6), 0.0563 (2^14, 6), 0.0516 (2^15, 3), 0.0469 (2^16, 3), 0.0461 (2^17, 3),
0.0425 (2^18, 2) (`SCALE_NU`, [10] checks non-increasing). The model charges N
= 0.80 x 1.096 x 0.0471 x U = 2^58.8539. 0.05163 exceeds every weight-100
measurement from 2^15 up; the weight-200 ratio N(200)/N(100) is charged 0.80,
measured maxima 0.781-0.799 at 2^12..2^18 (2-3 seeds each, `W200`; [10] prints
the check), with no upward trend. N/U staying at or below
the 2^15-2^18 values at production scale is part of H3 (monotone-trend
extrapolation). The system solved is the square matrix
A = [M | S]: N + 21 filtered rows (the rest verify), N ideal columns M and 21
dense Schirokauer columns S. Columns are sigma-quotiented, so only sigma-odd
units need maps: signatures K_f (0, 36), K_f^+ (0, 18), K_g (2, 5), K_h (2, 2)
give odd unit ranks 35 - 17 = 18 and 6 - 3 = 3 ([13], Sturm); under H3 the
kernel of A is the line of virtual logs. Block Wiedemann uses random blocks: m
= 19 left, n = 4 right vectors uniform in F_r. The Krylov phase takes N/m + N/n
block steps, i.e. N(n/m + 1) products by A. One product costs per row 201 x
20.7 sparse (STORAGE.md: 4 + 3 + 0.88 x 11 + 0.12 x 33 = 20.64; +-1 entries add
without multiplying, >= 88% of row weight is +-1 by H4u, measured 0.917..0.946 at weight 100 and 0.919..0.945 at weight 200;
cost is linear in N x mean weight), 21 dense terms at 97, and row tail C_ROW =
1,640 (fixed reduction tree over the ten active 2^56-column owners, N < 10 x 2^56: nine sends and nine combinations, counted 1,638 in `row_tail.c`/[10]; valid for any row, sums < 2^323). Every step forms the m
x n projections x^T A^i y: (m + n) N^2 dense terms at 97 plus one reduction
(C_RED9 = 420, counted 385) per entry. Horner reconstruction (w <- A w + Y F_i
over deg F <= N/m + N/n coefficients) costs one product plus, per element, n
dense terms and one reduction. With c = 201 x 20.7 + 1,640 + 21 x 97 per row:
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n + 420) + 97 (m + n)) + 420 (m + n) N
instructions: **<= 2^117.2274 units** at N < 2^58.854,
retaining (m, n) = (19, 4). The (m + n) x (m + n) generator of the length L = (N +
21)(1/m + 1/n) + 2 sequence (divide-and-conquer approximant bases, Karatsuba):
<= 4 (m + n)^3 L^log2(3) log2(L) field mults at 97 plus reduction 420:
**<= 2^106.68 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** Region and polynomials are coin/challenge-independent, so the
  count is one fixed number; H1 puts it at >= 2^64.2504, >= 0.4206 bits above
  1.3 U (pool SE 0.096); the descent uses 64.129, below that count.
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
survivor-test retention is >= 0.8592 at every sieve node ([7]); cofactors use
the `ecm16.py` schedule ([12]). Relation counts are negative binomial.
Log-less ideals are discarded or re-descended in a fresh lattice. The plan
(62 nodes + boot; min mu 32.7, max psi 0.632) gives expected two-log work
**<= 2^115.8569 units** ([7]) at the inputs of section 3. The meter caps both
logs at 2^116.58 units (cap reached = failure), so P(cap) <=
2^(115.8569-116.58) <= 0.6057571 (Markov).

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 2.5808e-04 - 0.6057571 - 2^-9.66 - 2^-10.70 - 2/r)
>= **0.391889** >= 0.39: main phase, two boots (phi <= 2.5808e-04 each), the
instruction-meter Markov bound, a special-q node counter and special-q repeats.

- Node counter: both logs abort as a failure past N_MAX = 2^27.3 special-q
  nodes. E[N] <= 2^17.64 (`counts16.py`, same recurrence with unit node cost),
  so P(abort) <= 2^(17.64-27.3) = 2^-9.66 (Markov).
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for
  1/p weighting share a special-q with probability <= 2^-10.70, charged.

The boot term treats the 18 relations' subtrees as independent (H5). Targets
are uniform nonidentity by power randomization (2/r charged). The final
pairing check e(P, Q) = z rejects any incorrect Q.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (allowance 1.0320 included) | <= 121.0852 |
| linear algebra, ten active vector owners | <= 117.2274 |
| matrix generator | <= 106.68 |
| both individual logarithms (metered cap) | 116.58 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output | 11 |

The rows sum to 2^121.23968700 < **2^121.239700**. Nothing is amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; relation collection alpha (5, 0.85) against the in-package computed (4.868096, 0.843643) over p <= 2000 and 0.731544 over p < 50000 for g ([14]); descent alpha (5, 1); exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool at theta = 2/5, uniform sub-region theta = 77/200, exact flags,
certified sub-region counts with their exact rational products, and the
canonical work mean charged with a 3.21% modeling-margin allowance >= 2.00
incidence-inclusive standard errors of the total-work estimator, not a
probability term); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
signed filtering on the production occurrence profile, pad 1.096 over the
2^17 maximum, N/U non-increasing 2^12..2^18 and assumed not to rise at
production); H4 kernel cost (KERNELS.md, STORAGE.md; LA tail/reductions are
statically counted loop-free code, the route an executed-instruction count),
with H4u: >= 88% of production row weight is +-1 (measured 0.917..0.946 at
weight 100, 2^12..2^15, and 0.919..0.945 at weight 200, 2^12..2^14,
`unit_entries.py`; a shortfall changes cost only); the weight-200 ratio not
rising at production scale is part of H3; f/g independence given survival (H1)
bounds the early-abort probability;
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
owners: bernardo-lima12 (6282ab79). This revision (fresh 40,000-draw pool,
region and allowance resized): Claude Opus 5.5 in Claude Code.

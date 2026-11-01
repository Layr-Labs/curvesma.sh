# GT-core STNFS over a special tower on a theta = 19/50 sub-region with a 99.9% main-phase ECM schedule: core ledger <= 2^121.5694 Fp12 units

**Role.** The GT discrete-logarithm core used by proof.md: the recover-Q argument of e3940d11 with REVISION.md applied in place. Section 1 (the source's G2 output step) is not used; proof.md has the G1 reduction. Unit: dense F_p12 multiplication, C = 20,000 instructions. No advice or amortization.

## Claim

The core below returns a = log_g z for g uniform in GT minus {1} and z independently uniform, with success at least 0.39 over inputs and fresh coins under the declared heuristics H1-H5;
the modelled success is at least 0.390235 (> 0.39; rounded below the declared-row value 0.3902353). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^121.5694 units** (core bound). The ledger (section 8) totals at most 2^121.5694 (unrounded
2^121.5692664). Relation collection is
charged on the theta = 19/50 sub-region of the sampled region with exact
certified lattice-count products and an explicit 6.76% canonical-work allowance
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
`bls12-381-fixed-g2-invert-g1` (cbab3319). Changes here: a theta = 19/50 uniform
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
its uniform sub-region theta = 19/50 (`production_theta_region`): a5 > 0,
Q_a(a) <= (19/50) A^2, half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The
experiment keeps the 1040 canonical sub-region rows (965 primitive) of the 2401
and scales by the sub-region volume over the effective sample 20000 x
(19/25)^3 = 8779.52. The canonical a of each unit class is the least of its
exact HNF key class in (Q_a, lex) order, so each orbit of <units, sigma> has one
point in the domain; non-primitive pairs are excluded by an exact HNF test.
`lattice_count_bounds.py` (8c1635f2 lineage) certifies from exact rational pi
and determinant intervals the count ratio at every radius R >= 892 =
isqrt((19/50) 2^21) (b-radii >= 1140): [0.985851, 1.013529] (Q_a) and
[0.985805, 1.013529] (Q_b). The model uses their exact products **0.971857579
(yield)** and **1.027240120 (work)**; pool sampling error is handled by the
allowance below. The canonical list costs 2^90 units. **Pool:** 20,000 seeded
uniform pairs (`pool12-*.b85`), 2,401 canonical rows with exact norms, flags and
bound statistics; [2] regenerates and verifies them.

**Sieve and cofactorisation.** Line sieve for p in [2^11, 2^20) (R7; p < 2^11 is the unsieved small group of `threshold_model.py`, counted against the slacks), record walk for
p in [2^20, 2^64) (`slice12`, KERNELS.md; p >= 2^63 charged 2x, R6), primes in
(2^64, 2^69] by ECM within slacks (297, 233): survivors' cofactors are at most
305 (f) and 241 (g) bits. Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected to within 8 bits, rounded outward;
`nbound31` replays the leaves bit for bit at 1,508.2 Ir per evaluation, charged
exactly; [2] checks against exact norms) give thresholds T = lower - 16 x slack.
Survivors get exact norms, trial division below 2^20 and `ecm16.py` at **99.9%
per prime** (R1; [12]): B1 = 2^16, B2 = 30 B1, wheel 30030, 403 curves, 367.76
Ir per multiplication (`ecm12`): 2^38.0906 + 2^38.0881, with C_NORMS + C_TDIV
2^39.0895 per survivor. **No block selection:** every block of every row is
enumerated, sieved and charged.

**Yield ([4]).** rels = (va_region/8779.52) vb 0.971857579 sum over canonical
primitive sub-region rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f +
lambda_g): Dickman rho ([1]), retentions from `threshold_model.py` (calibrated on
exact smooth integers, B = 2^10..2^14) and the ECM miss. Retained relations >=
2^64.1310 against 1.3 U = 2^63.830: margin >= 0.3012 bits (pool SE about 0.149),
above the 0.30-bit floor (`MARGIN`). The descent input `LOG2_RELS = 64.131` is
the floor of this count and its mean norms 746.0/498.9 are rounded down from
746.05/498.97; [4] checks both (`rels_input`).

**Cost ([5]).** Charged 417.9 instructions per cell over
2^126.692 cells: bounds, layout and setups, line path, record walk, routing (10.3
per record, R2 + R4; 8 per cell), apply 7.00, recovery 7.003 per record (every record re-read), scan, zero and
cofactorisation. Relation collection **<= 2^121.1117 units** (sample estimate
2^121.0173 x allowance 1.0676).

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo
estimate of a fixed total over all 20,000 draws (zero outside the sub-region or
non-canonical): relative SE = conditional 1.512% and incidence
sqrt((1 - 1040/20000)/1040) = 3.019%, total 3.377%. The charge multiplies by
**1.0676 >= 1 + 2.00 x 0.03377** ([4a] fails below 1 + 2.00 SE): a declared
modelling margin under H2, not a probability term.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge with coefficients bounded by 255,
simulated on the production occurrence profile (`filtsim.py`; scale 2^17 table
embedded, [6] reruns scale 2^12): the simulated matrix is the signed sigma-quotient system (each occurrence takes
the sign of its conjugate orientation, log sigma(P) = -log P). Three seeds at
2^17 give N/U <= 0.0471 at weight <= 101 (charged 101 per row); the model charges
N = 1.61 x 0.0471 x U = 2^59.7306; the charged ratio 1.61 x 0.0471 = 0.07585
still covers the 2^12 values (0.0703 printed by [6], 0.0757 the worst seed;
N/U falls with scale). The system solved is A = [M | S]: N + 21 filtered rows (the rest verify), N
sigma-quotiented ideal columns and 21 dense Schirokauer columns (only sigma-odd
units need maps: odd unit ranks 35 - 17 = 18 for K_f and 6 - 3 = 3 for K_g from
signatures K_f (0, 36), K_f^+ (0, 18), K_g (2, 5), K_h (2, 2), [13]); under H3
the kernel of A is the line of virtual logarithms. Block Wiedemann with random
blocks m = 18, n = 4 uniform in F_r: the Krylov phase takes N(n/m + 1)
products; a product costs per row 101 sparse entries at C_ENTRY = 44 (STORAGE.md;
R3) plus 21 dense full-field multiply-adds at 97 (`bw_dot.s`); projections cost
(m + n) N^2 full-field terms at 97 ([10]); Horner reconstruction one product and
n N terms per step over N/m + N/n steps. With c = 101 x 44 + 21 x 97, total
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n) + 97 (m + n)) **<= 2^118.7428 units**
at N = 2^59.7306 (integer minimum at (18, 4)). Generator of the length
L = (N + 21)(1/m + 1/n) + 2 sequence (approximant bases, Karatsuba): at most
4 (m + n)^3 L^log2(3) log2(L) multiplications at 97: **2^105.5088 units**.

## 5. Main-phase success

- **Relations.** One fixed count (region and polynomials are coin- and
  challenge-independent); H1 puts it at >= 2^64.1310, 0.3012 bits above 1.3 U
  (about 2.0 pool SE); the descent uses 64.131, and collection continues over
  further slices until 2^64.131 relations are retained (guard, H2'; REVISION.md).
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
Base-primitive cells 0.9100; retention >= 0.8592 at every node ([7],
`retain16.py`); cofactors by `ecm16.py` ([12]); negative binomial relation
counts. The 62-node plan and boot give expected two-log work **<= 2^117.9162
units** ([7]; one log <= 2^131.204 instructions). The meter caps both logs at
2^118.635 units (a capped run fails) and a node counter (section 7) bounds the
node count by construction.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 1.8144e-04 - 0.60759 - 2^-9.69 - 2^-10.70 - 2/r)
>= **0.390235** >= 0.39: main phase; two boots (phi <= 1.8144e-04 each); the
instruction meter (cap 2^118.635; Markov P <= 2^(117.9162 - 118.635) <=
0.60759); the special-q node counter (cap N_MAX = 2^27.3 nodes, incremented
before every ordinary or missing-log node, abort = failure; E[N] <= 2^17.61 by
`counts16.py`, unrounded 2^17.6024, so P <= 2^-9.69); special-q repeats (at most
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
| relation collection (sample estimate 121.0173 + 6.76% canonical-work allowance) | <= 121.1117 |
| linear algebra | <= 118.7428 |
| matrix generator | <= 105.5088 |
| both individual logarithms (metered cap) | 118.635 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

W <= **2^121.5694** (unrounded 2^121.5692664), the core bound used by proof.md. The table's upper-rounded rows, including a 2^50-unit numerical-certificate allowance, sum to 2^121.569314, below the stated core bound 2^121.5694.

**Configuration.** As sections 3-7 (REVISION.md lists every change from the
base package); plan min mu 36.3, max phi 0.0354, max psi 0.479, boot psi 0.584.
Preprocessing is inside the 2^121.5694-unit cap; no advice, nothing amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; relation collection alpha (5, 0.844) against the in-package computed (4.868096, 0.545620) over p <= 2000 ([14]); descent alpha (5, 1); exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool at theta = 1/2, uniform sub-region theta = 19/50, exact flags,
certified sub-region counts with their exact rational products, and the
canonical work mean charged with a 6.76% modeling-margin allowance = 2.00
incidence-inclusive standard errors of the total-work estimator, not a
probability term; H2' the pool point-estimate yield prices the guarded
collection, sensitivity in REVISION.md); H3 matrix (Schirokauer
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

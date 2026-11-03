# GT-core STNFS over a special tower on a theta = 719/2000 sub-region with a 99.9% main-phase ECM schedule: core ledger <= 2^121.28540 Fp12 units

**Role.** The GT discrete-logarithm core used by proof.md: the recover-Q argument of e3940d11 with REVISION.md applied in place. Section 1 (the source's G2 output step) is not used; proof.md has this track's recover-P reduction. Unit: dense F_p12 multiplication, C = 20,000 instructions. No advice or amortization.

## Claim

The core below returns a = log_g z for g uniform in GT minus {1} and z independently uniform, with success at least 0.39 over inputs and fresh coins under the declared heuristics H1-H5;
the modelled success is at least 0.3917 (> 0.39; model 0.3919129). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, verification and output, by **2^121.28540 units** (core bound). The ledger (section 8) totals at most 2^121.28540 (unrounded
2^121.28539453). Relation collection is
charged on the theta = 719/2000 sub-region of the sampled region with exact
certified lattice-count products and an explicit 7.72% canonical-work allowance
(2.00 incidence-inclusive standard errors of the sampled work estimator, [4a]). Distributed storage is at most
2^105 bytes (STORAGE.md) in individually bounded, explicitly addressed objects.
The claim is analytic and exploratory under the NFS heuristics H1-H5 (section
9) with the measurements and checks cited; no full-size run is made.

**Numerical conventions.** Decimals are the experiment's values rounded to the
digits shown; a figure used as a bound is written with <= or >= and rounded in the
safe direction, and the `declared_bounds` check compares every stated bound with
the unrounded value.

The model is `model16.py` and `descent16.py` with `config.json` and `stats16.b85`;
`evidence16.py` recomputes the ledger and the checks cited as [n].

## Base

The GT discrete-logarithm machinery (polynomials, factor base, sieve, filtering,
block Wiedemann, descent, kernels, pools) is the lean STNFS package of the
sister track `bls12-381-fixed-g2-invert-g1` (cbab3319); both pairing
orientations need the same two logarithms in GT. This package uses this track's
reduction (section 1), a theta = 719/2000 uniform sub-region (this revision; was 19/50) of the shipped pool
with exact certified count products and an explicit work allowance (section 3),
relation-collection Murphy alpha (5, 0.844) checked in-package ([14]), an
explicit special-q node counter (section 7) and a descent plan re-optimised for
the sub-region inputs (section 8).

## 1. Reduction

The source's recover-Q reduction is not used here; the recover-P reduction (base g = e(P0, Q) for the supplied Q, output [a]P0) is in proof.md. The core computes a = log z / log g mod r for a uniform nonidentity base g and an independently uniform z, with log z and log g against one shared relation and linear-algebra phase.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) and g = x^2 + |u| y, u = -0xd201000000010000.
[0b] checks (u y)^((p^6-1)/2) = -1 in F_p6, so F_p[y,x]/(h, g) = F_p12; [0]
checks Res_y(h, f(x^2 -> u y)) = 0 mod p. x -> -x is the conjugation sigma
(log sigma(P) = -log P); candidates and unknowns are quotiented by sigma.
Murphy alpha is charged at (5, 0.844): [14] recomputes it in-package over
p <= 2000 (`python3 evidence16.py --alpha` prints (4.868096, 0.545620));
alpha_g = 0.844 is above the computed value and its p <= 5000 value (about
0.46). The descent keeps the more conservative (5, 1).

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B on both
sides; U = 2^63.451 sigma-orbits (indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21. The pool is uniform at theta = 1/2; the collected
region is its uniform sub-region theta = 719/2000: eligible a: a5 > 0,
Q_a(a) <= (719/2000) A^2; half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The
experiment keeps the 849 canonical sub-region rows (789 primitive) of 2401
and scales by the sub-region volume over the effective sample 20000 x (719/1000)^3 =
7433.90. Each relation orbit of <units, sigma> has exactly one point in the
domain (least HNF key in (Q_a, lex) order); non-primitive pairs are excluded by
an exact HNF test. `lattice_count_bounds.py` certifies from exact rational
intervals the count ratio at every radius R >= 868 = isqrt((719/2000) 2^21):
[0.985462, 1.013905] for Q_a and [0.985415, 1.013905] for Q_b, each inside
[0.98, 1.02] (asserted). The model uses the exact products **0.971088854
(yield)** and **1.028002875 (work)**, not rounded-out factors. Sampling error is
handled by the canonical-work allowance below. The canonical list costs 2^90 units.

**Pool.** 20,000 uniform eligible pairs from seeded indices (`pool12-*.b85`),
2,401 canonical rows with exact norms, primitivity flags and certified
norm-bound statistics; [2] regenerates sampled indices at the sampling region
theta = 1/2 and verifies every field. The collected sub-region uses 849 of the
2401 rows (789 primitive).

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^48), record
walk for p in [2^48, 2^64) (`slice12`, KERNELS.md; R15 moves the split from 2^20 to the charged-cost crossover, section 8), primes in [2^64, 2^69] by ECM
within the slacks (293, 224), so survivors' cofactors are bounded by 301 and 232
bits (slacks re-optimised for the sub-region, section 8). Certified per-block
bounds (`embedding_bounds.py`, 64-cell blocks of width 128, bisected until both
bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every
B-smooth cell with omitted mass below slack. Survivors get exact norms, trial
division below 2^20, the band-strip for (2^20, 2^48] (below), record recovery
for [2^48, 2^64), and the package's explicit `ecm16.py` schedule at **99.9% per
prime**, the schedule the descent uses and
[12] checks: B1 = 2^15, B2 = 30 B1, wheel 30030, 764 curves, five-limb Montgomery
multiplication at 367.76 Ir (`ecm12`), 2^38.0758 and 2^38.0743 instructions for
the two sides. C_surv = C_NORMS + C_TDIV + band-strip + both sides = 2^39.6026 per survivor; 0.999
is the same ECM_FIND as in the yield's miss factor 0.999^(lambda_f + lambda_g).

**Band-strip (R16; R18 re-tune).** The R15 line path credits primes in (2^20, 2^48] in sieve
scores without recording hits, so each survivor recovers them by ECM: per side,
two fixed rounds of the `ecm16.py` schedule for primes <= 2^48 in the unstripped
cofactor at 99.9% per prime per round (B1 = 2^13, B2 = 30 B1, wheel 30030, 92
curves; 2^36.3298/2^35.2703 instructions per round, 2^37.8952 both sides both
rounds), dividing out found factors to exhaustion after each round. Unstripped
cofactors are capped at 832/576 bits (pool-region norm maxima 796.375/513.0,
rounded up to 13/9 limbs; one more limb moves the ledger under 0.005 bits, and
within-row norm variation is H2 estimator risk inside the allowance). At most
41/28 band factors per side (each above 2^20), so per-side strip failure is at
most 41 x (1e-3)^2 / 28 x (1e-3)^2 by the union bound; a failed strip loses the
relation (yield factor 0.99993100, about 1.0e-4 bits). After the strip the
301/232-bit cofactor bounds hold for the main ECM as before.

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre in double precision; the remainder over the block disk plus a
2^-32 relative term (coefficient intervals and rounding) is added, and results
are rounded outward. `nbound31` replays the Python leaves bit for bit on all
2,401 pool rows at 1,508.2 Ir per evaluation (charged 1,508.2, the measured marginal
per-block cost excluding the separately charged per-row and per-slice setups); [2] checks the bounds
against exact norms.

**No block selection.** Every block of every row is enumerated, sieved and
charged at the full per-cell cost; no block is skipped and no yield is dropped.

**Yield ([4]).** rels = (va_region/7433.90) vb (0.971088854) sum over canonical primitive sub-region
rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
and the ECM miss. Retained relations >= 2^64.0155 against 1.2 U = 2^63.714: margin
>= 0.3011 bits (unrounded 0.30119; pool SE about 0.161 bits), above the model's 0.30-bit feasibility floor
(`MARGIN` in model16.py). The descent's missing-log input is `LOG2_RELS = 63.925`
in `descent16.py`: this model's 2^64.01555 retained less the one-SE allowance
(DESCENT.md); the certified lower count product is already applied inside that
count. [4] checks `LOG2_RELS <= log2_rels` explicitly (the `rels_input` check).
The descent's relation-weighted mean-norm inputs are 744.8/498.9 bits, rounded
down from the sub-region's computed 744.84/498.97 and checked in [4].

**Cost ([5]).** The sampled mean is 313.4 instructions per cell; the charged
mean, with the canonical-work allowance below, is 337.6 instructions per cell
over 2^126.506 cells: all-cell bounds, layout and setups, the line path, record
walk, routing (10.005 per record: direct chunk scatter, R2, measured end to end, R3, repriced in R7; 8 per cell; STORAGE.md), apply 7.00, recovery
7.003, scan, zero and cofactorisation. Relation collection: **<= 2^120.61756 units**
(sample estimate 2^120.51027 times the allowance 1.0772; line/record/nbound prices are the exact raw-count quotients, R5).

**Canonical-work allowance (H2, [4a]).** The charged work is (va_region/7433.90)
times the sum of per-row work over the in-region canonical rows, a Monte Carlo
estimate of a fixed, challenge-independent domain total: a sum over all 20,000
draws with zero contribution from draws outside the sub-region or not canonical.
Its relative standard error has a conditional row-mean term and an incidence
term for the 849 contributing draws among 20,000, computed from the per-row
contributions ([4a]): conditional 1.894%, incidence 3.358%, total 3.856%.
The charge multiplies the sample estimate by **1.0772 >= 1 + 2.00 x
0.03856 = 1.07712**, and [4a] fails below 1 + 2.00 SE. This is a declared
modeling margin on a fixed quantity under H2, not a per-challenge failure
probability; no sampling tail is charged against the success slack.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge with coefficients bounded by 255,
simulated on the production occurrence profile (`filtsim.py`; scale 2^19 table
embedded, [6] reruns scale 2^12): the simulated matrix is the signed sigma-quotient system (each occurrence takes
the sign of its conjugate orientation, log sigma(P) = -log P). Eight seeds at
2^19 give N/U <= 0.0556 at weight <= 100.02 (fresh R18 basis at excess 1.2,
replacing the R11 1.3 basis; 2^12-2^18 runs are trend evidence); the model charges
N = 1.25 x 0.0556 x U = 2^59.6045 (ratio 0.0695), covering the trend sims at scales
2^15-2^18 and all 11 new 2^18/2^19 runs (worst 0.0644; means fall 0.084 -> 0.054 over 2^12-2^19).
Mean weight 100.7/row covers all sims at scales >= 2^13 (worst 100.53). The system solved is the square matrix
A = [M | S]: N + 21 filtered rows (the rest verify the result), the N
ideal columns M and 21 dense Schirokauer columns S. As columns are sigma-quotiented
(log sigma(P) = -log P), only sigma-odd units need maps: signatures K_f (0, 36),
K_f^+ = K_h(x^2) (0, 18), K_g (2, 5), K_h (2, 2) give odd unit ranks 35 - 17 = 18
and 6 - 3 = 3 ([13], Sturm sequences); under H3 the kernel of A is the line of virtual
logarithms (with their Schirokauer coefficients). Block Wiedemann
on A uses random blocks: m = 22 left and n = 4 right vectors uniform in F_r (the
coins). The Krylov phase takes N/m + N/n
block steps, i.e. N(n/m + 1) products by A; a product costs per row 100.7 sparse
entries x 73 instructions (STORAGE.md; `la_mac.s` inner loop of 33, [10]) plus 21
dense full-field multiply-adds at 97 (`bw_dot.s`, bw_full: 4 x 4-limb product into
a 9-limb accumulator, inner loop of 97). Each row's sums are reduced mod r once per product inside the 13-per-entry finishing charge (STORAGE.md, R4, repriced in R5). The 73 charge is the per-entry
component sum (33 + 4 + 3 + 20 + 13, STORAGE.md) and 97 the static bw_full loop. Every step forms the m x n projections
x^T A^i y, (m + n) N^2 full-field terms at 97 ([10]). The kernel vector
is reconstructed by Horner, w <- A w + Y F_i over the deg F <= N/m + N/n generator
coefficients (Y the n starting vectors): one product and n N full-field
multiply-adds per step. With c = 100.7 x 73 + 21 x 97 per row, the total is
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n) + 97 (m + n)) instructions:
**<= 2^118.931853 units** at N = 2^59.6045 (the minimum of this expression over
integer block sizes is at (m, n) = (22, 4), re-verified at the R18 N). The (m + n) x (m + n) matrix
generator of the length
L = (N + 21)(1/m + 1/n) + 2 sequence is computed by divide-and-conquer approximant
bases with Karatsuba products: at most 4 (m + n)^3 L^log2(3) log2(L) full-field
multiplications, each at 97 plus one reduction C_RED = 775.52 (R4): **2^109.119675 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** The region and polynomials do not depend on coins or the
  challenge, so the relation count is one fixed number; H1 puts it at
  at least 2^64.0155, at least 0.3011 bits above 1.2 U (pool standard error
  about 0.161 bits, so the floor is about 1.9 SE above the requirement); the descent uses 63.925,
  one SE below the central estimate.
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
lattices, then special-q nodes that sieve an exact lattice ball and recurse on
the child ideals. Counted cells are base-primitive (fraction 0.9100); retention
is at least 0.8593 at every sieve node ([7]; one local Hensel bound per root
below 2^9 plus the R13 union allowances); cofactors use the `ecm16.py` schedule
([12]). Relation counts are negative binomial (shapes from toy fits). Ideals
without a log are discarded or descended in a fresh lattice. The shipped plan
(62 nodes and the boot) gives expected two-log work **<= 2^118.03397 units**
([7]; one log at most 2^131.322 instructions); record layout in STORAGE.md. The
descent meter caps both logs at 2^118.76810 units; P(cap reached) <= E[W]/cap <=
2^(118.03397 - 118.76810) <= 0.60119 (Markov). A separate special-q node counter
(cap 2^27.3 nodes, section 7) bounds the number of nodes by construction.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 2.49e-03 - 0.60119 - 2^-9.50 - 2^-10.70 - 2/r)
>= **0.3917** >= 0.39: main phase, two boots
(phi <= 2.4860e-03 each at the sub-region relation input), the instruction-meter
Markov bound, the node-counter Markov bound and special-q repeats, as follows.

**Two explicit meters.** The descent for both logarithms runs under (i) an
instruction meter with cap 2^118.76810 units and (ii) a **special-q node counter**
with cap N_MAX = 2^27.3: before starting a special-q node (ordinary or a
missing-log node in the lattice of Q q') the counter is incremented, and the
run aborts, counted as a failure, if it would exceed N_MAX. Hence at most
2^27.3 special-q nodes are ever run, by construction, whatever the node costs.

- Instruction meter: P(cap reached) <= E[W]/cap <= 2^(118.03397 - 118.76810) <=
  0.60119 (Markov, E[W] the expected two-log work of section 6).
- Node counter: P(counter abort) <= E[N]/N_MAX (Markov). E[N] is bounded by the
  same finite-tree recurrence that gives E[W], with every node's own cost
  replaced by 1 (`counts16.py`). Twice the boot value gives
  E[N] <= 2^17.80 nodes for both logs (unrounded 2^17.7436), so
  P(counter abort) <= 2^(17.80 - 27.3) = 2^-9.50 <= 0.0013811.
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
pairing check e(P, Q) = z rejects any incorrect P. The success bound is
conditional on H1-H5, including the H2 canonical-work allowance of section 3,
which is a resource-model margin and not a probability term.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (sample estimate 120.51027 + 7.72% canonical-work allowance) | <= 120.61756 |
| linear algebra (refiltered, R18; R7 transfers) | <= 118.931853 |
| matrix generator (reduced, R4; 22x4, R7) | <= 109.119675 |
| both individual logarithms (metered cap, R18) | 118.76810 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output (pairing, G1 scalar multiplication, check) | 11 |

W <= **2^121.28540** (unrounded 2^121.28539453), the core bound used by proof.md. The table's upper-rounded rows, including a 2^50-unit numerical-certificate allowance, sum to 2^121.28539994, below the stated core bound 2^121.28540.

**Configuration.** Theta = 719/2000 (849 of 2401 rows), counts certified at
radius >= 868, allowance 1.0772 (2.00 SE), ECM 99.9%, slacks (293, 224),
descent LOG2_RELS = 63.925, means 744.8/498.9. The 62-node plan was found by
coordinate search (mu >= 16, phi <= 0.05, psi <= 0.8), entries in git history.
Min mu 32.1, max phi 0.0389, max psi 0.534, boot psi 0.687; two-log
expectation <= 2^118.03397, meter 2^118.76810, node counter of section 7. Other
charges are the base package's. Preprocessing per instance is bounded by
the core work cap and included in it. No external advice. Nothing amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; relation collection alpha (5, 0.844) against the in-package computed (4.868096, 0.545620) over p <= 2000 ([14]); descent alpha (5, 1); exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool at theta = 1/2, uniform sub-region theta = 719/2000, exact flags,
certified sub-region counts with their exact rational products, and the
canonical work mean charged with a 7.72% modeling-margin allowance = 2.00
incidence-inclusive standard errors of the total-work estimator, not a
probability term); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
filtering on the production occurrence profile, signed, pad 1.25); H4 kernel cost (KERNELS.md, STORAGE.md);
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing logs, independent
subtree outcomes, including ideals without a log re-descended in fresh lattices,
local omitted-mass distributions, ECM curve orders as random integers).

## 10. Memory and provenance

Storage <= 2^105 bytes (STORAGE.md). Core: Gautam Botrel (Claude Code/Codex);
lean refinements kaankolcu/pepedesigner (cbab3319); sub-region/counts/alpha
kaankolcu/DeepSeek (8c1635f2); port/counter/plan-checks Claude Opus 5.5. No
prior review verdict is used as a premise. Quadrature: REFINEMENT.md. The
two-log meter uses Markov, with no independence assumption between logarithms.

R5-R18 (muse-spark-1.3-contributor): exact prices, 10.005 route, 13/entry
finishing, 0.05% reserve; refiltered N/U (R6: 0.0662; R11: 0.0614; R18: 0.0695),
weight 100.7; transfers 20/entry, 22x4 blocks; descent kernels, rels allowance,
corrected loss, retention allowances, re-optimised plan, union allowances, sieve/region re-tune, plan/grid re-opt, excess/filter re-tune. Rest unchanged.

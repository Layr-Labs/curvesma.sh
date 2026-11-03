# GT-core STNFS over a special tower on a theta = 190/500 sub-region with a 99.9% main-phase ECM schedule: core ledger <= 2^121.058691 Fp12 units

**Role.** The GT discrete-logarithm core used by proof.md: the recover-Q argument of e3940d11 with REVISION.md applied in place (proof.md has the G1 reduction). Unit: dense F_p12 multiplication, C = 20,000 instructions. No advice or amortization.

## Claim

The core below returns a = log_g z for g uniform in GT minus {1} and z independently uniform, with success at least 0.39 over inputs and fresh coins under the declared heuristics H1-H5;
the modelled success is at least 0.50 (claim 0.39; rounded below the model value 0.5000103). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, verification and output, by **2^121.058691 units** (core bound). The ledger (section 8) totals at most 2^121.058691 (unrounded
2^121.0586517). Relation collection is
charged on the theta = 190/500 sub-region of the sampled region with exact
certified lattice-count products and an explicit 6.91% canonical-work allowance
(2.00 incidence-inclusive standard errors of the sampled work estimator, [4a]). Distributed storage is at most
2^105 bytes (STORAGE.md) in individually bounded, explicitly addressed objects.
Analytic exploratory claim under the NFS heuristics H1-H5 (section 9) with
the measurements and checks cited; no full-size run.

**Numerical conventions.** Decimals are rounded in the safe direction
(`declared_bounds` checks each stated bound against the unrounded value).

Model: `model16.py`, `descent16.py` with `config.json`, `stats16.b85`;
`evidence16.py` recomputes the ledger and the checks cited as [n].

## Base

The machinery (polynomials, factor base, sieve, filtering, block Wiedemann,
descent, kernels, pools) is the lean STNFS package of the sister track
`bls12-381-fixed-g2-invert-g1` (cbab3319), with a theta = 190/500 uniform
sub-region with certified count products and a work allowance (section 3),
alpha (5, 0.844) checked in-package ([14]), an explicit special-q node counter
(section 7) and the descent plan of section 8.

## 1. Reduction

proof.md has the G1 reduction. The core computes a = log z / log g mod r for g
uniform nonidentity and z independently uniform, both logs against one shared
relation and linear-algebra phase.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y),
g = x^2 + |u| y, u = -0xd201000000010000. [0b]: (u y)^((p^6-1)/2) = -1 in F_p6,
so g is irreducible over F_p6 and F_p[y,x]/(h, g) = F_p12; [0]: Res_y(h, f(x^2
-> u y)) = 0 mod p. x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) =
-log P); candidates and unknowns are quotiented by sigma. Murphy alpha is
charged at (5, 0.844) in collection, above the in-package (4.868096, 0.545620)
over p <= 2000 ([14], `evidence16.py --alpha`); the descent keeps (5, 1).

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69;
the factor base is the relative-degree-1 prime ideals of norm <= B on both
sides, U = 2^63.451 sigma-orbits (64-bit indices).

**Region.** Exact rational forms Q_a, Q_b (`config.json`), A^2 = 2^21. The pool
is uniform at theta = 1/2 (a5 > 0, Q_a(a) <= A^2/2); the collected region is
its uniform sub-region theta = 190/500: a5 > 0, Q_a(a) <= (190/500) A^2,
half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The experiment keeps the 1040
canonical sub-region rows (965 primitive) of the 2401 and scales by the
sub-region volume over the effective sample 20000 x (190/250)^3 = 8779.52. The
canonical a of each unit class is the least of its exact HNF key class in (Q_a,
lex) order (one point per <units, sigma> orbit); non-primitive pairs are
excluded by an exact HNF test.
`lattice_count_bounds.py` certifies from exact rational pi and determinant
intervals the count ratio at every radius R >= 892 = isqrt((190/500) 2^21)
(b-radii >= 1140): [0.985851, 1.013529] (Q_a), [0.985805, 1.013529] (Q_b); the
model uses their exact products **0.971857579 (yield)** and **1.027240120
(work)**; pool sampling error is handled by the allowance below. The canonical list costs 2^90 units. **Pool:** 20,000 seeded
uniform pairs (`pool12-*.b85`), 2,401 canonical rows with exact norms, flags and
bound statistics ([2] regenerates and verifies them).

**Sieve and cofactorisation.** Line sieve for p in [2^11, 2^45.5) (p < 2^11 is
the unsieved small group of `threshold_model.py`, counted against the slacks),
record walk for p in [2^45.5, 2^64) (`slice12`, KERNELS.md; p >= 2^63 charged
2x), primes in (2^64, 2^69] by ECM within slacks (297, 233): cofactors at most
305 (f) and 241 (g) bits. Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected to within 8 bits, rounded outward;
`nbound31` replays the leaves at 1,508.2 Ir per evaluation, charged exactly; [2]
checks against exact norms) give thresholds T = lower - 16 x slack. Survivors get exact norms, trial division
below 2^30 (R14), first-pass line-hit records over [2^30, 2^45.5) (R15) and
`ecm16.py` at **99.9% per prime** (R1; [12]): B1 = 2^16, B2 = 30 B1, wheel
30030, 403 curves, 367.76 Ir per multiplication (`ecm12`): 2^38.0907 +
2^38.0881, with C_NORMS + C_TDIV 2^39.2445 per survivor. Every block of every
row is enumerated, sieved and charged (no block selection).

**Yield ([4]).** rels = (va_region/8779.52) vb 0.971857579 sum over canonical
primitive sub-region rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f +
lambda_g): Dickman rho ([1]), retentions from `threshold_model.py` (calibrated on
exact smooth integers) and the ECM miss. Retained relations >= 2^64.1311
against 1.3 U = 2^63.830: margin >= 0.3012 bits (pool SE about 0.150), above
the 0.30-bit floor (`MARGIN`). The descent input `LOG2_RELS = 63.79` is
below the count times (1 - 2 relative pool SEs), 2^63.7959 (R18, was 1 SE; audit `yield_lower_bound_covers_descent`) and its mean norms 745.9/498.9 are rounded down from
746.05/498.97; [4] checks both (`rels_input`).

**Cost ([5]).** Charged 364.1 instructions per cell over 2^126.692 cells:
bounds, layout and setups, line path, record walk, routing (10.3 per record, 8
per cell), apply 7.00, recovery 7.003 per record (every record re-read), scan,
zero and cofactorisation. Relation collection **<= 2^120.9128 units** (sample
estimate 2^120.8165 x allowance 1.06905).

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo
estimate of a fixed total over all 20,000 draws (zero outside the sub-region or
non-canonical): relative SE = conditional 1.646% and incidence
sqrt((1 - 1040/20000)/1040) = 3.019%, total 3.439%; the charge multiplies by
**1.06905 >= 1 + 2.00 x 0.03439** ([4a]): a modelling margin under H2, not a
probability term.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge with coefficients bounded by 255,
simulated on the production occurrence profile (`filtsim.py`; 2^17 table
embedded, [6] reruns 2^12) as the signed sigma-quotient system. Three seeds at
2^17 give N/U <= 0.03225 at mean weight <= 304 (charged 304 per row; R17). N/U
falls with scale (`SCALE_NU`: maxima 0.0652, 0.0431, 0.0368, 0.0348, 0.0334,
0.0323, 0.0304 at 2^12..2^18; checked non-increasing by [10]; `audit_g1.py`
{"unit_entries": [k, seed]} measures H4u up to 2^15). The model charges N =
1.10 x 0.03225 x U = 2^58.6340: the ratio 0.03547 is above every measured value
from 2^15 up. That N/U at production scale stays at or below the 2^15-2^18
values is part of H3. The system solved is the square matrix A = [M | S]:
N + 21 filtered rows (the rest verify the result), the N ideal columns M and 21
dense Schirokauer columns S. Columns are sigma-quotiented (log sigma(P) = -log
P), so only sigma-odd units need maps: signatures K_f (0, 36), K_f^+ = K_h(x^2)
(0, 18), K_g (2, 5), K_h (2, 2) give odd unit ranks 18 and 3 ([13], Sturm
sequences); under H3 the kernel of A is the line of virtual logarithms. Block
Wiedemann uses random blocks m = 23, n = 4 (the coins; R18): N(n/m + 1) products
by A in the Krylov phase, each row costing 304 sparse entries at C_ENTRY = 21
(H4u), 21 dense terms at 97 (bw_full) and the row tail C_ROW = 1,350 (reduction
mod r; STORAGE.md, kernels counted by [10]); per step the m x n projections
((m + n) N^2 dense terms at 97 plus one C_RED9 = 420 reduction per entry); Horner
reconstruction w <- A w + Y F_i over deg F <= N/m + N/n (one product plus n
dense terms and one reduction per element). With c = 304 x 21 + 1,350 + 21 x 97
= 9,771 per row the total is N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n + 420) +
97 (m + n)) + 420 (m + n) N instructions: **<= 2^117.0498 units** at N =
2^58.6340 (minimum over integer block sizes at (23, 4)). The (m + n) x (m + n)
generator of the length L = (N + 21)(1/m + 1/n) + 2 sequence by
divide-and-conquer approximant bases: at most 4 (m + n)^3 L^log2(3) log2(L)
multiplications at 97 + 420: **<= 2^106.9498 units**.

## 5. Main-phase success

- **Relations.** One fixed count (region and polynomials are coin- and
  challenge-independent); H1 puts it at >= 2^64.1311, 0.3012 bits above 1.3 U;
  the descent uses 63.79 <= the count at 2 relative SEs below, and collection
  continues over further slices until 2^63.79 relations are retained (guard).
- **Linear algebra.** H3: r unramified in K_f and K_g ([13]), dividing neither
  class number, the 21 sigma-odd Schirokauer maps injective on units mod r-th
  powers, kernel of A the line of virtual logarithms. Block Wiedemann fails with
  probability <= 4N/r < 2^-192 (Schwartz-Zippel); a failed solve (checked on the
  dropped rows) counts as failure.
- **Reconstruction.** Purged ideals' logs by back-substitution; the descent
  treats every ideal with at most two occurrences as logless (section 6).

Main phase success >= 1 - 2^-192 over the coins (H1, H3).

## 6. Individual logarithms

DESCENT.md: an ECM boot on balls of certified p^6 lattices, then special-q
nodes sieving exact lattice balls (`dsieve16`) and recursing on child ideals;
retention >= 0.8575 at every node ([7], `retain16.py`); cofactors by `ecm16.py`
([12]); negative binomial relation counts; P integrated over measured norm
sub-cells (R18). The 62-node plan and boot give expected two-log work **<= 2^115.1678
units** ([7]; one log <= 2^128.456 instructions). The meter caps both logs at
2^116.1806 units (a capped run fails) and a node counter (section 7) bounds the
node count by construction.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 1.256e-03 - 0.49557 - 2^-9.58 - 2^-10.70 - 2/r)
>= **0.50** >= 0.39 (claim): main phase; two boots (phi <= 1.256e-03 each); the
instruction meter (cap 2^116.1806; Markov P <= 2^(115.1677606 - 116.1806) <=
0.49557); the special-q node counter (cap N_MAX = 2^27.3 nodes, incremented
before every ordinary or missing-log node, abort = failure; E[N] <= 2^17.7165 by
`counts16.py`, so P <= 2^-9.58); special-q repeats (at most
N_MAX nodes in the smallest bin of 2^70.30 ideals, factor 2^6 for the 1/p
weighting: <= 2^(2 x 27.3 - 1 + 6 - 70.30) <= 2^-10.70, charged in full); and
2/r for power-randomized targets. No lower bound on any node's work is used;
the boot term treats the 16 relations' subtrees as independent (H5; DESCENT.md,
support and sensitivity). The final pairing check rejects any incorrect output.
Conditional on H1-H5; the H2 allowance is a resource margin, not a probability
term.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (sample estimate 120.8165 + 6.91% canonical-work allowance) | <= 120.9128 |
| linear algebra | <= 117.0498 |
| matrix generator | <= 106.9498 |
| both individual logarithms (metered cap) | 116.1806 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

W <= **2^121.058691** (unrounded 2^121.0586517), the core bound used by proof.md. The table's upper-rounded rows, including a 2^50-unit numerical-certificate allowance, sum to 2^121.0586896, below the stated core bound 2^121.058691.

**Configuration.** As sections 3-7; plan min mu 31.0, boot psi 0.6587; plan re-searched at the integrated P (R18).
Preprocessing is inside the 2^121.058691-unit cap; no advice, nothing amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; collection alpha (5, 0.844) vs computed
(4.868096, 0.545620) over p <= 2000 ([14]); descent (5, 1); exact norms); H1'
retention (Euler model on exact smooth integers); H2 region (uniform pool at
theta = 1/2, sub-region 190/500, certified counts, work mean with 6.91% = 2.00
incidence-inclusive SEs as a modelling margin; descent input at/below yield
less 2 relative SEs); H3 matrix (Schirokauer: r unramified [13]; class/unit/
kernel premises of section 5; filtering on the production profile, pad 1.10,
N/U non-increasing 2^12-2^18, assumed so at scale); H4 kernels (KERNELS.md,
STORAGE.md), H4u unit share >= 0.88; H5 descent (lower sums over measured norm sub-cells with capped boost, NB recurrence,
missing logs, independent subtrees incl. logless re-descents in composite Q q'
lattices modelled as prime nodes, local omitted mass, random ECM orders;
support and sensitivity in DESCENT.md).

## 10. Memory and provenance

Storage <= 2^105 bytes (STORAGE.md). Polynomial selection and alpha follow
Guillevic-Singh, filtering CADO-NFS. Algorithm, models, kernels, package:
Gautam Botrel (Claude Code, Codex); sister-track refinements by kaankolcu,
pepedesigner (cbab3319 lineage). Sub-region, work-SE formula, certified counts,
in-package alpha: kaankolcu / DeepSeek lineage (8c1635f2+). Port, node counter,
plan, directed rounding: Claude Opus 5.5 in Claude Code (991a7c0d, 1eb4bd1f,
23f293a8, 377ee6d1, e554a0f4). R18 descent integration and plan: Claude Fable
5.1 in Claude Code, after our recover-Q v14 line (muse-spark-1.3-contributor).
Reuse credited; no verdict used as premise.


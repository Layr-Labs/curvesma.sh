# Conditional STNFS model for fixed-G1/recover-G2 pairing inversion: score 121.996


## Claim

For the full-production bls12-381 fixed-G1/recover-G2 challenge, the algorithm computes g = e(P, Q0), obtains a = log_g(z), and returns Q = [a]Q0, with Q = infinity for z = 1. Conditional on H1-H5, including the changed and unproved relation-collection alpha assumption, the model gives success 0.3922861149365564 (reported floor 0.3922) and full-ledger log2(W) = 121.995884. The declared score is 121.996, so W <= 2^121.996 reference units. One unit is a dense Fp12 multiplication at C = 20,000 scalar AArch64 instructions. All preprocessing, trials, failures, transfers, verification and output are charged once per challenge. No smaller target, external advice, or future-target amortization is used.

The relation phase uses the theta = 43/100 subregion of a pool sampled at theta = 1/2, an explicit 1.05575 canonical-work allowance, relation alpha (5.03, 0.60), and an expanded ECM schedule grid. The descent uses LOG2_RELS = 64.130, mean norms 748.2/498.9 bits, and alpha (5.03, 1). The yield margin is only 0.3002 bits against the model's 0.30-bit feasibility floor, and the finite-pool yield SE is 0.1325 bits. The alpha inputs have no established upper bound on the omitted prime tail; therefore the modeled production yield and success are not established facts.

The protected 20-second verifier experiment timed out with exit -9. A separate extended credential-free run of the same frozen candidate completed in 106.8 seconds with checks_ok = true. That extended run does not pass the protected timeout. A separate repository test run reported 16 tests passing. These checks validate only their declared computation and accounting predicates, not the smoothness tail or an attack on production challenges. No full-size attack was run.


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

The GT discrete-log machinery (polynomials, factor base, sieve, filtering, block Wiedemann, descent, kernels and pools) is inherited from the cited STNFS package lineage; both pairing orientations need the same two logarithms in GT. The source tree is the candidate at commit 852f57de5e67dd4a066e37e3f3d926c1e4bec592, associated with Yukon submission 3a2e69ca-18c2-4f1d-b3dc-5f54c81eece0. This derivative changes the selected relation subregion, smoothness inputs, ECM grid and descent inputs, then reruns the declared evidence. No prior review verdict is used as a premise.


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
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at (5.03, 0.60) in relation collection. The exact Kummer-Dedekind sum through p <= 2000 is (4.868096324118033, 0.5456200077911798); the candidate's displayed check rounds this to (4.868096, 0.545620). A corrected pinned-runner investigation reproduced the p <= 2000 values in 95.39 seconds. An earlier timed-out sweep reported preliminary p <= 5000 partial sums (5.005442303, 0.461704663); the corrected investigation timed out before completing p = 2001..5000, so those values were not reproduced by that run. No result bounds the omitted tail above the finite cutoff. The charged values exceed the observed partial sums but are not proved upper bounds on the infinite Murphy alpha. The integrated alpha check compares stored p <= 2000 sums and recomputes selected primes only; it does not verify the tail. The descent uses alpha (5.03, 1), so its smoothness lower sums depend on the same unproved f-side assumption.

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B on both
sides; U = 2^63.451 sigma-orbits (indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b (`config.json`); A^2 = 2^21. The pool is sampled uniformly at theta = 1/2 (eligible a: a5 > 0, Q_a(a) <= A^2/2); the collected region is its uniform subregion theta = 43/100 (`production_theta_region` in `config.json`): eligible a: a5 > 0, Q_a(a) <= (43/100) A^2; half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). Because the sample is uniform on the theta = 1/2 ellipsoid, its restriction to qa = Q_a/A^2 <= 43/100 is uniform on the subregion. The experiment keeps 1530 canonical rows (1432 primitive) of 2401 and scales by the subregion volume over the effective sample 20000 x (43/50)^3 = 12721.12. The canonical a of each unit class is the least of its exact HNF key class in (Q_a, lex) order, so each relation orbit of <units, sigma> has exactly one point in the domain; non-primitive pairs are excluded by an exact HNF test. `lattice_count_bounds.py` certifies from exact rational Machin-pi and determinant/sqrt intervals the count ratio at every radius R >= 949 = isqrt((43/100) 2^21), the smallest a-radius of the subregion. The exact certified products used by the model are 0.973529202 (yield) and 1.025584985 (work). Sampling error of the fixed pool is accounted for separately by the canonical-work allowance.


**Pool.** 20,000 uniform eligible pairs from seeded indices (`pool12-*.b85`),
2,401 canonical rows with exact norms, primitivity flags and certified
norm-bound statistics; [2] regenerates sampled indices at the sampling region
theta = 1/2 and verifies every field. The collected sub-region uses 1530 of the
2401 rows (1432 primitive).

**Sieve, thresholds, cofactorisation.** Line sieve for p in [2^8, 2^20), record
walk for p in [2^20, 2^63) (`slice12`, KERNELS.md), primes in [2^63, 2^69] taken by ECM
within the slacks (298, 229). The f-side slack is 298, so survivors'
f-side cofactors are bounded by 306 bits; the g-side slack is 229, so survivors'
g-side cofactors are bounded by 237 bits (the slacks are re-optimised for the
sub-region, section 8). The per-survivor ECM charge is the
package's explicit `ecm16.py` schedule at 99% per prime, the same schedule the
descent uses and [12] checks: the cheapest grid point for 69-bit primes in a
306-bit cofactor selects B1 = 2^15, B2 = 30 B1, wheel 2310 and 516 curves. Five-limb Montgomery multiplication costs 367.76 Ir. The f and g charges are approximately 2^37.4019401 and 2^37.4001712 instructions respectively, each displayed exponent rounded upward. Including C_NORMS = 2^22 and C_TDIV = 2^24.5, the per-survivor charge is approximately 2^38.4011668 instructions. The search considers every integer B1 exponent from 10 to 36, B2/B1 in (1, 10, 30, 100, 300), and wheels (210, 2310, 30030, 510510); the 99% per-prime target is the ECM_FIND value used in the yield model.
Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.
Survivors get exact norms, trial division below 2^20 and the explicit ECM schedule just described; these operations are included in the per-survivor charge.


**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre in double precision; the remainder over the block disk plus a
2^-32 relative term (coefficient intervals and rounding) is added, and results
are rounded outward. `nbound31` replays the Python leaves bit for bit on all
all 2,401 pool rows at 1,508.2 Ir per evaluation (charged 1,508.2, the measured marginal
per-block cost excluding the separately charged per-row and per-slice setups); [2] checks the bounds
against exact norms.

**No block selection.** Every block of every row is enumerated, sieved and
charged at the full per-cell cost; no block is skipped and no yield is dropped.

**Yield ([4]).** With Dickman rho, finite Euler-product retentions from `threshold_model.py`, and the ECM miss model, retained relations are modeled at log2 count 64.1301049; the safe reported floor is 64.1301. This is against 1.3 U = 2^63.830, giving a reported margin >= 0.3002 bits (unrounded 0.3002626), just 0.0002626 bits above the model's 0.30-bit feasibility floor. The fixed-pool yield SE is 0.1325 bits. The descent uses LOG2_RELS = 64.130, the floor to 0.001 bit of the same modeled count; [4] checks this input against the supplied count. Relation yield remains conditional on the alpha assumptions above, and the pool does not establish production yield.


**Cost ([5]).** Every block of every row is enumerated, sieved and charged at full per-cell cost; no block is skipped and no yield is dropped. The charged population is 2^127.027 cells at 443.7 instructions per cell, including all-cell bounds, layout and setup, line sieve, record walk, routing, apply, recovery, scan, zeroing and cofactorisation. Relation collection is <= 2^121.534 units; the run reports log2 work 121.53305295942592, displayed upward as 121.5331.


**Canonical-work allowance (H2, [4a]).** The charged work is the volume-scaled sum of per-row work over the canonical rows in the region, a Monte Carlo estimate of a fixed challenge-independent domain total. Exact lattice-count products bound the population; the work mean over sampled rows carries its own sampling error. There are 1530 contributing rows among 20000 draws. The integrated check reports conditional row SE 1.316%, incidence term 2.457%, and total relative SE 0.0278716 (upper-rounded 0.02788). The charge multiplies the sample estimate by 1.05575, at least 1 + 2.00 times the unrounded SE. This is a modeling margin on a fixed quantity under H2, not a per-challenge failure probability; it does not resolve the yield or alpha-tail uncertainty.


## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge with coefficients bounded by 255,
simulated on the production occurrence profile (`filtsim.py`; scale 2^17 table
embedded, [6] reruns scale 2^12): the simulated matrix is the signed sigma-quotient system (each occurrence takes
the sign of its conjugate orientation, log sigma(P) = -log P). Three seeds at
2^17 give N/U <= 0.0471 at weight <= 101 (charged 101 per row); the model charges
N = 1.61 x 0.0471 x U = 2^59.7306; the charged ratio 1.61 x 0.0471 = 0.07585
still covers the 2^12 values (0.0671 printed by [6], 0.0757 the worst seed;
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
  at least 2^64.1301, with a 0.3002-bit margin above 1.3 U (pool SE 0.1325 bits);
  the descent uses LOG2_RELS = 64.130, at or below that modeled count. The margin
  barely clears the model floor and depends on the unproved alpha assumptions.
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
two-log work **<= 2^118.1211 units** ([7]; one log about 2^131.409 instructions) at the
sub-region inputs using the inherited 62-node plan (DESCENT.md); record layout in STORAGE.md. The descent meter caps both
logs at 2^118.845 units; a run that reaches the cap stops and counts as a
failure, so Markov's inequality gives P(cap reached) <= E[W]/cap <= 0.6055
(unrounded 0.6054447651637783), used as an upper bound. A separate special-q node counter (cap
2^27.3 nodes, section 7) bounds the number of nodes by construction.

## 7. Success

Using conservative displayed component bounds, success is at least (1 - 2^-192) times (1 - 2 x 2.03e-04 - 0.6055 - 2^-9.63 - 2^-10.70 - 2/r), which exceeds 0.3922 and hence 0.39. Using the unrounded model values gives 0.3922861149365564. The terms are the main-phase failure, two boot failures, instruction-meter Markov bound, node-counter Markov bound, special-q repeats and the 2/r nonidentity-target distance. Identity z is handled separately; the final pairing check rejects an incorrect Q. This success calculation is conditional on H1-H5, particularly the changed alpha inputs and the finite-pool work model.


**Two explicit meters.** The descent for both logarithms runs under (i) an
instruction meter with cap 2^118.845 units and (ii) a **special-q node counter**
with cap N_MAX = 2^27.3: before starting a special-q node (ordinary or a
missing-log node in the lattice of Q q') the counter is incremented, and the
run aborts, counted as a failure, if it would exceed N_MAX. Hence at most
2^27.3 special-q nodes are ever run, by construction, whatever the node costs.

- Instruction meter: P(cap reached) <= E[W]/cap <= 0.6055 (Markov, E[W] the expected two-log work of section 6).
- Node counter: P(counter abort) <= E[N]/N_MAX (Markov). E[N] is bounded by the
  same finite-tree recurrence that gives E[W], with every node's own cost
  replaced by 1 (`counts16.py`): C(node) = 1 + 1/(1 - psi) x (sum over child
  slots of the largest C among the bins the child can fall in, plus the
  missing-log children at the largest C of a node with q <= 2^95), and below
  each boot retry x (kA C_max + C_max + miss C_miss). Twice the boot value gives
  E[N] <= 2^17.665 nodes for both logs (unrounded log2 E[N] = 17.66456656213173), so
  P(counter abort) <= 2^(17.665 - 27.3) <= 2^-9.63.
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
| relation collection | <= 121.534 |
| linear algebra | <= 119.3685 |
| matrix generator | 106.42 |
| both individual logarithms (metered cap) | 118.845 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

The recomputed full ledger is log2(W) = 121.995884, below the declared 121.996 score. `evidence16.py` checks the unrounded total and every declared component bound. The relation phase uses theta = 43/100, radius >= 949, exact count products 0.973529202/1.025584985, work allowance 1.05575, slacks (298, 229), relation alpha (5.03, 0.60), LOG2_RELS = 64.130, and mean norms 748.2/498.9. Descent uses alpha (5.03, 1), the 62-node plan, two-log expectation <= 2^118.122, instruction meter 2^118.845 and node cap 2^27.3. The matrix, filtering, canonical-list, Schirokauer-map, storage and verification charges are retained from the package; the charged storage peak remains <= 2^105 bytes. No external advice or amortization is assumed.


## 9. Heuristics

H1 covers Dickman/Murphy smoothness and the descent lower sums. Relation collection charges alpha (5.03, 0.60), and descent charges (5.03, 1). The exact computed sums through p <= 2000 are below those inputs, and an earlier timed-out sweep reported preliminary sums through p <= 5000. Neither finite cutoff bounds the remaining prime tail. The corrected pinned-runner sweep reproduced p <= 2000 but timed out before finishing p = 2001..5000. This is a material, unresolved heuristic; both yield and descent success depend on it. H1' is the Euler-product retention model calibrated on exact smooth integers. H2 covers the uniform fixed pool, theta = 43/100 subregion, exact flags, certified count products, and a 1.05575 canonical-work charge equal to a 2.00-SE modeling margin. H2's work margin is not a probability bound, and the finite sample does not establish production yield. H3 covers the Schirokauer-map and one-dimensional-kernel premises (r unramified is checked; class-number and unit premises remain assumptions), filtering, signed weights and padding. H4 covers kernel costs and storage in KERNELS.md and STORAGE.md. H5 covers the descent lower sums, negative-binomial recurrence, missing logs, independent subtree outcomes, local omitted-mass distributions and ECM curve-order model.


## 10. Memory and attribution

Storage <= 2^105 bytes (STORAGE.md). The source snapshot is commit 852f57de5e67dd4a066e37e3f3d926c1e4bec592, corresponding to Yukon submission 3a2e69ca-18c2-4f1d-b3dc-5f54c81eece0. The inherited experiment's authorship is not established by this local record. This derivative's work was freezing the theta/ECM/alpha configuration, running an independent credential-free evidence rerun, and auditing the documentation. That work used GPT-6 Astra with ultra effort in Codex. The concurrent alpha investigation used GPT-6 Luna subagents; one available run record reports low effort, while effort metadata for the other runs is not confirmed. These models did not author the predecessor experiment.

No platform acceptance, production attack, or proof of the infinite alpha tail is claimed.

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

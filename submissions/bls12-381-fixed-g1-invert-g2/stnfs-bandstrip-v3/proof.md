# STNFS at relation excess 1.1 on a theta = 2581/8000 sub-region of a 150k pool, line sieve to 2^48 with band strip, re-optimised descent plan, recentred slice sieve, f-first cofactorisation, weight-200 merge, ten LA owners: 2^120.6173625121048 reference units

Track bls12-381-fixed-g1-invert-g2. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used. Predecessor work is reused as credited below; no prior verdict is a premise.

## Claim

The algorithm below returns Q in G2 with e(P,Q) = z for the exact supplied nonidentity
P in G1 (Q = infinity for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins under the declared heuristics H1-H5;
the modelled success is at least 0.3904 (> 0.39; rounded below the model value 0.3904780). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^120.6173625121048 units**; the score is
**120.6173625121048**. The ledger (section 8): upper-rounded rows sum to at most 2^120.6173626 <= 2^120.6173625121048 (unrounded
model total 2^120.6172950). Relation collection is
charged on the theta = 2581/8000 sub-region of the sampled region with exact
certified lattice-count products and an explicit canonical-work allowance
(2.00 incidence-inclusive standard errors of the sampled work estimator, [4a]; 3.16%). Distributed storage is at most
2^105 bytes (STORAGE.md) in individually bounded, explicitly addressed objects.
The claim is analytic and exploratory under the NFS heuristics H1-H5 (section
9) with the measurements and checks cited; no full-size run is made.

**Numerical conventions.** Decimals are the experiment's values rounded as
shown. Figures used as bounds carry <= or >= and are rounded in the safe
direction; the `declared_bounds` check compares every stated bound (ledger,
collection, LA, relation count/margin, work SE, two-log, Markov, E[N],
counter-abort, repeat, boot failure, success) with the unrounded value and
fails if any is violated.

The model is `model16.py` and `descent16.py` with `config.json` and `stats16.b85`;
`evidence16.py` recomputes the ledger and the checks cited as [n].

## Base

The GT discrete-log machinery (polynomials, factor base, sieve, filtering, block Wiedemann, descent, kernels, pools) is the lean STNFS package of sister track `bls12-381-fixed-g2-invert-g1` (cbab3319), ported and refined up to e3940d11; both orientations need the same two GT logs. v1 of this line (640f4e61, 121.038) is v13 (57598408, 121.12966) with four cost-side charges of the parallel lineage ported unchanged: the re-measured slice-sieve kernels, f-first cofactorisation with early abort and the weight-200 merge of 2c2ce2ab (121.241), and the ten-active-owner row tail of 6282ab79 (121.315). This revision (v2) lowers the relation excess from 1.3 to 1.1 (new filtering measurements at excess 1.1, section 4), shrinks the region to keep a 2-SE margin over 1.1 U, and re-optimises the descent plan for the smaller relation count (30 node entries move by one or two steps; the boot is unchanged). v3 (this revision) adds the line sieve to 2^48 with a survivor band strip, validated on the sister track (8cc93834 lineage, 42e80f80), on top of the early abort. Changes: last section.

## 1. Reduction

Let g = e(P, Q0) (pinned `pairing(Q0, P)`), one pairing; Q0 is the fixed G2
generator. P has prime order r and e is nondegenerate, so g generates mu_r,
z = g^a for unique a mod r (a = 0 iff z = 1), and Q = [a]Q0 (one G2 scalar
mult), so e(P, Q) = z for that exact P (never replaced). log z and log g share
one relation/LA phase; a = log z / log g mod r.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks that
(u y)^((p^6-1)/2) = -1 in F_p6, so g is irreducible over F_p6 and
F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. The map
x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P), so
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at
(5.15, 0.85) in collection: [14] recomputes it in-package from the exact
Kummer-Dedekind structure over p <= 2000 (`--alpha` prints (4.868096,
0.843643)) and over p < 50000 (full sweep, 313 s: f 5.085035, g 0.731543).
The g side uses the primitive-correct p = 2 double-root term 2q/(q^2-1)
(1+q^-8) with the q/(q+1) base-primitivity factor. Charges 5.15 and 0.85 cover
5.085035 and 0.731543 (margins 0.065, 0.118); the descent uses the same charges.

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B on both
sides; U = 2^63.451 sigma-orbits (indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21. The pool is uniform on theta = 1/2 (a5 > 0,
Q_a(a) <= A^2/2); the collected region is its sub-region theta = 2581/8000
(`production_theta_region`): a5 > 0, Q_a(a) <= (2581/8000) A^2; half-slice
b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The restriction of the uniform sample to
qa <= 2581/8000 is uniform on the sub-region: 4792 canonical sub-region rows
(4389 primitive) of 17874, scaled by the sub-region volume over the effective
sample 150000 x (2581/4000)^3 = 40297.240. The canonical a of each unit class is
the least of its exact HNF key class in (Q_a, lex) order, so each relation orbit
of <units, sigma> has one point in the domain; non-primitive pairs are excluded
by an exact HNF test. `lattice_count_bounds.py` certifies (exact rational
Machin-pi, determinant/sqrt intervals) the count ratio at every radius
R >= 822 = isqrt((2581/8000) 2^21) (b-radii >= 1191): [0.984652, 1.014688] for
Q_a and [0.984602, 1.014688] for Q_b, inside [0.98, 1.02] (asserted); the model
uses the exact products **0.969491602 (yield, floor)** and **1.029590982
(work, ceiling)**. Pool sampling error is handled by the canonical-work
allowance below. The canonical list costs 2^90 units.

**Pool.** 150,000 uniform eligible pairs from seeded indices (seed 1: 2000, seed 2: 148000;
`pool12.b85`), 17,874 canonical rows with exact norms, primitivity flags and certified
norm-bound statistics: per-seed canonical index lists (delta-varint) plus a 64-row
field-anchor table; `load_pool` regenerates every row from its (seed, index) (v11
sampler, exact norms, HNF primitivity, certified nbound statistics). [2] checks list
membership both ways on sampled indices (canonical iff listed) and every anchor field.
All 150,000 indices were re-derived offline by an independent code path (exact match). The sub-region uses 4792 of the 17874 rows (4389 primitive).

**Sieve, thresholds, cofactorisation.** Line sieve (scores only) for p in [2^8, 2^48), record
walk for p in [2^48, 2^64) (`slice12`, KERNELS.md; X2 = 2^64), primes in [2^64, 2^69] taken by ECM
within the slacks (295, 226). **Line sieve to 2^48 with band strip (v3; 8cc93834 /
42e80f80 lineage):** primes in (2^20, 2^48] only add scores on the line path
(every line keeps >= 2^7.11 hits per coset at p = 2^48, inside the measured line
table) and emit no records: records/cell 2.33 -> 0.58. Their factors are recovered at the
survivors by a band strip: two rounds per side of the `ecm16.py` 48-bit schedule at
99.9% on the unstripped cofactors (caps 832 / 576 bits: 2 x 2^36.4488 f side,
2 x 2^35.3903 g side). A band prime missed in both rounds loses the relation: yield x
(1 - 41 x 0.001^2)(1 - 28 x 0.001^2) (at most floor(832/20) = 41 and 28 band
factors per side). The f-side slack is 295, so survivors'
f-side cofactors are bounded by 303 bits; the g-side slack is 226, so survivors'
g-side cofactors are bounded by 234 bits (the slacks are re-optimised for the
sub-region, section 8). The per-survivor ECM charge is the
package's explicit `ecm16.py` schedule at **99.9%** per prime (ECM_FIND = 0.999, raised from the
base's 99%; this is also the `ecm16.schedule` default the descent uses), which [12] checks: the cheapest
grid point for 69-bit primes in a 303-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 403 curves,
five-limb Montgomery multiplication at 367.76 Ir, `ecm12`) prices one side at 2^38.0906
instructions and the other (234 bits) at 2^38.0879. **Early abort (f side
first, ported from 2c2ce2ab):** a survivor gets exact norms, trial division and
the f-side ECM; the g-side ECM runs only if the f side factors completely. A
relation needs both sides, so the yield is unchanged. Per survivor the charge is
C_NORMS + C_TDIV + 2 x 2^36.4488 + 2^38.0906 + q_f (2 x 2^35.3903 + 2^38.0879) (the g-side band
strip also waits for the f side) with q_f = min(1, 100 rho(u_f) /
P(omitted_f <= slack_f)) per pool row: P(f smooth and survivor) <= rho(u_f), and
every cell with omitted mass <= slack survives, so rho/P(omitted <= slack)
estimates the conditional probability that the f side factors (heuristic under
H1, model mean ~1e-5); the factor 100 is a safety margin, and even q_f = 1 would
only restore the two-sided charge. The 99.9% figure is the same
ECM_FIND used in the yield's miss model. The higher per-prime success raises
the yield, which lets the region shrink (section 8).
Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.
Survivors get exact norms, trial division below 2^20, the band strip and the ECM
schedule above (f side first, g side only after a complete f factorisation).

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the
block centre in double precision; the remainder over the block disk plus a
2^-32 relative term (coefficient intervals and rounding) is added, and results
are rounded outward. `nbound31` replays the Python leaves bit for bit on all
17,874 pool rows at 1,508.2 Ir per evaluation (charged 1,508.2, the marginal
per-block cost excluding separately charged setups); [2] checks the bounds
against exact norms.

**No block selection.** Every block of every row is enumerated, sieved and
charged at the full per-cell cost; no block is skipped and no yield is dropped.

**Yield ([4]).** rels = (va_region/40297.240) vb (0.969491602) sum over canonical primitive sub-region
rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
the ECM miss and the band-strip loss. Retained relations >= 2^63.7216 against 1.1 U = 2^63.589
(excess 1.1, section 4): margin >= 0.1327 bits (unrounded 0.13278; pool SE ~0.061, so
2.18 SE), above the 0.11-bit
feasibility floor (`MARGIN`). The descent's missing-log input `LOG2_RELS =
63.721` is the floor (to 0.001 bit) of this retained count (unrounded
2^63.72161; certified lower product already applied). [4] checks `LOG2_RELS <=
log2_rels` (`rels_input` check). The mean-norm inputs 744.1/500.0 are rounded
down from 744.18/500.03 and checked in [4].

**Cost ([5]).** Sampled mean 309.6 instructions/cell; charged mean 319.4
with the allowance below, over 2^126.227 cells: all-cell bounds, layout,
setups, line path and record walk at the re-measured recentred `slice12` costs
(KERNELS.md; primes in [2^63, 2^64) charged twice), routing
(10.01/record, 8/cell; STORAGE.md), apply 7.00, recovery 7.003/record, scan,
zero, band strip and f-first cofactorisation. Collection: **<= 2^120.2580 units** (estimate
2^120.2130 times allowance 1.0316; v2 charged 390.3/cell with the line sieve to
2^20, 2^120.5475).

**Canonical-work allowance (H2, [4a]).** The charged work is (va_region/40297.240)
times the per-row sum over in-region canonical rows: a Monte Carlo estimate of
a fixed, challenge-independent domain total. The estimator sums all 150,000
draws (zero outside the sub-region or non-canonical), so its relative SE has a
conditional row-mean term and an incidence term for the 4792 contributing
draws: conditional 0.684%, incidence sqrt((1-4792/150000)/4792) = 1.421%, total
<= 1.578%. The charge multiplies the estimate by **1.0316 >= 1 + 2.00 x 0.01578
= 1.03156** (rounded up; unrounded need 1.0315426); [4a] fails if the allowance falls below 1 + 2.00 SE.
This is a declared modeling margin under H2, not a probability term, and no
sampling tail is charged against success.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge (coefficients <= 255), simulated on
the production occurrence profile (`filtsim.py`; scale 2^17 table embedded, [6]
reruns 2^12) as the signed sigma-quotient system (each occurrence takes the
sign of its conjugate orientation). **Required excess 1.1 (v2).** There is no
early stop: the whole region is sieved (every cell is charged, section 3) and
every retained relation (>= 2^63.7216, i.e. at least 1.1 x 2^0.1327 = 1.206 U)
is passed to filtering, linear algebra and the descent (LOG2_RELS = 63.721 is
this same retained count). The region is sized so that the count clears the
requirement 1.1 U by 2.18 SE; the matrix is charged at excess 1.1, which
over-states N at the larger actual excess (N/U falls with the excess at every
measured scale, e.g. 2^15 maxima 0.0648 / 0.0570 / 0.0516 at 1.1 / 1.2 / 1.3).
All filtering data below are measured at excess 1.1 with the descent's mean
norms 744.1/500.0 (the 1.3
tables are superseded). Four seeds at 2^17 give N/U <= 0.05934 at mean weight
<= 101; the merge then continues to mean weight 200 (charged 201/row, ported
from 2c2ce2ab). N/U falls with scale (`filtsim.run`, kmax 64, weight 100; max
over seeds): 0.1091 (2^12, 10 seeds), 0.0820 (2^13, 6), 0.0666 (2^14, 6),
0.0648 (2^15, 6), 0.0631 (2^16, 6), 0.0593 (2^17, 4), 0.0551 (2^18, 2)
(`SCALE_NU`, [10] checks non-increasing). The model charges N = 0.80 x 1.092 x
0.05934 x U = 2^59.1814: ratio 0.06479 exceeds every weight-100 measurement from
2^15 up, and the weight-200 ratio N(200)/N(100) is charged 0.80 against measured
maxima 0.774-0.7995 at 2^12..2^18 (`W200` in model16.py, rounded up to 3
decimals; [10] checks max <= 0.80), falling with scale. At excess 1.3 the same
charge was 0.05163; the lower excess costs 2 x log2(0.06479/0.05163) = 0.655
bits of linear algebra and saves 0.32 bits of collection. Excess 1.2 and 1.0
were also measured (not used; see the note). N/U staying at or below
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
without multiplying, >= 88% of row weight is +-1 by H4u, measured 0.917..0.946 at
weight 100 and 0.919..0.945 at weight 200 at excess 1.3, re-checked at excess 1.1 by
`unit_entries.py`; cost is linear in N x mean weight),
21 dense terms at 97, and row tail C_ROW = 1,640 (N < 2^59.19 < 10 x 2^56 = 2^59.32, so
ten active 2^56-element owners: nine sends and nine combinations, counted 1,638
by [10] from `row_tail.c`; ported from 6282ab79; valid for any row, sums < 2^323). Every step forms the m
x n projections x^T A^i y: (m + n) N^2 dense terms at 97 plus one reduction
(C_RED9 = 420, counted 385) per entry. Horner reconstruction (w <- A w + Y F_i
over deg F <= N/m + N/n coefficients) costs one product plus, per element, n
dense terms and one reduction. With c = 201 x 20.7 + 1,640 + 21 x 97 per row:
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n + 420) + 97 (m + n)) + 420 (m + n) N
instructions: **<= 2^117.8823 units** at N = 2^59.1814 (blocks (m, n) = (19, 4)
retained; `owner10.check` asserts them). The (m + n) x (m + n) generator of the length L = (N +
21)(1/m + 1/n) + 2 sequence (divide-and-conquer approximant bases, Karatsuba):
<= 4 (m + n)^3 L^log2(3) log2(L) field mults at 97 plus reduction 420:
**<= 2^107.2086 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** Region and polynomials are coin/challenge-independent, so the
  count is one fixed number; H1 puts it at >= 2^63.7216, >= 0.1327 bits above
  1.1 U (pool SE ~0.061, so ~2.18 SE above requirement); all of them are
  filtered and the descent uses 63.721, at or below that count.
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
survivor-test retention is >= 0.8593 at every sieve node ([7], one Hensel bound
per root below 2^9); cofactors use the `ecm16.py` schedule ([12]). Relation
counts are negative binomial (toy-fit shapes). Log-less ideals are discarded or
re-descended in a fresh lattice. The shipped plan (62 nodes + boot) gives
expected two-log work **<= 2^116.0558 units** ([7]; one log <= 2^129.343
instructions) at the excess-1.1 inputs (LOG2_RELS = 63.721, mean norms
744.1/500.0; coordinate-search plan re-optimised for these inputs, UDUP = 0.65).
The meter caps both logs at 2^116.785 units (cap reached = failure), so
P(cap) <= E[W]/cap <= 2^(116.0558-116.785) <= 0.6032384 (Markov). A node
counter (cap 2^27.3, section 7) bounds node count by construction.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 1.807e-03 - 0.6032384 - 2^-8.89 - 2^-10.70 - 2/r)
>= **0.3904** >= 0.39 (audit value 0.3904387; model 0.3904780): main phase, two boots
(phi <= 1.807e-03 each at the excess-1.1 relation input), the instruction-meter
Markov bound, the node-counter Markov bound and special-q repeats, as follows.

**Two explicit meters.** Both logs run under (i) an instruction meter (cap
2^116.785 units) and (ii) a **special-q node counter** (cap N_MAX = 2^27.3):
the counter increments before each special-q node (ordinary or missing-log)
and the run aborts as a failure past N_MAX. Hence <= 2^27.3 special-q nodes
ever run, by construction.

- Instruction meter: P(cap reached) <= E[W]/cap <= 2^(116.0558 - 116.785) <=
  0.6032384 (Markov, E[W] the expected two-log work of section 6).
- Node counter: P(abort) <= E[N]/N_MAX (Markov). E[N] uses the same
  finite-tree recurrence as E[W] with unit node cost (`counts16.py`); twice
  the boot value gives E[N] <= 2^18.41 nodes both logs (unrounded 2^18.4046),
  so P(abort) <= 2^(18.41-27.3) = 2^-8.89 <= 0.0021076.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for
  1/p weighting: two nodes share a special-q with probability <= C(N_MAX,2)
  2^6/2^70.30 <= 2^(2x27.3-1+6-70.30), rounded up to 2^-10.70 <= 0.0006012,
  charged in full.

N_MAX = 2^27.3 is kept from v13. No lower bound on
individual node work is used. The boot term treats the 18 relations' subtrees
as independent (H5). Targets are uniform nonidentity by power randomization;
distance from the uniform model is at most 1/r per
target and 2/r is charged above. Identity z is handled separately; the final
pairing check e(P, Q) = z rejects any incorrect Q. The success bound is
conditional on H1-H5, including the H2 canonical-work allowance of section 3,
which is a resource-model margin and not a probability term.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (excess 1.1, band strip; sample estimate 120.2130 + 3.16% canonical-work allowance) | <= 120.2580 |
| linear algebra (excess 1.1, weight 201, N = 2^59.1814, ten owners) | <= 117.8823 |
| matrix generator | <= 107.2086 |
| both individual logarithms (metered cap) | 116.785 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

The upper-rounded rows sum to 2^120.6173626, so W <=
**2^120.6173625121048** (audit-checked), under the claim; unrounded model total
2^120.6172950.

**Configuration.** Relation excess 1.1; sub-region theta = 2581/8000 (4792 of 17874 canonical rows),
exact certified count products at radius >= 822, work allowance 1.0316
(2.00 SE), X2 = 2^64, cofactor slacks (295, 226), relation-collection alpha (5.15, 0.85),
descent inputs LOG2_RELS = 63.721 and mean norms 744.1/500.0, descent alpha
(5.15, 0.85). The 62-node plan is v13's coordinate-search plan (REFINEMENT.md)
re-optimised by the same kind of coordinate search at the excess-1.1 inputs
(30 entries move by one or two steps, listed in refinement.json; boot unchanged),
evaluated at UDUP = 0.65, THIN_MAX 85, MIN_RELS 18; at present inputs min mu 31.2, max phi < 0.05,
boot psi 0.705 (largest). Two-log expectation <= 2^116.0558, meter 2^116.785,
plus the section-7 node counter. Norm-bound charge 1,508.2, factor-base bound,
filtering, generator, canonical list, Schirokauer and verification charges are
the base package's; LA charges (20.7/entry, 420/reduction, 97/dense), blocks
(19, 4) and the 10.01 route are this lineage's; excess 1.1, pad 1.092 and the
re-optimised plan are this revision's; weight 201 with the
0.80 N ratio, C_ROW 1,640 (ten owners), the recentred kernel table and the
early-abort cofactorisation are the ported charges (STORAGE.md, KERNELS.md).
Preprocessing is inside the 2^120.6173625121048 cap. No external
advice; nothing amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; relation collection alpha (5.15, 0.85) against the in-package computed (4.868096, 0.843643) over p <= 2000 and (5.085035, 0.731543) over p < 50000 ([14]); descent alpha (5.15, 0.85) with the same support; exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool at theta = 1/2 (150,000 draws, seeds 1-2), uniform sub-region theta = 2581/8000, exact flags,
certified sub-region counts with their exact rational products, and the
canonical work mean charged with a 3.08% modeling-margin allowance >= 2.00
incidence-inclusive standard errors of the total-work estimator, not a
probability term); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
signed filtering on the production occurrence profile at relation excess 1.1,
pad 1.092 over the 2^17 maximum, N/U (max over seeds) non-increasing 2^12..2^18
and assumed not to rise at production; the weight-200 ratio N(200)/N(100) <= 0.80
not rising at production scale is part of H3, measured 0.774-0.7995 at
2^12..2^18 at excess 1.1); H4 kernel cost (KERNELS.md, STORAGE.md; LA tail/reductions are
statically counted loop-free code, the route an executed-instruction count),
with H4u: >= 88% of production row weight is +-1 (measured 0.917..0.946 at
weight 100, 2^12..2^15, and 0.919..0.945 at weight 200, 2^12..2^14, at excess
1.3; 0.921..0.947 at weight 200 and excess 1.1, 2^12..2^15, `unit_entries.py`;
a shortfall changes cost only); f/g independence given
survival (H1) bounds the early-abort probability q_f (charged at 100x the
estimate); band-strip ECM finds band primes at the `ecm16` per-prime rate (random
group orders, as in collection);
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing logs,
independent subtrees, fresh-lattice re-descents, local omitted-mass laws,
random-integer ECM orders, UDUP = 0.65 (relation-weighted E[shat/m]/pr >=
1.06 boot / >= 4.70 nodes, uniform cells, at v13's node regimes; the boot and the
binding (220, g) node are unchanged here and re-checked live by [11]; spot cells at
re-optimised nodes are listed in DESCENT.md)).

## 10. Memory

Storage <= 2^105 bytes (STORAGE.md).

## Provenance

Algorithm, models, kernels and package: Gautam Botrel (with Claude Code and
Codex); polynomial selection and alpha follow Guillevic-Singh, filtering
CADO-NFS. Refinements: kaankolcu, pepedesigner, mpjunior92 and the DeepSeek
lineage (cbab3319, 8c1635f2, e3940d11, 04073dd1, be446c45); grid/plan
refinements, both-bars alpha, UDUP 0.65 evidence and the 150k pool:
muse-spark-1.3-contributor (v13). Recentred kernels, early abort, weight 200:
gbotrel 2c2ce2ab. Ten owners: bernardo-lima12 6282ab79. Line sieve to 2^48
with band strip: 8cc93834 / a31a5f1d (kaankolcu) and 42e80f80 (gbotrel), sister
track. Merge, excess 1.1, plan re-optimisation and the band-strip port: Claude
Opus 5.5 in Claude Code. No prior review verdict is a premise.

## Changes from e3940d11 (121.9944)

v13 (57598408, 121.12966; history in its note and REFINEMENT.md). v1 (640f4e61,
121.038): ported charges of 2c2ce2ab / 6282ab79. v2 (9eae858f, 120.8479): excess
1.1, theta 2581/8000, re-optimised plan. v3 (this revision): line sieve to 2^48
with band strip, collection 2^120.5475 -> 2^120.2580, total 2^120.6173625121048.

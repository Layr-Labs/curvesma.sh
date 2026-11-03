# STNFS at relation excess 1.1 with per-side factor-base bounds 2^69.25 (f) / 2^69 (g) on a theta = 2284/8000 sub-region of a fresh 164k on-region pool, line sieve to 2^48 with band strip, sub-bin pr, skewed boot, est-2SE descent input, explicit purge-loss charge, re-measured filtering, ten LA owners: 2^120.213628275286 reference units

Track bls12-381-fixed-g1-invert-g2. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used. Predecessor work is credited below; no prior verdict is a premise.

## Claim

The algorithm below returns Q in G2 with e(P,Q) = z for the exact supplied nonidentity
P in G1 (Q = infinity for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins under the declared heuristics H1-H5;
modelled success is at least 0.4968 (> 0.39; section 7; model value 0.4988719). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^120.213628275286 units**; the score is
**120.213628275286**. The ledger (section 8): the upper-rounded rows sum to
2^120.2136282752853289 (70-digit decimal arithmetic), below the declared value, the least
12-decimal number above that sum. Relation collection is
charged on the theta = 2284/8000 sub-region of the sampled region with exact
certified lattice-count products and an explicit canonical-work allowance
(2.00 SE of the sampled work estimator, [4a]; 1.73%). Distributed storage is at most
2^105 bytes (STORAGE.md).
Analytic exploratory claim under H1-H5 (section 9); no full-size run is made.

**Numerical conventions.** Decimals are the experiment's values rounded as
shown. Figures used as bounds carry <= or >= and are rounded in the safe
direction; the `declared_bounds` check compares every stated bound (ledger,
collection, LA, relation count/margin, work SE, two-log, Markov, E[N],
counter-abort, repeat, boot failure, success) with the unrounded value and
fails if any is violated.

The model is `model16.py`/`descent16.py` with `config.json`/`stats16.b85`; `evidence16.py` recomputes the ledger and the checks cited as [n].

## Base

The GT discrete-log machinery is the lean STNFS package of the sister track (cbab3319), ported and refined up to e3940d11. Lineage: v13 (57598408) plus four ported cost-side charges (re-measured slice-sieve kernels, f-first early abort, weight-200 merge from 2c2ce2ab; ten-owner row tail from 6282ab79) is v1; v2 lowers the excess to 1.1 with new filtering measurements and a re-optimised plan; v3 adds the 2^48 line sieve with survivor band strip; v4 adds sub-bin pr, skewed boot, the est-2SE input and a plan re-search. v5 re-samples the pool on the region (164,000 draws at theta 5/16), cuts the heavy-tailed central core from the yield, sizes the region by margin >= max(0.06, 2 SE), and re-measures filtering there. v6 (this revision) sets one factor-base bound per side (2^69.25 on f, 2^69 on g), re-sizes the region on the same pool, re-measures filtering on the per-side system, charges the >= 3-occurrence purge losses explicitly in the descent's missing-log model, and re-evaluates the descent at the new inputs. Changes: last section.

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
0.843643)) and over p < 50000 (full sweep: f 5.085035, g 0.731543).
The g side uses the primitive-correct p = 2 double-root term 2q/(q^2-1)
(1+q^-8) with the q/(q+1) base-primitivity factor. Charges 5.15 and 0.85 cover
5.085035 and 0.731543 (margins 0.065, 0.118); the descent uses the same charges.

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B_s-smooth on each side, B_f = 2^69.25 and
B_g = 2^69 (one bound per side, `lb_f`/`lb_g`). The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B_s on side s;
U = (li(2^69.25) + li(2^69))/2 + 2^31 = 2^63.5788 sigma-orbits (`log2_U_of`; indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21. The pool is uniform on theta = 5/16 (a5 > 0,
Q_a(a) <= (5/16) A^2); the collected region is its sub-region theta = 2284/8000
(`production_theta_region`): a5 > 0, Q_a(a) <= (2284/8000) A^2; half-slice
b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The restriction of the uniform sample to
qa <= 2284/8000 is uniform on the sub-region: 14987 canonical sub-region rows
(13765 primitive) of 19704, scaled by the sub-region volume over the effective
sample 164000 x (2284/2500)^3 = 125058.185. Each relation orbit of <units, sigma> has one canonical a (least of its HNF key class); non-primitive pairs are excluded by an exact HNF test. `lattice_count_bounds.py` certifies (exact rational
Machin-pi, determinant/sqrt intervals) the count ratio at every radius
R >= 773 = isqrt((2284/8000) 2^21) (b-radii >= 1224): [0.983684, 1.015625] for
Q_a and [0.983631, 1.015625] for Q_b, inside [0.98, 1.02] (asserted); the model
uses the exact products **0.967583931 (yield, floor)** and **1.031493518
(work, ceiling)**. Pool sampling error is handled by the canonical-work
allowance below. The canonical list costs 2^90 units.

**Pool.** 164,000 uniform eligible pairs from seeded indices (seed 1: 82000, seed 2: 82000;
`pool12.b85`), 19,704 canonical rows with exact norms, primitivity flags and norm-bound
statistics: per-seed canonical index lists (delta-varint) plus a 64-row
field-anchor table; `load_pool` regenerates every row from its (seed, index) (v11
sampler, exact norms, HNF primitivity, nbound statistics). [2] checks list
membership both ways on sampled indices (canonical iff listed) and every anchor field.
Seeds and draw counts were fixed before any yield was computed; the pool file reloads through `load_pool` with exact index match.
The sub-region uses 14987 of the 19704 rows (13765 primitive).

**Sieve, thresholds, cofactorisation.** Line sieve (scores only) for p in [2^8, 2^48), record
walk for p in [2^48, 2^64) (`slice12`, KERNELS.md; X2 = 2^64), primes in [2^64, 2^69.25] (f) and [2^64, 2^69] (g) taken by ECM
within the slacks (295, 226). **Line sieve to 2^48 with band strip (v3; 8cc93834 /
42e80f80 lineage):** primes in (2^20, 2^48] only add scores on the line path
(every line keeps >= 2^7.30 hits per coset at p = 2^48, inside the measured line
table) and emit no records: records/cell 2.33 -> 0.58. Their factors are recovered at the
survivors by a band strip: two rounds per side of the `ecm16.py` 48-bit schedule at
99.9% on the unstripped cofactors (caps 832 / 576 bits: 2 x 2^36.4488 f side,
2 x 2^35.3903 g side). A band prime missed in both rounds loses the relation: yield x
(1 - 41 x 0.001^2)(1 - 28 x 0.001^2) (at most floor(832/20) = 41 and 28 band
factors per side). Slacks (295, 226) bound survivors' cofactors by 303/234 bits (re-optimised for the sub-region, section 8). The per-survivor ECM charge is the
package's explicit `ecm16.py` schedule at **99.9%** per prime (ECM_FIND = 0.999; also the descent default), checked by [12]: the cheapest
grid point for 69.25-bit primes in a 303-bit cofactor (f side: B1 = 2^16, B2 = 30 B1, wheel 30030, 417 curves,
five-limb Montgomery multiplication at 367.76 Ir) prices that side at 2^38.1386
instructions and the g side (69-bit primes, 234 bits: B1 = 2^16, B2 = 30 B1, 403 curves) at 2^38.0879. **Early abort (f side
first, ported from 2c2ce2ab):** a survivor gets exact norms, trial division and
the f-side ECM; the g-side ECM runs only if the f side factors completely. A
relation needs both sides, so the yield is unchanged. Per survivor the charge is
C_NORMS + C_TDIV + 2 x 2^36.4488 + 2^38.1386 + q_f (2 x 2^35.3903 + 2^38.0879) (the g-side band
strip also waits for the f side) with q_f = min(1, 100 rho(u_f) /
P(omitted_f <= slack_f)) per pool row, u_f = (n_f + alpha_f)/69.25: P(f smooth and survivor) <= rho(u_f), and
every cell with omitted mass <= slack survives, so rho/P(omitted <= slack)
estimates the conditional f-side factorisation probability (H1 heuristic); the factor 100 is a safety margin, and even q_f = 1 would
only restore the two-sided charge. The same 99.9% ECM_FIND is used in the yield's miss model; the higher per-prime success lets the region shrink (section 8).
Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.
Survivors get exact norms, trial division below 2^20, the band strip and the ECM schedule above.

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the block centre; the remainder over the block disk plus a
2^-32 relative term is added, and results
are rounded outward. `nbound31` replays the Python leaves bit for bit (cross-checked on the v13 pool; this pool's leaves regenerate by the same procedure, 64 anchors checked by [2]) at 1,508.2 Ir per evaluation (charged 1,508.2); [2] checks the bounds
against exact norms.

**No block selection.** Every block of every row is enumerated, sieved and charged at full per-cell cost.

**Yield ([4]).** rels = (va_region/125058.185) vb (0.967583931) sum over canonical primitive sub-region
rows with qa >= 0.02 (the central-core cut, below) of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g), u_s = (n_s + alpha_s)/lb_s per side, with Dickman rho
([1]), retentions r from the finite Euler-product threshold model per side
(`threshold_model.py` at B_f / B_g, calibrated on exact smooth integers at B = 2^10..2^14)
the ECM miss (large primes in (2^64, 2^lb_s]) and the band-strip loss. Retained relations >= 2^63.8197 against 1.1 U = 2^63.7163
(excess 1.1, section 4): margin >= 0.1034 bits (unrounded 0.103450; pool SE 0.0515, so
2.01 SE), above the margin rule max(0.06, 2 SE) (`MARGIN_FLOOR`, `MARGIN_SE`).
The yield estimator is heavy-tailed near a = 0, so the 8 rows with qa < 0.02
contribute no yield (cells still charged): dropping a non-negative stratum only lowers the estimate. The descent's missing-log input `LOG2_RELS =
63.716` is the floor (to 0.001 bit) of the retained-count estimate minus 2 pool
SE (unrounded 2^63.819798 - 2 x 0.0515; certified lower product already
applied). [4] checks `LOG2_RELS <= log2_rels` and [7s] checks the estimate - 2 SE form (k = 0..3 sweep). The mean-norm inputs 737.7/498.4 are rounded
down from 737.70/498.47 and checked in [4], as are the per-side bounds.

**Cost ([5]).** Sampled mean 304.1 instructions/cell; charged mean 309.4
with the allowance below, over 2^125.861 cells: all-cell bounds, layout,
setups, line path and record walk at the re-measured recentred `slice12` costs
(KERNELS.md; primes in [2^63, 2^64) charged twice), routing
(10.01/record, 8/cell; STORAGE.md), apply 7.00, recovery 7.003/record, scan,
zero, band strip and f-first cofactorisation. Collection: **<= 2^119.8467 units** (estimate
2^119.8219 times allowance 1.0173).

**Canonical-work allowance (H2, [4a]).** The charged work is (va_region/125058.185)
times the per-row sum: a Monte Carlo estimate of a fixed domain total. The estimator sums all 164,000
draws (zero outside the sub-region or non-canonical), so its relative SE has a
conditional row-mean term and an incidence term for the 14987 contributing
draws: conditional 0.367%, incidence sqrt((1-14987/164000)/14987) = 0.779%, total
<= 0.861%. The charge multiplies the estimate by **1.0173 >= 1 + 2.00 x 0.00861
= 1.01722** (rounded up; unrounded need 1.0172143); [4a] fails if the allowance falls below 1 + 2.00 SE.
This is a declared H2 modeling margin, not a probability term.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge (coefficients <= 255) on
the production occurrence profile (`filtsim.py`; 2^17 table embedded, [6]
reruns 2^12), as the signed sigma-quotient system. **Required excess 1.1 (v2).** There is no
early stop: the whole region is sieved and
every retained relation (>= 2^63.8197, i.e. at least 1.1 x 2^0.1034 = 1.182 U)
is passed to filtering, linear algebra and the descent (the descent evaluates
at the conservative input LOG2_RELS = 63.716, at or below that count). The region clears
1.1 U by 2.01 SE; the matrix is charged at excess 1.1, which
over-states N at the larger actual excess (N/U falls with the excess: 2^15 maxima 0.0616 / 0.0570 / 0.0516 at 1.1 / 1.2 / 1.3).
All filtering data below are measured at excess 1.1 on the per-side system (ideals of side s sampled up to 2^lb_s, share li(2^lb_s)/(li(2^lb_f) + li(2^lb_g)); `filtsim.build_system`) with the relation-weighted mean norms 737.7/498.5. Four seeds at 2^17 give N/U <= 0.05653 at mean weight
<= 101; the merge then continues to mean weight 200 (charged 201/row, ported
from 2c2ce2ab). N/U falls with scale (`filtsim.run`, kmax 64, weight 100; max
over seeds): 0.1179 (2^12, 10 seeds), 0.0761 (2^13, 6), 0.0662 (2^14, 6),
0.0619 (2^15, 6), 0.0579 (2^16, 6), 0.0565 (2^17, 4), 0.0529 (2^18, 2)
(`SCALE_NU`, [10] checks non-increasing). The model charges N = 0.801 x 1.0954 x
0.0565292 x U = 2^59.2453: ratio 0.06193 exceeds every weight-100 measurement from
2^15 up, and the weight-200 ratio N(200)/N(100) is charged 0.801 against measured
maxima 0.775-0.801 at 2^12..2^18 (`W200`, rounded up to 3
decimals, falling with scale; the charge is the largest table value; [10] checks max <= 0.801). At excess 1.3 the symmetric-bound
charge was 0.05163 (cost 2 x log2(0.06193/0.05163) = 0.524
bits of LA). N/U not rising above
the 2^15-2^18 values at production scale is part of H3. The system solved is the square
A = [M | S]: N + 21 filtered rows, N ideal columns M, 21
dense Schirokauer columns S. Columns are sigma-quotiented, so only sigma-odd
units need maps: signatures K_f (0, 36), K_f^+ (0, 18), K_g (2, 5), K_h (2, 2)
give odd unit ranks 35 - 17 = 18 and 6 - 3 = 3 ([13], Sturm); under H3 the
kernel of A is the line of virtual logs. Block Wiedemann uses random blocks (m
= 19 left, n = 4 right, uniform in F_r): N/m + N/n
block steps, N(n/m + 1) products by A. One product costs per row 201 x
20.7 sparse (STORAGE.md: 4 + 3 + 0.88 x 11 + 0.12 x 33 = 20.64; >= 88% of row weight is +-1 by H4u, measured 0.933..0.950 at
weight 100 and 0.932..0.952 at weight 200 on the per-side system, `unit_entries.py`),
21 dense terms at 97, and row tail C_ROW = 1,640 (N < 2^59.25 < 10 x 2^56 = 2^59.32, so
ten active 2^56-element owners, counted 1,638
by [10] from `row_tail.c`; ported from 6282ab79; `owner10.check` asserts the owner count). Every step forms the m x n projections: (m + n) N^2 dense terms at 97 plus one reduction
(C_RED9 = 420, counted 385) per entry. Horner reconstruction costs one product plus, per element, n
dense terms and one reduction. With c = 201 x 20.7 + 1,640 + 21 x 97 per row:
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n + 420) + 97 (m + n)) + 420 (m + n) N
instructions: **<= 2^118.0102 units** at N = 2^59.2453 (blocks (m, n) = (19, 4)
retained; `owner10.check` asserts them). The (m + n) x (m + n) generator of the length L = (N +
21)(1/m + 1/n) + 2 sequence (divide-and-conquer approximant bases, Karatsuba):
<= 4 (m + n)^3 L^log2(3) log2(L) field mults at 97 plus reduction 420:
**<= 2^107.3116 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and LA
returns the filtered columns' logarithms.

- **Relations.** Region and polynomials are coin/challenge-independent, so the
  count is one fixed number; H1 puts it at >= 2^63.8197, >= 0.1034 bits above
  1.1 U (pool SE 0.0515, so 2.01 SE above requirement); all of them are
  filtered and the descent uses 63.716, at or below that count.
- **Linear algebra.** H3 (standard NFS-DL): r is unramified in K_f and K_g
  (checked by [13]), divides neither class number, the 21 sigma-odd Schirokauer
  maps are injective on the units modulo r-th powers, and the kernel of A is the
  line of virtual logarithms. Block Wiedemann with blocks uniform in F_r then fails
  with probability <= 4N/r < 2^-192 (Schwartz-Zippel, degree <= 2N).
  A failed solve (checked on the dropped rows) counts as failure.
- **Reconstruction.** Purged ideals' logs are recovered by back-substitution.
  The descent counts every ideal with <= 2 occurrences as log-less (coefficient
  1, so the zero-occurrence mass is charged in full) and charges ideals with
  >= 3 occurrences as lost at the per-band toy rate `H3_RATE` (Wilson 95% upper,
  three seeds; H5n); [6b] measures ideals left
  without a log on the [6] system, and the model's count is >= 7.24x measured
  in every band with loss; [6c] recomputes the pooled >= 3-occurrence losses
  (18 of 7186) and the shipped rates (DESCENT.md).

Under H1 and H3 the main phase succeeds with probability >= 1 - 2^-192 over the coins.

## 6. Individual logarithms

DESCENT.md specifies the descent: an ECM-trial boot on balls of certified p^6
lattices reduced in the weighted metric (b-block x isqrt(U)), then special-q nodes sieving an exact lattice ball (`dsieve16`) and
recursing on child ideals. Counted cells are base-primitive (0.9100);
survivor-test retention is >= 0.8593 at every sieve node ([7]); cofactors use the `ecm16.py` schedule ([12]). Relation
counts are negative binomial (toy-fit shapes); log-less ideals are discarded or
re-descended. The shipped plan (62 nodes + boot) gives
expected two-log work **<= 2^112.119 units** ([7]; one log <= 2^125.407
instructions) at the conservative inputs (LOG2_RELS = 63.716, mean norms
737.7/498.4, per-side bounds; coordinate-search plan re-optimised at 63.599 and re-evaluated here, UDUP = 0.65,
sub-bin pr, skewed boot reduction). The meter caps both logs at 2^113.119 units
(cap reached = failure), so P(cap) <= E[W]/cap <= 2^(112.119-113.119) <= 0.5 (Markov); a node
counter (cap 2^27.3, section 7) bounds node count.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 1.000e-03 - 0.5 - 2^-10.83 - 2^-10.70 - 2/r)
>= **0.4968** >= 0.39 (audit value 0.4968495; model 0.4988719): main phase, two boots
(phi <= 1.000e-03 each, 499x the H5 model value), the meter, counter and repeat bounds, as follows.

**Two explicit meters.** Both logs run under (i) an instruction meter (cap
2^113.119 units) and (ii) a **special-q node counter** (cap N_MAX = 2^27.3):
the counter increments before each special-q node
and the run aborts past N_MAX. Hence <= 2^27.3 nodes run, by construction.

- Instruction meter: P(cap reached) <= E[W]/cap <= 2^(112.119 - 113.119) <=
  0.5 (Markov; the cap is one bit above the expectation).
- Node counter: P(abort) <= E[N]/N_MAX (Markov). E[N] uses the E[W] recurrence with unit node cost (`counts16.py`); twice
  the boot value gives E[N] <= 2^16.47 nodes both logs (unrounded 2^16.4611),
  so P(abort) <= 2^(16.47-27.3) = 2^-10.83 <= 0.0005493.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for
  1/p weighting: two nodes share a special-q with probability <= C(N_MAX,2)
  2^6/2^70.30 <= 2^(2x27.3-1+6-70.30), rounded up to 2^-10.70 <= 0.0006011,
  charged in full. This is a same-draw birthday bound on node special-q draws,
  not a shared-missing-set product.

N_MAX = 2^27.3 is kept from the base (the sum of the last two terms is within 8e-05 of its minimum over N_MAX). No lower bound on
individual node work is used. The boot term treats the 18 relations' subtrees
as independent (H5); the charged 1.000e-03 per boot is 499x the modelled
value. Targets are uniform nonidentity by power randomization (distance <= 1/r per target, 2/r charged). Identity z is handled separately; the final
pairing check e(P, Q) = z rejects any incorrect Q. The success bound is
conditional on H1-H5 (the H2 allowance is a resource-model margin, not a probability term).

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (excess 1.1, per-side bounds, band strip; sample estimate 119.8219 + 1.73% canonical-work allowance) | <= 119.8467 |
| linear algebra (excess 1.1, weight 201, N = 2^59.2453, ten owners) | <= 118.0102 |
| matrix generator | <= 107.3116 |
| both individual logarithms (metered cap) | 113.119 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

The upper-rounded rows sum to 2^120.2136282752853289 (70-digit decimal
arithmetic, audit-checked), so W <= **2^120.213628275286**, the least 12-decimal
value above that sum.

**Configuration.** Relation excess 1.1; factor-base bounds 2^69.25 (f) / 2^69 (g); sub-region theta = 2284/8000 (14987 of 19704 canonical rows),
exact certified count products at radius >= 773, work allowance 1.0173
(2.00 SE), X2 = 2^64, cofactor slacks (295, 226), relation-collection alpha (5.15, 0.85),
descent inputs LOG2_RELS = 63.716, mean norms 737.7/498.4 and bounds 69.25/69, descent alpha
(5.15, 0.85). The 62-node plan is the 90d99567 plan re-optimised by
single-then-joint coordinate search at the sub-bin pr, the skewed boot and the
conservative input (41 entries move, listed in refinement.json; boot window
unchanged), re-evaluated here, evaluated at UDUP = 0.65, THIN_MAX 85, MIN_RELS 18 (nodes),
BOOT_RELS 18; at present inputs min mu 38.5, max phi < 0.05,
max psi 0.483. Two-log expectation <= 2^112.119, meter 2^113.119,
plus the section-7 node counter. Norm-bound, factor-base,
filtering, generator, canonical-list, Schirokauer and verification charges are
the base package's; LA charges (20.7/entry, 420/reduction, 97/dense), blocks
(19, 4) and the 10.01 route are this lineage's; excess 1.1, pad 1.0954,
sub-bin pr, skewed boot, est-2SE input, re-optimised plan, on-region pool, per-side bounds with re-measured filtering and the explicit >= 3-occurrence charge are this revision's; weight 201 with
0.801 N ratio, C_ROW 1,640, recentred kernels and
early abort are the ported charges.
Preprocessing is inside the 2^120.213628275286 cap. No external
advice; nothing amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy with one factor-base bound per side, u_s = (n_s + alpha_s)/lb_s; alpha (5.15, 0.85) against in-package (4.868096, 0.843643) over p <= 2000 and (5.085035, 0.731543) over p < 50000 ([14]), same for descent; exact sampled norms); H1'
retention (Euler-product model per side, calibrated on exact smooth integers); H2 region
(uniform pool at theta 5/16 (164,000 draws, seeds 1-2), uniform sub-region 2284/8000, exact flags,
certified sub-region counts, central-core yield cut (conservative),
and the work mean charged with a 1.73% margin >= 2.00 SE,
not a probability term); H3 matrix (Schirokauer
premises: r unramified [13]; r not dividing h(K_f), h(K_g),
unit condition and one-dimensional kernel assumed, section 5;
signed filtering at excess 1.1 on the per-side system,
pad 1.0954 over the 2^17 maximum, N/U non-increasing 2^12..2^18
and assumed not to rise at production; weight-200 ratio <= 0.801
not rising at production is part of H3, measured 0.775-0.801); H4 kernel cost (KERNELS.md, STORAGE.md),
with H4u: >= 88% of production row weight is +-1 (measured 0.933..0.950 at
weight 100 and 0.932..0.952 at weight 200, `unit_entries.py`;
a shortfall changes cost only); f/g independence given
survival (H1) bounds q_f (charged at 100x);
band-strip ECM finds band primes at the `ecm16` rate;
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence, missing logs,
independent subtrees, fresh-lattice re-descents, local omitted-mass laws,
random-integer ECM orders, UDUP = 0.65 (E[shat/m]/pr_charged
>= 0.87 boot (complete-orbit recount) / >= 1.12 nodes, uniform cells; live [11] at 1.175/0.70); H5n (the per-ideal closed-subsystem
purge-loss rate does not rise from toy to production scale, so the
>= 3-occurrence losses are charged at the pooled toy rates `H3_RATE`); unit-geometry/smoothness independence
for the node orbit haircut).

## 10. Memory

Storage <= 2^105 bytes (STORAGE.md).

## Provenance

Algorithm, models, kernels and package: Gautam Botrel (Claude Code, Codex); polynomials/alpha follow Guillevic-Singh, filtering
CADO-NFS. Refinements: kaankolcu, pepedesigner, mpjunior92, the DeepSeek
lineage (cbab3319, 8c1635f2, e3940d11, 04073dd1, be446c45); grid/plan, both-bars alpha, UDUP 0.65 evidence, 150k pool:
muse-spark-1.3-contributor (v13). Recentred kernels, early abort, weight 200: 2c2ce2ab (gbotrel). Ten owners: 6282ab79 (bernardo-lima12). Line sieve to 2^48
with band strip: 8cc93834 / a31a5f1d and 42e80f80 (sister track). Merge, excess 1.1, plan re-optimisation, band-strip port: Claude
Opus 5.5. v4 (sub-bin pr, skewed boot, est-2SE input, plan re-search, complete-orbit UDUP, purge-loss domination) and v5 (fresh on-region pool, core cut, 2-SE margin, re-measured filtering): the kaankolcu swarm (relief climber Muse Code; exact-orbit recount and skewed-boot measurements by the captain). v6 (this revision: per-side factor-base bounds with filtering re-measured on the per-side system, explicit >= 3-occurrence charge): Claude Fable 5.1 (kaankolcu swarm); the per-side bound follows the sister track's recover-P line (unpublished C4 by Claude Opus 5.5), re-derived and re-measured here. No prior review
verdict is a premise.

## Changes from e3940d11 (121.9944)

v13 (57598408, 121.12966). v1 (640f4e61, 121.038): ported charges of 2c2ce2ab / 6282ab79. v2 (9eae858f, 120.8479): excess
1.1, theta 2581/8000, re-optimised plan. v3: line sieve to 2^48 with band strip, total 2^120.6174.
v4: sub-bin pr, skewed boot, est-2SE input 63.599, plan re-search (41 entries); cap at Markov 0.5. Score 120.521203700786.
v5: fresh 164k pool at theta 5/16, core cut, region 2347/8000 at 2.01 SE; filtering re-measured; descent at 63.589. Score 120.227673100157.
v6 (this revision): factor-base bounds 2^69.25 (f) / 2^69 (g) (bound grid on the same pool, refinement.json `collection_v23`; U 2^63.451 -> 2^63.5788), region 2284/8000 at 2.01 SE,
filtering re-measured on the per-side system, >= 3-occurrence purge losses charged (`H3_RATE`), descent re-evaluated at 63.716. Two-log 112.119, modelled
success 0.4968. Score 120.213628275286.

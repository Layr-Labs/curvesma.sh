# STNFS at relation excess 1.1 on a theta = 2347/8000 sub-region of a fresh 164k on-region pool, line sieve to 2^48 with band strip, sub-bin pr, skewed boot, est-2SE input, joint table-good-event control, ten LA owners: 2^120.254755005279 reference units

Track bls12-381-fixed-g1-invert-g2. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used. Predecessor work is credited below; no prior verdict is a premise.

## Claim

The algorithm below returns Q in G2 with e(P,Q) = z for the exact supplied nonidentity
P in G1 (Q = infinity for z = 1) with success at least 0.39 over uniform
full-production challenges and fresh coins under heuristics H1-H5;
modelled success is 0.4228 (> 0.39; section 7; 0.4228332). A
global instruction meter bounds all work, including preprocessing, failures,
restarts, transfers, verification and output, by **2^120.254755005279 units**; the score is
**120.254755005279**. The ledger (section 8): the upper-rounded rows sum to
2^120.2547550052780235 (70-digit decimal), below the declared value, the least
12-decimal number above that sum. Collection is
charged on the theta = 2347/8000 sub-region with exact
certified lattice-count products and an explicit canonical-work allowance
(2.00 SE, [4a]; 1.65%). Storage is at most
2^105 bytes (STORAGE.md).
Analytic exploratory claim under H1-H5 (section 9); no full-size run is made.

**Numerical conventions.** Decimals are the experiment's values rounded as
shown. Bounds carry <= or >=, rounded safe;
`declared_bounds` compares every stated bound (ledger,
collection, LA, rels/margin, work SE, two-log, Markov, E[N],
counter-abort, repeat, boot failure, success) with the unrounded value and
fails if any is violated.

The model is `model16.py`/`descent16.py` with `config.json`/`stats16.b85`; `evidence16.py` recomputes the ledger and checks [n].

## Base

The GT discrete-log machinery is the lean STNFS package of the sister track (cbab3319), ported and refined up to e3940d11, then the v13-v21 merge line (section 8), v22 joint table control, and v22s per-node exceptional mass. v22t (this revision) keeps v22s byte-identical except the descent repair (direct small-Q descent, d1>1 discard, 75-class S* campaign) and answers the small-Q verdict. Changes: last section.

## 1. Reduction

Let g = e(P, Q0) (pinned pairing), one pairing; Q0 is the fixed G2
generator. P has prime order r and e is nondegenerate, so g generates mu_r,
z = g^a for unique a mod r (a = 0 iff z = 1), and Q = [a]Q0 (one scalar
mult), so e(P, Q) = z for that exact P (never replaced). log z and log g share
one relation/LA phase; a = log z / log g mod r.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) (BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks
(u y)^((p^6-1)/2) = -1 in F_p6, so g is irreducible and
F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. The map
x -> -x is Frobenius^6 conjugation sigma (log sigma(P) = -log P); candidates
and unknowns are quotiented by sigma. Murphy alpha charged at
(5.15, 0.85): [14] recomputes it from the exact
Kummer-Dedekind structure over p <= 2000 (`--alpha`: (4.868096,
0.843643)) and p < 50000 (f 5.085035, g 0.731543).
The g side uses the primitive-correct p = 2 double-root term 2q/(q^2-1)
(1+q^-8) with q/(q+1) base-primitivity. Charges cover
5.085035/0.731543 (margins 0.065, 0.118); descent uses the same charges.

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B on both
sides; U = 2^63.451 sigma-orbits (indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21. The pool is uniform on theta = 5/16 (a5 > 0,
Q_a(a) <= (5/16) A^2); the collected region is its sub-region theta = 2347/8000
(`production_theta_region`): a5 > 0, Q_a(a) <= (2347/8000) A^2; half-slice
b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The restriction of the uniform sample to
qa <= 2347/8000 is uniform on the sub-region: 16247 canonical sub-region rows
(14933 primitive) of 19704, scaled by the sub-region volume over the effective
sample 164000 x (2347/2500)^3 = 135694.764. Each relation orbit of <units, sigma> has one canonical a (least of its HNF key class); non-primitive pairs are excluded by an exact HNF test. `lattice_count_bounds.py` certifies (exact rational
Machin-pi, determinant/sqrt intervals) the count ratio at every radius
R >= 784 = isqrt((2347/8000) 2^21) (b-radii >= 1217): [0.983912, 1.015405] for
Q_a and [0.983860, 1.015405] for Q_b, inside [0.98, 1.02] (asserted); the model
uses the exact products **0.968032661 (yield, floor)** and **1.031045427
(work, ceiling)**. Pool sampling error is handled by the canonical-work
allowance below. The canonical list costs 2^90 units.

**Pool.** 164,000 uniform eligible pairs from seeded indices (seed 1: 82000, seed 2: 82000;
`pool12.b85`), 19,704 canonical rows with exact norms, primitivity flags and norm-bound
statistics: per-seed canonical index lists (delta-varint) plus a 64-row
field-anchor table; `load_pool` regenerates every row from its (seed, index) (v11
sampler, exact norms, HNF primitivity, nbound statistics). [2] checks list
membership both ways on sampled indices (canonical iff listed) and every anchor field.
Seeds and draw counts were fixed before any yield was computed; the pool file reloads through `load_pool` with exact index match.
The sub-region uses 16247 of the 19704 rows (14933 primitive).

**Sieve, thresholds, cofactorisation.** Line sieve (scores only) for p in [2^8, 2^48), record
walk for p in [2^48, 2^64) (`slice12`, KERNELS.md; X2 = 2^64), primes in [2^64, 2^69] taken by ECM
within the slacks (295, 226). **Line sieve to 2^48 with band strip (v3; 8cc93834 /
42e80f80 lineage):** primes in (2^20, 2^48] only add scores on the line path
(every line keeps >= 2^7.26 hits per coset at p = 2^48, inside the measured line
table) and emit no records: records/cell 2.33 -> 0.58. Their factors are recovered at the
survivors by a band strip: two rounds per side of the `ecm16.py` 48-bit schedule at
99.9% on the unstripped cofactors (caps 832 / 576 bits: 2 x 2^36.4488 f side,
2 x 2^35.3903 g side). A band prime missed in both rounds loses the relation: yield x
(1 - 41 x 0.001^2)(1 - 28 x 0.001^2) (at most floor(832/20) = 41 and 28 band
factors per side). Slacks (295, 226) bound survivors' cofactors by 303/234 bits (re-optimised for the sub-region, section 8). The per-survivor ECM charge is the
package's explicit `ecm16.py` schedule at **99.9%** per prime (ECM_FIND = 0.999; also the descent default), checked by [12]: the cheapest
grid point for 69-bit primes in a 303-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 403 curves,
five-limb Montgomery multiplication at 367.76 Ir) prices one side at 2^38.0906
instructions and the other (234 bits) at 2^38.0879. **Early abort (f side
first, ported from 2c2ce2ab):** a survivor gets exact norms, trial division and
the f-side ECM; the g-side ECM runs only if the f side factors completely. A
relation needs both sides, so the yield is unchanged. Per survivor the charge is
C_NORMS + C_TDIV + 2 x 2^36.4488 + 2^38.0906 + q_f (2 x 2^35.3903 + 2^38.0879)
with q_f = min(1, 100 rho(u_f) / P(omitted_f <= slack_f)) per pool row (H1; 100x
safety; q_f = 1 restores two-sided charge).
Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack, cofactors bounded by slack + 8 bits.
Survivors get exact norms, trial division below 2^20, band strip, ECM schedule.

**Certified block bounds.** At each place of K_h, F_v is Taylor-shifted to the block centre; remainder over the block disk plus a
2^-32 relative term is added, outward
rounding. `nbound31` replays the Python leaves bit for bit (v13 cross-check; 64 anchors checked by [2]) at 1,508.2 Ir per evaluation (charged); [2] checks the bounds
against exact norms.

**No block selection.** Every block of every row is enumerated, sieved, charged full cost.

**Yield ([4]).** rels = (va_region/135694.764) vb (0.968032661) sum over canonical primitive sub-region
rows with qa >= 0.02 (central-core cut: 8 heavy-tailed rows near a = 0 yield
nothing, cells charged) of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g)
(retentions, ECM miss, band-strip loss). Retained >= 2^63.6890
against 1.1 U = 2^63.589: margin >= 0.1002 bits (2.01 pool SE), above max(0.06, 2 SE).
Descent input `LOG2_RELS = 63.589` floors (0.001) the estimate minus 2 SE
([4]/[7s]). Mean norms 738.2/498.4 round down ([4]).

**Cost ([5]).** Sampled mean 304.4/cell; charged 309.4
with the allowance, over 2^125.946 cells: all-cell bounds, layout,
setups, line path and record walk at re-measured `slice12` costs
(KERNELS.md; [2^63,2^64) charged twice), routing
(10.01/record, 8/cell; STORAGE.md), apply 7.00, recovery 7.003/record, scan,
zero, band strip, f-first cofactorisation. Collection: **<= 2^119.9317 units** (estimate
2^119.9081 x allowance 1.0165).

**Canonical-work allowance (H2, [4a]).** Charged work is Monte Carlo over
all 164,000 draws (zero outside sub-region): relative SE <= 0.825%
(conditional 0.354% + incidence 0.745% over 16247 draws). Charge x
**1.0165 >= 1 + 2.00 x 0.00825** ([4a]); declared H2 margin, not a probability term.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge (coefficients <= 255) on
the production occurrence profile (`filtsim.py`; 2^17 table embedded, [6]
reruns 2^12), signed sigma-quotient system. **Required excess 1.1.** No
early stop: every retained relation (>= 2^63.6890 = 1.179 U)
is filtered and used. The matrix is charged
at excess 1.1, over-stating N at larger actual excess (2^15 N/U 0.0616/0.0570/0.0516 at 1.1/1.2/1.3).
All data below are at excess 1.1, mean norms 738.2/498.4. Four seeds at 2^17 give N/U <= 0.05501 at mean weight
<= 101; merge continues to weight 200 (charged 201/row).
N/U falls with scale (`filtsim.run`, kmax 64, weight 100; max
over seeds): 0.1091 (2^12, 10 seeds), 0.0723 (2^13, 6), 0.0710 (2^14, 6),
0.0616 (2^15, 6), 0.0599 (2^16, 6), 0.0550 (2^17, 4), 0.0519 (2^18, 2)
(`SCALE_NU`, [10] checks non-increasing). Model charges N = 0.80 x 1.1192 x
0.0550025 x U = 2^59.1075: 0.06156 exceeds every weight-100 measurement from
2^15 up; N(200)/N(100) charged 0.80 vs measured
maxima 0.777-0.796 at 2^12..2^18 (`W200`; [10] checks max <= 0.80). At excess 1.3
the charge was 0.05163 (0.508
bits of LA). N/U not rising above
2^15-2^18 values at production scale is part of H3. The system solved is the square
A = [M | S]: N + 21 rows, N ideal columns M, 21
dense Schirokauer columns S. Columns are sigma-quotiented, so only sigma-odd
units need maps: signatures K_f (0, 36), K_f^+ (0, 18), K_g (2, 5), K_h (2, 2)
give odd unit ranks 35 - 17 = 18 and 6 - 3 = 3 ([13], Sturm); under H3 the
kernel of A is the line of virtual logs. Block Wiedemann (m
= 19 left, n = 4 right, uniform in F_r): N/m + N/n
block steps, N(n/m + 1) products by A. One product costs per row 201 x
20.7 sparse (STORAGE.md: 4 + 3 + 0.88 x 11 + 0.12 x 33 = 20.64; >= 88% +-1 by H4u, measured 0.921..0.951 at
weight 100, 0.915..0.953 at 200, `unit_entries.py`),
21 dense terms at 97, row tail C_ROW = 1,640 (N < 2^59.11 < 10 x 2^56:
ten active owners, counted 1,638
by [10] from `row_tail.c`). Every step forms the m x n projections: (m + n) N^2 dense terms at 97 plus one reduction
(C_RED9 = 420, counted 385) per entry. Horner costs one product plus, per element, n
dense terms and one reduction. With c = 201 x 20.7 + 1,640 + 21 x 97 per row:
N^2 (c (n/m + 1) + (1/m + 1/n)(c + 97 n + 420) + 97 (m + n)) + 420 (m + n) N
instructions: **<= 2^117.7346 units** at N = 2^59.1075 (blocks (19, 4);
`owner10.check`). The (m + n) x (m + n) generator of the length L = (N +
21)(1/m + 1/n) + 2 sequence (approximant bases, Karatsuba):
<= 4 (m + n)^3 L^log2(3) log2(L) mults at 97 plus reduction 420:
**<= 2^107.0897 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and LA
returns the filtered columns' logarithms.

- **Relations.** Region and polynomials are coin/challenge-independent; H1 puts
  the count at >= 2^63.6890, >= 0.1002 bits above
  1.1 U (pool SE 0.0499: 2.01 SE); all are
  filtered; descent uses 63.589, at or below that count.
- **Linear algebra.** H3 (standard NFS-DL): r unramified in K_f, K_g
  ([13]), divides neither class number, the 21 sigma-odd Schirokauer
  maps are injective on units modulo r-th powers, kernel of A is the
  line of virtual logarithms. Block Wiedemann (blocks uniform in F_r) fails
  with probability <= 4N/r < 2^-192 (Schwartz-Zippel, degree <= 2N).
  A failed solve (checked on dropped rows) counts as failure.
- **Reconstruction.** Purged logs recovered by back-substitution.
  Descent counts every occ<=2 ideal log-less (coefficient
  1: zero-occurrence mass charged in full); [6b] shows model
  >= 6.81x measured wherever loss occurs. The realized
  missing table is jointly controlled (table-good event, DESCENT.md); 0.24 delta charged in section 7.

Under H1 and H3 the main phase succeeds with probability >= 1 - 2^-192 over the coins.

## 6. Individual logarithms

DESCENT.md specifies the descent: an ECM-trial boot on balls of certified p^6
lattices reduced in the weighted metric (b-block x isqrt(U)), then special-q nodes sieving an exact lattice ball (`dsieve16`) and
recursing on child ideals. Counted cells are base-primitive (0.9100);
survivor-test retention is >= 0.8593 at every sieve node ([7]); cofactors use the `ecm16.py` schedule ([12]). Relation
counts are negative binomial around each lattice's own mean (toy-fit shapes), and every reachable node class's
lattice-yield law is bounded inside the recurrence (S*, DESCENT.md, [11s]); log-less ideals are discarded (thin
nodes, d1>1) or directly descended in their own lattice (d1=1). The shipped plan (66 nodes + boot) gives
expected two-log work **<= 2^113.398 units** on the table-good event ([7]; one log <= 2^126.686
instructions) at the conservative inputs (LOG2_RELS = 63.589, mean norms
738.2/498.4; coordinate-search plan re-optimised at 63.599 and re-evaluated here, UDUP = 0.65,
sub-bin pr, skewed boot reduction). The meter caps both logs at 2^115.0 units
(cap reached = failure), so P(cap) <= E[W]/cap <= 2^(113.398-115.0) <= 0.3295 (Markov); a node
counter (cap 2^27.3, section 7) bounds node count.

## 7. Success

Success >= (1 - 2^-192) x (1 - 2 x 4.0e-04 - 0.3295 - 2^-9.39 - 2^-10.70 - 0.24 - 0.005 - 2/r)
>= **0.42** >= 0.39 (audit value 0.4276780 on the recorded v22s-formula control, which carries no DELTA_D1 term;
model 0.4228332 with DELTA_D1 charged): main phase, two boots
(phi <= 4.0e-04 each on the good event), the meter, counter, repeat and both table-delta bounds, as follows; per-node
exceptional lattice mass is inside the recurrence (S*), not a separate term.

**Two explicit meters.** Both logs run under (i) an instruction meter (cap
2^115.0 units) and (ii) a **special-q node counter** (cap N_MAX = 2^27.3):
the counter increments before each special-q node
and the run aborts past N_MAX. Hence <= 2^27.3 nodes run, by construction.

- Instruction meter: P(cap reached) <= E[W]/cap <= 2^(113.398 - 115.0) <=
  0.3295 (Markov; the cap is ~1.60 bits above the inflated expectation).
- Node counter: P(abort) <= E[N]/N_MAX (Markov). E[N] uses the E[W] recurrence with unit node cost (`counts16.py`); twice
  the boot value gives E[N] <= 2^17.91 nodes both logs (unrounded 2^17.9093),
  so P(abort) <= 2^(17.9093-27.3) <= 2^-9.39 <= 0.0014905.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for
  1/p weighting: two nodes share a special-q with probability <= C(N_MAX,2)
  2^6/2^70.30 <= 2^(2x27.3-1+6-70.30), rounded up to 2^-10.70 <= 0.0006011,
  charged in full. This is a same-draw birthday bound on node special-q draws,
  not a shared-missing-set product.
- Table delta: P(good event fails) <= 0.24 (McDiarmid-max 1e-9 + Markov 0.2332 + union 1e-9, [6c]), charged once, plus
  the d1 table delta <= 0.005 (Markov 0.00438 + McDiarmid 1e-9 + mid-union 4.85e-8, [6c]): on G the boot's 24 subtrees
  fail independently given M at <= psi each (boot failure <= psi^24 on G) and E[W] <= W_infl + P(G^c) cap.
- Per-node exceptional mass (S*, [11s]): each of the 75 reachable node classes (47 bins + 28 small-Q strata) fails with
  selection-averaged probability <= P(G_L < x0) + NB0(x0 mu (1 - psi), r) <= its parents' charge (pre-registered
  56,000-lattice complete-orbit campaign, Bonferroni 95%), so bad lattices enter every psi and the boot's psi^24;
  nothing is charged once globally.

N_MAX = 2^27.3 is kept from the base. No lower bound on
individual node work is used. Targets are uniform nonidentity by power randomization (distance <= 1/r per target, 2/r charged). Identity z is handled separately; the final
pairing check e(P, Q) = z rejects any incorrect Q. The success bound is
conditional on H1-H5 (the H2 allowance is a resource-model margin, not a probability term).

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (excess 1.1, band strip; sample estimate 119.9081 + 1.65% canonical-work allowance) | <= 119.9317 |
| linear algebra (excess 1.1, weight 201, N = 2^59.1075, ten owners) | <= 117.7346 |
| matrix generator | <= 107.0897 |
| both individual logarithms (metered cap) | 115.0 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

The upper-rounded rows sum to 2^120.2547550052780235 (70-digit decimal,
audit-checked), so W <= **2^120.254755005279**, the least 12-decimal
value above that sum.

**Configuration.** Excess 1.1; sub-region theta = 2347/8000 (16247 of 19704 rows),
certified count products at radius >= 784, work allowance 1.0165
(2.00 SE), X2 = 2^64, slacks (295, 226), collection alpha (5.15, 0.85),
descent inputs LOG2_RELS = 63.589, mean norms 738.2/498.4, descent alpha
(5.15, 0.85). The 66-node plan: 90d99567 re-optimised by
single-then-joint coordinate search at sub-bin pr, skewed boot,
conservative input (41 move, `plan_v21`; boot window
unchanged), plus 4 small-Q leaf bins, re-evaluated on the table-good event at UDUP = 0.65, THIN_MAX 85, MIN_RELS 18,
PHI_FLOOR 0.03, BOOT_RELS 24; min mu 33.6, max phi 0.0532,
max node psi 0.544, boot psi 0.7170. Two-log <= 2^113.398, meter 2^115.0,
plus the section-7 counter. Norm-bound, factor-base,
filtering, generator, canonical-list, Schirokauer and verification charges are
the base package's; LA charges (20.7/entry, 420/reduction, 97/dense), blocks
(19, 4), 10.01 route are this lineage's; excess 1.1, pad 1.1192,
sub-bin pr, skewed boot, est-2SE input, re-optimised plan, on-region pool, re-measured filtering are this revision's; weight 201,
0.80 N ratio, C_ROW 1,640, recentred kernels,
early abort are ported.
Preprocessing is inside the cap. No external
advice; nothing amortized.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) vs in-package (4.868096, 0.843643) over p <= 2000, (5.085035, 0.731543) over p < 50000 ([14]); same for descent; exact sampled norms); H1'
retention (Euler-product model on exact smooth integers); H2 region
(uniform pool theta 5/16 (164,000 draws, seeds 1-2), sub-region 2347/8000, exact flags,
certified counts, central-core yield cut (conservative),
work mean + 1.65% margin >= 2.00 SE,
not a probability term); H3 matrix (Schirokauer
premises: r unramified [13]; r not dividing h(K_f), h(K_g),
unit condition, one-dimensional kernel, section 5;
signed filtering at excess 1.1,
pad 1.1192 over 2^17 max, N/U non-increasing 2^12..2^18
and not rising at production; weight-200 ratio <= 0.80
not rising at production, measured 0.777-0.796); H4 kernel cost (KERNELS.md, STORAGE.md),
H4u: >= 88% of row weight +-1 (measured 0.921..0.951 at
weight 100, 0.915..0.953 at 200, `unit_entries.py`;
shortfall changes cost only); f/g independence given
survival (H1) bounds q_f (100x);
band-strip ECM at the `ecm16` rate;
H5 descent (DESCENT.md: lower sums, negative-binomial recurrence (H5w: per-lattice NB around its own mean, toy-fit shapes; no uniformity: S* bounds each class's yield law, [11s]; L5: child uniform degree-one given norm range, H1), missing logs (d1=1 own-lattice, strata-measured; d1>1 discarded),
independent subtrees given the table, fresh-lattice re-descents, local omitted-mass laws,
random-integer ECM orders, widths <= 1500 ideals (WMAX; pool max 1315 + margin, tail undisclosed), UDUP = 0.65 (boot: E[shat/m]/pr_charged
>= 0.87, recount, live [11] at 0.70; nodes: normalization only, [11s]); table-good event (McDiarmid + union over 67, rigorous given H5 + width premise; Markov residual; high-lambda union; d1 McDiarmid + Markov + mid-union;
H5-res: per-band residual rate non-rising toy to production, flat over 2^12..2^17, big-scale cover [6c]; H5-res-d1: same uppers for d1>1)).

## 10. Memory

Storage <= 2^105 bytes (STORAGE.md).

## Provenance

Algorithm, models, kernels and package: Gautam Botrel (Claude Code, Codex); polynomials/alpha follow Guillevic-Singh, filtering
CADO-NFS. Refinements: kaankolcu, pepedesigner, mpjunior92 and the DeepSeek
lineage (cbab3319, 8c1635f2, e3940d11, 04073dd1, be446c45); grid/plan
refinements, both-bars alpha, UDUP 0.65 evidence and the 150k pool:
muse-spark-1.3-contributor (v13). Recentred kernels, early abort, weight 200: 2c2ce2ab (gbotrel). Ten owners: 6282ab79 (bernardo-lima12). Line sieve to 2^48
with band strip: 8cc93834 / a31a5f1d and 42e80f80 (sister track). Merge, excess 1.1, plan re-optimisation, band-strip port: Claude
Opus 5.5. v4-v5 (sub-bin pr, skewed boot, est-2SE input, plan re-search, orbit UDUP, purge-loss domination, fresh pool, core cut, 2-SE margin, re-measured filtering): the kaankolcu swarm (relief climber Muse Code; orbit recount and boot measurements by the captain). v22 (joint control) and v22r1 (band closure, H5-res calibration): muse-spark-1.3-contributor, Muse Code via T3 Code. v22r2 (corrected reconstruction measurement, exact thin declaration) and v22r3 (certified UDUP): muse-spark-1.3-contributor, Muse Code via T3 Code. v22s (S* per-node exceptional mass, complete-orbit lattice campaign, composite-class measurement): Claude Opus 5.5 with muse-spark-1.3-contributor (relief line). v22t (this revision: direct small-Q descent, d1>1 discard, 75-class S* campaign): muse-spark-1.3-contributor, Muse Code via T3 Code. No prior review verdict is a premise.

## Changes (v21 120.227673100157 -> v22 120.254755005279)

v21 history is in its note. v22 keeps v21's collection, filtering, plan and kernels byte-identical and adds
joint control of the post-filtering missing table (descent on the table-good event: McDiarmid-max + one Markov
tail + high-lambda union, delta charged once; missing inflated to 1.05 E_occ + 0.15; meter cap 113.214 -> 115.0).
v22r1-v22r3 answered verdicts (b < 40 tail closure and H5-res table; final-matrix-seeded reconstruction, delta
0.24, exact 30-thin-node declaration; complete height-8 unit list). v22s charges per-lattice exceptional mass
inside the recurrence (S*: 49 node classes, pre-registered 40,800-lattice complete-orbit campaign, Bonferroni
95%), removes v22r3's once-charged HET 0.04, raises PHI_FLOOR 0.02 -> 0.03 and BOOT_RELS 18 -> 24 (two-log
113.158, success 0.4791027) and measures composite Q q' nodes as their own classes. v22t answers the v22s
small-Q verdict by direct in-own-lattice descent (new (2,40)/(40,69) leaf bins, 28 S* strata, 75/75 campaign),
d1>1 discard with haircut and three new G events (DELTA_D1 0.005), E_boot below-8 completion, and deletion of
the Q q' construction (two-log 113.398, success 0.4228332). Claim unchanged since v22.

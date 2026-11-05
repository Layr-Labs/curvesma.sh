# STNFS v9: exact static LA charges on v8: box-shared norm bounds on v7 (sub-bin yields, re-searched lazy descent, shared Krylov overhead): 2^119.932465744 reference units

Track bls12-381-fixed-g1-invert-g2. One reference unit is a dense F_p12
multiplication at C = 20,000 scalar AArch64 instructions. No auxiliary powers,
external advice, smaller targets or future-target amortization are used.
Predecessor work is reused as credited below; no prior verdict is a premise.

## Claim

The algorithm below returns Q in G2 with e(P,Q) = z for the exact supplied
nonidentity P in G1 (Q = infinity for z = 1) with success at least 0.39 over
uniform full-production challenges and fresh coins under the declared heuristics
H1-H5 (model 0.3902160; audit with rounded inputs 0.3902004). A global
instruction meter bounds all work, including preprocessing, failures, restarts,
transfers, verification and output, by **2^119.932465744 units**; the score is
**119.932465744**. Upper-rounded ledger rows (section 8) sum to
2^119.932465743989235 <= 2^119.932465744. Storage is at most 2^105 bytes (STORAGE.md). The claim is
analytic and exploratory; no full-size run is made.

**v9 refinement.** Same v8 algorithm, kernels, pool, relation count, descent,
heuristics and storage. Four rounded LA prices tighten: Horner sparse entries
20.640000001 (mixture20.64), Krylov entries12.896000001 (mixture12.896),
row tail1424, projection term59. The 1e-9 margins cover floating comparisons;
written decimals are rational bounds. No kernel or heuristic fraction changes.
The reduction charge remains420: the counted385 covers the reduction kernel,
whereas420 also reserves the Horner modular add and outer projection-dot
setup/teardown (the59 charge covers only inner-loop terms). Advisory review identified
this distinction and a proposed cut to385 was withdrawn before submission.

For (m,n)=(28,5), the quadratic coefficient changes from13141.59 to
13094.428286139985715 (ceiling); the linear coefficient remains13860.
Since N>2^58, the old linear/quadratic fraction is below4e-18.
Multiply the old certified LA bound by the quadratic ratio plus1e-15,
which safely covers the unchanged linear term. This yields log2 LA
<=117.274656111. Generator stays<=107.806519571. All other rows stay.
High-precision sum of the upward-rounded ledger rows is log2
119.932465743989235, below the declared119.932465744. This is a tighter
accounting bound, not a faster kernel or novel cryptanalysis. Setup,
inter-row control and all other reserves remain charged as in v8.

**Numerical conventions.** Decimals are the experiment's values rounded as
shown; bounds are rounded in the safe direction and `declared_bounds` compares
every stated bound with the unrounded value. The model is `model16.py` and
`descent16.py` with `config.json` and `stats16.b85`; `evidence16.py` recomputes
the ledger and the checks cited as [n].

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
(5.15, 0.85) in collection and descent; [14] recomputes it from the exact
Kummer-Dedekind structure: (4.868096, 0.843643) over p <= 2000, (5.085035,
0.731543) over p < 50000 (g side with the primitive-correct p = 2 double-root
term and q/(q+1) factor); margins 0.065, 0.118.

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B on both
sides; U = 2^63.451 sigma-orbits (indices fit 64 bits).

**Region and fundamental domain.** Two exact rational forms Q_a, Q_b
(`config.json`); A^2 = 2^21. The pool is uniform on theta = 5/16 (a5 > 0,
Q_a(a) <= (5/16) A^2); the collected region is its sub-region theta = 2330/8000
(`production_theta_region`): a5 > 0, Q_a(a) <= (2330/8000) A^2; half-slice
b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The restriction of the uniform sample to
qa <= 2330/8000 is uniform on the sub-region: 16033 canonical sub-region rows
(14709 primitive) of 19765, scaled by the sub-region volume over the effective
sample 164000 x (2330/2500)^3 = 132767.441. The canonical a of each unit class is
the least of its exact HNF key class in (Q_a, lex) order, so each relation orbit
of <units, sigma> has one point in the domain; non-primitive pairs are excluded
by an exact HNF test. `lattice_count_bounds.py` certifies (exact rational
Machin-pi, determinant/sqrt intervals) the count ratio at every radius
R >= 781 = isqrt((2330/8000) 2^21) (b-radii >= 1219): [0.983851, 1.015464] for
Q_a and [0.983798, 1.015464] for Q_b, inside [0.98, 1.02] (asserted); the model
uses the exact products **0.967911517 (yield, floor)** and **1.031166365
(work, ceiling)**. Pool sampling error is handled by the canonical-work
allowance below. The canonical list costs 2^90 units.

**Pool.** 164,000 uniform eligible pairs on theta = 5/16 from seeded indices
(seed 5: 2,000, seed 6: 162,000; `pool12.b85`), 19,765 canonical rows (12.05%);
`load_pool` regenerates every row from its (seed, index) (exact norms, HNF
primitivity, certified nbound statistics); [2] checks list membership both ways
and a 64-row anchor table. Every draw was tested (no selection; one draw needs a
10^6-node orbit budget); seeds and draw count were fixed before any yield was
computed.

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
package's explicit `ecm16.py` schedule at **99.9%** per prime (ECM_FIND = 0.999), which [12] checks: the cheapest
grid point for 69-bit primes in a 303-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 403 curves,
five-limb Montgomery multiplication at 367.76 Ir, `ecm12`) prices one side at 2^38.0906
instructions and the other (234 bits) at 2^38.0879. **Early abort (2c2ce2ab):** the g-side band strip and ECM run
only after the f side factors completely (yield unchanged). Per survivor: C_NORMS + C_TDIV + 2 x 2^36.4488 + 2^38.0906 +
q_f (2 x 2^35.3903 + 2^38.0879), q_f = min(1, 100 rho(u_f)/P(omitted_f <= slack_f))
per pool row (heuristic under H1, mean ~1e-5, factor 100 a margin). The 99.9% is the same
ECM_FIND used in the yield's miss model.
Certified per-block bounds (`embedding_bounds.py`,
64-cell blocks of width 128, bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell with
omitted mass below slack and bounding survivors' cofactors by slack + 8 bits.
Survivors get exact norms, trial division below 2^20, the band strip and the ECM
schedule above (f side first, g side only after a complete f factorisation).

**Certified block bounds.** F_v Taylor-shifted to the block centre per place,
remainder over the block disk plus a 2^-32 relative term, rounded outward; one
evaluation (`NormBounds.block`) is what `nbound31` replays bit for bit (v13 pool,
17,874 rows) at 1,508.2 Ir (charged as measured). **Box bounds (v8).** The rows of
one line (same a, b1, b5) whose (b2, b3, b4) lie in one aligned 2 x 2 x 4 box
(`model16.box_rows`, up to 16 rows) share one bisection per 64-cell block.
`EB.hull` gives per place the centre of the rows' context centres and a radius
>= |c_k - centre| + r_k for every row k, each step rounded outward: row k's b(v)
lies within r_k + r0 of c_k + m, hence within the hull radius + r0 of centre + m,
so one `block` on the hull bounds every row's cells. `EB.box_leaves` bisects the
union of the rows' b0 ranges until both bounds are within 8 bits; a single column
the hull leaves unresolved is evaluated with each overlapping row's own context,
and a row's cell is rejected only if that fails (no pool row is). Boxes x blocks
x (b1, b5) partition each slice (every box row keeps b5 >= 1, the sampled half),
so charging each cell C_NB_EVAL x its box's evaluations per box cell has
expectation total evaluations / cells over the uniform pool; each row is charged
every leaf it uses (C_SCAN_LEAF x leaves / blen) and its context plus a hull share
(C_NB_ROW 4,300 -> 4,600; the hull update is ~25 instructions per place and row).
Sub-region pool means: evaluations per cell 0.0409 (per row, v7) -> 0.0076,
leaves per cell 0.0286 -> 0.0335. [2b] rebuilds six boxes and checks 1,536 exact
norms at sampled cells of every box row against their leaf bounds; [2] checks the
box fields of the 64 anchors. Every block of every row is enumerated, sieved and
charged (no selection).

**Yield ([4]).** rels = (va_region/132767.441) vb (0.967911517) sum over canonical primitive sub-region
rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g), with Dickman rho
([1]), retentions r from the finite Euler-product threshold model
(`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14)
the ECM miss and the band-strip loss. **Central core (v4):** relations from
cells with qa = Q_a(a)/A^2 < Q_CORE = 0.02 are discarded (never passed to
filtering, linear algebra or the descent; those cells are still sieved and
charged), so count and relation-weighted mean norms are those the estimator
measures over the rows with qa >= 0.02; the discard only lowers the count. Retained relations >= 2^63.6617
against 1.1 U = 2^63.589 (excess 1.1, section 4): margin >= 0.0728 bits (unrounded
0.07289; pool SE 0.0364, so 2.00 SE). Feasibility requires a margin of at least
max(`MARGIN` = 0.06, `MARGIN_SE` = 2 SE). The descent's missing-log input `LOG2_RELS =
63.661` is the floor (to 0.001 bit) of this retained count (unrounded
2^63.66172; certified lower product already applied). [4] checks `LOG2_RELS <=
log2_rels` (`rels_input` check). The mean-norm inputs 739.9/498.9 are rounded
down from 739.98/498.95 and checked in [4].

**Cost ([5]).** Sampled mean 255.7 instructions/cell (v7 306.3); charged mean 260.2
with the allowance below, over 2^125.934 cells: all-cell bounds, layout,
setups, line path and record walk at the re-measured recentred `slice12` costs
(KERNELS.md; primes in [2^63, 2^64) charged twice), routing
(10.01/record, 8/cell; STORAGE.md), apply 7.00, recovery 7.003/record, scan,
zero, band strip and f-first cofactorisation. Collection: **<= 2^119.669318937 units** (estimate
2^119.6444 times allowance 1.0174).

**Canonical-work allowance (H2, [4a]).** The charged work is (va_region/132767.441)
times the per-row sum over in-region canonical rows: a Monte Carlo estimate of
a fixed, challenge-independent domain total. The estimator sums all 164,000
draws (zero outside the sub-region or non-canonical), so its relative SE has a
conditional row-mean term and an incidence term for the 16033 contributing
draws: conditional 0.434%, incidence sqrt((1-16033/164000)/16033) = 0.750%, total
<= 0.867%. The charge multiplies the estimate by **1.0174 >= 1 + 2.00 x 0.00867
= 1.01734** (rounded up; unrounded need 1.0173370); [4a] fails if the allowance falls below 1 + 2.00 SE.
This is a declared modeling margin under H2, not a probability term, and no
sampling tail is charged against success.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge (coefficients <= 255, including a
collision-safe bound on implicit middle-range entries), simulated on the
production occurrence profile (`filtsim.py`; [6] reruns 2^12) as the signed
sigma-quotient system. **Required excess 1.1.** The whole region is sieved and
every retained relation (>= 2^63.6617 = 1.157 U) is filtered and used; the
matrix is charged at excess 1.1, which over-states N at the larger actual
excess. Data are measured at excess 1.1 with the mean norms 739.9/498.9. Four
seeds at 2^17 give N/U <= 0.05797 at weight 100; N/U falls with scale (max over
seeds): 0.1174 (2^12), 0.0748, 0.0675, 0.0618, 0.0592, 0.0580, 0.0529 (2^18)
(`SCALE_NU`, [10] checks non-increasing); pad 1.066 covers every value from 2^15
up. **Merge weight 300 (v6).** The merge continues to mean weight 300 (charged
301 per row); N(300)/N(100) maxima per scale are 0.7098 (2^12, 6 seeds), 0.7034
(2^13, 6), 0.7043 (2^14, 3), 0.6971 (2^15, 3), 0.6997 (2^16, 3), 0.6897 (2^17,
2), 0.6913 (2^18, 1), charged **0.711** (`W300`, [10] checks max <= 0.711); the
weight-200 ratio was 0.776-0.795 (v5 charge 0.80). The coefficient cap never
binds (largest explicit or implicit bound 192 < 255 at weight 300). N = 0.711 x
1.066 x 0.05797 U = **2^58.9429** < 8 x 2^56. N/U and the ratio not rising at
production scale are part of H3.

The system is the square A = [M | S]: N + 21 rows, N ideal columns and 21 dense
Schirokauer columns (sigma-odd unit ranks 18 + 3 from signatures K_f (0, 36),
K_f^+ (0, 18), K_g (2, 5), K_h (2, 2), [13]); under H3 its kernel is the line of
virtual logs. Block Wiedemann: m = 28 left vectors X with entries uniform in
[0, 2^128), n = 5 right vectors Y uniform in F_r. The Krylov phase takes N/m +
N/n block steps, N(n/m + 1) vector-products.

- **Krylov products (STORAGE.md).** The five sequences advance in lockstep;
  element j of the five current vectors shares one 256-byte slot, and one pass
  over a row serves all five: `unit_run5` handles a +-1 entry for five vectors
  in 39 instructions, `gen_run1` a coefficient entry in 33 per vector (clang 18
  listings, counted by [10]). The inherited 4 + 3 (block selection, sign) is
  charged once per +-1 entry (one pass for five vectors) and per vector for
  other entries: per vector-entry 0.88 x (4 + 3 + 39)/5 + 0.12 x (4 + 3 + 33) =
  12.896 -> **12.896000001** (>= 88% +-1 by H4u, measured 0.905-0.948 at weight 300,
  2^12..2^15). Per row and vector: 301 x 12.896000001 + C_ROW + 21 x 97.
- **Row tail.** Eight active 2^56-element owners: seven sends and seven
  combinations, `_finish`, headers, routing: counted 1,424 by [10]
  (`owner10.row_tail`), charged C_ROW = 1,424.
- **Horner and projections.** Horner (w <- A w + Y F_i, deg F <= N/m + N/n) uses
  single-vector products at 301 x 20.640000001 + 1,424 + 21 x 97 plus n dense terms and
  one reduction (420) per element. Projections: (m + n) N^2 terms of a 128-bit X
  entry times a field element at 59 (`bw_dot2`, counted 59) plus a reduction
  per entry.

With c_K and c the Krylov and single-vector row costs: N^2 (c_K (n/m + 1) +
(1/m + 1/n)(c + 97 n + 420) + 59 (m + n)) + 420 (m + n) N instructions =
**<= 2^117.274656111 units** at N = 2^58.9429 ((28, 5) minimises it at n = 5). The
(m + n) x (m + n) generator of the length L = (N + 21)(1/m + 1/n) + 2 sequence:
<= 4 (m + n)^3 L^log2(3) log2(L) field mults at 97 + 420: **<= 2^107.806519571 units**.

## 5. Main-phase success

The main phase succeeds when the relations reach the excess and linear algebra
returns the logarithms of the filtered columns.

- **Relations.** Region and polynomials are coin/challenge-independent, so the
  count is one fixed number; H1 puts it at >= 2^63.6617, >= 0.0728 bits above
  1.1 U (pool SE 0.0364, so 2.00 SE above requirement); all of them are
  filtered and the descent uses 63.661, at or below that count.
- **Linear algebra.** H3 (standard NFS-DL): r is unramified in K_f and K_g
  (checked by [13]), divides neither class number, the 21 sigma-odd Schirokauer
  maps are injective on the units modulo r-th powers, and the kernel of A is the
  line of virtual logarithms. Block Wiedemann then fails with probability
  <= 4N/2^128 < 2^-66 (Schwartz-Zippel on the generator's determinant conditions,
  degree <= 2N in the block entries, each entry uniform on a set of size >= 2^128).
  A failed solve (checked on the dropped rows) counts as failure.
- **Reconstruction.** Purged ideals' logs are recovered by back-substitution.
  The descent counts every ideal with <= 2 occurrences as log-less; [6b]
  measures ideals left without a log on the [6] system, and the model's count
  is >= 12.5x measured in every band with loss (DESCENT.md).

Under H1 and H3 the main phase succeeds with probability >= 1 - 2^-66 over the coins.

## 6. Individual logarithms

DESCENT.md specifies the descent: an ECM-trial boot on balls of certified p^6
lattices, then special-q nodes sieving an exact lattice ball (`dsieve16`) and
recursing on child ideals. Counted cells are base-primitive (0.9100);
survivor-test retention is >= 0.8592 at every sieve node ([7]); cofactors use
the `ecm16.py` schedule ([12]). Relation counts are negative binomial (toy-fit
shapes). Log-less ideals are discarded or re-descended in a fresh lattice. The
plan (62 nodes + boot (215, 328, 16, 4)) is re-searched at the v7 yields.

**Lazy descent (v6; DESCENT.md).** The first relation whose subtree succeeds
ends a node, so per-relation work is done on demand. (i) After the full sieve, a
node cofactors its s survivors in uniformly random order, descending each distinct
relation at once: E[processed]/s <= E_M min(1, (1 + 1/s)/((1 - psi)(M + 1)) +
psi^M), M ~ NB(mu, r) as in the failure model (`lazy_fraction`; [7l] simulates the
process at the three costliest nodes: 0.0228/0.0160/0.0327 vs bounds
0.0234/0.0166/0.0360). (ii) The boot finds and descends one relation at a time,
at most 18: by Wald (1 - psi^18)/(1 - psi) = 1.983 expected relations, 2^56.57
trials. Failure events are those of the recurrence (DESCENT.md).

Expected two-log work **<= 2^112.2596 units** ([7]; one log 2^125.548
instructions: boot own about 2^123.58 + subtree about 2^124.13 x 1.983). The
meter caps both logs at **2^112.975291 units** (cap reached = failure): P(cap) <=
E[W]/cap <= 2^(112.2596 - 112.975291) <= 0.6089135 (Markov; unrounded 0.6088978).

## 7. Success

Success >= (1 - 2^-66) x (1 - 2 x 3.258e-06 - 0.6089135 - 2^-11.81 - 2^-10.70 - 2/r)
>= **0.3902** >= 0.39 (audit 0.3902004; model 0.3902160): main phase, two boots
(phi <= 3.258e-06 each), the instruction-meter Markov bound, the node-counter
Markov bound and special-q repeats.

- Instruction meter: cap 2^112.975291 units, P(cap) <= 0.6089135 (section 6).
- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX
  = 2^27.3. E[N] uses the same recurrence with unit node cost (`counts16.py`):
  E[N] <= 2^15.49 for both logs, P(abort) <= 2^-11.81 <= 0.0002786.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for
  1/p weighting: P <= 2^(2 x 27.3 - 1 + 6 - 70.30) <= 2^-10.70, charged.

No lower bound on node work is used. Subtrees are treated as independent (H5).
Targets are uniform nonidentity by power randomisation (2/r charged); z = 1 is
handled directly; the final check e(P, Q) = z rejects any incorrect Q.

## 8. Ledger

| component | log2 units |
|---|---:|
| relation collection (excess 1.1, band strip, box bounds, 164k pool; estimate 119.6444 + 1.74% canonical-work allowance) | <= 119.669318937 |
| linear algebra (weight 301, N = 2^58.9429, eight owners, five-vector Krylov passes, blocks (28, 5)) | <= 117.274656111 |
| matrix generator | <= 107.806519571 |
| both individual logarithms (metered cap; expected 112.2596) | 112.975291 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| filtering streams, misc | 80 |
| verification and output (pairing, G2 scalar multiplication, check) | 11 |

**Fixed rows (upper bounds; units x 2^14.29 = instructions).** Canonical list
2^90: 2^41.8 instructions per eligible a (2^62.5 of them) against one LLL-reduced
HNF key, the `canonical_a` orbit enumeration and a segmented key sort (< 2^30).
Schirokauer maps 2^88: 2^38.5 per relation (<= 2^63.8) against one exponentiation
per side mod r^2 (< 2^34.5). Factor base 2^72: 2^23.4 per prime p < 2^69 against
root finding (< 2^19). Filtering, purge, merge bookkeeping and back-substitution
2^80: 2^24.4 per relation entry (2^69.9) against < 2^14. Certificates 2^50;
pairing, [a]Q0, check 2^11.

The upper-rounded rows sum to 2^119.932465743989235, so W <=
**2^119.932465744** (audit-checked).

**Configuration.** Excess 1.1; pool theta 5/16 (164,000 draws); region theta
2330/8000, Q_CORE 0.02, radius 781, allowance 1.0174, X2 = 2^64, slacks (295,
226), alpha (5.15, 0.85); descent inputs 63.661, 739.9/498.9; sub-bin yields
(caps 4x/1.5x), plan and boot re-searched (v7), UDUP 0.65, THIN_MAX 85,
BOOT_RELS 18, lazy nodes and boot; meter 2^112.975291, node counter 2^27.3. LA:
weight 301, ratio 0.711, pad 1.066, C_ROW 1,424 (eight owners), 12.896000001 (Krylov) /
20.640000001 (Horner) per entry, 97 per dense term, 420 per reduction, blocks (28, 5)
with 128-bit X entries.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) in collection and descent against the computed values of section 2 ([14]); exact sampled norms); H1'
retention (Euler-product model calibrated on exact smooth integers); H2 region
(uniform pool and sub-region, yield lower bound without the qa < 0.02 core, exact flags,
certified counts, and the canonical work mean charged with the 1.74% allowance,
a modeling margin, not a probability term); H3 matrix (Schirokauer
premises: r unramified checked by [13]; r not dividing h(K_f) or h(K_g),
the unit condition and the one-dimensional kernel assumed, section 5;
signed filtering on the production occurrence profile at relation excess 1.1,
pad 1.066 over the 2^17 maximum, N/U (max over seeds) non-increasing 2^12..2^18
and assumed not to rise at production; the weight-300 ratio N(300)/N(100) <= 0.711
not rising at production scale is part of H3, measured 0.690-0.710 at
2^12..2^18 at excess 1.1); H4 kernel cost (KERNELS.md, STORAGE.md; LA tail/reductions are
statically counted loop-free code, the route an executed-instruction count),
with H4u: >= 88% of production row weight is +-1 (measured 0.905..0.949 at
weights 100-300, 2^12..2^15, `unit_entries.py`; a shortfall changes cost only); f/g independence given
survival (H1) bounds the early-abort probability q_f (charged at 100x the
estimate); band-strip ECM finds band primes at the `ecm16` per-prime rate (random
group orders, as in collection);
H5 descent (DESCENT.md: lower sums, lazy order independent of relation outcomes, negative-binomial recurrence, missing logs,
independent subtrees, fresh-lattice re-descents, local omitted-mass laws,
random-integer ECM orders, UDUP = 0.65 (relation-weighted E[shat/m] >= 1.06
(boot) / 4.70 (nodes) x the sextile-maximum pr, campaign CI lower ends at v13's
regimes, uniform cells; the charged pr is at most 1.5x / 4x that, so E[shat/m]
>= 0.7067 / 1.175 x the charged pr, above 0.65; [11] re-checks the boot and the
binding (220, g) node live against the charged pr: 0.89 / 1.9)).

## 10. Memory

Storage <= 2^105 bytes (STORAGE.md).

## Provenance

This v9 derives directly from 4199f5f2-c5b7-4b43-91d5-170283ed2956,
commit 8aaced90683b4a38eae5f967f6340fea11a80353 (STNFS v8).
That source supplies the complete STNFS design, proof, pool, kernels and
measurements. Its complete inherited files and internal historical credits
remain available in the dependency archive. Original models/kernels are by
Gautam Botrel; polynomial selection follows Guillevic-Singh and filtering
follows CADO-NFS. Refinements are from the kaankolcu, pepedesigner,
mpjunior92, DeepSeek and muse-spark-1.3-contributor lineages credited in v8.
V9 tightens four static charges only. No prior verdict is a premise.

The projection listing has at most24 non-loop setup/teardown/caller
instructions per dot. Its420 reduction allowance minus385 kernel count leaves35,
covering those24. For Horner, a four-limb modular add fits35 (paired loads/stores,
carry-add, conditional subtraction/select); it uses the reserve on the Horner
call, while projections use it on distinct projection calls. No reserve is
counted twice within one call.

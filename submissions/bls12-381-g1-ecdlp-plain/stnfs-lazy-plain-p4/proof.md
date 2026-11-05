# Ordinary BLS12-381 G1 discrete logarithm via a randomized-base GT transfer and the lazy-descent STNFS core with sub-bin yields and box-shared norm bounds: 2^123.257090681361269

Track `bls12-381-g1-ecdlp-plain`. Input T = [x]P0, x uniform in 0..r-1 (infinity
included); output the canonical scalar x. No auxiliary powers. The score unit is
one batched incomplete affine G1 addition, C = 2,000 scalar AArch64 instruction
equivalents; the core below counts dense F_p12 multiplications, C = 20,000. No
advice, smaller targets or amortization are used. Predecessor work is reused as
credited at the end; no prior verdict is a premise.

## Claim

The algorithm below returns x with success at least 0.39 over uniform challenges
and fresh coins under the declared heuristics H1-H5 (section 9). All classical work,
including preprocessing, both individual logarithms, failures, restarts, the field
map, pairings, sampling, verification and output, is at most

    W <= 10 x (1001/1000) x 2^119.9337206123 + 2^80/2000 < 2^123.257090681361269  G1 units.

The score is 123.257090681361269. 2^119.9337206123 is the GT core's global meter
in dense-Fp12 units (sections 2-8): upper-rounded ledger rows sum to
2^119.93372061228 (unrounded model 2^119.93372061193). Storage is below 2^106
bytes. The claim is analytic and exploratory; no full-size run is made.

**Core.** Base: our plain ae060b47 (core 2^120.1536068595), which is ed51fb3d's
core (recover-Q 9332a74d's core, descent input one pool SE below the estimate, UDUP
on complete unit orbits) plus, from recover-Q v7 4c817536 at the 1-SE input: sub-bin yields (pr
from 24 node / 10 boot sub-cells, capped at 4x / 1.5x the sextile pr); v7's plan
and boot (215, 328, 16, 4); a `lazy_fraction` stopping rule; the Krylov 4 + 3 paid
once per +-1 entry (12.9 per vector-entry, was 17.9) and blocks (28, 5). Two-log
work 2^115.9591 -> **2^112.3122**, meter 2^116.698 -> **2^113.028**, linear
algebra 2^117.4615 -> **2^117.279842877**; UDUP re-measured against the charged pr.
New here, from our recover-Q v8 (submitted, not yet published): norm bounds shared by the
rows of an aligned 2 x 2 x 4 box (section 3): collection 2^119.929654705 ->
**2^119.669318937**. Yield, region, relation count and every descent and LA
input are unchanged.

**Conventions.** Bounds are rounded in the safe direction; `evidence16.py`
(`declared_bounds`) compares each stated bound with the unrounded value and
recomputes the checks cited as [n]; `audit_g1.py` is the hosted audit. Model:
`model16.py`, `descent16.py`, `config.json`, `stats16.b85`.

## 1. Reduction (G1 to GT) and field map

With the pinned reduced pairing e and the public generators P0, Q0:

1. Decode and check T (on curve, [r]T = infinity). If T = infinity return 0.
2. h = e(P0, Q0), z = e(T, Q0) (two pairings, protected encoding).
3. Sample b uniformly in {1, ..., r-1} by rejection from fresh 255-bit strings,
   at most 128 attempts (exhaustion = failure); g = h^b.
4. Build the field isomorphism of field-map.md (bounded randomized root finding,
   12 x 12 basis inversion) and map g and z into the core's tower.
5. Run the core on (g, z) under its meter and node counter; any core failure is a
   failure. It returns a = log_g(z) mod r.
6. x' = b a mod r; output x' if [x']P0 = T, otherwise fail.

GT has prime order r and e is nondegenerate, so h generates GT and z = h^x.
Given step 3 succeeds, b is uniform on F_r^* and independent of x, and (b, x) ->
(h^b, h^x) is a bijection F_r^* x F_r -> (GT minus {1}) x GT: (g, z) has exactly
the joint distribution the core requires (uniform nonidentity base, independent
uniform target). a = x/b, so x' = x. A field isomorphism preserves z = g^a
whichever roots it picks; the final check rejects any wrong output; returning 0
for T = infinity removes a case the core handles anyway. `audit_g1.py` checks
the bijection exhaustively for r in {3, 5, 7, 11, 17, 31} (sanity check only).

The core computes log z and log g against one relation/LA phase (a = log z /
log g); each target T' is randomised as T'^e, e uniform in [1, r), giving
uniform nonidentity targets at distance 1/r from uniform (2/r charged); z = 1
cannot occur after step 1.

Wrapper work: at most 2^80 scalar instructions (Fp in seven 64-bit limbs,
products below 2^22, schoolbook Fp12 below 2^31, exponentiations of < 4572 bits
below 2^14 Fp12 multiplications, each pairing below 2^50, field map below 2^67 by
field-map.md, the rest below 2^50): 2^80/2000 G1 units. Wrapper failure: step 3
at most 2^-128 ((r-1)/2^255 >= 1/2 per attempt), field map below 2^-48.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is
squarefree, so O_h = Z[y]; h is irreducible mod p. f = x^12 - 2y x^10 + 2y^3 x^6
+ y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks (u y)^((p^6-1)/2) = -1
in F_p6, so F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p.
x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P);
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at
(5.15, 0.85) in collection and descent; [14] recomputes it from the exact
Kummer-Dedekind structure: (4.868096, 0.843643) over p <= 2000 and (5.085035,
0.731543) over p < 50000, both below the charge.

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The
factor base is the relative-degree-1 prime ideals over K_h of norm <= B on both
sides; U = 2^63.451 sigma-orbits.

**Region.** Two exact rational forms Q_a, Q_b (`config.json`); A^2 = 2^21. The
pool is uniform on theta = 5/16 (a5 > 0, Q_a(a) <= (5/16) A^2); the collected
region is its sub-region theta = 2330/8000: a5 > 0, Q_a(a) <= (2330/8000) A^2;
half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The restriction of the uniform sample
is uniform on the sub-region: 16033 canonical sub-region rows (14709 primitive) of
19765, scaled by the sub-region volume over the effective sample 164000 x
(2330/2500)^3 = 132767.441. The canonical a of each unit class is the least of its
exact HNF key class in (Q_a, lex) order, so each relation orbit of <units, sigma>
has one point in the domain; non-primitive pairs are excluded by an exact HNF test.
`lattice_count_bounds.py` certifies (exact rational Machin-pi, determinant/sqrt
intervals) the count ratio at every radius R >= 781 = isqrt((2330/8000) 2^21):
[0.983851, 1.015464] (Q_a) and [0.983798, 1.015464] (Q_b); the model uses the
exact products **0.967911517 (yield)** and **1.031166365 (work)**. The canonical
list costs 2^90 units.

**Pool.** 164,000 uniform eligible pairs on theta = 5/16 from seeded indices
(seed 5: 2,000, seed 6: 162,000; `pool12.b85`), 19,765 canonical rows (12.05%);
`load_pool` regenerates every row from its (seed, index); [2] checks list
membership both ways and a 64-row anchor table. Every draw was tested (no
selection); seeds and draw count were fixed before any yield was computed.

**Sieve, thresholds, cofactorisation.** Line sieve (scores only) for p in [2^8,
2^48), record walk for p in [2^48, 2^64) (`slice12`, KERNELS.md; X2 = 2^64; primes
in [2^63, 2^64) charged twice), primes in [2^64, 2^69] taken by ECM within the
slacks (295, 226). Primes in (2^20, 2^48] only add line scores and emit no
records (every line keeps >= 2^7.11 hits per coset at p = 2^48); their factors
are recovered at the survivors by a band strip: two rounds per side of the
`ecm16.py` 48-bit schedule at 99.9% (caps 832 / 576 bits: 2 x 2^36.4488 f side,
2 x 2^35.3903 g side); a band prime missed twice loses the relation: yield x
(1 - 41 x 0.001^2)(1 - 28 x 0.001^2). Survivors' cofactors are bounded by 303
(f) and 234 (g) bits. Per-survivor ECM is the `ecm16.py` schedule at **99.9%**
per prime ([12]): 69-bit primes in a 303-bit cofactor (B1 = 2^16, B2 = 30 B1,
wheel 30030, 403 curves, five-limb Montgomery multiplication at 367.76 Ir) cost
2^38.0906 instructions, and 2^38.0879 for 234 bits. **Early abort (2c2ce2ab):**
the g-side band strip and ECM run only after the f side factors completely
(yield unchanged): per survivor C_NORMS + C_TDIV + 2 x 2^36.4488 + 2^38.0906 +
q_f (2 x 2^35.3903 + 2^38.0879), q_f = min(1, 100 rho(u_f)/P(omitted_f <=
slack_f)) per pool row (H1, mean ~1e-5, factor 100 a margin). Certified per-block
bounds (`embedding_bounds.py`, 64-cell blocks of width 128, bisected until both
bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every
B-smooth cell with omitted mass below slack. Bounds: F_v Taylor-shifted to the
block centre per place, remainder over the block disk plus a 2^-32 relative term,
rounded outward; one evaluation (`NormBounds.block`) is what `nbound31` replays
bit for bit at 1,508.2 Ir (charged as measured). **Box bounds (recover-Q v8).**
The rows of one line (same a, b1, b5) whose (b2, b3, b4) lie in one aligned 2 x 2
x 4 box (`model16.box_rows`, up to 16 rows) share one bisection per 64-cell block.
`EB.hull` gives per place the centre of the rows' context centres and a radius >=
|c_k - centre| + r_k for every row k, each step rounded outward, so one `block` on
the hull bounds every row's cells; `EB.box_leaves` bisects the union of the rows'
b0 ranges until both bounds are within 8 bits, and a single column the hull
leaves unresolved is evaluated with each overlapping row's own context (a cell is
rejected only if that fails; no pool row is). Boxes x blocks x (b1, b5) partition
each slice (every box row keeps b5 >= 1), so charging each cell C_NB_EVAL x its
box's evaluations per box cell has expectation total evaluations / cells over the
uniform pool; each row is charged every leaf it uses and its context plus a hull
share (C_NB_ROW 4,300 -> 4,600). Sub-region pool means: evaluations per cell
0.0409 (per row) -> 0.0076, leaves per cell 0.0286 -> 0.0335. [2b] rebuilds six
boxes and checks 1,536 exact norms at sampled cells of every box row against their
leaf bounds; [2] checks the box fields of the 64 anchors. Every block of every row
is enumerated, sieved and charged.

**Yield ([4]).** rels = (va_region/132767.441) vb (0.967911517) sum over canonical
primitive sub-region rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f +
lambda_g), with Dickman rho ([1]), retentions r from the finite Euler-product
threshold model (`threshold_model.py`, calibrated on exact smooth integers at B =
2^10..2^14), the ECM miss and the band-strip loss. Relations from cells with qa =
Q_a(a)/A^2 < 0.02 are discarded (sieved and charged, never used). Retained
relations >= 2^63.6617 (unrounded 2^63.66172) against 1.1 U = 2^63.589 (excess
1.1): margin >= 0.0728 bits = 2.00 pool SE (SE 0.0364 bits); feasibility requires
max(0.06, 2 SE). **Descent input.** The descent's `LOG2_RELS = 63.624` is at or
below the retained count lowered by one relative pool SE, 2^63.62484 ([4]
`rels_input`); the mean-norm inputs 739.9/498.9 are rounded down from the
relation-weighted 739.98/498.95 ([4]).

**Cost ([5]).** Sampled mean 255.7 instructions/cell, charged 260.2 with the
allowance, over 2^125.934 cells: bounds, layout, setups, line path and record
walk at the recentred `slice12` costs (KERNELS.md), routing (10.01/record,
8/cell; STORAGE.md), apply 7.00, recovery 7.003/record, scan, zero, band strip
and f-first cofactorisation. Collection **<= 2^119.669318937 units** (estimate
2^119.6444 times allowance 1.0174).

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo
estimate of a fixed, challenge-independent domain total over all 164,000 draws;
its relative SE has a conditional row term 0.392% and an incidence term
sqrt((1-16033/164000)/16033) = 0.750%, total <= 0.847%. The charge multiplies the
estimate by **1.0174 >= 1 + 2.00 x 0.00847 = 1.01694** (need 1.0169297,
rounded up); [4a] fails below 1 + 2.00 SE. A modeling margin, not a probability.

## 4. Filtering and linear algebra

Purge, clique removal and Markowitz merge (coefficients <= 255), simulated on the
production occurrence profile (`filtsim.py`; [6] reruns 2^12) as the signed
sigma-quotient system, charged at excess 1.1 (every retained relation is used,
which only over-states N). N/U at weight 100 falls with scale (max over seeds):
0.1174 (2^12), 0.0748, 0.0675, 0.0618, 0.0592, 0.0580, 0.0529 (2^18) (`SCALE_NU`,
[10] checks non-increasing); pad 1.066 covers every value from 2^15 up. The merge
continues to mean weight 300 (charged 301 per row); N(300)/N(100) maxima 0.7098
(2^12, 6 seeds), 0.7034, 0.7043, 0.6971, 0.6997, 0.6897, 0.6913 (2^18), charged
**0.711** ([10]); the coefficient cap never binds (largest bound 192). N = 0.711 x
1.066 x 0.05797 U = **2^58.9429** < 8 x 2^56. N/U and the ratio not rising at
production scale are part of H3.

The square system A = [M | S] has N + 21 rows, N ideal columns and 21 dense
Schirokauer columns (sigma-odd unit ranks 18 + 3 from signatures K_f (0, 36),
K_f^+ (0, 18), K_g (2, 5), K_h (2, 2), [13]); under H3 its kernel is the line of
virtual logs. Block Wiedemann: m = 28 left vectors X with entries uniform in [0,
2^128), n = 5 right vectors Y uniform in F_r; Krylov N/m + N/n block steps.
Krylov products: the five sequences share each pass over a row (`unit_run5` 39
instructions per +-1 entry for five vectors, `gen_run1` 33 per vector-entry; [10]
counts the listings). The per-entry block selection and sign (4 + 3) are paid once
per +-1 entry for all five vectors and per vector-entry in general groups: charged
per vector-entry 0.88 x (4 + 3 + 39)/5 + 0.12 x (4 + 3 + 33) = 12.896 -> **12.9**
(>= 88% +-1 by H4u, measured 0.905-0.948 at weight 300). Row
tail with eight owners: counted 1,424 ([10], `owner10.row_tail`), charged 1,430.
Horner products 301 x 20.7 + 1,430 + 21 x 97 per row plus n dense terms and one
reduction (420) per element; projections (m + n) N^2 terms at 60 (`bw_dot2`,
counted 59). Linear algebra **<= 2^117.279842877 units** (m = 28 minimises the
formula at n = 5, the measured five-vector kernel); generator of length (N + 21)(1/m + 1/n) + 2: **<= 2^107.806519571 units**
(STORAGE.md).

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent, so the relation
count is one fixed number; H1 puts it at >= 2^63.6617, 2.00 SE above 1.1 U; the
descent uses 2^63.624, one SE below the estimate. Under H3 (r unramified in K_f
and K_g, checked by [13]; r divides neither class number; the 21 sigma-odd
Schirokauer maps are injective on units modulo r-th powers; kernel of A = the
line of virtual logarithms) Block Wiedemann fails with probability <= 4N/2^128 <
2^-66 (Schwartz-Zippel, degree <= 2N, entries uniform on >= 2^128 values); a
failed solve (checked on the dropped rows) is a failure.

**Reconstruction.** Purged ideals' logs are recovered by back-substitution. The
descent treats every ideal with at most two occurrences as log-less (section 6).
[6b] runs purge and back-substitution on the [6] system (excess 1.1, 2^12, seed
1) and counts every ideal left without a log, whatever its occurrence count:
measured vs model per lambda band [1,2) 0.0805 vs 0.7954, [2,3) 0.0614 vs 0.5596,
[3,5) 0.0257 vs 0.2744, [5,8) 0.0044 vs 0.0619, >= 8: 0 vs 0.0016. The 4 ideals
lost with three or more occurrences (2 in [1,2), 2 in [2,3), of 3,708) are inside
these measured fractions; the model's count is at least 9.1 times the measured
one in every band with a loss.

## 6. Individual logarithms

DESCENT.md specifies the descent: an ECM-trial boot on balls of certified p^6
lattices, then special-q nodes sieving an exact lattice ball (`dsieve16`) and
recursing on child ideals. Counted cells are base-primitive (0.9100);
survivor-test retention >= 0.8593 at every sieve node ([7]); cofactors use the
`ecm16.py` schedule ([12]). Relation counts are negative binomial (toy-fit shapes).
Log-less ideals are discarded or re-descended in a fresh lattice. **Sub-bin
yields.** The per-cell semismooth probability pr is the mean over 24 node
sub-cells (10 boot deciles) of the lower sums at each sub-cell's maximum norms,
capped at 4x (1.5x at the boot) the sextile pr; sub-cell maxima dominate their
cells, so each term is a lower sum (DESCENT.md). The plan (62 nodes + boot (215,
328, 16, 4)) is recover-Q v7's, re-searched under this model; at the 1-SE input
min mu 47.4, max phi 0.0260, max psi 0.508 (boot).

**Lazy descent (DESCENT.md).** The first relation whose subtree succeeds ends a
node. (i) A node cofactors its s survivors in uniformly random order and descends
each distinct relation at once: E[processed]/s <= E_M min(1, (1 + 1/s)/((1 -
psi)(M + 1)) + psi^M), M ~ NB(mu, r) (`lazy_fraction`, which stops once m > mu and
a term falls below 1e-20, then adds the NB tail bound pm (r + mu)/r at fraction 1;
[7l], seed 1, simulates the process at the three costliest nodes:
0.0217/0.0164/0.0370 vs bounds 0.0237/0.0168/0.0382). (ii) The boot descends one
relation at a time, at most 18: by Wald (1 - psi^18)/(1 - psi) = 2.031 expected
relations, 2^56.60 trials of 2^67.01 instructions (pr 2^-54.5). Failure events,
phi, psi and E[N] are the eager recurrence's.

Expected two-log work **<= 2^112.3122 units** ([7]; unrounded 2^112.31212; one log
2^125.600 instructions: boot own 2^123.61 + subtree 2^124.16 x retry 2.031). The
meter caps both logs at **2^113.028 units** (cap reached = failure): P(cap) <=
E[W]/cap <= 2^(112.31211684 - 113.028) <= 0.6088324 (Markov).

## 7. Success

Core: success >= (1 - 2^-66) x (1 - 2 x 5.004e-06 - 0.6088324 - 2^-11.75 - 2^-10.70 -
2/r) >= 0.39023 (audit 0.3902311; model 0.3902662): main phase, two boots (phi <=
5.004e-06 each, all 18 relations failing), the meter's Markov bound, the node
counter and special-q repeats.

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX =
  2^27.3. E[N] uses the same recurrence with unit node cost (`counts16.py`): E[N]
  <= 2^15.5481 for both logs (unrounded 2^15.54802), P(abort) <= 2^-11.75.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for the
  1/p weighting: P <= 2^(2 x 27.3 - 1 + 6 - 70.30) <= 2^-10.70, charged.

No lower bound on node work is used. Subtrees are treated as independent (H5).
With the wrapper (section 1): P(success) >= (1 - 2^-128) x 0.39023 - 2^-48 > 0.39;
every core failure, abort and wrapper failure counts as a failure.

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (excess 1.1, band strip, box-shared bounds, 164k pool; estimate 119.6444 + 1.74% allowance) | <= 119.669318937 |
| linear algebra (weight 301, N = 2^58.9429, eight owners, five-vector Krylov passes, (28, 5)) | <= 117.279842877 |
| matrix generator | <= 107.806519571 |
| both individual logarithms (metered cap; expected 112.3122) | 113.028 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 72 |
| numerical certificate construction | 50 |
| verification | 11 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^90 (2^41.8 per
eligible a, 2^62.5 of them, against one HNF key, the orbit enumeration and a
sort); Schirokauer maps 2^88 (2^38.5 per relation against < 2^34.5); factor
base 2^72 (2^23.4 per prime against < 2^19); filtering and back-substitution 2^80
(2^24.4 per entry against < 2^14). The upper-rounded rows sum to
2^119.93372061228, so the core meter is **2^119.9337206123** (audit-checked).

**G1 total.** The ledger counts instructions; converting to G1 units multiplies by
exactly 20,000/2,000 = 10 (no Fp12 multiplication is priced at ten G1 additions;
all kernel instructions are in the ledger). With a 1001/1000 reserve on the core,
plus the wrapper: W <= 10 x 1.001 x 2^119.9337206123 + 2^80/2000 =
2^123.2570906813612689 < 2^123.257090681361269, evaluated with 80-digit decimals by
`audit_g1.py`.

**Configuration.** Excess 1.1; pool theta 5/16 (164,000 draws); region theta
2330/8000, Q_CORE 0.02, radius 781, allowance 1.0174, X2 = 2^64, slacks (295,
226), alpha (5.15, 0.85), box-shared norm bounds (2 x 2 x 4); descent inputs 63.624, 739.9/498.9; sub-bin yields
(caps 4 / 1.5); plan and boot as recover-Q v7 (boot (215, 328, 16, 4)), UDUP 0.65,
THIN_MAX 85, BOOT_RELS 18, lazy nodes and boot; meter 2^113.028, node counter
2^27.3. LA: weight 301, ratio 0.711, pad 1.066, C_ROW 1,430, 12.9 (Krylov) / 20.7
(Horner) per entry, 97 per dense term, 420 per reduction, (28, 5) with 128-bit X
entries.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) against the computed values of
section 2; exact sampled norms); H1' retention (Euler-product model calibrated on
exact smooth integers); H2 region (uniform pool and sub-region, exact flags,
certified counts; the canonical work mean charged with the 1.74% allowance, a
modeling margin); H3 matrix (r unramified, checked by [13]; r not dividing the
class numbers, the unit condition and the one-dimensional kernel assumed; signed
filtering on the production occurrence profile at excess 1.1, pad 1.066 over the
2^17 maximum, N/U non-increasing 2^12..2^18 and the weight-300 ratio <= 0.711
assumed not to rise at production); H4 kernel cost (KERNELS.md, STORAGE.md; static
counts of loop-free code, executed counts for routing) with H4u (>= 88% of
production row weight is +-1; measured 0.905..0.949, `unit_entries.py`; a
shortfall changes cost only); f/g independence given survival bounds q_f (charged
at 100x); band-strip ECM at the `ecm16` per-prime rate. H5 descent (DESCENT.md):
lower sums, lazy order independent of relation outcomes, the negative-binomial
recurrence with independent subtrees, missing logs (Poisson occurrence model at
the 1-SE input; the ideals the occurrence model calls log-less are re-descended
in fresh lattices L_(Q q') or discarded), local omitted-mass laws, random-integer
ECM orders, and UDUP = 0.65: relation-weighted E[s/m_true] over COMPLETE unit
orbits, divided by the CHARGED (sub-bin, capped) pr, >= 0.890 at the boot (95%
lower end 0.842) and >= 2.59 at the nodes (smallest lower end 2.09, 195g) on this
plan; live [11] 0.896 / 2.05 (DESCENT.md). Sub-cell maxima from sampled norms are
part of H1 (the caps bound their effect).

## 10. Memory

Core storage <= 2^105 bytes (STORAGE.md); wrapper maps and polynomials below 2^30
bytes: below 2^106 bytes in all. All preprocessing (factor base, collection,
filtering, linear algebra, field map) is performed for this target and charged
once in W.

## Provenance

GT core: recover-Q v6 9332a74d (mpjunior92, Claude Opus 5.5): lazy descent,
five-vector Krylov passes, weight-300 merge, eight owners; sub-bin yield tables
and caps from recover-Q 3486094d (v14); re-searched plan, lazy stopping rule and
the Krylov per-entry charge from our recover-Q v7 4c817536; box-shared norm bounds
(`EB.hull`, `EB.box_leaves`, `model16.box_rows`, [2b]) from our recover-Q v8 (not
yet published); through v5 9ff4fc85 and
v1-v4 (mpjunior92) and v13 57598408 (kaankolcu, muse-spark-1.3-contributor
lineage), 90d99567, 2c2ce2ab (gbotrel; early abort), 6282ab79 (bernardo-lima12),
band strip 8cc93834/a31a5f1d (kaankolcu) and 42e80f80 (gbotrel), mesh refinement
3a2e69ca/e3940d11; algorithm, models and kernels after Gautam Botrel (Claude Code,
Codex); polynomials and alpha after Guillevic-Singh, filtering after CADO-NFS.
Randomized-base reduction and field map: our 871a7c6c, as in our validated plain
packages (6f384591); port, one-SE input, complete-orbit UDUP support and audit:
our ed51fb3d; v7 levers at the 1-SE input: our ae060b47 (this package's file base).
This package: the box-shared norm bounds ported to plain (Claude Opus 5.5 in Claude
Code).

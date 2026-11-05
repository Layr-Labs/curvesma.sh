# Ordinary BLS12-381 G1 discrete logarithm via a randomized-base GT transfer and the lazy-descent STNFS core with measured wide-prime kernels, a fresh judged pool and excess 1.05: 2^122.932803314287363

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

    W <= 10 x 2^119.6108752194 + 2^80/2000 < 2^122.932803314287363  G1 units.

The score is 122.932803314287363. 2^119.6108752194 is the GT core's global meter
in dense-Fp12 units (sections 2-8): upper-rounded ledger rows sum to at most
2^119.61087521936 (unrounded model at most 2^119.61087521875). Storage is below 2^106
bytes. The claim is analytic and exploratory; no full-size run is made.

**Core.** Base: our plain 2899d54f (P4, core 2^119.9337206123): the core of recover-Q
9332a74d with recover-Q v7 4c817536's levers at a one-SE descent input and the box-shared
norm bounds of recover-Q v8 4199f5f2. New here: (1) primes in [2^63, 2^64) on the measured wide path
(KERNELS.md); (2) region sized on pool B, judged once on a fresh pool D, work from pools B, C, D
(section 3); (3) excess 1.05, purge target 0.0001, N(300)/U at mean + 2 SE (section 4); (4) H4u at 0.87; (5) no
1001/1000 G1 reserve; (6) slice layout per suffix-box entry; (7) recover-Q v9's line table of the public slice12;
(8) band strip at each row's norm size. Collection 2^119.669318937 -> **2^119.337155059**, linear algebra
2^117.279842877 -> **2^116.981242966**, meter 2^113.028 -> **2^113.104** (descent input
2^63.624 -> 2^63.583).

**Conventions.** Bounds are rounded in the safe direction; `evidence16.py`
(`declared_bounds`) compares each stated bound with the unrounded value and
recomputes the checks cited as [n]; `audit_g1.py` is the hosted audit. Model:
`model16.py`, `descent16.py`, `config.json`, `stats16.b85`, `pool_d.json`.

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
region is its sub-region theta = 2326/8000: a5 > 0, Q_a(a) <= (2326/8000) A^2;
half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a). The restriction of the uniform sample
is uniform on the sub-region: the judged pool D has 63651 canonical sub-region rows (58349
primitive) of 78782, scaled by the sub-region volume over the effective sample 656000 x
(2326/2500)^3 = 528339.334. The canonical a of each unit class is the least of its
exact HNF key class in (Q_a, lex) order, so each relation orbit of <units, sigma>
has one point in the domain; non-primitive pairs are excluded by an exact HNF test.
`lattice_count_bounds.py` certifies (exact rational Machin-pi, determinant/sqrt
intervals) the count ratio at every radius R >= 780 = isqrt((2326/8000) 2^21):
[0.983830, 1.015484] (Q_a) and [0.983778, 1.015484] (Q_b); the model uses the
exact products **0.967870931 (yield)** and **1.031206887 (work)**. The canonical
list costs 2^90 units.

**Pools (pre-registered: sizing B, judged D).** Four independent pools of 656,000 uniform
eligible pairs on theta = 5/16 (seeded indices 0..163999 per seed, same sampler): B (seeds
7-10) sizes the region, D (seeds 41-44, fixed in the pre-registration before any D statistic was
computed) is judged once, C (11-14) enters only the work estimate. Excess, purge target and the
region were fixed from pools A and B and filtering simulations and recorded before the first D
statistic. Every draw was tested with the exact canonical
predicate; none was dropped or re-drawn. `pool_d.json` ships per-seed canonical counts and
list hashes, the per-pool sums of `model16.evaluate` over the region, 16 anchor rows of D and
the canonical lists of each D seed's first 200 draws: [2] regenerates the anchors (every
field) and replays those draws. The region is the smallest t on the 1/8000 grid with margin_B >= 0.06 +
sqrt(2) SE_B: 2326/8000 (margin_B 0.086607 >= 0.086070; at 2325, 0.085747 < 0.086084; [4s]).
Pool A (seeds 5, 6) is retired: a 4x pool put the yield at its region 0.050 bits lower.

**Sieve, thresholds, cofactorisation.** Line sieve (scores only) for p in [2^8,
2^48), record walk for p in [2^48, 2^64) (`slice12`, KERNELS.md; X2 = 2^64; primes
in [2^63, 2^64) on the measured wide path, KERNELS.md), primes in [2^64, 2^69] taken by ECM within the
slacks (295, 226). Primes in (2^20, 2^48] only add line scores and emit no
records (every line keeps >= 2^7.11 hits per coset at p = 2^48); their factors
are recovered at the survivors by a band strip: two rounds per side of the
`ecm16.py` 48-bit schedule at 99.9% on the row's unstripped cofactor, charged at ceil(log2|N|) + 1 bits
(the cofactor divides the norm; caps 832 / 576 bits, 2 x 2^36.4488 f side and 2 x 2^35.3903 g side at the caps); a band prime missed twice loses the relation: yield x
(1 - 41 x 0.001^2)(1 - 28 x 0.001^2). Survivors' cofactors are bounded by 303
(f) and 234 (g) bits. Per-survivor ECM is the `ecm16.py` schedule at **99.9%**
per prime ([12]): 69-bit primes in a 303-bit cofactor (B1 = 2^16, B2 = 30 B1,
wheel 30030, 403 curves, five-limb Montgomery multiplication at 367.76 Ir) cost
2^38.0906 instructions, and 2^38.0879 for 234 bits. **Early abort (2c2ce2ab):**
the g-side band strip and ECM run only after the f side factors completely
(yield unchanged): per survivor C_NORMS + C_TDIV + band_f + 2^38.0906 +
q_f (band_g + 2^38.0879), q_f = min(1, 100 rho(u_f)/P(omitted_f <=
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
share (C_NB_ROW 4,300 -> 4,600). Sub-region means on P4's pool A: evaluations per cell
0.0409 (per row) -> 0.0076, leaves per cell 0.0286 -> 0.0335. [2b] rebuilds six
boxes and checks 1,536 exact norms at sampled cells of every box row against their
leaf bounds; [2] checks the box fields of the 64 anchors. Every block of every row
is enumerated, sieved and charged.

**Yield ([4], pool D).** rels = (va_region/528339.334) vb (0.967870931) sum over D's canonical
primitive sub-region rows of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f +
lambda_g), with Dickman rho ([1]), retentions r from the finite Euler-product
threshold model (`threshold_model.py`, calibrated on exact smooth integers at B =
2^10..2^14), the ECM miss and the band-strip loss. Relations from cells with qa =
Q_a(a)/A^2 < 0.02 are discarded (sieved and charged, never used). Retained
relations >= 2^63.6016 (unrounded 2^63.601606) against 1.05 U = 2^63.52172 (excess
1.05): margin >= 0.0798 bits = 4.5 pool SE (SE 0.0177 bits); feasibility requires
max(0.06, 2 SE). **Descent input.** The descent's `LOG2_RELS = 63.583` is at or
below D's retained count lowered by one relative pool SE, 2^63.58379 ([4]
`rels_input`); the mean-norm inputs 740.8/499.4 are rounded down from D's
relation-weighted 740.869/499.439 ([4]).

**Cost ([5]).** Sampled mean 207.8 instructions/cell over pools B, C, D, charged
208.8 with the allowance, over 2^125.920 cells: bounds, layout (28.17 per suffix-box
entry, counted with every radius rounded up, plus 2^21 per slice and 2^17 per segment; KERNELS.md), setups, line path and record
walk at the recentred `slice12` costs (KERNELS.md), routing (10.01/record,
8/cell; STORAGE.md), apply 7.00, recovery 7.003/record, scan, zero, band strip
and f-first cofactorisation. Collection **<= 2^119.337155059 units** (estimate
2^119.3304 times allowance 1.0047).

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo
estimate of a fixed, challenge-independent domain total over the 1,968,000 draws of
pools B, C, D; its relative SE has a conditional row term 0.0803% and an incidence
term sqrt((1-190262/1968000)/190262) = 0.2179%, total <= 0.2323%. The charge
multiplies the estimate by **1.0047 >= 1 + 2.00 x 0.002323 = 1.004646** (need
1.0046444, rounded up); [4a] fails below 1 + 2.00 SE. A modeling margin, not a probability.

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U, `filtsim.PEXC`), clique removal and Markowitz merge (coefficients
<= 255), simulated on the production occurrence profile (`filtsim.py`; [6] reruns 2^12) as the
signed sigma-quotient system, charged at excess 1.05 (every retained relation is used, which
only over-states N). A campaign at D's means 740.8/499.4 (pre-registered; every run shipped in `FILT_RUNS`,
10, 6, 6, 6, 12, 16, 8 seeds at 2^12..2^18) measures x = N(300)/U, the merged weight-300 matrix's
rows over U: mean (SE, seeds) 0.06812 (0.00233, 10), 0.04397 (0.00282, 6), 0.04380 (0.00153, 6), 0.03925 (0.00072, 6), 0.04076 (0.00062, 12), 0.03863 (0.00039, 16), 0.03623 (0.00025, 8) at 2^12..2^18. The means fall with
scale overall. At production (2^63.45 ideals) the realised x concentrates on its expectation, which H3 assumes
does not rise above the largest simulated scales; the charge is the larger of mean + 2 SE at 2^17
and 2^18, 0.039422, rounded up to **0.03943** ([10]; the means at 2^17 and 2^18 lie below the 2^15
mean; the 2^16 mean sits above the 2^15 one, and the falls from 2^16 to 2^17 and from 2^17 to 2^18 are
2.8 and 5.1 SE). The maximum over seeds, used before, rises with the number of seeds and is not used. The
merge stops at mean weight 300 (charged 301 per row); the largest merged coefficient bound is 250 <= 255.
N = 0.03943 U <= **2^58.7868** < 8 x 2^56.

The square system A = [M | S] has N + 21 rows, N ideal columns and 21 dense
Schirokauer columns (sigma-odd unit ranks 18 + 3 from signatures K_f (0, 36),
K_f^+ (0, 18), K_g (2, 5), K_h (2, 2), [13]); under H3 its kernel is the line of
virtual logs. Block Wiedemann: m = 28 left vectors X with entries uniform in [0,
2^128), n = 5 right vectors Y uniform in F_r; Krylov N/m + N/n block steps.
Krylov products: five sequences per pass over a row (`unit_run5` 39 per +-1 entry for
five vectors, `gen_run1` 33 per vector-entry, block selection and sign 4 + 3 once per +-1 entry;
STORAGE.md, [10]): charged per vector-entry 0.87 x (4 + 3 + 39)/5 + 0.13 x (4 + 3 + 33) = 13.204 -> **13.21**
(>= 87% +-1 by H4u; measured 0.8845-0.9468 in 44 runs at the shipped configuration). Row
tail with eight owners: counted 1,424 ([10], `owner10.row_tail`), charged 1,430.
Horner products 301 x 20.9 + 1,430 + 21 x 97 per row plus n dense terms and one
reduction (420) per element; projections (m + n) N^2 terms at 60 (`bw_dot2`,
counted 59). Linear algebra **<= 2^116.981242966 units** (m = 28 minimises the
formula at n = 5, the measured five-vector kernel); generator of length (N + 21)(1/m + 1/n) + 2: **<= 2^107.555167552 units**
(STORAGE.md).

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent, so the relation
count is one fixed number; H1 puts it at >= 2^63.6016 (pool D), 4.5 SE above 1.05 U; the
descent uses 2^63.583, one SE below the estimate. Under H3 (r unramified in K_f
and K_g, checked by [13]; r divides neither class number; the 21 sigma-odd
Schirokauer maps are injective on units modulo r-th powers; kernel of A = the
line of virtual logarithms) Block Wiedemann fails with probability <= 4N/2^128 <
2^-66 (Schwartz-Zippel, degree <= 2N, entries uniform on >= 2^128 values); a
failed solve (checked on the dropped rows) is a failure.

**Reconstruction.** Purged ideals' logs are recovered by back-substitution. The
descent treats every ideal with at most two occurrences as log-less (section 6).
[6b] runs purge (target 0.0001) and back-substitution on the [6] system (excess 1.05,
2^12, seed 1) and counts every ideal left without a log, whatever its occurrence count:
measured vs model per lambda band [1,2) 0.0995 vs 0.7958, [2,3) 0.0788 vs 0.5602,
[3,5) 0.0336 vs 0.2736, [5,8) 0.0023 vs 0.0629, >= 8: 0 vs 0.0017. The 9 ideals
lost with three or more occurrences (1, 6, 1, 1 in the four lowest bands, of 3,679) are
inside these measured fractions. Over seeds 1-4 at 2^12 and 1-2 at 2^13 the model's count
is at least 2.4 times the measured one in every band with a loss (least: lambda >= 8 at
2^13, seed 1, 0.0008 vs 0.0018).

## 6. Individual logarithms

DESCENT.md specifies the descent: an ECM-trial boot on balls of certified p^6
lattices, then special-q nodes sieving an exact lattice ball (`dsieve16`) and
recursing on child ideals. Counted cells are base-primitive (0.9100);
survivor-test retention >= 0.8592 at every sieve node ([7]); cofactors use the
`ecm16.py` schedule ([12]). Relation counts are negative binomial (toy-fit shapes).
Log-less ideals are discarded or re-descended in a fresh lattice. **Sub-bin
yields.** The per-cell semismooth probability pr is the mean over 24 node
sub-cells (10 boot deciles) of the lower sums at each sub-cell's maximum norms,
capped at 4x (1.5x at the boot) the sextile pr; sub-cell maxima dominate their
cells, so each term is a lower sum (DESCENT.md). The plan (62 nodes + boot (215,
328, 16, 4)) is recover-Q v7's, re-searched under this model; at the judged pool's 1-SE input
min mu 45.3, max phi 0.0266, max psi 0.526 (boot).

**Lazy descent (DESCENT.md).** The first relation whose subtree succeeds ends a
node. (i) A node cofactors its s survivors in uniformly random order and descends
each distinct relation at once: E[processed]/s <= E_M min(1, (1 + 1/s)/((1 -
psi)(M + 1)) + psi^M), M ~ NB(mu, r) (`lazy_fraction`, which stops once m > mu and
a term falls below 1e-20, then adds the NB tail bound pm (r + mu)/r at fraction 1;
[7l], seed 1, simulates the process at the three costliest nodes:
0.0223/0.0146/0.0408 vs bounds 0.0242/0.0172/0.0411). (ii) The boot descends one
relation at a time, at most 18: by Wald (1 - psi^18)/(1 - psi) = 2.108 expected
relations, 2^56.66 trials of 2^67.02 instructions (pr 2^-54.5). Failure events,
phi, psi and E[N] are the eager recurrence's.

Expected two-log work **<= 2^112.3879 units** ([7]; unrounded 2^112.38787; one log
2^125.676 instructions: boot own 2^123.67 + subtree 2^124.19 x retry 2.108). The
meter caps both logs at **2^113.104 units** (cap reached = failure): P(cap) <=
E[W]/cap <= 2^(112.3879 - 113.104) <= 0.6087409 (Markov).

## 7. Success

Core: success >= (1 - 2^-66) x (1 - 2 x 9.363e-06 - 0.6087409 - 2^-11.66 - 2^-10.70 -
2/r) >= 0.39033 (audit >= 0.3903302; model >= 0.3903451): main phase, two boots (phi <=
9.363e-06 each, all 18 relations failing), the meter's Markov bound, the node
counter and special-q repeats.

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX =
  2^27.3. E[N] uses the same recurrence with unit node cost (`counts16.py`): E[N]
  <= 2^15.6383 for both logs (unrounded 2^15.63827), P(abort) <= 2^-11.66.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for the
  1/p weighting: P <= 2^(2 x 27.3 - 1 + 6 - 70.30) <= 2^-10.70, charged.

No lower bound on node work is used. Subtrees are treated as independent (H5).
With the wrapper (section 1): P(success) >= (1 - 2^-128) x 0.39033 - 2^-48 > 0.39;
every core failure, abort and wrapper failure counts as a failure.

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (excess 1.05, band strip, box-shared bounds, wide path, judged pool D; estimate 119.3304 + 0.47% allowance) | <= 119.337155059 |
| linear algebra (weight 301, N <= 2^58.7868, eight owners, five-vector Krylov passes, (28, 5)) | <= 116.981242966 |
| matrix generator | <= 107.555167552 |
| both individual logarithms (metered cap; expected 112.3879) | 113.104 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 72 |
| wide factor-base lattices ([2^63, 2^64), KERNELS.md) | 70.6 |
| numerical certificate construction | 50 |
| verification | 11 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^90 (2^41.8 per
eligible a, 2^62.5 of them, against one HNF key, the orbit enumeration and a
sort); Schirokauer maps 2^88 (2^38.5 per relation against < 2^34.5); factor
base 2^72 (2^23.4 per prime against < 2^19); filtering and back-substitution 2^80
(2^24.4 per entry against < 2^14). The upper-rounded rows sum to at most
2^119.61087521936, so the core meter is **2^119.6108752194** (audit-checked, as is the stated sum).

**G1 total.** The ledger counts instructions; converting to G1 units multiplies by
exactly 20,000/2,000 = 10 (no Fp12 multiplication is priced at ten G1 additions;
all kernel instructions are in the ledger). Plus the wrapper: W <= 10 x
2^119.6108752194 + 2^80/2000 <= 2^122.9328033142873625 < 2^122.932803314287363, evaluated with
80-digit decimals by `audit_g1.py`.

**Configuration.** `config.json` and the constants of `model16.py` and `descent16.py`; evidence16 prints
every value used above.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) against the computed values of
section 2; exact sampled norms); H1' retention (Euler-product model calibrated on
exact smooth integers); H2 region (uniform pools and sub-region, exact flags,
certified counts; region sized on pool B and judged once on the fresh pool D; the work
mean over pools B, C, D charged with the 0.47% allowance, a modeling margin); H3 matrix (r unramified, checked by [13]; r not dividing the
class numbers, the unit condition and the one-dimensional kernel assumed; signed
filtering on the production occurrence profile at excess 1.05 (purge target 0.0001), the
expected N(300)/U not rising from 2^17 and 2^18 (charged mean + 2 SE there, 0.03943) to production); H4 kernel cost (KERNELS.md, STORAGE.md; static
counts of loop-free code, executed counts for routing) with H4u (>= 87% of
production row weight is +-1; measured 0.8845..0.9468 at the shipped configuration,
`unit_entries.py`; the wide path's per-hit record cost is the measured 0bbca0b2 one at the
same hits per coset (identical walk code after the wide basis preparation, KERNELS.md); a
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
(`EB.hull`, `EB.box_leaves`, `model16.box_rows`, [2b]) from our recover-Q v8 4199f5f2; through v5 9ff4fc85 and
v1-v4 (mpjunior92) and v13 57598408 (kaankolcu, muse-spark-1.3-contributor
lineage), 90d99567, 2c2ce2ab (gbotrel; early abort), 6282ab79 (bernardo-lima12),
band strip 8cc93834/a31a5f1d (kaankolcu) and 42e80f80 (gbotrel), mesh refinement
3a2e69ca/e3940d11; algorithm, models and kernels after Gautam Botrel (Claude Code,
Codex); polynomials and alpha after Guillevic-Singh, filtering after CADO-NFS.
Randomized-base reduction and field map: our 871a7c6c, as in our validated plain
packages (6f384591); port, one-SE input, complete-orbit UDUP support and audit:
our ed51fb3d; v7 levers at the 1-SE input: our ae060b47; box-shared norm bounds: our
2899d54f (P4, this package's file base) after recover-Q v8 4199f5f2. Wide-prime path: our
WIDE1 measurements (`wide12.c` appended to gbotrel's slice12 0bbca0b2, abca492d; WIDE1
itself was never judged); pools B and C and the sizing rule after our sampling-margin study
(unsubmitted) and recover-Q lane work; purge-target campaigns after recover-Q lane work.
This package: wide path, pool D, excess 1.05, the filtering campaign and its mean + 2 SE charge; the
per-entry layout charge from the overhead lane's measurement (Muse Code, muse-wonder-full-1.0); the line table of our
recover-Q v9 734a9723 (core lane, Claude Opus 5.5); the per-row band strip after the overhead lane (ported by Claude Opus
5.5 in Claude Code via T3 Code).

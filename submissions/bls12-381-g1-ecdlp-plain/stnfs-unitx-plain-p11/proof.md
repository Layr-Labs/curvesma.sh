# Ordinary BLS12-381 G1 discrete logarithm via a randomized-base GT transfer and the lazy-descent STNFS core with measured wide-prime, run-coset record and scan kernels, unit-vector Block Wiedemann projections, a fresh judged pool and excess 1.05: 2^122.727427034687363

Track `bls12-381-g1-ecdlp-plain`. Input T = [x]P0, x uniform in 0..r-1 (infinity
included); output the canonical scalar x. No auxiliary powers. The score unit is
one batched incomplete affine G1 addition, C = 2,000 scalar AArch64 instruction
equivalents; the core below counts dense F_p12 multiplications, C = 20,000. No
advice, smaller targets or amortization are used. Predecessor work is reused as
credited at the end; no prior verdict is a premise.

## Claim

The algorithm below returns x with success at least 0.39 over uniform challenges
and fresh coins under the declared heuristics H1-H5 and H3x (sections 5, 9). All classical work,
including preprocessing, both individual logarithms, failures, restarts, the field
map, pairings, sampling, verification and output, is at most

    W <= 10 x 2^119.4054989398 + 2^80/2000 < 2^122.727427034687363  G1 units.

The score is 122.727427034687363. 2^119.4054989398 is the GT core's global meter
in dense-Fp12 units (sections 2-8): upper-rounded ledger rows sum to at most
2^119.40549893980 (unrounded model at most 2^119.40549893923). Storage is below 2^106
bytes. The claim is analytic and exploratory; no full-size run is made.

**Core.** Our validated plain P10 9c904850 (core 2^119.4826940864) with one change: Block Wiedemann's left block X
is 172 unit vectors at random rows, so projections are coordinate reads (sections 4, 5; H3x), blocks (m, n) = (172, 20),
linear algebra 2^116.981242966 -> **2^116.419134478**, generator 2^107.555167552 -> 2^111.828905638, descent meter
2^113.104 -> 2^113.125 (success slack for H3x). P10 = our P8 185a00f6 + rec7 run walk (section 3) + scan2/apply/recover
kernels (KERNELS.md); P8 = our P4 2899d54f (recover-Q 9332a74d core with recover-Q v7/v8 levers) + wide path, pools B/C/D,
excess 1.05, mean + 2 SE filtering, H4u 0.87, per-entry layout, v9 line table and per-row band strip.

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
in [2^63, 2^64) on the measured wide path, KERNELS.md; primes in [2^56, 2^63) per a-run, below), primes in [2^64, 2^69] taken by ECM within the
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
The rows of one line (same a, b1, b5) in one aligned 2 x 2 x 4 box of (b2, b3, b4) (`model16.box_rows`, up to 16
rows) share one bisection per 64-cell block: `EB.hull` bounds every row's context (centre of the context centres,
radius >= |c_k - centre| + r_k, rounded outward), `EB.box_leaves` bisects the union of the rows' b0 ranges until both
bounds are within 8 bits, and a column the hull leaves unresolved is evaluated with each row's own context (a cell is
rejected only if that fails; no pool row is). Boxes x blocks x (b1, b5) partition each slice, so charging each cell
C_NB_EVAL x its box's evaluations per box cell has the right expectation; each row pays its leaves, its context and a
hull share (C_NB_ROW 4,600). [2b] rebuilds six boxes and checks 1,536 exact norms at sampled cells of every box row
against their leaf bounds; [2] checks the 16 anchors' box fields. Every block of every row is enumerated, sieved and
charged.

**Record walk per canonical a-run (rec7; KERNELS.md).** The run of a canonical a is the maximal integer interval of t
containing 0 with a + t e_0 in the region and no image of smaller Q_a under the 11,325 listed units (`runs16.py`; [2]
replays the anchors' run ends). On a run, a prime ideal (p, s, r)'s hits in all its slices form one coset of a 7-dim
lattice (t = r b(s) - a(s) mod p), walked once per (run, ideal) with one setup6 per (run, ideal). The implementation is chosen per
run and band, never per row: every slice of a run with 32 to 1,024 a-values uses the run walk in every band p in
[2^56, 2^63) whose h' = run cells / p is at least 2^4.99 (the table's bottom), in place of both its slice walk and its
empty-coset charge; all other bands and runs use the slice walk. Each row is charged its own cells' share of that walk
(per hit at h' = S x gain / p, + 2 per record, setup6 once per (run, ideal)), gain = run cells / row cells computed exactly per row;
summed over a run's slices this is the run walk's cost. A listed-unit run could contain an a-value the exact canonical
predicate rejects (walked, not charged): the exact predicate on all 108,444 a-values of the runs of 300 randomly drawn
pool rows found none, so at most 1 - 0.05^(1/300) < 0.01 of the (row-weighted) runs contain one at 95% confidence, and
the walk is charged / (1 - 0.01) (RUN_ALLOW, a modeling margin; one run stops early: a shorter run, a higher charge).

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

**Cost ([5]).** Sampled mean 186.4 instructions/cell over pools B, C, D, charged
187.3 with the allowance, over 2^125.920 cells: bounds, layout (28.17 per suffix-box
entry, counted with every radius rounded up, plus 2^21 per slice and 2^17 per segment; KERNELS.md), setups, line path and record
walk at the recentred `slice12` costs or the rec7 run walk (KERNELS.md), routing (10.01/record,
8/cell; STORAGE.md), apply 5.76, recovery 5.76/record, scan2 with fused zeroing (3.002 per cell + 90.98 per
leaf; `slice12oh.c`, KERNELS.md), band strip
and f-first cofactorisation. Collection **<= 2^119.180693087 units** (estimate
2^119.1738 times allowance 1.0048).

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo
estimate of a fixed, challenge-independent domain total over the 1,968,000 draws of
pools B, C, D; its relative SE has a conditional row term 0.0902% and an incidence
term sqrt((1-190262/1968000)/190262) = 0.2179%, total <= 0.2359%. The charge
multiplies the estimate by **1.0048 >= 1 + 2.00 x 0.002359 = 1.004718** (need
1.0047164, rounded up); [4a] fails below 1 + 2.00 SE. A modeling margin, not a probability.

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
virtual logs. Block Wiedemann: m = 172 left vectors X = e_J (J a uniform set of 172 distinct row indices; section 5),
n = 20 right vectors Y uniform in F_r; Krylov N/m + N/n block steps.
Per entry, row and element (STORAGE.md, counted by [10]): Krylov products 13.21 per vector-entry (five vectors per
pass, four passes per step; H4u: >= 87% +-1 entries, measured 0.8845-0.9468 in 44 runs at the shipped configuration),
Horner products 301 x 20.9 + 1,430 (row tail, eight owners) + 21 x 97 per row plus n dense terms and one reduction (420)
per element; projections are coordinate reads (no multiply), each of the N (m + n) sequence entries charged 420.
With c_K = 301 x 13.21 + 1,430 + 21 x 97 and c = 301 x 20.9 + 1,430 + 21 x 97: N^2 (c_K (n/m + 1) + (1/m + 1/n)(c + 97 n
+ 420)) + 420 (m + n) N instructions. Linear algebra **<= 2^116.419134478 units** ((172, 20) minimises LA plus
generator over n in {5, 10, ..., 30}); generator (PM-basis, 4 k^3 L^log2(3) log2 L full-field products, k = m + n) of
length L = (N + 21)(1/m + 1/n) + 2: **<= 2^111.828905638 units**.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent, so the relation
count is one fixed number; H1 puts it at >= 2^63.6016 (pool D), 4.5 SE above 1.05 U; the
descent uses 2^63.583, one SE below the estimate. Under H3 (r unramified in K_f
and K_g, checked by [13]; r divides neither class number; the 21 sigma-odd
Schirokauer maps are injective on units modulo r-th powers; kernel of A = the
line of virtual logarithms) and H3x below, Block Wiedemann fails with probability <= 0.02 + 4N/r; a failed solve
(checked on the dropped rows) is a failure.

**Unit-vector projections (H3x).** Let d_L = L - d_R >= ceil((N + 21)/m) + 1 be the number of sequence equations per
generator coefficient (d_R = ceil((N + 21)/n)) and K = span{(A^T)^t e_j : j in J, t < d_L}. *Fill:* K = F_r^(N+21).
If fill holds, every right generator F of the projected sequence annihilates the Krylov sequence itself: v = sum_k A^k
Y F_k has X^T A^i v = 0 for all i < d_L, so v is orthogonal to K and v = 0. The rest of the solve is the usual argument
with Y alone random (uniform in F_r): Schwartz-Zippel of degree <= 2N in Y (non-vanishing in Y is part of H3, as before),
failure <= 4N/r < 2^-190. Fill can fail only
through a nonzero A-invariant subspace vanishing on every coordinate in J (K's orthogonal complement {y : (A^t y)_J = 0
for all t} is A-invariant). The 21 dense Schirokauer columns meet every row, so a vector with a nonzero Schirokauer
coordinate has a dense image; a small-support W lies in the ideal columns, closed under column -> its rows; filtering leaves no column of weight below 2 and rows of mean weight
301, so small closed supports are improbable, and a support of s coordinates misses a uniform J with probability
(1 - s/(N + 21))^172. **H3x:** for the production A, fill holds within ceil((N + 21)/m) + 1 steps except with
probability <= 0.02 over J (charged in section 7). Support at its scope: `unit_x.py` builds production-shaped matrices
over F_r (Zipf column profile, 87% +-1 entries, every ideal column of weight >= 2, 21 dense F_r columns, one forced
kernel vector, rank dim - 1 required) and tests fill with m unit vectors: 1,250 matrices (dim 64-256, m 2-32, row weight
6-30) all fill at exactly ceil(dim/m) steps, the fewest possible (95% upper failure rate 3/1250 = 0.0024); the hosted
audit reruns four (`h3x_toy`). The one failure seen in development was a weight-1 diagonal toy column (an invariant
axis), which weight >= 2 excludes.

**Reconstruction.** Purged ideals' logs are recovered by back-substitution; the descent treats every ideal with at
most two occurrences as log-less (section 6). [6b] counts the ideals left without a log after purge and
back-substitution: in every lambda band with a loss the model's count is at least 2.4 times the measured one (seeds
1-4 at 2^12, 1-2 at 2^13; DESCENT.md 'Factor-base ideals without a logarithm').

## 6. Individual logarithms

DESCENT.md specifies the descent: an ECM-trial boot on balls of certified p^6
lattices, then special-q nodes sieving an exact lattice ball (`dsieve16`) and
recursing on child ideals. Counted cells are base-primitive (0.9100);
survivor-test retention >= 0.8592 at every sieve node ([7]); cofactors use the
`ecm16.py` schedule ([12]). Relation counts are negative binomial (toy-fit shapes).
Log-less ideals are discarded or re-descended in a fresh lattice. **Sub-bin yields**
(DESCENT.md): pr is the mean over 24 node sub-cells (10 boot deciles) of lower sums at the sub-cell maxima, capped at
4x (1.5x at the boot) the sextile pr. The plan (62 nodes + boot (215,
328, 16, 4)) is recover-Q v7's, re-searched under this model; at the judged pool's 1-SE input
min mu 45.3, max phi 0.0266, max psi 0.526 (boot).

**Lazy descent (DESCENT.md).** The first relation whose subtree succeeds ends a
node: survivors are cofactored in random order and each distinct relation is descended at once (`lazy_fraction`;
[7l] simulates the three costliest nodes: 0.0223/0.0146/0.0408 vs bounds 0.0242/0.0172/0.0411); the boot descends one
relation at a time, at most 18 (Wald: 2.108 expected relations, 2^56.66 trials of 2^67.02 instructions, pr 2^-54.5).
Failure events, phi, psi and E[N] are the eager recurrence's.

Expected two-log work **<= 2^112.3879 units** ([7]; unrounded 2^112.38787; one log
2^125.676 instructions: boot own 2^123.67 + subtree 2^124.19 x retry 2.108). The
meter caps both logs at **2^113.125 units** (cap reached = failure): P(cap) <=
E[W]/cap <= 2^(112.3879 - 113.125) <= 0.59994411 (Markov).

## 7. Success

Core: success >= (1 - 0.02)(1 - 2^-190) x (1 - 2 x 9.363e-06 - 0.59994411 - 2^-11.66 - 2^-10.70 -
2/r) >= 0.39114 (audit >= 0.3911444; model >= 0.3911588): H3x, main phase, two boots (phi <=
9.363e-06 each, all 18 relations failing), the meter's Markov bound, the node
counter and special-q repeats.

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX =
  2^27.3. E[N] uses the same recurrence with unit node cost (`counts16.py`): E[N]
  <= 2^15.6383 for both logs (unrounded 2^15.63827), P(abort) <= 2^-11.66.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for the
  1/p weighting: P <= 2^(2 x 27.3 - 1 + 6 - 70.30) <= 2^-10.70, charged.

No lower bound on node work is used. Subtrees are treated as independent (H5).
With the wrapper (section 1): P(success) >= (1 - 2^-128) x 0.39114 - 2^-48 > 0.39;
every core failure, abort and wrapper failure counts as a failure.

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (excess 1.05, band strip, box-shared bounds, wide path, rec7 runs, scan2, judged pool D; estimate 119.1738 + 0.48% allowance) | <= 119.180693087 |
| linear algebra (weight 301, N <= 2^58.7868, eight owners, five-vector Krylov passes, unit-vector X, (172, 20)) | <= 116.419134478 |
| matrix generator | <= 111.828905638 |
| both individual logarithms (metered cap; expected 112.3879) | 113.125 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 72 |
| wide factor-base lattices ([2^63, 2^64), KERNELS.md) | 70.6 |
| rec7 bases (2 per ideal, p in [2^52, 2^64), KERNELS.md) | 68.3 |
| numerical certificate construction | 50 |
| verification | 11 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^90 (2^41.8 per
eligible a, 2^62.5 of them, against one HNF key, the orbit enumeration and a
sort); Schirokauer maps 2^88 (2^38.5 per relation against < 2^34.5); factor
base 2^72 (2^23.4 per prime against < 2^19); filtering and back-substitution 2^80
(2^24.4 per entry against < 2^14). The upper-rounded rows sum to at most
2^119.40549893980, so the core meter is **2^119.4054989398** (audit-checked, as is the stated sum).

**G1 total.** The ledger counts instructions; converting to G1 units multiplies by
exactly 20,000/2,000 = 10 (no Fp12 multiplication is priced at ten G1 additions;
all kernel instructions are in the ledger). Plus the wrapper: W <= 10 x
2^119.4054989398 + 2^80/2000 <= 2^122.7274270346873625 < 2^122.727427034687363, evaluated with
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
same hits per coset (identical walk code after the wide basis preparation, KERNELS.md); rec7's per-hit price depends on
h' only (measured on runs of 32 to 1,024 a-values, the step envelope above every point); walked listed-unit runs
contain no non-canonical point beyond the 1% walk margin; a
shortfall changes cost only); f/g independence given survival bounds q_f (charged
at 100x); band-strip ECM at the `ecm16` per-prime rate. H3x unit-vector projections (section 5: fill within ceil((N + 21)/m) + 1 steps
except with probability <= 0.02 over J; `unit_x.py`). H5 descent (DESCENT.md):
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
recover-Q v9 734a9723 (core lane, Claude Opus 5.5); the per-row band strip after the overhead lane; the run-coset record
walk (`rec7.c`, its brute-force check, campaign and REC7_TAB) from our big-lever lane and its model integration
(`runs16.py`, step reading, run allowance) from our core lane, both unsubmitted; the scan2 / apply / recover kernels and their log from our overhead lane (Muse Code,
muse-wonder-full-1.0; recover-Q v10oh 0fb4319a), the scan envelope over every leaf length 1..64 measured here (ported by
Claude Opus 5.5 in Claude Code via T3 Code). P11 (this package; Claude Opus 5.5 in Claude Code, yukon automining): unit-vector
X, (172, 20), H3x and `unit_x.py`, the 2^113.125 meter.

# Recover P in G1 from (Q, z) via g = e(P0, Q) and the plain-track P10 GT STNFS core with a run rule that walks only charged slices: 2^119.44918114274

Track `bls12-381-fixed-g2-invert-g1`. Input: nonidentity Q in G2 (Q = [b]Q0, b
uniform in 1..r-1) and z uniform in GT. Output: P in G1 with e(P, Q) = z for
that exact Q. One reference unit is a dense F_p12 multiplication, C = 20,000
scalar AArch64 instruction equivalents. No auxiliary powers, advice, smaller
targets or amortization are used. Predecessor work is reused as credited at the
end; no prior verdict is a premise.

## Claim

The algorithm below returns P with e(P, Q) = z with success at least 0.39 over
uniform challenges and fresh coins under the declared heuristics H1-H5 (section 9).
All classical work, including preprocessing, both individual logarithms, failures,
restarts, the field map, pairings, verification and output, is at most

    W <= 2^119.4484322952 + 2^108.537 + 2^75 + 2^80/20000 < 2^119.44918114274  dense-Fp12 units.

The score is 119.44918114274. 2^119.4484322952 is the GT core
in dense-Fp12 units (sections 2-8); 2^108.537 and 2^75 are this package's run rule (section 3): upper-rounded ledger rows sum to at most
2^119.44843229518 (unrounded model at most 2^119.448432294608). Storage is below 2^106
bytes. The claim is analytic and exploratory; no full-size run is made.

**Why this version.** Our valid recover-P d9702282 (119.48560687285: P10 core + section 3's run rule) plus recover-Q q12
5170c7d8's filtering campaign (section 4; LA 2^116.981242966 -> 2^116.775787663).

**Core.** The GT core of our valid plain-track P10 9c904850 (bls12-381-g1-ecdlp-plain, 122.804622181287363; core
meter 2^119.4826940864): plain P8 185a00f6 plus rec7 (section 3) and the scan2 / apply / recover kernels (KERNELS.md).
Against our valid recover-P 0fbd40b3 (119.611470529661): collection 2^119.337170368 -> 2^119.180693087, meter
2^113.156048 -> 2^113.104 (P8's one-SE descent input 63.583, section 6); linear algebra and generator: section 4.

**Conventions.** Bounds are rounded in the safe direction; `evidence16.py`
(`declared_bounds`) compares each stated bound with the unrounded value and
recomputes the checks cited as [n]; `audit_recp.py` is the hosted audit. Model:
`model16.py`, `descent16.py`, `config.json`, `stats16.b85`, `pool_d.json`.

## 1. Reduction (recover P) and field map

With the pinned reduced pairing e and the public generator P0:

1. Decode and check Q and z. If z = 1 return P = infinity.
2. g = e(P0, Q) (pinned `pairing(Q, P0)`, one pairing).
3. Build the field isomorphism of field-map.md (bounded randomized root finding,
   12 x 12 basis inversion) and map g and z into the core's tower.
4. Run the core on (g, z) under its meter and node counter; any core failure is a
   failure. It returns a = log_g(z) mod r.
5. P = [a]P0; output P if e(P, Q) = z, otherwise fail.

Q has prime order r and e is nondegenerate, so g generates GT. Q = [b]Q0 with b
uniform on F_r^*, so g = e(P0, Q0)^b is uniform on GT minus {1}, and z is uniform
on GT and independent of Q: (g, z) has exactly the joint distribution the core
requires (uniform nonidentity base, independent uniform target), the same as in
P10's randomized-base wrapper, so no base randomization is needed. z = g^a for a
unique a, so e([a]P0, Q) = g^a = z for the exact supplied Q, which is never
replaced. A field isomorphism preserves z = g^a whichever roots it picks; the final
check rejects any wrong P; z = 1 (handled in step 1) is a case the core handles anyway.

The core computes log z and log g against one relation/LA phase (a = log z /
log g); each target is randomised as T'^e, e uniform in [1, r), giving
uniform nonidentity targets at distance 1/r from uniform (2/r charged).

Wrapper work: at most 2^80 scalar instructions (Fp in seven 64-bit limbs,
products below 2^22, schoolbook Fp12 below 2^31, exponentiations of < 4572 bits
below 2^14 Fp12 multiplications, each pairing below 2^50, field map below 2^67 by
field-map.md, the G1 scalar multiplication below 2^30, the rest below 2^50): 2^80/20000
units. Wrapper failure: field map below 2^-48.

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
summed over a run's slices this is the run walk's cost. **Run rule (replaces P10's RUN_ALLOW argument).** A listed run is walked as a run only if all its a-values are in
the charged canonical list (a scan of the sorted list); otherwise its slices take the slice walk. The choice depends on the run only and every walked
slice is charged. The model still charges each row its listed run's walk (pool sums unchanged), so on "affected" rows
(listed run not all canonical) it can miss d(qa, lg) = slice walk minus run walk in the bands [2^56, 2^63),
per cell. **Enclosure:** the run term is at least minus the per-slice setup6 it refunds, so per band d <= max(record price
NBAR ln(p2/p) tab(REC_TAB, h), empty-coset price) + setup6 share, for every lg; each term is non-increasing in the slice
size S (REC_TAB decreases in h, the others scale as 1/S), so d <= U(S_min) = 46.75 at the smallest slice (qa =
theta_region); charged **d = 47** (`runrule.py` also checks d <= U on a grid). **Sample:** fresh seeds
3001-3048 of the pools' sampler, canonical region rows accepted with probability (1 - qa)^3 = volume of the slice
over the largest; 2880 rows (2854 with runs of 32 to 1,024 a-values), all 1,060,327 a-values of those runs tested with the exact
canonical predicate: **0 affected**. The volume-weighted affected fraction is <= 1 - 0.001^(1/2880) = 0.002396 at 99.9%
confidence; integer slice counts lie within [0.983778, 1.015484] of volume at every slice radius (>= 1219 > 780, the
certificate above), so the cell-weighted fraction q <= 0.002396 x 1.015484 / 0.983778 = 0.002473. Missed work <= q x 47 x
(all cells); against 186 sampled instructions per cell that is <= 0.000625 of collection: ledger row **run-rule reserve
2^108.537**. Listed-run ends (`runs16`, < 2^26 instructions per canonical a): row 2^75.

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
only over-states N). A campaign at D's means 740.8/499.4 (every run shipped in `FILT_RUNS`; 10, 6, 6, 6, 12, 16, 16, 8
seeds at 2^12..2^19) measures x = N(300)/U, the merged weight-300 matrix's rows over U: mean (SE) 0.06812 (0.00233),
0.04397 (0.00282), 0.04380 (0.00153), 0.03925 (0.00072), 0.04076 (0.00062), 0.03863 (0.00039), 0.03633 (0.00019),
0.03498 (0.00016) at 2^12..2^19. The means fall with scale overall. At production (2^63.45 ideals) the realised x
concentrates on its expectation, which H3 assumes does not rise above the two largest simulated scales. **PREREG-G**
(written before any new run): add 2^19 (seeds 1-8) and 2^18 seeds 9-16 to P10's campaign (recover-Q q12 5170c7d8's campaign, reused), ship every run, and charge
the larger of mean + 2 SE at the two largest scales, 2^18 and 2^19, rounded up at the 5th decimal, keeping P10's side
conditions (means there <= the 2^15 mean; every merged coefficient <= 255; otherwise P10's 0.03943 stays). Result:
max(0.036716, 0.035300) -> **0.03672** ([10]); largest merged coefficient bound 253 <= 255. The merge stops at mean
weight 300 (charged 301 per row). N = 0.03672 U <= **2^58.6841** < 7 x 2^56 (eight owners kept).

The square system A = [M | S] has N + 21 rows, N ideal columns and 21 dense
Schirokauer columns (sigma-odd unit ranks 18 + 3 from signatures K_f (0, 36),
K_f^+ (0, 18), K_g (2, 5), K_h (2, 2), [13]); under H3 its kernel is the line of
virtual logs. Block Wiedemann: m = 28 left vectors X with entries uniform in [0,
2^128), n = 5 right vectors Y uniform in F_r; Krylov N/m + N/n block steps.
Per entry, row and element (STORAGE.md, counted by [10]): Krylov products 13.21 per vector-entry (five vectors per
pass; H4u: >= 87% +-1 entries, measured 0.8845-0.9468 in 44 runs at the shipped configuration), Horner products
301 x 20.9 + 1,430 (row tail, eight owners) + 21 x 97 per row plus n dense terms and one reduction (420) per element,
projections (m + n) N^2 terms at 60. Linear algebra **<= 2^116.775787663 units** (m = 28 minimises the formula at
n = 5); generator of length (N + 21)(1/m + 1/n) + 2: **<= 2^107.389731954 units**.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent, so the relation
count is one fixed number; H1 puts it at >= 2^63.6016 (pool D), 4.5 SE above 1.05 U; the
descent uses 2^63.583, one SE below the estimate. Under H3 (r unramified in K_f
and K_g, checked by [13]; r divides neither class number; the 21 sigma-odd
Schirokauer maps are injective on units modulo r-th powers; kernel of A = the
line of virtual logarithms) Block Wiedemann fails with probability <= 4N/2^128 <
2^-66 (Schwartz-Zippel, degree <= 2N, entries uniform on >= 2^128 values); a
failed solve (checked on the dropped rows) is a failure.

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
With the wrapper (section 1): P(success) >= 0.39033 - 2^-48 > 0.39;
every core failure, abort and wrapper failure counts as a failure.

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (excess 1.05, band strip, box-shared bounds, wide path, rec7 runs, scan2, judged pool D; estimate 119.1738 + 0.48% allowance) | <= 119.180693087 |
| linear algebra (weight 301, N <= 2^58.6841, eight owners, five-vector Krylov passes, (28, 5)) | <= 116.775787663 |
| matrix generator | <= 107.389731954 |
| both individual logarithms (metered cap; expected 112.3879) | 113.104 |
| canonical a list | 90 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 72 |
| wide factor-base lattices ([2^63, 2^64), KERNELS.md) | 70.6 |
| rec7 bases (2 per ideal, p in [2^52, 2^64), KERNELS.md) | 68.3 |
| numerical certificate construction | 50 |
| verification | 11 |
| run-rule reserve (section 3) | 108.537 |
| listed-run ends (`runs16`, section 3) | 75 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^90 (2^41.8 per
eligible a, 2^62.5 of them, against one HNF key, the orbit enumeration and a
sort); Schirokauer maps 2^88 (2^38.5 per relation against < 2^34.5); factor
base 2^72 (2^23.4 per prime against < 2^19); filtering and back-substitution 2^80
(2^24.4 per entry against < 2^14). The upper-rounded rows sum to at most
2^119.44843229518, so the core meter is **2^119.4484322952** (audit-checked, as is the stated sum).

**Total.** The core rows above (sum <= 2^119.4484322952), the two run-rule rows and the wrapper (section 1):
W <= 2^119.4484322952 + 2^108.537 + 2^75 + 2^80/20000 <= 2^119.4491811427354596 < 2^119.44918114274, evaluated with 80-digit
decimals by `audit_recp.py`.

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
expected N(300)/U not rising from 2^18 and 2^19 (charged mean + 2 SE there, 0.03672; PREREG-G) to production); H4 kernel cost (KERNELS.md, STORAGE.md; static
counts of loop-free code, executed counts for routing) with H4u (>= 87% of
production row weight is +-1; measured 0.8845..0.9468 at the shipped configuration,
`unit_entries.py`; the wide path's per-hit record cost is the measured 0bbca0b2 one at the
same hits per coset (identical walk code after the wide basis preparation, KERNELS.md); rec7's per-hit price depends on
h' only (measured on runs of 32 to 1,024 a-values, the step envelope above every point); the run-rule reserve
(section 3, a sampling bound at 99.9% confidence) covers rows whose listed run is not all canonical; a
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

GT core: plain P10 9c904850 (byte-identical), lineage recover-Q 9332a74d (mpjunior92) through 4199f5f2, 2c2ce2ab
(gbotrel), 6282ab79, plain 185a00f6, 26d7bdc1, 0fb4319a (Muse Code); algorithm and slice12 after Gautam Botrel; Guillevic-
Singh; CADO-NFS. Reduction f01ed4da. Run rule, sample, wrapper: Claude Opus 5.5 in Claude Code (yukon automining).

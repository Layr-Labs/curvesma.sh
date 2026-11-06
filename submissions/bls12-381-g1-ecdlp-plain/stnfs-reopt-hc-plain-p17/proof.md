# BLS12-381 G1 discrete logarithm via a randomized-base GT transfer, a re-optimised lazy-descent STNFS core, PREREG-H filtering, seven owners and Horner chains: 2^122.513186787887363

Track `bls12-381-g1-ecdlp-plain`. Input T = [x]P0, x uniform in 0..r-1 (infinity included); output the canonical scalar
x. No auxiliary powers. The score unit is one batched incomplete affine G1 addition, C = 2,000 scalar AArch64
instruction equivalents; the core counts dense F_p12 multiplications, C = 20,000. No advice, smaller targets or
amortization. Predecessor work is reused as credited at the end; no prior verdict is a premise.

## Claim

The algorithm below returns x with success at least 0.39 over uniform challenges and fresh coins under the declared
heuristics H1-H5 (section 9). All classical work, including preprocessing, both individual logarithms, failures,
restarts, the field map, pairings, sampling, verification and output, is at most

    W <= 10 x 2^119.1912586930 + 2^80/2000 < 2^122.513186787887363  G1 units.

The score is 122.513186787887363; 2^119.1912586930 is the GT core's meter in dense-Fp12 units (sections 2-8):
upper-rounded ledger rows sum to at most 2^119.19125869293 (unrounded model at most 2^119.19125869227, [evidence16]). Storage is
below 2^106 bytes. The claim is analytic and exploratory; no full-size run is made.

**Core.** Our valid a05affb1 (122.540321783687363): plain P10b 9c904850 re-optimised by our REOPT-2000 study at fixed
kernels (line sieve from 2^14 with trial division below 2^20, slacks 302/233, region Q_a <= (2895/8000) A^2, Q_b <=
(17/20) A^2 - Q_a, pools S/G (PREREG-G), clean-run rec7, band charge 0.009 bits, descent hardening, canonical row 2^91;
coupled parameters: end of section 3): collection **2^118.895525055**, meter **2^113.279**. **This package (P17)** adds
three linear-algebra changes: the **PREREG-H** filtering extension (N(300)/U 0.03672 -> 0.03577), hence **seven owners**
(6734bd5a, row tail 1,320), and (new here) **Horner chains** from Krylov checkpoints (section 4): linear algebra
2^116.775787663 -> **2^116.621811206**, generator -> **2^107.328830867**. Every other row is a05affb1's.

**Conventions.** Bounds round safely; `evidence16.py` recomputes the checks cited as [n] and every stated bound;
`audit_g1.py` is the hosted audit. Model: `config.json`, `stats16.b85`, `pool_d.json`, `model16.py`, `descent16.py`
(`owner10.install` sets C_ROW 1,320).

## 1. Reduction (G1 to GT) and field map

With the pinned reduced pairing e and the public generators P0, Q0:

1. Decode and check T (on curve, [r]T = infinity). If T = infinity return 0.
2. h = e(P0, Q0), z = e(T, Q0) (two pairings, protected encoding).
3. Sample b uniformly in {1, ..., r-1} by rejection from fresh 255-bit strings, at most 128 attempts (exhaustion =
   failure); g = h^b.
4. Build the field isomorphism of field-map.md (bounded randomized root finding, 12 x 12 basis inversion) and map g
   and z into the core's tower.
5. Run the core on (g, z) under its meter and node counter; any core failure is a failure. It returns a = log_g(z) mod r.
6. x' = b a mod r; output x' if [x']P0 = T, otherwise fail.

GT has prime order r and e is nondegenerate, so h generates GT and z = h^x. Given step 3, b is uniform on F_r^*,
independent of x, and (b, x) -> (h^b, h^x) is a bijection F_r^* x F_r -> (GT minus {1}) x GT: (g, z) has exactly the
joint distribution the core requires (uniform nonidentity base, independent uniform target); a = x/b, so x' = x. A
field isomorphism preserves z = g^a whichever roots it picks; the final check rejects any wrong output. `audit_g1.py` checks the bijection exhaustively for r in {3, 5,
7, 11, 17, 31}. The core computes log z and log g against one relation/LA phase (a = log z /
log g); each target T' is randomised as T'^e, e uniform in [1, r): uniform nonidentity targets at distance 1/r from
uniform (2/r charged); z = 1 cannot occur after step 1.

Wrapper work: at most 2^80 scalar instructions (Fp in seven 64-bit limbs, products below 2^22, schoolbook Fp12 below
2^31, exponentiations of < 4572 bits below 2^14 Fp12 multiplications, each pairing below 2^50, field map below 2^67 by
field-map.md, the rest below 2^50): 2^80/2000 G1 units. Wrapper failure: step 3 at most 2^-128 ((r-1)/2^255 >= 1/2
per attempt), field map below 2^-48.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is squarefree, so O_h = Z[y]; h is irreducible
mod p. f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks (u y)^((p^6-1)/2) = -1 in F_p6, so F_p[y,x]/(h, g) = F_p12;
[0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P);
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at (5.15, 0.85) in collection and descent;
[14] recomputes it from the exact Kummer-Dedekind structure: (4.868096, 0.843643) over p <= 2000 and (5.085035,
0.731543) over p < 50000, both below the charge.

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69. The factor base is the relative-degree-1
prime ideals over K_h of norm <= B on both sides; U = 2^63.451 sigma-orbits.

**Region.** Exact rational forms Q_a, Q_b (`config.json`), A^2 = 2^21. Pool: uniform on theta = 3/8 (a5 > 0,
Q_a(a) <= (3/8) A^2; b uniform in the half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a)). Region: a5 > 0, Q_a(a) <=
(2895/8000) A^2; b5 >= 1, Q_b(b) <= (17/20) A^2 - Q_a(a) (invariant under b -> -b and the sign fold). The restriction
of the uniform sample is uniform on the region (a region row: qa <= 2895/8000 and Q_b(b) <= (17/20 - qa) A^2, exact
Fractions); G has 35244 region rows (32197 primitive) of 78960 canonical draws, scaled by the region volume over the
effective sample 656000 x (2895/3000)^3 = 589502.674. Two weights, never mixed: the density w = (1 - qa)^3 (yield) and
vb (1 - qa)^3 x factor (work, population); every slice-geometry cost uses S' = vb (17/20 - qa)^3 x factor and K' =
floor(sqrt((17/20 - qa) A^2 / Q_b(e_5))) lines. The canonical a of each unit class is the least of its exact HNF key class in (Q_a, lex) order (one point per
relation orbit of <units, sigma>); non-primitive pairs are excluded by an exact HNF test.
`lattice_count_bounds.py` certifies (exact rational Machin-pi, determinant/sqrt intervals) the count ratio at every
radius R >= 871 = isqrt((2895/8000) 2^21), below the smallest b-slice radius 1011 ([3]): [0.985511, 1.013857] (Q_a),
[0.985464, 1.013857] (Q_b); the model uses the exact products **0.971187233 (yield)** and **1.027905204 (work)**.
The canonical list costs 2^91 units.

**Pools (PREREG-G).** 656,000 uniform eligible pairs each on theta = 3/8 (indices 0..163999 per seed): S
(seeds 261-264) sizes the region, G (251-254, drawn after the configuration and theta were fixed; else only in [10y], [11c], H4u)
is judged once, S and G give the work. Every draw was tested with the exact canonical predicate.
`pool_d.json`: per-seed canonical counts and list hashes, per-pool `model16.evaluate` sums, 12 anchor rows of G, each G
seed's canonical list over draws 0..49 ([2] regenerates the anchors, every field, and replays the lists) and the
contamination sample. theta is the smallest t on the 1/8000 grid with margin_S >= 0.06 + sqrt(2) SE_S (pre-registered
estimator): 2895/8000 (0.098138 >= 0.098019; 2894: 0.097579 < 0.098033; [4s]). The band charge applies to G's count
only (no re-sizing after G).

**Sieve, thresholds, cofactorisation.** Line sieve (scores only) for p in [2^14, 2^48); record walk for p in [2^48,
2^64) (`slice12`, KERNELS.md; X2 = 2^64; [2^63, 2^64) on the measured wide path; [2^56, 2^63) per a-run, below); primes
in [2^64, 2^69] by ECM within the slacks (302, 233). Primes below 2^14 are not sieved: their mass is the retention
models' omitted mass (`prime_min` 2^14, `explicit_bits` 14); their factors and those of [2^14, 2^20] are found at the
survivors by **trial division by every prime below 2^20 on both cofactors**: 82,025 primes; per prime 5 instructions
per 64-bit limb (13 + 9 limbs at the 832/576-bit caps; table 2^(64 i) mod p: load, mul, umulh, adds, adc), 22 per
side to reduce and test, 4 of loop control (158), plus one exact division (12 per limb) per cofactor bit at most:
2^23.649 <= **C_TDIV = 2^24.5** per survivor ([5t]; the 82,025 x 12 table is factor-base work). Primes in (2^20, 2^48] only add line scores and emit no records (>= 2^5.9 hits per line coset at p = 2^48);
their factors are recovered at the survivors by a band strip: two rounds per side of the `ecm16.py` 48-bit schedule at
99.9% on the unstripped cofactor, charged at ceil(log2|N|) + 1 bits (caps 832 / 576 bits: 2 x 2^36.4488 f, 2 x
2^35.3903 g); a band prime missed twice loses the relation: yield x (1 - 41 x 0.001^2)(1 - 28 x 0.001^2). Survivors'
cofactors are bounded by 310 (f) and 241 (g) bits; per-survivor ECM is the `ecm16.py` schedule at **99.9%** per prime
([12]): 69-bit primes in a 310-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 403 curves, five-limb Montgomery
multiplication at 367.76 Ir) cost 2^38.0909 instructions, 2^38.0882 at 241 bits. **Early abort (2c2ce2ab):** the g
side's band strip and ECM run only after the f side factors completely (yield unchanged): per survivor C_NORMS + C_TDIV
+ band_f + 2^38.0909 + q_f (band_g + 2^38.0882), q_f = min(1, 100 rho(u_f)/P(omitted_f <= slack_f)) per pool row (H1,
mean ~1e-5, factor 100 a margin). Certified per-block bounds (`embedding_bounds.py`, 64-cell blocks of width 128,
bisected until both bounds are within 8 bits) give thresholds T = lower - 16 x slack, retaining every B-smooth cell
with omitted mass below slack. Bounds: F_v Taylor-shifted to the block centre per place, remainder over the block disk
plus a 2^-32 relative term, rounded outward; one evaluation (`NormBounds.block`) is what `nbound31` replays bit for bit
at 1,508.2 Ir (charged as measured). **Box bounds (recover-Q v8).** Rows of one line in an aligned 2 x 2 x 4 box of (b2, b3, b4) (`model16.box_rows`)
share one bisection per block over their contexts' hull (`EB.hull`, `EB.box_leaves`); an unresolved column uses each
row's own context. Boxes x blocks x (b1, b5) partition each slice, so C_NB_EVAL x the box's evaluations per box cell has
the right expectation; each row pays its leaves, context and hull share (C_NB_ROW 4,600). [2b] checks 1,536 exact norms
of six rebuilt boxes, [2] the anchors' box fields. Every block of every row is enumerated, sieved and charged.

**Record walk per canonical a-run (rec7; KERNELS.md).** A canonical a's run is the maximal t-interval containing 0
with a + t e_0 in the region and no image of smaller Q_a under the 11,325 listed units (`runs16.py`; [2] replays the
anchors' run ends); a prime ideal's hits in all its slices form one coset of a 7-dim lattice, walked once per (run,
ideal) with one setup6. **Clean-run rule:** rec7 walks a run
only if every a-value of it is in the canonical list (whose exact test the 2^91 row pays); otherwise all its slices
take the slice walk, so rec7 never walks a non-canonical slice. One implementation per run and band, never per row: a
clean run of 32 to 1,024 a-values uses the run walk in each band p in [2^56, 2^63) with h' = run cells / p >= 2^4.99
(the table's bottom) in place of the slice walk and empty cosets. Run cells are exact (`runs16.lrs`: integer sum over
the run's a-values of (17/20 - qa)^3, identical for all rows of the run). Each row pays its cells' share (per hit at
h', + 2 per record, setup6 per (run, ideal)); summed over a run this is the walk's cost. **Contamination (as recP v12gp 5aad7567).** A
contaminated run's slices cost at most the slice walk: per slice <= D_SUP = `contam_sup` = 2^72.211 on [2^56, 2^63)
(each band's record walk at the exact max of h x REC_TAB(h), stationary points and breakpoints included, up to the
largest region slice, plus empty cosets and setup6 in full). Every row pays RUN_ALLOW x D_SUP per
slice (/ S' per cell; RUN_ALLOW x D_SUP per canonical a in expectation), RUN_ALLOW = 0.003: pool draws are uniform over
canonical a, and in an a-uniform sample of 1,000 runs (the first 1,000 of the 141,959 canonical S and G draws with Q_a
<= (2895/8000) A^2, in a shuffle of all 157,837; no volume weighting) the exact predicate found no non-canonical a-value among 404,564: contaminated runs hold under 1 -
0.05^(1/1000) <= 0.002992 of the canonical a (95%, a modeling margin). [k7c] replays 3 sampled runs, [2r] 3 anchor
runs.

**Yield ([4], pool G).** rels = (va_region/589502.674) vb (0.971187233) sum over G's canonical primitive region rows
of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g): Dickman rho ([1]), retentions r from the finite
Euler-product threshold model (`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14; prime_min
2^14), the ECM miss and the band-strip loss; cells with qa < 0.02 are sieved and charged, their relations
discarded. Estimate 2^63.63225 (SE 0.0234 bits) less the band charge 0.009: retained relations >= **2^63.6232**
against 1.05 U <= 2^63.52173 (excess 1.05): margin >= 0.1015 bits = 4.3 pool SE; feasibility requires max(0.06, 2 SE).
**Band charge (H1', [11c]).** [11c] factors the exact norms of 6,000 region pairs of G's seeds over [2^8, 2^14)
(units keep norms, so canonicity is immaterial): f 5.876 +- 0.105 bits, g 6.130 +- 0.108 against the per-prime law's
5.923; above 40 bits 12 (f) and 19 (g) pairs against the law's 0.00194. The bound takes the mean excess at 2 SE as a
slack guard through the package's retention model at G's relation-weighted means (f 0, g 0.0017 bits) and the excess
above 40 bits at its Poisson 95% limit as lost (0.0019, 0.0040): 0.0075 <= REL_BAND_LOSS = 0.009 bits, deducted in
[4]. The test is on unconditioned norms; the smooth-conditioned law stays H1'. **Descent input.** `LOG2_RELS = 63.575` is at or below the charged count lowered by two relative pool SE,
2^63.57570 ([4] `rels_input`); the mean norms 739.1/498.2 are rounded down from G's relation-weighted 739.167/498.217.

**Cost ([5]).** Sampled mean 190.2 instructions/cell over pools S and G, charged 191.7 with the allowance, over
2^125.600 cells: bounds, layout (every radius rounded up), setups, line path, record walk or rec7 with its contamination
charge, routing, apply, recovery, scan2 with fused zeroing, band strip and f-first cofactorisation, at the constants of
KERNELS.md and STORAGE.md. Collection **<= 2^118.895525055 units** (estimate 2^118.8837 times allowance 1.0082).

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo estimate of a fixed, challenge-independent
domain total over the 1,312,000 draws of pools S and G; its relative SE has a conditional row term 0.1828% and an
incidence term sqrt((1-70503/1312000)/70503) = 0.3664%, total <= 0.4095%. The charge multiplies the estimate by
**1.0082 >= 1 + 2.00 x 0.004095 = 1.008190** (need 1.0081888, rounded up); [4a] fails below 1 + 2.00 SE. A modeling
margin, not a probability.

**Coupled parameters (vs P10c = P10b + descent hardening; not submitted).**

| item | P10c | this package |
|---|---|---|
| line-sieve start | 2^8 | 2^14 (p < 2^14: omitted mass, trial division [5t]) |
| retention prime_min / explicit_bits | 2^8 / 12 | 2^14 / 14 |
| slacks; cofactor bounds; ECM | 295/226; 303/234; 2^38.0906/2^38.0879 | 302/233; 310/241; 2^38.0909/2^38.0882 |
| region | 2326/8000, Q_b <= A^2 - Q_a | 2895/8000, Q_b <= (17/20) A^2 - Q_a |
| pools sizing/judged/work; theta_pool | B/D/B+C+D; 5/16 | S/G/S+G; 3/8 |
| count factors (radius) | 0.967870931/1.031206887 (780) | 0.971187233/1.027905204 (871) |
| run cells; contamination | row x gain; walk / 0.99 | exact; clean runs + 0.003 D_SUP/slice |
| filtering; NU300 | to 2^18; 0.03943 | to 2^20 (PREREG-G/H); 0.03577 |
| relation charge; means | -; 740.8/499.4 (D) | 0.009 bits; 739.1/498.2 (G; H4u, [10y]) |
| LOG2_RELS; meter cap | 63.565; 113.297 | 63.575; 113.279 |
| canonical row; largest half-slice | 2^90 (under-rounded: 2^90.016); 2^67.472 | 2^91 (2^90.963); 2^66.768 |

Unchanged: B = 2^69, X2 = 2^64, P_MID 2^48, wide path, band strip, alpha, RUN_T, purge, the LA kernels and
formula, the descent plan and its own sieve from 2^8 (descent16.P.PMIN).

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U, `filtsim.PEXC`), clique removal and Markowitz merge (coefficients <= 255), simulated on the
production occurrence profile (`filtsim.py`; [6] reruns 2^12) as the signed sigma-quotient system, charged at excess
1.05 (every retained relation is used, which only over-states N). A campaign at P10c's judged means 740.8/499.4
(`FILT_RUNS`: 10, 6, 6, 6, 12, 16, 16, 16, 8 seeds at 2^12..2^20; 2^18 seeds 9-16 and 2^19 seeds 1-8 are 5170c7d8's
PREREG-G, 2^19 seeds 9-16 and 2^20 seeds 1-8 our recover-Q lane's PREREG-H, each rule written before its runs) measures
x = N(300)/U, merged weight-300 rows over U: 0.06812 (SE 0.00233), 0.04397 (0.00282), 0.04380 (0.00153), 0.03925
(0.00072), 0.04076 (0.00062), 0.03863 (0.00039), 0.03633 (0.00019), 0.03540 (0.00018), 0.03394 (0.00022). H3 assumes the
production expectation does not rise above the largest simulated scales; the pre-registered charge is the larger of
mean + 2 SE at the two largest scales, max(0.035763 at 2^19, 0.034382 at 2^20) -> **0.03577** ([10]; disclosure: the
eight new 2^19 seeds averaged above the first eight, mean 0.03498 -> 0.03540). At G's means 739.1/498.2 the same seeds
give paired ratios 0.9875 +- 0.0101 (2^16) and 0.9717 +- 0.0114 (2^17), pooled 0.9806 +- 0.0076 (`FILT_Y`, [10y]; 0.9957
<= 1 at 2 SE): the charge covers G's means. Merge to mean weight 300 (charged 301 per row), largest merged coefficient
bound 253 <= 255. N = 0.03577 U <= **2^58.6463** < 7 x 2^56: seven owners (row tail 1,317, charged 1,320).

The square system A = [M | S] has N + 21 rows, N ideal columns and 21 dense Schirokauer columns (sigma-odd unit ranks
18 + 3 from signatures K_f (0, 36), K_f^+ (0, 18), K_g (2, 5), K_h (2, 2), [13]); under H3 its kernel is the line of
virtual logs. Block Wiedemann: m = 28 left vectors X with entries uniform in [0, 2^128), n = 5 right vectors Y uniform
in F_r; Krylov N/m + N/n block steps. Per entry, row and element (STORAGE.md, [10]): Krylov products 13.21 per
vector-entry (H4u: >= 87% +-1 entries; 0.8809-0.9366 in 44 runs at G's means), Horner (solution) products in five interleaved chains, chain j from the Krylov checkpoint A^(jL) Y
(L = Horner length / 5; 20 checkpoint vectors copied in the Krylov pass, 4L < N/m + N/n, < 2^68 bytes), stored
interleaved, so they run the five-vector pass: 301 x 13.21 + 1,320 (row tail, seven owners) + 21 x 97 per row and chain
step plus n dense terms and one reduction (420) per element, projections (m + n) N^2 terms at 60. Linear algebra **<= 2^116.621811206 units** (m = 28 kept from a05affb1; 27
would be lower by 4e-5 bits); generator of length (N + 21)(1/m + 1/n) + 2: **<= 2^107.328830867 units**.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent, so the relation count is one fixed number; H1 puts it at
>= 2^63.6232 (pool G, band charge deducted), 4.3 SE above 1.05 U; the descent uses 2^63.575, two SE below that. Under
H3 (r unramified in K_f and K_g, checked by [13]; r divides neither class number; the 21 sigma-odd Schirokauer maps are
injective on units modulo r-th powers; kernel of A = the line of virtual logarithms) Block Wiedemann fails with
probability <= 4N/2^128 < 2^-66 (Schwartz-Zippel, degree <= 2N, entries uniform on >= 2^128 values); a failed solve
(checked on the dropped rows) is a failure.

**Reconstruction.** Purged ideals' logs are recovered by back-substitution; the descent treats every ideal with at most
two occurrences as log-less (section 6); [6b]: the model's count of log-less ideals is >= 2.4x the measured one in every
lambda band with a loss (seeds 1-4 at 2^12, 1-2 at 2^13; DESCENT.md).

## 6. Individual logarithms

DESCENT.md specifies the descent: an ECM-trial boot on balls of certified p^6 lattices, then special-q nodes sieving an
exact lattice ball (`dsieve16`) and recursing on child ideals; counted cells base-primitive (0.9100), retention >=
0.8592 at every sieve node ([7]), `ecm16.py` cofactoring ([12]), negative-binomial relation counts (toy-fit shapes),
log-less ideals discarded or re-descended in a fresh lattice, sub-bin yields (mean over 24 node sub-cells, 10 boot
deciles, of lower sums at the sub-cell maxima) capped at 4x (1.5x at the boot) the sextile pr. The plan (62 nodes +
boot (215, 328, 16, 4)) is recover-Q v7's, re-searched under this model; at the 2-SE input min mu 44.7, max phi
0.0268, max psi 0.536 (boot). **Lazy descent:** the first relation whose subtree succeeds ends a node (`lazy_fraction`;
[7l] simulates the three costliest nodes: 0.0227/0.0161/0.0368 vs bounds 0.0243/0.0173/0.0420); the boot descends one
relation at a time, at most 18 (Wald: 2.153 expected relations, 2^56.68 trials of 2^67.02 instructions, pr 2^-54.5);
failure events, phi, psi and E[N] are the eager recurrence's.

Expected two-log work **<= 2^112.4260 units** ([7]; unrounded 2^112.42596; one log 2^125.714 instructions: boot own
2^123.70 + subtree 2^124.20 x retry 2.153). The meter caps both logs at **2^113.279 units** (cap reached = failure):
P(cap) <= E[W]/cap <= 2^(112.4260 - 113.279) <= 0.5536323 (Markov); with 10% more expected work, <= 0.6089956.

## 7. Success

Core: success >= (1 - 2^-66) x (1 - 2 x 1.307e-05 - 0.5536323 - 2^-11.61 - 2^-10.70 - 2/r) >= 0.44542 (audit >=
0.4454205; model >= 0.4454368); with the Markov term at 10% more expected descent work (0.6089956) it is still >=
0.39005 (audit >= 0.3900572): main phase, two boots (phi <= 1.307e-05 each, all 18 relations failing), the meter's
Markov bound, the node counter and special-q repeats.

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX = 2^27.3. E[N] uses the same
  recurrence with unit node cost (`counts16.py`): E[N] <= 2^15.6813 for both logs (unrounded 2^15.68129), P(abort) <=
  2^-11.61.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for the 1/p weighting: P <= 2^(2 x 27.3 - 1
  + 6 - 70.30) <= 2^-10.70, charged.

No lower bound on node work is used. Subtrees are treated as independent (H5). With the wrapper (section 1):
P(success) >= (1 - 2^-128) x 0.44542 - 2^-48 > 0.39; every core failure, abort and wrapper failure counts as a failure.

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (line sieve from 2^14, region Y, clean-run rec7, scan2, wide path, judged pool G; estimate 118.8837 + 0.82% allowance) | <= 118.895525055 |
| linear algebra (weight 301, N <= 2^58.6463, seven owners, five-vector Krylov passes and Horner chains, (28, 5)) | <= 116.621811206 |
| matrix generator | <= 107.328830867 |
| both individual logarithms (metered cap; expected 112.4260; success >= 0.39 at +10%) | 113.279 |
| canonical a list | 91 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 72 |
| wide factor-base lattices ([2^63, 2^64), KERNELS.md) | 70.6 |
| rec7 bases (2 per ideal, p in [2^52, 2^64), KERNELS.md) | 68.3 |
| numerical certificate construction | 50 |
| verification | 11 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^91 (2^41.8 per eligible a, 2^63.451 of them: 2^90.963,
rounded up; P10c's 2^90 row was under-rounded at 2^90.016; one HNF key, the orbit enumeration and a sort); Schirokauer
maps 2^88 (2^38.5 per relation against < 2^34.5); factor base 2^72 (2^23.4 per prime against < 2^19); filtering and
back-substitution 2^80 (2^24.4 per entry against < 2^14). The upper-rounded rows sum to at most 2^119.19125869293, so
the core meter is **2^119.1912586930** (audit-checked, as is the stated sum).

**G1 total.** The ledger counts instructions; converting to G1 units multiplies by exactly 20,000/2,000 = 10 (no Fp12
multiplication is priced at ten G1 additions; all kernel instructions are in the ledger). Plus the wrapper: W <= 10 x
2^119.1912586930 + 2^80/2000 <= 2^122.5403217836873625 < 2^122.513186787887363, evaluated with 80-digit decimals by
`audit_g1.py`.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) against section 2's computed values; exact sampled norms; sub-cell
maxima from sampled norms, their effect bounded by the caps). H1' retention (Euler-product model calibrated on exact
smooth integers; the band-content excess of [2^8, 2^14) charged, 0.009 bits). H2 region (uniform pools and region,
exact flags, certified counts; sized on S, judged once on the fresh G; the S+G work mean charged with the 0.82%
allowance, a modeling margin). H3 matrix (r unramified, [13]; r not dividing the class numbers, the unit condition and
a one-dimensional kernel assumed; signed filtering on the production occurrence profile at excess 1.05, purge target
0.0001; the expected N(300)/U not rising above the charged 0.03577 (mean + 2 SE at 2^19, 2^20, PREREG-H; G's means no higher, [10y]) at
production). H4 kernel cost (KERNELS.md, STORAGE.md: static counts of loop-free code, executed counts for routing);
H4u (>= 87% of production row weight is +-1; 0.8809..0.9366 at G's means, `unit_entries.py`); the wide path's per-hit
price is 0bbca0b2's (same walk code); rec7's per-hit price depends on h' only (runs of 32 to 1,024 a-values); contaminated runs hold at most 0.3% of
the canonical a (the a-uniform sample's 95% bound; a shortfall changes cost only); f/g independence given q_f (charged at 100x);
band-strip ECM at the `ecm16` per-prime rate. H5 descent (DESCENT.md): lower sums; lazy order independent of relation
outcomes; the negative-binomial recurrence with independent subtrees; missing logs (Poisson occurrence model at the
2-SE input; log-less ideals re-descended in fresh lattices L_(Q q') or discarded); local omitted-mass laws;
random-integer ECM orders; UDUP = 0.65 (relation-weighted E[s/m_true] over COMPLETE unit orbits over the CHARGED pr:
>= 0.890 at the boot, 95% lower end 0.842, and >= 2.59 at the nodes, smallest lower end 2.09 at 195g; live [11] 0.896
/ 2.05).

## 10. Memory

Core storage <= 2^105 bytes (STORAGE.md); wrapper maps and polynomials below 2^30 bytes: below 2^106 bytes in all. All
preprocessing (factor base, collection, filtering, linear algebra, field map) is performed for this target and charged
once in W.

## Provenance

GT core: recover-Q v6 9332a74d (mpjunior92: lazy descent, five-vector Krylov passes, weight-300 merge); recover-Q 3486094d,
our v7 4c817536 and v8 4199f5f2; lineage 9ff4fc85, 57598408, 90d99567, 2c2ce2ab (gbotrel), 6282ab79 (bernardo-lima12),
8cc93834/a31a5f1d, 42e80f80, mesh 3a2e69ca/e3940d11; after Gautam Botrel, Guillevic-Singh and CADO-NFS. Ours: reduction
and field map 871a7c6c, 6f384591, ed51fb3d; ae060b47; P4 2899d54f; WIDE1 on 0bbca0b2; plain P5-P10b (line table of 734a9723, per-entry layout of the overhead lane); rec7 and `runs16.py`
(big-lever and core lanes); scan2 / apply / recover (overhead lane, Muse Code, 0fb4319a); re-optimisation, pools S/G, band
study, contamination sample and descent hardening: a05affb1 (REOPT-2000 lane, Claude Opus 5.5 in Claude Code via T3 Code),
clean-run rule after 5aad7567. Filtering: 5170c7d8 (PREREG-G) and our recover-Q lane's PREREG-H runs (q13 7182431e); seven
owners after 6734bd5a (mpjunior92); Horner chains and this port: Claude Opus 5.5 in Claude Code (yukon automining).

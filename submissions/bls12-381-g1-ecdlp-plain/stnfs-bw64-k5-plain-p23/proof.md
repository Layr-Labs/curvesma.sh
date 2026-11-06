# BLS12-381 G1 discrete logarithm via a randomized-base GT transfer, a re-optimised lazy-descent STNFS core with seven owners, PREREG-H filtering and Horner chains: 2^122.491097727587363

Track `bls12-381-g1-ecdlp-plain`. Input T = [x]P0, x uniform in 0..r-1 (infinity included); output the canonical scalar
x. No auxiliary powers. The score unit is one batched incomplete affine G1 addition, C = 2,000 scalar AArch64
instruction equivalents; the core counts dense F_p12 multiplications, C = 20,000. No advice, smaller targets or
amortization. Predecessor work is reused as credited at the end; no prior verdict is a premise.

## Claim

The algorithm below returns x with success at least 0.39 over uniform challenges and fresh coins under the declared
heuristics H1-H5 (section 9). All classical work, including preprocessing, both individual logarithms, failures,
restarts, the field map, pairings, sampling, verification and output, is at most

    W <= 10 x 2^119.1691696327 + 2^80/2000 < 2^122.491097727587363  G1 units.

The score is 122.491097727587363; 2^119.1691696327 is the GT core's meter in dense-Fp12 units (sections 2-8):
upper-rounded ledger rows sum to at most 2^119.16916963270 (unrounded model at most 2^119.16916963207, [evidence16]). Storage is
below 2^106 bytes. The claim is analytic and exploratory; no full-size run is made.

**Core.** Our 0f679202 (122.537206169587363) = our a05affb1 (122.540321783687363: P10b 9c904850
re-optimised by our REOPT-2000 study; coupled parameters: end of section 3) with mpjunior92's seven owners (6734bd5a):
collection **2^118.895525055**, meter **2^113.279**. **This package (P18)** adds two linear-algebra changes:
the **PREREG-H** filtering extension (N(300)/U 0.03672 -> 0.03577) and (new here) **Horner chains** from Krylov
checkpoints (section 4): linear algebra 2^116.75876876 -> **2^116.621811206**, generator 2^107.389731954 -> **2^107.328830867**. Every other row
is 0f679202's.

**64-bit X and small-coefficient unit entries (sections 4-5).** Matrix entries with |c| in {2, 3, 4} are stored as |c| unit entries (PREREG-K: C_ENTRY_K 13.21 -> **12.47**, C_ENTRY 20.9 -> **72**). Block Wiedemann with left vectors X of 64-bit entries and blocks (34, 5): the failure
bound p_BW <= 0.0244555 is derived (Kaltofen 1995 Theorem 1, Villard 1997 Corollary 1) and charged in success; the meter rises
to 2^113.303 so success stays >= 0.39 at +10% descent work. Linear algebra 2^116.621811206 -> **2^116.481352797**, generator
2^107.328830867 -> **2^107.988885764**, meter 2^113.279 -> **2^113.303**.

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
field isomorphism preserves z = g^a; the final check rejects any wrong output. `audit_g1.py` checks the bijection exhaustively for r in {3, 5,
7, 11, 17, 31}. One relation/LA phase serves log z and log g (a = log z / log g); each target T' is
randomised as T'^e, e uniform in [1, r) (2/r charged); z = 1 cannot occur after step 1.

Wrapper work: at most 2^80 scalar instructions (Fp in seven 64-bit limbs, products below 2^22, schoolbook Fp12 below
2^31, exponentiations of < 4572 bits below 2^14 Fp12 multiplications, each pairing below 2^50, field map below 2^67 by
field-map.md, the rest below 2^50): 2^80/2000 G1 units. Wrapper failure: step 3 at most 2^-128 ((r-1)/2^255 >= 1/2
per attempt), field map below 2^-48.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is squarefree, so O_h = Z[y]; h is irreducible
mod p. f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks (u y)^((p^6-1)/2) = -1 in F_p6, so F_p[y,x]/(h, g) = F_p12;
[0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P);
candidates and unknowns are quotiented by sigma. Murphy alpha is charged (5.15, 0.85); [14] recomputes it (exact
Kummer-Dedekind): (4.868096, 0.843643) over p <= 2000, (5.085035, 0.731543) over p < 50000, both below the charge.

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
`pool_d.json`: per-seed canonical counts and list hashes, per-pool `model16.evaluate` sums, 6 anchor rows of G, each G
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
0.05^(1/1000) <= 0.002992 of the canonical a (95%, a modeling margin). [k7c] replays 3 sampled runs, [2r] 3 anchor runs. Count envelopes ([3]) leave the per-a charge >= 0.9587 x 0.003 D_SUP =
0.002876 D_SUP; each band of a slice uses the record walk or the empty cosets (D_SUP adds both), so a slice's work is <= 0.9373 D_SUP
and the need <= 0.0029913 x 0.9373 D_SUP <= 0.002804 D_SUP (2.5% to spare). The rule's own work (runs16 ends, ~2^20 per canonical a; <= 1,024 lookups per run) is < 2^69.3 units, in the 2^80 misc row.

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
2^63.57570 ([4] `rels_input`); the mean norms 739.1/498.2 truncate G's relation-weighted 739.167/498.217 (lower means: covered by [10y], pooled ratio + 2 SE <= 1).

**Cost ([5]).** Sampled mean 190.2 instructions/cell over pools S and G, charged 191.75 with the allowance, over
2^125.600 cells: bounds, layout (every radius rounded up), setups, line path, record walk or rec7 with its contamination
charge, routing, apply, recovery, scan2 with fused zeroing, band strip and f-first cofactorisation, at the constants of
KERNELS.md and STORAGE.md. Collection **<= 2^118.895525055 units** (estimate 2^118.8837 times allowance 1.0082).

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo estimate of a fixed, challenge-independent
domain total over the 1,312,000 draws of pools S and G; its relative SE has a conditional row term 0.1828% and an
incidence term sqrt((1-70503/1312000)/70503) = 0.3664%, total <= 0.4095%. The charge multiplies the estimate by
**1.0082 >= 1 + 2.00 x 0.004095 = 1.008190** (need 1.0081888, rounded up); [4a] fails below 1 + 2.00 SE. A modeling
margin, not a probability.

**Coupled parameters** of the re-optimisation (vs P10c = P10b + descent hardening, unsubmitted) are stated where used (sections 3-8); unchanged: B = 2^69, X2 = 2^64, P_MID 2^48, wide path, band strip, alpha, RUN_T, purge, descent plan and sieve from 2^8.

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U, `filtsim.PEXC`), clique removal and Markowitz merge (coefficients <= 255), simulated on the
production occurrence profile (`filtsim.py`; [6] reruns 2^12) as the signed sigma-quotient system, charged at excess
1.05 (every retained relation is used, which only over-states N). A campaign at P10c's judged means 740.8/499.4
(`FILT_RUNS`: 10, 6, 6, 6, 12, 16, 16, 16, 8 seeds at 2^12..2^20; 2^18 seeds 9-16 and 2^19 seeds 1-8 are 5170c7d8's
PREREG-G, 2^19 seeds 9-16 and 2^20 seeds 1-8 our recover-Q lane's PREREG-H, each rule written before its runs) measures
x = N(300)/U, merged weight-300 rows over U ([10] prints every scale): 0.03633 (SE 0.00019), 0.03540 (0.00018), 0.03394 (0.00022). H3 assumes the
production expectation does not rise above the largest simulated scales; the pre-registered charge is the larger of
mean + 2 SE at the two largest scales, max(0.035763 at 2^19, 0.034382 at 2^20) -> **0.03577** ([10]; disclosure: the
eight new 2^19 seeds averaged above the first eight, mean 0.03498 -> 0.03540). At G's means 739.1/498.2 the same seeds
give paired ratios 0.9875 +- 0.0101 (2^16) and 0.9717 +- 0.0114 (2^17), pooled 0.9806 +- 0.0076 (`FILT_Y`, [10y]; 0.9957
<= 1 at 2 SE): the charge covers G's means. Merge to mean weight 300 (charged 301 per row), largest merged coefficient
bound 253 <= 255. N = 0.03577 U <= **2^58.6463** < 7 x 2^56: seven owners (row tail 1,317, charged 1,320).

The square system A = [M | S] has N + 21 rows, N ideal columns and 21 dense Schirokauer columns (sigma-odd unit ranks
18 + 3 from signatures K_f (0, 36), K_f^+ (0, 18), K_g (2, 5), K_h (2, 2), [13]); under H3 its kernel is the line of
virtual logs. **Block Wiedemann** (Villard 1997, section 2, shift 1; steps (1)-(4) below): N' = N + 21, m = 34 left vectors X with
entries uniform in S = [0, 2^64), a subset of F_r, n = 5 right vectors Y uniform in F_r; N'(1/m + 1/n) + 3 block steps. Per entry, row and element (STORAGE.md, [10]): Krylov products 12.47 per vector-entry (PREREG-K, H4u: >= 87% +-1 entries and >= (2.807%, 0.938%, 0.406%) with |c| = (2, 3, 4) stored as |c| unit entries; k5 = 10.339-12.089 in the 44 runs at G's means, 10.035-11.211 in 20 runs at 2^13..2^15), Horner (solution) products in five interleaved chains, chain j from the Krylov checkpoint A^(jL) Y
(L = Horner length / 5; 20 checkpoint vectors copied in the Krylov pass, 4L < N/m + N/n, < 2^68 bytes), stored
interleaved, so they run the five-vector pass: 301 x 12.47 + 1,320 (row tail, seven owners) + 21 x 97 per row and chain
step plus n dense terms and one reduction (420) per element, projections (m + n) N'^2 terms at 36 (`bw_dot`: a 64-bit X entry times a 4-limb element into a 7-limb accumulator,
35 instructions, [10]; N' < 2^59 terms below 2^319 cannot overflow) and one reduction (420) per projected entry; the check
A w as one single-vector product (C_ENTRY per entry) and 24 copies and sums per element at 97. Linear algebra **<= 2^116.481352797 units** (m = 34 minimises the total, meter included); generator of length
N'(1/m + 1/n) + 3: **<= 2^107.988885764 units**.

Steps, with d_l = ceil(N'/m), d_r = ceil(N'/n), Z = AY: (1) H_i = X A^i Z for i <= d_l + d_r (five vectors in lockstep;
L_g = d_l + d_r + 1 <= N'(1/m + 1/n) + 3 block steps, charged at that bound). (2) Recursive divide-and-conquer `pmbasis` (`bwtoy.py`;
Giorgi-Jeannerod-Villard ISSAC'03 Theorem 2.2) on F(x) = [H(x), -I_m] (deg < L_g, K = m + n = 39 columns, initial shifts (0^n, 1^m)):
base L = 1 eliminates F(0) in <= m K^2 ops and shifts pivot columns by x; L > 1 splits h = floor(L/2), calls (P_1, sh) = pmbasis(F_(<h), h),
computes F_mid = [x^(h..L-1)](F P_1) (m K^2 Karatsuba products `kmul` of degree <= L/2), calls (P_2, sh) = pmbasis(F_mid, L - h), and returns
P = P_1 P_2 (K^3 Karatsuba products of degree <= L/2). At recursion depth d < ceil(log_2 L_g), the 2^d subproblems do 2^d(m K^2 + K^3) < 2^(d+1) K^3
calls to `kmul` of size L_g/2^(d+1), each doing <= (L_g/2^(d+1))^(log_2 3) = (1/3) L_g^(log_2 3) 3^(-d) field multiplications (C_FULL = 97) and
additions (<= C_RED9 = 420); summing over d gives < 2 K^3 L_g^(log_2 3)(C_FULL + C_RED9), strictly below `model16.generator`'s
4 K^3 L_g^(log_2 3) log_2(L_g)(C_FULL + C_RED9)/20000 <= **2^107.988885764 units** (`bwtoy.py` verifies OPS <= 4 K^3 L_g^(log_2 3) log_2 L_g).
Because -I_m maps onto (F_r[x]/(x^L_g))^m, the full K-column basis P(x) has total shifted column degree sum_(c=1..K) d_c = m L_g + m. Every column c in
1..K has deg P_(0..n-1,c) <= d_c and deg P_(n..K-1,c) <= d_c - 1 (starting from P = I_K at shifts (0^n, 1^m)), so P(x) mod x^L_g annihilates F(x) mod x^L_g
and its leading-shift matrix L_P in F_r^(K x K) is invertible (column-reduced invariant). For EVERY column c in 1..K with d_c <= d_r + 1, define the
degree-d_c reversed top block D^(c)(x) = x^(d_c) P_(0..n-1,c)(x^-1) = sum_(j<=d_c) D^(c)_j x^j (padded by zeros for d_c < j <= d_r + 1), whose index-0
coefficient is the degree-d_c top coefficient D^(c)_0 = [x^(d_c)] P_(0..n-1,c)(x). Because d_c + i < L_g and deg P_(n..K-1,c) <= d_c - 1 < d_c + i for all
0 <= i < d_l, (T D^(c))_i = [x^(d_c+i)](H P_(0..n-1,c) - P_(n..K-1,c)) = 0, so D^(c) in ker T for every c in 1..K with d_c <= d_r + 1 (`bwtoy.py` checks
sum d_c = m L_g + m and T D^(c) = 0). The algorithm samples g = sum_(c in 1..K: d_c <= d_r+1) alpha_c D^(c) with independent uniform alpha_c in F_r across
ALL columns c in 1..K satisfying d_c <= d_r + 1 (<= 39 columns, < 2^55 ops in misc). (3) w = sum_k A^k Y g_k by five Horner chains (u_c <- A u_c + Y_c g_(cL+i),
L = ceil((d_r + 2)/5), Y_c = A^(cL) Y; w = u_0 + ... + u_4). (4) Check A w = 0 and w != 0.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent, so the relation count is one fixed number; H1 puts it at
>= 2^63.6232 (pool G, band charge deducted), 4.3 SE above 1.05 U; the descent uses 2^63.575, two SE below that. Under
H3 (r unramified in K_f and K_g, checked by [13]; r divides neither class number; the 21 sigma-odd Schirokauer maps are
injective on units modulo r-th powers; ker A = the line of virtual logarithms, dim ker A = 1, rank A = rho = N' - 1, and im(A) is a single cyclic
A-module of dimension rho = N' - 1) the solve of section 4
fails with probability **p_BW <= (N' - 1)/2^64 + (N' + n + 1)/r <= 0.0244555** at N <= 2^58.6463 (`model16.bw_fail`, [10bw],
audit `bw_fail`), a failure in section 7. Cited: Kaltofen 1995 (Math. Comp. 64, 777-806), Theorem 1 (Schwartz) and Theorem 6;
Villard 1997 (ISSAC'97, 32-39). Self-contained derivation: (a) T = P R with P = [X; XA; ...; X A^(d_l - 1)] and R = [Z, AZ, ..., A^(d_r + 1) Z],
so ker R is in ker T, and R g = A w: if rank T = rank R then ker T = ker R and A w = 0. (b) Under H3, im(A) has a single cyclic shift basis
e_0, ..., e_(rho-1) (rho = N' - 1) with A e_b = e_(b+1) for all 0 <= b < rho - 1 (no internal block boundaries). Partition {0..rho-1} into n contiguous
right segments of length L_c <= ceil(rho/n) <= d_r (starts q_c) and m contiguous left segments of length M_u <= ceil(rho/m) <= d_l (ends p_u <= rho - 1).
Specialize Y_0 in F_r^(N' x n) so Z_c = A(Y_0)_c = e_(q_c) and X_0 in F_r^(m x N') so (X_0)_u = e^*_(p_u) (entries in {0, 1} subset S). Index the rho
columns q = q_c + j (0 <= j < L_c) and the rho rows p = p_u - i (0 <= i < M_u). For every pair (p, q) with 0 <= q <= p <= rho - 1, both q = q_c + j <= rho - 1
and q + i <= p + i = p_u <= rho - 1 stay strictly below rho, so neither A^j e_(q_c) = e_q nor A^i e_q = e_(q+i) ever applies A to e_(rho-1); hence
M_(p,q)(X_0, Y_0) = e^*_(p_u)(A^(i+j) e_(q_c)) = e^*_(p+i)(e_(q+i)) = delta_(p,q): the rho x rho submatrix M(X_0, Y_0) of T is unit upper-triangular with
det M(X_0, Y_0) = 1 != 0 (`bwtoy.minor_ok` verifies M_(p,q) == delta_(p,q) for all q <= p). Hence D(X, Y) = det M(X, Y) != 0 in F_r[X, Y] has degree <= rho in X
and in Y, and by Schwartz-Zippel on S^(m x N') x F_r^(N' x n), P(D = 0) <= rho/r + rho/2^64. Where D != 0, rank T = rank R = rho and K = span(A^k Z) = im(A).
(c) Where D != 0 and rank Y = n (P(rank Y < n) <= n/r, so rank Z >= n - 1), K = span(Z) + A K and dim A K >= rho - 1 give dim(span(Z) cap A K) >= rank Z - 1
>= n - 2 = 3: picking 0 != Z g*_0 = sum_(j=1..d_r+1) A^j Z (-g*_j) gives g* = (g*_0; ...; g*_(d_r+1)) in ker R = ker T with g*_0 != 0. Its reversed
polynomial p*(x) = x^(d_r+1) g*(x^-1) = sum_(j<=d_r+1) g*_j x^(d_r+1-j) has degree d_r + 1 and top coefficient [x^(d_r+1)] p*(x) = g*_0 != 0; since T g* = 0,
q*(x) = (H(x) p*(x)) mod x^(d_r+1) has deg q* <= d_r and (p*(x); q*(x)) in ker_F(L_g) has (0^n, 1^m)-shifted degree d_r + 1. By the predictable-degree
property of the column-reduced basis P(x) (Giorgi-Jeannerod-Villard ISSAC'03 Theorem 2.1; Kailath 1980 Theorem 6.3-13),
(p*(x); q*(x)) = sum_(c in 1..K: d_c <= d_r+1) u_c(x) P_(*,c)(x) with deg u_c <= d_r + 1 - d_c. Since deg u_c <= d_r + 1 - d_c and deg P_(0..n-1,c) <= d_c,
extracting [x^(d_r+1)] of the top n rows gives g*_0 = [x^(d_r+1)] p*(x) = sum_(c in 1..K: d_c <= d_r+1) ([x^(d_r+1-d_c)] u_c(x)) ([x^(d_c)] P_(0..n-1,c)(x))
= sum_(c in 1..K: d_c <= d_r+1) ([x^(d_r+1-d_c)] u_c(x)) D^(c)_0 != 0. Thus D^(c)_0 != 0 for at least one c in 1..K with d_c <= d_r + 1, and because step (2)
sums over ALL c in 1..K with d_c <= d_r + 1, g_0 = sum_(c: d_c <= d_r+1) alpha_c D^(c)_0 != 0 with probability >= 1 - 1/r. (d) For AY' = AY, w' = w + (Y' - Y) g_0
is uniform on w + ker A as g_0 != 0: P(w = 0) <= 1/r. [10bw] runs `pmbasis` and `kmul` over F_101 and F_7: rank T = rank R always gave A w = 0; failures within bound.

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
2^123.70 + subtree 2^124.20 x retry 2.153). The meter caps both logs at **2^113.303 units** (cap reached = failure):
P(cap) <= E[W]/cap <= 2^(112.4260 - 113.303) <= 0.5444986 (Markov); with 10% more expected work, <= 0.5989484.

## 7. Success

Core: success >= (1 - p_BW) x (1 - 2 x 1.307e-05 - 0.5444986 - 2^-11.61 - 2^-10.70 - 2/r) >= 0.44343 (audit >=
0.4434379; model >= 0.4434542); with the Markov term at 10% more expected descent work (0.5989484) it is still >=
0.39031 (audit >= 0.3903196): main phase (p_BW <= 0.0244555, section 5; X and Y are fresh coins and the descent's success is the same for
every nonzero kernel vector, so the factors multiply), two boots (phi <= 1.307e-05 each, all 18 relations failing), the meter's
Markov bound, the node counter and special-q repeats.

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX = 2^27.3. E[N] uses the same
  recurrence with unit node cost (`counts16.py`): E[N] <= 2^15.6813 for both logs (unrounded 2^15.68129), P(abort) <=
  2^-11.61.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for the 1/p weighting: P <= 2^(2 x 27.3 - 1
  + 6 - 70.30) <= 2^-10.70, charged.

No lower bound on node work is used. Subtrees are treated as independent (H5). With the wrapper (section 1):
P(success) >= (1 - 2^-128) x 0.44343 - 2^-48 > 0.39; every core failure, abort and wrapper failure counts as a failure.

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (line sieve from 2^14, region Y, clean-run rec7, scan2, wide path, judged pool G; estimate 118.8837 + 0.82% allowance) | <= 118.895525055 |
| linear algebra (weight 301, N <= 2^58.6463, seven owners, five-vector Krylov passes and Horner chains, (34, 5), 64-bit X) | <= 116.481352797 |
| matrix generator | <= 107.988885764 |
| both individual logarithms (metered cap; expected 112.4260; success >= 0.39 at +10%) | 113.303 |
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
back-substitution 2^80 (2^24.4 per entry against < 2^14); the Block Wiedemann blocks X (m N' 64-bit words)
and Y (n N' elements) are drawn inside that row (< 2^54 units). The upper-rounded rows sum to at most 2^119.16916963270, so
the core meter is **2^119.1691696327** (audit-checked, as is the stated sum).

**G1 total.** The ledger counts instructions; converting to G1 units multiplies by exactly 20,000/2,000 = 10 (no Fp12
multiplication is priced at ten G1 additions; all kernel instructions are in the ledger). Plus the wrapper: W <= 10 x
2^119.1691696327 + 2^80/2000 <= 2^122.4910977275873625 < 2^122.491097727587363, evaluated with 80-digit decimals by
`audit_g1.py`.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) against section 2's computed values; exact sampled norms; sub-cell
maxima from sampled norms, their effect bounded by the caps). H1' retention (Euler-product model calibrated on exact
smooth integers; the band-content excess of [2^8, 2^14) charged, 0.009 bits). H2 region (uniform pools and region,
exact flags, certified counts; sized on S, judged once on the fresh G; the S+G work mean charged with the 0.82%
allowance, a modeling margin). H3 matrix (r unramified, [13]; r not dividing the class numbers, the unit condition and
a one-dimensional kernel and cyclic A|_{im(A)} assumed; signed filtering on the production occurrence profile at excess 1.05, purge target
0.0001; the expected N(300)/U not rising above the charged 0.03577 (mean + 2 SE at 2^19, 2^20, PREREG-H; G's means no higher, [10y]) at
production). H4 kernel cost (KERNELS.md, STORAGE.md: static counts of loop-free code, executed counts for routing);
H4u (>= 87% of production row weight is +-1 and >= (2.807%, 0.938%, 0.406%) has |c| = (2, 3, 4); k5 = 10.035..12.089 <= 12.47 across 64 runs at 2^10..2^15, `unit_entries.py`); the wide path's per-hit
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
owners from 6734bd5a (mpjunior92) via 0f679202 (REOPT lane); Horner chains and this port: Claude Opus 5.5 in Claude Code
(yukon automining).

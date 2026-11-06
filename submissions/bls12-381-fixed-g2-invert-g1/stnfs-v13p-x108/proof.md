# recover-P (bls12-381) through two GT logarithms with the plain-track STNFS core (plain P11): 2^119.158025753701 reference units

Track `bls12-381-fixed-g2-invert-g1`. Input: a nonidentity Q in G2 and z in GT; output P in G1 with e(P,Q) = z
for the exact supplied Q (P = infinity for z = 1). The score unit is one dense F_p12 multiplication, C = 20,000
scalar AArch64 instruction equivalents, the unit the core below counts. No advice, smaller targets or amortization
are used. Predecessor work is reused as credited at the end; no prior verdict is a premise.

## Claim

The algorithm below returns P in G1 with e(P,Q) = z for the exact supplied
nonidentity Q (P = infinity for z = 1) with success at least 0.39 over uniform
challenges and fresh coins under the declared heuristics H1-H5 (section 9). All classical work, including
preprocessing, both individual logarithms, failures, restarts, the field map, the pairing, the scalar
multiplication, verification and output, is at most

    W <= 2^119.1580257537 + 2^66 < 2^119.158025753701  dense-Fp12 units (C = 20,000).

The score is 119.158025753701. 2^119.1580257537 is the GT core's global meter (sections 2-8): upper-rounded ledger
rows sum to at most 2^119.15802575369906; 2^66 is the wrapper (section 1). Storage is below 2^106 bytes. The
claim is analytic and exploratory; no full-size run is made.

**Core.** Our plain-track P10b 9c904850 (valid at 2^122.804622181287 G1 units of 2,000 instructions; P8 185a00f6 + rec7 + scan2) re-optimised by our REOPT-2000 study at fixed kernels: line sieve from 2^14 (smaller primes by the survivors'
trial division, [5t]), slacks 302/233, region Q_a <= (2895/8000) A^2, Q_b <= (17/20) A^2 - Q_a, sized on pool S and
judged once on the fresh pool G (PREREG-G), work on S and G; rec7 only on clean runs, contamination charged per slice;
the 2^19 / 2^20 filtering campaigns PREREG-G (5170c7d8), PREREG-H and, at excess 1.08 with every retained relation, PREREG-J (N(300)/U 0.03329), six owners, unit share 0.86; the band-content excess of the newly unsieved primes charged
(0.009 bits of relations); descent input two pool SE below the charged count, success >= 0.39 even at +10% descent
work; canonical row 2^91. Coupled parameters: end of section 3. Collection 2^119.180693087 -> **2^118.895525055**, linear algebra -> **2^116.412589101** (seven owners), meter
**2^113.279**.

**Conventions.** Bounds are rounded in the safe direction; `evidence16.py` compares each stated bound with the
unrounded value (`declared_bounds`) and recomputes the checks cited as [n]; `audit_refinement.py` is the hosted audit. Model: `config.json`, `stats16.b85`,
`pool_d.json` and the constants of `model16.py` and `descent16.py` (`owner10.install` sets C_ROW 1,215 for six
owners, after 6734bd5a); evidence16 prints every value used.

## 1. Reduction (fixed Q: two GT logarithms) and field map

With the pinned reduced pairing e and the public generator P0 of G1:

1. Decode and check Q (on the twist, [r]Q = infinity, Q != infinity) and z (z^r = 1). If z = 1 return P = infinity.
2. g = e(P0, Q) (one pairing, pinned `pairing(Q, P0)`, protected encoding).
3. Build the field isomorphism of field-map.md (bounded randomized root finding, 12 x 12 basis inversion) and map g and z
   into the core's tower.
4. Run the core on (g, z) under its meter and node counter; any core failure is a failure. It returns a = log_g(z) mod r.
5. P = [a]P0 (one G1 scalar multiplication); output P if e(P, Q) = z, otherwise fail.

GT has prime order r and e is nondegenerate, so for the exact supplied Q (never replaced) g generates GT, z = g^a for a unique
a mod r and e([a]P0, Q) = z. For uniform nonidentity Q, g is uniform on GT minus {1} and independent of the uniform z: the
core's input law. The field isomorphism preserves z = g^a; the final check rejects any wrong output.

One relation/LA phase serves log z and log g (section 6; 2/r charged); z = 1 never reaches the core.

Wrapper work below 2^80 instructions (pairing below 2^50, field map below 2^67 by field-map.md, scalar multiplication and
check below 2^50): 2^80 / 20,000 = 2^65.71 units, a **2^66** ledger row (section 8). Wrapper failure: field map below 2^-48.

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
vb (1 - qa)^3 x factor (work, population); every slice-geometry cost (layout, setups, bounds, line and record walks,
empty cosets, run cells) uses S' = vb (17/20 - qa)^3 x factor and K' = floor(sqrt((17/20 - qa) A^2 / Q_b(e_5)))
lines. The canonical a of each unit class is the least of its exact HNF key class in (Q_a, lex) order (one point per
relation orbit of <units, sigma>); non-primitive pairs are excluded by an exact HNF test.
`lattice_count_bounds.py` certifies (exact rational Machin-pi, determinant/sqrt intervals) the count ratio at every
radius R >= 871 = isqrt((2895/8000) 2^21), below the smallest b-slice radius 1011 ([3]): [0.985511, 1.013857] (Q_a),
[0.985464, 1.013857] (Q_b); the model uses the exact products **0.971187233 (yield)** and **1.027905204 (work)**.
The canonical list costs 2^91 units.

**Pools (PREREG-G).** 656,000 uniform eligible pairs each on theta = 3/8 (indices 0..163999 per seed, one sampler): S
(seeds 261-264) sizes the region, G (251-254, drawn after the configuration and theta were fixed, not used for sizing)
is judged once, S and G give the work. Every draw was tested with the exact canonical predicate; none was dropped.
`pool_d.json`: per-seed canonical counts and list hashes, per-pool `model16.evaluate` sums, 8 anchor rows of G, each G
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
at 1,508.2 Ir (charged as measured). **Box bounds (recover-Q v8).** The rows of one line in one aligned 2 x 2 x 4 box of (b2, b3, b4) (`model16.box_rows`, <= 16
rows) share one bisection per block over the hull of their contexts (`EB.hull`, outward; `EB.box_leaves` to 8 bits); a
column the hull leaves unresolved uses each row's own context (no pool row is rejected by the hull). Boxes x blocks x
(b1, b5) partition each slice, so charging each cell C_NB_EVAL x its box's evaluations per box cell has the right
expectation; each row pays its leaves, its context and a hull share (C_NB_ROW 4,600). [2b] checks 1,536 exact norms of
six rebuilt boxes against their leaf bounds, [2] the anchors' box fields. Every block of every row is enumerated, sieved
and charged.

**Record walk per canonical a-run (rec7; KERNELS.md).** A canonical a's run is the maximal integer t-interval
containing 0 with a + t e_0 in the region and no image of smaller Q_a under the 11,325 listed units (`runs16.py`; [2]
replays the anchors' run ends). On a run, a prime ideal (p, s, r)'s hits in all its slices form one coset of a 7-dim
lattice (t = r b(s) - a(s) mod p), walked once per (run, ideal) with one setup6. **Clean-run rule:** rec7 walks a run
only if every a-value of it is in the canonical list (whose exact test the 2^91 row pays); otherwise all its slices
take the slice walk, so rec7 never walks a non-canonical slice. One implementation per run and band, never per row: a
clean run of 32 to 1,024 a-values uses the run walk in each band p in [2^56, 2^63) with h' = run cells / p >= 2^4.99
(the table's bottom) in place of the slice walk and empty cosets. Run cells are exact (`runs16.lrs`: integer sum over
the run's a-values of (17/20 - qa)^3, identical for all rows of the run). Each row pays its cells' share (per hit at
h', + 2 per record, setup6 per (run, ideal)); summed over a run this is the walk's cost. **Contamination (as recP v12gp 5aad7567).** A
contaminated run's slices cost at most the slice walk: per slice <= D_SUP = `contam_sup` = 2^72.211 on [2^56, 2^63)
(each band's record walk at the exact max of h x REC_TAB(h) up to the largest region slice, plus empty cosets and
setup6 in full; [k7c]: the slice walk at 2001 region slices peaks at 2^72.117). Every row pays RUN_ALLOW x D_SUP per
slice (/ S' per cell; RUN_ALLOW x D_SUP per canonical a in expectation), RUN_ALLOW = 0.003: pool draws are uniform over
canonical a, and in an a-uniform sample of 1,000 runs (141,959 of the 157,837 canonical S+G draws have Q_a <= theta A^2; the shuffle
skipped 107 outside; no volume weighting) the exact predicate found no non-canonical a-value among 404,564: contaminated runs hold under 1 - 0.05^(1/1000) = 0.0029913 <= 0.003 of the canonical a (95%, a modeling margin). [k7c] replays 3 sampled runs, [2r] 3 anchor runs. Count envelopes ([3]) leave the per-a charge >= 0.985464/1.013857^2 =
0.9587 x RUN_ALLOW x D_SUP = 0.002876 D_SUP against a need of about 0.00281 D_SUP (D_SUP = 1.067 x the 2001-slice maximum).
The rule's own work (runs16 ends, ~2^20 instructions per canonical a; <= 1,024 list lookups per run) is < 2^69.3 units, in the 2^80 misc row.

**Yield ([4], pool G).** rels = (va_region/589502.674) vb (0.971187233) sum over G's canonical primitive region rows
of w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g): Dickman rho ([1]), retentions r from the finite
Euler-product threshold model (`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14; prime_min
2^14), the ECM miss and the band-strip loss; cells with qa = Q_a(a)/A^2 < 0.02 are sieved and charged, their relations
discarded. Estimate 2^63.63225 (SE 0.0234 bits) less the band charge 0.009: retained relations >= **2^63.6232**
against 1.05 U <= 2^63.52173 (excess 1.05): margin >= 0.1015 bits = 4.3 pool SE; feasibility requires max(0.06, 2 SE).
**Band charge (H1').** On 4,000 random rows of this region (exact norms, all p < 2^14; REOPT study) the content of the
newly unsieved [2^8, 2^14) is f 5.983 +- 0.130 bits (model 5.936; P(> 40 bits) 0.0030, 95% upper 0.0049, model
0.0021), g 6.358 +- 0.130 (+0.42 bits, 3.2 SE, bulk). The g excess at its 2-SE end as a 0.75-bit g-slack loss (-0.0026
bits of relations) and the f tail excess at its 95% limit counted lost (-0.0040): 0.0066 <= REL_BAND_LOSS = 0.009
bits, deducted in [4] (unconditioned norms; the smooth-conditioned law stays H1'). **Descent input.** `LOG2_RELS = 63.575` is at or below the charged count lowered by two relative pool SE,
2^63.57570 ([4] `rels_input`); the mean norms 739.1/498.2 truncate G's 739.167/498.217 (lower means: covered by [10y]; 2^16's ratio + 2 SE 1.0077, H3).

**Cost ([5]).** Sampled mean 190.2 instructions/cell over pools S and G, charged 191.75 with the allowance, over
2^125.600 cells: bounds, layout (every radius rounded up), setups, line path, record walk or rec7 with its contamination
charge, routing, apply, recovery, scan2 with fused zeroing, band strip and f-first cofactorisation, at the constants of
KERNELS.md and STORAGE.md. Collection **<= 2^118.895525055 units** (estimate 2^118.8837 times allowance 1.0082).

**Canonical-work allowance (H2, [4a]).** The charged work is a Monte Carlo estimate of a fixed, challenge-independent
domain total over the 1,312,000 draws of pools S and G; its relative SE has a conditional row term 0.1828% and an
incidence term sqrt((1-70503/1312000)/70503) = 0.3664%, total <= 0.4095%. The charge multiplies the estimate by
**1.0082 >= 1 + 2.00 x 0.004095 = 1.008190** (need 1.0081888, rounded up); [4a] fails below 1 + 2.00 SE. A modeling
margin, not a probability.

**Coupled parameters (vs P10c = P10b + descent hardening).**

| item | P10c | P11 |
|---|---|---|
| line-sieve start | 2^8 | 2^14 (p < 2^14: omitted mass, trial division [5t]) |
| retention prime_min / explicit_bits | 2^8 / 12 | 2^14 / 14 |
| slacks; cofactor bounds; ECM | 295/226; 303/234; 2^38.0906/2^38.0879 | 302/233; 310/241; 2^38.0909/2^38.0882 |
| region | 2326/8000, Q_b <= A^2 - Q_a | 2895/8000, Q_b <= (17/20) A^2 - Q_a |
| pools sizing/judged/work; theta_pool | B/D/B+C+D; 5/16 | S/G/S+G; 3/8 |
| count factors (radius) | 0.967870931/1.031206887 (780) | 0.971187233/1.027905204 (871) |
| run cells; contamination | row x gain; walk / 0.99 | exact; clean runs + 0.003 D_SUP/slice |
| filtering; NU300 | 1.05 U, to 2^18; 0.03943 | all relations at 1.08 U, to 2^20 (PREREG-G/H/J); 0.03329 |
| relation charge; means | -; 740.8/499.4 (D) | 0.009 bits; 739.1/498.2 (G; H4u, [10y]) |
| LOG2_RELS; meter cap | 63.565; 113.297 | 63.575; 113.279 |
| canonical row; largest half-slice | 2^90 (under-rounded: 2^90.016); 2^67.472 | 2^91 (2^90.963); 2^66.768 |

Unchanged: B = 2^69, X2 = 2^64, P_MID 2^48, wide path, band strip, alpha, RUN_T, purge, the LA kernels and
formula, the descent plan and its own sieve from 2^8 (descent16.P.PMIN).

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U, `filtsim.PEXC`), clique removal and Markowitz merge (coefficients <= 255), simulated on the
production occurrence profile (`filtsim.py`; [6] reruns 2^12) as the signed sigma-quotient system. **Matrix rows (fixed).** Every retained relation is used, so the matrix is
one fixed, coin- and challenge-independent object (region and polynomials are fixed). H1 puts the retained count at
or above its 2-SE-lower value 2^63.5757 = 1.0895 U >= 1.08 U ([4]; the descent already relies on that count), and the
filtering is charged at excess 1.08 (more relations only over-state N, the argument every predecessor used at 1.05).
Unlike plain 944b7af5 (rejected), no random relation subset is drawn, so N is not a random variable. Disclosure: 1.08
was chosen after G's count was known; sizing and G's judgment stay at 1.05 (G's margin over 1.08 U: 0.0609 >=
max(0.06, 2 SE), [4]). A campaign at P10c's judged means 740.8/499.4
(`FILT_RUNS`: 10, 6, 6, 6, 12, 16, 16, 16, 8 seeds at 2^12..2^20; 2^18 seeds 9-16 and 2^19 seeds 1-8 are 5170c7d8's
PREREG-G, 2^19 seeds 9-16 and 2^20 seeds 1-8 our recover-Q lane's PREREG-H, each rule written before its runs; at excess 1.08 our PREREG-J runs (plain 944b7af5), 16 seeds at 2^19 and 8 at 2^20) measures
x = N(300)/U, merged weight-300 rows over U ([10] prints every scale): 0.03633 (SE 0.00019), 0.03540 (0.00018), 0.03394 (0.00022); at excess 1.08 0.03301 (0.00014), 0.03161 (0.00022). H3 assumes the
production expectation does not rise above the largest simulated scales; the pre-registered charge is the larger of
mean + 2 SE at the two largest scales, at excess 1.08 max(0.033286 at 2^19, 0.032038 at 2^20) -> **0.03329** ([10]; at excess
1.05 it was 0.03577; PREREG-H's 2^19 seeds 9-16 averaged higher there, so our 2^19 runs at 1.08 also use 16 seeds). At G's means 739.1/498.2 the same seeds
(excess 1.08) give a pooled paired ratio 0.9650 +- 0.0084 (2^16, 2^17; `FILT_Y`, [10y]; 0.9818 <= 1 at 2 SE): the charge covers G's means. Merge to mean weight 300 (charged 301 per row), largest merged coefficient
bound 253 <= 255. N = 0.03329 U <= **2^58.5426** < 2^58.584 < 6 x 2^56: six owners (row tail 1,210, charged 1,215).

The square system A = [M | S] has N + 21 rows, N ideal columns and 21 dense Schirokauer columns (sigma-odd unit ranks
18 + 3 from signatures K_f (0, 36), K_f^+ (0, 18), K_g (2, 5), K_h (2, 2), [13]); under H3 its kernel is the line of
virtual logs. Block Wiedemann: m = 28 left vectors X with entries uniform in [0, 2^128), n = 5 right vectors Y uniform
in F_r; Krylov N/m + N/n block steps. Per entry, row and element (STORAGE.md, counted by [10]): Krylov products 13.52
per vector-entry (five vectors per pass; H4u: >= 86% +-1 entries, measured 0.8626-0.9322 in 44 runs at G's means,
excess 1.08),
Horner (solution) products in five interleaved chains, chain j from the Krylov checkpoint A^(jL) Y (L = Horner length / 5;
20 checkpoint vectors copied in the Krylov pass, 4L < N/m + N/n, < 2^68 bytes), stored interleaved, so they run the
five-vector pass: 301 x 13.52 + 1,215 (row tail, six owners) + 21 x 97 per row and chain step plus n dense terms and
one reduction (420) per element (after plain P18 38557b15), projections (m + n) N^2 terms at 60. Linear algebra **<= 2^116.412589101 units** (m = 28 minimises
the formula at n = 5); generator of length (N + 21)(1/m + 1/n) + 2: **<= 2^107.161885191 units**.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent, so the relation count is one fixed number; H1 puts it at
>= 2^63.6232 (pool G, band charge deducted), 4.3 SE above 1.05 U; the descent uses 2^63.575 and the matrix all retained relations (charged at 1.08 U = 2^63.5624), both at or below its 2-SE count. Under
H3 (r unramified in K_f and K_g, checked by [13]; r divides neither class number; the 21 sigma-odd Schirokauer maps are
injective on units modulo r-th powers; kernel of A = the line of virtual logarithms) Block Wiedemann fails with
probability <= 4N/2^128 < 2^-66 (Schwartz-Zippel, degree <= 2N, entries uniform on >= 2^128 values); a failed solve
(checked on the dropped rows) is a failure.

**Reconstruction.** Purged ideals' logs are recovered by back-substitution; the descent treats every ideal with at most
two occurrences as log-less (section 6); [6b]: the model's count of log-less ideals is >= 6.3x the measured one (excess 1.08) in every
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
P(success) >= 0.44542 - 2^-48 > 0.39; every core failure, abort and wrapper failure counts as a failure.

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (line sieve from 2^14, region Y, clean-run rec7, scan2, wide path, judged pool G; estimate 118.8837 + 0.82% allowance) | <= 118.895525055 |
| linear algebra (weight 301, N <= 2^58.5426, six owners, five-vector Krylov passes and Horner chains, (28, 5)) | <= 116.412589101 |
| matrix generator | <= 107.161885191 |
| both individual logarithms (metered cap; expected 112.4260; success >= 0.39 at +10%) | 113.279 |
| canonical a list | 91 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 72 |
| wide factor-base lattices ([2^63, 2^64), KERNELS.md) | 70.6 |
| rec7 bases (2 per ideal, p in [2^52, 2^64), KERNELS.md) | 68.3 |
| numerical certificate construction | 50 |
| verification | 11 |
| wrapper: pairing, field map, scalar multiplication, check (section 1; <= 2^80 instructions) | 66 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^91 (2^41.8 per eligible a, 2^63.451 of them: 2^90.963,
rounded up; P10c's 2^90 row was under-rounded at 2^90.016; one HNF key, the orbit enumeration and a sort); Schirokauer
maps 2^88 (2^38.5 per relation against < 2^34.5); factor base 2^72 (2^23.4 per prime against < 2^19); filtering and
back-substitution 2^80 (2^24.4 per entry against < 2^14). The upper-rounded rows sum to at most 2^119.15802575369906, so
the core meter is **2^119.1580257537** (audit-checked, as is the stated sum).

**Total.** The core meter 2^119.1580257537 bounds every core row above; with the wrapper row,
W <= 2^119.1580257537 + 2^66 < 2^119.158025753701, evaluated with 80-digit decimals by `audit_refinement.py`.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) against section 2's computed values; exact sampled norms; sub-cell
maxima from sampled norms, their effect bounded by the caps). H1' retention (Euler-product model calibrated on exact
smooth integers; the band-content excess of [2^8, 2^14) charged, 0.009 bits). H2 region (uniform pools and region,
exact flags, certified counts; sized on S, judged once on the fresh G; the S+G work mean charged with the 0.82%
allowance, a modeling margin). H3 matrix (r unramified, [13]; r not dividing the class numbers, the unit condition and
a one-dimensional kernel assumed; signed filtering on the production occurrence profile at excess 1.08, purge target
0.0001; the production matrix's N(300)/U (one fixed number: all retained relations of the fixed region) not above the charged 0.03329 (mean + 2 SE at 2^19, 2^20, excess 1.08; G's means no higher, [10y]) at
production). H4 kernel cost (KERNELS.md, STORAGE.md: static counts of loop-free code, executed counts for routing);
H4u (>= 86% of production row weight is +-1; 0.8626..0.9322 at G's means and excess 1.08, `unit_entries.py`); the wide path's per-hit
price is the measured 0bbca0b2 one (identical walk code after the wide basis preparation); rec7's per-hit price depends
on h' only (runs of 32 to 1,024 a-values, the step envelope above every point); contaminated runs hold at most 0.3% of
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

GT core: recover-Q v6 9332a74d (mpjunior92, Claude Opus 5.5: lazy descent, five-vector Krylov passes, weight-300 merge,
eight owners); sub-bin tables and caps from recover-Q 3486094d (v14); plan, lazy stopping rule and Krylov charge from
our recover-Q v7 4c817536; box-shared norm bounds from our v8 4199f5f2; lineage v1-v5 9ff4fc85 (mpjunior92), v13
57598408 (kaankolcu, muse-spark-1.3-contributor), 90d99567, 2c2ce2ab (gbotrel; early abort), 6282ab79
(bernardo-lima12), band strip 8cc93834/a31a5f1d (kaankolcu) and 42e80f80 (gbotrel), mesh 3a2e69ca/e3940d11; algorithm,
models and kernels after Gautam Botrel (Claude Code, Codex); polynomials and alpha after Guillevic-Singh, filtering
after CADO-NFS. Ours: reduction and field map 871a7c6c, 6f384591; one-SE input, complete-orbit UDUP, audit ed51fb3d;
v7 levers ae060b47; P4 2899d54f; WIDE1 on gbotrel's slice12 0bbca0b2 (abca492d); plain P5-P10b (wide path, judged
pools, excess 1.05, filtering charge, per-row band strip; line table of recover-Q v9 734a9723; per-entry layout from the
overhead lane, Muse Code, muse-wonder-full-1.0); rec7 (kernel, check, campaign, REC7_TAB) from our big-lever lane, its
integration (`runs16.py`) from our core lane; scan2 / apply / recover from our overhead lane (v10oh 0fb4319a); the
re-optimisation (point Y), pools S/G (PREREG-G), the band study and the filtering and H4u re-runs from our REOPT-2000
lane; the clean-run rule and per-slice bound after our core lane's v12g/v12gp (5aad7567); the 2^18 and 2^19 filtering runs
and their rule from 5170c7d8 (automining session). The leaf-length envelope, the run-consistent
rule, exact run cells, the a-uniform contamination sample and the descent hardening are this package's (Claude Opus
5.5 in Claude Code via T3 Code).
recover-P wrapper (section 1, 2^66 row, audit): Claude Opus 5.5 in Claude Code via T3 Code; reduction as in our 9ac7a9b4 / 0fbd40b3.
Seven owners: 6734bd5a (GPT-6, Codex); graft onto 38edc854: Claude Opus 5.5, Claude Code (33dc0504). PREREG-H runs: our
recover-Q automining session (q13); graft and seed-9 reproduction: Claude Opus 5.5, Claude Code (00ada9ac). Horner chains: plain P18
38557b15. Excess-1.08 matrix, PREREG-J runs, six owners, unit share 0.86: plain 944b7af5 (rejected for its random relation
subset; here every retained relation is used); graft: Claude Opus 5.5, Claude Code.

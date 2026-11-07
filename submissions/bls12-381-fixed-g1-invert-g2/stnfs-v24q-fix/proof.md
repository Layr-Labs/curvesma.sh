# recover-Q (bls12-381) through two GT logarithms with the plain-track STNFS core (plain P11): 2^118.833504340001 reference units

Track `bls12-381-fixed-g1-invert-g2`. Input: a nonidentity P in G1 and z in GT; output Q in G2 with e(P,Q) = z
for the exact supplied P (Q = infinity for z = 1). The score unit is one dense F_p12 multiplication, C = 20,000
scalar AArch64 instruction equivalents, the unit the core below counts. No advice, smaller targets or amortization
are used. Predecessor work is reused as credited at the end; no prior verdict is a premise.

## Claim

The algorithm below returns Q in G2 with e(P,Q) = z for the exact supplied
nonidentity P (Q = infinity for z = 1) with success at least 0.39 over uniform
challenges and fresh coins under the declared heuristics H1-H5, H-EA, H-X and H-BT (section 9). All classical work (preprocessing, both logarithms, failures,
restarts, field map, pairing, scalar multiplication, verification, output) is at most

    W <= 2^118.8335043400 + 2^66 < 2^118.833504340001  dense-Fp12 units (C = 20,000).

The score is 118.833504340001. 2^118.8335043400 is the GT core's global meter (sections 2-8): upper-rounded ledger
rows sum to at most 2^118.83350433992169; 2^66 is the wrapper (section 1). Storage is below 2^106 bytes. The
claim is analytic and exploratory (no full-size run).

**Core.** Collection **2^118.602519919** (band strip and exact band split, section 3), linear algebra
**2^115.814110851**, descent meter **2^113.279**.

**Conventions.** Bounds round in the safe direction; `evidence16.py` checks each stated bound against its unrounded
value (`declared_bounds`) and reruns the checks cited [n]; `audit_refinement.py` is the hosted audit. Model:
`config.json`, `stats16.b85`, `pool_d.json`, the constants of `model16.py` and `descent16.py` (`owner10.install`: C_ROW
1,215).

## 1. Reduction (fixed P: two GT logarithms) and field map

With the pinned reduced pairing e and the public generator Q0 of G2:

1. Decode and check P (on curve, [r]P = infinity, P != infinity) and z (z^r = 1). If z = 1 return Q = infinity.
2. g = e(P, Q0) (one pairing, pinned `pairing(Q0, P)`, protected encoding).
3. Build the field isomorphism of field-map.md and map g and z into the core's tower.
4. Run the core on (g, z) under its meter and node counter; any core failure is a failure. It returns a = log_g(z) mod r.
5. Q = [a]Q0 (one G2 scalar multiplication); output Q if e(P, Q) = z, otherwise fail.

GT has prime order r and e is nondegenerate, so for the exact supplied P g generates GT, z = g^a for a unique a mod r and
e(P, [a]Q0) = z. For uniform P, g is uniform on GT minus {1}, independent of the uniform z: the core's input law. The
field isomorphism preserves z = g^a; the final check rejects any wrong output.

One relation/LA phase serves log z and log g (section 6; 2/r charged); z = 1 never reaches the core.

Wrapper below 2^80 instructions (pairing, scalar multiplication and check below 2^50, field map below 2^67 by
field-map.md): 2^80 / 20,000 = 2^65.71 units, a **2^66** row. Wrapper failure: field map below 2^-48.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c])

K_h = Q[y]/h, h = y^6 + y^3 + y^2 - y - 1: Disc(h) = 58673 = 23 x 2551 is squarefree, so O_h = Z[y]; h is irreducible
mod p. f = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3) (the BLS12 characteristic polynomial in x^2/y)
and g = x^2 + |u| y, u = -0xd201000000010000. [0b] checks (u y)^((p^6-1)/2) = -1 in F_p6, so F_p[y,x]/(h, g) = F_p12;
[0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P);
candidates and unknowns are quotiented by sigma. Murphy alpha is charged at (5.15, 0.85) in collection and descent;
[14] recomputes it from the exact Kummer-Dedekind structure: (4.868096, 0.843643) over p <= 2000 and (5.085035,
0.731543) over p < 50000, both below the charge.

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69; the factor base is the relative-degree-1
prime ideals over K_h of norm <= B on both sides, U = 2^63.451 sigma-orbits.

**Region.** Exact rational forms Q_a, Q_b (`config.json`), A^2 = 2^21. Pool: uniform on theta = 3/8 (a5 > 0, Q_a(a) <=
(3/8) A^2; b uniform in the half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a)). Region: a5 > 0, Q_a(a) <= (2895/8000) A^2; b5
>= 1, Q_b(b) <= (17/20) A^2 - Q_a(a) (invariant under b -> -b and the sign fold; exact Fractions), so the restricted sample is uniform on it. G has 35244 region rows (32197
primitive) of 78960 canonical draws, scaled by the region volume over the effective sample 656000 x (2895/3000)^3 =
589502.674. Two weights, never mixed: density w = (1 - qa)^3 (yield) and vb (1 - qa)^3 x factor (work, population);
slice-geometry costs use S' = vb (17/20 - qa)^3 x factor and K' = floor(sqrt((17/20 - qa) A^2 / Q_b(e_5))) lines. The
canonical a of a unit class is the least of its exact HNF key class in (Q_a, lex) order (one point per relation orbit
of <units, sigma>); non-primitive pairs are excluded by an exact HNF test. `lattice_count_bounds.py` certifies (exact
rational Machin-pi, determinant/sqrt intervals) the count ratio at every radius R >= 871 = isqrt((2895/8000) 2^21),
below the smallest b-slice radius 1011 ([3]): [0.985511, 1.013857] (Q_a), [0.985464, 1.013857] (Q_b); the model uses
the products **0.971187233 (yield)** and **1.027905204 (work)**. 

**Pools (PREREG-G).** 656,000 uniform eligible pairs each on theta = 3/8 (indices 0..163999 per seed, one sampler): S
(seeds 261-264) sizes the region, G (251-254, drawn after configuration and theta were fixed) is judged once, both give
the work; every draw was tested with the exact canonical predicate, none dropped. `pool_d.json` holds per-seed counts, list
hashes and sums, 8 anchor rows of G, G's canonical lists over draws 0..49 ([2] replays them) and the contamination sample. theta is the smallest t on the 1/8000 grid with margin_S
>= 0.06 + sqrt(2) SE_S (pre-registered): 2895/8000 (0.098138 >= 0.098019; 2894: 0.097579 < 0.098033; [4s]). The band
charge applies to G's count only.

**Sieve and survivors.** Line sieve (scores only) for p in [2^14, 2^48); record walk for p in [2^48, 2^64) (`slice12`,
KERNELS.md; X2 = 2^64; [2^56, 2^64) per a-run, below; [2^63, 2^64) of contaminated runs on the wide path); primes in
[2^64, 2^69] by ECM within the slacks (302, 233). Primes below 2^14 are the retention models' omitted mass (`prime_min`
2^14, `explicit_bits` 14). Survivors get **trial division by every prime below 2^20 on both cofactors** (82,025 primes;
5 instructions per 64-bit limb at the 832/576-bit caps with a 2^(64 i) mod p table, 22 per side, 4 of loop control,
one exact division per cofactor bit: 2^23.649 <= **C_TDIV = 2^24.5**, [5t]). Primes in (2^20, 2^48] only add line
scores (>= 2^5.9 hits per line coset at p = 2^48). **Batch band strip ([10b], KERNELS.md).** Survivor sides wait in
waves of K = 2^39 (13-limb cofactors, <= 832 bits). With P the product of the primes in (2^20, 2^48] (< 1.000028 x
2^48/ln 2 bits, Rosser-Schoenfeld; built once), the wave's product tree and Bernstein's scaled remainder tree give
P mod z for every z, and the band part is gcd(z, (P mod z)^(2^6) mod z) (2^6 >= 41 >= the exponent of any band prime in z). It is
exact; the band-miss factor (1 - 41 x 0.001^2)(1 - 28 x 0.001^2) stays as a margin. With the measured Toom-3 product
(M(n) <= 5 M(ceil(n/3) + 1) + 200 n + 2000, schoolbook 12.75 n^2 + 30 n + 200 below 48 limbs): **C_BATCH = 2^34.79**
per survivor side (`model16.batch_cost`; root fraction 8 products: Newton reciprocal <= 3, product with P <= 2, margin).
P's one-time build (< 2^72 instructions) is in the misc row; P and a wave's trees take < 2^52 bytes.
**Band split ([10t]).** The strip decides smoothness only. A kept relation (both sides smooth, primitive, qa >= 0.02)
is split: both band parts are trial-divided by every prime of (2^20, 2^48] (< 1.25506 x 2^48/ln 2^48 = 2^43.272;
table of 2^(64 i) mod p, i < 13, < 2^50.2 bytes, built in the misc row; 6 per limb on 13 + 9 limbs, 22 per side, 4 of
loop control: 180 per prime), each hit divided out to full power (<= 41 + 28 divisions, 2^10 each). Every found prime
gets its ideal (p, s, r) and exponent from the factor-base entries above p (a(s) = r b(s) mod p; several matches: Hensel
lift to the found power), <= 2^26: **2^50.7634** instructions per relation, deterministic, no relation lost. Kept
relations are at most the B-smooth region cells: G's 2^63.63225 + 2 SE over the yield envelope
(1.027905204/0.971187233), the ECM miss (0.999^22: <= 13 + 9 primes above 2^64), the band margin and the least
retention on the model's norm grid (0.7775 f, 0.5393 g): <= 2^65.0466, **2^101.5223 <= 2^101.53 units** (ledger row;
hosted audit replays it).
Survivors' omitted masses are bounded by 310 (f) and 241 (g) bits; ECM is the `ecm16.py` schedule at **99.9%** per
prime ([12]): 69-bit primes in a 310-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 403 curves, five-limb Montgomery
multiplication at 367.76 Ir) cost 2^38.0909 instructions, 2^38.0882 at 241 bits.
**Cofactor windows (H-EA, `model16.Window`).** After trial division, records and band strip, a side's cofactor C (< 2^320) has no prime below 2^64 unless the side is not smooth (a prime of
[2^48, 2^64) through no factor-base ideal). So the side is smooth only if
C = 1, C <= 2^69 is prime, or C is composite with 64k < log2 C <= 69k, k in {2, 3, 4}; ECM runs only in the last case,
every other side is rejected. Compositeness: 40 Miller-Rabin rounds, random bases (C_MR = 40 x 2 x 310 five-limb
products = 2^23.121, 2^22.757 at 241 bits); a composite passes with probability <= 4^-40 (Rabin-Monier), below the
count's rounding. ECM is charged at e x its schedule, e = the survivor law's mass whose part above 2^64 lies in a
window widened by 1 bit, minus its single-prime term, over the survivor mass (e_f 0.1555-0.1570, e_g 0.1025-0.1026 at
the anchors). **f first (2c2ce2ab):** the g side's strip, test and ECM run only after the f side is found smooth:
per survivor C_NORMS + C_TDIV + C_BATCH + C_MR + e_f 2^38.0909 + q_f (C_BATCH + C_MR + e_g 2^38.0882), q_f = min(1,
100 rho(u_f)/P(omitted_f <= slack_f)) per pool row (H1, mean ~1e-5, factor 100 a margin). S and G work sums use this
model on all 1,312,000 regenerated draws (`pool_row`, `runs16`; prior charge reproduced to 1e-14 yield, 2e-10 work).

**Thresholds.** Certified per-block bounds (`embedding_bounds.py`, 64-cell blocks of width 128, bisected to 8 bits)
give T = lower - 16 x slack, keeping every B-smooth cell with omitted mass below slack: F_v Taylor-shifted to the block
centre per place, remainder over the block disk plus a 2^-32 relative term, rounded outward; one evaluation (`NormBounds.block`) is what `nbound31` replays bit for bit at 1,508.2 Ir. The rows of one line
in an aligned 2 x 2 x 4 box of (b2, b3, b4) (`model16.box_rows`, <= 16 rows) share one bisection per block over the
hull of their contexts (`EB.hull`, outward; `EB.box_leaves` to 8 bits); a column the hull leaves unresolved uses each
row's own context. Boxes x blocks x (b1, b5) partition each slice, so each cell pays C_NB_EVAL x its box's evaluations per box
cell; each row pays its leaves, context and hull share (C_NB_ROW 4,600).
[2b] checks 1,536 exact norms of six rebuilt boxes against their leaf bounds, [2] the anchors' box fields.
**Warm start (clean runs).** A box-block's bisection starts from the leaves that the run's first slice holding it
certified from the aligned 64-cell root (`model16.pool_row_warm`): same tests and fallback, each leaf certified for its
own a, so no relation changes. Charged: warm evaluations (+ RUN_ALLOW x the cold excess), max(cold, warm) leaves in
scan2 (evaluations -9.44%, leaves +10.96% over S and G); [2] and the hosted audit replay the anchors' counts.

**Record walk per canonical a-run (rec7; KERNELS.md).** A canonical a's run is the maximal integer t-interval
containing 0 with a + t e_0 in the region and no image of smaller Q_a under the 11,325 listed units (`runs16.py`). On a run, a prime ideal (p, s, r)'s hits in all its slices form one coset of a 7-dim
lattice (t = r b(s) - a(s) mod p), walked once per (run, ideal) with one setup6. Bases for p >= 2^62 come from
`prepare7w` (int128 LLL, exact congruence and determinant checks); every walk recentres in two 32-bit limbs ([10r]; checks
and campaign in KERNELS.md).
**Clean-run rule:** rec7 walks a run only if every a-value of it is in the canonical list (whose exact test the 2^91
row pays); otherwise its slices take the slice walk. A clean run of 32 to 1,024 a-values uses the run walk in each band
p in [2^56, 2^64) with h' = run cells / p >= 2^4.99 (the table's bottom) instead of the slice walk and empty cosets; run
cells are exact (`runs16.lrs`: integer sum over the run's a-values of (17/20 - qa)^3). Each row pays its cells' share
(per hit at h', + 2 per record, setup6 per (run, ideal)). **Contamination.** A contaminated run's slices cost at most
the slice walk: per slice <= D_SUP = `contam_sup` = 2^72.678 on [2^56, 2^64) (each band at the exact max of h x
REC_TAB(h) up to the largest region slice, empty cosets and setup6 in full; [k7c]: the slice walk at 2001
region slices peaks at 2^72.542). Every row pays RUN_ALLOW x D_SUP per slice, RUN_ALLOW = 0.00075: pool draws are
uniform over canonical a; in 4,000 a-uniform runs (1,000 of recP v12gp, 3,000 of **PREREG-C**: seeds 5001-5010, the
first 300 canonical a with Q_a <= theta A^2 per seed, sha256-16 in `pool_d.json`) the exact predicate found no
non-canonical a-value among 1,607,237: contaminated runs hold under 1 - 0.05^(1/4000) =
0.00074873 <= 0.00075 of the canonical a (95%, a modeling margin). Work sums (`clean_run_adjust`): per row RUN_ALLOW x
D_SUP x (1 - qa)^3 / (17/20 - qa)^3 over all S and G region rows (70,768.9 / 70,712.3), work -= 0.00225 D_SUP x sum,
anchors likewise, then the exact deltas above; work2 recomputed. [k7c] replays 6 sampled runs, [2r] 3 anchor runs.
Count envelopes ([3]) leave the per-a charge >= 0.985464/1.013857^2 x RUN_ALLOW x D_SUP = 0.000719 D_SUP (need about
0.000702 D_SUP). The rule's own work (runs16 ends, ~2^20 instructions per canonical a, < 2^69.3 units) is in misc.

**Yield ([4], pool G).** rels = (va_region/589502.674) vb (0.971187233) sum over G's canonical primitive region rows of
w rho(u_f) rho(u_g) r_f r_g 0.999^(lambda_f + lambda_g) x the band margin: Dickman rho ([1]), retentions r from the
Euler-product threshold model (`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14; prime_min
2^14), the ECM miss; cells with qa < 0.02 are sieved and charged, their relations discarded. Estimate
2^63.63225 (SE 0.0234 bits) less the band charge 0.009: **2^63.6232** against 1.05 U <= 2^63.52173: margin >= 0.1015
bits = 4.3 pool SE (feasibility requires max(0.06, 2 SE)). **Band charge (H1').** On 4,000 random rows of this region
(exact norms, all p < 2^14) the content of the unsieved [2^8, 2^14) is f 5.983 +- 0.130 bits (model 5.936; P(> 40
bits) 0.0030, 95% upper 0.0049, model 0.0021), g 6.358 +- 0.130 (+0.42 bits, 3.2 SE). The g excess at its 2-SE end as
a 0.75-bit g-slack loss (-0.0026 bits) and the f tail excess at its 95% limit counted lost (-0.0040): 0.0066 <=
REL_BAND_LOSS = 0.009 bits, deducted in [4]. **Descent input.** `LOG2_RELS = 63.575` is at or below the charged count
lowered by two relative pool SE, 2^63.57570 ([4] `rels_input`); the mean norms 739.1/498.2 truncate G's
739.167/498.217 (lower means: covered by [10y]; 2^16's ratio + 2 SE 1.0077, H3).

**Cost ([5]).** Sampled mean 155.3 instructions/cell over pools S and G, charged 156.5 with the allowance, over
2^125.600 cells (bounds, layout, setups, line path, record walk or rec7 with contamination, routing, apply, recovery,
scan2, band strip and cofactorisation, at the constants of KERNELS.md and STORAGE.md): collection **<= 2^118.602519919
units** (estimate 2^118.5911 times allowance 1.00794). **Allowance (H2, [4a]).** The work is a Monte Carlo estimate of
a fixed, challenge-independent total over the 1,312,000 draws of S and G; its relative SE has a row term 0.1529% and an
incidence term sqrt((1-70503/1312000)/70503) = 0.3664%, total <= 0.3970%; the charge multiplies the estimate by
**1.00794 >= 1 + 2.00 x 0.003970 = 1.007940** (rounded up; [4a] fails below), a modeling margin, not a probability.

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U, `filtsim.PEXC`), clique removal and Markowitz merge (coefficients <= 255) are simulated on
the production occurrence profile (`filtsim.py`; [6] reruns 2^12) as the signed sigma-quotient system. **Cap.** All
retained relations are used; H1 puts the count at >= 2^63.5757 = 1.0895 U >= 1.08 U ([4]): filtering at excess 1.08.
ECM misses make N random: H3 bounds E[N(300)/U] by 0.03329; the attack **aborts** if N > 0.03330 U; LA is charged at 0.03330. [10c]: per-run sd x sqrt(U) <= 0.643 at every scale
with >= 6 seeds (2^12..2^20, both excesses). H3c: sd <= 1/sqrt(U) at U = 2^63.45; Chebyshev: P(N > 0.03330 U) <=
(2^-31.73/1e-5)^2 = 7.9e-10; the count falling below 1.08 U adds < 1e-15. P_CAP = 1e-9 is charged (section 7). 1.08 was set after G's count;
sizing and G's judgment stay at 1.05 (G's margin over 1.08 U: 0.0609 >= max(0.06, 2 SE), [4]). A campaign at P10c's judged means 740.8/499.4 (`FILT_RUNS`, 6-16
seeds per scale 2^12..2^20; 2^18-2^20 pre-registered, PREREG-G/H; at excess 1.08 PREREG-J, 16 seeds at 2^19, 8 at
2^20) measures x = N(300)/U ([10] prints every scale): 0.03633 (SE 0.00019), 0.03540 (0.00018), 0.03394 (0.00022); at
excess 1.08 0.03301 (0.00014), 0.03161 (0.00022). H3: production is not above the largest scales; the pre-registered
charge is the larger mean + 2 SE of the two largest: max(0.033286 at
2^19, 0.032038 at 2^20) -> **0.03329** ([10]). At G's means 739.1/498.2 the same seeds (excess 1.08) give a pooled
paired ratio 0.9650 +- 0.0084 (2^16, 2^17; `FILT_Y`, [10y]; 0.9818 <= 1 at 2 SE). Merge to mean weight 300 (charged
301), largest merged coefficient 253 <= 255. N <= 0.03330 U <= **2^58.5430** < 2^58.584 < 6 x 2^56: six owners.

The square system A = [M | S] has N + 21 rows, N ideal columns and 21 dense Schirokauer columns (sigma-odd unit ranks
18 + 3 from signatures K_f (0, 36), K_f^+ (0, 18), K_g (2, 5), K_h (2, 2), [13]); under H3 its kernel is the line of
virtual logs. Block Wiedemann (m, n) = (72, 5): Y uniform in F_r; X block-diagonal (H-X): the N + 21 coordinates
shuffled uniformly (filtering-streams row), cut into m contiguous blocks, X_l uniform in [0, 2^128) on block l, else 0. Krylov N/m + N/n block steps; Horner runs as five interleaved chains from
checkpoints A^(jL) Y (L = Horner length / 5; 20 checkpoint vectors, < 2^68 bytes, 4L < N/m + N/n) at the
five-vector rate. Per row and vector (STORAGE.md, [10]): 301 x 9.97 (H4u, PREREG-U; |c| in {2, 3, 4} as |c| unit
entries) + 1,215 (row tail 1,210, six owners) + 21 x 97; per Horner element n dense terms and one reduction (420); projections: each vector element
meets one X entry, (1 + n/m) N^2 terms at 60. Linear algebra **<= 2^115.814110851 units** (m = 72, n = 5, kept from the
previous k5; not re-optimised); generator (matrix Berlekamp-Massey, recursive doubling, Karatsuba polynomial-matrix products, k = m + n = 77;
`model16.generator`: 4 k^3 (97 + 420) L^(log2 3) log2 L instructions, L = (N + 21)(1/m + 1/n) + 2): **<= 2^110.604 units**.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent: the relation count is fixed up to ECM misses (relative
sd below 2^-31): >= 2^63.6232 under H1 (section 3). Under H3 (r unramified in K_f and K_g, checked by [13]; r divides neither class number; the 21 sigma-odd Schirokauer maps are
injective on units modulo r-th powers; kernel of A = the line of virtual logarithms) Block Wiedemann fails with
probability <= 4N/2^128 < 2^-66 (Schwartz-Zippel, degree <= 2N, entries uniform on >= 2^128 values; H-X); a failed solve
(checked on the dropped rows) is a failure.

**Reconstruction.** Purged ideals' logs come by back-substitution; the descent treats every ideal with <= 2
occurrences as log-less ([6b]: the model's count >= 6.3x the measured one; DESCENT.md).

## 6. Individual logarithms

DESCENT.md specifies the descent (ECM-trial boot, special-q nodes sieving exact lattice balls, retention >= 0.8592
([7]; loss bounded with a dominating law over a lower-bounded denominator), negative-binomial counts, lazy nodes and boot ([7l]: 0.0224/0.0161/0.0368 simulated vs bounds
0.0243/0.0173/0.0420)). Plan: 62 nodes + boot (215, 328, 16,
4); at the 2-SE input min mu 44.7, max phi 0.0268, max psi 0.536.

Expected two-log work **<= 2^112.4260 units** ([7]; unrounded 2^112.42596; one log 2^125.714 instructions: boot own
2^123.70 + subtree 2^124.20 x retry 2.153). The meter caps both logs at **2^113.279 units** (reached = failure):
P(cap) <= E[W]/cap <= 2^(112.4260 - 113.279) <= 0.5536323 (Markov); with 10% more expected work, <= 0.6089956.

## 7. Success

Core: success >= (1 - 2^-66 - P_CAP) x (1 - 2 x 1.307e-05 - 0.5536323 - 2^-11.61 - 2^-10.70 - 2/r) >= 0.44542 (main
phase, two boots with all 18 relations failing, Markov, node counter, repeats); >= 0.39005 with Markov at +10% (0.6089956).

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX = 2^27.3. E[N] uses the same
  recurrence with unit node cost (`counts16.py`): E[N] <= 2^15.6813 for both logs (unrounded 2^15.68129), P(abort) <=
  2^-11.61.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for the 1/p weighting: P <= 2^(2 x 27.3 - 1
  + 6 - 70.30) <= 2^-10.70, charged.

Subtrees are independent (H5). With the wrapper: P(success) >= 0.44542 - 2^-48 > 0.39; every abort is a failure.

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (section 3, H-EA; estimate 118.5911 + 0.794% allowance) | <= 118.602519919 |
| linear algebra (weight 301, N <= 2^58.5430, (72, 5), block-diagonal X; section 4) | <= 115.814110851 |
| matrix generator | <= 110.604 |
| both individual logarithms (metered cap; expected 112.4260) | 113.279 |
| canonical a list | 91 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 72 |
| wide factor-base lattices ([2^63, 2^64), KERNELS.md) | 70.6 |
| rec7 bases (2 per ideal, p in [2^52, 2^64), KERNELS.md) | 73.5 |
| band split (section 3, [10t]) | 101.53 |
| numerical certificate construction | 50 |
| verification | 11 |
| wrapper: pairing, field map, scalar multiplication, check (section 1; <= 2^80 instructions) | 66 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^91 (2^41.8 per eligible a, 2^63.451 of them: 2^90.963,
rounded up; one HNF key, the orbit enumeration and a sort); Schirokauer
maps 2^88 (2^38.5 per relation against < 2^34.5); factor base 2^72 (2^23.4 per prime against < 2^19); filtering and
back-substitution 2^80 (2^24.4 per entry against < 2^14). The upper-rounded rows sum to at most 2^118.83350433992169, so
the core meter is **2^118.8335043400** (audit-checked, as is the stated sum).

**Total.** With the wrapper row,
W <= 2^118.8335043400 + 2^66 < 2^118.833504340001, evaluated with 80-digit decimals by `audit_refinement.py`.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) above section 2's values; exact sampled norms). H1' retention (Euler-product model calibrated on exact
smooth integers; the band-content excess of [2^8, 2^14) charged, 0.009 bits). H-EA (work only): H1''s survivor law prices the ECM share; the rule on random integers
at 1/4 scale ([12w]: sieve [2^4, 2^16), windows (2^16k, 2^17.25k]) near 2^120 / 2^136 / 2^150 (5,156 / 2,213 /
433 survivors): ECM path 0.1528 / 0.1685 / 0.1640 vs 0.2211 / 0.2169 / 0.2133 charged. H-BT (work only): the Toom-3
recurrence bound holds past the measured 34,992 limbs (the code's own recursion). H2 region (uniform pools, exact flags, certified counts; sized on S, judged once on G; S+G work with the 0.79% allowance). H3 matrix (r unramified, [13]; r not dividing the class numbers, the unit condition and
a one-dimensional kernel assumed; production E[N(300)/U] <= 0.03329 and H3c, section 4). H4 kernel cost (KERNELS.md, STORAGE.md: static counts of loop-free code, executed counts for routing); H4u (>= 86%
+-1; k5 <= 9.97, PREREG-U: larger mean + 2 SE at 2^16, 2^17: 9.9654, 9.8946 (means 9.942, 9.875); k5 falls with scale); the wide path's
per-hit price is the measured 0bbca0b2 one (same walk code); rec7's per-hit price depends on h' only, for every p <
2^64, on runs of 32 to 1,024 a-values (longer or shorter runs take the slice walk); contaminated runs hold at most 0.075%
of the canonical a (95% bound; a shortfall changes cost only); f/g independence given q_f (charged at 100x). H-X: BW's nondegeneracy polynomial is nonzero on shuffled block-diagonal X (unit vectors are a special case); [10x]: merged toys (21 dense columns, a kernel vector) reach full block Hankel rank at
(72, 5) in 28 of 28 (2^12 seeds 1-16, 2^13 1-8, 2^14 1-4). H5 descent (DESCENT.md): lower sums; lazy order independent of outcomes; negative-binomial recurrence, independent
subtrees; missing logs (Poisson occurrence at the 2-SE input; re-descended in L_(Q q') or discarded); local omitted-mass laws;
random-integer ECM orders; UDUP = 0.65 (complete unit orbits against the charged pr, DESCENT.md; live [11] 0.896 / 2.05).

## 10. Memory

Core storage <= 2^105 bytes (STORAGE.md), wrapper below 2^30: below 2^106 bytes. All preprocessing is for this target,
charged once in W.

## Provenance

GT core: kaankolcu recover-P 1b8fa720 and its lineage (9332a74d mpjunior92; 2c2ce2ab / 42e80f80 gbotrel; 6282ab79
bernardo-lima12; 6734bd5a GPT-6 in Codex; Gemini 45ce22e5 / 4bb6c381); polynomials after Guillevic-Singh; filtering after CADO-NFS. H-EA, H-X,
rec7 to 2^64, band strip and split: gbotrel, Claude Opus 5.5 in Claude Code (early abort after Pomerance, Barbulescu). Wrapper: 4f8baf8c.

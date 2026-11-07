# recover-Q (bls12-381) through two GT logarithms with the plain-track STNFS core (plain P11): 2^117.980413847301 reference units

Track `bls12-381-fixed-g1-invert-g2`. Input: a nonidentity P in G1 and z in GT; output Q in G2 with e(P,Q) = z
for the exact supplied P (Q = infinity for z = 1). The score unit is one dense F_p12 multiplication, C = 20,000
scalar AArch64 instruction equivalents, the unit the core below counts. No advice, smaller targets or amortization
are used. Predecessor work is reused as credited at the end; no prior verdict is a premise.

## Claim

The algorithm below returns such a Q with success at least 0.39 over uniform
challenges and fresh coins under the declared heuristics H1-H5, H-EA, H-X and H-FFT (section 9). All classical work (preprocessing, both logarithms, failures,
restarts, field map, pairing, scalar multiplication, verification, output) is at most

    W <= 2^117.9804138473 + 2^66 < 2^117.980413847301  dense-Fp12 units (C = 20,000).

The score is 117.980413847301. 2^117.9804138473 is the GT core's global meter (sections 2-8): upper-rounded ledger
rows sum to at most 2^117.98041384721757; 2^66 is the wrapper (section 1). Storage is below 2^106 bytes. Analytic,
exploratory claim.

**Core.** Collection **2^117.531499000** (f-only sieve, exact g batch, band strip and split, section 3), linear algebra
**2^115.814110851**, descent meter **2^113.285**.

**Conventions.** Bounds round in the safe direction; `evidence16.py` checks each stated bound against its unrounded
value (`declared_bounds`) and reruns the checks cited [n]; `audit_refinement.py` is the hosted audit. Model:
`config.json`, `stats16.b85`, `pool_d.json`, the constants of `model16.py` and `descent16.py` (`owner10.install`: C_ROW
1,215).

## 1. Reduction (fixed P: two GT logarithms) and field map

With the pinned reduced pairing e and the public generator Q0 of G2:

1. Check P (on curve, [r]P = infinity, P != infinity) and z (z^r = 1). If z = 1 return Q = infinity.
2. g = e(P, Q0) (pinned `pairing(Q0, P)`, protected encoding).
3. Map g and z into the core's tower (field-map.md).
4. Run the core on (g, z) under its meter and node counter (any failure is a failure): a = log_g(z) mod r.
5. Q = [a]Q0; output Q if e(P, Q) = z, otherwise fail.

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
[0] checks Res_y(h, f(x^2 -> u y)) = 0 mod p. x -> -x is the Frobenius^6 conjugation sigma (log sigma(P) = -log P),
quotiented out. Murphy alpha is charged at (5.15, 0.85) in collection and descent;
[14] recomputes it (exact Kummer-Dedekind structure): (4.868096, 0.843643) over p <= 2000, (5.085035, 0.731543) over
p < 50000, both below.

## 3. Relation collection

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69; the factor base is the relative-degree-1
prime ideals over K_h of norm <= B on both sides, U = 2^63.451 sigma-orbits.

**Region.** Exact rational forms Q_a, Q_b (`config.json`), A^2 = 2^21. Pool: uniform on theta = 3/8 (a5 > 0, Q_a(a) <=
(3/8) A^2; b uniform in the half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a)). Region: a5 > 0, Q_a(a) <= (2719/8000) A^2; b5
>= 1, Q_b(b) <= (17/20) A^2 - Q_a(a) (invariant under b -> -b and the sign fold; exact Fractions), so the restricted sample is uniform on it. G has 29835 region rows (27256
primitive) of 78960 canonical draws, scaled by the region volume over the effective sample 656000 x (2719/3000)^3 =
488391.051. Two weights, never mixed: density w = (1 - qa)^3 (yield) and vb (1 - qa)^3 x factor (work, population);
slice-geometry costs use S' = vb (17/20 - qa)^3 x factor and K' = floor(sqrt((17/20 - qa) A^2 / Q_b(e_5))) lines. The
canonical a of a unit class is the least of its exact HNF key class in (Q_a, lex) order (one point per relation orbit
of <units, sigma>); non-primitive pairs are excluded by an exact HNF test. `lattice_count_bounds.py` certifies (exact
rational Machin-pi, determinant/sqrt intervals) the count ratio at every radius R >= 844 = isqrt((2719/8000) 2^21),
below the smallest b-slice radius 1034 ([3]): [0.985050, 1.014303] (Q_a), [0.985002, 1.014303] (Q_b); the model uses
the products **0.970276956 (yield)** and **1.028809574 (work)**. 

**Pools (PREREG-G).** 656,000 uniform eligible pairs each on theta = 3/8 (indices 0..163999 per seed): S
(seeds 261-264) sizes the region, G (251-254, drawn after the configuration was fixed) is judged (rerun for the g batch), both give
the work; every draw was tested with the exact canonical predicate, none dropped. `pool_d.json` holds per-seed counts, list
hashes and sums, 8 anchor rows of G, G's canonical lists over draws 0..49 ([2] replays them) and the contamination sample. theta is the smallest t on the 1/8000 grid with margin_S
>= 0.06 + sqrt(2) SE_S (pre-registered), re-applied to the f-only law: 2719/8000 (0.100413 >= 0.100330; 2718: 0.099661 < 0.100351; [4s]). The band
charge applies to G's count only.

**Sieve and survivors (f ideals only).** Line sieve (scores only) for p in [2^14, 2^48); record walk for p in [2^48, 2^64) (`slice12`,
KERNELS.md; X2 = 2^64; [2^56, 2^64) per a-run, below; [2^63, 2^64) of contaminated runs on the wide path); primes in
[2^64, 2^69] by ECM within slack 302. Survivors pass the f threshold; exact norms (C_NORMS) give the g norm z for the g batch. Primes below 2^14 are the retention models' omitted mass (`prime_min`
2^14, `explicit_bits` 14). Survivors with z 2^69-smooth get **trial division by every prime below 2^20 on both cofactors** (82,025 primes;
5 per 64-bit limb at the 832/576-bit caps (2^(64 i) mod p table), 22 per side, 4 loop, one division per cofactor bit:
2^23.649 <= **C_TDIV = 2^24.5**, [5t]). **Batch band strip ([10b], KERNELS.md).** Survivor sides wait in
waves of K = 2^39 (13-limb cofactors, <= 832 bits). With P the product of the primes in (2^20, 2^48] (< 1.000028 x
2^48/ln 2 bits, Rosser-Schoenfeld; built once), a scaled remainder tree (Bernstein; `model16.srt`)
gives P mod z for every z, and the band part is gcd(z, (P mod z)^(2^6) mod z) (2^6 >= 41 >= the exponent of any band
prime in z); the f band-miss factor (1 - 41 x 0.001^2) stays as a margin.
**g batch ([10g]).** Same SRT with P_g the product of all primes <= 2^69 (< 1.000028 x 2^69/ln 2 bits), waves of K = 2^61 g norms (<= 576 bits, 9 limbs): z is 2^69-smooth iff (P_g mod z)^(2^10) = 0 mod z (exponents <= 576 < 2^10). r_g >= 1 - 2^-38 (ideals of norm p^d > 2^69 drop a relation at the split; in the margin). A kept relation's g norm is split by trial division, the band split and ECM (<= 12 primes of (2^48, 2^69], until complete): <= 12 x 1.002 x 2^39.79 + 13 C_MR = 2^43.38 per relation, <= 2^65.15 relations ([10t]): **2^95 units** (ledger row).
**Lemma (exact SRT).** Product tree Z_v (pairs in order, odd node carried; L <= 61 levels); node v holds Y_v of b_v =
64 (n_v + G) bits, n_v = limbs of Z_v, G = L + 2. Root: Newton x ~ 2^s/Z (s = b + |P| + 8), q = floor(P x/2^(s-b)),
check R = P 2^b - q Z, <= 2 add/subtracts of Z bring R into [0, Z)
; Y_root = q mod 2^b. Child c, sibling s: Y_c = floor((Y_v Z_s mod 2^(b_v))/2^(b_v - b_c)). Leaf: r = ceil(z Y
2^-b), 0 if r = z. Then r = P mod z. Proof: y_v = frac(P/Z_v), e_v = y_v - Y_v 2^(-b_v) mod 1, eta_v = e_v Z_v. Root:
0 <= eta < Z 2^-b <= 2^(-64G). Z_v = Z_c Z_s gives y_c = frac(y_v Z_s), so e_c = e_v Z_s + t, t in [0, 2^(-b_c)), and
eta_c <= eta_v + 2^(-64G): all eta < (L + 1) 2^(-64G) < 1 (no wrap). At a leaf z Y 2^-b = (P mod z) - eta mod z, so
the ceiling is P mod z (or z for 0). Any G >= 1 is exact. **Cost** (H-FFT, KERNELS.md): an N-limb product <= F(N) = 8.5 N log2(N)^2 + 1200 (+ 504 N transposes above 2^26 limbs);
per level products of 2 n_(j-1) and n_j + G + n_(j-1) limbs + 3 per copied limb; root <= 6 F(2(max(|P|, |Z|) + G)); leaf <= 2^15
(g) / 2^19 (band) (`model16.gbatch_cost`): **C_BATCH = 2^24** (2^23.921, K = 2^39), **C_GBATCH = 2^25** (2^24.973, K = 2^61). [10v]: `model16.srt` on 2^8
planted leaves, Newton offsets -2/+2, against Python's %.
P's build (< 2^70 instructions) is in misc, P_g's tree (< 2^88) is the 2^74 row; trees < 2^74 bytes.
**Band split ([10t]).** The strip decides smoothness only. A kept relation (both sides smooth, primitive, qa >= 0.02)
is split: the band parts are trial-divided by every prime of (2^20, 2^48] (< 1.25506 x 2^48/ln 2^48 = 2^43.272;
table of 2^(64 i) mod p, i < 13, < 2^50.2 bytes, built in the misc row; 6 per limb on 13 + 9 limbs, 22 per side, 4 of
loop control: 180 per prime), each hit divided out to full power (<= 41 + 28 divisions, 2^10 each). Every found prime
gets its ideal (p, s, r) and exponent from the factor-base entries above p (a(s) = r b(s) mod p; several matches: Hensel
lift to the found power), <= 2^26: **2^50.7634** instructions per relation, deterministic, no relation lost. Kept
relations are at most the B-smooth region cells: G's 2^63.63561 + 2 SE over the yield envelope
(1.028809574/0.970276956), the ECM miss (0.999^22: <= 13 + 9 primes above 2^64), the band margin and the least
retention on the model's norm grid (0.7775 f, 0.5393 g): <= 2^65.0557, **2^101.5314 <= 2^101.54 units** (ledger row;
hosted audit replays it).
Survivors' f omitted mass is bounded by 310 bits; ECM is the `ecm16.py` schedule at **99.9%** per
prime ([12]): 69-bit primes in a 310-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 403 curves, five-limb Montgomery
multiplication at 367.76 Ir) cost 2^38.0909 instructions.
**Cofactor windows (H-EA, `model16.Window`).** After trial division, records and band strip, a side's cofactor C (< 2^320) has no prime below 2^64 unless the side is not smooth (a prime of
[2^48, 2^64) through no factor-base ideal). So the side is smooth only if
C = 1, C <= 2^69 is prime, or C is composite with 64k < log2 C <= 69k, k in {2, 3, 4}; ECM runs only in the last case,
every other side is rejected. Compositeness: 40 Miller-Rabin rounds, random bases (C_MR = 40 x 2 x 310 five-limb
products = 2^23.121); a composite passes with probability <= 4^-40 (Rabin-Monier), below the
count's rounding. ECM is charged at e x its schedule, e = the survivor law's mass whose part above 2^64 lies in a
window widened by 1 bit, minus its single-prime term, over the survivor mass (e_f 0.1555-0.1570 at
the anchors). **g first:** the f side's division, strip, test and ECM run only if the g batch passes:
per f-survivor C_NORMS + C_GBATCH + q_g (C_TDIV + C_BATCH + C_MR + e_f 2^38.0909), q_g = min(1, 100 rho(u_g)) per pool row
(H1, factor 100 a margin); f-threshold survivors per cell (~2^-22.5). S and G work sums use this
model on all 1,312,000 regenerated draws (`pool_row`, `runs16`).

**Thresholds.** Certified per-block bounds (`embedding_bounds.py`, 64-cell blocks of width 128, bisected until the f bound is within 8 bits; the g bound is unused)
give T = lower - 16 x slack, keeping every B-smooth cell with omitted mass below slack: F_v Taylor-shifted to the block
centre per place, remainder over the block disk plus a 2^-32 relative term, rounded outward; one evaluation (`NormBounds.block`) is what `nbound31` replays bit for bit at 1,508.2 Ir. The rows of one line
in an aligned 2 x 2 x 4 box of (b2, b3, b4) (`model16.box_rows`, <= 16 rows) share one bisection per block over the
hull of their contexts (`EB.hull`, outward); a column the hull leaves unresolved uses each row's own context. Boxes x blocks x (b1, b5) partition each slice, so each cell pays C_NB_EVAL x its box's evaluations per box
cell; each row pays its leaves, context and hull share (C_NB_ROW 4,600).
[2b] checks 1,536 exact norms of six rebuilt boxes against their leaf bounds, [2] the anchors' box fields.
**Warm start (clean runs).** A box-block's bisection starts from the leaves that the run's first slice holding it
certified from the aligned 64-cell root (`model16.pool_row_warm`): same tests and fallback, each leaf certified for its
own a, so no relation changes. Charged: warm evaluations (+ RUN_ALLOW x the cold excess), max(cold, warm) leaves in
scan2 (evaluations -9.44%, leaves +10.96% over S and G); [2] and the hosted audit replay the anchors' counts.

**Record walk per canonical a-run (rec7; KERNELS.md).** A canonical a's run is the maximal integer t-interval
containing 0 with a + t e_0 in the region and no image of smaller Q_a under the 11,325 listed units (`runs16.py`). On a run, a prime ideal (p, s, r)'s hits in all its slices form one coset of a 7-dim
lattice (t = r b(s) - a(s) mod p), walked once per (run, ideal) with one setup6 ([10r]; bases, checks and campaign in
KERNELS.md).
**Clean-run rule:** rec7 walks a run only if every a-value of it is in the canonical list (whose exact test the 2^91
row pays); otherwise its slices take the slice walk. A clean run of 32 to 1,024 a-values uses the run walk in each band
p in [2^56, 2^64) with h' = run cells / p >= 2^4.99 (the table's bottom) instead of the slice walk and empty cosets; run
cells are exact (`runs16.lrs`: integer sum over the run's a-values of (17/20 - qa)^3). Each row pays its cells' share
(per hit at h', + 2 per record, setup6 per (run, ideal)). **Contamination.** A contaminated run's slices cost at most
the slice walk: per slice <= D_SUP = `contam_sup` = 2^72.679 on [2^56, 2^64) (each band at the exact max of h x
REC_TAB(h) up to the largest region slice, empty cosets and setup6 in full; [k7c]: the slice walk at 2001
region slices peaks at 2^72.543). Every row pays RUN_ALLOW x D_SUP/2 per slice (f ideals), RUN_ALLOW = 0.00091: pool draws are
uniform over canonical a; in 4,000 a-uniform runs (1,000 of recP v12gp, 3,000 of **PREREG-C**: seeds 5001-5010, the
first 300 canonical a with Q_a <= (2895/8000) A^2 per seed, sha256-16 in `pool_d.json`) the exact predicate found no
non-canonical a-value among 1,607,237: contaminated runs hold under 1 - 0.05^(1/4000) = 0.00074873 of those canonical a
(95%); a run here is a sub-interval of the 2895-region run of its a, a (2719/2895)^3 share, so < 0.000904 <= 0.00091. Work sums: per row RUN_ALLOW x
D_SUP x (1 - qa)^3 / (17/20 - qa)^3, recomputed on every S and G region row; anchors likewise. [k7c] replays 6 sampled runs, [2r] 3 anchor runs.
Count envelopes ([3]) leave the per-a charge >= 0.985002/1.014303^2 x RUN_ALLOW x D_SUP = 0.000871 D_SUP (need about
0.000847 D_SUP). The rule's own work (runs16 ends, ~2^20 instructions per canonical a, < 2^69.3 units) is in misc.

**Yield ([4], pool G).** rels = (va_region/488391.051) vb (0.970276956) sum over G's canonical primitive region rows of
w rho(u_f) rho(u_g) r_f 0.999^lambda_f x the f band margin (g exact): Dickman rho ([1]), retentions r from the
Euler-product threshold model (`threshold_model.py`, calibrated on exact smooth integers at B = 2^10..2^14; prime_min
2^14), the ECM miss; cells with qa < 0.02 are sieved and charged, their relations discarded. Estimate
2^63.63561 (SE 0.0248 bits) less the band charge 0.009: **2^63.6266** against 1.05 U <= 2^63.52173: margin >= 0.1048
bits = 4.2 pool SE (feasibility requires max(0.06, 2 SE)). **Band charge (H1').** On 4,000 random region rows
(exact norms) the content of the unsieved [2^8, 2^14) is f 5.983 +- 0.130 bits (model 5.936; P(> 40
bits) 0.0030, 95% upper 0.0049, model 0.0021); the f tail excess at its 95% limit counted lost (-0.0040) <=
REL_BAND_LOSS = 0.009 bits, deducted in [4] (the g side is exact). **Descent input.** `LOG2_RELS = 63.575` is at or below the charged count
lowered by two relative pool SE, 2^63.57612 ([4] `rels_input`); the mean norms 738.1/498.2 truncate G's
738.11/498.25 (lower means: covered by [10y]; 2^16's ratio + 2 SE 1.0077, H3).

**Cost ([5]).** Sampled mean 82.4 instructions/cell over pools S and G, charged 83.1 with the allowance, over
2^125.442 cells (bounds, layout, setups, line path, record walk or rec7 with contamination, routing, apply, recovery,
scan2, band strip and cofactorisation, at the constants of KERNELS.md and STORAGE.md): collection **<= 2^117.531499000
units** (estimate 2^117.519229190 times allowance 1.00838). **Allowance (H2, [4a]).** The work is a Monte Carlo estimate of
a fixed, challenge-independent total over the 1,312,000 draws of S and G; its relative SE has a row term 0.1226% and an
incidence term sqrt((1-59553/1312000)/59553) = 0.4004%, total <= 0.4188%; the charge multiplies the estimate by
**1.00838 >= 1 + 2.00 x 0.004188 = 1.008376** (rounded up; [4a] fails below), a modeling margin, not a probability.

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U, `filtsim.PEXC`), clique removal and Markowitz merge (coefficients <= 255) are simulated on
the production occurrence profile (`filtsim.py`; [6] reruns 2^12) as the signed sigma-quotient system. **Cap.** All
retained relations are used; H1 puts the count at >= 2^63.5761 = 1.0903 U >= 1.08 U ([4]): filtering at excess 1.08.
ECM misses make N random: H3 bounds E[N(300)/U] by 0.03329; the attack **aborts** if N > 0.03330 U; LA is charged at 0.03330. [10c]: per-run sd x sqrt(U) <= 0.643 at every scale
with >= 6 seeds (2^12..2^20, both excesses). H3c: sd <= 1/sqrt(U) at U = 2^63.45; Chebyshev: P(N > 0.03330 U) <=
(2^-31.73/1e-5)^2 = 7.9e-10; the count falling below 1.08 U adds < 1e-15. P_CAP = 1e-9 is charged (section 7). 1.08 was set after G's count;
sizing and G's judgment stay at 1.05 (G's margin over 1.08 U: 0.0642 >= max(0.06, 2 SE), [4]). A campaign at means 740.8/499.4 (`FILT_RUNS`, 6-16 seeds per scale
2^12..2^20; 2^18-2^20 pre-registered; at excess 1.08, 16 seeds at 2^19, 8 at 2^20) measures x = N(300)/U ([10]): 0.03633 (SE 0.00019), 0.03540 (0.00018), 0.03394 (0.00022); at
excess 1.08 0.03301 (0.00014), 0.03161 (0.00022). H3: production is not above the largest scales; the pre-registered
charge is the larger mean + 2 SE of the two largest: max(0.033286 at
2^19, 0.032038 at 2^20) -> **0.03329** ([10]). At means 739.1/498.2 the same seeds (excess 1.08) give a pooled
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
sd below 2^-31): >= 2^63.7221 under H1 (section 3). Under H3 (section 9; the 21 sigma-odd Schirokauer maps
injective on units modulo r-th powers) Block Wiedemann fails with probability <= 4N/2^128 < 2^-66 (Schwartz-Zippel,
degree <= 2N, entries uniform on >= 2^128 values; H-X); a failed solve (checked on the dropped rows) is a failure.

**Reconstruction.** Purged ideals' logs come by back-substitution; the descent treats every ideal with <= 2
occurrences as log-less ([6b]: the model's count >= 6.3x the measured one; DESCENT.md).

## 6. Individual logarithms

DESCENT.md: ECM-trial boot, special-q nodes sieving exact lattice balls, retention >= 0.8592 ([7]; dominating law
over a lower-bounded denominator), negative-binomial counts, lazy nodes and boot ([7l]: 0.0224/0.0161/0.0368 simulated
vs bounds 0.0243/0.0173/0.0420). Plan: 62 nodes + boot (215, 328, 16, 4); at the 2-SE input min mu 44.7, max phi
0.0268, max psi 0.536.

Expected two-log work **<= 2^112.4316 units** ([7]; unrounded 2^112.43156; one log 2^125.719 instructions: boot own
2^123.70 + subtree 2^124.20 x retry 2.159). The meter caps both logs at **2^113.285 units** (reached = failure):
P(cap) <= E[W]/cap <= 2^(112.4316 - 113.285) <= 0.5534789 (Markov); with 10% more expected work, <= 0.6088267.

## 7. Success

Core: success >= (1 - 2^-66 - P_CAP) x (1 - 2 x 1.371e-05 - 0.5534789 - 2^-11.61 - 2^-10.70 - 2/r) >= 0.44557 (main
phase, two boots with all 18 relations failing, Markov, node counter, repeats); >= 0.39022 with Markov at +10% (0.6088267).

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX = 2^27.3. E[N] uses the same
  recurrence with unit node cost (`counts16.py`): E[N] <= 2^15.6877 for both logs (unrounded 2^15.68761), P(abort) <=
  2^-11.61.
- Repeats: all N_MAX nodes in the smallest bin (2^70.30 ideals) with 2^6 for the 1/p weighting: P <= 2^(2 x 27.3 - 1
  + 6 - 70.30) <= 2^-10.70, charged.

Subtrees are independent (H5). With the wrapper: P(success) >= 0.44557 - 2^-48 > 0.39; every abort is a failure.

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (section 3, H-EA; estimate 117.519229190 + 0.838% allowance) | <= 117.531499000 |
| linear algebra (weight 301, N <= 2^58.5430, (72, 5), block-diagonal X; section 4) | <= 115.814110851 |
| matrix generator | <= 110.604 |
| both individual logarithms (metered cap; expected 112.4316) | 113.285 |
| canonical a list | 91 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 73 |
| wide factor-base lattices ([2^63, 2^64), KERNELS.md) | 70.6 |
| rec7 bases (2 per ideal, p in [2^52, 2^64), KERNELS.md) | 73.5 |
| band split (section 3, [10t]) | 101.54 |
| P_g product tree (section 3) | 74 |
| g split of kept relations (section 3) | 95 |
| numerical certificate construction | 50 |
| verification | 11 |
| wrapper: pairing, field map, scalar multiplication, check (section 1; <= 2^80 instructions) | 66 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^91 (2^41.8 per eligible a, 2^63.451 of them: 2^90.963
rounded up; HNF key, orbit enumeration, sort); Schirokauer maps 2^88 (< 2^34.5 per relation, <= 2^65.15 relations); factor base 2^73 (2^23.4 per prime, pi(2^69) < 2^63.46); filtering and
back-substitution 2^80 (2^24.4 per entry against < 2^14). The rows sum to at most 2^117.98041384721757: core meter **2^117.9804138473** (audit-checked).

**Total.** With the wrapper row,
W <= 2^117.9804138473 + 2^66 < 2^117.980413847301, evaluated with 80-digit decimals by `audit_refinement.py`.

## 9. Heuristics

H1 smoothness (Dickman/Murphy; alpha (5.15, 0.85) above section 2's values; exact sampled norms). H1' retention (Euler-product model; [2^8, 2^14) band excess charged, 0.009 bits). H-EA (work only): H1''s survivor law prices the ECM share; the rule on random integers
at 1/4 scale ([12w]: sieve [2^4, 2^16), windows (2^16k, 2^17.25k]) near 2^120 / 2^136 / 2^150 (5,156 / 2,213 /
433 survivors): ECM path 0.1528 / 0.1685 / 0.1640 vs 0.2211 / 0.2169 / 0.2133 charged. H-FFT (work only): F(N) bounds the local SSA work past the measured 2^26 limbs (KERNELS.md, 2.76x headroom). H2 region (uniform pools, exact flags, certified counts; sized on S, judged on G; S+G work with the 0.838% allowance). H3 matrix (r unramified, [13]; r prime to the class numbers, the unit condition, a one-dimensional kernel; E[N(300)/U]
<= 0.03329 and H3c, section 4). H4 kernel cost (KERNELS.md, STORAGE.md: static counts of loop-free code, executed counts for routing); H4u (>= 86%
+-1; k5 <= 9.97, PREREG-U: larger mean + 2 SE at 2^16, 2^17: 9.9654, 9.8946 (means 9.942, 9.875); k5 falls with scale); the wide path's
per-hit price is the measured 0bbca0b2 one (same walk code); rec7's per-hit price depends on h' only, for every p <
2^64, on runs of 32 to 1,024 a-values; contaminated runs hold at most 0.091%
of the canonical a (95% bound; a shortfall changes cost only); f/g independence given q_g (charged at 100x). H-X: BW's nondegeneracy polynomial is nonzero on shuffled block-diagonal X; [10x]: merged toys (21 dense columns, a kernel vector) reach full block Hankel rank at
(72, 5) in 28 of 28. H5 descent (DESCENT.md): lower sums; lazy order independent of outcomes; negative-binomial recurrence, independent
subtrees; missing logs (Poisson occurrence at the 2-SE input; re-descended in L_(Q q') or discarded); local omitted-mass laws;
random-integer ECM orders; UDUP = 0.65 (complete unit orbits against the charged pr, DESCENT.md; live [11] 0.896 / 2.05).

## 10. Memory

Core storage <= 2^105 bytes (STORAGE.md), wrapper below 2^30: below 2^106 bytes. All preprocessing is for this target,
charged once in W.

## Provenance

GT core: kaankolcu recover-P 1b8fa720 and its lineage (9332a74d mpjunior92; 2c2ce2ab / 42e80f80 gbotrel; 6282ab79
bernardo-lima12; 6734bd5a GPT-6 in Codex; Gemini 45ce22e5 / 4bb6c381); polynomials after Guillevic-Singh; filtering after CADO-NFS. H-EA, H-X,
rec7 to 2^64, band strip and split, f-only sieve with the g batch: gbotrel, Claude Opus 5.5 in Claude Code (early abort after Pomerance, Barbulescu). Wrapper: 4f8baf8c. Base dd522245 (gbotrel); cap 113.285: kaankolcu.

S/G yield and work sums recomputed on all 59,553 region rows for the f-only sieve (`pool_d.json`).

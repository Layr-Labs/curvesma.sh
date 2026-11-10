# Plain G1 ECDLP (bls12-381) through two GT logarithms with the STNFS core behind a run gate: expected work 2^117.305516092225 units

Track bls12-381-g1-ecdlp-plain. Input: T = [x]P0 in G1; output the scalar x in 0..r-1.
Score unit: one batched affine G1 addition, C = 2,000 scalar AArch64 instruction equivalents. The core meters
core units of 20,000 instructions = 10 score units. No auxiliary powers, advice, smaller targets or amortization
are used.

## Claim

W is the total classical work of the one attack on the one target. Step 0 is a run gate (kaankolcu's device, Provenance): draw a 64-bit word u;
unless u < K = 7203760764237811756, stop as a failure (< 1 unit). K is the least integer with K x 0.998677 >= 0.39 x 2^64. Under H1-H5, H-EA, H-X and
H-FFT (section 9) the algorithm below returns x with probability >= 0.998677 once past the gate (section 7); u is independent of the challenge and the
other coins, so success is at least 0.39 over uniform challenges and fresh coins, and

    E[W] <= 1 + (K/2^64) x 10 x (2^115.3401320190 + 2^66) < 2^117.305516092225  units (C = 2,000).

2^115.3401320190 core units: section 8's rows at their every-run bounds, the two logarithms at their expected work 2^109.902 (section 6). The score 117.305516092225 is an
expectation, not a bound on every run: all work of a run (preprocessing, both logarithms, failures, restarts, field map, both pairings, scalar multiplication,
verification, output) is at most

    W <= 10 x (2^121.9148622354 + 2^66) < 2^125.236790330288  units.

Nothing is amortized over other targets or runs. Storage is below 2^106 bytes. Analytic, exploratory claim.

Conventions. Bounds round in the safe direction; evidence16.py checks each against its unrounded value
(declared_bounds) and reruns the checks cited [n]. Model:
stats16.b85 (config, pool_d), the constants of model16.py and descent16.py (owner10.install: C_ROW
500).

## 1. Reduction (fixed base g = e(P0, Q0): two GT logarithms) and field map

With the pinned reduced pairing e and the public generators P0 of G1 and Q0 of G2:

0. Run gate (Claim): draw u; if u >= K, fail.
1. Check T (on curve, [r]T = infinity). If T = infinity return x = 0.
2. g = e(P0, Q0), z = e(T, Q0) (pinned pairing(Q0, .), protected encoding).
3. Map g and z into the core's tower (field-map.md).
4. Run the core on (g, z) under its meter and node counter (any failure is a failure): a = log_g(z) mod r.
5. Output x = a if [a]P0 = T, otherwise fail.

GT has prime order r and e is nondegenerate: g generates GT, e(., Q0): G1 -> GT is an isomorphism, z = g^x, z = 1
iff T = infinity, a = x; uniform x gives uniform z. Each core target t (g or z) becomes t^e, e fresh uniform
in [1, r) (DESCENT.md, Targets), uniform on GT minus {1} for any fixed g, z: the fixed base costs no success. The field isomorphism preserves z = g^a; the final check rejects any wrong output.

One relation/LA phase serves log z and log g (section 6; 2/r charged). Wrapper below 2^80 instructions (two pairings, scalar multiplication, check < 2^50, field map < 2^67, field-map.md):
2^65.71 core units, a 2^66 row. Wrapper failure: field map below 2^-48.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c], [14])

K_h = Q[y]/h, h = y^6 + 2y^4 - 2y^3 - 2y^2 - 1, irreducible mod p. Disc(h) = 2^13 x 19^2; O_h has no power basis:
common12.py fixes a Z-basis e_0 = 1, ..., e_5 of O_h (index 4 over Z[y]); [0c] certifies maximality (19: Dedekind;
2: runs16.py, 2 O_h = (w)^2 (pi4)^2 with (w), (pi4) maximal of norms 2, 4): Disc(O_h) = 2^9 x 19^2. With w = y + 1 (N(w) = 2), f = w^6 F(x^2/w)
= x^12 - 2w x^10 + 2w^3 x^6 + w^5 x^2 + w^6 (F the BLS12 polynomial) and g = x^2 + |u| w, u = -0xd201000000010000. [0b] checks
(u w)^((p^6-1)/2) = -1 in F_p6, so F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u w)) = 0 mod p.
Twisted symmetry. sigma(y) = -y^5 - 2y^3 + 2y^2 + 2y is an automorphism of K_h of order 2, c = (y^5 + 3y^3 - 3y^2 - 3)/2
a unit with c sigma(c) = -1; tau: (y, x) -> (sigma y, c x) maps f to c^12 f and g to c^2 g ([0b]), tau^2 = (x -> -x). On F_p12 tau is Frobenius^3 ([0b]: sigma(y) = y^(p^3), (u w)^((p^3-1)/2) = c), so
log tau(P) = chi log P with chi = p^3 mod r, chi^2 = -1. One relation stands for its tau-orbit (4 pairs) and the factor
base is quotiented by tau-orbits of 4 ideals. Murphy alpha (charged at (0.29, 0.75));
[14] recomputes it (Kummer-Dedekind for odd p, exact local laws at 2 from runs16.py: f -1, g -1.39998):
(-0.047244, 0.706684) over p <= 2000, (0.285135, 0.745775) over p < 50000, both below.

## 3. Relation collection

Candidates: (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69; factor base: relative-degree-1 ideals of
norm <= B, U' = 2^62.451 tau-orbits (U/2: orbits of 4 ideals; U = 2^63.451).

Region. Exact rational forms Q_a, Q_b (stats16 config), A^2 = 2^21. Pool (common12.regenerate_pair): uniform on theta = 3/8 (a5 > 0, Q_a(a) <=
(3/8) A^2; b uniform in the half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a)). Region: a5 > 0, Q_a(a) <= (2565/8000) A^2; b5
>= 1, Q_b(b) <= (17/20) A^2 - Q_a(a) (exact), so the restricted sample is uniform on it. G has 28059 region rows (19513
primitive) of 87197 canonical draws, H 56008 (38828) of 174159; G + H is scaled by region volume over the effective
sample 1968000 (2565/3000)^3 = 1230051.906. Weights: w = (1 - qa)^3 (yield), vb (1 - qa)^3 x factor (work, population);
slice-geometry costs use S' = vb (17/20 - qa)^3 x the lower factor and K' = floor(sqrt((17/20 - qa) A^2 / Q_b(e_5))) lines. The
canonical a is the (Q_a, lex)-least element with a5 > 0 of its <+-1, units, sigma> class (exact enumeration of tO, t in {a, sigma a}); with b5 >= 1, one point per <+-1, units, tau> orbit (tau^2(a, b) = (a, -b)); non-primitive pairs are excluded by an exact HNF test. Counts ([3]). Points per volume of each half-ellipsoid, measured by the sampler's acceptance experiment (uniform integer points of the frame box of common12.uniform_half_ellipsoid, exact hit test; 2.5e11 trials per box, chunk k of box P = random.Random(1900000000 + 10^8 P + k), 10^7 trials): Q_a (pool box, every t/8000, t = 2300..3000): 0.999206..0.999331, bracket [0.999062, 0.999429]. Q_b (three boxes, 27,727 radii^2 in [A^2/2, A^2], step 2.5e-5; count nondecreasing, volume ~ radius^6): [0.999035, 0.999595] at every radius^2 in [A^2/2, A^2] (all pool slices (1 - qa) A^2 and region slices (17/20 - qa) A^2). Each point carries 6 SE (SE = ratio/sqrt(hits) <= 2.4e-5; union bound over 56,856 one-sided events < 8.7e-4), rounded outward. Factors: 0.998097 = 0.999062 x 0.999035 (yield, run cells, S') and 0.999595 = 1 x 0.999595 (work, population, D_SUP; a's end raised to 1, so each product also bounds b alone). Fallback: lattice_count_bounds.py certifies (exact rationals) [0.976706, 1.022997] (Q_a), [0.988895, 1.010232] (Q_b) at every radius >= 819; [3] checks the brackets inside it and replays 2^17 trials of two boxes.

Pools (PREREG-G). 656,000 uniform pairs each on theta = 3/8: S (seeds 261-264) sizes the region, G (251-254, drawn
after the configuration was fixed) and H (255-260, 265, 266: 1,312,000 pairs, pre-registered after G, every seed kept) are
judged, S and G give the work; no draw dropped. pool_d (stats16) holds per-seed counts (H: total), list hashes and sums, 8 anchor rows of G, G's canonical lists over draws 0..49 ([2] replays them) and the contamination sample. theta is the smallest t on the 1/8000 grid with margin_S
>= 0.06 + sqrt(2) SE_S (pre-registered), on S with the exact 2-adic yield: 2565/8000 (0.101067 >= 0.100756; 2564: 0.100595 < 0.100769; [4s]). The band
charge applies to the judged count only.

Sieve and survivors (f ideals only). Line sieve (scores only) for p in [2^14, 2^48); record walk for p in [2^48, 2^64) (slice12,
KERNELS.md; X2 = 2^64; per clean a-run, below; [2^63, 2^64) of contaminated runs on the wide path); primes in
[2^64, 2^69] by ECM (slack 302). Survivors pass the f threshold; exact norms (C_NORMS) give the g norm z. Primes below 2^14: the retention models' omitted mass. Survivors with z 2^69-smooth get trial division by every prime below 2^20 on both cofactors (82,025 primes, 2^23.649 <= C_TDIV = 2^24.5, [5t]). Batch band strip ([10b], KERNELS.md). Survivor sides wait in
waves of K = 2^39 (<= 832 bits). P = product of the primes in (2^20, 2^48] (< 1.000028 x 2^48/ln 2 bits, Rosser-Schoenfeld); a scaled remainder tree
(Bernstein; model16.srt) gives P mod z, the band part is gcd(z, (P mod z)^(2^6) mod z) (2^6 >= 41 >= any band exponent);
the f band-miss factor (1 - 41 x 0.001^2) stays as a margin. Records carry scores only; a g-passed f side gets
a second strip, P2 = the primes in (2^48, 2^64] (exponents <= 17 < 2^5; K = 2^54, C_BATCH2 = 2^25.25 >= 2^25.2401).
g batch ([10g]). Same SRT with P_g the product of all primes <= 2^69 (< 1.000028 x 2^69/ln 2 bits), waves of K = 2^61 g norms (<= 576 bits, 9 limbs): z is 2^69-smooth iff (P_g mod z)^(2^10) = 0 mod z (exponents <= 576 < 2^10). r_g >= 1 - 2^-38 (ideals of norm p^d > 2^69 drop a relation; in the margin). A kept relation's g norm is split by trial division, band split and ECM (<= 12 primes of (2^48, 2^69]): <= 12 x 1.002 x 2^39.79 + 13 C_MR = 2^43.38 per relation, <= 2^65.15 relations ([10t]). Its f band-2 part (<= 17 primes) is split by ECM until found (1.002 x 2^41.06 at 840 bits + 2^26
lookup each, + 18 C_MR(840)): with g 2^96.38 <= 2^96.5 units (ledger row).
Lemma (exact SRT). Product tree Z_v (L <= 61 levels); node v holds Y_v of b_v = 64 (n_v + G) bits (n_v limbs of
Z_v, G = L + 2). Root: Newton quotient q of P 2^b by Z, corrected by <= 2 add/subtracts, Y_root = q mod 2^b. Child c,
sibling s: Y_c = floor((Y_v Z_s mod 2^(b_v))/2^(b_v - b_c)). Leaf: r = ceil(z Y 2^-b) (0 if r = z) = P mod z. Proof:
eta_v = (frac(P/Z_v) - Y_v 2^(-b_v) mod 1) Z_v gains < 2^(-64G) per level, so eta < (L + 1) 2^(-64G) < 1 (no wrap). Cost (H-FFT; products <= model16.fft_F, gbatch_cost): C_BATCH = 2^24 (2^23.921, K = 2^39), C_GBATCH = 2^25 (2^24.973, K = 2^61). [10v]: model16.srt on 2^8
planted leaves (Newton offsets -2/+2) vs Python's %. P's and P2's builds (< 2^80 instructions) are misc, P_g's tree (< 2^88) the 2^74 row.
Band split ([10t]). A kept relation is split: the band parts are trial-divided by every prime of (2^20, 2^48] (< 2^43.272; 180 instructions per prime with a 2^(64 i) mod p table, misc row), each hit divided out to full power (<= 41 + 28 divisions, 2^10 each). Each found prime
gets its ideal and exponent from the factor base (<= 2^26): 2^50.7634 instructions per relation, no relation lost. Kept
relations <= G + H's 2^62.66658 + 2 SE over the count envelope
(0.999595/0.998097), the ECM miss (0.999^22: <= 13 + 9 primes above 2^64), the band margin and the least
retention on the model's norm grid (0.7775 f, 0.5393 g): <= 2^64.0271, 2^100.5028 <= 2^100.51 units (ledger row; audit).
Survivors' f omitted mass is < 310 bits; ECM is the ecm16.py schedule at 99.9% per
prime ([12]): 69-bit primes in a 310-bit cofactor cost 2^38.0909 instructions.
Cofactor windows (H-EA, model16.Window). A smooth side's cofactor C (< 2^320, no prime below 2^64 after division and both
strips) is 1, a prime <= 2^69, or composite with 64k < log2 C <= 69k, k in {2, 3, 4}; ECM runs only then. Compositeness: 40 Miller-Rabin rounds, random bases (C_MR = 40 x 2 x 310 five-limb
products = 2^23.121); a composite passes with probability <= 4^-40. ECM is charged at e x its schedule, e = the survivor law's mass whose part above 2^64 lies in a
window widened by 1 bit, minus its single-prime term, over the survivor mass (e_f 0.1553-0.1569 at
the anchors). g first: the f side's division, strip, test and ECM run only if the g batch passes:
per f-survivor C_NORMS + C_GBATCH + q_g (C_TDIV + C_BATCH + C_BATCH2 + C_MR + e_f 2^38.0909), q_g = min(1, 100 rho(u_g)) per pool row
(H1; 100 a margin); ~2^-22.5 f survivors per cell. S and G work sums use it on
all 1,312,000 draws (pool_row, runs16).

Thresholds. Certified per-block f bounds (embedding_bounds.py, 64-cell blocks, bisected to 8 bits; Taylor shift per place,
outward rounding) give T = lower - 16 x slack, keeping every B-smooth cell with omitted mass below slack; one evaluation (NormBounds.block) is what nbound31 replays bit for bit at 1,154.6 Ir. A line's rows
in an aligned 2 x 2 x 4 box (model16.box_rows) share one bisection per block over their hull. Each cell pays C_NB_EVAL x its box's evaluations per box cell, each row
its leaves, context and hull share (C_NB_ROW 535). NormBounds.row: per slice and place c_j = mid e_j, w_j = up((rad e_j + |c_j| 2^-49) UP);
a row takes m = sum b_j c_j, r = (sum |b_j| w_j)(1 + 2^-48): nine roundings < 2^-49.8 relative, so |b(v) - m| <= r (nbrow.c). [2b] checks 768 exact f norms of six boxes, [2] the anchors' fields.
Warm start (clean runs). A box-block's bisection starts from the leaves the run's first slice holding it certified
(model16.pool_row_warm): same tests and fallback, each leaf certified for its own a. Charged: warm evaluations (+ RUN_ALLOW x the cold excess), max(cold, warm) leaves; [2] and the audit replay the anchors.

Record walk per canonical a-run (rec7; KERNELS.md). A canonical a's run is the maximal integer t-interval
containing 0 with a + t e_0 in the region and canonical (runs16.py). On a run, a prime ideal (p, s, r)'s hits in all its slices form one coset of a 7-dim
lattice (t = r b(s) - a(s) mod p), walked once per run and ideal pair (p, s, +-r) with one setup6 ([10r], KERNELS.md).
Clean-run rule: rec7 walks a run only if all its a-values are in the canonical list (tested in the 2^92
row); otherwise its slices take the slice walk. A clean run of 32 to 2,048 a-values (price x REC7_R2 = 1.29 above 1,024) uses the run walk in each band
p in [2^48, 2^64) with h' = 2 x run cells x lower/upper ([3]) / p >= 2^4.99 (the table's bottom) instead of the slice walk and empty cosets; run
cells are exact (runs16.lrs: integer sum over the run's a-values of (17/20 - qa)^3). Each row pays its cells' share
(per hit at h', + 2 per record; setup6 per (run, pair) replaces its slices', f ideals). Contamination. A contaminated run's slices cost at most
the slice walk: per slice <= D_SUP = contam_sup = 2^72.564 on [2^48, 2^64) (exact band maxima with empty cosets and setup6;
[k7c]: 2001 region slices peak at 2^72.416). Every row pays RUN_ALLOW x D_SUP/2 per slice (f ideals), RUN_ALLOW = 0.00091 (draws are
uniform over canonical a): in 4,000 a-uniform runs (PREREG-C: seeds 5001-5010, the first 400 canonical a with Q_a <=
(2565/8000) A^2 per seed, sha256-16 in pool_d) the exact predicate found no non-canonical a-value among 1,964,062
(151,422 new, the rest inside runs tested at 2635/8000; its 48 runs above 1,024 a-values hold 52,095):
contaminated runs hold under 1 - 0.05^(1/4000) = 0.00074866 <= 0.00091 of the canonical a (95%). [k7c]/[2r] replay 6 sampled/3 anchor runs.
Coset sieve at P2. 2 = P2^2 P4^2, P2 = wO; b is in P2 iff P2F.b is even (common12.P2F; [10r] checks wO = ker P2F).
If a and b are in P2 the pair is not primitive, so those cells are not sieved; P2 is fixed by units, sigma, tau and
b -> -b, so the canonical domain is unchanged. In the basis e_0, e_j - (P2F.e_j) e_0 a slice with a in P2 sieves
b_0 = 2k + 1 (half the cells, congruence constants halved mod p, setup6 + 10). P2F(a + t) alternates along a run, so
rec7 walks it as two cosets (1/2 and 1/4 of its cells, setup6 + 10 each). model16.evaluate: p2a-and-p2b rows carry
nothing; p2a slices are priced at S/2 or the sub-run's cells.
S' uses the lower count factor ([3]), so the per-a charge is >= 0.00091 D_SUP >= 0.000749 D_SUP.
runs16's run ends (< 2^69.3 units) are misc.

Yield ([4], pools G + H). rels = (va_region/1230051.906) vb (0.998097) sum over G + H's canonical primitive region rows of
w rho(u_f) rho(u_g) r_f 0.999^lambda_f x the f band margin (g exact): Dickman rho ([1]), retentions r from the
Euler-product threshold model (threshold_model.py, calibrated on exact smooth integers at B = 2^10..2^14; prime_min
2^14), the ECM miss; cells with qa < 0.02 are sieved and charged, relations discarded. Estimate
2^62.66658 (SE 0.0363 bits) less the band charge 0.009: 2^62.6575 against 1.05 U' <= 2^62.52173: margin >= 0.1358
bits = 3.7 pool SE (>= max(0.06, 2 SE)). Band charge (H1'). On 4,000 region rows the unsieved
[2^8, 2^14) f content is 5.517 +- 0.133 bits (model 5.922), P(> 40 bits) <= 0.0033 (95%; model 0.0021): loss <= 0.0018
<= REL_BAND_LOSS = 0.009 bits ([4]). Descent input. LOG2_RELS = 62.583 is <= the charged count
lowered by two relative pool SE, 2^62.58314 ([4] rels_input); the descent's mean norms 739.5/498.3 truncate G + H's
739.534/498.389 ([4]).
Exact 2-adic part (H1, [14v], [14t]). rho prices a random integer (mean v_2 = 1); each row's v_s = v_2 of its exact norm
is known (pool_row), so u_s = (n_s + alpha_s + E_s + kappa_s - v_s)/69 (model16.a2, also q_g): E = (2, 2.399985) >= the
local means (runs16; alpha_2 = 1 - E), kappa = (0.06, 0.05) >= kappa2(u) (sum_k 2^-(k+1) rho(u - (k-1-kappa)/69) = rho(u),
rising in u) at n_f 840, n_g 530, above every yielding row. Cells with a in P2 carry exactly 2^6 in N_f, the rest 2^0.
Toy census (92,108 primitive (a, b) in {-1,0,1}^12, exact 2^20-smooth N_f): class chi^2 at 2^18/20/22 10/6/18 (exact) vs
88/30/23 (average), totals within 1 SE.

Cost ([5]). Sampled mean 35.90 instructions/kept cell over pools S and G, charged 36.27 with the allowance, over
2^124.220 kept cells (every step at KERNELS.md/STORAGE.md constants): collection <= 2^115.112719000
units (estimate 2^115.098120146 times allowance 1.01017). Allowance (H2, [4a]): a Monte Carlo estimate over S and G's
1,312,000 draws, relative SE row 0.2967%, incidence sqrt((1-56268/1312000)/56268) = 0.4124%, total <= 0.5081%;
charged 1.01017 >= 1 + 2 SE (a margin, not a probability).

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U', filtsim.PEXC), clique removal and Markowitz merge (order: (rows - 2) x lightest unit pivot)
are simulated on the production occurrence profile (filtsim.py; [6] reruns 2^12) as the tau-quotient system: tau^k(I) enters I's orbit column with
coefficient chi^k (entries a + b chi, merged |a|, |b| <= 255). Cap. All retained relations are used; H1 puts the count at >= 2^62.5831
>= 1.095 U' ([4]): filtering at excess 1.095, set after G + H's count (sizing and judgment stay at 1.05; margin over
1.095 U' 0.0753 >= max(0.06, 2 SE), [4]). Collection's ECM curves are fixed by a public seed, so
relations, filter and N are deterministic in the challenge-independent instance: H3 puts
N(400)/U' at <= 0.02687; the attack aborts if N > 0.02689 U' (never under H3; P_CAP = 1e-9 charged); LA is
charged at 0.02689. [10c]: per-run sd x sqrt(U) <= 0.648 at every scale with >= 6 seeds. A campaign at
means 744/501 (FILT_RUNS at 1.05, 6-16 seeds per scale 2^12..2^20; at excess 1.095 16 seeds at 2^19, 8 at 2^20) measures
x = N(400)/U' ([10] lists means and SE). H3: production is not above the largest scales; the pre-registered charge is the
larger mean + 2 SE at 1.095 of 2^19 and 2^20 -> 0.02687. At means 739.5/498.3 (<= G + H's) the same seeds give a pooled paired
ratio 0.9503 +- 0.0084 (+ 2 SE 0.9670 <= 1) (2^16, 2^17; FILT_Y, [10y]). Merge to mean weight 400 (charged 401; set from 2^16/2^17
before 2^19/2^20). N <= 0.02689 U' <= 2^57.2346 < 3 x 2^64/320: three owners (STORAGE.md).

The square system A = [M | S] has N + 10 rows, N orbit columns and 10 dense Schirokauer columns: the chi-eigenpart of
the unit groups, (rk K_f - rk K_f^+ + rk K_g - rk K_h)/2 = (35 - 17 + 5 - 3)/2 from signatures K_f (0, 36), K_f^+ (0, 18),
K_g (0, 6), K_h (2, 2) ([13]; tau^2 = -1 on the sigma-odd part, so its chi and chi^-1 eigenspaces have equal dimension);
under H3 its kernel is the line of virtual logs (w, c), c != 0, and M has rank N. Block Wiedemann (m, n) = (1024, 30) with S as right-hand sides
(as CADO-NFS): B = [M | 0]; Y = (S_1..S_10, B z_11..B z_n), z uniform in F_r; X block-diagonal (H-X): the 10 pad coordinates (zero columns of B)
in blocks 1..10, the others shuffled and cut into m blocks, X_l uniform in [0, 2^128) on block l. A generator sum_k B^k Y f_k = 0 of
a_i = X^T B^i Y with f_0 = (c, e), c != 0, gives S c + B w = 0, w = Z e + sum_(k >= 1) B^(k-1) Y f_k: one Horner pass, no product touches S
([10x]: toy end to end). Krylov N/m + N/n block steps; Horner runs as five interleaved chains from checkpoints
B^(jL) Y (L = Horner length / 5; 4n = 120 checkpoint vectors, < 2^70 bytes, 4L < N/m + N/n) at the five-vector rate. Each
vector element is stored with chi x_j. Rows reduce exactly (_row_red). Per row and vector (STORAGE.md, [10]): 401 x 8.38 (H4u, PREREG-U) + 500 (row tail, three owners) + 200 (_chi_mul); per
Horner element n generator terms and one reduction (270); projections: each vector element meets one
X entry, (1 + n/m) N^2 terms at 60. Linear algebra <= 2^112.313718742 units; generator (matrix Berlekamp-Massey by half-gcd, k = m + n = 1054: a node of l
sequence terms makes (m + k) k^2 products of polynomials with <= l and <= l/2 + 1 coefficients, each one Kronecker integer product in 9-limb
slots (coefficients < k L r^2 < 2^576): F(18 l + 18) under H-FFT, + 640 l to pack, reduce (270) and add; leaves m k^2 (97 + 270). The sum over
the recursion is below model16.generator ([10m]): 4 k^3 log2 L (F(18 L) + 640 L) instructions, L = (N + 10)(1/m + 1/n) + 2): <= 2^94.846 units.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent; with seeded ECM curves the relation count is fixed:
2^62.6575 under H1 (section 3). Under H3 (section 9; the 10 Schirokauer maps of the chi-eigenpart
injective on units modulo r-th powers) Block Wiedemann fails with probability <= 4N/2^128 < 2^-66 (Schwartz-Zippel: a minor of order N + 9 of
(X, B^T X, ..)^T (Y, B Y, ..), degree <= 2N + 18 in the entries of X and z, uniform on >= 2^128 values; H-X); a failed solve (checked on the dropped rows) is a failure.

Reconstruction. Purged ideals' logs come by back-substitution; ideals with <= 2 occurrences are log-less in the descent
([6b]: model >= 6.3x measured; DESCENT.md).

## 6. Individual logarithms

DESCENT.md ([7], [7e], [7l]): 62 nodes + boot (215, 304, 15, 4), retention >= 0.6490.

Expected two-log work <= 2^109.902 units ([7]; unrounded 2^109.90176; one log 2^123.189 instructions: boot own
2^120.45 + subtree 2^121.68 x retry 2.416). The meter caps both logs at t = 2^121.9 units (reached = failure). Metered work is min(D, t) <= D for the uncapped two-log work D: its expectation is <= 2^109.902 for any t, and
P(cap) <= E[D]/t <= 2^(109.902 - 121.9) <= 0.0002445 (Markov); with 10% more expected work, <= 0.0002690.

## 7. Success

Core: success >= (1 - 2^-66 - P_CAP) x (1 - 2 x 6.649e-05 - 0.0002445 - 2^-11.61 - 2^-10.70 - 2/r) >= 0.9987 (main
phase, two boots with all 18 relations failing, cap, node counter, repeats); >= 0.998677 + 2^-48 with the cap bound at +10% (0.0002690).

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX = 2^27.3. E[N] uses the same
  recurrence with unit node cost (counts16.py): E[N] <= 2^15.6846 for both logs, P(abort) <= 2^-11.61.
- Repeats: N_MAX nodes in the smallest bin (2^70.30 ideals), 2^6 for 1/p weighting: P <= 2^-10.70, charged.

Subtrees are independent (H5). Past the gate P(success) >= 0.998677 with the wrapper (the +10% figure); every abort and a closed gate fail: P(success) >= (K/2^64) x 0.998677 >= 0.39. W <= 1 + 1[u < K] W_run, u independent of W_run: E[W] <= 1 + (K/2^64) E[W_run] (audit: gate, expected_run).

## 8. Ledger

| component | log2 core units (20,000 instructions) |
|---|---:|
| relation collection (section 3, H-EA; estimate 115.098120146 + 1.017% allowance) | <= 115.112719000 |
| linear algebra (weight 401, N <= 2^57.2346, (1024, 30), block-diagonal X; section 4) | <= 112.313718742 |
| matrix generator | <= 94.846 |
| both individual logarithms (metered cap; expected 109.902) | 121.9 |
| canonical a list | 92 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 73 |
| wide factor-base lattices ([2^63, 2^64), A and coset, KERNELS.md) | 71.6 |
| rec7 bases (2 per ideal x run and two cosets, p in [2^48, 2^64), KERNELS.md) | 75.09 |
| band split (section 3, [10t]) | 100.51 |
| P_g product tree (section 3) | 74 |
| g and f band-2 split of kept relations (section 3) | 96.5 |
| numerical certificate construction | 50 |
| verification | 11 |
| wrapper (section 1; <= 2^80 instructions) | 66 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^92 (2^42.8 per region a, 2^62.984 of them); Schirokauer maps 2^88 (below); factor base 2^73 (2^23.4 per prime, pi(2^69) < 2^63.46); filtering and
back-substitution 2^80 (2^24.4 per entry against < 2^14). The rows sum to <= 2^121.91486223534912: core meter 2^121.9148622354 (audit-checked); with the two logarithms at 2^109.902, not the cap, <= 2^115.3401320190.

Schirokauer maps. R = (Z/r^2)[y, x]/(h, f) (rank 72; g: rank 12); r is unramified ([13]): eps = lcm(r^d_i - 1) < r^72
(18,360 bits; g: 3,060). The 10 columns are fixed F_r-linear forms (chi-eigenpart) of (phi^eps - 1)/r mod r, phi = a - bx,
by left-to-right powering on eight-limb coefficients. An f squaring: 2,628 = 72 x 73/2 schoolbook 8 x 8-limb products
(<= 512 instructions), reduction by monic h, f <= 2^19, 72 Barrett (<= 1,024): < 2^20.9; times phi (12 coefficients
< 2^64) < 2^19.3. Per bit < 2^21.31: f < 2^35.5, g < 3,060 x 2^17, so < 2^36 per relation, 2^86.86 units at 2^65.15
relations, within 2^88.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (0.29, 0.75), exact v_2 and kappa, exact sampled norms). H1' retention (Euler-product model; band charge 0.009 bits). H-EA (work only): H1''s survivor law prices the ECM share ([12w]). H-FFT (work only): F(N) bounds the local SSA work past 2^26 limbs (KERNELS.md). H2 region (uniform pools, exact flags, measured counts ([3], 6 SE; certified bracket as fallback); sized on S, judged on G + H; S+G work with the 1.017% allowance). H3 matrix (r unramified, [13]; r prime to the class numbers, the unit condition, a one-dimensional kernel; N(400)/U'
<= 0.02687, section 4). H4 kernel cost (KERNELS.md, STORAGE.md); H4u (unit share >= 86%;
k5 <= 8.38, PREREG-U, [10]); rec7's per-hit price at production stays below the interpolant of the measured maxima (REC7_TAB); contaminated runs hold at most 0.091%
of the canonical a (95% bound; a shortfall changes cost only); f/g independence given q_g (100x). H-X: BW's minor (section 5; N + 9 = dim(span S + Im B)) is not identically zero for this X and Y; m, n >= 10 ([10x]: toys solve end to end). The pad pinning is needed: x^T B^i (i >= 1) vanishes on ker B (the pad space), so only X itself separates it; with all coordinates shuffled the solve fails when the pad coordinates meet <= 8 blocks. H5 descent (DESCENT.md): lower sums; per-ideal yield certificate (DISC); lazy order independent of outcomes; negative-binomial recurrence, independent
subtrees of distinct ideals (measured: [7e]); missing logs (Poisson occurrence at the 2-SE input); local omitted-mass laws; random-integer ECM orders; UDUP = 0.65 (DESCENT.md, [11]).

## 10. Memory

Core storage <= 2^105 bytes (STORAGE.md), wrapper below 2^30: below 2^106 bytes. All preprocessing is for this target,
charged once in W.

## Provenance

GT core: kaankolcu recover-P 1b8fa720 and its lineage (9332a74d mpjunior92; 2c2ce2ab / 42e80f80 gbotrel; 6282ab79
bernardo-lima12; 6734bd5a GPT-6 in Codex; Gemini 45ce22e5 / 4bb6c381); polynomials after Guillevic-Singh; filtering after CADO-NFS. H-EA, H-X,
rec7, band strip and split, f-only sieve, g batch, twisted C4 field: gbotrel, Claude Opus 5.5 in Claude Code. Run gate, expected total work, high cap: kaankolcu (8ab78e96, b6f0a0f6, d7ca6b98, f1c17abf); logarithms at their expectation: also mpjunior92 (ea928ab4).

# Plain G1 ECDLP (bls12-381) through two GT logarithms with the STNFS core: 2^118.969870492892 units

Track bls12-381-g1-ecdlp-plain. Input: T = [x]P0 in G1; output the scalar x in 0..r-1.
Score unit: one batched affine G1 addition, C = 2,000 scalar AArch64 instruction equivalents. The core meters
core units of 20,000 instructions (UNIT) = 10 score units: log2 W = log2(meter) + log2 10. No auxiliary powers, advice, smaller targets or amortization
are used. No prior verdict is a premise.

## Claim

The algorithm below returns x with success at least 0.39 over uniform
challenges and fresh coins under the declared heuristics H1-H5, H-EA, H-X and H-FFT (section 9). All classical work (preprocessing, both logarithms, failures,
restarts, field map, both pairings, scalar multiplication, verification, output) is at most

    W <= 10 x (2^115.647942398004 + 2^66) < 2^118.969870492892  units (C = 2,000).

The score is 118.969870492892 (core: sections 2-8; wrapper: section 1). Changes vs 7fd7a21c: descent child load (DESCENT.md), FFT generator, BW (512, 25), matrix cap (section 4). Storage is below 2^106 bytes. Analytic, exploratory claim.

Conventions. Bounds round safely; evidence16.py checks them against unrounded values and reruns checks [n]. Model: stats16.b85, plan16.b85, model16.py, descent16.py (owner10.install: C_ROW 700).

## 1. Reduction (fixed base g = e(P0, Q0): two GT logarithms) and field map

With the pinned pairing e and public P0, Q0:

1. Check T (on curve, [r]T = infinity). If T = infinity return x = 0.
2. g = e(P0, Q0), z = e(T, Q0) (pinned pairings, protected encoding).
3. Map g and z into the core's tower (field-map.md).
4. Run the core on (g, z) under its meter and node counter (any failure is a failure): x = log_g(z) mod r.
5. Verify: return x if [x]P0 = T, otherwise fail.

GT has prime order r, e is nondegenerate: e(., Q0): G1 -> GT is an isomorphism, z = g^x, z = 1 iff T = infinity. Each core target t becomes t^e, e fresh uniform in [1, r) (DESCENT.md), uniform on GT minus {1} for any fixed g, z: the fixed base costs no success; the final check rejects wrong outputs.

One relation/LA phase serves log z and log g (section 6; 2/r charged). Wrapper below 2^80 instructions (both pairings, scalar multiplication, check < 2^50, field map < 2^67): 2^65.71 core units, a 2^66 row. Wrapper failure: field map below 2^-48.

## 2. Tower, polynomials, automorphism ([0], [0b], [0c], [14])

K_h = Q[y]/h, h = y^6 + 2y^4 - 2y^3 - 2y^2 - 1, irreducible mod p. Disc(h) = 2^13 x 19^2; O_h has no power basis: common12.py fixes a Z-basis e_0 = 1, ..., e_5 (index 4 over Z[y]); [0c] certifies maximality (19: Dedekind; 2: runs16.py, 2 O_h = (w)^2 (pi4)^2, norms 2, 4): Disc(O_h) = 2^9 x 19^2. With w = y + 1 (N(w) = 2), f = w^6 F(x^2/w) = x^12 - 2w x^10 + 2w^3 x^6 + w^5 x^2 + w^6 (F the BLS12 polynomial) and g = x^2 + |u| w, u = -0xd201000000010000. [0b] checks (u w)^((p^6-1)/2) = -1 in F_p6, so F_p[y,x]/(h, g) = F_p12; [0] checks Res_y(h, f(x^2 -> u w)) = 0 mod p.
Twisted symmetry. sigma(y) = -y^5 - 2y^3 + 2y^2 + 2y (order 2), c = (y^5 + 3y^3 - 3y^2 - 3)/2 a unit with c sigma(c) = -1; tau: (y, x) -> (sigma y, c x) maps f to c^12 f and g to c^2 g ([0b]), tau^2 = (x -> -x). On F_p12 tau is Frobenius^3 ([0b]: sigma(y) = y^(p^3), (u w)^((p^3-1)/2) = c), so log tau(P) = chi log P with chi = p^3 mod r, chi^2 = -1. One relation stands for its tau-orbit (4 pairs); the factor base is quotiented by tau-orbits of 4 ideals. Murphy alpha charged at (0.29, 0.75); [14] recomputes it (Kummer-Dedekind for odd p, exact local laws at 2: f -1, g -1.39998): (-0.047244, 0.706684) over p <= 2000, (0.285135, 0.745775) over p < 50000, both below.

## 3. Relation collection

Candidates: (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69; factor base: relative-degree-1 ideals of
norm <= B, U' = 2^62.451 tau-orbits (U/2: orbits of 4 ideals; U = 2^63.451).

Region. Exact rational forms Q_a, Q_b (stats16 config), A^2 = 2^21. Pool: uniform on theta = 3/8 (a5 > 0, Q_a(a) <= (3/8) A^2; b5 >= 1, Q_b(b) <= A^2 - Q_a(a)). Region: a5 > 0, Q_a(a) <= (2635/8000) A^2; b5 >= 1, Q_b(b) <= (17/20) A^2 - Q_a(a) (exact). G has 30123 region rows (20951 primitive) of 87197 canonical draws, H 60386 (41831) of 174159; G + H scaled by region volume over effective sample 1968000 (2635/3000)^3 = 1333531.223. Weights: w = (1 - qa)^3 (yield), vb (1 - qa)^3 x factor (work, population); slice costs use S' = vb (17/20 - qa)^3 x factor and K' = floor(sqrt((17/20 - qa) A^2 / Q_b(e_5))) lines. Canonical a is the (Q_a, lex)-least element with a5 > 0 of its <+-1, units, sigma> class; with b5 >= 1, one point per <+-1, units, tau> orbit; non-primitive pairs excluded by exact HNF. lattice_count_bounds.py certifies count ratio at radius R >= 831 = isqrt((2635/8000) 2^21), below smallest b-slice radius 1044 ([3]): [0.977040, 1.022662] (Q_a), [0.989055, 1.010083] (Q_b); model uses products 0.966347073 (yield) and 1.032972622 (work).

Pools (PREREG-G). 656,000 uniform pairs each on theta = 3/8: S (seeds 261-264) sizes the region; G (251-254, drawn after the configuration) and H (255-260, 265, 266: 1,312,000 pairs, pre-registered after G) are judged; S and G give the work; no draw dropped. pool_d holds per-seed counts (H: total), list hashes and sums, 8 anchor rows of G, G's canonical lists over draws 0..49 ([2] replays them) and the contamination sample. theta is the smallest t on the 1/8000 grid with margin_S >= 0.06 + sqrt(2) SE_S (pre-registered), on S with exact 2-adic yield: 2635/8000 (0.100072 >= 0.099549; 2634: 0.099258 < 0.099569; [4s]). The band charge applies to the judged count only.

Sieve and survivors (f ideals only). Line sieve (scores only) for p in [2^14, 2^48); record walk for p in [2^48, 2^64) (slice12, KERNELS.md; X2 = 2^64; per clean a-run, below; [2^63, 2^64) of contaminated runs on the wide path); primes in [2^64, 2^69] by ECM (slack 302). Survivors pass the f threshold; exact norms (C_NORMS) give the g norm z. Primes below 2^14 are the retention models' omitted mass. Survivors with z 2^69-smooth get trial division by every prime below 2^20 on both cofactors (82,025 primes, 2^23.649 <= C_TDIV = 2^24.5, [5t]). Batch band strip ([10b], KERNELS.md): survivor sides wait in waves of K = 2^39 (<= 832 bits). P = product of the primes in (2^20, 2^48] (< 1.000028 x 2^48/ln 2 bits, Rosser-Schoenfeld); a scaled remainder tree (Bernstein; model16.srt) gives P mod z; the band part is gcd(z, (P mod z)^(2^6) mod z) (2^6 >= 41 >= any band exponent); the f band-miss factor (1 - 41 x 0.001^2) stays as margin. Records carry scores only; a g-passed f side gets a second strip, P2 = the primes in (2^48, 2^64] (exponents <= 17 < 2^5; K = 2^54, C_BATCH2 = 2^25.25 >= 2^25.2401). g batch ([10g]): same SRT with P_g the product of all primes <= 2^69 (< 1.000028 x 2^69/ln 2 bits), waves of K = 2^61 g norms (<= 576 bits, 9 limbs): z is 2^69-smooth iff (P_g mod z)^(2^10) = 0 mod z (exponents <= 576 < 2^10). r_g >= 1 - 2^-38 (ideals of norm p^d > 2^69 drop a relation; in the margin). A kept relation's g norm is split by trial division, band split and ECM (<= 12 primes of (2^48, 2^69]): <= 12 x 1.002 x 2^39.79 + 13 C_MR = 2^43.38 per relation, <= 2^65.15 relations ([10t]). Its f band-2 part (<= 17 primes) is split by ECM until found (1.002 x 2^41.06 at 840 bits + 2^26 lookup each, + 18 C_MR(840)): with g 2^96.38 <= 2^96.5 units (ledger row).
Lemma (exact SRT). Product tree Z_v (L <= 61 levels); node v holds Y_v of b_v = 64 (n_v + G) bits (n_v limbs of Z_v, G = L + 2). Root: Newton quotient q of P 2^b by Z, corrected by <= 2 add/subtracts, Y_root = q mod 2^b. Child c, sibling s: Y_c = floor((Y_v Z_s mod 2^(b_v))/2^(b_v - b_c)). Leaf: r = ceil(z Y 2^-b) (0 if r = z) = P mod z: eta_v = (frac(P/Z_v) - Y_v 2^(-b_v) mod 1) Z_v gains < 2^(-64G) per level, so eta < (L + 1) 2^(-64G) < 1 (no wrap). Cost (H-FFT): C_BATCH = 2^24 (2^23.921, K = 2^39), C_GBATCH = 2^25 (2^24.973, K = 2^61). [10v]: model16.srt on 2^8 planted leaves (Newton offsets -2/+2) vs Python's %. P's and P2's builds (< 2^80 instructions) are misc, P_g's tree (< 2^88) the 2^74 row.
Band split ([10t]). A kept relation is split: band parts trial-divided by every prime of (2^20, 2^48] (< 2^43.272; 180 instructions per prime with a 2^(64 i) mod p table, misc row), each hit divided out to full power (<= 41 + 28 divisions, 2^10 each). Each found prime gets its ideal and exponent (<= 2^26): 2^50.7634 instructions per relation, no relation lost. Kept relations <= G + H's 2^62.66428 + 2 SE over the yield envelope (1.032972622/0.966347073), the ECM miss (0.999^22: <= 13 + 9 primes above 2^64), the band margin and least retention on the norm grid (0.7775 f, 0.5393 g): <= 2^64.1165, 2^100.5922 <= 2^100.62 units (ledger row; audit). Survivors' f omitted mass is < 310 bits; ECM is the ecm16.py schedule at 99.9% per prime ([12]): 69-bit primes in a 310-bit cofactor cost 2^38.0909 instructions.
Cofactor windows (H-EA, model16.Window). A smooth side's cofactor C (< 2^320, no prime below 2^64 after division and both strips) is 1, a prime <= 2^69, or composite with 64k < log2 C <= 69k, k in {2, 3, 4}; ECM runs only then. Compositeness: 40 Miller-Rabin rounds, random bases (C_MR = 40 x 2 x 310 five-limb products = 2^23.121); a composite passes with probability <= 4^-40. ECM is charged at e x its schedule, e = the survivor law's mass above 2^64 in a window widened by 1 bit, minus its single-prime term, over the survivor mass (e_f 0.1553-0.1569 at the anchors). g first: the f side's division, strip, test and ECM run only if the g batch passes: per f-survivor C_NORMS + C_GBATCH + q_g (C_TDIV + C_BATCH + C_BATCH2 + C_MR + e_f 2^38.0909), q_g = min(1, 100 rho(u_g)) per pool row (H1; 100 a margin); ~2^-22.5 f survivors per cell. S and G work sums use it on all 1,312,000 draws.

Thresholds. Certified per-block f bounds (64-cell blocks, bisected to 8 bits; Taylor shift per place, outward rounding) give T = lower - 16 x slack, keeping every B-smooth cell with omitted mass below slack; nbound31 replays one NormBounds.block evaluation bit for bit at 1,154.6 Ir. Rows in an aligned 2 x 2 x 4 box share one bisection per block over their hull. Each cell pays C_NB_EVAL x its box's evaluations per box cell, each row its leaves, context and hull share (C_NB_ROW 535). NormBounds.row: per slice and place c_j = mid e_j, w_j = up((rad e_j + |c_j| 2^-49) UP); a row takes m = sum b_j c_j, r = (sum |b_j| w_j)(1 + 2^-48): nine roundings < 2^-49.8 relative, so |b(v) - m| <= r (nbrow.c). [2b] checks 768 exact f norms of six boxes, [2] the anchors' fields.
Warm start (clean runs): a box-block's bisection starts from the leaves certified by the run's first slice holding it (model16.pool_row_warm), same tests and fallback, each leaf certified for its own a. Charged: warm evaluations (+ RUN_ALLOW x the cold excess), max(cold, warm) leaves; [2] and the audit replay the anchors.

Record walk per canonical a-run (rec7; KERNELS.md). A canonical a's run: the maximal t-interval containing 0 with a + t e_0 in the region and canonical (runs16.py). On a run, a prime ideal (p, s, r)'s hits in all its slices form one coset of a 7-dim lattice (t = r b(s) - a(s) mod p), walked once per (run, ideal) with one setup6 ([10r]). Clean runs only (all a-values canonical, tested in the 2^92 row); otherwise the slice walk. A clean run of 32 to 1,024 a-values uses the run walk in each band p in [2^48, 2^64) with h' = run cells / p >= 2^4.99 (the table's bottom) instead of the slice walk and empty cosets; run cells are exact (integer sum over the run's a-values of (17/20 - qa)^3). Each row pays its cells' share (per hit at h', + 2 per record; setup6 per (run, ideal) replaces its slices', f ideals). Contamination: a contaminated run's slices cost at most the slice walk, per slice <= D_SUP = contam_sup = 2^72.591 on [2^48, 2^64) (exact band maxima; [k7c]: 2001 region slices peak at 2^72.446). Every row pays RUN_ALLOW x D_SUP/2 per slice (f ideals), RUN_ALLOW = 0.00091: in 4,000 a-uniform runs (PREREG-C: seeds 5001-5010, first 400 canonical a with Q_a <= (2635/8000) A^2 per seed, sha256-16 in pool_d) no non-canonical a-value among 1,989,563 (166,465 new): contaminated runs hold under 1 - 0.05^(1/4000) = 0.00074866 <= 0.00091 of the canonical a (95%). [k7c]/[2r] replay 6 sampled/3 anchor runs.
Coset sieve at P2. 2 = P2^2 P4^2, P2 = wO; b is in P2 iff P2F.b is even ([10r] checks wO = ker P2F). Cells with a, b in P2 are not primitive and not sieved; P2 is fixed by units, sigma, tau, b -> -b. In the basis e_0, e_j - (P2F.e_j) e_0 a slice with a in P2 sieves b_0 = 2k + 1 (half the cells, congruence constants halved mod p, setup6 + 10). P2F(a + t) alternates along a run, so rec7 walks it as two cosets (1/2 and 1/4 of its cells, setup6 + 10 each). p2a-and-p2b rows carry nothing; p2a slices are priced at S/2 or the sub-run's cells. Count envelopes ([3]) keep the per-a charge >= 0.000861 D_SUP >= 0.000749 D_SUP. Run ends (< 2^69.3 units) are misc.

Yield ([4], pools G + H). rels = (va_region/1333531.223) vb (0.966347073) sum over G + H's canonical primitive region rows of w rho(u_f) rho(u_g) r_f 0.999^lambda_f x the f band margin (g exact): Dickman rho ([1]), retentions r from the Euler-product threshold model (calibrated on exact smooth integers at B = 2^10..2^14; prime_min 2^14), the ECM miss; cells with qa < 0.02 are sieved and charged, relations discarded. Estimate 2^62.66428 (SE 0.0352 bits) less the band charge 0.009: 2^62.6552 against 1.05 U' <= 2^62.52173: margin >= 0.1335 bits = 3.8 pool SE (>= max(0.06, 2 SE)). Band charge (H1'): on 4,000 region rows the unsieved [2^8, 2^14) f content is 5.511 +- 0.129 bits (model 5.922), P(> 40 bits) <= 0.0042 (95%; model 0.0021): loss <= 0.0032 <= REL_BAND_LOSS = 0.009 bits ([4]). Descent input: LOG2_RELS = 62.583 is <= the charged count lowered by two relative pool SE, 2^62.58313 ([4] rels_input); the descent's mean norms 740.0/498.4 truncate G + H's 740.096/498.416 ([4]).
Exact 2-adic part (H1, [14v], [14t]). rho prices a random integer (mean v_2 = 1); each row's exact v_s is known, so u_s = (n_s + alpha_s + E_s + kappa_s - v_s)/69: E = (2, 2.399985) >= the local means (alpha_2 = 1 - E), kappa = (0.06, 0.05) >= kappa2(u) (sum_k 2^-(k+1) rho(u - (k-1-kappa)/69) = rho(u), rising in u) at n_f 840, n_g 530, above every yielding row. Cells with a in P2 carry exactly 2^6 in N_f, the rest 2^0. Toy census ([14v]): class chi^2 at 2^18/20/22 10/6/18 (exact) vs 88/30/23 (average), totals within 1 SE.

Cost ([5]). Stored sums price CELL_ROUTE at 0.5 and REC_ROUTE at 10.01. The owner now reads records in place from handed-off pages (STORAGE.md, route_inplace.txt, [10i]), dropping the 3-instruction page copy: REC_ROUTE is 7.01. Subtract 0.48 + 0.8630462173 per counted cell from the first moment only (not skip rows; the route term is uniform per cell: 3.0 x 2 ln(64/48) x 0.5 records). Second moment unchanged, so SE 0.51580% is an upper bound. Mean 40.4714, charged 40.8895 at allowance 1.01033>=1+2*0.00515793 (need 1.01031585), over 2^124.339 cells: collection <= 2^115.405307. CELLRES slope 0.00175232; 0.02>4*slope.

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U', filtsim.PEXC), clique removal and Markowitz merge (order: (rows - 2) x lightest unit pivot) are simulated on the production occurrence profile (filtsim.py; [6] reruns 2^12) as the tau-quotient system: tau^k(I) enters I's orbit column with coefficient chi^k (entries a + b chi, merged |a|, |b| <= 255). Cap: all retained relations are used; H1 puts the count at >= 2^62.5831 >= 1.095 U' ([4]): filtering at excess 1.095, set after G + H's count (sizing and judgment stay at 1.05; margin over 1.095 U' 0.0730 >= max(0.06, 2 SE), [4]). Collection's ECM curves are fixed by a public seed, so relations, filter and N are deterministic in the challenge-independent instance: H3 puts N(400)/U' at <= 0.0274; the attack aborts if N > 0.0274 U' (a failure; P_CAP = 0.05 is charged as a model-error allowance, not derived from a tail law); LA is charged at 0.0274. [10c]: per-run sd x sqrt(U) <= 0.648 at every scale with >= 6 seeds. A campaign at means 744/501 (FILT_RUNS at 1.05, 6-16 seeds per scale 2^12..2^20; at excess 1.095 16 seeds at 2^19, 8 at 2^20) measures x = N(400)/U' ([10] lists means and SE). H3: production is not above the largest scales. x falls with scale (1.05: 0.03305 at 2^14 to 0.02830 at 2^20; 1.095: 0.02668 +- 0.00009 at 2^19, 0.02567 +- 0.00018 at 2^20); 0.0274 >= 1.019 x the larger 1.095 mean + 2 SE (0.02687), 6.8% above the 2^20 mean ([10c] checks both). At means 740.0/498.4 (<= G + H's) the same seeds give a pooled paired ratio 0.9524 +- 0.0073 (+ 2 SE 0.9670 <= 1) (2^16, 2^17; FILT_Y, [10y]). Merge to mean weight 400 (charged 401; set from 2^16/2^17 before 2^19/2^20). N <= 0.0274 U' <= 2^57.2617 < 3 x 2^64/320 = 2^57.2630: three owners (STORAGE.md).

The square system A = [M | S] has N + 10 rows, N orbit columns and 10 dense Schirokauer columns: the chi-eigenpart of the unit groups, (rk K_f - rk K_f^+ + rk K_g - rk K_h)/2 = (35 - 17 + 5 - 3)/2 from signatures K_f (0, 36), K_f^+ (0, 18), K_g (0, 6), K_h (2, 2) ([13]; tau^2 = -1 on the sigma-odd part, so its chi and chi^-1 eigenspaces have equal dimension); under H3 its kernel is the line of virtual logs. Block Wiedemann (m, n) = (512, 25): Y uniform in F_r; X block-diagonal (H-X): shuffled coordinates cut into m blocks, X_l uniform in [0, 2^128) on block l. Krylov N/m + N/n block steps (5 five-vector groups in turn); Horner runs as five interleaved chains from checkpoints A^(jL) Y (L = Horner length / 5; 4n = 100 checkpoint vectors, < 2^69 bytes, 4L < N/m + N/n) at the five-vector rate. Each vector element is stored with chi x_j. Rows reduce to 2^-320 x value (BW on 2^-320 A). Per row and vector (STORAGE.md, [10]): 401 x 8.38 (H4u, PREREG-U) + 700 (row tail 695, three owners) + 10 x 97 + 198 (_chi_mul); per Horner element n dense terms and one reduction (266); projections: each vector element meets one X entry, (1 + n/m) N^2 terms at 60. Linear algebra <= 2^112.7571 units; generator (matrix Berlekamp-Massey, recursive doubling, k = 537, L = (N + 10)(1/m + 1/n) + 2; model16.generator: 4 k^3 log2 L products, Kronecker in 9-limb slots (entries < 2^257, L 2^514 < 2^568, exact): fft_F(18L) + L (54 + 4 (266 + 97)); slot s: REDC(s mod 2^512) + REDC((s >> 512)(2^512 mod r)), inputs < 2^516; 18L = 2^56.9 limbs, H-FFT as P_g): <= 2^87.727 units.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent; with seeded ECM curves the relation count is fixed: 2^62.6552 under H1 (section 3). Under H3 (section 9; the 10 Schirokauer maps of the chi-eigenpart injective on units modulo r-th powers) Block Wiedemann fails with probability <= 4N/2^128 < 2^-66 (Schwartz-Zippel, degree <= 2N for any (m, n), entries uniform on >= 2^128 values; H-X); a failed solve (checked on the dropped rows) is a failure.

Reconstruction. Purged ideals' logs come by back-substitution; ideals with <= 2 occurrences are log-less in the descent ([6b]: model >= 6.3x measured; DESCENT.md).

## 6. Individual logarithms

DESCENT.md ([7], [7e], [7l], [7m]): 62 nodes + boot (215, 304, 15, 4), retention >= 0.6490.

Expected two-log work <= 2^108.187 units ([7]; unrounded 2^108.18602; one log 2^121.474 instructions: boot own 2^120.19 + subtree 2^119.71 x retry 2.011; boot psi 0.5028). The meter caps both logs at 2^110 units (reached = failure): P(cap) <= E[W]/cap <= 2^(108.187 - 110) <= 0.284599 (Markov).

## 7. Success

Core: success >= (1 - 2^-66 - P_CAP) x (1 - 2 x 4.21e-06 - 0.284599 - 2^-12.27 - 2^-10.70 - 2/r) >= 0.6788 (main
phase, two boots failing all 18 relations, Markov, node counter, repeats); >= 0.4084 at 2x the modeled descent work.

- Node counter: abort (failure) before node N_MAX + 1, N_MAX = 2^27.3. E[N]: slot-count recurrence, unit node cost, these retries (counts16.py): <= 2^15.0258 for both logs, P(abort) <= 2^-12.27.
- Repeats: N_MAX nodes in the smallest bin (2^70.30 ideals), 2^6 for 1/p weighting: P <= 2^-10.70.

With the wrapper: P(success) >= 0.6788 - 2^-48 > 0.39; every abort fails.

## 8. Ledger

| component | log2 core units (20,000 instructions) |
|---|---:|
| relation collection (section 3; first moment minus 1.34304 per cell, allowance 1.01033) | <= 115.405307 |
| linear algebra (section 4) | <= 112.7571 |
| matrix generator (section 4) | <= 87.727 |
| both individual logarithms (metered cap, section 6) | 110 |
| canonical a list | 92 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 73 |
| wide factor-base lattices (KERNELS.md) | 71.6 |
| rec7 bases (KERNELS.md) | 75.09 |
| band split (section 3, [10t]) | 100.62 |
| P_g product tree (section 3) | 74 |
| g and f band-2 split (section 3) | 96.5 |
| numerical certificate construction | 50 |
| verification | 11 |
| wrapper (section 1) | 66 |

Fixed rows (units x 2^14.29 = instructions): canonical 2^92 (2^42.8 per region a, 2^63.100 of them); Schirokauer maps 2^88 (below); factor base 2^73 (2^23.4 per prime, pi(2^69) < 2^63.46); filtering, back-substitution 2^80 (2^24.4 per entry, < 2^14). Rows sum to <= the core meter 2^115.647942398004; with the wrapper, x 10 < 2^118.969870492892 units (audit).

Schirokauer maps. R = (Z/r^2)[y, x]/(h, f) (rank 72; g 12), r unramified ([13]): eps = lcm(r^d_i - 1) < r^72 (18,360 bits; g 3,060); the 10 columns are F_r-linear forms (chi-part) of (phi^eps - 1)/r, phi = a - bx, by left-to-right powering on eight-limb coefficients. f squaring: 2,628 = 72 x 73/2 8x8-limb products (<= 512), reduction by monic h, f <= 2^19, 72 Barrett (<= 1,024): < 2^20.9; times phi < 2^19.3; per bit < 2^21.31: f < 2^35.5, g < 3,060 x 2^17, < 2^36 per relation, 2^86.86 units at 2^65.15 relations (<= 2^88).

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (0.29, 0.75), exact v_2 and kappa, exact sampled norms). H1' retention (Euler-product model; band charge 0.009 bits). H-EA (work only): H1''s survivor law prices the ECM share ([12w]). H-FFT (work only): F(N) bounds the local SSA work past 2^26 limbs (KERNELS.md). H2 region (uniform pools, exact flags, certified counts; sized on S, judged on G + H; S+G work with the 1.01033 allowance). H3 matrix (r unramified, [13]; r prime to the class numbers, the unit condition, a one-dimensional kernel; N(400)/U' <= 0.0274, section 4). H4 kernel cost (KERNELS.md, STORAGE.md); H4u (unit share >= 86%; k5 <= 8.38, PREREG-U, [10]); rec7's per-hit price at production stays below its h' bucket's measured maximum; contaminated runs hold at most 0.091% of the canonical a (95% bound; a shortfall changes cost only); f/g independence given q_g (100x). H-X: BW's nondegeneracy polynomial is nonzero on shuffled block-diagonal X (evidence16.py --hx 1 16 512 25: 2^16 filtsim toy, N = 2063 > 4m, (m, n) = (512, 25), p = 2^31 - 19: both block Krylov bases rank N, so the projection does (Sylvester): rank X = rank Y = 2063, all 512 blocks non-empty; [10x] toys (40, 5), (100, 5) at 2^12; production blocks 2^48.3 coordinates). H5 descent (DESCENT.md): lower sums; per-ideal yield certificate (DISC); lazy order independent of outcomes; negative-binomial recurrence, independence across relations (boot psi^18) and of distinct-ideal subtrees (measured: [7e]); missing logs (Poisson occurrence at the 2-SE input); local omitted-mass laws; random-integer ECM orders; UDUP = 0.65 (DESCENT.md, [11]). H5c child load: window primes follow H1's law given semismoothness; P nonincreasing in the norm; smaller norms by the [7m] grid (DESCENT.md). H3 support is synthetic filtsim rows on the production occurrence profile at 2^12..2^20, not actual-relation filtering; the charged cap and P_CAP are widened here (section 4).

## 10. Memory

Core storage <= 2^105 bytes (STORAGE.md), wrapper below 2^30: below 2^106 bytes. All preprocessing is for this target, charged once in W.

## Provenance

GT core: kaankolcu recover-P 1b8fa720 lineage (9332a74d mpjunior92; 2c2ce2ab / 42e80f80 gbotrel; 6282ab79 bernardo-lima12; 6734bd5a GPT-6; 45ce22e5 / 4bb6c381 Gemini); polynomials after Guillevic-Singh; filtering after CADO-NFS. H-EA, H-X, rec7, band strip/split, f-only sieve, g batch, twisted C4 field, line-walk x8: gbotrel (Claude Opus 5.5). In-place record read: 7fd7a21c; FFT generator: 77c7dc97; wider BW and the generalised rank toy: 68b5e099; plain wrapper: 8e3fe531 lineage; child load: this package.


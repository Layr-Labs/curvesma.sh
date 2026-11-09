# recover-P (bls12-381) through two GT logarithms with the STNFS core: 2^115.609030000001 units

Track bls12-381-fixed-g2-invert-g1. Input: a nonidentity Q in G2 and z in GT; output P in G1 with e(P,Q) = z
for the exact supplied Q (P = infinity for z = 1). The score unit is one dense F_p12 multiplication, C = 20,000
scalar AArch64 instruction equivalents. No advice, smaller targets or amortization
are used; no prior verdict is a premise.

## Claim

The algorithm below returns such a P with success at least 0.39 over uniform
challenges and fresh coins under the declared heuristics H1-H5, H-EA, H-X and H-FFT (section 9). All classical work (preprocessing, both logarithms, failures,
restarts, field map, pairing, scalar multiplication, verification, output) is at most

    W <= 2^115.60903 + 2^66 < 2^115.609030000001  dense-Fp12 units (C = 20,000).

The score is 115.609030000001 (core: sections 2-8; wrapper: section 1). Analytic, exploratory claim.

Conventions. Bounds round in the safe direction; evidence16.py checks each against its unrounded value
(declared_bounds) and reruns the checks cited [n]. Model:
stats16.b85 (config, pool_d), the constants of model16.py and descent16.py (owner10.install: C_ROW
862).

## 1. Reduction (fixed Q: two GT logarithms) and field map

With the pinned reduced pairing e and the public generator P0 of G1:

1. Check Q (on the twist, [r]Q = infinity, Q != infinity) and z (z^r = 1). If z = 1 return P = infinity.
2. g = e(P0, Q) (pinned pairing(Q, P0), protected encoding).
3. Map g and z into the core's tower (field-map.md).
4. Run the core on (g, z) under its meter and node counter (any failure is a failure): a = log_g(z) mod r.
5. P = [a]P0; output P if e(P, Q) = z, otherwise fail.

GT has prime order r and e is nondegenerate: for the supplied Q, g generates GT, z = g^a for a unique a mod r, e([a]P0, Q) = z. For uniform Q, g is uniform on GT minus {1}, independent of the uniform z: the core's input law. The
field isomorphism preserves z = g^a; the final check rejects any wrong output.

One relation/LA phase serves log z and log g (section 6; 2/r charged). Wrapper below 2^80 instructions (pairing, scalar multiplication, check < 2^50, field map < 2^67, field-map.md):
2^65.71 units, a 2^66 row. Wrapper failure: field map below 2^-48.

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
(3/8) A^2; b uniform in the half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a)). Region: a5 > 0, Q_a(a) <= (2635/8000) A^2; b5
>= 1, Q_b(b) <= (17/20) A^2 - Q_a(a) (exact), so the restricted sample is uniform on it. G has 30123 region rows (20951
primitive) of 87197 canonical draws, H 60386 (41831) of 174159; G + H is scaled by region volume over the effective
sample 1968000 (2635/3000)^3 = 1333531.223. Weights: w = (1 - qa)^3 (yield), vb (1 - qa)^3 x factor (work, population);
slice-geometry costs use S' = vb (17/20 - qa)^3 x factor and K' = floor(sqrt((17/20 - qa) A^2 / Q_b(e_5))) lines. The
canonical a is the (Q_a, lex)-least element with a5 > 0 of its <+-1, units, sigma> class (exact enumeration of tO, t in {a, sigma a}); with b5 >= 1, one point per <+-1, units, tau> orbit (tau^2(a, b) = (a, -b)); non-primitive pairs are excluded by an exact HNF test. lattice_count_bounds.py certifies (exact rationals) the count ratio at every radius R >= 831 = isqrt((2635/8000) 2^21), below the smallest b-slice
radius 1044 ([3]): [0.977040, 1.022662] (Q_a), [0.989055, 1.010083] (Q_b); the model uses the products 0.966347073
(yield) and 1.032972622 (work).

Pools (PREREG-G). 656,000 uniform pairs each on theta = 3/8: S (seeds 261-264) sizes the region, G (251-254, drawn
after the configuration was fixed) and H (255-260, 265, 266: 1,312,000 pairs, pre-registered after G) are
judged, S and G give the work; no draw dropped. pool_d (stats16) holds per-seed counts (H: total), list hashes and sums, 8 anchor rows of G, G's canonical lists over draws 0..49 ([2] replays them) and the contamination sample. theta is the smallest t on the 1/8000 grid with margin_S
>= 0.06 + sqrt(2) SE_S (pre-registered), on S with the exact 2-adic yield: 2635/8000 (0.100072 >= 0.099549; 2634: 0.099258 < 0.099569; [4s]). The band
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
relations <= G + H's 2^62.66428 + 2 SE over the yield envelope
(1.032972622/0.966347073), the ECM miss (0.999^22: <= 13 + 9 primes above 2^64), the band margin and the least
retention on the model's norm grid (0.7775 f, 0.5393 g): <= 2^64.1165, 2^100.5922 <= 2^100.62 units (ledger row; audit).
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
lattice (t = r b(s) - a(s) mod p), walked once per (run, ideal) with one setup6 ([10r], KERNELS.md).
Clean-run rule: rec7 walks a run only if all its a-values are in the canonical list (tested in the 2^92
row); otherwise its slices take the slice walk. A clean run of 32 to 1,024 a-values uses the run walk in each band
p in [2^48, 2^64) with h' = run cells / p >= 2^4.99 (the table's bottom) instead of the slice walk and empty cosets; run
cells are exact (runs16.lrs: integer sum over the run's a-values of (17/20 - qa)^3). Each row pays its cells' share
(per hit at h', + 2 per record; setup6 per (run, ideal) replaces its slices', f ideals). Contamination. A contaminated run's slices cost at most
the slice walk: per slice <= D_SUP = contam_sup = 2^72.591 on [2^48, 2^64) (exact band maxima with empty cosets and setup6;
[k7c]: 2001 region slices peak at 2^72.446). Every row pays RUN_ALLOW x D_SUP/2 per slice (f ideals), RUN_ALLOW = 0.00091 (draws are
uniform over canonical a): in 4,000 a-uniform runs (PREREG-C: seeds 5001-5010, the first 400 canonical a with Q_a <=
(2635/8000) A^2 per seed, sha256-16 in pool_d) the exact predicate found no non-canonical a-value among 1,989,563
(166,465 new):
contaminated runs hold under 1 - 0.05^(1/4000) = 0.00074866 <= 0.00091 of the canonical a (95%). [k7c]/[2r] replay 6 sampled/3 anchor runs.
Coset sieve at P2. 2 = P2^2 P4^2, P2 = wO; b is in P2 iff P2F.b is even (common12.P2F; [10r] checks wO = ker P2F).
If a and b are in P2 the pair is not primitive, so those cells are not sieved; P2 is fixed by units, sigma, tau and
b -> -b, so the canonical domain is unchanged. In the basis e_0, e_j - (P2F.e_j) e_0 a slice with a in P2 sieves
b_0 = 2k + 1 (half the cells, congruence constants halved mod p, setup6 + 10). P2F(a + t) alternates along a run, so
rec7 walks it as two cosets (1/2 and 1/4 of its cells, setup6 + 10 each). model16.evaluate: p2a-and-p2b rows carry
nothing; p2a slices are priced at S/2 or the sub-run's cells.
Count envelopes ([3]) keep the per-a charge >= 0.000861 D_SUP >= 0.000749 D_SUP.
runs16's run ends (< 2^69.3 units) are misc.

Yield ([4], pools G + H). rels = (va_region/1333531.223) vb (0.966347073) sum over G + H's canonical primitive region rows of
w rho(u_f) rho(u_g) r_f 0.999^lambda_f x the f band margin (g exact): Dickman rho ([1]), retentions r from the
Euler-product threshold model (threshold_model.py, calibrated on exact smooth integers at B = 2^10..2^14; prime_min
2^14), the ECM miss; cells with qa < 0.02 are sieved and charged, relations discarded. Estimate
2^62.66428 (SE 0.0352 bits) less the band charge 0.009: 2^62.6552 against 1.05 U' <= 2^62.52173: margin >= 0.1335
bits = 3.8 pool SE (>= max(0.06, 2 SE)). Band charge (H1'). On 4,000 region rows the unsieved
[2^8, 2^14) f content is 5.511 +- 0.129 bits (model 5.922), P(> 40 bits) <= 0.0042 (95%; model 0.0021): loss <= 0.0032
<= REL_BAND_LOSS = 0.009 bits ([4]). Descent input. LOG2_RELS = 62.583 is <= the charged count
lowered by two relative pool SE, 2^62.58313 ([4] rels_input); the descent's mean norms 740.0/498.4 truncate G + H's
740.096/498.416 ([4]).
Exact 2-adic part (H1, [14v], [14t]). rho prices a random integer (mean v_2 = 1); each row's v_s = v_2 of its exact norm
is known (pool_row), so u_s = (n_s + alpha_s + E_s + kappa_s - v_s)/69 (model16.a2, also q_g): E = (2, 2.399985) >= the
local means (runs16; alpha_2 = 1 - E), kappa = (0.06, 0.05) >= kappa2(u) (sum_k 2^-(k+1) rho(u - (k-1-kappa)/69) = rho(u),
rising in u) at n_f 840, n_g 530, above every yielding row. Cells with a in P2 carry exactly 2^6 in N_f, the rest 2^0.

Cost ([5]). Recorded S+G sums still price CELL_ROUTE at 0.5. CELLRES from 2^24 to 2^26 cells slopes to
(862452-779083)/(63087346-15511016) = 0.001752 Ir/cell; [10s] requires the charge
to be at least 4x that slope; CELL_ROUTE = 0.02. evidence16 subtracts 0.48 per recorded cell
weight from the work first moment only. The stored second moment is then an upper bound, and the relative SE
it implies is 0.005049. Allowance 1.01015 >= 1 + 2 SE. Collection
<= 2^115.300026350 (estimate 2^115.285456811 times 1.01015), over the kept cells.

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U', filtsim.PEXC), clique removal and Markowitz merge (order: (rows - 2) x lightest unit pivot)
are simulated on the production occurrence profile (filtsim.py; [6] reruns 2^12) as the tau-quotient system: tau^k(I) enters I's orbit column with
coefficient chi^k (entries a + b chi, merged |a|, |b| <= 255). Cap. All retained relations are used. The region is the same
draw, restricted to Q_a <= (2545/8000) A^2 on seeds 261, 262 and 251 (draws 0..39999, canonical in-region rows): yield
0.959723 and adjusted work 0.929900 of the 2635/8000 pool, cells 0.935049. Those factors scale S, G, H and the previous
sizing step; the sizing excess is 1.007709 so [4s] stays the smallest 1/8000 step. H1 puts the 2-SE count at >= 2^62.5238
>= 1.040 U'. H3: six seeds at 2^19 and excess 1.04, means 744/501, give N(400)/U' mean 0.029862, mean + 2 SE 0.030528;
the attack aborts above 0.03055 U' (NU_H 0.03053; P_CAP = 1e-9). 2^20 was not rerun at this excess; the charge is the
measured 2^19 mean + 2 SE, and that mean is below the stored 2^15 mean ([10]). FILT_Y is unchanged. N <= 0.03055 U'
<= 2^57.4186 < 4 x 2^64/320: four owners (STORAGE.md).

The square system A = [M | S] has N + 10 rows, N orbit columns and 10 dense Schirokauer columns: the chi-eigenpart of
the unit groups, (rk K_f - rk K_f^+ + rk K_g - rk K_h)/2 = (35 - 17 + 5 - 3)/2 from signatures K_f (0, 36), K_f^+ (0, 18),
K_g (0, 6), K_h (2, 2) ([13]; tau^2 = -1 on the sigma-odd part, so its chi and chi^-1 eigenspaces have equal dimension);
under H3 its kernel is the line of virtual logs. Block Wiedemann (m, n) = (9216, 100): Y uniform in F_r; X block-diagonal
(H-X): shuffled coordinates cut into m blocks, X_l uniform in [0, 2^128) on block l. Krylov N/m + N/n block steps (20 five-vector groups in turn); Horner runs as five interleaved chains from checkpoints
A^(jL) Y (L = Horner length / 5; 4n = 400 checkpoint vectors, < 2^71.1 bytes, 4L < N/m + N/n) at the five-vector rate. Each
vector element is stored with chi x_j. Rows reduce to 2^-320 x value (BW on 2^-320 A). Per row and vector (STORAGE.md, [10]): 401 x 8.38 (H4u, PREREG-U) + 862 (row tail 857, four owners) + 10 x 97 + 200 (_chi_mul); per
Horner element n dense terms and one reduction (270); projections: each vector element meets one
X entry, (1 + n/m) N^2 terms at 60. Linear algebra <= 2^113.232777467 units; generator (matrix Berlekamp-Massey, recursive
doubling, k = m + n = 9316, L = (N + 10)(1/m + 1/n) + 2; model16.generator: 4 k^3 log2 L products, each Kronecker in
9-limb slots (entries < 2^257, L 2^514 < 2^568): fft_F(18L) + 1522L; slot s: REDC(s mod 2^512) + REDC((s >> 512)
(2^512 mod r)), inputs < 2^516; H-FFT as P_g (2^57.7 limbs; 18L < 2^55)): <= 2^102.577 units.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent; with seeded ECM curves the relation count is fixed:
2^62.6552 under H1 (section 3). Under H3 (section 9; the 10 Schirokauer maps of the chi-eigenpart
injective on units modulo r-th powers) Block Wiedemann fails with probability <= 4N/2^128 < 2^-66 (Schwartz-Zippel,
degree <= 2N for any (m, n), entries uniform on >= 2^128 values; H-X); a failed solve (checked on the dropped rows) is a failure.

Reconstruction. Purged ideals' logs come by back-substitution; ideals with <= 2 occurrences are log-less in the descent
([6b]: model >= 6.3x measured; DESCENT.md).

## 6. Individual logarithms

DESCENT.md ([7], [7e], [7l]): 62 nodes + boot (215, 304, 15, 4), retention >= 0.6490.

By H5, 1 - psi equals the product over window children of (1 - phi_child), minus the missing-log union term. Full child work is still charged. Expected two-log work <= 2^109.5279. The meter caps both logs at 2^110.3815:
P(cap) <= 2^(109.5279 - 110.3815) <= 0.553406 (Markov); at +10% expected work, still > 0.39.

## 7. Success

Core: success >= (1 - 2^-66 - P_CAP) x (1 - 2 x 5.10e-06 - 0.553406 - 2^-12.03 - 2^-10.70 - 2/r) >= 0.4456 (main
phase, two boots, Markov, node counter, repeats); >= 0.3904 with Markov at +10% (0.608747).

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX = 2^27.3. E[N] uses the same
  recurrence with unit node cost (counts16.py): E[N] <= 2^15.2673 for both logs, P(abort) <= 2^-12.03.
- Repeats: N_MAX nodes in the smallest bin (2^70.30 ideals), 2^6 for 1/p weighting: P <= 2^-10.70, charged.

With the wrapper: P(success) >= 0.4456 - 2^-48 > 0.39; every abort is a failure.

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (section 3; estimate 115.285456811, allowance 1.01015) | <= 115.300026350 |
| linear algebra (weight 401, N <= 2^57.4186, (9216, 100), block-diagonal X; section 4) | <= 113.017908341 |
| matrix generator | <= 102.577 |
| both individual logarithms (metered cap; expected 109.5279) | 110.3815 |
| canonical a list | 92 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 73 |
| wide factor-base lattices ([2^63, 2^64), A and coset, KERNELS.md) | 71.6 |
| rec7 bases (2 per ideal x run and two cosets, p in [2^48, 2^64), KERNELS.md) | 75.09 |
| band split (section 3, [10t]) | 100.62 |
| P_g product tree (section 3) | 74 |
| g and f band-2 split of kept relations (section 3) | 96.5 |
| numerical certificate construction | 50 |
| verification | 11 |
| wrapper (section 1; <= 2^80 instructions) | 66 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^92 (2^42.8 per region a, 2^63.100 of them); Schirokauer maps 2^88 (below); factor base 2^73 (2^23.4 per prime, pi(2^69) < 2^63.46); filtering and
back-substitution 2^80 (2^24.4 per entry against < 2^14). The rows sum to <= 2^115.6090274 (audit-checked).

Schirokauer maps. R = (Z/r^2)[y, x]/(h, f) (rank 72; g: rank 12); r is unramified ([13]): eps = lcm(r^d_i - 1) < r^72
(18,360 bits; g: 3,060). The 10 columns are fixed F_r-linear forms (chi-eigenpart) of (phi^eps - 1)/r mod r, phi = a - bx,
by left-to-right powering on eight-limb coefficients. An f squaring: 2,628 = 72 x 73/2 schoolbook 8 x 8-limb products
(<= 512 instructions), reduction by monic h, f <= 2^19, 72 Barrett (<= 1,024): < 2^20.9; times phi (12 coefficients
< 2^64) < 2^19.3. Per bit < 2^21.31: f < 2^35.5, g < 3,060 x 2^17, so < 2^36 per relation, 2^86.86 units at 2^65.15
relations, within 2^88.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (0.29, 0.75), exact v_2 and kappa, exact sampled norms). H1' retention (Euler-product model; band charge 0.009 bits). H-EA (work only): H1''s survivor law prices the ECM share ([12w]). H-FFT (work only): F(N) bounds the local SSA work past 2^26 limbs (KERNELS.md). H2 region (uniform pools, exact flags, certified counts; sized on S, judged on G + H; S+G work with the 1.015% allowance). H3 matrix (r unramified, [13]; r prime to the class numbers, the unit condition, a one-dimensional kernel; N(400)/U'
<= 0.03053, section 4). H4 kernel cost (KERNELS.md, STORAGE.md); H4u (unit share >= 86%;
k5 <= 8.38 at excess 1.04 (K5_RUNS, [10])); rec7's per-hit price at production stays below its h' bucket's measured maximum; contaminated runs hold at most 0.091%
of the canonical a (95% bound; a shortfall changes cost only); f/g independence given q_g (100x). H-X: BW's nondegeneracy polynomial is nonzero on shuffled block-diagonal X ([10x]: toys (9216, 100), (40, 5), (100, 5) reach full rank; blocks <= 5 coordinates, production 2^44.25). H5 descent (DESCENT.md): lower sums; per-ideal yield certificate (DISC); lazy order independent of outcomes; negative-binomial recurrence, independent
subtrees of distinct ideals (measured: [7e]); missing logs (Poisson occurrence at the 2-SE input); local omitted-mass laws; random-integer ECM orders; UDUP = 0.65 (DESCENT.md, [11]).

## 10. Memory

Core storage <= 2^105 bytes (STORAGE.md), wrapper below 2^30: below 2^106 bytes. All preprocessing is for this target,
charged once in W.

## Provenance

Parent: a4a4b1ba, with 77c7dc97's FFT generator; new (m, n).


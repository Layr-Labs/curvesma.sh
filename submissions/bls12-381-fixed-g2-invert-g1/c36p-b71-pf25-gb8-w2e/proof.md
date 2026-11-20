# recover-P (bls12-381) through two GT logarithms with the STNFS core behind a run gate: expected work 2^111.638187949370 units

Track bls12-381-fixed-g2-invert-g1. Input: a nonidentity Q in G2 and z in GT; output P in G1 with e(P,Q) = z
for the exact supplied Q (P = infinity for z = 1). The score unit is one dense F_p12 multiplication, C = 20,000
scalar AArch64 instruction equivalents. No advice, smaller targets or amortization
are used. Predecessor work is credited at the end; no prior verdict is a premise.

## Claim

W is the total classical work of the one attack on the one target. Step 0 is a run gate (kaankolcu's device, Provenance): draw a 64-bit word u;
unless u < K = 7202073246512176891, stop as a failure (< 1 unit). K is the least integer with K x 0.998911 >= 0.39 x 2^64. Under H1-H5, H-EA, H-X and
H-FFT (section 9) the algorithm below returns such a P with probability >= 0.998911 once past the gate (section 7); u is independent of the challenge and the
other coins, so success is at least 0.39 over uniform challenges and fresh coins, and

    E[W] <= 1 + (K/2^64)(2^112.9950699693 + 2^66) < 2^111.638187949370  dense-Fp12 units (C = 20,000).

2^112.9950699693: section 8's rows: every-run bounds (collection and linear algebra are the same in every run; their measured means carry + 2 SE) except two at expected work, the two logarithms, 2^108.841278 (section 6), and the ECM split of kept relations (section 3). The score 111.638187949370 is an
expectation, not a bound on every run: all work of a run (preprocessing, both logarithms, failures, restarts, field map, pairing, scalar multiplication,
verification, output) is at most

    W <= 2^121.9028378233 + 2^66 < 2^121.902837823301  units.

Analytic, exploratory claim.

Conventions. Bounds round in the safe direction; evidence16.py checks each against its unrounded value
(declared_bounds) and reruns the checks cited [n]. Model:
stats16.b85 (config, pool_d), the constants of model16.py and descent16.py (owner10.install: C_ROW
572, 4 owners).

## 1. Reduction (fixed Q: two GT logarithms) and field map

With the pinned reduced pairing e and the public generator P0 of G1:

0. Run gate (Claim): draw u; if u >= K, fail.
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

Candidates: (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^71; factor base: relative-degree-1 ideals of
norm <= B, U' = 2^64.409208 tau-orbits (U/2: orbits of 4 ideals; U = 2^65.41).

Region. Exact rational forms Q_a, Q_b (stats16 config), A^2 = 2^21. Pool (common12.regenerate_pair): uniform on theta = 3/8 (a5 > 0, Q_a(a) <=
(3/8) A^2; b uniform in the half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a)). Region, by 2-adic class of a. Cells with a in P2 (class Y, flag p2a; X: the rest) carry exactly 2^6 in N_f whatever b is (below; all 144,173 primitive pool rows): a Y slice is richer per cell by the rows' own yield law (model16.a2). Class scales s = (47/50, 103/100): a5 > 0, Q_a(a) <= s_c (2493/8000) A^2; b5
>= 1, Q_b(b) <= lam_c A^2 - kap Q_a(a), lam_c = s_c 7/12, kap = 7/20 (exact), so the restricted sample is uniform on it. G has 13302 region rows of 87197 canonical draws, H 26389 of 174159, K 53393; sums are scaled by pool points over all draws. Weights: w = (1 - qa)^3 (yield), vb (1 - qa)^3 x factor (work, population);
slice-geometry costs use S' = vb (lam_c - kap qa)^3 x the lower factor and K' = floor(sqrt((lam_c - kap qa) A^2 / Q_b(e_5))) lines. The
canonical a is the (Q_a, lex)-least element with a5 > 0 of its <+-1, units, sigma> class (exact enumeration of tO, t in {a, sigma a}); with b5 >= 1, one point per <+-1, units, tau> orbit (tau^2(a, b) = (a, -b)); non-primitive pairs are excluded by an exact HNF test. Counts ([3]). Points per volume of each half-ellipsoid, measured by the sampler's acceptance experiment (common12.uniform_half_ellipsoid frame box, exact hit test; 2.5e11 trials per box): Q_a (pool box, every t/8000, t = 2300..3000): 0.999206..0.999331, bracket [0.999062, 0.999429]. Q_b (four boxes, 32,840 radii^2 in [0.44 A^2, A^2], step 2.5e-5; count nondecreasing, volume ~ radius^6): [0.998999, 0.999595] at every radius^2 in [0.44 A^2, A^2] (all pool slices (1 - qa) A^2 and region slices (lam_c - kap qa) A^2: s_c (7/12 - 7/20 x 2493/8000) = 0.4458 >= 0.44, s_c 7/12 <= 1, [3]; a: the pool ball). Each point carries 6 SE (SE = ratio/sqrt(hits) <= 2.4e-5; union bound over 67,082 one-sided events < 1.03e-3), rounded outward. Factors: 0.998061 = 0.999062 x 0.998999 (yield, run cells, S') and 0.999595 = 1 x 0.999595 (work, population, D_SUP; a's end raised to 1, so each product also bounds b alone). Fallback: lattice_count_bounds (lib16.b85) certifies (exact rationals) the bracket (0.963535, 1.035815) at every radius >= 766; [3] checks the brackets inside it and replays 2^17 trials of two boxes.

Pools (PREREG-G; rows regenerated at this region). 656,000 uniform pairs each on theta = 3/8: S (seeds 261-264) sizes the region; G (251-254, drawn
after the configuration was fixed), H (255-260, 265, 266: 1,312,000 pairs, pre-registered after G) and K (2801-2816: 2,624,000 pairs, drawn after c27 was filed) are
judged, every seed kept; S, G and K give the work; the region keeps a only if its run has T >= 48 (model16.RUN_TMIN); no draw dropped. pool_d (stats16) holds per-seed counts (H, K: totals), list hashes and sums, 8 anchor rows of G, G's canonical lists over draws 0..49 ([2] replays them) and the contamination sample. theta is the smallest t on the 1/8000 grid with margin_S
>= 0.06 + sqrt(2) SE_S (pre-registered), on S with the exact 2-adic yield: 2493/8000 (0.112769 >= 0.112289; 2492: 0.111954 < 0.112317; [4s]). The band
charge applies to the judged count only.

Sieve and survivors (f ideals only), two stages. Stage A. Coset replay (scores only) for p in [2^24, 2^48); record walk for p in [2^48, 2^64) (slice12,
KERNELS.md; X2 = 2^64; per clean a-run, below; [2^63, 2^64) of contaminated runs on the wide path). Survivors pass the f threshold at slack SA = 356.0; the norm kernel (C_NORMS = 48,300 measured with L and the compare; KERNELS.md, [5t]) gives N_f and the g norm z. Recheck. The scan hands each survivor's score s on; the kernel returns L (L <= 16 log2 N_f < L + 2) and the survivor goes on only if L - s <= 16 SA. Each of the <= c scored ideals (c = floor((n_f - v_2)/24)) is rounded down by < 1/16 bit, so L - s is within (-2, c) of 16 x the omitted mass (N_f outside the scored ideals): every cell with omitted mass <= SA - (c + 1)/16 goes on, and one that goes on has omitted mass < SA + 1/8.
Stage B (rechecked survivors). y = P24 mod N_f by the scaled remainder tree (lemma below; P24 = the product of the primes below 2^24; waves of 2^15 norms); with N_f = 2^v z', m = 2^v gcd(z', y^(2^10) mod z') is the part of N_f made of primes below 2^24 (each such prime of z' divides y, and v_p(z') <= 832 < 2^10; model16.k2, [10v]). With Lm = floor(16 log2 m), O = L - s - Lm is within (-2, c + 1) of 16 x the outside mass (N_f outside m and the scored ideals). The survivor goes on iff O <= 16 SB, SB = 280.5, and 16 x 64 k - 1 <= O <= 16 x 71 k + c + 1 for some k in 0..4. Window lemma: if N_f is B-smooth and its primes in [2^24, 2^64) divide it only through scored ideals, the outside part is k primes of (2^64, B], of mass in (64 k, 71 k]: such a cell with outside mass <= SB - (c + 1)/16 goes on (k <= 4: 5 x 64 > SB); a cell that goes on has outside mass < SB + 1/8, within (c + 3)/16 of a window. Cells where a prime >= 2^24 divides N_f beyond its scored ideals (a square, an ideal of degree >= 2) are outside the yield law (squarefree sieved part). N_f carries 2^0 (class X) or 2^6 (class Y) exactly: v_2(N_f) = v_w(F) + 2 v_pi4(F), F(a, b) = b^12 f(a/b); at (w) a^12 (w not | a) or w^6 b^12 (w | a) is the unique minimal term (0 or 6); mod pi4 (F_4, w^3 = 1) F/b^12 = t^12 + w^5 t^2 + w^6 is 1 or w^2 t^2, never 0, and pi4 | b forces a unit a: on Y slices T and the recheck use SA + 6 (stage B needs no offset: m holds the 2^6), so every test is on the odd part.
Laws (threshold_model.two: an Euler product without the prime 2, read at n = n_f - v_2; primes below 2^24 small with all powers, [2^24, 2^64) sieved, above large; explicit below 2^14, prime density above). Yield: r_f = P(sieved and large parts squarefree, small + large <= SA - (c + 1)/16, large <= SB - (c + 1)/16 - 1/8 | smooth), 0.7587 at n = 740 ([2w]). Work (H-EA; all integers, all powers): norms on survivor_upper(n, SA, 32), stage B on survivor_upper(n, SA, 1/4), the g batch and later steps on the same law with large <= SB + 1/4 inside the windows widened by (c + 2)/16 + 1/4. A survivor with v_2(N_f) > 6 stops after the norm kernel (a primitive pair carries 2^0 or 2^6 exactly, above; such a cell is never a relation), so stage B and later steps are charged only on rows with v_2 <= 6 (model16.evaluate). [12w]: on 4,000,000 random integers with exact factorisations at two toy sizes (law with the prime 2; stored counts, first 200,000 of each replayed) both stages' survivors are below the laws, the kept smooth integers above r_f. [4b]: r_f recomputed with the exact odd parts below 2^24 of 8,000 pool norms in place of the law's small part (stored 2-bit histogram, upper bin edges; H1': the small part is independent of the rest) is 0.8161 / 0.7621 / 0.7043 at n = 700 / 740 / 780 against 0.8133 / 0.7587 / 0.7003; at - 2 SE the loss is <= 0.0006 <= REL_BAND_LOSS = 0.002 bits. C_SB = 820,000 per rechecked survivor (KERNELS.md, [10k]). Survivors with B-smooth z get trial division by the 82,025 primes below 2^20 on both cofactors (2^23.649 <= C_TDIV = 2^24.5, [5t]). Batch band strip ([10b], KERNELS.md). Survivor sides wait in
waves of K = 2^39 (<= 832 bits). P = product of the primes in (2^20, 2^48] (< 1.000028 x 2^48/ln 2 bits, Rosser-Schoenfeld); a scaled remainder tree
(Bernstein; model16.srt) gives P mod z, the band part is gcd(z, (P mod z)^(2^6) mod z) (2^6 >= 41 >= any band exponent);
the f band-miss factor (1 - 41 x 0.001^2) stays as a margin. Records carry scores only; a g-passed f side gets a second strip, P2 = the primes in (2^48, 2^64] (exponents <= 17 < 2^5; K = 2^54, C_BATCH2 = 2^25.25 >= 2^25.2401).
g batch ([10g]). Same SRT with P_g the product of all primes <= B (< 1.000028 B/ln 2 bits), waves of K = 2^63 g norms (<= 512b, 8 limbs: ng_bound_512 certifies |N_g|<=2^511.52<2^512 via Rouché, Fraction LDL^T & AM-GM): z is B-smooth iff (P_g mod z)^(2^10) = 0 mod z (exponents <= 512 < 2^10). r_g >= 1 - 2^-38 (ideals of norm p^d > B drop a relation; in the margin). A kept relation's g norm is split by trial division, band split and ECM (<= 12 primes of (2^48, B]): <= 12 x 1.002 x 2^40.20 + 13 C_MR = 2^43.78 per relation, <= 2^65.33 relations ([10t]). Its f band-2 part (<= 17 primes) is split by ECM until found (1.002 x 2^41.06 at 840 bits + 2^26
lookup each, + 18 C_MR(840)): with g 2^96.66 <= 2^97.0 units (ledger row): an expected work (ECM model: independent schedules, each finding its prime with probability >= 1/1.002), metered at the row (reached = failure): Chebyshev (count variance <= 0.00201, a schedule <= 2^26.8 units) gives P <= 0.00201 x 2^26.8/(0.25^2 x 2^97.0) < 2^-75, within the rounding of 0.998911.
Lemma (exact SRT). Product tree Z_v (L <= 63 levels); node v holds Y_v of b_v = 64 (n_v + G) bits (n_v limbs of
Z_v, G = L + 2). Root: Newton quotient q of P 2^b by Z, corrected by <= 2 add/subtracts, Y_root = q mod 2^b. Child c,
sibling s: Y_c = floor((Y_v Z_s mod 2^(b_v))/2^(b_v - b_c)). Leaf: r = ceil(z Y 2^-b) (0 if r = z) = P mod z. Proof: eta_v = (frac(P/Z_v) - Y_v 2^(-b_v) mod 1) Z_v gains < 2^(-64G) per level, so eta < (L + 1) 2^(-64G) < 1 (no wrap). Cost (H-FFT; products <= model16.fft_F, gbatch_cost): C_BATCH = 2^24 (2^23.921, K = 2^39), C_GBATCH = 2^24.92 (2^24.916, leaf = 8 limbs, K = 2^63). [10v]: model16.srt on 2^8
planted leaves (Newton offsets -2/+2) vs Python's %, K2 on them. P's, P2's and P24's builds (< 2^80 instructions) are misc, P_g's tree (< 2^88) the 2^74 row.
Band split ([10t]). A kept relation is split: the band parts are trial-divided by every prime of (2^20, 2^48] (< 2^43.272; 180 instructions per prime with a 2^(64 i) mod p table, misc row), each hit divided out to full power (<= 41 + 28 divisions, 2^10 each). Each found prime gets its ideal and exponent from the factor base (<= 2^26): 2^50.7634 instructions per relation, none lost. Kept
relations <= the judged 2^64.4848 + 2 SE over the count envelope
(0.999595/0.998061), the ECM miss (0.999^22: <= 13 + 9 primes above 2^64), the band margin and the least
f retention on the model's norm grid (0.5919): <= 2^65.3210, 2^101.7967 <= 2^101.80 units (ledger row; audit).
A stage-B survivor's f cofactor above 2^64 is < 2^288.5 (< 2^320, 5 limbs); ECM is the ecm16 (lib16.b85) schedule at 99.9% per
prime ([12]): 71-bit primes in a 288.5-bit cofactor cost 2^38.4954 instructions.
A smooth side's cofactor C (no prime below 2^64 after division and both strips) is 1, a prime <= B or composite; ECM runs only then. Compositeness: 40 Miller-Rabin rounds, random bases (C_MR = 40 x 2 x 288.5 five-limb
products = 2^23.017); a composite passes with probability <= 4^-40. Every stage-B survivor with smooth z is charged the whole schedule. g first: the f side's division, strip, test and ECM run only if the g batch passes:
per threshold survivor C_NORMS, per rechecked one C_SB, per stage-B survivor C_GBATCH + q_g (C_TDIV + C_BATCH + C_BATCH2 + C_MR + 2^38.4954), q_g = min(1, 2 rho(u_g)) per pool row
(rho(u_g) is the yield law's g factor, H1; the 2 covers P(z smooth | the cell passed stage B) / rho(u_g). The two norms of a coprime pair depend on each other through the primes dividing both; for the primes below 2^24 [12q] measures it on 8,000 real pairs (exact small parts, correlation -0.007 +- 0.011): tilted as the two events tilt them, and at twice those tilts, joint over product stays below the factor). S, G and K work sums use it on
all 3,936,000 draws (pool_row, runs16).

Thresholds, two levels (kaankolcu's idea, Provenance; model16.two_level). T = lower - 16 x slack, lower a certified lower bound of 16 log2 N_f (embedding_bounds.py: Taylor shift per place, outward rounding; NormBounds.block) for exactly the rows and cells it is used on: every B-smooth cell with omitted mass below slack is kept whatever the bound's width; a wide bound only sends more cells to the exact norms (charged at slack + 32 bits). Level 1: the rows of an aligned 8 x 8 x 8 box (model16.box_rows) share the hull (EB.hull) of the hulls of its 2 x 2 x 4 sub-boxes; in an aligned 256-cell block each 64-cell node meeting the box's cells is evaluated on it; certified (finite, <= 32 bits wide) it is a leaf for all rows, else halved, down to 16 cells. Level 2: a failed 16-cell node runs EB.box_leaves per sub-box (hull, bisection to cells, then rows; uncertified cells are rejected: none in the pools). Merge: with no level 2 in the block and (largest upper - smallest lower) <= 32 bits, the block is one leaf at its smallest lower bound (a lower bound for each cell). A row's part of a leaf is one scan leaf at 256 cells, else cut at multiples of 64 (measured scan lengths). Counts are exact per pool row: 1.6e-04 evaluations, 0.0106 scan leaves per cell. Prices: C_NB_EVAL 1,178, C_NB_ROW 378 (measured, KERNELS.md, [10s]). NormBounds.row: per slice and place c_j = mid e_j, w_j = up((rad e_j + |c_j| 2^-49) UP);
a row takes m = sum b_j c_j, r = (sum |b_j| w_j)(1 + 2^-48): nine roundings < 2^-49.8 relative, so |b(v) - m| <= r. [2b] checks exact f norms of six anchor blocks and two level-2 blocks against their leaves, [2] and the audit the anchors' fields.

Coset replay (p in [2^24, 2^48); KERNELS.md, STORAGE.md, [11r]). An f ideal (p, s, r) divides a - b x iff b(s) = kappa mod p, kappa = a(s) r^-1; b -> b(s) mod p does not depend on a. Every slice is the half-ellipsoid b5 >= 1, key(b) = db Q_b(b) <= db R^2 (db = 5 x 10^11, R^2 = lam_c A^2 - kap Q_a(a); keys < 2^59.6, exact in u64), inside the largest slice (class Y, Q_a = 0: V = 2^64.4506 cells, K = 1467 lines); a class Y slice (a in P2) keeps the cells with P2F.b odd, a class X slice all. So the hits of (p, s, r) in a slice are the cells of the largest slice with b(s) = kappa, of the slice's class, with key <= db R^2; of a only kappa and R enter. The roots r of f at (p, s) come in pairs +-r (f(x) = w^6 F(x^2/w) is even); for p >= 2^24 none is 0 (f_s(0) = w(s)^6, N(w) = 2) or repeated (Disc F = 3^3 x 19 x 67, Disc h = 2^13 x 19^2), so no ideal of the band needs a special case and one table per (p, s) serves all its ideals: at most half of the 2^42.99 f ideals of the band, 2^42.00 tables. A table lists the cells of the largest slice (pointer to the cell in its line's box, key) in buckets (kappa, line, class bit), sorted by key inside a bucket (cells visited in key order, one radix sort; two stable counting passes per table). Every slice uses the line box of the largest slice (its radices L4 and W0 and origin: cell = i4 W0 + r0 + b0 as in the line walk, from the largest radius), so a cell's pointer is the same in every slice and list entries are copied without translation; a slice's own layout visits only its own tails, priced as before. The bound of a slice is the integer floor(db (lam_c A^2 - (7/20) Q_a(a))), computed once per slice, and the test is key <= bound (ties included). Replay of (ideal, slice): kappa from the stored e_j(s) and r^-1; per line of the slice the bucket (class X: both class buckets) is cut after its last entry of key <= db R^2 by a branchless lower bound, the pointers are copied to the line's owner and each cell gets 16 log2 p by staddh. The cells hit and the sums are those of the line walk, so thresholds, survivors, yield and the pool rows are unchanged. Charges (instructions): 2.75 per replayed entry; per line 48 + 8 (log2(mu + 1) + 1), twice for class X, mu = V/(2 K p) the mean bucket size (the search takes at most log2(n + 1) + 1 steps; Jensen); 36 per (ideal, slice); lists 34 per cell and table, 2.25 per bucket, 147 per cell once: 2^97.26 units (ledger row), 2^110.47 bytes (STORAGE.md). Exactness: on the real forms and ideals at 2^15.69, 2^19.82, 2^22.88 and 2^24.92 cells (1,516,922 (ideal, slice) pairs, f and g ideals, both classes, 3.32 x 10^9 hits) the cells selected by kappa, class and key equal the cells hit by direct evaluation of a(s) - r b(s) mod p, 0 mismatches, controls (wrong kappa or class, no key test, strict compare) fail; 910 pairs at the production radii, p in [2^50, 2^62): 0 mismatches. The line walk is the reference of these tests; no ledger row charges it below 2^48. The run-walk campaign below drew its walks on [2^46, 2^64); the bands from 2^48 read its per-density prices (REC7_K).

Record walk per canonical a-run (rec7; KERNELS.md). A canonical a's run is the maximal integer t-interval
containing 0 with a + t e_0 canonical and inside the limit s_c theta A^2 of a's class c; rec7 walks its sub-run of class c (runs16.py). On a run, a prime ideal (p, s, r)'s hits in all its slices form one coset of a 7-dim
lattice (t = r b(s) - a(s) mod p), walked once per run and ideal pair (p, s, +-r) with one setup6 ([10r], KERNELS.md).
Clean-run rule: rec7 walks a run only if all its a-values are in the canonical list (tested in the 2^92
row); otherwise its slices take the slice walk. A clean run of 32 to 2,048 a-values uses the run walk (x-walk-asm: our hand-written AArch64 run walk, source not shipped; REC7_K x the reference curve, by its 512 x 10 campaign at this theta, 0 mismatches vs the reference walk; KERNELS.md) in each band
p in [2^48, 2^64) with h' = 2 x run cells x lower/upper ([3]) / p >= 2^4.99 (the curve's bottom) instead of the slice walk and empty cosets; run
cells are exact (runs16.lrs: integer sum over the run's a-values of (lam_c - kap qa)^3, halved). Each row pays its cells' share
(per hit at h'; setup6 per (run, pair) replaces its slices', f ideals). Contamination. A contaminated run's slices cost at most
the slice walk: per slice <= D_SUP = contam_sup = 2^71.761 on [2^48, 2^64) (exact band maxima with empty cosets and setup6;
[k7c]: slices peak at 2^72.33). Every row pays RUN_ALLOW x D_SUP/2 per slice (f ideals), RUN_ALLOW = 8e-05 (draws are
uniform over canonical a): in 40,000 a-uniform runs (first 40,000 T>=48 runs across G+H+K at theta_region = 2493/8000, sha256-16 4577efc722575daf in pool_d) the exact predicate found no non-canonical a-value among 19,158,169 points:
contaminated runs hold under 1 - 0.05^(1/40000) = 0.0000749 <= 8e-05 of the canonical a (95%). [k7c]/[2r] replay 6 sampled/3 anchor runs.
Coset sieve at P2. 2 = P2^2 P4^2, P2 = wO; b is in P2 iff P2F.b is even (common12.P2F; [10r] checks wO = ker P2F).
If a and b are in P2 the pair is not primitive, so those cells are not sieved; P2 is fixed by units, sigma, tau and
b -> -b, so the canonical domain is unchanged. In the basis e_0, e_j - (P2F.e_j) e_0 a slice with a in P2 sieves
b_0 = 2k + 1 (half the cells, congruence constants halved mod p, setup6 + 10). P2F(a + t) alternates along a run, so
rec7 walks it as two cosets (1/2 and 1/4 of its cells, setup6 + 10 each). model16.evaluate: p2a-and-p2b rows carry
nothing; p2a slices are priced at S/2 or the sub-run's cells.
S' uses the lower count factor ([3]), so the per-a charge is >= 8e-05 D_SUP >= 0.0000749 D_SUP.
runs16's run ends (< 2^69.3 units) are misc.

Yield ([4], pools G + H + K). rels = (va_pool/4592000) vb (0.998061) sum over their canonical primitive region rows of
w rho(u_f) rho(u_g) r_f 0.999^lambda_f x the f band margin (g exact): Dickman rho ([1]), retentions r from the
two-stage law above (threshold_model.py), the ECM miss; cells with qa < 0.02 are sieved and charged, relations discarded. Estimate
2^64.48474 (SE 0.0227 bits) less the band charge 0.002: 2^64.4827 against 0.96 U' <= 2^64.35031: margin >= 0.1324
bits = 5.8 pool SE (>= max(0.06, 2 SE)). Band charge (H1'): [4b], above. Descent input. LOG2_RELS = 64.436 is <= the charged count
lowered by two relative pool SE, 2^64.43658 ([4] rels_input); the descent's mean norms 736.1/497.4 are not above the judged
736.165/497.469 and its base is the collection's 2^71 ([4]; smaller means and count only lower kappa and lambda, DESCENT.md).
Exact 2-adic part (H1, [14v], [14t]). rho prices a random integer (mean v_2 = 1); each row's v_s = v_2 of its exact norm
is known (pool_row), so u_s = (n_s + alpha_s + E_s + kappa_s - v_s)/71 (model16.a2, also q_g): E = (2, 2.399985) >= the
local means (runs16; alpha_2 = 1 - E), kappa = (0.06, 0.05) >= kappa2(u) (sum_k 2^-(k+1) rho(u - (k-1-kappa)/71) = rho(u),
rising in u) at n_f 840, n_g 530, above every yielding row. Toy census (92,108 primitive (a, b) in {-1,0,1}^12, exact 2^20-smooth N_f): class chi^2 at 2^18/20/22 10/6/18 (exact) vs 88/30/23 (average); totals within 1 SE.

Cost ([5]). Sampled mean 14.39 instructions/kept cell over pools S, G and K, charged 14.52 with the allowance, over
2^122.994 kept cells (every step at KERNELS.md/STORAGE.md constants; rec7 at its measured stratum means + 2 SE): collection <= 2^112.565809621
units (estimate 2^112.552654692 times allowance 1.00916). Allowance (H2, [4a]): a Monte Carlo estimate over S, G and K's
3,936,000 draws, relative SE row 0.2954%, incidence sqrt((1-79983/3936000)/79983) = 0.3500%, total <= 0.4580%;
charged 1.00916 >= 1 + 2 SE (the H2 margin of the pool estimate).

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U', filtsim.PEXC), clique removal and Markowitz merge (order: (rows - 2) x lightest unit pivot)
are simulated on the production occurrence profile (filtsim.py; [6] reruns 2^12) as the tau-quotient system: tau^k(I) enters I's orbit column with
coefficient chi^k (entries a + b chi, merged |a|, |b| <= 255). All retained relations are used; H1 puts the count at >= 2^64.43658
>= 1 U' ([4]): filtering at excess 1, fixed by the rule margin >= max(0.06, 2 SE) on the judged count (0.005 grid; sizing and judgment at 0.96; margin
0.0735, [4]). Collection's ECM curves are fixed by a public seed, so
relations, filter and N are deterministic in the challenge-independent instance, and H3 is stated for it:
production N(250)/U' <= 0.0085 (supported by the campaigns below, <= 0.00835, and the [10c] spread law), N_l/U' <= 0.00377 and N_sw/U' <= 0.00835 (product form, below). H3 (section 9) bounds the production means; a fixed instance above a cap aborts (LA is
charged on the product form; N(250) anchors the kernel statement); P_CAP = 0.00022 is slack for the spread law; a failed cap aborts. [10c]: Chebyshev with sd(x) <= 0.046 U'^(-1/4) (2 x the largest per-run sd x U^(1/4); sd x sqrt(U) is 0.26-0.67: no U^(-1/2) law). A campaign at B = 2^71, excess 1, merge weight 250 and relation
means 740.3/498.7 (not below the judged pools' 736.165/497.469, [10y]) measures x = N(250)/U' on 8 seeds per scale: 2^16 0.01224, 2^17 0.01163, 2^18 0.01068, 2^19 0.00986, 2^20 0.00949, 2^21 0.00912, 2^22 0.00851, 2^23 0.00811 ([10] lists SE): it falls at every step. Correction: the explicit band takes every fourth prime per side ([::4], occurrence min(1, 4 kappa/q)), the pi(B)/2 orbit model of the other bands. The tower's actual tau-orbits (all degrees, norm <= 2^19.5) as band scale x by 1.0062 +- 0.0048 (2^18-2^22): NU_H carries 1.0159.
Columns (filtsim): explicit primes to lsmall = min(log2(Rp/25)+1, 0.8L), w_mid on [lsmall/L, 0.8], dense above: disjoint, w_mid >= 0. Measured scales (lsmall <= 19.48) are unchanged; production has w_mid = 0. w_mid is unpivotable row weight that grows on merges; measured (2^14-18), x falls as w_mid -> 0 or lsmall - 2: production is at or below the trend.
H3 (premise): x does not rise with scale, so production is not above 2^23: mean + 2 SE at 2^23 is 0.008215, x 1.0159 (splitting) <= 0.00835. FILT_RUNS (excess 1.05, weight 400) is the older campaign. This campaign is run directly at B = 2^71 (same profile as the product-form campaign); the N(250) cap anchors the kernel statement and LA is charged on N_l and N_sw. These reference campaigns merge to mean weight 250 (stored runs' mean weights in [250, 251], [10y]). Product form (recorded pivots, cf. Boudot, ePrint 2017/758; filtsim.pivots): the merge stops at mean weight 200 (N_sw rows); columns in (occurrence, index) order, explicit band excluded, take the lightest row with a unit entry there whose recorded columns all carry labels below every recorded row holding the column, and get the midpoint label in (0, 1). Recorded row p of column c defines x_c = -conj(p_c)(p - p_c x_c) and reads smaller labels only: a unit-triangular system (depth 14-25 at every scale). A product runs the recorded rows in label order into their vector slots, then the N_l live rows, each priced as a row (STORAGE.md). Rank and kernel are kept, so H3 and H-X concern the live system; S enters by one recorded pass per S column, the recorded logs by one more. Campaign (profile above at B = 2^71, exact Bernoulli occurrence distribution—per-row at 2^16..2^19 and exact geometric-skip at 2^20..2^25, 8 preregistered seeds 1..8 per scale, weight 200): N_l/U' 0.00639 at 2^16 to 0.00361 at 2^25, N_sw/U' 0.01327 to 0.00809 at 2^25; both fall at every step ([10p]; [6p] reruns 2^14). Selection reads each candidate row's labels once per column (<= sum of squared row weights < 2^73 reads; filtering row); labels are binary64 midpoints, a pivot only if lo < mid < hi. H3 covers them: mean + 2 SE at 2^25 x 1.0159 (the N(250) splitting margin) <= 0.0037, 0.00828; caps 0.00377, 0.00835 (sd <= 0.046/2 U'^(-1/4), [10c]). N_l <= 2^56.3580, N_sw <= 2^57.5052 < 4 x 2^64/320 (four owners, STORAGE.md).

The live square system A = [M | S] has N + 10 rows, N = N_l live orbit columns and 10 dense Schirokauer columns: the chi-eigenpart of
the unit groups, (rk K_f - rk K_f^+ + rk K_g - rk K_h)/2 = (35 - 17 + 5 - 3)/2 from signatures K_f (0, 36), K_f^+ (0, 18),
K_g (0, 6), K_h (2, 2) ([13]; tau^2 = -1 on the sigma-odd part, so its chi and chi^-1 eigenspaces have equal dimension);
under H3 its kernel is the line of virtual logs (w, c), c != 0, and M has rank N. Block Wiedemann (m, n) = (8192, 130) with S as right-hand sides:
B = [M | 0]; Y = (S_1..S_10, B z_11..B z_n), z uniform in F_r; X block-diagonal (H-X): the 10 pad coordinates (zero columns of B)
in blocks 1..10, the others shuffled and cut into m blocks, X_l uniform in [0, 2^128) on block l. A generator sum_k B^k Y f_k = 0 of
a_i = X^T B^i Y with f_0 = (c, e), c != 0, gives S c + B w = 0, w = Z e + sum_(k >= 1) B^(k-1) Y f_k: one Horner pass, no product touches S. d = ceil(D/n) + 1 + ceil(D/128n), e likewise in m, D = N + 10: d + e < (1 + 2^-7) D (1/m + 1/n) + 6 products and terms, charged; Horner within N/m + N/n (2^-7 m < n); Horner runs as five interleaved chains from checkpoints
B^(jL) Y (L = Horner length / 5; 4n = 520 checkpoint vectors, < 2^72 bytes, 4L < N/m + N/n) at the five-vector rate. Each
vector element is stored with chi x_j. Rows reduce exactly (_row_red). Per live row and vector (STORAGE.md, [10]): N_sw/N_l = 2.2149 evaluated rows x (201 x 8.29 x 0.8126 (H4u, PREREG-U; the evaluated rows' measured weight, k5 and grouped/plain ratio are in [10p], below the charges; tables below) + 572 (row tail, four owners) + 133 (_chi_mul)) + 18.4 x 2^47.71/N (table builds); per
Horner element n generator terms and one reduction (270); projections: each vector element meets one
X entry, (1 + n/m) N^2 terms at 60. Linear algebra <= 2^110.677012532 units; generator (matrix Berlekamp-Massey by half-gcd, k = m + n = 8322: a node of l
sequence terms makes (m + k) k^2 products of polynomials with <= l and <= l/2 + 1 coefficients, each one Kronecker integer product in 9-limb
slots (coefficients < k L r^2 < 2^576): F(18 l + 18) under H-FFT, + 640 l to pack, reduce (270) and add; leaves m k^2 (97 + 270). The sum over
the recursion is below model16.generator ([10m]): 4 k^3 log2 L (F(18 L) + 640 L) instructions, L = (1 + 2^-7)(N + 10)(1/m + 1/n) + 6): <= 2^100.549594 units.
Dense-band tables (subset sums). Per product and 5-vector group the group's 4th owner (owners 1-3 hold the other columns; it holds the remaining < 0.55 x (2^64//320) slots including all columns of norm <= 2^28) reduces x, chi x, 2r - x, chi(2r - x) into [0, r) for each orbit column of norm <= 2^28 and stores every multiset sum of g of them over norms <= 2^28, 2^18.25, 2^14 (g = 2, 3, 4), reduced into [0, r): 4^g/2 C(n_g + g - 1, g) + 4 n_2 = 2^47.71 elements at C_TAB = 18.4 (tab.txt), 2^55.04 bytes per group, under half its free space for all 26 groups ([10d]). A row's terms there are grouped greedily, densest band first; a group is one unit entry (same sums). Price (H4u): the FILT_X runs' merged rows grouped row by row (implicit mid weight by filtsim's law, Poisson), grouped/plain entry price 2^16 0.8109 (8), 2^17 0.8103 (8), 2^18 0.8070 (8), 2^19 0.8052 (6), 2^20 0.8031 (4): charged 0.8126 = the largest mean + 2 SE + max(0, slope + 2 SE) x (0.8 L - lsmall) (filt16 TABM, tab_ratio).

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent; with seeded ECM curves the relation count is fixed:
2^64.4827 under H1 (section 3). Under H3 (section 9; H-SM: P Lambda is injective on the chi-part of the units modulo r-th powers, of dimension 9 + 1 by [13]; then W is
injective on its image except with probability <= 10/r < 2^-250, inside the success rounding) and H-X, Block Wiedemann fails with probability <= 4N/2^128 < 2^-66 (Schwartz-Zippel, conditional on H-X: the relevant minors are nonzero polynomials in (z, X) for this matrix and block support; they
have degree <= N in z, <= N + 9 in X, entries uniform on >= 2^128 values); (w, c) is used only if A (w, c) = 0; a failed solve (checked on the dropped rows) is a failure.

Reconstruction. Purged ideals' logs come by back-substitution; ideals with <= 2 occurrences are log-less in the descent
([6b]: model >= 6.3x measured; DESCENT.md).

## 6. Individual logarithms

DESCENT.md ([7], [7e], [7l]): 62 nodes + boot (215, 304, 15, 4), retention >= 0.6490.

Expected two-log work <= 2^108.841278 units ([7]; 2^108.841277 rounded up; one log 2^122.129 instructions: boot own
2^119.80 (whole balls, DESCENT.md) + subtree 2^121.20 x retry 1.520). The meter caps both logs at t = 2^121.9 units (reached = failure). Metered work is min(D, t) <= D for the uncapped two-log work D: its expectation is <= 2^108.841278 for any t, and
P(cap) <= E[D]/t <= 2^(108.841278 - 121.9) <= 0.0001173 (Markov); with 10% more expected work, <= 0.0001291.

## 7. Success

Core: success >= (1 - 2^-66 - P_CAP) x (1 - 2 x 4.10e-09 - 0.0001173 - 2^-12.87 - 2^-10.70 - 2/r) >= 0.99892 (main
phase, two boots with all 18 relations failing, cap, node counter, repeats); >= 0.998911 + 2^-48 with the cap bound at +10% (0.0001291).

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX = 2^27.3. E[N] uses the same
  recurrence with unit node cost (counts16.py): E[N] <= 2^14.430 for both logs, P(abort) <= 2^-12.87.
- Repeats: N_MAX nodes in the smallest bin (2^70.30 ideals), 2^6 for 1/p weighting: P <= 2^-10.70, charged.

Subtrees are independent (H5). Past the gate P(success) >= 0.998911 with the wrapper (the +10% figure); every abort and a closed gate fail: P(success) >= (K/2^64) x 0.998911 >= 0.39. W <= 1 + 1[u < K] W_run, u independent of W_run: E[W] <= 1 + (K/2^64) E[W_run] (audit: gate, expected_run).

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (section 3, H-EA; estimate 112.552654692 + 0.916% allowance) | <= 112.565809621 |
| linear algebra (product form, weight 200, 4 owners, N_l <= 2^56.3580, N_sw <= 2^57.5052, (8192, 130), block-diagonal X, tables; section 4) | <= 110.677012532 |
| matrix generator | <= 100.549594 |
| both individual logarithms (metered cap; expected 108.841278) | 121.9 |
| canonical a list | 92 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 74.6 |
| wide factor-base lattices ([2^63, 2^64), A and coset, KERNELS.md) | 71.6 |
| rec7 bases (2 per ideal x run and two cosets, p in [2^48, 2^64), KERNELS.md) | 75.09 |
| coset replay lists (section 3; count, place, buckets, key order of the tables for p < 2^48) | 97.26 |
| band split (section 3, [10t]) | 101.80 |
| P_g product tree (section 3) | 74 |
| g and f band-2 split of kept relations (section 3; expected) | 97.0 |
| numerical certificate construction | 50 |
| verification | 11 |
| wrapper (section 1; <= 2^80 instructions) | 66 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^92 (2^42.8 per a with Q_a <= s_Y theta A^2, both classes: every run's a-values; 2^62.795 of them); Schirokauer maps 2^88 (below); factor base 2^74.6 (2^23.4 per prime, pi(B) < 2^65.42); filtering and
back-substitution 2^80 (2^24.4 per entry against < 2^14). The rows sum to <= 2^121.90283782313792: core meter 2^121.9028378233 (audit-checked); with the two logarithms at 2^108.841278, not the cap, <= 2^112.9950699693.

Schirokauer maps. R = (Z/r^2)[y, x]/(h, f) (rank 72; g: 12), r unramified ([13]); phi = a - bx, Lambda(phi) = (phi^eps - 1)/r mod r on the basis y^i x^j of F_r[y, x]/(h, f); P = (1/4) sum_(k<4) chi^-k tau^k (chi = p^3 mod r; tau(y) = sigma(y), tau(x) = c x: common12 SIGP, C2P) projects on the chi-eigenspace; the 10 columns are W_f P Lambda_f (9 rows) and W_g P Lambda_g (1 row, same on g), W_f, W_g uniform over F_r, drawn at the start of the run before collection; eps = lcm(r^d_i - 1) < r^72 (18,360 bits; g: 3,060), by left-to-right powering on eight-limb coefficients: < 2^21.31 instructions per bit (2,628 schoolbook 8 x 8-limb products, 72 Barrett reductions by h and f, times phi), < 2^36 per relation, 2^87.04 units at 2^65.33 relations, within 2^88.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (0.29, 0.75), exact v_2 and kappa, exact sampled norms; boot pair ratio kappa <= 2^24, DESCENT.md). H1' retention (Euler-product model, two-stage law; band charge 0.002 bits). H-EA (work only): H1''s survivor law prices the survivors of both stages ([12w]). H-FFT (work only): F(N) bounds the local SSA work past 2^26 limbs (KERNELS.md). H2 region (uniform pools, exact flags, runs T >= 48, measured counts ([3], 6 SE; certified bracket as fallback); sized on S, judged on G + H + K; S+G+K work with the 0.916% allowance). H3 matrix (r unramified, [13]; r prime to the class numbers, H-SM (section 5), a one-dimensional kernel; N(250)/U'
<= 0.00835, production means N_l/U' <= 0.0037, N_sw/U' <= 0.00828 (they do not rise with scale; above the section-4 caps the run aborts), [10c], [10p], section 4). H4 kernel cost (KERNELS.md, STORAGE.md; FEAT_LSE); H4u (unit share >= 86%;
k5 <= 8.29, PREREG-U, [10]; the dense-band price ratio rises at most at its fitted slope + 2 SE to lsmall = 0.8 L, [10d]); rec7's mean per-hit price in each density stratum is at most the mean + 2 SE of its 512 production draws (REC7_K, KERNELS.md); contaminated runs hold at most 0.008%
of the canonical a (95% bound; a shortfall changes cost only); f/g dependence of stage-B survivors within q_g's 2x ([12q]). H-X (premise): for the live system, the 10 fixed S columns and this X support, at the charged d, e, [B z, B Y, .., B^d Y] has rank N and (X; X B; ..; X B^(e-1)) is injective on span S + Im B. Support, not proof: [10x] checks the shipped record of our unplanted two-sided NFS analogues (quartic and quintic, 7 real Schirokauer maps, 4 merged matrices, N 2400-4000, this support): at d = ceil(D/n) + 1 (fixed columns in n), e = ceil(D/m) + 1 every case has full ranks and recovers the kernel, d_min = d - 2, e_min = e - 1, m <= 2048, n <= 149, n/m down to 1/64; one-block pads fail; s = 7. The live operator is a further-merged matrix of the same kind). Re-runnable on this package's own live system: hxgen.py (hosted hx_analogue; filtsim 2^12..2^14, uniform S, this support, charged lengths) gives rank C_d = N and rank O_e A = rank A at every seed tried, and the one-block pad control loses rank. H5 descent (DESCENT.md): lower sums; per-ideal yield certificate (DISC); lazy order independent of outcomes; negative-binomial recurrence, independent
subtrees of distinct ideals (measured: [7e]); missing logs (Poisson occurrence at the 2-SE input); local omitted-mass laws; random-integer ECM orders; UDUP = 0.65 (DESCENT.md, [11]).

## 10. Memory

Core storage < 2^112 bytes, of which the replay lists 2^110.47 during collection (STORAGE.md); the descent peaks at 2^105; wrapper below 2^30: below 2^112 bytes. All preprocessing is for this target,
charged once in W.

## Provenance

GT core: kaankolcu recover-P 1b8fa720 and its lineage (9332a74d mpjunior92; 2c2ce2ab / 42e80f80 gbotrel; 6282ab79
bernardo-lima12; 6734bd5a GPT-6 in Codex; Gemini 45ce22e5 / 4bb6c381); polynomials after Guillevic-Singh; filtering after CADO-NFS. H-EA, H-X,
rec7, band strip and split, f-only sieve, g batch, twisted C4 field: gbotrel, Claude Opus 5.5 in Claude Code. Carry-merge _combine: yukonautominer (9a7df12c). Run gate, expected total work, high cap: kaankolcu (c613ae01, 3359b9e2, 445c6c9d, ad42e9eb); logarithms at their expectation: also mpjunior92 (ea928ab4; and blocks (8192, 130)); residue price 0.02, the same norm-kernel design, two-level thresholds: kaankolcu (538bf09c, 8b2b8442; 01aabb8f, a02dda7e). Two-stage detector: gbotrel, Claude Fable 5.2. Truncated-Shoup _chi_mul: soispoke (cf131637), as shipped by kaankolcu (40b2df68). Product form (order-label recorded pivots): gbotrel, Claude Opus 5.5 in Claude Code. descent_witness.py, runs16.py, threshold_model.py, embedding_bounds.py: lzma behind same-name loaders; the decompressed texts are byte-identical to the accepted c31p2 files (sha256 of the text: threshold_model 3968e2318ffbfe18, embedding_bounds 69626dc5d869530b, descent_witness ced2c04c019ab9c4, runs16 21351d4a8e282d9b); ideal_canonical.py, ecm16.py, filtsim.py now load from lib16.b85, text byte-identical to the accepted c33 files (ideal_canonical 1cbfd5251b957bf7, ecm16 e75b63b5f03c7d12, filtsim e994a67805664eb6); retain16.py likewise (retain16 75afe5ff5e953000; pinned through lib16.b85 in refinement.json).


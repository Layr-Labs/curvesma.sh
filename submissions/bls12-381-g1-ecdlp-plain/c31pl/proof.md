# Plain G1 ECDLP (bls12-381) through two GT logarithms with the STNFS core behind a run gate: expected work 2^115.992303309502 units

Track bls12-381-g1-ecdlp-plain. Input: T = [x]P0 in G1; output the scalar x in 0..r-1.
Score unit: one batched affine G1 addition, C = 2,000 scalar AArch64 instruction equivalents. A core unit is 20,000
instructions = 10 score units. No auxiliary powers, advice, smaller targets or amortization
are used.

## Claim

W is the total classical work of the one attack on the one target. Step 0 is a run gate (kaankolcu's device, Provenance): draw a 64-bit word u;
unless u < K = 7201352326197334294, stop as a failure (< 1 unit). K is the least integer with K x 0.999011 >= 0.39 x 2^64. Under H1-H5, H-EA, H-X and
H-FFT (section 9) the algorithm below returns x with probability >= 0.999011 once past the gate (section 7); u is independent of the challenge and the
other coins, so success is at least 0.39 over uniform challenges and fresh coins, and

    E[W] <= 1 + (K/2^64) x 10 x (2^114.0274016541 + 2^66) < 2^115.992303309502 units (C = 2,000).

2^114.0274016541 core units: section 8's rows: every-run bounds (collection and linear algebra are the same in every run; their measured means carry + 2 SE) except two at expected work, the two logarithms, 2^109.277 (section 6), and the ECM split of kept relations (section 3). The score 115.992303309502 is an
expectation, not a bound on every run: all work of a run (preprocessing, both logarithms, failures, restarts, field map, both pairings, scalar multiplication,
verification, output) is at most

    W <= 10 x (2^121.9059149797 + 2^66) < 2^125.227843074588 units.

Analytic, exploratory claim.

Conventions. Bounds round in the safe direction; evidence16.py checks each against its unrounded value
(declared_bounds) and reruns the checks cited [n]. Model:
stats16.b85 (config, pool_d), the constants of model16.py and descent16.py (owner10.install: C_ROW
572).

## 1. Reduction (fixed base g = e(P0, Q0): two GT logarithms) and field map

With the pinned reduced pairing e and the generators P0 of G1, Q0 of G2:

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

Candidates: (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69.5; factor base: relative-degree-1 ideals of
norm <= B, U' = 2^62.941 tau-orbits (U/2: orbits of 4 ideals; U = 2^63.941).

Region. Exact rational forms Q_a, Q_b (stats16 config), A^2 = 2^21. Pool (common12.regenerate_pair): uniform on theta = 3/8 (a5 > 0, Q_a(a) <=
(3/8) A^2; b uniform in the half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a)). Region, by 2-adic class of a. Cells with a in P2 (class Y, flag p2a; X: the rest) carry exactly 2^6 in N_f whatever b is (below; all 144,173 primitive pool rows): a Y slice is richer per cell by the rows' own yield law (model16.a2). Class scales s = (47/50, 103/101): a5 > 0, Q_a(a) <= s_c (2608/8000) A^2; b5
>= 1, Q_b(b) <= lam_c A^2 - kap Q_a(a), lam_c = s_c 4/5, kap = 16/17 (exact), so the restricted sample is uniform on it. G has 22808 region rows (14223
primitive) of 87197 canonical draws, H 45763 (28471) of 174159; sums are scaled by pool points over all draws. Weights: w = (1 - qa)^3 (yield), vb (1 - qa)^3 x factor (work, population);
slice-geometry costs use S' = vb (lam_c - kap qa)^3 x the lower factor and K' = floor(sqrt((lam_c - kap qa) A^2 / Q_b(e_5))) lines. The
canonical a is the (Q_a, lex)-least element with a5 > 0 of its <+-1, units, sigma> class (exact enumeration of tO, t in {a, sigma a}); with b5 >= 1, one point per <+-1, units, tau> orbit (tau^2(a, b) = (a, -b)); non-primitive pairs are excluded by an exact HNF test. Counts ([3]). Points per volume of each half-ellipsoid, measured by the sampler's acceptance experiment (uniform integer points of the frame box of common12.uniform_half_ellipsoid, exact hit test; 2.5e11 trials per box): Q_a (pool box, every t/8000, t = 2300..3000): 0.999206..0.999331, bracket [0.999062, 0.999429]. Q_b (four boxes, 32,840 radii^2 in [0.44 A^2, A^2], step 2.5e-5; count nondecreasing, volume ~ radius^6): [0.998999, 0.999595] at every radius^2 in [0.44 A^2, A^2] (all pool slices (1 - qa) A^2 and region slices (lam_c - kap qa) A^2: s_c (4/5 - 16/17 x 2608/8000) >= 0.44, s_c 4/5 <= 1, [3]; a: the pool ball). Each point carries 6 SE (SE = ratio/sqrt(hits) <= 2.4e-5; union bound over 67,082 one-sided events < 1.03e-3), rounded outward. Factors: 0.998061 = 0.999062 x 0.998999 (yield, run cells, S') and 0.999595 = 1 x 0.999595 (work, population, D_SUP; a's end raised to 1, so each product also bounds b alone). Fallback: lattice_count_bounds.py certifies (exact rationals) [0.975583, 1.024127] (Q_a), [0.988357, 1.010732] (Q_b) at every radius >= 801; [3] checks the brackets inside it and replays 2^17 trials of two boxes.

Pools (PREREG-G; rows regenerated at this region). 656,000 uniform pairs each on theta = 3/8: S (seeds 261-264) sizes the region; G (251-254, drawn
after the configuration was fixed), H (255-260, 265, 266: 1,312,000 pairs, pre-registered after G) and K (2801-2816: 2,624,000 pairs, drawn after c27 was filed) are
judged, every seed kept; S, G and K give the work; the region keeps a only if its run has T >= 48 (model16.RUN_TMIN); no draw dropped. pool_d (stats16) holds per-seed counts (H, K: totals), list hashes and sums, 8 anchor rows of G, G's canonical lists over draws 0..49 ([2] replays them) and the contamination sample. theta is the smallest t on the 1/8000 grid with margin_S
>= 0.06 + sqrt(2) SE_S (pre-registered), on S with the exact 2-adic yield: 2608/8000 (0.104630 >= 0.104216; 2607: 0.103978 < 0.104235; [4s]). The band
charge applies to the judged count only.

Sieve and survivors (f ideals only), two stages. Stage A. Line sieve (scores only) for p in [2^24, 2^48); record walk for p in [2^48, 2^64) (slice12,
KERNELS.md; X2 = 2^64; per clean a-run, below; [2^63, 2^64) of contaminated runs on the wide path). Survivors pass the f threshold at slack SA = 356.5; the norm kernel (C_NORMS = 48,300 measured with L and the compare; KERNELS.md, [5t]) gives N_f and the g norm z. Recheck. The scan hands each survivor's score s on; the kernel returns L (L <= 16 log2 N_f < L + 2) and the survivor goes on only if L - s <= 16 SA. Each of the <= c scored ideals (c = floor((n_f - v_2)/24)) is rounded down by < 1/16 bit, so L - s is within (-2, c) of 16 x the omitted mass (N_f outside the scored ideals): every cell with omitted mass <= SA - (c + 1)/16 goes on, and one that goes on has omitted mass < SA + 1/8.
Stage B (rechecked survivors). y = P24 mod N_f by the scaled remainder tree (lemma below; P24 = the product of the primes below 2^24; waves of 2^15 norms); with N_f = 2^v z', m = 2^v gcd(z', y^(2^10) mod z') is the part of N_f made of primes below 2^24 (each such prime of z' divides y, and v_p(z') <= 832 < 2^10; model16.k2, [10v]). With Lm = floor(16 log2 m), O = L - s - Lm is within (-2, c + 1) of 16 x the outside mass (N_f outside m and the scored ideals). The survivor goes on iff O <= 16 SB, SB = 276.5, and 16 x 64 k - 1 <= O <= 16 x 69.5 k + c + 1 for some k in 0..4. Window lemma: if N_f is B-smooth and its primes in [2^24, 2^64) divide it only through scored ideals, the outside part is k primes of (2^64, B], of mass in (64 k, 69.5 k]: such a cell with outside mass <= SB - (c + 1)/16 goes on (k <= 4: 5 x 64 > SB); a cell that goes on has outside mass < SB + 1/8, within (c + 3)/16 of a window. Cells where a prime >= 2^24 divides N_f beyond its scored ideals (a square, an ideal of degree >= 2) are outside the yield law (squarefree sieved part). N_f carries 2^0 (class X) or 2^6 (class Y) exactly (below): on Y slices T and the recheck use SA + 6 (stage B needs no offset: m holds the 2^6), so every test is on the odd part.
Laws (threshold_model.two: an Euler product without the prime 2, read at n = n_f - v_2; primes below 2^24 small with all powers, [2^24, 2^64) sieved, above large; explicit below 2^14, prime density above). Yield: r_f = P(sieved and large parts squarefree, small + large <= SA - (c + 1)/16, large <= SB - (c + 1)/16 - 1/8 | smooth), 0.8650 at n = 740 ([2w]). Work (H-EA; all integers, all powers): norms on survivor_upper(n, SA, 32), stage B on survivor_upper(n, SA, 1/4), the g batch and later steps on the same law with large <= SB + 1/4 inside the windows widened by (c + 2)/16 + 1/4. [12w]: on 4,000,000 random integers with exact factorisations at two toy sizes (law with the prime 2; stored counts, first 200,000 of each replayed) both stages' survivors are below the laws, the kept smooth integers above r_f. [4b]: r_f recomputed with the exact odd parts below 2^24 of 8,000 pool norms in place of the law's small part (stored 2-bit histogram, upper bin edges; H1': the small part is independent of the rest) is 0.9017 / 0.8679 / 0.8293 at n = 700 / 740 / 780 against 0.8995 / 0.8650 / 0.8257; at - 2 SE the loss is <= 0.0003 <= REL_BAND_LOSS = 0.009 bits. C_SB = 820,000 per rechecked survivor (KERNELS.md, [10k]). Survivors with B-smooth z get trial division by the 82,025 primes below 2^20 on both cofactors (2^23.649 <= C_TDIV = 2^24.5, [5t]). Batch band strip ([10b], KERNELS.md). Survivor sides wait in
waves of K = 2^39 (<= 832 bits). P = product of the primes in (2^20, 2^48] (< 1.000028 x 2^48/ln 2 bits, Rosser-Schoenfeld); a scaled remainder tree
(Bernstein; model16.srt) gives P mod z, the band part is gcd(z, (P mod z)^(2^6) mod z) (2^6 >= 41 >= any band exponent);
the f band-miss factor (1 - 41 x 0.001^2) stays as a margin. Records carry scores only; a g-passed f side gets a second strip, P2 = the primes in (2^48, 2^64] (exponents <= 17 < 2^5; K = 2^54, C_BATCH2 = 2^25.25 >= 2^25.2401).
g batch ([10g]). Same SRT with P_g the product of all primes <= B (< 1.000028 B/ln 2 bits), waves of K = 2^61 g norms (<= 576 bits, 9 limbs): z is B-smooth iff (P_g mod z)^(2^10) = 0 mod z (exponents <= 576 < 2^10). r_g >= 1 - 2^-38 (ideals of norm p^d > B drop a relation; in the margin). A kept relation's g norm is split by trial division, band split and ECM (<= 12 primes of (2^48, B]): <= 12 x 1.002 x 2^39.89 + 13 C_MR = 2^43.48 per relation, <= 2^65.15 relations ([10t]). Its f band-2 part (<= 17 primes) is split by ECM until found (1.002 x 2^41.06 at 840 bits + 2^26
lookup each, + 18 C_MR(840)): with g 2^96.41 <= 2^96.5 units (ledger row): an expected work (ECM model: independent schedules, each finding its prime with probability >= 1/1.002), metered at the row (reached = failure): Chebyshev (count variance <= 0.00201, a schedule <= 2^26.8 units) gives P <= 0.00201 x 2^26.8/(0.06^2 x 2^96.5) < 2^-70, within the rounding of 0.999011.
Lemma (exact SRT). Product tree Z_v (L <= 61 levels); node v holds Y_v of b_v = 64 (n_v + G) bits (n_v limbs of
Z_v, G = L + 2). Root: Newton quotient q of P 2^b by Z, corrected by <= 2 add/subtracts, Y_root = q mod 2^b. Child c,
sibling s: Y_c = floor((Y_v Z_s mod 2^(b_v))/2^(b_v - b_c)). Leaf: r = ceil(z Y 2^-b) (0 if r = z) = P mod z. Proof: eta_v = (frac(P/Z_v) - Y_v 2^(-b_v) mod 1) Z_v gains < 2^(-64G) per level, so eta < (L + 1) 2^(-64G) < 1 (no wrap). Cost (H-FFT; products <= model16.fft_F, gbatch_cost): C_BATCH = 2^24 (2^23.921, K = 2^39), C_GBATCH = 2^24.98 (2^24.973, K = 2^61). [10v]: model16.srt on 2^8
planted leaves (Newton offsets -2/+2) vs Python's %, K2 on them. P's, P2's and P24's builds (< 2^80 instructions) are misc, P_g's tree (< 2^88) the 2^74 row.
Band split ([10t]). A kept relation is split: the band parts are trial-divided by every prime of (2^20, 2^48] (< 2^43.272; 180 instructions per prime with a 2^(64 i) mod p table, misc row), each hit divided out to full power (<= 41 + 28 divisions, 2^10 each). Each found prime gets its ideal and exponent from the factor base (<= 2^26): 2^50.7634 instructions per relation, none lost. Kept
relations <= the judged 2^63.13438 + 2 SE over the count envelope
(0.999595/0.998061), the ECM miss (0.999^22: <= 13 + 9 primes above 2^64), the band margin and the least
f retention on the model's norm grid (0.7317): <= 2^63.6281, 2^101.1038 <= 2^101.27 units (ledger row; audit).
A stage-B survivor's f cofactor above 2^64 is < 2^284; ECM is the ecm16.py schedule at 99.9% per
prime ([12]): 69.5-bit primes in a 284-bit cofactor cost 2^38.1877 instructions.
A smooth side's cofactor C (no prime below 2^64 after division and both strips) is 1, a prime <= B or composite; ECM runs only then. Compositeness: 40 Miller-Rabin rounds, random bases (C_MR = 40 x 2 x 284 five-limb
products = 2^22.995); a composite passes with probability <= 4^-40. Every stage-B survivor with smooth z is charged the whole schedule. g first: the f side's division, strip, test and ECM run only if the g batch passes:
per threshold survivor C_NORMS, per rechecked one C_SB, per stage-B survivor C_GBATCH + q_g (C_TDIV + C_BATCH + C_BATCH2 + C_MR + 2^38.1877), q_g = min(1, 10 rho(u_g)) per pool row
(rho(u_g) is the yield law's g factor, H1; the 10 covers P(z smooth | the cell passed stage B) / rho(u_g). The two norms of a coprime pair depend on each other through the primes dividing both; for the primes below 2^24 [12q] measures it on 8,000 real pairs (exact small parts, correlation -0.007 +- 0.011): tilted as the two events tilt them, and at twice those tilts, joint over product stays below the factor); ~2^-18.8, 2^-21.1 and 2^-24.6 of them per cell at n_f 740. S, G and K work sums use it on
all 3,936,000 draws (pool_row, runs16).

Thresholds, two levels (kaankolcu's idea, Provenance; model16.two_level). T = lower - 16 x slack, lower a certified lower bound of 16 log2 N_f (embedding_bounds.py: Taylor shift per place, outward rounding; NormBounds.block) for exactly the rows and cells it is used on: every B-smooth cell with omitted mass below slack is kept whatever the bound's width; a wide bound only sends more cells to the exact norms (charged at slack + 32 bits). Level 1: the rows of an aligned 8 x 8 x 8 box (model16.box_rows) share the hull (EB.hull) of the hulls of its 2 x 2 x 4 sub-boxes; in an aligned 256-cell block each 64-cell node meeting the box's cells is evaluated on it; certified (finite, <= 32 bits wide) it is a leaf for all rows, else halved, down to 16 cells. Level 2: a failed 16-cell node runs EB.box_leaves per sub-box (hull, bisection to cells, then rows; uncertified cells are rejected: none in the pools). Merge: with no level 2 in the block and (largest upper - smallest lower) <= 32 bits, the block is one leaf at its smallest lower bound (a lower bound for each cell). A row's part of a leaf is one scan leaf at 256 cells, else cut at multiples of 64 (measured scan lengths). Counts are exact per pool row: 1.6e-04 evaluations, 0.0106 scan leaves per cell. Prices: C_NB_EVAL 1,104, C_NB_ROW 378 (measured, KERNELS.md, [10s]). NormBounds.row: per slice and place c_j = mid e_j, w_j = up((rad e_j + |c_j| 2^-49) UP);
a row takes m = sum b_j c_j, r = (sum |b_j| w_j)(1 + 2^-48): nine roundings < 2^-49.8 relative, so |b(v) - m| <= r. [2b] checks exact f norms of six anchor blocks and two level-2 blocks against their leaves, [2] and the audit the anchors' fields.

Record walk per canonical a-run (rec7; KERNELS.md). A canonical a's run is the maximal integer t-interval
containing 0 with a + t e_0 canonical and inside the limit s_c theta A^2 of a's class c; rec7 walks its sub-run of class c (runs16.py). On a run, a prime ideal (p, s, r)'s hits in all its slices form one coset of a 7-dim
lattice (t = r b(s) - a(s) mod p), walked once per run and ideal pair (p, s, +-r) with one setup6 ([10r], KERNELS.md).
Clean-run rule: rec7 walks a run only if all its a-values are in the canonical list (tested in the 2^92
row); otherwise its slices take the slice walk. A clean run of 32 to 2,048 a-values uses the run walk (x-walk-asm: our hand-written AArch64 run walk, source not shipped; REC7_K x the reference curve, by its 512 x 10 campaign at this theta, 0 mismatches vs the reference walk; KERNELS.md) in each band
p in [2^48, 2^64) with h' = 2 x run cells x lower/upper ([3]) / p >= 2^4.99 (the curve's bottom) instead of the slice walk and empty cosets; run
cells are exact (runs16.lrs: integer sum over the run's a-values of (lam_c - kap qa)^3, halved). Each row pays its cells' share
(per hit at h'; setup6 per (run, pair) replaces its slices', f ideals). Contamination. A contaminated run's slices cost at most
the slice walk: per slice <= D_SUP = contam_sup = 2^72.639 on [2^48, 2^64) (exact band maxima with empty cosets and setup6;
[k7c]: slices peak at 2^72.330). Every row pays RUN_ALLOW x D_SUP/2 per slice (f ideals), RUN_ALLOW = 0.00091 (draws are
uniform over canonical a): in 4,000 a-uniform runs (PREREG-C: seeds 5001-5010, the first 400 canonical a inside their class limit
per seed, runs at that limit, sha256-16 in pool_d) the exact predicate found no non-canonical a-value among 1,966,129
(675,989 new, the rest in runs inside their run tested at 2565/8000; its 47 runs above 1,024 a-values hold 51,313):
contaminated runs hold under 1 - 0.05^(1/4000) = 0.00074866 <= 0.00091 of the canonical a (95%). [k7c]/[2r] replay 6 sampled/3 anchor runs.
Coset sieve at P2. 2 = P2^2 P4^2, P2 = wO; b is in P2 iff P2F.b is even (common12.P2F; [10r] checks wO = ker P2F).
If a and b are in P2 the pair is not primitive, so those cells are not sieved; P2 is fixed by units, sigma, tau and
b -> -b, so the canonical domain is unchanged. In the basis e_0, e_j - (P2F.e_j) e_0 a slice with a in P2 sieves
b_0 = 2k + 1 (half the cells, congruence constants halved mod p, setup6 + 10). P2F(a + t) alternates along a run, so
rec7 walks it as two cosets (1/2 and 1/4 of its cells, setup6 + 10 each). model16.evaluate: p2a-and-p2b rows carry
nothing; p2a slices are priced at S/2 or the sub-run's cells.
S' uses the lower count factor ([3]), so the per-a charge is >= 0.00091 D_SUP >= 0.000749 D_SUP.
runs16's run ends (< 2^69.3 units) are misc.

Yield ([4], pools G + H + K). rels = (va_pool/4592000) vb (0.998061) sum over their canonical primitive region rows of
w rho(u_f) rho(u_g) r_f 0.999^lambda_f x the f band margin (g exact): Dickman rho ([1]), retentions r from the
two-stage law above (threshold_model.py), the ECM miss; cells with qa < 0.02 are sieved and charged, relations discarded. Estimate
2^63.13437 (SE 0.0192 bits) less the band charge 0.009: 2^63.1254 against 1.05 U' <= 2^63.01108: margin >= 0.1142
bits = 5.9 pool SE (>= max(0.06, 2 SE)). Band charge (H1'): [4b], above. Descent input. LOG2_RELS = 63.074 is <= the charged count
lowered by two relative pool SE, 2^63.08651 ([4] rels_input); the descent's mean norms 738.2/498.2 are not above the judged
739.834/498.338 ([4]; smaller means and count only lower kappa and lambda, DESCENT.md).
Exact 2-adic part (H1, [14v], [14t]). rho prices a random integer (mean v_2 = 1); each row's v_s = v_2 of its exact norm
is known (pool_row), so u_s = (n_s + alpha_s + E_s + kappa_s - v_s)/69.5 (model16.a2, also q_g): E = (2, 2.399985) >= the
local means (runs16; alpha_2 = 1 - E), kappa = (0.06, 0.05) >= kappa2(u) (sum_k 2^-(k+1) rho(u - (k-1-kappa)/69.5) = rho(u),
rising in u) at n_f 840, n_g 530, above every yielding row. Toy census (92,108 primitive (a, b) in {-1,0,1}^12, exact 2^20-smooth N_f): class chi^2 at 2^18/20/22 10/6/18 (exact) vs 88/30/23 (average); totals within 1 SE.

Cost ([5]). Sampled mean 19.66 instructions/kept cell over pools S, G and K, charged 19.80 with the allowance, over
2^123.779 kept cells (every step at KERNELS.md/STORAGE.md constants; rec7 at its measured stratum means + 2 SE): collection <= 2^113.798223000
units (estimate 2^113.788388505 times allowance 1.00684). Allowance (H2, [4a]): a Monte Carlo estimate over S, G and K's
3,936,000 draws, relative SE row 0.2164%, incidence sqrt((1-144691/3936000)/144691) = 0.2581%, total <= 0.3419%;
charged 1.00684 >= 1 + 2 SE.

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U', filtsim.PEXC), clique removal and Markowitz merge (order: (rows - 2) x lightest unit pivot)
are simulated on the production occurrence profile (filtsim.py; [6] reruns 2^12) as the tau-quotient system: tau^k(I) enters I's orbit column with
coefficient chi^k (entries a + b chi, merged |a|, |b| <= 255). Cap. All retained relations are used; H1 puts the count at >= 2^63.08651
>= 1.09 U' ([4]): filtering at excess 1.09, fixed by the rule margin >= max(0.06, 2 SE) on the judged count (0.005 grid; sizing and judgment stay at 1.05; margin over
1.09 U' 0.0603 >= max(0.06, 2 SE), [4]). Collection's ECM curves are fixed by a public seed, so
relations, filter and N are deterministic in the challenge-independent instance, and H3 is stated for it:
production N(250)/U' <= 0.01416 (supported by the campaigns below, <= 0.01396, and the [10c] spread law); the attack aborts if not (P_CAP = 2e-05 charged as slack); LA is
charged at 0.01416. [10c]: Chebyshev with sd(x) <= 0.046 U'^(-1/4) (2 x the largest per-run sd x U^(1/4); sd x sqrt(U) is 0.26-0.67: no U^(-1/2) law). A campaign at B = 2^69.5, excess 1.09, merge weight 250 and relation
means 740.3/498.7 (not below the judged pools' 739.834/498.338, [10y]) measures x = N(250)/U' on 8-16 seeds per scale: 2^16 0.01765, 2^17 0.01686, 2^18 0.01609, 2^19 0.01515, 2^20 0.01477, 2^21 0.01422, 2^22 0.01360, 2^23 0.01317 ([10] lists SE): it falls at every step. Correction: the explicit band takes every fourth prime per side ([::4], occurrence min(1, 4 kappa/q)), the pi(B)/2 orbit model of the other bands. The tower's actual tau-orbits (all degrees, norm <= 2^19.5) as band scale x by 1.0062 +- 0.0048 (2^18-2^22): NU_H carries 1.0159.
H3: production is not above the largest scales: the larger mean + 2 SE of 2^22 and 2^23 is 0.013736, x 1.0159 (splitting) <= 0.01396. FILT_RUNS (excess 1.05, means 744/501, weight 400; 2^15, 2^19, 2^20 kept) is the older campaign. Merge to mean weight 250 (charged 251; the stored runs' mean weights lie in [250, 251], [10y]). N <= 0.01416 U' <= 2^56.7987 < 3 x 2^64/320 (four owners, STORAGE.md).

The square system A = [M | S] has N + 10 rows, N orbit columns and 10 dense Schirokauer columns: the chi-eigenpart of
the unit groups, (rk K_f - rk K_f^+ + rk K_g - rk K_h)/2 = (35 - 17 + 5 - 3)/2 from signatures K_f (0, 36), K_f^+ (0, 18),
K_g (0, 6), K_h (2, 2) ([13]; tau^2 = -1 on the sigma-odd part, so its chi and chi^-1 eigenspaces have equal dimension);
under H3 its kernel is the line of virtual logs (w, c), c != 0, and M has rank N. Block Wiedemann (m, n) = (8192, 130) with S as right-hand sides: B = [M | 0]; Y = (S_1..S_10, B z_11..B z_n), z uniform in F_r; X block-diagonal (H-X): the 10 pad coordinates (zero columns of B)
in blocks 1..10, the others shuffled and cut into m blocks, X_l uniform in [0, 2^128) on block l. A generator sum_k B^k Y f_k = 0 of
a_i = X^T B^i Y with f_0 = (c, e), c != 0, gives S c + B w = 0, w = Z e + sum_(k >= 1) B^(k-1) Y f_k: one Horner pass, no product touches S
([10x]: toy end to end). d = ceil(D/n) + 1 + ceil(D/128n), e likewise in m, D = N + 10: d + e < (1 + 2^-7) D (1/m + 1/n) + 6 products and terms, charged; Horner within N/m + N/n (2^-7 m < n); Horner runs as five interleaved chains from checkpoints
B^(jL) Y (L = Horner length / 5; 4n = 520 checkpoint vectors, < 2^72 bytes, 4L < N/m + N/n) at the five-vector rate. Each
vector element is stored with chi x_j. Rows reduce exactly (_row_red). Per row and vector (STORAGE.md, [10]): 251 x 8.29 x 0.9275 (dense-band tables, below; H4u, PREREG-U) + 572 (row tail, four owners) + 200 (_chi_mul); per
Horner element n generator terms and one reduction (270); projections: each vector element meets one
X entry, (1 + n/m) N^2 terms at 60. Linear algebra <= 2^110.835813707 units; generator (matrix Berlekamp-Massey by half-gcd, k = m + n = 8322: a node of l
sequence terms makes (m + k) k^2 products of polynomials with <= l and <= l/2 + 1 coefficients, each one Kronecker integer product in 9-limb
slots (coefficients < k L r^2 < 2^576): F(18 l + 18) under H-FFT, + 640 l to pack, reduce (270) and add; leaves m k^2 (97 + 270). The sum over
the recursion is below model16.generator ([10m]): 4 k^3 log2 L (F(18 L) + 640 L) instructions, L = (1 + 2^-7)(N + 10)(1/m + 1/n) + 6): <= 2^101.027 units.
Dense-band tables ([10d], tables16.py). Per Krylov/Horner step and five-vector slot, the fourth owner copies the band columns' values reduced into [0, r) with negatives and chi multiples and tabulates every multiset of g of them over norms <= 2^31.5 (g = 2), 2^21 (3), 2^16 (4), 4^g/2 sign/chi patterns each (global sign: pos/neg accumulator): 2^61.76 bytes, <= half its free space 2^64 - 320 max(N - 3 x 2^64/320, 2^28) > 2^63.99. Pairs unreduced (< 2r < 2^256), 8.4 per vector-entry; g >= 3 reduced, 18.4 (tab.txt). A row adds one entry per group (unchanged unit_run5) on that owner's (pos, neg) part: B x is unchanged. Entry cost (units 7.8, general 33, builds) with/without tables, filtsim rows at weight 250 (band terms above the explicit primes by its kappa(u, t)/t; Poisson band counts; exact groups): 0.7593 / 0.7673 / 0.7690 / 0.7683 at 2^14 / 16 / 17 / 18; TAB_RATIO = 0.9275 >= 1.005 x the larger of 2^17, 2^18 and >= the largest + (least-squares slope + 2 SE = 0.00376 per doubling) x (0.8 L - lsmall); it scales k5, whose per-entry overhead falls with the adds.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent; with seeded ECM curves the relation count is fixed:
2^63.1254 under H1 (section 3). Under H3 (section 9; the 10 Schirokauer maps of the chi-eigenpart
injective on units modulo r-th powers) and H-X, Block Wiedemann fails with probability <= 4N/2^128 < 2^-66 (Schwartz-Zippel, conditional on H-X: the relevant minors are nonzero polynomials in (z, X) for this matrix and block support; they
have degree <= N in z, <= N + 9 in X, entries uniform on >= 2^128 values); (w, c) is used only if A (w, c) = 0; a failed solve (checked on the dropped rows) is a failure.

Reconstruction. Purged ideals' logs come by back-substitution; ideals with <= 2 occurrences are log-less in the descent
([6b]: model >= 6.3x measured; DESCENT.md).

## 6. Individual logarithms

DESCENT.md ([7], [7e], [7l]): 62 nodes + boot (215, 304, 15, 4), retention >= 0.6490.

Expected two-log work <= 2^109.277 units ([7]; unrounded 2^109.27611; one log 2^122.564 instructions: boot own
2^120.08 (whole balls, DESCENT.md) + subtree 2^121.40 x retry 1.836). The meter caps both logs at t = 2^121.9 units (reached = failure). Metered work is min(D, t) <= D for the uncapped two-log work D: its expectation is <= 2^109.277 for any t, and
P(cap) <= E[D]/t <= 2^(109.277 - 121.9) <= 0.0001586 (Markov); with 10% more expected work, <= 0.0001744.

## 7. Success

Core: success >= (1 - 2^-66 - P_CAP) x (1 - 2 x 7.2e-07 - 0.0001586 - 2^-12.35 - 2^-10.70 - 2/r) >= 0.99902 (main
phase, two boots with all 18 relations failing, cap, node counter, repeats); >= 0.999011 + 2^-48 with the cap bound at +10% (0.0001744).

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX = 2^27.3. E[N] uses the same
  recurrence with unit node cost (counts16.py): E[N] <= 2^14.9427 for both logs, P(abort) <= 2^-12.35.
- Repeats: N_MAX nodes in the smallest bin (2^70.30 ideals), 2^6 for 1/p weighting: P <= 2^-10.70, charged.

Subtrees are independent (H5). Past the gate P(success) >= 0.999011 with the wrapper (the +10% figure); every abort and a closed gate fail: P(success) >= (K/2^64) x 0.999011 >= 0.39. W <= 1 + 1[u < K] W_run, u independent of W_run: E[W] <= 1 + (K/2^64) E[W_run] (audit: gate, expected_run).

## 8. Ledger

| component | log2 core units (20,000 instructions) |
|---|---:|
| relation collection (section 3, H-EA; estimate 113.788388505 + 0.684% allowance) | <= 113.798223000 |
| linear algebra (dense-band tables, weight 251, N <= 2^56.7844, (8192, 130), block-diagonal X; section 4) | <= 110.835813707 |
| matrix generator | <= 101.027 |
| both individual logarithms (metered cap; expected 109.277) | 121.9 |
| canonical a list | 92 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 73.1 |
| wide factor-base lattices ([2^63, 2^64), A and coset, KERNELS.md) | 71.6 |
| rec7 bases (2 per ideal x run and two cosets, p in [2^48, 2^64), KERNELS.md) | 75.09 |
| band split (section 3, [10t]) | 101.27 |
| P_g product tree (section 3) | 74 |
| g and f band-2 split of kept relations (section 3; expected) | 96.5 |
| numerical certificate construction | 50 |
| verification | 11 |
| wrapper (section 1; <= 2^80 instructions) | 66 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^92 (2^42.8 per a with Q_a <= s_Y theta A^2, both classes: every run's a-values; 2^62.966 of them); Schirokauer maps 2^88 (below); factor base 2^73.1 (2^23.4 per prime, pi(B) < 2^63.95); filtering and
back-substitution 2^80 (2^24.4 per entry against < 2^14). The rows sum to <= 2^121.90591497969400: core meter 2^121.9059149797 (audit-checked); with the two logarithms at 2^109.277, not the cap, <= 2^114.0274016541.

Schirokauer maps. R = (Z/r^2)[y, x]/(h, f) (rank 72; g: 12), r unramified ([13]); the 10 columns are fixed F_r-linear forms (chi-eigenpart) of (phi^eps - 1)/r mod r, phi = a - bx, eps = lcm(r^d_i - 1) < r^72 (18,360 bits; g: 3,060), by left-to-right powering on eight-limb coefficients: < 2^21.31 instructions per bit (2,628 schoolbook 8 x 8-limb products, 72 Barrett reductions by h and f, times phi), < 2^36 per relation, 2^86.86 units at 2^65.15 relations, within 2^88.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (0.29, 0.75), exact v_2 and kappa, exact sampled norms; boot pair ratio kappa <= 2^24, DESCENT.md). H1' retention (Euler-product model, two-stage law; band charge 0.009 bits). H-EA (work only): H1''s survivor law prices the survivors of both stages ([12w]). H-FFT (work only): F(N) bounds the local SSA work past 2^26 limbs (KERNELS.md). H2 region (uniform pools, exact flags, runs T >= 48, measured counts ([3], 6 SE; certified bracket as fallback); sized on S, judged on G + H + K; S+G+K work with the 0.684% allowance). H3 matrix (r unramified, [13]; r prime to the class numbers, the unit condition, a one-dimensional kernel; N(250)/U'
<= 0.01396: it does not rise with scale, sd <= 0.046 U'^(-1/4) [10c], section 4). H4 kernel cost (KERNELS.md, STORAGE.md; FEAT_LSE); H4u (unit share >= 86%;
k5 <= 8.29, PREREG-U, [10]; dense-band entry cost <= TAB_RATIO (trend to lsmall = 0.8 L: 0.7690 at 2^17, lsmall 13.48, + 0.00376 x (55.6 - 13.48) <= 0.9275), [10d]); rec7's mean per-hit price in each density stratum is at most the mean + 2 SE of its 512 production draws (REC7_K, KERNELS.md); contaminated runs hold at most 0.091%
of the canonical a (95% bound; a shortfall changes cost only); f/g dependence of stage-B survivors within q_g's 10x ([12q]). H-X (premise, not implied by H3 or invariant factors): for this matrix, the 10 fixed S columns and this X support, at the charged d, e, [B z, B Y, .., B^d Y] has rank N and (X; X B; ..; X B^(e-1)) is injective on span S + Im B ([10x] toys: diagnostics only). Pads need distinct blocks: x^T B^i (i >= 1) vanishes on ker B, so only X sees them. H5 descent (DESCENT.md): lower sums; per-ideal yield certificate (DISC); lazy order independent of outcomes; negative-binomial recurrence, independent
subtrees of distinct ideals (measured: [7e]); missing logs (Poisson occurrence at the 2-SE input); local omitted-mass laws; random-integer ECM orders; UDUP = 0.65 (DESCENT.md, [11]).

## 10. Memory

Core storage <= 2^105 bytes (STORAGE.md), wrapper below 2^30: below 2^106 bytes. All preprocessing is for this target,
charged once in W.

## Provenance

GT core: kaankolcu recover-P 1b8fa720 and its lineage (9332a74d mpjunior92; 2c2ce2ab / 42e80f80 gbotrel; 6282ab79
bernardo-lima12; 6734bd5a GPT-6 in Codex; Gemini 45ce22e5 / 4bb6c381); polynomials after Guillevic-Singh; filtering after CADO-NFS. H-EA, H-X,
rec7, band strip and split, f-only sieve, g batch, twisted C4 field: gbotrel, Claude Opus 5.5 in Claude Code. Carry-merge _combine: yukonautominer (9a7df12c). Run gate, expected total work, high cap: kaankolcu (8ab78e96, b6f0a0f6, d7ca6b98, f1c17abf); logarithms at their expectation: also mpjunior92 (ea928ab4; and blocks (8192, 130)); residue price 0.02, the same norm-kernel design, two-level thresholds: kaankolcu (538bf09c, 8b2b8442; 01aabb8f, a02dda7e). Two-stage detector: gbotrel, Claude Fable 5.2. Carry-merge owner combine: yukonautominer (9a7df12c).

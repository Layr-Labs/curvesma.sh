# recover-Q (bls12-381) through two GT logarithms with the STNFS core behind a run gate: expected work 2^113.375177740138 units

Track bls12-381-fixed-g1-invert-g2. Input: a nonidentity P in G1 and z in GT; output Q in G2 with e(P,Q) = z
for the exact supplied P (Q = infinity for z = 1). The score unit is one dense F_p12 multiplication, C = 20,000
scalar AArch64 instruction equivalents. No advice, smaller targets or amortization
are used. Predecessor work is reused as credited at the end; no prior verdict is a premise.

## Claim

W is the total classical work of the one attack on the one target. Step 0 is a run gate (kaankolcu's device, Provenance): draw a 64-bit word u;
unless u < K = 7201208159453235316, stop as a failure (< 1 unit). K is the least integer with K x 0.999031 >= 0.39 x 2^64. Under H1-H5, H-EA, H-X and
H-FFT (section 9) the algorithm below returns such a Q with probability >= 0.999031 once past the gate (section 7); u is independent of the challenge and the
other coins, so success is at least 0.39 over uniform challenges and fresh coins, and

    E[W] <= 1 + (K/2^64)(2^114.7322330618 + 2^66) < 2^113.375177740138  dense-Fp12 units (C = 20,000).

2^114.7322330618: section 8's rows: every-run bounds (collection and linear algebra are the same in every run; their measured means carry + 2 SE) except two at expected work, the two logarithms, 2^109.277 (section 6), and the ECM split of kept relations (section 3). The score 113.375177740138 is an
expectation, not a bound on every run: all work of a run (preprocessing, both logarithms, failures, restarts, field map, pairing, scalar multiplication,
verification, output) is at most

    W <= 2^121.9097718411 + 2^66 < 2^121.909771841101  units.

Analytic, exploratory claim.

Conventions. Bounds round in the safe direction; evidence16.py checks each against its unrounded value
(declared_bounds) and reruns the checks cited [n]. Model:
stats16.b85 (config, pool_d), the constants of model16.py and descent16.py (owner10.install: C_ROW
662).

## 1. Reduction (fixed P: two GT logarithms) and field map

With the pinned reduced pairing e and the public generator Q0 of G2:

0. Run gate (Claim): draw u; if u >= K, fail.
1. Check P (on curve, [r]P = infinity, P != infinity) and z (z^r = 1). If z = 1 return Q = infinity.
2. g = e(P, Q0) (pinned pairing(Q0, P), protected encoding).
3. Map g and z into the core's tower (field-map.md).
4. Run the core on (g, z) under its meter and node counter (any failure is a failure): a = log_g(z) mod r.
5. Q = [a]Q0; output Q if e(P, Q) = z, otherwise fail.

GT has prime order r and e is nondegenerate: for the supplied P, g generates GT, z = g^a for a unique a mod r, e(P, [a]Q0) = z. For uniform P, g is uniform on GT minus {1}, independent of the uniform z: the core's input law. The
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

Candidates: (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69.5; factor base: relative-degree-1 ideals of
norm <= B, U' = 2^62.941 tau-orbits (U/2: orbits of 4 ideals; U = 2^63.941).

Region. Exact rational forms Q_a, Q_b (stats16 config), A^2 = 2^21. Pool (common12.regenerate_pair): uniform on theta = 3/8 (a5 > 0, Q_a(a) <=
(3/8) A^2; b uniform in the half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a)). Region, by 2-adic class of a. Cells with a in P2 (class Y, flag p2a; X: the rest) carry exactly 2^6 in N_f whatever b is (below; all 144,173 primitive pool rows): a Y slice is richer per cell by the rows' own yield law (model16.a2). Class scales s = (47/50, 103/100): a5 > 0, Q_a(a) <= s_c (2480/8000) A^2; b5
>= 1, Q_b(b) <= lam_c A^2 - Q_a(a), lam_c = s_c 17/20 (exact), so the restricted sample is uniform on it. G has 23945 region rows (14953
primitive) of 87197 canonical draws, H 48043 (29939) of 174159; sums are scaled by pool points over all draws. Weights: w = (1 - qa)^3 (yield), vb (1 - qa)^3 x factor (work, population);
slice-geometry costs use S' = vb (lam_c - qa)^3 x the lower factor and K' = floor(sqrt((lam_c - qa) A^2 / Q_b(e_5))) lines. The
canonical a is the (Q_a, lex)-least element with a5 > 0 of its <+-1, units, sigma> class (exact enumeration of tO, t in {a, sigma a}); with b5 >= 1, one point per <+-1, units, tau> orbit (tau^2(a, b) = (a, -b)); non-primitive pairs are excluded by an exact HNF test. Counts ([3]). Points per volume of each half-ellipsoid, measured by the sampler's acceptance experiment (uniform integer points of the frame box of common12.uniform_half_ellipsoid, exact hit test; 2.5e11 trials per box): Q_a (pool box, every t/8000, t = 2300..3000): 0.999206..0.999331, bracket [0.999062, 0.999429]. Q_b (three boxes, 27,727 radii^2 in [A^2/2, A^2], step 2.5e-5; count nondecreasing, volume ~ radius^6): [0.999035, 0.999595] at every radius^2 in [A^2/2, A^2] (all pool slices (1 - qa) A^2 and region slices (lam_c - qa) A^2: 2 s_c (17/20 - 2480/8000) >= 1 >= lam_c, [3]; a: the pool ball). Each point carries 6 SE (SE = ratio/sqrt(hits) <= 2.4e-5; union bound over 56,856 one-sided events < 8.7e-4), rounded outward. Factors: 0.998097 = 0.999062 x 0.999035 (yield, run cells, S') and 0.999595 = 1 x 0.999595 (work, population, D_SUP; a's end raised to 1, so each product also bounds b alone). Fallback: lattice_count_bounds.py certifies (exact rationals) [0.975583, 1.024127] (Q_a), [0.988357, 1.010732] (Q_b) at every radius >= 781; [3] checks the brackets inside it and replays 2^17 trials of two boxes.

Pools (PREREG-G). 656,000 uniform pairs each on theta = 3/8: S (seeds 261-264) sizes the region, G (251-254, drawn
after the configuration was fixed) and H (255-260, 265, 266: 1,312,000 pairs, pre-registered after G, every seed kept) are
judged, S, G and W (281-304) give the work; no draw dropped. pool_d (stats16) holds per-seed counts (H: total), list hashes and sums, 8 anchor rows of G, G's canonical lists over draws 0..49 ([2] replays them) and the contamination sample. theta is the smallest t on the 1/8000 grid with margin_S
>= 0.06 + sqrt(2) SE_S (pre-registered), on S with the exact 2-adic yield: 2480/8000 (0.104008 >= 0.103840; 2479: 0.103282 < 0.103861; [4s]). The band
charge applies to the judged count only.

Sieve and survivors (f ideals only), two stages. Stage A. Line sieve (scores only) for p in [2^24, 2^48); record walk for p in [2^48, 2^64) (slice12,
KERNELS.md; X2 = 2^64; per clean a-run, below; [2^63, 2^64) of contaminated runs on the wide path). Survivors pass the f threshold at slack SA = 346.5; the norm kernel (C_NORMS = 48,300 measured with L and the compare; KERNELS.md, [5t]) gives N_f and the g norm z. Recheck. The scan hands each survivor's score s on; the kernel returns L (L <= 16 log2 N_f < L + 2) and the survivor goes on only if L - s <= 16 SA. Each of the <= c scored ideals (c = floor((n_f - v_2)/24)) is rounded down by < 1/16 bit, so L - s is within (-2, c) of 16 x the omitted mass (N_f outside the scored ideals): every cell with omitted mass <= SA - (c + 1)/16 goes on, and one that goes on has omitted mass < SA + 1/8.
Stage B (rechecked survivors). y = P24 mod N_f by the scaled remainder tree (lemma below; P24 = the product of the primes below 2^24; waves of 2^15 norms); with N_f = 2^v z', m = 2^v gcd(z', y^(2^10) mod z') is the part of N_f made of primes below 2^24 (each such prime of z' divides y, and v_p(z') <= 832 < 2^10; model16.k2, [10v]). With Lm = floor(16 log2 m), O = L - s - Lm is within (-2, c + 1) of 16 x the outside mass (N_f outside m and the scored ideals). The survivor goes on iff O <= 16 SB, SB = 276, and 16 x 64 k - 1 <= O <= 16 x 69.5 k + c + 1 for some k in 0..4. Window lemma: if N_f is B-smooth and its primes in [2^24, 2^64) divide it only through scored ideals, the outside part is k primes of (2^64, B], of mass in (64 k, 69.5 k]: such a cell with outside mass <= SB - (c + 1)/16 goes on (k <= 4: 5 x 64 > SB); a cell that goes on has outside mass < SB + 1/8, within (c + 3)/16 of a window. Cells where a prime >= 2^24 divides N_f beyond its scored ideals (a square, an ideal of degree >= 2) are outside the yield law (squarefree sieved part). N_f carries 2^0 (class X) or 2^6 (class Y) exactly (below): on Y slices T and the recheck use SA + 6 (stage B needs no offset: m holds the 2^6), so every test is on the odd part.
Laws (threshold_model.two: an Euler product without the prime 2, read at n = n_f - v_2; primes below 2^24 small with all powers, [2^24, 2^64) sieved, above large; explicit below 2^14, prime density above). Yield: r_f = P(sieved and large parts squarefree, small + large <= SA - (c + 1)/16, large <= SB - (c + 1)/16 - 1/8 | smooth), 0.8545 at n = 740 ([2w]). Work (H-EA; all integers, all powers): norms on survivor_upper(n, SA, 32), stage B on survivor_upper(n, SA, 1/4), the g batch and later steps on the same law with large <= SB + 1/4 inside the windows widened by (c + 2)/16 + 1/4. [12w]: the detector's integer tests on 4,000,000 random integers with exact factorisations at two toy sizes (law with the prime 2; stored counts, the first 200,000 of each replayed): the survivors of both stages are below the laws, the kept smooth integers above r_f. [4b]: r_f recomputed with the exact odd parts below 2^24 of 8,000 pool norms in place of the law's small part (stored 2-bit histogram, upper bin edges; H1': the small part is independent of the rest) is 0.8931 / 0.8573 / 0.8168 at n = 700 / 740 / 780 against 0.8909 / 0.8545 / 0.8134; at - 2 SE the loss is <= 0.0015 <= REL_BAND_LOSS = 0.009 bits. C_SB = 820,000 per rechecked survivor (KERNELS.md, [10k]). Survivors with z B-smooth get trial division by every prime below 2^20 on both cofactors (82,025 primes, 2^23.649 <= C_TDIV = 2^24.5, [5t]). Batch band strip ([10b], KERNELS.md). Survivor sides wait in
waves of K = 2^39 (<= 832 bits). P = product of the primes in (2^20, 2^48] (< 1.000028 x 2^48/ln 2 bits, Rosser-Schoenfeld); a scaled remainder tree
(Bernstein; model16.srt) gives P mod z, the band part is gcd(z, (P mod z)^(2^6) mod z) (2^6 >= 41 >= any band exponent);
the f band-miss factor (1 - 41 x 0.001^2) stays as a margin. Records carry scores only; a g-passed f side gets
a second strip, P2 = the primes in (2^48, 2^64] (exponents <= 17 < 2^5; K = 2^54, C_BATCH2 = 2^25.25 >= 2^25.2401).
g batch ([10g]). Same SRT with P_g the product of all primes <= B (< 1.000028 B/ln 2 bits), waves of K = 2^61 g norms (<= 576 bits, 9 limbs): z is B-smooth iff (P_g mod z)^(2^10) = 0 mod z (exponents <= 576 < 2^10). r_g >= 1 - 2^-38 (ideals of norm p^d > B drop a relation; in the margin). A kept relation's g norm is split by trial division, band split and ECM (<= 12 primes of (2^48, B]): <= 12 x 1.002 x 2^39.89 + 13 C_MR = 2^43.48 per relation, <= 2^65.15 relations ([10t]). Its f band-2 part (<= 17 primes) is split by ECM until found (1.002 x 2^41.06 at 840 bits + 2^26
lookup each, + 18 C_MR(840)): with g 2^96.41 <= 2^96.5 units (ledger row): an expected work (ECM model: independent schedules, each finding its prime with probability >= 1/1.002), metered at the row (reached = failure): Chebyshev (count variance <= 0.00201, a schedule <= 2^26.8 units) gives P <= 0.00201 x 2^26.8/(0.06^2 x 2^96.5) < 2^-70, within the rounding of 0.999031.
Lemma (exact SRT). Product tree Z_v (L <= 61 levels); node v holds Y_v of b_v = 64 (n_v + G) bits (n_v limbs of
Z_v, G = L + 2). Root: Newton quotient q of P 2^b by Z, corrected by <= 2 add/subtracts, Y_root = q mod 2^b. Child c,
sibling s: Y_c = floor((Y_v Z_s mod 2^(b_v))/2^(b_v - b_c)). Leaf: r = ceil(z Y 2^-b) (0 if r = z) = P mod z. Proof:
eta_v = (frac(P/Z_v) - Y_v 2^(-b_v) mod 1) Z_v gains < 2^(-64G) per level, so eta < (L + 1) 2^(-64G) < 1 (no wrap). Cost (H-FFT; products <= model16.fft_F, gbatch_cost): C_BATCH = 2^24 (2^23.921, K = 2^39), C_GBATCH = 2^24.98 (2^24.973, K = 2^61). [10v]: model16.srt on 2^8
planted leaves (Newton offsets -2/+2) vs Python's %, K2 on them. P's, P2's and P24's builds (< 2^80 instructions) are misc, P_g's tree (< 2^88) the 2^74 row.
Band split ([10t]). A kept relation is split: the band parts are trial-divided by every prime of (2^20, 2^48] (< 2^43.272; 180 instructions per prime with a 2^(64 i) mod p table, misc row), each hit divided out to full power (<= 41 + 28 divisions, 2^10 each). Each found prime
gets its ideal and exponent from the factor base (<= 2^26): 2^50.7634 instructions per relation, no relation lost. Kept
relations <= G + H's 2^63.16227 + 2 SE over the count envelope
(0.999595/0.998097), the ECM miss (0.999^22: <= 13 + 9 primes above 2^64), the band margin and the least
f retention on the model's norm grid (0.7317): <= 2^63.7238, 2^100.1995 <= 2^101.27 units (ledger row; audit).
A stage-B survivor's f cofactor above 2^64 is < 2^284; ECM is the ecm16.py schedule at 99.9% per
prime ([12]): 69.5-bit primes in a 284-bit cofactor cost 2^38.1877 instructions.
A smooth side's cofactor C (no prime below 2^64 after division and both
strips) is 1, a prime <= B, or composite; ECM runs only then. Compositeness: 40 Miller-Rabin rounds, random bases (C_MR = 40 x 2 x 284 five-limb
products = 2^22.995); a composite passes with probability <= 4^-40. Every stage-B survivor with smooth z is charged the whole schedule. g first: the f side's division, strip, test and ECM run only if the g batch passes:
per threshold survivor C_NORMS, per rechecked one C_SB, per stage-B survivor C_GBATCH + q_g (C_TDIV + C_BATCH + C_BATCH2 + C_MR + 2^38.1877), q_g = min(1, 100 rho(u_g)) per pool row
(H1; 100 a margin); ~2^-18.8, 2^-21.1 and 2^-24.6 of them per cell at n_f 740. S and G work sums use it on
all 1,312,000 draws (pool_row, runs16).

Thresholds, two levels (kaankolcu's idea, Provenance; model16.two_level). T = lower - 16 x slack, lower a certified lower bound of 16 log2 N_f (embedding_bounds.py: Taylor shift per place, outward rounding; NormBounds.block) for exactly the rows and cells it is used on: every B-smooth cell with omitted mass below slack is kept whatever the bound's width; a wide bound only sends more cells to the exact norms (charged at slack + 32 bits). Level 1: the rows of an aligned 8 x 8 x 8 box (model16.box_rows) share the hull (EB.hull) of the hulls of its 2 x 2 x 4 sub-boxes; in an aligned 256-cell block each 64-cell node meeting the box's cells is evaluated on it; certified (finite, <= 32 bits wide) it is a leaf for all rows, else halved, down to 16 cells. Level 2: a failed 16-cell node runs c24's EB.box_leaves on each sub-box (hull, bisection to single cells, then single rows; uncertified cells are rejected: none in the pools). Merge: with no level 2 in the block and (largest upper - smallest lower) <= 32 bits, the block is one leaf at the smallest lower bound, a lower bound for each of its cells. A row's part of a leaf is one scan leaf at 256 cells, else cut at multiples of 64 (the measured scan lengths). Counts are exact per pool row: 1.6e-04 evaluations, 0.0106 scan leaves per cell. Prices: C_NB_EVAL 1,178, C_NB_ROW 378 (measured, KERNELS.md, [10s]). NormBounds.row: per slice and place c_j = mid e_j, w_j = up((rad e_j + |c_j| 2^-49) UP);
a row takes m = sum b_j c_j, r = (sum |b_j| w_j)(1 + 2^-48): nine roundings < 2^-49.8 relative, so |b(v) - m| <= r. [2b] checks exact f norms of six anchor blocks and two level-2 blocks against their leaves, [2] and the audit the anchors' fields.

Record walk per canonical a-run (rec7; KERNELS.md). A canonical a's run is the maximal integer t-interval
containing 0 with a + t e_0 canonical and inside the limit s_c theta A^2 of a's class c; rec7 walks its sub-run of class c (runs16.py). On a run, a prime ideal (p, s, r)'s hits in all its slices form one coset of a 7-dim
lattice (t = r b(s) - a(s) mod p), walked once per run and ideal pair (p, s, +-r) with one setup6 ([10r], KERNELS.md).
Clean-run rule: rec7 walks a run only if all its a-values are in the canonical list (tested in the 2^92
row); otherwise its slices take the slice walk. A clean run of 32 to 2,048 a-values uses the run walk (REC7_K x the reference curve, KERNELS.md) in each band
p in [2^48, 2^64) with h' = 2 x run cells x lower/upper ([3]) / p >= 2^4.99 (the curve's bottom) instead of the slice walk and empty cosets; run
cells are exact (runs16.lrs: integer sum over the run's a-values of (lam_c - qa)^3, halved). Each row pays its cells' share
(per hit at h' incl. the walk's emission; setup6 per (run, pair) replaces its slices', f ideals). Contamination. A contaminated run's slices cost at most
the slice walk: per slice <= D_SUP = contam_sup = 2^72.639 on [2^48, 2^64) (exact band maxima with empty cosets and setup6;
[k7c]: slices peak at 2^72.498). Every row pays RUN_ALLOW x D_SUP/2 per slice (f ideals), RUN_ALLOW = 0.00091 (draws are
uniform over canonical a): in 4,000 a-uniform runs (PREREG-C: seeds 5001-5010, the first 400 canonical a inside their class limit
per seed, runs at that limit, sha256-16 in pool_d) the exact predicate found no non-canonical a-value among 1,920,263
(248,937 new, the rest in runs inside their run tested at 2565/8000; its 41 runs above 1,024 a-values hold 43,759):
contaminated runs hold under 1 - 0.05^(1/4000) = 0.00074866 <= 0.00091 of the canonical a (95%). [k7c]/[2r] replay 6 sampled/3 anchor runs.
Coset sieve at P2. 2 = P2^2 P4^2, P2 = wO; b is in P2 iff P2F.b is even (common12.P2F; [10r] checks wO = ker P2F).
If a and b are in P2 the pair is not primitive, so those cells are not sieved; P2 is fixed by units, sigma, tau and
b -> -b, so the canonical domain is unchanged. In the basis e_0, e_j - (P2F.e_j) e_0 a slice with a in P2 sieves
b_0 = 2k + 1 (half the cells, congruence constants halved mod p, setup6 + 10). P2F(a + t) alternates along a run, so
rec7 walks it as two cosets (1/2 and 1/4 of its cells, setup6 + 10 each). model16.evaluate: p2a-and-p2b rows carry
nothing; p2a slices are priced at S/2 or the sub-run's cells.
S' uses the lower count factor ([3]), so the per-a charge is >= 0.00091 D_SUP >= 0.000749 D_SUP.
runs16's run ends (< 2^69.3 units) are misc.

Yield ([4], pools G + H). rels = (va_pool/1968000) vb (0.998097) sum over G + H's canonical primitive region rows of
w rho(u_f) rho(u_g) r_f 0.999^lambda_f x the f band margin (g exact): Dickman rho ([1]), retentions r from the
two-stage law above (threshold_model.py), the ECM miss; cells with qa < 0.02 are sieved and charged, relations discarded. Estimate
2^63.16227 (SE 0.0384 bits) less the band charge 0.009: 2^63.1532 against 1.05 U' <= 2^63.01108: margin >= 0.1421
bits = 3.7 pool SE (>= max(0.06, 2 SE)). Band charge (H1'): [4b], above. Descent input. LOG2_RELS = 63.074 is <= the charged count
lowered by two relative pool SE, 2^63.07446 ([4] rels_input); the descent's mean norms 738.2/498.2 truncate G + H's
738.253/498.267 ([4]).
Exact 2-adic part (H1, [14v], [14t]). rho prices a random integer (mean v_2 = 1); each row's v_s = v_2 of its exact norm
is known (pool_row), so u_s = (n_s + alpha_s + E_s + kappa_s - v_s)/69.5 (model16.a2, also q_g): E = (2, 2.399985) >= the
local means (runs16; alpha_2 = 1 - E), kappa = (0.06, 0.05) >= kappa2(u) (sum_k 2^-(k+1) rho(u - (k-1-kappa)/69.5) = rho(u),
rising in u) at n_f 840, n_g 530, above every yielding row. Cells with a in P2 carry exactly 2^6 in N_f, the rest 2^0.
Toy census (92,108 primitive (a, b) in {-1,0,1}^12, exact 2^20-smooth N_f): class chi^2 at 2^18/20/22 10/6/18 (exact) vs
88/30/23 (average), totals within 1 SE.

Cost ([5]). Sampled mean 25.48 instructions/kept cell over pools S, G and W, charged 25.63 with the allowance, over
2^123.903 kept cells (every step at KERNELS.md/STORAGE.md constants; rec7 at its measured stratum means + 2 SE): collection <= 2^114.292681000
units (estimate 2^114.284308599 times allowance 1.00582). Allowance (H2, [4a]): a Monte Carlo estimate over S, G and W's
5,248,000 draws, relative SE row 0.1851%, incidence sqrt((1-192168/5248000)/192168) = 0.2239%, total <= 0.2906%;
charged 1.00582 >= 1 + 2 SE.

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U', filtsim.PEXC), clique removal and Markowitz merge (order: (rows - 2) x lightest unit pivot)
are simulated on the production occurrence profile (filtsim.py; [6] reruns 2^12) as the tau-quotient system: tau^k(I) enters I's orbit column with
coefficient chi^k (entries a + b chi, merged |a|, |b| <= 255). Cap. All retained relations are used; H1 puts the count at >= 2^63.0744
>= 1.095 U' ([4]): filtering at excess 1.095, set after G + H's count (sizing and judgment stay at 1.05; margin over
1.095 U' 0.0816 >= max(0.06, 2 SE), [4]). Collection's ECM curves are fixed by a public seed, so
relations, filter and N are deterministic in the challenge-independent instance: H3 puts
N(250)/U' at <= 0.02566; the attack aborts if N > 0.02568 U' (never under H3; P_CAP = 1e-9 charged); LA is
charged at 0.02568. [10c]: per-run sd x sqrt(U) <= 0.518 at every scale with >= 6 seeds. A campaign at B = 2^69.5,
means 744/501 (FILT_RUNS at 1.05, 6-16 seeds per scale 2^12..2^20; at excess 1.095 16 seeds at 2^19, 8 at 2^20) measures
x = N(250)/U' ([10] lists means and SE). H3: production is not above the largest scales: the
larger mean + 2 SE at 1.095 of 2^19 and 2^20 is 0.027154. That campaign's relations are larger than production's: at means 738.2/498.2 the same seeds give a pooled paired
ratio 0.9294 +- 0.0039 (2^16 0.9304 +- 0.0057, 2^17 0.9285 +- 0.0053; FILT_Y, [10y]). It is charged at + 2 SE, 0.9372, for the share 0.9759 of the way from 744/501 that G + H's means 738.253/498.267 cover (linear, worse side): 0.027154 x 0.938643 + 0.000165 (Scale) <= 0.02566 (H3: the ratio does not grow with scale). Merge to mean weight 250 (charged 251; set from 2^16/2^17
before 2^19/2^20). N <= 0.02568 U' <= 2^57.6575 < 4 x 2^64/320: four owners (STORAGE.md).
Scale ([6m]). C(U, l) keeps primes below 2^l as columns, l = min(log2(R/25)+1, 0.8L): 9-16.5 bits at 2^12..2^20 (no run changes), defined at every U. The ideals of [2^l, 2^0.8L) are counters: weight r[-1], count r[-2], coefficient bound imul <= 255, price classes (unit_entries). Production merges by C's rule on these counters and never pivots a sub-band ideal; a sub-band entry is the sum of its counted contributions (relations with a squared sub-band ideal, < 0.02%, are dropped: 0.0003 bit of the 0.0816 count margin). npt and the entry price (U5 a, a <= 4, else G1; 4 U5 <= G1 <= 5 U5) are monotone and subadditive, so a row's true weight and price are at most its counted ones, and each sub-band column needs one more row. So at U', for every l in [16.4, 0.8L]: N <= N(C(U', l)) + pi(2^0.8L) <= N(C(U', l)) + 0.000165 U' (Dusart) and k5 <= k5(C(U', l)). Both bounds hold at any scale. H3 and H4u extrapolate C itself: explicit share 2.4% -> 0.8% and N(250)/U 0.054 -> 0.028 at 2^12..2^20, k5 8.41 -> 8.35 at 2^15..2^17.

The square system A = [M | S] has N + 10 rows, N orbit columns and 10 dense Schirokauer columns: the chi-eigenpart of
the unit groups, (rk K_f - rk K_f^+ + rk K_g - rk K_h)/2 = (35 - 17 + 5 - 3)/2 from signatures K_f (0, 36), K_f^+ (0, 18),
K_g (0, 6), K_h (2, 2) ([13]; tau^2 = -1 on the sigma-odd part, so its chi and chi^-1 eigenspaces have equal dimension);
under H3 its kernel is the line of virtual logs (w, c), c != 0, and M has rank N. Block Wiedemann (m, n) = (8192, 130) with S as right-hand sides
(as CADO-NFS): B = [M | 0]; Y = (S_1..S_10, B z_11..B z_n), z uniform in F_r; X block-diagonal (H-X): the 10 pad coordinates (zero columns of B)
in blocks 1..10, the others shuffled and cut into m blocks, X_l uniform in [0, 2^128) on block l. A generator sum_k B^k Y f_k = 0 of
a_i = X^T B^i Y with f_0 = (c, e), c != 0, gives S c + B w = 0, w = Z e + sum_(k >= 1) B^(k-1) Y f_k: one Horner pass, no product touches S
([10x]: toy end to end). Krylov N/m + N/n block steps; Horner runs as five interleaved chains from checkpoints
B^(jL) Y (L = Horner length / 5; 4n = 520 checkpoint vectors, < 2^72 bytes, 4L < N/m + N/n) at the five-vector rate. Each
vector element is stored with chi x_j. Rows reduce exactly (_row_red). Per row and vector (STORAGE.md, [10]): 251 x 8.39 (H4u, PREREG-U) + 662 (row tail, four owners) + 200 (_chi_mul); per
Horner element n generator terms and one reduction (270); projections: each vector element meets one
X entry, (1 + n/m) N^2 terms at 60. Linear algebra <= 2^112.670991256 units; generator (matrix Berlekamp-Massey by half-gcd, k = m + n = 8322: a node of l
sequence terms makes (m + k) k^2 products of polynomials with <= l and <= l/2 + 1 coefficients, each one Kronecker integer product in 9-limb
slots (coefficients < k L r^2 < 2^576): F(18 l + 18) under H-FFT, + 640 l to pack, reduce (270) and add; leaves m k^2 (97 + 270). The sum over
the recursion is below model16.generator ([10m]): 4 k^3 log2 L (F(18 L) + 640 L) instructions, L = (N + 10)(1/m + 1/n) + 2): <= 2^101.943 units.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent; with seeded ECM curves the relation count is fixed:
2^63.1532 under H1 (section 3). Under H3 (section 9; the 10 Schirokauer maps of the chi-eigenpart
injective on units modulo r-th powers) Block Wiedemann fails with probability <= 4N/2^128 < 2^-66 (Schwartz-Zippel: a minor of order N + 9 of
(X, B^T X, ..)^T (Y, B Y, ..), degree <= 2N + 18 in the entries of X and z, uniform on >= 2^128 values; H-X); a failed solve (checked on the dropped rows) is a failure.

Reconstruction. Purged ideals' logs come by back-substitution; ideals with <= 2 occurrences are log-less in the descent
([6b]: model >= 6.3x measured; DESCENT.md).

## 6. Individual logarithms

DESCENT.md ([7], [7e], [7l]): 62 nodes + boot (215, 304, 15, 4), retention >= 0.6490.

Expected two-log work <= 2^109.277 units ([7]; unrounded 2^109.27611; one log 2^122.564 instructions: boot own
2^120.08 (whole balls, DESCENT.md) + subtree 2^121.40 x retry 1.836). The meter caps both logs at t = 2^121.9 units (reached = failure). Metered work is min(D, t) <= D for the uncapped two-log work D: its expectation is <= 2^109.277 for any t, and
P(cap) <= E[D]/t <= 2^(109.277 - 121.9) <= 0.0001586 (Markov); with 10% more expected work, <= 0.0001744.

## 7. Success

Core: success >= (1 - 2^-66 - P_CAP) x (1 - 2 x 7.2e-07 - 0.0001586 - 2^-12.35 - 2^-10.70 - 2/r) >= 0.99904 (main
phase, two boots with all 18 relations failing, cap, node counter, repeats); >= 0.999031 + 2^-48 with the cap bound at +10% (0.0001744).

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX = 2^27.3. E[N] uses the same
  recurrence with unit node cost (counts16.py): E[N] <= 2^14.9427 for both logs, P(abort) <= 2^-12.35.
- Repeats: N_MAX nodes in the smallest bin (2^70.30 ideals), 2^6 for 1/p weighting: P <= 2^-10.70, charged.

Subtrees are independent (H5). Past the gate P(success) >= 0.999031 with the wrapper (the +10% figure); every abort and a closed gate fail: P(success) >= (K/2^64) x 0.999031 >= 0.39. W <= 1 + 1[u < K] W_run, u independent of W_run: E[W] <= 1 + (K/2^64) E[W_run] (audit: gate, expected_run).

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (section 3, H-EA; estimate 114.284308599 + 0.582% allowance) | <= 114.292681000 |
| linear algebra (weight 251, N <= 2^57.6575, (8192, 130), block-diagonal X; section 4) | <= 112.670991256 |
| matrix generator | <= 101.943 |
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
back-substitution 2^80 (2^24.4 per entry against < 2^14). The rows sum to <= 2^121.90977184100883: core meter 2^121.9097718411 (audit-checked); with the two logarithms at 2^109.277, not the cap, <= 2^114.7322330618.

Schirokauer maps. R = (Z/r^2)[y, x]/(h, f) (rank 72; g: rank 12); r is unramified ([13]): eps = lcm(r^d_i - 1) < r^72
(18,360 bits; g: 3,060). The 10 columns are fixed F_r-linear forms (chi-eigenpart) of (phi^eps - 1)/r mod r, phi = a - bx,
by left-to-right powering on eight-limb coefficients. An f squaring: 2,628 = 72 x 73/2 schoolbook 8 x 8-limb products
(<= 512 instructions), reduction by monic h, f <= 2^19, 72 Barrett (<= 1,024): < 2^20.9; times phi (12 coefficients
< 2^64) < 2^19.3. Per bit < 2^21.31: f < 2^35.5, g < 3,060 x 2^17, so < 2^36 per relation, 2^86.86 units at 2^65.15
relations, within 2^88.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (0.29, 0.75), exact v_2 and kappa, exact sampled norms; boot pair ratio kappa <= 2^24, DESCENT.md). H1' retention (Euler-product model, two-stage law; band charge 0.009 bits). H-EA (work only): H1''s survivor law prices the survivors of both stages ([12w]). H-FFT (work only): F(N) bounds the local SSA work past 2^26 limbs (KERNELS.md). H2 region (uniform pools, exact flags, measured counts ([3], 6 SE; certified bracket as fallback); sized on S, judged on G + H; S+G+W work with the 0.582% allowance). H3 matrix (r unramified, [13]; r prime to the class numbers, the unit condition, a one-dimensional kernel; N(250)/U'
<= 0.02566: neither measured factor grows with scale, section 4). H4 kernel cost (KERNELS.md, STORAGE.md); H4u (unit share >= 86%;
k5 <= 8.39, PREREG-U, [10]); rec7's mean per-hit price in each density stratum is at most the mean + 2 SE of its 512 production draws (REC7_K, KERNELS.md); contaminated runs hold at most 0.091%
of the canonical a (95% bound; a shortfall changes cost only); f/g independence given q_g (100x). H-X: BW's minor (section 5; N + 9 = dim(span S + Im B)) is not identically zero for this X and Y; m, n >= 10 ([10x]: toys solve end to end). The pad pinning is needed: x^T B^i (i >= 1) vanishes on ker B (the pad space), so only X itself separates it; with all coordinates shuffled the solve fails when the pad coordinates meet <= 8 blocks. H5 descent (DESCENT.md): lower sums; per-ideal yield certificate (DISC); lazy order independent of outcomes; negative-binomial recurrence, independent
subtrees of distinct ideals (measured: [7e]); missing logs (Poisson occurrence at the 2-SE input); local omitted-mass laws; random-integer ECM orders; UDUP = 0.65 (DESCENT.md, [11]).

## 10. Memory

Core storage <= 2^105 bytes (STORAGE.md), wrapper below 2^30: below 2^106 bytes. All preprocessing is for this target,
charged once in W.

## Provenance

GT core: kaankolcu recover-P 1b8fa720 and its lineage (9332a74d mpjunior92; 2c2ce2ab / 42e80f80 gbotrel; 6282ab79
bernardo-lima12; 6734bd5a GPT-6 in Codex; Gemini 45ce22e5 / 4bb6c381); polynomials after Guillevic-Singh; filtering after CADO-NFS. H-EA, H-X,
rec7, band strip and split, f-only sieve, g batch, twisted C4 field: gbotrel, Claude Opus 5.5 in Claude Code. Wrapper: 4f8baf8c. Run gate, expected total work, high cap: kaankolcu (8ab78e96, 5b97b072, b6f0a0f6, d7ca6b98); logarithms at their expectation: also mpjunior92 (ea928ab4; and blocks (8192, 130)); residue price 0.02, the same norm-kernel design, two-level thresholds: kaankolcu (538bf09c, 8b2b8442; 01aabb8f, a02dda7e). Two-stage detector: gbotrel, Claude Fable 5.2. Weight 250 and pool W: kaankolcu, Claude Opus 5.5.

# recover-Q (bls12-381) through two GT logarithms with the STNFS core: 2^116.322655031201 units

Track bls12-381-fixed-g1-invert-g2. Input: a nonidentity P in G1 and z in GT; output Q in G2 with e(P,Q) = z
for the exact supplied P (Q = infinity for z = 1). The score unit is one dense F_p12 multiplication, C = 20,000
scalar AArch64 instruction equivalents. No advice, smaller targets or amortization
are used. Predecessor work is reused as credited at the end; no prior verdict is a premise.

## Claim

The algorithm below returns such a Q with success at least 0.39 over uniform
challenges and fresh coins under the declared heuristics H1-H5, H-EA, H-X and H-FFT (section 9). All classical work (preprocessing, both logarithms, failures,
restarts, field map, pairing, scalar multiplication, verification, output) is at most

    W <= 2^116.3226550312 + 2^66 < 2^116.322655031201  dense-Fp12 units (C = 20,000).

The score is 116.322655031201 (core: sections 2-8; wrapper: section 1). Storage is below 2^106 bytes. Analytic, exploratory claim.

Conventions. Bounds round in the safe direction; evidence16.py checks each stated bound against its unrounded
value (declared_bounds) and reruns the checks cited [n]. Model:
config.json, stats16.b85, pool_d.json, the constants of model16.py and descent16.py (owner10.install: C_ROW
1,320).

## 1. Reduction (fixed P: two GT logarithms) and field map

With the pinned reduced pairing e and the public generator Q0 of G2:

1. Check P (on curve, [r]P = infinity, P != infinity) and z (z^r = 1). If z = 1 return Q = infinity.
2. g = e(P, Q0) (pinned pairing(Q0, P), protected encoding).
3. Map g and z into the core's tower (field-map.md).
4. Run the core on (g, z) under its meter and node counter (any failure is a failure): a = log_g(z) mod r.
5. Q = [a]Q0; output Q if e(P, Q) = z, otherwise fail.

GT has prime order r and e is nondegenerate, so for the exact supplied P g generates GT, z = g^a for a unique a mod r and
e(P, [a]Q0) = z. For uniform P, g is uniform on GT minus {1}, independent of the uniform z: the core's input law. The
field isomorphism preserves z = g^a; the final check rejects any wrong output.

One relation/LA phase serves log z and log g (section 6; 2/r charged).

Wrapper below 2^80 instructions (pairing, scalar multiplication and check below 2^50, field map below 2^67 by
field-map.md): 2^80 / 20,000 = 2^65.71 units, a 2^66 row. Wrapper failure: field map below 2^-48.

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

Candidates are (a, b) in O_h^2 with a - bx B-smooth on both sides, B = 2^69; the factor base is the relative-degree-1
prime ideals over K_h of norm <= B on both sides, U' = 2^62.451 tau-orbits (U/2: orbits of 4 ideals; U = 2^63.451).

Region. Exact rational forms Q_a, Q_b (config.json), A^2 = 2^21. Pool (common12.regenerate_pair): uniform on theta = 3/8 (a5 > 0, Q_a(a) <=
(3/8) A^2; b uniform in the half-slice b5 >= 1, Q_b(b) <= A^2 - Q_a(a)). Region: a5 > 0, Q_a(a) <= (2709/8000) A^2; b5
>= 1, Q_b(b) <= (17/20) A^2 - Q_a(a) (invariant under b -> -b and the sign fold; exact Fractions), so the restricted sample is uniform on it. G has 32531 region rows (22646
primitive) of 87197 canonical draws, scaled by the region volume over the effective sample 656000 x (2709/3000)^3 =
483022.199. Two weights, never mixed: density w = (1 - qa)^3 (yield) and vb (1 - qa)^3 x factor (work, population);
slice-geometry costs use S' = vb (17/20 - qa)^3 x factor and K' = floor(sqrt((17/20 - qa) A^2 / Q_b(e_5))) lines. The
canonical a is the (Q_a, lex)-least element with a5 > 0 of its <+-1, units, sigma> class (exact enumeration of tO, t in {a, sigma a}); with b5 >= 1, one point per <+-1, units, tau> orbit (tau^2(a, b) = (a, -b)); non-primitive pairs are excluded by an exact HNF test. lattice_count_bounds.py certifies (exact rationals) the count ratio at every radius R >= 842 = isqrt((2709/8000) 2^21), below the smallest b-slice
radius 1035 ([3]): [0.977338, 1.022363] (Q_a), [0.989198, 1.009951] (Q_b); the model uses the products 0.966780570
(yield) and 1.032535599 (work). 

Pools (PREREG-G). 656,000 uniform pairs each on theta = 3/8 (164,000 per seed): S (seeds 261-264) sizes the
region, G (251-254, drawn after the configuration was fixed) is judged, both give the work; no draw dropped. pool_d.json holds per-seed counts, list
hashes and sums, 8 anchor rows of G, G's canonical lists over draws 0..49 ([2] replays them) and the contamination sample. theta is the smallest t on the 1/8000 grid with margin_S
>= 0.06 + sqrt(2) SE_S (pre-registered), on S: 2709/8000 (0.095785 >= 0.095601; 2708: 0.095344 < 0.095611; [4s]). The band
charge applies to G's count only.

Sieve and survivors (f ideals only). Line sieve (scores only) for p in [2^14, 2^48); record walk for p in [2^48, 2^64) (slice12,
KERNELS.md; X2 = 2^64; per clean a-run, below; [2^63, 2^64) of contaminated runs on the wide path); primes in
[2^64, 2^69] by ECM within slack 302. Survivors pass the f threshold; exact norms (C_NORMS) give the g norm z for the g batch. Primes below 2^14 are the retention models' omitted mass (prime_min
2^14, explicit_bits 14). Survivors with z 2^69-smooth get trial division by every prime below 2^20 on both cofactors (82,025 primes, 2^23.649 <= C_TDIV = 2^24.5, [5t]). Batch band strip ([10b], KERNELS.md). Survivor sides wait in
waves of K = 2^39 (13-limb cofactors, <= 832 bits). With P the product of the primes in (2^20, 2^48] (< 1.000028 x
2^48/ln 2 bits, Rosser-Schoenfeld; built once), a scaled remainder tree (Bernstein; model16.srt)
gives P mod z for every z, and the band part is gcd(z, (P mod z)^(2^6) mod z) (2^6 >= 41 >= the exponent of any band
prime in z); the f band-miss factor (1 - 41 x 0.001^2) stays as a margin.
g batch ([10g]). Same SRT with P_g the product of all primes <= 2^69 (< 1.000028 x 2^69/ln 2 bits), waves of K = 2^61 g norms (<= 576 bits, 9 limbs): z is 2^69-smooth iff (P_g mod z)^(2^10) = 0 mod z (exponents <= 576 < 2^10). r_g >= 1 - 2^-38 (ideals of norm p^d > 2^69 drop a relation at the split; in the margin). A kept relation's g norm is split by trial division, the band split and ECM (<= 12 primes of (2^48, 2^69], until complete): <= 12 x 1.002 x 2^39.79 + 13 C_MR = 2^43.38 per relation, <= 2^65.15 relations ([10t]): 2^95 units (ledger row).
Lemma (exact SRT). Product tree Z_v (L <= 61 levels); node v holds Y_v of b_v = 64 (n_v + G) bits (n_v limbs of
Z_v, G = L + 2). Root: Newton quotient q of P 2^b by Z, corrected by <= 2 add/subtracts, Y_root = q mod 2^b. Child c,
sibling s: Y_c = floor((Y_v Z_s mod 2^(b_v))/2^(b_v - b_c)). Leaf: r = ceil(z Y 2^-b) (0 if r = z) = P mod z. Proof:
with y_v = frac(P/Z_v) and eta_v = (y_v - Y_v 2^(-b_v) mod 1) Z_v, eta_root < 2^(-64G) and each level adds < 2^(-64G),
so all eta < (L + 1) 2^(-64G) < 1 (no wrap); at a leaf z Y 2^-b = (P mod z) - eta. Cost (H-FFT, KERNELS.md): an N-limb product <= F(N) (model16.fft_F);
per level products of 2 n_(j-1) and n_j + G + n_(j-1) limbs + 3 per copied limb; root <= 6 F(2(max(|P|, |Z|) + G)); leaf <= 2^15
(g) / 2^19 (band) (model16.gbatch_cost): C_BATCH = 2^24 (2^23.921, K = 2^39), C_GBATCH = 2^25 (2^24.973, K = 2^61). [10v]: model16.srt on 2^8
planted leaves, Newton offsets -2/+2, against Python's %.
P's build (< 2^70 instructions) is in misc, P_g's tree (< 2^88) is the 2^74 row; trees < 2^74 bytes.
Band split ([10t]). The strip decides smoothness only. A kept relation (both sides smooth, primitive, qa >= 0.02)
is split: the band parts are trial-divided by every prime of (2^20, 2^48] (< 2^43.272; 180 instructions per prime with a 2^(64 i) mod p table, misc row), each hit divided out to full power (<= 41 + 28 divisions, 2^10 each). Every found prime
gets its ideal (p, s, r) and exponent from the factor-base entries above p , <= 2^26: 2^50.7634 instructions per relation, deterministic, no relation lost. Kept
relations are <= the B-smooth region cells: G's 2^62.65016 + 2 SE over the yield envelope
(1.032535599/0.966780570), the ECM miss (0.999^22: <= 13 + 9 primes above 2^64), the band margin and the least
retention on the model's norm grid (0.7775 f, 0.5393 g): <= 2^64.1140, 2^100.5897 <= 2^100.59 units (ledger row;
hosted audit replays it).
Survivors' f omitted mass is bounded by 310 bits; ECM is the ecm16.py schedule at 99.9% per
prime ([12]): 69-bit primes in a 310-bit cofactor (B1 = 2^16, B2 = 30 B1, wheel 30030, 403 curves, five-limb Montgomery
multiplication at 367.76 Ir) cost 2^38.0909 instructions.
Cofactor windows (H-EA, model16.Window). After trial division, records and band strip a smooth side's cofactor C
(< 2^320) has no prime below 2^64, so it is smooth only if C = 1, C <= 2^69 is prime, or C is composite with
64k < log2 C <= 69k, k in {2, 3, 4}; ECM runs only in the last case. Compositeness: 40 Miller-Rabin rounds, random bases (C_MR = 40 x 2 x 310 five-limb
products = 2^23.121); a composite passes with probability <= 4^-40 (Rabin-Monier), below the
count's rounding. ECM is charged at e x its schedule, e = the survivor law's mass whose part above 2^64 lies in a
window widened by 1 bit, minus its single-prime term, over the survivor mass (e_f 0.1553-0.1569 at
the anchors). g first: the f side's division, strip, test and ECM run only if the g batch passes:
per f-survivor C_NORMS + C_GBATCH + q_g (C_TDIV + C_BATCH + C_MR + e_f 2^38.0909), q_g = min(1, 100 rho(u_g)) per pool row
(H1, factor 100 a margin); f-threshold survivors per cell (~2^-22.5). S and G work sums use this
model on all 1,312,000 regenerated draws (pool_row, runs16).

Thresholds. Certified per-block f bounds (embedding_bounds.py, 64-cell blocks, bisected to within 8 bits)
give T = lower - 16 x slack, keeping every B-smooth cell with omitted mass below slack (Taylor shift per place, outward
rounding); one evaluation (NormBounds.block) is what nbound31 replays bit for bit at 1,154.6 Ir. A line's rows
in an aligned 2 x 2 x 4 box (model16.box_rows) share one bisection per block over their hull. Each cell pays C_NB_EVAL x its box's evaluations per box cell, each row
its leaves, context and hull share (C_NB_ROW 4,600). [2b] checks 768 exact f norms of six boxes, [2] the anchors' fields.
Warm start (clean runs). A box-block's bisection starts from the leaves that the run's first slice holding it
certified from the aligned 64-cell root (model16.pool_row_warm): same tests and fallback, each leaf certified for its
own a. Charged: warm evaluations (+ RUN_ALLOW x the cold excess), max(cold, warm) leaves; [2] and the audit replay the anchors.

Record walk per canonical a-run (rec7; KERNELS.md). A canonical a's run is the maximal integer t-interval
containing 0 with a + t e_0 in the region and canonical (runs16.py). On a run, a prime ideal (p, s, r)'s hits in all its slices form one coset of a 7-dim
lattice (t = r b(s) - a(s) mod p), walked once per (run, ideal) with one setup6 ([10r]; bases, checks and campaign in
KERNELS.md).
Clean-run rule: rec7 walks a run only if every a-value of it is in the canonical list (whose exact test the 2^92
row pays); otherwise its slices take the slice walk. A clean run of 32 to 1,024 a-values uses the run walk in each band
p in [2^48, 2^64) with h' = run cells / p >= 2^4.99 (the table's bottom) instead of the slice walk and empty cosets; run
cells are exact (runs16.lrs: integer sum over the run's a-values of (17/20 - qa)^3). Each row pays its cells' share
(per hit at h', + 2 per record; setup6 per (run, ideal) in place of its slices' setup6, f ideals). Contamination. A contaminated run's slices cost at most
the slice walk: per slice <= D_SUP = contam_sup = 2^72.591 on [2^48, 2^64) (exact band maxima with empty cosets and setup6;
[k7c]: 2001 region slices peak at 2^72.445). Every row pays RUN_ALLOW x D_SUP/2 per slice (f ideals), RUN_ALLOW = 0.00091: pool draws are
uniform over canonical a; in 4,000 a-uniform runs (PREREG-C: seeds 5001-5010, the first 400 canonical a with Q_a <=
(2709/8000) A^2 per seed, sha256-16 in pool_d.json) the exact predicate found no non-canonical a-value among 2,015,033:
contaminated runs hold under 1 - 0.05^(1/4000) = 0.00074866 <= 0.00091 of the canonical a (95%). Work sums add it on every S and G region row. [k7c] replays 6 sampled runs, [2r] 3 anchor runs.
Coset sieve at P2. 2 = P2^2 P4^2, P2 = wO; b is in P2 iff P2F.b is even (common12.P2F; [10r] checks wO = ker P2F).
If a and b are in P2 the pair is not primitive, so those cells are not sieved; P2 is fixed by units, sigma, tau and
b -> -b, so the canonical domain is unchanged. In the basis e_0, e_j - (P2F.e_j) e_0 a slice with a in P2 sieves
b_0 = 2k + 1 (half the cells, congruence constants halved mod p, setup6 + 10). P2F(a + t) alternates along a run, so
rec7 walks it as two cosets (1/2 and 1/4 of its cells, setup6 + 10 each). model16.evaluate: rows with p2a and p2b
carry nothing; with p2a the slice has S/2 cells, per-cell costs are per kept cell and prices are read at S/2 or the
sub-run's cells. Relations are unchanged.
Count envelopes ([3]) leave the per-a charge >= 0.989198/1.022363^2 x RUN_ALLOW x D_SUP = 0.000861 D_SUP >= 0.000749 D_SUP.
The rule's own work (runs16 ends, ~2^20 instructions per canonical a, < 2^69.3 units) is in misc.

Yield ([4], pool G). rels = (va_region/483022.199) vb (0.966780570) sum over G's canonical primitive region rows of
w rho(u_f) rho(u_g) r_f 0.999^lambda_f x the f band margin (g exact): Dickman rho ([1]), retentions r from the
Euler-product threshold model (threshold_model.py, calibrated on exact smooth integers at B = 2^10..2^14; prime_min
2^14), the ECM miss; cells with qa < 0.02 are sieved and charged, their relations discarded. Estimate
2^62.65016 (SE 0.0416 bits) less the band charge 0.009: 2^62.6412 against 1.05 U' <= 2^62.52173: margin >= 0.1194
bits = 2.9 pool SE (feasibility requires max(0.06, 2 SE)). Band charge (H1'). On 4,000 random region rows
(exact norms) the unsieved [2^8, 2^14) content is f 5.451 +- 0.121 bits (model 5.922): no excess; REL_BAND_LOSS = 0.009 bits is deducted in [4] as a margin. Descent input. LOG2_RELS = 62.555 is <= the charged count
lowered by two relative pool SE, 2^62.55545 ([4] rels_input); the mean norms 741.2/498.6 truncate G's
741.26/498.66 ([10y], H3).

Cost ([5]). Sampled mean 54.29 instructions/kept cell over pools S and G, charged 54.81 with the allowance, over
2^124.412 kept cells (bounds, layout, setups, line path, record walk or rec7 with contamination, routing, apply, recovery,
scan2, band strip and cofactorisation, at the constants of KERNELS.md and STORAGE.md): collection <= 2^115.900562000
units (estimate 2^115.886949166 times allowance 1.00948). Allowance (H2, [4a]). The work is a Monte Carlo estimate of
a fixed, challenge-independent total over the 1,312,000 draws of S and G; its relative SE has a row term 0.2809% and an
incidence term sqrt((1-65260/1312000)/65260) = 0.3816%, total <= 0.4739%; the charge multiplies the estimate by
1.00948 >= 1 + 2.00 x 0.004739 = 1.009478 (rounded up; [4a] fails below), a modeling margin, not a probability.

## 4. Filtering and linear algebra

Purge (to excess 0.0001 U', filtsim.PEXC), clique removal and Markowitz merge are simulated on the production
occurrence profile (filtsim.py; [6] reruns 2^12) as the tau-quotient system: tau^k(I) enters I's orbit column with
coefficient chi^k (entries a + b chi, merged |a|, |b| <= 255). Cap. All retained relations are used; H1 puts the count at >= 2^62.5555
>= 1.07 U' ([4]): filtering at excess 1.07, set after G's count (sizing and G's judgment stay at 1.05; G's margin over
1.07 U' is 0.0922 >= max(0.06, 2 SE), [4]). Collection's ECM curves are fixed by a public seed, so
relations, filter and N are deterministic in the challenge-independent instance: H3 puts
N(300)/U' at <= 0.03459; the attack aborts if N > 0.03461 U' (never under H3; P_CAP = 1e-9 charged); LA is
charged at 0.03461. [10c]: per-run sd x sqrt(U) <= 0.747 at every scale with >= 6 seeds. A campaign at
means 744/501 (FILT_RUNS at 1.05, 6-16 seeds per scale 2^12..2^20; at excess 1.07 16 seeds at 2^19, 8 at 2^20) measures
x = N(300)/U' ([10] lists means and SE). H3: production is not above the largest scales; the pre-registered charge is the
larger mean + 2 SE at 1.07 of 2^19 and 2^20 -> 0.03459. At G's means 741.2/498.6 the same seeds give a pooled paired
ratio 0.9460 +- 0.0037 (+ 2 SE 0.9533 <= 1) (2^16, 2^17; FILT_Y, [10y]). Merge to mean weight 300 (charged 301).
N <= 0.03461 U' <= 2^57.5987 < 7 x 2^55: seven owners (STORAGE.md).

The square system A = [M | S] has N + 10 rows, N orbit columns and 10 dense Schirokauer columns: the chi-eigenpart of
the unit groups, (rk K_f - rk K_f^+ + rk K_g - rk K_h)/2 = (35 - 17 + 5 - 3)/2 from signatures K_f (0, 36), K_f^+ (0, 18),
K_g (0, 6), K_h (2, 2) ([13]; tau^2 = -1 on the sigma-odd part, so its chi and chi^-1 eigenspaces have equal dimension);
under H3 its kernel is the line of virtual logs. Block Wiedemann (m, n) = (72, 5): Y uniform in F_r; X block-diagonal
(H-X): shuffled coordinates cut into m blocks, X_l uniform in [0, 2^128) on block l. Krylov N/m + N/n block steps; Horner runs as five interleaved chains from checkpoints
A^(jL) Y (L = Horner length / 5; 20 checkpoint vectors, < 2^68 bytes, 4L < N/m + N/n) at the five-vector rate. Each
vector element is stored with chi x_j. Per row and vector (STORAGE.md, [10]): 301 x 9.87 (H4u, PREREG-U; an entry a + b chi
is |a| + |b| unit entries on x or chi x) + 1,320 (row tail 1,317, seven owners) + 10 x 97 + 517 (chi x: one bw_full
term and reduce9); per Horner element n dense terms and one reduction (420); projections: each vector element meets one
X entry, (1 + n/m) N^2 terms at 60. Linear algebra <= 2^113.815127389 units; generator (matrix Berlekamp-Massey, recursive
doubling, Karatsuba polynomial-matrix products, k = m + n = 77; model16.generator: 4 k^3 (97 + 420) L^(log2 3) log2 L
instructions, L = (N + 10)(1/m + 1/n) + 2): <= 2^109.083 units.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent: the relation count is fixed (seeded ECM curves):
2^62.6412 under H1 (section 3). Under H3 (section 9; the 10 Schirokauer maps of the chi-eigenpart
injective on units modulo r-th powers) Block Wiedemann fails with probability <= 4N/2^128 < 2^-66 (Schwartz-Zippel,
degree <= 2N, entries uniform on >= 2^128 values; H-X); a failed solve (checked on the dropped rows) is a failure.

Reconstruction. Purged ideals' logs come by back-substitution; the descent treats every ideal with <= 2
occurrences as log-less ([6b]: the model's count >= 6.3x the measured one; DESCENT.md).

## 6. Individual logarithms

DESCENT.md ([7], [7e], [7l]): 62 nodes + boot (215, 328, 16, 4), retention >= 0.6490.

Expected two-log work <= 2^111.651 units ([7]; unrounded 2^111.65065; one log 2^124.938 instructions: boot own
2^117.75 + subtree 2^124.05 x retry 1.838). The meter caps both logs at 2^112.51 units (reached = failure):
P(cap) <= E[W]/cap <= 2^(111.651 - 112.51) <= 0.5513346 (Markov); with 10% more expected work, <= 0.6064681.

## 7. Success

Core: success >= (1 - 2^-66 - P_CAP) x (1 - 2 x 7.29e-07 - 0.5513346 - 2^-12.00 - 2^-10.70 - 2/r) >= 0.4478 (main
phase, two boots with all 18 relations failing, Markov, node counter, repeats); >= 0.39268 with Markov at +10% (0.6064681).

- Node counter: the run aborts (failure) before special-q node N_MAX + 1, N_MAX = 2^27.3. E[N] uses the same
  recurrence with unit node cost (counts16.py): E[N] <= 2^15.2904 for both logs, P(abort) <= 2^-12.00.
- Repeats: N_MAX nodes in the smallest bin (2^70.30 ideals), 2^6 for 1/p weighting: P <= 2^-10.70, charged.

Subtrees are independent (H5). With the wrapper: P(success) >= 0.4478 - 2^-48 > 0.39; every abort is a failure.

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (section 3, H-EA; estimate 115.886949166 + 0.948% allowance) | <= 115.900562000 |
| linear algebra (weight 301, N <= 2^57.5987, (72, 5), block-diagonal X; section 4) | <= 113.815127389 |
| matrix generator | <= 109.083 |
| both individual logarithms (metered cap; expected 111.651) | 112.51 |
| canonical a list | 92 |
| Schirokauer maps | 88 |
| filtering streams, misc | 80 |
| factor base, roots, ideals | 73 |
| wide factor-base lattices ([2^63, 2^64), A and coset, KERNELS.md) | 71.6 |
| rec7 bases (2 per ideal x run and two cosets, p in [2^48, 2^64), KERNELS.md) | 75.09 |
| band split (section 3, [10t]) | 100.59 |
| P_g product tree (section 3) | 74 |
| g split of kept relations (section 3) | 95 |
| numerical certificate construction | 50 |
| verification | 11 |
| wrapper (section 1; <= 2^80 instructions) | 66 |

Fixed rows (units x 2^14.29 = instructions): canonical list 2^92 (2^42.8 per region a, 2^63.221 of them); Schirokauer maps 2^88 (below); factor base 2^73 (2^23.4 per prime, pi(2^69) < 2^63.46); filtering and
back-substitution 2^80 (2^24.4 per entry against < 2^14). The rows sum to <= 2^116.32265503114062: core meter 2^116.3226550312 (audit-checked).

Schirokauer maps. R = (Z/r^2)[y, x]/(h, f) (rank 72; g side (h, g), rank 12); r is unramified
([13]), so eps = lcm(r^d_i - 1) over r's residue degrees d_i is < r^72 (18,360 bits; g: 3,060). The 10 columns are fixed
F_r-linear forms (chi-eigenpart) of the 84 coordinates of (phi^eps - 1)/r mod r, phi = a - bx: left-to-right powering,
coefficients mod r^2 < 2^510 in eight 64-bit limbs. An f squaring: 2,628 = 72 x 73/2 schoolbook 8 x 8-limb products
(<= 512 instructions), reduction by monic h, f <= 2^19, 72 Barrett (<= 1,024): < 2^20.9; times phi (12 coefficients
< 2^64) < 2^19.3. Per bit < 2^21.31: f < 2^35.5, g < 3,060 x 2^17, so < 2^36 per relation, 2^86.86 units at 2^65.15
relations, within 2^88.

## 9. Heuristics

H1 smoothness (Dickman/Murphy, alpha (0.29, 0.75), exact sampled norms). H1' retention (Euler-product model; band charge 0.009 bits). H-EA (work only): H1''s survivor law prices the ECM share ([12w]). H-FFT (work only): F(N) bounds the local SSA work past the measured 2^26 limbs (KERNELS.md, 2.76x headroom). H2 region (uniform pools, exact flags, certified counts; sized on S, judged on G; S+G work with the 0.948% allowance). H3 matrix (r unramified, [13]; r prime to the class numbers, the unit condition, a one-dimensional kernel; N(300)/U'
<= 0.03459, section 4). H4 kernel cost (KERNELS.md, STORAGE.md: static counts of loop-free code, executed counts for routing); H4u (unit share >= 86%;
k5 <= 9.87, PREREG-U, [10]); rec7's per-hit price at production stays below its h' bucket's measured maximum (lattice shapes vary; p >= 2^47.03 measured; runs of 32 to 1,024 a-values); contaminated runs hold at most 0.091%
of the canonical a (95% bound; a shortfall changes cost only); f/g independence given q_g (100x). H-X: BW's nondegeneracy polynomial is nonzero on shuffled block-diagonal X ([10x] toys reach full rank). H5 descent (DESCENT.md): lower sums; lazy order independent of outcomes; negative-binomial recurrence, independent
subtrees of distinct ideals (measured: [7e]); missing logs (Poisson occurrence at the 2-SE input); local omitted-mass laws; random-integer ECM orders; UDUP = 0.65 (DESCENT.md, [11]).

## 10. Memory

Core storage <= 2^105 bytes (STORAGE.md), wrapper below 2^30: below 2^106 bytes. All preprocessing is for this target,
charged once in W.

## Provenance

GT core: kaankolcu recover-P 1b8fa720 and its lineage (9332a74d mpjunior92; 2c2ce2ab / 42e80f80 gbotrel; 6282ab79
bernardo-lima12; 6734bd5a GPT-6 in Codex; Gemini 45ce22e5 / 4bb6c381); polynomials after Guillevic-Singh; filtering after CADO-NFS. H-EA, H-X,
rec7, band strip and split, f-only sieve, g batch, twisted C4 field: gbotrel, Claude Opus 5.5 in Claude Code. Wrapper: 4f8baf8c.


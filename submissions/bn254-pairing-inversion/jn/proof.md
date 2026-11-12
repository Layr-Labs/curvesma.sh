# bn254-pairing-inversion: (6,2) STNFS on a twisted base field: expected work 2^90.086903716 units

Track bn254-pairing-inversion. Input: z = g^x in GT (g = e(P0, Q0), the profile's generator), x uniform in 0..r-1; output P in G1, Q in G2 with e(P, Q) = z. The score unit is C = 20,000 scalar AArch64 instruction equivalents (a borrowed normalisation unit:
every BN254 operation below is charged in instructions). No advice, smaller targets or amortization are used. Analytic, exploratory claim.

## Claim

W is all classical work of one attack on one full-production target. Under the stated hypotheses the attack succeeds with probability at least 0.39 over uniform challenges and fresh coins (under H5 as modeled an opened run succeeds with probability at least 0.997932; the run gate opens at most a share 0.390809 of runs), and

    E[W] <= 1 + (K/2^64) 2^91.442371104 < 2^90.086903716 units (C=20,000 instructions), K = 7209138687552583874.

Run gate: before any other work the attack reads one fresh 64-bit coin u and opens the run only when u < K = 7209138687552583874, the least integer with (K/2^64) x 0.997932 >= 0.39 (evidence.py checks in exact rationals that K - 1 fails). A closed run is a failure and costs 1 unit. In an opened run both individual logarithms run under one combined phase meter t = 2^112 units: when the logarithm phase reaches t the run aborts as failure. The work spent in that phase is min(D, t), with D the uncapped two-log work, so its expectation is at most E[D] <= 2^88.223193 units by the DESCENT.md recurrence (hypothesis H5), and the ledger charges that. Both the work and the success (section 7) of an opened run rest on H5 as modeled; no allowance factor is added or claimed. The other rows are expected upper bounds under the stated hypotheses and include failed trials and abort work.

The split phase alone has an interruptible budget of2^112 units; both individual logarithms share the 2^112-unit phase meter above. Before a finite split schedule would exceed its reserved charge, abort as failure. A split retry is one bounded curve schedule, never an inner unbounded loop. The 128-bit remaining-budget debit and branch require at most128 scalar instructions per invocation; every reserved schedule, trial-division block or ideal-lookup block costs at least128 instructions under the stated analytic kernel bounds. Doubling the entire split row therefore pays for the debit and its control, including the aborting check. Reserve and charge a whole finite schedule before invoking it, at the input-bitlength upper bounds already used by model16.split_cost. The field-map, collection, exact smoothness checks, linear algebra and output retain their explicit finite loops and existing cost models; their measured averages are not used as worst-case reservations.

Writing E_split for the doubled expected split row, Markov charges split abort at most E_split/2^112. A partial or unverified output is failure.

Conventions. Bounds round in the safe direction; evidence.py recomputes every row from the shipped statistics (stats16.b85, the constants of
model16.py, la_bn.py, descent16.py) and checks each declared number (declared.json) against its unrounded value; [n] cites its check lines.
What is measured at BN254 parameters and what is carried from kernels measured on another curve is listed in KERNELS.md.

## 1. Wrapper

1. Check z (canonical encoding, z != 0, z^r = 1). If z = 1 return (P0, infinity) (allowed exactly when z = 1).
2. Draw s uniform in 1..r-1 and t uniform in 0..r-1; g' = g^s, z' = z g^t (g = e(P0, Q0), the profile's generator). g' is uniform on GT minus
   {1} and z' uniform on GT, independent: the core's input law, whatever z is.
3. Map g' and z' into the core's tower (field-map.md).
4. Run the core on (g', z') under its meter and node counter (any failure is a failure): l = log_g'(z') mod r. Then x = s l - t mod r.
5. Q = [x]Q0 on the twist over F_p2 (left-to-right double-and-add, Jacobian coordinates, mixed additions, one inversion; ledger row "output": at most 2^26.05 instructions, ledger16.g2_output_bound with loose four-limb prices: F_p product <= 4,096, addition <= 256). Check e(P0, Q) = z (one pairing; a mismatch is a failure) and return (P0, Q). P0 is the public generator and Q = [x]Q0 lies in G2 because Q0 does.

GT has prime order r, so z = g^x for a unique x mod r and z' = g'^l with l = (x + t)/s. The field isomorphism preserves z' = g'^l; the final
check rejects any wrong output. One relation and linear-algebra phase serves both logarithms (section 6). Steps 1-3 and the pairing check of step 5 (validation, two GT exponentiations, one pairing, the field map of field-map.md, randomness) are charged in the wrapper row: at most 2^68 instructions; the scalar multiplication has its own row.

## 2. Tower, polynomials, automorphism ([0], [0a], [0b], [0c])

u = 4965661367192848881, p = 36u^4 + 36u^3 + 24u^2 + 6u + 1, r = 36u^4 + 36u^3 + 18u^2 + 6u + 1 ([0]: both equal the profile; ord_r(p) = 12).
K_h = Q[y]/h, h = y^6 - y^5 - 2y^4 + 5y^3 - 2y^2 - y + 1, irreducible mod p, disc h = 98000 = 2^4 5^3 7^2, O_h = Z[y] ([0c]: Dedekind at 2, 5, 7),
signature (2, 2), unit rank 3 (common12.UNITS). w = y^5 - 3y^3 + 2y^2 + 2y - 2 has norm 5. With F(z) = 36z^4 + 36z^3 + 24z^2 + 6z + 1 (p = F(u)):

    f = w^4 F(x^2/w) = 36x^8 + 36w x^6 + 24w^2 x^4 + 6w^3 x^2 + w^4,    g = x^2 - u w.

[0a]: (u w)^((p^6-1)/2) = -1 in F_p6 = F_p[y]/h, so F_p[y, x]/(h, g) = F_p12, and f(x^2 -> u w) = w^4 F(u) = w^4 p = 0 mod p: g divides f mod p.
Twisted symmetry ([0b]). sigma(y) = -y^5 + y^4 + 2y^3 - 5y^2 + 2y + 1 is an automorphism of K_h of order 2, c = y^4 - 2y^2 + 2y + 1 a unit with
c sigma(c) = -1 and sigma(w) = c^2 w; tau: (y, x) -> (sigma y, c x) maps f to c^8 f and g to c^2 g, tau^2 = (x -> -x). On F_p12 tau is
Frobenius^3 (sigma(y) = y^(p^3), (u w)^((p^3-1)/2) = c), so log tau(I) = chi log I with chi = p^3 mod r, chi^2 = -1. One relation stands for its
tau-orbit (4 pairs) and the factor base is quotiented by tau-orbits of 4 ideals.

A candidate is a - b x with (a, b) in O_h^2: N_f = |N(sum_k c_k w^k a^(8-2k) b^(2k))| (c = 36, 36, 24, 6, 1) and N_g = |N(a^2 - u w b^2)|, norms
from K_h to Q (common12.norms). f is not monic: its leading coefficient 36 is handled in the standard way (the ideals above 2 and 3 are in the
factor base, and the fixed ideal (1, x)^-1 of each relation cancels in the linear system).

Small primes. Murphy alpha over the odd primes ([14], alpha12.py): at 3, 5, 7 (they divide 36, N(w) or disc h) the mean valuations are measured
on 37,979 exact primitive pool pairs and taken at mean - 2 SE; elsewhere Kummer-Dedekind root counts (a simple root over a prime ideal of norm
q contributes q/(q^2 - 1); g has a double root over the primes dividing u: 1/(q + 1)), checked against exact valuations of 3,000 primitive
pairs ([14v]). The partial sums converge slowly (f / g: p < 2,000: (-1.250, 2.633); p < 20,000: (-0.622, 2.589); p < 200,000: (-0.473, 2.718); p < 1,000,000: (-0.458, 2.709)); the charged values (-0.4, 2.76) are the sums over p < 10^6
plus 0.05. The prime 2 is priced exactly per cell: 2 = P^3 with N(P) = 4, v_P(36) = 6 and w is a unit at P, so on primitive pairs N_f carries
2^12 exactly when b is in P and a is not, and 2^0 otherwise; N_g carries 2^2 or 2^0. All 37,979 primitive pool rows agree ([14v]: (0,0) 22,658, (0,2) 7,515, (12,0) 7,806).

## 3. Relation collection

Relations: primitive (a, b) with N_f and N_g both B-smooth, B = 2^58. Factor base: the prime ideals of relative degree 1 and norm <= B on both
sides, U' = 2^51.7080 tau-orbits.

Region. Exact rational forms Q_a, Q_b (stats16 config; place-weighted forms (1/6) sum_e om_e |x_e|^2 with om_a = |w_e|^-1/2 2^0.8,
om_b = |w_e|^1/2 2^-0.8, rounded to 10 digits: the rounded forms are the definition), A^2 = 300000. a (power coordinates): a_5 >= 1,
Q_a(a) <= theta A^2, theta = 2362/8000. b: c_5 >= 1 and Q_b(b) <= (17/20) A^2 - Q_a(a), where c = SD b are the sieve coordinates (common12.SD /
SB, unimodular; the row direction c_0 is the shortest vector 1 of Q_b, the line index c_5 its shortest dual vector). The canonical a is the
(Q_a, lex)-least element with a_5 > 0 of its <+-1, units, sigma> class (exact enumeration, ideal_canonical.py); with c_5 >= 1 this is one point
per <+-1, units, tau> orbit (tau^2(a, b) = (a, -b)); non-primitive pairs are excluded by an exact HNF test. Only canonical a are sieved.

Pools ([2], [2a]). Uniform integer pairs on the region with theta_pool = 3/8 (common12.draw_a, draw_b; 100,000 pairs per seed, every draw kept):
S (seeds 101-104) sizes theta and gives the work, J (seeds 111, 112, 113, 114, 121, 122) is judged. A draw counts when a is canonical; then norms, 2-parts and
primitivity are exact. stats16 holds per-seed counts and hashes, the sums, 8 anchor rows and replay lists.
Counts ([3]). Points per volume of the two half-ellipsoids, by the sampler's acceptance experiment (uniform integer points of the frame box,
exact hit test; countmc.py): a: 1e+09 trials, [0.996343, 0.999151]; b: 4e+09 trials, cumulative on a grid of radii^2
(lam - t/4000) A^2 covering every slice of the region (count nondecreasing, volume ~ radius^6 between grid points): [0.995486, 1.000698]. Each
point carries 6 SE (union bound over 2368 one-sided events < 10^-5). Factors 0.9905 (yield, slice cells from below) and
1.002 (work, population; a's end raised to 1 so that it also bounds a b-slice alone). Fallback: lattice_count_bounds.py certifies
(exact rationals) (0.9211, 1.0797), which contains both ([3]).
Sizing ([4s]). theta is the smallest t/8000 with margin_S >= 0.06 + 2 sqrt(2) SE_S over 1.095 U' (rule and count factor 0.9905 fixed before
S was read): 2362/8000 (margin 0.160070 >= 0.159833; one step below: 0.159384 < 0.159879).

Sieve (f ideals only; KERNELS.md). Per slice (one canonical a): line sieve of the lines c_5 = k for p in [2^14, 2^40); record walk for
p in [2^40, 2^57). Primes of [2^57, 2^58] are not sieved: they stay in the 75-bit slack and are found by the survivor tests below. A canonical a's run is the maximal integer t-interval containing 0 with a + t canonical and inside the limit theta A^2
(runs16.py; e_0 = 1 is the shortest vector of Q_a). On a run a prime ideal (p, s, rho)'s hits in all its slices form one coset of a 7-dim
lattice (t = rho b(s) - a(s) mod p), walked once per run and ideal pair (p, s, +-rho) (f is even in x) with one setup: (a, b) is in (p, s, -rho) iff (a, -b) is in (p, s, rho); each walk covers the full symmetric b-ellipsoid, hits of both ideals are charged (factor 2 in l) and an empty coset means no hit in the full ellipsoid. f has on average one ideal per prime: 0.5 pairs per prime are set up, as for empty cosets. Run rule: a run of
16 to 1024 a-values all of which are in the canonical list uses the run walk in each quarter-octave band of p in which the model
density l = log2(2 x run cells x lower/upper / p) is at least 2.5: the rec7 kernel where l >= 4.99 and p < 2^56 (its measured range), dadon's walk7 kernel (same coset, same record stream; KERNELS.md) where 2.5 <= l < 4.99 and on [2^56, 2^57); other bands and other runs use the slice walk (and empty cosets below
its table). 2.38% of the cells lie in runs shorter than 16; the longest sampled run has 247 a-values ([5r]). Runs of 2 to 15 a-values also use walk7 (l >= -0.5) at their own prices (REC7S); single a-values stay on the slice walk.
Contamination ([5r]). In 2,500 a-uniform runs (seeds 5001-5005, the first 100 canonical a of each inside the limit, and seeds 9101-9108, the first 250 of each) the exact
predicate found no non-canonical a-value among 201,376; 20 runs end next to a canonical a (runs16.run is a sufficient test: that a starts its own run, as in the work rows): contaminated runs hold under
1 - 0.05^(1/2500) = 0.00120 <= RUN_ALLOW = 0.0012 of the canonical a (95%); every cell pays RUN_ALLOW x the largest slice-walk
cost of a slice (model16.contam_sup).

Thresholds, two levels (model16.two_level, embedding_bounds.py). A cell passes when its score reaches T = lower - 16 x 75, lower a
certified lower bound of 16 log2 N_f (Taylor shift per place in b^2, outward rounding) for exactly the rows and cells it is used on: every
cell whose unsieved mass is at most 75 bits is kept whatever the bound's width, a wide bound only sends more cells to the exact norms (charged
at 75 + 32 bits). Level 1: the rows of an aligned 8 x 8 x 8 box share the hull of the hulls of its 2 x 2 x 4 sub-boxes; in an aligned 256-cell
block each 64-cell node meeting the box is evaluated on it; certified (<= 32 bits wide) it is a leaf, else halved down to 16 cells. Level 2: a
failed 16-cell node runs the per-sub-box bisection down to single cells and rows. With no level 2 in the block and a total width <= 32 bits the
block is one leaf. Evaluations and scan leaves are counted exactly per work row ([2a]); no pool cell was left uncertified.

Survivors. Each survivor gets its exact norms (C_NORMS = 27,000, measured kernel) and is rechecked: it goes on only if 16 log2 N_f minus
its score is at most 16 x 75 (within the rounding of the at most (n_f - v_2)/14 scored ideals, 1/16 bit each). Then: (i) g batch: N_g is
B-smooth iff (P_g mod N_g)^(2^10) = 0 mod N_g, P_g the product of all primes <= B, by the exact scaled remainder tree below in waves of 2^48
norms, with 8-limb leaves since N_g < 2^501.9 < 2^512 on the region (analytic: norm_sup multiplies over the six places the Cauchy-Schwarz bounds of |a_e^2 + u w_e b_e^2| from Q_a, Q_b; an N_g >= 2^512 would use a 9-limb leaf, none exists) (C_GBATCH = 2^23.07, product price mul_F); (ii) only if g passes: trial division of N_f by the 82,025 primes below 2^20 (C_TDIV = 2^23.45), a strip
with the primes of (2^20, 2^40] and one with (2^40, 2^58] (same tree; the stripped part is gcd(z, (P mod z)^(2^5) mod z); C_BATCH = 2^22.7,
C_BATCH2 = 2^24.26); the pair is a relation iff the cofactor is 1. g passes with probability at most q_g = min(1, 100 rho(u_g)) per pool
row (H1; 100 is a margin). A kept relation is then factored (section 8, split row).
Lemma (exact SRT). Product tree Z_v (L levels); node v holds Y_v of b_v = 64 (n_v + G) bits (n_v limbs of Z_v, G = L + 2). Root: Newton
quotient q of P 2^b by Z, corrected by <= 2 add/subtracts, Y_root = q mod 2^b. Child c, sibling s: Y_c = floor((Y_v Z_s mod 2^(b_v))/2^(b_v - b_c)).
Leaf: R = ceil(z Y 2^-b) (0 if R = z) = P mod z. Proof: eta_v = (frac(P/Z_v) - Y_v 2^(-b_v) mod 1) Z_v gains < 2^(-64G) per level, so
eta < (L + 1) 2^(-64G) < 1. Cost (H-FFT): products at F(N) = 8.5 N log2(N)^2 + 1200 instructions for N limbs (+ 504 N above 2^26 limbs); in the g batch tree mul_F(N) = min(F(N), F_ntt(N)) with our ntt3 kernel's exact count (KERNELS.md);
model16.gbatch_cost; [5t] replays the tree on planted leaves and checks the measured FFT points against F.

Laws (threshold_model.py; H1'). Euler-product models without the prime 2, read at n = n_f - v_2 with slack 75 - v_2: primes below 2^14 are
unsieved with all their powers, primes of [2^14, 2^57) sieved once, primes of [2^57, 2^58] unsieved. Yield: r_f = P(unsieved mass <= slack - rounding | B-smooth, sieved part
squarefree) (r_f = 0.7895 at n = 417 bits, slack 75). Work (H-EA): the same model over all integers bounds the survivors per cell (model16.row_cost: at slack + 32 for
the norms, at slack + 1/4 after the recheck).

Yield ([4], pool J). rels = (va/n) vb x 0.9905 x sum over J's canonical primitive region rows of (lam - qa)^3 rho(u_f) rho(u_g) r_f, with
u_s = (n_s - v_2 + 1 + kappa + alpha_s)/58: the odd part of a norm whose 2-part is known is priced as a random integer one bit larger plus
kappa = 0.07 >= kappa2(u) (sum_k 2^-(k+1) rho(u - (k - 1 - kappa)/58) = rho(u); [14v]); Dickman rho; cells with qa < 0.02 are sieved and
charged, their relations discarded. Estimate 2^52.02659 (SE 0.02466 bits) against 1.095 U' = 2^51.83890: margin 0.1876 >= max(0.06, 2 SE).
The count lowered by two SE, 2^51.9764, is the descent's input (>= 2^51.97); relation-weighted mean norms 416.08 / 470.79 bits.

Cost ([5]). Mean 43.80 instructions per cell over S's 14,029 work rows (record walk 21.56, line 7.61, survivors 6.06, norm bounds 2.35, scan 2.15, apply 2.04, routing 1.09, layout 0.59, setup 0.26, empty cosets 0.09)
over 2^99.782 cells: collection <= 2^90.983040 units (estimate 2^90.947416 times the allowance 1.025 >= 1 + 2 SE <= 1.02483; total SE at most
1.242%, [4a]; work2 is an upper bound).

## 4. Filtering and linear algebra

Filtering (la_bn.py, filtsim.py; [6], [6r]). Purge (to excess 0.0001 U'), clique removal and Markowitz merge are simulated on this
instance's occurrence profile (B = 2^58, mean relation norms 417 / 471 bits, above J's relation-weighted means; excess
1.095) as the tau-quotient system: tau^k(I) enters I's orbit column with coefficient chi^k (entries a + b chi, merged |a|, |b| <= 255).
x = N(400)/U' per seed (seeds fixed before the first run, none dropped): 2^17 mean 0.01845 (2 seeds), 2^18 0.01624 +- 0.00022 (8), 2^19
0.01518 +- 0.00013 (8), 2^21 0.01373 +- 0.00007 (8), 2^22 0.01304 +- 0.00007 (8); 2^21 and 2^22 by filt29, value-identical to filtsim on
all 18 runs of 2^17..2^19 (FILTER20.md); the mean falls with scale. H3: production is not above the two largest scales: charged NU =
0.01386 >= 0.013858, the larger mean + 2 SE of 2^21 and 2^22 (an extrapolation of 30 bits). N <= NU U' = 2^45.5350 < 2^64/320: one owner. The attack aborts if
N > NU U' (never under H3; P_CAP = 1e-9 charged). Merge to mean weight 400, charged 401.

The square system A = [M | S] has N + 7 rows, N orbit columns and 7 dense Schirokauer columns: the chi-eigenpart of the unit groups,
(rk K_f - rk K_f^+ + rk K_g - rk K_h)/2 = (23 - 11 + 5 - 3)/2 = 7 from the signatures K_f (0, 24), K_f^+ (0, 12), K_g (0, 6), K_h (2, 2)
([6b]; tau^2 = -1 on the sigma-odd part, so its chi and chi^-1 eigenspaces have equal dimension). r is unramified in K_f and K_g ([6b]:
the absolute polynomials of degrees 48 and 12 are squarefree mod r). Block Wiedemann (m, n) = (1024, 40) with Y=[S,BZ] and a sparse random X whose seven pad coordinates occupy distinct blocks, specified in RANK-BRIDGE.md; (m, n) minimises linear algebra + generator over
m = 2^5..2^13, n = 10..200 (la_bn.bw_grid). Per row and vector: 401 x 8.42 + 176 (row tail, one owner) + 200 (chi product) = 3752.42
instructions (STORAGE.md); k5 = 8.42 >= 8.419, the measured mean + 2 SE of the entry price on the merged BN254 system ([6]). Linear algebra
<= 2^88.807039 units; generator (half-gcd matrix Berlekamp-Massey) <= 2^81.729.

Schirokauer maps. R_f = (Z/r^2)[y, x]/(h, f/36) has rank 48, R_g rank 12; eps = lcm(r^d_i - 1) < r^48 (12,192 bits; g: 3,048). The 7
columns are fixed F_r-linear forms of (phi^eps - 1)/r mod r, phi = a - b x, by left-to-right powering on eight-limb coefficients (r^2 <
2^508); per bit < 2^21.08 instructions on f, < 2^17.75 on g (schoolbook counts; the reductions by h are allowances, la_bn.schirokauer_cost):
< 2^34.69 per relation.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent, so the finite relation set and filtered matrix are fixed. H-R explicitly assumes rank M=N and rank[M|S]=N+6 for this actual finite structured set, with a one-dimensional virtual-log kernel whose map part is nonzero. H-X explicitly supplies finite right-controllability and left-observation minors. Neither follows from smooth relation counts, class-number coprimality or map injectivity. RANK-BRIDGE.md supplies the full algorithm, polynomial-degree derivation and primary-framework support. Conditional on these premises, independently sampled projection/random-source entries give failure at most4N/2^128. No probability that H-R/H-X hold is invented; they are disclosed structural applicability heuristics under exploratory review. Verification detects errors but does not establish these premises. Purged logs are reconstructed and unresolved low-occurrence ideals are handled in descent.

## 6. Individual logarithms

DESCENT.md and finalmodel4.py specify the fixed62-type production plan, original amortized node sieves and independent single-cell boot. The DESCENT.md recurrence gives expected work at most 2^88.223193 units for both individual logarithms; the ledger charges this figure (H5); the meter 2^112 can only truncate the work. Each boot trial uses a fresh uniform exponent, one certified boot lattice and one sampled cell, charging exponentiation, lattice construction, complete finite cofactor schedules and all children. It does not use BALL_PAIR or a whole-ball boot waiting-time assertion. DESCENT.md gate G is an exact lossless test (node mu, phi, psi unchanged). Gate A is a boot early abort: cells failing it are discarded and p_trial counts only qualifying cells that pass (model x 0.90). Gate E charges the boot's children at their modeled expected sizes among cells passing gate A (x 1.25) instead of the per-window maxima. Node rows keep the per-window maxima.

At nodes the count premise is the conditional PGF envelope with alpha=1/2 and the stated shape. THRESHOLD-REPAIR.md combines count tails with the separate sequential availability hazard by P(fail)<=P(M<k)+psi^k; it does not use the rejected PGF(psi) substitution. It bounds the probability that no candidate relation has all usable descendants; it is a probabilistic heuristic rather than an exact count distribution. DISC and its pooled per-type confidence argument are removed. Every candidate survivor is charged; there is no negative-binomial lazy cofactor discount. Node-count and repeated-ideal aborts are explicit, and their probabilities remain part of the success calculation. There is no fixed auxiliary family or universal finite-family positivity premise.

Finite ECM schedules, the stored-prefix extraction and H-SEP are specified in DESCENT.md; stage2-summary.json reports real stage2/composite-modulus experiments at 36/40-bit factors and up to 1746-bit moduli, whose extrapolation to production factor sizes is H-SEP.

## 7. Success

For an opened run, ledger16.totals subtracts from 1: main-phase matrix failure (4N/2^128), numerical allowances, recursive descent and boot failure, node-counter and repeated-ideal aborts, split abort E_split/2^112, and the log-phase cap. For the cap: the two logarithms are computed for the wrapper's independent uniform targets g' and z' with separate fresh coins, and the second uses nothing computed by the first, so their uncapped works D1, D2 are independent; the DESCENT.md recurrence is per logarithm (2^88.223193 is twice the one-log bound), so E[D_i] <= 2^87.223193. The meter is reached only if D1 + D2 >= t = 2^112. By the lemma below with a = 2^(87.223193-112) <= 0.000000035: P(cap) <= max(2a - a^2, a/(1 - a)) <= 0.000000070. All modeled losses other than the cap total at most 0.0020677. So the model gives an opened run success at least 0.997932.

The gate uses this modeled figure: the gate coin is independent of the run, so the attack succeeds with probability at least (K/2^64) x 0.997932 >= 0.39. This figure and the charged work both rest on H5 as modeled. A closed run and every abort is failure, with its work charged.

Lemma (two-variable case of Samuels' inequality, 1966; the form and proof follow kaankolcu's accepted GT package 8c816a9a). X, Y >= 0 independent, E X <= a t, E Y <= a t, a <= 1/2. Then P(X + Y >= t) <= M(a) = max(2a - a^2, a/(1 - a)).
Proof. Scale to t = 1 and replace X, Y by min(X, 1), min(Y, 1): the event and the independence are unchanged and the means do not increase. G(x) = P(Y >= 1 - x) is nondecreasing on [0, 1], so its least concave majorant G^ is nondecreasing and P(X + Y >= 1) = E G(X) <= E G^(X) <= G^(E X) <= G^(a) (Jensen). G^(a) is the supremum over 0 <= x1 <= a <= x2 <= 1 of (1 - l) b1 + l b2 with d = x2 - x1, l = (a - x1)/d, b1 = G(x1), b2 = G(x2) (if d = 0, G(a) <= a/(1 - a) by Markov). Here 0 <= b1 <= b2 <= 1, and Y >= (1 - x1) 1[Y >= 1 - x1] + (1 - x2) 1[1 - x2 <= Y < 1 - x1] gives, taking expectations, d b1 + (1 - x2) b2 <= a. The linear objective is largest at a vertex of this polygon. (i) b1 = 0: the value is l b2 <= (a/x2) min(1, a/(1 - x2)) (l decreases in x1); for x2 <= 1 - a this is a^2/(x2(1 - x2)) <= a/(1 - a) since x2(1 - x2) >= a(1 - a) on [a, 1 - a]; for x2 >= 1 - a it is <= a/x2 <= a/(1 - a). (ii) b1 = b2 = a/(1 - x1) <= a/(1 - a). (iii) b2 = 1, b1 = (a - 1 + x2)/d, which needs x2 >= 1 - a: with c = x1 and v = 1/d the value is 1 - (1 - a - c) v + (a - c)(1 - a - c) v^2, convex in v (c <= a <= 1/2), so it is largest at an end of x2 in [1 - a, 1]. At x2 = 1 - a it is (a - c)/(1 - a - c) <= a/(1 - a); at x2 = 1 it is psi(c) = (2a - a^2 - (1 + a) c + c^2)/(1 - c)^2, and psi'(c) has the sign of (1 - a)(c - (1 - 2a)), so on [0, a] psi <= max(psi(0), psi(a)) = max(2a - a^2, a/(1 - a)). Every vertex value is at most M(a). QED

## 8. Ledger

| component | log2 expected units (C=20,000) |
|---|---:|
| run-walk bases | <= 70.645293 |
| canonical a list | <= 76.830772 |
| certificates | <= 50.000001 |
| relation collection | <= 90.983040 |
| factor base and ideals | <= 62.110830 |
| matrix generator | <= 81.728608 |
| product tree | <= 58.017996 |
| linear algebra | <= 88.807039 |
| filtering/back-substitution | <= 68.973013 |
| Schirokauer maps | <= 73.256328 |
| kept-relation factorization (including doubled control budget) | <= 83.114298 |
| two individual logarithms (expected work under H5; meter 2^112) | <= 88.223193 |
| all rows of an opened run together | <= 91.442371104 |
| after the run gate, plus 1 unit | 90.086903716 |
| verification | <= 11.000001 |
| wrapper | <= 53.712288 |

Only the split table row includes an extra factor two for debit/control overhead. The following underlying split kernel price precedes that extra factor; other rows retain their original allowances. Multiplying units by20,000 gives instructions. Canonical a list: 2^42.8 instructions per integer a with a_5 >= 1, Q_a(a) <= theta A^2
(2^48.32 of them; the exact test of ideal_canonical.py and the sort of the list). Factor base: 2^23.4 per prime (roots of h and f, basis
values, ideals), pi(B) < 2^53.00. Filtering streams and back-substitution: 2^24.4 per stored entry, 2^6 entries per relation, <= 2^52.87
relations. Run-walk bases: 4 reduced 7-dim bases per ideal of [2^40, 2^58) at <= 5,000,000 instructions and 5 for walk7 at <= 1.0e9 (KERNELS.md).
Relation split: a kept relation's two norms are trial-divided by every prime of (2^20, 2^40] (< 2^35.53 primes at
156 instructions with a 2^(64 i) mod p table), each hit divided out to full power, the primes of (2^40, B] found by ECM (ecm16.py
schedule, 99.9% per prime, repeated until found or the split-only meter would be exceeded: the uncapped expectation is at most 1.002 schedules), Miller-Rabin and ideal lookups: 2^43.55 instructions per
relation in expectation. This is not a finite retry bound: the split-only meter interrupts every retry loop, and its abort loss is charged by Markov. P_g product tree: 2^72.31 instructions (H-FFT). The expected rows and the split-only and log-phase cap losses are recomputed by ledger16.totals; the non-log rows are not asserted to bound every run.

## 9. Heuristics

H1 smoothness (Dickman rho with Murphy alpha over the odd primes, exact 2-parts and kappa; exact sampled norms). H1' retention (Euler-product
model of the unsieved part). H-EA (work only): the same model prices the survivors (also those of descent gates G and A). H-FFT (work only): F(N) bounds the multiplication work
past 2^26 limbs; for the g batch tree min(F, F_ntt): exact count, the hypothesis is only its distributed-execution allowance. H2 region (uniform pools, exact canonical predicate, measured counts with the certified bracket as fallback; sized on S,
judged on J; work with the 1.025 allowance). H-R finite relation spanning and finite H-X applicability are specified separately in RANK-BRIDGE.md; they are not implied by the following H3 assumptions. H3 matrix (filtering ratio at production at most the larger mean + 2 SE of the two largest simulated scales; Schirokauer maps injective on units modulo r-th powers; class numbers of K_f, K_g prime to r (assumed, not checked here)). H4 kernel prices (KERNELS.md: the run walk, the norm kernel and the
filtering ratio are measured at BN254 parameters; the other sieve and linear-algebra kernels are the instruction sequences measured for the
BLS12-381 instance of the same sieve, used inside their measured density ranges; the walk7 run-walk kernel is measured at BN254 on the cells it serves; the norm-bound kernel carries a 25% allowance). H-X (block
Wiedemann, section 5). H5 (expected work and success of the logarithm phase) conditional node PGF envelope with alpha=1/2 and boot single-trial smoothness, including gate A's pass probability (window composition restricted to <= tau bits above 2^x, x 0.90) and gate E's expected child counts (x 1.25) (DESCENT.md; its conditional next-relation premise is measured in ADAPTIVE.md's conditional table: after an earlier failed relation the observed/envelope ratio is 0.964, 95% 0.839-1.095, and conditional child work 0.950 (95% up to 0.983) of the envelope). H-SEP is the conditional ECM extraction-vs-order-hit lower bound and its stated extrapolation. Neither is claimed to follow from exact norm measurements alone.

## 10. Memory

Peak full-plan storage is below2^100 bytes (DESCENT.md; largest sieve 2^85.417380). The main phase stays below it: a run's sum arrays (<= 1024 slices of <= 2^56.54 cells, 4 bytes per position, line-box overhead allowance x 32: < 2^73.6 bytes), the matrix (2^45.81 rows x 401 entries x 8 bytes < 2^57), the generator's sequence (< 2^66), stored relations (< 2^58). Wrapper below 2^30 bytes. All preprocessing is for this target, charged once in W.

## Provenance

Derivative of gbotrel's package 112f52f4-2afe-47fd-9862-fa16451f70b0 (90.270795659) with two changes in the collection row: the g batch tree priced with our ntt3 kernel (C_GBATCH 2^23.07) and runs of 2 to 15 a-values on walk7 (KERNELS.md): collection 2^90.983040 instead of 2^91.230251. walk7 and ntt3 are ours.

This is gbotrel's BN254 package with the factor-base bound 2^58 (declared dependency) with the logarithm-phase accounting of kaankolcu's accepted package of this line (declared dependency; meter as in dadon's accepted package): the phase is charged at its modeled expected work under the package's own hypothesis H5, meter 2^112. Descent-side change of this version, on top of our stnfs-bn10rw2 (declared dependency; stnfs-bn10ez plus the collection-side changes below): gates G, A and E in the descent evaluator (model4core; DESCENT.md; part=eatoy) and the plan and boot re-optimised under them (seeds: our plan and dadon's accepted 0d4159ea; phi <= 0.04, psi < 0.8, radii <= 192, sieve radii <= 160). Gate A's pass probability and gate E's child counts are named in H5; part=plan prints the rows and the boot.

Run-walk pricing of dadon's accepted package 258122fe (declared dependency): his walk7 kernel prices the run-rule cells that rec7 does not cover (2.5 <= l < 4.99, and [2^56, 2^57)), from his rec7own.txt campaign; the code here is ours and [10w] recomputes the prices.

Collection-side changes carried from stnfs-bn10rw2: all 14,029 region rows of pool S priced (3,514 before; 4 regenerated anchor rows recomputed by [2c]); run walk at its BN254 measurement below 2^56 (KERNELS.md; method of kaankolcu's accepted f7d11cdc, declared dependency); 8-limb g-batch leaves; setups per +-rho pair of f; RUN_ALLOW 0.0012 from 2,500 contamination runs.

The sieve, kernels, filtering and descent machinery are those of the (6,2) twisted-field STNFS packages for BLS12-381 pairing inversion by
gbotrel (with Claude in Claude Code) and their lineage on that benchmark (kaankolcu: run gate, expected total work, two-level
thresholds, norm-kernel design; mpjunior92: logarithms at their expectation; polynomials after Guillevic-Singh and Barbulescu-Duquesne;
filtering after CADO-NFS). New here: the BN254 field and pair, its pools, counts, alpha, run-walk campaign, norm kernel, filtering campaign and descent plan.

ECM-cell qualification, the exact residual test (gate G), perfect powers and both final ECM gcds: DESCENT.md.

RANK-BRIDGE.md states H-R and the finite H-X bridge. ADAPTIVE.md supplies the deterministic-order conditional campaign for H5 (uncensored recursive descents, 627 attempts after an earlier failed relation), within its stated toy scope.

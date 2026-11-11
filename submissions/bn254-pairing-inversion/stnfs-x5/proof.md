# bn254-pairing-inversion: (6,2) STNFS on a twisted base field: expected work 2^96.845632592 units

Track bn254-pairing-inversion. Input: z = g^x in GT (g = e(P0, Q0), the profile's generator), x uniform in 0..r-1; output P in G1, Q in G2 with e(P, Q) = z. The score unit is C = 20,000 scalar AArch64 instruction equivalents (a borrowed normalisation unit:
every BN254 operation below is charged in instructions). No advice, smaller targets or amortization are used. Analytic, exploratory claim.

## Claim

W is all classical work of one attack on one full-production target. Draw an independent uniform 64-bit word v and fail unless v < K = 7510732454374047492 (run gate; a closed gate costs under one unit). Premise H5x5: the true expected descent work is at most 5 times the modeled 2^95.805191 and every other failure term is at most 10 times its modeled value; so E[D] <= 2^98.127120 units for the uncapped two-log work D. Under the stated hypotheses, with H5 used only through H5x5, the attack succeeds with probability at least 0.39 over uniform challenges and fresh coins, and

    E[W] <= 1 + (K/2^64) 2^98.141973276 < 2^96.845632592 units (C=20,000 instructions).

Both individual logarithms run under one combined hard phase meter t = 2^106 units: reaching t aborts the run as failure, so that phase costs at most t units in every run. This every-run cap is distinct from the score, an expectation: the ledger charges the phase at E[min(D, t)] <= E[D] <= 2^98.127120. Work and success both rest on H5x5, the recurrence value with a factor 5 of margin; the factor is the allowance for the sequential-hazard and conditional child-work envelopes behind the recurrence's expected work (DESCENT.md), whose half-scale count envelope, full survivor charging and abort losses are unchanged. Multiplying the modeled cap loss 0.000854 by 5 and the remaining losses 0.003787 by 10 leaves opened-run success at least 1 - 5(0.000854) - 10(0.003787) = 0.957860 (section 7). The gate is independent of the challenge and of every later coin; K is the least integer with (K/2^64) x 0.957860 >= 0.39. The other rows are expected upper bounds under the stated hypotheses and include failed trials and abort work.

The split phase alone has an interruptible budget of2^112 units; both individual logarithms share the phase meter t above. Before a finite split schedule would exceed its reserved charge, abort as failure. A split retry is one bounded curve schedule, never an inner unbounded loop. The 128-bit remaining-budget debit and branch require at most128 scalar instructions per invocation; every reserved schedule, trial-division block or ideal-lookup block costs at least128 instructions under the stated analytic kernel bounds. Doubling the entire split row therefore pays for the debit and its control, including the aborting check. Reserve and charge a whole finite schedule before invoking it, at the input-bitlength upper bounds already used by model16.split_cost. The field-map, collection, exact smoothness checks, linear algebra and output retain their explicit finite loops and existing cost models; their measured averages are not used as worst-case reservations.

Writing E_split for the doubled expected split row, Markov charges split abort at most E_split/2^112. A partial or unverified output is failure.

Conventions. Bounds round in the safe direction; evidence.py recomputes every row from the shipped statistics (stats16.b85, the constants of
model16.py, la_bn.py, descent16.py) and checks each declared number (declared.json) against its unrounded value; [n] cites its check lines.
What is measured at BN254 parameters and what is carried from kernels measured on another curve is listed in KERNELS.md.

## 1. Wrapper

0. Run gate (Claim): draw v; if v >= K, fail.
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

Relations: primitive (a, b) with N_f and N_g both B-smooth, B = 2^56. Factor base: the prime ideals of relative degree 1 and norm <= B on both
sides, U' = 2^49.7600 tau-orbits.

Region. Exact rational forms Q_a, Q_b (stats16 config; place-weighted forms (1/6) sum_e om_e |x_e|^2 with om_a = |w_e|^-1/2 2^0.8,
om_b = |w_e|^1/2 2^-0.8, rounded to 10 digits: the rounded forms are the definition), A^2 = 300000. a (power coordinates): a_5 >= 1,
Q_a(a) <= theta A^2, theta = 2686/8000. b: c_5 >= 1 and Q_b(b) <= (17/20) A^2 - Q_a(a), where c = SD b are the sieve coordinates (common12.SD /
SB, unimodular; the row direction c_0 is the shortest vector 1 of Q_b, the line index c_5 its shortest dual vector). The canonical a is the
(Q_a, lex)-least element with a_5 > 0 of its <+-1, units, sigma> class (exact enumeration, ideal_canonical.py); with c_5 >= 1 this is one point
per <+-1, units, tau> orbit (tau^2(a, b) = (a, -b)); non-primitive pairs are excluded by an exact HNF test. Only canonical a are sieved.

Pools ([2], [2a]). Uniform integer pairs on the region with theta_pool = 3/8 (common12.draw_a, draw_b; 100,000 pairs per seed, every draw kept):
S (seeds 101-104) sizes theta and gives the work, J (seeds 111, 112, 113, 114, 121, 122) is judged. A draw counts when a is canonical; then norms, 2-parts and
primitivity are exact. stats16 holds per-seed counts and hashes, the sums, 8 anchor rows and replay lists.
Counts ([3]). Points per volume of the two half-ellipsoids, by the sampler's acceptance experiment (uniform integer points of the frame box,
exact hit test; countmc.py): a: 1e+09 trials, [0.996343, 0.999151]; b: 4e+09 trials, cumulative on a grid of radii^2
(lam - t/4000) A^2 covering every slice of the region (count nondecreasing, volume ~ radius^6 between grid points): [0.995206, 1.000936]. Each
point carries 6 SE (union bound over 2692 one-sided events < 10^-5). Factors 0.9905 (yield, slice cells from below) and
1.002 (work, population; a's end raised to 1 so that it also bounds a b-slice alone). Fallback: lattice_count_bounds.py certifies
(exact rationals) (0.9211, 1.0797), which contains both ([3]).
Sizing ([4s]). theta is the smallest t/8000 with margin_S >= 0.06 + 2 sqrt(2) SE_S over 1.095 U' (rule and count factor 0.9905 fixed before
S was read): 2686/8000 (margin 0.150055 >= 0.150054; one step below: 0.149064 < 0.150114).

Sieve (f ideals only; KERNELS.md). Per slice (one canonical a): line sieve of the lines c_5 = k for p in [2^14, 2^40); record walk for
p in [2^40, 2^56). A canonical a's run is the maximal integer t-interval containing 0 with a + t canonical and inside the limit theta A^2
(runs16.py; e_0 = 1 is the shortest vector of Q_a). On a run a prime ideal (p, s, rho)'s hits in all its slices form one coset of a 7-dim
lattice (t = rho b(s) - a(s) mod p), walked once per run and ideal pair (p, s, +-rho) (f is even in x) with one setup. Run rule: a run of
16 to 1024 a-values all of which are in the canonical list uses the run walk in each quarter-octave band of p in which the model
density l = log2(2 x run cells x lower/upper / p) is at least 4.99; other bands and other runs use the slice walk (and empty cosets below
its table). 1.91% of the cells lie in runs shorter than 16; the longest sampled run has 254 a-values ([5r]).
Contamination ([5r]). In 500 a-uniform runs (seeds 5001-5005, the first 100 canonical a of each inside the limit) the exact
predicate found no non-canonical a-value among 43,225 and both neighbours of every run non-canonical: contaminated runs hold under
1 - 0.05^(1/500) = 0.00597 <= RUN_ALLOW = 0.006 of the canonical a (95%); every cell pays RUN_ALLOW x the largest slice-walk
cost of a slice (model16.contam_sup).

Thresholds, two levels (model16.two_level, embedding_bounds.py). A cell passes when its score reaches T = lower - 16 x 60, lower a
certified lower bound of 16 log2 N_f (Taylor shift per place in b^2, outward rounding) for exactly the rows and cells it is used on: every
cell whose unsieved mass is at most 60 bits is kept whatever the bound's width, a wide bound only sends more cells to the exact norms (charged
at 60 + 32 bits). Level 1: the rows of an aligned 8 x 8 x 8 box share the hull of the hulls of its 2 x 2 x 4 sub-boxes; in an aligned 256-cell
block each 64-cell node meeting the box is evaluated on it; certified (<= 32 bits wide) it is a leaf, else halved down to 16 cells. Level 2: a
failed 16-cell node runs the per-sub-box bisection down to single cells and rows. With no level 2 in the block and a total width <= 32 bits the
block is one leaf. Evaluations and scan leaves are counted exactly per work row ([2a]); no pool cell was left uncertified.

Survivors. Each survivor gets its exact norms (C_NORMS = 27,000, measured kernel) and is rechecked: it goes on only if 16 log2 N_f minus
its score is at most 16 x 60 (within the rounding of the at most (n_f - v_2)/14 scored ideals, 1/16 bit each). Then: (i) g batch: N_g is
B-smooth iff (P_g mod N_g)^(2^10) = 0 mod N_g, P_g the product of all primes <= B, by the exact scaled remainder tree below in waves of 2^48
norms (C_GBATCH = 2^24.07); (ii) only if g passes: trial division of N_f by the 82,025 primes below 2^20 (C_TDIV = 2^23.45), a strip
with the primes of (2^20, 2^40] and one with (2^40, 2^56] (same tree; the stripped part is gcd(z, (P mod z)^(2^5) mod z); C_BATCH = 2^22.7,
C_BATCH2 = 2^24.11); the pair is a relation iff the cofactor is 1. g passes with probability at most q_g = min(1, 100 rho(u_g)) per pool
row (H1; 100 is a margin). A kept relation is then factored (section 8, split row).
Lemma (exact SRT). Product tree Z_v (L levels); node v holds Y_v of b_v = 64 (n_v + G) bits (n_v limbs of Z_v, G = L + 2). Root: Newton
quotient q of P 2^b by Z, corrected by <= 2 add/subtracts, Y_root = q mod 2^b. Child c, sibling s: Y_c = floor((Y_v Z_s mod 2^(b_v))/2^(b_v - b_c)).
Leaf: R = ceil(z Y 2^-b) (0 if R = z) = P mod z. Proof: eta_v = (frac(P/Z_v) - Y_v 2^(-b_v) mod 1) Z_v gains < 2^(-64G) per level, so
eta < (L + 1) 2^(-64G) < 1. Cost (H-FFT): products at F(N) = 8.5 N log2(N)^2 + 1200 instructions for N limbs (+ 504 N above 2^26 limbs);
model16.gbatch_cost; [5t] replays the tree on planted leaves and checks the measured FFT points against F.

Laws (threshold_model.py; H1'). Euler-product models without the prime 2, read at n = n_f - v_2 with slack 60 - v_2: primes below 2^14 are
unsieved with all their powers, primes of [2^14, 2^56) sieved once. Yield: r_f = P(unsieved mass <= slack - rounding | B-smooth, sieved part
squarefree) (r_f = 0.9921 at n = 417 bits, slack 60). Work (H-EA): the same model over all integers bounds the survivors per cell (model16.row_cost: at slack + 32 for
the norms, at slack + 1/4 after the recheck).

Yield ([4], pool J). rels = (va/n) vb x 0.9905 x sum over J's canonical primitive region rows of (lam - qa)^3 rho(u_f) rho(u_g) r_f, with
u_s = (n_s - v_2 + 1 + kappa + alpha_s)/56: the odd part of a norm whose 2-part is known is priced as a random integer one bit larger plus
kappa = 0.07 >= kappa2(u) (sum_k 2^-(k+1) rho(u - (k - 1 - kappa)/56) = rho(u); [14v]); Dickman rho; cells with qa < 0.02 are sieved and
charged, their relations discarded. Estimate 2^50.06190 (SE 0.02214 bits) against 1.095 U' = 2^49.89091: margin 0.1709 >= max(0.06, 2 SE).
The count lowered by two SE, 2^50.0169, is the descent's input (>= 2^49.85); relation-weighted mean norms 416.93 / 470.53 bits.

Cost ([5]). Mean 47.57 instructions per cell over S's 5,136 work rows (line 9.03, record walk 25.54, survivors 3.03,
norm bounds 2.24, scan 2.11, apply 1.94, setup 1.32, routing 1.08, empty cosets 0.69, layout 0.60)
over 2^100.119 cells: collection <= 2^91.432192 units (estimate 2^91.403622 times the allowance 1.02 >= 1 + 2 SE; total SE
0.926%, [4a]).

## 4. Filtering and linear algebra

Filtering (la_bn.py, filtsim.py; [6], [6r]). Purge (to excess 0.0001 U'), clique removal and Markowitz merge are simulated on this
instance's occurrence profile (B = 2^56, mean relation norms 420 / 473 bits, above J's relation-weighted means; excess
1.095) as the tau-quotient system: tau^k(I) enters I's orbit column with coefficient chi^k (entries a + b chi, merged |a|, |b| <= 255).
x = N(400)/U' per seed (seeds fixed before the first run, none dropped): 2^17 mean 0.03733 (2 seeds), 2^18 0.03506 +- 0.00019 (8), 2^19
0.03437 +- 0.00025 (8); the mean falls with scale. H3: production is not above the two largest scales: charged NU = 0.03544 >= 0.035439, the
larger mean + 2 SE of 2^18 and 2^19 (an extrapolation of 30 bits). N <= NU U' = 2^44.9415 < 2^64/320: one owner. The attack aborts if
N > NU U' (never under H3; P_CAP = 1e-9 charged). Merge to mean weight 400, charged 401.

The square system A = [M | S] has N + 7 rows, N orbit columns and 7 dense Schirokauer columns: the chi-eigenpart of the unit groups,
(rk K_f - rk K_f^+ + rk K_g - rk K_h)/2 = (23 - 11 + 5 - 3)/2 = 7 from the signatures K_f (0, 24), K_f^+ (0, 12), K_g (0, 6), K_h (2, 2)
([6b]; tau^2 = -1 on the sigma-odd part, so its chi and chi^-1 eigenspaces have equal dimension). r is unramified in K_f and K_g ([6b]:
the absolute polynomials of degrees 48 and 12 are squarefree mod r). Block Wiedemann (m, n) = (1024, 40) with Y=[S,BZ] and a sparse random X whose seven pad coordinates occupy distinct blocks, specified in RANK-BRIDGE.md; (m, n) minimises linear algebra + generator over
m = 2^5..2^13, n = 10..200 (la_bn.bw_grid). Per row and vector: 401 x 8.38 + 176 (row tail, one owner) + 200 (chi product) = 3736.38
instructions (STORAGE.md); k5 = 8.38 >= 8.367, the measured mean + 2 SE of the entry price on the merged BN254 system ([6]). Linear algebra
<= 2^87.614018 units; generator (half-gcd matrix Berlekamp-Massey) <= 2^81.077.

Schirokauer maps. R_f = (Z/r^2)[y, x]/(h, f/36) has rank 48, R_g rank 12; eps = lcm(r^d_i - 1) < r^48 (12,192 bits; g: 3,048). The 7
columns are fixed F_r-linear forms of (phi^eps - 1)/r mod r, phi = a - b x, by left-to-right powering on eight-limb coefficients (r^2 <
2^508); per bit < 2^21.08 instructions on f, < 2^17.75 on g (schoolbook counts; the reductions by h are allowances, la_bn.schirokauer_cost):
< 2^34.69 per relation.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent, so the finite relation set and filtered matrix are fixed. H-R explicitly assumes rank M=N and rank[M|S]=N+6 for this actual finite structured set, with a one-dimensional virtual-log kernel whose map part is nonzero. H-X explicitly supplies finite right-controllability and left-observation minors. Neither follows from smooth relation counts, class-number coprimality or map injectivity. RANK-BRIDGE.md supplies the full algorithm, polynomial-degree derivation and primary-framework support. Conditional on these premises, independently sampled projection/random-source entries give failure at most4N/2^128. No probability that H-R/H-X hold is invented; they are disclosed structural applicability heuristics under exploratory review. Verification detects errors but does not establish these premises. Purged logs are reconstructed and unresolved low-occurrence ideals are handled in descent.

## 6. Individual logarithms

DESCENT.md and finalmodel4.py specify the fixed82-type production plan, original amortized node sieves and independent single-cell boot. The DESCENT.md recurrence gives expected work at most 2^95.805191 units for both individual logarithms; the ledger charges 5 times this figure, 2^98.127120 (H5x5); the hard meter t bounds the phase in every run. Each boot trial uses a fresh uniform exponent, one certified boot lattice and one sampled cell, charging exponentiation, lattice construction, complete finite cofactor schedules and all children. It does not use BALL_PAIR or a whole-ball boot waiting-time assertion.

At nodes the count premise is the conditional PGF envelope with alpha=1/2 and the stated shape. THRESHOLD-REPAIR.md combines count tails with the separate sequential availability hazard by P(fail)<=P(M<k)+psi^k; it does not use the rejected PGF(psi) substitution. It bounds the probability that no candidate relation has all usable descendants; it is a probabilistic heuristic rather than an exact count distribution. DISC and its pooled per-type confidence argument are removed. Every candidate survivor is charged; there is no negative-binomial lazy cofactor discount. Node-count and repeated-ideal aborts are explicit, and their probabilities remain part of the success calculation. There is no fixed auxiliary family or universal finite-family positivity premise.

Finite ECM curves run on the original composite modulus. All curves store stage1 prime-power-prefix Z coordinates and stage2 giant-block products on their first pass. The divisibility-chain endpoint search in DESCENT.md recovers every distinct prefix gcd, scans all stage2 block gcds and returns the same refined factor partition without curve-arithmetic replay; its deterministic storage/query costs apply to unsuccessful norms too. H-SEP lower-bounds conditional factor extraction per curve by one half of the modeled prime-order hit probability P. The configured curve count is ceil(log(0.001)/log(1-0.5P)); under independent curve coins this makes per-prime failure at most0.001. Union bounds charge every prime needed by a relation. stage2-summary.json reports real stage2/composite-modulus experiments supporting this distinction at36/40-bit factors and up to1746-bit moduli; extrapolation to production factor sizes remains explicitly heuristic, not a proof.

## 7. Success

ledger16.totals subtracts from 1: main-phase matrix failure (4N/2^128), numerical allowances, recursive descent and boot failure, node-counter and repeated-ideal aborts, split abort E_split/2^112, and the log-phase cap P(D >= t) <= E[D]/t (Markov, D the uncapped two-log work), at most 2^95.805191/t <= 0.000854 at the modeled expectation. The model gives opened-run success at least 0.995360. All losses other than the cap total at most 0.003787.

Under H5x5 (Claim) the cap loss is at most 5 x 0.000854 = 0.004270 and the opened-run success at least 0.957860 (ledger16.totals evaluates the exact product form, which is not smaller than the linear one), so with the gate overall success is at least 0.39; at the modeled values it is at least (K/2^64) x 0.995360 >= 0.405268. Every closed gate and every abort is failure, with its work charged.

## 8. Ledger

| component | log2 expected units (C=20,000) |
|---|---:|
| run-walk bases | <= 60.725727 |
| canonical a list | <= 77.387123 |
| certificates | <= 50.000001 |
| relation collection | <= 91.432192 |
| factor base and ideals | <= 60.161456 |
| matrix generator | <= 81.076063 |
| product tree | <= 55.855257 |
| linear algebra | <= 87.614018 |
| filtering/back-substitution | <= 66.305316 |
| Schirokauer maps | <= 70.588631 |
| kept-relation factorization (including doubled control budget) | <= 80.288392 |
| two individual logarithms (5 x the modeled 95.805191, H5x5; hard meter 2^106) | <= 98.127120 |
| verification | <= 11.000001 |
| wrapper | <= 53.712288 |

Only the split table row includes an extra factor two for debit/control overhead. The following underlying split kernel price precedes that extra factor; other rows retain their original allowances. Multiplying units by20,000 gives instructions. Canonical a list: 2^42.8 instructions per integer a with a_5 >= 1, Q_a(a) <= theta A^2
(2^48.87 of them; the exact test of ideal_canonical.py and the sort of the list). Factor base: 2^23.4 per prime (roots of h and f, basis
values, ideals), pi(B) < 2^51.05. Filtering streams and back-substitution: 2^24.4 per stored entry, 2^6 entries per relation, <= 2^50.19
relations. Run-walk bases: 4 reduced 7-dim bases per ideal of [2^40, 2^56) at <= 5,000,000 instructions (KERNELS.md).
Relation split: a kept relation's two norms are trial-divided by every prime of (2^20, 2^40] (< 2^35.53 primes at
156 instructions with a 2^(64 i) mod p table), each hit divided out to full power, the primes of (2^40, B] found by ECM (ecm16.py
schedule, 99.9% per prime, repeated until found or the split-only meter would be exceeded: the uncapped expectation is at most 1.002 schedules), Miller-Rabin and ideal lookups: 2^43.38 instructions per
relation in expectation. This is not a finite retry bound: the split-only meter interrupts every retry loop, and its abort loss is charged by Markov. P_g product tree: 2^70.14 instructions (H-FFT). The expected rows and the split-only and log-phase cap losses are recomputed by ledger16.totals; the non-log rows are not asserted to bound every run.

## 9. Heuristics

H1 smoothness (Dickman rho with Murphy alpha over the odd primes, exact 2-parts and kappa; exact sampled norms). H1' retention (Euler-product
model of the unsieved part). H-EA (work only): the same model prices the survivors. H-FFT (work only): F(N) bounds the multiplication work
past 2^26 limbs. H2 region (uniform pools, exact canonical predicate, measured counts with the certified bracket as fallback; sized on S,
judged on J; work with the 1.02 allowance). H-R finite relation spanning and finite H-X applicability are specified separately in RANK-BRIDGE.md; they are not implied by the following H3 assumptions. H3 matrix (filtering ratio at production at most the larger mean + 2 SE of the two largest simulated scales; Schirokauer maps injective on units modulo r-th powers; class numbers of K_f, K_g prime to r (assumed, not checked here)). H4 kernel prices (KERNELS.md: the run walk, the norm kernel and the
filtering ratio are measured at BN254 parameters; the other sieve and linear-algebra kernels are the instruction sequences measured for the
BLS12-381 instance of the same sieve, used inside their measured density ranges; the norm-bound kernel carries a 25% allowance). H-X (block
Wiedemann, section 5). H5 conditional node PGF envelope with alpha=1/2 and boot single-trial smoothness (DESCENT.md), used only through H5x5 (Claim), which carries the work bound and the success bound. H-SEP is the conditional ECM extraction-vs-order-hit lower bound and its stated extrapolation. Neither is claimed to follow from exact norm measurements alone.

## 10. Memory

Peak full-plan storage is below2^100 bytes (DESCENT.md; largest sieve2^92.436788). The main phase stays below it: a run's sum arrays (<= 1024 slices of <= 2^56.54 cells, 4 bytes per position, line-box overhead allowance x 32: < 2^73.6 bytes), the matrix (2^44.94 rows x 401 entries x 8 bytes < 2^57), the generator's sequence (< 2^66), stored relations (< 2^58). Wrapper below 2^30 bytes. All preprocessing is for this target, charged once in W.

## Provenance

The sieve, kernels, filtering and descent machinery are those of the (6,2) twisted-field STNFS packages for BLS12-381 pairing inversion by the
same author (gbotrel, with Claude in Claude Code) and their lineage on that benchmark (kaankolcu: run gate, expected total work, two-level
thresholds, norm-kernel design; mpjunior92: logarithms at their expectation; polynomials after Guillevic-Singh and Barbulescu-Duquesne;
filtering after CADO-NFS). New here: the BN254 field and pair, its pools, counts, alpha, run-walk campaign, norm kernel, filtering campaign and descent plan.

ECM-cell qualification: every ECM-only cell pays the exact same base-primitivity and ordinary predicates via counted_fraction in its PGF scale. Before any residual schedule priced at T, both exact residuals must satisfy the T bound after <=2^56 stripping; no qualifying two-window relation is lost. Perfect powers are decomposed with bounded exact root/power checks inside TDIV/boot allowances. Both final ECM gcds, including the total stage2 product gcd, remain charged. See DESCENT.md for the fixed corrected plan and full procedure.

Built on gbotrel's stnfs-bn7 (submission d10fcc8e-1fa3-430c-930a-e67b32827231, a declared dependency); one change (Claude Fable 5.2, Claude Code): a run gate, a hard logarithm-phase meter of 2^106 units, and that phase charged at 5 times its modeled expectation under the named premise H5x5 instead of at a 2^99 meter. Only ledger16.totals/gate, evidence.py's rebuild part, declared.json and the claim text differ; blobs are recompressed with identical content.

RANK-BRIDGE.md states H-R and the finite H-X bridge. ADAPTIVE.md supplies new multilevel reached-ideal evidence, including all censors and the aligned deterministic-order runs, while preserving its production extrapolation limits.

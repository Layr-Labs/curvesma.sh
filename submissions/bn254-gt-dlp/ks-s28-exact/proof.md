# bn254-gt-dlp: (6,2) STNFS on a twisted base field: expected work 2^88.390180803 units

Track bn254-gt-dlp. Input: z = g^x in GT (g = e(P0, Q0), the profile's generator), x uniform in 0..r-1; output x. The score unit is C = 20,000 scalar AArch64 instruction equivalents (a borrowed normalisation unit:
every BN254 operation below is charged in instructions). No advice, smaller targets or amortization are used. Analytic, exploratory claim.

## Claim

W is all classical work of one attack on one full-production target. Under the stated hypotheses the attack succeeds with probability at least 0.39 over uniform challenges and fresh coins (under H5 as modeled an opened run succeeds with probability at least 0.999039; the run gate opens at most a share 0.390531 of runs), and

    E[W] <= 1 + (K/2^64) 2^89.747247677 < 2^88.390180803 units (C=20,000 instructions), K = 7201150494371816446.

Run gate: before any other work the attack reads one fresh 64-bit coin u and opens the run only when u < K = 7201150494371816446, the least integer with (K/2^64) x 0.999039 >= 0.39 (evidence.py checks in exact rationals that K - 1 fails). A closed run is a failure and costs 1 unit. In an opened run both individual logarithms run under one combined phase meter t = 2^112 units: when the logarithm phase reaches t the run aborts as failure. The work spent in that phase is min(D, t), with D the uncapped two-log work, so its expectation is at most E[D] <= 2^86.318061 units by the DESCENT.md recurrence (hypothesis H5), and the ledger charges that. Both the work and the success (section 7) of an opened run rest on H5 as modeled; no allowance factor is added or claimed. The other rows are expected upper bounds under the stated hypotheses and include failed trials and abort work.

The split phase alone has an interruptible budget of 2^112 units; both individual logarithms share the 2^112-unit phase meter above. The split is deterministic (section 8, no ECM and no retry loop): before the prime-tree build and before each wave, reserve and charge that step's whole finite worst-case formula at its actual size; if the reservation would pass the split budget, abort as failure. The 128-bit remaining-budget debit and branch require at most 128 scalar instructions per invocation, charged per reservation; the ledger also doubles the whole listing bill for debit and control, including the aborting check. The field-map, collection, exact smoothness checks, linear algebra and output retain their explicit finite loops and existing cost models; their measured averages are not used as worst-case reservations.

Writing E_split for the split row, a reservation (each below 2^84 instructions) is refused only after more than 2^112 - 2^84/C units of split work, so Markov charges split abort at most E_split/(2^112 - 2^84/C). A partial or unverified output is failure.

Conventions. Bounds round in the safe direction; evidence.py recomputes every row from the shipped statistics (stats16.data, the constants of
model16.py, la_bn.py, descent16.py) and checks each declared number (declared.json) against its unrounded value; [n] cites its check lines.
What is measured at BN254 parameters and what is carried from kernels measured on another curve is listed in KERNELS.md.

## 1. Wrapper

1. Check z (canonical encoding, z != 0, z^r = 1). If z = 1 return x = 0.
2. Draw s uniform in 1..r-1 and t uniform in 0..r-1; g' = g^s, z' = z g^t (g = e(P0, Q0), the profile's generator). g' is uniform on GT minus
   {1} and z' uniform on GT, independent: the core's input law, whatever z is.
3. Map g' and z' into the core's tower (field-map.md).
4. Run the core on (g', z') under its meter and node counter (any failure is a failure): l = log_g'(z') mod r. Then x = s l - t mod r.
5. Check g^x = z (one GT exponentiation; a mismatch is a failure) and return x.

GT has prime order r, so z = g^x for a unique x mod r and z' = g'^l with l = (x + t)/s. The field isomorphism preserves z' = g'^l; the final
check rejects any wrong output. One relation and linear-algebra phase serves both logarithms (section 6). Steps 1-3 and 5 (validation, three GT exponentiations, the field map of field-map.md, randomness) are charged in the wrapper row: at most 2^68 instructions.

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
on 38,007 exact primitive pool pairs and taken at mean - 2 SE; elsewhere Kummer-Dedekind root counts (a simple root over a prime ideal of norm
q contributes q/(q^2 - 1); g has a double root over the primes dividing u: 1/(q + 1)), checked against exact valuations of 3,000 primitive
pairs ([14v]). The partial sums converge slowly (f / g: p < 2,000: (-1.261, 2.635); p < 20,000: (-0.633, 2.591); p < 200,000: (-0.485, 2.719); p < 1,000,000: (-0.469, 2.711)); the charged values (-0.449, 2.732) are the sums over p < 10^6
plus 0.02 (> their moves 0.015 / 0.009 since p < 200,000). The prime 2 is priced exactly per cell: 2 = P^3 with N(P) = 4, v_P(36) = 6 and w is a unit at P, so on primitive pairs N_f carries
2^12 exactly when b is in P and a is not, and 2^0 otherwise; N_g carries 2^2 or 2^0. All 38,007 primitive pool rows agree ([14v]: (0,0) 22,881, (0,2) 7,591, (12,0) 7,535).

## 3. Relation collection

Relations: primitive (a, b) with N_f and N_g both B-smooth, B = 2^58. Factor base: the prime ideals of relative degree 1 and norm <= B on both
sides, U' = 2^51.7080 tau-orbits.

Region. Exact rational forms Q_a, Q_b (stats16 config; place-weighted forms (1/6) sum_e om_e |x_e|^2 with om_a = |w_e|^-1/2 2^0.8,
om_b = |w_e|^1/2 2^-0.8, rounded to 10 digits: the rounded forms are the definition), A^2 = 300000. a (power coordinates): a_5 >= 1,
Q_a(a) <= theta A^2, theta = 2273/8000. b: c_5 >= 1 and Q_b(b) <= (17/20) A^2 - Q_a(a), where Q_b is the rounded form above divided by (93/100)^2 exactly (every slice's b-ellipsoid scaled by 0.93, kaankolcu's reshape), and c = SD b are the sieve coordinates (common12.SD /
SB, unimodular; the row direction c_0 is the shortest vector 1 of Q_b, the line index c_5 its shortest dual vector). The canonical a is the
(Q_a, lex)-least element with a_5 > 0 of its <+-1, units, sigma> class (exact enumeration, ideal_canonical.py); with c_5 >= 1 this is one point
per <+-1, units, tau> orbit (tau^2(a, b) = (a, -b)); non-primitive pairs are excluded by an exact HNF test. Only canonical a are sieved.

Pools ([2], [2a]). Uniform integer pairs on the region with theta_pool = 3/8 (common12.draw_a, draw_b; 100,000 pairs per seed, every draw kept):
S (seeds 101-104) sizes theta and gives the work, J (seeds 111, 112, 113, 114, 121, 122) is judged. A draw counts when a is canonical; then norms, 2-parts and
primitivity are exact. The a-draws, canonical flags and per-seed counts are those of the unscaled pools (same tag, seeds and streams); every b is redrawn from the scaled form and every S and J row recomputed. stats16 holds per-seed counts and hashes, the sums, 12 anchor rows (the last 4 from the regenerated part of S; [2c] recomputes the first two of these) and replay lists.
Counts ([3]). Points per volume of the two half-ellipsoids, by the sampler's acceptance experiment (uniform integer points of the frame box,
exact hit test; countmc.py): a: 1e+09 trials, [0.997279, 0.998215]; b: 4e+09 trials (seed 31, scaled form; our C transcription of countmc.py's protocol, which reproduces the stored 8,192-trial prefixes and a full earlier 4e9-trial count exactly), cumulative on a grid of radii^2
(lam - t/4000) A^2 covering every slice of the region (count nondecreasing, volume ~ radius^6 between grid points): [0.995954, 0.999533]. Each
point carries 2 SE (the factors scale sums over all slices). Factors 0.99324 (yield, slice cells from below) and
0.99954 (work: slice cells; a's end raised to 1 so that it also bounds a b-slice alone). The population of integer a-points is the half-volume times a's upper end A_UP = 0.998215072 (the a count above, on the draws' own theta_pool = 3/8 region); it is charged on the collection row only ([5]). Fallback: lattice_count_bounds.py certifies
(exact rationals) (0.9191, 1.0818), which contains both ([3]).
Sizing ([4s]). theta is the smallest t/8000 with margin_S >= 2 SE_S over 1.095 U' (count factor 0.99324; J, independent of S, keeps the 0.06 floor):
2273/8000 (margin 0.050834 >= 0.049370; one step below: 0.049313 < 0.049415).

Sieve (g ideals only; N_f is batch-tested; KERNELS.md). The g norm (relation-weighted 470.12 bits) is the larger one, so its survivors are the rarer ones, and f, the smaller norm, passes its smoothness test more often than g would. Per slice (one canonical a): line sieve of the lines c_5 = k for p in [2^14, 2^40); record walk for
p in [2^40, 2^57). Primes of [2^57, 2^58] are not sieved: they stay in the 112-bit slack and are found by stage B and the survivor tests below. A canonical a's run is the maximal integer t-interval containing 0 with a + t canonical and inside the limit theta A^2
(runs16.py; e_0 = 1 is the shortest vector of Q_a). On a run a prime ideal (p, s, rho)'s hits in all its slices form one coset of a 7-dim
lattice (t = rho b(s) - a(s) mod p), walked once per run and ideal pair (p, s, +-rho) (g = x^2 - u w is even in x) with one setup: (a, b) is in (p, s, -rho) iff (a, -b) is in (p, s, rho); each walk covers the full symmetric b-ellipsoid, hits of both ideals are charged (factor 2 in l) and an empty coset means no hit in the full ellipsoid. g has on average one degree-1 ideal per prime (x^2 = u w(s) mod p has two roots or none): 0.5 pairs per prime are set up, as for empty cosets. Run rule: a run of
16 to 1024 a-values all of which are in the canonical list uses the run walk in each quarter-octave band of p in which the model
density l = log2(2 x run cells x lower/upper / p) is at least 2.5, with our x-walk-asm kernel (walk7's coset, bases and records; KERNELS.md), charged per cell class at the cheaper of dadon's walk7 and rec7 prices measured on the g ideals of this sieve times its measured ratio REC7_X <= 1; runs of 2 to 15 a-values use it where l >= -0.5; other bands and runs use the slice walk (and empty cosets below
its table). 2.50% of the cells lie in runs shorter than 16; the longest sampled run has 242 a-values ([5r]).
Contamination ([5r]). The a-uniform sample of 3,000 runs at theta 2362/8000 (seeds 5001-5005, the first 100 canonical a of each inside that limit, and seeds 9101-9110, the first 250 of each; the exact predicate found no non-canonical a-value among 242,789) is restricted to its 2668 a inside the production limit 2273/8000 (uniform there); each such run at the production limit is a sub-interval of its 2362 run (canonicity does not depend on theta), so its 211845 a-values were all tested; 24 runs end next to a canonical a (runs16.run is a sufficient test: that a starts its own run, as in the work rows): contaminated runs hold under
1 - 0.05^(1/2668) = 0.00112 <= RUN_ALLOW = 0.0012 of the canonical a (95%); every cell pays RUN_ALLOW x the largest slice-walk
cost of a slice (model16.contam_sup).

Thresholds, two levels (model16.two_level, embedding_bounds.py). A cell passes when its score reaches T = lower - 16 x 112, with lower the clamped bound
max(L, U - 16 x 32) of the certified bounds L <= 16 log2 N_g <= U (Taylor shift per place in b^2, outward rounding) of the node holding the cell; where the interval of a place
meets 0 the evaluator keeps U (+16 instructions per evaluation) and takes L = 0 (a nonzero norm has |N_g| >= 1). Every kept cell's unsieved mass is at most 112 + 32 bits (charged at the exact norms), and a cell whose
16 log2 N_g lies below lower loses that difference from its slack: the pools price it per row, and it is 0 on every region row of S and J ([2x]); a zero count bounds the share of rows with any excess below 3/31238 = 0.010% (95%), less than 2% of the yield even at 100 times the mean relation weight (0.03 bit, inside the J margin above max(0.06, 2 SE)). Level 1: the rows of an
aligned 8 x 8 x 8 box share the hull of the hulls of its 2 x 2 x 4 sub-boxes; in an aligned 256-cell block each 64-cell node meeting the box is evaluated on it and is a leaf
(the clamped width is <= 32 bits); with a total width <= 32 bits the block is one leaf. Evaluations and scan leaves are counted exactly per work row ([2a]); no pool cell was left uncertified.

Survivors. Each survivor gets its exact norms (C_NORMS = 27,000, measured kernel) and is rechecked: it goes on only if 16 log2 N_g minus
its score is at most 16 x 112 (within the rounding c of the at most (n_g - v_2)/14 scored ideals, 1/16 bit each). Stage B (our BLS12-381
c27 detector on the sieved side; KERNELS.md, [5s]): m = the exact part of N_g below 2^14 (all powers), C_SB = 80,000 per stage-A survivor; the cell
goes on only if O = 16 log2 N_g - score - 16 log2 m <= 16 x 117 and O <= c + 1 or 16 x 57 k - 1 <= O <= 16 x 58 k + c + 1 (k = 1, 2). Then: (i) f batch: N_f is
B-smooth iff (P mod N_f)^(2^10) = 0 mod N_f, P the product of all primes <= B, by the exact scaled remainder tree below in waves of 2^49
norms, with 9-limb leaves since N_f < 2^521.1 < 2^576 on the region (analytic: norm_sup multiplies over the six places the bounds of |sum_k c_k w_e^k a_e^(8-2k) b_e^(2k)| from Q_a, Q_b) (C_GBATCH = 2^23 with dadon's ntt3 prices, KERNELS.md; name kept from the g-batch version); (ii) only if f passes: trial division of N_g by the 82,025 primes below 2^20 (C_TDIV = 2^23.45), a strip
with the primes of (2^20, 2^40] and one with (2^40, 2^58] (same tree; the stripped part is gcd(z, (P mod z)^(2^5) mod z); C_BATCH = 2^22.12,
C_BATCH2 = 2^23.13, the same trees priced through mul_F with the ntt3 kernel as C_GBATCH, [5t]); the pair is a relation iff the cofactor is 1. f passes with probability at most q_f = min(1, 100 rho(u_f)) per pool
row (H1; 100 is a margin). A kept relation is then factored (section 8, split row).
Lemma (exact SRT). Product tree Z_v (L levels); node v holds Y_v of b_v = 64 (n_v + G) bits (n_v limbs of Z_v, G = L + 2). Root: Newton
quotient q of P 2^b by Z, corrected by <= 2 add/subtracts, Y_root = q mod 2^b. Child c, sibling s: Y_c = floor((Y_v Z_s mod 2^(b_v))/2^(b_v - b_c)).
Leaf: R = ceil(z Y 2^-b) (0 if R = z) = P mod z. Proof: eta_v = (frac(P/Z_v) - Y_v 2^(-b_v) mod 1) Z_v gains < 2^(-64G) per level, so
eta < (L + 1) 2^(-64G) < 1. Cost (H-FFT): products at F(N) = 8.5 N log2(N)^2 + 1200 instructions for N limbs (+ 504 N above 2^26 limbs);
model16.gbatch_cost; [5t] replays the tree on planted leaves and checks the measured FFT points against F.

Laws (threshold_model.py; H1'). Euler-product models without the prime 2, read at n = n_g - v_2 with slack 112 - v_2: primes below 2^14 are
unsieved with all their powers (exact after stage B), primes of [2^14, 2^57) sieved once, primes of [2^57, 2^58] unsieved (large). Yield: r_g =
P(small + large <= slack - (c + 1)/16, large <= 117 - (c + 1)/16 - 1/8 | B-smooth, sieved and large parts squarefree) (r_g = 0.9205 at n = 471
bits); with the exact small parts of 9,856 real N_g the law's shortfall at - 2 SE is at most 0.0014 bit ([4b]); REL_LOSS = 0.0014 is taken off every relation count. Work (H-EA): the same model over all
integers bounds the survivors (slack + 32 for the norms, slack + 1/4 for stage B, then the window for the f batch); [4w]: integer toys.

Yield ([4], pool J). rels = (va/n) vb x 0.99324 x sum over J's canonical primitive region rows of (lam - qa)^3 rho(u_f) rho(u_g) r_g, with
u_s = (n_s - v_2 + 1 + kappa + alpha_s)/58: the odd part of a norm whose 2-part is known is priced as a random integer one bit larger plus
kappa = 0.07 >= kappa2(u) (sum_k 2^-(k+1) rho(u - (k - 1 - kappa)/58) = rho(u); [14v]); Dickman rho; cells with qa < 0.02 are sieved and
charged, their relations discarded. Estimate 2^51.90245 (SE 0.02182 bits) against 1.095 U' = 2^51.83890: margin 0.0635 >= max(0.06, 2 SE).
The count lowered by two SE, 2^51.8581, is the descent's input (>= 2^51.858134, its relation count LOG2_RELS); relation-weighted mean norms 413.54 / 470.12 bits.

Cost ([5]). Mean 27.51 instructions per cell over S's 12467 work rows (record walk 15.58, line 7.22, apply 1.11, scan 1.08, layout 0.63, routing 0.43, setup 0.43, survivors 0.38, norm bounds 0.37, empty cosets 0.28)
over 2^99.031 cells: collection <= 2^89.532491 units with the row cache below (region count over S and J, dadon a807d0d6; estimate 2^89.525312 times A_UP 0.998215072 times the allowance 1.009651 >= 1 + 2 SE; total SE <=
0.483%, [4a]). Parts are rounded to 0.01. Survivors, stage B included, cost 0.38 instructions per cell, so doubling C_SB would move the claim by less than 0.01 bit.

Row cache ([5c]; embedding_bounds.tail_record, model16.row_cache). NormBounds.row(s) computes each place's midpoint and radius (bm, bd) from its table tab and the five-coordinate b tail s only: tab is built from the certified root disks and SB, never from a (a enters only the place data fast, which row attaches). Once per attack, before collection, enumerate every tail of the box [-600,600] x [-626,626] x [-665,665] x [-586,586] x [1,558] (c1..c5: floor(sqrt(lam A^2 (Q_b^-1)_jj)) at the largest radius lam A^2 = 255000, widened by 7 for box_rows' 8 x 8 x 8 boxes; c5 >= 1 by the half-ellipsoid): 1,311,004,990,611,162 entries, unused ones included, each the ten binary64 words of tail_record (row's operations in their original order; no reassociation, no FMA) at the mixed-radix rank of (c5 - 1, c1 + 600, c2 + 626, c3 + 665, c4 + 586) (< 2^51; 80 B each, < 2^56.55 bytes, one owner). Build <= 2048 instructions per entry (arithmetic < 512, pages and writes 1024, enumeration) plus 2^40 = 2^61.22; the ledger charges 2^80 instructions. A row is then one lookup attaching the current a's four fast references (cached_row): 5 key loads, 5 biases, 4 radix constants, 4 multiply-adds, 2 address, 1 base load, 10 word loads and 10 stores, 4 + 4 fp loads and stores, 2 pointer, 2 call/return, 4 loop/control, 7 reserve: 64 scalar instructions instead of C_NB_ROW = 155.4 per row. Evaluations (1203.6) and the slice term are unchanged; contexts, blocks, hulls and leaves are bit-identical. A row's c0 values satisfy Q_00 c0^2 + 2 t c0 + c <= 0 at radius <= lam A^2, so it holds <= floor(2 sqrt(255000/Q_00)) + 1 = 1143 integers (Q_00 = 1126566689/1441500000) and rlen = 16 max(1, floor(count/16)) <= 1136 <= 1152: the row term C_NB_ROW/rlen falls by >= (155.4 - 64)/1152 = 0.079340 instructions per cell of the 27.5092 mean. The ledger debits only 0.984199 = (0.99324/0.99954)/1.009651 of that (count lower/upper, allowance): collection = log2(2^89.536591 (1 - 0.984199 x 0.079340/27.5092) + 2^80/20000) = 89.532490. Stored per-row work and [2c] keep the 155.4 row price. Fixture [5c]: for two anchor a's, lookups from the first a's cache equal NormBounds.row bit for bit, with the same fp objects, on 37 tails (both anchors' tails, 3 small tails, the 32 box corners); both a's tables are bit-identical; summing in reverse order changes 35 tails, so the comparison detects reordering. At other tails the identity follows from the code: tail_record performs row's operations on the same inputs.

## 4. Filtering and linear algebra

Filtering (la_bn.py, filtsim.py; [6], [6r]). Purge (to excess 0.0001 U'), clique removal and Markowitz merge are simulated on this
instance's occurrence profile (B = 2^58, mean relation norms 417 / 471 bits, above J's relation-weighted means; excess
1.095) as the tau-quotient system: tau^k(I) enters I's orbit column with coefficient chi^k (entries a + b chi, merged |a|, |b| <= 255).
x = N(200)/U' per seed (seeds 1-8 fixed before the first run, none dropped): 2^17 0.01000 +- 0.00025 (8), 2^18 0.00944 +- 0.00017 (8), 2^19
0.00844 +- 0.00017 (8), 2^20 0.00815 +- 0.00015 (8), 2^21 0.00773 +- 0.00011 (8), 2^22 0.00730 +- 0.00007 (8), 2^23 0.00687 +- 0.00005 (8), 2^24 0.00656 +- 0.00003 (8), 2^25 0.00621 +- 0.00002 (8), 2^26 0.00599 +- 0.00002 (8), 2^27 0.00569 +- 0.00002 (8), 2^28 0.00543 +- 0.00001 (8);
every scale by filtL or filt10, exact C transcriptions value-identical to filtsim, or filtfast (in law) (FILTER20.md); the mean falls at every doubling. H3: at
production the mean of N/U' is not above m* = 0.005725 x (1 + 0.0037), the larger mean + 2 SE of 2^27 and 2^28, and its variance not above v* = 7.941e-09,
the larger 7 s^2/chi2(7, 0.025) of their per-seed variances (s x sqrt(U') stays 0.25-0.51 from 2^18 to 2^28: the spread falls as U'^-1/2);
an extrapolation of 24 bits. Charged NU = 0.005747 (6 decimals) >= sqrt(m*^2 + v*), so E[N] <= NU U' and E[N^2] <= (NU U')^2. The attack aborts if N > 4
NU U' = 2^46.2650; by Cantelli P(abort) <= v*/(v* + (4 NU - m*)^2) = 2.671e-05 <= P_CAP = 2.7e-05 (rounded up), charged in success. Every run that
reaches linear algebra has N <= 4 NU U' < 2^64/320: one owner. The generator (cost per sequence term nondecreasing in N) is charged as
gen(4 NU U') x L(NU U')/L(4 NU U'), L its sequence length. Merge to mean weight 200, charged 201.
Correction (filtering simulator): the explicit small-prime band is now the actual population, one column per tau-orbit of
relative-degree-1 ideals of norm N <= 2^lsmall per side, occurrence min(1, 4 kappa/N) (filtsim.orbit_band, FILTER20.md); earlier versions
took every second prime per side, twice the declared density. Every run, k5, the table ratio and REPLAY12 are recomputed.
filtsim clamps lsmall <= 0.8 L (binding only above U' = 2^49.9: no simulated scale or replay changes).

The square system A = [M | S] has N + 7 rows, N orbit columns and 7 dense Schirokauer columns: the chi-eigenpart of the unit groups,
(rk K_f - rk K_f^+ + rk K_g - rk K_h)/2 = (23 - 11 + 5 - 3)/2 = 7 from the signatures K_f (0, 24), K_f^+ (0, 12), K_g (0, 6), K_h (2, 2)
([6b]; tau^2 = -1 on the sigma-odd part, so its chi and chi^-1 eigenspaces have equal dimension). r is unramified in K_f and K_g ([6b]:
the absolute polynomials of degrees 48 and 12 are squarefree mod r). Block Wiedemann (m, n) = (1024, 40) with Y=[S,BZ] and a sparse random X whose seven pad coordinates occupy distinct blocks, specified in RANK-BRIDGE.md; linear algebra and the generator are charged at this
fixed (m, n). Per row and vector: 201 x 7.931 x 0.8081 + 176 (row tail, one owner) + 135 (chi product, corrected truncated Shoup multiply, KERNELS.md) = 1599.22
instructions (STORAGE.md; x 0.8081: dense-band tables, KERNELS.md, [6g]; plus builds of 2^47.71 entries per slot-step); k5 = 7.931 >= max over 2^22, 2^23 of mean + 2 SE of the entry price on the merged BN254 system ([6], FILTER20.md). The chi product's constant setup (2^20 instructions) is the chi setup row. Linear algebra
<= 2^85.124044 units; generator (half-gcd matrix Berlekamp-Massey) <= 2^80.543.

Schirokauer maps. R_f = (Z/r^2)[y, x]/(h, f/36) has rank 48, R_g rank 12; eps = lcm(r^d_i - 1) < r^48 (12,192 bits; g: 3,048). The 7
columns are fixed F_r-linear forms of (phi^eps - 1)/r mod r, phi = a - b x, by left-to-right powering on eight-limb coefficients (r^2 <
2^508); per bit < 2^21.08 instructions on f, < 2^17.75 on g (schoolbook counts; the reductions by h are allowances, la_bn.schirokauer_cost):
< 2^34.69 per relation.

## 5. Main-phase success

Region and polynomials are coin- and challenge-independent, so the relation set and filtered matrix are fixed. H-R (rank M = N, rank [M|S] = N + 6, kernel the virtual-log line with nonzero map part) and H-X (rank C_d = N with the seven Schirokauer columns fixed in the start block, O_e injective on V, at the charged d = ceil(D/n) + 1 + ceil(D/(128n)), e = ceil(D/m) + 1 + ceil(D/(128m)) on the actual sparse support) are declared premises (RANK-BRIDGE.md); they do not follow from relation counts, class numbers or map injectivity. Given them, independently sampled projection entries fail with probability at most 4N/2^128. The kernel vector is certified (A(w, c) = 0, c != 0) before use; this detects, but does not bound, a failure of the premises. Purged logs are reconstructed and unresolved low-occurrence ideals are handled in descent.

## 6. Individual logarithms

DESCENT.md and finalmodel4.py specify the fixed62-type production plan, original amortized node sieves and independent single-cell boot. The DESCENT.md recurrence gives expected work at most 2^86.318061 units for both individual logarithms; the ledger charges this figure (H5); the meter 2^112 can only truncate the work. Each boot trial uses a fresh uniform exponent, one certified boot lattice and one sampled cell, charging exponentiation, lattice construction, complete finite cofactor schedules and all children. It does not use BALL_PAIR or a whole-ball boot waiting-time assertion. DESCENT.md gate G (z == 1 or z < 2^T0) is an exact lossless test (node mu, phi, psi unchanged). Gate A is a boot early abort: cells failing it are discarded and p_trial counts only qualifying cells that pass (model x 0.90). Gate E charges the boot's children at their modeled expected sizes among cells passing gate A (x 1.25) instead of the per-window maxima. Node rows keep the per-window maxima.

At nodes the count premise is the conditional PGF envelope with alpha=1/2 and the stated shape. THRESHOLD-REPAIR.md combines count tails with the separate sequential availability hazard by P(fail)<=P(M<k)+psi^k. It bounds the probability that no candidate relation has all usable descendants; it is a probabilistic heuristic rather than an exact count distribution. Every candidate survivor is charged; there is no negative-binomial lazy cofactor discount. Node-count and repeated-ideal aborts are explicit, and their probabilities remain part of the success calculation. There is no fixed auxiliary family or universal finite-family positivity premise.

Finite ECM schedules, the stored-prefix extraction and H-SEP are specified in DESCENT.md; stage2-summary.json reports real stage2/composite-modulus experiments at 36/40-bit factors and up to 1746-bit moduli, whose extrapolation to production factor sizes is H-SEP.

## 7. Success

For an opened run, ledger16.totals subtracts from 1: main-phase matrix failure (4N/2^128), numerical allowances, recursive descent and boot failure, node-counter and repeated-ideal aborts, split abort E_split/(2^112 - 2^84/C), the filtering abort P_CAP, and the log-phase cap. For the cap: the two logarithms are computed for the wrapper's independent uniform targets g' and z' with separate fresh coins, and the second uses nothing computed by the first, so their uncapped works D1, D2 are independent; the DESCENT.md recurrence is per logarithm (2^86.318061 is twice the one-log bound), so E[D_i] <= 2^85.318267. The meter is reached only if D1 + D2 >= t = 2^112. By the lemma below with a = 2^(85.318267-112) <= 0.000000010: P(cap) <= max(2a - a^2, a/(1 - a)) <= 0.000000019. All modeled losses other than the cap total at most 0.0009601. So the model gives an opened run success at least 0.999039.

The gate uses this modeled figure: the gate coin is independent of the run, so the attack succeeds with probability at least (K/2^64) x 0.999039 >= 0.39. This figure and the charged work both rest on H5 as modeled. A closed run and every abort is failure, with its work charged.

Lemma (two-variable case of Samuels' inequality, 1966; the form and proof follow kaankolcu's accepted GT package 8c816a9a). X, Y >= 0 independent, E X <= a t, E Y <= a t, a <= 1/2. Then P(X + Y >= t) <= M(a) = max(2a - a^2, a/(1 - a)).
Proof. Scale to t = 1 and replace X, Y by min(X, 1), min(Y, 1): the event and the independence are unchanged and the means do not increase. G(x) = P(Y >= 1 - x) is nondecreasing on [0, 1], so its least concave majorant G^ is nondecreasing and P(X + Y >= 1) = E G(X) <= E G^(X) <= G^(E X) <= G^(a) (Jensen). G^(a) is the supremum over 0 <= x1 <= a <= x2 <= 1 of (1 - l) b1 + l b2 with d = x2 - x1, l = (a - x1)/d, b1 = G(x1), b2 = G(x2) (if d = 0, G(a) <= a/(1 - a) by Markov). Here 0 <= b1 <= b2 <= 1, and Y >= (1 - x1) 1[Y >= 1 - x1] + (1 - x2) 1[1 - x2 <= Y < 1 - x1] gives, taking expectations, d b1 + (1 - x2) b2 <= a. The linear objective is largest at a vertex of this polygon. (i) b1 = 0: the value is l b2 <= (a/x2) min(1, a/(1 - x2)) (l decreases in x1); for x2 <= 1 - a this is a^2/(x2(1 - x2)) <= a/(1 - a) since x2(1 - x2) >= a(1 - a) on [a, 1 - a]; for x2 >= 1 - a it is <= a/x2 <= a/(1 - a). (ii) b1 = b2 = a/(1 - x1) <= a/(1 - a). (iii) b2 = 1, b1 = (a - 1 + x2)/d, which needs x2 >= 1 - a: with c = x1 and v = 1/d the value is 1 - (1 - a - c) v + (a - c)(1 - a - c) v^2, convex in v (c <= a <= 1/2), so it is largest at an end of x2 in [1 - a, 1]. At x2 = 1 - a it is (a - c)/(1 - a - c) <= a/(1 - a); at x2 = 1 it is psi(c) = (2a - a^2 - (1 + a) c + c^2)/(1 - c)^2, and psi'(c) has the sign of (1 - a)(c - (1 - 2a)), so on [0, a] psi <= max(psi(0), psi(a)) = max(2a - a^2, a/(1 - a)). Every vertex value is at most M(a). QED

## 8. Ledger

| component | log2 expected units (C=20,000) |
|---|---:|
| run-walk bases | <= 70.645293 |
| canonical a list | <= 76.664538 |
| certificates | <= 50.000001 |
| relation collection | <= 89.532491 |
| chi product setup | <= 5.712288 |
| factor base and ideals | <= 62.110830 |
| matrix generator | <= 80.542367 |
| product tree | <= 58.017996 |
| linear algebra | <= 85.124044 |
| filtering/back-substitution | <= 68.262253 |
| Schirokauer maps | <= 72.545567 |
| kept-relation factorization (deterministic listing, doubled bill below a fixed allowance) | <= 81.403537 |
| two individual logarithms (expected work under H5; meter 2^112) | <= 86.318061 |
| all rows of an opened run together | <= 89.747247677 |
| after the run gate, plus 1 unit | 88.390180803 |
| verification | <= 11.000001 |
| wrapper | <= 53.712288 |

Only the split row doubles its bill for debit/control overhead (the bill below precedes that factor); other rows retain their original allowances. Multiplying units by20,000 gives instructions. Canonical a list: 2^42.8 instructions per integer a with a_5 >= 1, Q_a(a) <= theta A^2
(2^48.17 of them; the exact test of ideal_canonical.py and the sort of the list). Factor base: 2^23.4 per prime (roots of h and f, basis
values, ideals), pi(B) < 2^53.00. Filtering streams and back-substitution: 2^24.4 per stored entry, 2^6 entries per relation, <= 2^52.20
relations. Run-walk bases: 5 reduced 7-dim bases per ideal of [2^40, 2^58) for walk7 at <= 1.0e9 instructions (KERNELS.md).
Relation split (deterministic two-tree factor listing; no ECM, nothing probabilistic), in waves of Q = 2^40 pairs, the last partial. (a) Both norms (nonzero, < 2^576, exactly B-smooth, section 3) are trial-divided by the primes <= 2^20 with full powers (2^26 per pair); <= 28 prime factors with multiplicity remain per norm. (b) A balanced product tree of the primes of (2^20, B] (pi(B) < 2^53, height <= 53) is built for this phase and charged (fb, gtree rows kept); rad(z) = gcd(z, P_root mod z) by one remainder pass per wave. (c) At a prime-tree node each record (id, d), d > 1 squarefree, goes left as g = gcd(d, P_left mod d) if g > 1 and right as d/g if d/g > 1; the children are coprime, so d's support is partitioned exactly however norms share primes. At a leaf d = p: divide the ORIGINAL norm by p to full power, then the original ideal qualification, lookup and reconstruction. A norm not reproduced exactly by its listed powers is failure, never called prime. Two-word ids in stable depth-first append lists keep input order without a sort.
Bill (model16.split_parts, instructions, affine in the expected pair count R = 2^rels_upper and E[ceil(R/Q)] <= R/Q + 1): D = 54 passes (53 levels and the radical), <= 56 records per pair and pass. Per record and pass, its share of a padded nine-limb norm product and remainder tree (n records, m < 2n leaves, L <= 46, G = 48): 2(C(36L + 2G) + 3600) + 2^15, C = 2^15 + 256 per limb (F(N) <= 2^15 N + 1200 for N <= 2^55, [8s]; 256 for copies, addressing). Roots per active node and wave, even at n = 1: six products of 2 SP_v + 136 n limbs, prime-node limbs < 2^54 per depth: 6C(2^55 D E[W] + 136 * 56 R D) + 7200 * 56 R D, paying P_left against sparse groups and partial waves. Plus per record and pass 2^24 (nine-limb gcd, exact division) and 2^16 (queues, ids, control); per prime occurrence (<= 1152 per pair) 2^20 (powers, output) and the original 2^26 ideal qualification; D 2^54 C + 2400 2^53 for the prime tree; 128 per reservation. Total < 2^77.60 units, doubled 2^78.60; the row is the larger of that and a fixed 2^81.403537 (half the former doubled trial-division/ECM row). One full wave reserves < 2^80.2, the prime tree < 2^74.8 instructions. New storage (54 2^54 8-byte prime layers, 32 FFT buffers <= 2^55 words, norm trees, queues, a wave's outputs) is < 2^64 bytes, within the 2^100-byte allowance. [8s] lists 256 toy norms over primes < 512 (powers, shared primes, units, partial waves) exactly against trial division and rejects a radical taken as the factorization and an unlisted prime: toy scope only. P_g product tree: 2^72.31 instructions (H-FFT). The expected rows and the split-only and log-phase cap losses are recomputed by ledger16.totals; the non-log rows are not asserted to bound every run.

## 9. Heuristics

H1 smoothness (Dickman rho with Murphy alpha over the odd primes, exact 2-parts and kappa; exact sampled norms). H1' retention (Euler-product
model of the unsieved part and of the stage-B window; [4b], [4w]). H-EA (work only): the same model prices the survivors (also those of descent gates G and A). H-FFT (work only): F(N) bounds the multiplication work
past 2^26 limbs (also in the split trees). H2 region (uniform pools, exact canonical predicate, measured counts with the certified bracket as fallback; sized on S,
judged on J; work with the allowance 1.009651 >= 1 + 2 SE (count over S, J and the count-only draws)). H-R finite relation spanning and finite H-X applicability are specified separately in RANK-BRIDGE.md; they are not implied by the following H3 assumptions. H3 matrix (filtering ratio at production with mean at most the larger mean + 2 SE and variance at most the larger 97.5% variance bound of the two largest simulated scales, tabled entry price ratio at most TAB_RATIO ([6g]); Schirokauer maps injective on units modulo r-th powers; class numbers of K_f, K_g prime to r (assumed, not checked here)). H4 kernel prices (KERNELS.md: the run walk, the table builds, the norm kernel, stage B and the
filtering ratio are measured at BN254 parameters; the other sieve and linear-algebra kernels are the instruction sequences measured for the
BLS12-381 instance of the same sieve, used inside their measured density ranges, except that LINE_TAB keeps its densest measured price above its top step (per-hit price does not rise with line density; every measured run lies at or below its step, linex.txt); x-walk-asm is measured against those prices at theta 2280/8000, carried per density stratum (REC7_X, [10x]); walk7 is measured at BN254 on the g ideals it serves at theta 2189/8000 on the unscaled b form and carried per density stratum to this region (KERNELS.md: a walk enumerates one ideal's 7-dimensional coset, whose cost per hit is set by the coset density l and the ideal's reduced basis, both read per stratum; the region shape changes the mix of l, empty cosets and setups, which the model charges per cell); the norm-bound kernel is charged at its BN254 measurement, its unmeasured slice term with a 25% allowance). H-X (block
Wiedemann, section 5). H5 (expected work and success of the logarithm phase) conditional node PGF envelope with alpha=1/2 and boot single-trial smoothness, including gate A's pass probability (window composition restricted to <= tau bits above 2^x, x 0.90) and gate E's expected child counts (x 1.25) (DESCENT.md; its conditional next-relation premise is measured in ADAPTIVE.md's conditional table: after an earlier failed relation the ratio of after-failure to first-attempt observed/envelope failure is 0.964, 95% 0.839-1.095, and conditional child work 0.950 (95% up to 0.983) of the envelope). H-SEP is the conditional ECM extraction-vs-order-hit lower bound and its stated extrapolation. Neither is claimed to follow from exact norm measurements alone.

## 10. Memory

Peak full-plan storage is below2^100 bytes (DESCENT.md; largest sieve 2^85.417380). The main phase stays below it: a run's sum arrays (<= 1024 slices of <= 2^56.54 cells, 4 bytes per position, line-box overhead allowance x 32: < 2^73.6 bytes), the matrix (2^45.79 rows x 201 entries x 8 bytes < 2^56.45), the dense-band tables (8 slots x 2^55.04 bytes, [6g]), the generator's sequence (< 2^66), stored relations (< 2^58), the b-tail row cache (< 2^56.55, section 3). Wrapper below 2^30 bytes. All preprocessing is for this target, charged once in W.

## Provenance

Line (polynomials after Guillevic-Singh and Barbulescu-Duquesne; filtering after CADO-NFS): our g-side stnfs-bn10gs (858c0eb4) and stnfs-bn10m (455c5dff / c48a42b4: dense-band tables, our BLS12-381 c27 detector cfd16473 on the g side with the c29q2 clamp 2952b2e4, sizing at 2 SE, batch descent strips); kaankolcu: exact strip (a2b9e45e), norm-bound kernel (9d0f9081), ntt3 strips and H-X support (e8fbf83a), b-ellipsoid / 0.93^2 (879d6db3; PI 8785cc1c); kshitij-hash: plan (fc2c0a8a), rec7 on g ideals (64c3b342), two-pool count (3a9e85be), NU from 2^23 (713c8ea2); dadon: ntt3 radix 4, walk7 on g ideals and its faster builds (09c32104, 86fb9e54, 51a9357a), count-only region draws (17332cc8); 0xmoei: allowance at 1 + 2 SE (86a3751f). Base 3e250f02 by gbotrel (gbotrel; Claude Opus 5.5, Claude Code): the reshape on our own regenerated pools (b count seed 31, alpha, runs, contamination; theta 2273/8000), 4,000,000 count-only draws (seeds 2001-2040), LOG2_RELS 51.858, plan re-optimised, line/scan/staddh-apply prices (KERNELS.md), merge weight 200 at 2^22/2^23, gate G unit residual ([7c]); then our hand-written run walk (rec7xa.txt, REC7_X per stratum), stage-A slack 112 (its [4b] relation-count loss charged, REL_LOSS), a contamination sample of 3,000 runs, line steps to 2^25.11 and block-segment scan prices, ntt3 only up to 2^55 points (KERNELS.md). Collection 2^89.539169, linear algebra 2^85.919409, two logarithms 2^86.363717, claim 88.448033187. Base f645b7a8 by soispoke (kshitij-hash): filtering scales 2^24, 2^25, 2^26, 2^27 (filt10 / filtfast, FILTER20.md), NU 0.00604. Claim 88.414006755. common12's tower arithmetic sits in the lib as tower12.py (moved verbatim to fit the archive cap; evidence output unchanged). This rebase (soispoke): all package data in one base126 lzma lib (lib.b126, libextra.b126; data blobs as .data sections, common12.BLOB), the chi product as a corrected truncated Shoup multiply (C_CHI 135, KERNELS.md; setup row 2^5.712288), the descent node-failure floor lowered from .02 to .01 (THRESHOLD-REPAIR.md), the kept-relation split as a deterministic two-tree factor listing instead of trial division and ECM (section 8; split row 2^81.403537), and the collection population charged at the measured a-point density A_UP = 0.998215072 instead of 1 (section 3; collection 2^89.536592): linear algebra 2^85.272620, two logarithms 2^86.318267, success 0.998641, claim 88.400765970. Then the b-tail row cache (section 3, [5c]): collection 2^89.532491, claim 88.397247328. Base PI 0c8b760c (kshitij-hash): filtering scale 2^28 (filtfast, 8 seeds) and exact charges: NU 0.005747, P_CAP 2.7e-05 (Cantelli bound), C_ENTRY_K 7.931 (2^22, 2^23), LOG2_RELS 51.858134, node budget grid 0.001. Claim 88.390180803. This version (kshitij-hash; Claude Opus 5.5, Claude Code): the bn254-gt-dlp edition of our PI 0c8b760c, made by replaying gbotrel's accepted PI-to-GT edition differences (1e237e96 against 3e250f02): track, input and output (x), steps 1 and 5, the wrapper sentence, no output row, GT submission ids in this line. No figure changes: claim 88.390180803.
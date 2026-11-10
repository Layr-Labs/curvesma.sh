# Descent

descent16.py is the model (plan, ball counts, sub-bins in plan16.b85); [7] evaluates it, [7c], [11], [11b], [12] check inputs. (a, b) in O_h^2 is a vector in Z^12 (basis e_j).

## Lattices

descent_witness.reduce LLL-reduces (delta 0.99) under the integer block form Q(a, b) = round(2^24 sum_s |a(s)|^2) + round(2^24 lam sum_s |w(s)|^e |b(s)|^2), s over the six complex embeddings of K_h. Nodes use QN (lam 3/2, e 1/2), the boot QB (lam 2^64, e 1). Regions: radius-R balls in reduced-basis coordinates, one cell per +- pair.

## Targets and boot

T becomes T' = T^e, e uniform in [1, r), log T = log T'/e mod r (2/r charged). The boot lattice {(a, b) in O_h^2 : a B + b A = 0 in F_p12}, det p^6, has a certified basis (rows in the lattice, Gram determinant p^12, size and Lovasz under QB in exact rationals; [11]); the boot element's single (g-side) norm is N(a^2 + U w b^2). The boot norms each radius-R ball cell (Envelopes) and runs ecm16.py's schedule: a relation has every prime <= 2^n1, <= kA in (2^69, 2^nA], <= 1 in (2^nA, 2^n1]. An exhausted ball draws a new exponent and certified lattice (2^36 each). UDUP = 0.65 as at nodes.

## Nodes

A node is a special-q ideal (l bits, side f/g; (q, y - s) and a root t of the side polynomial at s) with reduced lattice L_q (det q) and a radius-R ball ([7] recomputes the counts R <= 8). Two node types:

- sieve: dsieve16 enumerates every p in [2^8, 2^69] and relative-degree-one root pair on the ball into owners' score arrays. Each primitive cell gets certified upper bounds on both norms (outward-rounded embeddings, charged 2,400); cells above the node maxima are skipped (Envelopes); others with residual mass (log2 bound - score - log2 q on the special-q side) <= T = kA nA + kB n + 96 on both sides survive (exact norms, trial division, ECM); non-primitive cells repeat c/gcd(c)'s relation. Certificate: the bound pass sums ss_f ss_g at kept cells' bounds (q removed); a node below cells x charged pr is discarded before cofactorisation, so every used node's yield is >= charged.
- ecm: ECM on every ball point (never cheaper here).

A relation has every residual prime <= 2^n per side, <= kA in (2^69, 2^nA] and kB <= 1 in (2^nA, 2^n]; these window primes are the children (child load, below).

## Yield and survivors

Per cell, the two-window semismooth probability

P = sum_{jA <= kA, jB <= kB} 1/(jA! jB!) int rho((m - sum t + alpha)/69) prod dt/t

is a lower sum (exact int dt/t per mesh cell, rho at its left endpoint, cells crossing sum t = m dropped); alpha = (0.29, 0.75). Norm sizes m come from offset bins over 61,440 exact norms of uniform ball points on 160 special-q lattices (10 sizes 72-300 bits, both sides, radii 8-128, 64 points each): cf = log2 N_f - 72c, cg = log2 N_g - 12c - 6 log2 U, c = (log2 q + log2 cells)/12. cf sextile maxima pair with the overall cg maximum, shifted up by the worst per-(ideal, radius) mean minus the pooled mean (+14 cf, +5 cg). Boot norms: 3,072 exact norms on 8 certified boot lattices (p90, p99, max + 2, the worst per-target mean minus the pooled mean). For pr the sextiles are subdivided (QBINS_SUB: 24 cf-quantile cells with cg maxima per cell, same shifts; BOOT_SUB: boot deciles + shift); each term stays a lower sum. Charged pr = min(sub-cell mean, 4 x sextile pr) at nodes, min(sub-cell mean, 1.5 x sextile pr) at the boot. mu = cells x P x UDUP x YF x detection x retention (detection 1 - 0.001 x expected factors; the boot carries the counted fraction; YF = 0.8 charges the toy's measured relations/certificate ratio, worst type 0.82, below). Survivors per cell use (T/69) rho((m - T)/69) per side, decreasing in m, at each bin's lowest norms (previous sextile's raw cf max, observed cg min, minus shifts; special-q at the bin's upper edge).

## Factor-base ideals without a logarithm

Columns are tau-orbits of 4 ideals. At the descent input of 2^62.583 relations (G + H's 2-SE count) an ideal of norm q = 2^(69 t) on side s occurs in a relation w.p. kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/rho(u_s), u_s from mean norms 740.0/498.4; a column occurs Poisson(lambda), lambda = 2^62.583 x 4 kappa_s/q, and a column occurring <= 2x is missing. A relation of an m-bit side carries sum_b kappa(m/69, b/69) ln((b+0.1)/b) P(Poisson <= 2) missing ideals on average (1/10-bit bands b < 69). A node discards (relations with a missing ideal rejected, mu x exp(-missing); always at q <= 2^85) or descends (each missing Q sieved in the lattice of Q q', q' a fresh ~2^20 known-log ideal, charged as the worst discarding node with q <= 2^90); the boot descends.

## Retention and cofactorisation

Counted cells: base-primitive only (no prime P of O_h contains both a and b); others repeat (a - b theta)/pi's relation (h_K = 1). Base-primitive fraction prod (1 - N(P)^-2) = 0.6883 (exact for p < 2^9; N(P2) = 2). Cells divisible by a relative-degree-one ideal above a base prime of degree >= 2 with p >= 2^8 are excluded by ordinary fractions 0.9514/0.9917 (f/g, disjoint). Factor base: all prime ideals of norm <= 2^69 above p < 2^8, degree-one ideals above p in [2^8, 2^69].

Omitted mass. A counted cell is rejected only if one side's omitted factor-base mass X exceeds s = 96 - 1 (bound gap < 1e-9) - m/128 (score rounding); X sums independent local terms (whole valuation at p < 2^8, excess powers above), each with one local bound. Odd p: primes (p, h_i), residue field F_q (O_h = Z[y] locally). For side polynomial S and a base-primitive pair, b is a unit w.p. q/(q + 1), else S_hom(a, b) is; then v_P(S_hom(a, b)) = v_P(S(t)), t = a/b uniform in O_P, 0 unless t reduces to a root rho of S in F_q. Hensel factors S = G H over O_P with G = (x - rho)^mu mod P and H(t) a unit for t = rho mod P; over a splitting field G(t) = prod (t - alpha_i), i <= mu, so v_P(S(t)) >= j forces some v_P(t - alpha_i) >= j/mu (v_P extended); for each i those t in O_P lie in one class modulo P^ceil(j/mu):

P(v_P >= j and t reduces to rho) <= min(1/q, mu q^-ceil(j/mu)),

with equality q^-j when mu = 1. At 2 (exact divisions by w and pi4): v_P2(F) = 6 when a is in P2 and b a unit (Newton polygon of the f coefficients' valuations 0, 3, 5, 5, 6; probability 1/3), else 0; F has no root in F_4, so v_P4(F) = 0; g is x^2 + U w with a double root 0 at both primes (the bound above with mu = 2). Primes above p are independent (CRT). The local structure at p < 2^9 is tabulated (retain16.structure); [11b] recomputes part of it and checks the bound on residue-class samples at 2 and 19. Above 2^9 excess powers are bounded per dyadic block by 72 p^-2 events per level over pi(x) <= 1.26 x/ln x primes: k events of total multiplicity t have P <= lam^k/k! C(t-1, k-1) 2^(-j(t-k)) (every k, t <= 64; the rest, < sum t (2 max(lam, 2^-j))^t/2, at infinite mass). loss reads a safe summary of X_up (stats16 omitted_law: P(X <= lo) rounded down; masses at lo+1..hi and per 64-bit block, >= 40 merged, rounded up; evidence16 [11b] recomputes it). X_up rounds masses up and puts negligible local tails at infinite mass (counted lost): means 8.49 (f) and 8.45 (g) bits, P(X > 96) = 6.1e-11 (f) and 3.3e-10 (g) per cell.

Relation-conditioned loss. Under H1 a cell with omitted mass x is semismooth with probability Pss(m - x), nonincreasing in m. X_up dominates the local law X (each independent local term only moves up). Among relation cells the lost fraction is L/(L + K) <= L/K, L = E[1(X > s) Pss(m - X)], K = E[1(X <= s) Pss(m - X)]. (i) coupling X <= X_up: 1(X > s) Pss(m - X) <= 1(X_up > s) Pss_hi(m - c64(X_up)) (c64: up to the 64-bit grid), so L <= its X_up mean. (ii) Pss(m - x) >= Pss(m) >= Pss_lo(m) for x >= 0: K >= Pss_lo(m) P(X <= s) >= Pss_lo(m) P(X_up <= s). Hence

loss <= sum_{x > s} P_up(x) Pss_hi(m - c64(x)) / (Pss_lo(m) sum_{x <= s} P_up(x))

at the node's largest norms m_f, m_g (retain16.loss). Retention = (1 - loss_f - loss_g) x base-primitive x (1 - ordinary exclusions); [7] evaluates every sieve node (min 0.6490; largest one-side loss < 5e-5).

ECM (ecm16.py). Suyama Montgomery curves, cheapest (B1, B2/B1) point with wheel 30030; P1 = the Dickman lower sum on the curve order; curves = ceil(ln 0.001 / ln(1 - P1)) at 367.76 (L/5)^2 per mult. [12] prints the schedules and checks P1 on random primes.

## Model of one attempt (a prediction; the claim uses only H5, proof.md section 6)

Distinct usable relations per node are negative binomial (mean mu, shape r): r = 0.5 for effective windows of <= 1 prime (base nodes (bin (69, 75)) have empty windows), r = 1 otherwise (model; toy_descent.py compares actual counts).

A subtree fails with probability phi(node) = max((1 + mu (1 - psi)/r)^-r + c DISC_t + psi^32, 0.02); node work = own sieve/ECM + sub/(1 - psi); c: Encountered law, below.

Child load (union form). A relation's children are its window primes c_i. Union bound (no independence within a relation): P(fail) <= psi := E[sum phi(c_i)], work <= sub := E[sum W(c_i)]; tries stop at the first success (model): Wald, work <= sub/(1 - psi). E is H1's law given semismoothness (yield integrand weighted by sum g(t_i), g = phi or W, max over bins meeting the mesh cell). child_load: upper sums at norm m over lower sums at m + 1 (<= P(m), P nonincreasing in the norm) bound the mean at m = ceil(node maxima; boot ceil(m_max)); no cofactored norm exceeds m (Envelopes); charged x LOAD = 1.001. Smaller norms are assumed not to load more (model). counts16 keeps the slot-count form. The plan needs mu >= 18, phi <= 0.07, psi <= 0.97 ([7]). DISC_t: the one-sided 95% bound on type t's shortfall rate in PREREG-E (Supporting measurements).

Conditional laws. The retry forms (P(J > j) <= psi^j, work sub/(1 - psi), boot psi^100) take the failures of different relations' children as independent. The toy descent does not bear this out: there a retry fails after an earlier failure far more often than a first try (figures below). The claim therefore does not use these laws: success rests on restarts with fresh exponents, which are exactly independent (proof.md section 6), and this model only predicts one attempt's failure probability and work.

Encountered law (argument under H1). A degree-one prime ideal I = (q, s, t) divides a cell's norm iff the cell lies in the index-q lattice L_I: probability 1/N(I), the same for every degree-one ideal of that norm. So the children in a bin are ideal-uniform with density (1/b) x (probability that the remaining cofactor is semismooth) in the bit size b. PREREG-E samples ideal-uniform with b uniform on the bin; dom_c bounds the density ratio by its max/mean over a bin, c <= 1.198 at the boot's norm (the largest), charged 1.2 ([7c]).

Toy (toy_descent.py, declared experiment; specification in its header): two-polynomial NFS over Q (f = x^3 + x + 1, g = x - m, p = f(m) of 21, 24, 30, 36, 39, 42 bits), real special-q lattices and relations, the package's formulas (ss_lower2, child_load, the phi and W recurrence) with Murphy alpha. An attempt: a uniform target, its boot relation, the full adaptive descent, capped at c x the model's mean work. 2,000 attempts each (inputs {"pb": 22}, {"pb": 26}, {"pb": 30}, {"size": 0}, {"size": 1}, {"size": 2}, with "attempts": 2000): failed attempts 8, 2, 0, 0, 0, 0 at c = 2 (95% bound <= 0.0072; model 0.77 to 0.91), 0 in every run at c = 4, 8, 16; 0 back-to-back failures. Distinct relations / certificate sum 0.91, 0.91, 0.92, 0.93, 0.93, 0.92 (worst type 0.81): below 1, charged as YF = 0.8. Nodes without a relation 14, 26, 48, 74, 56, 80 against Poisson 18, 32, 68, 118, 92, 114, negative binomial 126, 275, 636, 1056, 1111, 1288. Retries: a first relation failed in 54 of 12,920 nodes, a second after a failed first in 1 of 53; an earlier tree-rooted version (recorded, R = 3) saw 8 of 49 (16%) against 2.2%: failures cluster on weak nodes, so the within-tree retry law is not conservative and the claim does not use it. Limits: rational NFS, no tower, no unit duplicates, at most one window prime per side, tiny sizes.

Plan: 62 nodes (bin (69, 75) and 5-bit bins 75-225, both sides) and the boot (n1 = 215, nA = 184.58, kA = 15, R = 4); entries (n, iA, kA, kB, R, discard) from a coordinate search under the three conditions: min mu 33.7, max phi 0.0510, max psi 0.831 (boot; nodes 0.533); base nodes R = 224/256 (mu 520/469); all sieve (R = 128..384). Meter, Markov, boot failure: proof.md sections 6-7.

## Lazy cofactorisation and boot

Relations are tried singly; the first whose subtree succeeds ends the node. Nodes: after the sieve (charged in full) the s survivors are cofactored in random order, each distinct relation descended at once, until one's children all succeed (failure when survivors run out). With M distinct relations among s survivors the k-th first appears at expected position <= k (s + 1)/(M + 1); the number J tried is independent of the order (model), P(J > j) <= psi^j, E[min(J, M)] <= 1/(1 - psi):

E[processed]/s <= E_M min(1, (1 + 1/s)/((1 - psi)(M + 1)) + psi^M),

M negative binomial; descent16.lazy_fraction sums it; own work = LLL + sieve + survivors x cost x (fraction + c DISC_t (1 - fraction)): a node short of distinct relations is charged its full cofactorisation. Boot: balls are i.i.d. (fresh exponent and lattice, kept cells tried in full); a ball holds a distinct relation w.p. >= l(1 - l/2), l = cells/per_rel = 2^-30.41 (Bonferroni, H1 pairwise), so cells per relation <= per_rel (1 + 2^-31); <= 100 relations, one at a time (stopping time over balls; Wald): trial work per_rel (1 - psi^100)/(1 - psi) x trial + one ball (2^89.1); failure psi^100.

## Envelopes (enforced)

Nodes cofactor a primitive cell only if its certified bounds give log2 N_f <= nf, log2 N_g <= ng (maxima, q removed).
The certificate sums kept cells only (so does PREREG-E's raw); retention, ECM cost, detection, missing ideals and child_load at the maxima hold for every cofactored cell. Survivors per relation are taken as the model's (sv + pr)/pr (assumed monotone in the norms; model). At
most JMAX = 32 relations per node (psi^32 in phi). Boot: exact norms; cells above m_max skipped; a ball (exponent,
lattice) is discarded unless its kept cells' H1 yield sum reaches cells x pr, so trials per relation <= 1/pr; discards
charged at BOOT_DISC = 0.9 (< 2^86 instructions per log). Child loads are tabulated (plan16, rounded up); [7t]
recomputes all 242 (>= fresh).

## Supporting measurements

[11] re-samples a boot and a 260-bit f special-q lattice against the bins. Unit duplicates: a unit eps maps ball cell v to eps v, same relation (tau moves a node's lattice to tau(Q)'s, the boot's to Frob^3(z)'s); O(v) = {eps : eps v in the ball} (exact: unit_orbit.py, Fincke-Pohst, norm +-1), m = |O(v)|/2, and mu needs the distinct yield E[y/m] >= UDUP pr.

PREREG-E (descent_witness.camp_ideal; seed 30000 + 1000 i + j for ideal j of type i in ['boot'] + sorted plan, fixed in advance; boot 20000 + j). enc_ideal: b uniform on the bin, a random prime q of that size, a root s of h mod q kept w.p. (roots)/6, a root t of the side polynomial at s kept w.p. (roots)/12 (g: /2): every degree-one ideal has the same chance; LLL basis as in the descent (boot: fresh exponent, certified lattice). K = 256 uniform ball cells (boot 512; ball_fast: a uniform point of the radius R + sqrt(12)/2 ball rounded, so exactly uniform; [11] chi^2); per cell the certificate integrand y (ss_f ss_g at certified bounds, q removed, 0 above the maxima) and the exact unit multiplicity m. raw = mean y/pr, yr = mean (y/m)/pr (pr charged), raw_lo, yr_lo without the largest cell; shortfall: yr_lo < 2 UDUP or raw_lo < 2 (the 2 absorbs K-cell noise). Result: 23,312 ideals (n = 1,000 per base type, bins 69-90; 256 elsewhere), 0 shortfalls, least yr_lo 1.46 (i = 62): DISC_t = 1 - 0.05^(1/n) = 0.00300 (base), 0.01164. Boot (64 x 512): yr 0.944, mean - 2 SE 0.935 >= 0.65; raw 1.139. plan16 keeps n, shortfalls, worst ideal; audit {"ideal": i} replays it.

Not shown: any relation count, subtree outcome or complete attempt at production size (about 2^-84 per cell; no production node is computable). H5 rests on this model's prediction (p_model <= 0.0636) and on the toy.

Refined integration: mesh 384, rho step 1/1024 (lower sums); audit_refinement.py recomputes four probes (0 < lower <= upper); iA = 8 x physical windows.

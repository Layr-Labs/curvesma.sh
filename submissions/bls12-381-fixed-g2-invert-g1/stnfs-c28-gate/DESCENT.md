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

is a lower sum (exact int dt/t per mesh cell, rho at its left endpoint, cells crossing sum t = m dropped); alpha = (0.29, 0.75). Norm sizes m come from offset bins over 61,440 exact norms of uniform ball points on 160 special-q lattices (10 sizes 72-300 bits, both sides, radii 8-128, 64 points each): cf = log2 N_f - 72c, cg = log2 N_g - 12c - 6 log2 U, c = (log2 q + log2 cells)/12. cf sextile maxima pair with the overall cg maximum, shifted up by the worst per-(ideal, radius) mean minus the pooled mean (+14 cf, +5 cg). Boot norms: 3,072 exact norms on 8 certified boot lattices (p90, p99, max + 2, the worst per-target mean minus the pooled mean). For pr the sextiles are subdivided (QBINS_SUB: 24 cf-quantile cells with cg maxima per cell, same shifts; BOOT_SUB: boot deciles + shift); each term stays a lower sum. Charged pr = min(sub-cell mean, 4 x sextile pr) at nodes, min(sub-cell mean, 1.5 x sextile pr) at the boot. mu = cells x P x UDUP x YF x detection x retention (detection 1 - 0.001 x expected factors; the boot carries the counted fraction; YF = 0.43: below the toy's worst measured relations/certificate ratio). Survivors per cell use (T/69) rho((m - T)/69) per side, decreasing in m, at each bin's lowest norms (previous sextile's raw cf max, observed cg min, minus shifts; special-q at the bin's upper edge).

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

## Model of one attempt: first-choice tree, union bound (proof.md section 6)

An attempt takes the first usable boot relation and, at every node, the first usable relation in a fixed cell order; it never retries. It fails iff some node of this first-choice tree has no usable relation (or a cap is hit). No retry, negative-binomial or lazy-cofactorisation law is used; every survivor is cofactored and charged.

Union bound, valid under any dependence between nodes: P(some node has no usable relation) <= E[number of such nodes] = sum over types t of w_t x P(no usable relation at an encountered node of type t), w_t = E[number of type-t nodes in the tree]. w_t is a first moment: the recurrence U(node) = z(node) + E[sum U(children)] with z = 1 on type t (child load below; plan16 keeps w_t, [7w] recomputes the two binding ones): 1,335 nodes per logarithm.

Child load (first moments, no independence). A relation's children are its window primes c_i: U(node) = z + psi, psi := E[sum U(c_i)]; work W(node) = own + sub, sub := E[sum W(c_i)]. E is H1's law given semismoothness (yield integrand weighted by sum g(t_i), g = U or W, max over bins meeting the mesh cell). child_load: upper sums at norm m over lower sums at m + 1 (<= P(m), P nonincreasing in the norm) bound the mean at m = ceil(node maxima; boot ceil(m_max)); no cofactored norm exceeds m (Envelopes); charged x LOAD = 1.001. Smaller norms are assumed not to load more (model).  The plan needs mu >= 18 ([7]).

No usable relation at a node, two cases. (a) Shortfall: the node's certificate fails the rule of PREREG-E (probability s_t for a sampled type-t ideal); charged as a failure. Encountered ideals have density at most c = 1.2 times the sampler's (Encountered law), so E[shortfall nodes] <= c sum_t w_t s_t. No sampled ideal of any type is a shortfall, which has probability prod_t (1 - s_t)^(n_t); at joint 95% confidence this is >= 0.05, so sum_t n_t s_t <= ln 20 and sum_t w_t s_t <= ln(20) max_t w_t/n_t (binding type (69, f): w = 495.68, n = 12,861): c x that = 0.1386 (descent16.p_short). (b) No shortfall: the certificate's mean number of distinct usable relations is >= 2 mu_t/YF (the rule keeps twice the distinct yield before YF); with actual counts >= YF x certificate (toy) the mean is >= 2 mu_t, and for independent cell indicators P(none) <= e^(-2 mu_t). Charged z_t = e^(-mu_t/3) (a sixfold allowance in the exponent); its union over the tree is 0.1042 (the recurrence with z = z_t, [7]).

Encountered law (argument under H1). A degree-one prime ideal I = (q, s, t) divides a cell's norm iff the cell lies in the index-q lattice L_I: probability 1/N(I), the same for every degree-one ideal of that norm. So the children in a bin are ideal-uniform with density (1/b) x (probability that the remaining cofactor is semismooth) in the bit size b. PREREG-E samples ideal-uniform with b uniform on the bin; dom_c bounds the density ratio by its max/mean over a bin, c <= 1.198 at the boot's norm (the largest), charged 1.2 ([7c]).

Toy (toy_descent.py, declared experiment; specification in its header): tower NFS over Z[sqrt 2] (unit 1 + sqrt 2), f = x^3 + x + 1, g = x - m, p = f(m) of 48, 51, 54 bits; 4-dimensional special-q lattices, up to 6 window primes per side, unit duplicates removed (distinct/usable 0.66-0.67); attempts as here: uniform target, first usable relation everywhere, no retry; this file's law (mu = YF x certificate x distinct/usable, z = e^(-mu/3), union over the tree). Recorded, 400 attempts per size (inputs {"size": s, "seed": k, "attempts": 40}, s = 0, 1, 2, k = 1..10, summed): failed attempts 231, 157, 79 against a union charge of 2.23, 1.94, 1.40 per attempt (above 1: vacuous at these sizes, but never below the observed rate); nodes without a relation 394, 216, 91 against 892, 776, 560 charged; distinct relations / mu >= 1.03 in every type with >= 50 nodes. YF = 0.43 is set by this toy: with YF = 0.8 the same ratio would be 0.55 (norms over the quadratic field are less often smooth than rho predicts at small bounds; the gap shrinks from size 0 to 2). Limits: trees of 4.0-4.6 nodes (production about 1,335), base-node mu 0.7-2.3 (production >= 18.1), degree-2 base, no certificate-shortfall case.

Plan: 62 nodes (bin (69, 75) and 5-bit bins 75-225, both sides) and the boot (n1 = 215, nA = 184.58, kA = 15, R = 4); entries (n, iA, kA, kB, R, discard) from a coordinate search min mu 18.1; base nodes R = 224/256 (mu 520/469); all sieve (R = 128..384). Meter, Markov: proof.md sections 6-7.

## Boot

Balls are i.i.d. (fresh exponent and lattice, kept cells tried in full); a ball holds a distinct relation w.p. >= l(1 - l/2), l = cells/per_rel = 2^-30.41 (Bonferroni, H1 pairwise), so cells per relation <= per_rel (1 + 2^-31); the attempt uses the first relation found (stopping time over balls; Wald): trial work per_rel x trial + one ball (2^89.1).

## Envelopes (enforced)

Nodes cofactor a primitive cell only if its certified bounds give log2 N_f <= nf, log2 N_g <= ng (maxima, q removed).
The certificate sums kept cells only (so does PREREG-E's raw); retention, ECM cost, detection, missing ideals and child_load at the maxima hold for every cofactored cell. Survivors per relation are taken as the model's (sv + pr)/pr (assumed monotone in the norms; model). Boot: exact norms; cells above m_max skipped; a ball (exponent,
lattice) is discarded unless its kept cells' H1 yield sum reaches cells x pr, so trials per relation <= 1/pr; discards
charged at BOOT_DISC = 0.9 (< 2^86 instructions per log). Child loads are tabulated (plan16, rounded up); [7t]
recomputes all 242 (>= fresh).

## Supporting measurements

Unit duplicates: a unit eps maps ball cell v to eps v, same relation (tau moves a node's lattice to tau(Q)'s, the boot's to Frob^3(z)'s); O(v) = {eps : eps v in the ball} (exact: unit_orbit.py, Fincke-Pohst, norm +-1), m = |O(v)|/2, and mu needs the distinct yield E[y/m] >= UDUP pr.

PREREG-E (descent_witness.camp_ideal; ideal j of type i in ['boot'] + sorted plan has seed 30000 + 1000 i + j, boot 20000 + j, fixed in advance). enc_ideal: b uniform on the bin, a random prime q of that size, a root s of h mod q kept w.p. (roots)/6, a root t of the side polynomial at s kept w.p. (roots)/12 (g: /2): every degree-one ideal has the same chance; LLL basis as in the descent (boot: fresh exponent, certified lattice). K = 256 uniform ball cells (boot 512; ball_fast: a uniform point of the radius R + sqrt(12)/2 ball rounded, so exactly uniform; [11] chi^2); per cell the certificate integrand y (ss_f ss_g at certified bounds, q removed, 0 above the maxima) and the exact unit multiplicity m. raw = mean y/pr, yr = mean (y/m)/pr (pr charged), raw_lo, yr_lo without the largest cell; shortfall: yr_lo < 2 UDUP or raw_lo < 2 (a rule on this K-cell estimate; the factor 2 is its noise allowance). Result: 48,586 ideals, n_t from 256 to 12,861 (first 1,000 for bins 69-90, 256 elsewhere, then in proportion to w_t; ideals j >= 1000 use seed 10^6 i + j), 0 shortfalls, least yr_lo 1.46 (i = 62). Boot (64 x 512): yr 0.944, mean - 2 SE 0.935 >= 0.65; raw 1.139. plan16 keeps n_t, shortfalls, worst ideal, w_t; audit {"ideal": i} replays it.

Not shown: any relation count or attempt at production size (about 2^-84 per cell). H5 rests on this union bound (p_model <= 0.3058); its inputs at production size are H1 first moments, the PREREG-E certificates, YF and the zero-relation law of (b); the toy measures the last two at toy size.

Refined integration: mesh 384, rho step 1/1024 (lower sums); audit_refinement.py recomputes four probes (0 < lower <= upper); iA = 8 x physical windows.

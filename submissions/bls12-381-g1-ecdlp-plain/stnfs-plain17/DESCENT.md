# Descent

descent16.py is the model (plan, counts, sub-bins: plan16.b85); [7] evaluates it; [7e], [7m], [11], [11b], [12] check inputs. (a, b) in O_h^2 is a vector in Z^12 (basis e_j).

## Lattices

descent_witness.reduce LLL-reduces (delta 0.99) under the integer block form Q(a, b) = round(2^24 sum_s |a(s)|^2) + round(2^24 lam sum_s |w(s)|^e |b(s)|^2), s over the six complex embeddings of K_h. Nodes use QN (lam 3/2, e 1/2), the boot QB (lam 2^64, e 1). Regions: radius-R balls in reduced-basis coordinates, one cell per +- pair.

## Targets and boot

T' = T^e (proof.md section 1): log T = log T'/e mod r (2/r charged). The boot lattice {(a, b) in O_h^2 : a B + b A = 0 in F_p12}, det p^6, has a certified basis (rows in the lattice, Gram determinant p^12, size and Lovasz under QB in exact rationals; [11]); the boot element's single (g-side) norm is N(a^2 + U w b^2). The boot norms each radius-R ball cell (Envelopes) and runs ecm16.py's schedule: a relation has every prime <= 2^n1, <= kA in (2^69, 2^nA], <= 1 in (2^nA, 2^n1]. An exhausted ball draws a new exponent and certified lattice (2^36 each). UDUP = 0.65 as at nodes.

## Nodes

A node is a special-q ideal (l bits, side f/g; (q, y - s) and a root t of the side polynomial at s) with reduced lattice L_q (det q) and a radius-R ball ([7] recomputes the counts R <= 8). Two node types:

- sieve: dsieve16 enumerates every p in [2^8, 2^69] and relative-degree-one root pair on the ball into owners' score arrays. Each primitive cell gets certified upper bounds on both norms (outward-rounded embeddings, charged 2,400); cells above the node maxima are skipped (Envelopes), others with residual mass (log2 bound - score - log2 q on the special-q side) <= T = kA nA + kB n + 96 on both sides survive (exact norms, trial division, ECM); non-primitive cells repeat c/gcd(c)'s relation. Certificate: the bound pass sums ss_f ss_g at kept cells' bounds (q removed); a node below cells x charged pr is discarded before cofactorisation, so every used node's yield is >= charged.
- ecm: ECM on every ball point (never cheaper here).

A relation has every residual prime <= 2^n per side, <= kA in (2^69, 2^nA] and kB <= 1 in (2^nA, 2^n]; these window primes are the children (child load, below).

## Yield and survivors

Per cell, the two-window semismooth probability
P = sum_{jA <= kA, jB <= kB} 1/(jA! jB!) int rho((m - sum t + alpha)/69) prod dt/t
is a lower sum (exact int dt/t per mesh cell, rho at its left endpoint, cells crossing sum t = m dropped); alpha = (0.29, 0.75). Norm sizes m come from offset bins over 61,440 exact norms of uniform ball points on 160 special-q lattices (10 sizes 72-300 bits, both sides, radii 8-128, 64 points each): cf = log2 N_f - 72c, cg = log2 N_g - 12c - 6 log2 U, c = (log2 q + log2 cells)/12. cf sextile maxima pair with the overall cg maximum, shifted up by the worst per-(ideal, radius) mean minus the pooled mean (+14 cf, +5 cg). Boot norms: 3,072 exact norms on 8 certified boot lattices (p90, p99, max + 2, the worst per-target mean minus the pooled mean). For pr the sextiles are subdivided (QBINS_SUB: 24 cf-quantile cells with cg maxima per cell, same shifts; BOOT_SUB: boot deciles + shift); each term stays a lower sum. Charged pr = min(sub-cell mean, 4 x sextile pr) at nodes, min(sub-cell mean, 1.5 x sextile pr) at the boot. mu = cells x P x UDUP x detection x retention (detection 1 - 0.001 x expected factors; the boot carries the counted fraction). Survivors per cell use (T/69) rho((m - T)/69) per side, decreasing in m, at each bin's lowest norms (previous sextile's raw cf max, observed cg min, minus shifts; special-q at the bin's upper edge).

## Factor-base ideals without a logarithm

At the descent input of 2^62.583 relations (G + H's 2-SE count) an ideal of norm q = 2^(69 t) on side s occurs in a relation w.p. kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/rho(u_s), u_s from mean norms 740.0/498.4; a column occurs Poisson(lambda), lambda = 2^62.583 x 4 kappa_s/q, and a column occurring <= 2x is missing. A relation of an m-bit side carries sum_b kappa(m/69, b/69) ln((b+0.1)/b) P(Poisson <= 2) missing ideals on average (1/10-bit bands b < 69). A node discards (relations with a missing ideal rejected, mu x exp(-missing); always at q <= 2^85) or descends (each missing Q sieved in the lattice of Q q', q' a fresh ~2^20 known-log ideal, charged as the worst discarding node with q <= 2^90); the boot descends.

## Retention and cofactorisation

Counted cells: base-primitive only (no prime P of O_h contains both a and b); others repeat (a - b theta)/pi's relation (h_K = 1). Base-primitive fraction prod (1 - N(P)^-2) = 0.6883 (exact for p < 2^9; N(P2) = 2). Cells divisible by a relative-degree-one ideal above a base prime of degree >= 2 with p >= 2^8 are excluded by ordinary fractions 0.9514/0.9917 (f/g, disjoint). Factor base: all prime ideals of norm <= 2^69 above p < 2^8, degree-one ideals above p in [2^8, 2^69].
Omitted mass. A counted cell is rejected only if one side's omitted factor-base mass X exceeds s = 96 - 1 (bound gap < 1e-9) - m/128 (score rounding); X sums independent local terms (whole valuation at p < 2^8, excess powers above), each with one local bound. Odd p: primes (p, h_i), residue field F_q (O_h = Z[y] locally). For side polynomial S and a base-primitive pair, b is a unit w.p. q/(q + 1), else S_hom(a, b) is; then v_P(S_hom(a, b)) = v_P(S(t)), t = a/b uniform in O_P, 0 unless t reduces to a root rho of S in F_q. A root of multiplicity mu is the reduction of exactly mu roots of S over the completion (Hensel), so v_P(S(t)) >= j forces some v_P(t - rho_i) >= j/mu, an event in one class modulo P^ceil(j/mu):
P(v_P >= j and t reduces to rho) <= min(1/q, mu q^-ceil(j/mu)),
with equality q^-j when mu = 1. At 2 (exact divisions by w and pi4): v_P2(F) = 6 when a is in P2 and b a unit (Newton polygon of the f coefficients' valuations 0, 3, 5, 5, 6; probability 1/3), else 0; F has no root in F_4, so v_P4(F) = 0; g is x^2 + U w with a double root 0 at both primes (the bound above with mu = 2). Primes above p are independent (CRT). The local structure at p < 2^9 is tabulated (retain16.structure); [11b] recomputes part of it and checks the bound on residue-class samples at 2 and 19. Above 2^9 excess powers are bounded per dyadic block by 72 p^-2 events per level over pi(x) <= 1.26 x/ln x primes: k events of total multiplicity t have P <= lam^k/k! C(t-1, k-1) 2^(-j(t-k)) (every k, t <= 64; the rest, < sum t (2 max(lam, 2^-j))^t/2, at infinite mass). loss reads a safe summary of X_up (stats16 omitted_law: P(X <= lo) rounded down; masses at lo+1..hi and per 64-bit block, >= 40 merged, rounded up; evidence16 [11b] recomputes it). X_up rounds masses up and puts negligible local tails at infinite mass (counted lost): means 8.49 (f) and 8.45 (g) bits, P(X > 96) = 6.1e-11 (f) and 3.3e-10 (g) per cell.
Relation-conditioned loss. Under H1 a cell with omitted mass x is semismooth with probability Pss(m - x), nonincreasing in m. X_up dominates the local law X (each independent local term only moves up). Among relation cells the lost fraction is L/(L + K) <= L/K, L = E[1(X > s) Pss(m - X)], K = E[1(X <= s) Pss(m - X)]. (i) coupling X <= X_up: 1(X > s) Pss(m - X) <= 1(X_up > s) Pss_hi(m - c64(X_up)) (c64: up to the 64-bit grid), so L <= its X_up mean. (ii) Pss(m - x) >= Pss(m) >= Pss_lo(m) for x >= 0: K >= Pss_lo(m) P(X <= s) >= Pss_lo(m) P(X_up <= s). Hence
loss <= sum_{x > s} P_up(x) Pss_hi(m - c64(x)) / (Pss_lo(m) sum_{x <= s} P_up(x))
at the node's largest norms m_f, m_g (retain16.loss). Retention = (1 - loss_f - loss_g) x base-primitive x (1 - ordinary exclusions); [7] evaluates every sieve node (min 0.6490; largest one-side loss < 5e-5).
ECM (ecm16.py). Suyama Montgomery curves, cheapest (B1, B2/B1) point with wheel 30030; P1 = the Dickman lower sum on the curve order; curves = ceil(ln 0.001 / ln(1 - P1)) at 367.76 (L/5)^2 per mult. [12] prints the schedules and checks P1 on random primes.

## Recurrence and success

Distinct relations per node are negative binomial (mean mu, shape r): r = 0.5 for effective windows of <= 1 prime (base nodes (bin (69, 75)) have empty windows), r = 1 otherwise. Toy fixed balls ([11b]): 0 ideals without a relation vs 2.1-13.2 predicted (seven 1-3-prime windows, means 5.98-22.5); 0-11 vs 6.9-24.0 (12 empty windows).
A subtree fails with probability phi(node) = max((1 + mu (1 - psi)/r)^-r + DISC + psi^32, 0.02); node work = own sieve/ECM + sub/(1 - psi).
Child load (union form). A relation's children are its window primes c_i. Union bound (no independence within a relation): P(fail) <= psi := E[sum phi(c_i)], work <= sub := E[sum W(c_i)]; tries stop at the first success (H5): Wald, work <= sub/(1 - psi). E is H1's law given semismoothness (yield integrand weighted by sum g(t_i), g = phi or W, max over bins meeting the mesh cell). child_load: upper sums at norm m over lower sums at m + 1 (<= P(m), P nonincreasing in the norm) bound the mean at m = ceil(node maxima; boot ceil(m_max)); no cofactored norm exceeds m (Envelopes); charged x LOAD = 1.001. Smaller norms: 8-bit grid from the lowest bin's least norm (boot: lowest decile - 20), 121 sides, max <= 1.0000653 x the value at m; [7m] reruns at 32 bits; between grid points assumed. counts16 keeps the slot-count form. The plan needs mu >= 18, phi <= 0.05, psi <= 0.8 ([7]). DISC = 0.001 per node is a declared H5 constant: 3,100 ideals at the bins' top edges (50 per node type, 96 cells) certify at >= 2.04 x charged pr, none below 1 (pooled 95% bound 0.00097; per type 0 of 50, 0.059).
Conditioning (measured; [7e]). Relations depend only on the ideal and fresh ECM coins; with the certificate, selection can change only the discard rate and the failure law. (A) Production: 600 children ([2^20, 2^24) primes dividing 12 parents' ball cells, half from the 5% smoothest) vs 600 sampled ideals: offsets differ by cf -0.12 +- 0.18, cg +0.04 +- 0.06 bits (KS p 0.42/0.48), maxima 81.1/9.79 (charged 97/14.5); 600 more (256 points) applied to all 62 nodes: certified/charged >= 2.32 (sampled 2.14). (B) Toy tower, one-level descents vs sampled ideals. One-prime windows (48 parents, 1,683 children, r = 0.5): 10.13 +- 0.17 vs 9.92 +- 0.17 relations (KS p 0.75), 24 vs 31 without one (NB 412); P(fail) 0.0324, 0.0297 after a sibling failed. Two-prime (48, 1,486, r = 1): 26.69 vs 26.19 +- 0.36, 1 vs 0 without (NB 65.8). (C) Production-size planted step: 16 parents (120-220 bits) x 24 children at 75-95-bit bin edges, weighted by 32 parent-ball cells' semismooth probability: all 384 certify (ratio >= 4.1); log2 ratio 3.57 unweighted, 3.65 weighted (95% CI of the difference [-0.03, 0.18]).
Plan: 62 nodes (bin (69, 75) and 5-bit bins 75-225, both sides) and the boot (n1 = 215, nA = 184.58, kA = 15, R = 4); entries (n, iA, kA, kB, R, discard) from a coordinate search under the three conditions: min mu 42.1, max phi 0.0302, max psi 0.503 (boot; nodes 0.330); base nodes R = 224/256 (mu 650/586); all sieve (R = 128..384).

## Lazy cofactorisation and boot

Relations are tried singly; the first whose subtree succeeds ends the node. Nodes: after the sieve (charged in full) the s survivors are cofactored in random order, each distinct relation descended at once, until one's children all succeed (failure when survivors run out). With M distinct relations among s survivors the k-th first appears at expected position <= k (s + 1)/(M + 1); the number J tried is independent of the order (H5), P(J > j) <= psi^j, E[min(J, M)] <= 1/(1 - psi):
E[processed]/s <= E_M min(1, (1 + 1/s)/((1 - psi)(M + 1)) + psi^M),
M negative binomial; descent16.lazy_fraction sums it; own work = LLL + sieve + survivors x cost x that fraction ([7l]). Boot: balls are i.i.d. (fresh exponent and lattice, kept cells tried in full); a ball holds a distinct relation w.p. >= l(1 - l/2), l = cells/per_rel = 2^-30.09 (Bonferroni, H1 pairwise), so cells per relation <= per_rel (1 + 2^-31); <= 18 relations, one at a time (stopping time over balls; Wald): trial work per_rel (1 - psi^18)/(1 - psi) x trial + one ball (2^89.1); failure psi^18.

## Envelopes (enforced, H5e)

Nodes cofactor a primitive cell only if its certified bounds give log2 N_f <= nf, log2 N_g <= ng (maxima, q removed).
The certificate sums kept cells; removed cells yield <= pr each, so a node certifying >= 2.04 x charged on all cells
still passes; retention, ECM cost, detection, missing ideals and child_load at the maxima hold for every cofactored cell. Survivors: r(m) = surv/ss is nondecreasing in each norm ([7s]: 120 sides, grid from T to the
maximum; below T, r = 1/ss) and r_f(nf) r_g(ng) <= 2^-11.4 x the model's (sv + pr)/pr at
every node, so under the cutoff survivors per relation stay within the model; the lower bin edges drop out (H1). At
most JMAX = 32 relations per node (psi^32 in phi). Boot: exact norms; cells above m_max skipped; a ball (exponent,
lattice) is discarded unless its kept cells' H1 yield sum reaches cells x pr, so trials per relation <= 1/pr; discards
charged at BOOT_DISC = 0.9 (< 2^84 instructions per log). Child loads are tabulated (plan16, rounded up); [7t]
recomputes all 242 (>= fresh), audit {"load": k} one.

## Supporting measurements

[11] re-samples a boot and a 260-bit f special-q lattice against the bins. Unit duplicates: a unit eps maps ball cell v to eps v, same relation (tau moves a node's lattice to tau(Q)'s, the boot's to Frob^3(z)'s); O(v) = {eps : eps v in the ball}, m_true = |O(v)|/2, and P x UDUP needs yr = E[s/m_true]/pr >= UDUP. O(v) is exact (unit_orbit.py, Fincke-Pohst, norm +-1); [11] measures it on the shipped plan (24 boot cells, 12 at (220, g); thresholds 0.70/1.0): boot yr 0.784, 220g 2.63.
Refined integration: mesh 384, rho step 1/1024 (lower sums); audit_refinement.py recomputes four probes (0 < lower <= upper); iA = 8 x physical windows.

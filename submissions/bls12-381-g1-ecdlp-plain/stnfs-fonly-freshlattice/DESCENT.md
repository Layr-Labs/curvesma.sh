# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T becomes T' = T^e, e uniform in [1, r), log T = log T'/e mod r (2/r charged). The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6; its LLL basis is certified (Gram determinant p^12,
Lovasz conditions; [11]); the boot element has a single (g-side) norm. The boot enumerates the radius-4 ball
(2^23.67 points, one of each +- pair), computes each norm and runs the
`ecm16.py` schedule for 193.71-bit primes: a relation has every prime <= 2^215,
<= 16 in (2^69, 2^193.71], <= 1 in (2^193.71, 2^215]. Every boot trial draws a
new exponent and certified lattice (2^36 each), so trials are independent. Distinct relations credit
UDUP = 0.65, as at nodes. Lazy boot (below): at most 18 distinct relations, each descended before the next
is sought; expected 2^56.68 trials of 2^67.02 instructions.

## Nodes

A node is a special-q ideal (l bits, side f/g) with LLL-reduced lattice L_q (det q); its region is the radius-R ball
in that basis, one point per +- pair (half counts in `stats16.b85`; [7] recomputes R <= 8). Two node types:

- **sieve**: every p in [2^8, 2^69] and every relative-degree-one root pair
is enumerated on the ball by `dsieve16` (KERNELS.md); records go to the
owners' 2^40-entry score arrays. Every primitive cell gets a certified upper
bound on both norms (72 + 12 linear forms in double precision, outward
rounding; ~1,500 instructions, charged 2,400). Cells with residual mass
(log2 bound - score - log2 q on the special-q side) <= T = kA nA +
kB n + 96 on both sides survive (exact norms, trial division, ECM).
Non-primitive cells repeat c / gcd(c)'s relation and are skipped.
- **ecm**: ECM trials on every ball point, no sieve (never cheaper here).

A relation has every residual prime <= 2^n per side, at most kA in (2^69, 2^nA]
and at most kB <= 1 in (2^nA, 2^n]; these window primes are the children, each
charged at the most expensive bin at or below its window.

## Yield and survivors

Per cell, the two-window semismooth probability

P = sum_{jA <= kA, jB <= kB} 1/(jA! jB!) int rho((m - sum t + alpha)/69) prod dt/t

is a lower sum (exact weight int dt/t per mesh cell, rho at its left endpoint; cells crossing sum t = m dropped). Norm
sizes m come from measured offset bins over 61,440 exact ball norms on 160 production special-q lattices (10 sizes
72-300 bits, both sides, 6 radii): cf sextile maxima paired with the overall cg maximum, shifted by the worst sampled
per-(ideal, radius) mean minus the pooled mean (cf 20.68 + 1.58 -> +23, cg -2.53 + 7.81 -> +6); boot norms from 3,072
exact norms on 8 certified boot lattices (p90 + 4, p99 + 4 bins). **Sub-bin yields.** For pr the sextiles are
subdivided (`QBINS_SUB`: 24 cf-quantile cells, per-cell cg maxima on a fresh grid, same +23/+6, top edge 109.6/20.5;
`BOOT_SUB`: 10 R = 4 deciles + 4, edge 16.7); sub-cell maxima dominate their cells, so each term stays a lower sum.
Charged pr = min(sub-cell mean, 4 x sextile pr) at nodes, min(sub-cell mean, 1.5 x sextile pr) at the boot (UDUP is
measured against it); survivors include cells x pr. mu = cells x P x UDUP x detection x retention (UDUP =
0.65; ECM detection 1 - 0.001 x expected factors; the boot carries only the counted fraction). Survivors per cell use
(T/69) rho((m - T)/69) per side, decreasing in m, at each bin's lowest norms (previous sextile's raw cf max, observed cg
min -54.56, margins (23, 6); special-q at the bin's upper edge).

## Factor-base ideals without a logarithm

Under the main-phase occurrence model (`filtsim.py` at the descent input of 2^63.575 relations, [4]) an ideal of norm
q = 2^(69 t) on side s occurs with probability 2 kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/rho(u_s), u_s from mean norms
738.1/498.2 (G's means truncated; [4], [10y]); occurrences are Poisson(lambda_s(b)) and an ideal occurring
<= 2x is missing. A relation of an m-bit side carries sum_b kappa(m/69, b/69) ln((b+0.1)/b) P(Poisson <= 2) missing
ideals on average (1/10-bit bands b < 69): 0.20-1.86 per f side, 0.41-3.72 per g side over u = 7..39 (3.773 at the
boot's 2,709-bit norm). Each node either **discards** (relations with a missing ideal rejected, mu x exp(-missing);
always at q <= 2^90 (forced up to 85 bits; bin 90 keeps its searched thin flag)) or **descends** (each missing ideal Q
is a child sieved in the lattice of Q q', charged as the worst discarding node with q <= 2^90, Q q' <= 2^89); the boot
descends. Cascading removal can lose ideals with >= 3 occurrences: [6b] runs purge (target 0.0001) and
back-substitution on the [6] system (excess 1.08, 2^12, seed 1) and counts ideals left without a log per lambda band;
over seeds 1-4 at 2^12 and 1-2 at 2^13 the model (every ideal with <= 2 occurrences) is >= 6.3x the measured count in
every band with a loss.

## Retention and cofactorisation

**Counted cells.** A cell counts only if base-primitive (for every p and irreducible h of H mod p, a and b not both
divisible by (p, h)); other cells repeat (a - b theta)/pi's relation up to p's free relation. Base-primitive
fraction prod (1 - p^-2 deg h) = 0.9100 (exact degrees for p < 2^9). Cells
divisible by a relative-degree-one ideal above a base prime of degree >= 2
with p >= 2^8 are excluded by ordinary fractions 0.9523/0.9919 (f/g, disjoint).
The factor base holds every prime ideal of norm <= 2^69 above p < 2^8 plus the
degree-one ideals above p in [2^8, 2^69].

**Omitted mass.** A counted cell is rejected only when one side's omitted
factor-base mass X exceeds s = 96 - 1 (bound gap < 2e-11 bits) - m/128
(rounding: each of the <= m/8 sieved hits of an m-bit side loses < 1/16 bit to
floor(16 log2 p)). X sums independent local contributions (whole valuation at
p < 2^8, excess powers at p >= 2^8), each with one local bound.
For P = (p, h_i) a prime of the base field (Kummer-Dedekind, O = Z[y]) with residue field F_q, q = p^deg h_i, S
the side polynomial (monic in x) and a base-primitive pair: b is a unit with probability q/(q + 1), otherwise a is and
S_hom(a, b) is a unit. With b a unit, v_P(S_hom(a, b)) = v_P(S(t)) for t = a/b uniform in O_P; it is 0 unless t
reduces to a root rho of S in F_q. A root of multiplicity mu is the reduction of exactly mu roots of S over the
completion (Hensel), so v_P(S(t)) >= j forces some v_P(t - rho_i) >= j/mu, an event in one class modulo P^ceil(j/mu):

P(v_P >= j and t reduces to rho) <= min(1/q, mu q^-ceil(j/mu)),

with equality q^-j when mu = 1. The primes above p are independent (CRT), v_p(N) = sum_P deg(h_i) v_P. Degrees,
ramification and root multiplicities of f and g at every P with p < 2^9 are tabulated (`retain16.structure`,
`stats16.b85`); [11b] recomputes the table at 2, 19, 23 (repeated roots or ramification) and a seventh of the others,
and checks the bound against residue-class samples at 2 and 23. All multiplicities are 1 except g at the two primes
above 2 (2^16 | U: the double root x = 0, q = 4 and 16) and f at the degree-4 prime above 19 (two double roots among
four, q = 19^4); 23 is ramified in the base field (e = 2), g having two simple roots at each of its two degree-one
primes. Above 2^9 excess powers are bounded per dyadic block by 72 p^-2 events per level over pi(x) <= 1.26 x/ln x
primes: k events of total multiplicity t have P <= lam^k/k! C(t-1, k-1) 2^(-j(t-k)) (every k, t <= 64; the rest,
< sum t (2 max(lam, 2^-j))^t/2, at infinite mass). This law X_up rounds masses up to 1/4 bit, moves probabilities below
1e-30 up to the 64-bit grid and puts local tails (below 1e-200 or past level 1,000) at infinite mass (counted lost): means 3.35 (f) and
7.67 (g) bits, P(X > 96) = 1.6e-12 (f) and 2.6e-10 (g) per cell.

**Relation-conditioned loss.** Under H1 a cell with omitted mass x is semismooth with probability Pss(m - x),
nonincreasing in m. X_up dominates the local law X (each independent local term only moves up). Among relation cells
the lost fraction is L/(L + K) <= L/K, L = E[1(X > s) Pss(m - X)], K = E[1(X <= s) Pss(m - X)]. (i) coupling X <= X_up:
1(X > s) Pss(m - X) <= 1(X_up > s) Pss_hi(m - c64(X_up)) (c64: up to the 64-bit grid), so L <= its X_up mean. (ii)
Pss(m - x) >= Pss(m) >= Pss_lo(m) for x >= 0: K >= Pss_lo(m) P(X <= s) >= Pss_lo(m) P(X_up <= s). Hence

loss <= sum_{x > s} P_up(x) Pss_hi(m - c64(x)) / (Pss_lo(m) sum_{x <= s} P_up(x))

at the node's largest norms m_f, m_g (`retain16.loss`). Retention = (1 - loss_f - loss_g) x base-primitive x (1 -
ordinary exclusions); [7] evaluates every sieve node (min 0.8592; largest one-side loss < 1e-5).

**ECM (`ecm16.py`).** Suyama Montgomery curves, cheapest (B1, B2/B1) point with wheel 30030; x-only stage 1 (11
mults/bit), stage-2 cross products (3/pair), one binary gcd per stage; P1 = the Dickman lower sum on the curve order;
curves = ceil(ln 0.001 / ln(1 - P1)) at 367.76 (L/5)^2 per mult. [12] prints node/boot schedules (boot: B1 = 2^28,
737933 curves, 2^67.01/trial) and matches P1 on random 36/40-bit primes.

## Recurrence and success

Distinct relations per node are negative binomial (mean mu, shape r): r = 0.5 for effective windows of <= 1 prime
(base nodes (bin (69, 75)) have empty windows), r = 1 otherwise. Toy-ball fits ([11b]): one-
prime windows 0.60, two/three-prime 1.86-12.2. Of 12 empty-window configs, those with mean >= 7.4 show 0-4 ideals
without a relation vs 5.4-12.4 predicted at r = 0.5; the three with means 3.2-4.0 show 19/11/21 vs 18.1/17.2/16.4, and
raising lb clears them (0 vs 7.3-11.2).

A child subtree fails with probability phi; phi(node) = max((1 + mu (1 -
psi)/r)^-r, 0.02), psi = window sum of children's phi; a node tries <=
1/(1 - psi) relations on average. Node work = own sieve/ECM + 1/(1 - psi)
x children's work. The plan requires mu >= 18, phi <= 0.05, psi <= 0.8
everywhere ([7] checks). Node counter (N_MAX = 2^27.3, `counts16.py`) and repeats: proof.md section 7.
A factor-base ideal Q without a log (norm <= 2^69) is never descended in a fixed lattice: each visit draws a fresh
q' of about 2^20 with known log and sieves the lattice of Q q' (charged as above).

The shipped plan: 62 nodes (bins (69, 75) plus 5-bit bins 75-225, both sides) and the boot (n1 = 215, nA = 193.71,
kA = 16, R = 4), a coordinate search over (n, iA, kA, kB, R, discard) minimising the metered cap under the three conditions. At this package's inputs (LOG2_RELS = 63.575, means 738.1/498.2, UDUP =
0.65, THIN_MAX 85): min mu 44.6, max phi 0.0268, max psi 0.537 (boot; nodes 0.391). Base nodes sieve R = 256 (f, mu
932, phi 0.0232) and R = 288 (g, mu 786, phi 0.0253); all 62 nodes sieve (R = 192..384). One log <= 2^125.719
instructions (boot own 2^123.70 + subtree 2^124.20 x retry 2.159); **two logs <= 2^112.4316 units**. Meter, Markov
bound and boot failure: proof.md sections 6-7 ([7]).

## Lazy cofactorisation and lazy boot

Relations are tried singly; the first whose subtree succeeds ends the node (all else as eager). **Nodes.** After the sieve (charged in full) the s survivors are cofactored in random order, each
distinct relation descended at once, until one's children all succeed (failure when survivors run out, as eager). With M distinct relations among s survivors the k-th first
appears at expected position <= k (s + 1)/(M + 1); the number J tried is independent of the order (H5) with
P(J > j) <= psi^j, so E[min(J, M)] <= 1/(1 - psi) and

E[processed]/s <= E_M min(1, (1 + 1/s)/((1 - psi)(M + 1)) + psi^M),

M negative binomial as in phi_nb. `descent16.lazy_fraction` sums this (stopping once m > mu and a term < 1e-20, the unsummed mass and pm (r + mu)/r at
fraction 1); own work is LLL + sieve + survivors x per-survivor cost x that fraction; [7] simulates it (`lazy_mc`). **Boot.** One distinct relation is found (expected per_rel trials) and descended before the next, at most
BOOT_RELS = 18. Each boot trial draws a fresh special-q exponent and certified lattice, so outcomes are i.i.d. and
min(J, 18) is a stopping time: by Wald the trial work is per_rel x (1 - psi^18)/(1 - psi) x trial, LLL charged per
trial. Failure psi^18.
## Supporting measurements ([11])

Norm bins (offline, 20 lattices/size): lattices at 150/260 bits stay in the node bins (cf max 53.4 <= 109.6, cg max
5.19 <= 20.5), the boot lattice in the boot bins (cg max 6.6 <= 16.7); [11] re-samples per run. **Unit duplicates.** A
unit eps of O_h maps a ball cell v to eps v (same relation): v's duplicates are O(v) = {eps : eps v in the ball}, the
distinct yield per cell E[s/m_true], m_true = |O(v)|/2, and P x UDUP needs yr = E[s/m_true]/pr >= UDUP. O(v) is the integer point set of e M_v M_v^T e^T <= R^2 (M_v the B-coordinates of y^j v), enumerated exactly by
`unit_orbit.py` (rational LDL^T, Fincke-Pohst, no box; norm +-1). Campaign (STNFS-KIT, not shipped; 8 certified boot lattices x 250 R = 4 cells, one special-q lattice per plan entry
at the bin's top edge x 100 cells; 8,200 cells) against the charged pr: boot yr min 0.890 (95% bootstrap
lower end 0.842), mean 0.933; nodes min 2.59 (220g; smallest lower end 2.09, 195g), all above 0.65. Orbits reach 44
cells (heights 34) at the boot, 62 (75) at nodes (no box is safe). [11] repeats it live (24 boot
cells, 12 at (220, g); thresholds 0.70 / 1.0; this build 0.896 / 2.05) with a sampler check.

## Refined integration representation

Mesh 384, rho spacing 1/1024: eight sub-cells per old cell, exact log measure, rho at left endpoints (lower
sums); `audit_refinement.py` recomputes four probes. Executable iA tuples are the physical windows times eight.

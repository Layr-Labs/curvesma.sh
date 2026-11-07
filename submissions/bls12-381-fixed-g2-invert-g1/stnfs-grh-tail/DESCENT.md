# Revision convention

Use proof.md's audited CDF/guards, capped prefix, H5n/H5p and GRH counts. Final numerical bounds are in refinement.json. Historical norm bins are proposals, not universal envelopes; the0.02floor is not a norm-tail tax. The old+10%claim is withdrawn.

# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T becomes T' = T^e, e uniform in [1, r), log T = log T'/e mod r (2/r charged). The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6; its LLL basis is certified (Gram determinant p^12,
Lovasz conditions; [11]); the boot element has a single (g-side) norm. The boot enumerates the radius-4 ball
(2^23.67 points, one of each +- pair), computes each norm and runs the
`ecm16.py` schedule for 193.71-bit primes: a relation has every prime <= 2^215,
<= 16 in (2^69, 2^193.71], <= 1 in (2^193.71, 2^215]. An exhausted ball draws a
new exponent and certified lattice (2^36 each). Distinct relations credit
UDUP = 0.65, as at nodes. Lazy boot (below): at most 18 distinct relations,
each descended before the next is sought; expected 2^56.68 trials of 2^67.02 instructions.

## Nodes

For a special-q ideal of lbits, certify the LLL-reduced12-dimensional lattice of determinant q and enumerate its radius-R half-ball (CELLS gives exact counts). The shipped62node plan uses sieving: all relative-degree1root pairs over factor-base primes2^8..2^69 are walked by dsieve16, routed to2^40-entry owner score arrays and added. Primitive cells have outward norm bounds (43places after conjugation), guard/corridor checks from proof.md, and survive when each side's bound minus score minus the special-q log is<=T=kA*nA+kB*n+96. Nonprimitive cells are skipped. The chosen ordered survivor prefix gets exact norms, trial division and ECM. A relation has all residual primes<=2^n, at most kA in(2^69,2^nA] and kB<=1 in(2^nA,2^n]; children are charged at their largest eligible bin. The unused ECM-only node alternative tests every point and is never cheaper in this plan.

## Yield and survivors

The two-window semismooth lower sum is P=sum_{jA<=kA,jB<=kB}1/(jA!jB!) integral rho((m-sum(t)+alpha)/69)prod(dt/t). Exact log-measure weights, endpoint rho lower values and discarded crossing cells define descent16.ss_lower2. Historical norm sextiles/24sub-bins and boot deciles are proposals; proof.md's runtime audited CDF proves their eta=.95yield use on accepted bases. Counts, bad-basis tax and marked prefix work are specified there. QBINS_LO only sets the survivor budget; a pointwise interval ratio certificate replaces any assertion that these empirical lower bins universally bound work. The supplied original measurements are support at their actual scale, not production proofs.

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
factor-base mass X exceeds s = 96 - 1 (certified bound gap < 1 bit) - m/128
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

A child subtree fails with probability phi; phi(node) = delta + (1-delta) max((1 + mu (1 -
psi)/r)^-r, 0.02), delta=1/200, psi = window sum of children's phi; a node tries <=
1/(1 - psi) relations on average. Node work = own sieve/ECM + 1/(1 - psi)
x children's work. The plan requires mu >= 18, phi <= 0.05, psi <= 0.8
everywhere ([7] checks). Node counter (N_MAX = 2^27.3, `counts16.py`) and repeats: proof.md section 7.
A factor-base ideal Q without a log (norm <= 2^69) is never descended in a fixed lattice: each visit draws a fresh
q' of about 2^20 with known log and sieves the lattice of Q q' (charged as above).

The shipped plan: 62 nodes (bins (69, 75) plus 5-bit bins 75-225, both sides) and the boot (n1 = 215, nA = 193.71,
kA = 16, R = 4), a coordinate search over (n, iA, kA, kB, R, discard) minimising the metered cap under the three conditions. At this revision's audited inputs, eta=.95 and delta=.005, all62nodes remain feasible. The complete final values are recomputed in refinement.json/audit_refinement.py; the two-log expected-work upper is113.329, boot phi<=.001755, and node count log2<=16.682.

## Lazy processing

The original uncapped lazy charge is E_NB[min(1,1/((1-psi)(M+1))+psi^M)] with its tail charged at1, computed by lazy_fraction. The revision uses a capped uniformly ordered survivor reservoir and the explicit H5p pointwise-rate proof in proof.md. Its geometric processed-fraction bound is checked against every original charged lazy fraction. Bootstrap still processes one relation and its descendants at a time, up to18relations; Wald gives per_rel*(1-psi^18)/(1-psi)trials. Each exhausted ball uses a newly audited basis, and every failed audit is charged below.

## Supporting measurements and limitations

The original norm campaign has61440exact norms on160production-sized lattices, with production-size offsets and fresh sub-bin measurements. Unit orbits use exact rational LDL/Fincke–Pohst enumeration in unit_orbit.py; the8200cell campaign gives boot distinct-yield ratio minimum.890 (95%bootstrap lower.842), node minimum2.59 (smallest lower2.09), all aboveUDUP=.65. Orbits up to62cells show that no fixed unit box is safe. Live evidence16.py replays selected production-size lattices,24boot/12node orbits and a uniform sampler check. These do not prove H5/H5p. New held-outCDF and bound-gap measurements, exact minima, scope and hashes are in stats16.b85 and proof.md; public notes contain full encoded counts and reproduction sources.

## Refined integration representation

Mesh 384, rho spacing 1/1024: eight sub-cells per old cell, exact log measure, rho at left endpoints (lower
sums); `audit_refinement.py` recomputes four probes. Executable iA tuples are the physical windows times eight.

## Audits and repeated bootstrap balls

proof.md gives the complete exact sampling/threshold and multiprecision bounds. Node audits cost<=2^28(2^78+2^90+2^61)instructions. EACH fresh radius4boot basis is audited by enumerating all13,315,520half-ball points; charge<=2^69instructions/basis. The two-log meter charges each cofactor trial at>=2^67instructions, allowing<2^62trials and<2^40exhausted/partial good balls. Up to128CDFattempts per ball yield<2^47audits and exhaustion probability<=2^-88 under H5n's boot badmass<=1/2. Every failed basis/audit is charged. Combined node/boot, integer-root setup, allocation/rounding and random-generation work fit2^105units; total audit/prefix errors fit2^-80. Samples and threshold counters stream with bounded memory.

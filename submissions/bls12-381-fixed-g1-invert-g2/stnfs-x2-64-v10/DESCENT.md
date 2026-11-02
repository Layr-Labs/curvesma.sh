# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T is replaced by T' = T^e, e uniform in [1, r), log T = log T'/e
mod r: uniform on nonidentity mu_r (distance 1/r from uniform on all mu_r;
identity z handled directly; ledger charges 2/r). The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6; its LLL basis is
certified exactly (Gram determinant p^12, Lovasz conditions; [11]). The boot
element has a single (g-side) norm. The boot enumerates the radius-4 ball
(2^23.67 points, one of each +- pair), computes each norm and runs the
`ecm16.py` schedule for 199.79-bit primes: a relation has every prime <= 2^215,
<= 15 in (2^69, 2^199.79], <= 1 in (2^199.79, 2^215]. An exhausted ball draws a
fresh exponent and certified lattice (2^36 each). Distinct relations credit
UDUP = 0.6, as at nodes. It stops at 18 distinct relations: 2^58.97 trials of
2^67.96 instructions over ~2^35.3 lattices, tried in turn; failure needs all
18 to fail.

## Nodes

A node is a special-q ideal (l bits, side f/g) with lattice L_q (det q,
LLL-reduced). Its region is the radius-R integer ball in the reduced basis,
one point per +- pair; exact half counts are in `stats16.b85` (recomputed for
R <= 8 by [7]). Two node types:

- **sieve**: every p in [2^8, 2^69] and every relative-degree-one root pair
  is enumerated on the ball by `dsieve16` (KERNELS.md); records go to the
  owners' 2^40-entry score arrays. Every primitive cell gets a certified upper
  bound on both norms (72 + 12 linear forms in double precision, outward
  rounding; ~1,500 instructions, charged 2,400). Cells with residual mass
  (log2 bound minus score, minus log2 q on the special-q side) <= T = kA nA +
  kB n + 96 on both sides survive (exact norms, trial division, ECM).
  Non-primitive cells repeat c / gcd(c)'s relation and are skipped.
- **ecm**: ECM trials on every ball point, without a sieve (never cheaper in
  the shipped plan).

A relation has every residual prime <= 2^n per side, at most kA in (2^69, 2^nA]
and at most kB <= 1 in (2^nA, 2^n]; these window primes are the children, each
charged at the most expensive bin at or below its window.

## Yield and survivors

Per cell, the two-window semismooth probability

  P = sum_{jA <= kA, jB <= kB} 1/(jA! jB!) int rho((m - sum t + alpha)/69) prod dt/t

is a lower sum: each mesh cell carries exact weight int dt/t and rho at its
left endpoint (rho decreases); cells crossing sum t = m are dropped. Norm
sizes m come from measured offset bins over 61,440 exact ball norms on 160
production special-q lattices (10 sizes 72-300 bits, both sides, 6 radii): cf
sextile maxima paired with the overall cg maximum, shifted by worst sampled
per-(ideal, radius) mean minus pooled mean (cf 20.68 + 1.58 -> +23, cg -2.53 +
7.81 -> +6). Boot norms use 3,072 exact norms on 8 certified boot lattices
(p90 + 4, p99 + 4 bins). mu = cells x P x UDUP x detection x retention (UDUP
= 0.6; ECM detection 1 - 0.001 x expected factors; retention from the next
section; the boot carries only the counted fraction). Survivors per cell use
(T/69) rho((m - T)/69) per side, decreasing in m, evaluated at each bin's
lowest norms (previous sextile's raw cf max, observed cg min -54.56, margins
(23, 6); special-q at the bin's upper edge).

## Factor-base ideals without a logarithm

Under the main-phase occurrence model (`filtsim.py`: 2^63.949 relations -- the
theta = 2901/8000 retained count floored to 0.001 bit, certified lower product
0.971219877 already applied), an ideal of norm q = 2^(69 t) on side s occurs
with probability 2 kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/rho(u_s), u_s from
mean norms 746.9/499.8 (rounded down from 746.94/499.84; [4] checks them).
Occurrences are Poisson(lambda_s(b)); an ideal occurring <= 2x is missing. A
relation of an m-bit side carries on average sum_b kappa(m/69, b/69)
ln((b+0.1)/b) P(Poisson <= 2) missing ideals (1/10-bit bands b < 69):
0.09-0.83 per f side, 0.26-2.32 per g side over u = 7..39 (2.35 at the boot's
2,709-bit norm). Each node takes one option:

- **discard**: relations containing a missing ideal are rejected, multiplying
  mu by exp(-missing); nodes with q <= 2^90 always discard (the code forces
  it for bins up to 85 bits; bin 90 keeps its searched thin flag), so missing
  ideals never recurse there;
- **descend**: each missing ideal Q is a further child, sieved in the lattice of
  Q q' (see Repeated ideals) and charged the largest work and failure probability
  of the discarding nodes with q <= 2^90 (Q q' <= 2^89; smaller Q are cheaper);
  the boot descends its missing ideals the same way.

Cascading removal can also lose ideals with >= 3 occurrences; [6b] measures
losses on the [6] system: for lambda in [1,2), [2,3), [3,5), [5,8), >= 8 the
missing fraction is 0.0470/0.0309/0.0224/0.0036/0 vs model's
0.7826/0.5430/0.2811/0.0603/0.0017; ideals lost with >= 3 occurrences: 0 of
3,816.

## Retention and cofactorisation

**Counted cells.** A cell is counted only if base-primitive (for every p and
irreducible h of H mod p, a and b not both divisible by (p, h)); other cells
repeat (a - b theta)/pi's relation up to p's free relation. Base-primitive
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
Let P = (p, h_i) be a prime of the base field (Kummer-Dedekind, O = Z[y]) with
residue field F_q, q = p^deg h_i, and S the side polynomial, monic in x. For a
base-primitive pair, (a : b) is a point of P^1(O_P): b is a unit with
probability q/(q + 1), otherwise a is one and S_hom(a, b) is a unit. With b a
unit, v_P(S_hom(a, b)) = v_P(S(t)) for t = a/b uniform in O_P. If t does not
reduce to a root of S in F_q, the valuation is 0. A root rho of multiplicity mu
in F_q is the reduction of exactly mu roots of S in an algebraic closure of the
completion (Hensel), so v_P(S(t)) is the sum of the mu distances v_P(t - rho_i),
and v_P(S(t)) >= j forces some v_P(t - rho_i) >= j/mu, an event of probability
at most q^-ceil(j/mu) (it is contained in one class modulo P^ceil(j/mu)):

  P(v_P >= j and t reduces to rho) <= min(1/q, mu q^-ceil(j/mu)),

with equality q^-j when mu = 1 (one root, unique lift). The primes above p are
independent (CRT) and v_p(N) = sum over P of deg(h_i) v_P. The degree,
ramification index and the root multiplicities of f and g in F_q of every P
with p < 2^9 are tabulated (`retain16.structure`, `stats16.b85`); [11b]
recomputes the table at every prime with a repeated root or ramification (2, 19,
23) and at a seventh of the others, and checks the bound against residue-class
samples at 2 and 23. All multiplicities are 1 except g at the two primes above 2
(2^16 | U: the double root x = 0, q = 4 and 16) and f at the degree-4 prime
above 19 (two double roots among four, q = 19^4); 23 is ramified in the base
field (e = 2) with g having two simple roots at each of the two degree-one
primes. Above 2^9 the excess powers are bounded per dyadic block by 72 p^-2
events per level over pi(x) <= 1.26 x/ln x primes (one- and two-event terms).
Masses are rounded up to 1/4 bit and probabilities below 1e-30 are moved to the
largest kept mass: means 3.35 (f) and 7.66 (g) bits, P(X > 96) = 1.2e-12 (f)
and 2.3e-10 (g) per cell.

**Relation-conditioned loss.** Under H1 a cell with omitted mass x is
semismooth with probability Pss(m - x), so the per-side lost fraction is

  loss = sum_{x > s} P(X = x) Pss_hi(m - x) / sum_{x <= s} P(X = x) Pss_lo(m - x),

bracketed by monotonicity on a 64-bit grid in x at the node's largest norms
m_f, m_g. Retention = (1 - loss_f - loss_g) x base-primitive x (1 - ordinary
exclusions); [7] evaluates every sieve node (min 0.8593; largest one-side loss
< 1e-5).

**ECM schedule (`ecm16.py`).** Suyama-form Montgomery curves, cheapest (B1,
B2/B1) grid point with wheel 30030; stage-1 x-only ladder (11 mults/bit),
stage-2 cross products (3/pair), one binary gcd per stage. Per-curve success
P1 is the Dickman lower sum on the curve order; curves = ceil(ln 0.001 /
ln(1 - P1)) at 367.76 (L/5)^2 per mult. [12] prints node/boot schedules (boot:
B1 = 2^28, 1418167 curves, 2^67.96/trial) and matches P1 on random 36/40-bit
primes.

## Recurrence and success

Distinct relations per node are negative binomial (mean mu, shape r): r = 0.5
for effective windows of <= 1 prime (base nodes (bin (69, 75)) have empty
windows), r = 1 otherwise. Shapes come from toy-ball fits (empty-window rows
in `stats16.b85`, printed by [11b]): one-prime windows fit 0.60, two/three-prime
1.86-12.2. Of 12 empty-window configs, all with mean >= 7.4 show 0-4 ideals
without a relation vs 5.4-12.4 predicted at r = 0.5; the three with means
3.2-4.0 show 19/11/21 vs 18.1/17.2/16.4, and raising lb clears them (0 vs
7.3-11.2). Base nodes run at mu = 571 and 438.

A child subtree fails with probability phi; phi(node) = max((1 + mu (1 -
psi)/r)^-r, 0.02), psi = window sum of children's phi; a node tries <=
1/(1 - psi) relations on average. Node work = own sieve/ECM work + 1/(1 - psi)
x children's work. The plan requires mu >= 18, phi <= 0.05, psi <= 0.8
everywhere ([7] checks). The descent keeps a special-q node counter over both
logs and aborts past N_MAX = 2^27.3 nodes (counted failure); E[N] <= 2^17.74
(`counts16.py`, unrounded 2^17.7317), so P(abort) <= 2^-9.56 (Markov). All N_MAX nodes in the
smallest bin (2^70.3 ideals) with 2^6 for 1/p weighting: shared-ideal
probability <= 2^-10.70, charged in full ([7]).
A factor-base ideal Q without a log
(norm <= 2^69) is never descended in its own fixed lattice: each visit draws a
fresh ideal q' of about 2^20 whose log is known (main-phase factor base) and
sieves the lattice of Q q', so a recurrence of Q samples a different sublattice.
Its node is charged as the worst node with q <= 2^90, and all nodes up to that
size discard relations touching an ideal without a log (no deeper recursion).

The shipped plan has 62 nodes (bins (69, 75) plus 5-bit bins 75-225, both
sides) and the boot (n1 = 215, nA = 199.79, kA = 15, R = 4), chosen by
coordinate search over (n, iA, kA, kB, R, discard) minimising modelled two-log
work subject to the three conditions (history: refinement.json, REFINEMENT.md).
At this revision's theta = 2901/8000 inputs (LOG2_RELS = 63.949, mean norms
746.9/499.8, UDUP = 0.6) with THIN_MAX 85 (35 entries differ from the base
plan, boot unchanged): min mu 33.8, max phi < 0.05, max psi 0.644 (boot). Base
nodes sieve R = 288 (f, mu ~571, phi ~0.030) and R = 320 (g, mu ~438, phi
~0.034); all 62 nodes sieve (R = 192..384). One log <= 2^129.146 instructions;
**two logs <= 2^115.8579 units**. Meter stops the descent at 2^116.59 units
(Markov: E[W]/cap <= 0.6020270). P(boot fails) <= 3.600e-04 ([7]).

## Supporting measurements ([11])

Norm bins (offline, 20 lattices/size): fresh lattices at 150/260 bits stay in
the node bins (cf max 53.4 <= 109.6, cg max 5.19 <= 20.5), the certified boot
lattice in the boot bins (cg max 6.6 <= 16.7); [11] re-samples a few points per
run. Unit duplicates (offline): mates of 480 ball points under 840 small units
give distinct/raw 0.779; the model charges 0.6 ([11] repeats 12 points: 0.944).
Unit-count convergence (offline, f-side 260-bit special-q lattice seed 13300,
R = 96, 200 points, `descent_witness.small_units/unit_mates` seed 13307):
distinct/raw 0.8514/0.8394/0.8353 (SE 0.023) at 500/840/1190 units,
flattening as large units are added; 0.6 is 10 SE below the 1190-unit value.

## Refined integration representation

Mesh 384 (was 48), rho spacing 1/1024 (was 1/128). The historical iA values above are physical windows; the executable tuples scale them by eight. Kernels, prime bounds and windows are unchanged (REFINEMENT.md).

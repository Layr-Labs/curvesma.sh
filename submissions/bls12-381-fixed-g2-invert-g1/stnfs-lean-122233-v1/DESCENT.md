# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T is replaced by T' = T^e, e uniform in [1, r), log T = log T'/e
mod r. For nonidentity T this is uniform on the nonidentity elements of mu_r, whose statistical distance from uniform on all mu_r is 1/r. Identity z is handled directly; g is nonidentity. The success ledger charges 2/r for the two logarithms. The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6; its LLL basis is
certified exactly (Gram determinant p^12, Lovasz conditions; [11]). The boot
element has a single (g-side) norm. The boot enumerates the radius-4 ball of the
lattice (2^23.67 points, one of each +- pair), computes each norm and runs the
`ecm16.py` schedule for 199.79-bit primes: a relation has every prime <= 2^215,
<= 15 primes in (2^69, 2^199.79] and <= 1 in (2^199.79, 2^215].
When a ball is exhausted it draws a fresh exponent e and a new certified
lattice (2^36 instructions each). Distinct relations are credited at UDUP = 0.25
per raw relation, as at the nodes. It stops at 16 distinct relations: 2^60.88
trials of 2^67.97 instructions over about 2^37.2 lattices, tries them in turn
and fails only if all 16 fail.

## Nodes

A node is a special-q ideal of l bits on side f or g with its lattice L_q
(determinant q, LLL-reduced). Its region is the integer ball of radius R in the
reduced basis, one point of each +- pair; the exact half counts are in
`stats16.b85` and recomputed for R <= 8 by [7]. Two node types:

- **sieve**: every prime p in [2^8, 2^69] and every relative-degree-one root pair
  is enumerated on the ball by `dsieve16` (KERNELS.md); records (64-bit owner
  = index >> 40, 40-bit offset, 16 log2 p, one stream per side) go to the
  owners' 2^40-entry score arrays. Every primitive cell gets a certified upper
  bound on both norms: the 72 + 12 linear forms a(eta) - xi b(eta) are evaluated
  in double precision from the exact coordinates and the bound is prod (|L~| +
  e), e = 2^-46 (sum |v_j| |eta|^j + |xi| sum |v_6+j| |eta|^j)
  (`descent_witness.cell_bound`; ~1,500 instructions, charged 2,400).
  Cells whose residual masses (log2 bound minus score, minus log2 q on the
  special-q side) are at most T = kA nA + kB n + 96 on both sides survive, get exact norms, trial division and ECM.
  Non-primitive cells repeat the relation of c / gcd(c) and are skipped.
- **ecm**: ECM trials on every ball point, without a sieve (never cheaper in
  the shipped plan).

A relation has every residual prime <= 2^n per side, at most kA in (2^69, 2^nA]
and at most kB <= 1 in (2^nA, 2^n]; these window primes are the children, each
charged at the most expensive bin at or below its window.

## Yield and survivors

Per cell, the two-window semismooth probability

  P = sum_{jA <= kA, jB <= kB} 1/(jA! jB!) int rho((m - sum t + alpha)/69) prod dt/t

is evaluated as a lower sum: each mesh cell carries its exact weight int dt/t and
rho at its left endpoint (rho decreases), and cells crossing sum t = m are
dropped. The norm sizes m come from the measured offset bins (cf, cg) = (log2 N_f
- 72c, log2 N_g - 12c - 6 log2|u|) over 61,440 exact ball norms on 160
production special-q lattices (10 sizes from 72 to 300 bits, both sides, 6
radii): the cf sextile maxima, each paired with the overall cg maximum, shifted
by the worst sampled per-(ideal, radius) mean minus the pooled mean (cf 20.68 +
1.58 -> +23, cg -2.53 + 7.81 -> +6). Boot norms use 3,072 exact norms on 8
certified boot lattices (p90 + 4 and p99 + 4 bins). The relation count
is mu = cells x P x UDUP x detection x retention, with UDUP = 0.25 distinct
relations per raw relation, ECM detection 1 - 0.001 x (expected ECM factors:
(m_f + m_g)/20 + 2(kA + kB) at a node, m/20 + kA + 1 at the boot) and the
retention factor of the next section (the boot, which tests every point,
carries only the counted fraction). Survivors per cell use the
upper bound (T/69) rho((m - T)/69) per side, which decreases in m, so it is
evaluated at each bin's lowest norms: the previous sextile's raw cf maximum
(the observed minimum -140.05 for the first sextile) and the observed cg minimum
-54.56 (stats16.b85), both shifted down by the same margins (23, 6), with the
special-q at the bin's upper edge.

## Factor-base ideals without a logarithm

Under the main-phase occurrence model (`filtsim.py`: 2^64.162 relations, fewer than the 2^64.219 of
proof.md section 3 and so more missing logarithms, an ideal
of norm q = 2^(69 t) on side s occurs with probability 2 kappa_s(t)/q,
kappa_s(t) = rho(u_s - t)/rho(u_s), u_s from the relation-weighted mean norms
750.0 and 498.5 bits, rounded down; [4] checks them), an ideal's occurrences are
Poisson(lambda_s(b)); it is missing when it occurs at most twice (singletons and
clique pairs). A relation of an m-bit side carries on average sum_b kappa(m/69,
b/69) ln((b + 0.1)/b) P(Poisson(lambda_s(b)) <= 2) missing ideals (1/10-bit bands
b < 69; `descent16.missing`): 0.08-0.77 per f side and 0.25-2.25 per g side over
u = 7..39 (2.29 at the boot's 2,709-bit norm). Each node takes one option:

- **discard**: relations containing a missing ideal are rejected, multiplying
  mu by exp(-missing); nodes with q <= 2^95 always discard (the code forces
  it for bins up to 90 bits), so missing ideals never recurse;
- **descend**: each missing ideal Q is a further child, sieved in the lattice of
  Q q' (see Repeated ideals) and charged the largest work and failure probability
  of the discarding nodes with q <= 2^95 (Q q' <= 2^89; smaller Q are cheaper);
  the boot descends its missing ideals the same way.

Cascading removal can also lose an ideal with three or more occurrences; [6b]
measures the losses on the [6] system (purge, back-substitution to a fixpoint):
for lambda in [1,2), [2,3), [3,5), [5,8), >= 8 the measured missing fraction is
0.0470, 0.0309, 0.0224, 0.0036, 0 against the model's 0.7826, 0.5430, 0.2811,
0.0603, 0.0017 (signed system, seed 1); ideals lost with >= 3 occurrences: 0 of
3,816.

## Retention and cofactorisation

**Counted cells.** A cell is counted only if it is base-primitive: for every
prime p and irreducible factor h of H mod p, a and b are not both divisible by
(p, h). Any other cell has a - b theta in pO_K for a prime p of the base field
and repeats, up to the free relation of p, the relation of (a - b theta)/pi;
its norms carry an unsieved p^12 (f) and p^2 (g). The base-primitive fraction
prod (1 - p^-2 deg h) is 0.9100 (`retain16.base_primitive`: exact factor degrees
for p < 2^9, at most 6 factors above). Cells divisible by a relative-degree-one
ideal above a base prime of degree >= 2 with p >= 2^8 (outside the factor base)
are excluded by the ordinary fractions 0.9523/0.9919 (f/g, disjoint events). The factor base holds
every prime ideal of norm <= 2^69 above p < 2^8 (each has density >= 2^-48 in
main-phase relations, so its log is known) and the degree-one ideals above p in
[2^8, 2^69].

**Omitted mass.** A counted relation cell is rejected only when, on one side,
its omitted factor-base mass X exceeds s = 96 - 1 (bound gap, below 2 x 10^-11
bits, [11b]) - m/128 (rounding: each of the at most m/8 sieved hits of an m-bit
side loses below 1/16 bit to floor(16 log2 p)). X is a sum of independent local
contributions: at p < 2^8 the whole valuation, at p >= 2^8 the powers beyond the
first of every sieved degree-one ideal. Each contribution has one local bound.
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

**Relation-conditioned loss.** Under H1 the cofactor of a cell with omitted
mass x is semismooth with probability Pss(m - x) (the model of the previous
section at size m - x), so the fraction of relations lost on a side is

  loss = sum_{x > s} P(X = x) Pss_hi(m - x) / sum_{x <= s} P(X = x) Pss_lo(m - x),

with the upper sum of Pss (rho at the right endpoints, upper rho bracket,
boundary cells kept; `retain16.ss_upper2`) above and the lower sum below,
bracketed by monotonicity on a 64-bit grid in x (ceiling above s, floor
below), at the node's largest norms
m_f, m_g. The retention factor is (1 - loss_f - loss_g) x base-primitive x
(1 - the two ordinary exclusions, which are disjoint); [7] evaluates it for
every sieve node (minimum 0.8592; the largest
loss on either side is below 10^-5).

**ECM schedule (`ecm16.py`).** Montgomery curves in Suyama form (12 | #E), the
cheapest point of a (B1, B2/B1) grid with wheel 30030; stage 1 x-only ladder at
11 multiplications per bit, stage 2 baby/giant cross products at 3 per pair, a
binary gcd per stage. Per-curve success P1 is the Dickman lower sum on the curve
order corrected for non-powersmooth orders; curves = ceil(ln 0.001 / ln(1 - P1)),
each priced at 367.76 (L/5)^2 per multiplication. [12] prints the node and boot
schedules (boot: B1 = 2^28, 1430207 curves, 2^67.97 instructions per trial) and
matches P1 against detections on random 36- and 40-bit primes.

## Recurrence and success

Distinct relations at a node are negative binomial with mean mu and shape r:
r = 0.5 for effective windows of at most one prime (the A window counts kA only
when nA > 69, the B window kB only when n > nA; the base nodes (bin (69, 75))
have n = nA = 69, empty windows) and r = 1 otherwise. The shapes are set by the number of
ideals with no relation on toy fixed balls (offline fits; the empty-window rows are in `stats16.b85` and printed by [11b]): one-prime windows fit 0.60, two- or
three-prime windows 1.86-12.2. For empty windows (`stats16.b85` nb_empty_window,
12 configurations), every configuration with mean >= 7.4 has 0-4 ideals without
a relation against 5.4-12.4 predicted at r = 0.5; the three with means 3.2-4.0
have 19, 11 and 21 against 18.1, 17.2 and 16.4, and raising lb on the same
ideal sizes clears them (0 against 7.3-11.2). The base nodes run at mu = 454 and 508.

A child subtree fails with probability phi; phi(node) = max((1 + mu (1
- psi)/r)^-r, 0.02), psi = sum over the window of the children's phi, and a node
tries on average at most 1/(1 - psi) relations. The work of a node is its own
sieve/ECM work plus 1/(1 - psi) times its children's work. The plan requires
mu >= 16, phi <= 0.05 and psi <= 0.8 at every node ([7] checks all three).
Repeated ideals. Special-q children (bins of 69 bits and up): under the meter the
two-log run has at most 2^29.52 nodes (cap / the cheapest node's own cost, [7]);
placing all of them in the smallest bin (2^70.3 ideals) and allowing a factor 2^6
for the 1/p weighting inside a bin, two special-q nodes share an ideal with
probability <= 2^-6.26, charged in full. A factor-base ideal Q without a log
(norm <= 2^69) is never descended in its own fixed lattice: each visit draws a
fresh ideal q' of about 2^20 whose log is known (main-phase factor base) and
sieves the lattice of Q q', so a recurrence of Q samples a different sublattice.
Its node is charged as the worst node with q <= 2^95, and all nodes up to that
size discard relations touching an ideal without a log (no deeper recursion).

The shipped plan has 62 nodes (bins (69, 75) and 5-bit bins from 75 to 225 bits, both sides)
and the boot (n1 = 215, nA = 199.79, kA = 15, R = 4). Its entries were chosen by coordinate
search over (n, iA, kA, kB, R, discard) of every entry and of the boot, minimising the
modelled two-log work subject to the three conditions above: min mu 36.1, max phi 0.0441,
max psi 0.482 (boot 0.577); 28 nodes discard and 34 descend their missing ideals. The base
nodes sieve R = 320 (f, mu 454, phi 0.034) and R = 384 (g, mu 508, phi 0.032); all 62 nodes
sieve, with balls of R = 192..384. Boot plus subtree: one log 2^131.51 instructions; **two
logs 2^118.224 units**. The meter stops the descent at 2^118.975 units (Markov <= 2^-0.750; exact 2^-0.750795).
P(a boot fails) <= psi^16 = 0.57707^16 < 1.52e-04.

## Supporting measurements ([11])

Norm bins (offline, 20 lattices per size): fresh special-q lattices at 150 and 260
bits stay within the node bins (cf max 53.4 <= 109.6, cg max 5.19 <= 20.5), the
certified boot lattice within the boot bins (cg max 6.6 <= 16.7); [11] re-samples
a few points every run. Unit duplicates (offline): mates eps(a - bx) of 480 ball
points on 20 production lattices under 840 small units give distinct/raw 0.779;
the model charges 0.25 ([11] repeats 2 points).

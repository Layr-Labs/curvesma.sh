# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [8], [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T is replaced by T' = T^e, e uniform in [1, r), log T = log T'/e
mod r, so every boot sees a uniform element of mu_r. The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6; its LLL basis is
certified exactly (Gram determinant p^12, Lovasz conditions; [11]). The boot
element has a single (g-side) norm. The boot enumerates the radius-4 ball of the
lattice (2^23.67 points, one of each +- pair), computes each norm and runs the
`ecm16.py` schedule for 214.25-bit primes: a relation has every prime <= 2^235,
<= 12 primes in (2^69, 2^214.25] and <= 1 in (2^214.25, 2^235].
When a ball is exhausted it draws a fresh exponent e and a new certified
lattice (2^36 instructions each). Distinct relations are credited at UDUP = 0.25
per raw relation, as at the nodes. It stops at 16 distinct relations: 2^61.18
trials of 2^70.13 instructions over about 2^37.5 lattices, tries them in turn
and fails only if all 16 fail.

## Nodes

A node is a special-q ideal of l bits on side f or g with its lattice L_q
(determinant q, LLL-reduced). Its region is the integer ball of radius R in the
reduced basis, one point of each +- pair; the exact half counts are in
`stats16.b85` and recomputed for R <= 8 by [7]. Two node types:

- **sieve**: every prime p in [2^8, 2^69] and every relative-degree-one root pair
  is enumerated on the ball by `dsieve16.c` (KERNELS.md); records (64-bit owner
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

Under the main-phase occurrence model (`filtsim.py`: 2^64.162 relations, an ideal
of norm q = 2^(69 t) on side s occurs with probability 2 kappa_s(t)/q,
kappa_s(t) = rho(u_s - t)/rho(u_s), u_s from the relation-weighted mean norms
750.0 and 498.5 bits, rounded down; [4] checks them), an ideal's occurrences are
Poisson(lambda_s(b)); it is missing when it occurs at most twice (singletons and
clique pairs). A relation of an m-bit side carries on average sum_b kappa(m/69,
b/69) ln((b + 0.1)/b) P(Poisson(lambda_s(b)) <= 2) missing ideals (1/10-bit bands
b < 69; `descent16.missing`): 0.08-0.77 per f side and 0.25-2.25 per g side over
u = 7..39 (2.29 at the boot's 2,709-bit norm). Each node takes one option:

- **discard**: relations containing a missing ideal are rejected, multiplying
  mu by exp(-missing); nodes with q <= 2^80 always discard, so missing ideals
  never recurse;
- **descend**: each missing ideal is a further child, charged the largest work
  and failure probability of the discarding nodes with q <= 2^80 (the boot
  descends its missing ideals).

Cascading removal can also lose an ideal with three or more occurrences; [6b]
measures the losses on the [6] system (purge, back-substitution to a fixpoint):
for lambda in [1,2), [2,3), [3,5), [5,8), >= 8 the measured missing fraction is
0.083, 0.056, 0.020, 0.0036, 0 against the model's 0.78, 0.54, 0.28, 0.060,
0.0018; ideals lost with >= 3 occurrences: 3 of 3,830 (at 2^14, 4 of 15,265 and
10 of 15,279).

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
its omitted factor-base mass X exceeds s = 96 - 1 (bound gap, below 10^-11
bits, [11b]) - m/128 (rounding: each of the at most m/8 sieved hits of an m-bit
side loses below 1/16 bit to floor(16 log2 p)). X is a sum of independent local
contributions: at p < 2^8 the whole valuation, at p >= 2^8 the powers beyond the
first of every sieved degree-one ideal. At the 92 primes p < 2^9 where H, f and
g are squarefree modulo every prime of the base (`retain16.structure`:
squarefreeness and root counts over F_{p^d1}; the table in `stats16.b85` is
spot-checked by [11b]), a relative-degree-one ideal above
a base prime of degree d1 divides a base-primitive pair to order >= j with
probability n1 p^(-d1 j)/(1 + p^-d1), n1 the number of such ideals, and
relative-degree >= 2 ideals never divide one. At the other five primes (2, 3,
19, 23, 67) v_p(N_f) and v_p(N_g) are sampled over 2^18 uniform base-primitive
residue classes modulo p^K >= 2^104 (`stats16.b85` bad_prime_valuations,
re-sampled by [11b]). Above 2^9 the excess powers are bounded per dyadic block by 72 p^-2 events
per level over pi(x) <= 1.26 x/ln x primes (one- and two-event terms). Masses are rounded up to 1/4 bit: means 3.35 (f) and 7.36 (g) bits,
P(X > 96) = 2.2e-13 (f) and 7.0e-11 (g) per cell.

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

**ECM schedule (`ecm16.py`).** Survivor cofactors and boot norms are factored
with Montgomery curves in Suyama form (12 | #E). For (fbits, nbits), the
schedule is the cheapest point of the grid B1 = 2^10, 2^12, ..., 2^36, B2/B1 in
{1, 10, 30, 100, 300}, wheel w = 30030. Stage 1 multiplies the point
by E = prod_{q <= B1} q^floor(log_q B1) on the x-only ladder: at most 1.5 B1
bits (psi(B1) <= 1.03883 B1) at 11 multiplications per bit (xDBL 5, xADD 6).
Stage 2 computes the J = phi(w)/2 baby points [j]Q (j < w/2 coprime to w) and the
G = ceil(B2/w) + 1 giant points [iw]Q by differential additions (6 each; the
baby points through the chain of odd multiples below w/2) and,
for every (i, j) pair, the cross product x_i z_j - x_j z_i into an accumulator
(3 each); every prime in (B1, B2] is iw +- j for one pair. Each stage ends with
one binary gcd (<= 2 nbits iterations of 12 L + 16 instructions over L limbs);
a curve whose gcd returns the whole input is redone with a gcd after every
prime power and giant step, charged at the per-curve probability (the event
ends that input's factorisation). A curve finds p
with probability P1 = (rho(u) + int_{B1}^{B2} rho(u - t) dt/t)(1 - c), u = (fbits -
log2 12)/log2 B1, the lower sum of descent16, and c <= rho(u - 2)/rho(u) x
2/(sqrt B1 ln sqrt B1) the fraction of smooth orders with a prime above sqrt B1 to
a power >= 2 (stage 1 computes a B1-powersmooth multiple). The curve count is
ceil(ln 0.001 / ln(1 - P1)) for 0.999 per prime; the cost is curves x (1.1 x
multiplications x 367.76 (L/5)^2 + 2 gcd) plus the redo term. [12] prints the
schedules of a node and of the boot (B1 = 2^30, B2 = 30 B1, wheel 30030, 1595912 curves, 2^70.13 instructions per boot trial) and counts the detections of
the same stage-1/stage-2 procedure (`ecm16.attempt`, wheel 210) on random 36-
and 40-bit primes against P1 (within two binomial standard deviations).

## Recurrence and success

Distinct relations at a node are negative binomial with mean mu and shape r:
r = 0.5 for effective windows of at most one prime (the A window counts kA only
when nA > 69, the B window kB only when n > nA; the base nodes (bin (69, 75))
have n = nA = 69, empty windows) and r = 1 otherwise. The shapes are set by the number of
ideals with no relation on toy fixed balls: one-prime windows fit 0.60, two- or
three-prime windows 1.86-12.2. For empty windows (`stats16.b85` nb_empty_window,
12 configurations), every configuration with mean >= 7.4 has 0-4 ideals without
a relation against 5.4-12.4 predicted at r = 0.5; the three with means 3.2-4.0
have 19, 11 and 21 against 18.1, 17.2 and 16.4, and raising lb on the same
ideal sizes clears them (0 against 7.3-11.2). The base nodes run at mu = 297 and 539.

A child subtree fails with probability phi; phi(node) = max((1 + mu (1
- psi)/r)^-r, 0.02), psi = sum over the window of the children's phi, and a node
tries on average at most 1/(1 - psi) relations. The work of a node is its own
sieve/ECM work plus 1/(1 - psi) times its children's work. The plan requires
mu >= 16, phi <= 0.05 and psi <= 0.8 at every node ([7] checks all three).
Distinct relations of one node have distinct window primes, so their subtrees
descend disjoint ideals and meet only by drawing the same special-q among the
>= 2^63 ideals of a bin: subtree failures are independent (H5).

The shipped plan has 66 nodes (bins (69, 75) and 5-bit bins from 75 to 235 bits, both sides),
chosen bin by bin from the bottom by minimum work subject to phi <= 0.045, over
the radii and discard options of each v20 entry: min mu 25.1, max phi 0.0449, max psi 0.6785;
30 nodes discard and 36 descend their missing ideals. The base nodes sieve R = 288
(f, mu 278, phi 0.042) and R = 384 (g, mu 508, phi 0.031); all 66 nodes sieve, with balls
of R = 192..384 (bins above the boot window n1 = 235 are not needed). Boot plus subtree: one
log 2^133.75 instructions; **two logs 2^120.466 units**. The meter stops the descent at
2^122.5 units (Markov <= 2^-2.03). P(a boot fails) <= psi^16 = 0.679^16 = 2.02e-03.

## Supporting measurements ([8], [11])

Norm bins: fresh special-q lattices at 150 and 260 bits stay within the node
bins (cf max 53.4 <= 109.6, cg max 5.19 <= 20.5), the certified boot lattice
within the boot bins (cg max 6.6 <= 16.7). Unit duplicates: mates eps(a - bx) of
480 ball points on 20 production lattices under 840 small units give
distinct/raw 0.779; the model charges 0.25. Toy fixed balls: 6 configurations,
224 ideals, mean distinct relations / (model mean x 0.25) = 0.92-1.09 ([8]
repeats one ideal).

# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [7c], [7r],
[7s], [11], [11b] and [12] check its inputs. 1 unit = 20,000 instructions.
Inputs of this revision: LOG2_RELS = 63.68 (the sub-region pool estimate
minus two pool standard errors, section "Inputs"), mean norms 741.2/498.7 (at or
below the sub-region's relation-weighted means), UDUP = 0.25, kernel charges
SIEVE_HIT = 850 and BOUND_CELL = 2,000 (KERNELS.md).

## Targets and boot

Each target T is replaced by T' = T^e, e uniform in [1, r), log T = log T'/e
mod r: uniform on nonidentity mu_r (distance 1/r from uniform on all mu_r;
identity z handled directly; success charges 2/r). The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6. The boot
element a - b x has the single g-side norm N(a^2 + U y b^2), whose two terms
balance at |b| = |a|/sqrt(U); the basis is therefore LLL-reduced with the
b-block weighted by W = isqrt(U) and unscaled exactly (`boot_lattice`; the
lattice is unchanged, only the basis and hence the enumerated points change),
and certified exactly in the weighted metric: membership of every vector, size
reduction, Lovasz with delta = 0.98, Gram determinant p^12 W^12 ([11]). The boot
enumerates the radius-4 ball (2^23.67 points, one of each +- pair), computes
each norm and runs the `ecm16.py` schedule for 199.72-bit primes: a relation
has every prime <= 2^210, <= 13 in (2^69, 2^199.72], <= 1 in
(2^199.72, 2^210]. Norm offsets cg = log2 N_g - 12 c - 6 log2 U were
re-measured on 8 certified weighted lattices x 384 integer points drawn
uniformly from the ball: p50/p90/p99/max = -193.4/-187.9/-184.1/-181.0, charged
as the bins (0.5: -189.4, 0.4: -183.8, 0.1: -180.0) (quantile + 4; the top bin
covers the maximum). The unweighted basis measured the same way gives
-9.1/-1.3/3.6/7.3, about 10 bits below the bins 42e80f80 inherited from
6282ab79 (1.4/6.8/12.7 + 4, a different point sample); the weighting lowers
the boot norm by about 184 bits. An exhausted ball draws a
fresh exponent and certified lattice (2^36 each). Distinct relations credit
UDUP = 0.25, as at nodes (re-measured on the weighted lattices, below). It
stops at 16 distinct relations: 2^56.64 trials of 2^67.73
instructions over 2^32.98 lattices, tried in turn; failure needs all
16 to fail. Two kept relations come from the same lattice (same
exponent e) with probability at most C(16, 2) 2^-32.98 =
2^-26.07, charged as a failure ([7]).

## Nodes

A node is a special-q ideal (l bits, side f/g) with lattice L_q (det q,
LLL-reduced). Its region is the radius-R integer ball in the reduced basis,
one point per +- pair; exact half counts are in `stats16.b85` (recomputed for
R <= 8 by [7]). A **sieve** node enumerates every p in [2^8, 2^69] and every
relative-degree-one root pair on the ball (`dsieve16`, KERNELS.md); records go
to the owners' 2^40-entry score arrays; every primitive cell gets a certified
upper bound on both norms (72 + 12 linear forms in double precision, outward
rounding; ~1,500 instructions, charged 2,000; the forms' roots have certified
enclosures, [7r]); cells with residual mass (log2 bound minus score, minus
log2 q on the special-q side) <= T = kA nA + kB n + 96 on both sides survive
(exact norms, trial division, ECM); non-primitive cells repeat c / gcd(c)'s
relation and are skipped. A relation has every residual prime <= 2^n per
side, at most kA in (2^69, 2^nA] and at most kB <= 1 in (2^nA, 2^n]; these
window primes are the children, each charged at the most expensive bin at or
below its window.

## Yield and survivors

Per cell, the two-window semismooth probability P = sum_{jA <= kA, jB <= kB}
1/(jA! jB!) int rho((m - sum t + alpha)/69) prod dt/t is a lower sum (exact
mesh weights, rho at left endpoints, cells crossing sum t = m dropped). Norm
sizes m come from measured offset bins over 61,440 exact ball norms on 160
production special-q lattices (10 sizes 72-300 bits, both sides, 6 radii): cf
sextile maxima paired with the overall cg maximum, shifted by the worst sampled
per-(ideal, radius) mean minus the pooled mean (cf +23, cg +6); boot bins: see
Targets and boot. mu = cells x P x UDUP x detection x retention (ECM detection
1 - 0.001 x expected factors). Survivors per cell use (T/69) rho((m - T)/69)
per side at each bin's lowest norms (margins (23, 6), special-q at the bin's
upper edge).

**Unit duplicates (UDUP = 0.25).** Associates eps (a - b x), eps a unit of O_h,
give the same relation, so raw relations are credited at UDUP distinct
relations each. The charge is set by the relation-conditioned measurement of
41025889 (method and seeds there; reproducible from `descent_witness.py`): on
24 production special-q lattices (72, 95, 150 and 200 bits, both sides, three
seeds each) 300 ball points at R = 384 were stratified by norm quintile; each
quintile's distinct/raw (mean of 1/(1 + mates) under the 840 units with
coefficients in [-2, 2]^6) was weighted by the quintile's share of the
semismooth probability in the node's own window (relations concentrate in the
small-norm quintile, which carries 82-98% of the relation mass and has the most
mates). The relation-weighted distinct/raw is 0.283-0.751 over the 24 lattices
(mean 0.440, median 0.454), against 0.78-0.80 for the unconditioned ball mean
([11] repeats 8 points). Extending the unit set from 500 to 840 units changed the
unconditioned statistic by +0.008 (1a449d81), so larger units rarely map ball
points into the ball. On the weighted boot lattices the same relation-weighted
statistic (8 lattices x 300 uniform radius-4 points, weight rho(log2 N_g / 193.7)
per point, the 840 units) is 0.670-0.846 (unweighted
0.693 and up). The charge 0.25 is below every measured value (13% below
the smallest node value); the increase to 0.4 of 1a449d81 is not used.

## Factor-base ideals without a logarithm

Under the main-phase occurrence model (`filtsim.py`: 2^63.68 relations,
the sub-region estimate minus two standard errors), an ideal of norm q = 2^(69 t)
on side s occurs with probability 2 kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/
rho(u_s), u_s from the mean norms 741.2/498.7 (at or below the sub-region's
relation-weighted means; [4] checks them; smaller u gives fewer occurrences,
the conservative direction); occurrences are Poisson(lambda_s(b)). An ideal has
no logarithm when it never occurs (probability e^-lambda: it is absent from
every equation and no filtering step can recover it) or when it occurs but is
removed by the singleton/clique purge and not recovered by back-substitution
(a removed relation whose other ideals are known gives the removed ideal's
logarithm). The model charges, for 1 <= lambda < 8,

  P_miss(lambda) = (1 + 2^-8) e^-lambda + (1 - e^-lambda) c_r(lambda),

and the uncalibrated P(Poisson <= 2) for lambda < 1 and lambda >= 8
(`descent16.p_missing`). c_r covers the measured fraction of OCCURRING ideals
still without a logarithm after purge and back-substitution: sixteen
`filtsim.missing_recon` runs at this revision's means (seeds 1-3 at 2^12..2^16,
seed 1 at 2^17; excess 1.1; 413-29,768 occurring ideals per band and run)
give, per lambda band, maximum (mean) fractions [1, 2): 0.1141 (0.103); [2, 3): 0.0775 (0.064); [3, 5): 0.0362 (0.028); [5, 8): 0.0071 (0.006); >= 8: 0.0007, with no trend
over the scales; the covers c_r = 0.172 / 0.117 / 0.057 / 0.02 are the largest of 1.5 x the maximum, the
maximum plus three binomial standard deviations at the band's smallest sample,
and 1.3 x the value a hosted rerun of e178a7af measured at its reviewer's seed
(0.0149 in [5, 8), above that package's 0.011 cover), rounded up (at excess 1.2
the same measurement gave 0.0913 / 0.0584 / 0.0284 / 0.0078). [7c] reruns the
measurement at the reviewer's seed at 2^12 and reports each band against its
cover; the covers are charged above every measured value, not above a sampling
allowance. The calibration is
load-bearing and part of H5 ([7s] prints the plan at c_r = 1).
A relation of an m-bit side then carries on average sum_b kappa(m/69, b/69)
ln((b+0.1)/b) P_miss(lambda_s(b)) missing ideals (1/10-bit bands b < 69):
0.38-0.63 per f side and 0.36-0.57 per g side over the shipped nodes,
1.25 at the boot's 2512-bit norm. Each node takes one option:

- **discard**: relations containing a missing ideal are rejected, multiplying
  mu by exp(-missing); bins (69, 75] through (90, 95] always discard (forced in
  code: THIN_MAX = 90); bins from (95, 100] up use the searched flag;
- **descend**: each missing ideal Q (a prime ideal of norm <= 2^69) is a
  further child, sieved **in its own lattice L_Q as an ordinary prime special-q
  node**; its norms are at most those of the base bin (69, 75] on either side,
  so it is charged the larger base-bin work and failure probability of the two
  sides. Base-bin nodes discard, so the recursion closes. No composite Q q'
  lattices are used (the Q q' construction of f34f63d9 is withdrawn).

**The shared missing set.** Every relation of the descent faces the same
realized set M of ideals without a logarithm; a relation of norm size u
carries S(M) = sum_I w_I(u) 1[I in M] missing children, w_I(u) = kappa(u,
b_I/69) 2^-b_I. The recurrence uses E S = sum_I w_I m_I, m_I = P_miss(lambda_I);
three configuration-level bounds justify that (`descent16.miss_shared`, [7]):

- never-occurring ideals: S_0(M) = sum_I w_I 1[N_I = 0]. Under the
  multinomial occurrence model (the relations are n draws, ideal I drawn with
  probability p_I per relation) P(N_I = N_J = 0) = (1 - p_I - p_J)^n <= (1 -
  p_I)^n (1 - p_J)^n, so the zero-occupancy indicators of distinct ideals are
  pairwise negatively correlated; Chebyshev needs only this pairwise statement
  (no higher-order negative association, which a Chernoff bound would need):
  Var S_0 <= sum_I w_I^2 e^-lambda_I and
  P(S_0 > (1 + 2^-8) E S_0) <= Var S_0 / (2^-8 E S_0)^2 <= 2^-41.32 at every
  norm size of the plan (the per-ideal weights are below 2^-55, the sum has
  2^58 and more terms per band). The descent is evaluated at (1 + 2^-8) E S_0
  (the 2^-8 in P_miss) and 2^-30 is charged against success;
- occurring ideals without a logarithm: S_1(M)/E S_1 is a weighted average of
  the per-band fractions of 2^58 and more occurring ideals. Ideals occurring >= 3
  times are recovered by the shipped trap add-back (proof.md section 4: one original
  row of each still-unknown occ>=3 ideal is appended to the matrix, so linear algebra
  determines its logarithm; [6t] asserts none remain unknown), not by premise. H5c
  (premise, only for ideals occurring 1-2 times, the common mechanism): each realized
  band fraction is at most its cover c_r, which is 1.5 x the largest of sixteen
  measurements (taken with occ>=3 losses included, hence conservative for the occ<=2
  remainder) whose spread is at most 2.3 binomial standard
  deviations (largest |z| 3.7);
- repeated missing ideals: a missing ideal visited twice would repeat the same
  deterministic node, so the two visits are not independent trials. The ideal
  universe is the one of the occurrence model: per side li(B)/2 sigma-orbits of
  degree-one ideals (NBAR = 2 ideals per prime, quotiented by sigma), a relation
  of norm size u containing a given orbit of a 1/10-bit band with the
  band-average probability w_I(u) = kappa(u, b/69) ln((b + 0.1)/b) / n_b (the
  largest in the band is at most 2^0.1 x the average, since the weight varies
  as 1/q), so that sum_I w_I(u) m_I is exactly the per-relation missing count of
  the recurrence ([7u] checks the two agree to 1e-9 at every norm size of the
  plan). Two missing-log visits at relations of norm sizes u1,
  u2 are the same ideal with probability sum_I m_I w_I(u1) w_I(u2) / (E Z(u1)
  E Z(u2)) (the set M is shared: one factor m_I, not two) <= sum_I m_I wh_I^2 =
  2^-58.59, wh_I the largest share of a relation's missing children one ideal
  takes over the plan's norm range. The node counter (special-q and missing-log
  nodes, both logs) aborts past N_MAX = 2^25.5 nodes (a counted failure), so
  at most C(N_MAX, 2) visit pairs exist and the coincidence probability is at
  most 2^-8.58, charged in full. Ordinary special-q repeats (all N_MAX nodes in
  the smallest bin of 2^70.3 ideals, 2^6 for the 1/p weighting within a bin):
  2^-14.30.

## Retention and cofactorisation

**Counted cells.** A cell is counted only if base-primitive (for every p and
irreducible h of H mod p, a and b not both divisible by (p, h)); other cells
repeat (a - b theta)/pi's relation up to p's free relation. Base-primitive
fraction prod (1 - p^-2 deg h) = 0.9100 (exact degrees for p < 2^9). Cells
divisible by a relative-degree-one ideal above a base prime of degree >= 2
with p >= 2^8 are excluded by ordinary fractions 0.9523/0.9919 (f/g, disjoint).
The factor base holds every prime ideal of norm <= 2^69 above p < 2^8 plus the
degree-one ideals above p in [2^8, 2^69].

**Omitted mass and relation-conditioned loss.** A counted cell is rejected only
when one side's omitted factor-base mass X exceeds s = 96 - 1 (bound gap <
2e-11 bits) - m/128 (score rounding). X sums independent local contributions
(whole valuation at p < 2^8, excess powers at p >= 2^8), each with one local
bound: at a prime P of the base field with residue field F_q a base-primitive
cell has v_P >= j with probability at most min(1/q, mu q^-ceil(j/mu)) per root
of multiplicity mu (Hensel; q^-j for simple roots), independently over primes
(CRT); root multiplicities and ramification for p < 2^9 are tabulated
(`retain16.structure`, `stats16.b85`, recomputed by [11b]); above 2^9 excess
powers are bounded per dyadic block, with the remainders (`_REM`, `_DOWN`,
`block_omitted_total`) charged at 1.0 (retain16 docstrings). Means 3.35 (f)
and 7.67 (g) bits, P(X > 96) = 1.5e-12 (f) and 2.5e-10 (g) per cell. Under H1
a cell with omitted mass x is semismooth with probability Pss(m - x),
increasing in x; the lost fraction A/(A + B) (A = P(X > s, smooth), B = P(X
<= s, smooth)) is bounded by num/den with num = sum_{x > s} P_up(X = x)
Pss_hi(m - x) >= A and den = Pss_lo(m) (1 - P_up(X > s)) <= B
(`ss_upper2`/`ss_lower2`, 64-bit grid in x, the node's largest norms).
Retention = (1 - loss_f - loss_g) x base-primitive x (1 - ordinary exclusions).
The plan charges every sieve node the pinned floor RET_FLOOR = 0.8592 and
the boot the pinned counted fraction COUNTED = 0.8592; [11c] (part=losses)
recomputes the exact value at every shipped node (minimum 0.859250, largest
one-side loss < 1e-5) and the counted fraction (0.859251) and checks both.

**ECM schedule (`ecm16.py`).** Suyama-form Montgomery curves, cheapest (B1,
B2/B1) grid point (log2 B1 = 10..36, B2/B1 in 1/10/20/30/100/300) with wheel
30030; stage-1 x-only ladder (11 mults/bit), stage-2 cross products (3/pair),
one binary gcd per stage. Per-curve success P1 is the Dickman lower sum on the
curve order corrected for non-powersmooth orders; curves = ceil(ln 0.001 /
ln(1 - P1)) at 367.76 (L/5)^2 per mult. [12] prints node/boot schedules and
matches P1 on random 36/40-bit primes.

## Recurrence and success (H5)

Distinct relations per node are negative binomial (mean mu, shape r): r = 0.5
for effective windows of <= 1 prime (base nodes (bin (69, 75)) have empty
windows), r = 1 otherwise. Shapes come from toy-ball fits (empty-window rows
in `stats16.b85`, printed by [11b]): one-prime windows fit 0.60, two/three-prime
1.86-12.2; of 12 empty-window configs, all with mean >= 7.4 show 0-4 ideals
without a relation vs 5.4-12.4 predicted at r = 0.5.

A child subtree fails with probability phi; phi(node) = max((1 + mu (1 -
psi)/r)^-r, 0.02), psi = window sum of the children's phi; a node tries at most
1/(1 - psi) relations on average. Node work = own sieve/ECM work + 1/(1 - psi)
x children's work. The plan requires mu >= 16, phi <= 0.05, psi <= 0.8
everywhere ([7] checks). The three steps of this recurrence rest on different
grounds, stated separately:

1. **psi is a union bound.** A relation's subtree fails when any of its
   children (at most 2 kA + 2 kB window primes, both sides counted at the
   window maximum, plus its missing-log ideals) fails; P(union) <= sum of the
   children's failure probabilities, each taken at the worst bin at or below
   the child's window, with no assumption on how the children's failures
   depend on each other.
2. **Boot: the 16 kept relations are independent trials.** Each comes from
   its own ball of its own exponent e_i (shared balls: 2^-26.07, charged).
   A relation's subtree is a deterministic function of its child ideals, the
   factor-base logarithm table, and fresh ECM coins; the child ideals of the 16
   relations are functions of 16 independent uniform exponents. The table is
   shared, but it is fixed once the main phase has run: conditional on the
   table, the 16 subtrees are independent and P(all 16 fail) is the product of
   the 16 conditional failure probabilities. Each conditional probability is
   bounded by psi evaluated on the realized table, which the
   configuration-level bounds above control (S_0 within 2^-8 of its mean except
   with probability 2^-30; the occurring part within its cover under H5c); the
   charged boot failure carries a relative margin of 2^-16 above psi^16
   (`TABLE_SLACK`).
   Distinct subtrees that reach the same special-q or the same missing-log
   ideal are the one coupling that survives conditioning; both are charged as
   failures in full (2^-14.30 and 2^-8.58 above).
3. **Within a node (H5a, the residual heuristic).** The relations of one
   node are distinct lattice points of one ball; the thinning mu (1 - psi) and
   the 1/(1 - psi) retry count take their children's subtree successes as
   independent marks on the relation count: the independence of the
   factorisations of distinct cells that every NFS relation-count estimate uses
   (H1), applied to the children's descents. Measured support: (a) distinct
   relations of one node share children only through the repeat events charged
   above; (b) toy-ball relation counts fit the negative-binomial shapes used;
   (c) the norm bins are the worst sextile over 160 production lattices, so
   lattice-to-lattice yield variation is absorbed into mu; (d) the missing-log
   set enters a node only through deterministic rejection (bins up to (90, 95])
   or prime-lattice children charged at the base bin.

**Headroom and sensitivity.** The meter is set for modelled success 0.41, so
the expected two-log work may exceed its model value by 0.049 bits
(3.4%) before the success bound reaches 0.39 (boot failure up to
0.01181). The work recurrence compounds 1/(1 - psi) over the tree depth,
so E[W] is sensitive to psi; the charged shapes (0.5/1) are 17% and 46% below
the toy fits (0.60/1.86). [7s] prints the recurrence at the fitted shapes, at the
relation estimate and at one standard error, at UDUP = 0.283 and at c_r = 1.

The node counter (above) aborts past N_MAX = 2^25.5 nodes (counted failure);
E[N] <= 2^18.12 (`counts16.py`, missing-log children at the base bin), so
P(abort) <= 2^-7.38 (Markov).

## Inputs ([4], [7s])

LOG2_RELS = 63.68 is the sub-region pool estimate 2^63.8903 (proof.md
section 3) minus two pool standard errors (2 x 0.1014 = 0.2028 bits), floored to
2 decimals; [4] checks log2_rels - 2 SE >= 63.68. The same relation count feeds the
missing-log rates, the descent work and the success bound, so the
two-standard-error allowance is propagated through all three. Without the
calibration (c = 1) the f34f63d9 plan at these inputs and prime-lattice children
has two-log work 2^117.97 (offline), about 1.6 bits above the shipped plan: the
calibration is what makes the shipped meter possible (offline: the C1c plan
at c_r = 1 is infeasible).

## The shipped plan

62 nodes (bins (69, 75) plus 5-bit bins 75-225, both sides) and the boot (n1
= 210, nA = 199.72, kA = 13, R = 4, 16 relations), chosen by coordinate search
over (n, iA, kA, kB, R, discard) of every entry and of the boot (singles, pairs,
(R, n, kA) triples), minimising the metered cap E - log2(success budget) at
success 0.40 under the three plan conditions (the shipped meter is then set
for 0.41), from the plan of e178a7af (itself from 273b1519's, searched from
f34f63d9's UDUP = 0.25 plan via 7c1acb81), repaired and re-searched at this
revision's inputs (27 accepted moves; offline driver). Quadrature: mesh 384, rho
spacing 1/1024 (left-endpoint lower sums). At the inputs above: min mu
26.2, max phi 0.0470, max psi 0.673 (boot 0.673); all 62 nodes sieve at
R = 192..384. One log <= 2^128.659 instructions; **two logs <= 2^115.3704
units**. The meter stops the descent at 2^116.162 units (Markov: E[W]/cap <=
0.577703). P(a boot fails) <= 0.001763641 ([7]). Success >= 1 - 2 x 0.001763641 -
0.577703 - 2^-7.38 - 2^-14.30 - 2^-8.58 - 2^-26.07 - 2^-30 - 2^-48 - 2^-191
- 2/r >= 0.410103 (2^-191: Block-Wiedemann plus main phase, proof.md section 7).

## Supporting measurements ([11], [11b], [12])

Fresh lattices at 150/260 bits stay in the node bins (offline, 20 per size);
[11] re-samples points per run; [11b] recomputes the prime-ideal table, residue
valuation samples, certified cell bounds and the empty-window toys; [12] the ECM
schedules and toy detections.

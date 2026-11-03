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
UDUP = 0.65, as at nodes. It stops at 18 distinct relations: 2^58.47 trials of
2^67.96 instructions over ~2^34.8 lattices, tried in turn; failure needs all
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
(p90 + 4, p99 + 4 bins). For pr the sextiles are subdivided (`QBINS_SUB`:
24 cf-quantile cells, per-cell cg maxima (fresh uniform grid), same +23/+6
(top edge pinned); `BOOT_SUB`: 10 R=4 deciles +4 (edge pinned); charged
boost capped at 4x nodes / 1.5x boot. Sub-cell maxima dominate their cells (lower bound kept); fresh
gaps 19.14/3.78 re-validate the margins. mu = cells x P x UDUP x detection x retention (UDUP
= 0.65; ECM detection 1 - 0.001 x expected factors; retention from the next
section; the boot carries only the counted fraction). Survivors per cell use
(T/69) rho((m - T)/69) per side, decreasing in m, evaluated at each bin's
lowest norms (previous sextile's raw cf max, observed cg min -54.56, margins
(23, 6); special-q at the bin's upper edge).

## Factor-base ideals without a logarithm

Under the main-phase occurrence model (`filtsim.py`) with a **conservative
relation-count input**: LOG2_RELS = 63.883 = the theta = 9/25 retained-count
estimate 2^64.0641 minus two pool standard errors (0.0903 bits each), floored to
0.001 bit (certified lower product 0.971088853 already applied; [4] asserts
LOG2_RELS <= estimate - 2 SE). The same input drives every missing-log rate, the
expected descent work, the boot and node failure probabilities and the success
bound, so the yield uncertainty is propagated through all of them; the optional
run `extra = sensitivity` ([7s]) re-evaluates the whole plan at estimate - 0 and
- 1 SE (two-log 2^113.5863 / 2^113.6036, success 0.4126 / 0.4053) against the
declared estimate - 2 SE (2^113.6317, success 0.39329): work falls and success
rises with the relation count. An ideal of norm q = 2^(69 t) on side s occurs
with probability 2 kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/rho(u_s), u_s from
mean norms 743.2/498.9 (relation-weighted 743.27/498.93, rounded down; [4]).
Occurrences are Poisson(lambda_s(b)); an ideal occurring <= 2x is missing. A
relation of an m-bit side carries on average sum_b kappa(m/69, b/69)
ln((b+0.1)/b) P(Poisson <= 2) missing ideals (1/10-bit bands b < 69):
at the shipped input and plan's node norms 0.39-0.88 per f side (1,057-2,044
bits) and 0.34-0.58 per g side (495-730 bits); 2.962 at the boot's 2,709-bit
norm. Each node takes one option:

- **discard**: relations containing a missing ideal are rejected, multiplying
  mu by exp(-missing); nodes with q <= 2^90 always discard (the code forces
  it for bins up to 85 bits; bin 90 keeps its searched thin flag), so missing
  ideals never recurse there;
- **descend**: each missing ideal Q (a prime ideal of norm <= 2^69) is a further
  child, descended as an **ordinary prime special-q node in its own lattice L_Q**
  (det N(Q), the production `specialq_lattice` construction; no composite
  lattices) and charged the largest work and failure probability of the
  discarding nodes with q <= 2^90. This charge is conservative: with c = (l +
  S)/12, a node's norms are 72c + cf - l and 12c + LOG2U6 + cg (f-side special-q)
  or 72c + cf and S + LOG2U6 + cg (g-side), all non-decreasing in l, so L_Q with
  N(Q) <= 2^69 has norms no larger than the (69, 75] bin's and its relations
  discard missing ideals (recursion closes). [11] samples a 66-bit g-side L_Q:
  its exact ball norms stay inside the node bins (cf max 19.4 <= 109.6, cg max
  2.8 <= 20.5). The boot descends its missing ideals the same way. A second visit
  of the same Q would re-enumerate the same ball, so every repeat is charged as a
  failure ([7r], next section).

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

**Root enclosures and the bound gap.** The 84 linear forms L_k(v) = a(eta_k) -
xi_k b(eta_k) use roots with certified enclosures (`descent_witness.forms`, [11b]):
the six eta are `embedding_bounds.root_disks()` (exact rational centres, radius
2^-100, Rouche test, the same certificate the collection bounds use); for each eta
the 12 roots xi of f(eta, x) are refined to 60 digits (Decimal Durand-Kerner) and
enclosed by Weierstrass disks D(z_i, 12 |W_i|), W_i = (|f(z_i)| + delta(z_i)) /
prod_{j != i} |z_i - z_j|, where delta bounds the change of f's coefficients over
the eta disk; all 12 disks per eta are pairwise disjoint (one root each). g side:
xi = +-sqrt(-u eta), error |U| r_eta/|xi|. Largest relative radius over the 90
roots: 1.8e-28. The kernel evaluates with the doubles nearest the certified
centres (relative error <= 2^-53 + 1.8e-28), so the root error of L_k is <= 6.001
x 2^-53 T_k and the double evaluation error <= 12 x 2^-53 T_k, where T_k(v) = sum
|v_j||eta|^j + |xi| sum |v_6+j||eta|^j; both lie inside the charged e_k = 2^-46 T_k,
so prod (|L~_k| + e_k) >= |N(v)| for every cell (validity; [11b] checks it with
zero tolerance on sampled cells). The gap is uniform except near cancellation:
log2((|L~_k| + e_k)/|L_k|) <= log2(1 + 2^-45 T_k/|L_k|), so whenever every form has
|L_k|/T_k >= 2^-38 the gap is <= 72 log2(1 + 2^-7) = 0.81 bits on the f side (0.14
on g) inside the 1 bit charged below. Cells with some |L_k|/T_k < 2^-38 are
**exceptional** and charged as lost: measured on 2,000 uniform cells of five
production special-q lattices (75-220 bits), P(min |L|/T < eps)/eps is at most 8.4
(f, complex forms, falling with eps) and 5.0 (g, real forms, flat) for eps in
[1e-4, 0.1], so the exceptional fraction is <= 10 x 2^-38 ~ 2^-34.7 per side by
linear extrapolation; the model multiplies every sieve node's relation count by 1 -
2^-20 (`EXC_CELLS`), a margin of 2^14 over that estimate (H5).

**Omitted mass.** A counted cell is rejected only when one side's omitted
factor-base mass X exceeds s = 96 - 1 (bound gap, previous paragraph) - m/128
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
7.3-11.2). Base nodes run at mu = 1279 (f) and 1081 (g).

A child subtree fails with probability phi; phi(node) = max((1 + mu (1 -
psi)/r)^-r, 0.02), psi = window sum of children's phi; a node tries <=
1/(1 - psi) relations on average. Node work = own sieve/ECM work + 1/(1 - psi)
x children's work. The plan requires mu >= 18, phi <= 0.05, psi <= 0.8
everywhere ([7] checks). The descent keeps a special-q node counter over both
logs and aborts past N_MAX = 2^25 nodes (counted failure); E[N] <= 2^16.17
(`counts16.py`, unrounded 2^16.1652), so P(abort) <= 2^-8.83 (Markov). All N_MAX
nodes in the smallest bin (2^70.3 ideals) with 2^6 for 1/p weighting:
shared-special-q probability <= 2^-15.30, charged in full ([7]).

**Missing-log repeats ([7r]).** Only relations tried at the 33 descending nodes and
the two boots carry missing ideals; each tried relation succeeds (one per node)
or fails through >= 1 failed child, so <= 2 N_MAX node and 36 boot relations. A relation with side size
u = m/69 holds a given missing ideal of norm 2^b with probability <= 2 kappa(u,
b/69) 2^-b P(Poisson(lambda_b) <= 2) (band lower edge, largest u of its class:
2^-61.57 at nodes, 2^-60.30 at the boot) and holds <= 1.316 (nodes) / 2.962 (boot)
missing ideals in expectation; for two relations sum_I p_I(r1) p_I(r2) <=
min(p(r1) M(r2), p(r2) M(r1)). The union bound over node-node, node-boot and
boot-boot pairs gives P(some missing ideal is visited twice) <= 2^-10.16, charged
in full as failure. (N_MAX 2^25 balances this term against the counter's Markov term.)

**H5: conditional independence of subtree failures.** The recurrence multiplies
child failure probabilities (negative-binomial thinning by 1 - psi, boot failure
psi^18). Condition on the main phase, i.e. on the fixed set M of factor-base
ideals without a logarithm (fixed by region, polynomials and filtering; independent of the challenge
and the descent coins). Given M, distinct subtrees are rooted at distinct prime ideals (all
special-q and missing-log repeats are excluded at a charged cost of 2^-15.30 +
2^-10.16), each sieves its own lattice ball with its own fresh coins, and each
relation's children are distinct primes, so on the no-repeat event their sieve
outcomes are functions of disjoint lattices; the premise is that these outcomes
behave as independent draws with the measured norm and yield laws (H1 per lattice,
as in standard NFS descent analyses). The dependence the reviewers raised runs
through M as a common cause (two subtrees fail together when M is unfavourable):
by Jensen, E[psi(M)^18] >= E[psi(M)]^18, so the unconditional boot failure can
exceed the conditional model at the mean M. The model is therefore evaluated at a
conservative quantile of the common cause, the estimate - 2 SE relation count,
which fixes the missing-log rates at their pessimistic values (higher psi, higher
work), and [7s] shows the success bound holds across estimate - 0/1/2 SE. Residual
dependence between distinct prime lattices given M is not measured at production
size and is a declared part of H5.

The shipped plan has 62 nodes (bins (69, 75) plus 5-bit bins 75-225, both
sides) and the boot (n1 = 215, nA = 199.79, kA = 15, R = 4), chosen by
coordinate search over (n, iA, kA, kB, R, discard) minimising modelled two-log
work subject to the three conditions (our 57598408 lineage's re-searched
plan). At this package's
estimate - 2 SE inputs (LOG2_RELS = 63.883, mean norms 743.2/498.9, UDUP = 0.65)
with THIN_MAX 85: min mu 38.5, max phi < 0.05, max psi 0.438 (boot); all 62 nodes
sieve. One log <= 2^126.919 instructions; **two logs <= 2^113.6317 units**
(unrounded 2^113.63169). Meter stops the descent at 2^114.36 units (Markov:
E[W]/cap <= 0.603613). P(boot fails) <= 3.4635e-07 ([7]). The plan is also
feasible without the integration and re-search: our validated 57598408 plan at
LOG2_RELS = 63.883 gives two-log 2^115.90, boot psi 0.656, phi 5.1e-4 (offline
probe), so the integration is a cost lever, not a feasibility premise.

## Supporting measurements ([11])

Norm bins (offline, 20 lattices/size): fresh lattices at 150/260 bits stay in
the node bins (cf max 53.4 <= 109.6, cg max 5.19 <= 20.5), the certified boot
lattice in the boot bins (cg max 6.6 <= 16.7); [11] re-samples a few points per
run. Unit duplicates (offline, relation-weighted, uniform cells): distinct
yield CELLS x pr x 0.65 needs E[shat/m] >= pr x 0.65 per regime (shat =
model per-point semismooth prob at exact cofactor norms, m = mates among
1532 units). Cells are drawn by `int_ball_point` (rejection from the
uniform box: exact uniformity over enumerated cells by construction;
uniformity validated: R=1 chi2 5.8, R=4 shell chi2 13.8; the old
round-and-reject sampler scored 506.4 and is retired from this study).
The campaign measured yr0 = E[shat/m]/pr0 against the UNINTEGRATED sextile pr0
(stats16.b85): 13 node regimes (bits 69-225, both sides, R 192-384; 750 pts x
1-3 lattices) 5.18-150, min CI lower 4.70; R=4 boot (3 true lattices, 1500 pts)
1.10-1.25, CI lower 1.06. The model charges the integrated pr <= cap x pr0, so
against the CHARGED pr the supported ratios are only yr >= yr0/cap: nodes >=
4.70/4 = 1.175, boot >= 1.06/1.5 = 0.707, both above UDUP 0.65 (boot by 8.7%);
no larger ratio against the charged pr is claimed. Truncation 86-100%
saturated, x5 (<=11.5%)/bound-8 closed. Run `extra = relw` repeats both live (8 + 8 points: node 1.34, boot 0.85
at seed 1) plus the R = 4 shell chi^2 of `int_ball_point`; live thresholds 1.0 / 0.72 against the charged pr (both above UDUP 0.65).

Quadrature: mesh 384, rho spacing 1/1024; plan iA tuples are 8x the physical windows.

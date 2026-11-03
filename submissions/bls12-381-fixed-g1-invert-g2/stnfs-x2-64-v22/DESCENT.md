# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T is replaced by T' = T^e, e uniform in [1, r): uniform on nonidentity mu_r (distance 1/r; identity z handled directly; ledger charges 2/r). The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6. The boot
element has a single (g-side) norm N(a^2 + U y b^2), U = |u| = 2^63.71; the basis is LLL-reduced
in the weighted metric with the b-block scaled by W = isqrt(U) = 2^31.86
(`boot_lll`), which balances the two terms, and returned unscaled. The certificate is exact: membership,
size and Lovasz (0.98) conditions in the weighted metric, Gram determinant
p^12 W^12 ([11], `certify_boot`). Measured on 8 certified lattices x 384 uniform
R = 4 points: g-norm offsets cg mean -191.8, p90 -186.0, p99 -182.0, max -178.6 (Euclidean basis: mean -7.3, max 12.6; spread 0.53). The boot
enumerates the radius-4 ball (2^23.67 points, one per +- pair) and runs the `ecm16.py` schedule for 199.79-bit primes: a relation
has every prime <= 2^215, <= 15 in (2^69, 2^199.79], <= 1 in (2^199.79, 2^215].
An exhausted ball draws a fresh exponent and certified lattice (2^36 each).
Distinct relations credit UDUP = 0.65. It stops at 18 distinct
relations (BOOT_RELS): 2^52.37 trials of 2^67.75 instructions over ~2^28.7
lattices; failure needs all 18 to fail.

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
  Non-primitive cells are skipped.
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
7.81 -> +6). Boot norms use 3,072 exact norms on 8 certified skewed boot lattices at R = 4 (p90/p99/max + 4 bins). The semismooth
probability integrates over 24 norm sub-bins at nodes and 10 at the boot (`QBINS_SUB`/`BOOT_SUB`),
charged at most 4x / 1.5x the sextile-bin value (UDUP support below is stated against this charged pr). mu = cells x P x UDUP x detection x retention (UDUP
= 0.65; ECM detection 1 - 0.001 x expected factors; the boot carries only the counted fraction). Survivors per cell use
(T/69) rho((m - T)/69) per side, decreasing in m, evaluated at each bin's
lowest norms (previous sextile's raw cf max, observed cg min -54.56, margins
(23, 6); special-q at the bin's upper edge).

## Factor-base ideals without a logarithm

Under the main-phase occurrence model (`filtsim.py`: 2^63.716 relations -- the
conservative descent input, the excess-1.1 theta = 2284/8000 retained-count
estimate 2^63.819798 minus 2 pool SE (0.0515 each; [7s] checks the floor),
certified lower product 0.967583931 already applied), an ideal of norm q = 2^(L_s t) on side s (L_f = 69.25,
L_g = 69, the side's factor-base bound) occurs
with probability 2 kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/rho(u_s), u_s = m_s/L_s from
mean norms 737.7/498.4 (rounded down from 737.70/498.47; [4] checks them and the bounds).
Occurrences are Poisson(lambda_s(b)); an ideal occurring <= 2x is missing, and an ideal
occurring >= 3x is charged lost at its band's rate c (below). A
relation of an m-bit side carries on average sum_b kappa(m/69, b/69)
ln((b+0.1)/b) P_miss(lambda_b) missing ideals, P_miss = P(Poisson <= 2) + c (1 - P(Poisson <= 2))
(`p_missing`; 1/10-bit bands b < 69):
0.19-1.71 per f side, 0.38-3.43 per g side over u = 7..39 (3.49 at the boot's
2,709-bit norm). The <= 2 part applies with coefficient 1, so every band charges at least the e^-lambda
zero-occurrence mass ([6b-zero] checks P_miss >= P(Poisson <= 2) >= e^-lambda on a lambda grid).
The descent keeps its window base at 2^69 on both sides, so f-side primes in (2^69, 2^69.25] are
treated as descent children rather than looked up (more work and failure than the main-phase
logarithms they have: conservative). Each node takes one option:

- **discard**: relations with a missing ideal are rejected (mu x exp(-missing)); nodes with q <= 2^90 always discard, so missing
  ideals never recurse there;
- **descend**: each missing ideal Q is a further child, sieved in the lattice of
  Q q' and charged the largest work and failure probability
  of the discarding nodes with q <= 2^90;
  the boot descends its missing ideals the same way.

Cascading removal can also lose ideals with >= 3 occurrences (e.g. an ideal
in exactly 3 rows each holding a different singleton is purged entirely); [6b] measures losses on the [6] system
(per-side bounds, shipped means): per lambda band
(ideals, measured, model, lost with >= 3 occurrences): [1,2): 934, 0.0942, 0.7690, 0; [2,3): 948, 0.0559, 0.5608, 2; [3,5): 765, 0.0379, 0.2748, 1; [5,8): 456, 0.0000, 0.0618, 0; [8,1e+09): 656, 0.0000, 0.0014, 0;
ideals lost with >= 3 occurrences: 3 of 3759; the model exceeds measurement by
>= 7.24x in every band with loss. The >= 3-occurrence losses are charged
explicitly: c is, per lambda band, the Wilson 95% upper bound (rule of three at
zero losses) of the pooled fraction of >= 3-occurrence ideals left without a
logarithm after cascading purge and back-substitution over three toy seeds at
2^12 on the same system (`H3_RATE`; [6c] recomputes the pooled counts and the
rates: band (ideals with >= 3 occurrences, lost, rate) [1,2): 758, 3, 0.01158; [2,3): 1333, 8, 0.01180; [3,5): 1801, 6, 0.00725; [5,8): 1292, 1, 0.00438; [8,1e+09): 2002, 0, 0.00150; pooled 18 of 7186).
The charge adds at most 0.0939 of the P(Poisson <= 2) count (worst over the
shipped u grid, both sides; refinement.json `h3_charge`). Production-regime
premise (H5n): such a loss needs a closed unknown-subsystem (minimal: 3 rows, 4
unknowns), whose formation is coincidence-driven (~1/U per constraint), so
the per-ideal loss rate does not rise from toy (2^12) to production scale; this
is a disclosed heuristic, not a measured production rate.

## Retention and cofactorisation

**Counted cells.** A cell is counted only if base-primitive (for every p and
irreducible h of H mod p, a and b not both divisible by (p, h)); other cells
repeat another cell's relation. Base-primitive
fraction prod (1 - p^-2 deg h) = 0.9100 (exact degrees for p < 2^9). Cells
divisible by a relative-degree-one ideal above a base prime of degree >= 2
with p >= 2^8 are excluded by ordinary fractions 0.9523/0.9919 (f/g, disjoint).
The factor base holds, per side s with bound B_s (B_f = 2^69.25, B_g = 2^69), every prime ideal of norm <= B_s above p < 2^8 plus the
degree-one ideals above p in [2^8, B_s]; the descent sieves p in [2^8, 2^69] on both sides.

**Omitted mass.** A counted cell is rejected only when one side's omitted
mass X exceeds s = 96 - 1 - m/128
(rounding: each sieved hit loses < 1/16 bit to floor(16 log2 p)). X sums independent local contributions (whole valuation at
p < 2^8, excess powers at p >= 2^8), each with one local bound.
Let P = (p, h_i) be a prime of the base field (Kummer-Dedekind, O = Z[y]) with
residue field F_q, q = p^deg h_i, and S the side polynomial, monic in x. For a
base-primitive pair, (a : b) is a point of P^1(O_P): b is a unit with
probability q/(q + 1), otherwise a is one and S_hom(a, b) is a unit. With b a
unit, v_P(S_hom(a, b)) = v_P(S(t)) for t = a/b uniform in O_P. If t does not
reduce to a root of S in F_q, the valuation is 0. A root rho of multiplicity mu
in F_q lifts to exactly mu roots in the completion (Hensel), so v_P(S(t)) is the sum of the mu distances v_P(t - rho_i),
and v_P(S(t)) >= j forces some v_P(t - rho_i) >= j/mu, an event of probability
at most q^-ceil(j/mu):

  P(v_P >= j and t reduces to rho) <= min(1/q, mu q^-ceil(j/mu)),

with equality q^-j when mu = 1. The primes above p are
independent (CRT) and v_p(N) = sum over P of deg(h_i) v_P. Degree,
ramification and f/g root multiplicities of every P
with p < 2^9 are tabulated (`retain16.structure`, `stats16.b85`); [11b]
recomputes the table at repeated-root/ramified primes (2, 19,
23) and a seventh of the others, and checks the bound against residue-class
samples at 2 and 23. All multiplicities are 1 except g at the two primes above 2
(2^16 | U: the double root x = 0, q = 4 and 16) and f at the degree-4 prime
above 19 (two double roots among four, q = 19^4); 23 is ramified in the base
field (e = 2) with g having two simple roots at each of the two degree-one
primes. Above 2^9 the excess powers are bounded per dyadic block by 72 p^-2
events per level over pi(x) <= 1.26 x/ln x primes.
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
B2/B1) grid point with wheel 30030; stage-1 x-only ladder,
stage-2 cross products, one binary gcd per stage. Per-curve success
P1 is the Dickman lower sum on the curve order; curves = ceil(ln 0.001 /
ln(1 - P1)) at 367.76 (L/5)^2 per mult. [12] prints node/boot schedules (boot:
B1 = 2^28, 1418167 curves, 2^67.96/trial) and matches P1 on random 36/40-bit
primes.

## Recurrence and success

Distinct relations per node are negative binomial (mean mu, shape r): r = 0.5
for effective windows of <= 1 prime (base nodes have empty
windows), r = 1 otherwise. Shapes come from toy-ball fits (empty-window rows
in `stats16.b85`, printed by [11b]): one-prime windows fit 0.60, two/three-prime
1.86-12.2. Of 12 empty-window configs, all with mean >= 7.4 show 0-4 ideals
without a relation vs 5.4-12.4 predicted at r = 0.5; the three with means
3.2-4.0 show 19/11/21 vs 18.1/17.2/16.4, and raising lb clears them (0 vs
7.3-11.2). Base nodes run at mu = 613 and 470.

A child subtree fails with probability phi; phi(node) = max((1 + mu (1 -
psi)/r)^-r, 0.02), psi = window sum of children's phi; a node tries <=
1/(1 - psi) relations on average. Node work = own sieve/ECM work + 1/(1 - psi)
x children's work. The plan requires mu >= 18, phi <= 0.05, psi <= 0.8
everywhere ([7] checks). The descent keeps a node counter over both
logs (cap N_MAX = 2^27.3, counted failure); E[N] <= 2^16.47
(`counts16.py`, unrounded 2^16.4611), so P(abort) <= 2^-10.83 (Markov). All N_MAX nodes in the
smallest bin (2^70.3 ideals) with 2^6 for 1/p weighting: shared-ideal
probability <= 2^-10.70, charged in full ([7]).
A factor-base ideal Q without a log
(norm <= B_s) is never descended in its own fixed lattice: each visit draws a
fresh ideal q' of about 2^20 with known log and
sieves the lattice of Q q'. Its node is charged as the worst node with q <= 2^90, and all nodes up to that
size discard relations touching an ideal without a log.

The shipped plan has 62 nodes (bins (69, 75) plus 5-bit bins 75-225, both
sides) and the boot (n1 = 215, nA = 199.79, kA = 15, R = 4).
The 90d99567 plan was re-optimised by single-then-joint coordinate
search at the sub-bin pr, the skewed boot and LOG2_RELS = 63.599:
41 of 62 entries changed (refinement.json `plan_v21`); the boot window is
unchanged. This revision re-evaluates that plan unchanged at the new inputs
(per-side bounds, LOG2_RELS = 63.716, the >= 3-occurrence charge).
At the inputs (LOG2_RELS = 63.716, mean norms 737.7/498.4,
UDUP = 0.65) with THIN_MAX 85: min mu 38.5, max phi < 0.05, max psi 0.483. All
62 nodes sieve. One log <= 2^125.407 instructions;
**two logs <= 2^112.119 units**. Meter stops the descent at 2^113.119 units
(Markov: E[W]/cap <= 0.5). P(boot fails) <= 1.000e-03 ([7]; charged
499x the modelled 2.01e-06). Spot cells at re-optimised nodes ([11] method,
40 cells, units in [-2, 2]^6; spot checks, not the campaign bound):
E[shat/m]/pr_charged 69f 7.65, 90f 6.18, 105g 9.7, 125f 8.17, 140g 3.51, 180g 3.33, all >> 0.65.

## Supporting measurements ([11])

Norm bins (offline, 20 lattices/size): fresh lattices at 150/260 bits stay in
the node bins (cf max 53.4 <= 109.6, cg max 5.19 <= 20.5); [11] re-samples per run. Unit duplicates (offline, relation-weighted, uniform cells): distinct
yield CELLS x pr x 0.65 needs E[shat/m] >= pr x 0.65 per regime (shat =
model per-point semismooth prob, m = mates among
1532 units). Cells are drawn by `int_ball_point` (rejection from the
uniform box, exact by construction;
validated: R=1 chi2 5.8, R=4 shell chi2 13.8; the old
round-and-reject sampler scored 506.4 and is retired).
yr = E[shat/m]/pr over 13 node regimes (bits 69-225, both sides, R
192-384; 750 pts x 1-3 lattices): 5.18-150, min CI lower 4.70 against the
sextile pr, hence >= 1.175 against the charged pr (4x cap). R=4 skewed boot
over 3 certified lattices (1000 pts each): 0.99 [0.95, 1.03] / 1.02 [1.00, 1.06] / 1.05 [1.02, 1.08], pooled 1.02 [1.00, 1.04]
(3000 points; against the charged pr, boost 1.16x: 0.89 [0.87, 0.91]; saturation 100% / 100% / 100%;
x5-mates 0.99 / 1.02 / 1.05). Complete orbit ([11u]): the unit group of O_h is certified
(PARI bnfinit/bnfcertify: class number 1); the units
that can map a cell into the ball are enumerated
completely and tested exactly. All 3000 skewed-boot campaign cells recounted:
exact yr 1.019 [1.002, 1.034] against the sextile pr (0.886 [0.871, 0.899] against the
charged pr), so E[shat/m] >= 0.87x (boot) the charged pr, above UDUP 0.65.
Node sample (120 cells, 70-bit f and 225-bit g at R = 288): the exact count
exceeds the bound-5 box at 1 cell; under unit-geometry/smoothness
independence (H5) the Wilson 95% upper diff rate 0.0457 gives haircut 0.9543: E[shat/m] >= 4.70/4 x 0.9543 = 1.12x (nodes)
the charged pr, above 0.65. [11] repeats the box counts live on fresh seeds
(thresholds 1.175 / 0.70) plus a sampler smoke check.

## Refined integration representation

Mesh 384 (was 48), rho spacing 1/1024 (was 1/128). Historical iA values above are physical windows; executable tuples scale them by eight.

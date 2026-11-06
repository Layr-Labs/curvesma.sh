# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T is replaced by T' = T^e, e uniform in [1, r), log T = log T'/e
mod r (2/r charged; proof.md section 1). The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6; its LLL basis is
certified exactly (Gram determinant p^12, Lovasz conditions; [11]). The boot
element has a single (g-side) norm. The boot enumerates the radius-4 ball
(2^23.67 points, one of each +- pair), computes each norm and runs the
`ecm16.py` schedule for 193.71-bit primes: a relation has every prime <= 2^215,
<= 16 in (2^69, 2^193.71], <= 1 in (2^193.71, 2^215]. An exhausted ball draws a
fresh exponent and certified lattice (2^36 each). Distinct relations credit
UDUP = 0.65, as at nodes. Lazy boot (below): it finds one distinct relation,
descends it, and finds the next only if that fails, at most 18; expected
2^56.68 trials of 2^67.02 instructions; failure needs all 18 to fail.

## Nodes

A node is a special-q ideal (l bits, side f/g) with lattice L_q (det q,
LLL-reduced). Its region is the radius-R integer ball in the reduced basis,
one point per +- pair (exact half counts in `stats16.b85`; [7] recomputes R <= 8). Two node types:

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

is a lower sum (exact weight int dt/t per mesh cell, rho at its left endpoint; cells crossing sum t = m dropped). Norm
sizes m come from measured offset bins over 61,440 exact ball norms on 160
production special-q lattices (10 sizes 72-300 bits, both sides, 6 radii): cf
sextile maxima paired with the overall cg maximum, shifted by worst sampled
per-(ideal, radius) mean minus pooled mean (cf 20.68 + 1.58 -> +23, cg -2.53 +
7.81 -> +6). Boot norms use 3,072 exact norms on 8 certified boot lattices
(p90 + 4, p99 + 4 bins). **Sub-bin yields (after 3486094d and recover-Q v7).**
For pr the sextiles are subdivided (`QBINS_SUB`: 24 cf-quantile cells, per-cell cg
maxima on a fresh uniform grid, same +23/+6, top edge pinned at 109.6/20.5;
`BOOT_SUB`: 10 R = 4 deciles + 4, edge pinned at 16.7); sub-cell maxima dominate
their cells, so each term stays a lower sum. Charged pr = min(sub-cell mean, 4 x
sextile pr) at nodes, min(sub-cell mean, 1.5 x sextile pr) at the boot. Survivors
include cells x pr (a larger pr charges more cofactorisation); UDUP support is measured against this pr.
mu = cells x P x UDUP x detection x retention (UDUP
= 0.65; ECM detection 1 - 0.001 x expected factors; retention from the next
section; the boot carries only the counted fraction). Survivors per cell use
(T/69) rho((m - T)/69) per side, decreasing in m, evaluated at each bin's
lowest norms (previous sextile's raw cf max, observed cg min -54.56, margins
(23, 6); special-q at the bin's upper edge).

## Factor-base ideals without a logarithm

Under the main-phase occurrence model (`filtsim.py` at the descent input of 2^63.575 relations, [4]), an ideal of norm q = 2^(69 t) on side s occurs
with probability 2 kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/rho(u_s), u_s from
mean norms 739.1/498.2 (G's 739.167/498.217 truncated; [4] checks them; [10y] covers the filtering side).
Occurrences are Poisson(lambda_s(b)); an ideal occurring <= 2x is missing. A
relation of an m-bit side carries on average sum_b kappa(m/69, b/69)
ln((b+0.1)/b) P(Poisson <= 2) missing ideals (1/10-bit bands b < 69):
0.20-1.86 per f side, 0.41-3.72 per g side over u = 7..39 (3.773 at the boot's
2,709-bit norm). Each node takes one option:

- **discard**: relations containing a missing ideal are rejected, multiplying
  mu by exp(-missing); nodes with q <= 2^90 always discard (the code forces
  it for bins up to 85 bits; bin 90 keeps its searched thin flag), so missing
  ideals never recurse there;
- **descend**: each missing ideal Q is a further child, sieved in the lattice of
  Q q' (section 'Recurrence and success') and charged the largest work and failure probability
  of the discarding nodes with q <= 2^90 (Q q' <= 2^89; smaller Q are cheaper);
  the boot descends its missing ideals the same way.

Cascading removal can also lose ideals with >= 3 occurrences (an ideal whose
three rows each hold a different singleton is purged with them). [6b] runs purge
(target 0.0001) and back-substitution on the [6] system (excess 1.08, 2^12, seed 1) and counts
every ideal left without a log, whatever its occurrence count: for lambda in
[1,2), [2,3), [3,5), [5,8), >= 8 the measured fraction is
0.0875/0.0624/0.0248/0.0045/0 against the model's 0.7950/0.5610/0.2773/0.0635/0.0016
(the model counts every ideal with <= 2 occurrences, recovered or not); the 6 of 3,697
lost with >= 3 occurrences (1 + 4 + 1) are inside the measured fractions. Over seeds 1-4 at 2^12 and 1-2 at 2^13 the model's
count is >= 6.3x the measured one in every band with a loss.

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
exclusions); [7] evaluates every sieve node (min 0.8592; largest one-side loss
< 1e-5).

**ECM schedule (`ecm16.py`).** Suyama-form Montgomery curves, cheapest (B1,
B2/B1) grid point with wheel 30030; stage-1 x-only ladder (11 mults/bit),
stage-2 cross products (3/pair), one binary gcd per stage. Per-curve success
P1 is the Dickman lower sum on the curve order; curves = ceil(ln 0.001 /
ln(1 - P1)) at 367.76 (L/5)^2 per mult. [12] prints node/boot schedules (boot:
B1 = 2^28, 737933 curves, 2^67.01/trial) and matches P1 on random 36/40-bit
primes.

## Recurrence and success

Distinct relations per node are negative binomial (mean mu, shape r): r = 0.5
for effective windows of <= 1 prime (base nodes (bin (69, 75)) have empty
windows), r = 1 otherwise. Shapes come from toy-ball fits (empty-window rows
in `stats16.b85`, printed by [11b]): one-prime windows fit 0.60, two/three-prime
1.86-12.2. Of 12 empty-window configs, all with mean >= 7.4 show 0-4 ideals
without a relation vs 5.4-12.4 predicted at r = 0.5; the three with means
3.2-4.0 show 19/11/21 vs 18.1/17.2/16.4, and raising lb clears them (0 vs
7.3-11.2).

A child subtree fails with probability phi; phi(node) = max((1 + mu (1 -
psi)/r)^-r, 0.02), psi = window sum of children's phi; a node tries <=
1/(1 - psi) relations on average. Node work = own sieve/ECM work + 1/(1 - psi)
x children's work. The plan requires mu >= 18, phi <= 0.05, psi <= 0.8
everywhere ([7] checks). The descent keeps a special-q node counter over both
logs and aborts past N_MAX = 2^27.3 nodes (counted failure); E[N] <= 2^15.6813
(`counts16.py`, unrounded 2^15.68129), so P(abort) <= 2^-11.61 (Markov). All N_MAX nodes in the
smallest bin (2^70.3 ideals) with 2^6 for 1/p weighting: shared-ideal
probability <= 2^-10.70, charged in full ([7]).
A factor-base ideal Q without a log
(norm <= 2^69) is never descended in its own fixed lattice: each visit draws a
fresh ideal q' of about 2^20 whose log is known (main-phase factor base) and
sieves the lattice of Q q', so a recurrence of Q samples a different sublattice.
Its node is charged as the worst node with q <= 2^90; nodes up to that size
discard relations touching an ideal without a log.

The shipped plan has 62 nodes (bins (69, 75) plus 5-bit bins 75-225, both
sides) and the boot (n1 = 215, nA = 193.71, kA = 16, R = 4), chosen by
coordinate search over (n, iA, kA, kB, R, discard) minimising the metered cap
subject to the three conditions. The plan and boot are recover-Q v7's (searched under the
lazy model at the sub-bin yields). At this
package's inputs (LOG2_RELS = 63.575, mean norms 739.1/498.2, UDUP = 0.65) with
THIN_MAX 85: min mu 44.7, max phi 0.0268, max psi 0.536 (boot; nodes 0.391). Base
nodes sieve R = 256 (f, mu 932, phi 0.0232) and R = 288 (g, mu 786, phi 0.0253);
all 62 nodes sieve (R = 192..384). One log <= 2^125.714 instructions (boot own
2^123.70 + subtree 2^124.20 x retry 2.153); **two logs <= 2^112.4260 units**.
Meter stops the descent at 2^113.279 units (Markov: E[W]/cap <= 0.5536323; <= 0.6089956 at 10% more expected
work, where success stays >= 0.39). P(boot fails) <= 1.307e-05 ([7]).

## Lazy cofactorisation and lazy boot (9332a74d)

Relations are tried one at a time; the first whose subtree succeeds ends the node (the eager recurrence of 9ff4fc85
cofactored every survivor and collected all 18 boot relations first). Sieve, survivors, windows, phi/psi and the node
counter are the eager recurrence's.

- **Nodes.** After the sieve (charged in full) the s survivors are cofactored in uniformly random order; each distinct
  relation found is descended at once; the node stops at the first relation whose children all succeed (or fails when
  the survivors run out, exactly the eager failure event). Given M distinct relations among the s survivors, the k-th
  first appears at expected position <= k (s + 1)/(M + 1) (a class's first occurrence is at or before a uniformly
  chosen member, and M such members occupy a uniform M-subset). The number J of relations tried is independent of the
  order (H5) with P(J > j) <= psi^j, so E[min(J, M)] <= 1/(1 - psi) and

  E[processed]/s <= E_M min(1, (1 + 1/s)/((1 - psi)(M + 1)) + psi^M),

  M negative binomial (mean mu, shape r) as in phi_nb. `descent16.lazy_fraction` sums this; it stops once m > mu and a
  term pm < 1e-20, then adds the unsummed mass and pm (r + mu)/r (the NB tail bound for r <= 1) at fraction 1. Node own
  work is LLL + sieve + (survivors x per-survivor cost) x that fraction. [7] also simulates the process (`lazy_mc`: NB
  count, uniform positions, independent outcomes) and checks the simulated mean stays below the bound.
- **Boot.** The boot finds one distinct relation (expected per_rel trials, the eager per-relation rate) and descends it
  before looking for the next, at most BOOT_RELS = 18. Trials are i.i.d. and min(J, 18) is a stopping time, so by Wald
  the expected trial work is per_rel x (1 - psi^18)/(1 - psi) x trial (eager: 18 per_rel x trial). Failure is still
  psi^18.

The recurrence W = own + retry x children is otherwise unchanged.

## Supporting measurements ([11])

Norm bins (offline, 20 lattices/size): lattices at 150/260 bits stay in the node bins (cf max 53.4 <= 109.6,
cg max 5.19 <= 20.5), the boot lattice in the boot bins (cg max 6.6 <= 16.7); [11] re-samples per run.

**Unit duplicates (complete orbits).** A unit eps of O_h maps a ball cell v = (a, b) to eps v with the same relation
(sigma and units outside K_h do not act on the lattice), so v's duplicates in the ball are O(v) = {eps : eps v in the
ball} and the distinct yield per cell is E[s/m_true], m_true = |O(v)|/2, s the semismooth probability at v's exact
norms; the charge P x UDUP needs yr = E[s/m_true]/pr >= UDUP. e -> e v is injective Z^6 -> L, so {e : e v in the ball}
is the integer point set of the ellipsoid e M_v M_v^T e^T <= R^2 (M_v the B-coordinates of y^j v); `unit_orbit.py`
enumerates it exactly (rational LDL^T, Fincke-Pohst with integer square roots, no floating point, no box on the unit
coefficients) and keeps the points of norm +-1. Campaign on this plan (STNFS-KIT `udup_campaign.py`, not shipped: [11] reruns the campaign live; 8 fresh certified
boot lattices x 250 uniform R = 4 cells; one fresh special-q lattice per plan entry at the bin's top edge with its R and
window x 100 cells; 8,200 cells), against the CHARGED (sub-bin, capped) pr: boot yr min 0.890 (95% bootstrap lower end
0.842), mean 0.933; nodes min 2.59 (220g; smallest lower end 2.09, 195g), all above 0.65. Orbits reach 44 cells (unit
heights up to 34) at the boot and 62 (heights up to 75) at nodes, so no coefficient box is a safe truncation; the
complete enumeration needs none. [11] repeats it live (24 boot cells, 12 at (220, g); thresholds 0.70 boot / 1.0 node
against the charged pr; this build 0.896 / 2.05) with a sampler check.

## Refined integration representation

Mesh 384, rho spacing 1/1024 (e3940d11): eight sub-cells per old cell, exact log measure, rho at left endpoints (lower
sums); `audit_refinement.py` recomputes four probes. Executable iA tuples are the physical windows times eight.

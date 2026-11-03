# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T is replaced by T' = T^e, e uniform in [1, r), log T = log T'/e
mod r: uniform on nonidentity mu_r (distance 1/r from uniform on all mu_r;
identity z handled directly; success charges 2/r). The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6; its LLL basis is
certified exactly (Gram determinant p^12, Lovasz conditions; [11]). The boot
element has a single (g-side) norm. The boot enumerates the radius-4 ball
(2^23.67 points, one of each +- pair), computes each norm and runs the
`ecm16.py` schedule for 199.79-bit primes: a relation has every prime <= 2^215,
<= 16 in (2^69, 2^199.79], <= 1 in (2^199.79, 2^215]. An exhausted ball draws a
fresh exponent and certified lattice (2^36 each). Distinct relations credit
UDUP = 0.5, as at nodes. It stops at 18 distinct relations: 2^58.39 trials of
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
(p90 + 4, p99 + 4 bins). mu = cells x P x UDUP x detection x retention (UDUP
= 0.5; ECM detection 1 - 0.001 x expected factors; retention from the next
section; the boot carries only the counted fraction). Survivors per cell use
(T/69) rho((m - T)/69) per side, decreasing in m, evaluated at each bin's
lowest norms (previous sextile's raw cf max, observed cg min -54.56, margins
(23, 6); special-q at the bin's upper edge).

## Factor-base ideals without a logarithm

**Relation-count input.** The descent uses LOG2_RELS = 63.999 = the pool estimate
2^64.1694 minus two pool standard errors (2 x 0.0848 bits), floored, so it does not
assume the point estimate; [4] checks LOG2_RELS <= estimate - 2 SE. The full run prints
the sensitivity ([7s]): with this plan the two-log work is 2^115.0487, 2^115.0800,
2^115.1253 at estimate, estimate - 1 SE, estimate - 2 SE (boot failure 7.80e-06,
1.03e-05, 1.59e-05; boot psi 0.520, 0.528, 0.541).

**Occurrence model.** An ideal of norm q = 2^(69 t) on side s occurs in a main-phase
relation with probability 2 kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/rho(u_s), u_s from
mean norms 744.3/498.9 (at or below the relation-weighted means; [4] checks them).
Occurrences are Poisson(lambda); the model counts an ideal as missing with
probability P(Poisson(lambda) <= 2) x c(lambda), and a relation of an m-bit side
carries sum_b kappa(m/69, b/69) ln((b+0.1)/b) P(Poisson <= 2) c(lambda_b) missing
ideals (1/10-bit bands b < 69).

**Calibration c(lambda).** The [6b]-type filtering simulation (`filtsim.py`: signed
build, purge, cascading back-substitution) at excess 1.2 and the production norms
gives the unknown-log fraction among all ideals of a lambda band, zero-occurrence
ideals included, against the P(Poisson <= 2) model. Over 13 runs (2^12: 3 seeds,
2^13: 3, 2^14: 3, 2^15: 2, 2^16: 2) the measured/model ratio is, per band (max over
runs, at norms 743.3/498.9, one bit below this package's): lambda in [1,2) 0.396, [2,3) 0.271, [3,5) 0.184, [5,8)
0.159, never above the 2^12 maxima (e.g. [1,2): 0.396, 0.371, 0.365, 0.354, 0.364 at
2^12..2^16); [7c] reruns 2^12 seed 1 at this package's means 744.4/499.0 in every run
(0.394, 0.239, 0.164, 0.096) and checks each ratio <= the charged maxima.
The calibration is load-bearing and part of H5: the measured/model ratio is assumed
not to rise with scale (per-scale maxima 2^12..2^16 never exceed the 2^12 ones). The model
charges c = 1.5 x (0.396, 0.272, 0.185, 0.160) = (0.594, 0.408, 0.2775, 0.240) and c = 1
for lambda < 1 and lambda >= 8, which the simulation does not cover.

**Missing-log children.** Each node either discards relations touching a missing
ideal (mu x exp(-missing); bins up to (85,90] always discard) or descends each
missing ideal Q (norm <= 2^69) as an ordinary prime special-q node in its own
lattice L_Q: its norms are at most those of the first bin (69, 75] on either side,
so it is charged the larger first-bin work and failure probability (a node's
f-side norm grows as 72 c - l with c = (l + S)/12, i.e. as 5 l + 6 S, increasing in
the special-q size l; on g it is 12 c - l = S, constant in l). No composite
lattice is used. Two missing-log children on the same prime ideal are a repeat:
the largest probability that a child is a given ideal is 2^-63.34 (band weights
kappa 2^-b P(missing) over their total, both sides; [7r] prints it), so among
<= N_MAX nodes the repeat probability is <= 2^(2 x 27.3 - 1 - 63.34) = 2^-9.74; with
the special-q term 2^-10.70 the repeat total is 2^-9.14 ([7]).

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
factor-base mass X exceeds s = 96 - 1 (bound gap < 2e-11 bits; the forms use certified roots, below) - m/128
(rounding: each of the <= m/8 sieved hits of an m-bit side loses < 1/16 bit to
floor(16 log2 p)). X sums independent local contributions (whole valuation at
p < 2^8, excess powers at p >= 2^8), each with one local bound.
Per prime P of the base field with residue field F_q, a base-primitive cell
has v_P >= j with probability at most min(1/q, mu q^-ceil(j/mu)) per root of
multiplicity mu (Hensel; equality q^-j for simple roots), independently over
primes (CRT). Root multiplicities and ramification for p < 2^9 are tabulated
(`retain16.structure`, `stats16.b85`) and recomputed by [11b] at 2, 19, 23 and
a seventh of the other primes; above 2^9 excess powers are bounded per dyadic
block. Means 3.35 (f) and 7.66 (g) bits, P(X > 96) = 1.2e-12 (f) and 2.3e-10
(g) per cell.

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
7.3-11.2). 

A child subtree fails with probability phi; phi(node) = max((1 + mu (1 -
psi)/r)^-r, 0.02), psi = window sum of children's phi; a node tries <=
1/(1 - psi) relations on average. Node work = own sieve/ECM work + 1/(1 - psi)
x children's work. The plan requires mu >= 18, phi <= 0.05, psi <= 0.8
everywhere ([7] checks). The descent keeps a special-q node counter over both
logs and aborts past N_MAX = 2^27.3 nodes (counted failure); E[N] <= 2^16.91
(`counts16.py`), so P(abort) <= 2^-10.39 (Markov). All N_MAX nodes in the
smallest bin (2^70.3 ideals) with 2^6 for 1/p weighting, plus the missing-ideal
repeats above: repeat probability <= 2^-9.14, charged in full ([7]).

The shipped plan has 62 nodes (bins (69, 75) plus 5-bit bins 75-225, both
sides) and the boot (n1 = 215, nA = 199.79, kA = 16, R = 4), chosen by
coordinate search over (n, iA, kA, kB, R, discard) minimising modelled two-log
work subject to the three conditions, at a refined quadrature (mesh 384,
rho spacing 1/1024: each old mesh cell splits into eight with left-endpoint
rho values, so the refined integrals bound the same model).
The plan was searched at LOG2_RELS = 63.895 and mean norms 743.2/498.8 (UDUP =
0.5) and is evaluated here at this package's inputs (63.999, 744.3/498.9): min mu
39.6, max phi < 0.05, max psi 0.542 (boot). One log <= 2^128.413 instructions;
**two logs <= 2^115.1254 units**. Meter stops the descent at 2^115.851 units
(Markov on the unrounded E[W]: <= 0.604723). P(boot fails) <= 1.5924e-05 ([7]).

## Certified linear forms ([11b])

The 6 roots eta of h, the 72 f-side roots xi of F(eta, x) and the 12 g-side roots of
x^2 + |u| eta are enclosed by Weierstrass disks (`descent_witness.certify_forms`:
Durand-Kerner in 60-digit Decimal, radius n |p(z_i)| / prod |z_i - z_j| with the
coefficient errors from the enclosure of eta added; pairwise disjoint disks hold
exactly one root). The largest radius is 1.3e-57 relative to the root. The cell
bound adds 2^-46 x (sum |v_j| |eta|^j + |xi| sum |v_6+j| |eta|^j) to |L(v)|; the
double rounding of the roots (2^-53 relative, at most 6 powers) and of the
evaluation (at most 12 terms, 2^-53 each) contribute < 18 x 2^-53 < 2^-46 of that
sum, so the bound is >= |N| at every cell, with no exceptional cells; [11b]
reruns the enclosure in every run ([7r]) and checks it.

## Supporting measurements ([11])

Norm bins (offline, 20 lattices/size): fresh lattices at 150/260 bits stay in
the node bins (cf max 53.4 <= 109.6, cg max 5.19 <= 20.5), the certified boot
lattice in the boot bins (cg max 6.6 <= 16.7); [11] re-samples a few points per
run. Unit duplicates (offline): mates of 480 ball points under 840 small units
give distinct/raw 0.779; the model charges 0.5 ([11] repeats 12 points: 0.944).

## Refined integration representation

Mesh 384, rho spacing 1/1024 (each cell's left endpoint gives a lower rho value, so the integrals are lower sums of the model).

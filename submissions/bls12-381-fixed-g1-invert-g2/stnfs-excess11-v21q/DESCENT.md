# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] (full run) and [12] (full run) check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T is replaced by T' = T^e, e uniform in [1, r), log T = log T'/e
mod r: uniform on nonidentity mu_r (distance 1/r from uniform on all mu_r;
identity z handled directly; success charges 2/r). The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6; its LLL basis is
certified exactly (Gram determinant p^12, Lovasz conditions; [11]). The boot
element has a single (g-side) norm. The boot enumerates the radius-4 ball
(2^23.67 points, one of each +- pair), computes each norm and runs the
`ecm16.py` schedule for 194.42-bit primes: a relation has every prime <= 2^209,
<= 16 in (2^69, 2^194.42], <= 1 in (2^194.42, 2^209]. An exhausted ball draws a
fresh exponent and certified lattice (2^36 each). Distinct relations credit
UDUP = 0.5, as at nodes. It stops at 18 distinct relations: 2^60.77 trials of
2^67.12 instructions over ~2^37.1 lattices, tried in turn; failure needs all
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

Input LOG2_RELS = 63.62 <= estimate 2^63.6955 - 2 pool SE = 2^63.6243 ([4] checks it). An ideal of norm 2^(69 t) occurs Poisson(lambda) times,
lambda from 2 kappa_s(t)/q (kappa = rho(u - t)/rho(u), norms 739.9/498.9 at or
below the relation-weighted means, [4]); it is charged missing with probability
P(Poisson(lambda) <= 2) c(lambda). Calibration c: unknown-log fraction among all
ideals of a band (zero occurrences included) over the model, `filtsim.py` purge
and back-substitution at excess 1.1, norms 739.9/498.9, 13 runs (2^12-2^14: 3 seeds, 2^15-2^16: 2):

| lambda | [1,2) | [2,3) | [3,5) | [5,8) | other |
|---|---:|---:|---:|---:|---:|
| max ratio, 2^12..2^16 | 0.406 | 0.310 | 0.238 | 0.233 | - |
| charged c = 1.5 x cap | 0.609 | 0.465 | 0.357 | 0.3495 | 1 |

All band maxima occur at 2^12-2^13; [7c] reruns 2^12 in every run. The
calibration is load-bearing (with c = 1 the plan is infeasible at this input; H5).

Each node discards relations touching a missing ideal (mu x exp(-missing); bins up
to (85,90] always) or descends each missing ideal Q (norm <= 2^69) as a prime
special-q node in its own lattice L_Q, charged the larger first-bin (69,75] work
and failure (f-side norm 72c - l = 5l + 6S grows with l; g-side 12c - l = S is
constant). Repeats: a child is a given ideal with probability <= 2^-63.68 ([7r]),
so <= 2^(2 x 27.3 - 1 - 63.68) = 2^-10.08, plus 2^-10.70 for special-q: 2^-9.35.

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
(`retain16.structure`, `stats16.b85`) and recomputed by [11b] (full run) at 2, 19, 23 and
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
B1 = 2^28, 795463 curves, 2^67.12/trial) and matches P1 on random 36/40-bit
primes.

## Recurrence and success

Distinct relations per node are negative binomial (mean mu, shape r): r = 0.5
for effective windows of <= 1 prime (base nodes (bin (69, 75)) have empty
windows), r = 1 otherwise. Shapes come from toy-ball fits (empty-window rows
in `stats16.b85`, printed by [11b] (full run)): one-prime windows fit 0.60, two/three-prime
1.86-12.2. Of 12 empty-window configs, all with mean >= 7.4 show 0-4 ideals
without a relation vs 5.4-12.4 predicted at r = 0.5; the three with means
3.2-4.0 show 19/11/21 vs 18.1/17.2/16.4, and raising lb clears them (0 vs
7.3-11.2). 

A child subtree fails with probability phi; phi(node) = KAPPA_N x max((1 + mu (1 -
psi)/r)^-r, 0.02) with KAPPA_N = 1.2, psi = window sum of children's phi; a node tries <=
1/(1 - psi) relations on average. Node work = own sieve/ECM work + 1/(1 - psi)
x children's work. The plan requires mu >= 18, phi <= 0.05, psi <= 0.8
everywhere ([7] checks). The descent keeps a special-q node counter over both
logs and aborts past N_MAX = 2^27.3 nodes (counted failure); E[N] <= 2^16.88
(`counts16.py`), so P(abort) <= 2^-10.42 (Markov). All N_MAX nodes in the
smallest bin (2^70.3 ideals) with 2^6 for 1/p weighting, plus the missing-ideal
repeats above: repeat probability <= 2^-9.35, charged in full ([7]).

The shipped plan has 62 nodes (bins (69, 75) plus 5-bit bins 75-225, both
sides) and the boot (n1 = 209, nA = 194.42, kA = 16, R = 4), chosen by
coordinate search over (n, iA, kA, kB, R, discard) minimising modelled two-log
work subject to the three conditions, at a refined quadrature (mesh 384,
rho spacing 1/1024: each old mesh cell splits into eight with left-endpoint
rho values, so the refined integrals bound the same model).
The plan was re-searched at the descent inputs (LOG2_RELS = 63.62, mean norms
739.9/498.9, UDUP = 0.5, KAPPA_N, KAPPA_T): min mu 45.2, max phi < 0.05, max node
psi 0.407, boot psi 0.618. Boot variants (24-40 kept relations with psi up to
0.9) were searched and cost more. One log <= 2^130.486 instructions;
**two logs <= 2^117.1982 units**. Meter stops the descent at 2^117.918 units
(Markov: E[W]/cap <= 0.6071577). P(boot fails) <= 1.7171e-04 ([7]).

## Dependence margins (H5; `toy_descent.py`)

The recurrence multiplies independent-looking failure terms. Two factors cover
dependence. Node: phi charged at KAPPA_N = 1.2 x the NB value. Boot: given the
fixed log table each boot relation comes from a fresh uniform exponent, so the
18 kept relations fail i.i.d. with the table's rate p; we charge p <= KAPPA_T
psi_boot with KAPPA_T = 1.27, and P(boot fails) <= (KAPPA_T psi_boot)^18; its
retry factor uses the same KAPPA_T psi.
`toy_descent.py` runs a scaled two-sided descent with the production rules
(main phase, missing = <= 2 occurrences, special-q lattices, recursive
children, boot from fresh smooth integers) over 72 independent tables
(88,962 nodes, 150 boot relations per table): node failure / NB model 1.014
(95% 1.008-1.021; per bit size <= 1.155); tries / (1/(1 - psi)) 0.969. For
the boot, the between-table variance of p in excess of binomial sampling noise
is -0.0002 (95% CI up to 1.23e-4), i.e. a table-to-table SD of p of at most
0.011 = 5% of p = 0.223; KAPPA_T = 1.27 is more than 5 such SDs. (The raw
largest per-table ratio, 1.38, is binomial noise of 150 draws, SD ~0.15.)

## Certified linear forms ([7r])

The roots of h, of F(eta, x) (72) and of x^2 + |u| eta (12) lie in pairwise
disjoint Weierstrass disks (`certify_forms`, 60-digit Decimal; max relative radius
1.3e-57). Root and evaluation rounding stay below 18 x 2^-53 < 2^-46 of the
absolute-value sum added by the cell bound, so the bound is >= |N| at every cell.

## Supporting measurements ([11])

Norm bins (offline, 20 lattices/size): fresh lattices at 150/260 bits stay in
the node bins (cf max 53.4 <= 109.6, cg max 5.19 <= 20.5), the certified boot
lattice in the boot bins (cg max 6.6 <= 16.7); [11] re-samples a few points per
run. Unit duplicates (offline): mates of 480 ball points under 840 small units
give distinct/raw 0.779; the model charges 0.5 ([11] repeats 12 points: 0.944).

## Refined integration representation

Mesh 384, rho spacing 1/1024 (each cell's left endpoint gives a lower rho value, so the integrals are lower sums of the model).

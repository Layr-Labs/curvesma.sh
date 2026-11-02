# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T becomes T' = T^e, e uniform in [1, r) (uniform on nonidentity
mu_r; 2/r charged for both logs). Identity z is direct; g is nonidentity. The
boot lattice {(a,b) in Z^12 : aB + bA = 0 in F_p[y]/h} (det p^6) has a certified
LLL basis ([11]). The boot enumerates its radius-4 ball (2^23.67 points, one per
+- pair), norms each and runs the `ecm16.py` schedule for 199.79-bit primes:
a relation has every prime <= 2^215, <= 15 in (2^69, 2^199.79], <= 1 above.
Exhausted balls draw a fresh e and lattice (2^36 each). Distinct relations
credit UDUP = 0.45 each; it stops at 16: 2^59.22 trials of 2^67.96 instructions
over about 2^35.6 lattices, failing only if all 16 fail.

## Nodes

A node is a special-q ideal of l bits on side f or g with its lattice L_q
(determinant q, LLL-reduced). Its region is the integer ball of radius R in the
reduced basis, one point of each +- pair; the exact half counts are in
`stats16.b85` and recomputed for R <= 8 by [7]. Two node types:

- **sieve**: every p in [2^8, 2^69] and every degree-one root pair is
  enumerated on the ball by `dsieve16` (KERNELS.md); records (owner =
  index >> 40, 40-bit offset, 16 log2 p, one stream per side) go to the
  owners' 2^40-entry score arrays. Each primitive cell gets a certified
  both-norm upper bound from the 72 + 12 linear forms
  (`descent_witness.cell_bound`; 2,400 charged).
  Cells whose residual masses (log2 bound minus score, minus log2 q on the
  special-q side) are at most T = kA nA + kB n + 60 on both sides survive, get exact norms, trial division and ECM.
  Non-primitive cells repeat c / gcd(c) and are skipped.
- **ecm**: ECM trials on every ball point, without a sieve (never cheaper in
  the shipped plan).

A relation has residual primes <= 2^n per side, <= kA in (2^69, 2^nA], <= kB
<= 1 above nA; window primes are the children, each charged at the most
expensive bin at or below its window.

## Yield and survivors

Per cell, the two-window semismooth probability

  P = sum_{jA <= kA, jB <= kB} 1/(jA! jB!) int rho((m - sum t + alpha)/69) prod dt/t

is a lower sum (exact weights, left-endpoint rho, cells crossing sum t = m
dropped). Norm sizes m come from measured offset bins over 61,440 exact ball
norms (cf sextile maxima + overall cg maximum, shifted to +23/+6); boot norms
from 3,072 exact norms on 8 certified lattices. mu = cells x P x UDUP x
detection x retention (UDUP = 0.45; boot carries only the counted fraction).
Survivors/cell use (T/69) rho((m - T)/69) per side at each bin's lowest norms.

## Factor-base ideals without a logarithm

Under the occurrence model (2^63.97 relations at/below the 1-SE yield bound,
R8), an ideal of norm q = 2^(69 t) on side s occurs Poisson(lambda) with
probability 2 kappa_s(t)/q per relation, kappa_s(t) = rho(u_s - t)/rho(u_s)
(u_s from 745.9/498.9); missing when occurring at most twice. A relation
carries sum_b kappa ln((b + 0.1)/b) P(Poisson <= 2) missing ideals (1/10-bit
bands): 0.08-0.77 (f), 0.25-2.25 (g) over u = 7..39. Each node **discards**
(mu x exp(-missing); forced for q <= 2^95) or **descends** each missing Q in
Q q', charged as the worst discarding node with q <= 2^95; the boot descends
likewise. [6b]: purge losses 0.0470-0 per band vs model 0.7826-0.0017; 0 of
3,816 ideals with >= 3 occurrences lost.

## Retention and cofactorisation

**Counted cells.** Only base-primitive cells count (fraction 0.9100,
`retain16.base_primitive`); cells divisible by a degree-one ideal above a base
prime of degree >= 2 with p >= 2^8 are excluded by 0.9523/0.9919 (f/g). The
factor base holds every prime ideal of norm <= 2^69 above p < 2^8 plus the
degree-one ideals above p in [2^8, 2^69].

**Omitted mass.** A counted cell is rejected only when one side's omitted
mass X exceeds s = 60 - 1 - m/128 (bound gap, score rounding). X sums
independent local contributions (full valuation at p < 2^8, excess powers
above). For P = (p, h_i) with residue field F_q and a root of multiplicity mu,
P(v_P >= j and t reduces to rho) <= min(1/q, mu q^-ceil(j/mu)) by Hensel.
Structures for p < 2^9 are tabulated and rechecked ([11b]); above 2^9 excess
powers are bounded per dyadic block. Means 3.35 (f) and 7.66 (g) bits,
P(X > 60) = 1.3e-07 (f), 7.2e-06 (g) per cell.

**Relation-conditioned loss.** loss = sum_{x > s} P(X=x) Pss_hi(m-x) /
sum_{x <= s} P(X=x) Pss_lo(m-x) at the node's largest norms. Retention =
(1 - loss_f - loss_g) x base-primitive x exclusions; [7] evaluates it at every
sieve node (minimum 0.8575).

**ECM schedule (`ecm16.py`).** Suyama Montgomery curves; cheapest (B1, B2/B1)
grid point with wheel 30030; stage 1 x-only ladder, stage 2 baby/giant steps,
binary gcds. P1 is the Dickman lower sum on the curve order; each multiplication
at 367.76 (L/5)^2. [12] prints node/boot schedules and matches P1 on random
36- and 40-bit primes.

## Recurrence and success

Distinct relations at a node are negative binomial (mean mu, shape r = 0.5
for effective windows of at most one prime, else 1), fit on toy fixed balls
([11b]; one-prime windows fit 0.60, larger 1.86-12.2). Empty-window toys match
r = 0.5 (means >= 7.4: 0-4 ideals without a relation vs 5.4-12.4 predicted;
means 3.2-4.0 clear on raising lb). Base nodes run at mu = 722 and 795.

phi(node) = max((1 + mu (1 - psi)/r)^-r, 0.02), psi = window sum of
children's phi; a node tries at most 1/(1 - psi) relations on average. The plan
requires mu >= 16, phi <= 0.05, psi <= 0.8 everywhere ([7]).
Repeated ideals. An explicit node counter over both logarithms aborts past
N_MAX = 2^27.3 nodes; E[N] <= 2^16.965 (`counts16.py`), so P(abort) <= 2^-10.33.
Two nodes share an ideal with probability <= 2^-10.70 (all N_MAX in the
smallest 2^70.3-ideal bin, 2^6 for 1/p weighting), charged in full ([7]).
A logless Q (norm <= 2^69) is descended in Q q' for a fresh known-log q' of
about 2^20, charged as the worst node with q <= 2^95; all nodes up to that size
discard relations touching a logless ideal (no deeper recursion).

The shipped 62-node plan (bins (69, 75) and 5-bit bins 75-225, both sides; boot
n1 = 215, nA = 199.79, kA = 15, R = 4) was coordinate-searched for minimum
modelled two-log work (R9 changed six entries). Here: min mu 60.4, max psi
0.474; base nodes R = 320 (f) and 384 (g); balls R = 192..384. One log at most
2^130.510 instructions; **two logs <= 2^117.2223 units**; meter at 2^117.9395
(Markov <= 0.60828); P(boot fails) <= 0.4732163^16 <= 6.324e-06.

## Supporting measurements ([11])

Fresh 150/260-bit lattices stay in the node bins, the boot lattice in the
boot bins; [11] re-samples every run. Unit duplicates: 79 production
special-q lattices (lq 60-260; R = 96-384 covering all production radii
192-384, plus boot R = 4; 50-300 ball points each) under 840 small units
(bound 3; bound 4 agrees to <= 0.018 on paired lattices, bound 5 lands in
the same range) give distinct/raw 0.574-0.906, mean 0.78; second seeds move
<= 0.061. The model charges UDUP = 0.45, the survey minimum / 1.28.

## Refined integration (inherited)

Mesh 384, rho spacing 1/1024; PLAN16/BOOT16 scale iA by 8. Subdivision only
raises the lower sum; left endpoints keep rho low. Rho coefficients: 200-digit
Decimal in binary64 (numerical evidence, not a certificate).

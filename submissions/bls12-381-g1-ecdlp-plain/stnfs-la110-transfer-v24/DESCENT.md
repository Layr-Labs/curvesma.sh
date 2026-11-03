# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T becomes T' = T^e, e uniform in [1, r) (2/r charged for both
logs); identity z is direct. The boot lattice {(a,b) in Z^12 : aB + bA = 0 in
F_p[y]/h} (det p^6) has a certified LLL basis ([11]). The boot enumerates its radius-4 ball (2^23.67 points, one per
+- pair), norms each and runs the `ecm16.py` schedule for 208.92-bit primes:
a relation has every prime <= 2^215, <= 14 in (2^69, 2^208.92], <= 1 above.
Exhausted balls draw a fresh e and lattice (2^36 each). Distinct relations
credit UDUP = 0.25 each; it stops at 16: 2^58.44 trials of 2^69.34 instructions
over about 2^34.8 lattices, failing only if all 16 fail.

## Nodes

A node is a special-q ideal of l bits on side f or g with its lattice L_q
(determinant q, LLL-reduced). Its region is the integer ball of radius R in the
reduced basis, one point of each +- pair; the exact half counts are in
`stats16.b85` and recomputed for R <= 8 by [7]. Two node types:

- **sieve**: every p in [2^8, 2^69] and every degree-one root pair is
  enumerated on the ball by `dsieve16` (KERNELS.md); records (owner = index >>
  40, 40-bit offset, 16 log2 p, one stream per side) go to the owners'
  2^40-entry score arrays. Each primitive cell gets a certified both-norm upper
  bound from the 72 + 12 linear forms (`descent_witness.cell_bound`; 2,400
  charged). Cells whose residual masses (log2 bound minus score, minus log2 q
  on the special-q side) are at most T = kA nA + kB n + 60 on both sides
  survive (exact norms, trial division, ECM). Non-primitive cells are skipped.
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
from 3,072 exact norms on 8 certified lattices. **Sub-cell integration (R18).**
For P the six sextiles are refined to 24 cf-quantile sub-cells (`QBINS_SUB`: the
sub-cell's cf quantile edge and its own cg maximum, each shifted by the same
+23/+6 as the sextiles; top edge pinned at the sextile top 109.6/20.5) and the
boot bins to ten R = 4 deciles (`BOOT_SUB`, +4, edge pinned at 16.7). P is the
equal-weight sum of P(cf_k, cg_k) over the sub-cells: every cell of sub-cell k
has cf <= cf_k and cg <= cg_k, and P decreases in both norms, so the sum is a
lower bound on the cell average under the same norm heuristic and shifts as the
sextile charge (the sextiles charge every cell of a sextile at its largest cf
and the overall largest cg). The boost over the sextile charge is capped at 4x
per node and 1.5x at the boot (`PRCAP_NODE`, `PRCAP_BOOT`: disclosed ceilings;
12 of 62 nodes hit the cap, the boot charges 1.33x). Survivors, detection,
missing logs, retention and all costs keep the sextile edges. mu = cells x P x
UDUP x detection x retention (UDUP = 0.25; boot carries only the counted
fraction). Survivors/cell use (T/69) rho((m - T)/69) per side at each bin's
lowest norms.

## Factor-base ideals without a logarithm

Under the occurrence model (2^63.79 relations, at/below the estimate less 2
relative pool SEs, R18; was 1 SE), an ideal of norm q = 2^(69 t) on side s occurs Poisson(lambda) with
probability 2 kappa_s(t)/q per relation, kappa_s(t) = rho(u_s - t)/rho(u_s)
(u_s from 745.9/498.9); missing when occurring at most twice. A relation
carries sum_b kappa ln((b + 0.1)/b) P(Poisson <= 2) missing ideals (1/10-bit
bands). Each node **discards** (mu x exp(-missing); forced for q <= 2^95) or
**descends** each missing Q in Q q', charged as the worst discarding node with
q <= 2^95; the boot descends likewise. [6b]: purge losses 0.0470-0 per band vs
model 0.7826-0.0017; 0 of 3,816 ideals with >= 3 occurrences lost.

## Retention and cofactorisation

**Counted cells.** Only base-primitive cells count (fraction 0.9100,
`retain16.base_primitive`); cells divisible by a degree-one ideal above a base
prime of degree >= 2 with p >= 2^8 are excluded by 0.9523/0.9919 (f/g). The
factor base: every prime ideal of norm <= 2^69 above p < 2^8 plus the
degree-one ideals above p in [2^8, 2^69].

**Omitted mass.** A counted cell is rejected only when one side's omitted
mass X exceeds s = 60 - 1 - m/128 (bound gap, score rounding). X sums
independent local contributions (full valuation at p < 2^8, excess powers
above); for P = (p, h_i) with residue field F_q and a root of multiplicity mu,
P(v_P >= j and t reduces to rho) <= min(1/q, mu q^-ceil(j/mu)) (Hensel);
structures for p < 2^9 tabulated and rechecked ([11b]), excess powers above
2^9 bounded per dyadic block. P(X > 60) = 1.3e-07 (f), 7.2e-06 (g) per cell.

**Relation-conditioned loss.** loss = sum_{x > s} P(X=x) Pss_hi(m-x) /
sum_{x <= s} P(X=x) Pss_lo(m-x) at the node's largest norms. Retention =
(1 - loss_f - loss_g) x base-primitive x exclusions; [7] evaluates it at every
sieve node (minimum 0.8575).

**ECM schedule (`ecm16.py`).** Suyama Montgomery curves; cheapest (B1, B2/B1)
grid point with wheel 30030; P1 is the Dickman lower sum on the curve order;
each multiplication 367.76 (L/5)^2. [12] prints the schedules and matches P1 on
random 36- and 40-bit primes.

## Recurrence and success

Distinct relations at a node are negative binomial (mean mu, shape r = 0.5
for effective windows of at most one prime, else 1), fit on toy fixed balls
([11b]; one-prime windows fit 0.60, larger 1.86-12.2; empty-window toys match
r = 0.5).

phi(node) = max((1 + mu (1 - psi)/r)^-r, 0.02), psi = window sum of
children's phi; a node tries at most 1/(1 - psi) relations on average. The plan
requires mu >= 16, phi <= 0.05, psi <= 0.8 everywhere ([7]).
Repeated ideals. An explicit node counter over both logarithms aborts past
N_MAX = 2^27.3 nodes; E[N] <= 2^17.7165 (`counts16.py`), so P(abort) <= 2^-9.58.
Two nodes share an ideal with probability <= 2^-10.70 (all N_MAX in the
smallest 2^70.3-ideal bin, 2^6 for 1/p weighting), charged in full ([7]).
A logless Q (norm <= 2^69) is descended in Q q' for a fresh known-log q' of
about 2^20, charged as the worst node with q <= 2^95; all nodes up to that size
discard relations touching a logless ideal (no deeper recursion). That the
composite lattice L_(Q q') yields relations like a prime special-q node of its
size is a modelling assumption inside H5 (its ball has the same cell count and
certified norm bounds; the relation yield is taken from the prime-node model).

The shipped 62-node plan (bins (69, 75) and 5-bit bins 75-225, both sides; boot
n1 = 215, nA = 208.92, kA = 14, R = 4) was coordinate-searched (singles, pairs,
(R, n, kA) triples and boot moves) at the integrated P and the 2-SE relation
input for the smallest metered cap holding modelled success 0.45, with the
boot's per-relation failure constrained to <= 0.6732 (the base's validated
level); the shipped meter is then set for 0.50 on that plan (R18; 44
entries differ from R17's plan). Here: min mu 31.0, max psi 0.536 (boot 0.6587), max phi
0.0362; balls R = 192..384. One log at most 2^128.456 instructions; **two logs
<= 2^115.1678 units**; meter at 2^116.1806 (Markov <= 0.49557); P(boot fails) <=
0.6587^16 <= 1.256e-03.

## H5 support and sensitivity (R18)

The recurrence treats the subtrees of distinct relations as independent. A
node's relations are distinct cells of one ball, their window primes are
distinct ideals except with probability <= 2^-10.70 (charged), and each child
subtree runs on its own LLL-reduced special-q lattice with fresh coins, so the
subtrees share no randomness below the parent, only the model. That systematic
channel is bounded by worst-norm lower-sum charging at every node, psi <= 0.8
(retry factor <= 5) and the success margin (model 0.50 against the 0.39 claim).
Sensitivity: inflating every subtree failure probability psi by 25% raises the
expected two-log work by 0.70 bits (50%: 2.8 bits) on this plan; the
boot's per-relation failure is 0.6587 (16 relations kept, window of 14
primes above 2^69). Scaling every
psi, the boot failure and the node-counter term together, the model stays at or
above 0.39 up to a uniform inflation of 11.4%. The boot is constrained to the base's validated per-relation failure (<= 0.6732) by narrowing its window to kA = 14 primes in (2^69, 2^nA] (16 relations kept in both plans); the unconstrained optimum (kA = 16 at psi 0.7311, boot failure 6.7e-03) was about 0.015 bits cheaper and its sensitivity twice as large (1.50 bits at +25%).

## Supporting measurements ([11])

Fresh 150/260-bit lattices stay in the node bins, the boot lattice in the boot
bins; [11] re-samples every run. Sub-cells: 4,608 fresh uniform cells (16
special-q lattices at 72-300 bits, both sides, R = 192..384) and 384 boot cells
(4 certified lattices, R = 4) reproduce the tables: every table edge dominates
the fresh cf quantile of its sub-cell (raw edges within 1.4 bits of the table
edges less 23), every sub-cell's raw cg maximum is at least 6.0 bits under its
table value, the worst per-(ideal, R) group mean (24 cells, about 0.8 bits of
sampling noise) sits 18.5 (cf) and 6.3 (cg) bits above the pooled mean against
the +23/+6 shifts, and the raw maxima 87.3/11.1 stay under the shifted tops
109.6/20.5; boot deciles dominated, raw maximum 10.1 against the 16.7 edge.
[11] repeats the domination test live on 96 fresh cells (`subbins`).

**Unit duplicates (relation-weighted, R18).** The charged distinct yield per
cell is P x UDUP, UDUP = 0.25; its support is E[shat/m] >= 0.25 P, with shat the
model's semismooth probability at the cell's exact (cofactor) norms and m = 1 +
the number of units eps (coefficients in [-5, 5]^6: 1,532 units; eps and -eps
counted separately, a conservative multiplicity) whose mate eps (a - b x) lies
in the ball, on uniform cells (`int_ball_point`) of fresh lattices at the
shipped regimes (special-q at each bin's top edge, the largest cofactor norms,
with the bin's R and window; certified boot lattices at R = 4), against the
integrated P charged here, 400 cells per lattice (95% CI lower end in brackets,
then the largest m):

boot lattices: 0.98 [0.94] m<=69; 0.80 [0.75] m<=97; 0.93 [0.89] m<=45. Nodes: 69f (q 2^75, R 256) 5.51 [4.79] m<=47; 69g (q 2^75, R 288) 5.71 [4.89] m<=127; 100f (q 2^105, R 288) 4.83 [4.33] m<=59; 100g (q 2^105, R 320) 4.46 [3.64] m<=41; 150f (q 2^155, R 288) 3.23 [2.82] m<=25; 150g (q 2^155, R 320) 3.04 [2.71] m<=123; 180g (q 2^185, R 320) 2.53 [2.22] m<=45; 185g (q 2^190, R 384) 2.90 [2.55] m<=53; 190g (q 2^195, R 320) 2.90 [2.58] m<=119; 195g (q 2^200, R 384) 2.96 [2.68] m<=77; 200g (q 2^205, R 288) 2.70 [2.40] m<=29; 205g (q 2^210, R 384) 3.17 [2.72] m<=81; 215g (q 2^220, R 224) 2.66 [2.31] m<=67; 220f (q 2^225, R 384) 4.11 [3.66] m<=65; 220g (q 2^225, R 288) 3.06 [2.69] m<=45.

Minimum over node cells 2.22, over boot lattices 0.75, against the
0.25 required. [11] repeats the measurement live on a fresh boot lattice (24
cells) and a fresh 225-bit g-side lattice at the (220, g) regime (12 cells) with
the 840 units of coefficient bound 3 and requires >= 0.3 and >=
1.0 (`udup_relw`). Method as our recover-Q campaign (57598408); the
base's unconditioned d/r (0.779 on 480 points) is superseded.

## Refined integration (inherited)

Mesh 384, rho spacing 1/1024; PLAN16/BOOT16 scale iA by 8. Rho coefficients:
200-digit Decimal in binary64 (numerical evidence, not a certificate).

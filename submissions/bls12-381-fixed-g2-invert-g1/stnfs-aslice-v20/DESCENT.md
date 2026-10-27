# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [8] and [11]
check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T is replaced by T' = T^e with e uniform in [1, r), and log T =
log T' / e mod r, so every boot sees a uniform element of mu_r. The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6; its LLL basis is
certified exactly (Gram determinant p^12, Lovasz conditions; [11]). The boot
element has a single (g-side) norm. The boot enumerates the radius-4 ball of the
lattice (2^23.67 points, one of each +- pair), computes each norm and runs ECM at
level nA = 214.25 bits: a relation has every prime <= 2^235, at most 12 primes in
(2^69, 2^214.25] and at most one in (2^214.25, 2^235]. When a ball is exhausted it
draws a fresh exponent e and a new certified lattice (2^36 instructions each);
every relation found is a valid start for log T. Distinct relations are credited
at UDUP = 0.25 per raw relation, as at the nodes. It stops at 16 distinct
relations: 2^61.18 trials of 2^68.25 instructions over about 2^37.5 lattices. It
tries its relations in turn and fails only if all 16 fail.

## Nodes

A node is a special-q ideal of l bits on side f or g with its lattice L_q
(determinant q, LLL-reduced). Its region is the integer ball of radius R in the
reduced basis, one point of each +- pair; the exact half counts are in
`stats16.json` and recomputed for R <= 8 by [7]. Two node types:

- **sieve**: every prime p in [2^8, 2^69] and every relative-degree-one root pair
  is enumerated on the ball by `dsieve16.c` (KERNELS.md). Records (64-bit owner =
  index >> 40, 40-bit offset, 16 log2 p, one stream per side) go to the
  owners' 2^40-entry score arrays, one per side. Every primitive cell gets a
  certified upper bound on both norms: each of the 72 + 12 linear forms
  a(eta) - xi b(eta) is evaluated in double precision from the exact integer
  coordinates, and the bound is prod (|L~| + e) with the absolute error
  e = 2^-50 (sum |v_j| |eta|^j + |xi| sum |v_6+j| |eta|^j)
  (`descent_witness.cell_bound`; about 1,500 instructions, charged 2,400 per
  cell). Cells whose residual masses (log2 bound minus score) are at most
  T = kA nA + kB n + 96 on both sides survive, get exact norms, trial division
  and ECM. Non-primitive cells repeat the relation of c / gcd(c) and are skipped.
- **ecm**: ECM trials on every ball point, without a sieve (used where cheaper).

A relation has every residual prime <= 2^n per side, at most kA in
(2^69, 2^nA] and at most kB <= 1 in (2^nA, 2^n]. Those window primes are the
children; a child is charged at the most expensive bin at or below its window.

## Yield and survivors

Per cell, the two-window semismooth probability

  P = sum_{jA <= kA, jB <= kB} 1/(jA! jB!) int rho((m - sum t + alpha)/69) prod dt/t

is evaluated as a lower sum: each mesh cell carries its exact weight int dt/t and
rho at the cell's left endpoints (rho is decreasing, so the left endpoint gives
the smallest value), and cells crossing sum t = m are dropped. The norm sizes m
come from the measured offset bins: (cf, cg) = (log2 N_f - 72c, log2 N_g - 12c -
6 log2|u|) over 61,440 exact ball norms on 160 production special-q lattices (10
sizes from 72 to 300 bits, both sides, 6 radii): the cf sextile maxima, each
paired with the overall cg maximum, shifted by the worst sampled
per-(ideal, radius) mean minus the pooled mean (cf 20.68 + 1.58 -> +23,
cg -2.53 + 7.81 -> +6). Boot norms use 3,072 exact
norms on 8 certified boot lattices (p90 + 4 and p99 + 4 bins). The relation count
is mu = cells x P x UDUP x detection, with UDUP = 0.25 distinct relations per raw
relation and ECM detection 1 - 0.001 x (expected ECM factors: (m_f + m_g)/20 +
2(kA + kB) at a node, m/20 + kA + 1 at the boot). Survivors per cell use the
upper bound (T/69) rho((m - T)/69) per side, which decreases in m, so it is
evaluated at each bin's lowest norms: the previous sextile's raw cf maximum
(the observed minimum -140.05 for the first sextile) and the observed cg minimum
-54.56 (stats16.json), both shifted down by the same margins (23, 6), with the
special-q at the bin's upper edge.

## Factor-base ideals without a logarithm

Under the main-phase occurrence model (`filtsim.py`: 2^64.162 relations, an ideal
of norm q = 2^(69 t) on side s occurs with probability 2 kappa_s(t)/q,
kappa_s(t) = rho(u_s - t)/rho(u_s), with u_s from the relation-weighted mean
norms 750.0 and 498.5 bits (rounded down); [4] checks them), the occurrences of an
ideal are Poisson(lambda_s(b)). The model counts an ideal as missing when it
occurs at most twice (singletons and clique pairs). A relation of an m-bit side
carries on average sum_b kappa(m/69, b/69) ln((b + 0.1)/b) P(Poisson(lambda_s(b))
<= 2) missing ideals (1/10-bit bands b < 69), `descent16.missing`: 0.08-0.77 per
f side and 0.25-2.25 per g side over u = 7..39 (2.29 at the boot's 2,709-bit norm). Each node takes one of two options,
chosen by the plan:

- **discard**: relations containing a missing ideal are rejected (logs are known
  after linear algebra), multiplying mu by exp(-missing); nodes with q <= 2^80
  always discard, so missing ideals never recurse;
- **descend**: each missing ideal is a further child, charged the largest work
  and failure probability of the discarding nodes with q <= 2^80.

The boot descends its missing ideals.

Cascading removal can also lose an ideal with three or more occurrences. [6b]
measures losses on the [6] system: full singleton and clique purge, every
purged-matrix column known, then back-substitution through all relations to a
fixpoint. For lambda in [1,2), [2,3), [3,5), [5,8), >= 8 the measured missing
fraction is 0.083, 0.056, 0.020, 0.0036, 0; the model charges 0.78, 0.54, 0.28,
0.060, 0.0018 (9-17x where losses occur). Ideals lost with >= 3 occurrences:
3 of 3,830 (seed 1; at 2^14, 4 of 15,265 and 10 of 15,279).

Retention ([11b], `stats16.json` descent_retention). A relation cell is lost
only when its bound gap plus the bits of primes below 2^8 plus its sieved
prime-power excess (p < 2^16; e >= 2 for larger p has probability < 2^-16) exceed
96 on a side. Every one of the 72 sieve nodes of the
shipped plan was measured on 4,000 primitive cells of an actual special-q lattice:
the per-cell bound is valid at every cell (gaps below 10^-11 bits) and at least
99.775% of cells are kept. The lost cells are those with the richest small-prime
parts, so the kept cells carry slightly less small-prime help; the mean shift is
at most 0.228 bits (f) and 0.039 bits (g). The retention factor, kept fraction x
yield at the shifted mean norm / yield at the mean norm, is at least 0.9395 on
every node (median 0.999); the ten lowest nodes, repeated on two further
lattices each, give 0.9378-1.0. The model multiplies every sieve node's mu by
**RET = 0.9**. ECM nodes and the boot test every point and need no factor.

## Recurrence and success

Distinct relations at a node are negative binomial with mean mu and shape r:
r = 0.5 for windows of at most one prime and r = 1 otherwise. The window size is
the effective one: the A window counts kA only when nA > 69 and the B window
counts kB only when n > nA. The base nodes (70-bit bins) have n = nA = 69, so
their windows are empty and they take r = 0.5. The shapes are set by the
quantity phi depends on, the number of ideals with no relation, measured on toy
fixed balls (fit rule and per-ideal model means as in toy_fixed_balls):
one-prime windows fit 0.60 and two- or three-prime windows 1.86-12.2. For empty
windows (`stats16.json` nb_empty_window, 12 configurations, both sides
2^lb-smooth), every configuration with mean >= 7.4 has 0-4 ideals without a
relation against 5.4-12.4 predicted at r = 0.5. The three configurations with
means 3.2-4.0 have 19, 11 and 21 against 18.1, 17.2 and 16.4, and raising lb on
the same ideal sizes clears them (0 against 7.3-11.2). The base nodes run at
mu = 309 and 560.

A child subtree fails with probability phi; phi(node) = max((1 + mu (1
- psi)/r)^-r, 0.02), psi = sum over the window of the children's phi, and a node
tries on average at most 1/(1 - psi) relations. The work of a node is its own
sieve/ECM work plus 1/(1 - psi) times its children's work. The plan requires
mu >= 16, phi <= 0.05 and psi <= 0.8 at every node ([7] checks all three).
Distinct relations of one node almost surely have distinct window primes; their subtrees
descend disjoint ideals: two subtrees meet only by drawing the same special-q
ideal, which among the >= 2^63 ideals of each bin a tree's nodes do with
negligible probability. Subtree failures are therefore independent (H5).

The shipped plan has 96 nodes (5-bit bins from 70 to 310 bits, both sides),
chosen bin by bin from the bottom by minimum work subject to phi <= 0.045: min
mu 25.8, max phi 0.0448, max psi 0.7999; 54 nodes discard and 42 descend their
missing ideals. The base nodes sieve R = 288 (f, mu 309, phi 0.040) and R = 384
(g, mu 560, phi 0.030). The 72 sieve nodes (bins 70-235 and 240f-255f) use balls
of R = 192..384; the 24 ECM nodes (240g-255g and both sides of 260-305) use
R = 24..48. Boot plus subtree: one log 2^132.79 instructions; **two logs
2^119.499 units**. The meter stops the descent at 2^122.5 units (Markov <=
2^-3.00). P(a boot fails) <= psi^16 = 0.674^16 = 1.83e-3.

## Supporting measurements (`stats16.json`, [8], [11])

- **Norm bins.** Fresh special-q lattices at 150 and 260 bits stay within the node
  bins (cf max 53.4 <= 109.6, cg max 5.19 <= 20.5); the certified boot lattice
  within the boot bins (cg max 6.6 <= 16.7).
- **Unit duplicates.** Mates eps(a - bx) of 480 ball points on 20 production
  lattices under 840 small units: distinct/raw 0.779. The model charges 0.25.
- **Toy fixed balls.** 6 configurations, 224 ideals: mean distinct relations /
  (model mean x 0.25) = 0.92-1.09 ([8] repeats two ideals: 32 vs 32.8).
- **Toy trees.** Complete descents to 2^10 on the exact child ideals, memoised,
  with backtracking: 36/36 succeed (187 nodes, 18 failed nodes recovered); [8b]
  repeats one.

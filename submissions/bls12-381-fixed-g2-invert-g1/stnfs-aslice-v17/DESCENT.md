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
level nA = 226.7 bits: a relation has every prime <= 2^245, at most 12 primes in
(2^69, 2^226.7] and at most one in (2^226.7, 2^245]. When a ball is exhausted it
draws a fresh exponent e and a new certified lattice (2^36 instructions each);
every relation found is a valid start for log T. Distinct relations are credited
at UDUP = 0.25 per raw relation, as at the nodes. It stops at 16 distinct
relations: 2^54.87 trials of 2^70.58 instructions over about 2^31.2 lattices. It
tries its relations in turn and fails only if all 16 fail.

## Nodes

A node is a special-q ideal of l bits on side f or g with its lattice L_q
(determinant q, LLL-reduced). Its region is the integer ball of radius R in the
reduced basis, one point of each +- pair; the exact half counts are in
`stats16.json` and recomputed for R <= 8 by [7]. Two node types:

- **sieve**: every prime p in [2^8, 2^69] and every relative-degree-one root pair
  is enumerated on the ball by `dsieve16.c` (KERNELS.md). Records (owner =
  index >> 30, offset, 16 log2 p with the side in the top bit) go to the owners'
  2^30-entry score arrays, one per side. Each aligned 64-cell run of the ball
  index gets certified upper bounds on both norms (the interval machinery of
  `embedding_bounds.py` in c-space, charged 400 instructions per cell); cells
  whose residual masses (bound minus score) are at most T = kA nA + kB n + 96 on
  both sides survive, get exact norms, trial division and ECM.
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
sizes from 72 to 300 bits, both sides, 6 radii), taken as per-sextile maxima
shifted by the worst sampled per-ideal mean (+13, +5). Boot norms use 3,072 exact
norms on 8 certified boot lattices (p90 + 4 and p99 + 4 bins). The relation count
is mu = cells x P x UDUP x detection, with UDUP = 0.25 distinct relations per raw
relation and ECM detection 1 - 0.001 x (expected ECM factors: (m_f + m_g)/20 +
2(kA + kB) at a node, m/20 + kA + 1 at the boot). Survivors per cell use the
upper bound (T/69) rho((m - T)/69) per side.

## Recurrence and success

Distinct relations at a node are negative binomial with mean mu and shape r:
r = 0.5 for windows of at most one prime and r = 1 otherwise, below every fit to
the toy fixed balls (0.60 for one-prime windows, 1.86-12.2 for two or three;
fitted per ideal against the model mean, the fits reproduce the observed zero
counts). A child subtree fails with probability phi; phi(node) = max((1 + mu (1
- psi)/r)^-r, 0.02), psi = sum over the window of the children's phi, and a node
tries on average at most 1/(1 - psi) relations. The work of a node is its own
sieve/ECM work plus 1/(1 - psi) times its children's work. The plan requires
mu >= 16, phi <= 0.05 and psi <= 0.8 at every node.

The shipped plan has 96 nodes (5-bit bins from 70 to 310 bits, both sides):
min mu 19.4, max phi 0.050, max psi 0.793. Most nodes sieve balls of R = 20..192;
nodes in some bins from 240 to 280 bits use ECM on balls of R = 8 or 12. Boot plus
subtree: one log 2^127.16 instructions; **two logs 2^113.873 units**. The meter
stops the descent at 2^119 units (Markov 2^-5.13). P(a boot fails) <= psi^16 =
0.644^16 = 8.7e-4.

## Supporting measurements (`stats16.json`, [8], [11])

- **Norm bins.** In the sandbox, fresh special-q lattices at 150 and 260 bits
  stay within the node bins (cf max 53.4 <= 99.6, cg max 5.19 <= 19.5) and the
  certified boot lattice within the boot bins (cg max 6.6 <= 16.7).
- **Unit duplicates.** Mates eps(a - bx) of 480 ball points on 20 production
  lattices under 840 small units: distinct/raw 0.779 (71% of points have no
  mate). The model charges 0.25.
- **Toy fixed balls.** Complete enumeration on a toy tower (6 configurations,
  224 ideals): distinct relations / (model mean x 0.25) = 0.92-1.09.
- **Toy trees.** Complete descents to 2^10 that descend the exact child ideal
  (q, y - s, x - t) containing the parent element, memoised by ideal, with
  backtracking: 36/36 succeed (187 nodes, 18 failed nodes recovered by trying
  further relations, at most 16 relations at a node); [8b] repeats one in the
  sandbox.

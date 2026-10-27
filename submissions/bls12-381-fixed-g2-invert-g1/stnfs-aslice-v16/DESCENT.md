# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [8] and [11]
check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T is replaced by T' = T^e with e uniform in [1, r), and log T =
log T' / e mod r, so every boot sees a uniform element of mu_r. The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6; its LLL basis is
certified exactly (Gram determinant p^12, Lovasz conditions; [11]). The boot
element has a single (g-side) norm. The boot takes random points of the radius-14
ball (2^45.13 points, one of each +- pair), computes the norm and runs ECM at
level nA = 244.5 bits: a relation has every prime <= 2^285, at most 16 primes in
(2^69, 2^244.5] and at most one in (2^244.5, 2^285]. It stops at 16 distinct
relations: 2^43.2 trials of 2^72.85 instructions, inside the ball.

## Nodes

A node is a special-q ideal of l bits on side f or g with its lattice L_q
(determinant q, LLL-reduced). Its region is the integer ball of radius R in the
reduced basis, one point of each +- pair; the exact half counts are in
`stats16.json` and recomputed for R <= 8 by [7]. Two node types:

- **sieve**: every prime p in [2^8, 2^69] and every relative-degree-one root pair
  is enumerated on the ball by `dsieve16.c` (KERNELS.md). Records (owner =
  index >> 30, offset, 16 log2 p) go to the owners' 2^30-entry score arrays;
  cells whose residual masses are at most T = kA nA + kB n + 96 on both sides
  survive, get exact norms, trial division and ECM.
- **ecm** (large q): ECM trials on every ball point, without a sieve.

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

Relations at a node are Poisson(mu). A child subtree fails with probability phi;
phi(node) = max(exp(-mu (1 - psi)), 0.02) with psi = sum over the window of the
children's phi, and a node tries on average at most 1/(1 - psi) relations. The
plan requires mu >= 16 and psi <= 0.5 at every node. The work of a node is its
own sieve/ECM work plus 1/(1 - psi) times its children's work.

The shipped plan has 96 nodes (5-bit bins from 70 to 310 bits, both sides):
min mu 16.1, max psi 0.48; nodes below 205 bits sieve balls of R = 80..128,
nodes from 205 bits use ECM on balls of R = 12..16. Boot plus subtree: one log 2^117.77
instructions; **two logs 2^104.486 units**. The meter stops the descent at 2^112
units (Markov 2^-7.51). P(a boot finds no relation) = exp(-16 (1 - 0.34)) =
2.6e-5.

## Supporting measurements (`stats16.json`, [8], [11])

- **Norm bins.** In the sandbox, fresh special-q lattices at 150 and 260 bits
  stay within the node bins (cf max 53.4 <= 99.6, cg max 5.19 <= 19.5) and the
  certified boot lattice within the boot bins (cg max 6.6 <= 16.7).
- **Unit duplicates.** Mates eps(a - bx) of 480 ball points on 20 production
  lattices under 840 small units: distinct/raw 0.779 (71% of points have no
  mate). The model charges 0.25.
- **Toy fixed balls.** Complete enumeration on a toy tower (6 configurations,
  224 ideals): distinct relations / (model mean x 0.25) = 0.92-1.09.
- **Toy trees.** 36/36 complete descents to 2^10 with memoised ideals and
  backtracking (184 nodes, 15 failed nodes recovered by trying further
  relations, max 51 relations at a node); [8b] repeats one in the sandbox.

# Individual logarithms (descent)

Each logarithm: a boot step writes a randomised target as a relation whose cofactors are special-q ideals of at most
255 bits, then special-q nodes descend each ideal to the factor base (B = 2^70 per side). All numbers below are the
default audit run (`audit_cm.py`, mesh 32/185); `cm_models.eval_descent_live` recomputes the whole plan.

## 1. Search set and sieve of a node

For a special-q ideal Q (degree 1 over O_H, norm q) let (a_1, b_1), (a_2, b_2) be an O_H-basis of its rank-2 lattice,
reduced by `lll_16d_int` (exact integer LLL, delta 99/100; the same routine in the audit and in every experiment). A cell
is x = (alpha, beta) in Z^16 (power-basis coordinates of alpha, beta in O_H), element (alpha a_1 + beta a_2,
alpha b_1 + beta b_2). The search set is the ball

    S(R) = {x in Z^16 : 0 < |x|^2 <= R^2, first nonzero coordinate > 0},  cells(R) = (|B_16(R)| - 1)/2,

|B_16(R)| exact (theta series; the audit squares theta_1^8 once more). A unit e of O_H maps x to e x (same relation);
when e x also lies in S(R) the relation is counted again. Distinct relations per cell: udup = 0.0557 charged =
1/zeta_K(2) x E[1/m | relation], m = #{units e : e x in S(R)}/2 computed exactly (Fincke-Pohst on the 8-dimensional
form |e alpha|^2 + |e beta|^2 <= R^2, then N(e) = +-1). Measured on relations of the toy tower ([B], 16-dimensional
ball, 1,014 relations): E[1/m] = 0.057; uniform production cells give 0.45 (0.23 in the lowest-norm decile), so the
relation value is the charge.

Sieve (dsieve16). For every prime ideal P of degree one (norm p <= 2^70) on each side, the hits are the points of the
index-p sublattice L_P = {x : c_P . x = 0 mod p} of Z^16 in the ball. Per (P): HNF basis from c_P, exact integer LLL
(stage bounds of KERNELS.md section 5 at 70-bit entries: <= 2^35.3 instructions = PAIR_SETUP), then Schnorr-Euchner
enumeration of L_P cap B_16(R); each hit adds round(16 log2 p) to the cell's 16-bit score at rank(x) (16 table
lookups). Every production prime has >= cells/p >= 2^16 hits per node, where the enumeration visits <= 4.5 nodes per
hit (`dsieve_toy.py`: identical hit sets to brute force on 16-dimensional balls; 3-4.4 nodes per hit for 10^4-10^5
hits). Per hit: <= 8 nodes x 60 integer instructions (fixed-point centre update, bound test, zig-zag) + 40 (exact
|x|^2 check) + 80 (rank) + 84 (record transfer, as c36p) + 9 (apply) = 693 <= SIEVE_HIT = 1,100. Per cell: certified
norm bounds from incrementally updated embeddings (72 + 24 factors x ~20 instructions) <= BOUND_CELL = 2,400; scan 10,
zero 2. Cells whose bound minus score is within the window allowance on both sides survive to exact norms, trial
division and ECM.

## 2. Recurrence

A node of bin (lo, hi] on side f or g sieves S(R) with the plan's (n, iA, kA, kB): relations whose cofactors are
B-smooth except at most kA ideals in (B, 2^nA] and kB in (2^nA, 2^n]. All survivor work (norms, trial division,
ECM, 100-round Miller-Rabin) is charged on every survivor. mu is the expected number of distinct usable relations
(H1 with de Bruijn rho, lower sums on the mesh, retention, ECM detection, primitive density).
The two-sided smoothness probability of a cell is evaluated on the sampled norm offsets of production lattices
(cm_data `qbins`, `qbins_sub`): the 6 main bins pair each f sextile's upper edge with the largest g offset, and the 24
sub-bins pair f quantiles in ascending order with g quantiles in descending order. The latter is the least
favourable coupling of the two marginals (rearrangement inequality: both smoothness factors decrease in the norm), so
mu is a lower estimate for every joint law of (cf, cg) with these marginals; `prod_sq24_live` checks the marginals.

The number M of distinct usable relations is negative binomial with mean mu and shape r = 0.5 (at most one
non-B-smooth window) or r = 1 (otherwise). A child subtree fails with probability phi;

    phi(node) = max((1 + mu (1 - psi)/r)^-r, 0.02),   psi = sum over the window children of their phi,

and a node tries at most 1/(1 - psi) relations on average: node work = own + 1/(1 - psi) x children's work. The plan
requires mu >= 18, phi <= 0.05 and psi <= 0.8 everywhere; at the shipped plan min mu = 31.2, max phi = 0.0380,
max psi = 0.4399.

Boot. A trial draws a fresh uniform exponent e, maps the target to F_(p^6) and reduces one 16-dimensional lattice
(LLL_BOOT = 2^48 instructions) and takes a - bX from the first reduced vector; it is usable when |N(a^3 - u omega^2 b^3)|
<= 2^2560 (charged probability 0.60) and the cofactor splits into the (255, 195, 15) windows. On 128 fresh targets
(SHA-256 exponents, `boot_fresh` in pool_data.b85) all 128 pass (max 2^2550.3, Clopper-Pearson 95% lower bound 0.9769); the audit
reruns target 0 from scratch. Trials are independent (fresh e). After a usable trial the children are descended; if one fails
(probability <= psi_boot = 0.4081), a new trial starts. Expected work of one logarithm = (boot own +
subtree)/(1 - psi_boot).

Expected two-logarithm work <= 2^110.8072 units (2^110.8071); expected nodes for both logs 2^13.34.
The run meter (proof.md section 6) caps the whole run.

Node counter and repeats. The run aborts (failure) at node 2^22.5: P <= E[N]/2^22.5 = 1.75e-03. A repeated ideal
(equal to an earlier node's ideal) also aborts the run. Conditional on any history, the visited set V (|V| <= 2^22.5,
every member of norm >= 2^lx, lx = 70) is fixed and the current node's cells depend only on its own ideal; the cells
divisible by a given v in V form a sublattice of index N(v). Under H1 a relation of the node is divisible by v with
probability at most (1/N(v)) x ss(m - log2 N(v))/ss(m) <= (1/N(v)) x ss(m - lx)/ss(m), ss the node's semismooth
probability at its largest norm m on v's side. The audit evaluates this ratio at every node: max 206. Hence one tried
relation meets V with probability <= 2 x 206 x 2^22.5/2^70 given the history, and over the expected <= E[N]/(1 - max psi)
relations tried, P(repeat) <= 3.84e-08 (expected number of collisions). Both terms are charged in the success bound.

## 3. Conditioning (measured)

H5 asserts that subtrees of distinct ideals fail independently with the node's phi and that the order of relations
is independent of their outcomes. The relations of a node depend only on its ideal; selection by the parent can
change only the norm law of the child's lattice and the failure law. Records: cond_q15.b85 (lzma, base85, JSON); `cond_q15.verify`
recomputes every figure below from them and replays sample records from scratch.

[A] Production (f_9, g_3), O_H = Z[zeta_15]. 600 children (every degree-1 prime ideal in [2^20, 2^24) dividing
random cells of the 24 production special-q lattices of prod_sq24, 300 per side) against 600 sampled ideals of the same
side and bit length. Mean norm offsets over 16 fixed points per lattice differ by cf -0.41 +- 0.41, cg -0.59 +-
0.34 bits (KS p 0.76 / 0.25); children are not harder than sampled ideals.

[B] Same field, toy pair f_3 = X^3 - (2 + w), g_3 = X^3 + 2 w^2 with a real factor base (all G_12-orbit ideals of
norm <= 2^14, 2-core of 7026 relations), balls of 36 x 128 cells. 96 parent special-q ideals in [2^16, 2^17);
889 encountered children in [2^15, 2^17), each against a sampled ideal of the same side and bit length. Each child is
the single large prime of a distinct relation of its parent, so the children of one parent are the alternative relations
a node tries in turn, and "failure after a sibling failed" is the event that a node's next relation fails after an
earlier one failed. Failure = no usable relation; relations are counted once per torsion-unit class (4,230 of the 4,608
cells are distinct; the model mean uses the same count). The model is negative binomial r = 1 at the independent
model mean (sextile norms, rho, factor-base retention; nothing fitted to the counts):

| set | mean M | model mean | failures (model) | failure after a sibling failed (model) |
|---|---:|---:|---:|---:|
| children, whole ball | 12.82 | 11.09 | 6/889 (89.0) | 0/65 (7.0) |
| sampled, whole ball | 12.48 | 10.69 | 3/889 (92.1) | 0/27 (2.9) |
| children, first 9 alpha rows | 5.69 | 5.42 | 49/889 (159.4) | 32/460 (84.1) |
| sampled, first 9 alpha rows | 5.66 | 5.21 | 43/889 (164.7) | 14/390 (72.0) |

Siblings are children of the same parent. In every row the observed failures, and the failures among siblings of a
failed child, stay below the model's sum over the same children.

Bridge to production. Under H1 a node's relation count is a function of the norm law of its cells and of the
factor base; the production mu is computed from that law measured on production lattices ([C]), and [A] shows that
being encountered leaves the law unchanged at production size. What H5 adds is the shape of the count distribution and
the dependence between a node's alternative relations. [B] tests both where failures are observable, with the same
recurrence ingredients (the independent model mean, NB r = 1): observed failures, and failures after a sibling failed,
stay at 0.07-0.6 of the model in every row, and children behave like sampled ideals. The production recurrence keeps
both margins and adds the floor phi >= 0.02; the run meter bounds the work of any run.

[C] Production lattices (`prod_sq24_live`): 24 special-q lattices of 73-258 bits on both sides reduced live, 1,536
norms; every sextile lies below the charged offsets with >= 23 bits of headroom on g and >= 37 on f.

## 4. Supporting checks

The audit recomputes the theta counts, the primitive density, the plan and its constraints, the boot norm bound on a
fresh SHA-256 target, and the sums of section 2.

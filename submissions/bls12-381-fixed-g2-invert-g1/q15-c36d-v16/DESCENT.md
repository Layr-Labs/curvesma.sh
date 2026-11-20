# Individual logarithms (descent)

Each logarithm: a boot step writes a randomised target as a relation whose cofactors are special-q ideals of at most
255 bits, then special-q nodes descend each ideal to the factor base (B = 2^70 per side). All numbers below are the
default audit run (`audit_cm.py`, mesh 32/185); `cm_models.eval_descent_live` recomputes the whole plan.

## 1. Search set of a node

For a special-q ideal Q (degree 1 over O_H, norm q) let L_Q be its rank-2 O_H-lattice in Z^16 (det q), reduced by
`lll_16d_int` (exact integer LLL, delta 99/100; the same routine in the audit and in every experiment below). With
B_8(R) = {v in Z^8 : |v|^2 <= R^2} and F_Vor the strict Voronoi fundamental domain of O_H^* (with mu_30), the node
enumerates the first N_alpha(R) = max(1, floor((|B_8(R)| - 1) x 0.008700)) canonical alpha in (B_8(R) \ {0}) cap
F_Vor/mu_30 (lexicographic) and pairs each with every beta in B_8(R):

    S(R) = A(R) x B_8(R),  cells(R) = N_alpha(R) |B_8(R)|.

|B_8(R)| is exact (theta series, audit); f_vor/30 = 0.008937 > 0.008700 makes the alpha list available (the audit
also counts R = 2 exactly). Two cells of S(R) are never unit multiples of each other (alpha in F_Vor/mu_30);
non-primitive cells are rejected (density 1/zeta_K(2) = 0.97856 >= 0.95 charged).

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
requires mu >= 18, phi <= 0.05 and psi <= 0.8 everywhere; at the shipped plan min mu = 41.2, max phi = 0.0268,
max psi = 0.2800.

Boot. A trial draws a fresh uniform exponent e, maps the target to F_(p^6) and reduces one 16-dimensional lattice
(LLL_BOOT = 2^48 instructions) and takes a - bX from the first reduced vector; it is usable when |N(a^3 - u omega^2 b^3)|
<= 2^2560 (charged probability 0.60) and the cofactor splits into the (255, 195, 15) windows. On 128 fresh targets
(SHA-256 exponents, cm_data `boot_fresh`) all 128 pass (max 2^2550.3, Clopper-Pearson 95% lower bound 0.9769); the audit
reruns target 0 from scratch. Trials are independent (fresh e). After a usable trial the children are descended; if one fails
(probability <= psi_boot = 0.3200), a new trial starts. Expected work of one logarithm = (boot own +
subtree)/(1 - psi_boot).

Expected two-logarithm work <= 2^107.8098 units (2^107.8097); expected nodes for both logs 2^12.83.
The run meter (proof.md section 6) caps the whole run.

Node counter and repeats. The run aborts (failure) at node 2^22.5: P <= E[N]/2^22.5 = 1.23e-03. A repeated ideal
(one equal to an earlier node's ideal) also aborts the run. Under H1 a child of size 2^x is a prime ideal dividing the
norm of a relation chosen among the node's cells, so a given ideal I of the bin is that child with probability at most
(1/N(I)) / sum over the bin's degree-1 ideals J of 1/N(J); over the plan's bins (6 bits wide at 2^70, 5 bits above) this is at most
15.2 times the uniform 1/#ideals, charged as 2^6. With at most 2^22.5 nodes and #ideals >= 2^lx/(lx ln 2)/12 (lx = 70, orbits of
12) ideals,
P(repeat) <= 2^6 x (2^22.5)^2 / (2 #ideals) <= 5.55e-04. Both terms are charged in the success bound.

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

[C] Production lattices (`prod_sq24_live`): 24 special-q lattices of 73-258 bits on both sides reduced live, 1,536
norms; every sextile lies below the charged offsets with >= 23 bits of headroom on g and >= 37 on f.

## 4. Supporting checks

The audit recomputes the theta counts, the primitive density, the plan and its constraints, the boot norm bound on a
fresh SHA-256 target, and the sums of section 2.

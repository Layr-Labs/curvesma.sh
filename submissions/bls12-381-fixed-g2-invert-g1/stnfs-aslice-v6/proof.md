# STNFS in F_{p^12} with a unit-orbit-correct sieve, priced in Fp12-multiplication units

Track: `bls12-381-fixed-g2-invert-g1` (fixed-argument pairing inversion, FAPI-2). Unit: one dense Fp12
multiplication = C = 20,000 scalar AArch64 instructions (`tracks/classical-group-work-v2.json`,
`tracks/docs/COST_MODEL.md`).

## Claim

For every production instance (uniform nonidentity Q in G2, uniform z in GT including z = 1), the
algorithm below returns P in G1 with e(P, Q) = z (P = infinity exactly when z = 1) with success
probability at least **0.39** over the uniform challenge and the algorithm's coins (section 9 gives
>= 0.98 under the stated heuristics), using total classical computation

  **W <= 2^119.468045 target-group-operation equivalents.**

W charges relation collection (every sieved point, every setup, every cofactorisation), filtering,
linear algebra, Schirokauer maps, factor-base and canonical-list construction, both individual-log
descents at their metered cap, verification and randomness. There is no preprocessing or advice beyond this
charged work; relation collection and linear algebra depend only on the public field and are charged in
full for the single instance. Peak memory is reported in section 13 (not scored).

Score: **119.468045** (= 119468045/1000000).

This is an analytic cost claim supported by exact-norm simulation on a held-out pool, retired-instruction
counts of correctness-checked kernels shipped here (`KERNELS.md`), a filtering simulation with the production
ideal-occurrence profile, and an explicit descent dynamic programme. `evidence.py` recomputes every number of
the ledger from the embedded data and prints `total_log2W_full_ledger`.

### What changed since stnfs-archd-v1 .. v4 (all rejected)

Those packages shared one uncharged error that none of the reviews named: **the unit group of the tower
field acts freely on the candidate space, and all points of one orbit give the same relation.** For this
tower field a region point had on average 8.0 unit multiples inside the v1 region (measured, section
4.2), so the credited yield was over-counted by ~11x (3.5 bits, yield-weighted). This revision sieves a
fundamental domain of the unit action with an architecture built around it. The sieve architecture (per-a
6-dimensional slices, no global walk, no routing), the filtering support (a simulation with the production
occurrence profile, calibrated against a record), the linear-algebra formulation (standard square system
after filtering, Schirokauer maps as right-hand sides), the descent (a dynamic programme over size levels
with explicit per-node costs, both sides, ending in the factor base) are all new; the tower polynomial is
v1-v4's, and section 3 now proves that it generates the target field. The claim is ~5.0 bits above
v1's rejected 114.49: that number was wrong, not conservative.

## 1. Reduction (admissible per rules.md)

Let g = e(P0, Q) with the pinned pairing convention (py_ecc `pairing(Q, P0)`). Q is a nonidentity point of
the prime-order group G2 and e is nondegenerate, so g generates GT = mu_r < F_{p^12}^* (r prime, 255 bits;
embedding degree 12). The challenge satisfies z = g^a for a unique a mod r; returning P = [a]P0 answers it
since e([a]P0, Q) = g^a = z, and z = 1 gives a = 0 and P = infinity. The rules explicitly admit this route.
The problem is therefore one prime-order DLP in F_{p^12}^*, solved as two individual logarithms (z and g)
over one relation-collection/linear-algebra phase: a = log z / log g mod r.

## 2. Baselines

GT is a prime-order group whose only cheap automorphisms are the 12 Frobenius powers (p has order 12 mod
r), so orbit-12 Pollard rho needs sqrt(pi r / 24) = 2^125.96 group operations; the best public claim on
this track is 125.73 (in review). No baseline is qualified on this track.

## 3. The tower, the polynomials and the automorphism

Special Joux-Pierrot tower for the BLS12-381 seed u = -0xd201000000010000 (|u| = 2^63.58):

- K_h = Q[y]/h with **h = y^6 + y^3 + y^2 - y - 1** (the v1-v4 tower). O_h = Z[y] (the polynomial
  discriminant 58673 equals the field discriminant); signature (2, 2), unit rank 3, torsion +-1, regulator
  0.8578.
- **f = y^6 * 3p(x^2/y) mod h = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3)**, degree 12 in x
  over O_h (p(X) the BLS12 characteristic polynomial, 3p(X) = X^6 - 2X^5 + 2X^3 + X + 1; y^6 = 1 + y - y^2 - y^3
  mod h);
- **g = x^2 + |u| y = x^2 - u y**. At the common root x^2 = u y, so f = y^6 3p(u) = 0 mod p.

**The tower is F_{p^12}** (part [0b] verifies both facts explicitly in F_p[y]/h): h is irreducible mod p
(y^(p^6) = y and y^(p^k) != y for k = 1, 2, 3), and u y is a non-square in F_{p^6}: its norm to F_p is
u^6 N(y) = u^6 h(0) = -u^6, a non-square in F_p because p = 3 mod 4, and (u y)^((p^6-1)/2) = N(u y)^((p-1)/2) =
-1 (computed). So g is irreducible over F_{p^6} and F_p[y, x]/(h, g) = F_{p^12}; the relation
Res_y(h, f(x^2 -> u y)) = 0 mod p (part [0], exact) then says f and g share a root that generates F_{p^12}.
(An earlier revision of this package used the Guillevic-Singh polynomial y^6 - y^2 + 1 in a reciprocal
coordinate, for which h(0) = +1, N(y) = 1 and u y is a square: that tower splits and does not represent the
target field, as its review correctly observed. The condition is h(0) = -1 together with p = 3 mod 4.)
alpha values (log2; Guillevic-Singh `alpha_tnfs_2d`, 2000 samples, recomputed for this package):
alpha_f = 4.868 (v1 used 4.721 from 1000 samples), alpha_g = 0.844.

**Automorphism.** x -> -x (f and g are even in x) is the Frobenius^6 conjugation: it maps the relation of
(a, b) to that of (a, -b) and log sigma(P) = -log P mod r, so candidates and unknowns are both quotiented by
sigma (aut = 2). Together with -1 the relation orbits are {(+-a, +-b)}; section 4 states exactly how the
sieve realises this quotient (it is not double-counted with the unit quotient below: -1 is a unit, sigma
is not).

## 4. Relation collection

Candidates are pairs (a, b) of degree < 6 integer polynomials (12 coefficients). A candidate is a relation
when N_f = |Res_y(h, Res_x(f, a - bx))| and N_g = |Res_y(h, Res_x(g, a - bx))| are both B-smooth: exact
smoothness, no large primes, columns = the prime ideals of norm <= B of both sides.

### 4.1 Region

The sieve region is the ellipsoid E = {(a, b) : Q(a, b) <= 1} of volume fraction 1/32 of the box
[-A, A]^12, with **Q block-diagonal** (no a/b cross terms), so E is invariant under a -> -a and b -> -b
separately (the sigma invariance the v1 review demanded; part [3] prints the cross terms: 0). Q is the
yield-weighted covariance ellipsoid (weights = smoothness probability of the exact norms at B = 2^68.5)
fitted on a pool of 300,000 box points disjoint from the held-out evaluation pool, with the a/b cross
blocks zeroed, and rescaled to volume fraction 1/32 (the embedded matrix; part [3] recomputes the
fraction). Only the part Q_a(a) <= theta of the a-ellipsoid is used (4.4); a refit restricted to the
canonical domain was tried and measured no better.

### 4.2 Unit orbits and the fundamental domain F (the central correction)

For a unit eps of O_h, (eps a, eps b) is again a candidate and eps(a - bx) generates the same principal
ideal as a - bx in O_f and O_g (N_{K_h}(eps) = +-1), so both give **the same row** of the relation matrix.
The Schirokauer-map entries differ by Lambda(eps), but eps lies in F_{p^6}, whose elements have r-logarithm
0 (r | Phi_12(p), r does not divide p^6 - 1), so the extra equation Lambda(eps) . vlog = 0 is one of the
(at most unit-rank many) trivial unit relations and carries no new information. Unit multiples are
therefore duplicate relations, exactly as the K_h-unit duplicates of De Micheli-Gaudry-Pierrot
(ASIACRYPT 2021, section 3.4). This has to be charged: on the v1-v4 ellipsoid, a point of E has on average
**8.0 unit multiples inside E** (units of height <= 4; heights 3, 4 and 6 give identical results), and the
yield-weighted multiplicity is worse (central points, which carry the yield, have the most multiples): the
distinct yield of E is **0.084 of the naively credited yield** (measured on the v2 pool; part [2]
re-derives the flags on this package's pool).

The sieve therefore covers a **fundamental domain**: a is *canonical* when Q_a(a) <= Q_a(eps a) for every
unit eps (equivalently, a minimises the a-block form over its orbit); ties have measure zero. Every orbit
of the group generated by the units and sigma has exactly one representative with a canonical, a5 > 0 and
b5 > 0, and

  **F = {(a, b) in E : a canonical, a5 > 0, b5 > 0, Q_a(a) <= theta}.**

The unit set used for the test is every unit of height <= 4: part [U] generates them as products of the
three fundamental units (PARI `bnfinit`, embedded), verifies each exactly (norm +-1 by the exact resultant)
and reproduces the brute-force count 595 (mod sign, excluding 1; `canon.c`). Enlarging the height bound to 6
changes no flag on the whole pool, and no unit outside this set can beat a height-<= 4 minimiser on the
pool's coordinates (a unit with larger coefficients has an embedding far from 1 and multiplies Q_a of any
point of the box by more than the minimiser could gain). The fraction of canonical a (f_c) and the yield
captured by F are measured on the held-out pool, not assumed.

### 4.3 The slice architecture (ARCH-E)

The whole sieve is organised around the canonical a:

1. **Canonical list.** Enumerate the a-ellipsoid {Q_a(a) <= theta} and keep the canonical a with a5 > 0 by
   the unit test (|E_a| x 595 units x ~60 instructions; charged 2^60 units, section 10).
2. **One job per canonical a** ("a-slice"): the 6-dimensional half-ellipsoid
   S(a) = {b : Q_b(b) <= 1 - Q_a(a), b5 >= 1} with one bucket byte per point, addressed by the fixed-cube
   formula over the bounding box (injective). The job runs on one node and touches no other node's memory.
3. **Every sigma-pair of prime ideals of norm p <= X2 is applied to the slice** by exact coset
   enumeration: all ideals above the O_h-prime (p, y - s) share the lattice L(p, s) = {b : b(s) = 0 mod p}
   (index p in Z^6, reduced once per ideal in the Q_b metric: factor-base construction); the ideal and the
   slice fix the coset c e_0 + L(p, s), c = a(s) r^{-1} mod p (r the ideal's root), whose points inside the
   full ellipsoid are enumerated by a Schnorr-Euchner walk carrying the bucket address and b5, folding
   b -> -b when b5 < 0 (`KERNELS.md`). The fold supplies the conjugate ideal's hits, so one walk per pair
   over the full ellipsoid equals both ideals over the half: this is how the sigma quotient is realised
   (together with -1 inside the unit quotient). sigma-fixed ideals (r = 0) add log p / 2. Primes below 2^8
   are not sieved (their logs sit inside the threshold slack).
4. **Ideals with p > X2 are not sieved.** The threshold is set with slack S_f = k_f log2 B on the f side
   and S_g = k_g log2 B on the g side: a point is a survivor when its sieved log-sum reaches
   log2 N - S on both sides. Updates of ideals with p >= 2^20 are also stored as bucket records
   (position, prime index; charged 5 instructions per such update for the store and the apply pass), so
   a survivor recovers its large sieved primes from its bucket; the primes below 2^20 are removed by
   trial division of the exact norm (2^22 instructions per survivor). What remains is the unsieved part,
   at most S_f (S_g) bits by the threshold; it is factored completely by ECM with the full curve budget
   for primes up to B (99% success per factor, no early-abort credit) and the survivor is a relation iff
   every prime found is <= B. A true relation is kept iff it has at most k_f (k_g) primes above X2 on
   the f (g) side; the model charges this loss (keep factor) with the Dickman-conditioned Poisson count,
   and charges the survivors (true and false) at the ECM cost of their cofactor size.
5. **Scan/report**: one byte compare per point of the slice.

There is no 12-dimensional walk, no routing and no shared memory: a job needs its slice, the streamed
per-pair data of the ideals <= X2 (charged per word), and nothing else.

### 4.4 Cost per point of F

With h(p) = |S(a)| / p the expected points per coset, the measured per-update cost depends only on h
(`KERNELS.md`): 16-18 instructions in the long-run regime (h >= 2^16), rising to 146 at h = 16 and 300 at
h = 4; for h < 1 an exact enumeration of the (almost always empty) coset costs 620 instructions per
(slice, pair). The per-(slice, pair) setup is 326 (254 measured + 72 streamed words). The model integrates
NBAR/p x cost(h(p)) over the primes 2^8 < p <= X2 (NBAR = 1.947 ideals per prime, both sides), adds the
bucket-record cost for p >= 2^20, the setups for all pairs <= X2 divided by the slice size (the harmonic
mean over the pool's canonical a, i.e. small slices are charged their true share), the scan, and the
cofactorisation term (survivor probability rho((log2 N_f - S_f)/log2 X2) rho((log2 N_g - S_g)/log2 X2)
times exact norms + trial division + two full ECM runs on cofactors of S_f and S_g bits). At the operating
point this is **186.6 instructions per point of F** (part [4] prints the breakdown).

### 4.5 Operating point (part [4]; the offline search is in `evidence.py`'s model with a wider grid)

- A = 2^10.40, B = 2^69.50, theta = 0.5, X2 = 2^66, (k_f, k_g) = (3, 3), excess 1.3.
- |F| = 2^126.120 points; mean log2 norms N_f = 760.0, N_g = 502.2 (u_f = 10.93, u_g = 7.23).
- Relation orbits (distinct, after the keep factor 0.842): 2^64.753, against the required
  excess x U = 2^64.319 (U = Li(B) ideal orbits), i.e. a margin of 0.43 bits, above the required
  0.35-bit margin (the pool's bootstrap standard deviation of the yield at this point is 0.054 bits,
  so the margin is 6.5 standard deviations).
- Relation collection: **2^119.377 units**.

## 5. Filtering and linear algebra

**Unknowns.** The sigma-orbits of the degree-1 prime ideals of norm <= B of both sides: U = Li(B) =
2^63.94. Relation orbits are the equations; the sieve produces excess x U of them.

**Filtering (H3).** Purge (singleton removal, then clique removal of the heaviest components of the
weight-2 graph until the excess is spent) and Markowitz merge with unit pivots and coefficients bounded by
255, as in CADO-NFS. Its effect is **simulated** on a random relation hypergraph with the production
occurrence profile: a column of norm q = B^t on side s appears in a relation orbit with probability
2 kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/rho(u_s) (the Dickman conditioning of B-smooth integers, the
same identity the yield model uses), columns distributed over t with the prime density B^t/t; columns
too heavy to be eliminated (t < 0.8) are charged as a dense per-row weight. This keeps the production
(column weight, row weight) profile exactly at any scale. Results (`filtsim.py`; part [5] runs scale 2^15
in the sandbox, the embedded table has scales 2^15-2^19 and two seeds): at this profile and excess, the
matrix has **N/U = 0.0076** rows per unknown at row weight <= 150 (the ratio is stable from 2^17 to
2^19). The model charges **N = 2 x (simulated N/U) x U = 2^57.89** (a factor-2 pad on the dimension,
i.e. x4 on the linear algebra) and row weight 150. Calibration against the only comparable record with the
same methodology (De Micheli-Gaudry-Pierrot, F_{p^6} 521-bit, exact B-smoothness, no excess): the
simulator reproduces the merge ratio (26-35% of the purged rows vs 29-33% reported) and the fraction of
ideals appearing (82% vs 87%), but keeps 61% of the rows after singleton removal where the record kept
38% (part [5]). The difference is attributable to their special-q sieving (which adds structure the random
model lacks) and to our norm-size estimates for their computation; it affects the purge step, whose
output our model does not use directly (the clique step spends the excess anyway), and it is the reason
for the factor-2 pad. The dependence on the smoothness parameters is strong (the embedded grid spans
N/U = 0.002 to 0.06), which is why the table is indexed by (u_f, u_g, excess).

**Square system.** After filtering the matrix is N x (N + k) with k = the number of Schirokauer maps
(unit ranks of the two absolute fields, 41 in v2-v4's count; <= 71 + 11 in any case). The sparse part S'
is made square by keeping N rows (the remaining rows are verification rows); the k dense SM columns are
treated as right-hand sides: block Wiedemann solves S' x_j = d_j for the k columns d_j simultaneously
(blocking factor n >= k, no extra products), and the k dropped rows give a k x k dense system for the SM
unknowns. No A^T A, no rectangular black box, no row scaling: the sparse matrix is multiplied as stored.

**Cost.** 3 N single-vector products of N x 150 entries at **35 instructions per entry** (shipped kernel
`la_mac.s`, 33 per entry branch-free, part [7] re-counts it; the per-row Barrett finish is amortised):
**2^115.442 units**. Rank: a random sparse square matrix over F_r with every column of weight >= 2 (purge
guarantees it) and small, non-degenerate coefficients is singular with probability O(1/r) beyond the
structural dependencies counted in H3 (identical-support weight-2 column pairs and the like are
eliminated by merge as 0 = 0 rows and cost one unit of excess each; their expected number at N = 2^58
is < 2^-4). The verification rows detect a wrong kernel vector; a failure is repaired by collecting 10%
more relations, which is inside the excess 1.3 chosen at the operating point (the yield margin is not
used for this).

**Schirokauer maps.** The SM values of 2^64.753 relations at <= 2^37 instructions each (exponent
<= 72 x 255 bits in O_f / r^2): <= 2^88 units.

## 6. Why no large primes

Large primes above B would raise the yield per point but put columns of expected weight < 1 into the
system; the simulation of section 5 shows the matrix size is governed by the top of the factor base
already, and the cofactorisation tier of 4.3 is the cheaper way to use the primes just below B. This
package therefore keeps exact B-smoothness.

## 7. Why nothing pairing-specific helps (and is not used)

Final-exponentiation preimages are polynomial-time (Kim-Cheon 2012/657) but the fibre has ~2^4314
elements and contains one Miller value; Miller inversion needs the unreduced Miller value (Satoh
2019/385), which the challenge does not reveal; hidden-root algebraic attacks have degree >= 2^359
(Vercauteren 2008/261); torus index calculus is L(1/2); Cheon-type restriction needs auxiliary powers,
which this track does not supply.

## 8. Individual logarithms: boot and descent (part [6])

Two logarithms (z and g); their ratio is a. The target is represented in F_p[y, x]/(h, g) (the isomorphism
from the pinned representation is computed once).

**Boot (initial splitting, g side).** For a randomiser t, T' = T gamma^t. With T'^{-1} = u + v theta,
(a - b theta) T'^{-1} has zero theta-component iff a v - b u = 0 in F_{p^6}: six F_p-linear conditions on
the 12 coefficients, so L_T' = {(a, b) : a - b theta in T' F_{p^6}^*} is a lattice of index p^6 (part [8]
verifies membership of every basis vector and |det| = p^6 exactly for the embedded example, the v2-v4 one
for this h). Every nonzero point of it lifts T' times a subfield element (r-log 0). After LLL (2^30 instructions,
once per t) the points of a ball of 2^37 candidates have g-side norms of ~2^2690 bits (measured over
random targets: mean 2661 at the smallest radius to 2723 at x256; the model uses 2725). A candidate is
tested by exact norm computation and a **full ECM run** for primes up to B1 (curves = expected x ln 100
for 99% success per factor; the GMP-ECM table up to 80 digits, B1 <= 2^265 so the table is not
extrapolated): expected candidates 1/rho(2725/log2 B1). The boot's g-side prime ideals above B become the
roots of the tree.

**Tree.** A node is a prime ideal of norm 2^l on side s. Its special-q lattice (12-dim, det 2^l) is reduced
once (2^30). The node targets a bound 2^n < 2^l: among 2^S lattice points of a ball (coordinates
2^(l/12 + S/12), norms N_f = 72 c + 6, N_g = 12 c + 6 log2 |u| + 2 with the q-side norm divided by 2^l), it
needs one with both norms 2^n-smooth; S is the least with 2^S P >= 4 (P from Dickman with alpha), so a
node fails with probability <= e^-4 (charged by the x2 node pad). Candidates are tested either by exact
norms + full ECM to 2^n on both sides (the only option when n >> 65: the factor base carries no
information about 2^n-smoothness then), or, when n <= 80, by sieving with the ideals <= 2^65 in the
special-q lattice (per-(node, ideal) 6-dim reduction 2^13 charged for every pair, 250 instructions per point,
slack k x n per side with the kept fraction and the survivor rate computed as in 4.3) and ECM on survivors.
**Every prime factor of either norm with norm in (B, 2^n] is a child** (expected count per 1-bit size bin
from the Dickman-conditioned density, x1.5 pad), and a child of size 2^l' costs the cost of an l'-node: the
dynamic programme over 1-bit levels from B up to 2^265 computes, for each (level, side), the cheapest
target and the expected total cost of the subtree, and then the cheapest boot bound B1. Children are
descended on **their own side**: an f-side child's special-q divides N_f, a g-side child's divides N_g;
no side is ever chosen after the fact. The tree ends at B exactly: a node's target n >= B and every child
is a factor-base ideal when n = B. Terminal ideals must have a determined logarithm: ideals of norm <= B
that appear in no relation have none; their expected fraction at this excess is < 2^-10 (the Poisson
occurrence rate at the top of the factor base is >= 1.5 x excess at this B), and the descent re-descends
such an ideal from level B + 1 (one more node, inside the x2 node pad).

Result (part [6]): per logarithm **2^96.09 units**, both **2^97.09**. The attack runs at most
**2^104 units** of descent per challenge, else aborts: by Markov, P(abort) <= 2^(97.09 - 104) =
0.0083.

## 9. Success probability

Relation collection and linear algebra do not involve the challenge. Given the relation count (H1-H3 with
the margin) they succeed up to the rank failure of section 5 (repaired by the excess). The descent succeeds
within its meter with probability >= 1 - 0.0083; z = 1 is free. Overall success >= 0.98 >= **0.39** over
uniform (Q, z) and the algorithm's coins (the challenge is exactly uniform over GT: P(ab = x) = 1/r for
every x).

## 10. Total work

| component | log2 units |
|---|---|
| relation collection (F sieved: updates, empty cosets, setups, scan, cofactorisation) | 119.377 |
| linear algebra (block Wiedemann, N = 2 x simulated, weight 150, 35 instr/entry) | 115.442 |
| descent, two logarithms, metered cap | <= 104 |
| Schirokauer maps | <= 88 |
| factor-base construction (2^62 ideals: roots 2^20 + 6-dim reduction 2^13 + tables) | <= 70 |
| canonical a-list (|E_a| x 595 units x ~60 instructions) | <= 60 |
| filtering passes, relation/matrix storage, ellipsoid fit, challenge pairing | <= 60 |
| verification (one pairing, subgroup checks, [a]P0, serialization) | <= 11 |
| **total** | **119.468044** |

The submitted bound **119.468045** is the full-ledger logarithm rounded up at 6 decimals; part [4]/ledger
prints `total_log2W_full_ledger = 119.468044` from the embedded data.

## 11. Sensitivity (printed by evidence.py)

- all sieve constants x1.5: total 120.014;
- matrix pad 4 instead of 2: total 119.712;
- the neighbourhood A +- 0.1, B +- 0.5 of the operating point (part [4]) is within 0.1 bit;
- the 0.35-bit yield margin is 6.5 bootstrap standard deviations of the pool estimate;
- without the unit correction (crediting every point of E, as v1-v4 did) the same model gives ~3.5 bits
  less, which is what v1 claimed and what this package shows to be an over-count.

## 12. Heuristics (all disclosed, with support and limitations)

- **H1 (smoothness).** Norms are exact resultants (part [2] re-derives them from seeds). B-smoothness
  probability: Murphy/Dickman rho(u) + (1 - gamma) rho(u-1)/ln N with GS21 alpha values recomputed for this
  package; rho from Marsaglia's exact series (part [1]). A model, not an execution; v1-v4 reviews accepted
  it at this scope.
- **H1' (yield margin).** The operating point is accepted only if the pool's distinct yield exceeds
  excess x U by >= 0.35 bit (2.4 bootstrap standard deviations of the embedded pool's estimate).
- **H2 (region and orbit capture).** f_c and the yield of F are measured on a held-out pool uniform over
  the region (seeds embedded, regenerated in part [2]); the unit test is exact; the only extrapolation is
  from pool statistics (2^13.5 canonical points) to the region.
- **H3 (filtering).** Random hypergraph with the production occurrence profile (section 5): the standard
  heuristic behind every NFS filtering prediction; simulated here, calibrated against one record with
  the stated discrepancy on the singleton step, padded x2 on the dimension.
- **H4 (kernel constants).** Every sieve constant is a retired count of a shipped kernel checked against
  brute force; the production slice (2^63 points) is larger than the benchmark slice (2^26-2^27), which
  the model handles by indexing the cost on points-per-coset, the quantity that fixes the walk's regime;
  the instruction metric prices loads/stores at one each (no cache or latency effects, as the metric
  prescribes). The full-size orchestration has not been executed, as for every NFS claim at this size.
- **H5 (descent).** Explicit dynamic programme with the stated ECM costs (GMP-ECM curve table within its
  range; 99% per-factor success; mulmod = 2.2 L^2 + 40 instructions for L limbs) and x2 node pads; the
  boot norm is measured on random targets of this field. The 0.39 requirement tolerates an expectation
  6.9 bits above the modelled one.
- **H6 (pairing transfer).** e(P0, Q) generates GT for nonidentity Q; standard.

Limitations: no NFS-DLP at this field size has been executed; the claim is an analytic bound grounded in
exact-norm simulation, measured kernels and simulated filtering, all of which are reproducible from the
package.

## 13. Memory (reported, not scored)

A job holds one slice: the bounding box of the largest half-slice is <= 2^66 bytes (one byte per point);
jobs are independent, so the peak is (number of concurrent jobs) x 2^66 bytes plus the streamed ideal
data (2^60.5 pairs x 288 bytes = 2^68.7 bytes, read-only, replicated or streamed). The
relation store is 2^64.753 x ~48 bytes = 2^70.3 bytes; the filtered matrix 2^57.89 x 150 x 5
bytes = 2^67.4 bytes; block-Wiedemann vectors < 2^67.4 bytes. Nothing is addressed across
nodes during sieving.

## Attribution

Research, validation and this package: Claude Fable 5.1 (harness Claude Code), building on the v1-v4
packages (Kimi K3 / Claude Opus 5.5, harness Claude Code) for the exact-norm pool machinery, the LA kernel
and the boot-lattice construction. Polynomial selection and alpha follow Guillevic-Singh (Math.
Cryptology 2021) and their public code; the filtering methodology follows CADO-NFS and De
Micheli-Gaudry-Pierrot (ASIACRYPT 2021); the FAPI-2 context follows Galbraith-Hess-Vercauteren
(2007/256), Kim-Cheon (2012/657), Satoh (2019/385) and Vercauteren (2008/261). No credentials or private
data are included.

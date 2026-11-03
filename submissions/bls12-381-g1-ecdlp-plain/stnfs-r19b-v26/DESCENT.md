# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Boot and nodes

Each target T becomes T' = T^e, e uniform in [1, r) (2/r charged for both
logs; identity direct). The boot lattice {(a,b) in Z^12 : aB + bA = 0 in
F_p[y]/h} (det p^6) has a certified LLL basis ([11]); the boot enumerates its
radius-4 ball (2^23.67 points), norms each and runs `ecm16.py` for
208.92-bit primes: a relation has every prime <= 2^215, <= 16 in (2^69,
2^208.92], <= 1 above. Exhausted balls draw a fresh e and lattice (2^36 each).
Distinct relations credit UDUP = 0.25; it stops at 16: 2^56.53 trials of
2^69.34 instructions over about 2^32.9 lattices.

A node is a special-q ideal of l bits on side f or g with its LLL-reduced
lattice L_q; its region is the integer ball of radius R (half counts in
`stats16.b85`). Sieve nodes enumerate every p in [2^8, 2^69] and root pair on
the ball (`dsieve16`, KERNELS.md); each primitive cell gets a certified
both-norm bound (2,400 charged); cells with residual masses <= T = kA nA + kB
n + 60 survive (exact norms, trial division, ECM). A relation has residual
primes <= 2^n, <= kA in (2^69, 2^nA], <= kB <= 1 above; these window primes are
children, charged at the most expensive bin at or below their window.
Retention (omitted mass, Hensel bounds, [11b]) >= 0.8575 at every sieve node.

## Yield

P = sum_{jA <= kA, jB <= kB} 1/(jA! jB!) int rho((m - sum t + alpha)/69) prod
dt/t is a lower sum (left-endpoint rho; mesh 384, rho spacing 1/1024, 200-digit
Decimal coefficients: numerical, not certified). Norms come from measured
offset bins (61,440 exact ball norms: cf sextile maxima, overall cg maximum,
+23/+6; boot 3,072 norms, +4). **Sub-cell integration (R18):** P is the
equal-weight sum over 24 cf-quantile sub-cells (`QBINS_SUB`, own cg maxima, top
edge pinned at 109.6/20.5) and ten boot deciles (`BOOT_SUB`, pinned at 16.7):
every cell of sub-cell k has cf <= cf_k, cg <= cg_k and P decreases in both, so
the sum lower-bounds the cell average. Boost over the sextile charge capped at
4x per node and 1.5x at the boot (14 of 62 nodes at the cap; boot
1.24x). Survivors, missing logs, retention and costs keep the sextile edges.
[11] repeats the sub-cell domination test live.

**Unit duplicates (red lines e, j; complete orbits).** Charged distinct yield
per cell is P x UDUP, UDUP = 0.25; support: E[shat/m] >= 0.25 P over uniform
ball cells v of fresh lattices at the shipped regimes, shat the semismooth
probability at v's exact norms, m(v) the number of units eps of O_h with eps v
in the ball (eps, -eps counted separately: conservative). m(v) is complete, not
a truncated box: eps v = c'B with ||c'|| <= R gives |sigma_i(eps)| <= M_i /
|sigma_i(a)| (M_i = R ||(sigma_i(B_k^a))_k||, Cauchy-Schwarz; also with b), so
every coefficient of eps is <= h(v) = max_j sum_i |Emb^-1[j][i]| T_i
(`mate_height_bound`, 2% float margin). All 2,442 units of height <= 8 are
enumerated exhaustively (`all_units(8)`; `gen_units(8)` reproduces them), so
cells with h(v) <= 8 get the exact orbit count; cells above contribute 0 to the
certified ratio, a lower bound. Our plain campaign (R18b, R18's regimes, 300
node / 400 boot cells per lattice, 95% CI lower ends of certified ratios vs the
integrated P): 69f 2.10, 69g 2.58, 100f 2.63, 100g 1.89, 150f 1.52, 150g 1.47,
180g 1.14, 185g 1.27, 190g 1.43, 195g 1.33, 200g 1.33, 205g 1.20, 215g 1.28,
220f 2.17, 220g 1.11; boot 0.79, 0.67, 0.77 (full counts >= 2.21 nodes, >= 0.74
boot; height-5 counts equal height-8 counts on >= 95% of yield mass). [11]
re-measures at this plan with `gen_units(8)` (2,442 asserted) on 48 fresh boot
cells and 16 cells of a fresh 225-bit (220, g) lattice and requires full ratios
>= 0.3 (boot) and >= 1.0 (node), certified >= 0.25 (`udup_relw`).

## Missing logs (red lines a, c, h)

Occurrence model: 2^63.97 relations (2-SE lower bound); an ideal I of norm q =
2^(69 t) on side s occurs Poisson(lambda_I) times, lambda_I = 2^63.97 x 2
kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/rho(u_s), u_s from the 2-SE lower
means 738.0/497.9; I is in the missing set M when it occurs at most twice:
m_I = P(Poisson <= 2), zero occurrences (e^-lambda) included in full, no
calibration factor (red line i). [6b]: purge losses 0.0470-0 per band vs model
0.7826-0.0017; no ideal with >= 3 occurrences lost. Ideals below 2^40 are
missing with probability < 2^-1000 (ignored).

**One shared M (red line h).** All descent relations see the same realised
set M. A relation of m = 69u bits contains I with probability w_I = kappa(u,
t)/q, so given M its missing count is S_u(M) = sum_{I in M} w_I (per
sigma-orbit C = {I, sigma I}, log sigma I = -log I: weight 2 kappa/q). The
recurrence uses (1 + 2^-8) mu_hi(u), mu_hi(u) >= E S_u = sum w_I m_I an upper
sum over 1/50-bit bands 2^40..2^69 (kappa at the band top, Poisson rate at its
minimum; within 1% of the integral; the base's 1/10-bit lower-edge sums were
5-8% low). **Chebyshev:** the indicators 1[N_C <= 2] are decreasing in the
occurrence counts, so under the Poisson/multinomial occurrence model (a fixed
number of relation draws) they are pairwise negatively correlated and Var S_u
<= sum_C w_C^2 m_C (`miss_pair(u, u)`); P(S_u > (1 + 2^-8) mu_hi) <= Var/(2^-8
mu_hi)^2. S_u grows with u, so one event per relation type (118, at their
largest norms) covers every relation: union <= 2^-37.14, charged 2^-30. No
association premise beyond the occurrence model is used. Off that event the
charged rate bounds every missing count; a thin node keeps a relation with
probability P(no missing ideal | M) >= prod_sides 1 - S + r (S^2 - c S)/2 -
S^3/6 (Bonferroni; r = kappa(u, 2)/kappa(u, 1)^2 bounds the pair ratio by
log-concavity of rho, grid-checked by [7]; `p_none`, 2-4% below exp(-S); taken
at the node's largest norms, where S is largest: S falls much faster than r as
u drops; this within-relation negative correlation is H1/Dickman, while
across ideals only the occurrence model is used). **H5c:** an ideal occurring at
least three times keeps its logarithm after purge and back-substitution ([6b],
simulation at 2^12: 0 of 3,816 such ideals lost);
P(Poisson <= 2) is charged at every lambda (uncalibrated where it is small,
lambda >= 8).

A node **discards** relations touching a missing ideal (forced for q <
2^95) or **descends** each. **Prime special-q (red line c):** a missing I is
an ordinary prime special-q node in its own lattice L_I (no composite lattice),
run with the first-bin entry (all primes <= 2^69, thin, so no recursion) and
charged per 5-bit sub-bin of [2^40, 2^69) and side at the sub-bin's largest
norms: work 2^110.261, failure 0.0200 (the composite charge was 2^112.214). [11]
adds 96 cells of fresh 2^45 (f) and 2^68 (g) lattices to the domination test.
**Repeats:** two visits to one I reuse its lattice, so a repeat is a failure.
Relations i, j share a missing orbit with probability sum_C w_C(i) w_C(j) m_C
(m_C once: one M), bounded per side by band sums (sum_band 1/N^2 <= 2^-b
sum_band 1/N) at the largest norms. A node tries at most one relation per
child plus one, so <= 2 N_MAX node and 2 kA boot relations: P <= C(2 N_MAX, 2)
(T_f + T_g) + 2 kA (2 N_MAX + 2 kA) T_boot = 2^-7.81 ([7]).

**Input propagation (red line a).** Fewer relations or smaller means raise
missing rates, work and failure monotonically (shipped plan):

| LOG2_RELS | means f/g | two-log | boot phi |
|---:|---:|---:|---:|
| 64.16 (estimate) | 744.3/499.0 | 2^114.1035 | 4.35e-04 |
| 64.07 (1 SE) | 741.2/498.4 | 2^114.1500 | 5.57e-04 |
| 63.97 (2 SE, shipped) | 738.0/497.9 | 2^114.2817 | 1.26e-03 |
| 63.88 (3 SE) | 734.8/497.4 | 2^114.6418 | 7.36e-03 |

## Recurrence and success (H5)

Distinct relations per node are negative binomial (mean mu, shape 0.5 for
windows of at most one prime, else 1; toy-ball fits [11b]: one-prime windows
0.60, larger 1.86-12.2). phi(node) = max((1 + mu (1 - psi)/r)^-r, 0.02), psi =
window sum of the children's phi; a node tries 1/(1 - psi) relations on
average; work = own + 1/(1 - psi) x children's work. The plan requires mu >=
16, phi <= 0.05, psi <= 0.8 ([7]). The three steps rest on different grounds:

1. **psi is a union bound.** A relation's subtree fails when any child (window
   primes, both sides at the window maximum, plus its missing-log ideals)
   fails; P(union) <= the sum of the children's failure probabilities, each at
   the worst bin at or below its window. No assumption on how the children's
   failures depend on each other.
2. **Boot: the 16 kept relations are independent trials given the
   main-phase table.** Each comes from the ball of its own uniform exponent;
   two from one lattice has probability <= C(16, 2)/2^32.9 = 2^-25.95,
   charged as failure. A relation's subtree is a deterministic function of its
   child ideals, the factor-base logarithm table and fresh ECM coins; the child
   ideals of the kept relations are functions of independent exponents. The
   table is shared but fixed once the main phase has run: conditional on it the
   subtrees are independent and P(all fail) is the product of the conditional
   failure probabilities. Each is at most the charged psi: psi is a union bound
   over the relation's children (step 1), and off the Chebyshev event (2^-30,
   charged) and under H5c the realised table's missing counts are at most the
   charged rates (relative margin 2^-16, `TABLE_SLACK`). Distinct subtrees reaching the same special-q or the
   same missing-log ideal are the coupling that survives conditioning; both are
   charged as failures in full.
3. **Within a node (H5a, the residual heuristic).** The relations of one node
   are distinct points of one ball; the thinning mu (1 - psi) and the
   1/(1 - psi) retries take their children's subtree successes as independent
   marks on the relation count: the independence of distinct cells'
   factorisations that every NFS relation-count estimate uses (H1), applied to
   their descents. Support: (a) distinct relations share children only through
   the repeat events charged above; (b) toy-ball counts fit the
   negative-binomial shapes used, and the charged shapes are below the fits;
   (c) the norm bins are worst-sextile envelopes over production lattices, so
   lattice-to-lattice variation is inside mu; (d) missing logs enter a node only
   by deterministic rejection or as prime-lattice children at the base bin.

**Headroom.** The meter is set for modelled success 0.41: the expected
two-log work may exceed its model value by 0.048 bits before the bound
reaches 0.39.

Node counter over both logs (missing-log nodes included): abort past N_MAX =
2^25.7 (balancing abort against repeats); E[N] <= 2^17.735 (`counts16.py`), so
P(abort) <= 2^-7.96. Window repeats: <= 2^-13.90 (all N_MAX nodes in the
smallest 2^70.3-ideal bin, 2^6 for 1/p). Repeat term (window and missing-log
pairs, missing-set tail charged 2^-30, boot shared lattice): 2^-7.78, charged in
full.

The 62-node plan (bins (69, 75) and 5-bit bins to 225; boot n1 = 215, nA =
208.92, kA = 16, R = 4) was coordinate-searched on this model (boot
per-relation failure <= 0.6732, all nodes sieve, R >= 192; 47 entries
differ from R17). Min mu 32.1, max psi 0.495 (boot 0.6588), max phi
0.0404. One log <= 2^127.570 instructions; **two logs <= 2^114.28169
units**; meter **2^115.0703** (Markov <= 0.578902); P(boot fails) <= 1.2578e-03.

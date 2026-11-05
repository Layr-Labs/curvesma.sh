# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [7c], [7e],
[7j], [7l], [7p], [7r], [7s], [7t], [7u], [7v], [11], [11b] and [12] check its
inputs.
1 unit = 20,000 instructions. Inputs: LOG2_RELS = 63.62 (the sub-region
pool estimate minus two pool standard errors, "Inputs"), mean norms
740.3/498.7 (at or below the sub-region's relation-weighted means), UDUP =
0.25, kernel charges SIEVE_HIT = 850 and BOUND_CELL = 2,000 (KERNELS.md).

## Targets and boot

Each target T is replaced by T' = T^e, e uniform in [1, r), log T = log T'/e
mod r: uniform on nonidentity mu_r (identity handled directly; success charges
2/r). The boot lattice {(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has
determinant p^6. The boot element a - b x has the single g-side norm
N(a^2 + U y b^2), whose terms balance at |b| = |a|/sqrt(U); the basis is
LLL-reduced with the b-block weighted by W = isqrt(U) and unscaled exactly
(`boot_lattice`; the lattice is unchanged) and certified exactly in the
weighted metric: membership, size reduction, Lovasz with delta = 0.98, Gram
determinant p^12 W^12 ([11]). The boot enumerates the radius-4 ball (2^23.67
points, one of each +- pair), computes each norm and runs the `ecm16.py`
schedule for 189.42-bit primes: a relation has every prime <=
2^205, <= 13 in (2^69, 2^189.42], <= 1 in
(2^189.42, 2^205].

**Uniform ball points.** `descent_witness.ball_point` keeps uniform points of
[-R, R]^12 in the ball, one per +- pair ([7v]: chi^2/df about 1; rounding a
real-ball point gives about 17). Norm offsets cg = log2 N_g - 12 c - 6 log2 U,
48 certified weighted lattices x 96 points/radius: p50/p90/p99/max
-191.6/-186.0/-182.1/-178.9 (R = 4), -191.5/-186.1/-182.0/-178.3 (R = 96);
esconc, 32 other lattices, 12,288 points at R = 4: -191.63/-186.01/-182.17/
-175.07. Charged bins (0.5: -187.5, 0.4: -181.9, 0.1: -175.0): p50/p90 of the
worse radius ceiled to 1 decimal plus 4; top decile at the largest of all
21,504 norms; stochastically above both samples, so a lower sum. Unweighted
(24 x 48 at R = 4): -7.2/1.0/5.9/8.5, about 184 bits higher. An exhausted ball draws
a fresh exponent and certified lattice (2^36 each). Distinct relations credit UDUP = 0.25 (below). The boot finds one
distinct relation (expected per_rel trials) and descends it before looking
for the next, at most 28 relations (lazy boot, below): 2^59.35 trials of
2^66.15 instructions over 2^35.68 lattices; failure needs all 28 to fail.
Two kept relations come from the same lattice with probability at most
C(28, 2) 2^-35.68 <= 2^-27.11, charged as a failure ([7]).

## Nodes

A node is a special-q ideal (l bits, side f/g) with lattice L_q (det q,
LLL-reduced). Its region is the radius-R integer ball in the reduced basis, one
point per +- pair; exact half counts are in `stats16.b85` (recomputed for R <= 8
by [7]). A **sieve** node enumerates every p in [2^8, 2^69] and every
relative-degree-one root pair on the ball (`dsieve16`, KERNELS.md); every
primitive cell gets a certified upper bound on both norms (72 + 12 linear
forms, outward rounding, charged 2,000; the forms' roots have certified
enclosures, [7r]); cells with residual mass (log2 bound minus score, minus log2
q on the special-q side) <= T = kA nA + kB n + 96 on both sides survive (exact
norms, trial division, ECM); non-primitive cells are skipped. A relation has
every residual prime <= 2^n per side, at most kA in (2^69, 2^nA] and at most kB <= 1 in (2^nA, 2^n]; these
window primes are the children, each charged at the most expensive bin at or
below its window.

## Yield and survivors

Per cell, the two-window semismooth probability P = sum_{jA <= kA, jB <= kB}
1/(jA! jB!) int rho((m - sum t + alpha)/69) prod dt/t is a lower sum (exact
mesh weights, rho at left endpoints). Norm sizes m come from offset bins over
61,440 exact ball norms on 160 production lattices (10 sizes 72-300 bits, both
sides, 6 radii): cf sextile maxima with the overall cg maximum, shifted by the
worst sampled per-(ideal, radius) mean minus the pooled mean (cf +23, cg +6); a
60-lattice sweep stays inside the bins by >= 7.25 bits.
mu = cells x P x UDUP x detection x retention (detection 1 - 0.001 x expected
factors). Survivors/cell use (T/69) rho((m - T)/69) per side at each bin's
lowest norms (margins (23, 6), special-q at the bin's upper edge).

**Unit duplicates (UDUP = 0.25).** Associates give the same relation.
Relation-weighted distinct/raw (41025889: 840 units in [-2, 2]^6, stratified,
24 production lattices at R = 384): 0.283-0.751 (mean 0.440); weighted boot
lattices 0.670-0.846; 500 -> 840 units: +0.008; [11] repeats 8 points. 0.25 is
below every measured value.

## Factor-base ideals without a logarithm: one ideal universe

Each side has the sigma-orbits of degree-one prime ideals of norm q <= 2^69,
counted per 1/100-bit sub-band [b, b + 0.01) as n = (Li(2^(b+.01)) -
Li(2^b))/2 (sum over b >= 40: (Li(2^69) - Li(2^40))/2 per side, [7u]). An
orbit of norm q = 2^(69 t) occurs in a relation of norm size u (in units of
69 bits) with probability w(u) = 2 kappa(u, t)/q, kappa(u, t) = rho(u - t)/
rho(u), and in the main phase Poisson(lambda) times, lambda = R 2 kappa(u_s,
t)/q, R = 2^63.62, u_s the side's mean norm/69 (smaller u_s gives fewer
occurrences: conservative). The same n, w and lambda are used everywhere
(missing children, the table event, the repeat bound). Within a sub-band lambda
lies in [lo, hi] (q and kappa at the edges) and every bound takes the worst
edge. The edge bounds read 1.2-1.4% above a 1/100-bit midpoint reference at
these inputs.

An ideal has no logarithm when it never occurs (probability e^-lambda; no
filtering step can recover it) or when it occurs but is purged and not
recovered by back-substitution. The charge is

  P_miss = (1 + 2^-8) e^-lo + (1 - e^-hi) c_occ(lo),

c_occ(lambda) = c_r for 1 <= lambda < 8 (0.172 / 0.117 / 0.057 / 0.02 on
[1, 2) / [2, 3) / [3, 5) / [5, 8)); for lambda >= 8, max(P(N <= 2 | N >= 1),
0.5 e^-0.45 lambda), so P_miss >= P(Poisson <= 2) there; 1 below
lambda = 1 (none occurs: lambda >= 1.12 here); c_occ is non-increasing.
c_r covers the measured unknown fraction of OCCURRING ideals after purge and
back-substitution: 16 `filtsim.missing_recon` runs (seeds 1-3 at 2^12..2^16,
seed 1 at 2^17; excess 1.1; 413-29,768 occurring ideals per band and run),
maximum (mean)
[1, 2): 0.1141 (0.103); [2, 3): 0.0775 (0.064); [3, 5): 0.0362 (0.028); [5, 8):
0.0071 (0.006), no trend over scales; c_r is the largest of 1.5 x the maximum,
the maximum plus three binomial SD at the smallest sample, and 1.3 x a hosted
rerun value (0.0149 in [5, 8)), rounded up. [7c] reruns at the reviewer's seed
at 2^12 against the covers with no allowance. Occ>=3 purge losses count as
unknown, so the covers include them. Above lambda = 8 the occurring losses
were measured at excess 1.1 on 200 + 120 fixed-universe draws at 2^12 and 2^13:
[8, 10) 66, [10, 12) 12, [12, 16) 3, [16, 24) 0, [24, inf) 0 lost among about
64,000 / 42,000 / 49,000 / 55,000 / 97,000 occurring ideal-draws, against 612 /
156 / 52.8 / 6.59 / 0.07 charged; below lambda = 24 charge/1.5 exceeds lost + 2
sqrt(lost) + 2 in every class. An m-bit side carries on average sum n w(m/69)
P_miss missing ideals (`miss_side`): 0.43-0.72 per f side, 0.41-0.66 per g
side, 1.40 at the boot's norm. Each node takes one option:

- **discard**: relations containing a missing ideal are rejected (mu x
  exp(-missing)); bins (69, 75] through (90, 95] always discard (THIN_MAX = 90);
- **descend**: each missing ideal Q (a prime ideal of norm <= 2^69) is a further
  child, sieved **in its own lattice L_Q as an ordinary prime special-q node**;
  its norms are at most those of the base bin (69, 75] on either side, so it is
  charged the larger base-bin work and failure probability of the two sides.
  Base-bin nodes discard, so the recursion closes. No composite lattices.

## The realized table (configuration-level control)

The main phase leaves one table: the set M of ideals without a logarithm,
shared by every relation of both descents. For a relation of size u on side s
the conditional number of missing children is S(M) = S_0(M) + S_1(M), S_0 =
sum_I w_I(u) 1[I never occurs], S_1 = sum_I w_I(u) 1[I occurs, I in M]. The
recurrence charges the band sums of (1 + 2^-8) e^-lo and (1 - e^-hi) c_occ(lo):
A_0(u) >= E S_0 and C_1(u). All its quantities (node yields through
exp(-missing), psi, children, work, E[N]) are non-decreasing in the missing
counts, so on the event

  E = {for both sides and every u = u4/4, u4 = 24..160: S_0 <= (1 + 2^-8) A_0,
       S_0 >= B_0 - 2^-8 A_0, S_1 <= C_1}

(B_0 <= E S_0: the lower-edge band sum) every conditional quantity is at most
its charged value. u4 = 24..160 contains every node, child and boot norm of the
plan ([7] checks). P(not E) is charged in full (DELTA_CONC = 2^-30).

**Model.** The ideal universe is fixed (the factor base is deterministic); the
main-phase relations X_1..X_n are independent (H1, the occurrence model above);
filtering (singleton removal, clique rounds that remove the heaviest 1% of
two-ideal components, back-substitution) is a deterministic function of the
relations. `filtsim` implements exactly this procedure.

**S_0.** With n independent relations P(N_I = N_J = 0) = (1 - p_I - p_J +
p_IJ)^n, p_IJ = P(one relation contains both orbits). Under H1 p_IJ <= p_I p_J
([7j]). Divisibility: an ideal over (q, t), t a root of h mod q, contains
a - b x iff a(t) = b(t) r mod q; residues at distinct (q, t) are independent
and uniform (CRT on O_h/q), so the density is the product; two orbits over one
(q, t) never both contain a primitive pair (one side: b(t) = 0 forces a(t) = 0;
across sides f, g need a common root mod (q, y - t), i.e. q | N(Res_x(f, g)) =
p^12 3^12, none in [2^40, 2^69]). Smoothness: the yield model is the product
rho(u_f) rho(u_g) (factor 1 across sides); on one side the factor is
rho(u-a-b) rho(u)/(rho(u-a) rho(u-b)) <= 1, log rho being concave ([7j]:
u in [1.5, 12], a, b in [40/69, 1]). So 1[N_I = 0] are pairwise negatively
correlated: Var S_0 <= V_0 = sum n w+^2 e^-lo, and Chebyshev gives
P(|S_0 - E S_0| > 2^-8 A_0) <= V_0/(2^-8 A_0)^2. Union over 274 pairs: 2^-34.29.

**S_1: Efron-Stein.** S_1 depends on the shared relation graph (purge cascades
and clique removal move several ideals together), so no correlation sign is
assumed. Instead:

*Lemma 1 (Efron-Stein, delete-one form).* If Z = f(X_1..X_n) with independent
X_i and each Z_i is any function of the X_j, j != i, then Var Z <= sum_i
E(Z - Z_i)^2. Proof: with E_i = E[. | X_1..X_i], Z - EZ = sum_i Delta_i,
Delta_i = E_i Z - E_{i-1} Z, orthogonal, so Var Z = sum E Delta_i^2. By
independence E_{i-1} Z = E_i[E^(i) Z] (E^(i) integrates out X_i only), so
Delta_i = E_i[Z - E^(i) Z] and Jensen gives E Delta_i^2 <= E(Z - E^(i) Z)^2 <=
E(Z - Z_i)^2, the conditional mean being the X^(i)-measurable L^2-projection.

Take Z_i = S_1 of the n - 1 relations other than X_i (exact delete-one rerun).
Every cross-ideal covariance, including purge cascades and changed clique
rounds, is inside Z - Z_i.

*Lemma 2 (influence).* Write w_I(u) = lambda_I g(t_I)/n with g(t) =
kappa(u, t)/kappa(u_s, t), monotone in t. Then |S_1 - Z_i| <= (G/n) Phi_i,
Phi_i = sum_I lambda_I |1[I in M_1] - 1[I in M_1 without X_i]|, G = max g, and
E S_1 >= (g_lo/n) Lambda, Lambda = E sum_I lambda_I 1[I in M_1], g_lo = min g
(over t in [40/69, 1]).

*Lemma 3 (tail).* With the premises below, E S_1 <= C_1/1.5, so P(S_1 > C_1)
<= Var S_1/(C_1/3)^2 <= 4 (G/g_lo)^2 n E[Phi^2]/Lambda^2 <= 4 (G/g_lo)^2
XI/N_s, N_s = 2^62.45 the side's orbits with t >= 0.8. Union over the 274
pairs: 2^-36.73 at XI = 2^14.

Then P(not E) <= 2^-34.29 + 2^-36.73 = 2^-34.04 <= 2^-30 ([7p]); the charge
holds for XI up to 2^20.6.

**Premises (H5c), each with a falsifying test.**

- *H5c-mean*: per lambda class, the production mean fraction of occurring
  ideals left unknown is at most the calibration maximum (c_r/1.5; for lambda >=
  8, c_occ/1.5). Support: the sixteen runs above (means 0.103 / 0.064 / 0.028 /
  0.006 below the maxima 0.1141 / 0.0775 / 0.0362 / 0.0071) and the
  fixed-universe draws above lambda = 8; the production relation count gives the
  lambdas of toy excess 1.19. Within a class w and lambda both scale as
  1/q, so the weighting favours the higher-lambda, lower-loss ideals. [7c]
  checks c_r >= 1.5 x the premise means and reruns the measurement at the
  reviewer's seed against the covers.
- *H5c-ES*: the deletion influence is local: XI = N_s n E[Phi_s^2]/Lambda_s^2
  does not grow with the system. Exact delete-one reruns at excess 1.1 (a C
  transcription, checked equal on 70 deletions): single systems 2^11..2^16
  (13 systems): 288-880 per side; pooled over 8 systems at 2^8, 40 seeds:
  82-3,350 (median 508). XI = 2^14 is 18.6 x the largest single-system value
  and 4.9 x the largest pooled one; [7e] measures XI at the reviewer's seed
  and fails above 2^14. Production ideals with t < 0.8 have lambda > 2^13
  (charged < 2^-1000 each); in 11 excess-1.1 systems no high-lambda column
  flipped and the largest flipped lambda was 10.8.

**Teeth.** On fixed universes at excess 1.1 the realized missing mass's variance
exceeds the independent-indicator value 1.09-2.37 times (0.9-5.8 SE); the
Efron-Stein estimate is above the variance in every row ([7t]).

**Boot and recurrence under the table.** On any table in E, the 28 boot
relations come from independent exponents and distinct lattices (shared
lattices charged above), each subtree a deterministic function of its child
ideals, the table and fresh ECM coins, so P(all fail | table) <= psi^28;
E[q(M)^28] <= psi^28 + P(not E), P(not E) charged once for both logs; the
charged boot failure keeps relative margin 2^-16 above psi^28 (`TABLE_SLACK`).
The meter and node counter use E[W | table] and E[N | table], bounded on E by
the charged recurrence.

## Repeats and the node counter

A node counter aborts past N_MAX = 2^28.5 nodes (both logs; a counted
failure). V <= N_MAX, so E[pairs] <= N_MAX E[N]/2, E[N] <= 2^17.48
(`counts16.py`): P(abort) <= 2^-11.02. Special-q repeats (2^6/|smallest bin|
per pair, bin 2^70.30 ideals) <= 2^-19.33. Two visits from distinct relations
pick the same ideal with probability sum_{I in M} pi_1(I) pi_2(I), pi(I) =
w_I/S(M) <= w_I/(B_0 - 2^-8 A_0) on E (M is shared: one factor m_I, not two);
averaged over M at most sum_I m_I wh_I^2 = 2^-55.71, so repeats <= 2^-10.73.

## Retention and cofactorisation

**Counted cells.** A cell is counted only if base-primitive (for every p and
irreducible h of H mod p, a and b not both divisible by (p, h)); base-primitive
fraction prod (1 - p^-2 deg h) = 0.9100; cells divisible by a relative-degree-one
ideal above a base prime of degree >= 2 with p >= 2^8 are excluded
(0.9523/0.9919). The factor base holds every prime ideal of norm <= 2^69 above
p < 2^8 plus the degree-one ideals above p in [2^8, 2^69].

**Omitted mass.** A counted cell is rejected only when one side's omitted
mass X exceeds 96 - 1 - m/128. X sums independent local Hensel bounds
(`retain16`, tables in `stats16.b85`, recomputed by [11b]): P(X > 96) =
1.5e-12 (f) / 2.5e-10 (g) per cell; the lost fraction is bounded by num/den at
the node's largest norms. Retention = (1 - loss_f - loss_g) x base-primitive x
(1 - ordinary exclusions). Every sieve node charges RET_FLOOR = 0.8592, the
boot COUNTED = 0.8592; [11c] recomputes every node's exact value (minimum
0.859250, largest one-side loss < 1e-5) and the counted fraction (0.859251).

**ECM (`ecm16.py`).** Suyama-form Montgomery curves, cheapest (B1, B2/B1) grid
point with wheel 30030; curves = ceil(ln 0.001 / ln(1 - P1)) at 367.76 (L/5)^2
per multiplication, P1 the Dickman lower sum on the curve order ([12] prints
the schedules and matches P1 on random 36/40-bit primes).

## Lazy cofactorisation and lazy boot

Relations are tried one at a time and the first whose subtree succeeds ends the
node; the predecessor still cofactored every survivor and collected all 16 boot
relations first. This revision does that work on demand; sieve, survivors,
windows, phi/psi and the node counter are unchanged, and every failure event is
the same.

- **Nodes.** After the sieve (charged in full), the s survivors are cofactored
  in uniformly random order (a fresh random permutation, independent of the
  relations' contents and outcomes by construction); each distinct relation
  found is descended at once and the node stops at the first relation whose
  children all succeed (or fails when the survivors run out, exactly the old
  failure event). Given M distinct relations among the s survivors, the k-th
  one first appears at expected position <= k (s + 1)/(M + 1) (the first
  occurrence of a class is at or before a uniformly chosen member, and M such
  members occupy a uniform M-subset). The number J of relations tried is
  independent of the order with P(J > j) <= psi^j (the retry premise of step 3
  below, unchanged), so E[min(J, M)] <= 1/(1 - psi) and

  E[processed]/s <= E_M min(1, (1 + 1/s)/((1 - psi)(M + 1)) + psi^M),

  with M negative binomial (mean mu, shape r) as in phi_nb.
  `descent16.lazy_fraction` evaluates this sum; node own work is LLL + sieve +
  (survivors x per-survivor cost) x that fraction. [7l] also runs the process
  itself (`lazy_mc`: NB count, uniform positions, independent outcomes) at the
  three costliest nodes and checks that the simulated mean is at most the
  bound plus a Monte Carlo tolerance of 0.01 (simulated 0.1391/0.0966/0.0706
  vs bounds 0.1563/0.1146/0.0740): it checks the arithmetic of the bound under
  H5a, not H5a itself.
- **Boot.** The boot finds one distinct relation (expected per_rel trials, the
  old per-relation rate) and descends it before looking for the next, at most
  28 relations. The trials are i.i.d. and the number of relations found,
  min(J, 28), is a stopping time, so by Wald's identity the expected trial
  work is per_rel x (1 - psi^28)/(1 - psi) x trial (3.558 relations instead of
  28). Failure is psi^28 (all 28 fail).

The recurrence W = own + retry x children is otherwise unchanged.

## Recurrence and success (H5)

Distinct relations per node are negative binomial (mean mu, shape r): r = 0.5
for effective windows of <= 1 prime (base nodes have empty windows), r = 1
otherwise. Shapes come from toy-ball fits (empty-window rows in `stats16.b85`,
printed by [11b]): one-prime windows fit 0.60, two/three-prime 1.86-12.2; of 12
empty-window configs, all with mean >= 7.4 show 0-4 ideals without a relation
vs 5.4-12.4 predicted at r = 0.5.

A child subtree fails with probability phi; phi(node) = max((1 + mu (1 -
psi)/r)^-r, 0.02), psi = window sum of the children's phi; a node tries at most
1/(1 - psi) relations on average. Node work = own work (LLL + sieve + the lazy
fraction of the survivor cofactorisation) + 1/(1 - psi) x children's work. The
plan requires mu >= 16, phi <= 0.05, psi <= 0.8 everywhere ([7] checks). The
three steps rest on different grounds:

1. **psi is a union bound.** A relation's subtree fails when any of its
   children (at most 2 kA + 2 kB window primes, both sides counted at the
   window maximum, plus its missing-log ideals) fails; P(union) <= sum of the
   children's failure probabilities, each taken at the worst bin at or below
   the child's window, with no assumption on how they depend on each other.
2. **Boot: independent trials given the table:** the kept relations' child
   ideals are functions of independent uniform exponents; the table is
   controlled by E (charged 2^-30); shared lattices and repeats are charged.
3. **Across the relations of one node (H5a, the residual heuristic).** The
   thinning mu (1 - psi) and the retry count take the subtree outcomes of a
   node's distinct relations as independent marks. Established: a subtree is a
   deterministic function of its child ideals, the table and fresh ECM coins;
   distinct relations share a child only through charged repeat events; each
   child is descended in its own prime lattice; the shared table is controlled
   by E. Not established: that distinct children found in one parent lattice
   have uncorrelated outcomes. H5a is a heuristic; its consequences are bounded
   below.

**Headroom and robustness to correlated outcomes.** The meter is set for
modelled success 0.700122 (the claim is 0.39): E[W] may exceed its model value
by 1.028 bits (a factor 2.04) before the bound reaches 0.39. [7s] also
evaluates the plan under a common-shock alternative to H5a: at every node, with
probability delta all relations fail together and every survivor is then
cofactored and every relation tried (retry (1 - delta)/(1 - psi) + delta mu,
lazy fraction (1 - delta) x bound + delta, phi = delta + (1 - delta) phi_nb);
success at the shipped cap stays >= 0.39 for delta <= 0.0023 at all nodes at
once. [7l] checks the lazy bound's arithmetic under H5a, not H5a. The charged
shapes (0.5/1) are 17% and 46% below the toy fits (0.60/1.86); [7s] also prints
the recurrence at the fitted shapes, at the estimate and one standard error,
at UDUP = 0.283 and at c_r = 1.

## Inputs ([4], [7s])

LOG2_RELS = 63.62 is the sub-region pool estimate 2^63.6922 (proof.md
section 3) minus two pool standard errors (2 x 0.0356 bits), floored to 2
decimals; [4] checks log2_rels - 2 SE >= 63.62. The same relation count
feeds the missing-log rates, the descent work and the success bound. At c_r =
1 the shipped plan is infeasible ([7s]): the calibration is load-bearing (H5).

## The shipped plan

62 nodes (bins (69, 75) plus 5-bit bins 75-225, both sides) and the boot (n1
= 205, nA = 189.42, kA = 13, R = 4, 28 relations). The node plan is
4ff4138b's (coordinate-searched there) with two entries moved for these inputs:
(155, g) iA 236 -> 232 and (205, g) n 160 -> 164, which restores phi <= 0.05;
the boot was re-searched under the lazy model (n1 210 -> 205, iA 356 -> 340,
16 -> 28 relations). [7] checks mu >= 16, phi <= 0.05, psi <= 0.8 everywhere.
Quadrature: mesh 384, rho spacing 1/1024 (left-endpoint lower sums). At the
inputs above: min mu 25.4, max phi 0.0485, max psi 0.719; all 62 nodes sieve at
R = 192..384. One log <= 2^128.040 instructions; **two logs <= 2^114.7514
units**. The meter stops the descent at 2^116.495 units (Markov: E[W]/cap <=
0.298612). P(a boot fails) <= 0.00009717 ([7]). Success >= 1 - 2 x 0.00009717
- 0.298612 - 2^-11.02 - 2^-19.33 - 2^-10.73 - 2^-27.11 - 2^-30 - 2^-48 -
2^-65.54 - 2/r >= 0.700121 (2^-65.54: Block Wiedemann plus main phase with
128-bit X entries, proof.md section 7).

## Supporting measurements

[11] re-samples uniform boot and special-q points per run against the bins;
[11b] the retention inputs and empty-window toys; [12] the ECM schedules.

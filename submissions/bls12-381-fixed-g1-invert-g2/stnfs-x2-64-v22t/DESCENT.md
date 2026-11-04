# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T is replaced by T' = T^e, e uniform in [1, r): uniform on nonidentity mu_r (distance 1/r; identity z direct; ledger charges 2/r). The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6, one (g-side)
norm N(a^2 + U y b^2), U = |u| = 2^63.71; the basis is LLL-reduced
in the weighted metric (b-block x W = isqrt(U) = 2^31.86, `boot_lll`) and returned unscaled. Exact certificate: membership,
size and Lovasz (0.98) in the weighted metric, Gram determinant
p^12 W^12 ([11], `certify_boot`). On 8 certified lattices x 384 uniform
R = 4 points: g-norm offsets cg mean -191.8, p90 -186.0, p99 -182.0, max -178.6 (Euclidean: mean -7.3, max 12.6; spread 0.53). The boot
enumerates the radius-4 ball (2^23.67 points, one per +- pair) with the `ecm16.py` schedule for 199.79-bit primes: a relation
has every prime <= 2^215, <= 15 in (2^69, 2^199.79], <= 1 in (2^199.79, 2^215].
An exhausted ball draws a fresh exponent and lattice (2^36 each).
Distinct relations credit UDUP = 0.65. It stops at 24 distinct
relations (BOOT_RELS): 2^52.84 trials of 2^67.75 instructions over ~2^29.1
lattices; failure needs all 24 to fail.

## Nodes

A node is a special-q ideal (l bits, side f/g) with LLL-reduced lattice L_q
(det q). Its region is the radius-R integer ball in the reduced basis, one
point per +- pair; exact half counts are in `stats16.b85` (recomputed for
R <= 8 by [7]):

- **sieve**: every p in [2^8, 2^69] and every degree-one root pair is
  enumerated by `dsieve16` (KERNELS.md); records go to the owners' 2^40-entry
  score arrays. Every primitive cell gets a certified upper bound on both norms
  (72 + 12 linear forms, outward rounding; ~1,500 instructions, charged 2,400).
  Cells with residual mass (log2 bound minus score, minus log2 q on the
  special-q side) <= T = kA nA + kB n + 96 on both sides survive (exact norms,
  trial division, ECM). Non-primitive cells are skipped.
- **ecm**: ECM trials on every ball point (never cheaper in the shipped plan).

A relation has every residual prime <= 2^n per side, at most kA in (2^69, 2^nA]
and at most kB <= 1 in (2^nA, 2^n]; window primes are the children, each
charged at the most expensive bin at or below its window.

## Yield and survivors

Per cell, the two-window semismooth probability

  P = sum_{jA <= kA, jB <= kB} 1/(jA! jB!) int rho((m - sum t + alpha)/69) prod dt/t

is a lower sum: each mesh cell carries exact weight int dt/t and rho at its
left endpoint (rho decreases); cells crossing sum t = m are dropped. Norms m
come from measured offset bins over 61,440 exact ball norms on 160 special-q
lattices (10 sizes 72-300 bits, both sides, 6 radii): cf sextile maxima with
the overall cg maximum, shifted by worst sampled per-(ideal,radius) mean minus
pooled mean (cf 20.68+1.58 -> +23, cg -2.53+7.81 -> +6). Boot norms use 3,072
exact norms on 8 certified
skewed lattices at R = 4 (p90/p99/max + 4 bins). The probability integrates
over 24 norm sub-bins at nodes, 10 at the boot (`QBINS_SUB`/`BOOT_SUB`),
charged at most 4x / 1.5x the sextile value (UDUP support is against this pr).
mu = cells x P x UDUP x detection x retention (UDUP = 0.65; ECM detection
1 - 0.001 x expected factors). Survivors per cell use (T/69) rho((m-T)/69) per
side at each bin's lowest norms (previous sextile's cf max, cg min -54.56,
margins (23,6); special-q at the bin's upper edge).

## Factor-base ideals without a logarithm

Under the occurrence model (`filtsim.py`: 2^63.589 relations, the conservative
input: excess-1.1 theta-2347/8000 estimate 2^63.689099 minus 2 pool SE, [7s];
certified product 0.968032661 applied), an ideal of norm q = 2^(69 t) on side s
occurs with probability 2 kappa_s(t)/q, kappa_s(t) = rho(u_s-t)/rho(u_s), u_s
from mean norms 738.2/498.4 ([4]). Occurrences are Poisson(lambda_s(b)); occ
<= 2 is missing. An m-bit side carries on average sum_b kappa(m/69,b/69)
ln((b+0.1)/b) P(Poisson <= 2) missing ideals (1/10-bit bands over the full base;
b < 40 terms are exactly 0.0 in binary64, lambda >= 2^23; [6c] checks the tail
is 0.0 at the shipped input and 63.5):
0.20-1.84 f-side, 0.40-3.68 g-side over u = 7..39 (3.74 at the boot's
2,709-bit norm). The proxy (coefficient 1) charges at least the e^-lambda
mass per band ([6b-zero]). Each node takes one option:

- **discard**: relations with a missing ideal are rejected (mu x exp(-missing)); the 34 thin nodes discard, so missing
  ideals never recurse there: all 12 with q < 2^90 (THIN_MAX 85 forces thin: the 8 v22s ones plus the 4 new
  (2,40)/(40,69) leaf bins below) plus the 22 plan-flagged thin=1 nodes
  (90f/g, 95g, 100f/g, 105f/g, 110f/g, 115f, 120f/g, 130f/g, 135f/g, 145f/g, 170f, 175f, 205f, 220f; exact list in
  refinement.json `declared_thin`, audit-bound to the `node_cost` flags);
- **descend (v22t: direct)**: each missing d1=1 ideal Q is a further child at the
  other 32 nodes, sieved in its own lattice L_Q with its side's leaf-bin entry ((2,40) R96 for
  N(Q) < 2^40, (40,69) R256/288 otherwise; routed by norm, exact), failure bounded by pm
  (S* charge, [11s]); the boot descends likewise. Missing d1>1 ideals are never descended
  (d1>1-discard haircut below). No auxiliary ideal exists.

Cascading removal can also lose occ>=3 ideals (e.g. an ideal
in 3 rows each holding a different singleton is purged); [6b] measures losses on the [6] system after purge + clique + merge (weight 200) + final-matrix-seeded reconstruction (soundness below): per lambda band
(ideals, measured, model, R1-form/other/merge): [1,2): 897, 0.0970, 0.7938, 0/2/0; [2,3): 961, 0.0822, 0.5602, 1/1/0; [3,5): 757, 0.0343, 0.2761, 0/0/0; [5,8): 454, 0.0044, 0.0628, 1/0/0; [8,1e+09): 667, 0.0000, 0.0016, 0/0/0;
model exceeds measurement by >= 6.81x wherever loss occurs (merge loses nothing).

**Reconstruction soundness.** `missing_recon` seeds known logs from
the FINAL matrix only, then applies
elimination pivots in reverse order and original rows to fixpoint. Pivot rows
are true relations (elimination preserves row space); triangularity holds (an
eliminated column is deleted everywhere, never reintroduced), so
each pivot row carries at most one unknown in reverse order (complete, no
circularity; machine-checked on every campaign system). Original rows are true
relations (one-unknown rule sound); deleted rows are on-disk truth, usable.

## Table-good event (v22 joint control + v22t d1 control)

S_j(M) = missing mass of node j's relations given the realized post-filtering
table M (67 functions: 66 nodes + boot); S_j = S_occ,j + S_res,j (occ<=2 pre-filtering; residual
= filtering-lost occ>=3). Descent runs on G = G_occ AND G_res AND G_hi AND G_d1:
- G_occ (McDiarmid + union): max_j (S_occ,j - 1.05 E_j) <= 0.01. Per-draw
  difference <= 2*WMAX*wmax (WMAX = 1500, pool bits <= 1400 live; wider widths
  are the disclosed H5 width premise; wmax = 400/R covers kappa ratios <= 4;
  lambda <= 100): P(S_j - 1.05E_j > 0.01) <= exp(-2(.01)^2/sum_c2); union over
  67: [6c] -1926.9 + log67 <= ln 1e-9, charged 1e-9. Rigorous under H5.
- G_res (Markov): boot-weighted residual C(M) <= 0.14 (boot weights dominate,
  kappa ratio <= 0.865, so all S_res,j <= C(M)).
  E[C] = 0.03265180 (Wilson-95 rates x boot masses over 8 <= b < 69, plus exact
  below-2^8 d1=1 sum 0.000327, v22t); P(G_res^c) <= 0.2332, charged 0.24.
  H5-res (disclosed rate transfer): pooled rates over 2^12..2^17 are
  0.00318/0.00205/0.00230/0.00234/0.00236/0.00231 (flat over 32x; counts in
  refinement.json, bands in the note); each charged upper covers its 2^16/2^17
  point rates (`ok_cover`). Tree obstructions (3 singleton-tainted rows) are
  scale-free given occurrences, included in the measured residual.
- G_hi (union): no lambda > 100 ideal occ<=2 (<= 2^(b+6)/band x
  e^-lambda(1+lambda+lambda^2/2); [6c]: 8.78e-20, charged 1e-9).
- G_d1 (v22t): relations touching a missing d1>1 ideal are dropped (yield x
  exp(-B), B = 1.05e-6+0.022+0.02 = 0.042). G_occ_d1: max_j(S_d1,j-E) <= 0.022
  (McDiarmid [6c] -28.9+log67 <= ln1e-9, charged 1e-9; E_occ_d1 <= 1e-6 proven).
  G_res_d1: d1 residual <= 0.02 (Markov; E_d1 = 8.768e-5 = below-8 6.771e-5 +
  recorded sieve 5.4e-6 + tail 1.456e-5; delta 0.00438 <= 0.005). G_hi_d1mid: no
  occ<=3 d1>1 ideal with q <= 2^59 ([6c] 4.85e-8 <= 1e-7, charged at value;
  pi(x) <= 1.26x/ln x). d1>1 rates are H5-res-d1 (disclosed: H5-res uppers apply
  by degree; filtering acts on columns). The sieve is a RECORDED upper (live
  redo infeasible, ~25 CPU-min; script + output + PARI xcheck ship in relief).
P(G^c) <= 0.245 + 3e-9. On G, miss = 1.05 E_occ + 0.15, every mu carries
exp(-0.042), boot failure <= psi^24 (fresh exponents given M); E[W] <= W_infl +
P(G^c) cap. Success subtracts both deltas (proof 7).

## Retention and cofactorisation

**Counted cells.** A cell is counted only if base-primitive (for every p and
irreducible h of H mod p, a and b not both divisible by (p, h)); other cells
repeat another cell's relation. Base-primitive
fraction prod (1 - p^-2 deg h) = 0.9100 (exact degrees for p < 2^9). Cells
divisible by a relative-degree-one ideal above a base prime of degree >= 2
with p >= 2^8 are excluded by ordinary fractions 0.9523/0.9919 (f/g; joint survival is the Frechet lower bound ordf+ordg-1, valid under overlap).
The factor base holds every prime ideal of norm <= 2^69 above p < 2^8 plus the
degree-one ideals above p in [2^8, 2^69].

**Omitted mass.** A counted cell is rejected only when one side's omitted
mass X exceeds s = 96 - 1 - m/128
(rounding: each sieved hit loses < 1/16 bit to floor(16 log2 p)). X sums independent local contributions (whole valuation at
p < 2^8, excess powers at p >= 2^8), each with one local bound.
Let P = (p, h_i) be a prime of the base field (Kummer-Dedekind, O = Z[y]) with
residue field F_q, q = p^deg h_i, and S the side polynomial, monic in x. For a
base-primitive pair, (a : b) is a point of P^1(O_P): b is a unit with
probability q/(q + 1), otherwise a is one and S_hom(a, b) is a unit. With b a
unit, v_P(S_hom(a, b)) = v_P(S(t)) for t = a/b uniform in O_P. If t does not
reduce to a root of S in F_q, the valuation is 0. A root rho of multiplicity mu
in F_q lifts to exactly mu roots in the completion (Hensel), so v_P(S(t)) is the sum of the mu distances v_P(t - rho_i),
and v_P(S(t)) >= j forces some v_P(t - rho_i) >= j/mu, an event of probability
at most q^-ceil(j/mu):

  P(v_P >= j and t reduces to rho) <= min(1/q, mu q^-ceil(j/mu)),

with equality q^-j when mu = 1. The primes above p are
independent (CRT) and v_p(N) = sum over P of deg(h_i) v_P. Degree,
ramification and f/g root multiplicities of every P
with p < 2^9 are tabulated (`retain16.structure`, `stats16.b85`); [11b]
recomputes the table at repeated-root/ramified primes (2, 19,
23) and a seventh of the others, and checks the bound against residue-class
samples at 2 and 23. All multiplicities are 1 except g at the two primes above 2
(2^16 | U: the double root x = 0, q = 4 and 16) and f at the degree-4 prime
above 19 (two double roots among four, q = 19^4); 23 is ramified in the base
field (e = 2) with g having two simple roots at each of the two degree-one
primes. Above 2^9 the excess powers are bounded per dyadic block by 72 p^-2
events per level over pi(x) <= 1.26 x/ln x primes.
Masses are rounded up to 1/4 bit and probabilities below 1e-30 are moved to the
largest kept mass: means 3.35 (f) and 7.66 (g) bits, P(X > 96) = 1.2e-12 (f)
and 2.3e-10 (g) per cell.

**Relation-conditioned loss.** Under H1 a cell with omitted mass x is
semismooth with probability Pss(m - x), so the per-side lost fraction is

  loss = sum_{x > s} P(X = x) Pss_hi(m - x) / sum_{x <= s} P(X = x) Pss_lo(m - x),

bracketed by monotonicity on a 64-bit grid in x at the node's largest norms
m_f, m_g. Retention = (1 - loss_f - loss_g) x base-primitive x (1 - ordinary
exclusions); [7] evaluates every sieve node (min 0.8593; largest one-side loss
< 1e-5).

**ECM schedule (`ecm16.py`).** Suyama-form Montgomery curves, cheapest (B1,
B2/B1) grid point with wheel 30030; stage-1 x-only ladder,
stage-2 cross products, one binary gcd per stage. Per-curve success
P1 is the Dickman lower sum on the curve order; curves = ceil(ln 0.001 /
ln(1 - P1)) at 367.76 (L/5)^2 per mult. [12] prints node/boot schedules (boot:
B1 = 2^28, 1418167 curves, 2^67.96/trial) and matches P1 on random 36/40-bit
primes.

## Recurrence and success

Distinct relations per node lattice are negative binomial around the lattice's
own mean (H5w, proof 9): r = 0.5 for windows of <= 1 prime (base nodes: empty),
r = 1 otherwise. Toy-ball fits (empty-window rows in `stats16.b85`, [11b]):
one-prime 0.60, two/three-prime 1.86-12.2. Of 12 empty-window configs, all with
mean >= 7.4 show 0-4 relation-less ideals vs 5.4-12.4 at r = 0.5; the three at
3.2-4.0 show 19/11/21 vs 18.1/17.2/16.4, cleared by raising lb (0 vs 7.3-11.2).
Pooling adds only dispersion. Base nodes run at mu = 769.6 and 649.5.

A child subtree fails with probability phi; phi(node) = max((1 + mu (1 -
psi)/r)^-r, 0.03), psi = window sum of the children's charges (union bound, no sibling independence); a node tries <=
1/(1 - psi) relations on average. Node work = own sieve/ECM work + 1/(1 - psi)
x children's work. The plan requires mu >= 18, phi <= 0.054, psi <= 0.8
everywhere ([7] checks; the phi cap is a self-imposed feasibility gate, raised v22s 0.05 -> v22t 0.054 for the
evaluated max 0.05323 + margin — no soundness bound consumes it). The descent keeps a node counter over both
logs (cap N_MAX = 2^27.3, counted failure); E[N] <= 2^17.91
(`counts16.py`, unrounded 2^17.9093), so P(abort) <= 2^-9.39 (Markov). All N_MAX nodes in the
smallest bin (2^70.3 ideals) with 2^6 for 1/p weighting: shared-ideal
probability <= 2^-10.70, charged in full ([7]).
A factor-base ideal Q without a log
(norm <= 2^69, relative degree 1) is descended directly in its own lattice L_Q with the leaf-bin entry of its
side ((2,40) or (40,69)); its node is charged pm, the worst thin-bin failure probability ([11s]). Ideals of
relative degree > 1 are never descended (d1>1-discard haircut, table-good event above). No auxiliary ideal exists.

**Per-node exceptional mass (S\*).** For a node class c (bin/side or small-Q
stratum) let charge_c be the least parent-slot charge (phimax, or pm), lambda_L the
lattice's own model mean and G_L = lambda_L/mu_c. A class-c node fails with
selection-averaged probability at most P(G_L < x0) + NB0(x0 mu_c (1 - psi_c), r_c)
for every x0 > 0 (NB0 decreases in its mean); when this is at most charge_c per
class, each psi and the boot's psi^24 hold as computed (no Jensen step). [11s]
checks the 75 classes reachable from the boot (47 bins + 28 strata; 160f-220f,
215g, 220g never visited) on a pre-registered campaign: 600-1200 lattices per
class (56,000), each a uniform degree-one ideal at its bin top (strata: top
1-bit slice, exact enumeration at top <= 2^11), `lll`-reduced, K = 200 uniform
cells, Y = s/(1+m/2)/(pr x 0.65) over the complete PARI-certified unit group.
Flagged when mean of min(Y,b) < m; a G_L < x0 lattice escapes with probability
<= beta = exp(-3K(m-x0)^2/(2b(2x0+m))) (Bernstein), so P(G_L < x0) <=
CP(k,n)/(1-beta) at 1 - 0.05/75 (Bonferroni 95%). Result: 13 of 56,000 flagged;
every class passes; tightest 69g (+0.0049). Transfer: given its norm range a
child is uniform degree-one (H1); top-of-bin/slice sampling is conservative
(residual norms grow with q; uniform overweights harder ideals vs 1/q). Boot
spread changes work only (linear in mean yield; 3000-cell recount).

The shipped plan has 66 nodes (bins (2,40), (40,69), (69,75) plus 5-bit bins
75-225, both sides) and the boot (n1 = 215, nA = 199.79, kA = 15, R = 4): the
90d99567 plan re-optimised by single-then-joint coordinate search at the sub-bin
pr, skewed boot and LOG2_RELS = 63.599 (41 of 62 move, `plan_v21`; boot window
unchanged) plus the 4 small-Q leaf bins. At LOG2_RELS = 63.589, mean norms
738.2/498.4, UDUP = 0.65, THIN_MAX 85, PHI_FLOOR 0.03, BOOT_RELS 24, on the
table-good event: min mu 33.6, max phi 0.0532, max node psi 0.544, boot psi
0.7170. All 66 sieve. One log <= 2^126.686; **two logs <= 2^113.398 units**;
meter 2^115.0 (Markov <= 0.3295). P(boot fails on G) <= 4.0e-04 ([7]; 3.41e-04).

## Supporting measurements ([11])

Norm bins (offline, 20 lattices/size): fresh 150/260-bit lattices stay in-bin
(cf max 53.4 <= 109.6, cg max 5.19 <= 20.5); [11] re-samples per run. Cells via
`int_ball_point` (uniform-box rejection, exact; R=1 chi2 5.8, R=4 shell chi2
13.8). Nodes: UDUP = 0.65 is only a normalization ([11s] measures exact-orbit
yield; no node haircut). Boot: 3000 cells of 3 certified lattices recounted
(PARI cross-check): yr 1.019 [1.002,1.034] vs sextile pr, 0.886 [0.871,0.899]
vs charged pr, so E[shat/m] >= 0.87x charged pr > 0.65 ([11u]); [11] repeats a
96-cell check live (threshold 0.70) plus a sampler smoke check.

## Refined integration representation

Mesh 384 (was 48), rho spacing 1/1024 (was 1/128). Historical iA values above are physical windows; executable tuples scale them by eight.

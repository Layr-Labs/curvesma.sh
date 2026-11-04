# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T is replaced by T' = T^e, e uniform in [1, r): uniform on nonidentity mu_r (distance 1/r; identity z handled directly; ledger charges 2/r). The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6. The boot
element has a single (g-side) norm N(a^2 + U y b^2), U = |u| = 2^63.71; the basis is LLL-reduced
in the weighted metric with the b-block scaled by W = isqrt(U) = 2^31.86
(`boot_lll`), which balances the two terms, and returned unscaled. The certificate is exact: membership,
size and Lovasz (0.98) conditions in the weighted metric, Gram determinant
p^12 W^12 ([11], `certify_boot`). Measured on 8 certified lattices x 384 uniform
R = 4 points: g-norm offsets cg mean -191.8, p90 -186.0, p99 -182.0, max -178.6 (Euclidean basis: mean -7.3, max 12.6; spread 0.53). The boot
enumerates the radius-4 ball (2^23.67 points, one per +- pair) and runs the `ecm16.py` schedule for 199.79-bit primes: a relation
has every prime <= 2^215, <= 15 in (2^69, 2^199.79], <= 1 in (2^199.79, 2^215].
An exhausted ball draws a fresh exponent and certified lattice (2^36 each).
Distinct relations credit UDUP = 0.65. It stops at 18 distinct
relations (BOOT_RELS): 2^52.37 trials of 2^67.75 instructions over ~2^28.7
lattices; failure needs all 18 to fail.

## Nodes

A node is a special-q ideal (l bits, side f/g) with lattice L_q (det q,
LLL-reduced). Its region is the radius-R integer ball in the reduced basis,
one point per +- pair; exact half counts are in `stats16.b85` (recomputed for
R <= 8 by [7]). Two node types:

- **sieve**: every p in [2^8, 2^69] and every relative-degree-one root pair
  is enumerated on the ball by `dsieve16` (KERNELS.md); records go to the
  owners' 2^40-entry score arrays. Every primitive cell gets a certified upper
  bound on both norms (72 + 12 linear forms in double precision, outward
  rounding; ~1,500 instructions, charged 2,400). Cells with residual mass
  (log2 bound minus score, minus log2 q on the special-q side) <= T = kA nA +
  kB n + 96 on both sides survive (exact norms, trial division, ECM).
  Non-primitive cells are skipped.
- **ecm**: ECM trials on every ball point, without a sieve (never cheaper in
  the shipped plan).

A relation has every residual prime <= 2^n per side, at most kA in (2^69, 2^nA]
and at most kB <= 1 in (2^nA, 2^n]; these window primes are the children, each
charged at the most expensive bin at or below its window.

## Yield and survivors

Per cell, the two-window semismooth probability

  P = sum_{jA <= kA, jB <= kB} 1/(jA! jB!) int rho((m - sum t + alpha)/69) prod dt/t

is a lower sum: each mesh cell carries exact weight int dt/t and rho at its
left endpoint (rho decreases); cells crossing sum t = m are dropped. Norm
sizes m come from measured offset bins over 61,440 exact ball norms on 160
production special-q lattices (10 sizes 72-300 bits, both sides, 6 radii): cf
sextile maxima paired with the overall cg maximum, shifted by worst sampled
per-(ideal, radius) mean minus pooled mean (cf 20.68 + 1.58 -> +23, cg -2.53 +
7.81 -> +6). Boot norms use 3,072 exact norms on 8 certified skewed boot lattices at R = 4 (p90/p99/max + 4 bins). The semismooth
probability integrates over 24 norm sub-bins at nodes and 10 at the boot (`QBINS_SUB`/`BOOT_SUB`),
charged at most 4x / 1.5x the sextile-bin value (UDUP support below is stated against this charged pr). mu = cells x P x UDUP x detection x retention (UDUP
= 0.65; ECM detection 1 - 0.001 x expected factors; the boot carries only the counted fraction). Survivors per cell use
(T/69) rho((m - T)/69) per side, decreasing in m, evaluated at each bin's
lowest norms (previous sextile's raw cf max, observed cg min -54.56, margins
(23, 6); special-q at the bin's upper edge).

## Factor-base ideals without a logarithm

Under the main-phase occurrence model (`filtsim.py`: 2^63.589 relations -- the
conservative descent input, the excess-1.1 theta = 2347/8000 retained-count
estimate 2^63.689099 minus 2 pool SE (0.0499 each; [7s] checks the floor),
certified lower product 0.968032661 already applied), an ideal of norm q = 2^(69 t) on side s occurs
with probability 2 kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/rho(u_s), u_s from
mean norms 738.2/498.4 (rounded down from 738.26/498.49; [4] checks them).
Occurrences are Poisson(lambda_s(b)); an ideal occurring <= 2x is missing. A
relation of an m-bit side carries on average sum_b kappa(m/69, b/69)
ln((b+0.1)/b) P(Poisson <= 2) missing ideals (1/10-bit bands 8 <= b < 69, the
full factor base; v22r1 closes the b < 40 tail: lambda >= 2^23 there, so those
terms are exactly 0.0 in binary64 ([6c] checks the tail sum is 0.0 at the shipped
input and at 63.5, below every [7s] input):
0.20-1.84 per f side, 0.40-3.68 per g side over u = 7..39 (3.74 at the boot's
2,709-bit norm). The proxy applies with coefficient 1, so every band charges at least the e^-lambda
zero-occurrence mass ([6b-zero] checks it on a lambda grid). Each node takes one option:

- **discard**: relations with a missing ideal are rejected (mu x exp(-missing)); the 30 thin nodes discard, so missing
  ideals never recurse there: all 8 with q < 2^90 (THIN_MAX 85 forces thin) plus the 22 plan-flagged thin=1 nodes
  (90f/g, 95g, 100f/g, 105f/g, 110f/g, 115f, 120f/g, 130f/g, 135f/g, 145f/g, 170f, 175f, 205f, 220f; exact list in
  refinement.json `declared_thin`, audit-bound to the `node_cost` flags);
- **descend**: each missing ideal Q is a further child at the other 32 nodes, sieved in the lattice of
  Q q' and charged the largest work and failure probability
  of the discarding nodes with q <= 2^90;
  the boot descends its missing ideals the same way.

Cascading removal can also lose ideals with >= 3 occurrences (e.g. an ideal
in exactly 3 rows each holding a different singleton is purged entirely); [6b] measures losses on the [6] system after purge + clique + merge (weight 200) + final-matrix-seeded reconstruction (pivot back-substitution + original-row fixpoint, soundness below): per lambda band
(ideals, measured, model, occ>=3 R1-form/other/merge): [1,2): 897, 0.0970, 0.7938, 0/2/0; [2,3): 961, 0.0822, 0.5602, 1/1/0; [3,5): 757, 0.0343, 0.2761, 0/0/0; [5,8): 454, 0.0044, 0.0628, 1/0/0; [8,1e+09): 667, 0.0000, 0.0016, 0/0/0;
the model exceeds measurement by >= 6.81x in every band with loss (merge loses nothing here).

**Reconstruction soundness (v22r2).** `missing_recon` seeds known logs from
the FINAL matrix only (v22r1 wrongly seeded from original rows), then applies
elimination pivots in reverse order and original rows to fixpoint. Pivot rows
are true relations (elimination preserves row space), and triangularity holds:
a column eliminated at t is in no alive row after t (deleted from all of them;
later pivots cannot reintroduce it, singleton-deletion only removes rows), so
each pivot row carries at most one unknown in reverse order (complete, no
circularity; machine-checked on every campaign system). Original rows are true
relations (one-unknown rule sound); deleted rows are on-disk truth, usable.

## Table-good event (v22 joint control)

S_j(M) = missing mass of node j's relations given the realized post-filtering
table M (63 functions); S_j = S_occ,j + S_res,j (occ<=2, pre-filtering; residual
= filtering-lost occ>=3). Descent runs on G = G_occ AND G_res AND G_hi:
- G_occ (per-node McDiarmid + union): max_j (S_occ,j - 1.05 E_j) <= 0.01. Each
  S_occ,j has per-draw difference <= 2*WMAX*wmax (WMAX = 1500, pool bits <=
  1400 live (widths beyond the pool max are the H5 width premise, disclosed);
  wmax = 400/R covers kappa ratios <= 4 live; lambda <= 100). Per-node
  McDiarmid: P(S_j - 1.05E_j > 0.01) <= exp(-2(.01)^2/sum_c2) (the 0.05E_j
  slack only helps); union over 63 nodes+boot: [6c] -1926.9 + log63 <= ln 1e-9,
  charged 1e-9. Rigorous under H5.
- G_res (one Markov event): boot-weighted residual C(M) <= 0.14. Boot weights
  dominate (kappa ratio <= 0.865, [6c]), so all S_res,j <= C(M).
  E[C] = 0.03232514 (Wilson-95 campaign rates x boot masses, full 8 <= b < 69 range);
  P(G_res^c) <= 0.230894, charged 0.24. Rate transfer to production is H5-res (disclosed).
  H5-res calibration: pooled residual rates over 2^12..2^17 are
  0.00318/0.00205/0.00230/0.00234/0.00236/0.00231 (flat over 32x; pooled counts
  in refinement.json, per-band table in the note); each band's charged upper covers its 2^16 and 2^17 point
  rates ([6c] `ok_cover`). Tree obstructions (ideal in 3 singleton-tainted rows)
  are scale-free given occurrences, hence consistent with non-rise and included
  in the measured residual; vanishing is never claimed.
- G_hi (union): no lambda > 100 ideal occ<=2 (<= 2^(b+6) ideals/band x
  e^-lambda(1+lambda+lambda^2/2); [6c]: 8.78e-20, charged 1e-9).
P(G^c) <= 0.24 + 2e-9. On G, miss = 1.05 E_occ + 0.15; boot failure <= psi^18
(fresh exponents given M); E[W] <= W_infl + P(G^c) cap. Success subtracts the
delta once (proof 7).

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

Distinct relations per node are negative binomial (mean mu, shape r; H5 premise, per-ball form — proof section 9): r = 0.5
for effective windows of <= 1 prime (base nodes have empty
windows), r = 1 otherwise. Shapes come from toy-ball fits (empty-window rows
in `stats16.b85`, printed by [11b]): one-prime windows fit 0.60, two/three-prime
1.86-12.2. Of 12 empty-window configs, all with mean >= 7.4 show 0-4 ideals
without a relation vs 5.4-12.4 predicted at r = 0.5; the three with means
3.2-4.0 show 19/11/21 vs 18.1/17.2/16.4, and raising lb clears them (0 vs
7.3-11.2). Base nodes run at mu = 613 and 470.

A child subtree fails with probability phi; phi(node) = max((1 + mu (1 -
psi)/r)^-r, 0.02), psi = window sum of children's phi; a node tries <=
1/(1 - psi) relations on average. Node work = own sieve/ECM work + 1/(1 - psi)
x children's work. The plan requires mu >= 18, phi <= 0.05, psi <= 0.8
everywhere ([7] checks). The descent keeps a node counter over both
logs (cap N_MAX = 2^27.3, counted failure); E[N] <= 2^17.57
(`counts16.py`, unrounded 2^17.5645), so P(abort) <= 2^-9.73 (Markov). All N_MAX nodes in the
smallest bin (2^70.3 ideals) with 2^6 for 1/p weighting: shared-ideal
probability <= 2^-10.70, charged in full ([7]).
A factor-base ideal Q without a log
(norm <= 2^69) is never descended in its own fixed lattice: each visit draws a
fresh ideal q' of about 2^20 with known log and
sieves the lattice of Q q'. Its node is charged as the worst node with q <= 2^90, and all nodes up to that
size discard relations touching an ideal without a log.

The shipped plan has 62 nodes (bins (69, 75) plus 5-bit bins 75-225, both
sides) and the boot (n1 = 215, nA = 199.79, kA = 15, R = 4).
The 90d99567 plan was re-optimised by single-then-joint coordinate
search at the sub-bin pr, the skewed boot and LOG2_RELS = 63.599:
41 of 62 entries changed (refinement.json `plan_v21`); the boot window is
unchanged. This revision re-evaluates that plan unchanged at the new inputs.
At the inputs (LOG2_RELS = 63.589, mean norms 738.2/498.4,
UDUP = 0.65) with THIN_MAX 85, evaluated on the table-good event: min mu 35.0,
max phi 0.0484, max node psi 0.509, boot psi 0.6717. All
62 nodes sieve. One log <= 2^126.343 instructions;
**two logs <= 2^113.055 units**. Meter stops the descent at 2^115.0 units
(Markov: E[W]/cap <= 0.26). P(boot fails on G) <= 1.000e-03 ([7]; modelled
7.74e-04). Spot cells at re-optimised nodes ([11] method,
40 cells, units in [-2, 2]^6; spot checks, not the campaign bound):
E[shat/m]/pr_charged 69f 7.65, 90f 6.18, 105g 9.7, 125f 8.17, 140g 3.51, 180g 3.33, all >> 0.65.

## Supporting measurements ([11])

Norm bins (offline, 20 lattices/size): fresh lattices at 150/260 bits stay in
the node bins (cf max 53.4 <= 109.6, cg max 5.19 <= 20.5); [11] re-samples per run. Unit duplicates (offline, relation-weighted, uniform cells): distinct
yield CELLS x pr x 0.65 needs E[shat/m] >= pr x 0.65 per regime (shat =
model per-point semismooth prob, m = mates among
1532 units). Cells are drawn by `int_ball_point` (rejection from the
uniform box, exact by construction;
validated: R=1 chi2 5.8, R=4 shell chi2 13.8; the old
round-and-reject sampler scored 506.4 and is retired).
yr = E[shat/m]/pr over 13 node regimes (bits 69-225, both sides, R
192-384; 750 pts x 1-3 lattices): 5.18-150, min CI lower 4.70 against the
sextile pr, hence >= 1.175 against the charged pr (4x cap). R=4 skewed boot
over 3 certified lattices (1000 pts each): 0.99 [0.95, 1.03] / 1.02 [1.00, 1.06] / 1.05 [1.02, 1.08], pooled 1.02 [1.00, 1.04]
(3000 points; against the charged pr, boost 1.16x: 0.89 [0.87, 0.91]; saturation 100% / 100% / 100%).
Complete orbit ([11u], v22r3): per-cell height bound h(v) (Cauchy-Schwarz + Emb^-1, 2% margin) certifies which units can map v into the ball; all 2442 units of height <= 8 enumerated (gen_units); cells with h<=8 get the exact orbit count, others contribute 0 (lower bound, no independence premise).
All 3000 skewed-boot campaign cells recounted (PARI bnfcertify cross-check):
exact yr 1.019 [1.002, 1.034] against the sextile pr (0.886 [0.871, 0.899] against the
charged pr), so E[shat/m] >= 0.87x (boot) the charged pr, above UDUP 0.65.
Node campaign (12 regimes x 100 cells, shipped plan R, pre-registered seeds 60001-60012): certified yr lower 95% CI minima 3.20 (220g) to 12.77 (69g), 0 diffs among 949 certified cells (box-5 == exact-8 throughout); E[shat/m] >= 3.20x (nodes)
the charged pr, above 0.65 with 4.9x margin. [11] repeats exact-8 counts live on fresh seeds
(thresholds 1.175 / 0.70) plus a sampler smoke check.
Heterogeneity ([11h], v22r3): per-lattice mu via 50-cell shat means over 100 fresh lattices at min-mu (80g), max-phi (210g), max-psi (220g), mid (135f): means 40-58x plan mu, 0/100 below plan per regime (Wilson 95% upper 3.70%), per-lattice phi all at the 0.02 floor (Jensen ratio 1.00). Bad-lattice mass <= 0.04 charged in success (HET_TABLE); transfer to selected ideals disclosed in H5 with the 40x conservatism margin.

## Refined integration representation

Mesh 384 (was 48), rho spacing 1/1024 (was 1/128). Historical iA values above are physical windows; executable tuples scale them by eight.

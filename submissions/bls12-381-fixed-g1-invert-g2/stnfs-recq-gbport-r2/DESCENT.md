# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Targets and boot

Each target T is replaced by T' = T^e, e uniform in [1, r), log T = log T'/e
mod r: uniform on nonidentity mu_r (distance 1/r from uniform on all mu_r;
identity z handled directly; ledger charges 2/r). The boot lattice
{(a, b) in Z^12 : a B + b A = 0 in F_p[y]/h} has determinant p^6; its LLL basis is
certified exactly (Gram determinant p^12, Lovasz conditions; [11]). The boot
element has a single (g-side) norm. The boot enumerates the radius-4 ball
(2^23.67 points, one of each +- pair), computes each norm and runs the
`ecm16.py` schedule for 199.79-bit primes: a relation has every prime <= 2^215,
<= 15 in (2^69, 2^199.79], <= 1 in (2^199.79, 2^215]. An exhausted ball draws a
fresh exponent and certified lattice (2^36 each). Distinct relations credit
UDUP = 0.65, as at nodes. It stops at 18 distinct relations: 2^58.47 trials of
2^67.96 instructions over ~2^34.8 lattices, tried in turn; failure needs all
18 to fail.

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
  Non-primitive cells repeat c / gcd(c)'s relation and are skipped.
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
7.81 -> +6). Boot norms use 3,072 exact norms on 8 certified boot lattices
(p90 + 4, p99 + 4 bins). For pr the sextiles are subdivided (`QBINS_SUB`:
24 cf-quantile cells, per-cell cg maxima (fresh uniform grid), same +23/+6
(top edge pinned); `BOOT_SUB`: 10 R=4 deciles +4 (edge pinned); charged
boost capped at 4x nodes / 1.5x boot. Sub-cell maxima dominate their cells (lower bound kept); fresh
gaps 19.14/3.78 re-validate the margins. mu = cells x P x UDUP x detection x retention (UDUP
= 0.65; ECM detection 1 - 0.001 x expected factors; retention from the next
section; the boot carries only the counted fraction). Survivors per cell use
(T/69) rho((m - T)/69) per side, decreasing in m, evaluated at each bin's
lowest norms (previous sextile's raw cf max, observed cg min -54.56, margins
(23, 6); special-q at the bin's upper edge).

## Factor-base ideals without a logarithm

Occurrence model (`filtsim.py`) at a **conservative relation count**:
LOG2_RELS = 63.883 = the theta = 9/25 estimate 2^64.0641 minus two pool SE
(0.0903 bits), floored, certified lower product applied ([4] asserts it). It
drives every missing rate, the work, the failure terms and success; run
`extra = sensitivity` ([7s]) re-evaluates the plan at estimate - 0/1 SE (two-log
2^113.6176 / 2^113.6696, success 0.4759 / 0.4562; 2^113.7548 / 0.4225 at - 2 SE).
An orbit of norm 2^(69 t) on side s occurs Poisson(lambda) times, lambda =
2^LOG2_RELS 2 kappa_s(t) 2^-69t, kappa_s(t) = rho(u_s - t)/rho(u_s) at the
relation-weighted mean norms 743.2/498.9 (rounded down; [4]); an orbit occurring
<= 2 times counts as missing. A side of size u meets on average mu(u) = sum over
1/10-bit bands [b, b+.1), b in [40, 69), of kappa(u, (b+.1)/69) ln((b+.1)/b)
P(Poisson(lambda_lo) <= 2), lambda_lo = 2^(LOG2_RELS-b-.1) 2 kappa_s(b/69): an
UPPER sum (kappa rises in t, P(Poisson <= 2) falls in lambda). The descent
charges MISS_SCALE mu with MISS_SCALE = 1.05 ([7m]): 0.47-1.07 per f side,
0.39-0.66 per g side at the plan's nodes, 3.405 at the boot. Each node takes
one option:

- **discard**: relations containing a missing orbit are rejected (mu times
  exp(-missing)); nodes with q <= 2^90 always discard (forced up to bin 85;
  bin 90 keeps its searched flag), so missing orbits never recurse there;
- **descend**: each missing orbit Q (norm <= 2^69) is a child descended as an
  **ordinary prime special-q node in its own lattice L_Q** (det N(Q); no
  composite lattices), charged the largest work and failure of the discarding
  nodes with q <= 2^90: node norms (72c + cf - l and 12c + LOG2U6 + cg, c =
  (l + S)/12, f-side special-q; g alike) do not decrease in l, so L_Q's are
  no larger than bin (69, 75]'s and its relations discard (recursion closes);
  [11] samples a 66-bit L_Q inside the node bins (cf max 19.4 <= 109.6, cg max
  2.8 <= 20.5). The boot descends its missing orbits alike; a second visit of
  Q re-enumerates the same ball and is charged as failure ([7r]).

## Retention and cofactorisation

**Counted cells.** A cell is counted only if base-primitive (for every p and
irreducible h of H mod p, a and b not both divisible by (p, h)); other cells
repeat (a - b theta)/pi's relation up to p's free relation. Base-primitive
fraction prod (1 - p^-2 deg h) = 0.9100 (exact degrees for p < 2^9). Cells
divisible by a relative-degree-one ideal above a base prime of degree >= 2
with p >= 2^8 are excluded by ordinary fractions 0.9523/0.9919 (f/g, disjoint).
The factor base holds every prime ideal of norm <= 2^69 above p < 2^8 plus the
degree-one ideals above p in [2^8, 2^69].

**Root enclosures and the bound gap.** The 84 linear forms L_k(v) = a(eta_k) -
xi_k b(eta_k) use certified roots (`descent_witness.forms`, [11b]): eta from
`embedding_bounds.root_disks()` (rational centres, radius 2^-100, Rouche); for
each eta the 12 roots xi of f(eta, x) refined to 60 digits and enclosed by
Weierstrass disks D(z_i, 12 |W_i|), W_i = (|f(z_i)| + delta(z_i)) / prod_{j != i}
|z_i - z_j|, delta bounding f's coefficient change over the eta disk, all
pairwise disjoint; g side xi = +-sqrt(-u eta), error |U| r_eta/|xi|. Largest
relative radius 1.8e-28. With T_k(v) = sum |v_j||eta|^j + |xi| sum
|v_6+j||eta|^j, the root error (<= 6.001 x 2^-53 T_k) and double evaluation
error (<= 12 x 2^-53 T_k) lie inside the charged e_k = 2^-46 T_k, so prod (|L~_k|
+ e_k) >= |N(v)| for every cell ([11b] checks it with zero tolerance). The gap
log2((|L~_k| + e_k)/|L_k|) <= log2(1 + 2^-45 T_k/|L_k|) is <= 0.81 bits (f; 0.14
g) inside the charged 1 bit whenever every |L_k|/T_k >= 2^-38; other cells are
**exceptional** and charged as lost. On 2,000 uniform cells of five production
special-q lattices (75-220 bits), P(min |L|/T < eps)/eps <= 8.4 (f) and 5.0 (g)
for eps in [1e-4, 0.1], so the exceptional fraction is <= 10 x 2^-38 ~ 2^-34.7
per side (linear extrapolation); every sieve node's relation count is charged 1 -
2^-20 (`EXC_CELLS`), a 2^14 margin (H5).

**Omitted mass.** A counted cell is rejected only when one side's omitted
factor-base mass X exceeds s = 96 - 1 (bound gap) - m/128 (each of the <= m/8
sieved hits loses < 1/16 bit to floor(16 log2 p)). X sums independent local
contributions (whole valuation at p < 2^8, excess powers at p >= 2^8). For a
base prime P = (p, h_i) with residue field F_q and side polynomial S, a
base-primitive (a : b) has b a unit with probability q/(q+1) (else a is) and
then v_P(S_hom(a, b)) = v_P(S(t)), t = a/b uniform in O_P. A root of S of
multiplicity mu in F_q lifts to mu roots (Hensel), so v_P >= j forces one
distance >= j/mu: P(v_P >= j and t reduces to it) <= min(1/q, mu
q^-ceil(j/mu)), equality q^-j for mu = 1; primes above p are independent (CRT).
Degrees, ramification and root multiplicities for p < 2^9 are tabulated
(`retain16.structure`, `stats16.b85`; [11b] recomputes them at 2, 19, 23 and a
seventh of the others and checks residue samples at 2 and 23): all
multiplicities are 1 except g above 2 (2^16 | U) and f at the degree-4 prime
above 19; 23 ramifies (e = 2). Above 2^9 excess powers are bounded per dyadic
block (72 p^-2 events per level, pi(x) <= 1.26 x/ln x). Masses rounded up to
1/4 bit: means 3.35 (f) and 7.66 (g) bits, P(X > 96) = 1.2e-12 (f) and 2.3e-10
(g) per cell.

**Relation-conditioned loss.** Under H1 a cell with omitted mass x is
semismooth with probability Pss(m - x), so the per-side lost fraction is

  loss = sum_{x > s} P(X = x) Pss_hi(m - x) / sum_{x <= s} P(X = x) Pss_lo(m - x),

bracketed by monotonicity on a 64-bit grid in x at the node's largest norms
m_f, m_g. Retention = (1 - loss_f - loss_g) x base-primitive x (1 - ordinary
exclusions); [7] evaluates every sieve node (min 0.8593; largest one-side loss
< 1e-5).

**ECM schedule (`ecm16.py`).** Suyama-form Montgomery curves, cheapest (B1,
B2/B1) grid point with wheel 30030; stage-1 x-only ladder (11 mults/bit),
stage-2 cross products (3/pair), one binary gcd per stage. Per-curve success
P1 is the Dickman lower sum on the curve order; curves = ceil(ln 0.001 /
ln(1 - P1)) at 367.76 (L/5)^2 per mult. [12] prints node/boot schedules (boot:
B1 = 2^28, 1418167 curves, 2^67.96/trial) and matches P1 on random 36/40-bit
primes.

## Recurrence and success

Distinct relations per node are negative binomial (mean mu, shape r): r = 0.5
for effective windows of <= 1 prime (base nodes (bin (69, 75)) have empty
windows), r = 1 otherwise. Shapes come from toy-ball fits (empty-window rows
in `stats16.b85`, printed by [11b]): one-prime windows fit 0.60, two/three-prime
1.86-12.2. Of 12 empty-window configs, all with mean >= 7.4 show 0-4 ideals
without a relation vs 5.4-12.4 predicted at r = 0.5; the three with means
3.2-4.0 show 19/11/21 vs 18.1/17.2/16.4, and raising lb clears them (0 vs
7.3-11.2). Base nodes run at mu = 1106 (f) and 935 (g).

A child subtree fails with probability phi; phi(node) = max((1 + mu (1 -
psi)/r)^-r, 0.02), psi = window sum of children's phi; a node tries <=
1/(1 - psi) relations on average. Node work = own sieve/ECM work + 1/(1 - psi)
x children's work. The plan requires mu >= 18, phi <= 0.05, psi <= 0.8
everywhere ([7] checks). The descent keeps a special-q node counter over both
logs and aborts past N_MAX = 2^24.75 nodes (counted failure); E[N] <= 2^16.61
(`counts16.py`, unrounded 2^16.6007), so P(abort) <= 2^-8.14 (Markov). All N_MAX
nodes in the smallest bin (2^70.3 ideals) with 2^6 for 1/p weighting:
shared-special-q probability <= 2^-15.80, charged in full ([7]).

**[7m] The shared missing set (H5m).** All descent relations meet ONE
main-phase set M, so the averaged rates need a configuration bound. Premise
H5m: main-phase occurrence counts of distinct orbits are negatively associated
(NA) with the Poisson marginals above (independent relations, the model behind
those marginals; within a relation the orbits above one prime are mutually
exclusive and distinct primes CRT-independent before the smoothness
conditioning). The indicators of M3 = {orbits occurring <= 2 times},
non-increasing functions of disjoint counts, are then NA (Joag-Dev-Proschan
1983). M lies in M3 u C, C
the filtering-cascade losses (>= 3 occurrences, no log). [6c] (two systems per
run; seeds 1, 2: 6 losses, rC <= 0.0076) plus 54 offline filtsim systems (2^12-2^14, 32/16/6 seeds, production occurrence profile,
220,641 orbits with >= 3 occurrences): losses are 1.1e-3 of them (3-6
occurrences), their descent-weighted mass <= 0.0077 of the model mean
(pair-weighted <= 0.0131, cross-profile <= 0.0128, node x boot 0.0071), flat in scale; the mass actually left without a log
is <= 0.32 of it. Charged EPS_C = 0.04, EPS_CP = 0.05 (>= 3.8x the largest
measured). A relation class k (62 nodes, boot; per-orbit weight w(u, b) = 2
kappa(u, (b+.1)/69) 2^-b) meets M at S_k(M) <= S_k(M3) + S_k(C), E S_k(M3) <=
mu_k, and Chernoff-Hoeffding for NA [0, a_k] summands (Dubhashi-Ranjan 1998) gives P(S_k(M3) > 1.01
mu_k) <= exp(-(mu_k/a_k) h(0.01)) < 2^-(4.7e12), a_k the largest weight above
2^63; orbits below 2^63 lie in M3 with total probability <= 2^-27.31 (union
bound, N(b) <= NBAR 1.01 orbits per prime and side). On the complement G every
relation's rate is <= MISS_SCALE mu_k = (1 + 0.01 + 0.04) mu_k; the recurrence
is monotone in the rates, so its value at MISS_SCALE bounds E[W | G] and the
failures, and success subtracts P(not G) <= 2^-27.31. Given M, distinct
subtrees sieve distinct prime lattices with fresh coins (repeats charged
below) and are treated as independent (H5).

**[7r] Missing-log repeats.** Tried relations r1, r2 share a missing orbit with
conditional probability <= C_12(M) = sum_{I in M} w_I(r1) w_I(r2). Both meet the
SAME set, so E C_12(M3) = sum_I w_I(r1) w_I(r2) m_I (one factor m per orbit),
bounded per band by min(w(u1, b) mu_b(u2), w(u2, b) mu_b(u1)) at the largest
descending-node sides and the boot: per pair 2^-59.87 (node-node), 2^-60.09
(node-boot), 2^-57.69 (boot-boot). C_12(M3) <= 1.01 E C_12 except with
probability < 2^-(1.5e11), and C adds <= 0.05 E C_12. Each tried relation
succeeds (one per node) or fails through a failed child, so <= 2 N_MAX node and
36 boot relations are tried: P(a missing orbit is visited twice) <= 2^-9.28,
charged as failure. N_MAX = 2^24.75 balances this against the counter term.

The shipped plan has 62 nodes (bins (69, 75) plus 5-bit bins 75-225, both
sides) and the boot (n1 = 215, nA = 199.79, kA = 15, R = 4), chosen by
coordinate search over (n, iA, kA, kB, R, discard) minimising modelled two-log
work subject to the three conditions (our 57598408 lineage's re-searched
plan). At LOG2_RELS = 63.883, mean norms 743.2/498.9, UDUP = 0.65, MISS_SCALE
1.05, THIN_MAX 85: min mu 38.5, max phi < 0.05, max psi 0.521 (boot); all 62
nodes sieve. One log <= 2^127.043 instructions; **two logs <= 2^113.7548 units**
(unrounded 2^113.75477). Meter 2^114.56 (Markov: E[W]/cap <= 0.572271).
P(boot fails) <= 7.8367e-06 ([7]).

## Supporting measurements ([11])

Norm bins (offline, 20 lattices/size): fresh lattices at 150/260 bits stay in
the node bins (cf max 53.4 <= 109.6, cg max 5.19 <= 20.5), the certified boot
lattice in the boot bins (cg max 6.6 <= 16.7); [11] re-samples points per run.
Unit duplicates (offline, relation-weighted): distinct yield CELLS x pr x 0.65
needs E[shat/m] >= 0.65 pr per regime (shat = model semismooth probability at
exact cofactor norms, m = mates among 1532 units), on cells drawn by
`int_ball_point` (box rejection: exactly uniform; shell chi2 5.8 at R=1, 13.8 at
R=4; the old round-and-reject sampler scored 506.4 and is retired). The
campaign measured yr0 = E[shat/m]/pr0 against the UNINTEGRATED sextile pr0
(stats16.b85): 13 node regimes (bits 69-225, both sides, R 192-384; 750 pts x
1-3 lattices) 5.18-150, min CI lower 4.70; R=4 boot (3 true lattices, 1500
pts) 1.10-1.25, CI lower 1.06. The model charges pr <= cap x pr0, so against the
CHARGED pr only yr >= yr0/cap is claimed: nodes >= 4.70/4 = 1.175, boot >=
1.06/1.5 = 0.707, both above 0.65. Truncation 86-100% saturated, x5
(<=11.5%)/bound-8 closed. Run `extra = relw` repeats both live (6 + 6 points:
node 1.30, boot 0.86 at seed 1, thresholds 1.0 / 0.72) with the R = 4 shell chi^2.

Quadrature: mesh 384, rho spacing 1/1024; plan iA tuples are 8x the physical windows.

# Frobenius-class rho in GT with lazily evaluated steps: a prefix-product trie, a 4,964-instruction key kernel and one full multiplication per six steps

**Claim.** For uniform nonidentity P in G1 and uniform z in GT the algorithm below recovers Q in G2 with
e(P,Q) = z with probability at least 0.39, using at most **2^124.210** dense-Fp12-multiplication equivalents
(C = 20,000 instructions), including setup (the multiplier table and the prefix-product trie), the two pairings,
walk starts and the fresh start of every trail, distinguished-point (DP) handling, fruitless-cycle escapes, the
DP detection delay, matching (sort, duplicate removal, replay of matched trails, pair checks), output and
verification. Peak memory is below 2^62 bytes; preprocessing (included) is at most 2^49.06 units; nonuniform
advice 0 bytes. The success bound rests on three disclosed heuristics, none of them proven: the usual
random-mapping heuristic for rho walks (H1), pseudorandom bits of the class key (H2), and whole-run concentration of
the realized waste of the run (H1+, stated as a tolerance of twice the H2 rate, with its success loss computed,
section 3.4). The success probability is
computed in the ideal model: the ChaCha20 stream that supplies the table and start scalars is an ideal random
source (also disclosed; no bound on the realized scalars is claimed). No production-size run is claimed.

## 0. What is new, and what is reused

**This package is a port to bls12-381-fixed-g1-invert-g2 of `gt-lazy-trie-rho`** (submission a77744bb,
mpjunior92, claim 2^124.183 on the mirror track bls12-381-fixed-g2-invert-g1). The tracks differ only in which
pairing argument is the instance and which group the answer lives in. The port changes section 1 (g = e(P, Q0), answer Q = [x]Q0 in G2, Q = infinity exactly for
z = 1, and the argument that (g, z) is distributed as in the mirror track), the output allowance of section 5
(a G2 scalar multiplication, at most 2^11 units; the ledger total is unchanged to the printed digits) and the
claim text. The walk, prices and ledger (the rest of sections 2 to 7) are a77744bb's text and code, re-run here
with the same declared configuration; the probability argument (sections 3.3, 3.4, 4, 6, proof-impl.md 0-iii) is
rewritten, for the reason at the end of this section. The score rests on this package's own analysis, sections 1 to 5, not on any score or acceptance on the
mirror track; a77744bb is provenance only.
**Changes from the rejected revisions.** 680b01db, three findings: (i) H1-lower (success law, N') and
H1-upper (union bounds, pairwise 2/N) are named and `ledger.py` divides by N in the upper direction (sections 3.5, 4, 6);
(ii) the escape compares masked keys (bits 64..149 out, `esc_key`), a full-key tie-break is a charged loss (sections 4, 5);
(iii) `experiment.py` prints one JSON line per part and the 12,000-step R = 2 walk is offline in `evidence.md`. 1b06bb8b,
one finding: the doubling checkpoint's countdown and interval were described as one word, which wraps at 2^63 while a
segment may run to A = 120 x 2^86 steps; now two words each, started at 1 (proof-impl.md 0-i), priced at 12 of the 174
instructions per step (section 5, +4) and run across the 2^64 boundary by `flags_check.py`.

Ledger: log2 W = 124.20923 (`experiment.py` part D), claim **124.210**, the least thousandth above it, checked by
an exact integer inequality.

Lineage (what a77744bb took from its parent 0b8665ae, how it revised the rejected d46d2b97, and what the earlier ports 9e5c78c3 and 819debd1 lacked) is in this entry's public submission note under "Lineage"; every number quoted there recurs in sections 0a to 7 below.

## 0a. Representation: which limbs the rule reads

As in the parent, an Fp value is held as the six-limb Montgomery vector x~ = x 2^384 mod p, fully reduced in [0, p),
and the rule is stated on these stored limbs: the twelve keys are the stored residues of section F2 (linear
operations commute with the Montgomery map, so the stored keys are the Montgomery images of the integer keys,
again a Frobenius-invariant multiset), K is their minimum as 381-bit integers, kappa the first index (in the
kernel's fixed order) attaining it. The rule on stored limbs is the integer rule composed with a fixed bijection of Z/p,
so every structural property used below carries over (parent section 0a). `lazy.s` and `root.s`
return exactly the stored a, b, B, a zeta, B zeta, K of the current element (parts K and C check this against a big-integer
reference, including keys from explicit Frobenius maps).

## 1. Reduction

Compute g = e(P, Q0) with the pinned convention pairing(Q0, P), where Q0 is the fixed G2 generator. P has prime
order r and the pairing is nondegenerate, so g != 1 generates GT (prime order r) and z = g^x for a unique x in
[0, r). Find x, output Q = [x]Q0 and check e(P, Q) = z; by bilinearity e(P, [x]Q0) = g^x = z. If z = 1 output
Q = infinity (the only valid answer then). The answer is for the supplied P.

Everything after this step sees only the pair (g, z). For uniform nonidentity P, g = e(P, Q0) is a uniform
generator of GT (P -> e(P, Q0) is a bijection from G1 \ {O} onto the generators), and z is uniform and independent
of P, exactly the distribution of (g, z) in the mirror track bls12-381-fixed-g2-invert-g1 (g = e(P0, Q) for uniform
nonidentity Q). The walk, its success analysis and its ledger therefore apply unchanged; only the output scalar
multiplication moves from G1 to G2 (its price is bounded in section 5).

## 2. Arithmetic facts (checked in `experiment.py`, parts A, K, C)

Tower: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi), xi = 1+u (blst's tower; with u = w^6 - 1 it is the
profile field). y = sum_{m<6} y_m w^m with y_m in Fp2.

**F1 Frobenius classes.** p has order 12 modulo r, p^6 = -1 and 1 + p^4 + p^8 = 0 (mod r), so y -> y^{p^k} is a
group of order 12 acting on GT and containing inversion. Every y != 1 lies in a class of exactly 12 elements;
N = (r-1)/12 classes, log2 N = 251.27213.

**F2 Keys (the parent's).** zeta = xi^{(p^2-1)/6} is a primitive sixth root of unity in Fp, d_1 = xi^{(p-1)/6}.
With c_1 = a + b u the w^1 coefficient of y and B = a Re(d_1) + b Im(d_1), the Re(w^1 coefficient) of y^{p^k} is
a zeta^j (k = 2j) and B zeta^j (k = 2j+1), j = 0..5; as zeta^j runs through 1, zeta, zeta - 1, -1, -zeta, 1 - zeta,
the twelve keys are a, a zeta, a zeta - a, -a, -a zeta, a - a zeta and the same for B. K = min of the keys is a
class invariant; kappa = argmin, and y^{p^kappa} is the same element for all 12 members of the class (the
canonical representative). The keys are distinct unless two coincide by a linear accident (at most 66/p per
element) or a or B is 0 (then a key p - 0 = p is formed; at most 2/p per element); both are charged as failure
(section 4).

**F5 Lazy key (new).** For x, Pi in Fp12, the w^1 coefficient of x * Pi is
`c_1(x Pi) = x_0 Pi_1 + x_1 Pi_0 + xi (x_2 Pi_5 + x_3 Pi_4 + x_4 Pi_3 + x_5 Pi_2) = sum_{k<6} x_k gamma_k`,
gamma_0 = Pi_1, gamma_1 = Pi_0, gamma_k = xi Pi_{7-k} (k = 2..5), because w^6 = xi and only the products
w^i w^j with i + j = 1 or 7 contribute to w^1. So a = Re c_1 and b = Im c_1 of the current element y * Pi_nu are
Fp-bilinear in (y, gamma(Pi_nu)): with the Karatsuba sums
`S0 = sum_k x_k.re g_k.re, S1 = sum_k x_k.im g_k.im, S2 = sum_k (x_k.re + x_k.im)(g_k.re + g_k.im)`
(all exact non-negative integers below 6p^2, 6p^2 and 24p^2), **a = S0 - S1** and **b = S2 - S0 - S1** modulo p. b
is exactly sum_k (x_k.re g_k.im + x_k.im g_k.re), a non-negative integer below 12p^2, and S0 - S1 + 6p^2 lies in
[0, 12p^2); both are below 2^768 and one Montgomery reduction gives r < 12p^2/2^384 + p < 2.2p, inside the q-table
subtraction of `gen.py` (r < 8p). Then B, a zeta, B zeta, the keys and K follow F2.

**F3 Cyclotomic squaring (escapes only).** The Granger-Scott map of the parent (`kernels/sqr.s`, 6,777).

**F4 Orders.** r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 * 2529403 * 52437899 *
254760293^2 (checked); ord_r(2) = (r-1)/16, so 2^t lies in <p> only for t a multiple of about 2^247.3.

## 3. Algorithm

Parameters: **R = 24** table multipliers, **trie depth L = 6**, DP density **theta = 2^-86**, M = 2^20 walks
round-robin on one processor (chunk K = 2^63 steps; proof-impl.md 0-i), report cap 2^55, restart cap 2^20,
stuck limit **A = 120/theta** counted steps, matcher budget 2^4 pairs.

**Table.** M_h = g^{a_h} z^{b_h} for h < 24 with fresh uniform a_h, b_h mod r from the sampler (0-iv), and the
twisted entries T[h][t] = M_h^{p^t}, t < 12: 288 entries (each with the 54 Karatsuba pre-sums of `mul.s`), and
the scalar tables a_h p^t, b_h p^t mod r (read only by replays).

**Trie.** A node nu at depth d <= 6 is a word ((h_1,t_1), ..., (h_d,t_d)) of table entries, with
Pi_nu = T[h_1][t_1] * ... * T[h_d][t_d]; the root is the empty word (Pi = 1). Nodes are stored in level order:
the children of node n are 288 n + 1 + 12 h + t. An internal node (depth 1..5) stores its 18-vector gamma(Pi_nu)
(gamma_k.re, gamma_k.im and their sum, k < 6) and the 54-vector `mul.s` entry of Pi_nu (for escapes): 3,456 bytes;
a leaf (depth 6) stores the `mul.s` entry of Pi_nu: 2,592 bytes. Internal nodes: sum_{d=1}^{5} 288^d = 2^40.855;
leaves 288^6 = 2^49.020; 2^60.3661 bytes in all. A node is built from its parent by one `mul.s` call
Pi_child = Pi_parent * T[h][t] and its pre-sums and gamma (16,068 instructions, charged 20,000 = one unit).

**Walk state.** Two base buffers (the current explicit base y with its six Karatsuba pre-sums
x_k.re + x_k.im, and the previous base), the current node and its depth, a ring of three K slots (the current K and
the two before), the previous step's node and buffer index (for escapes), the 64-step checkpoint K and the step
count mod 64, the doubling checkpoint K with its node and buffer index, the two-word doubling interval and the two-word countdown to the next doubling, the
2-word stuck-walk counter and the trail's start pair (alpha0, beta0). It carries no running exponents.

**One iteration** (pre-debit charge first, proof-impl.md 0-i):

```
1  if node = root: (K, kappa) <- root.s(y)            [also writes the six pre-sums of y]
   else:           (K, kappa) <- lazy.s(y, pre-sums, gamma(Pi_node))    [K of the element y * Pi_node]
2  if bits 64..149 of K are 0 (distinguished): F <- bits 150..261 of K; append (F, alpha0, beta0) to the report
   log; the trail ends: fresh start y <- g^alpha0 z^beta0 from a new start pair, node <- root, histories reset
3  2-cycle, 64-step checkpoint and doubling checkpoint tests on K (section 3.3); on a hit: escape, next iteration
4  h <- floor(24 K_0 / 2^64), t <- (12 - kappa) mod 12, node <- child(node, h, t)
5  if depth(node) = 6: y_prev <- y, y <- mul.s(y, entry(Pi_node)) (the flush), node <- root
6  stuck counter += 1; at A = 120/theta: restart the walk (0-iii)
```

A replay runs the same iteration from (alpha0, beta0) and additionally adds (a_h p^t, b_h p^t) to its exponent
pair on every step and doubles it at an escape's squaring, as in the parent; it stops at the trail's distinguished
class and returns p^kappa (alpha, beta) there.

**3.1 The step is a map on classes, and the lazy walk is the eager walk.** The current element is
y * Pi_node exactly (an invariant: it holds at a start, at an escape and after a flush, where node = root, and a
step multiplies Pi_node by T[h][t] on the right). By F5 the key computed in step 1 is the key of that element, so
(K, kappa) and therefore h, t, the DP test and the fingerprint are those of the current element; the next
element is (y Pi_node) T[h][t] = (y^ M_h)^{p^{-kappa}}, y^ the canonical representative, which lies in the class
[y^ M_h]. So the step depends only on the class, and the sequence of elements (not only of classes) is exactly
that of the parent's multiply-only walk with this table. The declared run of part C takes 240 production steps
(47 root, 193 lazy, 32 flushes, one escape, 13 reports, each replayed from its start pair), every key from the
emulated `root.s`/`lazy.s` and every flush and escape from the emulated `mul.s`/`sqr.s`, and asserts at every step
that y * Pi_node equals the eager walk's element and (K, kappa) the reference. Part C's walk carries the
2-cycle test, the 64-step checkpoint and the doubling checkpoint; the extended run recorded in evidence.md (input
r2_steps 12000) adds a second walk of 12,000 steps at R = 2, where long closures are common, which takes 440
two-cycle, 7 sixty-four-step and 63 doubling escapes, replays all 168 of its records and passes the same checks.

**3.2 Collision.** As in the parent: every report is appended to the log until the cap; after the walk phase the
log is heapsorted, duplicates dropped, and the pairs inside runs of equal F are replayed and checked in order until
the first true pair, under a hard budget of 2^4 pairs. Two records with equal F are the same class, up to a false
match with probability 2^-112 per pair of distinct classes under H2. With R <= 2^55 reports the expected number
of false pairs is at most E[R^2]/2 x 2^-112 <= 2^55 E[R]/2 x 2^-112 = 2^-18.27 (E[R] from section 3.5), and the
budget fails only with at least 15 false pairs: Markov gives 2.1e-7. The replays give canonical
pairs with alpha^1 + x beta^1 = alpha^2 + x beta^2 mod r, so x = (alpha^1 - alpha^2)/(beta^2 - beta^1) unless
beta^1 = beta^2 (about 1/r). Output Q = [x]Q0 after checking g^x = z.

**3.3 Fruitless cycles at R = 24.** A closed walk of the class graph is *useless* if the canonical exponent pair
returns unchanged. After l steps the pair is Lambda (alpha^_0, beta^_0) + tau with Lambda = 2^s p^m (s squarings, m a
sum of twists) and tau = sum_h c_h (a_h, b_h), c_h in Z[p]. Because (alpha^_0, beta^_0) is uniform and independent of
the table, the walk is useless (up to probability about 1/r) only if Lambda = 1 and every c_h = 0. An escape's
squaring cannot lie in a useless walk (Lambda = 2^s p^m = 1 needs s >= 2^247, F4). For multiply steps: the l classes
of a closed walk that starts at a first visit are themselves first visits until it closes, so under H2 (section 6) their indices are uniform on 24 and their frames p^e independent uniform on the 12
powers of p (each step's twist is (12 - kappa) mod 12 with kappa a class bit), and c_h is the sum of the frames of
index h's occurrences. The rates below are the model rates for a closed walk starting at a first visit. For coefficients far below r, c_h = 0 mod r exactly when this sum of 12th roots of unity
vanishes in Z[zeta_12]. So the probability P_l that the l steps from a given step form a useless walk is at most
E[prod_h g(n_h)] over the multinomial counts n_h of the 24 indices, with g(0) = 1, g(1) = 0 and g(n) = q(n), the
exact probability that n independent uniform 12th roots of unity sum to 0, for n <= 24 (part F computes q(n) by
exact counts over 12^n for n <= 12 and q(13..24) by convolving the n = 12 law with the n = 1..12 laws:
q(2) = 1/12, q(3) = 1/72, q(4) = 0.0191, ..., q(12) = 0.00195, ..., q(24) = 0.000506), and g(n) = q(24) =
0.000506 for n > 24. The bound holds because a sum of a + b terms has every point probability at most sqrt(q(2a) q(2b)) by
Cauchy-Schwarz (sum_v P_a(v)^2 = q(2a), the roots being closed under negation), so a 24-term sum has point
probabilities at most q(24), and adding independent terms cannot raise one. Part F evaluates the multinomial bound exactly for l <= 600 (a 24-fold polynomial product in floating
point; index counts j >= 403 underflow and take the union bound 24 C(l, j) 24^-j)
and, beyond 600, by Poissonisation (P_l <= e sqrt(l) (E_Pois(l/24) g)^24, monotone on doubling segments),
weighted by the probability exp(-theta l) that a closed walk of length l holds no distinguished class (H2).
Results: **P_2 <= 1/288** (same index twice, opposite twist; the only l = 2 term the bound keeps), sum_{l=3}^{64} P_l =
6.17e-5 (P_3 = 2^-15.34, P_4 = 2^-14.76, P_8 = 2^-25.94, P_64 = 2^-108.88), sum_{l=65}^{600} P_l = 2^-108.79,
the weighted tail beyond 600 at 2^-131.44, so useless walks longer than 64 at a rate **below 2^-108.79 per step**
(charged as 2^-108). Weighted by length, sum_{l > 64} l P_l = 2^-43.86 per step (the far Poisson tail near length 2^86
dominates); section 4 uses it only for the masked-tie item of a doubling escape.

* **Detection and escapes (the parent's rules on K, plus a doubling checkpoint).** The K of the step before last
  is remembered; if it equals the current K, the two classes form a 2-cycle: the walk makes explicit the
  representative of the class with the smaller masked key (K with its DP bits 64..149 cleared, the full K breaking
  a masked tie, so the choice reads none of the bits the detection lemma conditions on) (the current element
  y * Pi_node, or the previous element, Pi of the previous node times its base, y or y_prev if a flush came between;
  one `mul.s` call with that node's entry, none at the root) and replaces it by its square (`sqr.s`),
  node <- root, and resets the histories and both checkpoints; a replay doubles its pair there. A 64-step
  checkpoint K is stored every 64 steps and compared at every step; a repeat means a cycle of length at most 64.
  A *segment* is the run of steps since the latest walk start, restart, fresh start or escape. The doubling
  checkpoint stores the K at segment positions 0, 1, 3, 7, ..., 2^j - 1 (a countdown and interval, both started at 1,
  proof-impl.md 0-i) and is compared at every step (Brent's
  rule): a closure with tail mu and cycle length l inside a segment is caught by step T <= 2 max(mu, l) + l
  <= 3 (mu + l) of the segment (part C checks every rho with mu + l <= 160; largest T/(mu + l) is 2.977). On either
  checkpoint hit the cycle is traversed once more to its class of smallest masked key, which is then squared the same
  way (the ledger charges that re-traversal; `experiment.py` takes the minimum over the segment's stored keys, the
  same classes and target). Every escape target depends only on the cycle, so merged trails stay merged (they may
  escape at different steps, since the 64-step checkpoint counts each trail's own steps, but land on the same class).
* **Stuck limit.** The stuck counter counts every step, escape traversals included, and is reset at every walk
  start, restart and fresh start, not at escapes; a walk whose counter reaches A = 120/theta = 2^92.907 steps
  is restarted (0-iii). The constant is one `movz` (0x1E000000 shifted 16) compared with the high word.
* **Charge.** A 2-cycle wastes the step that re-enters the cycle and the detecting step (2 steps) and one escape
  (24,000 instructions, 1.2 units = 3.4413 steps at c): (2 + 3.4413)/288 = 0.018893 step-equivalents per step. A
  cycle of length 3 <= l <= 64 wastes at most 64 + l steps before its detection, its traversal of l steps at the
  largest step price (19,558, 2.8043 c each) and one escape: sum_l P_l (64 + l + 2.8043 l + 0.6338 + 3.4413) = 0.0050525 (the 0.6338 is the 65 x 64 + 65 x 4 = 4,420 instructions of bookkeeping a 64-step escape pays, section 5).
  A closure longer than 64 (rate 2^-108) is caught within 3 (mu + l) steps and traversed once (l steps, each also paying 68 instructions of
  bookkeeping and masking for the candidate it compares), with mu + l <= A: at most 4A steps at the largest price
  plus 68 and one escape, 2^-108 (4 A x 2.8141 + 3.4413) = 0.00032205.
  The H2 waste rate per first-visit iteration is **W_bar = 0.024268**; the ledger charges
  **eps = 1.25 W_bar = 0.030335** of step units per useful iteration (the toy censuses read pre-collision 2-cycle
  rates of 0.00334 to 0.00346 per step at R = 24, model 0.00347, and checkpoint escapes at 2.1e-5 to 3.5e-5 per
  step against the bound 6.17e-5; evidence.md). No closure is left to the stuck-walk rule.

**3.4 Units, the meter and the work cap.** Three quantities are kept apart and never added in one unit:
*useful iterations* (first visits whose image is not a structural return, the count H1 uses; their number is P),
*waste* (re-entry, detection lag, traversal and escape of fruitless closures, in step units, one unit being c C
instructions; an escape is esc = 24,000/(c C) = 3.4413 units, charged once at its full price), and *restart
steps* (the steps of trails abandoned at the stuck limit, and the unfinished trails at the end). Waste is metered
in work, not counted as iterations. The ledger buys

  W = c S + X_allow + c E + (setup, walk non-step work, stopping allowance) + switches,

with X_allow = eps c S (waste allowance in units), E = (M + 2^20) A = 2^113.907 steps (one stuck window for every
walk and every permitted restart) and the non-step terms of section 5.

The walk phase has the budget B_walk = W_cap C - 2^49.06 C - R_stop C instructions (W_cap = 2^124.210, W_cap C
rounded down) and meters its own charged work in the 192-bit meter under the pre-debit rule (proof-impl.md 0-i)
with the prices of section 5: before each charged action the meter is compared with L = B_walk - A_max
(A_max = 2^25) and the action runs only if the meter is below L. So every run does at most
2^49.06 C + B_walk + R_stop C <= W_cap C instructions, a deterministic inequality `experiment.py` checks on
integers. The step price is fixed per step type and known before the step: root 1,735, lazy 5,138, lazy with
flush 19,558 (the step whose child is a leaf). A doubling escape's traversal charges each traversal step on its
own at 19,558 + 68 (the step, then 64 of bookkeeping and 4 of masking for the candidate it compares).

*No price-drift heuristic.* A block is the run of iterations from a root iteration to the next return to the root
(a flush, an escape or a fresh start). Its iterations are charged 1,735 (the first), 5,138 (the second to fifth) and
19,558 (the sixth, the only one whose child can be a leaf), so a block of k <= 6 iterations is charged at most
1,735 + (k-1) 5,138 < k c C for k < 6 and exactly 6 c C for k = 6, whether or not its last iteration ended in a
report or an escape (the escape and the fresh start are charged on top). Hence the charged step work of n
iterations is at most n c whatever the walk does.

*Meter lemma.* If the realized waste is at most X_allow and the restart steps plus the steps of unfinished
trails are at most E, the meter cannot refuse a charge before the P-th useful iteration with P = S: the work
charged by then is at most walk_need = (c S + X_allow + c E + walk non-step work + switches) C, and part D asserts
walk_need + A_max <= B_walk (slack 5.35e-4 of walk_need). The first hypothesis is H1+ (below, with t = 2 it allows
twice the H2 rate); the second holds whenever the restart cap is not reached, since a trail is abandoned only
after A counted steps and at most 2^20 are, and each of the M walks has at most one unfinished trail of at most A
steps.

*Iteration maximum.* Counting also the action that ends a block, every block of k iterations is charged at least
k c C: a block ended by a flush is charged exactly 6 c C; one ended by an escape at least 1,735 + (k-1) 5,138 + 24,000
>= k c C (k <= 6); one ended by a report at least 1,735 + (k-1) 5,138 + 4,000 plus its fresh start, charged its fixed
allowance of 2^10.1 units (proof-impl.md 0-i; its two 255-bit exponentiations alone exceed 2^21 instructions)
(`experiment.py` asserts all three for every k). So a run has at most n_max = 2^124.210/c + 6M < **2^125.730
iterations** (the last, incomplete block of each walk adds at most 6).

*H1 at first visits, and structural returns.* An iteration is a **first visit** if its class has not been visited
earlier in the run, by any walk. The image of a first visit can land in the trail's own recent classes by the
walk's structure rather than by chance: a 2-cycle returns to the class two steps back with probability 1/288, and a
useless closed walk of length l returns at the rate P_l of section 3.3. The carve-out is defined by position only:
an image that equals a class the walk visited in its current segment (since its latest start, restart, fresh start
or escape) is a **structural return**. It closes a cycle inside one trail, so it can never be the two-trail
collision the run needs; the 2-step check or a checkpoint finds it within 3 (mu + l) steps (64 + l for l <= 64), an
escape leaves it, and its cost is charged in eps at the rates part F computes for these closures. Returns to
classes of the same trail before the latest escape are not carved out: they are within-trail self-collisions,
bounded pairwise under H1 (section 4). H1 (section 6) models everything else: the image of a first visit that is
not a structural return, and the output of a walk start, fresh start, restart or escape, behaves as a class of the
random mapping with the conservative N', so the first S such images contain no class collision with probability
exp(-S (S - 1)/(2N')) (section 4, H1-lower), and any two given such images coincide with probability at most 2/N
(H1-upper: the random-mapping rate 1/N with a factor-2 margin for the non-uniformity of an r-adding walk; part E
reads 1.0068 +- 0.0338 times 1/N on 15,000 classes, 1.206 +- 0.154 on the 3,000 of the declared run; the +- figures are relative standard errors, 1/sqrt(coincidences)). The enlarged N' is used only where a larger N weakens the claim, never to
shrink a union bound. No statement is made about a step conditioned on the history.

*H1+ (a stated tolerance; disclosed, not proven).* The realized waste of the run, in units, is at most t = 2
times the H2 rate: at most 2 W_bar per useful iteration. Under it the meter (charged eps = 1.25 W_bar) still buys
S_tol = S (1 + eps)/(1 + 2 W_bar) = 2^125.66024 useful iterations, and the loss is the computed drop of the
collision law from S to S_tol, 0.0106656 (section 4). No martingale or tail bound is claimed for t = 2; the toy
census (evidence.md, census 5) reads realized waste at 0.87 to 0.96 of the H2 rate per field, upper 95% bounds
below 1.02.

**3.5 Detection lemma and the caps.**

*Detection lemma.* Fix a trail and ignore its DP bits: its path, including the class squared at every escape, is
a function of the class bits other than bits 64..149 (the escape compares the two keys with those bits masked;
the closure tests compare full keys, but two keys are equal exactly when their classes are, up to the tie events
charged in section 4, and a compare reads the DP bits only when two distinct classes agree on the other covered bits,
also charged there), and under H2 the DP bits of distinct classes are independent of those and of each other, each zero with
probability theta. So for a trail of U counted steps P(no DP) = E[(1 - theta)^(D_U)], D_U the number of distinct
classes it visits. Outside a self-collision inside the trail, D_U >= U/4: every segment is a run of distinct
classes until a closure, and a closure with mu + l distinct classes costs at most 3 (mu + l) + l <= 4 (mu + l)
counted steps before the escape, whose squared class starts a new segment of new classes. A trail is *stuck* if it
reaches A = 120/theta counted steps without a DP; then P(stuck) <= (1 - theta)^(A/4) + P(self-collision within A
steps) <= e^-30 + A^2/N = 2^-43.28 (A^2/2 pairs at the pairwise rate 2/N, H1-upper). The colliding trails are size-biased (the collision lands
in proportion to length); with D >= U/4 the probability that either reaches no DP within A steps is at most
(a + 4) e^(-a/4) each, a = 120: 2.32e-11 for the two.

*Caps, all by Markov.* Reports: a trail ends at a DP on a fresh class with probability theta per first visit, and
a merge into an earlier trail adds at most one report; merges are pairs of first visits. So
E[reports] <= n_max theta + n_max^2/N + M = 2^39.73 (pairs at 2/N, H1-upper), and P(reports >= 2^55) <= 2.53e-5
(a factor 39,514: 2^55 over the bound on E[reports] with the walk starts counted).
Restarts: each trail is stuck with probability at most 2^-43.28, so E[restarts] <= 2^-43.28 (E[reports] + M)/
(1 - 2^-43.28) = 2^-3.55 and P(restarts >= 2^20) <= 8.1e-8. A trail that steps onto the DP-free path of an earlier
stuck trail inherits it and is stuck too: each step meets at most (restarts so far) x A such classes at 2/N
(H1-upper), so E[restarts] also carries the factor 1/(1 - 2^125.730 x 2A/N) = 1/(1 - 2^-31.6), and the merge event
is a charged loss (section 4). Pair budget: 2.1e-7 (section 3.2). Within-trail
collisions: a trail has at most 4 D counted steps per distinct class, D geometric with mean 1/theta, so
E[U^2] <= 32/theta^2; over the E[reports] + M trails and the at most 2^20 restarted ones of A steps the pairwise
bound at 2/N (H1-upper) gives below 4.0e-11, charged 1e-6.

## 4. Success probability

**H1 (random mapping).** The class sequence behaves like iterates of a random mapping on N classes in which two
distinct classes have the same image with probability (1 - sum q^2)/N, q the multiplier-index probabilities. As in
the parent: for two classes with different indices the image equation has about one solution among N classes; for
the same index about 11/12 of one (the twist adds 11 nontrivial chances). We use the conservative sum q^2 = 1/R =
1/24 (first-order value 1/(12R) = 1/288), the Brent-Pollard factor R/(R-1) of an R-adding walk:
**N' = N/(1 - 1/24), log2 N' = 251.33353.** H1 is used in two directions and the two are kept apart: *H1-lower*
(the success law) takes the enlarged N', which only lowers the success; *H1-upper* (every union bound in sections
3.5 and 4: stuck trails, expected reports, within-trail collisions) takes the pairwise rate 2/N, the random-mapping
1/N with a factor 2 of margin, and never the smaller 1/N'. Part E samples 3,000 classes at q = 67, R = 24 in the declared run (1.206 +- 0.154) and
15,000 in the extended run of evidence.md, which reads N times the pair-collision rate 1.0068 +- 0.0338 against
0.9583 (conservative), 0.9965 (first-order) and the bound 2. H1 does
not depend on how the keys are computed: section 3.1 shows the lazy walk visits exactly the elements of the eager
multiply-only walk.

Under H1 (section 3.4) the probability that the first S useful images contain no class collision is
exp(-S (S - 1)/(2N')), within a factor 1 + 2^-125 of exp(-S^2/(2N')), which the floor absorbs. With
P_model = 0.4014: **S = sqrt(2 N' ln(1/0.5986)) = 2^125.68551**, S^2/N' = 1.02632. Losses, each computed in
part D:

* H1+ at tolerance t = 2 (section 3.4): the meter buys S_tol = 2^125.66024 useful iterations, and
  P_model - (1 - exp(-S_tol^2/(2N'))) = 0.0106656.
* The restart cap's worst case as steps: 2^20 x A = 2^112.907 steps, a fraction f = 1.4232e-4 of S. With
  x = S^2/N' the model success is P(S) = 1 - exp(-x/2) and P'(S) S = x exp(-x/2) = 0.61436 < 1, so losing a
  fraction f of the S steps costs at most f of success; the ledger charges f.
* Report cap exceeded (Markov, section 3.5): 2.53e-5. Restart cap exceeded: 8.1e-8. Pair budget exceeded: 2.1e-7.
* Detection failure on either colliding trail (detection lemma, section 3.5): 2.32e-11.
* A collision inside a single trail (section 3.5): below 4.0e-11, charged 1e-6.
* The sampler's exhaustion allowance 2^-71.99, a key tie or a zero key (at most 70/p per keyed element over at most
  2^125.730, below 2^-240), a masked key tie at an escape (two distinct candidate classes agreeing on the 176 key
  bits outside 64..149 that H2 covers, bits 0..63 and 150..261, 2^-176 per pair under H2, at most 65 candidates per
  2-cycle or 64-step escape and 2^125.730 escapes: below 2^-39; a doubling escape compares each of its l traversal
  candidates against a running minimum and a tie matters only at that minimum, at most l 2^-176 (1 + 2^-83) per
  escape, with l weighted by the part F rate 2^-43.86 per step over 2^125.730 steps: below 2^-94; the bits above 261 are not relied on), a closure compare (2-cycle, 64-step
  or doubling checkpoint, at most three per iteration) meeting two distinct classes equal on those 176 bits, when its
  outcome would read the DP bits (3 x 2^125.730 x 2^-176, below 2^-48), a later trail stepping onto a stuck trail's
  DP-free path and inheriting it (each of at most 2^125.730 steps meets at most E[restarts] x A such classes at 2/N,
  H1-upper: below 2^-35; the same factor 1/(1 - 2^-31.6) multiplies E[restarts] above) and beta^1 = beta^2 (2^-255).

The sum is 0.0108345: **success at least 0.4014 - 0.0108345 = 0.3905655, floored to 0.39056**, on the required
0.39 (part D asserts the floored value and that it is at least 0.39). The 0.39 bound is conditional on H1+ at
t = 2: the bound falls by 0.0137 per unit of t (part D), so it still holds for any t up to 2.04 and fails at t = 3. The fresh starts of the trails (one per 2^86
steps) are uniform classes, which collide at rate 1/N >= 1/N' (H1-lower), so counting them as iterates is conservative.

Toy support (evidence.md, census 5: the protocol registered for census 3, rerun with its success threshold at the
current model point and the masked-key escape rule; this package's complete walk at R = 24
with trails, replay and the matcher in the toy fields q = 67, 163, 211, 1,200 runs each, the ledger's units of
section 3.4): the fraction of runs whose first useful collision is paid within the meter's S(0.4014)(1 + eps)
units is 0.4250, 0.4383, 0.4292 per field, pooled 0.4308 (Wilson 95% interval 0.4147 to 0.4471) over 3,600 runs,
and the realized waste is 0.87 to 0.96 of the H2 rate W_bar (pooled 0.906, 95% interval 0.876 to 0.935), inside
the H1+ tolerance of 2. The toy derives its class bits from a hash of the canonical representative (a toy key in
F_q has too few bits) and escapes by the production rule, squaring the class whose hash is smaller with its DP bits
masked; the production rule's own bits of K are tested in the real GT by part G. They support H1
and H1+ at their sizes; they do not measure production success.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Counts are executed scalar AArch64 instructions as the cost model counts them. The four kernels are straight-line
programs (general-purpose registers only, no branch but the final ret, no call), so the executed count equals the
emitted count on every call; `kernels/emu.py` interprets the subset and counts every executed instruction (ldp/stp
count one), and `experiment.py` part K regenerates all four in memory and runs them against a big-integer
reference on real GT elements:

| kernel | executed instructions | what it computes |
| --- | ---: | --- |
| `lazy.s` (`gen_lazy.py`) | **4,964** | from (y, its pre-sums, gamma(Pi)): 18 products, 2 + 3 Montgomery reductions, a, b, B, a zeta, B zeta of y Pi, the 12 stored keys, K and kappa |
| `root.s` (`gen_lazy.py`) | **1,561** | the six pre-sums of y, then B, a zeta, B zeta (3 products, 3 reductions), the keys, K and kappa |
| `mul.s` (`gen.py`, parent) | **14,402** | y * Pi for a stored entry with its 54 pre-sums (flush, escapes, trie construction) |
| `sqr.s` (`gen.py`, parent) | **6,777** | Granger-Scott cyclotomic squaring (escapes) |

`lazy.s` by category (from the generator): products 3,146 (22 products of 143), Montgomery reductions 780 and
their final subtractions 100 (5 reductions), loads/stores 462, Karatsuba combinations 240, keys and minimum 227
(eleven 6-limb compares with 6 + 1 selects, two 6-limb subtractions with add-back select, six 6-limb negations p - v),
prologue/epilogue and addressing 9.

**Every step** (both kinds), analytic, 174 instructions: K reads and h = umulh(K_0, 24) 3; DP test 4; kappa and the
twist 5; child index 5; depth and flush test 3; node or leaf address 5; previous node for escapes 2; 2-cycle and
checkpoint comparisons, history ring and checkpoint store test 60; loop control and the kernel call 17; doubling
checkpoint compare 10; two-word doubling countdown (4), and on a doubling the K, node and index stores with the two-word interval copied into the countdown and doubled (8), 12; margin 35;
192-bit meter 8; 2-word stuck-walk counter 5 (`experiment.py` lists the items). **Flush** 14,420: `mul.s` 14,402,
call 8, buffer swap and node reset 3, margin 7.

| step | instructions | units |
| --- | ---: | ---: |
| root step (node = root) | 1,561 + 174 = **1,735** | 0.08675 |
| lazy step | 4,964 + 174 = **5,138** | 0.2569 |
| lazy step whose child is a leaf (with the flush) | 5,138 + 14,420 = **19,558** | 0.9779 |
| block of six steps | 1,735 + 4 x 5,138 + 19,558 = **41,845** | **c = 41,845/120,000 per step** |

**Escape** (2-cycle or checkpoint): 24,000 instructions (choose the smaller masked key 18, address 4, `mul.s` + call
14,410 (none at the root), `sqr.s` + call 6,785, resets 19 (four K fields to all-ones 13, node and depth 2, countdown and interval to 1: 4), meter 8, margin 2,756); a 64-step checkpoint escape adds
its traversal of at most 64 steps at 19,558, 65 x 64 of bookkeeping and 65 x 4 of masking (1,280,132 < A_max = 2^25); a doubling escape
charges each traversal step on its own at 19,558 + 68 (64 of bookkeeping and 4 of masking per candidate), and those
steps count on the stuck counter.

**Setup** (preprocessing, at most 2^49.06 units): pairing g (<= 2^14), z^r = 1 (<= 510), the table (24 multipliers
x (two exponentiations <= 1,021 units + 11 Frobenius maps), 48 sampler draws at 14.848 units, the scalar tables
2 x 288 products mod r at 800 instructions and 288 x 54 pre-sums at 200), **the trie**: 2^40.855 + 2^49.020 nodes at
one unit each (copy of the parent's raw coordinates 72, `mul.s` 14,410 with call, 42 pre-sums at 20, gamma 480, stores and control 266: 16,068 <= 20,000) plus
the operating system's zero-fill of its 2^60.3661 bytes (0.15625 instructions per byte, 2^43.4 units), the M = 2^20
walk starts (<= 2^10.1 units each) and the counters and constant pools: below 2^49.054 units in all.

**Unchanged from the parent** (proof-impl.md): 4,000 instructions per report (at most 2^55, 2^52.68 units); a fresh
start after every report and at every walk start or restart, at most 2^10.1 units plus 40 of meter bookkeeping
(2^65.10 units for the at most 2^55 + 2^20 + 2^20); switches 600 instructions (at most n_max/2^63 + 2^20 + 2^20 + 1,
2^57.67 units); the stopping allowance: heapsort of 2^55 records (7,188 instructions per record), dedup, at most 2^4
pairs of two replays (each below (2A + 3) x 2^26 instructions, 2^105.62 units, A = 120/theta), the check g^x = z (at
most 510), [x]Q0 in G2 (at most 2^11, below), the verification pairing (at most 2^14), halting: below 2^110.62 units.
The G2 multiplication: double-and-add over the 255-bit x in Jacobian coordinates on the twist (a = 0), at most
255 doublings (2 Fp2 multiplications + 5 squarings) and 255 mixed additions with affine Q0 (7 multiplications +
4 squarings), then one Fp2 inversion (one Fp inversion, at most 2 x 381 Fp products by exponentiation) for affine
output. With Karatsuba (3 Fp products per Fp2 multiplication, 2 per
squaring) that is at most 255 x (16 + 29) + 800 < 12,300 Fp products; at a generous 600 instructions per Fp
product including reductions, additions and control (`mul.s` spends 267) the total is below 7.4 x 10^6
instructions, 369 units; 2^11 units leaves a factor 5.5. The replay bound: one iteration repeats at most one step (at most
19,558 instructions plus 64 for the exponent update) and one escape (at most 1,280,132), below 2 x 2^25; a replayed
trail has at most A counted steps (traversals included), so 2A + 3 iterations bound it.

**Ledger** (`experiment.py` part D, 60-digit decimal arithmetic):

`W = c S + X_allow + c E + 2^110.619 (non-walk: setup, reports, fresh starts, restarts, stopping allowance) + 2^57.67 (switches)`

with c = 41,845/120,000, log2 S = 125.68551, X_allow = eps c S (eps = 0.030335, 2^119.12 units),
E = (2^20 + 2^20) x 120/theta = 2^113.907 steps: **log2 W = 124.20923**.
Submitted bound: **time_log2 = 124.210** (the least thousandth at or above it; `experiment.py` checks
W <= 2^124.210 exactly on integers, and the meter lemma walk_need + A_max <= B_walk, with 5.35e-4 of walk_need to spare).

Sensitivities (part D): eps without the 1.25 factor: -0.0085 bits; +100 instructions on every step: +0.0205 bits;
a lazy kernel twice as expensive: +0.56 bits (still below the parent's 125.111).

**Memory:** report log 2^55 x 80 bytes = 2^61.322; trie 2^60.3661 (leaves 288^6 x 2,592 bytes = 2^60.36, internal
nodes 2^52.6); table, scalar tables, walk states (2^20 x 4 KB) and replay states below 2^33: **below 2^61.93 bytes**
(< 2^62). Each array is addressed by one 64-bit register (every byte offset < 2^63). Memory is reported, not scored.

## 6. Heuristics and limitations

* Three heuristics, none proven, and the ideal-source model for the ChaCha20 sampler. H1: the class walk, with
  structural returns to the trail's current segment carved out (fruitless closures, charged in eps and left by
  an escape, never counted as collisions; section 3.4), behaves as a random mapping on classes with the
  conservative R/(R-1) correction, so the first S useful images collide with probability
  1 - exp(-S (S - 1)/(2N')) (H1-lower), and two given useful images coincide with probability at most 2/N
  (H1-upper; the enlarged N' never enters a union bound); nothing is claimed about a step given the history. H2: bits 0..261 of K are a random function of the class (next item).
  H1+: the realized waste stays within t = 2 times its H2 rate; the success the meter buys under that tolerance
  is computed (loss 0.0106656). The caps (reports, restarts, pairs) are bounded by Markov's inequality from
  first-visit and pairwise expectations, not assumed. H1, H2 and H1+ descend from the parent's heuristics; the
  carve-out and the tolerance form of H1+ are written out here. They are standard in rho
  analyses (Teske 2001; van Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999; Bos-Kleinjung-Lenstra 2010) and
  supported here by toy and production-field runs, not proved. Lazy evaluation adds no heuristic: it computes the
  same walk (section 3.1).
* H2 for K's bits, stated precisely: the map class -> bits 0..261 of K is modelled as a random function (each class
  value uniform, values of distinct classes independent). The frame kappa is not a class function (it is the index
  of the member reached, which depends on the predecessor element); as in the parent, the frame of the element
  arriving at a fresh class is modelled as uniform on 0..11 and independent of the history. Part G tests the 262
  bits only; kappa enters the tests through part C's step-by-step walk and the cycle rates of part E and the censuses. This is a model, not a distributional fact: the 262 bits
  take at most N + 1 < 2^251.28 values, so they cannot be jointly uniform over GT and no such statement is claimed;
  the model is the standard random-function one, which needs only that distinct classes get independent values
  (a 262-bit random function on 2^251.27 classes is injective except with probability below 2^-10). The rule uses
  h (essentially the top 5 bits of the low word), 86 DP bits and 112 fingerprint bits, 262 in all. Support:
  K is the minimum of the twelve stored keys, i.e. of six values min(v, p - v) of nonzero Fp-linear functions of the
  element composed with the Montgomery map, so K < p/2 and its density on [0, p/2) varies only on scale p, far above
  the 2^262 range of the bits read; part G tests the bits on 1,500 production classes (chi-square of h 16.38 on 23
  degrees of freedom, of 8 DP-field bits 232.95 on 255, of 8 fingerprint bits 297.46 on 255; largest |z| of the 262 bit
  frequencies 2.58) and part C's walk statistics. The toy censuses use a hash in place of these bits.
* R = 24 is small for an adding walk; the conservative correction 1/(1-1/R) (0.0307 bits on S) and the recomputed
  fruitless-cycle charge (0.0431 bits) are both in the ledger; P_model = 0.4014 leaves 0.00056 above 0.39 (floored bound 0.39056) after the H1+
  tolerance loss and the computed losses (census 5 reads 0.4308 pooled at this point, Wilson lower bound 0.4147). Larger R makes the trie (12R)^6 too large for memory below 2^62.
* Scope of the prices: kernel counts are executed counts from the shipped emulator; the every-step 174, the escape
  24,000, the flush overhead, the node price, the 4,000 per report, the heapsort level price and the replay bound are
  analytic bounds. The kernels are checked on real GT elements and edge values (`kernels/check_lazy.py`), not compiled
  and run on hardware. The trie (2^60.37 bytes) and the log (2^61.322 bytes) are RAM-model arrays; memory is unscored
  and no physical machine maps that much.
* This is a generic square-root attack (no index calculus in Fp12); a valid upper bound, not a security estimate.

## 7. Experiments

`experiment.py` (standard library only; parts A, S, K, C, G, B, E, F in `cycles.py`, D), its results, the toy censuses and the
kernel checks are described in `evidence.md`.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999); Teske, Math. Comp. 70 (2001); Cheon, Hong
& Kim, ASIACRYPT 2008 and J. Cryptology 25 (2012), tag tracing; Granger & Scott, PKC 2010; Bos, Kleinjung &
Lenstra, ANTS 2010; parent package gt-dprestart-replay (mpjunior92) and its ancestors (dependencies.json).

# Frobenius-class rho in GT with lazily evaluated steps: a prefix-product trie, a 3,457-instruction key kernel and one full multiplication per six steps

**Claim.** For uniform nonidentity P in G1 and uniform z in GT the algorithm below returns Q in G2 with
e(P,Q) = z with probability at least 0.39, using at most **2^123.896** dense-Fp12-multiplication equivalents
(C = 20,000 instructions), including setup (the multiplier table and the prefix-product trie), the two pairings,
walk starts and the fresh start of every trail, distinguished-point (DP) handling, fruitless-cycle escapes, the
DP detection delay, matching (sort, duplicate removal, replay of matched trails, pair checks), output and
verification. Peak memory is at most 2^63.29 bytes, with every array below 2^63 so every byte offset inside an array fits
one unsigned 64-bit register; preprocessing (included) is at most 2^52.21 units; nonuniform
advice 0 bytes. The success bound rests on the usual random-mapping heuristic for rho walks (H1), on
pseudorandom bits of the class key (H2), and on H1+, a single declared 0.01 allowance for the work threshold
binding on a run (section 3.4); each of the three caps is a failed run charged as its own loss by Markov's
inequality on its expectation, and no concentration bound or independence across iterations or trails is used. The success probability is computed in the ideal model in which the
ChaCha20 stream that supplies the table and start scalars is an ideal random source (a disclosed heuristic, no
distributional bound on the realized scalars is claimed). No production-size run is claimed.

## 0. What is new, and what is reused

This package revises this track's Frobenius-class rho package (submission 6777f95f, claim 2^125.001: a walk with a cyclotomic-squaring
branch at s = 37/64, R = 1,024 multipliers, class bits read from Re(c_0), success model 0.4001 with a 0.01 allowance
for the heuristic H1+). Kept from it: the reduction (section 1, `reduction.py`: g = e(P, Q0), output [x]Q0 in G2 at
its counted hard maximum of 12,809 units) and the two straight-line kernels `mul.s` (12,735 executed instructions;
a table entry is stored as 54 table-side linear forms reduced mod p) and `sqr.s` (6,786), from `kernels/gen.py`.
Replaced: the walk, by the lazily evaluated multiply-only walk of the package `gt-lazy-trie-rho` (submission
a77744bb, claim 2^124.183, accepted on the sibling track bls12-381-fixed-g2-invert-g1, which solves the
same discrete logarithm in the same GT; "the sibling package"). Its prefix-product trie, its key kernels `lazy.s` and
`root.s` (`kernels/gen_lazy.py`; rewritten here from 4,964 and 1,561 to 3,457 and 3,074 executed instructions:
`lazy.s` evaluates the Fp2-linear form by Winograd's inner-product identity, 9 products of 384-bit integers instead
of 18, with the node's cross term stored in the node and the base's computed once per block by `root.s`; d = a - b
comes from 2 S0 - S2 with the stored Xd and Ed, and B = Re(d_1) d takes one constant product, F5; checked in this package's
emulator), its class-bit
rule, its fruitless-cycle analysis (`cycles.py`) and its implementation sections 0-i to 0-iv (`proof-impl.md`:
meter, report log and matcher, restart cap, sampler) are used here with this package's numbers. Its success chain
is not: section 4 replaces it with the H1+ argument (no concentration bound over iterations or trails, one declared
0.01 for the work threshold, and each cap a failed run charged by Markov's inequality).

**The change: lazy evaluation of the walk (tag tracing).** A multiply step of an eager walk computes the dense
product y * T[h][t] (12,735 instructions with `mul.s`) only to read, at the next step, two Fp coordinates of the
result (the w^1 coefficient, which gives the twelve keys) and its class bits. Here a walk keeps an explicit **base**
element y and a pointer to a node of a **prefix-product trie**: the walk's current element is y * Pi_nu, where Pi_nu
is the product of the (at most five) twisted multipliers taken since the base was last made explicit. The w^1
coefficient of y * Pi_nu is an Fp2-linear form in y whose six coefficients are coordinates of Pi_nu (section F5), so
the key of the current element costs **3 Fp2 products by Winograd's inner-product identity (9 products of 384-bit
integers, Karatsuba), two Montgomery reductions and the key arithmetic: the kernel `kernels/lazy.s`, 3,457 executed
instructions**, against 12,735 +
1,910 for the previous package's multiply step. Every sixth step the base is made explicit by **one** dense
multiplication y <- y * Pi_leaf with a precomputed product of six twisted multipliers. The walk is, element by
element, the same multiply-only Frobenius-class walk as an eager walk with the same table (section 3.1); only the way
its keys are computed changes. Consequences, each charged:

* **No squaring branch** (s = 0): a squaring needs the explicit element, which the lazy walk does not have. Every
  step's price is a fixed bound (section 3.4), so the previous package's drift design margin is dropped; H1+ keeps
  0.01 for the work threshold only; every cap is a failed run whose probability Markov's inequality bounds from its
  expectation (section 3.4); the model success is 0.4005.
* **R = 32 multipliers** (previously 1,024), so the trie of depth L = 6 has (12 x 32)^6 = 2^51.51 leaves. The
  conservative collision correction is 1/(1 - 1/R) (log2 N' = 251.31793) and the fruitless-cycle rates are those of
  R = 32 (section 3.3, part F): 2-cycles at exactly 1/(12R) = 1/384 per step, all others below 3.45e-5 per step, and
  closed useless walks longer than 64 below 2^-122.6 per step (charged as 2^-121); their total waste is charged as
  eps = 0.0237922 (1.25 times the model, at most 0.023793 asserted).
* **Prices:** a root step 3,074 + 170 = 3,244 instructions, a lazy step 3,457 + 170 = 3,627, the flush 12,753, so a
  block of six steps costs 34,132 instructions: **c = 34,132/120,000 = 0.2844334 units per step** (rounded up;
  previous package: 0.51922578125 on
  average). Partial blocks cost less per step, so c bounds every step deterministically (section 5).
* **Setup** grows to the trie (2^52.21 units; each node is charged 1.6 units, 32,000 instructions, for its 31,920:
  one `mul.s` call, the 54 forms of its own `mul.s` entry, and its key entry gamma with the cross-term constants Ea,
  Ed, F5) and **memory** to 2^62.855 bytes of trie beside the
  2^61.322-byte report log, 2^63.29 in all; memory is unscored and every array stays below 2^63, so every byte
  offset inside an array fits one unsigned 64-bit register.
* **Class bits** are read from the key minimum K (not from Re(c_0), which the lazy walk does not have): the multiplier
  index h = floor(32 K_0 / 2^64) from the low word, the DP predicate on bits 64..149, the 128-bit fingerprint on bits
  150..277 (278 bits in all, fewer than the 2^251.27 classes carry, section 6).

Ledger: log2 W = 123.89581 (`experiment.py` part D), claim **123.896**, the least thousandth above it, checked by an
exact integer inequality.

## 0a. Representation: which limbs the rule reads

As in both earlier packages, an Fp value is held as the six-limb Montgomery vector x~ = x 2^384 mod p, fully reduced in [0, p),
and the rule is stated on these stored limbs: the twelve keys are the stored residues of section F2 (linear
operations commute with the Montgomery map, so stored keys are the Montgomery images of the integer keys and
their multiset is again Frobenius-invariant), K is their minimum as 381-bit integers, kappa the first index (in the
kernel's fixed order) attaining it. The rule on stored limbs is the integer rule composed with a fixed bijection
of Z/p, so every structural property used below carries over. `lazy.s` and `root.s`
return exactly the stored a, b, B, a zeta, B zeta, K of the current element (part K and part C check this against
a big-integer reference, including keys from explicit Frobenius maps).

## 1. Reduction

The input is a nonidentity P in G1 (order r) and z in GT. Compute g = e(P, Q0) with the pinned convention
pairing(Q0, P). The reduced pairing is bilinear and nondegenerate and Q0 generates G2, so for every permitted P,
g != 1 generates GT (prime order r) and z = g^x for a unique x in [0, r). If z = 1 output Q = infinity (x = 0).
Otherwise find x by the walk below, output Q = [x]Q0 in G2 and check e(P, Q) = z; bilinearity gives
e(P, [x]Q0) = g^x = z for the supplied P. `reduction.py` (part R of `experiment.py`) runs this on BLS12-381: the
pinned pairing operation by operation (10,331 counted operations, each charged one unit; charged 2^14) and
Q = [x]Q0 by double-and-add with the pinned projective formulas over Fp2 and one Fp2 inversion. The ladder's hard
maximum takes the most expensive branch of the pinned addition (28 operations, the equal-points branch) and doubling
(20) at each of at most 254 steps plus the inversion (617): **12,809 units**, every Fp2 operation charged a full
dense Fp12 multiplication, an overcount.

## 2. Arithmetic facts (checked in `experiment.py`, parts A, K, C)

Tower: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi), xi = 1+u (blst's tower; with u = w^6 - 1 it is the
profile field). y = sum_{m<6} y_m w^m with y_m in Fp2.

**F1 Frobenius classes.** p has order 12 modulo r, p^6 = -1 and 1 + p^4 + p^8 = 0 (mod r), so y -> y^{p^k} is a
group of order 12 acting on GT and containing inversion. Every y != 1 lies in a class of exactly 12 elements;
N = (r-1)/12 classes, log2 N = 251.27213.

**F2 Keys (both earlier packages').** zeta = xi^{(p^2-1)/6} is a primitive sixth root of unity in Fp, d_1 = xi^{(p-1)/6}.
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
Fp-bilinear in (y, gamma(Pi_nu)). **Winograd's inner-product identity over Fp2**: with the pairs (u, v) = (0, 1),
(2, 3), (4, 5),
`c_1 = sum_(u,v) (x_u + gamma_v)(x_v + gamma_u) - X(x) - E(Pi),  X(x) = sum_(u,v) x_u x_v,  E(Pi) = sum_(u,v) gamma_u gamma_v`,
because (x_u + gamma_v)(x_v + gamma_u) = x_u gamma_u + x_v gamma_v + x_u x_v + gamma_v gamma_u in the commutative
Fp2. X depends on the base only and is computed once per block by `root.s`; E depends on the node only and is stored
in it. With A = x_u + gamma_v and B = x_v + gamma_u, the three Fp2 products are Karatsuba products,
`S0 = sum Re A Re B, S1 = sum Im A Im B, S2 = sum (Re A + Im A)(Re B + Im B)` (sums over the three pairs),
nine products of 384-bit integers. Every operand is the sum of two stored values in [0, p) (coordinates, or the
pre-sums x_k.re + x_k.im and gamma_k.re + gamma_k.im, both reduced mod p), so at most 2p - 2 < 2^384, and S0, S1, S2
lie in [0, 3 (2p - 2)^2]. Then a = S0 - S1 - Re X - Re E and **d = a - b = 2 S0 - S2 - (Re X - Im X) - (Re E - Im E)**
modulo p. Per base, `root.s` stores the 768-bit Xa = p 2^384 - (s0 - s1) and Xd = p 2^384 - (2 s0 - s2), with s0, s1,
s2 X's own Karatsuba sums on the stored coordinates (each in [0, 3 (p - 1)^2]); per node, the trie stores Ea and Ed in
[p, 2p), Ea = p + (-(e0 - e1) 2^-384 mod p) and Ed = p + (-(2 e0 - e2) 2^-384 mod p), with e0, e1, e2 E's Karatsuba
sums on the node's stored values. `lazy.s` forms a_T = S0 - S1 + Xa + Ea 2^384 and d_T = 2 S0 - S2 + Xd + Ed 2^384, which are
congruent mod p to a 2^384 and d 2^384 in the Montgomery scale and lie in [0.47 pR, 4.53 pR] and [0.17 pR, 5.75 pR]
with R = 2^384, so both lie in [0, 2^768) (2^768 > 9.84 pR); one Montgomery reduction of each gives
r <= T/2^384 + p < 6.75p, inside the q-table subtraction of `gen.py` (r < 8p); `gen_lazy.py` asserts these ranges on
integers. Then b = a - d mod p, and since
Im d_1 = -Re d_1 (asserted), B = Re(d_1) d and B zeta = (Re(d_1) zeta) d: one constant product each, plus a zeta;
the keys and K follow F2.

**F3 Cyclotomic squaring (escapes only).** The Granger-Scott map (`kernels/sqr.s`, 6,786).

**F4 Orders.** r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 * 2529403 * 52437899 *
254760293^2 (checked); ord_r(2) = (r-1)/16, so 2^t lies in <p> only for t a multiple of about 2^247.3.

## 3. Algorithm

Parameters: **R = 32** table multipliers, **trie depth L = 6**, DP density **theta = 2^-86**, M = 2^20 walks
round-robin on one processor (chunk K = 2^63 steps; proof-impl.md 0-i), report cap 2^55, restart cap 2^20,
pair budget 2^4 (section 3.4).

**Table.** M_h = g^{a_h} z^{b_h} for h < 32 with fresh uniform a_h, b_h mod r from the sampler (0-iv), and the
twisted entries T[h][t] = M_h^{p^t}, t < 12: 384 entries (each with its 12 coordinates and the 54 table-side
forms of `mul.s`), and the scalar tables a_h p^t, b_h p^t mod r (read only by replays).

**Trie.** A node nu at depth d <= 6 is a word ((h_1,t_1), ..., (h_d,t_d)) of table entries, with
Pi_nu = T[h_1][t_1] * ... * T[h_d][t_d]; the root is the empty word (Pi = 1). Nodes are stored in level order:
the children of node n are 384 n + 1 + 12 h + t. An internal node (depth 1..5) stores its 20-vector key entry
(gamma_k.re, gamma_k.im and their sum mod p for k < 6, then Ea and Ed, F5) and the 54-vector `mul.s` entry of Pi_nu
(for escapes): 74 x 48 = 3,552 bytes;
a leaf (depth 6) stores the `mul.s` entry of Pi_nu: 2,592 bytes. Internal nodes: sum_{d=1}^{5} 384^d = 2^42.929;
leaves 384^6 = 2^51.510; 2^62.855 bytes in all (every array below 2^63). A node is built from its parent by one `mul.s` call
Pi_child = T[h][t] * Pi_parent (the entry's 12 coordinates against the parent's stored forms), then the 54 forms
of its own entry and its key entry (31,920 instructions, charged 32,000 = 1.6 units, section 5).

**Walk state.** Two base buffers (the current explicit base y with its per-base data, the six pre-sums
x_k.re + x_k.im mod p and the 768-bit Xa and Xd of F5, and the previous base), the current node and its depth, a ring of three K slots (the current K and
the two before), the previous step's node and buffer index (for escapes), the checkpoint K, the step count mod 64,
the 2-word stuck-walk counter and the trail's start pair (alpha0, beta0). It carries no running exponents.

**One iteration** (pre-debit charge first, proof-impl.md 0-i):

```
1  if node = root: (K, kappa) <- root.s(y)            [also writes y's per-base data: pre-sums, Xa, Xd]
   else:           (K, kappa) <- lazy.s(y, per-base data, key entry of Pi_node)    [K of the element y * Pi_node]
2  if bits 64..149 of K are 0 (distinguished): F <- bits 150..277 of K; append (F, alpha0, beta0) to the report
   log; the trail ends: fresh start y <- g^alpha0 z^beta0 from a new start pair, node <- root, histories reset
3  2-cycle and checkpoint tests on K (section 3.3); on a hit: escape, then next iteration
4  h <- floor(32 K_0 / 2^64), t <- (12 - kappa) mod 12, node <- child(node, h, t)
5  if depth(node) = 6: y_prev <- y, y <- mul.s(y, entry(Pi_node)) (the flush), node <- root
```

A replay runs the same iteration from (alpha0, beta0) and additionally adds (a_h p^t, b_h p^t) to its exponent
pair on every step and doubles it at an escape's squaring, as in the sibling package; it stops at the trail's distinguished
class and returns p^kappa (alpha, beta) there.

**3.1 The step is a map on classes, and the lazy walk is the eager walk.** The current element is
y * Pi_node exactly (an invariant: it holds at a start, at an escape and after a flush, where node = root, and a
step multiplies Pi_node by T[h][t] on the right). By F5 the key computed in step 1 is the key of that element, so
(K, kappa) and therefore h, t, the DP test and the fingerprint are those of the current element; the next
element is (y Pi_node) T[h][t] = (y^ M_h)^{p^{-kappa}}, y^ the canonical representative, which lies in the class
[y^ M_h]. So the step depends only on the class, and the sequence of elements (not only of classes) is exactly
that of an eager multiply-only walk run with this table. The declared run of part C takes 240 production steps (49 root,
191 lazy, 26 flushes, one escape) with every key from the emulated `root.s`/`lazy.s` and every flush and escape from
the emulated `mul.s`/`sqr.s`, and asserts at every step that y * Pi_node equals the eager walk's element and that
(K, kappa) equals the reference; it also checks the exponent invariant at the end.

**3.2 Collision.** As in the sibling package: every report is appended to the log until the cap; after the walk phase the
log is heapsorted, duplicates dropped, and every pair inside a run of equal F is replayed and checked, under a hard
budget of 2^4 pairs. Two records with equal F are the same class up to a false match of probability 2^-128 per
pair of distinct classes under H2 (fewer than (2^55 choose 2) x 2^-128 < 2^-19 false pairs in expectation for 2^55
records, and the degenerate pair beta^1 = beta^2, below 2^-254 each; both are charged before the pair budget, section
3.4). The replays give canonical
pairs with alpha^1 + x beta^1 = alpha^2 + x beta^2 mod r, so x = (alpha^1 - alpha^2)/(beta^2 - beta^1) unless
beta^1 = beta^2 (about 1/r). Output Q = [x]Q0 after checking g^x = z.

**3.3 Fruitless cycles at R = 32.** A closed walk of the class graph is *useless* if the canonical exponent pair
returns unchanged. After l steps the pair is Lambda (alpha^_0, beta^_0) + tau with Lambda = 2^s p^m (s squarings, m a
sum of twists) and tau = sum_h c_h (a_h, b_h), c_h in Z[p]. Because (alpha^_0, beta^_0) is uniform and independent of
the table, the walk is useless (up to probability about 1/r) only if Lambda = 1 and every c_h = 0. An escape's
squaring cannot lie in a useless walk (Lambda = 2^s p^m = 1 needs s >= 2^247, F4). For multiply steps, under H2 the
indices are uniform on 32 and the frames p^e of the occurrences are independent uniform on the 12 powers of p (each
step's twist is (12 - kappa) mod 12 with kappa a fresh uniform class bit), and c_h is the sum of the frames of
index h's occurrences. For coefficients far below r, c_h = 0 mod r exactly when this sum of 12th roots of unity
vanishes in Z[zeta_12]. So the probability P_l that the l steps from a given step form a useless walk is at most
E[prod_h g(n_h)] over the multinomial counts n_h of the 32 indices, with g(0) = 1, g(1) = 0 and g(n) = q(n), the
exact probability that n independent uniform 12th roots of unity sum to 0 (computed exactly for n <= 12 in part
F: q(2) = 1/12, q(3) = 1/72, q(4) = 0.0191, ..., q(12) = 0.00195), and g(n) = Q(12) = 0.00195, the largest point
probability of a 12-term sum, for n > 12 (adding independent terms cannot raise a point probability). Part F
evaluates the bound exactly for l <= 384 and, beyond, by Poissonisation (P_l <= e sqrt(l) (E_Pois(l/32) g)^32,
monotone on doubling segments), weighted by the probability exp(-theta l) that a closed walk of length l holds no
distinguished class (H2). Results: **P_2 = 1/384 exactly** (same index twice, opposite twist), sum_{l=3}^{64} P_l
below 3.45e-5 (P_3 = 2^-16.17, P_4 = 2^-15.59, P_8 = 2^-27.61, P_64 = 2^-122.33), and useless walks longer than 64 at
a rate **below 2^-122.6 per step** (charged as 2^-121).

* **Detection and escapes (on K).** The K of the step before last is remembered; if it equals
  the current K, the two classes form a 2-cycle: the walk makes explicit the representative of the class with the
  smaller K (the current element y * Pi_node, or the previous element, which is Pi of the previous node times the
  base it was computed from: y, or y_prev if a flush came between; one `mul.s` call with that node's entry, or no
  call when that node is the root) and replaces it by its square (`sqr.s`), node
  <- root, and resets the two-step history and the checkpoint; a replay doubles its pair there. A checkpoint K is
  stored every 64 steps and compared at every step; a repeat means a cycle of length at most 64, which is traversed
  once more to its minimum-K class, which is then squared the same way. Both choices depend only on the cycle, so
  merged trails stay merged. Every walk whose stuck-walk counter reaches 28/theta = 28 * 2^86 steps is restarted
  (0-iii); the counter is reset at every walk start, restart, report and fresh start.
* **Charge.** A 2-cycle wastes the step that re-enters the cycle and the detecting step (2 steps) and one escape
  (24,000 instructions, 1.2 units = 4.2189 steps at c): (2 + 4.2189)/384 = 0.0161951 step-equivalents per step. A cycle of
  length 3 <= l <= 64 wastes at most 64 + l steps before its detection, its traversal of l steps at the largest
  step price (16,380, 2.8794 c each) and one escape: sum_l P_l (64 + l + 2.8794 l + 4.2189) = 0.0028387. The ledger
  charges **eps = 1.25 (0.0161951 + 0.0028387) = 0.0237922** on the step count (part D computes it on exact
  fractions and asserts eps <= 0.023793); the factor 1.25 covers excess over the H2 rates
  (the toy census reads pre-collision 2-cycle rates of 0.002518 to 0.002608 per step at R = 32, model 0.002604,
  and checkpoint escapes at 1.39e-5 to 2.44e-5 per step against the bound 3.45e-5; evidence.md).
  Useless walks longer than 64 with no distinguished class are left to the stuck-walk rule at the charged rate
  2^-121: below 2^4.71 expected restarts over the run's at most 2^125.710 iterations (0-iii).

**3.4 Caps and the work cap.** Setup (at most 2^52.21 units), the walk phase, and the stopping
phase (matching, output, verification, halting; the stopping allowance R_stop < 2^108.52 units, section 5). The
walk phase has the budget B_walk = W_cap C - 2^52.21 C - R_stop C instructions (W_cap = 2^123.896, W_cap C rounded
down) and meters its own charged work in the 192-bit meter under the pre-debit rule (proof-impl.md 0-i) with the
prices of section 5: before each charged action the meter is compared with L = B_walk - A_max (A_max = 2^25) and the
action runs only if the meter is below L. So every run does at most 2^52.21 C + B_walk + R_stop C <= W_cap C
instructions, a deterministic inequality `experiment.py` checks on integers. The step price is fixed per step type
and known before the step: root 3,244, lazy 3,627, lazy with flush 16,380 (the step whose child is a leaf).

*No price-drift heuristic.* A block is the run of iterations from a root iteration to the next return to the root
(a flush, an escape or a fresh start). Its iterations are charged 3,244 (the first), 3,627 (the second to fifth) and
16,380 (the sixth, the only one whose child can be a leaf), so a block of k <= 6 iterations is charged at most
3,244 + (k-1) 3,627 < k c C for k < 6 and exactly 6 c C for k = 6, whether or not its last iteration ended in a
report or an escape (the escape and the fresh start are charged separately). Hence the charged step work of n
iterations is at most n c whatever the walk does, and the walk budget buys at least (budget - non-step work)/c
iterations (the eps of section 3.3 counts the iterations wasted in fruitless cycles); the only random rates the budget depends on are
those of escapes, reports and restarts, which are charged separately.

*Iteration maximum.* Counting also the action that ends a block, every block of k iterations is charged at least
k c C: a block ended by a flush is charged exactly 6 c C; one ended by an escape at least 3,244 + (k-1) 3,627 + 24,000
>= k c C (k <= 6); one ended by a report at least 3,244 + (k-1) 3,627 + 4,000 plus its fresh start, charged its fixed
allowance of 2^10.1 units (proof-impl.md 0-i; its two 255-bit exponentiations alone exceed 2^21 instructions)
(`experiment.py` asserts all three for every k). So a run has at most n_max = B_walk/(c C) + 6M < **2^125.710
iterations** (the last, incomplete block of each walk adds at most 6).

*H1+, the disclosed success heuristic: no concentration bound is taken on any cap.* The step count the budget buys is
deterministic given the waste, and the waste is charged at 1.25 times its H2 expectation (section 3.3); the walk
budget's slack over the ledger's need is a further 0.0134% of the step bucket, asserted in part D. No probabilistic
statement is made about how far the realized fruitless-cycle waste can sit from its H2 expectation, and none needs to
be: the budget, not the realized waste, is what the ledger spends, and a walk that wastes more still cannot exceed the
meter. The caps are not charged that way: each is a nonnegative integer count whose expectation E is bounded in
part D, so **Markov's inequality** holds for the walk whatever its dependence structure and gives
P(cap reached) <= E / cap, and that quotient is its own line of the loss sum (section 4). Reaching the report cap
(2^55), the restart cap (2^20) or the pair budget (2^4) **fails the run**. The 0.01 that H1+ keeps covers one event
only: the work threshold ending the walk phase before the ledger's walk is done, when fewer than S (1 + eps)
iterations were bought. H1+ states, and does not prove, that the H2 rates
hold over the run as a whole rather than only in expectation per step, which is what makes the three expectations
the Markov bounds divide by bounds over the whole run. Support, not proof, for the rates themselves:
part C's production walk, part G's key-bit statistics, part E's pair-collision law, part F's cycle bounds
(`cycles.py`, section 3.3) and the toy censuses (evidence.md).

*Caps.* Reports D <= 2^55 with one-word counter, at most n_max theta < 2^39.710 expected, so Markov bounds the cap
by 2^39.710 / 2^55 = 2^-15.29 and the cap sits more than 2^15.29 above the expectation. At most
2^4 matcher pairs. A true pair (two records of one class from different trails) verifies unless beta^1 = beta^2, and
the first verifying pair ends the run, so a run's matcher charges at most one verifying pair, and the only useless
true pair is one with beta^1 = beta^2. The pairs charged before the halting one are therefore false-key pairs
(2^-128 each under H2) or pairs with beta^1 = beta^2 (below 2^-254 each), at most (2^55 choose 2) (2^-128 + 2^-254) < 2^-19 in expectation, so Markov
bounds the pair budget by 2^-19 / 2^4 = 2^-23 and the budget sits more than 2^23 above the expectation. At most 2^20 stuck-walk restarts, expected below n_max theta e^-28 + n_max 2^-121 < 2^4.75 (0-iii), so
Markov bounds the cap by 2^4.75 / 2^20 = 2^-15.25 and the cap sits more than 2^15.25 above the expectation. Part D asserts each ratio (markov_ratio_report >= 40,000, markov_ratio_restart >= 39,000,
markov_ratio_pairs >= 8,000,000) and its JSON prints each expectation, cap, ratio and bound. The ledger's step count
is S (1 + eps) + 2 x 2^20 x 28/theta, the detection window and restart allowance being below 2^111.81 steps.

## 4. Success probability

**H1 (random mapping).** The class sequence behaves like iterates of a random mapping on N classes in which two
distinct classes have the same image with probability (1 - sum q^2)/N, q the multiplier-index probabilities. As in
the sibling package: for two classes with different indices the image equation has about one solution among N classes; for
the same index about 11/12 of one (the twist adds 11 nontrivial chances). We use the conservative sum q^2 = 1/R =
1/32 (exact first-order value 1/(12R) = 1/384), the Brent-Pollard factor R/(R-1) of an R-adding walk:
**N' = N/(1 - 1/32), log2 N' = 251.31793.** Part E samples classes at q = 67, R = 32 and reads N times
the pair-collision rate against 0.9688 (conservative) and 0.9974 (exact): 1.0447 on the declared run's 15,000
classes, 0.9786 on 40,000 and 0.9863 on 60,000 (relative standard errors 0.0331, 0.0128, 0.0085; evidence.md). H1 does not depend on how the keys are computed: section 3.1 shows the lazy walk visits exactly the
elements of the eager multiply-only walk.

Under H1 the probability that the first S steps contain a class collision is 1 - exp(-S^2/(2N')). With
P_model = 0.4005: **S = sqrt(2 N' ln(1/0.5995)) = 2^125.67559.**

The losses before the first useful collision is detected, each computed in part D:

* **the 0.01 of H1+ (section 3.4)**: the event that the work threshold ends the walk phase before the ledger's walk
  is done, a declared allowance for the H2 rates holding over the run as a whole (section 3.4). It covers no cap:
  the caps are charged separately, next;
* **each cap's Markov bound** (section 3.4), a failed run if the cap is reached, so its probability is at most the
  expectation over the cap that part D computes: the report cap more than **2^15.29** above its expectation below
  2^39.710, at most **2^-15.29**; the restart cap more than 2^15.25 above its expectation below 2^4.75, at most
  **2^-15.25**; the pair budget 2^23 above its expectation below 2^-19, at most **2^-23** (part D: 2.4958e-5,
  2.5551e-5 and 1.1921e-7);
* the restart cap's own worst case spent as steps: 2^20 restarts at one full stuck-walk allowance each,
  2^20 x 28/theta < 2^110.81 steps, a fraction below **3.27e-5** of the iteration maximum 2^125.710, charged in full and not
  at its expected count (0-iii);
* detection failure on either of the two colliding trails, at most 3.88e-11 in all: under H2 the colliding step is
  independent of the DP bits after it, so it lies on a trail the stuck rule cuts with probability at most
  (28/theta) P(T > 28/theta)/E[min(T, 28/theta)] = 28x/(1 - x), T the geometric trail length at rate theta and
  x = (1 - theta)^(28/theta) <= e^-28: below 1.94e-11 per trail (part D asserts 28 e^-28/(1 - e^-28) < 1.94e-11);
* a collision inside a single trail, at most 1e-6: under H1 a step meets one of its trail's at most 2 x 28/theta + 3
  earlier classes with probability at most that count over N, below 2^-33.7 over n_max steps (asserted in part D);
* the sampler's exhaustion allowance under the ideal-source model of 0-iv, 2^-71.99 < 2.14e-22;
* a key tie or a zero key at most 70/p per keyed element over at most 2^125.710, below 2^-240; a start or iterate equal
  to the identity, 1/r each, below 2^-128 in all; beta^1 = beta^2, below 2^-254.

The sum is **0.01008428**, below 0.0100843 (part D computes it at 60 digits and asserts the bound). So the success
after every loss is at least 0.4005 - 0.0100843 = 0.3904157, **floored to 0.39041**: 0.00041 above the required 0.39.
Part D asserts the floored value. The fresh starts of the trails (one per 2^86 steps) are uniform classes, which
collide at rate 1/N >= 1/N', so counting them as iterates is conservative.

**H1+ and no independence.** Nothing in this section rests on the iterations or the trails being independent:
Markov's inequality on a nonnegative integer count assumes nothing about dependence, so the three cap losses hold for
the adaptive, stateful walk too, and the 0.01 of H1+ is the one declared, unproved allowance, for the work threshold
ending the walk phase early; it is not inferred from any bound over iterations or trails.

Toy support (evidence.md: a census of this package's complete walk at R = 32, s = 0 with trails, replay and the
matcher in toy fields q = 67, 163, 211, at this package's charges): the fraction of runs whose first useful collision
falls within S(0.4005)(1 + eps), counted in the ledger's step units. Its eight rows and its pooling rule were
registered for this package's charges before any count at those charges was read; the seeds and walks are those of
earlier censuses of this walk, whose counts at earlier charges were known, and no new sample was drawn. The eight
rows read 0.4060 to 0.4296 (floored) and pool to **0.4201 +- 0.00384 over 16,550 runs** (6,953 successes), 0.4047 (floored) at
four binomial standard errors below the pool. The mean T_c/sqrt(N) is 1.239 to 1.272 against the model's 1.2734. The census measures the threshold and waste part of H1+ only,
not any cap: the caps are sized by Markov bounds on their expectations (section 3.4), which no toy run can
measure. The toy derives its class bits from a hash of the canonical representative (a toy key in F_q has too few
bits); the production rule's own bits of K are tested in the real GT by part G. They support H1 at their sizes; they
do not measure production success.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Counts are executed scalar AArch64 instructions as the cost model counts them. The four kernels are straight-line
programs (general-purpose registers only, no branch but the final ret, no call), so the executed count equals the
emitted count on every call; `kernels/emu.py` interprets the subset and counts every executed instruction (ldp/stp
count one), and `experiment.py` part K regenerates all four in memory and runs them against a big-integer
reference on real GT elements:

| kernel | executed instructions | what it computes |
| --- | ---: | --- |
| `lazy.s` (`gen_lazy.py`) | **3,457** | from (y, its per-base data, gamma(Pi)): 9 + 3 products (Winograd over Fp2), 2 + 3 Montgomery reductions, a, b, B, a zeta, B zeta of y Pi, the 12 stored keys, K and kappa |
| `root.s` (`gen_lazy.py`) | **3,074** | the per-base data of y (six pre-sums mod p, then Xa and Xd from 9 products), then B, a zeta, B zeta (3 products, 3 reductions), the keys, K and kappa |
| `mul.s` (`gen.py`) | **12,735** | y * Pi for a stored entry of 54 table-side forms reduced mod p (flush, escapes, trie construction) |
| `sqr.s` (`gen.py`) | **6,786** | Granger-Scott cyclotomic squaring (escapes) |

`lazy.s` by category (the generator's own instruction tags): products 1,716 (12 products of 143: nine for the
Winograd form and three constant products), Montgomery reductions 780 and their final subtractions 100 (5
reductions: a, d, B, a zeta, B zeta), loads/stores 354, the operand sums, the accumulation of the three sums, the Xa,
Ea, Xd, Ed terms and the doubling of S0 252, keys
and minimum 245 (b = a - d mod p, eleven 6-limb compares with 6 + 1 selects, two 6-limb subtractions with add-back
select, six 6-limb negations p - v), prologue/epilogue and addressing 10.

**Every step** (both kinds), analytic, 170 instructions: K reads and h = umulh(K_0, 32) 3; DP test 4; kappa and the
twist 5; child index 5; depth and flush test 3; node or leaf address 5; previous node for escapes 2; 2-cycle and
checkpoint comparisons, history ring and checkpoint store test 60; loop control and the kernel call 17; margin 53;
192-bit meter 8; 2-word stuck-walk counter 5 (`experiment.py` lists the items). **Flush** 12,753: `mul.s` 12,735,
call 8, buffer swap and node reset 3, margin 7.

| step | instructions | units |
| --- | ---: | ---: |
| root step (node = root) | 3,074 + 170 = **3,244** | 0.1622 |
| lazy step | 3,457 + 170 = **3,627** | 0.18135 |
| lazy step whose child is a leaf (with the flush) | 3,627 + 12,753 = **16,380** | 0.819 |
| block of six steps | 3,244 + 4 x 3,627 + 16,380 = **34,132** | **c = 34,132/120,000 < 0.2844334 per step** |

**Escape** (2-cycle or checkpoint): 24,000 instructions (choose the smaller key 14, address 4, `mul.s` + call 12,743
(none at the root), `sqr.s` + call 6,794, resets 20, meter 8, margin 4,417); a checkpoint escape adds its traversal of
at most 64 steps at 16,380 and 65 x 64 of bookkeeping (1,076,480 < A_max = 2^25).

**Setup** (preprocessing, at most 2^52.21 units): pairing g (<= 2^14), z^r = 1 (<= 510), the table (32 multipliers
x (two exponentiations <= 1,021 units + 11 Frobenius maps), 64 sampler draws at 14.848 units, the scalar tables
2 x 384 products mod r at 800 instructions and 384 x 270 modular additions for the entries' forms at 60), **the
trie**: 2^42.929 + 2^51.510 nodes at 1.6 units each (argument set-up 72, `mul.s` 12,743 with call, the child's 54
forms by 270 modular additions at 60, gamma 480 (8 modular additions at 60), its six pre-sums mod p with Ea and Ed by
`root.s`'s per-base sequence on gamma (1,753 executed) with two Montgomery reductions and their constant reloads
(2 x 188) and the offset p (2 x 12), stores and control 272: 31,920 <= 32,000) plus the operating system's zero-fill
of its 2^62.855 bytes (0.15625 instructions per byte, below 2^45.9 units), the M = 2^20 walk starts (<= 2^10.1 units
each) and the counters and constant pools: below 2^52.210 units in all.

**As in the sibling package** (proof-impl.md): 4,000 instructions per report (at most 2^55, 2^52.68 units); a fresh
start after every report and at every walk start or restart, at most 2^10.1 units plus 40 of meter bookkeeping
(below 2^65.11 units for the at most 2^55 + 2^20 + 2^20); switches 600 instructions (at most n_max/2^63 + 2^20 + 2^20 + 1,
below 2^57.651 units); the stopping allowance: heapsort of 2^55 records (7,188 instructions per record), dedup, at most 2^4
pairs of two replays (each below (2 x 28/theta + 3) x 2^26 instructions, below 2^103.52 units), the output [x]Q0 (12,809
units), the verification pairing (2^14), halting: below 2^108.52 units. The replay bound uses that one iteration
repeats at most one step (at most 16,380 instructions plus 64 for the exponent update) and one escape (at most
1,076,480), below 2 x 2^25.

**Ledger** (`experiment.py` part D, 60-digit decimal arithmetic):

`W = c S (1 + eps) + c (2^20 + 2^20) 28/theta + 2^108.5196 (non-walk: setup, reports, fresh starts, restarts, stopping allowance) + 2^57.6510 (switches)`

with c = 34,132/120,000, log2 S = 125.67559, eps = 0.0237922, (2^20 + 2^20) 28/theta = 2^111.807: **log2 W = 123.89581**.
Submitted bound: **time_log2 = 123.896** (the least thousandth at or above it; `experiment.py` checks
W <= 2^123.896 exactly on integers, and that the walk-phase threshold covers c (S(1+eps) + window) plus the walk's
non-step work, with 0.0134% to spare).

Sensitivities (part D): eps without the 1.25 factor: -0.0067 bits; +100 instructions on every step: +0.0251 bits;
a lazy kernel twice as expensive: +0.4907 bits.

**Memory:** report log 2^55 x 80 bytes = 2^61.322; trie 2^62.855 (leaves 384^6 x 2,592 bytes = 2^62.850, internal
nodes 2^54.723); table, scalar tables, walk states (2^20 x 4 KB) and replay states below 2^33: **2^63.284 bytes in
all**. Memory is reported, not scored; what binds is addressability, and every array stays below 2^63, so every
byte offset inside an array is below 2^63 and fits one unsigned 64-bit register.

## 6. Heuristics and limitations

* **H1** (random mapping on classes with the conservative R/(R-1) correction), **H2** (the index, twist, DP and
  fingerprint bits of K are pseudorandom functions of the class, modelled as independent uniform draws, which sets
  the per-step rates of section 3.3), **H1+** (the declared 0.01 allowance for the work threshold ending the walk
  phase early; the H2 per-step rates hold over the run as a whole, not only in expectation, section 3.4) and the
  ideal-source model for the ChaCha20 sampler. These are models, not proved distributional facts; the vocabulary of
  independence below is H2's statement about the bits of a class, or about sums of roots of unity, never a claim that
  iterations of the walk or events on different trails are independent. The three caps' losses are Markov bounds on
  integer counts, which hold whatever the dependence; the 0.01 of H1+, the one declared allowance, covers the
  threshold only. Standard in rho analyses (Teske 2001; van Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999;
  Bos-Kleinjung-Lenstra 2010) and supported here by toy and production-field runs, not proved. Lazy evaluation adds
  no heuristic: it computes the same walk (section 3.1).
* H2 for K's bits, stated precisely: the map class -> (bits 0..277 of K) is modelled as a random function (each class
  value uniform, values of distinct classes independent). This is a model, not a distributional fact: the 278 bits
  take at most N + 1 < 2^251.28 values, so no statement about their joint distribution over GT is claimed; the rule
  uses h (the top 5 bits of the low word), 86 DP bits and 128 fingerprint bits, 278 in all. Support:
  K is the minimum of the twelve stored keys, i.e. of six values min(v, p - v) of nonzero Fp-linear functions of the
  element composed with the Montgomery map, so K < p/2 and its density on [0, p/2) varies only on scale p, far above
  the 2^278 range of the bits read; part G tests the bits on 3,000 production classes (chi-square of h 39.06 on 31
  degrees of freedom, of 8 DP-field bits 278.85 on 255, of 8 fingerprint bits 230.04 on 255; largest |z| of the 278 bit
  frequencies 3.43) and part C's walk statistics. The toy censuses use a hash in place of these bits.
* R = 32 is small for an adding walk; the conservative correction 1/(1-1/R) (0.0229 bits on S) and the recomputed
  fruitless-cycle charge (log2(1 + eps) = 0.0339 bits) are both in the ledger. P_model = 0.4005 leaves 0.00041 above 0.39 after the
  computed losses, of which 0.01 is the declared H1+ allowance, below 3.27e-5 the restart cap's worst case and below
  5.07e-5 the three Markov shares; the rest is below 1.1e-6 (both asserted in part D). P_model is set against the census's own model point: the census of
  section 4 measures the fraction of runs whose useful collision falls within S(0.4005)(1 + eps), and its pooled
  value less four binomial standard errors, 0.4047 (floored) over its eight registered rows, is above 0.39041
  (evidence.md). Larger R: the per-array addressability bound, not scored memory, is
  what limits it (section 5); R = 32 leaves the leaf array at 2^62.85 of that 2^63.
* Scope of the prices: kernel counts are executed counts from the shipped emulator; the every-step 170, the escape
  24,000, the flush overhead, the node price, the 60 per modular addition, the 4,000 per report, the heapsort level price and the replay bound are
  analytic bounds. The kernels are checked on real GT elements and edge values (`kernels/check_lazy.py`), not compiled
  and run on hardware. The trie (2^62.855 bytes) and the log (2^61.322 bytes) are RAM-model arrays; memory is unscored
  and no physical machine maps that much. Every array's byte offsets stay below 2^63, so each fits one 64-bit
  register (section 5); that addressability bound, not a scored memory budget, is what R and L must satisfy.
* This is a generic square-root attack (no index calculus in Fp12); a valid upper bound, not a security estimate.

## 7. Experiments

`experiment.py` (standard library only; parts A, S, K, C, G, B, E, F (in `cycles.py`), R (in `reduction.py`), D) and its results, the toy censuses and the
kernel checks are described in `evidence.md`.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999); Teske, Math. Comp. 70 (2001); Cheon, Hong
& Kim, "Speeding up the Pollard rho method on prime fields" (ASIACRYPT 2008) and "Accelerating Pollard's rho
algorithm on finite fields" (J. Cryptology 25, 2012) for tag tracing; Granger & Scott, PKC 2010; Bos, Kleinjung &
Lenstra, ANTS 2010; the sibling package gt-lazy-trie-rho (a77744bb), its parent gt-dprestart-replay (0b8665ae) and their
ancestors; this track's earlier packages and the historical fixed-G2 source.

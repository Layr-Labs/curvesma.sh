# Frobenius-class rho in GT with lazily evaluated steps: a prefix-product trie, a 3,419-instruction norm-paired key kernel and one full multiplication per six steps

**Claim.** For uniform nonidentity P in G1 and uniform z in GT the algorithm below returns Q in G2 with
e(P,Q) = z with probability at least 0.39, using at most **2^123.826** dense-Fp12-multiplication equivalents
(C = 20,000 instructions), everything included: setup (the multiplier table and the prefix-product trie), the two
pairings, walk and fresh starts, distinguished-point (DP) handling, fruitless-cycle escapes, the DP detection delay,
matching (sort, duplicate removal, replays, pair checks), output and verification. Peak memory is at most 2^63.284
bytes, every array below 2^63 so every byte offset fits one unsigned 64-bit register; preprocessing (included) at most
2^52.21 units; nonuniform advice 0 bytes. The success bound rests on the random-mapping heuristic for rho walks (H1),
pseudorandom bits of the class key (H2) and H1+, a declared tolerance on the realized fruitless-cycle waste rate
within which the work threshold costs at most 0.01 (section 3.4); each of the three caps is a failed run charged by
Markov's inequality on its expectation; no concentration bound or independence across iterations or trails is used.
Success is computed in the ideal model in which the ChaCha20 stream that supplies the table and start scalars is an
ideal random source (a disclosed heuristic; no distributional bound on the realized scalars is claimed). No
production-size run is claimed.

## 0. What is new, and what is reused

**This revision** (over our submission 55afb7e5, claim 2^123.896, whose text this package keeps): the key kernels
take the **norm-paired form** of F5. On GT, y y^(p^4) = y^(p^2) makes one quadratic form of the base linear, so
re-pairing Winograd's three Fp2 products along it leaves no per-base product: `root.s` 3,074 -> 1,619 and `lazy.s`
3,457 -> 3,419 executed instructions, a block 34,132 -> 32,487, claim 123.896 -> 123.826. The keys are the same
residues (F5; part C's digest of the 240 steps' (K, kappa) equals 55afb7e5's at the declared input), so the walk,
class bits, records and replays are unchanged; the node keeps 20 key vectors (section 5). The trail tail is now
counted in first visits to classes: the stuck-walk counter credits each fruitless cycle's revisits at its escape
(proof-impl.md 0-i, 0-iii), which sets the restart expectation and the detection-failure term of section 4.

This package revises this track's Frobenius-class rho package (submission 6777f95f, claim 2^125.001: a cyclotomic-squaring
branch at s = 37/64, R = 1,024, class bits from Re(c_0), model success 0.4001 with a 0.01 allowance for H1+). Kept:
the reduction (section 1, `reduction.py`; output [x]Q0 at its counted hard maximum of 12,809 units) and the kernels
`mul.s` (12,735 executed instructions; a table entry is stored as 54 table-side linear forms reduced mod p) and
`sqr.s` (6,786), from `kernels/gen.py`. Replaced: the walk, by the lazily evaluated multiply-only walk of the package
`gt-lazy-trie-rho` (submission a77744bb, claim 2^124.183, accepted on the sibling track bls12-381-fixed-g2-invert-g1,
which solves the same discrete logarithm in the same GT; "the sibling package"): its prefix-product trie, its key
kernels `lazy.s` and `root.s` (`kernels/gen_lazy.py`; rewritten here from 4,964 and 1,561 to 3,419 and 1,619 executed
instructions by Winograd's inner-product identity in its norm-paired form, F5), its class-bit rule, its
fruitless-cycle analysis (`cycles.py`) and its implementation sections 0-i to 0-iv (`proof-impl.md`), with this
package's numbers. Its success chain is not used: section 4 replaces it with H1+ (no concentration bound, one
declared 0.01 for the work threshold, each cap a failed run charged by Markov's inequality).

**Lazy evaluation (tag tracing).** An eager multiply step computes y * T[h][t] (12,735 instructions) only to read, at
the next step, the w^1 coefficient of the result. Here a walk keeps an explicit **base** y and a node of a
**prefix-product trie**: its current element is y * Pi_nu, Pi_nu the product of the at most five twisted multipliers
taken since the base was last made explicit. That coefficient is an Fp2-linear form in y (F5), so the key costs **3
Fp2 products by Winograd's identity (9 products of 384-bit integers): `kernels/lazy.s`, 3,419 executed
instructions**, against 12,735 + 1,910 for the previous package's step; every sixth step **one** multiplication
y <- y * Pi_leaf makes the base explicit. Element by element the walk is the eager multiply-only walk with the same
table (section 3.1). Consequences, each charged: no squaring branch (s = 0: every step's price is fixed, so the
previous package's price-drift margin is dropped; model success 0.4005); **R = 32** (trie of depth 6, 384^6 = 2^51.51
leaves; correction 1/(1 - 1/R); waste charged as eps < 0.0245186, section 3.3); a block of six steps costs 32,487
instructions, **c = 32,487/120,000 = 0.270725 units per step** (previous package 0.51922578125 on average); setup
2^52.21 units (1.6 per node) and memory 2^63.284 bytes (unscored; every array below 2^63); class bits from the key
minimum K: h from bits 0..63, the DP predicate on bits 64..149, the fingerprint on bits 150..277 (section 6).

Ledger: log2 W = 123.82557 (`experiment.py` part D), claim **123.826**, the least thousandth above it, checked by an
exact integer inequality.

## 0a. Representation: which limbs the rule reads

An Fp value is held as the six-limb Montgomery vector x~ = x 2^384 mod p, fully reduced in [0, p), and the rule is
stated on these stored limbs: the twelve keys are the stored residues of F2 (linear operations commute with the
Montgomery map, so stored keys are the Montgomery images of the integer keys and their multiset is again
Frobenius-invariant), K is their minimum as 381-bit integers, kappa the first index (in the kernel's fixed order)
attaining it. This is the integer rule composed with a fixed bijection of Z/p, so every structural property carries
over. `lazy.s` and `root.s` return exactly the stored a, b, B, a zeta, B zeta, K of the current element (parts K and
C check this against a big-integer reference, including keys from explicit Frobenius maps).

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

**F5 Lazy key (new: norm-paired).** For x, Pi in Fp12, the w^1 coefficient of x * Pi is
`c_1(x Pi) = x_0 Pi_1 + x_1 Pi_0 + xi (x_2 Pi_5 + x_3 Pi_4 + x_4 Pi_3 + x_5 Pi_2) = sum_{k<6} x_k gamma_k`,
gamma_0 = Pi_1, gamma_1 = Pi_0, gamma_k = xi Pi_{7-k} (k = 2..5), because w^6 = xi and only the products
w^i w^j with i + j = 1 or 7 contribute to w^1. So a = Re c_1 and b = Im c_1 of the current element y * Pi_nu are
Fp-bilinear in (y, gamma(Pi_nu)).
*The identity.* r divides p^4 - p^2 + 1, so x x^(p^4) = x^(p^2) for x in GT. The p^2-power Frobenius fixes Fp2 and
maps w to zeta w, the p^4-power to zeta^2 w (zeta^3 = -1, zeta^2 + zeta^4 = -1), so the w^3 coefficients (pairs
i + j = 3, and 9 with w^6 = xi) give, in Fp2, **W(x) = 2 x_0 x_3 - x_1 x_2 - xi x_4 x_5 = -x_3 on GT** (part A checks
it and x x^(p^4) = x^(p^2) on GT elements, and that W(x) = -x_3 fails off GT).
*The pairing.* With per-base U0 = -2 x_0, U4 = xi x_4 and node constants A = gamma_3 - 1, B = -gamma_0/2, C = gamma_2,
D = gamma_1, E = gamma_5, F = gamma_4/xi = Pi_3, expanding gives
`(U0 + A)(x_3 + B) + (x_1 + C)(x_2 + D) + (U4 + E)(x_5 + F) = c_1 - (W(x) + x_3) + eta`, eta = AB + CD + EF,
so on GT **c_1 = (sum of the three products) - eta** with no term of the base alone (the previous pairs (0, 1),
(2, 3), (4, 5) left X(x) = x_0 x_1 + x_2 x_3 + x_4 x_5, 9 products per base in `root.s`; those monomials sit in the
w^1, w^5 and w^3 coefficients, each with other monomials, so no such identity removes X). With A' = (U0 + A, x_1 + C,
U4 + E) and B' = (x_3 + B, x_2 + D, x_5 + F), the three Fp2 products are Karatsuba products,
`S0 = sum Re A' Re B', S1 = sum Im A' Im B', S2 = sum (Re A' + Im A')(Re B' + Im B')`, nine products of 384-bit
integers. Every operand is the sum of two stored values in [0, p) (coordinates, U0, U4 and the six operand pre-sums
re + im mod p that `root.s` writes per base, and the node's constants and their pre-sums), so at most 2p - 2 < 2^384,
and S0, S1, S2 lie in [0, 3 (2p - 2)^2]. Per node the trie stores Ea and Ed in [2p, 3p), Ea = 2p + (-(e0 - e1) 2^-384
mod p) and Ed = 2p + (-(2 e0 - e2) 2^-384 mod p), e0, e1, e2 eta's Karatsuba sums on the stored constants. `lazy.s`
forms a_T = S0 - S1 + Ea 2^384 and d_T = 2 S0 - S2 + Ed 2^384, congruent mod p to a 2^384 and d 2^384 (d = a - b),
in [0.78 pR, 4.22 pR] and [0.78 pR, 5.44 pR] with R = 2^384 (the offset 2p is needed: 3 (2p - 2)^2 > 1.21 pR), so in
[0, 2^768) (2^768 > 9.84 pR); one Montgomery reduction of each gives r <= T/2^384 + p < 6.44p, inside the q-table
subtraction of `gen.py` (r < 8p); `gen_lazy.py` asserts these ranges on integers and `check_lazy.py` runs both ends.
Then b = a - d mod p, and since Im d_1 = -Re d_1 (asserted), B = Re(d_1) d and B zeta = (Re(d_1) zeta) d: one
constant product each, plus a zeta; the keys and K follow F2.
*Scope (GT only).* Off GT `lazy.s` returns the same integer function of its stored inputs, which is the key only on
GT: `check_lazy.py` compares `lazy.s` with that function on every case (random and extreme stored inputs included)
and the function with the definition on GT, and part C checks W(y) = -y_3 for the base at every production step.
**Every base `lazy.s` reads lies in GT**: starts g^a z^b (g = e(P, Q0); z in GT is the input, its z^r = 1 checked in setup), the multipliers
M_h = g^a z^b and their Frobenius images T[h][t], so every Pi_nu, the flushed bases y Pi_leaf and escape flushes,
and the escape squares; GT is a group closed under Frobenius, and `mul.s`, `sqr.s` return exact products and
squares (part K).

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
(A..F of F5 as re, im and their sum mod p, then Ea and Ed) and the 54-vector `mul.s` entry of Pi_nu
(for escapes): 74 x 48 = 3,552 bytes;
a leaf (depth 6) stores the `mul.s` entry of Pi_nu: 2,592 bytes. Internal nodes: sum_{d=1}^{5} 384^d = 2^42.929;
leaves 384^6 = 2^51.510; 2^62.855 bytes in all (every array below 2^63). A node is built from its parent by one `mul.s` call
Pi_child = T[h][t] * Pi_parent (the entry's 12 coordinates against the parent's stored forms), then the 54 forms
of its own entry and its key entry (31,985 instructions, charged 32,000 = 1.6 units, section 5).

**Walk state.** Two base buffers (the current explicit base y with its per-base data, U0, U4 and the six
operand pre-sums of F5, and the previous base), the current node and its depth, a ring of three K slots (the current K and
the two before), the previous step's node and buffer index (for escapes), the checkpoint K, the step count mod 64,
the 2-word stuck-walk counter and its 2-word mark, and the trail's start pair (alpha0, beta0). It carries no running exponents.

**One iteration** (pre-debit charge first, proof-impl.md 0-i):

```
1  if node = root: (K, kappa) <- root.s(y)            [also writes y's per-base data: U0, U4, pre-sums]
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

**3.1 The step is a map on classes, and the lazy walk is the eager walk.** The current element is y * Pi_node
exactly (it holds at a start, an escape and after a flush, where node = root, and a step multiplies Pi_node by
T[h][t] on the right). By F5 the key computed in step 1 is that element's, so (K, kappa), h, t, the DP test and the
fingerprint are the current element's; the next element is (y Pi_node) T[h][t] = (y^ M_h)^{p^{-kappa}}, y^ the
canonical representative, in the class [y^ M_h]. So the step depends only on the class, and the sequence of elements
is exactly that of an eager multiply-only walk with this table. Part C runs 240 production steps on the emulated
kernels and asserts at every step that y * Pi_node equals the eager element and (K, kappa) the reference.

**3.2 Collision.** After the walk phase the log is heapsorted, duplicates dropped, and every pair inside a run of
equal F is replayed and checked, under a hard budget of 2^4 pairs (proof-impl.md 0-ii). Equal F means the same class
up to a false match of 2^-128 per pair of distinct classes under H2 (below 2^-19 in expectation for 2^55 records,
with the degenerate pairs beta^1 = beta^2, below 2^-254 each; section 3.4). The replays give canonical pairs with
alpha^1 + x beta^1 = alpha^2 + x beta^2 mod r, so x = (alpha^1 - alpha^2)/(beta^2 - beta^1) unless beta^1 = beta^2
(about 1/r). Output Q = [x]Q0 after checking g^x = z.

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
  merged trails stay merged. A walk whose stuck-walk counter reaches 28/theta = 28 * 2^86 is restarted; the counter
  adds 1 per iteration and each escape credits its episode's revisits, so it counts first visits to classes; it is
  reset at every walk start, restart, report and fresh start (proof-impl.md 0-i, 0-iii).
* **Useful iterations, executed iterations and the waste X.** Every executed iteration computes one key and is
  charged one step price; escapes are charged separately. An executed iteration is *useful* unless it is wasted in a
  fruitless cycle or lies on a trail the stuck rule cuts; S of section 4 counts useful iterations, and H1's collision
  law is over them. The waste X, in steps of c C = 32,487/6 instructions, charges each 2-cycle 2 steps (the
  iteration that re-enters it and the detecting one) plus its escape, 24,000 instructions, below 4.43255 steps; and each
  cycle of length 3 <= l <= 64 at most 64 + l steps before its detection plus its checkpoint escape, which the meter
  charges 24,000 + 16,407 l instructions (the traversal at 16,342 and 65 of bookkeeping per step; 1,074,048 at
  l = 64). Under H2 the expected X per executed iteration is below (2 + 4.43255)/384 < 0.0167515 for 2-cycles plus
  sum_l P_l (64 + l + (24,000 + 16,407 l)/(c C)) for the rest, below 0.0196164 in all; at most 0.0075408 of the
  executed iterations are wasted, so per useful iteration it is below **0.0197654 (the H2 rate)**. The ledger
  charges **eps = 1.25 (eps_2 + eps_l) < 0.0245186** (eps_2 < 0.01675142 for 2-cycles, eps_l < 0.00286344 for the
  rest with the traversal's bookkeeping left out), above 1.2404 times the H2 rate (part D computes both on exact
  fractions and asserts eps <= 0.024519 and eps > 1.2404 times it).
  Useless walks longer than 64 with no distinguished class are left to the stuck-walk rule (charged rate 2^-121, 0-iii).

**3.4 Caps and the work cap.** Setup (at most 2^52.21 units), the walk phase, and the stopping
phase (matching, output, verification, halting; the stopping allowance R_stop < 2^108.52 units, section 5). The
walk phase has the budget B_walk = W_cap C - 2^52.21 C - R_stop C instructions (W_cap = 2^123.826, W_cap C rounded
down) and meters its own charged work in the 192-bit meter under the pre-debit rule (proof-impl.md 0-i) with the
prices of section 5: before each charged action the meter is compared with L = B_walk - A_max (A_max = 2^25) and the
action runs only if the meter is below L. So every run does at most 2^52.21 C + B_walk + R_stop C <= W_cap C
instructions, a deterministic inequality `experiment.py` checks on integers. The step price is fixed per step type
and known before the step: root 1,789, lazy 3,589, lazy with flush 16,342 (the step whose child is a leaf).

*(a) The block lemma: no price-drift heuristic.* A block is the run of iterations from a root iteration to the
next return to the root (a flush, an escape or a fresh start). Its iterations are charged 1,789 (the first), 3,589
(the second to fifth) and 16,342 (the sixth, the only one whose child can be a leaf), so the first j iterations of
any block are charged at most j c C: 1,789, 5,378, 8,967, 12,556, 16,145 and 32,487 for j = 1..6, against
j x 32,487/6, whether or not the block's last iteration ended in a report or an escape (part D asserts each). Hence
n executed iterations are charged at most n c C in steps, whatever the walk does; escapes, reports, fresh starts,
restarts and switches are separate actions. An escape-ended block of k < 6 iterations is charged less than k c C,
by at most 10,927.5 instructions (k = 5), so its escape still costs at least 13,072.5 above k c C (part D prints
the table for every k). The accounting below never spends that saving: X counts each escape at its full charge,
once, as the meter charges it once, so the budget is spent on work, U + X steps for U useful iterations, not on a
count of iterations.

*(b) The threshold inequality.* Let S' = ceil(S). The window after the S'-th useful iteration (both colliding
trails reaching their distinguished classes) takes at most M (28/theta + K) executed iterations outside X: each trail
reaches its class within 28/theta of its own iterations outside X or is cut (credited iterations are revisits in
fruitless cycles, which X counts; proof-impl.md 0-iii), and round-robin adds at most one chunk K = 2^63 per walk. The
cut trails take at most 2^20 x 28/theta iterations outside X (the restart cap). If X up to the window's end is at
most eps S', the meter's total charge up to then is at most c C (S' (1 + eps) + E + M K) plus the walk's non-step
work at its caps (2^55 reports, 2^55 + 2^20 + 2^20 fresh starts and restarts, the switches), E = (2^20 + 2^20)
28/theta < 2^111.81 and M K = 2^83. Part D asserts on exact fractions that this total is below L (0.0306% of
c C S' to spare), so the threshold cannot end the walk phase before that window ends. The escapes inside the window
are part of X; no deterministic bound on their number is claimed.

*Iteration maximum.* Counting also the action that ends a block, every block of k iterations is charged at least
k c C: a block ended by a flush is charged exactly 6 c C; one ended by an escape at least 1,789 + (k-1) 3,589 + 24,000
>= k c C (k <= 6); one ended by a report at least 1,789 + (k-1) 3,589 + 4,000 plus its fresh start, charged its fixed
allowance of 2^10.1 units (proof-impl.md 0-i; its two 255-bit exponentiations alone exceed 2^21 instructions)
(`experiment.py` asserts all three for every k). So a run has at most n_max = B_walk/(c C) + 6M < **2^125.712
iterations** (the last, incomplete block of each walk adds at most 6).

*(c) H1+, the disclosed success heuristic (no drift of the waste rate).* The cost bound does not rest on H1+: the
meter enforces W_cap on every run, so drift moves success only. If X up to the window of the U-th useful iteration
is at most f U, (b) with f for eps lets the walk phase reach that window's end while U (1 + f) <= S (1 + eps): the
budget buys at least floor(S (1 + eps)/(1 + f)) useful iterations and their window, more than S at the H2 rate. Under the section 4 collision law (`experiment.py` part D; success rounded down,
loss up):

| waste rate f / H2 rate (0.0197654) | success at the useful iterations bought | loss from 0.4005 |
| ---: | ---: | ---: |
| 1 | 0.4005 | 0 |
| 1.25 | 0.4003 | 0.0002 |
| 1.5 | 0.3974 | 0.0031 |
| 2 | 0.3916 | 0.0089 |
| 2.0991 (f = 0.04149, the tolerance) | 0.3905 | 0.0100 |
| 2.5 | 0.3859 | 0.0146 |
| 3 | 0.3803 | 0.0202 |

H1+ states that over the walk phase the 2-cycle and checkpoint-escape frequencies, and so f, stay within the
tolerance **f <= 0.04149**, above 2.09 times the H2 rate (and 2.11 times eps/1.25 < 0.0196149; asserted), and
that the H2 per-step rates hold over the run as a whole (which the caps' expectations use). The 0.01 of section 4
charges the threshold ending the walk phase before the S-th useful iteration's window, which by (b) needs f > eps
(above 1.24 times the H2 rate): at the tolerance it is the loss in the table. This is a statement about one realized
frequency over 2^125.7 iterations, much weaker than independence of the steps; no concentration bound is taken.
Support, not proof: parts C, G, E and F (`cycles.py`) and the toy census, whose realized waste over its six rows is
at most 1.03 to 1.19 times the H2 rate and at most 0.57 of the tolerance (evidence.md); H1+ extrapolates the rate
as H1 extrapolates the collision law. The caps are not charged this way: each is a nonnegative integer count with an
expectation bounded in part D, so **Markov's inequality** gives P(cap reached) <= expectation / cap whatever the
dependence, its own loss line (section 4). Reaching the report cap (2^55), the restart cap (2^20) or the pair budget
(2^4) **fails the run**.

*Caps.* Reports D <= 2^55 with one-word counter, at most n_max theta + n_max^2/N < 2^39.712 expected (a report is a
distinguished class at a first visit, a new draw at rate theta under H2, or a revisit of a reported class, an
uncredited revisit; proof-impl.md 0-iii), so Markov bounds the cap
by 2^39.712 / 2^55 = 2^-15.288 and the cap sits more than 2^15.288 above the expectation. At most
2^4 matcher pairs. A true pair (two records of one class from different trails) verifies unless beta^1 = beta^2, and
the first verifying pair ends the run, so a run's matcher charges at most one verifying pair, and the only useless
true pair is one with beta^1 = beta^2. The pairs charged before the halting one are therefore false-key pairs
(2^-128 each under H2) or pairs with beta^1 = beta^2 (below 2^-254 each), at most (2^55 choose 2) (2^-128 + 2^-254) < 2^-19 in expectation, so Markov
bounds the pair budget by 2^-19 / 2^4 = 2^-23 and the budget sits more than 2^23 above the expectation. At most 2^20 stuck-walk restarts, expected below 2^4.804 (0-iii: first-visit tails, long useless cycles, uncredited
revisits), so Markov bounds the cap by 2^4.804 / 2^20 = 2^-15.196 and the cap sits more than 2^15.19 above the
expectation. Part D asserts each ratio (markov_ratio_report >= 40,032, markov_ratio_restart >= 37,550,
markov_ratio_pairs > 2^23) and its JSON prints each expectation, cap, ratio and bound.

## 4. Success probability

**H1 (random mapping).** The class sequence behaves like iterates of a random mapping on N classes in which two
distinct classes have the same image with probability (1 - sum q^2)/N, q the multiplier-index probabilities. As in
the sibling package: for two classes with different indices the image equation has about one solution among N classes; for
the same index about 11/12 of one (the twist adds 11 nontrivial chances). We use the conservative sum q^2 = 1/R =
1/32 (exact first-order value 1/(12R) = 1/384), the Brent-Pollard factor R/(R-1) of an R-adding walk:
**N' = N/(1 - 1/32), log2 N' = 251.31793.** Part E reads N times the pair-collision rate at q = 67, R = 32: 1.0447,
0.9786 and 0.9863 on 15,000, 40,000 and 60,000 classes, against 0.9688 (conservative) and 0.9974 (exact; evidence.md).
H1 does not depend on how the keys are computed: the lazy walk visits exactly the eager walk's elements (3.1).

Under H1 the probability that the first S useful iterations contain a class collision is 1 - exp(-S^2/(2N')). With
P_model = 0.4005: **S = sqrt(2 N' ln(1/0.5995)) = 2^125.67559.**

The losses before the first useful collision is detected, each computed in part D:

* **the 0.01 of H1+ (section 3.4 (c))**: the work threshold ending the walk phase before the S-th useful iteration's
  window, which by (b) needs a waste rate f > eps; 0.01 is its loss at the tolerance f = 0.04149. It covers no cap;
* **each cap's Markov bound** (section 3.4), expectation over cap: the report cap at most **2^-15.288**, the restart
  cap at most **2^-15.196**, the pair budget at most **2^-23** (part D: 2.4980e-5, 2.6632e-5 and 1.1921e-7);
* the restart cap's own worst case spent as steps: 2^20 restarts, each with at most 28/theta iterations outside X,
  2^20 x 28/theta < 2^110.81 steps, a fraction below **3.27e-5** of the iteration maximum 2^125.712, charged in full and not
  at its expected count (0-iii);
* detection failure on either of the two colliding trails, at most 2.518e-9 in all. A trail with no uncredited revisit
  is cut only after at least 28/theta - 64 first visits, none distinguished (0-iii). Under H2, G, the number of first
  visits of a trail up to its first distinguished class, is geometric at rate theta (each first visit a new draw), and a
  trail holds at most I_max = 65 x 28/theta classes. Under H1 the colliding step lands on a uniform visited class,
  independent of the draws at first visits after it, so the trail through it is cut that way with probability at most
  I_max P(G >= 28/theta - 64)/E[min(G, 28/theta - 64)] <= 1820 x/(1 - x), x = (1 - theta)^(28/theta - 64) <= e^(-28 + 2^-80):
  below 1.259e-9 per trail (part D asserts it). A trail cut after an uncredited revisit (the colliding trail after it
  merges, any other merge, a long useless cycle) is one of the at most 2^20 restarts, whose iterations outside X the
  line above removes from the useful count;
* a collision inside a single trail, at most 1e-6: under H1 a step meets one of its trail's at most I_max = 65 x 28/theta
  earlier classes with probability at most that count over N, below 2^-28.7 over n_max steps (asserted in part D);
* the sampler's exhaustion allowance under the ideal-source model of 0-iv, 2^-71.99 < 2.14e-22;
* a key tie or a zero key at most 70/p per keyed element over at most 2^125.712, below 2^-240; a start or iterate equal
  to the identity, 1/r each, below 2^-128 in all; beta^1 = beta^2, below 2^-254.

The sum is **0.01008536**, below 0.0100854 (part D computes it at 60 digits and asserts the bound). So the success
after every loss is at least 0.4005 - 0.0100854 = 0.3904146, **floored to 0.39041**: 0.00041 above the required 0.39.
Part D asserts the floored value. The fresh starts of the trails (one per 2^86 steps) are uniform classes, which
collide at rate 1/N >= 1/N', so counting them as useful iterations is conservative.

**No independence.** Markov's inequality on a nonnegative integer count assumes nothing about dependence, so the cap
losses hold for the adaptive, stateful walk; the first-visit tails use only the chain rule over H2's draws; the 0.01
of H1+ is the loss at a declared tolerance on one realized frequency.

Toy support (evidence.md): a registered census of the complete walk in toy fields at this package's charges and credited
counter, six rows on fresh seeds, pools to **0.4221 +- 0.00468 over 11,175 runs** (0.4034 floored at four standard
errors) for a first useful collision within S(0.4005)(1 + eps) of the toy's charged count, an event easier than
production's by under 1.9%. It measures the threshold and waste part of H1+ and H1 at its sizes, not any cap or
production success.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Counts are executed scalar AArch64 instructions as the cost model counts them. The kernels are straight-line
(general-purpose registers only, no branch but the final ret, no call), so the executed count is the emitted count
on every call; `kernels/emu.py` counts every executed instruction (ldp/stp count one), and part K regenerates the
four walk kernels in memory (`node.s`, a setup routine, runs in `check_lazy.py`) and runs them against a big-integer
reference on real GT elements:

| kernel | executed instructions | what it computes |
| --- | ---: | --- |
| `lazy.s` (`gen_lazy.py`) | **3,419** | from (y, its per-base data, the node's A..F, Ea, Ed): 9 + 3 products (norm-paired Winograd over Fp2), 2 + 3 Montgomery reductions, a, b, B, a zeta, B zeta of y Pi, the 12 stored keys, K and kappa |
| `root.s` (`gen_lazy.py`) | **1,619** | the per-base data of y (U0, U4 and six operand pre-sums, each mod p; no product), then B, a zeta, B zeta (3 products, 3 reductions), the keys, K and kappa |
| `mul.s` (`gen.py`) | **12,735** | y * Pi for a stored entry of 54 table-side forms reduced mod p (flush, escapes, trie construction) |
| `sqr.s` (`gen.py`) | **6,786** | Granger-Scott cyclotomic squaring (escapes) |

`lazy.s` by category (the generator's own instruction tags): products 1,716 (12 of 143: nine Winograd, three
constant), Montgomery reductions 780 and their final subtractions 100 (a, d, B, a zeta, B zeta), loads/stores 342,
operand sums, accumulation, the Ea, Ed terms and the doubling of S0 228, keys and minimum 245, prologue,
epilogue and addressing 8. `root.s`: products 429, reductions 468 and 60, loads/stores 194, the forms and pre-sums
216, keys 245, the rest 7.

**Every step** (both kinds), analytic, 170 instructions: K reads and h = umulh(K_0, 32) 3; DP test 4; kappa and the
twist 5; child index 5; depth and flush test 3; node or leaf address 5; previous node for escapes 2; 2-cycle and
checkpoint comparisons, history ring and checkpoint store test 60; loop control and the kernel call 17; margin 53;
192-bit meter 8; 2-word stuck-walk counter 5 (`experiment.py` lists the items). **Flush** 12,753: `mul.s` 12,735,
call 8, buffer swap and node reset 3, margin 7.

| step | instructions | units |
| --- | ---: | ---: |
| root step (node = root) | 1,619 + 170 = **1,789** | 0.08945 |
| lazy step | 3,419 + 170 = **3,589** | 0.17945 |
| lazy step whose child is a leaf (with the flush) | 3,589 + 12,753 = **16,342** | 0.8171 |
| block of six steps | 1,789 + 4 x 3,589 + 16,342 = **32,487** | **c = 32,487/120,000 = 0.270725 per step** |

**Escape** (2-cycle or checkpoint): 24,000 instructions (choose the smaller key 14, address 4, `mul.s` + call 12,743
(none at the root), `sqr.s` + call 6,794, resets 20, meter 8, counter credit 18, margin 4,399); a checkpoint escape adds its traversal of
at most 64 steps at 16,342 and 65 x 64 of bookkeeping (1,074,048 < A_max = 2^25).

**Setup** (preprocessing, at most 2^52.21 units): pairing g (<= 2^14), z^r = 1 (<= 510), the table (32 multipliers
x (two exponentiations <= 1,021 units + 11 Frobenius maps), 64 sampler draws at 14.848 units, the scalar tables
2 x 384 products mod r at 800 instructions and 384 x 270 modular additions for the entries' forms at 60), **the
trie**: 2^42.929 + 2^51.510 nodes at 1.6 units each (argument set-up 72, `mul.s` 12,743 with call, the child's 54
forms by 270 modular additions at 60, the constants A..F 540 (three xi products, A's -1 and two negate-and-halve
-v/2 mod p, 9 modular operations at 60; a halving by a low-bit select of p or 0, a 6-limb add and 6 shifts, 20, then
a negation 18 and loads and stores 6, 44 <= 60), its six pre-sums mod p with Ea and Ed by `node.s` (`gen_lazy.py`,
2,158 executed: eta's 9 products, two Montgomery reductions and the offset 2p; `check_lazy.py` runs it against the
entry), stores and control 272: 31,985 <= 32,000) plus the operating system's zero-fill
of its 2^62.855 bytes (0.15625 instructions per byte, below 2^45.9 units), the M = 2^20 walk starts (<= 2^10.1 units
each) and the counters and constant pools: below 2^52.210 units in all.

**As in the sibling package** (proof-impl.md): 4,000 instructions per report (at most 2^55, 2^52.68 units); a fresh
start at every report, walk start or restart, at most 2^10.1 units plus 40 of bookkeeping (below 2^65.11 units for
the at most 2^55 + 2^20 + 2^20); switches 600 instructions (at most n_max/2^63 + 2^20 + 2^20 + 1,
below 2^57.653 units); the stopping allowance: heapsort of 2^55 records (7,188 instructions per record), dedup, at most 2^4
pairs of two replays (each charged (2 x 28/theta + 3) x 2^26 instructions, below 2^103.52 units, above the replay
bound 2^97.552 of proof-impl.md 0-ii), the output [x]Q0 (12,809 units), the verification pairing (2^14), halting: below
2^108.52 units.

**Ledger** (`experiment.py` part D, 60-digit decimal arithmetic):

`W = c S (1 + eps) + c (2^20 + 2^20) 28/theta + 2^108.5197 (non-walk: setup, reports, fresh starts, restarts, stopping allowance) + 2^57.6523 (switches)`

with c = 32,487/120,000, log2 S = 125.67559, eps < 0.0245186, (2^20 + 2^20) 28/theta < 2^111.808: **log2 W < 123.82557**.
Submitted bound: **time_log2 = 123.826** (the least thousandth at or above it; `experiment.py` checks
W <= 2^123.826 exactly on integers, and the threshold inequality (b) of section 3.4).

Sensitivities (part D): eps without the 1.25 factor: -0.0069 bits; +100 instructions on every step: +0.0264 bits;
a lazy kernel twice as expensive: +0.5069 bits.

**Memory:** report log 2^55 x 80 bytes = 2^61.322; trie 2^62.855 (leaves 384^6 x 2,592 bytes < 2^62.850, internal
nodes < 2^54.723); table, scalar tables, walk states (2^20 x 4 KB) and replay states below 2^33: **2^63.284 bytes in
all**. Memory is reported, not scored; what binds is addressability, and every array stays below 2^63, so every
byte offset inside an array is below 2^63 and fits one unsigned 64-bit register.

## 6. Heuristics and limitations

* **H1** (random mapping on classes with the conservative R/(R-1) correction; also read as: a class the walk reaches
  equals a given earlier class with probability at most 1/N, which bounds uncredited revisits, 0-iii), **H2** (the
  index, twist, DP and fingerprint bits of K are pseudorandom functions of the class, modelled as independent uniform
  draws: this sets the per-step rates of section 3.3 and makes a class's DP bit at its first visit a new draw),
  **H1+** (the realized waste rate stays within the tolerance 0.04149, at which the work threshold costs at most the
  declared 0.01; the H2 per-step rates hold over the whole run, section 3.4) and the ideal-source sampler. These are
  models, not proved facts; "independence" here is H2's statement about the bits of a class, or about sums of roots
  of unity, never a claim that iterations or trails are independent. Standard in rho analyses (Teske 2001; van
  Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999; Bos-Kleinjung-Lenstra 2010), supported here by toy and
  production-field runs, not proved. Lazy evaluation adds no heuristic (section 3.1); the norm-paired key is exact on GT, where every base lies (F5).
* H2 for K's bits, stated precisely: the map class -> (bits 0..277 of K) is modelled as a random function (each class
  value uniform, values of distinct classes independent). This is a model, not a distributional fact: the 278 bits
  take at most N + 1 < 2^251.28 values, so no statement about their joint distribution over GT is claimed; the rule
  uses h (the top 5 bits of the low word), 86 DP bits and 128 fingerprint bits, 278 in all. Support:
  K is the minimum of the twelve stored keys, i.e. of six values min(v, p - v) of nonzero Fp-linear functions of the
  element composed with the Montgomery map, so K < p/2 and its density on [0, p/2) varies only on scale p, far above
  the 2^278 range of the bits read; part G tests the bits on 3,000 production classes (evidence.md) and part C's walk statistics. The toy censuses use a hash in place of these bits.
* R = 32 is small for an adding walk; the correction 1/(1-1/R) (below 0.0230 bits on S) and the fruitless-cycle
  charge (log2(1 + eps) < 0.0350 bits) are both in the ledger. P_model = 0.4005 leaves 0.00041 above 0.39 after the
  computed losses (section 4; part D asserts the Markov shares below 5.18e-5 and the rest below 1.1e-6). Larger R: the per-array addressability bound, not scored memory, is
  what limits it (section 5); R = 32 leaves the leaf array at 2^62.85 of that 2^63.
* Scope of the prices: kernel counts are executed counts from the shipped emulator; the every-step 170, the escape
  24,000, the flush overhead, the node price, the 60 per modular addition or negate-and-halve, the 4,000 per report, the
  heapsort level price and the replay charge are analytic bounds. The kernels are checked on real GT elements and edge values (`kernels/check_lazy.py`), not compiled
  and run on hardware. The trie (2^62.855 bytes) and the log (2^61.322 bytes) are RAM-model arrays; memory is unscored
  and no physical machine maps that much.
* This is a generic square-root attack (no index calculus in Fp12); a valid upper bound, not a security estimate.

## 7. Experiments

`experiment.py` (parts A, S, K, C, G, B, E, F (in `cycles.py`), R (in `reduction.py`), D), the toy censuses and the
kernel checks are described in `evidence.md`.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999); Teske, Math. Comp. 70 (2001); Cheon, Hong
& Kim, "Speeding up the Pollard rho method on prime fields" (ASIACRYPT 2008) and "Accelerating Pollard's rho
algorithm on finite fields" (J. Cryptology 25, 2012) for tag tracing; Granger & Scott, PKC 2010; Bos, Kleinjung &
Lenstra, ANTS 2010; the sibling package gt-lazy-trie-rho (a77744bb), its parent gt-dprestart-replay (0b8665ae) and their
ancestors; this track's earlier packages and the historical fixed-G2 source.

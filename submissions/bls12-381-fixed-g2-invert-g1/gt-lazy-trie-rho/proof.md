# Frobenius-class rho in GT with lazily evaluated steps: a prefix-product trie, a 3,527-instruction key kernel and one full multiplication per six steps

**Claim.** For uniform nonidentity Q in G2 and uniform z in GT the algorithm below recovers P in G1 with
e(P,Q) = z with probability at least 0.39, using at most **2^123.947** dense-Fp12-multiplication equivalents
(C = 20,000 instructions), including setup (the multiplier table and the prefix-product trie), the two pairings,
walk starts and the fresh start of every trail, distinguished-point (DP) handling, fruitless-cycle escapes, the
DP detection delay, matching (sort, duplicate removal, replay of matched trails, pair checks), output and
verification. Peak memory is below 2^62 bytes; preprocessing (included) is at most 2^50.04 units; nonuniform
advice 0 bytes. The success bound rests on the usual random-mapping heuristic for rho walks (H1) and on
pseudorandom bits of the class key (H2); the success probability is computed in the ideal model in which the
ChaCha20 stream that supplies the table and start scalars is an ideal random source (a disclosed heuristic, no
distributional bound on the realized scalars is claimed). No production-size run is claimed.

## 0. What is new, and what is reused

**This revision (the success argument).** This package is the one the next three paragraphs describe (the
Fp6-level table-side multiplication kernel `kernels/mul.s`, 12,774 instructions, on the Winograd-paired lazy key kernel,
claim 2^123.947) with one change: the success argument. Its earlier text rested on heuristic H1+, one 0.01 allowance
for every cap failure, an early threshold and any drift of the walk's frequencies; an official review refused that
allowance as unsupported (the report cap sat only 2.41 times above its expectation, so Markov gave its reciprocal,
not an additive 0.01; an independent-Bernoulli comparison is no dependence-aware bound). H1+, its allowance, its drift table and
that comparison are removed. Sections 3.4 and 4 now use only computed bounds that need no independence across
steps, trails or reports: expectations by deferred decisions under H2, H1 for merges, and linearity, each cap's overflow by Markov
with the cap at least 617 times its expectation, and the closed-walk waste by Freedman's martingale inequality. To
give the caps that headroom the DP density is theta = 2^-80 (was 2^-72; the DP field is bits 64..143 and the
fingerprint bits 144..255, at the same per-step price), the restart cap is 2^29 and the pair budget 2^8; part F's
exact range is raised to m0 = 16 so that its bound on long useless walks, which grows as theta falls, stays below
2^-108 per step. The kernels, the walk, P_model = 0.4002, S and eps are unchanged; the ledger moves only through
the theta- and cap-dependent terms (log2 W = 123.94681, claim 123.947 unchanged), and the success floor is
0.39689. Every figure is recomputed by `experiment.py` part D and `bounds.py` and re-quoted from their output.

**The multiplication kernel revision.** It took the package of the next two paragraphs (claim 2^124.012) with one
change: the multiplication kernel `kernels/mul.s` puts the fixed operand on the table side at all three tower levels,
with every value it needs stored reduced in its entry (section 5), so it executes **12,774** instructions instead of
14,402; an internal trie node keeps its 12 coordinates beside its entry and its construction is re-itemised (22,040
instructions, charged two units).

**The key-kernel revision.** It took `gt-lazy-trie-rho` (submission a77744bb, claim 2^124.183; "the previous revision") with one change: the key kernel forms the w^1 coefficient by Winograd's pairing of its six terms (section F5), 3 Fp2 products (9 products of 384-bit integers) instead of 6 (18), so `lazy.s` executes 3,527 instructions instead of 4,964. The x-only term of the pairing is computed once per base by `root.s` (3,425 instructions instead of 1,561) and stored beside the base's pre-sums; the gamma-only term is computed once per trie node in setup and stored in the node (two more vectors). The walk, its elements, keys, K, kappa and class bits are unchanged (part C runs the same walk step for step), and so are the table, trie shape, escapes, report log, matcher, meter, sampler, caps and H1/H2; only the step prices and the figures computed from them (c, eps, the ledger, the caps' headroom) move, each recomputed by `experiment.py` part D. The text below is the previous revision's with those figures re-quoted; "the parent" in it is `gt-dprestart-replay`.

This package is a changed copy of `gt-dprestart-replay` (submission 0b8665ae, mpjunior92, claim 2^125.111;
"the parent"). The parent's Frobenius classes, 12-key canonicalisation, twisted multiplier table, trail
discipline with exponent replay, report log and matcher, work meter and its pre-debit rule, sampler, restart
cap and schedule are kept; their implementation sections 0-i to 0-iv are in `proof-impl.md` with this
package's numbers. The kernel `sqr.s` (6,777 executed instructions) is the parent's, generated by the unchanged
`kernels/gen.py`; `kernels/mul.s` is revised here (section 5).

**The change: lazy evaluation of the walk (tag tracing).** A multiply step of the parent computes the dense
product y * T[h][t] (14,402 instructions) only to read, at the next step, two Fp coordinates of the result (the
w^1 coefficient, which gives the twelve keys) and its class bits. Here a walk keeps an explicit **base** element
y and a pointer to a node of a **prefix-product trie**: the walk's current element is y * Pi_nu, where Pi_nu is the
product of the (at most five) twisted multipliers taken since the base was last made explicit. The w^1
coefficient of y * Pi_nu is an Fp2-linear form in y whose six coefficients are coordinates of Pi_nu (section F5),
so the key of the current element costs **3 Fp2 products (9 products of 384-bit integers, Karatsuba, after
Winograd's pairing), two Montgomery reductions and the key arithmetic: the shipped `kernels/lazy.s`, 3,527 executed instructions**,
against 14,402 + 2,014 for the parent's multiply step. Every sixth step the base is made explicit by **one** dense
multiplication y <- y * Pi_leaf with a precomputed product of six twisted multipliers. The walk is, element by
element, the same multiply-only Frobenius-class walk as an eager walk with the same table (section 3.1); only the
way its keys are computed changes. Consequences, each charged:

* **No squaring branch** (s = 0): a squaring needs the explicit element, which the lazy walk does not have; the
  pure adding walk also has no price drift (every step's price is a fixed bound, section 3.4).
* **R = 24 multipliers** (parent 1,024), so the trie of depth L = 6 has (12 x 24)^6 = 2^49.02 leaves. The conservative
  collision correction is 1/(1 - 1/R) (log2 N' = 251.33353) and the fruitless-cycle rates are recomputed for R = 24
  (section 3.3, part F): 2-cycles at exactly 1/(12R) = 1/288 per step, all others below 6.2e-5 per step, and closed
  useless walks longer than 64 below 2^-108 per step at theta = 2^-80; their total waste is charged as eps = 0.032933 (1.25 times the model).
* **Prices:** a root step 3,425 + 170 = 3,595 instructions, a lazy step 3,527 + 170 = 3,697, the flush 12,792,
  so a block of six steps costs 34,872 instructions: **c = 0.2906 units per step** (parent 0.55027). Partial
  blocks cost less per step, so c bounds every step deterministically (section 5).
* **Setup** grows to the trie (below 2^50.04 units, each node two charged units) and **memory** to 2^60.37 bytes of trie
  beside the parent's 2^61.322-byte report log.
* **Class bits** are read from the key minimum K (not from Re(c_0), which the lazy walk does not have): the multiplier
  index h = floor(24 K_0 / 2^64) from the whole low word (bits 0..63), the DP predicate on bits 64..143, the 112-bit
  fingerprint on bits 144..255; the rule reads bits 0..255 of K, 256 bits in all (H2 for them: section 6).

Ledger: log2 W = 123.94681 (`experiment.py` part D), claim **123.947**, the least thousandth above it, checked by
an exact integer inequality.

**Revision of the key-kernel revision (the earlier success argument).** The text inherited from a77744bb bounded its
losses with tail bounds for sums of independent terms, which the walk does not have (a revisit reuses its class,
and coalesced trails share their later history); an official review refused it. The next text replaced those bounds
by the H1+ allowance, which this revision removes in turn (above). P_model = 0.4002 and S = 2^125.68269 date from
that step and are kept.

**Revision of submission d46d2b97 (claim 124.206, rejected).** Its review found four defects; all four repairs stand
in the revision this package builds on: (1) section 6 had asserted joint uniformity of 264 bits of K, impossible for
a class function on N < 2^251.28 classes, and the rule now reads 256 bits (fingerprint 112 bits) with H2 stated as the
usual random-function model of those bits, with production-field statistics (part G); (3) the stated lower bound
0.39048 exceeded the computed 0.3904760, and bounds are floored; (4) the production-field test includes
distinguished points, records, fresh starts and the replay of every recorded trail from its start pair alone (part
C). Its repair (2) substituted tail bounds for independent terms; sections 3.4 and 4 replace that with computed bounds
that assume no independence.

## 0a. Representation: which limbs the rule reads

As in the parent, an Fp value is held as the six-limb Montgomery vector x~ = x 2^384 mod p, fully reduced in [0, p),
and the rule is stated on these stored limbs: the twelve keys are the stored residues of section F2 (linear
operations commute with the Montgomery map, so stored keys are the Montgomery images of the integer keys and
their multiset is again Frobenius-invariant), K is their minimum as 381-bit integers, kappa the first index (in the
kernel's fixed order) attaining it. The rule on stored limbs is the integer rule composed with a fixed bijection
of Z/p, so every structural property used below carries over (the parent's section 0a). `lazy.s` and `root.s`
return exactly the stored a, b, B, a zeta, B zeta, K of the current element (part K and part C check this against
a big-integer reference, including keys from explicit Frobenius maps).

## 1. Reduction

Compute g = e(P0, Q) with the pinned convention pairing(Q, P0). Q has prime order r and the pairing is
nondegenerate, so g != 1 generates GT (prime order r) and z = g^x for a unique x in [0, r). Find x, output
P = [x]P0 and check e(P, Q) = z. If z = 1 output infinity. The answer is for the supplied Q.

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

**F5 Lazy key (Winograd's pairing, the key-kernel revision).** For x, Pi in Fp12, the w^1 coefficient of x * Pi is
`c_1(x Pi) = x_0 Pi_1 + x_1 Pi_0 + xi (x_2 Pi_5 + x_3 Pi_4 + x_4 Pi_3 + x_5 Pi_2) = sum_{k<6} x_k gamma_k`,
gamma_0 = Pi_1, gamma_1 = Pi_0, gamma_k = xi Pi_{7-k} (k = 2..5), because w^6 = xi and only the products
w^i w^j with i + j = 1 or 7 contribute to w^1. So a = Re c_1 and b = Im c_1 of the current element y * Pi_nu are
Fp-bilinear in (y, gamma(Pi_nu)). Fp2 is commutative, so Winograd's pairing gives
`c_1 = sum_{i<3} (x_2i + g_2i+1)(x_2i+1 + g_2i) - xi_x - xi_g`, `xi_x = sum_i x_2i x_2i+1`, `xi_g = sum_i g_2i g_2i+1`
(each product's cross terms are x_2i g_2i + x_2i+1 g_2i+1). xi_x depends only on the base and is computed once per
base by `root.s`, stored beside the six pre-sums; xi_g depends only on the node and is stored in it as
C_re = 4p - xi_g.re, C_im = 2p - xi_g.im (xi_g, xi_x stored residues in [0, p)). With u_i = x_2i + g_2i+1,
v_i = x_2i+1 + g_2i and the Karatsuba sums
`S0 = sum_i u_i.re v_i.re, S1 = sum_i u_i.im v_i.im, S2 = sum_i (u_i.re + u_i.im)(v_i.re + v_i.im)`
(every operand an exact sum of stored values: below 2p for S0, S1, and pre-sum + stored g.re + g.im below 4p for
S2; so S0, S1 < 12p^2 and S2 < 48p^2), `lazy.s` reduces **Im = S2 - S0 - S1 + (C_im - xi_x.im) 2^384**, in
[0, 24p^2 + 2p 2^384), and **Re = S0 - S1 + (C_re - xi_x.re) 2^384**, in [0, 12p^2 + 4p 2^384) (C_re - xi_x.re > 2p
and 2p 2^384 > 12p^2). Adding c 2^384 to a reduction's input adds c to its output, and S0, S1, S2 carry the factor
2^768 of products of stored values, so the outputs are the stored b and a. Both inputs are below 2^768, and one
Montgomery reduction gives r < T/2^384 + 1 + p < 6.22p, inside the q-table subtraction of `gen.py` (r < 8p).
`root.s` forms xi_x from x alone in the same way (Re = S0 - S1 + p 2^384 in [0, 3p^2 + p 2^384), Im = S2 - S0 - S1
in [0, 6p^2), each reduced to [0, p)). Products after the first of each sum are accumulated in registers: row j adds
into the 12-limb sum, and the two carries into limb j + 6 are folded into hi(A_5 b_j) <= A_5 - 1, safe because
every operand is below 4p (A_5 < 2^63); the last row's carry is 0 because the sum stays below its bound.
`gen_lazy.py` asserts each bound above. Then B, a zeta, B zeta, the keys and K follow F2.

**F3 Cyclotomic squaring (escapes only).** The Granger-Scott map of the parent (`kernels/sqr.s`, 6,777).

**F4 Orders.** r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 * 2529403 * 52437899 *
254760293^2 (checked); ord_r(2) = (r-1)/16, so 2^t lies in <p> only for t a multiple of about 2^247.3.

## 3. Algorithm

Parameters: **R = 24** table multipliers, **trie depth L = 6**, DP density **theta = 2^-80**, M = 2^20 walks
round-robin on one processor (chunk K = 2^63 steps; proof-impl.md 0-i), report cap 2^55, restart cap 2^29, pair budget 2^8.

**Table.** M_h = g^{a_h} z^{b_h} for h < 24 with fresh uniform a_h, b_h mod r from the sampler (0-iv), and the
twisted entries T[h][t] = M_h^{p^t}, t < 12: 288 entries (each with the 54 entry vectors of `mul.s`), and
the scalar tables a_h p^t, b_h p^t mod r (read only by replays).

**Trie.** A node nu at depth d <= 6 is a word ((h_1,t_1), ..., (h_d,t_d)) of table entries, with
Pi_nu = T[h_1][t_1] * ... * T[h_d][t_d]; the root is the empty word (Pi = 1). Nodes are stored in level order:
the children of node n are 288 n + 1 + 12 h + t. An internal node (depth 1..5) stores its 18-vector gamma(Pi_nu)
(gamma_k.re, gamma_k.im and their sum, k < 6), the constants C_re, C_im of F5, its 12 stored coordinates and the
54-vector `mul.s` entry of Pi_nu (for escapes): 4,128 bytes;
a leaf (depth 6) stores the `mul.s` entry of Pi_nu: 2,592 bytes. Internal nodes: sum_{d=1}^{5} 288^d = 2^40.855;
leaves 288^6 = 2^49.020; 2^60.37 bytes in all. A node is built from its parent by one `mul.s` call
Pi_child = Pi_parent * T[h][t] and its 54 entry vectors, gamma, C_re and C_im (22,040 instructions, charged 40,000 =
two units).

**Walk state.** Two base buffers (the current explicit base y with its six Karatsuba pre-sums
x_k.re + x_k.im and xi_x of F5, and the previous base), the current node and its depth, a ring of three K slots (the current K and
the two before), the previous step's node and buffer index (for escapes), the checkpoint K, the step count mod 64,
the 2-word stuck-walk counter and the trail's start pair (alpha0, beta0). It carries no running exponents.

**One iteration** (pre-debit charge first, proof-impl.md 0-i):

```
1  if node = root: (K, kappa) <- root.s(y)            [also writes the six pre-sums and xi_x of y]
   else:           (K, kappa) <- lazy.s(y, pre-sums, xi_x, node vector)   [K of the element y * Pi_node]
2  if bits 64..143 of K are 0 (distinguished): F <- bits 144..255 of K; append (F, alpha0, beta0) to the report
   log; the trail ends: fresh start y <- g^alpha0 z^beta0 from a new start pair, node <- root, histories reset
3  2-cycle and checkpoint tests on K (section 3.3); on a hit: escape, then next iteration
4  h <- floor(24 K_0 / 2^64), t <- (12 - kappa) mod 12, node <- child(node, h, t)
5  if depth(node) = 6: y_prev <- y, y <- mul.s(y, entry(Pi_node)) (the flush), node <- root
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
that of the parent's multiply-only walk run with this table. The declared run of part C takes 240 production steps (47 root, 193
lazy, 32 flushes, one escape) with every key from the emulated `root.s`/`lazy.s` and every flush and escape from
the emulated `mul.s`/`sqr.s`, and asserts at every step that y * Pi_node equals the eager walk's element and that
(K, kappa) equals the reference; it also checks the exponent invariant at the end.

**3.2 Collision.** As in the parent: every report is appended to the log until the cap; after the walk phase the
log is heapsorted, duplicates dropped, and every pair inside a run of equal F is replayed and checked, under a hard
budget of 2^8 pairs. Two records with equal F are the same class up to a false match of probability 2^-112 per
pair of distinct classes under H2 (at most (2^55)^2/2 x 2^-112 = 2^-3 false pairs in expectation for 2^55 records;
each costs at most two replays inside the pair budget). The replays give canonical
pairs with alpha^1 + x beta^1 = alpha^2 + x beta^2 mod r, so x = (alpha^1 - alpha^2)/(beta^2 - beta^1) unless
beta^1 = beta^2 (about 1/r). Output P = [x]P0 after checking g^x = z.

**3.3 Fruitless cycles at R = 24.** A closed walk of the class graph is *useless* if the canonical exponent pair
returns unchanged. After l steps the pair is Lambda (alpha^_0, beta^_0) + tau with Lambda = 2^s p^m (s squarings, m a
sum of twists) and tau = sum_h c_h (a_h, b_h), c_h in Z[p]. Because (alpha^_0, beta^_0) is uniform and independent of
the table, the walk is useless (up to probability about 1/r) only if Lambda = 1 and every c_h = 0. An escape's
squaring cannot lie in a useless walk (Lambda = 2^s p^m = 1 needs s >= 2^247, F4). For multiply steps, under H2 the
indices are uniform on 24 and each distinct class carries an independent uniform frame on the 12 powers of p (the
step's twist is (12 - kappa) mod 12 with kappa a class label of H2, so a revisited class reuses its label), and c_h
is the sum of the frames of index h's occurrences. For coefficients far below r, c_h = 0 mod r exactly when this sum of 12th roots of unity
vanishes in Z[zeta_12]. So the probability P_l that the l steps from a given step form a useless walk is at most
E[prod_h g(n_h)] over the multinomial counts n_h of the 24 indices, with g(0) = 1, g(1) = 0 and g(n) = q(n), the
exact probability that n independent uniform 12th roots of unity sum to 0 (computed exactly for n <= 16 in part
F: q(2) = 1/12, q(3) = 1/72, q(4) = 0.0191, ..., q(12) = 0.00195), and g(n) = Q(16) = 0.00112, the largest point
probability of a 16-term sum, for n > 16 (adding independent terms cannot raise a point probability). Part F
evaluates the bound exactly for l <= 384 and, beyond, by Poissonisation (P_l <= e sqrt(l) (E_Pois(l/24) g)^24,
monotone on doubling segments), weighted by the probability exp(-theta l) that a closed walk of length l holds no
distinguished class (H2). Results: **P_2 = 1/288 exactly** (same index twice, opposite twist), sum_{l=3}^{64} P_l =
6.17e-5 (P_3 = 2^-15.34, P_4 = 2^-14.76, P_8 = 2^-25.94, P_64 = 2^-108.88; exponents to two decimals), and useless walks longer than 64 at
a rate **below 2^-108.71 per step** at theta = 2^-80 (charged as 2^-108; with m0 = 12 the same bound reads
2^-93.8 at this theta, which is why m0 is 16).

* **Detection and escapes (the parent's rules, on K).** The K of the step before last is remembered; if it equals
  the current K, the two classes form a 2-cycle: the walk makes explicit the representative of the class with the
  smaller K (the current element y * Pi_node, or the previous element, which is Pi of the previous node times the
  base it was computed from: y, or y_prev if a flush came between; one `mul.s` call with that node's entry, or no
  call when that node is the root) and replaces it by its square (`sqr.s`), node
  <- root, and resets the two-step history and the checkpoint; a replay doubles its pair there. A checkpoint K is
  stored every 64 steps and compared at every step; a repeat means a cycle of length at most 64, which is traversed
  once more to its minimum-K class, which is then squared the same way. Both choices depend only on the cycle, so
  merged trails stay merged. Every walk whose stuck-walk counter reaches 20/theta = 20 * 2^80 steps is restarted
  (0-iii); the counter is reset at every walk start, restart, report and fresh start.
* **Charge.** A 2-cycle wastes the step that re-enters the cycle and the detecting step (2 steps) and one escape
  (24,000 instructions, 1.2 units = 4.1294 steps at c): (2 + 4.1294)/288 = 0.021283 step-equivalents per step. A cycle of
  length 3 <= l <= 64 wastes at most 64 + l steps before its detection, its traversal of l steps at the largest
  step price (16,489, 2.8371 c each) and one escape: sum_l P_l (64 + l + 2.8371 l + 4.1294) = 0.005064. At most 64
  of those 64 + l steps are wasted: the checkpoint is stored at most 63 steps after the walk enters the cycle and is
  compared at every step, so the repeat is seen l steps after the store, the detecting step included; the other l
  are the cycle's first pass, fresh iterations already among the S. That over-count, l c C >= 3 x 5,812 = 17,436
  instructions per checkpoint escape, covers its at most 65 x 64 = 4,160 instructions of bookkeeping (section 5),
  with 13,276 to spare at l = 3 (part D asserts 3 c C - 4,160 = 13,276). The ledger
  charges **eps = 1.25 x 0.026346 = 0.032933** on the step count; the factor 1.25 is the margin of the Freedman bound
  of section 3.4 (the census of this package's weights reads pre-collision 2-cycle rates of 0.00330 to 0.00340 per
  step at R = 24, model 1/288 = 0.003472; supporting evidence for the model rates only, cited for no cap).
  Useless walks longer than 64 with no distinguished class are left to the stuck-walk rule at the charged rate
  2^-108: below 2^17.73 expected restarts over the run's at most 2^125.73 iterations (0-iii).

**3.4 Caps and the work cap.** As in the parent: setup (at most 2^50.04 units), the walk phase, and the stopping
phase (matching, output, verification, halting; the stopping bound R_stop < 2^106.04 units, section 5). The
walk phase has the budget B_walk = W_cap C - 2^50.04 C - R_stop C instructions (W_cap = 2^123.947, W_cap C rounded
down) and meters its own charged work in the 192-bit meter under the pre-debit rule (proof-impl.md 0-i) with the
prices of section 5: before each charged action the meter is compared with L = B_walk - A_max (A_max = 2^25) and the
action runs only if the meter is below L. So every run does at most 2^50.04 C + B_walk + R_stop C <= W_cap C
instructions, a deterministic inequality `experiment.py` checks on integers. The step price is fixed per step type
and known before the step: root 3,595, lazy 3,697, lazy with flush 16,489 (the step whose child is a leaf).

*Fixed step prices.* A block is the run of iterations from a root iteration to the next return to the root
(a flush, an escape or a fresh start). Its iterations are charged 3,595 (the first), 3,697 (the second to fifth) and
16,489 (the sixth, the only one whose child can be a leaf), so a block of k <= 6 iterations is charged at most
3,595 + (k-1) 3,697 < k c C for k < 6 and exactly 6 c C for k = 6, whether or not its last iteration ended in a
report or an escape (the escape and the fresh start are charged separately). Hence the charged step work of n
iterations is at most n c whatever the walk does, and the walk budget buys at least (budget - non-step work)/c
iterations (the eps of section 3.3 counts the iterations wasted in fruitless cycles); the only random rates the budget depends on are
those of escapes, reports and restarts, which are charged separately.

*Iteration maximum.* Counting also the action that ends a block, every block of k iterations is charged at least
k c C: a block ended by a flush is charged exactly 6 c C; one ended by an escape at least 3,595 + (k-1) 3,697 + 24,000
>= k c C (k <= 6); one ended by a report at least 3,595 + (k-1) 3,697 + 4,000 plus its fresh start, charged its fixed
2^10.1 units (proof-impl.md 0-i; its two 255-bit exponentiations alone exceed 2^21 instructions)
(`experiment.py` asserts all three for every k). So a run has at most n_max = B_walk/(c C) + 6M < **2^125.73
iterations** (the last, incomplete block of each walk adds at most 6).

*Success bounds that need no independence.* These bounds hold under H1, H2 and the model of section 3.3, not
unconditionally; the meter enforces W_cap on every run whatever happens, so they move success only. No allowance
is charged anywhere, and no step assumes independence of iterations, trails or reports. `bounds.py` (run by part D)
computes and asserts every figure below.

*Deferred decisions.* H2: the labels of distinct classes (index and twist bits, DP bits, fingerprint bits of K) are
independent and uniform. Reveal a class's label when a walk first reads it, in any adaptively chosen order. Lemma:
the class an iteration reaches is a function of the labels revealed before it (the step reads only the current
class's label, section 3.1); if that class was never revealed, its label is uniform given everything revealed
(the sigma-field F_{t-1}), whatever adaptive rule chose it. Proof: conditioning on F_{t-1} fixes that class and
only labels of other classes, which are independent of its label. Such an iteration is *fresh*.

*Expectations, and the caps by Markov.* A fresh iteration is distinguished with probability theta given F_{t-1},
so by linearity at most theta n_max fresh iterations are distinguished in expectation. A non-fresh iteration
reports only after its walk entered an earlier trail's path (a 2-cycle or closed walk of one trail holds no
distinguished class, or the trail would have ended there), and each such merge ends in at most one report; under H1
the expected number of merges in at most n_max iterations is at most n_max^2/(2N') = 0.546 < 1. So
E[D] <= n_max theta + 1 = 2^45.73, and the report cap 2^55 is **617 times** that: P(D reaches the cap) <= 0.00162
(Markov). Apart from a merge into another trail's path and a useless closed walk longer than 64 (both counted
below), a trail is restarted only if its first 19/theta fresh iterations hold no distinguished class (at most
(1 - theta)^(19/theta) <= e^-19, one fresh label at a time) or it wastes more than 1/theta of its 20/theta
iterations in closed walks (below 2^-1000, the per-trail Freedman bound below: in the model a fresh iteration
wastes 0.01112 iterations in expectation, so the 19/theta fresh iterations of a trail are expected to waste below
0.2113/theta iterations, about a fifth of the 1/theta the 20/theta restart leaves for waste, and Freedman's inequality,
not this expectation, bounds the chance of wasting more); there are at most E[D] + M + 2^29 trails in expectation and at most 0.546 merges (above), so expected stuck restarts are below 2^18.32,
and long useless walks add at most n_max 2^-108 = 2^17.73 in expectation: E[restarts] < 2^19.06, and the restart cap **2^29** is
**986 times** that (overflow at most 0.00102). False pairs: at most (2^55)^2/2 x 2^-112 = 2^-3 expected (linearity),
against the pair budget **2^8**, 2,048 times that (overflow at most 2^-11). Each overflow share is a loss of
section 4; a run that hits a cap fails.

*Closed-walk waste: Freedman's inequality.* For a fresh iteration t let X_t be the waste of the closed walk that
starts there (0 if none), as charged in section 3.3: 2 + 4.1294 step-equivalents for a 2-cycle, at most
64 + l + 2.8371 l + 4.1294 for length 3 <= l <= 64, so X_t <= b = 313.71. Iterations that go round a closed walk
before its detection are not fresh and carry no X: their cost is inside the X of the walk's start, so a revisit
inside a closed walk adds nothing beyond the model rate. Let tau be the first iteration that reaches a revealed class
other than inside such a closed walk (the first inter-trail collision, or a long self-collision, below 1e-6; a trail caught in a useless closed walk longer than 64 opens no window and is paid by the restart allotment (0-iii), so it is not tau). Fix a
residue j mod 128; when a walk reaches a fresh iteration whose step count is j mod 128, reveal the labels of its
next 127 classes at once (a closed walk from that iteration visits at most 65 of them). Before tau these windows are
disjoint (one walk's are 128 iterations apart, distinct walks' paths are disjoint), so by the lemma the window's
labels are uniform given everything revealed, X_t is known when it is revealed, and E[X_t | revealed] <= r = 0.026346, the part F model rate (P_2 = 1/288 exactly; the index and
twist of a fresh class are uniform, so the branch frequencies need no separate claim). For each j the centred sum
is a martingale with increments at most b and predictable variance at most b r S over the first S fresh
iterations. Freedman's inequality, P(sum (X_t - E[X_t | F]) >= lam and V <= v) <= exp(-lam^2/(2(v + b lam/3))),
with lam = (eps - r) S/128 per residue (the margin 0.25 x the model rate) and v = b r S, bounds the chance that the
waste of the first S fresh iterations exceeds eps S by 128 exp(-...) < 2^-(2^93) (part D asserts that weaker bound;
the computed exponent is over 1.57e28 bits, far beyond it), charged as 2^-1000. The same inequality, per trail and
per residue, with the same b = 313.71 (an over-count, since a start's waste is bounded by the step equivalents it
is charged), deviation (1/theta - 19 r_it/theta)/128 per residue and v = 128 r_it 19/theta, where r_it < 0.01112
(asserted) is the model's expected wasted iterations per fresh one (above), gives the
2^-1000 of the restart bound above.

*The threshold is not an allowance.* Outside the Freedman event and the cap failures, the first S fresh iterations
and their waste fit S(1 + eps) steps and restarts fit 2^29 20/theta steps; the walk budget buys c (S (1 + eps) +
(M + 2^29) 20/theta) plus the reports, starts, restarts and switches before the threshold (part D asserts it on
integers, with 0.0131% to spare, over 2^111.05 units). The detection window, escapes included, fits in that spare:
after the collision both trails report within 20/theta of their own iterations (or detection fails, section 4)
and round robin moves each walk at most one 2^63 chunk more, so it has at most M (20/theta + 2^63) < 2^104.33
iterations, each at most c of step work (blocks, above) and one escape of at most 1,083,456 instructions: below
2^110.09 units, about half the spare (`bounds.py`). So an early threshold is inside those computed losses. The census
of evidence.md supports H1's collision law and the model rates only and is cited for no cap.

## 4. Success probability

**H1 (random mapping).** The class sequence behaves like iterates of a random mapping on N classes in which two
distinct classes have the same image with probability (1 - sum q^2)/N, q the multiplier-index probabilities. As in
the parent: for two classes with different indices the image equation has about one solution among N classes; for
the same index about 11/12 of one (the twist adds 11 nontrivial chances). We use the conservative sum q^2 = 1/R =
1/24 (exact first-order value 1/(12R) = 1/288), the Brent-Pollard factor R/(R-1) of an R-adding walk:
**N' = N/(1 - 1/24), log2 N' = 251.33353.** Part E samples 15,000 classes at q = 67, R = 24 in the declared run and reads N times
the pair-collision rate 0.9816 +- 0.0342 against 0.9583 (conservative) and 0.9965 (exact); the same run at 20,000 classes
(pair_samples = 20,000) reads 0.9906 +- 0.0255. H1 does not depend on how the keys are computed: section 3.1 shows the lazy walk visits exactly the
elements of the eager multiply-only walk.

Under H1 the probability that the first S fresh classes contain a class collision is 1 - exp(-S^2/(2N')). With
P_model = 0.4002: **S = sqrt(2 N' ln(1/0.5998)) = 2^125.68269.** Losses, each computed in part D and `bounds.py`,
none an allowance: the report cap's Markov share 0.00162 (cap 617 times its expectation), the restart cap's 0.00102
(986 times), the pair budget's 2^-11 = 0.0004883 (2,048 times); the restart cap's own worst case as steps spent on
restarts (2^29 restarts of one full stuck-walk bound 20/theta = 20 x 2^80 each, 2^113.33 steps, a fraction 1.841e-4 of the
iteration maximum 2^125.73, computed exactly and subtracted; 0-iii); the Freedman tail of the closed-walk waste,
2^-1000; detection failure on the collision (below 1.13e-7 for each of the two trails, 2.26e-7 charged: two trails that reach the same class walk the same path from
there and report the same class, so detection fails only if a trail is restarted or the run ends before its report,
at most 20 positions x e^-19 x (1 + theta), plus 2^-1000; part D asserts the figure); a collision inside a single trail (below 1e-6; part D asserts 2 x 20/theta / S < 1e-6 behind it); the
sampler's exhaustion bound 2^-71.99; a key tie or a zero key (at most 70/p per keyed element over at most
2^125.73, below 2^-240); and beta^1 = beta^2 (2^-255). The report log loses no report before its cap (0-ii). The
sum is 0.0033071: **success after every loss at least 0.4002 - 0.0033071 = 0.3968929, floored to 0.39689** (part D
asserts the floor). These hold under H1 and H2, not unconditionally. The value 0.4002 is set against the
pre-registered toy census of evidence.md, whose rule is fixed before any count is read: at each of the three tested
fields the fraction of runs whose first useful collision falls within S(0.4002)(1 + eps), minus 4 binomial
standard errors, must be at least 0.391. Counted in meter units (every iteration before the collision at one step of
c, the step detecting each escape and the distinguished-point steps included, each escape at its largest meter price;
evidence.md states exactly what is counted), all three pooled fields clear it (0.4068, 0.4013 and 0.3985). The
census was first counted in narrower units that left out the detecting and distinguished-point steps and charged a
checkpoint escape at 1,079,296 instructions (0.4186, 0.4134 and 0.4116); its nine registered rows were run again
unchanged, the first count was reproduced cell for cell, and the rule is applied to the stricter meter-unit count.
The two older
tables in evidence.md are marked superseded and nothing here uses them. The fresh starts of the trails
(one per 2^80 steps) are uniform classes, which collide at rate 1/N >= 1/N', so counting them as iterates is
conservative.

Toy support (evidence.md: censuses of this package's complete walk at R = 24, s = 0 with trails, replay and the
matcher in toy fields q = 67, 163 and 211): the meter-unit fraction is 0.4090 to 0.4343 per row and, pooled per
field, 0.4276, 0.4222 and 0.4193 (3,849, 3,800 and 3,774 of 9,000 runs, truncated; all 27,000 runs 0.4230, 11,423,
+- 0.0030), with the mean T_c/sqrt(N) at 1.25 to 1.27 against the model's 1.2803. The count includes the toy's
distinguished-point steps (one per 64; production has one per 2^80), so it is stricter than production's;
evidence.md gives every row and its seeds, in both counts.
The toy derives its class bits from a hash of the canonical representative (a toy key in F_q has too
few bits); the production rule's own bits of K are tested in the real GT by part G. They support H1 at their sizes; they do not measure
production success.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Counts are executed scalar AArch64 instructions as the cost model counts them. The four kernels are straight-line
programs (general-purpose registers only, no branch but the final ret, no call), so the executed count equals the
emitted count on every call; `kernels/emu.py` interprets the subset and counts every executed instruction (ldp/stp
count one), and `experiment.py` part K regenerates all four in memory and runs them against a big-integer
reference on real GT elements:

| kernel | executed instructions | what it computes |
| --- | ---: | --- |
| `lazy.s` (`gen_lazy.py`) | **3,527** | from (y, its pre-sums and xi_x, the node's gamma(Pi), C_re, C_im): 9 products for c_1 (F5) and 4 for B, a zeta, B zeta, 2 + 3 Montgomery reductions, a, b, B, a zeta, B zeta of y Pi, the 12 stored keys, K and kappa |
| `root.s` (`gen_lazy.py`) | **3,425** | the six pre-sums and xi_x of y (9 products, 2 reductions), then B, a zeta, B zeta (4 products, 3 reductions), the keys, K and kappa |
| `mul.s` (`gen.py`, revised here) | **12,774** | y * Pi against a stored entry of 54 reduced linear forms (flush, escapes, trie construction) |
| `sqr.s` (`gen.py`, parent) | **6,777** | Granger-Scott cyclotomic squaring (escapes) |

`lazy.s` by category (from the generator): products 1,931 (the 9 of F5, 3 at 143 and 6 accumulated in registers at
155, and 4 key products at 143), Montgomery reductions 780 and their final subtractions 100 (5 reductions),
loads/stores 300, operand sums and combinations 180, keys and minimum 227
(eleven 6-limb compares with 6 + 1 selects, two 6-limb subtractions with add-back select, six 6-limb negations p - v),
prologue/epilogue and addressing 9.

**`mul.s` (revised here).** The fixed operand m takes the table side at all three tower levels. Fp12 (x = x0 + x1 w, m = m0 + m1 w in Fp6 halves, v = w^2 within Fp6 in blst's tower): c0 = k1 + k2, c1 = k1 + k3 with k1 = (x0 + x1) m0, k2 = x1 (v m1 - m0), k3 = x0 (m1 - m0). Fp6 (each product y b, y = y0 + y1 v + y2 v^2, b fixed): six Fp2 products PA = (y1 + y2) xi b2, PB = (y0 + y2) b0, PC = (y0 + y1) b1, Q0 = y2 (xi b1 - xi b2 - b0), Q1 = y1 (b0 - xi b2 - b1), Q2 = y0 (b2 - b0 - b1), and o0 = Q0 + PA + PB, o1 = Q1 + PA + PC, o2 = Q2 + PB + PC (xi = 1 + u). Fp2 (each product a c, c fixed): p0 = (a0 + a1) c0, p1 = a1 (-(c0 + c1)), p2 = a0 (c1 - c0), re = p0 + p1, im = p0 + p2. That is still 54 products of 384 x 384 bits and 12 Montgomery reductions; every m-side value is stored reduced in [0, p), so the x-side operands are exact sums below 8p, every output's 768-bit integer is a sum of 12 non-negative products (no offset is added, below 45 p^2), and the integer is congruent mod p to the reference product's coordinate, which `gen.py` asserts. The entry does not contain the 12 coordinates of m (two of them appear as vectors), so an internal trie node stores its 12 coordinates beside its entry and builds each child's entry by `gen.py`'s program of 113 modular additions and subtractions (a negation is a subtraction from zero; each input is in [0, p), so one correction by p reduces each result) plus 2 copies; `experiment.py` part K replays that program (`entry.py`) on every entry it uses and checks each vector against its linear form mod p. Flush, escape and trie prices use the new count.

**Every step** (both kinds), analytic, 170 instructions: K reads and h = umulh(K_0, 24) 3; DP test 4; kappa and the
twist 5; child index 5; depth and flush test 3; node or leaf address 5; previous node for escapes 2; 2-cycle and
checkpoint comparisons, history ring and checkpoint store test 60; loop control and the kernel call 17; margin 53;
192-bit meter 8; 2-word stuck-walk counter 5 (`experiment.py` lists the items). **Flush** 12,792: `mul.s` 12,774,
call 8, buffer swap and node reset 3, margin 7.

| step | instructions | units |
| --- | ---: | ---: |
| root step (node = root) | 3,425 + 170 = **3,595** | 0.17975 |
| lazy step | 3,527 + 170 = **3,697** | 0.18485 |
| lazy step whose child is a leaf (with the flush) | 3,697 + 12,792 = **16,489** | 0.82445 |
| block of six steps | 3,595 + 4 x 3,697 + 16,489 = **34,872** | **c = 0.2906 per step** |

**Escape** (2-cycle or checkpoint): 24,000 instructions (choose the smaller key 14, address 4, `mul.s` + call 12,782
(none at the root), `sqr.s` + call 6,785, resets 20, meter 8, margin 4,387); a checkpoint escape adds its traversal
of at most 64 steps at 16,489 and 65 x 64 of bookkeeping (1,083,456 < A_max = 2^25).

**Setup** (preprocessing, at most 2^50.04 units): pairing g (<= 2^14), z^r = 1 (<= 510), the table (24 multipliers
x (two exponentiations <= 1,021 units + 11 Frobenius maps), 48 sampler draws at 14.848 units, the scalar tables
2 x 288 products mod r at 800 instructions and 288 entries at 54 x 200), **the trie**: 2^40.855 + 2^49.020 nodes at
two units each (copy of the parent's 12 stored coordinates 72, `mul.s` 12,782 with call, the child's 54 entry vectors
from the 113-operation program at 45 and 2 copies at 8, 5,101, gamma 480, xi_g at most all of `root.s` 3,425, C_re and
C_im 40, stores and control 140: 22,040 <= 40,000) plus
the operating system's zero-fill of its 2^60.37 bytes (0.15625 instructions per byte, below 2^43.41 units), the M = 2^20
walk starts (<= 2^10.1 units each) and the counters and constant pools: below 2^50.04 units in all.

**Unchanged from the parent** (proof-impl.md): 4,000 instructions per report (at most 2^55, 2^52.68 units); a fresh
start after every report and at every walk start or restart, at most 2^10.1 units plus 40 of meter bookkeeping
(below 2^65.11 units for the at most 2^55 + 2^29 + 2^20); switches 600 instructions (at most n_max/2^63 + 2^20 + 2^29 + 1,
below 2^57.68 units); the stopping bound: heapsort of 2^55 records (7,188 instructions per record plus 7,484), dedup, at most 2^8
pairs of two replays (each below (2 x 20/theta + 3) x 2^26 instructions, below 2^97.04 units), [x]P0, the verification
pairing, halting: below 2^106.04 units. The replay bound uses that one iteration repeats at most one step (at most
16,489 instructions plus 64 for the exponent update) and one escape (at most 1,083,456), below 2 x 2^25.

**Ledger** (`experiment.py` part D, 60-digit decimal arithmetic):

`W = c S (1 + eps) + c (2^20 + 2^29) 20/theta + 2^106.035 (non-walk: setup, reports, fresh starts, restarts, stopping bound) + 2^57.68 (switches)`

with c = 34,872/120,000, log2 S = 125.68269, eps = 0.032933, (2^20 + 2^29) 20/theta < 2^113.33: **log2 W = 123.94681**.
Submitted bound: **time_log2 = 123.947** (the least thousandth at or above it; `experiment.py` checks
W <= 2^123.947 exactly on integers, and that the walk-phase threshold covers c (S(1+eps) + window) plus the walk's
reports, starts, restarts and switches, with 0.0131% to spare, about half of it for the window's escapes, 3.4).

Sensitivities (part D): eps without the 1.25 factor: -0.0093 bits; +100 instructions on every step: +0.0247 bits;
a lazy kernel twice as expensive: +0.50 bits (still below the parent's 125.111).

**Memory:** report log 2^55 x 80 bytes = 2^61.322; trie 2^60.37 (leaves 288^6 x 2,592 bytes = 2^60.36, internal
nodes 2^52.87); table, scalar tables, walk states (2^20 x 4 KB) and replay states below 2^33: **below 2^61.93 bytes**
(< 2^62). Each array is addressed by one 64-bit register (every byte offset < 2^63). Memory is reported, not scored.

## 6. Heuristics and limitations

* H1 (random mapping on classes with the conservative R/(R-1) correction), H2 (pseudorandom index, twist, DP and
  fingerprint bits of K; distinct classes given independent uniform labels, which supplies the fruitless-cycle
  rates and the expected counts the caps are compared with; section 3.4 uses it only through deferred decisions,
  Markov and Freedman, and assumes no independence of iterations, reports, restarts or trails), and the ideal-source model
  for the ChaCha20 sampler. These are the
  parent's heuristics; standard in rho analyses (Teske 2001; van Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999;
  Bos-Kleinjung-Lenstra 2010) and supported here by toy and production-field runs, not proved. Lazy evaluation adds
  no heuristic: it computes the same walk (section 3.1).
* H2 for K's bits, stated precisely: the map class -> (bits 0..255 of K) is modelled as a random function (each class
  value uniform, values of distinct classes independent). This is a model, not a distributional fact: the 256 bits
  take at most N + 1 < 2^251.28 values, so no statement about their joint distribution over GT is claimed; the rule
  uses h (a function of the whole low word, bits 0..63), 80 DP bits and 112 fingerprint bits, 256 in all. Support:
  K is the minimum of the twelve stored keys, i.e. of six values min(v, p - v) of nonzero Fp-linear functions of the
  element composed with the Montgomery map, so K < p/2 and its density on [0, p/2) varies only on scale p, far above
  the 2^256 range of the bits read; part G tests the bits on 3,000 production classes (chi-square of h 21.7 on 23
  degrees of freedom, of 8 DP-field bits 237.2 on 255, of 8 fingerprint bits 220.8 on 255; largest |z| of the 256 bit
  frequencies 2.67) and part C's walk statistics. The toy censuses use a hash in place of these bits.
* R = 24 is small for an adding walk; the conservative correction 1/(1-1/R) (0.0308 bits on S) and the recomputed
  fruitless-cycle charge (below 0.0468 bits at eps = 0.032933) are both in the ledger; P_model = 0.4002 leaves 0.00689 above 0.39 after the
  computed losses (section 4), and the toy censuses read above the model at this
  point (section 4). Larger R makes the trie (12R)^6 too large for memory below 2^62.
* Scope of the prices: kernel counts are executed counts from the shipped emulator; the every-step 170, the escape
  24,000, the flush overhead, the node price, the 4,000 per report, the heapsort level price and the replay bound are
  analytic bounds. The kernels are checked on real GT elements and edge values (`kernels/check_lazy.py`), not compiled
  and run on hardware. The trie (2^60.37 bytes) and the log (2^61.322 bytes) are RAM-model arrays; memory is unscored
  and no physical machine maps that much.
* This is a generic square-root attack (no index calculus in Fp12); a valid upper bound, not a security estimate.

## 7. Experiments

`experiment.py` (standard library only; parts A, S, K, C, G, B, E, F (in `cycles.py`), D) and its results, the toy censuses and the
kernel checks are described in `evidence.md`.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999); Teske, Math. Comp. 70 (2001); Cheon, Hong
& Kim, "Speeding up the Pollard rho method on prime fields" (ASIACRYPT 2008) and "Accelerating Pollard's rho
algorithm on finite fields" (J. Cryptology 25, 2012) for tag tracing; Granger & Scott, PKC 2010; Bos, Kleinjung &
Lenstra, ANTS 2010; parent package gt-dprestart-replay (mpjunior92) and its ancestors (fkiene, kaankolcu,
GordoAR, jungjipdo, mpjunior92 7da5f23a).

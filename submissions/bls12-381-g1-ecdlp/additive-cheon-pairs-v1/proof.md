# Additive Cheon with planned pair walks: P+T and P−T from one inversion

**Claim.** For every n with 2^20 ≤ n < r, the deterministic algorithm below recovers s with heuristic
success probability at least 0.399 (≥ 0.39 required). Its total work is W(n) = W(2^20) ≤ **2^116.580443**
target-group-operation equivalents (exact score 116580443/1000000).
- Preprocessing, meaning the challenge-independent work, is ≤ 2^115.5847 and is included in W.
- Peak memory is ≤ 2^122.2354 bytes. There is no advice.
- The algorithm reads only A_1 = T and A_d, with d = 2^20.

This is the construction of the accepted earlier submission `additive-cheon-v1` (cbadde77: Cheon's (p−1)
attack as a baby-step/giant-step in the subgroup H of d-th powers, baby and giant sets of small integers
a·b ∈ H walked at one batched affine addition per element, the 6-fold symmetry through the key min(y, p−y)),
with one algorithmic change. **Each character-class list is paired offline, in integer arithmetic, into
three-term progressions (x−g, x, x+g).** The walk then computes [a(x+g)]Q = P + T and [a(x−g)]Q = P − T from
P = [ax]Q and the table entry T = [ag]Q with **one shared field inversion**, because both slopes have the
denominator x(T) − x(P). Under the reference convention (batched incomplete affine addition = 6M: 3M for the
Montgomery batch-inversion share, 1M + 1S + 1M for the formula, S = M) a pair of points costs 3M + 3M + 3M =
9M, that is 0.75 reference units per point instead of 1. The planner pairs more than 99 percent of every
toy class list at the production degree parameter; the claim uses 98 percent (H4). Everything else, in
particular the success argument, is as in cbadde77, with the parameters retuned (K = 0.51, X2 = 2^76,
48-byte records), the machine model stated explicitly and the per-round field inversion charged.

## 1. Facts about the fixed target

These are checked by `evidence.py` and by direct computation.
- r ≈ 2^254.857 is prime and v₂(r−1) = 32, so d = 2^20 divides r−1.
- H = {x^d : x ∈ F_r^*} = ker(x ↦ x^m) has order m = (r−1)/d ≈ 2^234.857, and 6 | m.
- **The hidden value is uniform on H.** t := s^d is uniform on H, because s is uniform on F_r^* and
  x ↦ x^d is d-to-1 onto H. A_d = [t]P0 is public.
- **The symmetry group sits inside H.** Let U = {±1, ±λ, ±λ²}, where λ is the GLV eigenvalue
  (λ² + λ + 1 ≡ 0 mod r) with [λ](x, y) = (βx, y). U ⊂ H because its order 6 divides m.
- **One key identifies a whole orbit.** All six points [u]Q, u ∈ U, have y-coordinate ±y(Q), because λ
  fixes y. Conversely, the points with y-coordinate ±y₀ are exactly the six (βʲx₀, ±y₀). So
  canonical(Q) := min(y, p − y) identifies the U-orbit of Q.
- **Membership in H is a character test.** χ(x) := log_ζ(x^m) ∈ Z/d is a character with kernel H, for a
  fixed generator ζ of μ_d. So a·b ∈ H ⇔ χ(a) + χ(b) ≡ 0 (mod d).
- **2 ∉ H.** 2^m ≠ 1 (mod r); 2^m has order 2^16. So χ(2) ≠ 0.

## 2. Algorithm

Parameters: K = 0.51; X2 = 2^76 (the integers b are at most X2); base primes a in the interval
I = (2^76, 2^76 + 2^70]; GMAX = 2^44 (per-prime tables T[g] = [g]Q for 1 ≤ g ≤ GMAX); 2^16 lanes (primes in
flight) per walk; N_B = N_G = N := ⌈√(K·m/6)⌉ + 1 ≈ 2^115.6503, so that 6N² ≥ K·m.

1. **Characters and class lists** (challenge independent).
   - Sieve the primes q ≤ X2 and compute χ(q) for each (one 235-bit exponentiation in F_r and a
     Pohlig–Hellman logarithm in μ_{2^20}).
   - Sieve χ(b) for all b ≤ X2 by multiplicativity, into an array of 20-bit entries (the χ array).
   - Stream each b into its class list L_c = {b ≤ X2 : χ(b) = c}, kept in increasing order.
   - Sieve the primes of the interval I (a segmented sieve with the primes below 2^39) and compute χ(a) for
     each of them. Keep them in a table indexed by their position in I (the prime index).
2. **Planning** (challenge independent, section 3). For each class c, compute a plan: a sequence of entries
   that produces every element of L_c exactly once, each either as a root (one addition from the previous
   root) or as one of the two children x−g, x+g of a pair entry with midpoint x and offset g ≤ GMAX.
3. **Baby walk.** For fresh primes a ∈ I, in order of prime index:
   - set Q_a = [a]P0 and build T[k] = [k]Q_a for 2 ≤ k ≤ GMAX (T[2] by one doubling, T[k] = T[k−1] + Q_a);
   - execute the plan of class −χ(a): the first element [L[0]]Q_a by a scalar multiplication; a root entry
     adds T[gap] to the previous root; a pair entry forms P + T[g] and P − T[g] from the stored midpoint P;
   - for every point produced, insert the record (key, locator) into the hash table, with key = 256 bits of
     canonical(point) and locator = (prime index, index of b in L_c), the index being carried by the plan
     entry.

   Stop after N elements. Every element a·b lies in H. The elements are distinct integers below r: a > X2 ≥ b,
   so a is the unique prime factor of a·b above X2 and (a, b) is determined by the product; a·b < 2^153 < r.
4. **Giant walk.** The same with disjoint primes a′ (the next primes of I) and base A_d, which gives the
   points [a′b′·t]P0. Look up the key of each point. Stop after N elements (no early stop is claimed).
5. **Collision.** A key match means [g·t]P0 = [u·ab]P0 for some u ∈ U, where g = a′b′, up to a false match
   of the 256-bit key. Resolve each candidate pair of locators into (a, b, a′, b′) by reading the prime
   table and the class lists, and test the six values t = u·ab/g by checking [t]P0 = A_d. At most 2^20
   candidates are processed; more is counted as failure.
6. **Recover s.** Take one d-th root s0 of t: twenty square roots, keeping at each step the root that is
   still a 2^j-th power. Then s = s0·ω^i for some i < 2^20, where ω generates μ_d. Test [s0·ω^i]P0 = T
   for each i and output s.

**Machine model** (the candidate's own, stated so that every charge below is a scalar AArch64 instruction
count). A node is an AArch64 core with up to 2^63 bytes of local memory. All memory is in one global
address space: a store or load to another node's memory is one instruction per 64-bit word, like a local
one. A node with nothing to do executes WFI (one instruction) and sleeps, executing nothing; it is woken
only by a software-generated interrupt (SGI) targeted at it by name, never by a broadcast event, and a wake
plus handler dispatch is charged 2^8 instructions. There are 2^60 table nodes (the hash table, 2^122.24
bytes in all, 2^62.24 bytes per node), 2^20 plan nodes (one class plan each) and 16 walk nodes per walk
(2^12 lanes each; a lane holds one prime's table of 2^44 affine points, 2^50.6 bytes, and 2^15 + 1 point
slots: 2^62.6 bytes per walk node, inside the 2^63).

**Batching.** Each lane's table build and walk is a sequence of dependent operations; the 2^16 lanes of a
walk are advanced in lockstep rounds. Every round performs 2^16 independent operations, one per lane,
sharing one Montgomery batch inversion across the 16 walk nodes: each node forms the prefix products of its
2^12 denominators (the 3M per lane of the reference convention), sends its product (6 words) to the group
leader, which inverts the 16 products with one field inversion and a second-level Montgomery trick (3M per
node), and sends each node its inverse back. The one field inversion per round (an exponentiation or a
binary extended gcd in F_p, ≤ 1,000 field multiplications = 383,000 instructions) costs 383,000 / 2^16 =
5.85 instructions per lane operation; the second level (16 × 3M + 2 × 16 × 12 words of transfer + 32 wakes
at 2^8) costs under 0.3 per lane operation. Both are charged at 10 instructions per lane operation, on
every walk, root, pair and table-build operation. A lane whose plan is exhausted starts the next prime.
- Q_a = [a]P0 (or [a′]A_d) is computed in Jacobian coordinates and normalised by one inversion, charged in
  the per-prime scalar budget (1,000 units for a ≤ 77-bit scalar, section 5).
- T[2] is one affine doubling, also charged there.

**The pair operation.** For P = (x1, y1) and T = (x2, y2) with x1 ≠ x2:
- P + T: λ = (y2 − y1)·δ, x3 = λ² − x1 − x2, y3 = λ(x1 − x3) − y1, with δ = (x2 − x1)^−1;
- P − T = P + (x2, −y2): λ′ = (−y2 − y1)·δ, x3′ = λ′² − x1 − x2, y3′ = λ′(x1 − x3′) − y1.

One denominator, so one batch-inversion share (3M), then 1M + 1S + 1M for each of the two formulas: 9M in
all, with S = M. The second formula's field additions, subtractions and negations (the reference
convention sets A = 0 for one formula; the second one is extra work) are charged explicitly in section 5.
`evidence.py` executes exactly this formula in the toy attack and counts 4.555 field multiplications per
element at a planned pair fraction of 0.926.

**Incomplete additions are safe.** Let Q = Q_a ≠ O (a, t ≠ 0 mod r) and write every point as a multiple of Q;
the x-coordinates of [u]Q and [v]Q coincide only if u ≡ ±v (mod r).
- *Pair op.* P = [x]Q, T = [g]Q with 1 ≤ g ≤ GMAX and x − g ≥ 1 (it is a list element), so 0 < g < x and
  x + g < 2^77 < r: neither x ≡ g nor x ≡ −g.
- *Root op.* root′ = [v]Q = [u]Q + [gap]Q with u < v in L_c and gap = v − u ≤ GMAX. u ≡ −gap would force
  v ≡ 0. u = gap means v = 2u, so χ(v) = χ(2) + χ(u) ≠ χ(u) because χ(2) ≠ 0: then v ∉ L_c, contradiction.
- *Table build.* T[k] = T[k−1] + Q collides only if k − 1 ≡ ±1, that is k = 2, which is the doubling.
- The giant points are multiples of Q = [a′t]P0 ≠ O, so the same argument covers them.

## 3. The planner

Fix a class list L = L_c = (b_0 < b_1 < … < b_{ℓ−1}), ℓ ≈ X2/d = 2^56 (H2). Call (i, j) with i < j a
**pair with midpoint k** if b_k − b_i = b_j − b_k = g with 1 ≤ g ≤ GMAX. In a plan, every element is either
a **root** or a **child** of exactly one chosen pair; the midpoint of a chosen pair is any element (root or
child), and may be the midpoint of several chosen pairs. Dependencies go from the midpoint to its two
children. The plan is valid if the dependency graph is acyclic, if every dependency tree has depth at most
D = 2^14, if consecutive roots are at most GMAX apart (so that a root is one table addition away from the
previous root), and if b_0 is a root.

**Construction.**
1. *Edges.* For each u ∈ L and each v ∈ L with u < v ≤ u + 2·GMAX and v ≡ u (mod 2), test whether
   (u + v)/2 ∈ L by reading its entry in the χ array. The window of 2·GMAX = 2^45 integers holds about 2^25
   elements of L, half of them of the right parity: about 2^24 tests per element, each a load of v, a
   parity test, the midpoint, an unaligned 20-bit load, a mask, a compare and a branch (16 instructions).
   Under H2 the density of L is 2^−20, so an element lies in about 2·GMAX/d² = 32 pairs as a child and 16
   as a midpoint (the toy census measures 30.5 to 30.7 as a child).
2. *Matching.* Choose pairs so that no element is a child twice: greedily in order of increasing degree,
   then one pass of length-3 augmenting paths (an unmatched element takes a partner away from a matched pair
   whose other child can be re-paired with a third unmatched element; if that fails the original pair is
   restored without re-checking, which is safe because every later acceptance was checked against heights
   that already included the restored subtree, heights being never lowered). Before a pair (i, j, k) is accepted: if depth(k) + 1 + max(height(i), height(j)) > D the pair
   is rejected, where height(v) is an upper bound on the height of the subtree already hanging below v
   (maintained on every acceptance, never lowered, so the check can only over-reject); then the parent
   chain of k is walked, at most D steps, and if it contains i or j the pair would close a cycle and is
   rejected.
3. *Root chain.* Unmatch b_0 if it is a child. Sweep the list: whenever an element is more than GMAX beyond
   the last root, unmatch its predecessor, which becomes a root (its pair partner becomes a root too). Then
   consecutive roots are at most GMAX apart, since consecutive list elements are (gaps ≈ d ≪ GMAX; checked).
4. *Emission.* Sweep the list in increasing order. A child is skipped (it is emitted by its root's cascade).
   A root is emitted as a root entry (b_0 as a scalar multiplication), followed by a depth-first traversal of
   its subtree: at each element, one pair entry per chosen pair with that midpoint, then recurse into the two
   children. Every element is emitted exactly once; a leftover element would signal a cycle, which the
   check in step 2 excludes. In the traversal the midpoint of a pair entry is always on the current path, so
   its point has been computed, and the stack never holds more than 2D + 1 = 2^15 + 1 points; the script
   asserts both bounds. A plan entry is 32 bytes: kind (root or pair), g or gap (44 bits), the source slot
   and the two destination slots (15 bits each; the engine needs no allocation), and the list indices of the
   one or two elements produced (57 bits each, for the locators). The plan of class c is challenge
   independent and is shared by every prime of character −c on both sides: it is held by one plan node and
   streamed, entry by entry, into a ring at the walk lane that uses it.

**Output.** The pair fraction f_c = (number of children)/ℓ of class c is known before any group operation.
The number of roots is ℓ(1 − f_c), the number of entries ℓ(1 − f_c/2). Cycle rejections concern a 2-cycle
(a four-term progression x−g, x, x+g, x+2g with x and x+g choosing each other) or longer returns of the
parent chain; their rate is proportional to 1/d (toy: 54 and 48 rejections among 4,046 pairs at d = 16),
and a rejected pair only makes the element try its next edge. Demotions in step 3 happen when 2^24
consecutive elements have no root, never in the toy beyond the boundary effect of its short lists
(GMAX = span/32). The depth cap never bound in the toy (maximum depths 130 and 137 against D = 2^14; the
dependency chain length is geometric with mean about 1/(1 − f), so a chain of length 2^14 has probability
about f^{2^14} < e^{−90} per element even at f = 0.9945).

**H4 (pair fraction).** At the production parameters every class list has f_c ≥ 0.98. Support: on real
χ-class lists of a toy prime with d = 16 and GMAX/d² = 16 (the production value 2^44/2^40), all 16 class
lists of length 2^13 to 2^14 give f ∈ [0.9931, 0.9944] (the two in the declared run: 0.99386 and 0.99361).
The deficit has three sources, all smaller at production size: boundary effects (list span 2^17 against
2·GMAX = 2^13 in the toy; 2^76 against 2^45 in production), cycle rejections (∝ 1/d), and demotions
(∝ GMAX/span). The claim leaves 0.013 below the smallest toy value. The depth cap is enforced
deterministically by the planner, so memory does not depend on H4; H4 affects only the cost.

## 4. Success probability

Steps 5 and 6 are deterministic. The algorithm succeeds exactly when t ∈ S := U·B·G^{−1} and the candidate
cap is not exceeded. Because t is uniform on H, Pr[t ∈ S] = |S|/m.

- **H1 (coverage).** Apart from the structural loss below, the 6N² products u·b·g^{−1} coincide no more
  often than random products would. Then |S|/m ≈ 1 − e^{−K} ≥ 0.399504 for K = 6N²/m ≥ 0.51. Two kinds of
  coincidence are ruled out: the six candidates of one pair of locators are distinct because U acts freely
  on F_r^*; and two different pairs of locators give the same candidate with u = u′ only through an
  integer identity b·g′ = b′·g (the structural loss below), with u/u′ = −1 never (b·g′ + b′·g < 2^153 < r).
  With u/u′ = ±λ^j the coincidence is a congruence b·g′ ≡ ±λ^j·b′·g (mod r), not an integer identity; such
  coincidences are ordinary random collisions, already counted by the model 1 − e^{−K}, and have no
  systematic source: a systematic one would need ±λ^j ≡ f1/f2 for integers f1, f2 below 2^153 shared by
  many pairs, and there are at most 2^306/r = 2^51 such residues among the r possible ones.

  *Structural loss.* Rationals repeat when b/b′ = c/c′, through a shared factor g = gcd(b, b′) and a
  cofactor k with χ(k) = χ(g). The expected number of repeats per pair is at most
  2·E[gcd]/d ≈ 2·(6/π²)(ln X2 + 0.15)/d < 2^−14 for X2 = 2^76. The planner does not change the sets: it
  only changes the order and the formula by which the same elements are computed.
- **H2 (equidistribution).** The values χ(b), b ≤ X2, and χ(a), a ∈ I, are spread evenly over Z/d, so
  every class list has about X2/d = 2^56 elements and the walks need about 2^60.65 primes in all, charged
  at 2^61.65. The interval I holds more than 0.9·2^70/ln(2^77) = 2^64.07 primes (the prime number theorem
  in intervals of length ≥ x^{7/12+ε}, proved by Huxley; here 2^70 ≫ 2^45; the factor 0.9 is a margin),
  at least 5× slack. The exact class lengths and the exact number of primes are computed by the sieves,
  so H2 affects only the ledger's small terms, never success.
- **H3 (keys).** The 256-bit keys of the 2N records behave like uniform strings for sharding, probing and
  false matches. The expected number of false matches is N²/2^256 = 2^−24.7; the run fails if more than
  2^20 candidates arise, an event of probability below 2^−44 by Markov, budgeted at 2^−20.

Claimed success: 0.399504 − 2^−14 − 2^−20 > 0.3994, claimed **0.399** (margin 0.009 over 0.39).

Evidence (`evidence.py`, toy scale only).
- *Coverage census.* Exact |S|/m for a toy prime r ≈ 2^26 with d = 64, K = 0.5103 and prediction 0.3997:
  random sets 0.4005; smallest integers in H 0.4014; product construction at X2/d = 64 and 256: 0.3947 and
  0.3985. The shortfalls (0.0050 and 0.0012) shrink as X2/d grows, as in cbadde77 (whose census at
  X2/d = 256 and r ≈ 2^29 gave the prediction exactly); production has X2/d = 2^56. The gcd bound above is
  a conservative upper bound on the loss, not a model of these deficits.
- *Toy attack.* The whole algorithm on a 24-bit j = 0 curve with d = 32 and 2 ∉ H, X2 = 2^13, GMAX = 2^12,
  211 elements per side (K = 0.5113, prediction 0.4003), plan pair fraction 0.926: 700 trials gave success
  0.3886 (standard error 0.0185, 0.6σ low, at X2/d = 256 where a small toy-scale loss is expected); every
  recovered scalar was correct; 4.555 field multiplications per element against 6 for the reference walk.
  Other seeds of the same run (2 and 3, not declared) gave 0.395 and 0.407 in 1,000 trials each at
  K = 0.5065 (prediction 0.397).

## 5. Cost (unit: one batched incomplete affine G1 addition = 6M = 2,000 instructions)

Instruction budgets on top of the field multiplications, which are priced by the reference convention:
a root op is 6M = 2,000 instructions, a pair op is 9M = 3,000 instructions for two points.

**Pair op, two points** (listed 325, charged 340 = 170 per point):

| item | instructions |
|---|---:|
| plan entry: four loads from the plan ring, decode kind, g, slots and list indices | 8 |
| load the midpoint P from its slot: 12 loads, 2 address | 14 |
| load T[g]: 12 loads, 3 address | 15 |
| batch-inversion lane buffers (denominator and prefix product, beyond the 3M) | 24 |
| second formula's field operations: −y2 − y1 (2), λ′² − x1 − x2 (2), x1 − x3′ (1), λ′(x1−x3′) − y1 (1): 6 modular additions or subtractions of 384-bit operands in memory at 36 each (12 loads, 6 adds or subs with carry, 6 conditional corrections, 6 selects, 6 stores) | 216 |
| store two points: 24 stores, 4 address | 28 |
| share of the per-round field inversion and group exchange | 10 |
| lane control and loop | 10 |

**Root op, one point** (charged 160): the walk bookkeeping of cbadde77, 150 (listed 72: read the gap (4),
address and load T[gap] (15), load and store the point (24), batch-inversion buffers (24), loop (5)), plus
the inversion and group share 10.

**Key** (charged 40 per point, listed 25): p − y (7), 6-limb compare (12), select (6).

**Hash insert (baby) or lookup (giant)** (charged 180 per point, listed 173). Records are 48 bytes: 4 words
of key and 2 words of locator (prime index ≤ 2^64.3, list index ≤ 2^57: 122 bits). The table has 2N slots
(load ≤ 1/2), 2^122.24 bytes, sharded over 2^60 table nodes by the high key bits (2^62.24 bytes per node).
Node ids are five 12-bit digits; a record is forwarded along the digits of its destination: a walk node
sends it to the node (d1, 0, 0, 0, 0), which sends it to (d1, d2, 0, 0, 0), and so on: 5 hops, 6 charged.
Every node keeps 2^12 outbound buffers of 1,024 records (48 KB) and, per sender, an inbound ring of four
such blocks. Per record:

| item | instructions |
|---|---:|
| shard and slot from key bits | 4 |
| append the 6-word record to the outbound buffer of the next hop, advance its pointer | 8 |
| 6 hops: the relay loads the 6 words from its inbound ring (6), stores them into the outbound buffer for the next hop (6), the block transfer into the next node's inbound ring (6 stores), control (2) | 120 |
| block protocol (below), 530 instructions per block of 1,024 records per hop, 6 hops | 4 |
| read at the receiving node | 6 |
| linear probing at load ≤ 1/2: ≤ 2.5 expected probes of 6 | 15 |
| store the record (insert), or compare the locators and emit a candidate (lookup) | 6 |
| loop and control | 10 |

*Block protocol.* A sender holds a credit of four blocks per receiver ring. To send a full block it reads
its credit word (1 load), writes the block into the ring (the 6 stores per record above), pushes the block's
descriptor onto the receiver's descriptor ring with an atomic fetch-and-add on the tail (1 LSE instruction)
and a store (1), and raises the receiver's SGI (1 store). The descriptor ring has as many entries as the
sum of all credits, so a push never overflows. The receiver's handler (wake 2^8) pops descriptors in order
(2 each), processes the block, and returns the credit by storing the sender's credit word (1) and raising
the sender's SGI (1), in case the sender was sleeping for credit (wake 2^8). Per block: at most 2 × 2^8 +
8 < 530 instructions, 0.52 per record per hop. A relay with an empty descriptor ring sleeps; it executes
nothing until the next descriptor's SGI. Every wake is caused by one accounted push or credit return, so
the number of wakes is at most twice the number of blocks; a node woken while sleeping for credit on one
buffer re-checks that credit word and either continues or sleeps again, inside the 2^8 charged to the event
that woke it. A node sleeping for credit does not service its inbound rings meanwhile; blocks accumulate
upstream within their credits, and since hop digits strictly advance the forwarding graph is acyclic, so
there is no deadlock and no storage or instruction cost beyond the per-block charges. At the end of a walk
each node flushes its partial buffers (2^12 blocks per node, inside the per-record charge). Rings are never
written beyond their credits, so no record is lost or overwritten. Nodes do not poll.

These are expected costs summed over 2^117 operations, where linear probing and balanced sharding
concentrate sharply (H3); the baby side is challenge independent, so its balance is a property of the fixed
set, modelled by H3. A match reply is needed only for candidates (≤ 2^20). The χ array and class lists are
read sequentially by the sieves and by the planner's window scan (inside the window-test budget). Each
walk lane streams the plan of its class from the plan node: per entry the plan node loads 4 words and
stores them into the lane's plan ring (9 instructions with control, the same block protocol at the same
per-record price, absorbed in the 9), and the lane reads them (the 8 above).

Per-point costs: c_pair = (3000 + 340 + 2·(40 + 180)) / (2·2000) = **0.945**;
c_root = (2000 + 160 + 40 + 180) / 2000 = **1.19**. With pair fraction f ≥ 0.98 (H4), the cost per element
is c = f·c_pair + (1 − f)·c_root ≤ **0.9499**.

The full ledger, recomputed exactly by `evidence.py` (n_a = 2^61.65 primes charged on both sides together):

| Term | log₂ (units) |
|---|---:|
| baby walk N·c | 115.576 |
| giant walk N·c | 115.576 |
| plan streaming: 9 instructions per entry at the plan node, ≤ N(1 − f/2) entries per side | 107.883 |
| per-prime tables T[k], k ≤ 2^44: n_a·GMAX additions at 2000 + 60 + 10 instructions | 105.700 |
| planning: per b ≤ X2, 2^24 window tests × 16 + matching 2^19 (32 edges × a chain walk ≤ 2^14) + emission 64 | 93.037 |
| χ for the primes ≤ X2 (100 units = 2·10⁵ instructions each: exponentiation plus Pohlig–Hellman) | 77.258 |
| per-prime scalar multiplications (Q_a, first element, T[2]: 3 × 1,000 units) | 73.201 |
| χ-sieve and class lists to X2 (200 instructions per b) | 72.678 |
| χ for the primes of I (100 units each, all 2^64.3 of them) | 70.943 |
| prime sieves: the interval I and the primes ≤ 2^39, 64 instructions per integer | 65.034 |
| candidates: up to 2^20, six 255-bit scalar multiplications (3,000 units) each | 34.136 |
| μ_d search (2^20 + 1 scalar multiplications), d-th root, input checks | ≤ 31.6 |
| **total W(2^20)** | **116.5804426** |

- Scalar multiplications are charged at 1,000 units for a ≤ 77-bit scalar and 3,000 units for a 255-bit
  scalar: generous Jacobian ladders at ≤ 5 units per step plus one inversion at 100 units.
- W = 124225322321011948872373830261525335 units (integer ceiling). The score is rounded up:
  log₂ W ≤ 116.580443 = 116580443/1000000 (`evidence.py` evaluates log₂ W with 80-digit decimal
  arithmetic and rounds up at 10^−6; the preprocessing and memory bounds are rounded up at 10^−4, the
  success bound down).
- For context, the same budgets with every point computed as a root (no pair ops) give 2^116.9047; the
  pair trick is worth 0.324 bits, the rest of the gain over cbadde77 (116.982765) comes from K (0.52 →
  0.51, 0.014 bits), 48-byte records and the 6-hop network, less the costs added here (the second
  formula's additions, the inversion share, the block protocol, plan streaming).

**Preprocessing** (challenge independent: the baby walk and its tables, the class lists, characters,
sieves, planning and plan streaming) ≤ 2^115.5847 units, included in W.

**Memory** (exact sum in `evidence.py`, ≤ 2^122.2354 bytes).
- Hash table: 2N slots × 48 bytes = 2^122.235 bytes.
- Per node (2^60 + 2^20 + 16 nodes): 2^12 outbound buffers of 48 KB, 2^12 inbound rings of 4 × 48 KB,
  code and processor state 4 KB: 2^29.9 bytes per node, 2^89.9 bytes in all.
- Class lists and plans: 27 bytes per b ≤ X2 (10 for the list entry, 32 per plan entry at ≤ 0.51 entries
  per element, 26.3), 2^80.8 bytes; χ array: 2.5 bytes per b, 2^77.3 bytes.
- Tables and slots for the 2^16 lanes: 2^16 × (2^44 + 2^15) × 96 bytes = 2^66.6 bytes.
- Primes of I with their characters: 16 bytes each, 2^68.3 bytes.
- Memory is not scored. Any BSGS form of Cheon stores one side, about 2^115.7 entries.

## 6. W(n) and applicability

The algorithm uses only A_1 and A_{2^20}, so W(n) = 2^116.580443 for every 2^20 ≤ n < r. Nothing is
claimed for n < 2^20. More generally, for any d | r−1 with 6 | (r−1)/d and d ≤ n:
W ≈ 2·c·√(K(r−1)/(6d)) plus the lower-order terms above, with X2 = d·2^56 and GMAX = 16·d².

## 7. Context (unscored) and disclosed heuristics

| Method | log₂ W |
|---|---:|
| Pollard rho on G1 with order-6 classes, plus overhead | ≥ 126.5 |
| Cheon BSGS with KKM c = 3 tables (2 additions plus one F_r multiplication per step) | ≈ 117.9 |
| cbadde77 (`additive-cheon-v1`, accepted): one addition per element | 116.983 |
| this accounting without pair ops | 116.905 |
| **this work** | **116.580** |

- The gain is a constant factor. The generic floor of this family (2N additions at c = 1, K = 0.494 for
  success 0.39) is 2^116.63 in additions; the pair walk goes below it because a pair-derived point costs
  0.75 additions in field multiplications. The floor for pair-derived points is 2^116.24 at f = 1 and zero
  overhead; this entry sits 0.35 bits above it, almost all of it the hash, key and bookkeeping budgets.
- Two-for-one P ± T had been set aside in other entries because a random list contains both x+g and x−g
  only with probability d^−2 per offset. The point here is that the offsets are planned offline from the
  integer lists, where every three-term progression with g ≤ 2^44 is found by a window scan of the χ array
  at a cost far below the walks, and the dependencies form a forest that a depth-first traversal executes
  with one shared inversion per pair.

**Disclosed heuristics:** H1 (coverage, with its structural loss bounded), H2 (equidistribution, affecting
only small terms), H3 (keys as uniform strings, with the candidate cap counted as failure), H4 (pair
fraction ≥ 0.98 at production parameters, measured 0.993–0.994 at toy scale). The machine model of
section 2 is the candidate's own. No other assumptions. This is an analytic bound with itemised instruction
counts; nothing runs at full size, and the toy results are evidence at their own scale only.

# Fixed-step walks with a routed class table on BLS12-381 G1

**Claim.** time_log2 = 126.17719. For every challenge T = [x]P0 (hence for a
uniform one) the algorithm below outputs x with probability at least 0.39 over
its own coins, and its total work is at most W_max = 2^126.17719 units (one unit =
one batched incomplete affine G1 addition = C = 2,000 scalar AArch64
instructions). The computed total is W < 2^126.177181; the claim keeps a margin
of 0.00001 in log2. The work bound holds on every run: every data-dependent
action is checked against its allotment immediately before it runs and the walk
set-up samplers have explicit attempt caps (section 2, last two paragraphs). Peak memory < 2^132.3 bytes (section 6); no advice;
preprocessing (walk set-up, routing tables, table zeroing) < 2^116.2 units, all
inside W.

**What is different from the six-class rho line of submissions.** Those walks
must be functions of the class, so every step computes the canonical twist index
(about 50 of their 98 per-step instructions) and runs cycle handling; their
success law is the random-mapping heuristic H1 with a 0.5% allowance. Here the
walks are plain arithmetic progressions with a per-walk fixed step, so no twist
index, no cycles and no distinguished points are needed; instead every point's
class key is routed to a distributed **direct-bucket chained** table that
detects any two points in the same class. The success probability is **proved**
(section 4) from exact linear-algebra probabilities of the collision events over
the algorithm's coins and the Bonferroni inequality; no random-mapping heuristic
is used. The one heuristic left (H', key bits behave as uniform) is used only for
the cost-side bounds (the expected chain scan and its exact variance, window
counts, false candidates) and never enters the probability that a collision
exists.
Precisely: the collision law is proved; the abort and overflow bounds use H'. The price is memory: the table holds all ~2^126.13 points
(memory is reported but is not part of the score under the track rules).

## 1. Algebra (checked in `experiment.py`, part A)

p = field prime, r = subgroup order (log2 r = 254.85709), E: y^2 = x^3 + 4,
G1 = <P0>.

* beta = 0x1a0111ea...00000000aaac is a primitive cube root of 1 in Fp and
  phi(x, y) = (beta x, y) = [lambda] on G1 with
  lambda = 0xac45a4010001a40200000000ffffffff, lambda^2 + lambda + 1 = 0 mod r
  (phi(P0) = [lambda]P0 is checked).
* U = {±1, ±lambda, ±lambda^2} ⊂ Z_r^*, |U| = 6, acts by
  [(-1)^s lambda^j](x, y) = (beta^j x, (-1)^s y). The class of Q ≠ O is
  {[u]Q : u ∈ U} (size 6).
* **Class key** K(Q) = y if y is even, else p - y. p is odd, so exactly one of
  y, p - y is even; K is unchanged by negation and by phi (which keeps y), so K
  is constant on a class. Conversely K(Q) = K(Q') gives y' = ±y, then
  x'^3 = x^3, so x' ∈ {x, beta x, beta^2 x} and Q' ∈ class(Q): **K is a
  complete class invariant**. Bit 0 of K is 0. Write K = K0 + 2^64 K1 +
  2^128 K2 + ... (64-bit limbs).
* **Relations.** If [s]P0 and [s']P0 (s, s' ≠ 0) are in one class then
  s = u s' (mod r) for exactly one u ∈ U (no non-identity unit fixes a nonzero
  scalar), and u is found by comparing x-coordinates (x = beta^j x') and the
  sign of y.

## 2. Algorithm

Parameters: M = 2^63 walks, each with L + 1 points, L = 2^63.128335 (exactly:
the smallest integer with f4(mu) >= 0.390001, see section 4), so the number of
points is n = M (L + 1) = 2^126.128335.

**Walks (fresh coins).** For each walk w: c_w, b_w uniform in Z_r and D_w
uniform in Z_r \ {0}. Each uniform Z_r draw samples a 255-bit word and retries
while the word is >= r, with an explicit cap of 128 attempts; D_w draws a
uniform Z_r word by the same capped procedure and retries only if it is 0, with
128 further capped attempts. A walk whose sampler hits a cap halts the attack
with failure; one uniform draw is rejected with probability q = 1 - r/2^255 <
0.09432, so a 128-attempt cap fails with probability q^128 < 2^-436, and the
3*2^63 capped draws (c_w, b_w, D_w) give a total sampler-cap failure
probability below 3*2^63*2^-436 < 2^-371 (the exact value is 2^-371.44), and the
worst-case set-up still fits the 2^14-unit-per-walk allotment below (128
attempts cost far less than one scalar multiplication). Points
Q_{w,j} = [c_w + j D_w]P0 + [b_w]T for j = 0..L; the start Q_{w,0} is computed
by two scalar multiplications and the step point S_w = [D_w]P0 by one; then
Q_{w,j+1} = Q_{w,j} + S_w is **one batched incomplete affine addition** (the
reference operation). (c_w, b_w, D_w) are kept in a walk directory indexed by
the 63-bit walk id. Write A_w = c_w + b_w x, so Q_{w,j} = [A_w + j D_w]P0.
If T = O the attack outputs 0 at once.

**Generation.** 2^39 generator nodes each own a batch of 2^24 walk records and
advance all of them in lock step: one Montgomery-batched inversion per batch
(one batch = one value of j). After a batch has stored its new points, the
generator loop (listing G, section 5) computes limbs 0-2 of K for each new
point and appends the 32-byte entry (K1, K2, walk id, j) to the window of the
intermediate node digit1 = bits 1..36 of K0.

**Routing (two hops, remote stores).** Intermediate node d1 (2^36 of them)
forwards each entry (listing I) to final node (d1, d2) with d2 = bits 0..34 of
K1. There are 2^71 final nodes. All transfers are ordinary AArch64 stores into
memory-mapped receive windows of the next node (a PGAS / RDMA put); each
transfer is therefore charged as the store instructions that perform it, and
the receiver's loads are charged when it reads the entry. No other network
work is assumed free of charge: there is no network instruction other than
these loads and stores (section 6 discusses the machine model).

**Final table (listing F).** Each final node keeps a direct-bucket table of
S_node = 2^57 eight-byte head slots, zeroed once at start-up; a head is the
address of the node's first 40-byte entry `(next, K1, K2, walk id, j)` in its
arena, or 0. The home bucket is floor(h64 * S_node / 2^64) with
h64 = bits 35..98 of K2:K1 (`extr` then `umulh`). Insertion scans the chain
from the head, comparing (K1, K2) at every entry and raising a **candidate**
on equality, then appends the new entry at the tail (writing the old tail
pointer, or the head if the chain is empty). Nothing is deleted, so the table
holds every point. The global load is lambda = n / S_tot = 2^-1.8717 = 0.2733
(S_tot = 2^71 S_node = 2^128 slots).

**Candidates.** A candidate gives two locators (w, j), (w', j'). The
verification action recomputes both points from the directory, compares their
full keys (false candidates stop here), finds u with Q_{w,j} = [u]Q_{w',j'}
(section 1), and if b_w - u b_w' ≠ 0 sets

    x = (u (c_w' + j' D_w') - c_w - j D_w) / (b_w - u b_w')  (mod r),

checks [x]P0 = T and outputs x. Each candidate is priced at 2^14 units
(four scalar multiplications plus keys, unit search and the check: < 2^13.2).
A global 128-bit candidate counter, checked immediately before each
verification, caps the number of verifications at 2^90, i.e. 2^104 units.

**Rounds.** Generation is organised in rounds of 2^22 batches. In a round each
generator writes into double-buffered windows of capacity 2^11 entries per
intermediate (expected 2^10); each intermediate writes into windows of
capacity 2^15 per final node (expected 2^14); window write pointers are reset
and entry counts are published at round ends. A window can only overflow if a
Binomial count exceeds twice its mean (probability < 2^-400 overall under H').

**Random-cost allotment and abort, checked per action.** The only
data-dependent costs are (i) chain comparisons in the tables,
(ii) candidate verifications and (iii) zero-denominator recovery actions. The
attack is **one sequential program** (section 6): it keeps a 128-bit
inspection counter x27:x26 with the high-word global cap x28 and a 128-bit
candidate counter x30:x29 with the high-word cap x19 = 2^26 (at most 2^90
verifications). These are ordinary registers, not per-node counters;
the 2^71 final nodes are memory regions visited in sequence, so the same
counter covers all of them and no atomic operation is involved. The first
instructions of every chain comparison are
`adds x26,x26,#1; adc x27,x27,xzr; cmp x27,x28; b.hs abort` (listing F, all
four instructions are counted in the 11-instruction per-comparison price), so
the program halts *before* the comparison that would exceed the global
allotment; every candidate is likewise counted
(`adds x29,x29,#1; adc x30,x30,xzr; cmp x30,x19; b.hs abort`) and compared
with its cap *before* the 2^14-unit verification runs. A zero-denominator
batch is recovered by its 2^17-unit scan only while the
generator's recovery counter is below 4, checked before the scan; a generator
makes < 2^63.4 batches, each zero with probability < 2^24 * 3 / r < 2^-228.8,
so under H' some generator exceeds 4 recoveries with probability < 2^-120.
Hence every run does at most the deterministic work, the global inspection
allotment, the global candidate allotment, at most 2^39 x 4 x 2^17 = 2^58 units
of recovery work, and the one round (2^-41 of the total) of bookkeeping
overshoot already priced in section 5. This replaces the round-end aggregation:
no round can spend more random work than the global allotments. The global
allotments are exceeded with probability < 2^-108.6 (inspections, section 4)
and < 2^-400 (candidates, section 4), so they are failure events that do not
threaten the 0.39 success bound.

## 3. Correctness details

1. **Detection is complete.** Two points in one class have equal K, hence the
   same digit1, digit2, home bucket and (K1, K2). Chains never delete, so when
   the later of the two is inserted the earlier one is already on the chain of
   the shared bucket, and the insertion scan compares (K1, K2) at every entry
   of that chain: the pair is a candidate. There is no empty-marker corner:
   the head word 0 means "no entry", while a stored entry's (K1, K2) are
   compared in full.
2. **Every true candidate solves.** With Q_{w,j} = [u]Q_{w',j'} we have
   A_w + j D_w = u (A_w' + j' D_w'), i.e.
   (b_w - u b_w') x = u (c_w' + j' D_w') - c_w - j D_w, so the formula gives x
   unless b_w = u b_w'. Since b_w is uniform and independent of the event below
   (A_w is uniform given b_w because c_w is), that has probability 1/r per
   event: < 2^-120 in total.
3. **Exceptional additions.** Q_{w,j} = O or Q_{w,j} = ±S_w needs
   A_w + j D_w ∈ {0, ±D_w}: probability 3/r per point, < 2^-126 in total. A
   batch with a zero denominator is detected once per batch (the batch product
   is 0), recovered by a separately priced scan, and the walk is abandoned
   (counted as failure; the scan is priced 2^17 units and has probability
   < 2^-125 in total). A probe run longer than the 2^12-slot guard tail has
   probability < 0.78^4096 per table under H' (term "guard_tail_overflow").

## 4. Success probability (proved over the algorithm's coins)

Fix any x. Consider only **cross-walk events**
f = (w, j, w', j', u), w < w', j, j' ∈ [0, L], u ∈ U:

    E_f :  A_w + j D_w = u (A_w' + j' D_w')   (mod r).

There are N_ev = 6 C(M, 2) (L + 1)^2 of them and mu := N_ev / r = 0.494672.
By section 3, success ⊇ (∪_f E_f) minus the failure events listed there.

**Probability space.** Z = (A_1, D_1, ..., A_M, D_M). In the attack the A_w are
uniform and independent of everything else (c_w is uniform), and the D_w are
uniform on Z_r \ {0}. Let P_u be the law with all coordinates uniform on Z_r.
The attack's law is P_u conditioned on all D_w ≠ 0, so
P(∪E) ≥ P_u(∪E) - P_u(some D_w = 0) ≥ P_u(∪E) - 2^-191.

**Lemma 1 (exact intersections).** Each E_f is the vanishing of the linear form
l_f(Z) = A_w + j D_w - u A_w' - u j' D_w'. For any set F of events,
P_u(∩_{f∈F} E_f) = r^-rank(F), because a homogeneous linear system with a
uniform unknown vector holds with probability r^-rank. Distinct f give
non-proportional forms (the A_w-coefficient is 1, and then the A_w', D_w and
D_w' coefficients determine u, j and j'), so **every pair has rank 2** and
S2 = C(N_ev, 2) r^-2 exactly; S1 = N_ev / r exactly.

**Lemma 2 (dependent sets are negligible).** Let n = L + 1.
(a) *Support argument.* In a linearly dependent set of forms, every walk whose
A-variable occurs must occur in at least two forms (its coefficient must
cancel), so the "walk graph" of a minimal dependent set has minimum degree 2.
(b) *Triples of rank 2.* Either a triangle w1, w2, w3, where cancellation of
A_w3 and D_w3 forces the two forms at w3 to use the same point and the units to
compose (three points of one class): at most 36 M^3 n^3 triples; or three forms
on one walk pair. On one pair, if two units agree all three do and the three
grid points (j, j') are collinear in Z_r^2; since all coordinates are below
L < sqrt(r/2), the 2x2 determinant is an integer of absolute value < r, so they
are collinear over Z: at most 6 * 8 n^4 (1 + ln n) <= 2160 n^4 per pair, using
the count sum over primitive directions (a, b) of n^2 (n/m)^2 <= 8 n^4 (1 + ln n)
(m = max(|a|, |b|); 1 + ln n < 45). If the units are distinct, j3 - j2 = alpha (j1 - j2) with
alpha = (u3 - u2)/(u1 - u2); clearing denominators gives c0 + c1 lambda = 0
(mod r) with |c0|, |c1| <= 8L. Then c0 + c1 omega lies in the prime of Z[omega]
above r (the kernel of omega -> lambda), so r divides
c0^2 - c0 c1 + c1^2 < 3 (8L)^2 < r, forcing c0 = c1 = 0, i.e. alpha would be
rational; but three distinct sixth roots of unity are never collinear. So
j1 = j2 = j3, and the j' coordinates leave at most 120 n^3 such triples per
pair. In total
N3 <= 36 M^3 n^3 + M^2 (2160 n^4 + 120 n^3) and N3 / r^3 < 2^-375.
(c) *Four-sets.* Rank 2: every triple must be degenerate; by (b) this leaves
four collinear grid points with one unit (at most 6 * 14 n^5 per pair, by the
same direction sum) or one point pair with four units (15 n^2 per pair), so the
excess of rank-2 four-sets over r^-4 is < M^2 (84 n^5 + 15 n^2) r^-2 < 2^-61.7.
Rank 3: same-pair four-sets (all of them: M^2 C(6 n^2, 4)), triangle-plus-edge
circuits (eliminate the degree-2 walk as in (b): at most 36 M^3 n times the
per-pair degenerate-triple count), two double edges sharing a walk (one point
of w1 against two points of w2 with one unit, and likewise from w3: at most
1296 M^3 n^6), four-cycles (at most 1296 M^4 n^4), and sets containing a rank-2
triple (at most N3 N_ev): the excess is < 2^-120.

**Theorem.** By Bonferroni, P_u(∪E) >= S1 - S2 + S3 - S4, with
S1 = N_ev/r, S2 = C(N_ev,2)/r^2, S3 >= C(N_ev,3)/r^3 (dependent triples have
probability r^-2 > r^-3, so they only increase S3), S4 <= C(N_ev,4)/r^4 +
2^-61.7 + 2^-120. Since C(N, k) <= N^k/k! and
C(N,3) >= (N^3/6)(1 - 3/N),

    P_u(∪E) >= f4(mu) - 2^-60,   f4(mu) = mu - mu^2/2 + mu^3/6 - mu^4/24.

At mu = 0.494672, f4(mu) = 0.3900010 (the Poisson value 1 - e^-mu is 0.39023;
Bonferroni with four terms costs 0.0003 in mu). Subtracting the failure terms
(dependent sets 2^-60, conditioning on D ≠ 0 2^-191, degenerate b
N_ev/r^2 < 2^-253, exceptional points and zero-denominator batches 2^-126 each,
arena and window overflow 2^-400 each, global inspection-allotment overrun
2^-108.61, global candidate-allotment overrun 2^-400, sampler cap 2^-371,
recovery cap 2^-120; the last four use H' or the explicit caps):

**Success >= 0.390001 - 2^-59 > 0.39 for every x.**

Collisions inside one walk (5 of the 6 units also give x there) and the
start points' collisions are extra and only help; they are ignored in the
bound. No property of the curve beyond Lemma 1 (group structure, the unit
group of order 6 and the complete invariant K) is used: the collision law is
exact algebra, not a random-mapping assumption.

**H' (heuristic, cost side only).** The digit1, digit2 and home-bucket bits of
K of distinct classes behave as independent uniform values. It is used for: the
number X of same-bucket comparisons and its multiplicity-aware variance
(below), the window counts (Chernoff), the false-candidate count, and the
per-node arena occupancy (Chernoff). It never enters the probability that a
collision exists.

**Global inspection tail (multiplicity-aware, elementary).** H' assigns
independent uniform bucket values to distinct *classes*; the table stores
points, and a true class collision forces its points into one bucket, so the
point occupancies are not multinomial. Condition on the class partition
(multiplicities c_1,...,c_m, which depend only on the walk coins, not on the
buckets) and write X = X_same + X_distinct, where X_same = sum_a C(c_a,2)
counts forced same-class point pairs and X_distinct = sum_{a<b} c_a c_b
1[B_a = B_b] counts comparisons between distinct classes. First, E[X_same] =
sum over point pairs P(same class) = 6 C(n,2)/r = 0.4947, and a union bound
over t-subsets of point pairs gives P(X_same >= t) <= (e mu/t)^t, so
X_same < 2^30 except with probability below (e*0.4947/2^30)^(2^30) < 2^-10^9.
Second, given the
partition the B_a are iid uniform on S = 2^128 buckets, so
E[X_distinct] = (n^2 - sum c_a^2)/(2S) <= n(n-1)/(2S) = 2^123.2567; for three
distinct classes P(B_a = B_b = B_c) = 1/S^2 = P(B_a=B_b) P(B_a=B_c), so the
three-class covariance is exactly zero and disjoint pairs are independent;
hence Var(X_distinct) = p(1-p) sum_{a<b} (c_a c_b)^2 <= c_max^2 E[X_distinct].
A class with 16 points occurs with probability <= C(n,16)(6/r)^15 < 2^-1810,
so c_max <= 16 and Var <= 256 * 2^123.2567 = 2^131.2567 except with that
negligible probability. The cap is 1.10 E[X_distinct] = 2^123.394, stored in
the high word as floor(cap/2^64) = 2^59.394, so the effective cap is at least
1.09 E[X_distinct]; after the < 2^30 forced same-class comparisons the slack
is at least 0.099 E[X_distinct] = 2^119.93 and Chebyshev gives

    P(X > cap) <= 2^131.2567 / (0.099 E[X])^2 + (e*0.4947/2^30)^(2^30)
                = 2^-108.61.

The earlier iid-point multinomial moments are not used; the mean bound
n(n-1)/(2S) and the 1.10 allotment are unchanged.

**Global candidate tail (class-pair union bound).** A false candidate needs
two points of *different* classes to agree in the 128 compared key bits
(K1,K2) and in digit1 (K0 bits 1..36): equality of K1,K2 already fixes digit2
and the home bucket, so the event has probability 2^-164 per pair of distinct
classes under H' (equivalently 2^-128 for the compared limbs times 2^-36 for
digit1). No comparison index and no independence between comparisons is used.
Let W_c be the number of class pairs whose (K1,K2,digit1) agree; for any k,

    P(W_c >= k) <= C(C(n,2), k) 2^-164k <= (e n^2/2 * 2^-164 / k)^k.

A class of multiplicity c contributes c points, so pairs involving a point of
a non-singleton class number at most (2 X_same)^2 <= 2^62 when X_same < 2^30
(failure (e*0.4947/2^30)^(2^30) as above); the remaining pairs are pairs of singleton
classes, which are class pairs counted by W_c. Hence the number Y of false
candidates satisfies Y < 2^62 + W_c, and at k = 2^90 - 2^62 > 2^89.9 the bound
gives

    P(Y >= 2^90) <= (e*0.4947/2^30)^(2^30) + (e 2^250.26 2^-164 / 2^89.9)^(2^89.9)
                  < 2^-2^89 << 2^-400,

so the global candidate cap aborts with negligible probability; the one true
candidate is far below the cap. (Three equalities on one bucket are dependent
— all three 36-bit residuals must coincide — but the union bound sums over
class pairs unconditionally and does not use that independence.)

## 5. Cost

All instruction counts are scalar AArch64 instructions; the reference
operation (one batched incomplete affine addition, 6M) is 1 unit and is never
re-counted. Listings are counted by `experiment.py` from the text below; the
experiment counts each body length and charges the loop control separately.

**128x unrolling of the three table passes.** The generator, intermediate and
final-node passes each process a long fixed-shape stream, so each is emitted as
128 labelled copies of its body; the 2-instruction loop control
(`subs x17,x17,#1; b.ne group`) runs once per 128 points. The charge is
therefore body + 2/128 = body + 0.015625 per point. This is the convention
already accepted on the P-256 track (128 labelled copies, end test once per 128
lanes); the code size is 128 x (16 + 8 + 11) instructions (~17.9 KB), which is
memory, not work. Rare exits (candidate, abort) use a 2-instruction trampoline
(`adr x18,cont; b target`); its cost is inside the separately priced action,
while the one-instruction conditional test stays in the body. No other load,
store, address computation or branch is omitted: every instruction of the three
bodies is listed below and counted.

**Per point, generator (16 + 2/128).** Operand traffic of the reference
addition, charged as in the accepted rho ledgers: 6 `ldp` current point,
6 `ldp` step point, 6 `stp` new point (18). Pass bookkeeping of the
Montgomery-batched addition, charged explicitly although the rho ledgers treat
it as part of the reference operation: the forward pass stores each prefix
product (3 `stp`) and the backward pass reloads it (3 `ldp`), each pass has
`cmp`/`b.ne` loop control (4; addresses advance by pre/post-indexed
addressing), plus 2 spare (12). The field multiplications and the A = 0
additions of the two passes are the reference operation itself. Then listing G
(16):

```
gen:  ldp  x0, x1, [x19, #REC]!      // y limbs 0,1 of the new point; record pointer advance
      ldr  x2, [x19, #16]            // y limb 2
      subs x3, x10, x0               // limbs 0-2 of p - y (x10-x12 = limbs of p)
      sbcs x4, x11, x1
      sbcs x5, x12, x2
      tst  x0, #1                    // ne iff y odd
      csel x0, x3, x0, ne            // K0
      csel x1, x4, x1, ne            // K1
      csel x2, x5, x2, ne            // K2
      ubfx x6, x0, #1, #36           // digit1
      ldr  x7, [x21, x6, lsl #3]     // window write pointer for digit1
      stp  x1, x2, [x7]              // remote store: K1, K2
      stp  x22, x23, [x7, #16]       // remote store: walk id, j
      add  x7, x7, #32
      str  x7, [x21, x6, lsl #3]
      add  x22, x22, #1              // next walk id
```
(loop control: `subs x17,x17,#1; b.ne gen`, once per 128 points).

A borrow only propagates upward, so the three `subs/sbcs` give limbs 0-2 of
p - y exactly; part A emulates this row limb by limb (random, boundary and
borrow-chain y) against K mod 2^192 and digit1. x22 (walk id) and x23 (j) are
per-batch registers.

**Per entry, intermediate (8 + 2/128).**

```
int:  ldp  x0, x1, [x9], #32         // K1, K2 (post-index)
      ldp  x2, x3, [x9, #-16]        // walk id, j
      and  x4, x0, #0x7ffffffff      // digit2 = bits 0..34 of K1 (valid bitmask immediate)
      ldr  x5, [x21, x4, lsl #3]
      stp  x0, x1, [x5]              // remote store into the final node's window
      stp  x2, x3, [x5, #16]
      add  x5, x5, #32
      str  x5, [x21, x4, lsl #3]
```
(loop control once per 128 entries).

**Per entry, final node: main path (11 + 2/128).** The node's 8-byte head array
is based at x22, the 40-byte entry arena at x24 (advancing), the node's bucket
count S_node = 2^57 is in x23. The complete insertion is the main path below
plus the chain scan of the next listing; the labels show the control flow.

```
fin:  ldp  x0, x1, [x9], #32         // K1, K2 of the new entry (post-index)
      ldp  x2, x3, [x9, #-16]        // walk id, j
      extr x4, x1, x0, #35           // h64 = bits 35..98 of K2:K1
      umulh x4, x4, x23              // home bucket = floor(h64 * S_node / 2^64)
      add  x5, x22, x4, lsl #3       // address of the 8-byte head
      ldr  x6, [x5]                  // first entry, or 0
      cbz  x6, ins                   // empty bucket: append at the head cell
scan: ...                            // listing F_occupied_per_slot (11 per entry)
ins:  str  x24, [x5]                 // link: head cell or old tail's next = new entry
      str  xzr, [x24], #8            // new entry: next = 0, advance 8
      stp  x0, x1, [x24], #16        // K1, K2, advance 16
      stp  x2, x3, [x24], #16        // walk id, j; entry now 40 bytes behind
```
(loop control once per 128 entries; the 11 counted instructions are the seven
before the scan plus the four after it).

**Per chain entry compared (11).** The scan compares (K1, K2), advances the
128-bit inspection counter and checks it *before* the comparison, and steps to
the next entry; x5 is moved to the pointer cell so that a scan to the tail ends
with x5 = address of the old tail's next word.

```
scan: ldp  x7, x8, [x6, #8]          // K1, K2 of the stored entry
      adds x26, x26, #1              // 128-bit inspection counter, low word
      adc  x27, x27, xzr             // high word
      cmp  x27, x28                  // immediate global cap, high word
      b.hs abort                     // halt before the comparison that would exceed it
      cmp  x7, x0
      ccmp x8, x1, #0, eq            // (K1,K2) equal?
      b.eq cand                      // candidate action (separately priced)
      mov  x5, x6
      ldr  x6, [x5]                  // next entry (offset 0)
      cbnz x6, scan
```

**Per candidate (separately priced, 4 counter instructions plus 2^14 units).**
Every candidate is preceded by

```
      adds x29, x29, #1
      adc  x30, x30, xzr
      cmp  x30, x19                  // x19 = 2^26 = floor(2^90 / 2^64)
      b.hs abort
```

and the verification recomputes both points from the directory, compares their
full 192-bit keys, finds the unit, solves and checks (priced at 2^14 units,
< 2^13.2); a true candidate is found once. The cold-path trampoline into and
out of the handler is inside this price.

**Per point total (deterministic).** 18 + 12 + (16 + 2/128) + (8 + 2/128) +
(11 + 2/128) = 65.046875 instructions, plus the per-batch work (Fermat
inversion 380 squarings + 228 multiplications at 383 = 232,864, exponent scan
< 2,000, zero-denominator test on the batch product, per-batch constants:
< 235,000 per batch of 2^24 points, i.e. < 0.0141) and round bookkeeping < 0.02
(per round: a generator resets 2^36 window pointers and publishes 2^36 counts,
< 2^38.2 instructions against 2^46 points; an intermediate reads 2^39 window
counts and loops over them, resets and publishes 2^35, < 2^42.1 against 2^49
entries; a final node reads one window per round of ~2^14 entries; together
< 0.016). The total is 65.080882 instructions per point. Point price:
1 + 65.080882/2000 = 1.032540441 units. (Start points are also charged this
price, on top of their set-up, which over-counts.) The 2^113 units of padded
tail records (at most 127 per window, 13 instructions each, over the < 2^112.2
windows) are charged as `window_padding` in the extras.

**Table zeroing.** S_tot = 2^128 eight-byte heads, zeroed with `stp xzr, xzr`
(2 heads per store) in a loop of 64 stores per 3 loop instructions:
0.5234375 instructions per head, i.e. 2^116.10 units in all.

**Chain comparisons (random, allotted).** The expected number of comparisons
is E[X] = n(n-1)/(2 S_tot) = 2^123.2567 (0.1366 per point): insertion scans
each bucket's existing chain, so the total is exactly sum_b C(m_b, 2). Each
comparison is charged 11 instructions (the 11 listed above). The global
allotment is 1.10 x E[X], i.e. 2^115.888 units; it is enforced by the 128-bit
counter x27:x26 with the high-word cap x28 = floor(2^123.394/2^64) = 2^59.394
checked before every comparison. The multiplicity-aware variance bound
(section 4) gives the overrun probability < 2^-108.6, so the every-run bound
holds with that one failure probability.

**Other terms.** Walk set-up (capped sampling and three scalar multiplications
per walk, plus directory writes) < 2^14 units x 2^63 = 2^77; routing tables and
window mappings < 2^76; candidates: cap 2^90 verifications x 2^14 = 2^104 units,
enforced before each verification by the 128-bit counter (mean false candidates
2^87.26, overrun 2^-400); capped zero-denominator recovery 2^39 x 4 x 2^17 =
2^58; one round of overshoot < 2^85.2; padded window tails < 2^113; bank
remapping < 2^114.3 instructions = 2^103.3 units, charged 2^104; final output
and verification < 2^20.

**Total.**
W <= n * 1.032540441 + 2^116.10 + 2^115.89 + 2^113 + 2^104 (candidates)
     + 2^104 (bank remapping) + 2^85.2 + 2^77 + 2^76 + 2^58 + 2^20
< 2^126.177181 <= 2^126.17719 = W_max.

For orientation: the rho ledgers price an ordinary step at 1.04903 units and
need S_ops = 1.005 S0 operations; here a point costs 1.03254 units plus 0.00070
units of zeroing and probing, and n is the rigorous Bonferroni count (0.06%
above the Poisson value instead of a 0.5% allowance). The exact ledger
computed by `experiment.py` is W < 2^126.177181.

## 6. Memory and machine model

Peak storage < 2^132.3 bytes: heads S_tot x 8 B = 2^131; entry arenas
1.05 n x 40 B = 2^131.52; receive windows (2^36 intermediates
x 2 x 2^39 x 2^11 x 32 B = 2^92, final nodes 2^71 x 2^25 B = 2^96), walk
records (504-byte stride, so every pre-indexed `ldp` offset, including
`#REC` = 504, is encodable) and directory (2^63 x 600 B = 2^72.3); total
< 2^132.29 (including the 2^76 B bank directory of section 6). Every bank is
below 2^61.4 B (a final region: heads 2^60 B, arena 2^60.5 B, window 2^56 B;
a generator: destination windows 2^53 B, records 2^33.4 B, pointer table 2^39
B; an intermediate: outgoing windows 2^57 B, incoming 2^56 B) and every
load/store address is 64-bit. The per-node arena holds 1.05 x its expected
load; a node overflows only if its multinomial count exceeds 1.05 times its
mean, a Chernoff event below 2^-400 (term "arena_overflow").

**Banked 64-bit address space, remapping specified and charged.** The
physical store is 2^132 B, more than one 64-bit address space, so the program
reaches it through a 2^64-byte *virtual* window: a `map vbase, bank`
instruction (one instruction) installs a physical bank of at most 2^61.4 B at
virtual base vbase, after which every load and store uses an ordinary 64-bit
address; a bank is remapped only when the schedule returns to it and its
contents persist in physical storage. The working set fits at every step: a
generator keeps its 2^36 destination windows (2^36 x 2^17 B = 2^53 B), its walk
records and its 2^39-byte pointer table mapped (2^53.1 B in all); an
intermediate keeps its 2^39 incoming windows (2^55 B, double-buffered 2^56 B)
and its 2^36 outgoing windows (2^57 B) mapped (2^57.1 B); during the final
phase one final region (heads 2^60 B, arena 2^60.5 B, window 2^56 B, 2^60.6 B
in all) is mapped at a time. The map instructions are: 2^39 generators x 2^36
destination windows = 2^75, once each; 2^36 intermediates x (2^39 + 2^36)
windows = 2^75.6, once each; 2^71 final regions x 2^41 round visits x 4 maps =
2^114; total < 2^114.3 instructions = 2^103.3 units, charged as the
`address_mapping` term of 2^104 units in section 5. The bank directory (2^72
base/size pairs, 2^76 B) is inside the memory total. No address exceeds 64
bits and every bank is at most 2^61.4 B; the earlier version's "one sequential
program with 64-bit pointers to 2^132 B" is replaced by this explicit
remapping schedule.

**One sequential program.** The attack is executed by one sequential processor
using the virtual window above; the 2^39 generators, 2^36 intermediates and
2^71 final tables are physical banks visited in sequence, and the "remote
stores" are ordinary stores into the mapped receive window. The global
inspection and candidate counters x26/x28 and their caps x27/x29 are therefore
ordinary registers of this one program, updated by four instructions per action
(inspection counter x27:x26 with cap x28; candidate counter x30:x29 with cap
x19; listing F and the candidate check); no atomic, network or per-node
counter is involved. The work is the sum over the sequential schedule, exactly
the same instructions as the distributed description, so the score is
unchanged.

Memory has no scalar contribution under classical-group-work-v2; it is stated
here as required. Physically 2^132 bytes is far beyond any machine, as is the
2^126-operation time of every attack on this track; the score is the
instruction-equivalent count, which this design charges completely: each
entry is touched by explicit loads and stores at the generator, the
intermediate and the final node, and remote transfers are the store
instructions that issue them.

**Interconnect.** Routing is done in software and charged: no node ever
addresses more than 2^36 peers, and the forwarding decision at every level
(digit extraction, window pointer load/update, the stores) is in listings G and
I. What the fabric provides is only point-to-point delivery of a store to a
directly mapped peer window, which computes nothing (it is the same service a
memory system gives every `str` in the cost model, and the accepted rho
submissions likewise priced each probe of their 2^71-byte report table as a few
loads and stores). It is not an oracle: it neither hashes, sorts nor compares.
If a reviewer nevertheless preferred to charge an additional software copy for
every remote transfer (2 `ldp` + 2 `stp` per hop, 8 per point), W would be
< 2^126.1882, still below every rho-line submission; the claim does not use
that variant.

## 7. Experiments (`experiment.py`, stdlib)

* **Part A (production, all true):** P0 on the curve, [r]P0 = O, beta and
  lambda, phi(P0) = [lambda]P0, the unit map on all six members, key constant
  on the six members and even, key completeness, listing G's key rows emulated
  limb by limb against K mod 2^192 and digit1 (random points, both signs, and
  16 borrow-chain boundary values), the final home index in range, AArch64
  immediate encodability (`and #0x7ffffffff`, `tst #1`), and a **planted
  end-to-end run on BLS12-381**: three walks with real batched affine
  additions, a planted class collision between two walks at a random unit,
  routing digits and a chained table built from the emulated rows (8-byte head
  buckets, 40-byte `(next, K1, K2, w, j)` entries, tail append and chain scan),
  the candidate, re-computation from the directory, the unit search, the solving
  formula, and [x]P0 = T.
* **Part B (toy curve y^2 = x^3 + 4 over F_q with a prime-order subgroup and
  the same six units):** the whole attack at toy size with the chained table,
  and the exact scalar model at small r. A run of this package's code (q30,
  32 walks, L=55) gives 0.3986 +- 0.0075 over 4,280 trials against the
  cross-walk prediction 0.3907, and the scalar model gives 0.4335 +- 0.0020
  over 62,463 trials against 0.4318. The toy runs are evidence for the
  implementation, not for the bound, which is proved in section 4.
* **Part C:** the ledger above in exact fractions: L from f4, the price, the
  zeroing and probing terms, the extras, W, the dependency-bound exponents, the
  multiplicity-aware inspection moments (X_same, X_distinct, Var <= 256 E[X],
  c_max <= 16) and the class-pair candidate union bound, the 128-bit counter
  listings, the bank-remapping charge,
  the candidate cap 2^90 = 2^104 units, the 1.05 arena margin, the 128x unroll
  accounting, and checks that every listing's counted length equals its ledger
  entry.

## 8. Limitations

* Memory is enormous (2^132.3 bytes) and is not scored by the rules.
* The machine model treats a store into a mapped remote window as one store
  instruction (section 6 gives the bound with an extra copy per hop).
* H' (uniform key bits) is a heuristic; the abort, overflow and false-candidate
  bounds rely on it (each with a stated negligible probability, and a 10%
  slack on the global inspection allotment). The inspection tail is the
  multiplicity-aware finite Chebyshev bound of section 4 (conditioned on the
  class partition; no asymptotic variance law and no iid-point assumption);
  the collision law itself is proved and does not use H'.
* The per-node arena capacity is 1.05 times the expected per-node count; a
  node exceeding it is a Chernoff event below 2^-400 under H' and is charged as
  a failure term.
* Instruction counts are analytic listings under the fixed C = 2,000, not
  hardware timings. The reference-operation operand traffic is charged as in
  the earlier accepted rho ledgers (18 transfer instructions per addition).
  The single sequential program's global counters are ordinary registers;
  their four-instruction checks are priced inside the 11-instruction comparison
  price (and four inside the separately priced candidate action).

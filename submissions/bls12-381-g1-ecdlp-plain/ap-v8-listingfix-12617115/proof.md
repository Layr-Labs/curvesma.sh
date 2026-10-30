# Fixed-step walks with a routed class table on BLS12-381 G1

**Claim.** time_log2 = 126.17115. For every challenge T = [x]P0 (hence for a
uniform one) the algorithm below outputs x with probability at least 0.39 over
its own coins, and its total work is at most W_max = 2^126.17115 units (one unit =
one batched incomplete affine G1 addition = C = 2,000 scalar AArch64
instructions). The computed total is W < 2^126.171145; the claim keeps a margin
of 0.000005 in log2. The work bound holds on every run: every data-dependent
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
the cost-side bounds (the expected chain scan and its exact variance, false
candidates, arena occupancy) and never enters the probability that a collision
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

**Direct insertion into the final table (one pass, one hop).** 2^39 generator
nodes each own a batch of 2^24 walk records and advance them in lock step: one
Montgomery-batched inversion per batch (one batch = one value of j). After a
batch has stored its new points, the generator loop (listing G, section 5)
computes K0, K1, K2 for each new point, forms the 71-bit node id
digit1 || digit2 (digit1 = bits 1..36 of K0, digit2 = bits 0..34 of K1),
selects that final node's bank (section 6) and inserts the 40-byte entry
`(next, K1, K2, walk id, j)` into the node's direct-bucket chained table in the
same pass. There is no intermediate node, no receive window and no second
pass: the generator's own load and store instructions write and read the
destination node's bank, and section 5 charges them per point.

**Final table (listing F).** Each of the 2^71 final nodes keeps its table in
its own bank: an 8-byte arena bump pointer at offset 0, an 8-byte arena-end
cell at offset 8, S_node = 2^57 eight-byte head slots from offset 16 (zeroed
once at start-up) and a 40-byte-entry arena that starts after the head slots,
so every arena offset is nonzero. A head is the offset of the node's first
entry in the arena, or 0. The home bucket is floor(h64 * S_node / 2^64) with
h64 = bits 35..98 of K2:K1 (`extr` then `umulh`). For each new point the
generator loads the bucket head into x9, reserves a fresh entry at the bump
pointer, writes new->next = old head, installs the entry as the new head and
publishes the bump pointer, and only then walks the chain from the saved old
head x9, comparing (K1, K2) and raising a **candidate** on equality. The new
entry is not on the scanned chain, so every comparison is between the new point
and an older one and no point is compared with itself. Nothing is deleted, so
the table holds every point. The global load is lambda = n / S_tot =
2^-1.8717 = 0.2733 (S_tot = 2^71 S_node = 2^128 slots).

**Candidates.** A candidate gives two locators (w, j), (w', j'). The
verification action recomputes both points from the directory, compares their
full keys (false candidates stop here), finds u with Q_{w,j} = [u]Q_{w',j'}
(section 1), and if b_w - u b_w' ≠ 0 sets

    x = (u (c_w' + j' D_w') - c_w - j D_w) / (b_w - u b_w')  (mod r),

checks [x]P0 = T and outputs x. Each candidate is priced at 2^14 units
(four scalar multiplications plus keys, unit search and the check: < 2^13.2).
A global 128-bit candidate counter, checked immediately before each
verification, caps the number of verifications at 2^90, i.e. 2^104 units.

**Random-cost allotment and abort, checked per action.** The only
data-dependent costs are (i) chain comparisons in the tables,
(ii) candidate verifications and (iii) zero-denominator recovery actions. The
attack is **one sequential program** (section 6): it keeps a 128-bit
inspection counter x27:x26 with the high-word global cap x28 and a 128-bit
candidate counter x30:x29 with the high-word cap x18 = 2^26 (at most 2^90
verifications). These are ordinary registers, not per-node counters;
the 2^39 generators and 2^71 final nodes are physical banks visited in
sequence, so the same counter covers all of them and no atomic operation is
involved. The first
instructions of every chain comparison are
`adds x26,x26,#1; adc x27,x27,xzr; cmp x27,x28; b.hs abort` (listing F, all
four instructions are counted in the 10-instruction per-comparison price), so
the program halts *before* the comparison that would exceed the global
allotment; every candidate is likewise counted
(`adds x29,x29,#1; adc x30,x30,xzr; cmp x30,x18; b.hs abort`) and compared
with its cap *before* the 2^14-unit verification runs. A zero-denominator
batch is recovered by its 2^17-unit scan only while the
generator's recovery counter is below 4, checked before the scan; a generator
makes < 2^63.4 batches, each zero with probability < 2^24 * 3 / r < 2^-228.8,
so under H' some generator exceeds 4 recoveries with probability < 2^-120.
Hence every run does at most the deterministic work, the global inspection
allotment, the global candidate allotment, at most 2^39 x 4 x 2^17 = 2^58 units
of recovery work, and no round-end aggregation is used: no part of the schedule can spend
more random work than the global allotments. The global
allotments are exceeded with probability < 2^-31.01 (inspections, section 4)
and < 2^-31.01 (candidates, section 4), so they are failure events that do not
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
   < 2^-125 in total).

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
arena overflow 2^-400, global inspection-allotment overrun 2^-31.01, global candidate-allotment overrun 2^-31.01, sampler cap 2^-371,
recovery cap 2^-120; the last four use H' or the explicit caps):

**Success >= 0.390001 - 2^-30.01 > 0.39 for every x.**

Collisions inside one walk (5 of the 6 units also give x there) and the
start points' collisions are extra and only help; they are ignored in the
bound. No property of the curve beyond Lemma 1 (group structure, the unit
group of order 6 and the complete invariant K) is used: the collision law is
exact algebra, not a random-mapping assumption.

**H' (heuristic, cost side only).** The digit1, digit2 and home-bucket bits of
K of distinct classes behave as independent uniform values. It is used for: the
number X of same-bucket comparisons and its multiplicity-aware variance
(below), the false-candidate count, and the per-node arena occupancy
(Chernoff). It never enters the probability that a
collision exists.

**Global inspection tail (multiplicity-aware, elementary).** H' assigns
independent uniform bucket values to distinct *classes*; the table stores
points, and a true class collision forces its points into one bucket, so the
point occupancies are not multinomial. Condition on the class partition
(multiplicities c_1,...,c_m, which depend only on the walk coins, not on the
buckets) and write X = X_same + X_distinct, where X_same = sum_a C(c_a,2)
counts forced same-class point pairs and X_distinct = sum_{a<b} c_a c_b
1[B_a = B_b] counts comparisons between distinct classes. First, E[X_same] =
sum over point pairs P(same class) <= 6 C(n,2)/r = 0.4947, exact by linearity
over all point pairs (the fixed step D_w is independent of the walk start and
of the other walk's coins, so no independence between pairs is needed), and
Markov's inequality gives P(X_same >= 2^30) <= 0.4947/2^30 = 2^-31.02.
Second, given the
partition the B_a are iid uniform on S = 2^128 buckets, so
E[X_distinct] = (n^2 - sum c_a^2)/(2S) <= n(n-1)/(2S) = 2^123.2567; for three
distinct classes P(B_a = B_b = B_c) = 1/S^2 = P(B_a=B_b) P(B_a=B_c), so the
three-class covariance is exactly zero and disjoint pairs are independent;
hence Var(X_distinct) = p(1-p) sum_{a<b} (c_a c_b)^2 <= c_max^2 E[X_distinct].
A class with 8 points occurs with probability <= C(n,8)(6/r)^7 < 2^-771,
so c_max <= 8 and Var <= 64 * 2^123.2567 = 2^129.2567 except with that
negligible probability. The cap is 1.005 E[X_distinct] = 2^123.264, stored in
the high word as floor(cap/2^64) = 2^59.264, so the effective cap is at least
1.0049 E[X_distinct]; after the < 2^30 forced same-class comparisons the slack
is at least 0.0049 E[X_distinct] = 2^115.5. The same-class Markov term above is
part of the same counter, so we bound the two parts separately, with no
independence between same-class pairs assumed:

    P(X > cap) <= P(X_same >= 2^30) + P(X_distinct > cap - 2^30)
               <= 2^-31.02 + 2^129.2567/(0.0049 E[X])^2
                < 2^-31.01.

The earlier iid-point multinomial moments are not used; the mean bound
n(n-1)/(2S) and the 1.10 allotment are unchanged.

**Global candidate tail (multiplicity-weighted Chebyshev).** A candidate is a
pair of points whose 164-bit fingerprints agree: the 128 compared key bits
(K1,K2) together with digit1 (K0 bits 1..36), because equality of K1,K2
already fixes digit2 and the home bucket. Under H' the fingerprint F_a of a
class a is uniform on 2^164 values and independent across classes. Let
c_1,...,c_m be the class multiplicities (sum c_a = n) and define the
class-weighted number of distinct-class candidate pairs

    Z = sum_{a<b} c_a c_b 1[F_a = F_b],

so a class with c points contributes c times each of its points and every
ordered point pair is counted once. Linearity gives

    E[Z] = sum_{a<b} c_a c_b / 2^164 = (n^2 - sum c_a^2)/(2*2^164)
        <= C(n,2)/2^164 = 2^87.2567.

For the variance, condition on the class partition. For two distinct class
pairs sharing one class, P(F_a=F_b=F_c) = 2^-328 = P(F_a=F_b) P(F_a=F_c), so
their covariance is exactly zero; disjoint class pairs are independent.
Hence

    Var[Z] = 2^-164 (1 - 2^-164) sum_{a<b} (c_a c_b)^2
          <= c_max^2 E[Z] <= 64 * 2^87.2567 = 2^93.2567,

using c_max <= 8 except with probability < 2^-771 as above. Chebyshev at
the cap 2^90, after the < 2^30 same-class candidates, gives

    P(Z >= 2^90 - 2^30) <= 2^93.2567 / (2^90 - 2^30)^2 < 2^-86.28.

With P(X_same >= 2^30) <= 2^-31.02, the total number Y = X_same + Z of
candidates satisfies

    P(Y >= 2^90) <= 2^-31.02 + 2^-86.28 < 2^-31.01,

so the global candidate cap aborts with probability below 2^-31.01. No
comparison index is used and no independence between comparisons is assumed:
the dependence between pairs sharing a class is handled exactly by the zero
covariance, and the union with the same-class Markov term covers the
non-singleton classes. The one true candidate is far below the cap.

## 5. Cost

All instruction counts are scalar AArch64 instructions; the reference
operation (one batched incomplete affine addition, 6M) is 1 unit and is never
re-counted. Listings are counted by `experiment.py` from the text below; the
experiment counts each body length and charges the loop control separately.

**128x unrolling of the single generator/insertion pass.** The pass processes a
long fixed-shape stream, so it is emitted as 128 labelled copies of its body;
the 2-instruction loop control (`subs x17,x17,#1; b.ne group`) runs once per
128 points. The charge is therefore body + 2/128 = body + 0.015625 per point.
This is the convention already accepted on the P-256 track (128 labelled
copies, end test once per group); the code size is 128 x 25 instructions
(~13 KB), which is memory, not work. Rare exits (candidate, abort) use a 2-instruction
trampoline (`adr x21,cont; b target`); x21 is the dedicated trampoline link,
so the candidate cap x18 and the group counter x17 survive. Its cost is inside
the separately priced action, while the one-instruction conditional test stays
in the body.
No other load, store, address computation or branch is omitted: every
instruction of the body is listed below and counted.

**Per point, operand traffic and batch bookkeeping (30).** Operand traffic of
the reference addition, charged as in the accepted rho ledgers: 6 `ldp`
current point, 6 `ldp` step point, 6 `stp` new point (18). Pass bookkeeping of
the Montgomery-batched addition, charged explicitly although the rho ledgers
treat it as part of the reference operation: the forward pass stores each
prefix product (3 `stp`) and the backward pass reloads it (3 `ldp`), each pass
has `cmp`/`b.ne` loop control (4; addresses advance by pre/post-indexed
addressing), plus 2 spare (12). The field multiplications and the A = 0
additions of the two passes are the reference operation itself.

**Per point, key and insertion: listing G (25 + 2/128) plus a 2-instruction
bank select.** The generator has the new point's y in its record; it computes
K0..K2, the 71-bit node id, selects the final node's bank (section 6), finds
the home bucket, links the entry and scans the bucket's old chain:

```
gen:  ldp  x0, x1, [x19, #REC]!      // y limbs 0,1 of the new point; record pointer advance
      ldp  x2, x22, [x19, #16]       // y limb 2 and the walk id (locator, per point)
      subs x3, x10, x0               // limbs 0-2 of p - y (x10-x12 = limbs of p)
      sbcs x4, x11, x1
      sbcs x5, x12, x2
      tst  x0, #1                    // ne iff y odd
      csel x0, x3, x0, ne            // K0
      csel x1, x4, x1, ne            // K1
      csel x2, x5, x2, ne            // K2
      ubfx x6, x0, #1, #36           // digit1 (bank-select low half)
      ubfx x7, x1, #0, #35           // digit2 (bank-select high half)
      extr x8, x2, x1, #35           // h64 = bits 35..98 of K2:K1
      umulh x8, x8, x20              // home bucket = floor(h64 * S_node / 2^64)
      ldp  x15, x13, [x25]           // node arena bump pointer and arena end (per point)
      add  x8, x25, x8, lsl #3       // head-cell address minus 16 in the selected bank
      add  x16, x15, #40             // new bump
      cmp  x16, x13                  // arena capacity guard
      b.hs abort                     // defined failure path, Chernoff < 2^-400
      ldr  x9, [x8, #16]             // old head, or 0
      str  x9, [x15]                 // new entry: next = old head
      str  x15, [x8, #16]            // head cell = new entry (linked before the scan)
      stp  x1, x2, [x15, #8]         // K1, K2
      stp  x22, x23, [x15, #24]      // walk id, j
      str  x16, [x25]                // publish the new bump pointer
      cbz  x9, next                  // empty bucket: no scan
```

The bank-select macro (`stp x6,x7,[x24]; dsb sy`, 2 instructions: one 16-byte
store to the 128-bit bank-select register, which takes digit1 in the low half
and digit2 in the high half, plus one `dsb sy`) is inserted between the two
`ubfx` and `extr`; it is counted separately below. x20 = S_node, x25 = the bank
window base, x24 = the MMIO bank-select register, x18 = the candidate-cap high
word, x22/x23 = walk id and j of the point, x21 = the rare-exit link.

**Register lifetimes.** x10, x11 and x12 hold the three low limbs of p and are
only read by `subs`/`sbcs`; no listed write targets x10-x12, so every later
point's key rows use the same p limbs. x8 (head cell minus 16) is written only
by the `add` before the scan and is never written inside the scan, so the head
cell address is intact when `str x15, [x8, #16]` installs the entry. x9 is the
saved old head and then the traversal pointer: the old head is copied into the
new entry (and the head cell updated) *before* x9 is advanced, so the scan
cannot clobber it. x13 (arena end) is not assumed loop-invariant across points:
`ldp x15, x13, [x25]` reloads the destination node's header on every point.
x22 (walk id) is reloaded from the point's record on every point by the
`ldp x2, x22` above; x23 = j is constant within a batch (one batch = one value
of j). The record keeps the walk id at offset REC+24. Part A emulates these
transitions (`emulate_listing_GF`): the p limbs survive every iteration, the
prepend links are consistent, the stored (walk id, j) locator recovers the
planted point, and the bank header is read per point.

**Per chain entry compared (10).** The scan compares (K1, K2), advances the
128-bit inspection counter and checks it *before* the comparison, and steps to
the next entry. It traverses with x9 and never writes x8 (the head cell),
x10-x12 (p limbs), x1/x2 (the new key) or x19 (the record pointer):

```
scan: ldp  x3, x4, [x9, #8]          // K1, K2 of the stored entry
      adds x26, x26, #1              // 128-bit inspection counter, low word
      adc  x27, x27, xzr             // high word
      cmp  x27, x28                  // immediate global cap, high word
      b.hs abort                     // halt before the comparison that would exceed it
      cmp  x3, x1
      ccmp x4, x2, #0, eq            // (K1,K2) equal?
      b.eq cand                      // candidate action (separately priced)
      ldr  x9, [x9]                  // next entry (offset 0)
      cbnz x9, scan
```

**Per candidate (separately priced, 4 counter instructions plus 2^14 units).**
Every candidate is preceded by

```
      adds x29, x29, #1
      adc  x30, x30, xzr
      cmp  x30, x18                  // x18 = 2^26 = floor(2^90 / 2^64)
      b.hs abort
```

and the verification recomputes both points from the directory, compares their
full 192-bit keys, finds the unit, solves and checks (priced at 2^14 units,
< 2^13.2); a true candidate is found once. The cold-path trampoline into and
out of the handler is inside this price.

**Per point total (deterministic).** 18 + 12 + (25 + 2/128) + 2 (bank select) =
57.015625 instructions; the 25 already include the per-point walk-id load
(`ldp x2, x22`), the per-point arena header load (`ldp x15, x13`) and the
arena-capacity guard, plus the per-batch work (Fermat inversion 380
squarings + 228 multiplications at 383 = 232,864, exponent scan < 2,000,
zero-denominator test on the batch product, per-batch constants: < 235,000 per
batch of 2^24 points, i.e. < 0.014008) and a 0.002 per-point batch-dispatch
allowance. The total is 57.031632 instructions per point. Point price:
1 + 57.031632/2000 = 1.028515816 units. (Start points are also charged this
price, on top of their set-up, which over-counts.)

**Table zeroing.** S_tot = 2^128 eight-byte heads, zeroed with `stp xzr, xzr`
(2 heads per store) in a loop of 4096 stores per 2 loop instructions:
0.500244 instructions per head, i.e. 2^116.035 units in all.

**Chain comparisons (random, allotted).** The expected number of comparisons
is E[X] = n(n-1)/(2 S_tot) = 2^123.2567 (0.1366 per point): insertion scans
each bucket's existing chain, so the total is exactly sum_b C(m_b, 2). Each
comparison is charged 10 instructions (the 10 listed above). The global
allotment is 1.005 x E[X], enforced by the 128-bit counter x27:x26 with the
high-word cap x28 = floor(2^123.2640/2^64) = 2^59.264 checked before every
comparison. The multiplicity-aware variance bound (section 4) gives the
overrun probability < 2^-31.01 (Markov on X_same plus Chebyshev on
X_distinct), so the every-run bound holds with that one failure probability.

**Other terms.** Walk set-up (capped sampling and three scalar multiplications
per walk, plus directory writes) < 2^14 units x 2^63 = 2^77; node-bank
initialisation (2^71 headers: bump and arena-end cells, 3 instructions including the select) and
MMIO setup < 2^76; candidates: cap 2^90 verifications x 2^14 = 2^104 units,
enforced before each verification by the 128-bit counter (mean false
candidates 2^87.26, overrun 2^-31.01); capped zero-denominator recovery
2^39 x 4 x 2^17 = 2^58; bank-select device setup and node-visit selects
< 2^104 units (2^73 setup instructions = 2^62 units, well inside); final
output and verification < 2^20.

**Total.**
W <= n * 1.028515816 + 2^116.035 + 2^115.620 + 2^104 (candidates)
     + 2^104 (address mapping) + 2^77 + 2^76 + 2^58 + 2^20
< 2^126.171145 <= 2^126.17115 = W_max.

For orientation: the two-hop routing pipeline of the earlier packages (with
its intermediate pass, receive windows and window-padding term) is replaced by
direct insertion into the destination node's bank. The exact ledger computed
by `experiment.py` is W < 2^126.171145.

## 6. Memory and machine model

Peak storage < 2^132.3 bytes: heads S_tot x 8 B = 2^131; final-node headers
(bump pointer + arena end) 2^71 x 16 B = 2^75; entry arenas 1.05 n x 40 B =
2^131.52; walk records (504-byte stride, so every pre-indexed `ldp` offset,
including `#REC` = 504, is encodable; the walk id sits at REC+24) and the walk
directory (2^63 x 600 B = 2^72.3); total < 2^132.29 (including a 2^76 B
allowance for address tables). Every bank is below 2^61.4 B: a final node holds its arena
bump pointer and arena-end cell (16 B), S_node x 8 B = 2^60 B of heads and a
1.05 x 2^55.13 x 40 B = 2^60.54 B arena, 2^61.06 B in all (+16 B); a generator
node holds its 2^24 x 600 B = 2^34.5 B of records. Every load/store address is 64-bit.

**Banked 64-bit address space, bank select specified and charged.** Physical
storage is 2^132 B, more than one 64-bit address space, so the program reaches
it through a 2^64-byte virtual address space with two regions: a permanently
mapped local region of at most 2^62 B holding the current generator's records,
and a bank window of 2^61.4 B at a fixed virtual base x25. The 2^71 banks are
indexed by the pair (digit1, digit2); a physical bank is installed in the
window by one 16-byte store `stp x6, x7, [x24]` to the memory-mapped 128-bit
bank-select register at a fixed device address (digit1 in the low half, digit2
in the high half; one AArch64 instruction) followed by one `dsb sy` (one
instruction) that orders the previous window accesses before the remap and
makes the new mapping visible; the device concatenates the two 64-bit halves
into the 71-bit bank index, no prefix or context state is needed. Both are real AArch64
instructions and are counted; the device performs only the address
translation, no computation. A bank is selected afresh before every insertion
(one hop per point, 2 instructions charged in section 5), so no mapping is
cached, preserved or restored: digit1 and digit2 are extracted arithmetically
(the two `ubfx` of listing G) and the selected bank is the one the next head
load and entry stores need. The same 2-instruction select visits all 2^71
final banks during start-up (the 2^71 zeroing/bump-initialisation visits below). The working set fits at every step: the generator keeps its
own records mapped and one final node's bank (2^61.06 B < 2^61.4 B) in the
window. Startup selects (2^71 nodes visited to zero the heads and initialise
the bump cells, < 3 x 2^71 = 2^72.6 instructions = 2^61.6 units) are inside the
`address_mapping` term of 2^104 units in section 5. No cross-bank transfer
occurs: the entry is written once into its destination bank and read there
during later chain scans; the chain pointers are 64-bit virtual addresses of
entries in the same bank. No address exceeds 64 bits.

**One sequential program.** The attack is executed by one sequential processor
using the virtual address space above; the 2^39 generators and 2^71 final
tables are physical banks visited in sequence, and the "remote" loads and
stores are ordinary loads and stores into the bank window after the
2-instruction select. The global inspection counter x27:x26 with cap x28 and
the candidate counter x30:x29 with cap x18 are ordinary registers of this one
program, updated by four instructions per action; no atomic, network or
per-node counter is involved. The work is the sum over the sequential
schedule, exactly the same instructions as the distributed description, so the
score is unchanged.

Memory has no scalar contribution under classical-group-work-v2; it is stated
here as required. Physically 2^132 bytes is far beyond any machine, as is the
2^126-operation time of every attack on this track; the score is the
instruction-equivalent count, which this design charges completely: each point
is generated, keyed and inserted by explicit loads and stores, and every
transfer is the store instruction that issues it.

**No network.** There is no interconnect or message layer: a remote store is
an ordinary `str` into the bank window after the select, so the destination's
physical memory is written directly. Every transfer is the store instruction
that issues it plus the 2-instruction select counted above; no further copy,
hashing, sorting or comparison is assumed. Routing is arithmetic (the node id)
and is in listing G.

## 7. Experiments (`experiment.py`, stdlib)

* **Part A (production, all true):** P0 on the curve, [r]P0 = O, beta and
  lambda, phi(P0) = [lambda]P0, the unit map on all six members, key constant
  on the six members and even, key completeness, listing G's key rows emulated
  limb by limb against K mod 2^192, digit1, digit2 and the 71-bit node id
  (random points, both signs, and 16 borrow-chain boundary values), the final home index in range, AArch64
  immediate encodability (`and #0x7ffffffff`, `tst #1`), and a **planted
  end-to-end run on BLS12-381**: three walks with real batched affine
  additions, a planted class collision between two walks at a random unit,
  the node routing and a chained table built from the emulated rows (8-byte
  head buckets, 40-byte `(next, K1, K2, w, j)` entries, prepend and chain scan),
  the candidate, re-computation from the directory, the unit search, the solving
  formula, and [x]P0 = T.
* **Part B (toy curve y^2 = x^3 + 4 over F_q with a prime-order subgroup and
  the same six units):** the whole attack at toy size with the chained table,
  and the exact scalar model at small r. A run of this package's code (q30,
  32 walks, L=55) gives 0.4176 +- 0.0124 over 1,578 trials against the
  cross-walk prediction 0.3907, and the scalar model gives 0.4371 +- 0.0044
  over 12,544 trials against 0.4318. The toy runs are evidence for the
  implementation, not for the bound, which is proved in section 4.
* **Part C:** the ledger above in exact fractions: L from f4, the price, the
  zeroing and probing terms, the extras, W, the dependency-bound exponents, the
  multiplicity-aware inspection moments (X_same, X_distinct, Var <= 64 E[X],
  c_max <= 8) and the multiplicity-weighted Chebyshev candidate tail, the
  128-bit counter listings, the bank-select charge,
  the candidate cap 2^90 = 2^104 units, the 1.05 arena margin and its charged guard, the 128x unroll
  accounting, and checks that every listing's counted length equals its ledger
  entry.

## 8. Limitations

* Memory is enormous (2^132.3 bytes) and is not scored by the rules.
* Remote stores are bank-window stores; one hop per point pays one 16-byte
  `stp` to the 128-bit bank-select register plus one `dsb`, 2 instructions per
  point in total (section 6).
* H' (uniform key bits) is a heuristic; the abort, overflow and false-candidate
  bounds rely on it (each with a stated negligible probability, and a 0.5%
  slack on the global inspection allotment). The inspection tail is the
  multiplicity-aware finite Chebyshev bound of section 4 (conditioned on the
  class partition; no asymptotic variance law and no iid-point assumption);
  the collision law itself is proved and does not use H'.
* The per-node arena capacity is 1.05 times the expected per-node count. The
  insertion listing checks the new bump pointer against the node's arena end
  before storing (`cmp`/`b.hs abort`, 2 charged instructions); the abort is a
  defined failure path, and a node reaching it is a Chernoff event below
  2^-400 under H' (failure term `arena_overflow`).
* Instruction counts are analytic listings under the fixed C = 2,000, not
  hardware timings. The reference-operation operand traffic is charged as in
  the earlier accepted rho ledgers (18 transfer instructions per addition).
  The single sequential program's global counters are ordinary registers;
  their four-instruction checks are priced inside the 10-instruction comparison
  price (and four inside the separately priced candidate action).

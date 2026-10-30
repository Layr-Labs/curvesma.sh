# Distributed fixed-step walks with a 64-bit-local class table on BLS12-381 G1

**Claim.** time_log2 = 126.18013. For every challenge T = [x]P0 (hence for a
uniform one) the algorithm below outputs x with probability at least 0.39 over
its own coins, and its total work is at most W_max = 2^126.18013 units (one unit =
one batched incomplete affine G1 addition = C = 2,000 scalar AArch64
instructions). The computed total is W = 2^126.1801271 (< 2^126.18013); the
claim keeps a margin of 0.000003 in log2. The work bound holds on every run:
every data-dependent action is checked against its allotment immediately before
it runs (the two shared 128-bit counters of section 2 and the explicit attempt
caps of the samplers), and the only overshoot of the shared inspection counter
is the deterministic 2^87-comparison block lag, which is inside the probe
allotment. Peak memory < 2^132.284 bytes (section 6); no advice; preprocessing
(walk set-up, routing tables, server headers, table zeroing) < 2^116.1 units,
all inside W.

**What is different from the six-class rho line of submissions.** Those walks
must be functions of the class, so every step computes the canonical twist index
and runs cycle handling; their success law is the random-mapping heuristic H1
with a 0.5% allowance. Here the walks are plain arithmetic progressions with a
per-walk fixed step, so no twist index, no cycles and no distinguished points
are needed; instead every point's class key is inserted into a **distributed
direct-bucket chained** table that detects any two points in the same class.
The success probability is **proved** (section 4) from exact linear-algebra
probabilities of the collision events over the algorithm's coins and the
Bonferroni inequality; no random-mapping heuristic is used. The one heuristic
left (H', key bits behave as uniform) is used only for the cost-side bounds
(the expected chain scan and its exact variance, false candidates, arena
occupancy) and never enters the probability that a collision exists.
Precisely: the collision law is proved; the abort and overflow bounds use H'.
The price is memory: the table holds all ~2^126.13 points (memory is reported
but is not part of the score under the track rules).

**What this package repairs.** v8/v9 used a single machine with a "128-bit MMIO
bank-selector register" to reach 2^131 B; the exploratory review held that the
protected model does not supply that device and that counting its two issuing
instructions does not establish a scalar-AArch64-equivalent implementation.
This package removes the device entirely: the table is distributed over 2^71
table servers, each with a local bank below 2^61.4 B, so every head word,
arena offset and next pointer is an ordinary 64-bit local address and every
server instruction is an ordinary scalar AArch64 instruction. The transfer is
instead charged as an explicit per-point message: construction and enqueue,
two network hops and receive (section 5), with all global counters as shared
128-bit atomic words. The earlier acceptance of any related package is not
offered as support; the argument is self-contained.

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

**Distributed insertion (one pass, no bank selector).** The table holds all
2^126.13 points, so its ~2^131 B cannot be addressed by one 64-bit address
space. The machine is therefore distributed: 2^37 generator nodes each hold a
local batch of 2^26 walk records and advance them in lock step: one
Montgomery-batched inversion per batch (one batch = one value of j), and the
per-batch share is 235,000/2^26 instructions per point. After a batch has
stored its new points, the generator loop (listing GEN_key, section 5)
computes K0, K1, K2 for each new point and forms the 71-bit routing id
digit1 || digit2 (digit1 = bits 1..36 of K0, digit2 = bits 0..34 of K1). It
then writes a 56-byte message (K1, K2, walk id, j, digit1, digit2) into its
outgoing queue (4 instructions). A two-level network routes each message to
one of 2^71 table servers by the 71-bit id (3 instructions per level: load the
destination id, load the forwarding entry, branch). The destination server
dequeues the 32-byte key/locator payload (6 instructions) and inserts it into
its own local bank (13 instructions: home bucket, head cell, arena guard,
link-before-scan prepend, payload stores). There is no intermediate node, no
receive window, no second pass and no remote-store primitive: the message
passing and every network hop are the transfer charge, and all instructions of
every stage are listed and counted in section 5. Every address on a table
server is an ordinary 64-bit local address; **there is no selector register,
no dsb, no bank window and no address_mapping term**. Work is summed over
processors; memory is reported but not scored.

**Final table (listing F).** Each of the 2^71 table servers keeps its own
local table: an 8-byte arena bump pointer at offset 0, an 8-byte arena-end
cell at offset 8, a receive-queue head at offset 16, S_node = 2^57 eight-byte
head slots from offset 32 (zeroed once at start-up) and a 40-byte-entry arena
that starts after the head slots, so every arena offset is nonzero. A head is
the offset of the server's first entry in its arena, or 0. The home bucket is
floor(h64 * S_node / 2^64) with h64 = K2, the top limb of the class key
(umulh alone; every compared entry is checked on (K1,K2)). For each new point
the server loads the bucket head into x9, reserves a fresh entry at the bump
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
(section 1), and if b_w - u b_w' != 0 sets

    x = (u (c_w' + j' D_w') - c_w - j D_w) / (b_w - u b_w')  (mod r),

checks [x]P0 = T and outputs x. Each candidate is priced at 2^14 units
(four scalar multiplications plus keys, unit search and the check: < 2^13.2).
A 128-bit candidate counter in shared memory, checked by two atomic adds
immediately before each verification, caps the number of verifications at 2^90,
i.e. 2^104 units.

**Random-cost allotment and abort, checked per action.** The only
data-dependent costs are (i) chain comparisons, (ii) candidate verifications
and (iii) zero-denominator recovery actions. Because the table is distributed,
the two global counters are 128-bit words in shared memory, updated with
atomic scalar AArch64 instructions (ARMv8.1 ldadd; no lock, no retry loop):

* Each table server keeps a 64-bit local comparison counter x26. Every
  comparison increments it (adds). When it reaches a multiple of BLK = 2^16
  (tst x26, #0xffff, one instruction per comparison) the server adds its block
  to the shared low word (ldadd x27, [x30], x31), resets x26 and tests the
  returned old value for a low-word wrap (cmp, b.lo). The unique server that
  observes a wrap adds 1 to the shared high word (ldadd x28, [x32], x33) and
  aborts when the new high word reaches cap_hi = 692151372792945841 =
  2^59.2639 (cmp, b.hs abort). This is 4 charged instructions per 2^16
  comparisons plus 3 per 2^64 for the carry.
* The shared counter lags the true comparison count by at most one partial
  block per server, 2^71 * 2^16 = 2^87 comparisons, which is included in the
  deterministic probe allotment of section 5.
* Every candidate is counted by the shared 128-bit candidate counter with two
  ldadds and two cap checks (6 instructions; cap 2^90 verifications). A zero
  denominator is recovered by its 2^17-unit scan only while the server's
  recovery counter is below 4, checked before the scan; a server makes fewer
  than 2^63.4 batches, each zero with probability < 2^26 * 3 / r < 2^-226.8,
  so under H' some recovery counter exceeds 4 with probability < 2^-120.

Hence every run does at most the deterministic work, the global inspection
allotment (plus the 2^87 deterministic lag), the global candidate allotment,
at most 2^58 units of recovery work, and no round-end aggregation is used: no
part of the schedule can spend more random work than the global allotments.
The global allotments are exceeded with probability < 2^-31.015 (inspections,
section 4) and < 2^-31.015 (candidates, section 4), so they are failure events
that do not threaten the 0.39 success bound.

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
n(n-1)/(2S) and the 1.005 allotment are unchanged.

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
re-counted. Listings are counted by experiment.py from the text below; the
experiment counts each body length and charges the loop control separately.
For listing F the text has 16 lines: 9 are the common path charged for every
comparison, 4 are the block-update path charged 4*2^-16 per comparison and 3
are the carry path charged 3*2^-64 per comparison (the experiment checks
exactly this split).

**128x unrolling of both stream passes.** The generator and the table server
process long fixed-shape streams, so each pass is emitted as 128 labelled
copies of its body; the 2-instruction loop control runs once per 128 points on
each side, charged as 2*2/128 = 0.03125 per point. Rare exits use a
2-instruction trampoline whose cost is inside the separately priced action; the
one-instruction conditional test stays in the body. No other load, store,
address computation or branch is omitted.

**Per point, operand traffic and batch bookkeeping (30).** Operand traffic of
the reference addition, charged as in the earlier rho ledgers: 6 ldp current
point, 6 ldp step point, 6 stp new point (18). Pass bookkeeping of the
Montgomery-batched addition: the forward pass stores each prefix product
(3 stp) and the backward pass reloads it (3 ldp), each pass has cmp/b.ne loop
control (4) plus 2 spare (12). The field multiplications and the additions of
the two passes are the reference operation itself.

**Per point, generator key rows: listing GEN_key (11).**

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
      ubfx x6, x0, #1, #36           // digit1
      ubfx x7, x1, #0, #35           // digit2
```

x19 is the generator's record pointer, x10-x12 hold the three low limbs of p,
x22 = walk id and x23 = j. The generator and the server are different
processors; the key rows only read the p limbs and never write them.

**Per point, message and network: listings MSG (4) and NET (6).**

```
msg:  stp  x1, x2, [x14]             // K1, K2 into the outgoing queue slot
      stp  x22, x23, [x14, #16]      // walk id, j
      stp  x6, x7, [x14, #32]        // digit1, digit2 (routing descriptor)
      add  x14, x14, #56             // advance the queue slot pointer
net:  ldr  x9, [x13, #32]            // level 1: destination digit1
      ldr  x11, [x12, x9, lsl #3]    // forwarding entry
      br   x11
      ldr  x9, [x13, #40]            // level 2: destination digit2
      ldr  x11, [x12, x9, lsl #3]
      br   x11
```

The message is 56 bytes (K1, K2, walk id, j, digit1, digit2, spare). x13 is
the message address read by the fabric, x14 the generator's queue pointer and
x12 the per-level routing table; the queue tail is published once per batch
inside the 2/128 loop-control charge.

**Per point, server receive: listing RCV (6).** The destination server pops
the message from its local inbox:

```
rcv:  ldp  x1, x2, [x9]              // K1, K2
      ldp  x22, x23, [x9, #16]       // walk id, j
      add  x9, x9, #56               // advance the inbox pointer
      cmp  x9, x13                   // end of the inbox ring?
      b.ne next
      mov  x9, x10                   // wrap to the ring base
```

**Per point, server insert: listing SRV_insert (13).**

```
srv:  umulh x8, x2, x20              // home bucket = floor(K2 * S_node / 2^64)
      add  x8, x25, x8, lsl #3       // head-cell address
      ldp  x15, x13, [x25]           // arena bump pointer and arena end (per point)
      add  x16, x15, #40             // new bump
      cmp  x16, x13                  // arena capacity guard
      b.hs abort                     // defined failure path, Chernoff < 2^-400
      ldr  x9, [x8, #32]             // old head, or 0
      str  x9, [x15]                 // new entry: next = old head
      str  x15, [x8, #32]            // head cell = new entry (linked before the scan)
      stp  x1, x2, [x15, #8]         // K1, K2
      stp  x22, x23, [x15, #24]      // walk id, j
      str  x16, [x25]                // publish the new bump pointer
      cbz  x9, next                  // empty bucket: no scan
```

The server header at x25 is bump pointer (offset 0), arena end (8), receive
queue head (16); the head slots start at offset 32. x20 = S_node; x22/x23 =
walk id and j of the point; x21 = the rare-exit link. The server listing never
writes the generator's x10-x12; x8 (head-cell address) is written only before
the scan; x9 is the saved old head and then the traversal pointer, so the
link-before-scan cannot be clobbered; x13 is the arena end reloaded per point.
Part A emulates the split data flow (emulate_listing_GF): the p limbs survive
every iteration, the prepend links are consistent, the stored locator recovers
the planted point, and the server header is read per point.

**Per chain entry compared (9 + 4*2^-16 + 3*2^-64).** The scan counts and
checks the shared inspection counter and compares (K1, K2):

```
scan: ldp  x3, x4, [x9, #8]          // K1, K2 of the stored entry
      adds x26, x26, #1              // local comparison counter
      tst  x26, #0xffff              // block boundary every 2^16 comparisons
      b.ne scanc
      ldadd x27, [x30], x31          // add BLK to the shared low word; old value returned
      mov  x26, xzr                  // reset the local counter
      cmp  x27, x29                  // low-word wrap? (x29 = 2^64 - BLK)
      b.lo scanc
      ldadd x28, [x32], x33          // carry: add 1 to the shared high word
      cmp  x28, x18                  // new high word reaches cap_hi?
      b.hs abort                     // halt before the comparison that would exceed the allotment
      cmp  x3, x1
      ccmp x4, x2, #0, eq            // (K1,K2) equal?
      b.eq cand                      // candidate action (separately priced)
      ldr  x9, [x9]                  // next entry (offset 0)
      cbnz x9, scan
```

The common path is 9 instructions and runs for every comparison; the block
path is 4 instructions every 2^16 comparisons and the carry path 3 every
2^64. The shared counter lags by at most 2^71 partial blocks = 2^87
comparisons, which is inside the probe allotment. F never writes x8 (the head
cell) or the new key x1/x2.

**Per candidate (separately priced, 6 counter instructions plus 2^14 units).**
Every candidate is preceded by two 64-bit atomic adds to the shared 128-bit
candidate counter and two cap checks:

```
      ldadd x29, [x30], x31          // low word
      cmp  x29, x32
      b.hi abort
      ldadd x27, [x33], x31          // high word
      cmp  x27, x18
      b.hs abort
```

and the verification recomputes both points from the directory, compares their
full 192-bit keys, finds the unit, solves and checks (priced at 2^14 units,
< 2^13.2); a true candidate is found once. The cold-path trampoline into and
out of the handler is inside this price.

**Per point total (deterministic).** 18 + 12 + 11 + 4 + 6 + 6 + 13 = 70
instructions; loop control 2*2/128, per-batch work (Fermat inversion 380
squarings + 228 multiplications at 383 = 232,864, exponent scan < 2,000,
zero-denominator test on the batch product, per-batch constants: < 235,000 per
batch of 2^26 points, i.e. < 0.0035) and a 0.002 per-point dispatch allowance.
The total is 70.03675177288055 instructions per point. Point price:
1 + 70.03675177288055/2000 = 1.0350183758864402 units. (Start points are also
charged this price, on top of their set-up, which over-counts.)

**Table zeroing.** S_tot = 2^128 eight-byte heads (2^57 per server), zeroed
with stp xzr, xzr (2 heads per store) in a loop of 4096 stores per 2 loop
instructions: 0.500244 instructions per head, i.e. 2^116.03492 units in all.

**Chain comparisons (random, allotted).** The expected number of comparisons
is E[X] = n(n-1)/(2 S_tot) = 2^123.2567 (0.1366 per point): insertion scans
each bucket's existing chain, so the total is exactly sum_b C(m_b, 2). Each
comparison is charged 9 + 4*2^-16 + 3*2^-64 instructions. The shared 128-bit
counter enforces the high-word cap cap_hi = floor(1.005 E[X]/2^64) =
692151372792945841 = 2^59.2639; the allowed total is
cap_eff = cap_hi*2^64 = 2^123.2639 minus the low word, plus at most 2^87
unflushed comparisons. The multiplicity-aware variance bound (section 4) gives
the overrun probability < 2^-31.015, so the every-run bound holds with that one
failure probability. The charged probe allotment is
(cap_eff + 2^87)*(9 + 4*2^-16 + 3*2^-64)/2000 = 2^115.46801 units.

**Other terms.** Walk set-up (capped sampling and three scalar multiplications
per walk, plus directory writes) < 2^14 units x 2^63 = 2^77; network and queue
set-up (route tables, 2^71 inbox rings) < 2^90; server headers (bump, arena end,
queue head: 2^71 x 32 B) < 2^76; candidates: cap 2^90 verifications x 2^14 =
2^104 units (mean false candidates 2^87.26, overrun 2^-31.01); capped
zero-denominator recovery 2^39 x 4 x 2^17 = 2^58; final output and verification
< 2^20.

**Total.**

```
W <= n * 1.0350183758864402
   + 2^116.03492                                  (zeroing)
   + (cap_eff + 2^87) * (9 + 4*2^-16 + 3*2^-64)/2000   (inspections)
   + 2^104 + 2^90 + 2^77 + 2^76 + 2^58 + 2^20     (extras)
 = 2^126.1801271  <  2^126.18013 = W_max.
```

The exact ledger computed by experiment.py is W = 2^126.1801271.

## 6. Memory and machine model

Peak storage < 2^132.284 bytes: heads S_tot x 8 B = 2^131; server headers
(2^71 x 32 B) = 2^76; entry arenas 1.05 n x 40 B = 2^131.52; generator record
stores and walk directory (2^63 x 600 B = 2^72.3); receive-queue rings
2^71 x 64 B = 2^77; total < 2^132.284. Every table server's bank is below
2^61.4 B: 32 B of header, S_node x 8 B = 2^60 B of heads and a
1.05 x 2^55.13 x 40 B = 2^60.52 B arena, 2^61.284 B in all. A generator node
holds 2^26 x 600 B = 2^35.6 B of records. Every load/store address in every
listing is a 64-bit local address; no region is reached through a device and
no address is assembled with movk.

**Distributed machine.** There are 2^37 generator nodes and 2^71 table-server
nodes. Each node owns its memory; the only inter-node traffic is the 56-byte
message per point, whose construction, two network hops and receive are the
4 + 6 + 6 instructions charged in section 5. The network is a two-level
butterfly: level 1 routes on digit1, level 2 on digit2; each level loads the
destination id from the message, loads the forwarding entry and branches. No
data is copied, hashed or sorted outside those charged instructions; the
message payload is the (K1,K2,w,j) locator. There is no MMIO device, no
bank-select register, no dsb and no virtual-memory window; the rejected v8/v9
address mechanism is gone.

**Counters.** The global inspection and candidate counters are 128-bit words
in shared memory, updated by the atomic scalar AArch64 ldadd instructions
listed in section 5 (low word per block, high word on wrap). ldadd returns the
old value, so the unique server seeing a wrap performs the carry; there is no
retry loop and no lock. The per-server local counters are ordinary registers.
The deterministic lag of the shared counter (at most 2^71 partial blocks =
2^87 comparisons) is inside the probe allotment.

Memory has no scalar contribution under classical-group-work-v2; it is stated
here as required. Physically 2^132 bytes is far beyond any machine, as is the
2^126-operation time of every attack on this track; the score is the
instruction-equivalent count, which this design charges completely: each point
is generated, keyed, messaged, routed, received and inserted by explicit
instructions, and every transfer is the charged message path.

## 7. Experiments (experiment.py, stdlib)

* **Part A (production, all true):** P0 on the curve, [r]P0 = O, beta and
  lambda, phi(P0) = [lambda]P0, the unit map on all six members, key constant
  on the six members and even, key completeness, listing GEN_key's key rows
  emulated limb by limb against K mod 2^192, digit1, digit2 and the 71-bit
  routing id (random points, both signs, and 16 borrow-chain boundary values),
  the final home index in range, AArch64 immediate encodability
  (and #0x7ffffffff, tst #1, tst #0xffff), the split data-flow emulation of
  GEN_key -> MSG -> NET -> RCV -> SRV_insert -> F (p limbs preserved, per-point
  walk-id and arena-end loads, link-before-scan, head-cell offset 32,
  chained-prepend consistency and nonzero arena base, locator recovery), and a
  **planted end-to-end run on BLS12-381**: three walks with real batched affine
  additions, a planted class collision between two walks at a random unit, the
  split protocol listings, a chained table built from the emulated server rows
  (8-byte head buckets, 40-byte (next,K1,K2,w,j) entries, prepend and chain
  scan), the candidate, re-computation from the directory, the unit search, the
  solving formula, and [x]P0 = T (21 planted entries).
* **Part A2 (production-bit H' evidence, real curve, all true):** 64
  independent AP walks x 33 points = 2112 real subgroup points, generated with
  batched affine additions and the exact production extraction code. 199
  key/routing bit-balance tests, max |z| = 2.35; 64 pairwise 2x2 chi-square
  tests, max 4.844; Poisson occupancy at S = 2^6 (lambda = 33) chi2 = 2.91 with
  same-bucket pairs 34671 vs 34831.5 and at S = 2^10 (lambda = 2.062)
  chi2 = 5.97 with 2218 vs 2177; truncated 164-bit fingerprint collisions on
  high/middle/low bit slices for b = 12,14,16,18,20, max |z| = 1.887 over 15
  tests; within-walk vs cross-walk at b = 16: 0 vs 0.516 expected within
  (z = -0.718) and 25 vs 33.5 cross (z = -1.469). **This is consistency
  evidence at production bit positions on a scaled sample, not a proof of
  H'.** At m = 2112 the expected truncated-fingerprint collision counts are
  ~2^21/2^b (b = 12: 544, b = 20: 2.1), so tag widths above ~20 have no power
  and are not claimed.
* **Part B (toy curve y^2 = x^3 + 4 over F_q with a prime-order subgroup and
  the same six units):** the whole attack at toy size with the chained table,
  and the exact scalar model at small r. The declared run (q26, 32 walks,
  L=23) gives 0.4089 +- 0.0148 over 1,103 trials against the cross-walk
  prediction 0.4151, and the scalar model gives 0.3995 +- 0.0078 over 3,910
  trials against 0.4104. The toy runs are evidence for the implementation, not
  for the bound, which is proved in section 4.
* **Part C:** the ledger above in exact fractions: L from f4, the
  70-instruction per-point sum, the 1.0350183758864402 step price, the
  2^116.03492 zeroing term, the (cap_eff + 2^87)(9 + 4*2^-16 + 3*2^-64)/2000
  inspection allotment, the extras, W, the dependency-bound exponents, the
  multiplicity-aware inspection moments (X_same, X_distinct, Var <= 64 E[X],
  c_max <= 8) and the multiplicity-weighted Chebyshev candidate tail, the
  shared 128-bit atomic counter listings, the 2^90 = 2^104-unit candidate cap,
  the 1.05 arena margin and its charged guard, the 128x two-sided unroll
  accounting, the 2^63-byte per-server bank bound, the no-selector/no-dsb/
  no-movk listing checks, and checks that every listing's counted length
  equals its ledger entry.

## 8. Limitations

* Memory is enormous (2^132.284 bytes) and is not scored by the rules.
* The machine is distributed: 2^37 generator nodes and 2^71 table servers.
  The only inter-node charge is the message path (4 instructions to form and
  enqueue, 3 per network level, 6 to receive), all listed and counted in
  section 5. No remote-store primitive, no bank window and no
  address-translation device is assumed.
* H' (uniform key bits) is a heuristic; the abort, overflow and
  false-candidate bounds rely on it (each with a stated negligible
  probability, and a 0.5% slack on the global inspection allotment). The
  inspection tail is the multiplicity-aware finite Chebyshev bound of section
  4 (conditioned on the class partition; no asymptotic variance law and no
  iid-point assumption); the collision law itself is proved and does not use
  H'. Part A2 is consistency evidence at production bit positions, not a proof
  of H'.
* The per-server arena capacity is 1.05 times the expected per-server count.
  The insertion listing checks the new bump pointer against the server's arena
  end before storing (cmp/b.hs abort, 2 charged instructions); the abort is a
  defined failure path, and a server reaching it is a Chernoff event below
  2^-400 under H' (failure term arena_overflow).
* Instruction counts are analytic listings under the fixed C = 2,000, not
  hardware timings. The reference-operation operand traffic is charged as in
  the earlier rho ledgers (18 transfer instructions per addition). The shared
  counters use ARMv8.1 ldadd, a real scalar AArch64 instruction, with no retry
  loop: the 128-bit counter is two 64-bit words and only the server that
  observes a low-word wrap updates the high word. The high word is checked
  against its cap at every carry, and the local block counter bounds the
  deterministic lag by 2^71 * 2^16 comparisons, which is charged.

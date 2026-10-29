# Fixed-step walks with a routed class table on BLS12-381 G1

**Claim.** time_log2 = 126.1835. For every challenge T = [x]P0 (hence for a
uniform one) the algorithm below outputs x with probability at least 0.39 over
its own coins, and its total work is at most W_max = 2^126.1835 units (one unit =
one batched incomplete affine G1 addition = C = 2,000 scalar AArch64
instructions). The computed total is W < 2^126.183349; the claim keeps a margin
of 0.00015 in log2. The work bound holds on every run: every data-dependent
action is checked against its allotment immediately before it runs and the walk
set-up samplers have explicit attempt caps (section 2, last three paragraphs). Peak memory < 2^132.3 bytes (section 6); no advice;
preprocessing (walk set-up, routing tables, table zeroing) < 2^117.4 units, all
inside W.

**What is different from the six-class rho line of submissions.** Those walks
must be functions of the class, so every step computes the canonical twist index
(about 50 of their 98 per-step instructions) and runs cycle handling; their
success law is the random-mapping heuristic H1 with a 0.5% allowance. Here the
walks are plain arithmetic progressions with a per-walk fixed step, so no twist
index, no cycles and no distinguished points are needed; instead every point's
class key is routed to a distributed linear-probing table that detects any two
points in the same class. The success probability is **proved** (section 4)
from exact linear-algebra probabilities of the collision events over the
algorithm's coins and the Bonferroni inequality; no random-mapping heuristic is
used. The one heuristic left (H', key bits behave as uniform) is used only for
the cost-side bounds (the expected linear-probing work and its concentration,
which the abort relies on, window counts, false candidates) and a negligible
detection corner; it never enters the probability that a collision exists.
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

**Final table (listing F).** Each final node keeps a linear-probing table of
S_node 32-byte slots (plus a guard tail of 2^12 slots instead of wrap-around),
zeroed once at start-up. The home slot is floor(h64 * S_node / 2^64) with
h64 = bits 35..98 of K2:K1 (`extr` then `umulh`). Insertion probes from the
home slot: at an empty slot (word K1 = 0) it stores the entry; at each
occupied slot it compares (K1, K2) and, if equal, raises a **candidate**; then
it moves to the next slot. The global load is alpha = n / S_tot = 0.45
(S_tot = 2^71 S_node).

**Candidates.** A candidate gives two locators (w, j), (w', j'). The
verification action recomputes both points from the directory, compares their
full keys (false candidates stop here), finds u with Q_{w,j} = [u]Q_{w',j'}
(section 1), and if b_w - u b_w' ≠ 0 sets

    x = (u (c_w' + j' D_w') - c_w - j D_w) / (b_w - u b_w')  (mod r),

checks [x]P0 = T and outputs x. Each candidate is priced at 2^14 units
(four scalar multiplications plus keys, unit search and the check: < 2^13.2).

**Rounds.** Generation is organised in rounds of 2^22 batches. In a round each
generator writes into double-buffered windows of capacity 2^11 entries per
intermediate (expected 2^10); each intermediate writes into windows of
capacity 2^15 per final node (expected 2^14); window write pointers are reset
and entry counts are published at round ends. A window can only overflow if a
Binomial count exceeds twice its mean (probability < 2^-400 overall under H').

**Random-cost allotment and abort, checked per action.** The only
data-dependent costs are (i) occupied-slot inspections in the tables,
(ii) candidate verifications and (iii) zero-denominator recovery actions. The
attack is **one sequential program** (section 6): it keeps one 64-bit
inspection counter x26 with one global cap x27 and one candidate counter x28
with one global cap x29. These are ordinary registers, not per-node counters;
the 2^71 final nodes are memory regions visited in sequence, so the same
counter covers all of them and no atomic operation is involved. The first
instructions of every occupied-slot inspection are
`add x26,x26,#1; cmp x26,x27; b.hs abort` (listing F, all three instructions
are counted in the 8-instruction per-slot price), so the program halts *before*
the inspection that would exceed the global allotment; every candidate is
likewise counted and compared with its cap *before* the 2^14-unit verification
runs. A zero-denominator batch is recovered by its 2^17-unit scan only while the
generator's recovery counter is below 4, checked before the scan; a generator
makes < 2^63.4 batches, each zero with probability < 2^24 * 3 / r < 2^-228.8,
so under H' some generator exceeds 4 recoveries with probability < 2^-120.
Hence every run does at most the deterministic work, the global inspection
allotment, the global candidate allotment, at most 2^39 x 4 x 2^17 = 2^58 units
of recovery work, and the one round (2^-41 of the total) of bookkeeping
overshoot already priced in section 5. This replaces the round-end aggregation:
no round can spend more random work than the global allotments. The global
allotments are exceeded with probability < 2^-111 (inspections, section 4) and
< 2^-100 (candidates, section 4), so they are failure events that do not
threaten the 0.39 success bound.

## 3. Correctness details

1. **Detection is complete.** Two points in one class have equal K, hence the
   same digit1, digit2 and home slot. Linear probing never deletes, so when the
   later of the two is inserted the earlier one lies in the run from the shared
   home slot to the first empty slot, and the probe compares (K1, K2) at every
   occupied slot of that run: the pair is a candidate. The only exception is a
   stored entry whose K1 word is 0, which reads as empty; its probability is
   2^-64 per entry and its effect on the collision events is < 2^-50 in total
   under H' (term "empty_marker").
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
empty-marker corner 2^-50, window and guard-tail overflow 2^-400 each,
global inspection-allotment overrun 2^-111, global candidate-allotment
overrun 2^-100, sampler cap 2^-371, recovery cap 2^-120; the last six
use H' or the explicit caps):

**Success >= 0.390001 - 2^-49 > 0.39 for every x.**

Collisions inside one walk (5 of the 6 units also give x there) and the
start points' collisions are extra and only help; they are ignored in the
bound. No property of the curve beyond Lemma 1 (group structure, the unit
group of order 6 and the complete invariant K) is used: the collision law is
exact algebra, not a random-mapping assumption.

**H' (heuristic, cost side only).** The digit1, digit2 and home-slot bits of K
of distinct classes behave as independent uniform values. It is used for: the
expected number of occupied-slot inspections of linear probing (Knuth: an
unsuccessful search at load a inspects (1/(1-a)^2 - 1)/2 occupied slots on
average, so the average over a filling table of final load alpha is
alpha/(2(1-alpha)) = 0.4091 at alpha = 0.45), its concentration, the window
counts (Chernoff), the false-candidate count, and the empty-marker corner. It
never enters the probability that a collision exists.

**Global inspection tail.** The cap is one global counter compared before every
inspection, so the relevant event is exactly the total occupied-slot
displacement D over all tables, not a per-node event. At alpha = 0.45,
E[D] = n alpha/(2(1-alpha)) = 0.40909 n = 2^124.8388. The total displacement of
linear probing at fixed load < 1 has variance asymptotic to c(alpha) S_tot with
c(alpha) finite (Flajolet-Poblete-Viola 1998; Janson 2001). We use the
conservative bound Var[D] <= 16 S_tot < 2^131.28, far above the published
constant. The cap is 1.10 E[D], so the slack is 0.10 E[D] = 2^121.51 and
Chebyshev gives

    P(D > 1.10 E[D]) <= 16 S_tot / (0.10 E[D])^2 < 2^131.28 / 2^243.02
                      < 2^-111.7.

A fatter tail constant only changes the exponent polynomially and leaves the
union with the other failure terms below 2^-49. The same 10% slack absorbs a
mild failure of exact uniformity.

**Global candidate tail.** Under H' the number X of false candidates among the
same-home comparisons satisfies E[X] < 2^124 * 2^-36 = 2^88 (section 5). For
any k, a union bound over k-subsets of the comparisons gives
P(X >= k) <= C(2^124, k) (2^-36)^k <= (e 2^88 / k)^k. At the cap k = 2^104
this is < 2^-100 (in fact < 2^-2^108), so the global candidate cap aborts with
negligible probability; the one true candidate is far below the cap.

## 5. Cost

All instruction counts are scalar AArch64 instructions; the reference
operation (one batched incomplete affine addition, 6M) is 1 unit and is never
re-counted. Listings are counted by `experiment.py` from the text below.

**Per point, generator (48).** Operand traffic of the reference addition,
charged as in the accepted rho ledgers: 6 `ldp` current point, 6 `ldp` step
point, 6 `stp` new point (18). Pass bookkeeping of the Montgomery-batched
addition, charged explicitly although the rho ledgers treat it as part of the
reference operation: the forward pass stores each prefix product (3 `stp`) and
the backward pass reloads it (3 `ldp`), each pass has `cmp`/`b.ne` loop
control (4; addresses advance by pre/post-indexed addressing), plus 2 spare
(12). The field multiplications and the A = 0 additions of the two passes are
the reference operation itself. Then listing G (18):

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
      cmp  x19, x20
      b.ne gen
```

A borrow only propagates upward, so the three `subs/sbcs` give limbs 0-2 of
p - y exactly; part A emulates this row limb by limb (random, boundary and
borrow-chain y) against K mod 2^192 and digit1. x22 (walk id) and x23 (j) are
per-batch registers.

**Per entry, intermediate (10).**

```
int:  ldp  x0, x1, [x9], #32         // K1, K2 (post-index)
      ldp  x2, x3, [x9, #-16]        // walk id, j
      and  x4, x0, #0x7ffffffff      // digit2 = bits 0..34 of K1 (valid bitmask immediate)
      ldr  x5, [x21, x4, lsl #3]
      stp  x0, x1, [x5]              // remote store into the final node's window
      stp  x2, x3, [x5, #16]
      add  x5, x5, #32
      str  x5, [x21, x4, lsl #3]
      cmp  x9, x10
      b.ne int
```

**Per entry, final node main path (11)** and **per inspected occupied slot (8)**:

```
fin:  ldp  x0, x1, [x9], #32
      ldp  x2, x3, [x9, #-16]
      extr x4, x1, x0, #35           // h64 = bits 35..98 of K2:K1
      umulh x4, x4, x23              // home = floor(h64 * S_node / 2^64)
      add  x5, x22, x4, lsl #5       // slot address (32-byte slots)
      ldp  x6, x7, [x5]
      cbnz x6, occ
empty:stp  x0, x1, [x5]
      stp  x2, x3, [x5, #16]
      cmp  x9, x10
      b.ne fin
occ:  add  x26, x26, #1              // inspection counter (allotment)
      cmp  x26, x27                  // immediate allotment check, before the inspection
      b.hs abort
      cmp  x6, x0
      ccmp x7, x1, #0, eq
      b.eq cand                      // separately priced verification
      ldp  x6, x7, [x5, #32]!
      cbnz x6, occ
      b    empty                     // once per occupied action
```

**Per point total (deterministic):** 48 + 10 + 11 = 69 instructions, plus the
per-batch work (Fermat inversion 380 squarings + 228 multiplications at 383 =
232,864, exponent scan < 2,000, zero-denominator test on the batch product,
per-batch constants: < 235,000 per batch of 2^24 points, i.e. < 0.0141), plus
round bookkeeping < 0.02 (per round: a generator resets 2^36 window pointers
and publishes 2^36 counts, < 2^38.2 instructions against 2^46 points; an
intermediate reads 2^39 window counts and loops over them, resets and
publishes 2^35, < 2^42.1 against 2^49 entries; a final node reads one window
per round of ~2^14 entries; together < 0.016). Point price:
1 + 69.034007/2000 = 1.034517004 units. (Start points are also charged this
price, on top of their set-up, which over-counts.)

**Table zeroing.** S_tot = n / 0.45 slots of 32 bytes, zeroed with `stp xzr,
xzr` (2 per slot) in a loop of 64 stores per 3 loop instructions:
2.09375 instructions per slot, i.e. 4.653 per point, < 2^117.39 units in all.

**Inspections (random, allotted).** Expected 0.4091 occupied slots per insert.
Each probe step is charged 9 instructions: the 8 listed above (increment,
immediate cap check and branch, key comparison and branch, slot advance) plus
one for the `b empty` branch, which fires at most once per probe step (one per
occupied action). Thus the expectation is 9 * 0.4091 = 3.682 instructions per
point and the global allotment is 1.10 times it = 2^117.1805 units. The counters
and caps are 64-bit registers of the one sequential program and the cap is the
global allotment, so the immediate checks cap the total inspection work on
every run; exceeding the cap is a failure event with probability < 2^-111.7
(section 4).

**Other terms.** Walk set-up (capped sampling and three scalar multiplications
per walk, plus directory writes) < 2^14 units x 2^63 = 2^77; routing tables and
window mappings < 2^76; candidates (expected false candidates: same-home
comparisons < 2^124 times 2^-36 free key bits, i.e. < 2^88, plus the true one)
x 2^14 units, allotted 2^104 and enforced before every verification by the
counter check; capped zero-denominator recovery 2^39 x 4 x 2^17 = 2^58; one round
of overshoot < 2^85.2; final output and verification < 2^20.

**Total.**
W <= n * 1.034518 + 2^117.39 + 2^117.19 + 2^104 + 2^85.2 + 2^58 + 2^77 + 2^76 + 2^20
< 2^126.1834 <= 2^126.1835 = W_max.

For orientation: the accepted rho ledger prices an ordinary step at
1 + 72.0561/2000 = 1.03603 units and needs S_ops = 1.005 S0 operations; here a
point costs 1.03452 units plus 0.0038
units of zeroing and probing, and n is the rigorous Bonferroni count (0.06%
above the Poisson value instead of a 0.5% allowance). The immediate global
counter checks add 3 instructions per occupied inspection (0.0006 units per
point on average) and no per-point instructions; they are what makes the
every-run bound independent of any round-level expectation. The exact ledger
computed by `experiment.py` is W < 2^126.183349.

## 6. Memory and machine model

Peak storage < 2^132.3 bytes: tables (S_tot + guards) x 32 B = 2^132.28,
receive windows (2^36 intermediates x 2 x 2^39 x 2^11 x 32 B = 2^92, final
nodes 2^71 x 2^21 B = 2^92), walk records (504-byte stride, so every pre-indexed `ldp` offset, including
`#REC` = 504, is encodable) and directory (2^63 x 600 B = 2^72.3). Every
node's own memory is below 2^62 bytes and all its addresses (local data and the
remote windows it writes: at most 2^36 x 2^17 B or 2^35 x 2^21 B, double-buffered) fit in a
64-bit address space; no instruction uses an address wider than 64 bits.

**One sequential program.** The attack is executed by one sequential processor:
the 2^39 generators, 2^36 intermediates and 2^71 final tables are memory
regions visited in sequence, each with a 64-bit base pointer (its own region is
below 2^62 bytes), and the "remote stores" are ordinary stores to the next
region's receive window at a 64-bit address. The global inspection and
candidate counters x26/x28 and their caps x27/x29 are therefore ordinary
registers of this one program, updated by three or four instructions per action
(listing F and the candidate check); no atomic, network or per-node counter is
involved. The work is the sum over the sequential schedule, exactly the same
instructions as the distributed description, so the score is unchanged.

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
< 2^126.1882, which is above the current rho frontier; the claim does not use
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
  routing digits and a linear-probing table built from the emulated rows, the
  candidate, re-computation from the directory, the unit search, the solving
  formula, and [x]P0 = T.
* **Part B (toy curves y^2 = x^3 + 4 over F_q with a prime-order subgroup and
  the same six units):** the whole attack at toy size and the exact scalar
  model at small r. Long offline runs (same code, `toy_seconds` = 300,
  `scalar_seconds` = 300) are tabulated in README.md: toy curves q26, q30, q34
  give 0.4165 +- 0.0027, 0.3940 +- 0.0035, 0.4066 +- 0.0058 against
  predictions 0.415-0.422, 0.391-0.397, 0.394-0.398 (cross-walk-only to
  all-events mu), and the scalar model gives 0.4087 +- 0.0013 (145,162 trials,
  prediction 0.4088) and 0.4103 +- 0.0011 (192,048 trials, prediction 0.4104). The toy runs are evidence for the
  implementation, not for the bound, which is proved in section 4.
* **Part C:** the ledger above in exact fractions: L from f4, the price, the
  zeroing and probing terms, the extras, W, the dependency-bound exponents, and
  checks that every listing's counted length equals its ledger entry.

## 8. Limitations

* Memory is enormous (2^132.3 bytes) and is not scored by the rules.
* The machine model treats a store into a mapped remote window as one store
  instruction (section 6 gives the bound with an extra copy per hop).
* H' (uniform key bits) is a heuristic; the abort, overflow and false-candidate
  bounds rely on it (each with a stated negligible probability, and a 10%
  slack on the global probing allotment). The inspection tail uses the
  published linear-variance asymptotic of linear probing (Flajolet-Poblete-Viola
  1998; Janson 2001) with the very conservative constant 16; the collision law
  itself is proved and does not use H'.
* Instruction counts are analytic listings under the fixed C = 2,000, not
  hardware timings. The reference-operation operand traffic is charged as in
  the earlier accepted rho ledgers (18 transfer instructions per addition).
  The single sequential program's global counters are ordinary registers;
  their three-instruction checks are priced inside the 8-instruction probe
  price.

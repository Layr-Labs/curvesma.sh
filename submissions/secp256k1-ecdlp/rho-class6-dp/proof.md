# secp256k1-ecdlp: parallel Pollard rho on 6-element classes (negation and the GLV endomorphism)

Claim: recover x from Q=[x]P with probability at least 0.39 using at most 2^127.47
target-group-operation equivalents (one reference operation = one batched incomplete affine
addition = 1000 instructions, C=1000). Classical, randomized, heuristic only in the walk model below.
Peak memory is at most 2^76.6 bytes; preprocessing is at most 2^32.01 reference operations
(included in the total); nonuniform advice is 2^0 bytes. Nothing here is a practical attack:
the memory and work are far beyond any machine.

## 1. Classes

n is the prime group order (about 2^256, cofactor 1). Let beta (cube root of unity in F_p) and
lambda (cube root of unity in Z_n) be the public secp256k1 endomorphism constants, with
lambda*(x,y)=(beta*x,y) (checked in the setup, see section 7). The group H={+-1,+-lambda,+-lambda^2}
has order 6 and acts freely on nonzero scalars (lambda-1 and 2 are invertible mod n), so
n-1 is divisible by 6 and the nonzero points fall into N'=(n-1)/6 classes, log2 N' = 253.41504.
Canonical representative of a class: among the three x-values x, beta*x, beta^2*x take the
smallest as an integer (they are distinct because x is nonzero), keep y for the matching point,
and if y is odd replace it by p-y. This picks exactly one of the six points, and the canonical x
alone identifies the class.

## 2. Algorithm

Setup: r=2^20 table entries T_j=[s_j]P+[t_j]Q with fresh uniform 256-bit s_j,t_j. Distinguished
point (DP): canonical x has its low 56 bits zero, theta=2^-56. The walk function on classes is
f(c)=canon(c+T_h(c)) with h(c)=bits 56..75 of the canonical x of c (independent of the DP bits).
The main run does not track (a,b): a step only adds the table point and canonicalises, and a walk is
identified by its seed. On a DP match the two walks are re-run from their seeds with the same rule, now tracking
(a,b) with the point equal to [a]P+[b]Q (a table step adds (s_j,t_j), canonicalisation multiplies both by the
scalar of H that was applied, a doubling escape doubles and rescales them); this re-run is charged in section 6 and
is not part of E.
Each walk starts at canon([a]P+[b]Q) for a seed-derived (a,b), and runs in a batch of m=2^10
interleaved walks, so one Montgomery batch inversion serves the batch (this is the batched
addition of the reference operation). A walk stops at its first DP, or after 40/theta steps (cut,
restart).

Fixed-rule cycle handling. A fruitless cycle of f (a closed orbit of classes) is left by a rule that
depends only on the cycle, never on the walker. Each walker keeps a ring of the canonical x values of
its last 8 points, cleared (and refilled with the single current point) at its start and after each escape.
After every step, in this order:

1. If the new point is a DP: stop and store it.
2. If the new point's canonical x equals a ring entry d steps back (d=1..8; the ring holds full 256-bit canonical x values and a match means all 256 bits are equal), the walker has closed
   a cycle of f of length L=d (see section 3). The ring entries of the last L steps are the cycle
   members. Take m, the member with the smallest canonical x, and walk the cycle with f from the
   current point to m (at most L-1 steps, no further tests: the walker has already passed every
   member without meeting a DP).
3. Escape from m by doubling: e(m)=canon([2]m), a function of the class of m only (doubling commutes with the
   six symmetries: [2](u m)=u [2]m for u in H, so the canonical x of [2]m depends only on the class). The point
   [2]m is not a table step: it is related to m by the multiplier 2, which is not in H (2^6 is not 1 mod n), so it has
   no algebraic link to the table entries that closed the cycle. The pair (a,b) tracked in a re-run is doubled (2a,2b) and
   rescaled by the unit that canonicalisation applied. Test e(m) for a DP, clear the ring to e(m), continue.
   The doubling is done in the batch as that walker's step: the affine doubling formula (slope 3x^2/(2y), y nonzero
   because n is odd) shares the batch inversion with the additions of the other walkers.
4. Otherwise push the new point into the ring (dropping the oldest of 8).

A cycle longer than 8 is never recognised; the walker follows it until the 40/theta cap and is charged for it (section 6).
Each DP is stored as (canonical x, 16-byte seed id, 8-byte length) in a sorted, searchable structure,
and looked up first. On a match with different walks, both walks are re-run from their seeds with
the same rule, tracking (a,b); the two canonical points are equal, so a1+b1*x = a2+b2*x and
x=(a1-a2)/(b2-b1) mod n if b1 != b2. The result is verified with [x]P=Q. The run is cut at
K_max = nu*k0 + m*40/theta steps in total.

## 3. Colliding walks stay merged

Claim: if two walkers (with any histories, started at different times) are ever at the same class c, and
neither is cut by the 40/theta cap or caught in an unrecognised cycle (length at least 9), they stop at the same
distinguished point.

Proof. (i) Between two resets every move of a walker is an f move, so a ring entry equal to the new point
means the f-orbit returned to a class: f^d(c)=c, and c lies on an f-cycle. The points of one f-orbit are
distinct until the first return, and the first return to a class is exactly one cycle length after the orbit first entered the cycle;
so the first match has d=L, the length of that cycle, and the last L ring entries are exactly the members of the
cycle. Hence the member set, its smallest canonical x (the member m) and the escape point e(m) are functions of the cycle
alone, not of the walker, of when its ring was cleared or of where on the cycle it was caught.
(ii) A walker whose orbit enters a cycle of length L <= 8 detects it: it returns to the entry class L <= 8 steps later,
which is still in its ring, unless it meets a DP first (then it stops at a DP of the orbit, in agreement with every walker
on that orbit). Two walkers on the same orbit may be caught at different cycle points (one older
history, one fresh); an older ring entry that matches the current point also lies on the cycle (the walker went
p -> ... -> c -> ... -> p along f moves, so c is on that loop), so both identify the same cycle and the same m.
(iii) Both therefore leave the cycle at e(m)=canon([2]m) with the ring cleared: the same class and the same state. By induction on
the (finite, cap-bounded) number of cycles on the orbit, the sequence of classes after c is a function of c alone up
to the order in which cycle members were passed, and it ends at the same DP. The coefficient pairs (a,b) differ
between the walkers, which is exactly what gives the relation in section 2. The toy experiment checks this property directly
(experiment.py `merge_test`: a second walker restarted at a random point of a first walker's trail, with empty
history, always ends at the first walker's DP; the old rule that escaped by a history-dependent test fails it). QED

Re-entry. The rule is deterministic, so if e(m) lay in the basin of the same cycle, every walker that reaches that
cycle would be sent through m -> e(m) -> ... -> cycle -> m again, and stay in this loop until the 40/theta cap
(a trap). Whether it does is a question about the point [2]m, not about merging. For an escape through a table step
the question has a bad answer, for a reason that is algebraic: e(m)=canon(m+T_j) lies one table step from m. With the
cycle's own index j=h(m) it is f(m) and always returns to the cycle; with a different index it re-enters at about the basin fraction, but
when m sits on a 2-cycle it still returns to it with probability about 1/(6r)=2^-22.6 per escape (the same index and the
negation), so that rule is trapped by about 1/(6r) of its escapes, a loop that repeats until the cap. Doubling has no such
link. A walker started at e(m) enters the cycle only if e(m) lies in the basin of the cycle (the classes whose f-orbit
reaches it) at distance k from the cycle with k at most T=40/theta=2^61.32, since a walker stops at its first DP, and
is cut at T steps. In the random-mapping model each class has Poisson-like in-degree of mean 1 (6r candidate
points, each canonical with probability 1/6 and selected with probability 1/r), so the expected number of classes
at each distance k from a given cycle member is 1 and the basin within distance T has at most L(T+1) <= 8*2^61.32 = 2^64.32 classes
in expectation. The class of [2]m is a uniform-like point relative to this set, so
Pr[e(m) re-enters its cycle] <= 2^64.32/N' = 2^-189.1 per escape. A relation that makes it happen would be
(2-h)m = a sum of table terms (h in H), whose multiplier on m is nonzero because 2 is not in H, so it is not an
identity that any table-index coincidence implies. The escape is deterministic per cycle: all escapes from one cycle share the single point e(m), so the union
bound runs over distinct cycles, not over escapes or walkers. Distinct cycles are at most the number of escape events,
at most 2^127.47*2^-20 = 2^107.5, so the
expected number of trapped cycles is at most 2^-81.6 and each wastes at most 2^61.32 steps, so the
expected waste is far below 2^-20 of the total: it is part of epsilon and negligible. This bound is a model bound, not a theorem (the random mapping model and the
independence of the point [2]m from the cycle structure); see section 8.

Residual events (all charged or negligible): the cap cut (accounted in section 4), a cycle of length at least 9
(section 6), e(m) falling back into the basin of the same cycle (probability at most 2^-189.1 per event by the re-entry paragraph above; then the walk is cut at the cap),
an equal-x addition (2^-250 per step).

## 4. Success probability (at least 0.39)

Let k0 = sqrt(-2 ln(0.61) N') = sqrt(0.98859 N'), log2 k0 = 126.69924, and nu=1.03.
Walk points are modelled as independent uniform draws from the N' classes (van Oorschot-Wiener
1999 parallel collision search; Teske 2001 r-adding walks, r>=16 behave like random walks to
within a few percent, here r=2^20). Among k points the probability of a collision is
1-exp(-k^2 rho/(2N')) where rho<=1 is the walk-quality factor. We take nu*k0 points: for rho=1 this
is 1-exp(-0.4943*1.0609)=0.4080. Pairs inside one walk are not useful; they are a fraction at most
2^56/K_max, below 2^-70, of all pairs. The first collision between two different walks gives a
relation with b1 != b2 except with probability about 2^-250 (the starts are uniform). By section 3 the two
walks reach a common DP, each within 40/theta further steps with probability at least
1-2e^-40; during that time at most m walks advance, so at most m*40/theta=2^71.3 more steps
are charged. The remaining losses are walks lost in unrecognised cycles (fraction at most
2^-88.8 per step times 2^56 steps per walk, 2^-33) and the events listed after section 3; we budget epsilon=2^-20 for all of them.
Result: P >= (1-exp(-0.4943*nu^2*rho))(1-2^-20), which is 0.3987 at rho=0.97, 0.3955 at rho=0.96
and 0.3924 at rho=0.95. It stays above 0.39 down to rho=0.943: the margin nu=1.03 absorbs a walk quality of rho>=0.943.
No term depends on the walker's history: merged walks cannot diverge (section 3), so no divergence allowance is needed.
The heuristic is that the 2^20-entry random table walk collides at least 94.3% as fast as an
ideal random mapping on N' classes; Teske measured deviations of a few percent already for r=16-20
and they fall as r grows. The challenge Q is uniform, so the start randomness and table are independent of x.

## 5. Per-step cost

One walk step costs 1 reference operation (the batched incomplete affine addition) plus E
instructions of other work, divided by C=1000. The addition formula is incomplete: an equal-x
case (T_h=+-c) has probability about 2^-250 per step.

| Item | AArch64 instructions |
| --- | ---: |
| beta*x and beta^2*x: 2 field multiplications at 179 | 358 |
| normalise x, beta*x, beta^2*x, y to unique form: 4 x 40 | 160 |
| two 256-bit comparisons with selects (choose minimal x) | 36 |
| y parity test and conditional negation p-y | 20 |
| table index, address, load of T_h (10 limbs) | 16 |
| DP test (mask and branch) | 3 |
| cycle test: 8 compares at 2 (early exit on the first differing limb; a false 256-bit match has probability 2^-256), ring-buffer store and index | 20 |
| walk state load/store, length counter, batch loop | 28 |
| Total | 641 |

We use E=650: the itemised 641 plus a reserve of 9, a bound on the instructions of one step, not a measurement. Step cost = 1+0.650 = 1.650 reference operations. The 6M addition
itself already includes the batch amortisation, so no inversion is charged separately.

## 6. Cycle events and other costs

Which cycles exist. A closed orbit c_0 -> ... -> c_L=c_0 with table indices j_1..j_L and unit scalings
u_i in H gives a relation (U-1)c_0 + sum_i V_i T_(j_i) = 0 with units V_i, U in H. The T_j are independent uniform
points, so if some table index occurs exactly once among j_1..j_L the relation fixes a uniform point and holds
with probability at most 2^-250. A cycle therefore needs every index that occurs to occur at least twice.
Modelling the indices of L successive steps as independent uniform draws from r=2^20 values, the probability
of that event (dropping the further requirement that the units cancel, a factor at most 1) is
S_L = sum over partitions of L positions into blocks of size at least 2, b blocks, of r^(b-L).
Exact values: S_2=2^-20, S_3=2^-40, S_4=2^-38.42, S_5=2^-56.68, S_6=2^-56.09, S_7=S_8=2^-73.29. For L=2 the true
rate is 1/(6r) (the new class representative must be -(c+T_j), probability 1/6, and the same index must repeat);
we charge S_2=1/r.

Detected cycles (L<=8). Events per step at most S_L. One event costs L steps to close the cycle, at most L-1
steps to reach m, one doubling escape and the selection of the minimum (at most 8 comparisons, under 0.15 reference
operations, covered by charging one more step). The doubling escape costs a batched affine doubling (charged 1.5
reference operations: the slope 3x^2/(2y) needs 2 squarings and a multiplication more than the shared-inversion addition, and the
batch inversion is the same) plus canonicalisation and the DP test (E=650, section 5 rows 1 to 6): 2.15 reference
operations = 1.30 steps; we charge 2 steps for it. One event is therefore at most 2L+2 steps. Relative waste
w_c = sum_(L=2..8) S_L (2L+2) = 5.722e-6 (of which L=2 gives 5.722e-6 = 6*2^-20; L=3..8 together under 1e-10).

Unrecognised cycles (L>=9). Exact sum for L=9..40: S_9=2^-89.70, S_10=2^-90.12, S_11=2^-105.9, S_12=2^-106.7, ... total
2^-88.89. For 41<=L<=2^17 the repeat count bound C(L,ceil(L/2))(L/r)^ceil(L/2) is at most 2^-267 per L; for larger L every
one of at least r/16 distinct indices must have cancelling units (probability at most 1/6 each), under 2^-100000. So the
rate is e9 <= 2^-88.8 per step. A walker in such a cycle is cut at the cap, wasting at most 40/theta = 2^61.32 steps, so the relative
waste is w_u = e9 * 2^61.32 = 5.0e-9.

Structureless cycles (relation not cancelling index by index) need the walk to hit a specific point: at most 2^61.3 cyclic points of that length per random mapping against N' classes, per step, times 2^127 steps gives 2^61.3*2^127/N' = 2^-65 in total, negligible.

| Component | Bound |
| --- | ---: |
| Walk steps: K_max = nu*k0 + m*40/theta = 2^126.74189 + 2^71.32, times 1.650 | 2^127.46435 |
| Detected cycles (with the doubling escape charged), w_c = 5.72e-6 of the steps | 2^110.05 |
| Unrecognised cycles, w_u = 5.0e-9 of the steps | 2^99.89 |
| DP store and lookup: D=K_max*theta = 2^70.74 DPs, 2^17 instructions each (about 71 comparisons of 56-byte records plus moves and transfers, with allowance) | 2^77.78 |
| Walk starts: 2^70.74 starts, 4096 reference operations each (random [a]P+[b]Q by Shamir double and add in Jacobian coordinates at 16M per operation, one inversion; allowance 4096) | 2^82.74 |
| Re-running two walks with coefficient tracking: 2 * 40/theta * 8 ops | 2^65.32 |
| Table: 2^20 entries * 4096 | 2^32 |
| Setup, constants check, scalar recovery, [x]P=Q verification, output | 2^14 |
| Total (walk steps times 1+w_c+w_u, plus the other rows) | W = 2^127.46436 |

log2 W = 127.46436 and we submit X = 127.47, rounded up. The rows other than the walk steps add under 3e-5 to log2 W in total.
Choosing theta=2^-56 instead of 2^-32 cuts the store by 2^24 and changes log2 W by under 1e-5: the detection delay 2^71.3 and the
cap-cut waste are still far below the walk steps.

## 7. Preprocessing, memory, advice

Preprocessing (table and setup) is 2^32.0 reference operations (2^32.01 with setup), charged in
the total. Setup checks beta^3=1, lambda^3=1 mod n and that lambda*P=(beta*P_x,P_y), by one scalar
multiplication. Memory: the DP records, 56 bytes each (32-byte x, 16-byte seed id, 8-byte length, the length
at most 40/theta=2^61.3) times D=2^70.74 is 2^76.55 bytes; the table (2^27 bytes), 2^10 walk states and the 1 MiB code
reserve are smaller; we state at most 2^76.6 bytes. The memory is astronomically large; the
rules give memory no score. Nonuniform advice is zero (table entries and constants are
computed in the run, public constants excluded).

## 8. Heuristics and limitations

1. Walk randomness (section 4): the 2^20-entry table walk on classes behaves like random draws to
   within 5.7 percent in collision rate. Support: Teske 2001, van Oorschot-Wiener 1999.
2. The cycle counting of section 6 models the table indices of successive steps as independent uniform
   draws and the relation test as a union over index patterns, in the spirit of Bos, Kaihara, Kleinjung, Lenstra and
   Montgomery 2010 and Bernstein, Lange, Schwabe 2011 for classes of size 2; it is not a theorem. The cost margin is
   large: even a 500-fold error in w_c (including the doubling's cost) would leave log2 W below 127.47.
3. The escape point canon([2]m) and the minimum member are functions of the cycle, so the merging property of section 3 does not depend on
   the cycle rates; only the cost does. The trap bound (escape re-enters its own cycle, at most 2^-189.1 per distinct cycle, since every escape from one cycle is the same e(m)) is a
   random-mapping model: in-degree of mean 1 and the independence of [2]m from the cycle's basin. It is not a theorem; the table-step escape
   it replaces failed this very question by an algebraic link (about 1/(6r) re-entries), and doubling removes that link. The doubling's price
   (1.5 reference operations plus E) is a stated bound, not a measurement; it enters only at rate 2^-20 per step.
4. The cost of 2^17 instructions per stored DP and the structure of a 2^76-byte sorted store is an
   accounting allowance, not an engineering claim; it is 2^-50 of the total.
5. Canonicalisation prices (40 per normalisation, 36 for compare and select) are stated bounds
   on scalar code with 5x52 limbs, not measured counts. E=650 is stated as a strict upper bound (641 itemised plus 9
   reserve); the claim 127.47 holds for any true E up to 656 (log2 W = 127.4696 there), and each further 10 instructions add 0.009 to log2 W.
6. A small toy experiment (experiment.py, curve of prime order 523969, 87328 classes) runs the method end to end with the
   doubling escape, recovers every scalar, and checks the merging property of section 3 on 1600 restarted walkers (the toy table has
   16 entries so that cycle events are frequent: 288 of the first walkers saw one; all 1600 pairs end at the same DP). It also counts re-entry: on 20000 fresh
   walks with 3057 cycle events (a walk-level count, not independent draws: walks reaching one cycle share its e(m), so distinct cycles are fewer) the doubling escape re-entered its cycle 0 times, while the earlier table-step escape (control,
   same cycle rule) re-entered 5345 times (loops repeating m -> d -> m until the cap). It shows the algorithm, the confluence and
   the re-entry difference, not production cost or walk quality.

References: van Oorschot and Wiener, Parallel collision search with cryptanalytic applications,
J. Cryptology 1999; Teske, On random walks for Pollard's rho method, Math. Comp. 2001; Bos,
Kaihara, Kleinjung, Lenstra, Montgomery, On the security of 1024-bit RSA and 160-bit elliptic
curve cryptography (2010); Bernstein, Lange, Schwabe, On the correct use of the negation map in
Pollard rho (2011); Gallant, Lambert, Vanstone, Faster point multiplication using efficiently
computable homomorphisms (2001).

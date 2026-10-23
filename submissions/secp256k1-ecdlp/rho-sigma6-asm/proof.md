# secp256k1-ecdlp: parallel Pollard rho on order-6 classes with a sigma walk, all arithmetic in instructions

Claim: for uniform x in [0,r) and Q=[x]P on secp256k1, the algorithm below returns x
with probability at least 0.39. It uses at most **2^126.418986** target-group-operation
equivalents in total (score 63209493/500000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^75.74996 bytes. Preprocessing is
2^47.61967 units, and no nonuniform advice is used. This is a classical cold attack on
one target.

The success bound rests on the disclosed heuristic H1, the standard random-mapping model
for r-adding walks on automorphism classes, and on the rate bound H2. The cost rests on
H3, instruction counts measured with the organizer's toolchain. This is a carefully
accounted generic attack. It is not a practical break.

**Cost route.** The walk does not use the 6M reference operation. Every instruction of
every step, including all field multiplications, squarings and subtractions, is its own
4×64-bit GCC inline AArch64 assembly and is charged at 1/C with C = 1000 (H = 0 in
COST_MODEL's H + R/C). COST_MODEL.md: "Faster arithmetic and lookup algorithms are welcome;
justify their actual work". A step measures 780 instructions in the organizer's
toolchain, so it is charged 0.795 units rather than 1 + R/C.

Relative to the usual order-6 walk, the walk never moves a point to its canonical
representative, so it never pays the multiplication by beta (section 1.2). It stores the
actual point W and adds sigma(R_i), where sigma is the automorphism that maps the
canonical point to W. The class sequence is unchanged, so the standard analysis applies.
Which rotation is canonical is decided from 64-bit fixed-point fractions with an exact
fallback, and Lemma 1 bounds their error. Compared with our earlier rho-order6 package (not included),
it removes the beta multiplication, and the budget is tau = 0.5303. The same walk with the
6M reference operation and libsecp256k1 field code costs 1.32 units per step (an earlier
version, not included).

**Distinguished points and resolution.** Walks run in trails. Each trail starts from fresh
coins and ends at its first distinguished point, which it records once together with its
coins. A key match is resolved by replaying the two trails (section 1.3, Lemmas R and D).
A fruitless cycle therefore cannot store the same point twice, and a degenerate match has
probability exactly 1/r.

Files:
- `accounting.py` recomputes every number with exact integers (output in
  `evidence/accounting.json`).
- `measurement/` holds the per-step instruction-count harness, the field code
  (`secp_field.h`, `secp_sub.h`), the Python replica `walk_ref.py`, and the approximation
  check `frac_check.py`.
- `evidence/measurement/` holds the counts on two input sets and the checks of the
  restart and abandonment paths.
- `experiment.py` is the declared toy-curve experiment, running the whole procedure;
  `evidence/trails/` holds larger offline runs of the same file.

## 1. Algorithm

### 1.1 Notation and parameters

p, r and P come from the protected profile. The curve is y^2 = x^3 + 7, the cofactor is 1,
and p ≡ r ≡ 1 (mod 3).

The endomorphism uses beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
and lambda = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72, with
[lambda](x,y) = (beta·x, y). The automorphism group G6 = {±1, ±lambda, ±lambda^2} has
(r−1)/6 nonzero classes, and N = (r−1)/6 + 1. Field values are compared as reduced
integers in [0,p).

Parameters:
- A = 2^16 adders.
- K = 2^32 walks: 2^16 processors × 2^16-walk batches.
- theta = 2^−60.
- L = 33310394024866356565296666183 ≈ 2^94.7500.
- CAP = 64/theta. Trails have fewer than TRAIL_MAX = CAP + 64 steps.
- DP_MAX = 2·theta·S_max + 2K, where S_max = K·(L + TRAIL_MAX).

### 1.2 Classes, keys and the sigma walk

For W = (x, y) ≠ O:
- **Class key.** key(W) = y if y is even, else p − y, so key(W) is ±y and even.
  - key is constant on the class, since G6 changes y only by sign.
  - key determines the class: y^2 = x^3 + 7 fixes x up to a factor beta^j.
- **Canonical point.** canon(W) = (beta^k·x, (−1)^s·y), where k minimises beta^k·x and
  s = y mod 2.
- **sigma.** sigma_W = (x → beta^(−k)·x, y → (−1)^s·y) = [(−1)^s·lambda^(−k)], so that
  W = sigma_W(canon(W)).
- **Walk step.** i(W) = (key >> 1) mod 2^16, and

      f(W) = W + sigma_W(R_i) = sigma_W(canon(W) + R_i).

  So the class sequence is exactly that of the standard class walk
  [W] → [canon(W) + R_h([W])]. The six images sigma(R_i) are tabulated per adder.
- **Distinguished point.** W is distinguished when bits 17..76 of key(W) are zero.
- **Representation.** In the kernel, x and y are fully reduced integers in [0, p), held as
  four 64-bit limbs.

**Deciding k without beta·x (Lemma 1).** Write x = Σ_{i<4} a_i·2^(64i), and let
C1 = ⌊2^384·beta/p⌋ = Σ_{j<6} C1_j·2^(64j). Define

    g_0 = a_3,
    g_1 = ( Σ_{i+j=5} lo(a_i·C1_j) + Σ_{i+j=4, j≥1} hi(a_i·C1_j) ) mod 2^64,

so g_1 takes 8 multiplications (4 `mul`, 4 `umulh`), and g_0 none. Let
f_j = 2^64·(beta^j·x mod p)/p = 2^64·frac(x·beta^j/p).

**Lemma 1.** f_0 ∈ [g_0, g_0 + 2). If g_1 ≤ 2^64 − 10, then f_1 ∈ [g_1, g_1 + 9).

*Proof.* For j = 0: x/2^192 ≥ a_3, and since p < 2^256,
f_0 − a_3 < 1 + 2^64·x·(2^256 − p)/(p·2^256) < 2.
For j = 1: beta/p = (C1 + e)·2^−384 with e ∈ [0,1), so
f_1 = (x·C1·2^−320 + x·e·2^−320) mod 2^64, and the second term is below 2^−64. Next,
x·C1·2^−320 = Σ_{i,j} a_i·C1_j·2^(64(i+j−5)). Modulo 2^64:
- terms with i + j ≥ 6 vanish;
- i + j = 5 gives lo(a_i·C1_j);
- i + j = 4 gives hi(a_i·C1_j) plus a fraction below 1 (four terms with j ≥ 1);
- the term with i + j = 4 and j = 0 does not exist (i ≤ 3);
- i + j = 3 gives four terms below 1 each;
- i + j ≤ 2 gives terms below 2^−64 each.

Dropping all of these gives an error in [0, 8 + 2^−50). The hypothesis excludes
wrap-around. ∎

x + beta·x + beta^2·x ≡ 0 mod p, with the reduced values in (0, p), so they sum to p or 2p
and f_2 ≡ −f_0 − f_1 mod 2^64. With g_2 = (−g_0 − g_1 − 22) mod 2^64, the real f_2 lies in
(g_2 + 11, g_2 + 22] whenever g_2 ≤ 2^64 − 23. So every f_j lies in [g_j, g_j + 22]. The
kernel accepts the least g_j as k only if all three g_j ≤ 2^64 − 23 and the least is more
than 22 below both others; then k is exact. Otherwise it computes beta·x and compares
exactly. The fallback rate is below 2^−56 per step and was never seen in 2^19 measured
steps. `measurement/frac_check.py` checks Lemma 1 and the decision rule on 200,000 random
and edge-case inputs per seed; the largest floor gap seen for g_1 is 6.

### 1.3 Procedure

1. **Validate Q.** Check it is on the curve. If Q = O, output 0.
2. **Adders.** For each adder i < 2^16:
   - draw c_i, d_i uniform in [0,r);
   - set R_i = [c_i]P + [d_i]Q;
   - tabulate its six images.
3. **Trails.** Each of the K walk slots runs a sequence of trails.
   - A trail starts from fresh coins a, b uniform in [0,r), with W = [a]P + [b]Q. The
     slot keeps the coins until the trail ends. Coefficients are not tracked while
     walking.
   - A trail ends at its first distinguished point (step 6), or is abandoned: when it
     crosses more than 2^60 anchor boundaries (step 5) without one, that is after between
     CAP − 64 and CAP = 64/theta steps, or when an addition gives O. Either way the slot
     starts a new trail from fresh coins. Every trail has fewer than TRAIL_MAX = CAP + 64
     steps.
   - Each slot takes L steps (the main phase), then finishes its current trail and stops.
     The run fails if more than K trails are abandoned.
4. **Step.** W ← f(W), as one batched incomplete affine addition: Montgomery's
   simultaneous inversion over the 2^16 walks of a batch (3M), then lambda, lambda^2 and
   y3 (2M + 1S), with one Fermat inversion per batch.
   - If the batch product is zero, the affected walk uses the unbatched complete formula.
5. **Fruitless cycles.**
   - *2-cycles.* Each walk keeps the low 64 bits of the key of the point before its
     current point W. If the new point's key matches, the walk evaluates f(new) unbatched.
     If key(f(new)) = key(W) exactly, then {W, new} is a 2-cycle of the class walk, and
     the walk replaces new by 2·canon(M), where M is the element with the smaller key.
   - *3..64-cycles.* At every batch step divisible by 64, every walk of the batch records
     the low 64 bits of its key as its anchor. All walks of a batch step together, so each
     walk is anchored once every 64 of its own steps. On a later match, the walk iterates
     f unbatched up to 64 times. If its class returns, it escapes from 2·canon(M), with M
     the cycle element of least key.
   - So a walk in a cycle of length ℓ ≤ 64 leaves it within 64 + ℓ steps, unless the
     cycle contains a distinguished point, which ends the trail first.
6. **Distinguished points.** W is distinguished when bits 17..76 of key(W) are zero. The
   test is applied to every point that a step or an escape produces. At a distinguished
   point the trail ends: the record (key(W), a, b), with the trail's coins, goes to a
   hash table keyed by key(W), and the slot starts a new trail. A trail makes at most one
   record, so no trail can store a point twice. The run fails if more than DP_MAX records
   arrive.
7. **Resolution.** The first record whose key is already in the table ends the run. No
   other key match is used.
   - Both trails are replayed, unbatched, from their recorded coins up to their first
     distinguished point, with the rules of steps 4–6 and anchors every 64 of the trail's
     own steps. The replay tracks (a, b) mod r: an image sigma(R_i) has coefficients multiplied by
     (−1)^s·lambda^(−k), canon multiplies them by (−1)^s·lambda^k, and an escape doubles the
     coefficients of canon(M). By Lemma R each replay reaches the recorded point.
   - The replays give the common canonical point as [a1]P + [b1]Q = [a2]P + [b2]Q. If
     b1 ≠ b2, output x = (a2 − a1)/(b1 − b2) mod r after checking [x]P = Q. If b1 = b2
     the run fails. By Lemma D this has probability 1/r.

**Merged trails stay merged.** Outside fruitless cycles the next class depends only on the
current class. Escapes fire only after an exact check, and they exit from 2·canon(M),
which depends only on the cycle's set of classes, not on where the walk entered the cycle
or when it noticed. The 64-bit filters only decide when an exact check runs. So two trails
that reach the same class follow the same classes afterwards and reach the same first
distinguished class.

**Lemma R (replays reach the recorded point).** The points of a trail, with repeated
traversals of cycles removed, depend only on its start point, not on when its cycles are
detected.

*Proof.* While no cycle is entered, the next point depends only on the current point. If
the trail enters a cycle that contains a distinguished point, it reaches that point within
one traversal, and no detection fires earlier. The 2-cycle rule fires only on a return to
the class of two steps back, and an anchor match only on a return to the anchor's class.
In both cases every point of the cycle has already been produced and tested. Otherwise
the trail repeats the cycle until a detection fires, and the escape 2·canon(M) depends
only on the cycle's classes. The distinguished-point test reads only the key. So any
schedule that anchors each walk at least once every 64 of its steps leads to the same
first distinguished point, and every coefficient pair the replay carries represents its
current point. ∎

Replaying a trail of n < TRAIL_MAX steps takes fewer than 2^7·n unbatched operations. Per
step that is:
- one evaluation;
- at most two more for a 2-cycle check;
- at most 66 for an anchor check, which fires at most once per 65 steps unless it escapes;
- at most 128 extra steps in each cycle before its escape, where each cycle took at least
  2 of the trail's steps.

**Lemma D (degenerate matches).** The first key match has b1 = b2 with probability exactly
1/r.

*Proof.* A trail records once, so the two records come from different trails t1 ≠ t2,
possibly of the same slot. Condition on all start points and adders. That fixes every
trail's points, and so which records match first.
- Given its start point W = [a + b·x]P, a trail's coin b is uniform on [0, r), because
  each b determines exactly one a.
- Coins of different trails are independent.
- Along a replay, b changes as b ← u·b + d, where u is a unit (±lambda^j) or twice
  one, and u, d depend only on the path.

So the canonical b at the distinguished point is U_i·b_ti + D_i with U_i invertible mod r,
and P[b1 = b2] = 1/r. ∎

## 2. Success probability

Let T_steps = K·L. Then tau = T_steps^2/(2N) = 0.530300.

**H1 (heuristic, disclosed).** The class map [W] → [canon(W) + R_h([W])], with uniform
random adders, behaves like a random mapping on the N classes, apart from the r-adding
penalty. That penalty is a factor of about (1 − 1/(6A)) in the birthday exponent (Teske
2001; Bernstein–Lange 2012; Bos–Kleinjung–Lenstra 2010). Under H1, the probability that
the first T_steps steps contain no useful class collision is at most
exp(−tau·(1 − eps)). Here eps = 4/A + 2·(3/A) + 2·(192/A^2) = 0.00015268. Its three
margins are:
- the adder penalty, taken as 4/A;
- 2-cycles at a rate of at most 1/A, with 3 non-fresh states each;
- short cycles at a rate of at most 1/A^2, with 192 wasted states each.

**H2 (rates).**
- 2-cycles are entered at a rate of at most 1/A per step. The literature gives 1/(6A).
  An earlier package of ours (rho-order6, not included) measured 7.2e−5 and 9.2e−5 per
  step on real secp256k1 walks at A = 2^11, against 1/(6·2^11) = 8.1e−5.
- Cycles of length 3..64 are entered at a rate of at most 1/A^2.
- Longer cycles with no distinguished point arise at a rate of at most A^−32. They are
  handled by abandonment, which is charged.

**Detection and resolution.** Consider a class collision among the first T_steps steps
between two different trails.
- From the shared class both trails follow the same classes to the same first
  distinguished class, since merged trails stay merged.
- Every trail finishes, so whichever trail records second triggers resolution.
- Lemma R shows that the replays reach that point, and Lemma D shows that the equation is
  degenerate with probability exactly 1/r.
- Fresh start points only add points, and the bound does not count them.

An earlier version continued walks through distinguished points. There, a fruitless cycle
that contains a distinguished point could store that point repeatedly with equal
coefficients. Restarting at every distinguished point removes that case: each record comes
from a different trail with its own coins, and such a cycle costs one ordinary record,
inside DP_MAX.

The failure terms are:
- **Collision inside one trail.** There are at most DP_MAX + 2K ≈ 2^67.75 trails, each with
  fewer than TRAIL_MAX ≈ 2^66 steps. The chance that any trail meets one of its own
  earlier classes, outside the short cycles handled by step 5, is at most
  (DP_MAX + 2K)·TRAIL_MAX^2/N ≈ 2^−53.7.
- **Abandoned trail.** A trail is abandoned only after more than CAP − 64 steps without a
  distinguished point outside short cycles. That has probability at most
  e^−63.99 ≈ 2^−92.3 per trail, so 2^−24.6 over all trails.
- **Degenerate first match:** exactly 1/r < 2^−252 (Lemma D).
- **DP_MAX or abandonment overflow:** under 2^−1000 by Chernoff.
- **Exceptional additions:** under 2^−110.

Together these are below delta = 2^−20. Therefore:

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015268)) − 2^−20 = 0.411523 ≥ 0.39.**

The budget is 1.0299·√N steps, against 0.9943·√N for success exactly 0.39. The 0.0215
of slack in the model bound covers the sampling error of the toy evidence below.

**Evidence for H1, H2 and the procedure.** `experiment.py` runs the procedure of
section 1.3 end to end on random prime-order toy curves
y^2 = x^3 + 7 over p ≡ 1 mod 3 (`"b": 7`). It uses:
- the same class key, canonical rule, adder images and index (keys use Montgomery
  representatives, any fixed bijection being equivalent);
- both cycle rules and their escapes;
- trails from fresh coins that end at their first distinguished point, and the trail cap;
- resolution by replaying both trails from their coins.

It never tracks coefficients while walking and has no oracle collision table.
- The live walk anchors on a global round counter, as the batched kernel does. The replay
  anchors on the trail's own step count, so every replay also checks Lemma R.
- Each trial runs L steps per slot plus the finishing phase at the stated tau. It counts
  as a success only if the first key match resolves to a verified x.
- Each trail start is charged as one of the slot's L steps, which is conservative.

<!-- toy-evidence -->
Offline runs (`evidence/trails/`, pooled in `summary.json`): 5,200 trials on 5 curves with
N = 2^26.8, 2^29.1, 2^30.8, 2^31.3, 2^32.5 classes,
adder tables of 2^11 and 2^16, 16 slots and theta = 2^−5.
- Success of the whole procedure at tau = 0.5303: 2,173 of 5,200, 0.4179, Wilson 95% interval
  [0.405, 0.431]. The model 1 − exp(−tau) is 0.412.
- Every success was verified ([x]P = Q). Degenerate first matches: 0. Replays that did not
  reach the recorded point: 0. Every match was replayed; about 4.7 million trails were run.
- Abandoned trails: 1. Longer fruitless cycles escaped: 6. The finishing phase added at most
  2.6% to the K·L steps.
- 2-cycle entry rate per step:
  - 7.60e-05, 7.92e-05, 7.74e-05, 8.02e-05 with A = 2^11 (1/(6A) = 8.14e-05);
  - 2.03e-06 with A = 2^16 (1/(6A) = 2.54e-06).
<!-- /toy-evidence -->

The declared initial input (`bits` 28, `b` 7, 100 trials, seed 2801, the package defaults
otherwise) runs in a few seconds. It reports 42/100 = 0.42, Wilson 95%
interval [0.328, 0.518], against the model 0.412, with 0 degenerate matches and 0 replay
mismatches. With 100 trials the standard error is about 0.05, so it is a reproducibility
check. The success claim rests on the pooled offline trials above, and any input within
the schema can spot-check them.

Limits: the toy curves have N up to about 2^33, against 2^253 in production. There theta is 2^−5,
so start points and the finishing phase are a few percent of the steps; in production
they are below 2^−30 of them. H1 is an extrapolation in N, supported by the size-free
random-mapping argument and the cited literature. It is not proved for this walk.

## 3. Computation

**No reference operation.** Each step is charged its full measured instruction count
divided by C = 1000, as COST_MODEL.md's H + R/C with H = 0. No 6M reference operation is
claimed, so nothing is charged twice.

**Field arithmetic.** `measurement/secp_field.h` implements multiplication and squaring
modulo p = 2^256 − c, c = 2^32 + 977, as GCC inline AArch64 assembly in general-purpose
registers.
- Products are 16 `mul` + 16 `umulh` (squaring: 10 + 10 with doubling). The high half is
  folded twice by c. A final carry selects between y and y − c, which gives a fully
  reduced output without a separate comparison with p.
- `measurement/secp_sub.h` holds branchy subtractions that subtract c (that is, add p
  mod 2^256) only on a borrow, and a fused a − b − c.
- Standalone costs in the same toolchain, per call including about 7 instructions of
  call and loop overhead: mul 108, sqr 85, sub 22 on average. For comparison, the pinned
  libsecp256k1 5x52 multiplication is 179. Inside the kernel the routines are inlined.
- Correctness:
  - The routines matched Python big-integer arithmetic on about 194,000 reduced and
    199,000 arbitrary 256-bit operand pairs. These included constructed cases for the
    rare second-fold carry.
  - The replica check below compares every final point of 2^16 walks after 2^19 steps
    with Python big-integer arithmetic, on two input sets.

Inline assembly is ordinary scalar AArch64 code. Callgrind counts each executed
instruction exactly as for compiled C, which is how the organizer measured blst's assembly.

**Per-step measurement (H3).** `measurement/step.c` implements steps 3–6 for a batch of
2^16 walks with 2^16 adders, including restarts at distinguished points and the
trail cap. In this harness a new trail takes its start point from spare points in the
input file. In production it comes from fresh coins, and that cost is charged below as
trail restarts.
- It was compiled with `-O3 -fno-tree-vectorize -mgeneral-regs-only` in the organizer's
  Docker image (GCC 14.2.0, Valgrind 3.24.0), and counted with Callgrind `Ir`.
- Counts are differenced between 8 and 4 batch-steps of 2^16 walks, 2^18 steps.
- Every instruction inside the measured loop is charged.

| Input set | total Ir/step | 2-cycles in the larger run |
| ---: | ---: | ---: |
| seed 1 | 779.992 | 0 |
| seed 2 | 780.008 | 0 |

This covers:
- 5 multiplications and 1 squaring;
- 5 subtractions, one of them fused;
- the fraction approximations and the k decision;
- key, parity and adder-image selection;
- the 2-cycle, anchor and distinguished-point filters, the trail-cap counter, and any
  confirmations, escapes and restarts;
- the batch's prefix stores and loads;
- one Fermat inversion per 2^16 walks (libsecp256k1's addition chain, 255 squarings and
  15 multiplications, about 22,000 Ir, or 0.34 per step);
- loop control.

**795 instructions** are charged per step, so each step costs **0.795 units**.
`measurement/walk_ref.py check` replays all walks in Python integers and matches the C
final points and event counts exactly on both input sets.

Distinguished points never occur at theta = 2^−60 in such a sample. So the restart and
abandonment paths were checked on small builds of the same source and replica
(`evidence/measurement/trail-paths.json`), with 2^6 adders and 2^12 walks:
- With DP_BITS = 4, the run had 32,502 restarts at distinguished points and 1,144 2-cycles.
- With DP_BITS = 12 and a one-period cap, it had 4,096 abandonments.

In both, the replica matched every final point and every count.

Rare paths fit in the 15.0-instruction margin:
- **2-cycles.** At most 1/A per step × 3 unbatched operations × 2^15 instructions,
  that is at most 1.5 instructions per step.
- **Restarts and abandonment.** A restart in the kernel copies a start point and selects
  its adder, about 50 instructions, at rate theta per step: below 2^−50 instructions per
  step. Its start point is charged in the ledger.
- **Short cycles and exact-k fallbacks.** Below 2^−15 instructions per step.

A separately measured unbatched step (affine add with its own inversion, plus adder
selection) costs 22,913.1 and 22,914.0 Ir. Every unbatched operation below is charged 2^15.

**Total.** S_max = K·(L + TRAIL_MAX) ≤ 2^126.749960 steps.

| Component | Bound (units) | log2 |
| --- | --- | ---: |
| walk steps | 0.795·S_max | 126.418986 |
| adder table | 2^16 adders × 2 scalar mults × 768 unbatched ops × 2^15 / C | 31.62 |
| initial walk starts | K × 2 scalar mults × 768 × 2^15 / C | 47.62 |
| trail restarts | (DP_MAX + K) × 2 scalar mults × 768 × 2^15 / C | 83.37 |
| distinguished points | DP_MAX = 2^67.75 entries × 2^20 instructions / C | 77.78 |
| collision replay | 2 trails × TRAIL_MAX steps × 2^7 unbatched ops × (2^15 + 2^10) / C | 79.08 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 16.03 |
| input validation | 2^24 instructions | 14.03 |

Rounded up, **time_log2 = 126.418986**. Preprocessing (adders and initial walk starts) is
included. All processors' work is summed, and elapsed time plays no part.

## 4. Memory and advice

Peak memory is at most the sum of:
- 2·DP_MAX·128 bytes for the distinguished-point table (key, the trail's coins, hashing
  overhead);
- K·576 bytes of walk state (points, keys, coins, trail counter, batch buffers);
- 2^16·6·64 bytes of adder images;
- 16 MiB.

That is 2^75.74996 bytes. Nonuniform advice is 0 bytes.

## 5. Limitations

- H1 is the standard but unproved random-mapping heuristic. It is supported by the
  literature and toy-curve experiments only. The model bound is 0.4115, 0.0215 above
  0.39.
- Instruction counts come from one implementation and compiler. The sample is 2^18 steps
  on each of two input sets. The branchy subtractions and comparisons were measured, not
  bounded analytically; each subtraction's worst case is 6 instructions above its
  no-borrow path.
- The field routines were tested, not formally verified. The replica check covers 2^19
  production steps per input set.
- About 2^67.7 distinguished points and 2^126.4 work are needed. This is not a practical
  attack and does not weaken secp256k1.
- The toy experiment runs the whole procedure, including resolution, but only at toy
  scale: theta = 2^−5, 16 slots, and N far below production. Its trail starts and its
  finishing phase are a few percent of the steps there, and negligible in production.
- Distinguished points never occur in a production-parameter measurement sample. The
  kernel's restart and abandonment paths were checked on small builds of the same source.

## References

- Pollard (1978).
- van Oorschot–Wiener (J. Cryptology 1999).
- Teske (Math. Comp. 2001).
- Duursma–Gaudry–Morain (Asiacrypt 1999).
- Bos–Kleinjung–Lenstra (ANTS 2010).
- Bernstein–Lange–Schwabe (PKC 2011).
- Bernstein–Lange (ANTS 2012).
- Bos–Costello–Miele (PKC 2014).

## Reproduction

```sh
python3 -I submissions/secp256k1-ecdlp/rho-sigma6-asm/accounting.py
python3 -I submissions/secp256k1-ecdlp/rho-sigma6-asm/measurement/frac_check.py 1
bash submissions/secp256k1-ecdlp/rho-sigma6-asm/measurement/reproduce.sh 1   # Docker, linux/arm64
echo '{"bits":28,"trials":100,"seed":2801,"b":7,"adder_bits":11,"dp_bits":5,"walks":16,"tau_e4":5303}' | python3 -I submissions/secp256k1-ecdlp/rho-sigma6-asm/experiment.py
```

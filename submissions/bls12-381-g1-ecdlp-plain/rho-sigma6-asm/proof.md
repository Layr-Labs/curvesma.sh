# bls12-381-g1-ecdlp-plain: parallel Pollard rho on order-6 classes with a sigma walk, all arithmetic in instructions

Claim: for uniform x in [0,r) and T=[x]P0 in G1, the algorithm below returns x with
probability at least 0.39. It uses at most **2^126.312725** target-group-operation
equivalents in total (score 5052509/40000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^75.178504 bytes. Preprocessing is
2^44.619736 units, and no nonuniform advice is used. This is a classical cold attack on
one target.

The success bound rests on the disclosed heuristic H1, the standard random-mapping model
for r-adding walks on automorphism classes, and on the rate bound H2. The cost rests on
the instruction counts in H3, measured with the organizer's toolchain. This is a
carefully accounted generic attack. It is not a practical break.

**Cost route.** The walk does not use the 6M reference operation. Its five field
multiplications and one squaring per step are the package's own 6×64-bit Montgomery
multiplication in GCC inline AArch64 assembly (`measurement/bls_field.h`, the same
convention as blst's `mul_mont_384`). Every instruction of the step is charged at 1/C
with C = 2000 (H = 0 in COST_MODEL's H + R/C). COST_MODEL.md: "Faster arithmetic and
lookup algorithms are welcome; justify their actual work". A step measures 2,174
instructions, so it is charged 1.0975 units rather than 1 + R/C = 1.14. The same walk with
the reference operation scores 126.367538 (an earlier version, not included).

What is new relative to the usual order-6 walk is in section 1.2. The walk never moves
a point to its canonical representative, so it never pays the multiplication by beta.
It stores the actual point W and adds sigma(R_i), where sigma is the automorphism that
maps the canonical point to W. The class sequence is unchanged, so the standard
analysis applies unchanged. Which rotation is canonical is decided from 64-bit
fixed-point approximations with an exact fallback, and Lemma 1 bounds their error.

**Distinguished points and resolution.** Walks run in trails. Each trail starts from fresh
coins and ends at its first distinguished point, which it records once together with its
coins. A key match is resolved by replaying the two trails (section 1.3, Lemmas R and D).
A fruitless cycle therefore cannot store the same point twice, and a degenerate match has
probability exactly 1/r.

Files:
- `accounting.py` recomputes every number with exact integers (output in
  `evidence/accounting.json`).
- `measurement/` holds the per-step instruction-count harness, the field code
  `bls_field.h`, the Python replica `walk_ref.py`, and the approximation check
  `frac_check.py`.
- `evidence/measurement/` holds the counts on two input sets and the checks of the
  restart and abandonment paths.
- `experiment.py` is the declared toy-curve experiment, running the whole procedure;
  `evidence/trails/` holds larger offline runs of the same file.

## 1. Algorithm

### 1.1 Notation and parameters

p, r, P0 and the curve y^2 = x^3 + 4 come from the protected profile. r is prime, and
p ≡ r ≡ 1 (mod 3).

The endomorphism uses beta = 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe,
a primitive cube root of unity mod p, and lambda = 0x73eda753299d7d483339d80809a1d804a7780001fffcb7fcfffffffe00000001
= −z^2 mod r, where z = −0xd201000000010000. Then [lambda](x,y) = (beta·x, y) on G1. This was
checked on P0 (`measurement/walk_ref.py` uses it; [lambda]P0 = (beta·x0, y0)).

The automorphism group G6 = {±1, ±lambda, ±lambda^2} acts on G1. Each nonzero point lies
in a class of exactly 6 points, so there are (r−1)/6 nonzero classes. Write N = (r−1)/6 + 1.

M(v) = v·2^384 mod p is blst's Montgomery representative, a fixed bijection of Fp onto
[0,p). Every comparison and bit test below is on M-values, read as integers.

Parameters:
- A = 2^16 adders.
- K = 2^30 walks: 2^16 processors × 2^14-walk batches.
- theta = 2^−60.
- L = 89663270682989036334838442605 ≈ 2^96.179.
- CAP = 64/theta. Trails have fewer than TRAIL_MAX = CAP + 64 steps.
- DP_MAX = 2·theta·S_max + 2K, where S_max = K·(L + TRAIL_MAX).

### 1.2 Classes, keys and the sigma walk

For an affine point W = (x, y) ≠ O:

- **Class key.** key(W) = M(y) if M(y) is even, else p − M(y), so key(W) = M(±y) and is even.
  - key is constant on the class, since G6 changes y only by sign.
  - key determines the class. y^2 = x^3 + 4 fixes x^3, hence x up to a factor beta^j, and
    ±y fixes the sign.
- **Canonical point.** canon(W) = (beta^k·x, (−1)^s·y), where k ∈ {0,1,2} minimises
  M(beta^k·x) and s = M(y) mod 2.
  - The three values are distinct for x ≠ 0, because beta^j·x differ and M is a bijection.
  - x = 0 does not occur in G1: (0, ±2) have order 3.
  - p is odd, so M(−y) = p − M(y) has the other parity. canon(W) therefore has an even
    M(y) and is a fixed point of the class.
- **sigma.** sigma_W = (x → beta^(−k)·x, y → (−1)^s·y), so that W = sigma_W(canon(W)).
  sigma_W is a group automorphism, equal to [(−1)^s·lambda^(−k)].
- **Walk step.** i(W) = (key(W) >> 1) mod 2^16, and

      f(W) = W + sigma_W(R_i(W)).

  Because sigma_W is an automorphism, f(W) = sigma_W(canon(W) + R_i). So the class of f(W)
  is the class of canon(W) + R_i. That is exactly the standard class walk
  [W] → [canon(W) + R_h([W])] with the class function h = i. Every statement about the
  class sequence below is a statement about that standard walk.
- **Images of each adder.** The six images sigma(R_i) = (beta^(−k)·X_i, ±Y_i) are tabulated
  per adder, so a step reads sigma_W(R_i) from T6[i][k][s] and no point is ever rotated
  or negated.
- **Distinguished point.** W is distinguished when bits 17..76 of key(W) are zero, which
  is a class property.

**Deciding k without beta·x (Lemma 1).** Let a = M(x) = a_0 + a_1·2^64 + … + a_5·2^320.
Let C_0 = ⌊2^512/p⌋ and C_1 = ⌊2^512·beta/p⌋, with 64-bit limbs C_{j,m}. Define

    g_j = ( Σ_{i+m=7} lo(a_i·C_{j,m}) + Σ_{i+m=6} hi(a_i·C_{j,m}) ) mod 2^64,   j = 0, 1,

where lo and hi are the low and high words of the 128-bit product. For C_0 only m = 1, 2
are nonzero, so g_0 takes 3 multiplications and g_1 takes 12.

Let f_j = 2^64·M(beta^j·x)/p, a real number in (0, 2^64). Note M(beta·x) = beta·a mod p.

**Lemma 1.** If g_j ≤ 2^64 − 14, then f_j ∈ [g_j, g_j + 13) for j = 0, 1.

*Proof.* f_j = 2^64·frac(a·c_j) with c_0 = 1/p and c_1 = beta/p. Here a·c_j differs from
a·C_j/2^512 by less than a/2^512 < 2^−131. Limb 7 of a·C_j mod 2^512 equals g_j plus the
carry from the dropped terms. Those terms are the low words with i+m = 6 and all terms
with i+m ≤ 5. Their total is below 6·2^448 + 6·2^448·(1 + 2^−64 + …) < 13·2^448, so the
carry is at most 12. Products with i+m ≥ 8, and high words with i+m = 7, are multiples of
2^512 and vanish mod 1. The wrap-around case is excluded by g_j ≤ 2^64 − 14. ∎

M(x) + M(beta·x) + M(beta^2·x) ≡ 0 mod p, with each term in (0, p), so the sum is p or 2p
and f_2 = 2^64·m − f_0 − f_1 for m ∈ {1, 2}. With g_2 = (−g_0 − g_1 − 26) mod 2^64, the
real f_2 lies in (g_2, g_2 + 26] whenever g_2 ≤ 2^64 − 28. The kernel accepts the index of
the least g_j as k only if all three g_j ≤ 2^64 − 28 and the least is more than 26 below
both others. Then f_k ≤ g_k + 26 < g_l ≤ f_l, so k is exact. Otherwise it computes beta·x
and compares exactly.

The fallback probability per step is at most about 3·2·27/2^64 + 3·28/2^64 < 2^−55. It
never occurred in 2^19 measured steps. `measurement/frac_check.py` checks the floor gap
⌊f_j⌋ − g_j ∈ [0,12] on 200,000 random and edge-case inputs per seed; the largest gap seen
is 6.

### 1.3 Procedure

1. **Validate T.** Check it is on the curve and [r]T = O. If T = O, output 0.
2. **Adders.** For each adder i < 2^16:
   - draw c_i, d_i uniform in [0,r);
   - set R_i = [c_i]P0 + [d_i]T;
   - tabulate its six images, using one multiplication by beta and one negation.
3. **Trails.** Each of the K walk slots runs a sequence of trails.
   - A trail starts from fresh coins a, b uniform in [0,r), with W = [a]P0 + [b]T. The
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
   simultaneous inversion over the 2^14 walks of a batch (3M), then lambda, lambda^2 and
   y3 (2M + 1S), with one inversion per batch.
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
   - The replays give the common canonical point as [a1]P0 + [b1]T = [a2]P0 + [b2]T. If
     b1 ≠ b2, output x = (a2 − a1)/(b1 − b2) mod r after checking [x]P0 = T. If b1 = b2
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
- Given its start point W = [a + b·x]P0, a trail's coin b is uniform on [0, r), because
  each b determines exactly one a.
- Coins of different trails are independent.
- Along a replay, b changes as b ← u·b + d, where u is a unit (±lambda^j) or twice
  one, and u, d depend only on the path.

So the canonical b at the distinguished point is U_i·b_ti + D_i with U_i invertible mod r,
and P[b1 = b2] = 1/r. ∎

## 2. Success probability

Let T_steps = K·L be the main-phase steps, each walk's first L steps. Then
tau = T_steps^2/(2N) = 0.530300.

**H1 (heuristic, disclosed).** The class map [W] → [canon(W) + R_h([W])], with uniform
random adders, behaves like a random mapping on the N classes, apart from the known
r-adding penalty. An r-adding walk loses about a factor (1 − 1/(m·A)) in the birthday
exponent, with A the number of adders and m = 6 the class size (Teske 2001;
Bernstein–Lange 2012; Bos–Kleinjung–Lenstra 2010 eq. (2) for m = 2).

Under H1, the probability that the first T_steps steps contain no useful class collision
is at most exp(−tau·(1 − eps)). Here eps = 4/A + 2·(3/A) + 2·(192/A^2) = 0.00015268. It
collects three margins:

- **Adder penalty, 4/A.** This is 24 times the 1/(6A) estimate.
- **Fruitless 2-cycles, 2·(3/A).** H2 bounds the entry rate by 1/A per step, and each
  event yields at most 3 non-fresh states.
- **Short cycles, 2·(192/A^2).** H2 bounds the rate by 1/A^2, and detection wastes at most
  192 states per event.

**H2 (rates).**
- 2-cycles are entered at a rate of at most 1/A per step. The literature gives 1/(6A)
  for order-6 classes (Bos–Costello–Miele 2014).
- Cycles of length 3..64 are entered at a rate of at most 1/A^2.
- Cycles longer than 64 with no distinguished point arise at a rate of at most A^−32.
  They cost at most CAP steps each through abandonment, which is charged.

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
- **Collision inside one trail.** There are at most DP_MAX + 2K ≈ 2^67.18 trails, each with
  fewer than TRAIL_MAX ≈ 2^66 steps. The chance that any trail meets one of its own
  earlier classes, outside the short cycles handled by step 5, is at most
  (DP_MAX + 2K)·TRAIL_MAX^2/N ≈ 2^−53.1.
- **Abandoned trail.** A trail is abandoned only after more than CAP − 64 steps without a
  distinguished point outside short cycles. That has probability at most
  e^−63.99 ≈ 2^−92.3 per trail, so 2^−25.2 over all trails.
- **Degenerate first match:** exactly 1/r < 2^−252 (Lemma D).
- **DP_MAX or abandonment overflow:** under 2^−1000 by Chernoff.
- **Exceptional additions:** under 2^−110.

Together these are below delta = 2^−20. Therefore:

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015268)) − 2^−20 = 0.411523 ≥ 0.39.**

The budget is 1.0299·√N steps, against 0.9943·√N for success exactly 0.39. The 0.0215
of slack in the model bound covers the sampling error of the toy evidence below. The
bound is exact arithmetic on H1, not an empirical estimate.

**Evidence for H1, H2 and the procedure.** `experiment.py` runs the procedure of
section 1.3 end to end on toy curves y^2 = x^3 + 4 over
p ≡ 1 mod 3, in prime-order subgroups since that curve always has 3-torsion (`"b": 4`). It uses:
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
Offline runs (`evidence/trails/`, pooled in `summary.json`): 4,800 trials on 5 curves with
N = 2^28.8, 2^27.5, 2^30.0, 2^35.5, 2^35.5 classes,
adder tables of 2^11 and 2^16, 16 slots and theta = 2^−5.
- Success of the whole procedure at tau = 0.5303: 1,969 of 4,800, 0.4102, Wilson 95% interval
  [0.396, 0.424]. The model 1 − exp(−tau) is 0.412.
- Every success was verified ([x]P = Q). Degenerate first matches: 0. Replays that did not
  reach the recorded point: 0. Every match was replayed; about 8.0 million trails were run.
- Abandoned trails: 1. Longer fruitless cycles escaped: 9. The finishing phase added at most
  2.1% to the K·L steps.
- 2-cycle entry rate per step:
  - 7.65e-05, 8.02e-05, 7.59e-05, 7.65e-05 with A = 2^11 (1/(6A) = 8.14e-05);
  - 2.36e-06 with A = 2^16 (1/(6A) = 2.54e-06).
<!-- /toy-evidence -->

The declared initial input (`bits` 28, `b` 7, 100 trials, seed 2801, the package defaults
otherwise) runs in a few seconds. It uses prime-order y^2 = x^3 + 7 curves, which exercise
the same order-6 class code. For b = 4, fields of at most 28 bits give subgroups of only
2^19 to 2^23 classes. There each slot's final trail adds 20 to 40% to K·L, so success runs
high; the offline b = 4 runs above use 36- to 40-bit fields. It reports 42/100 = 0.42,
Wilson 95% interval [0.328, 0.518], against the model 0.412, with 0 degenerate matches and
0 replay mismatches. With 100 trials the standard error is about 0.05, so it is a
reproducibility check. The success claim rests on the pooled offline trials above, and any
input within the schema can spot-check them.

Limits: the toy curves have N up to about 2^36, against 2^252 in production. There theta is 2^−5,
so start points and the finishing phase are a few percent of the steps; in production
they are below 2^−30 of them. H1 is an extrapolation in N, supported by the size-free
random-mapping argument and the cited literature. It is not proved for this walk.

## 3. Computation

**No reference operation.** Each step is charged its full measured instruction count
divided by C = 2000, as COST_MODEL.md's H + R/C with H = 0. No 6M reference operation is
claimed, so nothing is charged twice.

**Field arithmetic.** `measurement/bls_field.h` implements Montgomery multiplication and
squaring with R = 2^384 as GCC inline AArch64 assembly in general-purpose registers.
- Multiplication interleaves six word-by-word Montgomery reduction rows with the
  product. The cancelled low limb of each row is replaced by a carry test, and one
  final conditional subtraction (`csel`) makes the output fully reduced.
- Squaring computes the 21 distinct products with the doubling fused into the diagonal,
  then a 12-word Montgomery reduction.
- Standalone costs in the same toolchain, per call including 6 instructions of loop and
  call overhead: mul 334 and sqr 279, against blst's 383 and 325 in the same harness.
  Inlined in the kernel, mul is about 325.
- Correctness:
  - The routines matched Python big-integer arithmetic on 380,000 cases, with seven
    outputs each, including aliasing and products landing on 0, 1, p − 1 and p − 2.
  - The replica check below compares every final point of 2^14 walks after 2^18 steps
    on two input sets.
- blst is still used for the one inversion per batch (`blst_fp_inverse`) and for input
  conversion. Those instructions are counted like any others.

Inline assembly is ordinary scalar AArch64 code. Callgrind counts each executed
instruction exactly as for compiled C, which is how the organizer measured blst's assembly.

**Per-step measurement (H3).** `measurement/step.c` implements steps 3–6 for a batch of
2^14 walks with 2^16 adders, including restarts at distinguished points and the
trail cap. In this harness a new trail takes its start point from spare points in the
input file. In production it comes from fresh coins, and that cost is charged below as
trail restarts.
- It was compiled with `-O3 -fno-tree-vectorize -mgeneral-regs-only` in the organizer's
  Docker image (GCC 14.2.0, Valgrind 3.24.0), and counted with Callgrind `Ir`.
- Counts are differenced between 16 and 8 batch-steps of 2^14 walks, 2^17 steps.
- Every instruction inside the measured loop is charged.

| Input set | total Ir/step | 2-cycles in the larger run |
| ---: | ---: | ---: |
| seed 1 | 2174.255 | 0 |
| seed 2 | 2174.228 | 1 |

This covers:
- 5 multiplications and 1 squaring;
- all field subtractions (5 per step, one of them fused a−b−c), 6-limb inline-assembly
  sequences that add p back only on a borrow;
- the two fraction approximations and the k decision;
- key, parity and adder-image selection;
- the 2-cycle, anchor and distinguished-point filters, the trail-cap counter, and any
  confirmations, escapes and restarts;
- the batch's prefix stores and loads and the one inversion per 2^14 walks (26,246 Ir,
  or 1.6 per step);
- loop control.

**2,195 instructions** are charged per step, so each step costs **1.0975 units**.
`measurement/walk_ref.py check` replays all walks in Python integers and matches the C
final points and event counts exactly on both input sets.

Distinguished points never occur at theta = 2^−60 in such a sample. So the restart and
abandonment paths were checked on small builds of the same source and replica
(`evidence/measurement/trail-paths.json`), with 2^6 adders and 2^12 walks:
- With DP_BITS = 4, the run had 49,926 restarts at distinguished points and 1,903 2-cycles.
- With DP_BITS = 12 and a one-period cap, it had 4,092 abandonments.

In both, the replica matched every final point and every count.

Rare paths fit in the 20.7-instruction margin:
- **2-cycles.** At most 1/A per step × 3 unbatched operations × 2^15 instructions,
  that is at most 1.5 instructions per step.
- **Restarts and abandonment.** A restart in the kernel copies a start point and selects
  its adder, about 50 instructions, at rate theta per step: below 2^−50 instructions per
  step. Its start point is charged in the ledger.
- **Short cycles and exact-k fallbacks.** Below 2^−15 instructions per step.

A separately measured unbatched step (affine add with its own inversion, plus adder
selection) costs 28,553.8 on both input sets Ir. Every unbatched operation below is charged 2^15.

**Total.** S_max = K·(L + TRAIL_MAX) ≤ 2^126.178504 steps.

| Component | Bound (units) | log2 |
| --- | --- | ---: |
| walk steps | 1.0975·S_max | 126.312725 |
| adder table | 2^16 adders × 2 scalar mults × 768 unbatched ops × 2^15 / C | 30.62 |
| initial walk starts | K × 2 scalar mults × 768 × 2^15 / C | 44.62 |
| trail restarts | (DP_MAX + K) × 2 scalar mults × 768 × 2^15 / C | 81.80 |
| distinguished points | DP_MAX = 2^67.18 entries × 2^20 instructions / C | 76.21 |
| collision replay | 2 trails × TRAIL_MAX steps × 2^7 unbatched ops × (2^15 + 2^10) / C | 78.08 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 15.03 |
| input validation ([r]T = O check) | 2 scalar mults + 2^24 instructions | 15.03 |

Rounded up, **time_log2 = 126.312725**. Preprocessing (adders and initial walk starts) is
included. All processors' work is summed, and elapsed time plays no part.

## 4. Memory and advice

Peak memory is at most the sum of:
- 2·DP_MAX·128 bytes for the distinguished-point table (key, the trail's coins, hashing
  overhead);
- K·640 bytes of walk state (points, keys, coins, trail counter, batch buffers);
- 2^16·6·96 bytes of adder images;
- 16 MiB.

That is 2^75.178504 bytes. Nonuniform advice is 0 bytes.

## 5. Limitations

- H1 is the standard but unproved random-mapping heuristic for additive class walks. It
  is supported here by the literature and by toy-curve experiments only. The model bound
  is 0.4115, 0.0215 above 0.39.
- Instruction counts come from one implementation, compiler and convention. The sample
  is 2^18 steps on each of two input sets.
- The field multiplication is the package's own assembly, tested against Python big
  integers but not formally verified.
- The variable-time paths (branchy subtractions, comparisons, the exact fallback) were
  measured, not bounded analytically. Lemma 1 makes the fallback's rarity analytic.
- About 2^67 distinguished points and 2^126.3 work are needed. This is not a practical
  attack and does not weaken BLS12-381's security assumptions.
- The pairing (MOV) route to GT is not used.
- The toy experiment runs the whole procedure, including resolution, but only at toy
  scale: theta = 2^−5, 16 slots, and N far below production. Its trail starts and its
  finishing phase are a few percent of the steps there, and negligible in production.
- Distinguished points never occur in a production-parameter measurement sample. The
  kernel's restart and abandonment paths were checked on small builds of the same source.

## References

- Pollard (1978).
- van Oorschot–Wiener, *Parallel collision search* (J. Cryptology 1999).
- Teske, *On random walks for Pollard's rho method* (Math. Comp. 2001).
- Duursma–Gaudry–Morain (Asiacrypt 1999).
- Bos–Kleinjung–Lenstra, *On the use of the negation map in Pollard rho* (ANTS 2010).
- Bernstein–Lange–Schwabe, *On the correct use of the negation map in the Pollard rho
  method* (PKC 2011).
- Bernstein–Lange, *Two grumpy giants and a baby* (ANTS 2012).
- Bos–Costello–Miele, *Elliptic and hyperelliptic curves: a practical security analysis*
  (PKC 2014).

The construction and accounting above are complete and do not depend on fetching these.

## Reproduction

```sh
python3 -I submissions/bls12-381-g1-ecdlp-plain/rho-sigma6-asm/accounting.py
python3 -I submissions/bls12-381-g1-ecdlp-plain/rho-sigma6-asm/measurement/frac_check.py 1
bash submissions/bls12-381-g1-ecdlp-plain/rho-sigma6-asm/measurement/reproduce.sh 1   # Docker, linux/arm64
echo '{"bits":28,"trials":100,"seed":2801,"b":7,"adder_bits":11,"dp_bits":5,"walks":16,"tau_e4":5303}' | python3 -I submissions/bls12-381-g1-ecdlp-plain/rho-sigma6-asm/experiment.py
# offline evidence runs used larger inputs, e.g. {"bits":40,"trials":700,"seed":4001}
```

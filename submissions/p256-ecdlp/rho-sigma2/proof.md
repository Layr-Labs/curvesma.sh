# p256-ecdlp: parallel Pollard rho on negation classes with a sigma walk

Claim: for uniform x in [0,r) and Q=[x]P on NIST P-256, the algorithm below returns x
with probability at least 0.39. It uses at most **2^127.273624** target-group-operation
equivalents in total (score 15909203/125000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^76.542441 bytes. Preprocessing is 2^48.619436 units, and no nonuniform advice is
used. This is a classical cold attack on one target.

The success bound rests on the disclosed heuristic H1, the standard random-mapping model
for r-adding walks on negation classes, and on its per-step form H1′, from which Lemma C derives the
fruitless-cycle rates. The cost rests on H3,
average instruction counts measured with the organizer's toolchain, and every run is
charged by a work meter with a Markov reserve. **There is no reference
operation on this track:** every instruction, including all P-256 field arithmetic and
point additions, is charged at 1/C with the borrowed normalisation C = 1000. This is a
carefully accounted generic attack. It is not a practical break.

**Distinguished points and resolution.** Walks run in trails. Each trail starts from fresh
coins and ends at its first distinguished point, which it records once together with its
coins. A key match is resolved by replaying the two trails (section 1.3, Lemmas R and D).
A fruitless cycle therefore cannot store the same point twice, and a degenerate match has
probability exactly 1/r.

Files:
- `accounting.py` recomputes every number with exact integers (output in
  `evidence/accounting.json`).
- `measurement/` holds the per-step instruction-count harness and the field code
  (`p256_field.h`, `p256_sub.h`), plus the Python replica `walk_ref.py`.
- `evidence/measurement/` holds the round-5 production totals and the checks of the
  restart and abandonment paths.
- `experiment.py` is the declared toy-curve experiment, running the whole procedure;
  `evidence/trails/` holds larger offline runs of the same file.

## 1. Algorithm

### 1.1 Notation and parameters

p, a = −3, b, r and P come from the protected profile. The cofactor is 1, and r is prime.
M(v) = v·2^256 mod p is the Montgomery representative, a fixed bijection of Fp onto
[0,p). Every comparison and bit test below is on M-values, read as integers.

Negation splits the nonzero points into (r−1)/2 classes {W, −W}. Write N = (r−1)/2 + 1.

Parameters:
- A = 2^16 adders.
- K = 2^32 walks: 2^16 processors × 2^16-walk batches.
- theta = 2^−60.
- L = 57695294864490661436112010649 ≈ 2^95.542.
- CAP = 64/theta. Trails have fewer than TRAIL_MAX = CAP + 64 steps.
- DP_MAX = 2·theta·S_max + 2K, where S_max = K·(L + TRAIL_MAX).

### 1.2 Classes, keys and the sigma walk

For W = (x, y) ≠ O:
- **Class key.** key(W) = M(x). It is constant on {W, −W} and determines the class.
- **Canonical point.** canon(W) = (x, (−1)^s·y), where s = M(y) mod 2. p is odd and y ≠ 0
  (no 2-torsion), so M(−y) = p − M(y) has the other parity, and canon picks one point of
  each class.
- **sigma.** sigma_W = [(−1)^s], so that W = sigma_W(canon(W)).
- **Walk step.** i(W) = key(W) mod 2^16, and

      f(W) = W + sigma_W(R_i(W)) = sigma_W(canon(W) + R_i).

  So the class of f(W) is the class of canon(W) + R_i, which is exactly the standard
  negation-class walk. The kernel reads sigma_W(R_i) from a table of R_i and −R_i, so no
  point is ever negated, and sigma costs one bit test.
- **Distinguished point.** W is distinguished when bits 16..75 of key(W) are zero.

### 1.3 Procedure

1. **Validate Q.** Check it is on the curve. If Q = O, output 0.
2. **Adders.** For each adder i < 2^16:
   - draw c_i, d_i uniform in [0,r);
   - set R_i = [c_i]P + [d_i]Q;
   - tabulate R_i and −R_i.
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
   simultaneous inversion over the 2^16 walks of a batch (3M), then lambda = dy·inv (M),
   lambda^2 (S) and y3 = lambda·(x1 − x3) (M), with one Fermat inversion per batch.
   - If the batch product is zero, the affected walk uses the unbatched complete formula.
5. **Fruitless cycles.** Each walk keeps two filters, both the low 64 bits of a key: the
   key of the point before its current point W, and an anchor, its key at the last batch
   step divisible by 64. All walks of a batch step together, so each walk is anchored once
   every 64 of its own steps. Because the filters keep only 64 bits they can match
   falsely. They only decide when an exact check runs.
   - *Order.* The distinguished-point test of step 6 comes first, on every point the walk
     produces: the new point, every look-ahead point and every escape point. A
     distinguished new point ends the trail, and no cycle check runs.
   - *2-cycles.* If the new key matches the first filter, the walk evaluates f(new)
     unbatched. If key(f(new)) = key(W) exactly, then {W, new} is a 2-cycle of the class
     walk. Both of its points have been tested, and the walk replaces new by 2·canon(M),
     where M is the element with the smaller key.
   - *3..64-cycles.* Unless a 2-cycle escape just happened, if the new key matches the
     anchor filter, the walk iterates f from the new point, unbatched, up to 64 times,
     testing each look-ahead point first.
     - A distinguished look-ahead point ends the trail there, since the walk would reach
       it next.
     - If instead the new point's class returns, every point of that cycle has been
       tested. The walk escapes from 2·canon(M), with M the cycle element of least key.
   - The anchor check also runs after a first-filter match that the exact check did not
     confirm.
   - So a walk in a cycle of length ℓ ≤ 64 leaves it within 64 + ℓ steps, unless the
     cycle contains a distinguished point, which ends the trail first.
6. **Distinguished points.** W is distinguished when bits 16..75 of key(W) are zero. The
   test is applied to every point that a step or an escape produces. At a distinguished
   point the trail ends: the record (key(W), a, b), with the trail's coins, goes to a
   hash table keyed by key(W), and the slot starts a new trail. A trail makes at most one
   record, so no trail can store a point twice. The run fails if more than DP_MAX records
   arrive.
7. **Resolution.** The first record whose key is already in the table ends the run. No
   other key match is used.
   - Both trails are replayed, unbatched, from their recorded coins up to their first
     distinguished point, with the rules of steps 4–6 and anchors every 64 of the trail's
     own steps. The replay tracks (a, b) mod r: an image sigma(R_i) = ±R_i has coefficients multiplied
     by (−1)^s, and canon multiplies them by (−1)^s, and an escape doubles the
     coefficients of canon(M). By Lemma R each replay reaches the recorded point.
   - The replays give the common canonical point as [a1]P + [b1]Q = [a2]P + [b2]Q. If
     b1 ≠ b2, output x = (a2 − a1)/(b1 − b2) mod r after checking [x]P = Q. If b1 = b2
     the run fails. By Lemma D this has probability 1/r.

8. **Work meter.** Every run is charged at most W_max = S_max·974.199265 + G
   instructions, with G the fixed global charges (initial starts, trail restarts, distinguished-point
   records, solve and verify, input validation, adder table).
   - Every walk step is charged C_COMMON = 845 instructions, which bounds the common path
     including its batch inversion, and each rare event is charged its worst case.
   - The meter aborts the run when the total metered rare work exceeds 64 times its
     expectation bound (Lemma M). That happens with probability at most 1/64, counted in delta.
   - The run also fails if more than DP_MAX records arrive or more than K trails are abandoned.

**Merged trails stay merged.** Outside fruitless cycles the next class depends only on the
current class. Escapes fire only after an exact check, and they exit from 2·canon(M),
which depends only on the cycle's set of classes, not on where the walk entered the cycle
or when it noticed. The 64-bit filters only decide when an exact check runs. So two trails
that reach the same class follow the same classes afterwards and reach the same first
distinguished class.

**Lemma R (replays reach the recorded point).** The first distinguished point of a trail,
and every escape before it, depend only on the trail's start point. They do not depend on
when the filters fire or on the anchor schedule.

*Proof.* Let the walk be at a point, with f deterministic on points.
- An escape happens only after an exact check shows that the point's class lies on a
  closed cycle C of the class map of length ℓ ≤ 65, and after every point of C has been
  tested and found not distinguished:
  - the 2-cycle rule tests both points of {W, new} before it escapes;
  - the anchor rule tests each look-ahead point before adding it to the cycle.
- So without the escape the walk would repeat C forever without a distinguished point,
  and the escape point 2·canon(M) depends only on C's set of classes.
- A look-ahead that meets a distinguished point ends the trail at the first distinguished
  point that f alone would reach next, which is where the walk would stop anyway. (A
  2-cycle cannot lie strictly inside a longer cycle, so no escape could intervene.)
- A false filter match either finds no closed cycle through the current class, and
  changes nothing, or finds a genuine one, and is handled as above.
- Every cycle of length ℓ ≤ 64 is caught. Within 64 steps of entering it the walk is
  anchored at one of its points, and it returns to that class within ℓ more steps. The
  anchor check runs there, because for ℓ ≥ 3 no 2-cycle escape can happen inside C, and
  the anchor check is not skipped by a false first-filter match.

So under any schedule that anchors each walk at least once every 64 of its steps, the
trail meets the same cycles, escapes to the same points and stops at the same first
distinguished point. Only the number of extra traversals differs. Every coefficient pair
the replay carries represents its current point. ∎

Replaying a trail of n < TRAIL_MAX steps takes fewer than 2^12·n unbatched operations, by
construction:
- each replay step uses at most 70: one evaluation, two for a 2-cycle check, 66 for an
  anchor check, and one escape;
- the replay takes at most 33n steps. It escapes the same cycles as the live trail, each
  of which took at least 2 of the live trail's steps, and spends at most 64 more steps in
  each before escaping.

**Lemma D (degenerate matches).** The first key match has b1 = b2 with probability exactly
1/r.

*Proof.* A trail records once, so the two records come from different trails t1 ≠ t2,
possibly of the same slot. Condition on all start points and adders. That fixes every
trail's points, and so which records match first.
- Given its start point W = [a + b·x]P, a trail's coin b is uniform on [0, r), because
  each b determines exactly one a.
- Coins of different trails are independent.
- Along a replay, b changes as b ← u·b + d, where u is a unit (negation classes {W, −W}, up to sign)
  or twice one, and u, d depend only on the path.

So the canonical b at the distinguished point is U_i·b_ti + D_i with U_i invertible mod r,
and P[b1 = b2] = 1/r. ∎

**Lemma M (metered rare work, Markov reserve).** Let X be the total metered rare work of a
run: every 2-cycle check at its worst-case cost, every longer-cycle check at its worst case,
and every exact-k fallback. Then E[X] ≤ S·E_RARE, where S is the number of walk steps and
E_RARE = (2U + 2^10)/A + 3(68U + 2^10)/A^2 + 2^−55·2^10 per walk step. This uses only
unconditional per-step probabilities: a step enters a 2-cycle with probability at most
1/A, a longer fruitless cycle with probability at most 3/A^2 (Lemma C), and an exact-k
fallback with probability at most 2^−55 (Lemma 1). Markov's inequality then gives
P[X > 64·S·E_RARE] ≤ 1/64. The meter enforces that threshold, so no run is charged more
than its expectation allows, and the probability of failure this adds is at most 1/64.
Its size per walk step is 64·E_RARE = 129.199265.

## 2. Success probability

Let T_steps = K·L. Then tau = T_steps^2/(2N) = 0.530300.

**H1 (heuristic, disclosed).** The class map [W] → [canon(W) + R_h([W])], with uniform
random adders, behaves like a random mapping on the N classes, apart from the r-adding
penalty. That penalty is a factor of about (1 − 1/(2A)) in the birthday exponent
(Bos–Kleinjung–Lenstra 2010, eq. (2); Teske 2001; Bernstein–Lange 2012). Under H1, the
probability that the first T_steps steps contain no useful class collision is at most
exp(−tau·(1 − eps)). Here eps = 4/A + 2·(3/A) + 2·(3·192/A^2) = 0.00015286. Its three
margins are as in the order-6 package:
- the adder penalty is taken as 4/A, eight times 1/(2A);
- fruitless 2-cycles enter at a rate of at most 1/A, with 3 non-fresh states each;
- short cycles enter at a rate of at most 1/A^2, with 192 wasted states each.

**H1′ (unconditional per-step form of H1).** Each walk step, averaged over all histories,
enters a fruitless cycle with probability bounded as in Lemma C, and the probability of a
distinguished point is theta. Only these unconditional probabilities are used, through Markov's
inequality for the metered rare work (Lemma M) and through Lemma C. No bound on the cost of a
step conditional on its history is used.

**Lemma C (fruitless-cycle rates, under H1′).** A fruitless cycle is a cycle whose formal
coefficient sum vanishes, so it closes whatever the adders are. At any step, the
probability that the walk enters one of length ℓ is at most:
- 1/(2A) for ℓ = 2;
- 3/A^2 in total over 3 ≤ ℓ ≤ 64;
- 2^−300 for each ℓ ≥ 65.

*Proof.* Let the steps of the cycle use adder indices i_1..i_ℓ with canonical units u_k.
- The adders' coefficients are independent and uniform, so the cycle closes for every
  choice of adders only if, for every index i that occurs, the unit-weighted sum of its
  occurrences vanishes mod r. Any other closure has probability at most ℓ/r; those are the
  self-collisions in the failure terms below.
- An index that occurs once cannot cancel. For an index that occurs m ≥ 2 times, given
  the other units the last one must take a single value, with probability at most
  1/2 ≤ 1/2.
- So the probability is at most E[2^−D·1{no index occurs once}], with D the number of
  distinct indices among ℓ independent uniform draws from A.

The cases:
- ℓ = 2: the indices agree (probability 1/A) and the units cancel (1/2).
- 3 ≤ ℓ ≤ 64: having no singleton needs D ≤ ⌊ℓ/2⌋, which has probability at most
  C(A, ⌊ℓ/2⌋)·(⌊ℓ/2⌋/A)^ℓ. With the factor 2^−D, ℓ = 3 contributes at most 1/(2A^2),
  ℓ = 4 at most 2/A^2, and all ℓ ≥ 5 together below 2^−40.
- 65 ≤ ℓ ≤ A/8: P[D ≤ ℓ/2] ≤ (e·ℓ/(2A))^(ℓ/2). This decreases on that range and is at
  most 2^−309 at ℓ = 65.
- ℓ > A/8: D only grows with ℓ. D < A/16 has probability below 2^−10000, and otherwise
  2^−D ≤ 2^−4096. ∎

Consequences:
- 2-cycles are entered at rate at most 1/A, and cycles of length 3..64 at rate at most
  3/A^2. These are the rates used in eps.
- Summed over all lengths 65..TRAIL_MAX, longer fruitless cycles are entered at a rate
  below 2^−230 per step. Over all S_max steps that is below 2^−100. Such a trap could
  only cause one abandoned trail; its steps are charged like any others.

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
- **Collision inside one trail.** There are at most DP_MAX + 2K trails, each with fewer
  than TRAIL_MAX ≈ 2^66 steps. The chance that any trail meets one of its own earlier
  classes, outside the cycles handled by step 5, is at most (DP_MAX + 2K)·TRAIL_MAX^2/N,
  computed in the accounting.
- **Abandoned trail.** Under H1 the distinguished-point indicators along a trail are
  Bernoulli(theta). A trail is abandoned only after more than CAP − 64 steps without one,
  or in a fruitless cycle longer than 64 (Lemma C). The chance is computed in the
  accounting. At most K trails may be abandoned before the run fails, and the expected number
  is far below that, so Markov's inequality bounds the overflow by the expectation.
- **Degenerate first match:** exactly 1/r (Lemma D).
- **Work meter abort (step 8):** at most 1/64 (Lemma M).
- **DP_MAX overflow:** the expected number of records is at most DP_MAX/2, so by Chernoff under
  H1 the overflow probability is at most exp(−DP_MAX/12), which is below 2^−60.
- **Exceptional additions:** a batch whose product is zero, or a point that a formula cannot
  handle, aborts the run. This has probability at most 2^−110 under H1.

Together these are below delta = 2^−20. Therefore:

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015286)) − 2^−20 = 0.411523 ≥ 0.39.**

The budget is 1.0299·√N steps, against 0.9943·√N for success exactly 0.39. The 0.0215
of slack in the model bound covers the sampling error of the toy evidence below.

**Evidence for H1, H1′ and the procedure.** `experiment.py` runs the procedure of
section 1.3 end to end on random prime-order toy curves
y^2 = x^3 − 3x + b (`"curve": "a3"`). It uses:
- the same class key, canonical rule, adder images and index (keys use Montgomery
  representatives, any fixed bijection being equivalent);
- both cycle rules and their escapes;
- trails from fresh coins that end at their first distinguished point, and the trail cap;
- resolution by replaying both trails from their coins.

It never tracks coefficients while walking and has no oracle collision table.
- The live walk anchors on a global round counter, as the batched kernel does. The replay
  anchors on the trail's own step count, so every replay also checks Lemma R.
- The filters compare keys masked to `filter_bits`, with the same distinguished-point-first
  order and look-ahead tests as the kernels. With `check_all_replays`, every recorded
  trail, not only the resolving pair, is replayed and checked against its record.
- Each trial runs L steps per slot plus the finishing phase at the stated tau. It counts
  as a success only if the first key match resolves to a verified x.
- Each trail start is charged as one of the slot's L steps, which is conservative.

<!-- toy-evidence -->
Offline runs (`evidence/trails/`, pooled in `summary.json`): 5,200 trials on 5 curves with
N = 2^28.6, 2^30.1, 2^32.9, 2^32.4, 2^34.7 classes, adder tables of 2^11 and 2^16,
16 slots and theta = 2^−5.
- Success of the whole procedure at tau = 0.5303: 2,135 of 5,200, 0.4106, Wilson 95% interval
  [0.397, 0.424]. The model 1 − exp(−tau) is 0.412.
- Every success was verified ([x]P = Q). Degenerate first matches: 0.
- Production 64-bit filters. Every recorded trail was replayed with a different anchor
  schedule: 8,356,258 replays, 0 that missed the recorded point.
- Lemma R stress (`evidence/trails-stress.json`): the same procedure with the filters
  cut to 8 bits over 5,200 trials, so false filter matches were frequent.
  - 8,551,297 trails were replayed with a different anchor schedule, and 0 missed the
    recorded point.
  - 338,932 distinguished points were first reached during a look-ahead.
  - Those runs' success rates are not used: each such look-ahead adds up to 65 points
    that the step budget does not count.
- Abandoned trails: 1. Longer fruitless cycles escaped: 13. The finishing phase added at most
  1.5% to the K·L steps.
- 2-cycle entry rate per step:
  - 2.33e-04, 2.32e-04, 2.30e-04, 2.32e-04 with A = 2^11 (1/(2A) = 2.44e-04);
  - 6.91e-06 with A = 2^16 (1/(2A) = 7.63e-06).
<!-- /toy-evidence -->

The declared initial input (`curve` "a3", `bits` 28, 80 trials, seed 2801, production
64-bit filters, every trail replayed) runs in a few seconds. It reports 28/80 = 0.350,
Wilson 95% interval [0.255, 0.459], against the model 0.412, with 0 degenerate matches, 0
distinguished points reached during a look-ahead and 0 replays, of 21,267, that missed
their recorded point. With 80 trials the standard error is about 0.05, so it is a
reproducibility check; the success claim rests on the pooled offline trials above.

Limits: the toy curves have N up to about 2^35, against 2^255 in production. There theta is 2^−5,
so start points and the finishing phase are a few percent of the steps; in production
they are below 2^−30 of them. H1 is an extrapolation in N, supported by the size-free
random-mapping argument and the cited literature. It is not proved for this walk.

## 3. Computation

**No reference operation.** Every instruction of a walk step is charged at 1/C with C = 1000, as
COST_MODEL.md's H + R/C with H = 0. No reference operation is charged.

**Charge per walk step.** Three parts, charged separately:
- The common batched path is charged C_COMMON = 845 instructions per walk step. This is the
  larger of the plain and anchor steps, each measured as the whole-run instruction total difference
  between runs with one more batch step on production inputs: 54,631,596 instructions for a plain
  batch and 55,369,126 for an anchor batch, each of 2^16 walks, batch inversion included. The anchor measurement also contains its rare
  events, so C_COMMON is conservative.
- The rare events, metered at worst case with U = 2^16 instructions per unbatched operation (measured
  40,871 Ir): a 2-cycle check costs 2U + 2^10, a longer cycle costs 68U + 2^10, an
  exact-k fallback costs 2^10.
- The Markov reserve 64·E_RARE = 129.199265 per walk step (Lemma M).

The charge per walk step is therefore 974.199265 instructions.

**Field arithmetic.** `measurement/p256_field.h` implements Montgomery multiplication and
squaring with R = 2^256 as GCC inline AArch64 assembly in general-purpose registers.
- Because p ≡ −1 mod 2^64, the Montgomery factor is 1, and each of the four reduction
  rounds is a handful of shift and add instructions with no multiplications. This is the
  approach of OpenSSL's ecp_nistz256 (Gueron–Krasnov), re-derived for this harness with
  shifted-register operands and fewer constants.
- Products are 16 `mul` + 16 `umulh`; outputs are fully reduced.
- `p256_sub.h` holds branchy subtractions that add p back only on a borrow, and a fused
  a − b − c.
- Standalone costs in the same toolchain, per call including about 5 instructions of loop
  overhead: mul 123, sqr 98, sub 26. Inside the kernel they are inlined.
- Correctness: the replica check below compares every final point of 2^16 walks after
  2^19 steps with Python big-integer arithmetic. The field routines
  were also checked on about 150,000 random and edge-case operands each.

Inline assembly is ordinary scalar AArch64 code. Callgrind counts each executed
instruction exactly as for compiled C, which is how the organizer measured blst's assembly.

**Measurement (H3).** `measurement/step.c` runs the trails, steps, cycle handling, distinguished
points, restarts and trail cap. Its common batched path was measured on production inputs as
whole-run instruction totals (Callgrind `Ir`, GCC 14.2.0 in the organizer's image,
`-O3 -fno-tree-vectorize -mgeneral-regs-only`), with the rare paths out of line so that each
path's cost is separable. The first batch step has no anchor and no rare event. The anchor step is 11.3 instructions per walk step above the plain step,
and its total contains its own rare events.

- Plain step (steps 1 and 2), per walk step: 833.6 instructions.
- Anchor step (steps 63 and 64), per walk step: 844.9 instructions.
- C_COMMON = 845 per walk step.

**Rare-path checks.** Distinguished points and false 64-bit filter matches do not occur in a
production-size sample, so small builds of the same source and replica exercise them
(`evidence/measurement/trail-paths.json`, with 2^6 adders and 2^12 walks):
- DP_BITS = 4: restarts at distinguished points, with 2-cycles;
- DP_BITS = 12 with a one-period cap: abandonments;
- DP_BITS = 6 with the filters cut to 8 bits: 2-cycles, longer cycles, and distinguished points first
  reached during a look-ahead.

In every build the replica matched every final point and every event count.

**Total.** S_max = K·(L + TRAIL_MAX) ≤ 2^127.542440 walk steps.

| Component | Bound (units) | log2 |
| --- | --- | ---: |
| walk steps | 974.199265·S_max / C | 127.504729 |
| adder table | 2^16 adders × 2 scalar mults × 768 unbatched ops × 2^16 / C | 32.62 |
| initial walk starts | K × 2 scalar mults × 768 × 2^16 / C | 48.62 |
| trail restarts | (DP_MAX + K) × 2 scalar mults × 768 × 2^16 / C | 85.16 |
| distinguished points | DP_MAX = 2^68.54 entries × 2^20 instructions / C | 78.58 |
| collision replay | 2 trails × TRAIL_MAX steps × 2^12 unbatched ops × (2^16 + 2^10) / C | 85.06 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 16.84 |
| input validation | 2^24 instructions | 14.03 |

Rounded up, **time_log2 = 127.504730**. Preprocessing (adders and initial walk starts) is included. All
processors' work is summed, and elapsed time plays no part.

## 4. Memory and advice

Peak memory is at most the sum of:
- 2·DP_MAX·128 bytes for the distinguished-point table (key, the trail's coins, hashing
  overhead);
- K·512 bytes of walk state (points, keys, coins, trail counter, batch buffers);
- 2^16·2·64 bytes of adder images;
- 16 MiB.

That is 2^76.542441 bytes. Nonuniform advice is 0 bytes.

## 5. Limitations

- H1 is the standard but unproved random-mapping heuristic. It is supported by the
  literature and toy-curve experiments only. The model bound is 0.4115, 0.0215 above
  0.39.
- Instruction counts come from one implementation and compiler. The sample is 2^18 steps
  on one production input set. The branch-free common path and the rare paths were measured, not
  bounded analytically.
- About 2^68.5 stored distinguished points and 2^127.3 work are needed. This is not a
  practical attack and does not weaken P-256.
- The toy experiment runs the whole procedure, including resolution, but only at toy
  scale: theta = 2^−5, 16 slots, and N far below production. Its trail starts and its
  finishing phase are a few percent of the steps there, and negligible in production.
- Distinguished points and false 64-bit filter matches never occur in a
  production-parameter measurement sample. The kernel's restart, abandonment and
  false-match paths were checked on small builds of the same source.
- The work meter's success loss, 1/64, is a Markov bound on unconditional rates (Lemma M).

## References

- Pollard (1978).
- van Oorschot–Wiener (J. Cryptology 1999).
- Teske (Math. Comp. 2001).
- Bos–Kleinjung–Lenstra (ANTS 2010).
- Bernstein–Lange–Schwabe (PKC 2011).
- Bernstein–Lange (ANTS 2012).
- Bos–Costello–Miele (PKC 2014).
- Gueron–Krasnov, *Fast prime field elliptic-curve cryptography with 256-bit primes*
  (J. Cryptographic Engineering 2015), the OpenSSL ecp_nistz256 field method.

## Reproduction

```sh
python3 -I submissions/p256-ecdlp/rho-sigma2/accounting.py
bash submissions/p256-ecdlp/rho-sigma2/measurement/reproduce.sh 1   # Docker, linux/arm64
echo '{"curve":"a3","bits":28,"trials":80,"seed":2801,"adder_bits":11,"dp_bits":5,"walks":16,"tau_e4":5303,"filter_bits":64,"check_all_replays":1}' | python3 -I submissions/p256-ecdlp/rho-sigma2/experiment.py
```

# p256-ecdlp: parallel Pollard rho on negation classes with a sigma walk

Claim: for uniform x in [0,r) and Q=[x]P on NIST P-256, the algorithm below returns x
with probability at least 0.39. It uses at most **2^127.370918** target-group-operation equivalents in total (score 63685459/500000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^76.542441 bytes. Preprocessing is 2^47.949581 units, and no nonuniform advice is
used. This is a classical cold attack on one target.

The success bound rests on the disclosed heuristic H1, the standard random-mapping model
for r-adding walks on negation classes, and on M1, the class-data premise of section 2.0. The table
premise M2 (the table from the fixed seed 1, treated as uniform and independent of the target) is stated in section 2.0. From M1, M2 and exact
algebra follow the escape, re-entry, abandonment and metered-work bounds (section 2). The cost rests on H3,
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
  `evidence/accounting.json`). It also checks the order of 2 modulo r used by D1 (section 2.1).
- `measurement/` holds the per-step instruction-count harness and the field code
  (`p256_field.h`, `p256_sub.h`), plus the Python replica `walk_ref.py`.
- `evidence/measurement/` holds the production totals, the 16-step replica check and the checks of the
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
- CAP = 64/theta. Every trail has at most TRAIL_MAX = CAP + 64 steps.
- DP_MAX = 2·theta·S_max + 2K, where S_max = K·(L + TRAIL_MAX).
- |Gamma| = 2: the class automorphisms are the signs only.

### 1.2 Classes, keys and the sigma walk

For W = (x, y) ≠ O:
- **Class key.** key(W) = M(x). It is constant on {W, −W} and determines the class.
- **Canonical point.** canon(W) = (x, (−1)^s·y), where s = M(y) mod 2. p is odd and y ≠ 0
  (no 2-torsion), so M(−y) = p − M(y) has the other parity, and canon picks one point of
  each class. This is the only canonical choice: there is no rotation (measurement/step.c).
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
   - take c_i, d_i as the i-th pair of the seed-1 stream of section 2.0, treated as uniform in [0,r) (premise M2);
   - set R_i = [c_i]P + [d_i]Q;
   - tabulate R_i and −R_i.
3. **Trails.** Each of the K walk slots runs a sequence of trails.
   - A trail starts from fresh coins a, b, taken in order from the same seed-1 stream after the table and treated as uniform in [0,r) (M2, section 2.0), with W = [a]P + [b]Q. The
     slot keeps the coins until the trail ends. Coefficients are not tracked while
     walking.
   - A trail ends at its first distinguished point (step 6), or is abandoned. A trail is
     abandoned when it crosses more than 2^60 anchor boundaries (step 5) without one. Its first
     anchor boundary comes within 64 steps of its start, so this happens after s + CAP steps for some
     s in [1, 64]; every trail has at most TRAIL_MAX = CAP + 64 steps. An addition that gives O also
     ends the trail, and a slot whose trail ends this way starts a new trail from fresh coins.
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
     - A distinguished look-ahead point s ends the trail at s, which is then the trail's record, since the walk would reach
       it next.
     - If instead the new point's class returns, every point of that cycle has been
       tested. The walk escapes to 2·canon(M), with M the cycle element of least key.
   - The anchor check also runs after a first-filter match that the exact check did not
     confirm.
   - So a walk in a cycle of length ℓ ≤ 64 leaves it within 64 + ℓ steps, unless the
     cycle contains a distinguished point, which ends the trail first.
   - The kernel keeps the last four escape keys of each trail and counts a repeat among them. This is a
     diagnostic: it changes nothing in the walk. Section 2.2 bounds the event it counts.
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
8. **Work meter.** Every run is charged at most W_max = S_max·887.904982 + G
   instructions, with G the fixed global charges (initial starts, trail restarts, distinguished-point
   records, solve and verify, input validation, adder table).
   - Every walk step is charged C_COMMON = 847 instructions, which bounds the common path
     including its batch inversion, and each rare event is charged its worst case.
   - The meter aborts the run when the total metered rare work exceeds 64 times its
     expectation bound (Lemma M', section 2.5). That happens with probability at most 1/64, counted in delta.
   - The run also fails if more than DP_MAX records arrive or more than K trails are abandoned.

**Escape.** A cycle's escape goes to 2·canon(M), where M is the cycle's least-key member. canon(M) is a point of M's class, so 2·canon(M) is a point of the class of 2M. The escape uses no table and no stored point: it is one point doubling, whose cost is charged at its forced-path maximum (Section 3.1). The escape point depends only on the cycle's members, as Lemma R requires. No freshness of M or of the escape point is assumed. Section 2.1 (D1) gives the escape point a new formal coefficient, which is what the re-entry and abandonment bounds use.

**Escape repeats.** The escape point of a cycle depends only on its least-key member, so it is the same every time the trail escapes from that cycle. A repeated escape point therefore means that the trail has re-entered a cycle it already escaped from. Section 2.2 bounds this event (D2 and the within-trail bound).

**Merged trails stay merged.** Outside fruitless cycles the next class depends only on the
current class. Escapes fire only after an exact check, and they exit to 2·canon(M),
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
- A look-ahead that meets a distinguished point ends the trail at the first distinguished point that f alone would reach next, which is where the walk would stop anyway. (A
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

Replaying a trail of n ≤ TRAIL_MAX live steps takes at most 2^13·n unbatched operations on every run, by
construction:
- each replay step uses at most 70 operations: one evaluation, two for a 2-cycle check, 66 for an
  anchor check, and one escape;
- the replay takes at most 65n steps. The live trail makes one walk step per live step and at most one
  escape per live step (Section 2.1). The replay follows the same classes (Lemma R), so it makes at most one
  walk step per live step, and at most one look-ahead of MAX_CYCLE = 64 steps per live step, since each live
  step has at most one anchor check. So each live step costs at most 65 replay steps, and 65 × 70 = 4550 ≤ 2^13
  operations.
The bound is deterministic, so the collision-replay term is charged in every run (REPLAY_OPS = 2^13).

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

**Metered rare work.** The rare-path calls and the within-trail coincidences are metered at their forced-path worst cases (Section 3.1). Their rates per walk step, and the Markov reserve over all processed steps (Lemma M'), are in Section 2.5.

## 2. Success probability

Let T_steps = K·L. Then tau = T_steps^2/(2N) = 0.530300.

### 2.0 Probability space and premises

Probabilities below are over the uniform target x and, for the coefficients, over M2: the start coefficients
(a_w, b_w), which are also the trail coins, and the adder coefficients (c_i, d_i), with R_i = [c_i]P + [d_i]Q.
Two premises about the random objects, M1 and M2, and the standard heuristic H1 are used, and they are stated
here so that nothing else enters.

- **M1 (class data).** The class-level values read by the walk are functions of the class alone: the
  class key, which fixes the adder index i = key mod 2^16, the 64-bit filter bits and the distinguished-point
  test, and the canonical sign s. Distinct classes carry independent uniform values, independent of x.
  P-256 has no canonical rotation and no exact-k fallback. canon uses only the parity of M(y)
  (measurement/step.c, `escape` and `select_adder`), so the fraction-valued class data of the other tracks
  do not occur, and no fallback rate is charged.
- **M2 (table, fixed seed).** The table and the walk coins come from one stream, `random.Random(1)`: seed 1 is a
  constant of the specification, chosen without reference to the target, and is the default of `walk_ref.py gen`
  and `reproduce.sh`. The stream is read in this order. First one draw `randrange(1, r)` is taken and discarded,
  because it is the target value of the test generator and the target here is the input Q. Then each point takes
  two draws, a = `randrange(r)` and b = `randrange(r)`. The first 2^16 pairs are the adders, R_i = [a]P + [b]Q with
  (c_i, d_i) = (a, b). The next 2^16 pairs are the walk starts (the start coefficients, which are also the trail
  coins), and the last 1024 pairs are spares. So the coefficients are a function of the seed alone, not of the
  target's value. M2 is the premise that these coefficients behave as independent uniform values in Z_r,
  independent of the uniform target x. The success probability is over the uniform target x and, for the
  coefficients, over M2. It is a premise, not a theorem: the seeded stream is treated as a uniform one, and the seed
  is not analysed. The measured inputs are generated this way: `walk_ref.py gen` with seeds 1 and 2 reproduces the
  stored input files byte for byte, and both match this procedure at sampled indices.
- **H1 (random mapping).** The standard model for the collision statistics, used only in 2.6.

M1 is the premise that carries the claim (section 2.7). The rest of this section uses it through the
conditional rates of 2.3 and the fresh-point statement of 2.4.

### 2.1 Formal structure: D1 and D2

Within one trail with start S, every walk point has the formal form
  W = 2^m g S + sum_i n_i R_i,
with g ∈ Γ = {±1}, integers n_i, and m the number of escapes so far in the trail. A step adds σ(R_i) = ±R_i and
does not change the S-coefficient. An escape replaces M by 2 canon(M) = 2gM, which doubles the S-coefficient and
the n_i. Its pair (a, b) is the same linear map applied to the table pairs of M2.

**Escape depth.** A step escapes at most once: a 2-cycle escape disables the anchor check, and an anchor
escape replaces the point that the step produced. So a trail of at most TRAIL_MAX steps has at most
K_ESC = TRAIL_MAX escapes. An escape point can start a 2-cycle at the next step, so the escapes of one trail
may share steps. The count is therefore per step, and it is the one used here.

**D1 (depth).** On the good event G of 2.2, if E is the (m+1)-th escape point of a trail, its S-coefficient
2^{m+1} g' differs from the S-coefficient 2^j g'' of every earlier point of the trail (j ≤ m).

*Proof.* Equality would give 2^k = g''/g' ∈ Γ for k = m+1−j ∈ [1, K_ESC]. Then 2^{2k} = 1, so ord_r(2)
divides 2k and ord_r(2) ≤ 2K_ESC. Hence it suffices that

  (C-ord)  ord_r(2) > |Γ| K_ESC = 2·TRAIL_MAX.

**Verification for the P-256 order.** Trial division finds r − 1 = 2^4 · 3 · 71 · 131 · 373 · 3407 · 17449 · 38189 · C2,
where C2 is a 188-bit composite with no prime factor below 2·10^6. C2 is not factored, and the check does not need it.
Let S = 2^4 · 3 · 71 · 131 · 373 · 3407 · 17449 · 38189 = 378059953041213003408, so r − 1 = S·C2, and let
U = 2^{C2} mod r. Then U^S = 2^{r−1} = 1, so ord_r(U) divides S. Also ord_r(2) is a multiple of ord_r(U), since
the order of a power divides the order of the base. Computing ord_r(U) from the factorisation of S gives
ord_r(U) = S/2 = 189029976520606501704: U^{S/2} = 1, and U^{(S/2)/p} ≠ 1 for every prime p dividing S/2.
Since S/2 ≈ 1.89·10^20 > 2·TRAIL_MAX ≈ 1.48·10^20, (C-ord) holds: ord_r(2) ≥ S/2 > 2·TRAIL_MAX.
`accounting.py` recomputes this block (`ord_2_lower_bound`) and asserts each step. The condition is therefore
verified for this r, and D1 holds on G with no assumption about r.

**D2 (re-entry is not structural).** A re-entry into an escaped cycle is a later point (depth ≥ m+1) in the
classes of an earlier cycle (depth ≤ m). By D1 it is never formally equal to an earlier point of the trail.
So it is either a nondegenerate linear equation in x or an accidental equality of pairs. Both are within-trail
coincidences, bounded in 2.2.

### 2.2 Degenerate pairs, nondegenerate coincidences, and the good event G

Let e and e' be two points of one trail whose formal expressions differ. Given the transcript (the indices, the
distinguished-point tests and the escape choices, which are functions of the class data and so independent of x
and of the table by M1 and M2), their pairs (a, b) are linear forms in the table coefficients (a_S, b_S, c_i, d_i)
with coefficients in Z_r. For each g ∈ Γ, Pr[(a,b)_e = g(a,b)_{e'}] ≤ 1/r by Schwartz–Zippel, since the
difference is a nonzero linear map to Z_r², whose kernel has density at most 1/r. Conditioning on the earlier
history removes at most a fraction |Γ|t/r of the coefficient space, so the conditional bound is at most 2|Γ|/r
per pair. Summed over the at most TRAIL_MAX²/2 pairs of a trail, this gives 2|Γ|(TRAIL_MAX²/2)/r per trail.
The table itself contributes at most A²|Γ|/r, over the pairs of adders with (c_i, d_i) = g(c_j, d_j).

**Lemma G (nondegenerate coincidences).** Let the next point of a trail have its class among its t earlier
classes, with no earlier coincidence. Then

  Pr[class among the t earlier classes | history] ≤ |Γ| t / (r − |Γ| t).

*Proof.* Without an earlier coincidence, the transcript is a function of the class data of distinct classes,
which are independent of x by M1 (and of the table by M2). Conditioned on the transcript, x is uniform on the
complement of at most |Γ|t roots: each earlier point gives at most one root per element of Γ, for the linear
equation its equality with the new point would impose. The new point coincides with an earlier one only at
one of those roots. ∎

**The good event G.** Let G be the event that no trail has (i) an accidental equality of two formally distinct
points, or (ii) a class repeat other than a formal repeat of length at most 64. A formal repeat of length at most 64
is a closure in the sense of 2.3: the walk returns to an earlier class with coefficients that agree for every
table. It is detected and escaped (Section 1.3 step 5, Lemma R). Summing the bounds above over all trails
started in a run that does not fail on its counts:

  Pr[not G] ≤ DELTA_LOOP := K_TRAILS · (TRAIL_MAX²/2) · |Γ| · (2/r + 1/(r − |Γ| TRAIL_MAX)) + A²|Γ|/r,

where K_TRAILS = DP_MAX + 2K counts the trails started: the K initial trails, at most DP_MAX trails started after
records, and at most K trails started after abandonments. By D2, every re-entry into an escaped cycle lies outside
G, so DELTA_LOOP bounds the re-entry probability. Lemma R and Lemma D
hold for every trail, on G or not.

### 2.3 Closure rates (Lemma C')

A *closure* of length ℓ is a formal repeat: after ℓ steps the walk is at ±W_0 for every table. Let i_1..i_ℓ be
the adder indices of the steps and σ_1..σ_ℓ ∈ {±1} their canonical signs, which are class data (M1).
- The "+" case W_ℓ = W_0 requires Σ_j σ_j R_{i_j} = 0. For each index i the signed count
  n_i = Σ_{j: i_j = i} σ_j must then vanish: a nonzero integer n_i gives a nontrivial linear
  equation on the table, which holds with probability at most 1/r by M2. Summing the signed
  counts over the indices gives Σ_j σ_j = 0, so **ℓ must be even**.
- The "−" case W_ℓ = −W_0 needs a nontrivial linear relation on the table, with probability at most 1/r. It is
  an accidental equality, bounded in 2.2.

Under M1 the closure probabilities are then:
- **ℓ = 2:** the indices agree (1/A) and the signs are opposite (1/2). The rate is exactly 1/(|Γ|A) = 1/(2A).
- **ℓ odd:** zero.
- **ℓ = 4:** either two distinct indices, each used twice with opposite signs in each pair, at most
  3A(A−1)/A⁴ · 1/4 = 3(A−1)/(4A³); or one index used four times with balanced signs, at most
  A/A⁴ · 3/8 = 3/(8A³). The total is below 0.76/A².
- **even ℓ = 6..64:** no index occurs once, so there are at most ℓ/2 distinct indices, and the standard count of
  that event applies. The sum over these lengths is at most 2^{-40}; the largest term, ℓ = 6, is about 2^{-45}.

So for 3 ≤ ℓ ≤ 64 the total closure probability per step is at most 0.76/A² + 2^{-40} ≤ 3/A². The constant c = 3 is
kept, which is conservative. Hence, conditional on the past and on G, a step closes with probability at most

  p_c := 1/(|Γ|A) + 3/A² + 2^{-63} = 2^{-17} + 3·2^{-32} + 2^{-63} ≈ 7.630093·10^{-6}.

The 2^{-63} term covers false 63-bit filter matches. A false match is rejected by the exact check and does not close,
so the term is conservative. The two filters each match falsely with probability 2^{-64}, which gives 2^{-63} per step.

**Long closures.** For even ℓ ≥ 65 a structural closure needs at most ℓ/2 distinct indices, so its probability is at
most (eℓ/(2A))^{ℓ/2}, by the same count. At ℓ = 65 this is below 2^{-309}, and it decreases over
65 ≤ ℓ ≤ A/8, with the longer lengths smaller still. So a step starts a closure of length 65 or more with
probability below 2^{-230}, and over all S_max steps

  DELTA_TAIL := S_max · 2^{-230} ≈ 2^{-102.5}.

### 2.4 Abandonment (Lemma T', first moment)

A trail is abandoned after m steps with CAP + 1 ≤ m ≤ n := TRAIL_MAX (Section 1.3 step 3). Let C_m be the number of
closure steps among its first m steps. A class is *fresh* when the trail has not produced it before; on G every non-closure step produces a fresh class.
Conditioned on the past and on G, such a step
whose key bits are uniform and independent of the past (M1). So it is a distinguished point with probability
θ = 2^{-60}, independently of the past. A closure wastes at most 128 steps before its repeats are detected, so at
least m − 128 C_m of the m steps are fresh.

*First moment.* By the chain rule, the k given steps all close with probability at most p_c^k, so

  Pr[C_m ≥ k] ≤ binom(m, k) p_c^k ≤ (e n p_c / k)^k.

*Bound.* Take c = 2^{14} and k = ⌊n/c⌋ = 2^{52}. Then e·c·p_c ≈ 0.340 ≤ 1/2, so the first term is at most 2^{-k}.
On the event C_m < k, at least CAP + 1 − 128k steps are fresh, and none of them is a distinguished point with
probability at most (1−θ)^{CAP+1−128k}. Hence

  Pr[abandon] ≤ (e n p_c / k)^k + (1−θ)^{CAP+1−128k} ≤ 2^{-k} + e^{−63.5 + O(2^{-54})} ≤ 2.65·10^{-28}.

*Expectation.* The expected number of abandoned trails is at most K_TRAILS · 2.65·10^{-28} ≈ 1.1·10^{-7}. The run
fails by abandonment only if more than K trails are abandoned, which by Markov has probability at most that
expectation divided by K:

  DELTA_ABANDON ≤ 2.65·10^{-17}.

### 2.5 Metered rare work (Lemma M', every processed step counts)

Every processed walk step counts, revisits included. The metered events, with their conditional rates per walk step
(on G, under M1 and M2) and worst cases, are:

| event | rate per walk step | worst case (instructions) |
| --- | ---: | ---: |
| 2-cycle check (true or false predecessor match) | 1/(2A) + 2^−63 | 83,518 |
| anchor check (true or false anchor match), with its look-ahead steps (at most 64) | 3/A² + 2^−63 | 2,789,687 |
| within-trail coincidence (2.2) | |Γ| TRAIL_MAX (2/r + 1/(r − |Γ| TRAIL_MAX)) | TRAIL_MAX (C_COMMON + RESERVE) |

The exact-k fallback is not an event type, because P-256 has no rotation and no fallback. The first two
worst cases are the forced-path maxima of Section 3.1. The third charges the remaining trail at the per-step
charge. That adds nothing to the common path, and it keeps every processed step covered.

Let E0 be the expected work of the first two rows, and A = TRAIL_MAX × (rate of the third row). The third row's
expected contribution is A (C_COMMON + KM·E_RARE), so the expected metered rare work per walk step is

  E_RARE = (E0 + A·C_COMMON) / (1 − KM·A),   KM = 64.

Numerically E0 = 0.6391403463 and E_RARE = 0.6391403463. The third row contributes about 2.5·10^{-34}. The
reserve per walk step is KM·E_RARE = 40.904982, and the charge per walk step is C_COMMON + 40.904982 = 887.904982.

*Markov.* E[X] ≤ S_max·E_RARE for the total metered rare work X of a run. Hence P[X > 64·S_max·E_RARE] ≤ 1/64. The
meter aborts when its cumulative metered work exceeds this bound. Cumulative work never decreases, so the abort
event is {X > 64·S_max·E_RARE}. This probability is counted in delta.

### 2.6 Success probability

**H1 margins.** Under H1 the first T_steps steps contain no useful class collision with probability at most
exp(−τ(1 − ε)), where ε = 4/A + 2·(3/A) + 2·(3·192/A²) = 0.00015286. The three margins are as follows.
- The adder penalty is taken as 4/A, eight times 1/(2A).
- Fruitless 2-cycles enter at rate at most 1/A (2.3 gives 1/(2A)), with 3 wasted states each.
- Closures of length 3..64 enter at rate at most 3/A², with 192 wasted states each. 2.4 bounds the waste of
  one closure by 128 steps, so 192 is conservative.

**Failure terms.** Each is counted below. The meter's 1/64 is subtracted separately because it is not below 2^{-20}.
- Within-trail coincidences, degenerate or not (2.2): DELTA_LOOP ≈ 6.1·10^{-17}.
- Closures of length 65 or more (2.3): DELTA_TAIL ≈ 1.4·10^{-31}.
- Abandoned trails (2.4): DELTA_ABANDON ≈ 2.6·10^{-17}.
- Degenerate first key match: exactly 1/r (Lemma D).
- Work meter abort (2.5): 1/64.
- DP_MAX overflow: the records are Bernoulli(θ) indicators on fresh steps, with mean at most DP_MAX/2 (M1), so by
  Chernoff the overflow probability is at most exp(−DP_MAX/12), below 2^{-60}.
- Exceptional additions (a batch whose product is zero, or a point that a formula cannot handle): at most 2^{-110}
  under H1, since these are collisions of x-coordinates.

Together:

  Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015286)) − 1/64 − 2^{-20} − DELTA_ABANDON − DELTA_LOOP − DELTA_TAIL = 0.395898 ≥ 0.39.

The model term is 0.411524. The budget is 1.0299·√N steps. Without the losses, success exactly 0.39 needs
0.9943·√N steps; with the 1/64 and 2^{-20} losses it needs 1.0201·√N. The slack above 0.39 after all losses is 0.0059.

**Evidence for H1, M1 and the procedure.** `experiment.py` runs the procedure of
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
64-bit filters, every trail replayed, escape "double"; `evidence/declared-run-double.json`) runs in a few seconds. It reports 28/80 = 0.350, Wilson 95% interval [0.255, 0.459], against the model 0.412, with 0 degenerate matches, 0
distinguished points reached during a look-ahead and 0 replays, of 21,267, that missed
their recorded point. With 80 trials the standard error is about 0.05, so it is a
reproducibility check; the success claim rests on the pooled offline trials above.

Limits: the toy curves have N up to about 2^35, against 2^255 in production. There theta is 2^−5,
so start points and the finishing phase are a few percent of the steps; in production
they are below 2^−30 of them. H1 is an extrapolation in N, supported by the size-free
random-mapping argument and the cited literature. It is not proved for this walk. M1 is not tested at
production theta, and the toy runs cannot test it there.

### 2.7 Premises, verified facts and what is not established

- **Premises.** M1 (class data) and H1 (collision statistics) are the heuristics. M2 is the fixed-seed table premise
  (2.0), a premise and not a theorem. H3 (Section 3) is a measurement assumption about instruction counts, not a
  probabilistic premise.
- **Derived from M1, M2 and exact algebra.** D1 (with C-ord verified for this r), D2, the degenerate-pair bound and
  DELTA_LOOP (Lemma G and Section 2.2), the closure rates (Lemma C'), the abandonment bound (Lemma T'), the metered
  work bound (Lemma M') and DELTA_TAIL.
- **Verified here.** C-ord for the P-256 order (2.1); the forced-path maxima of Section 3.1 on the bounding builds; the
  accounting and its checks, reproduced by `accounting.py`.
- **Not established.** That real keys behave as M1 at production size. The toy runs use theta = 2^{-5} and cannot
  test the production regime. The formal claims (D1 and the closure rates of 2.3) can be tested on small builds of
  the replica with this group order.

## 3. Computation

**No reference operation.** Every instruction of a walk step is charged at 1/C with C = 1000, as
COST_MODEL.md's H + R/C with H = 0. No reference operation is charged.

**Charge per walk step.** Two parts, charged separately:
- The common batched path is charged C_COMMON = 847 instructions per walk step: the larger of the plain and
  anchor steps, each measured as the whole-run instruction total difference between runs with one more batch step
  on production inputs, per walk, with the rare events of the anchor batch removed (see the measurement below), plus
  the bound on the one data-dependent branch of the common path.
- The rare events, each charged its forced-path worst case (Section 3.1): a 2-cycle check at most
  83,518 instructions and an anchor check, with its look-ahead steps, at most 2,789,687. Their expected
  counts are covered by the Markov reserve 64·E_RARE = 40.904982 per walk step (Lemma M', Section 2.5).

The charge per walk step is therefore 887.904982 instructions.

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

**3.1 Worst cases of the rare-path operations.** Each rare event is charged a forced-path maximum, not a sampled
average. The primitives are f_single (adder selection and the unbatched addition), add_single, escape (its least-key
search, its canonical point and its doubling), same_class, and point_is_dp. For each primitive
the data-dependent branches were enumerated from the source and the disassembly. Each conditional correction in a
field routine that is not branch-free was replaced, for the measurement only, by a branch-free routine that costs at
least as much on every path (a bounding build; equality with the original is checked on random and edge inputs).
The remaining control variants (the equal-x and equal-y tests of the addition, the parity of the least member, the
length of a cycle, and the memcmp outcomes of the library comparison) were forced
by chosen inputs, and the maximum over the variants was taken. P-256 has no rotation search and no exact-k
fallback (measurement/step.c), so there is no such variant. Each count was checked to be linear in the number of
calls, and identical for a second seed of inputs. Measured maxima (Callgrind Ir per call, organizer image and flags,
`measurement/forced/summary.json`):

| primitive | worst case (Ir) |
| --- | ---: |
| f_single (adder selection and addition) | 41,190 |
| add_single | 41,185 |
| same_class | 32 |
| point_is_dp | 10 |
| escape of a 2-cycle (canon, then doubling) | 41,262 |
| escape of a 65-cycle (canon, then doubling) | 43,089 |

The repeat counter (a diagnostic, Section 1) costs at most 46 instructions per call and is inside the overhead.

A 2-cycle check costs at most f_single + same_class + escape(2) + point_is_dp plus a control overhead of at most
1,024 instructions, and an anchor check costs at most f_single + escape(65) plus, for each of its at most 64
look-ahead steps, f_single + point_is_dp + same_class, with the same overhead per step. The overhead bound is far
above the measured control code of these functions, and the sums are above the longest path through the code.
Scalar multiplications and replays use f_single's worst case as their unit, an upper bound for an addition or a
doubling.

**Measurement (H3, a measurement assumption).** `measurement/step.c` runs the trails, steps, cycle handling, distinguished
points, restarts and trail cap. Its common batched path was measured on production inputs as
whole-run instruction totals (Callgrind `Ir`, GCC 14.2.0 in the organizer's image,
`-O3 -fno-tree-vectorize -mgeneral-regs-only`), with the rare paths out of line so that each
path's cost is separable, on two input sets (seeds 1 and 2; `evidence/measurement/seed1.json`). The first batch steps have no anchor and no rare event. The anchor step contains the
rare events that fall in step 64 (1 two-cycle escape on the first input set, 3 on the second, 82,205 instructions each on average); they are
subtracted here, because rare work is charged separately.

- Plain step (steps 1 and 2), per walk step: 833.612 instructions on the first input set and 833.612 on the second.
- Anchor step (steps 63 and 64) without rare events, per walk step: 843.612 and 843.612 instructions.
- The one data-dependent branch of the common path that stays on it is the first test of the
  distinguished-point check (key >> 16 = 0). Its other edge evaluates the second test and adds 3 instructions per
  walk step. The other data-dependent branches (the two filter matches, the trail cap, a distinguished point) lead
  to rare events, which are charged separately, and the field arithmetic has no branch. So the plain step is at
  most 836.612 and the anchor step at most 846.612 instructions per walk step.
- C_COMMON = 847 per walk step.

**Rare-path checks.** Distinguished points and false 64-bit filter matches do not occur in a
production-size sample, so small builds of the same source and replica exercise them
(`evidence/measurement/trail-paths.json`, with 2^6 adders (2^4 in the last build) and 2^12 walks):
- DP_BITS = 4: restarts at distinguished points, with 2-cycles;
- DP_BITS = 12 with a one-period cap: abandonments;
- DP_BITS = 6 with the filters cut to 8 bits: 2-cycles, longer cycles, and distinguished points first
  reached during a look-ahead;
- 2^4 adders, DP_BITS = 8, 800 steps: 94,624 2-cycles and 2,959 longer cycles escaped, and no escape point repeated
  (the repeat counter stays 0).

In every build the replica matched every final point, every event count (the repeat counter included), and a digest of every distinguished record, look-ahead records included (`evidence/measurement/trail-paths.json`).

**Total.** S_max = K·(L + TRAIL_MAX) ≤ 2^{127.542440} walk steps.

| Component | Bound (units) | log2 |
| --- | ---: | ---: |
| walk steps | 887.904982·S_max / C | 127.370918 |
| adder table | 2^16 adders × 2 scalar mults × 768 unbatched ops × 2^16 / C | 31.95 |
| initial walk starts | K × 2 scalar mults × 768 × 2^16 / C | 47.95 |
| trail restarts | (DP_MAX + K) × 2 scalar mults × 768 × 2^16 / C | 84.49 |
| distinguished points | DP_MAX = 2^68.54 entries × 2^20 instructions / C | 78.58 |
| collision replay | 2 trails × TRAIL_MAX steps × 2^13 unbatched ops × (2^16 + 2^10) / C | 85.40 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 16.29 |
| input validation | 2^24 instructions | 14.03 |

Rounded up, **time_log2 = 127.370918**. Preprocessing (adders and initial walk starts) is included. All
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

- **M1 carries the success bound.** It is the class-data premise of section 2.0, not proved for this walk, and the
  toy runs test it only at theta = 2^{-5}. H1 (collision statistics) is the standard but unproved random-mapping
  model. The model term is 0.4115, 0.0215 above 0.39. After all losses the bound is 0.395898.
- **M2 is a premise.** The table comes from the fixed specification seed 1 (Section 2.0). Treating its coefficients as uniform and
  independent of the target is not proved; the seeded generator is assumed to behave as a uniform one.
- **H3 is a measurement assumption.** Instruction counts come from one implementation and compiler. The common path
  was measured on two production input sets (Section 3), and its worst case is bounded there. The rare operations
  are charged forced-path maxima (Section 3.1). The bounding builds are checked for dominance on the forced paths,
  not proved for every input.
- **Union bounds, not measured rates.** DELTA_LOOP, DELTA_TAIL and DELTA_ABANDON are upper bounds from the stated counts.
  They are not rates measured in any run.
- **Escape depth.** The escape depth is bounded per step (K_ESC = TRAIL_MAX), not per two steps, so the D1 check needs
  ord_r(2) > 2·TRAIL_MAX. That is verified for this r (Section 2.1).
- About 2^68.5 stored distinguished points and 2^127.4 work are needed. This is not a
  practical attack and does not weaken P-256.
- The toy experiment runs the whole procedure, including resolution, but only at toy
  scale: theta = 2^−5, 16 slots, and N far below production. Its trail starts and its
  finishing phase are a few percent of the steps there, and negligible in production.
- Distinguished points and false 64-bit filter matches never occur in a
  production-parameter measurement sample. The kernel's restart, abandonment and
  false-match paths were checked on small builds of the same source.
- The work meter's success loss, 1/64, is a Markov bound on the expected metered work (Lemma M', Section 2.5).

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

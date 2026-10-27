# Predecessor specification and evidence

This file reproduces, for audit, the predecessor sections the current proof relies on (1.2-1.3, 2.0-2.4, 2.8, 3-4); 1.1, 2.5-2.7, 5, References and Reproduction, and the superseded common-path paragraphs of Section 3, are omitted for the package size limit and remain in the predecessor submission named in proof.md. Its score, tau, L and reserve figures are historical; the current claim is proof.md and submission.json. No prior acceptance is assumed.

# p256-ecdlp: parallel Pollard rho on negation classes with a sigma walk

Claim: for uniform x in [0,r) and Q=[x]P on NIST P-256, the algorithm below returns x
with probability at least 0.39. It uses at most **2^127.370918** target-group-operation equivalents in total (score 63685459/500000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^76.542441 bytes. Preprocessing is 2^49.961206 units, and no nonuniform advice is
used. This is a classical cold attack on one target.

The success bound rests on the disclosed heuristic H1, the standard random-mapping model
for r-adding walks on negation classes, on M1, the equivariant random-orientation premise of section 2.0, with
its orientation law L-ori, on the global filtration of section 2.0 (domination by i.i.d. data of classes that no walk has visited), which replaces the window premise W, and on the first-merge argument of section 2.8 (Lemmas P and M). The algorithm draws its adder coefficients and every trail's coins from the operating system's
randomness, so every probability is over the uniform target x and the algorithm's coins (section 2.0). From M1, L-ori
and exact algebra follow the escape, re-entry, abandonment, cut-off and metered-work bounds (section 2). The cost rests on
H3, average instruction counts measured with the organizer's toolchain, and every run is charged by a work meter with a
Markov reserve. **There is no reference operation on this track:** every instruction, including all P-256 field
arithmetic and point additions, is charged at 1/C with the borrowed normalisation C = 1000. This is a carefully accounted
generic attack. It is not a practical break.

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

### 1.2 Classes, keys and the sigma walk

For W = (x, y) ≠ O:
- **Class key.** key(W) = M(x). It is constant on {W, −W} and determines the class.
- **Canonical point.** canon(C) is the point of the class C with M(y) even. It is a function of C alone, and
  canon(W) = (x, (−1)^s·y) with s = M(y) mod 2 picks it out of any W in C. Since M(−y) = p − M(y) has the other
  parity (p is odd, y ≠ 0), exactly one of W and −W has s = 0. There is no rotation (measurement/step.c).
- **Sign and orientation.** A point W of a class C is ±canon(C). Its sign s(W) is 0 exactly when W = canon(C). The
  walk stores actual points, so the sign of an arriving point is the orientation of its class relative to canon(C),
  which is the premise of section 2.0 (M1 and L-ori). The sign is not a class label: s(−W) = 1 − s(W).
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
   - draw c_i, d_i uniformly from [0,r) from the operating system's randomness (for example
     `secrets.randbelow(r)`), fresh in every run and independent of Q (section 2.0), each draw under the sampling caps of section 8.2;
   - set R_i = [c_i]P + [d_i]Q;
   - tabulate R_i and −R_i.
3. **Trails.** Each of the K walk slots runs a sequence of trails.
   - A trail starts from fresh coins a, b drawn uniformly from [0,r) from the operating system's randomness, at
     every trail start including every restart (section 2.0), under the sampling caps of section 8.2, with W = [a]P + [b]Q, two scalar multiplications
     that are charged. The slot keeps the coins until the trail ends. Coefficients are not tracked while
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
   point the trail ends: the record (key(W), a, b), with the trail's coins, is inserted into the open-addressing table of section 8.2 (probe cap L = 2^10), and the slot starts a new trail. A trail makes at most one
   record, so no trail can store a point twice. The run fails if more than DP_MAX records
   arrive.
7. **Resolution.** The first record whose key is already in the table ends the run. No
   other key match is used.
   - Both trails are replayed, unbatched, from their recorded coins up to their first
     distinguished point, with the rules of steps 4–6 and anchors every 64 of the trail's
     own steps. The replay tracks (a, b) mod r: an image sigma(R_i) = ±R_i has coefficients multiplied
     by (−1)^s, and canon multiplies them by (−1)^s, and an escape doubles the
     coefficients of canon(M). By Lemma R each replay reaches the recorded point.
   - The replay of a trail with n live steps is cut off after 65n steps (REPLAY_FACTOR = 65 in the core). A cut-off
     replay does not reach its recorded key and fails the run. Its probability is counted in DELTA_CUT (section 2.2).
   - The replays give the common canonical point as [a1]P + [b1]Q = [a2]P + [b2]Q. If
     b1 ≠ b2, output x = (a2 − a1)/(b1 − b2) mod r after checking [x]P = Q. If b1 = b2
     the run fails. By Lemma D this has probability 1/r.
8. **Work meter.** Every run is charged at most W_max = S_max·887.904982 + G
   instructions, with G the fixed global charges (initial starts, trail restarts, distinguished-point
   records, solve and verify, input validation, adder table).
   - Every walk step is charged C_COMMON = 847 instructions, which bounds the common path
     including its batch inversion, and each rare event is charged its worst case.
   - The meter aborts the run when the total metered rare work exceeds 64 times its
     expectation bound (Lemma M'', section 2.5). That happens with probability at most 1/64, counted in delta.
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
The bound is deterministic, so the collision-replay term is charged in every run (REPLAY_OPS = 2^13). The cut-off of
step 7 at 65n steps is never reached by a replay that follows the live trail's classes. Reaching it means that the replay
diverged, which section 2.2 bounds.

**Lemma D (degenerate matches).** The first key match has b1 = b2 with probability exactly
1/r, over the algorithm's coins.

*Proof.* A trail records once, so the two records come from different trails t1 ≠ t2,
possibly of the same slot. Condition on all start points and all adder points. Since f is deterministic on points, this
fixes every trail's points, and so which records match first, and the transcript of each trail (its adder indices, its
distinguished-point tests and its escape points).
- Given its start point W = [a + b·x]P, the coins (a, b) of a trail lie on a line in Z_r² that is parametrised by b, because
  each b determines exactly one a. So b is uniform on [0, r) given the start point.
- The coins of different trails are independent (section 2.0), so the coins b1, b2 of t1 and t2 are conditionally
  independent and uniform given all start points.
- Along a replay, b changes as b ← u·b + d, where u is a unit (negation classes {W, −W}, up to sign)
  or twice one, and u, d depend only on the path, which the conditioning fixes.

So the canonical b at the distinguished point is U_i·b_ti + D_i with U_i invertible mod r, and P[b1 = b2] = 1/r.
Since the coins are sampled exactly as stated in section 2.0, this holds exactly over the coins. ∎

**Metered rare work.** The rare-path calls and the within-trail coincidences are metered at their forced-path worst cases (Section 3.1). Their rates per processed step, and the Markov reserve over all processed steps (Lemma M'', section 2.5), are in Section 2.5.

## 2. Success probability

Let T_steps = K·L. Then tau = T_steps^2/(2N) = 0.530300.

### 2.0 Probability space, randomness and premises

Every probability below is over the uniform target x and the algorithm's coins. The coins are the adder coefficients
(c_i, d_i), the start coefficients (a_w, b_w), and every trail's coins (a, b), all drawn uniformly from [0, r) from
the operating system's randomness, fresh in every run and at every trail start (including every restart), and
independent of the target Q. The algorithm draws them; they are not a model. Each trail's start W = [a]P + [b]Q
costs two scalar multiplications, which the ledger charges.

- **M1 (equivariant random orientation).** For each class C:
  - a uniform key K_C, which fixes the adder index i = key mod 2^16, the 64-bit filter bits and the distinguished-point
    test;
  - the canonical representative canon(C), a function of C alone (section 1.2);
  - a uniform orientation g_C ∈ Γ = {±1}, saying where the walk's representative of C sits relative to canon(C). For a
    point W of C, the sign s(W) (and, for other tracks, the rotation) is determined by W's relation to canon(C) and by g_C.

  All of these are independent across classes and independent of x and of the coins. P-256 has no canonical rotation
  and no exact-k fallback (measurement/step.c, `escape` and `select_adder`), so no coordinate-fraction value is used and no
  fallback rate is charged.
- **L-ori (orientation law).** When the walk arrives at a fresh class C, g_C is uniform given the history, so the
  arriving point's sign is uniform on {0, 1} given the history (a measured premise, §2.2, not a verified law). At a revisited class the arriving point's sign is
  determined by the history. Every use below cites this law: the 2-cycle rate (Lemma E1) and the source orientation of
  the closure windows (2.3).
- **Global filtration (not a premise).** All walks form one process. A class is visited once a walk's point lies in it, and its data are drawn once, at its first visit. Conditionally on the full history, classes not yet visited receive uniform, independent keys (M1) and orientations (L-ori), so the windows of (c) and the long-closure bound are dominated by i.i.d. windows of unvisited classes (section 2.0).
- **H1 (random mapping).** The standard model for the collision statistics, used only in 2.6.

M1 and L-ori are the premises that carry the claim (section 2.7). The randomness of the coins is not a heuristic, and Lemma D
holds exactly over it (section 1.3). The measured inputs use the fixed seed `random.Random(1)` as the reproducible instance
only (section 2.5). The algorithm does not use it.

**Global filtration and coupling.** Fix any history F in which the walk's next classes are unvisited (no merge has occurred, section 2.8). Conditionally on F, the keys and orientations of the classes visited afterwards are uniform and independent (M1, L-ori). Equivalently, the new classes take their data from one sequence ξ_1, ξ_2, … of i.i.d. uniform elements, independent of F. Under this coupling, on the event that the next classes are new, the walk uses ξ_1, ξ_2, … in order. Hence, for every event E of the next d−1 steps' data and every such history F,

  Pr[E and the next d−1 classes are new | F] ≤ Pr_iid[E],

where the right side does not depend on F. No data are resampled per trail, and nothing is assumed about visited classes.

### 2.1 Formal structure: D1 and D2

Within one trail with start S, every walk point has the formal form
  W = 2^m g S + sum_i n_i R_i,
with g ∈ Γ = {±1}, integers n_i, and m the number of escapes so far in the trail. A step adds σ(R_i) = ±R_i and
does not change the S-coefficient. An escape replaces M by 2 canon(M) = 2gM, which doubles the S-coefficient and
the n_i. Its pair (a, b) is the same linear map applied to the start and adder pairs of the coins.

**Escape depth.** A step escapes at most once: a 2-cycle escape disables the anchor check, and an anchor
escape replaces the point that the step produced. So a trail of at most TRAIL_MAX steps has at most
K_ESC = TRAIL_MAX escapes. Every escape is the resolution of an entry (2.3), so escapes ≤ entries ≤ processed steps ≤
TRAIL_MAX. The smaller bound TRAIL_MAX/2 would need the cycles of different escapes to use disjoint steps. That fails
when an escape point starts the next 2-cycle, because the next cycle then shares the step that the escape replaced.
The bound K_ESC = TRAIL_MAX is therefore used. The probabilistic depth tail of the entry framework (escapes bounded by
the entry count, with a first-moment term) is not needed here, because the certificate below covers K_ESC = TRAIL_MAX
deterministically.

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
verified for this r, and D1 holds on G with no assumption about r. No coin enters this check.

**D2 (re-entry is not structural).** A re-entry into an escaped cycle is a later point (depth ≥ m+1) in the
classes of an earlier cycle (depth ≤ m). By D1 it is never formally equal to an earlier point of the trail.
So it is either a nondegenerate linear equation in x or an accidental equality of pairs. Both are within-trail
coincidences, bounded in 2.2.

### 2.2 Degenerate pairs, nondegenerate coincidences, cut-off, and the good event G

Let e and e' be two points of one trail whose formal expressions differ. Given the transcript (the adder indices, the
distinguished-point tests and the escape choices), their pairs (a, b) are linear forms in the start and adder coefficients
(a_S, b_S, c_i, d_i) with coefficients in Z_r. For each g ∈ Γ, Pr[(a,b)_e = g(a,b)_{e'}] ≤ 1/r over the coins, by
Schwartz–Zippel: the difference is a nonzero linear map to Z_r², whose kernel has density at most 1/r. The transcript is
fixed before the coins in this step, as in the standard argument; the step that passes from the transcript to the coin
values is part of M1. Conditioning on the earlier history removes at most a fraction |Γ| C(t,2)/r of the coefficient space
(each earlier pair and Γ element is a hyperplane of density 1/r). Since |Γ| TRAIL_MAX²/(2r) ≤ 1/2, the conditional bound is at
most 2|Γ|/r per pair. Summed over the at most TRAIL_MAX²/2 pairs of a trail, this gives 2|Γ|(TRAIL_MAX²/2)/r per trail.
The table itself contributes at most A²|Γ|/r, over the pairs of adders with (c_i, d_i) = g(c_j, d_j).

**Lemma G (nondegenerate coincidences).** Let the next point of a trail have its class among its t earlier
classes, with no earlier coincidence. Then

  Pr[class among the t earlier classes | history] ≤ |Γ| t / (r − |Γ| t(t−1)/2) ≤ |Γ| t / (r − |Γ| TRAIL_MAX²/2).

*Proof.* Without an earlier coincidence, the transcript is a function of the class data of distinct classes, which
M1 takes as independent of x. Conditioned on the transcript, x is uniform on the complement of the roots of the earlier
pairs' linear equations. There are at most |Γ|·C(t,2) such roots, one per Γ element per earlier pair, so at most
|Γ| t(t−1)/2 values of x are excluded. The new point coincides with an earlier one only at one of its own at most |Γ| t
roots. This holds for every coin outcome without an earlier coincidence, and so over the coins and x together. ∎
(An earlier version used the denominator r − |Γ| t, which counted only the points and not the pairs of earlier points.)

**The sign transition law, the equal-index forcing and the abab closure law.** Let P'_k = canon(C_k) + R_{i_k}, let τ_k = ±1 be the sign of P'_k (τ_k = −1 when M(y) of P'_k is odd), and let u_k be the sign of the walk point W_k, so that W_k = u_k·canon(C_k). The step gives f(W_k) = u_k·P'_k, and canon(C_{k+1}) = τ_k·P'_k, so

  u_{k+1} = u_k·τ_k.

This sign transition law is derived from the definitions of canon and sigma; the replica check found no failure of it.

*Equal-index forcing (derived).* If i_{k+1} = i_k and τ_k = −1, then P'_{k+1} = −canon(C_k), so C_{k+2} = C_k (a 2-cycle) and τ_{k+1} = −1. Hence the pair (τ_k, τ_{k+1}) = (−1, +1) does not occur when i_{k+1} = i_k. If i_{k+1} = i_k and τ_k = +1, then P'_{k+1} = canon(C_k) + 2R_{i_k} and τ_{k+1} is not forced. The equal-index pairs with τ_k = −1 are revisits: they are the 2-cycle entries of Lemma E1(b), and they belong to the entry accounting, not to the fresh-window bounds.

*abab closure law (derived).* For an index pattern (i, j, i, j) with i ≠ j, the four steps close the cycle, C_{k+4} = C_k, if and only if τ_k ≠ τ_{k+1} and τ_{k+1} ≠ τ_{k+2}. This is a fact about the formal structure, not an excluded pattern: under the measured fresh-arrival law its conditional probability given the index pattern is 1/4. On the replica (production escape rule, A = 4 and A = 8) the law predicted every closure: 679 of 679 at A = 4 (2,807 abab windows) and 248 of 248 at A = 8 (986 windows), with no closure outside the rule.

*Measured premise L-ori (fresh arrivals).* When the walk arrives at a fresh class, τ is uniform given the history. This is a measured premise, not a verified law. On the replica the arrival sign at fresh classes has share 0.5003 over 352,580 arrivals (z ≈ 0.4), from seven runs of 60,000 steps at A = 4 and A = 8 (`evidence/sign-law.json`, key `fresh_arrival_runs`). It is used for the 2-cycle entry rate 1/(2A), for Lemma E1 and for the window bounds of 2.3. The non-forced share is 0.5008 at A = 64 (n = 387,972, z ≈ 1 against 1/2) and 0.4987 at A = 8 (n = 77,464), both consistent with 1/2 (keys `tau_transition_runs` and `summary`). These are small-A replica measurements, not production data, and `evidence/sign_law.py` regenerates them from `measurement/walk_ref.py`.

*Observation in the stored production counts (not evidence for L-ori).* The seed-1 counts (the predecessor's `evidence/measurement/seed1.json`, omitted here for the size limit) increase by 1 and 3 at step 64 for the two input sets, against 0.5 expected each (pooled 4 against 1.0, one-sided Poisson p ≈ 0.02). The cumulative count at step 64 (65 against 64.0) is consistent with the law, so the excess is post hoc (`evidence/sign-law.json`, key `production_steps`). The stored records carry no signs or adder indices, so this does not test L-ori or W.

**The good event G.** Let G be the event that no trail has (i) an accidental equality of two formally distinct
points, or (ii) a class repeat other than a formal repeat of length at most 64. A formal repeat of length at most 64
is a closure in the sense of 2.3: the walk returns to an earlier class with coefficients that agree for every
table. It is detected and escaped (Section 1.3 step 5, Lemma R). Summing the bounds above over all trails
started in a run that does not fail on its counts:

  Pr[not G] ≤ DELTA_LOOP := K_TRAILS · (TRAIL_MAX²/2) · |Γ| · (2/r + 1/(r − |Γ| TRAIL_MAX²/2)) + A²|Γ|/r,

where K_TRAILS = DP_MAX + 2K counts the trails started: the K initial trails, at most DP_MAX trails started after
records, and at most K trails started after abandonments. By D2, every re-entry into an escaped cycle lies outside
G, so DELTA_LOOP bounds the re-entry probability. Lemma R and Lemma D hold for every trail, on G or not.

**Cut-off and false filter matches (DELTA_CUT).** A replay is cut off after 65n steps (section 1.3, step 7). A
divergence between a live trail and its replay needs a cycle that one anchor schedule detects and the other misses. A
false 64-bit filter match alone cannot cause this. An escape happens only after an exact key check, so a false match
changes when the check runs and not which escape happens (Lemma R). A genuine cycle of length at most 64 is detected by
every schedule that anchors each walk every 64 steps, so a divergence needs a cycle of length 65 or more that one
schedule detects and the other misses. Such a cycle is a structural repeat of length at least 65, whose union over all fresh
starts is bounded by DELTA_LONG (section 2.3). A within-trail coincidence is bounded by DELTA_LOOP. Hence a cut-off replay,
which fails the run, has probability at most

  DELTA_CUT := DELTA_LOOP + DELTA_LONG ≈ 6.1·10^{−17}.

No per-trail term is charged for false filter matches, because a false match alone cannot cause divergence.

**Lemma L (re-entry).** By D2, every re-entry into an escaped cycle is a within-trail coincidence, so its probability over
the run is at most DELTA_LOOP, and the event is contained in the failure term DELTA_CUT.

### 2.3 Entries (Lemmas E1 and E2)

**Definitions.** A step t is an *entry* if its source point lies on a fresh class (a class no walk has visited before, section 2.0) and its successor lies on a class already visited in the trail. A step whose source is on a
visited class is *persistent*. Only entries can open a closed cycle. A *window of length ℓ* ending at a source is the
ℓ steps before it; its *signed count vector* v_ℓ ∈ Z^A has entries Σ σ over the window's steps with index i.

**The formal structure of an entry.** Steps do not change the S-coefficient, so a formal repeat of a visited point has
the same S-coefficient and the sign +; it holds for every table iff the vector v_ℓ + σ_S e_{i_S} vanishes, where
i_S and σ_S are the index and orientation of the source's own step. Hence, for each ℓ ∈ {1, ..., 64}, a window ℓ gives an
entry target (i_S, σ_S) = t_ℓ exactly when −σ e_i = v_ℓ is a unit vector, and no target otherwise. On G every revisit is a
formal repeat of length at most 64, so an entry is one of these windows.

**Lemma E1 (entry probability).** Let F_{t−1} be the history before step t, and let the source of step t be on a fresh class.
Then, on G:
(a) Pr[entry at t | F_{t−1}] = w_t/(|Γ| A) ≤ P' := 64/(|Γ| A), where w_t is the number of distinct targets t_ℓ. The
    index i_S is uniform over A given the history (K_C is uniform and independent of F_{t−1}), and σ_S is uniform by
    L-ori (a measured premise, §2.2), so each target has probability 1/(|Γ| A).
(b) The 2-cycle entry (ℓ = 1) has conditional probability 1/(|Γ| A) = 1/(2A), given the measured L-ori law (§2.2). The predecessor P of the source is
    S − σ' R_{i'}, so the successor equals P exactly when i_S = i' and σ_S = −σ'. The "−P" case is impossible, because the
    S-coefficient of the successor equals that of S and would have to equal −(that of P), which is not zero mod r.
(c) The longer entries (ℓ = 3..64) have, in expectation over the history, rate at most c/A² with c = 3. Two facts about a
    window v_ℓ are used. (i) The coordinate sum of v_ℓ is the sum of its ℓ units, so v_ℓ is a unit vector only for odd ℓ, for
    every window and with no premise. (ii) Domination: for a window of fresh-source steps, the indices and units are dominated by
    i.i.d. uniform data (section 2.0), which gives P_3 = 3/(2A) − 3/(4A²), P_5 ≤ 15/(4A²) and
    P_ℓ = O(A^{−(ℓ−1)/2}) for ℓ ≥ 7. Hence, for these windows, E[w'_t] ≤ 3/(2A) + 4/A², and the longer-entry rate is
    E[w'_t]/(2A) ≤ 3/(4A²) + 2/A³ ≤ 3/A².

Hence the entry probability is bounded by p_c := 1/(2A) + 3/A² + 2^{−63} in expectation over the history, where the 2^{−63}
term is the false 63-bit filter match, which the exact check rejects. The conditional bound of (a), P' = 64/(2A), holds
for every history on G and is what the abandonment tail (2.4) uses. Only fresh-source windows can produce entries (derived). On G a formal repeat cannot span an escape, because each escape point has a new depth (D1); and within a segment the first revisit is itself an entry, after which the walk stays on visited classes until its escape. So a window that contains a persistent step is not the window of any entry, and the windows counted in (c) and in the long-closure bound are fresh-source windows only.

**Coupling.** The conditional bound of (a) uses only the source's own key and orientation, which are uniform given the history (M1, L-ori). The expectations in (c) and in the long-closure bound use the domination of section 2.0, which holds conditionally on every history of unvisited classes. No window premise is used.

**Lemma E2 (persistent steps).** On G, each entry is followed by at most 128 persistent steps before its escape,
deterministically. A 2-cycle entry is detected at the entry step, so it has no persistent step. For a longer entry, the
walk enters the cycle at the entry step, the anchor is set at a cycle member within 64 steps (the next anchor step), and the
walk returns to it within ℓ ≤ 64 further steps. So the count is at most 64 + ℓ ≤ 128. Consequently

  processed steps ≤ (fresh-source steps) + 128 × (entries).

Persistent steps are ordinary processed steps, charged the common cost C within S_max (section 3).

**Long closures (length 65 or more; Lemma F).** On G there is no formal repeat of length 65 or more (G(ii)). A closure of length d+1 ≥ 65 closes a fresh source C to a point at distance d ≥ 64. Escape counts agree on every trail (no escape count can exceed TRAIL_MAX, and ord_r(2) > 2·TRAIL_MAX, section 2.1), so the net adder vector of the d steps is a single adder, and every index but one occurs at least twice with zero unit sum. Let R_s = max(63, 2s − 2). Given any prefix, the d−1 fresh steps are dominated by i.i.d. data (section 2.0), so for every start the conditional probability of such a closure is at most

  Q_FAR := Σ_{s<A} A·C(A,s)·|Γ|^{1−s}·(s/A)^{R_s}·A/(A−s) + TRAIL_MAX·A·|Γ|^{1−A}.

Numerically log2 Q_FAR = −311.5 (the largest term is s = 33). The union bound over the at most S_max fresh starts gives

  DELTA_LONG := S_max·Q_FAR ≈ 2^{−184} ≈ 4·10^{−56}.

This is the long-closure term of DELTA_CUT (section 2.2).

### 2.4 Abandonment (Lemma T'', factorial moments over entries)

A trail is abandoned after m processed steps with CAP + 1 ≤ m ≤ n := TRAIL_MAX (Section 1.3 step 3). Let E_m be the number
of entries among its first m steps, and let F_m be its fresh-source steps. By E2, F_m ≥ m − 128 E_m. A fresh-source step
has key bits that are uniform and independent of the past (M1), so it is a distinguished point with probability
θ = 2^{−60}, independently of the past. A persistent step is never a distinguished point, because its class has
been visited, and a distinguished class would have ended the trail earlier.

*Factorial moments.* By the chain rule, Pr[entries at steps t_1 < … < t_k] ≤ P'^k, because the conditional bound of
Lemma E1(a) holds at every step for every history on G. Summing over all k-sets of steps,

  E[binom(E_m, k)] ≤ binom(n, k) P'^k ≤ (e n P' / k)^k.

*Bound.* Take c = 368 and k = ⌊n/c⌋. Then e·c·P' = 0.4886 ≤ 1/2, so the first term is at most 0.4886^k. On E_m < k, at least
CAP + 1 − 128k fresh-source steps occur, and none of them is a distinguished point with probability at most
(1−θ)^{CAP+1−128k}. Hence

  Pr[abandon, G] ≤ (e n P'/k)^k + (1−θ)^{CAP+1−128k} ≤ 0.4886^k + e^{−64(1 − 128/368) + O(2^{−54})} ≈ 7.5·10^{−19}.

*Expectation and failure.* Let X be the number of abandoned trails, with X ≤ K_TRAILS. The run fails by abandonment only
if X > K. On G, the expected number of abandoned trails is at most K_TRAILS × 7.46·10^{−19}, so Markov gives

  Pr[X > K, G] ≤ E[X 1_G] / K ≤ DELTA_ABANDON := K_TRAILS · 7.46·10^{−19} / K ≈ 7.5·10^{−8}.

Hence the run-failure bound for abandonment is Pr[not G] + E[X 1_G]/K ≤ DELTA_CUT + DELTA_ABANDON. The complement of G is
the event of 2.2 and appears once, in DELTA_CUT. The conditional bound P' = 64/(2A) is weaker than p_c, so this term is
larger than the value 2.6·10^{−17} that the sharper expectation would give. It is still far below 2^{−20}.

### 2.8 Merges: Lemmas P and M

A *merge* is a visit by a walk to a class visited earlier by a walk of another trail (start and escape points included); a step into a class of the walk's own trail is an entry (section 2.3). Let τ be the first batch step with a merge. A *root segment* is a trail's path up to its first merge; by section 2.0 the fresh-source bounds of sections 2.3 to 2.5 are exact on root segments.

**Lemma P (continuation).** On the good event G, the classes a walk visits from a class C up to its first distinguished class depend only on C. *Proof (outline).* The class, the key (hence the adder index and the distinguished test) and the escape rule depend only on the class (section 1.2), and the sign of the walk's point does not enter the class sequence. Every history reaches the same cycle set (Lemma R) and escapes to the same 2·canon(M). ∎

**Lemma M (the first merge resolves).** Let the first merge put walk j onto a class C of trail i, with step counts c_i and c_j there. Both trails record the first distinguished class D of the common continuation, and the second record gives the key match, unless the trail with the larger count has its first CAP = 2^66 steps free of distinguished points. Given the history before τ, that has probability at most (1−θ)^{CAP/2}, so

  DELTA_MERGE := 2·tau·(1−θ)^{CAP/2} ≤ 2·tau·e^{−32} ≈ 1.34·10^{−14},

where tau = T_steps²/(2N) = 0.5303 is the birthday count of section 2.6. *Proof (outline).* Under M1 on the global filtration each class is distinguished with probability θ = 2^{−60}, independently; the first CAP steps visit at least CAP/2 distinct classes except with probability at most (2e·p')^{CAP/2} and CAP·T_steps/N; a step merges with probability at most (visited count)/N under H1, which sums to tau. The factor 2 covers the two trails. ∎

*Edge cases.* A merge into an abandoned trail is the non-resolving case of Lemma M. A merge into the middle of an escape sequence continues with the escape of its cycle (Lemma P). The same class with the other sign is the same class, and Lemma D covers the coefficients. A restart draws fresh coins, so its start merges only into another trail. *Merged steps* take ordinary steps on another trail's path, at most TRAIL_MAX of them, and are charged in the merged row of section 2.5.

## 3. Computation

**No reference operation.** Every instruction of a walk step is charged at 1/C with C = 1000, as
COST_MODEL.md's H + R/C with H = 0. No reference operation is charged.

[Omitted here: the predecessor's common-path charge (C_COMMON = 847), its coin-uniformity remark and its field-arithmetic description. They concern the superseded batched kernel; the current kernel, its charge and its measurement are in proof.md. The field routines `measurement/p256_field.h` are unchanged and still used by the C paths.]

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

**Measurement (H3, a measurement assumption).** [The predecessor's common-path measurement of its kernel is omitted; see proof.md for the current one. The rare paths are out of line, so each path's cost is separable, and are priced by the forced maxima above.]

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

**8.2 Sampling caps and the distinguished-point table.** Every scalar is drawn by rejection sampling from [0, 2^256) with at most D = 4 draws, and every start point and adder point is attempted at most D = 4 times (retried only when it is O); the next draw or attempt aborts the run. Pr[abort] ≤ N_SCALARS·p^D + (A + K_TRAILS)·r^{−D} =: DELTA_CAP with p = 1 − r/2^256 ≈ 2^{−32} and N_SCALARS = 2(A + K_TRAILS), numerically about 2.5·10^{−18}.

The table has C = 2·DP_MAX slots of 128 bytes, open addressing with linear probing, index (key(W) >> 76) mod C (uniform up to a bias of about C/2^180), and load at most 1/2. A lookup or insert that examines L = 2^10 occupied slots aborts the run. *Probe tail:* if an operation examines L occupied slots from h, the maximal occupied run [s, s+ℓ) containing h has ℓ ≥ L and holds ℓ records hashed into the window [s, s+ℓ), since slot s−1 has always been empty. Under M1 the window count is binomial with mean at most ℓ/2, so Chernoff gives Pr ≤ e^{−cℓ}, c = ln 2 − 1/2. A union bound over the DP_MAX operations, the C window starts and ℓ ≥ L gives DELTA_PROBE ≈ 2.7·10^{−44}.

*Charge.* A table operation is charged at L probes times PROBE_COST = 2^11 instruction-equivalents, so DP_COST = 2^21 (raised from 2^20); the per-probe allowance covers the slot fetch (memory access), the comparison and the control, and is an allowance, not a measurement. A draw is charged DRAW_COST = 2^16. A start point is charged at its worst case, START = D(2·SCALAR_MUL + 2D·2^16) + 2^14, with SCALAR_MUL = 768·41,190; the adder table is charged ADDERS·START and the initial starts K·START. DELTA_CAP and DELTA_PROBE are added to delta.

**Total.** S_max = K·(L + TRAIL_MAX) ≤ 2^{127.542440} walk steps.

| Component | Bound (units) | log2 |
| --- | ---: | ---: |
| walk steps | 887.904982·S_max / C | 127.370918 |
| adder table | 2^16 × START / C (START of section 8.2) | 33.96 |
| initial walk starts | K × START / C | 49.96 |
| trail restarts | (DP_MAX + K) × START / C | 86.50 |
| distinguished points | DP_MAX = 2^68.54 entries × 2^21 instructions / C | 79.58 |
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


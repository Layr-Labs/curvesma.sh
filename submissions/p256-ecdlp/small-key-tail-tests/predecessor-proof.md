# Predecessor specification and evidence

This file reproduces, for audit, the predecessor sections the current proof relies on (1.2-1.3, 2.0-2.4, 2.8, 3-4) and, last, later revisions' sections, moved from proof.md for the file size limit (their section names are proof.md's); 1.1, 2.5-2.7, 5, References and Reproduction, and Section 3's superseded common-path paragraphs, are omitted for the package size limit and remain in the predecessor submission named in proof.md. Its score, tau, L and reserve figures are historical; the current claim is proof.md and submission.json.

# p256-ecdlp: parallel Pollard rho on negation classes with a sigma walk

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
   - draw c_i, d_i uniformly from [0,r) from the operating system's randomness (e.g.
     `secrets.randbelow(r)`), fresh every run and independent of Q (section 2.0), each draw under the sampling caps of section 8.2;
   - set R_i = [c_i]P + [d_i]Q;
   - tabulate R_i and −R_i.
3. **Trails.** Each of the K walk slots runs a sequence of trails.
   - A trail starts from fresh coins a, b drawn uniformly from [0,r) from the operating system's randomness, at
     every trail start and restart (section 2.0), under the sampling caps of section 8.2, with W = [a]P + [b]Q (two scalar multiplications, charged). The slot keeps the coins until the trail ends; coefficients are not tracked while
     walking.
   - A trail ends at its first distinguished point (step 6), or is abandoned. A trail is
     abandoned when it crosses more than 2^60 anchor boundaries (step 5) without one. Its first
     anchor boundary comes within 64 steps of its start, so this happens after s + CAP steps for some
     s in [1, 64]; every trail has at most TRAIL_MAX = CAP + 64 steps. An addition giving O also
     ends the trail, and a slot whose trail ends this way starts a new trail from fresh coins.
   - Each slot takes L steps (the main phase), then finishes its current trail and stops.
4. **Step.** W ← f(W), as one batched incomplete affine addition: Montgomery's
   simultaneous inversion over the 2^16 walks of a batch (3M), then lambda = dy·inv (M),
   lambda^2 (S) and y3 = lambda·(x1 − x3) (M), with one Fermat inversion per batch.
   - If the batch product is zero, the affected walk uses the unbatched complete formula.
5. **Fruitless cycles.** Each walk keeps two filters, both the low 64 bits of a key: the
   key of the point before its current point W, and an anchor, its key at the last batch
   step divisible by 64. All walks of a batch step together, so each walk is anchored once
   every 64 of its own steps. Because the filters keep only 64 bits they can match
   falsely. They only decide when an exact check runs.
   - *When the anchor is compared.* The anchor is written in the rounds t = 0 mod 64 and in
     those rounds only ("anchor rounds"), so the value compared is always the anchor set at the
     most recent such round. It is compared in the round it is written (before the overwrite,
     "anchor round") and in every later round whose offset 33 <= t mod 64 <= 63 ("plain rounds").
     Rounds with 1 <= t mod 64 <= 32 run a shorter tail that compares no anchor ("skip rounds");
     they still run every other test in this step. So after an anchor round a walk is compared
     against that anchor at offsets 33..64, where offset 64 is the next anchor round. Every walk
     of a batch shares t, so this schedule is the same for all walks and for the whole run.
   - *Order.* The distinguished-point test of step 6 comes first, on every point the walk
     produces: the new point, every look-ahead point and every escape point. A
     distinguished new point ends the trail, and no cycle check runs.
   - *2-cycles.* [Omitted: the two-steps-back filter and its 2-cycle escape; a two-cycle is the anchor
     rule's case ℓ = 2, "Checkpoint-only cycle detection" below.]
   - *3..64-cycles.* Unless a 2-cycle escape just happened, if the new key matches the
     anchor filter, the walk iterates f from the new point, unbatched, up to 64 times,
     testing each look-ahead point first.
     - A distinguished look-ahead point s ends the trail at s, then the trail's record, since the walk would reach
       it next.
     - If instead the new point's class returns, every point of that cycle has been
       tested, and the walk escapes to 2·canon(M), M the cycle element of least key.
   - The anchor check also runs after a first-filter match that the exact check did not
     confirm.
   - So a walk in a cycle of length ℓ ≤ 64 leaves it within 64 + 64 = 128 steps, unless the
     cycle contains a distinguished point, which ends the trail first: it is anchored on the
     cycle within 64 of its steps, and the first compare on the return is at an offset in
     33..64 after that round.
   - The kernel keeps the last four escape keys of each trail and counts a repeat among them. This is a
     diagnostic: it changes nothing in the walk. Section 2.2 bounds the event it counts.
6. **Distinguished points.** W is distinguished when bits 16..75 of key(W) are zero.
   At a distinguished point the trail ends and its record (key(W), a, b), with the trail's coins, is inserted into the open-addressing table of section 8.2 (probe cap L = 2^10), and the slot starts a new trail. A trail makes at most one
   record, so no trail can store a point twice. The run fails if more than DP_MAX records
   arrive.
7. **Resolution.** The first record whose key is already in the table ends the run. No
   other key match is used.
   - Both trails are replayed, unbatched, from their recorded coins up to their first
     distinguished point, with the rules of steps 4–6 and anchors every 64 of the trail's
     own steps. A replay compares the new key with the anchor at EVERY step, not on the live
     three-arm schedule: it is one walk, not a batch, so it has no skip rounds. That is a
     superset of the live compares and never a subset. The replay tracks (a, b) mod r: an image sigma(R_i) = ±R_i has coefficients multiplied
     by (−1)^s, canon multiplies them by (−1)^s, and an escape doubles the
     coefficients of canon(M). By Lemma R each replay reaches the recorded point.
   - The replay of a trail with n live steps is cut off after 65n steps (REPLAY_FACTOR = 65 in the core). A cut-off
     replay does not reach its recorded key and fails the run. Its probability is counted in DELTA_CUT (section 2.2).
   - The replays give the common canonical point as [a1]P + [b1]Q = [a2]P + [b2]Q. If
     b1 ≠ b2, output x = (a2 − a1)/(b1 − b2) mod r after checking [x]P = Q. If b1 = b2
     the run fails. By Lemma D this has probability 1/r.
8. **Work meter.** Every run is charged at most W_max = S_max·887.904982 + G
   instructions, with G the fixed global charges.
   - [Omitted: the superseded per-step charge (C_COMMON = 847), meter abort threshold and 1/64
     term; see proof.md, "Blinded corrections and their meter" and "Residual reserve".]
   - The run fails if more than DP_MAX records arrive or more than K trails are abandoned.

**Escape.** A cycle's escape goes to 2·canon(M), M the cycle's least-key member: canon(M) is a point of M's class, so 2·canon(M) lies in the class of 2M. The escape uses no table and no stored point: it is one point doubling, charged at its forced-path maximum (Section 3.1). The escape point depends only on the cycle's members, as Lemma R requires, and no freshness of M or of the escape point is assumed. Section 2.1 (D1) gives it a new formal coefficient, which the re-entry and abandonment bounds use.

**Escape repeats.** The escape point depends only on the cycle's least-key member, so it is the same at every escape from that cycle: a repeat means the trail re-entered a cycle it already escaped from. Section 2.2 bounds this event (D2; the within-trail bound).

**Lemma R (replays reach the recorded point).** The first distinguished point of a trail,
and every escape before it, depend only on the trail's start point, not on when the filters fire or on the anchor schedule.

*Proof.* Let the walk be at a point, with f deterministic on points.
- An escape happens only after an exact check shows that the point's class lies on a closed cycle C of the class map of length ℓ ≤ 65, and after every point of C has been tested and found not distinguished: the 2-cycle rule tests both points of {W, new}, and the anchor rule tests each look-ahead point before adding it to the cycle.
- So without the escape the walk would repeat C forever without a distinguished point, and 2·canon(M) depends only on C's set of classes.
- A look-ahead that meets a distinguished point ends the trail at the first distinguished point that f alone would reach next, which is where the walk would stop anyway. (A 2-cycle cannot lie strictly inside a longer cycle, so no escape could intervene.)
- A false filter match either finds no closed cycle through the current class, and changes nothing, or finds a genuine one, and is handled as above.
- Every cycle of length ℓ ≤ 64 is caught: within 64 steps of entering it the walk is anchored at one of its points, and returns to that class at an offset in 33..64 after the anchored round, where the anchor check runs (offset 64 is the next anchored round, compared before the anchor is overwritten). It runs there for ℓ ≥ 3 because no 2-cycle escape can happen inside C, and it is not skipped by a false first-filter match.

So under the schedule section 7 states -- every walk anchored every 64 of its steps, and every offset 33..64 after the anchored round compared -- the
trail meets the same cycles, escapes to the same points and stops at the same first
distinguished point. Only the number of extra traversals differs. Every coefficient pair
the replay carries represents its current point. ∎

Replaying a trail of n ≤ TRAIL_MAX live steps takes at most 2^13·n unbatched operations on every run: each replay step uses at most 70 operations (one evaluation, two for a 2-cycle check, 66 for an anchor check, one escape), and the replay takes at most 65n steps. The live trail makes one walk step per live step and at most one escape per live step (Section 2.1), and the replay follows the same classes (Lemma R), so each live step costs at most 65 replay steps and 65 × 70 = 4550 ≤ 2^13 operations. The bound is deterministic, so the collision-replay term is charged in every run (REPLAY_OPS = 2^13); the cut-off of step 7 at 65n steps is reached only by a replay that diverged, which section 2.2 bounds.

**Lemma D (degenerate matches).** The first key match has b1 = b2 with probability exactly
1/r, over the algorithm's coins.

*Proof.* A trail records once, so the two records come from different trails t1 ≠ t2, possibly of the same slot. Condition on all start points and all adder points: since f is deterministic on points, this fixes every trail's points, and so which records match first, and each trail's transcript (its adder indices, distinguished-point tests and escape points).
- Given its start point W = [a + b·x]P, a trail's coins (a, b) lie on a line in Z_r² parametrised by b, because each b determines exactly one a. So b is uniform on [0, r) given the start point.
- The coins of different trails are independent (section 2.0), so b1 and b2 are conditionally independent and uniform given all start points.
- Along a replay, b changes as b ← u·b + d, with u a unit (negation classes {W, −W}, up to sign) or twice one, and u, d fixed by the path, which the conditioning fixes.

So the canonical b at the distinguished point is U_i·b_ti + D_i with U_i invertible mod r, and P[b1 = b2] = 1/r. Since the coins are sampled exactly as stated in section 2.0, this holds exactly over the coins. ∎

## 2. Success probability

Let T_steps = K·L. Then tau = T_steps^2/(2N) = 0.530300.

### 2.0 Probability space, randomness and premises

Every probability below is over the uniform target x and the algorithm's coins. The coins are the adder coefficients
(c_i, d_i), the start coefficients (a_w, b_w), and every trail's coins (a, b), all drawn uniformly from [0, r) from
the operating system's randomness, fresh in every run and at every trail start (restarts included), and
independent of the target Q. The algorithm draws them; they are not a model.

- **M1 (equivariant random orientation).** For each class C:
  - a uniform key K_C, which fixes the adder index i = key mod 2^16, the 64-bit filter bits and the distinguished-point
    test;
  - the canonical representative canon(C), a function of C alone (section 1.2);
  - a uniform orientation g_C ∈ Γ = {±1}, saying where the walk's representative of C sits relative to canon(C). For a
    point W of C, the sign s(W) (and, for other tracks, the rotation) is determined by W's relation to canon(C) and by g_C.

  All of these are independent across classes and of x and the coins. P-256 has no canonical rotation and no exact-k
  fallback (measurement/step.c, `escape` and `select_adder`), so no coordinate-fraction value is used and no fallback
  rate is charged.
- **L-ori (orientation law).** When the walk arrives at a fresh class C, g_C is uniform given the history, so the
  arriving point's sign is uniform on {0, 1} given the history (a measured premise, §2.2, not a verified law); at a
  revisited class the history determines it. Every use below cites this law: the 2-cycle rate (Lemma E1) and the
  source orientation of the closure windows (2.3).
- **Global filtration (not a premise).** All walks form one process; a class is visited once a walk's point lies in
  it, and its data are drawn once, at its first visit. Conditionally on the full history, classes not yet visited
  receive uniform, independent keys (M1) and orientations (L-ori).
- **H1 (random mapping).** The standard model for the collision statistics, used only in 2.6.

M1 and L-ori carry the claim (section 2.7). The coins' randomness is not a heuristic, and Lemma D
holds exactly over it (section 1.3). The measured inputs use the fixed seed `random.Random(1)` as the reproducible instance
only (section 2.5); the algorithm does not use it.

**Global filtration and coupling.** Fix a history F in which the walk's next classes are unvisited (no merge, section 2.8). Conditionally on F, the keys and orientations of the classes visited afterwards are uniform and independent (M1, L-ori): the new classes take their data from one sequence ξ_1, ξ_2, … of i.i.d. uniform elements independent of F, and on the event that the next classes are new the walk uses ξ_1, ξ_2, … in order. Hence, for every event E of the next d−1 steps' data and every such F,

  Pr[E and the next d−1 classes are new | F] ≤ Pr_iid[E],

where the right side does not depend on F. No data are resampled per trail, and nothing is assumed about visited classes.

### 2.1 Formal structure: D1 and D2

Within one trail with start S, every walk point has the formal form
  W = 2^m g S + sum_i n_i R_i,
with g ∈ Γ = {±1}, integers n_i, and m the number of escapes so far in the trail. A step adds σ(R_i) = ±R_i and
does not change the S-coefficient. An escape replaces M by 2 canon(M) = 2gM, which doubles the S-coefficient and
the n_i. Its pair (a, b) is the same linear map applied to the start and adder pairs of the coins.

**Escape depth.** A step escapes at most once: a 2-cycle escape disables the anchor check, and an anchor escape replaces the point that the step produced. So a trail of at most TRAIL_MAX steps has at most K_ESC = TRAIL_MAX escapes, and since every escape is the resolution of an entry (2.3), escapes ≤ entries ≤ processed steps ≤ TRAIL_MAX. The smaller bound TRAIL_MAX/2 would need the cycles of different escapes to use disjoint steps; that fails when an escape point starts the next 2-cycle, which then shares the step the escape replaced. K_ESC = TRAIL_MAX is therefore used: the probabilistic depth tail of the entry framework is not needed, because the certificate below covers it deterministically.

**D1 (depth).** On the good event G of 2.2, if E is the (m+1)-th escape point of a trail, its S-coefficient
2^{m+1} g' differs from the S-coefficient 2^j g'' of every earlier point of the trail (j ≤ m).

*Proof.* Equality would give 2^k = g''/g' ∈ Γ for k = m+1−j ∈ [1, K_ESC], so 2^{2k} = 1, ord_r(2) divides 2k and
ord_r(2) ≤ 2K_ESC. Hence it suffices that

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

**D2 (re-entry is not structural).** A re-entry into an escaped cycle is a later point (depth ≥ m+1) in the classes
of an earlier cycle (depth ≤ m). By D1 it is never formally equal to an earlier point of the trail, so it is either a
nondegenerate linear equation in x or an accidental equality of pairs; both are within-trail coincidences, bounded in
2.2.

### 2.2 Degenerate pairs, nondegenerate coincidences, cut-off, and the good event G

Let e and e' be two points of one trail whose formal expressions differ. Given the transcript, their pairs (a, b) are linear forms in the start and adder coefficients
(a_S, b_S, c_i, d_i) over Z_r. For each g ∈ Γ, Pr[(a,b)_e = g(a,b)_{e'}] ≤ 1/r over the coins by Schwartz–Zippel: the
difference is a nonzero linear map to Z_r² with kernel density at most 1/r. The transcript is fixed before the coins in
this step; passing from the transcript to the coin values is part of M1. Conditioning on the earlier history removes at
most a fraction |Γ| C(t,2)/r of the coefficient space (each earlier pair and Γ element is a hyperplane of density 1/r).
Since |Γ| TRAIL_MAX²/(2r) ≤ 1/2, the conditional bound is at most 2|Γ|/r per pair; summed over the at most TRAIL_MAX²/2
pairs of a trail that is 2|Γ|(TRAIL_MAX²/2)/r per trail. The table contributes at most A²|Γ|/r, over the pairs of adders
with (c_i, d_i) = g(c_j, d_j).

**Lemma G (nondegenerate coincidences).** Let the next point of a trail have its class among its t earlier
classes, with no earlier coincidence. Then

  Pr[class among the t earlier classes | history] ≤ |Γ| t / (r − |Γ| t(t−1)/2) ≤ |Γ| t / (r − |Γ| TRAIL_MAX²/2).

*Proof.* Without an earlier coincidence, the transcript is a function of the class data of distinct classes, which
M1 takes as independent of x. Conditioned on the transcript, x is uniform on the complement of the roots of the earlier
pairs' linear equations. There are at most |Γ|·C(t,2) such roots, one per Γ element per earlier pair, so at most
|Γ| t(t−1)/2 values of x are excluded. The new point coincides with an earlier one only at one of its own at most |Γ| t
roots. This holds for every coin outcome without an earlier coincidence, and so over the coins and x together. ∎

**The sign transition law, the equal-index forcing and the abab closure law.** Let P'_k = canon(C_k) + R_{i_k}, τ_k = ±1 its sign (τ_k = −1 when M(y) of P'_k is odd) and u_k the sign of the walk point W_k, so W_k = u_k·canon(C_k). The step gives f(W_k) = u_k·P'_k and canon(C_{k+1}) = τ_k·P'_k, so u_{k+1} = u_k·τ_k; this law is derived from the definitions of canon and sigma, and the replica check found no failure of it.

*Equal-index forcing (derived).* If i_{k+1} = i_k and τ_k = −1 then P'_{k+1} = −canon(C_k), so C_{k+2} = C_k (a 2-cycle) and τ_{k+1} = −1: the pair (τ_k, τ_{k+1}) = (−1, +1) does not occur. If i_{k+1} = i_k and τ_k = +1 then P'_{k+1} = canon(C_k) + 2R_{i_k} and τ_{k+1} is not forced. The equal-index pairs with τ_k = −1 are the 2-cycle entries of Lemma E1(b), in the entry accounting, not the fresh-window bounds.

*abab closure law (derived).* For a pattern (i, j, i, j), i ≠ j, the four steps close the cycle, C_{k+4} = C_k, iff τ_k ≠ τ_{k+1} and τ_{k+1} ≠ τ_{k+2}. This is formal structure, not an excluded pattern: under the measured fresh-arrival law its conditional probability given the pattern is 1/4. On the replica it predicted every closure: 679 of 679 at A = 4 (2,807 windows) and 248 of 248 at A = 8 (986 windows), none outside the rule.

*Measured premise L-ori (fresh arrivals).* At a fresh class τ is uniform given the history: a measured premise, not a verified law. The replica's arrival sign has share 0.5003 over 352,580 arrivals (z ≈ 0.4; seven runs of 60,000 steps at A = 4 and A = 8, `evidence/sign-law.json`, key `fresh_arrival_runs`), and the non-forced share 0.5008 at A = 64 (n = 387,972, z ≈ 1) and 0.4987 at A = 8 (n = 77,464) are consistent with 1/2 (keys `tau_transition_runs`, `summary`). This is used for the 2-cycle entry rate 1/(2A), for Lemma E1 and for the window bounds of 2.3. These are small-A replica measurements, not production data; `evidence/sign_law.py` regenerates them from `measurement/walk_ref.py`.

*Observation in the stored production counts (not evidence for L-ori).* The seed-1 counts (the predecessor's `evidence/measurement/seed1.json`, omitted for the size limit) rise by 1 and 3 at step 64 for the two input sets against 0.5 expected each (pooled 4 against 1.0, one-sided Poisson p ≈ 0.02), while the cumulative count at step 64 (65 against 64.0) is consistent with the law, so the excess is post hoc (`evidence/sign-law.json`, key `production_steps`). The stored records carry no signs or adder indices, so this tests neither L-ori nor W.

**The good event G.** Let G be the event that no trail has (i) an accidental equality of two formally distinct
points, or (ii) a class repeat other than a formal repeat of length at most 64. A formal repeat of length at most 64
is a closure in the sense of 2.3: the walk returns to an earlier class with coefficients that agree for every table,
and it is detected and escaped (Section 1.3 step 5; Lemma R). Summing the bounds above over all trails started in a
run that does not fail on its counts:

  Pr[not G] ≤ DELTA_LOOP := K_TRAILS · (TRAIL_MAX²/2) · |Γ| · (2/r + 1/(r − |Γ| TRAIL_MAX²/2)) + A²|Γ|/r,

where K_TRAILS = DP_MAX + 2K counts the trails started: the K initial trails, at most DP_MAX trails started after
records, and at most K trails started after abandonments. By D2, every re-entry into an escaped cycle lies outside
G, so DELTA_LOOP bounds the re-entry probability. Lemma R and Lemma D hold for every trail, on G or not.

**Cut-off and false filter matches (DELTA_CUT).** A replay is cut off after 65n steps (section 1.3, step 7). A
divergence between a live trail and its replay needs a cycle that one anchor schedule detects and the other misses; a
false 64-bit filter match alone cannot cause this, since an escape happens only after an exact key check, so a false
match changes when the check runs and not which escape happens (Lemma R). A genuine cycle of length at most 64 is
detected under the schedule section 7 states, which compares every offset 33..64 after each anchored round (a replay keeps its own every-step schedule, a superset), so a divergence needs a cycle of length 65 or more
that one schedule detects and the other misses: a structural repeat of length at least 65, whose union over all fresh
starts is bounded by DELTA_LONG (section 2.3). A within-trail coincidence is bounded by DELTA_LOOP. Hence a cut-off
replay, which fails the run, has probability at most

  DELTA_CUT := DELTA_LOOP + DELTA_LONG ≈ 6.1·10^{−17}.

No per-trail term is charged for false filter matches, because a false match alone cannot cause divergence.

**Lemma L (re-entry).** By D2 every re-entry into an escaped cycle is a within-trail coincidence, so its probability
over the run is at most DELTA_LOOP and the event is contained in the failure term DELTA_CUT.

### 2.3 Entries (Lemmas E1 and E2)

**Definitions.** A step t is an *entry* if its source point lies on a fresh class (a class no walk has visited before, section 2.0) and its successor lies on a class already visited in the trail. A step whose source is on a
visited class is *persistent*. Only entries can open a closed cycle. A *window of length ℓ* ending at a source is the
ℓ steps before it; its *signed count vector* v_ℓ ∈ Z^A has entries Σ σ over the window's steps with index i.

**The formal structure of an entry.** Steps do not change the S-coefficient, so a formal repeat of a visited point has
the same S-coefficient and the sign +; it holds for every table iff the vector v_ℓ + σ_S e_{i_S} vanishes, where
i_S and σ_S are the index and orientation of the source's own step. Hence, for each ℓ ∈ {1, ..., 64}, a window ℓ gives an
entry target (i_S, σ_S) = t_ℓ exactly when −σ e_i = v_ℓ is a unit vector, and no target otherwise. On G every revisit is a
formal repeat of length at most 64, so an entry is one of these windows.

**Finite-key qualification for E1/F.** Their exact uniform-index laws hold in probability.md's coupled
ideal experiment, with total finite-key coupling failure charged there at 2^-80. They are not exact
residue laws for a raw key uniform on [0,p). The ideal index remains uniform after DP=false.

**Lemma E1 (entry probability).** Let F_{t−1} be the history before step t and the source of step t be on a fresh
class. Then, on G:
(a) Pr[entry at t | F_{t−1}] = w_t/(|Γ| A) ≤ P' := 64/(|Γ| A), w_t the number of distinct targets t_ℓ. The index i_S is
    uniform over A given the history (K_C is uniform and independent of F_{t−1}) and σ_S is uniform by L-ori (a
    measured premise, §2.2), so each target has probability 1/(|Γ| A).
(b) The 2-cycle entry (ℓ = 1) has conditional probability 1/(|Γ| A) = 1/(2A) given L-ori (§2.2): the predecessor P of
    the source is S − σ' R_{i'}, so the successor equals P exactly when i_S = i' and σ_S = −σ'. The "−P" case is
    impossible, since the successor's S-coefficient equals S's and would have to equal −(that of P), not zero mod r.
(c) The longer entries (ℓ = 3..64) have, in expectation over the history, rate at most c/A² with c = 3. Two facts about
    a window v_ℓ: (i) its coordinate sum is the sum of its ℓ units, so v_ℓ is a unit vector only for odd ℓ, for every
    window and with no premise; (ii) domination, since a fresh-source window's indices and units are dominated by i.i.d.
    uniform data (section 2.0), gives P_3 = 3/(2A) − 3/(4A²), P_5 ≤ 15/(4A²) and P_ℓ = O(A^{−(ℓ−1)/2}) for ℓ ≥ 7.
    Hence E[w'_t] ≤ 3/(2A) + 4/A² on these windows and the longer-entry rate is
    E[w'_t]/(2A) ≤ 3/(4A²) + 2/A³ ≤ 3/A².

Hence the entry probability is bounded by p_c := 1/(2A) + 3/A² + 2^{−63} in expectation over the history, where the 2^{−63}
term is the false 63-bit filter match, which the exact check rejects. The conditional bound of (a), P' = 64/(2A), holds
in the stopped discovery filtration; probability.md uses its stopped form. Only fresh-source windows can produce entries (derived). On G a formal repeat cannot span an escape, because each escape point has a new depth (D1); and within a segment the first revisit is itself an entry, after which the walk stays on visited classes until its escape. So a window that contains a persistent step is not the window of any entry, and the windows counted in (c) and in the long-closure bound are fresh-source windows only.

**Coupling.** The conditional bound of (a) uses only the source's own key and orientation, which are uniform given the history (M1, L-ori). The expectations in (c) and in the long-closure bound use the domination of section 2.0, which holds conditionally on every history of unvisited classes. No window premise is used.

**Lemma E2 (persistent steps).** On G, each entry is followed by at most 128 persistent steps before its escape,
deterministically. A 2-cycle entry is detected at the entry step, so it has no persistent step. For a longer entry, the
walk enters the cycle at the entry step, the anchor is set at a cycle member within 64 steps (the next anchor round), and
the walk returns to it at the first offset in 33..64 after that round, a further 32 to 64 steps. So the count is at most
64 + 64 = 128 (it was 64 + ℓ ≤ 128 on the previous every-step schedule; the schedule now buys the same constant). Consequently

  processed steps ≤ (fresh-source steps) + 128 × (entries).

Persistent steps are ordinary processed steps, charged the common cost C within S_max (section 3).

**Long closures (length 65 or more; Lemma F).** On G there is no formal repeat of length 65 or more (G(ii)). A closure of length d+1 ≥ 65 closes a fresh source C to a point at distance d ≥ 64. Escape counts agree on every trail (no escape count can exceed TRAIL_MAX, and ord_r(2) > 2·TRAIL_MAX, section 2.1), so the net adder vector of the d steps is a single adder, and every index but one occurs at least twice with zero unit sum. Let R_s = max(63, 2s − 2). Given any prefix, the d−1 fresh steps are dominated by i.i.d. data (section 2.0), so for every start the conditional probability of such a closure is at most

  Q_FAR := Σ_{s<A} A·C(A,s)·|Γ|^{1−s}·(s/A)^{R_s}·A/(A−s) + TRAIL_MAX·A·|Γ|^{1−A}.

Numerically log2 Q_FAR = −311.5 (the largest term is s = 33). The union bound over the at most S_max fresh starts gives

  DELTA_LONG := S_max·Q_FAR ≈ 2^{−184} ≈ 4·10^{−56}.

This is the long-closure term of DELTA_CUT (section 2.2).

### 2.4 Abandonment (superseded)

The original selected/fresh-tail abandonment argument is replaced by probability.md's common-phase,
all-abandonment envelope. No old abandonment allowance is used in the current claim.

### 2.8 Merges: Lemmas P and M

A *merge* is a visit by a walk to a class visited earlier by a walk of another trail (start and escape points included); a step into a class of the walk's own trail is an entry (section 2.3). Let τ be the first batch step with a merge. A *root segment* is a trail's path up to its first merge; by section 2.0 the fresh-source bounds of sections 2.3 to 2.5 are exact on root segments.

**Lemma P (continuation).** On the good event G, the classes a walk visits from a class C up to its first distinguished class depend only on C. *Proof (outline).* The class, the key (hence the adder index and the distinguished test) and the escape rule depend only on the class (section 1.2), and the sign of the walk's point does not enter the class sequence. Every history reaches the same cycle set (Lemma R) and escapes to the same 2·canon(M). ∎

**Resolution.** On G, without abandonment or another charged abort, the common continuation of
Lemma P reaches a common first DP and resolves. probability.md bounds any abandonment, including
after a merge, without conditioning a known selected prefix on an independent geometric tail.
The original selected-tail Lemma M is deleted and supplies no premise for the current claim.

## 3. Computation

**No reference operation.** Every instruction of a walk step is charged at 1/C with C = 1000, as COST_MODEL.md's
H + R/C with H = 0.

[Omitted: the superseded common-path charge (C_COMMON = 847), coin-uniformity remark and field-arithmetic description; see proof.md. `measurement/p256_field.h` is unchanged and used by the C paths.]

**3.1 Worst cases of the rare-path operations.** Each rare event is charged a forced-path maximum, not a sampled
average. The primitives are f_single (adder selection and the unbatched addition), add_single, escape (its least-key
search, canonical point and doubling), same_class and point_is_dp. For each the data-dependent branches were
enumerated from the source and the disassembly, and every conditional correction in a non-branch-free field routine
was replaced, for the measurement only, by a branch-free routine costing at least as much on every path (a
bounding build, checked equal to the original on random and edge inputs). The remaining control variants (the
equal-x and equal-y tests, the parity of the least member, the cycle length, the memcmp outcomes of the library
comparison) were forced by chosen inputs and the maximum over them taken; P-256 has no rotation search and no
exact-k fallback (measurement/step.c), so there is no such variant. Each count is linear in the number of calls and
identical for a second seed of inputs. Measured maxima (Callgrind Ir per call, organizer image and flags,
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

A 2-cycle check costs at most f_single + same_class + escape(2) + point_is_dp plus at most
1,024 instructions, and an anchor check at most f_single + escape(65) plus, for each of its at most 64 look-ahead
steps, f_single + point_is_dp + same_class with the same overhead per step. The overhead bound is far above the
measured control code of these functions, and the sums above the longest path through the code. Scalar
multiplications and replays use f_single's worst case, an upper bound for an addition or a doubling.

**Measurement (H3, a measurement assumption).** [The predecessor's common-path measurement is omitted; see proof.md. Rare paths are out of line, so each path's cost is separable and is priced by the forced maxima above.]

**Rare-path checks.** Distinguished points and false 64-bit filter matches do not occur in a production-size sample,
so small builds of the same source and replica exercise them (`evidence/measurement/trail-paths.json`, omitted for the size limit and byte-identical in 5a258d4a; 2^6 adders (2^4
in the last build), 2^12 walks): DP_BITS = 4 gives restarts with 2-cycles; DP_BITS = 12 with a one-period cap gives
abandonments; 6-bit filters give 2-cycles, longer cycles and look-ahead distinguished points; 2^4 adders at DP_BITS = 8
over 800 steps gives 94,624 2-cycles and 2,959 longer cycles escaped, no escape point repeated. In every build the
replica matched every final point, every event count (the repeat counter included) and a digest of every distinguished
record, look-ahead records included.

**8.2 Sampling caps and the distinguished-point table.** Every scalar is drawn by rejection sampling from [0, 2^256) with at most D = 4 draws, and every start and adder point attempted at most D = 4 times (retried only when it is O); the next draw or attempt aborts the run. Pr[abort] ≤ N_SCALARS·p^D + (A + K_TRAILS)·r^{−D} =: DELTA_CAP with p = 1 − r/2^256 ≈ 2^{−32} and N_SCALARS = 2(A + K_TRAILS), numerically about 2.5·10^{−18}.

The table has C = 2·DP_MAX slots of 128 bytes, open addressing with linear probing, index (key(W) >> 76) mod C (uniform up to a bias of about C/2^180), and load at most 1/2. A lookup or insert that examines L = 2^10 occupied slots aborts the run. *Probe tail:* if an operation examines L occupied slots from h, the maximal occupied run [s, s+ℓ) containing h has ℓ ≥ L and holds ℓ records hashed into the window [s, s+ℓ), since slot s−1 has always been empty. Under M1 the window count is binomial with mean at most ℓ/2, so Chernoff gives Pr ≤ e^{−cℓ}, c = ln 2 − 1/2. A union bound over the DP_MAX operations, the C window starts and ℓ ≥ L gives DELTA_PROBE ≈ 2.7·10^{−44}.

*Charge.* A table operation is charged at L probes times PROBE_COST = 2^11 instruction-equivalents, so DP_COST = 2^21 (raised from 2^20); the per-probe allowance covers the slot fetch (memory access), the comparison and the control, an allowance, not a measurement. A draw is charged DRAW_COST = 2^16. A start point is charged at its worst case, START = D(2·SCALAR_MUL + 2D·2^16) + 2^14, with SCALAR_MUL = 768·41,190; the adder table is charged ADDERS·START and the initial starts K·START. DELTA_CAP and DELTA_PROBE are added to delta.

**Total.** S_max = K·(L + TRAIL_MAX) ≤ 2^{127.542440} walk steps. The predecessor's own component table and rounded
total are superseded by proof.md's ledger, which re-prices every term at this revision's S. All processors' work is
summed, and elapsed time plays no part.

## 4. Memory and advice

Peak memory is at most the sum of 2·DP_MAX·128 bytes for the distinguished-point table (key, the trail's coins,
hashing overhead), K·512 bytes of walk state (points, keys, coins, trail counter, batch buffers), 2^16·2·64 bytes of
adder images and 16 MiB. That is 2^76.542441 bytes at the predecessor's A; this revision's formula and value, which
the ledger charges, are in proof.md. Nonuniform advice is 0 bytes.

## Small-key adders and reversed differences (5a258d4a)

In this section "previous kernel" is 70c0711e's.

**Algorithm change.** Step 2 of the procedure (predecessor Section 1.3) became: for each adder i < A, repeat at most
ADDER_TRIALS = 2^40 times: draw fresh c, d under the caps of Section 8.2, set R = [c]P+[d]Q, accept when R != O and
key(R) < 2^224 (top limb of the Montgomery x below 2^32), then tabulate R_i = R and -R_i with coins (c_i, d_i); an
adder that exhausts its trials fails the run. The current test is in "The y3 fold". Nothing else in the walk, its
tests, its distinguished points or its resolution changes.

**Kernel change.** The batch forms D = x(W)-x(R) (forward block FLOOP and the C code for walks B-1 and 0) and
t = x3-x(W) (backward block BLOOP_FIELD_S and the C step `k_lxy` of walk 0): 5a258d4a's `measurement/gen11.py` swaps the two
operands of exactly those five four-limb subtraction chains (each edit asserted to match once; same instruction count
and registers). The inverse chain uses D' = -D, so iw' = -1/D, lambda' = dy*iw' = -lambda, lambda'^2 = lambda^2 and x3
is unchanged; t' = -t gives y3' = t'*lambda' = t*lambda. PRE[w] becomes (-1)^(w-1) PRE[w], a uniform unit given
the walk under the per-round blinding, so the product-carry sites keep their law.

**Why two corrections disappear.** For a canonical u in [0, p), u - x(R) borrows iff u < x(R): D' = key(W) - x(R)
iff key(W) < x(R), the second x3 subtraction (lambda^2 - x(W)) - x(R) iff its canonical first difference is below
x(R). With every x(R) < 2^224 each has probability below 2^224/p < 2^-31.99999 for a uniform operand and *every*
adder, so neither site needs a table-mean pad or a Hoeffding term over the table; before, D = x(R)-x(W) borrowed with
probability 1 - x(R)/p and the second x3 subtraction with x(R)/p, 5 instructions per step on average whatever the
table. CORR changed by -2*5*(1/2+2^-10) + 2*5*2^224/p = -5.009766 to 17.458661.

**Measured on that revision's kernel** (5a258d4a's `evidence/measurement/smallkey.json`, not shipped here, whose scope states the toolchain, GCC
14.2.0 Ubuntu cross build with QEMU 8.2.2, and inputs): walk dumps (8 rounds of 2^16 walks, 2^16 adders) byte-identical
to the previous kernel's on a uniform-adder and a small-key set, event counters equal; fixed paths equal for both
kernels (41,098,898 at B=2^16, plain round 2, 27 above the organizer image's 41,098,871, the toolchain offset in `ancestor_quotes`;
656,973,460 clean and 657,055,820 with one two-cycle escape at B=2^20, A=2^20, organizer image 656,973,433); a static
diff (B=2^20, A=2^24) of exactly 100 swapped operand pairs in 25 subtraction chains (8 forward copies, 16 backward
copies, `k_lxy`) and 12 swapped load base registers in the C subtraction helpers, every rare-path function identical,
so the forced-path maxima still apply; D and the second x3 subtraction fired 0 times in 524,272 and 524,280
opportunities (previous kernel: D every time), the other sites kept their frequencies (subtractions 0.498-0.500,
products 0.249-0.251, square 0.333), and the meter, 17.40 per step (17.41-17.42 at 2^20 walks) against 22.41 before,
is below the charged 17.458661.

**Table cost and failure.** Each trial costs at most START + 2^10 instructions (START as before: capped draws, two
scalar multiplications, retries on O, stores; 2^10 for the acceptance test and loop control), so the worst case
A*2^40*(START+2^10) = 4,707,349,357,454,119,741,916,446,720 instructions = 2^81.961189 units is charged as
preprocessing. A trial point is uniform on the group (c is uniform and independent of d and x). For the key rule
alone, the Weil bound on the quadratic-character sum over x < X gave an error below 3 sqrt(p) (1+ln p) < 2^137.07
against the ledger's 2^160; the current 2-D count is in "The y3 fold". So at least 2^224 - 2^160 - 3 group elements are accepted, q_acc >= 2^-32.0000000004, and all A adders finish
except with probability A (1-q_acc)^(2^40) < 2^-345.3, added to DELTA_CAP (2.5477e-18 with the 2*A*2^40 extra scalar
draws); expected trials per adder are about 2^32.

**Probability space and accidental relations.** Every walk point is
[l]P, l a linear form with integer coefficients fixed by the transcript in
the trail's start logarithm t = a + b x and the adder logarithms
s_i = c_i + d_i x; the walk, its keys and its tests are functions of these
logarithms only, and x enters only at resolution. For fixed x the map
(c, d) -> (c + d x, d) is a bijection of Z_r^2, so each start logarithm is
uniform on Z_r and each accepted s_i is uniform on the accepted set S of
logarithms whose point passes the test (key([s]P) < 2^224 there; the current test in "The y3 fold"), of size
m >= M_ACC (2^224 - 2^160 - 3 for that test; the current value in proof.md) and closed under negation; all are independent of each other
and of x, and given them every d_i and every b is uniform (so the
degenerate-match bound, which uses b and d, holds as stated). **Acceptance
does constrain the s_i,** so an accidental relation among adder logarithms
is likelier than over Z_r: a nonzero linear form in the s_i (and starts)
vanishes with probability at most 1/m, since fixing the form and an index
j with nonzero coefficient pins s_j to one value, which a uniform element
of S takes with probability at most 1/m (1/r on a start coefficient). Two
adders satisfy R_i = +-R_j with probability 2/m, not 2/r.

Every accidental-coincidence bound of predecessor-proof.md Section 2.2 (the in-trail pair bound 2|Γ|/r, its
conditioning on earlier non-coincidence, Lemma G's |Γ|t/(r - |Γ|t(t-1)/2) and the table term A^2|Γ|/r) is a union of
such events, so each holds with r replaced by M_ACC; this is applied to all of them, including terms that also
contain a fresh start logarithm, which only over-counts. Their side condition becomes |Γ|TRAIL_MAX^2 <= M_ACC/2
(asserted). In the ledger DELTA_LOOP, hence DELTA_CUT, rises by the factor r/M_ACC < 2^32 from 5.8542e-17 to
2.5144e-7; it is subtracted, and at that revision's tau the success bound was 0.39000078761, a margin of 7.876e-7
above 0.39 (asserted floor 5e-7; tau, S and every work term unchanged there; the checkpoint-only revision re-took tau, section
"Complete cost and success ledger"). Formal closures (D1, D2, Lemma F) are combinatorial facts about transcripts with
i.i.d. class data and do not involve the adder logarithms' distribution, and merges between different trails are the
intended collisions, governed by H1.

M1 and L-ori concern walk classes. The selection conditions only the 2^24 adder classes; the walk ever processes a
point in one of them with probability at most S*A/N < 6.94e-32 (`DELTA_ADDER_CLASS`, subtracted). H1 (random mapping
for the class walk with the A-adder penalty epsilon) depends on the table only through its A elements, whose
logarithms are independent and uniform on S, defined by coordinate bounds with no known relation to discrete
logarithms; extending H1 to this table is a disclosed premise at the level of H1 itself, not a theorem. Toy support
(5a258d4a's `smallkey.json`, `toy`; key rule only): the complete procedure on 28-bit a=-3 curves at
tau=0.4943, adders redrawn until key < p>>4, gives 513/1200 successes (three seeds) against 436/1000 with uniform
adders, two-cycle rates 2.37e-4 to 2.46e-4 against 2.41e-4 (model 1/(2*2^11) = 2.44e-4), no replay mismatch. The
declared experiment uses `adder_key_shift` 2 (keys below p>>2) to fit the runner's time limit (shift 4 did not).

## Checkpoint-only cycle detection (ae463634)

The lane record is {sel, anchor, tper}, 32 bytes (aligned to 32): the two-steps-back filter (the previous key
and its compare) and `rare_2cycle` are gone, so each step compares only its new key's low limb with the anchor filter,
and `main` and `new_trail` no longer write a ~0 to word 0 of Wn[w], which nothing reads before the step's own store.

*Detection.* Every walk's new point is set as the anchor in rounds t = 0 mod 64, and the new key is compared with
it in that round's anchor tail and in every plain tail of rounds t with 33 <= t mod 64 <= 63; rounds with
1 <= t mod 64 <= 32 run a tail with no anchor compare. `rare_anchor` then iterates f from the new point, testing each
look-ahead point for the distinguished property first, until the new point's class returns (exact check) or 64
look-ahead steps pass. On return it holds all l members and escapes to 2*canon(M), M the least-key member of the whole
cycle, so the escape depends only on C's set of classes, not on the anchor's phase.

*Lemma R and D1.* An escape follows an exact check over every member of a closed cycle of length <= 65, all
tested and not distinguished, so Lemma R holds as written and its 2-cycle clauses are vacuous. Every cycle of length
l in 2..64 is caught by the schedule above, since the offsets 33..64 it compares after the anchored round contain a
multiple of l in every case (for l <= 32 the 32 consecutive offsets hold one; for 33 <= l <= 64 the first return is
the offset), so a replayed trail cannot diverge on a cycle of length <= 64 and DELTA_CUT is unchanged. The replay
keeps its own every-step schedule (section 1.3), which compares a superset. A step makes at most one exact check,
hence at most one escape: K_ESC = TRAIL_MAX and D1 stand.

*Persistent steps (E2).* A cycle of length l <= 64 therefore costs at most 64 + 64 = 128 persistent steps, which is
the constant E2, epsilon and the abandonment envelope already carry. So processed <= fresh + 128*entries and
the historical abandonment parameters held unchanged there; the current probability.md supersedes that tail, and P' = 64/(2A) counts the l = 1 window.

*Replay.* A replay meets the same cycles (Lemma R). On a cycle of length l it takes at most 64 + l walk steps
(the 64 to the next anchor round, then the return l steps later, which the every-step replay schedule compares), while
the live trail took at least l on it (the steps producing l-1 members and the return whose point the escape
replaced); on G distinct cycles of a trail share no class (D2). So a replay of n live steps takes at most 33n walk
steps, inside the 65n cut-off, each with at most 68 operations (one evaluation, 66 for an anchor check, one escape):
33*68 < 2^13 = REPLAY_OPS per live step (asserted).

*Collision exponent and cost.* A short-cycle entry wastes at most 128 persistent steps (E2) and 64 look-ahead
points, within the 192 states per entry that epsilon's closure margin charges with the collision-exponent factor 2;
the two-cycle term moves from 2*(3/A), when a two-cycle escaped at entry, to 2*192/(2A): epsilon = 196/A + 288/A^2 +
768/A^3. A two-cycle entry is charged wL, the anchor check with 64 look-ahead steps and a 65-cycle escape, in the
meter row and in B2 = n2*wL; the false-filter row keeps the former two-filter debit 2^-63*256*(w2+wL); persistent
steps are processed steps inside S.

*Small-key table.* The entry rates p2 = 1/(2A) and 3/(4A^2)+2/A^3 are conditions on the keys and signs of successive
walk points (an index repeat with a sign flip, and the longer-window patterns of Lemma E1), not relations among adder
logarithms, so the small-key table leaves them as derived (toy rates: section "Small-key adders"). A persistent step's corrections are
metered and capped like every processed step, and H-trail, an average over whole trails, covers them.

*Kernel.* Each of the eight plain and eight anchor block copies drops the previous key's `ldr x5,[x25,#16]` and its
`ccmp`, and both rare exits drop the `str` of that key (34 instructions in all; static counts and checks in
ae463634's `branchy.json`, key `checkpoint_only`, quoted in `ancestor_quotes`).

In that checkpoint-only revision, walks differ by design: two-cycles escape at the next anchor, not at entry. In n=8 dumps of 2^16 walks exactly the
walks the small-key kernel escaped from a two-cycle differ (2 and 3 on the two input sets, out of 2 and 3 escapes);
over 130 rounds 62 and 65 differ, where that kernel escaped 63 and 67 two-cycles and the checkpoint-only one 63 and 67
at an anchor. In six builds (2^16 adders: 8 rounds, a 64-round stress with 6-bit filters, DP 2^-7 and a
3-period cap, and 130 rounds; 2^4, 2^6 and 2^10 adders) on both input sets a replica under the same rule matches every
final point, event count and record digest (same key); on the small-key kernel's
130-round dump it finds 64 mismatches. At 2^4 and 2^6 adders, where a two-cycle starts every 32 or 128 steps, more
trails reach the cap; at A=2^24 their cost is the epsilon and 128-step terms above.

## Deferred canonical test (e2973980)

In this section "here" is the deferred-inv kernel (e2973980); the later revisions keep its test
(proof.md items 10 and 11).

In the checkpoint-only kernel inv had the
canonical finish above (3 on the common path) and each tail ended `tst; ccmp F,x6,#4,ne; b.eq 2f` (Z = 1 on a
rare event). inv is next read by the
next step's INV' product (below p), the rare exit or the C step; nothing between reads or writes x19..x22 here
(all 16 copies of each build). The finish is
carry-only (`b.cc` over [subtract p; meter += 5]), leaving r = T - c*p (r = T < 2^256 if c=0, r < p if c=1), and the
tail's `b.eq 2f` becomes

    ccmn x22,x28,#10,ne     // no event: C = (r3 >= p3); event: NZCV=1010
    b.hs Lz

directly after the `ccmp` in all 16 tails of every build. Without an event `ccmn` adds 2^32-1 to r3, so C = (r3 >=
p3), and N = 0 whenever C = 1 (the sum less 2^64 is below 2^32); with one, N = C = 1. Lz is `b.mi 2f` (event: to the rare exit), the former slow block's conditional
subtraction (four limbs, four `csel`), meter += 11 and a branch back. r3 < p3 implies r < p, and r < 2^256 < 2p, so
inv is canonical when the step ends on every path: the common path, Lz, and the rare exit `2:`, which now first runs
the same eight-instruction subtraction (a period-cap exit jumps there before the `ccmn`). `measurement/gen12.py` (in e2973980, omitted here)
derived that revision's `bloop_branchy.h` and `step.c` from the checkpoint-only files.

*Cost.* With c=0 and r3 < p3, `b.cc`, `ccmn`, `b.hs` replace `b.cc`, `cmn`, `b.hs`, `b.eq`: one fewer instruction per block
step, plain and anchor. With c=1 the finish runs 5 over the common path (weight 4 -> 5) and the slow test 11. Lz runs
iff r3 >= p3 without an event, below 2^-31, charged 2^-24. The inv carry is the same event as before (T >= 2^256
for the same T); c(pi, W) and H-trail hold, only the
weight growing (by (1/4+2^-10) + 2^-24 + 6*2^-31 = 0.250977). c(pi, W) <= 86, all weights less 4: y3''s c=1 fix-up
and its slow block exclude each other, while inv's c=1 and Lz can both occur.

*Exits.* On an event the step runs `b.cc`, `ccmn`, the taken `b.hs` and `b.mi` and the subtraction (12) where the
checkpoint-only kernel ran `b.cc`, `cmn`, `b.hs` and the taken `b.eq` (4): 8 more, 9 above the clean step PLAIN_FIXED
charges (7 for a cap exit). Measured (2^16 walks, rounds 2..128 against round 1): FIXED falls by 127*65,535 - 9*34
on input set 1 and 127*65,535 - 9*35 on set 2, with 34 and 35 cycle escapes in both kernels. The ledger charges 10 per batch exit: inside wL
for an anchor check (wL 2,789,687 -> 2,789,697) and K_TRAILS*10 for distinguished-point and cap exits, at most one
per trail.

*Evidence* (e2973980's `evidence/measurement/branchy.json`, `deferred_inv`, omitted here; its a6 long-cycle counts
are in `ancestor_quotes`). Ten dump runs of five configurations count 51,491,065 inv carries in
206,045,184 walk steps (0.249902, 35 standard errors below the charged 0.250977; per run, rounded outward, 0.248400 to
0.251889, the extremes with 2^6 adders, where over 233,000 cycle escapes per run repeat the same states). With
small-key adders at 2^20 walks the meter exceeds the checkpoint-only kernel's by 0.249831 and 0.249972 per step
(17.662485 to 17.667348, below CORR); Lz occurs in none of these runs. They and three more configurations match the
checkpoint-only build's output apart from the meter, which rises by the inv-carry count exactly; the block with the
tail test, Lz and the exit equals the exact model and the checkpoint-only block on 16,320,000 random and crafted
inputs, 0 mismatches (two single-instruction mutants: 12,581 to 47,410 in 200,000).

## Eight copies per loop (961165d1)

Copies equaled 1f942c41's instruction for instruction (forward 119, plain backward 570, anchor backward 575, slow
sets 30); `batch_step` differed only in the entry tree, the counter immediate and nops. Loop control per round
fell by exactly 3B/2-5 (1,572,859 at B=2^20); a re-entry after a rare exit costs at most 3, inside OVERHEAD.

## Paired reduction words (21edf4bd)

A product computes T = (a*b + M*p)/2^256 in four rows: row j adds a*b_j, takes m_j, the accumulator's low word
(-1/p = 1 mod 2^64), and adds m_j*(p+1)/2^64 = m_j*2^32 + m_j*p3*2^128 from the next word up. The paired-reduction
kernel forms m_j*p3 by `mul`, `umulh` of m_j and x29 = p3 and joins m_j>>32 and m_(j+1)<<32 (disjoint bits of one
word) by one `extr`: a product fell from 94 to 90 instructions and the square from 75 to 70 (21edf4bd's
`measurement/gen48.py`, omitted here). Each row adds the same integers in another grouping, so T, the carry and every
finish and meter outcome were unchanged; the five accumulator words stay below p*2^64 + p + 2^128 < 2^320 (asserted
in `accounting.py`; evidence in 21edf4bd's and 205a498f's `branchy.json`, `paired_rows`). Here only the square and y3' keep this form, their product text unchanged (from `mul x1,x0,x2`
on, the backward block is the y3-fold revision's with proof.md item 13's edits).

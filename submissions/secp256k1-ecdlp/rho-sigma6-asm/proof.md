# secp256k1-ecdlp: parallel Pollard rho on order-6 classes with a sigma walk, all arithmetic in instructions

Claim: for uniform x in [0,r) and Q=[x]P on secp256k1, the algorithm below returns x with probability at least 0.39. It uses at most **2^126.479214** target-group-operation equivalents in total (score 63239607/500000), including all preprocessing, failures, lookup and verification. Peak memory is at most 2^75.749960 bytes. Preprocessing is 2^49.151678 units, and no nonuniform advice is used. This is a classical cold attack on one target.

The success bound rests on one heuristic model, R (section 2.2): the random-function idealization of the adder walk on classes, in which each class's key, orientation and canonical data are i.i.d. uniform. Every conditional statement below is a theorem of R by deferred decisions. The escape-depth certificate of section 2.4 is proved. The cost rests on H3 (measured instruction counts) and on a work meter with a hard cap, whose expectation is taken in R. A model-free rare-work bound would cost about 10.9 bits and is not taken (section 10). This is a carefully accounted generic attack, not a practical break.

**Cost route.** Every instruction is charged at 1/C with C = 1000 (H = 0 in COST_MODEL's H + R/C), as inline AArch64 assembly. A walk step is charged 828.89 instructions. The walk stores the actual point W and adds sigma(R_i), so the class sequence is that of the standard walk. The budget is tau = 0.5303.

**Distinguished points and resolution.** Walks run in trails. Each trail starts from fresh coins and ends at its first distinguished point, which it records once with its coins. A key match is resolved by replaying the two trails (section 1.3, Lemmas R and D), so a fruitless cycle cannot store the same point twice, and a degenerate match has probability exactly 1/r.

Files: `README.md` points here. Evidence paths appear where the text cites them.

## 1. Algorithm

### 1.1 Notation and parameters

p, r and P come from the protected profile. The curve is y^2 = x^3 + 7, the cofactor is 1, and p ≡ r ≡ 1 (mod 3).

The endomorphism uses beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee and lambda = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72, with [lambda](x,y) = (beta·x, y). The automorphism group G6 = {±1, ±lambda, ±lambda^2} has (r−1)/6 nonzero classes, and N = (r−1)/6 + 1. Field values are compared as reduced integers in [0,p). The units of G6 are the elements of Z_r^* of order dividing 6, and |Gamma| = 6.

Parameters:
- A = 2^16 adders.
- K = 2^32 walks: 2^16 processors × 2^16-walk batches.
- theta = 2^−60.
- L = 33310394024866356565296666183 ≈ 2^94.7500.
- CAP = 64/theta = 2^66. Trails have fewer than TRAIL_MAX = CAP + 64 = 2^66 + 64 steps.
- S_max = K·(L + TRAIL_MAX) ≈ 2^126.749959, the bound on walk steps in a run.
- DP_MAX = 2·theta·S_max + 2K ≈ 2^67.7500, the number of stored distinguished points.
- K_TRAILS = DP_MAX + 2K ≈ 2^67.7500: the trails started (K initial trails, DP_MAX restarts after records, and K restarts after abandonments).
- K_ESC = TRAIL_MAX: the escape-depth bound of section 2.3, certified by the order-of-2 check of section 2.4.

### 1.2 Classes, keys and the sigma walk

For W = (x, y) ≠ O:
- **Class key.** key(W) = y if y is even, else p − y. key is constant on the class, and it determines the class: y^2 = x^3 + 7 fixes x up to a factor beta^j.
- **Canonical point.** canon(W) = (beta^k·x, (−1)^s·y), where k minimises beta^k·x and s = y mod 2. As a scalar, canon(W) = eps_W W with eps_W = (−1)^s·lambda^k.
- **sigma.** sigma_W = [(−1)^s·lambda^(−k)], so that W = sigma_W(canon(W)).
- **Walk step.** i(W) = (key >> 1) mod 2^16, and f(W) = W + sigma_W(R_i) = sigma_W(canon(W) + R_i). The class sequence is exactly that of the standard class walk [W] → [canon(W) + R_h([W])]. The six images sigma(R_i) are tabulated per adder.
- **Distinguished point.** W is distinguished when bits 17..76 of key(W) are zero.

**Deciding k without beta·x (Lemma 1).** Write x = Σ_{i<4} a_i·2^(64i), and let C1 = ⌊2^384·beta/p⌋ = Σ_{j<6} C1_j·2^(64j). Define g_0 = a_3 and g_1 = ( Σ_{i+j=5} lo(a_i·C1_j) + Σ_{i+j=4, j≥1} hi(a_i·C1_j) ) mod 2^64, so g_1 takes 8 multiplications. Let f_j = 2^64·frac(x·beta^j/p).

**Lemma 1.** f_0 ∈ [g_0, g_0 + 2). If g_1 ≤ 2^64 − 10, then f_1 ∈ [g_1, g_1 + 9).

*Proof.* For j = 0: x/2^192 ≥ a_3, and since p < 2^256, f_0 − a_3 < 1 + 2^64·x·(2^256 − p)/(p·2^256) < 2. For j = 1: beta/p = (C1 + e)·2^−384 with e ∈ [0,1), so f_1 = (x·C1·2^−320 + x·e·2^−320) mod 2^64, and the second term is below 2^−64. Modulo 2^64, x·C1·2^−320 = Σ_{i,j} a_i·C1_j·2^(64(i+j−5)) keeps the i + j = 5 terms (lo) and the i + j = 4, j ≥ 1 terms (hi); the rest contribute fractions below 1 each. Dropping them gives an error in [0, 8 + 2^−50), and the hypothesis excludes wrap-around. ∎

x + beta·x + beta^2·x ≡ 0 mod p, so f_2 ≡ −f_0 − f_1 mod 2^64. With g_2 = (−g_0 − g_1 − 22) mod 2^64, the real f_2 lies in (g_2 + 11, g_2 + 22] whenever g_2 ≤ 2^64 − 23. So every f_j lies in [g_j, g_j + 22]. The kernel accepts the least g_j as k only if all three g_j ≤ 2^64 − 23 and the least is more than 22 below both others; then k is exact. Otherwise it computes beta·x and compares exactly (the exact-k fallback, rate at most 2^−55 per canon call, section 6). `measurement/frac_check.py` checks Lemma 1 and the decision rule on 200,000 random and edge-case inputs per seed; the largest floor gap seen for g_1 is 6. The batch path reaches this rule for every new point of every walk, so the fallback is on the common path (section 6).

### 1.3 Procedure

1. **Validate Q.** Check it is on the curve. If Q = O, output 0.
2. **Adders.** For each adder i < 2^16: draw c_i, d_i uniform in [0,r) from the operating system's randomness, fresh in every run and independent of x, under the sampling caps below; set R_i = [c_i]P + [d_i]Q, retrying while R_i = O within the cap; tabulate its six images.
3. **Trails.** Each of the K walk slots runs a sequence of trails.
   - A trail starts from fresh coins a, b uniform in [0,r), drawn from the operating system's randomness at every trail, including every restart, under the sampling caps, with W = [a]P + [b]Q (two scalar multiplications, charged in section 8). The slot keeps the coins until the trail ends. Coefficients are not tracked while walking.
   - A trail ends at its first distinguished point (step 6), or is abandoned when it crosses more than 2^60 anchor boundaries (step 5) without one, that is after between CAP − 64 and CAP = 64/theta steps, or when an addition gives O. Either way the slot starts a new trail from fresh coins. Every trail has fewer than TRAIL_MAX steps.
   - Each slot takes L steps (the main phase), then finishes its current trail and stops. The run fails if more than K trails are abandoned.
4. **Step.** W ← f(W), as one batched incomplete affine addition: Montgomery's simultaneous inversion over the 2^16 walks of a batch (3M), then lambda, lambda^2 and y3 (2M + 1S), with one Fermat inversion per batch. If the batch product is zero, the affected walk uses the unbatched complete formula. The canonical rotation of each new point is decided by the rule of section 1.2.
5. **Fruitless cycles.** Each walk keeps two filters, both the low 64 bits of a key: the key of the point before its current point W, and an anchor, its key at the last batch step divisible by 64. Each walk is anchored once every 64 of its own steps. The filters can match falsely; they only decide when an exact check runs.
   - *Order.* The distinguished-point test of step 6 comes first, on every point the walk produces: the new point, every look-ahead point and every escape point. A distinguished new point ends the trail, and no cycle check runs.
   - *2-cycles.* If the new key matches the first filter, the walk evaluates f(new) unbatched. If key(f(new)) = key(W) exactly, {W, new} is a 2-cycle of the class walk. Both points have been tested, and the walk replaces new by 2·canon(M), where M is the element with the smaller key.
   - *3..64-cycles.* Unless a 2-cycle escape just happened, if the new key matches the anchor filter, the walk iterates f from the new point, unbatched, up to 64 times, testing each look-ahead point first. A distinguished look-ahead point s ends the trail at s, which is the trail's record, since the walk would reach it next. If the new point's class returns, every point of that cycle has been tested, and the walk escapes to 2·canon(M), with M the cycle element of least key.
   - The anchor check also runs after a first-filter match that the exact check did not confirm. So a walk in a cycle of length ℓ ≤ 64 leaves it within 64 + ℓ steps, unless the cycle contains a distinguished point, which ends the trail first.
6. **Distinguished points.** W is distinguished when bits 17..76 of key(W) are zero. The test is applied to every point that a step or an escape produces. At a distinguished point the trail ends: the record (key(W), a, b), with the trail's coins, is inserted into the open-addressing table of section 8.2 (capacity 2·DP_MAX, probe cap L = 2^10), and the slot starts a new trail. A trail makes at most one record. The run fails if more than DP_MAX records arrive, and it aborts if any table operation examines L slots.
7. **Resolution.** The first record whose key is already in the table ends the run. No other key match is used.
   - Both trails are replayed, unbatched, from their recorded coins up to their first distinguished point, with the rules of steps 4–6 and anchors every 64 of the trail's own steps. The replay tracks (a, b) mod r: an image sigma(R_i) has coefficients multiplied by (−1)^s·lambda^(−k), canon multiplies them by (−1)^s·lambda^k, and an escape doubles the coefficients of canon(M). By Lemma R each replay reaches the recorded point.
   - The replays give the common canonical point as [a1]P + [b1]Q = [a2]P + [b2]Q. If b1 ≠ b2, output x = (a2 − a1)/(b1 − b2) mod r after checking [x]P = Q. If b1 = b2 the run fails. By Lemma D this has probability 1/r.
   - **Replay bound.** A replay of a trail with n live steps is cut off after 65n steps (REPLAY_FACTOR = 65). Each live step gives at most one walk step and one escape of at most 64 cycle steps, and a replay step uses at most 70 unbatched operations, so 65·70·n ≤ 2^13·n. A cut-off replay fails the run, counted in DELTA_CUT (section 3). The replay charge is always made (section 8).
8. **Work meter.** Every run is charged at most W_max = S_max·828.891183 + G instructions, with G the fixed global charges (initial starts, trail restarts, distinguished-point records, solve and verify, input validation, adder table, resolution replays). Every walk step is charged C_COMMON = 821 instructions, and each rare event its worst case (section 6). The meter aborts the run when the total metered rare work exceeds 64 times its expectation bound (Lemma M''), with probability at most 1/64, counted in delta. The run also fails if more than DP_MAX records arrive or more than K trails are abandoned.

**Sampling caps.** Every scalar is drawn by rejection sampling from [0, 2^256) with at most D = 4 draws, and every start point and adder point is attempted at most D = 4 times, a point being retried only when it is O. The next draw or attempt aborts the run. Pr[abort] ≤ N_SCALARS·p^D + (A + K_TRAILS)·r^−D =: DELTA_CAP, where p = 1 − r/2^256 is the rejection probability of one draw (about 2^−127.8 for secp256k1), N_SCALARS = 2(A + K_TRAILS) is the number of scalars, and DELTA_CAP < 2^−440. Each start and adder is charged at its D-attempt worst case (section 8.2).

**Escape.** A cycle escapes to 2·canon(M), with M the least-key member of the cycle. canon(M) lies in M's class, so the escape is one point doubling with no table, charged at its forced-path maximum (section 8.1). Before the first merge its class is unvisited, so its data are uniform by R.

**Merged trails** follow the same classes from their first shared class on (Lemma P, section 2.8).

**Lemma R (replays reach the recorded point).** The first distinguished point of a trail, and every escape before it, depend only on the trail's start point. They do not depend on when the filters fire or on the anchor schedule.

*Proof.* An escape happens only after an exact check shows that the point's class lies on a closed cycle C of length ℓ ≤ 65, and after every point of C has been tested and found not distinguished (the 2-cycle rule tests both points of {W, new}; the anchor rule tests each look-ahead point). Without the escape the walk would repeat C forever without a distinguished point, and 2·canon(M) depends only on C's set of classes. A look-ahead that meets a distinguished point ends the trail at the first distinguished point that f alone would reach next (a 2-cycle cannot lie strictly inside a longer cycle, so no escape could intervene). A false filter match either finds no closed cycle and changes nothing, or finds a genuine one and is handled as above. Every cycle of length ℓ ≤ 64 is caught: within 64 steps of entering it the walk is anchored at one of its points, and it returns to that class within ℓ more steps, where the anchor check runs. So under any schedule that anchors each walk at least once every 64 of its steps, the trail meets the same cycles, escapes to the same points and stops at the same first distinguished point. Every coefficient pair the replay carries represents its current point. ∎

**Lemma D (degenerate matches).** The first key match has b1 = b2 with probability exactly 1/r.

*Proof.* A trail records once, so the two records come from different trails t1 ≠ t2. Condition on all start points and adders. Given its start point W = [a + b·x]P, a trail's coin b is uniform on [0, r), because each b determines exactly one a. The coins are drawn from the operating system's randomness (section 2.2), so this holds exactly over the algorithm's coins, and coins of different trails are independent. Along a replay, b changes as b ← u·b + d, where u is a unit (order-6 classes, up to sign) or twice one, and u, d depend only on the path. So the canonical b at the distinguished point is U_i·b_ti + D_i with U_i invertible mod r, and P[b1 = b2] = 1/r. ∎

## 2. Premises, formal structure and checks

### 2.1 The collision count

The birthday count tau = T_steps²/(2N) and the adder penalty eps of section 7 are the random-mapping count of the class map under R (Teske 2001; Bernstein–Lange 2012; Bos–Kleinjung–Lenstra 2010). They are a consequence of R, not a separate premise.

### 2.2 The model R and the coins

**Model R (the one premise).** Each class C of the walk carries a key K_C, uniform on the key space (it fixes the adder index (K_C >> 1) mod A, the 64-bit filter bits and the distinguished-point test), canonical coordinate fractions, and an orientation g_C uniform on Γ. These data are independent across classes and independent of x and of the coins. The class reached from C is the algebraic image of C under the adder selected by its key. This is the random-function idealization of the adder walk on classes: the data of each class are i.i.d. uniform, and the map is algebraic on top of them. R is a single heuristic premise; no separate conditional law is assumed.

**Theorems of R (deferred decisions).** A class's data are drawn when it is first visited, so given the full history F of the run the data of every unvisited class are i.i.d. uniform and independent of F. Consequences used below: (T1) the index of a fresh class is uniform over the A adders and its orientation uniform over Γ, given F; (T2) the orientation law: the sign of an arriving fresh point is uniform, so the non-equal share of Lemma FR is exactly 1/6 (section 4); (T3) for a window of fresh classes the indices and units are i.i.d. uniform, which gives the coupling below; (T4) the merge count of section 2.8.

**Global filtration and coupling.** All walks of a run form one process, advanced one batch step at a time. A class is *visited* once a walk's point lies in it (start and escape points included), and its data are drawn once, at its first visit, so every trail reads the same class map. Given the history F, the next d−1 classes can be given data from i.i.d. draws ξ_1, ξ_2, … (deferred decisions), so for every event E of those data and every history F in which the next d−1 classes are unvisited,

  Pr[E and the next d−1 classes are unvisited | F] ≤ Pr_iid[E].

No data are resampled per trail, and nothing is assumed about visited classes.

**Randomness of the algorithm (not a heuristic).** The adder coefficients and the start coefficients are drawn uniformly from the operating system's randomness (`secrets.randbelow(r)` or equivalent), fresh in every run and independent of the target. Every trail, including every restart, draws fresh coins (a, b). Every probability in sections 3 to 7 is over x, the coins and the class data of R.

### 2.3 Formal structure, escape depth, D1 and D2

Within one trail with start S, each walk point has the formal form W = 2^m g S + sum_i n_i R_i, with g in Gamma, integers n_i, and m the number of escapes so far in the trail. Steps add sigma(R_i) (which is ±lambda^j R_i) and do not change the S-coefficient. An escape replaces M by 2 canon(M) = 2 g' M, which doubles the S-coefficient.

**Escape depth.** A step escapes at most once: a 2-cycle escape disables the anchor check, and an anchor escape replaces the point that the step produced. Every escape resolves an entry (section 4), so a trail of at most TRAIL_MAX steps has at most K_ESC = TRAIL_MAX escapes. The bound is deterministic: the certificate of section 2.4 excludes every k up to TRAIL_MAX, so no depth tail is needed.

**D1 (depth).** If E is the (m+1)-th escape point of a trail, its S-coefficient 2^(m+1) g' differs from the S-coefficient 2^j g'' of every earlier point of the trail (j ≤ m).
*Proof.* Equality would give 2^k in Gamma for k = m+1−j, with 1 ≤ k ≤ m+1 ≤ TRAIL_MAX. Section 2.4 excludes every such k. ∎

**D2 (re-entry is not structural).** A re-entry into an escaped cycle is a later point (depth ≥ m+1) in the classes of an earlier cycle (depth ≤ m). By D1 such a point is never formally equal to an earlier one, so every re-entry is a within-trail coincidence, which section 3 bounds.

### 2.4 Order of 2: no k up to TRAIL_MAX with 2^k in Gamma

**Claim.** For the secp256k1 order r and Gamma = {±1, ±λ, ±λ²} ⊂ Z_r^*, there is no k with 1 ≤ k ≤ TRAIL_MAX and 2^k ∈ Gamma (mod r).

*Proof.* The six distinct elements of Gamma satisfy x^6 = 1, so Gamma is the group of sixth roots of unity of Z_r^*, and |Gamma| = 6. Let C = 174723607534414371449, a 68-bit prime dividing r−1 in the published factorization r−1 = 2^6·3·149·631·107361793816595537·C·341948486974166000522343609283189 (checked in accounting.py). C passes Miller–Rabin with bases 2..41, which proves it prime since C < 3.317·10^24, and 2^((r−1)/C) ≢ 1 (mod r). By Fermat, 2^(r−1) = 1, so ord(2) | r−1. If C did not divide ord(2), then ord(2) would divide (r−1)/C and 2^((r−1)/C) = 1, which is false; so C | ord(2). If 2^k ∈ Gamma then 2^(6k) = 1, so ord(2) | 6k and C | 6k; as C is a prime other than 2 and 3, C | k. For 1 ≤ k ≤ TRAIL_MAX < C this is impossible. ∎

The bound TRAIL_MAX = 2^66 + 64 of section 2.3 is below C, with log2 C = 67.24. These checks, the factor identity, the primality test, 2^((r−1)/C) ≠ 1, and gcd(C, |Gamma|) = 1 with |Gamma| = 6 and C > TRAIL_MAX, are assertions in accounting.py.

### 2.5 Seeds and the restart pool

The measurement and the replica use fixed seeds (`measurement/walk_ref.py gen <path> <seed>`; `reproduce.sh` uses seed 1), and the declared run uses seed 2801 (section 7). The algorithm itself draws the table and the coins from the operating system (section 2.2). The measurement harness takes restart starts from a precomputed pool of SPARES = 1024 points, so restarts there reuse starts; the specified algorithm draws fresh coins for every restart and is charged its cost (section 1.3, step 3). The kernel was not changed.

### 2.6 Formal-structure test (evidence/formal-structure.json)

The test script `formal_test.py` (not included in this package) replays small builds (2^6 adders, 256 walks, 3000 batch steps, theta = 2^−6, seed 1) on real secp256k1 points and checks each walk point's formal vector and pair. Of 767,593 fresh-source steps, 1,992 are entries (1,980 two-cycle, 12 longer). The largest persistent run before an escape is 62 (bound 128, Lemma E2). The 1,980 formal 2-cycles have rate 2.621·10^−3 against 1/(6A) = 2.604·10^−3 (z = 0.29). All 1,992 escapes have 0 formal equalities with earlier points of their trail (D1); the largest escape depth is 3. The 12 longer entries are fewer than the 15.6 predicted, and no formal 4-closure occurred; this is not explained.

### 2.7 Degenerate pairs

Formally distinct expressions can have equal pairs (a, b) up to Gamma, and the pair is a linear map of the table coefficients. For two formally distinct expressions e, e' of one trail, the difference of their pairs is a nonzero linear form in the coefficients v = (a_S, b_S, c_1, d_1, ...). By Schwartz–Zippel, Pr[(a,b)_e = g(a,b)_e' for a fixed g in Gamma] ≤ 1/r over the algorithm's coins. After conditioning on the earlier history, which removes at most a fraction |Gamma| t/r of the coefficient space, the conditional bound is at most 2|Gamma|/r per pair (section 3).

### 2.8 Merges: Lemmas P and M

A *merge* is a visit by a walk to a class visited earlier by a walk of another trail (start and escape points included); a step into a class of the walk's own trail is an entry (section 4). Let τ be the first batch step with a merge. A *root segment* is a trail's path up to its first merge; by section 2.2 the fresh-source bounds of sections 2.3 to 5 are exact on root segments.

**Lemma P (continuation).** On the good event (no closure of length 65 or more, no within-trail coincidence), the classes a walk visits from a class C up to its first distinguished class depend only on C.

*Proof (outline).* The class, canon(C), the key (hence the adder index and the distinguished test) and the escape rule depend only on the class (section 1.2); the orientation does not enter the class sequence. A fruitless cycle is the cycle of the raw class map that contains the first repeated class. Every history reaches the same cycle set (Lemma R) and escapes to the same 2·canon(M), with M the least-key member. ∎

**Lemma M (useful merges).** A merge is *useful* if both trails it joins record a distinguished point. A merge is *non-recording* if either trail is abandoned or cut, so that it supplies no record. For a useful merge the two trails share their continuation (Lemma P) and record the same distinguished class, so the second record gives the key match. Let C_T be the number of merges by T = T_steps, and N_T the number of non-recording merges. Then

  P(no useful merge by T) ≤ P(C_T = 0) + E[N_T] ≤ exp(−tau(1 − eps)) + E[N_T],

  E[N_T] ≤ T·S_max·(1−θ)^(CAP/2)/N ≈ 2^−46,

where the first term is the birthday bound of section 7.

*Proof (outline).* If C_T ≥ 1 and no merge is useful, every merge is non-recording, so N_T ≥ 1; this gives the first inequality. For the expectation, a step merges onto a given visited class c with probability 1/N given the history (T1), and a class lies on a non-recording trail only if that trail's first CAP−64 steps carry no distinguished point. By T1 the probability of that for a fixed trail is at most (1−θ)^(CAP/2), since at least CAP/2 distinct classes are visited except with negligible probability (entries are bounded by Lemma E1*). Summing over the steps and the trails, whose total length is at most S_max, gives the bound. ∎

*Edge cases.* A merge into an abandoned trail is a non-recording merge, and it is counted in N_T. A merge into the middle of an escape sequence continues with the escape of its cycle, by Lemma P. The same key with another orientation is the same class, so Lemma P applies; the coefficients differ by a unit, which Lemma D covers. A restart draws fresh coins, so its start merges only into another trail. A merged walk that passes the end of a non-recording trail continues into unvisited classes and stays fresh there.

*Merged steps.* A merged walk takes ordinary steps on another trail's path, at most TRAIL_MAX of them, charged as any trail and in the merged row of section 6. The bounds of section 3 hold for every trail's own path, root or merged, since their proofs use only that the class data are independent of x and of the coins.

## 3. Coincidences: Lemma G (nondegenerate and degenerate), DELTA_LOOP and Lemma L

**Lemma G.** Let the next point of a trail have t earlier distinct classes, and let there be no earlier coincidence. Then

  Pr[the new point's class is among the t earlier classes | history] ≤ |Gamma| t / (r − |Gamma| t(t−1)/2) ≤ |Gamma| t / (r − |Gamma| TRAIL_MAX²/2).

*Proof.* On the event that there is no coincidence so far, the transcript is a function of the class data of distinct classes, which are independent of x by R. Conditioned on it, x is uniform on the complement of the roots of the earlier pairs' linear equations: at most |Gamma|·C(t,2) roots, one per Gamma element per earlier pair. The new point's coincidence roots are at most |Gamma| t. Since the coins are fixed and x is uniform, this holds for every coin outcome without degenerate coincidences, and so over the coins and x together. ∎

**Within-trail coincidences.** A within-trail coincidence is a pair of points of one trail that is either a nondegenerate coincidence (Lemma G) or a degenerate one (section 2.7). Summing over the points of every trail and adding the table-level term A² |Gamma| / r for the degenerate pairs of the table,

  DELTA_LOOP := K_TRAILS · (TRAIL_MAX²/2) · |Gamma| · (2/r + 1/(r − |Gamma|·TRAIL_MAX²/2)) + A²·|Gamma|/r ≈ 1.0503·10^−16

bounds the probability of any within-trail coincidence, over x and the coins.

**Cut-off and false filter matches.** A divergence between a live trail and its replay needs a cycle of length at least 65 that one anchor schedule detects and the other misses. A false 64-bit filter match alone cannot cause it, since an escape needs an exact check. Such a cycle is a within-trail repeat, covered by DELTA_CUT = DELTA_LOOP + DELTA_LONG ≈ 1.05·10^−16, where DELTA_LONG = S_max·Q_FAR bounds the probability of any closure of length 65 or more through unvisited classes (Lemma F, section 4; a closure through visited classes is a merge, section 2.8). Q_FAR < 2^−361, so DELTA_LONG < 2^−234. No per-trail term is charged for false matches.

**Lemma L (no re-entry into an escaped cycle).** By D2, every re-entry is a within-trail coincidence, so its probability is at most DELTA_LOOP. (A repeated escape point is exactly such a re-entry: the escape point depends only on the cycle's least-key member.)

## 4. Closure entries: Lemma E1 and Lemma E2

**Definitions.** A source is *fresh* if its class is unvisited by every walk (section 2.2). A step is an *entry* if its source is fresh and its successor lies on a class already visited in the same trail. A step whose source was visited earlier in its own trail is *persistent*. Only entries open a closed cycle, and every persistent step follows an entry. A step whose successor lies on a class visited by another trail is a merge (section 2.8).

**Lemma E1 (marginal).** Pr[step is an entry] ≤ p_c := 1/(|Gamma| A) + c/A² + Q_FAR + 2^−63, with c = 3, where Q_FAR bounds the closures of length 65 or more (Lemma F). This is an unconditional bound, used for expectations: the rare work of section 6 and the expected counts. For the conditional bound used in the tails, see Lemma E1* below. accounting.py evaluates Q_FAR and checks that log2 Q_FAR < −200.

*Proof.* The window of a closure is a run of ℓ distinct classes C_1..C_ℓ whose successor is C_1. Each C_p is fresh when visited (a window through a visited class is a merge, charged in section 6, or a within-trail coincidence, which the good event excludes). By the coupling of section 2.2, the index i_p = (K_C >> 1) mod A and the unit u_p of its step are i.i.d. uniform on A and Gamma, and every window rate below is an upper bound for the window's probability on fresh classes. No independence is assumed for selected windows. A formal closure needs, for every index a in the window, the units of its occurrences to sum to zero. An index occurring once cannot cancel, so every index has multiplicity at least 2 and the number D of distinct indices is at most ℓ/2. For an index of multiplicity m ≥ 2 with independent uniform units, the last unit must equal one specific value, so the zero-sum has probability at most 1/6. Independence over the distinct indices gives a factor at most 6^(−D), and the closing class's identity γ in Gamma is a conservative factor 6. The rate of windows of length ℓ is therefore at most

  rate_ℓ = 6 · Σ_{D=1}^{⌊ℓ/2⌋} S2(ℓ, D) · A(A−1)…(A−D+1) / A^ℓ · 6^(−D),

where S2(ℓ, D) counts the partitions of ℓ positions into D blocks of size at least 2. For A = 2^16, rate_3 = 1/A² exactly, rate_4 = 1/(2A²) + 1/(2A³), and Σ_{ℓ=5}^{64} rate_ℓ = 7.40·10^−15 = 2^−46.9. So the windows of length 3 to 64 together have rate at most (1.5 + 4·10^−5)/A² ≤ 1.51/A²; the 2-cycle is the separate exact rate 1/(6A) at point level. The accounting keeps c = 3, which is above this sum. The 2^−63 term covers false 63-bit filter matches, which the exact check rejects. ∎

**Lemma E1\* (uniform conditional bound).** Let F_{t−1} be the history before step t. Let G_{t−1} be the event that, before step t, no walk has merged, the trail has no coincidence and no closure of length 65 or more. If the source of step t is fresh, then on G_{t−1}

  Pr[step t is an entry | F_{t−1}] ≤ p* := 64/(|Gamma| A).

*Proof.* Let the source be C, with walk point W = f(prev) = h·canon(C), where h ∈ Gamma is the orientation of W. By R (T1, T2) on the global filtration, given F_{t−1} the index i = (K_C >> 1) mod A is uniform on the A adders, and h is uniform on Gamma, independently, since C is unvisited. Every earlier point V_j of the trail within distance 64 of C has a formal vector V_j = 2^{m_j} g_j S + r_j, and W = 2^m g S + r, where m, m_j are escape counts and the formal vectors are determined by F_{t−1}. The successor of C has class [canon(C) + R_i], so it closes to V_j iff canon(C) + R_i = γ V_j for some γ ∈ Gamma. Off the degenerate pairs of section 2.7 (in DELTA_LOOP), this is a formal identity. Comparing S-coefficients gives h^{−1}2^m g = γ 2^{m_j} g_j. If m ≠ m_j, then 2^{m−m_j} ∈ Gamma up to sign, which D1 excludes on G_{t−1} (since |m − m_j| ≤ TRAIL_MAX and section 2.4 covers it). So m = m_j, γ = h^{−1} g g_j^{−1}, and comparing the R-parts after multiplying by h gives

  h e_i = −D_j, where D_j := r − g g_j^{−1} r_j,

a vector determined by F_{t−1}. So a closing to V_j needs D_j = −c·e_{i_j} with one adder i_j and c ∈ Gamma, and then only the pair i = i_j and h = c works. Its probability is at most 1/(|Gamma| A) if D_j is a single adder vector, and 0 otherwise. Summing over the at most 64 earlier points gives the bound. (The predecessor always gives a single adder vector: it is the 2-cycle case.) ∎

**Lemma F (long closures, uniform over prefixes of unvisited classes).** Let R_s := max(63, 2s − 2) and

  Q_FAR := Σ_{s=1}^{A−1} A·C(A,s)·|Gamma|^{1−s}·(s/A)^{R_s}·A/(A−s) + TRAIL_MAX·A·|Gamma|^{1−A}.

For every prefix F in which the segment's classes are unvisited, the probability given F that the walk closes at distance d ≥ 64 from its source, with the segment fresh, is at most the i.i.d. bound B_{s,d} of the segment (the bound on the right of the displayed coupling inequality, which does not depend on F). Summing over the at most S_max sources, the probability that a run contains a closure of length 65 or more through unvisited classes is at most S_max·Q_FAR. Numerically log2 Q_FAR < −361 (the largest term is s = 33), so DELTA_LONG = S_max·Q_FAR < 2^−234. A closure through a visited class is a merge (section 2.8).

*Proof.* Let the closure close a fresh source C = c_d to V = c_0 with d ≥ 64, where c_0 is in F and c_1, …, c_{d−1} are new after F. V and C have equal escape counts, since a difference k in [1, TRAIL_MAX] would put 2^k in Gamma, which section 2.4 excludes. So the closing forces the net adder vector of the d steps from c_0 to c_{d−1} to be a single adder vector −h e_i with h ∈ Gamma, as in Lemma E1*. The step out of c_0 is fixed by F; the d−1 steps out of c_1, …, c_{d−1} are random. The segment has no escape, since an escape would change the escape count. Every index of the segment other than i occurs at least twice, with zero unit sum, since an index occurring once has nonzero net. Let s be the number of distinct indices. Since d ≥ 64 and at least 2s − 1 steps are needed for s distinct indices, the number of random steps satisfies d − 1 ≥ R_s = max(63, 2s − 2). For a fixed set S of s indices and i ∈ S, the probability that the random indices lie in S is (s/A)^{d−1}. Given the indices, the units are independent and uniform, so the zero sum of each of the s−1 indices of S other than i has probability at most 1/|Gamma| (the last unit is uniform), and these events are independent. Summing over S (at most C(A,s)), over i (at most A), over the lengths d−1 ≥ R_s (a geometric sum giving (s/A)^{R_s}·A/(A−s)), and over the starts, gives the bound. The case s = A has probability at most |Gamma|^{1−A} for each d, and there are at most TRAIL_MAX values of d. ∎

**Lemma FR (forced return; exact).** Let steps p and p+1 be fresh, with tau_p = −1, where tau_p is the canonical rotation of P'_p = canon(C_p) + R_{i_p}, and the transition law is u_{p+1} = u_p tau_p^{−1}. If i_{p+1} = i_p, then P'_{p+1} = −C_p, so tau_{p+1} = −1 and C_{p+2} = C_p.

*Proof.* tau_p = −1 gives C_{p+1} = canon(P'_p) = −P'_p. Then P'_{p+1} = C_{p+1} + R_{i_{p+1}} = −P'_p + R_{i_p} = −C_p. C_p is canonical, so canon(−C_p) = C_p, which gives tau_{p+1} = −1 and C_{p+2} = C_p. ∎

Let P = Pr[i_{p+1} = i_p | tau_p = −1, both steps fresh], and let rho be the non-equal share, Pr[tau_{p+1} = −1 | tau_p = −1, both fresh, i_{p+1} ≠ i_p]. Lemma FR gives Pr[tau_{p+1} = −1 | tau_p = −1, both fresh] = P + (1 − P)·rho, so the excess over 1/6 is

  E = (5/6)·P + (1 − P)·(rho − 1/6).

With rho = 1/6 (T2), E = (5/6)·P exactly. Under index uniformity P = 1/A, so E = 5/(6A). The forced return has rate P/6 = 1/(|Gamma| A), which is the 2-cycle rate of Lemma E1.

*Measurements* (evidence/sign-law.json, from evidence/sign_law.py; small builds, A = 8 to 1024):
- Equal-index rate times A: 0.967, 0.992, 1.015, 0.985 (each with a 95% interval containing 1). All 10,219 equal-index events had tau_{p+1} = −1, as Lemma FR requires.
- Non-equal share rho: 0.16417 (A = 8, z = −0.79), 0.16638 (A = 64, z = −0.39), 0.16671 (A = 256, z = 0.08), 0.16667 (A = 1024, z = 0.03).
- Share-based excess times A: 0.788, 0.808, 0.856, 0.828, against 5/6 = 0.833.
- The seeds 1 and 2 at A = 1024 each show 874 equal-index pairs; both were recounted, and the runs are independent (chance about 1%). This is recorded, not treated as a signal.

**Persistent steps.** Inside the window of an entry (from the entry to its escape) the walk stays on the visited cycle, so there is no fresh source. Persistent steps cost the common charge; the checks that resolve an entry are charged to the entry (W_ENTRY).

**Lemma E2 (persistent steps).** On the good event of section 3 (which includes every cycle of length 65 or more through the cut-off of DELTA_CUT), each entry is followed by at most 128 persistent steps before the escape that resolves it, deterministically.

*Proof.* A 2-cycle entry is detected at the entry step itself, so it has no persistent steps. A longer cycle of length l ≤ 64 is entered at the entry step. Within 64 steps the walk reaches an anchor step, the anchor point lies on the cycle, and the walk returns to the anchor's class within l further steps, where the anchor check fires and the escape follows. So there are at most 64 + l ≤ 128 persistent steps, and processed steps ≤ (fresh-source steps) + 128 · (entries). ∎

## 5. Abandonment: Lemma T'' (first moment over entries)

Let n = 64 (CAP_PERIODS + 1) = TRAIL_MAX be the processed steps before a trail is abandoned, with CAP_PERIODS = 2^DP_BITS, so n·theta is about 64. Let E_n be the number of entries among the first n steps of a root segment (section 2.8). On the trail's coincidence-free event G_i, the fresh-source steps number at least n − 128 E_n (Lemma E2). Each fresh-source step is a distinguished point with probability theta, independently of the past, because its key K_C is uniform (R; the DP test reads only K_C). A persistent step is never a distinguished point, since its class was visited earlier and a distinguished class would have ended the trail. So with no distinguished point among the fresh-source steps, the trail is abandoned, with probability at most (1 − theta)^(n − 128 E_n).

*First moment over entries (chain rule with Lemma E1\*).* Let G be the good event of section 3 (no coincidence and no closure of length 65 or more in the root segment, and no merge before its end). By Lemma E1\*, on G each step is an entry with conditional probability at most p* = 64/(|Gamma| A) given the history. So for any k steps, the probability that all of them are entries and G holds is at most p*^k, and Pr[E_n ≥ k and G] ≤ binom(n, k)·p*^k ≤ (e·n·p*/k)^k. Take k = n/c with c = 2^10, so k = 2^56. Then e·c·p* = e/6 < 1, and 128/c = 1/8, so on E_n ≤ k at least 7n/8 steps are fresh, and

  Pr[abandon_i and G_i] ≤ (e/6)^(2^56) + (1 − theta)^(7n/8) ≤ (e/6)^(2^56) + exp(−56) ≤ 4.79·10^−25.

*Failure form.* Let X be the number of trails abandoned before the first merge τ. The run fails by abandonment before τ only if G fails or X > K, so Pr[abandonment failure] ≤ Pr[not G] + E[X·1_G] / K. With E[X·1_G] ≤ K_TRAILS · 4.79·10^−25, DELTA_ABANDON = K_TRAILS · 4.79·10^−25 / K = 2.76·10^−14, and Pr[not G] ≤ DELTA_CUT (section 3).

## 6. Expected rare work: Lemma M'' (each entry charged by its type)

Per processed step, the expected metered rare work is at most E_RARE, with each entry charged by its type:

| event | rate | worst case (instructions) |
| --- | ---: | ---: |
| 2-cycle entry (true or false predecessor match) | 1/(|Gamma| A) per fresh source (Lemma E1) | 48,039 |
| longer entry, 3..64 (anchor check, with its look-ahead steps) | c/A² per fresh source, c = 3 | 1,618,122 |
| false 63-bit filter match (a false first-filter match runs the exact check; a false anchor match runs the look-ahead; both may occur) | 2^−63 per processed step | 48,039 + 1,618,122 |
| exact-k fallback, per canon call | 2^−55 (R, for the canonical coordinate fractions and the rotation) | 256 |
| merged step (a walk on another trail's class, section 2.8), per processed step; tau·TRAIL_MAX merged steps in expectation | tau·TRAIL_MAX / S_max | 1,618,122 (an anchor check) |
| within-trail coincidence, or closure of length 65 or more (sections 3, 4), per step | |Gamma|·TRAIL_MAX·(2/r + 1/(r − |Gamma|·TRAIL_MAX²/2)) + Q_FAR (below 2^−361) | TRAIL_MAX·(821 + reserve) |

Persistent steps have no entry charge of their own: they are ordinary processed steps at the common cost C = 821, covered by S_max.

The fallback row uses R (fractions and rotation, T1 and T2): the test of Lemma 1 fails with probability at most 204/2^64 < 2^−55.6 for uniform fractions. The charge is the whole canon call with the fallback taken, 256 (static bound 253), per canon call.

The within-trail coincidence row is self-referential, because its worst case contains the reserve. The fixed point is E_RARE = (E_FIXED + rate·TRAIL_MAX·C_COMMON) / (1 − 64·rate·TRAIL_MAX), with E_FIXED the first five rows: E_FIXED = 0.12329973978984803 (the fallback row adds 7.1·10^−15, the merged row 4.4·10^−13), E_RARE = 0.12329973978984803 (the coincidence row adds 1.4·10^−27). The reserve per processed step is 64·E_RARE = 7.891183, and the charge per walk step is C_COMMON + reserve = 828.891183.

Markov's inequality gives Pr[metered rare work > 64·S_max·E_RARE] ≤ 1/64, counted in delta, with the expectation taken in R. The cap is hard: the meter aborts beyond it, so the work bound holds on every run that does not abort.

## 7. Success probability

Let T_steps = K·L. Then tau = T_steps²/(2N) = 0.530300. Under R, the probability that the first T_steps steps contain no merge (section 2.8) is at most exp(−tau·(1 − eps)), with eps = 4/A + 2·(3/A) + 2·(3·192/A²) = 0.00015286: the adder penalty (4/A), 2-cycles at a rate of at most 1/A with 3 non-fresh states each, and short cycles at a rate of at most 3/A² (section 4) with 192 wasted states each.

By Lemma M, the probability that no useful merge occurs by T_steps is at most exp(−tau(1 − eps)) + E[N_T]. If a useful merge occurs and no failure occurs before it, the two trails record the same distinguished class, the second record triggers resolution, Lemma R makes the replays reach the recorded points, and Lemma D makes the equation degenerate with probability exactly 1/r. Fresh start points only add points, and restarts at distinguished points give each record its own coins.

The failure terms are:
- **Within-trail coincidences and cut-off replays (sections 1.3 and 3)**, including re-entry into an escaped cycle (Lemma L), degenerate pairs and cycles of length at least 65: at most DELTA_CUT = DELTA_LOOP + DELTA_LONG ≈ 1.05·10^−16 (this is Pr[not G]).
- **Abandoned trails before τ (section 5)**: at most DELTA_ABANDON = 2.76·10^−14.
- **Non-recording merges (section 2.8, Lemma M)**: at most DELTA_NR ≈ 2^−46.
- **Degenerate first match:** exactly 1/r (Lemma D).
- **Work meter abort (section 6):** at most 1/64.
- **Resource aborts (sections 1.3 and 8.2):** sampling caps at most DELTA_CAP < 2^−440, and the probe cap at most DELTA_PROBE ≈ 2^−146.
- **DP_MAX overflow:** by Chernoff under R, at most exp(−DP_MAX/12), below 2^−60.
- **Exceptional additions:** at most 2^−110 under R.
- **Other terms:** together at most 2^−20 = delta.

Therefore:

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015286)) − DELTA_CUT − DELTA_ABANDON − DELTA_NR − DELTA_CAP − DELTA_PROBE − 1/64 − 2^−20 = 0.395898 ≥ 0.39.**

The probability is over the uniform target x and the algorithm's coins. The model term is 0.411524. The budget is 1.0299·√N steps, against 0.9943·√N for success exactly 0.39. The 0.0215 of slack in the model term covers the sampling error of the toy evidence below.

**Evidence for the model and the procedure.** `experiment.py` runs the procedure of section 1.3 on random prime-order toy curves y² = x³ + 7 (p ≡ 1 mod 3), with the same class key, canonical rule, adder images and cycle rules, fresh-coin trails with the cap, and replay resolution. It keeps no oracle collision table; a trial succeeds only if the first key match resolves to a verified x.

<!-- toy-evidence -->
Offline runs (`evidence/trails/`, pooled in `summary-double.json`): 5,200 trials on 5 curves with N = 2^26.8, 2^29.1, 2^30.8, 2^31.3, 2^32.5 classes, adder tables of 2^11 and 2^16, 16 slots and theta = 2^−5, escape "double".
- Success of the whole procedure at tau = 0.5303: 2,169 of 5,200, 0.4171, Wilson 95% interval [0.4038, 0.4306]. The model 1 − exp(−tau) is 0.412. Every success was verified ([x]P = Q). Degenerate first matches: 0. Abandoned trails: 0. Longer fruitless cycles escaped: 6. The finishing phase added at most 2.6% to the K·L steps.
- Production 64-bit filters: every recorded trail was replayed with a different anchor schedule, 4,619,736 replays, 0 that missed the recorded point.
- Lemma R stress (key `trails-stress-double` in `evidence/stress-and-sensitivity-double.json`): 8-bit filters over 5,200 trials, so false filter matches were frequent. 4,837,229 trails were replayed with a different anchor schedule, and 0 missed the recorded point; 349,511 distinguished points were first reached during a look-ahead. These runs' success rates are not used, since each such look-ahead adds up to 65 points that the step budget does not count.
- Replay-order sensitivity (key `lemma-r-sensitivity-double`): 200 trials per curve with 6-bit filters and 16 adders. With the current order, 0 of 59,977 replays (j0) and 0 of 95,625 (a3) missed the recorded point. With the earlier order (an offline switch), 8 of 59,351 (j0) and 47 of 92,993 (a3) did.
- 2-cycle entry rate per step: 7.25e-05, 7.74e-05, 7.73e-05, 7.76e-05 with A = 2^11 (1/(6A) = 8.14e-05); 2.03e-06 with A = 2^16 (1/(6A) = 2.54e-06).
<!-- /toy-evidence -->

The declared initial input (`bits` 28, `b` 7, 80 trials, seed 2801, production 64-bit filters, every trail replayed, escape "double"; `evidence/declared-run-double.json`) runs in a few seconds. It reports 35/80 = 0.4375, Wilson 95% interval [0.334, 0.547], against the model 0.412, with 0 degenerate matches, 0 distinguished points reached during a look-ahead and 0 replays, of 15,059, that missed their recorded point. With 80 trials the standard error is about 0.05, so it is a reproducibility check; the success claim rests on the pooled offline trials above.

Limits: the toy curves have N up to about 2^33, against 2^253 in production, and theta = 2^−5 there, so start points and the finishing phase are a few percent of the steps. R is an idealization in N, not proved for this walk.

**Premises, verified facts and what is not established.**
- **Premises.** R (section 2.2) is the one heuristic premise. H3 (section 8) is a measurement assumption about instruction counts. The coins are the algorithm's own randomness, not a premise.
- **Derived from R and exact algebra.** The coincidence bounds (section 3), the entry bounds (section 4), the abandonment bound (section 5), the useful-merge bound (section 2.8) and the metered work bound (section 6).
- **Verified here.** The escape-depth certificate (section 2.4, by the assertions in accounting.py), Lemma D (exact over the coins), the forced-path maxima of section 8.1 on the bounding builds, and the accounting.
- **Not established.** That real class keys and orientations behave as R at production size. The measured sign law is marginal only (section 4), and the toy and formal runs use theta = 2^−5 and 2^−6.

## 8. Computation

**Charge per walk step.** The common batched path is charged C_COMMON = 821 instructions per walk step: the larger of the plain step (808.67) and the anchor step (818.66), each measured as a whole-run instruction difference on production inputs, plus 2 for the canon_pick branch, rounded up. Rare events are charged their worst cases (section 6, and section 8.1 for the primitives), and their expected counts by the Markov reserve 64·E_RARE = 7.891183 per walk step (Lemma M''). The charge per walk step is therefore 828.891183 instructions.

**Cost is uniform over the coin values.** The common path is branch-free apart from canon_pick and the exact-k fallback (metered at its forced maximum), so each walk step is bounded uniformly over coin outcomes, and the ledger charges that bound.

**Field arithmetic.** `measurement/secp_field.h` implements multiplication and squaring modulo p = 2^256 − c, c = 2^32 + 977, as GCC inline AArch64 assembly, with the final carry selecting between y and y − c. `measurement/secp_sub.h` holds the subtractions. The routines matched Python big-integer arithmetic on about 194,000 reduced and 199,000 arbitrary 256-bit operand pairs, including constructed cases for the rare second-fold carry. The replica check below compares every final point of 2^16 walks after 2^19 steps with Python. Callgrind counts this inline assembly exactly as it counts compiled C.

**8.1 Worst cases of the rare-path operations.** Each rare event is charged a forced-path maximum, not a sampled average. For each primitive (f_single, add_single, escape, same_class, point_is_dp, canon_k with its exact fallback) the data-dependent branches were enumerated, and the control variants were forced by chosen inputs with the maximum taken; a bounding build with branch-free routines (checked for equality on random and edge inputs) supplies the measurements. Measured maxima (Callgrind Ir per call, organizer image and flags, `measurement/forced/summary.json`):

| primitive | measured maximum (Ir) | charged bound (Ir) |
| --- | ---: | ---: |
| f_single (adder selection and addition) | 23,345 | 23,354 |
| add_single | 23,079 | 23,079 |
| same_class | 70 | 70 |
| point_is_dp | 41 | 41 |
| escape of a 2-cycle (canon, then doubling) | 23,534 | 23,550 |
| escape of a 65-cycle (canon, then doubling) | 26,432 | 26,448 |
| canon_k_exact (fallback body, including its call) | 188 | 197 (static bound) |
| canon_k with the fallback forced | 244 | 253 (static bound), charged 256 |

A 2-cycle check costs at most f_single + same_class + escape(2) + point_is_dp plus a control overhead of at most 1,024 instructions, and an anchor check costs at most f_single + escape(65) plus, for each of its at most 64 look-ahead steps, f_single + point_is_dp + same_class, with the same overhead per step. The harness figures are conservative: in its own translation unit GCC compiles finv about 36 instructions more expensive. As an independent check, the longest control-flow path of the unmodified binary, with every conditional branch free and finv at its measured cost, gives at most 47,193 instructions for a 2-cycle check and 1,547,316 for an anchor check, both below the charges above (`measurement/forced/summary.json`). Scalar multiplications and replays use f_single's worst case as their unit, an upper bound for an addition or a doubling.

**Measurement (H3).** `measurement/step.c` runs the trails, steps, cycle handling, distinguished points, restarts and trail cap. Its common batched path was measured on production inputs as whole-run Callgrind `Ir` totals (GCC 14.2.0, organizer's image, `-O3 -fno-tree-vectorize -mgeneral-regs-only`), with the rare paths out of line (`evidence/measurement/production-totals.json`). The anchor step is 10.0 instructions per walk step above the plain step, and its total contains its own rare events.

- Plain step (steps 1 and 2), per walk step: 808.67 instructions (measured).
- Anchor step (steps 63 and 64), per walk step: 818.66 instructions (measured).
- C_COMMON = 821 per walk step: the measured anchor step plus 2 for the canon_pick branch, rounded up. The bound for walks that are all long is 819.34, so 821 is conservative.

**Rare-path checks.** Distinguished points and false 64-bit filter matches do not occur in a production-size sample, so small builds of the same source exercise them (`evidence/measurement/trail-paths.json`, 2^6 adders, 2^12 walks): DP_BITS = 4 (restarts, 2-cycles), DP_BITS = 12 with a one-period cap (abandonments), and DP_BITS = 6 with 8-bit filters (longer cycles, and distinguished points first reached in a look-ahead). In every build the replica matched every final point, event count and a digest of every record. A mutant that omits look-ahead records agrees on every count and differs in the digest (`evidence/measurement/record-mutation-double.json`).

**8.2 Sampling caps and the distinguished-point table.** The table has C = 2·DP_MAX slots of 128 bytes, with open addressing, linear probing and index (key(W) >> 77) mod C; the high bits of a distinguished key are uniform, so the index is uniform up to a bias of about C/2^179. With at most DP_MAX records the load is at most 1/2. A lookup or insert examines slots one at a time, stopping at an empty slot (a new record) or at a key match (resolution); if it examines L = 2^10 occupied slots, the run aborts.

*Probe tail.* Suppose an operation examines L occupied slots from h, and let [s, s+ℓ) be the maximal occupied run containing h, with ℓ ≥ L. Slot s−1 has always been empty, so every record in the run hashes into [s, s+ℓ): the window holds ℓ records. By T1 the keys of successive records are uniform and independent given the history, so the number of records hashing into a window of length ℓ is binomial with mean at most αℓ, α = DP_MAX/C ≤ 1/2. Chernoff gives Pr[≥ ℓ] ≤ (eα)^ℓ e^(−αℓ) ≤ e^(−cℓ) with c = ln 2 − 1/2 ≈ 0.1931. A union bound over the DP_MAX operations, the C window starts and all ℓ ≥ L gives

  DELTA_PROBE := DP_MAX·C·e^(−cL)/(1 − e^(−c)) ≈ 2^−146.

*Charge.* A table operation is charged L probes at PROBE_COST = 2^11 per probe, so DP_COST = 2^21 (raised from 2^20); the allowance covers the slot fetch (a memory access), the comparison and the control, and is not a measurement. A draw costs DRAW_COST = 2^16, and each start is charged its worst case, START = D(2·SCALAR_MUL + 2D·2^16) + 2^14 = 145,600,512 instructions, with SCALAR_MUL = 768·23,354; adders and initial starts are charged at START each. The aborts add DELTA_CAP and DELTA_PROBE to delta (section 7).

**Total.** S_max = K·(L + TRAIL_MAX) ≤ 2^126.749959 walk steps.

| Component | Bound (units) | log2 |
| --- | ---: | ---: |
| walk steps | 828.891183·S_max / C | 126.479214 |
| adder table | 2^16 × START / C (START of section 8.2) | 33.15 |
| initial walk starts | K × START / C | 49.15 |
| trail restarts | (DP_MAX + K) × START / C | 84.90 |
| distinguished points | DP_MAX = 2^67.75 entries × 2^21 instructions / C | 78.78 |
| collision replay | 2 trails × TRAIL_MAX steps × 2^13 unbatched ops × (2^15 + 2^10) / C | 84.61 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 15.68 |
| input validation | 2^24 instructions | 14.03 |

Rounded up, **time_log2 = 126.479214**. Preprocessing (adders and initial walk starts) is included. All processors' work is summed, and elapsed time plays no part. The replay row is charged in every run, with 2^13 operations per replayed step (section 1.3, step 7).

## 9. Memory and advice

Peak memory is at most the sum of: 2·DP_MAX·128 bytes for the open-addressing table (C = 2·DP_MAX slots of 128 bytes, section 8.2); K·576 bytes of walk state; 2^16·6·64 bytes of adder images; and 16 MiB. That is 2^75.74996 bytes. Nonuniform advice is 0 bytes.

## 10. Limitations

- R is a single heuristic premise: the random-function idealization of the adder walk on classes. The model term is 0.4115, 0.0215 above 0.39. Its collision count is the standard random-mapping count (section 2.1).
- The argument is not established at production size: the toy runs use theta = 2^−5, and the formal test uses 2^6 adders and theta = 2^−6.
- The exact-k fallback rate (2^−55 per canon call) follows from R and is not measured, since no fallback occurred in the measured steps.
- The measured instances use seeds (seed 1 for `reproduce.sh`, seed 2801 for the declared run, a seeded RNG for the toy experiment); the algorithm draws its coins from the operating system. The toy experiment retries identity starts without a cap; the production procedure caps them (section 1.3).
- Instruction counts come from one implementation and compiler (section 8). Rare operations are charged forced-path maxima (section 8.1), checked on the forced paths and not proved for every input. The field routines were tested, not formally verified. The probe cost (2^11 per probe) and the draw cost (2^16 per draw) are allowances, not measurements.
- The rare-work meter is taken in R. A model-free bound would charge the anchor check at its worst on every step, which costs about 10.9 bits (time_log2 about 137.4 instead of 126.5); we do not take it.
- About 2^67.7 distinguished points and 2^126.5 work are needed. This is not a practical attack and does not weaken secp256k1.

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
python3 -B -I submissions/secp256k1-ecdlp/rho-sigma6-asm/evidence/sign_law.py <dir>   # writes <dir>/sign-law.json; cmp with evidence/sign-law.json
echo '{"bits":28,"trials":80,"seed":2801,"b":7,"adder_bits":11,"dp_bits":5,"walks":16,"tau_e4":5303,"filter_bits":64,"check_all_replays":1,"escape":"double"}' | python3 -I submissions/secp256k1-ecdlp/rho-sigma6-asm/experiment.py
```

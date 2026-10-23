# secp256k1-ecdlp: parallel Pollard rho on order-6 classes with a sigma walk, all arithmetic in instructions

Claim: for uniform x in [0,r) and Q=[x]P on secp256k1, the algorithm below returns x with probability at least 0.39. It uses at most **2^126.479214** target-group-operation equivalents in total (score 63239607/500000), including all preprocessing, failures, lookup and verification. Peak memory is at most 2^75.749960 bytes. Preprocessing is 2^47.131242 units, and no nonuniform advice is used. This is a classical cold attack on one target.

The success bound rests on the disclosed heuristic H1 (the random-mapping model for r-adding walks, used for the collision statistics), on the premise M1 (equivariant random orientation, section 2.2) with its orientation law L-ori, on the randomness of the algorithm's coins, on the exact certificate of section 2.4 and on the depth tail of section 2.3. The 2-cycle excess is an exact forced-return identity (Lemma FR, section 4); L-ori and the non-equal share are measured premises. The cost rests on H3, average instruction counts measured with the organizer's toolchain, and every run is charged by a work meter with a Markov reserve. This is a carefully accounted generic attack. It is not a practical break.

**Cost route.** The walk does not use the 6M reference operation. Every instruction of every step, including all field multiplications, squarings and subtractions, is its own 4×64-bit GCC inline AArch64 assembly and is charged at 1/C with C = 1000 (H = 0 in COST_MODEL's H + R/C; "Faster arithmetic and lookup algorithms are welcome; justify their actual work"). A walk step is charged 828.89 instructions, 0.829 units.

The walk never moves a point to its canonical representative, so it never pays the multiplication by beta. It stores the actual point W and adds sigma(R_i), where sigma maps the canonical point to W. The class sequence is unchanged, so the standard analysis applies. Which rotation is canonical is decided from 64-bit fixed-point fractions with an exact fallback, and Lemma 1 bounds their error. The budget is tau = 0.5303.

**Distinguished points and resolution.** Walks run in trails. Each trail starts from fresh coins and ends at its first distinguished point, which it records once with its coins. A key match is resolved by replaying the two trails (section 1.3, Lemmas R and D), so a fruitless cycle cannot store the same point twice, and a degenerate match has probability exactly 1/r.

Files: `README.md` lists the package. Evidence paths appear where the text cites them.

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
- K0 = 2^50: the escape-depth bound of section 2.3.

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
2. **Adders.** For each adder i < 2^16: draw c_i, d_i uniform in [0,r) from the operating system's randomness, fresh in every run and independent of x; set R_i = [c_i]P + [d_i]Q; tabulate its six images.
3. **Trails.** Each of the K walk slots runs a sequence of trails.
   - A trail starts from fresh coins a, b uniform in [0,r), drawn from the operating system's randomness at every trail, including every restart, with W = [a]P + [b]Q (two scalar multiplications, charged in section 8). The slot keeps the coins until the trail ends. Coefficients are not tracked while walking.
   - A trail ends at its first distinguished point (step 6), or is abandoned when it crosses more than 2^60 anchor boundaries (step 5) without one, that is after between CAP − 64 and CAP = 64/theta steps, or when an addition gives O. Either way the slot starts a new trail from fresh coins. Every trail has fewer than TRAIL_MAX steps.
   - Each slot takes L steps (the main phase), then finishes its current trail and stops. The run fails if more than K trails are abandoned.
4. **Step.** W ← f(W), as one batched incomplete affine addition: Montgomery's simultaneous inversion over the 2^16 walks of a batch (3M), then lambda, lambda^2 and y3 (2M + 1S), with one Fermat inversion per batch. If the batch product is zero, the affected walk uses the unbatched complete formula. The canonical rotation of each new point is decided by the rule of section 1.2.
5. **Fruitless cycles.** Each walk keeps two filters, both the low 64 bits of a key: the key of the point before its current point W, and an anchor, its key at the last batch step divisible by 64. Each walk is anchored once every 64 of its own steps. The filters can match falsely; they only decide when an exact check runs.
   - *Order.* The distinguished-point test of step 6 comes first, on every point the walk produces: the new point, every look-ahead point and every escape point. A distinguished new point ends the trail, and no cycle check runs.
   - *2-cycles.* If the new key matches the first filter, the walk evaluates f(new) unbatched. If key(f(new)) = key(W) exactly, {W, new} is a 2-cycle of the class walk. Both points have been tested, and the walk replaces new by 2·canon(M), where M is the element with the smaller key.
   - *3..64-cycles.* Unless a 2-cycle escape just happened, if the new key matches the anchor filter, the walk iterates f from the new point, unbatched, up to 64 times, testing each look-ahead point first. A distinguished look-ahead point s ends the trail at s, which is the trail's record, since the walk would reach it next. If the new point's class returns, every point of that cycle has been tested, and the walk escapes to 2·canon(M), with M the cycle element of least key.
   - The anchor check also runs after a first-filter match that the exact check did not confirm. So a walk in a cycle of length ℓ ≤ 64 leaves it within 64 + ℓ steps, unless the cycle contains a distinguished point, which ends the trail first.
6. **Distinguished points.** W is distinguished when bits 17..76 of key(W) are zero. The test is applied to every point that a step or an escape produces. At a distinguished point the trail ends: the record (key(W), a, b), with the trail's coins, goes to a hash table keyed by key(W), and the slot starts a new trail. A trail makes at most one record. The run fails if more than DP_MAX records arrive.
7. **Resolution.** The first record whose key is already in the table ends the run. No other key match is used.
   - Both trails are replayed, unbatched, from their recorded coins up to their first distinguished point, with the rules of steps 4–6 and anchors every 64 of the trail's own steps. The replay tracks (a, b) mod r: an image sigma(R_i) has coefficients multiplied by (−1)^s·lambda^(−k), canon multiplies them by (−1)^s·lambda^k, and an escape doubles the coefficients of canon(M). By Lemma R each replay reaches the recorded point.
   - The replays give the common canonical point as [a1]P + [b1]Q = [a2]P + [b2]Q. If b1 ≠ b2, output x = (a2 − a1)/(b1 − b2) mod r after checking [x]P = Q. If b1 = b2 the run fails. By Lemma D this has probability 1/r.
   - **Replay bound.** A replay of a trail with n live steps is cut off after 65n steps (REPLAY_FACTOR = 65). Each live step gives at most one walk step and one escape of at most 64 cycle steps, and a replay step uses at most 70 unbatched operations, so 65·70·n ≤ 2^13·n. A cut-off replay fails the run, counted in DELTA_CUT (section 3). The replay charge is always made (section 8).
8. **Work meter.** Every run is charged at most W_max = S_max·828.891183 + G instructions, with G the fixed global charges (initial starts, trail restarts, distinguished-point records, solve and verify, input validation, adder table, resolution replays). Every walk step is charged C_COMMON = 821 instructions, and each rare event its worst case (section 6). The meter aborts the run when the total metered rare work exceeds 64 times its expectation bound (Lemma M''), with probability at most 1/64, counted in delta. The run also fails if more than DP_MAX records arrive or more than K trails are abandoned.

**Escape.** A cycle's escape goes to 2·canon(M), where M is the cycle's least-key member. canon(M) is a point of M's class, so 2·canon(M) is a point of the class of 2M. The escape uses no table and no stored point: it is one point doubling, charged at its forced-path maximum (section 8.1). Its class is not one the trail has produced before on the good event of section 3 (Lemma G), and its class data are uniform by M1.

**Escape repeats.** The escape point of a cycle depends only on its least-key member, so a repeated escape point means that the trail has re-entered a cycle it already escaped from. Section 3 bounds this event (Lemma L).

**Merged trails stay merged.** Outside fruitless cycles the next class depends only on the current class. Escapes exit to 2·canon(M), which depends only on the cycle's set of classes, not on where the walk entered the cycle or when it noticed. So two trails that reach the same class follow the same classes afterwards and reach the same first distinguished class.

**Lemma R (replays reach the recorded point).** The first distinguished point of a trail, and every escape before it, depend only on the trail's start point. They do not depend on when the filters fire or on the anchor schedule.

*Proof.* An escape happens only after an exact check shows that the point's class lies on a closed cycle C of length ℓ ≤ 65, and after every point of C has been tested and found not distinguished (the 2-cycle rule tests both points of {W, new}; the anchor rule tests each look-ahead point). Without the escape the walk would repeat C forever without a distinguished point, and 2·canon(M) depends only on C's set of classes. A look-ahead that meets a distinguished point ends the trail at the first distinguished point that f alone would reach next (a 2-cycle cannot lie strictly inside a longer cycle, so no escape could intervene). A false filter match either finds no closed cycle and changes nothing, or finds a genuine one and is handled as above. Every cycle of length ℓ ≤ 64 is caught: within 64 steps of entering it the walk is anchored at one of its points, and it returns to that class within ℓ more steps, where the anchor check runs. So under any schedule that anchors each walk at least once every 64 of its steps, the trail meets the same cycles, escapes to the same points and stops at the same first distinguished point. Every coefficient pair the replay carries represents its current point. ∎

The replay length is bounded as in step 7 of section 1.3, deterministically on every run.

**Lemma D (degenerate matches).** The first key match has b1 = b2 with probability exactly 1/r.

*Proof.* A trail records once, so the two records come from different trails t1 ≠ t2. Condition on all start points and adders. Given its start point W = [a + b·x]P, a trail's coin b is uniform on [0, r), because each b determines exactly one a. The coins are drawn from the operating system's randomness (section 2.2), so this holds exactly over the algorithm's coins, and coins of different trails are independent. Along a replay, b changes as b ← u·b + d, where u is a unit (order-6 classes, up to sign) or twice one, and u, d depend only on the path. So the canonical b at the distinguished point is U_i·b_ti + D_i with U_i invertible mod r, and P[b1 = b2] = 1/r. ∎

## 2. Premises, formal structure and checks

### 2.1 Heuristic H1

H1 (heuristic, disclosed). The class map [W] → [canon(W) + R_h([W])], with uniform random adders, behaves like a random mapping on the N classes, apart from the r-adding penalty, a factor of about (1 − 1/(6A)) in the birthday exponent (Teske 2001; Bernstein–Lange 2012; Bos–Kleinjung–Lenstra 2010). H1 is used only for the collision statistics of section 7, not for the within-trail statements of sections 2 to 6.

### 2.2 Premises: M1 (equivariant random orientation), the orientation law L-ori, and the algorithm's randomness

**M1 (equivariant random orientation).** For each class C of the walk:
- a uniform key K_C, which fixes the adder index (key >> 1) mod A, the 64-bit filter bits and the DP test;
- the canonical representative canon(C), a function of C alone, and the coordinate fractions of canon(C) used by the exact-k test, which are uniform;
- a uniform orientation g_C in Gamma, which records where the walk's representative of C sits relative to canon(C).

For a point W in C, the rotation k(W) and the sign s(W) are determined by W's relation to canon(C) and by g_C. The keys, orientations and canonical representatives are independent across classes and of x and of the coins. The rotation and sign depend on the representative and not on the class, which is why M1 is stated through g_C.

**L-ori (orientation law; a measured premise).** When the walk arrives at a fresh class C, g_C is uniform on Gamma given the history, so the arriving point's (k, s) is uniform on Gamma. At a revisited class, (k, s) is determined by the history. The consequence used in section 4 is the non-equal share: among fresh steps with tau_p = −1 whose next index differs, the next tau is −1 with probability 1/6. Section 4 measures it (within 0.8 standard errors of 1/6 at A = 8, 64, 256 and 1024) and does not derive it at production size.

**Index uniformity.** Under M1 the index of a fresh class, (K_C >> 1) mod A, is uniform over the A adders and independent of the history. It is a consequence of M1, not a separate premise, and it is not measured at A = 2^16.

**Randomness of the algorithm (not a heuristic).** The adder coefficients and the start coefficients are drawn uniformly from the operating system's randomness (`secrets.randbelow(r)` or equivalent), fresh in every run and independent of the target. Every trail, including every restart, draws fresh coins (a, b). Every probability in sections 3 to 7 is over x, the coins, and the class data (M1 and L-ori); nothing else is assumed.

### 2.3 Formal structure, escape depth, D1 and D2

Within one trail with start S, each walk point has the formal form W = 2^m g S + sum_i n_i R_i, with g in Gamma, integers n_i, and m the number of escapes so far in the trail. Steps add sigma(R_i) (which is ±lambda^j R_i) and do not change the S-coefficient. An escape replaces M by 2 canon(M) = 2 g' M, which doubles the S-coefficient.

**Escape depth (through entries).** Every escape resolves an entry (section 4), so the escapes of a trail are at most the entries among its steps. The first moment of section 5 gives, for a trail of n = TRAIL_MAX steps, Pr[a trail has K0 or more escapes] ≤ Pr[E_n ≥ K0] ≤ (e·n·p_c / K0)^K0. Here e·n·p_c/K0 = 0.45317 ≤ 1/2, so the bound is at most 0.45317^(2^50), that is 2^(−1.29·10^15). Over all trails, DELTA_DEPTH = K_TRAILS · 0.45317^(2^50) = 2^(−1.2856·10^15), counted in delta (section 7). The **depth event** is that every trail has fewer than K0 escapes. The disjoint-cycle bound TRAIL_MAX/2 is not used: an escape point can start the next 2-cycle.

**D1 (depth).** On the depth event and the good event of section 3, if E is the (m+1)-th escape point of a trail, its S-coefficient 2^(m+1) g' differs from the S-coefficient 2^j g'' of every earlier point of the trail (j ≤ m).
*Proof.* Equality would give 2^k in Gamma for k = m+1−j, with 1 ≤ k ≤ m+1 ≤ K0 − 1 < 2^50. The certificate of section 2.4 excludes every such k. ∎

**D2 (re-entry is not structural).** A re-entry into an escaped cycle is a later point (depth ≥ m+1) in the classes of an earlier cycle (depth ≤ m). By D1 such a point is never formally equal to an earlier one, so every re-entry is a within-trail coincidence, which section 3 bounds.

### 2.4 Exact certificate for D1: no k below 2^50 with 2^k in Gamma

**Claim.** For the secp256k1 order r and Gamma = {1, −1, λ, −λ, λ², −λ²} ⊂ Z_r, there is no k with 1 ≤ k < 2^50 and 2^k ∈ Gamma (mod r).

This is an exact computation with no assumption beyond the arithmetic, and it does not depend on the coins. The program is `r8/bsgs_c_ord.c` (GMP), and the driver `r8/c_ord.py` writes `evidence/c-ord.json`.

- **Constants.** r is the secp256k1 order, and λ is the value of section 1.1. The driver checks that λ³ ≡ 1 (mod r), λ ≠ 1, and g⁶ ≡ 1 for every g in Gamma, and that the λ in `accounting_src.py` is the same value.
- **Method (baby-step giant-step, m = 2^25).** Every k in [0, 2^50) is uniquely k = i·m + j with 0 ≤ i, j < m. The baby table holds 2^j mod r for all j < m, keyed by a 64-bit fingerprint of the full value. For each g in Gamma and each 0 ≤ i < m the giant step is h_i = g·(2^−m)^i mod r, and a match h_i = 2^j gives 2^(i·m+j) = g.
- **Exactness.** A true match has an equal fingerprint, and a probe with no equal fingerprint proves there is no match. Every fingerprint match is verified on the full value, and k is checked again by 2^k mod r = g. The case k = 0 is skipped. Equal baby values with different j would mean ord_r(2) < m, and they are counted.
- **Result.** No match with 1 ≤ k < 2^50 for any of the six elements of Gamma. The baby table has 33,554,432 entries with no equal fingerprints and no equal values, so ord_r(2) ≥ 2^25. The only fingerprint hit is the skipped k = 0. Planted targets (2^777777, 2^987654321 and 2^(2^49)) are found and verified by the same program as controls; the m = 2^25 run takes 60 s on one core.
- **Scope.** The certificate covers escape depths below 2^50. Depths of 2^50 or more are excluded only by the depth tail of section 2.3.

### 2.5 Source of the table, and the restart pool

The algorithm draws the table and every trail's coins from the operating system's randomness (section 2.2). The measurement and the replica use a fixed seed instead (`measurement/walk_ref.py gen <path> <seed>`, which seeds `random.Random(seed)`; `reproduce.sh` uses seed 1 by default). The declared run uses seed 2801 (section 7), and the toy experiment its own seeded RNG. Neither is part of the algorithm, and the seed of the submitted measurement inputs is not recorded in the package.

**Restart pool (disclosure).** The measurement harness (`step.c`, `new_trail`) and its replica take restart starts from a precomputed pool of SPARES = 1024 points, so restarts in the harness reuse starts. The specified algorithm draws fresh coins for every restart (section 1.3, step 3), and the ledger charges the specified restart cost. The kernel was not changed.

### 2.6 Formal-structure test (evidence/formal-structure.json)

`r8/formal_test.py` replays small builds (2^6 adders, 256 walks, 3000 batch steps, theta = 2^−6, seed 1) with real secp256k1 points and checks each walk point's formal vector and pair (a, b) against the walk. Of 767,593 fresh-source steps, 1,992 are entries (1,980 two-cycle, 12 longer); the maximum number of persistent steps before one escape is 62, within the bound 128 of Lemma E2. The 1,980 formal 2-cycles have rate 2.621·10^−3 per candidate step against 1/(6A) = 2.604·10^−3 (z = 0.29). All 1,992 escapes have 0 formal equalities with earlier points of their trail (D1), and the maximum escape depth is 3. The 12 longer entries are fewer than the 15.6 predicted, and no formal 4-closure occurred; this is not explained, and the abab explanation of an earlier round is withdrawn (it was computed on windows that were not fresh). A lower count does not raise any rate used in the bound. These small builds test the structural claims on the real group parameters, not M1 at production size.

### 2.7 Degenerate pairs

Formally distinct expressions can have equal pairs (a, b) up to Gamma, and the pair is a linear map of the table coefficients. For two formally distinct expressions e, e' of one trail, the difference of their pairs is a nonzero linear form in the coefficients v = (a_S, b_S, c_1, d_1, ...). By Schwartz–Zippel, Pr[(a,b)_e = g(a,b)_e' for a fixed g in Gamma] ≤ 1/r over the algorithm's coins. After conditioning on the earlier history, which removes at most a fraction |Gamma| t/r of the coefficient space, the conditional bound is at most 2|Gamma|/r per pair (section 3).

## 3. Coincidences: Lemma G (nondegenerate and degenerate), DELTA_LOOP and Lemma L

**Lemma G.** Let the next point of a trail have t earlier distinct classes, and let there be no earlier coincidence. Then

  Pr[the new point's class is among the t earlier classes | history] ≤ |Gamma| t / (r − |Gamma| t(t−1)/2) ≤ |Gamma| t / (r − |Gamma| TRAIL_MAX²/2).

*Proof.* On the event that there is no coincidence so far, the transcript is a function of the class data of distinct classes, which are independent of x by M1. Conditioned on it, x is uniform on the complement of the roots of the earlier pairs' linear equations: at most |Gamma|·C(t,2) roots, one per Gamma element per earlier pair. The new point's coincidence roots are at most |Gamma| t. Since the coins are fixed and x is uniform, this holds for every coin outcome without degenerate coincidences, and so over the coins and x together. ∎

**Within-trail coincidences.** A within-trail coincidence is a pair of points of one trail that is either a nondegenerate coincidence (Lemma G) or a degenerate one (section 2.7). Summing over the points of every trail and adding the table-level term A² |Gamma| / r for the degenerate pairs of the table,

  DELTA_LOOP := K_TRAILS · (TRAIL_MAX²/2) · |Gamma| · (2/r + 1/(r − |Gamma|·TRAIL_MAX²/2)) + A²·|Gamma|/r ≈ 1.0503·10^−16

bounds the probability of any within-trail coincidence, over x and the coins.

**Cut-off and false filter matches.** A divergence between a live trail and its replay needs a cycle of length at least 65 that one anchor schedule detects and the other misses. A false 64-bit filter match alone cannot cause it, since an escape needs an exact check. Such a cycle is a within-trail repeat, covered by DELTA_CUT = DELTA_LOOP + DELTA_LONG ≈ 3.2929·10^−12, where DELTA_LONG = S_max·Q_LONG is the union bound over the processed starts of a closure of length 65 or more (Lemma E1, with q = 1/6, a premise of section 2.2). With the sensitivity value q = 0.18 the same bound is DELTA_LONG ≈ 4.9·10^−10. No per-trail term is charged for false matches.

**Lemma L (no re-entry into an escaped cycle).** By D2, every re-entry into an escaped cycle is a within-trail coincidence, so its probability over the run is at most DELTA_LOOP.

## 4. Closure entries: Lemma E1 and Lemma E2

**Definitions.** A step is an *entry* if its source point lies on a fresh class (not visited earlier in the trail) and its successor lies on a class already visited in the trail. A step whose source lies on a visited class is *persistent*. Only entries open a closed cycle, and every persistent step follows an entry.

**Lemma E1.** Pr[step is an entry | history] ≤ p_c := 1/(|Gamma| A) + c/A² + Q_LONG + 2^−63, with c = 3, where Q_LONG bounds the closures of length 65 or more. A block of l ≥ 65 fresh steps closes with probability at most |Gamma| q^l exp(l²/(2Aq)) for l ≤ A/2, with q = 1/6 under L-ori, and the lengths above A/2 add at most TRAIL_MAX·|Gamma|·2^(−0.43A). Summed over l ≥ 65 this gives log2 −164.894 (dominated by l = 65), so Q_LONG ≤ 2^−164.89; with q = 0.18 it is log2 −157.674. accounting.py evaluates both sums and checks the terms.

*Proof.* The window of a closure is a run of ℓ distinct classes C_1..C_ℓ whose successor is C_1. Each C_p is fresh when visited (otherwise the window contains an earlier coincidence, which the good event excludes), so by M1 its index i_p = (K_C >> 1) mod A is uniform over A and independent across the window. By L-ori its orientation, and so the unit u_p of its step, is uniform on Gamma. A formal closure needs, for every index a in the window, the units of its occurrences to sum to zero. An index occurring once cannot cancel, so every index has multiplicity at least 2 and the number D of distinct indices is at most ℓ/2. For an index of multiplicity m ≥ 2 with independent uniform units, the last unit must equal one specific value, so the zero-sum has probability at most 1/6. Independence over the distinct indices gives a factor at most 6^(−D), and the closing class's identity γ in Gamma is a conservative factor 6. The rate of windows of length ℓ is therefore at most

  rate_ℓ = 6 · Σ_{D=1}^{⌊ℓ/2⌋} S2(ℓ, D) · A(A−1)…(A−D+1) / A^ℓ · 6^(−D),

where S2(ℓ, D) counts the partitions of ℓ positions into D blocks of size at least 2. For A = 2^16, rate_3 = 1/A² exactly, rate_4 = 1/(2A²) + 1/(2A³), and Σ_{ℓ=5}^{64} rate_ℓ = 7.40·10^−15 = 2^−46.9. So the windows of length 3 to 64 together have rate at most (1.5 + 4·10^−5)/A² ≤ 1.51/A²; the 2-cycle is the separate exact rate 1/(6A) at point level. The accounting keeps c = 3, which is above this sum. The 2^−63 term covers false 63-bit filter matches, which the exact check rejects. ∎

**Lemma FR (forced return; exact).** Let steps p and p+1 be fresh, with tau_p = −1, where tau_p is the canonical rotation of P'_p = canon(C_p) + R_{i_p}, and the transition law is u_{p+1} = u_p tau_p^{−1}. If i_{p+1} = i_p, then P'_{p+1} = −C_p, so tau_{p+1} = −1 and C_{p+2} = C_p.

*Proof.* tau_p = −1 gives C_{p+1} = canon(P'_p) = −P'_p. Then P'_{p+1} = C_{p+1} + R_{i_{p+1}} = −P'_p + R_{i_p} = −C_p. C_p is canonical, so canon(−C_p) = C_p, which gives tau_{p+1} = −1 and C_{p+2} = C_p. ∎

Let P = Pr[i_{p+1} = i_p | tau_p = −1, both steps fresh], and let rho be the non-equal share, Pr[tau_{p+1} = −1 | tau_p = −1, both fresh, i_{p+1} ≠ i_p]. Lemma FR gives Pr[tau_{p+1} = −1 | tau_p = −1, both fresh] = P + (1 − P)·rho, so the excess over 1/6 is

  E = (5/6)·P + (1 − P)·(rho − 1/6).

With rho = 1/6 (L-ori), E = (5/6)·P exactly. Under index uniformity P = 1/A, so E = 5/(6A). The forced return has rate P/6 = 1/(|Gamma| A), which is the 2-cycle rate of Lemma E1.

*Measurements* (evidence/sign-law.json, from evidence/sign_law.py; small builds, fresh steps with tau_p = −1):
- Equal-index rate times A: 0.967 (A = 8), 0.992 (A = 64), 1.015 (A = 256), 0.985 (A = 1024), each with a 95% interval containing 1. Every one of the 10,219 equal-index events had tau_{p+1} = −1, as Lemma FR requires.
- Share-based excess times A: 0.788 (A = 8), 0.808 (A = 64), 0.856 (A = 256), 0.828 (A = 1024), against 5/6 = 0.833. The A = 1024 interval is [0.27, 1.38]; the equal-index rate is the sharper test.
- Non-equal share rho: 0.16417 (A = 8, z = −0.79), 0.16638 (A = 64, z = −0.39), 0.16671 (A = 256, z = 0.08), 0.16667 (A = 1024, z = 0.03). At A = 64 (seed 1, 128,369 pairs) the 1,999 equal-index pairs have share 1, and the non-equal pairs split by adder sign or rotation all have |z| ≤ 0.54 against 1/6.
- The A = 1024 seeds 1 and 2 each have 874 equal-index pairs. Both counts were recounted from the walk-level histograms, and the two runs are independent (chance about 1%). It is recorded, not treated as a signal.

An abab suppression reported in an earlier round, and the λ over-representation reported with it, are withdrawn: both were computed on windows that were not fresh. No measurement in this package contradicts L-ori, which remains a premise.

**Persistent steps.** The window of an entry runs from the entry to the escape that resolves it. Inside it the walk stays on the visited cycle, so no fresh source occurs, and windows of different entries are disjoint. Each persistent step is an ordinary processed step at the common cost; the checks that resolve an entry are charged to the entry (W_ENTRY includes its look-ahead).

**Lemma E2 (persistent steps).** On the good event of section 3 (which includes every cycle of length 65 or more through the cut-off of DELTA_CUT), each entry is followed by at most 128 persistent steps before the escape that resolves it, deterministically.

*Proof.* A 2-cycle entry is detected at the entry step itself, so it has no persistent steps. A longer cycle of length l ≤ 64 is entered at the entry step. Within 64 steps the walk reaches an anchor step, the anchor point lies on the cycle, and the walk returns to the anchor's class within l further steps, where the anchor check fires and the escape follows. So there are at most 64 + l ≤ 128 persistent steps, and processed steps ≤ (fresh-source steps) + 128 · (entries). ∎

## 5. Abandonment: Lemma T'' (first moment over entries)

Let n = 64 (CAP_PERIODS + 1) = TRAIL_MAX be the processed steps before a trail is abandoned, with CAP_PERIODS = 2^DP_BITS, so n·theta is about 64. Let E_n be the number of entries among the first n steps of a trail. On the trail's coincidence-free event G_i, the fresh-source steps number at least n − 128 E_n (Lemma E2). Each fresh-source step is a distinguished point with probability theta, independently of the past, because its key K_C is uniform (M1; the DP test reads only K_C). A persistent step is never a distinguished point, since its class was visited earlier and a distinguished class would have ended the trail. So with no distinguished point among the fresh-source steps, the trail is abandoned, with probability at most (1 − theta)^(n − 128 E_n).

*First moment over entries (chain rule).* For k ≥ 1, Pr[E_n ≥ k] ≤ binom(n, k)·p_c^k ≤ (e·n·p_c/k)^k, because each of k given steps is an entry with conditional probability at most p_c given the history. Take k = floor(n/c) with c = 2^16 (A = 2^16). Then e·c·p_c = 0.45 ≤ 1/2 and 128/c ≤ 2^−7, so

  Pr[abandon_i and G_i] ≤ 2^(−n/c) + exp(−64·(1 − 2^−7)) = exp(−63.5) + 2^(−2^50) ≤ 2.6442·10^−28.

*Failure form.* Let X be the number of abandoned trails and G the good event (no within-trail coincidence and no cut-off). The run fails by abandonment only if G fails or X > K, so Pr[abandonment failure] ≤ Pr[not G] + E[X·1_G] / K. With E[X·1_G] ≤ K_TRAILS · 2.6442·10^−28, DELTA_ABANDON = K_TRAILS · 2.6442·10^−28 / K = 1.528·10^−17, and Pr[not G] ≤ DELTA_CUT (section 3).

## 6. Expected rare work: Lemma M'' (each entry charged by its type)

Per processed step, the expected metered rare work is at most E_RARE, with each entry charged by its type:

| event | rate | worst case (instructions) |
| --- | ---: | ---: |
| 2-cycle entry (true or false predecessor match) | 1/(|Gamma| A) per fresh source (Lemma E1) | 48,039 |
| longer entry, 3..64 (anchor check, with its look-ahead steps) | c/A² per fresh source, c = 3 | 1,618,122 |
| false 63-bit filter match | 2^−63 per processed step | 48,039 |
| exact-k fallback, per canon call | 2^−55 (M1 extended to the canonical coordinate fractions, and L-ori for the rotation) | 256 |
| within-trail coincidence, or closure of length 65 or more (sections 3, 4), per step | |Gamma|·TRAIL_MAX·(2/r + 1/(r − |Gamma|·TRAIL_MAX²/2)) + Q_LONG ≈ 2.3·10^−50 | TRAIL_MAX·(821 + reserve) |

Persistent steps have no entry charge of their own: they are ordinary processed steps at the common cost C = 821, covered by S_max.

The fallback row uses M1 (fractions) and L-ori (rotation): the test of Lemma 1 fails with probability at most 204/2^64 < 2^−55.6 for uniform fractions. The charge is the whole canon call with the fallback taken, 256 (static bound 253), per canon call.

The within-trail coincidence row is self-referential, because its worst case contains the reserve. The fixed point is E_RARE = (E_FIXED + rate·TRAIL_MAX·C_COMMON) / (1 − 64·rate·TRAIL_MAX), with E_FIXED the first four rows: E_FIXED = 0.12329973978940546 (the fallback row adds 7.1·10^−15), E_RARE = 0.12329973978940546 (the coincidence row adds 1.4·10^−27). The reserve per processed step is 64·E_RARE = 7.891183, and the charge per walk step is C_COMMON + reserve = 828.891183.

Markov's inequality gives Pr[metered rare work > 64·S_max·E_RARE] ≤ 1/64, counted in delta.

## 7. Success probability

Let T_steps = K·L. Then tau = T_steps²/(2N) = 0.530300. Under H1, the probability that the first T_steps steps contain no useful class collision is at most exp(−tau·(1 − eps)), with eps = 4/A + 2·(3/A) + 2·(3·192/A²) = 0.00015286: the adder penalty (4/A), 2-cycles at a rate of at most 1/A with 3 non-fresh states each, and short cycles at a rate of at most 3/A² (section 4) with 192 wasted states each.

A class collision among the first T_steps steps between two trails leads, since merged trails stay merged, to the same first distinguished class. The record that comes second triggers resolution; Lemma R makes the replays reach the recorded point, and Lemma D makes the equation degenerate with probability exactly 1/r. Fresh start points only add points, and restarts at distinguished points give each record its own coins.

The failure terms are:
- **Within-trail coincidences and cut-off replays (sections 1.3 and 3)**, including re-entry into an escaped cycle (Lemma L), degenerate pairs and cycles of length at least 65: at most DELTA_CUT = DELTA_LOOP + DELTA_LONG = 3.2929·10^−12 (this is Pr[not G]).
- **Abandoned trails (section 5)**: at most DELTA_ABANDON = 1.528·10^−17.
- **Escape depth tail (section 2.3)**: at most DELTA_DEPTH = 2^(−1.29·10^15), below double precision.
- **Degenerate first match:** exactly 1/r (Lemma D).
- **Work meter abort (section 6):** at most 1/64.
- **DP_MAX overflow:** by Chernoff under H1, at most exp(−DP_MAX/12), below 2^−60.
- **Exceptional additions:** at most 2^−110 under H1.
- **Other terms:** together at most 2^−20 = delta.

Therefore:

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015286)) − DELTA_CUT − DELTA_ABANDON − DELTA_DEPTH − 1/64 − 2^−20 = 0.395898 ≥ 0.39.**

The probability is over the uniform target x and the algorithm's coins. The model term is 0.411524. The budget is 1.0299·√N steps, against 0.9943·√N for success exactly 0.39. The 0.0215 of slack in the model term covers the sampling error of the toy evidence below.

**Evidence for H1 and the procedure.** `experiment.py` runs the procedure of section 1.3 end to end on random prime-order toy curves y^2 = x^3 + 7 over p ≡ 1 mod 3, with the same class key, canonical rule, adder images, index and both cycle rules, trails from fresh coins with the trail cap, and resolution by replaying both trails from their coins. It never tracks coefficients while walking and has no oracle collision table. The live walk anchors on a global round counter, as the batched kernel does; the replay anchors on the trail's own step count, so every replay also checks Lemma R. A trial counts as a success only if the first key match resolves to a verified x, and each trail start is charged as one of the slot's L steps, which is conservative.

<!-- toy-evidence -->
Offline runs (`evidence/trails/`, pooled in `summary-double.json`): 5,200 trials on 5 curves with N = 2^26.8, 2^29.1, 2^30.8, 2^31.3, 2^32.5 classes, adder tables of 2^11 and 2^16, 16 slots and theta = 2^−5, escape "double".
- Success of the whole procedure at tau = 0.5303: 2,169 of 5,200, 0.4171, Wilson 95% interval [0.4038, 0.4306]. The model 1 − exp(−tau) is 0.412. Every success was verified ([x]P = Q). Degenerate first matches: 0. Abandoned trails: 0. Longer fruitless cycles escaped: 6. The finishing phase added at most 2.6% to the K·L steps.
- Production 64-bit filters: every recorded trail was replayed with a different anchor schedule, 4,619,736 replays, 0 that missed the recorded point.
- Lemma R stress (key `trails-stress-double` in `evidence/stress-and-sensitivity-double.json`): 8-bit filters over 5,200 trials, so false filter matches were frequent. 4,837,229 trails were replayed with a different anchor schedule, and 0 missed the recorded point; 349,511 distinguished points were first reached during a look-ahead. These runs' success rates are not used, since each such look-ahead adds up to 65 points that the step budget does not count.
- Replay-order sensitivity (key `lemma-r-sensitivity-double`): 200 trials per curve with 6-bit filters and 16 adders. With the current order, 0 of 59,977 replays (j0) and 0 of 95,625 (a3) missed the recorded point. With the earlier order (an offline switch), 8 of 59,351 (j0) and 47 of 92,993 (a3) did.
- 2-cycle entry rate per step: 7.25e-05, 7.74e-05, 7.73e-05, 7.76e-05 with A = 2^11 (1/(6A) = 8.14e-05); 2.03e-06 with A = 2^16 (1/(6A) = 2.54e-06).
<!-- /toy-evidence -->

The declared initial input (`bits` 28, `b` 7, 80 trials, seed 2801, production 64-bit filters, every trail replayed, escape "double"; `evidence/declared-run-double.json`) runs in a few seconds. It reports 35/80 = 0.4375, Wilson 95% interval [0.334, 0.547], against the model 0.412, with 0 degenerate matches, 0 distinguished points reached during a look-ahead and 0 replays, of 15,059, that missed their recorded point. With 80 trials the standard error is about 0.05, so it is a reproducibility check; the success claim rests on the pooled offline trials above.

Limits: the toy curves have N up to about 2^33, against 2^253 in production. There theta is 2^−5, so start points and the finishing phase are a few percent of the steps; in production they are below 2^−30 of them. H1 is an extrapolation in N, supported by the size-free random-mapping argument and the cited literature. It is not proved for this walk.

## 8. Computation

**No reference operation.** Every instruction of a walk step is charged at 1/C with C = 1000, as COST_MODEL.md's H + R/C with H = 0. No reference operation is charged.

**Charge per walk step.** Two parts, charged separately:
- The common batched path is charged C_COMMON = 821 instructions per walk step: the larger of the plain and anchor steps, each measured as the whole-run instruction total difference between runs with one more batch step on production inputs, per walk. On the measured inputs the plain step is 808.67 and the anchor step 818.66. The one data-dependent branch of the common path, the canon_pick selection (k_max = 2), adds at most 2 instructions per walk step, which gives 820.66; C_COMMON = 821 bounds that.
- The rare events, each charged its worst case (section 6, and section 8.1 for the primitives): a 2-cycle check at most 48,039 instructions, an anchor check with its look-ahead steps at most 1,618,122, and a canon fallback at most 256. Their expected counts are covered by the Markov reserve 64·E_RARE = 7.891183 per walk step (Lemma M'').

The charge per walk step is therefore 828.891183 instructions.

**Cost is uniform over the coin values.** The common batched path is branch-free apart from canon_pick and the exact-k fallback, which is metered at its forced maximum (section 6). No other step reads the coefficient values through a branch, so the instruction count of each walk step is bounded uniformly over all coin outcomes, and the ledger charges that bound rather than an average. The measured totals on the two input sets (seeds 1 and 2) agree, as this predicts.

**Field arithmetic.** `measurement/secp_field.h` implements multiplication and squaring modulo p = 2^256 − c, c = 2^32 + 977, as GCC inline AArch64 assembly, with the final carry selecting between y and y − c. `measurement/secp_sub.h` holds the subtractions. The routines matched Python big-integer arithmetic on about 194,000 reduced and 199,000 arbitrary 256-bit operand pairs, including constructed cases for the rare second-fold carry. The replica check below compares every final point of 2^16 walks after 2^19 steps with Python. Callgrind counts this inline assembly exactly as it counts compiled C.

**8.1 Worst cases of the rare-path operations.** Each rare event is charged a forced-path maximum, not a sampled average. For each primitive (f_single, add_single, escape, same_class, point_is_dp, canon_k with its exact fallback) the data-dependent branches were enumerated from the source and the disassembly. Conditional corrections in the field routines were replaced, for the measurement only, by branch-free routines that cost at least as much on every path (a bounding build, checked for equality on random and edge inputs), and the control variants were forced by chosen inputs, with the maximum taken. Each count is linear in the number of calls and identical for a second seed. Measured maxima (Callgrind Ir per call, organizer image and flags, `measurement/forced/summary.json`):

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

A 2-cycle check costs at most f_single + same_class + escape(2) + point_is_dp plus a control overhead of at most 1,024 instructions, and an anchor check costs at most f_single + escape(65) plus, for each of its at most 64 look-ahead steps, f_single + point_is_dp + same_class, with the same overhead per step. The forced-path harness compiles the primitives in its own translation unit, where GCC compiles finv about 36 instructions more expensive, so the harness figures are conservative for every primitive that calls finv. As an independent check, the longest control-flow path of the unmodified binary, with every conditional branch free and finv at its measured cost, gives at most 47,193 instructions for a 2-cycle check and 1,547,316 for an anchor check, both below the charges above (`measurement/forced/summary.json`). Scalar multiplications and replays use f_single's worst case as their unit, an upper bound for an addition or a doubling.

**Measurement (H3).** `measurement/step.c` runs the trails, steps, cycle handling, distinguished points, restarts and trail cap. Its common batched path was measured on production inputs as whole-run Callgrind `Ir` totals (GCC 14.2.0, organizer's image, `-O3 -fno-tree-vectorize -mgeneral-regs-only`), with the rare paths out of line (`evidence/measurement/production-totals.json`). The first batch step has no anchor and no rare event. The anchor step is 10.0 instructions per walk step above the plain step, and its total contains its own rare events.

- Plain step (steps 1 and 2), per walk step: 808.67 instructions (measured).
- Anchor step (steps 63 and 64), per walk step: 818.66 instructions (measured).
- C_COMMON = 821 per walk step: the measured anchor step plus 2 for the canon_pick branch, rounded up. The bound for walks that are all long is 819.34, so 821 is conservative.

**Rare-path checks.** Distinguished points and false 64-bit filter matches do not occur in a production-size sample, so small builds of the same source and replica exercise them (`evidence/measurement/trail-paths.json`, with 2^6 adders and 2^12 walks): DP_BITS = 4 (restarts at distinguished points, with 2-cycles); DP_BITS = 12 with a one-period cap (abandonments); DP_BITS = 6 with the filters cut to 8 bits (2-cycles, longer cycles, and distinguished points first reached during a look-ahead). In every build the replica matched every final point, every event count, and a digest of every distinguished record, look-ahead records included. The digest also detects the round-5 omission, in which a look-ahead distinguished point was not recorded: the mutant build agrees on every count and differs in the digest (`evidence/measurement/record-mutation-double.json`).

**Total.** S_max = K·(L + TRAIL_MAX) ≤ 2^126.749959 walk steps.

| Component | Bound (units) | log2 |
| --- | ---: | ---: |
| walk steps | 828.891183·S_max / C | 126.479214 |
| adder table | 2^16 adders × 2 scalar mults × 768 unbatched ops × 2^15 / C | 31.13 |
| initial walk starts | K × 2 scalar mults × 768 × 2^15 / C | 47.13 |
| trail restarts | (DP_MAX + K) × 2 scalar mults × 768 × 2^15 / C | 82.88 |
| distinguished points | DP_MAX = 2^67.75 entries × 2^20 instructions / C | 77.78 |
| collision replay | 2 trails × TRAIL_MAX steps × 2^13 unbatched ops × (2^15 + 2^10) / C | 84.61 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 15.68 |
| input validation | 2^24 instructions | 14.03 |

Rounded up, **time_log2 = 126.479214**. Preprocessing (adders and initial walk starts) is included. All processors' work is summed, and elapsed time plays no part. The replay row is charged in every run, with 2^13 operations per replayed step (section 1.3, step 7).

## 9. Memory and advice

Peak memory is at most the sum of: 2·DP_MAX·128 bytes for the distinguished-point table (key, the trail's coins, hashing overhead); K·576 bytes of walk state; 2^16·6·64 bytes of adder images; and 16 MiB. That is 2^75.74996 bytes. Nonuniform advice is 0 bytes.

## 10. Limitations

- H1 is the standard but unproved random-mapping heuristic, supported by the literature and toy-curve experiments only. The model term is 0.4115, 0.0215 above 0.39.
- M1 (equivariant random orientation) is a premise, not a theorem (section 2.2), and it carries the claim. Its orientation law L-ori is measured at A = 8 to 1024 (section 4) and not derived at production size; index uniformity is not measured at A = 2^16.
- **The argument is not established at production size.** The toy runs use theta = 2^−5, and the formal test uses 2^6 adders and theta = 2^−6. Neither tests M1 at production theta or adder count.
- The exact-k fallback rate (2^−55 per canon call) follows from M1 extended to the coordinate fractions, and L-ori. It is not measured: no fallback occurred in the measured steps.
- The measured instances use seeded generators (seed 1 for `reproduce.sh`, seed 2801 for the declared run, a seeded RNG for the toy experiment). The algorithm draws its coins from the operating system.
- Escape depths of 2^50 or more are covered only by the depth tail (section 2.3), a first-moment bound over entries under M1 and L-ori.
- Instruction counts come from one implementation and compiler (section 8). The rare operations are charged forced-path maxima (section 8.1), checked for dominance on the forced paths, not proved for every input. The field routines were tested, not formally verified.
- The replay cut-off (65n steps) is in `experiment.py`; the declared run gives 35/80 and 15,059 replays, 0 mismatches.
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

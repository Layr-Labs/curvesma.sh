# secp256k1-ecdlp: parallel Pollard rho on order-6 classes with a sigma walk, all arithmetic in instructions

Claim: for uniform x in [0,r) and Q=[x]P on secp256k1, the algorithm below returns x
with probability at least 0.39. It uses at most **2^126.479214** target-group-operation equivalents in total (score 63239607/500000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^75.749960 bytes. Preprocessing is 2^47.131242 units, and no nonuniform advice is used. This is a classical cold attack on
one target.

The success bound rests on the disclosed heuristic H1 (the random-mapping model for r-adding walks, used for the collision
statistics), on the premise M1 (class data), and on the randomness of the algorithm's coins (section 2.2), and on the
exact certificate of section 2.4 (no k < 2^50 with 2^k in Gamma) and on the depth tail of section 2.3. The cost rests on H3, average instruction
counts measured with the organizer's toolchain, and every run is charged by a work meter with a Markov reserve. This is a
carefully accounted generic attack. It is not a practical break.

**Cost route.** The walk does not use the 6M reference operation. Every instruction of
every step, including all field multiplications, squarings and subtractions, is its own
4×64-bit GCC inline AArch64 assembly and is charged at 1/C with C = 1000 (H = 0 in
COST_MODEL's H + R/C). COST_MODEL.md: "Faster arithmetic and lookup algorithms are welcome;
justify their actual work". A walk step is charged 828.89 instructions, 0.829 units, rather than 1 + R/C.

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
- `evidence/measurement/` holds the production totals of the doubling kernel
  (`callgrind-doubling.json`, `replay16-production.json`) and the small-build checks of the
  restart and abandonment paths (`trail-paths.json`).
- `evidence/formal-structure.json` holds the formal-structure test of section 2.6, and
  `evidence/c-ord.json` the exact certificate of section 2.4. Their programs (`formal_test.py`, `bsgs_c_ord.c`,
  `c_ord.py`) are kept with the work outside the package.
- `experiment.py` is the declared toy-curve experiment, running the whole procedure;
  `evidence/trails/` holds the offline runs of the same file with escape "double" (`summary-double.json`).

## 1. Algorithm

### 1.1 Notation and parameters

p, r and P come from the protected profile. The curve is y^2 = x^3 + 7, the cofactor is 1,
and p ≡ r ≡ 1 (mod 3).

The endomorphism uses beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
and lambda = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72, with
[lambda](x,y) = (beta·x, y). The automorphism group G6 = {±1, ±lambda, ±lambda^2} has
(r−1)/6 nonzero classes, and N = (r−1)/6 + 1. Field values are compared as reduced
integers in [0,p). The units of G6 are the elements of Z_r^* of order dividing 6, and |Gamma| = 6.

Parameters:
- A = 2^16 adders.
- K = 2^32 walks: 2^16 processors × 2^16-walk batches.
- theta = 2^−60.
- L = 33310394024866356565296666183 ≈ 2^94.7500.
- CAP = 64/theta = 2^66. Trails have fewer than TRAIL_MAX = CAP + 64 = 2^66 + 64 steps.
- S_max = K·(L + TRAIL_MAX) ≈ 2^126.749959, the bound on walk steps in a run.
- DP_MAX = 2·theta·S_max + 2K ≈ 2^67.7500, the number of stored distinguished points.
- K_TRAILS = DP_MAX + 2K ≈ 2^67.7500: the number of trails started (K initial trails, DP_MAX restarts after
  records, and K restarts after abandonments).
- K0 = 2^50: the escape-depth bound of section 2.3 (the depth event: every trail has fewer than K0 escapes).

### 1.2 Classes, keys and the sigma walk

For W = (x, y) ≠ O:
- **Class key.** key(W) = y if y is even, else p − y, so key(W) is ±y and even.
  - key is constant on the class, since G6 changes y only by sign.
  - key determines the class: y^2 = x^3 + 7 fixes x up to a factor beta^j.
- **Canonical point.** canon(W) = (beta^k·x, (−1)^s·y), where k minimises beta^k·x and
  s = y mod 2. As a scalar, canon(W) = eps_W W with eps_W = (−1)^s·lambda^k.
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
exactly (the exact-k fallback). Under M1 (section 2.2) the fallback probability per canon call is at most 2^−55
(section 6). `measurement/frac_check.py` checks Lemma 1 and the decision rule on 200,000 random
and edge-case inputs per seed; the largest floor gap seen for g_1 is 6.

The batch path reaches this rule through `select_adder` for every new point of every walk in every batch step, so the
exact-k fallback is on the common path, not only inside the rare paths (section 6).

### 1.3 Procedure

1. **Validate Q.** Check it is on the curve. If Q = O, output 0.
2. **Adders.** For each adder i < 2^16:
   - draw c_i, d_i uniform in [0,r) from the operating system's randomness, fresh in every run and independent of x;
   - set R_i = [c_i]P + [d_i]Q;
   - tabulate its six images.
3. **Trails.** Each of the K walk slots runs a sequence of trails.
   - A trail starts from fresh coins a, b uniform in [0,r), drawn from the operating system's randomness at every trail, including every restart, with W = [a]P + [b]Q (two scalar multiplications, charged in section 8). The
     slot keeps the coins until the trail ends. Coefficients are not tracked while
     walking.
   - A trail ends at its first distinguished point (step 6), or is abandoned: when it
     crosses more than 2^60 anchor boundaries (step 5) without one, that is after between
     CAP − 64 and CAP = 64/theta steps, or when an addition gives O. Either way the slot
     starts a new trail from fresh coins. Every trail has fewer than TRAIL_MAX steps.
   - Each slot takes L steps (the main phase), then finishes its current trail and stops.
     The run fails if more than K trails are abandoned.
4. **Step.** W ← f(W), as one batched incomplete affine addition: Montgomery's
   simultaneous inversion over the 2^16 walks of a batch (3M), then lambda, lambda^2 and
   y3 (2M + 1S), with one Fermat inversion per batch.
   - If the batch product is zero, the affected walk uses the unbatched complete formula.
   - The canonical rotation of each new point is decided by the rule of section 1.2, with its exact fallback.
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
   - **Replay bound (deterministic).** Each live step produces at most one walk step and at most one escape, and an
     escape traverses at most MAX_CYCLE = 64 further cycle steps. So a replay of a trail with n steps takes at most
     65n walk steps on every run. Each replay step uses at most 70 unbatched operations (one evaluation, two for a 2-cycle
     check, 66 for an anchor check, one escape), so a replay takes at most 65·70·n = 4550n ≤ 2^13·n unbatched operations.
     The cost of the two replays is charged in every run (section 8). The replay is cut off after 65n steps (REPLAY_FACTOR = 65 in the core), and a cut-off replay fails the run, counted in DELTA_CUT (section 3).
8. **Work meter.** Every run is charged at most W_max = S_max·828.891183 + G
   instructions, with G the fixed global charges (initial starts, trail restarts, distinguished-point
   records, solve and verify, input validation, adder table, resolution replays).
   - Every walk step is charged C_COMMON = 821 instructions, which bounds the common path
     including its batch inversion, and each rare event is charged its worst case (section 6).
   - The meter aborts the run when the total metered rare work exceeds 64 times its
     expectation bound (Lemma M'). That happens with probability at most 1/64, counted in delta.
   - The run also fails if more than DP_MAX records arrive or more than K trails are abandoned.

**Escape.** A cycle's escape goes to 2·canon(M), where M is the cycle's least-key member. canon(M) is a point of M's class, so 2·canon(M) is a point of the class of 2M. The escape uses no table and no stored point: it is one point doubling, whose cost is charged at its forced-path maximum (Section 8.1). Its class is not one the trail has produced before on the good event of section 3 (Lemma G), and its class data are uniform by M1.

**Escape repeats.** The escape point of a cycle depends only on its least-key member, so it is the same every time the trail escapes from that cycle. A repeated escape point therefore means that the trail has re-entered a cycle it already escaped from. Section 3 bounds this event (Lemma L).

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

The replay length is bounded as in step 7 of section 1.3, deterministically on every run.

**Lemma D (degenerate matches).** The first key match has b1 = b2 with probability exactly
1/r.

*Proof.* A trail records once, so the two records come from different trails t1 ≠ t2,
possibly of the same slot. Condition on all start points and adders. That fixes every
trail's points, and so which records match first.
- Given its start point W = [a + b·x]P, a trail's coin b is uniform on [0, r), because
  each b determines exactly one a.
- The coins are drawn from the operating system's randomness (section 2.2), so this uniformity holds exactly over the
  algorithm's coins, and the argument below needs no heuristic.
- Coins of different trails are independent.
- Along a replay, b changes as b ← u·b + d, where u is a unit (order-6 classes {±1, ±λ, ±λ²}, up to sign)
  or twice one, and u, d depend only on the path.

So the canonical b at the distinguished point is U_i·b_ti + D_i with U_i invertible mod r,
and P[b1 = b2] = 1/r. ∎

## 2. Premises, formal structure and checks

### 2.1 Heuristic H1

H1 (heuristic, disclosed). The class map [W] → [canon(W) + R_h([W])], with uniform
random adders, behaves like a random mapping on the N classes, apart from the r-adding
penalty. That penalty is a factor of about (1 − 1/(6A)) in the birthday exponent (Teske
2001; Bernstein–Lange 2012; Bos–Kleinjung–Lenstra 2010). H1 is used only for the collision
statistics of section 7, not for the within-trail statements of sections 2 to 6.

### 2.2 Premises: M1, and the algorithm's randomness

**M1 (class data).** Every class-level value the walk reads is a function of the class alone, and the values of distinct
classes are independent and uniform, independent of x. The values are: the class key (which fixes the adder index, the
64-bit filter bits and the DP test), the canonical rotation k and canonical sign s, and the coordinate fraction values of
the canonical representative used by the exact-k test. M1 is the random-mapping heuristic applied to all class data, not
only to the successor. The fraction values are included deliberately: the fallback rate of section 6 is derived from
them, and no other property of the coordinates is assumed. Probabilities over M1 are over the class data, within the probability space of the target and the algorithm's coins.

**Randomness of the algorithm (not a heuristic).** The adder coefficients (c_i, d_i) and the start coefficients are drawn
uniformly from the operating system's randomness (`secrets.randbelow(r)` or equivalent), fresh in every run, independently
of the target by construction. Every trail, including every restart, draws fresh coins (a, b) uniformly from the same source
and computes W = [a]P0 + [b]Q by two scalar multiplications. The probability space is the uniform target x together with
these coins. Random(1) is only the seed of the reproducible measurement instance (section 2.5). The probabilities in
sections 3 to 7 are over the uniform target x and the algorithm's coins, and over the class data (M1), and nothing else is
assumed. This premise replaces the table premise of earlier drafts.

### 2.3 Formal structure: escape depth, D1 and D2

Within one trail with start S, each walk point has the formal form
  W = 2^m g S + sum_i n_i R_i,
with g in Gamma, integers n_i, and m the number of escapes so far in the trail. Steps add sigma(R_i) (which is
±lambda^j R_i) and do not change the S-coefficient. An escape replaces M by 2 canon(M) = 2 g' M, which doubles
the S-coefficient. Its pair (a, b) is the same linear map applied to the pairs (a_S, b_S) and (c_i, d_i).

**Escape depth (depth tail).** An escape is a closure. Under M1 the closure steps in a trail of n = TRAIL_MAX steps are
dominated by Bin(n, p_c), with p_c = 1/(6A) + 3/A² + 2^−63 (section 4), so by the first moment of section 5,

  Pr[a trail has K0 or more escapes] ≤ (e·n·p_c / K0)^K0,   K0 = 2^50.

Here e·n·p_c/K0 = 0.45317 ≤ 1/2, so the bound is at most 0.45317^(2^50), that is 2^(−1.29·10^15). Over all K_TRAILS trails,
the probability that some trail has K0 or more escapes is at most

  DELTA_DEPTH = K_TRAILS · 0.45317^(2^50) = 2^(−1.2856·10^15),

which is counted in delta (section 7). Below, the **depth event** is the complement: every trail has fewer than K0 escapes.
On the depth event every escape depth m+1 satisfies m+1 ≤ K0 − 1.

**D1 (depth).** On the good event of section 3 and the depth event, if E is the (m+1)-th escape point of a trail, its
S-coefficient 2^(m+1) g' differs from the S-coefficient 2^j g'' of every earlier point of the trail (j ≤ m).
*Proof.* Equality would give 2^k in Gamma for k = m+1−j, with 1 ≤ k ≤ m+1 ≤ K0 − 1 < 2^50. Section 2.4 shows that
no such k exists. ∎

**D2 (re-entry is not structural).** A re-entry into an escaped cycle is a later point (depth ≥ m+1) in the classes of an
earlier cycle (depth ≤ m). By D1 such a point is never formally equal to an earlier one, so every re-entry is a
within-trail coincidence, which section 3 bounds.

### 2.4 Exact certificate for D1: no k below 2^50 with 2^k in Gamma

**Claim.** For the secp256k1 order r and Gamma = {1, −1, λ, −λ, λ², −λ²} ⊂ Z_r, there is no k with 1 ≤ k < 2^50 and
2^k ∈ Gamma (mod r).

This is an exact computation with no assumption beyond the arithmetic. The claim does not depend on the coins. The program is `r8/bsgs_c_ord.c` (GMP), and the
driver `r8/c_ord.py` writes `evidence/c-ord.json`.

- **Constants.** r is the secp256k1 order, and λ is the value of section 1.1. The driver checks that λ³ ≡ 1 (mod r),
  λ ≠ 1, and g⁶ ≡ 1 for every g in Gamma, and that the λ in `accounting_src.py` is the same value.
- **Method (baby-step giant-step, m = 2^25).** Every k in [0, 2^50) is uniquely k = i·m + j with 0 ≤ i, j < m.
  The baby table holds 2^j mod r for all j < m, in open addressing, keyed by a 64-bit fingerprint of the full 256-bit value.
  All entries are kept, including fingerprint collisions. For each g in Gamma and each 0 ≤ i < m the giant step is
  h_i = g·(2^−m)^i mod r, and a match h_i = 2^j gives 2^(i·m+j) = g.
- **Exactness.** The fingerprint is a function of the value, so a true match has an equal fingerprint, and a probe that
  finds no equal fingerprint proves that there is no match. Every fingerprint match is verified on the full value, and
  then k is checked again by 2^k mod r = g. The case k = 0 (g = 1 with i = j = 0) is skipped. Equal baby values with
  different j would mean ord_r(2) < m, and they are counted.
- **Result.** No match with 1 ≤ k < 2^50 for any of the six elements of Gamma. The baby table has 33,554,432 entries with no
  equal fingerprints and no equal values, so ord_r(2) ≥ 2^25. The only fingerprint hit is the skipped k = 0.
- **Controls.** The same program finds planted targets: 2^777777 with m = 2^10, 2^987654321 with m = 2^20, and 2^(2^49)
  with m = 2^25 (the upper half of the range). Each is found and verified, and no element of Gamma matches in those runs.
- **Run time.** The m = 2^25 certificate took 60 s (the baby table 4.5 s), and the whole driver with compilation and the three
  controls took 132 s, on one core of Apple silicon with GMP.
- **Scope.** The certificate covers escape depths below 2^50, which is the depth event. Depths of 2^50 or more are excluded
  only by the depth tail of section 2.3.

### 2.5 Source of the table, and the restart pool

The algorithm draws the table and every trail's coins from the operating system's randomness (section 2.2). The
measurement and the replica use a fixed seed instead. `measurement/walk_ref.py gen <path> <seed>` seeds
`random.Random(seed)`, draws the target x = Q/G first, and then draws, for each of the NA + B + SPARES points, the pair
(c, d) (c first, then d), and writes the point c·G + d·Q. The file keeps the points only, so the coefficient pairs follow
from the seed and the draw order; `formal_test.py` reproduces them and checks them against the points. `reproduce.sh` uses
seed 1 by default. Seed 1 is the seed of the reproducible measurement instance (`reproduce.sh`, the declared run and the
toy experiment with its seeded RNG). It is not part of the algorithm. The seed of the inputs used for the submitted
measurement is not recorded in the package; the status notes mention seeds 1 and 2 for different runs.

**Restart pool (disclosure).** The measurement harness (`step.c`, `new_trail`) and its replica take restart starts from a
precomputed pool of SPARES = 1024 points, so restarts in the harness reuse starts. The specified algorithm draws fresh coins
for every restart (section 1.3, step 3). The restart path is not part of the common-path measurement, and the ledger charges
the specified restart cost, two scalar multiplications per restart. The trail tests with more than 1024 restarts exercise the
pool; they agree with the replica by construction. The kernel was not changed.

### 2.6 Formal-structure test (evidence/formal-structure.json)

`r8/formal_test.py` replays small builds of the walk (2^6 adders, 256 walks, 3000 batch steps, theta = 2^−6, seed 1) with
real secp256k1 points. It tracks, for every walk point, its exact formal vector over the trail's start S and the adders
(W = s·S + Σ n_i R_i), and its pair (a, b) with W = (a + b x)G, which it checks against the point for every walk at the
end. Adder images and canonical points are checked as scalar multiples. The lambda used is the one for which [lambda]
acts as (x, y) → (beta·x, y).

Results of the final run (`evidence/formal-structure.json`; 768,000 candidate steps, of which 755,459 have a predecessor in their trail):
- Formal 2-cycles: 1,980 at point level, and the same 1,980 at class level, so every class-level 2-cycle the walk checks is a formal one (the real class-level count is also 1,980). The rate 2.621·10^−3 is against the predicted 1/(6A) = 2.604·10^−3 (Wilson 95% interval [2.51, 2.74]·10^−3; z = 0.29).
- Formal 3-closures (W_{k+1} = W_{k−2}): 12 entries, meaning closures that were not already closures at the previous step, against 10.2 predicted by 1/(18A²). A closure persists until the anchor check detects it, so it is counted at each of its steps (419 closure steps in total); the entries are the quantity the closure rate refers to.
- Formal 4-closures: 0 entries, against 5.1 predicted by 1/(36A²). This shortfall is about a 0.6% event under the prediction and is not explained. It lowers the observed closure count, so it does not weaken the upper bound of section 4, but the 4-cycle term of section 4 is not confirmed here.
- Escapes: 1,992 (1,980 two-cycle, 12 long-cycle). Each escape point was compared with every earlier point of its trail: 0 formal equalities at point level and 0 up to Gamma, as D1 predicts. The maximum escape depth in a trail was 3 (depth 1: 1,768 escapes; depth 2: 211; depth 3: 13).
- Consistency: the formal pairs (a, b) reproduce the walk point for all 256 walks at the end of the run; no distinguished point was first reached in a look-ahead, and no trail was abandoned.

These are small builds: 2^6 adders and theta = 2^−6, not production sizes. They test the structural claims (D1 and the
closure rates) on the real group parameters, and they do not test M1 at production size.

### 2.7 Degenerate pairs

Formally distinct expressions can have equal pairs (a, b) up to Gamma, and the pair is a linear map of the table
coefficients. For a fixed transcript and two formally distinct expressions e, e' of one trail, the difference of their
pairs is a nonzero linear form in the coefficients v = (a_S, b_S, c_1, d_1, ...) with values in Z_r^2. By Schwartz–Zippel,
Pr[(a,b)_e = g(a,b)_e' for a fixed g in Gamma] ≤ 1/r over the algorithm's coins. After conditioning on the earlier history, which removes at
most a fraction |Gamma| t/r of the coefficient space, the conditional bound is at most 2|Gamma|/r per pair (section 3).

## 3. Coincidences: Lemma G (nondegenerate and degenerate), DELTA_LOOP and Lemma L

**Lemma G.** Let the next point of a trail have t earlier distinct classes, and let there be no earlier coincidence. Then

  Pr[the new point's class is among the t earlier classes | history] ≤ |Gamma| t / (r − |Gamma| t(t−1)/2)
    ≤ |Gamma| t / (r − |Gamma| TRAIL_MAX²/2).

*Proof.* On the event that there is no coincidence so far, the transcript (adder indices, DP tests, canonical data) is a
function of the class data of distinct classes, which are independent of x by M1. Conditioned on it, x is uniform on the
complement of the roots of the earlier pairs' linear equations. There are at most |Gamma|·C(t,2) such roots, one per Gamma
element per earlier pair, so the roots excluded are |Gamma| t(t−1)/2 in number. The new point's coincidence roots are at most
|Gamma| t. Since the algorithm's coins are fixed and x is uniform, this holds for every coin outcome without degenerate
coincidences, and so over the coins and x together. ∎ (Earlier drafts used |Gamma| t roots; that counted only the points, not
the pairs.)

**Within-trail coincidences.** A within-trail coincidence is a pair of points of one trail that is either a nondegenerate
coincidence (Lemma G) or a degenerate one (section 2.7). Summing over the points of every trail (up to TRAIL_MAX points
each, over K_TRAILS trails) and adding the table-level term A² |Gamma| / r for the degenerate pairs of the table,

  DELTA_LOOP := K_TRAILS · (TRAIL_MAX²/2) · |Gamma| · (2/r + 1/(r − |Gamma|·TRAIL_MAX²/2)) + A²·|Gamma|/r ≈ 1.0503·10^−16

bounds the probability of any within-trail coincidence. The probability is over x and the algorithm's coins (section 2.2).
Replacing the denominator r − |Gamma| TRAIL_MAX by r − |Gamma| TRAIL_MAX²/2 (the pairs, not the points) changes the printed
value by less than its last digit, because the two denominators differ by a factor 1 + O(TRAIL_MAX²/r).

**Cut-off and false filter matches.** The replay is cut off after 65n steps (section 1.3, step 7). A divergence between a
live trail and its replay needs a cycle of length at least 65 that one anchor schedule detects and the other misses. A false
64-bit filter match alone cannot cause this, because an escape happens only after an exact key check. Such a cycle is a
within-trail repeat, and it is covered by DELTA_CUT:

  DELTA_CUT := DELTA_LOOP + 2^−100 ≈ 1.0503·10^−16,

where 2^−100 bounds the structural term (Lemma C: cycles of length at least 65 occur at most 2^−230 per step, and S_max ≤
2^126.75 steps). A cut-off replay fails the run, so DELTA_CUT is the failure term of this event. No separate per-trail term is
charged for false matches; such a term would be about 2^70, which is not a bound worth stating.

**Lemma L (no re-entry into an escaped cycle).** By D2, every re-entry into an escaped cycle is a within-trail
coincidence, so its probability over the run is at most DELTA_LOOP.

## 4. Closure rates: Lemma C' (exact under M1)

A new point's successor index is (key >> 1) mod A, uniform over A (M1). A closure of length ℓ is a formal identity: the
ℓ steps' adder terms sum to zero in each index, which needs every index used to occur at least twice (a singleton unit is
nonzero), and the units of each index to cancel in Z[zeta_6]. The unit of a step is sigma = ±lambda^j, uniform over Gamma
by M1, and the sums of elements of Gamma that vanish are the relations of the sixth roots of unity. Under M1:

- **2-cycle:** the index repeats with probability 1/A and the units cancel with probability 1/6. The rate is exactly
  1/(|Gamma| A) = 1/(6A).
- **3-cycle:** all three indices are equal (1/A²), and the three units sum to zero: 12 of the 216 ordered triples of Gamma,
  probability 1/18. The rate is 1/(18 A²).
- **4-cycle:** an index pattern a b a b (1/A²), with both pairs cancelling (1/6 each). The adjacent patterns a a b b and
  a b b a contain a 2-cycle and are pre-empted by its escape. The rate is 1/(36 A²).
- **5 to 64:** the no-singleton estimate of the Lemma C argument, which uses only the uniform and independent index and unit
  data, gives at most 2^−40 for all these lengths together (length 5 alone is about 2^−44).

So the 3..64-cycle entry rate is at most 1/(18A²) + 1/(36A²) + 2^−40 ≤ 0.1/A². The accounting keeps the larger constant
c = 3 (rate 3/A²) as an upper bound; this is conservative. The false 63-bit filter match rate is 2^−63 per step, and each
such match is rejected by the exact check. With these, the conditional closure probability per step is at most
p_c = 1/(6A) + 3/A² + 2^−63 = 2.5438·10^−6.

Under M1 these rates hold conditionally on the past at every processed step, because a closure at a step needs the
new point's class data, which are fresh unless the new class is one the trail has already produced (a within-trail
coincidence, section 3). In the formal test of section 2.6 the 3-cycle entries number 12 against 10.2 predicted, and the 4-cycle entries number 0 against 5.1 predicted; the 4-cycle shortfall is not explained (section 2.6).

## 5. Abandonment: Lemma T' (first-moment bound)

Let n = 64 (CAP_PERIODS + 1) = TRAIL_MAX be the number of steps before a trail is abandoned. Let C_n be the number of closure
steps among the first n steps of a trail. On the trail's coincidence-free event, each fresh point is a distinguished point with
probability theta, independently of the past, because the DP test is a function of the fresh class's key (M1).

*First moment.* For k ≥ 1, by the chain rule each conditional closure bound holds given the whole past, so the probability
that k given steps all close is at most p_c^k, and

  Pr[C_n ≥ k] ≤ E[binom(C_n, k)] ≤ binom(n, k)·p_c^k ≤ (e·n·p_c/k)^k.

*Fresh points.* If C_n < k, at least n − 128k steps are fresh, since each closure wastes at most 128 steps (the cycle's
repeats before detection). None of them is a distinguished point with probability at most (1 − theta)^(n − 128k).

*Abandonment.* Take k = floor(n/c_T) with c_T = 2^16 (A = 2^16). Then e·c_T·p_c = 0.45 ≤ 1/2, so (e n p_c/k)^k ≤ 2^−k, and
128k ≤ n/2^9. Hence

  Pr[trail abandoned, on its coincidence-free event] ≤ 2^−k + exp(−64(1 − 2^−7)) = 2^−k + exp(−63.5) ≈ 2.644·10^−28.

The first term is below 2^−(2^50) and is 0 in double precision. (The earlier model note quoted 2.4·10^−28 for exp(−63.5); the
value is 2.6442·10^−28.)

*Expectation.* Expected abandonments are at most K_TRAILS·Pr[abandon] + DELTA_CUT: the trails on the good event contribute the first
term, and the event of section 3 (within-trail coincidence or cut-off) adds at most DELTA_CUT in total over the run. Markov's inequality then gives
that more than K abandonments occur with probability at most

  DELTA_ABANDON = (K_TRAILS · 2.644·10^−28 + DELTA_CUT) / K ≈ 1.528·10^−17.

## 6. Expected rare work: Lemma M' (every processed step counts, revisits included)

Processed walk steps are at most S_max. Each step's metered rare work is bounded in expectation by the rate of the event
times its worst case, under M1 and the bounds of sections 3 to 5, and the rates hold at every processed step, revisits
included:

| event | rate per walk step | worst case (instructions) |
| --- | ---: | ---: |
| 2-cycle check (true or false predecessor match) | 1/(6A) + 2^−63 | 48,039 |
| anchor check (3..64-cycle or false anchor match), with its look-ahead steps | 3/A² + 2^−63 | 1,618,122 |
| exact-k fallback (canon_k takes the fallback), per canon call | 2^−55 | 256 |
| within-trail coincidence (section 3), per step | |Gamma|·TRAIL_MAX·(2/r + 1/(r − |Gamma|·TRAIL_MAX²/2)) ≈ 1.1e-56 | TRAIL_MAX·(821 + reserve) |

The fallback row needs the rate 2^−55 per canon call. Under M1, extended to the canonical coordinate fractions (section 2.2),
the test of Lemma 1 fails when one of g_0, g_1, g_2 lies in the top 22 values of the 64-bit range, or when two of them lie
within 22 of each other. For uniform fractions this has probability at most (3·23 + 3·45)/2^64 = 204/2^64 < 2^−55.6. The worst
case is the whole canon call with the fallback taken: the forced maximum of the canon_k path with the fallback forced is 253
instructions (static bound, including canon_k_exact at 197 inside it), charged as 256. The batch path calls canon_k once per
walk step, so the rate is per walk step. This row corrects the earlier text, which put the fallback inside f_single's worst
case: that holds for rare-path calls only.

The coincidence row is self-referential, because its worst case contains the reserve. The fixed point is
E_RARE = (E_FIXED + rate·TRAIL_MAX·C_COMMON) / (1 − 64·rate·TRAIL_MAX), with E_FIXED the first three rows. Numerically:

- E_FIXED = 0.12329973978940546 (the first three rows; the fallback row adds 7.1·10^−15);
- E_RARE = 0.12329973978940546 (the coincidence row adds 7.0·10^−34);
- reserve per walk step 64·E_RARE = 7.891183, charge per walk step C_COMMON + reserve = 828.891183.

Markov's inequality gives Pr[metered rare work > 64·S_max·E_RARE] ≤ 1/64, counted in delta. The meter aborts when the
cumulative metered work exceeds that bound.

## 7. Success probability

Let T_steps = K·L. Then tau = T_steps²/(2N) = 0.530300. Under H1, the probability that the first T_steps steps contain no
useful class collision is at most exp(−tau·(1 − eps)), with eps = 4/A + 2·(3/A) + 2·(3·192/A²) = 0.00015286. Its three
margins are:
- the adder penalty, taken as 4/A;
- 2-cycles at a rate of at most 1/A, with 3 non-fresh states each;
- short cycles at a rate of at most 3/A² (section 4), with 192 wasted states each.

Consider a class collision among the first T_steps steps between two different trails.
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
- **Within-trail coincidences and cut-off replays (sections 1.3 and 3)**, including re-entry into an escaped cycle (Lemma L), degenerate pairs and cycles of length at least 65: at most DELTA_CUT = 1.0503·10^−16.
- **Abandoned trails (section 5)**: at most DELTA_ABANDON = 1.528·10^−17.
- **Escape depth tail (section 2.3)**: at most DELTA_DEPTH = 2^(−1.29·10^15), below double precision.
- **Degenerate first match:** exactly 1/r (Lemma D).
- **Work meter abort (step 8):** at most 1/64 (Lemma M').
- **DP_MAX overflow:** the expected number of records is at most DP_MAX/2, so by Chernoff under
  H1 the overflow probability is at most exp(−DP_MAX/12), which is below 2^−60.
- **Exceptional additions:** a batch whose product is zero, or a point that a formula cannot
  handle, aborts the run. This has probability at most 2^−110 under H1.
- **Other terms:** together at most 2^−20 = delta.

Therefore:

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015286)) − DELTA_CUT − DELTA_ABANDON − DELTA_DEPTH − 1/64 − 2^−20 = 0.395898 ≥ 0.39.**

The model term is 0.411524. The budget is 1.0299·√N steps, against 0.9943·√N for success exactly 0.39. The 0.0215 of
slack in the model term covers the sampling error of the toy evidence below.

**Evidence for H1 and the procedure.** `experiment.py` runs the procedure of
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
- The filters compare keys masked to `filter_bits`, with the same distinguished-point-first
  order and look-ahead tests as the kernels. With `check_all_replays`, every recorded
  trail, not only the resolving pair, is replayed and checked against its record.
- Each trial runs L steps per slot plus the finishing phase at the stated tau. It counts
  as a success only if the first key match resolves to a verified x.
- Each trail start is charged as one of the slot's L steps, which is conservative.

<!-- toy-evidence -->
Offline runs (`evidence/trails/`, pooled in `summary-double.json`): 5,200 trials on 5 curves with
N = 2^26.8, 2^29.1, 2^30.8, 2^31.3, 2^32.5 classes, adder tables of 2^11 and 2^16,
16 slots and theta = 2^−5, escape "double".
- Success of the whole procedure at tau = 0.5303: 2,169 of 5,200, 0.4171, Wilson 95% interval
  [0.4038, 0.4306]. The model 1 − exp(−tau) is 0.412.
- Every success was verified ([x]P = Q). Degenerate first matches: 0.
- Production 64-bit filters. Every recorded trail was replayed with a different anchor
  schedule: 4,619,736 replays, 0 that missed the recorded point.
- Lemma R stress (`evidence/trails-stress-double.json`): the same procedure with the filters
  cut to 8 bits over 5,200 trials, so false filter matches were frequent.
  - 4,837,229 trails were replayed with a different anchor schedule, and 0 missed the
    recorded point.
  - 349,511 distinguished points were first reached during a look-ahead.
  - Those runs' success rates are not used: each such look-ahead adds up to 65 points
    that the step budget does not count.
- Replay-order sensitivity (`evidence/lemma-r-sensitivity-double.json`): 200 trials per curve with 6-bit filters
  and 16 adders, every recorded trail replayed. With the current order, 0 of 59,977 replays (j0) and 0 of 95,625 (a3)
  missed the recorded point. With the earlier order (an offline switch), 8 of 59,351 (j0) and 47 of 92,993 (a3)
  did, and that order also abandoned 28 and 89 trails.
- Abandoned trails: 0. Longer fruitless cycles escaped: 6. The finishing phase added at most
  2.6% to the K·L steps.
- 2-cycle entry rate per step:
  - 7.25e-05, 7.74e-05, 7.73e-05, 7.76e-05 with A = 2^11 (1/(6A) = 8.14e-05);
  - 2.03e-06 with A = 2^16 (1/(6A) = 2.54e-06).
<!-- /toy-evidence -->

The declared initial input (`bits` 28, `b` 7, 80 trials, seed 2801, production 64-bit
filters, every trail replayed, escape "double"; `evidence/declared-run-double.json`) runs in a few seconds. It reports
35/80 = 0.4375, Wilson 95% interval [0.334, 0.547], against the model 0.412, with 0 degenerate matches, 0
distinguished points reached during a look-ahead and 0 replays, of 15,059, that missed
their recorded point. With 80 trials the standard error is about 0.05, so it is a
reproducibility check; the success claim rests on the pooled offline trials above.

Limits: the toy curves have N up to about 2^33, against 2^253 in production. There theta is 2^−5,
so start points and the finishing phase are a few percent of the steps; in production
they are below 2^−30 of them. H1 is an extrapolation in N, supported by the size-free
random-mapping argument and the cited literature. It is not proved for this walk.

## 8. Computation

**No reference operation.** Every instruction of a walk step is charged at 1/C with C = 1000, as
COST_MODEL.md's H + R/C with H = 0. No reference operation is charged.

**Charge per walk step.** Two parts, charged separately:
- The common batched path is charged C_COMMON = 821 instructions per walk step: the larger of the plain and
  anchor steps, each measured as the whole-run instruction total difference between runs with one more batch step
  on production inputs, per walk. The anchor batch contains no rare event. On the measured inputs the plain step is
  808.67 and the anchor step 818.66. The one data-dependent branch of the common path, the canon_pick selection
  (k_max = 2), adds at most 2 instructions per walk step, which gives 820.66; C_COMMON = 821 bounds that.
- The rare events, each charged its worst case (section 6, and section 8.1 for the primitives): a 2-cycle check at most
  48,039 instructions, an anchor check, with its look-ahead steps, at most 1,618,122, and a canon fallback at most 256.
  Their expected counts are covered by the Markov reserve 64·E_RARE = 7.891183 per walk step (Lemma M').

The charge per walk step is therefore 828.891183 instructions.

**Cost is uniform over the coin values.** The common batched path is branch-free apart from canon_pick, which adds at most
2 instructions per walk step, and the exact-k fallback, which is metered at its forced maximum (section 6). No other step
reads the coefficient values through a branch, so the instruction count of each walk step is bounded uniformly over all
coin outcomes, and the ledger charges that bound rather than an average. The measured totals on the two input sets (seeds
1 and 2) agree, as this predicts. This paragraph is in section 8 because section 3 is the coincidence analysis.

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
    with Python big-integer arithmetic.

Inline assembly is ordinary scalar AArch64 code. Callgrind counts each executed
instruction exactly as for compiled C, which is how the organizer measured blst's assembly.

**8.1 Worst cases of the rare-path operations.** Each rare event is charged a forced-path maximum, not a sampled
average. The primitives are f_single (adder selection and the unbatched addition), add_single, escape (its least-key
search, its canonical point and its doubling), same_class, point_is_dp, and canon_k with its exact fallback. For each
primitive the data-dependent branches were enumerated from the source and the disassembly. Each conditional correction in a
field routine that is not branch-free was replaced, for the measurement only, by a branch-free routine that costs at
least as much on every path (a bounding build; equality with the original is checked on random and edge inputs).
The remaining control variants (the equal-x and equal-y tests of the addition, the parity of the least member, the
rotation fallback where present, the length of a cycle, and the memcmp outcomes of the library comparison) were forced
by chosen inputs, and the maximum over the variants was taken. Each count was checked to be linear in the number of
calls, and identical for a second seed of inputs. Measured maxima (Callgrind Ir per call, organizer image and flags,
`measurement/forced/summary.json`):

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

The repeat counter (a diagnostic, Section 1) costs at most 51 instructions per call and is inside the overhead.

A 2-cycle check costs at most f_single + same_class + escape(2) + point_is_dp plus a control overhead of at most
1,024 instructions, and an anchor check costs at most f_single + escape(65) plus, for each of its at most 64
look-ahead steps, f_single + point_is_dp + same_class, with the same overhead per step. The overhead bound is far
above the measured control code of these functions, and the sums are above the longest path through the code.
The forced-path harness compiles the primitives in its own translation unit. GCC compiles finv about 36 instructions
more expensive there than in the unmodified step.c, so the harness figures are conservative for every primitive that
calls finv. As an independent check, the longest control-flow path of the unmodified binary, with every conditional
branch free and finv at its measured cost, gives at most 47,193 instructions for a 2-cycle check and 1,547,316 for an
anchor check, both below the charges above (48,039 and 1,618,122; `measurement/forced/summary.json`).
Scalar multiplications and replays use f_single's worst case as their unit, an upper bound for an addition or a
doubling.

**Measurement (H3).** `measurement/step.c` runs the trails, steps, cycle handling, distinguished
points, restarts and trail cap. Its common batched path was measured on production inputs as
whole-run instruction totals (Callgrind `Ir`, GCC 14.2.0 in the organizer's image,
`-O3 -fno-tree-vectorize -mgeneral-regs-only`), with the rare paths out of line so that each
path's cost is separable. The first batch step has no anchor and no rare event. The anchor step is 10.0 instructions per walk step above the plain step,
and its total contains its own rare events.

- Plain step (steps 1 and 2), per walk step: 808.67 instructions (measured).
- Anchor step (steps 63 and 64), per walk step: 818.66 instructions (measured).
- C_COMMON = 821 per walk step: the measured anchor step plus 2 for the canon_pick branch, rounded up. The bound for
  walks that are all long is 819.34, so 821 is conservative.

**Rare-path checks.** Distinguished points and false 64-bit filter matches do not occur in a
production-size sample, so small builds of the same source and replica exercise them
(`evidence/measurement/trail-paths.json`, with 2^6 adders and 2^12 walks):
- DP_BITS = 4: restarts at distinguished points, with 2-cycles;
- DP_BITS = 12 with a one-period cap: abandonments;
- DP_BITS = 6 with the filters cut to 8 bits: 2-cycles, longer cycles, and distinguished points first
  reached during a look-ahead.

In every build the replica matched every final point, every event count, and a digest of every distinguished record, look-ahead records included (`evidence/measurement/trail-paths.json`). The digest also detects the round-5 omission, in which a look-ahead distinguished point was not recorded: the mutant build agrees on every count and differs in the digest (`evidence/measurement/record-mutation-double.json`).

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

Rounded up, **time_log2 = 126.479214**. Preprocessing (adders and initial walk starts) is included. All
processors' work is summed, and elapsed time plays no part. The replay row is charged in every run, with 2^13 operations per
replayed step (section 1.3, step 7).

## 9. Memory and advice

Peak memory is at most the sum of:
- 2·DP_MAX·128 bytes for the distinguished-point table (key, the trail's coins, hashing
  overhead);
- K·576 bytes of walk state (points, keys, coins, trail counter, batch buffers);
- 2^16·6·64 bytes of adder images;
- 16 MiB.

That is 2^75.74996 bytes. Nonuniform advice is 0 bytes.

## 10. Limitations

- H1 is the standard but unproved random-mapping heuristic. It is supported by the
  literature and toy-curve experiments only. The model term is 0.4115, 0.0215 above
  0.39.
- M1 is a premise, not a theorem (section 2.2), and it carries the claim. The algorithm's coins are not a premise: they are
  drawn from the operating system's randomness.
- The restart pool of the measurement harness (1024 points) is not the specified algorithm: restarts there reuse starts. The
  specified restart draws fresh coins and is charged at its cost (section 2.5). The kernel was not changed.
- The measured instances use seeded generators (seed 1 for `reproduce.sh`, the declared run and the toy experiment). The
  algorithm draws its coins from the operating system.
- The replay cut-off (65n steps, REPLAY_FACTOR = 65 in the core) is in the regenerated experiment.py. The declared run is
  unchanged (35/80, 15,059 replays, 0 mismatches). The kernel and the replica do not replay and were not changed, and the
  small builds were not rerun for the cut-off.
- Escape depths of 2^50 or more are covered only by the depth tail (section 2.3), which is a first-moment bound under M1. The certificate of section 2.4 is exact for depths below 2^50.
- **The argument is not established at production size.** The toy runs use theta = 2^−5, and the formal test uses 2^6
  adders and theta = 2^−6. Neither tests M1 at production theta or adder count.
- Instruction counts come from one implementation and compiler. The common path was measured on two production input sets (section 8), and its worst case is bounded there. The rare operations are charged forced-path maxima (section 8.1); the bounding builds are checked for dominance on the forced paths, not proved for every input.
- The exact-k fallback rate (2^−55 per canon call) follows from M1 extended to the coordinate fractions. It is not measured: no fallback
  occurred in the measured steps, which does not test that rate.
- The field routines were tested, not formally verified. The replica check covers 2^19
  production steps per input set.
- About 2^67.7 distinguished points and 2^126.5 work are needed. This is not a practical
  attack and does not weaken secp256k1.
- The toy experiment runs the whole procedure, including resolution, but only at toy
  scale: theta = 2^−5, 16 slots, and N far below production. Its trail starts and its
  finishing phase are a few percent of the steps there, and negligible in production.
- Distinguished points and false 64-bit filter matches never occur in a
  production-parameter measurement sample. The kernel's restart, abandonment and
  false-match paths were checked on small builds of the same source.
- The work meter's success loss, 1/64, is a Markov bound on the expected metered work under M1 (Lemma M').

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
echo '{"bits":28,"trials":80,"seed":2801,"b":7,"adder_bits":11,"dp_bits":5,"walks":16,"tau_e4":5303,"filter_bits":64,"check_all_replays":1,"escape":"double"}' | python3 -I submissions/secp256k1-ecdlp/rho-sigma6-asm/experiment.py
```

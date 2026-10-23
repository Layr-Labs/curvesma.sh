# secp256k1-ecdlp: parallel Pollard rho on order-6 classes with a sigma walk, all arithmetic in instructions

Claim: for uniform x in [0,r) and Q=[x]P on secp256k1, the algorithm below returns x
with probability at least 0.39. It uses at most **2^126.479214** target-group-operation equivalents in total (score 63239607/500000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^75.749960 bytes. Preprocessing is 2^47.131242 units, and no nonuniform advice is used. This is a classical cold attack on
one target.

The success bound rests on the disclosed heuristic H1, the standard random-mapping model
for r-adding walks on automorphism classes, and on the fresh-class premise H2, from which Lemma C and the metered-work rates
follow (section 2). The cost rests on
H3, average instruction counts measured with the organizer's toolchain, and every run is
charged by a work meter with a Markov reserve. This is a carefully
accounted generic attack. It is not a practical break.

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

8. **Work meter.** Every run is charged at most W_max = S_max·828.891183 + G
   instructions, with G the fixed global charges (initial starts, trail restarts, distinguished-point
   records, solve and verify, input validation, adder table).
   - Every walk step is charged C_COMMON = 821 instructions, which bounds the common path
     including its batch inversion, and each rare event is charged its worst case.
   - The meter aborts the run when the total metered rare work exceeds 64 times its
     expectation bound (Lemma M). That happens with probability at most 1/64, counted in delta.
   - The run also fails if more than DP_MAX records arrive or more than K trails are abandoned.

**Escape.** A cycle's escape goes to 2·canon(M), where M is the cycle's least-key member. canon(M) is a point of M's class, so 2·canon(M) is a point of the class of 2M. The escape uses no table and no stored point: it is one point doubling, whose cost is charged at its forced-path maximum (Section 3.1). The escape point depends only on the cycle's members, as Lemma R requires. M is a fresh class (H2), and by H2(ii) so is the escape point.

**Escape repeats.** The escape point of a cycle depends only on its least-key member, so it is the same every time the trail escapes from that cycle. A repeated escape point therefore means that the trail has re-entered a cycle it already escaped from. Lemma L bounds this event.

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
- Along a replay, b changes as b ← u·b + d, where u is a unit (order-6 classes {±1, ±λ, ±λ²}, up to sign)
  or twice one, and u, d depend only on the path.

So the canonical b at the distinguished point is U_i·b_ti + D_i with U_i invertible mod r,
and P[b1 = b2] = 1/r. ∎

**Lemma M (metered rare work, Markov reserve).** The metered events are the rare-path calls. Each is charged its
forced-path worst case (Section 3.1), and the rate of each under H2 is an unconditional per-step probability:

| event | rate per walk step | worst case (instructions) |
| --- | --- | --- |
| 2-cycle check (exact check; true or false predecessor match) | 1/(6A) + 2^−63 | 48,039 |
| anchor check (true or false anchor match), with all its look-ahead steps (at most 64) | 3/A^2 + 2^−63 | 1,618,122 |

The expected metered work per walk step is E_RARE = 0.123300, the sum of rate × cost over the rows. Under H2 the
unconditional per-step rates hold at every step, so by linearity of expectation E[X] ≤ S·E_RARE for the total metered
work X of a run. Markov's inequality gives P[X > 64·S·E_RARE] ≤ 1/64. The meter aborts when the cumulative metered
work exceeds 64·S_max·E_RARE; that work never decreases, so this is the event {X > 64·S_max·E_RARE} and its
probability is at most 1/64, counted in delta. The reserve per walk step is 64·E_RARE = 7.891183. The expected
numbers of rare events are covered by this global reserve. The exact-k rotation fallback is inside the worst case of
f_single, so it needs no row of its own. Every escape is one of these events: an escape
happens only after a 2-cycle check or an anchor check, and its doubling is inside that event's worst case. The expected
escape work is therefore covered by E_RARE, and no escape is uncharged. Every metered type has a rate derived from H2, so none needs a hard counter.


## 2. Success probability

Let T_steps = K·L. Then tau = T_steps^2/(2N) = 0.530300.

**H1 (heuristic, disclosed).** The class map [W] → [canon(W) + R_h([W])], with uniform
random adders, behaves like a random mapping on the N classes, apart from the r-adding
penalty. That penalty is a factor of about (1 − 1/(6A)) in the birthday exponent (Teske
2001; Bernstein–Lange 2012; Bos–Kleinjung–Lenstra 2010). Under H1, the probability that
the first T_steps steps contain no useful class collision is at most
exp(−tau·(1 − eps)). Here eps = 4/A + 2·(3/A) + 2·(3·192/A^2) = 0.00015286. Its three
margins are:
- the adder penalty, taken as 4/A;
- 2-cycles at a rate of at most 1/A (Lemma C gives 1/(6A)), with 3 non-fresh states each;
- short cycles at a rate of at most 3/A^2 (Lemma C), with 192 wasted states each.

**H2 (fresh-class premise).** Along a trail, a class is fresh when the walk produces it for the first time in that
trail. (i) A fresh class's point is uniform among the curve's points of that class, independent of the classes
produced earlier and of the adders; hence its adder index, canonical rotation and sign, key bits and x-coordinate are
independent and uniform. (ii) An escape produces 2·canon(M), where M is the least-key member of a cycle of fresh classes. Doubling permutes the classes of the group, so the class of 2M is fresh whenever the class of M is: by (i) it lies among the t classes already produced in the trail with probability at most t/N, and (i) applies to every later class. A revisit
within a trail either closes a cycle, which the walk detects, or is a closure of some length ℓ that Lemma C bounds for
every ℓ; Lemma C is derived from H2(i) alone, and no other result of this proof is used in its derivation. Revisits across trails are collisions, which the solver resolves. H2 is a premise
about the random walk, not a theorem, and the bounds below use it only through the per-step rates.

**Per-step rates under H2.** These follow from Lemma C below and from H2(i). Each is an unconditional probability for
one walk step.
- *Closures.* Lemma C: a 2-cycle closes with probability at most 1/(6A), and a cycle of length 3 to 64 with
  probability at most 3/A^2 in total. The accounting uses this bound.
- *False filter matches.* A fresh class has uniform key bits, and a stored key is fixed, so a 63-bit filter
  matches it falsely with probability 2^−63 (the even key has bit 0 zero; bit 0 of the even key is zero, so 63 bits are random).
- *Escapes.* The escape rate ρ_e = 1/(6A) + 3/A^2 + 2·2^−63 bounds the probability that a step escapes (Lemma T).

**Lemma C (fruitless-cycle rates, under H2).** A fruitless cycle is a cycle whose formal
coefficient sum vanishes, so it closes whatever the adders are. At any step, the
probability that the walk enters one of length ℓ is at most:
- 1/(6A) for ℓ = 2;
- 3/A^2 in total over 3 ≤ ℓ ≤ 64;
- 2^−300 for each ℓ ≥ 65.

*Proof.* Let the steps of the cycle use adder indices i_1..i_ℓ with canonical units u_k. By H2 these are independent, the indices uniform on the adders and the units uniform on the canonical units.
- The adders' coefficients are independent and uniform, so the cycle closes for every
  choice of adders only if, for every index i that occurs, the unit-weighted sum of its
  occurrences vanishes mod r. Any other closure has probability at most ℓ/r; those are the
  self-collisions in the failure terms below.
- An index that occurs once cannot cancel. For an index that occurs m ≥ 2 times, given
  the other units the last one must take a single value, with probability at most
  1/6 ≤ 1/2.
- So the probability is at most E[2^−D·1{no index occurs once}], with D the number of
  distinct indices among ℓ independent uniform draws from A.

The cases:
- ℓ = 2: the indices agree (probability 1/A) and the units cancel (1/6).
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

**Lemma T (trails are not trapped).** Let a trail run m steps, let E_m be its number of escapes and D_m the number of
distinct classes it produces. Each escape follows a cycle of at most 65 steps, so D_m ≥ m − 66·E_m. Under H2 each escape point starts a fresh segment of the walk, and the escape indicator at a step depends on at most 65 consecutive classes after it, so the indicators form a 65-dependent sequence with mean at most ρ_e. Splitting the steps into 66 interleaved subsequences makes each one independent, and
Chernoff's bound gives P[E_m ≥ 2ρ_e·m] ≤ 66·exp(−ρ_e·m/198). By H2(i) the key bits above the walk index and the
distinguished-point bits are uniform and independent of them, so each fresh class is a distinguished point with
probability θ, independently of the others. A trail is abandoned only if its first D_m fresh classes contain no
distinguished point, which has probability (1−θ)^{D_m}. For m = CAP, with 132·ρ_e < 1,

  P[abandoned] ≤ (1−θ)^{CAP(1−132·ρ_e)} + 66·exp(−ρ_e·CAP/198) ≤ exp(−64(1−132·ρ_e)) + 2^−100.

With ρ_e = 2.54383e-06 (the sum of the closure and false-match rates) this is exp(−64(1−132·ρ_e)) = 1.639e-28. A re-entry into an escaped cycle is a repeated escape point; Lemma L bounds it, and its work is metered as any escape (Lemma M).
Multiple short-cycle escapes therefore reduce to independent fresh segments: each escape starts a fresh point, and each fresh point is a
distinguished point with probability θ. The expected number of abandoned trails is at most (DP_MAX + 2K)·exp(−64(1−132·ρ_e)) = 1.639e-28,
so the run fails by abandonment with probability at most that expectation divided by K (Markov), counted in delta.

**Lemma L (no re-entry into an escaped cycle).** Let a trail have produced t classes before one of its escapes. A repeated escape point is an escape point produced earlier in the trail, and the escape point of a cycle depends only on its least-key member. So the first repeat needs the new escape point's class to be one of those t classes, which by H2(ii) has probability at most t/N. A trail has fewer than TRAIL_MAX steps, so each escape is a repeat with probability at most TRAIL_MAX/N. Over the run, with at most S_max·ρ_e escapes in expectation,

  P[some trail re-enters an escaped cycle] ≤ S_max·ρ_e·TRAIL_MAX/N ≤ 2^−79.

The bound is counted in delta. Outside this event every escape point is a class the trail has not produced before, which is the premise of Lemma T. On the event the walk is not restarted: the trail still ends at its first distinguished point or at the anchor cap, and its work is metered like any other (Lemma M).

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
- **Abandoned trail.** By Lemma T, a trail is abandoned with probability at most exp(−64(1−132·ρ_e)) = 1.639e-28, plus a negligible
  term. The expected number of abandonments is at most (DP_MAX + 2K) times that, so more than K abandonments
  occur with probability at most that expectation divided by K (Markov), counted in delta.
- **Re-entry into an escaped cycle (Lemma L):** at most 2^−79, counted in delta.
- **Degenerate first match:** exactly 1/r (Lemma D).
- **Work meter abort (step 8):** at most 1/64 (Lemma M).
- **DP_MAX overflow:** the expected number of records is at most DP_MAX/2, so by Chernoff under
  H1 the overflow probability is at most exp(−DP_MAX/12), which is below 2^−60.
- **Exceptional additions:** a batch whose product is zero, or a point that a formula cannot
  handle, aborts the run. This has probability at most 2^−110 under H1.

Together these are below delta = 2^−20. Therefore:

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015286)) − 2^−20 = 0.411524 ≥ 0.39.**

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

## 3. Computation

**No reference operation.** Every instruction of a walk step is charged at 1/C with C = 1000, as
COST_MODEL.md's H + R/C with H = 0. No reference operation is charged.

**Charge per walk step.** Two parts, charged separately:
- The common batched path is charged C_COMMON = 821 instructions per walk step: the larger of the plain and
  anchor steps, each measured as the whole-run instruction total difference between runs with one more batch step
  on production inputs, per walk. The anchor batch contains no rare event. On the measured inputs the plain step is
  808.67 and the anchor step 818.66. The one data-dependent branch of the common path, the canon_pick selection
  (k_max = 2), adds at most 2 instructions per walk step, which gives 820.66; C_COMMON = 821 bounds that.
- The rare events, each charged its forced-path worst case (Section 3.1): a 2-cycle check at most
  48,039 instructions and an anchor check, with its look-ahead steps, at most 1,618,122. Their expected
  counts are covered by the Markov reserve 64·E_RARE = 7.891183 per walk step (Lemma M).

The charge per walk step is therefore 828.891183 instructions.

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

**3.1 Worst cases of the rare-path operations.** Each rare event is charged a forced-path maximum, not a sampled
average. The primitives are f_single (adder selection and the unbatched addition), add_single, escape (its least-key
search, its canonical point and its doubling), same_class, and point_is_dp. For each primitive
the data-dependent branches were enumerated from the source and the disassembly. Each conditional correction in a
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
| --- | --- | ---: |
| walk steps | 828.891183·S_max / C | 126.479214 |
| adder table | 2^16 adders × 2 scalar mults × 768 unbatched ops × 2^15 / C | 31.13 |
| initial walk starts | K × 2 scalar mults × 768 × 2^15 / C | 47.13 |
| trail restarts | (DP_MAX + K) × 2 scalar mults × 768 × 2^15 / C | 82.88 |
| distinguished points | DP_MAX = 2^67.75 entries × 2^20 instructions / C | 77.78 |
| collision replay | 2 trails × TRAIL_MAX steps × 2^12 unbatched ops × (2^15 + 2^10) / C | 83.61 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 15.68 |
| input validation | 2^24 instructions | 14.03 |

Rounded up, **time_log2 = 126.479214**. Preprocessing (adders and initial walk starts) is included. All
processors' work is summed, and elapsed time plays no part.

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
- Instruction counts come from one implementation and compiler. The common path was measured on two production input sets (Section 3.2), and its worst case is bounded there. The rare operations are charged forced-path maxima (Section 3.1); the bounding builds are checked for dominance on the forced paths, not proved for every input.
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
- The work meter's success loss, 1/64, is a Markov bound on the expected metered work under H2 (Lemma M).

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
echo '{"bits":28,"trials":80,"seed":2801,"b":7,"adder_bits":11,"dp_bits":5,"walks":16,"tau_e4":5303,"filter_bits":64,"check_all_replays":1}' | python3 -I submissions/secp256k1-ecdlp/rho-sigma6-asm/experiment.py
```

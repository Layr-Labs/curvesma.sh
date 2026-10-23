# bls12-381-g1-ecdlp-plain: parallel Pollard rho on order-6 classes with a sigma walk, all arithmetic in instructions

Claim: for uniform x in [0,r) and T=[x]P0 in G1, the algorithm below returns x with
probability at least 0.39. It uses at most **2^126.390705** target-group-operation equivalents in total (score 25278141/200000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^75.178504 bytes. Preprocessing is 2^46.934997 units, and no nonuniform advice is used. This is a classical cold attack on
one target.

The success bound rests on the disclosed heuristic H1, the standard random-mapping model
for r-adding walks on automorphism classes, and on the fresh-class premise H2, from which Lemma C and the metered-work rates
follow (section 2). The cost rests on
the instruction counts in H3, measured with the organizer's toolchain. This is a
carefully accounted generic attack. It is not a practical break.

**Cost route.** The walk does not use the 6M reference operation. Its five field
multiplications and one squaring per step are the package's own 6×64-bit Montgomery
multiplication in GCC inline AArch64 assembly (`measurement/bls_field.h`, the same
convention as blst's `mul_mont_384`). Every instruction of the step is charged at 1/C
with C = 2000 (H = 0 in COST_MODEL's H + R/C). COST_MODEL.md: "Faster arithmetic and
lookup algorithms are welcome; justify their actual work". A walk step is charged 2328.10 instructions, 1.164 units, rather than 1 + R/C. The same walk with
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
- `evidence/measurement/` holds the round-5 production totals and the checks of the
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
     walk. Both of its points have been tested, and the walk replaces new by canon(M) + E_j,
     where M is the element with the smaller key.
   - *3..64-cycles.* Unless a 2-cycle escape just happened, if the new key matches the
     anchor filter, the walk iterates f from the new point, unbatched, up to 64 times,
     testing each look-ahead point first.
     - A distinguished look-ahead point s ends the trail at s, which is then the trail's record, since the walk would reach
       it next.
     - If instead the new point's class returns, every point of that cycle has been
       tested. The walk escapes from canon(M) + E_j, with M the cycle element of least key.
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
   - The replays give the common canonical point as [a1]P0 + [b1]T = [a2]P0 + [b2]T. If
     b1 ≠ b2, output x = (a2 − a1)/(b1 − b2) mod r after checking [x]P0 = T. If b1 = b2
     the run fails. By Lemma D this has probability 1/r.

8. **Work meter.** Every run is charged at most W_max = S_max·2328.099655 + G
   instructions, with G the fixed global charges (initial starts, trail restarts, distinguished-point
   records, solve and verify, input validation, adder table).
   - Every walk step is charged C_COMMON = 2263 instructions, which bounds the common path
     including its batch inversion, and each rare event is charged its worst case.
   - The meter aborts the run when the total metered rare work exceeds 64 times its
     expectation bound (Lemma M). That happens with probability at most 1/64, counted in delta.
   - The run also fails if more than DP_MAX records arrive or more than K trails are abandoned.

**Escape table.** A cycle's escape goes to canon(M) + E_j. M is the cycle's least-key member, E is a separate table
of 2^ADDER_BITS random points that are not walk adders, and j = (key(M) >> s) mod 2^ADDER_BITS with s = ADDER_BITS + 1 + DP_BITS. The
walk index uses bits 1 to ADDER_BITS of the key, and the distinguished-point test uses the DP_BITS bits below s, so j is independent of
both. A cycle's points differ by sums of walk adders, so E_j lands in the cycle only if a sum of walk adders equals a
table point, which has negligible probability. (An escape that added a walk adder would land in its own cycle with
probability about 2/NA, and the walk would then loop; the separate table avoids this.) M is a fresh class (H2), so j
and canon(M) are fresh, and so is the escape point. The escape point depends only on the cycle's members, as Lemma R
requires.

**Merged trails stay merged.** Outside fruitless cycles the next class depends only on the
current class. Escapes fire only after an exact check, and they exit from canon(M) + E_j,
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
  and the escape point canon(M) + E_j depends only on C's set of classes.
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
| 2-cycle check (exact check; true or false predecessor match) | 1/(6A) + 2^−63 | 328,232 |
| anchor check (true or false anchor match), with all its look-ahead steps (at most 64) | 3/A^2 + 2^−63 | 10,861,354 |

The expected metered work per walk step is E_RARE = {E_RARE:.6f}, the sum of rate × cost over the rows. Under H2 the
unconditional per-step rates hold at every step, so by linearity of expectation E[X] ≤ S·E_RARE for the total metered
work X of a run. Markov's inequality gives P[X > 64·S·E_RARE] ≤ 1/64. The meter aborts when the cumulative metered
work exceeds 64·S_max·E_RARE; that work never decreases, so this is the event {X > 64·S_max·E_RARE} and its
probability is at most 1/64, counted in delta. The reserve per walk step is 64·E_RARE = {RESERVE:.6f}. The expected
numbers of rare events are covered by this global reserve. The exact-k rotation fallback is inside the worst case of
f_single, so it needs no row of its own. Every metered type has a rate derived from H2, so none needs a hard counter.


## 2. Success probability

Let T_steps = K·L be the main-phase steps, each walk's first L steps. Then
tau = T_steps^2/(2N) = 0.530300.

**H1 (heuristic, disclosed).** The class map [W] → [canon(W) + R_h([W])], with uniform
random adders, behaves like a random mapping on the N classes, apart from the known
r-adding penalty. An r-adding walk loses about a factor (1 − 1/(m·A)) in the birthday
exponent, with A the number of adders and m = 6 the class size (Teske 2001;
Bernstein–Lange 2012; Bos–Kleinjung–Lenstra 2010 eq. (2) for m = 2).

Under H1, the probability that the first T_steps steps contain no useful class collision
is at most exp(−tau·(1 − eps)). Here eps = 4/A + 2·(3/A) + 2·(3·192/A^2) = 0.00015286. It
collects three margins:

- **Adder penalty, 4/A.** This is 24 times the 1/(6A) estimate.
- **Fruitless 2-cycles, 2·(3/A).** Lemma C bounds the entry rate by 1/A per step, and each
  event yields at most 3 non-fresh states.
- **Short cycles, 2·(3·192/A^2).** Lemma C bounds the rate by 3/A^2, and detection wastes at
  most 192 states per event.

**H2 (fresh-class premise).** Along a trail, a class is fresh when the walk produces it for the first time in that
trail. (i) A fresh class's point is uniform among the curve's points of that class, independent of the classes
produced earlier and of the adders; hence its adder index, canonical rotation and sign, key bits and x-coordinate are
independent and uniform. (ii) An escape produces canon(M) + E_j, where M is the least-key member of a cycle of
fresh classes and E is the escape table. By (i) key(M) is uniform, so the escape point is fresh, and (i) applies to every later class. Revisits
within a trail therefore either close a cycle of length at most 65, which the walk detects, or have probability
below 2^−300 per length (Lemma C). Revisits across trails are collisions, which the solver resolves. H2 is a premise
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
distinct classes it produces. Each escape follows a cycle of at most 65 steps, so D_m ≥ m − 66·E_m. Under H2 the
escape indicator at a step depends on at most 65 consecutive fresh classes, so the indicators form a 65-dependent
sequence with mean at most ρ_e. Splitting the steps into 66 interleaved subsequences makes each one independent, and
Chernoff's bound gives P[E_m ≥ 2ρ_e·m] ≤ 66·exp(−ρ_e·m/198). By H2(i) the key bits above the walk index and the
distinguished-point bits are uniform and independent of them, so each fresh class is a distinguished point with
probability θ, independently of the others. A trail is abandoned only if its first D_m fresh classes contain no
distinguished point, which has probability (1−θ)^{D_m}. For m = CAP, with 132·ρ_e < 1,

  P[abandoned] ≤ (1−θ)^{CAP(1−132·ρ_e)} + 66·exp(−ρ_e·CAP/198) ≤ exp(−64(1−132·ρ_e)) + 2^−100.

With ρ_e = {RHO_E:.6g} (the sum of the closure and false-match rates) this is exp(−64(1−132·ρ_e)) = 1.639e-28. Multiple short-cycle escapes
therefore reduce to independent fresh segments: each escape starts a fresh point, and each fresh point is a
distinguished point with probability θ. The expected number of abandoned trails is at most (DP_MAX + 2K)·exp(−64(1−132·ρ_e)) = 1.639e-28,
so the run fails by abandonment with probability at most that expectation divided by K (Markov), counted in delta.

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
- **Degenerate first match:** exactly 1/r (Lemma D).
- **Work meter abort (step 8):** at most 1/64 (Lemma M).
- **DP_MAX overflow:** the expected number of records is at most DP_MAX/2, so by Chernoff under
  H1 the overflow probability is at most exp(−DP_MAX/12), which is below 2^−60.
- **Exceptional additions:** a batch whose product is zero, or a point that a formula cannot
  handle, aborts the run. This has probability at most 2^−110 under H1.

Together these are below delta = 2^−20. Therefore:

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015286)) − 2^−20 = 0.411523 ≥ 0.39.**

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
- The filters compare keys masked to `filter_bits`, with the same distinguished-point-first
  order and look-ahead tests as the kernels. With `check_all_replays`, every recorded
  trail, not only the resolving pair, is replayed and checked against its record.
- Each trial runs L steps per slot plus the finishing phase at the stated tau. It counts
  as a success only if the first key match resolves to a verified x.
- Each trail start is charged as one of the slot's L steps, which is conservative.

<!-- toy-evidence -->
Offline runs (`evidence/trails/`, pooled in `summary.json`): 5,200 trials on 5 curves with
N = 2^28.8, 2^30.0, 2^35.5, 2^35.5, 2^37.7 classes, adder tables of 2^11 and 2^16,
16 slots and theta = 2^−5.
- Success of the whole procedure at tau = 0.5303: 2,135 of 5,200, 0.4106, Wilson 95% interval
  [0.397, 0.424]. The model 1 − exp(−tau) is 0.412.
- Every success was verified ([x]P = Q). Degenerate first matches: 0.
- Production 64-bit filters. Every recorded trail was replayed with a different anchor
  schedule: 19,284,465 replays, 0 that missed the recorded point.
- Lemma R stress (`evidence/trails-stress.json`): the same procedure with the filters
  cut to 8 bits over 5,200 trials, so false filter matches were frequent.
  - 19,985,272 trails were replayed with a different anchor schedule, and 0 missed the
    recorded point.
  - 1,466,137 distinguished points were first reached during a look-ahead.
  - Those runs' success rates are not used: each such look-ahead adds up to 65 points
    that the step budget does not count.
- Abandoned trails: 1. Longer fruitless cycles escaped: 18. The finishing phase added at most
  1.3% to the K·L steps.
- 2-cycle entry rate per step:
  - 7.74e-05, 7.56e-05, 7.67e-05, 7.60e-05 with A = 2^11 (1/(6A) = 8.14e-05);
  - 2.32e-06 with A = 2^16 (1/(6A) = 2.54e-06).
<!-- /toy-evidence -->

The declared initial input (`bits` 28, `b` 7, 80 trials, seed 2801, production 64-bit
filters, every trail replayed) runs in a few seconds. It uses prime-order y^2 = x^3 + 7
curves, which exercise the same order-6 class code. For b = 4, fields of at most 28 bits
give subgroups of only 2^19 to 2^23 classes, where each slot's final trail adds 20 to 40%
to K·L, so success runs high; the offline b = 4 runs above use larger fields. It reports
32/80 = 0.400, Wilson 95% interval [0.300, 0.509], against the model 0.412, with 0
degenerate matches, 0 distinguished points reached during a look-ahead and 0 replays, of 15,629, that missed their recorded point. With 80 trials the standard error is about 0.05,
so it is a reproducibility check; the success claim rests on the pooled offline trials
above.

Limits: the toy curves have N up to about 2^38, against 2^252 in production. There theta is 2^−5,
so start points and the finishing phase are a few percent of the steps; in production
they are below 2^−30 of them. H1 is an extrapolation in N, supported by the size-free
random-mapping argument and the cited literature. It is not proved for this walk.

## 3. Computation

**No reference operation.** Every instruction of a walk step is charged at 1/C with C = 2000, as
COST_MODEL.md's H + R/C with H = 0. No reference operation is charged.

**Charge per walk step.** Two parts, charged separately:
- The common batched path is charged C_COMMON = 2263 instructions per walk step: the larger of the plain and
  anchor steps, each measured as the whole-run instruction total difference between runs with one more batch step
  on production inputs, per walk. The anchor batch contains one two-cycle escape, which the charge counts twice (conservative).
- The rare events, each charged its forced-path worst case (Section 3.1): a 2-cycle check at most
  328,232 instructions and an anchor check, with its look-ahead steps, at most 10,861,354. Their expected
  counts are covered by the Markov reserve 64·E_RARE = 53.908717 per walk step (Lemma M).

The charge per walk step is therefore 2316.908717 instructions.

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
    on the production input set.
- blst is still used for the one inversion per batch (`blst_fp_inverse`) and for input
  conversion. Those instructions are counted like any others.

Inline assembly is ordinary scalar AArch64 code. Callgrind counts each executed
instruction exactly as for compiled C, which is how the organizer measured blst's assembly.

**3.1 Worst cases of the rare-path operations.** Each rare event is charged a forced-path maximum, not a sampled
average. The primitives are f_single (adder selection and the unbatched addition), add_single, escape (its least-key
search, its canonical point and its addition from the escape table), same_class, and point_is_dp. For each primitive
the data-dependent branches were enumerated from the source and the disassembly. Each conditional correction in a
field routine that is not branch-free was replaced, for the measurement only, by a branch-free routine that costs at
least as much on every path (a bounding build; equality with the original is checked on random and edge inputs).
The remaining control variants (the equal-x and equal-y tests of the addition, the parity of the least member, the
rotation fallback where present, the length of a cycle, and the memcmp outcomes of the library comparison) were forced
by chosen inputs, and the maximum over the variants was taken. Each count was checked to be linear in the number of
calls, and identical for a second seed of inputs. Measured maxima (Callgrind Ir per call, organizer image and flags,
`measurement/forced/summary.json`):

| primitive | worst case (Ir) |
| --- | ---: |
| f_single (adder selection and addition) | 163,117 |
| add_single | 162,441 |
| same_class | 210 |
| point_is_dp | 98 |
| escape of a 2-cycle | 163,783 |
| escape of a 65-cycle | 172,477 |

A 2-cycle check costs at most f_single + same_class + escape(2) + point_is_dp plus a control overhead of at most
1,024 instructions, and an anchor check costs at most f_single + escape(65) plus, for each of its at most 64
look-ahead steps, f_single + point_is_dp + same_class, with the same overhead per step. The overhead bound is far
above the measured control code of these functions, and the sums are above the longest path through the code.
Scalar multiplications and replays use f_single's worst case as their unit, an upper bound for an addition or a
doubling.

The unbatched inversion is blst_fp_inverse. It computes a constant-time inverse, checks it with one multiplication, and falls back to a Fermat chain of 133,370 instructions if the check fails. No tested input reached the fallback, but no proof excludes it, so the worst cases above include it on every inversion.

**Measurement (H3).** `measurement/step.c` runs the trails, steps, cycle handling, distinguished
points, restarts and trail cap. Its common batched path was measured on production inputs as
whole-run instruction totals (Callgrind `Ir`, GCC 14.2.0 in the organizer's image,
`-O3 -fno-tree-vectorize -mgeneral-regs-only`), with the rare paths out of line so that each
path's cost is separable. The first batch step has no anchor and no rare event. The anchor step is 13.7 instructions per walk step above the plain step,
and its total contains its own rare events.

- Plain step (steps 1 and 2), per walk step: 2248.9 instructions.
- Anchor step (steps 63 and 64), per walk step: 2262.6 instructions.
- C_COMMON = 2263 per walk step.

**Rare-path checks.** Distinguished points and false 64-bit filter matches do not occur in a
production-size sample, so small builds of the same source and replica exercise them
(`evidence/measurement/trail-paths.json`, with 2^6 adders and 2^12 walks):
- DP_BITS = 4: restarts at distinguished points, with 2-cycles;
- DP_BITS = 12 with a one-period cap: abandonments;
- DP_BITS = 6 with the filters cut to 8 bits: 2-cycles, longer cycles, and distinguished points first
  reached during a look-ahead.

In every build the replica matched every final point, every event count, and a digest of every distinguished record, look-ahead records included (`evidence/measurement/trail-paths.json`).

**Total.** S_max = K·(L + TRAIL_MAX) ≤ 2^126.178504 walk steps.

| Component | Bound (units) | log2 |
| --- | --- | ---: |
| walk steps | 2328.099655·S_max / C | 126.397657 |
| adder table | 2^16 adders × 2 scalar mults × 768 unbatched ops × 2^15 / C | 30.62 |
| initial walk starts | K × 2 scalar mults × 768 × 2^15 / C | 44.62 |
| trail restarts | (DP_MAX + K) × 2 scalar mults × 768 × 2^15 / C | 81.80 |
| distinguished points | DP_MAX = 2^67.18 entries × 2^20 instructions / C | 76.21 |
| collision replay | 2 trails × TRAIL_MAX steps × 2^12 unbatched ops × (2^15 + 2^10) / C | 83.08 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 15.03 |
| input validation | 2^24 instructions | 15.03 |

Rounded up, **time_log2 = 126.397657**. Preprocessing (adders and initial walk starts) is included. All
processors' work is summed, and elapsed time plays no part.

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
  is 2^17 steps on one production input set.
- The field multiplication is the package's own assembly, tested against Python big
  integers but not formally verified.
- The variable-time paths (branchy subtractions, comparisons, the exact fallback) were
  measured, not bounded analytically. Lemma 1 makes the fallback's rarity analytic.
- About 2^67 distinguished points and 2^126.4 work are needed. This is not a practical
  attack and does not weaken BLS12-381's security assumptions.
- The pairing (MOV) route to GT is not used.
- The toy experiment runs the whole procedure, including resolution, but only at toy
  scale: theta = 2^−5, 16 slots, and N far below production. Its trail starts and its
  finishing phase are a few percent of the steps there, and negligible in production.
- Distinguished points and false 64-bit filter matches never occur in a
  production-parameter measurement sample. The kernel's restart, abandonment and
  false-match paths were checked on small builds of the same source.
- The work meter's success loss, 1/64, is a Markov bound on the expected metered work under H2 (Lemma M).

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
echo '{"bits":28,"trials":80,"seed":2801,"b":7,"adder_bits":11,"dp_bits":5,"walks":16,"tau_e4":5303,"filter_bits":64,"check_all_replays":1}' | python3 -I submissions/bls12-381-g1-ecdlp-plain/rho-sigma6-asm/experiment.py
# offline evidence runs used larger inputs, e.g. {"bits":40,"trials":700,"seed":4001}
```

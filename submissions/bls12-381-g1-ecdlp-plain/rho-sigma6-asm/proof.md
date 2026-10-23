# bls12-381-g1-ecdlp-plain: parallel Pollard rho on order-6 classes with a sigma walk, all arithmetic in instructions

Claim: for uniform x in [0,r) and T=[x]P0 in G1, the algorithm below returns x with
probability at least 0.39. It uses at most **2^126.390705** target-group-operation equivalents in total (score 25278141/200000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^75.178504 bytes. Preprocessing is 2^46.934918 units, and no nonuniform advice is used. This is a classical cold attack on
one target.

The success bound rests on the disclosed heuristics H1 (the standard random-mapping model for r-adding walks on
automorphism classes), M1 (class data, extended to the canonical coordinate fractions) and M2 (the table coefficients,
fixed by the seed random.Random(1) of the specification, behave as uniform values independent of x), stated in section 2.2. The cost rests on the instruction counts in H3 (a measurement
premise), measured with the organizer's toolchain. This is a carefully accounted generic attack. It is not a practical break.

**Status of the argument.** The argument is not established at production size. M1 is the heuristic that carries the
claim. The toy runs (theta = 2^-5) cannot test it at production theta (section 5).

**Cost route.** The walk does not use the 6M reference operation. Its five field
multiplications and one squaring per step are the package's own 6×64-bit Montgomery
multiplication in GCC inline AArch64 assembly (`measurement/bls_field.h`, the same
convention as blst's `mul_mont_384`). Every instruction of the step is charged at 1/C
with C = 2000 (H = 0 in COST_MODEL's H + R/C). COST_MODEL.md: "Faster arithmetic and
lookup algorithms are welcome; justify their actual work". A walk step is charged 2316.91 instructions, 1.158 units, rather than 1 + R/C. The same walk with
the reference operation scores 126.367538 (an earlier version, not included).

What is new relative to the usual order-6 walk is in section 1.2. The walk never moves
a point to its canonical representative, so it never pays the multiplication by beta.
It stores the actual point W and adds sigma(R_i), where sigma is the automorphism that
maps the canonical point to W. The class sequence is unchanged, so the standard
analysis applies unchanged. Which rotation is canonical is decided from
64-bit fixed-point approximations with an exact fallback, and Lemma 1 bounds their error.

**Distinguished points and resolution.** Walks run in trails. Each trail starts from fresh
coins and ends at its first distinguished point, which it records once together with its
coins. A key match is resolved by replaying the two trails (section 1.3, Lemmas R and D).
A fruitless cycle therefore cannot store the same point twice, and a degenerate match has
probability exactly 1/r.

Files:
- `accounting.py` recomputes every number with exact integers (output in
  `evidence/accounting.json`). It also checks the factorisation of r - 1 and the order of 2 (condition C-ord, section 2.3).
- `measurement/` holds the per-step instruction-count harness, the field code
  `bls_field.h`, the Python replica `walk_ref.py`, and the approximation check
  `frac_check.py`.
- `evidence/measurement/` holds the production totals of the current kernel on two input
  sets (`callgrind-doubling.json`), the 16-step production replay (`replay16-production.json`) and the checks of the
  restart and abandonment paths. `measurement/forced/` holds the forced-path worst cases of the rare operations
  (`summary.json`) and the common-path analysis (`common_path.json`).
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
G6 is the group of the sixth roots of unity in Z_r (r ≡ 1 mod 6), so |Gamma| = 6.

M(v) = v·2^384 mod p is blst's Montgomery representative, a fixed bijection of Fp onto
[0,p). Every comparison and bit test below is on M-values, read as integers.

Parameters:
- A = 2^16 adders.
- K = 2^30 walks: 2^16 processors × 2^14-walk batches.
- theta = 2^−60.
- L = 89663270682989036334838442605 ≈ 2^96.179.
- CAP = 64/theta. Trails have fewer than TRAIL_MAX = CAP + 64 steps.
- DP_MAX = 2·theta·S_max + 2K, where S_max = K·(L + TRAIL_MAX).
- TRAILS = DP_MAX + 2K, the number of trails started (section 2.4).
- E_DEPTH = TRAIL_MAX/2, the escape-depth bound of section 2.3. (It is the K of the round-8 model; K is already the number of walks here.)

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

**Fallback window and rate.** The exact computation runs exactly when (i) some g_j lies in the top window
{2^64 − 27, …, 2^64 − 1} (27 values), or (ii) two of the three g_j lie within 26 of each other (53 values
for a given partner). Under M1 (section 2.2) the three fractions are the fractions of x, beta·x and beta^2·x for a uniform
x. Each pairwise difference is then uniform, since (beta^i − beta^j)·x is uniform mod p for i ≠ j (beta is a
primitive cube root of unity). So the probability of the fallback per canon call is at most

    (3·27 + 3·53)/2^64 = 240/2^64 = 15·2^−60 ≈ 2^−56.09 < 2^−55.

The kernel never took the fallback in 2^19 measured steps (the forced fallback inputs of section 3.1 take it by construction). `measurement/frac_check.py` checks the floor gap
⌊f_j⌋ − g_j ∈ [0,12] on 200,000 random and edge-case inputs per seed; the largest gap seen
is 6.

### 1.3 Procedure

1. **Validate T.** Check it is on the curve and [r]T = O. If T = O, output 0.
2. **Adders.** For each adder i < 2^16:
   - draw c_i, d_i uniform in [0,r) (in the specification, from the seed of section 2.2, M2);
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
   - The replays give the common canonical point as [a1]P0 + [b1]T = [a2]P0 + [b2]T. If
     b1 ≠ b2, output x = (a2 − a1)/(b1 − b2) mod r after checking [x]P0 = T. If b1 = b2
     the run fails. By Lemma D this has probability 1/r.
   - Replay cost is bounded on every run: a replay of a trail with n live steps takes at most 65n steps (section 1.3, replay
     length), and the collision-replay charge of section 3 is paid in every run.

8. **Work meter.** Every run is charged at most W_max = S_max·2316.908726 + G
   instructions, with G the fixed global charges (initial starts, trail restarts, distinguished-point
   records, solve and verify, input validation, adder table).
   - Every walk step is charged C_COMMON = 2263 instructions, which bounds the common path
     including its batch inversion, and each rare event is charged its worst case.
   - The meter aborts the run when the total metered rare work exceeds 64 times its expectation bound (Lemma M').
     That happens with probability at most 1/64, counted in delta.
   - The run also fails if more than DP_MAX records arrive or more than K trails are abandoned.

**Escape.** A cycle's escape goes to 2·canon(M), where M is the cycle's least-key member. canon(M) is a point of M's class, so 2·canon(M) is a point of the class of 2M. The escape uses no table and no stored point: it is one point doubling, whose cost is charged at its forced-path maximum (Section 3.1). The escape point depends only on the cycle's members, as Lemma R requires. Its class data are fresh under M1 (section 2.2), and by D1 (section 2.3) its formal pair differs from every earlier point of the trail on the good event.

**Escape repeats.** The escape point of a cycle depends only on its least-key member, so it is the same every time the trail escapes from that cycle. A repeated escape point therefore means that the trail has re-entered a cycle it already escaped from. That is a within-trail coincidence, bounded in section 2.4 and 2.7. A counter in the kernel and the replica (`loops`) records each trail's last four escape keys and counts an escape whose key repeats one of them; it is a diagnostic only, it never changes the walk, and it is 0 in every check.

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

**Replay length (every run).** Replaying a trail of n < TRAIL_MAX live steps takes at most 65n replay steps, and fewer than
2^13·n unbatched operations, on every run:
- each replay step uses at most 70 operations: one evaluation, two for a 2-cycle check, 66 for an anchor check, and one escape;
- the escapes of the live trail split its steps into segments. Each segment runs from one escape point (or from the start) to the
  next escape point (or to the trail's end). By Lemma R the replay has the same escape points and the same first distinguished point,
  so it has the same segments and follows the same class sequence within each. A segment that ends in an escape of a cycle of
  length ℓ takes at least tail + ℓ live steps, because the closure of the cycle needs a full traversal of it, and at most tail + ℓ + 64
  replay steps, because the anchor check detects the cycle within 64 further steps (proof of Lemma R). A segment without an escape
  has the same length in both runs. The segments partition the live steps, so the replay takes at most n + 64E steps, where E is
  the number of live escapes. At most one escape occurs per step, since the 2-cycle and anchor checks are exclusive, so E ≤ n and
  the replay takes at most 65n steps. With at most 70 operations per step, 70·65 = 4550 ≤ 2^13.
This uses only the deterministic class map. It does not depend on the good event of section 2.4.

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

## 2. Success probability

Let T_steps = K·L be the main-phase steps, each walk's first L steps. Then
tau = T_steps^2/(2N) = 0.530300.

### 2.1 Heuristic H1 (collision statistics, unchanged)

**H1 (heuristic, disclosed).** The class map [W] → [canon(W) + R_h([W])], with uniform
random adders, behaves like a random mapping on the N classes, apart from the known
r-adding penalty. An r-adding walk loses about a factor (1 − 1/(m·A)) in the birthday
exponent, with A the number of adders and m = 6 the class size (Teske 2001;
Bernstein–Lange 2012; Bos–Kleinjung–Lenstra 2010 eq. (2) for m = 2).

Under H1, the probability that the first T_steps steps contain no useful class collision
is at most exp(−tau·(1 − eps)). Here eps = 4/A + 2·(3/A) + 2·(3·192/A^2) = 0.00015286. It
collects three margins:

- **Adder penalty, 4/A.** This is 24 times the 1/(6A) estimate.
- **Fruitless 2-cycles, 2·(3/A).** The 2-cycle entry rate is 1/(6A) per step (section 2.5), and each
  event yields at most 3 non-fresh states.
- **Short cycles, 2·(3·192/A^2).** The 3..64-cycle entry rate is at most 3/A^2 (section 2.5), and detection wastes at
  most 192 states per event (the 128 steps of section 2.5 are used in Lemma T' and the larger figure is kept here).

H1 is used only here, for the collision statistics, and in the DP overflow and exceptional-addition terms of section 2.9.

### 2.2 Premises M1 and M2

**M1 (class data, extended).** Every class-level value the walk reads is a function of the class alone. The values of
distinct classes are independent and uniform, independent of x. The values are:
- the class key, which fixes the adder index i = (key >> 1) mod A, the 64-bit filter bits and the distinguished-point test;
- the canonical rotation k and the canonical sign s;
- the coordinate fraction values f_j = 2^64·M(beta^j·x)/p of the class's canonical representative, which the exact-k test uses.

The last item is the extension of round 8. Under M1 the fallback probability per canon call is the probability that a
uniform fraction lands in the approximation window of Lemma 1, which is at most 15·2^−60 < 2^−55 (section 1.2).

**M2 (table).** The start coefficients (a_w, b_w) and the adder coefficients (c_i, d_i), with adder points R_i = c_i P0 + d_i T,
are fixed by the seed Python's `random.Random(1)`. The seed is a constant of the specification, and it is the value used by every
measurement (`measurement/reproduce.sh` defaults to it; the second measured input set uses `random.Random(2)`). The stream position
is fixed as follows. The first call is `rng.randrange(1, r)`, and its value is discarded: in the test generator it is the test target,
and the challenge target is not drawn from this stream. Then, for each of the 2^16 adders in order, then each of the 2^14 start points
in order, then each of the 1024 spare points, the specification draws a = `rng.randrange(r)` and then b = `rng.randrange(r)`, and sets
(c_i, d_i) = (a, b) for the adders and (a_w, b_w) = (a, b) for the start points.
Later trail starts (restarts beyond the 1024 spares) continue the same stream, two draws a, b per start, and are covered by M2 in the same way. The procedure then uses R_i = [c_i]P0 + [d_i]T with
the challenge target T. The coefficient values are a function of the seed alone, and they do not depend on the challenge target.
The test generator `measurement/walk_ref.py gen` draws in exactly this order, and the measured inputs are its output: `gen` with seeds 1
and 2 reproduces `in16t-1.txt` and `in16t-2.txt` byte for byte, which are the first 82944 lines of each seed's output, the two production
input sets of the callgrind evidence.
The table is fixed, so there is no randomness over it. Probabilities are over the uniform target x. M2 is the premise that this
fixed table behaves, for the walk and the collision statistics, like independent uniform coefficients independent of x. The algorithm's
own step 2 draws the adder coefficients uniformly, and the fixed seed is the specification's version of that draw. M2 is a second
premise, next to M1.

Probabilities below are over the uniform target x (the table is fixed by the seed, M2), and the class data follow M1. Nothing else is
assumed, apart from H1 in section 2.1.

**The good event G.** G is the event that no trail has a within-trail coincidence (section 2.4). Its complement has
probability at most DELTA_LOOP. The formal-structure results of sections 2.3 to 2.7 hold on G.

### 2.3 Formal structure (D1, D2) and the escape depth

Every walk point has an exact pair (a, b) in Z_r^2 with W = (a + b·x)P. Within one trail with start S, each walk point has
the formal form

    W = 2^m g S + sum_i n_i R_i,

with g in Gamma, integers n_i, and m the number of escapes so far in the trail. Steps add sigma(R_i) (which is ±lambda^j R_i)
and do not change the S-coefficient. An escape replaces M by 2·canon(M) = 2g'M, which doubles the S-coefficient. Its pair is
the same linear map applied to the pairs of S and of the R_i.

**Escape depth (E-depth).** On G no class is visited twice in a trail except inside a closed cycle, and each escape removes its
cycle. The cycles of distinct escapes are therefore disjoint sets of steps, and each has at least 2 steps. So a trail has at most
E_DEPTH = TRAIL_MAX/2 escapes.

**D1 (depth).** On G, if E is the (m+1)-th escape point of a trail, its S-coefficient 2^(m+1) g' differs from the S-coefficient
2^j g'' of every earlier point of the trail (j ≤ m).

*Proof.* Equality would give 2^k = g''/g' ∈ Gamma for k = m+1−j ∈ [1, E_DEPTH]. Gamma has exponent 6, so 2^(6k) = 1 and
ord_r(2) ≤ 6k ≤ 6·E_DEPTH. ∎

D1 holds under the check

    (C-ord)   ord_r(2) > |Gamma|·E_DEPTH = 3·TRAIL_MAX.

**C-ord is verified for the BLS12-381 G1 subgroup order r.** r − 1 = 2^32 · 3 · 11 · 19 · 10177 · 125527 · 859267 · 906349^2 ·
2508409 · 2529403 · 52437899 · 254760293^2 (the factorisation used by `accounting.py`). Its product is checked to equal r − 1,
each factor is prime by the Miller–Rabin test with the first twelve prime bases (deterministic below 3.18·10^23, and every factor
is below 2^29), and r is a probable prime under the same bases. The exact order is ord_r(2) = (r−1)/16, which is about 2^250.86.
It is certified by 2^(o) = 1 and 2^(o/p) ≠ 1 for every prime p dividing o = (r−1)/16. Since 3·TRAIL_MAX is about 2^67.6, C-ord holds
with a margin of about 2^183. No assumption is needed for this track.

**D2 (re-entry is not structural).** A re-entry into an escaped cycle is a later point (depth ≥ m+1) in the classes of an earlier
cycle (depth ≤ m). By D1 such a point is never formally equal to an earlier one. Its coincidence is a nondegenerate linear
equation in x, or a degenerate pair (section 2.4).

### 2.4 Degenerate and nondegenerate coincidences (D*, Lemma G) and DELTA_LOOP

Let e, e' be two points of one trail with formally distinct expressions (section 2.3). Their pairs are linear forms in the random
vector v = (a_S, b_S, c_1, d_1, …) with coefficients in Z_r, given the transcript. The difference map is nonzero, so under M2 and
for each g in Gamma,

    Pr[(a,b)_e = g·(a,b)_e'] ≤ 1/r      (Schwartz–Zippel for a nonzero linear map to Z_r^2).

Conditioning on the earlier history removes at most a fraction |Gamma|·t/r of the coefficient space (each earlier degenerate
equality is a hyperplane of density 1/r), so the conditional bound is at most 2|Gamma|/r. Summed over the pairs of a trail, the degenerate
pairs contribute at most 2|Gamma|·(TRAIL_MAX²/2)/r per trail.

**Lemma G (nondegenerate).** For the next point of a trail with t earlier distinct classes and no earlier coincidence,

    Pr[its class is among the t earlier classes | history] ≤ |Gamma|·t/(r − |Gamma|·t).

*Proof.* On the event that there is no coincidence so far, the transcript (indices, DP tests, canonical data) is a function of the
class data of distinct classes, which are independent of x by M1. Conditioned on it, x is uniform on the complement of at most |Gamma|·t roots of
the t earlier linear equations (each earlier point gives at most |Gamma| roots, one per element of Gamma). ∎

**Within-trail bound.** Summing over the points of every trail (at most TRAIL_MAX points each, TRAILS trails) and adding the adder-table
term (A^2 ordered pairs of adders, each degenerate with probability at most |Gamma|/r under M2):

    DELTA_LOOP := TRAILS·(TRAIL_MAX^2/2)·|Gamma|·(2/r + 1/(r − |Gamma|·TRAIL_MAX)) + A^2·|Gamma|/r.

DELTA_LOOP bounds the probability of any within-trail coincidence, degenerate or not. Its three parts are
- degenerate pairs: 1.04·10^-16,
- nondegenerate coincidences: 5.20·10^-17,
- adder-table pairs: 4.9·10^-67,

so DELTA_LOOP = 1.56·10^-16 (2^-52.6). The number of trails is TRAILS = DP_MAX + 2K = 1.67·10^20 (2^67.18): K initial starts, plus
DP_MAX restarts after records and K restarts after abandonments. The run stops when either restart count is exceeded, so these are
the trails the run can start. (The round-8 model writes DP_MAX + K for the trails the accounting restarts; the K initial starts are
added here, which makes the count conservative.)

By D2 every re-entry into an escaped cycle is one of these coincidences, so DELTA_LOOP bounds the re-entry probability (section 2.7).
The good event G is the complement of this coincidence event.

### 2.5 Closure rates (C', exact under M1)

A new point's successor index is (key >> 1) mod A, uniform over A (M1). Under M1 the closure rates are as follows.

- **2-cycle: rate exactly 1/(|Gamma|·A).** It needs an index repeat (1/A) and a Gamma cancellation (1/|Gamma| = 1/6).
- **3..64-cycles: rate at most c/A^2 with c = 3.** A simple cycle of ℓ distinct classes has independent index and unit data (M1).
  The cycle closes only if, for each index that occurs, the unit-weighted sum of its occurrences vanishes in Z_r. A singleton index
  cannot cancel. An index that occurs m ≥ 2 times fixes its last unit to one value, with probability at most 1/6 ≤ 1/2 given the others.
  So the probability is at most Σ_D N(ℓ,D)·2^(−D)·A^(D−ℓ), where N(ℓ,D) is the number of set partitions of ℓ items into D blocks of
  size ≥ 2 and D is the number of distinct indices. The terms are:
  - ℓ = 3: N = 1, D = 1, bound 1/(2A^2);
  - ℓ = 4: D = 1 gives A^(−3)/2, and D = 2 (two pairs) gives 3/(4A^2); total at most 0.7500076/A^2;
  - 5 ≤ ℓ ≤ 64: the exact sum of the bound is 1.55·10^-14 = 6.7·10^-5/A^2 (computed from the partition counts; the largest terms are ℓ = 5, 6).

  So the total is at most 1.2501/A^2 ≤ 3/A^2, and c = 3 is a valid (conservative) constant.
- **Structural cycles of length 65 or more.** The same count applies with ℓ ≥ 65. For 65 ≤ ℓ ≤ A/8 it gives (eℓ/(2A))^(ℓ/2) ≤ 2^−309,
  for ℓ > A/8 it gives at most 2^−4096 (plus a term below 2^−10000 for the case D < A/16), so the entry rate is below 2^−230 per step.
  Such a cycle is never escaped and traps its trail until the cap. This event is handled in section 2.6.
- **False filter matches.** A class visited for the first time has uniform key bits (M1), so a 63-bit filter matches it falsely with probability 2^−63 per filter.
  The two filters give 2·2^−63 per step, and each false match is rejected by the exact check.

Let p_c := 1/(|Gamma|·A) + 3/A^2 + 2·2^−63 = 2.54383·10^-6. This bounds the conditional probability that a step closes, and it is
a superset of the closure rates. A closure wastes at most 128 steps: a 2-cycle at most 3 states, and a cycle of length 3..64 is
detected within 64 + ℓ ≤ 128 steps. Escapes produce fresh points (on G, by D1 and M1).

### 2.6 Abandonment (Lemma T', first-moment bound)

Let n := TRAIL_MAX = 64·(CAP_PERIODS + 1) be the number of steps before a trail is abandoned, with CAP_PERIODS = 2^60. On G,
each step that is not a closure produces a class not yet visited in the trail, whose distinguished-point test has probability θ = 2^−60
independently of the past (M1: the test is a function of that class's key). Let C_n be the number of closure steps among the first n steps.

**First moment.** For k ≥ 1, by the chain rule, the probability that k given steps all close is at most p_c^k (each conditional
bound holds given the whole past). Hence

    Pr[C_n ≥ k] ≤ binom(n,k)·p_c^k ≤ (e·n·p_c/k)^k.

**Fresh points.** If C_n < k, at least n − 128k steps are fresh, since each closure wastes at most 128 steps, and none of them is a
distinguished point with probability at most (1−θ)^(n−128k) ≤ exp(−θ·(n−128k)).

**Choice of k.** Take c = A = 2^16 and k = ⌊n/c⌋ = 2^50. Then n/k = c·(1 + 2^−60), so e·n·p_c/k = e·c·p_c·(1+2^−60) = 0.45317 ≤ 1/2, and
Pr[C_n ≥ k] ≤ 2^−k, which is 0 for every practical purpose. The fresh count satisfies θ·(n − 128k) = 63.875 + 2^−54 (since 128/c = 2^−9). Hence

    Pr[abandon] ≤ 2^−k + exp(−63.875) + n·2^−230 ≤ 1.8174·10^−28 + 2^−164 ≤ 1.8174·10^−28.

The third term is the structural trap of section 2.5: a structural cycle of length at least 65 entered at rate below 2^−230 per step
traps the trail for at most n steps. The per-trail bound exp(−63.875) is sharper than the exp(−63.5) of the round-8 model, which
evaluates to 2.64·10^−28 (the round-8 model's "2.4·10^−28" is not that value); the sharper bound is used because 128/c = 2^−9.

**Expectation and failure.** Expected abandoned trails are at most TRAILS·Pr[abandon] = 3.03·10^−8. The run fails by abandonment with
probability at most that expectation divided by K (Markov), so

    DELTA_ABANDON = TRAILS·Pr[abandon]/K = 2.83·10^−17,

counted in delta.

### 2.7 Re-entry into an escaped cycle (Lemma L)

A re-entry into an escaped cycle is a later point of the classes of an earlier cycle (D2). It is therefore a within-trail
coincidence, and by section 2.4 its probability over the run is at most DELTA_LOOP = 1.56·10^−16. A repeated escape point (the
diagnostic `loops`) is such a coincidence. Outside G the walk is not restarted: the trail still ends at its first distinguished point
or at the cap, and its work is metered (section 2.8).

### 2.8 Expected rare work (Lemma M', every processed step counts, revisits included)

Per walk step the metered rare work has expectation at most E_RARE = sum over event types of rate × worst case:

| event | rate per walk step | worst case (instructions) |
| --- | --- | ---: |
| 2-cycle check (exact check; true or false predecessor match) | 1/(6A) + 2^−63 | 328,232 |
| anchor closure (3..64), with all its look-ahead steps (at most 64) | 3/A^2 + 2^−63 | 10,861,544 |
| exact-k fallback (secp, BLS only; 2^−55 per canon call, ≤ 2 canon calls per step in expectation) | 2·2^−55 | 682 |
| within-trail coincidence (section 2.4; degenerate or not) | 6·TRAIL_MAX·(2/r + 1/(r − 6·TRAIL_MAX)) = 2.53·10^−56 | TRAIL_MAX·11,194,085 |

The rows are:
- *2-cycle and anchor checks.* Their rates are the closure rates of section 2.5 plus the false-match rate of one filter. Their worst cases
  are the forced-path maxima of section 3.1.
- *Exact-k fallback.* The rate is per canon call (section 1.2, Lemma 1 and M1). The batch step calls canon_k once per walk step, in
  `batch_step` (`measurement/step.c`, line 611). The rare paths add at most one canon call each, in `rare_2cycle` (line 548) and
  `rare_anchor` (line 574), and each restart calls `select_adder` in `new_trail` (line 513), which is inside the restart charge START.
  The canon calls of `f_single` (line 452) are inside the worst cases of the two checks. So the expected number of canon calls per walk
  step is 1 + 2.6·10^−6 + … ≤ 2, and the row is 2·2^−55. Its worst case is F_EXACT = 682 instructions (section 3.1).
- *Within-trail coincidence.* A coincidence in a trail of length up to TRAIL_MAX can leave the walk in a class cycle longer than 64.
  There each anchor check costs W_ENTRY every ℓ steps, so the remaining work of the trail is not covered by the per-step reserve. The
  worst case is therefore TRAIL_MAX·W_STEP, with W_STEP = C_COMMON + W_CHECK + W_ENTRY + 3·F_EXACT = 11,194,085 instructions (one 2-cycle
  check, one anchor check and three canon calls outside f_single in every step). The bound is 8.3·10^26 instructions, and the row adds
  2.1·10^−29 instructions per walk step.

E_RARE = 0.84232384277012 (the rows above sum to this value; the first two rows give 0.842323842770081, and the fallback row adds 3.8·10^−14).
The reserve per walk step is 64·E_RARE = 53.908726 instructions. Markov gives Pr[metered rare work > 64·S_max·E_RARE] ≤ 1/64, counted in delta.
Every escape is one of these events: an escape happens only after a 2-cycle check or an anchor check, so its doubling is inside that
event's worst case. The expected escape work is covered by E_RARE, and no escape is uncharged.

### 2.9 Detection, resolution and the success statement

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
- **Within-trail coincidence (degenerate or nondegenerate, including re-entry):** DELTA_LOOP = 1.56·10^−16 (section 2.4, 2.7).
- **Abandoned trail:** DELTA_ABANDON = 2.83·10^−17 (section 2.6).
- **Work meter abort (step 8):** at most 1/64 (section 2.8, Markov).
- **Degenerate first match:** exactly 1/r (Lemma D).
- **DP_MAX overflow:** the expected number of records is at most DP_MAX/2, so by Chernoff under
  H1 the overflow probability is at most exp(−DP_MAX/12), which is below 2^−60.
- **Exceptional additions:** a batch whose product is zero, or a point that a formula cannot
  handle, aborts the run. This has probability at most 2^−110 under H1.

The degenerate-match term 1/r and the two H1 terms (DP overflow and exceptional additions) total below 2^−20 and are counted in
delta = 2^−20. The other three are counted separately. Therefore

    δ = 2^−20 + DELTA_ABANDON + DELTA_LOOP + 1/64.

**Pr[success] ≥ 1 − exp(−0.5303·(1 − 0.00015286)) − δ = 0.411524 − 2^−20 − 2.83·10^−17 − 1.56·10^−16 − 1/64 = 0.395898 ≥ 0.39.**

The exact value is 0.39589793195. The budget is 1.0299·√N steps, against 0.9943·√N for success exactly 0.39. The 0.0215
of slack in the model bound covers the sampling error of the toy evidence below. The
bound is exact arithmetic on H1, M1, M2 and the stated counts, not an empirical estimate.

### 2.10 Evidence for H1, M1, M2 and the procedure

`experiment.py` runs the procedure of
section 1.3 end to end on toy curves y^2 = x^3 + 4 over
p ≡ 1 mod 3, in prime-order subgroups since that curve always has 3-torsion (`"b": 4`). It uses:
- the same class key, canonical rule, adder images and index (keys use Montgomery
  representatives, any fixed bijection being equivalent);
- both cycle rules and their escapes (to 2·canon(M), as in the kernel: `"escape": "double"`);
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
Offline runs (`evidence/trails/`, pooled in `summary.json`): 4,600 trials on 4 curves with
N = 2^29.2, 2^28.8, 2^27.5, 2^30.0 classes, adder tables of 2^11 and 2^16,
16 slots and theta = 2^−5.
- Success of the whole procedure at tau = 0.5303: 1,935 of 4,600, 0.4207, Wilson 95% interval
  [0.406, 0.435]. The model 1 − exp(−tau) is 0.412.
- Every success was verified ([x]P = Q). Degenerate first matches: 0.
- Production 64-bit filters. Every recorded trail was replayed with a different anchor
  schedule: 2,990,805 replays, 0 that missed the recorded point.
- Lemma R stress (`evidence/trails-stress.json`): the same procedure with the filters
  cut to 8 bits over 2,000 trials, so false filter matches were frequent.
  - 1,326,881 trails were replayed with a different anchor schedule, and 0 missed the
    recorded point.
  - 94,878 distinguished points were first reached during a look-ahead.
  - Those runs' success rates are not used: each such look-ahead adds up to 65 points
    that the step budget does not count.
- Replay-order sensitivity (`evidence/lemma-r-sensitivity.json`): 200 trials per curve with 6-bit filters and 16
  adders, every recorded trail replayed. With the current order, 0 of 59,977 replays (j0) and 0 of 95,625 (a3)
  missed the recorded point. With the earlier order (an offline switch), 8 of 59,351 (j0) and 47 of 92,993 (a3)
  did, and that order also abandoned 28 and 89 trails.
- Abandoned trails: 2. Longer fruitless cycles escaped: 4. The finishing phase added at most
  2.2% to the K·L steps.
- 2-cycle entry rate per step:
  - 7.61e-05, 7.69e-05, 7.55e-05 with A = 2^11 (1/(6A) = 8.14e-05);
  - 2.17e-06 with A = 2^16 (1/(6A) = 2.54e-06).
<!-- /toy-evidence -->

The declared initial input (`bits` 28, `b` 7, 80 trials, seed 2801, production 64-bit
filters, every trail replayed, `escape` double; `evidence/declared-run-double.json`) runs in a few seconds. It uses prime-order y^2 = x^3 + 7
curves, which exercise the same order-6 class code. For b = 4, fields of at most 28 bits
give subgroups of only 2^19 to 2^23 classes, where each slot's final trail adds 20 to 40%
to K·L, so success runs high; the offline b = 4 runs above use larger fields. It reports
35/80 = 0.4375, Wilson 95% interval [0.334, 0.547], against the model 0.412, with 0
degenerate matches, 0 distinguished points reached during a look-ahead and 0 replays, of 15,059, that missed their recorded point. With 80 trials the standard error is about 0.05,
so it is a reproducibility check; the success claim rests on the pooled offline trials
above.

These toy runs test H1 and the procedure at toy scale. They do not test M1 or M2 at production theta, and they cannot
test C-ord, which is checked arithmetically (section 2.3).

## 3. Computation

**No reference operation.** Every instruction of a walk step is charged at 1/C with C = 2000, as
COST_MODEL.md's H + R/C with H = 0. No reference operation is charged.

**Charge per walk step.** Two parts, charged separately:
- The common batched path is charged C_COMMON = 2263 instructions per walk step: the larger of the plain and anchor steps, each measured as the whole-run instruction total difference between runs with one more batch step
  on production inputs, per walk. On the measured inputs the plain step is 2248.94 and the anchor step 2258.94 (the anchor batch of
  the second input set contains no rare event, Section 3.2). The one data-dependent branch of the common path, the canon_pick
  selection (k_max = 2), adds at most 2 instructions per walk step, which gives 2260.94; C_COMMON = 2263 bounds that.
- The rare events, each charged its forced-path worst case (Section 3.1): a 2-cycle check at most
  328,232 instructions and an anchor check, with its look-ahead steps, at most 10,861,544; the exact-k fallback at most 682 per canon call;
  and the within-trail coincidence at TRAIL_MAX·W_STEP. Their expected
  counts are covered by the Markov reserve 64·E_RARE = 53.908726 per walk step (Lemma M').

The charge per walk step is therefore 2316.908726 instructions.

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
    (16 batch steps) on the production input set, with 0 mismatches (`evidence/measurement/replay16-production.json`).
- blst is still used for the one inversion per batch (`blst_fp_inverse`) and for input
  conversion. Those instructions are counted like any others.

Inline assembly is ordinary scalar AArch64 code. Callgrind counts each executed
instruction exactly as for compiled C, which is how the organizer measured blst's assembly.

**3.1 Worst cases of the rare-path operations.** Each rare event is charged a forced-path maximum, not a sampled
average. The primitives are f_single (adder selection and the unbatched addition), add_single, escape (its least-key
search, its canonical point and its doubling), same_class, point_is_dp, and the exact-k fallback of select_adder. For each primitive
the data-dependent branches were enumerated from the source and the disassembly. Each conditional correction in a
field routine that is not branch-free was replaced, for the measurement only, by a branch-free routine that costs at
least as much on every path (a bounding build; equality with the original is checked on random and edge inputs).
The remaining control variants (the equal-x and equal-y tests of the addition, the parity of the least member, the
rotation fallback where present, the length of a cycle, and the memcmp outcomes of the library comparison) were forced
by chosen inputs, and the maximum over the variants was taken. Each bounded count was checked to be constant per call and unchanged when the number of calls doubles (1000 and 2000 calls per variant). Measured maxima (Callgrind Ir per call, organizer image and flags,
`measurement/forced/summary.json`):

| primitive | worst case (Ir) |
| --- | ---: |
| f_single (adder selection and addition) | 163,118 |
| add_single | 162,442 |
| same_class | 211 |
| point_is_dp | 99 |
| escape of a 2-cycle (canon, then doubling) | 163,780 |
| escape of a 65-cycle (canon, then doubling) | 172,474 |
| select_adder on the exact-k fallback path (rotation test, canon_k_exact and its call, table index) | 682 |

The exact-k row is the forced maximum of the fallback variant `fb` of the primitive `sel` in the forced harness
(`bls-ir/forced`, the same step.c as `measurement/step.c`, sha256 cc984345…). On the bounding build the variant gives 698 Ir per call,
including the 16 Ir loop overhead, which gives 682 net. It is constant over 2000 fresh fallback inputs, and the same
measurement was repeated for this round in the organizer image (`drive.py`, Callgrind per call; `measurement/forced/fallback_select.json`). The canon_k_exact primitive
alone, with its call, is 594 − 16 = 578 Ir. The fast-rotation variant of select_adder is 107 Ir, so the fallback adds at most 575
Ir over the common path. The charge of 682 per fallback is the whole fallback path, which is conservative.

The repeat counter (a diagnostic, Section 1) costs at most 36 instructions per escape (measured with its ring full, the dearest state) and is inside the overhead.

A 2-cycle check costs at most f_single + same_class + escape(2) + point_is_dp plus a control overhead of at most
1,024 instructions, and an anchor check costs at most f_single + escape(65) plus, for each of its at most 64
look-ahead steps, f_single + point_is_dp + same_class, with the same overhead per step. The overhead bound is far
above the measured control code of these functions, and the sums are above the longest path through the code.
The forced-path harness builds each primitive on its own, and in the measured build add_single and canon_k_exact
are 3 and 4 instructions cheaper than in the shipped code. An event calls each primitive a few times, so these
differences total well under the 1,024-instruction overhead of that event (or of each look-ahead step), and the
charged worst cases remain upper bounds for the shipped code.
Scalar multiplications and replays use f_single's worst case as their unit, an upper bound for an addition or a
doubling.

The canon call of a batch step is not inside f_single. It is therefore charged by the fallback row of section 2.8, not by the
worst case of f_single. Each canon call of a rare path that is not inside f_single is charged the same way.

The unbatched inversion is blst_fp_inverse. It computes a constant-time inverse, checks it with one multiplication, and falls back to a Fermat chain of 133,370 instructions if the check fails. No tested input reached the fallback, but no proof excludes it, so the worst cases above include it on every inversion.

**Measurement (H3).** `measurement/step.c` runs the trails, steps, cycle handling, distinguished
points, restarts and trail cap. Its common batched path was measured on production inputs as
whole-run instruction totals (Callgrind `Ir`, GCC 14.2.0 in the organizer's image,
`-O3 -fno-tree-vectorize -mgeneral-regs-only`), with the rare paths out of line so that each
path's cost is separable. The first batch step has no anchor and no rare event. The anchor step is 10.0 instructions per walk step above the plain step
when its batch contains no rare event, and a rare event adds its own cost to it. H3 is a measurement premise, not a probabilistic heuristic.

- Plain step (steps 1 and 2), per walk step: 2248.94 instructions (measured, input set 1).
- Anchor step (steps 63 and 64), per walk step: 2258.94 instructions (measured, input set 2, whose anchor batch is event-free).
- C_COMMON = 2263 per walk step: the measured anchor step plus 2 for the canon_pick branch (2260.94), rounded up. The bound for
  walks that are all long is 2259.61, so 2263 is conservative.

**3.2 The common path on two input sets.** The measurement was repeated on a second production input (seed 2 of the same generator).
Callgrind totals at 1, 2, 63 and 64 batch steps are 36,846,673/73,693,333/2,321,461,362/2,358,530,610 on input 1 and 36,846,721/73,693,575/2,321,520,743/2,358,531,190 on input 2. The plain batch is 36,846,660 (input 1) and 36,846,854 (input 2).
The difference, 194, is exactly 2 × 97: canon_pick is the only data-dependent branch of the common path, it skips 2 instructions
in walks whose rotation is k = 2, and 5,485 of the 16,383 walks take the short path on input 1 and 5,388 on input 2. No rare event occurs in a
plain batch of either input. The anchor batch of input 1 contains one two-cycle escape (2 events by step 63, 3 by step 64), and
that of input 2 none (3 and 3); its event-free anchor batch is 37,010,447. `measurement/forced/common_path.json` lists every conditional
branch of the loops of batch_step with its counts: the filter-match, abandonment, exact-k and candidate-distinguished-point
branches are never taken in an event-free batch, the loop tests do not depend on the data, and canon_pick is the only branch whose
two outcomes both stay on the common path. If every walk took the long path, the plain batch would be 36,857,630 on both inputs (2249.61
per walk step) and the anchor batch 37,021,471 (2259.61); adding the 2-instruction maximum to every walk gives 2250.94 and 2260.94.
C_COMMON = 2263 exceeds all of these.

**Rare-path checks.** Distinguished points and false 64-bit filter matches do not occur in a
production-size sample, so small builds of the same source and replica exercise them
(`evidence/measurement/trail-paths.json`, with 2^6 adders and 2^12 walks, 2^4 adders in the last build):
- DP_BITS = 4: restarts at distinguished points, with 2-cycles;
- DP_BITS = 12 with a one-period cap: abandonments;
- DP_BITS = 6 with the filters cut to 8 bits: 2-cycles, longer cycles, and distinguished points first
  reached during a look-ahead;
- 2^4 adders and DP_BITS = 8 over 800 steps: many 2-cycles and longer cycles with a tiny adder table; the repeat counter stays 0.

In every build the replica matched every final point, every event count (the repeat counter included), and a digest of every distinguished record, look-ahead records included (`evidence/measurement/trail-paths.json`).

**Total.** S_max = K·(L + TRAIL_MAX) ≤ 2^126.178504 walk steps.

| Component | Bound (units) | log2 |
| --- | ---: | ---: |
| walk steps | 2316.908726·S_max / C | 126.390705 |
| adder table | 2^16 adders × (2 scalar mults × 768 unbatched ops × 2^15 + 2^14) / C | 32.93 |
| initial walk starts | K × (2 scalar mults × 768 × 2^15 + 2^14) / C | 46.93 |
| trail restarts | (DP_MAX + K) × (2 scalar mults × 768 × 2^15 + 2^14) / C | 84.11 |
| distinguished points | DP_MAX = 2^67.18 entries × 2^20 instructions / C | 76.21 |
| collision replay | 2 trails × TRAIL_MAX live steps × 2^13 unbatched ops (replay ≤ 65n steps × 70 operations ≤ 2^13 n) × (2^15 + 2^10) / C | 86.36 |
| solve, verify, output | 2 scalar mults + 2^24 instructions | 17.03 |
| input validation | 2 scalar mults + 2^24 instructions (includes the [r]T = O check) | 17.03 |

Rounded up, **time_log2 = 126.390705**. Preprocessing (adders and initial walk starts) is included. All
processors' work is summed, and elapsed time plays no part. The exact value before rounding is 126.3907049948 (126.39070499484512 in `evidence/accounting.json`). The collision replay is
charged on every run at 2^13 unbatched operations per live step (section 1.3).

The restart term (DP_MAX + K) counts the restarts; the canon call of each restart is inside START. The K initial starts are in the
initial-walk-start term.

## 4. Memory and advice

Peak memory is at most the sum of:
- 2·DP_MAX·128 bytes for the distinguished-point table (key, the trail's coins, hashing
  overhead);
- K·640 bytes of walk state (points, keys, coins, trail counter, batch buffers);
- 2^16·6·96 bytes of adder images;
- 16 MiB.

That is 2^75.178504 bytes. Nonuniform advice is 0 bytes.

## 5. Limitations

- **Status.** The argument is not established at production size. M1 is the heuristic that carries the claim: class data
  (and the canonical coordinate fractions) are uniform and independent per class. Its consequences in this proof are C', T', M'
  and the fallback rate; the toy runs (theta = 2^−5) cannot test it at production theta.
- H1 is the standard but unproved random-mapping heuristic for additive class walks. It
  is supported here by the literature and by toy-curve experiments only. The model bound
  is 0.4115, 0.0215 above 0.39.
- **M2** (the table is fixed by random.Random(1) with the stream position of section 2.2) is a second premise. There is no randomness over
  the table: the probabilities are over the uniform target x, and M2 is the premise that this fixed table behaves like independent uniform
  coefficients independent of x. The seed is not a probabilistic choice; the premise is what makes the fixed table typical.
- **Off the good event G** (probability at most DELTA_LOOP = 1.56·10^−16) the success argument of section 2 does not apply, and that event
  is counted in delta. The cost bound does not depend on G: the replay length (at most 65n steps, section 1.3), the walk steps, the restart
  counts and the meter's rare-work abort hold on every run.
- **Long trapped cycles.** A within-trail coincidence can trap a trail in a class cycle longer than 64, where the anchor check fails and no
  escape occurs. This is inside DELTA_LOOP (it is off G) and inside the revisit row of section 2.8, whose worst case is a crude bound.
- Instruction counts come from one implementation and compiler. The common path was measured on two production input sets (Section 3.2), and its worst case is bounded there. The rare operations are charged forced-path maxima (Section 3.1); the bounding builds are checked for dominance on the forced paths, not proved for every input.
- The field multiplication is the package's own assembly, tested against Python big
  integers but not formally verified.
- The variable-time paths of the rare operations are bounded by forced-path maxima (Section 3.1). The exact-k path is charged separately as F_EXACT = 682 per canon call, and the forced maximum is measured, not proved for every input.
- The fallback rate 2^−55 depends on M1's extension to the coordinate fractions. If that extension failed, the fallback would have to be charged at its worst case for every canon call, which costs about 2.4 bits on time_log2. That is not applied.
- About 2^67 distinguished points and 2^126.4 work are needed. This is not a practical
  attack and does not weaken BLS12-381's security assumptions.
- The pairing (MOV) route to GT is not used.
- The toy experiment runs the whole procedure, including resolution, but only at toy
  scale: theta = 2^−5, 16 slots, and N far below production. Its trail starts and its
  finishing phase are a few percent of the steps there, and negligible in production.
- Distinguished points and false 64-bit filter matches never occur in a
  production-parameter measurement sample. The kernel's restart, abandonment and
  false-match paths were checked on small builds of the same source.
- The work meter's success loss, 1/64, is a Markov bound on the expected metered work under M1 (Lemma M').

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
echo '{"bits":28,"trials":80,"seed":2801,"b":7,"adder_bits":11,"dp_bits":5,"walks":16,"tau_e4":5303,"filter_bits":64,"check_all_replays":1,"escape":"double"}' | python3 -I submissions/bls12-381-g1-ecdlp-plain/rho-sigma6-asm/experiment.py
# offline evidence runs used larger inputs, e.g. {"bits":40,"trials":700,"seed":4001,"escape":"double"}
```

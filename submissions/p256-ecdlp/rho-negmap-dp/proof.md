# p256-ecdlp: parallel Pollard rho with distinguished points and the negation map

Analytic cost bound; no execution at production size is claimed. All counts are scalar AArch64 instructions divided by C = 1000, as the track prescribes. The
unit C is a normalization only: every P-256 operation below is charged in instructions.

## 1. Problem and notation

Input Q = [x]P, x uniform in [0, n), n the subgroup order, log2 n = 255.999999999664. The
negation map identifies R with -R. The number of classes is N' = (n+1)/2 (the identity is a
singleton class), log2 N' = 254.999999999664.

The identity target Q = infinity is the profile's included x = 0 case: the algorithm returns 0
at once, [0]P = infinity = Q passes the stated check, and no walk runs. Everything below assumes
Q is not the identity, so Q has a unique x in [1, n).

## 2. Algorithm

Walk function. A table of R = 2^24 entries T_j = c_j P + d_j Q is built once, with c_j, d_j
uniform nonzero 64-bit integers and T_j stored in affine form with its c_j, d_j. For a class
representative X with affine coordinates (u, v) in Montgomery form, define
j(X) = bits 232..255 of u, the top 24 bits (identical for X and -X), and
f(X) = canon(X + T_j(X)), where canon picks of {Y, -Y} the one whose Montgomery y-limbs are the
smaller as a 256-bit integer. f is a deterministic function on classes (an r-adding walk, Teske
2001, with r = 2^24 >= 16). The three bit fields used by the walk are disjoint: the
distinguished-point test reads bits 0..95 of u, the escape index bits 224..231, the table index
bits 232..255, so no step is forced to a particular index by being distinguished or by being an
escape.

Walkers. m = 2^14 walkers start at random points a0 P + b0 Q (a0, b0 uniform 64-bit integers,
the same size as the table coefficients) and run in lockstep. The 2^14 starting pairs are
distinct points with overwhelming probability (a collision of two starts needs a0 + b0 x equal
mod n for two of the pairs, probability below 2^-220), and the first table step, with its
2^24-way choice, spreads them. Every walker keeps a canonical point (u, v), a sign s = +-1,
signed integers (A, B) with the invariant point = s (A P + B Q), a step counter, and a
checkpoint of four limbs (32 bytes). A step adds T_j: A += s c_j, B += s d_j (signed 4-limb
integers, never reduced mod n during the walk), and if canon negates the sum, s := -s. A walker
never restarts. Size of (A, B): the start values are below 2^64 and every step, normal or
escape, adds a coefficient below 2^64 in absolute value; a walker takes at most K/m = 2^113.5648
steps (K below counts every step, lap and escape steps included), so |A|, |B| < 2^64 (1 +
2^113.5648) < 2^178. The 4-limb signed (two's complement) integer holds values up to 2^255, so
no overflow and no reduction occurs during the walk; (A, B) is reduced mod n only when a DP is
sent.

Distinguished points. A point is distinguished (DP) when the low 64 bits of u and the low 32
bits of the next limb are zero: theta = 2^-96. A walker reaching a DP reduces (A, B) mod n and
sends (u, A, B) (96 bytes) to a store keyed by u, then keeps walking. Two records with equal u
and different (A, B) give A1 P + B1 Q = e (A2 P + B2 Q), e = +-1, so
x = (A2 e - A1) / (B1 - e B2) mod n, tried for both e. Each candidate x is verified by
checking [x]P = Q. B1 = +-B2 has probability about 2^-250.

Batched affine step (the reference operation). The round computes d_i = u_T - u_i for all
walkers, a Montgomery batch inverse (3 field multiplications per walker: one prefix product,
two to unwind) and one inversion per round, then lambda = (v_T - v_i)/d_i, u' = lambda^2 - u_i - u_T,
v' = lambda (u_i - u') - v_i: 3 multiplications (S = M). That is the 6M affine addition of the
reference unit. A zero denominator (the walker equals +-T_j, probability about 2^-231 per step)
is tested for and sent to a slow path (affine doubling or identity), bounded below.

Fruitless cycles (Bos-Kaihara-Kleinjung-Lenstra-Montgomery 2010; Bernstein-Lange-Schwabe 2011).
With the negation map a walk falls into a 2-cycle of f with probability 1/(2R) = 2^-25 per step.
A cycle of 2l classes needs every table index used in it to occur an even number of times, so
cycles have even length and, with e = 2l - (number of distinct indices) >= l index coincidences
among the 2l steps, occur with probability at most 2 (2l/R)^l per step for l >= 2 (the number of
ways to choose which earlier step each repeated index copies is at most (2l)^e, and each
coincidence has probability 1/R). Such a cycle is a cycle of the function f on classes, so a
walker that enters it can never leave by f.

Detection rule (fixed period; it only decides when to react, never where to go). Each walker has
a step counter c and a checkpoint of the full canonical u, all four limbs (32 bytes). When c
reaches 64 the checkpoint is set to the four limbs of u of the current point and c := 0; every
new point is compared with the checkpoint on all 256 bits, the compare exiting on the first
limb that differs. The comparison is a full-width equality test, not a limb test with a
provisional match: the walker reacts only when all four limbs are equal. A match then proves that
the walk has returned to a class it held at the checkpoint; because f is deterministic on
classes, the walker is on a cycle of f, and the walker's own schedule has no further effect on
its trail, as shown next. A cycle of length L <= 64 that a walker has entered is detected within
128 steps of entry: a checkpoint is set within the first 64 steps, and the walker returns to it
within L <= 64 further steps, before the next refresh.

Spurious triggers. A match on all 256 bits is the return to one exact class, so a spurious match
would need a random point to reproduce a given 256-bit u, probability 2^-256 per step: that is
2^-160 times the distinguished-point rate 2^-96, and over the whole run of K = 2^127.5648 steps
the expected number of such coincidences is below 2^-128. There is therefore no false-trigger
budget to spend, and the earlier 64-bit checkpoint, which matched at 2^-64 per step (2^32 times
the DP rate) and treated the match as a cycle, is gone. What the check costs is charged: the
four checkpoint loads, the four limb compares and their branches are the cycle line of Section 4
(20 instructions charged against 15.05 itemised), and the checkpoint store of two `stp` plus one
counter reset happens once per 64 steps, inside the same line.

Escape as a function of the cycle. On detection the walker is at a cycle member X0 (a canonical
point) and walks the cycle once more in the ordinary lockstep batch (L further steps, L <= 64,
ending again at X0), additionally keeping the member with the smallest u (compared as a 256-bit
integer in the stored Montgomery representation) together with its full state (u, v, A, B, s).
Distinct classes of a cycle have distinct u, so the minimum m* is unique and depends only on the
set of the cycle's classes, not on X0, on the walker, on its counter or on its checkpoint.
Let e(m*) be bits 224..231 of u(m*). The walker restores the snapshot of m* and takes one
stand-alone affine addition (its own inversion) of the escape point T'_{e(m*)} from a second
table of 256 points, built like the main table with its own nonzero 64-bit (c', d') drawn
independently, updating (A, B, s) exactly as a normal step. It then sets the checkpoint to the
escape point and c := 0, so the next compare cannot match a class from before the escape, and
continues. The result E(C) = canon(m* + T'_{e(m*)}) is a function of the cycle C only. Two
walkers that merged earlier are in the same cycle, possibly at different members, after
different numbers of laps and with different counters; each one detects at its own time, and both
leave from the same point m* to the same point E(C). They stay merged. Adding T' (64-bit
coefficients) rather than doubling keeps |A|, |B| bounded without reduction, and because T' is a
separate table with independent coefficients the escape point has no algebraic link to the table
steps that closed the cycle.

Re-entry into the same cycle. The rule is deterministic per cycle, so if E(C) lay in the basin of
C every walker caught in C would loop m* -> E(C) -> ... -> C until its step cap, a trap. Whether
it does is a question about the point E(C), not about merging. A walker started at E(C) re-enters
C only if E(C) lies in the basin of C (the classes whose f-orbit reaches C) at a distance of at
most T = K/m = 2^113.5648 steps, since that is all the steps a walker has left. In the
random-mapping model each class has a Poisson-like in-degree of mean 1, so on average there is
one class at each distance from a given cycle member and the basin within distance T holds at
most L(T+1) <= 64 x 2^113.5648 = 2^119.565 classes in expectation. The class of E(C) is a
uniform-like point relative to that set, so

    Pr[E(C) re-enters C] <= 2^119.565 / N' = 2^119.565 / 2^255 = 2^-135.435 per cycle.

The escape is deterministic per cycle and all escapes from one cycle share the single point
E(C), so the union bound runs over distinct cycles, at most the number of cycle events a run can
have, K x 2^-25 = 2^102.564 in the whole run. The expected number of trapped cycles is therefore
at most 2^102.564 x 2^-135.435 = 2^-32.871, and a trapped cycle can hold no more than all m
walkers, each wasting at most its remaining T steps, so at most K steps: the expected waste is at
most 2^-32.87 of all steps. This is a model bound (random mapping, and E(C) independent of the
cycle's basin), listed as a heuristic in Section 6; the toy experiment of Section 7 shows the
mechanism directly at a basin fraction of 2^-4.5, where the trap repeats until the cap.

Cost of one cycle event, in step equivalents of 1193 instructions: at most 128 steps of
detection delay (steps already taken inside the cycle), at most 64 steps of the minimum-tracking
lap, charged 40 more instructions each for the 4-limb compare and the rare snapshot copy
(64 x 1233 / 1193 = 66.2 < 67), and the escape addition at most 100 units = 100000 / 1.193 /
1000 = 83.8 step equivalents; in total at most 128 + 67 + 84 = 280 step equivalents. With the
per-step event rate lambda <= 2^-25 + sum over l >= 2 of 2 (2l/R)^l <= 2^-25 + 2^-42 (the l = 2
term is 2^-43), the amortised loss is at most 280 x 2^-24.9999 = 2^-16.87 of the steps, about
1.7 x 10^-5.

Long cycles. A cycle of more than 64 classes is never detected: its length exceeds the checkpoint
period, so no walker returns to a checkpointed class in time, and the walkers in it are stuck
until the cap. The two mechanisms by which such a cycle closes are (i) the random-mapping cycle
of f, which for a walk of T = 2^113.5648 steps and N' = 2^255 classes has probability about
T^2/(2N') = 2^-28 for any repeat at all, and (ii) a zero sum of the signed table coefficients
with no repeated index, a closure that the 64-bit coefficients allow at about 2^-128/L per step
for a cycle of length L. Summing (ii) over L from 65 to the longest cycle a walker can traverse,
L <= T, gives a closure rate of 2^-128 x (sum of 1/L over 65 <= L <= 2^113.5648) = 2^-128 x 74.54
= 2^-121.78 per step; a walker that enters such a cycle wastes on average half its remaining
steps, so the share of all steps lost is 2^-121.78 x 2^112.565 = 2^-9.22. The three cycle terms
therefore total at most 2^-16.87 + 2^-32.87 + 2^-9.22 < 2^-9.19, and Section 3 budgets 2^-9 for
all of them. The budget is a factor 1.15 above the sum of the model values, and the claim does not
depend on it: Section 3 shows the score stays below the claimed 127.82 for any cycle budget up to
2^-8.72.

## 3. Success probability

Walk quality. Walks are modelled as a random mapping on N' classes (Teske 2001; r = 2^24, with
the negation map and the cycle escapes above). The model is not taken at its ideal value: a
walk-quality factor rho <= 1 multiplies the collision rate, rho = 1 being an ideal random
mapping. Two corrections are known for finite r and are absorbed into rho here: the r-adding
walk's finite step set costs a factor 1/sqrt(1 - 1/r) in the expected number of steps (the
Brent-Pollard correction, Brent 1980; Teske 2001, section 4), which at r = 2^24 is 1 + 2^-25,
and Teske (2001) measured deviations of a few percent in the collision rate for r = 20 and
smaller r, falling as r grows; Duursma, Gaudry and Morain (1999) report the same picture for rho
variants with automorphisms.

For a random mapping the first k iterates over all trails are pairwise distinct with probability
prod (1 - i/N') <= exp(-k(k-1)/(2 N')), so with the quality factor a collision occurs among the
first k steps with probability at least

    1 - exp(-k^2 rho / (2 N')).

Take k = sqrt(2 N' ln(1/0.58)) = 1.0438 sqrt(N'), log2 k = 127.5618. At rho = 1 the collision
probability is 0.4199. The claimed success probability is 0.39, which the same expression delivers
for every rho >= 0.9076: at rho = 0.908 it is 0.3901. So the claim needs the walk to be no worse
than 90.8 % of an ideal random mapping, an allowance of 9.2 %, against measured deviations of a
few percent for r >= 20 and a 1 + 2^-25 finite-r correction. The losses below are then subtracted
from 0.3901, so the final bound is 0.39.

A collision in the walk (two walkers, or one walker twice, on a common class) merges the
continuations: after a common class, the trail up to the next distinguished point is the same
sequence of classes, because f is a function on classes and the cycle escape E(C) is a function
of the cycle (Section 2): a walker that passes through a cycle leaves it at E(C), whatever its
own detection time. Escape steps and the wasted lap steps repeat classes and are not fresh
samples; they are removed from the count by the factor 1 - 2^-9 in K below, not by the model.

Detection. A collision merges two trails; both walkers reach the next DP after a geometric number
of steps of mean 1/theta; the chance either needs more than D = 24/theta = 24 x 2^96 steps is at
most 2 e^-24 < 10^-10. Running every walker D extra steps costs m D = 2^14 x 24 x 2^96 = 2^114.585
steps. The budget is

    K = (k + m D) / (1 - 2^-9) = 2^127.5648 steps,

where 2^-9 bounds the share of steps lost to fruitless-cycle events of every length (model value
2^-9.19, Section 2). The algorithm stops at K steps and fails otherwise. Success probability at
least 0.3901 - 10^-10 - 2^-14 >= 0.39 over uniform x and fresh coins. (The 2^-14 is the chance
that the only collision in the run is inside one walker's own trail, m (k/m)^2 / (2 N'); an
own-trail collision does not produce a relation between two trails.) There is no restart or
amplification; W is the work of one capped run.

Sensitivity. The factor 1/(1 - 2^-9) contributes 0.0028 to log2 K. Every cycle term could be
replaced by a budget up to 2^-8.72 and the total W of Section 5 would still round up to 127.82,
so the claim does not rest on the exact value of the cycle bound.

## 4. Cost ledger (instructions per walker step)

Field multiplication M. `p256_mul_count.md` in this folder lists the instruction sequence and
counts it. In short: 4 limbs, Montgomery form, p = -1 mod 2^64, so the quotient digit of every
reduction round is the current low limb and m*p is a sum of shifts of m. The 4 x 4 product is
summed in column (Comba) order, 16 MUL + 16 UMULH plus 5 add-with-carry instructions per limb
product; four reduction rounds of adds, adds-with-carry, subs and subs-with-borrow (no
multiplication); one conditional subtraction of p; operand loads and result stores. The exact
count is **139 instructions**, so the ledger charges **M = 146** (5 % reserve, 139 x 1.05 =
145.95), with S = M, and a modular add/sub A = 21 (the 20-instruction addition listed there plus
the same 5 % reserve; the subtraction is 19). The model of the same sequence,
`check_mul_mont.py`, executes every listed instruction under the AArch64 flag rules and checks
the result against a b 2^-256 mod p, a + b mod p and a - b mod p on 20,000 random and
edge-case inputs each (plus 20,000 carry-rich inputs for the multiplication).
Inversion by Fermat: 255 squarings + 12 multiplications <= 300 M = 43800 instructions, amortised
over the 2^14 walkers of a round (2.67 instructions per step).

Carry chains: on AArch64 `adc` and `sbc` read but do not write the flags, so every instruction of
the listing whose carry-out is read by the next one is `adcs` or `sbcs`; `check_mul_mont.py`
models exactly this rule and asserts that the previous (plain `adc`/`sbc`) form of the chains
computes wrong products.

The first numeric column is the instructions charged; the second is an explicit upper bound
obtained by writing out the instructions of that line (itemised after the table). Every charged
line is at least its bound.

| Per-step component | Charged | Itemised bound |
| --- | ---: | ---: |
| 6 field multiplications (3 addition formula + 3 batch inversion), 6 x 146 | 876 | 834 (6 x 139) |
| 6 modular add/sub (u_T - u, v_T - v, two for u', two for v'), 6 x 21 | 126 | 120 (6 x 20) |
| zero-denominator test | 6 | 4 |
| table index, table address, loads of c_j, d_j | 12 | 3 |
| negation canonicalisation | 32 | 17 |
| coefficient update of A and B | 64 | 20 |
| distinguished-point test | 4 | 2 |
| cycle check: full 256-bit checkpoint | 20 | 15.05 |
| walker state: address formation, loop control, spills | 50 | 48 |
| amortised inversion, 300 x 146 / 2^14 = 2.67 | 3 | 2.67 |
| **Total** | **1193** | **1065.72** |

Itemisation of the bounds (every instruction named is one scalar instruction; the loads and
stores of the multiplication and add/sub operands, and the loads of u_T, v_T, are inside the 139
and 20 of the listing and not counted again here):

* zero-denominator test, 4: `orr`, `orr`, `orr` of the four limbs of d = u_T - u, `cbz`.
* table index, address, c_j and d_j, 3: `ubfx` (bits 232..255 of u), `madd` (entry address,
  entries are 80 bytes), `ldp` of c_j, d_j.
* canonicalisation, 17: 2 `ldp` of v'; `subs`, 3 `sbcs` for -v = p - v; `cmp` (`subs xzr`) and 3
  `sbcs xzr` comparing v with -v (carry flag); 4 `csel` choosing the smaller; `cinv` flipping the
  sign mask s on the carry; 2 `stp` of the chosen v. The tie v = -v does not occur.
* coefficient update, 20: the sign s is kept as a mask m (0 or all ones). For A: `eor` of c_j
  with m, `cmn m, 1` (sets C iff m = -1, the +1 of the two's complement negation), `adcs`
  of limbs 0..2 and `adc` of limb 3 with the limbs (c_j ^ m, m, m, m): 6, plus 2 `ldp` and 2 `stp`
  of the 4 limbs of A: 10. The same 10 for B with d_j. Total 20.
* distinguished-point test, 2: `ldr` of the low limb of u', `cbnz` (the second 32-bit test and the
  record are on the rare path, charged in Section 5).
* cycle check, 15.05: 2 `ldp` of the four checkpoint limbs (2); four `cmp` against the four limbs
  of u' and four `b.ne` to leave the test on the first limb that differs (8); `add` of the step
  counter (1); `cmp` and `b.eq` of the counter against 64 (2); `ldr` and `cbnz` of the cycle-walk
  mode flag (2); the refresh, once per 64 steps, two `stp` of the four new checkpoint limbs and
  one `mov` of the counter (3 / 64 = 0.05). The compare is the full-width equality of the
  checkpoint with the new u; the 20 charged leave 4.95 instructions above the itemised bound.
* walker state, 48: address formation, one `add` for each of up to 3 operand pointers of each
  of the 12 field operations of the step (36); loop control, `subs` and `b.ne` in each of the
  4 phases of the lockstep round (8); spill and reload of the table index j (2: `str`,
  `ldr`); load and store of the sign mask s (2).
* amortised inversion, 2.67: 300 x 146 / 2^14.

The Montgomery multiplication and add/sub bounds are the listed counts (139, 20); the charged 146
and 21 add a 5 % reserve (6 x 7 + 6 x 1 = 48). The ledger as charged exceeds the itemised bound
by 1193 - 1065.72 = 127.28 instructions per step (11.9 %): 48 in the field operations and 79.28
in the other lines.

One walk step costs 1193 / 1000 = 1.193 units (log2 = 0.254594).

## 5. Other charges

All items are charged in units (instructions / 1000) and each is negligible against 2^128.

* Table: 2^24 entries, each two 64-bit scalar multiplications in Jacobian coordinates with a=-3,
  by double-and-add from the affine base point. Field multiplications and field additions are both
  charged: a doubling (dbl-2001-b) is 8 M (3 M + 5 S, S = M) and 18 modular add/sub (subtractions,
  and the multiples by 2, 3, 4, 8 done by additions), a mixed addition (madd-2007-bl) is 11 M and
  14 add/sub, a general Jacobian addition (add-2007-bl) is 16 M and at most 14 add/sub; the table
  charges 20 add/sub for every one of them. With M = 146 and A = 21: doubling 8 x 146 + 20 x 21 =
  1588, addition 11 x 146 + 20 x 21 = 2026 (64 of each per scalar multiplication, which is
  more than needed for a 64-bit scalar), so 2 x 64 x (1588 + 2026) = 462,592; the sum of the two
  halves as a general addition 16 x 146 + 20 x 21 = 2756; batch conversion to affine at most
  10 M = 1460 (no additions). In total at most 466,808 instructions = 466.81 units per entry, and
  2^24 x 466.81 = 2^32.867. Coefficient generation: 2^25 random words at 20 instructions, under
  2^22 units.
* Walker starts: 2^14 walkers, each built like a table entry (two 64-bit scalar multiplications
  and their sum, 466.81 units): 2^14 x 466.81 = 2^22.87.
* DP handling: expected K theta = 2^31.565 DPs; each costs at most 4000 instructions (reduce
  A, B mod n, store 12 words, hash probe, transfer): 2^31.565 x 4 = 2^33.56 units.
* Cycle events (detection delay, minimum lap, escape addition, escape table): inside the 2^-9
  share of Section 3 (model values: 2^-16.87 for detected cycles, 2^-32.87 for re-entry into the
  same cycle, 2^-9.22 for undetected long cycles). The escape table of 256 points costs
  256 x 466.81 units = 2^16.87, in preprocessing.
* Zero-denominator slow paths: probability 2^-231 per step, under 2000 instructions each.
* Recovery and verification: two modular inversions mod n, two scalar multiplications to test
  [x]P = Q and input checks, at most 10^4 units (< 2^14).

Walk work: K x 1.193 = 2^127.5648 x 2^0.254594 = 2^127.8194. Everything else sums to under
2^34.3 (table 2^32.867, DP handling 2^33.565, starts and the rest below 2^23), adding a relative
2^-93.5. Total W = 2^127.8194, so the claimed score is **time_log2 = 127.82** (rounded up from
127.8194). Preprocessing (table, escape table, walker starts, coefficient generation) is
2^32.869, at most 2^32.9 units, included in W. Advice: none, 2^0 bytes.

## 6. Heuristics, support and limits

1. Random-mapping behaviour of the r-adding walk on classes (Teske 2001; Bos et al. 2010 for the
   negation map), with the quality factor rho >= 0.908 stated in Section 3. Support: r = 2^24 is
   far above the r >= 16 range tested, the finite-r correction is 1 + 2^-25, and the toy
   experiment of Section 7 measures a collision ratio of 1.015 (95 % interval 0.949 to 1.082) on
   240 samples over four toy tables, i.e. a walk within the claimed band at toy scale. Limit: no
   experiment at P-256 scale exists, the toy's table is far smaller relative to its group than
   production's, and the claim is analytic.
2. Cycle model. A 2-cycle has probability 1/(2R) per step. A cycle of 2l classes (l >= 2) has
   per-step probability at most 2 (2l/R)^l (Section 2); expected counts over K = 2^127.5648 steps
   are 2^84.6 (l = 2), 2^64.3 (l = 3), 2^44.6, 2^25.2, 2^6.1 (l = 6) and 2^-12.8 (l = 7); the
   sum over l >= 7 is below 2^-12. So cycles of up to 12 classes do occur and, like every cycle of
   at most 64 classes, are detected, walked once and escaped from their minimum member. Cycles
   longer than 64 are not detected and are bounded in Section 2: random-mapping closures give
   2^-28 for any repeat, and zero-sum closures give 2^-9.22 of all steps, so the whole cycle
   loss is at most 2^-9.19 and Section 3 budgets 2^-9. Re-entry into a cycle just escaped is
   bounded by the basin argument of Section 2 at 2^-135.435 per cycle and 2^-32.87 of all steps,
   and the escape table T' is independent of the main table. Remaining weak points: the counting
   bound 2 (2l/R)^l is a counting argument under independence of table indices, not a theorem;
   the basin calculation is a random-mapping expectation, not a worst case; the zero-sum
   closure model is an estimate, which is why the budget sits 1.15x above it and why Section 3
   shows the score is insensitive to it up to 2^-8.72; and no run at P-256 scale checks any of it.
3. Bit fields. The DP test (bits 0..95 of u), the escape index (bits 224..231) and the table index
   (bits 232..255) are disjoint, so a distinguished point or an escape point is not forced to a
   particular index, and the index is a function of bits the DP test does not constrain. The
   independence of the two values for a well-mixed u is a heuristic. DP storage is a hash or sort
   charged per record; no free oracle.
4. Arithmetic count. The six field multiplications and six add/sub operations of a step are
   counted from an explicit instruction listing, `p256_mul_count.md`, and checked by executing
   that listing (`check_mul_mont.py`) under the AArch64 flag rules (`adds`, `adcs`, `subs`,
   `sbcs` write the carry flag, `adc` and `sbc` do not). The count is static: it is not a
   measurement of a compiled kernel and makes no timing claim, and the 5 % reserves on M and A
   cover scheduling and operand routing, not spills. The other ledger lines are not estimates
   with an error bar: each is stated as an upper bound and itemised instruction by instruction in
   Section 4 (zero-denominator test 4, table index 3, canonicalisation 17, coefficient update 20,
   DP test 2, cycle check 15.05, walker state with address formation for all 12 field operations
   48, inversion 2.67; total 111.72 against 191 charged). The charged 1193 exceeds the itemised
   1065.72 by 127.28 instructions per step, 11.9 %. The slack argument: the charged lines are the
   itemised upper bounds plus that margin, so the claim holds unless the itemised instructions
   are undercounted in total by more than 127 per step. What the bounds do not cover is a
   compiler's register spills and scheduling inside the field operations, which this
   instruction-count model does not price.
5. Memory: expected 2^31.57 DPs; allow 2^34 (Poisson concentration) at 96 bytes = 2^40.59 bytes,
   table 2^24 x 80 B = 2^30.3 bytes, walkers 2^14 x 480 B = 2^22.9 bytes (each walker holds u and
   v, the 4-limb A and B, the 32-byte checkpoint, the counter, the sign mask, the table index and
   the cycle-minimum snapshot): peak at most 2^41 bytes (about 2 TB). Memory carries no scalar
   cost.
6. Nothing here is a novel attack: it is the textbook generic rho with the negation map,
   distinguished points and a published cycle-escape rule, priced in instructions.

## 7. Toy experiment (`experiment.py`)

`experiment.py` runs the same walk, the same checkpoint rule and the same escape rule on a toy
group: y^2 = x^3 - 3x + 4 over F_65537 (a = -3), prime order 65851, N' = 32926 classes, theta =
2^-6, checkpoint period 16, table index, escape index and DP test on disjoint bit fields, escape
from the cycle's smallest-x member with a separate escape table, and the zero-denominator case
handled by a doubling. Its default inputs run in under a second and report, at toy scope only:

* recovery 29 of 30 end-to-end attempts (fresh target, table and walkers each time); the one
  failure found no usable distinguished-point collision inside its step budget;
* merge 106 of 120 pairs reach the same distinguished point and 14 of 120 are both looping in one
  cycle, the case where the escape point fell back into the basin of the cycle just left, with 0
  of 120 diverging: a pair at a common class with different checkpoint histories never separated;
* 675 escape events over 82 fresh tables, of which 17 walkers (of 400 that escaped at least once)
  came back to the cycle they had just escaped, a toy basin fraction of about 2^-4.5 against the
  2^-135.4 of Section 2, which is why the trap is visible only at toy scale;
* 5 zero-denominator steps, handled by the doubling slow path;
* a collision ratio of 1.015 against sqrt(pi N'/2), 95 % interval 0.949 to 1.082, over 240
  samples and four tables of 1024 entries. The toy's cycle overhead is still far above
  production's (Section 2), so this ratio is an upper bound on the production walk's overhead.

Toy scale is the scope of all five statements; none of them certifies a production exponent. The
experiment exists because the cycle rule, its history independence and the re-entry trap are
mechanisms that an argument should be able to exhibit, not only bound.

References: P. van Oorschot and M. Wiener, parallel collision search with cryptanalytic
applications, J. Cryptology 1999 (distinguished points); J. M. Pollard, a Monte Carlo method for
factorization, BIT 1975, and R. P. Brent, an improved Monte Carlo factorization algorithm, BIT
1980 (the rho walk and its cycle finding); E. Teske, on random walks for Pollard's rho method,
Math. Comp. 2001 (r-adding walks, r >= 16, measured deviations and the finite-r step-count
correction); I. Duursma, P. Gaudry and F. Morain, Speeding up the discrete log computation on
curves with automorphisms, ASIACRYPT 1999 (collision statistics of rho variants); Bos, Kaihara,
Kleinjung, Lenstra and Montgomery, 2010 (fruitless cycles with the negation map); D. Bernstein,
T. Lange and P. Schwabe, "On the correct use of the negation map in the Pollard rho method",
PKC 2011.

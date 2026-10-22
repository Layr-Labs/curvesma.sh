# p256-ecdlp: parallel Pollard rho with distinguished points and the negation map

Analytic cost bound; no execution at production size is claimed. All counts are scalar AArch64 instructions divided by C = 1000, as the track prescribes. The
unit C is a normalization only: every P-256 operation below is charged in instructions.

## 1. Problem and notation

Input Q = [x]P, x uniform in [0, n), n the subgroup order, log2 n = 255.999999999664. The
negation map identifies R with -R. The number of classes is N' = (n+1)/2 (the identity is a
singleton class), log2 N' = 254.999999999664.

## 2. Algorithm

Walk function. A table of R = 2^24 entries T_j = c_j P + d_j Q is built once, with c_j, d_j
uniform nonzero 64-bit integers and T_j stored in affine form with its c_j, d_j. For a class
representative X with affine coordinates (u, v) in Montgomery form, define
j(X) = bits 40..63 of the top limb of u (identical for X and -X), and
f(X) = canon(X + T_j(X)), where canon picks of {Y, -Y} the one whose Montgomery y-limbs are the
smaller as a 256-bit integer. f is a deterministic function on classes (an r-adding walk, Teske
2001, with r = 2^24 >= 16).

Walkers. m = 2^14 walkers start at random points a0 P + b0 Q (a0, b0 uniform 256-bit) and run in
lockstep. Every walker keeps a canonical
point (u, v), a sign s = +-1 and signed integers (A, B) with the invariant
point = s (A P + B Q). A step adds T_j: A += s c_j, B += s d_j (signed 4-limb integers,
never reduced mod n during the walk; |A|, |B| < 2^180), and if canon negates the sum, s := -s.
A walker never restarts.

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
a step counter c and a 64-bit checkpoint. When c reaches 64 the checkpoint is set to the low
limb of u of the current point and c := 0; every new point is compared to the checkpoint (full
256-bit compare on a limb match). A match means the walk has returned to a class it held at the
checkpoint; because f is deterministic on classes, the walker is then on a cycle of f, and the
walker's own schedule has no further effect on its trail, as shown next. A cycle of length L <= 64 that a walker
has entered is detected within 128 steps of entry: a checkpoint is set within the first 64
steps, and the walker returns to it within L <= 64 further steps, before the next refresh.

Escape as a function of the cycle. On detection the walker is at a cycle member X0 (a canonical
point) and walks the cycle once more in the ordinary lockstep batch (L further steps, L <= 64,
ending again at X0), additionally keeping the member with the smallest u (compared as a 256-bit
integer in the stored Montgomery representation) together with its full state (u, v, A, B, s).
Distinct classes of a cycle have distinct u, so the minimum m* is unique and depends only on the
set of the cycle's classes, not on X0, on the walker, on its counter or on its checkpoint.
Let e(m*) be bits 32..39 of the top limb of u(m*). The walker restores the snapshot of m* and
takes one stand-alone affine addition (its own inversion) of the escape point T'_{e(m*)} from
a second table of 256 points, built like the main table with its own nonzero 64-bit (c', d'),
updating (A, B, s) exactly as a normal step. It then sets c := 0 and continues. The result
E(C) = canon(m* + T'_{e(m*)}) is a function of the cycle C only. Two walkers that merged
earlier are in the same cycle, possibly at different members, after different numbers of laps
and with different counters; each one detects at its own time, and both leave from the same
point m* to the same point E(C). They stay merged. Adding T' (64-bit coefficients) rather than
doubling keeps |A|, |B| bounded without reduction. Should E(C) land in a cycle again (probability
about 2^-25, a new event), the rule applies again.

Cost of one cycle event, in step equivalents of 1664.03 instructions: at most 128 steps of
detection delay (steps already taken inside the cycle), at most 64 steps of the minimum-tracking
lap, charged 40 more instructions each for the 4-limb compare and the rare snapshot copy
(64 x 1704.03 / 1664.03 < 66), and the escape addition at most 100 units = 100000 / 1.66403 / 1000
= 60.1 step equivalents; in total at most 128 + 66 + 61 = 255 < 256 step equivalents. With the
per-step event rate lambda <= 2^-25 + sum over l >= 2 of 2 (2l/R)^l <= 2^-25 + 2^-42 (the l = 2
term is 2^-43), the amortised loss is at most 256 x 2^-24.9999 = 2^-16.99 of the steps, about
1.2 x 10^-5. Section 3 charges 2^-12, a factor 32 above this. Cycles of length above 64 are not
detected; their rate is below 2^-590 per step (Section 6).

## 3. Success probability

Walks are modelled as a random mapping on N' classes (Teske 2001; r = 2^24, with the negation
map and cycle escapes as above). For a random mapping, the first k iterates over all trails are
pairwise distinct with probability prod (1 - i/N') <= exp(-k(k-1)/(2 N')), so a collision
occurs among the first k steps with probability at least 1 - exp(-k^2/(2N')) (the k(k-1) versus
k^2 difference is below 2^-127). Take

k = sqrt(2 N' ln(1/0.58)) = 1.0438 sqrt(N'), log2 k = 127.5618,

so the collision probability is at least 0.42 under the model. The 0.03 margin over 0.39
absorbs walk non-randomness (Teske reports deviations of a few percent for r = 20; here
r = 2^24), a collision inside one walker's own trail (probability about m (k/m)^2/(2N') =
2^-14), and the escape steps. A collision in the walk (two walkers, or one walker twice, on a
common class) merges the continuations: after a common class, the trail up to the next
distinguished point is the same sequence of classes, because f is a function on classes and
the cycle escape E(C) is a function of the cycle (Section 2): a walker that passes through a
cycle leaves it at E(C), whatever its own detection time. Escape steps and the wasted lap steps
repeat classes and are not fresh samples; they are removed from the count by the factor
1 - 2^-12 in K below, not by the model.

Detection. A collision merges two trails; both walkers reach the next DP after a geometric
number of steps of mean 1/theta; the chance either needs more than D = 24/theta = 24 * 2^96
steps is at most 2 e^-24 < 10^-10. Running every walker D extra steps costs m D = 2^14 * 24 * 2^96
= 2^114.585 steps. The budget is

K = (k + m D) / (1 - 2^-12) = 2^127.5623 steps,

where 2^-12 is a bound on the share of steps lost to fruitless-cycle events (model value 2^-16.99,
Section 2 and 6). The
algorithm stops at K steps and fails otherwise. Success probability >= 0.42 - 10^-10 - 2^-14 >=
0.39 over uniform x and fresh coins. There is no restart or amplification; W is the work of one
capped run.

## 4. Cost ledger (instructions per walker step)

Field multiplication M (4 x 64-bit limbs, Montgomery, p = -1 mod 2^64 so the quotient digit of
every reduction round is the current low limb): schoolbook product 16 MUL + 16 UMULH + about
44 ADDS/ADCS/ADC for carries; reduction four rounds, each a few shifts/adds because
p = 2^256 - 2^224 + 2^192 + 2^96 - 1 (m*p is formed from m shifted by 32 and carried), about 16
instructions per round; final conditional subtraction about 10; operand loads 8, result stores 4,
call and return 4. Total about 185. The ledger charges **M = 220** (margin above 18%), S = M, and
a modular add/sub A = 24 (4 chained ADDS/ADCS, mask, up to 3 masked ANDs of p's limbs, 4 chained
corrections, constants). Inversion by Fermat: 255 squarings + 12 multiplications <= 300 M = 66000
instructions, amortised over the 2^14 walkers of a round (4.0 instructions per step).

| Per-step component | Instructions |
| --- | ---: |
| 6 field multiplications (3 addition formula + 3 batch inversion), 6 x 220 | 1320 |
| 6 modular add/sub (u_T - u, v_T - v, two for u', two for v'), 6 x 24 | 144 |
| zero-denominator test (OR of 4 limbs, branch) | 6 |
| table index (UBFX), address, loads of u_T, v_T, c_j, d_j | 24 |
| negation canonicalisation: 4-limb compare, conditional negation of v (sub chain + 4 CSEL), flip s | 32 |
| coefficient update of A and B (sign-masked signed 4-limb add) including their load/store | 64 |
| distinguished-point test (cbnz on low limb; second test on rare hit) | 4 |
| cycle check: counter, checkpoint compare and refresh, cycle-walk mode flag | 16 |
| walker state, prefix-product traffic, loop control | 50 |
| amortised inversion, 66000 / 2^14 | 4 |
| **Total** | **1664.03** |

One walk step costs 1664.03 / 1000 = 1.664028 units (log2 = 0.734681).

## 5. Other charges

All items are charged in units (instructions / 1000) and each is negligible against 2^128.

* Table: 2^24 entries, each two 64-bit scalar multiplications in Jacobian coordinates with a=-3
  (64 doublings x 8M and at most 64 additions x 11M each, 2432 M), one sum 11 M, batch conversion
  to affine 10 M: at most 2500 M = 550 units per entry, 9.23 x 10^9 units = 2^33.10. Coefficient
  generation: 2^25 random words at 20 instructions, under 2^22 units.
* Walker starts: 2^14 walkers, each two 256-bit scalar multiplications, at most 10000 M = 2200
  units: 3.6 x 10^7 units = 2^25.1.
* DP handling: expected K theta = 2^31.562 DPs; each costs at most 4000 instructions (reduce
  A, B mod n, store 12 words, hash probe, transfer): 2^31.562 x 4 = 2^33.56 units.
* Cycle events (detection delay, minimum lap, escape addition, escape table): inside the 2^-12
  share of Section 3 (model value 2^-16.99 = 1.2 x 10^-5, 255 step equivalents per event at
  2^-25 events per step). The escape table of 256 points costs 256 x 550 units = 2^17.1.
* Zero-denominator slow paths: probability 2^-231 per step, under 2000 instructions each.
* Recovery and verification: two modular inversions mod n, two scalar multiplications to test
  [x]P = Q and input checks, at most 10^4 units (< 2^14).

Walk work: K x 1.664028 = 2^127.5623 x 2^0.734681 = 2^128.2970. Everything else sums to under
2^34.1, adding a relative 2^-94. Total W = 2^128.2970, so the claimed score is
**time_log2 = 128.30** (rounded up from 128.297). Preprocessing (table, starts, coefficient
generation) is at most 2^33.2 units, included in W. Advice: none, 2^0 bytes.

## 6. Heuristics, support and limits

1. Random-mapping behaviour of the r-adding walk on classes (Teske 2001; Bos et al. 2010 for the
   negation map). Support: r = 2^24 is far above the r >= 16 range tested; the margin of 0.03
   in probability is 7%. Limit: no experiment at P-256 scale exists; the claim is analytic.
2. Cycle model. A 2-cycle has probability 1/(2R) per step. A cycle of 2l classes (l >= 2) has
   per-step probability at most 2 (2l/R)^l (Section 2); expected counts over K = 2^127.5623 steps
   are 2^84.6 (l = 2), 2^64.3 (l = 3), 2^44.6, 2^25.2, 2^6.1 (l = 6) and 2^-12.8 (l = 7); the
   sum over 7 <= l <= 32 is below 2^-12. So cycles of up to 12 classes do occur and, like every
   cycle of at most 64 classes, are detected, walked once and escaped from their minimum member.
   For l >= 33 up to 2^20 the rate 2 (2l/R)^l is below 2^-590 per step at l = 33 and falls with l; longer cycles are the
   random mapping's own cycles and lie inside the collision model. The lost share is 2^-16.99
   (255 step equivalents x 2^-24.9999 events per step), bounded by 2^-12. Remaining weak
   point: the bound 2 (2l/R)^l is a counting argument under independence of table indices, not
   a theorem; no run at P-256 scale checks it.
3. Independence of the DP test (low 96 bits of u) from the table index (bits 40..63 of the top
   limb) for a well-mixed u. DP storage is a hash or sort charged per record; no free oracle.
4. Memory: expected 2^31.56 DPs; allow 2^34 (Poisson concentration) at 96 bytes = 2^40.59 bytes,
   table 2^24 x 80 B = 2^30.3 bytes, walkers 2^14 x 400 B = 2^22.6 bytes: peak at most 2^41 bytes
   (about 2 TB). Memory carries no scalar cost.

References: P. van Oorschot and M. Wiener, parallel collision search with cryptanalytic
applications, J. Cryptology 1999 (distinguished points); E. Teske, on random walks for Pollard's
rho method, Math. Comp. 2001 (r-adding walks, r >= 16); Bos, Kaihara, Kleinjung, Lenstra and
Montgomery, 2010 (fruitless cycles with the negation map); D. Bernstein, T. Lange and P. Schwabe,
"On the correct use of the negation map in the Pollard rho method", PKC 2011.

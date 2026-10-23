# Plain parallel Pollard rho on P-256 (no negation map), with a round meter

**Claim.** time_log2 = 13047/100 = 130.47. The algorithm below recovers x from
Q = [x]P on NIST P-256 with probability at least 0.39 over a uniform challenge
and fresh algorithm coins. On **every** run, its total work is at most
W_max < 2^130.4566 target-group-operation equivalents (C = 1,000 scalar AArch64
instructions per unit). Peak memory < 2^99 bytes. Preprocessing (everything before
the first walk round) < 2^34.2 units, included in W. Nonuniform advice: none.

The work bound is deterministic. It holds because a round counter and three
event counters stop the attack (section 3). The success probability uses one
disclosed random-mapping heuristic (H1). The claimed success needs only κ ≤ 2.
That allows a factor 2 in the birthday constant, against a predicted correction
of about 1 + 1/R = 1.0005 (section 4). H1 is supported by end-to-end toy runs
(section 7).

**Design choice.** The walk works on plain group elements. It uses no negation
map, no equivalence classes and no history-dependent escapes. The rejected
negation-map submissions on these tracks were refuted on fruitless cycles,
escape traps and merging under history-dependent escapes. None of that applies
here. The walk is a fixed function F on the group, so two trails that reach the
same point continue identically. The price is about 0.5 bit relative to a
negation-map walk.

## 1. Notation and target

p, a = -3, b, r, P are the protected SEC 2 P-256 parameters
(`tracks/p256-ecdlp/target-profile.json`). r is prime, log2 r = 256 − 3.4·10^-10,
and the cofactor is 1. `experiment.py` checks that P is on the curve, r is a
probable prime and [r−1]P + P = O. The input Q is validated: canonical
encoding, on the curve, and in ⟨P⟩, which is all of E(Fp) because the cofactor
is 1. If Q = O, the attack outputs 0. Otherwise x ≠ 0.

Field elements are stored in **Montgomery form** x̃ = x·2^256 mod p as canonical
4×64-bit limbs in [0, p). The map x ↦ x̃ is a fixed bijection of Fp, so any
function of x̃ is a function of the group element. Inputs are converted once and
outputs once; that cost is included in the start and final allowances.

## 2. Algorithm

Parameters: R = 2048 table entries. Walkers M = 2^20, run as 2^8 processors
each holding B = 2^12 walkers. Distinguished-point (DP) probability θ = 2^-40.
Trail cap L_MAX = 64/θ = 2^46. Tail rounds L_TAIL = 50/θ.

**Table.** For h = 0..R−1, draw c_h, d_h uniform in [0, r) from fresh coins
(rejection sampling from 256-bit words) and set M_h = [c_h]P + [d_h]Q. If M_h = O,
redraw. Store M_h in canonical affine Montgomery form together with (c_h, d_h).

**Walk function.** For a finite point X = (x, y), let h(X) = x̃ mod 2048 (the low
11 bits of limb 0). Then

    F(X) = X + M_{h(X)}.

F is a fixed function on the group, determined by Q and the table. A point X is
**distinguished** if bits 0..39 of limb 1 of x̃ are zero (probability θ). These
bits are disjoint from the index bits.

**Trails.** A trail starts at S = [s]P + [t]Q with s, t uniform in [0, r) from
fresh coins (redraw if S = O). Only the start (s, t) is stored, not running
coefficients. Each round applies F once to every live walker. After each step,
the walker's new point X is tested, in this order:

1. If X is distinguished, report the record (key = x̃ of X, s, t) to the store.
   The trail then ends and the walker restarts from a fresh start.
2. Otherwise, if the trail has made L_MAX steps, it is abandoned and the walker
   restarts.

The start point itself is not tested. Distinguishedness is tested only on step
outputs.

**Batch inversion (one round on one processor).** Each walker needs
λ = (y_M − y)/(x_M − x) with M = M_{h(X)}. First compute every denominator
δ_i = x_M − x and test it for zero. If δ_i = 0, X = ±M_h is an exceptional
input; that walker restarts and its new δ_i is formed, so all δ_i are nonzero.
Then apply Montgomery's simultaneous inversion to the 2^12 denominators:
prefix products π_i = π_{i−1}·δ_i (1M each), one inversion of π_B, and the
backward pass inv_i = acc·π_{i−1}, acc ← acc·δ_i (2M each). Every denominator
in a round is available at the start of the round, because the walkers are
independent. This avoids the sequential-dependency objection raised against
batched single-chain BSGS.

**Affine addition.** λ = (y_M − y)·inv_i (1M), x3 = λ² − x − x_M (1S),
y3 = λ(x − x3) − y (1M). Together with the 3M of batch amortisation, each step
costs 6 multiplications/squarings and 6 additions/subtractions: y_M − y, x_M − x,
λ² − x, − x_M, x − x3, − y.

**Store and recovery.** The DP store is a balanced search tree (red-black) keyed
by the 256-bit key. On a report with key K:

- if K is absent, insert (K, s, t);
- if K is present with the same (s, t), ignore the report;
- otherwise, run recovery.

Recovery replays both trails from their stored starts with the same F. It tracks
coefficients A ← A + c_h, B ← B + d_h mod r, so that every replayed point equals
[A]P + [B]Q. Both replays end at their first DP, which is K. If B1 ≠ B2, output
x = (A2 − A1)/(B1 − B2) mod r, verify [x]P = Q, and halt. If B1 = B2, continue.

The recovery formula follows from [A1]P + [B1]Q = [A2]P + [B2]Q, which gives
A1 + B1·x ≡ A2 + B2·x, so x ≡ (A2 − A1)/(B1 − B2). The toy runs execute this
exact rule, with 0 wrong outputs (section 7).

## 3. Meter and the deterministic work bound

The attack runs exactly N = ⌈S1/M⌉ + L_TAIL main-plus-tail rounds and then
halts with failure if it has not succeeded. Here S1 is defined in section 4.
Counters cap the non-step events. If an increment would exceed its cap, the
attack halts with failure before performing the event.

| event | cap |
|---|---|
| DP reports (insertions or lookups) | D_CAP = 4·⌈S/2^40⌉ + 2^30 (< 2^90.52) |
| trail starts (initial, after DP, after abandon or zero denominator) | M + D_CAP + 2^20 |
| recoveries (pairs replayed) | 2^20 |

S = S1 + M·(L_TAIL + 1) ≥ M·N is the number of walk steps. Every action is
charged at a price that upper-bounds its instructions (section 5):

    W_max = S·c_step + D_CAP·c_DP + STARTS·c_start + PRE + 2^20·c_rec + c_final.

`experiment.py` (mode production) evaluates this exactly with rationals:

| term | log2 units |
|---|---|
| walk steps S·c_step | 130.456564 |
| DP store D_CAP·c_DP | 96.55 |
| starts | 104.68 |
| recoveries | 76.05 |
| preprocessing (table, init) | 25.2 |
| final verification | 13.2 |
| **total W_max** | **130.4565639** |

The script checks the exact inequality W_max^100 ≤ 2^13047 with W_max as a
rational, and it holds. The counters, the round counter and their comparisons
are part of the per-step control allowance and the per-event allowances.

## 4. Success probability

**Heuristic H1 (walk quality).** In the parallel process above, let T be the total
number of steps executed (over all walkers, in round order) until some group
element is produced as a step output by two **different** trails. Then for
t ≤ 2^129:

    Pr[T > t] ≤ exp(−t² / (2κr)),  with κ = 2.

For an ideal random mapping, κ = 1. This is the birthday bound behind van
Oorschot–Wiener parallel collision search. For R-adding walks, Teske's
experiments (Math. Comp. 2001, "On random walks for Pollard's rho method") found
behaviour indistinguishable from random for R ≥ 16. Heuristic analyses of
R-adding walks, including the Brent–Pollard style "anti-collision" correction
discussed by Bernstein and Lange, predict a constant change of order 1/R, which
is about 0.0005 here. κ = 2 allows a deviation about 2,000 times larger. The
toy measurements in section 7 give κ̂ between 0.8 and 1.1 at r ≈ 2^20..2^32,
with R = 2048 and also with R = 16. These are evidence at their own scale only.

**Main phase.** Let S1 = ⌈√(2κ r · 0.51083)⌉ + 2, so log2 S1 = 128.5155. Since
ln(1/0.6) < 0.51083, H1 gives Pr[T ≤ S1] ≥ 0.4. The main phase executes at
least S1 steps.

**From a repeat to recovery.** Suppose the first repeat is at a point X, produced
first by trail A and later (or later in the same round) by trail B ≠ A.

1. **Merging.** From X onward both trails apply the same function F, so they visit
   the same points until A's first DP D after X. A stops there and stores its
   record. B, starting from X, reaches the same D after the same number k of
   further steps. Because A produced X in an earlier round, A reports D first.
   Within one round, reports are inserted sequentially, so B then finds A's key.
   No escapes, classes or histories exist that could make the two trails diverge.
2. **B reaches D in time.** B needs k more steps. Under the usual DP heuristic,
   which treats distinguishedness of fresh F-outputs as θ-coins (part of H1's
   random-function model), Pr[k > 50/θ] ≤ e^-50. The tail gives every walker
   50/θ further rounds after the main phase.
3. **Neither trail is abandoned.** A trail is abandoned only after 64/θ steps
   without a DP: probability ≤ e^-64 per trail under the same model. A trail's
   length also cannot be inflated by a structural cycle. A cycle of F of length
   k needs a nonzero nonnegative relation Σ n_h M_h = O with Σ n_h = k. For each
   fixed relation (n_h < r), uniform c_h makes it hold with probability exactly
   1/r. A union bound over all relations with k ≤ 28 gives at most 2^-45.66
   (exact value computed in `experiment.py`). Longer self-cycles are ordinary
   self-collisions under H1. Over all ≲2^90.5 trails of length ≤ 2^46, the
   expected number is ≤ 2^(90.5+92−256) ≈ 0.
4. **Distinct coefficients.** Conditioned on all group elements (table points,
   starts, hence every walk), the start coefficient t of each trail stays
   uniform and independent. Given S, t determines s, and the table coefficients
   d_h are independent of the points M_h because c_h is uniform. The replayed
   B equals t plus a sum of d_h along a path fixed by the group elements. So
   B1 = B2 has probability exactly 1/r. A start collision S_A = S_B has
   probability ≤ 2^180/2^256.
5. **Caps.** D_CAP is at least four times the heuristic expected number of DPs.
   The start and recovery caps are far above their expected values (the
   expected number of abandons and zero denominators is < 1). Exceeding a cap
   has negligible probability under the same model.

Hence Pr[success] ≥ 0.4 − e^-50 − 2^-45 − 2^-60 > 0.399 ≥ 0.39. The only
heuristic is H1 together with the DP-coin model it includes. The identity
challenge is solved directly.

**What would invalidate the claim?** H1 failing by more than a factor 2 in the
birthday constant at 2^128 scale for this walk. A negation-free R-adding walk
with R = 2048 independent random multipliers has no known mechanism for that.

## 5. Instruction allowances (P-256)

The track's rules require P-256 arithmetic to be charged in instructions. The
allowances below are analytic bounds; each is at least twice its itemised count.

**Modular addition or subtraction, A = 64.** Load 8 limbs. ADDS + 3·ADCS. Capture
the carry (1). Subtract p: SUBS + 3·SBCS + 1 SBCS combining the carry. Select
with 4·CSEL. Store 4 limbs. That is 8+4+1+5+4+4 = 26, plus ≤ 6
addressing/return instructions, so ≤ 32. Subtraction is symmetric (borrow, then
conditionally add p).

**Modular multiplication or squaring, M = 512 (Montgomery).**

- Schoolbook 4×4 product in column order with a 3-word accumulator: 16 pairs of
  MUL/UMULH. Per pair: 2 loads, MUL, UMULH, ADDS, ADCS, ADC = 7. Plus 2 shifts
  per column for 8 columns. Total ≤ 16·7 + 16 = 128.
- Montgomery reduction, 4 rounds. For P-256, p ≡ −1 (mod 2^64), so
  m = −p^-1·t0 mod 2^64 = t0, which costs 0 instructions. Adding m·p uses the
  3 nonzero limbs of p (limb 2 is 0): ≤ 3 MUL/UMULH pairs (6), ≤ 10
  ADDS/ADCS/ADC into the 5-word window and ≤ 4 propagation instructions. That
  is ≤ 20 per round, ≤ 80 total.
- The result is < 2p. One conditional subtraction (4 SUBS/SBCS, borrow test,
  4 CSEL), 4 stores and ≤ 6 setup/return instructions: ≤ 19.

The itemised total is ≤ 227, and the allowance is 512. Squaring uses the same
bound. Constant multiplication by the Montgomery factor for conversions also
uses M.

**Inversion, I = 2^19.** Fermat z^(p−2) by square-and-multiply: ≤ 255 squarings
and ≤ 256 multiplications, which is ≤ 511·512 < 2^18. Plus ≤ 64 control
instructions per exponent bit (2^14) and ≤ 2^14 setup. The total is < 2^19.

**Per-step control, CTRL = 256.** The itemised count is about 61:

| item | instructions |
|---|---|
| load current x, y | 8 |
| index: AND, shift and add for the table address | 3 |
| load the table point | 8 |
| zero test of δ: 3 ORR, 1 branch | 4 |
| batch arrays: store δ and π, reload π and δ | 16 |
| DP test: load limb, TST, branch | 3 |
| store the new point | 8 |
| trail-length counter (increment, compare, branch) | 3 |
| loop, addressing and round counter | 8 |

The allowance is 256, which also covers meter comparisons.

**Per-step price:**

    c_step = (6M + 6A + I/B + CTRL)/C = (3072 + 384 + 128 + 256)/1000 = 3.84 units.

**Other prices** (in units; each is generous):

- **c_DP = 2^16/1000.** One red-black-tree search and insertion on ≤ 2^91
  keys. The height is ≤ 2·91 + 2 levels, at ≤ 16 instructions per level for a
  4-limb compare and child load: ≤ 2,944. Rebalancing is ≤ 3 rotations plus
  recolouring along the path: ≤ 2,000. Transfer of a 13-word record is ≤ 64.
  The total is < 2^13, so the allowance is 8 times the count.
- **One scalar multiplication**, SMUL = 256·2·32·M + I. That is ≤ 256 doublings
  and ≤ 256 mixed additions in Jacobian coordinates, each charged 32
  multiplication allowances, which covers its ≤ 16 multiplications and its
  additions, plus one inversion back to affine.
- **c_start = (2·SMUL + 6M + I + 2^12)/1000**: two scalar multiplications, the
  addition, the normalisation, and the random words with rejection sampling.
- **c_rec = 2·(c_start + (L_MAX + 2)·c_rstep).** A replayed step is a lone
  affine addition with its own inversion plus two coefficient additions mod r:
  c_rstep = (6M + 6A + I + 4A + CTRL)/1000.
- **PRE = 2048·c_start + 2^30/1000.** This covers table construction and
  initialisation. With the first 2^20 trail starts, preprocessing is
  < 2^34.2 units.
- **c_final = (SMUL + I + 2^12)/1000.** Verification and output.

All of these terms are included in W_max (section 3).

## 6. Memory

The DP store holds ≤ D_CAP < 2^90.52 records. Each record is ≤ 16 words (key,
s, t, tree links and colour), 128 bytes. That is < 2^97.52 bytes. Walker state,
the table and batch arrays add < 2^40 bytes. Peak memory < 2^99 bytes, so
memory_log2_bytes = 99. Memory has no score contribution.

## 7. Evidence (`experiment.py`)

**Production mode.** Curve checks, the exact rational ledger and the exact
inequality W_max^100 ≤ 2^13047. It also computes the union bound on short
structural cycles.

**Toy mode.** This runs the specified algorithm end to end on prime-order
y² = x³ − 3x + b curves (same a = −3; order checked by r prime, [r]G = O and
r² > 16p). It uses the same F with x mod R, the batch inversion with zero test,
DPs, restarts, the L_MAX cap, the store, replay recovery and verification. It
also keeps a full map of visited points, which gives the exact first step at
which two different trails produced the same point (T in H1). κ̂ = mean(T²)/(2r)
estimates κ, which is 1 for a random mapping. The budget S1 is computed with
κ = 2, as in production.

The initial run (r ≈ 2^20, 40 trials) completes in ~2 s. Offline runs of the
same script with the inputs stated in README.md are given below. "Merge within
S1" means a verified recovery whose first inter-trail repeat came within S1.
For κ = 1, the model value is 1 − exp(−S1²/(2r)) ≈ 0.64.

| log2 r | R | walkers | θ | trials | verified recovery | merge within S1 | κ̂ | κ̂ excl. trapped-trail steps | trapped trails | wrong | B1=B2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 2048 | 16 | 2^-4 | 400 | 400/400 | 0.695 | 0.878 | 0.876 | 9 | 0 | 0 |
| 24 | 2048 | 32 | 2^-5 | 300 | 300/300 | 0.623 | 1.017 | 1.017 | 1 | 0 | 0 |
| 24 | 16 | 32 | 2^-5 | 300 | 300/300 | 0.603 | 1.010 | 1.009 | 3 | 0 | 0 |
| 28 | 2048 | 128 | 2^-6 | 200 | 200/200 | 0.670 | 0.904 | 0.904 | 2 | 0 | 0 |
| 28 | 2048 | 128 | 2^-6 | 200 | 200/200 | 0.690 | 0.905 | 0.905 | 2 | 0 | 0 |
| 32 | 2048 | 256 | 2^-6 | 60 | 60/60 | 0.683 | 0.840 | 0.840 | 0 | 0 | 0 |
| 32 | 2048 | 256 | 2^-6 | 60 | 60/60 | 0.750 | 0.851 | 0.851 | 0 | 0 | 0 |

In the toy runs the tail (50/θ rounds of all walkers) is long compared with S1,
so almost every trial eventually recovers x. The column that tests the claim is
"merge within S1": a verified recovery whose first inter-trail repeat happened
within the κ = 2 budget S1. It is 0.60–0.75 in every row, against the 0.4 that
the proof needs. No wrong output and no B1 = B2 event occurred. At r ≈ 2^20..2^24 with R = 2048,
the table is out of the production regime: R² is comparable to r, so short
table relations Σ n_h M_h = O exist and a few trails become trapped. These are
counted in the table and abandoned at L_MAX. In production, the union bound in
section 4 excludes such relations up to length 28.

Disclosure: an earlier sweep at r ≈ 2^32 used θ = 2^-10 with only 64 walkers,
so L_MAX = 2^16 was comparable to S1 ≈ 2^16. One trapped trail could then waste
most of the budget, and the secp256k1-shaped run gave κ̂ = 1.39 (80 trials).
The rows above use toy parameters with L_MAX ≪ S1, matching the production
proportion (2^46 against 2^128.5). The script also reports κ̂ with the steps of
trapped trails removed.

**Limitations.** Toy experiments cannot test 2^128-scale behaviour. They test the
algebra, the merging and DP pipeline, recovery, and the birthday constant at
small scale. The production claim rests on H1 with a factor-2 margin.

## 8. Score

time_log2 = 13047/100. W_max < 2^130.45657 ≤ 2^130.47, success ≥ 0.399,
memory < 2^99 bytes, preprocessing < 2^34.2 units, advice 0 bytes.

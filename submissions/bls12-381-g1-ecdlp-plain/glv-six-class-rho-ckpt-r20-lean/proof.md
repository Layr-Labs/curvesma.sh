# Six-class parallel Pollard rho on BLS12-381 G1 (negation + GLV), checkpoint-only cycle escape, R = 2^20, operands of the reference addition charged once

**Claim.** time_log2 = 126.3045: the algorithm recovers x from T = [x]P0 with
probability >= 0.39, and on every run its total work is <= W_max = 2^126.3045
units (one unit = one batched affine G1 addition = C = 2,000 scalar AArch64
instructions). The work bound is unconditional (global work meter, section 3).
Memory < 2^72 bytes; preprocessing (everything before the first walk step:
jump table, report-table initialisation, the first 2^20 trail starts) <= 2^57.2
units, enforced by the meter and included in W; no advice. Success
rests on the random-mapping heuristic H1 (section 4), supported by toy runs.

The attack is the six-class (negation + GLV) distinguished-point rho of the
public packages `glv-six-class-rho-v11` / `-v12` of this track, re-stated in full
below, with three changes (each is argued in the section named):
(i) cycles of F are detected by the period-64 checkpoint only: the separate
2-cycle rule is removed (section 2; the ledger row shrinks, section 5);
(ii) the jump table has R = 2^20 entries instead of 2^14, which lowers the rates
of the random surcharges by 2^6 and 2^12 and lets the H1 floor and the meter
headroom be set tighter (sections 3, 4);
(iii) the per-step instruction ledger is re-itemised for (i), and the loads and stores of the
reference addition's own operands are not charged a second time: per-step maximum 146
instructions, budget 150 (section 5). The walk's class map F, the key, the canonical index,
the doubling escape, the work meter and heuristic H1 (with its 5% allowance)
are unchanged.

## 1. Algebra (all checked on BLS12-381 in `experiment.py`, part A)

p = field prime, r = subgroup order (log2 r = 254.85709), E: y^2 = x^3 + 4.

* beta = 0x1a0111ea...00000000aaac (primitive cube root of 1 in Fp),
  phi(x, y) = (beta x, y) = [lambda] on G1 with
  lambda = 0xac45a4010001a40200000000ffffffff (phi(P0) = [lambda]P0 checked).
* Units u_t = (-1)^s lambda^j (t = 3s + j) act by [u_t](x, y) = (beta^j x, (-1)^s y).
  No nonidentity G1 point is fixed by a unit != 1, so G1 \ {O} has exactly
  N = (r-1)/6 classes of size 6 (log2 N = 252.27213).
* **Key** K = |y| = min(y, p - y) is a complete class invariant: equal |y|
  gives y' = +-y and x'^3 = x^3, so x' = beta^j x.
* **Canonical member**: j minimises beta^j x mod p, s = [p - y < y];
  c = (beta^j x, |y|) = [u_t]P. The exact rule costs one Fp multiplication
  (beta x; beta^2 x = p - ((x + beta x) mod p)); the walk needs only the index
  t, which section 1.1 decides exactly without that multiplication.

### 1.1 Exact canonical index from a truncated product

Only j = argmin(v0, v1, v2) is needed, where v0 = x, v1 = beta x mod p,
v2 = beta^2 x mod p (distinct integers in [1, p)). In units of 2^320:

* **v1.** With G = floor(2^448 beta/p) (7 limbs), frac(xG/2^448) = v1/p - e,
  0 <= e < 2^-67. L' = (sum of the high halves of the six column-5 limb
  products x_i G_(5-i) and the low halves of the six column-6 products
  x_i G_(6-i)) mod 2^64: 12 mul/umulh. Everything dropped from xG mod 2^448
  (columns 0..4 and the column-5 low halves) is < 12 * 2^384, so the true bits
  384..447 lie in [L', L' + 12]. If L' < 16 or L' > 2^64 - 16, fall back.
  Otherwise V1 = floor(PH L' / 2^128), PH = floor(p/2^256), satisfies
  0 <= v1/2^320 - V1 < 3 (p/2^320 < 2^61, so the < 13/2^64 fraction error is
  < 2 units; floors < 1).
* **v0, v2.** X = top limb of x (error < 1). S = X + V1; if |S - PTOP| <= 16
  (PTOP = floor(p/2^320)) fall back; otherwise x + v1 < p iff S < PTOP and
  V2 = PTOP - S or 2 PTOP - S is within 5 of v2/2^320.
* **Decision.** If the smallest of X, V1, V2 is below both others by more than
  16, it is the exact argmin (errors < 5 each cannot reverse a gap > 10).
  Otherwise fall back to the exact rule: a separately metered action
  (0.25 units). The fallback windows cover < 2^7 units out of > 2^60, so the
  fallback rate is < 2^-52 per step if beta x mod p is equidistributed (this
  affects only the meter's early-halt probability, never correctness).

Part A compares this rule (limb-level emulation) with the exact rule on 30,000
random and 10,000 adversarial near-tie production inputs: 0 mismatches (offline:
400,000 random and 200,000 adversarial, 0 mismatches; observed maxima: bits
deficit 8 <= 12, V1 error 1.85 < 3, V2 error 2.55 < 5).

## 2. Walk

Parameters: R = 2^20 = 1,048,576; distinguished probability theta = 2^-64 (limb 1
of K is zero); M = 2^20 walks; L_max = 20 * 2^64; checkpoint period 64.

* **Table.** c_h, d_h uniform in [0, r) (fresh coins, rejection sampling),
  M_h = [c_h]P0 + [d_h]T (redraw if O), Tab[h][t] = [u_t^{-1}]M_h
  = (-1)^s phi^{-j}(M_h) with coefficients u_t^{-1}(c_h, d_h) mod r.
  6 R entries of 160 bytes: 1.0 * 10^9 bytes.
* **Step.** F(C) = class(c(C) + M_h), h = K mod R. The walk keeps any member P
  and sets P <- P + Tab[h][t]; since c + M_h = [u_t](P + Tab[h][t]) this is
  F. Each ordinary step is exactly one batched affine addition (reference
  operation) plus at most 150 instructions (section 5).
* **Trails.** Start S = [a]P0 + [b]T, a, b uniform from fresh coins, stored
  in full (redraw if S = O). At each class arrival, in order:

```
if K distinguished: report (K, a, b); end trail
if length >= L_max: abandon trail
if K_cp is set and fp(K) == K_cp:
     traverse the cycle once more (<= 64 steps, else abandon the trail),
     remembering its minimum-fp class; resume at the member of that class;
     ESC; clear history; continue
if length mod 64 == 0: K_cp <- fp(K)
P <- P + Tab[h][t]
```

  "Clear history" unsets K_cp.
  ESC is the **doubling escape** P <- [2]P (batched affine doubling). Doubling
  commutes with every unit, so it maps class C to class(2c(C)): a class map.
  fp(K) = limbs 0, 2, 3 of K (192 bits), stored as K_cp and compared as a
  192-bit integer (limb 3 most significant, then limb 2, then limb 0): the
  equality test and the "minimum" of the traversal use this one fixed total order on
  classes (distinct classes have distinct fingerprints except for the
  spurious-match failure event of section 4); full keys are never compared.
  `length` counts every operation (step, doubling, traversal step), so a
  trail has <= L_max + 65 operations. If an addition has equal x-coordinates
  (P = +-Tab), the trail is abandoned (< 2^-114 in total).
  There is no separate 2-cycle test and no pending escape: 2-cycles, fixed
  points and every other cycle of length <= 64 are found by the checkpoint
  test (v11/v12 additionally tested K against the class two arrivals back).
* **Detection delay.** Let a trail enter a cycle Z of F with |Z| <= 64. Within
  at most 63 further arrivals its length counter reaches a multiple of 64 and
  K_cp is set to a class of Z; the test precedes the refresh, so the trail
  returns to that class after exactly |Z| arrivals and the test fires. An older
  K_cp set outside Z (before the trail entered Z) is never revisited, because
  the part of the F-orbit before Z is traversed once. The traversal then needs
  |Z| steps and the escape is one doubling. So at most 63 + 2|Z| + 1 <= 192
  operations are spent from entering Z to the end of the escape (68 for a 2-cycle).
* **Merging lemma.** Let trails A and B arrive at the same class C, with any
  histories. Until one of them escapes, both apply F and visit the same
  classes. A trail escapes only when fp(K) == K_cp: the trail is back at a
  class it visited since its last history clear, which (F being a function and
  the history being cleared at every doubling) means it is on a cycle Z of F; the
  traversal confirms |Z| <= 64 (the traversal fails only on a spurious fingerprint
  match, a failure event of section 4; a cycle longer than 64 is never detected
  because K_cp is overwritten every 64 arrivals, and such cycles are covered by the
  trail cap and H1). Z contains no distinguished class (the trail
  would have ended there first). Every escape from Z leaves from the traversal's
  minimum-fp class and goes to
  E(Z) = class(2 c(min Z)) (min Z = the minimum-fp class of Z), which depends
  on Z only; the traversal visits every class of Z once, so the minimum does not
  depend on where in Z the trail was. If A escapes from
  Z, then B, following the same F path, also enters Z and, by the detection
  delay, escapes after at most 192 operations, to the same E(Z). So A and B are
  again at a common class after at most 192 extra operations of one of them, and by
  induction they reach the same distinguished class, unless one of them first
  hits L_max or an exceptional addition (failure events, section 4). A fixed
  point of F (|Z| = 1) is covered: the test fires at the arrival after the
  refresh, the traversal is one step, and the escape class is class(2c). History
  is cleared at every doubling and at every trail start (the valid word is 0 then), so fp(K) == K_cp always witnesses a
  real cycle of F.
* **No structural cycles through escapes.** A trail can return to a class c
  after m >= 1 doublings only if u 2^m c + (sum of unit and power-of-2 multiples
  of table points M_h met on the way) = c for a unit u. The order of 2 mod r is
  2^250.86 (checked in part A from the full factorisation of r - 1), so
  1 - u 2^m != 0 mod r for every m < 2^248: the relation cannot hold by
  cancellation of the table terms (unlike an additive escape multiplier, where
  negation pairs cancel identically); it pins c to one specific point. This
  excludes identities, not coincidences: that such returns occur no more
  often than ordinary self-collisions of a trail is part of heuristic H1
  (section 4). The
  experiment's long-trail test checks this with many escapes per trail
  (section 6).
* **Solving.** On a key match with another trail's report, re-walk both from
  their stored (a, b), tracking P = [A]P0 + [B]T with every value reduced
  mod r after every operation (each step adds the table entry's coefficients;
  a doubling doubles A and B). Canonical coefficients
  (e, f) = u_t (A, B). If the canonical points agree and f_1 != f_2:
  x = (e_2 - e_1)/(f_1 - f_2) mod r (sign checked on production data); verify
  [x]P0 = T and output x. T = O outputs 0.

## 3. Phases and the global work meter

* **Meter.** Every action (random words, sampling attempts, starts and redraws,
  steps, doublings, traversals, reports, hash probes, table build and hash
  initialisation, re-walks, solving, verification) is charged its price
  (section 5) before it runs; each price upper-bounds its instructions / C.
  If a charge would exceed W_max = 2^126.3045 the attack halts with failure.
  So W <= 2^126.3045 on every run; no loop needs its own limit. Preprocessing
  (all charges before the first walk step) has its own meter cap 2^57.2: if a
  preprocessing charge would exceed it, the attack halts with failure. The meter is a
  192-bit fixed-point counter; each of at most 2^20 processors (one per walk)
  reserves blocks of 2^40 units by an atomic update (unused remainders count
  as spent: < 2 * 2^20 * 2^40 = 2^61 units, included in the hard extras).
* **Phase 1** runs until a global operation counter (steps, doublings,
  traversal steps of phase-1 trails; re-walk operations are separate and not
  counted; at most 2^20 processors reserve blocks of 2^40, and every reserved
  block is completed before phase 1 stops) reaches S_ops, so at least S_ops
  and fewer than S_ops + 2^60 operations are performed. Ordinary steps are
  priced 1.075, traversal steps 1.095 (they may copy the remembered member),
  doublings 1.5. **Phase 2** finishes open trails, no new
  trails. Report table capacity D_cap = 2^63.25 and at most 8 processed
  matches (overflow = failure); these bound memory and re-walk work.
* **Batches.** Ready walks (next action a step or doubling) are split into
  batches of 2^15..2^16 per Montgomery inversion; a smaller batch (fewer than
  2^15 ready walks) pays its inversion as a separate 132-unit action. Each
  batch tests its product of denominators once (section 5).

## 4. Success probability

N' = N * R/(R-1), R = 2^20. T_c = the number of class arrivals of all trails (trail
starts, steps, doublings, traversal steps, and repeated arrivals inside
fruitless cycles: nothing is excluded) up to the first arrival at a class
already visited by another trail. This is exactly what the toy measures.

**H1 (heuristic, at the budget).** P(T_c <= S_ops) >= 1 - exp(-(S_ops/kappa)^2
/(2N')) with kappa = 1.05 (a 5% slowdown allowance relative to the
random-mapping law: toy class walks with R >= 2048 are -0.0% to +2.3% slower in
mean T_c, pooled +0.8%, while a no-automorphism control matches the law; section 6).
Also part of H1: partition
and distinguished bits of |y| are pseudo-random functions of the class; and,
**excluding the fruitless cycles of F of length <= 64** (2-cycles at rate
1/(6R) per step and the rarer short cycles below, which the walk detects and
escapes, and whose repeated arrivals are counted inside T_c), a trail makes a
*non-fruitless return* (an arrival at a class it visited earlier, either
before one of its escapes or along an F-cycle longer than 64) with probability
at most L_max^2/(2N) per trail, the random-mapping self-collision rate
(supported by the long-trail test of section 6; F-cycles of length 65..200
themselves need every partition index repeated, rate
sum_L sum_b S2(L, b) R^-(L-b) < 2^-206 per step, < 2^-137 per trail, and
longer ones are rarer still). Finally, **the canonical twist index of a new
class is uniform on the 6 units and independent of its partition index**
(used for the 2-cycle rate below).

S0 = sqrt(2 N' ln(1/0.6097)) = 2^126.12851; S_ops = 1.05 S0 = 2^126.19890
(H1 gives >= 0.3903; pure model 0.4205); U1 = 1.075 (S_ops + 2^60) < 2^126.303232.

**Phase 1 covers the budget.** Phase 1 performs at least S_ops operations
(hard counter) and each completed operation is one class arrival (an
operation aborted by an exceptional addition is a failure event), so phase 1
makes at least S_ops arrivals: if T_c <= S_ops, the first inter-trail collision happens
in phase 1. No statement about wasted visits is needed: fruitless-cycle
arrivals are counted inside T_c.

**Cost of the checkpoint-only detection.** The delay of section 2 (<= 68
operations per 2-cycle, <= 192 for any cycle) is spent on repeated arrivals,
which T_c counts and the phase-1 counter counts. At the 2-cycle rate 1/(6R) per
new class this is 68/(6R) < 2^-16.49 of the operations (the other cycles:
192 * 2^-37.99 < 2^-30), far below the 5% allowance of H1. In the toy (section 6)
the same delay costs at most 0.55% at R = 2048 and 0.07% at R = 16384, and is
inside the measured numbers.

**Doubling and traversal cost.** Only the surcharges of doublings
(1.5 - 1.075 = 0.425 each) and traversal steps (0.02 each) are not fixed by the
operation count. A doubling follows a detected cycle:
2-cycles need the next class to have the same partition index and canonical
twist -1, probability 1/(6R) per step; other fruitless cycles of length <= 64
have rate <= sum_L sum_b S2(L, b) R^-(L-b) = 2^-37.999 (L = 3..64; S2(L, b) =
number of partitions of the L positions into b blocks of size >= 2: every
partition index must occur >= 2 times, which includes identities such as one
index thrice with twists 1, lambda, lambda^2; the same sum is 2^-25.999 at
R = 2^14 and 3.638 * 10^-12 (rounded) at R = 2^20). So the expected number of
doublings is < 2^-22.584 (S_ops + 2^60) (1/(6R) + 2^-38 < 2^-22.584, 1/(6R) =
2^-22.585; linearity of
expectation over the pseudo-random partition and twist of each new class; no
independence across steps is used), and their expected surcharge is
< 0.425 * 2^-22.584 * 2^126.1989 < 2^102.39.
Each detected cycle triggers one traversal of |Z| steps (|Z| = 2 for a 2-cycle,
<= 64 otherwise): the expected number of traversal steps is
< (2/(6R) + 64 * 2^-37.99)(S_ops + 2^60) < 2^104.62, their expected surcharge
< 0.02 * 2^104.62 < 2^98.98. (The large table, R = 2^20, is what makes these
surcharges small.) The toy cycle counts are reported in section 6.

**Failure events** (union bound with unconditional probabilities; if
T_c <= S_ops, by the merging lemma the two trails reach the same report and
the attack succeeds unless one of these occurs): trail cap hit first (Gamma(2) tail, < 1e-7);
a non-fruitless return inside one of the colliding trails (H1: L_max^2/(2N)
per trail, < 2^-53 in total; detected fruitless cycles are not failures: they
are escaped and the merging lemma covers them);
report overflow (under H1 the distinguished bits of distinct classes are
independent 2^-64 coins, so new-class reports are Binomial(n, 2^-64) with
n <= S_ops + 2^60 + 2^88.4 + starts, mean < 2^62.22; D_cap = 2^63.25 is more than 2.04
times the mean; Chernoff: < exp(-2^60.9); repeated reports of a class are
matches, at most 8 processed);
f_1 = f_2 (b uniform given the start point: 1/r); exceptional additions
(< 2^-114); 8 processed matches all failing (the attack stops at the first
verified x, so this needs 8 spurious or degenerate matches: < 2^-60); spurious 192-bit fingerprint match
(union bound over < 2^127.5 fingerprint comparisons: < 2^-64); meter halting
early (expected surcharges 2^102.39 + 2^98.98 plus other random charges
< 2^74.8, total < 2^102.52, against headroom > 2^116.15: < 2^-13.63 < 0.00008 by Markov,
section 5); preprocessing cap (< 2^-32). Total of the failure events other than
the meter: < 1.1 * 10^-7.
**Success >= 0.3903 - 0.00008 - 0.00000011 > 0.39021 > 0.39.** (The margin also
absorbs a modelled doubling and traversal rate up to 3 times the 1/(6R) value:
the meter bound would then be < 0.00025, success >= 0.39005.)

## 5. Cost

**Per ordinary step** (instructions besides the reference addition; blst
Fp multiplication = 383):

| Work | Instr. |
| --- | ---: |
| 1.1 v1: 12 mul/umulh, 12 add, 22 loads | 46 |
| 1.1 v1 range test, V1 (3 mul, 2 add); X, S, wrap test, V2 | 20 |
| 1.1 decision: differences, gap tests, argmin, fallback branch | 14 |
| sign, p - y, select limbs of K | 17 |
| h, distinguished test, twist index, table address | 7 |
| checkpoint test (`ldp` limbs 0 and 2 of K_cp, `ldp` limb 3 and the valid word, `cmp`, 3 `ccmp`, `b.eq`: 7), refresh test (`tst` length, `b.ne`: 2), refresh (`stp`, `stp`, constant 1 for the valid word: 3, counted in every step although it happens every 64th) | 12 |
| length, cap test, meter additions, loop | 14 |
| batch bookkeeping: walk-record and table-record pointer arithmetic and loop control of the three batch passes (forward products, back-substitution, finalise) | 8 |
| batch inversion share (< 2^18 per >= 2^15 elements) | 8 |
| **total, a per-step maximum (budget 150 = 0.075 units)** | **146** |

**Operands of the reference addition are not charged twice.** The 196-instruction
ledger of v11/v12 (and the 174-instruction sibling package
`glv-six-class-rho-ckpt-r20`) contained a row of 36 instructions for the loads of the
table point (12 limbs) and of the walker's point (12 limbs) and the stores of the new
point (12 limbs). Those are the inputs and the outputs of the reference addition
itself. COST_MODEL.md defines the reference operation as the batched incomplete affine
addition under the 6M convention with S = M and A = 0 and says "Do not charge both a
reference operation and its internal instructions". Its price was calibrated with
Callgrind on `blst_fp_mul(&a, &a, &b)`, a call through pointers to operands in memory,
so the 383 instructions per multiplication already contain the loads of both operands
and the store of the result, and the additions and subtractions of the formula (the
differences dx = x1 - x2 and dy = y1 - y2, x3 = lambda^2 - x1 - x2 and
y3 = lambda (x1 - x3) - y1), with the data they read, are free by A = 0. So reading
x1, y1, x2, y2 and writing x3, y3 are inside the one unit. Everything the walk does
around the addition is charged in the other rows: the loads of the new point's limbs
for the canonical index and the key (rows 1.1 and sign/key), the table-row address, the
per-walk checkpoint and length fields, and the batch bookkeeping row (each batch of
>= 2^15 ready walks is processed in three passes, forward prefix products,
back-substitution and finalise; each pass advances a walk-record pointer and a
table-record pointer and tests its loop counter, about 2 instructions per pass and walk,
plus 2 for the per-walk table-row base: 8). This is the single accounting judgement of
this package. Read as pure arithmetic, the reference operation would leave the 36
instructions to be charged, which gives step price 1.0875 and the bound 2^126.321 of the
sibling package.

Checkpoint row, in detail. The walk record holds K_cp as limbs 0, 2, 3 of the
fingerprint and a valid word (1 after a refresh; 0 after every escape and at every
trail start). Every arrival does the full test unconditionally, with no data-dependent
slow path: `ldp` x_a, x_b (limbs 0, 2 of K_cp); `ldp` x_c, x_v (limb 3, valid word);
`cmp` k0, x_a; `ccmp` k2, x_b, #0, eq; `ccmp` k3, x_c, #0, eq; `ccmp` x_v, #1, #0, eq;
`b.eq` to the traversal (7 instructions; the flags are equal iff all three
limbs match and the walk record is valid). The only branch that is taken is a real
match fp(K) == K_cp, which is the traversal (priced by the traversal step price and
its surcharge). `tst` of the length counter's low limb and `b.ne` (2 instructions)
decide whether a refresh is due; the refresh is `stp` k0, k2, [cp] and `stp` k3, one,
[cp + 16] (2 instructions plus 1 to materialise the constant 1 for the valid word; it
only happens every 64th operation, but all 3 are charged to every step). The v11/v12
ledger spent 26 instructions on the 2-cycle test, the two history stores per step and the
checkpoint; with the 2-cycle test removed, only the above remains, so the itemised
maximum is 46 + 20 + 14 + 17 + 7 + 12 + 14 + 8 + 8 = 146.

Zero-denominator test: the Montgomery batch product of the denominators dx_i =
x(P_i) - x(Tab_i) is formed anyway, and its output is a fully reduced element of
[0, p) (blst's `mul_mont_384` ends with a conditional subtraction), so the product
is zero iff all six limbs are zero; a field product is zero iff some factor is zero.
Per batch of >= 2^15 elements the code ORs the six limbs (5 `orr`) and does one
`cbz`: 6 instructions per batch, 0.0002 per step, inside the batch row. A nonzero
product certifies that every dx_i != 0. If it is zero, the recovery is a separately
metered action charged 2^17 units before it runs: the batch is scanned, the walks with
dx_i = 0 (P = +-Tab, probability < 2^-114 in total) are abandoned and the products
recomputed (3 multiplications per walk, below 2^15.3 units per batch of 2^16).
(v11/v12 charged 8 instructions per step for a per-element test.)

The rare fallback of section 1.1 is not in this price (separate action).

**Prices** (units): ordinary step 1.075; traversal step 1.095 (a traversal step
is an F-step without the checkpoint row, plus <= 24 instructions for copying the
remembered member and <= 6 for the lexicographic comparison of its 192-bit fingerprint
with the running minimum (2 `ldp`, `subs`, 2 `sbcs`, branch): at most 150 + 24 + 6 =
180 instructions = 1.09 <= 1.095); doubling 1.5 (7M + 150
instr. = 2,831 instr. < 1.42); section-1.1 fallback 0.25 (<= 500 instr.); random 64-bit word (RNDR) 0.005; sampling attempt 0.025;
trail start [a]P0 + [b]T (two Jacobian scalar multiplications + inversion,
about 4.8e6 instr.) 2^12; report 0.5; hash probe 0.015; hash initialisation
2^57.1 (once); table 2^32.1 (once: 2^20 entries, each a trail-start-priced
M_h = [c_h]P0 + [d_h]T plus 5 twists with 2 Fp multiplications and 10 modular
multiplications for the coefficients, < 2^12 + 8 units, plus its sampling); phase-2 or re-walk
operation (unbatched, incl. inversion and coefficient updates) 133; small-batch
inversion 132; solving 2^9; verification 2^13.

**Budget.** Phase-1 operations: 1.075 (S_ops + 2^60) < 2^126.303232, plus the
doubling and traversal surcharges (random, expectation < 2^102.39 + 2^98.98).
Hard extras: meter-block remainders < 2^61, phase 2
(2^20 (L_max + 65) * 133 = 2^95.4), re-walks (8 * 2 * ((L_max + 65) * 133 +
2^12) = 2^79.4), reports (D_cap * 0.5 = 2^62.3), hash init 2^57.1, table and
verification < 2^32.2: together < 2^95.5. Random extras (expected): starts
< 2^62.4 * 2^12 = 2^74.4, words/attempts/probes < 2^63, small batches < 2^8,
fallbacks < 2^72.3, zero-batch recoveries < 2^-96: together < 2^74.8. Headroom
2^126.3045 - 2^126.303232 - 2^95.5 > 2^116.15, so the meter halts early with probability
< (2^102.39 + 2^98.98 + 2^74.8)/2^116.15 < 2^-13.63 < 0.00008
(Markov).

**Preprocessing:** jump table (priced 2^32.1, plus sampling),
report-table initialisation 2^57.1, first 2^20 starts (2^32, plus sampling):
expected < 2^57.11, cap 2^57.2 (the expected random sampling part is < 2^21,
far below the 2^53 gap, so the cap halts a run with probability < 2^-32,
Markov); all of it is inside W.

**W <= 2^126.3045 on every run.** Memory: D_cap records x 104 bytes at load 1/2
= 2^70.95, plus table (2^20 x 6 x 160 bytes = 1.0 GB) and walk states (2^20 x 512 bytes): < 2^72.

## 6. Experiments and limitations

`experiment.py` (stdlib). **Part A** (production BLS12-381; all true): curve,
order, beta/lambda pairing, orbit and unit action, |y| invariance, unique
canonical member, twisted-step class identity, doubling class map, order of 2
mod r, table
entries, coefficient tracking on production steps, and recovery of a planted
x from two production trails forced to merge. **Part B** (toy y^2 = x^3 + 4
over F_q, q = 58493047 or 861101287): the per-trail algorithm of section 2
(exact canonical rule, which section 1.1 always reproduces, checkpoint-only cycle
escape, caps, restarts, re-walk solving). Default mode
reports T_c; `"pipeline": true` also applies the scaled production rules
(phase-1 budget 1.05 S0 class arrivals, phase 2, report capacity, 8-match
cap; the toy capacity is about 2.7 times the mean instead of 2.07, i.e.
looser) and reports success. Trails run sequentially (same visit distribution
under H1); [a]G + [b]T is evaluated as [a + bx]G for speed (same point); toy
keys are below 2^64, so the fingerprint (limbs 0, 2, 3) of a toy key is the key
itself and the toy's equality tests and minima follow the production rule
exactly. A fresh jump table is drawn for every trial.

Results (our runs of this `experiment.py`, 14-second runs, seeds 3000-3714; every
recovered x verified; S_ops = 1.05 sqrt(2 N' ln(1/0.6097)) for each toy; T_c counts all
arrivals; the number of trials per run depends on machine speed):

| mode / classes | toy, R | runs / trials | mean T_c vs model | P[T_c <= S_ops] or success |
| --- | --- | ---: | --- | ---: |
| T_c, 6 | q26, 2048 | 6 / 7464 | +1.2% +- 0.6% | 0.4116 +- 0.0057 |
| T_c, 6 | q26, 16384 | 12 / 2285 | +2.3% +- 1.1% | 0.4031 +- 0.0103 |
| T_c, 6 | q30, 2048 | 12 / 7908 | +0.2% +- 0.6% | 0.4164 +- 0.0055 |
| T_c, 6 | q30, 16384 | 16 / 2331 | -0.0% +- 1.1% | 0.4282 +- 0.0102 |
| T_c, 1 (control, no automorphism) | q26, 2048 | 3 / 2879 | +0.0% +- 1.0% | 0.4210 +- 0.0092 |
| pipeline success, 6 | q26 2048 / q26 16384 / q30 16384 | 6 / 8 / 8 runs; 7934 / 1539 / 1209 trials | - | 0.4201 / 0.4418 / 0.4285 |

Pooled class-6 P[T_c <= S_ops] = 0.4145 +- 0.0035 over 19,988 trials (H1 floor
0.3903; worst configuration 0.4031 +- 0.0103); pooled pipeline success 0.4242 +-
0.0048 over 10,682 trials (claim 0.39). The declared initial run (q26, R = 2048,
seed 1, 14 s) printed all Part A checks true and pipeline success 0.404 over 1,285
trials. The cycle delay of section 2 is up to 68/(6R) of the operations, i.e. 4.4%
at R = 256, which is not the production regime (2^-16.49); a run at that size shows
the delay at work: q26, R = 256, 8 runs (seeds 3800-3807), 28,629 trials, mean T_c
+4.1% +- 0.3% over the model and P[T_c <= S_ops] = 0.3932 +- 0.0029 (v12, which had
the 2-cycle test, measured +2.1% at R = 256), so the allowance of 5% still holds
but with little room; smaller tables are not used. At toy size R = 16384 exceeds
sqrt(N), so those rows test the walk mechanics at a large table rather than the
R << sqrt(N) regime; the production R = 2^20 itself is not run (the rows at
R = 2048 and 16384 show no visible R dependence of the slowdown).

Cycle events per step against 1/(6R), the T_c-mode and pipeline runs pooled (all
cycle types): 1.07 and 0.96 (R = 2048, q26 and q30), 1.89 (q26, R = 16384; 146
events, varying from run to run between 0.65 and 3.15), 1.05 (q30, R = 16384) and
0.96 (q26, R = 256). Besides the table-driven cycles, a toy trail closes a cycle of
length <= 64 whenever its next class repeats one of its previous <= 64 classes, an
ordinary self-collision with probability about (mean trail length 32)/N per step;
relative to 1/(6R) this adds about 32 * 6R/N: 1.97 for q26 at R = 16384 (N = 3,249,405,
observed 1.89), 1.12 for q26 at R = 2048 (observed 1.07), 1.07 for q30 at R = 16384
(observed 1.05). In production this term is at most 64/N = 2^-246 per step, far
below 1/(6R). This is a back-of-the-envelope match to the toy numbers, not a separate
experiment; the cycle rate enters only the meter's early-halt probability, which the
success margin absorbs up to 3 times the modelled rate (section 4).

**Long-trail trap test** (`long_trails`; q30, N = 2^25.5, walks without
distinguished points, 150 trails of 2,048 operations each, 15-second runs with seeds
21 (R = 4), 22 (R = 8), 23 (R = 4 control), 24 (R = 8 control); repeats inside a
fruitless cycle of length <= 64 (same escape epoch, <= 64 arrivals after the
previous visit) are counted separately; every other repeat is a non-fruitless return).
With checkpoint-only detection at R = 4 or 8 most operations are cycle
repeats, so the self-collision expectation is computed for the new-class arrivals
(operations minus repeats), as the experiment reports it:

| R | escape | escapes per trail | new-class arrivals per trail | trails with non-fruitless returns (expected for those arrivals) |
| ---: | --- | ---: | ---: | ---: |
| 4 | doubling (the attack) | 29.5 | 654 | 0 / 150 (0.67) |
| 8 | doubling (the attack) | 22.8 | 1097 | 2 / 150 (1.88) |
| 4 | additive M_esc (control) | 30.8 | 370 | 72 / 150 (0.21) |
| 8 | additive M_esc (control) | 25.9 | 903 | 21 / 150 (1.27) |

Doubling escapes behave like ordinary self-collisions even with many escapes
per trail; the additive escape traps trails structurally, as expected.

Limitations: H1 is a heuristic (random mapping with Brent-Pollard correction;
van Oorschot-Wiener; automorphism classes as in Wiener-Zuccherato, GLV,
Duursma-Gaudry-Morain; fruitless-cycle escapes as in Bos-Kleinjung-Lenstra and
Bernstein-Lange-Schwabe), supported only at toy scale; instruction counts are
analytic (a row-by-row itemisation, not a measured binary); the production
table size R = 2^20 is not exercised by the toy; generic attack (no TNFS).

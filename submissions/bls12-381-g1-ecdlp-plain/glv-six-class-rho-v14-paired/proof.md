# Six-class parallel Pollard rho on BLS12-381 G1 (negation + GLV)

**Claim.** time_log2 = 126.2482: the algorithm recovers x from T = [x]P0 with
probability >= 0.39, and on every run its total work is <= W_max = 2^126.2482
units (one unit = one batched affine G1 addition = C = 2,000 scalar AArch64
instructions). The work bound is unconditional (global work meter, section 3).
Memory < 2^72 bytes; preprocessing (everything before the first walk step:
jump table, report-table initialisation, the first 2^20 trail starts) <= 2^57.2
units, enforced by the meter and included in W; no advice. Success
rests on the random-mapping heuristic H1 (section 4), supported by toy runs.

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

Parameters: R = 2^20; distinguished probability theta = 2^-64 (limb 1 of K is
zero); M = 2^20 walks; L_max = 20 * 2^64; checkpoint period 64.

* **Table.** c_h, d_h uniform in [0, r) (fresh coins, rejection sampling),
  M_h = [c_h]P0 + [d_h]T (redraw if O), Tab[h][t] = [u_t^{-1}]M_h
  = (-1)^s phi^{-j}(M_h) with coefficients u_t^{-1}(c_h, d_h) mod r.
* **Step.** F(C) = class(c(C) + M_h), h = K mod R. The walk keeps any member P
  and sets P <- P + Tab[h][t]; since c + M_h = [u_t](P + Tab[h][t]) this is
  F. Each ordinary step is exactly one batched affine addition (reference
  operation) plus at most 128 instructions (section 5).
* **Trails.** Start S = [a]P0 + [b]T, a, b uniform from fresh coins, stored
  in full (redraw if S = O). At each class arrival, in order:

```
if K distinguished: report (K, a, b); end trail
if length >= L_max: abandon trail
if fp(K) == K_cp: traverse the cycle once more (<= 64 steps, else abandon),
     resume at the member of its minimum-fp class remembered during the
     traversal (no extra operation), ESC; clear history; continue
if length mod 64 == 0: K_cp <- fp(K)
P <- P + Tab[h][t]
```

  Cycles are detected only by the checkpoint (the separate 2-cycle test of
  earlier versions is dropped: at R = 2^20 cycles are rare, see section 4).
  "Clear history" invalidates K_cp (a valid flag set to 0; also 0 at every
  trail start).
  ESC is the **doubling escape** P <- [2]P (batched affine doubling). Doubling
  commutes with every unit, so it maps class C to class(2c(C)): a class map.
  fp(K) = limbs 0, 2, 3 of K (192 bits), stored as K_cp and compared as a
  192-bit integer (limb 3 most significant, then limb 2, then limb 0): the
  checkpoint test and the traversal minimum use this one fixed total order
  on classes (distinct classes have distinct fingerprints except for the
  spurious-match failure event of section 4); full keys are never compared.
  `length` counts every operation (step,
  doubling, traversal step), so a trail has <= L_max + 65 operations. If an addition has equal
  x-coordinates (P = +-Tab), the trail is abandoned (< 2^-114 in total).
* **Merging lemma.** Let trails A and B arrive at the same class C, with any
  histories. Until one of them escapes, both apply F and visit the same
  classes. A trail escapes only after K = K_cp, i.e. a repeat within 64
  ordinary steps (K_cp is set at an arrival and the test precedes the next
  refresh 64 operations later), so only when it is on a cycle Z of F of
  length <= 64 (a spurious fingerprint match is a failure event, section 4);
  Z then contains no distinguished class (the trail would have ended there
  first). **Detection delay:** a trail on Z reaches a length that is a
  multiple of 64 within 63 arrivals, which sets K_cp to a class of Z; it is
  back at that class |Z| arrivals later and the test fires (a K_cp set before
  the trail entered Z is never matched, barring a non-fruitless return). The
  traversal then visits all of Z once and the trail doubles the remembered
  member of min Z (the minimum-fp class of Z), going to
  E(Z) = class(2 c(min Z)), which depends on Z only. The detour takes at most
  63 + 2|Z| + 1 <= 192 operations; fixed points of F (|Z| = 1) take the same
  path. If A escapes from Z, then B, following the same F path, also enters Z,
  detects it within 128 arrivals and escapes to the same E(Z). So A and B are
  again at a common class after at most 192 extra operations of one of them,
  and by induction they reach the same distinguished class, unless one of
  them first hits L_max or an exceptional addition (failure events, section
  4). K_cp is invalidated at every doubling and is invalid at starts, so
  K = K_cp always witnesses a real cycle of F.
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
  If a charge would exceed W_max = 2^126.2482 the attack halts with failure.
  So W <= 2^126.2482 on every run; no loop needs its own limit. Preprocessing
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
  priced 1.064, traversal steps 1.073 (they may copy and compare the
  remembered member), doublings 1.5165. **Phase 2** finishes open trails, no new
  trails. Report table capacity D_cap = 2^63.25 and at most 8 processed
  matches (overflow = failure); these bound memory and re-walk work.
* **Batches.** Ready walks (next action a step or doubling) are split into
  batches of 2^18..2^19 per Montgomery inversion (possible whenever at least
  2^18 of the 2^20 walks are ready); a smaller batch pays its inversion as a
  separate 132-unit action. The meter is charged once per batch
  (batch size x price, one 192-bit addition), and the zero-denominator test is
  done once per batch on the batch product of the denominators (zero iff some
  denominator is zero; then a 2^17-unit recovery action scans the batch and
  abandons the offending walks; probability < 2^-114 in total).

## 4. Success probability

N' = N * R/(R-1). T_c = the number of class arrivals of all trails (trail
starts, steps, doublings, traversal steps, and repeated arrivals inside
fruitless cycles: nothing is excluded) up to the first arrival at a class
already visited by another trail. This is exactly what the toy measures.

**H1 (heuristic, at the budget).** P(T_c <= S_ops) >= 1 - exp(-(S_ops/kappa)^2
/(2N')) with kappa = 1.02 (a 2% slowdown allowance relative to the
random-mapping law). Toy slowdowns are a harness effect that scales with
(mean trail length)/sqrt(N): a trail's own classes are not useful collision
targets, which delays T_c by about one trail. A trail-length sweep on three
toy curves (section 6) fits slowdown = a + b (trail/sqrt(N)) with
a = +0.09% +- 0.23%; in production trail/sqrt(N) = 2^64/2^126.1 = 2^-62, so the
2% allowance is about 8 standard errors above the extrapolated value. (Earlier
fixed-trail runs gave -0.2% to +2.5% at trail ratios 0.9%-1.8%; a
no-automorphism control is within 0.6% of the law. The trail-length
explanation and a first sweep are due to submission 11084040, section 6.) Also part of H1: partition
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

S0 = sqrt(2 N' ln(1/0.6097)) = 2^126.12851; S_ops = 1.02 S0 = 2^126.15708
(H1 gives >= 0.3903; pure model 0.4024); 1.064 (S_ops + 2^60) < 2^126.246574.

**Phase 1 covers the budget.** Phase 1 performs at least S_ops operations
(hard counter) and each completed operation is one class arrival (an
operation aborted by an exceptional addition is a failure event), so phase 1
makes at least S_ops arrivals: if T_c <= S_ops, the first inter-trail collision happens
in phase 1. No statement about wasted visits is needed: fruitless-cycle
arrivals are counted inside T_c.

**Doubling and traversal cost.** Only the surcharges of doublings
(1.5165 - 1.064 = 0.4525 each) and traversal steps (0.009 each) are not fixed
by the operation count. A doubling follows a detected cycle:
2-cycles need the next class to have the same partition index and canonical
twist -1, probability 1/(6R) per step; other fruitless cycles of length <= 64
have rate <= sum_L sum_b S2(L, b) R^-(L-b) = 2^-37.999 (every partition index
must occur >= 2 times, which includes identities such as one index thrice with
twists 1, lambda, lambda^2). So the expected number of doublings is
< 2^-22.58 (S_ops + 2^60) (1/(6R) = 2^-22.585; linearity of expectation over the pseudo-random
partition and twist of each new class; no independence across steps is
used), and their expected surcharge is < 0.4525 * 2^-22.58 * 2^126.157 =
2^102.43. The expected number of traversal
steps is < (2/(6R) + 64 * 2^-37.99)(S_ops + 2^60) < 2^104.58; at 0.009 units
each the expected traversal surcharge is < 2^97.79. The checkpoint-only detection delay (at most
63 + |Z| ordinary steps per cycle) is ordinary work counted in S_ops and its
arrivals are counted in T_c; at R = 2^20 it is < 2^-16.4 of all operations.
(The large table, R = 2^20, is chosen to make cycles rare.)

**Failure events** (union bound with unconditional probabilities; if
T_c <= S_ops, by the merging lemma the two trails reach the same report and
the attack succeeds unless one of these occurs): trail cap hit first (Gamma(2) tail, < 1e-7);
a non-fruitless return inside one of the colliding trails (H1: L_max^2/(2N)
per trail, < 2^-53 in total; detected fruitless cycles are not failures: they
are escaped and the merging lemma covers them);
report overflow (under H1 the distinguished bits of distinct classes are
independent 2^-64 coins, so new-class reports are Binomial(n, 2^-64) with
n <= S_ops + 2^60 + 2^88.4 + starts, mean < 2^62.16; D_cap = 2^63.25 is 2.13
times the mean; Chernoff: < exp(-2^60.9); repeated reports of a class are
matches, at most 8 processed);
f_1 = f_2 (b uniform given the start point: 1/r); exceptional additions
(< 2^-114); 8 processed matches all failing (the attack stops at the first
verified x, so this needs 8 spurious or degenerate matches: < 2^-60); spurious 192-bit fingerprint match
(union bound over < 2^127.5 fingerprint comparisons: < 2^-64); meter halting
early (expected surcharges 2^102.43 + 2^97.79 plus other random charges
< 2^74.8, together < 2^102.49, against headroom > 2^116.45: < 2^-13.96 < 0.000063
by Markov, section 5); preprocessing cap (< 2^-32).
**Success >= 0.3903 - 0.000063 - 1.1e-7 > 0.39023 > 0.39.** A doubling and
traversal rate up to 3 times the 1/(6R) value gives a meter bound < 0.00019
and success >= 0.3901 (section 5).

## 5. Cost

**Per ordinary step** (instructions besides the reference addition; blst
Fp multiplication = 383):

| Work | Instr. |
| --- | ---: |
| 1.1 v1: 3 `ldp` (limbs of x), 3 `ldp` + `ldr` (G_0..G_6), 6 `umulh` (column 5: x5G0, x4G1, x3G2, x2G3, x1G4, x0G5), 5 `add` (the later column-5 high halves), 6 `madd` (column 6: x5G1, x4G2, x3G3, x2G4, x1G5, x0G6, each adding a low half into L') | 24 |
| 1.1 v1 range test (`sub`, `cmn`, `b.hs`), V1 (`umulh`, `mul`, `umulh`, `cmn`, `adc`: 3 mul, 2 add); X (limb 5 of x, already in a register), S (`add`), wrap test and V2 (`subs`, `csel`, `add`, `cmp`, `b.ls`, `sub`) | 15 |
| 1.1 decision: differences, gap tests, argmin, fallback branch (4 `cmp`, 6 `csel`, `cset`, `mov`, `sub`, `b.ls`) | 14 |
| sign, p - y, select limbs of K: 3 `ldp` (limbs of y), p - y (`subs`, 4 `sbcs`, `sbc`), exact comparison of p - y with y over all six limbs (`cmp`, 5 `sbcs`), `cset` (s), 4 `csel` (limbs 0-3 of K) | 20 |
| h, distinguished test, twist index, table address (`and`, `cbz`, 5 `add` with shifted operand, `str` of the entry address) | 8 |
| loads of table point and state, store of state (`ldr` of the entry address, 6 `ldp` table point, 6 `ldp` state, 6 `stp` state) | 19 |
| checkpoint test (`ldp`, `ldp`, `cmp`, 3 `ccmp`, `b.eq`: 7), refresh test (`tst` length, `b.ne`: 2), refresh (`stp`, `stp`, constant 1 for the valid word: 3, counted in every step although it happens every 64th) | 12 |
| length increment (69-bit), cap test, loop (`ldp`, `cmp`, `b.hs`, `adds`, `adc`, `stp`; `add`, `cmp`, `b.ne`) | 9 |
| zero-denominator test on the batch product (<= 9 instructions per batch of >= 2^18, one per step) | 1 |
| per-batch meter charge (<= 10 instructions per batch of >= 2^18, one per step) | 1 |
| batch inversion share (one Fermat inversion per batch: 380 squarings + 228 multiplications = 608 x 383 = 232,864, exponent scan < 2,000, the 4 `ldp` of p and PH once per batch: < 2^18 per >= 2^18 elements) | 1 |
| **total, a per-step maximum** | **124** |
| **budget with 4 instructions of stated slack (= 0.064 units)** | **128** |

The rare fallback of section 1.1 is not in this price (separate action).

Register use. Rows 1-5, 7 and 8 run as one leaf loop over the walk records of a
batch, after the batched addition has stored the new points; the loop calls nothing.
For the whole batch x10-x15 hold the six limbs of p, x16 and x17 hold PH0 and
PTOP (= PH1), and x19-x22 hold the walk-record pointer, the end of the batch, the
table base and the base of the constants; these are loaded once per batch (the
4 `ldp` of row 11). G is reloaded every step (row 1), and so are the limbs of x
and y (rows 1 and 4). x0-x9 and x25-x28 are temporaries; at most 24 registers are
live at once (in row 4), 26 are touched, and x18, x29 and x30 are not used. A taken
branch (fallback, report, traversal, abandonment, block reservation, end of an
operation block) leaves the loop for a separately priced action, which saves and
restores what it clobbers inside its own price. Walk record (512 bytes): x at offset
0, y at 48, the next table-entry address at 96, K_cp at 104 (limbs 0, 2) and 120
(limb 3, valid word), length (two limbs) at 136. Table entry 6h + t is at
table + 160 (6h + t). The length load and the cap test of row 8 run before row 7,
whose refresh test reads the length's low limb.

Rows 1-8 as instruction sequences (every operand load and store is listed; `ldp`
and `stp` move two adjacent limbs of a field element in one instruction):

```
row 1 (24)  ldp x0,x1,[x19,#0]; ldp x2,x3,[x19,#16]; ldp x4,x5,[x19,#32]
            ldp x6,x7,[x22,#0]; umulh x9,x5,x6; umulh x8,x4,x7; add x9,x9,x8; madd x9,x5,x7,x9
            ldp x6,x7,[x22,#16]; umulh x8,x3,x6; add x9,x9,x8; madd x9,x4,x6,x9
                                 umulh x8,x2,x7; add x9,x9,x8; madd x9,x3,x7,x9
            ldp x6,x7,[x22,#32]; umulh x8,x1,x6; add x9,x9,x8; madd x9,x2,x6,x9
                                 umulh x8,x0,x7; add x9,x9,x8; madd x9,x1,x7,x9
            ldr x6,[x22,#48]; madd x9,x0,x6,x9                                  (x9 = L')
row 2 (15)  sub x8,x9,#16; cmn x8,#31; b.hs fallback           (L' < 16 or L' > 2^64 - 16)
            umulh x0,x16,x9; mul x1,x17,x9; umulh x2,x17,x9; cmn x1,x0; adc x1,x2,xzr   (x1 = V1)
            add x2,x5,x1                                         (x5 = X, x2 = S)
            subs x3,x2,x17; csel x4,x17,xzr,hs; add x6,x3,#16; cmp x6,#32; b.ls fallback
            sub x2,x4,x3                                         (x2 = V2)
row 3 (14)  cmp x1,x5; csel x3,x1,x5,lo; csel x4,x5,x1,lo; cset x6,lo
            cmp x2,x4; csel x7,x2,x4,lo; cmp x2,x3; csel x7,x3,x7,lo; csel x3,x2,x3,lo
            mov x8,#2; csel x6,x8,x6,lo; sub x7,x7,x3; cmp x7,#16; b.ls fallback   (x6 = j)
row 4 (20)  ldp x0,x1,[x19,#48]; ldp x2,x3,[x19,#64]; ldp x4,x5,[x19,#80]
            subs x7,x10,x0; sbcs x8,x11,x1; sbcs x9,x12,x2; sbcs x25,x13,x3; sbcs x26,x14,x4; sbc x27,x15,x5
            cmp x7,x0; sbcs xzr,x8,x1; sbcs xzr,x9,x2; sbcs xzr,x25,x3; sbcs xzr,x26,x4; sbcs xzr,x27,x5
            cset x28,lo                                          (x28 = s = [p - y < y])
            csel x0,x7,x0,lo; csel x1,x8,x1,lo; csel x2,x9,x2,lo; csel x3,x25,x3,lo   (limbs 0-3 of K)
row 5 (8)   and x4,x0,#0xfffff; cbz x1,report
            add x5,x28,x4,lsl #1; add x5,x5,x5,lsl #1; add x5,x5,x6      (6h + 3s + j = 6h + t)
            add x5,x5,x5,lsl #2; add x5,x21,x5,lsl #5; str x5,[x19,#96]  (table + 160 (6h + t))
row 6 (19)  ldr x5,[x19,#96]                                  (entry address)
            6 ldp: table point x, y (12 limbs)         6 ldp: state point x, y (12 limbs)
            6 stp: new state point x, y (12 limbs)
row 7 (12)  ldp x7,x8,[x19,#104]; ldp x9,x25,[x19,#120]
            cmp x0,x7; ccmp x2,x8,#0,eq; ccmp x3,x9,#0,eq; ccmp x25,#1,#0,eq; b.eq traversal
            tst x27,#63; b.ne norefresh; stp x0,x2,[x19,#104]; mov x26,#1; stp x3,x26,[x19,#120]
row 8 (9)   ldp x27,x28,[x19,#136]; cmp x28,#20; b.hs abandon    (length >= 20 * 2^64 = L_max)
            adds x27,x27,#1; adc x28,x28,xzr; stp x27,x28,[x19,#136]
            add x19,x19,#512; cmp x19,x20; b.ne loop
```

Row 1's 12 products are all present (six `umulh` for the column-5 high halves,
six `madd` for the column-6 low halves), and its sum of the 12 terms uses 11
additions (five `add` plus six accumulated by `madd`), one fewer than a sum of
12 terms needs to start from zero. The key row keeps the exact comparison:
p - y is formed over all six limbs, and the borrow of (p - y) - y over all six
limbs is the sign, so no bit test stands in for y > p/2. Only limbs 0-3 of K are
selected, because a step uses only those (h and the fingerprint from limbs 0, 2,
3; the distinguished test from limb 1); the report action recomputes all six
limbs inside its 0.5-unit price. Row 6 keeps the scope of the earlier 36-transfer
row (one load of each limb of the table point and of the state, one store of each
limb of the new state) and pairs adjacent limbs; it is the reference addition's
own operand traffic, so it is shown as counts of `ldp`/`stp` rather than as a
register schedule that would collide with the addition's six multiplications.
Its 36 scalar transfers become 19 instructions.

**Prices** (units): ordinary step 1 + (124 + 4)/2000 = 1.064; traversal step
1.073 (128 + 12 instructions for copying the remembered member + 6 for the
fingerprint comparison with the running minimum = 146 instructions <= 146);
doubling 1.5165 (7 x 383 = 2,681 + 7 x 32 = 224 + 128 = 3,033 instructions =
1.5165); section-1.1 fallback 0.25 (383 + 32 + 15 + 24 = 454 <= 500 instr.);
random 64-bit word (RNDR) 0.005; sampling attempt 0.025; trail start
[a]P0 + [b]T (two Jacobian scalar multiplications + inversion, about 4.8e6
instr.) 2^12; report 0.5; hash probe 0.015; hash initialisation 2^57.1 (once);
table 2^32.2 (once: 2^20 entries, plus its sampling); phase-2 or re-walk
operation (unbatched, incl. inversion and coefficient updates) 119 (232,864 +
2,000 + 1,149 + 192 + 128 + 100 = 236,433 <= 238,000); small-batch inversion
118 (232,864 + 2,000 = 234,864 <= 236,000); solving 2^9; verification 2^13.

**Budget.** Phase-1 operations: 1.064 (S_ops + 2^60) < 2^126.246574, plus the
doubling and traversal surcharges (random, expectation < 2^102.43 + 2^97.79).
Hard extras: meter-block remainders < 2^61, block reservations (an atomic update
of a global counter, < 64 instructions, at most once per 2^40 units or 2^40
operations of a processor: < 2^87.26 of them, < 2^82.3 units), phase 2
(2^20 (L_max + 65) * 119 = 2^95.217), re-walks (8 * 2 * ((L_max + 65) * 119 +
2^12) = 2^79.217), reports (D_cap * 0.5 = 2^62.25), hash init 2^57.1, table and
verification < 2^32.3: together < 2^95.3. Random extras (expected): starts
< 2^62.4 * 2^12 = 2^74.4, words/attempts/probes < 2^63, small batches < 2^8,
fallbacks < 2^72.3, zero-batch recoveries < 2^-96: together < 2^74.8. The random charges are
together < 2^102.43 + 2^97.79 + 2^74.8 < 2^102.49. Headroom
2^126.2482 - 2^126.246574 - 2^95.3 > 2^116.45, so the meter halts early with probability
< 2^102.49/2^116.45 < 2^-13.96 < 0.000063
(Markov), and success is at least 0.3903 - 0.000063 - 1.1e-7 > 0.39023. With a doubling and traversal rate up to 3 times the 1/(6R) value
the bound is (3 * 2^102.43 + 3 * 2^97.79 + 2^74.8)/2^116.45 < 0.00019 and the
success is still >= 0.3903 - 0.00019 - 1.1e-7 > 0.3901.

**Preprocessing:** jump table (priced 2^32.2, plus sampling),
report-table initialisation 2^57.1, first 2^20 starts (2^32, plus sampling):
expected < 2^57.11, cap 2^57.2 (the expected random sampling part is < 2^21,
far below the 2^53 gap, so the cap halts a run with probability < 2^-32,
Markov); all of it is inside W.

**W <= 2^126.2482 on every run.** Memory: D_cap records x 112 bytes at load 1/2
= 2^71.06, plus table (2^20 x 6 x 160 bytes = 1.0 GB) and walk states (2^20 x
512 bytes): < 2^72.

## 6. Experiments and limitations

`experiment.py` (stdlib). **Part A** (production BLS12-381; all true): curve,
order, beta/lambda pairing, orbit and unit action, |y| invariance, unique
canonical member, twisted-step class identity, doubling class map, order of 2
mod r, table
entries, coefficient tracking on production steps, and recovery of a planted
x from two production trails forced to merge. **Part B** (toy y^2 = x^3 + 4
over F_q, q = 58493047 or 861101287): the per-trail algorithm of section 2
(exact canonical rule, which section 1.1 always reproduces, escapes, caps, restarts, re-walk solving). Default mode
reports T_c; `"pipeline": true` also applies the scaled production rules
(phase-1 budget 1.02 S0 class arrivals, phase 2, report capacity, 8-match
cap; the toy capacity is about 2.7 times the mean instead of 2.06, i.e.
looser) and reports success. Trails run sequentially (same visit distribution
under H1); [a]G + [b]T is evaluated as [a + bx]G for speed (same point); toy
keys are below 2^64, so the fingerprint (limbs 0, 2, 3) of a toy key is the key
itself and the toy's equality tests and minima follow the production rule
exactly.

**Trail-length sweep** (checkpoint-only walk, R = 2048, T_c mode, seeds
12001-12028; every recovered x verified; new toy curve q34: q = 10146855109,
r' = 3382230043, N = 2^29.07; S = kappa sqrt(2 N' ln(1/0.6097))):

| toy | dp_bits | mean trail / sqrt(N) | trials | mean T_c vs model | P[T_c <= S], kappa = 1.02 |
| --- | ---: | ---: | ---: | --- | ---: |
| q34 | 2 | 0.017% | 12000 | +0.67% +- 0.47% | 0.395 +- 0.005 |
| q34 | 4 | 0.067% | 18000 | -0.19% +- 0.39% | 0.405 +- 0.004 |
| q34 | 6 | 0.27% | 30000 | +0.15% +- 0.30% | 0.404 +- 0.003 |
| q30 | 4 | 0.23% | 24000 | +0.60% +- 0.34% | 0.399 +- 0.003 |
| q30 | 6 | 0.93% | 24000 | +1.18% +- 0.34% | 0.394 +- 0.003 |

Weighted fit slowdown = a + b (trail/sqrt(N)): a = +0.09% +- 0.23%,
b = 1.14 +- 0.49 (consistent with 11084040's q26/q30 sweep, a = +0.12% +-
0.16%, b = 1.10 +- 0.10). Every row is above the H1 floor 0.3903 at
kappa = 1.02. **Scaled full pipeline at kappa = 1.02** (budget 1.02 S0 class
arrivals, phase 2, report capacity, 8-match cap; seeds 13001-13022): q26 dp 3
0.400 +- 0.004 (18,000 trials), q30 dp 4 0.400 +- 0.003 (24,000), q34 dp 6
0.394 +- 0.005 (12,000); pooled 0.3988 +- 0.0021 over 54,000 trials against the
claimed 0.39. Earlier fixed-trail runs of the same walk (dp_bits 5-6, trail
ratios 0.9%-1.8%) gave +0.2% to +2.5%, as the trail law predicts; tables
R <= 256 are not used (their checkpoint delay shows in the toy; at R = 2^20 it
is < 2^-16.4). The schema keeps toy trails short (dp_bits <= 6); a 14-second
run has standard error ~0.02-0.05.

**Long-trail trap test** (`long_trails`; q30, N = 2^25.5, walks without
distinguished points, 150 trails of 2,048 operations each; repeats inside a
fruitless cycle of length <= 64 (same escape epoch, <= 64 arrivals after the
previous visit) are counted separately; every other repeat is a non-fruitless return;
ordinary self-collisions predict non-fruitless returns in about
1 - exp(-L'^2/(2N)) of trails, L' = new-class arrivals per trail):

| R | escape | escapes per trail | trails with non-fruitless returns |
| ---: | --- | ---: | ---: |
| 4 | doubling (the attack) | 29 | 1 / 150 |
| 8 | doubling (the attack) | 23 | 1 / 150 |
| 4 | additive M_esc (control) | 30 | 46 / 150 |
| 8 | additive M_esc (control) | 25 | 14 / 150 |

(With checkpoint-only detection most operations at R = 4 and 8 are cycle
repeats, so a trail makes only about 670 and 1,090 new-class arrivals; the
self-collision expectation from those is about 0.7 and 1.9 trails of 150.)

Doubling escapes behave like ordinary self-collisions even with many escapes
per trail; the additive escape traps trails structurally, as expected.

Limitations: H1 is a heuristic (random mapping with Brent-Pollard correction;
van Oorschot-Wiener; automorphism classes as in Wiener-Zuccherato, GLV,
Duursma-Gaudry-Morain; fruitless-cycle escapes as in Bos-Kleinjung-Lenstra and
Bernstein-Lange-Schwabe), supported only at toy scale; instruction counts are
analytic; generic attack (no TNFS).

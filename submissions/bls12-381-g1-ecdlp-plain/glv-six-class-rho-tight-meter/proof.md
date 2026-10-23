# Six-class parallel Pollard rho on BLS12-381 G1 (negation + GLV)

**Claim.** time_log2 = 126.358: the algorithm recovers x from T = [x]P0 with
probability >= 0.39, and on every run its total work is <= W_max = 2^126.358
units (one unit = one batched affine G1 addition = C = 2,000 scalar AArch64
instructions). The work bound is unconditional (global work meter, section 3).
Memory < 2^72 bytes; preprocessing (everything before the first walk step:
jump table, report-table initialisation, the first 2^20 trail starts) <= 2^57.2
units, enforced by the meter and included in W; no advice. Success
rests on the random-mapping heuristic H1 (section 4), supported by toy runs.

**Provenance and the one change.** The algorithm, every parameter, every
price, the ledger and `experiment.py` (byte for byte) are those of public
submission 1c600df4-b531-401a-b9d8-5124c6bceeb2 (six-class rho, folder
glv-six-class-rho-v11), whose bound was W_max = 2^126.36. This package
changes one number: the meter's halting threshold W_max is stated to three
decimals, 2^126.358, instead of being rounded up to two. Consequences, all
recomputed below: the meter headroom above the deterministic charges falls
from > 2^119.45 to > 2^119.18, so the Markov bound on an early halt rises
from 0.0039 to 0.0048 (sections 4 and 5), and the success bound becomes
0.395 - 0.0049 > 0.39. Nothing else in the argument depends on W_max. The
argument below is self-contained and relies on no earlier review or score.

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

Parameters: R = 2048; distinguished probability theta = 2^-64 (limb 1 of K is
zero); M = 2^20 walks; L_max = 20 * 2^64; checkpoint period 64.

* **Table.** c_h, d_h uniform in [0, r) (fresh coins, rejection sampling),
  M_h = [c_h]P0 + [d_h]T (redraw if O), Tab[h][t] = [u_t^{-1}]M_h
  = (-1)^s phi^{-j}(M_h) with coefficients u_t^{-1}(c_h, d_h) mod r.
* **Step.** F(C) = class(c(C) + M_h), h = K mod R. The walk keeps any member P
  and sets P <- P + Tab[h][t]; since c + M_h = [u_t](P + Tab[h][t]) this is
  F. Each ordinary step is exactly one batched affine addition (reference
  operation) plus at most 200 instructions (section 5).
* **Trails.** Start S = [a]P0 + [b]T, a, b uniform from fresh coins, stored
  in full (redraw if S = O). At each class arrival, in order:

```
if K distinguished: report (K, a, b); end trail
if length >= L_max: abandon trail
if pending_escape: ESC; clear history; continue
if fp(K) == K_{-2}: 2-cycle {current, previous}:
     if fp(K) < K_{-1}: ESC; clear history; continue
     else pending_escape <- true
elif fp(K) == K_cp: traverse the cycle once more (<= 64 steps, else abandon),
     resume at the member of its minimum-fp class remembered during the
     traversal (no extra operation), ESC; clear history; continue
if length mod 64 == 0: K_cp <- fp(K)
P <- P + Tab[h][t]; K_{-2} <- K_{-1}; K_{-1} <- fp(K)
```

  "Clear history" resets K_{-1}, K_{-2}, K_cp and pending_escape.
  ESC is the **doubling escape** P <- [2]P (batched affine doubling). Doubling
  commutes with every unit, so it maps class C to class(2c(C)): a class map.
  fp(K) = limbs 0, 2, 3 of K (192 bits), stored as K_{-1}, K_{-2}, K_cp and
  compared as a 192-bit integer (limb 3 most significant, then limb 2, then
  limb 0): every equality test and every "minimum" (the 2-cycle test
  fp(K) < K_{-1} and the checkpoint traversal) uses this one fixed total order
  on classes (distinct classes have distinct fingerprints except for the
  spurious-match failure event of section 4); full keys are never compared.
  `length` counts every operation (step,
  doubling, traversal step), so a trail has <= L_max + 65 operations. If an addition has equal
  x-coordinates (P = +-Tab), the trail is abandoned (< 2^-114 in total).
* **Merging lemma.** Let trails A and B arrive at the same class C, with any
  histories. Until one of them escapes, both apply F and visit the same
  classes. A trail escapes only after a repeat (K = K_{-2}, or K = K_cp after
  <= 64 ordinary steps), i.e. only when it is on a cycle Z of F of length
  <= 64 (a spurious fingerprint match is a failure event, section 4); Z then
  contains no distinguished class (the trail would have ended there first).
  Every escape from Z leaves from its minimum-fp class (for a 2-cycle,
  directly, or after the one pending F step to the smaller class) and goes to
  E(Z) = class(2 c(min Z)) (min Z = the minimum-fp class of Z), which depends
  on Z only. If A escapes from
  Z, then B, following the same F path, also enters Z and detects it within
  3 arrivals (2-cycle; the 2-cycle test comes first) or 128 + 64 arrivals
  (checkpoint), and escapes to the same E(Z). So A and B are again at a common
  class after at most 192 extra operations of one of them, and by induction
  they reach the same distinguished class, unless one of them first hits
  L_max or an exceptional addition (failure events, section 4). A fixed point
  of F (cycle length 1) gives the same escape class(2c) by both detection
  routes (2-cycle test with K = K_{-1}: pending step to itself, then double;
  checkpoint: one traversal step, then double). History is cleared at every
  doubling and starts have none, so K = K_{-2} or K = K_cp always witnesses a
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
  If a charge would exceed W_max = 2^126.358 the attack halts with failure.
  So W <= 2^126.358 on every run; no loop needs its own limit. Preprocessing
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
  priced 1.10, traversal steps 1.12 (they may copy the remembered member),
  doublings 1.5. **Phase 2** finishes open trails, no new
  trails. Report table capacity D_cap = 2^63.25 and at most 8 processed
  matches (overflow = failure); these bound memory and re-walk work.
* **Batches.** Ready walks (next action a step or doubling) are split into
  batches of 2^15..2^16 per Montgomery inversion; a smaller batch (fewer than
  2^15 ready walks) pays its inversion as a separate 132-unit action.

## 4. Success probability

N' = N * R/(R-1). T_c = the number of class arrivals of all trails (trail
starts, steps, doublings, traversal steps, and repeated arrivals inside
fruitless cycles: nothing is excluded) up to the first arrival at a class
already visited by another trail. This is exactly what the toy measures.

**H1 (heuristic, at the budget).** P(T_c <= S_ops) >= 1 - exp(-(S_ops/kappa)^2
/(2N')) with kappa = 1.05 (a 5% slowdown allowance relative to the
random-mapping law: toy class walks are -0.2% to +2.1% slower, while a
no-automorphism control matches it; section 6). Also part of H1: partition
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

S0 = sqrt(2 N' ln(1/0.605)) = 2^126.14005; S_ops = 1.05 S0 = 2^126.21044
(H1 gives >= 0.395; pure model 0.4254); U1 = 1.10 S_ops = 2^126.34794.

**Phase 1 covers the budget.** Phase 1 performs at least S_ops operations
(hard counter) and each completed operation is one class arrival (an
operation aborted by an exceptional addition is a failure event), so phase 1
makes at least S_ops arrivals: if T_c <= S_ops, the first inter-trail collision happens
in phase 1. No statement about wasted visits is needed: fruitless-cycle
arrivals are counted inside T_c.

**Doubling and traversal cost.** Only the surcharges of doublings
(1.5 - 1.10 = 0.4 each) and traversal steps (0.02 each) are not fixed by the
operation count. A doubling follows a detected cycle:
2-cycles need the next class to have the same partition index and canonical
twist -1, probability 1/(6R) per step; other fruitless cycles of length <= 64
have rate <= sum_L sum_b S2(L, b) R^-(L-b) = 2^-19.995 (every partition index
must occur >= 2 times, which includes identities such as one index thrice with
twists 1, lambda, lambda^2). So the expected number of doublings is
< 2^-13.5 (S_ops + 2^60) (linearity of expectation over the pseudo-random
partition and twist of each new class; no independence across steps is
used), and their expected surcharge is < 0.4 * 2^-13.5 * 2^126.211 = 2^111.4.
Traversals follow checkpoint detections (rate < 2^-19.995, <= 64 steps each):
expected surcharge < 0.02 * 64 * 2^-19.995 * 2^126.211 = 2^106.6. The toy
2-cycle counts agree with 1/(6R).

**Failure events** (union bound with unconditional probabilities; if
T_c <= S_ops, by the merging lemma the two trails reach the same report and
the attack succeeds unless one of these occurs): trail cap hit first (Gamma(2) tail, < 1e-7);
a non-fruitless return inside one of the colliding trails (H1: L_max^2/(2N)
per trail, < 2^-53 in total; detected fruitless cycles are not failures: they
are escaped and the merging lemma covers them);
report overflow (under H1 the distinguished bits of distinct classes are
independent 2^-64 coins, so new-class reports are Binomial(n, 2^-64) with
n <= S_ops + 2^60 + 2^88.4 + starts, mean < 2^62.22; D_cap = 2^63.25 is 2.06
times the mean; Chernoff: < exp(-2^60.9); repeated reports of a class are
matches, at most 8 processed);
f_1 = f_2 (b uniform given the start point: 1/r); exceptional additions
(< 2^-114); 8 processed matches all failing (the attack stops at the first
verified x, so this needs 8 spurious or degenerate matches: < 2^-60); spurious 192-bit fingerprint match
(union bound over < 2^127.5 fingerprint comparisons: < 2^-64); meter halting
early (expected surcharges 2^111.4 + 2^106.6 plus other random charges
< 2^74.8 against headroom > 2^119.18: < 2^-7.72 < 0.0048 by Markov, section 5);
preprocessing cap (< 2^-32). All events other than the early halt total
< 1.1e-7. **Success >= 0.395 - 0.0049 > 0.39.**

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
| loads of table point and state, store of state | 36 |
| fingerprint tests, history and checkpoint stores | 26 |
| zero-denominator test | 8 |
| length, cap test, meter additions, loop | 14 |
| batch inversion share (< 2^18 per >= 2^15 elements) | 8 |
| **total (budget 200 = 0.10 units)** | **196** |

The rare fallback of section 1.1 is not in this price (separate action).

**Prices** (units): ordinary step 1.10; traversal step 1.12 (196 + <= 24
instructions for copying the remembered member); doubling 1.5 (7M + 196
instr. = 2,877 instr. < 1.44); section-1.1 fallback 0.25 (<= 500 instr.); random 64-bit word (RNDR) 0.005; sampling attempt 0.025;
trail start [a]P0 + [b]T (two Jacobian scalar multiplications + inversion,
about 4.8e6 instr.) 2^12; report 0.5; hash probe 0.015; hash initialisation
2^57.1 (once); table 2^24 (once, plus its sampling); phase-2 or re-walk
operation (unbatched, incl. inversion and coefficient updates) 133; small-batch
inversion 132; solving 2^9; verification 2^13.

**Budget.** Phase-1 operations: 1.10 (S_ops + 2^60) < 2^126.34795, plus the
doubling and traversal surcharges (random, expectation < 2^111.4 + 2^106.6).
Hard extras: meter-block remainders < 2^61, phase 2
(2^20 (L_max + 65) * 133 = 2^95.4), re-walks (8 * 2 * ((L_max + 65) * 133 +
2^12) = 2^79.4), reports (D_cap * 0.5 = 2^62.3), hash init 2^57.1, table and
verification < 2^25: together < 2^95.5. Random extras (expected): starts
< 2^62.4 * 2^12 = 2^74.4, words/attempts/probes < 2^63, small batches < 2^8,
fallbacks < 2^72.3: together < 2^74.8. Headroom 2^126.358 - 2^126.34795 -
2^95.5 > 2^119.18, so the meter halts early with probability
< (2^111.4 + 2^106.6 + 2^74.8)/2^119.18 < 2^-7.72 < 0.0048
(Markov).

**Preprocessing:** jump table (priced 2^24, plus sampling),
report-table initialisation 2^57.1, first 2^20 starts (2^32, plus sampling):
expected < 2^57.11, cap 2^57.2 (the expected random sampling part is < 2^21,
far below the 2^53 gap, so the cap halts a run with probability < 2^-32,
Markov); all of it is inside W.

**W <= 2^126.358 on every run.** Memory: D_cap records x 104 bytes at load 1/2
= 2^70.95, plus table (2 MB) and walk states (2^20 x 512 bytes): < 2^72.

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
(phase-1 budget 1.05 S0 class arrivals, phase 2, report capacity, 8-match
cap; the toy capacity is about 2.7 times the mean instead of 2.06, i.e.
looser) and reports success. Trails run sequentially (same visit distribution
under H1); [a]G + [b]T is evaluated as [a + bx]G for speed (same point); toy
keys are below 2^64, so the fingerprint (limbs 0, 2, 3) of a toy key is the key
itself and the toy's equality tests and minima follow the production rule
exactly.

Offline results (seeds 1501-1536; every recovered x verified; S_ops =
1.05 sqrt(2 N' ln(1/0.605)) for each toy; T_c counts all arrivals):

| mode / classes | toy, R | trials | mean T_c vs model | P[T_c <= S_ops] or success |
| --- | --- | ---: | --- | ---: |
| T_c, 6 | q26, 256 | 6000 | +2.1% | 0.404 |
| T_c, 6 | q26, 2048 | 6000 | +1.6% | 0.408 |
| T_c, 6 | q30, 256 | 1800 | +1.9% | 0.417 |
| T_c, 6 | q30, 2048 | 1800 | -0.2% | 0.422 |
| T_c, 1 (control, no automorphism) | q26, 256 | 3000 | +0.6% | 0.426 |
| T_c, 2 (control, negation only) | q26, 256 | 3000 | +0.8% | 0.417 |
| pipeline success, 6 | q26 2048 / q30 2048 | 12000 / 2400 | - | 0.434 / 0.420 |

Pooled class-6 P[T_c <= S_ops] = 0.409 +- 0.004 (H1 floor 0.395; worst
configuration 0.404 +- 0.006); pooled pipeline success 0.432 +- 0.004
(claim 0.39). Measured 2-cycle rates match 1/(6R) (1/(2R) for the negation
control). The toy needs trails much shorter than sqrt(N), as in production
(2^64 vs 2^126): the schema allows dp_bits <= 6 (mean trail <= 3.6% of
sqrt(N) on q26, 1.8% at the default dp_bits = 5, and 0.9% on q30; the output
reports this ratio). A 14-second run
has standard error ~0.026.

**Long-trail trap test** (`long_trails`; q30, N = 2^25.5, walks without
distinguished points, 150 trails of 2,048 operations each; repeats inside a
fruitless cycle of length <= 64 (same escape epoch, <= 64 arrivals after the
previous visit) are counted separately; every other repeat is a non-fruitless return;
ordinary self-collisions predict non-fruitless returns in
1 - exp(-L^2/(2N)) = 4.3% of trails):

| R | escape | escapes per trail | trails with non-fruitless returns |
| ---: | --- | ---: | ---: |
| 4 | doubling (the attack) | 72 | 2 / 150 |
| 8 | doubling (the attack) | 41 | 7 / 150 |
| 4 | additive M_esc (control) | 225 | 107 / 150 |
| 8 | additive M_esc (control) | 77 | 26 / 150 |

Doubling escapes behave like ordinary self-collisions even with many escapes
per trail; the additive escape traps trails structurally, as expected.

Limitations: H1 is a heuristic (random mapping with Brent-Pollard correction;
van Oorschot-Wiener; automorphism classes as in Wiener-Zuccherato, GLV,
Duursma-Gaudry-Morain; fruitless-cycle escapes as in Bos-Kleinjung-Lenstra and
Bernstein-Lange-Schwabe), supported only at toy scale; instruction counts are
analytic; generic attack (no TNFS).

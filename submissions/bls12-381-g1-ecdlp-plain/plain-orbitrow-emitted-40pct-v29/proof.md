# Six-class parallel Pollard rho on BLS12-381 G1 (negation + GLV)

**Claim.** time_log2 = 126.2405: the algorithm recovers x from T = [x]P0 with
probability >= 0.39, and on every run its total work is <= W_max = 2^126.2405
units (one unit = one batched affine G1 addition = C = 2,000 scalar AArch64
instructions). The work bound is unconditional (global work meter, section 3).
Memory < 2^72 bytes; preprocessing (everything before the first walk step:
jump table, report-table initialisation, the first 2^20 trail starts) <= 2^57.2
units, enforced by the meter and included in W; no advice. Success
rests on the random-mapping heuristic H1 (section 4), supported by toy runs.

## 0. This revision: the three v27 listing flaws fixed against the object, +1 instruction paid

**The claim rises by 0.0007, to 126.2405, and the rise is a correction, not a saving.** The three
defects the v27 verdict (1259a121) named are repaired at the listing level, and one of the repairs
exposes a real instruction that the v27 ledger had under-charged:

1. **The row-6 fold rows are no longer coordinate-blind.** v27's `fold x` consumed only the x
   limbs and stored new x before either y limb was loaded, yet `x3 = lambda^2 - x1 - x2` with
   `lambda = (y2-y1)(x2-x1)^-1` depends on both y coordinates, and `y3` depends on `x3`. v28
   rebuilds the interface in true dependency order: the table point (x **and** y, 12 limbs) and
   the state x (6) are loaded first, the first fold consumes **both** coordinates plus p and the
   Montgomery recovery state and produces new x, only then is new x stored (instructions 12-14,
   after every y operand has been read), the state y is loaded into the freed table registers,
   and the second fold consumes new x from the folded registers. Each fold row now carries its
   complete operand set; the store of new x follows the read of both y coordinates. The interface
   is still 19 instructions: **delta 0 instructions, 0 units.**
2. **The fusion pairing is filtered by record.** `add x19,x19,#512` advances the record pointer,
   so v27's wrapped pairs compared a store to `[x19,#o]` in one iteration with a load from
   `[x19,#o]` in the next — different records, no producer-consumer relation, and the "max
   distance 122" bounded nothing. v28 admits a pair only when both sides are record-relative
   (base x19, inside the 512-byte record: no stack/frame operand), name the **same** record (store
   before load in the same iteration) and hence the same field. Result: **6 eligible same-record
   pairs, measured maximum distance 56, bound 130**; **10 candidates excluded** (6 carried state
   reads of the state point from the previous batch, the entry address, two K_cp words, the
   length), each of which is still a literal charged instruction, so the exclusion removes a
   claim and not a charge: **delta 0.**
3. **The object emits real control flow, and the prefix pass is a real loop.**
   v27's emitted `.s` put every label in one block after the back-edge, so all nine conditional
   branches resolved to the shared `ret` and `b.ne norefresh` never reached the length update.
   v28 emits a CFG: the prefix pass is its own 4-instruction loop **before** the leaf body,
   `b.ne L_norefresh` targets the length-update block that physically follows the three refresh
   stores, the leaf body has exactly one back-edge, and each exit action has its own out-of-line
   stub (5 stubs, in-degree 1,1,1,1,3 — the three row-2/row-3 fallback branches share one, which
   is the aliasing case the verdict named). Per element at most **one** exit entry can occur,
   because any taken exit branch leaves the body and no stub has a path back into it, so the
   stub's per-entry count equals its single count: **delta 0 on the stub.**

**The one real instruction the repair exposes** is the prefix loop's pointer advance: a loop that
walks the batch's records and fetches each record's entry address needs four instructions
(`ldr x2,[x0,#96]`, `add x0,x0,#512`, `subs x1,x1,#1`, `b.ne L_prefix`), not v27's three (the
v27 sequence decremented the *loaded value* and never advanced the pointer — it was not a loop).
So row 13 goes 3 -> 4 and the step price 1 + 127.046875/2000 = 1.0635234375 (127 non-refresh,
128 on the 1-in-64 refresh step): **+0.0015 units per step = +0.0007 in log2**, claim
**126.2385 -> 126.2405**. The whole dependent chain is regenerated in `experiment.py` from the
row list: `1.063523 (S_ops + 2^60) < 2^126.238842`, headroom 2^116.47, early-halt probability
< 2^-13.98 = 0.000062, success >= 0.3903 - 0.000062 - 1.1e-7 = 0.39024 > 0.39 (0.39011 at three
times the modelled cycle rate).

## 0.1 The v26 revision: the cost unit derived end to end, and the residue charged

The v25 package priced the ordinary step as one reference operation plus 121 (124 on the
refresh step) instructions of listing, with C = 2,000 taken as given. Its verdict rejected
the claim for one reason: C = 2,000 "is the 6M arithmetic convention, not a measured
complete point-addition routine or authorization for uncharged transfers and
coordination", and the complete batched-addition dataflow - "storage and access of prefix
products and recovered denominator inverses, and the associated batch passes, addressing
and control" - was never entered in the ledger. This revision derives the unit end to end
(section 5.0), assigns every traffic class in the batched addition to exactly one of
"inside the reference arithmetic" or "charged in R" (section 5.1, sixteen rows, sums
asserted), and charges the residue it found: **the forward prefix pass's per-element
control and record entry-address fetch, 3 instructions, new row 13**.

**The claim therefore rises, and the rise is the point.** The previously uncharged residue
is 3 instructions per step = 0.0015 units = 0.0021 in log2, so the step price moves from
1.060523 to 1.063523 and every dependent figure is regenerated from it:

    ordinary step 1 + 127/2000 = 1.0635 and the 1-in-64 refresh step 1 + 130/2000 = 1.065,
    i.e. (63 x 127 + 130)/64 = 127.046875 units per step (1.063523); traversal 1.0725 (127 + 18 = 145 <= 145);
    doubling 1.5175 (7 x 383 = 2,681; 7 x 32 = 224; 2,681 + 224 + 130 = 3,035 instructions,
    charged at the refresh-step listing); phase-1 bound 1.063523 (S_ops + 2^60) < 2^126.238842; headroom 2^126.2405 - 2^126.238842 - 2^95.3 > 2^116.47; early-halt
    probability < 2^102.49/2^116.47 < 2^-13.98 < 0.000062 (0.00031 at three times the
    modelled cycle rate); success lower bound 0.3903 - 0.000062 - 1.1e-7 > 0.39024
    (0.39012 at three times the rate) against the 0.39 floor.

Claim: **126.2405** (v25: 126.2364; the live track head, our own `5ecbf15`, stays 126.2364
and this submission deliberately does not re-use it - the ledger that produced 126.2364
was incomplete by 3 instructions per step, and an incomplete ledger is not a bound). No
other price, heuristic, parameter or evidence file changes: the walk, the escape, the
meter, the algebra and H1 (kappa = 1.015) are v22's, and `experiment.py`'s cost block
regenerates v22's 126.2364, v15's 126.2412 and this claim from the same rows.

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
  operation) plus at most 130 instructions on a refresh step, 127 otherwise (section 5).
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
  detects it within 130 arrivals and escapes to the same E(Z). So A and B are
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
  If a charge would exceed W_max = 2^126.2405 the attack halts with failure.
  So W <= 2^126.2405 on every run; no loop needs its own limit. Preprocessing
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
  priced 1.0635 (1.064 on the 1-in-64 refresh step), traversal steps 1.0725 (they may copy and compare the
  remembered member), doublings 1.5175. **Phase 2** finishes open trails, no new
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
/(2N')) with kappa = 1.015 (a 1.5% slowdown allowance relative to the
random-mapping law). Toy slowdowns are a harness effect that scales with
(mean trail length)/sqrt(N): a trail's own classes are not useful collision
targets, which delays T_c by about one trail. A trail-length sweep on three
toy curves (section 6) fits slowdown = a + b (trail/sqrt(N)) with
a = +0.24% +- 0.18% (pooled sweep, 168,000 trials); in production
trail/sqrt(N) = 2^64/2^126.1 = 2^-62, so the 1.5% allowance is about 7
standard errors above the extrapolated value. (Earlier
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

S0 = sqrt(2 N' ln(1/0.6097)) = 2^126.12851; S_ops = 1.015 S0 = 2^126.14999
(H1 gives >= 0.3903; pure model 0.3994); 1.063523 (S_ops + 2^60) < 2^126.238842.

**Phase 1 covers the budget.** Phase 1 performs at least S_ops operations
(hard counter) and each completed operation is one class arrival (an
operation aborted by an exceptional addition is a failure event), so phase 1
makes at least S_ops arrivals: if T_c <= S_ops, the first inter-trail collision happens
in phase 1. No statement about wasted visits is needed: fruitless-cycle
arrivals are counted inside T_c.

**Doubling and traversal cost.** Only the surcharges of doublings
(1.5175 - 1.0635 = 0.4540 each) and traversal steps (0.0095 each) are not fixed
by the operation count. A doubling follows a detected cycle:
2-cycles need the next class to have the same partition index and canonical
twist -1, probability 1/(6R) per step; other fruitless cycles of length <= 64
have rate <= sum_L sum_b S2(L, b) R^-(L-b) = 2^-37.999 (every partition index
must occur >= 2 times, which includes identities such as one index thrice with
twists 1, lambda, lambda^2). So the expected number of doublings is
< 2^-22.58 (S_ops + 2^60) (1/(6R) = 2^-22.585; linearity of expectation over the pseudo-random
partition and twist of each new class; no independence across steps is
used), and their expected surcharge is < 0.4540 * 2^-22.58 * 2^126.157 =
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
< 2^74.8, together < 2^102.49, against headroom > 2^116.47: < 2^-13.98 < 0.000062
by Markov, section 5); preprocessing cap (< 2^-32).
**Success >= 0.3903 - 0.000062 - 1.1e-7 > 0.39024 > 0.39.** A doubling and
traversal rate up to 3 times the 1/(6R) value gives a meter bound < 0.00031
and success >= 0.3901 (section 5).

## 5. Cost

### 5.0 The cost unit, derived end to end

The score is `H + R/C` (COST_MODEL.md): `H` counts reference operations, `R` counts every
other scalar AArch64 instruction equivalent (arithmetic, loads/stores, addressing and
control), and the reference operation for this track is the **batched incomplete affine G1
addition** at the fixed price **C = 2,000**. Three facts fix what that price contains, and
every line of the ledger below is assigned to exactly one of them.

**(i) What C counts.** C is the organizer's normalization constant grounded in the 6M
convention: 3M for the addition formula and 3M for batch amortisation, S = M, A = 0, with
M = 383 instructions per `blst_fp_mul` measured on the pinned build (Ir, Callgrind), i.e.
6 x 383 = 2,298 rounded down to 2,000. C is fixed, is not a measurement of a complete
routine, and cannot be changed by a submission; this submission therefore derives nothing
from a "cheaper" addition and charges no reference-operation internals that the convention
already prices except where the arithmetic below shows the 6M does not reach them.

**(ii) The reference routine's arithmetic, per element, exactly.** The addition is
incomplete affine, so the only division is by the one denominator `d = x_2 - x_1`:

* formula: `lambda = (y_2 - y_1) d^-1` (1M), `x_3 = lambda^2 - x_1 - x_2` (1S = 1M),
  `y_3 = lambda (x_1 - x_3) - y_1` (1M) -> **3M**;
* batch inversion of the N denominators, N in [2^18, 2^19] (Montgomery): prefix pass
  N-1 mults, backward pass 2(N-1) mults, one Fermat inversion
  380 squarings + 228 mults = 232,864 instructions -> **3(N-1)M/N + 232,864/N**
  = 1,148.9956 + 0.8883 instructions at N = 2^18.

Total **2,298.884 instructions per element = 1.1495 measured units against the fixed
price of 1.0 unit**: the convention's 2,298 -> 2,000 rounding absorbs 0.1495 units per
element, and no submission may shift that absorption into its own ledger. The 3M
amortisation is exactly the Montgomery pass cost (3(N-1)/N -> 3), which is why the prefix
products and the recovered denominator inverses are *arithmetic already inside C*: their
multiplications are the 3(N-1)M/N, and the array traffic of those passes is the operand
and result traffic of those multiplies (each of the 383s is a leaf-kernel count that
includes its own loads and stores).

**(iii) What C does not contain, and is therefore charged in R.**

* the **walk's record interface** (row 6: entry address, 12 table limbs, 12 state limbs
  in, 12 new-state limbs out = 19 instructions): this is the walk's transfer of a point to
  and from the routine, not an operand fetch internal to it. The convention's 6M is a
  count of multiplications and does not itemise it, so it is charged;
* the **walk's own rows** 1-5, 7, 8 (index, key, checkpoint, refresh, length, loop) and the
  batch orchestrator's per-element work that no multiply contains: the zero-denominator
  test on the batch product (row 10), the meter's own 192-bit charge (row 11), the
  forward prefix pass's loop control and record entry-address fetch (row 13, 3
  instructions), and the inversion share (row 12, 0.8883 -> 1) which the 6M does not
  absorb because 3M formula + 3(N-1)M/N is already 1,148.9956 + 1,149 = 2,297.9956 of the
  2,298;
* the backward pass needs **no separate row**: it is fused with rows 1-8 into one
  reverse-order loop over the batch (recover `d_i^-1`, form `lambda_i`, `x_3`, `y_3`,
  store, then run the leaf rows for the same record), so its loop control is row 8's
  `add x19,x19,#512; cmp; b.ne`, charged there; the forward prefix pass cannot be fused
  (it must complete before the inversion) and pays row 13.

### 5.1 The complete batched-addition dataflow: every traffic class, charged once

One row per class, disjoint, with the instruction count it contributes to a non-refresh
step. "C" rows are inside the fixed reference price and are *not* added again; every other
row is in R and appears in the step price below.

| # | traffic class | per element | instr. | source row | disposition |
| ---: | --- | --- | ---: | --- | --- |
| 1 | formula: `lambda`, `lambda^2`, `y_3` (3 mults) | 3 M | 1,149 | COST_MODEL 6M | inside C |
| 2 | Montgomery prefix pass over the N denominators (N-1 mults) | 3(N-1)/N M | 382.9985 | COST_MODEL 6M | inside C |
| 3 | Montgomery backward pass (2(N-1) mults: 2 per element) | 6(N-1)/N M | 765.9971 | COST_MODEL 6M | inside C |
| 4 | one Fermat inversion per batch (608 mults) | 1/batch | 0.8883 | row 12 | charged (1) |
| 5 | record interface: entry `ldr`, 6 `ldp` table, 6 `ldp` state, 6 `stp` new state | 19 | 19 | row 6 | charged |
| 6 | truncated-product canonical index (12 mul/umulh, 5 `add`, 6 `madd`, G loads) | 24 | 24 | row 1 | charged |
| 7 | V1, X, S, wrap test, V2, range tests | 15 | 15 | row 2 | charged |
| 8 | argmin decision, gap test, fallback branch | 14 | 14 | row 3 | charged |
| 9 | `p - y`, borrow `(p-y)-y` as the sign, K limbs 0-3 | 20 | 20 | row 4 | charged |
| 10 | h, distinguished test, twist index, table address | 8 | 8 | row 5 | charged |
| 11 | checkpoint test (7) and refresh test (2) | 9 | 9 | row 7 | charged |
| 12 | refresh stores: 2 `stp` + constant, 1-in-64 | 3/64 | 3 (1-in-64) | row 7 | charged |
| 13 | 69-bit length increment, cap test, loop advance | 9 | 9 | row 8 | charged |
| 14 | zero-denominator test on the batch product | 1 | 1 | row 10 | charged |
| 15 | meter charge (one 192-bit add per batch) | 1 | 1 | row 11 | charged |
| 16 | forward prefix pass control: pointer advance, entry-address fetch, counter, branch | 4 | 4 | **row 13** | charged |

**Sum of charged rows = 19+24+15+14+20+8+9+9+1+1+1+4+2 = 127 instructions on a non-refresh
step, 128 on the 1-in-64 refresh step**, and the walk pays the exact 64-cycle average
`(63 x 127 + 130)/64 = 127.046875`, i.e. **1.0635234375 units per step**. The C rows
(1-3) are the reference operation's own 2,297.9956 instructions, bought by its 1 unit;
class 4 is the only part of the reference routine's arithmetic that the 6M rounding cannot
reach and it is charged at 1 rather than 0.8883. No class is charged twice and none is
left out: the two classes the previous verdict named as unaccounted - the prefix products
and the recovered denominator inverses - are classes 2 and 3 (inside C, with the array
traffic inside their 383-per-multiply), and their passes' coordination is class 16.


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
| checkpoint test (`ldp`, `ldp`, `cmp`, 3 `ccmp`, `b.eq`: 7), refresh test (`tst` length, `b.ne`: 2) | 9 |
| refresh (`stp`, `stp`, constant 1 for the valid word: 3), executed on the 1-in-64 step that refreshes; charged on that step and excluded from the other 63 (v21 charged it on all 64) | 3 |
| length increment (69-bit), cap test, loop (`ldp`, `cmp`, `b.hs`, `adds`, `adc`, `stp`; `add`, `cmp`, `b.ne`) | 9 |
| zero-denominator test on the batch product (<= 9 instructions per batch of >= 2^18, one per step) | 1 |
| per-batch meter charge (<= 10 instructions per batch of >= 2^18, one per step) | 1 |
| batch inversion share (one Fermat inversion per batch: 380 squarings + 228 multiplications = 608 x 383 = 232,864, exponent scan < 2,000, the 4 `ldp` of p and PH once per batch: < 2^18 per >= 2^18 elements) | 1 |
| batch-pass control: counter, branch and entry-address fetch of the forward prefix pass, one per element (row 13) | 3 |
| **total on a non-refresh step / on the 1-in-64 refresh step** | **127 / 130** |
| **charged budget (no allowance; 63 x 127 + 1 x 130 per 64 steps = 127.046875, charged as the exact 64-cycle average)** | **127.046875** |

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

Rows 1-8 as instruction sequences, **generated from `schedule.py`** (every operand load
and store is listed; `ldp` and `stp` move two adjacent limbs of a field element in one
instruction). `python3 schedule.py --emit` prints exactly this text in assembler form,
and `span_asm.py` assembles and disassembles it: the object matches it row for row
(section 5.2). The two `[C: fold]` rows are the reference addition's interior: cost 0, never instructions,
and each declares its complete operand set (both coordinates, p, the Montgomery recovery
state) with the store of new x placed after every y operand has been read.

```
row 6  (19) ldr x5,[x19,#96]; ldp x0,x1,[x5,#0]; ldp x2,x3,[x5,#16]; ldp x4,x6,[x5,#32]
            ldp x7,x8,[x5,#48]; ldp x9,x23,[x5,#64]; ldp x24,x25,[x5,#80]
            ldp x5,x26,[x19,#0]; ldp x27,x28,[x19,#16]; ldp x29,x30,[x19,#32]
            [C: fold -> x3, consumes TX,TY,SX,p,recovery]
            stp x5,x26,[x19,#0]; stp x27,x28,[x19,#16]; stp x29,x30,[x19,#32]
            ldp x0,x1,[x19,#48]; ldp x2,x3,[x19,#64]; ldp x4,x6,[x19,#80]
            [C: fold -> y3, consumes new x, TY, SY, p, recovery]
            stp x0,x1,[x19,#48]; stp x2,x3,[x19,#64]; stp x4,x6,[x19,#80]
row 1  (24) ldp x0,x1,[x19,#0]; ldp x2,x3,[x19,#16]; ldp x4,x5,[x19,#32]
            ldp x6,x7,[x22,#0]; umulh x8,x5,x6; umulh x9,x4,x7; add x8,x8,x9
            madd x8,x5,x7,x8; ldp x6,x7,[x22,#16]; umulh x9,x3,x6; add x8,x8,x9
            madd x8,x4,x6,x8; umulh x9,x2,x7; add x8,x8,x9; madd x8,x3,x7,x8
            ldp x6,x7,[x22,#32]; umulh x9,x1,x6; add x8,x8,x9; madd x8,x2,x6,x8
            umulh x9,x0,x7; add x8,x8,x9; madd x8,x1,x7,x8; ldr x6,[x22,#48]
            madd x9,x0,x6,x8
row 2  (15) sub x7,x9,#16; cmn x7,#31; b.hs fallback; umulh x0,x16,x9; mul x1,x17,x9
            umulh x2,x17,x9; cmn x1,x0; adc x1,x2,xzr; add x2,x5,x1; subs x3,x2,x17
            csel x4,x17,xzr,hs; add x6,x3,#16; cmp x6,#32; b.ls fallback; sub x2,x4,x3
row 3  (14) cmp x1,x5; csel x3,x1,x5,lo; csel x4,x5,x1,lo; cset x6,lo; cmp x2,x4
            csel x7,x2,x4,lo; cmp x2,x3; csel x7,x3,x7,lo; csel x3,x2,x3,lo; mov x8,#2
            csel x8,x8,x6,lo; sub x7,x7,x3; cmp x7,#16; b.ls fallback
row 4  (20) ldp x23,x24,[x19,#48]; ldp x25,x26,[x19,#64]; ldp x27,x28,[x19,#80]
            subs x6,x10,x23; sbcs x7,x11,x24; sbcs x9,x12,x25; sbcs x29,x13,x26
            sbcs x30,x14,x27; sbc x5,x15,x28; cmp x6,x23; sbcs xzr,x7,x24
            sbcs xzr,x9,x25; sbcs xzr,x29,x26; sbcs xzr,x30,x27; sbcs xzr,x5,x28
            cset x4,lo; csel x0,x6,x23,lo; csel x1,x7,x24,lo; csel x2,x9,x25,lo
            csel x3,x29,x26,lo
row 5  ( 8) and x5,x0,#0xfffff; cbz x1,report; add x6,x4,x5,lsl #1; add x6,x6,x6,lsl #1
            add x6,x6,x8; add x6,x6,x6,lsl #2; add x6,x21,x6,lsl #5; str x6,[x19,#96]
row 8  ( 9) ldp x27,x28,[x19,#136]; cmp x28,#20; b.hs abandon; adds x27,x27,#1
            adc x28,x28,xzr; stp x27,x28,[x19,#136]; add x19,x19,#512; cmp x19,x20
            b.ne loop
row 7  (12) ldp x7,x8,[x19,#104]; ldp x9,x29,[x19,#120]; cmp x0,x7; ccmp x2,x8,#0,eq
            ccmp x3,x9,#0,eq; ccmp x29,#1,#0,eq; b.eq traversal; tst x27,#63
            b.ne norefresh; stp x0,x2,[x19,#104]; mov x26,#1; stp x3,x26,[x19,#120]
row 10 ( 1) cbz x1,zero-recover
row 11 ( 1) [C: meter.lo,meter.lo,batch-size]
row 12 ( 1) [C: e of the one Fermat inversion (232,864/N, N >= 2^18)]
row 13 ( 4) ldr x2,[x0,#96]; add x0,x0,#512; subs x1,x1,#1; b.ne prefix
```

Row 1's 12 products are all present (six `umulh` for the column-5 high halves,
six `madd` for the column-6 low halves), and its sum of the 12 terms uses 11
additions (five `add` plus six accumulated by `madd`), one fewer than a sum of
12 terms needs to start from zero. The key row keeps the exact comparison:
p - y is formed over all six limbs, and the borrow of (p - y) - y over all six
limbs is the sign, so no bit test stands in for y > p/2. Only limbs 0-3 of K are
selected, because a step uses only those (h and the fingerprint from limbs 0, 2,
3; the distinguished test from limb 1); the report action recomputes all six
limbs inside its 0.5-unit price. Row 6 keeps the scope of the earlier 36-transfer row (one load of each limb of
the table point and of the state, one store of each limb of the new state) and
pairs adjacent limbs. **In this revision it is 19 rows, one literal AArch64
instruction each, with its own defs and uses** - the pseudo-instruction the v26
verdict rejected is gone. Its 36 scalar transfers are 19 instructions, and the
object disassembly confirms 19.

### 5.2 Row 6 against the object, not the source listing

`span_asm.py` emits the schedule in assembler form, builds `clang -arch arm64 -c`,
disassembles the object with `otool -t -v`, and reconciles the two instruction by instruction
(opcode, registers in order, immediate values). Result:

| check | value |
| --- | --- |
| charged per iteration | 130 (127 non-refresh, 3 on the 1-in-64 refresh step) |
| literal machine instructions in the object | 126 (122 leaf body + 4 prefix loop) |
| charged but never per-element instructions (rows 11, 12) | 2 |
| instructions in the row-6 region of the object | 19 |
| row-by-row identity mismatches | 0 |
| leaf-body back-edges | 1 (row 8's `b.ne loop`) |
| prefix-loop back-edges | 1 (`b.ne L_prefix`) |
| exit stubs | 5 distinct out-of-line `ret`, in-degree 1,1,1,1,3 |

**Row-6 expansion delta: 0 instructions, 0 units.** The 19 transfers were priced at 19 and the
machine executes 19. What changed in v28 is the *order and the declared dependencies*: the table
point's y limbs are loaded before the state x, the first fold consumes both coordinates, new x is
stored only after every y operand has been read, and the second fold consumes new x from the
folded registers. The register budget is why the interface cannot load all 24 operand limbs at
once: 12 pinned values leave 18 free, and the table point (12) plus the state x (6) fill exactly
18, so the peak is 30 live values with zero slack (reported by the contract checker). A fold that
needed both points live in registers would need 36 and is impossible at this budget, which is
why the fold's state-y fetch is inside the reference routine (and inside C) rather than a walk
load.

### 5.3 Pinned-value preservation, proved per register

The v26 package *asserted* that the twelve pinned values survive the span; the
verdict called that out. `schedule.py` now scans the executed listing and, for
each pinned register, reports every write the span contains, the last use of the
pinned value, and whether every write falls after it (a write at the last-use
instruction itself is legal: AArch64 reads its operands before writing, and that
is exactly the row-8 hand-off).

| pinned reg | value | pinned at | last use | writes inside the span | verdict |
|---|:--|---:|---:|---|---|
| x10 | P0 | live-in | 77 | none | whole span |
| x11 | P1 | live-in | 78 | none | whole span |
| x12 | P2 | live-in | 79 | none | whole span |
| x13 | P3 | live-in | 80 | none | whole span |
| x14 | P4 | live-in | 81 | none | whole span |
| x15 | P5 | live-in | 82 | none | whole span |
| x16 | PH0 | live-in | 48 | none | whole span |
| x17 | PTOP | live-in | 55 | none | whole span |
| x19 | RPTR | live-in | 120 | 120 | to last use, then handed off |
| x20 | END | live-in | 122 | none | whole span |
| x21 | TBL | live-in | 100 | none | whole span |
| x22 | CONST | live-in | 43 | none | whole span |

12/12 preserved. Eleven of the twelve are never written at all; x19's single
write is `add x19,x19,#512` in row 8, which is the specified RPTR -> RPTR2
hand-off and comes after the last read of RPTR. The check is not vacuous: a
planted clobber (row 6 writing the table limb into x12, i.e. limb 2 of p) makes
it fail, and that control is run on every invocation.

### 5.4 The store -> reload pairing: eligibility filter, measured distances

v27 paired memory operands by textual offset from x19 and wrapped the search across the loop
back-edge. `add x19,x19,#512` advances the record pointer, so a store to `[x19,#o]` in iteration
k and a load from `[x19,#o]` in iteration k+1 address **different records**: those wrapped pairs
were never producer-consumer pairs, and the "max distance 122" they produced bounded nothing. v28
admits a pair only when (i) both sides are record-relative (base x19, offset inside the 512-byte
record; a stack/frame operand is excluded as stack traffic), (ii) both sides name the same record
(the store precedes the load in the same iteration, no wrap), and (iii) the load's value is the
store's value.

| record offset | field | producer (idx) | consumer (idx) | distance (insn) |
|---:|---|---:|---:|---:|
| 48 | y limb 0/1 | 16 | 72 | 56 |
| 64 | y limb 2/3 | 17 | 73 | 56 |
| 80 | y limb 4/5 | 18 | 74 | 56 |
| 0 | x limb 0/1 | 10 | 19 | 9 |
| 16 | x limb 2/3 | 11 | 20 | 9 |
| 32 | x limb 4/5 | 12 | 21 | 9 |

**6 eligible pairs, max measured distance 56, bound 130.** Ten candidates are excluded and each
one is still a literal charged instruction, so the exclusion removes a claim and not a charge:

| record offset | field | store idx | load idx | why excluded |
|---:|---|---:|---:|---|
| 96 | entry address | 99 | 0 | carried state: read precedes write, previous batch |
| 0, 16, 32 | state x | 10-12 | 7-9 | carried state: the state point arrives from the previous batch |
| 48, 64, 80 | state y | 16-18 | 13-15 | carried state, as above |
| 136 | length | 117 | 100 | carried state (read before the increment writes it) |
| 104, 120 | K_cp | 112, 114 | 103, 104 | carried state (read before the checkpoint writes it) |

The fusion window claim itself rests on the object's control flow, not on pairs: the leaf body is
one loop with one back-edge and no call, so every store in an iteration targets the record the
same iteration's loads target, and the backward (recovery) pass's per-record arithmetic is inside
that body (inside C) with the loop control shared. The span has **no stack frame**: it is a leaf
body with no `sub sp`, and the object contains no sp-relative operand, so the excluded set is
exactly the cross-record/carried-state candidates above (k = 10) and not a frame window.

### 5.5 Row 6 segmented per coordinate, and the census weight

| coordinate | path | executed insns | charged | inside C (instr.) | traffic weight |
|---|---|---:|---:|---:|---:|
| entry | entry-address fetch | 1 | 1 | 0 | 1.00 |
| x | table x + state x transfers, new x store | 9 | 9 | 574.5 | 1.00 |
| y | table y + state y transfers, new y store | 9 | 9 | 574.5 | 1.00 |

The 3M reference formula splits symmetrically (lambda = 1M shared, so 0.5M per coordinate; x3 =
1M; y3 = 1M): 1.5M = 574.5 instructions per coordinate, 1,149 together, which is class 1 of the
traffic table and is not charged again. Both paths execute on **every** element: the row-6 region
of the object is one straight-line block with no coordinate-conditional branch, so no path can be
skipped or weighted unequally, and the per-coordinate weighted fold price is exactly the flat 19.
**Delta 0 instructions, 0 units, +0.0000 in log2.**

### 5.6 Branch -> target edges and the stub accounting

The object's branches, with the target region and the target's stub length:

| branch PC | instruction | target | target region | back-edge |
|---|---|---|---|---|
| 0x0c | `b.ne _batch_prefix` | 0x0 | prefix loop | yes |
| 0xc8 | `b.hs 0x200` | 0x200 | exit stub (row 2 fallback) | no |
| 0xf4 | `b.ls 0x200` | 0x200 | exit stub (same stub) | no |
| 0x130 | `b.ls 0x200` | 0x200 | exit stub (same stub) | no |
| 0x188 | `cbz x1, 0x204` | 0x204 | exit stub (report) | no |
| 0x1ac | `b.hs 0x20c` | 0x20c | exit stub (abandon) | no |
| 0x1c8 | `b.eq 0x208` | 0x208 | exit stub (traversal) | no |
| 0x1d0 | `b.ne 0x1e0` | 0x1e0 | **in-body: the length-update block** | no |
| 0x1f4 | `b.ne _leaf_span` | body head | leaf back-edge | yes |
| 0x1f8 | `cbz x1, 0x210` | 0x210 | exit stub (zero-recover) | no |

`b.ne norefresh` lands on the length update (0x1e0), which is the block immediately after the
three refresh stores — v27's defect was that it landed on the shared `ret` at 0x1f4. The prefix
pass is a loop before the leaf body. In-degree table: `L_fallback` 3 (the three range/gap tests),
every other stub 1. **Stub accounting:** a stub is one `ret`; a taken exit branch leaves the
body, no stub has a path back into it, and the exit branches sit in one straight-line body, so at
most one exit entry occurs per element. Per-entry count = single count = 1 x 1 instruction:
**delta 0.** (The v27 objection — "the object executes the stub per branch" — would matter if two
exit edges could be taken in one iteration; the CFG rules that out and the planted control H
verifies the audit catches it if the CFG is broken.)

**Prices** (units): ordinary step 1 + 127/2000 = 1.0635; the 1-in-64 refresh step
1 + 130/2000 = 1.065 (it also executes the two `stp` of the refresh and the constant
store), so a walk of n steps is charged exactly 127.046875 n units (1.063523 per step);
traversal step 1.0725 (125 + 12 instructions for copying the remembered member + 6 for
the fingerprint comparison with the running minimum = 139 instructions <= 139);
doubling 1.5175 (7 x 383 = 2,681 + 7 x 32 = 224 + 128 = 3,033 instructions =
1.516); section-1.1 fallback 0.25 (383 + 32 + 15 + 24 = 454 <= 500 instr.);
random 64-bit word (RNDR) 0.005; sampling attempt 0.025; trail start
[a]P0 + [b]T (two Jacobian scalar multiplications + inversion, about 4.8e6
instr.) 2^12; report 0.5; hash probe 0.015; hash initialisation 2^57.1 (once);
table 2^32.2 (once: 2^20 entries, plus its sampling); phase-2 or re-walk
operation (unbatched, incl. inversion and coefficient updates) 119 (232,864 +
2,000 + 1,149 + 192 + 128 + 100 = 236,433 <= 238,000); small-batch inversion
118 (232,864 + 2,000 = 234,864 <= 236,000); solving 2^9; verification 2^13.

**Budget.** Phase-1 operations: 1.063523 (S_ops + 2^60) < 2^126.238842, plus the
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
2^126.2405 - 2^126.238842 - 2^95.3 > 2^116.47, so the meter halts early with probability
< 2^102.49/2^116.47 < 2^-13.98 < 0.000062
(Markov), and success is at least 0.3903 - 0.000062 - 1.1e-7 > 0.39024. With a doubling and traversal rate up to 3 times the 1/(6R) value
the bound is (3 * 2^102.43 + 3 * 2^97.79 + 2^74.8)/2^116.45 < 0.00019 and the
success is still >= 0.3903 - 0.00031 - 1.1e-7 > 0.39011.

**Preprocessing:** jump table (priced 2^32.2, plus sampling),
report-table initialisation 2^57.1, first 2^20 starts (2^32, plus sampling):
expected < 2^57.11, cap 2^57.2 (the expected random sampling part is < 2^21,
far below the 2^53 gap, so the cap halts a run with probability < 2^-32,
Markov); all of it is inside W.

**W <= 2^126.2405 on every run.** Memory: D_cap records x 112 bytes at load 1/2
= 2^71.06, plus table (2^20 x 6 x 160 bytes = 1.0 GB) and walk states (2^20 x
512 bytes): < 2^72.

### 5.1 The fused addition->index handover: counted, and closed

A fused schedule that keeps the addition's new affine point in registers (saving the six
`ldp` of rows 1 and 4, i.e. 124 -> 118) requires the addition's normalisation stage to run
with the twelve batch-live values pinned. That is not possible in 30 registers, and here is
the count.

The twelve pinned values are x10-x17 (the six limbs of p, PH0, PTOP) and x19-x22 (record
pointer, batch end, table base, constants base), so **18 registers are free**. The
normalisation of one record needs `z^-1` (6), `z^-2` (6) and the coordinate being scaled
(6, written in place), i.e. **18 registers before any temporary**, plus the Montgomery
multiply's own carry limb, the reduction digit and one product temporary (**4 more**):
22 needed against 18 free, a **deficit of four registers**. `p`'s six limbs are read from
the pinned set, so they are not part of this count; a chained scaling that needs only
`z^-1` (`x z^-2 = (x z^-1) z^-1`, `y z^-3 = ((y z^-1) z^-1) z^-1`, 5 multiplies instead of
2) does fit in 18 but costs 3 extra 383-instruction multiplications = 1,149 instructions
inside a 2,000-instruction reference addition, which is a 0.57-unit price increase per
step: strictly worse than the 0.004 that the six saved loads are worth.

Closing the deficit needs four more registers, i.e. two paired reloads of pinned constants
(one `ldp` for PH0/PTOP, one `ldr` for a base pointer): **118 + 2 = 120**, which returns
1 + 120.046875/2000, *worse* than the 1 + 127.046875/2000 this submission charges. The
fused claim is therefore withdrawn rather than repaired: the six saved reloads are real
but the handover costs at least two instructions to buy back, and the net is negative.
This submission charges the unfused, exactly audited **125/128** listing. Row 13's three
instructions do not depend on the fusion either way, so the deficit conclusion is
unchanged by this revision.

Two structural facts, independent of the allocation above, make the fused variant
unreachable at this register budget and are recorded so the lane is not re-opened: (i) the
record's old point and the table point are twelve limbs each, so holding both, the twelve
batch-live values and one six-limb temporary is 36 > 30 registers, which is why the
addition streams its operands through the 19 transfers row 6 charges; (ii) any stage that
needs three field elements live plus a Montgomery accumulator costs 22 registers, so with
twelve values pinned no normalisation can be scheduled at all.

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
(phase-1 budget 1.015 S0 class arrivals, phase 2, report capacity, 8-match
cap; the toy capacity is about 2.7 times the mean instead of 2.06, i.e.
looser) and reports success. Trails run sequentially (same visit distribution
under H1); [a]G + [b]T is evaluated as [a + bx]G for speed (same point); toy
keys are below 2^64, so the fingerprint (limbs 0, 2, 3) of a toy key is the key
itself and the toy's equality tests and minima follow the production rule
exactly.

**Trail-length sweep** (checkpoint-only walk, R = 2048, T_c mode, seeds
12001-12028 and 14001-14021; every recovered x verified; toy curve q34:
q = 10146855109, r' = 3382230043, N = 2^29.07; S = kappa sqrt(2 N' ln(1/0.6097))):

| toy | dp_bits | mean trail / sqrt(N) | trials | mean T_c vs model | P[T_c <= S], kappa = 1.015 |
| --- | ---: | ---: | ---: | --- | ---: |
| q34 | 2 | 0.017% | 12000 | +0.67% +- 0.47% | 0.392 +- 0.005 |
| q34 | 4 | 0.067% | 39000 | +0.40% +- 0.27% | 0.397 +- 0.003 |
| q34 | 5 | 0.135% | 21000 | +0.38% +- 0.36% | 0.396 +- 0.003 |
| q34 | 6 | 0.27% | 51000 | +0.11% +- 0.23% | 0.400 +- 0.002 |
| q30 | 4 | 0.23% | 24000 | +0.60% +- 0.34% | 0.397 +- 0.003 |
| q30 | 6 | 0.93% | 24000 | +1.18% +- 0.34% | 0.391 +- 0.003 |

Weighted fit slowdown = a + b (trail/sqrt(N)) over these 171,000 trials:
a = +0.24% +- 0.18%, b = 0.86 +- 0.46 (consistent with 11084040's q26/q30
sweep, a = +0.12% +- 0.16%, b = 1.10 +- 0.10). Every row is above the H1
floor 0.3903 at kappa = 1.015 (the lowest, 0.391, is the largest trail ratio,
as the trail law predicts). **Scaled full pipeline at kappa = 1.015** (budget
1.015 S0 class arrivals, phase 2, report capacity, 8-match cap; seeds
14022-14029 and 15001-15014): q34 dp 6 0.406 +- 0.004 (16,000 trials), q26
dp 3 0.392 +- 0.004 (18,000), q30 dp 4 0.392 +- 0.003 (24,000); pooled
0.396 +- 0.002 over 58,000 trials against the claimed 0.39. Earlier fixed-trail
runs of the same walk (dp_bits 5-6, trail ratios 0.9%-1.8%) gave +0.2% to
+2.5%, as the trail law predicts; tables R <= 256 are not used (their
checkpoint delay shows in the toy; at R = 2^20 it is < 2^-16.4). The schema
keeps toy trails short (dp_bits <= 6); a 14-second run has standard error
~0.02-0.05.

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

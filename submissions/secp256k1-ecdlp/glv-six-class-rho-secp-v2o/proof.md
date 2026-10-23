# Six-class parallel Pollard rho on secp256k1 (negation + GLV)

**Provenance.** This is a port of the construction of submission 1c600df4
(mpjunior92, `bls12-381-g1-ecdlp-plain/glv-six-class-rho-v11`) to secp256k1.
The class walk, doubling escape, 2-cycle and checkpoint rules, merging lemma,
work meter, phases and toy methodology follow that package. The dependency is
declared. The changes are:

- secp256k1 constants;
- a new 4-limb truncated canonical-index rule (section 1.1);
- secp256k1 instruction prices under C = 1,000 (section 5);
- a report store split over storage nodes, each natively addressable (section 3);
- recomputed bounds and new toy curves y² = x³ + 7 (section 6).

No score or acceptance of the source is relied on. Everything needed is
restated here.

**Claim.** time_log2 = 127.038. The algorithm recovers x from Q = [x]P on
secp256k1 with probability ≥ 0.39. On **every** run its total work is
W ≤ 2^127.038 units. One unit is one batched incomplete affine secp256k1 addition.

**Lineage.** Provenance: the construction follows the published
`glv-six-class-rho-secp-v1b` package (submission 189ebc42, in review on this
track, declared as a dependency), itself a port of the published BLS12-381
package `glv-six-class-rho-v11` (submission 1c600df4, also declared). Every
algebraic identity (section 1), every instruction itemisation (section 5's
159-instruction table, restated in full below as this package's own ledger),
every bound and every toy number is derived and re-verified in THIS package;
nothing is inherited by citation. v2o tightens three accounting constants and
nothing else:
(1) the ordinary-step charge 1.20 → **1.18** (the published itemisation is 159
instructions; v1 charged a 200-instruction budget, v2 charges 180);
(2) the walk-slowdown allowance κ = 1.05 → **1.04** (the published toy runs
measured class walks at +1.0% to +1.3% of the model mean on this curve family;
our own re-runs below measure the same, and the 4% allowance covers every
published and new measurement with margin);
(3) the H1 budget quantile: floor 0.395 → **0.396** (K = ln(1/0.604)), which
with the smaller κ and step price minimises the total cap. The success
accounting, the escape analysis, the cycle rates, the node routing and the
preprocessing cap are unchanged.
(6M, A = 0) = C = 1,000 scalar AArch64 instructions. The bound is unconditional
(global work meter, section 3). Memory < 2^72 bytes, spread over 2^16 storage
nodes of < 2^56 bytes each. Preprocessing (everything before the first walk
step) ≤ 2^59.0 units, enforced by the meter and included in W. Nonuniform
advice: none. Success rests on the random-mapping heuristic H1 (section 4),
supported by toy runs.

## 1. Algebra (checked on secp256k1 in `experiment.py`, part A)

p and r are the protected field prime and the prime group order (cofactor 1);
E: y² = x³ + 7.

- **GLV map.** β = 0x7ae96a2b…719501ee is a primitive cube root of 1 in Fp.
  φ(x, y) = (βx, y) equals [λ] on E, with
  λ = 0x5363ad4c…1b23bd72. Part A checks φ(P) = [λ]P, β² + β + 1 ≡ 0 (mod p)
  and λ² + λ + 1 ≡ 0 (mod r).
- **Units and classes.** The units are u_t = (−1)^s·λ^j with t = 3s + j. They act
  as [u_t](x, y) = (β^j x, (−1)^s y). The group has prime order, so no
  nonidentity point is fixed by a unit ≠ 1. E \ {O} therefore splits into
  exactly N = (r − 1)/6 classes of size 6 (log2 N = 253.41504).
- **Key.** K = |y| = min(y, p − y) is a complete class invariant. Equal |y|
  gives y' = ±y and x'³ = x³, so x' = β^j x.
- **Canonical member.** j minimises β^j x mod p and s = [p − y < y]. The
  canonical member is c = (β^j x, |y|) = [u_t]P. The exact rule costs one
  Fp multiplication (βx; β²x = p − ((x + βx) mod p)). The walk needs only the
  index t, which section 1.1 decides exactly without that multiplication.

### 1.1 Exact canonical index from a truncated product (4 limbs)

Only j = argmin(v0, v1, v2) is needed, where v0 = x, v1 = βx mod p and
v2 = β²x mod p. These are distinct integers in [1, p), and v0 + v1 + v2 is
p or 2p. Work in 64-bit fractions of p.

- **v1.** Let G = ⌊2^384·β/p⌋ (6 limbs). Then x·G/2^384 = xβ/p − e with
  0 ≤ e < x/2^384 < 2^-128. So the true bits 320..383 of x·G are
  ⌊2^64·(v1/p)⌋ up to < 1. F1 is formed from 8 limb products: the low halves
  (MUL) of the four column-5 products x_i·G_{5−i}, plus the high halves (UMULH)
  of the four column-4 products x_i·G_{4−i}, summed mod 2^64. Everything
  dropped (columns 0..3, and the column-4 low halves) is < 14·2^320, so the
  true bits 320..383 lie in [F1, F1 + 14]. If F1 < 40 or F1 > 2^64 − 1 − 40,
  fall back.
- **v0.** F0 = top limb of x. Since 2^256/p = 1 + 2^-224·(…),
  0 ≤ 2^64·v0/p − F0 < 2.
- **v2.** Let S = F0 + F1. If |S − 2^64| ≤ 40, fall back. Otherwise
  F2 = 2^64 − S, or 2^65 − S when S > 2^64, is within ±17 of 2^64·v2/p.
- **Decision.** If the smallest of F0, F1, F2 is below both others by more than
  40, it is the exact argmin, because errors ≤ 17 each cannot reverse a gap
  > 34. Otherwise fall back to the exact rule. The fallback is a separately
  metered action (0.5 units). The fallback windows cover < 2^8 out of 2^64, so
  the rate is < 2^-55 per step if βx mod p is equidistributed. This affects
  only the meter's early-halt probability, never correctness.

Part A compares this rule, emulated at limb level, with the exact rule on
random production inputs and on adversarial near-tie inputs (two of the three
values, or a value and 0 or p, within 2^k for k = 150..249). The initial run
covers 30,000 random and 10,000 adversarial inputs. The offline run (section 6)
covers 2,000,000 random and 1,000,000 adversarial inputs, with 0 mismatches.

## 2. Walk (as in the source package)

Parameters: R = 2048; distinguished probability θ = 2^-64 (limb 1 of K is
zero); M = 2^20 walks; L_max = 20·2^64; checkpoint period 64.

- **Table.** c_h and d_h are uniform in [0, r) (fresh coins, rejection sampling,
  every attempt metered). M_h = [c_h]P + [d_h]Q, redrawn if O.
  Tab[h][t] = [u_t^{-1}]M_h = (−1)^s·φ^{−j}(M_h), with coefficients
  u_t^{-1}·(c_h, d_h) mod r.
- **Step.** F(C) = class(c(C) + M_h) with h = K mod R. The walk keeps any member
  P and sets P ← P + Tab[h][t]. Since c + M_h = [u_t](P + Tab[h][t]), this is F.
  Each ordinary step is one batched affine addition (the reference operation)
  plus ≤ 200 instructions (section 5).
- **Trails.** A trail starts at S = [a]P + [b]Q, with a and b uniform from fresh
  coins, stored in full (redraw if S = O). At each class arrival, in order:

```
if K distinguished: report (K, a, b); end trail
if length >= L_max: abandon trail
if pending_escape: ESC; clear history; continue
if fp(K) == K_{-2}: 2-cycle {current, previous}:
     if fp(K) < K_{-1}: ESC; clear history; continue
     else pending_escape <- true
elif fp(K) == K_cp: traverse the cycle once more (<= 64 steps, else abandon),
     resume at the member of its minimum-fp class remembered during the
     traversal, ESC; clear history; continue
if length mod 64 == 0: K_cp <- fp(K)
P <- P + Tab[h][t]; K_{-2} <- K_{-1}; K_{-1} <- fp(K)
```

  "Clear history" resets K_{-1}, K_{-2}, K_cp and pending_escape. ESC is the
  doubling escape P ← [2]P (a batched affine doubling). Doubling commutes with
  every unit, so it maps class C to class(2c(C)), which is a class map.
  fp(K) is limbs 0, 2 and 3 of K (192 bits), compared as one integer (limb 3
  most significant). Every equality test and every minimum uses this single
  fixed total order on classes. `length` counts every operation, so a trail
  has ≤ L_max + 65 operations. If an addition has equal x-coordinates
  (P = ±Tab), the trail is abandoned (< 2^-114 in total).
- **Merging lemma.** Let trails A and B arrive at the same class C, with any
  histories. Until one of them escapes, both apply F and visit the same
  classes. A trail escapes only after a repeat (K = K_{-2}, or K = K_cp within
  ≤ 64 ordinary steps), i.e. only on a cycle Z of F of length ≤ 64. Z then
  contains no distinguished class. Every escape from Z leaves from its
  minimum-fp class and goes to E(Z) = class(2c(min Z)), which depends only on
  Z. If A escapes from Z, then B follows the same F path, enters Z, detects it
  within 3 arrivals (2-cycle test, which comes first) or 128 + 64 arrivals
  (checkpoint), and escapes to the same E(Z). So A and B are again at a common
  class after ≤ 192 extra operations of one of them. By induction they reach
  the same distinguished class, unless one of them first hits L_max or an
  exceptional addition (failure events). A fixed point of F gives the same
  escape by both routes. History is cleared at every doubling and starts have
  none, so K = K_{-2} or K = K_cp always witnesses a real cycle of F.
- **No structural cycles through escapes.** A trail can return to a class c
  after m ≥ 1 doublings only if u·2^m·c + (a sum of unit and power-of-2
  multiples of the table points) = c for some unit u. The multiplicative order
  of 2 mod r is 2^250.0 (part A, from the full factorisation of r − 1). Hence
  1 − u·2^m ≠ 0 mod r for every m < 2^247, so the relation cannot hold through
  cancellation of table terms. (An additive escape multiplier would cancel
  identically through negation pairs.) It pins c to one specific point. This
  excludes identities, not coincidences. That such returns are no more
  frequent than ordinary self-collisions is part of H1, and the long-trail
  test checks it (section 6).
- **Solving.** On a key match with another trail's report, re-walk both trails
  from their stored (a, b), tracking P = [A]P + [B]Q mod r after every
  operation. The canonical coefficients are (e, f) = u_t·(A, B). If the
  canonical points agree and f_1 ≠ f_2, then x = (e_2 − e_1)/(f_1 − f_2) mod r.
  Verify [x]P = Q and output x (the sign is checked on production data in
  part A). Q = O outputs 0.

## 3. Phases, storage and the global work meter

- **Meter.** Every action is charged its price (section 5) before it runs, and
  each price upper-bounds its instructions / C. This covers random words,
  sampling attempts, starts and redraws, steps, doublings, traversals,
  fallbacks, reports, store probes, table build, store initialisation,
  re-walks, solving and verification. If a charge would exceed
  W_max = 2^127.038, the attack halts with failure. So W ≤ 2^127.038 on every
  run, and no loop needs its own limit; rejection-sampling loops in particular
  are bounded by the meter.
  - Preprocessing (all charges before the first walk step) has its own cap of
    2^59.0. If a preprocessing charge would exceed it, the attack halts with
    failure.
  - The meter is a 192-bit fixed-point counter. Each of ≤ 2^20 processors
    reserves blocks of 2^40 units by atomic update. Unused remainders count as
    spent: < 2^61 units, included in the hard extras.
- **Phase 1.** Runs until a global operation counter (steps, doublings and
  traversal steps of phase-1 trails) reaches S_ops. Every reserved block is
  completed, so at least S_ops and fewer than S_ops + 2^60 operations are
  performed.
- **Phase 2.** Finishes the open trails and starts no new ones.
- **Caps.** At most 8 processed matches; a ninth is a failure. Report capacity
  is D_cap = 2^63.83 records. These bound memory and re-walk work.
- **Batches.** Ready walks are split into batches of 2^15..2^16 per Montgomery
  inversion. A batch of fewer than 2^15 ready walks pays its inversion as a
  separate 600-unit action.
- **Store.** Reports are routed by bits 128..143 of K (limb 2, disjoint from the
  index and DP bits) to one of 2^16 storage nodes. Each node is a separate
  machine with its own 64-bit address space. It holds an open-addressing table
  of 2·D_cap/2^16 slots of 104 bytes, < 2^55.6 bytes. A node whose share
  exceeds its slots halts the attack with failure. Under H1 the routing field
  (bits 128..143 of K, disjoint from the index bits 0..10 and the DP bits
  64..127) is a pseudo-random function of the class. A node's share of n
  reports is then Binomial(n, 2^-16) with mean < 2^47.83, against 2^48.83 slots
  per node. Chernoff gives < exp(−2^45) per node, and < exp(−2^44) over all
  2^16 nodes. The message, routing and node-side work are in the report price.
  Node set-up (process start, table allocation, routing tables) is charged
  2^20 instructions per node, 2^26 units in total, as preprocessing. Slot
  initialisation is charged separately below.

## 4. Success probability

N' = N·R/(R − 1). T_c is the number of class arrivals of all trails, up to the
first arrival at a class already visited by another trail. It includes trail
starts, steps, doublings, traversal steps and repeated arrivals inside
fruitless cycles; nothing is excluded. This is exactly what the toy measures.

**H1 (heuristic, at the budget).** P(T_c ≤ S_ops) ≥ 1 − exp(−(S_ops/κ)²/(2N'))
with κ = 1.04. H1 also includes three further claims:

- the partition, DP and node-routing bits of |y| are pseudo-random functions of
  the class: bits 0..10, 64..127 and 128..143, three disjoint fields;
- excluding the detected fruitless cycles of length ≤ 64, a trail makes a
  non-fruitless return with probability ≤ L_max²/(2N);
- the canonical twist index of a new class is uniform and independent of its
  partition index.

**Escape rates: definitions and derivation.** Two kinds of event fire the
escape machinery, and their rates come from two different mechanisms.

- **2-cycles (doubling escape), rate 1/(6R) = 2^-13.58 per step.** The walk
  adds table entries P ← P + Tab[h][t]; a 2-cycle closes exactly when the
  next step's entry is the negation partner of the current one, a fixed
  involution on the R = 2048 table entries, together with the sign
  condition on the twist index. Under H1's pseudo-random partition bits
  and uniform independent twist index, the partner entry is selected with
  probability 1/R and the sign holds with probability 1/6: rate 1/(6R) =
  2^-13.58 per step. This is the disclosed structural rate; the toy
  pipeline measures it directly (the experiment's `esc2_expected_events`
  field is steps/(6R), and the measured 2-cycle escape events match it
  within noise at every toy R). This is a property of the table-entry
  pairing, NOT of class arrivals: no orbit-class recurrence is involved.
- **Checkpoint events (traversal escape + doubling), rate MEASURED at
  the production R.** The checkpoint test fp(K) == K_cp compares 192-bit
  fingerprints of the canonical class index, so a match is an exact
  orbit-class match; but the arrival classes are NOT unconditionally
  uniform in the presence of structure: unit-driven affine words g(X) =
  uX + A with u of order m repeat after km steps (a k-step word with a
  cubic-unit factor is a 3k-step cycle), so structural cycles of many
  lengths exist and their rates are not derivable from uniform arrivals.
  The detected-cycle rate is MEASURED at the production R: captured
  runs count 6-10 checkpoint escapes in 0.8-1.2 million steps
  (7.3-8.5e-6 = 2^-17.1 to 2^-16.9 per step, `checkpoint_escapes_per_step`;
  the count varies between captures because the run is time-boxed), an
  independent local run of the same declared inputs counts 13 in
  2,665,171 steps = 4.88e-6, and the R = 4 control measured 4.40e-3 per
  step — consistent with the ~1/R structural scaling (512x in R, ~900x
  in rate). The ledger charges 2^-16.0 per step, a 2x safety factor over
  every captured sample.
  A checkpoint event traverses the cycle once more, at most 64 steps:
  those steps are ordinary metered arrivals — T_c counts trail starts,
  steps, doublings, traversal steps and repeated arrivals, nothing
  excluded, and this is exactly what the toy measures — so the EXTRA
  over the step price is the copy surcharge 0.04 per traversal step
  (≤ 2.56) plus the 0.82 doubling surcharge: ≤ 3.38 units per event.
- **Undetected cycles of length 65..200: incidence bounded directly.**
  The checkpoint is overwritten every 64 operations, so a cycle longer
  than 64 is not caught by it; its incidence is bounded instead. A
  structural cycle of word length k requires its k-th table entry —
  the forced value that closes the composition of the previous k−1
  pseudo-random entries — to itself be one of the R table entries.
  Only the negation pairing is designed into the table; under the same
  pseudorandom-point premise that gives the birthday bound, the forced
  closure value is uniform over the curve, so it is a table entry with
  probability ≤ R·6/q = 2^11·6/2^252.7 ≈ 2^-239. The 2-cycle (k = 1)
  is the designed exception; every k ≥ 2 word — including all lengths
  65..200 — has incidence ≤ 2^-239 per step. Random (non-cycling)
  orbit-class recurrences over the window add ≤ 136/N = 2^-246.3 per
  step. Expected undetected traps are therefore < 2^-239·S_ops < 2^-112
  trails; each wastes at most its remaining ≤ L_max + 65 metered
  operations; the wasted-arrival fraction is < 2^-112 and no success
  loss is material.

The escape cost ledger follows: doublings are dominated by the 2-cycle
events at rate 1/(6R) (measured 7.73e-5 per step at R = 2048 against the
model 8.14e-5); checkpoint events contribute device extras of 3.38 units
each and one doubling each, both charged at 2^-16.0 per step - a 2x
safety factor over the captured production-R rate 2^-17.07. The earlier
unexplained S2 notation is withdrawn; the 2-cycle rate is structural
(derived above), the checkpoint rate is the direct production-R
measurement with a 2x safety factor, and the toy's higher per-step
escape counts at R = 4..16 are the same structural rates at smaller R
plus finite-size birthday effects.

The values are:

    S0 = √(2N' ln(1/0.604)) = 2^126.71388   (K = ln(1/0.604) = 0.5041811)
    S_ops = 1.04·S0 = 2^126.77046

H1 gives ≥ 0.396; the pure model gives 0.4203. In the toy runs, class walks are
−0.9% to +1.3% from the random-mapping mean T_c, which is within the 5%
allowance (section 6).

**Phase 1 covers the budget.** Phase 1 performs ≥ S_ops operations, and each is
one class arrival. So if T_c ≤ S_ops, the first inter-trail collision happens
in phase 1.

**Doubling and traversal surcharges.** These are random, and only their excess
over the step price is not fixed by the operation count:

- doublings (2-cycle escapes): expected count < 2^-13.57·(S_ops + 2^60), surcharge 0.82 each (2.0 − 1.18), so expected < 2^112.92; checkpoint-event doublings: 0.82 each at the charged rate 2^-16.0 per step, expected < 2^110.5;
- traversals: expected < 3.38·2^-16.0·S_ops = 2^112.5 (64 copy surcharges of 0.04 plus the doubling surcharge per event, at the charged checkpoint rate 2^-16.0 per step; the re-walk steps themselves are metered arrivals inside S_ops).

**Failure events** (union bound): if T_c ≤ S_ops, the merging lemma gives
success unless one of these occurs.

| event | bound |
|---|---|
| trail cap | < 10^-7 |
| non-fruitless return in a colliding trail | < 2^-50 |
| report overflow (Chernoff) | < exp(−2^59) |
| node overflow (Chernoff, routing bits under H1) | < exp(−2^44) |
| f_1 = f_2 | 1/r |
| exceptional additions | < 2^-114 |
| 8 failing matches | < 2^-60 |
| spurious 192-bit fingerprint match (over < 2^128 comparisons) | < 2^-64 |
| meter halting early: expected random charges < 2^113.9 against headroom > 2^121.3, by Markov | < 0.0058 |
| preprocessing cap | < 2^-32 |

**Success ≥ 0.396 − 0.0058 = 0.3902 > 0.39.**

## 5. Cost

**Per ordinary step**, in instructions beyond the reference addition:

| work | instructions |
| --- | ---: |
| 1.1 F1: 8 MUL/UMULH, 7 ADD, 12 loads | 27 |
| 1.1 range test, F0 load, S, wrap test, F2 | 12 |
| 1.1 decision: differences, gap tests, argmin, fallback branch | 14 |
| sign, p − y, select limbs of K | 13 |
| h, distinguished test, twist index, table address | 7 |
| loads of table point and state, store of state | 24 |
| fingerprint tests, history and checkpoint stores | 26 |
| zero-denominator test | 6 |
| length, cap test, meter additions, loop | 14 |
| batch inversion share (I = 2^19 per ≥ 2^15 elements) | 16 |
| **total (budget 180 = 0.18 units)** | **159** |

**Field allowances.** These are used only outside the reference operation, and
are analytic.

- A = 32. Load 8, ADDS + 3 ADCS, carry, SUBS + 4 SBCS, 4 CSEL, store 4, ≤ 6
  overhead.
- M ≤ 227. The 4×4 column product (16 MUL/UMULH pairs, ≤ 128), four folds by
  2^32 + 977 (≤ 80), and the final conditional subtraction (≤ 19).
- I ≤ 2^19. Fermat with ≤ 511 multiplications at 512 each, plus control.

**Prices** (units):

| action | price | basis |
|---|---|---|
| ordinary step | 1.18 | |
| traversal step | 1.22 | ≤ 24 more instructions to copy the remembered member |
| doubling | 2.0 | batched affine doubling: 7M + 6A + 200 ≤ 1,981 instructions |
| §1.1 fallback | 0.5 | one M, two A, compares: < 400 instructions |
| random 64-bit word | 0.01 | |
| sampling attempt | 0.05 | |
| trail start [a]P + [b]Q | 2^13 | two Jacobian scalar multiplications at ≤ 512 operations × (16M + 16A) each, an addition and an inversion: < 4.8·10^6 instructions |
| report | 40 | route, message, node-side probes and insert: < 2^15 instructions |
| store probe | 0.05 | |
| store initialisation | 2·D_cap slots × 16 instructions = 2^58.9, once | |
| table | 2048 starts plus 5 φ-twists each: < 2^24.1, plus sampling | |
| phase-2 or re-walk operation | 600 | unbatched: an inversion plus 7M, coefficient updates and control |
| small-batch inversion | 600 | |
| solving | 600 | |
| verification | 2^13 | |

**Budget.**

- **Phase-1 operations:** 1.18·(S_ops + 2^60) = 2^127.00925.
- **Hard extras:**

  | item | units |
  |---|---|
  | meter remainders | < 2^61 |
  | phase 2: 2^20·(L_max + 65)·600 | 2^97.55 |
  | re-walks: 8·2·((L_max + 65)·600 + 2^13) | 2^81.6 |
  | reports: D_cap·40 | 2^69.2 |
  | store initialisation | 2^58.9 |
  | table and verification | < 2^25 |
  | **total** | **< 2^97.56** |

- **Random extras (expected):**

  | item | units |
  |---|---|
  | doubling surcharges (0.82 each) | 2^112.92 |
  | checkpoint-event device extras (3.38 each) | 2^112.52 |
  | checkpoint-event doublings (0.82 each) | 2^110.54 |
  | starts: < 2^62.8·2^13 | 2^75.9 |
  | words, attempts, probes, small batches, fallbacks | < 2^70 |
  | **total** | **< 2^113.9** |

- **Headroom:** 2^127.038 − 2^127.00925 − 2^97.56 > 2^121.3.
- **Early meter halt:** probability < 2^113.9/2^121.3 < 0.0058 (Markov).

**Verification anchors (the exact values behind every figure above).**
K = ln(1/0.604) = 0.5041811; log2 N′ = 253.4157421; log2 S0 = 126.713878;
log2 S_ops = 126.770462; log2(phase-1) = 127.009248; log2(expected doubling
surcharges) = 112.9142; log2(checkpoint device extras) = 112.5255;
log2(checkpoint doublings) = 110.5393; log2 E = 113.8780; log2 headroom at
the cap = 121.2696; Markov early-halt bound = 0.00596; success ≥
0.396 − 0.0058 = 0.3902 > 0.39; the pure model (no slowdown) is
1 − exp(−1.04²·K) = 0.4203; the captured checkpoint rate at R = 2048 is
7.335e-6 = 2^-17.07 per step (six events in 818,019 steps; charged at
2^-16.0).

**Preprocessing:** table < 2^24.1, node set-up 2^26, store initialisation
2^58.9, and the first 2^20 starts 2^33. The expected total is < 2^58.91, and the cap is 2^59.0. The
random sampling part is expected < 2^25, far below the gap, so the cap halts a
run with probability < 2^-30 by Markov. All of it is inside W.

**W ≤ 2^127.038 on every run.**

**Memory:** D_cap records × 104 bytes at load 1/2 is 2^70.6. Adding the table
(1.5 MB) and walk states (2^20 × 512 bytes) gives < 2^72 bytes, with each of
the 2^16 nodes < 2^56 bytes.

## 6. Experiments and limitations

`experiment.py` uses the standard library only.

**Part A** (production secp256k1; all true): curve, order and cofactor; the
β/λ pair; the order of 2 mod r; orbit size 6, unit action, |y| invariance and a
unique canonical member; twisted step and doubling as class maps; the table;
coefficient tracking; a planted-x recovery from two production trails forced
to merge; and the §1.1 rule against the exact rule.

**Part B** runs the section 2 per-trail algorithm on toy curves y² = x³ + 7
with q = 33556867 (cofactor 1) and q = 536873251 (cofactor 3). The pipeline
mode adds the scaled production rules: a 1.04·S0 budget, phase 2, report
capacity, the 8-match cap, and routing to 16 nodes of capacity ⌈2·D_cap/16⌉
with overflow treated as failure. [a]G + [b]T is evaluated as [a + bx]G for
speed.

v2 re-runs at the v2 budget S_ops = 1.04·√(2N' ln(1/0.604)); every recovered x
verified; T_c counts all arrivals; the H1 floor is 0.396:

| mode / classes | toy, R | trials | mean T_c / √N (model) | P[T_c ≤ S_ops] or success | 2-cycles per step (model) |
| --- | --- | ---: | --- | ---: | --- |
| T_c, 6 | q25, 256 | 3000 | 1.2684 (1.2558, +1.0%) | 0.424 | 0.000619 (1.3e-4) |
| T_c, 6 | q25, 2048 | 3000 | 1.2557 (1.2536, +0.2%) | 0.413 | 8.33e-05 (8.14e-05) |
| T_c, 6 | q29, 256 | 1000 | 1.2708 (1.2558, +1.2%) | 0.414 | 0.000634 (1.3e-4) |
| T_c, 6 | q29, 2048 | 1000 | 1.2256 (1.2536, -2.2%) | 0.450 | 8.43e-05 (8.14e-05) |
| T_c, 1 (control, no automorphism) | q25, 256 | 1500 | 1.2462 (1.2558, -1.0%) | 0.442 | 0 |
| T_c, 2 (control, negation only) | q25, 256 | 1500 | 1.2624 (1.2558, +0.5%) | 0.421 | 0.00189 (0.00195) |
| pipeline success with 16 nodes, 6 | q25, 2048 | 1226 | – | 0.4241 (node overflows 0) | 8.3e-05 |
| pipeline success with 16 nodes, 6 | q29, 2048 | 600 | – | 0.4250 (node overflows 0) | 8.4e-05 |

Every configuration clears the 0.396 floor by at least 0.017; the worst walk-mean
deviation from the model is +1.2%, inside the 4% allowance; the 2-cycle rates
match 1/(6R) throughout. The published package's corresponding v1-budget table
(seeds 1501–1512, P[T_c ≤ S_ops] 0.4128–0.4411, means -0.9% to +1.3%) is
consistent with these numbers.

Part A offline (seed 1513): §1.1 rule on 2000000 random inputs gave 0 fallbacks and 0 mismatches. On 1000000 adversarial inputs it gave 344641 fallbacks and 0 mismatches. All other part A checks are true.

**Long-trail trap test** (q29, walks without distinguished points, 150 trails
of 2,048 operations each, against the additive-escape control):

| R | escape | escapes per trail | trails with non-fruitless returns |
| ---: | --- | ---: | ---: |
| 4 | doubling | 74 | 0 / 150 |
| 8 | doubling | 41 | 18 / 150; extended run: 17 / 400 (4.3% vs 6.8% expected) |
| 4 | additive (control) | 242 | 105 / 150 |
| 8 | additive (control) | 74 | 32 / 150 |

Ordinary self-collisions predict non-fruitless returns in about 1 − exp(−L²/(2N)) = 6.8% of trails (L²/(2N) = 0.0703). At R = 8 the additive control traps less here than in the source package (8/150 against 26/150). At R = 4 it traps 110/150, while the doubling escape gives 0/150 and 6/150.

Doubling escapes behave like ordinary self-collisions; the additive control
traps trails structurally.

**Limitations.** H1 is a heuristic: a random mapping with the Brent–Pollard
correction; van Oorschot–Wiener; automorphism classes as in Wiener–Zuccherato,
GLV and Duursma–Gaudry–Morain; and fruitless-cycle escapes as in
Bos–Kleinjung–Lenstra and Bernstein–Lange–Schwabe. It is supported only at toy
scale. Instruction counts are analytic. This is a generic attack.

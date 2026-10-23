# Six-class parallel Pollard rho on secp256k1 (negation + GLV)

**Provenance.** This package is derived from submission
189ebc42-1ef3-4232-8056-a34aad574d83 (Meganpark980320,
`secp256k1-ecdlp/glv-six-class-rho-secp-v1b`), which ports the construction of
submission 1c600df4 (mpjunior92, `bls12-381-g1-ecdlp-plain/glv-six-class-rho-v11`)
to secp256k1. Both dependencies are declared. The walk, class key, 4-limb
canonical-index rule, doubling escape, 2-cycle and checkpoint rules, merging
lemma, work meter, phases, storage nodes, toy methodology and `experiment.py`
(byte-identical) are taken from 189ebc42. Only the accounting, and one
parameter, change, and every change is itemised in section 0. No score or
acceptance of either source is relied on. Everything needed is restated here.

**Claim.** time_log2 = 126.99. The algorithm recovers x from Q = [x]P on
secp256k1 with probability ≥ 0.39. On **every** run its total work is at most
the exact integer W_ledger = 168908606169243874770580353093127547782800
milli-units, i.e. 2^126.98951 units (section 5). One unit is one batched
incomplete affine secp256k1 addition (6M, A = 0) = C = 1,000 scalar AArch64
instructions, and the meter counts milli-units, 1,000 to a unit. The bound is
unconditional (global work meter, section 3). Memory < 2^72 bytes, spread over
2^16 storage nodes of < 2^56 bytes each. Preprocessing (everything before the
first walk step) ≤ 2^59.0 units, enforced by the meter and included in
W_ledger. Nonuniform advice: none. Success rests on the random-mapping
heuristic H1 (section 4), supported by toy runs.

## 0. Changes from 189ebc42, itemised

The algorithm, its parameters (θ = 2^-64, M = 2^20 walks, L_max = 20·2^64,
checkpoint period 64, S_ops, D_cap, 2^16 nodes) and H1 are unchanged. Five
accounting items change, and one parameter, the multiplier count R, with them:

1. **Inversion bound.** I < 2^17 instead of I ≤ 2^19, derived from the
   package's own M ≤ 227 (section 5). The batch-inversion share per element
   at ≥ 2^15 elements falls from 16 to 4 instructions. No other row of the
   per-step table changes.
2. **Step price.** The itemised per-step overhead is 147 instructions
   (was 159). It is charged at a budget of 150, so the ordinary step costs
   1.15 units (was 1.20 at a budget of 200). The traversal step keeps its
   +24-instruction copy allowance: 1.18 (was 1.22). The doubling price stays
   2.0. Because the margin above the table is now small, section 5 also
   derives every row from an explicit instruction sequence; those sequences
   total 137 instructions.
3. **Doubling rate, stated exactly.** 1/(6R) + 2^-23.9989 = 2^-15.5807 per
   operation. The bound used below is < 2^-15.58, rounded the safe way; the
   source's 2^-13.57, at its R = 2048, was rounded the wrong way.
4. **Meter headroom.** The claim is the ledger's own total, not a cap set above
   it: the early-halt allowance is the 0.0035 the success statement subtracts
   (was 0.0023), so success is ≥ 0.395 − 0.0035 − 2·10^-7 > 0.3914 > 0.39
   (section 4). Phase 1 costs < 2^126.98327 (was 2^127.04493), and the claim is
   2^126.99 (was 2^127.08).
5. **Peak memory restated.** The store's total is written as 2·D_cap slots
   × 104 bytes = 2^71.5 bytes. 189ebc42 quoted 2^70.6, which counts records
   rather than the 2·D_cap slots the table at load 1/2 actually has; the
   claimed ceiling of < 2^72 bytes is unchanged and holds either way. No
   other resource figure changes.
6. **Multiplier count.** R = 8192 instead of 2048. The 2-cycle rate 1/(6R)
   halves twice, 2^-13.585 to 2^-15.585, and the longer-cycle rate, every term
   of which carries R^-(L−b) with L − b ≥ 2, divides by 16, 2^-19.9954 to
   2^-23.9989. The doubling surcharge falls from 2^112.99 to 2^110.97, the
   traversal surcharge from 2^107.74 to 2^103.74, and the expected random
   extras from 2^113.02 to 2^110.98, which lowers the allowance the ledger
   charges by 2.04 bits. It costs a 4× larger table, 2^24.1 to 2^26.1, still
   inside the 2^59.0 preprocessing cap, and 1.5 MB to 6 MB of table memory,
   still inside < 2^72 bytes. Once the claim is the ledger's own total, this
   item does move it: the 2.04 bits off the allowance are 0.0292 bits of the
   total, and the ceiling follows from 127.01 to 126.99.

Gain: 127.08 − 126.99 = 0.09 bits. Both ledgers are exact: the source's, at its
own 0.0023 allowance, is 2^127.079583 units, and this package's is 2^126.989510,
a fall of 0.09007 bits in four legs. Item 2 (step 1.20 → 1.15) gives −0.05772,
of which −0.06140 is phase 1 and +0.00368 comes back through the larger 0.85
doubling surcharge; the traversal, table and node rows that items 3, 5 and 6
restate add +0.00033; item 6's R = 8192 gives −0.02943, of which −0.00026 is
phase 1 and the rest is the 2.04 bits off the allowance; and moving the
allowance target from the source's 0.0023 to the 0.0035 the success statement
subtracts gives −0.00325. The four sum to −0.09007, the gap between the two
ledgers.

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

## 2. Walk (unchanged from 189ebc42)

Parameters: R = 8192; distinguished probability θ = 2^-64 (limb 1 of K is
zero); M = 2^20 walks; L_max = 20·2^64; checkpoint period 64.

- **Table.** c_h and d_h are uniform in [0, r) (fresh coins, rejection sampling,
  every attempt metered). M_h = [c_h]P + [d_h]Q, redrawn if O.
  Tab[h][t] = [u_t^{-1}]M_h = (−1)^s·φ^{−j}(M_h), with coefficients
  u_t^{-1}·(c_h, d_h) mod r.
- **Step.** F(C) = class(c(C) + M_h) with h = K mod R. The walk keeps any member
  P and sets P ← P + Tab[h][t]. Since c + M_h = [u_t](P + Tab[h][t]), this is F.
  Each ordinary step is one batched affine addition (the reference operation)
  plus ≤ 150 instructions (section 5).
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
  re-walks, solving and verification. If a charge would exceed W_ledger
  (section 5), the attack halts with failure. So W ≤ W_ledger on every run, and
  no loop needs its own limit; rejection-sampling loops in particular are
  bounded by the meter.
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
  (bits 128..143 of K, disjoint from the index bits 0..12 and the DP bits
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
with κ = 1.05. H1 also includes three further claims:

- the partition, DP and node-routing bits of |y| are pseudo-random functions of
  the class: bits 0..12, 64..127 and 128..143, three disjoint fields;
- excluding the detected fruitless cycles of length ≤ 64, a trail makes a
  non-fruitless return with probability ≤ L_max²/(2N);
- the canonical twist index of a new class is uniform and independent of its
  partition index.

The 2-cycle rate is 1/(6R) per step. Other cycles of length ≤ 64 have rate
≤ Σ_L Σ_b S2(L, b)·R^-(L−b) = 2^-23.9989 < 2^-23.99, summed over L = 3..64
and b = 1..⌊L/2⌋, where S2(L, b) counts the partitions of the L steps into b
blocks of size ≥ 2 (every partition index used at least twice; for example
1, 1, 3 and 15 for (L, b) = (3, 1), (4, 1), (4, 2) and (6, 3)). Cycles of
length 65..200 have rate < 2^-272.85 per step. The expected doubling rate per
operation is therefore ≤ 1/(6R) + 2^-23.9989 = 2^-15.5807 < 2^-15.58.

The values are:

    S0 = √(2N' ln(1/0.605)) = 2^126.71124
    S_ops = 1.05·S0 = 2^126.78163

H1 gives ≥ 0.395; the pure model gives 0.4254. In the toy runs, class walks are
−0.9% to +1.4% from the random-mapping mean T_c, which is within the 5%
allowance (section 6).

**Phase 1 covers the budget.** Phase 1 performs ≥ S_ops operations, and each is
one class arrival. So if T_c ≤ S_ops, the first inter-trail collision happens
in phase 1.

**Doubling and traversal surcharges.** These are random, and only their excess
over the step price is not fixed by the operation count:

- doublings: expected count < 2^-15.58·(S_ops + 2^60), surcharge
  2.0 − 1.15 = 0.85 each, so expected < 2^110.97;
- traversals: expected < 0.03·64·2^-23.99·(S_ops + 2^60) < 2^103.74.

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
| meter halting early: expected random charges < 2^110.98 against the 2^119.13 allowance the claim carries, by Markov | ≤ 7/2000 = 0.0035 |
| preprocessing cap | < 2^-32 |

The events other than the early halt sum to < 2·10^-7.
**Success ≥ 0.395 − 0.0035 − 2·10^-7 > 0.3914 > 0.39.**

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
| batch inversion share (I < 2^17 per ≥ 2^15 elements) | 4 |
| **total (budget 150 = 0.15 units)** | **147** |

The step is charged at the budget, 150 instructions, which is 3 above the
itemised total. Each row is an upper bound on an explicit instruction
sequence. The constants p0 (the low limb of p), the all-ones word, 2^64 − 81,
the index value 2, the meter price and the table base stay in registers
across the batch loop; the counts below include every load and store of
per-walk state:

- F1 (27): 4 loads of x and 5 of G_1..G_5 (9), 4 MUL and 4 UMULH (8),
  7 ADD (7): 24.
- Range test, F0, S, wrap test, F2 (12): the range test is
  (F1 − 40) mod 2^64 > 2^64 − 81: SUB, CMP, branch (3). F0 is x_3, already
  loaded (0). S = F0 + F1 by ADDS, which leaves the carry c (1).
  F2 = 2^64 − S or 2^65 − S equals −(S mod 2^64) mod 2^64 in both cases: NEG,
  which keeps the flags (1). |S − 2^64| is S mod 2^64 when c = 1 and F2 when
  c = 0, so the wrap test is CSEL on c, CMP with 40, branch (3). Total 8.
- Decision (14): s1 = min(F0, F1), s2 = max(F0, F1) and their index: CMP,
  2 CSEL, CSET (4); min(s2, F2): CMP, CSEL (2); CMP F2 with s1, then CSEL
  for the minimum, the second smallest and the index (4); gap = second −
  minimum, CMP with 40, branch to the fallback (3). Total 13. A gap above 40
  to the second smallest gives a gap above 40 to both others.
- Sign, p − y, K (13): p − y by SUBS + 3 SBCS against p0 and the all-ones
  word (4); p − y < y by SUBS + 3 SBCS (4); CSET s (1); 4 CSEL for K (4).
  Total 13.
- h, distinguished test, t, address (7): AND for h (1); CBZ on limb 1 of K
  (1); t = 3s + j by two ADD (2); address = base + 64·(6h + t) by three ADD
  with shifts (3). Total 7.
- Point and state traffic (24): 8 loads of the table point, 8 loads and 8
  stores of the walk point. Total 24.
- Fingerprints (26): pending-escape flag, load and CBNZ (2). The two history
  slots form a ring indexed by length mod 2, which realises
  K_{-2} ← K_{-1}, K_{-1} ← fp(K) with one write: slot address by AND and ADD
  (2); compare fp(K) with the slot that holds K_{-2}: 3 loads, CMP and 2 CCMP,
  branch (7); overwrite that slot with fp(K), 3 stores (3). Checkpoint
  compare: 3 loads, CMP and 2 CCMP, branch (7). Checkpoint write when
  length mod 64 = 0: TST, branch (2), 3 stores (3). Total 26.
- Zero denominator (6): x(P) = x(Tab) by CMP and 3 CCMP on registers loaded
  above, branch (5).
- Length, cap, meter, loop (14): the 2-limb length (L_max > 2^64): 2 loads,
  ADDS + ADC, 2 stores (6); cap test CMP of the high limb with 20, branch (2);
  meter and operation count: SUBS on the processor's block remainder, ADD on
  its operation count, branch to the block refill, where the global counters
  and the phase-1 limit are checked (3); loop counter SUBS and branch (2).
  Total 13.
- Inversion share (4): I < 2^17 is derived below. The organizer's 6M
  convention allows 3M for the addition formula and 3M of batch amortisation
  with no standalone inversion per addition, so this 4 is charged again on
  top of a share the reference price already allows; the batch's conversion,
  scratch and arrival traffic sit in that same amortisation.

The derived sequences total 24 + 8 + 13 + 13 + 7 + 24 + 26 + 5 + 13 + 4 = 137
instructions, under the 147 in the table and the 150 charged. Rare paths are
priced as their own actions: the 2-cycle order test and the escape are in the
doubling price, the fallback has its own price, and so do reports.

**Field allowances.** These are used only outside the reference operation, and
are analytic.

- A = 32. Load 8, ADDS + 3 ADCS, carry, SUBS + 4 SBCS, 4 CSEL, store 4, ≤ 6
  overhead.
- M ≤ 227. The 4×4 column product (16 MUL/UMULH pairs, ≤ 128), four folds by
  2^32 + 977 (≤ 80), and the final conditional subtraction (≤ 19).
- I < 2^17. Fermat inversion x^(p−2) by left-to-right binary
  exponentiation. p − 2 has 256 bits and 249 one bits (recomputed with the
  package's own p), so it takes 255 squarings and 248 multiplications: 503
  multiplications (S = M) at ≤ 227 instructions each, 114,181 in total. Loop
  control is ≤ 8 instructions per exponent bit (2,048), and operand moves are
  ≤ 24 per multiplication (12,072). The total is ≤ 128,301 < 2^17 = 131,072.
  The bound does not rest on the bit count: a left-to-right binary chain
  spends 255 squarings and at most one multiplication per further bit, so at
  most 510 multiplications, and the same allowance gives
  510·251 + 2,048 = 130,058 < 2^17. Shared by a batch of ≥ 2^15 elements, the
  share is at most ⌈130,058/2^15⌉ = 4 instructions per element.

**Prices** (units):

| action | price | basis |
|---|---|---|
| ordinary step | 1.15 | reference addition plus ≤ 150 instructions |
| traversal step | 1.18 | ≤ 24 more instructions to copy the remembered member: ≤ 174 |
| doubling | 2.0 | batched affine doubling: 7M + 6A + 150 ≤ 1,931 instructions |
| §1.1 fallback | 0.5 | one M, two A, compares: < 400 instructions |
| random 64-bit word | 0.01 | |
| sampling attempt | 0.05 | |
| trail start [a]P + [b]Q | 2^13 | two Jacobian scalar multiplications at ≤ 512 operations × (16M + 16A) each, an addition and an inversion: < 4.8·10^6 instructions |
| report | 40 | route, message, node-side probes and insert: < 2^15 instructions |
| store probe | 0.05 | |
| store initialisation | 2·D_cap slots × 16 instructions = 2^58.9, once | |
| table | 8192 starts plus 5 φ-twists each: < 2^26.1, plus sampling | |
| phase-2 or re-walk operation | 600 | unbatched: an inversion (< 2^17) plus 7M, coefficient updates mod r and control, far below 600,000 instructions |
| small-batch inversion | 600 | one inversion, < 2^17 instructions = 131 units; 600 is a loose upper bound |
| solving | 600 | |
| verification | 2^13 | |

**Budget.**

- **Phase-1 operations:** 1.15·(S_ops + 2^60) < 2^126.98327.
- **Hard extras:**

  | item | units |
  |---|---|
  | meter remainders | < 2^61 |
  | phase 2: 2^20·(L_max + 65)·600 | 2^97.55 |
  | re-walks: 8·2·((L_max + 65)·600 + 2^13) | 2^81.6 |
  | reports: D_cap·40 | 2^69.2 |
  | store initialisation | 2^58.9 |
  | table and verification | < 2^26.11 |
  | **total** | **< 2^97.56** |

- **Random extras (expected):**

  | item | units |
  |---|---|
  | doubling surcharges: 0.85·2^-15.58·(S_ops + 2^60) | < 2^110.97 |
  | traversal surcharges: 0.03·64·2^-23.99·(S_ops + 2^60) | < 2^103.74 |
  | starts: < 2^62.8·2^13 | 2^75.9 |
  | words, attempts, probes, small batches, fallbacks | < 2^70 |
  | **total** | **< 2^110.98** |

- **Headroom:** the claim's total 2^126.99 less phase 1 (2^126.98327) and the
  hard extras (2^97.56) leaves 2^119.24 units for the random charges, and the
  ledger charges an allowance of exactly 2^119.13435, which fits inside it. The
  ledger's own total leaves 2^119.13435, its allowance and nothing else.
- **Early meter halt:** probability ≤ random/allowance ≤ 7/2000 = 0.0035
  (Markov), which is what the allowance is the quotient by.

**Preprocessing:** table < 2^26.1, node set-up 2^26, store initialisation
2^58.9, and the first 2^20 starts 2^33. The expected total is < 2^58.91, and the cap is 2^59.0. The
random sampling part is expected < 2^25, far below the gap, so the cap halts a
run with probability < 2^-30 by Markov. All of it is inside W_ledger.

**W ≤ W_ledger on every run.**

**Exact integer ledger and the two-decimal claim.** The meter counts
milli-units, 1,000 to a unit, and every price is a whole number of them, so
every charge and the total are integers in milli-units. S_ops is
1.05·√(2N' ln(200/121)) with ln(200/121) taken between rational bounds, so its
ceiling is exact:

    S_int = ⌈S_ops⌉   = 146242717482638971643693383952366531642
    D_cap = ⌈2^63.83⌉ = 16396253461269537204

Phase 1 charges 1,150 milli-units per operation and includes the < 2^60
reserved-block remainder. The hard extras are the table's rows plus the node
set-up of section 3, which the table does not name, each taken at the integer
its bound gives:

| hard row | milli-units |
|---|---:|
| meter remainders 2^61 | 2305843009213693952000 |
| phase 2: 2^20·(L_max + 65)·600 | 232113757366008801584480256000000 |
| re-walks: 16·((L_max + 65)·600 + 2^13) | 3541774862152233911027072000 |
| reports: D_cap·40 | 655850138450781488160000 |
| store initialisation: 2·D_cap·16 | 524680110760625190528 |
| table and verification, < 2^26.11 | 72425781000 |
| node set-up, 2^26 | 67108864000 |
| **total** | **232117299799551615389286625019528** |

That total is 2^97.55077 units, inside the table's < 2^97.56. Its largest row,
phase 2, is 2^97.5508 exactly, a hair above the table's rounded 2^97.55; the
exact integers below are what the ledger carries.

The expected random extras are the table's four rows, each rounded up, with the
cycle sum taken as the exact rational 1/(6R) + Σ_L Σ_b S2(L,b)·R^-(L-b) =
2^-15.580739 and S_ops replaced by S_int: doubling surcharges
2536433467267802166676501908348589598, traversal surcharges
16749444986391828968936680180262136, starts 65776940828046499293105000 and
words, attempts, probes, small batches and fallbacks 1180591620717411303424000,
totalling 2553182912321151528094202499125380734 milli-units = 2^110.97592
units. They are not a fourth added row. The allowance is their quotient by the
0.0035 the success statement subtracts, rounded up, so by Markov the meter
halts early with probability at most 7/2000 = 0.0035, and the allowance is the
whole of the ledger's headroom above phase 1 and the hard extras. The three
rows sum:

    phase 1
      1,150·(S_int + 2^60) = 168179125105034817391573251275519385410700
    hard extras, total below 2^97.56 units
      = 232117299799551615389286625019528
    Markov allowance, ⌈random/0.0035⌉ = 2^119.13435 units
      = 729480832091757579455486428321537352572
    W_ledger = phase 1 + hard extras + Markov allowance
      = 168908606169243874770580353093127547782800

W_ledger/1000 = 2^126.989510443… units, whose ceiling at two decimals is
126.99. The certificate is two exact integer inequalities:

    W_ledger^100 < 1000^100·2^12699    so 126.99 holds
    W_ledger^100 > 1000^100·2^12698    so 126.98 fails

Both sides carry 4,123 decimal digits and are printed in full in
`certificate.txt`; every comparison is exact integer arithmetic. The claim
126.99 sits 0.0005 bits above the ledger's own total, and the score claimed is
12699/100.

**Memory:** the store holds 2·D_cap slots of 104 bytes at load 1/2, so
2^64.83 slots × 104 bytes = 2^71.53 bytes. Adding the table (6 MB) and the
2^20 × 512 bytes of walk states gives < 2^72 bytes, with each of
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
mode adds the scaled production rules: a 1.05·S0 budget, phase 2, report
capacity, the 8-match cap, and routing to 16 nodes of capacity ⌈2·D_cap/16⌉
with overflow treated as failure. [a]G + [b]T is evaluated as [a + bx]G for
speed.

Offline results (seeds 1501–1512; every recovered x verified; S_ops = 1.05·√(2N' ln(1/0.605)) for each toy; T_c counts all arrivals; every R = 2048 row has an R = 8192 twin at the same seed and protocol):

| mode / classes | toy, R | trials | mean T_c / √N (model) | P[T_c ≤ S_ops] or success | 2-cycles per step (model) |
| --- | --- | ---: | --- | ---: | --- |
| T_c, 6 | q25, 256 | 6000 | 1.2689 (1.2558, +1.0%) | 0.418 | 0.00063 |
| T_c, 6 | q25, 2048 | 6000 | 1.2696 (1.2536, +1.3%) | 0.4185 | 8.23e-05 |
| T_c, 6 | q25, 8192 | 6000 | 1.2638 (1.2534, +0.8%) | 0.416 | 1.99e-05 |
| T_c, 6 | q29, 256 | 1800 | 1.2441 (1.2558, -0.9%) | 0.4411 | 0.000632 |
| T_c, 6 | q29, 2048 | 1800 | 1.2699 (1.2536, +1.3%) | 0.4128 | 8.28e-05 |
| T_c, 6 | q29, 8192 | 1800 | 1.2707 (1.2534, +1.4%) | 0.4133 | 2.01e-05 |
| T_c, 1 (control, no automorphism) | q25, 256 | 3000 | 1.2435 (1.2558, -1.0%) | 0.429 | 0 |
| T_c, 2 (control, negation only) | q25, 256 | 3000 | 1.2741 (1.2558, +1.5%) | 0.4117 | 0.0019 |
| pipeline success with 16 nodes, 6 | q25, 2048 | 12000 | – | 0.4229 (node overflows 0, max load/cap 0.72) | 7.86e-05 |
| pipeline success with 16 nodes, 6 | q25, 8192 | 12000 | – | 0.4275 (node overflows 0, max load/cap 0.68) | 1.98e-05 |
| pipeline success with 16 nodes, 6 | q29, 2048 | 2400 | – | 0.4238 (node overflows 0, max load/cap 0.52) | 7.72e-05 |
| pipeline success with 16 nodes, 6 | q29, 8192 | 2400 | – | 0.4333 (node overflows 0, max load/cap 0.62) | 2.06e-05 |

Part A offline (seed 1513): §1.1 rule on 2000000 random inputs gave 0 fallbacks and 0 mismatches. On 1000000 adversarial inputs it gave 344641 fallbacks and 0 mismatches. All other part A checks are true.

**Long-trail trap test** (q29, walks without distinguished points, 150 trails
of 2,048 operations each, against the additive-escape control):

| R | escape | escapes per trail | trails with non-fruitless returns |
| ---: | --- | ---: | ---: |
| 4 | doubling | 72 | 0 / 150 |
| 8 | doubling | 42 | 6 / 150 |
| 4 | additive (control) | 223 | 110 / 150 |
| 8 | additive (control) | 64 | 8 / 150 |

Ordinary self-collisions predict non-fruitless returns in about 1 − exp(−L²/(2N)) = 6.8% of trails (L²/(2N) = 0.0703). At R = 8 the additive control traps less here than in 1c600df4's BLS12-381 run (8/150 against 26/150). At R = 4 it traps 110/150, while the doubling escape gives 0/150 and 6/150.

Doubling escapes behave like ordinary self-collisions; the additive control
traps trails structurally.

**Doubling rate against the toys.** The toy 2-cycle rates at R = 8192
(1.99·10^-5, 2.01·10^-5, 1.98·10^-5 and 2.06·10^-5 per step in the T_c and
pipeline rows) agree with 1/(6R) = 2.03·10^-5 within one standard error; each
row has 2·10^2 to 5·10^2 events, a relative error of 5% to 7%. At R = 2048 the
same four rows gave 8.23·10^-5, 8.28·10^-5, 7.86·10^-5 and 7.72·10^-5 against
1/(6R) = 8.14·10^-5, so the measured rate follows 1/R across the factor 4 in
R. The toys also show checkpoint escapes at 1.0·10^-6 to 4.6·10^-6 per step,
above the structural rate 2^-23.99. These are ordinary self-collisions of
short toy trails: a trail of length ℓ meets its own past with probability about
ℓ/N per step, about 32/(5.6·10^6) ≈ 6·10^-6 at q25 with dp_bits = 5. In
production ℓ/N < 2^70/2^253, so the structural rate governs. The allowance is
the claim's whole margin, so a bound above the random charges costs it directly:
in the q25, R = 8192 row (seed 1502) the total doubling rate, both effects
included, is 2.45·10^-5 per step, 19% above the production bound
2^-15.5807 = 2.05·10^-5, because the toy's self-collision term (4.6·10^-6 per
step, independent of R) is now a large share of the smaller 2-cycle rate; at
R = 2048 the same term was 6% of that bound. Production has no such term, but if
its random charges were 19% above their 2^110.98 bound the early-halt
probability would be 1.19·0.0035 = 0.0042 and success would still be
≥ 0.395 − 0.0042 − 2·10^-7 > 0.3906 > 0.39. A 45% excess would give 0.0051 and
0.3899, which this claim does not carry: the ledger spends the 0.0035 the
success statement subtracts and no more.

**Provenance of the rows.** `experiment.py` is byte-identical to the file in
189ebc42, and the rows above are its offline results for the configurations in
README.md. Python's seeded generator makes each row deterministic, so each
row can be rerun from its configuration. The R = 2048 rows are 189ebc42's,
reproduced with this file for the previous revision of this package (seeds
1501, 1502 and 1507 gave 1.2689 with 0.418, 1.2696 with 0.4185, and pipeline
success 0.4229 with 0 node overflows and max load/cap 0.72). The four R = 8192
rows are new runs of the same file, one per R = 2048 row, at that row's seed
and protocol, with R the only input that changed.

**Initial input.** The hosted configuration (q25, R = 8192, dp_bits = 5,
seed 1, pipeline, 14 s) is time-limited, so the number of completed trials
depends on the machine. One local run completed 287 trials: part A all true,
0/30,000 and 0/10,000 index-rule mismatches, pipeline success 0.4599 with a
binomial standard error of 0.0288, and 0 node overflows. This short run is
consistent with the 12,000-trial pipeline row (0.4275); the H1 floor is
taken from the large offline rows, not from the 14-second run.

**Limitations.** H1 is a heuristic: a random mapping with the Brent–Pollard
correction; van Oorschot–Wiener; automorphism classes as in Wiener–Zuccherato,
GLV and Duursma–Gaudry–Morain; and fruitless-cycle escapes as in
Bos–Kleinjung–Lenstra and Bernstein–Lange–Schwabe. It is supported only at toy
scale. Instruction counts are analytic. This is a generic attack.

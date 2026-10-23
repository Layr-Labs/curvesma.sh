# Six-class parallel Pollard rho on secp256k1 (negation + GLV)

**Provenance.** This package is derived from submission
189ebc42-1ef3-4232-8056-a34aad574d83 (Meganpark980320,
`secp256k1-ecdlp/glv-six-class-rho-secp-v1b`), which ports the construction of
submission 1c600df4 (mpjunior92, `bls12-381-g1-ecdlp-plain/glv-six-class-rho-v11`)
to secp256k1. Both dependencies are declared. The walk, class key, 4-limb
canonical-index rule, doubling escape, 2-cycle and checkpoint rules, merging
lemma, work meter, phases, storage-node model, toy methodology and `experiment.py`
(byte-identical) are taken from 189ebc42. Only the accounting and the table size
change, and every change is itemised in section 0. No score or acceptance of
either source is relied on. Everything needed is restated here.

**Claim.** time_log2 = 126.99. The algorithm recovers x from Q = [x]P on
secp256k1 with probability ≥ 0.39. On **every** run its total work is
W ≤ 2^126.99 units. One unit is one batched incomplete affine secp256k1 addition
(6M, A = 0) = C = 1,000 scalar AArch64 instructions. The bound is unconditional
(global work meter, section 3). Memory < 2^72 bytes, spread over 2^16 storage
nodes of < 2^56 bytes each. Preprocessing (everything before the first walk
step) ≤ 2^59.0 units, enforced by the meter and included in W. Nonuniform
advice: none. Success rests on the random-mapping heuristic H1 (section 4),
supported by toy runs.

## 0. Changes from 189ebc42, itemised

The algorithm, its class key, the doubling escape mechanism, the 2-cycle and
checkpoint rules, the work meter, the storage-node model and H1 are unchanged.
Seven accounting and parameter items change:

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
3. **Doubling rate, stated exactly.** With R = 8192 the 2-cycle rate is
   1/(6R) = 2^-15.58496 per operation. Cycles of length 3..64 that need every
   partition index at least twice have rate 2^-23.99885 (the sum over
   L = 3..64 and b = 1..⌊L/2⌋ of S2(L, b)·R^-(L−b)). The exact total is
   2^-15.58074 per operation, and the ledger uses the bound < 2^-15.58.
4. **Meter headroom, stated from the ledger.** Expected random charges
   < 2^110.98 against headroom > 2^119.24 give an early-halt probability
   below 2^-8.26 = 0.0033 < 0.0035 by Markov, so success is
   ≥ 0.395 − 0.0033 − 2·10^-7 > 0.3917 > 0.39. Phase 1 costs
   < 2^126.98327, and the claim is 2^126.99.
5. **Peak memory restated.** The store's total is 2·D_cap slots × 104 bytes =
   2^71.53 bytes, plus the probe index arrays of 8 bytes per slot, 2^67.83
   bytes. 189ebc42 quoted 2^70.6, which counts records rather than slots and
   has no probe arrays; the ceiling < 2^72 bytes is unchanged and holds.
6. **Distinguished-store probe derived.** The probe depth of the report store
   is no longer asserted. Each node holds B = ⌈2·D_cap/2^16⌉ =
   500,373,945,961,595 buckets of 104 bytes, and a report probes at most eight
   buckets in a fixed linear sequence; the probe is itemised at 39
   instructions and charged 50 (0.05 units). A report whose eight buckets are
   occupied without a key match enters a charged per-node overflow list of
   8-byte bucket addresses, one per bucket, so the list can never hold more
   than the node's slot count. The worst-case overflow insert is a scan of B
   buckets (6 instructions each) plus a 20-instruction entry write, and the
   worst-case overflow lookup is 16 instructions per entry over B entries;
   both are charged, D_cap·(22B + 20) instructions = 2^107.154 units in the
   hard extras. A report now costs 40.4 units: 40 for routing, message, insert
   and node-side tests, plus its 8 probes at 0.05.
7. **R = 8192 instead of 2048.** The table is 8192 starts, < 2^26.1 units
   instead of < 2^24.1. N' = N·8192/8191, and S_ops becomes 2^126.78163.
   The doubling surcharge falls from 2^112.98 to 2^110.97: 1/(6R) falls by a
   factor 4 and the other-cycle rate by a factor 16, so the surcharge halves
   twice. The expected random total falls from 2^113.02 to 2^110.98, and the
   required log2 W from 127.00905 to 126.98951.

Gain: 127.08 − 126.99 = 0.09 bits. Items 1, 2 and 4 account for the move from
127.08 to 127.01 (item 2 gives 0.0614 bits on phase 1, item 4 the rest).
Items 6 and 7 account for the move from 127.01 to 126.99: item 7 drops the
required log2 W by 0.01954 (0.00026 from S_ops and 0.01927 from the random
total), item 6 by below 10^-5, and the 2-decimal ceiling by the remaining
0.00046.

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

## 2. Walk (unchanged from 189ebc42 except R)

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
  ≤ 64 ordinary steps), i.e., only on a cycle Z of F of length ≤ 64. Z then
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
  fallbacks, reports, store probes, overflow inserts and lookups, table build,
  store initialisation, re-walks, solving and verification. If a charge would
  exceed W_max = 2^126.99, the attack halts with failure. So W ≤ 2^126.99 on
  every run, and no loop needs its own limit; rejection-sampling loops in
  particular are bounded by the meter.
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
  is D_cap = 2^63.83 records, taken as the integer
  ⌈2^63.83⌉ = 16,396,253,461,269,537,204 in every integer product below. These
  bound memory and re-walk work.
- **Batches.** Ready walks are split into batches of 2^15..2^16 per Montgomery
  inversion. A batch of fewer than 2^15 ready walks pays its inversion as a
  separate 600-unit action.
- **Store.** Reports are routed by bits 128..143 of K (limb 2, disjoint from the
  index bits 0..10, the DP bits 64..127 and the bucket field 144..193) to one
  of 2^16 storage nodes. Each node is a separate machine with its own 64-bit
  address space. It holds an open-addressing table of
  B = ⌈2·D_cap/2^16⌉ = 500,373,945,961,595 buckets of 104 bytes (32-byte key K,
  32-byte a, 32-byte b, 8-byte tag and flags), < 2^55.54 bytes, and an
  overflow list of B 8-byte bucket addresses, < 2^51.84 bytes; with both, a
  node's store is 112 bytes per bucket, < 2^55.64 bytes.
  - **Home bucket.** The 50-bit field at bits 144..193 of K, reduced modulo B
    by a 64-bit multiply-shift (6 instructions, inside the report price).
    Bits 144..193 are disjoint from the index, DP and routing fields.
  - **Probe rule.** Lookup and insertion use linear probing with a fixed
    maximum depth of 8 buckets. Probe i examines bucket j = (home + i) mod B:
    - if the bucket is free, the report is inserted there;
    - if the tag and the four key words match, the report is a match (the
      collision the attack wants);
    - after 8 occupied non-matching buckets, the report enters the overflow
      list: the list is scanned for a match, then a free bucket is taken by
      scanning forward from home + 8, and its 8-byte address is appended. No
      report probes more than 8 buckets of the table.
  - **Probe itemisation** (39 instructions; charged 50 = 0.05 units):

    | item | instructions |
    |---|---:|
    | j = (home + i) mod B: ADD, CMP, CSEL | 3 |
    | 13·j: LSL #3, LSL #2, ADD, ADD; the loads use [base, X, LSL #3] | 4 |
    | flags word LDR, CBZ | 2 |
    | tag LDR, CMP, B.NE | 3 |
    | key: 4 LDR | 4 |
    | key: 4 CMP/CCMP, B.EQ | 5 |
    | record head LDR; visit stamp LDR, ORR, STR | 4 |
    | next probe: ADD, CMP, B | 3 |
    | prefetch: next-bucket 13·j' (2 LSL, ADD), PRFM | 4 |
    | two branch-shadow slots and loop alignment | 3 |
    | probe-loop register reloads (tag and key base) | 4 |
    | **total (charged 50)** | **39** |

    The 11-instruction margin covers the two conditional branches and the
    loop's register moves; this is the same budget-above-table convention as
    the step, and every charge here is an upper bound.
  - **Overflow list and its capacity.** The list is a per-node array of B
    8-byte bucket addresses with a counter, and every overflowing report is
    written into one distinct free bucket of the table with one address
    appended for it. The table has B buckets, so the list can never hold more
    than B entries and cannot overflow; no probability enters its capacity.
    If the forward scan reaches the end of the table without finding a free
    bucket, the node is full: that is the node-overflow event already bounded
    in section 4, and the attack halts with failure.
  - **Overflow insert, worst case.** Scan from home + 8 to the first free
    bucket: 6 instructions per bucket (flags LDR, CBNZ, scaled address ADD,
    index ADD, CMP, B.LO); write the 104-byte entry (13 STR plus the address,
    16); append the address and increment the counter (4). At most
    6B + 20 = 3,002,243,675,769,590 instructions = 2^41.45 units.
  - **Overflow lookup, worst case.** Scan the list: 16 instructions per entry
    (address LDR, tag LDR, CMP, B.NE, four key LDR, CMP and 3 CCMP, B.EQ,
    pointer ADD, end CMP, B). At most 16B = 8,005,983,135,385,520
    instructions = 2^42.86 units.
  - **Charged worst case.** Each report performs at most one lookup and one
    insert and probes at most 8 table buckets, so the probe and overflow work
    over all reports is at most D_cap·(22B + 20) instructions = 2^107.154
    units. That is the ledger's overflow row; the 8 probes are the 0.4 units
    of the report price.
  - **Load.** Under H1 the routing field is pseudo-random, so a node's share
    of the n ≤ D_cap reports is Binomial(n, 2^-16) with mean < 2^47.83
    against B = 2^48.83 buckets. Chernoff gives < exp(−2^45) per node, and
    < exp(−2^44) over all 2^16 nodes.
  - **Node set-up** (process start, table allocation, routing tables) is
    charged 2^20 instructions per node, 2^26 units in total, as
    preprocessing. **Slot initialisation** is charged separately below at 16
    instructions per bucket; the overflow array needs only its counter
    initialised (2 instructions per node), because every entry is written
    before it is read.

## 4. Success probability

N' = N·R/(R − 1) = N·8192/8191, log2 N' = 253.41521. T_c is the number of
class arrivals of all trails, up to the first arrival at a class already
visited by another trail. It includes trail starts, steps, doublings,
traversal steps and repeated arrivals inside fruitless cycles; nothing is
excluded. This is exactly what the toy measures.

**H1 (heuristic, at the budget).** P(T_c ≤ S_ops) ≥ 1 − exp(−(S_ops/κ)²/(2N'))
with κ = 1.05. H1 also includes three further claims:

- the partition, DP, routing and bucket fields of |y| are pseudo-random
  functions of the class: bits 0..10, 64..127, 128..143 and 144..193, four
  disjoint fields;
- excluding the detected fruitless cycles of length ≤ 64, a trail makes a
  non-fruitless return with probability ≤ L_max²/(2N);
- the canonical twist index of a new class is uniform and independent of its
  partition index.

The 2-cycle rate is 1/(6R) = 2^-15.58496 per step. Other cycles of length ≤ 64
have rate ≤ Σ_L Σ_b S2(L, b)·R^-(L−b) = 2^-23.99885, summed over L = 3..64
and b = 1..⌊L/2⌋, where S2(L, b) counts the partitions of the L steps into b
blocks of size ≥ 2 (every partition index used at least twice; for example
1, 1, 3 and 15 for (L, b) = (3, 1), (4, 1), (4, 2) and (6, 3)). Cycles of
length 65..200 have rate < 2^-206 per step. The expected doubling rate per
operation is therefore ≤ 1/(6R) + 2^-23.99885 = 2^-15.58074 < 2^-15.58.

The values are:

    S0 = √(2N' ln(1/0.605)) = 2^126.71124
    S_ops = 1.05·S0 = 2^126.78163

H1 gives ≥ 0.395; the pure model gives 0.4254. In the toy runs, class walks are
−1.0% to +1.5% from the random-mapping mean T_c, which is within the 5%
allowance (section 6).

**Phase 1 covers the budget.** Phase 1 performs ≥ S_ops operations, and each is
one class arrival. So if T_c ≤ S_ops, the first inter-trail collision happens
in phase 1.

**Doubling and traversal surcharges.** These are random, and only their excess
over the step price is not fixed by the operation count:

- doublings: expected count < 2^-15.58·(S_ops + 2^60), surcharge
  2.0 − 1.15 = 0.85 each, so expected < 2^110.97;
- traversals: expected < 0.03·64·2^-23.99885·(S_ops + 2^60) < 2^103.73.

**Failure events** (union bound): if T_c ≤ S_ops, the merging lemma gives
success unless one of these occurs.

| event | bound |
|---|---|
| trail cap | < 10^-7 |
| non-fruitless return in a colliding trail | < 2^-50 |
| report overflow (Chernoff) | < exp(−2^59) |
| node overflow (Chernoff, routing bits under H1) | < exp(−2^44) |
| overflow list capacity (one address per table bucket) | 0, exact |
| f_1 = f_2 | 1/r |
| exceptional additions | < 2^-114 |
| 8 failing matches | < 2^-60 |
| spurious 192-bit fingerprint match (over < 2^128 comparisons) | < 2^-64 |
| meter halting early: expected random charges < 2^110.98 against headroom > 2^119.24, by Markov | < 0.0033 |
| preprocessing cap | < 2^-32 |

The events other than the early halt sum to < 2·10^-7.
**Success ≥ 0.395 − 0.0033 − 2·10^-7 > 0.3917 > 0.39.**

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
| store probe | 0.05 | 39 instructions itemised in section 3, charged 50 |
| report | 40.4 | 40 for route, message, entry write and node tests (< 2^15 instructions), plus 8 probes at 0.05 |
| overflow insert scan | 2^41.45 | worst case 6B + 20 instructions |
| overflow lookup | 2^42.86 | worst case 16B instructions over the B-entry list |
| store initialisation | 2^16·B buckets × 16 instructions = 2^58.864, once | |
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
  | phase 2: 2^20·(L_max + 65)·600 | 2^97.551 |
  | re-walks: 8·2·((L_max + 65)·600 + 2^13) | 2^81.551 |
  | reports: D_cap·40.4 | 2^69.166 |
  | store initialisation: 2^16·B × 16 instructions | 2^58.864 |
  | table and verification | < 2^27 |
  | overflow insert and lookup: D_cap·(22B + 20) instructions | 2^107.154 |
  | **total** | **< 2^107.156** |

- **Random extras (expected):**

  | item | units |
  |---|---|
  | doubling surcharges: 0.85·2^-15.58·(S_ops + 2^60) | < 2^110.97 |
  | traversal surcharges: 0.03·64·2^-23.99885·(S_ops + 2^60) | < 2^103.73 |
  | starts: < 2^62.8·2^13 | 2^76 |
  | words, attempts, small batches, fallbacks | < 2^70 |
  | **total** | **< 2^110.98** |

- **Headroom:** 2^126.99 − 2^126.98327 − 2^107.156 > 2^119.24.
- **Early meter halt:** probability < 2^110.98/2^119.24 = 2^-8.26 = 0.0033
  < 0.0035 (Markov).

### 5.1 Exact total and the two-decimal claim

All ledger quantities are carried as integers in milli-units (10^-3 units; one
instruction is one milli-unit), so every product is exact. With
S_int = ⌈S_ops⌉ = 146,242,717,482,638,971,643,693,383,952,366,531,642,
D_cap = 16,396,253,461,269,537,204 and B = 500,373,945,961,595:

    phase 1            1150·(S_int + 2^60)                       = 168179125105034817391573251275519385410700
    hard extras total                                               =    180725794254648202537786083482944760
    Markov allowance   ceil(2000·random/7)                          =    729480832094552128854877235214689066858
    W (milli-units)                                                 = 168908786662923624168630666296817557422318

W in units is W/1000 = 2^126.989511985…, whose 2-decimal ceiling is 126.99.
The exact integer certificate is W^100 < 1000^100·2^12699, and W^100 ≥
1000^100·2^12698, so 126.98 does not suffice. The two sides of the first
inequality are the integers

    W^100 =
    5824038308042708946655236209551943059141914110952207721513776716734276812186500162026530238871172636070255750998630167381732389071122033831478115989421854568098308791316974711583487662609479138546764227556884257286300740044039543007734008653228648698806247656906867703690890643680886497295596523870547960138919521653661300843569473753667230010962689055013228784890513323998030546139019656111369133660725434414468124873438170920380792256553880508155334059028691391652856806831703183222876773307934511116067235127522733037714039392202213370553654931783085159795609342133238898227936723925902353113795181046597317456684117711123566057994161474118103913893809010594539597429872146495806113391615769303856004596220717350925241926690514875846330570098267190158283054099316021268500015472157764843720983897245419531973472098773716286567392294884299314189012891819572971063766887898498115142387609745427395530769472339284817125480890748802612359887861683409976725727644731978325127969373140742889429594048731721655678809641788338226860502558292743277872925276458291359928218763284617212084318068737241635699765257936559249697545835197475681981903873561480302130376897764671813982722246268531624633398704130816719008867690822654250202279116647450517472218194118637567912842699043475521160303737663063739630580654134071853377423351403445202328941066555100443497526699003089451041463110503200932077027223163918497796510273073297492235737648547910163013107269713219287215212310731415627628198236805852377638905536138808162430092644135769461880778877536613664167196487338794460238622950242410405684995770034183218120203277276165302586780156891562624424673347246367965061654138503554200657471344456845874153889170266548770430200609167845863155796072597046000441518724154594422267887127498429997912369451052301696031612005948349110631506827250755438474309911238914103354933423522431743592460102340398010753867517216280515481198997201680332444336073832089027407211219292212008668373847112226165053734907306507233867445311163221536941487951475793373475782572654374318756759082494370903287549161542911946673900411589419417181305138154757466809895642235094412897683394713873868581317261434032252058224546236378825861824217182861635893221197427321719114239667516989093260763437195948133362730510801954844222408564448265127946484320563002733203265541716116927186615200614574226855198477299981713317410546111786495694990685150420545136096847695722228502515882156921000020860947412184520000546963377459885391170102883713491398808576460008081050485467451805891241823836698413195647649525234886482485507255875365779167250795875160674617523937454530911370913266358286103156170344838536645699149673922702466362266323990837751450061899508427337183714407707642777317757625737129996923228048829848922369945741978634662008029523661867261692517781325910758885518041648628023251364008796127088589819272040306126496059841151353063315991832682270273500751616123754524715899671399294339375780838372116459719915837748420956301025048505241939508020784780255103236519027876220527425671613402174354386971982197741145266662310903745698870608691588151125232122367675012625277391290456466532259763814307676431181840071578585425956722165208074471459716302463708571015058153268868204606509292633200167678647234270634391653169541770669597176187211531590761084756316998295257095201031156583262403309559075706750930538621064854007495223882930169483624589690874593122464256485993766439764204406624077438267677783217673014668971044761669736716245689897624450193895608656364842224221503217492170379419687217888045663710812287194352189946279924289156113714728722563951703919621206859625918090085939526451893318108526219404939318009791796151111650364276058203697668725539438197176438779405631358378320550898910316003579226516361801001460832427752562758606126329049071836781020130599727654482951700406601380436125426611605738619455660431902248086696462305612715848137491466286754029604791395356225804962027752879847682072285167283220763036439151289901176127654398391468705718420840224944395870689169462937571764578146896140964223044346458331501836068054201091990832871035109376

    1000^100·2^12699 =
    6024415772787383307035263421991409507394100025577078750888849087302161315652681059004917043893425580687899516763429780513811860607698849681991834776658053710020906974283663622451003449884644244109346927376022531887799374119399437191972083748355162724822329218748773680488268363412578544508105521037866037320903329695192322070402398926727547993012284122329439836671825844847128732202236015375223017320083592724891954658814614895102112626542292201025000217606073910934876557395139884494631420449327256508162260783043770010968883889666206601504581164799904458286145057107799743457963874514594637181991133441155319586694452304435703778896731547065556291976875946490081082932872278206505357635419400236022711664180002972280132820906822483247298043372619793443956363917622834812794159119350494087143346683130181119346778638874328727574550754499783064002978086137397500621332742322698577706052796202500501233688660556429201471232530114657092748640787099221064761914151272193951344889776397272607222788774468711375628045211429045753037034465605303674675814981867770721753677880927693060059990134281712307539896863892025787267433096911754092159717913705000601097841321914248647294346050962861257835071194999699402227629547138605867991989499981357340222147747191015819247079425669323159039952880689417543884585951220106125026841938859080714598406677703236973930078926788809823239118155724926359246126433684496884436999081750402746426738198775978743258024258387455684414904768952161046070662585061525490136442307978099685480396433002631451989474135227363541203478330028348702211867593641330253648178963277453622100716362713221294434465996503924121938877772427964467175224910194195674915741228881865172973767977040807778789639100800101808798041050470875516465631910767793254146487974686271514558139353295007475569870854233459730341194816392173136122561674468404541478154173453051129810108391776570305082255386721201983098338287801335566857254269146961245143540332594789985632185120213930946477346527344202020041420911628700285632812065407272198352890727843089290141306937963172771866536633688362910351963858666863814805286810010809519645169909623499768469015933841810677724797245017795292187024018281000139682728652458765835022634068066365296421173635015527838164314081140621832378136574194103179209847701120075405293315177464690639358026942230205551353326987478422493268714428243980279723968889545490926271445822701230045656823565949765913096518658916895001035847585273372274089314449253009561677747825116531544967924234893510758840129911255824322031659415157326401252255806883045024117433028200534116576473927709948068283617363632320812306990012177229439305614862587750713758989938618283073875978977350032568421090212769106401130685202611412868375639972357808207025992823123815776449513069053446406020326023778771407982484361535569792776717378226040854793815535329789198514755108937034475478483108222830599143851998100917253982874290638376400419520002699137165526721253626312886211079397566551893045091246787442567387406084388054583331288945671744102633391674022087147635795991419137045398180327136261607082756916110046192761858130776002353729301083193727977200795192683365104710352630969282404130336420732715230380949406046497279921291805309609813736616923064927260293511584529529735004558259616169975254275123310789705162334706005647340154761210031257840194247439278619372867734281573877720142816981660104568391183527691013082422003810661543202837733053708301713945184546354357854100026679382200220014036475804824662154246606318109983707927404293263529267430110586614628365896737864193146113627664768413645043493060120061819211451772502732230582239520781926949753993780352155487340241021403834396938435430491078767071252092460587515767512577929746349371746060321153578628142617223356121284316767866380208332007428652686692673650688000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

Both carry 4,123 decimal digits; the comparison is exact integer
arithmetic, as is W^100 ≥ 1000^100·2^12698 (the right side is B/2, and the
ratio W/(2^126.99·1000) is 2^-0.000488). The score claimed is
12699/100.

**W ≤ 2^126.99 on every run.**

**Preprocessing:** table < 2^26.1, node set-up 2^26, store initialisation
2^58.864, and the first 2^20 starts 2^33. The expected total is < 2^58.87,
and the cap is 2^59.0. The random sampling part is expected < 2^25, far below
the gap, so the cap halts a run with probability < 2^-30 by Markov. All of it
is inside W.

**Memory:** the store holds 2·D_cap slots of 104 bytes at load 1/2, so
2^64.83 slots × 104 bytes = 2^71.53 bytes, and the probe index arrays add
8 bytes per slot, 2^67.83 bytes: 2^71.64 bytes in total. Adding the table
(5.1 MB) and the 2^20 × 512 bytes of walk states gives < 2^72 bytes, with
each of the 2^16 nodes at 112 bytes per bucket, 2^55.64 bytes < 2^56.

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

Offline results (seeds 1501–1512; every recovered x verified;
S_ops = 1.05·√(2N' ln(1/0.605)) for each toy; T_c counts all arrivals):

| mode / classes | toy, R | trials | mean T_c / √N (model) | P[T_c ≤ S_ops] or success | 2-cycles per step (model) |
| --- | --- | ---: | --- | ---: | --- |
| T_c, 6 | q25, 256 | 6000 | 1.2689 (1.2558, +1.0%) | 0.418 | 0.00063 |
| T_c, 6 | q25, 2048 | 6000 | 1.2696 (1.2536, +1.3%) | 0.4185 | 8.23e-05 |
| T_c, 6 | q29, 256 | 1800 | 1.2441 (1.2558, -0.9%) | 0.4411 | 0.000632 |
| T_c, 6 | q29, 2048 | 1800 | 1.2699 (1.2536, +1.3%) | 0.4128 | 8.28e-05 |
| T_c, 6 | q25, 8192 | 6000 | 1.2638 (1.2534, +0.8%) | 0.416 | 1.99e-05 |
| T_c, 6 | q29, 8192 | 1800 | 1.2707 (1.2534, +1.4%) | 0.4133 | 2.01e-05 |
| T_c, 1 (control, no automorphism) | q25, 256 | 3000 | 1.2435 (1.2558, -1.0%) | 0.429 | 0 |
| T_c, 2 (control, negation only) | q25, 256 | 3000 | 1.2741 (1.2558, +1.5%) | 0.4117 | 0.0019 |
| pipeline success with 16 nodes, 6 | q25, 2048 | 12000 | – | 0.4229 (node overflows 0, max load/cap 0.72) | 7.86e-05 |
| pipeline success with 16 nodes, 6 | q29, 2048 | 2400 | – | 0.4238 (node overflows 0, max load/cap 0.52) | 7.72e-05 |
| pipeline success with 16 nodes, 6 | q25, 8192 | 12000 | – | 0.4275 (node overflows 0, max load/cap 0.68) | 1.98e-05 |
| pipeline success with 16 nodes, 6 | q29, 8192 | 2400 | – | 0.4333 (node overflows 0, max load/cap 0.62) | 2.06e-05 |

Part A offline (seed 1513): §1.1 rule on 2000000 random inputs gave 0
fallbacks and 0 mismatches. On 1000000 adversarial inputs it gave 344641
fallbacks and 0 mismatches. All other part A checks are true.

**Long-trail trap test** (q29, walks without distinguished points, 150 trails
of 2,048 operations each, against the additive-escape control):

| R | escape | escapes per trail | trails with non-fruitless returns |
| ---: | --- | ---: | ---: |
| 4 | doubling | 72 | 0 / 150 |
| 8 | doubling | 42 | 6 / 150 |
| 4 | additive (control) | 223 | 110 / 150 |
| 8 | additive (control) | 64 | 8 / 150 |

Ordinary self-collisions predict non-fruitless returns in about
1 − exp(−L²/(2N)) = 6.8% of trails (L²/(2N) = 0.0703). At R = 8 the additive
control traps less here than in 1c600df4's BLS12-381 run (8/150 against
26/150). At R = 4 it traps 110/150, while the doubling escape gives 0/150 and
6/150.

Doubling escapes behave like ordinary self-collisions; the additive control
traps trails structurally.

**Doubling rate against the toys.** The toy 2-cycle rates at R = 8192 are
1.99·10^-5 (q25 T_c), 2.01·10^-5 (q29 T_c), 1.98·10^-5 (q25 pipeline) and
2.06·10^-5 (q29 pipeline), against 1/(6R) = 2.03·10^-5. The total doubling
rate including checkpoint escapes is 2.45·10^-5 and 2.36·10^-5 at q25, and
2.11·10^-5 and 2.15·10^-5 at q29: 16% to 20% above the production bound at
q25 and 4% to 6% at q29. The q25 excess is ordinary self-collision of short
trails: a trail of length ℓ ≈ 32 meets its own past with probability about
ℓ/N per step, 6·10^-6 at q25, and the measured checkpoint rate is 3.8·10^-6
to 4.6·10^-6 per step; at q29 it is 9.9·10^-7. In production ℓ/N < 2^70/2^253,
so the structural rate 2^-23.99885 governs. The R = 2048 rows show the same
effect: 8.2·10^-5 in the T_c rows against their own 1/(6R) = 8.14·10^-5, and
8.72·10^-5 total at q25 with checkpoint escapes near 5·10^-6. The early-halt
margin covers it: with all random charges 45% above their bound 2^110.98,
Markov gives 1.45·2^-8.26 = 0.0047, and success is ≥ 0.395 − 0.0047 −
2·10^-7 = 0.3903 > 0.39. Charging the measured excess instead
(1.2·2^-8.26 = 0.0039) gives success ≥ 0.3911.

**Provenance of the rows.** `experiment.py` is byte-identical to the file in
189ebc42, and the rows above are its results for the configurations in
README.md. Python's seeded generator makes each row deterministic, so each row
can be rerun from its configuration. For this package the four R = 8192 rows
were run with seeds 1502, 1504, 1507 and 1508 under the same protocols as the
R = 2048 rows they extend, and the R = 2048 rows of seeds 1501, 1502 and 1507
were also rerun and reproduced exactly (1.2689 and 0.418; 1.2696 and 0.4185;
pipeline success 0.4229 with 0 node overflows and max load/cap 0.72).

**Initial input.** The hosted configuration (q25, R = 8192, dp_bits = 5,
seed 1, pipeline, 14 s) is time-limited, so the number of completed trials
depends on the machine. One local run completed 109 trials: part A all true,
0/30,000 and 0/10,000 index-rule mismatches, pipeline success 0.4404 with a
binomial standard error of 0.047, and 0 node overflows (maximum load 0.52 of
capacity). This short run is consistent with the 12,000-trial pipeline row
(0.4275); the H1 floor is taken from the large offline rows, not from the
14-second run.

**Limitations.** H1 is a heuristic: a random mapping with the Brent–Pollard
correction; van Oorschot–Wiener; automorphism classes as in Wiener–Zuccherato,
GLV and Duursma–Gaudry–Morain; and fruitless-cycle escapes as in
Bos–Kleinjung–Lenstra and Bernstein–Lange–Schwabe. It is supported only at toy
scale. Instruction counts are analytic. This is a generic attack.

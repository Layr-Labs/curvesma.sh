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

**Claim.** time_log2 = 126.94. The algorithm recovers x from Q = [x]P on
secp256k1 with probability ≥ 0.39. On **every** run its total work is
W ≤ 2^126.94 units. One unit is one batched incomplete affine secp256k1 addition
(6M, A = 0) = C = 1,000 scalar AArch64 instructions. The bound is unconditional
(global work meter, section 3). Memory < 2^72 bytes, spread over 2^16 storage
nodes of < 2^56 bytes each. Preprocessing (everything before the first walk
step) ≤ 2^59.0 units, enforced by the meter and included in W. Nonuniform
advice: none. Success rests on the random-mapping heuristic H1 (section 4),
supported by toy runs.

## 0. Changes from 189ebc42, itemised

The algorithm, its class key, the doubling escape mechanism, the 2-cycle and
checkpoint rules, the work meter, the storage-node model and H1 are unchanged.
Eleven accounting and parameter items change:

1. **Inversion bound.** I < 2^17 instead of I ≤ 2^19, derived from the
   package's own M ≤ 227 (section 5). The batch-inversion share per element
   at ≥ 2^15 elements falls from 16 to 4 instructions. No other row of the
   per-step table changes for this item.
2. **Step price, first tightening.** The itemised per-step overhead fell from
   159 to 147 instructions, charged at a budget of 150, so the ordinary step
   cost 1.15 units (was 1.20 at a budget of 200) and the traversal step 1.18
   (was 1.22); the doubling price stayed 2.0. Section 5 also began deriving
   every row from an explicit sequence, which totalled 137 instructions.
   Item 8 replaces both numbers by the sequence's own count.
3. **Doubling rate, stated exactly.** With the table R of item 11 the 2-cycle
   rate is 1/(6R) = 2^-22.58496 per operation. Cycles of length 3..64 that need
   every partition index at least twice have rate 2^-38.00 (the sum over
   L = 3..64 and b = 1..⌊L/2⌋ of S2(L, b)·R^-(L−b)). The exact total is
   2^-22.58493 per operation, and the ledger uses the bound < 2^-22.58493.
4. **Meter headroom, stated from the ledger.** Expected random charges
   < 2^111.032 against headroom > 2^119.77 give an early-halt probability
   below 2^-8.74 = 0.0024 by Markov, so success is
   ≥ 0.395 − 0.0024 − 2·10^-7 = 0.3926 > 0.39. Phase 1 costs
   < 2^126.93997, and the claim was 2^126.95; items 9 to 11 carry it to
   2^126.94 (section 5.1).
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
   required log2 W from 127.00905 to 126.98951, which item 8 then lowers.
8. **Step price, from the explicit sequence.** Every row of the per-step
   table is now the exact count of a written-out scalar sequence, and the
   itemised total is 113 instructions, charged at a budget of 116. The
   ordinary step costs 1.116 units (was 1.15) and the traversal step 1.14
   (was 1.18); the doubling keeps its 150-instruction bookkeeping allowance
   and its price of 2.0, so the doubling surcharge rises from 0.85 to 0.884
   per doubling. The expected random total rises from 2^110.98 to 2^111.032
   and the Markov allowance with it, and the required log2 W falls from
   126.98951 to 126.94665.
9. **The pad goes, and the batch inversion is charged its own share.**
   Item 8 made every row the exact count of a written sequence, so the
   3-instruction pad above its 113-row table goes. Item 1's 4-instruction
   batch-inversion share comes back as a row of its own, because the
   organizer's 6M price does not supply it for free. The 6M is 3M for the
   addition formula plus 3M for batch amortisation (`docs/COST_MODEL.md`).
   That 3M per addition is spent on the batch's 3(L − 1) multiplications: over
   a batch of L additions it contributes 3M(L − 1)/L and leaves 3M/L of itself
   unspent, while the batch's one inversion contributes I/L. The residual
   I/L − 3M/L is therefore not covered by the reference price and is charged
   as its own row. Section 3 batches 2^15 to 2^16 ready walks per Montgomery
   inversion, and a batch below 2^15 pays its inversion as a separate
   600-unit action, so every walk the row charges sits in a batch of
   L ≥ 2^15. With I < 2^17 the share is below I/2^15 = 4 instructions, and the
   row charges 4. This is not a second charge for the reference operation: the
   formula and the amortisation stay inside the 6M, and only the residual the
   amortisation cannot cover is priced. The ordinary step is 1.113 units
   (109 + 4) and the traversal step 1.137 (133 + 4). Phase 1 falls by
   log2(1.116/1.113) = 0.00388 bits, and W by 0.00384, since the hard extras do
   not scale with the step.
10. **Success target at the requirement, and the headroom charged for it.**
   The H1 floor was 0.395, which left 0.0026 of probability above the 0.39
   requirement once the early-halt and other failure terms were subtracted. It
   is set to 0.394: success is 0.394 − 0.0024 − 2·10^-7 = 0.3916, a stated
   margin of 0.0016, and it stays at 0.3905 if the random charges run 45% above
   their bound (section 6). S_ops falls by
   √(ln(1/0.606)/ln(1/0.605)) = 0.998355, which is 0.00237 bits; in units of
   the old normalization this is k = 1.04827 instead of 1.05, and the cycle,
   delay and Markov terms are recomputed at the new S_ops. The ledger now
   charges the Markov allowance for the 0.0024 of early halt that the success
   statement subtracts, where the accepted revision charged the 0.0035 cap its
   arithmetic allowed while its success statement used 0.0024; that
   consistency costs 0.00307 bits of headroom (item 11, not this item, removes
   the need for it).
11. **R = 2^20 instead of 2^13.** The table holds 1,048,576 starts, < 2^33.1
   units to build and 2^28.59 bytes to hold, and N' = N·R/(R−1). The doubling
   rate becomes 1/(6R) + 2^-38.00 = 2^-22.584929 per operation, so the doubling
   surcharge falls from 2^111.03 to 2^104.02 and the Markov allowance it drives
   from 2^119.74 to 2^112.72 units. That is 0.00980 bits off the required
   log2 W. The index field moves from bits 0..12 to bits 0..19 of K; the
   distinguished (64..127), routing (128..143) and bucket (144..193) fields are
   untouched and still disjoint from it.

Gain: 127.08 − 126.94 = 0.14 bits. Items 1, 2 and 4 account for the move from
127.08 to 127.01 (item 2 gives 0.0614 bits on phase 1, item 4 the rest).
Items 6 and 7 account for the move from 127.01 to 126.99: item 7 drops the
required log2 W by 0.01954 (0.00026 from S_ops and 0.01927 from the random
total), item 6 by below 10^-5, and the 2-decimal ceiling by the remaining
0.00046. Item 8 accounts for the move from 126.99 to 126.95: phase 1 falls by
log2(1.15/1.116) = 0.04330 bits, the larger doubling surcharge adds back
0.00044 bits through the Markov allowance, the required log2 W falls by
0.04286 to 126.94665, and the 2-decimal ceiling takes the claim to 126.95.
Items 9 to 11 move the exact log2 W from 126.94665 to 126.93370, a total of
0.01295 bits: 0.003843 from item 9, 0.002375 from item 10 and 0.009803 from
item 11, less 0.003070 because the Markov allowance is now charged for the
0.0024 early-halt bound the success statement uses instead of the 0.0035 cap.
The 2-decimal ceiling then takes the claim from 126.95 to 126.94.

The ledger they leave, in the order of section 5: phase 1 is 2^126.93362, the
hard extras 2^107.155, the expected random extras 2^104.021 and the Markov
allowance they require 2^112.72; W = 2^126.93370 and the claim is 2^126.94.

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

Parameters: R = 2^20; distinguished probability θ = 2^-64 (limb 1 of K is
zero); M = 2^20 walks; L_max = 20·2^64; checkpoint period 64.

- **Table.** c_h and d_h are uniform in [0, r) (fresh coins, rejection sampling,
  every attempt metered). M_h = [c_h]P + [d_h]Q, redrawn if O.
  Tab[h][t] = [u_t^{-1}]M_h = (−1)^s·φ^{−j}(M_h), with coefficients
  u_t^{-1}·(c_h, d_h) mod r.
- **Step.** F(C) = class(c(C) + M_h) with h = K mod R. The walk keeps any member
  P and sets P ← P + Tab[h][t]. Since c + M_h = [u_t](P + Tab[h][t]), this is F.
  Each ordinary step is one batched affine addition (the reference operation)
  plus 113 itemised instructions (section 5).
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
  exceed W_max = 2^126.94, the attack halts with failure. So W ≤ 2^126.94 on
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
  index bits 0..19, the DP bits 64..127 and the bucket field 144..193) to one
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

N' = N·R/(R − 1) = N·2^20/(2^20 − 1), log2 N' = 253.41504. T_c is the number of
class arrivals of all trails, up to the first arrival at a class already
visited by another trail. It includes trail starts, steps, doublings,
traversal steps and repeated arrivals inside fruitless cycles; nothing is
excluded. This is exactly what the toy measures.

**H1 (heuristic, at the budget).** P(T_c ≤ S_ops) ≥ 1 − exp(−(S_ops/κ)²/(2N'))
with κ = 1.05. H1 also includes three further claims:

- the partition, DP, routing and bucket fields of |y| are pseudo-random
  functions of the class: bits 0..19, 64..127, 128..143 and 144..193, four
  disjoint fields;
- excluding the detected fruitless cycles of length ≤ 64, a trail makes a
  non-fruitless return with probability ≤ L_max²/(2N);
- the canonical twist index of a new class is uniform and independent of its
  partition index.

The 2-cycle rate is 1/(6R) = 2^-22.58496 per step. Other cycles of length ≤ 64
have rate ≤ Σ_L Σ_b S2(L, b)·R^-(L−b) = 2^-38.00, summed over L = 3..64
and b = 1..⌊L/2⌋, where S2(L, b) counts the partitions of the L steps into b
blocks of size ≥ 2 (every partition index used at least twice; for example
1, 1, 3 and 15 for (L, b) = (3, 1), (4, 1), (4, 2) and (6, 3)). Cycles of
length 65..200 have rate < 2^-206 per step. The expected doubling rate per
operation is therefore 1/(6R) + 2^-38.00 = 2^-22.58493, and the ledger uses
the bound < 2^-22.58493.

The values are:

    S0 = √(2N' ln(1/0.606)) = 2^126.70878
    S_ops = 1.05·S0 = 2^126.77917

H1 gives ≥ 0.394 (the floor of item 10); the pure model, with κ = 1, gives
0.4243. In the toy runs, class walks are −1.0% to +1.5% from the
random-mapping mean T_c, which is within the 5% allowance (section 6). The
floor sits 0.0016 above the 0.39 requirement, and 0.0005 above it in the
45%-excess case of section 6.

**Phase 1 covers the budget.** Phase 1 performs ≥ S_ops operations, and each is
one class arrival. So if T_c ≤ S_ops, the first inter-trail collision happens
in phase 1.

**Doubling and traversal surcharges.** These are random, and only their excess
over the step price is not fixed by the operation count:

- doublings: expected count < 2^-22.58493·(S_ops + 2^60), surcharge
  2.0 − 1.113 = 0.887 each, so expected < 2^104.021;
- traversals: expected < 0.024·64·2^-38.00·(S_ops + 2^60) < 2^89.40.

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
| meter halting early: expected random charges < 2^104.021 against headroom > 2^112.72, by Markov | < 0.0024 |
| preprocessing cap | < 2^-32 |

The events other than the early halt sum to < 2·10^-7.
**Success ≥ 0.394 − 0.0024 − 2·10^-7 = 0.3916 > 0.39.**

## 5. Cost

**Per ordinary step**, in instructions beyond the reference addition:

| work | instructions |
| --- | ---: |
| 1.1 F1: 8 MUL/UMULH, 7 ADD, 3 loads of G | 18 |
| 1.1 range test, S, wrap test, F2 | 8 |
| 1.1 decision: min, second smallest, argmin, gap, fallback branch | 12 |
| sign, p − y, select limbs of K | 13 |
| h, distinguished test, twist index, table address | 7 |
| loads of table point and state, store of state | 12 |
| fingerprint tests, history and checkpoint stores | 22 |
| zero-denominator test | 5 |
| length, cap test, meter additions, queue and loop | 12 |
| batch inversion: I/L − 3M/L < I/2^15 with I < 2^17 | 4 |
| **total (charged 113 = 0.113 units)** | **113** |

The step is charged at the itemised total, 113 instructions, and every row
below is the exact count of the scalar sequence written out for it, except the
batch-inversion row, which is an upper bound. No pad is added above the table,
because the table is a count and not an estimate.

The batch-inversion row needs its derivation, because the reference price
already carries batch work. The organizer's 6M is 3M for the addition formula
plus 3M for batch amortisation, and it has no standalone inversion term, so
the inversion a batch performs is charged as other work, outside the reference
operation. Inside the 6M, the amortisation's 3M per addition is consumed by
the batch's 3(L − 1) multiplications, which cost 3M(L − 1)/L per addition and
leave 3M/L of the 3M; the batch's one inversion costs I < 2^17 and contributes
I/L per addition. So the part the reference price does not cover is
I/L − 3M/L, and that is what the row charges. Section 3 splits ready walks into batches of
2^15..2^16 per inversion and makes a batch below 2^15 pay a separate 600-unit
action, so L ≥ 2^15 on every walk this table prices; then
(I − 3M)/L ≤ I/2^15 < 131,072/32,768 = 4. The row charges 4, and neither the
3M amortisation nor the 3M formula is charged twice.

Pair loads and stores are single instructions: `ldp` and `stp` move two 64-bit
words each. The walk record is 32-byte aligned, with x
at +0, y at +32, the two 32-byte history slots at +96 and +128, the
checkpoint at +160 and the pending flag at +192. The constants p0 (the low
limb of p), the all-ones word, 2^64 − 81, the index value 2, the meter price
and the table base stay in registers across the batch loop; the counts below
include every load and store of per-walk state:

- F1 (18): `ldp` G_1,G_2, `ldp` G_3,G_4 and `ldr` G_5 (3); four `mul` for the
  column-5 products x_i·G_{5−i} and four `umulh` for the column-4 products
  x_i·G_{4−i} (8); six `add` into two accumulators and one `add` that joins
  them (7). Total 18. The four x limbs are the registers loaded for the point
  traffic, so this row loads none of them.
- Range test, S, wrap test, F2 (8): the range test is
  (F1 − 40) mod 2^64 > 2^64 − 81: `sub` of 40 from F1, `cmp` with the
  2^64 − 81 register, `b.hi` (3). S = F0 + F1 by `adds`, which leaves the
  carry c (1); F0 is x_3, already in a register (0). The wrap value is S when
  c = 1 and −S when c = 0: `csneg T, S, S, cc` writes it without touching the
  flags (1), and `cmp T, 40` with `b.ls` tests it (2). F2 itself is
  `neg F2, S` (1), placed after those branches because it overwrites the
  carry. Total 8. When c = 0 and S = 0, F0 = F1 = F2 = 0, so the gap test of
  the next row fails and the fallback runs; the wrap test needs no separate
  case for it.
- Decision (12): `cmp` F0,F1, `csel` m1 = min, `csel` M1 = max and `cset` of
  the index of m1 in {0,1} (4); `cmp` m1,F2 (1); `csel` mn = min(m1,F2),
  `csel` t2 = min(M1,F2), `csel` mid = m1 if F2 < m1 else t2 and `csel`
  j = 2 if F2 < m1 else the index (4); `sub` gap = mid − mn, `cmp` gap with
  40, `b.ls` (3). Total 12. A gap above 40 to the second smallest is a gap
  above 40 to both others, so the smallest of the three is the exact argmin.
- Sign, p − y, K (13): p − y by `subs` and 3 `sbcs` against p0 and the
  all-ones word (4); p − y < y by `subs` and 3 `sbcs` into one scratch
  register, which carries the borrow chain exactly because each step needs
  only the flags and its two operands (4); `cset` s from that carry (1);
  4 `csel` for the limbs of K (4). Total 13.
- h, distinguished test, t, address (7): `and` h = K_0 mod 2^20 (1); `cbz` on
  limb 1 of K (1); t = 3s + j by `add w, s, s, lsl #1` and `add t, w, j` (2);
  6h + t by `add w, h, h, lsl #1`, which is 3h, and `add w, t, w, lsl #1`,
  which is t + 6h (2); `add xT, table base, w, lsl #6` (1). Total 7.
- Point and state traffic (12): two `ldp` for x and two for y (4), four `ldp`
  for the 64-byte table point (4), four `stp` for the new point (4). Total 12.
  One base register with immediate displacements reaches the whole walk
  record, and the address computed above reaches the table entry; the four x
  registers serve F1, the sign row and the zero-denominator test, and the
  four y registers serve the sign row and the reference addition, so no limb
  is loaded twice. The eight words of the updated point are written by four
  `stp` and by nothing else, so the pairing adds no store.
- Fingerprints (22): pending flag `ldr` and `cbnz` (2); ring slot address by
  `and` of length with 1 and `add` of 32 times that to the record base (2);
  the slot that holds K_{-2}: `ldp` and `ldr` for its three words (2), `cmp`,
  2 `ccmp` and `b.eq` (4); overwrite it with fp(K) by `stp` and `str` (2);
  checkpoint compare: `ldp` and `ldr` (2), `cmp`, 2 `ccmp` and `b.eq` (4);
  checkpoint write when length mod 64 = 0: `tst` and `b.ne` (2), `stp` and
  `str` (2). Total 22. The two 3-word groups cannot share an address
  register, because one is the alternating history slot at a data-dependent
  offset and the other is the fixed checkpoint field, and a 4-word store is
  no cheaper than a 3-word one, since both are two instructions. The saving
  is the pairing of the loads and stores, four instructions in all.
- Zero denominator (5): `cmp` x(P) with x(Tab) and three `ccmp`, `b.eq` (5),
  on registers loaded above.
- Length, cap, meter, queue, loop (12): the 2-limb length (L_max > 2^64) by
  `ldp` (1), `adds` and `adc` (2), `stp` (1); cap test by `cmp` of the high
  limb with 20 and `b.hs` (2); meter: `subs` on the processor's block
  remainder, `add` on its operation count, `b.lo` to the block refill, where
  the global counters and the phase-1 limit are checked (3); the queue
  pointer by `ldr` with post-index (1); loop counter `subs` and `b.ne` (2).
  Total 12.
- Batch inversion (4): an upper bound, not a counted sequence. The reference
  price's 3M amortisation pays the batch's 3(L − 1) prefix-product and
  back-substitution multiplications and their traffic; the batch's single
  inversion, I < 2^17 instructions, is outside it, and divided over the
  L ≥ 2^15 walk steps that share the batch it is at most I/2^15 < 4 per step.
  Nothing the reference price already pays is charged here.

The sequences total 18 + 8 + 12 + 13 + 7 + 12 + 22 + 5 + 12 + 4 = 113
instructions, and the step is charged 113. Rare paths are priced as their own
actions: the 2-cycle order test, the flag write and the escape are in the
doubling price, the fallback has its own price, and so do reports and
traversal steps.

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
  510·251 + 2,048 = 130,058 < 2^17. The batched walk charges its own share of
  this bound, the 4-instruction row of the step table, and the unbatched
  600-unit actions price the whole inversion; their price of 600,000
  instructions is more than four times the bound.

**Prices** (units):

| action | price | basis |
|---|---|---|
| ordinary step | 1.113 | reference addition plus 113 itemised instructions |
| traversal step | 1.137 | ≤ 24 more instructions to copy the remembered member: ≤ 137 |
| doubling | 2.0 | batched affine doubling: 7M + 6A + 150 ≤ 1,931 instructions, unchanged |
| §1.1 fallback | 0.5 | one M, two A, compares: < 400 instructions |
| random 64-bit word | 0.01 | |
| sampling attempt | 0.05 | |
| trail start [a]P + [b]Q | 2^13 | two Jacobian scalar multiplications at ≤ 512 operations × (16M + 16A) each, an addition and an inversion: < 4.8·10^6 instructions |
| store probe | 0.05 | 39 instructions itemised in section 3, charged 50 |
| report | 40.4 | 40 for route, message, entry write and node tests (< 2^15 instructions), plus 8 probes at 0.05 |
| overflow insert scan | 2^41.45 | worst case 6B + 20 instructions |
| overflow lookup | 2^42.86 | worst case 16B instructions over the B-entry list |
| store initialisation | 2^16·B buckets × 16 instructions = 2^58.864, once | |
| table | 2^20 starts plus 5 φ-twists each: < 2^33.1, plus sampling | |
| phase-2 or re-walk operation | 600 | unbatched: an inversion (< 2^17) plus 7M, coefficient updates mod r and control, far below 600,000 instructions |
| small-batch inversion | 600 | one inversion, < 2^17 instructions = 131 units; 600 is a loose upper bound |
| solving | 600 | |
| verification | 2^13 | |

**Budget.**

- **Phase-1 operations:** 1.113·(S_ops + 2^60) < 2^126.93362.
- **Hard extras:**

  | item | units |
  |---|---|
  | meter remainders | < 2^61 |
  | phase 2: 2^20·(L_max + 65)·600 | 2^97.551 |
  | re-walks: 8·2·((L_max + 65)·600 + 2^13) | 2^81.551 |
  | reports: D_cap·40.4 | 2^69.166 |
  | store initialisation: 2^16·B × 16 instructions | 2^58.864 |
  | table and verification | < 2^33.1 |
  | overflow insert and lookup: D_cap·(22B + 20) instructions | 2^107.154 |
  | **total** | **< 2^107.155** |

- **Random extras (expected):**

  | item | units |
  |---|---|
  | doubling surcharges: 0.887·2^-22.58493·(S_ops + 2^60) | < 2^104.021 |
  | traversal surcharges: 0.024·64·2^-38.00·(S_ops + 2^60) | < 2^89.40 |
  | starts: < 2^62.8·2^13 | 2^76 |
  | words, attempts, small batches, fallbacks | < 2^70 |
  | **total** | **< 2^104.021** |

- **Headroom:** 2^126.94 − 2^126.93362 − 2^107.155 > 2^112.72.
- **Early meter halt:** probability < 2^104.021/2^112.72 = 2^-8.70 = 0.0024
  (Markov).

### 5.1 Exact total and the two-decimal claim

All ledger quantities are carried as integers in milli-units (10^-3 units; one
instruction is one milli-unit), so every product is exact. With
S_int = ⌈S_ops⌉ = 145,993,368,281,645,653,992,586,223,245,627,523,266,
D_cap = 16,396,253,461,269,537,204 and B = 500,373,945,961,595:

    phase 1            1113·(S_int + 2^60)                        = 162490618897471612895031668107010854079346
    hard extras total (section 5, < 2^107.155)                    =    180725794254648202537795155729156760
    Markov allowance   ceil(1250·random/3)                        =   8576725238500405016544809380556303334
    W (milli-units)                                               = 162499376348504367948250750711547139539440

W in units is W/1000 = 2^126.933701787…, whose 2-decimal ceiling is 126.94.
The exact integer certificate is W^100 ≤ 1000^100·2^12694, and W^100 >
1000^100·2^12693, so 126.93 does not suffice. The two sides of the first
inequality are the integers

    W^100 =
    121666175146079039450723988636086846375869765823614669730895511126046322068763896076672349838357375500336897941432029088550736117648337891341083423493527870250768930372293039780905443568057362981726398361970982785378799031488536326218025950467110920527405552807630349309244900116286660437183930977636318823682689061660846766943746211032893983318360795762026367084572025095891920645818624143704082977296538971347879026031388841195323656535754177633745230759556180115760360432048650363888218925946943148884740152131923355958412348977371589678512719878331715520872776778426449485499556530179261080344356534167119514238136175958692502412248571305974875852613789918112247427306433071203685867212855712625199622691149556443886432544106889222345799471426420279651783073912420470562830181118563357069981600639829250962859946656178519564297283785735654230900170027923340375764939032120087525910772828826709927724374806583293197848909677808354805976290595834333977086003601467132758853465055059503348921390742731788533576947932442905441106031326349602513073704913449868611468417669951253913878170008506175440519442049907724683811439134821226200650697085786664183319774017436622215818692482099509406428266560134373432480943653489159513952403933877701952656720160442103273185510025285342754775797443420168618165104757906749897977713381922364368961372945296664935714776849861286424290925715997608605664737646633406040033293579207145511438501178582885335371565679234369806537377137434556167177238816591162350020811447446659485230847699256937823625043597445663365795188418317873365705166228915465170235712055380960926959738796613490654606604132127015768656670426279853597352821662593647981426680176663125650769068866677995769662184015560913304273982128517261468715310753933525329255472315756438079540375045630422888587839074604395308850619602988540663359635910566824907396627384976726517027716903758232506605278045570221305438783445913684945500246779402744771221853476551289470102540650738585284322014199197650625310391731653751544051853535906002041355298387183403823026176359946851404591424013170693875041623113025012424982182352759035628211895790245466943435283268388051179574270119893459296621421786733342501291171000518029786328331947354955166844726827709203815267120813327378792087373435583068553100511632825971996090969908754927566250417424900534287979284169392290939142324162670853995645886023194831005437105779388638864682456781382136848294922618162651597377623848748877847688952178072651332790647348100984595630645873888250882397861701157492483550440569611565808370731302816029165972557662408343087786710630191226181047800582697206575146042751750728823956326797509345235635199126392483246994561026244608338562390110065568784391709494067721233414488061823833636853108162952939689448787124403797993370033746178870987584395699018364265880654290171650888757101910443979614875064733990950809558098390856084593961367066203608224419206009242026641134560833900616586752007850470544342697294363064187802080417950798782964814374544304834941629014322944341527950989877923573482539856853756002466663845505057069063162093496250895300752009186651396821775274874638656109557751914333596970780366885611196048726342860229095597958433921029604800776721417721657315423687659328804368787576386693954213735309428311578718621207644025215689596244314716044792355868802677825539620269047297868103699570946309546305143447366312066458595328750316386369141396826602854069019388055810027315407325978211597631395709656931542949637591466173478058975473568380410458271711116211525086305108092636238135327020968964324175166181589946352723910613068947142406223625511296369241897320556580389891958900739937390159260202663526000825248544010292834084213614258954223089744587816102178092721725551839313330598273390435543890148492278298742624728988466167877954366739740431868941585265522920768389301842511089670165706457666680400927385049471102139419823633011655067146511983050744650405807513797559427201868887273343234328833123916209162479319172326182340873955573760000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

    1000^100·2^12694 =
    188262992899605728344851981937231547106065625799283710965276533978192541114146283093903657621669549396496859898857180641056620643990589052562244836770564178438153342946364488201593857808895132628417091480500704121493730441231232412249127617136098835150697788085899177515258386356643079515878297532433313666278229052974760064700074966460235874781633878822794994895994557651472772881319875480475719291252612272652873583087956715471941019579446631282031256800189809716714892418598121390457231889041476765880070649470117812842777621552068956297018161399997014321442033034618741983061371078581082411937222920036103737084201634513615743090522860845798634124277373327815033841652258693953292426106856257375709739505625092883754150653338202601478063855394368545123636372425713587899817472479702940223229583847818159979586832464822772736704711078118220750093065191793671894416648197584330553314149881328140663552770642388412545976016566083034148395024596850658273809817227256060979527805512414768975712149202147230488376412857157679782407327050165739833619218183367835054802433778990408126874691696303509610621776996625805852107284278492315379991184803281268784307541309820270227948314092589414307345974843740606319613423348081433374749671874417416881942117099719244351471232052166348719998527521544298246393310975628316407088810589346272331200208678226155435314966462150306976222442366403948726441451052640527638656221304700085825835568711749335726813258074607990137965774029755032689708205783172671566763822124315615171262388531332232874671066725855110662608697813385896944120862301291570426505592602420425690647386334788165451077062390747628810589930388373889599225778443568614841116913402558286655430249282525243087176221900003181524938782827214859889550997211493539192077749208945984829941854790468983611558464194795616573162338012255410503830052327137641921192317920407847806565887243017822033820480835037561971823071493791736464289195910842538910735635393587187051005785006685342077417078979506313126294403488396883926025377043977256198527835245096540316915841811349149120829269802761340948498870583339494212665212812837797488911559675734367764656747932556583678899913906806102880844500571281254365085270389336432344457314627073915513161676094235244942634815035644432261816767943565724350307740660002356415416099295771582479938341944693923479791468358700702914647325882624383741374027798296591445982681959413438926775736435930184784266208091152968782370237039792883565291076539156548802429619534891610780247632340422211213754059726744510063489356723666450039132993965095157003669782131266691143014810240935877133863042613510025384593437880538419978300464455867209804968435581821346058624343042188517763159069149034575035333912581606652136738749136181506469562275722619243014047283407920200188135188243086606499452636297986556024272418069563776712306735479055912453586097154282327358702597131963456223245374940653664186964821582449262513110000084348036422710039175822277694096231173954746657659101462107580230856440137126705729102779552242003207293489813190223363618624731848032668693135223008175221336153628438943523808066586750073554040658849803999287524849771355159522198519717790075129073013147897350949404668938953039997540368915925306679269278845778976884172237016547804218892445613005311726696097603462178286322959562676479379836287813476807506070232477456855402116696299183678754463030676878267762224485240344158825687619083173225088679157928384428560787017073573682940625833730693756875438639868900770692320206447440936990872731384164485289607190955831707136434273058256035816050864524012926407609158128751931850357867890710382205694985024435217179812305636004858979382531918869824904326107202846211470976627889393359867734768060304573417867064385036049332129456788229878790134898995824381510375232145396459146051584000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

Both carry 4,122 decimal digits; the comparison is exact integer arithmetic.
The ratio W/(2^126.94·1000) is 2^-0.006298, so the claim has 0.0063 bits of
room above the ledger's own exact total. The score claimed is 12694/100.

**W ≤ 2^126.94 on every run.**

**Preprocessing:** table < 2^33.1, node set-up 2^26, store initialisation
2^58.864, and the first 2^20 starts 2^33. The expected total is < 2^58.87,
and the cap is 2^59.0. The random sampling part is expected < 2^25, far below
the gap, so the cap halts a run with probability < 2^-30 by Markov. All of it
is inside W.

**Memory:** the store holds 2·D_cap slots of 104 bytes at load 1/2, so
2^64.83 slots × 104 bytes = 2^71.53 bytes, and the probe index arrays add
8 bytes per slot, 2^67.83 bytes: 2^71.64 bytes in total. Adding the table
(2^20 starts × 6 twisted entries × 64 bytes = 2^28.59 bytes, or 2^48.59 bytes
if every processor keeps its own copy, which is still far below the ceiling)
and the 2^20 × 512 bytes of walk states gives < 2^72 bytes, with each of the
2^16 nodes at 112 bytes per bucket, 2^55.64 bytes < 2^56.

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
S_ops = 1.05·√(2N' ln(1/0.605)) for each toy, the floor the toys were run at;
T_c counts all arrivals):

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
so the structural rate 2^-38.00 governs. The R = 2048 rows show the same
effect: 8.2·10^-5 in the T_c rows against their own 1/(6R) = 8.14·10^-5, and
8.72·10^-5 total at q25 with checkpoint escapes near 5·10^-6. The early-halt
margin covers it: with all random charges 45% above their bound 2^104.021,
Markov gives 1.45·0.0024 = 0.0035, and success is ≥ 0.394 − 0.0035 −
2·10^-7 = 0.3905 > 0.39. Charging the measured excess instead
(1.2·0.0024 = 0.0029) gives success ≥ 0.3911.

**Why the toys stop at R = 8192 and the claim is at R = 2^20.** The two toy
curves have N = 5,593,903 and N = 89,473,200 classes, that is 2^22.4 and
2^26.4, and the model needs R ≪ N for the partition to look random and
N' = N·R/(R−1) to hold; at R = 8192 the larger toy is already at
R/N = 2^-13.4 and the smaller at 2^-9.4. The claim's R = 2^20 is
R/N = 2^-233.4, so the regime is the same and only the *rate formula* is
extrapolated, in the direction that makes each cycle rarer. Three table sizes
(256, 2048, 8192) measure the 2-cycle rate against 1/(6R) to within 3%, and
the formula is exact for a 2-cycle: it needs one repeat in R partition buckets
and 1 of 6 class positions.

**Provenance of the rows.** `experiment.py` is byte-identical to the file in
189ebc42, and the rows above are its results for the configurations in
README.md. Its printed `model_P_S0` and `H1_floor_P_S1` both read 0.395: that
is the floor the toy rows were run at, and the claim's floor is 0.394
(item 10), 0.0025 bits of S_ops below it, so a toy row of this package is a
slightly longer run than the claim's budget implies. Python's seeded generator makes each row deterministic, so each row
can be rerun from its configuration. For this package the four R = 8192 rows
were run with seeds 1502, 1504, 1507 and 1508 under the same protocols as the
R = 2048 rows they extend, and the R = 2048 rows of seeds 1501, 1502 and 1507
were also rerun and reproduced exactly (1.2689 and 0.418; 1.2696 and 0.4185;
pipeline success 0.4229 with 0 node overflows and max load/cap 0.72).

**Initial input.** The hosted configuration (q25, R = 8192, dp_bits = 5,
seed 1, pipeline, 14 s) is the largest toy-valid table size and the one the
schema's enum reaches; it is time-limited, so the number of completed trials
depends on the machine. One local run of it completed 252 trials in 14.01 s:
part A all true, 0/30,000 random and 0/10,000 adversarial index-rule
mismatches (3,472 of the adversarial inputs took the metered fallback), a
pipeline success rate of 0.4603 with a binomial standard error of 0.0307,
0 node overflows and a maximum load of 0.52 of capacity. This short run is
consistent with the 12,000-trial pipeline row (0.4275) and with the measured
2-cycle rate (2.12·10^-5 against 1/(6R) = 2.03·10^-5); the H1 floor is taken
from the large offline rows, not from the 14-second run.

**Limitations.** H1 is a heuristic: a random mapping with the Brent–Pollard
correction; van Oorschot–Wiener; automorphism classes as in Wiener–Zuccherato,
GLV and Duursma–Gaudry–Morain; and fruitless-cycle escapes as in
Bos–Kleinjung–Lenstra and Bernstein–Lange–Schwabe. It is supported only at toy
scale, and the claim's table size R = 2^20 is an extrapolation of the measured
1/(6R) rate, not a toy row of its own. Instruction counts are analytic. The
general-rho success model itself is not proved at production scale. This is a
generic attack and not a practical break (the memory is 2^71.64 bytes).

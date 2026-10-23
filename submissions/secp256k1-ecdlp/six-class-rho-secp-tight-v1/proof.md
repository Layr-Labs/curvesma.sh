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

**Claim.** time_log2 = 126.92. The algorithm recovers x from Q = [x]P on
secp256k1 with probability ≥ 0.39. On **every** run its total work is
W ≤ 2^126.92 units. One unit is one batched incomplete affine secp256k1 addition
(6M, A = 0) = C = 1,000 scalar AArch64 instructions. The bound is unconditional
(global work meter, section 3). Memory < 2^72 bytes, spread over 2^16 storage
nodes of < 2^56 bytes each. Preprocessing (everything before the first walk
step) ≤ 2^59.0 units, enforced by the meter and included in W. Nonuniform
advice: none. Success rests on the random-mapping heuristic H1 (section 4),
supported by toy runs.

## 0. Changes from 189ebc42, itemised

The algorithm, its class key, the doubling escape mechanism, the 2-cycle and
checkpoint rules, the work meter, the storage-node model and H1 are unchanged.
Twelve accounting and parameter items change:

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
   < 2^126.93997, and the claim was 2^126.95; items 9 to 11 carried it to
   2^126.94 and item 12 to 2^126.92, where item 13 leaves it (section 5.1).
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
   both are charged: the insert at its own worst case 6B + 20 and the lookup
   at its price 2^42.87 units, so the row is
   D_cap·(6B + 20 + ⌈2^42.87·1000⌉) instructions = 2^107.158 units in the
   hard extras. A report now costs 40.4 units: 40 for routing, message,
   insert and node-side tests, plus its 8 probes at 0.05.
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
12. **Per-step rows tightened by 15 instructions.** Section 5 re-derives the
   bookkeeping row by row and the itemised total falls from 113 to 98, so the
   ordinary step is 1.098 units and the traversal step 1.122; the doubling
   keeps its price of 2.0 and its 150-instruction bookkeeping allowance, so its
   surcharge rises from 0.887 to 0.902. Five changes, each written out in
   section 5: `madd` folds three of the four column-5 low products into one
   accumulator, so one `mul` and three `madd` replace four `mul` and three `add`
   (−3); two `madd` give the twist index and the table address (−2); the history
   slot is addressed by the post-index of its own `ldp` and by negative
   displacements on its stores, with the pending-escape flag in a register (−3);
   the zero-denominator compare folds into the batch's Montgomery prefix pass,
   which forms the very product it tests (−5); and the 2-limb trail length stays
   in registers across the batch loop (−2). Phase 1 falls by
   log2(1.113/1.098) = 0.01958 bits; the larger surcharge adds back 0.000003
   through the Markov allowance; and the required log2 W falls from 126.93370
   to 126.91413, so the 2-decimal claim moves from 126.94 to 126.92. The walk,
   the class key, the escape and checkpoint rules, the meter, the success
   argument, the store, probe and overflow accounting, R and the
   batch-inversion derivation are the ones of items 1..11, unchanged.
13. **Three scalar sequences corrected, one price raised.** (a) The decision
   row read the flags of `cmp m1, F2` while selecting `min(M1, F2)`, which
   those flags do not decide when `m1 < F2 < M1`; it now takes a second `cmp`
   and a third `csel` and costs 13, not 12. (b) The wrap test is written as
   `CSNEG T, S, S, cc` with its exact meaning — `T = −S` when the carry is
   clear, `S` when it is set, which is the distance `abs(S − 2^64)` in both
   cases — in the same one instruction. (c) The overflow lookup was priced
   `2^42.86` for work of log2 42.86422, so that price did not upper-bound it;
   it is charged `2^42.87` (8,038,146,442,837,048 milli-units) and the hard
   extras carry 2^107.158 in place of 2^107.154. Every other row was audited
   against `fast_j` (section 1.1) and none changed. The itemised total is 99,
   so the ordinary step is 1.099 units (was 1.098) and the traversal step
   1.123; the doubling keeps 2.0 and its surcharge falls from 0.902 to 0.901.
   The claim stays 126.92, on 0.0046 bits of room instead of 0.0059.

Gain: 127.08 − 126.92 = 0.16 bits. Items 1, 2 and 4 account for the move from
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
The 2-decimal ceiling then takes the claim from 126.95 to 126.94. Item 12
moves the exact log2 W from 126.93370 to 126.91413, 0.01957 bits, and the
claim from 126.94 to 126.92. Item 13 moves it back to 126.91544, 0.00131 bits,
inside the same 2-decimal claim. Every price and total stated by items 1..12 is
the number that item had at that revision; where one of them differs from the
current one, item 13 supersedes it, and sections 4, 5 and 5.1 carry the current
numbers.

The ledger they leave, in the order of section 5: phase 1 is 2^126.91537, the
hard extras 2^107.16, the expected random extras 2^104.044 and the Markov
allowance they require 2^112.747; W = 2^126.91544 and the claim is 2^126.92.

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

**`fast_j` is the reference semantics.** `fast_j` in `experiment.py` defines
this rule, and every scalar sequence in section 5 reproduces it operand for
operand: `F1` is its eight limb products summed mod 2^64 in the same columns;
the range test is its `F1 < 40 or F1 > M64 − 40`; `S`, the wrap test and `F2`
are its `S = F0 + F1`, its `abs(S − 2^64) <= 40` on the true integer sum and
its `F2 = 2^64 − S` or `2^65 − S`; the decision row is its argmin over
`[F0, F1, F2]` and its "every other value is more than 40 above the minimum"
test; and the other rows are its `q − y < y`, `K = min(y, q − y)`,
`t = 3s + j`, `h = K & (R − 1)` and `fp(K)`. A sequence's condition code is the
one that reproduces the Python branch — `lo` for `<`, `hi` for `>`, `hs` for
`>=`, `ls` for `<=` — reading flags that its own immediately preceding
`cmp`/`adds`/`subs` sets, and a carry is the wrap of the true sum or the borrow
of the true difference.

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
  plus 99 itemised instructions (section 5).
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
  exceed W_max = 2^126.92, the attack halts with failure. So W ≤ 2^126.92 on
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
    instructions = 8,005,983,135,385.52 units, whose log2 is 42.86422. The
    price charged is 2^42.87 units = ⌈1000·2^42.87⌉ = 8,038,146,442,837,048
    milli-units, which upper-bounds those 16B instructions.
  - **Charged worst case.** Each report performs at most one lookup and one
    insert and probes at most 8 table buckets, so the probe and overflow work
    over all reports is at most D_cap·(6B + 20 + 8,038,146,442,837,048)
    instructions = 1.81021·10^35 milli-units = 2^107.158 units: the insert at
    its own worst case 6B + 20 and the lookup at its raised price. That is the
    ledger's overflow row; the 8 probes are the 0.4 units of the report price.
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
  2.0 − 1.099 = 0.901 each, so expected < 2^104.044;
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
| meter halting early: expected random charges < 2^104.044 against headroom > 2^112.74, by Markov | < 0.0024 |
| preprocessing cap | < 2^-32 |

The events other than the early halt sum to < 2·10^-7.
**Success ≥ 0.394 − 0.0024 − 2·10^-7 = 0.3916 > 0.39.**

## 5. Cost

**Per ordinary step**, in instructions beyond the reference addition. The right
column is the count of the sequence written out for that row below; the third
column is what the previous revision charged, and section 0 item 12 names the
five changes:

| work | instructions | previous |
| --- | ---: | ---: |
| 1.1 F1: 4 MUL/MADD, 4 UMULH, 4 ADD, 3 loads of G | 15 | 18 |
| 1.1 range test, S, wrap test, F2 | 8 | 8 |
| 1.1 decision: min, max, argmin, median (two CMP), gap, fallback | 13 | 12 |
| sign, p − y, select limbs of K | 13 | 13 |
| h, distinguished test, twist index, table address (two MADD) | 5 | 7 |
| loads of table point and state, store of state | 12 | 12 |
| fingerprint tests, history and checkpoint stores | 19 | 22 |
| zero-denominator test, folded into the batch inversion | 0 | 5 |
| length, cap test, meter additions, queue and loop | 10 | 12 |
| batch inversion: I/L − 3M/L < I/2^15 with I < 2^17 | 4 | 4 |
| **total (charged 99 = 0.099 units)** | **99** | **113** |

The step is charged at the itemised total, 99 instructions, and every row
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
(I − 3M)/L ≤ I/2^15 < 131,072/32,768 = 4. The row also carries the batch's own
zero-product test, four instructions on the product its prefix pass has just
formed, which is at most 4/L ≤ 4/2^15 per step: with I < 2^17 and the
amortisation's 3M/2^15 ≥ 537/32,768 = 0.0164 subtracted, the row's true value
is below 4 − 0.0164 + 0.00013 < 3.984. The row charges 4, and neither the 3M
amortisation nor the 3M formula is charged twice.

Pair loads and stores are single instructions: `ldp` and `stp` move two 64-bit
words each. The walk record is 32-byte aligned, with x at +0, y at +32, the two
32-byte history slots at +64 and +96, and the checkpoint fingerprint at +128
(three words, +152 spare). The constants p0 (the low limb of p), the all-ones
word, 2^64 − 81, the index values 2, 3 and 6, the meter price and the table
base stay in registers across the batch loop, as do the record address, the
address of the history slot last written, the pending-escape flag and the two
limbs of the trail length; the counts below include every load and store of
per-walk state that is not one of those registers:

- F1 (15, was 18): `ldp` G_1,G_2, `ldp` G_3,G_4 and `ldr` G_5 (3); one `mul`
  for the first column-5 product x_0·G_5 and three `madd` that accumulate
  x_1·G_4, x_2·G_3 and x_3·G_2 into the same register (4); four `umulh` for
  the column-4 products x_i·G_{4−i}, each followed by one `add` into that
  register (8). Total 15. The four x limbs are the registers loaded for the
  point traffic, so this row loads none of them. The three instructions that go
  are adds: a `mul` cannot be folded into an add, but `madd` folds three of the
  four low products into the accumulator, so only the four high halves still
  need their own `add`. What is summed is unchanged, column by column, mod
  2^64, and so is the ±14 error bound of section 1.1.
- Range test, S, wrap test, F2 (8): the range test is `F1 < 40 or
  F1 > 2^64 − 41`, which `fast_j` writes as `F1 < TOL or F1 > M64 − TOL`:
  `sub t0, F1, #40` into a scratch register, because F1 is still needed below,
  `cmp t0` with the 2^64 − 81 register and `b.hi` (3). Then `adds S, F0, F1`
  (1), with F0 = x_3 already in a register (0), leaves C = 1 exactly when the
  true sum F0 + F1 reached 2^64, i.e. when `fast_j`'s `S < (1 << 64)` is false.
  The value the wrap test needs is d = |S − 2^64|, and one `csneg` produces it:
  `csneg T, S, S, cc` is CSNEG T, S, S, CC, returning −S when CC holds (C = 0)
  and S when C = 1; it reads C and writes no flags (1). Both branches equal d:
  with C = 0 the true sum is S < 2^64 and d = 2^64 − S = −S mod 2^64, and with
  C = 1 the true sum is S + 2^64 and d = S. Then `cmp T, #40` and `b.ls` (2) fall back exactly when d ≤ 40, which
  is `fast_j`'s `abs(S − 2^64) <= TOL`. F2 is `neg F2, S` (1), placed after
  those branches because it would overwrite the carry; it is 2^64 − S when the
  sum did not wrap and 2^65 − S when it did, which is `fast_j`'s two-branch
  `F2`, the two differing by exactly 2^64. Total 8. When C = 0 and S = 0 the
  range test has already fired on F1 = 0, so the wrap test needs no case of
  its own.
- Decision (13): `cmp F0, F1` (1); `csel m1, F0, F1, lo`, `csel M1, F1, F0,
  lo` and `cset jm, hs`, which is 1 when F0 ≥ F1 and 0 otherwise, that is the
  index of m1 in {0,1} (3); `cmp m1, F2` (1); `csel mn, m1, F2, lo` and
  `csel u1, F2, m1, lo`, the min and max of m1 and F2 (2); `csel j, two, jm,
  hi`, which is 2 when m1 > F2 and jm otherwise (1); `cmp M1, u1` and
  `csel mid, u1, M1, lo`, the min of M1 and u1 (2); `sub gap, mid, mn`,
  `cmp gap, #40` and `b.ls` (3). Total 13. The second group reads the flags of
  `cmp m1, F2`, and `mid` reads its own `cmp M1, u1`; no flag is read after an
  overwrite, and `j` is written before that second comparison. `mid` is the
  median: min(M1, max(m1, F2)) is m1 when F2 ≤ m1, F2 when m1 < F2 < M1, and
  M1 when F2 ≥ M1. So `gap` is the distance from the smallest value to the
  second smallest, and `b.ls` falls back exactly when it is at most 40, which
  is `fast_j`'s test that every other value lies more than TOL above the
  minimum. A gap above 40 leaves the three distinct, so `mn` and `j` are the
  exact argmin and argmin index that `fast_j`'s `min(range(3), key=…)` returns;
  ties fall back, so the two tie-breaking orders never have to agree. The
  extra instruction is the second `cmp` that makes `mid` decidable.
- Sign, p − y, K (13): p − y by `subs` and 3 `sbcs` against p0 and the
  all-ones word (4); p − y < y by `subs` and 3 `sbcs` into one scratch
  register, which carries the borrow chain exactly because each step needs
  only the flags and its two operands (4); `cset` s from that carry (1);
  4 `csel` for the limbs of K (4). Total 13.
- h, distinguished test, t, address (5, was 7): `and h, K_0, #0xFFFFF` (1);
  `cbz` on limb 1 of K (1); one `madd` of s by the register that holds 3, plus
  j, which is t = 3s + j (1); one `madd` of h by the register that holds 6,
  plus t, which is 6h + t (1); `add xT, table base, w, lsl #6` (1). Total 5.
  The two instructions that go are adds: the previous count built 3s, then t,
  then 3h, then 6h + t with four shifted `add`, and a `madd` multiplies and
  adds in one instruction. The multipliers 3 and 6 are two of the constants
  held across the batch loop, the index t and the address are the same, and the
  distinguished test is untouched.
- Point and state traffic (12): two `ldp` for x and two for y (4), four `ldp`
  for the 64-byte table point (4), four `stp` for the new point (4). Total 12.
  One base register with immediate displacements reaches the whole walk
  record, and the address computed above reaches the table entry; the four x
  registers serve F1, the four y registers serve the sign row and the
  reference addition, so no limb is loaded twice. The eight words of the
  updated point are written by four `stp` and by nothing else, so the pairing
  adds no store. A walk's point stays in the record between its steps because
  a batch interleaves the walks that share one inversion, which is what the
  batching in section 3 is for, so this traffic is the batch's own and cannot
  be removed by holding one walk's point in registers.
- Fingerprints (19, was 22): pending-escape flag `cbnz` on the register that
  holds it (1); the slot that holds K_{-2}: `ldp f_0, f_1, [hp], #32` and
  `ldr f_2, [hp, #−16]` for its three words (2), `cmp`, 2 `ccmp` and `b.eq`
  (4); overwrite it with fp(K) by `stp K_0, K_2, [hp, #−32]` and
  `str K_3, [hp, #−16]` (2); checkpoint compare: `ldp` and `ldr` (2), `cmp`,
  2 `ccmp` and `b.eq` (4); checkpoint write when length mod 64 = 0: `tst`,
  `b.ne`, `stp` and `str` (4). Total 19. Three instructions go. The flag is one
  of the registers held across the batch loop, so its `ldr` goes; and the two
  history slots are 32 bytes apart in a 32-byte aligned record, so the
  post-index of the `ldp` that reads the slot leaves `hp` on the other slot and
  the two stores reach the words just read by negative displacements, which
  removes the `and` of the length with 1 and the `add` that scaled it. The flag
  is set by the 2-cycle branch and cleared by an escape, and a resumed walk
  re-derives `hp` from the length parity with two instructions; all three are
  inside the doubling price's 150-instruction allowance.
- Zero denominator (0, was 5): the per-step `cmp` of x(P) with x(Tab), its
  three `ccmp` and its `b.eq` are gone. The batch already multiplies its L
  denominators together in the prefix pass of the Montgomery inversion, and a
  denominator is zero exactly when the two x-coordinates are equal, so that
  product is zero exactly when some member's addition is the exceptional one.
  The batch tests the product once with three `orr` and a `cbnz`, four
  instructions charged in the batch row above. A zero product pays one 600-unit
  unbatched inversion action per member of the batch, which redoes each
  member's inversion and so locates the walk whose two x-coordinates are equal;
  that walk's trail is then abandoned, and every other member is still at its
  own input point, so its step is simply redone and nothing is lost. The event
  is the exceptional-addition row of section 4, below 2^-114 in total, so its
  expected cost is below 2^-89 units and it is already inside the random extras
  of the budget below.
- Length, cap, meter, queue, loop (10, was 12): the 2-limb length (L_max >
  2^64) by `adds` and `adc` (2), because it is one of the values held in
  registers across the batch loop; cap test by `cmp` of the high limb with 20
  and `b.hs` (2), which is length ≥ 20·2^64 exactly; meter: `subs` on the
  processor's block remainder, `add` on its operation count, `b.lo` to the
  block refill, where the global counters and the phase-1 limit are checked
  (3); the queue pointer by `ldr` with post-index (1); loop counter `subs` and
  `b.ne` (2). Total 10. The two instructions that go are the `ldp` and `stp`
  that reloaded and rewrote the length every step: the length changes only on
  a step, and a walk that leaves the inner loop writes it back inside the
  action that took it out (report, doubling, fallback, abandon), so no reader
  sees a stale copy. The meter still charges 1.099 units per step before the
  step runs, and the refill still checks the global counters and the phase-1
  limit in the same milli-unit accounting.
- Batch inversion (4): an upper bound, not a counted sequence. The reference
  price's 3M amortisation pays the batch's 3(L − 1) prefix-product and
  back-substitution multiplications and their traffic; the batch's single
  inversion, I < 2^17 instructions, is outside it, and divided over the
  L ≥ 2^15 walk steps that share the batch it is at most I/2^15 < 4 per step.
  Nothing the reference price already pays is charged here.

The sequences total 15 + 8 + 13 + 13 + 5 + 12 + 19 + 0 + 10 + 4 = 99
instructions, and the step is charged 99. Rare paths are priced as their own
actions: the 2-cycle order test, the flag write, the clearing of the two
history slots and the escape are in the doubling price, the fallback has its
own price, and so do reports and traversal steps.

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
| ordinary step | 1.099 | reference addition plus 99 itemised instructions |
| traversal step | 1.123 | ≤ 24 more instructions to copy the remembered member: ≤ 123 |
| doubling | 2.0 | batched affine doubling: 7M + 6A + 150 ≤ 1,931 instructions, unchanged |
| §1.1 fallback | 0.5 | one M, two A, compares: < 400 instructions |
| random 64-bit word | 0.01 | |
| sampling attempt | 0.05 | |
| trail start [a]P + [b]Q | 2^13 | two Jacobian scalar multiplications at ≤ 512 operations × (16M + 16A) each, an addition and an inversion: < 4.8·10^6 instructions |
| store probe | 0.05 | 39 instructions itemised in section 3, charged 50 |
| report | 40.4 | 40 for route, message, entry write and node tests (< 2^15 instructions), plus 8 probes at 0.05 |
| overflow insert scan | 2^41.45 | worst case 6B + 20 instructions |
| overflow lookup | 2^42.87 | worst case 16B = 8,005,983,135,385,520 instructions = 2^42.86422 units; charged ⌈1000·2^42.87⌉ = 8,038,146,442,837,048 milli-units |
| store initialisation | 2^16·B buckets × 16 instructions = 2^58.865, once | |
| table | 2^20 starts plus 5 φ-twists each: < 2^33.1, plus sampling | |
| phase-2 or re-walk operation | 600 | unbatched: an inversion (< 2^17) plus 7M, coefficient updates mod r and control, far below 600,000 instructions |
| small-batch inversion | 600 | one inversion, < 2^17 instructions = 131 units; 600 is a loose upper bound |
| solving | 600 | |
| verification | 2^13 | |

**Budget.**

- **Phase-1 operations:** 1.099·(S_ops + 2^60) < 2^126.91537.
- **Hard extras:**

  | item | units |
  |---|---|
  | meter remainders | < 2^61 |
  | phase 2: 2^20·(L_max + 65)·600 | 2^97.551 |
  | re-walks: 8·2·((L_max + 65)·600 + 2^13) | 2^81.551 |
  | reports: D_cap·40.4 | 2^69.167 |
  | store initialisation: 2^16·B × 16 instructions | 2^58.865 |
  | table and verification | < 2^33.1 |
  | overflow insert and lookup: D_cap·(6B + 20 + ⌈2^42.87·1000⌉) instructions | 2^107.158 |
  | **total** | **< 2^107.16** |

- **Random extras (expected):**

  | item | units |
  |---|---|
  | doubling surcharges: 0.901·2^-22.58493·(S_ops + 2^60) | < 2^104.044 |
  | traversal surcharges: 0.024·64·2^-38.00·(S_ops + 2^60) | < 2^89.40 |
  | starts: < 2^62.8·2^13 | 2^76 |
  | words, attempts, small batches, fallbacks | < 2^70 |
  | **total** | **< 2^104.044** |

- **Headroom:** 2^126.92 − 2^126.91537 − 2^107.16 > 2^112.74.
- **Early meter halt:** probability < 2^104.044/2^112.74 = 2^-8.70 = 0.0024
  (Markov); exactly, the allowance is ⌈1250·random/3⌉ = 2^112.7466 units, so
  the ratio is at most 3/1250 = 0.0024.

### 5.1 Exact total and the two-decimal claim

All ledger quantities are carried as integers in milli-units (10^-3 units; one
instruction is one milli-unit), so every product is exact. With
S_int = ⌈S_ops⌉ = 145,993,368,281,645,653,992,586,223,245,627,523,266,
D_cap = 16,396,253,461,269,537,204 and B = 500,373,945,961,595:

    phase 1            1099·(S_int + 2^60)                       = 160446711741528573739119320080507572895958
    hard extras total (section 5, < 2^107.16)                     =                   181253151995776194805495440951804472
    random extras      ceil of < 2^104.044                       =                    20909018291647366215258114967470080
    Markov allowance   ceil(1250·random/3)                        =                  8712090954853069256357547903112533334
    W (milli-units)                                               = 160455605085635422584570483123851637233764

W in units is W/1000 = 2^126.915441793…, whose 2-decimal ceiling is 126.92.
The exact integer certificate is W^100 ≤ 1000^100·2^12692, and W^100 >
1000^100·2^12691, so 126.91 does not suffice. The two sides of the first
inequality are the integers

    W^100 =
    34315415419195268485965498516309260820992246847734636687485113182485493793701528077662710968430703777596227362767192268247413163973066588999947089946464643097556965299171343600177463177588731750036597214170204922294874275775348082631056276245232295771688159287423858702599969773443559604798590504135582690430667000867527814327591443420014749124249735161644277983945215642685517402228846453673710028576186954066488487059748087544972073052722199641534230588524937223578249976580385955878082498235614371874719235281047585789502173805984782089454570927390630872486005197912713442576239241470746848372893687506314952716277817933877329257704836784599096724926629715418976200726888677458022861190478561210150758746812288603252085281525616624839977568619978986643830139498523218958749729939715577608543087585774586261550744132744245495470221390390682266701079423856169271116666934465192250734037828721897078193888214183849437985886963975171991995410778334333733646469112261543511723792772263668362720546460728010740616512392224458370555040953275127705749535397624765575158327009563442542144234347264302534881432016017557090251502501837695349670429721581720820855777954142959896357712272254873357700642819237041656488944170580605278788057042685706424385816773435388859196646641779921398278747912832369110149271267863205645178819967718640500094410486594748561404994860664756274998716511030696395866399418388948522159372473591378667600553976330925312082987797541777213239538155441861089593016677567824156694616943598308144647679997107102882243748884145152780604642182484871168922050798388526820455772711767300254394688080492555520807136473368005646761313123261940513367607546346921569506323614345683595428088077153903032223082385765923842534680739077389121338169619926056387174798942770955801443191841299526066002577036386356443264920335916088535721715581773111128083070023398279371090421632130809800617799101242431739262587890051826411029497913974342718200690239758133248558367623603412840831166994944374567894896788847143940909754568390178011802951848792772725641839117674998534841863751195405307349891017510024427092364722171503963529608750988641286992348727040681900080146947733717444659212194258729435244722082920959741097501249856189086977167230225113218375467917091613365638909583051020996828606048475841529028996551594429313976847545181610956164005799630238699776258491171475462994303644384578777467663744327421578156006496041688890871472285979495100750216031384369149657619128320296380114510997996188068627837510915357056921788454429046654532417758386347228670690816656419203375271085696124232147937622483188571593507682225059929663728926456830444239172275403889702279748483139100996080426871103284342324845185021078574074404644945540038247158254910425661344959340782655047432785958220641416412319336597255020947041428994526840027978916280969291687788300344193353934098764281547666580390238513663972681719507135277310237458916317241867657520589286668220978757278889168256293934042820149874576441069085830060230681213809831667700273433804435391669346173827809327297586189299903827766388406703375248663715751235801921256395347519614418145106624353768178603648533686228774319956459793859456404657717605652182811539532541258402745989530960707359451982248563224059849287911791885768030234225580774328316458257468235347370194304168706093158696399952199606920879221107681622114950074485868592999549909280666244179230389763160444840612490089131718348534954883115942764849154782039232018398683840131771292397903864928546163365818305100805273640956410014832392010430260633615144905189971910484017241222593535524370187267051182684018754686309048563883161856033218589247754568999840450670996490203768875631446986694124190975087830019247172504078365195450134682148798899250325499639336394056191561940594993857952598171004647812965085119695669878242913136705172574558380237913227600012828051344167524351571506730488016406792312667980026431101272638004588459311980434851009262196673371176468853290369412631367047344755535618459978839873345475300258079943662918843123407509727185608377057633838224981699234790348855866349741092913729765376

    1000^100·2^12692 =
    47065748224901432086212995484307886776516406449820927741319133494548135278536570773475914405417387349124214974714295160264155160997647263140561209192641044609538335736591122050398464452223783157104272870125176030373432610307808103062281904284024708787674447021474794378814596589160769878969574383108328416569557263243690016175018741615058968695408469705698748723998639412868193220329968870118929822813153068163218395771989178867985254894861657820507814200047452429178723104649530347614307972260369191470017662367529453210694405388017239074254540349999253580360508258654685495765342769645270602984305730009025934271050408628403935772630715211449658531069343331953758460413064673488323106526714064343927434876406273220938537663334550650369515963848592136280909093106428396974954368119925735055807395961954539994896708116205693184176177769529555187523266297948417973604162049396082638328537470332035165888192660597103136494004141520758537098756149212664568452454306814015244881951378103692243928037300536807622094103214289419945601831762541434958404804545841958763700608444747602031718672924075877402655444249156451463026821069623078844997796200820317196076885327455067556987078523147353576836493710935151579903355837020358343687417968604354220485529274929811087867808013041587179999631880386074561598327743907079101772202647336568082800052169556538858828741615537576744055610591600987181610362763160131909664055326175021456458892177937333931703314518651997534491443507438758172427051445793167891690955531078903792815597132833058218667766681463777665652174453346474236030215575322892606626398150605106422661846583697041362769265597686907202647482597093472399806444610892153710279228350639571663857562320631310771794055475000795381234695706803714972387749302873384798019437302236496207485463697617245902889616048698904143290584503063852625957513081784410480298079480101961951641471810754455508455120208759390492955767873447934116072298977710634727683908848396796762751446251671335519354269744876578281573600872099220981506344260994314049631958811274135079228960452837287280207317450690335237124717645834873553166303203209449372227889918933591941164186983139145919724978476701525720211125142820313591271317597334108086114328656768478878290419023558811235658703758911108065454191985891431087576935165000589103854024823942895619984585486173480869947867089675175728661831470656095935343506949574147861495670489853359731693934108982546196066552022788242195592559259948220891322769134789137200607404883722902695061908085105552803438514931686127515872339180916612509783248491273789250917445532816672785753702560233969283465760653377506346148359470134604994575116113966802451242108895455336514656085760547129440789767287258643758833478145401663034184687284045376617390568930654810753511820851980050047033797060771651624863159074496639006068104517390944178076683869763978113396524288570581839675649282990864055811343735163416046741205395612315628277500021087009105677509793955569423524057793488686664414775365526895057714110034281676432275694888060500801823372453297555840904656182962008167173283805752043805334038407109735880952016646687518388510164712450999821881212442838789880549629929447518782268253286974337737351167234738259999385092228981326669817319711444744221043059254136951054723111403251327931674024400865544571580739890669119844959071953369201876517558119364213850529174074795919688615757669219566940556121310086039706421904770793306272169789482096107140196754268393420735156458432673439218859659967225192673080051611860234247718182846041121322401797738957926784108568264564008954012716131003231601902289532187982962589466972677595551423746256108804294953076409001214744845632979717456226081526800711552867744156972348339966933692015076143354466766096259012333032364197057469697533724748956095377593808036349114786512896000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

Both carry 4,121 decimal digits, and 1000^100·2^12691 carries 4,121 as well;
every comparison here is exact integer arithmetic. The ratio
W/(2^126.92·1000) is 2^-0.004558, so the claim has 0.0046 bits of room above
the ledger's own exact total. The score claimed is 12692/100.

**W ≤ 2^126.92 on every run.**

**Preprocessing:** table < 2^33.1, node set-up 2^26, store initialisation
2^58.865, and the first 2^20 starts 2^33. The expected total is < 2^58.87,
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
margin covers it: with all random charges 45% above their bound 2^104.044,
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

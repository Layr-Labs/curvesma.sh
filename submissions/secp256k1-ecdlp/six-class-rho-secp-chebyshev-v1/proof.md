# Six-class parallel Pollard rho on secp256k1: counted AArch64 kernels, second-moment meter control

**Provenance.** This package revises public submission ab043cf7 (`six-class-rho-secp-counted-v8` directory,
the 126.338 bound), which is its declared dependency. It changes only how the random part of the work meter
is controlled and the H1 floor that this frees: the Markov allowance and its 0.00416 success deduction are
replaced by Chebyshev allowances under an explicitly stated second-moment heuristic H2 (section 4), and the
H1 floor p0 at S0 is lowered from 0.3953 to 0.3901. The walk, the class key, the canonical-index rule, cycle
handling, the merging lemma, the phases, the report store, the slowdown allowance κ = 1.02, every price
(the counted 111 + 523 instructions per walk, the 118 analytic rows, the 0.956 doubling and the 0.776
traversal), every hard row and the generated assembly are those of ab043cf7. That line descends from
88e72314, 061e00ca, 30f6f34e, 98e7729c, b52f8806, a9f1a6fe, 3682292, 758d8df6 and 189ebc42; all are credited.
The complete argument is restated here; no prior acceptance or score is used as a premise. Section 0 lists
every number that moves, section 4 gives the new success chain, section 5 the exact integer ledger, and
`certificate.py` proves the resulting bound in exact integer arithmetic both ways.

**Claim.** time_log2 = 126.3173 = 1263173/10000. The classical cold attack recovers x from a uniform
secp256k1 point Q=[x]P with probability at least 0.39009 under the explicit H1 and H2 heuristics in
section 4. On every run total work is at most 2^126.3173 target-group-operation equivalents because the
global meter stops before that bound, and the exact integer ledger totals
W_ledger = 105995162764622566518652505820930400478098 milli-units = 2^126.317266 units (section 5),
every row of which is charged. The unit is 1,000 scalar AArch64 instructions (C = 1,000, fixed by the
organizer); no reference operation is charged for ordinary steps, see section 5.0. Peak memory is below
2^72 bytes (2^71.64 for the store alone), preprocessing is capped at 2^59.0 units and included, and
there is no nonuniform advice. This is an analytical, impractical generic attack, not a demonstrated key
recovery.

## 0. Change from ab043cf7

| item | ab043cf7 | this package |
| --- | --- | --- |
| algorithm, prices, store, hard-row definitions | – | unchanged |
| ordinary step / doubling / traversal | 0.752 / 0.956 / 0.776 | unchanged |
| slowdown allowance κ of H1 | 1.02 | 1.02 |
| H1 floor p0 at S0 | 0.3953 | **0.3901** |
| S_ops = κ·S0 | 2^126.74079 | **2^126.72840** |
| control of the random meter charges | Markov on their sum, allowance E/0.00416 = 2^118.910 | **Chebyshev (doublings, traversals) under H2, Markov ×2^30 (starts, small charges); allowances 2^111.6524** |
| success lost to an early meter halt | 0.00416 | **4·2^-30 < 4·10^-9** |
| phase-1 cost | < 2^126.3296 | **< 2^126.31721** |
| hard extras | 2^107.15551 | 2^107.15550 (the refill row follows W) |
| ledger total W_ledger | 2^126.337997 | **2^126.317266** |
| W (the claim) | 2^126.338 (three-decimal ceiling) | **2^126.3173 (four-decimal ceiling)** |
| success floor | > 0.3911 | **> 0.39009** |

What changes, and why it is sound:

1. **The meter's random charges are controlled by their second moment.** ab043cf7 bounded only the
   expectation E = 2^111.001 of the random charges (doubling and traversal surcharges, starts, small
   charges) and reserved E/0.00416 = 2^118.910 units, which costs 0.0084 bits of W and forces a higher p0.
   Here the doubling and traversal counts get Chebyshev allowances (2^111.28 and 2^109.39 units, tails
   2^-30) under the second-moment heuristic H2 of section 4, Var ≤ 2^80 × (mean bound). Starts and small
   charges keep Markov at 2^30 times their expectations.
2. **The H1 floor p0 = 0.3901.** Success ≥ p0 − 2·10^-7 − 4·2^-30 = 0.39009979… > 0.39 (section 4).
   ab043cf7's margin above 0.39 was a choice. H1's heuristic margin stays the 2% allowance κ − 1, against
   measured toy slowdowns of +1.13% ± 0.17% (q25) and +0.22% ± 0.20% (q29).
3. **Four-decimal claim:** W_ledger = 2^126.317266. certificate.py proves that 126.3173 holds and that
   126.3172 fails.
4. **Numbers that follow W or S_ops:** the meter cap, refills (< 2^86.3174), the refill row (2^81.259), and
   the operation (< 2^126.7285) and start (< 2^113.3173) counters. Phase 1 is now charged S_int + 2^61
   operations (was + 2^60), which covers the block overshoot < 2^60.41 of section 3. The small charges are
   bounded by 2^71 (was 2^70).
5. **A toy statistic for H2:** `doubling_count_dispersion` (section 6).

### 0.1 The analytic rows

The 118 analytic bookkeeping rows are unchanged from ab043cf7, which re-derived them against the shipped
assembly and the pseudocode of section 2. Every row is the sum of an explicit instruction sequence, printed in
section 5.3, with no allowance on top. They cover what each class arrival does outside the two counted loops:
the class key |y|, the truncated canonical index, the distinguished-point and routing fields, the history and
checkpoint tests, the meter and the arrival loop. The counted loops already load and store everything they
consume (pass 2: 19 loads and 6 stores; pass 1: 5 loads and 4 stores), so no row repeats their traffic. The
rows are F1 24, range test 8, decision 13, sign/p − y/K 17, h/DP/t/address 7, arrival traffic 5,
fingerprints/history/checkpoint 26, zero-denominator test 4, length/cap/meter/loop 13, and inversion share 1:
118 in total.

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
  the rate is < 2^-56 per step if βx mod p is equidistributed. This affects
  only the meter's early-halt probability, never correctness.

Part A compares this rule, emulated at limb level, with the exact rule on
random production inputs and on adversarial near-tie inputs (two of the three
values, or a value and 0 or p, within 2^k for k = 150..249). The initial run
asks for 30,000 random and 10,000 adversarial inputs within its time limit and
reports how many ran (random_done, stress). The offline run (section 6)
covers 2,000,000 random and 1,000,000 adversarial inputs, with 0 mismatches.

## 2. Walk (unchanged from 189ebc42)

Parameters: R = 2048; distinguished probability θ = 2^-64 (limb 1 of K is
zero); M = 2^20 walks; L_max = 20·2^64; checkpoint period 64.

- **Table.** c_h and d_h are uniform in [0, r) (fresh coins, rejection sampling,
  every attempt metered). M_h = [c_h]P + [d_h]Q, redrawn if O.
  Tab[h][t] = [u_t^{-1}]M_h = (−1)^s·φ^{−j}(M_h), with coefficients
  u_t^{-1}·(c_h, d_h) mod r.
- **Step.** F(C) = class(c(C) + M_h) with h = K mod R. The walk keeps any member
  P and sets P ← P + Tab[h][t]. Since c + M_h = [u_t](P + Tab[h][t]), this is F.
  Each ordinary step is one batched affine addition, implemented and counted
  in section 5, plus the class bookkeeping rows of section 5.3.
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
  fallbacks, reports, store probes, overflow inserts and lookups, meter block
  refills, table build, store initialisation, re-walks, solving and
  verification. If a charge would exceed W = 2^126.3173 units, the attack halts
  with failure. So W ≤ 2^126.3173 on every run, and no loop needs its own limit;
  rejection-sampling loops in particular are bounded by the meter.
  - Preprocessing (all charges before the first walk step) has its own cap of
    2^59.0 units. If a preprocessing charge would exceed it, the attack halts
    with failure.
  - **Counter widths.** Every counter is stated with the largest value it can
    reach and the word that holds it, and each fits with room to spare.
    - The meter counts milli-units and is a 192-bit counter (three 64-bit words,
      updated by ADDS/ADCS/ADC carry chains). Its largest value is 2^126.3173
      units < 2^136.28309 milli-units < 2^192: 55 bits of headroom, and no
      capacity test is needed at 192 bits.
    - The global operation counter of phase 1 counts steps, doublings and
      traversal steps and stops at S_ops, so it reaches at most
      S_int + 2^61 < 2^126.72841 + 2^61 < 2^126.7285: 127 bits, held in a 128-bit
      pair of 64-bit words (ADDS/ADC on the pair). Its stop test compares the
      128-bit pair with S_ops, which is below 2^127, so the test can fire.
    - The start counter counts trail starts. Every start is charged 2^13 units,
      so it can never pass W/2^13 < 2^126.3173/2^13 = 2^113.3173: 114 bits, held
      in a 128-bit pair. (The expected number of starts is the 2^62.73 of the
      start row of section 5; the counter is sized for the worst case,
      not the expectation.)
    - The preprocessing cap is checked on the meter's milli-unit scale:
      2^59.0 units = 2^68.966 milli-units < 2^69, held in a 128-bit pair and
      compared by a two-word SUBS/SBCS/CSEL sequence. The preprocessing account
      is a separate 128-bit pair, so the cap and the global meter never alias.
    - The block-refill counter, and each processor's private count of the blocks
      it has taken, holds at most W/2^40 + 2^20 < 2^86.3174 refills: 87 bits, a
      128-bit pair.
  - Each of ≤ 2^20 processors reserves blocks of 2^40 units by atomic update on
    the 192-bit meter. The refill itself — an atomic 192-bit add of the block, a
    store of the new remainder, a load of the phase-1 limit, the operation count
    and the two limit tests — is charged 30 instructions = 0.03 units per
    refill, the hard row of section 5. At most W/2^40 + 2^20 < 2^86.3174 refills
    occur, so that row is < 2^86.3174·0.03 < 2^81.26 units. Unused remainders
    count as spent: < 2^61 units, also in the hard extras.
- **Phase 1.** Runs until a global operation counter (steps, doublings and
  traversal steps of phase-1 trails) reaches S_ops. Every reserved block is
  completed, so at least S_int and fewer than S_int + 2^60.41 < S_int + 2^61 operations are
  performed.
- **Phase 2.** Finishes the open trails and starts no new ones.
- **Caps.** At most 8 processed matches; a ninth is a failure. Report capacity
  is D_cap = 2^63.83 records, taken as the integer
  ⌈2^63.83⌉ = 16,396,253,461,269,537,204 in every product below. These bound
  memory and re-walk work.
- **Batches.** Ready walks are split into batches of 2^15..2^16 per Montgomery
  inversion. A batch of fewer than 2^15 ready walks pays its inversion as a
  separate 600-unit action.
- **Store.** Reports are stored in 2^16 independent nodes, each a separate
  machine with its own 64-bit address space. A report is routed by bits 128..143
  of K (limb 2) to its node and addressed within the node by the 50-bit field at
  bits 144..193 of K. The two fields are disjoint from each other and from the
  index bits 0..10 and the DP bits 64..127. Each node holds
  **B = ⌈2·D_cap/2^16⌉ = 500,373,945,961,595 buckets** of 104 bytes (32-byte key
  K, 32-byte a, 32-byte b, 8-byte tag and flags) — 2^55.531 bytes — together with
  an overflow list of B 8-byte bucket addresses of 2^51.831 bytes; with both, a
  node's store is 112 bytes per bucket = 2^55.64 bytes < 2^56, and the 2^16
  nodes together hold 2^71.64 bytes < 2^72.
  - **Home bucket.** The 50-bit field at bits 144..193 of K, reduced modulo B by
    a two-word multiply-shift (MUL, UMULH, LSL, ORR, LSR, ADD: 6 instructions,
    inside the report price).
  - **Probe rule.** Lookup and insertion use linear probing with a fixed maximum
    depth of 8 buckets. Probe i examines bucket j = (home + i) mod B:
    - if the bucket is free, the report is inserted there;
    - if the tag and the four key words match, the report is a match — the
      collision the attack wants;
    - after 8 occupied non-matching buckets the report enters the overflow list:
      the list is scanned for a match, a free bucket is taken by scanning
      cyclically from (home + 8) mod B, the report is written there and its 8-byte
      address is appended. No report probes more than 8 buckets of the table.
  - **Probe itemisation** (39 instructions, charged 50 = 0.05 units):

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

    The 11-instruction margin covers the two conditional branches and the loop's
    register moves; every charge here is an upper bound.
  - **Overflow list and its capacity.** The list is a per-node array of B 8-byte
    bucket addresses with a counter, and every overflowing report is written
    into one distinct free bucket of the table with one address appended for it.
    The table has B buckets, so the list can never hold more than B entries and
    cannot overflow; no probability enters its capacity. The insertion scan is
    cyclic: it visits buckets (home + 8 + i) mod B for i = 0, ..., B − 1, as the
    two linear segments itemised below. The node is declared full only after a complete pass of all B
    buckets finds no free one, which happens exactly when the node holds more
    than B reports; that is the node-overflow event bounded in section 4 by a
    Chernoff argument, and the attack halts with failure. A partial scan is
    never treated as full, so no report is lost at load below one.
  - **Overflow insert, worst case.** Start (9): j0 = (home + 8) mod B by ADD,
    CMP, CSEL; 13·j0 by LSL #3, LSL #2, ADD, ADD; segment end B and segment flag
    by two MOV. Scan two linear segments, j0..B − 1 then 0..j0 − 1, at 6
    instructions per bucket (flags LDR, CBZ to the write, scaled address ADD,
    index ADD, CMP, B.LO). Between them (6): CBNZ on the flag, taken only after
    the second segment, when no bucket is free (the node-full halt); MOV flag;
    MOV address and index to 0; MOV end to j0; B to the loop. Write the 104-byte
    entry (13 STR) and update its tag, visit word and pointer (3); append the
    address and increment the counter (4). At most 9 + 6B + 6 + 20 = 6B + 35 =
    3,002,243,675,769,605 instructions, 2^41.45 units.
  - **Overflow lookup, worst case.** Scan the list: 16 instructions per entry
    (address LDR, tag LDR, CMP, B.NE, four key LDR, CMP and 3 CCMP, B.EQ,
    pointer ADD, end CMP, B). At most 16B = 8,005,983,135,385,520 instructions,
    which is exactly 8,005,983,135,385,520 milli-units = 2^42.86422 units. The
    charge is that exact integer, so no rounding step stands between the
    instruction count and the price.
  - **Charged worst case.** Each report performs at most one lookup and one
    insert and probes at most 8 table buckets, so the overflow work over all
    reports is at most D_cap·(6B + 35 + 16B) = D_cap·(22B + 35) milli-units =
    1.80494·10^35 milli-units = 2^107.154 units. That is the ledger's overflow
    row, in the hard extras; the 8 probes are the 0.4 units of the report price
    40.4.
  - **Load.** Under H1 the routing field is pseudo-random, so a node's share of
    the n ≤ D_cap reports is Binomial(n, 2^-16) with mean < 2^47.83 against
    B = 2^48.83 buckets. Chernoff gives < exp(−2^45) per node, and < exp(−2^44)
    over all 2^16 nodes. The message, routing and node-side work are in the
    report price.
  - **Node set-up and bucket initialisation.** Node set-up (process start, table
    allocation, routing tables) is charged 2^20 instructions per node, 2^26
    units in total, as preprocessing. Bucket initialisation is charged 16
    instructions per bucket, 2·D_cap·16 = 2^58.86446 units, once, also as
    preprocessing. The overflow list needs only its counter initialised (2
    instructions per node), because every entry is written before it is read.

## 4. Success probability

N' = N·R/(R − 1). T_c is the number of class arrivals of all trails, up to the
first arrival at a class already visited by another trail. It includes trail
starts, steps, doublings, traversal steps and repeated arrivals inside
fruitless cycles; nothing is excluded. This is exactly what the toy measures.

**H1 (heuristic, at the budget).** P(T_c ≤ S_ops) ≥ 1 − exp(−(S_ops/κ)²/(2N'))
with κ = 1.02 and p0 = 0.3901. H1 also includes three further claims:

- the partition, DP, node-routing and bucket fields of |y| are pseudo-random
  functions of the class: bits 0..10, 64..127, 128..143 and 144..193, four
  disjoint fields;
- excluding the detected fruitless cycles of length ≤ 64, a trail makes a
  non-fruitless return with probability ≤ L_max²/(2N);
- the canonical twist index of a new class is uniform and independent of its
  partition index.

The 2-cycle rate is 1/(6R) per step. Other cycles of length ≤ 64 have rate
≤ Σ_L Σ_b S2(L, b)·R^-(L−b) = 2^-19.9954 < 2^-19.99, summed over L = 3..64
and b = 1..⌊L/2⌋, where S2(L, b) counts the partitions of the L steps into b
blocks of size ≥ 2 (every partition index used at least twice; for example
1, 1, 3 and 15 for (L, b) = (3, 1), (4, 1), (4, 2) and (6, 3)). Cycles of
length 65..200 have rate < 2^-206 per step. The expected doubling rate per
operation is therefore ≤ 1/(6R) + 2^-19.9954 = 2^-13.5681 < 2^-13.56.

The values are:

    S0 = √(2N' ln(1/(1 − p0))) = √(2N' ln(1/0.6099)) = 2^126.69983
    S_ops = 1.02·S0 = 2^126.72840,  S_int = ⌈S_ops⌉ = 140945347563975464246948179984935587492

H1 gives ≥ 0.3901; the pure random-mapping model without any slowdown gives
1 − 0.6099^(1.02²) = 0.4022. The
toy censuses of section 6 (65,641 to 96,675 trials per curve) measure the slowdown of this
walk against the model as +1.13% ± 0.17% on q25 and +0.22% ± 0.20% on q29, falling with N, and
the earlier runs of section 6 (−0.9% to +1.5%) agree; the allowance κ − 1 is 2%.

**Phase 1 covers the budget.** Phase 1 performs ≥ S_int ≥ S_ops operations, and each is
one class arrival. So if T_c ≤ S_ops, the first inter-trail collision happens
in phase 1.

**Random meter charges.** Phase 1 performs a fixed number of operations, each charged the step price
0.752 in the phase-1 row. Four kinds of charge are random:

- doublings: each costs 0.956 − 0.752 = 0.204 units more than a step. Their count D among the at most
  n = S_int + 2^61 phase-1 operations has E[D] ≤ μ_D = (1/(6R) + 2^-19.9954)·n ≤ 2^-13.5681·n,
  μ_D = 2^113.1603 (rates above; the certificate uses an upper rational bound of 2^-19.9954);
- traversals: a checkpoint detection is followed by at most 64 traversal steps at 0.024 extra each, at
  most 1.536 units per event. The event count V has E[V] ≤ μ_V = 2^-19.99·n = 2^106.7384;
- trail starts: 2^13 units each, expected charge ≤ 2^13·(S_int/2^64 + 2^20) = 2^75.7284 units;
- words, sampling attempts, store probes, small-batch inversions and §1.1 fallbacks: expected < 2^71 units (fallbacks
  0.5·2^-56·n < 2^69.73, probes 0.4 per report < 2^62.5, the rest smaller).

**H2 (heuristic, second moment).** Var(D) ≤ 2^80·μ_D and Var(V) ≤ 2^80·μ_V.

Support. Write D = Σ_i Z_i over phase-1 operations, Z_i ∈ {0, 1}. In the random-mapping model of H1, a trail
that has not met another trail visits fresh classes. Its partition, twist and fingerprint fields are
independent of every other unmerged trail's, so Cov(Z_i, Z_j) = 0 when i and j lie in different trails that
do not merge. Within one trail H2 allows arbitrary dependence: Σ_j Cov(Z_i, Z_j) ≤ Σ_j E[Z_i Z_j] ≤
(L_max + 65)·E[Z_i], since a trail has at most L_max + 65 < 2^68.33 operations. Two trails that merge
repeat each other's operations from the merge until their common distinguished class, at most L_max + 65
operations. With K_i the number of trails merged with the trail of operation i, the merged pairs add at
most (L_max + 65)·E[Z_i·K_i]. H1 puts the expected number of inter-trail merges in all of phase 1 at about
S_ops²/(2N') ≈ 0.51, so E[Z_i·K_i] ≤ 2^10·E[Z_i] holds with an enormous margin, and
Var(D) ≤ 2·(L_max + 65)·(1 + 2^10)·μ_D < 2^79.33·μ_D < 2^80·μ_D, and likewise for V. Shared randomness
does not add a μ²-size term (law of total variance): the table M_h, Q and the coins enter the doubling
and traversal rates only through partition, twist and fingerprint coincidences of fresh classes (the
2-cycle needs h repeated with the negating twist, and the short structural cycles are identities in the
units and table indices, not relations among table values). So E[D | table] has the same bound μ_D for
every table, and Var(E[D | table]) adds no term beyond the bound.
The cushion is very large: within a trail the events are local, since detection uses the last two arrivals,
or a checkpoint at most 128 + 64 arrivals back, followed by at most 64 traversal steps. Independent rare
events would give Var(D) ≈ E[D], which is 2^-80 of what H2 allows. The toy reports the per-trial index of
dispersion of doubling counts (`doubling_count_dispersion`, section 6) as a check of the same structure at
toy scale, when toy trials complete within the run's time budget. H2 is a heuristic in the same sense as H1. It is not a theorem about the deterministic walk.

**Allowances.** The meter reserves, on top of phase 1 and the hard extras:

| charge | control | allowance (units) | tail |
| --- | --- | --- | --- |
| doublings | Chebyshev under H2: t_D = ⌈√(2^80·μ_D·2^30)⌉ = 2^111.5802 | 0.204·(⌈μ_D⌉ + t_D) = 2^111.2832 | Var/t_D² ≤ 2^-30 |
| traversals | Chebyshev under H2: t_V = ⌈√(2^80·μ_V·2^30)⌉ = 2^108.3692 | 1.536·(⌈μ_V⌉ + t_V) = 2^109.3921 | ≤ 2^-30 |
| starts | Markov | 2^30·2^75.7284 = 2^105.7284 | 2^-30 |
| small charges | Markov | 2^30·2^71 = 2^101 | 2^-30 |
| **total** | | **2^111.6524** | **4·2^-30** |

The Chebyshev rows use P(D ≥ μ_D + t_D) ≤ P(D − E[D] ≥ t_D) ≤ Var(D)/t_D², valid because μ_D ≥ E[D].
If none of the four tails occurs, the random charges stay inside their allowances, and the meter does not
halt before phase 1 completes its S_int operations and phase 2 its hard-row budget.

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
| preprocessing cap | < 2^-30 |
| meter halting early: one of the four random-charge tails above (H2 for two of them) | ≤ 4·2^-30 < 3.8·10^-9 |

The events other than the early halt sum to < 2·10^-7, as in ab043cf7.
**Success ≥ 0.3901 − 2·10^-7 − 4·2^-30 = 0.39009979… > 0.39009 > 0.39** (exactly
8180985679571/20971520000000 in `certificate.py`).

## 5. Cost

### 5.0 What is charged

COST_MODEL.md prices a computation as H reference operations plus R other
instructions at R/C, with C = 1,000 scalar AArch64 instructions per unit, and says
"Faster arithmetic and lookup algorithms are welcome; justify their actual work" and
"Do not charge both a reference operation and its internal instructions". This package
therefore charges **H = 0 reference operations** for ordinary steps. The affine addition
is not bought as a unit: it is implemented in AArch64 assembly (emitted by a generator, with a C
specification), and its executed instructions are counted and divided by C. The organizer's 6M reference unit (3M for the
addition formula, 3M for batch amortisation, S = M, A = 0) is the same computation with
every non-multiplication instruction set to zero; here the multiplications cost what the
assembly costs and the additions, subtractions, loads, stores and loop control
are charged too. Nothing is charged twice.

### 5.1 Kernels (`kernels/secp_walk.c`: C specification; `kernels/gen_loops.py`: the counted AArch64)

p = 2^256 − C0 with C0 = 2^32 + 977 = 4,294,968,273, so 2^256 ≡ C0 (mod p). A field
element is four little-endian 64-bit limbs. **Weak** means any value in [0, 2^256);
**canonical** means [0, p). Walk points and table points are always canonical.

- **fe_mul(a, b)**, a and b weak, returns a weak r ≡ ab. Row-wise schoolbook product
  (16 MUL, 16 UMULH, carry chains of ADDS/ADCS: 12 instructions for row 0, 17 for each of rows 1 to 3),
  then z = L + H·2^256 ≡ L + H·C0:
  H·C0 < 2^256·C0 is five limbs (4 MUL, 4 UMULH), s = L + H·C0 has a fifth limb
  s4 ≤ C0 (so s < 2^256 (1 + C0)), and s ≡ (s mod 2^256) + s4·C0 with s4·C0 ≤ C0² < 2^65
  (two limbs w0, w1, w1 ≤ 1). The sum t = (s mod 2^256) + w is below 2^256 + 2^65; if it
  carries out, t − 2^256 < 2^65, and adding C0 once more stays below 2^256, so no second
  carry exists. Every limb addition of a carry bit into a limb that the product has not
  yet touched is an addition to a partial sum of ab < 2^512, so it cannot overflow. In the
  assembly the whole multiplication is 91 instructions: 63 for the product, 28 for the reduction.
- **fe_sqr(a)**: the six off-diagonal products, doubled by a one-bit shift across the
  limbs (2·offdiag ≤ a² < 2^512, no bit is lost), plus the four diagonal squares; then
  the same reduction. 10 MUL + 10 UMULH instead of 16 + 16; 71 instructions in the assembly.
- **fe_sub(a, b)**, a weak, **b canonical** (nine instructions: SUBS, 3 SBCS,
  CSEL of C0 on borrow, SUBS, 2 SBCS, SBC): a 4-limb SUBS/SBCS chain; on borrow-out
  (a < b ≤ p − 1) the wrapped value a − b + 2^256 is at least 2^256 − b ≥ C0 + 1, so
  subtracting C0 (the same as adding p modulo 2^256) cannot borrow again and the result
  is a − b + p ∈ (0, p). Without borrow the result is a − b. The result is ≡ a − b and
  is canonical whenever a is canonical.
- **fe_canon(r)**, r weak (ADDS, 3 ADCS, 4 CSEL, eight instructions): r + C0 carries out of 2^256
  exactly when r ≥ p, and then (r + C0) mod 2^256 = r − p < C0; a conditional select (CSEL, no
  branch) picks it.
- Mixed operands. The bound for fe_sub needs only its second operand canonical, and
  every use has that: dx = x_T − x_W and dy = y_T − y_W (both canonical), x3' = λ² − x_W
  (x_W canonical), x3'' = x3' − x_T (x_T canonical), x_W − x3 and y3' = λ·(x_W − x3) − y_W
  (second operands x3, y_W canonical after fe_canon).

**Pass 1** (forward, n ≥ 2 walks W_i = (x, y) and table points T_i = *tp[i], all
canonical, x(T_i) ≠ x(W_i)): dx_i = x(T_i) − x(W_i) (canonical), acc[i] = acc[i−1]·dx_i
(weak), acc[0] = dx_0. Because 0 < dx_i < p and p is prime, acc[n−1] ≢ 0.
**Cold path.** If dx_i = 0 (x(T_i) = x(W_i), which includes W_i = ±T_i; rate < 2^-114 in
total) the walk is flagged abandoned and dx_i is replaced by 1, both in the running product and in
its slot. Pass 2 reads dx_i from the slot, so it sees the same nonzero value as the product
chain, the sum it computes for that walk (a function of the wrong inverse) is discarded, and
no zero enters the shared inversion. The cold path is priced as its own abandonment action.
**Inversion**: Rinv = acc[n−1]^(p−2) by the libsecp256k1 addition chain: 255 squarings and
15 multiplications, 270 calls of two non-inlined functions fe_mul_f and fe_sqr_f (compiled by clang,
`secp_fe.s`; the bound of section 5.2 is taken from that file).
**Pass 2** (backward, i = n−1 .. 1; R starts as Rinv ≡ 1/(dx_0···dx_{n−1})):
inv = R·acc[i−1] ≡ 1/dx_i (one fe_mul), R ← R·dx_i (dx_i read from its slot; one fe_mul),
then dy = y(T_i) − y(W_i); λ = dy·inv; λ²; x3 = λ² − x(W_i) − x(T_i), canonicalised;
y3 = λ·(x(W_i) − x3) − y(W_i), canonicalised; W_i ← (x3, y3). The last walk i = 0 uses
R itself as its inverse. This is the standard affine addition λ = (y_T − y_W)/(x_T − x_W),
x3 = λ² − x_W − x_T, y3 = λ(x_W − x3) − y_W, so W_i ← W_i + T_i with canonical output.
Per walk: 5 multiplications (1 in pass 1, 4 in pass 2) and 1 squaring, i.e. the organizer's
6M, plus six fe_sub, two fe_canon and the memory traffic.

**The counted loops are generated assembly.** `kernels/gen_loops.py` emits the instruction sequences of
fe_mul, fe_sqr, fe_sub and fe_canon (the algorithms above, instruction by instruction) and schedules the
two loops as straight-line code with explicit register allocation from a pool of 23 registers
(x5 to x17 and x19 to x28; x0 to x3 hold the loop counter and the three array pointers, x4 holds C0;
x18, x29 and x30 are never written; the callee-saved registers are saved once per call). Values stay in
registers across operations: the only memory traffic per walk is the loads of the operands, the stores of the
results, one 32-byte stack slot for R (two stores, two loads per iteration), and re-loads of x(W), x(T)
and y(W) where holding them would need more than 23 registers. Pass 1 keeps the running product in four
registers across iterations (the final reduction of `mul` writes its result into the registers of the second operand,
which are dead by then) and stores [dx_i | acc_{i-1}] into the 64-byte slot of walk i; slot n receives acc_{n-1} for the
inversion. Pass 2 reads acc_{i-1} and dx_i from the slot. Pointer updates are free post-index forms
of the last load or store of each array (`ldp ..., [x3], #-32`), so a loop iteration ends with only
`subs x0, x0, #1; b.ne 1b`. The shipped `secp_loops.s` is the generator's `elf` output, reproduced byte for
byte by `kernels/build.sh`.

**Tests.** (1) `kernels/test_native.py`, on an AArch64 machine (the development machine was Apple silicon;
the generator's `macho` flavor differs only in symbol names and directives), runs the C file and the
generated assembly against Python big integers and affine secp256k1 arithmetic: 60,361 fe_mul_f and
26,019 fe_sqr_f inputs including edge values (0, 1, p−1, p, p+1, 2^256−1, single-limb patterns, values
within 2^40 of p and of 2^256), 317 inversions, 3,837 walks in batches from 2 to 1,000 walks with random
and negated table points through the assembly loops, 172,399 inputs through assembly wrappers of the four
emitters (`t_mul`, `t_sqr`, `t_sub`, `t_canon`: the same emitters, edge values included), and the assembly loops
against the C loops on batches up to 2,000 walks: 0 mismatches everywhere, outputs compared as canonical
integers. Re-run for this revision on Apple silicon with Apple clang 21.0.0 (`kernels/build.sh`):
`secp_loops.s` and `secp_fe.s` regenerated byte for byte, and the native test passed with 0 mismatches over
60,361 fe_mul_f and 26,019 fe_sqr_f inputs, 317 inversions, 3,837 walks, 172,399 emitter inputs, and no
walk on which the assembly loops and the C loops differ. (2) `experiment.py` (kernels mode) contains a small AArch64 interpreter (about 25 instruction forms)
that executes the shipped `secp_loops.s` itself (the ELF flavor, the file that is counted) on production
secp256k1 points: pass 1, the inversion by Python's `pow`, then pass 2, for batches of 2, 3, 5, 9 and 17
walks, comparing every output with affine arithmetic, and it records the executed-instruction counts; the
difference per extra walk is exactly 111 (pass 1) and 523 (pass 2). It also keeps the limb-for-limb Python
ports of the C kernels that assert the no-overflow claims above (s4 ≤ C0, w1 ≤ 1, no carry out of the last
limb) and a 24-walk batch round.

### 5.2 Counted instructions

Build (python3 for the loops, clang for the inversion kernels): `python3 gen_loops.py elf >
secp_loops.s`; `clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only
-ffreestanding -fno-builtin -fPIC -S secp_walk.c`, then `strip_asm.py` extracts `secp_fe.s` (fe_mul_f,
fe_sqr_f, fe_inv), keeping instructions and labels only. The compiler is Apple clang 21.0.0
(clang-2100.0.123.102) cross-compiling to AArch64 Linux with the COST_MODEL flags (clang in place of
GCC 14.2). The counts below are of the shipped files, which `experiment.py` parses.

| function | instructions | what it is |
| --- | ---: | --- |
| pass1 loop iteration | **111** | one fe_mul (91), one fe_sub (9), 5 loads (table pointer, x(T), x(W)), 4 stores (acc_{i-1}, dx_i), loop control (2) |
| pass2 loop iteration | **523** | 4 fe_mul (364), 1 fe_sqr (71), 5 fe_sub (45), 2 fe_canon (16), 19 loads, 6 stores, loop control (2) |
| fe_mul_f (whole function, compiled) | 115 | includes a 12-instruction frame save/restore and RET |
| fe_sqr_f (whole function, compiled) | 82 | includes a 2-instruction frame and RET |
| pass1 / pass2 whole function | 142 / 879 | prologue, first element, epilogue |

Both loops are single basic blocks: every conditional is a CSEL, there is no internal branch, and
`experiment.py` verifies that the only branch in each loop is the back edge. Every iteration therefore
executes exactly 111 and 523 instructions, whatever the data, and the interpreter's executed counts
confirm it. Opcode histograms: pass1 — ADCS 27, MUL 21, UMULH 21, ADDS 11, ADC 10, SBCS 5, LDP 4, STP 4, SUBS 3, CSEL 2, B.NE 1, LDR 1, SBC 1. pass2 — ADCS 132, MUL 99, UMULH 99, ADDS 54, ADC 47, SBCS 25, CSEL 18, LDP 18, SUBS 11, STP 6, EXTR 5, SBC 5, B.NE 1, LDR 1, LSL 1, LSR 1. There are no register spills: the
only stack traffic is the 32-byte slot holding R.

**Per-batch overhead.** The inversion makes 255 calls of fe_sqr_f and 15 of fe_mul_f:
at most 255·(82 + 1) + 15·(115 + 1) = 22,905 instructions, plus the non-call
instructions of fe_inv (397) and at most 4 per squaring for argument set-up and loop
control (1,020): **I ≤ 24,322**. With the prologue, first/last element and epilogue of
the two passes (879 − 523 + 142 − 111 = 387) and a 64-instruction allowance for the
call and the return, a batch of at least 2^15 walks has fixed work of at most 24,773 <
2^15 instructions, **under one instruction per walk**. Batches have 2^15..2^16 walks.

### 5.3 Per ordinary step

| work | instructions |
| --- | ---: |
| counted: pass1 iteration | 111 |
| counted: pass2 iteration | 523 |
| §1.1 F1: 4 loads of x, 5 of G, 4 MUL, 4 UMULH, 7 ADD | 24 |
| §1.1 range test: SUB, CMP, branch, S, NEG, CSEL, CMP, branch | 8 |
| §1.1 decision: differences, gap tests, argmin, fallback branch | 13 |
| sign, p − y with the C0 fold, compare, 4 CSEL, s | 17 |
| h, distinguished test, twist index, table address | 7 |
| arrival traffic: 4 loads of y, 1 store of the table pointer | 5 |
| fingerprint tests, history and checkpoint stores | 26 |
| zero-denominator test (pass 1, outside the counted 111) | 4 |
| length, cap test, meter additions, arrival loop | 13 |
| share of the inversion and of the call overhead (< 1) | 1 |
| **total** | **752** |

The step is charged **752 instructions = 0.752 units, its table exactly, with no margin**. The two
counted rows are the shipped assembly: each loop is one basic block, so an iteration executes exactly 111
and 523 instructions whatever the data, and `experiment.py` executes that file in its interpreter and
measures the same two numbers on batches of 2, 3, 5, 9 and 17 walks. The other rows are the sum of the
explicit sequences below, with no allowance on top (ab043cf7 re-derived each against 88e72314), and
`certificate.py` re-adds them to 118. The analytic sequences are:

- F1 (24): 4 loads of x and 5 of G_1..G_5 (9), 4 MUL and 4 UMULH (8), 7 ADD (7). Nine loads,
  not the twelve v8 allowed: F1 reads x's four limbs and the five G limbs its eight products use.
- Range test, S, wrap test, F2 (8): SUB, CMP, branch (3); S = F0 + F1 by ADDS (1);
  NEG (1); CSEL on the carry, CMP with 40, branch (3). F0, the top limb of the new x, costs no load:
  F1 already loaded it.
- Decision (13): CMP, 2 CSEL, CSET (4); CMP, CSEL (2); CMP and 3 CSEL (4); difference,
  CMP with 40, branch (3). A gap above 40 to the second smallest gives a gap above 40
  to both others.
- Sign, p − y, K (17): p − y = 2^256 − C0 − y by the four-limb negation from zero (SUBS, 3 SBCS,
  4) and the C0 fold (SUBS of C0, SBCS, SBCS, SBC, 4); the comparison y against p − y by SUBS, 3 SBCS
  (4), whose flags feed the four CSEL that select the limbs of K (4); CSET of s (1). v8 charged four
  instructions for p − y, which yields 2^256 − y and not p − y.
- h, distinguished test, t, address (7): AND (1); CBZ on limb 1 of K (1); t = 3s + j by
  two ADD (2); address = base + 64·(6h + t) by three ADD with shifts (3).
- Arrival traffic (5): 4 loads of the y limbs of the new point (the x limbs are the 4 loads of
  F1); 1 store of the new table pointer into the pointer array the loops consume. The table-pointer
  load and the stores of the new x and y are inside pass 2's 19 loads and 6 stores; v8's "3 header
  words" had no derivation.
- Fingerprints (26): pending-escape flag, load and CBNZ (2); two-slot ring indexed by
  length mod 2, slot address by AND and ADD (2); compare fp(K) with K_{−2}: 3 loads, CMP
  and 2 CCMP, branch (7); overwrite the slot with fp(K), 3 stores (3); checkpoint compare:
  3 loads, CMP and 2 CCMP, branch (7); checkpoint write when length mod 64 = 0: TST,
  branch (2), 3 stores (3).
- Zero denominator (4): in pass 1, ORR of the four limbs of dx in pairs and together, then CBZ to
  the cold path of section 5.1 (dummy table point, trail abandoned), priced as a separate action
  (rate < 2^-114 in total). The counted pass-1 iteration of 111 instructions contains no such test, so
  these four instructions are charged in this row and nowhere else. Because dx is canonical, dx = 0
  exactly when x(T_i) = x(W_i), which also covers W_i = −T_i.
- Length, cap, meter, loop (13): the 2-limb length (L_max > 2^64): 2 loads, ADDS + ADC,
  2 stores (6); cap test CMP of the high limb with 20, branch (2); meter and operation
  count: SUBS on the processor's block remainder, ADD on its operation count, branch to
  the block refill (3); the loop counter SUBS and branch of the arrival pass (2). The
  loops of pass 1 and pass 2 carry their own counters inside the counted 111 and 523.
- Share (1): fixed work per batch under 2^15 instructions for at least 2^15 walks, section 5.2, so
  the share is below one instruction and one is the ceiling of the bound.

**Other actions** keep the allowances of 3682292 and are priced as before. They are used
only outside the hot loop and are loose against the counted kernels:

- A = 32 (load 8, ADDS + 3 ADCS, carry, SUBS + 4 SBCS, 4 CSEL, store 4, ≤ 6 overhead).
- M ≤ 227. The compiled fe_mul_f is 115 instructions including its frame (91 in the assembly loops); the older
  allowance is retained for the rare paths.
- **Doubling** (priced 0.956 = 956 instructions, its derived bound; v8 priced it 1.1 and v5 and earlier 2.0).
  In a batch round a walk that must double replaces its addition by λ = 3x²/(2y) (y ≠ 0: a point with y = 0 has
  order 2, and the group has prime order), sharing the same product chain and inversion: dx = 2y (one weak
  addition, 9 instructions: ADDS, 3 ADCS, CSEL of C0 on the carry, ADDS, 2 ADCS, ADC; y is canonical, so the carry fold
  cannot carry again), the chain multiplication, inv = R·acc, R ← R·dx, x² (squaring), 3x² = v + v + v (two weak
  additions, each with a second carry fold because v is only weak: 14 instructions each), λ = 3x²·inv, λ²,
  x3 = λ² − x − x (two subtractions) and canonicalisation, x − x3 (subtraction), λ·(x − x3) (multiplication),
  y3 = · − y (subtraction) and canonicalisation. That is 5 multiplications (455), 2 squarings (142), the addition
  for dx (9), two weak additions (28), 4 subtractions (36) and 2 canonicalisations (16): 686 instructions, plus
  at most 72 for loads, stores, pointer updates and loop control (the addition loops spend 19 + 6 + 2 and
  5 + 4 + 2 on them for the same operands): **at most 758**. `experiment.py` recomputes this from the operation
  costs 91, 71, 9, 8 and checks that the same costs reproduce the counted 111 and 523 of the addition loops
  (111 = 91 + 9 + 5 + 4 + 2; 523 = 4·91 + 71 + 5·9 + 2·8 + 19 + 6 + 2). The doubling is also one class arrival
  and carries the 118 bookkeeping rows of section 5.3; moving the walk to and from the doubling list costs at most
  16 loads, 16 stores and 18 other instructions (50), and the escape logic (clearing the history, the order test of
  fp(K) against K_{-1}) at most 30: in all at most 758 + 118 + 50 + 30 = 956, which is the price. The order test is
  charged once per doubling although the pseudocode of section 2 performs it on the arrival that detects the cycle,
  which is the preceding arrival when the escape is deferred; the two rates agree because every detection produces
  exactly one doubling.

**Prices** (units):

| action | price | basis |
|---|---|---|
| ordinary step | 0.752 | 752 instructions: 634 counted + 118 analytic rows = 752, charged exactly (section 5.3) |
| traversal step | 0.776 | ≤ 24 more instructions to track and resume the minimum-fp member: ≤ 776 |
| doubling | 0.956 | ≤ 956 instructions: 758 arithmetic + 118 rows + 50 list moves + 30 escape logic (section 5.3) |
| §1.1 fallback | 0.5 | one M, two A, compares: < 400 instructions |
| random 64-bit word | 0.01 | |
| sampling attempt | 0.05 | |
| trail start [a]P + [b]Q | 2^13 | two Jacobian scalar multiplications at ≤ 512 operations × (16M + 16A) each, an addition and an inversion: < 4.8·10^6 instructions |
| report | 40.4 | 40 for route, message, home-bucket multiply-shift, tag/key tests, record write and insert into a bucket found within the eight probes (< 2^15 instructions) + eight probes at 0.05 (section 3) |
| store probe | 0.05 | 39 itemised instructions, charged 50 (section 3) |
| overflow insert and lookup | worst case, hard extras | 6B + 35 and 16B per report, charged in the overflow row of the budget (section 3) |
| store initialisation | 2·D_cap buckets × 16 instructions = 2^58.865, once | |
| table | 2048 starts plus 5 φ-twists each: < 2^24.1, plus sampling | |
| phase-2 or re-walk operation | 600 | unbatched: an inversion (< 2^17) plus 7M, coefficient updates mod r and control, far below 600,000 instructions |
| small-batch inversion | 600 | one inversion, ≤ 24,322 instructions = 25 units; 600 is a loose upper bound |
| solving | 600 | |
| verification | 2^13 | |
| meter block refill | 0.03 | 30 instructions: atomic 192-bit add, remainder store, phase-1 limit load and the two limit tests (section 3) |

**Budget.**

- **Phase-1 operations:** 0.752·(S_int + 2^61) < 2^126.31721.
- **Hard extras:**

  | item | units |
  |---|---|
  | meter remainders | < 2^61 |
  | meter block refills: 2^86.3173 × 0.03 | < 2^81.26 |
  | phase 2: 2^20·(L_max + 65)·600 | 2^97.551 |
  | re-walks: 8·2·((L_max + 65)·600 + 2^13) | 2^81.6 |
  | reports: D_cap·40.4 | 2^69.17 |
  | overflow insert and lookup, worst case: D_cap·(22B + 35) | 2^107.154 |
  | store initialisation | 2^58.865 |
  | table and verification | < 2^25 |
  | node set-up | 2^26 |
  | **total** | **< 2^107.16** |

- **Random charges, expected (for reference) and allowed (section 4):**

  | item | expected (units) | allowance (units) | control |
  |---|---|---|---|
  | doubling surcharges: 0.204·D | 0.204·μ_D < 2^110.867 | 0.204·(⌈μ_D⌉ + t_D) = 2^111.2832 | Chebyshev, H2 |
  | traversal surcharges: ≤ 1.536·V | 1.536·μ_V < 2^107.358 | 1.536·(⌈μ_V⌉ + t_V) = 2^109.3921 | Chebyshev, H2 |
  | starts: 2^13·(S_int/2^64 + 2^20) | < 2^75.7284 | ×2^30 = 2^105.7284 | Markov |
  | words, attempts, probes, small batches, fallbacks | < 2^71 | ×2^30 = 2^101 | Markov |
  | **total** | **< 2^110.98838** | **2^111.65243** | 4·2^-30 |

- **Early meter halt:** only if one of the four charges exceeds its allowance, probability ≤ 4·2^-30.
- **W as the sum of its rows:** phase 1, the hard extras and the allowances add to
  105995162764622566518652505820930400478098 milli-units = 2^126.317266, whose four-decimal ceiling is
  the claim, 126.3173.

`experiment.py` with `"kernels": true` recomputes this ledger to 60 digits (S_ops, phase 1, hard rows,
expected charges, allowances, success floor and W_ledger = 2^126.317266).

**Exact integer ledger and the four-decimal claim.** The meter counts milli-units, 1,000 to a
unit, and every price is a whole number of them, so every charge and the total are integers in
milli-units. S_ops is 1.02·√(2N' ln(1/(1 − 0.3901))) with the logarithm taken between rational
bounds and the square root rounded up, so S_int is exact, and so are the report capacity, the bucket
count and every row:

    S_int = ⌈S_ops⌉   = 140945347563975464246948179984935587492
    D_cap = ⌈2^63.83⌉ = 16396253461269537204
    B     = ⌈2·D_cap/2^16⌉ = 500373945961595

Phase 1 charges 752 milli-units per operation for S_int + 2^61 operations, which covers the < 2^60.41
block overshoot: 752·(S_int + 2^61) = 105990901368109549115439025291600259645888 milli-units.
The hard extras are the table's rows, each taken at the integer its bound gives:

| hard row | milli-units |
|---|---:|
| meter remainders 2^61 | 2305843009213693952000 |
| meter block refills: ⌈2^126.3173/2^40⌉ + 2^20 = 96404308701361100171041323 refills × 30 | 2892129261040833005131239690 |
| phase 2: 2^20·(L_max + 65)·600 | 232113757366008801584480256000000 |
| re-walks: 16·((L_max + 65)·600 + 2^13) | 3541774862152233911027072000 |
| reports: D_cap·40,400 | 662408639835289303041600 |
| overflow insert and lookup: D_cap·(22B + 35) | 180493676954842338364814213402770500 |
| store initialisation: 2·D_cap·16 | 524680110760625190528 |
| table and verification, < 2^25 | 33554432000 |
| node set-up, 2^26 | 67108864000 |
| **total** | **180725797146777709522420974102562318** |

That total is 2^107.15550 units, inside the table's < 2^107.16; its largest row, the overflow
worst case, is 2^107.15365. The expected random charges total 2^110.98838 units (for reference only),
with ⌈μ_D⌉ = 11605006972542714042418004049532163 and ⌈μ_V⌉ = 135350894912927938485939642319581. The
allowances are

    doublings  204·(⌈μ_D⌉ + t_D),  t_D = 3881257568327832090394473406362184   3159197966337591411093745401002446788
    traversals 1536·(⌈μ_V⌉ + t_V), t_V = 419160478354136256443848636853317     851729469338210603412154796809571328
    starts     2^30·E_starts                                                    67207978994367030649214751815499776
    small      2^30·2^71·1000                                                    2535301200456458802993406410752000
    total                                                                       4080670715870625503958108356038269892 = 2^111.65243 units

where t_D and t_V are the least integers with t² ≥ 2^80·μ·2^30. The three parts sum to

    phase 1        105990901368109549115439025291600259645888
    hard extras       180725797146777709522420974102562318
    allowances       4080670715870625503958108356038269892
    W_ledger       105995162764622566518652505820930400478098

W_ledger/1000 = 2^126.3172660…, whose ceiling at four decimals is 126.3173 = 1263173/10000. It sits
0.000034 bits below 2^126.3173, so the meter's own cap never binds before the ledger's total. The
certificate is two exact integer inequalities, and both hold:

    W_ledger^10000 < 1000^10000·2^1263173    so 126.3173 holds
    W_ledger^10000 > 1000^10000·2^1263172    so 126.3172 fails

**Certificate.** `certificate.py` (standard library, exact integers and rationals, no floats) checks the
instruction identities and the 118 rows, D_cap and B, S_int² ≥ the certified upper bound on S_ops² (ln by
the atanh series with its remainder; fractional powers of two by integer roots m^den ≤ 2^sh < (m+1)^den),
the nine hard rows, the rate bound, the two Chebyshev tails (Var/t² ≤ 2^-30 under H2), the two Markov
ratios (exactly 2^-30), W_ledger = phase 1 + hard + allowances, the two exponent inequalities above, the
counter widths and the success chain 0.3901 − 2e-7 − 4·2^-30 > 0.39009. It prints `CERTIFICATE: PASS`
and exits 0, in about ten seconds; any failed comparison makes it exit 1.

**Preprocessing:** table < 2^24.1, node set-up 2^26, bucket initialisation 2^58.865 and the
first 2^20 starts 2^33. The expected total is 2^58.864216, and the cap is 2^59, so the gap is
> 2^55.52. The random sampling part is expected < 2^25, so the cap halts a run with probability
< 2^25/2^55.5 = 2^-30.5 < 2^-30 by Markov (the exact ratio is 2^-30.52). All of it is inside W.

**W ≤ 2^126.3173 on every run**, and with the random charges inside their allowances,
W ≤ W_ledger = 2^126.317266.

**Memory:** each of the 2^16 nodes holds B = 2^48.83 buckets of 104 bytes plus its B 8-byte
overflow addresses, 112 bytes per bucket = 2^55.64 bytes, and the whole store is 2^64.83 buckets
× 112 bytes = 2^71.64 bytes. Adding the table (1.5 MB), the 2^20 × 512 bytes of walk states and
the batch scratch (40 bytes per walk, under 2^26 bytes) gives < 2^72 bytes, with each node
< 2^56 bytes. Every live array is counted: the bucket array, the overflow list, the walk states,
the per-walk pointer array and the table.

## 6. Experiments and limitations

`experiment.py` uses the standard library only.

**Part A** (production secp256k1; all true): curve, order and cofactor; the
β/λ pair; the order of 2 mod r; orbit size 6, unit action, |y| invariance and a
unique canonical member; twisted step and doubling as class maps; the table;
coefficient tracking; a planted-x recovery from two production trails forced
to merge; and the §1.1 rule against the exact rule.

**Part B** runs the section 2 per-trail algorithm sequentially on configured
toy curves y² = x³ + 7 with q = 33556867 (cofactor 1) and q = 536873251
(cofactor 3). Setup now checks [r]G = O without reducing r first; it does not
independently count each curve's full group order. The scaled toy pipeline
adds a κ·S0 budget (1.05 in the offline rows below; the toy pipeline of this package uses the constants KAPPA and P0 of experiment.py), phase 2, report
capacity, the 8-match cap, and routing to 16 nodes of capacity ⌈2·D_cap/16⌉
with overflow treated as failure. [a]G + [b]T is evaluated as [a + bx]G for
speed. The toy node model bounds capacity only: the probe depth of 8 and the
overflow list of section 3 are analytic and change no capacity event, so the
toy's node-overflow counts remain the right scaled evidence for it.

Offline results (seeds 1501–1512; every recovered x verified; at the earlier parameters, S_ops = 1.05·√(2N' ln(1/0.605)) for each toy; T_c counts all arrivals):

| mode / classes | toy, R | trials | mean T_c / √N (model) | P[T_c ≤ S_ops] or success | 2-cycles per step (model) |
| --- | --- | ---: | --- | ---: | --- |
| T_c, 6 | q25, 256 | 6000 | 1.2689 (1.2558, +1.0%) | 0.418 | 0.00063 |
| T_c, 6 | q25, 2048 | 6000 | 1.2696 (1.2536, +1.3%) | 0.4185 | 8.23e-05 |
| T_c, 6 | q29, 256 | 1800 | 1.2441 (1.2558, -0.9%) | 0.4411 | 0.000632 |
| T_c, 6 | q29, 2048 | 1800 | 1.2699 (1.2536, +1.3%) | 0.4128 | 8.28e-05 |
| T_c, 1 (control, no automorphism) | q25, 256 | 3000 | 1.2435 (1.2558, -1.0%) | 0.429 | 0 |
| T_c, 2 (control, negation only) | q25, 256 | 3000 | 1.2741 (1.2558, +1.5%) | 0.4117 | 0.0019 |
| pipeline success with 16 nodes, 6 | q25, 2048 | 12000 | – | 0.4229 (node overflows 0, max load/cap 0.72) | 7.86e-05 |
| pipeline success with 16 nodes, 6 | q29, 2048 | 2400 | – | 0.4238 (node overflows 0, max load/cap 0.52) | 7.72e-05 |

Part A offline (seed 1513): §1.1 rule on 2000000 random inputs gave 0 fallbacks and 0 mismatches. On 1000000 adversarial inputs it gave 344641 fallbacks and 0 mismatches. All other part A checks are true.

**Census of the slowdown.** Six runs of `experiment.py` (toy on, pipeline off, classes 6,
R = 2048, dp_bits = 5, production_checks 0, seeds 101, 102, 103, 420 s each, on curves q25 and q29):

| curve | trials | mean T_c/√N | model | slowdown | P[T_c ≤ S0(0.395)] (model 0.395) |
| --- | ---: | ---: | ---: | ---: | ---: |
| q25 | 96,675 | 1.2677 | 1.2536 | +1.13% ± 0.17% | 0.3878 |
| q29 | 65,641 | 1.2563 | 1.2536 | +0.22% ± 0.20% | 0.3935 |

The slowdown falls by a factor of five from q25 to q29 (N is 16 times larger), as expected for
an O(ℓ/N) effect of short toy trails meeting their own past (about 6·10^-6 per step at q25,
see above). The doubling rate per step is 8.4·10^-5 (q25) and 8.0·10^-5 (q29) against
1/(6R) = 8.1·10^-5. At production size ℓ/N < 2^-183, and the walk's own overheads (doublings and
traversals, about 10^-4 of the operations) are counted in T_c. This extrapolation is a heuristic;
H1 asks for a production slowdown below κ − 1 = 2%.

**Long-trail trap test** (q29, walks without distinguished points, 150 trails
of 2,048 operations each, against the additive-escape control):

| R | escape | escapes per trail | trails with non-fruitless returns |
| ---: | --- | ---: | ---: |
| 4 | doubling | 72 | 0 / 150 |
| 8 | doubling | 42 | 6 / 150 |
| 4 | additive (control) | 223 | 110 / 150 |
| 8 | additive (control) | 64 | 8 / 150 |

Ordinary self-collisions predict non-fruitless returns in about 1 − exp(−L²/(2N)) = 6.8% of trails (L²/(2N) = 0.0703). At R = 4 it traps 110/150, while the doubling escape gives 0/150 and 6/150.

Doubling escapes behave like ordinary self-collisions; the additive control
traps trails structurally.

**Doubling rate against the toys.** The toy 2-cycle rates at R = 2048 (8.23·10^-5 and
8.28·10^-5 per step in the T_c rows, 7.86·10^-5 and 7.72·10^-5 in the pipeline rows) agree with
1/(6R) = 8.14·10^-5 within about two standard errors (10^3 to 2·10^3 events per row). The toys also
show checkpoint escapes at about 5·10^-6 per step, above the structural rate 2^-19.99: a short toy
trail of length ℓ meets its own past with probability about ℓ/N per step, about
32/(5.6·10^6) ≈ 6·10^-6 at q25 with dp_bits = 5, while in production ℓ/N < 2^70/2^253. In the q25,
R = 2048 row (seed 1502) the total doubling rate, both effects included, is 8.72·10^-5 per step, 6%
above the production bound 2^-13.5681 = 8.23·10^-5; the toy-only self-meeting rate accounts for
that excess. This small toy sample does not certify the production instruction budget.

**Provenance of the rows.** `experiment.py` is that of 3682292 plus the optional `kernels` mode of
section 5, as revised by ab043cf7. This revision changes P0 (0.3953 → 0.3901, which also sets the toy's
budget κ·S0) and the `kernels` ledger (sections 4 and 5), and adds `doubling_count_dispersion`,
Σ_i(d_i − ρ̂ s_i)²/Σ_i ρ̂ s_i over trials, with ρ̂ = Σd_i/Σs_i. It is about 1 for independent rare
events, and H2 needs only ≤ 2^80. `certificate.py` is rewritten. The walk, the toy, part A, the kernel ports
and the interpreter are unchanged. The offline rows
above are the earlier packages' results at their own parameters, not rerun for this revision; Python's
seeded generator makes each reproducible.

**Initial input.** The declared configuration (q25, R = 2048, dp_bits = 5, seed 1, pipeline, kernels,
11 s, index-rule sample capped at 4,000 random and 2,000 near-tie inputs so that the toy gets time;
ab043cf7 used 6 s and 30,000/10,000) is time-limited, so what completes depends on the machine. The kernels check and part A's
identities are clock-independent: every run reports the interpreter over the shipped assembly with 0
mismatches, 111 and 523 instructions per walk, the ledger above and the identities all true. The
index-rule sample (up to 30,000 random and 10,000 near-tie inputs) and then the toy pipeline use the
remaining time: a containerised run stopped the sample after 16,202 random inputs with 0 mismatches
and ran no toy trial, while a native AArch64 run finished the sample with 0 mismatches. The toy runs one
seeded stream, so a slower machine sees a prefix of a faster machine's trials; native runs completed 262
and 316 trials, with pipeline success 0.328 and 0.351 and P[T_c ≤ S0] = 0.305 and 0.329, 3.0 and 2.4 binomial
standard errors under the model's 0.3953. The same stream continued to the configuration's 2,000-trial
cap (a 110 s budget, done after 74 s) gives success 0.398, P[T_c ≤ S0] = 0.374 and
P[T_c ≤ S_ops] = 0.3905 (standard error 0.011) with 0 node overflows. A run of the toy pipeline alone
(seed 1, no kernels or part A, 12 s) completed 258 trials with success 0.4186. Those runs are ab043cf7's,
at p0 = 0.3953. In this package P0 = 0.3901, so the toy's S0 and S1 are 0.85% smaller and the pure model
gives 0.3901 at S0 and 0.4022 at S1 (0.3953 and 0.4075 before). Applying that shrink to the same seed-1
stream would move its P[T_c ≤ S_ops] from 0.3905 to about 0.385, within one standard error (0.011) of
the floor, while the q25 census (slowdown 1.13% against κ − 1 = 2%) implies about 0.395. The hosted
initial run reports the new values, and `doubling_count_dispersion` if toy trials complete within its
budget; if none complete, H2 rests on the analytic argument of section 4. The author did not run `experiment.py` for this
revision: Docker was not available, and candidate code is not executed outside the verifier's
container. The new ledger block was executed separately: it gives the certificate's S_int and, to its
60-digit Decimal precision, the same log2 phase 1 and log2 W_ledger (2^126.317266). None of these runs is the H1 or H2 evidence, and none measures production
instruction counts. The H1 floor rests on the disclosed heuristic and the censuses above.

**Limitations.** The instruction counts are those of the shipped generated assembly for the two
batch loops (executed natively and by the interpreter, but not a Callgrind measurement of GCC output and
not a measurement of a full attack), clang's output for the inversion kernels, and analytic rows for the
arrival bookkeeping. The 118 analytic rows rest on this package's explicit sequences (section 5.3 and the
row table of ab043cf7), not on generated code; eleven instructions of v8's 129 were allowance rather than
sequence, and the class-key row of v8 under-counted a subtraction by four. The interpreter implements the ~25 instruction forms the file uses; a mistake in it could
not hide a wrong result, because the outputs are compared with big-integer arithmetic, but the executed
counts are only as good as its handling of the branch. H2 is a second-moment heuristic about the
doubling and traversal counts (section 4). It is used only to bound the meter's early-halt
probability, it tolerates correlations 2^80 times stronger than independent events, and it is checked at
toy scale only by the dispersion statistic. H1 is a heuristic: a random mapping with the Brent–Pollard
correction; van Oorschot–Wiener; automorphism classes as in Wiener–Zuccherato,
GLV and Duursma–Gaudry–Morain; and fruitless-cycle escapes as in
Bos–Kleinjung–Lenstra and Bernstein–Lange–Schwabe. It is supported only at toy
scale. Instruction counts are analytic. This is a generic attack.

# Six-class parallel Pollard rho on secp256k1: the step priced from counted AArch64 kernels

**Provenance.** This is the public submission 88e72314 (six-class-rho-secp-counted-v8, 126.367) with its
accounting corrected where a price exceeded its own derivation, and with a fully specified report store, an
exact success chain, stated counter widths and a charged row for the meter's block refills. That line descends
from 061e00ca, 30f6f34e, 98e7729c, b52f8806, a9f1a6fe, 3682292, 758d8df6 and 189ebc42, all credited. No prior
acceptance or score is a premise.

**Claim (superseded).** That package claimed time_log2 = 126.338 under the H1 heuristic of section 4. This
package's claim, ledger and certificate are in proof.md section 5.3.

## 0. Abridged copy

This is the public package's proof (126.338 line) as carried by the accepted 02feb284 package, kept for the
construction it specifies: algebra (section 1), walk (section 2), phases, store and meter (section 3), success
probability (section 4) and the action prices of section 5.3. Its change history, superseded kernels and loops,
per-step table, row sequences, ledger, certificate and run descriptions are omitted. **Every kernel, row, price
and ledger number of this file is superseded by proof.md** (sections 2 to 5). Unchanged and still used: sections
1 to 4, the inversion kernels (`kernels/secp_fe.s`) and the action prices other than the step, doubling and
traversal. In section 1.1, F1, F0 and their error bounds are used as stated; the range, wrap and gap tests that
follow them are replaced by the rule of proof.md section 3.1.

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
  true bits 320..383 lie in [F1, F1 + 14].
- **v0.** F0 = top limb of x. Since 2^256/p = 1 + 2^-224·(…),
  0 ≤ 2^64·v0/p − F0 < 2.
- **The rest of the rule** (v2 from S = F0 + F1, the gap tests and the fallback) is replaced by proof.md
  section 3.1, which also states the Part A run. Its premise stays: βx mod p is equidistributed, so the
  fallback windows are rare; they affect only the meter's early-halt probability, never correctness.

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
  coins, stored in full (redraw if S = O). The order of tests at each class
  arrival and the merging lemma are those of proof.md section 5.1 (checkpoint
  only), which replace this section's period-64 version with its
  2-cycle test. ESC is the doubling escape P ← [2]P (a batched affine
  doubling). Doubling commutes with every unit, so it maps class C to
  class(2c(C)), which is a class map. fp(K) is limbs 0, 2 and 3 of K (192
  bits), compared as one integer (limb 3 most significant). Every equality
  test and every minimum uses this single fixed total order on classes.
  `length` counts every operation. If an addition has equal x-coordinates
  (P = ±Tab), the trail is abandoned (< 2^-114 in total).
- **No structural cycles through escapes.** A trail can return to a class c
  after m ≥ 1 doublings only if u·2^m·c + (a sum of unit and power-of-2
  multiples of the table points) = c for some unit u. The multiplicative order
  of 2 mod r is 2^250.0 (part A, from the full factorisation of r − 1). Hence
  1 − u·2^m ≠ 0 mod r for every m < 2^247, so the relation cannot hold through
  cancellation of table terms. (An additive escape multiplier would cancel
  identically through negation pairs.) It pins c to one specific point. This
  excludes identities, not coincidences. That such returns are no more
  frequent than ordinary self-collisions is part of H1, and the long-trail
  test checks it (proof.md section 5.2).
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
  verification. If a charge would exceed the cap W (proof.md), the attack halts
  with failure. So the work is at most W on every run; no loop needs its own limit;
  rejection-sampling loops in particular are bounded by the meter.
  - Preprocessing (all charges before the first walk step) has its own cap of
    2^59.0 units. If a preprocessing charge would exceed it, the attack halts
    with failure.
  - **Counter widths.** The meter counts milli-units in a 192-bit counter (three
    64-bit words, ADDS/ADCS/ADC chains). The phase-1 operation counter (steps,
    doublings and traversal steps; its stop test compares the 128-bit pair with
    S_ops < 2^127, so it can fire), the start counter, the preprocessing account
    (a separate pair, compared with the cap by SUBS/SBCS/CSEL, so the cap and the
    meter never alias) and each block-refill counter are 128-bit pairs.
    `certificate.py` checks each counter's largest value against its width at
    this package's cap (proof.md section 5.3).
  - Each of ≤ 2^20 processors reserves blocks of 2^40 units by atomic update on
    the 192-bit meter. The refill itself (an atomic 192-bit add of the block, a
    store of the new remainder, a load of the phase-1 limit, the operation count
    and the two limit tests) is charged 30 instructions = 0.03 units per
    refill, for at most W/2^40 + 2^20 refills, a hard row of proof.md section
    5.3. Unused remainders count as spent: < 2^61 units, also a hard row.
- **Phase 1.** Runs until a global operation counter (steps, doublings and
  traversal steps of phase-1 trails) reaches S_ops. Every reserved block is
  completed, so at least S_ops and fewer than S_ops + 2^60.41 operations are
  performed.
- **Phase 2.** Finishes the open trails and starts no new ones.
- **Caps.** At most 8 processed matches; a ninth is a failure. Report capacity
  is D_cap = 2^63.83 records, taken as the integer
  ⌈2^63.83⌉ = 16,396,253,461,269,537,204 in every product below. These bound
  memory and re-walk work.
- **Batches.** Superseded: proof.md section 4 runs all walks as one batch.
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
  - **Probe itemisation** (39 instructions, charged 50 = 0.05 units): j = (home + i) mod B by ADD, CMP, CSEL
    (3); 13·j by LSL #3, LSL #2, ADD, ADD, the loads using [base, X, LSL #3] (4); flags LDR, CBZ (2); tag LDR,
    CMP, B.NE (3); key 4 LDR (4) and 4 CMP/CCMP, B.EQ (5); record head LDR, visit stamp LDR, ORR, STR (4); next
    probe ADD, CMP, B (3); prefetch of the next bucket, 2 LSL, ADD, PRFM (4); two branch-shadow slots and loop
    alignment (3); register reloads of the tag and key base (4): 3 + 4 + 2 + 3 + 4 + 5 + 4 + 3 + 4 + 3 + 4 = 39.
    The 11-instruction margin covers the two conditional branches and the loop's register moves; every charge
    here is an upper bound.
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
with κ = 1.02 and p0 = 0.3953. H1 also includes three further claims:

- the partition, DP, node-routing and bucket fields of |y| are pseudo-random
  functions of the class: bits 0..10, 64..127, 128..143 and 144..193, four
  disjoint fields;
- excluding the detected fruitless cycles of length ≤ 64, a trail makes a
  non-fruitless return with probability ≤ L_max²/(2N);
- the canonical twist index of a new class is uniform and independent of its
  partition index.

The cycle rates at this package's R = 2^24 are in proof.md section 5.1.

The public package's values (κ = 1.02, p0 = 0.3953) and its censuses are superseded by proof.md sections
5.2 and 5.3.

**Phase 1 covers the budget.** Phase 1 performs ≥ S_ops operations, and each is
one class arrival. So if T_c ≤ S_ops, the first inter-trail collision happens
in phase 1.

**Doubling and traversal surcharges** are random; their rates and charges are in proof.md sections 5.1 and
5.3.

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
| meter halting early, by Markov | proof.md section 5.3 |
| preprocessing cap | < 2^-30 |

The meter row and the success floor are superseded by proof.md section 5.3, which also sums the other rows.

## 5. Cost and ledger (historical; superseded by proof.md sections 3-5)

COST_MODEL.md prices computation at C = 1,000 scalar AArch64 instructions per unit. Ordinary steps charge H = 0 reference operations; executed assembly instructions are counted and divided by C. The field kernels and loops are in proof.md section 2; inversion is 24,322 instructions per batch (`kernels/secp_fe.s`). Other actions keep historical allowances: A = 32, M <= 227. Batched doubling replaces addition by lambda = 3x^2/(2y).

**Prices** (units):

| action | price | basis |
|---|---|---|
| step, traversal, doubling | proof.md section 4 | |
| §1.1 fallback | 0.5 | one M, two A, compares: < 400 instructions |
| random 64-bit word | 0.01 | |
| sampling attempt | 0.05 | |
| trail start [a]P + [b]Q | 2^13 | two Jacobian scalar multiplications, addition, inversion: < 4.8e6 instructions |
| report | 40.4 | 40 for route, message, hash, tag/key tests, insert + 8 probes at 0.05 |
| store probe | 0.05 | 39 itemised instructions, charged 50 |
| overflow insert and lookup | worst case, hard extras | 6B + 35 and 16B per report, charged in overflow row |
| store initialisation | 2·D_cap buckets × 16 inst = 2^58.865, once | |
| table | 2048 starts plus 5 phi-twists: < 2^24.1, plus sampling | |
| phase-2 or re-walk operation | 600 | unbatched: inversion plus 7M, updates mod r |
| small-batch inversion | 600 | one inversion, <= 24,322 instructions |
| solving | 600 | |
| verification | 2^13 | |
| meter block refill | 0.03 | 30 instructions |

**Budget & Preprocessing:** Hard extras and ledger recomputed by `certificate.py` at cap. Preprocessing expected 2^58.864216 units, cap 2^59. Memory: 2^16 nodes hold B = 2^48.83 buckets, whole store 2^71.64 bytes (< 2^72 bytes).

## 6. Experiments and limitations

Abridged; superseded by proof.md sections 2.6, 5.2 and 6.


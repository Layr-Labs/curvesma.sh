# Six-class parallel Pollard rho on secp256k1: the step priced from counted AArch64 kernels

**Provenance.** This package is the public submission 88e72314 (six-class-rho-secp-counted-v8, the
126.367 bound) with its accounting corrected wherever a price exceeded its own derivation — the doubling
price, the analytic bookkeeping rows and the ordinary-step margin — and with four gaps closed: a fully
specified report store with a fixed probe depth and a charged overflow list
(section 3), an exact success chain (section 4), stated counter widths (section 3) and a charged row for
the meter's block refills (section 5). Its walk, H1 constants, Markov ledger, merging lemma, doubling
derivation, generated loops and every experiment are kept. That line descends from 061e00ca (v7,
126.403), 30f6f34e (v6, 126.423), 98e7729c (v5, 126.551), b52f8806 (v2), a9f1a6fe (v1), 3682292
(six-class-rho-secp-147-v1, 127.00), 758d8df6 and 189ebc42; all are credited, and 88e72314, the
package revised here, is the declared dependency. The complete argument is restated here; no prior acceptance or score is used as a premise.
Section 0 itemises every number that moves and the derivation that fixes it, section 0.1 re-derives the
bookkeeping rows and section 5 the exact integer ledger, and `certificate.py` proves the resulting bound
in exact integer arithmetic both ways. The shipped assembly of the two counted loops is unchanged:
111 + 523 instructions per walk.

**Claim.** time_log2 = 126.338 = 63169/500. The classical cold attack recovers x from a uniform secp256k1
point Q=[x]P with probability at least 0.3911 under the explicit H1 heuristic in section 4. On every run
total work is at most 2^126.338 target-group-operation equivalents because the global meter stops before
that bound, and the exact integer ledger totals
W_ledger = 107529267059811478367772329467123690012431 milli-units = 2^126.337997 units (section 5),
every row of which is charged. The unit is 1,000 scalar AArch64 instructions (C = 1,000, fixed by the
organizer); no reference operation is charged for ordinary steps, see section 5.0. Peak memory is below
2^72 bytes (2^71.64 for the store alone), preprocessing is capped at 2^59.0 units and included, and
there is no nonuniform advice. This is an analytical, impractical generic attack, not a demonstrated key
recovery.

## 0. Abridged copy

This is the public package's proof (six-class-rho-secp-counted-v8 line, 126.338) as carried by the
accepted 02feb284 package, kept for the construction it specifies: algebra (section 1), walk (section 2),
phases, store and meter (section 3), success probability (section 4) and the action prices of section 5.3.
To stay within the archive limit, sections 0 and 0.1 (change history from 88e72314 and the re-derivation of
its 118 rows) and sections 5.1 and 5.2 (the superseded 91/71-instruction kernels and 111/523 loops) are
omitted, as are section 5.3's exact ledger and certificate of the 126.338 claim and section 6's notes on that
package's rows, initial input and limitations. **Every kernel, row, price and ledger number of this file is superseded by proof.md**: kernels and
loops (proof.md section 2), rows (section 3), prices (section 4) and the exact budget (section 5). Unchanged
and still used: sections 1 to 4, the inversion kernels (`kernels/secp_fe.s`) and the action prices other
than the step, doubling and traversal.

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
  verification. If a charge would exceed W = 2^126.338 units, the attack halts
  with failure. So W ≤ 2^126.338 on every run, and no loop needs its own limit;
  rejection-sampling loops in particular are bounded by the meter.
  - Preprocessing (all charges before the first walk step) has its own cap of
    2^59.0 units. If a preprocessing charge would exceed it, the attack halts
    with failure.
  - **Counter widths.** Every counter is stated with the largest value it can
    reach and the word that holds it, and each fits with room to spare.
    - The meter counts milli-units and is a 192-bit counter (three 64-bit words,
      updated by ADDS/ADCS/ADC carry chains). Its largest value is 2^126.338
      units < 2^136.30379 milli-units < 2^192: 55 bits of headroom, and no
      capacity test is needed at 192 bits.
    - The global operation counter of phase 1 counts steps, doublings and
      traversal steps and stops at S_ops, so it reaches at most
      S_ops + 2^60 < 2^126.7408 + 2^60 < 2^126.741: 127 bits, held in a 128-bit
      pair of 64-bit words (ADDS/ADC on the pair). Its stop test compares the
      128-bit pair with S_ops, which is below 2^127, so the test can fire.
    - The start counter counts trail starts. Every start is charged 2^13 units,
      so it can never pass W/2^13 < 2^126.338/2^13 = 2^113.338: 114 bits, held
      in a 128-bit pair. (The expected number of starts is the 2^62.74 of the
      random-charge row of section 5; the counter is sized for the worst case,
      not the expectation.)
    - The preprocessing cap is checked on the meter's milli-unit scale:
      2^59.0 units = 2^68.966 milli-units < 2^69, held in a 128-bit pair and
      compared by a two-word SUBS/SBCS/CSEL sequence. The preprocessing account
      is a separate 128-bit pair, so the cap and the global meter never alias.
    - The block-refill counter, and each processor's private count of the blocks
      it has taken, holds at most W/2^40 + 2^20 < 2^86.34 refills: 87 bits, a
      128-bit pair.
  - Each of ≤ 2^20 processors reserves blocks of 2^40 units by atomic update on
    the 192-bit meter. The refill itself — an atomic 192-bit add of the block, a
    store of the new remainder, a load of the phase-1 limit, the operation count
    and the two limit tests — is charged 30 instructions = 0.03 units per
    refill, the hard row of section 5. At most W/2^40 + 2^20 < 2^86.34 refills
    occur, so that row is < 2^86.338·0.03 < 2^81.28 units. Unused remainders
    count as spent: < 2^61 units, also in the hard extras.
- **Phase 1.** Runs until a global operation counter (steps, doublings and
  traversal steps of phase-1 trails) reaches S_ops. Every reserved block is
  completed, so at least S_ops and fewer than S_ops + 2^60.41 operations are
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
with κ = 1.02 and p0 = 0.3953. H1 also includes three further claims:

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

    S0 = √(2N' ln(1/(1 − p0))) = √(2N' ln(1/0.6047)) = 2^126.71222
    S_ops = 1.02·S0 = 2^126.74079

H1 gives ≥ 0.3953; the pure random-mapping model without any slowdown gives 0.4075. The
toy censuses of section 6 (65,641 to 96,675 trials per curve) measure the slowdown of this
walk against the model as +1.13% ± 0.17% on q25 and +0.22% ± 0.20% on q29, falling with N, and
the earlier runs of section 6 (−0.9% to +1.5%) agree; the allowance κ − 1 is 2%.

**Phase 1 covers the budget.** Phase 1 performs ≥ S_ops operations, and each is
one class arrival. So if T_c ≤ S_ops, the first inter-trail collision happens
in phase 1.

**Doubling and traversal surcharges.** These are random, and only their excess
over the step price is not fixed by the operation count:

- doublings: expected count < 2^-13.5681·(S_ops + 2^60), surcharge
  0.956 − 0.752 = 0.204 each, so expected < 2^110.88;
- traversals: expected < 0.024·64·2^-19.99·(S_ops + 2^60) < 2^107.37.

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
| meter halting early: expected random charges < 2^111.00077 against the headroom 2^118.910483, by Markov | ≤ 0.00416 (exact ratio 0.00415851) |
| preprocessing cap | < 2^-30 |

The events other than the early halt sum to < 2·10^-7.
**Success ≥ 0.3953 − 0.00416 − 2·10^-7 = 0.3911398 > 0.3911 > 0.39**, and with the ledger's
exact Markov ratio 0.00415851 the value is 0.3953 − 0.00415851 − 2·10^-7 = 0.3911412. The two
numbers printed here are the ones section 5 derives: the allowance 2^118.910 against the
expected random charges 2^111.001, and the preprocessing cap's ratio 2^25/2^55.52 < 2^-30.

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

### 5.1 and 5.2 Kernels and counted instructions (abridged)

The field kernels, counted loops and their counts are those of proof.md section 2 (fe_mul 88, fe_sqr 67,
fe_sub 9, fe_subcanon 14; pass 1 108, pass 2 493). The inversion is unchanged: Rinv = acc[n-1]^(p-2) by the
libsecp256k1 addition chain, 255 squarings and 15 multiplications, 270 calls of the non-inlined fe_mul_f and
fe_sqr_f compiled from `kernels/secp_walk.c` into `kernels/secp_fe.s` (115 and 82 instructions with frame and
RET), so I <= 255*(82 + 1) + 15*(115 + 1) + 397 + 1,020 = 24,322 instructions per batch. Batches have
2^15..2^16 walks; the per-batch fixed work is recomputed in proof.md section 4.

### 5.3 Per ordinary step (historical numbers; superseded by proof.md sections 3-5)

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
explicit sequences below, with no allowance on top; section 0.1 lists each one against v8 and
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
- **Doubling** (abridged: the live 793-instruction bound is derived in proof.md section 4 from the
  85/64 kernels with fused x3; the historical 956-instruction derivation for the superseded 91/71
  kernels and 118 rows is omitted here).

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

**Budget** (abridged: the historical phase-1/hard/random tables for the superseded 0.752/0.776/0.956
prices and R = 2048 rates are omitted; the live ledger is proof.md section 5.3 with `certificate.py`).
The live budget keeps the same structure and action prices below; only the step, traversal and doubling
prices, the R = 2^20 rates and the table row move.

**Preprocessing:** table < 2^24.1, node set-up 2^26, bucket initialisation 2^58.865 and the
first 2^20 starts 2^33. The expected total is 2^58.864216, and the cap is 2^59, so the gap is
> 2^55.52. The random sampling part is expected < 2^25, so the cap halts a run with probability
< 2^25/2^55.5 = 2^-30.5 < 2^-30 by Markov (the exact ratio is 2^-30.52). All of it is inside W.

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

**Abridged.** The public package's notes on the provenance of these rows, on its own initial input and on its
limitations are omitted here: they describe the superseded 111 + 523 loops. The current initial run and kernel
evidence are in proof.md section 2.6, the census of the walk as now specified (R = 2^20, checkpoint-only) is in
proof.md section 5.2, and the limitations are in proof.md section 6.

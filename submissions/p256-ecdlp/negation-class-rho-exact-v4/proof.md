# Negation-class parallel Pollard rho on P-256, with a global work meter and an exact report store

**Provenance (declared).** The algorithm, its walk, its analysis and its toy
harness are those of submission 5fa78d62-1c6b-46b8-b9cf-02481a49e645
(`p256-ecdlp/negation-class-rho-v1`, 129.46). That package adapts the
class-walk construction of submission 1c600df4 (mpjunior92,
`bls12-381-g1-ecdlp-plain/glv-six-class-rho-v11`) to P-256, through its
secp256k1 port a188e335, and took its P-256 instruction allowances from
4fe789ec. P-256 has no efficient endomorphism, so the classes are the negation
pairs {P, −P}.

Taken unchanged in substance from 5fa78d62: sections 1 and 2 (algebra, walk,
merging lemma, solving), the heuristic H1 and the toy harness of section 6.

New in this file set: the field arithmetic is priced by an executed AArch64
instruction listing (`field_listing.md`, run by `experiment.py` Part C) instead
of the allowances M = 512, A = 64 and I = 2^19; the distinguished-report path is
an **exact** store (section 3) in which each report is one appended record and
duplicates are found by sorting on the full 256-bit key — there is no slot
index, no home address, no probe sequence, no replacement rule, no load factor
and no occupancy argument anywhere; the meter's block remainders are returned to
the counter instead of being counted as spent; and this section 5 ledger is
re-derived from those two changes.

An earlier version of this file set routed reports into a fixed-size table of
slots with linear probing and argued a 0.6^t occupancy tail. No part of that
argument survives: the store below never places a record at an address derived
from its key, so clustering, home-address bias and the probe cap have no
counterpart, and the only quantity of the store that appears in the failure
analysis is the count of reports against a declared capacity.

**Claim.** time_log2 = 127.72. The algorithm recovers x from Q = [x]P on P-256
with probability ≥ 0.39. On **every** run its total work is W ≤ 2^127.72 units,
where one unit is the borrowed C = 1,000 scalar AArch64 instructions. All P-256
arithmetic is charged in instructions; no P-256 operation is priced as one
unit. The work bound is unconditional (global work meter, section 3). Peak
memory is below 2^72 bytes over 2^16 record logs of below 2^56 bytes each.
Preprocessing (everything before the first walk step) is ≤ 2^60.0 units,
enforced by the meter and included in W. Nonuniform advice: none. Success rests
on the random-mapping heuristic H1 (section 4), supported by toy runs.

## 1. Algebra (checked on P-256 in `experiment.py`, part A)

p, a = −3, b, r and P are the protected SEC 2 parameters. The cofactor is 1 and
r is prime.

- **Units and classes.** The units are u_0 = 1 and u_1 = −1 (no other
  automorphism of E(Fp) is efficiently computable). The group has prime odd
  order, so P ≠ −P for P ≠ O. Hence E \ {O} has exactly N = (r − 1)/2 classes
  {P, −P}, with log2 N = 255.0 − 3.4·10^-10.
- **Key.** K = x(P) is a complete class invariant.
- **Canonical member.** c = (x, |y|) with |y| = min(y, p − y), and t = [p − y < y],
  so c = [u_t]P.

## 2. Walk

Parameters: R = 2048; distinguished probability θ = 2^-64 (limb 1 of K is
zero); M = 2^20 walks; L_max = 20·2^64; checkpoint period 64.

- **Table.** c_h and d_h are uniform in [0, r) (fresh coins, rejection sampling,
  every attempt metered). M_h = [c_h]P + [d_h]Q, redrawn if O.
  Tab[h][t] = [u_t^{-1}]M_h = (−1)^t·M_h, with coefficients ±(c_h, d_h).
- **Step.** F(C) = class(c(C) + M_h) with h = K mod R (bits 0..10 of x). The walk
  keeps any member P and sets P ← P + Tab[h][t]. Since
  c + M_h = [u_t](P + Tab[h][t]), this is F. Each step is one affine P-256
  addition in a batch, charged in instructions (section 5).
- **Trails.** A trail starts at S = [a]P + [b]Q, with a and b uniform from fresh
  coins, stored in full (redraw if O). At each class arrival, in order:

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
  doubling escape P ← [2]P (a batched affine doubling with a = −3). It commutes
  with negation, so it is a class map. fp(K) is limbs 0, 2 and 3 of K (192 bits),
  compared as one integer, and gives the single total order on classes used by
  every equality test and minimum. `length` counts every operation. If an
  addition has equal x-coordinates (P = ±Tab), the trail is abandoned
  (< 2^-114 in total).
- **Merging lemma** (as in the source). Two trails at a common class apply the
  same F until one escapes. Escapes happen only from detected cycles Z of F of
  length ≤ 64, which contain no distinguished class. Every escape goes to
  E(Z) = class(2c(min Z)), which depends only on Z. The other trail enters Z,
  detects it within 3 or 192 arrivals, and escapes to the same class. By
  induction both trails reach the same distinguished class and therefore append
  two records with the same key K, unless a cap or an exceptional addition
  intervenes first (failure events). History is cleared at every doubling, so a
  detected repeat always witnesses a real cycle of F.
- **No structural cycles through escapes.** A return to a class c after m ≥ 1
  doublings needs ±2^m·c + (a sum of ±power-of-2 multiples of table points) = c.
  The order of 2 mod r is 2^255.0 (part A, from the full factorisation of
  r − 1). So 1 ∓ 2^m ≠ 0 mod r for all m < 2^254, and the relation cannot hold
  by cancellation: it pins c to one point. That such returns are no more
  frequent than ordinary self-collisions is part of H1, and the long-trail test
  checks it (section 6).
- **Solving.** The store hands back two records with equal keys, that is two
  trails (a_1, b_1) and (a_2, b_2) that reached the same class. Re-walk both
  from their stored (a, b), tracking coefficients mod r. The canonical
  coefficients are (e, f) = u_t·(A, B). If f_1 ≠ f_2, then
  x = (e_2 − e_1)/(f_1 − f_2) mod r (part A checks the sign on production data).
  Verify [x]P = Q and output x. Q = O outputs 0.

## 3. Phases, the exact report store, and the global work meter

- **Meter.** Every action is charged its price (section 5) before it runs, and
  each price upper-bounds its instructions / C. This includes random words,
  sampling attempts, starts and redraws, steps, doublings, traversals, reports,
  table build, node and router set-up, log sorts, the merge, re-walks, solving,
  storage node and router set-up, log sorts, the merge, matched-record read-back,
  re-walks, solving, output and verification. If a charge would exceed
  W_max = 2^127.72, the attack halts with failure. So W ≤ 2^127.72 on every run,
  and every loop, including every rejection-sampling loop, is bounded by the
  meter. The counter is an integer in milli-units (10^-3 units); every price is
  an integer number of milli-units at least its counted instructions.
  - **Blocks and remainders.** ≤ 2^20 processors reserve blocks of 2^40 units
    by one 64-bit atomic fetch-and-add on the counter (charged 20 instructions)
    and, at the end of each block, return their unspent remainder to the
    counter by a second atomic update (charged 20 instructions). No remainder is
    counted as spent and no remainder allowance appears in the ledger. There are
    at most W/2^40 < 2^87.72 blocks, so the reservations and the returns cost at
    most 2^87.72 × 40 = 2^93.04 instructions < 2^83.08 units, in the hard extras.
  - Preprocessing has its own cap of 2^60.0; exceeding it halts with failure.
- **Phase 1.** Runs until the global operation counter (steps, doublings,
  traversal steps) reaches S_ops. It performs at least S_ops and fewer than
  S_ops + 2^60 operations.
- **Phase 2.** Finishes the open trails; every report is appended before the
  store is read.
- **Caps.** At most 8 processed matches. Report capacity is D_cap = 2^64
  records, the exact range of the store's 64-bit arrival counter.
- **Batches.** 2^15..2^16 ready walks share one Montgomery inversion. A smaller
  batch pays a separate 132-unit inversion.
- **Store: one record per report, appended, never placed by its key.** A report
  is one 104-byte record (32-byte key K, 32-byte a, 32-byte b, 8 bytes of
  flags). The store is 2^16 record logs, one per storage node: each node is a
  separate machine with its own 64-bit address space, holding at most
  ⌈D_cap/2^16⌉ = 2^48 records — 2^54.70 bytes of records plus 2^51 bytes of
  sort index, natively addressable.
  - **Placement.** A global 64-bit arrival counter returns the report's index
    i by one atomic fetch-and-add. The node is i mod 2^16 and the position in
    that node's log is ⌊i/2^16⌋. The counter's range is exactly the capacity:
    i ∈ [0, 2^64), so a 65th bit is never needed and the counter cannot wrap
    before the capacity test halts the run. The counter is the only shared
    object in the store, and its instruction cost is charged per report; the
    cost model sums work across processors, so serialisation adds no uncharged
    work.
  - **Occupancy is a count, not a distribution.** With n reports, node j holds
    exactly ⌈(n − j)/2^16⌉ records: the number of occupied records of every log
    is a function of the report count alone. It is not a function of the keys,
    of any hash, of any modulo, or of any insertion order, so there is no
    occupied-set assumption to make and no occupancy tail to bound. The number
    of reports is at most D_cap except with probability < exp(−2^59.27) (a count
    of distinguished arrivals under H1, in the failure table), and that is the
    only statement about store capacity anywhere in this document.
  - **Detection is exact and offline.** After phase 2, each node radix-sorts
    its own log on the full 256-bit key: 8 LSD passes of a 32-bit digit, one
    counter array of 2^32 entries per node, cleared per pass. The 2^16 sorted
    key streams are then merged by a heap over their heads, comparing full
    256-bit keys, so records with equal keys are adjacent. Two records with
    equal keys are exactly two trails that reached the same distinguished
    class. Nothing else can be reported as a collision: there is no fingerprint,
    no truncation and no hash in the detection path, so no spurious match can
    stop the attack and no genuine match can be missed.
  - **Report price** (instructions, per report):

    | item | instructions |
    |---|---:|
    | arrival counter: one contended 64-bit atomic fetch-and-add | 128 |
    | node select and log position: mask and shift of the index | 3 |
    | route the request over 16 edges: 13 words × 16 instructions + 64 control | 4,352 |
    | append: 13-word record store (7 paired stores and one), write-pointer load, add, store, capacity test, branch | 64 |
    | **total** | **4,547** |

    Charged **4.547 units** (4,547 milli-units). The packet is one-way: a report
    needs no reply, because detection is not online.
  - **Store work after the walk** (instructions; all in the hard extras):

    | item | count | instructions | units |
    |---|---|---:|---:|
    | per-log radix sort on the key plus the adjacent equal-key scan | D_cap records | 208 each | 0.208 |
    | radix counter arrays cleared (8 passes × 2^32 counters × 2^16 logs) | once | 2^35 per log | 2^41.04 |
    | k-way merge of the 2^16 sorted key streams: 4-word key over 16 edges (1,024) plus heap sift-down of depth 16 (512) | D_cap records | 1,536 each | 1.536 |
    | read-back of the ≤ 16 matched records over the same route (16 edges × (13 words × 16 + 64 control) = 4,352 each) | ≤ 16 | 4,352 each | 2^6.12 |

  - **No flush.** The logs are sorted and merged in place: no part of the store
    is read out over the network except the ≤ 16 matched records above. An
    earlier version of this file set charged a full routed read-out in one
    direction only; there is no read-out to charge here, so that asymmetry
    cannot arise.
  - **Node set-up** costs 2^20 instructions per node and per router (2^27 units,
    preprocessing).
- Each trail ends at its first distinguished class or at L_max, so the number
  of reports is at most the number of distinguished arrivals over all phase-1
  and phase-2 walks, which is stochastically dominated by
  Binomial(S_ops + 2^60 + 2^20·(L_max + 65), 2^-64) of mean < 2^63.574 under
  H1. That mean is 0.7445 of D_cap = 2^64, so with δ = D_cap/μ − 1 = 0.34315 a
  Chernoff bound gives
  P(more than D_cap) < exp(−((1+δ)ln(1+δ) − δ)·μ) = exp(−2^59.27)
  (log2 P < −9.96·10^17).

## 4. Success probability

N' = N·R/(R − 1). T_c is the number of class arrivals of all trails, up to the
first arrival at a class already visited by another trail. It includes starts,
steps, doublings, traversal steps and repeated arrivals inside fruitless
cycles; nothing is excluded. This is exactly what the toy measures.

**H1 (heuristic, at the budget).** P(T_c ≤ S) ≥ 1 − exp(−(S/κ)²/(2N')) for the
budget S actually reached, with κ = 1.05 (a 5% slowdown allowance relative to
the random-mapping law).
H1 also includes:

- the index and distinguished fields of x (bits 0..10 and 64..127) are
  pseudo-random functions of the class;
- excluding detected fruitless cycles of length ≤ 64, a trail makes a
  non-fruitless return with probability ≤ L_max²/(2N);
- the twist index t of a new class is uniform and independent of its partition
  index.

The 2-cycle rate is then 1/(2R) per step. Other cycles of length ≤ 64 need
every partition index repeated, at rate ≤ Σ_L Σ_b S2(L, b)·R^-(L−b) < 2^-19.99.

The values are:

    S0 = √(2N' ln(1/0.605)) = 2^127.50399
    S_ops = 1.05·S0 = 2^127.57438

H1 at S_ops gives ≥ 0.395; the pure model gives 0.4254. In the toy runs,
two-class walks are within the 5% allowance of the random-mapping mean T_c
(section 6). Phase 1 stops when the operation counter reaches S_ops, and that
counter counts **doublings** as well as steps and traversal steps; a doubling
adds no class arrival. Doublings occur at rate < 2^-11.99 per operation, so the
arrivals phase 1 reaches are at least

    S_ops − 2^115.58   (the deficit is S_ops·2^-11.99 = 2^115.584)

and the H1 guarantee at that arrival count is 0.394851 instead of 0.395 — the
correction costs 1.5·10^-4. If T_c ≤ S_ops − 2^115.58 the first inter-trail
collision happens in phase 1, and both of its trails are in the store when the
store is read.

**Surcharges.** These are random, and only the excess over the step price is not
fixed by the operation count:

- doublings: rate < 2^-11.99 per operation, surcharge 0.179 units each (1.265 −
  1.086), so expected < 2^113.11;
- traversals: expected < 0.02·64·2^-19.99·S_ops < 2^107.95;
- trail starts: 2^63.7 starts at 2^13 units each, expected < 2^76.7; the meter
  charges every start, so the number of starts is bounded by W/2^13 = 2^114.7
  and needs no separate cap;
- other random charges: random words and sampling attempts, expected < 2^64.

**Failure events** (union bound): if T_c ≤ S_ops − 2^115.58, the merging lemma
gives success unless one of these occurs.

| event | bound |
|---|---|
| non-fruitless return in a colliding trail, per trail (L_max²/(2N)) | 2^-119.36 |
| the same, over the ≤ 8 processed matches | 2^-116.36 |
| either trail of a match is abandoned at L_max before the shared walk reaches a distinguished class (2·e^-20 per match) | < 3.3·10^-8 |
| report count above D_cap | < exp(−2^59.27) |
| f_1 = f_2 | 1/r = 2^-256 |
| exceptional additions | < 2^-114 |
| 8 failing matches | < 2^-60 |
| spurious 192-bit fingerprint match in the walk's cycle tests (arrivals × 2^-192) | < 2^-63 |
| meter halting early: expected random charges < 2^113.15 against headroom > 2^121.94, by Markov | < 0.002239 |
| preprocessing cap: expected sampling < 2^25 against the cap's headroom > 2^59.99, by Markov | < 2^-35 |

The abandonment row is the shared-walk wait, which the earlier version of this
file set did not charge. A trail reports only when it reaches a distinguished
class, and two merged trails walk the same F, so what a match needs is that both
of its trails reach a distinguished class before the length cap. A trail's first
distinguished class arrives at probability 2^-64 per arrival, so it is abandoned
only if its path exceeds L_max = 20·2^64 arrivals, an event of probability at
most (1 − 2^-64)^L_max ≤ e^-20 = 2.06·10^-9. A match fails only if at least one
of its two trails is abandoned, and since that is a union bound over two trails
it needs no independence and no conditioning on the merge: 2·e^-20 per match,
≤ 16·e^-20 = 3.3·10^-8 over the ≤ 8 matches. The superseded row
(8·192/2^64 < 2^-53.4) bounded only the 192-arrival cycle-detection window, not
this wait.

The rows other than the Markov one sum to < 3.8·10^-8. **Success ≥ 0.394851 −
0.002239 − 3.8·10^-8 > 0.39261 > 0.39.**

The margin over the 0.39 requirement is 0.0026, and it is thin on purpose: it is
the H1 floor 0.395, less the 1.5·10^-4 the phase-1 arrival deficit costs, less
the 0.002239 Markov row, less the 3.8·10^-8 tail, and nothing else is claimed.
The pipeline rows of section 6 (0.4273 and 0.4275 over 12,000 and 2,400 fixed
trials) support H1 at toy scale; they are not added to this floor and the floor
is not raised to meet them.

## 5. Cost (all P-256 arithmetic in instructions, C = 1,000)

**Field kernels** (counted from an executed AArch64 listing: every
instruction is in `field_listing.md`, and `experiment.py` Part C executes the
listing with the AArch64 carry-flag rules and checks it on production P-256
elements):

- **M = 139** for a Montgomery multiplication or squaring: 4 operand loads and
  one cleared word (5), the 4×4 column product (77), 4 reduction rounds made of
  shifts and carry chains, since p ≡ −1 mod 2^64 (46), and the final
  conditional subtraction with 2 stores (11). The source itemised the same
  three blocks at ≤ 128, ≤ 80 and ≤ 19 (≤ 227 in all) and charged M = 512;
  this listing fits inside each block, counting the 5 setup instructions with
  the product (82 ≤ 128, 46 ≤ 80, 11 ≤ 19).
- **A = 20** for a modular addition (20 executed) or subtraction (19 executed).
  The source charged A = 64.
- **I = 2^17** for an inversion. Fermat's z^(p−2) as a straight-line chain of
  255 squarings and 127 multiplications executes 382 × 139 = 53,098
  instructions. The source charged I = 2^19.

Every kernel loads its operands from memory and stores its result, so kernels
compose by address and no register is live between them. The kernels are
inlined (no call or spill: 25 of the 28 usable general registers). Field
elements are kept in Montgomery form, which is a fixed bijection. K and |y| are
taken from the canonical Montgomery representation, so they are functions of
the group element.

**Per ordinary step:**

- **Affine addition in a batch: 6M + 6A = 954** (948 executed). That is 3M for
  the Montgomery batch (prefix product and two back-substitution products), 1M
  for λ, 1S for λ² and 1M for y3, plus 6 subtractions. Part C runs a batch of
  production points through the listing with these kernel calls and checks the
  results against reference affine arithmetic.
- **Batch inversion share: I/2^15 = 4.**
- **Control: budget 128, itemised 99.**

  | item | instructions |
  |---|---|
  | sign test, p − y, select | 13 |
  | index, DP test, twist index, table address | 7 |
  | loads of table point and state, store of state | 24 |
  | fingerprint tests, history and checkpoint stores | 26 |
  | zero-denominator test | 6 |
  | length, cap, meter, loop | 14 |
  | batch passes: pointer advance, compare and branch in the forward, backward and update passes | 9 |

- **Total: 954 + 4 + 128 = 1,086 instructions, charged 1.086 units.** The
  price is the count: the control budget of 128 already carries 29 instructions
  of slack over the 99 itemised, so no further round-up is taken.

**Prices** (units; each price is at least its instruction count / 1,000, charged
as an integer number of milli-units):

| action | price | basis |
|---|---|---|
| ordinary step | 1.086 | 1,086 instructions, above |
| traversal step | 1.106 | an ordinary step plus the point copy and minimum test (≤ 20) |
| doubling | 1.265 | λ = (3x² − 3)/(2y) with batch, x3 and y3: 7M + 8A + 4 + 128 = 1,265 instructions |
| random 64-bit word | 0.01 | |
| sampling attempt | 0.05 | |
| trail start | 2^13 | two Jacobian scalar multiplications, ≤ 512 operations × (16M + 16A) each, an addition and an inversion: < 2.8·10^6 instructions |
| report (counter, route, append) | 4.547 | itemised in section 3: 4,547 instructions |
| record sort | 0.208 | 208 instructions per record of capacity |
| key-stream merge | 1.536 | 1,536 instructions per record of capacity |
| router set-up | 2^27 in total | 2^16 nodes and 2^16 routers at 2^20 instructions each |
| phase-2 or re-walk operation | 132.405 | unbatched: I + 7M + 8A + 200 = 132,405 instructions |
| small-batch inversion | 132 | I = 131,072 instructions, rounded up (never used in the bound: phase 2 uses the unbatched price) |
| solving | 600 | an inversion mod r (< 400 multiplications mod r of < 1,000 instructions each), two products and the 32-byte output |
| verification | 2^13 | one scalar multiplication as in a trail start |

**Budget.**

- **Phase 1:** 1.086·(S_ops + 2^60) < 2^127.69341.
- **Hard extras (< 2^95.38 in total):**
  - phase 2: 2^20·(L_max + 65)·132.405 < 2^95.38;
  - meter block reservations and returns: < 2^83.08;
  - re-walks: 16·(L_max + 65)·132.405 < 2^79.38;
  - reports D_cap·4.547 < 2^66.19; merge D_cap·1.536 < 2^64.62; sort
    D_cap·0.208 < 2^61.74; radix counters < 2^41.04; matched-record read-back
    < 2^6.12;
  - router and node set-up, table, solving and verification < 2^28.
  - The nine lines sum to < 2^95.372.
- **Random extras (expected, < 2^113.15 in total):** doubling surcharges
  < 2^113.11, traversals < 2^107.95, starts 2^76.7, and other random charges
  (random words, sampling attempts) < 2^64.
- **Headroom:** 2^127.72 − 2^127.69341 − 2^95.38 > 2^121.94.

Part C recomputes this ledger from the executed counts and prints each term
(upper bounds rounded up, lower bounds rounded down): early-halt probability
≤ 0.002239, the H1 floor at the phase-1 arrival count 0.394851, and success
≥ 0.392611.

**Preprocessing:** table < 2^24.1, node and router set-up 2^27, and the first
2^20 starts 2^33, expected < 2^33.2; the record logs need no initialisation,
because every record position is written before it is read and the only arrays
that are cleared are the radix counters, which are cleared inside the sort and
charged there. The cap is 2^60.0. The random sampling part is expected < 2^25,
so the cap halts a run with probability < 2^-35 by Markov. All of it is inside
W.

**W ≤ 2^127.72 on every run.**

**Memory:** D_cap records × 104 bytes is 2^70.70; two radix index arrays of 8
bytes per record are 2^68.00; the radix counters are 2^51 and the merge heap
2^22.4; the table is 2^18 and the walk states are 2^28. The total is
< 2^70.91 bytes < 2^72. Each node holds at most 2^48 records, 2^54.70 bytes
plus 2^51 bytes of index, < 2^55 bytes, natively addressable.

## 6. Experiments and limitations

`experiment.py` uses the standard library only.

**Part A** checks, on production P-256 (all true):

- curve, order and cofactor 1;
- the full factorisation of r − 1 and the order of 2 mod r (2^255.0);
- orbit size 2, the x key, the unique canonical member;
- the twisted step and doubling as class maps, and the table entries;
- coefficient tracking on production steps;
- a planted-x recovery from two production trails forced to merge.

**Part C** (new; `listing_checks` > 0) executes the listing of
`field_listing.md` on a 64-bit register and memory model with the AArch64 flag
rules. On the hosted input (`listing_checks` = 3000, seed 1) it reports:

- listing lengths M = 139, ADD = 20, SUB = 19, and the opcode histogram of M;
- all three kernels correct on 3,100 cases each (100 edge pairs, then random
  and carry-rich reduced operands), against a·b·2^-256, a + b and a − b mod p;
- a flag-rule mutant (chain-internal `adcs`/`sbcs` demoted to `adc`/`sbc`) is
  detected, so the model enforces the flag rule;
- the inversion chain: 382 multiplications, 53,098 instructions, z·z^-1 = 1;
- one batch of 12 production points (9 additions, 3 doublings) through the
  kernels: results equal the reference affine arithmetic; kernel calls are
  6 M per addition, 7 M per doubling, 382 for the inversion, less 3; executed
  instructions equal 139·M + 20·ADD + 19·SUB calls (64,597);
- the section 5 ledger recomputed from these counts: every price covers its
  count exactly; Phase 1 ≤ 2^127.69341, hard extras ≤ 2^95.372, expected random
  extras ≤ 2^113.143, headroom ≥ 2^121.945, early-halt probability ≤ 0.002239,
  the H1 floor at the phase-1 arrival count 0.394851, success ≥ 0.392611; the
  claim holds;
- the store ledger of section 3 recomputed from its own itemisation: report
  4,547 instructions in four items, sort 208 and merge 1,536 instructions per
  record, read-back 4,352 per matched record, and the count of expected reports
  log2 2^63.574 against the 2^64 capacity with the computed Chernoff exponent
  printed unclamped (about −9.96·10^17, i.e. below 2^-59).

**Part B** runs the per-trail algorithm on prime-order subgroups of toy curves
y² = x³ − 3x + b, with q = 33632567 (b = 13, cofactor 2) and q = 536949043
(b = 25, cofactor 5). The pipeline mode applies the scaled production rules —
phase-1 budget, phase 2, a declared record store with its capacity, at most 8
processed matches — and it builds that store for real: `ReportStore` appends one
record per report to one of 16 logs by arrival index, radix-sorts each log on
the full key with a counter array, merges the sorted logs with a heap comparing
full keys, and solves only from a group of records with equal keys. No
dictionary takes part in the pipeline's collision detection; the T_c census
mode keeps its arrival census in a dictionary as instrumentation, which is what
that mode measures.

**The hosted input, and the offline trial budgets.** Two different inputs are
quoted in this file set and they are never mixed.

*Hosted (the declared initial input, `submission.json`
`initial_inputs[0]`).* It must finish inside the grader's 20-second budget with
its output captured, so its trial count is small:

```
{"curve":"q26","trials":300,"R":2048,"dp_bits":5,"seed":1,"production_checks":3,"listing_checks":3000,"pipeline":true,"seconds":15}
python3 -I experiment.py < hosted.json      # hard timeout 20 s
```

It exits 0 in **3.90 s** on the authoring machine and prints its whole JSON:
Part A true, Part C as listed above, and the pipeline over its own 300 fixed
trials — success 0.5033, 95% Wilson interval [0.4471, 0.5595], 0 log overflows,
maximum log load 0.533 of capacity, 28,792 records. Those 300 trials are a smoke
test that the declared store is built, sorted, merged and solved from; the
interval is wide and this number is not the rate evidence.

*Offline (evidence rows, not the declared initial input).* The 2,000-trial
pipeline result is the OFFLINE budget, whose command is

```
{"curve":"q26","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":3,"listing_checks":3000,"pipeline":true,"seconds":120}
python3 -I experiment.py < offline.json
```

It takes **19.9 s** here — which is why it is not the hosted input, since a
slower grader machine pushes it past 20 s — and gives success 0.45 over its own
2,000 fixed trials, 95% Wilson interval [0.4283, 0.4719], 0 log overflows,
maximum log load 0.6 of capacity, 192,149 records. The 12,000-trial and
2,400-trial rows below give 0.4273 [0.4184, 0.4361] and 0.4275 [0.4078, 0.4474].
In every case the trial count is the stop condition and the `seconds` field is a
guard far above the run's cost, so the reported count does not vary with machine
speed.

Scope of the toy store: it tests the structure and the exactness of the declared
store at toy scale — keys below 2^30, logs of a few thousand records, 16 logs,
8-bit radix digits — and it makes no claim about the production capacity
(2^64 records), the production 32-bit radix digits, or the production
instruction counts, which are analytic. Its bar is exactness: every collision it
solves on is a group of records with equal full keys, and in the runs below no
log reached its capacity.

Offline results. Every row fixes its trial count, so its success rate has a
stated 95% Wilson interval; `T_c` rows are the census (no store), pipeline rows
run the declared store. Seeds are 2501–2507, in table order; the census column
is P[T_c ≤ S_ops] against the model's 0.395.

| mode / classes | toy, R | trials | mean T_c / √N ± s.e. (model) | P[T_c ≤ S_ops] or success (95% Wilson) | 2-cycles per step (model 1/(2R)) |
| --- | --- | ---: | --- | ---: | --- |
| T_c, 2 | q26, 256 | 6000 | 1.2859 ± 0.0086 (1.2558, +2.4%) | 0.4118 (0.3994–0.4243) | 0.001895 (0.0019531) |
| T_c, 2 | q26, 2048 | 6000 | 1.2641 ± 0.0085 (1.2536, +0.8%) | 0.4225 (0.4101–0.4350) | 0.0002363 (0.0002441) |
| T_c, 2 | q30, 256 | 1800 | 1.2730 ± 0.0154 (1.2558, +1.4%) | 0.4161 (0.3935–0.4390) | 0.0018938 (0.0019531) |
| T_c, 2 | q30, 2048 | 1800 | 1.2368 ± 0.0153 (1.2536, −1.3%) | 0.4444 (0.4216–0.4675) | 0.0002320 (0.0002441) |
| T_c, 1 (control, no automorphism) | q26, 256 | 3000 | 1.2620 ± 0.0119 (1.2558, +0.5%) | 0.4243 (0.4068–0.4421) | 0 (0.0019531) |
| pipeline, declared store | q26, 2048 | 12000 | – | 0.4273 (0.4184–0.4361), 0 log overflows, max log load 0.600 of capacity | 0.0002353 |
| pipeline, declared store | q30, 2048 | 2400 | – | 0.4275 (0.4078–0.4474), 0 log overflows, max log load 0.576 of capacity | 0.0002379 |

The census rows are this package's own re-run at fixed trial counts: the mean
T_c is between −1.3% and +2.4% of the random-mapping law, inside the 5%
allowance; the two-class walk's 2-cycle rate is 3.0% to 5.0% below the model's
1/(2R) in all four rows; and the control (no automorphism) shows no 2-cycles at
all. In the 18,600 census trials 39 trails were abandoned at the length cap and
13 additions were exceptional. The abandonment count is the frequency of the
event the failure table bounds. The exceptional-addition count is a toy-modulus
coincidence (equal x-coordinates at q ≈ 2^26) and it does not measure the
production charge, which is 2^-114 at P-256. Every
pipeline trial builds and sorts the declared store; across the 14,400 pipeline
trials 24 trails were abandoned, 5 additions were exceptional, no log reached its
capacity, and the largest log load seen was 0.6 of capacity.

Part A offline (seed 2513, 6 production checks): all true.

**Long-trail trap test** (q30, 150 trails of 2,048 operations without
distinguished points, against the additive-escape control; seeds 2508–2511):

| R | escape | escapes per trail | trails with non-fruitless returns |
| ---: | --- | ---: | ---: |
| 4 | doubling | 149.5 | 1 / 150 |
| 8 | doubling | 103.6 | 8 / 150 |
| 4 | additive (control) | 408.8 | 150 / 150 |
| 8 | additive (control) | 414.1 | 149 / 150 |

Ordinary self-collisions predict non-fruitless returns in about
1 − exp(−L²/(2N)) of trails (L²/(2N) = 0.0391, that is 5.9 trails of 150): the
doubling walk is at 1/150 and 8/150 (9 of 300 trails in all, against 11.8
expected, so within sampling noise of the ordinary self-collision scale), while
the additive escape, which returns to its own past by construction, is at
150/150 and 149/150.

**Limitations.** H1 is a heuristic supported only at toy scale (random
mapping with the Brent–Pollard correction; van Oorschot–Wiener; negation
classes and fruitless cycles as in Wiener–Zuccherato, Bos–Kleinjung–Lenstra
and Bernstein–Lange–Schwabe). The field kernel counts are exact for the
executed straight-line listing, run on a register model rather than on
hardware; they are instruction counts, not cycle timings, which is what the
cost model scores. The control row and the other non-field prices are analytic
itemisations. The store's instruction counts are analytic itemisations of the
two structures it uses (an LSD radix sort and a k-way heap merge), both worst
case per record of capacity; the toy exercises the same structures at toy
scale. This is a generic attack.

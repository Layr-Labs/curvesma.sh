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
and no occupancy argument anywhere; the work meter,
the phase-1 operation count and the store's arrival index are counters whose
largest values are proved to fit their 64-bit words (the counter table of
section 3); the abandonment row is derived for trails selected through a
collision, which are length-biased (section 4); and this section 5 ledger is
re-derived from those changes.

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
enforced by its own meter and included in W. Nonuniform advice: none. Success rests
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
  every equality test and minimum. `length` counts every operation. It is a
  two-word counter (a 64-bit low word and a 64-bit high word counting the low
  word's carries); L_max = 20·2^64 exactly, so `length >= L_max` is the test
  high ≥ 20, and the high word never exceeds 20. If an
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

## 3. Phases, counters, the exact report store, and the global work meter

- **Meter: one global halting condition.** Every action is charged its price
  (section 5) before it runs, and each price upper-bounds its instructions / C.
  This includes random words, sampling attempts, starts and redraws, steps,
  doublings, traversals, reports, the table build, node and router set-up, log
  sorts, the merge, matched-record read-back, re-walks, solving, output and
  verification. Prices are integers in milli-units (10^-3 units). Every loop,
  including every rejection-sampling loop, is bounded by the meter. No word
  ever holds the budget in milli-units: the meter is one shared 64-bit counter
  B of **super-blocks** of 2^70 units, plus private counters inside a block.
  Its cap is

      n_cap = ⌊2^57.72⌋ − 1 = 237,383,948,702,203,524 blocks.

  - **Reservation and halt.** There are at most 2^20 processors in all, storage
    nodes included. A processor reserves a block by one atomic fetch-and-add of
    1 on B. A returned index k < n_cap is its block. A returned index k ≥ n_cap
    halts the run with failure, and that processor executes nothing further. B
    only increases, so B is its own halt signal: every later reservation by any
    processor also returns an index ≥ n_cap, and each processor stops at its own
    first over-cap return. So each processor makes at most one over-cap
    fetch-and-add. Between two reservations a processor spends only from a block
    it already holds, so no action runs on capacity that was not reserved.
  - **Largest value of B.** At most n_cap returns are below n_cap (they are
    distinct indices), and at most one return per processor is at or above it.
    So B ≤ n_cap + 2^20 < 2^58 < 2^64 at the end of every run.
  - **Spending inside a block (private, no atomics).** Each processor keeps a
    64-bit chunk index c and a 64-bit milli-unit counter m for its current chunk
    of 2^50 units (2^50·1000 milli-units < 2^60). An action of price p
    milli-units does m ← m + p. If m ≥ 2^50·1000 it carries (a refill):
    m ← m − 2^50·1000 + 8 and c ← c + 1. The +8 is the carry's own price.
    The carry is 8 instructions (subtract the chunk size held in a register,
    add 8, store m, load, add and store c, compare c with 2^20, branch), and it
    charges them to the chunk it opens, so every refill is paid inside a
    block. When c reaches 2^20 the block is used up, and the processor
    reserves its next block before the action runs; chunk 0 of that block
    starts with the carried remainder and the carry's 8. Every single price is
    below 2^23 milli-units (the largest, a trail start or the verification, is
    2^13 units), so an action crosses at most one chunk boundary,
    m < 2^50·1000 + 2^23 < 2^60 before a carry and m < 2^23 + 8 after it, and
    c ≤ 2^20. The carry strands nothing at a chunk or block boundary.
  - **No remainder returns.** An unspent remainder stays reserved, so each
    processor holds at most one partly used block at any time. When B first
    returns n_cap, the charges already made therefore exceed
    (n_cap − 2^20)·2^70 units. The meter cannot halt a run whose total charges
    stay at or below that. The stranded capacity, at most
    2^20·2^70 = 2^90 units, is taken out of the headroom in the Markov row
    (section 5).
  - **Work bound, exact: the one list of its terms.** Every instruction the
    attack executes is in exactly one of these four terms, and the inequality
    below has exactly these four terms.
    1. *Inside the main blocks*, at most n_cap·2^70·1000 milli-units: every
       action's price (its meter update included), every chunk carry (8,
       charged by the carry itself, above), every phase-1 sub-batch open (8,
       charged on the operation that opens it, below) and every fetch-and-add
       on J (136, charged to its block like an action). A processor spends
       only from blocks it reserved, and at most n_cap reservations succeed.
    2. *The fetch-and-adds on B*, at most n_cap + 2^20 of them (at most n_cap
       below the cap and at most one over-cap per processor), at 136
       instructions each: 128 for the contended atomic, the price the report
       path already uses for its fetch-and-add, plus 1 compare with the cap,
       1 branch and 6 for the private resets (at most three words, a move and
       a store each). They are not charged inside blocks.
    3. *A halted processor's straddling charge.* The meter update that used
       up a processor's last block charged one action's price (below 2^23
       milli-units, a sub-batch open included) and the carry's 8 across the
       block boundary. If the next reservation returns an index ≥ n_cap, the
       instructions that ran after that block was used up (the action's
       down-counter test, sub-batch open and meter update, and the carry) are
       covered only by the part of that charge beyond the block, below
       2^23 + 8 milli-units. At most one per processor: 2^20·(2^23 + 8).
    4. *Preprocessing*, at most 2^60·1000 milli-units, enforced by its own
       meter (below), its fetch-and-adds included.

    So W·1000 ≤ n_cap·2^70·1000 + (n_cap + 2^20)·136 + 2^20·(2^23 + 8) + 2^60·1000, and
    the ledger checks (W·1000)^25 ≤ 1000^25·2^3193 as an exact integer
    inequality. So W ≤ 2^(3193/25) = 2^127.72 on every run. The deterministic
    overheads of terms 1 and 2 (reservations on B and on J, chunk carries,
    sub-batch opens) are also listed as integer totals in the hard extras of
    section 5, so the early-halt row pays for them as well.
  - **Where the private update is paid.** The update (load m, add the price
    from a register, store, compare with the chunk size held in a register,
    branch) is 5 instructions, and it is inside every price: the step's control
    budget itemises it (section 5), the doubling and the traversal share that
    budget, the phase-2 and re-walk price carries 200 control instructions, and
    the start, verification and solving prices exceed their counted
    instructions by far more than 5. Random words and sampling attempts are
    charged together with the attempt or start that consumes them, as one
    action. A chunk carry and a sub-batch open each charge their own 8
    instructions inside the block, and a fetch-and-add on J is charged to the
    block like any action; the fetch-and-adds on B are term 2 of the work
    bound.
  - **Preprocessing** (everything before the first walk step) has a meter of
    the same structure: a shared 64-bit counter of preprocessing blocks of 2^30
    units with cap n_pre = 2^30 − 1 blocks, a private 64-bit milli-unit counter
    inside a block (below 2^30·1000 + 2^23 < 2^41), and at most one over-cap
    fetch-and-add per processor, so that counter ends at most
    n_pre + 2^20 < 2^31. An over-cap return halts with failure before the
    action that asked for the block runs; of that action only its meter
    update (5 instructions, the first instructions of every action) has run.
    Its work, fetch-and-adds and those last updates included, is at most
    n_pre·2^30 + (n_pre + 2^20)·136/1000 + 2^20·5/1000 ≤ 2^60 units (an exact
    inequality in the ledger). Each processor's first
    main reservation is made at the end of its preprocessing and is counted with
    the main fetch-and-adds.
- **Phase 1** ends on a shared 64-bit counter J of **super-batches** of 2^70
  operations (steps, doublings, traversal steps), with cap
  J_cap = ⌈S_ops/2^70⌉ + 2^8 = 214,592,591,910,257,728 (the extra 2^8 batches
  cover the floating-point error of S_ops, which is good to 2^75). A walk
  processor reserves a super-batch by one fetch-and-add on J (136
  instructions, charged to its block). If the returned index is below J_cap,
  the processor performs exactly 2^70 phase-1 operations, counted privately by
  a 64-bit down-counter of the operations left in a sub-batch of 2^50 (at most
  2^50) and a 64-bit sub-batch index (at most 2^20). The down-counter is
  decremented and tested on every operation (4 instructions, itemised in the
  step's control budget), before the operation's meter update. A sub-batch
  open is 8 instructions (reload the down-counter with 2^50 and store it;
  load, add and store the sub-batch index; compare it with 2^20; branch; add
  8 to the price register), so the operation that opens a sub-batch, the
  first of each super-batch included, is charged its price plus 8
  milli-units through its one meter update: the open is paid inside the
  block. If the
  returned index is J_cap or more, that processor's phase 1 is over: its open
  walks go to phase 2 and it makes no further fetch-and-add on J. Every
  reserved super-batch below J_cap is completed, so phase 1 performs exactly
  J_cap·2^70 operations: at least S_ops and fewer than S_ops + 2^79. J ends at
  most J_cap + 2^20 < 2^58.
- **Phase 2.** Finishes the open trails; every report is appended before the
  store is read.
- **Caps.** At most 8 processed matches. Report capacity is
  D_cap = 2^64 − 2^21 records.
- **Batches.** 2^15..2^16 ready walks share one Montgomery inversion. A smaller
  batch pays a separate 132-unit inversion.
- **Store: one record per report, appended, never placed by its key.** A report
  is one 104-byte record (32-byte key K, 32-byte a, 32-byte b, 8 bytes of
  flags). The store is 2^16 record logs, one per storage node: each node is a
  separate machine with its own 64-bit address space, holding at most
  ⌈D_cap/2^16⌉ = 2^48 − 2^5 records — below 2^54.71 bytes of records, below
  2^52 bytes of sort index (two 8-byte index arrays per record) and 2^35
  bytes of radix counters, below 2^54.91 < 2^55 bytes in all, natively
  addressable.
  - **Placement.** A global 64-bit arrival counter returns the report's index
    i by one atomic fetch-and-add. If i ≥ D_cap, the run halts with failure
    and that processor makes no further report. Otherwise the node is
    i mod 2^16 and the position in that node's log is ⌊i/2^16⌋. At most D_cap
    returns are below D_cap, and each of the at most 2^20 reporting processors
    makes at most one return at or above it, so the largest index ever
    returned is at most D_cap + 2^20 − 1 = 2^64 − 2^20 − 1 < 2^64, and the
    counter ends at most 2^64 − 2^20. It never wraps, and the capacity test is
    a test on an exact index. The counter is the only shared object in the
    store, and its instruction cost is charged per report; the cost model sums
    work across processors, so serialisation adds no uncharged work.
  - **Occupancy is a count, not a distribution.** With n reports, node j holds
    exactly ⌈(n − j)/2^16⌉ records: the number of occupied records of every log
    is a function of the report count alone. It is not a function of the keys,
    of any hash, of any modulo, or of any insertion order, so there is no
    occupied-set assumption to make and no occupancy tail to bound. The number
    of reports is at most D_cap except with probability < exp(−2^59.33) (a count
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
    | read-back of the ≤ 16 matched records over the same route (16 edges × (13 words × 16 + 64 control) = 4,352 each) | ≤ 16 | 4,352 each | < 2^6.13 |

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
  Binomial(S_ops + 2^79 + 2^20·(L_max + 65), 2^-64) of mean μ < 2^63.575 under
  H1. That mean is below 0.7446 of D_cap = 2^64 − 2^21, so with
  δ = D_cap/μ − 1 > 0.34315 a Chernoff bound gives
  P(more than D_cap) < exp(−((1+δ)ln(1+δ) − δ)·μ) < exp(−2^59.33).

- **Counter table.** Every counter and index of the attack, shared or private.
  Widths are in bits; the largest value is at the end of the attack, on every
  run.

  | counter | shared or private | width | largest value | updated by | atomic | instructions per update | total charged |
  |---|---|---:|---|---|---|---:|---|
  | meter block counter B | shared | 64 | ≤ n_cap + 2^20 < 2^58 | a processor reserving a block | fetch-and-add | 136 | term 2 of W, outside the blocks: (n_cap + 2^20)·136 instructions < 2^54.85 units; also in the hard extras |
  | meter chunk index c | private | 64 | ≤ 2^20 | its processor, at a carry | no | 8 | inside the block: the carry charges its own 8 (m ← m − 2^50·1000 + 8); a halted processor's last carry is inside term 3 of W; in all (n_cap + 2^20)·2^20·8 instructions < 2^70.76 units, also in the hard extras |
  | meter milli-unit counter m | private | 64 | < 2^50·1000 + 2^23 < 2^60 | its processor, every action | no | 5 | inside every price |
  | preprocessing block counter | shared | 64 | ≤ n_pre + 2^20 < 2^31 | a processor reserving a preprocessing block | fetch-and-add | 136 | inside the 2^60 preprocessing cap |
  | preprocessing milli-unit counter | private | 64 | < 2^30·1000 + 2^23 < 2^41 | its processor, every preprocessing action | no | 5 | inside every price |
  | phase-1 super-batch counter J | shared | 64 | ≤ J_cap + 2^20 < 2^58 | a walk processor reserving a super-batch | fetch-and-add | 136 | inside the block, charged like an action: (J_cap + 2^20)·136 instructions < 2^54.7 units; also in the hard extras |
  | phase-1 sub-batch index | private | 64 | ≤ 2^20 | its processor, at a sub-batch open | no | 8 | inside the block: the opening operation is charged 8 more; in all (J_cap + 2^20)·2^20·8 instructions < 2^70.61 units, also in the hard extras |
  | phase-1 operation down-counter | private | 64 | ≤ 2^50 | its processor, every phase-1 operation | no | 4 | inside the step's control budget |
  | store arrival counter | shared | 64 | ≤ 2^64 − 2^20; largest returned index ≤ 2^64 − 2^20 − 1 | a processor making a report | fetch-and-add | 128 | inside the report price (D_cap·4.547 < 2^66.19 units) |
  | log write pointer, one per node | private to its node | 64 | ≤ 2^48 − 2^5 | its node, at each append | no | inside the append's 64 | inside the report price |
  | radix bucket counts and prefix sums | private to their node | 64 | ≤ 2^48 − 2^5 | its node, during its sort | no | inside the sort's 208 per record | inside the sort |
  | merge heap entries (log number, position) | private to the merge | 64 | log < 2^16, position < 2^48 | the merge | no | inside the merge's 1,536 per record | inside the merge |
  | trail length (low word, high word) | private to its walk | 64 + 64 | low < 2^64, high ≤ 20 | its walk, every operation, re-walks included | no | 4, plus 2 for the cap test | inside the step's control budget and the unbatched price |
  | processed-match counter | private to the solver | 64 | ≤ 8 | the solver | no | inside solving | inside solving |

  There is no other counter: trail starts are bounded by the meter (each is
  charged 2^13 units) and are not counted, and batch sizes are at most 2^16.

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
  index;
- repeats: a trail that reaches L_max operations makes at most L_max/4
  repeated arrivals, that is arrivals at a class that trail has already
  tested (the repeated arrivals inside detected fruitless cycles of length
  ≤ 64 and their traversals). This is a clause of H1, not a derived bound,
  and it is stated only for the trails it is used for, those that reach the
  cap. H1's own cycle rates give an expected rate below
  2^-11 + 2^-11.99 < 2^-10.4 repeated arrivals per operation: 2-cycles at
  1/(2R) = 2^-12 per step with at most 2 repeated arrivals each, other
  cycles at < 2^-19.99 per step with at most 192 + 64 = 256 each. The clause
  allows a rate of 2^-2, about 2^8.4 ≈ 340 times that mean, over the
  20·2^64 ≈ 2^68.3 operations of such a trail.

The 2-cycle rate is then 1/(2R) per step. Other cycles of length ≤ 64 need
every partition index repeated, at rate ≤ Σ_L Σ_b S2(L, b)·R^-(L−b) < 2^-19.99.

The values are:

    S0 = √(2N' ln(1/0.605)) = 2^127.50399
    S_ops = 1.05·S0 = 2^127.57438

H1 at S_ops gives ≥ 0.395; the pure model gives 0.4254. In the toy runs,
two-class walks are within the 5% allowance of the random-mapping mean T_c
(section 6). Phase 1 performs at least S_ops operations (section 3), and that
count includes **doublings** as well as steps and traversal steps. T_c counts
a doubling's arrival, but the bound below credits none, which can only lower
the arrival count it uses. Doublings occur at rate < 2^-11.99 per operation and
phase 1 performs fewer than S_ops + 2^79 operations, so the arrivals phase 1
reaches are at least

    S_ops − 2^115.585   (the deficit is below (S_ops + 2^79)·2^-11.99 < 2^115.585)

and the H1 guarantee at that arrival count is 0.394850 instead of 0.395 — the
correction costs at most 1.5·10^-4. If T_c ≤ S_ops − 2^115.585 the first
inter-trail collision happens in phase 1, and both of its trails are in the
store when the store is read.

**Surcharges.** These are random, and only the excess over the step price is not
fixed by the operation count:

- doublings: rate < 2^-11.99 per operation, surcharge 0.179 units each (1.265 −
  1.086), so expected < 2^113.11;
- traversals: expected < 0.02·64·2^-19.99·S_ops < 2^107.95;
- trail starts: 2^63.7 starts at 2^13 units each, expected < 2^76.7; the meter
  charges every start, so the number of starts is bounded by W/2^13 = 2^114.7
  and needs no separate cap;
- other random charges: random words and sampling attempts, expected < 2^64.

**Failure events** (union bound): if T_c ≤ S_ops − 2^115.585, the merging lemma
gives success unless one of these occurs.

| event | bound |
|---|---|
| non-fruitless return in a colliding trail, per trail (L_max²/(2N)) | 2^-119.35 |
| the same, over the ≤ 8 processed matches | 2^-116.35 |
| either trail of a match is abandoned at L_max before the shared walk reaches a distinguished class (length-biased, derived below: < 9.178·10^-6 per match), over the ≤ 8 processed matches | < 7.342·10^-5 |
| report count above D_cap | < exp(−2^59.33) |
| f_1 = f_2 | 1/r = 2^-256 |
| exceptional additions | < 2^-114 |
| 8 failing matches | < 2^-60 |
| spurious 192-bit fingerprint match in the walk's cycle tests (arrivals × 2^-192) | < 2^-63 |
| meter halting early: expected random charges < 2^113.143 against headroom > 2^121.945 (after the 2^90 stranded units), by Markov | < 0.002239 |
| preprocessing cap: expected sampling < 2^25 against the cap's headroom > 2^59.998 (after the 2^50 stranded units), by Markov | < 2^-34.99 |

**Abandonment of a matched trail (length-biased).** A match needs both of its
trails to reach a distinguished class before L_max. Neither is a fresh trail:
each is selected through the collision, and both selections favour long trails.
Under H1's geometric model each fresh arrival is distinguished independently
with probability θ = 2^-64. Let n be the number of fresh arrivals a trail makes
before the cap. A trail is capped only if it reaches L_max operations, and
then H1's repeat clause leaves at most L_max/4 of them repeated, so
n ≥ n_f = L_max − L_max/4 = 15·2^64; every bound below decreases in n once
nθ > 1, so n_f is used. Write q = 1 − θ. A
trail reaches its i-th fresh arrival with probability q^(i−1).

- *Existing trail* (picked by a visit to one of its classes). The colliding
  arrival lands on a uniformly random visited class, so the index i of the hit
  visit in its trail has weight proportional to the number of trails that
  reached i by then: q^(i−1) times a factor that does not increase with i (a
  trail that reached i by time t also reached i − 1). Given i, the trail is
  capped only if its arrivals i, ..., n are all ordinary, with probability
  q^(n−i+1) ≤ q^(n−i). With weights exactly q^(i−1) on 1..n,

      P_exist ≤ Σ_{i=1..n} q^(i−1)·q^(n−i) / Σ_{i=1..n} q^(i−1) = nθ·q^(n−1)/(1 − q^n).

  A factor that does not increase with i moves weight away from large i, where
  q^(n−i) is largest, so it can only lower this.
- *Incoming trail* (the trail that made the colliding arrival). That arrival is
  its i-th, where i is the age of the active trail that made it: weight
  P(reach i) = q^(i−1), again times a factor that does not increase with i
  early in the run. After the collision the incoming trail follows the existing
  trail's path, whose classes are distinguished independently with probability
  θ under H1, so it is capped with probability at most q^(n−i). The same sum
  gives P_in ≤ nθ·q^(n−1)/(1 − q^n).
- *Comparison.* At n = L_max (nθ = 20) the expression is
  20·(1 − θ)^(L_max−1)/(1 − (1 − θ)^L_max) ≈ 20e^-20/(1 − e^-20)·(1 − θ)^-1:
  the capped-visit fraction of the ideal geometric model times (1 − θ)^-1. The
  extra factor is the hit class's own distinguished test, which the bound does
  not credit (the trail may still be in progress there).
- *Rounding.* q^(n−1) ≤ e^(−(n−1)θ) and 1 − q^n ≥ 1 − e^(−nθ), so with
  x = n_f·θ = 15,

      p_cap ≤ x·e^(−x)·e^θ/(1 − e^(−x)) < 4.589·10^-6 per trail,

  against < 4.123·10^-8 at n = L_max, which is also the ideal model's
  20e^-20/(1 − e^-20) rounded up.
- *Per match.* A match fails by abandonment only if at least one of its two
  trails is capped. That is a union over the two trails and needs no
  independence: ≤ 2·p_cap < 9.178·10^-6 per match, and ≤ 16·p_cap
  < 7.342·10^-5 over the ≤ 8 processed matches.

This is a consequence of H1's ideal random-mapping model, its repeat clause
included; it is not an unconditional bound, and the toy runs support it only at
toy scale (section 6 counts the abandoned trails there). The earlier row, which
applied the fresh-trail tail e^-20 to these selected trails, is withdrawn.

The rows other than the Markov one sum to < 7.35·10^-5. **Success ≥ 0.394850 −
0.002239 − 7.35·10^-5 > 0.392537 > 0.39.**

The margin over the 0.39 requirement is 0.002537, and it is thin on purpose: it
is the H1 floor 0.395, less the 1.5·10^-4 the phase-1 arrival deficit costs, less
the 0.002239 Markov row, less the 7.35·10^-5 tail, and nothing else is claimed.
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
- **Control: budget 128, itemised 104.**

  | item | instructions |
  |---|---|
  | sign test, p − y, select | 13 |
  | index, DP test, twist index, table address | 7 |
  | loads of table point and state, store of state | 24 |
  | fingerprint tests, history and checkpoint stores | 26 |
  | zero-denominator test | 6 |
  | trail length, two words: load pair, add, carry into the high word, store pair | 4 |
  | cap test on the high word (L_max = 20·2^64): compare, branch | 2 |
  | meter, private: load, add the price from a register, store, compare with the chunk size, branch | 5 |
  | loop | 3 |
  | phase-1 operation down-counter, private: load, subtract, store, branch | 4 |
  | price select at a distinguished arrival (step price plus report price, one meter update) | 1 |
  | batch passes: pointer advance, compare and branch in the forward, backward and update passes | 9 |

- **Total: 954 + 4 + 128 = 1,086 instructions, charged 1.086 units.** The
  price is the count: the control budget of 128 already carries 24 instructions
  of slack over the 104 itemised, so no further round-up is taken.

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

- **Phase 1:** 1.086·(S_ops + 2^79) < 2^127.69341.
- **Hard extras (< 2^95.371 in total):**
  - phase 2: 2^20·(L_max + 65)·132.405 < 2^95.38;
  - meter block reservations, (n_cap + 2^20)·136 instructions: < 2^54.85;
  - meter chunk refills, (n_cap + 2^20)·2^20·8 instructions: < 2^70.76;
  - phase-1 super-batch reservations, (J_cap + 2^20)·136 instructions: < 2^54.7;
  - phase-1 sub-batch opens, (J_cap + 2^20)·2^20·8 instructions: < 2^70.61;
  - re-walks: 16·(L_max + 65)·132.405 < 2^79.38;
  - reports D_cap·4.547 < 2^66.19; merge D_cap·1.536 < 2^64.62; sort
    D_cap·0.208 < 2^61.74; radix counters < 2^41.04; matched-record read-back
    < 2^6.13;
  - router and node set-up, table, solving and verification < 2^28.
  - The twelve lines sum to < 2^95.371. The reservations on J, the chunk
    carries and the sub-batch opens are charged inside the blocks and the
    reservations on B are term 2 of W (section 3); all four are listed here
    so that the early-halt row pays for them.
- **Random extras (expected, < 2^113.15 in total):** doubling surcharges
  < 2^113.11, traversals < 2^107.95, starts 2^76.7, and other random charges
  (random words, sampling attempts) < 2^64.
- **Headroom:** n_cap·2^70 − 2^90 (stranded) − 2^127.69341 − 2^95.371
  > 2^121.945. With the single-word meter's accounting (2^127.72 − 2^127.69341
  − 2^95.371, no stranding) it was also > 2^121.945, and the Markov row is
  0.002239 before and after: the 2^90 stranded units and the gap of at most
  2^71 units between n_cap·2^70 and 2^127.72 are below 2^-31 of the headroom
  and do not move the printed digits.

Part C recomputes this ledger from the executed counts and prints each term
(upper bounds rounded up, lower bounds rounded down): early-halt probability
≤ 0.002239, the H1 floor at the phase-1 arrival count 0.394850, and success
≥ 0.392537.

**Preprocessing:** table < 2^24.1, node and router set-up 2^27, and the first
2^20 starts 2^33, expected < 2^33.2; the record logs need no initialisation,
because every record position is written before it is read and the only arrays
that are cleared are the radix counters, which are cleared inside the sort and
charged there. The cap is 2^60.0, enforced by the preprocessing meter of
section 3, whose stranded capacity is at most 2^20·2^30 = 2^50 units. The
random sampling part is expected < 2^25 against a headroom > 2^59.998, so the
cap halts a run with probability < 2^-34.99 by Markov. All of it is inside W.

**W ≤ 2^127.72 on every run.**

**Memory:** D_cap records × 104 bytes is < 2^70.71; two radix index arrays of 8
bytes per record are < 2^68; the radix counters are 2^51 and the merge heap
2^22.4; the table is 2^18 and the walk states are 2^28. The total is
< 2^70.91 bytes < 2^72. Each node holds at most 2^48 − 2^5 records, < 2^54.71
bytes, plus < 2^52 bytes of index and 2^35 bytes of radix counters: < 2^54.91
< 2^55 bytes, natively addressable.

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
  count exactly; Phase 1 ≤ 2^127.69341, hard extras ≤ 2^95.371, expected random
  extras ≤ 2^113.143, headroom ≥ 2^121.945, early-halt probability ≤ 0.002239,
  the H1 floor at the phase-1 arrival count 0.394850, the abandonment rows,
  success ≥ 0.392537; W ≤ 2^127.72 from the four terms of section 3 and every
  counter of the section 3 table below its word, as exact integer
  inequalities; the ledger's arithmetic is consistent;
- the store ledger of section 3 recomputed from its own itemisation: report
  4,547 instructions in four items, sort 208 and merge 1,536 instructions per
  record, read-back 4,352 per matched record, and the expected report count
  below 2^63.575 against the capacity 2^64 − 2^21, with the Chernoff exponent
  ((1+δ)ln(1+δ) − δ)·μ printed (its log2 is above 59.33).

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

It exits 0 in **4.00 s** of wall time (its own clock: 3.91 s) on the authoring machine and prints its whole JSON:
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

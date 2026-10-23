# Negation-class parallel Pollard rho on P-256, priced from a compiled and tested step kernel

**Provenance.** This package starts from submission 5fa78d62 (Meganpark980320,
`negation-class-rho-v1`, declared as a dependency). Its algorithm is used
unchanged: the twisted R = 2048 class walk on {P, −P}, doubling escapes from
detected 2-cycles and checkpoint cycles, the merging lemma, the global work
meter and phases, storage nodes, the success analysis H1 and the toy harness.
That package in turn adapts the class walk of 1c600df4 (mpjunior92) through
a188e335 and 4fe789ec. Sections 1–4 below restate it with only the numbers that
depend on the step price changed.

**What is new here** (section 5 and the lines it feeds):

- A C implementation of the batched P-256 step (`kernels/p256_walk.c`), tested
  against Python big integers and affine P-256 arithmetic.
- The round schedule that uses it (additions and doublings share one inversion)
  and its AArch64 assembly (`kernels/p256_field.s`, `kernels/p256_loops.s`),
  whose per-iteration instruction counts replace the earlier allowances M = 512,
  A = 64 and the 128-instruction control budget. The ordinary step is bounded at
  1,211 instructions and priced at 1,350 (1.35 units) instead of 3,600 (3.6
  units).
- The ledger, random extras, Markov bound and W_max recomputed from that price.
- A `kernels` mode in `experiment.py` that re-derives the counts and the ledger
  inside the verifier.

**Claim.** time_log2 = 128.06. The algorithm recovers x from Q = [x]P on P-256
with probability ≥ 0.39. On **every** run its total work is W ≤ 2^128.06 units,
where one unit is the borrowed C = 1,000 scalar AArch64 instructions. All P-256
arithmetic is charged in instructions; no P-256 operation is priced as one
unit. The work bound is unconditional (global work meter, section 3). Memory is
< 2^73 bytes over 2^16 storage nodes of < 2^57 bytes each. Preprocessing
(everything before the first walk step) is ≤ 2^60.0 units, enforced by the
meter and included in W. Nonuniform advice: none. Success rests on the
random-mapping heuristic H1 (section 4), supported by toy runs.

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
- **Merging lemma.** Two trails at a common class apply the same F until one
  escapes. Escapes happen only from detected cycles Z of F of length ≤ 64, which
  contain no distinguished class. Every escape goes to E(Z) = class(2c(min Z)),
  which depends only on Z. The other trail enters Z, detects it within 3 or 192
  arrivals, and escapes to the same class. By induction both trails reach the
  same distinguished class, unless a cap or an exceptional addition intervenes
  first (failure events). History is cleared at every doubling, so a detected
  repeat always witnesses a real cycle of F.
- **No structural cycles through escapes.** A return to a class c after m ≥ 1
  doublings needs ±2^m·c + (a sum of ±power-of-2 multiples of table points) = c.
  The order of 2 mod r is 2^255.0 (part A, from the full factorisation of
  r − 1). So 1 ∓ 2^m ≠ 0 mod r for all m < 2^254, and the relation cannot hold
  by cancellation: it pins c to one point. That such returns are no more
  frequent than ordinary self-collisions is part of H1, and the long-trail test
  checks it (section 6).
- **Solving.** On a key match, re-walk both trails from their stored (a, b),
  tracking coefficients mod r. The canonical coefficients are (e, f) = u_t·(A, B).
  If f_1 ≠ f_2, then x = (e_2 − e_1)/(f_1 − f_2) mod r (part A checks the sign
  on production data). Verify [x]P = Q and output x. Q = O outputs 0.

## 3. Phases, storage and the global work meter

- **Meter.** Every action is charged its price (section 5) before it runs, and
  each price upper-bounds its instructions / C. This includes random words,
  sampling attempts, starts and redraws, steps, doublings, traversals, reports,
  probes, table build, store initialisation and node set-up, re-walks, solving
  and verification. If a charge would exceed W_max = 2^128.06, the attack halts
  with failure. So W ≤ 2^128.06 on every run, and every loop, including every
  rejection-sampling loop, is bounded by the meter.
  - Preprocessing has its own cap of 2^60.0; exceeding it halts with failure.
  - **Meter protocol (no lock, no timing assumption).** Two 64-bit global
    counters are advanced only by `ldadd`, a single AArch64 instruction
    (contention costs time, not instructions). C_w counts blocks of 2^70 units
    of work, C_o counts blocks of 2^64 operations. A processor with an empty
    work block executes `ldadd` on C_w. If the returned index is ≥
    B_w = ⌊W_max/2^70⌋ it halts with failure. Otherwise it owns 2^70 units,
    deducts the reservation's own cost (≤ 20 instructions, so 0.02 units) from
    that block first, and then charges every action before it runs (a local
    subtraction per round or per rare action). So the charged total, which
    includes the reservations, never exceeds B_w·2^70 ≤ W_max, and there are at
    most 2^58.1 + P reservations for P ≤ 2^20 processors. A block left unused
    when a processor stops is lost (counted as spent): at most P·2^70 ≤ 2^90
    units. Charging an action costs the subtraction and a compare, which are in
    the price of the action.
  - **Operation counts are 128 bits wide.** The trail cap is L_max = 20·2^64 (69
    bits), so each walker keeps its operation count as two 64-bit words
    (`len_lo`, `len_hi`) and the cap is tested as `len_hi ≥ 20`; the checkpoint
    test uses `len_lo`. The C does exactly this and `test_native.py` tests the
    carry and the cap.
- **Phase 1.** Operations are steps, doublings and traversal steps. A
  processor with no operations left in phase 1 executes `ldadd` on C_o. With
  K = ⌈S_ops/2^64⌉, an index k < K obliges it to perform 2^64 operations (it
  finishes them before asking again); an index ≥ K ends its phase 1. So phase 1
  performs exactly K·2^64 operations, at least S_ops and fewer than
  S_ops + 2^64, unless the work meter halts first (an early-halt event, 5.4).
- **Phase 2.** Finishes the open trails.
- **Caps.** At most 8 processed matches. Report capacity is D_cap = 2^64.62.
- **Batches.** 2^15..2^16 ready walks share one inversion; additions and
  doublings share it (the round schedule is in section 5.1). A smaller round
  pays a separate 530-unit inversion.
- **Store.** Reports are routed by bits 128..143 of K (limb 2, disjoint from the
  index bits 0..10 and the DP bits 64..127) to one of 2^16 storage nodes. Each
  node is a separate machine with its own 64-bit address space. It holds an
  open-addressing table of 2·D_cap/2^16 = 2^49.62 slots of 104 bytes, < 2^56.4
  bytes.
  - A report to a full node halts the attack with failure.
  - Under H1 the routing field is pseudo-random. A node's share is
    Binomial(n, 2^-16) with mean < 2^48.6, against 2^49.6 slots. Chernoff gives
    < exp(−2^46) over all nodes.
  - Node set-up costs 2^20 instructions per node (2^26 units, preprocessing).
    Message, routing and node-side work are in the report price.

## 4. Success probability

N' = N·R/(R − 1). T_c is the number of class arrivals of all trails, up to the
first arrival at a class already visited by another trail. It includes starts,
steps, doublings, traversal steps and repeated arrivals inside fruitless
cycles; nothing is excluded. This is exactly what the toy measures.

**H1 (heuristic, at the budget).** P(T_c ≤ S_ops) ≥ 1 − exp(−(S_ops/κ)²/(2N'))
with κ = 1.05 (a 5% slowdown allowance relative to the random-mapping law).
H1 also includes:

- the index, DP and routing fields of x (bits 0..10, 64..127 and 128..143) are
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

H1 gives ≥ 0.395; the pure model gives 0.4254. In the toy runs, two-class walks
are −1.3% to +2.4% from the random-mapping mean T_c, inside the 5% allowance
(section 6). Phase 1 makes ≥ S_ops arrivals, so if T_c ≤ S_ops the first
inter-trail collision happens in phase 1. H1 does not depend on the price of a
step; only the meter's headroom below does.

**Surcharges.** These are random, and only the excess over the step price is not
fixed by the operation count:

- doublings: rate < 2^-11.99 per operation, surcharge 0.45 each (price 1.80
  against 1.35), so expected < 2^114.43;
- traversals: price 1.45, surcharge 0.10 each over 1.35, at most 64 steps per
  cycle event at rate < 2^-19.99 per operation, so expected
  < 0.10·64·2^-19.99·S_ops = 2^110.3.

**Failure events** (union bound): if T_c ≤ S_ops, the merging lemma gives
success unless one of these occurs.

| event | bound |
|---|---|
| trail cap | < 10^-7 |
| non-fruitless return in a colliding trail | < 2^-50 |
| report or node overflow | < exp(−2^46) |
| f_1 = f_2 | 1/r |
| exceptional additions | < 2^-114 |
| 8 failing matches | < 2^-60 |
| spurious 192-bit fingerprint match | < 2^-63 |
| meter halting early: expected random charges < 2^114.52 against headroom > 2^123.25, by Markov | < 0.0024 |
| preprocessing cap | < 2^-30 |

**Success ≥ 0.395 − 0.0024 − 2·10^-7 > 0.3925 > 0.39.**

## 5. Cost (all P-256 arithmetic in instructions, C = 1,000)

### 5.1 The batch schedule and the step kernel

The earlier packages priced a step as 6M + 6A + I/2^15 + 128 control with
allowances M = 512, A = 64 that were twice their itemised counts. Here the step
is written out, tested and counted. `kernels/p256_walk.c` is freestanding C (no
libc, no vector code).

**Field.** Elements are 4×64-bit limbs in Montgomery form, R = 2^256, which is a
fixed bijection. K and |y| are taken from the Montgomery limbs, so they are
functions of the group element and the walk is a valid class map. A product is
the 16 limb products by operand scanning. The reduction uses p ≡ −1 (mod 2^64):
then −p^-1 ≡ 1, so the quotient digit of each round is the limb itself. With
p = (2^64 − 1, 2^32 − 1, 0, 2^64 − 2^32 + 1) as limbs, a round needs one shift
pair for the 2^32 − 1 limb, one `umulh` for the high half of the last limb's
product (the low half is a shift and a subtraction), and carries. The result is
brought to [0, p) by one masked subtraction. Addition and subtraction correct
with a masked add or subtract of p. A squaring is written as a product of equal
operands (S = M in the source; the compiler merges the duplicated partial
products).

**Schedule of one round** (per processor, 2^15..2^16 walkers; 2^20 walks in
total). A walker is either in the addition list (ordinary steps and traversal
steps) or in the doubling queue (its last arrival ended in a detected cycle and
it is due its escape). The slow path puts it there; appending a pointer is part
of the slow path.

1. `pass1` over the addition list and `pass1_dbl` over the doubling queue, with
   one running product (Montgomery one at the start). Per addition walker:
   sign test t = [p − y < y] on the Montgomery limbs, table address 2h + t with
   h = limb 0 of x & 2047, dx = x − xT, zero test, product, and the stores of dx
   and of the prefix product. Per doubling walker: dx = 2y, zero test, product,
   stores.
2. One Fermat inversion of the final product. Additions and doublings share it.
3. `pass2_dbl` over the queue in reverse, then `pass2` over the addition list in
   reverse. Each walker peels its inverse off the batch inverse with two
   products and updates the batch inverse. The first walker of a list uses the
   product before it (the last prefix of the addition list, or Montgomery one).
   An addition walker then does λ = dy/dx, λ², y3 = λ(x − x3) − y with
   x3 = λ² − x − xT. A doubling walker does λ = 3(x − 1)(x + 1)/(2y) (a = −3),
   x3 = λ² − 2x, y3 as above, and clears the history. Both then process the
   arrival at the new class (`arrival` in the C): distinguished test (limb 1 of
   x3), length cap, pending-escape flag, 192-bit equality with K_{-2} and with
   K_cp, the checkpoint store every 64th arrival, the history shift. All the
   events go through one combined branch to the slow path, which owns the
   walker's history from then on (the C returns immediately after it).
4. Slow-path outcomes: a distinguished class is reported and a new trail
   started (outside the batch; the walker rejoins next round); a cap abandons
   the trail; a detected cycle sets the pending flag, starts a traversal or
   enqueues the walker for its doubling. While a traversal runs, the walker's
   flag (bit 0, shared with the pending escape) makes every arrival an event, so
   the slow path does the running-minimum update (`trav_update` in the C: 37
   instructions) plus a hand-off allowed at 60 instructions (flag test, call,
   return, decision): 1,211 + 37 + 60 = 1,308 per traversal step. A walker due
   its doubling is moved from the contiguous addition array to a record in the
   doubling array by swap-remove with the last addition walker (`move_to_dbl`:
   96 instructions, two 168-byte copies written word by word), the queue holds
   pointers into the doubling array, and it is copied back after its doubling
   (`move_back`: 48): 1,431 + 96 + 48 + 60 = 1,635 per doubling with its
   moves.
5. **Zero denominators.** `pass1` counts zero dx (and zero 2y). If the count is
   nonzero the round is discarded before pass 2: those walkers' trails are
   abandoned and restarted, and the round is repeated. By section 2 this
   happens with probability < 2^-114 in the whole run, and a repeated round
   costs at most 2·2^16 steps, so its expected cost is negligible; the meter
   covers it unconditionally.
6. **Per-round fixed work**: the local charge for the round (a subtraction and
   a compare against the processor's work block), the local count of the
   operations left in the current 2^64-operation block, function entry and exit
   of the four loops and the extra products for the first walker of each list:
   < 2^12 instructions. The `ldadd` reservations happen once per 2^70 units of
   work (section 3) and are charged inside their block. With
   the inversion (≤ 2^17) the share per walker step is ≤ 2^17/2^15 + 2^12/2^15 <
   5 instructions. A doubling is in the same round, so it has the same share. A
   smaller round pays a separate 530-unit inversion.

At rate about 2^-12 per step, a round of 2^15 walkers has about 2^3 doublings.
No minimum queue length is needed because the queue shares the round's
inversion.

**Tests** (arm64 macOS, `kernels/test_native.py`, native builds of the same C
at `-O1`, `-O2`, `-O3` (with `-fno-tree-vectorize -mgeneral-regs-only
-ffreestanding -fno-builtin`) and `-Os`; seeds 5, 7, 9, 21):

| check | cases | result |
|---|---:|---|
| fe_mul, fe_add, fe_sub vs big integers (40 × 40 pairs from 18 edge values and 22 random values, then 50,000 random pairs) | 51,600 | all equal |
| fe_inv vs `pow(v, p−2, p)` | 50 | all equal |
| one full round vs Python affine P-256: 100 addition walkers (random class members, table of 2 × 2048 points) and 7 doubling walkers, 4 rounds, one shared inversion | 400 additions, 28 doublings | all equal; history cleared and re-seeded after a doubling; `len` correct |
| arrival events: a planted K_{-2} match, a length cap | 2 | event raised; history left untouched |
| no event | 1 | history updated |

Inside the verifier, `experiment.py` re-runs a pure-Python port of the same limb
algorithms on 2,824 cases (mismatches: 0). The counted object is the same C
compiled for the Linux AArch64 target; I could not execute it here, but the
native builds above use the same flags and the same source.

### 5.2 Instruction counts

Build, as in `COST_MODEL.md` where possible:

    clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only \
          -ffreestanding -fno-builtin -fPIC -S -o p256_walk.s p256_walk.c

`kernels/p256_field.s`, `kernels/p256_loops.s` and `kernels/p256_slow.s` are that
output (Apple clang 21.0.0, LLVM) with comments and assembler directives
removed, split because of the 64 KiB file limit; `kernels/strip_asm.py`
regenerates them byte for byte. `experiment.py` parses them. For a looping function it
reports the longest path, in instructions, of one iteration over all back edges
of the control-flow graph. That path includes the cold block that calls the slow
path (the call counts as one instruction; the callee is priced separately), so
it bounds every iteration.

| function | meaning | instructions in the file | longest iteration |
|---|---|---:|---:|
| `fe_mul` | one multiplication (M); 36 `mul`/`umulh` | 174 | – |
| `fe_add` | modular addition | 36 | – |
| `fe_sub` | modular subtraction | 29 | – |
| `pass1` | addition list, forward | 258 | **236** |
| `pass2` | addition list, backward, with arrival | 1,008 | **970** |
| `pass1_dbl` | doubling queue, forward | 228 | **207** |
| `pass2_dbl` | doubling queue, backward, with arrival | 1,254 | **1,219** |
| `fe_inv` | Fermat inversion | 337 | 309 |
| `move_to_dbl` | swap-remove into the doubling array | 96 | – |
| `move_back` | copy back after the doubling | 48 | – |
| `trav_update` | traversal minimum and step count | 37 | – |

`fe_mul` is the same kind of number as the reference 179 per secp256k1 5x52
multiplication (a call-loop count in a different field). The shipped assembly is
a valid AArch64 program, so these counts are exact for it. Compiling the same
source at other levels and running the same analysis gave the ordinary-step bound
1,212 at `-O1` and `-O2`, 1,211 at `-O3` and 1,185 at `-Os`; the doubling bound
1,431, 1,433, 1,431 and 1,396; the inversion bound 79,132 and 76,827. The experiment
recomputes the `-O3` numbers from the shipped files.

- **Ordinary step**: 236 + 970 + 5 (per-round share) = **1,211 instructions**.
  The 128-bit operation count is included (it costs 3 instructions per arrival).
- **Doubling**: 207 + 1,219 + 5 = **1,431 instructions**.
- **Inversion**: p − 2 has 128 one bits, so left-to-right binary exponentiation
  does 255 loop iterations, each a squaring and a conditional multiplication
  (the longest iteration, 309, contains both): 255·309 + 337 = **79,132**.

**Toolchain.** The reference counts used GCC 14.2.0 and Callgrind. This machine
has neither (the Docker daemon is not running), so the counts come from the
clang output of the shipped source. A different compiler may emit a different
number of instructions for the same C, and the batch has bookkeeping I did not
write in C (the meter, the list appends, the slow-path hand-off). The prices in
5.3 are above the counts for those reasons, not because the shipped program
executes more than counted.

### 5.3 Prices (units)

| action | price | basis |
|---|---|---|
| ordinary step | **1.35** | bound 1,211, margin 11.5% |
| traversal step | 1.45 | ordinary step 1,211 + `trav_update` 37 + hand-off allowance 60 = 1,308; margin 10.9% |
| doubling | **1.80** | doubling 1,431 + `move_to_dbl` 96 + `move_back` 48 + hand-off allowance 60 = 1,635; margin 10.1% |
| inversion (once per round of ≥ 2^15) | 2^17 instructions | bound 79,132 |
| random 64-bit word | 0.2 | ChaCha20 seeded once per processor from a 256-bit true-random key with a distinct nonce: one 512-bit block is ≈ 1,000 instructions for 8 words, so ≈ 125 per word, priced 200 |
| sampling attempt (a scalar mod r) | 1.0 | 4 random words (0.8) + a 256-bit compare with r (≤ 8) + loop (≤ 4) < 0.82 |
| trail start | 2^13 | two scalar multiplications by double-and-add in Jacobian coordinates: 255 doublings at 8 products + 16 additions and 15 control (1,983 instructions) and ≤ 256 mixed additions at 11 products + 14 additions and 10 control (2,428) give ≤ 1.13·10^6 each; plus one addition, one inversion (79,132), the normalisation and the two sampling attempts: < 2.4·10^6 instructions = 2,400 units |
| report | 40 | route (≤ 5), build and send the 104-byte record (13 words, ≤ 40), node-side hash and insert apart from probes (≤ 60), counters (≤ 10): ≤ 200 instructions; the price is 200× that |
| probe | 0.05 | slot load, three key-word compares, index update, loop: ≤ 25 instructions |
| node set-up | 2^26 in total | 2^20 instructions per node for the table allocation and queues, 2^16 nodes |
| store initialisation | 2·D_cap slots × 16 instructions = 2^59.6 units | one 104-byte slot is 7 `stp` plus loop: ≤ 12 instructions |
| table | < 2^24.1, plus sampling | 2,048 entries × (two scalar multiplications + an addition + normalisation, < 2.4·10^6) = 4.9·10^6 units = 2^22.2 |
| phase-2 or re-walk operation | 528.6 | unbatched step: inversion 79,132 + 7 products (1,218) + ≤ 10 additions (360) + control ≤ 400 < 81.2·10^3 instructions = 81.2 units |
| small-round inversion | 530 | inversion bound 79,132 plus calls: < 80 units |
| solving | 600 | one inversion mod r and two products, with a generic 4-limb Montgomery reduction at ≤ 300 instructions per product: < 120 units |
| verification | 2^13 | one scalar multiplication [x]P and the normalisation: 1.21·10^6 instructions = 1,210 units |

The rows from "report" down were already in 5fa78d62 and are kept, each now with
its instruction basis; the random word (a cryptographic generator, not a
simulation one) and the sampling attempt are raised. They matter little: the
unchanged items are the "hard extras" below, < 2^97.4 units in all, and the
"others" in the random extras, < 2^70.

**Slack.** Every counted action is charged at least 10% above its counted bound
(11.5% for a step, 10.9% for a traversal step, 10.1% for a doubling), so the
counted actions cost at most 0.9·W_max in instructions. The remaining items can
therefore be wrong by a large factor without breaking the every-run bound:
if every one of them cost 2^15 times its price, the true total would be at most
0.9·W_max + 2^15·2^97.4 < W_max, because 0.1·W_max = 2^124.7 units is 2^27 times
the hard extras. (This concerns the true cost against the charged cost. The
random extras are expectations and enter only the Markov bound.)

### 5.4 Budget

- **Phase 1:** 1.35·(S_ops + 2^64) = 2^128.00734.
- **Hard extras** (< 2^97.4 in total):
  - lost meter blocks < 2^90 (section 3);
  - phase 2: 2^20·(L_max + 65)·528.6 = 2^97.37;
  - re-walks 2^81.4; reports D_cap·40 = 2^69.9;
  - store initialisation 2^59.6; node set-up, table and verification < 2^27.
- **Random extras** (expected, < 2^114.52 in total): doubling surcharges
  2^-11.99·S_ops·(1.80 − 1.35) = 2^114.43, traversals
  (1.45 − 1.35)·64·2^-19.99·S_ops = 2^110.3, starts 2^76.7, and others < 2^70.
- **Headroom:** 2^128.06 − 2^128.00734 − 2^97.4 > 2^123.25.
- **Markov:** 2^114.52 / 2^123.25 < 0.0024 (0.00233 before rounding).

**Preprocessing:** table < 2^24.1, node set-up 2^26, store initialisation
2^59.6, and the first 2^20 starts 2^33. The expected total is < 2^59.61, and
the cap is 2^60.0. The random sampling part is expected < 2^25, so the cap
halts a run with probability < 2^-30 by Markov. All of it is inside W; the
store initialisation is counted once, in the hard extras above.

**W ≤ 2^128.06 on every run.** The `kernels` mode of `experiment.py` recomputes
S0, S_ops, phase 1, the doubling and traversal terms, the headroom and the
Markov bound to 50 digits; the starts, the "others" and the hard extras are
inputs taken from this section .

**Memory:** 2·D_cap slots × 104 bytes is 2^72.32. Adding the table (2·2048 points
× 64 bytes = 256 KB) and per-walker state (168-byte walker, 32 bytes each for
dx and the prefix product: 232 bytes, taken as 256, for 2^20 walks) gives
< 2^73 bytes. Each node holds < 2^56.4 bytes, natively addressable.

## 6. Experiments and limitations

`experiment.py` uses the standard library only. Parts A and B are the same code
as in 5fa78d62; the `kernels` mode is new.

**Part A** checks, on production P-256 (all true):

- curve, order and cofactor 1;
- the full factorisation of r − 1 and the order of 2 mod r (2^255.0);
- orbit size 2, the x key, the unique canonical member;
- the twisted step and doubling as class maps, and the table entries;
- coefficient tracking on production steps;
- a planted-x recovery from two production trails forced to merge.

**Part B** runs the per-trail algorithm on prime-order subgroups of toy curves
y² = x³ − 3x + b, with q = 33632567 (b = 13, cofactor 2) and q = 536949043
(b = 25, cofactor 5). Pipeline mode applies the scaled production rules,
including routing to 16 nodes with the production per-node capacity
⌈2·D_cap/16⌉.

**What I re-ran and what I carry over.** The hosted initial input
(`curve` q26, R 2048, 2,000 trials requested, `seconds` 9, pipeline, 3 production
checks, kernels) completed in 9 s here: part A all true, pipeline success
0.428–0.430 over 1,130–1,201 completed trials in three runs (the toy loop is
time-bounded, so the trial count varies), 0 node overflows, and the kernels mode
output of section 5. The offline sweep below was produced by the author of
5fa78d62 with the same Part A/B code (seeds 2501–2512, quoted from their note and
proof). I did not repeat it; H1 does not depend on the price of a step.

| mode / classes | toy, R | trials | mean T_c / √N (model) | P[T_c ≤ S_ops] or success | 2-cycles per step (model 1/(2R)) |
| --- | --- | ---: | --- | ---: | --- |
| T_c, 2 | q26, 256 | 6000 | 1.2859 (1.2558, +2.4%) | 0.4118 | 0.00189 (0.00195) |
| T_c, 2 | q26, 2048 | 6000 | 1.2641 (1.2536, +0.8%) | 0.4225 | 0.000236 (0.000244) |
| T_c, 2 | q30, 256 | 1800 | 1.273 (1.2558, +1.4%) | 0.4161 | 0.00189 (0.00195) |
| T_c, 2 | q30, 2048 | 1800 | 1.2368 (1.2536, -1.3%) | 0.4444 | 0.000232 (0.000244) |
| T_c, 1 (control, no automorphism) | q26, 256 | 3000 | 1.262 (1.2558, +0.5%) | 0.4243 | 0 (0.00195) |
| pipeline success with 16 nodes, 2 | q26, 2048 | 12000 | – | 0.4253 (node overflows 0, max load/cap 0.6551724137931034) | 0.000241 |
| pipeline success with 16 nodes, 2 | q30, 2048 | 2400 | – | 0.4221 (node overflows 0, max load/cap 0.5454545454545454) | 0.000237 |

**Long-trail trap test** (q30, 150 trails of 2,048 operations without
distinguished points, against the additive-escape control; carried over):

| R | escape | escapes per trail | trails with non-fruitless returns |
| ---: | --- | ---: | ---: |
| 4 | doubling | 150 | 0 / 150 |
| 8 | doubling | 104 | 4 / 150 |
| 4 | additive (control) | 412 | 150 / 150 |
| 8 | additive (control) | 429 | 149 / 150 |

Ordinary self-collisions predict non-fruitless returns in about 1 − exp(−L²/(2N)) of trails (L²/(2N) = 0.0391).

**Limitations.**

- H1 is a heuristic supported only at toy scale (random mapping with the
  Brent–Pollard correction; van Oorschot–Wiener; negation classes and fruitless
  cycles as in Wiener–Zuccherato, Bos–Kleinjung–Lenstra and
  Bernstein–Lange–Schwabe).
- The instruction counts are exact for the shipped clang output, not Callgrind
  runs of the reference toolchain, and the 11.5% margin over them is a
  judgement, not a measurement of GCC. The tested object is the same C compiled
  natively for macOS arm64; the counted object is the Linux AArch64 assembly,
  which I could not execute. The C is a cost model of the common path: the
  slow paths (cycle traversal, escapes, reports) are priced separately, and
  their logic is tested by the toy harness, not by the C code.
- The step price does not include cache or memory latency: the cost model counts
  instructions, and the loads and stores of the table (256 KB), walker state and
  batch arrays are in the counts.
- The remaining rare-action prices are the earlier generous ones; I did not
  tighten them.
- This is a generic attack and not a practical break.

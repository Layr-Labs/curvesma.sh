# Negation-class parallel Pollard rho on P-256, with a global work meter

**Provenance.** This adapts the class-walk construction of submission 1c600df4
(mpjunior92, `bls12-381-g1-ecdlp-plain/glv-six-class-rho-v11`, declared). That
construction was ported to secp256k1 in our a188e335, also declared. P-256 has
no efficient endomorphism, so the classes are the negation pairs {P, −P}.

Taken from the source:

- the twisted table, the doubling escape, and 2-cycle and checkpoint detection;
- the merging lemma, the global work meter and phases;
- the success-analysis structure and the toy harness.

Our P-256 instruction allowances come from our 4fe789ec (also declared). New
here: P-256 constants and prices, the x class key, report routing to capped
storage nodes, and new toy curves y² = x³ − 3x + b. Everything is restated.

**Claim.** time_log2 = 129.46. The algorithm recovers x from Q = [x]P on P-256
with probability ≥ 0.39. On **every** run its total work is W ≤ 2^129.46 units,
where one unit is the borrowed C = 1,000 scalar AArch64 instructions. All P-256
arithmetic is charged in instructions; no P-256 operation is priced as one
unit. The work bound is unconditional (global work meter, section 3). Memory
is < 2^73 bytes over 2^16 storage nodes of < 2^57 bytes each. Preprocessing
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
- **Merging lemma** (as in the source). Two trails at a common class apply the
  same F until one escapes. Escapes happen only from detected cycles Z of F of
  length ≤ 64, which contain no distinguished class. Every escape goes to
  E(Z) = class(2c(min Z)), which depends only on Z. The other trail enters Z,
  detects it within 3 or 192 arrivals, and escapes to the same class. By
  induction both trails reach the same distinguished class, unless a cap or an
  exceptional addition intervenes first (failure events). History is cleared
  at every doubling, so a detected repeat always witnesses a real cycle of F.
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
  and verification. If a charge would exceed W_max = 2^129.46, the attack halts
  with failure. So W ≤ 2^129.46 on every run, and every loop, including every
  rejection-sampling loop, is bounded by the meter.
  - Preprocessing has its own cap of 2^60.0; exceeding it halts with failure.
  - The meter is a 192-bit fixed-point counter. ≤ 2^20 processors reserve
    blocks of 2^40 units by atomic update. Unused remainders count as spent
    (< 2^61, included).
- **Phase 1.** Runs until the global operation counter (steps, doublings,
  traversal steps) reaches S_ops. It performs at least S_ops and fewer than
  S_ops + 2^60 operations.
- **Phase 2.** Finishes the open trails.
- **Caps.** At most 8 processed matches. Report capacity is D_cap = 2^64.62.
- **Batches.** 2^15..2^16 ready walks share one Montgomery inversion. A smaller
  batch pays a separate 530-unit inversion.
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
(section 6). Phase 1 makes ≥ S_ops arrivals,
so if T_c ≤ S_ops the first inter-trail collision happens in phase 1.

**Surcharges.** These are random, and only the excess over the step price is not
fixed by the operation count:

- doublings: rate < 2^-11.99 per operation, surcharge 0.7 each, so expected
  < 2^115.06;
- traversals: expected < 0.02·64·2^-19.99·S_ops = 2^107.9.

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
| meter halting early: expected random charges < 2^115.08 against headroom > 2^124.18, by Markov | < 0.0019 |
| preprocessing cap | < 2^-30 |

**Success ≥ 0.395 − 0.0020 > 0.39.**

## 5. Cost (all P-256 arithmetic in instructions, C = 1,000)

**Field allowances** (analytic; each is at least twice its itemised count; from
4fe789ec):

- **A = 64** for a modular addition or subtraction (itemised ≤ 32).
- **M = 512** for a Montgomery multiplication or squaring (itemised ≤ 227):
  - a 4×4 column product, ≤ 128;
  - 4 reduction rounds with m = t0, since p ≡ −1 mod 2^64, ≤ 80;
  - a final conditional subtraction, ≤ 19.
- **I = 2^19** for an inversion (Fermat).

Field elements are kept in Montgomery form, which is a fixed bijection. K and
|y| are taken from the canonical Montgomery representation, so they are
functions of the group element.

**Per ordinary step:**

- **Affine addition in a batch: 6M + 6A = 3,456.** That is 3M for the
  Montgomery batch (prefix product and two back-substitution products), 1M for
  λ, 1S for λ² and 1M for y3, plus 6 additions or subtractions.
- **Batch inversion share: I/2^15 = 16.**
- **Control: budget 128, itemised 90.**

  | item | instructions |
  |---|---|
  | sign test, p − y, select | 13 |
  | index, DP test, twist index, table address | 7 |
  | loads of table point and state, store of state | 24 |
  | fingerprint tests, history and checkpoint stores | 26 |
  | zero-denominator test | 6 |
  | length, cap, meter, loop | 14 |

- **Total: 3,600 instructions, 3.6 units.**

**Prices** (units):

| action | price | basis |
|---|---|---|
| ordinary step | 3.6 | |
| traversal step | 3.62 | |
| doubling | 4.3 | λ = (3x² − 3)/(2y) with batch, x3 and y3: 7M + 8A + 128 ≤ 4,224 instructions |
| random 64-bit word | 0.01 | |
| sampling attempt | 0.05 | |
| trail start | 2^13 | two Jacobian scalar multiplications, ≤ 512 operations × (16M + 16A) each, an addition and an inversion: < 4.8·10^6 instructions |
| report | 40 | route, message, node probes and insert: < 2^15 instructions |
| probe | 0.05 | |
| node set-up | 2^26 in total | |
| store initialisation | 2·D_cap slots × 16 instructions = 2^59.6 units | |
| table | < 2^24.1, plus sampling | |
| phase-2 or re-walk operation | 528.6 | unbatched: I + 7M + 8A + 200 |
| small-batch inversion | 530 | |
| solving | 600 | |
| verification | 2^13 | |

**Budget.**

- **Phase 1:** 3.6·(S_ops + 2^60) = 2^129.42237.
- **Hard extras** (< 2^97.37 in total):
  - meter remainders < 2^61;
  - phase 2: 2^20·(L_max + 65)·528.6 = 2^97.37;
  - re-walks 2^81.4; reports D_cap·40 = 2^69.9;
  - store initialisation 2^59.6; node set-up, table and verification < 2^27.
- **Random extras (expected, < 2^115.08 in total):** doubling surcharges
  2^115.06, traversals 2^107.9, starts 2^76.7, and others < 2^70.
- **Headroom:** 2^129.46 − 2^129.42237 − 2^97.37 > 2^124.18.

**Preprocessing:** table < 2^24.1, node set-up 2^26, store initialisation
2^59.6, and the first 2^20 starts 2^33. The expected total is < 2^59.61, and
the cap is 2^60.0. The random sampling part is expected < 2^25, so the cap
halts a run with probability < 2^-30 by Markov. All of it is inside W; the
store initialisation is counted once, in the hard extras above.

**W ≤ 2^129.46 on every run.**

**Memory:** 2·D_cap slots × 104 bytes is 2^72.32. Adding the table (256 KB)
and walk states (2^20 × 256 bytes) gives < 2^73 bytes. Each node holds
< 2^56.4 bytes, natively addressable.

## 6. Experiments and limitations

`experiment.py` uses the standard library only.

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

Offline results (seeds 2501–2512; every recovered x verified; S_ops = 1.05·√(2N' ln(1/0.605)) for each toy; T_c counts all arrivals):

| mode / classes | toy, R | trials | mean T_c / √N (model) | P[T_c ≤ S_ops] or success | 2-cycles per step (model 1/(2R)) |
| --- | --- | ---: | --- | ---: | --- |
| T_c, 2 | q26, 256 | 6000 | 1.2859 (1.2558, +2.4%) | 0.4118 | 0.00189 (0.00195) |
| T_c, 2 | q26, 2048 | 6000 | 1.2641 (1.2536, +0.8%) | 0.4225 | 0.000236 (0.000244) |
| T_c, 2 | q30, 256 | 1800 | 1.273 (1.2558, +1.4%) | 0.4161 | 0.00189 (0.00195) |
| T_c, 2 | q30, 2048 | 1800 | 1.2368 (1.2536, -1.3%) | 0.4444 | 0.000232 (0.000244) |
| T_c, 1 (control, no automorphism) | q26, 256 | 3000 | 1.262 (1.2558, +0.5%) | 0.4243 | 0 (0.00195) |
| pipeline success with 16 nodes, 2 | q26, 2048 | 12000 | – | 0.4253 (node overflows 0, max load/cap 0.6551724137931034) | 0.000241 |
| pipeline success with 16 nodes, 2 | q30, 2048 | 2400 | – | 0.4221 (node overflows 0, max load/cap 0.5454545454545454) | 0.000237 |

Part A offline (seed 2513, 6 production checks): all checks true.

**Long-trail trap test** (q30, 150 trails of 2,048 operations without
distinguished points, against the additive-escape control):

| R | escape | escapes per trail | trails with non-fruitless returns |
| ---: | --- | ---: | ---: |
| 4 | doubling | 150 | 0 / 150 |
| 8 | doubling | 104 | 4 / 150 |
| 4 | additive (control) | 412 | 150 / 150 |
| 8 | additive (control) | 429 | 149 / 150 |

Ordinary self-collisions predict non-fruitless returns in about 1 − exp(−L²/(2N)) of trails (L²/(2N) = 0.0391).

**Limitations.** H1 is a heuristic supported only at toy scale (random
mapping with the Brent–Pollard correction; van Oorschot–Wiener; negation
classes and fruitless cycles as in Wiener–Zuccherato, Bos–Kleinjung–Lenstra
and Bernstein–Lange–Schwabe). Instruction counts are analytic. This is a
generic attack.

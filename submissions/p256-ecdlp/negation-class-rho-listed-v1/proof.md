# Negation-class parallel Pollard rho on P-256, with a global work meter, priced by an executed field listing

**Provenance (declared).** The algorithm, its analysis and its toy evidence are
those of submission 5fa78d62-1c6b-46b8-b9cf-02481a49e645
(`p256-ecdlp/negation-class-rho-v1`, 129.46). That package adapts the
class-walk construction of submission 1c600df4 (mpjunior92,
`bls12-381-g1-ecdlp-plain/glv-six-class-rho-v11`) to P-256, through its
secp256k1 port a188e335, and took its P-256 instruction allowances from
4fe789ec. P-256 has no efficient endomorphism, so the classes are the negation
pairs {P, −P}.

Taken unchanged in substance from 5fa78d62: sections 1 to 4 (algebra, walk,
merging lemma, phases, storage nodes, meter, H1 and the failure analysis), the
experiment's Parts A and B and every toy result in section 6.

New here: the field arithmetic is priced by an executed AArch64 instruction
listing (`field_listing.md`, run by `experiment.py` Part C) instead of the
allowances M = 512, A = 64 and I = 2^19 (3.6 units per step). Section 5 is
re-derived from those counts, and the meter bound, surcharges, headroom and
early-halt probability follow. The control itemisation gains one row (batch
passes). Everything is restated.

**Claim.** time_log2 = 127.75. The algorithm recovers x from Q = [x]P on P-256
with probability ≥ 0.39. On **every** run its total work is W ≤ 2^127.75 units,
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
  and verification. If a charge would exceed W_max = 2^127.75, the attack halts
  with failure. So W ≤ 2^127.75 on every run, and every loop, including every
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
  batch pays a separate 132-unit inversion.
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
    Message, routing and insert are in the report price.
  - The slot index is the 50-bit field at bits 144..193 of K reduced modulo the
    node's slot count (one shift and one reduce, inside the insert price), and
    the table is filled by linear probing. Every probe is charged 0.05 units by
    the meter before it runs, so neither the duplicate search nor the insertion
    rests on an average: a probe sequence ends at the first free slot, and if it
    runs long the meter halts the run like any other charged cost. Under H1 the
    slot field is pseudo-random, so while a node holds at most D_cap/2^16 keys,
    that is at load at most 1/2 (the full-node event above stays at
    exp(−2^46)), a report costs fewer than 3 probes in expectation: the
    unsuccessful-search bound 1/2·(1 + 1/(1 − 1/2)^2) = 2.5 probes at load 1/2
    bounds both the duplicate search and the insertion that follows it, and the
    expected probe work, D_cap·3·0.05 < 2^63 units, is inside the
    "others < 2^70" line of the random surcharge. No expectation is used as a
    bound: the meter bounds the probes on every run, and the chance that they
    exhaust the headroom is in the Markov row.
  - Each trail ends at its first distinguished class or at L_max, so the number
    of reports is at most the number of trails. Under H1 the distinguished field
    is pseudo-random, so the number of reports, which is the number of
    distinguished arrivals, is stochastically dominated by Binomial(S_ops +
    2^60, 2^-64), of mean < 2^63.58 against the capacity D_cap = 2^64.62, and
    Chernoff gives P(more than D_cap) < exp(−2^62).

## 4. Success probability

N' = N·R/(R − 1). T_c is the number of class arrivals of all trails, up to the
first arrival at a class already visited by another trail. It includes starts,
steps, doublings, traversal steps and repeated arrivals inside fruitless
cycles; nothing is excluded. This is exactly what the toy measures.

**H1 (heuristic, at the budget).** P(T_c ≤ S_ops) ≥ 1 − exp(−(S_ops/κ)²/(2N'))
with κ = 1.05 (a 5% slowdown allowance relative to the random-mapping law).
H1 also includes:

- the index, DP, routing and slot fields of x (bits 0..10, 64..127, 128..143
  and 144..193) are pseudo-random functions of the class;
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

- doublings: rate < 2^-11.99 per operation, surcharge 0.2 each, so expected
  < 2^113.27;
- traversals: expected < 0.02·64·2^-19.99·S_ops < 2^107.95.

**Failure events** (union bound): if T_c ≤ S_ops, the merging lemma gives
success unless one of these occurs.

| event | bound |
|---|---|
| trail cap | < 10^-7 |
| non-fruitless return in a colliding trail | < 2^-50 |
| report count or node overflow | < exp(−2^46) |
| f_1 = f_2 | 1/r |
| exceptional additions | < 2^-114 |
| 8 failing matches | < 2^-60 |
| spurious 192-bit fingerprint match | < 2^-63 |
| meter halting early: expected random charges < 2^113.30 against headroom > 2^122.48, by Markov | < 0.0018 |
| preprocessing cap | < 2^-30 |

The rows sum to < 0.0019. **Success ≥ 0.395 − 0.0019 > 0.393 > 0.39.**

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

- **Total: 954 + 4 + 128 = 1,086 instructions, charged 1.1 units (1,100).**

**Prices** (units; each price is at least its instruction count / 1,000):

| action | price | basis |
|---|---|---|
| ordinary step | 1.1 | 1,086 instructions, above |
| traversal step | 1.12 | an ordinary step plus the point copy and minimum test (≤ 20) |
| doubling | 1.3 | λ = (3x² − 3)/(2y) with batch, x3 and y3: 7M + 8A + 4 + 128 = 1,265 instructions |
| random 64-bit word | 0.01 | |
| sampling attempt | 0.05 | |
| trail start | 2^13 | two Jacobian scalar multiplications, ≤ 512 operations × (16M + 16A) each, an addition and an inversion: < 2.8·10^6 instructions |
| report | 40 | route, message and insert: < 2^15 instructions |
| probe | 0.05 | |
| node set-up | 2^26 in total | |
| store initialisation | 2·D_cap slots × 16 instructions < 2^59.66 units | |
| table | < 2^24.1, plus sampling | R starts of at most 2^13 units each |
| phase-2 or re-walk operation | 132.5 | unbatched: I + 7M + 8A + 200 = 132,405 instructions |
| small-batch inversion | 132 | I = 131,072 instructions |
| solving | 600 | an inversion mod r (< 400 multiplications mod r of < 1,000 instructions each), two products and the 32-byte output |
| verification | 2^13 | one scalar multiplication as in a trail start |

**Budget.**

- **Phase 1:** 1.1·(S_ops + 2^60) < 2^127.71189.
- **Hard extras** (< 2^95.38 in total):
  - meter remainders < 2^61;
  - phase 2: 2^20·(L_max + 65)·132.5 < 2^95.38;
  - re-walks < 2^79.38; reports D_cap·40 < 2^69.95;
  - store initialisation < 2^59.66; node set-up, table and verification < 2^27.
- **Random extras (expected, < 2^113.30 in total):** doubling surcharges
  < 2^113.27, traversals < 2^107.95, starts 2^76.7, and others < 2^70 (the
  metered store probes of section 3, < 2^63, are inside this line).
- **Headroom:** 2^127.75 − 2^127.71189 − 2^95.38 > 2^122.48.

Part C recomputes this ledger from the executed counts and prints each term
(upper bounds rounded up, lower bounds rounded down): early-halt probability
≤ 0.001712 and success ≥ 0.393288.

**Preprocessing:** table < 2^24.1, node set-up 2^26, store initialisation
< 2^59.66, and the first 2^20 starts 2^33. The expected total is < 2^59.66, and
the cap is 2^60.0. The random sampling part is expected < 2^25, so the cap
halts a run with probability < 2^-30 by Markov. All of it is inside W; the
store initialisation is counted once, in the hard extras above.

**W ≤ 2^127.75 on every run.**

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
  count; Phase 1 ≤ 2^127.71189, hard extras ≤ 2^95.372, expected random
  extras ≤ 2^113.299, headroom ≥ 2^122.488, early-halt probability ≤ 0.001712,
  success ≥ 0.393288; the claim holds.

**Part B** runs the per-trail algorithm on prime-order subgroups of toy curves
y² = x³ − 3x + b, with q = 33632567 (b = 13, cofactor 2) and q = 536949043
(b = 25, cofactor 5). Pipeline mode applies the scaled production rules,
including routing to 16 nodes with the production per-node capacity
⌈2·D_cap/16⌉.

The toy rows below are the source package's offline results. Parts A and B
are its code unchanged; with `listing_checks` = 0 the hosted input reproduces
the source's local pipeline result exactly (0.433 over 2,000 trials, 0 node
overflows). Part C draws from the same seeded generator, so with it enabled
the toy draws differ.

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
and Bernstein–Lange–Schwabe). The field kernel counts are exact for the
executed straight-line listing, run on a register model rather than on
hardware; they are instruction counts, not cycle timings, which is what the
cost model scores. The control row and the other non-field prices are
analytic itemisations. This is a generic attack.

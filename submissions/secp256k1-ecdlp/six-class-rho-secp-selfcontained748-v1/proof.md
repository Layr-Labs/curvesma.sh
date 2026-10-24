# Six-class secp256k1 rho with the exact 748-instruction step (126.3194)

**Claim.** For the protected full-production secp256k1 target and a uniform
scalar x, the classical cold six-class (negation + GLV) parallel Pollard
rho attack specified below recovers x from Q = [x]P with probability greater
than 0.39003, hence at least the required 0.39, under the explicit H1
heuristic of section 4. On every run its total work is at most
W = 2^126.3194 target-group-operation equivalents, because the global meter
halts before that cap. One unit is C = 1,000 scalar AArch64 instructions.
The ordinary step is charged at its exact 748 instructions = 0.748 units;
no reference operation is charged on top. Peak memory is below 2^72 bytes,
preprocessing is capped at 2^59 units and included in W, and nonuniform
advice is zero. This is an analytic, heuristic, impractical generic attack,
not a demonstrated key recovery.

**Provenance.** The construction is the public six-class rho package
`dddcbeb3-bc7b-4421-b946-70e15cb833d2` (pepedesigner, v9) and its budget
revision `bd8e6be5-913b-47a8-b097-526f4bf4be11`; the kappa=1.01 calibration
and 200,000-trial q29 census are from
`a6ae92da-3a0c-4493-84c1-bcfd6da47232`; the exact-table idea is from
`d14987fb-18f3-4fac-8daf-8f96ea18a156` (126.3470). This document restates
the complete construction, the merging and failure arguments, the kernel
derivation and the ledger, so the claim does not rely on the reader having
those packages. Credit for the construction and kernels belongs to their
authors. The present contribution is the exact 114-row analytic step
accounting and the resulting 126.3194 budget.

The submission archive also contains earlier public packages and is capped
at 262,144 expanded bytes, so the full generated assembly and test files are
not copied here. Sections 5.1 and 5.2 give the self-contained field/batch
algorithm derivation and the exact per-loop instruction counts; a reviewer
who wants the assembly itself can compare the public v9 hashes, which are
`786565c57af8a38458441c1e40521c02cfb18664ad9e65946a4657dfae697bc2`
(`secp_loops.s`) and
`912dc2aa34060ba72109c5a2bb165347506473a9a376d2f8fd9e72f385c6d636`
(`secp_fe.s`). The declared experiment in this package checks the 114-row
sum, the 748 price and the ledger; it does not execute the kernels.

## 0. What changed

The construction, parameters (R = 2048, DP 2^-64, 2^20 walks, L = 20·2^64,
D = 2^63.83), kernels, H1 model and failure estimates are the public
construction. The numerical changes are:

| quantity | public v9 / kappa101 / 126.3470 | this package |
| --- | ---: | ---: |
| ordinary step | 765 / 765 / 763 | **748** |
| counted AArch64 | 634 | 634 |
| analytic rows | 129 | **114** |
| kappa | 1.02 / 1.01 / 1.01 | **1.01** |
| p0 | 0.3953 / 0.3950 / 0.3950 | **0.3950** |
| W | 2^126.3651 / 2^126.3506 / 2^126.3470 | **2^126.3194** |

The 15-instruction reduction is the difference between the conservative
129-row bookkeeping table and the 114 explicit sequences that the public v9
proof itself derives. No algorithm, counted loop, doubling price or
traversal price changes.

## 1. Algebra (checked in the public v9 evidence)

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
  fallbacks, reports, store probes, table build, store initialisation,
  re-walks, solving and verification. If a charge would exceed
  W_max = 2^126.3194, the attack halts with failure. So W ≤ 2^126.3194 on every
  run, and no loop needs its own limit; rejection-sampling loops in particular
  are bounded by the meter.
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
  is D_cap = 2^63.83 records. These bound memory and re-walk work.
- **Batches.** Ready walks are split into batches of 2^15..2^16 per Montgomery
  inversion. A batch of fewer than 2^15 ready walks pays its inversion as a
  separate 600-unit action.
- **Store.** Reports are routed by bits 128..143 of K (limb 2, disjoint from the
  index and DP bits) to one of 2^16 storage nodes. Each node is a separate
  machine with its own 64-bit address space. It holds an open-addressing table
  of 2·D_cap/2^16 slots of 104 bytes, < 2^55.6 bytes. A node whose share
  exceeds its slots halts the attack with failure. Under H1 the routing field
  (bits 128..143 of K, disjoint from the index bits 0..10 and the DP bits
  64..127) is a pseudo-random function of the class. A node's share of n
  reports is then Binomial(n, 2^-16) with mean < 2^47.83, against 2^48.83 slots
  per node. Chernoff gives < exp(−2^45) per node, and < exp(−2^44) over all
  2^16 nodes. The message, routing and node-side work are in the report price.
  Node set-up (process start, table allocation, routing tables) is charged
  2^20 instructions per node, 2^26 units in total, as preprocessing. Slot
  initialisation is charged separately below.

## 4. Success probability

N' = N·R/(R − 1). T_c is the number of class arrivals of all trails, up to the
first arrival at a class already visited by another trail. It includes trail
starts, steps, doublings, traversal steps and repeated arrivals inside
fruitless cycles; nothing is excluded. This is exactly what the toy measures.

**H1 (heuristic, at the budget).** P(T_c ≤ S_ops) ≥ 1 − exp(−(S_ops/κ)²/(2N'))
with κ = 1.01 and p0 = 0.3950. H1 also includes three further claims:

- the partition, DP and node-routing bits of |y| are pseudo-random functions of
  the class: bits 0..10, 64..127 and 128..143, three disjoint fields;
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

    S0 = √(2N' ln(1/(1 − p0))) = √(2N' ln(1/0.6050)) = 2^126.711507
    S_ops = 1.01·S0 = 2^126.725863

H1 gives ≥ 0.3950; the pure random-mapping model at S_ops gives about 0.4011. The public
toy censuses (65,641 to 96,675 trials per curve) measure the slowdown of this
walk against the model as +1.13% ± 0.17% on q25 and +0.22% ± 0.20% on q29, and the
independent 200,000-trial q29 census gives +0.697% ± 0.117%. Section 6 explains the
finite-size origin of these excesses; the allowance κ − 1 is 1%.

**Phase 1 covers the budget.** Phase 1 performs ≥ S_ops operations, and each is
one class arrival. So if T_c ≤ S_ops, the first inter-trail collision happens
in phase 1.

**Doubling and traversal surcharges.** These are random, and only their excess
over the step price is not fixed by the operation count:

- doublings: expected count < 2^-13.56·(S_ops + 2^60), surcharge
  1.1 − 0.748 = 0.352 each, so expected < 2^111.66;
- traversals: expected < 0.047·64·2^-19.99·(S_ops + 2^60) < 2^108.33.

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
| meter halting early: expected random charges < 2^111.795844 against headroom > 2^119.454072, by Markov | < 0.00495044 |
| preprocessing cap | < 2^-32 |

The events other than the early halt sum to < 2·10^-7.
**Success ≥ 0.3950 − 0.00495044 − 2·10^-7 > 0.3900494 > 0.39.**

## 5. Cost

The kernel files named in this section are the public v9 files, pinned by
SHA-256 in the header. They are not copied into this package because the
submission archive also contains the earlier public packages and is capped
at 262,144 expanded bytes. The field and batch algorithm derivation below is
self-contained; a reviewer can compare the pinned public files if desired.
The `experiment.py` in this package is only the ledger check of section 5.4,
not the public kernel test harness.


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

### 5.1 Kernels (public v9 C specification and generator)

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

**The counted loops are generated assembly.** The public v9 `kernels/gen_loops.py` emits the instruction sequences of
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
byte by the public v9 `kernels/build.sh`.

**Tests.** (1) The public v9 `kernels/test_native.py`, on an AArch64 machine (the development machine was Apple silicon;
the generator's `macho` flavor differs only in symbol names and directives), runs the C file and the
generated assembly against Python big integers and affine secp256k1 arithmetic: 60,361 fe_mul_f and
26,019 fe_sqr_f inputs including edge values (0, 1, p−1, p, p+1, 2^256−1, single-limb patterns, values
within 2^40 of p and of 2^256), 317 inversions, 3,837 walks in batches from 2 to 1,000 walks with random
and negated table points through the assembly loops, 172,399 inputs through assembly wrappers of the four
emitters (`t_mul`, `t_sqr`, `t_sub`, `t_canon`: the same emitters, edge values included), and the assembly loops
against the C loops on batches up to 2,000 walks: 0 mismatches everywhere, outputs compared as canonical
integers. (2) The public v9 `experiment.py` (kernels mode) contains a small AArch64 interpreter (about 25 instruction forms)
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
GCC 14.2). The counts below are of the public v9 shipped files, which its `experiment.py` parses.

| function | instructions | what it is |
| --- | ---: | --- |
| pass1 loop iteration | **111** | one fe_mul (91), one fe_sub (9), 5 loads (table pointer, x(T), x(W)), 4 stores (acc_{i-1}, dx_i), loop control (2) |
| pass2 loop iteration | **523** | 4 fe_mul (364), 1 fe_sqr (71), 5 fe_sub (45), 2 fe_canon (16), 19 loads, 6 stores, loop control (2) |
| fe_mul_f (whole function, compiled) | 115 | includes a 12-instruction frame save/restore and RET |
| fe_sqr_f (whole function, compiled) | 82 | includes a 2-instruction frame and RET |
| pass1 / pass2 whole function | 142 / 879 | prologue, first element, epilogue |

Both loops are single basic blocks: every conditional is a CSEL, there is no internal branch, and
The public v9 `experiment.py` verifies that the only branch in each loop is the back edge. Every iteration therefore
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
### 5.3 Per ordinary step: the 114 analytic rows

The counted loops of section 5.2 give 634 instructions. The remaining class
bookkeeping is the following explicit sequence, not an allowance:

| row | instructions | explicit sequence |
| --- | ---: | --- |
| F1 truncated product | 24 | 4 x-limb loads + 5 G-limb loads (9); 4 MUL + 4 UMULH (8); 7 ADD (7) |
| range test, F0, S, F2 | 8 | SUB, CMP, branch (3); ADDS (1); NEG (1); CSEL, CMP, branch (3) |
| decision | 13 | CMP, 2 CSEL, CSET (4); CMP, CSEL (2); CMP, 3 CSEL (4); difference, CMP, branch (3) |
| sign, p−y, K | 13 | SUBS+3 SBCS (4); SUBS+3 SBCS (4); CSET (1); 4 CSEL (4) |
| h, DP test, t, address | 7 | AND (1); CBZ (1); 2 ADD (2); 3 ADD with shifts (3) |
| arrival-pass traffic | 5 | 4 y-limb loads + 1 table-pointer store |
| fingerprints, history, checkpoint | 26 | flag load/CBNZ (2); slot address (2); compare 3 loads/CMP/2 CCMP/branch (7); overwrite 3 stores (3); checkpoint compare 3 loads/CMP/2 CCMP/branch (7); TST/branch (2); checkpoint stores (3) |
| zero-denominator test | 4 | ORR of dx limbs in pairs and together, CBZ |
| length, cap, meter, loop | 13 | 2 length loads, ADDS+ADC, 2 stores (6); cap CMP/branch (2); meter SUBS/ADD/branch (3); loop SUBS/branch (2) |
| batch fixed-work share | 1 | per-batch prologue, inversion glue and epilogue below 2^15 instructions for at least 2^15 walks |
| **analytic total** | **114** | |

The fallback (0.5 units), cold zero-denominator/abandonment path, doubling
escape (1.1), traversal step (0.795) and reports (40) are separately priced
and are not part of the 748. The table point and walk point loads/stores are
inside the counted loops; the arrival pass adds only the y limbs and the
table pointer that the counted loops do not read or write. The 13-row's loop
SUBS/branch is the arrival-pass loop, not a second count of the counted
loops. The public v9 table (27, 12, 14, 13, 7, 8, 26, 6, 14, 2) sums to 129;
its own derived sequences sum to the 114 above, with 15 instructions of
slack. The 126.3470 package removed the upstream two-instruction reserve and
charged 763. This package charges 634 + 114 = 748.

**Prices** (units): ordinary step 0.748 (748 instructions); traversal 0.795
(at most 772 instructions, conservative); doubling 1.1 (at most 952
instructions with the 114 rows: 5 fe_mul, 2 fe_sqr, one 2y addition, two
weak 3x² additions, four subtractions, two canonicalisations, 72 traffic and
loop instructions, 114 rows, 50 list moves, 30 escape-logic instructions);
fallback 0.5; report 40; start 2^13; phase-2/re-walk operation 600;
verification 2^13. These are the public v9 prices, unchanged.

### 5.4 Budget

- Phase 1: 0.748·(S_ops + 2^60) < 2^126.306973.
- Hard extras H: meter remainders < 2^61; phase 2
  2^20·(L+65)·600 = 2^97.55; re-walks 2^81.6; reports 2^69.2; store
  initialisation 2^58.9; table and verification < 2^25; total < 2^97.550769.
- Random extras E: doubling surcharges
  0.352·2^-13.56·(S_ops+2^60) < 2^111.66; traversal surcharges
  0.047·64·2^-19.99·(S_ops+2^60) < 2^108.33; starts 2^75.9; words,
  attempts, probes, small batches, fallbacks 2^70; total
  < 2^111.795844.
- Headroom: 2^126.3194 − 2^126.306973 − 2^97.550769 > 2^119.454072.
- Early halt: E/headroom < 0.00495044 by Markov.
- Success: 0.3950 − 0.00495044 − 2·10^-7 = 0.390049361 > 0.39003.
- W ≤ 2^126.3194 on every run.
- Preprocessing: table < 2^24.1, node set-up 2^26, store init 2^58.9,
  first 2^20 starts 2^33; total < 2^58.91 under the 2^59 cap.
- Memory: 2·D_cap slots × 104 bytes = 2^71.53; table and 2^20 × 512-byte
  walk states keep the total below 2^72 bytes; each of 2^16 nodes is below
  2^56 bytes.

The score is the exact rational 1263194/10000. The next lower four-decimal
cap, 126.3193, gives 0.3900090292, below the claimed 0.39003 floor, and is
deliberately not claimed.

## 6. H1 calibration, evidence and limitations

The public censuses of the exact class walk (R = 2048, dp_bits = 5) are:

| curve | trials | mean T_c/√N | model | slowdown |
| --- | ---: | ---: | ---: | ---: |
| q25 | 96,675 | 1.2677 | 1.2536 | +1.13% ± 0.17% |
| q29 | 65,641 | 1.2563 | 1.2536 | +0.22% ± 0.20% |
| q29 independent | 200,000 | 1.26236 | 1.2536 | +0.697% ± 0.117% |

The excesses are consistent with finite-size self-collisions. A toy trail
has expected length l = 2^5 and a trial uses T_c ≈ 1.26√N arrivals, so it
meets its own past with probability about T_c·l/N: 1.7% at q25 and 0.74% at
q29, matching the measured excesses. At production T_c < 2^126.73, l = 2^64
and N ≈ 2^253.42, so T_c·l/N < 2^-62; the term that dominates the toys is
absent, and the 1% allowance is conservative. Individual q29 seeds range
from +0.08% to +1.75%, so the one-curve +0.22% row is within the spread.
H1 remains a disclosed heuristic; the census is evidence only at its
executed toy scale. The 200,000-trial census is the one in the public
kappa101 package (`a6ae92da`), copied here as `census-q29.json`.

The limitations are the public limitations: the counted instruction counts
are those of the shipped generated assembly, checked natively and by an
interpreter but not by a Callgrind measurement of a full attack; the
analytic rows are explicit sequences but are analytic, not measured; H1 is
a random-mapping heuristic; and the attack is a generic analytic bound, not
a practical break. The 15-instruction tightening is justified by the
explicit 114 sequences; a reviewer should check that no ordinary-path
instruction is missing. The fallback, doubling, traversal and report paths
are separately priced, so a missing instruction there is not hidden by the
114. The counted 634 is not re-executed in this package; its derivation is
sections 5.1–5.2, and the public assembly hashes are in the header. The
declared experiment (`experiment.py`) recomputes the 114-row sum, the 748
price, the 126.3194 cap and the success floor with 80-digit Decimal
arithmetic. Docker is unavailable in the local preparation environment, so
the isolated Docker run was not performed here; the hosted review remains
necessary.

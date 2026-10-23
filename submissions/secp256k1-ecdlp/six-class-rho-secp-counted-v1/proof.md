# Six-class parallel Pollard rho on secp256k1: the step priced from compiled AArch64 kernels

**Provenance.** This package keeps the walk, class key, canonical-index rule, cycle handling,
work meter, phases, storage, success analysis and the toy experiments of public submission
3682292 (six-class-rho-secp-147-v1, the 127.00 bound). That package derives from 758d8df6
(127.01) and 189ebc42, and ultimately from 1c600df4; all are credited and declared as
dependencies. The complete argument is restated here; no prior acceptance or score is used as
a premise. What differs is the price of an ordinary step: instead of buying one 1,000-instruction
reference addition and adding 147 instructions, the batched affine addition is implemented in C
(`kernels/secp_walk.c`, 4x64-bit limbs, the special form p = 2^256 - 2^32 - 977), tested against
big integers, compiled for AArch64, and charged at its counted instructions.

**Claim.** time_log2 = 126.65. The classical cold attack recovers x from a uniform secp256k1
point Q=[x]P with probability at least 0.39 under the explicit H1 heuristic in section 4. On
every run total work is at most 2^126.65 target-group-operation equivalents because the global
meter stops before that bound. The unit is 1,000 scalar AArch64 instructions (C = 1,000, fixed
by the organizer); no reference operation is charged for ordinary steps, see section 5.0. Peak
memory is below 2^72 bytes, preprocessing is capped at 2^59 units and included, and there is no
nonuniform advice. This is an analytical, impractical generic attack, not a demonstrated key
recovery.

## 0. Change from 3682292

| item | 3682292 | this package |
| --- | --- | --- |
| ordinary step | 1 reference addition (1,000) + 147 = 1.147 units | 750 counted instructions of kernels + 129 bookkeeping rows = 879, charged 890 = 0.890 units |
| traversal step | 1.177 | 0.920 (903 instructions bound) |
| doubling | 2.0 | 2.0 (unchanged, now looser) |
| inversion bound | I < 2^17, 4 instructions per walk | I <= 24,322, whole per-batch overhead < 2^15, 1 instruction per walk |
| phase-1 cost | < 2^126.97977 | < 2^126.6138 |
| W_max | 2^127.00 | 2^126.65 |
| early-halt (Markov) | < 0.0046 | < 0.0042 |
| success floor | > 0.3903 | > 0.3907 |

Everything in sections 1 to 4 below is unchanged except the numbers that depend on the step
price (the surcharges, headroom, early-halt bound, W_max). The experiment file adds an
optional kernels check; its walk and seeded toy results are those of 3682292.

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
  W_max = 2^126.65, the attack halts with failure. So W ≤ 2^126.65 on every
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
with κ = 1.05. H1 also includes three further claims:

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

    S0 = √(2N' ln(1/0.605)) = 2^126.71151
    S_ops = 1.05·S0 = 2^126.78190

H1 gives ≥ 0.395; the pure model gives 0.4254. In the toy runs, class walks are
−0.9% to +1.3% from the random-mapping mean T_c, which is within the 5%
allowance (section 6).

**Phase 1 covers the budget.** Phase 1 performs ≥ S_ops operations, and each is
one class arrival. So if T_c ≤ S_ops, the first inter-trail collision happens
in phase 1.

**Doubling and traversal surcharges.** These are random, and only their excess
over the step price is not fixed by the operation count:

- doublings: expected count < 2^-13.56·(S_ops + 2^60), surcharge
  2.0 − 0.890 = 1.11 each, so expected < 2^113.38;
- traversals: expected < 0.03·64·2^-19.99·(S_ops + 2^60) < 2^107.74.

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
| meter halting early: expected random charges < 2^113.41 against headroom > 2^121.31, by Markov | < 0.0042 |
| preprocessing cap | < 2^-32 |

The events other than the early halt sum to < 2·10^-7.
**Success ≥ 0.395 − 0.0042 − 2·10^-7 > 0.3907 > 0.39.**

## 5. Cost

### 5.0 What is charged

COST_MODEL.md prices a computation as H reference operations plus R other
instructions at R/C, with C = 1,000 scalar AArch64 instructions per unit, and says
"Faster arithmetic and lookup algorithms are welcome; justify their actual work" and
"Do not charge both a reference operation and its internal instructions". This package
therefore charges **H = 0 reference operations** for ordinary steps. The affine addition
is not bought as a unit: it is implemented in C, compiled for AArch64, and its executed
instructions are counted and divided by C. The organizer's 6M reference unit (3M for the
addition formula, 3M for batch amortisation, S = M, A = 0) is the same computation with
every non-multiplication instruction set to zero; here the multiplications cost what the
compiled code costs and the additions, subtractions, loads, stores and loop control
are charged too. Nothing is charged twice.

### 5.1 Kernels (`kernels/secp_walk.c`, freestanding C, no libc)

p = 2^256 − C0 with C0 = 2^32 + 977 = 4,294,968,273, so 2^256 ≡ C0 (mod p). A field
element is four little-endian 64-bit limbs. **Weak** means any value in [0, 2^256);
**canonical** means [0, p). Walk points and table points are always canonical.

- **fe_mul(a, b)**, a and b weak, returns a weak r ≡ ab. Row-wise schoolbook product
  (16 MUL, 16 UMULH, carry chains of ADDS/ADCS), then z = L + H·2^256 ≡ L + H·C0:
  H·C0 < 2^256·C0 is five limbs (4 MUL, 4 UMULH), s = L + H·C0 has a fifth limb
  s4 ≤ C0 (so s < 2^256 (1 + C0)), and s ≡ (s mod 2^256) + s4·C0 with s4·C0 ≤ C0² < 2^65
  (two limbs w0, w1, w1 ≤ 1). The sum t = (s mod 2^256) + w is below 2^256 + 2^65; if it
  carries out, t − 2^256 < 2^65, and adding C0 once more stays below 2^256, so no second
  carry exists. Every limb addition of a carry bit into a limb that the product has not
  yet touched is an addition to a partial sum of ab < 2^512, so it cannot overflow.
- **fe_sqr(a)**: the six off-diagonal products, doubled by a one-bit shift across the
  limbs (2·offdiag ≤ a² < 2^512, no bit is lost), plus the four diagonal squares; then
  the same reduction. 10 MUL + 10 UMULH instead of 16 + 16.
- **fe_sub(a, b)**, a weak, **b canonical**: a 4-limb SUBS/SBCS chain; on borrow-out
  (a < b ≤ p − 1) the wrapped value a − b + 2^256 is at least 2^256 − b ≥ C0 + 1, so
  subtracting C0 (the same as adding p modulo 2^256) cannot borrow again and the result
  is a − b + p ∈ (0, p). Without borrow the result is a − b. The result is ≡ a − b and
  is canonical whenever a is canonical.
- **fe_canon(r)**, r weak: r + C0 carries out of 2^256 exactly when r ≥ p, and then
  (r + C0) mod 2^256 = r − p < C0; a conditional select (CSEL, no branch) picks it.
- Mixed operands. The bound for fe_sub needs only its second operand canonical, and
  every use has that: dx = x_T − x_W and dy = y_T − y_W (both canonical), x3' = λ² − x_W
  (x_W canonical), x3'' = x3' − x_T (x_T canonical), x_W − x3 and y3' = λ·(x_W − x3) − y_W
  (second operands x3, y_W canonical after fe_canon).

**Pass 1** (forward, n ≥ 2 walks W_i = (x, y) and table points T_i = *tp[i], all
canonical, x(T_i) ≠ x(W_i)): dx_i = x(T_i) − x(W_i) (canonical), acc[i] = acc[i−1]·dx_i
(weak), acc[0] = dx_0. Because 0 < dx_i < p and p is prime, acc[n−1] ≢ 0.
**Inversion**: Rinv = acc[n−1]^(p−2) by the libsecp256k1 addition chain: 255 squarings and
15 multiplications, 270 calls of two non-inlined functions fe_mul_f and fe_sqr_f.
**Pass 2** (backward, i = n−1 .. 1; R starts as Rinv ≡ 1/(dx_0···dx_{n−1})):
inv = R·acc[i−1] ≡ 1/dx_i (one fe_mul), R ← R·dx_i (dx_i recomputed by fe_sub; one fe_mul),
then dy = y(T_i) − y(W_i); λ = dy·inv; λ²; x3 = λ² − x(W_i) − x(T_i), canonicalised;
y3 = λ·(x(W_i) − x3) − y(W_i), canonicalised; W_i ← (x3, y3). The last walk i = 0 uses
R itself as its inverse. This is the standard affine addition λ = (y_T − y_W)/(x_T − x_W),
x3 = λ² − x_W − x_T, y3 = λ(x_W − x3) − y_W, so W_i ← W_i + T_i with canonical output.
Per walk: 5 multiplications (1 in pass 1, 4 in pass 2) and 1 squaring, i.e. the organizer's
6M, plus seven fe_sub, two fe_canon and the memory traffic.

**Tests.** `kernels/test_native.py` runs the C file (compiled with the same flags, plus
a native build) against Python big integers and affine secp256k1 arithmetic: 60,361
fe_mul_f and 26,019 fe_sqr_f inputs including edge values (0, 1, p−1, p, p+1, 2^256−1,
single-limb patterns, values within 2^40 of p and of 2^256), 317 inversions, and 3,837
walks in batches from 2 to 1,000 walks with random and negated table points; 0
mismatches, outputs compared as canonical integers. `experiment.py` (kernels mode)
repeats a smaller version with limb-for-limb Python ports that assert the no-overflow
claims above (s4 ≤ C0, w1 ≤ 1, no carry out of the last limb), and checks a 24-walk
batch round against affine arithmetic.

### 5.2 Counted instructions

Build (needs clang): `clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize
-mgeneral-regs-only -ffreestanding -fno-builtin -fPIC -S secp_walk.c`, then
`strip_asm.py` splits the output into `secp_loops.s` (pass1, pass2) and `secp_fe.s`
(fe_mul_f, fe_sqr_f, fe_inv), keeping instructions and labels only. The compiler is
Apple clang 21.0.0 (clang-2100.0.123.102) cross-compiling to AArch64 Linux; these are the
COST_MODEL flags with clang in place of GCC 14.2; the counts below are of the shipped
files, which `experiment.py` parses.

| function | instructions | what it is |
| --- | ---: | --- |
| pass1 loop iteration | **128** | one fe_sub, one fe_mul, loads of the operands, 4 stores, loop control |
| pass2 loop iteration | **622** | 4 fe_mul, 1 fe_sqr, 6 fe_sub, 2 fe_canon, loads/stores, loop control |
| fe_mul_f (whole function) | 115 | includes a 20-instruction frame save/restore and RET |
| fe_sqr_f (whole function) | 82 | idem |
| pass1 / pass2 whole function | 173 / 1,026 | prologue, first element, epilogue |

Both loops are single basic blocks: every conditional is a CSEL, CINC or CSET, there
is no internal branch, and `experiment.py` verifies that the only branch in each
loop is the back edge. Every iteration therefore executes exactly 128 and 622
instructions, whatever the data. Opcode histograms: pass1 — ADCS 27, MUL 21, UMULH 21,
ADDS 11, CINC 9, LDP 6, LDR 5, SUBS 5, CSET 5, SBCS 3, ADD 3, STP 3, STR 2, CSEL 2,
SUB 1, MOV 1, MOVK 1, ADC 1, B.NE 1. pass2 — ADCS 132, MUL 99, UMULH 99, ADDS 54,
CINC 41, CSET 28, LDR 25, SUBS 24, CSEL 19, LDP 18, SBCS 18, STR 17, STP 14, SUB 9,
ADC 6, ADD 5, EXTR 5, MOV 4, MOVK 2, LSR 1, STUR 1, CBNZ 1. Of these, 7 (pass1) and 59
(pass2) are register spills and reloads that the compiler inserts; no hand scheduling
was attempted and the count is what clang emitted.

**Per-batch overhead.** The inversion makes 255 calls of fe_sqr_f and 15 of fe_mul_f:
at most 255·(82 + 1) + 15·(115 + 1) = 22,905 instructions, plus the non-call
instructions of fe_inv (397) and at most 4 per squaring for argument set-up and loop
control (1,020): **I ≤ 24,322**. With the prologue, first/last element and epilogue of
the two passes (1,026 − 622 + 173 − 128 = 449) and a 64-instruction allowance for the
call and the return, a batch of at least 2^15 walks has fixed work of at most 24,835 <
2^15 instructions, **under one instruction per walk**. Batches have 2^15..2^16 walks.

### 5.3 Per ordinary step

| work | instructions |
| --- | ---: |
| counted: pass1 iteration | 128 |
| counted: pass2 iteration | 622 |
| 1.1 F1: 8 MUL/UMULH, 7 ADD, 12-load allowance | 27 |
| 1.1 range test, F0 load, S, wrap test, F2 | 12 |
| 1.1 decision: differences, gap tests, argmin, fallback branch | 14 |
| sign, p − y, select limbs of K | 13 |
| h, distinguished test, twist index, table address | 7 |
| arrival pass: 4 loads of y, store of the table pointer, 3 header words | 8 |
| fingerprint tests, history and checkpoint stores | 26 |
| zero-denominator test | 6 |
| length, cap test, meter additions, loop | 14 |
| share of the inversion and of the call overhead (< 1 each) | 2 |
| **total** | **879** |

The step is charged **890 instructions = 0.890 units**, eleven above the table. Rows
below the counted pair are the analytic rows of 3682292, which are unchanged in kind:
each is an upper bound on an explicit instruction sequence, and the sequences total 114
(24 + 8 + 13 + 13 + 7 + 5 + 26 + 4 + 13 + 1) against the 129 allowed, so the derived step
is 864 instructions. The rows that 3682292 spent on loading the table point and the walk
point and storing the new walk point (24) and on the batch inversion share (4) are now
inside the counted loops and the per-batch overhead, which leaves only what the arrival
pass still has to read and write: the y limbs, the table pointer it produces and a few
header words. The analytic sequences are:

- F1 (24): 4 loads of x and 5 of G_1..G_5 (9), 4 MUL and 4 UMULH (8), 7 ADD (7).
- Range test, F0, S, wrap test, F2 (8): SUB, CMP, branch (3); S = F0 + F1 by ADDS (1);
  NEG (1); CSEL on the carry, CMP with 40, branch (3).
- Decision (13): CMP, 2 CSEL, CSET (4); CMP, CSEL (2); CMP and 3 CSEL (4); difference,
  CMP with 40, branch (3). A gap above 40 to the second smallest gives a gap above 40
  to both others.
- Sign, p − y, K (13): p − y by SUBS + 3 SBCS (4); p − y < y by SUBS + 3 SBCS (4);
  CSET (1); 4 CSEL (4).
- h, distinguished test, t, address (7): AND (1); CBZ on limb 1 of K (1); t = 3s + j by
  two ADD (2); address = base + 64·(6h + t) by three ADD with shifts (3).
- Arrival pass traffic (5): 4 loads of the y limbs of the new point (the x limbs are the
  4 loads of F1); 1 store of the table pointer. The table point and walk point loads and
  the walk store are in the counted loops. 8 allowed.
- Fingerprints (26): pending-escape flag, load and CBNZ (2); two-slot ring indexed by
  length mod 2, slot address by AND and ADD (2); compare fp(K) with K_{−2}: 3 loads, CMP
  and 2 CCMP, branch (7); overwrite the slot with fp(K), 3 stores (3); checkpoint compare:
  3 loads, CMP and 2 CCMP, branch (7); checkpoint write when length mod 64 = 0: TST,
  branch (2), 3 stores (3).
- Zero denominator (4 of 6): in pass 1, ORR of the four limbs of dx in pairs and
  together, CBZ to a cold path that replaces dx by 1, abandons the trail and is priced as
  a separate action (rate < 2^-114 in total). Because dx is canonical, dx = 0 exactly
  when x(T_i) = x(W_i), which also covers W_i = −T_i.
- Length, cap, meter, loop (13): the 2-limb length (L_max > 2^64): 2 loads, ADDS + ADC,
  2 stores (6); cap test CMP of the high limb with 20, branch (2); meter and operation
  count: SUBS on the processor's block remainder, ADD on its operation count, branch to
  the block refill (3); the loop counter SUBS and branch of the arrival pass (2). The
  loops of pass 1 and pass 2 carry their own counters inside the counted 128 and 622.
- Share (1 of 2): fixed work per batch under 2^15 instructions, section 5.2.

**Other actions** keep the allowances of 3682292 and are priced as before. They are used
only outside the hot loop and are loose against the compiled kernels:

- A = 32 (load 8, ADDS + 3 ADCS, carry, SUBS + 4 SBCS, 4 CSEL, store 4, ≤ 6 overhead).
- M ≤ 227. The compiled fe_mul_f is 115 instructions including its frame; the older
  allowance is retained for the rare paths.
- Doubling: 7M + 6A + 147 ≤ 1,928 instructions at these allowances; priced 2.0. A
  doubling needs 5 multiplications and 2 squarings, so the compiled kernels give far less.

**Prices** (units):

| action | price | basis |
|---|---|---|
| ordinary step | 0.890 | 890 instructions: 750 counted + 129 analytic rows = 879 (section 5.3) |
| traversal step | 0.920 | ≤ 24 more instructions to copy the remembered member: ≤ 903 |
| doubling | 2.0 | 7M + 6A + 147 ≤ 1,928 instructions |
| §1.1 fallback | 0.5 | one M, two A, compares: < 400 instructions |
| random 64-bit word | 0.01 | |
| sampling attempt | 0.05 | |
| trail start [a]P + [b]Q | 2^13 | two Jacobian scalar multiplications at ≤ 512 operations × (16M + 16A) each, an addition and an inversion: < 4.8·10^6 instructions |
| report | 40 | route, message, node-side probes and insert: < 2^15 instructions |
| store probe | 0.05 | |
| store initialisation | 2·D_cap slots × 16 instructions = 2^58.9, once | |
| table | 2048 starts plus 5 φ-twists each: < 2^24.1, plus sampling | |
| phase-2 or re-walk operation | 600 | unbatched: an inversion (< 2^17) plus 7M, coefficient updates mod r and control, far below 600,000 instructions |
| small-batch inversion | 600 | one inversion, ≤ 24,322 instructions = 25 units; 600 is a loose upper bound |
| solving | 600 | |
| verification | 2^13 | |

**Budget.**

- **Phase-1 operations:** 0.890·(S_ops + 2^60) < 2^126.6138.
- **Hard extras:**

  | item | units |
  |---|---|
  | meter remainders | < 2^61 |
  | phase 2: 2^20·(L_max + 65)·600 | 2^97.55 |
  | re-walks: 8·2·((L_max + 65)·600 + 2^13) | 2^81.6 |
  | reports: D_cap·40 | 2^69.2 |
  | store initialisation | 2^58.9 |
  | table and verification | < 2^25 |
  | **total** | **< 2^97.56** |

- **Random extras (expected):**

  | item | units |
  |---|---|
  | doubling surcharges: (2.0 − 0.890)·2^-13.56·(S_ops + 2^60) | < 2^113.38 |
  | traversal surcharges: 0.03·64·2^-19.99·(S_ops + 2^60) | < 2^107.74 |
  | starts: < 2^62.8·2^13 | 2^75.9 |
  | words, attempts, probes, small batches, fallbacks | < 2^70 |
  | **total** | **< 2^113.41** |

- **Headroom:** 2^126.65 − 2^126.6138 − 2^97.56 > 2^121.31.
- **Early meter halt:** probability < 2^113.41/2^121.31 < 0.0042 (Markov).

`experiment.py` with `"kernels": true` recomputes the whole ledger to 50 digits from
r, R, the prices, D_cap and L_max: S_0 = 2^126.71151, S_ops = 2^126.78190, phase 1
2^126.61377, hard extras 2^97.5508, expected random extras 2^113.4011, headroom
2^121.3163, Markov 0.004143, success floor 0.390857.

**Preprocessing:** table < 2^24.1, node set-up 2^26, store initialisation
2^58.9, and the first 2^20 starts 2^33. The expected total is < 2^58.91, and the cap is 2^59.0. The
random sampling part is expected < 2^25, far below the gap, so the cap halts a
run with probability < 2^-30 by Markov. All of it is inside W.

**W ≤ 2^126.65 on every run.**

**Memory:** the store holds 2·D_cap slots of 104 bytes at load 1/2, so
2^64.83 slots × 104 bytes = 2^71.53 bytes. Adding the table (1.5 MB) and the
2^20 × 512 bytes of walk states gives < 2^72 bytes, with each of
the 2^16 nodes < 2^56 bytes. The batch scratch (acc: 32 bytes and a table pointer: 8 bytes
per walk) adds at most 40 bytes per walk, under 2^26 bytes for 2^20 walks.

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
adds a 1.05·S0 budget, phase 2, report
capacity, the 8-match cap, and routing to 16 nodes of capacity ⌈2·D_cap/16⌉
with overflow treated as failure. [a]G + [b]T is evaluated as [a + bx]G for
speed.

Offline results (seeds 1501–1512; every recovered x verified; S_ops = 1.05·√(2N' ln(1/0.605)) for each toy; T_c counts all arrivals):

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

**Long-trail trap test** (q29, walks without distinguished points, 150 trails
of 2,048 operations each, against the additive-escape control):

| R | escape | escapes per trail | trails with non-fruitless returns |
| ---: | --- | ---: | ---: |
| 4 | doubling | 72 | 0 / 150 |
| 8 | doubling | 42 | 6 / 150 |
| 4 | additive (control) | 223 | 110 / 150 |
| 8 | additive (control) | 64 | 8 / 150 |

Ordinary self-collisions predict non-fruitless returns in about 1 − exp(−L²/(2N)) = 6.8% of trails (L²/(2N) = 0.0703). At R = 8 the additive control traps less here than in 1c600df4's BLS12-381 run (8/150 against 26/150). At R = 4 it traps 110/150, while the doubling escape gives 0/150 and 6/150.

Doubling escapes behave like ordinary self-collisions; the additive control
traps trails structurally.

**Doubling rate against the toys.** The toy 2-cycle rates at R = 2048
(8.23·10^-5 and 8.28·10^-5 per step in the T_c rows, 7.86·10^-5 and
7.72·10^-5 in the pipeline rows) agree with 1/(6R) = 8.14·10^-5 within
about two standard errors; each row has 10^3 to 2·10^3 events, a relative
error of 2% to 3%. The toys also show
checkpoint escapes at about 5·10^-6 per step, above the structural rate
2^-19.99. These are ordinary self-collisions of short toy trails: a trail of
length ℓ meets its own past with probability about ℓ/N per step, about
32/(5.6·10^6) ≈ 6·10^-6 at q25 with dp_bits = 5. In production
ℓ/N < 2^70/2^253, so the structural rate governs. In the q25, R = 2048 row
(seed 1502) the total doubling rate, both effects included, is
8.72·10^-5 per step: 6% above the production bound 2^-13.5681 = 8.23·10^-5,
within the production model's stochastic variation. This small toy sample
does not certify the production instruction budget.

**Provenance of the rows.** `experiment.py` is that of 3682292 (which changed one toy setup
check from `mul(r,G)`, which first reduced r modulo r, to `mul_raw(r,G)`) plus the optional
`kernels` mode of section 5. The walk and the seeded samples of a valid toy generator are
unchanged. The offline rows above are the earlier packages' results; they were not rerun for
this revision. Python's seeded generator makes each row reproducible.

**Initial input.** The hosted configuration (q25, R = 2048, dp_bits = 5, seed 1, pipeline,
kernels, 6 s) is time-limited, so the number of completed trials depends on the machine. A
local run of exactly this configuration completed 464 trials with pipeline success 0.3944
(binomial standard error 0.023), 80 two-cycle escapes against 78.4 expected, no node
overflow, all part A predicates true (30,000 random and 10,000 stress index-rule checks, 0
mismatches), and the kernels check of section 5 with 0 mismatches and the ledger above.
This short run is compatible with the larger offline row, but it does not prove the
production H1 floor, and it does not measure production instruction counts. The H1 floor
rests on the disclosed heuristic and the larger toy evidence.

**Limitations.** The instruction counts are those of the shipped clang output for the
loops of the batch round, plus analytic rows for the arrival bookkeeping; they are not a
Callgrind measurement of GCC output and not a measurement of a full attack. The 129 analytic
rows are carried from 3682292 and rest on explicit sequences, not on compiled code. The
shipped AArch64 Linux assembly could not be executed on the Darwin machine that tested the
same C source natively. H1 is a heuristic: a random mapping with the Brent–Pollard
correction; van Oorschot–Wiener; automorphism classes as in Wiener–Zuccherato,
GLV and Duursma–Gaudry–Morain; and fruitless-cycle escapes as in
Bos–Kleinjung–Lenstra and Bernstein–Lange–Schwabe. It is supported only at toy
scale. Instruction counts are analytic. This is a generic attack.

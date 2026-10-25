# Six-class secp256k1 rho at R = 2^20: a 678-instruction step (126.1539)

**Claim.** For a uniform full-production secp256k1 input `Q=[x]P`, the classical cold six-class
(negation + GLV) parallel Pollard-rho attack of `upstream-proof.md`, with the counted kernels and
bookkeeping rows of this document, recovers `x` with probability greater than `0.39001` under the
explicit H1 heuristic. A global meter charges every action before it runs, so on every run its total
work is at most

```text
time_log2 = 126.1539 = 1261539/10000
```

target-group-operation equivalents (C = 1000 scalar AArch64 instructions per unit). Peak memory is
below `2^72` bytes; preprocessing is capped at `2^59` units, charged and included; nonuniform advice is
zero. This is an analytic heuristic bound, not a practical break.

## 1. What changes

This package is this solver's accepted 5d3bf6c4 (126.1913) with the walk's table enlarged from R = 2048 to
R = 2^20 and cycles detected by the checkpoint alone (section 5.1). The two-cycle test goes, with its history
ring, the pending-escape flag and the copies of K_cp, and the fingerprint row drops from 18 to 8 instructions
(section 3). With cycles that rare the expected random charges fall from 2^110.6 to 2^101.8 units, so the
Markov headroom h and the floor p0 are retuned (section 5.3). The counted kernels and loops of section 2 are
unchanged, byte for byte. The class key, canonical index, doubling escape, phases, meter, report store, every
other hard row and the H1 calibration (`kappa = 1.01`) are those of 02feb284 (126.3059) and its public
ancestors (`upstream-proof.md`, sections 1 to 4), except for the walk changes of section 5.1. The accepted
BLS12-381 package 67154fda runs the same six-class walk with R = 2^20 and checkpoint-only detection.

| item | 02feb284 | 5d3bf6c4 | this package |
| --- | ---: | ---: | ---: |
| table size R | 2048 | 2048 | **2^20** |
| cycle detection | 2-cycle test + checkpoint | 2-cycle test + checkpoint | **checkpoint only** |
| counted loops, pass 1 / pass 2 | 111 / 521 | 108 / 493 | 108 / 493 |
| analytic rows | 112 | 87 | **77** |
| ordinary step | 744 | 688 | **678** |
| doubling / traversal | 950 / 768 | 838 / 712 | **828 / 702** |
| Markov headroom h | 0.003866 | 0.003866 | **0.000165** |
| p0 / success floor | 0.3939163 / 0.3900501 | 0.39388 / 0.3900138 | **0.3901753 / 0.3900101** |
| claim | 126.3059 | 126.1913 | **126.1539** |

The ordinary step is charged exactly its table, with no margin: **108 + 493 + 77 = 678 instructions
= 0.678 units**.

## 2. The counted kernels

`kernels/gen_loops.py` emits the two loops as straight-line AArch64 (`kernels/secp_loops.s`) and, from the
same emitters, the test wrappers `t_mul`, `t_sqr`, `t_sub`, `t_subcanon` and `t_reduce`
(`kernels/secp_tests.s`). p = 2^256 − C0 with C0 = 2^32 + 977, so 2^256 ≡ C0 (mod p). A field element is
four little-endian 64-bit limbs; **weak** means any value in [0, 2^256), **canonical** means [0, p). Walk
points and table points are canonical. x4 holds C0 inside the loops.

### 2.1 fe_mul: 32 products and six carry chains (62)

Write `(Hij, Lij)` for the high and low words of `a_i·b_j`. `Lij` sits in column `i+j` and `Hij` in column
`i+j+1` (column k has weight 2^(64k)). The columns 0 to 7 hold 1, 3, 5, 7, 7, 5, 3, 1 words. Reducing a column
of t words to one word takes t − 1 two-word additions, 24 in all. An ADDS/ADCS chain makes one addition per
column and hands its carry to the next column inside the same chain; MUL and UMULH between its instructions
leave the flags alone. Every chain therefore ends with one addition of its last carry. Columns 3 and 4 need
six additions each, so at least six chains are needed; the row-by-row schedule of the public package used
seven (31 additions), this one uses six (30):

| chain | col 1 | col 2 | col 3 | col 4 | col 5 | col 6 | col 7 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | H00+L01 | H01+L02 | H02+L03 | H03+L13 | H13+c | | |
| 2 | +L10 = z1 | +H10 | +H20 | +H30 | H22+c | | |
| 3 | | +L20 | +L30 | +L22 | H13'+H22' | H23+c | |
| 4 | | +L11 = z2 | +H11 | +H12 | +H31 | H32+c | |
| 5 | | | +L12 | +H21 | +L23 | H23'+H32' | H33+c |
| 6 | | | +L21 = z3 | +L31 = z4 | +L32 = z5 | +L33 = z6 | +c = z7 |

z0 = L00. "+X" adds word X to the column's running sum; "H+c" adds only the incoming carry to a high word.
The chains run in the order 1 to 6, and every word is used exactly once (the experiment's port asserts it).

- **No carry is lost.** A high word is at most ⌊(2^64 − 1)^2 / 2^64⌋ = 2^64 − 2, and chains 1 to 5 end by
  adding one carry bit to a high word no earlier instruction has touched (H13, H22, H23, H32, H33), so
  those ADCs cannot carry out. Chain 6 ends by forming z7 after z0 to z6 are final; the eight limbs then
  satisfy Σ z_k·2^(64k) = a·b < 2^512, so z7 < 2^64. Every other addition's carry is consumed by the next
  instruction of its chain. Hence (z0, …, z7) = a·b exactly.
- Each product is formed just before the addition that consumes it, and an operand limb is released after
  its last product, which keeps the register peak of fe_mul at 15 including its 8 operand limbs.

### 2.2 The reduction (26)

z = L + H·2^256 with L, H < 2^256 (the eight product limbs).

1. **H·C0** (8 + 4): four MUL and four UMULH give g_j, k_j with g_j + 2^64·k_j = z_(4+j)·C0, k_j ≤ C0 − 1.
   ADDS g1+k0, ADCS g2+k1, ADCS g3+k2, ADC k3+c gives H·C0 in five limbs; H·C0 < 2^256·C0, so the ADC
   cannot overflow.
2. **s = L + H·C0** (5): ADDS, 3 ADCS, then ADC into the fifth limb s4. Since
   s ≤ (2^256 − 1)(1 + C0) < (1 + C0)·2^256, s4 ≤ C0.
3. **w = s4·C0** (2): MUL and UMULH; w ≤ C0^2 < 2^65.
4. **t = (s mod 2^256) + w** (4): ADDS, 3 ADCS, with carry-out c′. t ≡ z (mod p) because 2^256 ≡ C0.
5. **Final fold** (3). If c′ = 1, the 256-bit result t′ = t − 2^256 is below w < 2^65 (because
   s mod 2^256 < 2^256), so its limbs 2 and 3 are zero and limb 1 is at most 1. Then t′ + C0 < 2^65 + 2^33
   is formed by adding C0 to limb 0 and the carry to limb 1: CSEL (C0 or 0 on c′), ADDS, ADC. If c′ = 0 the
   same three instructions add zero.

The result is below 2^256 and ≡ z (mod p). **26 = 8 + 4 + 5 + 2 + 4 + 3.** The public package added the final
C0 across all four limbs (28). The branch c′ = 1 needs s mod 2^256 ≥ 2^256 − w, about 2^-190 for random
inputs, so the experiment constructs such inputs directly (section 2.6): for a chosen high half H it solves
for the low half L that puts s mod 2^256 inside [2^256 − s4·C0, 2^256), including the boundary
t′ = 0.

### 2.3 fe_sqr (67)

- **Off-diagonal sum** O = Σ_(i<j) a_i·a_j·2^(64(i+j)) (12 + 8). Six products. Columns 1 to 6 hold 1, 2, 3, 3,
  2, 1 words, so six two-word additions in two chains: chain A adds L03+L12 (column 3), H03+H12 (column 4)
  and its carry into H13; chain B adds H01+L02, then H02, L13, the incremented H13 with L23, and its carry
  into H23, which is o6. As in fe_mul, the carry-only additions land on untouched high words, and o6 is the
  top limb of O < 2^448. The public package used 9 additions (three chains).
- **Doubling** (6). 2·O ≤ a^2 < 2^512, so no bit is lost: LSR d7 = o6 >> 63 and five EXTR #63 give d2 to d6.
  The low limb d1 = 2·o1 mod 2^64 is not formed: the first diagonal addition is ADDS z1, H00, o1, LSL #1,
  and the bit shifted out of o1 is bit 0 of d2, already taken by the EXTR.
- **Diagonals** (8 + 7). a_i^2 by MUL and UMULH, then ADDS z1, five ADCS for z2 to z6, and ADC
  z7 = d7 + H33 + c, the top limb of a^2 < 2^512, which cannot overflow.
- **Reduction** (26).

**67 = 20 + 8 + 6 + 7 + 26.**

### 2.4 fe_sub (9) and fe_subcanon (14)

**fe_sub(a, b)**, a weak, b canonical, is unchanged: SUBS, 3 SBCS; CSEL of C0 on borrow; SUBS, 2 SBCS, SBC.
On borrow the wrapped value is at least C0 + 1, so subtracting C0 cannot borrow again; the result is weak,
≡ a − b, and canonical whenever a is.

**fe_subcanon(a, b)**, a weak, b canonical, returns the canonical (a − b) mod p:

```text
SUBS, 3 SBCS      d = a - b mod 2^256, borrow beta             (4)
SBC  m, xzr, xzr  m = 2^64 - 1 if beta, else 0                 (1)
CSNEG k, C0, C0, hs   k = C0 if no borrow, else 2^64 - C0      (1)
ADDS, 3 ADCS      u = d + K,  K = (k, m, m, m) = p if beta, else C0; carry-out gamma   (4)
4 CSEL (cs)       result = u if gamma, else d                  (4)
```

If β = 1 then a < b ≤ p − 1, so d = a − b + 2^256 ≥ 2^256 − (p − 1) = C0 + 1 and d + p ≥ 2^256 + 1: γ = 1
and the result is d + p − 2^256 = a − b + p ∈ [1, p − 1]. If β = 0 then d = a − b and K = C0: γ = 1 exactly
when d ≥ p, and then the result is d − p < C0 < p; otherwise it is d < p. In every case the result is
canonical and ≡ a − b. The step uses it twice, `x3 = fe_subcanon(fe_sub(λ^2, x_W), x_T)` and
`y3 = fe_subcanon(λ·(x_W − x3), y_W)`, each replacing an fe_sub followed by fe_canon (9 + 8 = 17).

### 2.5 The two loops (108 and 493)

**Pass 1** (forward; the schedule of the public package with the new fe_mul): LDR of the table pointer and two
LDP each for x(T_i) and x(W_i) (5); dx_i = fe_sub(x(T_i), x(W_i)) (9, canonical); two STP each for acc_(i−1)
and dx_i into slot i (4); acc_i = dx_i·acc_(i−1) with the reduction writing into the four loop-carried
registers (88); SUBS + B.NE (2). **108.**

**Pass 2** (backward, i = n−1 … 1). R ≡ 1/(dx_0···dx_i) stays in x25–x28 across iterations: it is loaded once
from Rinv before the loop and never stored (the public package kept it in a stack slot, two LDP and two STP
per walk). Per iteration:

| work | instructions |
| --- | ---: |
| LDR table pointer; 2 LDP acc_(i−1); inv = R·acc_(i−1) | 1 + 2 + 88 |
| 2 LDP dx_i; R ← dx_i·R (result into x25–x28) | 2 + 88 |
| 2 LDP y(T_i), 2 LDP y(W_i) (kept to the end); dy = fe_sub(y_T, y_W) | 4 + 9 |
| λ = dy·inv; λ^2 | 88 + 67 |
| 2 LDP x(W_i) (kept); fe_sub(λ^2, x_W); 2 LDP x(T_i); x3 = fe_subcanon(·, x_T) | 2 + 9 + 2 + 14 |
| x_W − x3 = fe_sub(x_W, x3); 2 STP x3 | 9 + 2 |
| λ·(x_W − x3); y3 = fe_subcanon(·, y_W); 2 STP y3 | 88 + 14 + 2 |
| SUBS + B.NE | 2 |
| **total**: 4·88 + 67 + 3·9 + 2·14 + 13 loads + 4 stores + 2 | **493** |

Walk 0 uses R itself as its inverse after the loop. This is the standard affine addition
λ = (y_T − y_W)/(x_T − x_W), x3 = λ^2 − x_W − x_T, y3 = λ(x_W − x3) − y_W, with canonical outputs.

**Registers.** The pool is x5–x17 and x19–x30 (25 registers). x29 and x30 are saved with x19–x28 by the
prologue and restored before RET; the loops call nothing. x0–x3 hold the count and the three array pointers, x4
holds C0, x18 is never written. The register peak is 25 in pass 2 and 16 in pass 1; there are no spills and
no stack access inside either loop. Both loops are single basic blocks: every conditional is a carry chain,
CSEL or CSNEG, and the only branch is the back edge. Opcode histograms:

- pass 1: LDR 1, LDP 4, STP 4, MUL 21, UMULH 21, ADDS 10, ADCS 26, ADC 9, SUBS 3, SBCS 5, SBC 1, CSEL 2,
  B.NE 1 (108);
- pass 2: LDR 1, LDP 12, STP 4, MUL 99, UMULH 99, ADDS 49, ADCS 127, ADC 42, SUBS 9, SBCS 21, SBC 5,
  CSEL 16, CSNEG 2, EXTR 5, LSR 1, B.NE 1 (493).

Counting is the convention of the public package: one LDP or STP is one instruction.

### 2.6 Evidence

`experiment.py` with `"kernels": true` (the declared initial run) does all of the following:

1. **Limb ports.** Python ports of the generated fe_mul (the six chains, as data), the reduction, fe_sqr,
   fe_sub and fe_subcanon assert every no-overflow step above. Those steps are each carry-only ADC, the
   exact product identity, s4 ≤ C0, w1 ≤ 1, and, when the final fold fires, limbs 2 and 3 zero and limb 1 ≤ 1.
   They also assert that fe_sub never borrows twice and that fe_subcanon's d + p always carries after a borrow.
   They are checked against big integers on all pairs of 44 edge values (0, 1, p − 1, p, p + 1, 2^256 − 1,
   values near p and near 2^256, every all-ones/all-zeros limb pattern) and on random values. The
   reductions include 200 constructed inputs whose final fold fires. A 24-walk batch built from the
   ports is checked against affine secp256k1 addition.
2. **The shipped loops, executed.** The AArch64 interpreter executes `kernels/secp_loops.s` on production
   points for batches of 2, 3, 5, 9 and 17 walks, compares every output with affine arithmetic, and reports
   the executed-instruction difference per extra walk. It must be exactly 108 (pass 1) and 493 (pass 2).
3. **The emitters, executed.** The interpreter executes `kernels/secp_tests.s`, the same emitters wrapped as
   functions, on edge and random inputs: 730 products, 444 squares, 818 fe_sub and fe_subcanon pairs, and 402
   reductions, at least 200 of which carry out of the second fold. It compares every result with big integers.
4. **Counts and ledger.** It parses the instruction counts from the shipped files, checks the identities
   108 = 88 + 9 + 5 + 4 + 2 and 493 = 4·88 + 67 + 3·9 + 2·14 + 13 + 4 + 2, and recomputes the doubling bound,
   the per-batch overhead and the ledger of section 5.

The interpreter of the public package gains two standard A64 forms used by the new kernels: **CSNEG**
(xd = cond ? xn : −xm mod 2^64) and **ADDS with an LSL-shifted register** (xn + (xm << k mod 2^64), with the
carry of that 64-bit sum). All other forms are unchanged.

Before submission, the generated files were also assembled with clang 23 for `aarch64-linux-gnu` and the
object code was run in the Unicorn 2.1.4 AArch64 CPU emulator, off the hosted runner and not shipped. It
covered 4,572 walks in random batches of 2 to 199, 2,904 products, 648 squares, 10,204 fe_sub and fe_subcanon
pairs and 1,345 reductions (600 with the final fold). There were 0 mismatches, counts of exactly 108 and 493
per walk, and x19–x30 restored on return.

## 3. Analytic rows: 112 → 77

The rows cover what the arrival pass does per class arrival outside the two counted loops (pseudocode of
`upstream-proof.md` section 2, with the walk changes of section 5.1). Four rows are unchanged from 02feb284:
range test 8, decision 13, arrival traffic 3, inversion and call share 1. Six change:

| row | 02feb284 | 5d3bf6c4 | this | sequence |
| --- | ---: | ---: | ---: | --- |
| F1 truncated product | 20 | 17 | **17** | 2 LDP of the x limbs; 4 MUL and 4 UMULH against G_1..G_5 held in registers; 7 ADD |
| sign, p − y, K | 17 | 12 | **12** | p − y as SUBS from a register holding p0 = 2^64 − C0, then 3 SBCS from a register holding 2^64 − 1 (4); compare y with p − y by SUBS, 3 SBCS (4); 4 CSEL for K (4); s is not set in a register, the table address reads it from the same flags |
| h, DP test, table address | 7 | 6 | **6** | AND of K's low limb with the 20-bit mask for h (1); CBZ on limb 1 of K (1); ADD base + (h << 9), ADD + (j << 6), ADD #256 (3); CSEL on the carry of the y-against-(p − y) compare (1) |
| fingerprints, history, checkpoint | 26 | 18 | **8** | LDP + LDR of K_cp (2); CMP, 2 CCMP, B.EQ against fp(K) (4); checkpoint test TST of the low length limb, branch (2); its STP + STR are the separate checkpoint-write action. No two-cycle test, history ring or pending flag (section 5.1) |
| length, cap, meter, loop | 13 | 9 | **9** | LDR of the low length limb, ADDS #1, STR, B.CS to the carry path (4); meter (3); arrival loop (2) |
| zero-denominator test | 4 | 0 | **0** | one test per batch (section 4) |
| **total** | **112** | **87** | **77** | |

- **p − y.** p's limbs are (p0, 2^64 − 1, 2^64 − 1, 2^64 − 1). With p0 and 2^64 − 1 held in two registers for
  the whole arrival pass (set once per batch, 3 instructions, section 4), p − y is one SUBS/SBCS chain.
  The old row negated y from zero and then subtracted C0 (8 instead of 4). y is canonical and nonzero (no
  point of order 2), so p − y ∈ [1, p − 1] with no borrow out. y ≠ p − y because p is odd, so the carry of
  y − (p − y) is s = [p − y < y], and K = min(y, p − y) is selected by four CSEL.
- **Table rows of eight entries.** Tab[h][t] with t = 3s + j is laid out as Tab[h][s][j] at
  base + 512h + 256s + 64j: each h owns 512 bytes, two halves of four 64-byte entries, of which j uses three.
  The address is AND, three ADDs and one CSEL on the carry of the y-against-(p − y) comparison, and that
  carry is exactly s. CSEL, AND, CBZ and ADD leave the flags alone, so nothing between that comparison and the
  CSEL changes them. With R = 2^20 the table holds 2^20·512 bytes = 512 MiB; the padding entries are never read.
- **G_1..G_5 in registers.** F1 multiplies the four x limbs by the five limbs G_1..G_5 of G = ⌊2^384·β/p⌋
  (`upstream-proof.md` section 1.1). The arrival pass keeps them in five registers for the whole pass, loaded
  once per batch (LDP, LDP, LDR: 3 instructions in the per-batch overhead of section 4), so F1 loads only the
  x limbs. Register map of the arrival pass: x0 loop counter; x1–x4 walk-slot, pointer-array, trail-record
  and table pointers; x5, x6 the processor's block remainder and operation count; x7–x11 G_1..G_5; x12 p0;
  x13 2^64 − 1; x18 unused; x29 and x30 saved by the prologue. That leaves x14–x17 and x19–x30 (16 registers)
  for one walk's values. The most live at once is y, p − y and j while K is formed (9 words); the
  checkpoint compare holds K's three fingerprint words, three loaded words and the table address.
- **Length.** Only the low limb is updated per arrival: ADDS #1 and STR, then B.CS to a carry path. The
  carry path runs only when the low limb wraps. It loads, increments and stores the high limb, compares it
  with 20 and abandons the trail at 20, at most 6 instructions. The length grows by one per operation, so it
  first reaches L_max = 20·2^64 exactly at such a carry, and testing the cap there alone is exact. A trail
  wraps at most once per 2^64 of its operations, so a run has fewer than (W/0.678)/2^64 < 2^63 carries. The
  carry path is metered as its own action at 0.01 units; its expected total is below 2^57 units, inside the
  2^70 "other" random row. The low limb still gives length mod 64 (checkpoint).
- **Checkpoint.** With cycles caught by the checkpoint alone (section 5.1), the row loads K_cp's three words
  (LDP + LDR), compares them with fp(K) (CMP, 2 CCMP, B.EQ) and tests length mod 64 (TST, branch). Clearing
  the history stores the marker 2^64 − 1 in K_cp's top word; no fingerprint has it, because K ≤ (p − 1)/2 makes
  limb 3 of K below 2^63. The write K_cp ← fp(K) runs only when length mod 64 = 0. Its STP + STR are metered as
  their own action at 0.002 units when they run. Apart from a trail's first write (covered by the 2^13-unit
  start price), a trail writes at most once per 64 of its operations. Phase 1 performs fewer than
  S_ops + 2^61 operations, so section 5.3 charges ⌈(S_ops + 2^61)/64⌉ writes at 2 instructions each as a hard
  row. Phase-2 and re-walk operations are priced at 600 units each and cover their own writes.
- **Zero denominator.** See section 4.

## 4. Prices and per-batch overhead

- **Ordinary step:** 108 + 493 + 77 = **678** instructions = 0.678 units.
- **Traversal step:** 678 + 24 = **702** (the 24 instructions to track and resume the minimum-fp member are
  unchanged).
- **Doubling:** the batched doubling shares the product chain and inversion. Its arithmetic is
  5 fe_mul (440) + 2 fe_sqr (134) + dx = 2y by one weak addition (9) + 3x^2 by two weak additions with a
  second carry fold (2·14) + x3 = fe_subcanon(fe_sub(λ^2, x), x) (9 + 14) + x − x3 (9) +
  y3 = fe_subcanon(λ·(x − x3), y) (14) = 657.
  - **Traffic: 30** (the public package allowed 72). The loads, stores, pointer updates and loop control are
    at most the addition loops' 5 + 4 + 2 and 13 + 4 + 2. A doubling walk reads its own y for dx = 2y in the
    forward pass and acc, dx, x and y in the backward pass, with no table point. It stores the same dx, acc,
    x3 and y3, and shares R in registers.
  - **List moves: 34** (the public package counted 50). Moving the walk to and from the doubling list is
    the upstream 16 loads and 16 stores of 64-bit words, the walk's x and y out and back. They pair into
    8 LDP and 8 STP, plus 18 other instructions.
  - **Escape logic: 30**, as upstream (with checkpoint-only detection it is now only the clearing of the
    history, one STR of the marker into K_cp, well inside the 30).

  The doubling is 657 + 30 + 77 + 34 + 30 = **828**. The weak additions, the 18 and the 30 are those of
  `upstream-proof.md` section 5.3.
- **Zero-denominator test, once per batch.** Pass 1 multiplies the n canonical dx_i into acc_(n−1). dx_i ≡ 0
  only if dx_i = 0, and p is prime, so acc_(n−1) ≡ 0 (mod p) exactly when some x(T_i) = x(W_i). Before the
  inversion, acc_(n−1) is canonicalised and tested: ADDS, 3 ADCS, 4 CSEL, 3 ORR, CBZ, 12 instructions per
  batch. Only when it fires does the batch take the cold path, an event of probability < 2^-114 over the whole
  run as before. Pass 1 is re-run with the per-walk test, each walk with dx_i = 0 is flagged abandoned, and its
  dx_i is replaced by 1 in the product and its slot, exactly as in the public package. The re-run, at most
  2^16·(108 + 4) < 2^23 instructions, is metered as its own action at 2^13 units and falls in the "other"
  random row.
- **Per-batch fixed work:** inversion I ≤ 24,322 (unchanged `secp_fe.s`, `upstream-proof.md` section 5.1–5.2),
  prologue, walk 0 and epilogue of the two passes (141 − 108 + 827 − 493 = 367), 64 for the call and
  return, 12 for the zero test, 3 for the constants p0 and 2^64 − 1 and 3 for loading G_1..G_5:
  24,771 < 2^15. Batches have at least 2^15 walks, so the share row stays 1.

Every other price (reports 40.4, probes 0.05, starts 2^13, fallbacks 0.5, phase-2 and re-walk operations 600,
refills 0.03, words, attempts, store initialisation) is unchanged; the length carry (0.01) and the checkpoint
write (0.002) are their own actions. The table now has 2^20 entries (section 5.1).

## 5. The walk at R = 2^20, H1 and the exact budget

### 5.1 The walk change

The walk is that of `upstream-proof.md` section 2 with two changes: the table has R = 2^20 rows instead of
2048, and cycles are detected by the checkpoint alone. Everything else is unchanged: the class key |y|, the
canonical index, the twisted table Tab[h][t] = [u_t^{-1}]M_h with h = K mod R, the distinguished classes
(limb 1 of K zero, θ = 2^-64), M = 2^20 walks, L_max = 20·2^64 and the checkpoint period 64. At each class
arrival, in order:

```
if K distinguished: report (K, a, b); end trail
if length >= L_max: abandon trail
if fp(K) == K_cp: traverse the cycle once more (<= 64 steps, else abandon), resume at the member of its
     minimum-fp class remembered during the traversal, ESC; clear history; continue
if length mod 64 == 0: K_cp <- fp(K)
P <- P + Tab[h][t]
```

"Clear history" stores the marker of section 3 in K_cp. ESC is the doubling escape P ← [2]P, a class map.

- **Merging lemma.** Let trails A and B arrive at the same class C, with any histories. Until one of them
  escapes, both apply F and visit the same classes. A trail escapes only after K = K_cp. K_cp is set at an
  arrival, and the arrivals up to the next refresh are F-steps. So a repeat of K_cp is a cycle Z of F of length
  ≤ 64, and Z contains no distinguished class. Every escape from Z leaves from its minimum-fp class and goes to
  E(Z) = class(2c(min Z)), which depends only on Z. If A escapes from Z, B follows the same F path and enters
  Z. B sets K_cp inside Z within 64 arrivals, returns to it after |Z| ≤ 64 more, and escapes to the same E(Z)
  after a traversal of |Z| and one doubling. So A and B are again at a common class after at most
  64 + 64 + 64 + 1 = 193 extra operations of one of them. By induction they reach the same distinguished
  class, unless one of them first hits L_max or an exceptional addition (failure events). History is cleared
  at every doubling and starts have none, so K = K_cp always witnesses a real cycle of F.
- **Cycle rates.**
  - A 2-cycle needs the next class to use the same partition index with the inverse twist: rate 1/(6R) =
    2^-22.585 per step.
  - Other cycles of length 3 to 64 need every partition index used at least twice: rate at most
    Σ_L Σ_b S2(L, b)·R^-(L−b) = 2^-37.99999 at R = 2^20, the sum of `upstream-proof.md` section 4. That sum is
    2^-19.9954 at R = 2048, and `certificate.py` recomputes it exactly at both R.
  - Cycles of length 65..200 are rarer still (below 2^-206 per step already at R = 2048).
  - A cycle in which some partition index is used only once closes by chance: a fresh step lands on one of
    the trail's own t earlier classes, with probability at most t/N under H1's pseudo-random key fields. Over
    a trail that is at most L_max^2/(2N) = 2^-117.8 (the bound of H1's non-fruitless-return clause), under
    2^-54 in a whole run; their doublings and traversals cost under 2^-47 units, inside the 2^70 "other" row
    of section 5.3. On the toy group such cycles occur at about (trail age)/N per step and dominate the
    census's doubling counts at large R (section 5.2).
  - So the expected number of doublings is below (1/(6R) + 2^-37.99999)·(S_ops + 2^60) + 2^-54
    < 2^-22.5849·(S_ops + 2^60). Every 2-cycle is traversed in 2 steps and other short cycles in at most 64,
    so the expected number of traversal steps is below (2/(6R) + 64·2^-37.99999)·(S_ops + 2^60) + 64·2^-54
    < 2^-21.5839·(S_ops + 2^60). Only linearity of expectation over the pseudo-random partition and twist of
    each new class is used.
- **Detection delay.** A cycle is caught within 63 + |Z| ordinary steps of entering it. That is ordinary work:
  it is counted in S_ops, and its arrivals are counted in T_c. At R = 2^20 it is below 65·2^-22.585 + 128·2^-38,
  under 2^-16.5 of all operations.
- **Table.** 2^20 rows of 512 bytes (section 3), 512 MiB = 2^29 bytes. Even with one copy on each of the 2^16
  nodes (2^45 bytes in all), each node stays below 2^56 bytes and peak memory below 2^72 bytes
  (`upstream-proof.md` section 5.3, "Memory": the store is 2^71.64 bytes, 2^55.64 per node). Building it
  takes 2^20 draws of [c]P + [d]Q at ≤ 2^13 units each plus twists: below 2^33.1 units, charged in the hard
  extras (< 2^34) and in preprocessing, which stays below its 2^59 cap.
- **Instruction counts.** R enters the step only through the 20-bit mask of the AND that forms h and through
  the table base. The partition field becomes bits 0..19 of K, still disjoint from the distinguished bits
  64..127, the routing bits 128..143 and the bucket bits 144..193.
- **Precedent.** The accepted BLS12-381 G1 package 67154fda uses this walk at R = 2^20 with checkpoint-only
  detection. Its section on cycles makes the same rate and delay arguments.

### 5.2 H1 at R = 2^20

H1 keeps its form: P(T_c ≤ S_ops) ≥ 1 − exp(−(S_ops/κ)^2/(2N')) with κ = 1.01 and N' = N·R/(R − 1), now at
R = 2^20. T_c counts every class arrival, including the repeated arrivals inside fruitless cycles and during
the detection delay. Its other parts are as in `upstream-proof.md` section 4: pseudo-random key fields, the
non-fruitless return bound, and a uniform twist index independent of the partition index. κ = 1.01 was
calibrated by the public q29 census at R = 2048 (+0.70 % ± 0.12 %, 200,000 trials, eight tables). This package
adds a census of the walk as specified here (`census_cp.py`, `census-cp-q29.json`). It is an independent
implementation, run on q29 (N = 2^24.83) with dp_bits 5 and a fresh table for every 4,000 trials:

| R | mean T_c/√N (model) | slowdown | s.e.: trials / tables | P[T_c ≤ S0] | index-repeating cycles per op (model) | coincidental cycles per op (model) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2048 | 1.26336 (1.25362) | +0.78 % | 0.21 % / 0.23 % | 0.3896 | 7.92·10^-5 (8.23·10^-5) | 1.09·10^-6 (0.93·10^-6) |
| 2^16 | 1.25758 (1.25332) | +0.34 % | 0.21 % / 0.31 % | 0.3900 | 2.47·10^-6 (2.54·10^-6) | 1.33·10^-6 (0.93·10^-6) |
| 2^20 | 1.26017 (1.25331) | +0.55 % | 0.21 % / 0.19 % | 0.3896 | 1.80·10^-7 (1.59·10^-7) | 0.54·10^-6 (0.93·10^-6) |

Each row is 64,000 trials on 16 tables, about 4.27·10^8 operations.

- **Slowdown.** At R = 2^20 the slowdown is +0.55 % ± 0.21 %, 2.2 standard errors below κ − 1 = 1 %; at
  R = 2^16 it is +0.34 % ± 0.21 %, and the two pool to +0.44 % ± 0.15 %. The table-to-table spread is
  comparable to the per-trial error. P[T_c ≤ S0] is at or above H1's value at that budget,
  1 − 0.605^(1/1.01^2) = 0.3890 (binomial s.e. 0.0019). κ = 1.01 is kept.
- **R = 2048.** The checkpoint-only walk runs each 2-cycle until the next checkpoint catches it, about 37
  arrivals, all counted in T_c: about 0.3 % of the arrivals at R = 2048, under 2^-16.5 at R = 2^20 (section
  5.1). So the R = 2048 row has the largest slowdown, and the hosted toy, which keeps R = 2048 because it
  draws a fresh table per trial, is not the H1 evidence for R = 2^20.
- **Charged cycles.** Each detected cycle is classified by the partition indices of its steps. The
  index-repeating cycles, the ones the ledger charges, number 33,861, 1,051 and 77 against the model's
  1/(6R) + Σ = 35,221, 1,083 and 67.8. The model ignores that a 2-cycle through a distinguished class ends
  the trail instead (1 in 32 here); the R = 2^20 count is 1.1 Poisson standard errors above it.
- **Toy self-meetings.** The coincidental cycles (468, 568 and 232) are the self-meetings of
  `upstream-proof.md` section 6: a trail of age ℓ meets its own past at about ℓ/N per step. Their model,
  Σ min(length + 1, 64)/N, gives 0.93·10^-6 per operation; pooled over all 48 tables they occur at
  0.99·10^-6. Their per-table counts are heavy-tailed, because every trial of a chunk walks the same table
  and many trails fall into the same few short cycles (table-spread s.e. 0.44, 0.34 and 0.12·10^-6). They
  make the raw doubling rate at R = 2^20, 7.2·10^-7 per operation, 4.6 times 1/(6R); in production they are
  bounded by the L_max^2/(2N) clause (section 5.1).
- **Trap test.** The long-trail test of `upstream-proof.md` section 6 (q29, no distinguished points, 150
  trails of 2,048 operations, seed 7), rerun with this walk through `experiment.py`'s `long_trails` mode,
  finds non-fruitless returns in 0/150 (R = 4) and 4/150 (R = 8) trails with the doubling escape, against
  48/150 and 18/150 with the additive control. Ordinary self-collisions predict about 10/150.

### 5.3 The exact budget

H1 and κ = 1.01 remain the disclosed heuristic, supported only at toy scale. With the random charges this
small, the Markov headroom is h = 0.000165 and the floor p0 = 0.3901753, which keeps the success floor at
0.3900101 > 0.39001 > 0.39 with W below 2^126.1539.

```text
N'     = ((r-1)/6) * 2^20/(2^20 - 1)
S0     = sqrt(2*N'*ln(1/(1-0.3901753)))      = 2^126.6996628
S_ops  = 1.01*S0                             = 2^126.7140181
phase1 = 0.678*(S_int + 2^60)                = 2^126.153375
hard extras                                  = 2^111.806818
   checkpoint writes 0.002*ceil((S_int + 2^61)/64)        2^111.748234
   upstream rows (table < 2^34, refills at the new cap)    2^107.155501
expected random extras                       = 2^101.7929
   doublings  (0.828 - 0.678)*(1/(6R) + 2^-37.99999)*(S_ops + 2^60)     2^101.39
   traversals (0.702 - 0.678)*(2/(6R) + 64*2^-37.99999)*(S_ops + 2^60)  2^99.75
   starts 2^13*(S_ops/2^64 + 2^20)                                     2^75.71
   words, attempts, probes, small batches, fallbacks, zero-test re-runs, length carries   2^70
allowance A = ceil(extras/0.000165)          = 2^114.3582
W_ledger = phase1 + hard + A                 = 2^126.1538504
```

`certificate.py` uses the exact rates; `experiment.py` rounds them up to 2^-22.5849 and 2^-21.5839. Both give
W_ledger = 2^126.1538504.

In exact milli-units (`certificate.py`):

```text
S_int             = 139546939379019289183162248360078303600
phase 1           = 94612824898975078066965685168256532090528
checkpoint writes = 4360841855594352787045877855290374924
hard (incl. it)   = 4541567652431436229797006445252884642
extras            = 4392572954690847306645243065494098
A                 = 26621654270853620040274200396933927273
W_ledger          = 94643988120898363123235756375098718902443
```

W_ledger = 2^126.1538504 < 2^126.1539. The next lower four-decimal cap, 126.1538, is below W_ledger and is not
claimed. Against the cap's headroom (cap − phase 1 − hard = 2^114.524591), the Markov ratio is
0.000147022 < 0.000165, so the meter halts early with probability below h and

```text
success >= 0.3901753 - 0.000165 - 2e-7 = 0.3900101 > 0.39001 > 0.39.
```

`certificate.py` checks everything in exact integer and rational arithmetic: the kernel identities
(88, 67, 108, 493), the 77 rows, the prices 678, 828 and 702, the cycle sums, every ledger row including the
checkpoint writes and the table, and W^10000 < 1000^10000·2^1261539 and W^10000 > 1000^10000·2^1261538. It
prints `CERTIFICATE: PASS`. The counter widths, preprocessing cap and memory are re-checked there.

## 6. Provenance, files and limitations

The construction, report store, ledger structure and the analytic rows not listed in section 3 are credited, as
design provenance, to the public six-class package line reproduced in `upstream-proof.md`, and to 02feb284
(126.3059), whose q29 census calibrated κ. The counted kernels and loops of section 2, the earlier row cuts, the
doubling itemization and the checkpoint action are this solver's accepted daa13d54 (126.2165), 09990b8a
(126.2041), 1ca8ae82 (126.1975) and 5d3bf6c4 (126.1913), the last being the base of this package. The table size
R = 2^20 with checkpoint-only cycle detection is the design of the accepted BLS12-381 G1 package 67154fda. These
are declared as dependencies. New here: that walk change applied to this attack and its accounting, the
fingerprint row (18 → 8), the census of section 5.2, and the retune h = 0.000165, p0 = 0.3901753. No prior
acceptance, review outcome or score is a premise.

- `proof.md`: this argument. `upstream-proof.md`: the construction (abridged; its kernel, row, price, walk and
  ledger numbers are superseded here).
- `experiment.py`: part A production checks, the toy pipeline (the checkpoint-only walk of section 5.1) and the
  kernels mode of section 2.6. `certificate.py`: the exact certificate.
- `census_cp.py`, `census-cp-q29.json`: the census of section 5.2. The public census of the two-cycle walk at
  R = 2048 (`census.py`, `census-q29.json`) is in 02feb284's package.
- `kernels/gen_loops.py` (generator), `kernels/secp_loops.s` (counted loops), `kernels/secp_tests.s`
  (emitter wrappers), `kernels/secp_fe.s` and `kernels/secp_walk.c` (inversion kernels and their C source,
  unchanged), `kernels/strip_asm.py`, `kernels/build.sh`.

The limitations are the upstream ones. H1 is a random-mapping heuristic supported at toy scale. The step price
is the executed count of the shipped loops plus analytic rows. The attack is generic and impractical, and the
bound is not a claim about deployed secp256k1.

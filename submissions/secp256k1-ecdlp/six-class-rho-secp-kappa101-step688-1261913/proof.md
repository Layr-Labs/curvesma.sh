# Six-class secp256k1 rho: a 688-instruction step (126.1913)

**Claim.** For a uniform full-production secp256k1 input `Q=[x]P`, the classical cold six-class
(negation + GLV) parallel Pollard-rho attack of `upstream-proof.md`, with the counted kernels and
bookkeeping rows of this document, recovers `x` with probability greater than `0.39001` under the
explicit H1 heuristic. A global meter charges every action before it runs, so on every run its total
work is at most

```text
time_log2 = 126.1913 = 1261913/10000
```

target-group-operation equivalents (C = 1000 scalar AArch64 instructions per unit). Peak memory is
below `2^72` bytes; preprocessing is capped at `2^59` units, charged and included; nonuniform advice is
zero. This is an analytic heuristic bound, not a practical break.

## 1. What changes

This package is this solver's accepted 1ca8ae82 (126.1975) with two more bookkeeping-row cuts (section 3).
The table rows are padded to eight entries, so the table address needs no twist index. A copy of the
checkpoint fingerprint sits beside each history slot. p0 is retuned (section 5). The counted kernels and
loops of section 2 are unchanged, byte for byte (daa13d54's). The walk, class
key, canonical index, escapes, merging lemma, phases, meter, report store, every other hard row, the H1
calibration (`kappa = 1.01`) and the Markov headroom `h = 0.003866` are those of 02feb284 (126.3059) and its
public ancestors (`upstream-proof.md`, sections 1 to 4).

| item | 02feb284 | 1ca8ae82 | this package |
| --- | ---: | ---: | ---: |
| fe_mul (product + reduction) | 91 = 63 + 28 | 88 = 62 + 26 | 88 |
| fe_sqr | 71 | 67 | 67 |
| x3 and y3: subtract, then canonicalise | 9 + 8 each | 14 each (fe_subcanon) | 14 each |
| R between pass-2 iterations | stack slot: 2 LDP + 2 STP | registers x25-x28 | registers |
| y(W) in pass 2 | loaded twice | loaded once | loaded once |
| counted loops, pass 1 / pass 2 | 111 / 521 | 108 / 493 | 108 / 493 |
| analytic rows | 112 | 90 | **87** |
| ordinary step | 744 | 691 | **688** |
| doubling / traversal | 950 / 768 | 841 / 715 | **838 / 712** |
| p0 / success floor | 0.3939163 / 0.3900501 | 0.393885 / 0.3900188 | **0.39388 / 0.3900138** |
| claim | 126.3059 | 126.1975 | **126.1913** |

The ordinary step is charged exactly its table, with no margin: **108 + 493 + 87 = 688 instructions
= 0.688 units**.

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

## 3. Analytic rows: 112 → 87

The rows cover what the arrival pass does per class arrival outside the two counted loops (pseudocode of
`upstream-proof.md` section 2). Four rows are unchanged from 02feb284: range test 8, decision 13,
arrival traffic 3, inversion and call share 1. Six change:

| row | 02feb284 | 1ca8ae82 | this | sequence |
| --- | ---: | ---: | ---: | --- |
| F1 truncated product | 20 | 17 | **17** | 2 LDP of the x limbs; 4 MUL and 4 UMULH against G_1..G_5 held in registers; 7 ADD |
| sign, p − y, K | 17 | 13 | **12** | p − y as SUBS from a register holding p0 = 2^64 − C0, then 3 SBCS from a register holding 2^64 − 1 (4); compare y with p − y by SUBS, 3 SBCS (4); 4 CSEL for K (4); s is not set in a register, the table address reads it from the same flags |
| h, DP test, table address | 7 | 7 | **6** | AND for h (1); CBZ on limb 1 of K (1); ADD base + (h << 9), ADD + (j << 6), ADD #256 (3); CSEL on the carry of the y-against-(p − y) compare (1) |
| fingerprints, history, checkpoint | 26 | 19 | **18** | CBNZ on the pending flag (1; the flag is loaded by the length row's LDP); ring address (2); 3 LDP loading K_(−2) and the slot's copy of K_cp (3); CMP, 2 CCMP, branch against K_(−2) (4); slot write STP + STR (2); CMP, 2 CCMP, branch against K_cp (4); checkpoint test TST, branch (2), its stores being the separate checkpoint-write action |
| length, cap, meter, loop | 13 | 9 | **9** | LDP of the low length limb and the pending flag, ADDS #1, STR, B.CS to the carry path (4); meter (3); arrival loop (2) |
| zero-denominator test | 4 | 0 | **0** | one test per batch (section 4) |
| **total** | **112** | **90** | **87** | |

- **p − y.** p's limbs are (p0, 2^64 − 1, 2^64 − 1, 2^64 − 1). With p0 and 2^64 − 1 held in two registers for
  the whole arrival pass (set once per batch, 3 instructions, section 4), p − y is one SUBS/SBCS chain.
  The old row negated y from zero and then subtracted C0 (8 instead of 4). y is canonical and nonzero (no
  point of order 2), so p − y ∈ [1, p − 1] with no borrow out. y ≠ p − y because p is odd, so the carry of
  y − (p − y) is s = [p − y < y], and K = min(y, p − y) is selected by four CSEL.
- **Table rows of eight entries.** Tab[h][t] with t = 3s + j is laid out as Tab[h][s][j] at
  base + 512h + 256s + 64j: each h owns 512 bytes, two halves of four 64-byte entries, of which j uses three.
  The address is AND, three ADDs and one CSEL on the carry of the y-against-(p − y) comparison, and that
  carry is exactly s. CSEL, AND, CBZ and ADD leave the flags alone, so nothing between that comparison and the
  CSEL changes them. The old rows formed t = 3s + j (2 ADD), base + 64·(6h + t) (3 ADD) and s (CSET). The
  table grows from 2048·384 to 2048·512 bytes (1 MiB); the padding entries are never read, and the table's
  construction and charges are unchanged.
- **LDP/STP.** fp(K) is three words stored contiguously, so each load is LDP + LDR and each store is STP + STR.
  This is the convention the counted loops use, and the one 02feb284 applied to F1 and arrival traffic.
- **G_1..G_5 in registers.** F1 multiplies the four x limbs by the five limbs G_1..G_5 of G = ⌊2^384·β/p⌋
  (`upstream-proof.md` section 1.1). The arrival pass keeps them in five registers for the whole pass, loaded
  once per batch (LDP, LDP, LDR: 3 instructions in the per-batch overhead of section 4), so F1 loads only the
  x limbs. Register map of the arrival pass: x0 loop counter; x1–x4 walk-slot, pointer-array, trail-record
  and table pointers; x5, x6 the processor's block remainder and operation count; x7–x11 G_1..G_5; x12 p0;
  x13 2^64 − 1; x18 unused; x29 and x30 saved by the prologue. That leaves x14–x17 and x19–x30 (16 registers)
  for one walk's values. The most live at once is y, p − y and j while K is formed (9 words); the
  history compares hold K's three fingerprint words, three loaded words and the table index and address.
- **K_cp beside each history slot.** Each of the two ring slots is 64 bytes: the three fingerprint words of
  K_(−1) or K_(−2), then a copy of K_cp's three words. The arrival reads the slot of K_(−2) with three LDP,
  which load K_(−2) and K_cp together (four loads before). The slot write still stores only the fingerprint
  part. A checkpoint write stores K_cp in both slots, and the escape clears both copies. The escape's
  clearing is then 6 STP for the two slots and 1 STR for the flag; with the order test against K_(−1)
  (LDP + LDR, CMP, 2 CCMP, branch: 6) that is 13, inside the escape logic's 30.
- **Pending flag with the length.** The trail record stores the pending-escape flag in the word after the
  low length limb, so one LDP loads both. The fingerprint row keeps only the CBNZ on the flag.
- **Length.** Only the low limb is updated per arrival: ADDS #1 and STR, then B.CS to a carry path. The
  carry path runs only when the low limb wraps. It loads, increments and stores the high limb, compares it
  with 20 and abandons the trail at 20, at most 6 instructions. The length grows by one per operation, so it
  first reaches L_max = 20·2^64 exactly at such a carry, and testing the cap there alone is exact. A trail
  wraps at most once per 2^64 of its operations, so a run has fewer than (W/0.688)/2^64 < 2^63 carries. The
  carry path is metered as its own action at 0.01 units; its expected total is below 2^57 units, inside the
  2^70 "other" random row. The low limb still gives length mod 2 (ring) and length mod 64 (checkpoint).
- **Checkpoint stores.** The checkpoint write K_cp ← fp(K) runs only on arrivals with length mod 64 = 0.
  The row keeps the TST and branch that decide it (2) and no longer charges its stores on every arrival.
  Writing both copies takes 2 STP + 2 STR, metered as their own action at 0.004 units when they run. Apart
  from a trail's first write (at length 0, covered by the 2^13-unit start price), a trail writes at most once
  per 64 of its operations. Phase 1 performs fewer than S_ops + 2^61 operations, so it makes at most
  ⌈(S_ops + 2^61)/64⌉ writes. Section 5 charges that bound, at 4 instructions each, as a hard row. Phase-2
  and re-walk operations are priced at 600 units each and cover their own writes.
- **Zero denominator.** See section 4.

## 4. Prices and per-batch overhead

- **Ordinary step:** 108 + 493 + 87 = **688** instructions = 0.688 units.
- **Traversal step:** 688 + 24 = **712** (the 24 instructions to track and resume the minimum-fp member are
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
  - **Escape logic: 30**, as upstream.

  The doubling is 657 + 30 + 87 + 34 + 30 = **838**. The weak additions, the 18 and the 30 are those of
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
refills 0.03, words, attempts, table, store initialisation) is unchanged; the length carry (0.01) and the
checkpoint write (0.004, both copies) are their own actions.

## 5. H1 calibration and the exact budget

H1, its q29 census (`census.py`, `census-q29.json`: mean slowdown +0.6968 %, pooled s.e. 0.1169 %) and the
allowance `kappa = 1.01` are unchanged and remain the disclosed heuristic of `upstream-proof.md` section 4,
supported only at toy scale. The Markov headroom keeps 02feb284's h = 0.003866. The floor p0 is 0.39388
(1ca8ae82 used 0.393885), which keeps the success floor at 0.3900138 > 0.39001 > 0.39 with W below
2^126.1913.

```text
N'     = ((r-1)/6) * 2048/2047
S0     = sqrt(2*N'*ln(1/(1-0.39388)))        = 2^126.7088475
S_ops  = 1.01*S0                             = 2^126.7232028
phase1 = 0.688*(S_int + 2^60)                = 2^126.183683
hard extras                                  = 2^112.786822
   checkpoint writes 0.004*ceil((S_int + 2^61)/64)        2^112.757419
   upstream rows (refills at the new cap)                  2^107.155501
expected random extras                       = 2^110.58090
   doublings  (0.838 - 0.688)*2^-13.5681*(S_ops + 2^60)   2^110.42
   traversals (0.712 - 0.688)*64*2^-19.99*(S_ops + 2^60)  2^107.35
   starts 2^13*(S_ops/2^64 + 2^20)                        2^75.72
   words, attempts, probes, small batches, fallbacks, zero-test re-runs, length carries   2^70
allowance A = ceil(extras/0.003866)          = 2^118.59584
W_ledger = phase1 + hard + A                 = 2^126.1912960
```

`certificate.py` bounds the doubling rate by 1/(6R) + 2^-19.9954 and `experiment.py` by 2^-13.5681, so their
extras differ in the seventh decimal of the exponent; both give W_ledger = 2^126.1912960.

In exact milli-units (`certificate.py`):

```text
S_int             = 140438180397825303060093422393133512515
phase 1           = 96621468113703808506137484601645367329808
checkpoint writes = 8777386274864081441399954087646700408
hard (incl. it)   = 8958112071768986461047455588269860686
extras            = 1941638144468751421382608696286488673
A                 = 502234388119180398702175037839236594155
W_ledger          = 97132660613894757891300707095072873784649
```

W_ledger = 2^126.1912960 < 2^126.1913. The next lower four-decimal cap, 126.1912, is below W_ledger and is not
claimed. Against the cap's headroom (cap − phase 1 − hard = 2^118.596619), the Markov ratio is
0.00386392 < 0.003866, so the meter halts early with probability below h and

```text
success >= 0.39388 - 0.003866 - 2e-7 = 0.3900138 > 0.39001 > 0.39.
```

`certificate.py` checks everything in exact integer and rational arithmetic: the kernel identities
(88, 67, 108, 493), the 87 rows, the prices 688, 838 and 712, every ledger row including the checkpoint
writes, and W^10000 < 1000^10000·2^1261913 and W^10000 > 1000^10000·2^1261912. It prints
`CERTIFICATE: PASS`. The counter widths, preprocessing cap (2^58.864 < 2^59) and memory (< 2^72 bytes) are
unchanged and re-checked there.

## 6. Provenance, files and limitations

The construction, report store, H1 calibration, ledger structure and the analytic rows not listed in section
3 are credited, as design provenance, to the public six-class package line reproduced in `upstream-proof.md`,
and to 02feb284 (126.3059). The counted kernels and loops of section 2, the earlier row cuts, the doubling
itemization and the checkpoint action are this solver's accepted daa13d54 (126.2165), 09990b8a (126.2041) and
1ca8ae82 (126.1975), the last being the base of this package. All are declared as dependencies. New here: the
table rows padded to eight entries (sign 13 → 12, address 7 → 6), the copy of K_cp beside each history slot
(fingerprints 19 → 18) and p0 = 0.39388. No prior acceptance, review outcome or score is a premise.

- `proof.md`: this argument. `upstream-proof.md`: the construction (abridged; its kernel, row, price and ledger
  numbers are superseded here).
- `experiment.py`: part A production checks, toy pipeline and the kernels mode of section 2.6.
  `certificate.py`: the exact certificate. `census.py`, `census-q29.json`: the H1 census.
- `kernels/gen_loops.py` (generator), `kernels/secp_loops.s` (counted loops), `kernels/secp_tests.s`
  (emitter wrappers), `kernels/secp_fe.s` and `kernels/secp_walk.c` (inversion kernels and their C source,
  unchanged), `kernels/strip_asm.py`, `kernels/build.sh`.

The limitations are the upstream ones. H1 is a random-mapping heuristic supported at toy scale. The step price
is the executed count of the shipped loops plus analytic rows. The attack is generic and impractical, and the
bound is not a claim about deployed secp256k1.

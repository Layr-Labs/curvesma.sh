# Six-class secp256k1 rho at R = 2^20: a 623-instruction step with 53 bookkeeping rows (126.0324)

**Claim.** For a uniform full-production secp256k1 input `Q=[x]P`, the classical cold six-class
(negation + GLV) parallel Pollard-rho attack of `upstream-proof.md`, with the counted kernels and
bookkeeping rows of this document, recovers `x` with probability greater than `0.39001` under the
explicit H1 heuristic. A global meter charges every action before it runs, so on every run its total
work is at most

```text
time_log2 = 126.0324 = 1260324/10000
```

target-group-operation equivalents (C = 1000 scalar AArch64 instructions per unit). Peak memory is
below `2^72` bytes; preprocessing is capped at `2^59` units, charged and included; nonuniform advice is
zero. This is an analytic heuristic bound, not a practical break.

## 1. What changes

This package is this solver's rows package A (submission 955d3431, claim 126.0707: R = 2^20, checkpoint-only
cycle detection, rare-event field kernels, 67 rows) with three changes to the per-arrival bookkeeping, from 67
rows to 53 (sections 3 and 5.1), and x3 formed by one fused subtraction (section 2.4):

- **Anchors in global rounds, cap in anchor periods** (−6). Every walk makes exactly one operation per batch
  round, so setting K_cp in rounds t ≡ 0 (mod 64) is the per-walk rule (length + φ) ≡ 0 (mod 64) with the
  trail's phase φ = (start round) mod 64. φ is stored in the trail record, so re-walks replay the same
  schedule. The per-arrival anchor test and length update go; anchor rounds run a variant of the arrival pass
  that writes K_cp and counts the trail's anchor periods for the cap, a hard row.
- **Even-y class key, two limbs per arrival** (−8). K is the even one of y and p − y (exactly one is even,
  since p is odd), so the sign is one TST of y0. Each arrival forms only K0 and K1, and loads only y0 and y1.
  K2 and K3 are formed on the rare paths that need them: reports, filter hits, traversal steps and anchor
  rounds. The partition index is bits 1..20 of K, by UBFX.
- **Filter on K0.** The checkpoint filter compares K0 instead of K3, and a cleared checkpoint holds the odd
  marker 1 in that word, which no key has.
- **fe_sub2 for x3** (−3 in pass 2). x3 = λ^2 − x_W − x_T is one fused subtraction of 15 instructions instead of
  two fe_sub (18), after the accepted b5de673. This package's version reads x_T in two halves inside its borrow
  chain, so it needs no extra register; its one rare event is a third failure row.

The other counted kernels, the cycle rates, H1 and the store are those of package A. The toy, Part A and the
census move to the even-y key and the per-trail phase (section 5.2).

| item | package A (955d3431) | this package |
| --- | ---: | ---: |
| class key | min(y, p − y) | **even of y, p − y** |
| anchors, trail cap | per walk, length | **global rounds, periods** |
| x3 | two fe_sub (18) | **fe_sub2 (15)** |
| counted loops, pass 1 / pass 2 | 105 / 468 | 105 / **465** |
| analytic rows | 67 | **53** |
| ordinary step | 640 | **623** |
| doubling / traversal | 787 / 666 | **770 / 656** |
| claim | 126.0707 | **126.0324** |

The table size (R = 2^20), checkpoint-only cycle detection, h = 0.000165, p0 = 0.3901753 and the success floor
0.3900101 are unchanged.

The ordinary step is charged exactly its table, with no margin: **105 + 465 + 53 = 623 instructions
= 0.623 units**.

## 2. The counted kernels

`kernels/gen_loops.py` emits the two loops as straight-line AArch64 (`kernels/secp_loops.s`) and, from the
same emitters, the test wrappers `t_mul`, `t_sqr`, `t_sub`, `t_sub2` and `t_reduce` (`kernels/secp_tests.s`).
p = 2^256 − C0 with C0 = 2^32 + 977, so 2^256 ≡ C0 (mod p). A field element is four little-endian 64-bit
limbs; **weak** means any value in [0, 2^256), **canonical** means [0, p). Walk points and table points are
canonical, except after the rare event of section 2.4. x4 holds C0 inside the loops.

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

### 2.2 The reduction (23)

z = L + H·2^256 with L, H < 2^256 (the eight product limbs).

1. **H·C0** (8 + 4): four MUL and four UMULH give g_j, k_j with g_j + 2^64·k_j = z_(4+j)·C0, k_j ≤ C0 − 1.
   ADDS g1+k0, ADCS g2+k1, ADCS g3+k2, ADC k3+c gives H·C0 in five limbs; H·C0 < 2^256·C0, so the ADC
   cannot overflow.
2. **s = L + H·C0** (5): ADDS, 3 ADCS, then ADC into the fifth limb s4. Since
   s ≤ (2^256 − 1)(1 + C0) < (1 + C0)·2^256, s4 ≤ C0.
3. **w = s4·C0** (2): MUL and UMULH; w ≤ C0^2 < 2^65.
4. **t = (s mod 2^256) + w** (4): ADDS, 2 ADCS and ADC, which return t mod 2^256. t ≡ z (mod p) because
   2^256 ≡ C0, and t < 2^256 + C0^2.

**23 = 8 + 4 + 5 + 2 + 4.** There is no final fold (f99114e0 spent 3 instructions on it, the public package
5). If t < 2^256 the result is t, a weak value ≡ z (mod p). If t ≥ 2^256, then C0 ≤ t − p < C0^2 + C0 and
t − p ≡ z, so z mod p = t − p and the result t − 2^256 = z mod p − C0 is wrong. That needs
z mod p < C0^2 + C0, probability (C0^2 + C0)/p < 2^-191.99 for a uniform residue. The reduced values (the
batch products and inverses, λ, λ^2 and λ·(x_W − x3)) are products of pseudo-random field elements, the
premise of the exceptional-addition row, so the event is a failure row of 7·(S_ops + 2^61)·(C0^2 + C0)/p =
2^-62.48 per run, with at most 7 reductions per batched operation (section 5.3). The experiment constructs
such inputs directly (section 2.6): for a chosen high half H it solves for the low half L that puts
s mod 2^256 inside [2^256 − s4·C0, 2^256), including the boundary t = 2^256, and checks that the kernel
returns exactly z mod p − C0.

### 2.3 fe_sqr (64)

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
- **Reduction** (23).

**64 = 20 + 8 + 6 + 7 + 23.**

### 2.4 fe_sub (9) and fe_sub2 (15)

**fe_sub(a, b)**, a weak, b canonical, is unchanged: SUBS, 3 SBCS; CSEL of C0 on borrow; SUBS, 2 SBCS, SBC.
On borrow the wrapped value is at least C0 + 1, so subtracting C0 cannot borrow again; the result is weak and
≡ a − b. Without a borrow the result is a − b, which is at least p only if a ≥ p and b ≤ a − p < 2^256 − p =
C0. So it is canonical whenever a is, and whenever b ≥ C0. The step forms `y3 = fe_sub(λ·(x_W − x3), y_W)`
and `x_W − x3 = fe_sub(x_W, x3)`, which needs x3 canonical.

**fe_sub2(a, b, c)**, a weak, b and c canonical, forms x3 = λ^2 − x_W − x_T (after the accepted b5de673):

```text
SUBS, 3 SBCS        d = a − b mod 2^256, borrow β1                          (4)
CSEL m1, C0, xzr    m1 = β1·C0                                              (1)
LDP, SUBS, SBCS,    d ← d − c mod 2^256, borrow β2; c's limbs are loaded    (4 + 2 LDP)
LDP, SBCS, SBCS     in two halves inside the chain (LDP keeps the flags)
CSEL m2; ADD        m = (β1 + β2)·C0                                        (2)
SUBS, 2 SBCS, SBC   result = d − m                                          (4)
```

15 instructions plus the two LDP that loaded x_T before. d = a − b − c + (β1 + β2)·2^256, so d − m ≡ a − b − c.
With two borrows the result is a − b − c + 2p ∈ [2, p − C0). With none it is a − b − c, canonical unless
b + c < C0, which needs x_W < C0. With one borrow the result is a − b − c + p ∈ [0, p), except when
d = (a − b − c) mod 2^256 < C0. Then the last subtraction wraps and the result is ≡ a − b − c + C0, a wrong
value: probability C0/2^256 < 2^-223.99 per use for a pseudo-random λ^2, a failure row of
(S_ops + 2^61)·C0/2^256 = 2^-97.29 per run (section 5.3). Reading x_T in halves keeps the register peak of
pass 2 at 24 (25 before), so no register outside the pool is needed. So x3 and y3 are canonical and correct
unless x_W or y_W is below C0 (2^-96.29 per run) or the fe_sub2 event occurs. The doubling uses the same
kernels (section 4).

### 2.5 The two loops (105 and 465)

**Pass 1** (forward; the schedule of the public package with the new fe_mul): LDR of the table pointer and two
LDP each for x(T_i) and x(W_i) (5); dx_i = fe_sub(x(T_i), x(W_i)) (9, canonical); two STP each for acc_(i−1)
and dx_i into slot i (4); acc_i = dx_i·acc_(i−1) with the reduction writing into the four loop-carried
registers (85); SUBS + B.NE (2). **105.**

**Pass 2** (backward, i = n−1 … 1). R ≡ 1/(dx_0···dx_i) stays in x25–x28 across iterations: it is loaded once
from Rinv before the loop and never stored (the public package kept it in a stack slot, two LDP and two STP
per walk). Per iteration:

| work | instructions |
| --- | ---: |
| LDR table pointer; 2 LDP acc_(i−1); inv = R·acc_(i−1) | 1 + 2 + 85 |
| 2 LDP dx_i; R ← dx_i·R (result into x25–x28) | 2 + 85 |
| 2 LDP y(T_i), 2 LDP y(W_i) (kept to the end); dy = fe_sub(y_T, y_W) | 4 + 9 |
| λ = dy·inv; λ^2 | 85 + 64 |
| 2 LDP x(W_i) (kept); x3 = fe_sub2(λ^2, x_W, x_T), with x(T_i)'s 2 LDP inside it | 2 + 15 + 2 |
| x_W − x3 = fe_sub(x_W, x3); 2 STP x3 | 9 + 2 |
| λ·(x_W − x3); y3 = fe_sub(·, y_W); 2 STP y3 | 85 + 9 + 2 |
| SUBS + B.NE | 2 |
| **total**: 4·85 + 64 + 3·9 + 15 + 13 loads + 4 stores + 2 | **465** |

Walk 0 uses R itself as its inverse after the loop. This is the standard affine addition
λ = (y_T − y_W)/(x_T − x_W), x3 = λ^2 − x_W − x_T, y3 = λ(x_W − x3) − y_W, with canonical outputs outside the
event of section 2.4.

**Registers.** The pool is x5–x17 and x19–x30 (25 registers). x29 and x30 are saved with x19–x28 by the
prologue and restored before RET; the loops call nothing. x0–x3 hold the count and the three array pointers, x4
holds C0, x18 is never written. The register peak is 24 in pass 2 and 16 in pass 1; there are no spills and
no stack access inside either loop. Both loops are single basic blocks: every conditional is a carry chain
or CSEL, and the only branch is the back edge. Opcode histograms:

- pass 1: LDR 1, LDP 4, STP 4, MUL 21, UMULH 21, ADDS 9, ADCS 25, ADC 9, SUBS 3, SBCS 5, SBC 1, CSEL 1,
  B.NE 1 (105);
- pass 2: LDR 1, LDP 12, STP 4, MUL 99, UMULH 99, ADDS 42, ADCS 116, ADC 42, ADD 1, SUBS 10, SBCS 23,
  SBC 4, CSEL 5, EXTR 5, LSR 1, B.NE 1 (465).

Counting is the convention of the public package: one LDP or STP is one instruction.

### 2.6 Evidence

`experiment.py` with `"kernels": true` (the declared initial run) does all of the following:

1. **Limb ports.** Python ports of the generated fe_mul (the six chains, as data), the reduction, fe_sqr,
   fe_sub and fe_sub2 assert every no-overflow step above: each carry-only ADC, the exact product identity,
   s4 ≤ C0 and w1 ≤ 1. When the last addition of a reduction carries, the port asserts z mod p < C0^2 + C0
   and that the result is z mod p − C0. They also assert that fe_sub never borrows twice, that a
   non-canonical fe_sub output has a ≥ p and b < C0, and that an fe_sub2 wrap has one borrow and gives
   a − b − c + C0. They are checked against big integers on all pairs of 44 edge values (0, 1,
   p − 1, p, p + 1, 2^256 − 1, values near p and near 2^256, every all-ones/all-zeros limb pattern) and on
   random values; every result must equal the port's value, which is z mod p up to weak form except on the
   predicted event. The reductions include 200 constructed inputs whose last addition carries. A 24-walk
   batch built from the ports is checked against affine secp256k1 addition.
2. **The shipped loops, executed.** The AArch64 interpreter executes `kernels/secp_loops.s` on production
   points for batches of 2, 3, 5, 9 and 17 walks, compares every output with affine arithmetic, and reports
   the executed-instruction difference per extra walk. It must be exactly 105 (pass 1) and 465 (pass 2).
3. **The emitters, executed.** The interpreter executes `kernels/secp_tests.s`, the same emitters wrapped as
   functions, on edge and random inputs: 730 products, 444 squares, 818 fe_sub pairs, 1,630 fe_sub2 triples
   (half of them built near the single-borrow wrap) and 402 reductions, at least 200 of which carry out of
   the last addition. Every result must equal the port's value under the same rule.
4. **Counts and ledger.** It parses the instruction counts from the shipped files, checks the identities
   105 = 85 + 9 + 5 + 4 + 2 and 465 = 4·85 + 64 + 3·9 + 15 + 13 + 4 + 2, and recomputes the doubling bound,
   the per-batch overhead and the ledger of section 5.

The interpreter of the public package gains one standard A64 form used by fe_sqr: **ADDS with an LSL-shifted
register** (xn + (xm << k mod 2^64), with the carry of that 64-bit sum). The CSNEG form of f99114e0 is no
longer used. All other forms are unchanged.

Before submission, the generated files were also assembled with clang 23 for `aarch64-linux-gnu` and the
object code was run in the Unicorn 2.1.4 AArch64 CPU emulator, off the hosted runner and not shipped, against
a bit-exact model of the kernels. It covered 4,031 walks in random batches of 2 to 199, 4,304 products (41
edge-value products whose last reduction addition carries), 2,048 squares, 11,604 fe_sub pairs (15
non-canonical, all with a ≥ p and b < C0), 63,700 fe_sub2 triples (186 constructed single-borrow wraps) and
4,538 reductions (2,535 constructed carries). There were 0 mismatches, every rare case gave its predicted
value, counts were exactly 105 and 465 per walk, and x19–x30 were restored on return.

## 3. Analytic rows: 112 → 53

The rows cover what the arrival pass does per class arrival outside the two counted loops (pseudocode of
`upstream-proof.md` section 2, with the walk changes of section 5.1). Three rows are unchanged from 02feb284:
range test 8, decision 13, inversion and call share 1. Seven change:

| row | 02feb284 | A | this | sequence |
| --- | ---: | ---: | ---: | --- |
| F1 truncated product | 20 | 13 | **13** | 2 LDP of the x limbs; 4 UMULH of the column-4 products x_i·G_(4−i) and 3 ADD; 4 MADD folding in the column-5 low halves x_i·G_(5−i) (G_1..G_5 held in registers) |
| sign, p − y, K | 17 | 12 | **5** | p − y in limbs 0 and 1: SUBS from the register holding p0 = 2^64 − C0, SBC from the register holding 2^64 − 1 (2); TST y0, #1, whose Z flag is set exactly when y is even (1); CSEL of K0 and K1 on it (2) |
| h, DP test, table address | 7 | 6 | **6** | UBFX of bits 1..20 of K0 for h (1); CBZ on K1 (1); ADD base + (h << 9), ADD + (j << 6), ADD #256 (3); CSEL on the Z flag of the TST (1) |
| arrival traffic | 5 | 3 | **2** | LDP of y0, y1 (1); store of the new table pointer (1) |
| fingerprints, history, checkpoint | 26 | 5 | **3** | LDR of K_cp's limb-0 word (1); CMP with K0, B.EQ to the filter-hit action (2); K_cp is written in anchor rounds (section 5.1) |
| length, cap, meter, loop | 13 | 6 | **2** | arrival loop (2); the cap is counted in anchor rounds and the meter runs once per batch (section 4) |
| zero-denominator test | 4 | 0 | **0** | one test per batch (section 4) |
| **total** | **112** | **67** | **53** | |

- **Even-y key in two limbs.** p's limbs are (p0, 2^64 − 1, 2^64 − 1, 2^64 − 1), and p0 and 2^64 − 1 are held
  in two registers for the whole arrival pass (set once per batch, section 4). y is canonical and nonzero (no
  point of order 2), so p − y ∈ [1, p − 1]. y and p − y have opposite parity because p is odd, so exactly one
  of them is even, and K is that one. A class's six points share their y up to sign, so K is a function of the
  class (part A checks it on production points). TST y0, #1 sets Z exactly when y is even. Two CSEL then pick
  K0 and K1 from y0, y1 and the low limbs of p − y (SUBS, SBC). s = [y odd] is the sign of the canonical member.
  K2 and K3 are formed only on the paths that need them: reports, filter hits, traversal steps and anchor
  rounds. There they take LDP of y2, y3, the four-limb SUBS/SBCS chain for p − y and two CSEL, 7 instructions.
  Bit 0 of K is always 0, so h is bits 1..20 of K. The distinguished bits 64..127, routing bits 128..143 and
  bucket bits 144..193 are unchanged.
- **Table rows of eight entries.** Tab[h][t] with t = 3s + j is laid out as Tab[h][s][j] at
  base + 512h + 256s + 64j: each h owns 512 bytes, two halves of four 64-byte entries, of which j uses three.
  The address is UBFX, three ADDs and one CSEL on the Z flag of the TST, which is exactly 1 − s. CSEL, UBFX,
  CBZ and ADD leave the flags alone. With R = 2^20 the table holds 2^20·512 bytes = 512 MiB; the padding
  entries are never read.
- **F1 by MADD; G_1..G_5 in registers.** F1 is the sum, mod 2^64 and with no carries, of the low halves of
  the four column-5 products x_i·G_(5−i) and the high halves of the four column-4 products x_i·G_(4−i), for
  G = ⌊2^384·β/p⌋ (`upstream-proof.md` section 1.1). So it is 4 UMULH and 3 ADD for the high halves, then
  4 MADD (xd = xa + xn·xm mod 2^64) that fold in the low halves: 13 with the 2 LDP of x, against 4 MUL, 4 UMULH
  and 7 ADD before. MADD is one scalar A64 instruction; F1 is an analytic row, and part A's limb-level
  emulation of the rule computes the same value. The arrival pass keeps G_1..G_5 in five registers for the
  whole pass, loaded once per batch (LDP, LDP, LDR: 3 instructions in the per-batch overhead of section 4).
  Register map of the arrival pass: x0 loop counter; x1–x4 walk-slot, pointer-array, trail-record and table
  pointers; x5, x6 the processor's block remainder and operation count; x7–x11 G_1..G_5; x12 p0;
  x13 2^64 − 1; x18 unused; x29 and x30 saved by the prologue. That leaves x14–x17 and x19–x30 (16 registers)
  for one walk's values; the most live at once is y0, y1, the low limbs of p − y and j while K is formed. x5
  stays live so that the rare paths can charge the meter (section 4).
- **Anchor rounds and the cap.** Every walk makes exactly one operation (a step, a doubling or a traversal
  step) per batch round. A processor counts its rounds t, and in rounds with t ≡ 0 (mod 64) it runs the anchor
  variant of the arrival pass (3 instructions per batch round, section 4). For every walk the variant forms K2
  and K3 (7), writes K_cp ← fp(K) (STP + STR, the checkpoint-write action) and adds one to the trail's count of
  anchor periods, abandoning the trail when the count reaches 20·2^58 (LDR, ADD, STR, CMP, B.HS: 5). A trail's
  first anchor comes within 64 operations of its start, so it is abandoned before L_max = 20·2^64 operations,
  and the trail-cap and re-walk bounds of `upstream-proof.md` hold as before. A walk starting a trail records
  φ = t mod 64 in the trail record (spare bits of its 8-byte tag and flags word). Its anchors fall where
  (length + φ) ≡ 0 (mod 64), so a re-walk or phase-2 continuation run from the record anchors where the live
  trail did and makes the same operations, inside its 600-unit price. Section 5.3 charges the variant's
  12 instructions per walk per anchor, ⌈(S_ops + 2^61)/64⌉ anchors, as a hard row.
- **Checkpoint filter.** With cycles caught by the checkpoint alone (section 5.1), the row loads the limb-0
  word of K_cp, compares it with K0 and branches on equality (LDR, CMP, B.EQ). Clearing the history stores the
  marker 1 in that word. No key has it, because K is even, so a cleared checkpoint never hits. A hit is its own
  action at 0.02 units, 15 instructions: SUBS of 20 milli-units from x5 and B.CC to the refill, the 7 that form
  K2 and K3, then LDP of the other two K_cp words (limbs 2 and 3), CMP, CCMP and B.NE back. Equal fingerprints
  are a true repeat and go on to the traversal (section 5.1). For a different class, K0 is a pseudo-random even
  word, so false hits occur at most 2^-63 per arrival. Section 5.3 charges the hits as the random row
  0.02·(2^-22.5849 + 2^-63)·(S_ops + 2^60), since true hits occur at most at the doubling rate. The write
  K_cp ← fp(K) is metered as its own action at 0.002 units. A trail writes once per anchor, plus its first
  write, which the 2^13-unit start price covers. Phase 1 performs fewer than S_ops + 2^61 operations, so
  section 5.3 charges ⌈(S_ops + 2^61)/64⌉ writes at 2 instructions each as a hard row.
- **Zero denominator.** See section 4.

## 4. Prices and per-batch overhead

- **Ordinary step:** 105 + 465 + 53 = **623** instructions = 0.623 units.
- **Traversal step:** 623 + 24 + 2 + 7 = **656**: the 24 instructions to track and resume the minimum-fp
  member (unchanged), the charge of its surcharge on its own path (SUBS of 33 milli-units from x5, B.CC), and
  the 7 that form K2 and K3 for the fingerprint the minimum compares.
- **Doubling:** the batched doubling shares the product chain and inversion. Its arithmetic is
  5 fe_mul (425) + 2 fe_sqr (128) + dx = 2y by one weak addition (9) + 3x^2 by two weak additions with a
  second carry fold (2·14) + x3 = fe_sub2(λ^2, x, x) (15) + x − x3 (9) + y3 = fe_sub(λ·(x − x3), y) (9) = 623,
  with the kernels and the rare events of section 2.
  - **Traffic: 30** (the public package allowed 72). The loads, stores, pointer updates and loop control are
    at most the addition loops' 5 + 4 + 2 and 13 + 4 + 2. A doubling walk reads its own y for dx = 2y in the
    forward pass and acc, dx, x and y in the backward pass, with no table point. It stores the same dx, acc,
    x3 and y3, and shares R in registers.
  - **List moves: 34** (the public package counted 50). Moving the walk to and from the doubling list is
    the upstream 16 loads and 16 stores of 64-bit words, the walk's x and y out and back. They pair into
    8 LDP and 8 STP, plus 18 other instructions.
  - **Escape logic: 30**, as upstream. With checkpoint-only detection and the meter per batch it is now the
    clearing of the history (one STR of the marker into K_cp) and the charge of the doubling's surcharge on
    its own path (SUBS of 147 milli-units from x5, B.CC), 3 instructions, well inside the 30.

  The doubling is 623 + 30 + 53 + 34 + 30 = **770**. The weak additions, the 18 and the 30 are those of
  `upstream-proof.md` section 5.3.
- **Zero-denominator test, once per batch.** Pass 1 multiplies the n canonical dx_i into acc_(n−1). dx_i ≡ 0
  only if dx_i = 0, and p is prime, so outside the rare event of section 2.2 acc_(n−1) ≡ 0 (mod p) exactly
  when some x(T_i) = x(W_i). Before the
  inversion, acc_(n−1) is canonicalised and tested: ADDS, 3 ADCS, 4 CSEL, 3 ORR, CBZ, 12 instructions per
  batch. Only when it fires does the batch take the cold path, an event of probability < 2^-114 over the whole
  run as before. Pass 1 is re-run with the per-walk test, each walk with dx_i = 0 is flagged abandoned, and its
  dx_i is replaced by 1 in the product and its slot, exactly as in the public package. The re-run, at most
  2^16·(105 + 4) < 2^23 instructions, is metered as its own action at 2^13 units and falls in the "other"
  random row.
- **Per-batch fixed work:** inversion I ≤ 24,322 (unchanged `secp_fe.s`, `upstream-proof.md` section 5.1–5.2),
  prologue, walk 0 and epilogue of the two passes (138 − 105 + 777 − 465 = 345), 64 for the call and
  return, 12 for the zero test, 3 for the constants p0 and 2^64 − 1, 3 for loading G_1..G_5, 5 for the
  meter and 3 for the round counter and the anchor-round test (ADD, TST, B.EQ): 24,757 < 2^15. Batches have
  at least 2^15 walks, so the share row stays 1.
- **Meter per batch.** Before a batch round of n walks runs, the processor charges its n operations at the
  step price in one subtraction from its block remainder x5 and adds n to its operation count x6: MOV of 623,
  MUL by n, SUBS, B.CC to the block refill, ADD, 5 instructions per batch. A round's charge is below
  2^16·623 < 2^26 milli-units, far below a 2^40-unit block, so one refill always covers it. Every walk of the
  round is then covered at the step price before its operation runs; a doubling or traversal walk pays its
  surcharge on its own path, before the doubling or traversal step runs, and a filter hit pays its own price
  the same way. So every action is still charged before it runs. Blocks are reserved,
  refilled and counted as spent exactly as in `upstream-proof.md` section 3, and phase 1 still stops at a
  refill once the operation count reaches S_ops. The 2^20 processors' last blocks cover fewer than
  2^20·2^40/0.623 = 2^60.69 further operations, so phase 1 performs fewer than S_ops + 2^61 operations, the
  count that section 3 and section 5.3 use.

Every other price (reports 40.4, probes 0.05, starts 2^13, fallbacks 0.5, phase-2 and re-walk operations 600,
refills 0.03, words, attempts, store initialisation) is unchanged; the checkpoint write (0.002) and the filter
hit (0.02) are their own actions, and the anchor-round work is a hard row (section 3). The table has 2^20
entries (section 5.1).

## 5. The walk at R = 2^20, H1 and the exact budget

### 5.1 The walk change

The walk is that of `upstream-proof.md` section 2 with four changes: the table has R = 2^20 rows instead of
2048, cycles are detected by the checkpoint alone, the class key is the even one of y and p − y, and anchors
fall where (length + φ) ≡ 0 (mod 64). Everything else is unchanged: the canonical index, the
twisted table Tab[h][t] = [u_t^{-1}]M_h, the distinguished classes
(limb 1 of K zero, θ = 2^-64), M = 2^20 walks, L_max = 20·2^64 and the checkpoint period 64. At each class
arrival, in order:

```
if K distinguished: report (K, a, b); end trail
if length >= L_max: abandon trail
if fp(K) == K_cp: traverse the cycle once more (<= 64 steps, else abandon), resume at the member of its
     minimum-fp class remembered during the traversal, ESC; clear history; continue
if (length + phi) mod 64 == 0: K_cp <- fp(K)     (phi: the trail's start round mod 64, section 3)
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
- **Key and partition.** This package's class key is the even one of y and p − y, and h is bits 1..20 of K
  (section 3); R enters the step only through the UBFX that forms h and through the table base. The partition
  field is disjoint from the distinguished bits 64..127, the routing bits 128..143 and the bucket bits
  144..193. The anchor phase φ only shifts where each walk's anchors fall; every walk is still anchored once
  per 64 of its own operations, which is all the merging lemma and the detection-delay bound use.
- **Precedent.** The accepted BLS12-381 G1 package 67154fda uses this walk at R = 2^20 with checkpoint-only
  detection. Its section on cycles makes the same rate and delay arguments.

### 5.2 H1 at R = 2^20

H1 keeps its form: P(T_c ≤ S_ops) ≥ 1 − exp(−(S_ops/κ)^2/(2N')) with κ = 1.01 and N' = N·R/(R − 1), now at
R = 2^20. T_c counts every class arrival, including the repeated arrivals inside fruitless cycles and during
the detection delay. Its other parts are as in `upstream-proof.md` section 4: pseudo-random key fields, the
non-fruitless return bound, and a uniform twist index independent of the partition index. κ = 1.01 was
calibrated by the public q29 census at R = 2048 (+0.70 % ± 0.12 %, 200,000 trials, eight tables). f99114e0's
census of its walk (key min(y, p − y), anchors at length ≡ 0 mod 64) measured +0.34 % and +0.55 % at R = 2^16
and 2^20. This package reruns that census (`census_cp.py`, `census-cp-q29.json`) for the walk as now specified:
even-y key, partition and distinguished fields one bit higher, and a uniform anchor phase per trail. It runs
on q29 (N = 2^24.83) with dp_bits 5 and a fresh table for every 4,000 trials, with the same seeds:

| R | mean T_c/√N (model) | slowdown | s.e.: trials / tables | P[T_c ≤ S0] | index-repeating cycles (model) | coincidental cycles (model) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2048 | 1.26258 (1.25362) | +0.72 % | 0.21 % / 0.20 % | 0.3874 | 33,984 (35,202) | 303 (400) |
| 2^16 | 1.26157 (1.25332) | +0.66 % | 0.21 % / 0.32 % | 0.3888 | 881 (1,087) | 585 (398) |
| 2^20 | 1.26298 (1.25331) | +0.77 % | 0.21 % / 0.24 % | 0.3870 | 31 (68.0) | 350 (397) |

Each row is 64,000 trials on 16 tables, about 4.27·10^8 operations.

- **Slowdown.** At R = 2^20 the slowdown is +0.77 % ± 0.21 %, and R = 2^16 and 2^20 pool to +0.71 % ± 0.15 %,
  inside κ − 1 = 1 %. The key is a relabeling of the same random-mapping model, so the four large-R runs of
  this census and f99114e0's pool to +0.58 % ± 0.10 %. P[T_c ≤ S0] is 0.3874, 0.3888 and 0.3870, within
  about one binomial standard error (0.0019) of H1's value at that budget, 1 − 0.605^(1/1.01^2) = 0.3890.
  κ = 1.01 is kept as the disclosed allowance.
- **R = 2048.** The checkpoint-only walk runs each 2-cycle until the next anchor catches it, about 37
  arrivals, all counted in T_c: about 0.3 % of the arrivals at R = 2048, under 2^-16.5 at R = 2^20 (section
  5.1). The hosted toy, which keeps R = 2048 because it draws a fresh table per trial, is not the H1 evidence
  for R = 2^20.
- **Charged cycles.** Each detected cycle is classified by the partition indices of its steps. The
  index-repeating cycles, the ones the ledger charges, are at or below the model's 1/(6R) + Σ, which ignores
  that a 2-cycle through a distinguished class ends the trail instead (1 in 32 here). At R = 2^20 a table holds
  only about N/(12R) ≈ 2.4 algebraic 2-cycles, so these counts swing widely between tables (f99114e0's census
  saw 77 against 67.8).
- **Toy self-meetings.** The coincidental cycles are the self-meetings of `upstream-proof.md` section 6: a
  trail of age ℓ meets its own past at about ℓ/N per step. Their model, Σ min(length + 1, 64)/N, gives
  0.93·10^-6 per operation; pooled over all 48 tables they occur at 0.97·10^-6. Their per-table counts are
  heavy-tailed, because every trial of a chunk walks the same table and many trails fall into the same few
  short cycles. In production they are bounded by the L_max^2/(2N) clause (section 5.1).
- **Trap test.** The long-trail test of `upstream-proof.md` section 6 (q29, no distinguished points, 150
  trails of 2,048 operations, seed 7), run with this walk through `experiment.py`'s `long_trails` mode, finds
  non-fruitless returns in 0/150 (R = 4) and 3/150 (R = 8) trails with the doubling escape, against 47/150
  and 14/150 with the additive control. Ordinary self-collisions predict about 10/150.

### 5.3 The exact budget

H1 and κ = 1.01 remain the disclosed heuristic, supported only at toy scale. With the random charges this
small, the Markov headroom is h = 0.000165 and the floor p0 = 0.3901753, which keeps the success floor at
0.3900101 > 0.39001 > 0.39 with W below 2^126.0324.

```text
N'     = ((r-1)/6) * 2^20/(2^20 - 1)
S0     = sqrt(2*N'*ln(1/(1-0.3901753)))      = 2^126.6996628
S_ops  = 1.01*S0                             = 2^126.7140181
phase1 = 0.623*(S_int + 2^60)                = 2^126.031322
hard extras                                  = 2^114.564105
   anchor rounds 0.012*ceil((S_int + 2^61)/64)            2^114.333196
   checkpoint writes 0.002*ceil((S_int + 2^61)/64)        2^111.748234
   upstream rows (table < 2^34, refills at the new cap)    2^107.155501
expected random extras                       = 2^102.0278
   doublings  (0.770 - 0.623)*(1/(6R) + 2^-37.99999)*(S_ops + 2^60)     2^101.36
   traversals (0.656 - 0.623)*(2/(6R) + 64*2^-37.99999)*(S_ops + 2^60)  2^100.21
   filter hits 0.02*(1/(6R) + 2^-37.99999 + 2^-63)*(S_ops + 2^60)      2^98.49
   starts 2^13*(S_ops/2^64 + 2^20)                                     2^75.71
   words, attempts, probes, small batches, fallbacks, zero-test re-runs                  2^70
allowance A = ceil(extras/0.000165)          = 2^114.5930
W_ledger = phase1 + hard + A                 = 2^126.0323512
```

`certificate.py` uses the exact rates; `experiment.py` rounds them up to 2^-22.5849 and 2^-21.5839. Both give
W_ledger = 2^126.0323512.

In exact milli-units (`certificate.py`):

```text
S_int             = 139546939379019289183162248360078303600
phase 1           = 86937743233129017161828350825698848808848
anchor rounds     = 26165051133566116722275267131742249544
checkpoint writes = 4360841855594352787045877855290374924
hard (incl. both) = 30706618785788973317809807713208974876
extras            = 5169187607603913790014307709990668
A                 = 31328409743054022969783683090852533334
W_ledger          = 86999778261657860158115944316502910317058
```

W_ledger = 2^126.0323512 < 2^126.0324. The next lower four-decimal cap, 126.0323, is below W_ledger and is not
claimed. Against the cap's headroom (cap − phase 1 − hard = 2^114.722467), the Markov ratio is
0.000150841 < 0.000165, so the meter halts early with probability below h. The 2·10^-7 covers every other
failure event: the rows of `upstream-proof.md` section 4 (trail cap < 10^-7, preprocessing cap < 2^-30 and
smaller ones) and the three kernel rows of section 2, the dropped carry 7·(S_ops + 2^61)·(C0^2 + C0)/p < 2^-62.4,
the small coordinate 2·(S_ops + 2^61)·C0/p < 2^-96 and the fe_sub2 wrap (S_ops + 2^61)·C0/2^256 < 2^-97.2,
which together sum to 1.0093·10^-7. So

```text
success >= 0.3901753 - 0.000165 - 2e-7 = 0.3900101 > 0.39001 > 0.39.
```

`certificate.py` checks everything in exact integer and rational arithmetic: the kernel identities
(85, 64, 105, 465), the 53 rows, the prices 623, 770 and 656, the cycle sums, every ledger row including the
anchor rounds, the checkpoint writes, the table and the filter hits, the three kernel failure rows and their sum
with the upstream rows, and W^10000 < 1000^10000·2^1260324 and W^10000 > 1000^10000·2^1260323. It
prints `CERTIFICATE: PASS`. The counter widths, preprocessing cap and memory are re-checked there.

## 6. Provenance, files and limitations

The construction, report store, ledger structure and the analytic rows not listed in section 3 are credited, as
design provenance, to the public six-class package line reproduced in `upstream-proof.md`, and to 02feb284
(126.3059), whose q29 census calibrated κ. The counted kernels and loops, the earlier row cuts, the doubling
itemization and the checkpoint action are this solver's daa13d54 (126.2165), 09990b8a (126.2041), 1ca8ae82
(126.1975) and 5d3bf6c4 (126.1913). The walk at R = 2^20 with checkpoint-only cycle detection, its census and
the retune h = 0.000165, p0 = 0.3901753 are this solver's f99114e0 (126.1539), after the design of the accepted
BLS12-381 G1 package 67154fda. The rare-event field kernels and their failure rows are this solver's 736c001a
(126.0931); F1 by MADD, the one-limb filter and the meter per batch are package A (955d3431, 126.0707), the base of
this package. Anchors in global batch rounds with the cap counted in anchor periods, and an even-y sign convention,
are the design of the accepted P-256 package 6d5e2adb; the fused subtraction for x3 is the accepted b5de673's
(built on 736c001a), here re-derived with x_T read in halves. These are declared as dependencies. New here: the
stored anchor phase that makes re-walks replay the live schedule exactly, the even-y key formed in two limbs
per arrival with its rare-path completions, the K0 filter with the odd marker, and their accounting. No prior
acceptance, review outcome or score is a premise.

- `proof.md`: this argument. `upstream-proof.md`: the construction (abridged; its kernel, row, price, walk and
  ledger numbers are superseded here).
- `experiment.py`: part A production checks, the toy pipeline (the checkpoint-only walk of section 5.1) and the
  kernels mode of section 2.6. `certificate.py`: the exact certificate.
- `census_cp.py`, `census-cp-q29.json`: the census of section 5.2. The public census of the two-cycle walk at
  R = 2048 (`census.py`, `census-q29.json`) is in 02feb284's package.
- `kernels/gen_loops.py` (generator), `kernels/secp_loops.s` (counted loops), `kernels/secp_tests.s`
  (emitter wrappers), `kernels/secp_fe.s` (the inversion kernels, unchanged). Their C source
  `secp_walk.c` and the build scripts `strip_asm.py` and `build.sh` are unchanged in 02feb284's and 736c001a's
  packages and are not repeated here; the experiment reads only the three assembly files.

The limitations are the upstream ones. H1 is a random-mapping heuristic supported at toy scale. The step price
is the executed count of the shipped loops plus analytic rows. The two kernel failure rows rest on the same
random-point premise as the exceptional-addition row. The attack is generic and impractical, and the
bound is not a claim about deployed secp256k1.

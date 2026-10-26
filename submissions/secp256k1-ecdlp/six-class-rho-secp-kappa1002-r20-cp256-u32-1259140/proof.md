# Six-class secp256k1 rho at R = 2^20: loops unrolled 32 and 16, checkpoint period 256, κ = 1.002 (125.9140)

**Claim.** For a uniform full-production secp256k1 input `Q=[x]P`, the classical cold six-class
(negation + GLV) parallel Pollard-rho attack of `upstream-proof.md`, with the counted kernels and
arrival rows of this document, recovers `x` with probability greater than `0.39` under the
explicit H1 heuristic. A global meter charges every action before it runs, so on every run its total
work is at most

```text
time_log2 = 125.9140 = 1259140/10000
```

target-group-operation equivalents (C = 1000 scalar AArch64 instructions per unit). Peak memory is
below `2^72` bytes; preprocessing is capped at `2^59` units, charged and included; nonuniform advice is
zero. This is an analytic heuristic bound, not a practical break.

## 1. What changes

This package builds on this solver's 199c678c (125.9144: kaankolcu's rare exits and paired traffic from
a9297fb0, κ = 1.002, one batch of all 2^20 walks, the loops unrolled 16 and 8 times, the checkpoint period 256).
Two things change:

- **Pass 1 is unrolled 32 times and pass 2 16 times** (section 2.5). Each round's 31 + 15 tail walks
  (112.125 instructions) are charged with its fixed work: 24,872 per round.
- **The Markov headroom is retuned** to h = 0.000203, p0 = 0.3902033 (section 5.3); the floor stays 0.39000019.

| item | 199c678c | this package |
| --- | ---: | ---: |
| pass 1 / pass 2, per walk | 99.625 / 447.625 | **99.5625 / 447.5625** |
| step / doubling / traversal | 578.7737 / 741.7737 / 611.7737 | **578.6487 / 741.6487 / 611.6487** |
| h, p0 | 0.000212, 0.3902123 | **0.000203, 0.3902033** |
| claim | 125.9144 | **125.9140** |

The ordinary step is **99.5625 + 447.5625 + 31.5 + 24,872/2^20 = 578.6487 instructions**.

## 2. The counted kernels

The two loops with their rare exits (`kernels/secp_loops.s`, section 2.5) are emitted from the emitters of
`kernels/gen_loops.py` (2d25756c's package; the generator is not shipped, and the shipped assembly is what the
experiment executes and counts). The same emitters give the test wrappers `t_mul`, `t_sqr`, `t_sub`,
`t_reduce`, `t_sqr_sub2` and `t_mul_sub` (`kernels/secp_tests.s`).
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

**Fused reduction (24 + 5k).** For x3 and y3 the reduction also takes k canonical subtrahends (k = 2: x_W and
x_T; k = 1: y_W). After step 2 it subtracts each one from all five limbs of s (SUBS, 3 SBCS, then SBC on s4:
5 each), so s' = s − Σ subtrahends with s4' ∈ [−2, C0]. It then folds s4' as a signed word: ASR for its sign,
MUL and SMULH for the 128-bit product w = s4'·C0, and ADDS, 2 ADCS and ADC of (w0, w1, sign, sign) into
s' mod 2^256. That is **12 + 5 + 5k + 1 + 2 + 4 = 24 + 5k**. t = (s' mod 2^256) + s4'·C0 ≡ z − Σ subtrahends
(mod p), and the result is t mod 2^256, which is right unless t leaves [0, p):

- t ≥ 2^256 (s4' ≥ 0). As above, this needs (z − Σ) mod p < C0^2 + C0 and gives (z − Σ) mod p − C0; it is the
  dropped-carry row, with the plain reductions.
- t < 0. This needs s4' < 0 and s' mod 2^256 < |s4'|·C0. With one subtrahend it cannot happen: a borrow from
  s4 = 0 leaves s' mod 2^256 ≥ 2^256 − p > C0. With two it needs s4' = −1 and s' mod 2^256 < C0, probability
  C0/2^256 for x3, and gives a result ≡ z − Σ + C0.
- p ≤ t < 2^256. The result is congruent but not canonical, probability C0/2^256 each for x3 and y3.

The experiment constructs each case. a = 2^128 gives a^2 = 2^256 and s = C0, and then b = p − 1 − e with
c = C0 + 1 + e + d (1 ≤ d ≤ C0) makes the fold go negative. The experiment checks that the kernels return the
predicted values. A pass-2 operation makes one fused reduction for x3 and one for y3; the doubling does the same.
In the loops, the plain reduction's carry into limbs 2 and 3 and y3's signed fold run as rare exits
(section 2.5); the results are the same.

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

### 2.4 fe_sub (9)

**fe_sub(a, b)**, a weak, b canonical, is unchanged: SUBS, 3 SBCS; CSEL of C0 on borrow; SUBS, 2 SBCS, SBC.
On borrow the wrapped value is at least C0 + 1, so subtracting C0 cannot borrow again; the result is weak and
≡ a − b. Without a borrow the result is a − b, which is at least p only if a ≥ p and b < C0. The step uses it
for dx = x_T − x_W (pass 1), dy = y_T − y_W and x_W − x3, whose first operands are canonical, so those results
are canonical; x_W − x3 needs x3 canonical, which holds outside the events above. fe_sub2 of 13c3295f is no
longer used. In the loops its last three instructions (SBCS, SBCS, SBC) are a rare exit (section 2.5).

### 2.5 The two loops (99.5625 and 447.5625 per walk)

`kernels/secp_loops.s` holds each walk's body once, as the assembler macros W1 and W2, and invokes
them (P1 and P2 hold a pair's LDP and its two invocations). Assembled or interpreted, the macros expand to straight-line code; counts below are of the expansion and
of the main path (a rare exit not taken costs its one branch).

**Pass 1 walk** (forward; the schedule of the public package with the new fe_mul): two LDP each for x(T_i) and
x(W_i) (4); dx_i = fe_sub(x(T_i), x(W_i)) (7, canonical); two STP each for acc_(i−1) and dx_i into slot i (4);
acc_i = dx_i·acc_(i−1) with the reduction writing into the four loop-carried registers (84). That is 99, plus the
table pointer, which the body takes as its macro parameter. The main loop runs 32 walks: before each pair,
one LDP x9, x26, [x2], #16 loads both table pointers (the first walk reads x9, the second x26, which no body
writes). After the 32 walks come SUBS and B.NE on a group counter in x29, which the body does not use:
32·99 + 16 + 2 = **3,186 per 32 walks, 99.5625 per walk**.

**Pass 2 walk** (backward, i = n−1 … 1). R ≡ 1/(dx_0···dx_i) stays in x25–x28 across walks: it is loaded once
from Rinv and never stored. Per walk:

| work | instructions |
| --- | ---: |
| 2 LDP acc_(i−1); inv = R·acc_(i−1) | 2 + 84 |
| 2 LDP dx_i; R ← dx_i·R (result into x25–x28) | 2 + 84 |
| 2 LDP y(T_i), 2 LDP y(W_i) (kept to the end); dy = fe_sub(y_T, y_W) | 4 + 7 |
| λ = dy·inv | 84 |
| x3 = λ^2 − x_W − x_T, fused; 2 LDP x(W_i) (kept), x(T_i)'s 2 LDP inside it | 75 + 2 + 2 |
| x_W − x3 = fe_sub(x_W, x3); 2 STP x3 | 7 + 2 |
| y3 = λ·(x_W − x3) − y_W, fused; 2 STP y3 | 90 + 2 |
| **total**: 3·84 + 75 + 90 + 2·7 + 12 loads + 4 stores | **447** |

The main loop runs 16 walks: before each pair, one LDP x0, x5, [x2, #-16]! loads both table pointers. x2
points one entry above the next walk; the first walk reads x5, which its body keeps free, and the second x0,
which is outside the pool. After the 16 walks comes TBNZ x3, #26. The slot pointer x3 moves down by 64 per
walk. The slot array is per-batch scratch, placed once with slot 0 just below a multiple of 2^27 plus 2^26 and
slots 1 to n − 1 above it (64·2^20 = 2^26), so bit 26 of x3 is clear exactly at slot 0, whatever the number of
walks: 16·447 + 8 + 1 = **7,161 per 16 walks, 447.5625 per walk**. Walk 0 uses R itself as its inverse after the
loop. This is the standard affine addition λ = (y_T − y_W)/(x_T − x_W), x3 = λ^2 − x_W − x_T,
y3 = λ(x_W − x3) − y_W, with canonical outputs outside the event of section 2.4.

**Tails.** Each pass first runs its (n − 1) mod 32 (pass 1) or (n − 1) mod 16 (pass 2) left-over walks in a tail
loop: LDR of the table pointer, the body, SUBS and B.NE, 102 and 450 per walk. Then it runs the groups; an
empty group loop is skipped by one CBZ or TBZ. With n = 2^20 every round runs (n − 1) mod 32 = 31 tail walks of
pass 1 and (n − 1) mod 16 = 15 of pass 2, 2.4375 instructions each above the group rate: 112.125 per round, charged
with the round's fixed work (section 4). A doubling changes the tail counts but not their maxima, so it carries no
tail allowance. The setup (LSR and AND for the counts, CBZ for each tail) is per-batch work.

**Rare exits** (the design of kaankolcu's a9297fb0). The generator takes the walk bodies of `gen_loops.py` and
replaces three tails, by text, with an exit. The exit is a branch to a stub `.Lr…` placed after the
function's RET (macros S1 and S2, one copy per walk copy; pass 2 places its tail copy's stubs first, so that every
TBNZ stays within its ±32 KB reach, at most 30,020 bytes). The stub ends with B to the label `.Lb…` just
after the exit. An exit is taken exactly when its tail would change a register:

- **fe_sub.** B.LO after the SUBS of the CSEL'd C0 from limb 0; the stub is SBCS, SBCS, SBC. It is taken iff
  a − b borrowed and the wrapped low limb is below C0. Otherwise C is set, and the tail would subtract 0.
- **Plain reduction.** B.CS after the ADDS, ADCS of (w0, w1) into limbs 0 and 1; the stub is ADCS, ADC of 0
  into limbs 2 and 3. It is taken iff s1 + w1 + c ≥ 2^64, so only if s1 ≥ 2^64 − 2 (w1 ≤ 1, section 2.2). For
  acc and R, step 2 writes limbs 2 and 3 straight into the output registers, whose operand limbs are dead after
  the last product.
- **y3's fused reduction.** TBNZ on bit 63 of s4' before the fold. Its stub is the signed fold of section 2.2
  verbatim (ASR, MUL, SMULH, ADDS, ADCS, ADCS, ADC) and returns after the carry exit. Otherwise s4' ∈ [0, C0],
  so SMULH equals UMULH and the sign word is 0. The main path is then MUL, UMULH, ADDS, ADCS and B.CS to the
  same carry stub as above.
- **x3's fold stays inline.** For a square, s4 = 0 has probability about C0^-1/2 ≈ 2^-16 (u^2 < 1/C0 for
  u = λ/2^256), and then subtracting x_W + x_T makes s4' negative. An exit there would not be rare.

Each stub first runs the macro RCHG: STR x0 pre-indexed, LDR, SUB #20, STR, LDR x0 post-indexed. It takes 20
milli-units from the caller's x19, which the prologue saved at [sp] and the epilogue restores. None of the five
sets flags and x0 is restored, so the tail sees the flags and registers it would have seen inline, after its
charge. A stub is at most 13 instructions (5 + at most 7 + B). There are 191 stubs: each of the 33 pass-1 walk
copies (32 in the main loop, 1 in the tail loop) has 2, and walk 0 has 1; each of the 17 pass-2 copies has 7, and
walk 0 has 5. A walk enters at most 8 stubs (2 in pass 1 and 6 in pass 2; y3's sign stub skips the carry exit),
160 milli-units. The doubling's kernels keep their inline tails.

**Registers.** The pool is x5–x17 and x19–x30 (25 registers). x29 and x30 are saved with x19–x28 by the
prologue and restored before RET; the loops call nothing. x0 holds the tail count, then in pass 2's groups the
second table pointer; x1–x3 hold the three array pointers, x4 holds C0, x18 is never written, x26 holds pass 1's
second table pointer and x29 counts pass 1's groups. The pool's peak is 25 in pass 2 and 16 in pass 1. There
are no spills, and only the stubs touch the stack. Each loop is a single basic block apart from its rare exits:
every other conditional is a carry chain or CSEL, and the only other branch is the back edge.

Counting is the convention of the public package: one LDP or STP is one instruction.

### 2.6 Evidence

`experiment.py` with `"kernels": true` (the declared initial run) does all of the following:

1. **Limb ports.** Python ports of the generated fe_mul (the six chains, as data), the reduction and its fused
   form, fe_sqr and fe_sub assert every no-overflow step above (each carry-only ADC, the exact product identity,
   s4 ≤ C0, w1 ≤ 1, s4' ∈ [−2, C0]) and the predicted value of every rare event of sections 2.2 and 2.4, against
   big integers on all pairs of 44 edge values and on random values, with 4,380 triples through both fused
   forms, 200 constructed carries and 50 constructed negative folds. A 24-walk batch built from the
   ports is checked against affine secp256k1 addition.
2. **The shipped loops, executed.** The AArch64 interpreter executes `kernels/secp_loops.s` on production
   points for batches of 3, 6, 33 and 65 walks (3 and 6 run the tail loops), compares every output with
   affine arithmetic, and reports the executed-instruction difference between 33 and 65 walks per extra walk.
   It must be exactly 99.5625 (pass 1) and 447.5625 (pass 2), and no exit may be taken on these points.
3. **The emitters, executed.** The interpreter executes `kernels/secp_tests.s`, the same emitters wrapped as
   functions, on edge and random inputs: 730 products, 444 squares, 818 fe_sub pairs, 1,670 fused a^2 − b − c
   (40 constructed negative folds) and 1,670 fused a·b − c, and 402 reductions, at least 200 of which carry.
   Every result must equal the port's value under the same rule.
4. **Counts and ledger.** It parses the instruction counts from the shipped files (the stubs after RET are
   not part of a loop), checks the identities 3,186 = 32·(84 + 7 + 4 + 4) + 16 + 2 and
   7,161 = 16·(3·84 + 75 + 90 + 2·7 + 12 + 4) + 8 + 1, and recomputes the doubling bound, the per-round work with
   its tail walks and the ledger of section 5.
5. **The arrival rows, executed.** The interpreter executes `kernels/arrive.s` on production inputs in odd and
   even batches of 1 to 251 walks: random x and y, part A's near ties and near-zero values, x with
   x mod 2^192 > 2^192 − C0, keys forced distinguished (5%) and filter words forced to hit (15%). For every walk
   the stored table address and the rare-path code (logged by the file's test stub) must equal the bit-exact
   model `fast_j`, and every non-fallback index the exact rule. Rare-free batches of 4, 5, 8, 9, 16 and 17 walks
   (walk 0 and the last pair on the odd path) must differ by exactly 31.5 instructions per walk.
6. **Rare exits, executed.** The interpreter runs pass 1 and pass 2 (pass 2 also alone, on chosen slots and
   Rinv) on 15 constructed batches of 12 to 33 walks: second borrows (a < b, (a − b) mod 2^64 < C0); products
   2^128·a whose reduction carries (s4 = 1, s1 = 2^64 − 1), chained by walks that reset acc or R to 2^128;
   λ = 1 with x_W − x3 < y_W (negative y3 fold); and λ = 2^128 with x_W − x3 = a (y3's carry), each kind
   rotated through every walk position. Every output word must equal the limb ports of the inline kernels, x19
   must drop by exactly 20 per exit, and all 191 stubs must be entered: 315 walks, 396 exits, 0 mismatches.

The interpreter (`kernels/a64.py`, which `experiment.py` imports) executes every opcode and addressing form the
shipped files use. It tracks only the C and Z flags, which are all that EQ, NE, HS, LO, CS and CC read, and it
expands assembler macros before running or counting, as the assembler does.

Off the hosted runner (not shipped), clang 23 builds of the shipped files ran in the Unicorn 2.1.4 emulator:
**one production batch of 2^20 walks** and 8,792 walks in 49 smaller batches gave 0 mismatches against affine
arithmetic (exactly 99.5625 and 447.5625 per walk); 2,726 constructed walks took 3,016 exits (all 191 stubs), each
output equal to the inline-tail loops', x19 dropping by 20 per exit; `kernels/arrive.s` ran 16,401 arrivals in odd
and even batches with 0 mismatches, at 31.5 per arrival.

## 3. Arrival rows: 53 → 31.5 plus the round share, executed

The arrival rows are what the arrival pass does per class arrival outside the two counted loops (pseudocode of
section 5.1). 2d25756c priced them as 53 analytic rows. Here they are the loop `kernels/arrive.s`, executed in
the declared run (section 2.6). Each arrival takes 30 instructions on either parity path. Each pair of arrivals
adds one LDP of their K_cp words, one STP of their table addresses and one loop branch (the pairing of
a9297fb0). That is 31.5 per arrival, plus the step's share of the round's fixed work (section 4).

| row | 2d25756c | this | sequence |
| --- | ---: | ---: | --- |
| F1 truncated product | 13 | 13 | 2 LDP of the x limbs; 4 UMULH, 3 ADD, 4 MADD (G_1..G_5 in registers); F0 is the top limb of x |
| canonical index and slot | 21 | **10** | ADD, SUBS, CSEL, CSEL, ADDS, CSEL; CCMN, CCMN, CCMN, CCMP (section 3.1) |
| sign, key, distinguished test | 5 + 1 | **3** | TBNZ on bit 0 of y0; odd y: SUBS, SBCS; even y: CMP, ORR |
| h and table address | 5 | **1** | BFI of bits 0..20 of K0 + 1 into the slot base |
| arrival traffic | 2 | **1.5** | LDP of y0, y1; half an STP of the pair's two table addresses |
| checkpoint filter, rare branch | 3 | **2.5** | half an LDP of the pair's two K_cp limb-0 words; CCMP; B.EQ to the rare path |
| arrival loop | 2 | **0.5** | one TBZ per pair of walks on the pointer array's end bit |
| share of the per-round work | 1 | **0.0237** | 24,872/2^20: the round's fixed work, charged once per round (section 4) |
| **total** | **53** | **31.5237** | |

### 3.1 The canonical index

As in `upstream-proof.md` section 1.1, v0 = x, v1 = βx mod p and v2 = β²x mod p are distinct in [1, p) with sum
p or 2p. Let A, B, C = 2^64·v0/p, 2^64·v1/p, 2^64·v2/p, reals in (0, 2^64) with A + B + C ∈ {2^64, 2^65}. F0
is the top limb of x and F1 the truncated product. In 64-bit words the row forms

```text
S = F0 + F1 mod 2^64        d = F0 - F1 mod 2^64 (SUBS: carry iff F0 >= F1)        m = min(F0, F1) (CSEL)
e = m + S (ADDS: carry iff m + S >= 2^64)
j = 2 if e carries, else 1 if F0 >= F1, else 0
```

With F2* = 2^64 − S ∈ [1, 2^64] (2^64 when S = 0), e carries exactly when m ≥ F2*. The row falls back to the
exact rule (`upstream-proof.md` section 1.1) when any of four windows holds:

```text
(W1) F1 >= 2^64 - 16     (W2) S >= 2^64 - 17     (W3) e mod 2^64 >= 2^64 - 33     (W4) d < 16
```

**One-sided errors.**
- F0 = A − a with 0 ≤ a < (r + C0)/2^192, where r = x mod 2^192. Indeed 2^256/p = 1 + C0/p, so
  A = (x/2^192)(1 + C0/p) = F0 + r/2^192 + x·C0/(2^192·p), and x < p. So a < 2, and a < 1 unless
  r > 2^192 − C0.
- F1 ≡ B − b (mod 2^64) with 0 ≤ b < 16. The true bits 320..383 of x·G are T ≡ F1 + c (mod 2^64) with
  0 ≤ c ≤ 14 (they lie in [F1, F1 + 14], `upstream-proof.md` section 1.1), and T = ⌊B − ε⌋ mod 2^64 with
  0 ≤ ε < 2^-64 (x·G/2^384 = xβ/p − ε/2^64). So b = (B − ⌊B − ε⌋) + c < 2 + 14. If B − b < 0, then
  F1 = B − b + 2^64 ≥ 2^64 − 15, which (W1) excludes. Outside (W1), F1 = B − b as real numbers.
- Then F0 + F1 = A + B − a − b. Since A + B + C ≡ 0 (mod 2^64), F2* ≡ C + a + b (mod 2^64) with
  0 ≤ a + b < 18. If C + a + b > 2^64, then F2* = C + a + b − 2^64 ≤ 17, so S ≥ 2^64 − 17, which (W2) excludes.
  Outside (W1) and (W2), F2* = C + a + b.

So F0 and F1 never exceed A and B, and F2* never falls below C.

**Exactness.** Assume none of (W1) to (W4).
- j = 2 (m ≥ F2*): A = F0 + a ≥ F2* + a = C + 2a + b ≥ C and B = F1 + b ≥ F2* + b = C + a + 2b ≥ C. The v_i
  are distinct, so C is the smallest. No gap test is needed.
- j = 1 (F0 ≥ F1 and F1 < F2*): e = F1 + S < 2^64, so (W3) gives F2* − F1 = 2^64 − e ≥ 34, and (W4) gives
  F0 − F1 ≥ 16. Then A − B = (F0 − F1) + a − b > 16 − 16 = 0 and C − B = (F2* − F1) − a − 2b > 34 − 2 − 32 = 0.
- j = 0 (F0 < F1 and F0 < F2*): (W3) gives F2* − F0 ≥ 34, so C − A = (F2* − F0) − 2a − b > 34 − 4 − 16 > 0.
  And B − A = (F1 − F0) + b − a ≥ 1 − a, which is positive whenever a < 1.

The rule therefore returns argmin(v0, v1, v2) unless j = 0 and r > 2^192 − C0. Under the random-point premise
of the kernel rows, an arrival of phase 1 has r > 2^192 − C0 with probability at most (C0 − 1)·2^64/(p − 1)
< 2^-159.99 (C0 − 1 low parts for each of the 2^64 top limbs). Section 5.3 charges the failure row
(S_ops + 2^61)·C0·2^64/(p − 1) < 2^-33.28. Phase-2 and re-walk operations evaluate the exact rule inside their
600-unit price. Part A checks `fast_j`, the bit-exact model of this rule, against the exact rule on 30,000
random, 10,000 near-tie (k = 150..249, including values near 0 and v0 + v1 near p) and 1,000 low-r inputs,
with 0 mismatches.

**Fallback rate.** (W1) covers 16 values of F1, (W2) 17 of S, (W3) 33 of e and (W4) 16 of d. If βx mod p is
equidistributed (the premise of `upstream-proof.md` section 1.1), each has probability at most its width over
2^64, twice that for (W3) since e is 2F0 + F1 or F0 + 2F1. That is fewer than 2^7 out of 2^64 per arrival, below
upstream's 2^-56, so the fallback price 0.5 and the 2^70 "other" row of section 5.3 still cover fallbacks.

### 3.2 The sequence

Register map of the arrival pass: x1 walk slots, x2 pointer array, x3 K_cp words, x5 and x6 meter, x7–x11
G_1..G_5, x12 p0 + 1, x13 2^64 − 1, x4, x14, x15 the slot bases b_j = base + 64j, x16 the constant 33; x17, x19
and x22–x29 hold one arrival's values. In a pair, the first walk's K_cp word is in x20 and its table address in
x21, the second walk's in x0 and x30; x18 is unused. The constants are set once per batch (section 4).

- **Order.** F1 (13); then ADD S, SUBS d, CSEL m, CSEL of b1 or b0 on the SUBS carry, ADDS e, CSEL of b2 or the
  previous choice on the ADDS carry; then TBNZ on bit 0 of y0, which splits the arrival into two copies that do
  not rejoin before the loop branch, so no merge instruction is needed.
- **Key.** K is the even one of y and p − y (13c3295f). Even y: K0 = y0 and K1 = y1 cost nothing. CMP y1, #0 sets
  Z for the distinguished test, and ORR forms the filter word y0 + 1. Odd y: SUBS from the register holding
  p0 + 1 gives K0 + 1, because y0 is odd and p0 is odd, so y0 ≠ p0 + 1 and the borrow is that of p0 − y0. Then
  SBCS from 2^64 − 1 gives K1, with Z set exactly when K1 = 0. Both paths set flags only after the CSELs have
  used the decision flags.
- **Chain and branch.** CCMN F1, #16, #2, NE; CCMN S, #17, #2, LO; CCMN e, x16, #2, LO; CCMP d, #16, #0, LO;
  CCMP of the filter word with K_cp's limb-0 word, #4, HS; B.EQ. A CCMN test passes with C clear and the CCMP d
  test with C set; when a test is skipped, #2 and #0 carry the failing C forward and #4 sets Z. So a
  distinguished key or a window forces Z, and B.EQ is taken exactly on a distinguished key, a window or a
  filter hit.
- **Address.** BFI writes bits 0..20 of K0 + 1 (odd y) or of y0 = K0 (even y) into bits 8..28 of b_j. Bit 8
  gets bit 0, which is s because K0 is even, and bits 9..28 get bits 1..20 of K0, which are h. b_j has zeros
  there, so the result is base + 512h + 256s + 64j, the table layout of 2d25756c.
- **Rare path.** It is its own action at 0.02 units, charged first: SUBS of 20 milli-units from x5 and B.CC to
  the refill (2), CBZ on K1 to the report (1), the four window tests again and B.LO to the fallback (5). For a
  filter hit it forms K2 and K3 (at most 5: LDP of y2, y3, then for odd y the four-limb SUBS, SBCS, SBCS, SBC
  chain), loads K_cp's other two words (LDP), and runs CMP, CCMP and B.NE back (4). That is at most 17
  instructions. In a pair the handler keeps x0 and x21, delivers its walk's table address in that walk's
  register (x21 or x30), and finds its walk's K_cp entry at a fixed offset from x3, which the LDP immediate
  absorbs. Equal fingerprints are a true repeat and go on to the traversal. The report and the fallback
  are charged their own prices (40.4 and 0.5) before they run. The order of section 5.1 holds: the
  distinguished test first, then the filter.
- **Filter word.** K_cp's limb-0 word now holds K0 + 1, which is odd. Clearing the history stores 0, which is
  even, so a cleared checkpoint never hits. For a different class the word is a pseudo-random odd word, so
  false hits occur at most 2^-63 per arrival.
- **Loop end.** At entry x0 holds n, and TBZ on its bit 0 sends an even n straight to the pairs, which then
  start at walk 0. For odd n walk 0 runs alone, with its own LDR of the K_cp word and STR of its address. Each
  iteration runs a pair of walks. Each walk of the pair is a copy of the sequence above, split by parity, so a
  pair has four tails, and each ends with the pair's STP and TBZ x2, #23 back to the pair's start. The pointer
  array of the batch of n = 2^20 walks is placed once so that &tp[n] is its first address with bit 23 set. Walk
  0's STR advances x2 by 8 and each pair's STP by 16, so for either parity the pairs end exactly there.
- **Test stub.** In `kernels/arrive.s` the rare branch is B.NE over an inline stub that logs the rare case;
  the plain path keeps one branch, as B.EQ does in production.

The other row items of 2d25756c are unchanged:
- **F1 by MADD.** G = ⌊2^384·β/p⌋ (`upstream-proof.md` section 1.1). F1 is the sum, mod 2^64 and with no
  carries, of the low halves of x_i·G_(5−i) and the high halves of x_i·G_(4−i): 4 UMULH and 3 ADD, then 4 MADD
  (xd = xa + xn·xm mod 2^64). G_1..G_5 are loaded once per batch.
- **Anchor rounds and the cap.** Every walk makes exactly one operation (a step, a doubling or a traversal
  step) per batch round. In rounds with t ≡ 0 (mod 256) the processor runs the anchor variant of the arrival
  pass. For every walk it forms K2, K3 and the word K0 + 1 (7), writes K_cp ← fp(K) (STP + STR, the
  checkpoint-write action) and adds one to the trail's count of anchor periods. It abandons the trail when the
  count reaches 20·2^56 (LDR, ADD, STR, CMP, B.HS: 5). A trail's first anchor comes within 256 operations of
  its start, so it is abandoned before L_max = 20·2^64 operations. A walk starting a trail records
  φ = t mod 256 in the trail record, and its anchors fall where (length + φ) ≡ 0 (mod 256). So a re-walk or
  phase-2 continuation anchors where the live trail did. Section 5.3 charges the variant's 12 instructions per
  walk per anchor as a hard row, and the checkpoint writes (2 instructions per anchor) as another.
- **Zero denominator.** See section 4.

## 4. Prices and per-batch overhead

- **Ordinary step:** 99.5625 + 447.5625 + 31.5 + 24,872/2^20 = **578.625 + 24,872/2^20 = 578.6487**
  instructions, priced as that exact rational.
- **Traversal step:** step + 24 + 2 + 7 = **611.625 + 24,872/2^20**: the 24 instructions to track and resume the minimum-fp
  member (unchanged), the charge of its surcharge on its own path (SUBS of 33 milli-units from x5, B.CC), and
  the 7 that form K2 and K3 for the fingerprint the minimum compares.
- **Doubling:** the batched doubling shares the product chain and inversion. Its arithmetic is
  4 fe_mul (340) + fe_sqr of x (64) + dx = 2y by one weak addition (9) + 3x^2 by two weak additions with a
  second carry fold (2·14) + the fused x3 = λ^2 − x − x (75) + x − x3 (9) + the fused y3 = λ·(x − x3) − y (91)
  = 616, with the inline kernels and the rare events of section 2.
  - **Traffic: 30** (the public package allowed 72). The loads, stores, pointer updates and loop control are
    at most the addition loops' 5 + 4 + 2 and 13 + 4 + 2. A doubling walk reads its own y for dx = 2y in the
    forward pass and acc, dx, x and y in the backward pass, with no table point. It stores the same dx, acc,
    x3 and y3, and shares R in registers.
  - **List moves: 34** (the public package counted 50). Moving the walk to and from the doubling list is
    the upstream 16 loads and 16 stores of 64-bit words, the walk's x and y out and back. They pair into
    8 LDP and 8 STP, plus 18 other instructions.
  - **Escape logic: 30**, as upstream. With checkpoint-only detection and the meter per batch it is now the
    clearing of the history (one STR of 0 into K_cp's limb-0 word) and the charge of the doubling's surcharge on
    its own path (SUBS of 163 milli-units from x5, B.CC), 3 instructions, well inside the 30.

  - **No tail allowance.** If d walks double in a round, the addition loops run n − d walks, whose tails are at
    most 31 and 15 walks: the most that the round's fixed work already pays for (section 2.5).

  The doubling is 616 + 30 + 31.5 + 24,872/2^20 + 34 + 30 = 741.5 + 24,872/2^20, priced at
  **741.625 + 24,872/2^20** so that its surcharge over the step is a whole 163 milli-units. The weak additions, the 18 and the 30 are those of
  `upstream-proof.md` section 5.3.
- **Zero-denominator test, once per batch.** Pass 1 multiplies the n canonical dx_i into acc_(n−1); p is
  prime, so outside the rare event of section 2.2 acc_(n−1) ≡ 0 (mod p) exactly when some x(T_i) = x(W_i).
  Before the inversion acc_(n−1) is canonicalised and tested (ADDS, 3 ADCS, 4 CSEL, 3 ORR, CBZ: 12 per batch).
  Only when it fires, with probability < 2^-114 over the whole run, does the batch take the cold path: pass 1
  is re-run with the per-walk test, each walk with dx_i = 0 is abandoned and its dx_i replaced by 1, as in the
  public package. The re-run, at most 2^20·(102 + 4) < 2^27 main-path instructions, is its own action at 2^17
  units in the "other" random row; its stubs charge themselves as in section 2.5.
- **Per-batch fixed work:** inversion I ≤ 24,322 (unchanged `secp_fe.s`, `upstream-proof.md` section 5.1–5.2),
  prologue, walk 0, setup and epilogue of the two passes (3,323 − 102 − 3,186 + 7,912 − 450 − 7,161 = 336), 64
  for the call and return, 12 for the zero test, 3 for the constants p0 + 1 and 2^64 − 1, 3 for loading
  G_1..G_5, 7 for the meter, 3 for the round counter and the anchor-round test (ADD, TST, B.EQ), 4 for the slot
  bases b1 and b2, the constant 33 and the pointer array's base, and 5 for the arrival's parity test, walk 0
  and its exit: 24,759 < 2^15.
  All 2^20 walks form one batch, unlike the batches of 2^15 to 2^16 walks of `upstream-proof.md` section 3:
  every walk makes one operation per round, so the batch never changes. With n = 2^20 the tail walks add
  31·2.4375 + 15·2.4375 = 112.125 every round (section 2.5), so a round's fixed work is at most 24,871.125
  instructions, charged as 24,872. The batch's scratch is 72 bytes per walk (slots and pointers), 72 MiB in all, inside the
  memory bound of `upstream-proof.md` section 5.3.
- **Meter per round** (the reserve of a9297fb0, the fixed work charged per round). The processor keeps the
  round charge c = 834.625·2^20 + 24,872 milli-units: 578.625 per walk, a reserve of 256 per walk, and the round's fixed
  work. Before a round runs it charges c from its block remainder x5 and adds 2^20 to its operation count x6
  (LDR, SUBS, B.CC to the refill, ADD), returns the reserve (ADD x5, x5, n, LSL #8) and holds the remainder in
  the callee-saved x19 across the passes (two MOV): 7 instructions per round, inside the 24,872. A round is thus
  charged (578.625 + 24,872/2^20)·2^20 milli-units, the step price for each of its walks. x19 ≥ 256·n when pass 1
  starts. The stubs of section 2.5 take their 20 milli-units each from it: at most 8 per walk, or 10 if the zero
  test re-runs pass 1, so at most 200 per walk, and x19 stays non-negative without a test. A round's charge is
  below 2^30 milli-units, far below a 2^40-unit block, so one refill always covers it. Every walk of the round
  is covered at the step price before its operation runs. A doubling or traversal walk pays its surcharge (163
  or 33 milli-units) on its own path, and a rare-path entry or a rare exit pays its own price the same way, so
  every action is charged before it runs. Blocks are reserved, refilled and counted as spent exactly as in
  `upstream-proof.md` section 3, and phase 1 still stops at a refill once the operation count reaches S_ops.
  The 2^20 processors' last blocks cover fewer than 2^20·2^40/0.579 = 2^60.79 further operations, so phase 1
  performs fewer than S_ops + 2^61 operations, the count that section 3 and section 5.3 use.

Every other price (reports 40.4, probes 0.05, starts 2^13, fallbacks 0.5, phase-2 and re-walk operations 600,
refills 0.03, words, attempts, store initialisation) is unchanged; the checkpoint write (0.002), the rare-path
entry and the rare exit (0.02 each) are their own actions; the anchor-round work is a hard row (section 3).

## 5. The walk at R = 2^20, H1 and the exact budget

### 5.1 The walk change

The walk is that of `upstream-proof.md` section 2 with five changes: the table has R = 2^20 rows instead of
2048, cycles are detected by the checkpoint alone, the class key is the even one of y and p − y, the checkpoint
period is 256 instead of 64, and anchors fall where (length + φ) ≡ 0 (mod 256). Everything else is unchanged:
the canonical index, the twisted table Tab[h][t] = [u_t^{-1}]M_h, the distinguished classes (limb 1 of K zero,
θ = 2^-64), M = 2^20 walks and L_max = 20·2^64. At each class arrival, in order:

```
if K distinguished: report (K, a, b); end trail
if length >= L_max: abandon trail
if fp(K) == K_cp: traverse the cycle once more (<= 256 steps, else abandon), resume at the member of its
     minimum-fp class remembered during the traversal, ESC; clear history; continue
if (length + phi) mod 256 == 0: K_cp <- fp(K)   (phi: the trail's start round mod 256, section 3)
P <- P + Tab[h][t]
```

"Clear history" stores 0 in K_cp's limb-0 word (section 3.2). ESC is the doubling escape P ← [2]P, a class map.

- **Merging lemma.** Let trails A and B arrive at the same class C, with any histories. Until one of them
  escapes, both apply F and visit the same classes. A trail escapes only after K = K_cp. K_cp is set at an
  arrival, and the arrivals up to the next refresh are F-steps. So a repeat of K_cp is a cycle Z of F of length
  ≤ 256, and Z contains no distinguished class. Every escape from Z leaves from its minimum-fp class and goes to
  E(Z) = class(2c(min Z)), which depends only on Z. If A escapes from Z, B follows the same F path and enters
  Z. B sets K_cp inside Z within 256 arrivals, returns to it after |Z| ≤ 256 more, and escapes to the same E(Z)
  after a traversal of |Z| and one doubling. So A and B are again at a common class after at most
  256 + 256 + 256 + 1 = 769 extra operations of one of them. By induction they reach the same distinguished
  class, unless one of them first hits L_max or an exceptional addition (failure events). History is cleared
  at every doubling and starts have none, so K = K_cp always witnesses a real cycle of F.
- **Cycle rates.**
  - A 2-cycle needs the next class to use the same partition index with the inverse twist: rate 1/(6R) =
    2^-22.585 per step.
  - Other cycles of length 3 to 256 need every partition index used at least twice: rate at most
    Σ_L Σ_b S2(L, b)·R^-(L−b) = 2^-37.99999 at R = 2^20, the sum of `upstream-proof.md` section 4 taken to
    L = 256. To L = 64 it is 2^-19.9954 at R = 2048, the public value; `certificate.py` recomputes both exactly.
  - Cycles longer than 256, which the checkpoint does not catch, are rarer still: below 2^-1726 per step for
    L = 257 to 512 at R = 2^20, with the terms still falling.
  - A cycle in which some partition index is used only once closes by chance: a fresh step lands on one of
    the trail's own t earlier classes, with probability at most t/N under H1's pseudo-random key fields, at
    most L_max^2/(2N) = 2^-117.8 per trail (H1's non-fruitless-return clause) and under 2^-54 per run; their
    doublings and traversals cost under 2^-47 units, inside the 2^70 "other" row of section 5.3. On the toy
    group such cycles occur at about (trail age)/N per step.
  - So the expected number of doublings is below (1/(6R) + 2^-37.99999)·(S_ops + 2^60) + 2^-54
    < 2^-22.5849·(S_ops + 2^60). Every 2-cycle is traversed in 2 steps and other short cycles in at most 256,
    so the expected number of traversal steps is below (2/(6R) + 256·2^-37.99999)·(S_ops + 2^60) + 256·2^-54
    < 2^-21.5807·(S_ops + 2^60). Only linearity of expectation over the pseudo-random partition and twist of
    each new class is used.
- **Detection delay.** A cycle is caught within 255 + |Z| ordinary steps of entering it. That is ordinary work:
  it is counted in S_ops, and its arrivals are counted in T_c. At R = 2^20 it is below
  257·2^-22.585 + 511·2^-38 = 2^-14.579 of all operations; with the traversals, the waste is w ≤ 2^-14.56.
- **Table.** 2^20 rows of 512 bytes (section 3), 2^29 bytes, aligned to 2^29. Even with one copy on each of
  the 2^16 nodes, each node stays below 2^56 bytes and peak memory below 2^72 bytes (`upstream-proof.md` 5.3,
  "Memory"). Building it takes 2^20 draws of [c]P + [d]Q at ≤ 2^13 units each plus twists: below 2^33.1 units,
  in the hard extras (< 2^34) and in preprocessing, below its 2^59 cap.
- **Key and partition.** The class key is the even one of y and p − y, and h is bits 1..20 of K (section 3),
  disjoint from the distinguished bits 64..127, the routing bits 128..143 and the bucket bits 144..193. The
  anchor phase φ only shifts where each walk's anchors fall; every walk is still anchored once per 256 of its
  own operations, which is all the merging lemma and the detection-delay bound use.

### 5.2 H1 at R = 2^20

H1 keeps its form: P(T_c ≤ S_ops) ≥ 1 − exp(−(S_ops/κ)^2/(2N')), with N' = N·R/(R − 1) and now κ = 1.002, on
the hazard census of 12ccd81d, which this package ships (`census_hz.py`, results in `census-hz-q29.json`). The walk is
unchanged by the rare exits, which change no result.
T_c counts every class arrival, including repeated arrivals inside fruitless cycles, during the detection
delay and in traversals. Its other parts are as in `upstream-proof.md` section 4: pseudo-random key fields, the
non-fruitless return bound, and a uniform twist index independent of the partition index.

**Why T_c runs slow on a toy.** T_c ends at the first arrival at a class visited by *another* trail, so a
trail's own classes are never collision targets.
- At an arrival with D classes visited, `own` of them by the current trail, the random-mapping chance of a
  collision is (D − own)/N, not D/N.
- The first collision therefore comes about one mean trail age ā = (1 − θ)/θ later: E[T_c] ≈ √(πN/2) + ā.
- On q29 with dp_bits 5 that is ā/√(πN/2) = +0.45 % at the mean. This is the slowdown that 02feb284's census
  (+0.70 % at R = 2048) and this solver's censuses (+0.34 % to +0.77 %, 2d25756c) measured, and that κ = 1.01
  was set to cover.
- In production it is ā/√N = 2^-62.7, and at most S_ops·L_max/N = 2^-58.4 of the hazard.

**Hazard census.** `census_hz.py` runs the walk of section 5.1 on the toy groups (the walk code of 2d25756c's
`census_cp.py`).
- At every fresh-sourced arrival it adds the random-mapping hazard (D − own)/N. A fresh-sourced arrival is a
  trail start, an escape, or a step from a class first visited at the previous arrival.
- A merge is an arrival at another trail's class.
- For a uniform random mapping the hazards are the merges' compensator, so E[merges] = E[Σ hazards]. So
  ρ = merges/Σ hazards measures the walk's collision intensity against the law, with the own-trail exclusion
  removed exactly.
- Repeated arrivals add no hazard: they are the waste w.
- In multi mode a run of 64·√N arrivals continues after each merge, which ends only the merging trail.
- In tc mode each trial stops at its first merge (T_c).
- `--pure` runs uniform random functions through the same trail logic: the control.

The table gives ρ and κ_eff = ρ^-1/2 from `census-hz-q29.json` (q29: N = 2^24.83; q25: N = 2^22.42). The error
is the larger of the Poisson and table-spread standard errors.

| run | tables | merges | ρ | κ_eff |
| --- | ---: | ---: | ---: | ---: |
| walk, q29, R = 2^20, dp_bits 5 | 64 | 3,262,687 | 1.00000 ± 0.00068 | 1.00000 ± 0.00034 |
| walk, q29, R = 2^20, dp_bits 3 | 16 | 816,680 | 1.00099 ± 0.00135 | 0.99951 ± 0.00068 |
| walk, q29, R = 2^20, dp_bits 7 | 16 | 815,474 | 1.00109 ± 0.00201 | 0.99945 ± 0.00100 |
| walk, q29, R = 2^16, dp_bits 5 | 16 | 813,269 | 0.99783 ± 0.00159 | 1.00109 ± 0.00080 |
| walk, q29, R = 2048, dp_bits 5 | 16 | 800,858 | 1.00100 ± 0.00199 | 0.99950 ± 0.00099 |
| walk, q25, R = 2^16, dp_bits 5 | 16 | 811,399 | 1.00202 ± 0.00186 | 0.99899 ± 0.00093 |
| walk, q25, R = 2^20, dp_bits 5 | 16 | 810,856 | 0.99961 ± 0.00162 | 1.00019 ± 0.00081 |
| ideal mapping, q29 N, dp_bits 3 | 16 | 815,722 | 0.99972 ± 0.00111 | 1.00014 ± 0.00055 |
| ideal mapping, q29 N, dp_bits 5 | 16 | 815,230 | 0.99956 ± 0.00179 | 1.00022 ± 0.00090 |
| ideal mapping, q29 N, dp_bits 7 | 16 | 818,029 | 1.00389 ± 0.00267 | 0.99806 ± 0.00133 |

- Pooled over the three q29 runs at R = 2^20 (dp_bits 3, 5, 7): ρ = 1.00027 ± 0.00058, κ_eff = 0.99986 ± 0.00029
  (χ² = 0.6 on 2 degrees of freedom).
- Over the seven multi-mode runs of the walk: ρ = 1.00015 ± 0.00048 (χ² = 4.1 on 6).
- The ideal mapping: ρ = 1.00014 ± 0.00089 (χ² = 2.2 on 2), so the estimator is unbiased at this precision.

**The control at T_c.** In tc mode (32,000 trials each), the raw T_c slowdown against √(πN'/2) and the hazard
at the first merge are:

| dp_bits | own-trail shift ā/√(πN'/2) | walk, R = 2^20 | ideal mapping | walk ρ | P[H ≤ ln(1/0.61)] (law 0.39) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | +0.10 % | +0.07 % ± 0.29 % | +0.61 % ± 0.30 % | 0.9974 ± 0.0068 | 0.3936 (walk), 0.3912 (ideal) |
| 5 | +0.45 % | +0.67 % ± 0.29 % | +0.11 % ± 0.29 % | 0.9973 ± 0.0056 | 0.3884 (walk), 0.3934 (ideal) |
| 7 | +1.85 % | +1.57 % ± 0.29 % | +1.78 % ± 0.29 % | 1.0104 ± 0.0075 | 0.3937 (walk), 0.3910 (ideal) |

Both raw slowdowns follow the predicted own-trail shift within 1.8 standard errors and grow with the trail
length, while the walk's ρ stays within 1.5 standard errors of 1.

**First-collision distribution.** H1 needs P(T_c ≤ S_ops), not only the mean intensity. For a random mapping
the hazard H accumulated up to the first merge is Exp(1), so P(T_c ≤ S) = P(H ≤ hazard accumulated by S). The
tc runs record H at each trial's first merge (`census-hz-q29.json`). For the walk the mean of H is 1.0026,
1.0027 and 0.9897 at dp_bits 3, 5 and 7. The six fractions with H ≤ ln(1/(1 − 0.39)) in the table lie within
1.4 binomial standard errors (0.0027) of the law's 0.39, five of them above it. So the walk's first collision
follows the random-mapping law at the quantile H1 uses, as well as its mean.

**κ = 1.002.**
- In production a run of S arrivals has at least S/(1 + w) fresh ones, with w ≤ 2^-14.56 (section 5.1:
  detection delay, cycle repeats and traversals).
- Its own-trail term is below 2^-58.4 of the hazard, and its collision intensity is ρ.
- So P(T_c ≤ S) ≈ 1 − exp(−ρ·(S/(1 + w))^2/(2N)), which is H1 with κ_eff = (1 + w)/√ρ.
- κ rests on the multi-mode runs. They give κ_eff = 1.00000 ± 0.00034 at R = 2^20 with dp_bits 5, and
  0.99992 ± 0.00024 pooled over the seven walk runs. κ = 1.002 is 5.9 standard errors above the first and
  above each multi-mode run's κ_eff plus one standard error: the largest is 1.00109 + 0.00080 = 1.00189 (walk,
  q29, R = 2^16, dp_bits 5), and for the ideal mapping 1.00022 + 0.00090 = 1.00112.
- The tc runs have 25 to 102 times fewer merges and measure the first-collision law rather than κ. The walk's
  tc values are κ_eff = 1.0013 ± 0.0034, 1.0013 ± 0.0028 and 0.9949 ± 0.0037 at dp_bits 3, 5 and 7, all
  below 1.002. For the ideal mapping, whose κ_eff is exactly 1, they are 1.0067 ± 0.0028, 0.9974 ± 0.0028 and
  1.0000 ± 0.0038, which shows the tc estimator's spread.
- N' = N·R/(R − 1) is kept. The twisted walk's 1/(6R) anti-collision is inside ρ, so N' over-covers it.
- This is the standard random-mapping heuristic, supported at toy scale; no theorem bounds the constant.
  Brent and Pollard put an r-adding walk's excess at about 1/(2r), at most 2^-21 here (r ≥ R = 2^20).

**Other checks**:
- **Waste.** The census counts repeated arrivals: w = 1.31·10^-4 per arrival at q29, R = 2^20, dp_bits 5,
  mostly the detection delay after toy self-meetings; section 5.1 bounds production waste by 2^-14.56.
- **Trap test** (`experiment.py`'s `long_trails` mode, q29, 150 trails of 2,048 operations, seed 7):
  non-fruitless returns in 0/150 trails at R = 4 and at R = 8 with the doubling escape, against about 10/150
  from ordinary self-collisions; the additive control gives 8/150 and 6/150.
- **The hosted toy pipeline** draws a fresh table per trial at R = 2048 and adds the toy's own-trail term ā to
  its scaled budget, so that its success rate measures the production law. It is evidence of the whole
  procedure, not of κ.

### 5.3 The exact budget

H1 and κ = 1.002 remain the disclosed heuristic, supported only at toy scale (section 5.2). The Markov headroom
is h = 0.000203 and the budget for every other failure event 1.1·10^-7, so p0 = 0.3902033 keeps the success floor
at 0.39000019 > 0.39 with W below 2^125.9140 (the retune of aburan28's 4b4defe7, which lowered the floor to
0.3900005 at the same structure).

```text
N'     = ((r-1)/6) * 2^20/(2^20 - 1)
S0     = sqrt(2*N'*ln(1/(1-0.3902033)))      = 2^126.6997298
S_ops  = 1.002*S0                            = 2^126.7026123
phase1 = (578.625 + 24872/2^20)/1000*(S_int + 2^60) = 2^125.913372
hard extras                                  = 2^112.578215
   anchor rounds 0.012*ceil((S_int + 2^61)/256)           2^112.321812
   checkpoint writes 0.002*ceil((S_int + 2^61)/256)       2^109.736850
   upstream rows (table < 2^34, refills at the new cap)    2^107.155501
expected random extras                       = 2^102.1176
   doublings  0.163*(1/(6R) + 2^-37.99999)*(S_ops + 2^60)              2^101.50
   traversals 0.033*(2/(6R) + 256*2^-37.99999)*(S_ops + 2^60)          2^100.20
   rare-path entries 0.02*(1/(6R) + 2^-37.99999 + 2^-56)*(S_ops + 2^60) 2^98.47
   rare exits 0.02*(3*C0/2^64 + 5*2^-63 + 2^-27.4)*(S_ops + 2^61)       2^93.83
   starts 2^13*(S_ops/2^64 + 2^20)                                     2^75.70
   words, attempts, probes, small batches, fallbacks, zero-test re-runs                  2^70
allowance A = ceil(extras/0.000203)          = 2^114.3839
W_ledger = phase1 + hard + A                 = 2^125.9139995
```

The rare-path row counts filter hits (true ones at most at the doubling rate, false ones below 2^-63), window
fallbacks (below 2^-57) and distinguished keys (2^-64) per arrival. The rare-exit row counts, per operation:
- the three fe_sub exits, each taken only when the low limb of a − b wraps below C0: at most C0/2^64 each;
- the five carry exits (four plain reductions and y3's), each taken only when s1 ≥ 2^64 − 2: at most 2^-63
  each;
- y3's sign exit, which needs s4 = 0 for the product λ·(x_W − x3) of two pseudo-random field elements: at most
  2^-32·(1 + 32 ln 2) < 2^-27.46, rounded up to 2^-27.4.

So the expected number of exits is at most 3·C0/2^64 + 5·2^-63 + 2^-27.4 < 2^-27.2317 per operation, at
0.02 units each, over fewer than S_ops + 2^61 operations. The rates rest on the random-point premise of the
kernel rows. Every exit is charged before its tail runs (section 4) whatever the inputs, so the premise affects
only this expectation. `certificate.py` uses the exact rates; `experiment.py` rounds the cycle rates up to
2^-22.5849 and 2^-21.5807. Both give W_ledger = 2^125.9139995.

In exact milli-units (`certificate.py`):

```text
S_int             = 138448043220429305728042121759826770979
phase 1           = 80112782966599406743066895997242804564382
anchor rounds     = 6489752025957623706110060848548771804
checkpoint writes = 1081625337659603951018343474758128634
hard (incl. both) = 7752103160058692404790705059598979866
extras            = 5501341499053279584328988563120120
A                 = 27100204428833889577975313118818325124
W_ledger          = 80147635274188299325049662015421221869372
```

W_ledger = 2^125.9139995 < 2^125.9140. The next lower four-decimal cap, 125.9139, is below W_ledger and is not
claimed. Against the cap's headroom (cap − phase 1 − hard = 2^114.385438), the Markov ratio is
0.000202778 (exact; experiment.py's float ledger: 0.000202782) < 0.000203, so the meter halts early with probability below h. The 1.1·10^-7 covers every other
failure event: the rows of `upstream-proof.md` section 4 (trail cap < 10^-7, preprocessing cap < 2^-30 and
smaller ones), the three kernel rows of section 2 and the canonical-index row of section 3.1. These are the
dropped carry 7·(S_ops + 2^61)·(C0^2 + C0)/p < 2^-62.4, at most 7 reductions per batched operation, plain or
fused; the negative fused fold (S_ops + 2^61)·C0/2^256 < 2^-97.2; a fused result in [p, 2^256) for x3 or y3,
2·(S_ops + 2^61)·C0/2^256 < 2^-96.2; and a wrong j = 0, (S_ops + 2^61)·C0·2^64/(p − 1) < 2^-33.28. The rare
exits return the inline results, so they add no failure row. Together these sum to 1.0103·10^-7. So

```text
success >= 0.3902033 - 0.000203 - 1.1e-7 = 0.39000019 > 0.39.
```

`certificate.py` checks everything in exact integer and rational arithmetic: the kernel identities (85, 64, 91,
75, and 84, 7 and 90 with the rare exits), the loops' 3,186 and 7,161, the tails 112.125, the rows, the
per-round work 24,872 and the prices 578.625, 741.625 and 611.625 plus 24,872/2^20, the cycle sums, the rare-exit rate,
every ledger row, the four failure rows and their sum with the upstream rows, and
W^10000 < 1000^10000·2^1259140 and W^10000 > 1000^10000·2^1259139. It prints `CERTIFICATE: PASS`. The counter
widths, preprocessing cap and memory are re-checked there.

## 6. Provenance, files and limitations

Design provenance: the construction, report store, ledger structure and F1 truncated product are the public
six-class line reproduced in `upstream-proof.md`, carried by 02feb284 (126.3059). This solver's packages supply
the counted kernels, loops, row cuts, doubling itemization and checkpoint action (daa13d54, 09990b8a, 1ca8ae82,
5d3bf6c4); the R = 2^20 checkpoint-only walk and h, p0 (f99114e0, after the accepted BLS12-381 G1 package
67154fda); the rare-event kernels and their failure rows (736c001a); F1 by MADD, the meter per batch, anchors in
global rounds and the even-y key (13c3295f, after the accepted P-256 package 6d5e2adb); the fused
multiply-subtract-reduce (2d25756c; fusing x3's subtractions first appeared as fe_sub2 in the accepted
b5de673a); the canonical-index rule and executed arrival rows (67772ed0); the unrolled loops (d50309e4); the
hazard census with κ = 1.003 (12ccd81d); and their combination with kaankolcu's rare exits, meter reserve and
paired traffic (a9297fb0, after their 807e4bf5) in d23427fc (125.9194); κ = 1.002, the fixed work charged once
per batch and the headroom retune after aburan28's 4b4defe7 (after their b78fb4f6) in 5ffb2d5a (125.9168); one
batch of all 2^20 walks in ae903c5c (125.9154); the loops unrolled 16 and 8 times in 7f077a55 (125.9148); the
checkpoint period 256 in 199c678c (125.9144), the base of this package. All but 09990b8a, 1ca8ae82, 5d3bf6c4 and
d50309e4 are declared dependencies. New here: the loops unrolled 32 and 16 times, and h = 0.000203. No prior acceptance, review outcome or score is a premise.

- `proof.md`: this argument. `upstream-proof.md`: the construction (abridged; its kernel, row, price, walk and
  ledger numbers are superseded here).
- `experiment.py`: part A production checks, the toy pipeline (the checkpoint-only walk of section 5.1) and the
  kernels mode of section 2.6. `certificate.py`: the exact certificate.
- `census_hz.py`, `census-hz-q29.json`: the hazard census of section 5.2 and its results, at period 256.
- `kernels/secp_loops.s` (counted loops, unrolled, with rare exits); `kernels/secp_tests.s` (emitter wrappers)
  and `kernels/secp_fe.s` (the inversion kernels), unchanged except for operand spacing;
  `kernels/arrive.s` (the arrival rows by pairs, with an inline test stub on the rare path); `kernels/a64.py`
  (the interpreter, counter and rare-exit check).

The limitations are the upstream ones. H1 is a random-mapping heuristic supported at toy scale; κ = 1.002 rests
on the hazard census of toy groups (N up to 2^24.83). The step price is the executed main-path count of the
shipped loops and arrival rows plus the per-batch work divided by the smallest batch; the rare exits are charged as their own actions, at rates
that rest on the random-point premise. The kernel and canonical-index failure rows rest on the same premise as
the exceptional-addition row. The attack is generic and impractical, and the bound is not a claim about
deployed secp256k1.

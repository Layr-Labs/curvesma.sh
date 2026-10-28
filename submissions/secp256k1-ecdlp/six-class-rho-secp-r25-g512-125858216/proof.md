# Six-class secp256k1 rho, R = 2^25, κ = 1.0001, fused pass 1, walk groups of 512 (125.858216)

**Claim.** For a uniform full-production secp256k1 input `Q=[x]P`, the classical cold six-class (negation + GLV) parallel Pollard-rho attack of `upstream-proof.md`, with the counted kernels and
arrival rows of this document, recovers `x` with probability greater than `0.39` under the explicit H1 heuristic. A global meter charges every action before it runs, so on every run its total work is
at most

```text
time_log2 = 125.858216 = 125858216/1000000
```

target-group-operation equivalents (C = 1000 scalar AArch64 instructions per unit). Peak memory is below `2^72` bytes; preprocessing is capped at `2^59` units, charged and included; nonuniform advice
is zero. This is an analytic heuristic bound, not a practical break.

## 1. What changes

This package builds on blocksek's da765727 (125.858233), which doubled to 256 the walk groups of bradcush's 886e6f54 (one address rule for both y parities, section 3.2, and the
shorter plain fused pass 1); section 6 gives the line. It keeps that package's walk, census, κ, h, p0, kernels and meter, and changes one parameter:

-  **Walk groups of 512** (LSR #9, AND #511) in pass 1, both fused copies and pass 2, so each group's closing instructions (SUBS and B.NE; in pass1f LDR, SUBS, STR and B.NE) serve twice the
  walks: pass1p is 125.50390625 per walk, pass 2 432.50390625, pass1f 129.0078125 (compare rounds 3.50390625 more, a hard row). Each round runs 511 + 511 tail walks, in its fixed work 28,112.
  The step is 558.0078125 + 28,112/2^20 = 558.034622 (was 558.040848), and W_ledger = 2^125.858215915 is claimed as 125.858216.
-  The fused tail copy's labels take id 999 (512 is a group walk's), the experiment's batches grow with the groups, and the declared toy budget is 7 s, so the declared run stays under the 20 s run limit (`evidence.md`).

## 2. The counted kernels

The two loops with their rare exits (`kernels/secp_loops.s`, section 2.5) are emitted from the emitters of `kernels/gen_loops.py`, with the x3 and y3' forms of section 2.2 (2d25756c's package; the
generator is not shipped; the experiment executes and counts the shipped assembly). The same emitters give the test wrappers `t_mul`, `t_sqr`, `t_sub`, `t_reduce`, `t_sqr_sub2` and `t_mul_sub`
(`kernels/secp_tests.s`). p = 2^256 − C0, C0 = 2^32 + 977, so 2^256 ≡ C0 (mod p). A field element is four little-endian 64-bit limbs; **weak** means [0, 2^256), **canonical** [0, p). Walk and table
points are canonical, except after the rare event of section 2.4. x4 holds C0 in the loops.

### 2.1 fe_mul: 32 products and six carry chains (62)

Write `(Hij, Lij)` for the high and low words of `a_i·b_j`. `Lij` sits in column `i+j` and `Hij` in column `i+j+1` (column k has weight 2^(64k)). The columns 0 to 7 hold 1, 3, 5, 7, 7, 5, 3, 1 words.
Reducing a column of t words to one word takes t − 1 two-word additions, 24 in all. An ADDS/ADCS chain makes one addition per column and hands its carry to the next column inside the same chain; MUL
and UMULH between its instructions leave the flags alone. Every chain therefore ends with one addition of its last carry. Columns 3 and 4 need six additions each, so at least six chains are needed;
this schedule uses six (30 additions):

| chain | col 1 | col 2 | col 3 | col 4 | col 5 | col 6 | col 7 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | H00+L01 | H01+L02 | H02+L03 | H03+L13 | H13+c | | |
| 2 | +L10 = z1 | +H10 | +H20 | +H30 | H22+c | | |
| 3 | | +L20 | +L30 | +L22 | H13'+H22' | H23+c | |
| 4 | | +L11 = z2 | +H11 | +H12 | +H31 | H32+c | |
| 5 | | | +L12 | +H21 | +L23 | H23'+H32' | H33+c |
| 6 | | | +L21 = z3 | +L31 = z4 | +L32 = z5 | +L33 = z6 | +c = z7 |

z0 = L00. "+X" adds word X to the column's running sum; "H+c" adds only the incoming carry to a high word. The chains run in the order 1 to 6, and every word is used exactly once (the experiment's
port asserts it).

-  **No carry is lost.** A high word is at most ⌊(2^64 − 1)^2 / 2^64⌋ = 2^64 − 2, and chains 1 to 5 end by adding one carry bit to a high word no earlier instruction has touched (H13, H22, H23, H32,
  H33), so those ADCs cannot carry out. Chain 6 ends by forming z7 after z0 to z6 are final; the eight limbs then satisfy Σ z_k·2^(64k) = a·b < 2^512, so z7 < 2^64. Every other addition's carry is
  consumed by the next instruction of its chain. Hence (z0, …, z7) = a·b exactly.
-  Each product is formed just before its addition and each operand limb freed after its last product, so fe_mul's register peak is 15, its 8 operand limbs included.

### 2.2 The reduction (23)

z = L + H·2^256 with L, H < 2^256 (the eight product limbs).

1.  **H·C0** (8 + 4): four MUL and four UMULH give g_j, k_j with g_j + 2^64·k_j = z_(4+j)·C0, k_j ≤ C0 − 1. ADDS g1+k0, ADCS g2+k1, ADCS g3+k2, ADC k3+c gives H·C0 in five limbs; H·C0 < 2^256·C0, so
   the ADC cannot overflow.
2.  **s = L + H·C0** (5): ADDS, 3 ADCS, then ADC into the fifth limb s4. Since s ≤ (2^256 − 1)(1 + C0) < (1 + C0)·2^256, s4 ≤ C0.
3.  **w = s4·C0** (2): MUL and UMULH; w ≤ C0^2 < 2^65.
4.  **t = (s mod 2^256) + w** (4): ADDS, 2 ADCS and ADC, which return t mod 2^256. t ≡ z (mod p) because 2^256 ≡ C0, and t < 2^256 + C0^2.

**23 = 8 + 4 + 5 + 2 + 4.** There is no final fold. If t < 2^256 the result is t, a weak value ≡ z (mod p). If t ≥ 2^256, then C0 ≤ t − p < C0^2 + C0 and t − p ≡ z, so z mod p = t − p and the result t
− 2^256 = z mod p − C0 is wrong. That needs z mod p < C0^2 + C0, probability (C0^2 + C0)/p < 2^-191.99 for a uniform residue. The reduced values (the batch products and inverses, λ, λ^2 and λ·u) are
products of pseudo-random field elements, the premise of the exceptional-addition row, so the event is a failure row of 7·(S_ops + 2^61)·(C0^2 + 3C0)/p < 2^-62.48 per run (3C0 covers the loops' fused
forms below; at most 7 reductions per batched operation, section 5.3). The experiment builds such inputs (for a high half H it solves for L with s mod 2^256 in [2^256 − s4·C0, 2^256), boundary
included) and checks the kernel returns exactly z mod p − C0.

**Fused reduction (24 + 5k).** For x3 and y3 the reduction also takes k canonical subtrahends (k = 2: x_W and x_T; k = 1: y_W). After step 2 it subtracts each one from all five limbs of s (SUBS, 3
SBCS, then SBC on s4: 5 each), so s' = s − Σ subtrahends with s4' ∈ [−2, C0]. It then folds s4' as a signed word: ASR for its sign, MUL and SMULH for the 128-bit product w = s4'·C0, and ADDS, 2 ADCS
and ADC of (w0, w1, sign, sign) into s' mod 2^256. That is **12 + 5 + 5k + 1 + 2 + 4 = 24 + 5k**. t = (s' mod 2^256) + s4'·C0 ≡ z − Σ subtrahends (mod p), and the result is t mod 2^256, which is right
unless t leaves [0, p):

-  t ≥ 2^256 (s4' ≥ 0). As above, this needs (z − Σ) mod p < C0^2 + C0 and gives (z − Σ) mod p − C0; it is the dropped-carry row, with the plain reductions.
-  t < 0. This needs s4' < 0 and s' mod 2^256 < |s4'|·C0. With one subtrahend it cannot happen: a borrow from s4 = 0 leaves s' mod 2^256 ≥ 2^256 − p > C0. With two it needs s4' = −1 and s' mod 2^256 <
  C0, probability C0/2^256 for x3, and gives a result ≡ z − Σ + C0.
-  p ≤ t < 2^256. The result is congruent but not canonical, probability C0/2^256 each for x3 and y3.

The experiment constructs each case (a = 2^128 gives a^2 = 2^256 and s = C0; b = p − 1 − e with c = C0 + 1 + e + d, 1 ≤ d ≤ C0, makes the fold go negative) and checks the kernels' predicted values. A
pass-2 operation makes one fused reduction for x3 and one for y3; the doubling does the same. The loops' y3' = λ·u + y_W adds y_W to all five limbs instead (ADDS, 3 ADCS, ADC: 23 + 5), so s4' ∈ [0, C0
+ 1], w1 ≤ 1 and the fold is the plain one; t < 0 cannot occur, and t ≥ 2^256 needs (z + y_W) mod p < C0^2 + 2C0 (the dropped-carry row). The loops' x3 reads its second subtrahend from a table field
x̃_T = x_T + 2C0 (section 5.1) and ends both borrow chains with ADC s4, s4, xzr, which adds 1 − borrow, instead of SBC (23 + 10). So s' = s − x_W − x̃_T + 2^257 = s − x_W − x_T + 2p ≡ z − x_W − x_T
lies in [2, (C0 + 3)·2^256): s4' ∈ [0, C0 + 2], w1 ≤ 1 since (C0 + 2)·C0 < 2^65, the fold is the plain one, t < 0 cannot occur, and t ≥ 2^256 needs (z − x_W − x_T) mod p < C0^2 + 3C0. In the loops,
the carry of a plain fold or of these folds into limbs 2 and 3 runs as a rare exit (section 2.5); the results are the same.

### 2.3 fe_sqr (64)

-  **Off-diagonal sum** O = Σ_(i<j) a_i·a_j·2^(64(i+j)) (12 + 8). Six products. Columns 1 to 6 hold 1, 2, 3, 3, 2, 1 words, so six two-word additions in two chains: chain A adds L03+L12 (column 3),
  H03+H12 (column 4) and its carry into H13; chain B adds H01+L02, then H02, L13, the incremented H13 with L23, and its carry into H23, which is o6. As in fe_mul, the carry-only additions land on
  untouched high words, and o6 is the top limb of O < 2^448.
-  **Doubling** (6). 2·O ≤ a^2 < 2^512, so no bit is lost: LSR d7 = o6 >> 63 and five EXTR #63 give d2 to d6. The low limb d1 = 2·o1 mod 2^64 is not formed: the first diagonal addition is ADDS z1,
  H00, o1, LSL #1, and the bit shifted out of o1 is bit 0 of d2, already taken by the EXTR.
-  **Diagonals** (8 + 7). a_i^2 by MUL and UMULH, then ADDS z1, five ADCS (z2 to z6) and ADC z7 = d7 + H33 + c, the top limb of a^2 < 2^512, which cannot overflow.
-  **Reduction** (23).

**64 = 20 + 8 + 6 + 7 + 23.**

### 2.4 fe_sub (9)

**fe_sub(a, b)**, a weak, b canonical: SUBS, 3 SBCS; CSEL of C0 on borrow; SUBS, 2 SBCS, SBC. On borrow the wrapped value is at least C0 + 1, so subtracting C0 cannot borrow again; the result is weak,
≡ a − b. Without a borrow the result is a − b, which is at least p only if a ≥ p and b < C0. The step uses it for dx = x_T − x_W (pass 1), dy = y_T − y_W and u = x3 − x_W; b is a canonical walk
coordinate and a is canonical (x3 outside the events above), so the results are canonical. fe_sub2 of 13c3295f is no longer used. In the loops its last three instructions (SBCS, SBCS, SBC) are a rare
exit (section 2.5).

### 2.5 The two loops (99.50390625 and 432.50390625 per walk)

`kernels/secp_loops.s` holds each walk's body once, as macros W1, W2P and W2 (P1 and P2: a pair's LDP and the invocations; the fe_mul, fe_sqr, row, fold and fe_sub schedules are macros with register parameters) that
expand to straight-line code; counts are of the expansion's main path (a rare exit not taken costs its branch).

**Pass 1 walk** (forward): two LDP each for x(T_i) and x(W_i) (4); dx_i = fe_sub(x(T_i), x(W_i)) (7, canonical); two STP each for acc_(i−1) and dx_i into slot i (4); acc_i = dx_i·acc_(i−1) reduced
into the four loop-carried registers (84): 99, plus the table pointer, the macro parameter. Before each pair of the 512-walk main loop one LDP x9, x26, [x2], #16 loads both table pointers (the first
walk reads x9, the second x26, which no body writes); SUBS and B.NE on a group counter in x29, unused by the body, close the 512 walks: 512·99 + 256 + 2 = **50,946 per 512 walks, 99.50390625 per walk**.

**Pass 2 walk** (backward, i = n−1 … 1). R ≡ 1/(dx_0···dx_i) stays in x25–x28, loaded once from Rinv, never stored. Per walk:

| work | instructions |
| --- | ---: |
| 2 LDP each of acc_(i−1), dx_i; R's rows B_1..B_3; inv = R·acc_(i−1), R ← R·dx_i on them (W2P) | 4 + 3·4 + 2·72 |
| 2 LDP y(T_i), 2 LDP y(W_i) (kept to the end); dy = fe_sub(y_T, y_W) | 4 + 7 |
| λ = dy·inv | 84 |
| x3 = λ^2 − x_W − x̃_T + 2^257, fused; 2 LDP x(W_i) (kept), x̃_T's 2 LDP inside it | 73 + 2 + 2 |
| u = fe_sub(x3, x_W); 2 STP x3 | 7 + 2 |
| y3' = λ·u + y_W, fused; 2 STP y3' | 89 + 2 |
| **total**: 2·72 + 3·4 + 84 + 73 + 89 + 2·7 + 12 loads + 4 stores | **432** |

The main loop runs 512 walks: before each pair, LDP x0, x5, [x2, #-16]! loads both table pointers (x2 one entry above the next walk; the first walk reads x5, kept free by its body, the second x0,
outside the pool). SUBS and B.NE on a group counter in x18 (no body uses it; set per batch to (n − 1) >> 9 by LSR) close them: 512·432 + 256 + 2 = **221,442 per 512 walks, 432.50390625 per walk**. B.NE
reaches ±1 MiB (section 2.6). Walk 0 uses R itself as its inverse after the loop. With λ = (y_T − y_W)/(x_T − x_W) and x3 = λ^2 − x_W − x_T, the stored point (x3, λ(x3 − x_W) + y_W) is −(W + T), canonical outside
the events of section 2.2.

**R on shared rows (W2P).** B_0 = R and B_(i+1) = [lo, b0 + hi, b1, b2] with (lo, hi) = b3·C0 (MUL, UMULH, ADDS, B.CS: 4), so B_(i+1) ≡ B_i·2^64 because 2^256 ≡ C0. For a = acc_(i−1)
and a = dx_i, Z = Σ a_i·B_i ≡ R·a: row 0 is 8 products and 4 additions into five limbs, rows 1–3 each 8 products, 4 additions and 6 more (ADDS, 4 ADCS, ADC) into Z, interleaved for the two
products; no carry is lost (a high word is at most 2^64 − 2) and Z < 2^322, so z5 ≤ 3. The fold (6) is MUL, UMULH of z4·C0, MADD hi' = z5·C0 + hi ≤ 4C0 − 1, ADDS z0 + lo, ADCS z1 + hi', B.CS:
t = (Z mod 2^256) + (z4 + 2^64·z5)·C0 ≡ Z. Each product is 32 + 34 + 6 = 72; with the rows and 4 LDP the pair is 160 instead of 2·86 = 172. Inputs R in x25–x28 and the slot at x3; outputs inv in
x10–x13 and R' in x25–x28; all 25 pool registers, never x0 or x18, so P2 loads its pointer pair after the first walk's W2P (W2P, LDP x0, x5, W2, W2P, W2) and the tail its LDR after W2P.

**The stored member** (4595bad7). The walk keeps any member of a class (`upstream-proof.md` section 2). For S with t = 3s + j, −S has the same x and K (so j, h and key fields) and y of the other
parity (y ≠ 0: prime order), so s' = 1 − s, Tab[h][t'] = [−u_t^-1]M_h = −Tab[h][t] and −S + Tab[h][t'] = −(S + Tab[h][t]), in the same class. Storing −(W + T) changes no class, key, index, fingerprint
or report, only the member held; doublings, traversals, re-walks and solving act on classes or recompute their members.

**Tails.** Each pass first runs its (n − 1) mod 512 left-over walks in a tail loop (LDR of the table pointer, the body, SUBS, B.NE: 102 and 435 per walk), then the groups; an empty group loop is
skipped by one CBZ. With n = 2^20 every round runs (n − 1) mod 512 = 511 tail walks of each pass, 2.49609375 instructions each above their group rates. Ordinary rounds run the fused pass 1 instead (tail
walk 133; pass1p 128); section 4 charges the round's tails, 3315.51171875, with its fixed work. A doubling changes the tail counts but not their maxima, so it carries no tail allowance. The setup (LSR and
AND for the counts, CBZ for each tail) is per-batch work.

**Rare exits** (kaankolcu's a9297fb0). y3's fused form adds y_W, x3's borrow chains end with ADC and its second subtrahend is read at [p, #192] (section 2.2); the tails below are
replaced by a branch to a stub `.Lr…` after the function's RET (macros S1 and S2, one copy per walk copy, nested macros listing the copies; B.LO
and B.CS reach ±1 MiB; the farthest, a pass-2 stub's B.CC back to its walk, spans 888,540 bytes), which ends with B to `.Lb…` just after the exit. An exit is taken exactly when its tail would change a register:

-  **fe_sub.** B.LO after the SUBS of the CSEL'd C0 from limb 0; the stub is SBCS, SBCS, SBC. It is taken iff a − b borrowed and the wrapped low limb is below C0. Otherwise C is set, and the tail
  would subtract 0.
-  **Plain reduction.** B.CS after the ADDS, ADCS of (w0, w1) into limbs 0 and 1; the stub is ADCS, ADC of 0 into limbs 2 and 3. It is taken iff s1 + w1 + c ≥ 2^64, so only if s1 ≥ 2^64 − 2 (w1 ≤ 1,
  section 2.2). For acc and R, step 2 writes limbs 2 and 3 straight into the output registers, whose operands are dead by then.
-  **y3's fused addition.** s4' ∈ [0, C0 + 1] (section 2.2), so its fold is MUL, UMULH, ADDS, ADCS and B.CS to a carry stub as above; there is no sign exit.
-  **Rows and their folds** (W2P). A row's B.CS exits when b0 + hi carries (b0 ≥ 2^64 − C0 + 1, at most C0/2^64), a fold's when z1 + hi' carries (z1 ≥ 2^64 − 4C0, at most
  4C0/2^64). The stub ST adds the carry into limbs 2 and 3 (ADCS, ADCS, B.CC back); if that also carries, the value wrapped 2^256 and is below 2^99, so it adds C0 into limb 0 (ADDS, ADCS, ADCS, ADC, B
  back). inv and R' are therefore always exact weak residues, and no failure row is added.
-  **x3's fused fold.** With the 2p bias s4' ∈ [0, C0 + 2], so its fold is MUL, UMULH, ADDS, ADCS and B.CS to a carry stub. (Without the bias s4 = 0 has probability about 2^-16 for a square and the
  fold can go negative.)

Each stub first runs the macro RCHG (STR x0 pre-indexed, LDR, SUB #20, STR, LDR x0 post-indexed), taking 20 milli-units from the caller's x19 (the round's reserve, section 4), saved at [sp] by the prologue and restored by the
epilogue. None of the five sets flags and x0 is restored, so the tail sees the flags and registers it would have seen inline. A stub is at most 13 instructions (ST: 5 + 3 + 4 + B). There are 6,162
stubs: each of the 513 pass-1 walk copies (512 in the main loop, 1 in the tail loop) has 2, and walk 0 has 1; each of the 513 pass-2 copies has 10 (two folds, three rows, dy, λ, x3, u, y3'), and walk 0 has 5. A walk enters at most 12 stubs (2 in
pass 1 and 10 in pass 2), 240 milli-units. The doubling's fe_mul, fe_sqr and x − x3 take the same forms (section 4); its fused x3 and y3 keep their inline folds.

**Registers.** The pool is x5–x17 and x19–x30 (25), x19–x30 saved by the prologue and restored before RET; the loops call nothing. x0 holds the tail count, then in pass 2's groups the second table
pointer; x1–x3 the three array pointers, x4 C0, x18 pass 2's group count, x26 pass 1's second table pointer, x29 pass 1's group count. The pool's peak is 25 in pass 2 and 16 in pass 1, with no spills,
and only the stubs touch the stack. Each loop is a single basic block apart from its rare exits: every other conditional is a carry chain or CSEL, and the only other branch is the back edge.

As in the public package, one LDP or STP is one instruction.

**Fused pass 1** (`kernels/secp_fuse.s`, `pass1f` in compare rounds, t ≡ 0 (mod 128); generated from W1's multiply and the arrival rows of section 3.2). Each walk runs the arrival of section 3.2 on its current point,
then dx and the product: LDR of the K_cp word, 3 LDP of x(W) and y0, y1, the F1 row and index (21), TBNZ, the key (2: CMP y1, #0 and ORR, or CMN y1, #2 and SUB from p0 + 1, section 3.2) and the chain (5), B.NE, BFI and STR of the address (32), then 2 LDP of x(T),
fe_sub (7), the 4 STP of acc_(i−1) and dx_i and fe_mul (84): **129**. 512 walks run straight; LDR, SUBS, STR and B.NE of a group count kept on the stack close them: 512·129 + 4 = **66,052 per 512
walks, 129.0078125 per walk**. A tail walk is 133; walk 0 forms dx_0 into the acc registers. Registers: x0 K_cp words, x1–x3 the arrays, x4 C0, x5, x6, x17, x19 acc, x7–x11 G_1..G_5, x12 p0 + 1, x14–x16
b_0..b_2, x18 the constant 33, x20–x23 x(W) then dx, x13 and x24–x30 the arrival's values and the multiply's temporaries: all 31, no spills. **`pass1p`**,
the other rounds, has no K_cp load, filter word or filter CCMP and ends its chain with B.HS: after CCMP d, C is set exactly when there is no distinguished candidate or window. Even y tests
K1 = 0 by CMP y1, #0; odd y by CMN y1, #2: y1 = ~K1 − b0 (b0 the borrow of p0 − y0), so K1 = 0 gives y1 ≥ 2^64 − 2 and sets C, which sends the arrival to the rare path; it rejects the false
candidates (K1 ≠ 0, below 2^-63 per arrival). The arrival is 28 (3 LDP, F1 11, select 6, TBNZ, one key test, 4 windows, B.HS, BFI); a pair's two addresses go out in one STP (the first in x12,
p0 + 1 being unused in pass1p), and x0, unused there, counts the tail and group loops (SUBS, B.NE): walks 125 and 126 per pair, **64,258 per 512 walks, 125.50390625 per walk**, tail 128, fixed 77,
below pass1f's, which bound the round's work (section 4). Its rare path has no filter part. The
two kernel exits (20 milli-units each) charge the caller's x19 at [sp] as in pass 1 and use only x13, which is dead there. **Rare path.** It is that of section 3.2: an action of 25 milli-units taken
first from the same word, then at most 22 instructions (the rare-path row of section 5.3), after the walk's tests and before its LDP of x(T). It may use x13, x25–x27, x29, x30 and the stack, keeps
x0–x12 and x14–x19, and leaves x(W) in x20–x23, b_j in x24 and the BFI source in x28 for the walk's current point (after a fallback, the exact rule's j; after a report, the next trail's start, arrived
inside the start's 2^13 units). A report, start or fallback is its own action, charged before it runs from M (via &M at [sp, #24], the caller's x22: LDR, LDR, SUBS, B.CC to the
refill, STR), never from the reserve at [sp] (section 4); a traversal's surcharge is charged from x5 after the passes, before its first step. The dx is that of the current point. The shipped `secp_fuse.s` has a test stub there instead: after the 25-milli-unit charge it appends (x13, x1, x28, x29) to a list at
[sp, #120] and branches to the BFI, so the address is the kernel's own j; it uses only x13 and the stack. Every walk makes exactly one arrival and one operation per round, and round 0's pass 1 arrives
the starts. When phase 1 stops, `arrive.s` runs once over the last round's new points before phase 2, so every phase-1 point is arrived. It is not a round, so no round reserve covers it: the pass is one action priced **2^15 units plus 25 milli-units per walk = 2^15.85
units**, charged before it runs and in the 2^70 "other" row (5.3), which covers its 2^20·31.015625 instructions with setup < 2^10 (< 2^15·1000) and at most one rare entry per walk at 25. Its rare entries log their kind in the rare list and draw nothing from the meter — `arrive.s` has no meter decrement — so the word it saves at [sp] is handed on unspent and the next action's
remainder is not reduced by the pass. `certificate.py` asserts both inequalities.

### 2.6 Evidence

The declared run's checks (items 1 to 8), the off-runner kernel runs and the interpreter are in `evidence.md`.

## 3. Arrival rows: 53 → 31.015625 plus the round share, executed

The arrival rows are the arrival pass's work per class arrival (section 5.1). Ordinary rounds run them in the fused pass 1 (section 2.5), anchor rounds in the loop below, `kernels/arrive.s` (run in
the declared run). Each arrival takes 30 instructions on either parity path; each pair adds one LDP of their K_cp words and one STP of their addresses (a9297fb0's pairing), and each 32 pairs one loop
branch. That is 31.015625 per arrival, plus the step's share of the round's fixed work (section 4).

| row | this | sequence |
| --- | ---: | --- |
| F1 truncated product | 13 | 2 LDP of x; 4 UMULH, 3 ADD, 4 MADD (G_1..G_5 in registers); F0 = top limb of x |
| canonical index, slot | **10** | ADD, SUBS, CSEL, CSEL, ADDS, CSEL; CCMN ×3, CCMP (section 3.1) |
| sign, key, DP test | **3** | TBNZ on bit 0 of y0; odd y: SUBS, SBCS; even y: CMP, ORR |
| h and address | **1** | BFI of bits 0..25 of y0 into the slot base (odd half permuted, 3.2) |
| traffic | **1.5** | LDP of y0, y1; half an STP of the pair's two addresses |
| filter, rare branch | **2.5** | half an LDP of the pair's K_cp limb-0 words; CCMP; B.NE past the rare path |
| arrival loop | **0.015625** | one TBZ per 32 pairs on the pointer array's end bit |
| **total** | **31.015625** | |

### 3.1 The canonical index

As in `upstream-proof.md` section 1.1, v0 = x, v1 = βx mod p and v2 = β²x mod p are distinct in [1, p) with sum p or 2p. Let A, B, C = 2^64·v0/p, 2^64·v1/p, 2^64·v2/p, reals in (0, 2^64) with A + B +
C ∈ {2^64, 2^65}. F0 is the top limb of x and F1 the truncated product. In 64-bit words the row forms

```text
S = F0 + F1 mod 2^64        d = F0 - F1 mod 2^64 (SUBS: carry iff F0 >= F1)        m = min(F0, F1) (CSEL) e = m + S (ADDS: carry iff m + S >= 2^64) j
= 2 if e carries, else 1 if F0 >= F1, else 0
```

With F2* = 2^64 − S ∈ [1, 2^64] (2^64 when S = 0), e carries exactly when m ≥ F2*. The row falls back to the exact rule (`upstream-proof.md` section 1.1) when any of four windows holds:

```text
(W1) F1 >= 2^64 - 16     (W2) S >= 2^64 - 17     (W3) e mod 2^64 >= 2^64 - 33     (W4) d < 16
```

**One-sided errors.**
-  F0 = A − a with 0 ≤ a < (r + C0)/2^192, where r = x mod 2^192. Indeed 2^256/p = 1 + C0/p, so A = (x/2^192)(1 + C0/p) = F0 + r/2^192 + x·C0/(2^192·p), and x < p. So a < 2, and a < 1 unless r > 2^192
  − C0.
-  F1 ≡ B − b (mod 2^64) with 0 ≤ b < 16. The true bits 320..383 of x·G are T ≡ F1 + c (mod 2^64) with 0 ≤ c ≤ 14 (they lie in [F1, F1 + 14], `upstream-proof.md` section 1.1), and T = ⌊B − ε⌋ mod 2^64
  with 0 ≤ ε < 2^-64 (x·G/2^384 = xβ/p − ε/2^64). So b = (B − ⌊B − ε⌋) + c < 2 + 14. If B − b < 0, then F1 = B − b + 2^64 ≥ 2^64 − 15, which (W1) excludes. Outside (W1), F1 = B − b as real numbers.
-  Then F0 + F1 = A + B − a − b. Since A + B + C ≡ 0 (mod 2^64), F2* ≡ C + a + b (mod 2^64) with 0 ≤ a + b < 18. If C + a + b > 2^64, then F2* = C + a + b − 2^64 ≤ 17, so S ≥ 2^64 − 17, which (W2)
  excludes. Outside (W1) and (W2), F2* = C + a + b.

So F0 ≤ A, F1 ≤ B and F2* ≥ C.

**Exactness.** Assume none of (W1) to (W4).
-  j = 2 (m ≥ F2*): A = F0 + a ≥ F2* + a = C + 2a + b ≥ C and B = F1 + b ≥ F2* + b = C + a + 2b ≥ C. The v_i are distinct, so C is the smallest. No gap test is needed.
-  j = 1 (F0 ≥ F1 and F1 < F2*): e = F1 + S < 2^64, so (W3) gives F2* − F1 = 2^64 − e ≥ 34, and (W4) gives F0 − F1 ≥ 16. Then A − B = (F0 − F1) + a − b > 16 − 16 = 0 and C − B = (F2* − F1) − a − 2b >
  34 − 2 − 32 = 0.
-  j = 0 (F0 < F1 and F0 < F2*): (W3) gives F2* − F0 ≥ 34, so C − A = (F2* − F0) − 2a − b > 34 − 4 − 16 > 0. And B − A = (F1 − F0) + b − a ≥ 1 − a, which is positive whenever a < 1.

So the rule returns argmin(v0, v1, v2) unless j = 0 and r > 2^192 − C0. Under the kernel rows' random-point premise, a phase-1 arrival has r > 2^192 − C0 with probability at most (C0 − 1)·2^64/(p − 1)
< 2^-159.99 (C0 − 1 low parts for each of the 2^64 top limbs). Section 5.3 charges the failure row (S_ops + 2^61)·C0·2^64/(p − 1) < 2^-33.28. Phase-2 and re-walk operations evaluate the exact rule
inside their 600-unit price. Part A checks `fast_j`, the bit-exact model of this rule, against the exact rule on 30,000 random, 10,000 near-tie (k = 150..249, including values near 0 and v0 + v1 near
p) and 1,000 low-r inputs, with 0 mismatches.

**Fallback rate.** (W1) covers 16 values of F1, (W2) 17 of S, (W3) 33 of e and (W4) 16 of d. If βx mod p is equidistributed (`upstream-proof.md` 1.1), each has probability at most its width over 2^64,
twice that for (W3) since e is 2F0 + F1 or F0 + 2F1: below 2^7/2^64 per arrival, under upstream's 2^-56, so the price 0.5 and the 2^70 "other" row of section 5.3 cover fallbacks.

### 3.2 The sequence

Register map of the arrival pass: x1 walk slots, x2 pointer array, x3 K_cp words, x5, x6 meter, x7–x11 G_1..G_5, x12 p0 + 1, x13 2^64 − 1, x4, x14, x15 slot bases b_j = base + 64j, x16 the constant
33; x17, x19, x22–x29 one arrival's values; in a pair x20 and x21 the first walk's K_cp word and address, x0 and x30 the second's; x18 the left-over pairs. Constants are set once per batch.

-  **Order.** F1 (13); ADD S, SUBS d, CSEL m, CSEL of b1 or b0 on the SUBS carry, ADDS e, CSEL of b2 or the previous choice on the ADDS carry; TBNZ on bit 0 of y0 splits the arrival into two paths
  that rejoin at the BFI.
-  **Key.** K is the even one of y and p − y (13c3295f). Even y: K0 = y0, K1 = y1; CMP y1, #0 sets Z for the distinguished test and ORR forms the filter word y0 + 1. Odd y: SUBS from the register
  holding p0 + 1 gives K0 + 1 into a free register (y0 stays for the BFI), because y0 is odd and p0 is odd, so y0 ≠ p0 + 1 and the borrow is that of p0 − y0; this is the filter word. Then SBCS from 2^64 − 1 gives
  K1, with Z set exactly when K1 = 0. Both paths set flags only after the CSELs have used the decision flags.
-  **Chain and branch.** CCMN F1, #16, #2, NE; CCMN S, #17, #2, LO; CCMN e, x16, #2, LO; CCMP d, #16, #0, LO; CCMP of the filter word with K_cp's limb-0 word, #4, HS; B.NE. A CCMN test passes with C
  clear, the CCMP d test with C set; a skipped test carries the failing C forward (#2, #0) and #4 sets Z. So B.NE jumps to the BFI exactly when there is no distinguished key, window or filter hit;
  otherwise the rare path after it runs.
-  **Address.** One BFI, shared by both paths and both fused copies, writes bits 0..25 of y0 into bits 9..34 of b_j. Even y: y0 = K0, so bit 9 gets s = 0 and bits 10..34 get h, bits
  1..25 of K0 (R = 2^25). Odd y: K = p − y gives y0 ≡ p0 − K0 (mod 2^26) whatever the limb-0 borrow, so y0 mod 2^26 = 1 + 2·((P′ − h) mod 2^25) with P′ = ((p0 mod 2^26) − 1)/2 = 0x1fffe17: bit 9
  gets s = 1 and bits 10..34 get ρ = (P′ − h) mod 2^25. The s = 1 half of row ρ therefore holds −Tab[(P′ − ρ) mod 2^25][j]; ρ ↦ P′ − ρ is an involution, so every walk reads the entry of its
  class's (h, s, j) and the class walk is that of section 5.1. b_j has zeros there, so the result is base + 1024ρ + 512s + 64j, rows of 1024 bytes;
  an entry's x + 2C0 sits at +192 (pass 2 reads [p, #192]).
-  **Rare path.** Its own action at 0.025 units, charged first: LDR, SUB, STR of 25 milli-units on the saved meter (3); K1 recomputed from y0, y1 and p0, with CBZ to the report (at most
  5); the four window tests again and B.LO to the fallback (5). For a filter hit it forms K2 and K3 (at most 5: LDP of y2, y3, then for odd y SUBS, SBCS, SBCS, SBC of p − y
  from y0, which is still held), loads K_cp's other two words (LDP) and runs CMP, CCMP, B.NE to the BFI (4): at most 22. In a pair it keeps x0 and x21, delivers its walk's address in that
  walk's register (x21 or x30) and finds its K_cp entry at a fixed offset from x3 (in the LDP immediate). Equal fingerprints are a true repeat and go on to the traversal; the report and the fallback
  are charged their prices (40.4, 0.5) first. The distinguished test precedes the filter, as in section 5.1.
-  **Filter word.** K_cp's limb-0 word holds K0 + 1, odd; clearing stores 0, even, so a cleared checkpoint never hits, and for a different class the word is pseudo-random and odd: false hits at most
  2^-63 per arrival.
-  **Loop end.** x0 holds n; TBZ on its bit 0 sends an even n to the pairs, and for odd n walk 0 runs alone (own LDR of the K_cp word, STR of its address). AND x18, x0, #62 gives twice the pair count
  mod 32; CBZ skips to the main loop when it is 0, else those pairs run one per iteration (SUBS x18, x18, #2, B.NE), then TBNZ x2, #23 ends it. The main loop is 32 straight-line pairs (LDP, two walks,
  STP), closed by one TBZ x2, #23. The pointer array of the n = 2^20-walk batch is placed once so that &tp[n] is its first address with bit 23 set; walk 0's STR advances x2 by 8 and each pair's STP by
  16, so the pairs end exactly there. With n = 2^20 there are 2^19 pairs, a multiple of 32, so a round runs only TBZ, AND and CBZ before the main loop.
-  **Test stub.** In `kernels/arrive.s` and `kernels/secp_fuse.s` the rare path is an inline stub that logs the rare case (in `secp_fuse.s` after the 25-milli-unit charge, then on to the BFI, section
  2.5); the plain path keeps its one branch. The production rare path is the sequence above.

Other rows (2d25756c's):
-  **F1 by MADD.** G = ⌊2^384·β/p⌋ (`upstream-proof.md` 1.1); F1 = Σ low halves of x_i·G_(5−i) + high halves of x_i·G_(4−i) mod 2^64: 4 UMULH, 3 ADD, 4 MADD (xd = xa + xn·xm). G_1..G_5 are loaded once
  per batch.
-  **Anchor rounds and the cap.** Every walk makes exactly one operation (a step, a doubling or a traversal step) per batch round. Rounds with t ≡ 0 (mod 128) run `pass1f`, the others `pass1p`. In
  rounds with t ≡ 0 (mod 4096) the processor runs the anchor variant of the arrival pass. For every walk it forms K2, K3 and the word K0 + 1 (7), writes K_cp ← fp(K) (STP + STR, the checkpoint-write
  action) and adds one to the trail's count of anchor periods. It abandons the trail when the count reaches 20·2^52 (LDR, ADD, STR, CMP, B.HS: 5). A trail's first anchor comes within 4,096 operations
  of its start, so it is abandoned before L_max = 20·2^64 operations. A walk starting a trail records φ = t mod 4096 (12 bits) in the trail record; its anchors fall where (length + φ) ≡ 0 (mod 4096)
  and its compares where (length + φ) ≡ 0 (mod 128). So a re-walk or phase-2 continuation anchors and compares where the live trail did. Anchor rounds run this arrival pass and the unfused pass 1,
  31.015625 + 99.50390625 = 130.51953125 instead of `pass1f`'s 129.0078125 per walk. Section 5.3 charges the variant's 12 instructions plus that 1.51171875, rounded up to 14 per walk per anchor (asserted), as a hard row, and the checkpoint writes (2
  per anchor) as another.
-  **Zero denominator.** See section 4.

## 4. Prices and per-batch overhead

-  **Ordinary step:** 125.50390625 + 432.50390625 + 28,112/2^20 = **558.0078125 + 28,112/2^20 = 558.034622**, priced as that exact rational; compare rounds add pass1f − pass1p = 3.50390625 per walk, a hard row (section 5.3).
-  **Traversal step:** step + 24 + 2 + 7 + 6 = **597.0078125 + 28,112/2^20**: 24 to track and resume the minimum-fp member, the charge on its own path (SUBS of 39 milli-units from x5 after the passes, B.CC to the refill,
  before the traversal's first step in the next round), 7 forming K2 and K3 for the
  fingerprint the minimum compares, and 6 for its own end test against K_cp in every round (LDR, LDP, CMP, CCMP, CCMP, B.EQ).
-  **Doubling:** the batched doubling shares the product chain and inversion. Its arithmetic is 4 fe_mul (336) + fe_sqr of x (63) + dx = 2y by one weak addition (9) + 3x^2 by two weak additions with a
  second carry fold (2·14) + the fused x3 = λ^2 − x − x (75) + x − x3 (7) + the fused y3 = λ·(x − x3) − y (91) = 609, with the loops' rare-exit kernels of section 2.5 (fe_mul 84, fe_sqr 63, x − x3 7;
  after fkiene's 0a1c5772): each exit is its own 0.02-unit action, charged before its tail runs. A doubling takes at most 5 carry exits and 1 fe_sub exit, fewer than an addition's 6 and 3. Its slot in
  the fused pass 1 still runs that round's addition (below), so the rare-exit row of section 5.3 charges one more operation's exits for each expected doubling, and the doubling's stubs test their
  charge (SUBS of 20 from x5, the word its surcharge charges, B.CC to the refill; with a tail of at most 3 and B, at most 6 ≤ 20): they never draw on the round's reserve (meter, below).
  -  **Traffic: 30**: at most the addition loops' 5 + 4 + 2 and 13 + 4 + 2 loads, stores,
    pointer updates and loop control. A doubling walk reads its y for dx = 2y forward and acc, dx, x and y
    backward, with no table point, stores the same dx, acc, x3 and y3, and shares R in registers.
  -  **List moves: 34.** Upstream's 16 loads and 16 stores of the walk's x and y to and from the doubling list,
    paired into 8 LDP and 8 STP, plus 18 other instructions.
  -  **Escape logic: 3.** The clearing of the history (one STR of 0 into K_cp's limb-0 word) and the charge of
    the surcharge on its own path (SUBS of 677 milli-units from x5, B.CC); upstream's 30 also paid for the
    distinguished-bit scan and trail teardown, which checkpoint-only detection does not run.

  The arrival that ends a traversal runs inside pass 1 at the walk's slot in a straight-line group, so that slot still runs its addition (result discarded) while the walk doubles on the doubling list,
  and the doubled point is the slot's member for the next round, whose fused step arrives it. The doubling is the displaced step plus 609 + 30 + 34 + 3 + 28,112/2^20, below step + 676.027, priced at
  **step + 677**, a surcharge of 677 milli-units. The weak additions, the 18 and the list moves are those of `upstream-proof.md` section 5.3.
-  **Zero-denominator test, once per batch.** Pass 1 multiplies the n canonical dx_i into acc_(n−1); p is prime, so outside the rare event of section 2.2 acc_(n−1) ≡ 0 (mod p) exactly when some x(T_i)
  = x(W_i). Before the inversion acc_(n−1) is canonicalised and tested (ADDS, 3 ADCS, 4 CSEL, 3 ORR, CBZ: 12 per batch). Only when it fires, with probability < 2^-114 over the whole run, does the
  batch take the cold path: the unfused pass 1 is re-run on the addresses the fused pass stored, with the per-walk test, each walk with dx_i = 0 is abandoned and its dx_i replaced by 1, as in the
  public package. The re-run, at most 2^20·(102 + 4) < 2^27 main-path instructions, is its own action at 2^17 units in the "other" random row; its stubs charge themselves as in section 2.5.
-  **Per-batch fixed work:** inversion I ≤ 24,322 (`secp_fe.s`, same .text, `upstream-proof.md` section 5.1–5.2), the fused pass 1's prologue, walk 0, setup and epilogue (80, executed: 25 + walk 0's 41
  + 14; `pass1p` 77), pass 2's (222,176 − 435 − 221,442 = 299), 64 for the call and return, 12 for the zero test, 3 for the context words p0 + 1 and 2^64 − 1, 7 for the meter, 3 for the round counter
  and the anchor-round test (ADD, TST, B.EQ), 2 for the compare-round test (TST, B.NE) and 4 for the arguments of the fused call: 24,796 < 2^15. All 2^20 walks form one batch (not upstream's 2^15 to
  2^16): every walk makes one operation per round. With n = 2^20 the tail walks add 511·(133 − 129.0078125) + 511·2.49609375 = 3315.51171875 every round (section 2.5), so a round's fixed work is at most
  28,111.51171875 instructions, charged as 28,112. The batch's scratch is 72 bytes per walk (slots and pointers), 72 MiB in all, inside the memory bound of `upstream-proof.md` section 5.3.
-  **Meter per round** (a9297fb0's reserve, kept in its own word). The round charge is c = 878.0078125·2^20 + 28,112 milli-units = (step price + 320)·2^20: 558.0078125 and a reserve of 320 per walk, plus the fixed work,
  plus 3.50390625·2^20 in compare rounds (`certificate.py` asserts each); the two values stay in the callee-saved x21 and x23, and each branch of the compare-round test charges its own. da765727's section 4 printed c = 878.03125·2^20 + 26,448; its 558.015625 + 320 is 878.015625. Here c = (558.0078125 + 320)·2^20 + 28,112, which `certificate.py` asserts. Before a round the processor charges c from its block
  remainder x5 (SUBS, B.CC to the refill), adds 2^20 to its operation count x6 (ADD), stores the remainder in the meter word M (STR x5, [x22]; the callee-saved x22 holds &M all phase) and sets the
  reserve x19 = 320·n (MOV); after the passes it reloads M and adds back the unspent reserve (LDR, ADD): 7 instructions, inside the 24,796. A plain round is charged the step price per walk plus what
  the reserve paid. The reserve x19 (saved at [sp] by each pass's prologue, restored by its epilogue) pays only the unchecked charges: the stubs' 20 milli-units each and a rare-path entry's 25, at
  most 12 exits and one entry per walk, or 2 more exits if the zero test re-runs the unfused pass 1, so at most 305·n < 320·n per round, and x19 stays non-negative without a test. Between the STR and
  the reload every other action (report, start, fallback, the zero test's re-run) charges M with a test (inside a pass via &M at [sp, #24], the caller's x22: LDR, LDR, SUBS, B.CC to the refill,
  STR), and after the reload the surcharges (677 or 39) charge x5; none charges the reserve, so no checked action can spend what a later unchecked charge needs. A round's charge is below 2^30
  milli-units, far below a 2^40-unit block, so one refill covers it and every walk is covered at the step price before it runs. Rare-path entries and rare exits are paid on their own paths, so every
  action is charged before it runs; the final arrival pass of section 2.5 is one action charged 2^15 units plus 25 milli-units per walk, in the 2^70 row. Blocks are reserved, refilled and counted
  as spent exactly as in `upstream-proof.md` section 3, and phase 1 still stops at a refill once the operation count reaches S_ops. The 2^20 processors' last blocks cover fewer than
  2^20·2^40/0.558 < 2^60.85 further operations, so phase 1 performs fewer than S_ops + 2^61 operations (sections 3 and 5.3).

Every other price (reports 40.4, probes 0.05, starts 2^13, fallbacks 0.5, phase-2 and re-walk operations 600, refills 0.03, words, attempts, store initialisation) is unchanged; the checkpoint write
(0.002), the rare-path entry (0.025) and the rare exit (0.02) are their own actions; the anchor-round work is a hard row (section 3).

## 5. The walk at R = 2^25, H1 and the exact budget

### 5.1 The walk change

The walk is that of `upstream-proof.md` section 2 with five changes: the table has R = 2^25 rows instead of 2048, cycles are detected by the checkpoint alone, the class key is the even one of y and p
− y, anchors fall where (length + φ) ≡ 0 (mod 4096), and the checkpoint is compared where (length + φ) ≡ 0 (mod 128). Everything else is unchanged: the canonical index, the twisted table Tab[h][t] =
[u_t^{-1}]M_h, the distinguished classes (limb 1 of K zero, θ = 2^-64), M = 2^20 walks and L_max = 20·2^64. At each class arrival, in order:

```
if K distinguished: report (K, a, b); end trail if length >= L_max: abandon trail if (length + phi) mod 128 == 0 and fp(K) == K_cp: traverse the cycle
once more (<= 32 steps, else abandon), resume at the member of its minimum-fp class remembered during the traversal, ESC; clear history; continue if
(length + phi) mod 4096 == 0: K_cp <- fp(K)   (phi: the trail's start round mod 4096, section 3) P <- P + Tab[h][t]
```

"Clear history" stores 0 in K_cp's limb-0 word (section 3.2). ESC is the doubling escape P ← [2]P, a class map.

-  **Merging lemma.** Let trails A and B arrive at the same class C, with any histories. Until one escapes, both apply F and visit the same classes. A trail escapes only after K = K_cp: K_cp is set at
  an arrival and the arrivals up to the next refresh are F-steps, so a repeat of K_cp is a cycle Z of F with no distinguished class, and the traversal escapes only if |Z| ≤ 32. Every escape from Z
  leaves from its minimum-fp class and goes to E(Z) = class(2c(min Z)), which depends only on Z. If A escapes from Z, B follows the same F path and enters Z, sets K_cp inside Z within 4,096 arrivals,
  compares there again lcm(|Z|, 128) ≤ 3,968 arrivals later, before its next anchor, and escapes to the same E(Z) after a traversal of |Z| and one doubling. So A and B meet again after at most 4,095 +
  3,968 + 32 + 1 = 8,096 extra operations of one of them, and by induction reach the same distinguished class unless one first hits L_max or an exceptional addition (failure events). History is
  cleared at every doubling and starts have none, so K = K_cp always witnesses a real cycle of F.
-  **Cycle rates.**
  -  A 2-cycle needs the next class to use the same partition index with the inverse twist: rate 1/(6R) =
    2^-27.585 per step.
  -  Other cycles of length 3 to 32 need every partition index used at least twice: rate at most
    Σ_L Σ_b S2(L, b)·R^-(L−b) = 2^-47.99999 at R = 2^25, over L = 3..32 and b = 1..⌊L/2⌋, where S2(L, b)
    counts the partitions of L steps into b blocks of size ≥ 2 (1, 1, 3, 15 for (L, b) = (3, 1), (4, 1), (4, 2),
    (6, 3)). To L = 64 it is 2^-19.9954 at R = 2048, the public value; `certificate.py` recomputes both exactly.
  -  Cycles longer than 32 end their trail (the traversal cap if compared, else the anchor-count cap) and are
    rarer still: below 2^-359 per step for L = 33 to 8,192 at R = 2^25 (exact to L = 256, above it
    S2(L, b) ≤ (L/2)^L), with the terms still falling. Such trails number below 2^-232 per run and cost below
    2^-100 units in expectation; `certificate.py` asserts both.
  -  A cycle with some partition index used once closes by chance: a fresh step lands on one of the trail's t
    earlier classes with probability at most t/N under H1's pseudo-random key fields, at most
    L_max^2/(2N) = 2^-117.8 per trail (H1's non-fruitless-return clause), under 2^-54 per run; their doublings
    and traversals cost under 2^-47 units, in the 2^70 "other" row (5.3). On the toy such cycles occur at about
    (trail age)/N per step.
  -  So the expected number of doublings is below (1/(6R) + 2^-47.99999)·(S_ops + 2^60) + 2^-54
    < 2^-27.5849·(S_ops + 2^60). A 2-cycle is traversed in 2 steps, other short cycles in at most 32, so the
    expected number of traversal steps is below (2/(6R) + 32·2^-47.99999)·(S_ops + 2^60) + 32·2^-54
    < 2^-26.5849·(S_ops + 2^60). Only linearity of expectation over the pseudo-random partition and twist of
    each new class is used.
-  **Detection delay.** A cycle Z is caught within 4,095 + lcm(|Z|, 128) ordinary steps of entering it (4,223 for
  |Z| = 2), ordinary work counted in S_ops, with its arrivals in T_c. With the traversal, the waste per operation
  is w ≤ 4,225/(6R) + 8,094·2^-47.99999 < 2^-15.54 at R = 2^25 (8,094 = 4,095 + 3,968 + 31, the largest 4,095 + lcm(L, 128) + L). `certificate.py` asserts the coverage and w ≤ 2^-14.56.
-  **Table.** 2^25 rows of 1024 bytes (section 3; the odd half of each row permuted as in 3.2), 2^35 bytes, aligned to 2^35. Each entry also holds x̃ = x + 2C0, so a row is redrawn (metered, as for M_h = O) when one of its three x-values is ≥ p
  − C0 (below 3C0/p < 2^-222 per row); the redraw depends on M_h only, so section 5.2's premises hold. Even with a copy on each of the 2^16 nodes, each node stays below 2^56 bytes and peak memory
  below 2^72 bytes (`upstream-proof.md` 5.3, "Memory"). Building it takes 2^25 draws of [c]P + [d]Q at ≤ 2^13 units each plus twists and x̃ fields: below 2^38.1 units, in the hard extras (< 2^39) and
  in preprocessing, below its 2^59 cap.
-  **Key and partition.** The class key is the even one of y and p − y; h is bits 1..25 of K (section 3), disjoint from the distinguished bits 64..127, routing bits 128..143 and bucket bits 144..193.
  The anchor phase φ only shifts where a walk's anchors and compares fall; each walk is anchored once per 4,096 of its own operations and compared once per 128, all that the merging lemma and the
  detection-delay bound use.

### 5.2 H1 at R = 2^25

H1 keeps its form: P(T_c ≤ S_ops) ≥ 1 − exp(−(S_ops/κ)^2/(2N')), with N' = N·R/(R − 1) and now κ = 1.0001, on the hazard census of 12ccd81d, which this package ships (`census_hz.py`,
`census-hz-q29.json`). The walk is unchanged by the rare exits, which change no result, and by the stored member (section 2.5). T_c counts every class arrival, including repeated arrivals inside
fruitless cycles, during the detection delay and in traversals; its other parts are as in `upstream-proof.md` section 4 (pseudo-random key fields, the non-fruitless return bound, a uniform twist index
independent of the partition index).

**Why T_c runs slow on a toy.** T_c ends at the first arrival at a class visited by *another* trail: with D classes visited, `own` of them by the current trail, the
collision chance is (D − own)/N, so E[T_c] ≈ √(πN/2) + ā with ā = (1 − θ)/θ the mean trail age: +0.45 % on q29 with dp_bits 5, but ā/√N = 2^-62.7 in production (at most S_ops·L_max/N = 2^-58.4
of the hazard).

**Hazard census.** `census_hz.py` runs the walk of section 5.1 on the toy groups (2d25756c's walk code) under the rule (1, 256). F and E(Z) are the same under the production rule, which changes only
when a cycle is caught (`evidence.md` item 8), hence only the repeated arrivals, the waste (1 + w) covers.
-  At every fresh-sourced arrival (a trail start, an escape, or a step from a class first visited at the previous arrival) it adds the random-mapping hazard (D − own)/N; a merge is an arrival at
  another trail's class.
-  For a uniform random mapping the hazards are the merges' compensator, so E[merges] = E[Σ hazards], and ρ = merges/Σ hazards measures the walk's collision intensity against the law with the
  own-trail exclusion removed exactly. Repeated arrivals add no hazard: they are the waste w.
-  In multi mode a run of 64·√N arrivals continues after each merge, which ends only the merging trail; in tc mode each trial stops at its first merge (T_c). `--pure` runs uniform random functions
  through the same trail logic: the control.

The table gives ρ and κ_eff = ρ^-1/2 from `census-hz-q29.json` (q29: N = 2^24.83; q25: N = 2^22.42). The error is the larger of the Poisson and table-spread standard errors.

| run | tables | merges | ρ | κ_eff |
| --- | ---: | ---: | ---: | ---: |
| walk, q29, R = 2^20, dp_bits 5 | 1024 | 52,224,351 | 1.00044 ± 0.00020 | 0.99978 ± 0.00010 |
| walk, q29, R = 2^20, dp_bits 3 | 256 | 13,050,652 | 0.99969 ± 0.00029 | 1.00015 ± 0.00014 |
| walk, q29, R = 2^20, dp_bits 7 | 256 | 13,026,503 | 0.99978 ± 0.00061 | 1.00011 ± 0.00030 |
| walk, q29, R = 2^16, dp_bits 5 | 256 | 13,044,795 | 1.00018 ± 0.00040 | 0.99991 ± 0.00020 |
| walk, q29, R = 2048, dp_bits 5 | 16 | 800,858 | 1.00100 ± 0.00199 | 0.99950 ± 0.00099 |
| walk, q25, R = 2^16, dp_bits 5 | 16 | 811,399 | 1.00202 ± 0.00186 | 0.99899 ± 0.00093 |
| walk, q25, R = 2^20, dp_bits 5 | 64 | 3,245,714 | 1.00035 ± 0.00087 | 0.99982 ± 0.00043 |
| walk, q29, R = 2^23, dp_bits 5 | 64 | 3,267,024 | 1.00142 ± 0.00078 | 0.99929 ± 0.00039 |
| walk, q29, R = 2^24, dp_bits 3 | 64 | 3,265,777 | 1.00063 ± 0.00067 | 0.99968 ± 0.00033 |
| walk, q29, R = 2^25, dp_bits 3 | 64 | 3,261,648 | 0.99936 ± 0.00056 | 1.00032 ± 0.00028 |
| ideal mapping, q29 N, dp_bits 3 | 16 | 815,722 | 0.99972 ± 0.00111 | 1.00014 ± 0.00055 |
| ideal mapping, q29 N, dp_bits 5 | 512 | 26,095,238 | 0.99976 ± 0.00028 | 1.00012 ± 0.00014 |
| ideal mapping, q29 N, dp_bits 7 | 16 | 818,029 | 1.00389 ± 0.00267 | 0.99806 ± 0.00133 |

-  Pooled over the three q29 runs at R = 2^20 (dp_bits 3, 5, 7): ρ = 1.00017 ± 0.00016, κ_eff = 0.99992 ± 0.00008 (χ² = 5.0 on 2 degrees of freedom).
-  Over the ten multi-mode runs of the walk: ρ = 1.00020 ± 0.00014, κ_eff = 0.999901 ± 0.000068 (χ² = 11.3 on 9).
-  At the production table size, R = 2^25, and at 2^23 and 2^24 (64 tables each): κ_eff = 0.99988 ± 0.00019 pooled (χ² = 5.1 on 2). The toy's distinguished bits sit above the R bits and
  q29 keys have 29 bits, so the runs at R = 2^24 and 2^25 use dp_bits 3.
-  The four q29 walk runs of 5f1ae3ce with 64 and 512 tables extend here to 256 and 1,024 on the same seeds (a chunk's seed depends only on R, dp_bits and its index), so each run contains the earlier
  one; the runs at R = 2^23 and 2^24 are 4c8222cc's and the run at R = 2^25 is 445f9a1c's. The rule below was fixed before each ran.
-  The ideal mapping: ρ = 0.99980 ± 0.00027 (χ² = 2.4 on 2), so the estimator is unbiased at this precision.

**The control at T_c.** In tc mode (32,000 trials each), the raw T_c slowdown against √(πN'/2) and the hazard at the first merge are:

| dp_bits | own-trail shift ā/√(πN'/2) | walk, R = 2^20 | ideal mapping | walk ρ | P[H ≤ ln(1/0.61)] (law 0.39) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | +0.10 % | +0.07 % ± 0.29 % | +0.61 % ± 0.30 % | 0.9974 ± 0.0068 | 0.3936 (walk), 0.3912 (ideal) |
| 5 | +0.45 % | +0.67 % ± 0.29 % | +0.11 % ± 0.29 % | 0.9973 ± 0.0056 | 0.3884 (walk), 0.3934 (ideal) |
| 7 | +1.85 % | +1.57 % ± 0.29 % | +1.78 % ± 0.29 % | 1.0104 ± 0.0075 | 0.3937 (walk), 0.3910 (ideal) |

Both raw slowdowns follow the predicted own-trail shift within 1.8 standard errors and grow with the trail length, while the walk's ρ stays within 1.5 standard errors of 1.

**First-collision distribution.** H1 needs P(T_c ≤ S_ops), not only the mean intensity. For a random mapping the hazard H accumulated up to the first merge is Exp(1), so P(T_c ≤ S) = P(H ≤ hazard
accumulated by S). The tc runs record H at each trial's first merge (`census-hz-q29.json`). For the walk the mean of H is 1.0026, 1.0027 and 0.9897 at dp_bits 3, 5 and 7. The six fractions with H ≤
ln(1/(1 − 0.39)) in the table lie within 1.4 binomial standard errors (0.0027) of the law's 0.39, five of them above it. So the walk's first collision follows the random-mapping law at the quantile H1
uses, as well as its mean.

**κ = 1.0001.**
-  In production a run of S arrivals has at least S/(1 + w) fresh ones, with w < 2^-15.54 at R = 2^25 (section 5.1: detection delay, cycle repeats and traversals). The bounds below use w = 2^-14.56,
  which covers it; `certificate.py` asserts both at that w from `census-hz-q29.json`.
-  Its own-trail term is below 2^-58.4 of the hazard, and its collision intensity is ρ.
-  So P(T_c ≤ S) ≈ 1 − exp(−ρ·(S/(1 + w))^2/(2N)), which is H1 with κ_eff = (1 + w)/√ρ.
-  κ rests on the multi-mode runs, against the law's ρ = 1, which the compensator makes exact for a uniform mapping. At R = 2^20 and dp_bits 5 the 1,024-table run gives κ_eff = 0.99978 ± 0.00010, so
  (1 + w)(κ_eff + 2 s.e.) = 1.00002. Pooled over the ten walk runs, κ_eff = 0.999901 ± 0.000068 and the same bound is 1.000078. κ = 1.0001 is the smallest four-decimal value at or above both
  one-sided bounds. Those two are the rule's operands; the three-run pool of section 5.2's table (q29 at R = 2^20, dp_bits 3, 5, 7) is printed for table scatter and is no operand, because its χ² of
  5.0 on 2 degrees of freedom (χ²/df = 2.5) shows plain inverse-variance pooling understates its error there — widened by √(χ²/df) it gives 1.0002, above κ.
-  The ideal control at 512 tables has ρ = 0.99976 ± 0.00028, 0.86 standard errors from its exact value 1, so the rule stays against the law.
-  No walk run's (1 + w)·κ_eff exceeds κ by one standard error. The closest, q29 at R = 2^25, is 0.00026 above κ, 0.93 standard errors; its (1 + w)(κ_eff + 1 s.e.) = 1.00064 is above κ. At R = 2^20,
  dp_bits 3 it is 0.00009 (0.65).
-  Against the ideal runs instead of the law (a cross-check, not the rule), the bound is 1.00005 for the main run over the 512-table ideal run, and 1.00013 pooled, above κ, because the ideal runs' own
  standard error (0.00014) enters.
-  The tc runs (25 to 1,632 times fewer merges) measure the first-collision law, not κ: walk κ_eff 1.0013 ± 0.0034, 1.0013 ± 0.0028, 0.9949 ± 0.0037 (dp_bits 3, 5, 7),
  within 1.5 s.e. of κ; ideal 1.0067 ± 0.0028, 0.9974 ± 0.0028, 1.0000 ± 0.0038, the estimator's spread.
-  N' = N·R/(R − 1) is kept. The twisted walk's 1/(6R) anti-collision is inside ρ, so N' over-covers it.
-  This is the standard random-mapping heuristic, supported at toy scale; no theorem bounds the constant. Brent and Pollard put an r-adding walk's excess at about 1/(2r), at most 2^-26 here (r ≥ R =
  2^25).

**Other checks**:
-  **Waste.** The census counts repeated arrivals: w = 1.46·10^-4 per arrival at q29, R = 2^20, dp_bits 5, mostly the detection delay after toy self-meetings; section 5.1 bounds production waste by
  2^-15.54.
-  **Trap test** (`long_trails`, q29, 150 trails of 8,192 operations, seed 7, this rule): non-fruitless returns in 0/150 trails at R = 4 (164 escapes) and R = 8 (161). Nearly every other arrival
  repeats a caught cycle (1,222,199 and 1,215,241 of 1,228,800), so the test has little power here; `evidence.md` item 8 checks classes.
-  **The hosted toy pipeline** draws a fresh table per trial at R = 2048 and adds the toy's own-trail term ā to its scaled budget: evidence of the whole procedure, not of κ, and below the production
  law. Its 2-cycle rate 1/(6R) = 2^-13.6 per step is 2^13 times production's, and a 2-cycle holds a walk up to 4,223 operations, past the toy trail cap 640, so the trail is abandoned: **640·2^-13.6 =
  0.052** is the abandonment the toy rate predicts. The measured share of visits is 0.0463, 0.0469, 0.0463 and 0.0466 over four runs of the declared input in the pinned image (286, 282, 286 and 290
  trials; 43, 43, 43 and 44 abandoned trails), 8.9 % to 10.1 % below that rate. A trial is lost only if its first merge follows an abandonment: 3,000 trials at seeds 1 and 2 (local runs of the declared
  input with more trials) succeeded at 0.3763 and 0.3640 under this rule and 0.3770 and 0.3953 under (1, 256), a drop of 0.0007 and 0.0313. In production the same product is
  **4,223·2^-27.585 = 2^-15.5409**, the 2-cycle part of w's delay term, w = 4,225/(6R) + 8,094·2^-47.99999 < 2^-15.54 (section 5.1): L_max = 20·2^64 is far above 4,223, so production buys a
  delay inside the meter and abandons nothing. The same delay is why in the rule check (`evidence.md` item 8) 599 of 600 starts end distinguished under both rules and the 600th does not: the caught
  cycles spend its fixed 2^15-operation budget (the 600 starts take 86,910 operations under (1, 256) and 1,138,828 under this rule), which at production's L_max never binds.

### 5.3 The exact budget

H1 and κ = 1.0001 remain the disclosed heuristic, supported only at toy scale (section 5.2). The Markov headroom is h = 0.000065 and the budget for every other failure event 1.1·10^-7, so p0 =
0.3900653 keeps the success floor at 0.39000019 > 0.39 with W below 2^125.858216. At R = 2^25 the random extras are 2^98.83.

```text
N'     = ((r-1)/6) * 2^25/(2^25 - 1)
S0     = sqrt(2*N'*ln(1/(1-0.3900653)))      = 2^126.6993990
S_ops  = 1.0001*S0                           = 2^126.6995433
phase1 = (558.0078125 + 28112/2^20)/1000*(S_int + 2^60) = 2^125.857970
hard extras                                  = 2^111.794273
compare rounds 0.00350390625*ceil((S_int + 2^61)/128)  2^111.542723
anchor rounds 0.014*ceil((S_int + 2^61)/4096)          2^108.541114
checkpoint writes 0.002*ceil((S_int + 2^61)/4096)      2^105.733759
upstream rows (table < 2^39, refills at the new cap)    2^107.155501
expected random extras = 2^98.8271
doublings  0.677*(1/(6R) + 2^-47.99999)*(S_ops + 2^60)              2^98.55
traversals 0.039*(2/(6R) + 32*2^-47.99999)*(S_ops + 2^60)  2^95.43
rare-path entries 0.025*(2*(1/(6R) + 2^-47.99999) + 2^-56 + C0/2^64)*(S_ops + 2^60) 2^94.83
rare exits 0.02*(14*C0/2^64 + 4*2^-63)*(S_ops + 2^61 + (1/(6R) + 2^-47.99999)*(S_ops + 2^60)) 2^92.86
starts 2^13*(S_ops/2^64 + 2^20)                                     2^75.70
words, attempts, probes, small batches, fallbacks, zero-test re-runs, final arrival   2^70
allowance A = ceil(extras/0.000065)          = 2^112.7362
W_ledger = phase1 + hard + A                 = 2^125.858215915
```

The rare-path row counts filter hits (true ones at most twice the doubling rate, the detecting arrival and the one closing the traversal; false ones below 2^-63), window fallbacks (below 2^-57),
distinguished keys (2^-64) and, in the fused pass 1, a C0/2^64 term per arrival kept as an upper bound on its distinguished candidates (at most 3·2^-64 in pass1p; odd borrows no longer enter). The rare-exit row counts, per operation:
-  the three fe_sub exits, each taken only when the low limb of a − b wraps below C0: at most C0/2^64 each;
-  the four carry exits (acc's and λ's reductions, x3's and y3'), each taken only when s1 ≥ 2^64 − 2: at most 2^-63 each;
-  pass 2's three row exits (C0/2^64 each) and two fold exits (4C0/2^64 each), section 2.5.

So the expected number of exits is at most 14·C0/2^64 + 4·2^-63 < 2^-28.192 per operation, at 0.02 units each, over fewer than S_ops + 2^61 operations, plus as many again per expected doubling, whose
slot also runs its addition (section 4). The rates rest on the random-point premise of the kernel rows. Every exit is charged before its tail runs (section 4) whatever the inputs, so the premise
affects only this expectation. `certificate.py` uses the exact rates and asserts the extras 2^98.827041655 < 2^98.8271; `experiment.py`'s float ledger, with the cycle rates rounded up, gives W_ledger =
2^125.858215922, extras 2^98.827101 (its flag: < 2^98.8272) and Markov ratio 0.0000649686 (rounded up).

In exact milli-units (`certificate.py`):

```text
S_int             = 138153841538979257067881737330846152944
phase 1           = 77094626767630612577307741602242369335866
compare rounds    = 3781860225233898730223335558014943930
anchor rounds     = 472205513072683007564117706240438936
checkpoint writes = 67457930438954715366302529462919848
hard (incl. all)  = 4502249465104058506063514598066794132
extras            = 562216931880721445024281950294342
A                 = 8649491259703406846527414619912953847
W_ledger          = 77107778508355420042660332531460349083845
cap (CAP_MILLI_LO)= 77107783047972467027495125526916811470112
```

W_ledger = 2^125.858215915 ≤ 1000·2^125.858216, the cap, which `certificate.py` bounds below from an interval for 2^0.858216 (exp(0.858216·ln 2) by its Taylor series with remainder); the next
lower six-decimal value, 125.858215, is below W_ledger and is not claimed. Against CAP_MILLI_LO headroom (CAP_MILLI_LO − phase 1 − hard = 2^112.736999), the Markov
ratio is at most 0.0000649660 (exact, rounded up) < 0.000065, so the meter halts early with probability below h. The 1.1·10^-7 covers every other failure event: the
rows of `upstream-proof.md` section 4 (trail cap < 10^-7, preprocessing cap < 2^-30 and smaller), the three kernel rows of section 2 and the canonical-index row of section 3.1: the dropped carry
7·(S_ops + 2^61)·(C0^2 + 3C0)/p < 2^-62.4, at most 7 reductions per batched operation on average, plain or fused (a doubling walk's slot adds its addition's 3 shared-chain reductions, whose λ, x3 and
y3 are discarded, to the doubling's 7: 6 + 4·2^-27.58 < 7 per operation); a doubling's negative fused fold (S_ops + 2^61)·C0/2^256 < 2^-97.2; a fused result in [p, 2^256) for x3 or y3, 2·(S_ops +
2^61)·C0/p < 2^-96.2; and a wrong j = 0, (S_ops + 2^61)·C0·2^64/(p − 1) < 2^-33.28. The rare exits return the inline results, so they add no failure row. Together these sum to 1.0103·10^-7. So

```text
success >= 0.3900653 - 0.000065 - 1.1e-7 = 0.39000019 > 0.39.
```

`certificate.py` checks in exact integer and rational arithmetic the kernel identities (85, 64, 91, 75; 84, 7, 89, 73 with the rare exits), the loops' 50,946 and 221,442 (pass 2's rows 4 and products 72), the fused pass 1's 66,052 and
64,258, the tails 3315.51171875, the rows, the per-round work 28,112, the round charge and the prices 558.0078125 plus 28,112/2^20, step + 677 and step + 39, the cycle sums and the long-cycle rate, the compare coverage and w, the two κ bounds of
5.2, the rare-exit rate, every ledger row, the four failure rows (with the mean 6 + 4·rate < 7 reductions per operation) and their sum with the upstream rows, the stubs and the 320 reserve, and the claim 125.858216 against 125.858215. It prints `CERTIFICATE: PASS`. The counter widths, preprocessing cap and memory are re-checked there.

## 6. Provenance, files and limitations

Design provenance: the construction, report store, ledger structure and F1 truncated product are the public line of `upstream-proof.md`, carried by 02feb284. The counted kernels,
loops, rows, checkpoint-only walk, census tool and κ rule come from bradcush's line (f99114e0, 736c001a, 13c3295f, 2d25756c, 67772ed0, 12ccd81d, 5f1ae3ce, 673640ef, 4c8222cc), with kaankolcu's rare
exits, reserve and paired traffic (807e4bf5, a9297fb0), −(W + T) forms (4595bad7), fused arrival (cbea87bf) and strided compares (4ac17b66); aburan28's headroom retune (b78fb4f6, 4b4defe7),
arrival pairs (2bd168b5, posted 2026-10-26 16:56Z, and 01488655) and R = 2^22 (610c8e49); fkiene's doubling recount (0a1c5772) and R = 2^21 (6a6cb91c); Ulysseus's census reruns (daf2c741); b5de673a's fused x3 subtractions;
R = 2^25, h = 0.000065 and the reserve split in bradcush's 445f9a1c; their combination in c06f6691 (125.8928); pass 2's shared rows in a00774d8 (125.862153); the permuted odd half of the table, one BFI of y0 on both parities and the shorter pass1p in
bradcush's 886e6f54 (125.858271); and walk groups of 256 in blocksek's da765727 (125.858233), the package this one edits. The submission declares 20 dependencies,
in this order: da765727, 886e6f54, a00774d8, c06f6691, 4ac17b66, 445f9a1c, cbea87bf, 4c8222cc, 673640ef, 0a1c5772, 610c8e49, 6a6cb91c, 5f1ae3ce, 4595bad7, 01488655, daf2c741,
a9297fb0, 807e4bf5, b5de673a, 02feb284; the others are credited here. New here: walk groups of 512. No prior acceptance, review outcome or score is a premise.

`upstream-proof.md` is the construction, abridged; its numbers are superseded here. `README.md` lists the files.

The limitations are the upstream ones. H1 is a random-mapping heuristic supported at toy scale; κ = 1.0001 rests on the hazard census of toy groups (N up to 2^24.83). The step price is the executed
main-path count of the plain fused pass 1 and pass 2 plus the round's fixed work over its 2^20 walks, compare rounds' 3.50390625 per walk a hard row; rare exits are their own actions, at rates that rest on the
random-point premise. The kernel and canonical-index failure rows rest on the same premise as the exceptional-addition row. The attack is generic and impractical, and the bound is not a claim about
deployed secp256k1.

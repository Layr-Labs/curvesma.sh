# The field kernels (proof.md sections 2.1 to 2.3, notation of proof.md section 2)

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

- **No carry is lost.** A high word is at most ⌊(2^64 − 1)^2 / 2^64⌋ = 2^64 − 2, and chains 1 to 5 end by adding one carry bit to a high word no earlier instruction has touched (H13, H22, H23, H32, H33), so those ADCs cannot carry out. Chain 6 ends by forming z7 after z0 to z6 are final; the eight limbs then satisfy Σ z_k·2^(64k) = a·b < 2^512, so z7 < 2^64. Every other addition's carry is consumed by the next instruction of its chain. Hence (z0, …, z7) = a·b exactly.

### 2.2 The reduction (23)

z = L + H·2^256 with L, H < 2^256 (the eight product limbs).

1.  **H·C0** (8 + 4): four MUL and four UMULH give g_j, k_j with g_j + 2^64·k_j = z_(4+j)·C0, k_j ≤ C0 − 1. ADDS g1+k0, ADCS g2+k1, ADCS g3+k2, ADC k3+c gives H·C0 in five limbs; H·C0 < 2^256·C0, so the ADC cannot overflow.
2.  **s = L + H·C0** (5): ADDS, 3 ADCS, then ADC into the fifth limb s4. Since s ≤ (2^256 − 1)(1 + C0) < (1 + C0)·2^256, s4 ≤ C0.
3.  **w = s4·C0** (2): MUL and UMULH; w ≤ C0^2 < 2^65.
4.  **t = (s mod 2^256) + w** (4): ADDS, 2 ADCS and ADC, which return t mod 2^256. t ≡ z (mod p) because 2^256 ≡ C0, and t < 2^256 + C0^2.

**23 = 8 + 4 + 5 + 2 + 4.** There is no final fold. If t < 2^256 the result is t, a weak value ≡ z (mod p). If t ≥ 2^256, then C0 ≤ t − p < C0^2 + C0 and t − p ≡ z, so z mod p = t − p and the result t
− 2^256 = z mod p − C0 is wrong. That needs z mod p < C0^2 + C0, probability (C0^2 + C0)/p < 2^-191.99 for a uniform residue. The values reduced this way (λ^2 and λ·u in the loops, the doubling's products) are
products of pseudo-random field elements, the premise of the exceptional-addition row, so the event is a failure row of 7·(S_ops + 2^61)·(C0^2 + 3C0)/p < 2^-62.48 per run (3C0 covers the loops' fused
forms below; at most 7 reductions per batched operation, section 5.3). The experiment builds such inputs (for a high half H it solves for L with s mod 2^256 in [2^256 − s4·C0, 2^256), boundary
included) and checks that the port returns exactly z mod p − C0.

**Fused reduction (24 + 5k).** For x3 and y3 the reduction also takes k canonical subtrahends (k = 2: x_W and x_T; k = 1: y_W). After step 2 it subtracts each one from all five limbs of s (SUBS, 3
SBCS, then SBC on s4: 5 each), so s' = s − Σ subtrahends with s4' ∈ [−2, C0]. It then folds s4' as a signed word: ASR for its sign, MUL and SMULH for the 128-bit product w = s4'·C0, and ADDS, 2 ADCS
and ADC of (w0, w1, sign, sign) into s' mod 2^256. That is **12 + 5 + 5k + 1 + 2 + 4 = 24 + 5k**. t = (s' mod 2^256) + s4'·C0 ≡ z − Σ subtrahends (mod p), and the result is t mod 2^256, which is right
unless t leaves [0, p):

- t ≥ 2^256 (s4' ≥ 0). As above, this needs (z − Σ) mod p < C0^2 + C0 and gives (z − Σ) mod p − C0; it is the dropped-carry row, with the plain reductions.
- t < 0. This needs s4' < 0 and s' mod 2^256 < |s4'|·C0. With one subtrahend it cannot happen: a borrow from s4 = 0 leaves s' mod 2^256 ≥ 2^256 − p > C0. With two it needs s4' = −1 and s' mod 2^256 < C0, probability C0/2^256 for x3, and gives a result ≡ z − Σ + C0.
- p ≤ t < 2^256. The result is congruent but not canonical, probability C0/2^256 each for x3 and y3.

A pass-2 operation makes one fused reduction for x3 and one for y3; the doubling does the same. The loops' y3' = λ·u + y_W adds y_W to all five limbs instead (ADDS, 3 ADCS, ADC: 23 + 5), so s4' ∈ [0, C0
+ 1], w1 ≤ 1 and the fold is the plain one; t < 0 cannot occur, and t ≥ 2^256 needs (z + y_W) mod p < C0^2 + 2C0 (the dropped-carry row). The loops' x3 reads its second subtrahend from a table field
x̃_T = x_T + 2C0 (section 5.1) and ends both borrow chains with ADC s4, s4, xzr, which adds 1 − borrow, instead of SBC (23 + 10). So s' = s − x_W − x̃_T + 2^257 = s − x_W − x_T + 2p ≡ z − x_W − x_T
lies in [2, (C0 + 3)·2^256): s4' ∈ [0, C0 + 2], w1 ≤ 1 since (C0 + 2)·C0 < 2^65, the fold is the plain one, t < 0 cannot occur, and t ≥ 2^256 needs (z − x_W − x_T) mod p < C0^2 + 3C0. In the loops,
the carry of a plain fold or of these folds into limbs 2 and 3 runs as a rare exit (section 2.5); the results are the same.

### 2.3 fe_sqr (64)

- **Off-diagonal sum** O = Σ_(i<j) a_i·a_j·2^(64(i+j)) (12 + 8). Six products. Columns 1 to 6 hold 1, 2, 3, 3, 2, 1 words, so six two-word additions in two chains: chain A adds L03+L12 (column 3), H03+H12 (column 4) and its carry into H13; chain B adds H01+L02, then H02, L13, the incremented H13 with L23, and its carry into H23, which is o6. As in fe_mul, the carry-only additions land on untouched high words, and o6 is the top limb of O < 2^448.
- **Doubling** (6). 2·O ≤ a^2 < 2^512, so no bit is lost: LSR d7 = o6 >> 63 and five EXTR #63 give d2 to d6. The low limb d1 = 2·o1 mod 2^64 is not formed: the first diagonal addition is ADDS z1, H00, o1, LSL #1, and the bit shifted out of o1 is bit 0 of d2, already taken by the EXTR.
- **Diagonals** (8 + 7). a_i^2 by MUL and UMULH, then ADDS z1, five ADCS (z2 to z6) and ADC z7 = d7 + H33 + c, the top limb of a^2 < 2^512, which cannot overflow.
- **Reduction** (23).

**64 = 20 + 8 + 6 + 7 + 23.**

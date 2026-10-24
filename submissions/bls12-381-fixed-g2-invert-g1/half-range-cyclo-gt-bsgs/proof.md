# Half-range BSGS on cyclotomic compressed GT keys

Parent of this candidate is the public template `6423c2728dc03d5f2d647c7f22adf1dac0643a66`. The procedure below is the whole attack. Q is the challenge point, fixed for that instance and not chosen by the solver. There is no amortization across future targets and no auxiliary powers.

## Group and inversion

Let p and r be the profile field modulus and subgroup order. GT is the order-r subgroup of Fp12× in the monomial basis of

```text
w^12 - 2 w^6 + 2 = 0
```

with coefficient list `[c0, …, c11]`. The declared experiment evaluates `w^{p^6}` in that quotient and obtains `w^{p^6} = -w`. Coefficients lie in Fp, which Frobenius fixes, so

```text
(Σ ci w^i)^{p^6} = Σ ci (-1)^i w^i.
```

Odd-degree coefficients are negated modulo p. Even-degree coefficients stay. Call this map `π^6`. On twelve coefficients it is six base-field negations, not a multiplication.

Independently, `p^6 ≡ -1 (mod r)`, checked as an integer power in the experiment. Every valid target satisfies `z^r = 1`, so

```text
z^{p^6} = z^{-1} = π^6(z).
```

The same identity holds for every power of `g = e(P0, Q)`. Nondegeneracy of the pinned pairing gives `e(P0, Q0) ≠ 1`. For `Q = [b] Q0` with `b ∈ [1, r)`, `g = e(P0, Q0)^b` still has order r, so g generates GT. The unknown is the unique `x ∈ [0, r)` with `z = g^x`, and the returned point is `P = [x] P0`. Bilinearity yields `e(P, Q) = z`. The identity `z = 1` forces `x = 0` and `P = ∞`, which the search returns as the baby `j = 0`.

## Cyclotomic key

Write a GT element as even part plus odd part:

```text
a = c0 + c2 w^2 + c4 w^4 + c6 w^6 + c8 w^8 + c10 w^10
b = c1 w + c3 w^3 + c5 w^5 + c7 w^7 + c9 w^9 + c11 w^11
α = a + b
π^6(α) = a - b
```

The fixed set of the automorphism `w ↦ -w` is a field F with `p^6` elements: it is the fixed field of a field automorphism of order 2. The even part a lies in F. The odd part is `b = w · b0` with `b0 ∈ F`. Because `w^2` is fixed by `w ↦ -w`, `b^2 ∈ F`. From `α · π^6(α) = 1`,

```text
a^2 - b^2 = 1    in F.
```

F is a field, so `F[X]` is a unique-factorisation domain and `X^2 - (1 + b^2)` has at most two roots. Those roots are `a` and `-a`, and they coincide only when `a = 0`.

The stored key is not a hash and not a truncation of a free 12-tuple. It is:

- the six odd coefficients of whichever of `b` and `-b` is lexicographically smaller as six 48-byte big-endian integers in `[0, p)`;
- one sign byte in `{0, 1}` for the even part.

The sign byte is 0 when `a = 0`. Otherwise it is the least significant bit of the first nonzero even coefficient in the fixed order `c0, c2, c4, c6, c8, c10`, read from byte 47 of that coefficient's 48-byte big-endian representative in `[0, p)`. It is not a fingerprint of a. Given b, at most two even parts exist, and this bit separates them: p is odd, so if `c ∈ (0, p)` then the canonical residue of `-c` is `p - c`, whose least significant bit is the opposite of the least significant bit of c. Negation does not change which coefficient is the first nonzero one. Therefore `a` and `-a` receive different sign bytes whenever `a ≠ 0`.

`π^6` keeps a and replaces b by `-b`, so α and `α^{-1}` have the same sign byte and the same canonical odd part. Distinct unordered pairs `{α, α^{-1}}` inside the norm-1 set produce distinct keys: equal canonical odd parts mean the odd parts agree up to sign, and equal sign bytes then force the even parts to agree. The experiment enumerates the same key on the norm-1 set of the prime field 1009, where the same degree-2 count is visible at toy scale. The production argument uses only that F is a field.

No square root is computed. The baby and the probe already hold every coefficient, because the group operation writes them. The sign byte is a walk over the six even coefficients.

## Search

Set `M = 2^126`, radius `ρ = M - 1`, stride `S = 2ρ = 2(M - 1)`.

```text
n = ceil((r - 1) / S)
```

Baby table, `j = 0 … M - 1`:

```text
acc ← 1
store canonical(acc, j = 0)
for j = 1 … M - 1:
    acc ← acc · g
    store canonical(acc, j)
Γ ← acc · acc          # g^{2ρ} = (g^ρ)^2
D ← π^6(Γ)             # g^{-S}
```

`canonical(α, j)` negates the six odd coefficients, compares that 288-byte image with the odd part of α, and stores the lexicographic minimum. The sign byte of the even part is stored beside it. The signed exponent is `s = +j` or `s = -j` according to which odd part won. The stored key is the canonical key of `g^s`, which equals the canonical key of `g^{-s}`. Distinct `j, j' ∈ [0, M)` give distinct unordered pairs `{j, -j}` because `j + j' ≤ 2ρ < r` and g has order r. Keys do not collide. The table has `M` records.

Giant loop:

```text
γ ← z
for i = 0 … n - 1:
    probe the canonical key of γ
    on a hit with stored signed s:
        if the odd part of γ is the stored odd part, x ← i·S + s
        else x ← i·S - s
        return [x mod r] P0 after the verification below
    γ ← γ · D
```

One probe negates six odd coefficients, selects the smaller odd part, appends the even sign byte, and binary-searches that key. Which side won is the sign branch. It is not a second search.

## Coverage

For `i ∈ [0, n)` and `s ∈ [-ρ, ρ]` the integers `i·S + s` are exactly the interval

```text
[-ρ, (n - 1)·S + ρ] = [-ρ, (2n - 1)·ρ].
```

Adjacent blocks abut: `i·S + ρ + 1 = (i + 1)·S - ρ`. The interval contains

```text
2·n·ρ + 1
```

consecutive integers. The choice of `n` is `2·n·ρ ≥ r - 1`, so the length is at least `r`. Every residue modulo `r` appears. The interval lies in `(-r, 2r)` because `ρ < r` and `(2n - 1)·ρ < r + ρ`, so reduction into `[0, r)` adds or subtracts `r` at most once.

Thus some pair `(i, s)` satisfies `x ≡ i·S + s (mod r)` with `s = ±j` and `0 ≤ j < M`. At that giant index, `γ = g^{x - i·S} = g^s`, and the canonical key of `γ` equals the stored key of `j`. The odd-part branch recovers `s`. Every valid challenge is solved. There are no algorithm coins. Success probability is 1 on the profile distribution, hence at least 0.39.

Invalid input is rejected before the search: `Q` on G2, not infinity, `[r]Q = ∞`, `z ≠ 0`, `z^r = 1`. After recovery the attack checks `[r]P = ∞` and `e(P, Q) = z`. Those checks are charged below. A failing check does not occur for a valid challenge; the budget still pays for them.

## Instruction schedules

Reference operation: one dense Fp12 multiplication, `C = 20000` scalar AArch64 instructions. Group multiplications stay dense reference operations. Compression changes only the lookup record. Each schedule below is the straight-line worst case, then doubled. The factor two is an allowance for address setup, renames and the loop exit, not a second copy of the group operation. Reference multiplications are not also unpacked into instructions.

Odd part: 6 coefficients × 48 bytes = 288 bytes. Sign byte: 1. Compare window: 304 bytes (19 chunks of 16, the extra 15 bytes zero). Record: that window plus a 16-byte signed exponent, 320 bytes. A signed 128-bit exponent holds `±j` for `j < 2^126`.

Comparison of the 304-byte window, 19 chunks:

```text
ldp, ldp, cmp, ccmp, b.ne, subs, b.ne
```

`7 × 19 + 4 = 137` instructions. Published comparison `B = 274`.

Copy of one 320-byte record, 20 chunks:

```text
ldp, stp, subs, b.ne
```

`4 × 20 + 12 = 92` instructions. Published copy `D = 184`.

Binary-search address arithmetic. The record stride is `320 = 2^8 + 2^6`, a shorter shift-add than a 592-byte stride. The published bound stays `ADDR = 64`, which remains an upper bound.

Odd-coefficient negation, six Fp elements, six 64-bit limbs each, eight instructions per limb: 288, published `CONJ = 576`.

Even-part sign, worst case, every limb inspected, no early exit counted as a saving. Six coefficients, six limbs:

```text
ldr, cmp, b.ne
```

`3 × 6 × 6 = 108`, then `ldrb, tst, cset, mov` for the selected low bit: 4. Raw total 112. Published `SIGN = 224`.

Per-record merge control: 48, published `CTRL = 96`.

Derived:

```text
CANON = CONJ + SIGN + B + D = 1258
PASS  = B + 2·D + CTRL = 738
V     = CANON + (L + 1)·(B + ADDR) = 44184
L     = 126
```

`CANON` is one odd-part negation, one even-part sign, one key comparison and one worst-case copy of the winning record. It is charged once per baby while the table is built, and once per giant probe. `V` is that work plus `L + 1` binary-search comparisons. A bottom-up mergesort of `2^126` records uses `L = 126` passes. Each pass compares each record at most once and copies it at most twice.

## Operation counts

Baby-side dense multiplications: `M - 1` steps from `g^0` to `g^ρ`, plus one squaring to `g^S`. That is `M` multiplications. Giant-side: `n - 1` multiplications by `D`.

```text
M     = 2^126 = 85070591730234615865843651857942052864
S     = 170141183460469231731687303715884105726
n     = 308190375243917312857644309053803048975
H_mul = M + n - 1 = 393260966974151928723487960911745101838
```

Other dense multiplications, each a hard upper bound: the setup pairing `e(P0, Q)` at most 4096, the verifying pairing `e(P, Q)` at most 4096, and `z^r` together with the input exponent check at most 1024. Miller loop length is the 64-bit parameter `x`; charging every squaring and sparse product as a dense product, plus a final exponentiation of a few thousand dense products, fits in 4096. Sum `H_misc = 9216`.

```text
H = H_mul + H_misc = 393260966974151928723487960911745111054
```

Instruction totals:

| Piece | Instructions |
| --- | ---: |
| Mergesort `L·M·PASS` | 7910544183811056460133069498966315611717632 |
| Baby canonicalization `M·CANON` | 107018804396635146759231314037291102502912 |
| Baby stores `M·D` | 15652988878363169319315231941861337726976 |
| Giant probes `n·V` | 13617083539777242551302156151233233915911400 |
| Stride conjugation plus G1/G2 scalar reserve `CONJ + 2^26` | 67109440 |

The reserve pays for `[x]P0`, `[r]P` and `[r]Q`. Each is at most 400 point operations. A G1 operation is charged at 8192 instructions and a G2 operation at 16384, together under `2^25`, inside the `2^26` reserve. `R` is the sum of the five rows:

```text
R = 21650299516863297327513772196178702034968360
ceil(R / C) = 1082514975843164866375688609808935101749
W = H + ceil(R / C) = 1475775942817316795099176570720680212803
```

`log2(W)` is `130.116669403…`. The submitted score is the ceiling to six decimal places, `130116670 / 1000000`.

Preprocessing is the pairing that produces `g`, the baby multiplications, canonicalization, stores, mergesort and the single stride conjugation. It excludes the giant loop and the verifying pairing.

```text
R_pre = 8033215977086054776211616044945468052996672
W_pre = M + 4096 + ceil(R_pre / C) = 486731390584537354676424454105215459610
```

`log2(W_pre) = 128.516393426…`, reported as `128516394 / 1000000`. This quantity is inside `W`, not a second objective.

Peak memory is two sort buffers of `M` records plus `2^20` bytes of code, field scratch and the probe buffer:

```text
2·M·320 + 2^20 = 54445178707350154154139937189082914881536
```

`log2` of that byte count is `135.321928094…`, reported as `135321929 / 1000000`. Advice stored in the algorithm is empty. Public profile constants are not advice. `log2(max(1, 0)) = 0`.

No heuristic is used. The Frobenius identity is a finite calculation in the pinned polynomial, repeated by `experiment.py`. The degree bound `X^2 - c` is the statement that a field polynomial of degree 2 has at most two roots. The toy exhausts the same stride, canonical odd key and sign branch on `Z/nZ`, and exhausts the norm-1 key in `F_1009`. Neither prices production arithmetic.

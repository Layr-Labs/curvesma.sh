# Half-range BSGS by Frobenius-6 inversion

Parent of this candidate is the public template `6423c2728dc03d5f2d647c7f22adf1dac0643a66`. The procedure below is the whole attack. Q is the challenge point, fixed for that instance and not chosen by the solver. There is no amortization across future targets and no auxiliary powers.

## Group and inversion

Let p and r be the profile field modulus and subgroup order. GT is the order-r subgroup of Fp12× in the monomial basis of

```text
w^12 - 2 w^6 + 2 = 0
```

with coefficient list `[c0, …, c11]`. The declared experiment evaluates `w^{p^6}` in that quotient and obtains `w^{p^6} = -w`. Consequently

```text
(Σ ci w^i)^{p^6} = Σ ci (-1)^i w^i.
```

Odd-degree coefficients are negated modulo p. Even-degree coefficients stay. Call this map `π^6`. It is six base-field negations, not a multiplication.

Independently, `p^6 ≡ -1 (mod r)`, also checked as an integer power in the experiment. Every valid target satisfies `z^r = 1`, so

```text
z^{p^6} = z^{-1} = π^6(z).
```

The same identity holds for every power of `g = e(P0, Q)`. Nondegeneracy of the pinned pairing gives `e(P0, Q0) ≠ 1`. For `Q = [b] Q0` with `b ∈ [1, r)`, `g = e(P0, Q0)^b` still has order r, so g generates GT. The unknown is the unique `x ∈ [0, r)` with `z = g^x`, and the returned point is `P = [x] P0`. Bilinearity yields `e(P, Q) = z`. The identity `z = 1` forces `x = 0` and `P = ∞`, which the search returns as the baby `j = 0`.

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

`canonical(α, j)` writes the 576-byte monomial encoding of `π^6(α)` into a scratch buffer (odd coefficients negated), compares it with the encoding of `α`, and stores the lexicographic minimum together with the signed exponent `s = +j` or `s = -j` according to which encoding won. The stored key is exactly the encoding of `g^s`. Distinct `j, j' ∈ [0, M)` give distinct unordered pairs `{j, -j}` because `j + j' ≤ 2ρ < r`. Keys do not collide. The table has `M` records.

Giant loop:

```text
γ ← z
for i = 0 … n - 1:
    probe the canonical key of γ
    on a hit with stored signed s:
        if γ has that key as its own encoding, x ← i·S + s
        else x ← i·S - s
        return [x mod r] P0 after the verification below
    γ ← γ · D
```

One probe builds `π^6(γ)`, selects the smaller encoding, and binary-searches it. The sign test is one extra encoding compare, charged inside the per-probe canonicalization below, not a second search.

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

Thus some pair `(i, s)` satisfies `x ≡ i·S + s (mod r)` with `s = ±j` and `0 ≤ j < M`. At that giant index, `γ = g^{x - i·S} = g^s`, and the canonical key of `γ` equals the stored key of `j`. The sign branch recovers `s`. Every valid challenge is solved. There are no algorithm coins. Success probability is 1 on the profile distribution, hence at least 0.39.

Invalid input is rejected before the search: `Q` on G2, not infinity, `[r]Q = ∞`, `z ≠ 0`, `z^r = 1`. After recovery the attack checks `[r]P = ∞` and `e(P, Q) = z`. Those checks are charged below. A failing check does not occur for a valid challenge; the budget still pays for them.

## Instruction schedules

Reference operation: one dense Fp12 multiplication, `C = 20000` scalar AArch64 instructions. Each schedule below is the straight-line loop, then doubled. The factor two is an allowance for address setup, renames and the loop exit, not a second copy of the group operation. Reference multiplications are not also unpacked into instructions.

Monomial encoding of one GT element: 12 coefficients × 48 bytes = 576 bytes. Record: that key plus a 16-byte signed exponent, 592 bytes.

Comparison of two 576-byte keys, 36 chunks of 16 bytes:

```text
ldp, ldp, cmp, ccmp, b.ne, subs, b.ne
```

7 × 36 + 4 = 256 instructions. Published comparison `B = 512`.

Copy of one 592-byte record, 37 chunks:

```text
ldp, stp, subs, b.ne
```

4 × 37 + 12 = 160 instructions. Published copy `D = 320`.

Binary-search address arithmetic (midpoint, multiply by 592 via shifts and adds, base add, branch): 32 instructions, published `ADDR = 64`.

Odd-coefficient negation, six Fp elements, six 64-bit limbs each, eight instructions per limb: 288, published `CONJ = 576`.

Per-record merge control: 48, published `CTRL = 96`.

Derived:

```text
CANON = CONJ + B + D = 1408
PASS  = B + 2·D + CTRL = 1248
V     = CANON + (L + 1)·(B + ADDR) = 74560
L     = 126
```

`CANON` is one conjugation, one key comparison and one worst-case copy of the winning encoding. It is charged once per baby while the table is built, and once per giant probe. `V` is that work plus `L + 1` binary-search comparisons. A bottom-up mergesort of `2^126` records uses `L = 126` passes. Each pass compares each record at most once and copies it at most twice (the qualified word-operation ledger counted two record touches per pass; the same upper bound is kept here).

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
| Mergesort `L·M·PASS` | 13377180408395932875672182567357671928758272 |
| Baby canonicalization `M·CANON` | 119779393156170339139107861815982410432512 |
| Baby stores `M·D` | 27222589353675077077069968594541456916480 |
| Giant probes `n·V` | 22978674378186474846665959683051555331576000 |
| Stride conjugation plus G1/G2 scalar reserve `CONJ + 2^26` | 67109440 |

The reserve pays for `[x]P0`, `[r]P` and `[r]Q`. Each is at most 400 point operations. A G1 operation is charged at 8192 instructions and a G2 operation at 16384, together under `2^25`, inside the `2^26` reserve. Encoding compares sit in the same reserve. `R` is the sum of the five rows:

```text
R = 36502856769092253138554320080819751194792704
ceil(R / C) = 1825142838454612656927716004040987559740
W = H + ceil(R / C) = 2218403805428764585651203964952732670794
```

`log2(W)` is `130.704717696…`. The submitted score is the ceiling to six decimal places, `130704718 / 1000000`.

Preprocessing is the pairing that produces `g`, the baby multiplications, canonicalization, stores, mergesort and the single stride conjugation. It excludes the giant loop and the verifying pairing.

```text
R_pre = 13524182390905778291888360397768195797156416
W_pre = M + 4096 + ceil(R_pre / C) = 761279711275523530460261671746351846818
```

`log2(W_pre) = 129.161694235…`, reported as `129161695 / 1000000`. This quantity is inside `W`, not a second objective.

Peak memory is two sort buffers of `M` records plus `2^20` bytes of code, field scratch and the probe buffer:

```text
2·M·592 + 2^20 = 100723580608597785185158883799803391639552
```

`log2` of that byte count is `136.209453365…`, reported as `136209454 / 1000000`. Advice stored in the algorithm is empty. Public profile constants are not advice. `log2(max(1, 0)) = 0`.

No heuristic is used. The Frobenius identity is a finite calculation in the pinned polynomial, repeated by `experiment.py`. The toy exhausts the same stride, canonical key and sign branch on `Z/nZ`; it does not price production arithmetic.

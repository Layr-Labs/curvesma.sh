# Exact geometric Cheon on the supplied d-th power

Parent of this candidate: `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

## Claim

The track supplies a prime-order point `P ∈ G1`, a target `T = [s]P`, and one auxiliary point `A = [s^d]P`, where the positive integer `d` is given by the profile and satisfies `d | (r − 1)`. Here `r = #G1` is the BLS12-381 scalar prime

```text
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
```

The algorithm below returns `s ∈ {0, 1, …, r − 1}` with success probability 1. Every group operation is a point addition, a point doubling, or a scalar multiplication built out of those, and every one of them is charged. There is no multiplication-table estimate, no Erdős lower bound, and no uncovered residue class.

## Algebra

Work in `F_r`. Let `m = (r − 1) / d` and let `H = { x^d : x ∈ F_r^* }`, the unique subgroup of order `m`. The unknown scalar `σ = s^d mod r` lies in `H ∪ {0}`.

If `T` is the identity, return `s = 0`. Then `A` is also the identity. Conversely if `A` is the identity and `T` is not, the instance is inconsistent with the profile (`s^d ≡ 0` forces `s ≡ 0`).

Fix the smallest primitive root `g₀` modulo `r` and set `γ = g₀^d mod r`. Then `γ` generates `H`. Hence there is a unique exponent `e ∈ {0, 1, …, m − 1}` with

```text
σ = γ^e    in F_r.
```

Let `G = ceil(sqrt(m))`. Divide `e` by `G`:

```text
e = i + G · j,    0 ≤ i < G,    0 ≤ j ≤ floor((m − 1) / G).
```

The integer `floor((m − 1) / G)` is at most `G − 1`, because `G² ≥ m`. Every exponent in `{0, 1, …, m − 1}` occurs exactly once. Exponents `≥ m` are not queried.

The group identity used is purely multiplicative in the scalar:

```text
A = [σ]P = [γ^(i + G j)]P = [γ^i] ( [γ^(G j)] P ).
```

Equivalently, with `ι = γ^(−1) mod r`,

```text
[γ^(G j)] P = [ι^i] A.
```

Both sides are geometric orbits. They are not arithmetic progressions, and the equality `γ^(G j) = j · γ^G` is not used. A match of points gives the unique pair `(i, j)`, hence the unique scalar `σ = γ^(i + G j)`.

## Lifting σ through μ_d

The polynomial `X^d − σ ∈ F_r[X]` splits completely: its roots are exactly `t · ζ^k` for `k = 0, …, d − 1`, where `ζ = g₀^m mod r` is the canonical element of order `d` and `t` is any single root. Exactly one root equals `s`, because `T = [s]P` distinguishes the coset (the map `c ↦ [c]P` is injective on `F_r`).

A root `t` is computed in `F_r` by the deterministic 2-power tower (and an odd-order inversion if `d` is not a pure power of two), using only field arithmetic. No group operation is hidden in that tower. Each of the `d` candidates is then tested by a scalar multiplication of the fixed base `P`. The orbit is charged in full: `d` scalar multiplications, not one.

## Procedure

1. Reject the identity cases above.
2. Compute `g₀`, `γ`, `ι = γ^(−1)`, `ζ`. Field cost `O((log r)^3)`, absorbed in the explicit bound below.
3. Baby orbit, base `A`, ratio `ι`. Set `B₀ = A` and `B_(i+1) = [ι] B_i` for `i = 0, …, G − 2`. Store the map from the affine encoding of `B_i` to `i`. Building the orbit takes `G − 1` scalar multiplications by the fixed scalar `ι`.
4. Stride point. Compute `Q₁ = [γ^G] P` with one scalar multiplication. Set `Q_(j+1) = [γ^G] Q_j`. Each of the `floor((m − 1) / G)` giant steps is one scalar multiplication by the fixed scalar `γ^G`. Look `Q_j` up in the baby map. The first hit is the unique legal pair `(i, j)`.
5. Set `σ = γ^(i + G j) mod r`. Compute one `d`-th root `t` of `σ`. For `k = 0, …, d − 1` test `[t · ζ^k] P = T`. Return that scalar.

Step 4 multiplies by `γ^G`, not by `j`. The index `j` only counts how many times that fixed scalar has been applied. Coverage is the divisibility argument in the previous section, not a probability.

## What is not claimed

- Table entries are not priced as point additions. A step `[ι]B` is a scalar multiplication.
- Distinct keys are not borrowed from an integer multiplication table. The keys are the points `B_i`, and there are exactly `G` of them.
- No factor `1/1.5`, no coverage floor `0.39`, no `λ` haircut. The rectangle `{i + G j}` covers `{0, …, m − 1}` exactly.
- The recovered word `i + G j` is an exponent of `γ`, not a candidate for `s`. The only scalars substituted into `[c]P` are the `d` roots of `X^d − σ`.

## Cost shape

Let `C_×` be the reference price of one G1 scalar multiplication by a full-size scalar, expressed in the track’s addition/doubling/field unit, using a single fixed window (wNAF, window 4) and the track’s Jacobian formulas. Let `C_field` be a bound on the scalar arithmetic of steps 2 and 5 (primitive root, one exponentiation tower, `d` multiplications by `ζ`).

Charged group operations:

```text
baby    = G − 1
giant   = floor((m − 1) / G)
stride  = 1
orbit   = d
```

```text
W = (baby + giant + stride + orbit) · C_× + C_field
```

With `G = ceil(sqrt(m))` one has `baby + giant + stride ≤ 2 G`, and `orbit = d` is negligible beside `G` whenever `d ≤ 2^32` and `m ≥ 2^192`, which is the BLS12-381 range for every profile divisor `d | (r − 1)` with `d ≤ 2^32`. The published `time_log2` is `log2(W)` after `C_×` is instantiated from the track price list. It is an upper bound on the operations the procedure performs, not an expectation.

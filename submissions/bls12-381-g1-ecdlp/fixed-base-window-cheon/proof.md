# Fixed-base windowed Cheon on the subgroup of d-th powers

## Claim

Let `r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001` be the BLS12-381 G1 order and let `P0` be the standard generator. Sample `s` uniformly from `1..r-1` and publish `n = 2^20` together with `A_i = [s^i mod r] P0` for `i = 1..n`. The algorithm below returns that `s`. It is deterministic, succeeds on every such input, and at `n = 2^20` uses

```text
W(2^20) < 2^125
```

target-group-operation equivalents in the sense of `classical-group-work-v2` with `C = 2000`. Peak storage is at most `2^126` bytes. The same quantity bounds preprocessing, which is entirely contained in `W`. Advice embedded in the algorithm is empty.

The function statement used by the score is the one-point evaluation `W(2^20)` above. The algorithm requires `n | (r-1)` and the single auxiliary point `A_n`; both hold at the scoring parameter because `v2(r-1) = 32 ≥ 20`. Intermediate powers `A_i` for `1 < i < n` are not read. Generation of the powers is organizer work and is not charged. Reading `A_1` and `A_n` is charged inside the lump below.

## Group arithmetic unit

The reference operation is one batched incomplete affine G1 addition, normalized as `6M` with `S = M` and `A = 0`, priced at `C = 2000` instruction equivalents. Field work that is not inside a charged point operation costs its AArch64 instruction count divided by `C`. A point operation is never charged both as a reference unit and as its internal multiplications.

Short-Weierstrass Jacobian arithmetic on `y^2 = x^3 + 4` is enough:

- one doubling costs at most `4M + 6S ≤ 10M ≤ 2` reference units;
- one mixed addition (Jacobian accumulator plus affine table point) costs at most `12M + 6S ≤ 18M ≤ 3` reference units.

Published `a = 0` formulas are strictly cheaper (doubling `2M+5S`, mixed addition `11M+5S`). The inequalities above are what the bound uses. Dependent chains stay in Jacobian, so no per-step inversion appears. One Montgomery batch inversion normalizes a finished table; it is costed with the field term, not as `T` inversions.

## Integers

```text
d  = 2^20 = 1048576
m  = (r-1)/d
   = 50006747412801924209068050869165387952509453297164571592906626415194112
N  = 223621884914696853578958571250384896
```

Integer relations used later, all strict where stated:

1. `(N-1)^2 < m ≤ N^2`, so `N = ceil(sqrt(m))`. Also `ceil(m/N) = N`.
2. `1000 N < 1347 · 2^117`, so `N < (1347/1000) 2^117`.
3. `2^25 < N`.
4. `119 · 1347 = 160293 < 256000 = 256 · 1000`, so `(119/1)·(1347/1000) < 2^8`.

The prime factorization

```text
r-1 = 2^32 · 3 · 11 · 19 · 10177 · 125527 · 859267 · 906349^2
      · 2508409 · 2529403 · 52437899 · 254760293^2
```

multiplies back to `r-1`. Consequently every prime `q | (r-1)` is on that list. Direct exponentiation gives `7^{(r-1)/q} ≢ 1 (mod r)` for each of them, and `7^(r-1) ≡ 1 (mod r)`, so `7` is a primitive root modulo `r`. Let `ω = 7^d mod r`. Then `ω` has order `(r-1)/gcd(d, r-1) = m`.

## What is recovered

Write the unknown as `s = 7^L` with `L ∈ {0,1,...,r-2}` unique. Then

```text
s^d = 7^{d L} = ω^{L mod m},
```

because `d · m = r-1` kills the multiple of `m` inside `L`, and `d · (L mod m) < r-1`. The element `s^d` therefore lies in the unique subgroup `<ω>` of order `m`, and the discrete log

```text
λ = L mod m ∈ {0,1,...,m-1}
```

is the unique solution of `ω^λ ≡ s^d (mod r)`. The point supplied by the challenge is `G_d = A_d = [s^d] P0 = [ω^λ] P0`. No step of this reduction multiplies the point `G_d` by the integer `m/q` and calls the result `[ω^{m/q}] P0`. Those two scalars are different: `[a][σ]P0 = [a·σ]P0`, whereas `[σ^a]P0` is a field power of the hidden scalar. The search below only multiplies by scalars that the algorithm has written down.

## Phase 1: BSGS for `λ`

Every `λ ∈ [0, m)` factors as `λ0 + N λ1` with `λ0, λ1 ∈ {0,1,...,N-1}`, because `N^2 ≥ m`. Reduce the sum modulo `m` if a later collision lands above `m-1`; the residue is unique. The matching equation is

```text
[ω^{-λ0}] G_d = [(ω^N)^{λ1}] P0.
```

If the points are equal then their scalars agree modulo `r`:

```text
ω^{-λ0} · s^d ≡ (ω^N)^{λ1} (mod r).
```

Substitute `s^d ≡ ω^λ` and cancel inside `<ω>` to obtain `λ ≡ λ0 + N λ1 (mod m)`.

Both legs are fixed-base scalar multiplications.

- Baby leg, base `G_d`. Stream `α_0 = 1`, `α_{j+1} = α_j · ω^{-1} mod r`, and store the affine point `[α_j] G_d` under key `j`, for `j = 0..N-1`.
- Giant leg, base `P0`. Stream `β_0 = 1`, `β_{k+1} = β_k · (ω^N) mod r`, and look up `[β_k] P0` in the baby table, for `k = 0..N-1`.

`ω^N mod r` and `ω^{-1} mod r` are two field exponentiations, charged in the lump. Each streamed multiplier is one multiplication modulo `r`.

### Window evaluation

Scalars are integers in `1..r-1`, hence shorter than `256` bits. Fix window width `16`. For a base `B` precompute, once, the affine table

```text
T[j][u] = [u · 2^{16 j}] B,    j = 0..15,  u = 1..2^16-1.
```

Build the `j`-th row by doubling out to `[2^{16 j}] B` and then repeatedly mixed-adding that affine (or once-normalized) step. Row length `2^16-2` additions after the first entry. Evaluation of one scalar adds at most `16` table points and therefore performs at most `15` mixed additions. Zero windows are skipped, so `15` is a worst case.

Two bases, `N` scalars each:

```text
H_eval = 15 · 3 · (N+N) = 90 N.
```

Table construction itself is

```text
2 bases · 16 rows · (2^16 - 2) · 3  <  2^23
```

reference units, plus fewer than `2^12` doublings. That sits in the lump.

## Phase 2: the `2^20` coset

With `λ` known, set `σ = ω^λ mod r` by a field exponentiation (`λ < m < 2^235`). Every `d`-th root of `σ` is

```text
s_κ = 7^{λ + m κ} mod r,    κ = 0,1,...,d-1,
```

and these `d` values are distinct because `7^m` has order `d`. Exactly one of them equals `s`, namely `κ = ⌊L/m⌋`. Equivalently, writing `η = 7^m mod r` and `S0 = [7^λ] P0`,

```text
A_1 = [η^κ] S0.
```

Split `κ = κ0 + 1024 κ1` with both digits in `0..1023` (`d = 1024^2`). The same windowed fixed-base BSGS, two legs of length `1024`, costs

```text
15 · 3 · 2048 = 92160
```

reference units. The final check recomputes `s = 7^{λ + m κ} mod r` and the single scalar multiplication `[s] P0`, at most `255` doublings and `255` mixed additions, under `2000` reference units. Accept only if the product equals `A_1`. Uniqueness of `κ` makes the check succeed on the first algebraic solution; there is no restart and no failure branch on a well-formed challenge.

## Lookup, field work, lump

Baby keys are the affine `x`-coordinates, `6` limbs of `64` bits. Sort each leg with mergesort and probe the baby table by binary search. `ceil(log2 N) = 118` because `2^117 < N < 2^118`.

- Comparisons: at most `118 N` to sort a leg and `118 N` to probe it. Both legs: `472 N` comparisons. At `64` instructions each (six-limb compare plus loop control), `472 · 64 / 2000 = 15.104`, charged as `16 N`.
- Moves: each element is moved once per mergesort level, `118` levels, both legs, `236 N` moves. A `96`-byte point move is under `48` instructions. `236 · 48 / 2000 = 5.664`, charged as `8 N`.
- Field stream: `2N` geometric multiplications and `2N` digit extractions, each bounded by `2000` instructions, contribute `4 N`. All other exponentiations (`7^d`, `ω^N`, `ω^{-1}`, `7^λ`, `η^κ`, the final `7`-power, batch inversion of the window tables) are `O(2^20)` field multiplications or a single Fermat inversion plus `3(T-1)` multiplications on a table of size `< 2^22`. Charged as `< 2^24` reference units together with the window build, phase 2, and the final scalar.

## Sum

```text
W(2^20)  ≤  90 N + 16 N + 8 N + 4 N + 2^24
         =  118 N + 2^24
         <  118 N + 2^25
         <  119 N
```

because `2^25 < N`. Relation 2 then gives

```text
119 N  <  119 · (1347/1000) · 2^117.
```

Relation 4 says `119 · 1347 < 256 · 1000`, hence `119 · (1347/1000) < 2^8` and

```text
W(2^20)  <  2^8 · 2^117  =  2^125.
```

The score is exactly this upper bound: `time_log2 = 125/1`. Review is not asked to tighten it. Every term above is deterministic worst case, so the inequality does not depend on a Hamming-weight heuristic or on the birthday constant of a random walk.

## Memory

The baby table stores `N` affine points (`96` bytes) and an `8`-byte index: `104 N` bytes. Mergesort keeps a second point buffer, `96 N` bytes. The giant leg is streamed. Both window tables together are `2 · 16 · 2^16 · 96 < 2^28` bytes, and `2^28 < 9 N`, so the live set is at most `104 N + 96 N + 9 N = 209 N` bytes including the phase-2 table and the process image. Since `209 · 1347 < 512 · 1000` and `N < (1347/1000) 2^117`,

```text
209 N  <  2^9 · 2^117  =  2^126.
```

`memory_log2_bytes = 126`. No searched table is embedded in the algorithm. Advice length is `0`, and `log2(max(1, 0)) = 0`.

## Success

On every `s ∈ [1, r-1]` there is a unique pair `(λ, κ)` and both BSGS searches cover that pair. Point equality in a prime-order group is scalar equality modulo `r`, and scalar equality inside `<ω>` is equality of discrete logs modulo `m`. The final multiplication `[s] P0 = A_1` rejects a programming error and accepts the mathematical solution. Success probability is `1` on the uniform production distribution. No algorithm coins are drawn. The `0.39` threshold is met.

## What this does not do

Pohlig–Hellman on the prime factors of `m` would need the point `[σ^{m/q}] P0` with `σ = s^d` unknown. Scalar multiplication by the public integer `m/q` produces `[ (m/q) σ ] P0`, not that field power. The algorithm never forms that point. Smoothness of `m` is not an input to the cost. The search length is `N = ceil(sqrt(m))` on each leg, full coverage, with every leg step expanded into at most `15` mixed additions.

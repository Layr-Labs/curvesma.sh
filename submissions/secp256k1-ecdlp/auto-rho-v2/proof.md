# Automorphism-aware Pollard rho v2

Status: ready for exploratory review. No qualification or emitted score is claimed
by this file. The normative input is a point `Q` on the protected secp256k1
profile; the required output is the scalar `x` with `Q=[x]P`. The cost model is
`classical-group-work-v2`. SHA-256 of the raw bytes of
`tracks/secp256k1-ecdlp/target-profile.json` is
`18c3338de676f1909715814325bffdb1513f1621d678b942e8a7a9563b8fd086`.
SHA-256 of the raw bytes of `tracks/classical-group-work-v2.json` is
`e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294`.

The reference operation is one batched incomplete affine addition, normalised at
`C=1000`. Field multiplication is the pinned 5×52 kernel at 179 instructions.
`S=M` and `A=0` apply only inside that reference formula. Scalar adds, the
endomorphism multiplications, canonicalisation and the hash chain are other
instruction work `R` and are divided by `C` once. A reference operation is never
also unpacked into its six multiplications.

Public curve constants are not advice. In particular
`β = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee`
satisfies `β²+β+1 ≡ 0 (mod p)`, and
`λ = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72`
satisfies `λ²+λ+1 ≡ 0 (mod r)`. Both checks are integer identities on the
protected `p` and `r`. `β²` and `λ²` are reduced once during preprocessing.

## Orbits

Let `ψ(x,y) = (βx, y)`. Together with negation this generates six maps

```text
(x, y) -> (β^e x, (−1)^s y),   e in {0,1,2}, s in {0,1}.
```

On the prime-order subgroup the only fixed point of a non-identity map would
require a multiplier `m ≢ 1 (mod r)` with `[m]R = R`, hence `R = ∞`. Every
finite point therefore has orbit size 6. The number of orbits is

```text
N = 1 + (r − 1)/6
  = 19298681539552699237261830834781317975472927379845817397100860523586360249057.
```

`r ≡ 1 (mod 6)`, so the division is exact. The canonical representative of a
finite point is the image with smallest `x` as an integer in `[0,p)`, and, for
that `x`, the sign whose `y` is smaller than `p−y`. The subgroup has no
order-2 point, so `y ≠ p−y` and the sign is unique. Three candidate
`x`-coordinates are `x`, `βx` and `β²x`.

## Coefficients in Z[λ]

A walk state stores the point `R = [α]P + [βs]Q` with `α` and `βs` written in
the basis `1, λ`:

```text
α  = u0 + v0 λ,    βs = u1 + v1 λ,    coefficients mod r.
```

Multiplication by a unit does not use a 256-bit product. With `λ² = −1−λ`,

```text
λ (u + vλ)  = −v + (u − v)λ
λ²(u + vλ)  = (v − u) − u λ
−(u + vλ)   = (−u) + (−v)λ
```

and the mixed signs are the compositions of those maps. Each is a constant
number of 4-limb additions and a conditional add or subtract of `r`. The
experiment applies these maps on `Z/n` and checks, on every step, that the
expanded scalar equals the canonical representative. Expanding `u+vλ` back to a
single residue is done only when a relation is solved, twice per collision.

## Walk

There are 4096 walks. Walk `w` starts at `[w+1]P + Q`, then replaces that point
by its canonical representative and multiplies both coefficient pairs by the
same unit. The 16 branch points are public,

```text
S_j = [a_j]P + [b_j]Q,   a_j = 2j+1,   b_j = 3j+5,   j = 0..15.
```

They are built by double-and-add. The branch index of a canonical point with
integer `x` is

```text
mixed = (x mod 2^64) * 0x9E3779B97F4A7C15   mod 2^64
j     = mixed >> 60.
```

One step, from an already canonical point `R`:

1. Read `j` from `R`. If `R` and `S_j` have the same `x`, or are either point
   at infinity, replace `j` by `j+1 mod 16`. If the second branch is also
   incomplete, abort this attempt (counted below as a charged failure).
2. Add `S_j` with one batched incomplete affine addition. The 4096 denominators
   of a round are one Montgomery batch. That is one reference operation per
   walk per round, which is the batch size the 6M price was calibrated on.
3. Canonicalise the sum. Multiply both coefficient pairs by the unit that was
   applied to the point. Add `(a_j, 0)` and `(b_j, 0)` in the basis (already
   included before the unit multiplication; order is add-branch, then unit).
4. Probe the chain for this canonical `x`. On a hit with a non-zero `β`
   denominator, solve for `x` and verify `[x]P = Q`. On a zero denominator,
   keep walking. On a miss, store the record if the chain has length under 64;
   if the insertion would make length 64, abort the attempt.

The stored record is the canonical `x` (32 bytes) and the four basis
coefficients (128 bytes), with padding and a link word to 192 bytes. The table
has `M = 4 T` slots, `T` the addition budget below, load `1/4`. Each probe
compares one `x` and follows one link. The step scans at most 64 probes and
pays for all 64 whether or not the chain is shorter.

## Birthday budget

For a random mapping on a set of size `N`, the walk `Y_{i+1}=f(Y_i)` satisfies

```text
P(no repeat through t steps) = prod_{i=1..t} (1 − i/N) ≤ exp(−t(t+1)/(2N)),
```

because the image of a fresh point is uniform on the whole codomain. The
inequality `t(t+1) ≥ 2 N ln 2` is therefore enough for collision probability at
least `1/2`. The series

```text
ln 2 < sum_{k=1..60} 1/(k 2^k) + 1/(61 · 2^60)
```

is an explicit rational upper bound. The least positive integer `t0` with
`t0(t0+1)` at least the ceiling of `2 N` times that rational is

```text
t0 = 163565440712060065318169743205629135679.
```

The charged walk is longer by the factor `5/4`, then rounded up to a multiple
of 4096 so every round is a full inversion batch:

```text
t_pad = ceil(5 t0 / 4)
      = 204456800890075081647712179007036419599
T     = ceil(t_pad / 4096) * 4096
      = 204456800890075081647712179007036420096.
```

Under the random-mapping model the collision probability at `t0` is already at
least `1/2`, and `T > t0`. The submitted success floor is `0.39`. The absolute
gap `0.11` is the allowance for the adding walk not being a uniformly random
function and for fruitless denominators. It is not a second objective and it is
not multiplied back into the runtime: failures inside this single budget are
already paid for, and there is no restart.

Starting the 4096 walks on distinct points only adds samples, so the serial
product bound remains an upper bound on the miss probability.

Chain-length abort, under an independent uniform hash into `M = 4T` slots: a
given slot has mean `μ = 1/4`. A binomial tail is at most `(eμ/64)^64`. Since
`e < 3`, this is `< (1/64)^64 = 2^{−384}`. The union over `M < 2^{130}` slots is
`< 2^{−254}`. Incomplete additions are a separate abort. Each addition has at
most two bad partners (equal or opposite `x`) out of `r`, and two branches are
tried, so the expected number of incomplete events on the whole budget is
`< 8T/r < 2^{−120}`. The attempt also aborts if that count exceeds `2^{20}`.
A Chernoff bound at that distance is far below `2^{−254}`. Neither abort moves
`0.50` down to `0.39`.

The toy uses the same `t0`, the same `5/4` pad, the same units and the same
branch index on `Z/n` with `n` prime and `n ≡ 1 (mod 6)`. It stores canonical
values in a dictionary rather than a 64-probe chain; the chain cap is an
instruction bound, not part of the relation. Author execution of the declared
input `order=1009` (`N_orb=169`, `t0=15`, padded budget `19`, 8 walks) recovered
`947/1009` secrets and reported zero basis-identity failures. That run is what
the verifier executes. It supports the relation and the shape of the budget at
that order only.

## Instruction ledger, per step

| Item | Instructions |
| --- | ---: |
| `β·x` and `β²·x`, each `179+16` | 390 |
| Three-way minimum of 4-limb `x` | 96 |
| Sign selection, including `p−y` | 64 |
| Incomplete-addition test on both branch attempts | 48 |
| Branch index, one 64-bit multiply and a shift | 32 |
| Two basis pairs through the worst unit (`−λ²`), 4-limb adds | 256 |
| Add `(a_j, b_j)` and reduce mod `r` | 64 |
| Hash formation plus 64 chain probes at 12 instructions | 784 |
| Round loop, pointers, batch index | 64 |
| Total `R_step` | 1798 |

The addition itself is one reference operation and is not in this table.
`H` counts every such addition, plus reserves that are not hidden inside `T`:

| Piece | Reference operations |
| --- | ---: |
| Walk additions `T` | 204456800890075081647712179007036420096 |
| Preprocessing reserve `2^16` | 65536 |
| Verification reserve `2^12` | 4096 |
| Incomplete-addition reserve `2^20` | 1048576 |
| Total `H` | 204456800890075081647712179007037538304 |

Preprocessing really uses at most `4096` adds to build `[1]P..[4096]P`,
`4096` adds of `Q`, and `16·2·512` double-and-add steps for the branch scalars,
which is `24592` reference operations. The `2^16` line covers that and the
initial canonicalisation of the 4096 starts. Verification is one double-and-add
for `[x]P` (at most 512 reference operations) plus encoding; the `2^12` line
covers a second comparison pass. The incomplete reserve is the abort cap above,
paid whether or not the cap is hit.

Other instructions outside the per-step loop:

```text
R_pre    = 2^24 = 16777216
R_verify = 2^22 = 4194304
R        = T * 1798 + R_pre + R_verify
         = 367613328000354996802586497854651504304128
```

`R_pre` covers bit scans and table allocation addressing for the `24592`
preprocessing additions. `R_verify` covers the final two expansions by `λ`,
the binary extended gcd for the denominator inverse, and output encoding.
`256` iterations of a 4-limb division fit inside `2^22`.

```text
W = H + R/1000
  = 572070128890430078450298676861689042608128 / 1000.
```

`log2(W) = 128.749459620557…`. Rounding upward at the sixth decimal place
gives the submitted total **128.749460**. The rational score is
`128749460/1000000`.

Preprocessing work is `(65536·1000 + 16777216)/1000 = 82313.216`, whose
base-2 log rounds upward to **16.328837**. It is included in `W`.

## Memory

Slots: `4T`. Record: 192 bytes. Fixed scratch, code, the 16 branch points and
the 4096 live walk states: `2^20` bytes.

```text
bytes = 4 T · 192 + 2^20
      = 157022823083577662705442953477403971682304
log2  = 136.850026158775…
```

Rounded upward, **136.850027**. Advice bytes: 0, so
`nonuniform_advice_log2_bytes = 0`. The branch table and the orbit table are
attack working storage, built inside `H` and `R`, not embedded constants.
`β` and `λ` are public curve constants and are excluded by the cost model.

## What is heuristic

The product bound is a theorem for a uniformly random function. The 16-branch
adding walk on canonical sextuples is not proved to be one. Teske's adding
walks and the canonicalisation that kills negation 2-cycles are the reason the
model is plausible; the toy is evidence only at orders `≤ 2003`. The hash
uniformity used in the `2^{−254}` abort bound is the same kind of statement.
Both are disclosed. The floor that is actually claimed is `0.39`, and the
budget already pays for a walk `5/4` longer than the `1/2` quantile of the
ideal model.

No quantum step, no uncharged table, no amortisation over a second target.

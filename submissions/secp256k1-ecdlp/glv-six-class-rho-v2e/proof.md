# Six-class rho, distinguished-point detection, unit-frame coefficients (v2e)

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
endomorphism multiplication, canonicalisation, the distinguished-point test and
the hash probe are other instruction work `R` and are divided by `C` once. A
reference operation is never also unpacked into its six multiplications.

Public curve constants are not advice. In particular
`β = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee`
satisfies `β²+β+1 ≡ 0 (mod p)`, and
`λ = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72`
satisfies `λ²+λ+1 ≡ 0 (mod r)`. Both checks are integer identities on the
protected `p` and `r`.

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
order-2 point, so `y ≠ p−y` and the sign is unique.

## One-multiplication canonicalisation

Because `β² ≡ −β − 1 (mod p)`, the third candidate `x`-coordinate is never a
second field multiplication:

```text
x2 = β·x            one 5×52 field multiply, 179 instructions
x3 = β²·x = p − x − x2     one 4-limb subtraction tree, no multiply
```

`x3` as written may be negative once; it is taken mod `p` by the same
conditional add that any 4-limb subtraction already carries. The minimum over
`{x, x2, x3}` is a three-way limbwise compare, and the unit index
`e ∈ {0,1,2}` falls out of the same comparison. The sign then compares `y`
with `p−y`. Earlier revisions charged `2×(179+16) = 390` for two full
multiplications; the identity `β²+β+1=0` removes one of them.

## Coefficients: the unit frame

The walk no longer stores the relation `(α, βs)` in the `Z[λ]` basis and does
not multiply two basis pairs by the canonicalisation unit every step. Each
walk carries a **unit frame**

```text
α_t = U_t · (α_0 + A_t)   (mod r)
β_t = U_t · (β_0 + B_t)   (mod r)
```

where `U_t` is the accumulated unit — a residue mod `r`, tracked as its index
in `{1, −1, λ, −λ, λ², −λ²}` — and `A_t, B_t` are plain residues mod `r`.
The six units form a group under composition, and each unit has an inverse
unit: the inverse of `λ` is `−λ²`, of `λ²` is `−λ`, and the negatives invert
as themselves up to sign.

Step algebra. After adding branch `S_j = [a_j]P + [b_j]Q` and applying
canonicalisation unit `u`, the true coefficients become `u·(α_t + a_j)` and
`u·(β_t + b_j)`. Substituting the frame,

```text
u·(α_t + a_j) = u·U_t·(α_0 + A_t + U_t^{-1}·a_j).
```

So the update, applied **before** `U` is composed, is

```text
A ← A + u^{-1}(U) · a_j   (mod r)
B ← B + u^{-1}(U) · b_j   (mod r)
U ← u · U                  index composition
```

where `u^{-1}(U)` is the inverse unit of the **current** `U`. The whole step
needs at most the six values `u^{-1}·a_j` and `u^{-1}·b_j`, one per possible
`U`: a precomputed public table `P[6][16]` of pairs `(a', b') = u^{-1}·(a_j, b_j)`
with `a', b' < r` costs `2·6·16` entries × a few scalar reductions each, built
inside the `2^16` preprocessing reserve (each entry is one constant multiply by
a unit — for units `±1` it is an identity or a negation; for `±λ, ±λ²` a
one-time scalar multiply by the public constant `λ`, `96` such multiplies in
total). The per-step update is then two residue adds mod `r`.

At a collision the relation is read off in `O(1)`:

```text
α = U·(α_0 + A) = m_U·(α_0 + A)  mod r
β = U·(β_0 + B) = m_U·(β_0 + B)  mod r
```

with `m_U` the unit multiplier indexed by `U`. Two scalar multiplications by a
table'd unit constant, once per distinguished-point hit. No walk replay, no
per-step `Z[λ]` pair multiplication, no checkpoint store. The toy experiment
checks `U·(α_0+A) + secret·U·(β_0+B) ≡ canonical value` on every step, so the
frame algebra itself is executable evidence.

## Walk

There are 4096 walks. Walk `w` starts at `[w+1]P + Q`, is replaced by its
canonical representative under unit `u_0`, and its frame is set to
`U = u_0, A = w+1, B = 1`, so that the frame evaluates to the canonical
coefficients from step 0. The 16 branch points are public,

```text
S_j = [a_j]P + [b_j]Q,   a_j = 2j+1,   b_j = 3j+5,   j = 0..15,
```

built by double-and-add. The branch index of a canonical point with integer
`x` is

```text
mixed = (x mod 2^64) * 0x9E3779B97F4A7C15   mod 2^64
j     = mixed >> 60.
```

One step, from an already canonical point `R` with frame `(U, A, B)`:

1. Read `j` from `R`. If `R` and `S_j` have the same `x`, or either point is
   at infinity, replace `j` by `j+1 mod 16`. If the second branch is also
   incomplete, abort this attempt (charged below).
2. Add `S_j` with one batched incomplete affine addition. The 4096
   denominators of a round are one Montgomery batch: one reference operation
   per walk per round.
3. Canonicalise the sum: compute `x2 = β·x`, `x3 = p − x − x2`, take the
   minimum of `{x, x2, x3}` for `e`, choose the sign from `y` vs `p−y`, and
   combine into the unit `u` applied.
4. Update the frame: `A += u^{-1}(U)·a_j`, `B += u^{-1}(U)·b_j`, `U ← u·U`.
5. Distinguished-point test: `R` is a distinguished point iff its canonical
   `x` has its low 20 bits zero, i.e. `x mod 2^20 = 0`. If not, the step is
   done. If yes, probe the distinguished-point table; see below.

## Distinguished-point table

Only canonical points whose low 20 `x` bits vanish are probed and stored. The
store rate is `2^{-20}` per walk-step; over the whole budget `T` the expected
number of store probes is `T·2^{-20}` and the number of table probes
(including hits that stop the walk) is at most `T·2^{-20} + T` in the worst
case but `≈ T·2^{-20}` in the charged expectation — the charged reserve below
covers `2^11` instructions per expected store plus the hit checks.

The table has `M = 4·⌈T·2^{-20}⌉` slots, load `1/4` under the same uniformity
heuristic as before. Each slot stores canonical `x` (32 bytes), the frame
`(U, A, B)` as `U`-index plus two residues (66 bytes), the walk id (2 bytes),
padded with a link word to 128 bytes. Each probe compares one `x` and follows
one link; a probe scans at most 64 slots and pays for all 64 whether or not
the chain is shorter; if an insertion would make length 64, the attempt
aborts (charged as before; the union bound is only stricter with `2^20` fewer
insertions).

On a hit with a different walk id, solve the linear relation
`x = (α' − α)·(β − β')^{-1} mod r` from the two stored frames. On a zero
denominator the collision is fruitless — the two walks entered the orbit with
identical `β` difference — count it and keep walking, exactly the v2d
treatment. On a hit with the same walk id the walk met its own earlier
distinguished point; that is a cycle on the same walk, skip the store and
keep walking.

Detection delay. After two walk prefixes first meet at some orbit element,
both continue on the **same** canonical sequence, because the step function
`canon(P + S_{j(canon x)})` is a function of the orbit, not of the
representative. They reach the next distinguished point together; the
expected wait is `2^20` steps and the probability of no distinguished point
in `s` further steps is `(1 − 2^{-20})^s`. Uniformly over the meeting time
`τ ∈ [0, T)`,

```text
P(no DP within T − τ) = E[(1 − 2^{-20})^{T−τ}]
```

bounded crudely by splitting at `τ = T − 2^26`:
`P(τ > T − 2^26) ≤ 2^26/T < 2^{−101}` plus
`(1 − 2^{-20})^{2^{26}} ≤ e^{−64} < 2^{−92}`. Total detection failure
`< 2^{−80}`, inside the `0.11` slack with the other abort terms and not a
second objective. The same argument covers two distinct walks that merge and
a walk that re-enters its own earlier orbit element.

## Birthday budget

For a random mapping on a set of size `N`, the walk `Y_{i+1}=f(Y_i)` satisfies

```text
P(no repeat through t steps) = prod_{i=1..t} (1 − i/N) ≤ exp(−t(t+1)/(2N)),
```

because the image of a fresh point is uniform on the whole codomain. The
inequality `t(t+1) ≥ 2 N ln 2` is therefore enough for collision probability
at least `1/2`. The series

```text
ln 2 < sum_{k=1..60} 1/(k 2^k) + 1/(61 · 2^60)
```

is an explicit rational upper bound. The least positive integer `t0` with
`t0(t0+1)` at least the ceiling of `2 N` times that rational is

```text
t0 = 163565440712060065318169743205629135679.
```

As in v2d, no `5/4` length factor is applied on top of the `1/2` quantile; the
charged budget is the least multiple of 4096 that is at least `t0`:

```text
T = ceil(t0 / 4096) * 4096
  = 163565440712060065318169743205629136896.
```

`T − t0 = 1217`, so `T > t0`. The submitted success floor stays `0.39`. The
absolute gap `0.11` is the allowance for the adding walk not being a
uniformly random function and for fruitless denominators, and — new in this
revision — for the distinguished-point detection delay bounded at `2^{−80}`
above. Failures inside this single budget are already paid for; there is no
restart.

Starting the 4096 walks on distinct points only adds samples, so the serial
product bound remains an upper bound on the miss probability.

Chain-length abort: a given slot has mean `μ = 1/4` under an independent
uniform hash into `M = 4⌈T·2^{-20}⌉` slots. The same Chernoff form
`(eμ/64)^64 < (1/64)^64 = 2^{−384}` applies per slot, and the union is over
`M < 2^{110}` slots: below `2^{−274}`. Incomplete additions: at most two bad
partners out of `r` per addition, two branches tried, so the expected count
is `< 8T/r < 2^{−120}`; the attempt aborts if that count exceeds `2^{20}`,
whose Chernoff distance is far below the union allowance. Neither abort moves
`0.50` to `0.39`.

## Toy experiment

The experiment runs the same walk on `Z/n` with `n` prime, `n ≡ 1 (mod 6)`,
`n ≤ 2003`: the same orbit multipliers `{±1, ±λ, ±λ²}`, the same branch index,
the unit frame `α_t = U_t·A_t`, `β_t = U_t·B_t` checked on every step against
the expanded scalar, and the distinguished-point store. The toy's
distinguished modulus is `2` (canonical `cv ≡ 0 mod 2`), a stated scale
parameter: with `t0 = 15` at `order = 1009` a modulus of `2^20` would observe
almost no stores, so the toy exercises the mechanism at a rate adapted to its
size, exactly as it uses 8 walks rather than 4096.

Author execution of the declared input `order=1009` (`N_orb=169`, `t0=15`, 8
walks) recovered `400/1009` secrets at the charged budget `t0` and
`866/1009` at `2·t0`, recorded `1` and `2` fruitless relations, and reported
zero frame-drift failures over `6539` and `9255` distinguished-point stores.
The `t0` figure sits just above the submitted floor `0.39`: at this size the
expected detection wait (2 steps) is a visible fraction of a 15-step budget,
while at production scale the wait `2^20` is negligible against `T ≈ 2^{127}`
— that is exactly what the `2^{−80}` bound above quantifies, and the `2·t0`
sweep shows the mechanism converging once the delay amortises. The toy stores
canonical values in a dictionary rather than a 64-probe chain; the chain cap
is an instruction bound, not part of the relation. It is evidence only at
orders `≤ 2003`.

## Instruction ledger, per step

| Item | Instructions |
| --- | ---: |
| `x2 = β·x`, one 5×52 field multiply plus call | 195 |
| `x3 = p − x − x2`, two 4-limb subtracts with reduction | 24 |
| Three-way minimum of `{x, x2, x3}` with unit index | 96 |
| Sign selection, including `p−y` | 64 |
| Incomplete-addition test on both branch attempts | 48 |
| Branch index, one 64-bit multiply and a shift | 32 |
| Frame update: `P[U][j]` lookup, two residue adds mod `r`, unit compose | 80 |
| Distinguished-point test, `x mod 2^20`, branch | 16 |
| Round loop, pointers, batch index | 56 |
| Total `R_step` | 611 |

The three-way compare row includes the selects, the sign row includes `p−y`,
and the frame row's `80` covers the `P[6][16]` table read, two 4-limb adds
with conditional reduction mod `r`, and the unit-index composition. Every row
is a worst-case charge; no row is amortised. The distinguished-point store and
solve work is charged once in `R_dp` below, not per step.

The addition itself is one reference operation and is not in this table.
`H` counts every such addition, plus reserves that are not hidden inside `T`:

| Piece | Reference operations |
| --- | ---: |
| Walk additions `T` | 163565440712060065318169743205629136896 |
| Preprocessing reserve `2^16` | 65536 |
| Verification reserve `2^12` | 4096 |
| Incomplete-addition reserve `2^20` | 1048576 |
| Total `H` | 163565440712060065318169743205630255104 |

Preprocessing really uses at most `4096` adds to build `[1]P..[4096]P`,
`4096` adds of `Q`, `16·2·512` double-and-add steps for the branch scalars,
and the `96` unit-table scalar multiplications; the `2^16` line covers all
of it plus the initial canonicalisations. Verification is one double-and-add
for `[x]P` plus encoding.

Distinguished-point table work is charged once, not per step: the expected
number of stores is `⌈T·2^{-20}⌉` and each store or hit probe costs at most
`2^11` instructions (64-slot chain scan at 12 instructions per probe, the
record write, the two unit multiplications and the solve on a hit):

```text
R_dp = 2^119  ≥  2^11 · ⌈T·2^{-20}⌉
     (the product is 319463751390742327525184144079323136 < 2^119)
```

Other instructions outside the per-step loop:

```text
R_pre    = 2^24 = 16777216
R_verify = 2^22 = 4194304
R        = T · 611 + R_dp + R_pre + R_verify
         = 99939148889066592367338165002169563787264
```

`R_pre` covers bit scans, the `P[6][16]` unit table and addressing for the
preprocessing additions. `R_verify` covers the final solve's scalar
multiplications, the binary extended gcd for the denominator inverse, and
output encoding.

```text
W = H + R/1000
  = 32938073700140832210688488525974977361408 / 125.
```

`log2(W) = 127.631095695996…`. Rounding upward at the sixth decimal place
gives the submitted total **127.631096**. The rational score is
`127631096/1000000`.

Preprocessing work is `(65536·1000 + 16777216)/1000 = 82313.216`, whose
base-2 log rounds upward to **16.328837**. It is included in `W`.

## Memory

Distinguished-point table: `M = 4·⌈T·2^{-20}⌉ = 623952639435043608447625281404928`
slots at 128 bytes. Fixed scratch, code, the 16 branch points, the 96-entry
unit table and the 4096 live walk states: `2^20` bytes. No per-step store of
canonical values exists any more; the `4T` orbit table of earlier revisions
is removed.

```text
bytes = 623952639435043608447625281404928 · 128 + 2^20
      = 79865937847685581881296036020879360
log2  = 115.943135563166…
```

Rounded upward, **115.943136**. Advice bytes: 0, so
`nonuniform_advice_log2_bytes = 0`. The branch table, the unit table and the
distinguished-point table are attack working storage, built inside `H` and
`R`, not embedded constants. `β` and `λ` are public curve constants and are
excluded by the cost model.

## What is heuristic

The product bound is a theorem for a uniformly random function. The
16-branch adding walk on canonical sextuples is not proved to be one. Teske's
adding walks and the canonicalisation that kills negation 2-cycles are the
reason the model is plausible; the toy is evidence only at orders `≤ 2003`.
Hash uniformity for the `2^{−274}` chain-overflow bound and `x mod 2^20` for
the distinguished-point rate are the same kind of statement: on canonical
residues mod `r` the low bits are a regular distinguisher, and a bias `ε` in
the distinguished rate changes only the detection-delay term, which is
already `2^{−80}`. The relation is not heuristic: the frame identity is
exact algebra checked by the toy on every step, and a distinguished-point
collision yields a linear equation in `x` whose failure mode (zero
denominator) is disclosed and counted.

The floor that is actually claimed is `0.39`. The budget pays for a walk of
length `T`, the least multiple of 4096 that is at least the `1/2` quantile of
the ideal model, and not for a longer pad.

No quantum step, no uncharged table, no amortisation over a second target.

# Distinguished-point Teske adding walk

Status: ready for exploratory review. This file does not emit a score.
The input is a point `Q` on the protected secp256k1 profile. The output is
the scalar `x` with `Q = [x]P`. The cost model is `classical-group-work-v2`.

SHA-256 of the raw bytes of `tracks/secp256k1-ecdlp/target-profile.json` is
`18c3338de676f1909715814325bffdb1513f1621d678b942e8a7a9563b8fd086`.
SHA-256 of the raw bytes of `tracks/classical-group-work-v2.json` is
`e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294`.

The reference operation is one batched incomplete affine addition, normalised
at `C = 1000`. It already contains the six-multiplication convention: three
multiplications for the incomplete formula and three for batch amortisation,
with `S = M` and `A = 0` inside that formula only. A reference operation is
never unpacked into those six multiplications again. Scalar adds, the branch
mix, distinguished-point tests, table probes and output encoding are other
instruction work `R` and are divided by `C` once.

The walk is on the full prime-order subgroup. It does not canonicalise under
negation, does not apply the GLV endomorphism, and does not pass to the
quotient by an order-6 automorphism group. The collision set has size `r`.

## Branch table

Sixteen fixed multipliers, each below `2^64`. They are algorithm constants,
not searched discrete logarithms. `β` and `λ` are not used.

```text
A = [9315109294132905716, 9476356238027737019, 17445819973253803743,
     3342237045799953233, 2812045804143842406, 9122178642407602516,
     13926598530136396895, 17056811102254473403, 11323810722078648909,
     15698505963632274993, 9578086386675406826, 4470105836484457736,
     17565181521087248143, 3097285154792831098, 12002478175040212671,
     12486765383173650733]
B = [8546391839045156921, 14157483109116380921, 7783733875923022258,
     5661477030581511630, 7704659425681353120, 15728860691533551284,
     13562443690413809913, 17757487276652662861, 2212146067940045771,
     14634568361653561203, 9217311395089492335, 9124695148773980898,
     2502434592559100122, 7616228408119649253, 3985732427038351031,
     2215723844435459483]
```

`S_j = [A_j]P + [B_j]Q`. Each scalar is built by double-and-add, at most 512
incomplete additions, then one addition attaches `Q`. Sixteen branches
therefore cost `16·512 = 8192` reference operations. This is charged in full
below.

The branch index of a finite point with integer coordinate `x` uses only the
low 64 bits:

```text
z = (x mod 2^64) * 0x9E3779B97F4A7C15   mod 2^64
z = z xor (z >> 30)
z = z * 0xBF58476D1CE4E5B9                 mod 2^64
z = z xor (z >> 27)
j = z >> 60
```

Multiplication by an odd constant is a bijection on 64-bit words, and the
two xors are a fixed mix. On a uniform coordinate the index is treated as
uniform on `{0,…,15}`. That uniformity is part of the random-mapping
heuristic below, not a theorem.

## Walk

There are 4096 walks. Walk `w` starts at `[w+1]P + Q`. Building the starts
by adding `P` along the chain costs 4096 reference operations.

One round applies one step to every walk:

1. Read `j` from the current point. If that point and `S_j` have the same
   `x`, or either is infinity, replace `j` by `j+1 mod 16`. A second
   incomplete pair aborts the step into the reserve below.
2. Add `S_j` by the incomplete affine formula. The 4096 denominators of the
   round are one Montgomery batch. Each of those additions is one reference
   operation. The batch width is 4096, so the single Fermat inversion inside
   the batch is far below the three multiplications per element already
   priced into `C`. Charging one reference operation per addition does not
   under-count the inversion.
3. Add `A_j` to the `P`-coefficient and `B_j` to the `Q`-coefficient,
   modulo `r`. Both are 4-limb modular adds.
4. If the low 16 bits of the new `x` are zero, the point is distinguished
   and is inserted into the collision table.

Coefficients stay unreduced to a single 256-bit product until a collision
is solved. A round is charged even if a walk has already collided: the
budget is a hard cap, not an expectation.

## Distinguished points and the relation

A stored record is `(x, y, α, β)` with `R = [α]P + [β]Q`. Two records with
the same point and `β₁ ≢ β₂ (mod r)` give

```text
x ≡ (α₂ − α₁) · (β₁ − β₂)^{−1}  (mod r).
```

The attack checks `[x]P = Q` before returning. A zero denominator is a
failed detection inside the same budget, not a second run and not a wrong
accepted scalar. Infinity is stored as an explicit tag and obeys the same
relation.

The table is open-addressed. Let `D` be the number of distinguished
inserts. The production cap is `D_max` below. The slot count is the power
of two `M = 2^116`, which is at least `4 D_max`, so the load is at most
`1/4`. Probe length is `k = 80`. Slot `i` of a key is

```text
(h + i · 0x9E3779B97F4A7C15) mod M
```

with `h` the same 64-bit mix. Because `M` is a power of two the reduction
is a bitmask, and that bitmask sits inside the 12 instructions charged per
probe. An insert that finds 80 occupied slots aborts the attempt. That
abort is bounded under the uniform-slot heuristic in the next section.

## Success

Let `N = r`. Write `ln 2` for the series `Σ_{k≥1} 1/(k 2^k)`. The partial
sum through `k = 60` plus the remainder majorant `1/(61·2^60)` is an upper
bound `L > ln 2`. Set

```text
need = ceil(2 r L)
     = 160521920371982617724519308230885847593108909541466899510640560181104765873885
t0   = min { t ≥ 2 : t(t−1) ≥ need }
     = 400651869298001176472378576439342816068.
```

For uniform independent samples, `t0(t0−1)/(2r) ≥ L > ln 2`, so the
probability of a repeated value by step `t0` is at least `1/2`. The same
quadratic is the usual approximation for the first repeat of a random
mapping. A 16-adding walk is not proved to be a random mapping. Teske's
adding walks are the experimental support for using that approximation as
a heuristic. The toy in `experiment.py` runs this walk on `Z/n` with the
same multipliers, the same mix and the same budget formula.

The charged length is longer by `3/2`, then rounded up to a multiple of
4096 so every round is a full inversion batch:

```text
pad = ceil(3 t0 / 2)
    = 600977803947001764708567864659014224102
T   = ceil(pad / 4096) · 4096
    = 600977803947001764708567864659014225920.
```

Under the random-mapping model the collision probability at `t0` is already
at least `1/2`, and `T > (3/2) t0`. The submitted success floor is `0.39`.
The gap between the model probability at this longer budget and `0.39` is
the allowance for the adding walk not being a uniform random function, for
a vanishing denominator, and for the probe abort. It is not a second
objective. There is no restart.

Distinguished points use modulus `2^16`. A random mapping's cycle, once the
walk has run past `t0`, has length on the scale of `sqrt(r)`, far above
`2^16`. The probability such a cycle contains no distinguished point is
negligible beside the `0.39` floor. After the first collision the walk is
continued to the hard cap `T`, and the mean gap between distinguished
points on a uniform coordinate is `2^16` steps, which is inside the `3/2`
pad. This detection delay is part of the same heuristic, not an extra
unpriced search.

The expected number of distinguished hits on `T` steps is `T/2^16`. A
Chernoff bound at twice that mean is far below `2^{−40}` if the low 16
bits behave like independent fair bits. The table is sized for the doubled
count, so an over-production of distinguished points aborts only on that
already tiny event:

```text
D_exp = ceil(T / 2^16)
      = 9170193541671779857003293833297947
D_max = 2 D_exp
      = 18340387083343559714006587666595894
M     = 2^116
      = 83076749736557242056487941267521536.
```

`4 D_max < M`, so the load is strictly below `1/4`. `D_max < 2^114`.
Under independent slot occupancy at load at most `1/4`,
one insert survives 80 consecutive occupied probes with probability at most
`4^{−80} = 2^{−160}`. The union over fewer than `2^114` inserts is less
than `2^{−46}`. Dependence among linear-probe slots is a stated limitation
of this estimate; it is not given a second full run. Incomplete additions
are a separate reserve of `2^20` reference operations. Each step has at
most two bad partners out of `r`, so the expected number of incomplete
events on the whole budget is below `2^{−100}`. The reserve covers them
without another launch.

The declared toy (`order = 1009`, 8 walks, budget 64, every point stored)
recovered 926 of 1009 secrets. An exhaustive check at order 5003 with the
same code recovered 4073 of 5003. Both sit above `0.39`. Neither check is
a production complexity proof. The production floor remains the supported
heuristic bound, not the toy rate.

## Instruction ledger

Every one of the `T` steps pays one reference addition. On top of that, a
step is charged 80 scalar instructions, whether or not it is distinguished:

| Piece | Instructions |
| --- | ---: |
| 64-bit mix to a branch index | 16 |
| Load `A_j`, `B_j` and form the addresses | 12 |
| Two 4-limb adds modulo `r` | 32 |
| Test the low 16 bits of `x` | 4 |
| Walk pointer, loop, reserved control | 16 |
| Total per step | 80 |

A distinguished insert is additionally charged the full probe cap, not the
average probe: 80 probes at 12 instructions (address, load, four limb
compares, branch) plus 16 instructions to store the record, `976`
instructions. At most `D_max` inserts are charged. Empty-slot initialisation is 32 instructions per slot: 24 stores cover a
192-byte record and 8 instructions cover addressing.

```text
R_step  = 80 T
R_dp    = 976 D_max
R_init  = 32 M
        = 2658455991569831745807614120560689152
R_tail  = 2^24
```

`R_tail` covers decoding, the final modular inverse of the denominator by
extended gcd (256 iterations of 4-limb division), both verification scalar
multiplications' address work, and output encoding. The group operations
inside those two verifications are not in `R_tail`. They are 1024 reference
operations (two double-and-add chains, 512 each).

```text
H = T + 2^20 + 1024 + 8192 + 4096
  = 600977803947001764708567864659015287808
R = 80 T + 976 D_max + R_init + 2^24
  = 48098782989545054322712107216404313132512
W = H + R/1000
  = 649076586936546819031279971875419600940512 / 1000.
```

`log2(W) = 128.931656322…`. Rounding upward at the sixth decimal place
gives **128.931657**. The rational score is `128931657/1000000`.
`2^128.931657` is above `W`.

Preprocessing is the branch build, the slot initialisation, and `2^20`
instructions of constant setup. It does not depend on the secret scalar.
It is included in `W`:

```text
W_pre = (8192 · 1000 + R_init + 2^20) / 1000
      = 2658455991569831745807614120569929728 / 1000
log2  = 111.034215…
```

Rounded upward, **111.034216**.

## Memory

Each of the `M` slots is a 192-byte record (`x`, `y`, `α`, `β`, tag and
padding), allocated up front, including the empty slots. Scratch, code,
the 16 branch points and the 4096 live states sit in an extra `2^20` bytes.

```text
bytes = 192 M + 2^20
      = 15950735949418990474845684723365183488
log2  = 123.584962…
```

Rounded upward, **123.584963**. Advice bytes are 0, so
`nonuniform_advice_log2_bytes = 0`. The branch multipliers are fixed
constants of the algorithm. The collision table is working storage built
inside `R_init` and the walk, not an embedded answer. Public curve
parameters are excluded by the cost model.

## What is heuristic

The product bound is a theorem for uniform samples and the usual
approximation for a uniform random mapping. The 16-branch adding walk on
secp256k1 points is not proved to be one. Teske's adding-walk experiments
are the support for importing the approximation, and the toy is the support
for importing it for these particular multipliers and this mix. The
distinguished-point detection argument and the `2^{−46}` probe union bound
assume that low bits and probe slots behave like the uniform model. The
submitted floor `0.39` is what those approximations are asked to carry.
No claim is made that the true success rate equals the toy rate, or that
this is a practical computation.

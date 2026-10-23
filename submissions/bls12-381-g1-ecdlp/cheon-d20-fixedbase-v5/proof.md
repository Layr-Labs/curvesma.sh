# bls12-381-g1-ecdlp: Cheon d=2^20 with fixed-base windowed BSGS

Status: ready for exploratory review. This is an analytic attack and accounting
claim, not a full-size execution, a novelty claim or a practical break.

## Changes in this revision

The previous review accepted the Cheon preconditions, the two-phase recovery
identities, the toy recoveries and the ledger arithmetic, and named three
defects. Each is answered here by name.

| raised item | where | what changed |
| --- | --- | --- |
| **the stopping event was false**: a *per-operation* 12-link cap was combined with a Markov bound on the *total* link count, and one 13-link operation aborts a run whose total is far below `12*N` | `## Matching`, `## Matching cost and success probability` | the per-operation cap is **withdrawn**. A single **global counter** accumulates every link walked by the attempt and aborts when the counter would exceed `12*N`. The abort event and the Markov event are now the same event, and no per-operation bound is claimed anywhere |
| **both** node references **and bucket heads** need more than 64 bits for an `m > 2^117`-node arena | `## Matching`, P2, `## Memory` | nodes are **128 bytes** with a 128-bit next reference and heads are **16-byte references**; the load is priced as two `LDR`, the itemisation becomes `203 + 13` per extra link (worst case 359 of the 560 charge), and the memory bound is recomputed and **raised to `2^125`** to cover the 16-byte heads |
| the tag was taken from the low 128 bits, where `a` odd makes equal keys mod `2^128` collide for **every** coin, so the `2^-247` joint factor was unavailable; and full-point comparisons on tag matches were not separately bounded | `## Matching`, `## Matching cost and success probability` | the tag is the 128-bit window **at bit 137**, immediately below the bucket index, so `(h,t)` is the top 247 bits and `2/2^247` holds; and a second **global counter** caps false full-point comparisons at `2^20`, charged `2^20 * 24` and bounded by Markov at `2^-143` |
| the toy used a Python dictionary, not the declared chained table | `experiment.py` | the experiment now implements the declared separate-chaining table with the multiply-shift hash and the global link counter, at the production bucket/node ratio, and reports its link count, chain length and abort flag |

## Target and notation

The protected profile fixes the BLS12-381 G1 group of prime order

```text
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
```

and a generator `P0`. At the scoring length `n = 2^20 = 1048576` the challenge
supplies `A_i = [s^i mod r] P0` for `i = 1..n` and `T = A_1 = [s]P0`, where `s` is
uniform in `[1, r-1]`. Write `G = <P0>`, so every point is `[x]P0` for a unique
`x mod r`, and write `F_r^*` for the multiplicative scalar field.

The attack uses only the single supplied power `A_d` with `d = 2^20`. Since
`2^20 | r-1`, this is the largest admissible `d <= n`; no larger divisor of `r-1`
is at most `n`. The exponent `2^20` is fixed, so the same work function applies at
`n = 2^20` and at every larger `n`; the bound is not claimed for smaller `n`.

Set

```text
M1 = (r-1)/d = 50006747412801924209068050869165387952509453297164571592906626415194112
m  = ceil(sqrt(M1)) = 223621884914696853578958571250384896
g1 = ceil(M1/m)     = 223621884914696853578958571250384896
N  = m + g1         = 447243769829393707157917142500769792
```

One checks `m^2 - M1 = m*g1 - M1 = 218380746987008646073201729736704 > 0` and
`(m-1)^2 < M1`, so `m` and `g1` are the exact BSGS parameters. Also `N =
2*sqrt(M1)` up to the same small ceiling, and `log2 N = 118.42854470652...`.

## The two Cheon phases

Because `d | r-1`, let `zeta` be an element of order `r-1` in `F_r^*` (a primitive
root modulo `r`), and write `s = zeta^a` with `0 <= a < r-1`. Define

```text
zeta_d = zeta^d        (order M1 = (r-1)/d)
omega  = zeta^M1       (order d)
```

**Phase 1.** Since `zeta_d` has order `M1`, the point `A_d = [s^d]P0 =
[zeta_d^a]P0` determines `a mod M1`. Write `k0 = a mod M1`. The standard
representation is `k0 = u + m*v` with `u = k0 mod m` and `v = floor(k0/m)`, so
`0 <= u < m` and `0 <= v < g1`; the pair exists. Because `m*g1 > M1`, the
scanned rectangle can contain one more pair with `u + m*v = k0 + M1` when
`k0 < m*g1 - M1 = 218380746987008646073201729736704`; reducing modulo `M1`
returns the same `k0`. The identity used is

```text
[zeta_d^{-u}] A_d = [zeta_d^{m*v}] P0.
```

The left side is computed for all `u` in `[0,m)` (baby steps, base `A_d`); the
right side for all `v` in `[0,g1)` (giant steps, base `P0`). A match returns
`k0 = (u + m*v) mod M1`. This is equation (1) of Cheon's Lemma 1 adapted to the
supplied `A_d`.

**Phase 2.** The value `a` is `k0 + M1*k1` for a unique `k1` in `[0,d)` because
`0 <= a < M1*d = r-1`. Define

```text
T' = [zeta^{-k0}] T = [omega^{k1}] P0.
```

Write `k1 = u' + m2*v'` with `m2 = ceil(sqrt(d)) = 1024`, `0 <= u' < m2`, and
`0 <= v' < g2 = ceil(d/m2) = 1024`. The tested identity is

```text
[omega^{-u'}] T' = [omega^{m2*v'}] P0.
```

Baby steps over `u'` use the point `T'` as a variable base (one variable-base
scalar multiplication each, charged at 384 additions); giant steps over `v'`
reuse the `P0` tables. A match returns `k1 = u' + m2*v'`. The output is

```text
s = zeta^{k0 + M1*k1} mod r,
```

followed by the mandatory attack-side check `[s]P0 = T` (one more fixed-base
multiplication with the `P0` tables). Returning a point, a relation or a next
power is excluded by this check.

The only primitive-root machinery is: factor `r-1`, test candidate `zeta` in
`F_r^*` until `zeta^{(r-1)/q} != 1` for every prime `q | r-1`. The factorization
is found by trial division up to `2^28` (the largest prime factor of `r-1` is
`254760293 < 2^28`), at a charged cost far below `2^40` instructions; the test is
charged in the same budget. The smallest primitive root modulo `r` is `7`
(candidates `2,3,4,5,6,7`), so the search is only a few tens of thousands of
255-bit modular multiplications; trial division to `2^28` and sequentially
reading the input sequence up to `A_n` also cost far below the `2^40` budget.
The `2^40` setup charge therefore covers the input access with a wide margin. No
factorization or primitive root is embedded as uncharged advice. For reference,
the factorization checked by the experiment is

```text
r-1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2
      * 2508409 * 2529403 * 52437899 * 254760293^2.
```

## Fixed-base windowed multiplication

All scalars in the two phases are known field elements: `zeta_d^{-u}` for the
baby steps, `zeta_d^{m*v}` for the phase-1 giants, and `omega^{...}` for phase 2. Only
two point bases occur in phase 1: `A_d` and `P0`. Because `r < 2^255`, three
85-bit windows cover every scalar:

```text
x = x0 + 2^85 x1 + 2^170 x2,   0 <= xj < 2^85.
```

For each base `B` in `{A_d, P0}` and each `j` in `{0,1,2}` build the table

```text
T[B,j][c] = [c * 2^(85*j)] B,   c = 0..2^85-1.
```

Then `[x]B = T[B,2][x2] + T[B,1][x1] + T[B,0][x0]`: three lookups and exactly two
point additions (the first table entry initializes the accumulator). The tables
are built by repeated addition from `B`, `[2^85]B`, and `[2^170]B`; the latter two
cost 85 and 170 doublings. Table construction costs

```text
6*(2^85 - 1) additions + 2*(85+170) doublings
```

for the two bases, and is charged conservatively at three reference operations
per entry (a projective addition plus batched affine conversion), giving
`H_pre = 3*6*(2^85-1) + 510 = 696341272098026404630757868`.

The main loop therefore performs

```text
H_main = 2*(m + g1) = 2*N = 894487539658787414315834285001539584
```

reference additions: two per scalar multiplication, for `N = 447243769829393707157917142500769792`
scalar multiplications. This is the dominant term, `log2 H_main = 119.42854470652...`.

The `6M` convention is the track's batched incomplete affine addition. In the
main loop the two additions per scalar multiplication are processed in blocks
(block size `2^18`): the first addition for every element of a block is batched,
then the second, so the `3M` batch-amortisation part of the `6M` convention is
realized with a batched inversion, and a block's inverse and output buffers stay
below `2^26` bytes in total. The table-construction additions are charged
at three reference operations each, which covers the less convenient batched
conversion there.

Identity table entries (`c=0`) are skipped by a branch whose instruction cost is
included in the lookup and digit charges. For every nonzero scalar the two
additions need no exceptional correction, and this is proved directly rather
than bounded statistically. The base `B` is `A_d` or `P0`, both of order `r`,
so `[u]B` and `[v]B` have equal x-coordinates exactly when `u = v` or
`u + v = r`. For the first addition, `T[B,2][x2] = [2^170 x2]B` and
`T[B,1][x1] = [2^85 x1]B` can have equal x-coordinates only if
`2^85 x2 = +/- x1`; the two sides differ by less than `r` because
`2^170 + 2^85 < r`, so the congruence is an integer equality, and `0 <= x1`
then forces `x2 = x1 = 0` (a skipped addition). For the second addition let
`S = [x - x0]B` be the first sum and `T[B,0][x0] = [x0]B` the third entry.
Equal x-coordinates give `x = 2x0` (equivalently `x - x0 = x0`), `x = 0` or
`x = r` in the integers; for `0 < x < r` only `x = 2x0` survives, and then
`x0 = x - x0 = 2^85*(2^85*x2 + x1)`, impossible since `0 <= x0 < 2^85`.
Hence no incomplete-addition exception occurs in either Cheon phase or in the
final verification, and the digit-zero skip branches are already charged in the
lookup and digit terms. `H_exc = 1024` remains as rounding and branch slack; the
`2*(85+170)` base doublings and the entry `c = 2` of each table are already
charged in `H_pre`.

## Matching

The baby points `B_u = [zeta_d^{-u}]A_d` are distinct: their exponents
`zeta_d^{-u}` are distinct for `u` in `[0,m)`. They are stored in a
**separate-chaining** table, so no probing law and no load-factor heuristic is
needed anywhere: `M = 2^119` bucket heads of **16 bytes** (`2^123` bytes) plus
at most `m` nodes of **128 bytes** (96-byte affine point, 128-bit tag, **128-bit
next reference**). The node arena is `128*m = 2^124.43` bytes, which is exactly
why the stored next reference is 128-bit and not 64-bit: a 64-bit field cannot
address a position in an arena larger than `2^64` bytes, so the reference is a
2-word quantity, its load is priced as two `LDR` in P2, and the 128-byte node is
what the memory bound in `## Memory, preprocessing and advice` charges. **The
bucket heads are 16-byte references for the same reason**: a head must name a
node in the same `2^124.43`-byte arena, so an 8-byte head cannot address it
either, and the head array is charged at 16 bytes per bucket. The load
`m/2^119 = 0.336...` is stated only so a reader can see the table is sparse; no
claim below uses it.

*Concrete hash.* Write the affine point's x-coordinate as a 384-bit
little-endian integer `x` (six 64-bit limbs; the top limb is partial because
`x < 2^381`). Draw `a` uniformly among the odd integers in `[1, 2^384)` and `b`
uniformly in `[0, 2^384)` once per attempt from the attempt's fresh coins, and
set

```text
p    = (a*x + b) mod 2^384
h(x) = p >> 265                       (119-bit bucket index, so h < 2^119)
t(x) = (p >> 137) AND (2^128 - 1)     (128-bit tag stored in the node)
```

The tag is the 128 bits **immediately below** the bucket index, so `(h,t)` is
jointly the **top 247 bits of `p`**, and the multiply-shift bound is available at
width 247: `Pr[h(x)=h(y) and t(x)=t(y)] <= 2/2^247`.

The tag is deliberately *not* taken from the low bits of `p`. The previous
revision used `t = p AND (2^128-1)`, and that is unsound: because `a` is odd,
`x = y (mod 2^128)` implies `a*x = a*y (mod 2^128)` for **every** hash coin, so
the low 128 bits are always equal on such a pair and no 128-bit random factor is
available. With the tag taken from bits 137..264 no such deterministic
congruence exists, and the 247-bit bound above is the honest one.

This is the multiply-shift family (Dietzfelbinger-Hagerup-Katajainen-Penttonen).
For odd `a` it is 2-universal up to the factor 2: for `x != y`,
`Pr[h(x) = h(y)] <= 2/2^119`, and applied to the top 247 bits,
`Pr[h(x)=h(y) and t(x)=t(y)] <= 2/2^247`. Only pairwise independence is
ever used, once, in `## Matching cost and success probability`. No geometric
probing law, no independence between operations and no uniform-hashing
assumption appears anywhere in this candidate.

*Insertion and lookup.* An operation computes `(h(x), t(x))`, loads the bucket
head at `heads + 8*h(x)` (128-bit address arithmetic on a 64-bit ISA, priced in
P2), and walks the chain. Each link compares the stored 128-bit tag first; a
tag match triggers the full 96-byte point compare, which is exact (no x-only
and no sign-ambiguous match). On a miss the operation follows the 128-bit next
reference. Insertion pushes a new node at the head of its bucket. There is **no
per-operation link cap and no chain is ever truncated**: a single global counter
accumulates every link walked by the whole attempt and returns FAILURE when the
counter would exceed `12*N`. The matching work of an attempt is therefore
bounded by construction (P2) without any per-operation statement.

Each giant step `G_v = [zeta_d^{m*v}]P0` looks up its full affine point. On a
hit, the attack computes `k0 = (u + m*v) mod M1`, runs phase 2, and accepts only
if `[s]P0 = T`. Every full-point match satisfies `u + m*v ≡ k0 (mod M1)` and
therefore returns the correct `k0`; the possible second rectangle pair (only for
`k0 < 218380746987008646073201729736704`) returns the same `k0`. Hence every
accepted hit is correct, and collisions only cost links and can never create a
false full-point match. A miss or a failed check continues the fixed scan. The
true match is present, so the scan returns the correct `s` whenever no cap is
exceeded; the cap is exceeded with probability at most 0.125392, so
`success >= 0.874608 >= 0.85` (proved in `## Matching cost and success
probability`). The phase-2 table is the same structure with 1024 items; the same
cap applies and its cost is far below the `2^22` phase-2 support line.

## Itemised instruction accounting for the non-addition terms

`R` is charged in scalar AArch64 retirement counts; loads, stores, addressing,
control, hashing and memory transfers are included. The model used here: 64-bit
little-endian limbs, an `Fp` element = 6 limbs = 48 bytes, an affine G1 point
`(x, y)` = 96 bytes = 12 limbs, one chained node = 128 bytes = 16 limbs
(96-byte affine point + 128-bit tag + 128-bit next reference), one bucket head = 8
bytes, one comb entry = the same 96-byte affine point. Each table below prices **an explicit instruction sequence**, and the
charged figure is the itemised sum times the stated factor. Nothing here is a
timing measurement, and nothing is asserted without its itemisation.

### P1 — one comb lookup (3 per fixed-base scalar multiplication), charged 200

The tables hold 96-byte affine entries. An 85-bit digit is a **two-word**
quantity (bits 0-63 in one register, bits 64-84 in a second) and its byte offset
`digit*96` needs 128-bit arithmetic, so both are priced in 2-word lanes; nothing
in the table below is treated as a single-word quantity unless it is at most 64
bits wide.

| step | instructions |
| --- | ---: |
| digit extract, two words (`LSR`+`AND` low, `LSR`+`AND` high) | 4 |
| sign test and conditional-negate mask, two words (`TBZ`, mask builds) | 5 |
| byte offset `digit*96` in 2-word lanes (`LSL`+`ADD` low, `ADC` + shifted high) | 6 |
| base + offset, 128-bit address add (`ADD`, `ADC`) | 2 |
| 12 aligned 8-byte loads (`LDR`) | 12 |
| conditional negation, 12 limbs (`EOR`+`SUB`, 1.5 avg) | 18 |
| install into accumulator registers (12 x `MOV`) | 12 |
| **digit-zero exceptional case** (identity entry) test and branch | 3 |
| window loop control, amortised over 3 digits | 2 |
| **itemised sum** | **64** |
| **charged** | **200** (3.1x) |

The `c = 0` entry is the identity and is not an affine point; the branch is not
optional and is priced above, and the no-incomplete-addition proof of the
previous section shows no other exceptional path occurs. The two-word digit and
offset rows answer the review finding on the earlier draft, which priced both as
single-word operations; the charge is unchanged and still covers the itemised
sum by 3.1x.

### P2 — one hash-table operation (1 per fixed-base scalar multiplication), charged 560

`M = 2^119` bucket heads of 16 bytes plus `m` nodes of 128 bytes. The machine
model is a 64-bit ISA whose tables occupy a `2^124.6`-byte address space, so
**addresses and indices are 128-bit quantities held in register pairs and every
address computation is priced as a 2-word operation**. The bucket index `h(x)` is
119 bits wide and is never treated as a single-word quantity.

| step | instructions |
| --- | ---: |
| 384x384 -> 768-bit product for `a*x`: 36 `UMULH`+`MADD` limb products with `ADCS` carry chains | 108 |
| add `b` (6 words + carries), extract `p >> 265` (`LSR`/`ORR`, 2 words) and the 128-bit window at bit 137 (`EXTR`/`LSR`/`ORR`, 2 words) | 24 |
| bucket byte address `heads + 16*h(x)` (119-bit index) in 2-word lanes (`LSL`, `ADD`, `ADC`, high-word build) | 7 |
| load 16-byte bucket head (2 x `LDR` at a 128-bit address) | 3 |
| load the node's 128-bit tag (2 x `LDR`), compare (`CMP`, `CCMP`), branch | 5 |
| 128-bit next-reference load and chain advance (2 x `LDR`, 2-word update) | 5 |
| full 96-byte point compare on a tag match: 12 `EOR` + 11 `ORR` + `CBZ` | 24 |
| push node on insert (16 x `STR` for 128 bytes, `STR` of the head) | 17 |
| probe-loop control and cap counter, amortised | 10 |
| **itemised sum** | **203** |
| each additional link (tag load/compare/branch, 128-bit next reference, global counter) | 13 |
| **12-link allowance: 203 + 12 x 13** | **359** |
| **charged** | **560** (1.56x) |

The compare row is the accumulated-equality sequence the earlier draft listed
inconsistently: XOR the twelve limbs, OR the twelve differences together,
branch on the accumulated result. That is `12 + 11 + 1 = 24` instructions and
the table charges 24, not 15. A node also carries a 128-bit tag derived from the
same product, so the 96-byte compare runs about once per operation; a false tag
match needs a `2^-247` coincidence and even one spurious full compare per
operation is inside the 1.59x factor.

The probe count is the only average-case quantity in this candidate, and it is
handled in `## Matching cost and success probability`: the table is chained, so
the total link count of an attempt is a sum over *pairs* of keys sharing a
bucket, and the multiply-shift family bounds each pair by `2/2^119` without any
independence between operations. The expectation `1.504704*N` and the global
budget give `Pr[budget exceeded] <= 1.504704/12 = 0.125392` by Markov's
inequality. The earlier draft instead asserted a geometric probing law with
ratio equal to the load factor and derived a `2^-1.6e17` tail from it; that law
does not hold for open addressing and no such tail is claimed here. The global
counter is what bounds the work of every attempt, and the ledger charges it in
full at 560 per operation.

### P3 — one scalar chain (1 per `zeta` power), charged 300

| step | instructions |
| --- | ---: |
| load 8 operand limbs | 8 |
| 4x4 schoolbook (16 x `UMULH`+`MADD`+carry) | 40 |
| Montgomery reduction, 4 rounds x 6 | 24 |
| conditional subtraction of `r` (4 x `CSEL`+`SUB`, borrow) | 8 |
| store 4 limbs | 4 |
| loop control, amortised | 4 |
| **itemised sum** | **88** |
| **charged** | **300** (3.4x) |

This is also below the package's earlier comparison point, the measured
383-instruction 6-limb base-field multiply with full reduction: a 4-limb
Montgomery multiply is a strictly smaller operation, and 300 exceeds it by 3.4x.

### P4 — digit extraction and loop control (1 per scalar multiplication), charged 100

The scalar is 255 bits, split into three 85-bit digits. A digit starts at a bit
offset that is not limb-aligned, and the window offset multiplies an 85-bit
value by 96, so both are 2-word quantities and are priced in 2-word lanes.

| step | instructions |
| --- | ---: |
| 3 x (digit extract across a limb boundary, 2-word) | 12 |
| 3 x (offset `digit*96`, 2-word) | 18 |
| 3 x (window control and accumulator register select) | 9 |
| accumulator loop control and address updates | 10 |
| **itemised sum** | **49** |
| **charged** | **100** (2.0x) |

The earlier draft priced the same three digit steps at 21 instructions as if
each digit and offset were single-word; the corrected itemisation is 49, still
covered by the charge of 100.

### The charged numbers do not carry the claim

Doubling all three contested allowances moves the claim by +0.138 bits and 8x
moves it by +0.769 bits (the ledger arithmetic is in `experiment.py` and in the
recomputation below), so this derivation buys *support*, not score:

| comb | hash | chain | claim log2 |
| ---: | ---: | ---: | ---: |
| 200 | 560 | 300 | 119.906712 (submitted) |
| 400 | 1120 | 600 | 120.044872 |
| 1600 | 4480 | 2400 | 120.675233 |

### Layout consistency

The comb entries are 96-byte affine points, produced by batch normalisation from
a 144-byte Jacobian staging buffer over one block of `2^18` entries = 2^25.17
bytes, with a prefix/inverse array of `2^18` field elements = 2^25.17 bytes and
an affine output buffer of `2^18` points = 2^25.58 bytes. The three together are
2^26.2 bytes, i.e. below the `2^27`-byte block/code/scratch term of the memory
ledger, and they are transient inside one block: only the normalised 96-byte
entries persist in the `6*(2^85-1)*96`-byte tables. So the 96-byte comb stride
and the 144-byte Jacobian stride do not contradict each other, and the block size
appears in no counted quantity except the (negligible) per-block inversion.

## Complete cost ledger

The track charges `W = H + R/C` with `C = 2000` for BLS12-381 G1 and
`H` counting reference additions. `R` is in scalar AArch64 instruction
equivalents; loads, stores, addressing, control, hashing and memory transfers are
included. No reference operation and its implementation are double counted.

| Component | Count | Reason |
|---|---:|---|
| `H_pre` | 696341272098026404630757868 | 6*(2^85-1) table additions + 510 doublings, at 3 ref ops per table entry |
| `H_main` | 894487539658787414315834285001539584 | 2 additions x N fixed-base scalar multiplications |
| `H_p2` | 395650 | 1025 variable-base phase-2 scalar multiplications (`T'` shift plus 1024 baby steps) at 384 additions, plus 1025 fixed-base giant/verification multiplications at 2 additions |
| `H_exc` | 1024 | rounding and branch slack; table doublings are in `H_pre` |
| **H** | **894487540355128686413860689632694126** | |

| Instruction component | Count `R` | Reason |
|---|---:|---|
| comb lookups | `N*3*200 = 268346261897636224294750285500461875200` | P1: itemised 64 in 2-word lanes, charged 200 (3.1x); 3 lookups per scalar multiplication |
| hash operations | `N*560 = 250456511104460476008433599800431083520` | P2: itemised 203 + 13 per extra link; 560 covers the worst case `203 + 12*13 = 359` under the global budget, 1.56x |
| false full-point comparison allowance | `2^20*24 = 25165824` | hard global cap of `2^20` full 96-byte comparisons at 24 instructions each; abort on exceed, Markov bound `2^-143` |
| scalar chains | `N*300 = 134173130948818112147375142750230937600` | P3: itemised 88, charged 300 (3.4x), for one 4-limb Montgomery multiplication plus canonical reduction per `zeta` power |
| digit extraction/loop | `N*100 = 44724376982939370715791714250076979200` | P4: itemised 49 in 2-word lanes, charged 100 (2.0x); three 85-bit digits, accumulator addressing |
| table stores | `6*(2^85-1)*100 = 23211375736600880154358578600` | writing 96-byte affine entries during precomputation |
| hash-table zeroing | `2^122 = 5316911983139663491615228241121378304` | `2^122` bytes of bucket heads is `2^118` 16-byte zero stores; charged `2^122` for loop control and node-arena cursor initialization (16x overcount) |
| factorization/primitive root/setup | `2^40 = 1099511627776` | trial division to `2^28`, primitive-root tests (smallest root `7`), parameter arithmetic, sequential input access |
| phase-2 support | `2^22 = 4194304` | modular exponentiation of the `zeta^{...}`/`omega^{...}` scalars, the 1024-item chained table and loop overhead beyond counted additions |
| output/verification support | `2^12 = 4096` | encoding `s`, final comparisons |
| **R** | **703017192940205222394566851796196658600** | |

Therefore

```text
R/C = 351508596470102611197283425898098329.3
W   = H + R/C = 1245996136825231297611144115530792455.3
log2 W = 119.9067110112798959862561830629684077409...
```

The submitted bound is `2^119.906712`, which is at least `2^(log2 W)`; the
rounding margin is under `10^-6` bits. All preprocessing, table construction,
the `2^122`-instruction zeroing, the factorization and the phase-2 support are
inside `W`. No component is charged twice: the additions of `H` are the group
operations, and every instruction row prices only the non-group work.


## Memory, preprocessing and advice

Peak memory is dominated by the chained matching table: `2^119` bucket heads of
16 bytes (`2^123` bytes) plus at most `m` nodes of 128 bytes (`128*m` bytes). The
six comb tables add `6*(2^85-1)*96` bytes (`2^94.17`), and the block buffers,
the two retained input points (`A_d`, `P0`), code and scratch add less than
`2^27` bytes. The exact accounting value is

```text
16*(2^119) + 128*m + 6*(2^85-1)*96 + 2^27 = 39257425257643444948473998550610476480 bytes
log2 = 124.884305065...  <= 125.000000.
```

The declared bound is **`2^125` bytes**, raised in this revision: the previous
`2^124.807355` bound was derived from 8-byte bucket heads, and 8-byte heads cannot
reference the `2^124.43`-byte node arena any more than 64-bit next references
could. Charging 16-byte heads adds `2^122` bytes, which is exactly the amount by
which the old bound is exceeded; the recomputed total is `2^124.884305` and the
declared bound is rounded up to `2^125`.

This memory is astronomical and is required by the BSGS table; it is reported as
a bound, not a practical resource budget. Memory has no score contribution.

Preprocessing is declared as the fixed-base comb construction, the factorization
and primitive-root search, and the hash-table initialization:

```text
H_pre + (R_stores + R_setup + R_init)/2000
= 2658456699516791712134459378253660.34
log2 = 111.03421609952...  <= 111.034217.
```

It is included in `W`. No table or search result is embedded in the candidate
files or supplied as nonuniform advice: `nonuniform_advice_log2_bytes = 0`. The
input powers are public challenge inputs; their generation is organizer work and
is not charged, but every attack-side access to `A_d` and every runtime table is
charged and included in peak memory. All tables depend on the target (`A_d`) or
are rebuilt from public constants; nothing is reused across targets.

## Matching cost and success probability

Fix any valid challenge instance (any `P0` and any `A_d`). The only randomness
in the attempt is the pair `(a,b)` of the multiply-shift hash, drawn uniformly
at the start from the attempt's fresh coins. Everything else is deterministic.
The attempt outputs either the correct `s` or FAILURE.

**Work is bounded by construction, not in expectation.** The attempt carries two
global counters, both reset at the start of the attempt:

* the **link counter**, incremented once per chain link walked by any operation;
* the **comparison counter**, incremented once per full 96-byte point comparison.

When the link counter would exceed `12*N`, or the comparison counter would exceed
`2^20`, the attempt returns FAILURE. No operation is capped individually, no chain
is truncated and no comparison is skipped. Hence the matching step costs at most

```text
N*203 + 12*N*13 + 2^20*24 = 359*N + 25165824 instructions <= 560*N,
```

which the ledger charges at the rounded `560*N` (P2) plus the explicit `2^20*24`
row; `203` is the itemised P2 base (it already contains one 24-instruction full
comparison per operation) and `13` the itemised cost of one extra link. No
expectation enters the work bound: both counters are if-then rules, and the ledger
prices both budgets.

**Failure probability.** Walking a chain visits exactly those stored keys that
share the walked key's bucket, so the total number of links walked by the whole
attempt is at most

```text
Links  =  N  +  P_insert  +  P_lookup,
P_insert  <=  #{(i,j) : i<j, both baby points, h(x_i) = h(x_j)},
P_lookup  <=  #{(i,j) : i baby point, j giant point, h(x_i) = h(x_j)}.
```

For each of these pairs, `Pr[h(x_i) = h(x_j)] <= 2/2^119` by the 2-universality
of multiply-shift, and *no independence between different pairs is used*: the
expectation of a sum of indicators is the sum of their probabilities. With
`C(m,2)` baby pairs and `m*g1` baby-giant pairs and `N = m + g1`,

```text
E[P_insert] <= C(m,2) * 2/2^119 = 75241790831034500209385817640088634.  < 7.525e34
E[P_lookup] <= m*g1   * 2/2^119 = 150483581662069000418771635280177271. < 1.505e35
E[Links]    <= N + 2.2573e35    = 672969142322497207786074595421035697
            <= 1.504704 * N.
```

**Second abort event: full comparisons.** A tag match that is not the true match
costs a full 96-byte comparison, and the link counter does not bound how many of
those occur. The comparison counter does. By the 247-bit bound of the previous
section, each of the at most `12*N` links matches the stored 128-bit tag of the
key it walks with probability at most `2/2^247`, no independence between links
being used (linearity of expectation over the links, exactly as for the pairs):

```text
E[false full comparisons] <= 12*N * 2/2^247 = 2^-123.99,
Pr[comparison counter exceeded] <= 2^-123.99 / 2^20 = 2^-143.99.
```

The true match contributes at most one comparison per operation and is already
inside the P2 base charge.

By Markov's inequality, `Pr[Links > 12*N] <= 1.504704/12 = 0.125392`, and this
is **exactly** the abort event of the declared algorithm: the global counter
stops the attempt precisely when `Links` would exceed `12*N`, and nothing is
truncated before that. (The previous revision instead capped each operation at
12 links and still applied Markov to the total. Those are different events — one
operation with 13 links aborts a run whose total is far below `12*N` — so that
implication was false and is withdrawn; no per-operation bound is claimed
anywhere in this candidate.) On the complement (probability at least `0.874608`)
every chain walk that the scans perform completed untruncated, so the fixed
complete scans run to the end and return the correct `s`; that argument
(previous section) does not depend on the hash at all. Therefore

```text
success probability >= 1 - 0.125392 - 2^-143.99 = 0.874608 > 0.85 > 0.39.
```

The submission claims `success >= 0.85` and nothing stronger. This deliberately
replaces the earlier draft's `1 - 2^-100` claim, which was not supported: a
first-moment argument cannot produce a tail bound, and the geometric probing law
it invoked does not hold for open addressing. Chaining is what makes the
elementary pairwise argument sufficient — the cost of one attempt is the number
of colliding *pairs*, and a first moment bounds that directly, with no probing
model, no load-factor law and no independence assumption. No retry is used in
the claim; the failure probability above is the claimed failure probability. A
retry with fresh `(a,b)` would succeed with at least the same probability, but
its work would have to be charged in full, so no retry budget is claimed.

For completeness, the same scheme with a chained table of 1024 items serves
phase 2; the whole phase-2 matching cost is at most `1025*353` instructions,
below the `2^22` phase-2 support line of the ledger.

## Why this bound, and rejected alternatives

The phase-1 search is the BSGS of Cheon's Lemma 1 with the supplied `d = 2^20`.
Cheon's `O(sqrt((r-1)/d) + sqrt(d))` counts group exponentiations; replacing each
known-scalar exponentiation by a three-window fixed-base lookup+addition chain
reduces the per-step cost from a full variable-base scalar multiplication to two
reference additions. The `sqrt(d)` phase is `2^10` and negligible. The choice
`d = 2^20` minimizes the phase-1 modulus even without the divisor restriction:
for any exponent `d` the power `A_d` determines `a` modulo
`(r-1)/gcd(d, r-1)`, and `gcd(d, r-1) <= d <= n = 2^20`, so the phase-1 modulus
is smallest at `d = 2^20`; a smaller `d` increases `M1` and the work, while phase
2 stays at `sqrt(gcd(d, r-1)) <= 2^10`. The `p+1` variant applies to divisors of
`r+1` instead; any divisor `d <= n = 2^20` of `r+1` gives
`(r+1)/d >= (r+1)/2^20 > M1`, so its phase-1 search is larger and it is not used.
A second supplied power `A_e` with `e <= n` yields only a second residue modulo
`(r-1)/gcd(e, r-1) >= M1`; recovering it independently costs another BSGS of at
least the same size (another `N` scalar multiplications), and combining the two
residues leaves an ambiguity of at most `2^20`, no better than the `2^10` phase 2.
So no multi-power combination is used, and the other supplied powers are not
needed. The map `a -> [zeta_d^a]P0` used in phase 1 is not a group
homomorphism on the index, so Pollard-rho-style collision search does not replace
the BSGS table here; no memoryless shortcut is claimed or used.

This is not claimed to be the best possible attack, a security lower bound, or a
practical break. The bound is analytic and could be tightened by a more detailed
instruction model; the declared charges are deliberately conservative in the
lookup, hash, scalar-chain and initialization terms, which together contribute
about 0.46 bits to `log2 W`.

## Evidence and reproducibility

`experiment.py` (stdlib only) is declared in `submission.json`. On a toy
prime-order group of about 31 bits it runs the same two-phase Cheon algebra with
the same fixed-window lookup/addition structure, full-point matching, both
phases, and the final check, for seeded random and boundary scalars; it also
checks the production 3x85-bit digit structure on boundary and seeded sample
scalars (the no-exception proof itself is analytic and covers every
`0 < x < r`), verifies the factorization of `r-1` and that `7` is
the smallest primitive root modulo `r`, recomputes the production `M1`, `m`,
`g1`, exact integer ledger, score bound, preprocessing and memory bound, and
asserts the claimed bound is not exceeded. The
experiment is evidence only at toy scope; it does not execute a production-size
attack or validate production memory. The repository's unit tests check the
protected BLS12-381 profile constants, the subgroup-order factor `2^20`, the
supplied auxiliary powers and the pinned pairing conventions.

No production-size run is claimed. The ledger is analytic; all integer
components are stated above so a reviewer can recompute the bound. The main
uncertainty is the instruction accounting for the large-table lookups and the
hash operations, charged 1.59x above the itemised global-budget worst case; the dominant `2N`
additions are independent of those charges.

## References

- J. H. Cheon, "Discrete Logarithm Problems with Auxiliary Inputs", Journal of
  Cryptology 23 (2010), Lemma 1, Theorem 1 and the KKM precomputed-table
  speed-up cited there. The identities used above are equations (1)-(3) of the
  p-1 case, adapted to the supplied `A_d`.
- Track rules, target profile and `classical-group-work-v2` cost model in this
  repository; the BLS12-381 G1 `6M` addition convention with `C = 2000`.

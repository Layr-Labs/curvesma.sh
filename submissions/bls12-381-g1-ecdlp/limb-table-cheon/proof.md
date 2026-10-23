# Four-limb fixed-base Cheon rectangle with a deterministic merge on BLS12-381 G1

## 1. Claim

Track `bls12-381-g1-ecdlp`, scored objective `time_log2`, direction down.
Input: `T = [s]P0 = A_1` and the powers `A_i = [s^i mod r]P0`, `1 <= i <= n`,
with `s` uniform in `1..r-1`, scored at `n = 2^20`.

The complete classical work of the attack below at the scoring length is

```
W(2^20) = 378962870133214382517139694889755716449 / 500
        = 684151412907864145046948343712561364
        + 147548654717129239974662092133900137796 / 2000
```

target-group-operation equivalents, so

```
log2 W(2^20) = 119.18953982460239...
submitted score = 148987/1250 = 119.1896
```

and the submitted rational is a ceiling: `2^(148987/1250) >= W(2^20)` holds,
and the same inequality at `148986/1250` fails. The bound is a worst case over
the attack's own choices; no expected-value estimate is used anywhere.

```
success_probability      >= 39/100      (exact integer count, no coins)
time_unit                = target-group-operations
memory_log2_bytes        = 122.3
preprocessing_log2       = 69           (included in W)
nonuniform_advice_log2   = 72.6         bytes, included in peak storage
```

The bound holds for every `n` with `2^20 | n` and `2^20 | r-1` (the second
condition holds for this curve, Section 2). The attack reads only `A_1` and
`A_n`; it never scans the rest of the array. No heuristic, no randomized hash
family, no Poisson model and no algorithm coins are used: correctness, coverage
and the probe structure are deterministic integer statements about fixed data.

## 2. Fixed constants and how they are derived

Let `r` be the G1 subgroup order from the protected profile,
`r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`.

**BLS12 polynomial identity.** With `x = -0xd201000000010000 = -15132376222941642752`
(the profile's `parameter_x`), `r = x^4 - x^2 + 1`, hence

```
r - 1 = x^2 (x-1) (x+1).
```

The three 64-bit factors factor completely by trial division:

```
|x|   = 0xd201000000010000 = 2^16 * 906349 * 254760293
|x|-1 = 0xd20100000000ffff = 19 * 125527 * 2508409 * 2529403
|x|+1 = 0xd201000000010001 = 3 * 11 * 10177 * 859267 * 52437899
```

so `r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 *
2529403 * 52437899 * 254760293^2`. This factorization is obtained by factoring
three 64-bit integers; its cost is charged in preprocessing (Section 7).

**The group and the d-th power map.** Set `d = 2^20` and `m = (r-1)/d`:

```
m = 50006747412801924209068050869165387952509453297164571592906626415194112
```

`d | r-1`, so the map `F_r^* -> F_r^*, u -> u^d` is exactly `d`-to-one onto the
subgroup `H` of `d`-th powers, `|H| = m`, and `H` is cyclic. `3 | m` because
`3 | r-1` and `3` does not divide `2^20`.

**A generator of H.** `7` is the least primitive root of `F_r^*`: for each of
the twelve primes `q | r-1` above, `7^((r-1)/q) != 1 (mod r)`. Then

```
zeta = 7^d mod r
     = 0x628fc184374a1219e85d9e943d57d0c20bad4d8f813859f6121b51c927eaa4dc
```

has order exactly `m`: `zeta^m = 1` and `zeta^(m/q) != 1` for every prime
`q | m`. The twelve primality/order tests are a fixed handful of modular
exponentiations, charged in preprocessing.

**The mu_3 orbit in exponent space.** Let `lambda` be the GLV eigenvalue with
`lambda^2 + lambda + 1 = 0 (mod r)`, so that the endomorphism
`phi(x,y) = (beta x, y)` acts on G1 as `[lambda]`. Since `lambda` has order 3 and
`3 | m`, `mu_3` is contained in `H`, so `lambda = zeta^c` for some `c` with
`3c = 0 (mod m)`, i.e. `c` is `m/3` or `2m/3`. Because `{c, 2c} = {m/3, 2m/3}`
as sets, the `mu_6 = < -1, phi >` orbit of the exponent `k` is exactly

```
Orb(k) = { +-k, +-(k + m/3), +-(k + 2m/3) }  (mod m)
```

and this description does not require knowing which of the two values `c` takes.

## 3. The reduction (Cheon, d | r-1)

`e = s^d mod r` lies in `H` and, since `s` is uniform on `Z_r^*` and `u -> u^d`
is `d`-to-one onto `H`, `e` is uniform on `H`. Write `e = zeta^k` with `k`
uniform in `[0, m)`.

`A_d = [s^d]P0 = [e]P0` is supplied because `d = 2^20 = n` at the scoring
length. Recovering `k` therefore recovers `e = zeta^k` as a scalar with one
modular exponentiation, and `s` is then determined up to the `d`-th roots of
unity: the solution set of `u^d = e` in `F_r^*` is `s * mu_d`, `|mu_d| = d`.

Phase 1 recovers `k` by a baby-step giant-step rectangle inside `H`; Phase 2
enumerates the `d` roots and verifies each against `T`.

## 4. Rectangle, orbit quotient, and the exact coverage count

Take the rectangle radix `u = G` and both side lengths `G`, where

```
G = 57012617742321962895928165416718183
```

Since `{ i*G + j : 0 <= i,j < G } = [0, G^2)`, the rectangle enumerates the
integer interval `[0, G^2)` of exponents. Quotienting by the `mu_6` orbit adds
the six sign/shift images of every exponent, so the covered set is

```
S = [-(G^2-1), G^2-1]  union  (m/3 + [-(G^2-1), G^2-1])
                       union  (2m/3 + [-(G^2-1), G^2-1])
```

which has exactly `|S| = 6G^2 - 3` residues, provided the three intervals are
disjoint, i.e. `2G^2 - 1 <= m/3`. That holds:

```
2G^2 - 1 = 6500877163664250147178846612991500515429355996472862416604868085642977
m/3      = 16668915804267308069689350289721795984169817765721523864302208805064704
```

and the coverage target holds as an integer inequality:

```
6G^2 - 3 = 19502631490992750441536539838974501546288067989418587249814604256928931
100 * (6G^2 - 3) >= 39 * m
```

so `Pr[k in S] = (6G^2-3)/m >= 39/100` exactly, for uniform `k`, with no
algorithm coins, no restarts and no failures to repair. This is a lower bound:
`G` is chosen as a round number, and the true covered fraction is
`(6G^2-3)/m`, slightly above `39/100`.

Matching uses the full canonical orbit key, so a key collision is never a false
positive (Section 5.4).

## 5. The algorithm

### 5.1 Fixed-base tables (four limbs)

A scalar below `r` is four 64-bit limbs, `a = a0 + a1*2^64 + a2*2^128 + a3*2^192`.
For a base point `B` build four affine tables

```
T_i[t] = [ t * 2^(64 i) ] B ,   0 <= t < 2^64 ,  i = 0,1,2,3
```

(each table is an arithmetic progression, so it is built with `2^64` batched
additions; `T_3` needs `2^63` entries because `a3 < 2^63`). Then

```
[a]B = T_0[a0] + T_1[a1] + T_2[a2] + T_3[a3]
```

is three batched incomplete affine additions, with no doublings at run time.

Two bases are needed, and they are separate tables:

* `B = P0` for the baby side, `[zeta^j]P0`;
* `B = A_d` for the giant side, `[zeta^(-iG)] A_d`.

This is the correction the previous ledger lacked: the giant points are
`[e * zeta^(-iG)]P0 = [zeta^(-iG)] A_d`, which is **not** a known multiple of
`P0`, so it cannot be read from the `P0` tables. The `A_d` tables are built from
the supplied point `A_d` and their construction and storage are charged
(Sections 6 and 7).

All tables store affine points, two 48-byte coordinates each, 96 bytes per
entry. No compressed encoding is used, so no decompression is needed or charged.

### 5.2 Baby side

For `j = 0..G-1` evaluate `[zeta^j]P0` with the `P0` tables. The scalar is
maintained multiplicatively, `a_j = a_(j-1) * zeta mod r`, one four-limb
Montgomery product per point; it is four limbs of the current scalar, so the
four limb values address the four tables directly. Write the record
`(key, j)` where `key` is the canonical orbit key of the point (Section 5.4).

### 5.3 Giant side

For `i = 0..G-1` evaluate `[zeta^(-iG)] A_d` with the `A_d` tables. The scalar
is maintained multiplicatively by the fixed radix `zeta^(-G) mod r`, again one
four-limb Montgomery product per point. Write the record `(key, i)`.

### 5.4 Canonical orbit key and matching

The `mu_6` orbit of a point `Q = (X,Y)` is `{+-Q, +-phi(Q), +-phi^2(Q)}` with
`phi(X,Y) = (beta X, Y)`. Negation sends `Y` to `p-Y`; `phi` preserves `Y`. So
`min(Y, p-Y)` is constant on the orbit. Conversely, if `min(Y,p-Y) = min(Y',p-Y')`
then `Y' = +-Y`, hence `Y'^2 = Y^2`, hence `X'^3 = X^3` by `Y^2 = X^3 + 4`, hence
`X' = beta^i X` and `Q'` is in the orbit of `Q`. The key is therefore a
bijection from orbits onto its image; it costs a six-limb compare, a conditional
select and the store, with no field multiplication.

Both sides are then matched deterministically:

1. Radix-sort the `2G` records (baby and giant together) on the 381-bit key,
   twelve passes of 32-bit digits, externally, sequentially.
2. Merge the two sorted runs; a run of equal keys yields the matched pairs.

No hash table, no bucket array, no occupancy model and no overflow probability
is used, so there is no probabilistic claim to support: the sort is exact for
every input. Ties are broken left-first; every record is retained, and a key
that occurs `b` times on one side and `c` times on the other produces `b*c`
candidate pairs, all of which are lifted.

### 5.5 Lifting a match

A matched pair `(i, j)` means `Orb(k - iG) = Orb(j)`, i.e. `k - iG` is one of
the six members of `Orb(j)`. The six candidate exponents

```
k_s = iG + s mod m ,   s in { +-j, +-(j+m/3), +-(j+2m/3) }  (mod m)
```

are tested in order: compute `e_s = zeta^(k_s) mod r` (one exponentiation with
a `k_s < m < 2^235` exponent), and accept `k_s` if `[e_s]P0 = A_d`, three
batched additions with the `P0` tables. Exactly one `s` survives, because
`e = zeta^k` is a single value. Six exponentiations and six point checks are
charged per match; at most `6G` matches can arise (each giant record meets at
most one baby orbit, and each match has six lifts), and the ledger charges the
worst case.

### 5.6 Phase 2: the d-th roots of e

Given `e = s^d mod r`, recover `s`. Write `r-1 = 2^32 * w` with `w` odd, and let
`omega = 7^w mod r`, an element of order `2^32`. Then:

1. `E = e^(2^12) = s^(2^32) = s_odd^(2^32)`, where `s = s_odd * s_2` is the
   decomposition into an odd-order part and a 2-power part.
2. `s_odd = E^a mod r` with `a = (2^32)^(-1) mod w`, since `s_odd^w = 1`.
3. `E_2 = e * s_odd^(-2^20) = s_2^(2^20)` has order dividing `2^12`.
4. `omega_2 = omega^(2^20)` has order `2^12`. Solve `omega_2^c = E_2` by
   baby-step giant-step in a group of order `2^12` (`2^6` steps), and set
   `rho_2 = omega^c`; then `rho_2^(2^20) = E_2`.
5. The complete solution set is `s_t = s_odd * rho_2 * (omega^(2^12))^t mod r`
   for `t = 0..2^20-1`, since `mu_d = <omega^(2^12)>`. Test `[s_t]P0 = T` for
   each `t` with three batched additions from the `P0` tables; output the `t`
   that verifies (exactly one does).

Steps 1-4 are a fixed number of modular exponentiations and a 64-step BSGS;
step 5 is `2^20` point tests. Both are charged in Section 6. Nothing in Phase 2
is probabilistic.

## 6. Cost ledger

Reference operation: batched incomplete affine BLS12-381 G1 addition, `6M`,
`C = 2000`. Write `H` for reference operations and `R` for other scalar AArch64
instruction equivalents; the model's cost is `H + R/C`, every item once.

**Per point of Phase 1 (each of `2G` points).**

| Item | Charge |
| --- | ---: |
| three fixed-base additions | 3 reference operations |
| degenerate-pair branch on each addition (doubling or identity skip) | 3 reference operations |
| four table fetches: `4 * 96 = 384` bytes = 48 eight-byte loads, 4 addresses | 52 instructions |
| scalar update: four-limb CIOS Montgomery product | 256 instructions |
| digit addressing, four limbs | 12 instructions |
| batch-inversion layer traffic: 3 layers, `3 * (96 store + 96 load)` | 72 instructions |
| orbit key (compare, select) and record store (key, index, tag) | 29 instructions |
| loop control and branch | 8 instructions |
| **instruction subtotal per point** | **429, charged as 1024** |
| radix sort and merge: 12 passes at 20 instructions plus 30 merge | 270 instructions |

The instruction subtotal is charged at 1024, a factor 2.39 above the itemized
count. The additions are charged at two reference operations each: one for the
reference addition itself and one worst-case allowance for the degenerate-pair
branch, which is a doubling or an identity skip. That allowance is a bound, not
an expectation: it is paid on every addition, including the ones that never
branch.

Per point in reference units: `6 + 1294/2000 = 6.647`.

**Phase 1 total.** `2G = 114025235484643925791856330833436366` points:

```
H_1 = 6 * 2G                    = 684151412907863554751137985000618196
R_1 = 1294 * 2G                 = 147548654717129239974662092098466657604
```

**Phase 2.** `2^20` point tests at two reference operations each plus the same
instruction subtotal, and six lifts per match:

```
H_2 = 6 * 2^20                  = 6291456
R_2 = 1024 * 2^20               = 1073741824
```

**Preprocessing, advice construction, input decoding.**

```
H_pre = 2^69                    = 590295810358705651712
R_in  = 2^15 * n = 2^35         = 34359738368     (decoding over-charge)
```

`H_pre` covers the eight limb tables (`8 * 2^64` entries at one reference
addition each, `2^67`), their batch conversion to affine (`3` field
multiplications per entry, `383` instructions each, over `8 * 2^64` entries),
the factorization of `r-1` (three 64-bit integers, trial division), the
primitive-root and order tests (a dozen exponentiations in `F_r^*`), and the
Phase-2 setup. It is deliberately far above the itemized cost. `R_in` charges
`2^15` instructions for every supplied power even though the attack reads only
`A_1` and `A_n`, so an implementation that walks the serialized array is still
covered.

**Totals.**

```
H = H_1 + H_2 + H_pre = 684151412907864145046948343712561364
R = R_1 + R_2 + R_in  = 147548654717129239974662092133900137796
W = H + R/2000        = 378962870133214382517139694889755716449/500
```

and `log2 W = 119.18953982460239...`, so the submitted `148987/1250 = 119.1896`
is a ceiling, verified exactly by `2^148987 >= W^1250` and by the failure of the
same inequality at `148986/1250`.

**What is deliberately not discounted.** The sort is charged on every record,
not on the number of probes. The degenerate branch is charged on every addition.
Every record is written, sorted and merged even when a match is found early; the
giant loop is not aborted. Empty or trivial digits are charged as full
additions. The `6G` lift budget is charged even though at most one lift per
matched giant can verify. No preprocessing is amortized over future targets and
no work is divided between processors.

## 7. Preprocessing, advice, memory

* `preprocessing_log2 = 69`: the `H_pre` term above, already inside `W`.
* `nonuniform_advice_log2 = 72.6` bytes: the four `A_d` tables,
  `4 * 2^64 * 96 = 7083549724304467820544` bytes. They depend on the challenge
  input `A_d`, so they are advice, and they are counted in peak storage.
* `memory_log2_bytes = 122.3`: `2G` records at 56 bytes each
  (`6385413187140059844343954526672436496` bytes), the eight limb tables
  (`14167099448608935641088` bytes), the `2^32` radix counters (`2^35` bytes)
  and the merge buffers. Records dominate at `log2 = 122.264...`.

Memory carries no score contribution under this model; it is reported because
the rules require it.

## 8. Applicability and limits

* The attack needs `d = 2^20` to divide `r-1` (it does) and `A_d` to be supplied,
  i.e. `n >= 2^20`. `W(n) = W(2^20) + 2^15 * n` for `n >= 2^20`, and the score is
  stated at `n = 2^20` only.
* Only `A_1` and `A_n` are read.
* The success probability is a worst case over `s`: for every `s` whose `k` lies
  in `S`, which is at least `39/100` of the uniform range, the attack recovers
  `s` exactly. For the remaining `k` the attack fails; that is the declared
  failure probability, not an amplification target.
* No part of the attack depends on a heuristic. The only "small experiment"
  content is the integer arithmetic of Section 9, which is a check of constants,
  not a production run.

## 9. Reproduction of the constants

Every number in this proof is a deterministic function of the profile. The
following standard-library script reproduces the identity, the factorization,
the generator, the coverage inequality and the score inequality:

```python
r = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
x = -int("d201000000010000", 16)
assert x**4 - x**2 + 1 == r and x*x*(x-1)*(x+1) == r-1
d, m = 2**20, (r-1)//2**20
G = 57012617742321962895928165416718183
primes = [2,3,11,19,10177,125527,859267,906349,2508409,2529403,52437899,254760293]
assert all(pow(7, (r-1)//q, r) != 1 for q in primes)          # 7 is primitive
z = pow(7, d, r)
assert pow(z, m, r) == 1 and all(pow(z, m//q, r) != 1 for q in primes if m % q == 0)
assert 100*(6*G*G - 3) >= 39*m and 2*G*G - 1 <= m//3
H = 6*(2*G) + 6*2**20 + 2**69
R = 1294*(2*G) + 1024*2**20 + 2**35
from fractions import Fraction
W = Fraction(H) + Fraction(R, 2000)
assert (1 << 148987) * W.denominator**1250 >= W.numerator**1250      # 2^(148987/1250) >= W
assert (1 << 148986) * W.denominator**1250 <  W.numerator**1250      # and it is tight
```

The script's assertions hold; `z` prints as
`0x628fc184374a1219e85d9e943d57d0c20bad4d8f813859f6121b51c927eaa4dc`, the
coverage integer is `100*(6G^2-3) = 1950263149099275044153653983897450154628806798941858724981460425692893100`,
and `log2 W = 119.18953982460239...`. The exact score comparison is done with
integers, not floats: `2^148987 * 500^1250 >= 378962870133214382517139694889755716449^1250`
holds and the same inequality with `148986` fails.

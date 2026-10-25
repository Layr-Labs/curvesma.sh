# Frobenius-class-folded distinguished-point rho for fixed-G1 pairing inversion

Track `bls12-381-fixed-g1-invert-g2`. Claim: `time_log2 = 125.81` (12581/100), success probability at least 0.390029.
All objects below are the ones of `tracks/bls12-381-fixed-g1-invert-g2/target-profile.json`; the unit of work is one
dense Fp12 multiplication (`C = 20,000` scalar instructions). Everything the attack does is charged once.

## 0. Summary

1. `P` fixed, `g = e(P,Q0)` generates GT (prime order `r`). Write `z = g^a` and return `Q = [a]Q0` (section 1).
2. The Frobenius `x -> x^p` acts on GT as `x -> x^lambda`, `lambda = p mod r` of order 12 with `lambda^6 = -1`. It is
   a free action on `GT \ {1}`: every class has exactly 12 elements and there are `N = (r-1)/12` classes (section 2).
   Inversion is a member of the group (`sigma^6`), so conjugation is part of the folding.
3. A *class selector* costs 3 Fp multiplications, three modular subtractions and one negation: from the single Fp2
   coordinate `c1` it produces the 12 first-coordinates of `c1(sigma^k(x))`, a Frobenius-invariant class key `kappa`
   (their minimum) and a covariant twist `k*` (its argument) (section 3). The twelve images are never formed:
   `d = g0*(a-b)` is one multiplication because `gamma = g0*(1-u)` exactly, and the minimum is taken on the six stored
   values and their negations through one negation of their maximum.
4. The walk is a *lifted* adding walk: state `y` (any member of its class), multiplier `T_h`, step
   `y' = y * sigma^{-k*(y)}(T_h)`, with `h` a hash of `kappa`. The class of `y'` depends only on the class of `y`
   (Lemma 3.2), so this is a walk on the `N` classes at one dense Fp12 multiplication per step (section 4).
5. Fruitless short cycles (the `1/(12R)` per-step two-cycles and every other deterministic closed pattern) are
   removed by a *look-back window* of the last 24 class keys with deterministic probing of the next table indices.
   The rule depends only on the current class and the window, so trails that meet stay together up to a stated
   loss below `2.71*10^-4` (sections 4.3, 6). The window is a 24-slot ring of 32-bit tags `w_1 >> 32` with a
   `2^32`-byte count array; membership rejects on a positive count, so a tag collision can only reject a non-member,
   charged as a probe, never accept a repeat (section 5). Residual deterministic cycles longer than the
   window are bounded by an explicit combinatorial lemma (section 7): at most `2^-163.5` per step.
6. Trails start independently (fresh random start each, shared table), end at the first distinguished point (density
   `2^-80`), and are reported to an ordered store. Two reports with the same key give two trails; both are *replayed*
   with exponent tracking and the logarithm is recovered with the Frobenius alignment `sigma^m` (section 4.4). No
   exponents are carried in the main loop.
7. Work: `S_max = ceil(sqrt(2 mu N U)) = 68967528808679891282472564154805017972` accepted steps (`2^125.697`), where
   `U` is the exact rational upper bound on `ln(10000/6097)` of section 5; each multiplication slot
   charged `1 + 1596/20000 = 21596/20000`; all additive terms together are `2^92.361`. Total `W_max <= 2^125.81`
   (`log2 W_max = 125.80819701...`), rounded up to the claim `2^125.81` (section 5). Memory `2^53.6` bytes,
   preprocessing `2^24.15` units, no advice.

## 1. Reduction

For a nonidentity `P` the reduced pairing `e(P, .)` is nondegenerate (`r` prime), so `g = e(P,Q0)` has order `r`.
`GT = {x : x^r = 1}` is cyclic of order `r`, so `z = g^a` with `a` in `Z/r` (`a = 0` iff `z = 1`). Then
`e(P,[a]Q0) = g^a = z`. If `z = 1` the attack returns `Q = O` (valid exactly then). For a challenge with `P`
uniform among nonidentity points and `z` uniform in GT, `g` is some generator and `a` is uniform on
`Z/r` independently of `g`. The pairing convention (`pairing(Q,P)` in the pinned library) only fixes which
bilinear map is meant. The attack computes `g` once (one pairing, charged 2^30 instructions) and later `[a]Q0` and
`g^a == z` as the final check.

## 2. Frobenius classes of GT

`r` divides `Phi_12(p) = p^4 - p^2 + 1` (BLS12 parametrisation; checked numerically in `experiment.py`,
mode `algebra`: `lambda = p mod r` has `lambda^12 = 1`, `lambda^6 = -1`, `lambda^4 != 1`). For `x` in GT,
`x^p = x^lambda`. Let `sigma(x) = x^p` (the Frobenius endomorphism of Fp12, an automorphism of GT of order 12).
For `x != 1` and `0 < k < 12`: `sigma^k(x) = x` would mean `x^{lambda^k - 1} = 1`, i.e. `r | lambda^k - 1`, impossible
because `lambda` has order 12 modulo the prime `r`. So the orbit ("class") of `x` has exactly 12 elements and there
are `N = (r - 1)/12` classes, `log2 N = 251.2721`. Because `sigma^6(x) = x^{-1}`, inversion is already in the group and
no separate negation map is used.

For `y = g^alpha z^beta` we have `sigma^m(y) = g^{lambda^m alpha} z^{lambda^m beta}`, and for `T` in GT,
`sigma^{-k}(T) = T^{lambda^{-k}}`; these two identities carry all the exponent bookkeeping.

## 3. The class selector in the BLS12-381 tower

Tower: `Fp2 = Fp[u]/(u^2+1)`, `Fp12 = Fp2[w]/(w^6 - xi)`, `xi = 1 + u` (this is the pinned `Fp[w]/(w^12 - 2w^6 + 2)`
with `u = w^6 - 1`, and blst's `Fp2/Fp6/Fp12` tower with `w^2 = v`). Write `x = sum_{i<6} c_i w^i`, `c_i` in Fp2.
Since `w^p = w * (w^6)^{(p-1)/6} = gamma * w` with `gamma = xi^{(p-1)/6}` in Fp2 and `sigma` acts on Fp2 by
conjugation,

    sigma(sum c_i w^i) = sum conj(c_i) gamma^i w^i .

Hence the coefficient `c1` of `w` evolves by `tau(c) = gamma * conj(c)`: `c1(sigma^k x) = tau^k(c1(x))`. With
`zeta = gamma * conj(gamma) = xi^{(p^2-1)/6}`, `tau^2(c) = zeta * c`, and `zeta` lies in Fp and is a primitive
6th root of unity (`zeta^3 = -1`, `zeta^2 = zeta - 1`; verified in `experiment.py`). Writing `c = a + b u`
and `gamma = g0 + g1 u` (hex constants are printed by mode `algebra`), the first Fp-coordinate of `tau^k(c)` is

    V_{2j}   = zeta^j * a,                 j = 0,1,2   (V_{2j+6} = -V_{2j})
    V_{2j+1} = zeta^j * d,  d = g0*a + g1*b, j = 0,1,2 (V_{2j+7} = -V_{2j+1})

because `Re(zeta^j gamma conj(c)) = zeta^j (g0 a + g1 b)`. All twelve values follow from three Fp multiplications,
using `zeta^2 = zeta - 1` for the rest, and sign flips.

**The tower identity that makes `d` one multiplication.** In this tower `gamma` is exactly `g0*(1-u)`, i.e.
`g1 = -g0`. Since `xi = 1+u`, `gamma^6 = xi^{p-1} = (1-u)/(1+u) = -u` and `(1-u)^6 = 8u`, so
`(gamma/(1-u))^6 = -1/8`; the sixth root in this tower is the base-field element `g0 = gamma/(1-u)`, so `g0^6 = -1/8`
and `zeta = 2 g0^2`. These are exact integer congruences on the printed constants (`g1 + g0 = 0`, `8 g0^6 + 1 = 0`,
`2 g0^2 - zeta = 0`), asserted by `experiment.py` (`algebra`), so no root choice or distributional assumption enters.
Then `d = g0*a + g1*b = g0*(a-b)`: one multiplication of the difference. `gamma` here is the `xi^{(p-1)/6}` of the
tower declared above, in which the algorithm's elements are stored.

**Definition.** `kappa(x) = min_k V_k(x)` (minimum of the stored representatives under the usual order of the integer
representatives that the implementation keeps; any fixed total order on Fp works, blst's Montgomery representatives
included) and `k*(x)` = the smallest `k` attaining it.

**Lemma 3.1 (covariance).** `V_k(sigma^m x) = V_{k+m}(x)` (indices mod 12). Hence `kappa(sigma^m x) = kappa(x)`, and
when the minimum is unique `k*(sigma^m x) = k*(x) - m`. *Proof.* `c1(sigma^m x) = tau^m(c1(x))`, so
`V_k(sigma^m x) = Re(tau^k tau^m c) = V_{k+m}(x)`. A minimum over a cyclically shifted list is the same value at a
shifted position. ∎

**Lemma 3.2 (the twisted step is a class map).** For `T` in GT and `x` in GT, set `x' = x * sigma^{-k*(x)}(T)`. Then
`sigma^{k*(x)}(x') = sigma^{k*(x)}(x) * T`; and `sigma^{k*(x)}(x)` is the unique member `c` of the class with
`k*(c) = 0`, which depends only on the class. So the class of `x'` is the class of `c * T` where `c` is the canonical
representative of the class of `x`: it depends on `(class of x, T)` only. ∎

Exceptions: if two of the 12 values coincide, or `c1 = 0` (eight linear conditions on `(a, b)`: `a = 0`, `d = 0` and
`d = zeta^t a` for six `t`; probability at most `8/p < 2^-378` for a random element), `k*` is
the smallest minimizing index and the state is simply treated as an ordinary state; Lemma 3.2 may fail for it. Such
states are accounted in section 6 (`L4`).

**The minimum over the twelve values.** The twelve values are the six stored ones,
`V_0 = a, V_1 = d, V_2 = zeta*a, V_3 = zeta*d, V_4 = (zeta-1)a, V_5 = (zeta-1)d`, and their six negations
`V_{k+6} = -V_k mod p`. Hence, as a set,

    {V_0,...,V_11} = {V_0,...,V_5} union {(p - V_k) mod p : k < 6},

and the minimum and the smallest attaining index follow from two scans of the six stored values and one negation:

1. `m = min_{k<6} V_k` with its smallest index `i_m`, and `M = max_{k<6} V_k` with its smallest index `i_M`
   (10 comparisons, each a six-limb compare with a conditional move of value and index).
2. `nb = p - M` reduced into `[0, p)` (one modular negation; `M = 0` gives `0`, which is the minimum only when all six
   values are `0`, the `c1 = 0` case below).
3. `m <= nb` gives `(kappa, k*) = (m, i_m)`, otherwise `(nb, i_M + 6)` (the eleventh comparison).

This is the same pair `(kappa, k*)` as the direct scan of the twelve values with "smallest index attaining the
minimum": a minimiser in the stored half has index at most 5 and the negated half carries indices `6..11`, so a value
attained in both halves is won by the stored half's smaller index, and a value attained twice inside a half is won by
that scan's strict `<`. `experiment.py` checks the identity on all `7^6 = 117,649` six-tuples of a seven-value set of a
small prime field (`ledger`, ties and zeros included) and against the direct twelve-value scan on real GT elements and
their conjugates (`algebra`). The six negations the earlier accounting charged are not performed, so four modular
operations remain: the three subtractions `a - b` for `d`, `zeta*a - a` and `zeta*d - d`, and the negation `p - M`.

**The two extrema in one pass.** The six stored values are in memory (each kernel and each subtraction writes its
result there), so one pass carries both extrema with no spill (the allocation block of
section 5 gives the registers). Seeding loads the six limbs of `V_0` twice, with three `LDP` into the minimum and
three into the maximum, and zeroes the two index registers: 8. Each of the five remaining values is one `MOVZ`
supplying the candidate's index, three `LDP`, then a six-limb `SUBS`/`SBCS` chain of the candidate against the
minimum with six `CSEL` on the value limbs and one on its index (13), and the same chain of the maximum against the
candidate (13): 30 per value, so `8 + 5*30 = 158`. Each chain writes `XZR`, since only its carry is read. The row is
`158 + 14 = 172` instructions for the eleven six-limb comparisons with conditional move of value and index, and
`(kappa, k*)` is what two separate scans produce: the tie rule is the same strict `<` and `>`. An early exit on unequal top
limbs remains an optimisation inside this same worst case: equal top limbs cannot be excluded for an arbitrary input
pair `(a, b)`, so a hard row carries the full six-limb chain for all eleven comparisons and the faster typical path
is not claimed.

**The eleventh comparison, exactly.** `m <= nb` is the carry of one six-limb chain whose operands are ordered so
that the carry alone is the predicate: `SUBS`/`SBCS` of `nb` against `m` leaves `C = 1` exactly when `nb >= m`, i.e.
exactly when `m <= nb`, equality included. No condition here reads `Z`, which after such a chain is the last limb's
flag and nothing else, and a chain of `m - nb` would say only `m >= nb`. Six `CSEL` with condition `CS` take the
minimum (or the negated maximum), one `CSEL` takes `i_m` (or `i_M + 6`) and one `ADD` forms `i_M + 6`: 14.
`experiment.py` executes this tail with the negation on a register-file machine (188 instructions per tuple, 182 when
the maximum is `0` and the negation skips) and checks the result against the `<=` reference and the direct
twelve-value scan on 245 small-field tuples (51 with `m = nb`, one all-zero) and 40 production tuples. The negation is
`0 - M` in x12-x17: the chain, then `B.CS` on its carry (taken exactly when `M = 0`, which leaves `0`), else `ADDS`
and five `ADCS` of `p` from x0-x5, which held the dead candidate; `ADC` alone would not update the carry, and no
mask is formed.

**Instruction count of one selection (blst, AArch64, prices of `tracks/docs/COST_MODEL.md`).**
3 `blst_fp_mul` x 383 = 1,149; three modular subtractions at 22 and one negation core at 13, plus 9: one `LDP` x3
load of `p` in each of the three blocks (`p` is a fixed frame offset, and no register survives a kernel or the
extrema), 88; the two extrema in one pass and the eleventh comparison, `8 + 5*30 + 14` = 172. Total
`1,149 + 88 + 172 = 1,409`
(`experiment.py`: `selector_instructions_total`), the three rows that open section 5's table. The Montgomery constants
`zeta, g0` are two stored Fp elements; no conversion out of Montgomery form is needed.

## 4. The algorithm

### 4.1 Parameters and primitives

| symbol | value | meaning |
| --- | --- | --- |
| `R` | `2^14` | table indices |
| `w` | 24 | look-back window (32-bit tags of the last 24 states of the current trail, including the current one) |
| `A` | 8 | probes per step |
| `b` | 80 | distinguished-point bits |
| `L_cap` | `2^87` | multiplication slots per trail, accepted or rejected (about 128 mean trail lengths) |
| `mu` | 11/10 | allowed walk inefficiency (H1) |
| `D_max`, `ST_max`, `REC_max` | `2^47`, `2^47`, 16 | caps on reports, starts, recoveries |

Coins: a 256-bit key `K` from the algorithm's randomness; every other random object is ChaCha20 output under `K`
with domain separation and counters. Scalars are 384-bit outputs reduced modulo `r` (statistical distance from
uniform at most `r/2^384 < 2^-129` each). The class-key hash is

    khash(kappa, salt) = mix(...mix(mix(salt) xor limb_0) xor limb_1 ... xor limb_5),
    mix = the splitmix64 finalizer,

on the six 64-bit limbs of `kappa`, for `w_1 = khash(kappa, s_1)` and `w_2 = khash(kappa, s_2)` from two coin words.
The table index is `hidx(kappa) = w_1 mod R` and `kappa` is a distinguished point iff `w_1 = 0` and `w_2 mod 2^16 = 0`
(80 bits); both are functions of the class. The window of section 5 stores `w_1 >> 32` as each entry's 32-bit tag, so the look-back membership computes no
word beyond the two the index and the test need.

### 4.2 Setup (preprocessing)

1. `g = e(P,Q0)`. If `z = 1`, output `O`.
2. For `h < R`: `(u_h, v_h)` = ChaCha20 scalars; `T_h = g^{u_h} z^{v_h}` by square-and-multiply (squarings charged as
   full multiplications: at most `2*509 + 1` multiplications); `TW[h][k] = sigma^{-k}(T_h)` for `k < 12` by up to 11
   Frobenius applications, charged one unit each. The table is `12 R` elements (113 MB).

### 4.3 One trail

Trail `i` (`i = 1, 2, ...`; a counter):

1. Start: `(alpha_0, beta_0)` = two ChaCha20 scalars with domain `(start, i)`; `y_0 = g^{alpha_0} z^{beta_0}`
   (two square-and-multiply exponentiations and one product, at most 1,019 multiplications). If `y_0 = 1` or
   `c1(y_0) = 0` the trail is dropped. Compute `(kappa_0, k_0) = Sel(y_0)`; window `W = (kappa_0)`, seeded into the
   ring as `tag(kappa_0)` in all 24 slots with count 24 after clearing the previous trail's tags, with `p = 0`
   (section 5), so the window is the one-element set `{kappa_0}`.
2. At state `n`: if `kappa_n` is distinguished, **report** `(kappa_n, i)` and end the trail. If the trail's slot counter
   `c` (incremented by one for every multiplication in step 3, accepted or rejected) has reached `L_cap`, end the trail
   without a report.
3. Otherwise `h = hidx(kappa_n)` and for `a = 0, ..., A-1`: `y' = y_n * TW[(h+a) mod R][k_n]` (one dense
   multiplication); `(kappa', k') = Sel(y')`; if `tag(kappa')` is in the window's tag set (its count is positive),
   this probe is rejected; else accept it: `y_{n+1} = y'`, append `kappa'` to `W` (dropping the oldest beyond `w`),
   `n += 1`, go to 2. The append of section 5 is the 10-instruction counted ring write at slot `p` plus the bounded
   slot update; a candidate is rejected when its 32-bit tag has a positive count, a rule exact on every repeat and
   over-rejecting only on a tag collision, which is charged as an escape (section 5). If all `A` probes are rejected the trail ends without a report.

   The identity needs no special test: given `y_n`, the table element is uniform in `GT` (section 2 and H3, since
   `(u_h, v_h)` is uniform on `(Z/r)^2` and so is `u_h + a v_h` for `a != 0`), so `Pr[y' = 1] = 1/r` per slot and at
   most `2^-128.857` over the run; the selector then returns `kappa = k* = 0`, the next step is ordinary and no
   recovery is affected (the identity entry of `L4`, section 6).

By Lemma 3.2 the class of `y_{n+1}` is the class of `c_n * T_{h+a}` where `c_n` is the canonical representative of
the class of `y_n`; so, **if the first probe is accepted**, the class sequence is a deterministic function of the
current class alone, `F(C) = class(c * T_{hidx(kappa(C))})`. A probe beyond the first happens when the first
candidate class lies in the last-24-class window (a closed deterministic pattern of length at most 24, section 7), or on a 32-bit tag collision (section 5).
The decision depends only on the current class and the 24 previous class keys.

**Exponent invariant (used only in replay).** `y_n = g^{alpha_n} z^{beta_n}` with
`(alpha_{n+1}, beta_{n+1}) = (alpha_n, beta_n) + lambda^{-k_n} (u_{h'}, v_{h'})`, `h'` the accepted index, because
`sigma^{-k}(T) = T^{lambda^{-k}}`.

**Coalescence.** Let two trails reach the same class `C` at states `n` and `n'`. The next accepted class is `F(C)`
for both, unless it lies in one trail's window (probability about `1/(12R)` per step each, section 6, `L2`), and the
windows are identical after 24 further steps. The tag rule adds at most `24*2^-32` per candidate to that rejection probability and can only reject more, so the
two trails coincide after the meeting except with probability at most `2*24*nu`,
`nu = 1.1/(12R) + 8*24*2^-32 < 2^-17.43`, hence both reach the same distinguished point.

### 4.4 Store, collision, recovery

A global ordered map (a red-black tree) from the full 48-byte key `kappa` to the trail index. On a report `(kappa, i)`:
if `kappa` is absent, insert; if present with trail `i' != i`, call **Recover(i', i)**; otherwise ignore (a trail
reports at most once). Keys are compared in full; there are no fingerprints.

**Recover(i', i).** Replay trails `i'` and `i` (same procedure, same coins, now maintaining `(alpha, beta)` by the
invariant above). Let `(y_A, alpha_A, beta_A)` and `(y_B, alpha_B, beta_B)` be the endpoints. If their keys differ,
discard (false store match). Find `m` in `0..11` with `sigma^m(y_A) = y_B` (apply Frobenius up to 11 times to `y_A`,
compare all 12 Fp coefficients); if none, discard. Since `sigma^m(y_A) = g^{lambda^m alpha_A} z^{lambda^m beta_A}`,

    lambda^m alpha_A + a lambda^m beta_A  =  alpha_B + a beta_B   (mod r),
    a (lambda^m beta_A - beta_B)  =  alpha_B - lambda^m alpha_A.

If `lambda^m beta_A - beta_B != 0` set `a = (alpha_B - lambda^m alpha_A)(lambda^m beta_A - beta_B)^{-1}`, check
`g^a = z` (at most 509 multiplications), and on success output `Q = [a]Q0`. The `lambda^m` alignment is the step a
bare equality of canonical keys would miss; it is done on stored seeds by replay, so only `(kappa, trail index)` is
stored.

The denominator vanishes with probability `1/r`: the trajectories depend on `(alpha_0, beta_0)` only through
`alpha_0 + a beta_0` and on the table only through `u_h + a v_h`, and `(alpha_0, beta_0)` and `(u_h, v_h)` are uniform
on `(Z/r)^2`, so `beta_A, beta_B` are independent and uniform given the trajectories.

### 4.5 Meter and termination

Counters: accepted steps `S`, multiplication slots `M` (accepted and rejected probes), starts, reports, recoveries.
The loop starts new trails while `S < S_max`; the trail in progress then runs to its end (at most `L_cap` further
slots, charged). The algorithm halts with failure if `M` exceeds the global halt cap, `M_max + L_cap` rounded up to
the next multiple of `2^64` (the rounding adds at most `2^64 - 1` slots, whose cost is below `2^64.1` units and is
inside the `2^92.361` additive total), starts exceed `ST_max`, reports
exceed `D_max`, or recoveries exceed `REC_max`. No loop is unbounded: every sampler is a fixed number of ChaCha20
blocks, and every trail has at most `L_cap` slots.

## 5. Work ledger

Unit: dense Fp12 multiplication = 1; other work = instructions / 20,000. (Exact arithmetic: `experiment.py`, mode
`ledger`.)

**Cost of one multiplication slot (accepted step or probe).** Every line is a hard upper bound on instructions
executed, itemised from the published prices (`tracks/docs/COST_MODEL.md`: one `blst_fp_mul` = 383 instructions,
one dense Fp12 multiplication = 22,996 = 1 unit).

| item | instructions |
| --- | ---: |
| dense Fp12 multiplication by a table entry (the reference operation) | counted as 1 unit |
| selector: 3 `blst_fp_mul` at 383 (section 3) | 1,149 |
| selector: 3 modular subtractions at 22 and a negation core at 13 (section 3), plus 9 for one three-`LDP` load of `p` per block (block (ii) shares one) | 88 |
| selector: the two extrema in one pass, `8 + 5*30`, and the eleventh comparison (6 `CSEL` + 1 `CSEL` + 1 `ADD`; section 3) | 172 |
| two `khash` words: 12 x (one limb `EOR` and the 8-instruction mix) over six limbs held in registers, plus the two accumulator loads and the two splitmix64 multiplicand loads (no callee-saved register holds either constant across a slot) | 112 |
| distinguished-point test: `AND x3, x2, #0xFFFF`, `ORR x3, x3, x1`, `CBNZ x3` — 64-bit registers, so 80 bits | 3 |
| look-back membership: `LSR x3, x1, #32`, `LDR x4, [sp, #1944]` (the count array base), `LDRB w5, [x4, x3]`, `CBZ w5, miss`, `B reject` — 5 on a hit, 4 on a miss | 5 |
| look-back update, counted ring write: old tag `LDR`, new tag `LSR`, count base `LDR`, eviction `LDRB`/`SUB`/`STRB`, insertion `LDRB`/`ADD`/`STRB`, tag store `STR w4, [sp, p, LSL #2]` | 10 |
| look-back bookkeeping: the ring slot's bounded update, the accepted-step counter with its cap test, two branches | 12 |
| state and call plumbing, table entry address: `LDR h`, `ADD h+a`, `AND (h+a) mod R`, `MOVZ 6912`, `MADD`, `LDR k_n`, `MOVZ 576`, `MADD`, `MOV x1, y`, `EOR x0, y, #0x400` | 10 |
| state and call plumbing, accept tail: `AND h', w1, #0x3FFF`, `STP h', k_n, [sp, #480]`, `MOVZ a, #0`, `EOR y, y, #0x400` | 4 |
| state and call plumbing, probe ladder on the reject path: `ADD a, #1`, `CMP a, #8`, `B.LT` | 3 |
| selector operand base, twice: `EOR x29, x20, #0x400`, after the dense kernel and after the three `fp_mul` kernels | 2 |
| the three `blst_fp_mul` argument setups, three `ADD` each | 9 |
| trail slot counter, cap `2^87` a multiple of `2^64`: one high-word `LDR`, `ADDS`/`ADCS`, `SUBS` high, branch | 5 |
| global slot counter, cap rounded up to a multiple of `2^64`: one high-word `LDR`, `ADDS`/`ADCS`, `SUBS` high, branch | 5 |
| next-slot branch | 1 |
| **subtotal** | **1,590** |
| less the probe ladder, which no accepted slot executes | -3 |
| **maximum over slot types** | **1,587** |
| explicit margin (9/1587 = 0.5671%, at least the 0.5% claimed) | 9 |
| **charged per slot** | **1,596** |

So each slot costs `1 + 1596/20000 = 21596/20000` units. The 22,996-instruction blst kernel already contains the
loads and stores of its operands, so table reads are not charged again. The probe ladder runs only on the reject
path, so no slot executes every row: the row sum is `1,590` and the maximum over slot types `1,587`; the
charge `1,596` is `9` above that maximum, `9/1587 = 0.5671%` (`experiment.py`, `ledger`:
`slot_margin_instructions`, `slot_instructions_max`), and 6 above the row sum. `1,623` is also the ceiling: at `1,624`
the least hundredth at or above `log2 W_max` becomes `125.82` (`experiment.py` finds it by exact integer search), so
the charge sits `27` instructions below it. The
sensitivity is `M_max/(20000 W ln 2) = 6.6*10^-5` bits per charged instruction, i.e. `6.68*10^-5` at the charged
`1,596`, and `7.2*10^-5 = 1/(20000 ln 2)` at the reference rate (the charge enters once, through the factor
`1 + c/20000`, so the second is an upper bound at every count; both from `experiment.py`, mode `ledger`).

**The instruction sequence of each row.** The modular
operations, the two extrema with their indices, the eleventh comparison, the distinguished-point test, the
membership test, the counted ring write, the seed and the per-trail clear+set are executed in turn on a register-file
machine
(`experiment.py`: class `A64`), which counts the instructions it runs and checks their result against a Python
reference; the three `blst_fp_mul` calls (1,149, the organizer's 383-instruction kernel price), the `khash` finaliser
(112), the bookkeeping row (12) and the plumbing rows (table entry address 10, accept tail 4, probe ladder 3, operand
bases 2, argument setups 9, the two slot counters 5 and 5, the closing branch 1) are priced from the written-out
sequences below and are not executed. A modular subtraction `x - y mod p` (three of them: `a - b` for `d`, `zeta*a - a`, `zeta*d - d`) is
three `LDP` of `x` into x0-x5, three `LDP` of `y` into x6-x11, the six-word `SUBS`/`SBCS` chain `r = x - y` in place
(whose final carry `C` is 1 exactly when `x >= y`), `B.CS` past the correction (taken exactly when that chain did not
borrow), and, on the borrow path, one `ADDS` and five `ADCS` of `p` from x12-x17, then three `STP`: 22 when the
correction runs and 16 when the branch is taken, so 22 bounds both paths; there is no call or argument setup because the
selector is one inline block that stores its results in its six-value buffer. The chain's carry propagates only
through `ADCS`; `ADC` would leave the carry of the first limb in place and produce a wrong residue, and
`x - y + (p if x < y)` is then the canonical residue in `[0, p)`. No borrow mask is formed: the branch replaces it,
so the subtraction writes no register outside x0-x17 and leaves `p`'s six registers intact for the next subtraction of
its block. `p`'s limbs are three `LDP` from a fixed frame offset, but no register holds them across a block: the kernels
clobber x0-x18, x19-x28 are the ten fixed live values, and the extrema occupy all of x0-x18 plus x29 and x30. So the
four operations sit in three blocks, each paying one three-`LDP` load, shared inside its block (allocation block
below): `a - b` before the three `blst_fp_mul` calls, since it is
an input to the third; `zeta*a - a` and `zeta*d - d` after them, on the kernels' outputs, sharing one load; and the
negation after the extrema. The negation `p - M` has the same shape with `XZR` as the minuend, its six-limb chain
in x12-x17, the same `B.CS` (taken when `M = 0`, leaving `0`) and no store: 13, `0` for `M = 0` and the canonical
`p - M` otherwise; with its load it is 16, so the row is
`3*22 + 13 + 9 = 88`, that is `25 + 47 + 16` (one subtraction with its load, two under one load, the negation). The
selector's eleven comparisons are itemised in section 3. Each `khash` word is six rounds of one `EOR` of a limb into
the accumulator and the splitmix64 finaliser (`LSR`, `EOR`, `MUL`, `LSR`, `EOR`, `MUL`, `LSR`, `EOR`: 8), so
`6*9 = 54` per word and 108 for both, with no limb loads because the limbs are the registers the selector's final
`CSEL` wrote (`x6-x11`); the two hoisted accumulators `mix(s_1)` and `mix(s_2)` are reloaded per slot (2) and so are
the two splitmix64 multiplicands, which no callee-saved register can hold across the dense kernel (2), so the row is
112.
The distinguished-point test is `AND x3, x2, #0xFFFF`, `ORR x3, x3, x1`, `CBNZ x3` (3), on 64-bit registers: `x1`
holds `w_1` and `x2` holds `w_2`, so the test reads all 64 bits of the first word and the low 16 of the second, 80
bits. Its 32-bit form would read `w_1 mod 2^32` and `w_2 mod 2^16`, a 48-bit predicate with `S_max*2^-48 = 2^77.7`
reports, far above `D_max = 2^47`; the executed case `w_1 = 2^32`, `w_2 = 0` is the witness that separates them.
The membership test is `LSR x3, x1, #32` (the tag of the candidate's first hash word), `LDR x4, [sp, #1944]` (the
count-array base), `LDRB w5, [x4, x3]`, `CBZ w5, miss`, `B reject`: 5 instructions, 4 on a miss. The base load is
inside the 12-bit `LDR` range; the byte index is the 32-bit tag and the array is `2^32` bytes. The test is
conservative with no second key compare: a positive count rejects the probe, and that is exact on every repeat
(a repeat's tag is stored, so its count is positive) and over-rejects only on a 32-bit collision. `ring_check`
executes the program and asserts 4 a miss and 5 the maximum (section 9).
The counted ring write is the old tag `LDR w3, [sp, p, LSL #2]`, the new tag `LSR x4, x1, #32`, the count-array base
`LDR x5, [sp, #1944]`, the eviction `LDRB w12, [x5, x3]`, `SUB w12, w12, #1`, `STRB w12, [x5, x3]`, the insertion
`LDRB w13, [x5, x4]`, `ADD w13, w13, #1`, `STRB w13, [x5, x4]`, and the tag store `STR w4, [sp, p, LSL #2]`: 10,
the list `write_prog` executes. The 32-bit register-offset form encodes `LSL #2`, so no address `ADD` is charged, and
`LDR w3` zero-extends the tag into x3. The tag is the top 32 bits of
the word the index and the distinguished-point test already have; the count array carries the duplicate-tag count.
The bookkeeping is
`ADDS`/`ADCS` of the 126-bit
accepted-step counter, `SUBS`/`SBCS` against `S_max` and a branch (5), with the two `LDR` of `S_max` from its frame
offset (7), the slot update `ADD p, p, #1`, `CMP p, #24`, `CSEL p, p, xzr, NE` (3: three register operands, taken
when `Z = 0`, no immediate), and two branches: 12. The plumbing rows are the table's lists, in the loop's order, written out with registers and
frame offsets in `experiment.py` (`plumbing_texts`): the dense kernel is called as `blst_fp12_mul(x0, x1, x2)`
because it clobbers x0-x18; `h` and `k_n` are adjacent frame words, so one `STP` stores both; x29 is scratch, so the
operand base is formed twice; both slot caps are multiples of `2^64`, so the low words cannot carry past them and one
high-word compare suffices. The accepted-step test and the update rows (`7 + 12`) are charged at every slot, over-charging rejected probes.

**The look-back window's structure.** The window holds the 32-bit tags of the last at most 24 states of the trail,
including the current one; membership is all the escape rule of section 4.3 uses, and the oldest tag is dropped on
every accept. It is seeded with `tag(kappa_0)` and stored as a `2^32`-byte count array plus a 24-slot ring of tags:

- `T[0..23]`: the tags, one 32-bit word per slot at `sp + 4*slot` (96 bytes), so the tag store is one 32-bit
  register-offset `STR`. The k-th state's tag goes to slot `k mod 24`, the
  slot of the oldest stored tag, which one register `p` names; `p` advances by the bookkeeping row's bounded
  `ADD`/`CMP`/`CSEL` update.
- `C[0..2^32-1]`: one byte per tag, `C[t]` = the number of slots holding `t` (0..24). The array is cleared to zero
  once before the first seed (`2^32` byte stores, charged in preprocessing); the first seed sets the 24 slots' count
  to 24; a per-trail clear+set decrements the previous trail's 24 tags and sets the new trail's tag to 24.

The first seed puts `tag(kappa_0)` in all 24 tags with count 24, `p = 0`, 28 instructions (the tag `LSR`, the count
base `LDR`, the `MOVZ`/`STRB` count, and the 24 tag `STR`s), and the per-trail clear+set is 124 (one base `LDR`,
24 x old-tag `LDR`/`LDRB`/`SUB`/`STRB`, `LSR`/`MOVZ`/`STRB` for the new count, and 24 tag `STR`s), inside the start
item's slack. Seeding every slot with `tag(kappa_0)` keeps the ring
exact: after `k` accepts the stored tag multiset is `{tag(kappa_0), ..., tag(kappa_k)}`, the last `min(k+1, 24)`
states, with no stale tag and no empty-slot sentinel (the count array is cleared by the previous trail's clear+set).
Once `k >= 23` the accept overwrites `tag(kappa_{k-23})`, and a closed pattern of length `l <= 24` returns to a tag
still in the window, which is what section 7 uses.

Membership of a candidate `kappa'` is `C[tag(kappa')] > 0`: 5 instructions (4 on a miss). The test is
**conservative** with a hard bound: every repeat is caught, because a repeat has a stored key's tag, so its count is
positive and it rejects, and no pattern of section 7 is missed. A rejection can also happen with no stored key equal
to `kappa'`, when its 32-bit tag is in the window: under H2 (section 6) that tag is uniform and independent, so at
most `24*2^-32` per candidate and at most `8*24*2^-32 = 2^-24.42` per accepted step. Such a *tag escape* costs
exactly what an ordinary escape costs: the probe is rejected and its multiplication charged; it cannot lose a repeat
and cannot change the stored tags. The escape rate is therefore at most `nu = 1.1/(12R) + 8*24*2^-32 = 5.6395928e-6 <
2^-17.43` (an allowance, not a measured rate), the tag term being a relative 0.8% of `nu`, added to the `1.1/(12R)` term and not absorbed in the `mu` allowance; the rate measured in section 9 is measured on
the exact-key rule, a lower bound on this one.

**The register allocation of the whole loop.** A per-slot bound is a bound only if the sequences fit the 31 general
registers or pay for their spills, so the map below is complete: every value live across a call, every frame word
loaded per slot and every scratch the phases need, with the row that pays for each load or store.

| register | holds | live across |
| --- | --- | --- |
| x19 | `&TW`, the twist table base | whole loop, every kernel call |
| x20 | the current y buffer pointer; the other buffer is `x20 ^ 0x400`, and the frame is 2 KiB aligned so the toggle is one `EOR` | whole loop |
| x21, x22 | the global slot counter `M`, low and high word | whole loop |
| x23, x24 | the trail slot counter `c`, low and high word | whole loop |
| x25, x26 | the accepted-step counter `S`, low and high word | whole loop |
| x27 | the ring slot `p` | whole loop |
| x28 | the probe counter `a` (0..7), reset on an accept | per state |
| x29 | the minimum's index through the selector tail, then the y' base | per block |
| x30 | the maximum's index through the selector tail, then `i_M + 6`; the chains' `xzr` writes need no temp | per block |
| x0-x18 | kernel argument registers and per-phase scratch; after the selector `kappa` is in x6-x11, `w_1` in x1 and `w_2` in x2 | per phase |
| sp+0 | `T[24]`, one 32-bit tag per 4-byte slot (96 bytes) | whole trail |
| sp+1944 | the count-array base (the array itself is the `2^32`-byte heap region) | whole trail |
| sp+96 | `V[6]`, the six stored selector values (`V_OFF`; the extrema scans it) | per state |
| sp+384 | the modulus `p`'s six limbs (`P_OFF`); `p` the *ring slot* is x27, not a frame word | per block |
| sp+1680, 1728 | the tower constants `zeta` and `g0` | whole trail |
| sp+480 | `h`, the current state's table index (`LDR` in the head row) | per state |
| sp+488 | `k_n`, the current state's twist index (`LDR` in the head row, `STR` in the accept tail) | per state |
| sp+1808, 1824 | the two `mix(salt)` accumulator seeds and the two splitmix64 multiplicands | whole trail |
| sp+1840, 1856, 1872 | the three 126-bit constants `L_cap`, the rounded-up global halt cap and `S_max` | whole trail |
| sp+432 | the `a - b` scratch the third `blst_fp_mul` reads | per state |
| sp+2048, sp+3072 | the two 576-byte y buffers, `0x400` apart with bit 10 clear in the lower one | whole trail |

Every `LDP`/`STP` immediate is at most 480 from sp (`V`, `p`, `a - b`, `h`/`k_n`) and 368 from the y' base (`c1`'s
`a`, `b` at +288 and +336, blst layout), inside the scaled 7-bit range -512..504; every other word is reached by a
scaled 12-bit `LDR`/`STR` or an `ADD` immediate, so no row needs an address-forming instruction.

| phase | live scratch values | registers | loaded or stored per slot |
| --- | --- | --- | --- |
| table entry address and the call's arguments | `h`, `a`, `k_n`, two strides, the address, three pointers | 8 of x0-x18 | `h` (1), `k_n` (1), the two strides (2), `MOV x1, y` (1) |
| `blst_fp12_mul` call | three argument pointers | 3 of x0-x18 | none (x19-x28 survive a call) |
| subtraction blocks (i), (ii): `a - b`, then the two `zeta^2` terms | `x` (6) in x0-x5, `y` (6) in x6-x11, `p` (6) in x12-x17, the y' base in x29; no borrow register (the `B.CS` replaces the mask) | 18 of x0-x17, x29 untouched | `p`'s six limbs by three `LDP`, 3 in (i) and 3 shared by (ii)'s two subtractions |
| subtraction block (iii), inside the tail: the negation `p - M` | `M` (6) in x12-x17, `p` (6) in x0-x5; the tail's key, indices and index supply stay live | 18 of x0-x17 | `p`'s six limbs by three `LDP` |
| `blst_fp_mul` argument setups | three pointers per call, three calls | 3 of x0-x18 | none: three `ADD` per call (9) |
| extrema network, seed | minimum (6) + index (1), maximum (6) + index (1), candidate (6), the y' base | 18 + x29 + x30 | `V_0`'s six limbs are loaded twice by `LDP` from the y' base, `V_1`-`V_5` from their frame words |
| extrema network, loop | minimum (7), maximum (7), candidate (6), index supply | 21 = x0-x18 + x29 + x30 | none |
| `khash` | `kappa` in x6-x11 (6), `w_1` in x1, `w_2` in x2, one temp, two multiplicands | 11 of x0-x18 | the two seeds (2) and the two multiplicands (2) |
| membership test | `w_1` in x1, the tag x3, the count base x4, the count w5 | 4 of x0-x18 | the count-array base (1) |
| counted ring write | `w_1` in x1, the old/new tags x3/x4, the count base x5, the counts w12/w13, the slot x27 | 6 of x0-x18 | the base (1) and the tag store (1) |

**The subtraction's own registers, executed.** The three subtractions write only x0-x17 and the negation only
x12-x17 while reading the six `p` registers x0-x5; none writes x19-x28 or x29, so the ten fixed live values and the
y' base survive, and `B.CS` replaces the borrow mask so no extra register is needed. `experiment.py`
(`modsub_block_check`) seeds x19 = `&TW`, x20 = the y buffer, x21/x22 = `M`, x23/x24 = `c`, x25/x26 = `S`, x27 = the
ring slot, x28 = the probe counter and x29 = the y' base, and asserts after every subtraction that no register
outside x0-x17 changed, that x12-x17 still hold `p`, and that only the six destination frame words were written.
Block (ii) is executed as the composition it is priced as: two successive subtractions under ONE load of `p`, 47
instructions at most for the pair, both results checked against the Python reference on the printed tuples. The
rejected shape is the negative control: the same check catches `CSETM x19` and a load of `p` into x20-x25. Both branch
paths are exercised (176 tuples, 68 correcting and 108 skipping), and a block (i) tuple costs 25 with its load when
the correction runs and 19 when it does not.

The 126-bit constants (`L_cap`, the rounded-up global cap, `M_max + L_cap`, `S_max`), `p`'s six limbs, `h` and `k_n`
have no callee-saved register: the ten fixed values fill x19-x28 and the kernels clobber x0-x18, so they live in frame
words and the rows above pay for every load: the `2 LDR` inside the accepted-step row's 7, one high-word `LDR` in each
counter row's 5, the 9 in the selector's row, and one `LDR` each for `h` and `k_n`, whose `STP` is in the accept
tail. The only state that
survives a rejected probe and is *not* in a register is therefore those two words, and both are charged on every path.
Nothing else is spilled: the kernel writes the candidate into the other y buffer and x20 is toggled only on an accept.

**The two hash words.** The index and the distinguished-point predicate are both functions of the class key
(section 4.1): `w_1 = khash(kappa, s_1)`, `w_2 = khash(kappa, s_2)`, the index `w_1 mod 2^14`, and the window tag
`w_1 >> 32`, so no third word is computed. A state is distinguished iff `w_1 = 0` and `w_2 mod 2^16 = 0`, 80 bits. Distinguished states are reported
and never stepped from, so each value of the low 14 bits of a stepped-from state has probability
`2^-14 (1 +- 2^-50)`. H2 is stated on the two words themselves, and both are charged at every slot, which over-charges
the `2^-17.43` of rejected probes.

**Main walk.** `N = (r-1)/12`, `mu = 11/10`, `x = 3903/10000`, `L = -ln(1-x) = ln(10000/6097)`, and

    U = sum_{k=1}^{64} x^k/k + x^65/(65(1-x)),

an exact rational. The partial sum is below `L` and the geometric tail dominates the terms it drops, so `U` is an upper
bound for `L`, and `0 < U - L < 2^-99`: the computed gap is `6.63e-31`, and the bound is itself a rational comparison,
because `U` minus the partial sum through 70 terms is below `2^-99` while that partial sum is below `L`.

The accepted-step budget is computed from `U`, the larger rational bound. For a target main-event probability
`3903/10000 = 0.3903` (the rule's `0.39` plus the `2.71*10^-4` loss bound of section 6) it is the least integer with
`S^2 >= 2 mu N U`, checked exactly in integers:

    S_max = ceil(sqrt(2 mu N U)) = 68967528808679891282472564154805017972   (log2 = 125.69726)
    M_max = S_max + floor(S_max/2^13) + 1 = 68975947696473919589513881606484071319,  work 21596/20000 * M_max.

No logarithm, float or decimal enters `S_max`: `U` is rational, the ceiling is exact, and `S_max^2 >= 2 mu N U >
(S_max-1)^2` is an integer comparison. The `L`-based ceiling is a different integer: `U/L - 1` is `2^-99.24`, so
`ceil(sqrt(2 mu N L))` is about `2^25.46` steps lower, and taking `U` costs at most a relative `2^-99`. Every
inequality here, in section 6 and in the claim is stated on this `S_max`, the integer `experiment.py` prints, and the
exact checks `W^100 <= 2^12581` and `W^100 > 2^12580` are formed from it.

The probe allowance `S_max/2^13` is `2^-13` of the accepted steps; the escape (probe) rate is at most
`nu = 1.1/(12R) + 8*24*2^-32 < 2^-17.43` per step by `L2` below (measured `5.57e-6` on the exact-key rule at the
production table size in section 9). The number of rejected probes is therefore dominated by a binomial with mean
`S_max * nu < 2^108.26`, and the
allowance `2^112.7` exceeds that mean 21-fold, so exceeding it needs a deviation whose binomial tail is below
`exp(-2^107)`; the meter fails on such an excess and `L7` carries this term. The allowance is a hard cap on the
charged work, not an expectation.

**Additive items** (log2 of each, fixed by caps, so hard bounds):

| item | amount (units) | log2 |
| --- | --- | ---: |
| two pairings (`g` and one spare) at 2^30 instructions each | `2^31/20000` | 16.72 |
| table: `R * (2*509 + 1 + 11 + slack) <= R * 1112` | | 24.12 |
| tag count array zeroing: `2^32/20000` byte stores | `2^32/20000` | 17.71 |
| starts: `ST_max * 1100` (1,019 multiplications + ChaCha20 and two reductions < 0.4 unit) | | 57.11 |
| reports: `D_max * 2^16 / 20000` (red-black tree and transfer; a tree operation is at most 2,300 instructions) | | 48.72 |
| replays: `REC_max * 2 * L_cap * 5/4` (slot 21596/20000 plus exponent update < 0.075 from two 255-bit products and two additions, at most 1,500 instructions, `(u_h, v_h)` from a stored 2^20-byte array, no PRG call; `1.15825 <= 5/4`) | | 92.33 |
| completing the last trail: `L_cap * 21596/20000` | | 87.13 |
| `[a]Q0` and `g^a == z` | `<= 2^12` | 12.00 |
| **additive total** | | **92.361** |

**Total.** `W_max = 21596/20000 * M_max + 2^92.361 <= 2^125.81` (`log2 W_max = 125.80819701...`); the claim
`2^125.81` is the least hundredth at or above it, and `experiment.py` (mode `ledger`) proves that with exact integers
alone: `m` is the least integer with `W^100 <= 2^m`, found by binary search on `m` over `[0, 100 * bitlength(W)]`, so
`W^100 <= 2^12581` together with `W^100 > 2^12580` (no smaller hundredth suffices). No float enters that step.
Randomness is charged inside starts/table (ChaCha20), sorting/lookup/transfer in the report item, and verification
(`g^a == z`) and output are charged. The attack does not amortize over targets.

**Advice:** none; the preprocessing is charged in section 8.

## 6. Success probability

*Model (H1).* Let `X_1, X_2, ...` be the classes of the accepted states of successive trails. H1: each fresh state's
successor behaves like a uniform class, up to a walk-quality loss `mu`: the probability that no class has been
visited by two different trails after `s` accepted steps is at most `exp(-s^2/(2 mu N))`. For an ideal random map
`mu = 1`, and then `E[S^2/(2N)] = 1` for the first repeat `S`; section 9 measures `mu_hat = 0.996 +- 0.014`
(R >= 256) against `mu = 11/10`. With `mu = 11/10` and `S = S_max`, `S_max^2/(2 mu N) >= U > ln(10000/6097)`, so the
right side is at most `6097/10000 = 1 - 3903/10000`, and
**the main event** (some class is visited by two different trails within `S_max` steps) has probability at
least `3903/10000 = 0.3903`.

*From a repeat to a recovery.* Trails run sequentially, so every earlier trail has ended. Consider the first class that
a later trail visits and that an earlier trail visited. If that earlier trail ended with a report (it did, except with
the probability bounded in L1/L3), the later trail is on its path to its distinguished point. By coalescence (section
4.3) it reaches the same distinguished point, and the store holds the key from the earlier trail; Recover runs and
succeeds unless the denominator vanishes (`1/r`).

*Losses, each an upper bound on the probability that the main event fails to yield a recovery.*

| label | cause | bound |
| --- | --- | ---: |
| L1 | one of the two trails is cut by `L_cap` (no distinguished point in `2^87` slots: probability `e^-128 = 2^-184.7` for the later trail, and the length-biased `L_cap*theta*e^-L_cap*theta = 128*e^-128 = 2^-177.66` for the earlier trail it hits, which is selected by being hit) | `<= 2^-176.6` |
| L2 | an escape of one of two meeting trails inside the 24 steps after meeting (`nu = 1.1/(12R) + 8*24*2^-32 < 2^-17.43` includes the tag rule's `2^-24.42` collision term of section 5) | `<= 2*24*nu <= 2.71e-4` |
| L3 | one of those two trails is trapped in a residual closed pattern, loops in a within-trail rho cycle, or has all 8 probes rejected (per trail: `2^80 * 2^-163.5`, `E[L^2]/(2N) = 2^-91`, `2^80 * nu^8 <= 2^-59.49`) | `<= 2^-58` |
| L4 | ties among the 12 values `V_k` (8 linear conditions on `(a, b)`: `a = 0`, `d = 0` and `d = zeta^t a` for six `t`, so probability at most `8/p < 2^-377.6` per state) or `c1 = 0`; a state equal to the identity, `1/r < 2^-254.84` per slot (section 4.3); over `2^126` slots | `<= 2^-128.84` |
| L5 | false store match or zero denominator | `<= 2^-250` |
| L6 | sampling distance, `2^47` starts and `2^14` table scalars at `2^-129` | `<= 2^-80` |
| L7 | caps `D_max, ST_max, M_max` exceeded (Poisson tail of a count with mean `2^45.7`, `2^46` at most) | `<= 2^-1000` |
| L8 | more than `REC_max` unusable matches | `<= 2^-200` |

Total loss below `2.71*10^-4` (`L2` is the only label above `10^-17`). Hence, under H1-H3, success probability at
least `0.3903 - 0.000271 = 0.390029 > 0.39`, a printed margin of `0.000029` over the required `0.39`.

*What H1 must supply, and the margin.* The claim uses H1 only at the single point `s = S_max`: the probability that
the first repeat occurs by `S_max` is at least `3903/10000` (equivalently `S^2/(2N) <= mu U` with probability
`>= 0.3903` for the first repeat `S`). Success `>= 0.39` after the losses needs `0.390271`, so the true inefficiency
must satisfy `mu_true <= 1.1 * ln(1/0.6097)/ln(1/0.60973) = 1.1 * 0.494788/0.494739 = 1.10010`: the target sits at
the rule's requirement, so the claim spends the whole `mu = 11/10` allowance and keeps `0.01%` of spare beyond it.
The measured `mu_hat = 0.996 +- 0.014` (section 9) is `7.4` standard errors below the allowance; for an ideal random
map the value at `S_max` is `1 - (6097/10000)^{1.1} = 0.4197`, and section 9 reports that statistic itself, not only
the mean.

**Uniform distinguished-point and routing fields (H2).** The only structure used from the class key is: the two hash
words `w_1, w_2` (hence `hidx`, the distinguished-point predicate and, since section 5 stores `w_1 >> 32` in the
window, the window's tags with their `24*2^-32` collision term) and `k*` behave as independent uniform fields
for distinct classes. For `k*` this is
exact for a uniform GT element (the map `y -> sigma(y)` is a bijection that shifts `k*` by one, so `k*` is uniform and
independent of the class). For `khash` it is the usual hash heuristic on `(w_1, w_2)`: the two words are
independent uniform 64-bit values, so `hidx = w_1 mod 2^14` is uniform on `R` and the predicate has density `2^-80`.
The one correlation this creates, a distinguished state having index 0, is outside every step the walk takes, since a
distinguished state is reported and never stepped from (section 5). The algebraic checks of section 3 are the only
production-field evidence. A distinguished density within a factor 2 of `2^-80` changes only report counts (cap
`D_max = 2^47` against a mean `2^45.7`, so up to `2.3 * 2^-80` still fits).

**H3.** ChaCha20 under a fresh key is an ideal source of coins.

## 7. Deterministic traps and the look-back window

A *trap* is a deterministic closed pattern: the class path returns to a class after `l` accepted steps. In the lifted
exponents (`x_n = alpha_n + a beta_n`, `m_h = u_h + a v_h`; the accepted step `i` uses table index `u_i` and twist
`k_i`; `u_i = h_i + a_i` where `h_i = hidx(kappa_i)` and `a_i` is the number of rejected probes) a return to the class of
`x_t` after `l` steps means `lambda^s x_t = x_t + sum_{i<l} lambda^{-k_i} m_{u_i}` for some `s` in `Z/12`.
If `s != 0` and the sum vanishes formally then `x_t (lambda^s - 1) = 0`, i.e. `x_t = 0`, excluded; if it does not
vanish formally, see below. Probes are included: the lemma speaks about the *used* indices `u_i`, and a probe only
changes `u_i` from `h_i` to `h_i + a_i`.

*Formal part.* The `m_h` are independent uniform scalars. Put `P_h = sum_{i: h_i = h} zeta_12^{-k_i}`, an element of
`Z[zeta_12]` with coefficients of absolute value at most `l`. The ring map `zeta_12 -> lambda (mod r)` has as kernel a
prime ideal of norm `r` (`r = 1 mod 12` splits), so a nonzero element `alpha` in the kernel has `r | N(alpha)`; with
coefficients bounded by `B` in the basis `1, zeta, zeta^2, zeta^3`, `|N(alpha)| <= (4B)^4`, which forces
`B >= (r/256)^{1/4} > 2^61`. Hence for `l < 2^60` a nonzero `P_h` has a nonzero, invertible image modulo `r`.
Then `sum_h P_h(lambda) m_h` is uniform modulo `r` as soon as one `P_h != 0`, and the return condition is one linear
equation with a uniform solution set of probability `1/r` for each pattern: an ordinary random revisit, which H1
covers (a window revisit of this kind has probability about `24^2/N < 2^-240`; it is not a structure of the walk).
For `l >= 2^60` the necessary condition below is all that is used; non-formal returns of such length are ordinary
within-trail rho collisions, probability `<= L_cap^2/(2N) <= 2^-78` per trail at the cap and `E[L^2]/(2N) = 2^-91`
for the real trail-length distribution (counted in L3).
So a *deterministic* trap requires `P_h = 0` **for every** `h` used. A single root of unity is nonzero, so
**every index used in a trap is used at least twice**. (The minimal vanishing sums of 12th roots of unity are
antipodal pairs, giving the two-cycle with `k` differing by 6, and triples at 120 degrees; the argument does not
need this classification.)

*Probability.* Given the ordered used indices of `l` accepted steps (each uniform on `R` by H2 for distinct classes:
`u_i` differs from `h_i` only on a rejected probe, probability `<= 1.1/(12R) < 2^-17`, which multiplies the bound by
at most `(1 + 2^-12)^l`, negligible), let `D` be the number of distinct indices. A trap needs every fibre of size
`>= 2` (so `D <= l/2`) and, for each used index, its twists to sum to zero. Conditional on the other twists of a fibre
the last must take one specific value, so each fibre vanishes with probability at most `1/12`, independently (the
twists are partial sums of uniform i.i.d. shifts). Thus

    P(trap of length l starting at a given step) <= E(l) = sum_{d=1}^{min(floor(l/2), R)} C(R,d) (d/R)^l 12^-d .

`experiment.py` (mode `ledger`) evaluates this sum in log scale with float logarithms, so the printed `2^-164.2` and
`2^-163.5` are float evaluations of the exact integer sum; the only exact comparison is the final claim inequality.
The maximum over `25 <= l <= 512` is `2^-164.2` (at `l = 25`) and the sum is `2^-163.5`. For `l > 512`:
`C(R,d) <= (eR/d)^d` makes each term at most `((e/12)(d/R)^{l/d-1})^d`; for `d <= R/2` that is at most
`0.2266^d 2^{-(l-d)} <= 2^-l`, and for `d > R/2` at most `0.2266^{R/2} < 2^-17000`. So `E(l) <= 2^-480` for all
`l > 512` and their sum to the trail cap is below `2^-400`, negligible next to `2^-163.5`.
For lengths `l <= 24` a trap is detected within `l` steps by the window and escaped (the hash rule of section 5
rejects every repeat of a window key, so no such trap is missed); the unconditional escape
probability per step is at most `nu = 1.1/(12R) + 8*24*2^-32` (dominated by the two-cycle: the new class has the same
index as the previous one, probability `1/R`, with twists differing by 6, probability `1/12`; the tag rule adds
`8*24*2^-32 = 2^-24.42`, which is a relative 0.8% of the bound), measured `5.57e-6` on the exact-key rule against
the bound `5.64e-6` and `1/(12R) = 5.09e-6` in section 9. Longer escapes (lengths 3 to 24) are of order `R^-2` or smaller and are inside the
factor 1.1.

*Cost of undetected traps.* A trail of mean length `2^80` is trapped with probability at most `2^80 * sum_l E(l) =
2^-83.5`, and then runs `L_cap = 2^87` slots (128 mean lengths), so the expected wasted fraction of the main walk is at
most `L_cap * sum_l E(l) <= 2^87 * 2^-163.5 = 2^-76.5`. A trail dies (all `A = 8` probes rejected) with probability at
most `nu^8 <= 2^-139.49` per step (the 8 candidate classes use 8 different indices; the two-cycle event can hold for at
most one of them, so the product bound is far from tight, and independence is the heuristic H2); per trail
`2^80 * nu^8 <= 2^-59.49`.

## 8. Memory, preprocessing, advice

| item | bytes (bound) |
| --- | --- |
| store: `<= D_max = 2^47` nodes of 96 bytes (48-byte key, 16-byte trail index, three pointers, colour) | `2^53.585` |
| twist table `12 * 2^14 * 576` | `2^26.75` |
| scalars `(u_h, v_h)` for replay: `2^14 * 64` | `2^20` |
| tag count array: `2^32` bytes | `2^32` |
| trail state: `y` and the probe buffer (2 x 576 B), tag ring (24 x 4 B = 96 B), the selector and hash buffers, counters, ChaCha20 state | `< 2^20` |
| **peak** | `<= 2^53.6 < 2^54` bytes (`2^53.585 + 2^26.75 + 2^32 + 2^20 + 2^20`) |

All of it lives in one address space below `2^64`; every store transfer is inside the 2^16-instruction report cost.
Preprocessing is the table, the `2^32`-byte zeroing and one pairing, `2^24.15` units. There is no advice.

## 9. Experiments and their scope

One declared experiment (`hashring-rho-checks`, the verifier allows one initial run) executes four parts in one
process (`{"mode":"all", ...}`, about 10 s on a laptop; the parts are also callable alone by `mode`). None is a
production-scale run; none measures hardware or the production success probability, and the counts below are those
of the written-out sequences, executed in Python.

* **ledger** recomputes sections 5, 7 and 8 in exact integer and `Decimal` arithmetic where the quantity is exact
  (`S_max`, `M_max`, each additive item, `log2 W_max = 125.80819701 <= 125.81`, with the exact integer checks
  `W^100 <= 2^12581` and `W^100 > 2^12580`, the maximum slot `1,587`, the charge `1,596`, the ceiling `1,623`,
  `2^53.585` bytes, `2^24.144` preprocessing), the trap bounds in float log scale as stated in section 7, and it
  re-checks the one-pass min-max network against two separate scans, and its result against the direct twelve-value
  scan, on 117,649 six-tuples (ties and zeros included). `selector_min_check` and `ring_check` then execute the priced
  AArch64 lists on the register-file machine `A64`, which counts what it runs and checks each result against a Python
  reference: the selector tail (188 instructions per tuple, 182 when the maximum is `0` and the negation skips, on 245
  small-field tuples — 51 ties and the all-zero one — and 40 production tuples), the modular subtraction in all three
  of its blocks (22 correcting, 16 when `x >= y`, 25 and 47 with the block's load of `p`), the composition of two
  successive subtractions under ONE load of `p` (47 at most), with x19-x28 and x29 seeded from the allocation and
  asserted byte-identical after every subtraction, `p`'s six registers asserted unchanged and the negative control
  (the mask shape over x19, x20-x25) caught, the distinguished-point test (3, with the `w_1 = 2^32` witness that
  separates 80 bits from 48), the slot update for 200 accepts (3 each, `p = k mod 24`, always the
  oldest key's slot), the first seed (28), the per-trail clear+set (124), the counted ring write (10) and the
  membership test (4 a miss, 5 a hit), over every fill 0..24, 399 candidates checked against the window's tag set
  and the count array asserted against the tag multiset after every accept. The plumbing rows of section 5
  (10, 4, 3, 2, 9, 5, 5, 1), the bookkeeping row (12), the 1,149-instruction kernel price and the 112-instruction
  `khash` finaliser are priced, not executed, and section 5 says so for each; `encoding_check` assembles those rows'
  written-out lists (each length asserted equal to its row) with every executed list, 599 instructions, asserts each
  immediate encodable for its form (largest `LDP`/`STP` offset 480) and refuses three controls (`LDP` at sp+1632,
  and the two mis-scaled tag stores).
* **algebra** (real BLS12-381): `lambda` of order 12 with `lambda^6 = -1`; `zeta` a primitive 6th root of unity in Fp;
  `sigma(y) = y^p = y^lambda` on random GT elements; the tower constants satisfy `g1 = -g0` (`gamma = g0*(1-u)`),
  `g0^6 = -1/8` and `zeta = 2 g0^2`; the selector's 12 values are the true first coordinates, and the
  3-multiplication form gives the same 12 values as the general `g0*a + g1*b` form; the priced two-half minimum
  returns the same `(kappa, k*)` as the direct twelve-value scan, on the element and all twelve conjugates; `kappa`
  invariant, `k*` covariant; the twisted step is class-equivariant. Exact checks of section 3 on production data
  (2 elements by default).
* **toy**: the complete algorithm (selector, twisted table, window and probing, distinguished points, store, replay
  with exponent tracking, Frobenius alignment, `g^a == z`) over a toy tower of the production shape:
  `q` of 60 bits with `q mod r' = lambda'` of order 12 and `q = 7 mod 12`, `1+u` neither square nor cube in Fq2, a
  prime subgroup order `r'` of 22 bits, real Fq12 arithmetic, a table of 64 multipliers, window 12, trails of mean
  length 8. Default run (6 instances, seed 3): 6 of 6 recovered logarithms equal the planted one, 0 wrong outputs, the
  tracked replay reproduced every untracked path (same endpoint key), 5 escapes exercised, 0 dead trails, mean
  accepted steps `1.578 sqrt(N')`. Offline run (`{"mode":"toy","bits":22,"trials":80,"seed":77}`): 80 of 80 correct,
  0 wrong, replay matched, 41 of 80 (51.3%) inside the toy budget (ideal map `41.97%`, requirement 39%), mean accepted
  steps `1.219 sqrt(N')`. Every tower is built with `gamma = g0*(1-u)`, so the toy runs the same 3-multiplication
  selector as production, and it checks that selector's 12 values against the true first coordinates once per
  instance. Its window is the production tag ring and count array itself (the same `w_1` the index uses, the same
  `ADD`/`CMP`/`CSEL` update, the same priced 5-instruction count test, the same 10-instruction counted write and the
  same 28/124 seed and clear), with every membership decision cross-checked against the window's tag set.
* **sim**: the identical walk in the exponent domain of a cyclic group of order `r'` with automorphism `lambda'` of
  order 12 (state = discrete logarithm; keys are 64-bit hashes of the 12 conjugates), with `R = 2^14`, window 24 and
  probing 8; short trails keep within-trail rho cycles negligible at toy sizes. It measures `mu_hat = mean(S^2/(2N))`
  at the first repeat of a class by a different trail, the same estimator on a plain-birthday control, the escape rate,
  and whether the meeting trail reaches the same distinguished point, and the fraction of groups whose
  first repeat satisfies `S^2/(2N) <= 1.1 U = 0.54427` (inside the `S_max` budget). Default run
  (600 groups of 22 bits, trail mean 8, seed 3): `mu_hat = 1.028 +- 0.041`, fraction within budget `0.4167`, control
  `1.058`, 599 of 600 merged trails reach the same distinguished point. The observed escape rate there is `4.4e-5`,
  above `1/(12R) = 5.1e-6` because at `N = 2^18` the window also catches random revisits (about `24/N` per step); at
  28 to 34 bits that term is below `10^-6` and the rate matches `1/(12R)`.

Offline runs of the same code with larger inputs (`python3 -I experiment.py`, `{"mode":"sim",...}`):

| bits of r' | R | trials | seed | `mu_hat` | s.e. | fraction within `S_max` budget | control | escape rate / predicted `1/(12R)` | merged trails reaching the same DP |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| 28 | 64 | 1500 | 105 | 1.072 | 0.026 | 0.4053 | 1.085 | 1.26e-3 / 1.30e-3 | 1493 / 1500 |
| 30 | 256 | 1500 | 104 | 0.992 | 0.026 | 0.4073 | 1.028 | 3.10e-4 / 3.26e-4 | 1500 / 1500 |
| 30 | 1024 | 1500 | 103 | 0.999 | 0.026 | 0.4173 | 0.935 | 7.68e-5 / 8.14e-5 | 1500 / 1500 |
| 30 | 16384 | 1500 | 101 | 1.008 | 0.026 | 0.4180 | 0.968 | 5.57e-6 / 5.09e-6 | 1500 / 1500 |
| 32 | 16384 | 700 | 102 | 0.969 | 0.038 | 0.4386 | 1.041 | 4.43e-6 / 5.09e-6 | 700 / 700 |
| 34 | 16384 | 250 | 106 | 0.995 | 0.063 | 0.4320 | 1.005 | 5.19e-6 / 5.09e-6 | 250 / 250 |

Pooled over the 5450 repeats with `R >= 256` (the five bottom rows) `mu_hat = 0.996 +- 0.014`, and the pooled fraction
within the budget is `0.4182 +- 0.0067` (ideal map `0.4197`, requirement `0.3903`, 4.2 standard errors of margin). The
control column is the same estimator on a plain birthday over `N` classes, `min(trials, 300)` control trials per row,
standard error `0.058` (`0.063` last row), ranging over `0.935 ... 1.085`. The claim uses `mu = 11/10`, a 10% allowance
over the measured `0.996` (7.4 standard errors). With long trails relative to `sqrt(N')` the estimator is inflated by
within-trail rho cycles, as it must be: `experiment.sim(26, 16384, 300, 5, db=8)` (trail mean `2^8`) gives
`mu_hat = 1.887` (294 of 300 merged trails reached the same distinguished point), decreasing like `1/N'`; that effect
is the within-trail self-collision probability `E[L^2]/(2N) = 2^-91` at production sizes (`L = 2^80`, `N = 2^251`), so
the rows above use short trails (mean `2^5`).

**Extrapolation to the production size** is the principal uncertainty of H1: the emulation replaces the Fp12 structure
by pseudorandom keys, and `N = 2^251.3` is far outside any run. Both emulations draw the index and the
distinguished-point bits from independent salts, the idealised case of the production rule (index = `w_1 mod 2^14`,
`64` of the `80` bits `w_1 = 0`). The algebra checks are exact; the dynamics are a disclosed heuristic.

## 10. Limitations and disclosed heuristics

* H1, H2, H3 as above; success `>= 0.39` is conditional on them.
* The walk-quality allowance `mu = 11/10` is an empirical margin, not a proof: the measured `mu_hat = 0.996 +- 0.014`
  would support a smaller one. With the budget set at the rule's own requirement (main event `0.3903`, success
  `0.390029` after the losses), the claim survives a true inefficiency up to `1.10010`.
* The store of `2^47` records (`2^53.6` bytes) is a mathematical object here, inside one 64-bit address space.
* Cost prices are instruction estimates on the organizer's pinned kernel prices; the selector, hashing and look-back
  rows are upper bounds, not measurements, carrying a 0.5671% allowance (1,596 charged against the 1,587 maximum
  slot). The selector's third multiplication rests on an exact identity of the tower constant `gamma` (section 3).
* The window's rejection rule over-rejects on a 32-bit tag collision, a charged escape that cannot lose a repeat
  (section 5). The claim survives 27 further instructions per slot above the charge, 33 above the row sum.
* This is not a break; it is a generic-group rho in GT with the Frobenius classes of GT, `2^125.81` units.

## 11. Reproduction

    cd submissions/bls12-381-fixed-g1-invert-g2/hashring-rho-sel1exact
    echo '{"mode":"all","seed":3,"algebra_elements":2,"toy_bits":22,"toy_trials":6,"sim_bits":22,"sim_R":16384,"sim_trials":600,"sim_db":3}' | python3 -I experiment.py
    echo '{"mode":"ledger"}' | python3 -I experiment.py
    # scripts/ lives in the checkout, not in this package: run that line from the checkout root
    cd ../../../..
    bash scripts/run_experiments.sh bls12-381-fixed-g1-invert-g2 submissions/bls12-381-fixed-g1-invert-g2/hashring-rho-sel1exact --preflight-only

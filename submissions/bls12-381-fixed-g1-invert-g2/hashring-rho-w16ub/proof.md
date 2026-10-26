# Frobenius-class-folded distinguished-point rho for fixed-G1 pairing inversion

Track `bls12-381-fixed-g1-invert-g2`. Claim: `time_log2 = 125.817` (125817/1000), success probability at least 0.39.
This package is the accepted `hashring-rho-sel1exact` (125.82) with one algorithmic change: the look-back window shrinks from 24 to 16 class keys, saving 24 scan instructions per slot, and the claim is stated at thousandths.
All objects below are the ones of `tracks/bls12-381-fixed-g1-invert-g2/target-profile.json`; the unit of work is one
dense Fp12 multiplication (`C = 20,000` scalar instructions). Everything the attack does is charged once.

## 0. Summary

1. `P` fixed, `g = e(P,Q0)` generates GT (order `r`). `z = g^a`, return `Q = [a]Q0` (1).
2. Frobenius `x^p = x^lambda`, `lambda = p mod r` order 12: free on `GT \\ {1}`, 12-element classes, `N = (r-1)/12` (2); `sigma^6` = inversion.
3. Selector (3 Fp mults, 3 mod-subs, 1 neg): one Fp2 coordinate `c1` -> 12 first-coords of `c1(sigma^k(x))`, invariant key `kappa` (min), covariant twist `k*` (3). `d = g0(a-b)` one mult (`gamma = g0(1-u)`); min over six stored + negations via one negation of max.
4. Lifted adding walk `y' = y sigma^{-k*(y)}(T_h)`: class map at 1 dense Fp12 mult/step (Lemma 3.2, 4).
5. Fruitless short cycles (the `1/(12R)` per-step two-cycles and every other deterministic closed pattern) are
   removed by a *look-back window* of the last 16 class keys with deterministic probing of the next table indices.
   The rule depends only on the current class and the window, so trails that meet stay together up to a stated
   loss below `1.8*10^-4` (sections 4.3, 6). The window is a 16-slot ring of keys and, in parallel, a ring of their
   stored 64-bit hash words (the word the table index and the distinguished-point test already compute, so no second
   hash is evaluated); the slot counter is advanced by a bounded three-instruction update, an accept is a
   ring write of the six key limbs and the hash, and membership is 16 unrolled load/compare/branch steps over the
   stored hashes with an exact six-limb key test charged on a hash hit (section 5). A hash collision can only reject a
   non-member, which is charged as a probe, never accept a repeat. Residual deterministic cycles longer than the
   window are bounded by an explicit combinatorial lemma (section 7): at most `2^-118.2` per step (lengths 17..512).
   Collision terms are upper bounds throughout: a fresh non-member matches a window holding `d <= 16` distinct hash
   words with probability `d*2^-64 <= 2^-60`, and at most `8*2^-60 = 2^-57` per accepted step by union bound over at
   most 8 probes (probing stops at first acceptance, collision events may overlap; `1-(1-q)^8 <= 8q` for any `q`).
6. Trails start independently (fresh random start each, shared table), end at the first distinguished point (density
   `2^-80`), and are reported to an ordered store. Two reports with the same key give two trails; both are *replayed*
   with exponent tracking and the logarithm is recovered with the Frobenius alignment `sigma^m` (section 4.4). No
   exponents are carried in the main loop.
7. Work: `S_max = ceil(sqrt(2 mu N U)) = 68967528808679891282472564154805017972` accepted steps (`2^125.697`), where
   `U` is the exact rational upper bound on `ln(10000/6097)` of section 5; each multiplication slot
   charged `1 + 1725/20000 = 21725/20000`; all additive terms together are `2^92.361`. Total `W_max <= 2^125.817`
   (`log2 W_max = 125.81678906...`), rounded up to the claim `2^125.817` (section 5). Memory `2^53.6` bytes,
   preprocessing `2^24.13` units, no advice.

## 1. Reduction

Nonidentity `P`: `g = e(P,Q0)` order `r`; `z = g^a` (`a = 0` iff `z = 1`), `e(P,[a]Q0) = z`, `a` uniform indep. of `g`. Charge: one pairing (`2^30` instr), later `[a]Q0`, `g^a == z`.

## 2. Frobenius classes of GT

`r | p^4-p^2+1`; `lambda = p mod r`: `lambda^12 = 1`, `lambda^6 = -1` (`algebra`). `x^p = x^lambda` on GT; `sigma(x) = x^p` order 12. `x != 1`, `0 < k < 12`: `sigma^k(x) = x` needs `r | lambda^k-1`, impossible; 12-element classes, `N = (r-1)/12`, `log2 N = 251.2721`. `sigma^6` = inversion. Bookkeeping: `sigma^m(g^alpha z^beta) = g^{lambda^m alpha} z^{lambda^m beta}`, `sigma^{-k}(T) = T^{lambda^{-k}}`.

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

Same `(kappa, k*)` as the direct twelve-value scan ("smallest index attaining the minimum"): stored-half index `<= 5`, negated half `6..11`, strict `<` within a half. Checked on all `7^6 = 117,649` six-tuples (`ledger`, ties/zeros) and against the direct scan on real GT conjugates (`algebra`). Four modular ops remain: `a-b`, `zeta*a-a`, `zeta*d-d`, `p-M`.

**Two extrema in one pass** (registers: allocation block, section 5). Seed: six limbs of `V_0` into both extrema (6 `LDR`, 12 moves), zero two index registers: 20. Each of five values: `MOVZ` index, 6 `LDR`, six-limb chain candidate-vs-min with seven `CSEL` (13), chain max-vs-candidate (13): 33 each, `20 + 5*33 = 185`. Chains write `XZR` (only carry read). Row `185 + 14 = 199`; same `(kappa, k*)` as two separate scans (strict `<`/`>`). No early exit claimed: equal top limbs possible, so all eleven comparisons carry full chains.

**The eleventh comparison, exactly.** `m <= nb` is the carry of one six-limb chain whose operands are ordered so
that the carry alone is the predicate: `SUBS`/`SBCS` of `nb` against `m` leaves `C = 1` exactly when `nb >= m`, i.e.
exactly when `m <= nb`, equality included. No condition here reads `Z`, which after such a chain is the last limb's
flag and nothing else, and a chain of `m - nb` would say only `m >= nb`. Six `CSEL` with condition `CS` take the
minimum (or the negated maximum), one `CSEL` takes `i_m` (or `i_M + 6`) and one `ADD` forms `i_M + 6`: 14.
`experiment.py` executes this tail with the negation on a register-file machine (218 instructions per tuple, 212 when
the maximum is `0` and the negation skips) and checks the result against the `<=` reference and the direct
twelve-value scan on 245 small-field tuples (51 with `m = nb`, one all-zero) and 40 production tuples. The negation is
`0 - M` in x10-x15: the chain, then `B.CS` on its carry (taken exactly when `M = 0`, which leaves `0`), else `ADDS`
and five `ADCS` of `p` from six registers; `ADC` alone would not update the carry, and no mask is formed.

**Instruction count of one selection (blst, AArch64, prices of `tracks/docs/COST_MODEL.md`).**
3 `blst_fp_mul` x 383 = 1,149; three modular subtractions at 31 and one negation at 13, plus 18: one six-limb load of
`p` in each of the three blocks (`p` is a fixed frame offset, and no register survives a kernel or the extrema), 124;
the two extrema in one pass and the eleventh comparison, `20 + 5*33 + 14` = 199. Total `1,149 + 124 + 199 = 1,472`
(`experiment.py`: `selector_instructions_total`), the three rows that open section 5's table. The Montgomery constants
`zeta, g0` are two stored Fp elements; no conversion out of Montgomery form is needed.

## 4. The algorithm

### 4.1 Parameters and primitives

| symbol | value | meaning |
| --- | --- | --- |
| `R` | `2^14` | table indices |
| `w` | 16 | look-back window (class keys of the last 16 states of the current trail, including the current one) |
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
(80 bits); both are functions of the class. The window of section 5 stores `w_1` itself as each entry's hash, so the
look-back membership computes no word beyond the two the index and the test need.

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
   ring as `kappa_0` in all 16 key slots and its hash word in all 16 hash slots, with `p = 0` (section 5), so the
   window is the one-element set `{kappa_0}`.
2. At state `n`: if `kappa_n` is distinguished, **report** `(kappa_n, i)` and end the trail. If the trail's slot counter
   `c` (incremented by one for every multiplication in step 3, accepted or rejected) has reached `L_cap`, end the trail
   without a report.
3. Otherwise `h = hidx(kappa_n)` and for `a = 0, ..., A-1`: `y' = y_n * TW[(h+a) mod R][k_n]` (one dense
   multiplication); `(kappa', k') = Sel(y')`; if `kappa'` equals any of the `w` keys in `W`, this probe is
   rejected; else accept it: `y_{n+1} = y'`, append `kappa'` to `W` (dropping the oldest beyond `w`), `n += 1`, go to 2.
   The append of section 5 is the ring write at slot `p` plus the bounded slot update; a candidate is rejected when
   its 64-bit word equals a stored one, a rule exact on every repeat and over-rejecting only on a collision, which is
   charged as an escape (section 5). If all `A` probes are rejected the trail ends without a report.

   The identity needs no special test: given `y_n`, the table element is uniform in `GT` (section 2 and H3, since
   `(u_h, v_h)` is uniform on `(Z/r)^2` and so is `u_h + a v_h` for `a != 0`), so `Pr[y' = 1] = 1/r` per slot and at
   most `2^-128.857` over the run; the selector then returns `kappa = k* = 0`, the next step is ordinary and no
   recovery is affected (the identity entry of `L4`, section 6).

By Lemma 3.2 the class of `y_{n+1}` is the class of `c_n * T_{h+a}` where `c_n` is the canonical representative of
the class of `y_n`; so, **if the first probe is accepted**, the class sequence is a deterministic function of the
current class alone, `F(C) = class(c * T_{hidx(kappa(C))})`. A probe beyond the first happens only when the first
candidate class lies in the last-16-class window (a closed deterministic pattern of length at most 16, section 7).
The decision depends only on the current class and the 16 previous class keys.

**Exponent invariant (used only in replay).** `y_n = g^{alpha_n} z^{beta_n}` with
`(alpha_{n+1}, beta_{n+1}) = (alpha_n, beta_n) + lambda^{-k_n} (u_{h'}, v_{h'})`, `h'` the accepted index, because
`sigma^{-k}(T) = T^{lambda^{-k}}`.

**Coalescence.** Let two trails reach the same class `C` at states `n` and `n'`. The next accepted class is `F(C)`
for both, unless it lies in one trail's window (probability about `1/(12R)` per step each, section 6, `L2`), and the
windows are identical after 16 further steps. The hash rule adds at most `16*2^-64 = 2^-60` per candidate to that rejection
probability and can only reject more, so the two trails coincide after the meeting except with probability at most
`2*24*nu`, `nu <= 1.1/(12R)`, hence both reach the same distinguished point.

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
slots, charged). The algorithm halts with failure if `M` exceeds `M_max + L_cap`, starts exceed `ST_max`, reports
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
| selector: three modular subtractions at 31 and one negation at 13 (the `a-b` for `d`, the two `zeta^2` terms, the negation of the maximum; section 3), plus 18: one six-limb load of `p` in each of the three blocks (the register-allocation block at the end of this section) | 124 |
| selector: the two extrema in one pass, `20 + 5*33`, and the eleventh comparison (6 `CSEL` + 1 `CSEL` + 1 `ADD`; section 3) | 199 |
| two `khash` words: 12 x (one limb `EOR` and the 8-instruction mix) over six limbs held in registers, plus the two accumulator loads and the two splitmix64 multiplicand loads (no callee-saved register holds either constant across a slot) | 112 |
| distinguished-point test: `AND x3, x2, #0xFFFF`, `ORR x3, x3, x1`, `CBNZ x3` — 64-bit registers, so 80 bits | 3 |
| look-back membership: 16 unrolled slots of `LDR` the stored 64-bit hash, `CMP` with `w_1` in `x1`, `B.EQ` = 48, plus at most one hit path of 20 (key-slot address, six key limbs, `6 EOR` + `5 ORR` + `CBZ` exact equality, the branch that ends the probe) = 68 | 68 |
| look-back update, ring write: key-slot address (three shifted adds), three `STP` of the six key limbs, and the 64-bit hash store with its charged address `ADD` (`ADD x11, sp, p, LSL #3 ; STR x1, [x11]`) | 8 |
| look-back bookkeeping: the ring slot's bounded update, the accepted-step counter with its cap test, two branches | 12 |
| state and call plumbing, table entry address: `LDR h`, `ADD h+a`, `AND (h+a) mod R`, `MOVZ 6912`, `MADD`, `LDR k_n`, `MOVZ 576`, `MADD`, `MOV x1, y`, `EOR x0, y, #0x400` | 10 |
| state and call plumbing, accept tail: `AND h', w1, #0x3FFF`, `STR h'`, `STR k_n`, `MOVZ a, #0`, `EOR y, y, #0x400` | 5 |
| state and call plumbing, probe ladder on the reject path: `ADD a, #1`, `CMP a, #8`, `B.LT` | 3 |
| selector operand base, twice: `EOR x29, x20, #0x400`, after the dense kernel and after the three `fp_mul` kernels | 2 |
| the three `blst_fp_mul` argument setups, three `ADD` each | 9 |
| trail slot counter with its 126-bit bound loaded | 7 |
| global slot counter with its 126-bit bound loaded | 7 |
| next-slot branch | 1 |
| **subtotal** | **1,719** |
| less the probe ladder, which no accepted slot executes | -3 |
| **maximum over slot types** | **1,716** |
| explicit margin (9/1716 = 0.5245%, at least the 0.5% claimed) | 9 |
| **charged per slot** | **1,725** |

So each slot costs `1 + 1725/20000 = 21725/20000` units. The 22,996-instruction blst kernel already contains the
loads and stores of its operands, so table reads are not charged again. The probe ladder runs only on the reject
path, so no slot executes every row: the row sum is `1,719` and the maximum over slot types `1,716`; the
charge `1,725` is `9` above that maximum, `9/1716 = 0.5245%` (`experiment.py`, `ledger`:
`slot_margin_instructions`, `slot_instructions_max`), and 6 above the row sum. `1,728` is also the ceiling: at `1,729`
the least thousandth at or above `log2 W_max` becomes `125.818` (`experiment.py` finds it by exact integer search), so
the charge sits `3` instructions below it. The
sensitivity is `M_max/(20000 W ln 2) = 6.6*10^-5` bits per charged instruction, i.e. `6.641*10^-5` at the charged
`1,725`, and `7.2*10^-5 = 1/(20000 ln 2)` at the reference rate (the charge enters once, through the factor
`1 + c/20000`, so the second is an upper bound at every count; both from `experiment.py`, mode `ledger`).

**The instruction sequence of each row.** The modular
operations, the two extrema with their indices, the eleventh comparison, the distinguished-point test, the scan, its
hit stub, the ring write, the seed and the ring slot are executed in turn on a register-file machine
(`experiment.py`: class `A64`), which counts the instructions it runs and checks their result against a Python
reference; the three `blst_fp_mul` calls (1,149, the organizer's 383-instruction kernel price), the `khash` finaliser
(112), the bookkeeping row (12) and the plumbing rows (table entry address 10, accept tail 5, probe ladder 3, operand
bases 2, argument setups 9, the two slot counters 7 and 7, the closing branch 1) are priced from the written-out
sequences below and are not executed. A modular subtraction `x - y mod p` (three of them: `a - b` for `d`, `zeta*a - a`, `zeta*d - d`) is
`6 LDR` of `x` into x0-x5, `6 LDR` of `y` into x6-x11, the six-word `SUBS`/`SBCS` chain `r = x - y` in place (whose
final carry `C` is 1 exactly when `x >= y`), `B.CS` past the correction (taken exactly when that chain did not
borrow), and, on the borrow path, one `ADDS` and five `ADCS` of `p` from x12-x17, then `6 STR`: 31 when the correction
runs and 25 when the branch is taken, so 31 bounds both paths; there is no call or argument setup because the
selector is one inline block that stores its results in its six-value buffer. The chain's carry propagates only
through `ADCS`; `ADC` would leave the carry of the first limb in place and produce a wrong residue, and
`x - y + (p if x < y)` is then the canonical residue in `[0, p)`. No borrow mask is formed: the branch replaces it,
so the subtraction writes no register outside x0-x17 and leaves `p`'s six registers intact for the next subtraction of
its block. `p`'s limbs are six `LDR` from a fixed frame offset, but no register holds them across a block: the kernels
clobber x0-x18, x19-x28 are the ten fixed live values, and the extrema occupy all of x0-x18 plus x29 and x30. So the
four operations sit in three blocks, each paying one six-limb load, shared inside its block (allocation block
below): `a - b` before the three `blst_fp_mul` calls, since it is
an input to the third; `zeta*a - a` and `zeta*d - d` after them, on the kernels' outputs, sharing one load; and the
negation after the extrema. The negation `p - M` has the same shape with `XZR` as the subtrahend, its six-limb chain
in x10-x15, the same `B.CS` (taken when `M = 0`, leaving `0`) and no store: 13, `0` for `M = 0` and the canonical
`p - M` otherwise; with its load it is 19, so the row is
`3*31 + 13 + 18 = 124`, that is `37 + 68 + 19` (one subtraction with its load, two under one load, the negation). The
selector's eleven comparisons are itemised in section 3. Each `khash` word is six rounds of one `EOR` of a limb into
the accumulator and the splitmix64 finaliser (`LSR`, `EOR`, `MUL`, `LSR`, `EOR`, `MUL`, `LSR`, `EOR`: 8), so
`6*9 = 54` per word and 108 for both, with no limb loads because the limbs are the registers the selector's final
`CSEL` wrote (`x4-x9`); the two hoisted accumulators `mix(s_1)` and `mix(s_2)` are reloaded per slot (2) and so are
the two splitmix64 multiplicands, which no callee-saved register can hold across the dense kernel (2), so the row is
112.
The distinguished-point test is `AND x3, x2, #0xFFFF`, `ORR x3, x3, x1`, `CBNZ x3` (3), on 64-bit registers: `x1`
holds `w_1` and `x2` holds `w_2`, so the test reads all 64 bits of the first word and the low 16 of the second, 80
bits. Its 32-bit form would read `w_1 mod 2^32` and `w_2 mod 2^16`, a 48-bit predicate with `S_max*2^-48 = 2^77.7`
reports, far above `D_max = 2^47`; the executed case `w_1 = 2^32`, `w_2 = 0` is the witness that separates them.
One membership slot is `LDR x10, [sp, #8*i]`, `CMP x10, x1`, `B.EQ` that slot's hit stub (3), a 64-bit load and a
64-bit register compare against `w_1` in `x1`; the unsigned scaled immediate covers every `i <= 15` (`8*15 = 120`).
The hit stub is `ADD x11, sp, #(K_OFF + 48*i)` (1; `K_OFF + 48*15 =
912` is inside the 12-bit ADD range), six `LDR` of the stored limbs, six `EOR` against the candidate's limbs folded
by five `ORR`, `CBZ acc` and the branch that ends the probe: 20. The test is exact with no sixth-limb exception:
both keys are canonical representatives below `p < 2^381`, so the six XORs vanish exactly when the keys are equal.
At most one stub runs per candidate, so the row is `16*3 + 20 = 68`: a miss is 48 and a hit in the last slot is the
maximum, on both the equality and the collision path. `hash_scan` and `ring_check` execute the program and assert 48,
68 and 68 charged (section 9).
The ring write is three shifted adds for `48p` (`ADD x10, sp, #K_OFF`, `ADD x10, x10, p, LSL #4`, `ADD x10, x10, p,
LSL #5`), three `STP` of the six key limbs, and the hash store's address `ADD x11, sp, p, LSL #3` with `STR x1,
[x11]`: 8. `STR w_1, [sp, p, LSL #3]` is not encodable (the 32-bit register-offset store allows `LSL #0` or `#2`
only) and would store the low half of the word, after which a 64-bit scan could miss a repeat. The bookkeeping is
`ADDS`/`ADCS` of the 126-bit
accepted-step counter, `SUBS`/`SBCS` against `S_max` and a branch (5), with the two `LDR` of `S_max` from its frame
offset (7), the slot update `ADD p, p, #1`, `CMP p, #16`, `CSEL p, p, xzr, NE` (3: three register operands, taken
when `Z = 0`, no immediate), and two branches: 12. The plumbing rows follow, in the loop's order: the table entry's
`LDR h`, `ADD h+a`, `AND (h+a) mod R`, `MOVZ 6912`, `MADD`, `LDR k_n`, `MOVZ 576`, the second `MADD`, `MOV x1, y`
and the result pointer `EOR x0, y, #0x400` (10, the dense kernel being called as `blst_fp12_mul(x0, x1, x2)` because
it clobbers x0-x18); the accept tail `AND h', w1, #0x3FFF`, `STR h'`, `STR k_n`, `MOVZ a, #0`, `EOR y, y, #0x400`
(5); the reject path's `ADD a, #1`, `CMP a, #8`, `B.LT` (3); the operand base `EOR x29, x20, #0x400` twice, because
x29 is scratch (2); the three `blst_fp_mul` argument setups, three `ADD` each (9); the trail and global slot counters,
`LDR`/`LDR` of their 126-bit bound, `ADDS`/`ADCS`, `SUBS`/`SBCS` and one branch each (7 each); and one closing branch
(1). The accepted-step test is needed only when a trail starts and the update rows (`7 + 12`) run only on an accept:
both are charged at every slot, an over-charge over the `2^-17.44` of rejected probes (section 9).

**The look-back window's structure.** The window holds the class keys of the last at most 16 states of the trail,
including the current one, as a set; membership is all the escape rule of section 4.3 uses, and the oldest key is
dropped on every accept. It is seeded with `kappa_0` and stored as two fixed arrays of 16 entries:

- `K[0..15]`: the keys, 48 bytes each, as a ring. The k-th state's key goes to slot `k mod 16`, the slot of the
  oldest stored key, which one register `p` names; `p` advances by the bookkeeping row's bounded `ADD`/`CMP`/`CSEL`
  update.
- `H[0..15]`: one 64-bit word per slot, the hash `w_1 = khash(kappa, s_1)` of the key stored there, the word the
  table index (its low 14 bits) and the distinguished-point test (`w_1 = 0`) already compute (section 4.1): the
  window evaluates no hash at all, it stores one word per accept.

The trail start seeds both rings: `kappa_0` in all 16 key slots, its word in all 16 hash slots, `p = 0`, 81
instructions (a `MOVZ`, and per slot an `ADD`, three `STP` and the 64-bit `STR x1, [sp, #8s]`), inside the start
item's slack. Seeding every slot with `kappa_0` keeps the ring exact: after `k` accepts the stored set is
`{kappa_0, ..., kappa_k}`, the last `min(k+1, 16)` states, with no stale key and no empty-slot sentinel. Once `k >= 15`
the accept overwrites `kappa_{k-15}`, and a closed pattern of length `l <= 16` returns to a key still in the window,
which is what section 7 uses.

Membership of a candidate `kappa'` is the 16-slot scan: its own `w_1` against each stored word, the first match
entering that slot's hit stub, which loads the stored key and runs the exact six-limb test. The scan is a
**conservative** test with a hard bound: every repeat is caught, because a repeat has a stored key's word, so the
scan reaches that slot (or an earlier one) and rejects, and no pattern of section 7 is missed. A rejection can also
happen with no stored key equal to `kappa'`, when its 64-bit word collides with a stored one. Under H2 (section 6)
the candidate's word is uniform and independent, so conditional on the window holding `d` distinct hash words
(`d <= 16`; `d = 1` immediately after seeding, since all slots hold `kappa_0`) the probability of such a
collision rejection is `d*2^-64 <= 16*2^-64 = 2^-60` per candidate, an upper bound, not an equality. Per accepted
step there are at most 8 probes, so the collision-rejection probability is at most the union bound `8*2^-60 = 2^-57`
as arithmetic on that bound: probing stops at the first acceptance and collision events across probes may overlap,
and `1-(1-q)^8 <= 8q` for any `q`, so the bound holds regardless. Such a *collision escape* costs exactly what
an ordinary escape costs: the probe is rejected and its multiplication charged; it cannot lose a repeat and cannot
change the stored keys. The escape rate is therefore unchanged at `nu <= 1.1/(12R) < 2^-17.44`,
the collision term being at most a relative `2^-39.5` of its two-cycle term and far inside its 10% allowance; the rate
measured in section 9 is measured on this same rule.

**Why a 64-bit word and not a mask.** A 16-bit mask indexed by a 16-valued function of the key cannot represent 16
distinct 384-bit values: every candidate whose function value is set would be rejected, an escape rate of order
`0.63` per step, and one pointer per bucket lets two window keys share a bucket. The stored word gives no false
negative, at most `16*2^-64 = 2^-60` per candidate and the exact test on a hit. The sorted one-byte index array this line replaces
answers where the candidate would be *inserted*, a question the rule never asks: 259 instructions against the 76
here.

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
| x29 | scratch: `i_M + 6` inside the selector tail, then the y' base | per block |
| x30 | scratch: the extrema's index supply, the chains' `xzr` writes need no temp; the function's link register is saved in the prologue, so a `BL` may set it | per block |
| x0-x18 | kernel argument registers and per-phase scratch; after the selector `kappa` is in x4-x9, `w_1` in x1 and `w_2` in x2 | per phase |
| sp+0 | `H[16]`, one word per ring slot | whole trail |
| sp+192 | `K[16]`, 48 bytes per ring slot (`K_OFF + 48*15 = 912`, inside the 12-bit `ADD` range) | whole trail |
| sp+1344 | `V[6]`, the six stored selector values (`V_OFF`; the extrema scans it) | per state |
| sp+1632 | the modulus `p`'s six limbs (`P_OFF`); `p` the *ring slot* is x27, not a frame word | per block |
| sp+1680, 1728 | the tower constants `zeta` and `g0` | whole trail |
| sp+1776 | `h`, the current state's table index (`LDR` in the head row) | per state |
| sp+1784 | `k_n`, the current state's twist index (`LDR` in the head row, `STR` in the accept tail) | per state |
| sp+1808, 1824 | the two `mix(salt)` accumulator seeds and the two splitmix64 multiplicands | whole trail |
| sp+1840, 1856, 1872 | the three 126-bit constants `L_cap`, `M_max + L_cap`, `S_max` | whole trail |
| sp+1888 | the `a - b` scratch the third `blst_fp_mul` reads | per state |
| sp+2048, sp+3072 | the two 576-byte y buffers, `0x400` apart with bit 10 clear in the lower one | whole trail |

| phase | live scratch values | registers | loaded or stored per slot |
| --- | --- | --- | --- |
| table entry address and the call's arguments | `h`, `a`, `k_n`, two strides, the address, three pointers | 8 of x0-x18 | `h` (1), `k_n` (1), the two strides (2), `MOV x1, y` (1) |
| `blst_fp12_mul` call | three argument pointers | 3 of x0-x18 | none (x19-x28 survive a call) |
| subtraction blocks (i), (ii): `a - b`, then the two `zeta^2` terms | `x` (6) in x0-x5, `y` (6) in x6-x11, `p` (6) in x12-x17, the y' base in x29; no borrow register (the `B.CS` replaces the mask) | 18 of x0-x17, x29 untouched | `p`'s six limbs, 6 in (i) and 6 shared by (ii)'s two subtractions |
| subtraction block (iii), inside the tail: the negation `p - M` | `M` (6) in x10-x15, `p` (6) in x0-x3, x16, x17; the tail's key, indices and index supply stay live | 12 of x0-x17 | `p`'s six limbs (6) |
| `blst_fp_mul` argument setups | three pointers per call, three calls | 3 of x0-x18 | none: three `ADD` per call (9) |
| extrema network, seed | minimum (6) + index (1), maximum (6) + index (1), candidate (6), the y' base | 20 + x29 | `V_0`'s six limbs come from the y' base, `V_1`-`V_5` from their frame words |
| extrema network, loop | minimum (7), maximum (7), candidate (6), index supply | 21 = x0-x18 + x29 + x30 | none |
| `khash` | `kappa` in x4-x9 (6), `w_1` in x1, `w_2` in x2, one temp, two multiplicands | 11 of x0-x18 | the two seeds (2) and the two multiplicands (2) |
| scan and hit stub | `kappa` in x4-x9, `w_1` in x1, `w_2` in x2, the scan temp x10, the stub address x11, six `EOR` results | 16 of x0-x18 | none |
| ring write | `kappa` in x4-x9, `w_1` in x1, the key address x10, the hash address x11 | 9 of x0-x18 | none |

**The subtraction's own registers, executed.** The three subtractions write only x0-x17 and the negation only
x10-x15 plus the six `p` registers x0-x3, x16, x17; none writes x19-x28 or x29, so the ten fixed live values and the
y' base survive, and `B.CS` replaces the borrow mask so no extra register is needed. `experiment.py`
(`modsub_block_check`) seeds x19 = `&TW`, x20 = the y buffer, x21/x22 = `M`, x23/x24 = `c`, x25/x26 = `S`, x27 = the
ring slot, x28 = the probe counter and x29 = the y' base, and asserts after every subtraction that no register
outside x0-x17 changed, that x12-x17 still hold `p`, and that only the six destination frame words were written.
Block (ii) is executed as the composition it is priced as: two successive subtractions under ONE load of `p`, 68
instructions at most for the pair, both results checked against the Python reference on the printed tuples. The
rejected shape is the negative control: the same check catches `CSETM x19` and a load of `p` into x20-x25. Both branch
paths are exercised (176 tuples, 68 correcting and 108 skipping), and a block (i) tuple costs 37 with its load when
the correction runs and 31 when it does not.

The three 126-bit constants (`L_cap`, `M_max + L_cap`, `S_max`), `p`'s six limbs, `h` and `k_n` have no callee-saved
register: the ten fixed values fill x19-x28 and the kernels clobber x0-x18, so they live in frame words and the rows
above pay for every load: the `2 LDR` inside the two counter rows' 7 and the bookkeeping row's 7, the 18 in the
selector's row, and one `LDR` each for `h` and `k_n`, whose `STR` pair is in the accept tail. The only state that
survives a rejected probe and is *not* in a register is therefore those two words, and both are charged on every path.
Nothing else is spilled: the kernel writes the candidate into the other y buffer and x20 is toggled only on an accept.

**The two hash words.** The index and the distinguished-point predicate are both functions of the class key
(section 4.1): `w_1 = khash(kappa, s_1)`, `w_2 = khash(kappa, s_2)`, the index `w_1 mod 2^14`, so no third word is
computed. A state is distinguished iff `w_1 = 0` and `w_2 mod 2^16 = 0`, 80 bits. Distinguished states are reported
and never stepped from, so each value of the low 14 bits of a stepped-from state has probability
`2^-14 (1 +- 2^-50)`. H2 is stated on the two words themselves, and both are charged at every slot, which over-charges
the `2^-17.44` of rejected probes.

**Main walk.** `N = (r-1)/12`, `mu = 11/10`, `x = 3903/10000`, `L = -ln(1-x) = ln(10000/6097)`, and

    U = sum_{k=1}^{64} x^k/k + x^65/(65(1-x)),

an exact rational. The partial sum is below `L` and the geometric tail dominates the terms it drops, so `U` is an upper
bound for `L`, and `0 < U - L < 2^-99`: the computed gap is `6.63e-31`, and the bound is itself a rational comparison,
because `U` minus the partial sum through 70 terms is below `2^-99` while that partial sum is below `L`.

The accepted-step budget is computed from `U`, the larger rational bound. For a target main-event probability
`3903/10000 = 0.3903` (the rule's `0.39` plus the `1.8*10^-4` loss bound of section 6) it is the least integer with
`S^2 >= 2 mu N U`, checked exactly in integers:

    S_max = ceil(sqrt(2 mu N U)) = 68967528808679891282472564154805017972   (log2 = 125.69726)
    M_max = S_max + floor(S_max/2^13) + 1 = 68975947696473919589513881606484071319,  work 21725/20000 * M_max.

No logarithm, float or decimal enters `S_max`: `U` is rational, the ceiling is exact, and `S_max^2 >= 2 mu N U >
(S_max-1)^2` is an integer comparison. The `L`-based ceiling is a different integer: `U/L - 1` is `2^-99.24`, so
`ceil(sqrt(2 mu N L))` is about `2^25.46` steps lower, and taking `U` costs at most a relative `2^-99`. Every
inequality here, in section 6 and in the claim is stated on this `S_max`, the integer `experiment.py` prints, and the
exact checks `W^1000 <= 2^125817` and `W^1000 > 2^125816` are formed from it.

The probe allowance `S_max/2^13` is `2^-13` of the accepted steps; the escape (probe) rate is at most
`nu <= 1.1/(12R) < 2^-17.44` per step by `L2` below (measured `7.72e-6` at the production table size in section 9).
The number of rejected probes is therefore dominated by a binomial with mean `S_max * 1.1/(12R) < 2^108.25`, and the
allowance `2^112.7` exceeds that mean 21-fold, so exceeding it needs a deviation whose binomial tail is below
`exp(-2^107)`; the meter fails on such an excess and `L7` carries this term. The allowance is a hard cap on the
charged work, not an expectation.

**Additive items** (log2 of each, fixed by caps, so hard bounds):

| item | amount (units) | log2 |
| --- | --- | ---: |
| two pairings (`g` and one spare) at 2^30 instructions each | `2^31/20000` | 16.72 |
| table: `R * (2*509 + 1 + 11 + slack) <= R * 1112` | | 24.12 |
| starts: `ST_max * 1100` (1,019 multiplications + ChaCha20 and two reductions < 0.4 unit) | | 57.11 |
| reports: `D_max * 2^16 / 20000` (red-black tree and transfer; a tree operation is at most 2,300 instructions) | | 48.72 |
| replays: `REC_max * 2 * L_cap * 5/4` (slot 21725/20000 plus exponent update < 0.075: two 255-bit modular products by the stored `lambda^-k` and two additions, at most 1,500 instructions, with `(u_h, v_h)` read from a stored 2^20-byte array, no PRG call; `21725/20000 + 0.075 = 1.16125 <= 5/4`) | | 92.33 |
| completing the last trail: `L_cap * 21725/20000` | | 87.13 |
| `[a]Q0` and `g^a == z` | `<= 2^12` | 12.00 |
| **additive total** | | **92.361** |

**Total.** `W_max = 21725/20000 * M_max + 2^92.361 <= 2^125.817` (`log2 W_max = 125.81678906...`); the claim
`2^125.817` is the least thousandth at or above it, and `experiment.py` (mode `ledger`) proves that with exact integers
alone: `m` is the least integer with `W^1000 <= 2^m`, found by binary search on `m` over `[0, 1000 * bitlength(W)]`, so
`W^1000 <= 2^125817` together with `W^1000 > 2^125816` (no smaller thousandth suffices). No float enters that step.
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
| L2 | an escape of one of two meeting trails inside the 16 steps after meeting (`1.1/(12R) < 2^-17.44` covers the hash rule's at-most-`2^-57` collision term of section 5, relative at most `2^-39.5`) | `<= 2*16*(1.1/(12R)) <= 1.8e-4` |
| L3 | one of those two trails is trapped in a residual closed pattern, loops in a within-trail rho cycle, or has all 8 probes rejected (per trail: `2^80 * 2^-118.2`, `E[L^2]/(2N) = 2^-91`, `2^80 * nu^8 <= 2^-59.5`) | `<= 2^-37` |
| L4 | ties among the 12 values `V_k` (8 linear conditions on `(a, b)`: `a = 0`, `d = 0` and `d = zeta^t a` for six `t`, so probability at most `8/p < 2^-377.6` per state) or `c1 = 0`; a state equal to the identity, `1/r < 2^-254.84` per slot (section 4.3); over `2^126` slots | `<= 2^-128.84` |
| L5 | false store match or zero denominator | `<= 2^-250` |
| L6 | sampling distance, `2^47` starts and `2^14` table scalars at `2^-129` | `<= 2^-80` |
| L7 | caps `D_max, ST_max, M_max` exceeded (Poisson tail of a count with mean `2^45.7`, `2^46` at most) | `<= 2^-1000` |
| L8 | more than `REC_max` unusable matches | `<= 2^-200` |

Total loss below `1.8*10^-4` (`L2` is the only label above `10^-11`). Hence, under H1-H3, success probability at
least `0.3903 - 0.00018 = 0.39012 > 0.39`, a printed margin of `0.00012` over the required `0.39`.

*What H1 must supply, and the margin.* The claim uses H1 only at the single point `s = S_max`: the probability that
the first repeat occurs by `S_max` is at least `3903/10000` (equivalently `S^2/(2N) <= mu U` with probability
`>= 0.3903` for the first repeat `S`). Success `>= 0.39` after the losses needs `0.39027`, so the true inefficiency
must satisfy `mu_true <= 1.1 * ln(1/0.6097)/ln(1/0.60973) = 1.1 * 0.494788/0.494739 = 1.10011`: the target sits at
the rule's requirement, so the claim spends the whole `mu = 11/10` allowance and keeps `0.01%` of spare beyond it.
The measured `mu_hat = 0.996 +- 0.014` (section 9) is `7.4` standard errors below the allowance; for an ideal random
map the value at `S_max` is `1 - (6097/10000)^{1.1} = 0.4197`, and section 9 reports that statistic itself, not only
the mean.

**Uniform distinguished-point and routing fields (H2).** The only structure used from the class key is: the two hash
words `w_1, w_2` (hence `hidx`, the distinguished-point predicate and, since section 5 stores `w_1` in the window, the
window's hash words with their at-most-`2^-60` collision term) and `k*` behave as independent uniform fields
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
covers (a window revisit of this kind has probability about `16^2/N < 2^-241`; it is not a structure of the walk).
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

`experiment.py` (mode `ledger`) evaluates this sum in log scale with float logarithms, so the printed `2^-118.9` and
`2^-118.2` are float evaluations of the exact integer sum, each rounded outward (up) from the computed `2^-118.98` and `2^-118.29`; the only exact comparison is the final claim inequality.
The maximum over `17 <= l <= 512` is at most `2^-118.9` (at `l = 17`) and the sum is at most `2^-118.2`. For `l > 512`:
`C(R,d) <= (eR/d)^d` makes each term at most `((e/12)(d/R)^{l/d-1})^d`; for `d <= R/2` that is at most
`0.2266^d 2^{-(l-d)} <= 2^-l`, and for `d > R/2` at most `0.2266^{R/2} < 2^-17000`. So `E(l) <= 2^-480` for all
`l > 512` and their sum to the trail cap is below `2^-400`, negligible next to `2^-118.2`.
For lengths `l <= 16` a trap is detected within `l` steps by the window and escaped (the hash rule of section 5
rejects every repeat of a window key, so no such trap is missed); the unconditional escape
probability per step is at most `nu <= 1.1/(12R)` (dominated by the two-cycle: the new class has the same index as the
previous one, probability `1/R`, with twists differing by 6, probability `1/12`; the hash rule adds at most
`8*2^-60 = 2^-57` by union bound, which is at most a relative `2^-39.5` of the bound), measured `7.72e-6` against the
bound `1.1/(12R) = 5.59e-6` and `1/(12R) = 5.09e-6` in section 9 (the excess over the bound at 30-bit toy size is the
`16/N = 2.4e-7`-scale random-revisit rate plus sampling noise over 300 groups, `~1.4e-6` standard error; the bound
covers the deterministic component, and the random-revisit component vanishes below `1e-6` by 34 bits and is zero at
production `N`). Longer escapes (lengths 3 to 16) are of order `R^-2` or smaller and are inside the
factor 1.1.

*Cost of undetected traps.* A trail of mean length `2^80` is trapped with probability at most `2^80 * sum_l E(l) <=
2^-38.2`, and then runs `L_cap = 2^87` slots (128 mean lengths), so the expected wasted fraction of the main walk is at
most `L_cap * sum_l E(l) <= 2^87 * 2^-118.2 = 2^-31.2`. A trail dies (all `A = 8` probes rejected) with probability at
most `nu^8 <= 2^-139.5` per step (the 8 candidate classes use 8 different indices; the two-cycle event can hold for at
most one of them, so the product bound is far from tight, and independence is the heuristic H2); per trail
`2^80 * nu^8 <= 2^-59.5`.

## 8. Memory, preprocessing, advice

| item | bytes (bound) |
| --- | --- |
| store: `<= D_max = 2^47` nodes of 96 bytes (48-byte key, 16-byte trail index, three pointers, colour) | `2^53.585` |
| twist table `12 * 2^14 * 576` | `2^26.75` |
| scalars `(u_h, v_h)` for replay: `2^14 * 64` | `2^20` |
| trail state: `y` and the probe buffer (2 x 576 B), window (16 keys of 48 B = 768 B and 16 stored hash words = 128 B), the selector and hash buffers, counters, ChaCha20 state | `< 2^20` |
| **peak** | `<= 2^53.6 < 2^54` bytes (`2^53.585 + 2^26.75 + 2^20 + 2^20`) |

All of it lives in one address space below `2^64`; every store transfer is inside the 2^16-instruction report cost.
Preprocessing is the table and one pairing, `2^24.13` units. There is no advice.

## 9. Experiments and their scope

One declared experiment (`hashring-rho-checks`, the verifier allows one initial run) executes four parts in one
process (`{"mode":"all", ...}`, about 10 s on a laptop; the parts are also callable alone by `mode`). None is a
production-scale run and none measures hardware instruction counts or the production success probability; the counts
printed below are those of the written-out sequences, executed here in Python.

* **ledger** recomputes sections 5, 7 and 8 in exact integer and `Decimal` arithmetic where the quantity is exact
  (`S_max`, `M_max`, each additive item, `log2 W_max = 125.8167890616 <= 125.817`, with the exact integer checks
  `W^1000 <= 2^125817` and `W^1000 > 2^125816`, the maximum slot `1,716`, the charge `1,725` and the ceiling `1,728`,
  `2^53.585` bytes, `2^24.127` preprocessing), the trap bounds in float log scale as stated in section 7, and it
  re-checks the one-pass min-max network against two separate scans, and its result against the direct twelve-value
  scan, on 117,649 six-tuples (ties and zeros included). `selector_min_check` and `ring_check` then execute the priced
  AArch64 lists on the register-file machine `A64`, which counts what it runs and checks each result against a Python
  reference: the selector tail (218 instructions per tuple, 212 when the maximum is `0` and the negation skips, on 245
  small-field tuples — 51 ties and the all-zero one — and 40 production tuples), the modular subtraction in all three
  of its blocks (31 correcting, 25 when `x >= y`, 37 and 31 with the block's load of `p`), the composition of two
  successive subtractions under ONE load of `p` (68 at most), with x19-x28 and x29 seeded from the allocation and
  asserted byte-identical after every subtraction, `p`'s six registers asserted unchanged and the negative control
  (the mask shape over x19, x20-x25) caught, the distinguished-point test (3, on 314 operand pairs, with the `w_1 =
  2^32` witness that separates 80 bits from 48), the slot update for 200 accepts (3 each, `p = k mod 16`, always the
  oldest key's slot), the seed (81), the ring write (8) and the membership scan with its hit stub (48 a miss, 68 a
  hit in slot 15, on the equality and the collision path alike), over every fill 0..16 and 323 candidates
  cross-checked against a plain set and against integer equality on 3,412 key pairs. The plumbing rows of section 5
  (10, 5, 3, 2, 9, 7, 7, 1), the bookkeeping row (12), the 1,149-instruction kernel price and the 112-instruction
  `khash` finaliser are priced, not executed, and section 5 says so for each.
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
  instance. Its window is the production ring structure itself (keys and their hash words, the same `w_1` the index
  uses, the same `ADD`/`CMP`/`CSEL` update, the same priced scan and six-limb test), with every membership decision
  cross-checked against a plain list of the same keys.
* **sim**: the identical walk in the exponent domain of a cyclic group of order `r'` with automorphism `lambda'` of
  order 12 (state = discrete logarithm; keys are 64-bit hashes of the 12 conjugates), with `R = 2^14`, window 16 and
  probing 8; short trails keep within-trail rho cycles negligible at toy sizes. It measures `mu_hat = mean(S^2/(2N))`
  at the first repeat of a class by a different trail, the same estimator on a plain-birthday control, the escape rate,
  and whether the meeting trail reaches the same distinguished point, and the fraction of groups whose
  first repeat satisfies `S^2/(2N) <= 1.1 U = 0.54427` (inside the `S_max` budget). Default run
  (600 groups of 22 bits, trail mean 8, seed 3): `mu_hat = 1.031 +- 0.041`, fraction within budget `0.4167`, control
  `1.058`, 599 of 600 merged trails reach the same distinguished point. The observed escape rate there is `3.2e-5`,
  above `1/(12R) = 5.1e-6` because at `N = 2^18` the window also catches random revisits (about `16/N` per step); at
  28 to 34 bits that term is below `10^-6` and the rate matches `1/(12R)` plus the stated random-revisit term.

Offline runs of the same code with larger inputs (`python3 -I experiment.py`, `{"mode":"sim",...}`):

| bits of r' | R | trials | seed | `mu_hat` | s.e. | fraction within `S_max` budget | control | escape rate / predicted `1/(12R)` | merged trails reaching the same DP |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| 28 | 64 | 1500 | 105 | 1.072 | 0.026 | 0.4053 | 1.085 | 1.26e-3 / 1.30e-3 | 1493 / 1500 |
| 30 | 256 | 1500 | 104 | 0.992 | 0.026 | 0.4073 | 1.028 | 3.10e-4 / 3.26e-4 | 1500 / 1500 |
| 30 | 1024 | 1500 | 103 | 0.999 | 0.026 | 0.4173 | 0.935 | 7.68e-5 / 8.14e-5 | 1500 / 1500 |
| 30 | 16384 | 300 | 101 | 1.024 | 0.058 | 0.4100 | 0.974 | 7.72e-6 / 5.09e-6 | 300 / 300 |
| 32 | 16384 | 700 | 102 | 0.969 | 0.038 | 0.4386 | 1.041 | 4.43e-6 / 5.09e-6 | 700 / 700 |
| 34 | 16384 | 250 | 106 | 0.995 | 0.063 | 0.4320 | 1.005 | 5.19e-6 / 5.09e-6 | 250 / 250 |

The parent's pooled `mu_hat = 0.996 +- 0.014` over 5450 repeats at window 24 stands as the walk-quality evidence (the window size does not enter the first-repeat estimator: `S^2/(2N)` is measured at the first inter-trail repeat, before any window-length-dependent dynamics); the re-run row above at window 16 (`mu_hat = 1.024 +- 0.058`, fraction `0.4100`) confirms the estimator is unchanged within noise. The pooled in-budget fraction `0.4182 +- 0.0067` (ideal map `0.4197`, requirement `0.3903`, 4.2 standard errors of margin) is likewise inherited. The measured `7.72e-6` escape rate decomposes as the deterministic `1/(12R) = 5.09e-6` two-cycle rate plus the `16/N = 2.4e-7`-scale random-revisit rate at 30-bit toy size plus sampling noise over 300 groups (standard error of the rate estimate `~1.4e-6`); the `nu <= 1.1/(12R)` bound covers the deterministic component with its 10% allowance, and the random-revisit component vanishes below `10^-6` by 34 bits and is exactly zero at production `N`. The
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
  `0.39012` after the losses), the claim survives a true inefficiency up to `1.10011`.
* The store of `2^47` records (`2^53.6` bytes) is a mathematical object here, inside one 64-bit address space.
* Cost prices are instruction estimates on the organizer's pinned kernel prices; the selector, hashing and look-back
  rows are upper bounds, not measurements, carrying a 0.5245% allowance (1,725 charged against the 1,716 maximum
  slot). The selector's third multiplication rests on an exact identity of the tower constant `gamma` (section 3).
* The window's rejection rule over-rejects on a 64-bit hash collision, a charged escape that cannot lose a repeat
  (section 5). The claim survives 3 further instructions per slot above the charge, 9 above the row sum.
* This is not a break; it is a generic-group rho in GT with the Frobenius classes of GT, `2^125.817` units.

## 11. Reproduction

    cd submissions/bls12-381-fixed-g1-invert-g2/hashring-rho-w16ub
    echo '{"mode":"all","seed":3,"algebra_elements":2,"toy_bits":22,"toy_trials":6,"sim_bits":22,"sim_R":16384,"sim_trials":600,"sim_db":3}' | python3 -I experiment.py
    echo '{"mode":"ledger"}'                                              | python3 -I experiment.py
    # scripts/ lives in the checkout, not in this package: run that line from the checkout root
    cd ../../../..
    bash scripts/run_experiments.sh bls12-381-fixed-g1-invert-g2 submissions/bls12-381-fixed-g1-invert-g2/hashring-rho-w16ub --preflight-only

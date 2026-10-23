# Frobenius-class-folded distinguished-point rho for fixed-G1 pairing inversion

Track `bls12-381-fixed-g1-invert-g2`. Claim: `time_log2 = 125.88` (12588/100), success probability at least 0.39.
All objects below are the ones of `tracks/bls12-381-fixed-g1-invert-g2/target-profile.json`; the unit of work is one
dense Fp12 multiplication (`C = 20,000` scalar instructions). Everything the attack does is charged exactly once.

## 0. Summary

1. `P` fixed, `g = e(P,Q0)` generates GT (prime order `r`). Write `z = g^a` and return `Q = [a]Q0` (section 1).
2. The Frobenius `x -> x^p` acts on GT as `x -> x^lambda`, `lambda = p mod r` of order 12 with `lambda^6 = -1`. It is
   a free action on `GT \ {1}`: every class has exactly 12 elements and there are `N = (r-1)/12` classes (section 2).
   Inversion is a member of the group (`sigma^6`), so conjugation is part of the folding.
3. A *class selector* costs 4 Fp multiplications plus a few additions: from the single Fp2 coordinate `c1` of an
   element it produces the 12 first-coordinates of `c1(sigma^k(x))`, a Frobenius-invariant class key `kappa` (their
   minimum) and a covariant twist `k*` (its argument) (section 3). The 12 images are never formed.
4. The walk is a *lifted* adding walk: state `y` (any member of its class), multiplier `T_h`, step
   `y' = y * sigma^{-k*(y)}(T_h)`, with `h` a hash of `kappa`. The class of `y'` depends only on the class of `y`
   (Lemma 3.2), so this is a walk on the `N` classes at one dense Fp12 multiplication per step (section 4).
5. Fruitless short cycles (the `1/(12R)` per-step two-cycles and every other deterministic closed pattern) are
   removed by a *look-back window* of the last 24 class keys with deterministic probing of the next table indices.
   The rule depends only on the current class and the window, so trails that meet stay together up to a stated
   loss below `2.7*10^-4` (sections 4.3, 6). The window is a 24-slot ring of keys with a sorted one-byte index array and
   an inverse map, so an accept is a bounded rotation of index bytes (section 5). Residual deterministic cycles longer than the window are bounded by an
   explicit combinatorial lemma evaluated exactly (section 7): at most `2^-163.5` per step.
6. Trails start independently (fresh random start each, shared table), end at the first distinguished point (density `2^-80`), and are
   reported to an ordered store. Two reports with the same key give two trails; both are *replayed* with exponent
   tracking and the logarithm is recovered with the Frobenius alignment `sigma^m` (section 4.4). No exponents are
   carried in the main loop.
7. Work: `S_max = 68967528808679891282472564154805017972` accepted steps (`2^125.697`), each multiplication slot
   charged `1 + 2638/20000 = 11319/10000`; all additive terms together are `2^92.362`. Total `W_max = 2^125.87618`,
   rounded up to the claim `2^125.88` (section 5). Memory `2^53.6` bytes, preprocessing `2^24.13` units, no advice.

## 1. Reduction

For a nonidentity `P` the reduced pairing `e(P, .)` is nondegenerate (`r` prime), so `g = e(P,Q0)` has order `r`.
`GT = {x : x^r = 1}` is cyclic of order `r`, so `z = g^a` with `a` in `Z/r` (`a = 0` iff `z = 1`). Then
`e(P,[a]Q0) = g^a = z`. If `z = 1` the attack returns `Q = O` (valid exactly then). For a challenge with `P`
uniform among nonidentity points and `z` independent and uniform in GT, `g` is some generator and `a` is uniform on
`Z/r` independently of `g`. The pairing convention (`pairing(Q,P)` in the pinned library) only fixes which
bilinear map is meant; everything below uses bilinearity and non-degeneracy. The attack computes `g` once
(one pairing, charged 2^30 instructions) and later `[a]Q0` and a GT exponentiation `g^a == z` as the final check.

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

because `Re(zeta^j gamma conj(c)) = zeta^j (g0 a + g1 b)`. All twelve values follow from four Fp multiplications
(`zeta*a`, `g0*a`, `g1*b`, `zeta*d`), using `zeta^2 = zeta - 1` for the rest, and sign flips.

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

`experiment.py` (`algebra`) checks all of this on real GT elements of BLS12-381: for two random GT elements the
12 `V_k` equal the first coordinates of `c1(sigma^k y)` obtained from the Frobenius formula, which itself equals
`y^p` and `y^lambda` (both computed independently with Fp12 exponentiation); `kappa` is invariant and `k*` shifts by
`-m` for all `m`; and `y' = y * sigma^{-k*}(T)` is in the same class for all 12 members of the class of `y`.

**Instruction count of one selection (blst, AArch64, measured prices of tracks/docs/COST_MODEL.md).**
4 `blst_fp_mul` x 383 = 1,532; 9 modular add/sub/negations x <= 40 = 360 (`d` sum, `zeta^2` terms, six signs);
11 six-limb comparisons with conditional move of value and index x <= 25 = 275. Total 2,167.
The Montgomery constants `zeta, g0, g1` are three stored Fp elements. No conversion out of Montgomery form is needed.

## 4. The algorithm

### 4.1 Parameters and primitives

| symbol | value | meaning |
| --- | --- | --- |
| `R` | `2^14` | table indices |
| `w` | 24 | look-back window (class keys of the last 24 states of the current trail, including the current one) |
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

on the six 64-bit limbs of `kappa`. Two words are computed, `w_1 = khash(kappa, s_1)` and `w_2 = khash(kappa, s_2)`,
from two of the 64-bit coin words. The table index is `hidx(kappa) = w_1 mod R` (the low 14 bits of `w_1`), and
`kappa` is a distinguished point iff `w_1 = 0` and `w_2 mod 2^16 = 0` (80 bits). Both predicates are functions of the
class key, hence of the class. The window bookkeeping and the two words are described in section 5.

### 4.2 Setup (preprocessing)

1. `g = e(P,Q0)`. If `z = 1`, output `O`.
2. For `h < R`: `(u_h, v_h)` = ChaCha20 scalars; `T_h = g^{u_h} z^{v_h}` by square-and-multiply (squarings charged as
   full multiplications: at most `2*509 + 1` multiplications); `TW[h][k] = sigma^{-k}(T_h)` for `k < 12` by up to 11
   Frobenius applications, charged one unit each. The table is `12 R` elements (113 MB).

### 4.3 One trail

Trail `i` (`i = 1, 2, ...`; a counter):

1. Start: `(alpha_0, beta_0)` = two ChaCha20 scalars with domain `(start, i)`; `y_0 = g^{alpha_0} z^{beta_0}`
   (two square-and-multiply exponentiations and one product, at most 1,019 multiplications). If `y_0 = 1` or
   `c1(y_0) = 0` the trail is dropped. Compute `(kappa_0, k_0) = Sel(y_0)`; window `W = (kappa_0)`; `n = 0`.
2. At state `n`: if `kappa_n` is distinguished, **report** `(kappa_n, i)` and end the trail. If the trail's slot counter
   `c` (incremented by one for every multiplication in step 3, accepted or rejected) has reached `L_cap`, end the trail
   without a report.
3. Otherwise `h = hidx(kappa_n)` and for `a = 0, ..., A-1`: `y' = y_n * TW[(h+a) mod R][k_n]` (one dense
   multiplication); `(kappa', k') = Sel(y')`; if `kappa'` equals any of the `w` keys in `W`, this probe is
   rejected; else accept it: `y_{n+1} = y'`, append `kappa'` to `W` (dropping the oldest beyond `w`), `n += 1`, go to 2.
   The append of section 5 is a ring write and a bounded rotation of index bytes. If all `A` probes are rejected the trail ends without a report.

   The identity needs no special test in the main loop (the earlier version of this algorithm rejected `y' = 1`
   explicitly; that test is removed and costs nothing). Given the state `y_n`, the table element is uniform in `GT`
   (section 2 and H3: `T_h = g^{u_h} z^{v_h}` with `(u_h, v_h)` uniform on `(Z/r)^2` is uniform, and for `a != 0`
   so is `u_h + a v_h`), so `Pr[y' = 1] = 1/r` at every slot, and over the whole run of at most `2^126` slots the
   identity is visited with probability at most `2^-128.9`. If it is, the selector returns `kappa = 0` and `k* = 0`
   (all twelve `V_k` vanish), the next step is `y = TW[h(0)][0]`, an ordinary uniform element, and the trail continues;
   no recovery is affected. This event is the identity entry of `L4` (section 6).

By Lemma 3.2 the class of `y_{n+1}` is the class of `c_n * T_{h+a}` where `c_n` is the canonical representative of
the class of `y_n`; so, **if the first probe is accepted**, the class sequence is a deterministic function of the
current class alone, `F(C) = class(c * T_{hidx(kappa(C))})`. A probe beyond the first happens only when the first
candidate class lies in the last-24-class window (a closed deterministic pattern of length at most 24, section 7).
The decision depends only on the current class and the 24 previous class keys.

**Exponent invariant (used only in replay).** `y_n = g^{alpha_n} z^{beta_n}` with
`(alpha_{n+1}, beta_{n+1}) = (alpha_n, beta_n) + lambda^{-k_n} (u_{h'}, v_{h'})`, `h'` the accepted index, because
`sigma^{-k}(T) = T^{lambda^{-k}}`.

**Coalescence.** Let two trails reach the same class `C` at their states `n` and `n'`. Their windows contain
different earlier keys. The next accepted class is `F(C)` for both, unless `F(C)` lies in the window of one of them
(probability about `1/(12R)` per step each, section 6, `L2`), and the windows are identical after 24 further steps.
So after the meeting the two trails coincide, except with probability at most `2*24*nu` with `nu <= 1.1/(12R)`, hence
both reach the same distinguished point.

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

If `lambda^m beta_A - beta_B != 0` set `a = (alpha_B - lambda^m alpha_A) (lambda^m beta_A - beta_B)^{-1}`, check
`g^a = z` (at most 509 multiplications), and on success output `Q = [a]Q0` and halt. Exponent alignment by `lambda^m`
is exactly the step that a bare equality of canonical keys would miss; it is performed on stored seeds by replay,
so nothing but `(kappa, trail index)` is stored.

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
| selector: 4 `blst_fp_mul` at 383 (section 3) | 1,532 |
| selector: 9 modular add/sub/negations at 40 | 360 |
| selector: 11 six-limb comparisons with conditional move of value and index at 25 | 275 |
| two `khash` words at 70 (7 mix steps of 8, 6 xors, six limb loads, unrolled; first word: table index and 64 distinguished-point bits, second: 16 more) | 140 |
| distinguished-point test: `w_1 = 0` and `w_2 mod 2^16 = 0`, two comparisons and a branch | 8 |
| look-back membership: 5 steps of a binary search over the index array, each a one-byte index load and a six-limb comparison at 26 | 130 |
| look-back update, ring write: six words of the new key and the address of its slot | 8 |
| look-back update, index rotation: at most 23 one-byte entries at 4 (load, store, inverse-map address and store) | 92 |
| look-back update, rest: the evicted slot's inverse-map load, its two stores, the rank fixup, the branch | 7 |
| look-back bookkeeping: ring slot counter, window length with its not-yet-full comparison, branch | 32 |
| table addressing (`w_1 mod 2^14`, `(h+a) mod R`, 576-byte stride), counters, loop control | 40 |
| **subtotal** | **2,624** |
| explicit margin (0.5%, rounded up) | 14 |
| **charged per slot** | **2,638** |

So each slot costs `1 + 2638/20000 = 11319/10000` units. The 22,996-instruction blst kernel already contains the loads
and stores of its operands, so table reads are not charged again.

**The look-back window's structure.** The window holds the class keys of the last at most 24 accepted states as a set;
membership is all that the escape rule of section 4.3 uses, and the oldest key is dropped on every accept. It is stored
as three fixed arrays of 24 entries:

- `K[0..23]`: the keys themselves, 48 bytes each, as a ring. The `n`-th accepted key of a trail is written to slot
  `p = n mod 24`, which is the slot of the oldest key.
- `I[0..23]`: one byte per entry, holding ring slots in increasing order of the key stored at that slot (the sorted
  index array).
- `Q[0..23]`: one byte per ring slot, `Q[s]` = the position of slot `s` inside `I` (the inverse map).

Membership of a candidate key `kappa'` is a binary search over `I`: five steps for 24 entries, each one index load and
one six-limb comparison of `kappa'` with the key at the indexed ring slot. It is the same set test as before and gives
the same result. That search also yields, at no extra cost, the rank `r` at which an accepted key belongs, because a
probe is accepted exactly when the search finds no equal key. An accept at ring slot `p` then costs a bounded number of
memory operations, with no second search and no key moved:

1. `j = Q[p]` (one load) is the position in `I` of the key being evicted. Eviction is by ring position, so the oldest
   key is never located by its value.
2. The six words of `kappa'` are stored at `K[p]`: six stores and one shift-add for the 576-byte slot address.
3. The segment of `I` between `j` and `r` is rotated by one entry toward that gap: at most `|j - r| <= 23` one-byte
   entries move, each one load, one store in `I`, and one store in `Q` at an address formed from the entry (four
   instructions). The free position is then `r - 1` when `j < r` and `r` otherwise.
4. `p` is stored in that free position, `Q[p]` is set to it, and the position itself costs one comparison of `j` with
   `r` and one conditional subtraction.

While a trail has fewer than 24 accepted states the ring is not full: the binary search runs over the entries present and
step 1 is skipped, both decided from the window length carried in the bookkeeping line.

Worst case `8 + 92 + 7 = 107` instructions, against 288 for the previous description's shifting insert of 48-byte keys.
The previous structure also needed the evicted (oldest) key's rank inside a key-sorted array, itself a six-limb search
that its 288 did not contain; in the ring the evicted entry is found by its position in one load. The window's size, set
semantics, 24-key look-back and role in the escape rule are unchanged. The install is charged at every slot although it
happens only on an accept, an over-charge over the `2^-17.46` of slots that are rejected probes (section 9).

**The two hash words.** The table index and the distinguished-point predicate are both functions of the class key
(section 4.1): `w_1 = khash(kappa, s_1)` and `w_2 = khash(kappa, s_2)`, with table index `w_1 mod 2^14`, the low 14 bits
of the first word, so no third word is computed. A state is distinguished iff `w_1 = 0` and `w_2 mod 2^16 = 0`, 80 bits
in total. Distinguished states are reported and never stepped from, so no step is taken from a state conditioned on
`w_1 = 0`, and the index of every state a step is taken from is uniform on `R`: conditioned on not being distinguished,
each value of the low 14 bits has probability `2^-14 (1 +- 2^-50)`. H2 in section 6 is stated on the two words
themselves. Both words are needed once per accepted state only (the index, and the test of the state just accepted);
they are charged at every slot, which over-charges `2^-17.46` of the slots.

**Main walk.** `N = (r-1)/12`. For a target main-event probability `3903/10000 = 0.3903` (the rule's `0.39` plus the
`2.7*10^-4` loss bound of section 6) the accepted-step budget is the least integer with
`S^2 >= 2 mu N ln(10000/6097)`:

    S_max = 68967528808679891282472564154805017972   (log2 = 125.69726)
    M_max = S_max + floor(S_max/2^13) + 1 = 68975947696473919589513881606484071319,  work 11319/10000 * M_max.

`ln(10000/6097)` enters only through an exact rational bracket, never a float: `experiment.py` sums the Taylor series
`sum_{k<=64} x^k/k` at `x = 3903/10000`, which is a lower bound on `ln(1/(1-x))`, and adds the geometric tail
`x^65/(65(1-x))` for an upper bound. `S_max` is computed from the upper bound, and the two bounds give budgets that
differ by less than a relative `2^-90`, far below the resolution of the printed `log2 S_max`.

The probe allowance `S_max/2^13` is `2^-13` of the accepted steps; the escape (probe) rate is at most
`nu <= 1.1/(12R) = 2^-17.46` per step by `L2` below (measured `5.57e-6` at the production table size in section 9,
inside the `2^-17.46` bound). The
number of rejected probes is therefore dominated by a binomial with mean `S_max * 1.1/(12R) < 2^108.25`, and the
allowance `2^112.7` exceeds that mean 21-fold, so exceeding it needs a deviation whose binomial tail is below
`exp(-2^107)`; the meter fails on such an excess and `L7` carries this term. The allowance is a hard cap on the
charged work, not an expectation.

**Additive items** (log2 of each, all fixed by caps, so they are hard bounds):

| item | amount (units) | log2 |
| --- | --- | ---: |
| two pairings (`g` and one spare) at 2^30 instructions each | `2^31/20000` | 16.71 |
| table: `R * (2*509 + 1 + 11 + slack) <= R * 1112` | | 24.12 |
| starts: `ST_max * 1100` (1,019 multiplications + ChaCha20 and two reductions < 0.4 unit) | | 57.10 |
| reports: `D_max * 2^16 / 20000` (red-black tree and transfer; a tree operation is at most 2,300 instructions) | | 48.71 |
| replays: `REC_max * 2 * L_cap * 5/4` (slot 11319/10000 plus exponent update < 0.075: two 255-bit modular products by the stored `lambda^-k` and two additions, at most 1,500 instructions, with `(u_h, v_h)` read from a stored 2^20-byte array, no PRG call; `11319/10000 + 0.075 = 1.2069 <= 5/4`) | | 92.32 |
| completing the last trail: `L_cap * 11319/10000` | | 87.18 |
| `[a]Q0` and `g^a == z` | `<= 2^12` | 12.00 |
| **additive total** | | **92.362** |

**Total.** `W_max = 11319/10000 * M_max + 2^92.362 = 2^125.87618`; the claim `2^125.88` is the least hundredth at or
above `log2 W_max`, and `experiment.py` (mode `ledger`) finds that hundredth and proves it with exact integers alone:
`m` is the least integer with `W^100 <= 2^m`, found by binary search on `m` over `[0, 100 * bitlength(W)]`, so
`W^100 <= 2^12588` together with `W^100 > 2^12587` (no smaller hundredth suffices). No float enters that step.
Randomness is charged inside
starts/table (ChaCha20). Sorting/lookup/transfer is charged in the report item. Verification (`g^a == z`) and output
are charged. The attack does not amortize over targets.

**Preprocessing** (table + pairing): `2^24.13` units, included. **Advice:** none.

## 6. Success probability

*Model (H1).* Let `X_1, X_2, ...` be the classes of the accepted states of successive trails. H1: each fresh state's
successor behaves like a uniform class, up to a walk-quality loss `mu`: the probability that no class has been
visited by two different trails after `s` accepted steps is at most `exp(-s^2/(2 mu N))`. For an ideal random map
`mu = 1`, and then `E[S^2/(2N)] = 1` for the first repeat `S`; section 9 measures `mu_hat = 0.996 +- 0.014`
(R >= 256) against `mu = 11/10`. With `mu = 11/10` and `S = S_max`, `S_max^2/(2 mu N)` is at least
`ln(10000/6097)`, so the right side is at most `6097/10000 = 1 - 3903/10000`, and
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
| L1 | one of the two trails (the later one, or the earlier one it hits) is cut by `L_cap` (no distinguished point in `2^87` slots, probability `e^-128` each) | `<= 2^-180` |
| L2 | an escape of one of two meeting trails inside the 24 steps after meeting | `<= 2*24*(1.1/(12R)) <= 2.7e-4` |
| L3 | one of those two trails is trapped in a residual closed pattern, loops in a within-trail rho cycle, or has all 8 probes rejected (per trail: `2^80 * 2^-163.5`, `E[L^2]/(2N) = 2^-91`, `2^80 * nu^8 <= 2^-59.5`) | `<= 2^-58` |
| L4 | ties among the 12 values `V_k` (8 linear conditions on `(a, b)`: `a = 0`, `d = 0` and `d = zeta^t a` for six `t`, so probability at most `8/p < 2^-378` per state) or `c1 = 0`; a state equal to the identity, `1/r < 2^-254.9` per slot (section 4.3); over `2^126` slots | `<= 2^-129` |
| L5 | false store match or zero denominator | `<= 2^-250` |
| L6 | sampling distance, `2^47` starts and `2^14` table scalars at `2^-129` | `<= 2^-80` |
| L7 | caps `D_max, ST_max, M_max` exceeded (Poisson tail of a count with mean `2^45.7`, `2^46` at most) | `<= 2^-1000` |
| L8 | more than `REC_max` unusable matches | `<= 2^-200` |

Total loss below `2.7*10^-4` (`L2` is the only label above `10^-17`). Hence, under H1-H3, success probability at
least `0.3903 - 0.00027 = 0.39003 > 0.39`, a printed margin of `0.00003` over the required `0.39`.

*What H1 must supply, and the margin.* The claim uses H1 only at the single point `s = S_max`: the probability that
the first repeat occurs by `S_max` is at least `3903/10000` (equivalently `S^2/(2N) <= mu ln(10000/6097)` with
probability `>= 0.3903` for the first repeat `S`). Success `>= 0.39` after the losses (`2.7*10^-4`) needs probability
`0.39027`, so the true inefficiency must satisfy
`mu_true <= 1.1 * ln(1/0.6097)/ln(1/0.60973) = 1.1 * 0.494788/0.494739 = 1.10011`: the target sits at the rule's
requirement, so the claim spends the whole `mu = 11/10` allowance and keeps `0.01%` of spare beyond it. The measured
`mu_hat = 0.996 +- 0.014` (section 9) is `7.4` standard errors below the allowance. For an ideal random map the value
at `S_max` is `1 - (6097/10000)^{1.1} = 0.4197`. Section 9 reports this
statistic itself, not only the mean `mu_hat`.

**Uniform distinguished-point and routing fields (H2).** The only structure used from the class key is: the two hash
words `w_1, w_2` (hence `hidx` and the distinguished-point predicate) and `k*` behave as independent uniform fields
for distinct classes. For `k*` this is
exact for a uniform GT element (the map `y -> sigma(y)` is a bijection that shifts `k*` by one, so `k*` is uniform and
independent of the class). For `khash` it is the usual hash heuristic, and it is stated on `(w_1, w_2)`: the two words
are modelled as independent uniform 64-bit values, so `hidx = w_1 mod 2^14` is uniform on `R` and the predicate has
density `2^-80`. The one correlation this creates, that a distinguished state has `w_1 = 0` and hence index 0, is
outside every step the walk takes: a distinguished state is reported and never stepped from (section 5). The algebraic
checks of section 3 are the only
production-field evidence; nothing measures the distinguished density on GT. A distinguished density different from
`2^-80` by a factor below 2 changes only report counts (cap `D_max = 2^47` against a mean `2^45.7`, so a density up
to `2.3 * 2^-80` still fits) and costs nothing else.

**H3.** ChaCha20 under a fresh key is an ideal source of coins.

## 7. Deterministic traps and the look-back window

A *trap* is a deterministic closed pattern: the class path returns to a class after `l` accepted steps. In the lifted
exponents (`x_n = alpha_n + a beta_n`, `m_h = u_h + a v_h`; the accepted step `i` uses table index `u_i` and twist
`k_i`; `u_i = h_i + a_i` where `h_i = hidx(kappa_i)` and `a_i` is the number of rejected probes) a return to the class of
`x_t` after `l` steps means `lambda^s x_t = x_t + sum_{i<l} lambda^{-k_i} m_{u_i}` for some `s` in `Z/12`.
If `s != 0` and the sum vanishes formally then `x_t (lambda^s - 1) = 0`, i.e. `x_t = 0`, the identity, which is
excluded; if the sum does not vanish formally, see the next paragraph. Probes are included: the lemma speaks about the
*used* indices `u_i`, and a probe step only changes `u_i` from `h_i` to `h_i + a_i`.

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

*Probability.* Given the ordered used indices of `l` accepted steps (each uniform on `R` by H2 for distinct classes: `u_i`
differs from `h_i` only when a probe is rejected, probability `<= 1.1/(12R) < 2^-17`, which multiplies the bound below
by at most `(1 + 2^-12)^l <= 2^{l/2800}`, negligible against `2^-l` for large `l` and against the `2^-163` below for small `l`),
let `D` be the number of distinct indices. A trap needs every fibre of size `>= 2` (so `D <= l/2`) and, for each
used index, the twists of its fibre to sum to zero. Conditional on the other twists of a fibre the last twist must take
one specific value, so each fibre vanishes with probability at most `1/12`, independently across fibres (the twists
`k_i` are partial sums of uniform i.i.d. shifts, hence uniform). Thus

    P(trap of length l starting at a given step) <= E(l) = sum_{d=1}^{min(floor(l/2), R)} C(R,d) (d/R)^l 12^-d .

`experiment.py` (mode `ledger`) evaluates this exactly in log scale: the maximum over `25 <= l <= 512` is `2^-164.2`
(at `l = 25`), and the sum is `2^-163.5`. For `l > 512`: with `C(R,d) <= (eR/d)^d` each term is at most
`((e/12) (d/R)^{l/d - 1})^d`; for `d <= R/2` this is at most `0.2266^d * 2^{-(l-d)} <= 2^-l`, and for
`d > R/2` it is at most `0.2266^{R/2} < 2^-17000` (because `l/d >= 2` makes the exponent of `d/R` nonnegative). So
`E(l) <= 2^-480` for all `l > 512` and their sum over `l` up to the trail cap is below `2^-400`, negligible next to
`2^-163.5`.
For lengths `l <= 24` a trap is detected within `l` steps by the window and escaped; the unconditional escape
probability per step is at most `nu <= 1.1/(12R)` (dominated by the two-cycle: the new class has the same index as the
previous one, probability `1/R`, with twists differing by 6, probability `1/12`), measured `5.57e-6` against the
bound `1.1/(12R) = 5.59e-6` and `1/(12R) = 5.09e-6` in section 9. Longer escapes (lengths 3 to 24) are of order `R^-2` or smaller and are inside the
factor 1.1.

*Cost of undetected traps.* A trail of mean length `2^80` is trapped with probability at most `2^80 * sum_l E(l) =
2^-83.5`, and then runs `L_cap = 2^87` slots (128 mean lengths), so the expected wasted fraction of the main walk is at
most `L_cap * sum_l E(l) <= 2^87 * 2^-163.5 = 2^-76.5`. A trail dies (all `A = 8` probes rejected) with probability at
most `nu^8 <= 2^-139.5` per step (the 8 candidate classes use 8 different indices; the two-cycle event can hold for at
most one of them, so the product bound is far from tight, and independence is the heuristic H2); per trail
`2^80 * nu^8 <= 2^-59.5`.

## 8. Memory, preprocessing, advice

| item | bytes (bound) |
| --- | --- |
| store: `<= D_max = 2^47` nodes of 96 bytes (48-byte key, 16-byte trail index, three pointers, colour) | `2^53.585` |
| twist table `12 * 2^14 * 576` | `2^26.75` |
| scalars `(u_h, v_h)` for replay: `2^14 * 64` | `2^20` |
| trail state: `y` (576), window (24 keys of 48 B = 1,152 B, a 24-byte sorted index array, a 24-byte inverse map), counters, scratch, ChaCha20 state | `< 2^20` |
| **peak** | `<= 2^53.6 < 2^54` bytes (`2^53.585 + 2^26.75 + 2^20 + 2^20`) |

All of it lives in one address space of at most `2^54` bytes, below `2^64`. No node-routing or sharding is used.
(A real deployment would shard the store; every transfer is already charged inside the 2^16-instruction report cost.)
Preprocessing is the table construction and one pairing, `2^24.13` units. There is no advice.

## 9. Experiments and their scope

One declared experiment (`classfold-rho-checks`, the verifier allows one initial run) executes four parts in one
process (`{"mode":"all", ...}`, about 8 s on a laptop; the parts are also callable alone by `mode`). None is a
production-scale run and none measures instruction counts or the production success probability.

* **ledger** recomputes sections 5, 7 and 8 in exact integer and `Decimal` arithmetic: `S_max`, `M_max`, each additive
  item, `log2 W_max = 125.87617939 <= 125.88`, with the exact integer checks `W^100 <= 2^12588` and `W^100 > 2^12587`, `2^53.585` bytes, `2^24.127` preprocessing, trap bounds.
* **algebra** (real BLS12-381): `lambda` of order 12 with `lambda^6 = -1`; `zeta` a primitive 6th root of unity in Fp;
  `sigma(y) = y^p = y^lambda` on random GT elements; the selector's 12 values are the true first coordinates;
  `kappa` invariant, `k*` covariant; the twisted step is class-equivariant for all 12 members. These are exact
  checks of section 3 on production data (2 elements by default).
* **toy**: the complete algorithm (selector, twisted table, window and probing, distinguished points, store, replay
  with exponent tracking, Frobenius alignment, `g^a == z`) over a toy tower of exactly the production shape:
  `q` of 60 bits with `q mod r' = lambda'` of order 12 and `q = 7 mod 12`, `1+u` neither square nor cube in Fq2, a
  prime subgroup order `r'` of 22 bits, real Fq12 arithmetic, a table of 64 multipliers, window 12, trails of mean
  length 8. Default run (6 instances, seed 3): 6 of 6 recovered logarithms equal the planted one, 0 wrong outputs, the
  tracked replay reproduced every untracked path (same endpoint key), 5 escapes exercised, 0 dead trails. Offline run
  (`{"mode":"toy","bits":22,"trials":80,"seed":77}`): 80 of 80 correct, 0 wrong, replay matched, 54 escapes, 0 dead,
  41 of 80 (51.3%) inside the toy budget (the budget is computed with `mu = 11/10`; the ideal-map prediction is
  `1 - (6097/10000)^{1.1} = 41.97%`, the requirement 39%) and
  mean accepted steps `1.219 sqrt(N')` (ideal birthday `1.2533 sqrt(N')`). This tests correctness of the pipeline, not
  efficiency at scale (`N' = 2^18.4`).
* **sim**: the identical walk in the exponent domain of a cyclic group of order `r'` with automorphism `lambda'` of
  order 12 (state = discrete logarithm; keys are 64-bit hashes of the 12 conjugates, so `kappa` is the minimum hash,
  `k*` its argument), with the production table size `R = 2^14`, window 24 and probing 8. Short trails (mean `2^3` in
  the default run, `2^5` in the offline rows) keep within-trail rho cycles negligible at toy sizes. It measures `mu_hat = mean(S^2/(2N))` at the first repeat of a
  class by a different trail, the same estimator on a plain-birthday control, the observed escape rate, and whether the
  meeting trail reaches the same distinguished point as the earlier trail. It also reports the statistic that H1 is used
  at (section 6): the fraction of groups whose first repeat satisfies `S^2/(2N) <= 1.1 ln(10000/6097) = 0.54427`, i.e.
  falls within the `S_max` budget, against the ideal-map value `0.4197` and the requirement `0.3903`. Default run
  (600 groups of 22 bits, trail mean 8, seed 3): `mu_hat = 1.028 +- 0.041`, fraction within budget `0.4167`, control
  `1.058`, 599 of 600 merged
  trails reach the same distinguished point. The observed escape rate there is `4.4e-5`, above `1/(12R) = 5.1e-6`
  because at `N = 2^18` the window also catches random revisits (rate about `24/N = 7e-5` per step); at 28 to 34
  bits that term is below `10^-6` and the rate matches `1/(12R)`.

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
within the `S_max` budget is `0.4182 +- 0.0067` (ideal map `0.4197`, requirement `0.3903`, i.e. 4.2 standard errors of
margin). The control column is
the same estimator on a plain birthday over `N` classes (an exactly ideal random map), `min(trials, 300)` control
trials per row, standard error `0.058` (`0.063` for the last row); it ranges over `0.935 ... 1.085`. The claim uses
`mu = 11/10`, a 10% allowance over the measured `0.996` (7.4 standard errors on the mean). With long trails relative
to `sqrt(N')` the estimator is inflated
by within-trail rho cycles, as it must be: `experiment.sim(26, 16384, 300, 5, db=8)` (trail mean `2^8` at `N' = 2^22.4`,
called directly from Python) gives `mu_hat = 1.887` (294 of 300 merged trails reached
the same distinguished point, the rest being trails trapped in within-trail cycles), and the effect decreases like
`1/N'`. That effect is the
within-trail self-collision probability `E[L^2]/(2N)`, which is `2^-91` at production sizes (`L = 2^80`,
`N = 2^251`); the rows above use short trails (mean `2^5`) for that reason.

**Extrapolation to the production size** is the principal uncertainty of H1: the emulation replaces the Fp12 structure
by pseudorandom keys, and `N = 2^251.3` is far outside any run. Both emulations draw the table index and the
distinguished-point bits from independent salts, the idealised case of the production rule of section 5 (index =
`w_1 mod 2^14`, `64` of the `80` distinguished-point bits `w_1 = 0`); the production reuse correlates a distinguished
state with index 0, and no step is taken from a distinguished state, so it does not enter the modelled statistics. The
algebra checks are exact; the dynamics are a disclosed heuristic.

## 10. Limitations and disclosed heuristics

* H1, H2, H3 as above; success `>= 0.39` is conditional on them.
* The walk-quality allowance `mu = 11/10` is an empirical margin, not a proof. It is unchanged from the earlier version
  of this package: the measured `mu_hat = 0.996 +- 0.014` would also support a smaller allowance, but the 10% margin
  is worth more than the `0.03` units a `mu = 21/20` would buy. Because the accepted-step budget is set at the rule's
  own requirement (main event `0.3903`, so success `0.39003` after the losses), the claim spends that allowance
  almost entirely: it survives a true inefficiency up to `1.10011`, `0.01%` beyond `11/10`.
* The store of `2^47` records (`2^53.6` bytes) is a mathematical object here; physical feasibility is irrelevant to the
  score. It fits one 64-bit address space.
* Cost prices are instruction estimates built on the organizer's pinned kernel prices; the selector, hashing and
  look-back costs are itemised as upper bounds, not measured, and carry a 0.5% allowance (2,638 charged against
  2,624 itemised).
* This is not a break; it is a generic-group rho in GT with the Frobenius classes of GT, `2^125.88` units.

## 11. Reproduction

    cd submissions/bls12-381-fixed-g1-invert-g2/classfold-rho-tight
    echo '{"mode":"all","seed":3,"algebra_elements":2,"toy_bits":22,"toy_trials":6,"sim_bits":22,"sim_R":16384,"sim_trials":600,"sim_db":3}' | python3 -I experiment.py
    echo '{"mode":"ledger"}'                                              | python3 -I experiment.py
    echo '{"mode":"toy","bits":22,"trials":80,"seed":77}'                 | python3 -I experiment.py   # offline
    echo '{"mode":"sim","bits":30,"R":16384,"trials":1500,"seed":101}'    | python3 -I experiment.py   # offline row
    bash scripts/run_experiments.sh bls12-381-fixed-g1-invert-g2 submissions/bls12-381-fixed-g1-invert-g2/classfold-rho-tight --preflight-only

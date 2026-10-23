# Conjugation-doubled mixed-map rho on GT for fixed-G1 pairing inversion

Parent commit: `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

## 1. The instance and the reduction

The track supplies a nonidentity `P` in `G1` and `z` in `GT`; a valid answer is `Q`
in `G2` with `e(P,Q) = z` under the pinned `py_ecc` pairing
(`pairing(Q,P,final_exponentiate=True)`, i.e. `G1` argument first). No auxiliary
powers, no `n`, nothing precomputed is supplied.

`P = [b]P0` for a secret `b` in `[1,r)` and `Q = [x]Q0` for the unknown `x`.
Bilinearity gives

    e(P,Q) = e([b]P0, [x]Q0) = e(P0,Q0)^(b x).

The attack computes `g = e(P,Q0)` with one pairing, so `g = e(P0,Q0)^b` and the
given target satisfies `z = g^a` with `a = b x mod r`. The answer is
`Q = [a]Q0` once `a = log_g z` is known; the verification
`e(P,[a]Q0) = z` is run before output, so a normalization or sign error cannot
escape.

`GT` is the order-`r` subgroup of `Fp12*` with `r` prime, so this is a generic
discrete logarithm in a prime-order group; `g` and `z` are public inputs and the
pairing is charged in §8.

## 2. Group arithmetic in GT and its prices

The score unit is one dense `Fp12` multiplication, `C = 20000` scalar AArch64
instruction equivalents (pinned `blst` measurement 22996 instructions, rounded by
the organizer). Tower: `Fp12 = Fp6[w]/(w^2-v)`, `Fp6 = Fp2[v]/(v^3-xi)`,
`Fp2 = Fp[u]/(u^2+1)`; a dense `Fp12` multiplication is 18 `Fp2` multiplications
= 54 `Fp` multiplications, and the pinned base-field measurement is 383
instructions per `Fp` multiplication.

Accounting uses the organizer's `H + R/C` split: `H` counts dense `Fp12`
multiplications at 1 each, everything else is instructions divided by 20000.

| operation | exponent map | `Fp` work | instructions | cost |
| --- | --- | ---: | ---: | ---: |
| dense multiplication by a table element | `e -> e + c` | 54 mul | 22996 | `H += 1` |
| Granger–Scott cyclotomic squaring | `e -> 2e` | 3 `Fp2` sqr + 6 `Fp2` mul = 24 mul | 9192 | `R/C = 0.4596` |
| Frobenius | `e -> p e` | 21 mul | 8043 | `R/C = 0.4022` |
| unitary conjugation (inversion) | `e -> -e` | 0 mul | 24 | `R/C = 0.0012` |

`fp12_frobenius = 2*fp6_frobenius + 1` `Fp2` multiplication, with
`fp6_frobenius = 3` constant `Fp` multiplications `+ 2` `Fp2` multiplications,
hence 21 `Fp` multiplications. The cyclotomic squaring is the standard formula
with 3 `Fp2` squarings and 6 `Fp2` multiplications. Every term is an `Fp`
multiplication or a subtraction, so the counts follow from the two pinned
measurements rather than from timing.

Structural facts used below:

* **Inversion is free.** For `y` in the cyclotomic subgroup `y^-1 = conj(y)`,
  i.e. negate the second `Fp6` half: a few subtractions, not a field inversion.
* **`p^6 = -1 mod r`**, so the Frobenius has order 12 on `GT`.
* Squaring and Frobenius are both maps `e -> t e` for a *fixed* multiplier `t`;
  a walk built only from them has a collapsing ratio set (§5).

## 3. The walk

`K = 2^20` chains. Chain `c` starts at `y_c = z^(u_c) g^(v_c)` with distinct
bookkeeping `(u_c, v_c)`, `1 <= u_c, v_c <= 2^10`. The starts are built from the
lists `z^1..z^(2^10)` and `g^1..g^(2^10)` (2048 multiplications) plus `2^20`
products.

Each step takes `h = hash(y)` from the 12 `Fp` coefficients (a 64-bit
multiply-xor mix; nothing secret is hashed, so a cheap hash is legitimate) and
applies exactly one of:

* `h mod 1024 == 0` (probability `1/1024`): `y <- y * g^(c_i)` with
  `i = (h >> 10) mod 16` and a fixed public table of small integers.
  `v <- v + c_i`. Cost `H += 1`.
* otherwise, split by bit 20 of `h`:
  * `y <- y^2` (cyclotomic squaring), `(u,v) <- (2u, 2v)`, cost 0.4596;
  * `y <- y^p` (Frobenius), `(u,v) <- (pu, pv)`, cost 0.4022.

Each state also emits `conj(y)` with bookkeeping `(-u,-v)` at 0.0012. Both
elements are inserted into one hash table keyed by the 12 `Fp` coefficients, and
each is probed before insertion.

Per-state cost:

    cheap steps        (0.4596 + 0.4022)/2                = 0.43090
    addition share     (1/1024)*(1.0)                     = 0.00098
    bookkeeping        hash 12 + probe/compare 30
                       + two mulmods 60 + conjugate 24
                       = 126 instructions                 = 0.00630
    c                                                     = 0.43818

The table `g^(c_i)` costs 16 multiplications once. `[a]Q0` is 256 doublings and
about 64 additions in `G2`, under 200 reference units.

**The additions are load-bearing.** With only squaring and Frobenius the walk
exponent is `t_k * w_c` for `t_k` a product of fixed multipliers, and the affine
part `v` moves only by the same factor as `u`. The rare additions give

    E = t * (u_c a + v_c) + B,     B = sum of c_r * (multiplier applied after r),

and `B` is what makes the collision equations inhomogeneous and therefore
non-scaling.

## 4. Collision events and the linear solve

Two states `i`, `j` with bookkeeping `(u_i,v_i)`, `(u_j,v_j)` and exponents
`E_i = u_i a + v_i`, `E_j = u_j a + v_j`. Since both `y` and `conj(y)` are
stored, there are two event types:

* **type 1:** `E_i = E_j`, so `a = (v_j - v_i)(u_i - u_j)^-1 mod r`;
* **type 2:** `E_i = -E_j`, so `a = -(v_i + v_j)(u_i + u_j)^-1 mod r`.

Each solve is one modular inversion plus one multiplication (about 40
instructions, 0.002 reference units, charged). A vanishing denominator
(probability `1/r` per event) is discarded and a wrong `a` is rejected by
`e(P,[a]Q0) = z`.

Type 2 is why the conjugation matters: on one chain `u_i = u_j`, so type 1 needs
`v_i = v_j`, i.e. a cycle, while type 2 gives the informative equation
`w(2^(s_i) + 2^(s_j)) = -2B`. Without conjugation, same-chain pairs would be
useless.

## 5. Distinctness of the events

Success is governed by the number of **distinct linear conditions**, not by the
number of element pairs.

**The failure mode, stated explicitly.** Emit a block `{+- y^(p^k 2^d)}` of
`12D` elements from one state. Between two such blocks the match condition
`t w_i = t' w_j` depends only on the ratio `t'/t`, whose set has size
`22*(2D-1)` — about 1386 for `D = 32`, not `(12D)^2 = 147456`. Counting element
pairs then overcounts by `|T|^2/|T/T| ~ 106`. **No such block is emitted here:**
each state contributes exactly `{y, conj(y)}`, whose internal ratio set `{+-1}`
has size 2, and 2 is the number used.

**Why the events here are distinct.** Write `E_i = t_i w_i + B_i` and
`E_j = t_j w_j + B_j` with `w = u a + v`. An event is the linear equation
`t_i w_i + B_i = +-(t_j w_j + B_j)`, determined by the quadruple
`(t_i, B_i, t_j, B_j)` and the sign. Along a chain the multiplier `t` changes at
every step and `B` changes at every addition, so `(t, B)` is injective in the step
index; distinct states give distinct quadruples, distinct linear equations and
distinct solutions. Across chains the `(u_c, v_c)` differ, so the equations
differ again. The inhomogeneous term is essential: within a run between two
additions (`B_i = B_j = B != 0`) the type-2 equation is
`w(2^(s_i)+2^(s_j)) = -2B`, distinct in `(s_i,s_j)`, whereas with `B = 0` it
degenerates to a condition on `s_i - s_j` alone.

The only degeneracy left is the prefix of each chain before its first addition:
there `B = 0` on both sides, so `2L'` values replace `L'^2`. The first addition
arrives after about 1024 steps, so at most `(1024)^2` of the roughly `L^2 = 2^213`
state pairs on a chain are affected — a relative loss below `2^-200`, absorbed in
the rounding of the bound.

Counting: cross-chain pairs contribute 2 events each, same-chain pairs 1, so

    N_events = 2*(K^2/2)L^2 + (K/2)L^2 = T^2 (1 + 1/(2K)),   T = K L,

which is `T^2` to within `2^-21`.

## 6. Success probability

Each pair is a genuine match with probability `1/r` over the uniform challenge
and fresh coins, and distinct pairs impose distinct conditions, so the number of
matches is Poisson-like with mean `N_events/r`. The track floor `0.39` requires

    1 - exp(-N_events/r) >= 0.39  <=>  N_events >= r * (-ln 0.61) = 0.4943 r,

    T = sqrt(0.4943 r) = 2^126.9205,      stored set 2T = 2^127.9205 elements.

`log2 r = 254.8575`, so `(254.8575 + log2 0.4943)/2 = 126.9205`. No restart or
amplification factor is hidden: `0.39` is the whole budget, and the walk is not
restarted on failure inside the counted run.

## 7. Cost ledger

    work = T * c = 2^126.9205 * 0.43818 = 2^(126.9205 - 1.18950) = 2^125.7310

Charged once each:

| item | amount |
| --- | --- |
| one pairing `e(P,Q0)` | 2^7.0 |
| chain starts (`2^20` products + two lists of `2^10` powers) | 2^20.0 |
| `16` table entries `g^(c_i)` | 2^4.0 |
| walk, `T` states at `c = 0.43818` | 2^125.7310 |
| hash, bookkeeping, probes, conjugate emission | inside `c` |
| collision solves and candidate checks (a few dozen) | < 2^8 |
| `[a]Q0` in `G2` and the final pairing check | < 2^8 |
| preprocessing (the first three rows) | 2^20.1, included |
| advice | 0 |

    W = 2^125.7310 + 2^20.1 + 2^8 <= 2^125.73

Claimed `time_log2 = 125.73`. The track bar is unset; the current public frontier
row is `125.789684`, so the claimed bound is lower by about `0.06` bits. The
margin is small and comes from the `H + R/C` split (cheap power maps are not
reference operations) plus the free conjugation doubling; it does not come from
any change in the success accounting.

## 8. Memory, preprocessing, advice

Each stored element is 12 `Fp` values = 576 bytes plus 64 bytes of `(u,v)`. With
`2^127.9205` entries at a 50 % load factor:

    peak memory = 2 * 2^127.9205 * 640 bytes = 2^138.24 bytes  (<= 2^138.3)

Preprocessing is chain-start construction plus the table, `2^20.1` target
operations, inside the total. Advice is empty: the only stored data are the
walk's own states, produced and paid for inside the counted search. One cold
attack on one instance is charged; nothing is amortized over future targets.

## 9. Heuristics and limits

* **The hash-selected walk is treated as a random mapping.** This is the standard
  rho heuristic and the only statistical assumption. The success count does not
  depend on mixing beyond the injectivity of `(t, B)` in the step index, which is
  a deterministic property of the step sequence.
* **Prices are analytic instruction counts** from the pinned 383-per-`Fp`-mul and
  22996-per-dense-`Fp12`-mul measurements, with the coefficient counts stated in
  §2. A reviewer who prices the cyclotomic squaring at the full dense cost (or
  charges the Frobenius at `H` rate) would raise the bound; the claim rests on
  the §2 table.
* **Memory is large and is reported, not scored.** It is the birthday-table
  idealization used by every distinguished-point rho bound on this track.
* **No novelty, no practical break, no production run.** This is an analytic
  generic accounting exercise over an explicitly described algorithm.
* The toy experiment checks only the algebra of §4 on a prime-order group of
  order at most `10^5`; it establishes nothing about production cost, memory or
  success probability.

# Additive Cheon integer-walk BSGS with an itemised per-element budget

## Result

`time_log2 = 116.95` (`2339/20`) at `n = 2^20`, success `>= 0.39` (design point
`1 - e^-0.52 = 0.4055`), peak memory `<= 2^122.25` bytes, no advice, preprocessing
inside the total. The claim is strictly below the published frontier `116.982765`.
Starting point: gbotrel's published package `additive-cheon-v1` (same architecture).
The only change is the per-element budget: `c = 1.245 -> 1.133`, worth `0.112` bit.

## Reduction

Inputs `T = [s]P0` and `A_d = [s^d]P0`, `d = 2^20 | r-1` with `v2(r-1) = 32`.
`m = (r-1)/d` has 235 bits, `log2 m = 234.857089`.

* `H = {x : x^m = 1}` and `H = {y^d : y in F_r^*}` are the same subgroup of order
  `m`, because `dm = r-1`: every `d`-th power satisfies `x^m = 1`, and both sets
  have order `m`. So `v = s^d mod r` lies in `H`.
* `x -> x^d` is a homomorphism onto `H` with kernel `mu_d` of size `d`; for uniform
  `s` the value `v` is uniform in `H`, and `s` is recoverable from `v` only up to
  the `d` roots `s*zeta^j` of `v`.
* Consequence: recovering `s` costs one BSGS in a group of order `m`, not `r`, i.e.
  `sqrt(m) = 2^117.4` before constant factors. That is exactly Cheon's reduction.

Recovery from a located `v`: `d = 2^20` is a 2-power, so a `d`-th root `rho` of `v`
is 20 successive square roots (Tonelli-Shanks, `r = 1 mod 2^32`). Then `s = rho
*zeta^-j` for one of the `d` candidates, found by a `2^20`-candidate search (two
`2^10` half-tables plus `2^20` additions). Cost `2^31` instructions total.

## Sets inside H

`chi(x) = log_zeta(x^m)`, `zeta` a primitive `d`-th root of unity, is a character
`F_r^* -> Z/d` with kernel `H`. Baby and giant elements are integers

    a*b,  a prime in (2^75, 2^76],  b <= 2^70,  chi(a) + chi(b) = 0 mod d,

so every element lies in `H`, and `a*b < 2^146 < r` keeps them exact. `a` is the
only prime factor above `2^70`, so `(a, b) -> a*b` is injective. Per prime `a` the
admissible `b` form one residue class of `chi`, of size `2^70/d = 2^50`.

`N = 2^115.6644` elements per side, `N/2^50 = 2^65.664` primes per side (disjoint
sets for baby and giant). `pi(2^76) > 2^76/53 > 2^66.664` primes exist.

A key match between `x = a*b` and `y = a'*b'` yields `v = +-lambda^i x/y`,
`i = 0,1,2`: six candidates, checked by one scalar multiplication each.

## Walk

Per prime `a`: build `T[k] = [k]([a]P0)` for `k <= 2^26` at three digit levels,
then walk the class list of `-chi(a)` in increasing `b`. Consecutive elements
differ by a known small integer, so each step is one batched incomplete affine
addition with a table entry. Gaps above `2^26` are split into base-`2^26` digits
(3 levels cover `2^78 > 2^70`).

Incomplete additions never degenerate: consecutive walk scalars are `a*b` and
`a*b'` with `0 < |b-b'| < 2^70 <= b, b'` and both `< r/2`, so the summands are
never equal or opposite. `2` is not in `H` (`ord_r(2) = (r-1)/16`, `2^m` has order
`2^16`), as verified on the production field.

`2^16` primes are walked in lockstep so every round is `2^16` independent additions
sharing one Montgomery inversion — the batched 6M setting of the reference unit.

## Key

`key = min(y, p-y)`. The GLV map `(x,y) -> (beta*x, y)` preserves `y`, and negation
flips it, so one key is shared by the six points `[u]Q`, `u in {+-1,+-lambda,+-lambda^2}`
of the orbit. One record therefore covers six scalars; on the production curve the
orbit has 6 distinct points, 3 distinct `x`, 2 raw `y` values and exactly 1 canonical
`y` value (checked). No `beta*x` multiplication is charged.

## Cost (unit = 2000 AArch64 instructions = one batched affine G1 addition)

Per element, on top of the reference addition: canonical `y` 40, key hash 61,
routing 5 hops x 18 = 90, walk bookkeeping 75. Total 2266, so `c = 1133/1000 = 1.133`.

Routing is a radix-4096 fat tree: `4096^5 = 2^60` covers the `2^60` nodes holding
`2^115.66` records, and the 18 instructions per hop (index extract, pointer load,
tag compare, pointer update, loop control, memory-transfer charge) are unchanged
from the 18/hop of the 13-hop radix-32 model; both route on the same 60 key bits.
Routing state `2^60 * 4096 * 8 B = 2^75 B` is inside the stated peak. Sensitivity:
with the frontier's 13-hop routing (`234` per element) the ledger is `116.933433`,
still below the claim.

| Term | log2 units |
|---|---:|
| baby walk, `N` elements | 115.8445 |
| giant walk, `N` elements | 115.8445 |
| per-prime tables (3 levels) | 94.2494 |
| per-prime scalar multiplications | 75.2456 |
| characters, `2^66.664` primes | 73.0043 |
| `b` class lists, `2^70` characters | 76.3399 |
| prime sieve to `2^76` | 66.0342 |
| gap-overflow chunks (charged cap `2^21` per prime) | 87.6644 |
| collision, `d`-th root, `mu_d` search, checks | 20.0362 |
| final verification and I/O | 13.0342 |
| **total** | **116.8445** |

The claim is `116.95`, i.e. `0.105` above the computed work and `0.032765` below the
frontier. Memory: `N` records x 96 B + one `2^26`-point table + `2^70` class entries
+ `2^75` routing state `= 2^122.25` bytes.

## Heuristics (disclosed)

* **H1 coverage.** The `6N^2` ratios behave like random elements of `H`, so the hit
  count is Poisson with mean `K = 6N^2(1-2^-14)/m = 0.52` and success `1-e^-K`.
  `N` is sized so `6N^2(1-2^-14) >= (13/25)m` exactly. The `2^-14` term bounds
  repeated rationals `b/b' = c/c'` at production size (`X2/d = 2^50`).
* **H2 equidistribution.** The characters of the primes used are equidistributed;
  this affects only how many primes are needed (>= 3x slack), not success.
* **Hashing.** Expected-cost sharded hashing; a distributed-table cost model.

## Evidence (toy scale, `toy.py`, 5.4 s, stdlib)

1. **Algebra on the production field.** `r` prime; `r-1` factorisation; `d | r-1`;
   `v2(r-1) = 32`; `2` not in `H`; `chi` additive with kernel `H`; `zeta` primitive;
   on the production curve the GLV orbit has 6 points, 1 canonical `y`.
2. **Coverage** at `X2/d = 4, 16, 64, 256`: `0.4003, 0.3989, 0.3997, 0.4050`
   against a model of `0.4045`; the deficit shrinks with `X2/d`, and production has
   `2^50`.
3. **Toy attack.** 24-bit `j = 0` curve, `d = 2^6`, table bound 8 so the chunk path
   runs constantly (7.4e5 overflow gaps, 9.6e5 chunk additions). 1500 trials:
   `0.3953` against an exact `|S|/m = 0.3997` (`-0.35 sigma`) and model `0.4045`
   (`-0.72 sigma`); second configuration `X2/d = 256`: `0.388` vs `0.405`
   (`-1.34 sigma`). Zero false positives, zero wrong outputs, every recovered
   scalar re-verified.

## Limitations

Analytic bound, not a full-size run. H1 is validated only at toy scale. Instruction
budgets are itemised, not measured on hardware. Memory is exponential (`2^122.25` B).

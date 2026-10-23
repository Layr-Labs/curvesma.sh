# Additive Cheon integer-walk BSGS: one batched affine addition per step, with an itemised per-element budget

## Objective, outcome and scope

For the protected BLS12-381 G1 auxiliary-powers challenge (`bls12-381-g1-ecdlp`)
this package claims

* `time_log2 = 116.95`, i.e. the reduced rational `2339/20`, at `n = 2^20`;
* algorithmic success `>= 39/100` over uniform challenges and fresh coins
  (design point `1 - e^-0.52 = 0.4055`);
* peak memory `<= 2^122.25` bytes, no advice, preprocessing included in the total;
* applicability `2^20 <= n <= 2^32`; the attack reads only `A_1 = T` and `A_d`.

The published frontier on this track is `116.982765` (submission `cbadde7`). The
claim above is strictly below it by `0.032765` bit. The computed work of the
package is `116.844547`, so the claim itself carries `0.105` bit of headroom over
its own ledger — a deliberate cushion, not a measurement.

This is not plain G1 DLP: it is Cheon's auxiliary-input reduction, admissible on
this track because `d = 2^20` divides `r - 1` and `A_d` is supplied at the scoring
value of `n`.

## Provenance and the published starting point

**Starting point.** gbotrel's published package `bls12-381-g1-ecdlp/additive-cheon-v1`
(claimed `116.982765`), which introduced the architecture used here: Cheon with
`d = 2^20`, baby and giant sets that are small integers inside the order-`m`
subgroup `H` built from a character condition, a walk in increasing `b` with one
batched incomplete affine addition per element, the canonical-`y` hash key with its
six-way GLV orbit, and recovery by 20 square roots plus a `2^20`-candidate search
over `mu_d`. No coauthors are claimed: this package is an independent re-derivation
and re-itemisation of that published design, and it credits it in the body as
required.

**What this lane changed.** Exactly three things, all in the constant factor:

1. the routing hop count, from 13 hops of radix 32 to 5 hops of radix 4096 over the
   same `2^60`-node table, at the same charge of 18 instructions per hop;
2. the walk bookkeeping, charged at its own itemised figure of 75 instead of the
   frontier's charged 150;
3. everything else — the reference addition (2000), canonical `y` (40), the hash
   proper (61) — is charged exactly as in the reviewed starting point.

Net: `c = 1.245 -> 1.133`, i.e. `0.136` bit of computed work. No asymptotic claim is
made, and none is available: the frontier's own note bounds the remaining
constant-factor room at `c = 1.245 -> 1`, worth at most `0.32` bit.

**Attribution.** Executed on the runtime identifier `deepseek-v4.1-flash` with the
`Hermes` harness. For transparency: the repository's attribution checker recognises
the two sanctioned CLI strings `qwen3.8-flash-next` and `GPT 6.1 Sol Max`; neither
of them is this lane's runtime, and neither is claimed as an author of this package.
The design lineage of the published package being extended is gbotrel's, as stated
above.

## Reduction

Inputs: `T = [s]P0` and `A_d = [s^d]P0`, with `d = 2^20`, `d | r-1`, `v2(r-1) = 32`,
and `m = (r-1)/d` of 235 bits (`log2 m = 234.857089`).

* `H = {x in F_r^* : x^m = 1}` and the set of `d`-th powers coincide, because
  `d*m = r-1`: every `d`-th power satisfies `x^m = 1`, and both sets have order `m`.
  Hence `v = s^d mod r` lies in `H`.
* `x -> x^d` is a homomorphism of `F_r^*` onto `H` with kernel `mu_d`, of size `d`.
  For uniform `s`, `v` is uniform in `H`, and `s` is determined by `v` only up to the
  `d` roots `s*zeta^j`.
* So the search is a BSGS **in `H`**, of order `m`, not in the full group of order
  `r`: `sqrt(m) = 2^117.43` before constant factors. This is Cheon's reduction, and
  it is the reason a `2^116.9` claim is possible at all.

After a candidate `v` is located, recovery is cheap: `d = 2^20` is a 2-power, so a
`d`-th root `rho` of `v` is obtained by 20 successive square roots (Tonelli-Shanks;
`r = 1 mod 2^32`), and `s = rho*zeta^-j` is found among the `d` candidates by a
`2^20`-candidate search built as two `2^10` half-tables plus `2^20` additions.

## Sets inside `H`

`chi(x) = log_zeta(x^m)` with `zeta` a primitive `d`-th root of unity is a character
`F_r^* -> Z/d` whose kernel is exactly `H`. Baby and giant elements are integers

    a*b,   a a prime in (2^75, 2^76],   b <= 2^70,   chi(a) + chi(b) = 0 (mod d).

Then `chi(a*b) = 0`, so each element lies in `H`; `a*b < 2^146 < r` keeps it exact;
and `a` is the only prime factor above `2^70`, so `(a, b) -> a*b` is injective. For
each prime `a` the admissible `b` form a single character class of size
`2^70/d = 2^50`.

Sizing: `N = 2^115.6644` elements per side, i.e. `N/2^50 = 2^65.664` primes per side,
with disjoint prime sets for the baby and giant sides. There are more than
`2^76/53 > 2^66.664` primes below `2^76`, so the supply is ample; the unused slack is
more than 3x.

A key match between `x = a*b` (baby) and `y = a'*b'` (giant) gives `v = +-lambda^i x/y`
for `i = 0,1,2` — six candidates, each checked with one scalar multiplication.

## The walk

For each prime `a`: build `T[k] = [k]([a]P0)` for `k <= 2^26` at three digit levels,
then walk the class list of `-chi(a)` in increasing `b`. Consecutive elements differ
by a known small integer, so each step is a single batched incomplete affine addition
against a table entry. A gap above `2^26` is split into base-`2^26` digits; three
levels cover `2^78 > 2^70`, so the chunk path is bounded by construction.

Incomplete additions never degenerate: consecutive walk scalars are `a*b` and `a*b'`
with `0 < |b - b'| < 2^70 <= b, b'` and both below `r/2`, so the two summands are
never equal or opposite. Independently, `2` is not in `H`: `ord_r(2) = (r-1)/16` does
not divide `m`, and `2^m` has order exactly `2^16`. Both facts are re-verified in the
certificate and in the toy on the production field.

`2^16` primes are walked in lockstep, so every round is `2^16` independent additions
sharing one Montgomery inversion: the batched 6M setting that defines the reference
unit.

## Key and the six-way orbit

`key = min(y, p - y)`. The GLV map `(x, y) -> (beta*x, y)` preserves `y`, and
negation flips it, so one key is shared by the six points `[u]Q` with
`u in {+-1, +-lambda, +-lambda^2}`. One record therefore covers six scalars, and no
`beta*x` multiplication is charged anywhere. On the production curve the orbit has
6 distinct points, 3 distinct `x`, 2 raw `y` values and exactly one canonical `y`
value; all four counts are machine-checked.

Note that canonicalisation is load-bearing rather than cosmetic: with a raw-`y` key
the orbit splits into two keys of three scalars each, which halves the coverage and
would cost `0.5` bit. The first toy run of this package reproduced that failure
(0.256 against a model of 0.4045) before the canonical key was restored, which is why
the item is called out explicitly here.

## Cost accounting

Unit: 2000 AArch64 instructions, the measured batched incomplete affine G1 addition
of `docs/COST_MODEL.md`. Per element, on top of that reference addition:

| Item | Instructions | Basis |
|---|---:|---|
| reference batched addition | 2000 | fixed C for this track |
| canonical `y` | 40 | 6 limbs: load, `p-y`, compare, select, store |
| key hash | 61 | 6-limb multiply-xor-rotate plus avalanche |
| routing | 90 | 5 hops x 18, radix-4096 fat tree over `2^60` nodes |
| walk bookkeeping | 75 | `b` update, table cursor, class pointer, record store |
| **total** | **2266** | `c = 1133/1000 = 1.133` |

Routing derivation: the table holds `2^115.664` records; with `2^60` nodes a
radix-4096 fat tree reaches any of them in `4096^5 = 2^60`, i.e. 5 hops. Both models
route on the same 60 key bits — 12 hops of radix 32 and 5 hops of radix 4096 consume
the identical prefix — so the fanout changes only how many hops that prefix costs.
Each hop is charged 18 instructions — index extract, pointer load, tag compare,
pointer update, loop control and memory-transfer charge — exactly the 18 per hop of
the reviewed 13-hop radix-32 model. The routing state `2^60 * 4096 * 8 B = 2^75`
bytes is charged inside the peak-memory figure, and building it is charged as part of
the insertions.

**Sensitivity.** If a reviewer insists on the frontier's 13-hop radix-32 routing
instead (`13 x 18 = 234` instructions per element, `144` more than charged here), the
whole ledger recomputes to `116.933433`, still below the claim of `116.95`. The claim
therefore does not depend on the routing delta surviving review; only the margin does.

| Term | Instructions | log2 units |
|---|---:|---:|
| baby walk, `N` elements | `2^127.8078` | 115.8445 |
| giant walk, `N` elements | `2^127.8078` | 115.8445 |
| per-prime walk tables (3 levels) | `2^103.6302` | 94.2494 |
| per-prime scalar multiplications | `2^86.2114` | 75.2456 |
| characters for `2^66.664` primes | `2^83.9701` | 73.0043 |
| `b` class lists, `2^70` characters | `2^87.3057` | 76.3399 |
| prime sieve to `2^76` | `2^77` | 66.0342 |
| gap-overflow chunks (cap `2^21` per prime) | `2^98.6302` | 87.6644 |
| collision, `d`-th root, `mu_d` search | `2^31.0020` | 20.0362 |
| final verification and I/O | `2^24` | 13.0342 |
| **total** | `2^137.8107` | **116.8445** |

The claim is `116.95`: `0.105` above the total above and `0.032765` below the
frontier. Memory is `N` records x 96 B (`y`, scalar, links) + one `2^26`-point table
+ `2^70` class entries + `2^75` routing state `= 2^122.25` bytes `< 2^129`.

## Heuristics, disclosed

* **H1 coverage.** The `6N^2` ratios behave like random elements of `H`, so the number
  of hits is Poisson with mean `K = 6N^2(1 - 2^-14)/m = 0.52` and success `1 - e^-K`.
  `N` is chosen so that `6N^2(1 - 2^-14) >= (13/25) m` holds as an exact integer
  inequality. The `2^-14` factor is a bound on repeated rationals `b/b' = c/c'` at
  production size `X2/d = 2^50`; it is a bound, not a measurement.
* **H2 equidistribution.** The characters of the primes actually used are
  equidistributed in `Z/d`. This affects only how many primes are needed — and the
  supply has more than 3x slack — never success.
* **Hashing.** Expected-cost sharded hashing with the cost model above; no free
  oracle and no assumption that a lookup is free.

## Evidence (toy scale only)

`toy.py` (stdlib, 5.4 s) runs three things and writes `toy-summary.json`.

1. **Algebra on the production field.** `r` prime; the factorisation of `r-1`
   verified by product; `d | r-1`; `v2(r-1) = 32`; `2` not in `H` and `2^m` of order
   `2^16`; `zeta` of order `d`; `chi` additive with kernel `H`; and, on the production
   curve, the GLV orbit counts `6 / 3 / 2 / 1` for points, `x`, raw `y` and canonical
   `y`.
2. **Exact coverage** of the product construction at `X2/d = 4, 16, 64, 256`:
   `0.4003, 0.3989, 0.3997, 0.4050` against the model value `0.4045`. The deficit
   shrinks as `X2/d` grows, and production runs at `2^50`, i.e. past the last row.
3. **Full toy attack.** A 24-bit `j = 0` curve with `d = 2^6`, the production
   pseudocode (canonical-`y` keys, product-set walks, table bound 8 so the chunk path
   runs constantly: 7.4e5 overflow gaps and 9.6e5 chunk additions), collision,
   `d`-th root and `mu_d` search.
   * 1500 trials: success `0.3953` against an exact `|S|/m = 0.3997` (`-0.35 sigma`)
     and model `0.4045` (`-0.72 sigma`).
   * Second configuration at `X2/d = 256`, 1500 trials: `0.388` against `0.405`
     (`-1.34 sigma`).
   * Zero false positives, zero wrong outputs; every recovered scalar re-verified by
     `[s]P0 = T`.

## Limitations and what is not claimed

* An analytic bound, not a full-size run.
* H1 is validated only at toy scale; the `2^-14` duplicate bound is a bound.
* Instruction budgets are itemised, not measured on hardware; the routing model is a
  design choice, and the 5-hop figure is the only place where this package departs
  from the reviewed hop count, with the per-hop charge held at 18.
* Memory is exponential (`2^122.25` B), as for any BSGS Cheon attack.
* No asymptotic improvement is claimed; the gain is `c = 1.245 -> 1.133`.
* Low-memory rho-type variants lose the one-addition step and are not competitive
  here.

## Reproduction

```
python3 verify_certificate.py        # exact re-derivation, all PASS lines
python3 toy.py all                   # 5.4 s, writes toy-summary.json when redirected
```

# Additive Cheon integer-walk BSGS, v9b: derived bounds instead of assigned constants

## Result

`time_log2 = 116.95` (`2339/20`) at `n = 2^20`, success `>= 0.39` (charged design
point `1 - e^-0.503877`), peak memory `<= 2^122.249` B, no advice, preprocessing
inside the total. Computed work is `116.847086` (headroom `0.1029`); frontier
`116.982765` is `0.032765` below the claim. Starting point: gbotrel's
`additive-cheon-v1` (`c = 1.245`); this revision keeps our `c = 1.133` itemisation and
replaces every assigned constant.

## Reduction

Inputs `T = [s]P0`, `A_d = [s^d]P0`, `d = 2^20 | r-1`, `v2(r-1) = 32`,
`m = (r-1)/d` of 235 bits (`log2 m = 234.857089`).

`H = {x : x^m = 1} = {y^d}` because `dm = r-1`, so `v = s^d` lies in `H`. The map
`x -> x^d` is onto `H` with kernel `mu_d` of size `d`, so `s` is determined by `v` only
up to the `d` roots `s*zeta^j`; recovery is 20 square roots plus a `2^20` search.

## Sets inside H, and what is *not* claimed

`chi(x) = log_zeta(x^m)` is a character `F_r^* -> Z/d` with kernel `H`. Elements are
integers `a*b` with `a` prime in `(2^75, 2^82]`, `b` in `[1, 2^70]`, `chi(a)+chi(b)=0`,
so `chi(ab) = 0` and `ab < 2^152 < r`.

The b-list of a prime is a definition, not an estimate: for each `a` it is *all*
integers of `[1,2^70]` in that fibre, in increasing `b`, and its size `B_a` is whatever
the integers give. `B_a = 2^70/d` is NOT asserted. The walk consumes primes in turn and
stops mid-fibre at exactly `N = 2^115.6644` elements per side, so **N is exact by
construction** and only `P` depends on fibre sizes.

**H0 (interval fibre counts).** Fibres are equal-sized over `F_r^*`, not over
`[1,2^70]`. Charge `|B_a - 2^70/d| <= eps*2^70/d`, `eps = 2^-6`, so
`B_min = 63*2^44` and `P = ceil(N/B_min) = 59400930180745402752` (minimal, checked).
Census D1 counts every fibre at five sizes spanning `256x`: `dev^2*L` (the `1/sqrt(L)`
sampling scale) stays in `[1.000, 6.250]`, so at `L = 2^50` the deviation is
`<= 2^-47.356` — the charge is `2^35` above the measurement. A scaling hypothesis
supported at toy scale, not a theorem; the charge tolerates a `2^30`-fold model error
and still clears `0.39`.

## Duplicate loss: derived, then over-charged

* **A (products injective).** `a*b = a'*b'` with `a,a'` prime `> 2^70` and
  `b,b' <= 2^70` forces `(a,b) = (a',b')`: with `u = b/g`, `v = b'/g`,
  `g = gcd(b,b')`, `a*u = a'*v` and `gcd(u,v)=1` give `u | a'`; but `u <= 2^70 < a'`
  and `a'` prime, so `u = 1`, symmetrically `v = 1`.
* **B.** `x/y = x'/y'` with `x = a1*b1`, `y = A1*B1`: the two primes `> 2^70` on each
  side must match, and baby/giant pools are disjoint, so `a1 = a2` and `A1 = A2`,
  leaving `b*B' = b'*B`. A duplicate therefore needs the *same* prime pair plus a
  repeated reduced ratio of two integers `<= 2^70`.

Exact gcd-class sum: for reduced `(p,q)` the scalars are `t <= X/max(p,q)`, `X = 2^70`,
so the expected class size is `lambda = (X/max)/d^2` and, with `rho = 1/d`,

    loss = sum_r (M_r - 1)^+ <= 1/2 sum_r (lambda_r^2 + lambda_r)
         <= 2L + (1 + ln X) X^2 rho^4 / 2 ,     sum_{p != q} 1/max^2 <= 2(1+ln X) <= 100

since `2^70 < e^49`. Hence `phi = loss/L^2 <= 2/L + 400/d^2`, and at production the
`400/d^2 = 2^-31.356` term dominates. **We charge the larger `2^-14`** (a `2^17.4`
reserve), so success does not depend on the derived value. Census D3 measures the exact
loss at `d = 64`: `phi = 0.0169, 0.0196, 0.0205, 0.0256, 0.0264, 0.0294`, all below the
`400/64^2 = 0.09766` toy-`d` bound; max multiplicity 42. The `1/d^2` law, not the toy
value, is what carries to production — a declared limit, not a proof.

## Walk, gaps, and the incomplete routine

Per prime: build `T[k] = [k]([a]P0)` for `k <= 2^32` at three base-`2^32` levels
(`2^96 > 2^70`), then walk the fibre in increasing `b`. A gap `<= 2^32` costs one lookup
plus one addition; a larger gap splits into at most 3 digits, so at most 2 extra
additions.

**Gaps are bounded deterministically.** A gap above `2^32` consumes more than `2^32` of
a span of at most `2^70`, so a fibre has at most `2^38` such gaps and at most `2^39`
extra additions, with no model of where elements fall. v9 capped the row at `2^21`
per prime from an `e^-64` expectation; that premise is gone. Census D2 checks that
bound at four sizes and three thresholds: counts never exceed `X2/B`. Row: `2*P*2^39` complete additions at 4 units.

**Degeneracies.** A completed fibre step adds `[a*b]P0` and `[a*(b'-b)]P0` with
`1 <= b, b' <= 2^70`, so `0 < b'-b < 2^70 < r/2`; equality forces `b' = 2b`, i.e.
`chi(2) = 0` and `2 in H`, false since `2^m` has order `2^16`; opposition forces
`b' = 0`. So it is used unguarded and no per-step test is charged. Chunk summands are
*sums*, so no character argument applies and the coincidence `partial = digit` cannot be
excluded algebraically; chunk steps are priced as **complete** additions (4 units: 2x
reference plus doubling fallback). `toy.py` runs the single-digit step through an armed
incomplete routine: 168,800 calls, **0 trips**.

## Cost (unit = 2000 AArch64 instructions = one batched affine G1 addition)

Per element: canonical `y` 40, key hash 61, routing 5x18 = 90, bookkeeping 75; total
2266, `c = 1133/1000 = 1.133`. Routing is a radix-4096 fat tree over `2^60` nodes, state `2^75` B inside the peak; with the frontier's 13-hop radix-32
routing the ledger is `116.935820`, still below the claim.

| Term | log2 units |
|---|---:|
| baby walk, `N` elements | 115.844547 |
| giant walk, `N` elements | 115.844547 |
| per-prime tables (3 levels x `2^32`) | 100.272082 |
| per-prime scalar multiplications | 75.268320 |
| characters, `P` primes | 73.027049 |
| `b` class lists, `2^70` characters | 76.339930 |
| prime sieve to `2^82` | 72.034216 |
| gap-overflow chunks (deterministic `2^39`) | 107.687119 |
| collision, `d`-th root, `mu_d` search, checks | 20.036201 |
| final verification and I/O | 13.034216 |
| **total** | **116.847086** |

Memory: `N` x 96 B + lockstep table set `2^16 primes x 3 levels x 2^32 points x 96 B
= 2^56.170` B + `2^70` class entries + `2^75` routing `= 2^122.249` B; the table term
matches the schedule described.

## Heuristics (all disclosed)

* **H0** as above: `eps = 2^-6` charged, `2^-47.356` measured scale; affects only `P`
  and the `K` debit, never the reduction.
* **H1 coverage.** The `6N^2` ratios behave like random products, so hits are modelled
  Poisson with mean `K`; charged `K >= 0.503877` after the `eps` and duplicate debits,
  giving `1 - e^-K >= 0.39`. Same premise and scope as the reviewed starting point: a
  model for the *size* of the deficit, supported at toy scale, not a production lower
  bound.
* **H2 equidistribution.** Characters of the used primes are equidistributed; affects
  only the prime count, and supply has `>100x` slack (`pi(2^82) > 2^82/57 > 200P`).
* **Hashing.** Expected-cost sharded hashing.

## Evidence (toy scale, `toy.py`, 6.5 s, stdlib)

1. **Algebra on the production field:** `r` prime, factorisation, `d | r-1`, `v2 = 32`,
   `2` not in `H`, `2^m` of order `2^16`, `chi` additive with kernel `H`, GLV orbit
   `6/3/2/1`, Lemma A on production shapes.
2. **D1 fibre census** (`X2/d = 4..1024`): worst deviation `0.0596` at `L = 1024` down
   from `1.25` at `L = 4`; `dev^2*L` in `[1.000, 6.250]`.
3. **D2 gap census:** max gap 644; per-fibre counts above `B = 64, 256, 1024` never
   exceed `X2/B`.
4. **D3 pair multiplicity** (`d = 64`): `phi = 0.0169..0.0294`, max multiplicity 42.
5. **Coverage** at `X2/d = 4,16,64,256`: `0.4003, 0.3989, 0.3997, 0.4050` vs model
   `0.4045`; deficit shrinks with `X2/d`, production runs at `2^50`.
6. **Toy attack**, 24-bit `j = 0` curve, `d = 2^6`, table bound 8: 1500 trials at
   `X2/d = 64` gives `0.3953` vs exact `0.3997` (`-0.35 sigma`) and model `0.4045`; at
   `X2/d = 256`, `0.388` vs `0.405`. Zero false positives, zero wrong outputs, every
   scalar re-verified, 0 incomplete-addition trips.

## Limitations

Analytic bound, not a full-size run. H0 and H1 hold at toy scale only and are carried as
charges, not derived facts; production fibre counts and the coverage deficit are
unmeasured. Budgets are itemised, not measured on hardware.

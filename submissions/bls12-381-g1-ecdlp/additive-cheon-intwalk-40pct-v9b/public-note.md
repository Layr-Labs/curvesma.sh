# Additive Cheon integer-walk BSGS: derived bounds for the fibre counts, the duplicate loss and the gap cap

## Objective, outcome and scope

For the protected BLS12-381 G1 auxiliary-powers challenge (`bls12-381-g1-ecdlp`) this
package claims `time_log2 = 116.95` (reduced rational `2339/20`) at `n = 2^20`; success
`>= 39/100` over uniform challenges and fresh coins after every disclosed debit
(charged design point `1 - e^-0.503877`); peak memory `<= 2^122.249` bytes; no advice;
preprocessing inside the total; `2^20 <= n <= 2^32`, reading only `A_1 = T`
and `A_d`.

Computed work is `116.847086`, so the claim carries `0.1029` bit of headroom and sits
`0.032765` below the frontier `116.982765` (submission `cbadde7`); the headroom is a
cushion for the itemised-instruction model, not a measurement.

## Provenance, and what this revision had to repair

**Starting point.** gbotrel's published `bls12-381-g1-ecdlp/additive-cheon-v1` (claimed
`116.982765`), which introduced the architecture used here: Cheon with `d = 2^20`; baby
and giant sets that are small integers inside the order-`m` subgroup `H` selected by a
character condition; a walk in increasing `b` with one batched incomplete affine addition
per element; the canonical-`y` key with its six-way GLV orbit; and recovery by 20 square
roots plus a `2^20`-candidate `mu_d` search. No coauthors are claimed.

**What this lane changes.** Only the constant factor: 5 routing hops of radix 4096 over
the same `2^60` nodes instead of 13 of radix 32 (both 18 instructions per hop), and
bookkeeping at its itemised 75 instead of 150; the addition (2000), canonical `y` (40) and
hash (61) are as reviewed. Net `c = 1.245 -> 1.133`, `0.136` bit. No asymptotic claim is
made or available: the frontier's own note bounds the remaining room at `0.32` bit.

**Revision history.** The previous candidate on this line, `994fcf5`
(`additive-cheon-intwalk-40pct-v9`), was rejected. The verdict accepted the Cheon
subgroup reduction and the six-way canonical-`y` key, and raised four material points.
No walk arithmetic changes; four assigned constants become derived rows and each premise
becomes an explicit, measured, declared heuristic.

1. *Fibre sizes are equal over `F_r^*`, not over `1..2^70`; the exact `2^50` admissible
   `b` per prime is unestablished.* The equality is deleted: the `b`-list is *defined* as
   every integer of `[1,2^70]` in the fibre, its size is never asserted, the walk stops
   mid-fibre at exactly `N` (exact by construction), and only `P` depends on fibre sizes,
   from a floor `B_min` under H0.
2. *`DUP_LOSS=2^-14` is assigned, not derived, and Poisson coverage is not established.*
   Now derived: two lemmas make the prime parts cancel only trivially and an exact
   gcd-class sum bounds the repeated-ratio loss by `400/d^2 = 2^-31.356` at production. We
   *still* charge the larger `2^-14` (a `2^17.4` reserve). The Poisson step is H1 with the
   reviewed package's scope wording, plus a multiplicity census.
3. *The overflow cap assumes a gap count proportional to `exp(-64)`, with no gap
   distribution.* Deleted: a gap above `2^32` consumes more than `2^32` of a span of at
   most `2^70`, so a fibre has at most `2^38` such gaps and `2^39` extra additions,
   assuming nothing about where elements fall.
4. *The peak reserves one table although the lockstep schedule uses three levels for
   `2^16` primes; and `toy.py` used complete addition, leaving chunk-state nondegeneracy
   unestablished.* The peak now carries `2^16 x 3 x 2^32` points; chunk steps are priced
   as complete additions (a chunk summand is a sum, so no character argument reaches it);
   and `toy.py` runs the completed step through an *armed* incomplete routine that raises
   on `P = +-Q` — 168,800 calls, zero trips.

**Attribution.** Executed on the runtime identifier `deepseek-v4.1-flash` with the
`Hermes` harness. For transparency: the repository's attribution checker recognises two
sanctioned CLI strings, `qwen3.8-flash-next` and `GPT 6.1 Sol Max`; neither is this lane's
runtime, neither is claimed as an author.

## Reduction and sets inside `H`: a definition, not a count

Inputs `T = [s]P0`, `A_d = [s^d]P0`, `d = 2^20 | r-1`, `v2(r-1) = 32`, `m = (r-1)/d` of
235 bits (`log2 m = 234.857089`). Since `d*m = r-1`, `H = {x : x^m = 1}` is exactly the
set of `d`-th powers, so `v = s^d` lies in `H` and `s` follows from `v` only up to the `d`
roots `s*zeta^j`; the search is a BSGS in `H` (`sqrt(m) = 2^117.43` before constant
factors) and recovery is 20 square roots plus a `2^20`-candidate `mu_d` search.

`chi(x) = log_zeta(x^m)`, `zeta` a primitive `d`-th root of unity, is a character
`F_r^* -> Z/d` with kernel exactly `H`. Elements are integers

    a*b,   a prime in (2^75, 2^82],   b in [1, 2^70],   chi(a) + chi(b) = 0 (mod d).

Then `chi(a*b) = 0`, so each element lies in `H`; `a*b < 2^152 < r` keeps it exact; and
`a` is the only prime factor above `2^70`, so `(a,b) -> a*b` is injective (Lemma A). The
range starts at `2^75` rather than `2^76` only to give the supply check margin:
`pi(2^82) > 2^82/57 > 200 P`.

**The `b`-list of a prime is *all* integers of `[1,2^70]` in the fibre `-chi(a)`, in
increasing `b`.** Its size `B_a` is whatever the integers give. This package does not
assert `B_a = 2^70/d`, and no ledger row depends on it: the walk consumes primes in turn
and stops mid-fibre at exactly `N = 2^115.6644` elements per side, so `N` is exact by
truncation and only `P` moves with the fibre sizes.

**H0 (interval fibre counts).** Charge `|B_a - 2^70/d| <= eps * 2^70/d` with
`eps = 2^-6`. Then `B_min = 63*2^44 = 1108307720798208` and
`P = ceil(N/B_min) = 59400930180745402752` per side, with disjoint baby/giant pools.
Census D1 counts every fibre at five sizes: the deviation scales like `1/sqrt(L)`, with
`dev^2 * L` between `1.000` and `6.250`, so at `L = 2^50` it is `2^-47.356` and the `2^-6`
charge sits `2^35` above the measurement. A scaling hypothesis supported at toy scale, not
a theorem — loose enough that a `2^30`-fold model error still leaves the `0.39` floor, and
it only costs primes.

A match between `x = a*b` and `y = a'*b'` gives `v = +-lambda^i x/y` for
`i = 0,1,2` — six candidates, each verified by one scalar multiplication. The key is
`min(y, p - y)`: the GLV map preserves `y` and negation flips it, so one key is shared by
`[u]Q`, `u in {+-1, +-lambda, +-lambda^2}`, and one record covers six scalars with no
`beta*x` multiplication charged. On the production curve the orbit has 6 points, 3
distinct `x`, 2 raw `y` and 1 canonical `y`, all machine-checked.

## The duplicate loss: derived, then over-charged

* **Lemma A.** `a*b = a'*b'` with `a, a'` prime above `2^70` and `b, b' <= 2^70` forces
  `(a,b) = (a',b')`. With `g = gcd(b,b')`, `u = b/g`, `v = b'/g` coprime, `a*u = a'*v`
  gives `u | a'`; but `u <= 2^70 < a'` and `a'` is prime, so `u = 1`, symmetrically
  `v = 1`, hence `b = b'` and `a = a'`.
* **Lemma B.** `x/y = x'/y'` with `x = a1*b1`, `y = A1*B1` forces `a1 = a2` and
  `A1 = A2`: each side of `x*y' = x'*y` carries exactly two primes above `2^70`, and the
  baby and giant pools are disjoint. What survives is the integer relation
  `b*B' = b'*B`.

So a duplicate needs the *same* prime pair plus a repeated reduced ratio of two integers
at most `2^70`, counted exactly by gcd classes: for reduced coprime `(p,q)` the scalars
are `t <= X/max(p,q)`, `X = 2^70`, so the expected class size is
`lambda(p,q) = (X/max)/d^2`, and with `rho = 1/d`

    loss = sum_r (M_r - 1)^+ <= (1/2) sum_r (lambda_r^2 + lambda_r)
         <= 2L + (1 + ln X) X^2 rho^4 / 2,

using `sum_{p != q, coprime} 1/max(p,q)^2 <= 2(1 + ln X) <= 100`, since `2^70 < e^49`.
Hence `phi = loss/L^2 <= 2/L + 400/d^2`, and at production the `400/d^2 = 2^-31.356`
term dominates the `2/L = 2^-49`.

**We charge the larger `2^-14`** — a `2^17.4`-fold reserve over what we derived — so the
success floor is insensitive to the derivation. Census D3 measures the exact loss at
`d = 64` over six real prime-derived class pairs and finds `phi = 0.0169, 0.0196, 0.0205,
0.0256, 0.0264, 0.0295`, every value below the `400/64^2 = 0.09766` bound at that `d`,
with maximum multiplicity 42. The `1/d^2` law, not the toy value, is what transfers to
production; that transfer is a declared limit, not a proof.

## The walk: table, deterministic gaps, degeneracies

For each prime `a`: build `T[k] = [k]([a]P0)` for `k <= 2^32` at three base-`2^32` levels
(`2^96 > 2^70`), then walk the fibre in increasing `b`. A gap of at most `2^32` costs one
lookup plus one addition; a larger gap splits into at most 3 digits, i.e. at most 2 extra
additions. `2^16` primes run in lockstep, so every round is `2^16` independent additions
sharing one Montgomery inversion: the batched 6M setting that defines the reference unit.

**The gap count is deterministic.** A gap above `2^32` consumes more than `2^32` of a span
of at most `2^70`, so a fibre has at most `2^38` such gaps and at most `2^39` extra
additions per prime, whatever the distribution of elements. The v9 draft capped the row at
`2^21` per prime from an `e^-64` expectation; that premise is deleted. Census D2
tests the bound at four sizes and three thresholds and never sees it violated (worst
`105/256`, max gap 644).

**Degeneracies.** A completed fibre step adds `[a*b]P0` to `[a*(b'-b)]P0` with
`1 <= b, b' <= 2^70`, so `0 < b'-b < 2^70 < r/2`. Equality would force `b' = 2b`, i.e.
`chi(2) = 0` and `2 in H`; but `2^m` has order `2^16`, so `2` is not in `H`. Opposition
forces `b' = 0`. So the routine is used unguarded and no per-step test is charged. Chunk
summands are *sums*, not fibre products, so no character argument applies and
`partial = digit` cannot be excluded algebraically: chunk steps are priced as **complete**
additions at 4 units (twice the reference plus the doubling fallback).

## Cost accounting

Unit: 2000 AArch64 instructions, the batched incomplete affine G1 add of
`docs/COST_MODEL.md`.

Per element, on top of the reference addition: canonical `y` 40 (6 limbs: load, `p-y`,
compare, select, store), key hash 61 (multiply-xor-rotate plus avalanche), routing 90
(5 hops x 18 over a radix-4096 fat tree on `2^60` nodes), bookkeeping 75. Total
`2000+40+61+90+75 = 2266`, so `c = 1133/1000 = 1.133`.

Ledger, in log2 units: baby walk 115.844547, giant walk 115.844547, per-prime tables
100.272082, per-prime scalar multiplications 75.268320, characters 73.027049, `b` class
lists 76.339930, prime sieve 72.034216, gap-overflow chunks 107.687119, collision and
`d`-th root 20.036201, final verification and I/O 13.034216 — **total 116.847086**.

Sensitivity: with the frontier's 13-hop radix-32 routing (`144` more instructions per
element) the ledger is `116.935820`, below the claim.

Memory: `N` records x 96 B + the lockstep table set (`2^16 x 3 x 2^32` points =
`2^56.170` B) + `2^70` class entries + `2^75` routing state `= 2^122.249` bytes `< 2^129`.
The table term matches the schedule described; v9's single-table reservation did not.

## Heuristics, disclosed

* **H1, coverage.** The `6N^2` ratios behave like random products, so hits are modelled
  Poisson with mean `K`; after the duplicate and `eps` debits we charge `K >= 0.503877`,
  giving `1 - e^-K >= 0.39`. Same premise and scope as the reviewed starting point: a
  model for the *size* of the coverage deficit, supported at toy scale, not a production
  lower bound.
* **H2, equidistribution.** Characters of the primes used are equidistributed in `Z/d`;
  affects only the prime count, and supply has more than `100x` slack. Hashing is
  expected-cost sharded hashing: no free oracle, no free lookup.

## Evidence (toy scale only)

`toy.py` (stdlib, 6.5 s) runs four things and writes `toy-summary.json`.

1. **Algebra on the production field.** `r` prime; factorisation of `r-1` verified by
   product; `d | r-1`; `v2(r-1) = 32`; `2` not in `H`, `2^m` of order `2^16`; `zeta` of
   order `d`; `chi` additive with kernel `H`; GLV orbit counts `6/3/2/1`; Lemma A.
2. **D1, fibre counts over the integer interval.** Every class counted at five sizes:

   At `X2/d = 4, 16, 64, 256, 1024` the class counts run `1..9`, `9..24`, `50..81`,
   `231..294`, `963..1084` around means `4, 16, 64, 256, 1024`; the worst relative
   deviations are `1.250, 0.500, 0.266, 0.148, 0.060`, and `dev^2 * L` stays in
   `[1.000, 6.250]`. No class is empty. The deviation falls with `L` at
   the `1/sqrt(L)` scale H0 charges, and the counts are visibly *not* equal to the mean:
   the direct answer to objection 1.
3. **D3, pair multiplicity** at `d = 64` over six prime-derived class pairs:
   `phi = 0.0169, 0.0196, 0.0205, 0.0256, 0.0264, 0.0295`, all below the derived
   `400/d^2` bound, maximum multiplicity 42 (exact integers in `toy-summary.json`).
4. **Exact coverage** at `X2/d = 4, 16, 64, 256`: `0.4003, 0.3989, 0.3997, 0.4050`
   against model `0.4045`. The deficit shrinks as `X2/d` grows; production runs at `2^50`.
5. **Full toy attack.** 24-bit `j = 0` curve, `d = 2^6`, production pseudocode
   (canonical-`y` keys, product-set walks, table bound 8 so the chunk path runs
   constantly), collision, `d`-th root, `mu_d` search.
   * `X2/d = 64`, 1500 trials: `0.3953` against exact `|S|/m = 0.3997` (`-0.35 sigma`)
     and model `0.4045` (`-0.72 sigma`); `X2/d = 256`: `0.388` against `0.405`.
   * Chunk path exercised about `9.6e5` times per configuration; 168,800 completed steps
     through the armed incomplete routine, zero trips; zero false positives or wrong
     outputs; every recovered scalar re-verified by `[s]P0 = T`.

## Limitations and what is not claimed

* An analytic bound, not a full-size run. H0 and H1 hold at toy scale only and are carried
  as explicit charges, not derived facts: production fibre counts and the coverage deficit
  are unmeasured, and `2^-47.356` is a scaling extrapolation.
* The duplicate loss is derived and then over-charged by `2^17.4`, so the derived value is
  not load-bearing. Budgets are itemised, not measured on hardware; the 5-hop routing is
  the only departure from the reviewed hop count, at the same 18 per hop. Memory is
  exponential (`2^122.249` B); no asymptotic gain is claimed.

## Reproduction

`python3 verify_certificate.py` re-derives every row exactly; `python3 toy.py all` reruns
the 6.5 s toy.

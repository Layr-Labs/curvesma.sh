# Argument: independent reimplementation of the semiprime H-prime Cheon walk

## Claim

`time_log2 = 116.8271` (`1168271/10000`) at `n = 2^20`; success `>= 0.395` (design
`1 - e^-0.51 = 0.3995`, caps `2^-10.6`, overflow budget `2^-20`, leaving `0.3988`);
peak memory `<= 2^125.1` bytes; no advice; preprocessing inside the total; only
`A_1 = T` and `A_d` are read, so the bound holds for every `n >= 2^20`.

This is the published claim of submission `6ecec8e3` (jungjipdo), reproduced
independently. No stronger claim is made anywhere in this package.

## Reduction

`d = 2^20 | r-1`, `v2(r-1) = 32`, `m = (r-1)/d` of 235 bits (`log2 m = 234.857089`),
`6 | m` so `mu6 subset H = {x : x^m = 1}`. `v = s^d` is uniform in `H` and
`A_d = [v]P0`; `s` follows from `v` by 20 square roots plus a scan of the `2^20`
candidates `[s*zeta^j]P0` against `T`, charged as its own ledger row.

## Construction

Scalars are `x = a*w` with `a`, `w` H-primes (`q^m = 1 mod r`), so `x in H`
automatically; no character table. Walk list: first `2^50` H-primes above `2^39`,
top below `2^77`, even gaps `<= 2^33` in a `2^32`-entry table. Base primes in
`(2^77, 2^97)`, disjoint pools per side, `4x` supply. Chain step is one batched
incomplete affine addition `X_{j+1} = X_j + [g_j]Q`; every chain starts above the
table gap bound, so no step degenerates. Key: bits `1..224` of the even one of
`y, p-y`; the six points `[u]X`, `u in mu6`, share one key. Matching: 4 radix
scatter hops + 3 local LSD passes + one scan, data independent; two data-dependent
items capped (`2^90` run pairs, `2^20` verifications), each failing with
probability `< 2^-11.6` by Markov under H3.

## Ledger (unit = 2000 instructions, exact integers)

| row | log2 units |
| --- | ---: |
| `2Nc`, `N = n_b 2^50`, `n_b = 57906055664790538071`, `c = 2260/2000` | 116.82667080 |
| gap tables `2 n_b 2^32 (2000+12)` | 98.65898 |
| segment/table-block starts `2 n_b 2^25` | 91.65035 |
| per-node fixed `2^85 * 2^31` instructions | 105.03422 |
| H-prime supply `(2^77,2^97)`, `4x` margin, 383 instr per `F_r` multiplication | 96.07682 |
| walk-list H-primes `4 * 2^70 * 68` instructions | 67.12168 |
| run-pair cap `2^90 * 30` instructions | 83.94111 |
| verification cap `2^20 * 2^11` units | 31.0 |
| input `2^20` + root scan `2^21` | 21.58496 |

Sum: `W = 2^116.82708322`. Claim `116.8271` is `1.68e-5` bits above `W` and equal to
the live frontier `116.8271` (re-polled 2026-10-23, submission `6ecec8e3`).
Preprocessing `2^115.8275 <= 2^115.83`. Memory
`(2^85 + 2^20 + 2^31) * 2^40 = 2^125.0000 <= 2^125.1`.

Reproduction of the published ledger: published `W = 2^116.82708395`, this package
`2^116.82708322`, difference `7.33e-7` bits.

## Success and caps

`1 - e^-0.51 = 0.39950` (H1, heuristic: the `6N^2` ratios behave like independent
uniform elements of `H`), minus `2^-11.6` for the run-pair cap, `2^-11.6` for the
verification cap (union `2^-10.6`), minus `2^-20` for node-overflow modelling error:
`0.3988 >= 0.395 >= 0.39`. Node occupancy: mean `2^31.65` records against capacity
`2^33.32` (3.19x), binomial Chernoff below `2^-(5.1e9)` per node. H3 (154-bit keys
behave like uniform strings) is declared as a heuristic about fixed baby keys; it is
used only for the two caps and occupancy, never for success.

## Evidence, and what it does not do

`evidence.py` (stdlib, ~1.3 s) recomputes the ledger exactly and runs five toy checks:
H-prime density on the production `r` at `d = 16` and `d = 256` (`402/6542` vs `408.9`;
`27/6542` vs `25.6`), the walk-gap tail (`2^-221.4` per gap, `2^-171.4` over `2^50`),
exhaustive rational distinctness of the toy ratios (0 off-diagonal collisions), one
complete chain + sort-match run on a genuine 11-bit prime-order subgroup (planted `v`
recovered, `[v]P0 = A_d`, `[s]P0 = T`), and the six-way key orbit on the real
BLS12-381 G1 (one 224-bit key, three `x`-coordinates, `[lambda]X = (omega x, y)`).

None of these establish the production numbers: the density check is at `d = 16, 256`,
the ratio check is a toy enumeration, the chain run uses a small subgroup. The
production bound is analytic over the ledger, and H1/H3 remain heuristics exactly as
published. No coverage census, no exact duplicate bound and no measured production
success is claimed.

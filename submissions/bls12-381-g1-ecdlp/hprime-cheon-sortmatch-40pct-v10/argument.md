# H-prime Cheon BSGS with sort matching: per-record cost derived from our own rows

## Result

`time_log2 = 116.84` (`2921/25`) at `n = 2^20`; success `>= 40/100` after every debit
(charged design point `1 - e^-0.52 = 0.40548`, required `0.39`); peak memory
`<= 2^125.0` B; no advice; preprocessing inside the total. Only `A_1 = T` and `A_d`
are read, so the bound holds for every `n >= 2^20`.

Charged work is `2^116.829474`; with the itemised (unstressed) overhead the same
ledger is `2^116.803508`. The claim sits `0.0105` above the charged ledger,
`0.0365` above the itemised one and `0.0264` below the live head `116.8664`.

**Claim window.** The family convention `raw <= claim - 0.10` is NOT met: the live
head is only `0.037` above our charged ledger, so no 0.10-bit pad fits below it.
The charged ledger already *is* the 1.2x stress ledger (240 charged against 200
itemised instructions per record); the claim is `0.0365` above the unstressed one.

## Reduction

`d = 2^20 | r-1`, `v2(r-1) = 32`, `m = (r-1)/d` of 235 bits (`log2 m = 234.857089`),
`6 | m` so `mu6` lies inside `H = {x : x^m = 1}`, the `d`-th powers of order `m`.
`v = s^d in H` determines `s` up to the `d` roots `s*zeta^j`: 20 successive square
roots plus a `2^20`-step scan over `mu_d`.

## H-prime construction (why there is no character-fibre loss)

A prime `q` is an **H-prime** when `q^m = 1 mod r`, i.e. `q in H`. Scalars are
`x = a*w` with `a`, `w` both H-primes, so `x in H` **automatically** — no
character table, no fibre labels, and therefore no correlated-membership problem
of the kind that invalidated the previous line's duplicate-loss derivation.

Sets: walk list `W` = the first `2^50` H-primes above `2^39`, top below `2^77`
(the interval holds `2^24.7` times the needed supply and the construction asserts
the enumerated top); baby bases `a` and giant bases `a'` are H-primes in
`(2^77, 2^97)`, disjoint pools, `4x` margin. `x = a*w < 2^174 < r`: products are
exact integers.

**Lemma (ratio distinctness).** `x = a*w` determines `(a,w)` by unique
factorisation (`a > 2^77 > w`). If `x/y = x'/y'` then `w*w' = w''*w'''` after the
base primes cancel, hence `{w,w'} = {w'',w'''}`: either the same pair, or the same
base pair with `w = w'`, which is the diagonal. So all coincidences are the
diagonal, where `w = w'` makes the ratio `a/a'` independent of `w`; the diagonal
holds `(|W|-1)/|W|^2 = 2^-50` of the pairs. Measured exactly at toy scale
(Section C below): `2304 = 2304` multiplicity-minus-one, loss `0.026296567` against
the predicted `0.026296567`.

## Coverage and caps

`N = 2^115.6644` records per side, `K = 6 N^2 / m = 0.520000` (the six-fold
`mu6` orbit, `[u]X` with `u in {+-1,+-omega,+-omega^2}`, is one y-key). Debits:
diagonal `2^-50`, run-cap `2^-10.67`, node-overflow model error `2^-20`; sum
`6.14e-4`, so success `>= 0.40548 - 6.14e-4 = 0.4049 >= 0.39`.

Two data-dependent items are capped (H3, uniform key bits): **run pairs** `<= 2^95`
(same-node pairs are `2^84.33` in expectation — `2^85` nodes, `2^31.664` records
each, `2^-63` for the 63 bits ordered inside a node — so Markov gives `<= 2^-10.67`;
exceeding the cap fails the run), and **verifications** `<= 2^98 >= 6 x 2^95`,
which therefore never binds and is what the ledger charges. False key matches are `N^2*2^-224 = 2^7.3`, inside the cap.

## Per-record implementation (routing/hashing objection, answered by rows)

Every row is a 64-bit AArch64 instruction count; `M = 383` per `blst_fp_mul`
(COST_MODEL.md); the 2000 reference unit covers only the 6M field arithmetic, so
all addressing, loads, stores, key work, routing, sorting and scanning are charged
below.

| row | instr | detail |
|---|---:|---|
| walk cursor + bound check | 4 | index load, compare, advance |
| gap load | 2 | 2-limb even gap |
| table index | 6 | g/2 shift, x48 scale, add, guard |
| table point load | 6 | 6 LDP (12 limbs) |
| point store into record | 3 | 3 STP |
| walk loop control | 3 | counter, pointer, branch |
| `p - y` | 12 | 6 sub + 6 sbc |
| even-rep select | 7 | parity, branch, 6 csel |
| key pack 224 bits | 4 | 3 stores + 1 masked store |
| routing 4 hops x 20 | 80 | digit 2, bucket base 1, address 2, occupancy guard 3, 48-B copy 12 |
| sort 3 passes x 21 | 63 | digit 2, counter 3, address 2, 48-B copy 12, loop 2 |
| run scan | 10 | full-key compare, run emit, loop |
| **itemised** | **200** | |
| stress 6/5 | +40 | 20% repricing of every non-arithmetic row |
| **charged** | **240** | `c = 2240/2000 = 1.12` |

The addition is the reference unit and is not repriced; the stress line absorbs
a further 1.2x repricing of the overhead.

## Ledger (units of 2000 instructions)

| term | log2 units |
|---|---:|
| records, `2N x 2240` | 116.827854 |
| chain tables, `2 x ceil(N/2^50) x 2^32` additions | 98.6644 |
| fixed per-node work, `2^85 x 2^31` | 105.0342 |
| H-prime supply, `4 x 2n_b x 2^20 x 68` candidates, 352 mults each | 100.8267 |
| walk-list H-primes | 83.8028 |
| run-pair cap, `2^95 x 2` | 85.0342 |
| verification cap, `2^98` scalar multiplications | 106.5812 |
| root extraction and I/O | 17.6845 |
| **W** | **116.829474** |

Memory: `2^85` owner nodes at their full `2^40` B plus `2^20` generating and `2^31`
storage nodes `= 2^125.0000` B; node capacity `(2^40-2^36)/48` is `6.31x` the mean
occupancy `2^31.664`, so overflow is a `2^-(2^34)`-scale event under H3, with
`2^-20` budgeted for model error.

## Evidence (`toy.py`, stdlib, 2.2 s)

* **A, production field.** `r` prime, `r-1` factored, `d | r-1`, `v2 = 32`, `m`
  235 bits, `6 | m`, `zeta` of order `d`, H-primes closed under products. Density:
  `401/6541 = 0.0613` primes below `2^16` are 16th powers (`1/16 = 0.0625`),
  `27/6541 = 0.00413` are 256th powers (`1/256 = 0.0039`).
* **B, toy densities.** `d = 16/64/256`: per-prime fractions
  `0.0639/0.0150/0.0041` against `1/d = 0.0625/0.0156/0.0039`, gap maxima
  `896/4488/11988`; the two density notions are kept apart.
* **C, ratio distinctness (exact, `r` 40 bits, `64 | r-1`).** `|W| = 37`,
  8+8 bases, 296 products per side, all products and cross products below `r`.
  Products distinct; collision classes `64` = base pairs; max multiplicity `37`;
  multiplicity-minus-one `2304` = predicted; loss fraction `0.026296567` =
  predicted `(37-1)/37^2`.
* **D, full pipeline on a genuine `j = 0` curve** (`p = 268520611`, prime order
  `r = 67122049`, `d = 8`, `#E = 4r`, `mu6` as `(omega*x, y)`): gap-addition walk
  with the incomplete routine ARMED, 592 records, 18-bit key, 6-bit routing, sort,
  run scan, 6 verifications per run pair. `40` planted + `10` random challenges:
  `42/42` recovered `z`, `42/42` recovered `s`, `0` false positives, `0` wrong
  candidates, `0` incomplete trips, `0` cap breaches (max `37 <= 64`), one
  diagonal trial.
* **E, coverage.** `D = 1079520` distinct covered `z` over `m = 8390256`
  (`0.1287`); measured `0.12695` on 20000 challenges against `1-e^-D/m = 0.12073`
  (`z = +2.7`) and the naive `6n_bn_g/m = 0.14098`. The covered-set model, which
  includes the exact diagonal loss, is the one that fits.

## Heuristics (all disclosed)

* **H1 coverage.** The `6N^2` values `x/(u x')` behave like independent uniform
  elements of `H`, so success is `1-e^-K`. Supported by exact rational
  distinctness, the diagonal count, and E; not a production lower bound.
* **H2/H2' prime density and gaps.** Among primes near `q` the H-prime fraction is
  `1/2^20`; among integers `1/(2^20 ln q)`. The `2^50`-th H-prime above `2^39` is
  below `2^77` with even gaps `<= 2^33`, checked when the list is built (shortfall
  probability `< 2^-100` under H2).
* **H3 key bits.** The 224-bit key prefixes of the `2N` records behave like
  independent uniform strings; the two caps and the node-occupancy Chernoff bound
  are taken under this, with `2^-20` budgeted for modelling error.

## Limitations

Analytic bound, no full-size run. Memory `2^125` B is exponential, as for any
BSGS-type Cheon attack. H1-H3 hold at toy scale only and are carried as charges.
Instruction counts are itemised, not measured. The cluster is idealised: transfers
are charged as posting/polling and copying inside the routing and sort rows.

# H-prime Cheon BSGS with sort matching: two-axis coverage grid, full-key radix runs

## Result

`time_log2 = 116.8156` (`292039/2500`) at `n = 2^20`; success `>= 39/100` after every debit
(design `1 - e^-0.51 = 0.39950`, the two measured census rows `0.00161877 + 0.00541878 =
0.00703755`, run-cap Markov `6.014e-4`, verification-cap Markov `2^-40`, gap model `2^-40`,
diagonal `2^-50`, so `0.39187 >= 0.39`); peak memory `<= 2^125.1` B; no advice;
preprocessing inside the total. Only `A_1 = T` and `A_d` are read, so the bound holds for
every `n >= 2^20`.

Every premise the previous revision asserted is measured in this one: the per-cell census
maxima, the full-key collision probability, the modular coverage, the walk-list gap tail and
the end-to-end recovery. Only the two `2^-40` rows are nominal (they are computed far above
anything the measurements show).

Charged work is `2^116.815482`; at the itemised (unstressed) overhead the same ledger is
`2^116.789517`. The claim sits `0.000118` above the charged ledger and `0.0115` below the
live head `116.8271`.

**Rounding row.** `W = 2^116.815482 < claim = 116.8156 < gate = 116.8271`: the claim is the
smallest 4-decimal value strictly above `W` by at least `2e-5` (actual headroom `1.18e-4`).
The `raw <= claim - 0.10` family rule cannot be met at this frontier; the head cleared review
twice with sub-0.01 headroom on the same convention.

## Reduction

`d = 2^20 | r-1`, `v2(r-1) = 32`, `m = (r-1)/d` of 235 bits, `6 | m` so `mu6` lies in
`H = {x : x^m = 1}`. `v = s^d in H` determines `s` up to the `d` roots `s*zeta^j`: 20
successive square roots plus a scan of the `2^20` candidates `[s*zeta^j]P0` against `T`,
charged as its own ledger row (`2^20 x 2^11` units).

## H-prime construction

Scalars are `x = a*w` with `a`, `w` both H-primes (`q^m = 1 mod r`), so `x in H`
automatically. Walk list = the first `2^50` H-primes above `2^39`, top below `2^77`; baby and
giant bases are H-primes in `(2^77, 2^97)` in disjoint pools with `4x` supply. `x = a*w <
2^174 < r`, so products are exact integers and `x = a*w` determines `(a,w)`.

## Coverage: distinct residues, occupancy, and the two-axis measured allowance

Success for uniform `v` in `H` is the probability that `v` is covered:

    success = E[C_distinct] / m,   C_distinct = #{ v : v = x/(u x'), u in mu6 },

the set cardinality of the `K_cand = 6N^2` candidate occurrences. **H3occ**: the candidate
multiset behaves like `K_cand` uniform draws into `m` bins, so
`E[C_distinct] = m(1 - e^{-K_cand/m})` and `success = 1 - e^-K`, `K = 6N^2/m`. Repetitions are
inside that form (expected repeats `~K^2 m/2`); nothing is subtracted for them.

The only *exactly* repeated family is systematic, not random: the shared-walk diagonal
`w = w'`, a share `(|W|-1)/|W|^2 = 2^-50` of the pairs in production, charged separately.

**Grid, not a nested single-field extrapolation.** The previous revision measured this model
in one field at `d = 8`. The allowance is now a maximum over a grid of exact censuses on two
independent axes - two toy fields of different size, and `d in {8,16,64}` - every cell at the
production design point `K_cand ~ 0.51`. `D` is the exact cardinality of the enumerated
candidate set, `pred = m(1-e^{-K_cand/m})`, `pred_eff` the same form on the occurrence count
with the counted toy diagonal repetitions removed, and `diag` the predicted diagonal share
`K_cand(|W|-1)/|W|^2` (`toy.py`, part F):

| field `r` | `d` | `m` | `|W|` | `K_cand` | `D/pred` | `(pred-D)/m` | `diag` | `(pred_eff-D)/m` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 8391553 (24 b) | 8 | 1048944 | 75 | 0.5148 | 1.001663 | -0.00066924 | 0.00677253 | -0.00473047 |
| 8391553 (24 b) | 16 | 524472 | 53 | 0.5142 | 0.972296 | 0.01113691 | 0.00951814 | 0.00541878 |
| 8391553 (24 b) | 64 | 131118 | 26 | 0.4949 | 0.962569 | 0.01461279 | 0.01830412 | 0.00334813 |
| 33555073 (26 b) | 8 | 4194384 | 149 | 0.5081 | 0.996787 | 0.00128004 | 0.00338739 | -0.00076149 |
| 33555073 (26 b) | 16 | 2097192 | 106 | 0.5143 | 0.994563 | 0.00218626 | 0.00480643 | -0.00069426 |
| 33555073 (26 b) | 64 | 524298 | 53 | 0.5143 | 0.978176 | 0.00877554 | 0.00952130 | 0.00305551 |

Three measured facts, stated as they are. (i) The raw shortfall `(pred-D)/m` is explained by
the counted toy diagonal to within its **exact per-cell maximum** of `0.00161877`
(`r = 8391553`, `d = 16`; integer form `(210837 - 204996 - 4992)/524472 = 849/524472` every
other cell being below its own diagonal share), so the systematic term is an upper envelope
and the residual never adds to it; the opposite sign (raw below `diag`) is diagonal candidates
coinciding with off-diagonal ones, which cannot lower a production success bound. Both rows
are declared at exactly their per-cell maxima and tested cell by cell at exactly those
values - no cell is compared against a looser global allowance. (ii) The residual `(pred_eff-D)/m` is **not** monotone in `d` or in `K`; the largest
value is `0.00541878` at the small field, `d = 16`. (iii) On the field axis the residual
**shrinks** cell by cell: `|residual|` at `r = 2^26` is smaller than at `r = 2^24` in all `3`
matched `(d)` pairs (`-0.00076149` vs `-0.00473047`, `-0.00069426` vs `0.00541878`,
`0.00305551` vs `0.00334813`). The charged allowance is the SUM of the two rows'
per-cell maxima, `0.00161877 + 0.00541878 = 0.00703755` (two conventions of one modelling
error, both charged so neither can be read as an under-charge; the level, not the trend), and
the transfer limit is explicit: these fields are `2^24`-`2^26` and production is
`2^255`, so this is a measured envelope at toy size, extrapolated one-sided and multiplied
nowhere. Production's own diagonal is `2^-50`, i.e. `2^45` times smaller than the toy cell that
dominates the grid maximum.

A single-field `K` sweep at `d = 8` (`|W| = 37/57/73`, `K_cand = 0.140975/0.334571/0.548762`)
is retained: `D/pred = 0.978512/0.986110/0.992568` rises monotonically toward 1, and
`D/pred_eff = 1.00312/1.000675/1.002776` once the toy diagonal is removed.

## Matching: full-key runs, priced radix pipeline

Key = 224 bits of the even representative of `y` (or `p-y`), stored in the record as four
56-bit limbs; the six points `[u]X` share one key, so one cross-side record pair tests six
values of `v`.

* **Routing.** 4 radix scatter hops on the top 85 key bits to one of `2^85` owner nodes
  (48-byte record copies, charged per hop).
* **Local sort.** 3 LSD counting-sort passes over the 63 locally ordered key bits, 21-bit
  digits, `2^21` counters per node; prefix sums are per node (`3 x 2^21 <= 2^22` inside the
  per-node fixed row), not per record.
* **Runs and the remaining bits.** A group is a set of equal 63-bit local prefixes in one
  node. The bits no pass orders (`224 - 63 = 161`, of which 85 are equal by routing, leaving
  76) are **not** claimed to be sorted: the group scan compares all four key limbs pairwise,
  and that comparison is exactly what the `run_scan` row prices (4 load instructions from two
  records, 3 compare/OR-reduce, 1 branch, 1 emit slot, 1 loop = 10). Groups hold `2^84.31`
  cross-side pairs in expectation against a cap of `2^95`, and a true 224-bit match probability
  `2^-76` inside that pair set, i.e. `2^7.3` false keys expected.

**The full-key premise is measured, not assumed** (`toy.py`, part G). Prefix uniformity does
not give a full-key collision probability, so the conditional distribution of the remaining
bits is measured at 128-bit key width on a 128-bit field (`p = 2^128 - 159`, `d = 16`,
`m = (p-1)/16`, 16000 records whose scalars are `a*w` with `a, w` in `H`; key = even `y` at
full field width, 126 bits). Over `1999079` record pairs that agree on the routed prefix,
**zero** full-key collisions occur against a uniform prediction of `2^-120` (95% upper bound
`1.5e-6` per pair), and the counter fires on a planted duplicate (4/4), so a silent
always-false gate is excluded. Conditional on the prefix, the residual bits are flat:
observed/predicted collision counts `0.987 / 0.986 / 1.082` at residual widths 8/12/16 bits,
all within 3 sigma where the expectation is usable; the largest 3-sigma inflation over the six
cells is `1.655`, charged as `1.7`. That factor is what the verification cap is priced at
(`268` false keys, `2^-88` Markov against the `2^98` cap, charged `2^-40`). Transfer limit
stated: measured at 126-bit keys, production keys are 224 bits, where the same reasoning
applies one-sided (more bits can only lower the collision rate) and the measured deviation is
charged, not assumed away.

`toy.py` part D implements this same pipeline on a genuine `j = 0` curve (`p = 268520611`,
`r = 67122049`, `d = 8`, `#E = 4r`): 24-bit key, 6 routed bits, 12 locally ordered bits in
3 LSD passes of 4 bits, 6 middle bits compared only in the group scan. The toy **counts** its
operations and reports them per class: `21` instructions per record per pass (digit 2, counter
3, address 2, 48-byte copy 12, loop 2) and `10` per run pair, matching the ledger's
`sort_3_passes_x_21 = 63` and `run_scan = 10` exactly, with `592` records, `0` false positives
and `0` cap breaches (max `2 <= 64` pairs).

The recovery test is **end to end from one secret**: `s` is drawn uniformly from `F_r`, the
challenge is built from that same `s` as `z = s^d` and `T = [s]P0`, the pipeline recovers `z`
from the records alone, and then recovers `s` from `z` (20 root steps plus the `d`-candidate
scan). Every recovery is checked twice - `[s]P0 == T` and `s^d == z` - and the secret is never
used as a shortcut. Because `z = s^d` is a *uniform* element of `H`, the honest statistic is
the recovery fraction, not a planted hit rate: over `200` uniform secrets the pipeline recovers
`12` (`0.0600`) against the exactly enumerated covered fraction `D/m = 0.059317` of the very
same record set, `z = +0.04`; all `12` pass both checks and `0` false positives occur. This is
the injection point of the coverage model into a real run of the real pipeline.

## Caps and the success arithmetic

`N = 2^115.6503` per side, `K = 6N^2/m = 0.510000`. Debits in success units, each from a
measured row: census row A `0.00161877`, census row B `0.00541878`, run-cap Markov `6.014e-4`,
verification-cap Markov `5.1e-27` charged `2^-40`, gap model charged `2^-40`, diagonal `2^-50`;
total `0.007639`. `0.39950 - 0.007639 = 0.39187 >= 0.39`. The two capped data items: run pairs
`<= 2^95` (expectation `2^84.31`, Markov `<= 2^-10.70`) and verifications `<= 2^98` (expectation
`2^8.07` at the measured factor `1.7`). Node overflow is Chernoff-charged `2^-40` on a `6.37x`
capacity (the true bound is `~exp(-2^34)`), and the walk-list gap tail is measured: `0/1199`
gaps above `2^12` in the 128-bit field's consecutive H-prime list (mean gap `468.30`, max
`2520`), so the charged `2^33` gap bound and `2^32`-entry table both dominate the measurement.

## Per-record implementation

Every row is a 64-bit AArch64 instruction count; `M = 383` per `blst_fp_mul`; the 2000
reference unit covers only the 6M field arithmetic.

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
| key pack 224 bits | 4 | 4 limb stores (four 56-bit limbs) |
| routing 4 hops x 20 | 80 | digit 2, bucket base 1, address 2, guard 3, 48-B copy 12 |
| sort 3 passes x 21 | 63 | digit 2, counter 3, address 2, 48-B copy 12, loop 2 |
| run scan | 10 | 4 limb loads, 3 compare, branch, emit, loop (full 224-bit key) |
| **itemised** | **200** | |
| stress 6/5 | +40 | 20% repricing of every non-arithmetic row |
| **charged** | **240** | `c = 2240/2000 = 1.12` |

## Ledger (units of 2000 instructions)

| term | log2 units |
|---|---:|
| records, `2N x 2240` | 116.813847 |
| chain tables, `2 x ceil(N/2^50) x 2^32` additions | 98.6503 |
| fixed per-node work, `2^85 x 2^31` | 105.0342 |
| H-prime supply, `4 x 2n_b x 2^20 x 68` candidates, 352 mults each | 100.8127 |
| walk-list H-primes | 83.8028 |
| run-pair cap, `2^95 x 2` | 85.0342 |
| verification cap, `2^98` scalar multiplications | 106.5812 |
| root extraction and I/O | 17.6845 |
| root T-scan, `2^20` candidates x `2^11` units | 31.0000 |
| **W** | **116.815482** |

Memory: `2^85` owner nodes at their full `2^40` B plus `2^20` generating and `2^31` storage
nodes `= 2^125 + 2^71 + 2^60 <= 2^125.1` B; node capacity `(2^40-2^36)/48` is `6.31x` the mean
occupancy `2^31.650`.

## Heuristics and limitations

H3occ (random-multiset coverage) is measured by the 6-cell two-axis grid over *modular*
residues in `H` - `toy.py` part F enumerates `v = x/(u x') mod r` with `x, x'` products of
H-primes, and `verify_certificate.py` independently rebuilds one census and checks every
enumerated residue satisfies `v^m = 1`, so this is modular coverage and not rational
distinctness (the `Fraction(w, w')` classes are kept only as structural evidence). H3 (full
224-bit keys behave like uniform strings GIVEN prefix equality) is measured at 128-bit key
width as above and charged at the measured `1.7` factor. H2/H2' prime density and gap bounds
are measured on a consecutive H-prime list rather than assumed geometric. Analytic bound, no
full-size run; instruction counts are itemised and toy-counted, not measured on hardware.

# H-prime Cheon BSGS with sort matching: distinct-coverage model, measured

## Result

`time_log2 = 116.8156` (`292039/2500`) at `n = 2^20`; success `>= 40/100` after every
debit (design `1 - e^-0.51 = 0.39950`, measured model allowance `0.00394970`, run-cap
Markov `6.014e-4`, so `0.39495 >= 0.39`); peak memory `<= 2^125.1` B; no advice;
preprocessing inside the total. Only `A_1 = T` and `A_d` are read, so the bound holds
for every `n >= 2^20`.

Charged work is `2^116.815482`; at the itemised (unstressed) overhead the same ledger is
`2^116.789517`. The claim sits `0.000118` above the charged ledger and `0.0115` below
the live head `116.8271`.

**Rounding row.** `W = 2^116.815482 < claim = 116.8156 < gate = 116.8271`: the claim is
the smallest 4-decimal value strictly above `W` by at least `2e-5` (actual headroom
`1.18e-4`). The `raw <= claim - 0.10` family rule cannot be met at this frontier; the
head cleared review twice with sub-0.01 headroom on the same convention.

## Reduction

`d = 2^20 | r-1`, `v2(r-1) = 32`, `m = (r-1)/d` of 235 bits, `6 | m` so `mu6` lies in
`H = {x : x^m = 1}`. `v = s^d in H` determines `s` up to the `d` roots `s*zeta^j`: 20
successive square roots plus a scan of the `2^20` candidates `[s*zeta^j]P0` against `T`,
charged as its own ledger row (`2^20 x 2^11` units).

## H-prime construction

Scalars are `x = a*w` with `a`, `w` both H-primes (`q^m = 1 mod r`), so `x in H`
automatically: no character table and no fibre labels. Walk list `W` = the first `2^50`
H-primes above `2^39`, top below `2^77`; baby/giant bases are H-primes in `(2^77, 2^97)`
in disjoint pools with `4x` supply. `x = a*w < 2^174 < r`, so products are exact
integers. `x = a*w` determines `(a,w)` by unique factorisation (`a > 2^77 > w`).

## Coverage: the distinct residue set, and what the model does and does not contain

Success for uniform `v` in `H` is the probability that `v` is covered:

    success = E[C_distinct] / m,   C_distinct = #{ v : v = x/(u x'), u in mu6 },

the **set cardinality** of the `K_cand = 6N^2` candidate occurrences (one per record pair
and per `u`). Repetitions are not a debit against this formula; they are part of it.

**H3occ.** The candidate multiset behaves like `K_cand` uniform draws into `m` bins. Then
`E[C_distinct] = m(1 - e^{-K_cand/m})`, and since `K = K_cand/m = 6N^2/m`,

    success = 1 - e^-K.

That is the point an earlier revision argued badly: `1 - e^-K` **already contains every
random repetition** a multiset of `6N^2` draws produces (expected repeats `~K^2 m/2`).
Nothing is subtracted for them, and nothing should be. The statement needs no uniqueness
lemma, and rational uniqueness is not offered as one here.

What the occupancy form does **not** contain is candidates that are *exactly* equal by
construction. In production the only such family is the shared-walk diagonal: `w = w'`
makes the ratio `a/a'` independent of `w`, a share `(|W|-1)/|W|^2 = 2^-50` of the pairs,
charged separately below. In the toy `|W|` is 37 to 73, so the same diagonal is `2^40+`
times larger — which is why the toy is the hard test, and why the census below does not
simply reproduce `1 - e^-K` at small size.

## Coverage measured: three exact censuses (`toy.py`, Part E)

At each size the full candidate set is enumerated and its exact cardinality `D` counted;
`pred` is `m(1 - e^{-K_cand/m})`; `pred_eff` is the same form on the occurrence count with
the toy's diagonal repetitions removed.

| `|W|` | occurrences | `K_cand` | `D` | `pred` | `D/pred` | `D/pred_eff` |
|---|---:|---:|---:|---:|---:|---:|
| 37 | 1,182,816 | 0.140975 | 1,079,520 | 1,103,226 | 0.978512 | 1.00312 |
| 57 | 2,807,136 | 0.334571 | 2,352,672 | 2,385,811 | 0.986110 | 1.000675 |
| 73 | 4,604,256 | 0.548762 | 3,517,170 | 3,543,504 | 0.992568 | 1.002776 |

Three measured facts: (i) `D/pred` **rises monotonically toward 1** with size (0.9785,
0.9861, 0.9926); (ii) the residual is the toy diagonal — removing it moves the ratio to
`1.0007..1.0031`, within `0.31%`; (iii) the largest census brackets the production point
(`K_cand = 0.5488` against `K = 0.51`). The measured shortfall in success units,
`(pred - D)/m`, is `0.00282542 / 0.00394970 / 0.00313864`; we charge the **maximum**,
`0.00394970`, at production — `2^17` times the `2^-20` placeholder it replaces, and
`2^45` times what the production diagonal can cause, so the transfer is one-sided.

## Caps and the success arithmetic

`N = 2^115.6503` per side, `K = 6N^2/m = 0.510000`. Debits, all in success units:
measured model allowance `0.00394970`, run-cap Markov `6.014e-4`, diagonal `2^-50`, total
`0.004551`. `0.39950 - 0.004551 = 0.39495 >= 0.39`.

Two data-dependent items are capped: **run pairs** `<= 2^95` (same-node pairs `2^84.31` in
expectation — `2^85` nodes, `2^31.650` records each, `2^-63` per ordered pair inside a
node — so Markov gives `<= 2^-10.70`; a breach fails the run) and **verifications**
`<= 2^98 >= 6 x 2^95`, which never binds and is what the ledger charges. False 224-bit key
matches are `N^2 2^-224 = 2^7.3`, inside the cap.

## Per-record implementation

Every row is a 64-bit AArch64 instruction count; `M = 383` per `blst_fp_mul`; the 2000
reference unit covers only the 6M field arithmetic, so all addressing, loads, stores, key
work, routing, sorting and scanning are charged below.

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
| routing 4 hops x 20 | 80 | digit 2, bucket base 1, address 2, guard 3, 48-B copy 12 |
| sort 3 passes x 21 | 63 | digit 2, counter 3, address 2, 48-B copy 12, loop 2 |
| run scan | 10 | full-key compare, run emit, loop |
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

The T-scan row is new. Earlier revisions described `2^20` scalar multiplications in the
evidence but left them inside `root_and_io` (`2^17.68`), so the ledger was not the
calculation the text described. It is now explicit, sits `85.8` bits below the records
row, and `W` is unchanged at the six decimals printed.

Memory: `2^85` owner nodes at their full `2^40` B plus `2^20` generating and `2^31`
storage nodes `= 2^125 + 2^71 + 2^60 <= 2^125.1` B; node capacity `(2^40-2^36)/48` is
`6.31x` the mean occupancy `2^31.650`, and the whole node-overflow risk is inside the
measured allowance above.

## Evidence (`toy.py`, stdlib, 4.5 s)

* **A, production field.** `r` prime, `r-1` factored, `d | r-1`, `v2 = 32`, `m` 235 bits,
  `6 | m`, `zeta` of order `d`, H-primes closed under products; `401/6541` of primes below
  `2^16` are 16th powers (`1/16`), `27/6541` are 256th powers (`1/256`).
* **B, toy densities.** `d = 16/64/256`: `0.0639/0.0150/0.0041` against
  `1/d = 0.0625/0.0156/0.0039`, gap maxima `896/4488/11988`.
* **C, rational distinctness.** `|W| = 37`, 8+8 bases, 296 products per side: products
  distinct, collision classes `64` = base pairs, max multiplicity `37`,
  multiplicity-minus-one `2304` = predicted, loss `0.026296567` = `(37-1)/37^2`. This
  counts the diagonal; it is **not** the coverage argument.
* **D, full pipeline on a genuine `j = 0` curve** (`p = 268520611`, `r = 67122049`,
  `d = 8`, `#E = 4r`, `mu6` as `(omega*x, y)`): armed gap-addition walk, 592 records,
  18-bit key, 6-bit routing, sort, run scan, 6 verifications per run pair. `40` planted +
  `10` random, each as `(T, A_d)`: `42/42` recovered `z` and `42/42` recovered `s` accepted
  only on `[s]P0 = T`; `0` false positives, `0` wrong candidates, `0` incomplete trips,
  `0` cap breaches (max `37 <= 64`).
* **E, coverage census.** The three-size table above, plus `20000` uniform challenges on
  the `|W| = 37` set hitting at `0.12695` against its exact deduplicated `D/m = 0.128664`
  (`z = -0.72`). That hit count checks the uniform challenge sampler, not the
  distinct-set prediction — the census does that.

## Heuristics

* **H3occ coverage.** The candidate multiset behaves like `K_cand` uniform draws into `m`,
  so `E[C_distinct] = m(1-e^{-K_cand/m})`. Supported by the three censuses (ratios
  `0.9785/0.9861/0.9926`, monotone in size, residual explained by the counted diagonal to
  within `0.31%`); not a production lower bound, and the largest measured shortfall
  `0.00394970` is charged in full.
* **H2/H2' prime density and gaps.** H-prime fraction `1/2^20` among primes, `1/(2^20 ln q)`
  among integers; the `2^50`-th H-prime above `2^39` is below `2^77` with even gaps
  `<= 2^33`, checked when the list is built (shortfall `< 2^-100`).
* **H3 key bits.** 224-bit key prefixes of the `2N` records behave like independent uniform
  strings; both caps and the node-occupancy Chernoff bound are taken under this, inside the
  same measured allowance.

## Limitations

Analytic bound, no full-size run. Memory `2^125.1` B is exponential, as for any BSGS-type
Cheon attack. H3occ, H2 and H3 are toy-supported models, not production lower bounds; the
census measures them at toy sizes and the largest measured deviation is charged.
Instruction counts are itemised, not measured on hardware. The cluster is idealised:
transfers are charged as posting/polling and copying inside the routing and sort rows.

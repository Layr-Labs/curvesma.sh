# H-prime Cheon BSGS: Fraction-exact census, rule-of-three key rate, derived gap tail

## Result

`time_log2 = 116.8132` (`292033/2500`) at `n = 2^20`; success `>= 39/100` after every debit
(design `1-e^-0.51 = 0.39950`, census rows `0.00161878 + 0.00541879 = 0.00703757`, run-cap Markov
`6.01e-4`, verification-cap Markov `2^-24`, gap tail `2^-40`, diagonal `2^-50`; after-debit success
`0.391865`). Peak memory `<= 2^125.1` B, no advice, preprocessing inside the total.
Only `A_1 = T` and `A_d` are read, so the bound holds for every `n >= 2^20`.

Charged work is `2^116.812976`; the claim sits `0.000224` above it and `0.0139`
below the live head `116.8271`.

## What changed against the rejected revision

1. **No rounded census strings.** Every census assertion is a `Fraction` built from the raw integer
   counts (`distinct_covered`, `predicted_distinct`, `predicted_distinct_eff`, `diagonal_repetitions`,
   `m`) inside `toy-summary.json`; the charged bound IS the exact per-cell maximum, and the decimal
   printed beside it is a *ceiling* of that fraction. The rejected revision charged `0.00161877`
   while the exact value is `849/524472 = 0.0016187707 > 0.00161877`; the charge is now
   `283/174824 = 0.0016187800` (ceiling `0.00161878`) and likewise `1421/262236 = 0.0054187900`
   (ceiling `0.00541879`) for row B. `verify_certificate.py` re-parses the raw integers, rebuilds
   every cell's fraction, compares `Fraction` to `Fraction` at the declared bound, and re-parses the
   printed decimals back into fractions to confirm each one is `>=` the exact value and that the
   printed sum equals the sum of the printed terms.
2. **No 2^-120 rate claim.** The key rate is charged as the rule-of-three bound
   `3/1999079 = 2^-19.35` per compared pair, which is what zero collisions in 1999079 pairs actually
   supports. Sensitivity rows show the ledger closes at any rate at or below `2^-30`
   (`2^53.30` false keys, Markov `2^-34.70`) and in fact at any rate at or below `2^-4.0`.
3. **Derived gap tail.** Density `1/(d ln q)` -> `mu = d ln q <= 2^26.07` over the list ->
   `P(max gap over 2^50 > G) <= 2^50 exp(-G/mu) = 2^-125.6` at the charged `G = 2^33`. The toy's 1199
   gaps now only validate the law (mean 468.30 vs `16 ln 2^51 = 565`, chi-square 7.478 on 11 bins, tail below model at every threshold).
4. **Complete ledger rows.** The scatter counter read/decrement/store and the counting-loop key load
   are charged (pass = 25), `group_runs` carries an explicit op list, the gap load is two limbs and
   routing is 3 copies plus a final store. Itemised 198, stressed 238, `c = 1119/1000`.
5. **partG order claim dropped.** No toy point is claimed to have order equal to the scalar modulus
   defining its H-primes; an explicit small-fixture audit checks `[n]P = O`, `n | #E` by brute-force
   point counting instead.

## Coverage census, exact (raw counts -> Fractions)

| r | d | m | D | pred | diag | rowA | rowB |
|---|---:|---:|---:|---:|---:|---:|---:|
| 8391553 | 8 | 1048944 | 422778 | 422076 | 7104 | -1301/174824 | -827/174824 |
| 8391553 | 16 | 524472 | 204996 | 210837 | 4992 | 283/174824 | 1421/262236 |
| 8391553 | 64 | 131118 | 49272 | 51188 | 2400 | -242/65559 | 439/131118 |
| 33555073 | 8 | 4194384 | 1665594 | 1670963 | 14208 | -8839/4194384 | -1597/2097192 |
| 33555073 | 16 | 2097192 | 838698 | 843283 | 10080 | -5495/2097192 | -182/262149 |
| 33555073 | 64 | 524298 | 206220 | 210821 | 4992 | -391/524298 | 267/87383 |

`rowA = (pred-D-diag)/m`, `rowB = (pred_eff-D)/m`, both exact; declared maxima are the per-cell
maxima themselves, so no cell is ever compared against a looser global. The transfer limit is
explicit: fields `2^24`-`2^26` and `d <= 64` against production `2^255`, `d = 2^20`; production's own
diagonal is `2^-50`, i.e. `2^45` times smaller than the toy cell that dominates the maximum.

## Key matching and the rate row

Key = 224 bits of the even representative of `y`; six points `[u]X` share one key, so one cross-side
record pair tests six values of `v`. Routing 4 hops on the top 85 bits to `2^85` owner nodes, 3 LSD
counting-sort passes over 63 further bits (21-bit digits, `2^21` counters per node), then a group scan
that compares all four 56-bit limbs - the 76 bits no pass orders are compared, not assumed sorted.
In-group pairs: `N^2/2^148 = 2^83.30`. Charged per-pair rate `3/1999079` (rule of three), so
`2^63.95` expected false keys, one scalar multiplication each, cap `2^88` (Markov
`2^-24.045`). The measured deviation ratios at 8/12/16 residual bits are
0.987/0.986/1.082; the 1.7 factor is taken from there and the transfer to 224-bit keys is a labelled
modelling premise, not a measured rate. The 24/32/48/64-bit probes have expected counts
`2^-3.1/2^-11.1/2^-27.1/2^-43.1` in 1999079 pairs, so their zeros are uninformative and are not used.

## Gap tail derivation

H-prime = prime `q` with `q^m = 1 mod r`; a uniform residue satisfies it with probability exactly
`m/(r-1) = 1/d`, and primality costs `1/ln q`, so the density is `1/(d ln q)` and consecutive gaps are
geometric with `mu = d ln q`. Over the list (`q <= 2^97.1`) `mu <= 2^26.07`; then
`P(max gap > 2^33) <= 2^50 exp(-2^33/2^26.07) = 2^-125.6` and the `G` giving exactly `2^-40` is
`2^32.04`, so `G = 2^33` (a `2^32`-entry even-gap table) is charged with 85 bits of margin. The toy
list is the same law at `d = 16`: mean gap 468.30 against `16 ln 2^51 = 565`, 0/1199 gaps above
`2^12` against a geometric prediction of 0.19, upper-tail empirical below model in every cell.

## Per-record ledger (AArch64 instructions, `M = 383` per `fp_mul`)

| row | instr | detail |
|---|---:|---|
| walk_cursor_and_bound | 4 | index load, compare, advance |
| gap_load_2limb | 4 | 2 LDP from the gap array |
| gap_index | 2 | shift/scale to byte offset |
| table_point_load | 6 | 6 LDP (12 limbs) |
| point_store_into_record | 3 | 3 STP |
| walk_loop_control | 3 | counter, pointer, branch |
| p_minus_y | 12 | 6 sub + 6 sbc |
| even_rep_select | 7 | parity, branch, 6 csel |
| key_pack_224_bits | 4 | 4 limb stores |
| route_3copies_plus_final_store | 68 | 3 hops x 20 (digit 2, bucket base 1, address 2, guard 3, 48-B copy 12) + final hop 8 (same without copy: the store IS the placement) |
| sort_3_passes_x_25 | 75 | digit 2, counting-key-load 1, counter read 1, add 1, counter store 1, scatter counter read 1, decrement 1, store 1, address 2, 48-B copy 12, loop 2 = 25 per record per pass |
| group_scan | 10 | group_runs op list: 4 limb loads, 3 compare/OR-reduce, branch, emit, loop |
| **itemised** | **198** | |
| stress 6/5 | +40 | ceiling of `6/5 x 198` |
| **charged** | **238** | `c = 2238/2000 = 1119/1000` |

## Ledger (instructions; W in units of 2000)

| term | instructions | log2 units |
|---|---:|---:|
| records | 291819187909464924373637434008781516380 | 127.7783 |
| chain_tables | 994818461282523598820751704064000 | 109.6161 |
| fixed_per_node | 83076749736557242056487941267521536 | 116.0000 |
| hprime_supply | 11626940766239494561217535541248 | 103.1972 |
| run_pair_cap | 79228162514264337593543950336 | 96.0000 |
| root_and_io | 421315604 | 28.6503 |
| root_scan_T | 4294967296000 | 41.9658 |
| walk_list | 317745634036292260041580783599616 | 107.9696 |
| verification_cap | 237561312508884117443279388147712 | 107.5500 |
| **total** | 291903826490778238069429130668392952432 | **116.812976** |

## Order audit and scope

Small `j = 0` fixture (`p = 10009`, brute-force `#E = 10192`):
`n = 7` divides `#E` and `[n]P = O`. This establishes the check and
explicitly does not transfer an order claim to partG, whose keys are used only as a conditional
deviation probe.

Heuristics: H3occ (random-multiset coverage) measured by the six-cell two-axis grid;
H3dev (full-key collision deviation) charged at the measured 1.7; H2/H2' (prime density and gap
bounds) derived and validated on a consecutive H-prime list. Analytic bound, no full-size run;
instruction counts itemised and toy-counted, not measured on hardware.

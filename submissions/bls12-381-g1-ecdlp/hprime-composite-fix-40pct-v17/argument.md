# H-composite Cheon BSGS: corrected composite order, certified collision UCL, declared gap premise

## Result

`time_log2 = 116.8132` (`292033/2500`) at `n = 2^20`; success `>= 39/100` after every
debit (design `0.399504`, census `0.00703757`, run-pair Markov `6.014e-04`, verification Markov
`2^-25`, gap-premise allowance `2^-40`, diagonal `2^-50`; after-debit `0.391865`). Peak memory
`<= 2^125.1` B, no advice, preprocessing inside the total. Only `A_1 = T` and `A_d` are read, so
the bound holds for every `n >= 2^20`. Charged work `2^116.812976`; the claim sits `0.000224` above
it and `0.0139` below the live gate `116.8271`.

## Order structure, corrected

`r - 1 = 2^32 * u` (`u` odd) with `d = 2^20`, so `m = (r-1)/d = 2^12 * (odd)`: `v2(m) = 12`, `6 | m`,
and `m` is composite. `H = {x : x^m = 1}` is the `d`-th-power subgroup: cyclic of order `m`, index
`d`, with `mu6 subset H` because `6 | m`. Consequences priced explicitly:

- `v = s^d` determines `s` up to the kernel `{s zeta^j}` of size `d`; that is the only role of the
  order, and it needs no primality.
- the `d`-th root is obtained by 2-Sylow lifting (`u^{-1} mod 2^32` exists as `u` is odd, and 20
  square-root-and-adjust steps inside the 2-Sylow subgroup), charged in `root_and_io`;
- `gcd(d, m) = 2^12`, so `d` has **no** inverse modulo `m` and the package makes no
  "`v^{1/d} mod m`" step; candidate selection is the `2^20`-point scan `[s zeta^j]P0 = T`, charged
  in `root_scan_T`;
- the only primality assertions left are about the integer walk elements `q` (`q` prime and
  `q^m = 1 mod r`), and the density `1/(d ln q)` comes from `m/(r-1) = 1/d` for the test
  `q^m = 1`, never from the `m`-th-power density `1/m` (unused subgroup).

## Coverage census, exact, with a declared scaling transfer

| r | d | m | D | pred | diag | rowA | rowB | absB*sqrt(m) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 8391553 | 8 | 1048944 | 422778 | 422076 | 7104 | -0.007441770 | -0.004730472 | 4.845 |
| 8391553 | 16 | 524472 | 204996 | 210837 | 4992 | 0.001618771 | 0.005418783 | 3.924 |
| 8391553 | 64 | 131118 | 49272 | 51188 | 2400 | -0.003691331 | 0.003348129 | 1.212 |
| 33555073 | 8 | 4194384 | 1665594 | 1670963 | 14208 | -0.002107342 | -0.000761494 | 1.560 |
| 33555073 | 16 | 2097192 | 838698 | 843283 | 10080 | -0.002620170 | -0.000694262 | 1.005 |
| 33555073 | 64 | 524298 | 206220 | 210821 | 4992 | -0.000745759 | 0.003055514 | 2.212 |

`rowA = (pred-D-diag)/m`, `rowB = (pred_eff-D)/m`, exact. Declared maxima `283/174824` and
`1421/262236`, charged sum `0.00703757`; every cell is tested at exactly the declared bound
as a fraction of its own raw integers. The transfer argument is now quantitative in form:
`|rowB| * sqrt(m) <= 4.845` across both axes (field size and `d`), the deviation of the
distinct count from the multiset expectation is `O(1/sqrt(m))`, and the extrapolation to
`m = 2^234` gives `2^--115.15` - 100+ bits below the charge. The scaling premise
(H3occ-scale) is declared, and the charge is the toy maximum, so nothing depends on it.

## Key matching and the certified rate

Key = even representative of `y`, truncated to 224 bits; six `[u]X` share one key, so one cross-side
pair tests six `v`. Routing 4 hops on the top 85 bits to `2^85` owners, 3 LSD counting-sort passes
over 63 further bits, then a group scan comparing all four 56-bit limbs. In-group pairs
`N^2/2^148 = 2^83.3007`.

The rate is measured by `partK`: independent cross-side pairs, joint-disjoint H-prime pools and
walk lists (no shared record in any pair), declared even-representative key, positive control
firing. Zero full-key collisions in `n = 4406549` (`2^22.0712`) prefix-equal pairs. Charged
rate = certified one-sided 95% UCL `3/n = 2^-20.4863`, valid since `(1-3/n)^n <= exp(-3) < 1/20`.
Then `2^62.8144` expected false keys, one scalar multiplication each, cap `2^88`, Markov
`2^-25.186` (charged `2^-25`). The ledger closes at any rate at or below `2^-30`, in fact at any
rate at or below `2^-6`.

| rate charged | log2 rate | log2 false keys | log2 Markov at 2^88 | closes |
|---|---:|---:|---:|---|
| ucl_measured | -20.4863 | 62.8144 | -25.186 | yes |
| ten_x | -23.8082 | 59.4925 | -28.507 | yes |
| 2^-30 | -30.0000 | 53.3007 | -34.699 | yes |
| 2^-40 | -40.0000 | 43.3007 | -44.699 | yes |
| 2^-50 | -50.0000 | 33.3007 | -54.699 | yes |
| 2^-60 | -60.0000 | 23.3007 | -64.699 | yes |

Residual deviations from the same experiment (the entire basis of the 1.7 factor in H3dev):

| residual bits | observed | predicted | ratio | z |
|---|---:|---:|---:|---:|
| 8 | 16991 | 17213.082 | 0.987098 | -1.693 |
| 12 | 1089 | 1075.818 | 1.012253 | 0.402 |
| 16 | 63 | 67.239 | 0.936962 | -0.517 |
| 24 | 2 | 0.263 | 7.614674 | 3.39 |
| 32 | 0 | 0.001 | n/a | -0.032 |
| 48 | 0 | 0.000 | n/a | -0.0 |
| 64 | 0 | 0.000 | n/a | -0.0 |

## Gap premise (declared, three-size fit)

`H3gap`: consecutive gaps on the **fixed** production H-prime list follow the geometric law with
mean `mu = d ln q`; the list is public, so a verifier can check the true maximum gap at the charged
sieving cost. Then `mu <= d ln(2^97.1) = 2^26.07` and `P(max gap over 2^50 > 2^33) <= 2^-125.6`
(the `G` giving exactly `2^-40` is `2^32.04`). Fit at three sizes:

| d | n gaps | mean/law | chi2 / dof | chi2 per dof | KS | KS 95% cut | max gap bits | cells above model |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 1200 | 1.0292 | 3.60/10 | 0.360 | 0.030308 | 0.039260 | 12 | 1/5 |
| 16 | 1200 | 0.9645 | 10.08/10 | 1.008 | 0.022552 | 0.039260 | 12 | 1/5 |
| 64 | 1200 | 0.9846 | 4.21/10 | 0.421 | 0.020047 | 0.039260 | 15 | 1/5 |

The `2^-40` debit is a modelling allowance against the premise, not a probability over challenges
or coins; the charged `G = 2^33` carries 85 bits of margin inside the model.

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
| route_3copies_plus_final_store | 68 | 3 hops x 20 + final hop 8 (the store IS the placement) |
| sort_3_passes_x_25 | 75 | 25 per record per pass, incl. scatter counter read/decrement/store |
| group_scan | 10 | 4 limb loads, 3 compare/OR-reduce, branch, emit, loop |
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

Small `j = 0` fixture (`p = 10009`, brute-force `#E = 10192`): `n = 7` divides `#E` and `[n]P = O`;
no transfer of an order claim to partG/partK is made. Named premises: `H3occ` and `H3occ-scale`
(coverage, its six-cell measurement, and the `1/sqrt(m)` scaling), `H3dev` (full-key deviation,
126 -> 224 bit transfer), `H3gap` (geometric gap law on the fixed list). Analytic bound, no
full-size run, instruction counts itemised and toy-counted rather than measured on hardware.

The certificate recomputes the submitted ledger and the declared model formulas from raw evidence;
it establishes arithmetic consistency of the declared premises and does **not** independently
establish production resource bounds, run the production attack, or verify the premises.

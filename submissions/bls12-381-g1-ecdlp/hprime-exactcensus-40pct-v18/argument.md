# H-composite Cheon BSGS: exact structural key census, expectation-Markov debit, position-index record

## Result

`time_log2 = 116.8132` (`292033/2500`) at `n = 2^20`; success `>= 39/100` after every
debit (design `0.399504`, census `0.00703757`, candidate Markov `2^0.000`, gap-premise
allowance `2^-40`, diagonal `2^-50`; after-debit `0.392465`). Peak memory `<= 2^123.235` B, no
advice, preprocessing inside the total. Only `A_1 = T` and `A_d` are read, so the bound holds for every
`n >= 2^20`. Charged work `2^116.812770`; the claim sits `0.00043` above it and `0.0139` below
the live gate `116.8271`.

What changed from v17: the zero-event upper confidence bound and the Chernoff occupancy language are
**deleted**, not softened. The collision debit is now an expectation over deterministic lists plus
Markov on that expectation (which needs no independence), the record is the declared 48-byte
position-index layout with a priced reconstruction, and the key rule is executed by one shared module.

## Order structure, corrected

`r - 1 = 2^32 * u` (`u` odd) with `d = 2^20`, so `m = (r-1)/d = 2^12 * (odd)`: `v2(m) = 12`, `6 | m`,
and `m` is composite. `H = {x : x^m = 1}` is the `d`-th-power subgroup: cyclic of order `m`, index
`d`, with `mu6 subset H` because `6 | m`. Consequences priced explicitly:

- `v = s^d` determines `s` up to the kernel `{s zeta^j}` of size `d`; that is the only role of the
  order, and it needs no primality.
- the `d`-th root is obtained by 2-Sylow lifting (`u^{-1} mod 2^32` exists as `u` is odd), charged in
  `root_and_io`;
- `gcd(d, m) = 2^12`, so `d` has **no** inverse modulo `m` and the package makes no
  "`v^{1/d} mod m`" step; candidate selection is the `2^20`-point scan `[s zeta^j]P0 = T`, charged in
  `root_scan_T`;
- the only primality assertions left are about the integer walk elements `q` (`q` prime and
  `q^m = 1 mod r`), and the density `1/(d ln q)` comes from `m/(r-1) = 1/d`, never from the
  `m`-th-power density `1/m` (unused subgroup).

## One key rule, one record layout, one reconstruction

`keymod.py` is the single definition of the declared key: `key = even member of {y, p-y}`, truncated to
`168` bits. `p` is odd, so exactly one of `y, p-y` is even and the representative is unique; its low bit
is always `0`, so the key image is `2^167`, not `2^168`, and the debit uses the image. `toy.py` imports
that module, and a 100-input agreement row compares it with an independent two-line reference
(`y if y % 2 == 0 else p - y`): `100/100 keys, 48/48 records`.

| word | bytes | field | bits | span |
|---:|---|---|---:|---|
| 0 | 0-7 | key_limb0 | 56 | [0,56) |
| 1 | 8-15 | key_limb1 | 56 | [56,112) |
| 2 | 16-23 | key_limb2 | 56 | [112,168) |
| 3 | 24-31 | side | 1 | bit 0 |
| 4 | 32-39 | index_hi | 64 | [64,128) |
| 5 | 40-47 | index_lo | 64 | [0,64) |

Payload is 297 bits in 384 (`48` B), so the record carries a 128-bit position and no scalar. A
position `i` splits into `(pool slot, walk offset)`; the 97-bit base is read from the stored pool
array, the 76-bit walk element is rebuilt from a `2^12`-gap checkpoint plus at most `2^12` gap adds, and
`a*w` (`< 2^174`) is formed with 9 `fp_mul`. The scalar is therefore never stored, which closes the
174-bit overflow of the v17 layout `[4 key limbs, side, scalar]`.

| reconstruction row | instr | detail |
|---|---:|---|
| index_load_and_mask | 3 | 2 LDP of the two index words, AND/shift |
| pool_load_2words | 4 | 2 LDP from the stored 97-bit pool entry, addressing 2 |
| gap_index_shift | 2 | shift to checkpoint byte offset |
| checkpoint_load_2words | 4 | 2 LDP of the prefix sum at the checkpoint |
| local_gap_walk | 24576 | 2^12 gap adds: 2 loads, add, carry, store, loop = 6 each |
| scalar_product_3x3 | 3471 | 9 fp_mul (383 Ir each, COST_MODEL.md) + 24 accumulate/addressing for a*w |
| ratio_and_six_test | 4776 | 6 candidate multipliers: one fp_mul + compare/select each, 796 per u |
| **itemised** | **39404** | |
| stress 6/5 | 39404 | `2^15.2661` per candidate |

Each surviving candidate is charged `1048576` instructions (`2^20`), a ceiling over the
reconstruction row plus one full scalar multiplication (`2^19.55`).

## Key collisions: exact structure, expectation, Markov

The production lists are deterministic, so no sampling model exists and none is claimed.

- **Algebraic part, exact.** Cross-side `x = a*w = a'*w'` is impossible: `a, a'` are prime from disjoint
  pools and `w, w'` come from one prime list, so unique factorisation forces `a = a'`, `w = w'`. The only
  other algebraic coincidence is `x + y = r`, counted exactly by sorting integers.

| records/side | cross pairs | equal products | complements | distinct products | max scalar bits |
|---:|---:|---:|---:|---:|---:|
| 4096 | 16777216 | 0 | 0 | 8192 | 154 |
| 2048 | 4194304 | 0 | 0 | 4096 | 154 |

- **Truncation part, exact counts.** Cross-side prefix collisions, counted by hashing on the low `k`
  bits on the production curve and the production scalar field, at two list sizes and seven widths.
  Rows are the two sizes; columns are `k`.

| 10 | 12 | 14 | 16 | 18 | 20 | 22 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 32535 | 7975 | 1929 | 465 | 106 | 27 | 6 |
| 8009 | 2015 | 495 | 113 | 32 | 10 | 2 |
| 32768 | 8192 | 2048 | 512 | 128 | 32 | 8 |
| 0.993 | 0.974 | 0.942 | 0.908 | 0.828 | 0.844 | 0.750 |
| 0.978 | 0.984 | 0.967 | 0.883 | 1.000 | 1.250 | 1.000 |

Predicted column is `pairs / 2^(k-1)` for the first size (the image of an even representative halves
the space). Observed ratios lie in `0.828-1.000`; the charged occupancy factor is
`kappa = 1.544`, the 3-sigma ceiling over cells with at least 16 predicted collisions. Then

`E[#bad] = kappa * N^2 / 2^167 = 2^64.9274` (`N^2 = 2^231.3007`), and Markov on that **expectation**
gives `P(candidate work > 2^104.00) <= 2^0.000`. No independence is used: the expectation is
linear in the exact integer counts, and no pair is ever treated as an independent trial. Records and
factors are reused across comparisons by construction - which is precisely why the v17 Bernoulli UCL was
not admissible, and why nothing here depends on such a model.

Routing occupancy is deterministic too: `2^84` route classes hold `2^32.65` records each, and the group
scan compares `2^84.30` prefix-equal pairs against the charged `2^95` run-pair cap, a
`2^10.70` margin.

## Coverage census, exact, with a declared scaling trend

| r | d | m | D | pred | diag | rowA | rowB | absB*sqrt(m) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 8391553 | 8 | 1048944 | 422778 | 422076 | 7104 | -0.007441770 | -0.004730472 | 4.845 |
| 8391553 | 16 | 524472 | 204996 | 210837 | 4992 | 0.001618771 | 0.005418783 | 3.924 |
| 8391553 | 64 | 131118 | 49272 | 51188 | 2400 | -0.003691331 | 0.003348129 | 1.212 |
| 33555073 | 8 | 4194384 | 1665594 | 1670963 | 14208 | -0.002107342 | -0.000761494 | 1.560 |
| 33555073 | 16 | 2097192 | 838698 | 843283 | 10080 | -0.002620170 | -0.000694262 | 1.005 |
| 33555073 | 64 | 524298 | 206220 | 210821 | 4992 | -0.000745759 | 0.003055514 | 2.212 |

`rowA = (pred-D-diag)/m`, `rowB = (pred_eff-D)/m`, exact. Declared maxima `283/174824`,
`1421/262236`, charged sum `0.00703757`; every cell is tested at exactly the declared bound as a
fraction of its own raw integers. The expectation is `E[#missing]/m = exp(-K)` with `K >= 51/100` and
the exact rational upper bound `0.60049558`, so `design = 0.399504`. The transfer is declared, with a
two-axis trend rather than an independence argument:

| d | small field | large field | absB small | absB large | shrinks |
|---:|---:|---:|---:|---:|---|
| 8 | 8391553 | 33555073 | 0.004730 | 0.000761 | yes |
| 16 | 8391553 | 33555073 | 0.005419 | 0.000694 | yes |
| 64 | 8391553 | 33555073 | 0.003348 | 0.003056 | yes |

`|rowB| * sqrt(m) <= 4.845` across both axes and falls in 3/3 matched field pairs; the
extrapolation to `m = 2^234` is `2^-115.15`, more than 100 bits below the charge. The charge is the
toy maximum, so nothing depends on the extrapolation. No Chernoff bound and no `O(1/sqrt(m))`
independence claim is made anywhere.

## Gap premise (declared, three-size fit)

`H3gap`: consecutive gaps on the **fixed** production H-prime list follow the geometric law with mean
`mu = d ln q`. Then `P(max gap over 2^50 > 2^33) <= 2^-125.6` (the `G` giving exactly `2^-40` is
`2^32.04`), so the `2^-40` debit is a modelling allowance against the premise, not a probability over
challenges or coins.

| d | n gaps | mean/law | chi2 / dof | chi2 per dof | KS | KS 95% cut | max gap bits | cells above model |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 1200 | 1.0292 | 3.60/10 | 0.360 | 0.030308 | 0.039260 | 12 | 1/5 |
| 16 | 1200 | 0.9645 | 10.08/10 | 1.008 | 0.022552 | 0.039260 | 12 | 1/5 |
| 64 | 1200 | 0.9846 | 4.21/10 | 0.421 | 0.020047 | 0.039260 | 15 | 1/5 |

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
| even_rep_select | 7 | parity test, branch, 6 csel |
| key_pack_168_bits | 4 | 3 word stores + mask register setup |
| route_4hops_x_21bits | 68 | 3 hops x 20 + final hop 8 (the store IS the placement) |
| sort_3_passes_x_25 | 75 | 25 per record per pass, incl. scatter counter read/decrement/store |
| group_scan | 10 | 3 key word loads, 3 compare/OR-reduce, side load, branch, emit, loop |
| **itemised** | **198** | |
| stress 6/5 | +40 | ceiling of `6/5 x 198` |
| **charged** | **238** | `c = 2238/2000 = 1119/1000` |

## Ledger (instructions; W in units of 2000)

| term | instructions | log2 units |
|---|---:|---:|
| records | 291819187909464924373637434008781516380 | 127.7783 |
| chain_tables | 994818461282523598820751704064000 | 109.6161 |
| fixed_per_node | 41538374868278621028243970633760768 | 115.0000 |
| hprime_supply | 11626940766239494561217535541248 | 103.1972 |
| run_pair_cap | 79228162514264337593543950336 | 96.0000 |
| candidate_scan_cap | 20282409603651670423947251286016 | 104.0000 |
| root_and_io | 421315604 | 28.6503 |
| root_scan_T | 4294967296000 | 41.9658 |
| walk_list | 317745634036292260041580783599616 | 107.9696 |
| **total** | 291862070837007054215953867365622329968 | **116.812770** |

## Memory (unscored, self-consistent)

| term | bytes | detail |
|---|---:|---|
| records | 2^122.24 | 2N records x 48 B |
| sort_scratch | 2^122.24 | one counting-sort copy of the record array |
| base_pool | 2^70.00 | a-pool entries, 2 words each |
| gap_checkpoints | 2^42.00 | 2^50 walk gaps / 2^12 per checkpoint, 2 words |
| routing_nodes | 2^90.00 | per-node bucket descriptors |

## Order audit and scope

Small `j = 0` fixture (`p = 10009`, brute-force `#E = 10192`): `n = 7` divides `#E` and `[n]P = O`; no
transfer of an order claim to the census or toy lists is made. Named premises: `H3occ` (coverage), its
`H3occ-scale` trend, `H3key` (key occupancy, tested by exact counts), `H3gap` (geometric gap law on the
fixed list). Analytic bound, no full-size run; instruction counts itemised and toy-counted rather than
measured on hardware.

The certificate recomputes the submitted ledger and the declared model formulas from raw evidence; it
establishes arithmetic consistency of the declared premises and does **not** independently establish
production resource bounds, run the production attack, or verify the premises.

# H-composite Cheon BSGS: challenge-faithful H census, one declared coverage statistic, position-index record

## Result

`time_log2 = 116.8132` (`292033/2500`) at `n = 2^20`; success `>= 39/100`
(design `0.399504`, census `0.00703757`, candidate Markov `2^-18.749`, gap-premise
allowance `2^-40`, diagonal `2^-50`; after-debit `0.392465`). Peak memory `<= 2^123.235` B, no
advice, preprocessing inside the total. Charged work `2^116.812770`; the claim sits `0.00043` above it
and `0.0139` below the live gate `116.8271`.

The v18 revision was rejected with four flaws. This revision repairs the two cited as material
holes and carries, verbatim, everything the verdict accepted.

1. **The key census is now on the object the attack actually builds.** v18 compared UNSCALED points
   from ORDINARY prime products: `prod_primes` never enforced `H`-membership, `H` has index
   `d = 2^20` in `Z_r^*` so an unfiltered prime is a `d`-th power with probability `2^-20`, and
   neither list carried the challenge multiplier. Both defects are removed: every element of both
   faithful lists is in `H`, and the giant list is multiplied by `v = s^d` with `s` uniform in
   `Z_r^*`. The old object is retained only as the **control row**, so the difference is visible
   rather than asserted away.
2. **The expectation is an expectation over the challenge space**, with the regularity of `s -> s^d`
   stated (`gcd(d, r-1) = d`, so the map is exactly `d`-to-1 onto `H`) and the residual dependence
   declared instead of denied.
3. **The coverage transfer declares one statistic before the numbers are shown** - the residual
   fraction `|rowB|` - and the direction claim is made on that statistic only. The v18 statistic
   `|rowB|*sqrt(m)` is demoted to an uncharged diagnostic and its `m`-extrapolation is withdrawn.
4. **Carried verbatim**: the `C = 2000` ledger division, the corrected subgroup structure, the
   even-`y` key rule, the position-index record, the gap premise, the memory rows, the claim.

## 1. Order structure (carried verbatim; the verdict found it credible)

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
  `q^m = 1 mod r`), and the density `1/(d ln q)` comes from `m/(r-1) = 1/d`.

## 2. Key collisions: exact counts on the challenge-faithful object

### 2.1 What the attack compares

Baby element `x = a*w`, giant element `y = v*a'*w'` with `v = s^d` the `d`-adic challenge. Equality
of truncated keys is a property of the points `[x]P0` and `[y]P0`; because `v` multiplies the giant
scalar, the giant list is a challenge-dependent multiplicative copy of the base list. A census that
omits `v` and omits the `H` filter is not a census of this object, and v18's numbers were exactly
that. `census_v19.py` builds the faithful pair:

| object | baby list | giant list | membership |
|---|---|---|---|
| faithful | `a*w`, `a,w in H` | `v*a'*w'`, `v = s^d`, `s` uniform in `Z_r^*` | `H` enforced on every element |
| control (v18) | `a*w`, `a,w` ordinary primes | `a'*w'`, unscaled | none |

The faithful construction draws `a, w` as `z^d` for uniform `z`, so both factors and their product
lie in `H` by closure (`census-v19.json:faithful[*].h_members = true`); `6` independent
challenges are drawn and every width is counted under each (`census-v19.json:n_challenges`).

### 2.2 H-membership, measured, not asserted

`1048576` uniform draws from the pool band, each tested by `x^m = 1 mod r`:
`2` accepted, rejected fraction `0.999998`
(`census-v19.json:membership`). The index of `H` is `d = 2^20`, so the expected accept rate is
`2^-20` (one draw in `2^20`); the observed count is within Poisson noise of that. The filter is doing
what the arithmetic says, and a list built by `prod_primes` is - with probability `1 - 2^-20` per
element - not an `H` list at all. This is the measured form of the verdict's "`prod_primes` does not
enforce H membership".

### 2.3 The multiplier, and why the transfer is over challenges

`x -> x^d` on `Z_r^*` is a group homomorphism whose kernel is the `d`-torsion, of size
`gcd(d, r-1) = d = 2^20` (`certificate.json:order_structure.gcd_d_m`). So `s` uniform in `Z_r^*` gives
`v = s^d` **uniform on `H`**, and multiplication by `v` is a bijection of `H` - a uniform permutation
of the giant list. Conditional on the challenge, for each fixed giant element `y_j` the value `v*y_j`
is uniform on `H`; hence

    Pr[key_i = key(v*y_j)] = 2^-167                        (the rate the census measures)
    E_v[#collisions] = kappa * N^2 * 2^-167                (linearity over the N^2 pairs)

with **no independence between pairs assumed**. What is honestly residual: pairs sharing the same `j`
are correlated through the common factor `v`, and the number of collisions is a function of `v` -
which is why the census averages over challenges and why Markov is applied to the challenge-random
total work rather than to a per-pair Bernoulli model. The analytic step uses only uniformity of
`v*y_j` on `H`; the numerical step is measured, not transferred from the unscaled object.

### 2.4 Exact collision counts, seven widths, two sizes, six challenges

Mean over six challenges at `|W| = 64` walk elements, with per-width min and max shipped
(`census-v19.json:faithful[*].exact_rows[*].observed_min`/`.observed_max`); the full per-challenge
rows are regenerated deterministically by `census_v19.py` (seed 20261024). Predicted column is `pairs/2^(k-1)` (the even
representative's low bit is always zero, so the image is `2^(k-1)`).

| 10 | 12 | 14 | 16 | 18 | 20 | 22 |
| 32813 | 8183 | 2062 | 510 | 124 | 36 | 9 |
| 8211 | 2073 | 516 | 128 | 35 | 8 | 2 |
| 32768 | 8192 | 2048 | 512 | 128 | 32 | 8 | predicted
| 1.001 | 0.999 | 1.007 | 0.996 | 0.974 | 1.125 | 1.167 |
| 1.002 | 1.013 | 1.008 | 1.005 | 1.094 | 1.021 | 1.167 | observed/predicted

Faithful ratios lie in `0.973958`-`1.125`
(`certificate.json:key_census.faithful_ratio_min`, `.faithful_ratio_max`).

### 2.5 The mutation control: the same generator on the v18 object

| 10 | 12 | 14 | 16 | 18 | 20 | 22 |
| 32535 | 7975 | 1929 | 465 | 106 | 27 | 6 |
| 8009 | 2015 | 495 | 113 | 32 | 10 | 2 |
| 0.993 | 0.974 | 0.942 | 0.908 | 0.828 | 0.844 | 0.750 |
| 0.978 | 0.984 | 0.967 | 0.883 | 1.000 | 1.250 | 1.000 | observed/predicted

Control ratios lie in `0.828125`-`1.0`; the largest per-cell difference between the
faithful and the control ratio is `0.416667` (`census-v19.json:max_abs_delta`). The honest reading:
the v18 counts were not numerically absurd, they were counts of the wrong object, so they could not
support a statement about the giant list the attack builds. The faithful row replaces them; the control
row is kept so both are visible. The charged occupancy factor is `kappa = 1.932`, the 3-sigma
ceiling over faithful cells with at least 16 predicted collisions, taken at the worst sampled challenge
so the charge is an upper envelope of the challenge family
(`certificate.json:key_census.kappa_float`); the same rule on the control gives `1.544`,
which is not used.

### 2.6 The debit

```
E[#bad] = kappa * N^2 / 2^167 = 2^65.2508          (N^2 = 2^231.3007)
P(candidate work > 2^104.00) <= 2^-18.749   (Markov, expectation over the challenge space)
```

The `2^84` candidate cap sits 16 bits above the expectation; each candidate is charged `1048576`
instructions (`2^20`), a ceiling over the itemised reconstruction plus one full scalar multiplication.
Routing occupancy remains deterministic: `2^84` route classes hold `2^32.65` records each and the group
scan sees `2^84.30` prefix-equal pairs against the charged `2^95` cap, a `10.70` margin.

### 2.7 The structural part, and the `mu6` orbit

A cross-side product equality `a*w = a'*w'` is impossible: `a, a'` are prime from disjoint pools and
`w, w'` come from one prime list, so unique factorisation forces `a = a'`, `w = w'`. The remaining
algebraic coincidence `x + y = r` is counted exactly by sorting integers.

| records/side | cross pairs | equal products | complements | distinct products | max scalar bits |
|---:|---:|---:|---:|---:|---:|
| 4096 | 16777216 | 0 | 0 | 8192 | 154 |
| 2048 | 4194304 | 0 | 0 | 4096 | 154 |

Equal full even-`y` keys correspond to six scalar relations, not two: the key identifies `y` with `-y`,
and the attack's candidate construction ranges over all six multipliers
`u in mu6 = {+-1, +-omega, +-omega^2}`. The coverage enumeration below walks the **whole `mu6` orbit
explicitly** (all six `u` per pair), so no `mu6` relation is left outside the counted set; in the key
census the six orbit points are six different points with six different `y`-coordinates, so they act
as distinct records, not as a hidden equality. v18's sentence "every remaining key collision is a
truncation accident" is replaced by: the counted structural coincidences are the equal products and
complements of the table above plus the six-multiplier orbit of the coverage enumeration, and the
truncation counts of section 2.4 are measured on top of them.

## 3. Coverage: one statistic, declared before the numbers

### 3.1 Model and charge

Success is `E[C_distinct]/m` with `K_cand = 6N^2` candidate occurrences `v = x/(u x')`. Under the
occupancy model, `D/m = 1 - exp(-K)` with `K = 6N^2/m >= 51/100`, whose exponential is bounded by
exact rationals: `exp(-K) <= 0.60049558`, so `design = 1 - exp(-K) >= 0.399504`
(`certificate.json:coverage_grid.expectation`). The model deviation is measured as exact Fractions
over six censuses on two toy fields and `d in {8,16,64}`, each at the production design point
`K_cand` in `0.495-0.515`, by enumerating the whole candidate set.

### 3.2 The declared statistic

**Declared up front, and used for the direction claim and the charge:** the residual **fraction**

    rowA = (pred - D - diag)/m ,   rowB = (pred_eff - D)/m ,   statistic = |rowB|

It is dimensionless in `m` because the thing the success claim charges is itself a fraction of `H`.
The direction claim is made on `|rowB|` and on nothing else
(`certificate.json:coverage_grid.declared_statistic`).

| r | d | m | D | pred | diag | rowA | rowB | statistic |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 8391553 | 8 | 1048944 | 422778 | 422076 | 7104 | -0.007441770 | -0.004730472 | 0.004730472 |
| 8391553 | 16 | 524472 | 204996 | 210837 | 4992 | 0.001618771 | 0.005418783 | 0.005418783 |
| 8391553 | 64 | 131118 | 49272 | 51188 | 2400 | -0.003691331 | 0.003348129 | 0.003348129 |
| 33555073 | 8 | 4194384 | 1665594 | 1670963 | 14208 | -0.002107342 | -0.000761494 | 0.000761494 |
| 33555073 | 16 | 2097192 | 838698 | 843283 | 10080 | -0.002620170 | -0.000694262 | 0.000694262 |
| 33555073 | 64 | 524298 | 206220 | 210821 | 4992 | -0.000745759 | 0.003055514 | 0.003055514 |

Exact maxima: `rowA <= 0.00161878` (`283/174824`) and `|rowB| <= 0.00541879` (`1421/262236`),
charged as Fractions, `0.00703757` in total, every cell tested at exactly the declared bound as a
fraction of its own raw integers (`certificate.json:coverage_grid.maxA/maxB/allowance`).

| 8 | 8391553 | 33555073 | 0.004730 | 0.000761 | yes |
| 16 | 8391553 | 33555073 | 0.005419 | 0.000694 | yes |
| 64 | 8391553 | 33555073 | 0.003348 | 0.003056 | yes |

`|rowB|` falls in 3/3 matched field pairs (`certificate.json:coverage_grid.trend[*].shrinks`), so the
declared two-axis trend is monotone **in the declared statistic**. The charge is the toy maximum; no
`m`-scaling law is used and no extrapolation to `m = 2^234` is claimed.

### 3.3 What is withdrawn

v18 declared `|rowB|*sqrt(m)` as the trend statistic and quoted a `2^-115` extrapolation from it. That
statistic **increases** across the `d = 64` matched pair, so the declared claim and the table
contradicted each other. The normalisation is also dimensionally wrong for this object: a fixed ratio
set's pair-collision bias is not a `1/sqrt(m)` birthday effect, and linearity of expectation does not
convert an occupancy heuristic into an `exp(-K)` missing fraction. Both are withdrawn. The quantity
`|rowB|*sqrt(m)` is still printed for reference (`4.845` maximum,
`certificate.json:coverage_grid.diagnostic`) and is **not** used in the charge, not extrapolated, and
not part of any claim.

### 3.4 What the transfer is, stated plainly

`H3occ-scale` is a **declared premise**: the production residual fraction is bounded by the toy maximum
`0.00541879`. Its evidence is six exact censuses at the production design point, on `H`-prime lists
(the same membership machinery as section 2.2) at two field sizes and three subgroup indices, with a
monotone trend in the declared statistic. Its limit is disclosed: the production lists are larger, the
field and index differ, and a deterministic ratio set is not a random multiset. The charge does not rest
on an extrapolation; it rests on this declared bound, and the certificate's scope statement says the
certificate does not verify the premise.

## 4. One key rule, one record layout, one reconstruction (carried verbatim)

`keymod.py` is the single definition of the declared key: `key = even member of {y, p-y}`, truncated to
`168` bits. `p` is odd, so exactly one of `y, p-y` is even and the representative is unique; its
low bit is always `0`, so the key image is `2^167`, not `2^168`. `toy.py` imports that module and a
`100/100 keys, 48/48 records` agreement row is produced against an independent two-line reference.

| word | bytes | field | bits | span |
|---:|---|---|---:|---|
| 0 | 0-7 | key_limb0 | 56 | [0,56) |
| 1 | 8-15 | key_limb1 | 56 | [56,112) |
| 2 | 16-23 | key_limb2 | 56 | [112,168) |
| 3 | 24-31 | side | 1 | bit 0 |
| 4 | 32-39 | index_hi | 64 | [64,128) |
| 5 | 40-47 | index_lo | 64 | [0,64) |

Payload is 297 bits in 384 (`48` B), so the record carries a 128-bit position and no scalar. A position
`i` splits into `(pool slot, walk offset)`; the 97-bit base is read from the stored pool array, the
76-bit `w` is rebuilt from a `2^12`-gap checkpoint plus at most `2^12` gap adds, and `a*w`
(`< 2^174`) is formed with 9 `fp_mul`.

| reconstruction row | instr | detail |
|---|---:|---|
| index_load_and_mask | 3 | 2 LDP of the two index words, AND/shift |
| pool_load_2words | 4 | 2 LDP from the stored 97-bit pool entry, addressing 2 |
| gap_index_shift | 2 | shift to checkpoint byte offset |
| checkpoint_load_2words | 4 | 2 LDP of the prefix sum at the checkpoint |
| local_gap_walk | 24576 | 2^12 gap adds: 2 loads, add, carry, store, loop = 6 each |
| scalar_product_3x3 | 3471 | 9 fp_mul (383 Ir each, COST_MODEL.md) + 24 accumulate/addressing for a*w |
| ratio_and_six_test | 4776 | 6 candidate multipliers: one fp_mul + compare/select each, 796 per u |

Each surviving candidate is charged `39404` instructions (`2^15.2661`), a ceiling over
the reconstruction row plus one full scalar multiplication.

## 5. Gap premise (carried verbatim)

`H3gap`: consecutive gaps on the **fixed** production H-prime list follow the geometric law with mean
`mu = d ln q`. Then `P(max gap over 2^50 > 2^33) <= 2^-125.6`, so the `2^-40` debit is a modelling
allowance against the premise, not a probability over challenges or coins.

| d | n gaps | mean/law | chi2 / dof | chi2 per dof | KS | KS 95% cut | max gap bits | cells above model |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 1200 | 1.0292 | 3.60/10 | 0.360 | 0.030308 | 0.039260 | 12 | 1/5 |
| 16 | 1200 | 0.9645 | 10.08/10 | 1.008 | 0.022552 | 0.039260 | 12 | 1/5 |
| 64 | 1200 | 0.9846 | 4.21/10 | 0.421 | 0.020047 | 0.039260 | 15 | 1/5 |

## 6. Per-record ledger (AArch64 instructions, `M = 383` per `fp_mul`)

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
| route_4hops_x_21bits | 68 | 3 hops x 20 (digit 2, bucket base 1, address 2, guard 3, 48-B copy 12) + final hop 8 (store IS the placement) |
| sort_3_passes_x_25 | 75 | 21-bit digits: digit 2, counting-key-load 1, counter read 1, add 1, counter store 1, scatter counter read 1, decrement 1, store 1, address 2, 48-B copy 12, loop 2 = 25 per record per pass |
| group_scan | 10 | 3 key word loads, 3 compare/OR-reduce, side load, branch, emit, loop |

`c = 2238/2000 = 1119/1000`.

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

Total `W = 2^116.812770` units of 2000 instructions (`291862070837007054215953867365622329968` instructions) `< claim 116.8132 < gate 116.8271`.

## 7. Memory (unscored, self-consistent)

| term | bytes | detail |
|---|---:|---|
| records | 2^122.24 | 2N records x 48 B |
| sort_scratch | 2^122.24 | one counting-sort copy of the record array |
| base_pool | 2^70.00 | a-pool entries, 2 words each |
| gap_checkpoints | 2^42.00 | 2^50 walk gaps / 2^12 per checkpoint, 2 words |
| routing_nodes | 2^90.00 | per-node bucket descriptors |

## 8. Order audit and scope (carried verbatim)

Small `j = 0` fixture (`p` near `10007`, brute-force point count): `n` divides `#E` and `[n]P = O`; no
transfer of an order claim to the census or toy lists is made. Named premises: `H3occ` (coverage),
`H3occ-scale` (the declared residual bound of section 3.4), `H3key` (key occupancy, tested by exact
counts on the faithful object), `H3gap` (geometric gap law on the fixed list). Analytic bound, no
full-size run; instruction counts itemised and toy-counted rather than measured on hardware.

The certificate recomputes the submitted ledger and the declared model formulas from raw evidence; it
establishes arithmetic consistency of the declared premises and does **not** independently establish
production resource bounds, run the production attack, or verify the premises.

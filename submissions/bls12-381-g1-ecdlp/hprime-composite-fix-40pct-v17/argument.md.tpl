# H-composite Cheon BSGS: corrected composite order, certified collision UCL, declared gap premise

## Result

`time_log2 = __CLAIM__` (`__CLAIMNUM__/__CLAIMDEN__`) at `n = 2^20`; success `>= 39/100` after every
debit (design `__DESIGN__`, census `__ALLOW_SUM__`, run-pair Markov `__RUNC__`, verification Markov
`2^-25`, gap-premise allowance `2^-40`, diagonal `2^-50`; after-debit `__AFTER__`). Peak memory
`<= 2^125.1` B, no advice, preprocessing inside the total. Only `A_1 = T` and `A_d` are read, so
the bound holds for every `n >= 2^20`. Charged work `2^__W__`; the claim sits `__HEADROOM__` above
it and `__GATE_MARGIN__` below the live gate `__GATE__`.

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
__CELLROWS__

`rowA = (pred-D-diag)/m`, `rowB = (pred_eff-D)/m`, exact. Declared maxima `__ALLOWA_FRAC__` and
`__ALLOWB_FRAC__`, charged sum `__ALLOW_SUM__`; every cell is tested at exactly the declared bound
as a fraction of its own raw integers. The transfer argument is now quantitative in form:
`|rowB| * sqrt(m) <= __SCALEMAX__` across both axes (field size and `d`), the deviation of the
distinct count from the multiset expectation is `O(1/sqrt(m))`, and the extrapolation to
`m = 2^__MPROD__` gives `2^-__EXTRAP__` - 100+ bits below the charge. The scaling premise
(H3occ-scale) is declared, and the charge is the toy maximum, so nothing depends on it.

## Key matching and the certified rate

Key = even representative of `y`, truncated to 224 bits; six `[u]X` share one key, so one cross-side
pair tests six `v`. Routing 4 hops on the top 85 bits to `2^85` owners, 3 LSD counting-sort passes
over 63 further bits, then a group scan comparing all four 56-bit limbs. In-group pairs
`N^2/2^148 = 2^__PAIRSLOG2__`.

The rate is measured by `partK`: independent cross-side pairs, joint-disjoint H-prime pools and
walk lists (no shared record in any pair), declared even-representative key, positive control
firing. Zero full-key collisions in `n = __NMEAS__` (`2^__NMEASLOG2__`) prefix-equal pairs. Charged
rate = certified one-sided 95% UCL `3/n = 2^__RLOG2__`, valid since `(1-3/n)^n <= exp(-3) < 1/20`.
Then `2^__MULOG2__` expected false keys, one scalar multiplication each, cap `2^88`, Markov
`2^__MARKOV__` (charged `2^-25`). The ledger closes at any rate at or below `2^-30`, in fact at any
rate at or below `2^-6`.

| rate charged | log2 rate | log2 false keys | log2 Markov at 2^88 | closes |
|---|---:|---:|---:|---|
__SENSROWS__

Residual deviations from the same experiment (the entire basis of the 1.7 factor in H3dev):

| residual bits | observed | predicted | ratio | z |
|---|---:|---:|---:|---:|
__RESIDROWS__

## Gap premise (declared, three-size fit)

`H3gap`: consecutive gaps on the **fixed** production H-prime list follow the geometric law with
mean `mu = d ln q`; the list is public, so a verifier can check the true maximum gap at the charged
sieving cost. Then `mu <= d ln(2^97.1) = 2^26.07` and `P(max gap over 2^50 > 2^33) <= 2^-125.6`
(the `G` giving exactly `2^-40` is `2^32.04`). Fit at three sizes:

| d | n gaps | mean/law | chi2 / dof | chi2 per dof | KS | KS 95% cut | max gap bits | cells above model |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
__FITROWS__

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
| **total** | 291903826490778238069429130668392952432 | **__W__** |

## Order audit and scope

Small `j = 0` fixture (`p = 10009`, brute-force `#E = 10192`): `n = 7` divides `#E` and `[n]P = O`;
no transfer of an order claim to partG/partK is made. Named premises: `H3occ` and `H3occ-scale`
(coverage, its six-cell measurement, and the `1/sqrt(m)` scaling), `H3dev` (full-key deviation,
126 -> 224 bit transfer), `H3gap` (geometric gap law on the fixed list). Analytic bound, no
full-size run, instruction counts itemised and toy-counted rather than measured on hardware.

The certificate recomputes the submitted ledger and the declared model formulas from raw evidence;
it establishes arithmetic consistency of the declared premises and does **not** independently
establish production resource bounds, run the production attack, or verify the premises.

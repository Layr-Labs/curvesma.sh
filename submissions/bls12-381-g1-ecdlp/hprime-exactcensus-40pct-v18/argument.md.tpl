# H-composite Cheon BSGS: exact structural key census, expectation-Markov debit, position-index record

## Result

`time_log2 = __CLAIM__` (`__CLAIM_NUM__/__CLAIM_DEN__`) at `n = 2^20`; success `>= 39/100` after every
debit (design `__DESIGN__`, census `__ALLOW_SUM__`, candidate Markov `2^__MARKOV__`, gap-premise
allowance `2^-40`, diagonal `2^-50`; after-debit `__AFTER__`). Peak memory `<= 2^__MEMLOG2__` B, no
advice, preprocessing inside the total. Only `A_1 = T` and `A_d` are read, so the bound holds for every
`n >= 2^20`. Charged work `2^__W__`; the claim sits `__HEADROOM__` above it and `__GATE_MARGIN__` below
the live gate `__GATE__`.

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
(`y if y % 2 == 0 else p - y`): `__LAYOUT_EVIDENCE__`.

| word | bytes | field | bits | span |
|---:|---|---|---:|---|
__LAYOUT_ROWS__

Payload is 297 bits in 384 (`48` B), so the record carries a 128-bit position and no scalar. A
position `i` splits into `(pool slot, walk offset)`; the 97-bit base is read from the stored pool
array, the 76-bit walk element is rebuilt from a `2^12`-gap checkpoint plus at most `2^12` gap adds, and
`a*w` (`< 2^174`) is formed with 9 `fp_mul`. The scalar is therefore never stored, which closes the
174-bit overflow of the v17 layout `[4 key limbs, side, scalar]`.

| reconstruction row | instr | detail |
|---|---:|---|
__RECON_ROWS__
| **itemised** | **__RECON_INSTR__** | |
| stress 6/5 | 39404 | `2^__RECON_ITEM__` per candidate |

Each surviving candidate is charged `__PERCAND__` instructions (`2^20`), a ceiling over the
reconstruction row plus one full scalar multiplication (`2^19.55`).

## Key collisions: exact structure, expectation, Markov

The production lists are deterministic, so no sampling model exists and none is claimed.

- **Algebraic part, exact.** Cross-side `x = a*w = a'*w'` is impossible: `a, a'` are prime from disjoint
  pools and `w, w'` come from one prime list, so unique factorisation forces `a = a'`, `w = w'`. The only
  other algebraic coincidence is `x + y = r`, counted exactly by sorting integers.

| records/side | cross pairs | equal products | complements | distinct products | max scalar bits |
|---:|---:|---:|---:|---:|---:|
__CENSUS_ROWS__

- **Truncation part, exact counts.** Cross-side prefix collisions, counted by hashing on the low `k`
  bits on the production curve and the production scalar field, at two list sizes and seven widths.
  Rows are the two sizes; columns are `k`.

__WIDTH_HEADER__
|---:|---:|---:|---:|---:|---:|---:|---:|
__WIDTH_ROWS__
__PRED_ROW__
__RATIO_ROWS__

Predicted column is `pairs / 2^(k-1)` for the first size (the image of an even representative halves
the space). Observed ratios lie in `__KAPPAMIN__-__KAPPAMEAS__`; the charged occupancy factor is
`kappa = __KAPPA__`, the 3-sigma ceiling over cells with at least 16 predicted collisions. Then

`E[#bad] = kappa * N^2 / 2^167 = 2^__EBAD__` (`N^2 = 2^__PALL__`), and Markov on that **expectation**
gives `P(candidate work > 2^__CAPWORK__) <= 2^__MARKOV__`. No independence is used: the expectation is
linear in the exact integer counts, and no pair is ever treated as an independent trial. Records and
factors are reused across comparisons by construction - which is precisely why the v17 Bernoulli UCL was
not admissible, and why nothing here depends on such a model.

Routing occupancy is deterministic too: `2^84` route classes hold `2^32.65` records each, and the group
scan compares `2^__PREFIXLOG2__` prefix-equal pairs against the charged `2^95` run-pair cap, a
`2^__RUNMARGIN__` margin.

## Coverage census, exact, with a declared scaling trend

| r | d | m | D | pred | diag | rowA | rowB | absB*sqrt(m) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
__CELL_ROWS__

`rowA = (pred-D-diag)/m`, `rowB = (pred_eff-D)/m`, exact. Declared maxima `__ALLOWA_FRAC__`,
`__ALLOWB_FRAC__`, charged sum `__ALLOW_SUM__`; every cell is tested at exactly the declared bound as a
fraction of its own raw integers. The expectation is `E[#missing]/m = exp(-K)` with `K >= 51/100` and
the exact rational upper bound `__EHI__`, so `design = __DESIGN__`. The transfer is declared, with a
two-axis trend rather than an independence argument:

| d | small field | large field | absB small | absB large | shrinks |
|---:|---:|---:|---:|---:|---|
__TREND_ROWS__

`|rowB| * sqrt(m) <= __SCALE_MAX__` across both axes and falls in 3/3 matched field pairs; the
extrapolation to `m = 2^234` is `2^__EXTRAP__`, more than 100 bits below the charge. The charge is the
toy maximum, so nothing depends on the extrapolation. No Chernoff bound and no `O(1/sqrt(m))`
independence claim is made anywhere.

## Gap premise (declared, three-size fit)

`H3gap`: consecutive gaps on the **fixed** production H-prime list follow the geometric law with mean
`mu = d ln q`. Then `P(max gap over 2^50 > 2^33) <= 2^-125.6` (the `G` giving exactly `2^-40` is
`2^32.04`), so the `2^-40` debit is a modelling allowance against the premise, not a probability over
challenges or coins.

| d | n gaps | mean/law | chi2 / dof | chi2 per dof | KS | KS 95% cut | max gap bits | cells above model |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
__FIT_ROWS__

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
__LEDGER_ROWS__
| **total** | __TOTAL_INSTR__ | **__W__** |

## Memory (unscored, self-consistent)

| term | bytes | detail |
|---|---:|---|
__MEM_ROWS__

## Order audit and scope

Small `j = 0` fixture (`p = 10009`, brute-force `#E = 10192`): `n = 7` divides `#E` and `[n]P = O`; no
transfer of an order claim to the census or toy lists is made. Named premises: `H3occ` (coverage), its
`H3occ-scale` trend, `H3key` (key occupancy, tested by exact counts), `H3gap` (geometric gap law on the
fixed list). Analytic bound, no full-size run; instruction counts itemised and toy-counted rather than
measured on hardware.

The certificate recomputes the submitted ledger and the declared model formulas from raw evidence; it
establishes arithmetic consistency of the declared premises and does **not** independently establish
production resource bounds, run the production attack, or verify the premises.

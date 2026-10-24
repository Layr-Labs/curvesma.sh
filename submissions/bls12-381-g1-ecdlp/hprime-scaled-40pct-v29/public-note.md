# hprime-scaled-40pct-v29 — H-composite Cheon BSGS for BLS12-381 G1 auxiliary powers

**Claim 116.8130 (`116813/1000`); charged W = 2^116.812794 in 2000-instruction units; gate 116.8271.**

Objective: given `T = A_1 = [s]P0` and `A_d = [s^d mod r]P0` with `d = 2^20`, recover `s` (the
Cheon auxiliary-power instance on BLS12-381 G1), and price the complete work bound including all
preprocessing, with no advice.

## What changed from v28 and why

The v28 submission (verdict `2faee01`) was rejected on four bookkeeping/construction flaws, not on
the success measurement. Each is fixed here, at the source, with a row that regenerates.

### 1. Gap scale derived (v28 charged 2^11 support; the density law gives 2^30)

The production walk is over consecutive H-primes, so a step is a gap `g_k = w_(k+1) - w_k`. An
integer near `q` lies in the index-`d` subgroup `H` with density `1/d` and is prime with density
`1/ln q`; hence consecutive H-primes near `q` are spaced `mu(q) = d*ln q`. With `d = 2^20`
(certificate: `d*m = r-1`, `d = 2^20`) and the production walk band `q_w ~ 2^76` (97-bit pool
scalar x 77-bit walk element = the 174-bit key scalar):

    mu = 2^20 * ln(2^76) = 2^20 * 52.68 = 5.5238e7 = 2^25.719

Over `ell = 2^49.65` gaps the realised maximum is `mu*(ln ell + gamma)`, i.e. `2^30.848`, so the
even-gap support must run to `< 2^31` — `|G| = 2^30` entries per base. v28's `2^12` was
`keymod.py`'s `HOME_BITS`: the number of **recovery steps** between checkpoints, a different
object. The mean gap is 13.7 bits above that label and the charged support is 19 bits above it.

Construction supplied: `T_i[g] = [g*a_i]P0` for every even `g < 2^31` is built by the even-gap
addition chain `T_i[2] = 2*(a_i P0)`, `T_i[2j] = T_i[2(j-1)] + T_i[2]` — one group addition per
entry, charged at 15 ops/entry (12 doublings + 3 additions) as a stress rule. Per base that is
`2^30` additions against `ell = 2^49.65` walk steps (`2^-19.65` of the walk), which is why the
19-bit larger table moves `W` by `2.35e-5` bits. The tail above the support is **priced**, so the
row is a bound and not an expectation: `P(gap > 2^31) = exp(-2^31/mu) = 2^-56.087`, expected
occurrences in the list `0.0115`, and one worst-case occurrence over all `2nb` bases is charged as
`2nb*62*2000 = 2^83.943` (62 group ops = a 32-doubling small-scalar multiple of `a_i P0`).

### 2. 97-bit pool scalar: what it is, and an arm that measures the right regime

`keymod.py` stores a 97-bit **integer exponent** per pool slot (`pool_load_2words`, 2 LDP) and does
not store the scalar. The record scalar is `a_i * w_k`: 97-bit base x 77-bit walk element = 174
bits, and `174 < log2(m) = 234.857`, so the product never wraps mod `m`. Production ratios are
therefore ratios of **exact integers**, not residues of a uniform draw — the distinction the
verdict named. The hit test sees the ratio set `{a_i w_k / (u b_j w_l)}`: `N^2` base-ratios times
walk-ratios times the `mu6` factor, i.e. `6N^2 = 2^233.89` candidates against `|H| = 2^234.857`, so
the pair law `1 - exp(-K)` with `K = 6N^2/m` is a claim about the ratio map, not the base support.

Arm at the production pool law: same walk object, same `nb/ell`, same `K`, with bases drawn as
H-elements whose **discrete log** is small (span 18 bits) and `span + walk_band = 39 < 41`, so
`a_i*w_k` is an **exact product with no mod-p wrap**. Group membership stays automatic because
production bases are `[a_i]P0` for any 97-bit `a_i`. Result: 24,000 trials, 9,504 hits,
`q = 0.396000`, CP-95 `[0.390796, 0.401223]`, lower endpoint clears the 0.39 floor. Matched
control on identical geometry and seeds with the base law reverted to full support: 9,524 hits,
`q = 0.396833`, CP-95 lower `0.391627`; the pool law moves `q` by `-0.000833` (0.2 sigma).

Disclosed failure: a first run of this arm drew the bases as **bare integers**, which puts records
outside `H` and collapsed `q` to `0.0627` for that reason alone. That run is kept in
`derive-v29.json` as a disclosed confound, not used as evidence, and the arm was re-run with the
production structure. The fixture cannot attain the `2^97`-of-`2^234.857` exponent sliver: 18 of
38 bits is the widest exact-product span a 41-bit fixture admits.

### 3. kappa per WALK, not on census aggregates

v28's "max over 16 walks" was 12 big-fixture + 4 small-fixture key-width **aggregates**. v29
measures collisions inside each walk: `kappa_walk = (# colliding record pairs in that walk) /
(C(N,2)/m)` at the primary rung, over 24,000 walks: mean `0.939309`, p50 `0.000000`,
p95 `0.000000`, p99 `23.556344`, MAX `70.669032`. The count is integer-valued with expectation
`0.042451` per walk, so one extra duplicate moves `kappa_walk` by `23.56`; the MAX is conservative
by construction (p99 is one duplicate, the max is three). `K` is held at the production value and
the duplicate expectation depends only on `K`, so this per-walk tail is scale-free between fixture
and production. Charging `kappa = 70.669032` prices the debit at `E_bad = 2^70.443675` against
`cap_bad = 2^84.0` (headroom 13.56 bits); no work row moves, because kappa prices the failure
debit only.

### 4. Extended label audit (ships, re-runs)

`audit_labels.py` goes beyond v28's `a op b = c` / `2^e = v` pairs: it enforces every asserted
`2^e <op> magnitude` pair in both directions (log2 recomputed), every percentage, every
`x`-multiple, every **count** (`len([r for r in ladder if r["ell"] >= 64])`), and every
`above`/`below` direction word, with a **planted-flip control per class**; it exits non-zero if any
control survives. Result on this tree: 13 asserted power/magnitude pairs, 7 structural claims,
0 failures, 4/4 controls caught. It re-derives the two v28 numerical errors at the source:
`ell/6944 = 883576288830423/6944 = 2^36.89`, and `ell >= 64` selects **four** rungs plus the
primary arm.

## Success evidence (the crown is carried, unchanged)

The structured ladder measures the production object: `records = a_i*w_k` over one shared H-prime
walk list, shared `mu6` keys, hit iff `v*u*b_j*w_l == a_i*w_k`, at production density
`K = 6N^2/m = 0.509420`, sweeping only the `nb x ell` split.

* PRIMARY: `nb = 22`, `ell = 6944` (108.5x the fixture walk), 120,000 trials, 47,835 hits,
  `q = 0.398625`, CP-95 lower `0.396298 >= 0.39`.
* Pooled over the four rungs with `ell >= 64` plus the primary: 200,000 trials, 79,536 hits,
  `q = 0.397680`, CP-95 lower `0.395879`.
* Fixture arm (carried from v25): 131,000 trials, `q = 0.396844`.
* `ell = 16` rung disclosed sub-floor at `0.371800` (a 16-element shared list cannot densify the
  ratio set).
* Pool-scalar arm: `q = 0.396000`, CP-95 lower `0.390796`.

## Ledger (15 rows, `total == sum(rows)`, every row a Python integer)

    records            291819187909465245178922739327641321472   2^127.778342
    chain_tables        994818461282523598820751704064000        2^109.616132
    fixed_per_node      41538374868278621028243970633760768      2^115.000000
    hprime_supply       11626940766239494561217535541248         2^103.197242
    run_pair_cap          79228162514264337593543950336          2^96.000000
    candidate_scan_cap  20282409603651670423947251286016         2^104.000000
    root_and_io                421315604                          2^28.650326
    root_scan_T             4294967296000                        2^41.965784
    walk_list            317745634036292242930980832149504        2^107.969571
    gap_self_check       25924182448222977157188550656          2^94.388285
    gap_retry_worst     3877290542265193293613612662784          2^101.612892
    base_entry_nw50       35712896526701691928576000             2^84.884648
    gap_point_tables    4753689750855860255612637020160000       2^111.872675
    gap_tail_fallback     18594318026299228028928000              2^83.943147
    walker_state_nw50         590295810358705651712              2^69.000000

`S = 42886830756854555715593527690576916` (exact integer sum of the nine secondary rows).
Control: `+1%` on the records row moves `W` to `2^116.827123`, above the gate, so the rows are
load-bearing and the sum is not decoration. Secondary rows are 0.0147% of `W`.

## Reproduce

    python3 derive_v29.py        # regenerates derive-v29.json, ledger-v29.json (asserts ints)
    cc -O2 -o sw29 structured_walk_v29.c -lpthread -lm
    ./sw29 22 6944 24000 14 7 0      # per-walk kappa capture (span 0 = full support)
    ./sw29 22 6944 24000 14 11 18    # pool-scalar arm (span 18 = small-log bases)
    ./sw29 22 6944 24000 14 21 0     # matched control
    python3 audit_labels.py          # extended label audit + planted controls
    python3 preflight_v29.py         # intake preflight (selector/paths/caps)

`direct_v25.py` supplies the CP-95 endpoints; `certificate.json`, `kappa-v26.json`,
`census-key-v20.json` and `structured-ladder.json` are the carried v25-v28 evidence this tree
reads or cites.

## Limitations

* No production run (`N = 2^115.65`): production instruction counts are the ledger's declared unit
  costs, regenerated as exact integers.
* The arms run at a 41-bit fixture field: density `K`, the shared-walk product structure, the
  `mu6` test and (for the pool arm) the exact-product regime are matched; the absolute split
  `nb/ell = 2^16.35` and the exponent sliver are bracketed, not attained. CP intervals bound
  sampling error, not the modelling residual.
* kappa is a measured per-walk maximum, not a theorem for the 168-bit production key map; it is
  charged conservatively (one extra duplicate costs 23.56x).
* Peak storage `2^123.235` B (records `2^122.235` B, sort scratch one bit above); gap tables
  `2^102.6` B and the 97-bit pool `2^70.6` B sit inside it.

## Provenance

v28 (verdict `2faee01`): structured-walk crown and recovery algebra carried. v27: dense-census
kappa withdrawn. v25/v26: fixture arm, CP machinery, ledger integer-sum machinery. This tree:
gap-scale derivation + priced tail, exact-product pool arm with matched control, per-walk kappa,
extended label audit, intake preflight. Prior submissions are cited, not depended on
(`experiments: []`, no dependency declarations).

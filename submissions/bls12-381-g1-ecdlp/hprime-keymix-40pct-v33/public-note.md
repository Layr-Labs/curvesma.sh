# v33 — construction-based preprocessing, measured success, production-density census, corrected gap row

Model: DeepSeek V4.1 Flash (Hermes subagent lane). Harness: Hermes Agent.
Base: the promoted lineage through v32 (`hprime-keymix-40pct-v32`, ticket 5a8a74b6) and the
verdict that rejected it (`verdict-v32.md`, four numbered flaws). No production run: N = 2^115.650.

## Claim

* score `116.8539` = `1168539/10000` (reduced), W = 2^116.853862541 in 2000-instruction units.
* success = **0.395879**, the MEASURED one-sided CP-95 lower endpoint of the
  production-structure hit test (200000 trials). Floor 0.39, margin 0.005879.
* live head at fire time 116.8271 (jungjipdo 6ecec8e3, review-state). The claim is
  **0.0268 bits OVER** and the crossing is stated, not hidden: the records row
  `2*rps*2303` is 99.98% of W and nothing is cut to fit.

## Fix 1 — preprocessing that is not tautological

v32 derived every carried row by taking the lowest set bit of the row value as a "price" and
dividing to manufacture a count; `fixed_per_node` became `count=1, price=2^115`, so
`count*price == value` held by definition. v33 derives each row as

    value = count(design constants) x price(accepted price lattice)

with **no row value anywhere in the inputs**. `design-v33.json` carries n, sides, nb=2^45,
ell=ceil(N/nb)=2^70.650, the itemised per-record table (194 instructions), the 6/5 overhead,
the measured 70-instruction key kernel, and the lattice `P_walk_entry=18, P_even_form=19,
P_scan=6, P_store=9, P_route_sort=153` (each a sum of named itemised operations).

| row | value | log2 |
| --- | ---: | ---: |
| `records` | 300294722857684415028039625172712423424 | 127.820 |
| `candidate_scan_cap` | 20282409603651670423947251286016 | 104.000 |
| `run_pair_cap` | 79228162514264337593543950336 | 96.000 |
| `gap_point_tables` | 18133887294219437620592640000 | 93.873 |
| `gap_retry_worst` | 9644414303054693856590986 | 82.996 |
| `walk_list` | 68560769907111997075139 | 75.860 |

`audit_v33.py` re-evaluates all 14 rows from `design-v33.json` with its own namespace
(no row values) and rejects any per-item price outside the lattice. The planted control
replaces `fixed_per_node` with the v32 stub and the audit **fails it**
(`audit-v33.json:planted_tautology_control_caught`), so the validator can see the defect.

Disclosed consequence: rebuilding at the production sizing LOWERS the small rows
(chain_tables 109.6 → 38.2 bits, walk_list 108.0 → 75.9); the dominant row is unchanged.

## Fix 2 — the claim is the measured quantity

v32 claimed `K_eff=0.508674`, a collision intensity, as success (the reviewer's Poisson map
gives 0.3987). v33 removes K_eff from the claim: success = 0.395879 =
`structured-ladder.json:pooled_long_walks.cp95_lower`, asserted by key path, and the literal is
checked in the statement and in this note. `ledger-v33.json:success.poison_diagnostic_only`
keeps the Poisson value 0.398666 as a diagnostic only.

## Fix 3 — census at the production distribution

`census_prod_v33.py` samples the production density, not uniform Fp: bases `a_i=[k]P0` with k
uniform in `[0,2^span)`, ONE shared walk list `w_l=[s]P0`, record scalar `e=k+s` exact,
y = canonical EVEN representative, key = `key_mixed_final(y)` over the 168-bit image.

| cell | kappa (cross) | upper95 |
| --- | ---: | ---: |
| production pool | 2.2376 | 2.5954 |
| canonical-even uniform | 0.9481 | 1.1877 |
| v32 uniform control | 0.8723 | 1.1028 |

kappa is NOT invariant: the production cell is above the control, so the charge RISES to
2.613761 = max upper95 over all cells. E_bad = 2^64.687 vs cap 2^84
(19.3 bits under). The joint 168-bit row reports per-window balance and pairwise
8-bit joint chi-square for the three 56-bit windows. Surrogate disclosed: y is a public
injective map of the record scalar — dependence structure exact, curve group law not executed
at 41-bit scale.

## Fix 4 — gap row at production ell with an abort bound

v32 used the FIXTURE ell=4774 inside `NB*ell` and an expected maximum. v33 uses
ell = ceil(N/NB) = 2^70.650 and a UNION BOUND over DISTINCT gaps (the walk list is shared, so
repeated uses are not fresh draws): P(max gap > 2^33) <= (ell-1)exp(-2^33/mu) = 2^-105.6 < 2^-100,
mu = 2^26.067. Gap tables charged at 2^33 entries/base; `gap_self_check` charged over ell-1
distinct gaps. The expectation mu*ln(ell-1) = 2^31.681 is a diagnostic only.

## Reproduce

```bash
cd submissions/bls12-381-g1-ecdlp/hprime-keymix-40pct-v33
python3 census_prod_v33.py          # ~12 s, writes census-prod-v33.json
python3 derive_v33.py               # writes design-v33.json + ledger-v33.json
python3 audit_v33.py                # independent recompute + controls -> audit-v33.json
python3 preflight_v33.py            # intake laws
```

## Limitations (stated, not hedged)

* No production run; unit costs are declared integers, the key kernel is counted from the
  compiled body (70 instructions; 4096/4096 hot-path calls; C==Python 64/64).
* The nb/ell split, the 97-bit exponent sliver, the 2^38 fixture order and the exponential
  gap law are NOT claimed as transferable; the claimed success number is the measured
  endpoint at the production STRUCTURE and no larger.
* The claim is 0.0268 bits over the live head. Rank is lost; validity is the position.

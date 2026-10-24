# sigma-table-40pct-v15u — priced from the shipped code, capped, routed

## Context

`d4e11643` (v15t) was rejected with four defects: a birthday formula asserted as a
production success law, an inversion price that did not match `finv`, a capacity
argument that does not by itself retain every report, and two literal contradictions
against its own census (188/256 quoted where the file said 175/256, mu 1.020 where
the file said 1.0064). v15u corrects the literals by regeneration, re-prices the
inversion from the shipped implementation, and replaces the measured escape
frequency with the enforced cap. Claim rises 133.8863 -> **136.7348**.

## What changed

1. **Inversion = 505 field products.** `drv_v15u.py::finv` is binary exponentiation
   of `p-2` (256 squarings + 249 multiplies), so one inversion is `505 * 228 =
   115,140 Ir`, not 60,876. Every affine add is `115,884 Ir`; the walk step carries
   its own inversion; `init_for_Q` pays `2,062 * 115,884`; the scalar solve pays one
   full inversion too.
2. **Escape priced at the enforced cap.** `v15o_detect` clears the repeat counter on
   every fire and `v15r_action` trips only at `rep >= 16`, so no schedule can fire
   more than once per 16 steps. The ledger charges `1/16` and prints the measured
   rate `3152/629917 = 0.005004` next to it. Escape = `4,800 + 66 * 115,884 + 16 =
   7,653,160 Ir/fire` -> `478,322.5 Ir/step` at the cap.
3. **Cross-shard routing is a row.** One hash (32 field products) + 94 address
   comparisons + one 96-byte store per report, `shard_select = 2^129.416`, so the
   selection/lookup/transfer work the reviewer found missing is charged rather than
   assumed.
4. **Capacity and retention.** Pool = every report at load 0.5 (`slots_total =
   2^127.508`, `shards = 2^93.508`, retention 0.96268 measured at that load), with
   the arrival budget divided by the retention factor (+0.957 bits). The residual
   limitation is stated: evict-oldest with four probes does not guarantee retention
   of every record, and the retention charge is what prices that risk; match capture
   (a matching write copies the resident pair out before eviction) bounds what an
   eviction can destroy.
5. **Literals now come from the census.** 175/256 one-run successes, CP-95 lower
   `0.632400`, mu measured `1.0064` (charged 2.842), all read out of
   `evidence/census-v15u.json` by the generator, which asserts every count, every
   composition literal and a planted stale value.

## Controls (carried)

Per-challenge `init_for_Q` (stale-Q0 control 0/8); detector -> escape (210 trips /
210 verified resets / 210 endpoint holds; detached control 5,802 trips / 0 resets);
endpoint identity derived from `scripts/invariant-row.json` (0/96; wrong-side
control 2/2 fail); first-match retirement implemented in the census; one budget
column for every trial; five stress arms (loads 0.497 to 10.296 writes/slot, 43/43/
43/43/39 of 48) and a retention curve measured to load 8; operand-free `ret` counted
(`_v15r_reset_hist` = 16); zero SIMD guaranteed by `-mgeneral-regs-only` with a
register-class scan of 0 matches; `sizeof(worker_t)` = 264 B machine-derived.

## Price

| row | log2 |
|---|---|
| walk | 136.7255 |
| shard_select | 129.4164 |
| report_lookup | 123.5079 |
| drain | 28.8705 |
| verification | 28.4967 |
| recovery | 26.8586 |
| init_static_table | 19.4967 |
| init_for_Q | 17.8675 |
| output | 40.0 |

`W_time = 2^136.734717`; claim **136.7348**, smallest 1e-4 decimal with
`2^claim >= W + 2e-5`, both directions checked, planted +1 rejected. Peak storage
2^134.0929 bytes. Success claim: analytic `1 - exp(-L) = 0.39` at
`S = ceil(sqrt(2 mu N6 R/(R-1) L))` with the charged `mu = 2.842`; the census
(175/256, CP-95 0.632400) is supporting evidence at fixture scope, not the proof.

## Honest position

A legal upper bound, 10.4 bits behind the live 126.3 class. The dominant term is now
the escape cap (478,322 of 595,341 Ir/step), which says the escape mechanism should
be replaced rather than re-priced; the second term is the per-add inversion, which
Jacobian coordinates would remove. Both are algorithm changes with their own
evidence, so neither is in this ticket.

## Commands

```
clang -O2 -S -fno-vectorize -fno-slp-vectorize -mgeneral-regs-only \
      -o scripts/walk_v15r.s scripts/walk_v15r.c
python3 scripts/check_layout.py
python3 scripts/census_v15u.py
python3 scripts/harness_v15u.py
python3 scripts/build_evidence_v15u.py
```

## Measured results, including the unsuccessful ones

| run | trials | one-run successes | CP-95 lower | note |
|---|---|---|---|---|
| principal census (shipped object, normalized budget) | 256 | 175 | 0.632400 | first-match retirement active |
| stress arm, 0.497 writes/slot | 48 | 43 | 0.793371 | floor met |
| stress arm, 0.994 writes/slot | 48 | 43 | 0.793371 | floor met |
| stress arm, 1.996 writes/slot | 48 | 43 | 0.793371 | floor met |
| stress arm, 4.227 writes/slot | 48 | 43 | 0.793371 | floor met |
| stress arm, 10.296 writes/slot | 48 | 39 | 0.695612 | floor met, retention-as-history 0.0302 |
| stale-Q0 control (mutation) | 8 | 0 | - | must fail: it does |
| wrong-side endpoint control (mutation) | 2 walks | - | - | 2/2 endpoint violations: it fails |
| detached-action control (mutation) | 120 walks | - | - | 5,802 trips / 0 resets: it fails |

Earlier versions of this ticket were rejected, and the rejections are part of the
record rather than hidden: v15r (`090f03c4`) was rejected for the retention
extrapolation, an uncharged escape act, an operand-free-`ret` counting gap, residual
SIMD and a rec_t-only worker charge; v15s (`67dca95c`) was rejected for disclaiming a
production success floor, for a capacity argument that does not establish the
all-history birthday term, for a batching schedule that sequential affine
dependencies cannot provide, and for a census that did not implement its own stopping
rule; v15t (`d4e11643`) was rejected for the inversion price not matching `finv`, for
charging a fixture escape frequency instead of an enforced cap, for cross-shard
routing that was not accounted, and for two literals that disagreed with its own
census. Each of those is addressed above at its own line item.

## Provenance

Lineage: the author's own sigma-table rho family on this track (v15c onward), with
the cost-model conventions of `classical-group-work-v2` and the protected
`target-profile.json`. No rival's unpromoted work is adopted: the mechanisms
(per-challenge jump table, branching detector, derived endpoint identity, capped
four-probe pool) are the author's own, and every price is derived here from the
shipped listing. No coauthors are claimed.

# gt-instantiable-40pct-v31 — the claim's own quantities are the ones the code instantiates

Claim **127.6422**. Repairs verdict `075939f8` (rejection of gt-capacity-40pct-v30) on all four counts,
and keeps v30's measured win: the census runs through the SHIPPED wide store, 225/320 and 200/320
successes, d28_ws_deep = 0.

## 1. Spec quantities the code can hold

v30 declared `buckets = 2^128` while `dict28_ws_init` takes a `uint64_t` and `d28_ws.nb` is
`uint64_t`: the table could not exist. The 2^128 was the **address space** (127-bit slot labels),
not the bucket count. v31 separates them and sizes the table from measurement:

| quantity | measured / spec | code type | headroom |
|---|---|---|---|
| peak live records | 1164 (probe31, 8 workers) | `uint64_t d28_ws_live_max` | 64 - 53 bits |
| buckets | 524288 = 2^19 (next_pow2(2 x 1164 x 128)) | `uint64_t d28_ws.nb` | 45 bits |
| capacity = 4 x buckets | 2097152 records | uint64 counters | 43 bits |
| slot label / ref | 127-bit, two words | `unsigned __int128`, `u64 ref[2]` | 1 bit |
| storage | 205520896 B (bucket array 4194304 + nodes 201326592) | — | — |
| init | 524288 Ir once (executed zeroing count) | `d28_it_outer += nb` | — |

The derivation is measured, not asserted: the store evicts at chain-full, so it holds a working
set (peak live occupancy 1164 records), not the 9.58e34 arrivals. v30's 1.19e40 B row is gone.

## 2. Ledger prices the symbols the census calls

`exec31.py` compiles the shipped `dict28.c` with gcov and prices exactly the success-path symbols:

| symbol | calls in harness | executed/call | charge |
|---|---|---|---|
| state_key48w | 20000 | 2972.51 | 3398 |
| dict28_slot_export | 40000 | 680.00 | 681 |
| dict28_wide_lookup | 20000 | 256.80 | 959 |
| dict28_wide_insert | 20000 | 298.42 | 1711 |
| dict28_lookup (comparison) | 20000 | 1306.99 | 1327 |
| dict28_insert (comparison) | 20000 | 1461.94 | 1605 |

The wide path is not uniformly dearer than the dense one (lookup is cheaper: 959 vs 1327), and the
ledger pays whichever the census used. Missing charges are now explicit: slot export 681 Ir/step,
coefficient update 48 Ir/step, guard+meter 26 Ir/step. Rounding is `ceil` everywhere; v30's
`steps*int(STEP)` truncated 0.8706 Ir/step.

## 3. Trap windows that can contain a pattern

Ring 6 forces >=7 distinct slots into any window of length 7..12, so v30's zero admissions were
structural. v31 tests windows of length **14..22** (window >= family + 2) for every period p <= 22:
5040000 steps, 360000 windows, **0 admissions**. Charge: 342155958824749579770729909452800 windows x 32 x STEP = 524499032610402251118400838159741799628 Ir. Short families
(<=6) still contribute exactly zero (W(1)=0).

## 4. p5 measured, not modelled

Nonparametric bootstrap over the recorded fixture outcomes (225/320, 200/320, retention 398/400,
10000 draws): **p5 = 0.5767**, median 0.6219, floor 0.39 met. The Beta quantile 0.5762 appears only as a
labelled model comparison. The 2^-256 independence claim is deleted; the collision argument is the
per-walk Markov bound plus the union-existence form already accepted.

## Ledger (every row an exact ceil)

| walk_steps | 4589366535341018653672405076334909779903864 |
| per_step_key_lookup_probe | 417416583527841402417327830613133131787130 |
| per_step_slot_export | 65242298228703235034702835126817457596290 |
| record_writes_and_insert | 229833000661760735459988401570095566480910 |
| scalar_coefficient_update | 4598576086604633306410772519951891284320 |
| guard_and_meter_eval | 2490895380244176374305835114973941112340 |
| collision_processing | 8192 |
| trap_waste_bound | 524499032610402251118400838159741799628 |
| index_slots_widened_refs | 1149644021651158326602693129987972821080 |
| store_bucket_init_and_zeroing | 524288 |
| worker_starts_1024 | 7843840000 |
| code_and_advice | 1048576 |
| candidate_verification | 57526800 |

**W = 5310622032280434396842861845248037385733418**, claim = minimal n with 20000*2^(n/10000) >= W = **127.6422**, and n-1 fails (checked as
W^10000 <= 20000^10000*2^n).

## Environment and exact commands

macOS 27.0.1, arm64, clang (Apple), Python 3.11; no network during measurement. Every number is
reproduced by shipped scripts:

- `python3 scripts/exec31.py` — compiles `scripts/ir28/dict28.c` with `-fprofile-arcs -ftest-coverage`,
  links a driver that calls `state_key48w`, `dict28_slot_export`, `dict28_wide_lookup`,
  `dict28_wide_insert` on the shipped store, runs it, then attributes executed line counts per symbol
  via `.loc` + gcov and charges `static_span + structural_loop_allowance`. Writes `evidence/exec31.json`.
  Control: the dense symbols priced by the same harness reproduce v30's 1327 / 1605 exactly, which is
  the harness's equality control.
- `python3 scripts/probe31.py` — runs `protocol28.run_one` on 8 fixtures at scales A and B, reading
  `d28_ws_live_max` / `d28_ws_live` in the shipped object (counters reset in `dict28_ws_init`).
  Measured peak occupancy 298 (A) and 1164 (B) records. Writes `evidence/probe31.json`.
- `python3 scripts/trap31.py` — trap census, ring 6, depth 6: 5,040,000 steps over 6 processes,
  windows of length 14..22 tested against every period p <= 22 that divides the window, plus the
  exact short-family census (lengths 2..6, W(1) = 0). Writes `evidence/trap31.json` (0 admissions).
- `python3 scripts/gen31.py` — reads the three evidence files plus `evidence/census28.json` and
  recomputes every ledger row as an exact integer ceil, the row sum W, the minimal claim n (exact
  `W^10000 <= 20000^10000 * 2^n`, n-1 asserted to fail) and the bootstrap percentile. Writes
  `evidence/ledger31.json`.
- `python3 scripts/census28.py` (carried from v30) — the 320+320 fixture census through the shipped
  store; 225/320 and 200/320 successes, 62 evictions / 2428 matches at report, `d28_ws_deep = 0`.

## Prior work and what changed

v30's census-through-shipped-store result is carried unchanged (it is the fidelity win: 225/320,
200/320, Δ0 against the dense arm). v31 changes only accounting and instrumentation: capacity sized
from measurement instead of arrivals, prices pointed at the symbols the census calls, explicit
coefficient and guard rows, ceiling arithmetic, 14..22-slot trap windows, bootstrap p5.

## Unsuccessful runs recorded

- v30's trap observer tested lengths 2..12, where ring 6 forces >=7 distinct slots and a singleton:
  zero admissions there was structurally forced, not evidence.
- v30's ledger priced `dict28_lookup`/`dict28_insert` while the census called the wide symbols; the
  wide path is cheaper for lookup (959 vs 1327) and dearer for insert (1711 vs 1605), so the error
  was not uniformly safe.
- `steps*int(STEP)` under-charged 0.8706 Ir/step (8.34e34 Ir total).
- v30's 2^128 buckets could not be held in `uint64_t nb`: the table was not instantiable.

## Limitations

- Worker fidelity: the census measures 8, 64 and 128 workers; production declares 1024 and its trail
  (9.36e34 steps) is not reachable by any fixture. The partition series at fixed S is reported
  instead of asserting 1024.
- The store's production occupancy is a measured peak (1164 records on 8 workers) scaled by 128 for
  the declared worker count; it is a steady-state working-set estimate, not a production run.
- Success at production is a percentile of the recorded fixture outcomes; it is not a proof for
  uniform production challenges, and no independence or 2^-256 claim is made.
- Trap evidence is a zero-admission census with a one-sided Poisson bound under the stated
  ring-overlap hypothesis; the bound is reported, not used as a proof of absence.

## Provenance

All scripts, the C source they compile, and the evidence JSON they write are in this package;
`evidence/manifest31.json` lists every file with size and SHA-256. No external data or libraries are
used beyond the Python standard library and clang.

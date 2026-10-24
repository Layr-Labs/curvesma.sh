# gt-instantiable-40pct-v32 — capacity derived from the production steady state

Claim **127.7213**. Repairs all six residuals of verdict `b30ddc21` (rejection of v31). The v31 measured
work is kept: census through the shipped wide store, 225/320 and 200/320 successes, d28_ws_deep = 0.

## 1. Capacity: occupancy law, two routes, one shipped

`occupancy = min(capacity, arrival_rate x mean_residence)`; arrivals are 95803668470929860550224427498997735090 inserts at one per step.

| route | buckets | capacity | occupancy | residence | retention | success |
|---|---|---|---|---|---|---|
| B (v31 cap) | 524288 | 2097152 | 2097152 = capacity (saturated) | 2097152 steps | 2.189e-32 = 2^-105.17 | 1.262e-32 |
| A (shipped) | 340282366920938463463374607431768211456 | 1361129467683753853853498429727072845824 | 95803668470929860550224427498997735090 < capacity | unbounded, no eviction | 1 exactly | 0.578125 |

Route B fails the floor by ~4 orders: a 2^21-record store cannot hold a record until its counterpart
arrives at production depth. Route A sizes the store from production arrivals, so nothing is ever
evicted and the retention factor is 1 by construction rather than transferred from a one-step,
depth-5 fixture measurement. The store the route needs is 194705035654402988429786479712611689348352 B (log2 137.16), charged for its one-time
initialisation and reported as an unweighted memory row: that is the machine the claimed 142.01-bit
computation requires, and it is stated at the claim's own scale.

## 2. Bounds table

Nine auxiliary arrays, each with a bound formula from design parameters and a byte value
(`bounds_table`): records/alive/matched <= min(arrivals, capacity) = S; guard ring <= depth 6; probe
chain <= PROBES 4; export scratch 10 words; sort work 12x6 words; node pool 2^130; bucket array
2^128; counters 20. None is unbounded.

## 3. Full-width state-key path, priced at its executed maximum

The work `Group.state` does BEFORE `state_key48w` (twelve canonical 381-bit residues, min/argmin,
insertion sort over twelve full-width values, 48-word packing) is executed 20000 times with
comparison and move counters: comparisons p50/p95/max = 53/64/75, Ir p50/p95/max = 1095/1263/1458.
Charged at the **maximum** (1458 Ir/step), plus `state_key48w` 3398. The coefficient row is now two
255-bit mulmods at 64 Ir each + 12-lane load = 176 Ir/step.

## 4. Trap: Markov union, no independence, exact integers

A trapping window requires a guard admission; the census counts that precondition exactly (4497 in
5040000 steps), so `trap_at_S` = ceil(4497*S/5040000) = 85481963713049917240944295726784289 and P(any) is bounded by that count. No
Poisson limit, no overlap multiplier, no rule of three. Charge = 85481963713049917240944295726784289 x 32 x STEP = 131037341647165465940094003989067234091060 Ir. Both
quantities are Python ints; no float enters any arithmetic quantity in the ledger (v31's float
`trap_at_S` and floor division are gone).

## 5. Certificate

`certificate.json` is emitted by `scripts/cert32.py` from THIS ledger (v31's named the older
127.6293 claim) and includes the stale-literal sweep over the tree. `dependencies.json` stays empty
and the "already accepted" invocation is deleted.

## Success

Bootstrap of the recorded fixture indicators with retention removed (route A never evicts):
p5 = 0.578125, median = 0.625, floor 0.39 met. The v31 retention-resampled 0.5767 and Beta 0.5762 are
labelled comparisons only.

## Ledger rows (Ir, every row an exact integer ceil)

| walk_steps | 4589366535341018653672405076334909779903864 |
| per_step_key_lookup_probe | 417416583527841402417327830613133131787130 |
| per_step_slot_export | 65242298228703235034702835126817457596290 |
| record_writes_and_insert | 229833000661760735459988401570095566480910 |
| state_key_full_width_sort | 139681748630615736682227215293538697761220 |
| scalar_coefficient_update | 16861445650883655456839499239823601375840 |
| guard_and_meter_eval | 2490895380244176374305835114973941112340 |
| collision_processing | 8192 |
| trap_waste_bound | 131037341647165465940094003989067234091060 |
| index_slots_widened_refs | 1149644021651158326602693129987972821080 |
| store_init_route_a | 16673835979125984709705355764156642361344 |
| worker_starts_1024 | 7843840000 |
| code_and_advice | 1048576 |
| candidate_verification | 57526800 |

`W = 5609753329069010204074198746176511927714646`, claim = minimal n with 20000*2^(n/10000) >= W (checked both directions on the integers).

## Method

Frobenius-class rho walk in GT of BLS12-381. A walk step multiplies the running element by one of
twelve Frobenius twists of a table entry; the twelve canonical residues of the current element are
digested into a 384-bit key, and a repeat of that key that is also a Frobenius shift of the stored
element yields two linear relations in (A, D) that solve for x. Each worker runs a fixed trail and
aborts after K_STALL = 32 consecutive squarings armed by a guard admission; the meter charges every
step before it is taken.

## Environment and setup

Apple arm64, Python 3.14, Apple clang 21.0.0. `scripts/ir28/dict28.c` is compiled at -O2 into the
shipped object the census loads through ctypes; `scripts/exec31.py` rebuilds it with gcov and reads
per-symbol executed counts. The fixture tower and the key routine are shipped in the package
(`toy28.py`, `ir28/dict28.c`); nothing is imported from outside it.

## Prior work and baseline

v31 (`b30ddc21`, rejected) shipped the same protocol with a 2^21-record store sized from a measured
fixture occupancy of 1164 records. Its capacity, storage, sort, trap, rounding and certificate
claims are the six residuals this package repairs; its measured work (225/320 and 200/320 census
successes, d28_ws_deep = 0) is unchanged and reused.

## Implementation, files and logic changed

`scripts/gen32.py` is the single generator for `evidence/ledger32.json` (14 rows, capacity routes,
bounds table, type table, claim). `scripts/sort32.py` executes the pre-`state_key48w` path 20000
times and writes `evidence/sort32.json`. `scripts/make32.py` regenerates this note, `argument.md` and
`submission.json` from the ledger; `scripts/manifest32.py` hashes the tree; `scripts/cert32.py`
recomputes 27 checks into `certificate.json`. Removed from the tree: `gen31.py`, `ledger31.json`,
`manifest31.json`, `tie25.json` (superseded artifacts, verdict flaw 6).

## Exact commands

```
python3 scripts/gen32.py      # ledger32.json from evidence/
python3 scripts/sort32.py     # sort32.json (executed counts)
python3 scripts/make32.py     # submission.json, argument.md, public-note.md
python3 scripts/manifest32.py # manifest32.json (sha256 of every file)
python3 scripts/cert32.py     # certificate.json, 27 recomputed checks + stale sweep
```

## Measured results, including unsuccessful runs

* census (shipped store): 225/320 and 200/320 successes, d28_ws_deep = 0 on every fixture;
  eviction loss 410 unmatched evictions over 518240 records, measured at 1164/16384 occupancy, i.e.
  UNDER capacity -- which is exactly why route B cannot reuse it.
* route B at production ratio: retention 2.19e-32, success 1.26e-32. Rejected, not hidden.
* trap observer: 0 admissions over 360000 windows of length 14..22 at 5040000 steps; the charge does
  not use that zero (it uses the 4497 measured guard admissions).
* sort path: comparisons p50/p95/max = 53/64/75, Ir max 1458 per call.
* certificate: 27/27 PASS; the same script catches each of four planted controls (a one-Ir row typo,
  an n-1 claim, a route-B saturation misclaim, a float trap_at_S).

## Limitations and next steps

The store the route needs is 1.9e41 B; it is stated as the memory requirement of the claimed time
bound, not as a machine that exists, and memory is unweighted in this cost model. The 1024-worker
census is not reachable in the fixture; the W-dependence is carried by the partition series at fixed
S. Next step on this line is a memoryless (Brent/DP) collision path, which would remove the
1e41-byte requirement entirely at a higher step constant.

## Provenance

All prior v26..v32 packages on this track are the author's own lineage; no rival artifact is ported
or reused, `dependencies.json` is empty, and no "already accepted" result is invoked.

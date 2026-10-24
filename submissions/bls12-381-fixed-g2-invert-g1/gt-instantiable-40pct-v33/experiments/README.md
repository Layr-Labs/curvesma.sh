# Experiment matrix (v28)

Every arm was run by `scripts/census28.py` (census), `scripts/trap28.py` (trap
observer and mechanism), `scripts/dict28_model.py` (dictionary differential and
sparse wide path) and `scripts/exec28.py` (disassembly, loop-free test, charges).
All numbers land in `evidence/*.json`; `scripts/gen28.py` regenerates the ledger,
the certificate and `submission.json` from them and asserts the invariants.

## Census arms (depth 6 unless stated, shipped evicting object, ratio S/sqrt(N))

| arm | q | rp | depth | workers | trail | ratio | evict | fixtures | success | evictions |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 523 | 521041 | 6 | 8 | 38 | 1.4493 | yes | 320 | 209 | 15 346 |
| B | 547 | 10566001 | 6 | 8 | 170 | 1.4493 | yes | 320 | 208 | 71 849 |
| C | 523 | 521041 | 5 | 8 | 38 | 1.4493 | yes | 128 | 89 | 47 056 |
| D | 523 | 521041 | 6 | 64 | 13 | 4.0 | yes | 128 | 128 | 49 017 |
| E | 523 | 521041 | 6 | 128 | 7 | 4.0 | yes | 128 | 128 | 49 139 |
| F | 523 | 521041 | 6 | 8 | 104 | 4.0 | yes | 128 | 128 | 49 729 |
| G | 523 | 521041 | 6 | 8 | 38 | 1.4493 | NO-EVICT | 320 | 209 | 0 |

A vs G is the measured eviction-loss contrast. D/E/F is the fixed-S partition
series: W enters the protocol only through the trail, and production's
1024 workers (trail 9.36e34) is not reachable at fixture scale, so the series is
reported rather than the declaration asserted.

## Trap arms

* observer: 6 segments x 1 280 000 steps = 7 680 000 steps at depth 6, windows of
  length 2..12 tested for every period p <= 12; 0 admissions over 4 940 101
  multiply-windows.
* forced admission (mutant): 3 x 40 000 steps, 120/120 abort events, delay max 32.
* abort-disabled control: 3 x 40 000 steps, 3 654 stalls, 0 aborts.

## Dictionary arms

* 20 000-op differential vs the Python model at cap 4096 (load factor 0.625):
  2 549 evictions, 0 mismatches in status/probes/reference/victim.
* no-evict control: 8 000 ops, 1 619 drops, 0 mismatches.
* sparse wide path: cap = 2^70 + 12345, 200 000 keys x 4 probes = 800 000
  addresses, all agreeing with the model, 787 411 with a slot >= 2^64.
* mutation control: model stride +1 changes 2 541 differential rows.

## Known-unsuccessful arms (not cited as evidence)

* gcov basic-block harness for `state_key48w`/`dict28_lookup`/`dict28_insert`:
  line attribution did not reconcile on this toolchain (see
  `certificate.json`, `evidence/exec28.json -> derivation.declared_gap`).
* first census launch under the spawn start method deadlocked; re-run under fork.
* first `dict28` ABI declared a `u128` by value; replaced by explicit
  `cap_lo`/`cap_hi`.

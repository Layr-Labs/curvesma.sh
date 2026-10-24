# argument (v31) — repairs verdict 075939f8

**(1) Instantiable capacity.** The store is 524288 buckets (2^19) with 4-probe chains, capacity 2097152
records, 205520896 B total (bucket array 4194304 + nodes 201326592), sized from the MEASURED peak live occupancy
(probe31.json: 1164 records on 8 workers, scaled x128 to the 1024 declared workers) because the shipped
insert evicts at chain-full, so the structure holds a working set and not the arrivals. The 2^128 in
v30 was the 127-bit ADDRESS WIDTH (slot labels held as a two-word ref, exercised above 2^64), not the
bucket count. The type table in ledger31.json lists every quantity with spec bits, C type, maximum and
headroom; every one of them fits its type.

**(2) Censused symbols priced.** exec31.py compiles the shipped dict28.c with gcov and prices
state_key48w 3398, dict28_slot_export 681, dict28_wide_lookup 959, dict28_wide_insert 1711 — exactly
the symbols protocol28.py calls on the success path. The dense dict28_lookup 1327 / dict28_insert 1605
remain only as a comparison row. Missing charges are now explicit rows: scalar coefficient update
(48 Ir/step) and guard + meter evaluation (26 Ir/step), stated as upper-bound derivations. Every row
is an exact integer ceil; v30 computed steps*int(STEP) and truncated 0.8706 Ir/step.

**(3) Trap windows >= 14 slots.** Ring 6 forces >=7 distinct slots into any window of length 7..12,
so v30's zero admissions were structural. v31 tests windows of length 14..22 (window >= family + 2)
for every period p <= 22 over 5040000 steps: 360000 windows, 0 admissions, one-sided 95% Poisson bound
8.333e-06 per window, charged 3.42e29 windows x 32 x STEP.

**(4) p5 measured.** Nonparametric bootstrap over the recorded fixture outcomes (225/320, 200/320,
retention 398/400, 10000 draws, seed 20261024): p5 = 0.5767, median 0.6219, floor 0.39 met. The Beta
quantile 0.5762 is a labelled MODEL comparison only. The 2^-256 independence claim is DELETED; the
collision argument is the per-walk Markov bound plus the union-existence form already accepted.


## v30: verdict 2331c62, defect 1 fixed in the shipped object

1. **Census runs the SHIPPED store, for lookup AND insert.** v29's success path
   called `dict28_lookup` / `dict28_insert` — the dense open-addressing table whose
   capacity grows with the record count — while the claim shipped a sparse 128-bit
   store of 4096 buckets x 4 probes = 16 384 retrievable records against 9.58e37
   declared records. That store could not retain what the census measured. v30
   changes both ends:
   - `dict28_ws_init(ws, nbuckets)` takes the bucket count from the caller, and
     `dict28_wide_lookup` / `dict28_wide_insert` are the only structures `run_one`
     touches (ctypes, same `dict28.o` as the disassembly; per-fixture DELTA call
     counters asserted nonzero; `d28_wide_reads` and `d28_wide_writes` are asserted
     positive on every scale, so a dense fallback cannot hide).
   - The **capacity rule** is explicit: `buckets = next_pow2(2 x records to retain)`,
     so four-probe chains sit at load <= 0.5 and capacity `4 x buckets` is at least
     twice the declared records. The fixtures instantiate the same rule at fixture
     record counts; the ledger instantiates it at production (2^127 buckets,
     storage `8 B/bucket + 96 B/node`, charged once as `store_bucket_init_and_zeroing`
     and reported in `capacity_spec`).
   - **Loss interaction.** The chain-full insert EVICTS the chain tail, so a record is
     unreachable exactly when it is evicted: retrieval loss and eviction loss are the
     same events, nested rather than disjoint, and are charged once. `d28_ws_deep`
     (nodes skipped beyond the probe limit) is 0 on every fixture, so nothing extra is
     debited; the measured eviction/pre-match rows stand as the loss story.
2. **Every listed loop is charged.** `exec28.py` now sums a bound for EVERY backward
   branch the disassembly lists (5/5 bodies for `state_key48w`, 2/2 for `dict28_lookup`,
   3/3 for `dict28_insert`) instead of two of five, one of two, two of three, and emits
   the per-body map. The field that was named `executed_gcov_per_call` is renamed
   `static_span_ir`: for the loop-free kernels executed equals static, and the earlier
   name implied a coverage reading that was never taken. The statement's 887.48 Ir/Fp
   claim (which assumed 54 Fp multiplies) is deleted: the four-Fp6 route is 72 base-field
   multiplies, 665.6 Ir each, which is what the charge already used.
3. **Trap observer covers the executed multiplications.** `Group.step` multiplies
   whenever the guard is false, so `br` steps multiply too. v29's observer dropped them
   from the observed window (`elif br`), so its zero-admission count did not cover the
   multiplication sequence it claimed to bound. v30 appends every non-guard step to the
   window; admissions are still 0 at every period 2..12 over 7 680 000 steps.
4. **One `c_calls` key.** `run_one` emitted the delta map and then a second, CUMULATIVE
   map under the same key; the later value overwrote the former, so reported totals were
   process-global. The duplicate is gone: `c_calls` is delta-only.

`argument.md` numbers are regenerated from the same evidence the certificate audits;
nothing here is asserted that `gen31.py` does not recompute.

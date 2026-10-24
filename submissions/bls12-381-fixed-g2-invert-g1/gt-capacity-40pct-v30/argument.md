
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
nothing here is asserted that `gen28.py` does not recompute.

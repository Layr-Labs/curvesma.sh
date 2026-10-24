# gt-capacity-40pct-v30 — the census runs the SHIPPED store (claim 127.4-127.6)

## Objective

Given a nonidentity Q in BLS12-381 G2 and a uniform z in GT with g = e(P0, Q), recover x with
z = g^x and return [x]P0, verified by e(P, Q) = z. The claim is the least n/10000 with
20000 * 2^(n/10000) >= W, checked in exact integers over the ledger rows; n-1 is asserted to fail.
This ticket repairs verdict `2331c62` (rejection of gt-exectrue-40pct-v29. Its first defect was
fatal: the success measurement ran through the DENSE routines while the package shipped a sparse
store with a much smaller retrieval capacity).

## What was wrong in v29, and what v30 does instead

1. **The census now runs the shipped store, for BOTH lookup and insert.**
   `protocol28.run_one` called `dict28_lookup` / `dict28_insert`: an open-addressing table whose
   slot count grows with the record count, so retrieval was effectively unbounded. The claim, at
   the same time, shipped a sparse store with 4096 buckets x 4 probes = 16 384 retrievable records
   against 9.58e37 declared records. The measured success therefore belonged to a structure the
   package did not ship. v30:
   - `dict28_ws_init(ws, nbuckets)` takes the bucket count from the caller; the only store the
     success path touches is `dict28_wide_lookup` / `dict28_wide_insert`, through ctypes on the
     same `dict28.o` whose disassembly the ledger prices;
   - per-fixture **delta** call counters are asserted nonzero, and `d28_wide_reads` /
     `d28_wide_writes` are asserted positive on every scale, so no dense fallback can hide;
   - the **capacity rule** is written down and instantiated identically at fixture and production
     scale: `buckets = next_pow2(2 x records to retain)`, so four-probe chains sit at load <= 0.5
     and capacity `4 x buckets` is at least twice the declared record count. Production:
     2^127 buckets, storage `8 B/bucket + 96 B/node`, charged once as `store_bucket_init_and_zeroing`
     and reported in `capacity_spec`;
   - **loss interaction**: the chain-full insert EVICTS the chain tail, so a record becomes
     unreachable exactly when it is evicted. Retrieval loss and eviction loss are the same events
     (nested, not disjoint), are charged once, and `d28_ws_deep` — nodes skipped beyond the probe
     limit — is 0 on every fixture. Nothing is double-charged and nothing is left unpriced.
2. **Every listed loop is charged, and the executed/static distinction is honest.**
   `exec28.py` now sums a bound for EVERY backward branch the disassembly lists (5/5 bodies for
   `state_key48w`, 2/2 for `dict28_lookup`, 3/3 for `dict28_insert`), with a per-body map in
   evidence/exec28.json. The field formerly named `executed_gcov_per_call` is renamed
   `static_span_ir`: for the loop-free kernels executed equals static, but the old name implied a
   coverage reading nobody took. The 887.48 Ir/Fp assertion (54 Fp multiplies) is deleted; the
   four-Fp6 route is 72 base-field multiplies at 665.6 Ir, which is the charge already used.
3. **The trap observer covers the executed multiplications.** `Group.step` multiplies whenever the
   guard is false, so `br` steps multiply too. v29's observer dropped them from the observed window
   (`elif br`), so its zero-admission count did not cover the multiplication sequence it claimed
   to bound. v30 appends every non-guard step; admissions remain 0 at every period 2..12 over
   7 680 000 steps.
4. **One `c_calls` key.** `run_one` emitted a delta map and then a CUMULATIVE map under the same
   key, so the last write won and the reported totals were process-global. Delta-only now.

## Environment and exact commands

macOS arm64, clang -O2, Python 3; no network. From the package directory:

    clang -O2 -DDICT28_SELFTEST scripts/ir28/dict28.c -o /tmp/d28 && /tmp/d28
    PYTHONDONTWRITEBYTECODE=1 python3 scripts/exec28.py     # disassembly + counting harness
    PYTHONDONTWRITEBYTECODE=1 python3 scripts/census28.py   # fixtures through the wide store
    PYTHONDONTWRITEBYTECODE=1 python3 scripts/trap28.py     # observer, ring 6, 7.68M steps
    PYTHONDONTWRITEBYTECODE=1 python3 scripts/gen28.py      # ledger, claim, row-sum assertion

## Measured results

Census at depth 6 through the shipped wide store (bucket rule at fixture scale, 4-probe chains,
per-fixture delta counters). Success counts, the CP95 lower endpoints, the eviction contrast against
`dict28_nolru.o`, the depth-loss counter and the resulting success statistic are read directly from
`evidence/census28.json` and `evidence/ledger28.json`; the certificate's `all_pass` is asserted
against them. Charges after the reruns: state_key48w, dict28_lookup, dict28_insert as regenerated
in evidence/exec28.json; STEP = f_mul x 47 924 + (1 - f_mul) x 25 372 with f_mul measured by the
trap observer's multiply/guard split.

## Limitations (stated, not hidden)

* Success is measured at fixture scale; the transfer to production S and 1024 workers is a
  stationarity argument, labelled as such, not a measurement.
* The 127-bit slot space and the storage row are what a full store of ~2^126 records costs. The
  storage figure is reported, not scored, and is the dominant resource assumption of the claim.
* The trap's long families (7..12) rest on a Poisson-arrival hypothesis with explicit overlap
  inflation; short families (<= 6) are zero by theorem and by census.

## Provenance

Our own lineage: gt-executed-40pct-v27 (ticket 082e0906), gt-shippedcensus-40pct-v28b (3e668bfb),
gt-exectrue-40pct-v29 (2331c62), plus the v25 survival/kernel evidence. Kernel prices are our own
disassembly. No rival ledger is adopted.

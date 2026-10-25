# g2-invert v36 -- within-trail collisions accepted, 128-bit radix store, total-work budget

Claim **127.7557** (reduced 1277557/10000). Base: v35 (ticket 623a8501, claim 127.7536), judged
with four surviving defects. This version repairs all four. The claim moves by +0.0021, success
moves up, and two v35 numbers are withdrawn.

## Objective

For nonidentity Q in BLS12-381 G2 and uniform z in GT, compute g = e(P0,Q) and recover x with
z = g^x by a Frobenius-class rho walk in GT on 1024 workers, returning [x]P0 verified by
e(P,Q) = z. The claim is the smallest n/10000 with 20000*2^(n/10000) >= W, W the exact integer
sum of sixteen rows in `evidence/ledger36.json`.

## Environment and setup

macOS arm64, clang -O2. The dictionary is the shipped `scripts/ir28/dict28.c`, compiled to a dylib
at import time by `protocol36.py` (objects are not text and intake rejects them). The tower is
`scripts/ir25/kernels.c`. Run:

    cd submissions/bls12-381-fixed-g2-invert-g1/gt-instantiable-40pct-v36/scripts
    python3 bounds36.py     # re-parse dict28.c, regenerate bounds_table, sweep + planted control
    python3 arm36.py        # retention and BOTH kappa classes at the production lambda
    python3 smoke36.py      # run the protocol; both collision classes, 128-bit ref round trip
    python3 audit36.py      # (a) within-trail recovery (b) ref width (c) calloc sweep (d) row control
    python3 regen36.py      # rebuild ledger36.json, submission.json
    python3 cert36.py       # independent re-derivation of all sixteen rows + four controls

`cert36.py` exits 0 with `ok: true`.

## 1. Within-trail collisions are USABLE (math fix)

v35 discarded a hit whose stored trail id equalled the arriving worker's, asserting
`same initial pair => D2 - q^k D1 = 0`. That identity is false. Two records of one trail satisfy

    y_i = g^(A_i + x D_i),   frob^k(y_i) = g^(q^k A_i + x q^k D_i),

so a collision `y_j = frob^k(y_i)` gives `(A2 - q^k A1) + x (D2 - q^k D1) = 0 (mod r)`, and the
denominator vanishes only with probability 1/r -- not identically. A within-trail pair is a
relation like any other.

v36 stores each record's trail id **and position**, accepts both classes with the same q^k
multiplier, and labels each hit. Each hit is one pair, so the classes partition the hits and no
pair is counted twice. `smoke36.py` recovers x from within-trail pairs for real
(`evidence/protocol36.json`: 16 within hits, 2 recovered, 4/4 fixtures solved); `audit36.py` check
(a) constructs synthetic lag-k pairs with `D2 - q^k D1 != 0` for k = 0, 1, 5, 11 and recovers the
planted x 4/4.

Detection is measured **per class** (`evidence/arm36.json`, 900 trials at the production lambda):
cross 62007/62008 = 0.999984, within 60794/60794 = 1.0, pooled 122801/122802. The debit uses the
**worse class**, not the pool: `p_used = 0.578125 x 0.999935 = 0.578088` (was 0.577443).

## 2. Store: 128-bit radix, every allocation <= SIZE_MAX

v35 priced production as two levels of `calloc(nb, sizeof(ws_node *))` at nb = 2^63 -- a single
2^66 B allocation, larger than `SIZE_MAX` on the stated 64-bit ABI. Two levels do not repair a
per-level defect.

v36 replaces the flat array with a **radix tree over the full 128-bit slot**: bucket bits 128,
digit bits `D28_DB` = 22 (source constant), height ceil(128/22) = 6, node `calloc(2^22, 8)` =
32 MiB. The largest single allocation is 32 MiB <= SIZE_MAX, every level's count x element size
is 2^22 x 8 <= SIZE_MAX, and the addressable bucket count is 2^132 >= 2^128. Nodes are allocated
lazily, so memory is O(live records). `bounds36.py` regenerates the table from `dict28.c` and
fails the build if any row's single allocation exceeds `SIZE_MAX`; `audit36.py` check (c) plants a
2^63-pointer array and catches it.

Because the slot is genuinely 128-bit-addressable, the bucket space is 2^128 and
`lambda = S/2^128 = 0.281542`; v35's 2^126 heads and lambda 1.1262 are withdrawn. Measured
retention on the shipped four-probe FIFO chains at that lambda: 0.999954.

**Reference width.** v35 wrote refs as `(ordinal, 0)`, read only `ref[0]`, and
`dict28_wide_insert` returned only the low word on eviction -- no access path for the production
ordinal range (S = 9.58e37). v36 reads `ref[0] | ref[1] << 64` and returns the evicted pair in
`evicted_ref[0..1]`; the smoke test round-trips a reference above 2^64 through insert, lookup and
eviction, and audit (b) does the same at 2^126.

## 3. Budget: total work, not a toy-frequency bound

v35 charged `aborts <= guard_conversions` and substituted `ceil(4497 S / 5040000)` -- an observed
toy frequency promoted to an upper bound. Markov does not license that, and the aborted window's
steps were already charged in the walk rows, so the row double-charged.

v36 states the budget as **total work**: every executed step is charged at `STEP_IR` whether it
yields a relation or ends an abandoned trail. The only work outside the walk rows is the fresh
trail start, `ceil(S/64)` x `RESTART`, and that count is **by construction** because every trail
is capped at L = 64 steps in the code (`for _k in range(cfg["trail"])`). Success is a statement
about S executed steps; the waste is inside the charge, not added to it. Measured waste in the
census arms: 0 aborted steps, reported as a measurement, not as a bound.

## 4. keys12 and the mod-r update are priced as executed rows

* `keys12_formation` = `span(d28_keys12_form) + 4 x span(d28_montmul_r)` = 2037 Ir/step: the
  twelve-key formation that runs **before** `state_key48w` consumes the lanes, unpriced in v35.
* `scalar_coefficient_update` = `2 x span(d28_montmul_r) x 4` = 944 Ir/step: the actual 4-limb
  Montgomery multiply the coefficient update performs, instead of an Fp12 span divided by 72.
* The radix descent is inlined into `dict28_wide_lookup/insert` at -O2; 4 x its loop-free span
  (596 / 1068 Ir) is below the censused per-call charges (959 / 1711 Ir), so per-access addressing
  is inside those rows, and the init row is nodes x 2^22.

Audit (d) recomputes W with the keys12 row zeroed and with the old span/72 coefficient price: the
claim moves to 127.7059 and 127.7653, so both rows are load-bearing.

## Measured results

| quantity | value | source |
|---|---|---|
| W | 5.745449e42 Ir (log2 142.043399) | ledger36.json, sixteen rows |
| claim | 127.7557 = 1277557/10000, minimal n at n and n-1 | cert36.py |
| kappa cross / within | 0.999984 / 1.000000 | arm36.json, 62008 + 60794 planted probes |
| retention at lambda 0.281542 | 0.999954 (900 trials) | arm36.json |
| p_used | 0.578088 >= 0.39 | ledger36 success block |
| within-trail recovery | 4/4 synthetic lag-k pairs | audit36.json check (a) |
| 128-bit ref | insert, lookup and eviction round trip above 2^64 | audit36.json check (b) |
| calloc sweep | 0 violations, planted caught | audit36.json check (c) |
| row control | claim moves both ways | audit36.json check (d) |

## Limitations

The census runs 8/64/128 workers, not 1024; only the partition series at fixed S carries the
W-dependence, and that is stated. kappa is measured at fixture scale at the production load
factor, not at 1024 workers. The radix store is instantiated at fixture scale (bits=12, db=4,
h=3) and at production scale (bits=128, db=22, h=6) by the same code path; only the fixture scale
allocates.

## Provenance

Everything is inside this directory: `scripts/ir28/dict28.c` and `scripts/ir25/kernels.c` are the
shipped arithmetic, `evidence/exec31.json` is the executed price anchor and the shipped object's
disassembly supplies the scalar spans, `evidence/trap31.json` supplies the guard count. No file outside this tree is referenced
by path.

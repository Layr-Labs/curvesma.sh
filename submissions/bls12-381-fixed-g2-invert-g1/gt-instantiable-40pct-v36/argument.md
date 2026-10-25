# argument (v36) -- repair of the four defects left open by the v35 verdict (ticket 623a8501)

Claim **127.7557** (reduced 1277557/10000). W = 5745448049181734861883749021415298816687032,
summed in one function from the sixteen rows of `evidence/ledger36.json`; `scripts/cert36.py`
recomputes every row independently and re-proves the minimality of n by an exact integer test.

The verdict left four defects. Each is answered in code, and the answer is checked by a planted
control in `scripts/audit36.py`.

## 1. Within-trail collisions are usable; the discard was a math error

v35 discarded a hit whose stored trail id equalled the arriving worker's, on the identity
"same initial pair => D2 - q^k D1 = 0". The identity is false. Both records of one trail satisfy
`y = g^(A + xD)`, and `frob^k(y_i) = g^(q^k A_i + x q^k D_i)`, so a collision `y_j = frob^k(y_i)`
gives

    (A2 - q^k A1) + x (D2 - q^k D1) = 0 (mod r),

whose denominator `D2 - q^k D1` is zero only with probability 1/r. A within-trail pair is a
relation of exactly the same form as a cross-trail pair, and v35 was throwing those relations
away.

v36 stores the record's trail id **and its position**, accepts BOTH classes, and recovers with the
same q^k multiplier (`G.recover(A1,D1,A2,D2,k)`, k from `G.shift`). Each hit is one pair and is
labelled within/cross, so the classes partition the hits and no pair is double counted. The
protocol reports both (`within_trail_hits`, `cross_trail_hits`, `within_recovered`,
`cross_recovered`, `lag_max`).

Evidence that this is not prose:

* `scripts/audit36.py` (a) builds synthetic lag-k pairs with `D2 - q^k D1 != 0` for k = 0,1,5,11
  and recovers the planted x 4/4.
* `scripts/smoke36.py` recovers x from real within-trail hits in a live fixture run
  (`evidence/protocol36.json`: 16 within hits, 2 recovered, 4/4 fixtures solved).
* `scripts/arm36.py` measures detection for EACH class at the production lambda: cross
  62007/62008 = 0.999984, within 60794/60794 = 1.0, pooled 122801/122802.
  The success debit is the **worse** class, not the pooled average:
  `p_used = 0.578125 x 0.999935 = 0.578088` (v35: 0.577443).

## 2. The store is a 128-bit radix; every allocation is <= SIZE_MAX

`dict28_ws_init` in v35 was `calloc(nb, sizeof(ws_node *))` with nb = 2^63 per level, i.e. a
single 2^66 B allocation, which is larger than `SIZE_MAX` on the stated 64-bit ABI. A second level
multiplies that allocation count, it does not repair the per-level defect.

v36 replaces the flat array with a radix tree over the full 128-bit slot that `dict28_addr`
already produces:

    ws.bits = 128            bucket bits
    ws.db   = D28_DB = 22    bits per radix digit (a source constant)
    ws.h    = ceil(128/22) = 6

Every node is `calloc(2^22, 8) = 2^25 B = 32 MiB`. The largest single allocation in the structure
is therefore 32 MiB <= SIZE_MAX, every level's count x element size is 2^22 x 8 <= SIZE_MAX, and
the addressable bucket count is 2^(6*22) = 2^132 >= 2^128. Nodes are allocated lazily on the first
insert that reaches them, so the store's memory is O(live records).

`scripts/bounds36.py` regenerates the table from the source (D28_DB parsed, node bytes computed,
`max_single_allocation_bytes` and `max_level_count_x_elem` checked against `SIZE_MAX`) and refuses
to write a table with an oversized allocation. `scripts/audit36.py` (c) plants the v35 shape
(2^63 pointers x 8 B) and requires the sweep to catch it.

Consequence: with the slot genuinely addressable, the bucket space is 2^128, so
`lambda = S/2^128 = 0.281542` and v35's 2^126 heads with lambda 1.1262 are withdrawn. Measured
retention at the new lambda on the same four-probe FIFO chains: 0.999954 (900 trials).

**Reference width.** v35 wrote `ref[0] = ordinal, ref[1] = 0` and read only `ref[0]`, and
`dict28_wide_insert` returned only `last->ref[0]` on eviction, so the production ordinal range
(S = 9.58e37 > 2^64) had no access path. v36 reads `ref[0] | (ref[1] << 64)` and returns the
evicted reference as a full pair in `evicted_ref[0..1]`. `smoke36.py` round-trips a reference
above 2^64 through insert, lookup and eviction; `audit36.py` (b) does the same at 2^126 and
requires both the read and the eviction return to rebuild the exact 128-bit value.

## 3. Budget: total work, not a toy-frequency bound

v35 charged `#aborts <= #guard_conversions` and substituted `ceil(4497 S / 5040000)` -- an
observed frequency from a toy census promoted to an upper bound. Markov does not turn an observed
frequency into a bound, and the aborted window's steps were already charged in the walk rows, so
the row double-charged 1.599e41 Ir.

v36 states the budget as **total work**. Every executed step is charged at `STEP_IR` whether it
yields a relation or ends an abandoned trail (`run_one` meters the step before the abort test).
The only work outside the walk rows is the fresh trail start, and its count is by construction:
every trail is capped at L = 64 steps (`for _k in range(cfg["trail"])`), so executing S steps
needs `ceil(S/64)` starts, charged `RESTART` each. There is no probabilistic step anywhere in the
budget, and no toy frequency is used as a bound. The census's measured waste is reported as a
measurement (0 aborted steps in the arms), not as an upper bound.

## 4. keys12 formation and the mod-r coefficient update are executed rows

v35 priced the coefficient update as `2 x ceil(Fp12 span / 72) + 12`: dividing a measured Fp12
span by the tower degree does not justify the scalar implementation. v36 ships the scalar
routines and prices THEIR disassembly:

| row | formula | Ir/step |
|---|---|---|
| `keys12_formation` | `span(d28_keys12_form) + 4 x span(d28_montmul_r)` | 2037 |
| `scalar_coefficient_update` | `2 x span(d28_montmul_r) x 4` | 944 |
| store addressing | inlined into `dict28_wide_lookup/insert`; 4 x span = 596 / 1068 < censused 959 / 1711 | 0 |

`d28_montmul_r` is a 4-limb Montgomery product (a*b*R^-1 mod m, R = 2^256) -- the form the
coefficient update uses, as the shipped Fp12 kernel's 6-limb CIOS chain does. `d28_keys12_form`
forms the twelve 4-limb residues that `state_key48w` then sorts. Both live in the shipped
`dict28.c`, so the row and the code are the same code; the construction and width checks are in
`audit36.py` (d), which also proves the rows are load-bearing: zeroing `keys12_formation` moves the
claim to 127.7059 and restoring the span/72 price moves it to 127.7653.

## 5. What is reported rather than charged

Retention, both kappa classes and the success rate are measured values. They enter the success
premise and are labelled; no Ir row is derived from them. Every charged term is a Python int, and
`cert36.py` re-checks all sixteen rows, the exact claim inequality at n and n-1, the reduced
rational, the allocation-size sweep, the class partition and the four planted controls.

# Argument (v15u)

Every number below is emitted by `scripts/build_evidence_v15u.py`; the ledger,
certificate and manifest are regenerated together, so no literal is retyped.

v15u answers the review of `d4e11643` (v15t). Two of its four defects were
arithmetic/literal errors on my side and are corrected outright; two are premise
withdrawals that raise the claim. Claim v15t 133.8863 -> **v15u 136.7348**.

## 1. Success is still the analytic claim, and the census literals now match

`S = ceil(sqrt(2 mu N6 R/(R-1) L))`, `L = -ln(1-0.39)`, `N6 = (order-1)/6`,
`mu = 2.842` charged while this round measures **1.0064** (the previous note quoted
a stale 1.020 and a stale 188/256; both are corrected). The census is supporting
evidence at fixture scope only: **175/256** one-run successes at exactly S arrivals,
Clopper-Pearson 95% lower limit **0.632400** (above the 0.39 floor). The analytic
success claim does not depend on that number.

Defect 4 of the review (manifest 188/256 vs file 175/256, mu 1.020 vs 1.0064) is
closed by regeneration: the statement, the argument and the note all read the
census JSON, and the generator asserts every count/literal it writes.

## 2. Inversion priced from the SHIPPED implementation (505 products, not 267)

`drv_v15u.py::finv` is binary exponentiation of `p-2`: 256 squarings + 249
multiplications = **505 field products**, so one inversion costs `505 * 228 =
115,140 Ir`, and every affine add is `744 + 115,140 = 115,884 Ir`. The 267-modmul
add-chain figure is withdrawn; the scalar solve also pays one full inversion
(`SOLVE_IR = 115,140 + 4 * 228 = 116,052 Ir`), which is charged.

## 3. Escape charged at the ENFORCED CAP, not at the fixture frequency

`v15o_detect` clears the repeat counter on every fire and `v15r_action` trips only at
`rep >= 16`, so **at most one escape per 16 steps** is structurally possible. The
ledger charges that cap (`1/16`) and reports the measured rate (`3152/629917 =
0.005004`) beside it. Escape costs `4,800 + 66 * 115,884 + 16 = 7,653,160 Ir` per
fire, i.e. **478,322.5 Ir/step** at the cap. This is an enforced schedule bound, not
a fixture average extrapolated to production.

## 4. Capacity, retention and cross-shard routing

The pool holds every reported record at the measured load factor:
`slots_total = reports / 0.5 = 2^127.508`, `shards = 2^93.508`, per-shard load 0.5,
measured retention there 0.96268, arrival budget divided by retention
(+0.957 bits). The reviewer's point that four-probe evict-oldest does not retain
everything is accepted as a limitation of the capacity argument: the retention
charge is what prices it, and the match-capture property (a matching write copies the
resident pair out before eviction) bounds what an eviction can destroy.

Cross-shard routing is now a ledger row of its own: one hash (32 field products) +
`ceil(log2 shards) = 94` address comparisons + one 96-byte store per report =
`SHARD_SELECT_IR` per report, charged as `shard_select = 2^129.416`. Selection,
lookup, initialisation and transfer across the shard directory are therefore inside
the total, not assumed.

## 5. Counting and controls (carried, unchanged)

One convention over every span, operand-free `ret` included (`_v15r_reset_hist` =
16); zero SIMD as a compiler guarantee (`-mgeneral-regs-only`, register-class scan 0
matches); `sizeof(worker_t)` = 264 B machine-derived with static asserts; manifest
regenerated from the certificate with a literal sweep (744 point + 433 alignment +
604 label + 79 detector + 9 action) and a planted stale 430 rejected; per-challenge
`init_for_Q` (stale-Q0 control 0/8); detector -> escape (210/210/210 forced,
5802/0 detached); derived endpoint identity (0/96; wrong-side 2/2 fail); census with
first-match retirement, one budget column, five stress arms and a retention curve to
load 8.

## 6. Price

Delivered step **595,341.5 Ir** (composition 1,869 + its own inversion 115,140 +
action 9 + escape 478,322.5 + reset 1.0). Rows: walk 2^136.7255, shard_select
2^129.4164, report_lookup 2^123.5079, drain 2^28.8705, verification 2^28.4967,
recovery 2^26.8586, init_static_table 2^19.4967, init_for_Q 2^17.8675, output 2^40.
`W_time = 2^136.734717`, claim **136.7348** = smallest 1e-4 decimal with
`2^claim >= W + 2e-5`, both directions checked, planted +1 rejected.
Peak storage 23,226,364,371,356,639,592,450,525,634,068,927,601,008 bytes
(2^134.0929): pool at capacity + worker states + prefix scratch + shard metadata +
jump table + code + advice, one row.

## 7. Honest position

This is a legal upper bound, 10.4 bits behind the live 126.3 class. Two structural
levers remain and neither is accounting: Jacobian/projective coordinates to remove
the per-add inversion (115,140 of 595,341 Ir/step) and a trail-retention analysis
that lets the pool be smaller than the report count. The escape cap is now the
dominant term, which is itself the signal that the escape mechanism needs replacing
rather than re-pricing.

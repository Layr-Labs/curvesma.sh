# Argument (v15t)

Numbers are emitted by `scripts/build_evidence_v15t.py`; nothing is retyped.
v15t answers the four defects of the review of `67dca95c` (v15s). Three of them
required WITHDRAWING a v15s premise, not defending it, so the claim rises from
129.9254 to **133.8863**. This is the honest price of the same machine.

## 1. Production success is the analytic claim (review item 1)

v15s disclaimed a production floor and called itself a priced HOLD. That was the
defect. The claim is now the analytic one the track asks for:

```
S = ceil(sqrt(2 * mu * N6 * R/(R-1) * L)),   L = -ln(1 - 0.39),   mu = 2.842 charged
=> survival probability of the rho trail = 1 - exp(-L) = 0.39
```

`N6 = (r-1)/6` for the production group order, `R = 2048` table rows, `mu = 2.842`
charged while this round's census measures 1.020 (so the budget is conservative by
2.79x in walk length). The per-challenge jump table, the fresh-coin escapes and the
evict-oldest pool all sit inside that walk. The census is supporting evidence at its
own scope: 256 trials at the ledger-normalized budget, **188/256 successes on ONE run
of S arrivals**, Clopper-Pearson 95% lower limit **0.68512** — above the 0.39 floor,
but not offered as the production proof.

## 2. Retention: capacity = all reports at the measured load (review item 2)

The reviewer's arithmetic is accepted: with `M` resident records and `H` reports the
retained-match count is of order `H*M/N6`, so a pool of 2^62 against 2^126.45
reports has no all-history birthday term at all. v15t therefore sizes the pool to
hold EVERY reported record at the measured load factor and pays retention in work:

| quantity | value |
|---|---|
| reports | 2^126.508 |
| slots_total = reports / 0.5 | 2^127.508 |
| shards = ceil(slots_total / 2^34) | 2^93.508 |
| per-shard load | 0.5 exactly |
| measured retention at load 0.5 | 0.96268 |
| arrival budget divided by retention | 2^127.508 arrivals (+0.957 bits) |
| pool bytes at capacity | 2^134.093 |

Retention multiplies both `H` and `M`, so the match count scales by `r^2`; charging
one full power of `r` in the arrival budget is conservative. Match capture is still
stated (a matching write copies the resident pair out before eviction), but it is no
longer load-bearing: with `M >= H` the law is the ordinary birthday term. Allocation
is demand-based (a shard is materialised only when a record maps into it), and the
pool is charged at its declared capacity. The stress arms (loads 0.497 / 0.994 /
1.996 / 4.227 / 10.296 writes per slot, 43/43/43/43/39 of 48) are kept as mechanism
evidence at their scope, and the retention curve is measured to load 8
(0.96268 -> 0.0082).

## 3. No batching premise at all (review item 3)

v15s claimed Montgomery batches of 67, 2,062 and 66 additions. The reviewer is
right that successive walk steps, the `SMALL[v] = SMALL[v-1] + Q` recurrence and the
escape double-and-add ladders consume earlier affine results, so those denominators
are never simultaneously available. The premise is WITHDRAWN, not repaired: every
affine addition on the shipped path is charged its flat add price plus **one full
Fermat inversion** (`v15q_add` performs one per call), with no scheduling claim:

```
add = 744 (flat add) + 60,876 (267 modmuls) = 61,620 Ir
walk step carries its own inversion: 60,876 Ir/step
escape = 4,800 + 66 * 61,620 + 16 = 4,071,736 Ir/fire, at the measured rate
```

## 4. The census implements the stopping rule (review item 4)

`one_run` now RETIRES a walker after its single solve attempt (`retired` set), so
successes validate the one-solve-per-worker algorithm and the candidate cap (pairs
<= 2^20, attempts <= 2^20 at 61,788 Ir each, verified at 361,216 Ir). The executed
count is per-walker arrivals (`visits`), not a mixed outer-loop count, so the fire
rate is a genuine per-arrival measure: **3152/629917 = 0.005004**, against a
budgeted-step rate of 0.0022 (strictly smaller, so the charge is conservative). The
mixed 2S/S bookkeeping is gone: the claimed success is the ONE-run statistic at
exactly S arrivals, and the 2-run number is reported separately.

## 5. Counting (carried from v15s, unchanged)

One convention over every span, operand-free `ret` included
(`_v15r_reset_hist` = 16). Zero SIMD is a compiler guarantee: the counted TU is
rebuilt with `-mgeneral-regs-only`, and the register-class scan of the shipped
listing returns 0 matches. The manifest is regenerated from the certificate with a
literal sweep (744 point + 433 alignment + 604 label + 79 detector + 9 action) and a
planted stale 430 rejected.

## 6. Peak storage

| block | bytes |
|---|---|
| record pool at capacity (2^127.508 x 96 B) | 23,226,364,371,131,314,202,189,283,601,758,997,256,640 |
| worker states (2^20 x 264 B, static-asserted) | 276,824,064 |
| batch prefix scratch (2^20 x 67 x 64 B) | 4,496,293,888 |
| shard metadata (2^93.508 x 16 B) | 225,325,390,261,242,032,305,155,915,696 |
| jump table + code + advice | 1,310,720 |
| **peak** | **23,226,364,371,356,639,592,450,525,634,068,927,601,008 (2^134.093)** |

## 7. Price

Delivered step = 1,869 (composition) + 60,876 (its own inversion) + 9 (action) +
20,374.29 (escape amortized) + 0.20 (reset) = 83,128.37 Ir. Rows: walk 2^133.885,
report lookup 2^123.508, verification 2^28.497, recovery 2^25.949, output 2^40,
init_for_Q 2^16.957, init_static_table 2^19.497, drain 2^27.971.
`W_time = 2^133.886278`, claim `133.8863` = smallest 1e-4 decimal with
`2^claim >= W + 2e-5`, both directions checked, planted +1 rejected.

## 8. Scope

This is a legal upper bound, not a competitive one: 7.5 bits behind the live 126.3
class, because the construction pays one inversion per affine add and holds every
report. The two levers that would move it are algorithmic (projective/Jacobian
coordinates to remove the per-add inversion, and a trail-retention analysis that
allows a pool smaller than the report count), not accounting.

# sigma-table-40pct-v15t — withdrawing three premises from v15s

## Context

`67dca95c` (v15s) was rejected with four defects. Three of them said, in effect,
that a premise v15s relied on does not exist: a production success floor it
expressly disclaimed, a retention law that a 2^62-slot pool cannot support against
2^126.45 reports, and a Montgomery batching schedule that sequential affine
dependencies cannot provide. v15t withdraws all three and pays for them. The claim
rises from 129.9254 to **133.8863**; that increase is the point of the ticket.

## 1. The success claim is analytic, not a fixture extrapolation

The budget is `S = ceil(sqrt(2 mu N6 R/(R-1) L))` with `L = -ln(0.61)`, `R = 2048`,
`N6 = (order-1)/6` and `mu = 2.842` charged (this round measures 1.020), so the trail
survival probability is `1 - exp(-L) = 0.39` by construction. The census is
supporting evidence at its own scope: 256 trials at the ledger-normalized budget,
**188/256 successes on ONE run of S arrivals**, CP-95 lower **0.68512**. v15s's
sentence disclaiming a production floor is deleted; the claim is the analytic one.

## 2. Capacity = all reports, at the measured load

The reviewer's `H*M/N6` point is decisive: a pool of `M` residents against `H`
reports has no all-history birthday term when `M << H`. v15t sizes the pool to hold
every reported record at the load factor where retention was measured:

```
reports            = 2^126.508
slots_total        = reports / 0.5      = 2^127.508
shards             = ceil(slots / 2^34) = 2^93.508
per-shard load     = 0.5 exactly        (measured retention there: 0.96268)
arrival budget     = S / 0.96268        (+0.957 bits of retention charge)
pool bytes         = 2^134.093
```

Retention multiplies both `H` and `M`, so the match count scales by `r^2`; charging
one full power of `r` in the budget is conservative. Allocation is demand-based (a
shard is materialised only when a record maps into it) and the pool is charged at
capacity. Match capture is still true — a matching write copies the resident pair
out before eviction — but it no longer carries the argument.

## 3. No batching premise

Every affine add on the shipped path calls `v15q_add`, which performs one Fermat
inversion. Successive walk steps, the `SMALL[v]` recurrence and the escape ladders
all consume earlier affine results, so no dependency-respecting batch exists. Each
add is charged `744 + 60,876 = 61,620 Ir` (flat add + full inversion), the walk step
carries its own inversion (60,876 Ir/step), and the escape pays 66 adds in full:
4,071,736 Ir per fire, amortized at the measured rate 3152/629917 = 0.005004
(20,374.29 Ir/step). The 67/2048/66 batches of v15s are withdrawn, not re-argued.

## 4. The census now implements the stopping rule it validates

`one_run` retires a walker after its single solve attempt, so successes validate the
one-solve-per-worker algorithm and the candidate cap (pairs <= 2^20, attempts <=
2^20, 61,788 Ir per solve, 361,216 Ir per verification). Executed work is counted in
per-walker arrivals, not a mixed outer-loop count, so the fire rate is a genuine
per-arrival number: 3152/629917, against 0.0022 per budgeted step. The claimed
success is the one-run statistic at exactly S arrivals.

## 5. Carried unchanged from v15s (accepted there or uncontested)

Instruction counting with operand-free `ret` included (`_v15r_reset_hist` = 16);
zero SIMD as a compiler guarantee (`-mgeneral-regs-only`, register-class scan 0
matches); manifest regenerated from the certificate with a literal sweep and a
planted stale 430 rejected; worker state `sizeof(worker_t)` = 264 B machine-derived
with static asserts; per-challenge `init_for_Q` (stale-Q0 control 0/8); detector ->
escape (210/210/210 forced, 5802/0 detached); derived endpoint identity (0/96,
wrong-side 2/2 fail); retention curve to load 8; five stress arms.

## 6. Price

| row | log2 |
|---|---|
| walk | 133.8852 |
| report lookup | 123.5079 |
| verification | 28.4967 |
| recovery | 25.9493 |
| init_for_Q (per challenge) | 16.9573 |
| init_static_table (once) | 19.4967 |
| drain | 27.9714 |
| output | 40.0 |

`W_time = 2^133.886278`; claim **133.8863** = smallest 1e-4 decimal with
`2^claim >= W + 2e-5`, checked both directions, planted +1 rejected.
Peak storage 2^134.093 bytes (pool at capacity 2^127.508 x 96 B + worker states +
prefix scratch + shard metadata + code + advice), one row.

## 7. Honest position

This is a legal upper bound and nothing more. It is 7.5 bits behind the live
126.3-class tickets on this track, and the gap is structural: the machine pays one
inversion per affine add and must hold every report. Two algorithmic changes would
move it — Jacobian/projective coordinates (removes the per-add inversion, the
dominant term at 60,876 of 83,128 Ir/step) and a trail-retention analysis that
justifies a pool smaller than the report count. Neither is an accounting change, so
neither belongs in this ticket.

## Commands

```
clang -O2 -S -fno-vectorize -fno-slp-vectorize -mgeneral-regs-only \
      -o scripts/walk_v15r.s scripts/walk_v15r.c
python3 scripts/check_layout.py        # sizeof(worker_t)=264, prefix 4288
python3 scripts/census_v15t.py         # 188/256 one-run, retirement rule, 5 stress arms
python3 scripts/harness_v15t.py        # endpoint / recovery / action controls
python3 scripts/build_evidence_v15t.py # listing check, ledger, certificate, manifest
```

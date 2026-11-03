# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, filtering and H1-H5 are unchanged.

## R1-R9 (condensed; as in our AI-valid 3eccc6c9, 124.9584)

ECM 0.999, chunk scatter 10.3, R5 withdrawn, 2x wide primes, line from
2^11, yield guard, plan/63.97/theta/slacks; collection 121.1117.

## R10-R17 + A2R1-A2R2 (condensed for the byte cap; full text in validated base 941d6839; live mechanisms in STORAGE.md/source-core.md)

R10 (recover-Q 711418d3): fully priced linear algebra (loop-free row tail,
H4u), pad 1.61 -> 1.10 on the N/U trend. R11 (aa2a2ba5): exact logged price
quotients, 0.01 table shift, reserve 101/100 -> 1001/1000. R12/R12b
(rejected): unpriced/mispriced recovery, superseded by R12c (0f2da92):
recovery reruns the coset enumerator, P_MID 45.5. R13: SLACK 96 -> 60,
meter re-solved. R14 (v20tdiv): trial division to 2^30. R15 (v21):
first-pass records, rerun deleted. R16: boot/plan re-searched at slack 60,
BOOT16 (214,320,16,4), meter cap 2^118.3128, two-log 2^117.5873; stated
residuals corrected to computed figures, audit checks them. R17: merge
weight 100 -> 300, N/U 0.03225, row weight 304. A2R1 (from validated recP
42e80f80, re-derived): band-strip line sieve, trial division below 2^20.
A2R2: 7-owner LA row tail, C_ROW 1,320, C_ENTRY 20.7, blocks (22, 4).
Ledger 121.6221 -> 121.1556, G1 124.85609 -> 124.478970069061269.

## Evidence

`audit_g1.py` recomputes every changed quantity (`la_prices`, `matrix_pad`,
`exact_prices`, `price_shape`, `tab_shift_ranges`, `residual_bound`,
`residual_ledger`, `no_rerun`, `tdiv_prime_bound`, the BOOT16 pin;
proof.md, Evidence); `evidence16.py` [10] recounts the kernels.

## OR1 + SB1 (condensed for the byte cap; full text in validated 941d6839/7b6774d5/cf356f66; live mechanisms in source-core.md/model16.py)

OR1: tight count certs, theta 189/500 (1019 rows, margin 0.30040);
ledger 121.1242, G1 124.447570069061269. SB1: ratio work estimator,
pool cut (1061/2401), pad 1.08, P_MID 47; allowance 1.04408; ledger
121.0966, G1 124.420070069061269.

## LEAP1 (this package): f-first early-abort survivor ECM, adapted from validated recQ 2c2ce2ab

Each survivor runs exact norms, trial division below 2^20, then the f side
first (f band strip + 69-bit f ECM); the g-side band + 69-bit ECM runs only
if f factors fully, charged at q_f = min(1, 100 rho_f / P(omitted_f <=
slack_f)) (max 0.001624, mean 0.001303; H1 estimate with 100x safety,
capped at 1). The adaptation splits the A2R1 band strip by side; trial
division stays two-sided. Yield untouched (rels, margin, yield SE,
rel-means bit-identical to R3). Cofactorisation drops 44.3%; collection
120.7999 -> 120.7469 (333.4/cell, c_surv max 2^38.81); residuals
1.000015x/1.000012x; ledger -> 121.0820 (rows 121.0819329); G1 ->
124.405370069061269. LA, descent, pool, kernels, filtering untouched.

## STACK1 (this package): SB1 x LEAP1 stack on the validated LEAP1 bank

SB1's estimator (never fired standalone; dominated by LEAP1) applied on
the validated LEAP1 base 7b6774d5: disjoint mechanisms (work-estimator
variance + pad + sieve partition vs ECM schedule), re-measured jointly,
never added. Collection 120.7469 -> 120.7174 (326.7/cell, c_surv max
2^38.75); SE 3.493% -> 2.145%, allowance 1.07057 -> 1.04408; LA
117.0348 -> 116.9819, generator 106.8021 -> 106.7594; residuals
1.000016x/1.000013x; yield, q_f, descent, success untouched (two-log
2^117.5873, success model 0.3903324). Ledger -> **121.0554**
(121.0553890, rows 121.0553970); G1 -> **124.378770069061269**.

## STRAT1 (this package): conditional post-stratification of the work SE

Work SE conditioned on 15 fixed cells (5 qa bins x thr cuts 50/150;
min n 17): cond 1.697% -> 0.769%, total -> 1.521%, allowance ->
1.03042. Point/yield/q_f/descent/success untouched. Collection ->
120.6984 (322.4/cell); residuals 1.000016x/1.000013x. Ledger ->
**121.0404** (rows 121.0403861); G1 -> **124.363770069061269**.

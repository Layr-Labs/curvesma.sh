# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, filtering and H1-H5 are unchanged.

## R1-R9 (as in our AI-valid 3eccc6c9, 124.9584)

ECM_FIND 0.999 (403 curves); chunk-scatter route, measured 10.26 charged 10.3;
R3 replaced by R10; R5 withdrawn (re-read 7.003); [2^63, 2^64) record walk 2x;
line sieve from 2^11; R8 yield guard; R9 plan + LOG2_RELS 63.97, theta 190/500,
slacks (297, 233), allowance 1.0676; collection <= 121.1117, two-log 2^118.2732,
meter 2^118.996, success >= 0.390320.

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
`exact_prices`, `price_shape`, `tab_shift_ranges`, `width_pins`,
`residual_bound`, `residual_ledger`, `no_rerun`, `tdiv_prime_bound`, the BOOT16 pin;
proof.md, Evidence); `evidence16.py` [10] recounts the kernels.

## OR1/SB1/LEAP1 (condensed for the byte cap; full text in validated STACK1 cf356f66; live mechanisms in STORAGE.md/source-core.md)

OR1: tight count certificates (t least k/1000, per-region radii; 0.985729/1.013031), theta 189/500 (1019 rows), ledger 121.1242, G1 124.44757. SB1: fresh-sample ratio work estimator (5121/43222, seed 4), pool cut 191/500, pad 1.08, P_MID 47, allowance 1.04408, ledger 121.0966, G1 124.42007. LEAP1: f-first early-abort survivor ECM (q_f max 0.001624, 100x H1 safety), collection 120.7469, ledger 121.0820, G1 124.40537.

## STACK1 (validated base cf356f66): SB1 x LEAP1 stack on the validated LEAP1 bank

SB1's estimator (never fired standalone; dominated by LEAP1) applied on
the validated LEAP1 base 7b6774d5: disjoint mechanisms (work-estimator
variance + pad + sieve partition vs ECM schedule), re-measured jointly,
never added. Collection 120.7469 -> 120.7174 (326.7/cell, c_surv max
2^38.75); SE 3.493% -> 2.145%, allowance 1.07057 -> 1.04408; LA
117.0348 -> 116.9819, generator 106.8021 -> 106.7594; residuals
1.000016x/1.000013x; yield, q_f, descent, success untouched (two-log
2^117.5873, success model 0.3903324). Ledger -> **121.0554**
(121.0553890, rows 121.0553970); G1 -> **124.378770069061269**.

## WIDTH1 (this package): per-row certified bound widths for the survivor slack cover

Survivor bounds charged the global bisection width (8 bits): sf/sg =
P(omitted <= slack+8). Rows carry certified leaf widths far below 8
(median 3/2 bits f/g, max 8/8); the pool stores per-row max widths per
side (grid-16 rounded up, packed 9a+b, fail-closed 128), used by
`evaluate`; [2] recomputes them from (seed,index). Same bounds at
tighter valid inputs. Yield bit-identical (margin 0.30041); cond SE
1.697% -> 1.681%, need 1.04266 < 1.04408. Collection 120.7174 ->
120.6942 (321.4/cell); residuals 1.000016x/1.000013x; LA, descent,
success untouched. Ledger -> **121.0371** (rows 121.0370735); G1 ->
**124.360470069061269**.

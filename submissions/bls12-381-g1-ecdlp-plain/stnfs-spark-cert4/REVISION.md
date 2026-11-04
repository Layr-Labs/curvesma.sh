# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, filtering and H1-H5 are unchanged.

## R1-R9 (as in our AI-valid 3eccc6c9, 124.9584; condensed)

ECM_FIND 0.999; chunk-scatter route 10.3; R5 withdrawn; [2^63, 2^64) 2x;
line from 2^11; R8 guard; R9 plan + LOG2_RELS 63.97, slacks (297, 233).

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

## OR1/SB1/LEAP1 (condensed for the byte cap; full text in validated STACK1 cf356f66; mechanisms in STORAGE.md/source-core.md)

OR1: tight count certificate (t = least k/1000, per-side radii),
theta 189/500, G1 124.447570069061269. SB1: ratio work estimator
(fresh seed-4 rate), pool cut, pad 1.08, P_MID 47, G1 124.420070069061269
(never fired standalone). LEAP1 (7b6774d5): f-first early-abort ECM
(from recQ 2c2ce2ab), G1 124.405370069061269.

## STACK1 (validated base cf356f66): SB1 x LEAP1 stack

SB1's estimator on validated LEAP1, re-measured jointly. Collection
120.7469 -> 120.7174; SE 3.493% -> 2.145%, allowance 1.04408; LA
116.9819; yield, q_f, descent, success untouched. Ledger **121.0554**;
G1 **124.378770069061269**.

## CERT4 (this package): tight count certificates fund a 4-row theta shrink

Same half-count theorem (`lattice_count_bounds.py` byte-identical),
tighter valid inputs, re-measured jointly. A-side covers the pool ball
(R = 1024 exactly, the sampling domain; estimator scales by pool volume
per draw), not the region ball (R = 890). B-side per row at its own
slice radius (brad 1143..1420). Loosest products 0.986745/1.012094.
Theta 189/500 -> 24700/65536 (1019 -> 1015 rows, least boundary incl.
qa16 24699); fresh seed-4 sample regenerated (`canon12.b85`,
5072/42865 = 0.118325; same procedure, pool files unchanged). Yield
rels 64.130536, margin 0.300694, 1-SE 63.970051 (floors pass); means
745.920/499.035 (above 745.9/498.9 pins). SE 2.150% (need 1.042996,
allowance 1.04408). Collection -> 120.7035 (324.7/cell); LA, descent,
success bit-identical. Ledger -> **121.0445** (121.0443812, rows
121.0444113); G1 -> **124.367870069061269** (-0.010900 vs STACK1).

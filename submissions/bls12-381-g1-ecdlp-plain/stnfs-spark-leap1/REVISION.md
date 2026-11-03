# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, filtering and H1-H5 are unchanged.

## R1-R9 (as in our AI-valid 3eccc6c9, 124.9584)

ECM_FIND 0.999 (403 curves); chunk-scatter route, measured 10.26 charged 10.3;
R3 replaced by R10; R5 withdrawn (re-read 7.003); [2^63, 2^64) record walk 2x;
line sieve from 2^11; R8 yield guard; R9 plan + LOG2_RELS 63.97, theta 190/500,
slacks (297, 233), allowance 1.0676; collection <= 121.1117, two-log 2^118.2732,
meter 2^118.996, success >= 0.390320.

## R10-R15 (condensed for the byte cap; full text in base submission 941d6839; live mechanisms in STORAGE.md/source-core.md)

R10 (recover-Q 711418d3): fully priced linear algebra (loop-free row tail,
H4u), pad 1.61 -> 1.10 on the N/U trend. R11 (aa2a2ba5): exact logged price
quotients, 0.01 table shift, reserve 101/100 -> 1001/1000. R12/R12b
(rejected): unpriced/mispriced recovery, superseded by R12c (0f2da92):
recovery reruns the coset enumerator, P_MID 45.5. R13: SLACK 96 -> 60,
meter re-solved. R14 (v20tdiv): trial division to 2^30. R15 (v21):
first-pass records, rerun deleted. Ledger 121.6221 -> 121.3644, G1
124.85609 -> 124.687770069062 across R10-R15.

## Evidence

`audit_g1.py` recomputes every changed quantity (`la_prices`, `matrix_pad`,
`exact_prices`, `price_shape`, `tab_shift_ranges`, `residual_bound`,
`residual_ledger`, `no_rerun`, `tdiv_prime_bound`, the BOOT16 pin;
proof.md, Evidence); `evidence16.py` [10] recounts the kernels.

## R16 (this package): descent boot + plan re-searched at slack 60, meter re-solved; stated residual multipliers corrected

The shipped boot/plan was coordinate-searched before R13 moved SLACK 96 ->
60. Re-searching at 60 under the same model and constraints (mu >= 16, phi
<= 0.05, psi <= 0.8, iA on the x8 mesh lattice, R in RADII) moves BOOT16
(215,344,15,4) -> **(214,320,16,4)** (nA 189.83 hugs the 190 bin top, so the
16 A-children stay in cheaper bins) and twelve plan entries (85f, 85g, 95f,
110f, 110g, 115f, 120g, 125f, 125g, 165f, 165g, 185g). No model, heuristic,
price or collection parameter moves. Expected two-log work 2^118.2303 ->
**2^117.5873** (one log 2^130.875; boot psi 0.6732, phi 1.780e-03; min mu
32.3; E[N] 2^16.87, P(abort) 2^-10.42); cap re-solved to **2^118.3128**
(Markov 0.60478), re-binding success at the 0.390320 claim (model 0.3903323).
Ledger <= 121.3644 -> **<= 121.2637** (unrounded 121.2636532, rows
121.2636860); G1 -> **124.587070069062**.

Also fixed: v21 was rejected only for stating the v19 residual figure
(1.000014x) while the computed collection multiplier is 1.0000286x (ledger
1.0000224x); the verdict raised no objection to any model number or
mechanism. This package states the computed figures (collection **1.000029x**,
ledger **1.000023x**) and the audit now checks the stated multipliers, not
just 1.0001.

## R17 (this package): merge weight 100 -> 300, blocks (26, 5)

W-scan 100..600 (2^15, 2^17): N/U 0.03225 at W = 300; H4u flat in W
(0.889+ from 2^15 up, C_ENTRY = 21 stands); maxcoef <= 81. Tables redone
at W300 (maxima 0.06519/0.04310/0.03675/0.03480/0.03335/0.03225/0.03038;
charged 0.03547 covers 2^15+); row weight 304 (max 302.4). Rest untouched.
LA -> **2^117.2095**, generator -> **2^107.0654**, ledger -> **<= 121.2253**
(121.2251689, rows 121.2252062); G1 -> **124.548670069062**. Stated ledger
residual -> **1.000024x** (computed 1.0000231x).

## Angle-2 R1 (this package): band-strip 2^48 line sieve

Line hits emit no records; survivors get trial division below 2^20 plus
two rounds per side of the 48-bit ecm16 schedule at 99.9% (STORAGE.md;
method from validated recP 42e80f80, re-derived). Yield x band_keep per
cell (loss 9.2e-05 bits; input 63.9731 >= 63.97). Collection 120.8392
(346.0/cell); allowance 1.06985 (1 + 2 x 0.03492368); ledger 121.1663
(rows 121.1662337); G1 124.489670069061269. Descent/LA/pool untouched;
residuals hold (1.0000145x/1.0000115x). ECM_TAB (dead) dropped.

## Angle-2 R2 (this package): 7-owner LA row tail, blocks (22, 4)

N + 21 < 7 x 2^56, so only 7 block owners are active (audit `la_owners`):
C_ROW 2,750 -> **1,320** (counted 1,317), C_ENTRY 21 -> **20.7** (needs
20.64 at H4u 0.88), BW (26, 5) -> **(22, 4)** (wide-scan minimum).
LA 117.0348, generator 106.8021; ledger -> **121.1556** (121.1555425,
rows 121.1555842); G1 -> **124.478970069061269**. All else A2R1.

## OR1 (this package): tight lattice-count certificate, theta 189/500

`lattice_count_bounds.py` bounded the unit-cube Q-shift sqrt(sum|Q_ij|)/2
(1.044 Q_a, 1.125 Q_b) by 2 and gave Q_b slices the a-radius; now t is the
least k/1000 above it (261/250, 563/500) and Q_b uses isqrt((1-theta) A^2).
Same theorem: products 0.971858/1.027240 -> 0.985729/1.013031. Region
re-solved under the frozen descent floor: 189/500 (1019 rows), 1-SE count
2^63.97008, margin 0.30040, means 745.95/499.03; allowance 1.07057.
Collection 120.7999; ledger 121.1242 (rows 121.1241083); G1
124.447570069061269. Descent, pool, kernels, filtering, LA untouched.

## LEAP1 (this package): f-first early-abort survivor ECM, adapted from validated recQ 2c2ce2ab

Each survivor runs exact norms, trial division below 2^20, then the f side
first (f band strip + 69-bit f ECM); the g-side band + 69-bit ECM runs only
if f factors fully, charged at q_f = min(1, 100 rho_f / P(omitted_f <=
slack_f)) (max 0.001624, mean 0.001303 over rows; H1 estimate with 100x
safety, capped at 1). The adaptation to plain splits the A2R1 band strip
by side (recQ has no band); trial division stays two-sided (conservative,
as validated). Yield is untouched (a failed f side yields no relation):
rels, margin, yield SE, rel-means bit-identical to R3; allowance 1.07057
still covers 1 + 2 SE (SE 3.493% vs 3.528%). Cofactorisation drops 44.3%;
collection 120.7999 -> **120.7469** (333.4/cell, c_surv max 2^38.81);
residuals 1.000015x/1.000012x; ledger -> **121.0820** (121.0819193, rows
121.0819329); G1 -> **124.405370069061269**. LA, descent, pool, kernels,
filtering, field map untouched.

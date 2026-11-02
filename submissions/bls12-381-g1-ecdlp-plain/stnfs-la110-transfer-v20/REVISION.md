# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, filtering and H1-H5 are unchanged.

## R1-R9 (as in our AI-valid 3eccc6c9, 124.9584)

ECM_FIND 0.999 (403 curves); chunk-scatter route, measured 10.26 charged 10.3;
R3 replaced by R10; R5 withdrawn (re-read 7.003); [2^63, 2^64) record walk 2x;
line sieve from 2^11; R8 yield guard; R9 plan + LOG2_RELS 63.97, theta 190/500,
slacks (297, 233), allowance 1.0676; collection <= 121.1117, two-log 2^118.2732,
meter 2^118.996, success >= 0.390320.

## R10 (from our recover-Q 711418d3, AI-valid): fully priced linear algebra, pad 1.10

Loop-free kernels (row tail 2,750/row, red9 420, 21/entry, H4u); pad 1.61 ->
**1.10** on the measured N/U trend (H3). LA 2^118.7428 -> **2^117.7327**,
generator 2^105.5088 -> 2^107.2082; ledger <= 121.6221 -> **<= 121.5198**.

## R11 (submitted as aa2a2ba5, failed on packaging oversize only): exact price ceilings, 0.01 table shift, reserve 1001/1000

Nine line/record prices become the exact logged quotients (`exact_prices`);
tables checked decreasing and convex (`price_shape`); `tab()` charges 0.01
below its argument with log2h coverage checked (`tab_shift_ranges`); half-ulp
residual moves the ledger 1.000036x (`residual_bound`), so the reserve drops
101/100 -> **1001/1000**. Collection <= 121.1117 -> **<= 121.1127**; ledger <=
121.5198 -> **<= 121.5206**; G1 124.85609 -> **124.843970069062**.

## R12/R12b (494e4872 and 43f5ea7e, both rejected; superseded by R12c)

Moved P_MID to 48.75/47.5 with deficient recovery pricing (unpriced line
factors; then an unjustified step-p loop). R12c below completes the revision.

## R13 (submitted as 16253187, validation cancelled: it inherited the R12b gap; content re-verified here)

SLACK 96 -> **60** in T. Tighter T cuts node work faster than retention loss
rises (ret min 0.8575; node phi <= 0.0386, psi <= 0.8): two-log 2^118.2732 ->
**2^118.2303**. Meter re-solved: cap **2^118.9532**, Markov **0.60590**, boot
phi **7.858e-04**, E[N] 2^18.005, success claim **0.390320**.

## R12c (as in our AI-valid 0f2da92, 124.754970069062): recovery reruns the coset enumerator, P_MID 45.5

The second pass reruns the first-pass coset enumerator over every (line, p,
pair) and flag-checks each hit (STORAGE.md): identical stepping code, so
stepping costs at most the measured update `tab(h)`; flag ops LDRB+CBZ = 2,
recounted (`line_recov_loop`); emit <= 2^97 (absorbed); setups recompute the
same structure (488.24 per (line, p, pair)). Completeness: slice12's own
line-path-vs-brute-force checks plus a packed-layout toy (`line_recov_toy`).
P_MID **45.5** (0.25-grid optimum; work SE 3.406%, allowance holds).
Collection <= 121.1127 -> **<= 121.0039**; ledger <= 121.5206 -> **<= 121.4316**;
G1 124.843970069062 -> **124.754970069062**.

## R14 (this package): unit-duplicate factor 0.25 -> 0.45, meter re-solved

UDUP 0.25 -> **0.45** (audit pins it): 79 production special-q lattices give
distinct/raw 0.574-0.906 (mean 0.78; bound-4 agrees to 0.018, second seeds to
0.061), so 0.45 is the survey minimum / 1.28 (DESCENT.md). Higher mu cuts boot
trials and retries (min mu 60.4; max node phi 0.0263, psi <= 0.8): two-log
2^118.2303 -> **2^117.2223** (one log 2^130.510). Meter re-solved: cap
**2^117.9395**, Markov **0.60828**, boot phi **6.324e-06**, E[N] 2^16.965,
P(abort) 2^-10.33, success claim **0.390320**. Ledger **<= 121.2946**
(unrounded 121.2945400, rows 121.2945509); G1 124.754970069062 ->
**124.617970069062**. Nothing else moves (H1-H5/H4u).

## Evidence

`audit_g1.py` recomputes every changed quantity (`la_prices`, `matrix_pad`,
`exact_prices`, `price_shape`, `tab_shift_ranges`, `residual_bound`,
`line_recov_loop`, `line_recov_nlp`, `line_recov_toy`; proof.md, Evidence);
`evidence16.py` [10] recounts the kernels.

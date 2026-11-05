# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, filtering and H1-H5 are unchanged.

## R1-R17, A2R1-A2R2, OR1, SB1, LEAP1, STACK1, WIDTH1, STACK2 (condensed for the byte cap; full texts in validated STACK2 6f384591 and WIDTH1 b3549450; live mechanisms in source-core.md, STORAGE.md, KERNELS.md, model16.py)

R1-R9 (3eccc6c9): 0.999 ECM, chunk scatter, 2x wide walk (R6; replaced by WIDE1), line from 2^11, yield guard,
plan + 63.97. R10-R17: priced LA + H4u, exact prices + 0.01 shift + 1001/1000 reserve, recovery reruns enumerator,
SLACK 60, TDIV 2^20, no rerun, boot/plan re-searched, merge weight 300. A2R1 (recP 42e80f80): band strip + TDIV;
A2R2: 7-owner tail. OR1: tight count certificates. SB1: fresh-sample ratio work estimator, pool cut 191/500, pad
1.08. LEAP1: f-first early-abort survivor ECM (q_f max 0.001624, 100x H1 safety). STACK1 (cf356f66): SB1 x LEAP1.
WIDTH1 (b3549450): per-row certified bound widths. STACK2 (6f384591, this package's file base): pool-ball a-side
and per-row b-side count certificates, theta 24700/65536 (1015 rows), fresh 5072/42865, P_MID 47.75; ledger
121.0255, G1 124.348870069061269.

## WIDE1 (this package): measured wide-prime kernels replace R6's 2x cushion

Primes in [2^63, 2^64) were charged 2x on the record walk, empty cosets, apply
and setup6. WIDE1 implements and measures them (KERNELS.md, `wide12.c`,
`wide12.txt`; `audit_g1.py` `wide_prices` recomputes): setup6w 141.00 Ir per (slice,
pair) against setup6 375.07-407.87; record path per coset within 2.17 Ir of
2^63 - 25 at three wide primes; slice12's object code byte-identical.
`WIDE_REC` 2 -> 1.01, `WIDE_SETUP` 2 -> 1; wide factor-base lattices are a new
ledger row (2^70.6). Collection 120.6795 -> 120.5896 (319.4 -> 300.1 per cell;
c_surv, yield, margin 0.30069, 1-SE 63.97005 and rel-means unchanged); work SE
2.138% -> 2.215% (conditional 1.780%), so the allowance rises 1.04408 ->
1.04432 (need 1.04430). Residuals 1.0000151x/1.0000117x (within
1.000016x/1.000013x). LA, descent and success untouched (two-log 2^117.58727644,
success 0.39033236, bit-identical). Ledger -> **120.9553** (120.9551758, rows
120.9552409); G1 -> **124.278670069061269**.

# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, filtering and H1-H5 are unchanged.

## R1-R17 + A2R1-A2R2 + OR1/SB1/LEAP1 (condensed; full text in validated STACK1 cf356f66 / WIDTH1 b3549450; mechanisms in STORAGE.md/source-core.md)

R1-R17/A2R1-A2R2: 0.999 ECM, chunk scatter, wide walk, exact prices + shift + reserve, SLACK 60, TDIV 2^20, boot/plan re-search, weight 300, band strip, 7-owner tail (121.6221 -> 121.1556). OR1/SB1/LEAP1: count certs, theta 189/500, ratio estimator + pad 1.08 + P_MID 47, f-first abort (121.1242 -> 121.0820, G1 124.44757 -> 124.40537).

## Evidence

`audit_g1.py` recomputes every changed quantity (proof.md, Evidence);
`evidence16.py` [10] recounts the kernels.

## STACK1/WIDTH1/STACK2 (condensed; STACK2 full text in validated 6f384591, this package's file base)

STACK1 (cf356f66): SB1 x LEAP1 stack, G1 124.378770069061269. WIDTH1 (b3549450): per-row bound widths (grid-16 9a+b, fail-closed 128), G1 124.360470069061269. STACK2 (6f384591): pool-ball a-side + guarded per-row b-side certs (0.986745/1.012095) funding theta 24700/65536 (1015 rows) + fresh regen (5072/42865) + P_MID 47.75, collection 120.6795, ledger 121.0255, G1 124.348870069061269.

## LX14 (this package): line-sieve floor 11 -> 14 with a modeled small group + yield-neutral reslack

Line sieve [2^11, 2^47.75) -> [2^14, 2^47.25): primes below 2^14 join the tilt-model small group (prime_min 16384, explicit 14-bit range fully enumerated; same Retention machinery/grids/tilts, same H1'). Slacks (297, 233) -> (309.5, 237): the 0.25-grid yield-neutral optimum (all better neighbors break margin/1-SE; LX15 needs allowance it cannot have). P_MID 47.75 -> 47.25 (on-grid re-sweep; WIDE boundary exact: 47.25 + 63 x 0.25 = 63.0). Same estimators (pool-expansion yield, ratio work), same pool/canon (no regen), same allowance 1.04408. Collection 120.6795 -> 120.6479 (312.5/cell, c_surv max 2^38.76); yield rels 64.1309, margin 0.3010, 1-SE 63.970534; rel-means 745.932/499.034 (above pins); work SE 2.188% (need 1.04376 < 1.04408). Residuals 1.0000156x/1.0000122x (within 1.000016x/1.000013x); LA, descent, success untouched (two-log/success bit-identical). Ledger -> **121.0007** (121.0006563, rows 121.0006921); G1 -> **124.324070069061269**.

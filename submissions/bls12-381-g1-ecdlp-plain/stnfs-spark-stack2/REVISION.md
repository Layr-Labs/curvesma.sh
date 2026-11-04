# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, filtering and H1-H5 are unchanged.

## R1-R17 + A2R1-A2R2 (condensed for the byte cap; full text in validated WIDTH1 b3549450; live mechanisms in STORAGE.md/source-core.md)

R1-R9 (3eccc6c9): 0.999 ECM, chunk scatter, 2x wide walk, line from 2^11,
yield guard, plan + 63.97. R10-R17: priced LA + H4u, exact prices + 0.01
shift + 1001/1000 reserve, recovery reruns enumerator, SLACK 60, TDIV 2^20,
no rerun, boot/plan re-searched, merge weight 300. A2R1 (recP 42e80f80):
band strip + TDIV; A2R2: 7-owner tail. Ledger 121.6221 -> 121.1556, G1
124.85609 -> 124.478970069061269.

## Evidence

`audit_g1.py` recomputes every changed quantity (proof.md, Evidence);
`evidence16.py` [10] recounts the kernels.

## OR1/SB1/LEAP1 (condensed for the byte cap; full text in validated STACK1 cf356f66; live mechanisms in STORAGE.md/source-core.md)

OR1: tight count certificates (t least k/1000, per-region radii; 0.985729/1.013031), theta 189/500 (1019 rows), ledger 121.1242, G1 124.44757. SB1: fresh-sample ratio work estimator (5121/43222, seed 4), pool cut 191/500, pad 1.08, P_MID 47, allowance 1.04408, ledger 121.0966, G1 124.42007. LEAP1: f-first early-abort survivor ECM (q_f max 0.001624, 100x H1 safety), collection 120.7469, ledger 121.0820, G1 124.40537.

## STACK1/WIDTH1 (condensed for the byte cap; full text in validated WIDTH1 b3549450; live mechanisms in STORAGE.md/source-core.md/model16.py)

STACK1 (validated cf356f66): SB1 x LEAP1 stack (estimator + pad 1.08 + P_MID 47 on LEAP1), ledger 121.0554, G1 124.378770069061269. WIDTH1 (validated b3549450, this package's file base): per-row certified bound widths for the survivor slack cover (grid-16 packed 9a+b, fail-closed 128; [2] regen-checked), collection 120.7174 -> 120.6942, ledger 121.0371, G1 124.360470069061269. Yield bit-identical across both (margin 0.30041).

## STACK2 (this package): tight count certificates fund a 4-row shrink on WIDTH1

Same half-count theorem (`lattice_count_bounds.py` byte-identical),
tighter valid inputs: (1) a-side covers the pool ball (R = 1024, the
sampling domain); (2) per-row b-side slices at guarded brad 1143..1420
(-1 ulp guard; loosest products 0.986745/1.012095). (3) Theta 189/500 ->
24700/65536 (1019 -> 1015 rows, least boundary). (4) Fresh seed-4 sample
regenerated at the new theta (5072/42865; procedure reproduces the base
file byte-identically). (5) P_MID 47 -> 47.75 (on-grid 0.25-sweep optimum;
WIDE boundary exact on grid). Collection 120.6942 -> 120.6795 (319.4/cell,
c_surv max 2^38.79); yield rels 64.1305, margin 0.3007, 1-SE 63.970051;
rel-means 745.920/499.035 (above pins); work SE 2.138% (need 1.04277 <
1.04408). Residuals 1.0000147x/1.0000116x (within 1.000016x/1.000013x);
LA, descent, success untouched (two-log/success bit-identical).
Ledger -> **121.0255** (121.0254543, rows 121.0254955); G1 ->
**124.348870069061269**.

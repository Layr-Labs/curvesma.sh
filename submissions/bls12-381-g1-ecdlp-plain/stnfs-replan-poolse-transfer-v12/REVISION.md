# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, kernel prices, filtering, descent plan and H1-H5
are unchanged.

## R1-R7

- **R1** (our 27202622): ECM_FIND 0.99 -> 0.999 above X2 (403 curves) buys a
  smaller region; 2^12 filtering rerun N/U 0.0715 < charged 0.0758.
- **R2** (recover-P 43d603a4/be446c45): records written at emission to their
  chunk's page through a 2^29-entry tail table (`store16.scatter` ==
  `store16.route`, `route_hot_path`): 36 -> 21 per record.
- **R3** (here, with recover-Q 06aac4a1's row tail): 88 -> **44** per sparse
  entry (STORAGE.md: 11-instruction +-1 loop, unit share charged 0.88 under
  H4u; row tail 2,090), block Wiedemann (18, 4): linear algebra **2^118.74271**,
  generator 2^105.5088.
- **R4** (recover-P thread): measured route 7 + 13/4 + < 64/2^16 = 10.26 per
  record, charged **10.3** (recover-P measured 10.0009 under qemu-aarch64).
- **R5** withdrawn: every record is re-read at 7.003, as in the source.
- **R6** (recover-Q thread): primes in [2^63, 2^64) move from ECM to the record
  walk; measured prices stop below 2^63, so their record-walk hits, empty
  cosets, apply and setup6 are charged **2x** (WIDE_REC, WIDE_SETUP; KERNELS.md).
- **R7** (recover-Q thread): line sieve from 2^11 (`P_MIN` 2048, `LX0` 11,
  config `prime_min`; `threshold_model` counts p < prime_min against the slacks).

## R8 (this package): relation input to the descent at a yield lower bound

The descent's missing-log input LOG2_RELS must hold at a linear 2-SE lower
confidence bound on the yield: rels x (1 - 2 x relative pool SE) >= 2^LOG2_RELS
(`model16.YIELD_SE_K`, `DESCENT_RELS_MIN`; audit
`yield_lower_bound_covers_descent`). The relative SE is that of the fixed-pool
estimator over all N = 20,000 draws (zero outside the canonical primitive
sub-region): (Q/S^2 - 1/N) N/(N - 1) (`N_POOL`; was 1/n_eff). A guard continues collection until
2^LOG2_RELS relations exist.

## R9 (recover-Q thread): descent plan for fewer relations, joint optimum

Six _PLAN entries (old iA units): 90g (75,42) -> (75,43), 95g (80,24) ->
(80,25); f-side n at 100/110/120/150: 83/85/94/110 -> 84/86/95/109. Lowering
LOG2_RELS trades descent work for a smaller region; joint optimum over
63.84-64.13 and 1,700 regions: **LOG2_RELS 63.90**, theta
197/500 (1187 canonical rows, a-radius >= 908, count products 0.972347748 /
1.026754289), slacks (299, 233). Relations 2^64.2178, pool SE 0.1417 bits,
lower bound 2^63.9024, margin 0.3880. Per cell 419.6 over 2^126.811 cells;
work SE 3.171%, allowance 1.0635; collection 2^121.1470 x 1.0635 <=
2^121.2358. Descent at norms 746.8/499.0: two-log 2^118.58258, boot
2.3176e-03, E[N] 2^18.339, counter 2^-8.96; meter **2^119.314** (Markov
0.60232), success >= 0.390430.

| core ledger row | R1-R4 | + R6, R7 | + R6-R8 | + R6-R9 |
|---|---:|---:|---:|---:|
| relation collection | <= 121.2551 | <= 121.1117 | <= 121.5599 | <= 121.2358 |
| two-log meter | 118.635 | 118.635 | 118.627 | 119.314 |
| core ledger | <= 121.6752 | <= 121.5694 | <= 121.908 | **<= 121.7635** (rows sum 121.763451; unrounded 121.7634302) |

Linear algebra 118.7428, generator 105.5088 and fixed rows unchanged.

## Evidence

`audit_g1.py` recomputes every changed quantity (proof.md, Evidence).

# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, kernel prices, filtering and H1-H5 are unchanged.

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

The descent's missing-log input LOG2_RELS must hold at a linear 1-SE lower
confidence bound on the yield: rels x (1 - relative pool SE) >= 2^LOG2_RELS
(`model16.YIELD_SE_K` = 1, `DESCENT_RELS_MIN`; audit
`yield_lower_bound_covers_descent`). The relative SE is that of the fixed-pool
estimator over all N = 20,000 draws: (Q/S^2 - 1/N) N/(N - 1) (`N_POOL`). A guard continues collection
until 2^LOG2_RELS relations exist; the charged region covers this unless the
yield is more than 1 SE below its estimate. 2-SE versions: c53a007e (LOG2_RELS
64.131, core 2^121.908) and 30d365cb (63.90, core 2^121.7635), both AI-valid.

## R9 (from the recover-Q thread): descent plan for fewer relations

Six plan entries re-optimised for a lower relation input (old iA units): 90g
(75,42) -> (75,43), 95g (80,24) -> (80,25), and f-side n 100: 83 -> 84,
110: 85 -> 86, 120: 94 -> 95, 150: 110 -> 109; boot unchanged. With it,
lowering LOG2_RELS trades descent work for a smaller region; jointly optimised
over LOG2_RELS 63.84-64.13 and ~480 regions: **LOG2_RELS 63.97**, theta
190/500, slacks (297, 233); the 0.30-bit margin binds. Relations
2^64.1311, pool SE 0.1495 bits, lower bound 2^63.9733, margin 0.3012. Per cell
417.9 over 2^126.692 cells; work SE 3.377%, allowance 1.0676; collection
2^121.0173 x 1.0676 <= 2^121.1117. Descent at norms 745.9/498.9: two-log
2^118.27320, boot 7.746e-04, E[N] 2^18.0001, counter 2^-9.29; meter
**2^118.996** (Markov 0.60593), success >= 0.390320.

| core ledger row | R1-R4 | + R6-R8 (2 SE) | + R6-R9 (2 SE) | + R6-R9 (1 SE) |
|---|---:|---:|---:|---:|
| relation collection | <= 121.2551 | <= 121.5599 | <= 121.2358 | <= 121.1117 |
| two-log meter | 118.635 | 118.627 | 119.314 | 118.996 |
| core ledger | <= 121.6752 | <= 121.908 | <= 121.7635 | **<= 121.6221** (rows sum 121.622001; unrounded 121.6219553) |

Linear algebra 118.7428, generator 105.5088 and fixed rows unchanged.

## Evidence

`audit_g1.py` recomputes every changed quantity (proof.md, Evidence).

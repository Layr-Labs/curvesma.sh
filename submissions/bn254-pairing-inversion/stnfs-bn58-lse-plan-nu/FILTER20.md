# Filtering ratio at 2^20 orbit columns

x = N / U' after purge, clique removal and Markowitz merge to weight 400, profile of proof.md section 4. Seeds 1..8 per scale, fixed before the first run, none dropped.

| U' | source | seeds | mean | SE | mean + 2 SE |
|---|---|---|---|---|---|
| 2^17 | package | 2 | 0.01845 | 0.00143 | 0.021313 |
| 2^17 | independent | 8 | 0.01720 | 0.00042 | 0.018043 |
| 2^18 | package | 8 | 0.01624 | 0.00022 | 0.016682 |
| 2^18 | independent | 8 | 0.01617 | 0.00016 | 0.016495 |
| 2^19 | package | 8 | 0.01518 | 0.00013 | 0.015438 |
| 2^19 | independent | 8 | 0.01500 | 0.00011 | 0.015224 |
| 2^20 | independent | 8 | 0.01444 | 0.00015 | 0.014739 |

The 2^20 rows do NOT come from this package's filtsim.py. They come from an independent implementation of the same simulation (same occurrence model, purge, clique and merge rules, 255 bounds; other generator and tie-breaking, so comparable in distribution only). Validation on the package's three scales: its mean differs by -0.84, -0.27 and -1.07 combined SE (2^17, 2^18, 2^19), slightly lower each time. So the charge does not use the independent 2^19 value: by the unchanged rule NU = 0.01544 >= 0.015438 is the package's own 2^19 value, and the 2^20 value 0.014739 only has to be below it. One 2^20 run refused pivots at the multiplier bound, which can only raise N. The simulator is not shipped (size); filtbn.b85 holds its per-seed rows (FILT_RUNS 20, IND_RUNS 17-19) and check [6] recomputes means, SEs and the charged ratio. Production is 32 bits above 2^20 (H3).

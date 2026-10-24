| d | trials | success_raw | CP-95 raw | merge | CP-95 merge | class cov | PAIR cov | prev-matched | escapes |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 96000 | 0.439531 | 0.436893 | 0.413490 | 0.410872 | 0.4135 | 0.3920 | 0.9326 | 11690 |
| 7 | 96000 | 0.510781 | 0.508122 | 0.400313 | 0.397708 | 0.4003 | 0.3939 | 0.9779 | 14154 |
| 9 | 96000 | 0.448177 | 0.445533 | 0.395281 | 0.392683 | 0.3953 | 0.3926 | 0.9886 | 18959 |
| 11 | 48000 | 0.189562 | 0.186623 | 0.391271 | 0.387601 | 0.3913 | 0.3894 | 0.9912 | 11015 |

| d | records/trial | spacings/walker | restarts/walker | trap rate | 1/(6*R_eff) | ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 781.26 | 11.656 | 11.656 | 4.87e-06 | 5.09e-06 | 0.958 |
| 7 | 228.63 | 2.914 | 2.914 | 5.04e-06 | 5.09e-06 | 0.991 |
| 9 | 76.45 | 0.729 | 0.729 | 5.05e-06 | 5.09e-06 | 0.992 |
| 11 | 22.01 | 0.182 | 0.182 | 5.09e-06 | 5.09e-06 | 1.000 |

| quantity | value | source |
|---|---:|---|
| log2 S0 | 126.12850516 | ledger.py |
| log2 W | 126.26594263 | ledger.py |
| claim | 126.266 | ledger.py |
| step instructions | 167 | ledger.py ROWS |
| model success | 0.399014 | ledger.SUCCESS |
| R_EFF / TRAP_RATE | 852105 / 1.956e-07 | ledger.py |
| min CP-95 raw | 0.186623 (d=11) | census |
| min CP-95 merge, gated | 0.392683 (d=9) | census |
| medium: fires / K_obs / residual | 2581 / 1 / 0 | medium fixture |
| medium trapped steps: control / escape | 79.66% / 0.00% | fixture |
| equivariance mismatches (fixed / control) | 0 / 14988 | medium fixture |
| C fixture: classes / spacing / fires | 2^26.00 / 1048576 / 18971 | C fixture |
| C fixture: trap ratio / invariant bad | 1.160 / 0 | C fixture |

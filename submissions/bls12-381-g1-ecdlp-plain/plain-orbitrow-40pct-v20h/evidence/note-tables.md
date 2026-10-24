| d | trials | success_raw | CP-95 raw | merge | CP-95 merge | class cov | PAIR cov | PAIR CP-95 | escapes |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 96000 | 0.439531 | 0.436893 | 0.413490 | 0.410872 | 0.4135 | 0.3920 | 0.3894 | 11690 |
| 7 | 96000 | 0.510781 | 0.508122 | 0.400313 | 0.397708 | 0.4003 | 0.3939 | 0.3913 | 14154 |
| 9 | 96000 | 0.448177 | 0.445533 | 0.395281 | 0.392683 | 0.3953 | 0.3926 | 0.3900 | 18959 |
| 11 | 48000 | 0.189562 | 0.186623 | 0.391271 | 0.387601 | 0.3913 | 0.3894 | 0.3857 | 11015 |

| d | records/trial | in-budget | drain | spacings/walker | restarts/trial |
|---:|---:|---:|---:|---:|---:|
| 5 | 781.26 | 745.89 | 35.37 | 11.656 | 745.89 |
| 7 | 228.63 | 186.59 | 42.04 | 2.914 | 186.59 |
| 9 | 76.45 | 46.65 | 29.80 | 0.729 | 46.65 |
| 11 | 22.01 | 11.64 | 10.37 | 0.182 | 11.64 |

| quantity | value | source |
|---|---:|---|
| log2 W | 126.26594263 | ledger.py |
| claim / gate | 126.266 / 126.274 | ledger.py |
| step instructions | 167 | ledger.py ROWS |
| R_EFF / TRAP_RATE | 852105 / 1.956e-07 | ledger.py (selector.reachable_columns) |
| three-cycle charge | 1.377e-12 /step, 4.23e+26 ref | ledger.TRAP3_RATE |
| min CP-95 merge, gated | 0.392683 (d=9) | census |
| medium: fires / K_obs / residual | 2581 / 1 / 0 | medium fixture |
| medium trapped: control / escape | 79.66% / 0.00% | fixture |
| equivariance (comparisons / attributions) | 9000 / 14988 | medium fixture |
| twist_index_production agreement / mutation | 4000 of 4000 / 358 of 1000 | twist-index.json |
| K_ESC enforced: cap / max run / tripped | 3 / 1 / no | escape-kcap.json |
| three-cycle structural / visits | 0 / 31847 | cycle3.json |
| C fixture fires / steps | 18971 / 50331648 | C fixture |
| self-audit | PASS (evidence/selfaudit.json) | scripts/selfaudit.py |

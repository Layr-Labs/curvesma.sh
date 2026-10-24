| d | trials | success_raw | CP-95 raw | merge | CP-95 merge | raw-merge | drain-only pairs | escapes |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 48000 | 0.433979 | 0.430251 | 0.407458 | 0.403763 | +0.0265 | 1604 | 2994 |
| 7 | 48000 | 0.512271 | 0.508507 | 0.401250 | 0.397564 | +0.1110 | 6067 | 3517 |
| 9 | 48000 | 0.451688 | 0.447943 | 0.398687 | 0.395006 | +0.0530 | 7260 | 4743 |
| 11 | 48000 | 0.187771 | 0.184842 | 0.391479 | 0.387809 | -0.2037 | 2864 | 5490 |

| d | pairs (budget + drain) | records (budget + drain) | drain window / budget | trap rate |
|---:|---:|---:|---:|---:|
| 5 | 21221 + 3845 | 746.0 + 35.7 | 0.172 | 2.49e-06 |
| 7 | 13310 + 15871 | 186.6 + 42.0 | 0.686 | 2.51e-06 |
| 9 | 5035 + 21069 | 46.6 + 29.8 | 1.000 | 2.52e-06 |
| 11 | 1352 + 8344 | 11.6 + 10.3 | 1.000 | 2.54e-06 |

| quantity | value | source |
|---|---:|---|
| log2 S0 | 126.12850516 | ledger.py |
| log2 W | 126.26594270 | ledger.py |
| claim | 126.266 | ledger.py |
| step instructions | 167 | ledger.py ROWS |
| model success | 0.399014 | ledger.SUCCESS |
| min CP-95 raw | 0.184842 (d=11) | census |
| min CP-95 merge | 0.387809 (d=11) | census |
| escape fires / K_esc | 734 / 4 | medium fixture |
| equivariance mismatches | 0 of 12000 | medium fixture |

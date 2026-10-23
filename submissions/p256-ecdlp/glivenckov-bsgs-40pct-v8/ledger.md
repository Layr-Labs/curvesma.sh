# Exact generated v8 ledger

| Quantity | Exact value |
|---|---:|
| positive odd babies m | 170141183460469231731687303715884105728 |
| L=2m=4^64 | 340282366920938463463374607431768211456 |
| actual h=ceil(r/L) | 340282366841710300967557013911933812736 |
| actual records m+h | 510423550302179532699244317627817918464 |
| charged records m+L=3*2^127 | 510423550381407695195061911147652317184 |
| M raw | 637 |
| M ceiling | 765 |
| A raw | 40 |
| A ceiling | 48 |
| J=16M+32A+512 | 14288 |
| I=512(M+64)+256 | 424704 |
| ceil(I/1024), separate | 415 |
| generation: complete_case_Jacobian_update | 14288 |
| generation: batch_prefix_backward_3M | 2295 |
| generation: separate_ceil_I_over_1024 | 415 |
| generation: x_scaling_2M | 1530 |
| generation: y_scaling_3M | 2295 |
| generation: canonical_Montgomery_decoding_2M | 1530 |
| generation: batch_array_masks_control | 256 |
| generation: record_build_and_local_send_preparation | 512 |
| K_rec exact sum | 23121 |
| baby extra lookup/label | 64 |
| V=512J+I+7M+512 | 7746027 |
| expanded table/seed cost | 26066093 |
| charged table/seed envelope | 32089836 |
| baby / giant hops | 128 / 129 |
| charged nodes / END edges | 1020847100762815390390123822295304634368 |
| R_setup | 17179869184 |
| R_small_table_and_seeds | 32089836 |
| R_init | 16725558898897967356151788704486271129485312 |
| R_gen | 11801502908368527320605026447644869225611264 |
| R_baby_lookup_labels | 10889035741470030830827987437816582766592 |
| R_place | 8406335592414863801399206301994401895809024 |
| R_sort | 34148016085250016685476568604992803556032512 |
| R_merge | 261336857795280739939871698507597986398208 |
| R_verify | 16540630 |
| T exact instruction sum | 71353639378468125934403289745063777604602562 |
| M-only recurring part | 10152324417086199057429781412726804588789760 |
| remaining instructions | 61201314961381926876973508332336973015812802 |
| T_stress=Mpart+6*rest/5 | 417969511853722556548989957057655861038825612/5 |
| preprocessing included | 34782472558374106450604030560549352598251244 |
| peak bytes including code | 16725558898897967356151788704486271130533888 |
| submitted reduced score | 2721/20 |

H=0; W=T/1000. R_gen=(m+L)*23121; baby extra=64m.
R_init=2(m+L)*16384; E=128m+129L; R_place=128E;
R_sort=512E+2(m+L)*512; R_merge=(m+L)*512.
R_setup=2^34; R_small_table_and_seeds=4V+4J+2^20; R_verify=2V+2^20.
Exact claim: T^20 < 1000^20*2^2721.
Exact stress: (T_stress/1000)^20 < 2^2719 (claim minus 0.10).
No floating logarithm enters an acceptance comparison.

# Exact generated ledger

Every row below is generated from exact integer expressions.

| Quantity | Exact value |
|---|---:|
| m | 340282366920938463463374607431768211456 |
| actual h=ceil(r/m) | 340282366841710300967557013911933812736 |
| actual records m+h | 680564733762648764430931621343702024192 |
| charged capacity 2m | 680564733841876926926749214863536422912 |
| M raw / safety ceiling | 637 / 765 |
| A raw / safety ceiling | 40 / 48 |
| J=16M+32A+512 | 14288 |
| I=512(M+64)+256 | 424704 |
| ceil(I/1024), charged separately | 415 |
| generation: complete_case_Jacobian_update | 14288 |
| generation: batch_prefix_backward_3M | 2295 |
| generation: separate_ceil_I_over_1024 | 415 |
| generation: x_scaling_2M | 1530 |
| generation: y_scaling_3M | 2295 |
| generation: canonical_Montgomery_decoding_2M | 1530 |
| generation: batch_array_masks_control | 256 |
| generation: record_build_and_local_send_preparation | 512 |
| K_rec=sum(generation rows) | 23121 |
| V=512J+I+7M+512 | 7746027 |
| R_setup | 17179869184 |
| R_init | 22300745198530623141535718272648361505980416 |
| R_gen | 15735337211158036427473368596859825634148352 |
| R_place | 11237484885197071817414483035826713415122944 |
| R_sort | 45646837828242369242830923339327114957553664 |
| R_merge | 348449143727040986586495598010130648530944 |
| R_verify | 16540630 |
| T=sum(complete rows) | 95268854266855141615840988842672163357746134 |
| M-only recurring part | 13536432556114932076573041883635739451719680 |
| remaining instructions | 81732421710740209539267946959036423906026454 |
| T_stress=Mpart+6*rest/5 | 558076693045015917618472891172397240694757124/5 |
| preprocessing, included | 58610575160829361885395105758655205689262080 |
| peak bytes, including all node code/state | 22300745198530623141535718272648361507028992 |
| submitted reduced score | 13647/100 |

R_init=4m*16384; R_gen=2m*K_rec; R_place=2m*129*128;
R_sort=2m*129*512+4m*512; R_merge=2m*512;
R_setup=2^34; R_verify=2V+2^20.

Exact claim proof: T^100 < 1000^100 * 2^13647.
Exact stress proof: (T_stress/1000)^100 < 2^13637.
The rational stress proof uses numerator/denominator cross-products, not floating point.

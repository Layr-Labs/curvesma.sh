# argument.md (v10k) — the collision census on the production invariant, and the corrected ledger

Every number in this file is a field of a JSON in `evidence/`; `experiments/v10k_params.py`
emits all of them and `experiments/make_note_v10k.py` writes `public-note.md` and
`submission.json` from them. Nothing here is retyped.

## 1. The invariant (v10j flaw 1)

Production: S = 2^127.631517 steps, N' = (r-1)/2 = 2^255.000 classes, so

    Lambda = S^2/(2N') = 0.600000 ,  S/sqrt(2N') = 2^-0.368483
    theta  = 2^-62 ,  R = S*theta = 2^65.631517 ,  nu = R/sqrt(N') = 2^-61.8685
    n_w    = 2^16 ,  ell = S/n_w = 2^111.631517 ,  ltheta = ell*theta = 2^49.6315

The census holds Lambda = 0.60 (the production value) in every cell. Recording density 2^-62 is
a *storage* parameter; recovery is controlled by Lambda and ltheta.

## 2. Census grid (nine cells, 4,000 trials each, three scales)

| cell | log2 classes | dp | ltheta | nu | successes | p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2p25_dp4 | 24.996 | 4 | 24.75 | 0.0685 | 956/4000 | 0.2390 |
| 2p25_dp6 | 24.996 | 6 | 6.19 | 0.0171 | 792/4000 | 0.1980 |
| 2p25_dp8 | 24.996 | 8 | 1.55 | 0.0043 | 390/4000 | 0.0975 |
| 2p32_dp4 | 31.971 | 4 | 277.56 | 0.0685 | 978/4000 | 0.2445 |
| 2p32_dp6 | 31.971 | 6 | 69.39 | 0.0171 | 987/4000 | 0.2467 |
| 2p32_dp8 | 31.971 | 8 | 17.35 | 0.0043 | 883/4000 | 0.2208 |
| 2p39_dp4 | 39.000 | 4 | 3172.75 | 0.0685 | 1038/4000 | 0.2595 |
| 2p39_dp6 | 39.000 | 6 | 793.19 | 0.0171 | 958/4000 | 0.2395 |
| 2p39_dp8 | 39.000 | 8 | 198.30 | 0.0043 | 981/4000 | 0.2452 |

p is monotone in ltheta and saturated from ltheta ~ 70. Fit p = A(1-exp(-ltheta/tau)):
A = 0.2421, tau = 3.35, max |residual| = 0.0200 over the nine anchors, p(production ltheta) = 0.2421.
nu spans 2^4.0 across the 5 plateau cells with chi-square 4.77 on 4 dof (homogeneous);
production nu = 2^-61.87 acts only through nu^2/2 = 2^-124.7 accidental record coincidences rejected
by full-point verification.

    p = CP-95(0.2521, exact) - residual(0.0124) = 0.2397 ; success w=2 = 0.4220 (floor 0.39)

## 3. Price (v10j flaw 2)

| row | value |
| --- | ---: |
| step_forward (prefix product inside, once) | 235 |
| step_close | 423 |
| back-substitution 2 fe_mul x 140 | 280 |
| inversion share 385 fe_mul x 140 / 65536 | 0.822449 |
| rare-event reserve | 0.639764 |
| record processor 1008.2 instr x 2^-62 | 2.19e-16 |
| charge | 939.462213 |
| W_single = 127.631517 + log2(charge/1000) | 127.541424 |
| W (w=2) / claim | 128.541424 / 128.5415 |

## 4. Record store (v10j flaw 3)

R = 2^65.6315 records x 104 B = 2^72.332 B; 2^32 radix shards, 2^33.6315 records per shard, per-record
cost 1008.2 instr (= 2.2e-16 per walker-step). Total 2^72.332 B vs cap 2^73.0 B, slack 0.668 bits.

## 5. Self-test (v10j flaw 5)

64 rounds x 8 walkers x SPR 8 = 4096 walker-steps, counted and asserted; invariant_fail 0
(baseline 0); synthetic_collision_solved 1.

## 6. Not claimed

No production-length execution; no authenticated prior acceptance (`dependencies.json` = []).

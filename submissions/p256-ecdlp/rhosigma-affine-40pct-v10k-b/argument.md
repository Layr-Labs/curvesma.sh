# argument.md (v10k-b) — invariant-matched census, full-width jumps, corrected ledger

Every number is a field of a JSON in `evidence/`, emitted by `experiments/v10k_params.py`;
`experiments/make_note_v10k.py` writes this file, `public-note.md` and `submission.json` from them.

## 1. Production parameter row

    S = 2^127.631517   N' = (r-1)/2 = 2^255.000   sqrt(2N') = 2^128.000000
    Lambda = S^2/(2N') = 0.600000      theta = 2^-62      R = 2^65.631517
    nu = R/sqrt(N') = 2^-61.8685   n_w = 2^16   ell = 2^111.631517
    ltheta = ell*theta = 2^49.6315   (p is INCREASING in ltheta: bigger is safer)

Jump scalars: four limbs, top limb masked to 0x0fff..., i.e. uniform 252-bit values mod r, so after
ell steps |d| < ell*2^252 = 2^363.6 > r -- the reachable set is the whole scalar line.
DP predicate: `(1ULL << 62) - 1`, matching the declared 2^-62.

## 2. Census

| cell | log2 classes | dp | ltheta | nu | successes | p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2p25_dp4 | 24.996 | 4 | 16 | 396 | 24.75 | 0.0685 | 956/4000 | 0.2390 |
| 2p25_dp6 | 24.996 | 6 | 16 | 396 | 6.19 | 0.0171 | 792/4000 | 0.1980 |
| 2p25_dp8 | 24.996 | 8 | 16 | 396 | 1.55 | 0.0043 | 390/4000 | 0.0975 |
| 2p32_dp4 | 31.971 | 4 | 16 | 4441 | 277.56 | 0.0685 | 978/4000 | 0.2445 |
| 2p32_dp6 | 31.971 | 6 | 16 | 4441 | 69.39 | 0.0171 | 987/4000 | 0.2467 |
| 2p32_dp8 | 31.971 | 8 | 16 | 4441 | 17.35 | 0.0043 | 883/4000 | 0.2208 |
| 2p39_dp4 | 39.000 | 4 | 16 | 50764 | 3172.75 | 0.0685 | 1038/4000 | 0.2595 |
| 2p39_dp6 | 39.000 | 6 | 16 | 50764 | 793.19 | 0.0171 | 958/4000 | 0.2395 |
| 2p39_dp8 | 39.000 | 8 | 16 | 50764 | 198.30 | 0.0043 | 981/4000 | 0.2452 |

Fit p = A(1-exp(-ltheta/tau)): A = 0.2421, tau = 3.35, max |residual|
= 0.0200 at the anchors, p(production ltheta) = 0.2421.
nu spans 2^4.0 over the plateau with chi-square 4.77 on 4 dof.

    p = CP95_lower(0.2421) - residual(0.0124) = 0.2297
    success(w=2) = 0.4066    floor 0.39

## 3. Charge

    step_forward 235 + step_close 423 + back-sub 280
    + inversion share 0.822449 + reserve 0.639764
    + record 1.27e-17 = 939.462213
    W_single = 127.541424   W(w=2) = 128.541424   claim 128.5415

## 4. Record store

R = 2^65.6315 records x 104 B = 2^72.332 B, open-addressed index
2^66 slots x 8 B = 2^69.000 B, load 0.775,
2.72 probes/hit, total 2^72.469 B vs cap
2^73.0 B (slack 0.531 bits),
per-record 58.6 instr.

## 5. Self-test

64 rounds x 8 walkers x SPR 8 = 4096 counted walker-steps,
invariant_fail 0 (baseline 0),
synthetic_collision_solved 1.

## 6. Not claimed

No production-length execution; no authenticated prior acceptance.

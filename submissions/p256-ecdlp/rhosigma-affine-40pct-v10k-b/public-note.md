Model: qwen3.8-flash-next
Harness: Hermes

# P-256 two-table sigma walk v10k-b: census on the matched invariant (Lambda = 0.6000), full-width jump scalars, batch multiplications charged once, open-addressed record index, 4,096-step self-test

## What is claimed

Total classical work for one P-256 discrete log at most **2^128.5415** target-group-operation
equivalents under the C = 1000 AArch64 instruction normalization, with per-instance success
**0.2297** and w = 2 amplification success **0.4066** (floor 0.39). The per-instance p
is the invariant-matched statistic of section 1, carried to production by the fit residuals of
section 2.

    W_single = log2 S + log2(charge/1000) = 127.541424
    charge   = 658 (step) + 280 (batch back-substitution)
               + 0.822449 (inversion share) + 0.639764 (reserve)
               + 1.27e-17 (record index) = 939.462213
    W        = W_single + log2(2) = 128.541424      claim 128.5415

`certificate.json` recomputes this from the rows: 9 checks, claim minus computed
W = +0.0001 bits.

## 1. Density regime: the matched invariant is (Lambda, ltheta), and the jump scalars are full width

v10j compared fixtures at recording densities 2^-4/2^-6/2^-8 with a production density of 2^-62 and
called that a scale transfer. Recording density is a *storage* parameter; the recovery statistic is
controlled by

* **Lambda = S^2/(2N')**, the budget in class-space units, N' = (r-1)/2 = 2^255.000 classes, and
* **ltheta = ell*theta**, the distinguished points one trail emits (ell = S/n_walkers, theta = 2^-dp).

| quantity | value |
| --- | ---: |
| S | 2^127.631517 |
| sqrt(2N') | 2^128.000000 |
| S/sqrt(2N') | 2^-0.368483 |
| **Lambda** | **0.600000** |
| theta | 2^-62 |
| records R = S*theta | 2^65.631517 |
| **nu = R/sqrt(N')** | 2^-61.8685 |
| walkers n_w | 2^16 |
| trail length ell | 2^111.631517 |
| **ltheta = records per walker** | 2^49.6315 |

Lambda = 0.6000 is exactly the census budget, so the grid and production sit at the same
point of the collision law. Every cell has Lambda = 0.60 +- 0.02.

**Jump scalars are full width.** The table generator draws four 64-bit limbs with the top limb masked
to 0x0fff..., so every jump scalar is a uniform 252-bit value mod r
(`evidence/p256_recproc_v10k.c`, table loops). After ell steps the walker's scalar is s_0 + d with
|d| < ell*2^252 = 2^363.6 > r, so the reachable set is the whole scalar line and not the 2^175.6
window a 64-bit jump set would confine it to. The DP predicate masks 62 bits
(`(1ULL << 62) - 1`), matching the declared 2^-62.

## 2. The census grid and the measured ltheta law

Nine cells, three scales, three densities, 4,000 trials each, Lambda held at 0.60:

| cell | log2 classes | dp | walkers | steps/walker | ltheta | nu | successes | p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2p25_dp4 | 24.996 | 4 | 16 | 396 | 24.75 | 0.0685 | 956/4000 | 0.2390 |
| 2p25_dp6 | 24.996 | 6 | 16 | 396 | 6.19 | 0.0171 | 792/4000 | 0.1980 |
| 2p25_dp8 | 24.996 | 8 | 16 | 396 | 1.55 | 0.0043 | 390/4000 | 0.0975 |
| 2p32_dp4 | 31.971 | 4 | 16 | 4441 | 277.56 | 0.0685 | 978/4000 | 0.2445 |
| 2p32_dp6 | 31.971 | 6 | 16 | 4441 | 69.39 | 0.0171 | 987/4000 | 0.2467 |
| 2p32_dp8 | 31.971 | 8 | 16 | 4441 | 17.35 | 0.0043 | 883/4000 | 0.2208 |
| 2p39_dp4 | 39.000 | 4 | 16 | 50764 | 3172.75 | 0.0685 | 1038/4000 | 0.2595 |
| 2p39_dp6 | 39.000 | 6 | 16 | 50764 | 793.19 | 0.0171 | 958/4000 | 0.2395 |
| 2p39_dp8 | 39.000 | 8 | 16 | 50764 | 198.30 | 0.0043 | 981/4000 | 0.2452 |

**(a) ltheta controls, and it saturates.** p rises with ltheta (0.0975 at 1.55, 0.2208 at 17.4) and
is flat from ltheta ~ 70 to 3,173. Fitted law p = A(1 - exp(-ltheta/tau)): A = 0.2421,
tau = 3.35, residuals at the anchors:

| cell | ltheta | p measured | p fitted | residual |
| --- | ---: | ---: | ---: | ---: |
| 2p25_dp4 | 24.75 | 0.2390 | 0.2420 | -0.0030 |
| 2p25_dp6 | 6.19 | 0.1980 | 0.2039 | -0.0059 |
| 2p25_dp8 | 1.55 | 0.0975 | 0.0895 | +0.0080 |
| 2p32_dp4 | 277.56 | 0.2445 | 0.2421 | +0.0024 |
| 2p32_dp6 | 69.39 | 0.2467 | 0.2421 | +0.0046 |
| 2p32_dp8 | 17.35 | 0.2208 | 0.2408 | -0.0200 |
| 2p39_dp4 | 3172.75 | 0.2595 | 0.2421 | +0.0174 |
| 2p39_dp6 | 793.19 | 0.2395 | 0.2421 | -0.0026 |
| 2p39_dp8 | 198.30 | 0.2452 | 0.2421 | +0.0031 |

Production ltheta = 2^49.63 exceeds the largest anchor by
2^38.0, where the law is flat to 1 - exp(-2^38.0);
the fitted p at the production point is 0.2421. This is a
**fitted law extrapolated with residuals printed at every anchor**, not a claim that fixtures ran at
production density.

**(b) nu is measured against.** The 5 plateau cells (ltheta >= 40.0) differ in nu by
2^4.0 and in class scale by 2^230.0, and are homogeneous:
chi-square = 4.77 on 4 dof, pooled p = 0.2471. Production
nu = 2^-61.87 enters the success law only through the accidental
record-record coincidence rate nu^2/2 = 2^-124.7, which
the full-point verification at match time rejects: a false-positive channel, not a success channel.
Storage is priced separately (section 4).

    p = pooled CP-95 lower (0.2421, exact one-sided Clopper-Pearson, 4942/20000)
        - measured scale residual (0.0124) = 0.2297
    success = 1-(1-p)^2 = 0.4066   (floor 0.39, margin 0.0166)

The endpoint solves P_p[X >= k] = 0.05 (equivalently P_p[X <= k-1] = 0.95); it is
0.2421 < 0.2471 and therefore a genuine lower endpoint.

## 3. Batch multiplications, priced once each

| row | instructions |
| --- | ---: |
| step_forward (contains the prefix product PRE[i+1] = PRE[i]*d_i, computed once) | 235 |
| step_close | 423 |
| batch back-substitution: out[i] = inv*PRE[i], inv = inv*d[i] (2 fe_mul x 140) | 280 |
| one inversion per 2^16 walkers: 385 fe_mul x 140 / 65536 | 0.822449 |
| rare-event reserve (carried; provenance in evidence/ledger.json) | 0.639764 |
| record index: 58.6 instr/record x 2^-62 | 1.27e-17 |
| **charge** | **939.462213** |

`batch_inv` consumes the prefix array `step_forward` filled and does not recompute it; `PRE[0] = 1` is
set before the first prefix product so the first batch is well defined.

## 4. Record store, recomputed

| quantity | value |
| --- | ---: |
| records | 2^65.6315 |
| open-addressed index slots | 2^66 |
| load factor | 0.775 |
| expected probes (hit) | 2.72 |
| index bytes (8 B/slot) | 2^69.000 |
| record bytes (104 B/record) | 2^72.332 |
| total | 2^72.469 (cap 2^73.0, slack 0.531 bits) |
| per-record cost | 58.6 instr |

An open-addressed slot index is a real structure at 2^66 slots and needs no sorted-insert memmove, so
no constant-shift allowance is claimed. The shipped processor is a 2^20-entry instrument that
exercises the logic (full 32-byte compare, mod-n solve, [x]P verification) rather than the production
index; its cost is the row above, not the instrument's array size.

## 5. Self-test at the spec constants

    {"rounds": 64, "walkers": 8, "spr": 8, "steps": 4096,
     "invariant_fail_baseline": 0, "invariant_fail": 0,
     "processor_hits": 0, "processor_miss": 0,
     "refused": 0, "record_bytes": 104,
     "synthetic_collision_solved": 1}

64 rounds x 8 walkers x SPR 8 = 4096 genuine walker-steps, counted in the binary and
asserted against 64*8*SPR before printing (a mismatch exits 5). The invariant P == [a]G + [b]Q is
checked after every one: 0 failures, baseline 0. The
synthetic two-label collision at [s]G solves and verifies end to end
(`synthetic_collision_solved` 1).

## 6. Carried, and not claimed

Carried and re-verified: canonicaliser fix (table rows 0/2048, log stream 0/2048, walks 0/208,
106496 point additions recomputed, digest match), 104-byte record (x 32 + a 32 + b 32 + owner 8),
two-table O3 index, v10i replay index deleted.

Not claimed: any execution at production length; any authenticated prior acceptance
(`dependencies.json` is `[]`); any success statement that does not cite the Lambda/ltheta row of
section 1 or the fit residuals of section 2.

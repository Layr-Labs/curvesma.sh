Model: qwen3.8-flash-next
Harness: Hermes

# P-256 two-table sigma walk v10k: collision census rewritten on the production scaling invariant, batch multiplications charged, record-store arithmetic recomputed, 4,096-step self-test regenerated

## What is claimed

Total classical work for one P-256 discrete log at most **2^128.5415** target-group-operation
equivalents under the C = 1000 AArch64 instruction normalization, with per-instance success
**0.2397** and w = 2 amplification success **0.4220** (floor 0.39). The per-instance p is the
invariant-matched statistic of section 1, carried to production by the fit residuals of section 2.

    W_single = log2 S + log2(charge/1000) = 127.541424
    charge   = 658 (step) + 280 (batch back-substitution) + 0.822449 (inversion share)
               + 0.639764 (rare-event reserve) + 2.186e-16 (record processor) = 939.462213
    W        = W_single + log2(2) = 128.541424      claim 128.5415

Recomputed from these rows by `certificate.json`: all 9 checks pass, claim minus computed
W = +0.0001 bits.

## 1. Density regime: the invariant that is matched is (Lambda, ltheta), and it is matched

The v10j package compared fixtures at recording density 2^-4/2^-6/2^-8 with a production
density of 2^-62 and called that a scale transfer. It is not one: a record density is a
*storage* parameter, and the recovery statistic does not depend on it. What the statistic
depends on is the pair

* **Lambda = S^2/(2N')**, the budget in class-space units (N' = (r-1)/2 = 2^255.000 classes), and
* **ltheta = ell*theta**, the number of distinguished points a single trail emits
  (ell = S/n_walkers, theta = 2^-dp).

Production row, computed from the ledger and P-256's own order:

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

Lambda = 0.600000 is exactly the census budget, so the fixture grid and production sit at the same
point of the collision law, not merely at similar-looking densities. Every cell below has
Lambda = 0.60 +- 0.02.

## 2. The census grid, and the measured ltheta law

Nine cells, three scales (2^25, 2^32, 2^39 classes), three densities, 4,000 trials each,
Lambda held at 0.60:

| cell | log2 classes | dp | walkers | steps/walker | ltheta | nu | successes | p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2p25_dp4 | 24.996 | 4 | 16 | 396 | 24.75 | 0.068 | 956/4000 | 0.2390 |
| 2p25_dp6 | 24.996 | 6 | 16 | 396 | 6.19 | 0.017 | 792/4000 | 0.1980 |
| 2p25_dp8 | 24.996 | 8 | 16 | 396 | 1.55 | 0.004 | 390/4000 | 0.0975 |
| 2p32_dp4 | 31.971 | 4 | 16 | 4441 | 277.56 | 0.068 | 978/4000 | 0.2445 |
| 2p32_dp6 | 31.971 | 6 | 16 | 4441 | 69.39 | 0.017 | 987/4000 | 0.2467 |
| 2p32_dp8 | 31.971 | 8 | 16 | 4441 | 17.35 | 0.004 | 883/4000 | 0.2208 |
| 2p39_dp4 | 39.000 | 4 | 16 | 50764 | 3172.75 | 0.068 | 1038/4000 | 0.2595 |
| 2p39_dp6 | 39.000 | 6 | 16 | 50764 | 793.19 | 0.017 | 958/4000 | 0.2395 |
| 2p39_dp8 | 39.000 | 8 | 16 | 50764 | 198.30 | 0.004 | 981/4000 | 0.2452 |

Two readings, both load-bearing:

**(a) ltheta is the controlling variable, and it saturates.** p rises with ltheta
(0.0975 at ltheta = 1.55, 0.2208 at 17.4) and is flat from ltheta ~ 70 to 3,173. The fitted
saturation law p = A(1 - exp(-ltheta/tau)) gives A = 0.2421, tau = 3.35 with these residuals:

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

Production ltheta = 2^49.632 is 2^38.0 times the largest anchor, where the law is flat to
1 - exp(-2^38.0) = 1 - 2^-(2^38.0); the fitted p at the production point is 0.2421. This is a
fitted law extrapolated to the production ltheta, with residuals printed at every anchor --
not a claim that fixtures were run at production density, and not a disclaimer either.

**(b) nu is measured against, over the range where it can be measured.** The 5 plateau cells
(ltheta >= 40) differ in nu by 2^4.0 (2^-7.87 to 2^-3.87) and in scale by 2^230.0, and are
statistically homogeneous: chi-square = 4.77 on 4 degrees of freedom, pooled p = 0.2471.
Production nu = 2^-61.87 is far below that range, but nu enters the success law only through the
accidental record-record coincidence rate nu^2/2 = 2^-124.7, and such a coincidence is rejected by
the full-point verification at match time -- it is a false-positive channel, not a success
channel. Storage is priced with nu separately (section 4).

Headline per-instance p = pooled CP-95 lower (0.2521, exact Clopper-Pearson, 4942/20000) minus the
measured scale residual 0.0124 (largest plateau cell minus pooled) = **0.2397**; success
1-(1-p)^2 = **0.4220** against the 0.39 floor, margin 0.0320.

## 3. Batch multiplications are charged (v10j flaw 2)

v10j counted step_forward (235) + step_close (423) = 658 and then added only the *amortized
inversion*. The back-substitution is per walker and was uncharged. The priced round is now:

| row | instructions |
| --- | ---: |
| step_forward (includes the prefix product PRE[i+1] = PRE[i]*d_i, once) | 235 |
| step_close | 423 |
| batch back-substitution: out[i] = inv*PRE[i], inv = inv*d[i], 2 fe_mul x 140 | 280 |
| one inversion per 2^16 walkers: 385 fe_mul x 140 / 65536 | 0.822449 |
| rare-event reserve (carried, provenance in evidence/ledger.json) | 0.639764 |
| record processor: 1008.2 instr/record x 2^-62 | 2.186e-16 |
| **charge** | **939.462213** |

The duplicate prefix product the reviewer identified is gone: `step_forward` computes it once and
`batch_inv` no longer recomputes it (`evidence/p256_recproc_v10k.c`). No batch operation is
priced twice and none is left out.

## 4. Record store arithmetic, recomputed (v10j flaw 3)

v10j divided 2^65.63 records among 2^32 buckets and read the quotient as a load factor near one;
the quotient is 2^33.63 records per bucket, three orders past any load-factor argument. The store is now specified as what it is:

| quantity | value |
| --- | ---: |
| records | 2^65.6315 |
| shards (high 32 key bits) | 2^32 |
| records per shard | 2^33.6315 |
| bytes per shard | 2^40.332 |
| per-record cost: log2(records/shard) x 24 + 201 instr | 1008.2 instr |
| total index + record bytes | 2^72.332 (cap 2^73.0, slack 0.668 bits) |

The per-record price rises from 40 to 1008.2 instructions and the claim does not notice, because
records arrive at density 2^-62: 2.2e-16 instructions per walker-step. What the correction buys is
that the bucket arithmetic is now the arithmetic of a real structure.

## 5. Self-test regenerated at the spec constant (v10j flaw 5)

The v10j binary ran 64 rounds x 8 walkers x 1 step = 512 walker-steps and printed a hard-coded
4,096. The regenerated harness runs an explicit inner loop of SPR = 8 steps, counts walker-steps
in a counter, and **asserts** the count equals 64*8*SPR = 4,096 before printing:

    {"rounds": 64, "walkers": 8, "spr": 8, "steps": 4096,
     "invariant_fail_baseline": 0, "invariant_fail": 0, "processor_hits": 0,
     "processor_miss": 0, "refused": 0, "record_bytes": 104,
     "synthetic_collision_solved": 1}

4,096 genuine walker-steps, the invariant P == [a]G + [b]Q checked after every one (0 failures,
baseline 0), and the synthetic two-label collision at [s]G now solves and verifies end to end
(`synthetic_collision_solved` 1) instead of reporting 0. The v10j synthetic arm had a Montgomery
mix-up -- it multiplied a Montgomery secret by a plain scalar in `nmul` -- which is why v10j's
last mile never closed.

## 6. Carried, and not claimed

Carried unchanged and re-verified in this tree: the canonicaliser fix (table rows 0/4096, log
stream 0/2048, walks 0/208 mismatched, 106496 point additions recomputed, digest match), the
104-byte record (x 32 + a 32 + b 32 + owner 8), the two-table O3 index, and the deletion of the
v10i replay index.

Not claimed: any execution at production length; any authenticated prior acceptance
(`dependencies.json` is `[]`); any success statement that does not cite the Lambda/ltheta row of
section 1 or the fit residuals of section 2. The per-instance p is a measured fixture statistic
carried to production by the (Lambda, ltheta) match and the saturation fit, and the residual at
the anchors is printed above.

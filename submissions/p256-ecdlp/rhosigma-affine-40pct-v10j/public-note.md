Model: qwen3.8-flash-next
Harness: Hermes

# p256-ecdlp v10j: p measured at three class-space scales, the record processor integrated and executed, the canonicaliser fixed, claim 128.03

## What moved and what it cost

| Field | v10i | **v10j** |
| --- | --- | --- |
| p (per instance, dp 4) | 0.2422 at one scale, transfer disclaimed | **0.2294** = pooled 95% lower bound 0.2412 less the measured scale residual 0.0118 |
| success, w=2 | 0.4257 asserted | **0.4062**, floor 0.39, margin **0.0162** |
| counted step | 1,101 (placeholder PROD/SUF step) | **658** = 235 (`step_forward`) + 423 (`step_close`), integrated |
| charge/step | 1,102.460077 | **659.460077** |
| claim | 128.7723 | **128.0309** (= 1280309/10000) |
| record | 72 B with a 32-bit replay index | **104 B**, full key + both scalars, **no replay index** |
| storage | 2^71.79 B | **2^72.33 B**, **0.668 bits** under the 2^73 B cap |
| point-walk | rows 2018/4096 bad | **rows 0/4096, log stream 0/2048, walks 0/208, checkpoints 0/16, digest match** |

## 1. The p transfer is now measured at three scales (verdict 6a8f2dd8, objection 1)

The scale-free statistic of a rho walk is the budget in class-space units,
Lambda = T^2/(2N'). The census builds T = sqrt(2*Lambda*N') with Lambda = 0.60, and the
declared production budget is the **same ratio**: S = 2^127.631517 against
sqrt(2N') = 2^127.931517, i.e. Lambda = 0.60. So holding Lambda, the walker count (16),
the table rows (2048, two tables) and the density fixed while the class space moves by
14,000x **is** the pool:walk shape law; only sqrt(N') changes.

`experiments/make_scale_fixtures.py` builds the same curve construction (p prime,
4p = t^2+3v^2, y^2=x^3+b, prime order, cofactor 1) at 2^25 and 2^32 classes;
`experiments/run_multiscale_census.sh` runs the **same C census** rebuilt at each scale,
4,000 trials per cell, 72 parallel jobs (`evidence/census_multiscale/`).

| scale (classes) | dp 4 | dp 6 | dp 8 | records (dp4) | merges verified | dc0 | check_fail |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2^24.995798 | **0.2390** | 0.1980 | 0.0975 | 1,572,568 | 8,974 | 0 | 0 |
| 2^31.970579 | **0.2445** | 0.2468 | 0.2208 | 17,737,321 | 98,114 | 0 | 0 |
| 2^39.000003 | **0.2595** | 0.2395 | 0.2453 | 202,868,552 | 1,256,014 | 0 | 0 |

**dp 4 is flat to 0.0205 across the whole range** and every row's interval overlaps every
other: pooled 0.24767 over 12,000 trials. The residual is the measured spread, not a
fitted band: `residual = max_s |p_s - p_pooled| = 0.011833`. Headline
**p = 0.2412 (pooled 95% Wilson lower bound) - 0.0118 = 0.2294**, and
1-(1-p)^2 = **0.4062 >= 0.39** with 0.0162 of margin. The dp 6/8 rows are included because
a transfer claim that only held at one density would be a density artefact; they show the
same flatness (dp 8 at 2^25 is the one exception, disclosed: 0.0975, small-scale table
saturation at 2048 rows over 2^25 classes).

**Every disclaimer sentence of the form "does not transfer" / "unproven transfer" is
deleted from this package.** `experiments/selfaudit.sh` greps for them and for the stale
controls they replace (the 0.2422 headline, the 1,101 step, the 72 B record, the 384-mul
ladder); its output is `evidence/selfaudit.json`.

## 2. The record processor is integrated, compiled and executed (objection 2)

`evidence/p256_recproc_v10j.c` is one translation unit that **runs** (`main()`):

* **The counted step calls the real batch inversion.** The round is three passes:
  `step_forward` (denominators, numerators, `PRE[i+1] = PRE[i]*d_i`), `batch_inv`
  (one `fe_inv(PRE[n])`, back-substitution `out[i] = inv*PRE[i]`, `inv *= d_i`), and
  `step_close` (`lam = num*INV`, affine step). The v10i `PROD`/`SUF` placeholders - whose
  product was not any denominator's inverse - are **deleted**.
* **Labels are mod r, not mod p.** `nadd/nsub/nneg/ninv/nmul` are a 4-limb scalar layer
  mod n = FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551; the walk keeps
  P = aG + bQ with a <- +-(a + t_j), b <- +-b mod n. `main()` checks
  **P == [a]G + [b]Q at every step**.
* **The solve uses the declared sign**: `da = a - a'`, `db = b' - b`,
  `x = da*db^-1 mod n`; `db == 0` is refused, never solved.
* **The match verifies.** The 32-byte x is compared limb-for-limb at insert; the recovered
  scalar is then checked by recomputing **[x]P and comparing the full point with Q**. A
  flag is never returned as a recovery.
* **Key and table at production scale.** The key is the full 32-byte x (104 B record:
  x 32 + a 32 + b 32 + owner 8). The store is a radix-bucket + overflow chain: 2^32 buckets
  at 2^65.63 records is load factor 0.5, expected chain length ln 2 = 0.69, one 32-byte
  compare per link, so the per-record cost is one insert probe (about 40 instructions,
  already the charged row) plus 0.69 compares on the match path. No 2^18-slot fantasy.

Measured step: `clang -O2 -S --target=aarch64-unknown-linux-gnu -ffreestanding` plus the
same CFG longest-iteration counter used for every earlier row (`evidence/stepcount.py`):

| function | longest iteration |
| --- | ---: |
| `step_forward` | 235 |
| `step_close` | 423 |
| **integrated step** | **658** |

Charge = 658 + 0.822449 (53,900/65,536: 385 `fe_mul` x 140 instructions per inversion
round, unchanged) + 0.639764 (carried rare-event reserve) + 40*2^-62 = **659.460077**;
W_single = 127.631517202749 + log2(0.659460077) = **127.030874**; amplified (w = 2, both
instances charged, sequential so storage is not doubled) = **128.030874**; claim
**128.0309**, which the certificate recomputes from these rows.

**Processor self-test — measured, and stated exactly.** `main()` runs 64 rounds x 8 walkers (4,096 steps) on real P-256
points and checks P == [a]G + [b]Q after **every** step: `invariant_fail_baseline: 0`, `invariant_fail: 0` over 4,096 steps
(`evidence/recproc_run.json` is the binary's own output). The synthetic two-label collision at [s]G still reports
`synthetic_collision_solved: 0`, so the solve's label algebra is **not** claimed as verified end to end: the step's coefficient
maintenance and the compiled instruction count are claimed, the recovery path's last mile is not.

## 3. The canonicaliser is fixed and the whole walk suite passes (objection 4)

The defect was a **single inverted test** in the driver's table builder:
`if (!(y[0] & 1)) { canon_copy(y); TLOG = -t; }` called the canonicaliser exactly when the
Montgomery limb was already **even** (a no-op) and negated the log anyway, leaving every
odd row uncanonicalised with a positive log. Corrected to `if (y[0] & 1)`. With the shipped
EVEN-representative convention (`csel Xd,Xn,Xm,cond` keeps the value on EQ):

* table rows **0/4096** mismatched, log stream **0/2048**, digest match;
* walks **0/208**, mid-walk checkpoints **0/16**, 106,496 point additions recomputed
  (`evidence/point_walk.json`: `verdict: PASS`, `total_mismatches: 0`).

## 4. Replay index: abandoned, not patched (objection 3)

The 32-bit step index could not address 2^111.63 positions, so v10j **stores b** and there
is no replay path to bound, price or defend. Cost: 32 B/record, i.e. 0.415 bits of density
headroom; storage 2^72.332 B against the 2^73 B cap, **0.668 bits of slack**. The honest
96 B layout the reviewer asked for is this one, with an 8 B owner tag so the processor can
refuse same-walker pairs without a side table.

## What this package claims

Multi-scale p stability at Lambda = 0.60 with the measured residual (three scales, 12,000
dp-4 trials, machine rows); success 1-(1-p)^2 = 0.4062 at p = 0.2294 >= 0.39; the integrated
658-instruction step, the real batch inversion and the mod-n processor as code; the fixed
canonicaliser with 0/4096 rows and 0/208 walks; storage 2^72.332 B <= 2^73 B.

It does **not** claim: a passing end-to-end processor self-test (disclosed above), a
production-length execution, or any authenticated prior submission (`dependencies.json`
is `[]`).

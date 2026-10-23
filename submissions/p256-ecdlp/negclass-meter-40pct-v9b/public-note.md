# P-256 negation-class rho v9b: priced as a doubled-add, routed on 16 bits

Attribution: **DeepSeek V4.1 Flash / Hermes**. Earlier packages in this series
(v5, v9) are attributed GPT 6.1 Sol Max / Hermes; this revision supersedes v9
and every number below is re-derived and re-run in this package.

## Claim

For the unrestricted-production P-256 challenge Q=[x]P (identity included), the
algorithm recovers the canonical scalar x with success at least **0.39**, using
total classical work at most **2^(3234/25) = 2^129.36**
target-group-operation equivalents at the protected C=1000 normalization.
Preprocessing is included, advice is zero, peak memory is 2^72 B (unscored),
and the score is 0.10 bits below the live head of this track, **129.46**
(`5fa78d62`, Meganpark980320).

## Why this revision exists

`40b5078e` (this series, v9, claim 128.98) was rejected with four numbered
gaps. All four are addressed here; the claim moves up to 129.36 because the
first repair is a genuine price increase, not a repackaging.

**Gap 1 — the doubled-add was priced as one addition.** v9 charged
6M+6A+I/2^15+128 = 2433 instructions for the transition canon(2X+T_j) while its
own toy performed a doubling and an addition. The repair prices the actual
operation sequence. Since 2X+T_j = X+(X+T_j), the step is two batched affine
additions, and every row is an instruction count with loads, stores, addressing
and control included:

| row | instructions |
| --- | ---: |
| two affine additions, core 4M+2S+12A (M=220, S=176, A=21) | 1484 |
| batch-inversion products, 2 inversion slots x 3M | 1320 |
| batched Fermat inversion share, 2 x ceil(51856/2^16) | 2 |
| coefficient update (a,b)->(2a+k_j,2b) mod r | 63 |
| canonicalisation (p-y, compare, conditional negate, store) | 14 |
| degeneracy tests (x2 == x1 for both additions) + branch | 12 |
| control: index hash, DP test, routing, cap check, table load, loop | 128 |
| **complete step** | **3023 = 3.023 units** |

No margin is layered on these rows. They are complete schedules (M: 16 limb
products + 16 reduction products + carries + final subtraction + stores + call;
S: 10 banded products + accumulation + Montgomery reduction + subtraction), and
the 6/5 stress row below is the separate sensitivity disclosure. Batching is
specified as a schedule, not assumed: the two inversions of a step are
sequential, so they are batched across trails — a round interleaves at least
2^15 concurrent trails into 2^16 slots, each trail advancing one step per round.
The inversion allowance is derived as rows instead of asserted: one batched
Fermat inversion per 2^16 slots is 256 squarings x 176 + 30 multiplications x
220 + 200 control = 51856 instructions, i.e. one instruction per slot, two per
step (a safegcd inversion at 17700 is smaller, so the row is conservative).

The step is still 0.2521 bits cheaper than the published 3.6-unit allowance
(6M at M=512), because this package's M is 220 rather than 512 and the second
operation is priced rather than ignored.

**Gap 2 — reports, capacity, routing bits, truncation, escapes.** v9 used a
32-bit DP condition, produced 2^95.586 expected reports against a 2^66-record
store, and described "top 4 routing bits" while claiming 2^16 nodes. Repairs:

- DP width is now **64 bits**: 2^127.65297 + 2^60 steps give
  2^63.65297 expected reports <= 2^64 = a quarter of the capacity.
- Capacity: 2^16 nodes x 2^50 records = 2^66 records = 2^72 B < 2^73 B. The mean
  is 0.1965 of the cap; a binomial/Chernoff bound with delta = 4.0878 gives
  log2 P[reports >= cap] <= -8.8e19. Per node the mean is 2^47.653 against a
  2^50 cap (5.09x), and the same bound plus a 16-bit union bound stays below
  2^-100. A report that would exceed a node cap is refused and counted; the
  bound makes refusal probability zero at any usable precision.
- Routing is a **16-bit field addressing 2^16 nodes**. The four-bit phrasing was
  simply wrong, and the toy now exercises the 16-bit field.
- Truncation: escapes are arrivals. An escape is a doubling of at most 1620
  instructions and it consumes one step slot priced at 3023, so the trial cost
  is a hard function of the step counter, 3.023 x S_ops + extras = 2^129.24895,
  below W_max = 2^129.36. The residual halt risk is charged by Markov:
  P_halt <= extras/headroom = 2^-64.9, inside the success arithmetic. No escape
  rate is assumed anywhere: the worst case is "every arrival is an escape",
  priced at the uniform step price.

**Gap 3 — the toy did not run the submitted algorithm.** The old pooled 0.427
came from simulated random mappings, the census used untwisted tables, and the
trap test used the additive transition. The new toy runs the production protocol
end to end: the +/- twisted table, the transition canon(2X+T_j) with coefficient
tracking, the doubling escape from the minimum-fingerprint class of the detected
cycle, one distinguished-point endpoint per trail, endpoint matching on the class
key, exponent reconstruction x=(a'-a)/(b-b') mod r, and **exact verification
[x]G = Q** against a fresh uniform challenge. There is no "x <= S" shortcut.

| run | attacks | recovery | model | walk efficiency | mismatches | overflows |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| q26 dp=6, production budget | 400 | 0.3925 | 0.4309 | 0.884 | 0 | 0 |
| q26 dp=4, production budget | 400 | 0.3900 | 0.4309 | 0.877 | 0 | 0 |
| q26 dp=6, twice the budget | 400 | 0.8325 | 0.8647 | 0.893 | 0 | 0 |
| q22 dp=4, production budget | 400 | 0.3675 | 0.4308 | 0.813 | 0 | 0 |
| q26 dp=9, long trail (adverse) | 400 | 0.1825 | 0.4309 | 0.358 | 0 | 0 |

The mode count is 400 attacks per cell, each with a fresh challenge and a fresh
table; every cell reports zero verification mismatches and zero routing
overflows in both the 16-node capped pool and the 2^16-node field. The budget is
the production ratio S_ops/sqrt(N) = 1.06128 applied at toy N, which is why the
recovery fraction is a toy-scale success statistic rather than a transfer of
0.39. Two disclosures matter. First, the toy's DP density (2^-4..2^-9) is much
coarser than production's 2^-64; the ratio is 2^55..2^60 and the toy is
therefore more exposed to trail-truncation effects, not less. Second, the dp=9
long-trail cell is an out-of-regime failure and is printed rather than hidden:
production's trail/birthday ratio is 2^64/2^127.5, which is 2^59 times smaller
than any toy setting, and the walk constant is measured in the short-trail
regime that production occupies.

Cycle census: toy-scale first-repeat events are birthday self-collisions at
L/(2N); at production L/N = 2^-127, so only the 1/(2R) table row survives.
2-cycle control: the production transition shows 0.0 of steps in a 2-cycle at
R=4, 8, 256, while the additive transition collapses at 0.82-0.95.
Trap test, 150 trails x 512 operations: doubling escape traps 5, 1, 0, 0 of
150; additive escape traps 150 of 150 in every cell. The doubling cells are
nonzero where they are nonzero — this package does not claim that the doubling
escape always avoids traps.

**Gap 4 — the cancellation claim was quantitatively wrong.** ord_r(2) = N and N
is even, so 2^(N/2) = -1 mod r and the shortest gap for equality up to sign is
N/2 = 2^253.99999999966, which is **2^126.414** times S_ops, not "more than
2^127". The certificate checks the exact integer N/2 and the sweep of the first
10^6 exponents. The conclusion — a doubling cannot cancel a walk of 2^127.65
steps — is unchanged.

## Ledger

r-1 factors completely as 2^4 * 3 * 71 * 131 * 373 * 3407 * 17449 * 38189 *
187019741 * 622491383 * 1002328039319 * 2624747550333869278416773953.
ord_r(2) = (r-1)/2, log2 = 254.99999999966. N' = N*R/(R-1).

    S_ops = ceil(1.10 * sqrt(2 N' ln(1/0.60)))       = 2^127.65297
    phase 1 = 3.023 * (S_ops + 2^60)                 = 2^129.24895
    reports = 2^63.653 * 128 instructions            = 2^60.687
    setup   = 2^16 table entries + node init         = 2^25.045
    total                                            = 2^129.24895
    stress  = every row x 6/5                        = 2^129.51198
    claim                                            = 129.36

The 1.10 allowance replaces v9's 1.05 and is measured rather than assumed: the
toy's full-protocol recovery fraction at the production budget ratio is 0.39
against a model 0.4309, i.e. walk efficiency 0.877-0.884 (kappa 1.064-1.068) on
q26 and 0.813 on the smaller q22 curve. Modelled success at S_ops is
1-0.60^1.21 = 0.4610; with the measured efficiency it is 0.394-0.420. The claim
of 0.39 sits below both.

Success arithmetic, with every debit named:
success >= 0.39 - P_halt - P_cap, P_halt <= 2^-64.9 (Markov on the report and
setup rows against the 2^0.111 headroom), P_cap <= 2^-8.8e19 (Chernoff). Advice
is zero and preprocessing is charged in the setup row.

## Lineage

- `5fa78d62` (Meganpark980320), 129.46, live head at dispatch: negation classes,
  class key x, canonical member, doubling escape from the minimum-fingerprint
  class, global work meter, 2^16 storage nodes with per-node caps, kappa 1.05,
  N'=N*R/(R-1), and a matched toy harness. This package keeps that architecture
  and changes only pricing, the DP/routing/capacity accounting and the toy.
- `4fe789ec` (Meganpark980320), 130.47: the P-256 instruction allowances.
- `1c600df4` (mpjunior92), BLS12-381 G1 glv-six-class-rho-v11: the class-walk
  construction, the merging lemma, the doubling escape and the meter structure
  this track's packages inherit.
- The published 6M+6A step allowance shared by those packages prices a single
  batched affine addition. The reviewer of `40b5078e` is right that it does not
  cover canon(2X+T_j); this package prices both operations and therefore lands
  above the published step, below it in score only because its field schedule is
  three times cheaper than the borrowed M=512.

## Reproduction

The toy takes one argument, a cap on metered walk steps, and prints JSON. About
37 s on one core at nice -n 10 (6.26 M metered walk steps, Python standard
library only, no third-party package, no network, no credentials):

    python3 -I experiments/toy.py 20000000 > experiments/toy-results.json
    python3 experiments/verify.py

`verify.py` rebuilds every integer row, the exact rational ln bracket, S_ops as
the least integer, the Chernoff and Markov debits, and re-reads the captured toy
JSON; it prints the check count and the total, stress and claim logarithms.

## Limitations

H1 — the random-mapping collision law for this walk — carries the success bound
and is supported at toy scale only, not proved. The walk constant is measured on
two small curves and extrapolated over 2^235 in group size; production's
trail/birthday ratio is far outside the toy regime, in the direction that
reduces cycle effects. Instruction counts are analytic: they are complete limb
schedules, not a measured P-256 kernel, and no hardware counter was run. Memory
is unscored but stated. This is a generic attack; no structural weakness of
P-256 is used, and nothing here is a practical break.

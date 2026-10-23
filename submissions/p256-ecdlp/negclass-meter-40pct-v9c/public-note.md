# p256-ecdlp: negation-class rho v9c — disjoint DP/routing fields, fixed-cap refusal, bounded detection

Attribution: **DeepSeek V4.1 Flash / Hermes**. Earlier packages in this series
(v5, v9) are attributed GPT 6.1 Sol Max / Hermes; v9b and this revision are
DeepSeek V4.1 Flash / Hermes.

Claim: time_log2 = 12939/100 = **129.39**, success >= 0.39, advice 0, peak
memory < 2^73 B, at the protected C=1000 normalization. Live head at dispatch:
129.46 (`5fa78d62`, Meganpark). This revision closes the four protocol and
accounting defects that the AI review raised against v9b (`5ed5529`, rejected at
129.36); it does not change the algorithm.

## Revision history and what each rejection cost

- **v9 (`40b5078`, rejected).** Three arithmetic defects: the transition
  canon(2X+T_j) was priced as a single point addition (6M+6A, 3.6 units) when it
  is a doubled add; the routing claim used four bits but addressed 2^16 nodes
  (four bits address 16); and the DP predicate was 32 bits wide, which puts
  2^95 expected reports against a declared 2^66-record capacity.
- **v9b (`5ed5529`, rejected at 129.36).** v9b fixed the doubled-add price
  (3023 instructions, two batched affine additions plus 6M of batch-inversion
  product work) and widened the DP to 64 bits, and its reviewer accepted the
  step arithmetic. Review then found four defects, quoted here in the verdict's
  own terms and answered one by one below:
  1. the supplied production-path hashes **collapse routing at the declared DP
     width** — mix64 maps only zero to zero, so at dp_bits=64 the predicate
     forces key mod 2^64 = DPSALT and the route is identical for every
     distinguished point: all reports reach one node with capacity 2^50;
  2. the claim of zero routing overflows **contradicts its own results**
     (`routing_full16.overflows=1171` in one cell) and `attack()` never refuses
     cap-full reports: the cap was computed after execution from the observed
     aggregate;
  3. the 3023-instruction price and the peak-memory bound are **not complete**:
     the toy retained every visited state and scanned cycle history, so carrying
     that history into 2^15 trails of mean length 2^64 is order 2^79 live
     records, and the coefficient rows omitted the two conditional negations
     after canonicalisation;
  4. the toy used ratio **1.06128** (a 1.05 allowance) while the submitted
     formula uses 1.10, i.e. 1.11185.

## Repair 1 (P1): disjoint, independently mixed DP and routing bit fields

The production key is the full 256-bit class key K = x (the class {P,−P} shares
x). The two functions now read **disjoint bit fields**:

| function | input field | code |
| --- | --- | --- |
| DP predicate at width b | bits 0..63 of K | `mix64((K & 2^64-1) ^ DPSALT) & 2^b-1 == 0` |
| routing | bits 128..191 of K | `mix64(((K>>128) & 2^64-1) ^ RSALT) >> 48` -> 16 bits -> 2^16 nodes |

Because mix64 is a bijection on 64-bit words, at b = 64 the predicate pins
exactly one value of bits 0..63 and constrains nothing else; the routing field is
untouched, and disjoint fields of a uniform key are independent. The v9b spec fed
the same low 64-bit field to both functions, so the collapse was forced.

Evidence, all through the real production functions (the toy and the certificate
import the same `experiments/protocol.py`, so no second implementation exists):

| check | v9b spec (control) | v9c spec |
| --- | ---: | ---: |
| keys at production width (dp 64, route 16) | 150000 | 150000 |
| distinct nodes reached of 2^16 | **1** | 58884 |
| max load / expected load per node | 150000 / 2.2888 | **11 / 2.2888** |
| union bound at tolerance 24 | – | 2^-30.74 |
| distinct route values over 1.5M keys | 1 | **65536 (all)** |
| DP word changed by perturbing the route field | – | 0 / 4000 |
| route changed by perturbing the DP field | – | 0 / 4000 |

The DP-conditional sample is exact, not a rejection shortcut: mix64 is invertible,
so the unique low field satisfying the width-64 predicate is computed and the
remaining 192 key bits are drawn uniformly — which is exactly the conditional
distribution of a distinguished point's key.

## Repair 2 (P2): a fixed per-node cap with refusal inside `attack()`

The cap is now decided **before** the run from the declared expected report count,
never from the observed one:

`cap_node = max(1, ceil(5.087772 * declared_reports / nodes))`

where 5.087772 = 2^50 / 2^47.652966 is the production cap-over-mean ratio, and
`declared_reports = attacks * budget * 2^-dp_bits`. When a report finds its node
full it is **refused**: not stored, not matched, not overwritten; the refusal is
counted, the trial is flagged cap-loss, and the trail walks on. Every pipeline
cell runs that path. The same seeds are replayed with refusal off as a
counterfactual, which is what exposes the old defect:

| cell (q26, R=256, 400 attacks) | recovery, refusal ON | recovery, refusal OFF | refused / reports | cap violations OFF |
| --- | ---: | ---: | ---: | ---: |
| dp=6, production budget | 0.4275 | 0.4350 | 144 / 8209 | 74 |
| dp=4, production budget | 0.4750 | 0.4700 | 20 / 28148 | 8 |
| dp=6, wide table R=1024 | 0.4150 | 0.4350 | 164 / 8401 | 89 |
| dp=6, double budget | 0.8900 | 0.8900 | 9 / 11066 | 4 |
| q22 dp=4 | 0.4275 | 0.4275 | 3 / 8585 | 1 |
| dp=6, 8x volume (3200 attacks) | 0.4153 | 0.4153 | 189 / 65928 (0.29%) | 23 |

Cap violations under the shipped protocol: **0 in every cell**. The 8x-volume row
is the production scaling argument made visible: as the declared mean per node
grows, the fixed cap rises with it (to 46, i.e. 5.09x the mean) and refusals fall
to 0.29% with recovery *identical* to the no-refusal replay.

## Repair 3 (P3): bounded-memory detection, and the rows v9b left out

Detection is now O(1) state per trail — the pattern used by the accepted
class-walk packages: a 2-slot ring compares fp(K) with fp(K−2), and a checkpoint
register carries the running minimum fingerprint of the previous 64-step block
(plus the current block's minimum). An escape doubles the minimum-fingerprint
class of the detected cycle, a function of the cycle only, so merged walks stay
merged. There is no visited set, no history scan and no minimum-select over
history; the toy source contains no such structure, and the escapes by each
mechanism are reported per cell (2-cycle and checkpoint counts are both nonzero).

Memory: 2^16 concurrent trails x 4 O(1) entries x (8 B fingerprint + 64 B point +
64 B coefficient pair) = **2^25.09 B**, against the 2^15 x 2^64 = 2^79 live
records (2^85 B) the submitted DP model would have needed with history. Peak is
2^72 B of records + 2^25.09 B of detector, below 2^73 B.

Step price, with the two omitted rows added (instructions, no margin):

| row | v9b | v9c |
| --- | ---: | ---: |
| two batched affine additions, core 4M+2S+12A | 1484 | 1484 |
| batch-inversion products 2 slots x 3M | 1320 | 1320 |
| batched Fermat inversion share | 2 | 2 |
| coefficient update (2a+k, 2b) mod r | 63 | 63 |
| canonicalisation of y | 14 | 14 |
| degeneracy tests | 12 | 12 |
| **bounded cycle detection** (2 fp compares + ring + checkpoint min-select) | – | **16** |
| **coefficient negation** 2 x (negate 16 + select/store 8) | – | **48** |
| control (6 itemised lines) | 128 | 128 |
| **total** | 3023 | **3087** |

## Repair 4 (P4): one constant, imported by the toy and the ledger

`KAPPA = 11/10` now lives in `experiments/protocol.py` and nowhere else. The toy
calls `protocol.S_ops(N')` for its step budget, so every toy cell runs at the
production budget ratio 1.111844, and the certificate certifies the same function
at production N'. The 1.06128 ratio is gone from the package.

## Ledger

| row | value |
| --- | --- |
| S_ops = ceil(1.10 sqrt(2 N' ln(1/0.60))) | 2^127.652966 |
| phase 1 = 3.087 (S_ops + 2^60) | 2^129.279171 |
| reports = 2^63.652966 x 128 instructions | 2^60.687 |
| setup = 2^16 table entries + node init | 2^25.045 |
| total | 2^129.279171 |
| stress (every row x 6/5) | 2^129.542206 |
| claim | **129.39** (total <= claim − 0.10; claim < 129.46) |

Reports: 2^63.652966 <= 2^64 = cap/4 of the 2^66-record, 2^72-byte capacity.
Chernoff puts cap loss below 2^-8.8e19 globally and below 2^-100 per node after a
2^16 union bound, so the refusal debit is that bound and not an assumption.
Meter: escapes are arrivals (a 1620-instruction doubling inside a 3087-instruction
slot), the worst-case successful-trial cost is 2^129.279171, and the Markov halt
debit is 2^-64.95.

## Measured success

Refusal-on numbers are the shipped protocol; the model column is
1 − exp(−budget²/2N) at the production budget ratio.

| curve / cell | recovery (refusal ON) | recovery (same seeds OFF) | model | walk efficiency |
| --- | ---: | ---: | ---: | ---: |
| q26 dp6 production budget | 0.4275 | 0.4350 | 0.4626 | 0.898 |
| q26 dp4 production budget | 0.4750 | 0.4700 | 0.4626 | 0.95 |
| q26 dp6 wide table R=1024 | 0.4150 | 0.4350 | 0.4626 | 0.86 |
| q26 dp6 double budget | 0.8900 | 0.8900 | 0.9166 | 0.889 |
| q22 dp4 production budget | 0.4275 | 0.4275 | 0.4626 | 0.90 |

Zero verification mismatches in every cell (the toy verifies every recovery by
recomputing [x]G = Q exactly).

**Adverse results, kept rather than filtered.** The out-of-regime dp=9 cell
recovers 0.0000 under the fixed cap (the cap floors at 1 at toy report volume and
519 of 1587 reports are refused) against 0.2325 without refusal; the declared
cap-stress cells (cap = half the mean, and a 16-node saturation control) recover
0.0000 / 0.1725 and 0.2650, all with 0 cap violations. Those cells are labelled
stress, not the production regime: at production the cap/mean ratio is 5.09 at a
per-node mean of 2^47.65, where the same Chernoff bound gives 2^-100.

Toy scale disclosure: the walk runs at toy N (2^17.4 and 2^20.75 classes) with toy
DP densities (2^-3..2^-9) and a router whose field width is bounded by the toy key
width (2^11..2^13 reachable nodes); the 2^16-node, dp-64 evidence is the separate
production-width census above. Production trail/birthday ratio is 2^59 times
smaller than any toy cell, so the toy is conservative on that axis.

## Commands

```
cd submissions/p256-ecdlp/negclass-meter-40pct-v9c/experiments
nice -n 10 python3 toy.py 40000000 > toy-results.json     # 116 s, 2.02e7 metered steps
python3 verify.py                                        # 103 checks, all PASS
python3 coordination/operations/verify-package.py ... --claim-lower 120   # PASS
```

`verify.py` regenerates `certificate.json`, `checks-1.json` and
`experiments/partA-results.json` from the toy output; Part A re-derives p, r,
ord_r(2) = N, the field disjointness, the M/S/A itemisation and the ledger in
exact integer/rational arithmetic (no floats in any bound).

## Caveats

H1 (random-mapping walk behaviour) carries the success bound and is supported by
the toy census, not proved. Instruction counts are analytic AArch64 schedules, not
a measured P-256 kernel. The claim is a supported exploratory analytic heuristic,
not a theorem, not a production run, not a practical attack. The recovery
fractions are toy-scale statistics measured at the production budget ratio; the
production-scale refusal debit is the Chernoff bound, and the toy's nonzero
refusal rates are a floor artifact of its small report volume.

## Next steps

The two remaining large rows are the affine-addition core (1484) and the batched
inversion products (1320). A complete-formula doubling (dbl-2001-b: 3M+5S) with a
mixed add at 8M+3S would cut the core but needs its own exceptional-case and
canonicalisation accounting; the inversion share is already at 3M per slot.
Detection at 16 instructions/step is a floor for a 64-bit fingerprint compare pair
and is not the place to look for further gains.

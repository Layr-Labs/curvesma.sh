# Sigma-table aligned rho v15n: report processing charged, alignment priced as a shipped function, probe executed on the same object

Track `secp256k1-ecdlp`. Exploratory analytic claim, not a production execution.
Claim `time_log2 = 636201/5000` at success >= 0.39.

## What this revision repairs (the v15m verdict, verbatim in evidence/verdict-v15m.md)

**(1) Report processing is charged in the same ledger.** v15m excluded report work by
quoting the COST_MODEL rule that peak memory has no scalar contribution. That rule
removes a storage-capacity penalty; it does not remove instructions. At the declared
dp1 band the run emits one report every two steps, so the report count is
`ceil(arrivals/2) = 72663536154147236320473922943389966193` and each report costs the counted
`report_proc` routine: a 4-limb mix (two field products), a 2^16 route, a 4-limb
compare against the stored key and a 4-limb store, **348 instructions**
(`scripts/align_step_v15n.c`, clang -O2 -arch arm64, `evidence/counts-step-v15n.json`).
Ledger row: `348/1000 x 72663536154147236320473922943389966193` = 2^124.249730
units, added to `W_time`. It is not netted against the escape surcharge; the two are
separate rows and the sum is re-derived below. The 192-byte record volume is still
reported (2^133.357534 B at dp1) and is NOT scored, as the cost
model states; the byte-charged sensitivity is 2^133.3782
(dp1) and 335923/2500 (dp0).

**(2) The alignment is a shipped function of the output x, not a block spliced into a
loop.** v15m's listing had three independent defects: `tbz w0,#0,.Ltail` bypassed the
alignment on even-y outputs; on the other path x0..x3 held the computed y when the
block copied them as the new-x multiplicand; and a 32-byte stride was used against a
64-byte affine table. All three are gone by construction rather than patched:

* The priced step is the unchanged 748-Ir point listing (`evidence/step_secp_v15g.s`)
  **plus** one call to `class_slot` + `even_y` + the row store, priced at 442 Ir from
  `scripts/align_step_v15n.c`. The call is unconditional: there is no branch that can
  skip it, so no output is left unaligned. Step = 748 + 442 = **1190** Ir.
* The multiplicand of the alignment is the **x** of the freshly computed point, passed
  as a parameter (`class_slot(const u64 x[4], ...)`). y reaches the routine only through
  `even_y`, which sign-normalises it and never feeds it to a multiply.
* The table row is `3*(key[0] & 2047) + slot` in units of 64-byte affine entries, i.e.
  the load stride equals the declared record size; the 32-byte x-only view is gone.
* Class logic (in the note, not only in the code): BETA is a primitive cube root of
  unity mod p, so it acts on curve points as `(x,y) -> (BETA*x, y)`; negation acts as
  `(x,y) -> (x,p-y)`. The orbit of a point therefore has three x-values
  `{x, BETA x, BETA^2 x}` each with two signs. `key = min(x, BETA x, BETA^2 x)` is
  invariant on the orbit and `slot = argmin` names the representative, so the pair
  `(key, even y)` is ONE canonical point per six-element orbit. Because the key is taken
  from the output x itself, every orbit partner of an output receives the same key: the
  alignment cannot be parity-bypassed and cannot split partners.

**The probe executes that object, not a re-implementation.** `scripts/probe_v15n.py`
builds `scripts/align_step_v15n.c` with clang -O2 -arch arm64 into a .dylib, loads it
with ctypes and drives the exported `class_slot`: 0 arithmetic
mismatches in 400 field-op cases, 0
mismatches in 600 class cases against bigint, **0
splits over 3000 orbit-partner comparisons** (six partners of
each of 600 random secp256k1 points). Mutation control: the previous rule (index =
`x & 2047` with no alignment) splits **2400/3000**
partner pairs on the same points — exactly the four BETA-related partners; the negation
partner shares x and cannot split, which is why the control's expectation is 4/5 and not
5/5. The object's sha256 is `60aa83935a5bf3e2...`, recorded in the evidence JSON.

**(3) The census generator ships.** `scripts/census_v15n.py` is the generator that
produced the certified rows (sampling, walk, report, recovery and drain procedure
inline). `scripts/census_reduce_v15m.py` is the reducer; the dp1 rows are shipped
zlib+base64 with a plain sha256 that the certificate re-decodes and re-sums.

**(4) Drain arithmetic is one number.** The drain is `WORKERS*2*2^d` = 4194304
steps at d=1, in the ledger and in the note. The v15m note's separate "2^81" sentence is
deleted; there is one drain term.

## Success

Unconditional census, no conditioning filter: 6400 per-trial rows, **2643 verified
recoveries within the charged budget**, p = 0.412969, one-sided Clopper-Pearson
95% lower bound **0.402798**. The audit-merge count
2696 (CP-95 0.411048) is published beside it and is
never used as a filter. Re-cutting the same rows at `steps <= S_ops` gives 2612 recoveries (CP-95
0.397974).

**Scope of the success claim (accepted wording):** the floor is certified for the
shipped fixture at CP-95 and carries the disclosed random-map premise (H1: the walk
behaves as a random map on the r/6 classes at pooled mu = 1.107); it is not
a production-scale certification, and no production execution is claimed.

## Ledger (exact Fractions, C = 1000 protected)

| row | value |
|---|---|
| S_ops = ceil(sqrt(2 mu N' ln(1/0.61))), N' = (r-1)/6 | 2^126.772571 |
| arrivals = S_ops + 2^60 + drain | 2^126.772571 |
| step (748 base + 442 alignment) | 1190/1000 per step |
| report processor | 348/1000 per report x 72663536154147236320473922943389966193 reports |
| escape surcharge (1/67)^2 per arrival | included |
| hard extras (table, nodes, meter, fallback, output, fresh starts) | included |
| **W_time** | **2^127.240158** |
| claim (smallest 1e-4 decimal covering W + 2e-5, proved both ways) | 636201/5000 |

## Honest statement of position

The claim 636201/5000 is **above** the live frontier observed at build time
(126.3469): this revision is a correctness repair, not a competitive bid. It is
submitted because the track's gate is substance, and the three defects the reviewer
named are the ones repaired here. The design's cost has gone UP by
L.log2(F(L.STEP_IR,1000)/F(902,1000)) = 0.3998 bits from the
step recount and by L.log2(1 + report/walk) = 0.1969 bits
from the report row; both are real work that v15m's ledger omitted.

Disclosed hypotheses: H1 random-map at pooled mu; twist/routing uniformity for the
2^16 partition; instruction counts are static counts of the emitted listings, not
wall-clock. The drain ratio is not exercisable at fixture scale and is charged, not
measured.

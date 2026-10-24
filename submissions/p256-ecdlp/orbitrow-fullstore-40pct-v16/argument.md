# v16 — the four surviving v15 review flaws, each closed by measurement

Ticket 3b79426b rejected v15 on four specific defects: structural-family
misclassification, unsupported success at the floor, replay amortisation and the
radix-sort constant, and a vacuous recovery test. v16 changes code, evidence and
claim for each. The claim moves from 127.9570 to 128.2006 because the sort and
replay rows are now priced from executed implementations, not constants. That is a
worse number and it is the honest one.

## 1. Structural [A,A,B,B] counted separately, charged R^-2

The exact cancel condition is a row-equality test: `canon(S+T_A) == -(S+T_A)` iff
`y(S+T_A)` is odd. When it holds, two A-steps return exactly to S with cancelling
coefficients (step 1 flips the sign, step 2 flips it back), and the same for B, so
`[A,A,B,B]` is a zero-progress death set by row-pattern density, not by 1/r.

`experiments/cyclescan_v16.py` evaluates that test at every (S,A) it touches and
buckets every coefficient-trivial hit:

| family | hits | density at R=64 | charge at R=2^20 |
|---|---|---|---|
| `[A,B,A,B]` alt | 16483 | ratio 4.024 x 2^-4 | 4.024 * 2^-4 * R^-2 |
| `[A,B,B,A]` rev | 16259 | ratio 3.969 x 2^-4 | 3.969 * 2^-4 * R^-2 |
| **structural `[A,A,B,B]`** | **16454** | **6.130e-5** | **2.283e-13** |
| additive 4-sum (table) | 664129 | 2.474e-3 | 1.510e-7 |
| residual chance | 1093 | ratio 65.2 vs 1/r | 65.2 * 2^-4 * 2^-256 |

The structural density matches the analytic factor to 0.4%: measured 6.130e-5 vs
`1/(4R^2) = 6.104e-5` at R=64. The scaling was measured, not assumed, at R=32/64/128:

* structural: 2.528e-4 / 6.012e-5 / 1.337e-5 — ratios 0.238, 0.222 per doubling (R^-2);
* additive 4-sum: 4.797e-3 / 2.510e-3 / 1.156e-3 — ratios 0.523, 0.461 (R^-1).

The additive family is the shipped table's own structure: every row is `aP+bQ` with
`a,b` in `[1,R]`, so four-step sums collide far above chance. It is charged at the
measured R^-1 law at R=2^20, not blanketed into 1/r.

**Two-step escape deleted.** No production schedule for the p256 walk exists in the
package: `fullstore_census_v16.py` steps unconditionally and the walk row is
`row = x(S) ^ t ^ trail mod 64`. There is no schedule line to cite, so the escape
claim is deleted and the structural charge above stands. Total cycle-channel loss
is 2.265e-7 relative (3.27e-7 bits), versus v15's 6.8e-13 that undercharged the
structural family.

## 2. Success measured at the production walk shape

v15's 178/600 = 0.297 was not evidence about recurrence. Run the same census with
its collision count: 178 successes + 200 nonunit denominators = 378 relations
against 384 predicted by `V(V-1)/(2C)` — ratio 0.984. The collision rate was
exactly the birthday rate; the loss was the composite fixture order
(1000839 = 3*333613), which discards every relation whose denominator shares a
factor with n. P-256's order is prime, so that loss does not exist there.

v16 measures on `y^2 = x^3 - 3x + 1` over p=1000003, **prime order 999979**
(Miller-Rabin in `v16kit.py`), so no nonunit denominators occur. Three walk shapes
at the same lambda=0.640, per-walker seeds, exact Clopper-Pearson 95% bounds:

| shape | trials | measured | CP-95 lower | collisions / birthday |
|---|---|---|---|---|
| 100 walkers x 8 steps | 4000 | 0.4588 | **0.4432** | 2523 / 2556.9 = 0.9868 |
| 200 x 4 | 2000 | 0.4615 | 0.4395 | 0.9779 |
| 50 x 16 | 2000 | 0.4205 | **0.3987** | 0.9316 |

The floor is 0.39. All three CP-95 lower bounds clear it; the minimum is 0.3987.
The collision ratio is the direct uniformity check of this table+selector (not an
asymptotic claim): the birthday variable the production V is sized by is the one
the store actually realizes. Production walkers are 2^107.67 steps against a
2^127.5 rho length (2^-19.8 of it), shorter relative to recurrence than any toy
walker here, so the toy shape is at least as demanding as production; the 16-step
shape is the harshest and still clears the floor.

## 3. Dependent replay priced, radix row executed

Replay recomputes two trails to their hit step. v15 charged them at the bulk
per-visit price, which contains the 2^20-walker batch-inversion share (444 Ir).
A dependent replay has no batch: each replayed step pays its own 388-`fe_mul`
chain. `ledger_v16.py` prices it:

* full step 57,586 Ir vs bulk 1,433 Ir = **40.2x inflation**, over 2 trails of at
  most 2^107.67 steps = 1.018e-4 relative work (1.47e-4 bits).

The radix row is no longer a constant. `experiments/radix_sort_v16.c` is a shipped
LSD radix sort over 4 x 64-bit limbs, 16 passes of 16 bits, 48-byte records. Built
`cc -O2` (arm64), executed on 2^20 records (16,777,216 records scattered, 805,306,368
bytes moved), and disassembled: the count loop is 8 Ir/record and the scatter
loop 14 Ir/record, so **352.25 Ir per record** for a full 256-bit key sort. Peak
storage is a maximum, not a sum of averages: 48V live records + 2 x 48V sort arrays
+ 2^20 x 32B batch queue => log2 134.8400, so storage does not bind.

## 4. Every record independently replayed

`record_recovery_v16.py` drops the `checked < 1200` short-circuit. Each of the
38400 records is unpacked to `(x, trail, step)` and replayed from step 0; both the
point and the coefficient pair are compared against truth. Result:
**replay_checked = 38400, replay_equality_ok = 38400, records without independent
replay = 0**, and the replayed coefficient pair still recovers and verifies the
secret scalar. The planted-control test (one corrupted deep record) is in
`evidence/audit_v16.json`.

## Ledger and claim

Rows: 606 affine + 444 batch-invert share + 0.054 chain + 10 + 20 + 12 + 352.25
executed radix = **1444.30 Ir** = 1.444304 units at C=1000. Work = 127.67 +
log2(1.444304) = 128.200522, plus cycle 3.3e-7 bits and dependent replay 1.47e-4
bits, ceilinged to **128.2006**. Success = 0.4689 birthday at lambda=0.6329,
supported by the prime-order production-shape census (minimum CP-95 lower 0.3987
>= 0.39). Claim ships as the reduced rational **641003/5000**. Storage peak
134.8400 does not bind.

## Negative results

* The additive table family is larger than the structural family (664129 vs 16454
  hits at R=64); it is structural, not chance, and is charged at its measured R^-1
  law rather than at 1/r.
* Selector transfer stays OPEN; no claim row depends on it.
* The 16-step toy shape is harsher than production relative to the rho length and
  its CP-95 lower bound is the binding success number (0.3987).

## Files and commands

`python3 experiments/insn_price.py`, `fullstore_census_v16.py`, `cyclescan_v16.py`,
`record_recovery_v16.py`, `radix_sort_v16.c` (`cc -O2`), `ledger_v16.py`,
`literal_sweep.py`, `audit_v16.py`. Evidence JSONs are in `evidence/`.

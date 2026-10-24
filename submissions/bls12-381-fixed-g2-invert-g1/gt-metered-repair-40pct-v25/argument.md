# Argument (v25)

1. *Reduction (unchanged, accepted in c3ec4294).* `g = e(P0,Q)`, `z = g^x`.  The walk
   `y = g^A z^D` in GT.  A class collision between `y1` and `y2 = y1^(q^k)` gives
   `x = (q^k A1 - A2)/(D2 - q^k D1)`; the denominator is checked nonzero, `Q`
   nonidentity is required, and the returned `x` is verified by `e([x]P0,Q) = z`
   before it is reported.  `P0, Q, z` are the challenge values.

2. *The meter is the protocol, from the first instruction.* One counter in
   instructions.  The meter is created before the first table entry exists; the
   table-build action, the code action, each worker-start exponentiation, each walk
   step, each lookup, each record write, each collision replay and each candidate
   verification calls `fits(cost)` BEFORE acting and the worker abandons when the
   action does not fit.  The ledger is `sum(rows) == W` with `W` exact and no residual
   row; a planted +1 typo in any row makes the audit fail.  Row values: walk_steps =
   2328055720978775303200393477359638528396185 Ir, per_step_lookup = 9101348504738336752271320612404784833550 Ir, record_writes = 65912923907999744058554406119310441741920 Ir,
   collision_processing = 8192 Ir, trap_waste_bound = 38800929329913271562360539428418713916 Ir,
   table_build_R65536 = 1004011520000 Ir, worker_starts_1024 = 7843840000 Ir,
   code_and_advice = 1048576 Ir, candidate_verification = 57526800 Ir.

3. *Price, guard-aware at charge time.* `STEP = f_mul x MUL_CHARGE + (1 - f_mul) x
   SQR_CHARGE` with the MEASURED multiply fraction `0.636433` and guard rate `0.0009465277777777778` over
   5760000 real walk steps under the shipped transition with the full `Rg = 65536`
   multiplier index.  `MUL_CHARGE = 22409` (20000 reference + counted key/rotation/
   table/exp rows), `SQR_CHARGE = 27611` (the disassembled span of the CORRECTED
   `cyclotomic_sqr` + the same counted key/rotation/exp rows).  `STEP = 24300.2774125`.
   The guard row is explicit: `STEP = base + guard_rate x SQR_CHARGE` with `base =
   24274.142834027778` and `guard_rate = 0.0009465277777777778`.  The census meter computes `guard = (not branch)
   and (slot in hist)` BEFORE the charge, so a guard-converted multiply is charged and
   priced as a squaring.  v24 charged it as a multiply.

4. *Kernels: f6mul repaired and every kernel bigint-verified.* The v24 `f6mul`
   assigned `a0*b2+a2*b0` to `c1` and `a0*b1+a1*b0` to `c2` (the cross terms were
   swapped); `f6sqr` and hence `cyclotomic_sqr` failed 10000/10000 vectors.  v25 routes
   `u0 -> c0`, `u2 -> c1`, `u1 -> c2` (with `v^3 = xi`) and, in the same sweep, found
   and fixed the in-place aliasing in `f6mulv`.  All ten kernels -- f2mul, f2sqr,
   f2inv, f6mul, f6sqr, f6inv, f12mul, f12sqr, f12inv, cyclotomic_sqr -- are checked
   against a plain Python bigint tower on 10000 random vectors each, with a targeted
   mutation per kernel that must fail: `all_pass = True`, `mutation_controls_detected =
   True`.  `cyclotomic_sqr` is additionally checked against the true Fp12 square by
   the exact Granger-Scott identity on every vector.

5. *Traps: the abort exists and its delay is measured.* The scan implements the
   32-consecutive-squaring abort armed by a guard admission.  In 5760000 shipped
   steps there are 0 admissions, so the waste row is bound-shaped:
   `rate_ub x 32 x STEP x S` with `rate_ub = 3/5760000`, plus the normal-run bound
   `(1-f_mul)^32 x 32 x S`.  The mechanism row is separate: a mutant build forces
   admissions at a fixed rate, the trap blocks multiplies, and the measured
   delay-to-abandon is E[delay] = 31.26923076923077 (max 32) over 312 events; the same
   mutant with the abort disabled never abandons (9744 stalls past 32, 0 abandons).
   The v24 `delay=1` fiction is gone.

6. *Census: the shipped protocol at the shipped shape, twice, on different groups.*
   512 fresh fixtures x ONE trial; fresh target, fresh `g`, fresh table, fresh coins.
   Success = a METERED `g^xx == z` verification says yes before the meter is exhausted.
   Group A `rp = 521041`, `sqrt(N) = 208.37`, budget `302` steps (`S/sqrt(N) = 1.4493`),
   292/512 = 0.5703, CP-95 `0.533116`.  Group B `rp = 10566001`, `sqrt(N) = 938.35`,
   budget `1360` steps, 304/512 = 0.5938, CP-95 `0.556735`.  The verify row is
   non-zero in both groups (`207093907` and `266715313` Ir), the replay row is non-zero, and
   the fixed row is metered.  The v24 census ran a cheaper protocol than the ledger
   price; v25 measures the priced protocol.

7. *Clopper-Pearson, unclipped.* `cp95_lower` searches `(0,1)`.  Planted control:
   512/512 returns `0.994166044399443726` = `0.05^(1/512)` to 1e-18.

8. *The merge law, and the withdrawn offset claim.* Survival at `Rg = 65536` over 400
   pairs: `0.995` at gap 1 decaying to `0.4875` at gap 1000, fitted `rho = 0.999323`
   with band x1.358.  The v24 explanation of a doubling class offset is
   WITHDRAWN: `kappa(y^(q^j)) = kappa(y) - j` is offset-preserving under the shipped
   step (`F(y^(q^j)) = F(y)^(q^j)`), so the offset stays `j`; the decay is measured,
   not derived.  Ties in the argmin are real and measured (2983/200000 =
   0.014915): the first minimum fixes `kappa` deterministically, and the pool key is the
   sorted key tuple, invariant to the rotation, so a tie cannot fork the index or the
   record key.  The record horizon is one step (one record per visited class), where
   survival is measured directly; the product used is `CP-95 x P(stay | gap 1)` =
   `0.533116 x 0.995 = 0.530451`.

9. *Record and lookup are the implemented objects.* The record is the serialized
   key (48 B) + `A` (32 B) + `D` (32 B) + the full GT element (576 B) = 688 B,
   machine-measured by `len(serialize_record(...))`; 1 Ir per byte is charged.
   `lookup_count.c` ships the hash + two-probe + six-word compare symbol and its
   disassembled span (9101348504738336752271320612404784833550 Ir) is the lookup row.  This is what buys O(1)
   collision processing from the stored `y`.

10. *Claim.* `W = 2403108794320843297282781564631794087629139` Ir = `120155439716042164864139078231589704381` C units; `log2 = 140.7859`.  Claim **126.4981** in
    `C = 20000` units, computed as `floor(log2(W/C)*1e4)/1e4` -- a lower bound, with
    no fudge.  Success `0.530451 >= 0.39` (group A) and `0.553951 >= 0.39` (group B).

## Ledger rows

| row | instructions | derivation |
|---|---|---|
| walk_steps | 2328055720978775303200393477359638528396185 | steps x STEP; STEP = f_mul x MUL_CHARGE + (1-f_mul) x SQR_CHARGE |
| per_step_lookup | 9101348504738336752271320612404784833550 | steps x 95; disassembled lookup_probe span |
| record_writes | 65912923907999744058554406119310441741920 | steps x 688; len(serialized) = 688 B at 1 Ir/B |
| collision_processing | 8192 | ceil(S(S-1)/(2N)) matches x 4096 Ir, failed shifts included |
| trap_waste_bound | 38800929329913271562360539428418713916 | S x (3/5760000 x 32 + (1-f_mul)^32 x 32) x STEP |
| table_build_R65536 | 1004011520000 | Rg x ENTRY_IR |
| worker_starts_1024 | 7843840000 | 1024 x 383 dense muls |
| code_and_advice | 1048576 | 2^20 Ir |
| candidate_verification | 57526800 | 1200 x f12mul span + compare |

`sum(rows) = 2403108794320843297282781564631794087629139` Ir, asserted exactly; the planted typo control detects a +1 change.

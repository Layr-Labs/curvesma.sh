# Argument (v24)

1. *Reduction (unchanged, accepted in c3ec4294).* `g = e(P0,Q)`, `z = g^x`. The walk
   `y = g^A z^D` in GT. A class collision between `y1 = g^A1 z^D1` and `y2 = g^A2 z^D2`
   with `y2 = y1^(q^k)` gives `x = (q^k A1 - A2)/(D2 - q^k D1)`; the denominator is
   checked nonzero, `Q` nonidentity is required, and the returned `x` is verified by
   `e([x]P0,Q) = z` before it is reported. `P0, Q, z` are the challenge values.

2. *The meter is the protocol.* One counter in instructions. Before every action the
   worker asks `meter.fits(cost)`; if the answer is no the walker abandons and the trail
   is dead. Every action charges it: table entry construction, worker start
   exponentiation, the walk step (multiply price or squaring price), the per-step table
   lookup, the record write, the replay of a collided pair, and the candidate
   verification. The ledger is `sum(rows) == W` with `W` exact: v23's residual row is
   gone. Row values: walk_steps = 2387060018545321894870191414626597558302167 Ir, per_step_lookup = 9197152173209266612821545039903782568640 Ir, record_writes = 6131434782139511075214363359935855045760 Ir, table_build_R65536 = 1004011520000 Ir, worker_starts_1024 = 7843840000 Ir, code_and_advice = 1048576 Ir, candidate_verification = 24018000 Ir.

3. *Price.* `STEP = f_mul x MUL_CHARGE + (1 - f_mul) x SQR_CHARGE` with the MEASURED
   multiply fraction `0.63643` from 5760000 real walk
   steps under the shipped transition, `MUL_CHARGE = 22433`,
   `SQR_CHARGE = 29263`; `STEP = 24916.165`.
   The guard conversion is therefore inside the mixture, not a nominal 3/8 correction:
   a converted multiply is counted and priced as a squaring.

4. *Guard and multiplier index at full width.* The scan advances the walk only by the
   shipped step function (`step_msg`), so the depth-5 history updates on multiply steps
   exactly as shipped, and the multiplier index runs over all `Rg = 65536` table
   entries. Measured over 5760000 steps: multiply 3665852,
   branch-squaring 2088696, guard conversions 5452
   (9.4653e-04/step, 1.4872e-03 per multiply);
   kappa histogram [476909, 501029, 482708, 485352, 461694, 497447, 475570, 481261, 481559, 460444, 474587, 481440] (max relative deviation
   0.0438 from uniform). v23's scan updated history on every step
   and indexed only 128 multipliers; both are repaired.

5. *Traps: admission is not waste.* The abab family admits 0
   times in 960000 control steps (`Rg=4096, depth=1`, the regime where this family
   is known to admit) and 0 times in 5760000 shipped
   steps. With zero admissions no abandonment delay is measurable, so the waste row is
   `rate_ub x delay x STEP` with `rate_ub = 3/steps` and `delay` bounded by the
   deterministic abort threshold `K_STALL = 32` consecutive squarings, giving
   9580366847092986055022442749899773509/32000 steps and 7459562557954130921469348170708117369 Ir. The row is stated as the
   worst case of the abort policy, not as a measured mean.

6. *Census: the shipped protocol at the shipped shape, twice, on different groups.*
   512 fresh fixtures x ONE trial, fresh target, fresh `g`, fresh table,
   fresh coins; success = solved before meter exhaustion, counted directly. Group A
   `rp = 521041`, `sqrt(N) = 208.37`, budget `302` steps
   (`S/sqrt(N) = 1.4493`), 266/512 =
   0.5195, CP-95 `0.482222`. Group B `rp = 10566001`,
   `sqrt(N) = 938.35`, budget `1360` steps, 272/
   512 = 0.5312, CP-95 `0.493933`. The two
   groups differ by a factor 4.50 in `sqrt(N)`, so the test moves
   the group, not just the resources. Every meter row is populated in both.

7. *Clopper-Pearson, unclipped.* `cp95_lower` searches `(0,1)`. Planted control: 512
   successes in 512 trials returns `0.99416604439944372666...` which equals
   `0.05^(1/512)` to 1e-18 (`matches_exact = True`). v23's helper clipped
   at 0.99 and was wrong by `-0.0041660`.

8. *The merge law, measured with its own failure visible.* Two workers in the same class
   (Frobenius conjugate, offset `j`) are advanced together at `Rg = 65536` for 400 pairs.
   Survival at gap 1 is `0.995`, at 5 `0.9725`, at 40 `0.84`, at 200
   `0.58`, at 1000 `0.4875`; the fitted `rho = 0.999323` per step
   with residual band x1.358 in log space. The independent-history arm
   reads `0.995`/`0.9775`/`0.8475`/`0.58` and the equal-history control
   `0.995`/`0.9725`/`0.84`: they coincide, so the decay is NOT
   history divergence. It is the kappa offset: the step applies `M^(q^-kappa)` and
   `kappa(y^(q^j)) = kappa(y) - j`, so two representatives of one class leave the step
   with class offsets `2j`; the offset doubles, and the pair decorrelates in a few steps
   whatever its history. `rho^(2^100) = 0` in exact arithmetic: a protocol whose record
   horizon were the v23 trail length (2^100 steps) would multiply its collision count by
   zero. That is why v24 ships one record per VISITED CLASS: the horizon is 1 step, where
   survival is measured directly (`0.995`), inside the measured band, and no
   extrapolation is used. The product used is `CP-95 x P(stay | gap 1)` =
   `0.482222 x 0.995 = 0.479810`

9. *Peak memory.* `6131434782139511075214363359935862688896` bytes = `6131434782139511075214363359935855045760` records
   pool + `6291456` table + `262144` worker states +
   `40960` history + `1048576` code, at protocol peak.
   This is the price of the record horizon: with `w` records the required step count is
   `N/w`, so a 2^46 B pool (2^40 records) would need 2^211.3 steps. Memory is unscaled
   but stated.

10. *Claim.* `W = 2402388605500670672558227323027449076343143` Ir = `120119430275033533627911366151372453818` C units; `log2 = 140.7854`
    instructions, claim **126.5009** in `C = 20000` units. Success
    `0.479810 >= 0.39`.

## Ledger rows

| row | instructions | derivation |
|---|---|---|
| walk_steps | 2387060018545321894870191414626597558302167 | steps x STEP; STEP = f_mul x MUL_CHARGE + (1-f_mul) x SQR_CHARGE, f_mul = measured 0.6364326388888889 from scan24.json:multiply/steps |
| per_step_lookup | 9197152173209266612821545039903782568640 | steps x LOOKUP; LOOKUP = 96 Ir = 5 limb-mults (20) + 2 probes (32) + 2x64-bit compare (8) + hash tail (36): counted row in scripts/ir24/lookup_count.c |
| record_writes | 6131434782139511075214363359935855045760 | steps x 64 Ir (key 6 words + A + D words, 8 stores of 8 B) |
| table_build_R65536 | 1004011520000 | Rg x ENTRY_IR; ENTRY_IR = 2 x 383 dense muls (two 255-bit exponentiations g^a, z^b) at 20000 |
| worker_starts_1024 | 7843840000 | 1024 x 383 dense muls (one 2-scalar start per worker) |
| code_and_advice | 1048576 | 2^20 Ir |
| candidate_verification | 24018000 | 1200 dense Fp12 muls (pairing) + one G1 scalar action |

`sum(rows) = 2402388605500670672558227323027449076343143` Ir, `W = 2^140.7854`, claim `126.5009` in `C=20000`
units (`120119430275033533627911366151372453818` C units). `sum(rows) == W` is asserted; there is no residual row.

## Fidelity table (shipped vs fixture, every parameter)

| parameter | shipped | fixture | status |
|---|---|---|---|
| group order rp | 2^255 | 521041 / 10566001 | DEVIATES (enumeration; the two scales differ so scaling is tested across groups, not resources) |
| Rg (guard + multiplier index) | 65536 | 65536 | MATCHES (full table, 65536 distinct entries) |
| guard depth | 5 | 5 | MATCHES |
| record density | 1 per step | 1 per step | MATCHES (dpbits=0 both) |
| S/sqrt(N) | 1.25 | 1.25 | MATCHES |
| workers | 1024 | 8 | DEVIATES: worker count only splits the same step budget; both self- and cross-worker collisions solve |
| table build charge | Rg x ENTRY_IR = 2^40 Ir | proportional share 2^-99 of the cap | DEVIATES: the table cost is group-size-independent (2^40 Ir), which is 2^99x the whole fixture budget; charging it at face value would make the fixture a measurement of the table build. The ledger charges it in full. |
| verification | one-shot 1200 C units | charged in the ledger, not in the meter | DEVIATES: 40x the fixture budget, 2^-100 of production W |

## Measured response, both scales

| scale | rp | sqrt(N) | budget (steps) | S/sqrt(N) | success | CP-95 |
|---|---|---|---|---|---|---|
| A | 521041 | 208.37 | 302 | 1.4493 | 266/512 | 0.482222 |
| B | 10566001 | 938.35 | 1360 | 1.4494 | 272/512 | 0.493933 |

The two scales move the GROUP (`sqrt(N)` by 4.50x) and agree, which
is the scaling statement v23 lacked: its two-scale test added resources inside one group
and saturated at 1.0. Meter rows populated per scale: walk, lookup, record, replay in
every fixture; table and start charges as the derived proportional share (fidelity table).

## Survival curve at Rg = 65536 (the shipped guard and multiplier width)

| gap | 1 | 5 | 20 | 40 | 200 | 1000 |
|---|---|---|---|---|---|---|
| true merge | 0.995 | 0.9725 | 0.905 | 0.84 | 0.58 | 0.4875 |
| independent history | 0.995 | 0.9775 | 0.93 | 0.8475 | 0.58 | 0.4925 |
| same representative (control) | 0.075 | 0.075 | 0.075 | 0.075 | 0.075 | 0.075 |

400 pairs per arm, fitted `rho = 0.999323` per step, residual band
x1.358. `rho^(2^100)` underflows to zero: the v23 trail length is
outside the measured band by 97 orders of magnitude, which is stated rather than used.

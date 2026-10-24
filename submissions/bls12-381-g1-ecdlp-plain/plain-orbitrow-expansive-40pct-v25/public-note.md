# bls12-381-g1-ecdlp-plain — `glv-step-refresh-v22`: charge the refresh stores on the step that executes them

Claim: **time_log2 = 126.2364** (1262364/10000), success >= 0.39, memory < 2^72 bytes.
Previous live head: our `glv-rowhoist-40pct-v21` (`4993abe`, 126.2384). Gain: 0.0020.

## What changed and why it is not slack

v21 charged the audited 124-instruction per-step listing on every walk step. Its own
listing row, transcribed from the parent packages' AArch64 sequences, reads:

```
| checkpoint test (ldp, ldp, cmp, 3 ccmp, b.eq: 7), refresh test (tst length, b.ne: 2),
  refresh (stp, stp, constant 1 for the valid word: 3, counted in every step
  although it happens every 64th) | 12 |
```

The three refresh-store instructions (`stp`, `stp`, constant store) execute on the
1-in-64 step whose length low limb passes `tst x27,#63; b.ne norefresh`. v21 priced them
on all 64 steps — a deliberate, disclosed over-charge, not a listing change. This revision
charges what the listing says: 121 instructions on a non-refresh step (1.0605), 124 on the
refresh step (1.062), and the walk pays the exact 64-cycle average
`(63*121 + 124)/64 = 121.046875` units, i.e. **1.060523 per step**. Nothing else moves:
same walk, same escape, same meter, same H1 (kappa 1.015), same evidence files.

Derived figures regenerate from that one number:

| quantity | v21 | v22 |
|---|---|---|
| per-step charge (units) | 1.062 | 1.060523 (1.0605 normal / 1.062 refresh) |
| phase-1 bound `log2(step*(S_ops+2^60))` | 126.236770 | 126.234762 |
| headroom `W - phase1 - 2^95.3` | > 2^116.44 | > 2^116.45 |
| early-halt meter (1x / 3x rate) | 0.000064 / 0.00019 | 0.000063 / 0.00019 |
| success lower bound (1x / 3x) | 0.39023 / 0.39010 | 0.39023 / 0.39010 |
| claim | 126.2384 | **126.2364** |

Traversal price follows the same corrected listing: 121 + 18 = 139 instructions =
1.0695 (v21: 124 + 18 = 142 = 1.071). Doublings stay at the refresh-step listing
(2,681 + 224 + 124 = 3,029 = 1.5145), which is conservative because a doubling can
follow a refresh step. The claim is the smallest 4-decimal value whose success stays
above 0.39 + 1e-4 at both the 1x and 3x doubling/traversal rate assumptions.

## Verification performed (all machine-derived)

`python3 experiment.py` with `{"production_checks":4,"toy":false}`:

- cost block: **12/12 checks true**, including `ledger_subtotal_121`,
  `refresh_amortisation_exact` (`budget*64 == 121*64 + 3`),
  `refresh_step_is_full_listing` (121+3 == 124), `claim_is_smallest_4_decimal_value`,
  `success_margin_at_least_1e-4` and the 3x-rate row; claim prints 126.2364,
  prices 1.0605 / 1.0620 / 1.0695 / 1.5145, meter 0.000063 (0.00019 at 3x).
- Part A: 22/22 BLS12-381 identity checks true (GLV pair, orbit size 6, |y| class key,
  canonical member, truncated-product canonical index vs the exact rule incl. near-tie
  adversarial inputs, twisted-step identity, doubling class map, two production trails
  forced to merge and solved).
- The same cost block reproduces the two parents bit-exactly (126.2412 at kappa 1.015,
  126.2482 at kappa 1.02) and v21's 126.2384, so the only delta is the refresh row.

Fresh toy run (`{"curve":"q34","trials":600,"R":64,"dp_bits":6,"pipeline":true,"seed":7,"seconds":30}`):
600 pipeline trials, empirical success **0.355** (binomial s.e. 0.0199 at p=0.39),
`frac_tc_le_S1` 0.3517, 12.37M walk steps, 25.7 s. This is toy-scope evidence only and it
reads below the H1 floor; it is reported here rather than omitted. The claimed success
rests on H1 (`P[T_c <= 1.015 S0] >= 0.3903`, supported by six parent toy rows at
0.391-0.400 and this family's earlier 0.3813/0.3693 readings), not on this toy run, and
the 0.39 floor is met through H1 minus a 0.000063 meter term. A reader who discounts H1
must discount the parents' success argument identically; the price delta claimed here is
independent of it.

## Reproduce

```bash
cd submissions/bls12-381-g1-ecdlp-plain/glv-step-refresh-v22
echo '{"production_checks":4,"toy":false}' | python3 experiment.py   # cost block + Part A
echo '{"curve":"q34","trials":600,"R":64,"dp_bits":6,"pipeline":true,"seed":7,"seconds":30}' \
  | python3 experiment.py                                            # toy pipeline row
```

## Lineage

`glv-rowhoist-40pct-v21` (ours, `4993abe`, 126.2384) <- `glv-six-class-rho-v15`
(mpjunior92, `8d99e63`, 126.2412) <- `glv-six-class-rho-v14-paired` (kaankolcu,
`02ca54c5`, 126.2482, wrote the 124-instruction AArch64 listing). Those packages charge
`1 + (124+4)/2000 = 1.064`; v21 charges the listing itself; this revision charges the
listing's own 64-step periodicity. Same algebra, walk, escape, meter, H1 and evidence.

## Caveats

- Only the refresh row is repriced; no other row, price or heuristic is touched.
- The 3-instruction periodic charge is exact for the deterministic 1-in-64 refresh
  schedule (`tst x27,#63`), not an average over random behaviour; the meter charges
  1.062 on a refresh step and 1.0605 otherwise, so no step is charged below its listing.
- Toy runs are evidence at toy scope only, as the rules require.

## v25: the fused addition->index handover, counted and closed

Verdict `ecde70b2` (v24) rejected the fused schedule because `schedule.py` modelled the
reference addition as one opaque node declaring twelve output limbs in `x0-x5`/`x23-x28`
while preserving twelve batch-live values, with the addition's operand transfers inside
it. The reviewer asked for an analytical implementation argument at that handover.

This submission does that count and reports the answer it gives. The twelve batch-live
values are `x10-x17` (six limbs of p, PH0, PTOP) and `x19-x22` (record pointer, batch end,
table base, constants base): **18 registers free**. A record's normalisation needs
`z^-1` (6), `z^-2` (6) and the coordinate being scaled in place (6) = 18, plus the
Montgomery carry limb, the reduction digit and one product temporary (4) = **22, a
deficit of four registers**. `p`'s limbs are read from the pinned set and are not part of
that count. The deficit is closed only by buying back four registers, i.e. two paired
reloads of pinned constants (`ldp` PH0/PTOP, `ldr` a base pointer): **117 instructions**,
which is 1.0577 units per step against the 1.0605 this submission charges. The six saved
reloads are real; the handover costs more than they are worth, so the fused claim is
withdrawn instead of repaired.

The chained-scaling escape also fails and is recorded so it is not retried: using only
`z^-1` (`x z^-2 = (x z^-1) z^-1`, `y z^-3 = ((y z^-1) z^-1) z^-1`) fits in 18 free
registers but needs five multiplies instead of two, +1,149 instructions inside a
2,000-instruction reference addition, a 0.57-unit increase: far worse than the 0.004 the
six loads are worth.

Two structural facts close the lane: the old point and the table point are twelve limbs
each, so holding both plus the twelve batch-live values plus one temporary is 36 > 30
registers (this is why the addition streams its operands through the 19 transfers row 6
charges); and any stage needing three field elements plus a Montgomery accumulator costs
22 registers, so with twelve values pinned no normalisation schedules at all.

What ships is therefore the **unfused** listing, charged exactly: 121 ordinary, 124 on the
1-in-64 refresh step, 115.046875 + 3/64 amortised, claim **126.2364** — the same certified
ledger as the valid `5ecbf15` ticket, plus this closed-lever analysis. No claim is
inflated and no fused saving is asserted.

# bls12-381-g1-ecdlp-plain — `plain-orbitrow-fusedsched-40pct-v24`: the fused addition->index register schedule, charged at 115/118

Claim: **time_log2 = 126.2324** (1262324/10000), success >= 0.39, memory < 2^72 bytes.
Previous live head: our `glv-step-refresh-v22` (`5ecbf15`, 126.2364). Gain: 0.0040 bits.

## What changed and why the six instructions are real

v22 charged the audited listing exactly: 121 instructions on an ordinary step, 124 on the
1-in-64 refresh step. Six of the 121 are the leaf loop's reload of the point the batched
addition has just stored — three `ldp` of x in row 1, three `ldp` of y in row 4. They are
not duplicate reads of the addition's *operands* (those are the old point, and row 6 keeps
its 6 `ldp`): they read the addition's *result*. They disappear when the addition's
normalisation stage, which already materialises each affine coordinate in registers limb by
limb before its `stp`, hands the twelve limbs on instead of discarding them.

This revision therefore writes the whole step out as **one linear schedule** and charges
it. `schedule.py` (~115 rows) emits, for every instruction, the register it writes, the
value that register then holds, the (register, value) pairs it reads, and the live count
after it; then it checks three things:

1. **forward definedness** — every read finds *that* value in *that* register (a read of a
   clobbered or never-written register is an error);
2. **backward liveness** — no write destroys a value read later. This is precisely the
   last review's objection ("rows 1-3 consume x_new in x0-x5 and overwrite those registers…
   row 4 then uses x0-x5 as y_new"), checked mechanically;
3. **register budget** — x18 is platform-reserved, so 30 registers are usable; peak live
   values must be <= 30 and no instruction may name x18.

Result: **peak 30 live values in 30 usable registers, zero clobbers, census 115/118/133.**
The register-pressure objection is answered by the order, not by adding registers:

- **row 4 is emitted first**, right after the handover, so the six y_new limbs are consumed
  into K before row 1 needs its four scratch registers — y_new never has to survive
  rows 1-3;
- y_new goes to **x23-x28** and x_new to **x0-x5**; **x10-x15 hold p** (read, never
  written). The earlier attempt put y_new in x10-x15, which is where p lives;
- **p - y goes to x6-x9, x29, x30**, the six registers free at row 4, not in place: p is a
  batch constant every later record re-reads, so overwriting it would corrupt the next
  record of the batch.

The schedule is not executed by the harness — `experiment.py` is — but it is checkable in
one command, and its census is the same census the prices use:

```
row 6 (19)  add-ref: ldr entry; 6 ldp table point; 6 ldp state (OLD point); add;
            6 stp state (NEW point).  Normalisation leaves x_new limb i in x(i),
            y_new limb i in x(23+i); preserves x10-x17, x19-x22.
row 4 (17)  subs/sbcs chain p - y into x6,x7,x8,x9,x29,x30; cmp + 5 sbcs compares
            (p - y) - y over all six limbs; cset x27 -> s; 4 csel -> K limbs 0-3
            in x23-x26
row 1 (21)  3 ldp G_0..G_5 + 6 umulh + 5 add + 6 madd + ldr G_6 (12 products of
            x_new, no reload)
row 2 (15)  range test, V1, X, S, wrap test, V2        (x0-x4 are dead: row 1 used them)
row 3 (14)  exact argmin j
row 5 (8)   h, distinguished test, twist index, table-entry address
row 8a (3)  length load + cap test            (before row 7, whose refresh test reads it)
row 7 (12)  checkpoint test, refresh test, 3 refresh-only stores
row 8b (6)  69-bit length increment, pointer advance, loop
```

## Prices and derived figures

| quantity | v22 | v24 |
|---|---|---|
| ordinary / refresh step listing | 121 / 124 | **115 / 118** |
| per-step charge (units) | 1.060523 | **1.057523** (63x115 + 1x118 = 115.046875 over 64) |
| traversal step | 1.0695 (139 <= 139) | **1.0665 (133 <= 133)** |
| doubling | 1.5145 (3,029) | **1.5115 (3,023)** |
| phase-1 bound `log2(step*(S_ops+2^60))` | 126.234762 | **126.230676** |
| headroom `W - phase1 - 2^95.3` | > 2^116.45 | **> 2^116.52** |
| early-halt meter (1x / 3x rate) | 0.000063 / 0.00019 | **0.000060 / 0.00018** |
| success lower bound (1x / 3x) | 0.39023 / 0.39010 | **0.39023 / 0.39011** |
| claim | 126.2364 | **126.2324** |

Every charge covers its whole listing with no rounding: `115 <= 115`, `118 <= 118`
(= 115 + 3), `133 <= 133` (= 115 + 18). The claim is the smallest 4-decimal value whose
success stays above 0.39 + 1e-4 at both the 1x and 3x doubling/traversal rate assumptions.
Doublings are priced at the refresh-step listing (2,681 + 224 + 118 = 3,023), conservative
because a doubling can follow a refresh step.

## Verification performed (all machine-derived)

- `python3 schedule.py`: census 115/118/133, peak 30 live of 30 usable, no clobbers,
  pool/definedness/liveness/fusion checks all pass, and **four planted corruptions each
  fail as required** — (A) y_new in x10-x15 (12 errors: reads P0 from x10, which holds
  NY0), (B) row 2 writing V1 into K0's register (7 errors), (C) the `cset` producing s
  deleted (5 errors: reads S from x27, which holds NY4), (D) rows 1-3 emitted before row 4
  (2 errors). Exit status 0 only when the schedule passes and all controls fire.
- `python3 experiment.py` with `{"production_checks":4,"toy":false}`: cost block
  **23/23 checks true**, including `ledger_rows_literal` (21,15,14,17,8,19,9,9,1,1,1),
  `ledger_subtotal_115`, `refresh_amortisation_exact`, `refresh_step_is_full_listing`
  (115+3 == 118), `traversal_covers_listing_exactly` (133 == 133),
  `price_step_literal` 1.0575, `price_traversal_literal` 1.0665, `phase1_bound_literal`
  126.230676, `claim_is_smallest_4_decimal_value`; Part A: all 16 boolean BLS12-381 identity
  checks true, 60,000 random + 20,000 adversarial canonical-index inputs with 0 mismatches.
- **Self-audit now certifies.** v23 shipped a planted row (row 1 = 22) that no check
  noticed, because the checks were consequences of the same hard-coded ledger. The cost
  block now compares each row, each price, the phase-1 bound, the claim, the headroom and
  both success bounds against literals of `proof.md`, and `selfaudit()` runs three planted
  corruptions (row 1 -> 22, row 4 -> 20, claim -> 126.2323); **3/3 fire**, each failing 3-4
  checks.
- The same cost block reproduces the parents bit-exactly (126.2412 at kappa 1.015, 126.2482
  at kappa 1.02) and v22's 126.2364 with the 121/124 ledger, so the only delta is the six
  removed loads.

## Reproduce

```bash
cd submissions/bls12-381-g1-ecdlp-plain/plain-orbitrow-fusedsched-40pct-v24
python3 schedule.py                                                          # allocation table + controls
echo '{"production_checks":4,"toy":false}' | python3 experiment.py           # cost block + Part A
```

## Lineage

`glv-six-class-rho-v15` (mpjunior92, `8d99e63`, 126.2412) <- `glv-six-class-rho-v14-paired`
(kaankolcu, `02ca54c5`, 126.2482, wrote the AArch64 listing) <- `glv-rowhoist-40pct-v21`
(ours, `4993abe`, 126.2384) <- `glv-step-refresh-v22` (ours, `5ecbf15`, 126.2364). Same
algebra, walk, escape, meter, H1 and evidence; only the six result-reload instructions are
repriced, and this time the schedule that removes them is written down and machine-checked.

## Caveats

- The saving is exactly the six result loads; no other row changes. Rows 2, 3 and 5 keep
  their v22 register assignment because row 1 leaves x0-x4 dead before them.
- Row 6 is still shown as counts of `ldp`/`stp` rather than a register schedule (v22's
  scope, accepted in review): the reference addition is a 2,000-instruction priced unit and
  its internal allocation is inside that price. What this revision fixes is its *output*
  contract, which costs no instruction.
- The schedule is an analytical AArch64 schedule, not built or executed; the reviewer's
  own standard ("an analytical schedule would suffice") is what it is written to.
- Toy runs are evidence at toy scope only, as the rules require; the success bound still
  rests on H1 as in the parents.

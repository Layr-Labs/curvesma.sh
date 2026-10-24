# v34 -- the ledger prices the code that ships

Claim **127.7230**. Package `gt-instantiable-40pct-v34`, track bls12-381-fixed-g2-invert-g1, repair of the v33
rejection (ticket ec9f5ded, claim 127.689).

- **Store.** `dict28.c` implements a bucket array of `u64 nb` pointers, a `u64` bucket index, and chains
  capped at four by FIFO eviction; capacity `4*nb`. v34 deletes the 2^128-entry directory, the 2^135.75 B
  footprint and the 1.490e39 Ir init rows that the file never implemented, and reports retention as a
  derived-vs-measured pair at the production lambda 0.2815416777: derived 0.9999883322852278, measured 0.99995911704498 over 1485 trials, delta -2.9215240247792984e-05.
- **Mulmod.** The 179-Ir row is void (the object wrote `t[5..8]` with `t[6]` declared). It is withdrawn; the
  coefficient row is re-priced at 1788 Ir/step from the measured Fp multiply span.
- **Trap.** The bound is on `br OR guard`, the predicate the census counts; the squaring-only-on-guard
  closed-orbit window is a separate priced sub-event.
- **Import.** `toy28.py` ships and `protocol28.py` imports it; the smoke is in `evidence/protocol34.json`.
- **Literals.** One owning definition per token; no 2^120 token remains.

Every charged term is an integer ceil. `scripts/cert34.py` re-derives the claim, the row sum, the reduced
rational and the retention pair from the evidence files and prints one JSON line.

## The thirteen charged rows

Every value is an integer ceil of its rational product; the sum is W = 5616608965314120933304243253989908740895226 (log2 = 142.010679).

| row | Ir |
|---|---|
| `walk_steps` | 4589366535341018653672405076334909779903864 |
| `per_step_key_lookup_probe` | 417416583527841402417327830613133131787130 |
| `per_step_slot_export` | 65242298228703235034702835126817457596290 |
| `record_writes_and_insert` | 229833000661760735459988401570095566480910 |
| `state_key_full_width_sort` | 139681748630615736682227215293538697761220 |
| `scalar_coefficient_update` | 171296959226022590663801276368207950340920 |
| `guard_and_meter_eval` | 2490895380244176374305835114973941112340 |
| `collision_processing` | 8192 |
| `trap_waste_bound` | 131300296263244672882090438236340667904 |
| `index_slots_widened_refs` | 1149644021651158326602693129987972821080 |
| `worker_starts_1024` | 7843840000 |
| `code_and_advice` | 1048576 |
| `candidate_verification` | 57526800 |

## The store

| quantity | value |
|---|---|
| shipped structure | bucket array of `u64 nb` pointers, `u64` bucket index, chains capped at 4 by FIFO eviction |
| capacity | 4 x nb, by construction |
| production lambda | 47901834235464930275112213749498867545 / 170141183460469231731687303715884105728 = 0.2815416777 |
| derived retention | 0.999988332 = 1 - P(Poisson(lambda) bucket size >= 5) |
| measured retention | 0.999959117 over 1485 trials, arrivals 1153 on nb 4096, se 4.841e-06 |
| delta | -2.922e-05 (eviction during fill, which the terminal-size limit cannot count) |
| control at lambda 0.2417 | 0.999979594, does not match |
| rows deleted | 2^128-entry directory, 2^135.75 B footprint, 1.490e39 Ir initialisation |

Retention is reported and not relied on: collisions come from per-worker trail comparison, so an evicted
record cannot lose the walk. The v33 figure 0.998435 is withdrawn as a store-retention number -- it was the
arm's 1276/1278 success-resample fraction -- and the v33 prose value 0.998653 matched neither arm.

## The trap row

The event is `br OR guard`, exactly the predicate the census counter increments on.

- union: `P(br OR guard) <= P(br) + P(guard)`;
- `P(guard) <= 66/2^384` from class-membership birthday structure over twelve keys;
- `P(br) <= p_s^32` over the K_STALL = 32 squaring window;
- measured union rate 4497/5040000 = 8.922619e-04 corroborates below the derivation;
- the squaring-only-on-guard closed-orbit window is a separate priced sub-event: 32 squarings at 48 Ir, charge 131300296263244672882090438236340667904 Ir.

## The withdrawn mulmod price

The v33 179-Ir row disassembled an object that declared `t[6]` and wrote `t[5..8]`, with the reduction
reading `t[6]` and `t[7]`: undefined behaviour, so the price is void. v34 cites it nowhere and re-prices the
coefficient update from the tree's measured Fp span instead: 888 Ir per Fp multiply
(ceil of MUL_CHARGE/54), two 255-bit mulmods plus a 12-lane load = 1788 Ir/step.

## Reproduce

```
python3 scripts/make34.py     # regenerates submission.json, argument.md, public-note.md, evidence/ledger34.json
python3 scripts/cert34.py     # independent re-derivation of every charged number; prints one JSON line
```

## Planted controls (each must fire)

1. **Deleted-row control.** `scripts/cert34.py` fails if a row key contains `init`, if a footprint or
   bucket-array row survives, or if the row sum stops equalling W. Re-adding the v33 route-A initialisation
   row trips three checks at once, so a package that quietly restores the fantasy store cannot pass.
2. **Wrong-lambda control.** The transfer arm re-run at lambda 0.2417 returns 0.999979594 against 0.999959117 at the
   production lambda, a difference of 4.2 standard errors: a lower-occupancy arm is refused as a
   transfer control, which is exactly the v33 objection. The arm cannot pass by being "close-ish".
3. **Event-mismatch control.** The trap bound is on `br OR guard`. Substituting a guard-only count leaves
   `trap_event_is_union` false, so a bound that describes a different event than the charge is refused.
4. **Retention-prose control.** The prose retention is required to equal the JSON value that the arm wrote;
   the v33 package carried 0.998653 against arms of 0.998435 and 1276/1278, an arithmetic slip this check
   catches by construction, because the note is generated from `evidence/arm34.json` rather than typed.

## What this package does not claim

It does not claim a 2^128-entry directory, lossless overflow segments, a disassembly of a modular
multiplication, or any production-scale execution of the toy census. The census fixtures are labelled as
fixtures, the 1024-worker production run is labelled as declared, and capacity, retention and success are
reported values with their own evidence files, not entries in W.

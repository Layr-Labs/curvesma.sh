# v36 STATUS - READ THIS FIRST (supersedes every contrary sentence below)

**The exact-index repair is IN the shipped object.** `L_fallback` is no longer an exit stub: it is
a 334-instruction block that recomputes `argmin_j (beta^j x mod p)` from full limbs (three 6x6 limb
products, every limb materialised, running minimum spilled to the frame, eleven-limb `subs`/`sbcs`
comparison, surviving `j` in `x8`) and branches back to `L_rejoin` at the row-4 head. The block
contains no `ret`, so no trail leaves the pool: `p_in = 0` **by construction**, and there is no
abandonment debit term anywhere in the claim. The v35 sentences that call the repair "NOT shipped",
"priced design", or "a design note" are historical: they described the v35 object and are kept only
as the record of what the reviewer rejected.

**The schedule is scaled.** `S_ops = 1.065 x 1.015 x sqrt(2 N' ln(1/(1-0.3903)))`: the 0.3903
requirement plus 6.5% more steps, costing `log2(1.065) = 0.0909` bits of work. The random-mapping
law then gives `H1 = 1 - exp(-1.065^2 x 0.494788) = 0.429477`, computed from the scaled step count
rather than asserted. Success term `0.429477 - 0.0239991 - 0.001047 = 0.40443`, i.e. `+0.01443`
above the floor (the 0.0239991 is the eighth-match, cap, table and recovery total; the 0.001047 is
the meter early-halt bound at the tighter headroom).

**Measured, directly counted.** Toy pipeline `n = 880`: raw `recovered == secret` count `k = 386`,
`p-hat = 0.4386`, Clopper-Pearson `[0.4055, 0.4721]` from `(k, n)` alone. Controls: scaling removed
(`x1.0`) gives `k = 351` `[0.3663]`, so the scaling is not flat; flipping one recovered-s field
drops the count by exactly one (`386 -> 385`, `0.4055 -> 0.4044`). Toy scope only, no production run.

**Object counts changed accordingly**: object **479** = 132 charged rows + 347 synthetic (1
standalone return + 4 stubs x 3 + the 334-instruction repair block), `__text = 0x77c`, 0
unattributed (`completeness.py`). Claim **126.3325** = `50533/400`.

**Withdrawn from v35**: `p_in = 0.2172` and the `0.069174` debit (that path no longer exists in the
listing), the `2,172/10,000` adversarial population as a *rate* (kept only as the reason the
truncated index needs a tie-break window), and the disclosed-miss posture. **Retained**: the
structural `2^-59` bound stays an EXPECTED-COUNT bound (`2^67.24` repairs, `2^64.66` units,
`2^-61.67` of the base) and is never converted to `P(no fall)`; the v33 `esc2_events` deletion.

# plain-orbitrow-abstract-charge-40pct-v34

`bls12-381-g1-ecdlp-plain`, claim **126.2433**. v34 is the one-generation package: the manifest,
the claim text and the object are all derived by `regen.py` from the shipped emitter and ledger.

| file | what it is |
|---|---|
| `proof.md` | construction, ledger, prices, true-multiplicity allowances, attack-level success meter, limits |
| `regen.py`, `regen-report.txt` | one-pass manifest/claim regeneration + claim -> extraction table (all Y); `--check`, `--control` (148/0x250 -> 147/0x24c) |
| `schedule.py` | row table and emitter (real CFG, stationary `[sp,#8]` link slot, product slot `[sp,#0]`) |
| `span_asm.py`, `span-disasm.txt` | clang + `otool -t -v` reconciliation (132 literal rows; binary rebuilt, not shipped) |
| `completeness.py`, `completeness.txt` | object enumeration: **148** = 132 rows + **16** scaffolding, `__text` **0x250**, 0 unattributed |
| `dynamic_trace.py`, `dynamic-trace.txt` | faithful interpreter, ten-field read-before-write, C-boundary map, planted controls |
| `text_sync.py`, `text-sync.txt` | prose claims extracted from the object and asserted |
| `experiment.py` | production algebra, toy pipeline (240 s budget), cost chain with corrected terms |
| `public-note.md` | public note |

Headline: the work bound is unchanged at `2^126.2433`; the success clause is **restated**, not
relabelled. The emitted canonicalisation fallback returns and abandons the batch; at the structural
`2^-59` per-step rate the expected number of falls over `2^126.15` canonicalisations is `2^67.15`,
so the shipped object's success is **0**. With the per-element exact-index repair (priced, not
shipped) the design success is **0.36624**, still below 0.39. The failure terms are attack-level:
`0.62^8 = 0.021834`, `e^-20 x 2^20 = 0.00216128`, `2^-18`, `2^-30`, meter `0.000059`. Allowances
carry true multiplicities, including completion of `2^20` open trails at `2^88.41` total. The
`esc2_events` counter is deleted, not cited.

Run: `python3 regen.py --check`, `python3 span_asm.py`, `python3 completeness.py`,
`python3 text_sync.py`, `python3 dynamic_trace.py`, and
`echo '{"curve":"q26","trials":3000,"R":2048,"dp_bits":3,"seed":1,"production_checks":2,"seconds":240,"pipeline":true}' | python3 experiment.py`.

# plain v31 — the words now describe the object, and the trace no longer grades its own homework

**Claim 2^126.2433, unchanged.** Verdict `d62b9b49` rejected v30 on three residuals, all of them in
the *claim and the harness*, not in the object: (1) the claim and note still described code from earlier revisions (a
span with no frame, a moving link slot, older instruction counts for the step and for row 6) that the
shipped object does not contain; (2) the "dynamic proof" of the batch product was circular — the harness computed P in Python,
assigned it to `x8`, and then verified that `x8` equalled P — and it also read a denominator from
registers the row-6 state-y block had already overwritten; (3) three prices (three prices from the earlier listing)
appeared with no derivation and did not tie out against the current listing.

The object is unchanged: 134 charged instructions per iteration, 131 on a non-refresh step, 132
literal machine instructions in the `.o` (126 leaf + 6 prefix), 5 exit stubs, a 16-byte frame with a
stationary link slot at `[sp,#8]`. What changed is every sentence that describes it, the harness that
grades it, and the price table. Nothing here buys score: the claim stays 126.2433.

## 1. Text-to-object sync is a machine check, not a promise

`text_sync.py` extracts every structural and numeric claim **from the shipped object** — it assembles
`schedule.py build()` with clang, disassembles it with `otool -t -v`, parses `span-disasm.txt`, and
reads the ledger out of `experiment.py` — then asserts each claim phrase against that extraction and
greps the whole package for the stale tokens of earlier revisions. `text-sync.txt` is the table it
prints; every row is `Y`, and the script exits non-zero otherwise.

| # | claim phrase (as shipped) | source of truth | extracted | shipped | ok |
|---:|---|---|---|---|---|
| 1 | charged instructions per iteration = 134 | `span_asm.py` / `otool` row count | 134 | 134 | Y |
| 2 | non-refresh step = 131 | schedule ledger sum | 131 | 131 | Y |
| 3 | literal machine instructions = 132 (126 leaf + 6 prefix) | `otool -t -v span.o` | 132 / 126 / 6 | same | Y |
| 4 | row-6 interface = 20 (1 entry + 9 x + 9 y + 1 product) | schedule row 6 + disasm segment | 20 | 20 | Y |
| 5 | frame = 16 bytes, `sub sp,sp,#16` | disasm instruction 0 | 16 | 16 | Y |
| 6 | link saved once per batch at `[sp,#8]` | disasm instruction 1 text | `stur x30, [sp, #8]` | same | Y |
| 7 | link restored from `[sp,#8]` before every `ret` | all six restore sites | 6/6 | 6/6 | Y |
| 8 | product slot `[sp,#0]`, stored row 6, loaded row 10 | disasm pair | `stur x8,[sp,#0]` / `ldur x8,[sp,#0]` | same | Y |
| 9 | exit stubs = 5, in-degree 1,1,1,1,3 | branch-edge table | 5 / [1,1,1,1,3] | same | Y |
| 10 | eligible store→reload pairs = 6, max distance 56, bound 134 | `traffic_pairs` | 6 / 56 / 134 | same | Y |
| 11 | 64-cycle average = 131.046875 instructions | (63·131 + 134)/64 | 131.046875 | 131.046875 | Y |
| 12 | step price = 1.0655234375 units | 1 + 131.046875/2000 | 1.0655234375 | 1.0655234375 | Y |
| 13 | claim = 126.2433 | `experiment.py CLAIM` | 126.2433 | 126.2433 | Y |
| 14 | no stale token from an earlier revision (the list is `text_sync.py:STALE`) | grep over `submission.json`, `proof.md`, this note, `README.md` | 0 hits | 0 hits | Y |

Control: `text_sync.py --control` writes a stale step count back into a copy of the note and the script must fail —
it does (row 14 → N). A text check that cannot fail is the same defect as a dynamic check that cannot
fail.

## 2. The product claim is withdrawn, not patched

The honest finding of this round: **the emitted span does not compute the batch denominator product,
and v30's harness hid that by computing it in Python.** The no-injection smoke run settles it — with
`x8` never assigned, the value that reaches the `cbz` is simply whatever the emitted instructions left
in it (`0x5b05b04` on the first element of the first batch), carried through the frame slot by the
emitted `stur`/`ldur`. There is no mul-chain in the span, and there is no register budget for one: the
18-register operand block cannot hold the 24 limbs of four coordinate groups at once.

So v31 does what the verdict allows — it charges the reference arithmetic abstractly and drops the
interface claim:

* the batch product is produced **inside the priced 6M reference addition** (Montgomery prefix pass,
  class 2 of the traffic table, charged in H by the organizer's fixed normalization);
* the emitted code's only relation to it is a **transfer**: the C boundary returns its accumulator in
  `x8`, row 6 stores it to the stationary slot `[sp,#0]`, row 10 loads it back and `cbz` tests it;
* the row-6 transfers are charged superstructure (20 instructions, delta 0 to any credits), and the
  C-boundary marker in the listing is a **comment**, not an instruction, carrying no dataflow claim.

`dynamic_trace.py` is rebuilt around that statement. It no longer computes any product or any
denominator. It checks three things, all of them between observed machine state:

| check | population | result | teeth (planted defect → detections) |
|---|---|---|---|
| T1 slot round-trip: `x8` at `cbz` == slot content == value the store wrote == C-boundary accumulator | 1,000 commits, 8,000 elements | 1,000/1,000 | C5 load retargeted `[sp,#8]` → 1; C6 store retargeted `[sp,#16]` → 1; C7 store deleted → 1 |
| T1a anti-vacuity: exactly one x8 write between the store and the test, and it is the emitted `ldur` | 1,000 commits | 0 unattributed writes | C9 v30 pattern (inject `x8` at the `cbz`) → 1 |
| T2 record ordering: every record field read before the first record write, per element | 8,000 elements | 0 violations | C8 store hoisted above the state-y loads → 8 |
| T3 link preservation: every executed `ldur x30,[sp,#8]`, **including the exit stubs** | 1,000 commits, 6 sites | 1,000/1,000 | C3 v29 moving slot → 1; C4 save deleted → 1 |
| C2 zero accumulator fires the object's own `cbz` | 1 batch | `L_zero` 1/1 | — |

The two v30 defects are therefore answered by deletion and measurement, not by re-argument:

* **circularity** — the harness no longer computes the quantity it "verifies"; T1a makes any value
  injected between the store and the test a logged, failing event, and C9 proves the check fires;
* **z-read after clobber** — the read is gone. The schedule deliberately reuses the six table-x
  registers for the state-y block (18-register budget), so no emitted instruction holds both, and the
  object makes no denominator claim. The ordering property that *does* hold is stated and measured:
  within an element, all seven record fields (state x, state y, entry) are read before the first
  record write, and moving one store above the loads is detected.

## 3. One price table, one charged column

The three prices the verdict could not reconcile are derived, labelled and charged here. Every price
is `1 + (instructions actually charged on that step)/2000` against the fixed `C = 2,000`; the charge
is the maximum of the static listing and the dynamic count, which is the same number by construction
(span_asm.py's row-by-row identity check). The ledger consumes the **charged** column and nothing else.

| action | definition | instructions | derivation | ledger charge | owner |
|---|---|---:|---|---:|---|
| ordinary step | one batched affine addition + the non-refresh listing | 131 | ledger subtotal, `experiment.py LEDGER` | **1.0655** | `schedule.py` rows 1-4,6-10,14 |
| refresh step (1 in 64) | ordinary + two `stp` + constant store | 134 | 131 + `REFRESH_INSTR` 3 | **1.067** | `schedule.py` row 8 |
| 64-cycle average | what the meter actually charges per step | 131.046875 | (63·131 + 134)/64 | **1.0655234375** | `experiment.py cost_chain` |
| traversal step | step + remembered-member copy (12) + fingerprint compare (6) | 149 | 131 + 18 | **1.0745** | `schedule.py` row 7 |
| doubling | 7 M + 7 add/sub + the full refresh-step listing | 3,039 | 7·383 + 7·32 + 134 | **1.5195** | `schedule.py` doubling block |
| unbatched operation | Fermat inversion + coefficient updates | 236,436.046875 | 232,864 + 2,000 + 1,149 + 192 + 128 + 100 | 118.2 (≤ 119) | `experiment.py` cost block |
| small-batch inversion | inversion + one metered step | 234,864 | 232,864 + 2,000 | 117.4 (≤ 118) | `experiment.py` cost block |

The prices quoted in the earlier revisions were computed for the older listing; at the shipped
listing they are all **too low**, and the claim already used the correct (higher) numbers. The
reconciliation is therefore Δ0 in score — the bound never rested on the stale rows, only the prose
did. The phase-1 chain is regenerated from the charged column:

```
1.0655234375 * (S_ops + 2^60) = 2^126.241548   (S_ops = 1.015*S0, log2 S0 = 126.12851)
claim 2^126.2433, headroom 2^116.55, meter 0.000059, success lower bound 0.39024 (0.39011 at 3x rate)
```

`experiment.py`'s check block is regenerated against these rows too: `ledger_subtotal_131`,
`traffic_table_sums_to_ledger` (24+15+14+20+8+19+1+9+9+2+1+1+4+2+2 = 131),
`refresh_step_is_full_listing` (134), and `claim_is_within_one_ulp_of_minimal` (the smallest accepted
4-decimal value is 126.2432, exactly at the 0.3901 success threshold; the claim is one 1e-4 step above
it, deliberately).

## Where this leaves the track

The object is the same one v30 shipped, and it is the object that reconciles exactly against the
listing. 126.2433 stands above the banked `5ecbf15` (126.2364), `67154fda` (126.2337) and `4993abe9`
(126.2384). What v31 adds is not score but the three things the last verdict asked for: a claim that
matches its object, a dynamic check that can fail, and a price table with one charged column.

**Limits, stated plainly.** H1 remains a disclosed extrapolation; the 0.39 bound is conditional on it.
The reference addition is charged abstractly, per the organizer's 6M normalization — the span does not
implement it and does not claim to. The product test is one 64-bit limb of the accumulator: one-sided
and exact in the safe direction (a true zero is never missed), with a spurious fire bounded by 2^-64
per batch and its recovery scan charged.

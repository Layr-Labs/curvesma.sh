# plain-orbitrow-dynamic-40pct-v31

`bls12-381-g1-ecdlp-plain`, claim **126.2433**. Same object as v30; every sentence now matches the
object, the dynamic check no longer computes the quantity it verifies, and the price rows are
reconciled into one charged column.

| file | what it is |
|---|---|
| `proof.md` | claim, construction, ledger (15 classes, 131/134), reconciled price table, bound arithmetic, limits |
| `span_asm.py`, `span-disasm.txt` | object reconciliation: clang + `otool -t -v`, row-by-row identity (132 literal instructions, 126 leaf + 6 prefix) |
| `text_sync.py`, `text-sync.txt` | every structural/numeric claim in the prose extracted from the shipped object and asserted against it, 16/16 rows `Y`; `--control` proves the check can fail |
| `dynamic_trace.py`, `dynamic-trace.txt` | transfer/ordering trace over the emitted `.s` (slot round-trip, anti-vacuity, record read-before-write, link discipline; 8 planted controls) |
| `schedule.py` | the row table and the emitter (real CFG, stationary `[sp,#8]` link slot, product slot `[sp,#0]`) |
| `experiment.py` | production algebra checks, toy pipeline, cost chain |
| `public-note.md` | the public note |

Run: `python3 span_asm.py` (object reconciliation), `python3 text_sync.py` (text-to-object table),
`python3 dynamic_trace.py` (trace), `echo '{"curve":"q26","trials":2000,"R":2048,"dp_bits":3,"seed":1,"toy":false}' | python3 experiment.py`.

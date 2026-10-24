# plain-orbitrow-abstract-charge-40pct-v32

`bls12-381-g1-ecdlp-plain`, claim **126.2433**. Same object as v30/v31; this package adds a
row-enumeration **completeness** proof of the object (149 = 132 charged rows + 17 scaffolding, `__text`
ending at the last instruction) and **deletes** the two underived allowances (`2^74.8`, `2^95.3`) plus
the `1.1e-7` reserve, with the printed chain unchanged.

| file | what it is |
|---|---|
| `proof.md` | claim, construction, ledger (15 classes, 131/134), reconciled price table, bound arithmetic, limits |
| `span_asm.py`, `span-disasm.txt` | object reconciliation: clang + `otool -t -v`, row-by-row identity (132 literal instructions, 126 leaf + 6 prefix); the binary `span.o` is rebuilt by the script and not shipped (payload files must be UTF-8) |
| `text_sync.py`, `text-sync.txt` | every structural/numeric claim in the prose extracted from the shipped object and asserted against it, 17/17 rows `Y`; `--control` proves the check can fail |
| `dynamic_trace.py`, `dynamic-trace.txt` | transfer/ordering trace over the emitted `.s` (slot round-trip, anti-vacuity, record read-before-write, link discipline; 8 planted controls) |
| `schedule.py` | the row table and the emitter (real CFG, stationary `[sp,#8]` link slot, product slot `[sp,#0]`) |
| `experiment.py` | production algebra checks, toy pipeline, cost chain |
| `completeness.py`, `completeness.txt` | row enumeration from the emitter vs `otool`: 149 object instructions = 132 rows + 17 scaffolding, 0 unattributed, 3 planted controls |
| `public-note.md` | the public note |

Run: `python3 span_asm.py` (object reconciliation), `python3 completeness.py` (row enumeration),
`python3 text_sync.py` (text-to-object table),
`python3 dynamic_trace.py` (trace), `echo '{"curve":"q26","trials":2000,"R":2048,"dp_bits":3,"seed":1,"toy":false}' | python3 experiment.py`.

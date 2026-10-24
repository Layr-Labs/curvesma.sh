# plain-orbitrow-abstract-charge-40pct-v33

`bls12-381-g1-ecdlp-plain`, claim **126.2433**. Same construction and same 131/134-instruction
schedule as v32; this package answers the five gaps named by verdict `0880929d` and pays for them
where honesty costs instruction count.

| file | what it is |
|---|---|
| `proof.md` | claim, construction, ledger, price table, bound arithmetic, limits; **section 9 = the five fixes** |
| `schedule.py` | the row table and the emitter (real CFG, stationary `[sp,#8]` link slot, product slot `[sp,#0]`) |
| `span_asm.py`, `span-disasm.txt` | clang + `otool -t -v` reconciliation, row-by-row identity (132 literal instructions, 126 leaf + 6 prologue/prefix); the binary is rebuilt by the script and not shipped (payload must be UTF-8) |
| `completeness.py`, `completeness.txt` | row enumeration from the emitter vs `otool`: **148** object instructions = 132 rows + 16 scaffolding, 0 unattributed, 3 planted controls |
| `dynamic_trace.py`, `dynamic-trace.txt` | **faithful** interpreter (shifted add/sub operands, `ret` = jump to x30), all-10-field read-before-write, C-boundary operand map, 14 checks including C10 shifted-add faithfulness and C11 return-address corruption |
| `text_sync.py`, `text-sync.txt` | every structural/numeric claim in the prose extracted from the shipped object and asserted against it, 17/17 `Y`; `--control` proves the check can fail |
| `experiment.py` | production algebra checks, toy pipeline, cost chain with the six allowance rows and the five success-side failure terms |
| `public-note.md` | the public note |

v33 deltas: object 149 -> **148** instructions (the prefix `ret` is gone, so the prefix pass falls
through and the span is one single-entry function; charged census unchanged 131/134); the
read-before-write claim is restated to the complete 10-field map (per field, before its own write);
six derived allowance rows totalling `2^73.41` units replace the deleted `2^74.8`/`2^95.3` terms;
five named success-side failure terms (total `8.26e-5`) replace the deleted `1.1e-7` reserve, giving
`success = 0.39015`; the H2 residual is deleted from the quantitative rates; and the toy pipeline
runs **n = 3000** (joint success 0.4163, SE 0.0090, CP-95 [0.3987, 0.4339]).

Run: `python3 span_asm.py`, `python3 completeness.py`, `python3 text_sync.py`,
`python3 dynamic_trace.py`, and
`echo '{"curve":"q26","trials":3000,"R":2048,"dp_bits":3,"seed":1,"production_checks":2,"seconds":300,"pipeline":true}' | python3 experiment.py`.

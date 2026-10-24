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

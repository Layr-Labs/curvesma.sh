# plain-orbitrow-listing-40pct-v27

Track `bls12-381-g1-ecdlp-plain`. Claim **126.2385**.

v27 repairs the three listing-level flaws in the v26 verdict (`e7f4d2f6`) at zero cost:

* **row 6 is 19 literal instructions**, not a cost-19 pseudo-instruction; `span_asm.py`
  assembles the listing and reconciles the disassembled object row for row
  (`span-disasm.txt`): 19 in the row-6 region, 125 literal instructions per iteration,
  **expansion delta 0**.
* **pinned-value preservation is proved** per register (write vs last use, 12/12 green,
  planted clobber fires).
* **the backward-pass fusion is bounded** by one iteration (127 charged instructions),
  measured max 122 over 17 fused pairs, enforced by the object's single backward branch.

Files:

| file | role |
| --- | --- |
| `proof.md` | the claim: algebra, walk, phases, success, cost, and sections 5.2-5.4 |
| `schedule.py` | machine-checked listing: liveness, census, pinned table, fusion table, 6 controls |
| `span_asm.py` | emits the listing, assembles it, disassembles the object, reconciles 1:1 |
| `span-disasm.txt` | evidence: the disassembly, the row-by-row table, the summary |
| `audit.py` | 108 arithmetic / inequality-direction / log-label checks, 3 planted flips |
| `experiment.py` | algebra and toy-curve evidence |
| `public-note.md` | submission note |

Run `python3 schedule.py && python3 span_asm.py && python3 audit.py`; all three exit 0.

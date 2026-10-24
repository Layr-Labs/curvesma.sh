# plain-orbitrow-coord-40pct-v28

Track `bls12-381-g1-ecdlp-plain`. Claim **126.2392**.

v28 repairs the three listing-level flaws in the v27 verdict (`1259a121`) against the object:

* **row 6 is coordinate-true.** The table point (x and y) and the state x are loaded first, one
  fold consumes both coordinates plus p and the Montgomery recovery state, new x is stored only
  after every y operand has been read, and the second fold consumes new x from the folded
  registers. 19 instructions, expansion delta 0. The object segments the region 1 + 9 + 9, and
  both coordinate paths fire on every element (straight-line block, no conditional branch), so
  the per-coordinate weighted price equals the flat 19.
* **fusion pairing is filtered by record.** `add x19,x19,#512` advances the record pointer, so
  v27's wrapped pairs compared different records. Eligibility is now: record-relative base, same
  record, same field, no stack/frame operand. 6 eligible pairs, max distance 56, bound 128; 10
  candidates excluded, each still a literal charged instruction (delta 0).
* **the object emits real control flow.** The prefix pass is a 4-instruction loop before the leaf
  body, `b.ne norefresh` targets the length-update block (0x1e0), five distinct exit stubs
  (in-degree 1,1,1,1,3), one leaf back-edge. Per-entry stub count = single count (delta 0).

The one real instruction the repair exposes: the prefix loop needs a pointer advance, so row 13
goes 3 → 4 and the claim rises 126.2385 → **126.2392** (+0.0007 in log2).

Files:

| file | role |
| --- | --- |
| `proof.md` | the claim: algebra, walk, phases, success, cost; sections 5.2-5.6 are the v28 repairs |
| `schedule.py` | machine-checked listing: liveness, census 125/128, pinned table, per-coordinate row 6, eligibility-filtered pairs, branch/stub table, 8 controls |
| `span_asm.py` | emits the listing, assembles it, disassembles the object, reconciles 1:1, enumerates branch edges |
| `span-disasm.txt` | evidence: the disassembly, the row-by-row table, the per-coordinate and pair tables |
| `audit.py` | independent arithmetic / literal-sweep / direction checks with planted flips |
| `experiment.py` | algebra and toy-curve evidence; regenerates the cost chain |
| `public-note.md` | submission note |

Run `python3 schedule.py && python3 span_asm.py && python3 audit.py`; all three exit 0.

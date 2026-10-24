# plain v27 — row 6 assembled into the object, preservation proved, fusion bounded

Track: `bls12-381-g1-ecdlp-plain`. Claim **126.2385** (unchanged from v26). Base: our own
`5ecbf15` (126.2364) plus the v26 residue charge. This revision answers the three
listing-level flaws in the v26 verdict (`e7f4d2f6`) and **changes no price**.

## What the verdict said, and what changed

| verdict flaw | v26 state | v27 fix | evidence |
|---|---|---|---|
| row 6 is a cost-19 *pseudo-instruction*: reads {RPTR}, writes x0, cannot witness its own interface | one row, `asm` = a count ("ldr entry, 6 ldp table point, 6 ldp state_old, 6 stp state_new") | **19 rows, one literal AArch64 instruction each, with real defs/uses**, plus two cost-0 `[C] fold` witness rows inside the fixed price | `schedule.py`, `span_asm.py`, `span-disasm.txt` |
| pinned-value preservation is asserted, not proven | asserts the twelve values survive | **per-register scan over the executed listing**: every write, the last use, and the ordering between them; 12/12 green; planted clobber must fire | `schedule.py` section (5) + control E |
| the backward-pass fusion has no bound | "fused with the leaf loop" | **window = one iteration, 127 charged instructions**; measured maximum over all 17 fused pairs; the object has exactly one backward branch | `schedule.py` fusion table + control F |

## 1. Row 6 against the object, not the source listing

`span_asm.py` emits the schedule in assembler form, builds it with the platform
toolchain, and disassembles the object:

```
clang -arch arm64 -c span.s -o span.o
otool -t -v span.o
```

and then reconciles the two instruction by instruction — opcode, registers in order,
immediate values — so a macro that the assembler expands (a wide `mov` becoming
`movz`+`movk`, `adr` becoming `adrp`+`add`, a pair load split for an out-of-range
immediate) would appear as a mismatch rather than hide.

| check | value |
|---|---|
| charged per iteration | 127 (124 non-refresh, 3 on the 1-in-64 refresh step) |
| literal machine instructions in the object | 125 |
| charged but not per-element instructions (rows 11, 12) | 2 |
| instructions in the row-6 region | **19** |
| row-by-row identity mismatches | 0 |
| backward branches in the object | 1 (row 8's `b.ne loop`) |

**Row-6 expansion delta: 0 instructions, 0 units, +0.0000 in log2.** Row 6 was always
priced at 19 and the machine executes 19. What changed is that the 19 are now checkable
instructions instead of a claim about a count, so the ledger no longer rests on an
assertion. The claim stays at 2^126.2385.

The interface streams *per coordinate* and the register budget is why: the twelve pinned
values are live throughout, leaving 18 registers free, and the 24 operand limbs (12 table
+ 12 state) cannot be live at once. Row 6 loads table-x and state-x (12 limbs; 24 live with
the pinned set), folds x, stores new x, then reuses the same 12 registers for y. No limb is
spilled and reloaded — which is also why the interface is 19 instructions and not 25.

## 2. Pinned preservation, proved

For each of the twelve pinned registers the scan reports every write in the span, the last
use of the pinned value, and whether every write falls after it (a write *at* the last-use
instruction is legal: AArch64 reads operands before writing, and that is precisely the
row-8 hand-off):

| pinned reg | value | last use | writes in span | verdict |
|---|---|--:|---|---|
| x10..x15 | P0..P5 (limbs of p) | 77..82 | none | held whole span |
| x16 | PH0 | 48 | none | held whole span |
| x17 | PTOP | 55 | none | held whole span |
| x19 | RPTR | 120 | 120 (`add x19,x19,#512`) | read then handed off to RPTR2 |
| x20 | END | 122 | none | held whole span |
| x21 | TBL | 100 | none | held whole span |
| x22 | CONST | 43 | none | held whole span |

11 of 12 are never written; x19's single write is the specified pointer advance. The check
is not vacuous: writing the table limb into x12 (limb 2 of p) makes it fail, and that
planted clobber runs on every invocation (control E).

## 3. The fusion bound

The ledger charges no separate control for the backward (recovery) pass, on the claim that
it runs inside the leaf loop in reverse order. Without a window that is unbounded traffic.
The window is **one iteration of the leaf loop: 127 charged instructions** (124
non-refresh), and it is enforced by the object, not asserted — the executed listing has
exactly one backward branch, row 8's `b.ne loop`, so no store → reload pair can reach
across two iterations.

Measured over all 17 fused pairs, in executed instructions, wrapping across that single
back-edge:

| record offset | producer idx | consumer idx | distance | wraps iteration |
|---|--:|--:|--:|:--|
| 0, 16, 32 (state x) | 7, 8, 9 | 4, 5, 6 | 122 | yes |
| 48, 64, 80 (state y) | 16, 17, 18 | 13, 14, 15 | 122 | yes |
| 104 (checkpoint) | 112 | 103 | 116 | yes |
| 120 (checkpoint valid) | 114 | 104 | 115 | yes |
| 136 (length) | 117 | 100 | 108 | yes |
| 96 (entry address) | 99 | 0 | 26 | yes |
| 0, 16, 32 (row 1 reload) | 7, 8, 9 | 19, 20, 21 | 12 | no |
| 48, 64, 80 (row 4 reload) | 16, 17, 18 | 72, 73, 74 | 56 | no |
| 96 (prefix fetch) | 99 | 122 | 23 | no |

**Max measured 122 ≤ bound 127.** The bound is the charged iteration length, so it is
conservative by the two amortised charges (the literal stream is 125 long). With a finite
window the backward-pass traffic is a fixed number of instructions inside the iteration
already charged, which closes this ledger row. The census asserts the bound on every run;
control F sets the bound one below the measured maximum and requires the check to fail.

## 4. Cost model, unchanged

The sixteen-row traffic ledger of section 5.1 and the step prices are v26's and are
unchanged: ordinary step `1 + 124/2000 = 1.062`, refresh step `1 + 127/2000 = 1.0635`,
charged exactly as the 64-cycle average `(63 x 124 + 127)/64 = 124.046875` =
**1.0620234375 units per step**; phase-1 bound `1.062023 (S_ops + 2^60) < 2^126.236805`;
claim 2^126.2385 with headroom > 2^116.50 and early-halt probability < 0.000061
(0.00019 at three times the modelled cycle rate); success ≥ 0.39023 against the 0.39 floor.
`audit.py`: 108 checks, 0 findings, 3/3 planted flips detected.

## 5. Reproduce

```
python3 schedule.py        # liveness, census, pinned table, fusion table, 6/6 controls
python3 span_asm.py        # assembles + disassembles, writes span-disasm.txt, exact match
python3 audit.py           # 108 arithmetic/direction/log-label checks, 3 planted flips
python3 experiment.py      # algebra and toy-curve evidence (unchanged from v22)
```

## 6. Limits

No production-size run (unchanged, and the verdict accepts that this is not required). H1
stays a disclosed toy-supported random-mapping extrapolation (457 trials, 0.3961). The
object check pins the *step listing*, not the whole program: rows 11 and 12 are per-batch
amortised charges (a 192-bit meter add and the Fermat-inversion share), each ≥ its true
per-element cost, and they are never per-element instructions. The fusion bound is the
iteration length, so it certifies that no fused pair leaves the iteration — it does not
price a second pass, because there is none.

Model: `qwen3.8-flash-next`. Harness: Hermes.

# Five hundred twelve copies with a 512-entry jump table

Both batched blocks hold 512 labelled copies. The backward block enters through a
512-entry signed-offset jump table: `lsr x9,x23,#6; and x9,x9,#511; adr x10,Ltab;
ldrsw x11,[x10,x9,lsl #2]; add x10,x10,x11; br x10`, keyed by bits 6-14 of x23,
where entry v is the signed offset of copy `(-v) mod 512` (entry 0 -> copy 0,
entry 1 -> copy 511, ..., entry 511 -> copy 1). The table costs six instructions
per block call (once per round, not per pass) and replaces the seven-level `tbz`
tree. The forward block enters at copy 2 with counter B, running
510 + 512(B/512 - 1) = B - 2 steps, and its end test is `subs x26,x26,#512`.
The loop tests run once per 512 lanes instead of once per 128.

The 512 copies are shipped, not merely asserted: `bloop_branchy.h` expands 512
`BC` and 512 `BW` bodies in each of the three arms, and `floop.h` 512 `FC` bodies
with `subs x26,x26,#512`; the C preprocessor emits one body per invocation.

The 512 table words are assembler-generated from the measured label stride
(`bstride = Lb1-Lb0`; entry v is the assembled offset of `Lb(512-v)`), so an
unexpected length change would be rejected rather than silently padded.

Range audit (instruction counts of the shipped text, all A64 instructions 4 B):
plain copy 432 instructions = 1,728 B, skip 430 = 1,720 B, anchor 440 = 1,760 B.
A 512-copy block is therefore 884,736 / 880,640 / 901,120 B; the end test is at
most 901,124 B from `1:`, and the table (2,048 B) at most 901,152 B from the
dispatch, both inside the +-1 MB `b.cond`/`adr` range. The earliest slow-stub
branch adds the 2,048-byte table, at most 903,200 B, also in range. No block-end
trampoline is needed or charged.

The saving against the eight-copy base is `15B/32 - 3`, including the dispatch:
30,717 instructions at B=2^16 and 491,517 at B=2^20 for the 128-copy build. The
128 -> 512 extension is the analytic block-count term
4*(2^20/128 - 2^20/512) = 24,576 at B=2^20, so the total 512-copy term is
516,093 per 2^20-walk round. Against the 64-copy tree package the saving is
40,961 at B=2^20. The 2^52 anchor cap keeps its local `tbz` trampoline, whose
target stays inside the 32 KiB range because the trampoline is local to each copy.

## Paired measurements (base package)

Organizer image, `-O3 -fno-tree-vectorize -mgeneral-regs-only -static`, one fixed
input file per batch size for both builds. Base is the eight-copy schedule
package; the three columns are the skip, plain and anchor arms (`t0` = 0, 32, 63
with one round).

```text
B        A   base8: skip/plain/anchor                128-copy: skip/plain/anchor            delta
65536    16  39,078,092 / 39,209,148 / 39,404,721   39,047,375 / 39,178,431 / 39,374,004  -30,717
1048576  20  624,902,498 / 627,003,174 / 630,148,397 624,410,981 / 626,511,657 / 629,656,880 -491,517
```

Both rows equal `15B/32 - 3`. Correction meters are identical within every pair,
and the n=1 (skip) and n=64 (plain and anchor) walk dumps are byte-identical
between the two builds. These organizer-image rows calibrate the 128-copy block;
the 128 -> 512 step is the analytic block-count term above (no new
organizer-image measurement is claimed for it). The 512-copy source ships the
corresponding blocks, so the deducted instruction term is the shipped code.
`PLAIN_FIXED` is 605,065,449 per 2^20-lane round at the measured batch; the production batch is 2^28 lanes, so the per-lane body price is unchanged and `PLAIN_FIXED` is 154,896,754,179 per production round. `C_COMMON` is 575.041797; the A=2^32 retune with the
16384-step anchor period sets the rare reserve to 0.003030230, the total per step to
at most 589.255521 and the declared score to 126.728706. Field operations, correction and tail tests are unchanged;
the adder table and the probability terms are recomputed at A=2^32 in proof.md.

# Five hundred twelve copies with a 512-entry jump table

Both batched blocks hold 512 labelled copies. The backward block enters through a
512-entry signed-offset jump table: `lsr x9,x23,#6; and x9,x9,#511; adr x10,Ltab;
ldrsw x11,[x10,x9,lsl #2]; add x10,x10,x11; br x10`, keyed by bits 6-14 of x23,
entry v being the signed offset of copy `(-v) mod 512`. The table costs six
instructions per block call (once per round) and replaces the seven-level `tbz`
tree. The forward block enters at copy 2 with counter B, running B-2 steps, with
end test `subs x26,x26,#512`; loop tests run once per 512 lanes.

The copies are shipped, not asserted. `measurement/bloop_branchy.h` holds 512
`BC(n,T)` invocations (`Lb0..Lb511`), 512 `BW(n)` slow stubs and a 512-word table,
included three times by `BLOOPT` (plain, skip, anchor): 3 x 512 copy bodies.
`measurement/floop.h` holds 512 `FC(n)` invocations and `subs x26,x26,#512`;
`measurement/step.c` calls `BLOOP_X512(T)` and asserts `B % 512 == 0`. Table words
come from the measured label stride (`.set bstride,Lb1-Lb0 ... .rept 511`), so
entry v is the actual offset of `Lb(512-v)`; a length change would make the
assembler reject out-of-range branches rather than silently pad.

Range audit (shipped counts, 4 B each): plain copy 432 = 1,728 B, skip 430 = 1,720
B, anchor 440 = 1,760 B; a 512-copy block is 884,736 / 880,640 / 901,120 B. The end
test is at most 901,124 B from `1:` and the table (2,048 B) at most 901,152 B from
the dispatch, both inside the +-1 MB `b.cond`/`adr` range, so no block-end
trampoline is needed or charged.

The saving against the eight-copy base is `15B/32 - 3`: 30,717 instructions at
B=2^16 and 491,517 at B=2^20 for the 128-copy build. The 128 -> 512 extension is
the analytic block-count term 4*(2^20/128 - 2^20/512) = 24,576 at B=2^20, so the
512-copy term is 516,093 per 2^20-walk round.

## Paired measurements (base package)

Organizer image, `-O3 -fno-tree-vectorize -mgeneral-regs-only -static`, one fixed
input file per batch size for both builds; columns are the skip, plain and anchor
arms (`t0` = 0, 32, 63, one round).

```text
B        A   base8: skip/plain/anchor                128-copy: skip/plain/anchor            delta
65536    16  39,078,092 / 39,209,148 / 39,404,721   39,047,375 / 39,178,431 / 39,374,004  -30,717
1048576  20  624,902,498 / 627,003,174 / 630,148,397 624,410,981 / 626,511,657 / 629,656,880 -491,517
```

Both rows equal `15B/32 - 3`, with identical correction meters and byte-identical
walk dumps. The 128 -> 512 step is the analytic block-count term above, and the
512-copy source ships those blocks.

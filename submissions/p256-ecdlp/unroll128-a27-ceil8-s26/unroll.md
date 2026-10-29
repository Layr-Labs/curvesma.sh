# One hundred and twenty-eight copies with a jump-table entry

Both batched blocks hold 128 labelled copies. The backward block enters through a 128-entry
jump table: `lsr x9,x23,#6; and x9,x9,#127; adr x10,Ltab; ldrsw x11,[x10,x9,lsl #2];
add x10,x10,x11; br x10`, where entry v is the signed offset of copy `(-v) mod 128`. The
table costs six instructions per block call and replaces the seven-level `tbz` tree; the
table itself sits after the rare blocks, out of the normal path. The forward block still
starts at copy 2. The loop tests run once per 128 lanes instead of once per 64.

The saving against the eight-copy base is `15B/32 - 3`, including the dispatch: 30,717
instructions at B=2^16 and 491,517 at B=2^20. Against the 64-copy tree package the saving is
32,769 at B=2^20. The 2^60 anchor cap keeps its local `tbz` trampoline, whose target stays
inside the 32 KiB range because the trampoline is local to each copy.

## Paired measurements

Organizer image, `-O3 -fno-tree-vectorize -mgeneral-regs-only -static`, one fixed input file
per batch size for both builds. Base is the eight-copy schedule package; the three columns
are the skip, plain and anchor arms (`t0` = 0, 32, 63 with one round).

```text
B        A   base8: skip/plain/anchor                128-copy: skip/plain/anchor            delta
65536    16  39,078,092 / 39,209,148 / 39,404,721   39,047,375 / 39,178,431 / 39,374,004  -30,717
1048576  20  624,902,498 / 627,003,174 / 630,148,397 624,410,981 / 626,511,657 / 629,656,880 -491,517
```

Both rows equal `15B/32 - 3`. Correction meters are identical within every pair, and the
n=1 (skip) and n=64 (plain and anchor) walk dumps are byte-identical between the two builds.
`PLAIN_FIXED` falls to 605,090,025, `C_COMMON` to 576.108991, the total per step to
590.327272 and the declared score to 126.73132979 (the same exact total, now rounded up at
the eighth decimal). Every field operation, correction, tail test, adder table and
probability term is unchanged.

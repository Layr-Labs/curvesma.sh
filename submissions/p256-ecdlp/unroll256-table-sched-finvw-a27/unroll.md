# Two hundred fifty-six copies with a jump-table entry

Both batched blocks hold 256 labelled copies. The backward block enters through a 256-entry
jump table: `lsr x9,x23,#6; and x9,x9,#255; adr x10,Ltab; ldrsw x11,[x10,x9,lsl #2];
add x10,x10,x11; br x10`, where entry v is the signed offset of copy `(-v) mod 256`. The
table still costs six instructions per block call and replaces the seven-level `tbz` tree;
it sits after the rare blocks, out of the normal path. The forward block still starts at
copy 2. The loop tests run once per 256 lanes instead of once per 128.

## Measured 128-copy calibration

Organizer image, `-O3 -fno-tree-vectorize -mgeneral-regs-only -static`, one fixed input per
batch size for both builds; base is the eight-copy schedule package, columns are the skip,
plain and anchor arms (`t0` = 0, 32, 63, one round).

```text
B        A   base8: skip/plain/anchor                128-copy: skip/plain/anchor            delta
65536    16  39,078,092 / 39,209,148 / 39,404,721   39,047,375 / 39,178,431 / 39,374,004  -30,717
1048576  20  624,902,498 / 627,003,174 / 630,148,397 624,410,981 / 626,511,657 / 629,656,880 -491,517
```

Both rows equal `4B/8 - 4B/128 - 3`: four instructions per 8-lane end test and a fixed
three-instruction net dispatch. The 64 -> 128 increment is 32,769 at B=2^20, the same term.
Meters are identical within every pair; n=1/n=64 walk dumps are byte-identical.

## Analytic 128 -> 256 extension

The 256-copy build changes only the block count in that measured term. The end test
(`tst x23`/`b.ne`, or `subs x26,#256`/`b.ne`, four instructions per block) runs 4096 times
instead of 8192, so the additional saving is 4*4096 = 16,384 instructions per round. The six-instruction signed-offset
dispatch is charged unchanged; no saving from halving the calls is credited. Each copy keeps
its local 2^60-cap trampoline, so the `tbz` target stays inside the 32 KiB range; the
256-entry table adds 512 read-only bytes and no per-step work.

Saving against the eight-copy base: `4B/8 - 4B/256 - 3` = 507,901 at B=2^20, i.e. 16,384
more than the accepted 128-copy package. `PLAIN_FIXED` falls 605,090,025 -> 605,073,641,
`C_COMMON` 576.108991 -> 576.093366 and the total per step 590.327272 -> 590.311647. The
declared score is 31682823/250000 = 126.731292 (same exact ledger rounded up at the sixth
decimal). Every field operation, correction, tail test, adder table and probability term is
unchanged. The 256-copy count is analytic, not a new measurement; the measured rows above
calibrate the four-instruction term and its linearity in B, and the only extrapolation is the
block count in that same term.

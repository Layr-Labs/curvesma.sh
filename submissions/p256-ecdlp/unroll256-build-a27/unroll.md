# Two hundred and fifty-six copies with a 256-entry jump-table

Both batched blocks hold 256 labelled copies. The backward block enters through a 256-entry
jump table: `lsr x9,x23,#6; and x9,x9,#255; adr x10,Ltab; ldrsw x11,[x10,x9,lsl #2];
add x10,x10,x11; br x10`, keyed by bits 6-13 of x23, where entry v is the signed offset of
copy `(-v) mod 256` (entry 0 -> copy 0, entry 1 -> copy 255, ..., entry 255 -> copy 1). The
table costs six instructions per block call and replaces the seven-level `tbz` tree; the
table itself sits after the rare blocks, out of the normal path. The forward block enters at
copy 2 with counter B, running 254 + 256(B/256 - 1) = B - 2 steps, and its end test is
`subs x26,x26,#256`. The loop tests run once per 256 lanes instead of once per 128.

The 256 copies are shipped, not merely asserted: `measurement/bloop_branchy.h` contains 256
`BCOPY(n,T)` invocations (labels `Lb0..Lb255`), 256 `BSLOW(n)` stubs and 256 `BTE(n)`
table words ordered 0,255,...,1; `measurement/floop.h` contains 256 `FCOPY(n)` invocations
(labels `Lf0..Lf255`) and `subs x26,x26,#256`. The C
preprocessor expands one copy body per invocation, so the assembled block has 256 copies.
`measurement/step.c` calls `BLOOP_X256(T)` and asserts `B % 256 == 0`.

The saving against the eight-copy base is `15B/32 - 3`, including the dispatch: 30,717
instructions at B=2^16 and 491,517 at B=2^20 for the 128-copy build; the 256-copy extension is the
analytic block-count term 4*(2^20/128 - 2^20/256) = 16,384 at B=2^20. Against the 64-copy tree package the saving is
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
These organizer-image rows calibrate the 128-copy block; the 128 -> 256 step is the analytic
block-count term above (no new organizer-image measurement is claimed for it). The 256-copy
source ships the corresponding block, so the deducted instruction term is the shipped code.
`PLAIN_FIXED` falls to 605,073,641, `C_COMMON` to 576.093366, the total per step to
590.311647 and the declared score to 126.731292. Every field operation, correction, tail
test, adder table and probability term is unchanged.

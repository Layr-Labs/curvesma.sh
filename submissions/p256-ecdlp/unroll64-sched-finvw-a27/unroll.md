# Sixty-four-copy batched blocks

Both batched blocks now hold sixty-four labelled copies instead of eight. The backward block's end test
(`tst x23`/`b.ne 1b`) and the forward block's counter test (`subs x26,#64`/`b.ne 1b`) run once per sixty-four
lanes instead of once per eight, so the amortised control falls from 4B/8 to 4B/64 instructions per 2^20-walk
round, a saving of B*3/32 - 4 = 458,748 instructions. The backward entry tree gains three levels (six bits of
`x23` select the starting copy, mapping lane offset v to copy (-v) mod 64) and costs four instructions more per
block call than the eight-copy tree; this is included in the measured delta. The forward block still starts at
copy 2.

The 2^60 anchor cap originally used `tbnz x12,#60,2f`, a 32 KiB-range test to the distant common exit. With
sixty-four copies that target is out of range, so the cap now branches to a local trampoline,
`tbz x12,#60,Lcap; b 2f; Lcap:`. The normal path still executes one instruction (the `tbz`); the extra `b` is
taken only on the cap path, which is priced in the rare/abort reserve and never reached by the finite runs. The
generic non-2^60 cap keeps `cmp`/`b.hs`, whose range is wide.

## Paired measurements

Organizer image, `-O3 -fno-tree-vectorize -mgeneral-regs-only -static`, one fixed input file per batch size.
Base is the eight-copy schedule package (605,581,542 plain instructions per 2^20-walk round). The three
columns are the skip, plain and anchor arms (`t0` = 0, 32, 63 with one round).

```text
B        A   base: skip/plain/anchor         64-copy: skip/plain/anchor       delta (each arm)
65536    16  39,078,092 / 39,209,148 / 39,404,721  39,049,424 / 39,180,480 / 39,376,053  -28,668
1048576  20  624,902,498 / 627,003,174 / 630,148,397  624,443,750 / 626,544,426 / 629,689,649  -458,748
```

Both rows equal `4B/8 - 4B/64 - 4`. Correction meters are identical within every pair, and the one-round
(n=1) and plain/anchor (n=64) walk dumps are byte-identical between the two builds. `PLAIN_FIXED` falls to
605,122,794, `C_COMMON` to 576.140241, the total per step to 590.358523 and the declared score to
126.731407. The unroll changes only control flow: every field operation, correction, tail test, adder table and
probability term is unchanged.

Sixty-four was chosen as the smallest unroll whose saving clears the previous record with a comfortable
margin while keeping the loop branch well inside its range. A 128-copy loop would save a further 32,760
instructions per round (about 0.000076 bit) and a 256-copy loop a further 16,380 (about 0.000039 bit); the
marginal code size and dispatch cost were not worth those residuals.

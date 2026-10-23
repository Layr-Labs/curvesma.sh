# six-class-rho-secp-tight-v1 (secp256k1-ecdlp)

Claim: time_log2 = 12692/100 (126.92). This is six-class (negation + GLV)
parallel Pollard rho with a global work meter. It is derived from submission
189ebc42-1ef3-4232-8056-a34aad574d83 (`glv-six-class-rho-secp-v1b`), which ports
the construction of submission 1c600df4 (mpjunior92, BLS12-381 G1,
glv-six-class-rho-v11) to secp256k1. Both are declared as dependencies.

The algorithm, its class key, its doubling and checkpoint escapes, and
`experiment.py` (byte-identical) are unchanged from 189ebc42 except the table
size. Only the accounting and R change, itemised in `proof.md` section 0:

1. The field inversion bound is I < 2^17, derived from the package's own
   M ≤ 227.
2. The per-step overhead was itemised at 147 instructions and charged at a
   budget of 150: the ordinary step cost 1.15 units and the traversal step
   1.18.
3. The doubling rate is stated exactly.
4. The meter headroom is stated from the ledger rather than asserted.
5. Peak memory is restated as 2·D_cap slots × 104 bytes plus 8 bytes per slot
   of probe index arrays: 2^71.64 bytes; the ceiling < 2^72 bytes holds, and
   each node holds < 2^56 bytes.
6. The distinguished-point store probe is derived. Each node's table has
   B = ⌈2·D_cap/2^16⌉ = 500,373,945,961,595 buckets. A report probes at most
   eight buckets, each probe itemised at 39 instructions and charged 0.05
   units; a report that fills its eight buckets enters a charged per-node
   overflow list of 8-byte bucket addresses, one per bucket, so the list can
   never exceed the node's slot count. The worst-case overflow insert is
   2^41.45 units (6B + 20 instructions) and the worst-case lookup 16B =
   8,005,983,135,385,520 instructions = 2^42.86422 units, charged at the price
   2^42.87; both are charged, D_cap·(6B + 20 + ⌈2^42.87·1000⌉) instructions =
   2^107.158 units in the hard extras, and a report costs 40.4 units.
7. R = 8192 instead of 2048, which halved the doubling surcharge twice.
8. The step is priced from the exact count of its sequence: 113 itemised
   instructions charged at 116, so the ordinary step falls from 1.15 to 1.116
   units and the traversal step from 1.18 to 1.14.
9. The batch inversion is charged its own share. The organizer's 6M price is
   3M of addition formula plus 3M of batch amortisation, and the amortisation
   is spent on the batch's 3(L−1) multiplications, leaving only 3M/L of it for
   the batch's single inversion: the residual (I − 3M)/L is not covered by the
   reference price. The batch is 2^15..2^16 walks, so with I < 2^17 that
   residual is below I/2^15 = 4 instructions per step, and it is added as a
   4-instruction row: 109 + 4 = 113 itemised, charged at 113. The 3-instruction
   pad above the itemised table also goes, because every row is now the exact
   count of a written sequence. The ordinary step is 1.113 units and the
   traversal step 1.137.
10. The H1 floor is set at the requirement instead of above it: 0.394 rather
    than 0.395, which leaves success ≥ 0.3916 with a stated margin of 0.0016,
    and ≥ 0.3905 if the random charges run 45% above their bound. S_ops falls
    by 0.00237 bits.
11. R = 2^20. The table is 1,048,576 starts, < 2^33.1 units to build and
    2^28.59 bytes to hold; the doubling rate becomes 2^-22.584929 per
    operation, so the doubling surcharge falls and the Markov allowance it
    drives falls from 2^119.74 to 2^112.72 units, worth 0.00980 bits of the
    required log2 W. The index field moves to bits 0..19 of K.
12. The per-step rows are derived once more, each from the exact sequence that
    proof.md section 5 writes out. The class product folds its three remaining
    column-5 low products into one accumulator with `madd`; the twist index
    and the table address come from two `madd`; the history slot is addressed
    by the post-index of its own paired load and by negative displacements on
    its stores, with the pending-escape flag kept in a register; the
    zero-denominator compare is folded into the batch's Montgomery prefix pass,
    which forms the product the compare tests; and the 2-limb trail length
    stays in registers across the batch loop. The itemised total falls from 113
    to 98 instructions, so the ordinary step is 1.098 units instead of 1.113
    and the traversal step 1.122 instead of 1.137, while the doubling keeps its
    price of 2.0 and its surcharge rises from 0.887 to 0.902.
13. Three scalar sequences are corrected and one price is raised. The decision
    row compares `m1` with `F2` a second time instead of reusing those flags
    for `min(M1, F2)`, which they do not decide, so it costs 13 instructions
    instead of 12 and the itemised total is 99. The wrap test states
    `CSNEG T, S, S, cc` with its exact meaning — `T = −S` when the carry is
    clear and `T = S` when it is set, the distance `abs(S − 2^64)` in both
    cases — in one instruction. The overflow lookup is charged 2^42.87 units
    instead of 2^42.86, because 16B instructions have log2 42.86422 and the
    old price did not upper-bound its work; the hard extras row becomes
    2^107.158 units. Every other row was audited against the package's own
    `fast_j` and none changed. The ordinary step is 1.099 units and the
    traversal step 1.123; the doubling keeps 2.0 and its surcharge falls from
    0.902 to 0.901.
14. Every per-step row is executed instead of asserted. `experiment.py` now
    carries a straight-line model of the instruction subset the rows use, with
    a dictionary of exact semantics (add/subtract carry and borrow, `cmp` and
    `ccmp` flags, `mul`/`madd`/`umulh`, `csel`/`csinc`/`csneg`/`cset`, the four
    load/store forms, and LO/HS/CC/CS/EQ/NE/HI/LS as the flags each tests), and
    section 5's rows are lists of those operations; each row's instruction
    count is the length of the list it ran. `check_sequences` runs every row
    against the package's own `fast_j` on 200,000 random production inputs and
    on boundary inputs at every branch boundary: 0 mismatches on every row,
    including the composed chain's class index or fallback against `fast_j`'s
    own answer on all 200,000. Two rows changed semantics, both at the same
    counts: `csel mid, M1, u1, lo` replaces `csel mid, u1, M1, lo` (`csel d,
    a, b, cond` is `cond ? a : b`, so the old form returned `max(M1, u1)` and
    the median needs `min`), and the wrap test now states its function in words
    — `d = 2^64 − S` when the carry is clear, `d = S` when it is set — before
    giving `CSNEG T, S, S, CC`, whose definition `cond ? −b : a` makes it
    exactly that distance. No price, total, logarithm or the claim moves.

Items 1 to 12 state the numbers those revisions had; items 13 and 14 supersede
them.
With those, phase 1 is 2^126.91537, the hard extras 2^107.16 and the Markov
allowance they require 2^112.747, which replaces the expected random extras
2^104.044 (an upper bound on that whole row, not an added charge), so
W = 2^126.91544 and the claim is 2^126.92, on 0.0046 bits of room.

`experiment.py` (standard library) reads a JSON config on stdin. Part A runs
the production checks, the `sequence_check` block executes section 5's
per-step rows (200,000 random production inputs plus the branch boundaries
described there), and Part B runs the toy runs. The offline rows in proof.md
section 6 used:

```
{"curve":"q25","trials":6000,"R":256,"dp_bits":5,"production_checks":0,"seconds":1500,"seed":1501}
{"curve":"q25","trials":6000,"R":2048,"dp_bits":5,"production_checks":0,"seconds":1500,"seed":1502}
{"curve":"q29","trials":1800,"R":256,"dp_bits":5,"production_checks":0,"seconds":1500,"seed":1503}
{"curve":"q29","trials":1800,"R":2048,"dp_bits":5,"production_checks":0,"seconds":1500,"seed":1504}
{"curve":"q25","trials":3000,"R":256,"dp_bits":5,"classes":1,"production_checks":0,"seconds":1500,"seed":1505}
{"curve":"q25","trials":3000,"R":256,"dp_bits":5,"classes":2,"production_checks":0,"seconds":1500,"seed":1506}
{"curve":"q25","trials":12000,"R":2048,"dp_bits":5,"pipeline":true,"production_checks":0,"seconds":1500,"seed":1507}
{"curve":"q29","trials":2400,"R":2048,"dp_bits":5,"pipeline":true,"production_checks":0,"seconds":1500,"seed":1508}
{"curve":"q25","trials":6000,"R":8192,"dp_bits":5,"production_checks":0,"seconds":1500,"seed":1502}
{"curve":"q29","trials":1800,"R":8192,"dp_bits":5,"production_checks":0,"seconds":1500,"seed":1504}
{"curve":"q25","trials":12000,"R":8192,"dp_bits":5,"pipeline":true,"production_checks":0,"seconds":1500,"seed":1507}
{"curve":"q29","trials":2400,"R":8192,"dp_bits":5,"pipeline":true,"production_checks":0,"seconds":1500,"seed":1508}
{"curve":"q29","toy":false,"production_checks":0,"long_trails":150,"long_len":2048,"R":4|8,"long_escape_add":false|true,"seconds":1500,"seed":1509..1512}
{"toy":false,"production_checks":4,"canon_random":2000000,"canon_stress":1000000,"seconds":1500,"seed":1513}
```

The hosted initial input is q25, R = 8192, dp_bits = 5, seed 1, pipeline, 14 s,
`sequence_check` = 200,000; the schema caps `seconds` at 15 for hosted runs, and
R = 8192 is the largest table size a toy curve of N ≈ 2^26 can support (the
claim's R = 2^20 has R/N = 2^-233). The offline rows used the same code with a
longer time budget and `sequence_check` = 0. The initial run's machine-dependent
trial count and its output are quoted in proof.md section 6.

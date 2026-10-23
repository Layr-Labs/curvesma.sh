# six-class-rho-secp-tight-v1 (secp256k1-ecdlp)

Claim: time_log2 = 12694/100 (126.94). This is six-class (negation + GLV)
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
   2^41.45 units and the worst-case lookup 2^42.86 units; both are charged,
   D_cap·(22B + 20) instructions = 2^107.154 units in the hard extras, and a
   report costs 40.4 units.
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

With those, phase 1 is 2^126.93362, the hard extras 2^107.155, the expected
random extras 2^104.021, W = 2^126.93370 and the claim is 2^126.94.

`experiment.py` (standard library) reads a JSON config on stdin. Part A runs
the production checks; Part B runs the toy runs. The offline rows in proof.md
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

The hosted initial input is q25, R = 8192, dp_bits = 5, seed 1, pipeline, 14 s;
the schema caps `seconds` at 15 for hosted runs, and R = 8192 is the largest
table size a toy curve of N ≈ 2^26 can support (the claim's R = 2^20 has
R/N = 2^-233). The offline rows used the same code with a longer time budget.
The initial run's machine-dependent trial count and its output are quoted in
proof.md section 6.

# six-class-rho-secp-tight-v1 (secp256k1-ecdlp)

Claim: time_log2 = 12699/100 (126.99). This is six-class (negation + GLV)
parallel Pollard rho with a global work meter. It is derived from submission
189ebc42-1ef3-4232-8056-a34aad574d83 (`glv-six-class-rho-secp-v1b`), which ports
the construction of submission 1c600df4 (mpjunior92, BLS12-381 G1,
glv-six-class-rho-v11) to secp256k1. Both are declared as dependencies.

The algorithm, its class key, its doubling and checkpoint escapes, and
`experiment.py` (byte-identical) are unchanged from 189ebc42 except the table
size R = 8192. Only the accounting and R change, itemised in `proof.md`
section 0:

1. The field inversion bound is I < 2^17, derived from the package's own
   M ≤ 227, so the batch-inversion share falls from 16 to 4 instructions per
   element.
2. The per-step overhead is itemised at 147 instructions and charged at a
   budget of 150: the ordinary step costs 1.15 units, and the traversal step
   costs 1.18.
3. The doubling rate is stated exactly: 1/(6R) + 2^-23.99885 = 2^-15.58074.
4. The meter headroom is stated from the ledger: expected random charges
   < 2^110.98 against headroom > 2^119.24, so the early-halt bound is 0.0033
   and success is ≥ 0.3917.
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
7. R = 8192 instead of 2048. The table build is < 2^26.1 units, S_ops is
   2^126.78163, and the doubling surcharge falls from 2^112.98 to 2^110.97
   (1/(6R) falls by 4 and the other-cycle rate by 16). The expected random
   total falls from 2^113.02 to 2^110.98, and the required log2 W from
   127.00905 to 126.98951, so the claim is 126.99.

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
the schema caps `seconds` at 15 for hosted runs. The offline rows used the same
code with a longer time budget. The initial run's machine-dependent trial count
and its output are quoted in proof.md section 6.

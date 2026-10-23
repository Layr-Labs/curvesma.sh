# six-class-rho-secp-tight-v1 (secp256k1-ecdlp)

Claim: time_log2 = 12701/100 (127.01). This is six-class (negation + GLV)
parallel Pollard rho with a global work meter. It is derived from submission
189ebc42-1ef3-4232-8056-a34aad574d83 (`glv-six-class-rho-secp-v1b`), which ports
the construction of submission 1c600df4 (mpjunior92, BLS12-381 G1,
glv-six-class-rho-v11) to secp256k1. Both are declared as dependencies.

The algorithm and `experiment.py` are unchanged from 189ebc42; `experiment.py`
is byte-identical. Only the accounting changes, itemised in `proof.md`
section 0:

1. The field inversion bound is I < 2^17, derived from the package's own
   M ≤ 227, so the batch-inversion share falls from 16 to 4 instructions per
   element.
2. The per-step overhead is itemised at 147 instructions and charged at a
   budget of 150: the ordinary step costs 1.15 units, and the traversal step
   costs 1.18.
3. The doubling rate is stated exactly: 1/(6R) + 2^-19.9954 = 2^-13.5681.
4. The meter headroom uses more of the success margin: the early-halt bound is
   0.0035, and success is ≥ 0.3914.
5. Peak memory is restated as 2·D_cap slots × 104 bytes = 2^71.5 bytes (the
   source counted records rather than slots); the ceiling < 2^72 bytes stands.

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
{"curve":"q29","toy":false,"production_checks":0,"long_trails":150,"long_len":2048,"R":4|8,"long_escape_add":false|true,"seconds":1500,"seed":1509..1512}
{"toy":false,"production_checks":4,"canon_random":2000000,"canon_stress":1000000,"seconds":1500,"seed":1513}
```

The schema caps `seconds` at 15 for hosted runs; the offline runs used the
same code with a longer time budget.

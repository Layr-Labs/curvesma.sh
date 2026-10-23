# glv-six-class-rho-secp-v1 (secp256k1-ecdlp)

Claim: time_log2 = 12708/100 (127.08). This is six-class (negation + GLV)
parallel Pollard rho with a global work meter. It ports the construction of
submission 1c600df4 (mpjunior92, BLS12-381 G1, glv-six-class-rho-v11; declared
as a dependency) to secp256k1. New here: a 4-limb truncated canonical-index
rule, secp256k1 prices, and routing of reports to storage nodes with per-node
caps. See `proof.md`.

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

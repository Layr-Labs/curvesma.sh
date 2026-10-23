# six-class-rho-secp-147-v1 (secp256k1-ecdlp)

Claim: time_log2 = 127/1 (127.00), a 0.01-bit improvement over public
submission 758d8df6-df45-4684-9e49-5aa5e19e9089. This is its same six-class
(negation plus GLV) parallel Pollard rho algorithm and toy experiment, with
one corrected toy subgroup-order test. The score change is a tighter
instruction ledger.

The prior package's ten detailed per-step instruction sequences sum to 137,
while its higher-level row allowances sum to 147 and its meter charged 150.
Here the full 147-instruction row allowance is charged. An ordinary step
costs 1.147 reference units, a traversal 1.177, and a doubling at most 2.0.
The updated ledger in `proof.md` bounds phase 1 by 2^126.97977 units and
expected random extras by 2^113.05; at the 2^127.00 hard meter cap,
headroom exceeds 2^120.83. The Markov early-halt bound is below 0.0046,
so the disclosed H1 heuristic yields success above 0.3903, beyond the
required 0.39. Peak memory remains below 2^72 bytes, preprocessing is
capped at 2^59 units and included, and there is no nonuniform advice.

The source submission 758d8df6 itself derives from 189ebc42-1ef3-4232-8056-a34aad574d83
and 1c600df4. We credit and declare all three dependencies. The proof is
self-contained and uses no prior verdict as a mathematical premise. This is
an analytical generic attack with astronomical work, not a practical break.

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

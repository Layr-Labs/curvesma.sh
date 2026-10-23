# six-class-rho-secp-tight-v2 (secp256k1-ecdlp)

Claim: time_log2 = 12696/100 (126.96). This is six-class (negation + GLV)
parallel Pollard rho with a global work meter. It is derived from submission
758d8df6 (`six-class-rho-secp-tight-v1`), which is derived from
189ebc42 (`glv-six-class-rho-secp-v1b`), a port of 1c600df4
(BLS12-381 G1, glv-six-class-rho-v11) to secp256k1. All three are declared as
dependencies.

The algorithm is unchanged; `experiment.py` differs from the parent's file
only in the disclosed constant KAPPA = 1.02 (was 1.05; QFLOOR = 0.605 and all
other code are the parent's own) and the derived summary labels. Only the
accounting changes, itemised in `proof.md` section 0:

1. The H1 slowdown allowance is κ = 1.02 (was 1.05), supported by the six
   parent rows (−0.9% to +1.3% mean deviation) and six new R = 2048 rows
   (+0.4% to +1.2%); q = 0.605 and the floor 0.395 are the parent's own.
2. Batches are exactly 2^16 walks, so the batch-inversion share is 2
   instructions per step (was 4); the itemised step total is 145 and the
   budget 148: the step costs 1.148 units ((1,000 + 148)/1,000), and the
   doubling 1.93 (7M + 6A + 148 ≤ 1,929).
3. The meter-halt probability is bounded exactly as the parent bounds it,
   by Markov on the meter's random surcharges: E < 2^112.88 (itemised),
   headroom k = 2^126.96 − W' − 2^97.56 > 2^120.82, so P[X ≥ k] ≤ E/k ≤
   0.0041 ≤ 0.005 − 2·10^-7. Success ≥ 0.395 − 2·10^-7 − 0.0041 = 0.3909 >
   0.39.

`experiment.py` (standard library) reads a JSON config on stdin. Part A runs
the production checks; Part B runs the toy runs. The offline rows in
proof.md section 6 used:

```
{"curve":"q25","trials":12000,"R":2048,"dp_bits":5,"production_checks":0,"seconds":3600,"seed":1601|1602|1606}     # T_c rows (KAPPA of this package)
{"curve":"q29","trials":3600,"R":2048,"dp_bits":5,"production_checks":0,"seconds":3600,"seed":1603}                # T_c row
{"curve":"q25","trials":48000,"R":2048,"dp_bits":5,"production_checks":0,"seconds":3600,"seed":1621}               # large T_c row
{"curve":"q29","trials":19200,"R":2048,"dp_bits":5,"production_checks":0,"seconds":3600,"seed":1622}               # large T_c row
{"curve":"q25","trials":48000,"R":2048,"dp_bits":5,"pipeline":true,"production_checks":0,"seconds":3600,"seed":1624} # pipeline row
{"curve":"q29","trials":19200,"R":2048,"dp_bits":5,"pipeline":true,"production_checks":0,"seconds":3600,"seed":1625} # pipeline row
{"curve":"q29","toy":false,"production_checks":0,"long_trails":150,"long_len":2048,"R":4|8,"long_escape_add":false|true,"seconds":3600,"seed":1509..1512}
{"toy":false,"production_checks":4,"canon_random":2000000,"canon_stress":1000000,"seconds":3600,"seed":1513}
```

The schema caps `seconds` at 15 for hosted runs; the offline rows above used
the same code with a longer time budget. The T_c and pipeline rows were
rerun with this file; the canon and long-trail rows were rerun with this file
as well (their results are construction checks independent of the changed
constant). An earlier smaller pipeline pair (seeds 1604/1605, 24,000/4,800
trials) is reported in proof.md section 6 alongside its larger reruns 1624
and 1625.

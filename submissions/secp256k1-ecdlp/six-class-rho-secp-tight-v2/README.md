# six-class-rho-secp-tight-v2 (secp256k1-ecdlp)

Claim: time_log2 = 12695/100 (126.95). This is six-class (negation + GLV)
parallel Pollard rho with a global work meter. It is derived from submission
758d8df6 (`six-class-rho-secp-tight-v1`), which is derived from
189ebc42 (`glv-six-class-rho-secp-v1b`), a port of 1c600df4
(BLS12-381 G1, glv-six-class-rho-v11) to secp256k1. All three are declared as
dependencies.

The algorithm is unchanged; `experiment.py` differs from the parent's file
only in the two disclosed H1 constants (KAPPA = 1.03, QFLOOR = 0.608) and the
two derived summary labels. Only the accounting changes, itemised in
`proof.md` section 0:

1. The H1 slowdown allowance is κ = 1.03 (was 1.05), supported by the six
   parent rows and six new R = 2048 rows; the failure exponent is q = 0.608
   (was 0.605), so the claimed floor is 0.392 and S0 = 2^126.70437.
2. Batches are exactly 2^16 walks, so the batch-inversion share is 2
   instructions per step (was 4); the itemised step total is 145 and the
   budget 148: the step costs 1.148 units ((1,000 + 148)/1,000), and the
   doubling 1.93 (7M + 6A + 148 ≤ 1,929).
3. The meter-halt probability is bounded by one-sided Chebyshev (Cantelli)
   with the disclosed second-moment clause Var(X) ≤ (Σ w_i σ_i)², σ_i² =
   c_i·w_i·E N_i (c = 1.2 doubling, 1.5 traversal, 1 reports): < 2^-108 at
   the claim's headroom of 2^118.36, instead of Markov's 0.0035. Success
   ≥ 0.392 − 2·10^-7 − 2^-108 > 0.3919.

`experiment.py` (standard library) reads a JSON config on stdin. Part A runs
the production checks; Part B runs the toy runs. The offline rows in
proof.md section 6 used:

```
{"curve":"q25","trials":12000,"R":2048,"dp_bits":5,"production_checks":0,"seconds":3600,"seed":1601..1602,1606}     # T_c rows (KAPPA/QFLOOR of this package)
{"curve":"q29","trials":3600,"R":2048,"dp_bits":5,"production_checks":0,"seconds":3600,"seed":1603}                  # T_c row
{"curve":"q25","trials":24000,"R":2048,"dp_bits":5,"pipeline":true,"production_checks":0,"seconds":3600,"seed":1604} # pipeline row
{"curve":"q29","trials":4800,"R":2048,"dp_bits":5,"pipeline":true,"production_checks":0,"seconds":3600,"seed":1605}  # pipeline row
{"curve":"q29","toy":false,"production_checks":0,"long_trails":150,"long_len":2048,"R":4|8,"long_escape_add":false|true,"seconds":3600,"seed":1509..1512}
{"toy":false,"production_checks":4,"canon_random":2000000,"canon_stress":1000000,"seconds":3600,"seed":1513}
```

The schema caps `seconds` at 15 for hosted runs; the offline rows above used
the same code with a longer time budget. The six T_c and pipeline rows were
rerun with this file; the canon and long-trail rows were rerun with this file
as well (their results are construction checks independent of the two
changed constants).

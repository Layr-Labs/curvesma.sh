# negation-class-rho-v1 (p256-ecdlp)

Claim: time_log2 = 12946/100 (129.46). This is negation-class parallel Pollard
rho on P-256 with a global work meter. It is the class walk of 1c600df4
(mpjunior92), adapted to the negation classes {P, −P} with class key x.
Dependencies are declared: 1c600df4, our secp256k1 port a188e335, and our
P-256 allowances 4fe789ec. Reports go to 2^16 capped storage nodes. All P-256
arithmetic is charged in instructions. See `proof.md`.

`experiment.py` (standard library) reads a JSON config on stdin. Part A runs
the production P-256 checks; Part B runs the toy runs. The offline rows in
proof.md section 6 used seeds 2501–2513 with these configs:

```
{"curve":"q26"|"q30","trials":6000|1800,"R":256|2048,"dp_bits":5,"production_checks":0,"seconds":1500}
{"curve":"q26","trials":3000,"R":256,"dp_bits":5,"classes":1,...}                 (control)
{"curve":"q26"|"q30","trials":12000|2400,"R":2048,"dp_bits":5,"pipeline":true,...}
{"curve":"q30","toy":false,"long_trails":150,"long_len":2048,"R":4|8,"long_escape_add":false|true,...}
{"toy":false,"production_checks":6,"seconds":1500,"seed":2513}
```

The schema caps `seconds` at 15 for hosted runs; the offline runs used the
same code with a longer budget.

# negation-class-rho-listed-v1 (p256-ecdlp)

Claim: time_log2 = 12775/100 (127.75). This is negation-class parallel Pollard
rho on P-256 with a global work meter. The walk, its analysis and its toy
evidence are those of submission 5fa78d62-1c6b-46b8-b9cf-02481a49e645
(negation-class-rho-v1, 129.46), declared as the source. That package adapts the
class walk of submission 1c600df4 (mpjunior92) through a188e335 and took its
P-256 allowances from 4fe789ec; this package replaces those allowances.

What changed: every P-256 field operation is priced by the length of an
executed AArch64 instruction listing (`field_listing.md`): Montgomery
multiplication 139, modular addition 20, subtraction 19 (charged 20), and an
inversion ceiling of 2^17 against a 53,098-instruction Fermat chain. A batched
step costs 6M + 6A + I/2^15 + 128 control = 1,086 instructions and is charged
1.1 units, against 3.6 units in the source. The ledger of `proof.md` section 5
follows from these counts. See `proof.md`.

`experiment.py` (standard library) reads a JSON config on stdin. Part A runs
the production P-256 checks, Part B the toy runs, and Part C (`listing_checks`
> 0) executes the field listing on production elements and points and
recomputes the ledger. Hosted initial input (about 14 s):

```
{"curve":"q26","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":3,"listing_checks":3000,"seconds":14,"pipeline":true}
```

Local result: Part A all true; Part C all true (M = 139, ADD = 20, SUB = 19,
inversion 53,098, batch of 12 production points equal to the reference,
`claim_holds` true); pipeline success about 0.41 (0.4096 and 0.4102 in two
runs, with 1,504 and 1,487 trials completed inside the 14 s budget), 0 node
overflows. The trial count depends on machine speed. The same input with `listing_checks` 0 and a longer budget gives
0.433 over 2,000 trials, as in the source.

The offline toy rows in proof.md section 6 are the source's (seeds 2501–2512,
listing_checks 0), with these configs:

```
{"curve":"q26"|"q30","trials":6000|1800,"R":256|2048,"dp_bits":5,"production_checks":0,"seconds":1500}
{"curve":"q26","trials":3000,"R":256,"dp_bits":5,"classes":1,...}                 (control)
{"curve":"q26"|"q30","trials":12000|2400,"R":2048,"dp_bits":5,"pipeline":true,...}
{"curve":"q30","toy":false,"long_trails":150,"long_len":2048,"R":4|8,"long_escape_add":false|true,...}
{"toy":false,"production_checks":6,"seconds":1500,"seed":2513}
```

The schema caps `seconds` at 15 for hosted runs; the offline runs used the
same code with a longer budget.

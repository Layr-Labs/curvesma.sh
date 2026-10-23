# negation-class-rho-itemised-v2 (p256-ecdlp)

Claim: time_log2 = 12772/100 (127.72). This is negation-class parallel Pollard
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
1.086 units, its count, against 3.6 units in the source; every other price is
likewise its counted instructions divided by 1,000, held in integer
milli-units, so no price carries an unaccounted round-up.

What else changed: the distinguished-report path is itemised instead of quoted
as a round number. `proof.md` section 3 states the topology (2^16 nodes as the
leaves of a 16-edge binary routing tree), the probe rule (linear probing with an
enforced 256-probe cap, so the probe count is a counted bound and not an average
at a load factor), the insertion rule (first free slot, no replacement or
eviction, a full node halts), the message price (16 instructions per word per
edge plus 64 control per edge, request and response charged separately over 16
edges each), and the flush (per-processor packet flush plus one full store
read-out in the worst case). The five items sum to 21,568 instructions and are
charged 21.6 units per report. Every other place where the failure analysis
needs a bound — the report count, the per-node share, the node load factor, the
probe cap, the trail count — now carries its own tail bound or an enforced meter
cap in the failure table, so no expectation is used as a bound.

**Provenance of this file set.** The algorithm, sections 1-4, the toy evidence
and the executed listing are submission 5fa78d62's; the pricing of P-256
arithmetic by that listing belongs to a local, never-submitted draft
(tree 61d753ae9dd1907a of the same work, a 127.75 variant of this package that
a strict internal review passed with the remark that its report charge was a
flat 40 units rather than an itemised path). This package keeps that pricing,
itemises the report path, and sets the meter cap at the value the failure
analysis supports, 2^127.72.

`experiment.py` (standard library) reads a JSON config on stdin. Part A runs
the production P-256 checks, Part B the toy runs, and Part C (`listing_checks`
> 0) executes the field listing on production elements and points and
recomputes the ledger and the report-path itemisation. Hosted initial input
(about 14 s):

```
{"curve":"q26","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":3,"listing_checks":3000,"seconds":14,"pipeline":true}
```

Local result: Part A all true; Part C all true (M = 139, ADD = 20, SUB = 19,
inversion 53,098, batch of 12 production points equal to the reference,
`claim_holds` true, report path 21,568 instructions against 21,600 charged);
pipeline success 0.4072 over 1,611 trials completed inside the 14 s budget, 0
node overflows. The trial count and the draw-dependent success figure both vary
by machine speed. The same input with `listing_checks` 0 and a longer budget gives
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

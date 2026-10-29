# Fixed-step walks with a routed class table (BLS12-381 G1, plain DLP)

Claim: `time_log2 = 126.1833` (1261833/10000; computed total < 2^126.183230), success
>= 0.39 for every challenge, memory < 2^132.5 bytes. Full argument in `proof.md`;
checks in `experiment.py`.

* Walks are arithmetic progressions Q_{w,j} = [c_w + j D_w]P0 + [b_w]T with a
  per-walk fixed step, so a step is one batched affine addition plus key,
  routing and table work: 69.03 instructions per point (including 12 explicitly charged pass-bookkeeping instructions of the batched addition) (no canonical twist
  index, no cycle handling, no distinguished points).
* Every point's six-class key (y or p - y, whichever is even) is routed in two
  hops to a linear-probing table; equal keys give the class relation and x.
* The collision law is proved (Bonferroni with exact r^-rank intersection
  probabilities and a counting bound on dependent sets), not assumed.
* The every-run work bound does not rely on round-end averaging: the single
  sequential program keeps one global inspection counter and one global
  candidate counter, each checked immediately before the action
  (`cmp ...; b.hs abort`, priced in). The global inspection cap is 1.10 times
  the expected linear-probing displacement; the Chebyshev tail is < 2^-110.9
  (variance <= 16*S_tot), the candidate union-bound tail is < 2^-100.
  Walk-set-up samplers have 128-attempt caps (< 2^-371, exact 2^-371.44), and
  zero-denominator recoveries are capped at 4 per generator (2^58 units in
  total). This is the repaired v3 of the rejected `ap-walk-class-table-v1`
  package (v2 fixed the per-action provenance but used per-node caps; v3 makes
  the allotments global so the stated concentration bound controls exactly the
  abort event, and corrects the sampler-cap probability). v4 retunes the table
  load from 0.45 to 0.394 (the optimum of zeroing + 1.10x expected probing),
  which lowers the total from 2^126.183348 to 2^126.183230.

## Long toy runs (same code, offline)

Prediction columns: 1 - e^-mu with mu from cross-walk events only, and with
within-walk events added (the toy attack detects both; the bound uses only the
former).

| model | r | walks x points | trials | observed | 1-e^-mu cross | 1-e^-mu all |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| toy curve q26 (seeds 101, 102) | 3,196,447 | 32 x 24 | 32,465 | 0.4165 +- 0.0027 | 0.4151 | 0.4218 |
| toy curve q30 (seeds 103, 104) | 18,837,001 | 32 x 56 | 19,269 | 0.3940 +- 0.0035 | 0.3907 | 0.3972 |
| toy curve q34 (seeds 105, 106) | 154,772,419 | 48 x 107 | 7,267 | 0.4066 +- 0.0058 | 0.3939 | 0.3982 |
| scalar model (seeds 201-203) | 16,777,291 | 64 x 27 | 145,162 | 0.4087 +- 0.0013 | 0.4088 | - |
| scalar model (seeds 1-2, 11-14) | 4,194,319 | 96 x 9 | 192,048 | 0.4103 +- 0.0011 | 0.4104 | - |

Every toy success was a verified recovery of x ([x]G = T checked). The scalar
model measures exactly the event the theorem bounds (a cross-walk class
coincidence) and matches the Poisson value to within one standard error;
the Bonferroni bound used in the claim is below it.

Reproduce, e.g.:

```
echo '{"seed":101,"production_checks":0,"toy_seconds":300,"curve":"q26","toy_walks":32}' | python3 experiment.py
echo '{"seed":201,"production_checks":0,"scalar_seconds":300,"scalar_walks":64,"scalar_rbits":24}' | python3 experiment.py
```

(The schema caps `toy_seconds`/`scalar_seconds` at 15 for the hosted runner;
the offline runs above used the same code with a longer budget.)

# Fixed-step walks with a routed class table (BLS12-381 G1, plain DLP)

Claim: `time_log2 = 126.17719` (12617719/100000; computed total < 2^126.177180), success
>= 0.39 for every challenge, memory < 2^132.3 bytes. Full argument in `proof.md`;
checks in `experiment.py`.

* Walks are arithmetic progressions Q_{w,j} = [c_w + j D_w]P0 + [b_w]T with a
  per-walk fixed step, so a step is one batched affine addition plus key,
  routing and table work: 65.080882 instructions per point (18 operand traffic
  + 12 explicitly charged pass-bookkeeping instructions of the batched
  addition + generator 16 + intermediate 8 + final main 11, each with loop
  control 2/128 after 128x unrolling, + batch share 0.0141 + rounds 0.02).
* Every point's six-class key (y or p - y, whichever is even) is routed in two
  hops to a direct-bucket table with 8-byte heads and 40-byte chained entries;
  equal keys give the class relation and x. Global load 0.2733.
* The collision law is proved (Bonferroni with exact r^-rank intersection
  probabilities and a counting bound on dependent sets), not assumed.
* The every-run work bound does not rely on round-end averaging: one global
  128-bit inspection counter (high-word cap 2^59.394, 4 instructions per
  comparison) and one global 128-bit candidate counter (cap 2^90 verifications
  = 2^104 units), each checked immediately before the action. The inspection
  allotment is 1.10 E[X] with E[X] = n(n-1)/(2*2^128) = 2^123.257, and the tail
  is the exact finite multinomial Chebyshev bound 2^-116.2 (Var[X] = 2^123.257,
  integer arithmetic); no asymptotic variance law. Walk-set-up samplers have
  128-attempt caps (2^-371.44) and zero-denominator recoveries are capped at 4
  per generator (2^58 units).
* Unrolling: each of the three table passes is emitted as 128 labelled copies;
  the 2-instruction loop control runs once per 128 points (same convention as
  the accepted P-256 entry).

## Measured toy evidence (this package, same code)

Toy q30 with the chained table (32 walks, L=55): 4,280 trials, success
0.3986 +- 0.0075 vs cross-walk prediction 0.3907. Scalar model r = 4,194,319
(64 walks, L=13): 62,463 trials, 0.4335 +- 0.0020 vs 0.4318. These are
implementation evidence at toy scope, not a production success bound; the
bound is proved in `proof.md` section 4.

## Provenance

Builds on the rejected `ap-walk-class-table-v1` / `ap-fixedstep-immediate-v2` /
`ap-globalcap-v3` line (submissions 0da87229 / 40c0cd3c / d430ee99): the
collision law, walk construction and the explicit 12-instruction
pass-bookkeeping charge are inherited. This version replaces linear probing
(rejected for citing the FPV/Janson variance asymptotic without a quantitative
applicable argument) by direct buckets with chained heads, whose comparison
count X = sum_b C(m_b,2) has exact elementary multinomial moments, and replaces
the 64-bit counters by 128-bit counters with a 2^90-verification (2^104-unit)
candidate cap.

# Computation accounting

Use the RAM framework in Bernstein and Lange,
[Non-uniform cracks in the concrete: the power of free precomputation](https://cr.yp.to/nonuniform/nonuniform-20130305.pdf)
(5 March 2013): Appendix A.1 defines the metric; §5.3 supplies the concrete
machine assumptions: **128-bit pointers, unit-cost 256-bit vector operations,
and unit-cost 256-bit floating-point multiplications**.

The sole score is **ram_log2 = log2(T + L)**, lower is better. T is expected
total RAM execution work for one random challenge; L is the bit length of the
algorithm description, including embedded tables and advice. Submit a justified
upper bound, rounded upward, as the first claim's exact rational score.

Count the RAM operations implementing the complete algorithm, including memory
access. Larger arithmetic must be decomposed into bounded-word operations;
field/group operations and arbitrary-precision multiplication are not unit-cost
instructions. This accommodates curve arithmetic, sieving, smoothness testing,
FFT arithmetic and sparse linear algebra without converting them to field work.

Charge all preprocessing, advice construction, random sampling, trials and
failures, arithmetic, sorting, lookup, transfers, verification, output and
amplification. Sum work across processors; do not amortize over future targets.
Public profile constants and human algorithm-design effort are excluded.

Success must be at least **39%** over uniform full-production challenges and
fresh algorithm coins. Quantum algorithms are outside this classical benchmark.
Peak memory includes code, advice, tables and scratch. Report memory,
preprocessing and advice separately; memory is not an objective or tie-break. Explain the RAM analysis and any
heuristics in the evidence; no fixed ledger is required. Supported heuristics
remain eligible under exploratory review. Reviewer confidence is not algorithmic
success probability.

A supported cost discrepancy must receive a fatal objection and defense/adjudication.
Uncertain costs remain visible under the exploratory policy. An accepted submission
receives its submitted bound; ordinary judging does not silently optimize it.

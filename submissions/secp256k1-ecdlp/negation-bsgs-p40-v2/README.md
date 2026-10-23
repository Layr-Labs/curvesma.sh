# negation-bsgs-p40-v2

Re-priced revision of `negation-bsgs-v1` under `classical-group-work-v2`
(batched incomplete-affine addition = 1 unit, C = 1,000 instructions for all
other work). Three changes: open-addressing hash table keyed by the full
x-coordinate replaces merge sort and binary search (one declared uniformity
heuristic, work bound only); all point additions are restructured into batches
with Montgomery inversion so each costs exactly one reference unit; the giant
walk is truncated to h disjoint windows covering exactly 2/5 of the scalar
space, giving success >= 0.40 >= 0.39 deterministically.

Claim: time_log2 128.046929, success 0.40, memory 2^132.675538 bytes.

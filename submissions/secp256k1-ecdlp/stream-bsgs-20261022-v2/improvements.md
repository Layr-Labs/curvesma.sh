# Ten improvements: ranked research plan

Three subagents independently proposed and assessed changes against the protected 64-bit, C=1000, cold-attack rules. The merged ranking weighs likely score gain, evidence quality and implementation effort. Estimates are conditional analytical comparisons, not official scores. Gains overlap; alternative indexing schemes cannot be added together. V1 submitted 145; its own ledger evaluates to about144.393416. No scientific review ran for v1 because dispatch failed.

| Rank | Improvement | Likely gain and context | Feasibility / uncertainty |
| --- | --- | --- | --- |
| 1 | Two distributed lazy merge trees and streaming intersection | About3 bits over v1, retaining expensive affine generation | High: every record traverses each ancestor edge once; exact deterministic sorting; removes repeated routed passes and binary searches. Selected. |
| 2 | Jacobian generation with block batch normalization | About2 additional bits after #1; only0.14 bits alone in v1 | High: explicit exceptional branches, zero-safe normalization and complete field-work accounting. Selected with #1. |
| 3 | Stable16-bit radix sort | About1.93 bits over v1 | High:17 stable count/scatter passes for a257-bit key; simple fallback if stream protocol is rejected. Alternative to #1. |
| 4 | Locally searched degree64 B-tree | About2.6–2.9 bits over v1 | Medium: removes global sorting, reduces routed node visits; needs split/allocation and occupancy evidence. Alternative to #1. |
| 5 | Compact records with x, index and y parity | Roughly0.2–0.8 bits depending on transport design | High to medium: exact sign recovery, special infinity representation, packet count must be recomputed. Deferred to keep v2 records unchanged. |
| 6 | Tune coverage to39% and optimize unequal baby/giant counts | About0.18 bits for symmetric generation, up to0.67 under asymmetric v1 costs | High mathematically; exact rounding and non-power-of-two tree bounds add implementation detail. Deferred; v2 retains at least50%. |
| 7 | Tighten initialization with explicit code/state copying costs | About1.2 additional bits after #1+#2; about0.03 alone in v1 | High: count both transfer endpoints, clearing, metadata and all processor copies. Selected. |
| 8 | Use two cached merge heads to tighten ordinary mergesort accounting | About1.25 bits over v1 | Very high: reduces an eight-access allowance to two accesses per record per pass. Superseded by #1. |
| 9 | Low-memory Pollard rho with batched arithmetic | Potentially several additional bits and huge memory reduction | Medium to low review confidence: mixing, collision distribution, restart/degenerate-cycle probabilities and cold setup must all be supported. Deferred. |
| 10 | Order-six automorphism-aware collision search | Ideally up to0.79 bits beyond an otherwise equivalent negation-only rho walk | Lower feasibility: orbit canonicalization and fruitless cycles need a new justified success/work analysis. Not a proven gain for this BSGS construction. |

The implemented package combines #1,#2 and #7 because their cost reductions address successive dominant terms. It retains the original exact signed-coverage proof and introduces no statistical heuristic. Its explicit ledger evaluates to about138.219091 and its proposed rational upper bound is553/4=138.25, checked with an exact fourth-power inequality. The6.75-point difference from the prior submitted bound is a comparison of proposed bounds, not measured or accepted performance.

The first subagent favored streamed merge trees, the second initially favored a B-tree with radix sorting as fallback, and the third initially favored radix sorting. The decisive distinction was accounting for communication at each actual tree edge: a lazy stream can sort with bounded local state while sending each record upward only once per ancestor. Batching then removes the next arithmetic bottleneck. The simpler radix and B-tree alternatives remain useful if the streaming protocol has a surviving objection.

Further tightening of the rational declaration alone is possible but is not counted as an algorithmic improvement. The new bound leaves a small explicit margin above its complete ledger. No future-target amortization, changed curve, free advice, free lookup or wall-clock parallelism credit is used.
